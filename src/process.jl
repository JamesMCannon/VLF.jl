# ============================================================================
# process.jl  —  Tier 2: raw counts -> calibrated/cleaned products.
#
# Everything here operates on the canonical full-day grid, so the legacy
# time-matching / deleteat! gymnastics are gone: gaps are explicit NaNs and the
# two channels are aligned by construction.
# ============================================================================

# --- calibration ------------------------------------------------------------
# Resolve `cal_num` for a single channel.
#   * A scalar applies to both channels; the sentinel -1.0 means "unset" so the
#     caller falls through to `cal_file` (returns `nothing`).
#   * A NamedTuple selects per channel by name (:EW / :NS). It is always treated
#     as "set" and must contain the key for the channel being calibrated.
_cal_num_for(cal_num::Real, ::Channel) = cal_num == -1.0 ? nothing : float(cal_num)
function _cal_num_for(cal_num::NamedTuple, ch::Channel)
    name = Symbol(chstr(ch))            # :EW or :NS
    haskey(cal_num, name) ||
        error("Per-channel cal_num has no entry for the $name channel (keys: $(keys(cal_num))).")
    return float(cal_num[name])
end

"""
    calibration_factor(amp_day::RawDay, params) -> Float64

Counts→pT scale factor for `amp_day`'s channel. Uses `params.cal_num` when set:
a scalar applies to both channels, while a NamedTuple `(EW = …, NS = …)` selects
the factor for `amp_day.channel`. A scalar `-1.0` (the unset sentinel) instead
reads `params.cal_file`, selects the NS/EW curve by `amp_day.channel`, and takes
the response at the frequency nearest `amp_day.Fc` (cal-file frequencies are in
kHz). A NamedTuple `cal_num` always takes precedence over `cal_file`.
"""
function calibration_factor(amp_day::RawDay, params::ProcessParams)
    cn = _cal_num_for(params.cal_num, amp_day.channel)
    if cn !== nothing
        return cn
    elseif !isempty(params.cal_file) && params.cal_file != "default"
        mf = matopen(params.cal_file)
        try
            curve = amp_day.channel === NS ? read(mf, "CalibrationNumberNS") :
                                              read(mf, "CalibrationNumberEW")
            freqs = abs.(curve[:, 1])           # kHz
            resp  = abs.(curve[:, 2])
            idx = argmin(abs.(freqs .* 1000 .- amp_day.Fc))
            return resp[idx]
        finally
            close(mf)
        end
    else
        error("No calibration information provided (set cal_num or cal_file).")
    end
end

"""
    calibrate(amp_day::RawDay, params) -> Vector{Float64}

Amplitude in pT. NaN gaps stay NaN.
"""
calibrate(amp_day::RawDay, params::ProcessParams) =
    amp_day.data .* calibration_factor(amp_day, params)

# --- amplitude units (pT <-> µV/m, linear <-> dB) ---------------------------
# Far-field plane-wave relation E = c·B, using the speed of light in sea-level
# air (c_vacuum / n_air, n_air ≈ 1.000293 at standard sea-level conditions).
const C_VACUUM        = 299_792_458.0          # m/s (exact)
const N_AIR_SEA_LEVEL = 1.000293               # radio refractive index, sea level
const C_AIR           = C_VACUUM / N_AIR_SEA_LEVEL

# B in pT (1e-12 T) → E in µV/m (1e-6 V/m): E = c_air·B[T]·1e6 = c_air·1e-6·B[pT].
"µV/m of E-field per pT of B-field for a plane wave in sea-level air (≈ 299.70)."
const PT_TO_UVM = C_AIR * 1e-6

"""
    pT_to_uVm(amp_pT)

Convert a magnetic-field amplitude in pT to electric field in µV/m for a plane
wave in sea-level air (`E = c_air · B`; [`PT_TO_UVM`](@ref) ≈ 299.70 µV/m per pT).
Broadcasts over scalars or arrays; `NaN` gaps are preserved.
"""
pT_to_uVm(amp_pT) = amp_pT .* PT_TO_UVM

"""
    to_db(x) -> Float64

Field amplitude in dB, `20·log10(x)`. `NaN` is passed through, and non-positive
values map to `NaN` (a true zero would be −∞ dB, which would poison plots/stats).
"""
to_db(x::Real) = (isnan(x) || x <= 0) ? NaN : 20 * log10(x)

"""
    apply_amplitude_units(amp_pT; efield=false, db=false) -> Vector{Float64}

Express a linear-pT amplitude vector in other units: pT or µV/m (`efield`), then
linear or dB (`db`). `NaN` gaps are preserved. Quadrature combination must be
done in linear units *before* this is applied.
"""
function apply_amplitude_units(amp_pT::AbstractVector; efield::Bool = false, db::Bool = false)
    lin = efield ? amp_pT .* PT_TO_UVM : amp_pT
    return db ? to_db.(lin) : collect(float.(lin))
end

"""
    amplitude(day::ProcessedDay, which=:combined; efield=false, db=false) -> Vector{Float64}

Return one of `day`'s amplitude products in the requested units, converted at
read time from the canonical linear-pT cache (so switching units never triggers
a recompute). `which` is `:EW`, `:NS`, or `:combined`; units follow `efield`
(pT→µV/m) and `db` (linear→dB). See [`apply_amplitude_units`](@ref).
"""
function amplitude(day::ProcessedDay, which::Symbol = :combined;
                   efield::Bool = false, db::Bool = false)
    v = which === :EW       ? day.EW_amp :
        which === :NS       ? day.NS_amp :
        which === :combined ? day.combined_amp :
        error("`which` must be :EW, :NS, or :combined (got :$which)")
    return apply_amplitude_units(v; efield = efield, db = db)
end

"""
    amplitudes(day::ProcessedDay; efield=false, db=false) -> NamedTuple

All three amplitude products converted at read time:
`(EW = …, NS = …, combined = …)`. See [`amplitude`](@ref).
"""
amplitudes(day::ProcessedDay; efield::Bool = false, db::Bool = false) =
    (EW       = amplitude(day, :EW;       efield = efield, db = db),
     NS       = amplitude(day, :NS;       efield = efield, db = db),
     combined = amplitude(day, :combined; efield = efield, db = db))

# --- ProcessedView: a units projection of a ProcessedDay --------------------
# The cache is canonical linear pT; conversion to µV/m / dB is a *presentation*
# layer that lives entirely on top of `ProcessedDay`. A `ProcessedView` owns its
# three converted amplitude arrays and borrows everything else (phase, metadata)
# from its parent. It is never cached and never an input to arithmetic: combine
# / calibrate in linear pT first, then view. Units only ever flow one way
# (pT -> µV/m / dB), and re-viewing always re-derives from the pT parent, so a
# conversion can never stack on itself.

"""
    AmplitudeUnits(efield, db)

The two orthogonal amplitude-unit axes, mirroring [`apply_amplitude_units`](@ref):
`efield` selects magnetic B (pT, `false`) vs electric E (µV/m, `true`); `db`
selects linear (`false`) vs `20·log10` decibels (`true`).
"""
struct AmplitudeUnits
    efield::Bool
    db::Bool
end

"Human-readable label for [`AmplitudeUnits`](@ref), e.g. `\"pT\"`, `\"µV/m (dB)\"`."
unit_label(u::AmplitudeUnits) = string(u.efield ? "µV/m" : "pT", u.db ? " (dB)" : "")

"""
    ProcessedView

A read-only units projection of a [`ProcessedDay`](@ref). The three amplitude
fields (`EW_amp`, `NS_amp`, `combined_amp`) are converted from the parent's
canonical linear pT and owned by the view; phase and all metadata (`time`, `Fc`,
`Fs`, `EW_pha`, `NS_pha`, `params`, …) are read through to the parent unchanged,
so `view.EW_pha` works exactly like `day.EW_pha` (and shares its storage — treat
phase as read-only here). `view.parent` recovers the underlying `ProcessedDay`.

Built with [`view_units`](@ref). Never cached (a [`ProcessedView`](@ref) passed to
the cache errors) and never fed back into [`combine_quadrature`](@ref) /
[`calibrate`](@ref), which require linear pT and accept only `ProcessedDay`.

See also [`view_units`](@ref), [`get_processed_view`](@ref), [`amplitudes`](@ref).
"""
struct ProcessedView
    parent::ProcessedDay
    units::AmplitudeUnits
    EW_amp::Vector{Float64}
    NS_amp::Vector{Float64}
    combined_amp::Vector{Float64}
end

# Names the view answers for itself; everything else forwards to the parent.
const _VIEW_OWN = (:parent, :units, :EW_amp, :NS_amp, :combined_amp)
Base.getproperty(v::ProcessedView, n::Symbol) =
    n in _VIEW_OWN ? getfield(v, n) : getproperty(getfield(v, :parent), n)
Base.propertynames(::ProcessedView) = (fieldnames(ProcessedDay)..., :units, :parent)

"""
    view_units(day::ProcessedDay; efield=false, db=false) -> ProcessedView
    view_units(view::ProcessedView; efield=false, db=false) -> ProcessedView

Project `day` into the requested amplitude units, reusing [`amplitudes`](@ref)
(which converts each stored pT field independently — `combined` is converted, not
recombined). Viewing a [`ProcessedView`](@ref) re-projects its parent, so units
never stack.
"""
function view_units(day::ProcessedDay; efield::Bool = false, db::Bool = false)
    a = amplitudes(day; efield = efield, db = db)
    return ProcessedView(day, AmplitudeUnits(efield, db), a.EW, a.NS, a.combined)
end
view_units(v::ProcessedView; efield::Bool = false, db::Bool = false) =
    view_units(getfield(v, :parent); efield = efield, db = db)

"""
    amplitude(v::ProcessedView, which=:combined) -> Vector{Float64}

Return one of the view's already-converted amplitude products (`:EW`, `:NS`, or
`:combined`) as a fresh copy. Passing `efield`/`db` is an error: the view is
already fixed to `v.units`, and re-converting would double-apply the scaling —
re-view the parent (`view_units(v.parent; …)`) to change units instead.
"""
function amplitude(v::ProcessedView, which::Symbol = :combined;
                   efield::Bool = false, db::Bool = false)
    (efield || db) && error("amplitude(::ProcessedView; efield/db) is not allowed — the view is already in $(unit_label(getfield(v, :units))). Re-view the parent with view_units(v.parent; …) to change units.")
    src = which === :EW       ? getfield(v, :EW_amp) :
          which === :NS       ? getfield(v, :NS_amp) :
          which === :combined ? getfield(v, :combined_amp) :
          error("`which` must be :EW, :NS, or :combined (got :$which)")
    return copy(src)
end

# Refuse to persist a presentation object: only canonical pT is cacheable.
_save_entry(::AbstractString, ::ProcessedView) =
    error("Refusing to cache a ProcessedView — only the canonical linear-pT ProcessedDay is cacheable. Cache the parent (view.parent) instead.")

function Base.show(io::IO, v::ProcessedView)
    print(io, "ProcessedView(", v.rx, "←", v.tx, " ", v.date, " | ",
          unit_label(getfield(v, :units)), ", amp ",
          round(_coverage(getfield(v, :combined_amp)); digits = 1), "% covered)")
end

function Base.show(io::IO, ::MIME"text/plain", v::ProcessedView)
    u = unit_label(getfield(v, :units))
    println(io, "ProcessedView: ", v.rx, " ← ", v.tx, "   ", v.date)
    println(io, "  Fc   = ", v.Fc, " Hz    Fs = ", v.Fs, " Hz    ",
            length(v.time), " samples    units = ", u)
    println(io, "  EW_amp        (", u, ")   ", _cov(getfield(v, :EW_amp)), " covered")
    println(io, "  NS_amp        (", u, ")   ", _cov(getfield(v, :NS_amp)), " covered")
    println(io, "  combined_amp  (", u, ")   ", _cov(getfield(v, :combined_amp)), " covered")
    println(io, "  EW_pha        (deg)  ", _cov(v.EW_pha), " covered")
    println(io, "  NS_pha        (deg)  ", _cov(v.NS_pha), " covered")
    print(io,   "  parent params: ", v.params)
end

# --- quadrature combine -----------------------------------------------------
"""
    combine_quadrature(ns, ew) -> Vector{Float64}

Combined amplitude `sqrt(ns^2 + ew^2)`, elementwise. NaN in either input
propagates to NaN (a sample is only "combined" where both channels exist).
"""
combine_quadrature(ns::AbstractVector, ew::AbstractVector) = sqrt.(ns .^ 2 .+ ew .^ 2)

# --- build / get ProcessedDay ----------------------------------------------
"""
    build_processed(cache, source_folder, rx, tx, date, params;
                    baseline_ns=nothing, baseline_ew=nothing) -> Union{ProcessedDay,Nothing}

Assemble a [`ProcessedDay`](@ref) from raw entries (pulled via [`get_raw`](@ref),
so they cache on first read). Amplitude channels are optional: a missing EW or NS
amplitude (and the resulting `combined_amp`) is filled with `NaN` and a warning is
emitted, rather than erroring — so a sweep over many days is not interrupted by a
single day whose instrument never produced that channel. Phase channels are
likewise optional (NaN if missing). Returns `nothing` only when *no* raw entry of
any kind exists for the `(date, rx, tx)`.

Amplitudes are stored in linear pT (use [`amplitude`](@ref) to view in µV/m or
dB). `baseline_ns`/`baseline_ew` are optional reference-receiver phase
[`RawDay`](@ref)s subtracted before stitching; where a baseline drops out,
`params.slope` extrapolates the reference forward (see [`detrend_phase`](@ref)).
With no baseline, `params.slope` linearly detrends.
"""
function build_processed(c::VLFCache, source_folder::AbstractString, rx, tx, date::Date,
                         params::ProcessParams;
                         baseline_ns::Union{Nothing,RawDay} = nothing,
                         baseline_ew::Union{Nothing,RawDay} = nothing,
                         baseline_amp::Union{Nothing,RawDay} = nothing,
                         dropout_ranges::Union{Nothing,Vector{UnitRange{Int}}} = nothing,
                         recompute::Bool = false)
    rxs, txs = Symbol(rx), Symbol(tx)
    key(ch, q) = DataKey(date, rxs, txs, ch, q)
    raw(ch, q) = get_raw(c, source_folder, key(ch, q); recompute = recompute)

    ew_amp = raw(EW, AMPLITUDE)
    ns_amp = raw(NS, AMPLITUDE)
    ew_pha = raw(EW, PHASE)
    ns_pha = raw(NS, PHASE)

    template = ew_amp !== nothing ? ew_amp :
               ns_amp !== nothing ? ns_amp :
               ew_pha !== nothing ? ew_pha :
               ns_pha !== nothing ? ns_pha : nothing
    template === nothing && return nothing
    Fs, Fc, tg = template.Fs, template.Fc, template.time
    n = length(tg)

    if ew_amp === nothing || ns_amp === nothing
        missing_amps = String[]
        ew_amp === nothing && push!(missing_amps, "EW")
        ns_amp === nothing && push!(missing_amps, "NS")
        @warn "Missing amplitude channel(s); amplitude/combined left as NaN \
               (phase still processed where present)." rx=rxs tx=txs date missing=missing_amps
    end

    ew_amp_pT = ew_amp === nothing ? fill(NaN, n) : calibrate(ew_amp, params)
    ns_amp_pT = ns_amp === nothing ? fill(NaN, n) : calibrate(ns_amp, params)

    ew_p = ew_pha === nothing ? fill(NaN, n) : copy(ew_pha.data)
    ns_p = ns_pha === nothing ? fill(NaN, n) : copy(ns_pha.data)

    # Copy baseline phase so masking/processing never mutates the caller's (or
    # the cached) RawDay arrays.
    bew = baseline_ew === nothing ? nothing : copy(baseline_ew.data)
    bns = baseline_ns === nothing ? nothing : copy(baseline_ns.data)

    # --- transmitter-dropout masking (raw untouched) ------------------------
    # Two evidence sources combined per-sample into one mask:
    #   (a) single-receiver detection on a trusted near-field baseline amplitude
    #       (gated on params.dropout_db); fires only where the baseline records.
    #   (b) externally supplied network-coincidence `dropout_ranges`, admitted
    #       ONLY at samples where the baseline is silent (NaN), so the trusted
    #       near-field detector is never overridden where it has evidence — the
    #       network only fills the baseline's gaps.
    # A baseline whose grid disagrees with the target is treated as absent.
    local usable_baseline
    if baseline_amp === nothing
        usable_baseline = nothing
    elseif isapprox(baseline_amp.Fs, Fs; rtol = 1e-9) && length(baseline_amp.data) == n
        usable_baseline = baseline_amp
    else
        @warn "baseline_amp grid ≠ target grid — treated as absent for dropout masking." rx=rxs tx=txs date
        usable_baseline = nothing
    end

    mask = falses(n)

    if params.dropout_db !== nothing
        if usable_baseline === nothing
            @warn "dropout_db set but no usable baseline_amp — single-receiver masking skipped." rx=rxs tx=txs date
        else
            for r in detect_dropouts(usable_baseline; drop_db = params.dropout_db,
                                     window_s = params.dropout_window,
                                     pad_s = params.dropout_pad,
                                     min_valid = params.dropout_min_valid)
                @views mask[r] .= true
            end
        end
    end

    if dropout_ranges !== nothing && !isempty(dropout_ranges)
        if any(r -> r.start < 1 || r.stop > n, dropout_ranges)
            @warn "dropout_ranges fall outside target grid 1:$n — external ranges skipped." rx=rxs tx=txs date
        else
            base_silent = usable_baseline === nothing ? trues(n) : isnan.(usable_baseline.data)
            @inbounds for r in dropout_ranges, i in r
                base_silent[i] && (mask[i] = true)
            end
        end
    end

    if any(mask)
        for v in (ew_amp_pT, ns_amp_pT, ew_p, ns_p)
            v[mask] .= NaN
        end
        bew === nothing || (bew[mask] .= NaN)
        bns === nothing || (bns[mask] .= NaN)
    end

    combined_pT = combine_quadrature(ns_amp_pT, ew_amp_pT)   # NaN where either masked

    ew_stitched = _clean_detrend_stitch(ew_p, bew, params, tg)
    ns_stitched = _clean_detrend_stitch(ns_p, bns, params, tg)

    return ProcessedDay(date, rxs, txs, Fc, Fs, tg,
                        ew_amp_pT, ns_amp_pT, combined_pT,
                        ew_stitched, ns_stitched, params)
end

"""
    get_processed(cache, source_folder, rx, tx, date, params;
                  recompute=false, baseline_ns=nothing, baseline_ew=nothing) -> Union{ProcessedDay,Nothing}

Load a cached [`ProcessedDay`](@ref) or build and cache one. Returns `nothing`
when no raw data exists for the `(date, rx, tx)`, so a loop over many days can
simply skip the empty ones. Days missing an amplitude channel are built (NaN
amplitude/combined, phase populated where present) and cached like any other
product — an absent channel usually reflects a real, permanent gap in what the
instrument recorded, so there is nothing to refresh.

Unlike the raw tier, the processed product depends on `params`. With
`recompute=false` a cached product is returned even if its stored provenance
differs from `params`, but a warning is emitted telling you to pass
`recompute=true` to rebuild — so a forgotten flag after a calibration/tolerance
change is visible rather than silent.
"""
function get_processed(c::VLFCache, source_folder::AbstractString, rx, tx, date::Date,
                       params::ProcessParams; recompute::Bool = false,
                       baseline_ns::Union{Nothing,RawDay} = nothing,
                       baseline_ew::Union{Nothing,RawDay} = nothing,
                       baseline_amp::Union{Nothing,RawDay} = nothing,
                       dropout_ranges::Union{Nothing,Vector{UnitRange{Int}}} = nothing)
    if !recompute
        cached = _load_entry(processed_path(c, date, rx, tx))
        if cached isa ProcessedDay
            provenance_matches(cached.params, params) ||
                @warn "Cached ProcessedDay was built with different parameters; \
                       returning the cached product. Pass recompute=true to rebuild." rx tx date
            return cached
        end
    end
    day = build_processed(c, source_folder, rx, tx, date, params;
                          baseline_ns = baseline_ns, baseline_ew = baseline_ew,
                          baseline_amp = baseline_amp, dropout_ranges = dropout_ranges, recompute = recompute)
    day === nothing && return nothing
    _save_entry(processed_path(c, date, rx, tx), day)
    _index_add_processed!(c, date, rx, tx)
    return day
end

"""
    get_processed_view(cache, source_folder, rx, tx, date, params;
                       efield=false, db=false, recompute=false,
                       baseline_ns=nothing, baseline_ew=nothing) -> Union{ProcessedView,Nothing}

Convenience wrapper: build/read the canonical linear-pT [`ProcessedDay`](@ref) via
[`get_processed`](@ref) (which caches it), then return a [`ProcessedView`](@ref) of
it in the requested units. Units are *not* part of the cache key — only the pT
parent is cached, and the projection happens after. Returns `nothing` when no raw
data exists for the `(date, rx, tx)`.
"""
function get_processed_view(c::VLFCache, source_folder::AbstractString, rx, tx, date::Date,
                            params::ProcessParams; efield::Bool = false, db::Bool = false,
                            recompute::Bool = false,
                            baseline_ns::Union{Nothing,RawDay} = nothing,
                            baseline_ew::Union{Nothing,RawDay} = nothing,
                            baseline_amp::Union{Nothing,RawDay} = nothing,
                            dropout_ranges::Union{Nothing,Vector{UnitRange{Int}}} = nothing)
    day = get_processed(c, source_folder, rx, tx, date, params;
                        recompute = recompute, baseline_ns = baseline_ns,
                        baseline_ew = baseline_ew, baseline_amp = baseline_amp, dropout_ranges=dropout_ranges)
    return day === nothing ? nothing : view_units(day; efield = efield, db = db)
end

# Default network settings, applied under any keys the caller omits from net_kw.
# Kept here so detect_dropouts_network and network_dropout_label always receive a
# complete, identical keyword set (the label is order- and value-sensitive).
const _NET_DEFAULTS = (; drop_db = 10.0, window_s = 300.0, pad_s = 1.0, 
                         min_valid = 30, min_coincident = nothing, min_receivers = 2)

# Copy of `p` with only `dropout_label` replaced. Field-preserving, so it survives
# future ProcessParams fields without edits here.
function _with_dropout_label(p::ProcessParams, label::AbstractString)
    base = (; (f => getfield(p, f) for f in fieldnames(ProcessParams))...)
    return ProcessParams(; merge(base, (; dropout_label = label))...)
end

"""
    NetworkJob(cache, src, rx, params; detect_channel=EW)

One receiver's inputs for [`get_processed_network`](@ref): its cache, source
folder, receiver code, per-receiver [`ProcessParams`](@ref), and the channel whose
amplitude votes in the coincidence detector. Bundling them prevents the
misalignment possible with parallel vectors. `params.dropout_label` is overwritten
by the orchestrator with the network label, so leave it unset.

`detect_channel` (`EW` or `NS`) selects only the *evidence* channel for dropout
detection. The detected mask is applied to every product of the receiver — both
amplitude channels and both phase channels — independent of this choice, because a
transmitter power-down drops all channels together. The selected channel is folded
into the network label (see [`network_dropout_label`](@ref)), so a product cached
under one detection channel does not satisfy a request under another.
"""
struct NetworkJob
    cache::VLFCache
    src::String
    rx::Symbol
    params::ProcessParams
    detect_channel::Channel
end
NetworkJob(cache::VLFCache, src::AbstractString, rx, params::ProcessParams;
           detect_channel::Channel = EW) =
    NetworkJob(cache, String(src), Symbol(rx), params, detect_channel)

"""
    get_processed_network(jobs::AbstractVector{NetworkJob}, tx, date;
                          baseline_amp=nothing, baseline_ns=nothing, baseline_ew=nothing,
                          net_kw=(;), recompute=false)
        -> Vector{Union{ProcessedDay,Nothing}}

Process a set of receivers of the same `tx`/`date` as a mutual transmitter-dropout
coincidence network. Detection runs once over each job's chosen detection-channel
amplitude (`NetworkJob.detect_channel`, default `EW`)
([`detect_dropouts_network`](@ref)); the resulting ranges are masked into every
target, but — per [`build_processed`](@ref)'s gating — only at samples where the
shared near-field `baseline_amp` is itself silent (NaN), so the trusted near-field
detector is never overridden where it has evidence.

`net_kw` overrides any of the network settings (`drop_db`, `window_s`, `pad_s`,
`min_valid`, `min_coincident`, `min_receivers`); omitted keys take `_NET_DEFAULTS`.
`baseline_amp`/`baseline_ns`/`baseline_ew` are shared across all targets.

Caching: a target is rebuilt when `recompute` is set, when it is uncached, or when
its stored provenance differs from the current parameters (which include the network
`dropout_label`); otherwise the cached product is honored. The coincidence detection
is skipped entirely when every target is already cached and provenance-matched — the
raw amplitudes are still gathered (a cache hit when raw is cached) to compute the
label that gates that decision. Returns one entry per job in input order; an entry is
`nothing` when that receiver has no raw data for `(date, rx, tx)`.

!!! note
    `dropout_label` records the network identity (contributing receiver codes plus
    settings) but not the `baseline_amp` used for gating, so honoring the cache on a
    label match assumes the same baseline across calls.
"""
function get_processed_network(jobs::AbstractVector{NetworkJob}, tx, date::Date;
                               baseline_amp::Union{Nothing,RawDay} = nothing,
                               baseline_ns::Union{Nothing,RawDay} = nothing,
                               baseline_ew::Union{Nothing,RawDay} = nothing,
                               net_kw = (;), recompute::Bool = false)
    txs = Symbol(tx)
    nk  = merge(_NET_DEFAULTS, net_kw)

    amps = RawDay[]
    for j in jobs
        a = get_raw(j.cache, j.src, DataKey(date, j.rx, txs, j.detect_channel, AMPLITUDE);
                    recompute = recompute)
        a === nothing || push!(amps, a)
    end
    have_network = length(amps) >= nk.min_receivers
    net_label = have_network ? network_dropout_label(amps; nk...) : ""

    jparams = [_with_dropout_label(j.params, net_label) for j in jobs]

    cached      = Vector{Union{ProcessedDay,Nothing}}(nothing, length(jobs))
    needs_build = falses(length(jobs))
    for i in eachindex(jobs)
        if recompute
            needs_build[i] = true
        else
            c = _load_entry(processed_path(jobs[i].cache, date, jobs[i].rx, txs))
            if c isa ProcessedDay && provenance_matches(c.params, jparams[i])
                cached[i] = c
            else
                needs_build[i] = true
            end
        end
    end

    net_ranges = (any(needs_build) && have_network) ?
                 detect_dropouts_network(amps; nk...) : UnitRange{Int}[]

    results = Vector{Union{ProcessedDay,Nothing}}(nothing, length(jobs))
    for i in eachindex(jobs)
        if !needs_build[i]
            results[i] = cached[i]
            continue
        end
        day = build_processed(jobs[i].cache, jobs[i].src, jobs[i].rx, txs, date, jparams[i];
                              baseline_ns = baseline_ns, baseline_ew = baseline_ew,
                              baseline_amp = baseline_amp, dropout_ranges = net_ranges,
                              recompute = recompute)
        if day !== nothing
            _save_entry(processed_path(jobs[i].cache, date, jobs[i].rx, txs), day)
            _index_add_processed!(jobs[i].cache, date, jobs[i].rx, txs)
        end
        results[i] = day
    end
    return results
end