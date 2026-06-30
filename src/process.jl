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
                    baseline=nothing, ref_channel=nothing,
                    dropout_ranges=nothing, recompute=false) -> Union{ProcessedDay,Nothing}

Assemble a [`ProcessedDay`](@ref) from raw entries. Amplitude and phase channels are
optional (missing → `NaN` with a warning); returns `nothing` only when no raw entry
of any kind exists for the `(date, rx, tx)`. Amplitudes are stored in linear pT.

`baseline` is a pre-cleaned [`ProcessedDay`](@ref) reference (typically a near-field
receiver of the same `tx`, built with `slope` = the carrier drift, `subtract_slope =
false`, and `dropout_db = nothing`). Its phase is subtracted per channel — NS←NS,
EW←EW when `ref_channel === nothing`, or one channel onto both when `ref_channel` is
`EW`/`NS` (the single-channel near-field convention). Its amplitude (the `ref_channel`
channel, else EW) drives transmitter-dropout detection. The reference is consumed
as-is; it is not re-cleaned here. `params.baseline` is taken as provenance, not
recomputed — stamp it via [`baseline_label`](@ref)/[`_with_baseline_label`](@ref) at
the `get_*` layer (which the `get_processed*` family does).

Warnings fire on the silent-corruption modes: a reference with `subtract_slope = true`
(differential keeps a residual `slope·t`), a slope mismatch (breaks ramp cancellation
and gap extrapolation), a dropout-masked reference amplitude (under-detects target
dropouts), or a tx/date/grid mismatch.
"""
function build_processed(c::VLFCache, source_folder::AbstractString, rx, tx, date::Date,
                         params::ProcessParams;
                         baseline::Union{Nothing,ProcessedDay} = nothing,
                         ref_channel::Union{Nothing,Channel} = nothing,
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

    # --- resolve the reference against the target grid ----------------------
    # The single ProcessedDay carries every baseline role: its NS/EW phase is the
    # reference subtracted below, and its amplitude drives dropout detection. A grid
    # disagreement drops it entirely; the convention warnings flag the silent-
    # corruption modes (ramp-removed reference, mismatched slope, masked reference
    # amplitude, or a tx/date mismatch).
    local base
    if baseline === nothing
        base = nothing
    elseif !(isapprox(baseline.Fs, Fs; rtol = 1e-9) && length(baseline.time) == n)
        @warn "baseline grid ≠ target grid — baseline ignored." rx=rxs tx=txs date
        base = nothing
    else
        baseline.tx == txs ||
            @warn "baseline tx $(baseline.tx) ≠ target tx $txs — source phase will not cancel." rx=rxs date
        baseline.date == date ||
            @warn "baseline date $(baseline.date) ≠ target date $date." rx=rxs tx=txs
        baseline.params.subtract_slope &&
            @warn "baseline built with subtract_slope=true; differential carries a residual slope·t." rx=rxs tx=txs date
        baseline.params.slope == params.slope ||
            @warn "baseline slope $(baseline.params.slope) ≠ target slope $(params.slope)." rx=rxs tx=txs date
        baseline.params.dropout_db === nothing ||
            @warn "baseline built with dropout_db set; its masked amplitude under-detects target dropouts." rx=rxs tx=txs date
        base = baseline
    end

    # Detection channel of the reference amplitude (pT; the dB-relative test is
    # scale-invariant). Tracks ref_channel when it names a channel, else EW.
    det_ch   = ref_channel === NS ? NS : EW
    base_amp = base === nothing ? nothing : (det_ch === NS ? base.NS_amp : base.EW_amp)

    # --- transmitter-dropout masking ----------------------------------------
    # (a) single-receiver detection on the reference amplitude (gated on dropout_db),
    # (b) external network-coincidence ranges, admitted only where the reference is
    #     silent (NaN) so the near-field detector is never overridden where it has
    #     evidence. The reference is immutable; the target phase is masked, which
    #     makes the differential NaN there regardless of the reference phase.
    mask = falses(n)
    if params.dropout_db !== nothing
        if base_amp === nothing
            @warn "dropout_db set but no usable baseline — single-receiver masking skipped." rx=rxs tx=txs date
        else
            for r in detect_dropouts(base_amp, Fs; drop_db = params.dropout_db,
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
            base_silent = base_amp === nothing ? trues(n) : isnan.(base_amp)
            @inbounds for r in dropout_ranges, i in r
                base_silent[i] && (mask[i] = true)
            end
        end
    end
    if any(mask)
        for v in (ew_amp_pT, ns_amp_pT, ew_p, ns_p)
            v[mask] .= NaN
        end
    end

    combined_pT = combine_quadrature(ns_amp_pT, ew_amp_pT)

    # --- per-channel phase referencing --------------------------------------
    # The reference is already anchored + stitched, so it is used directly. With
    # ref_channel set, one reference channel references BOTH target channels; with
    # nothing, NS←NS and EW←EW (which needs a two-channel reference).
    ref_ew = base === nothing ? nothing : (ref_channel === NS ? base.NS_pha : base.EW_pha)
    ref_ns = base === nothing ? nothing : (ref_channel === EW ? base.EW_pha : base.NS_pha)

    if base !== nothing
        ref_ew !== nothing && all(isnan, ref_ew) &&
            @warn "reference EW phase is entirely NaN — EW target left unreferenced." rx=rxs tx=txs date
        ref_ns !== nothing && all(isnan, ref_ns) &&
            @warn "reference NS phase is entirely NaN — NS target left unreferenced." rx=rxs tx=txs date
    end

    ew_stitched = _clean_detrend_stitch(ew_p, ref_ew, params, tg)
    ns_stitched = _clean_detrend_stitch(ns_p, ref_ns, params, tg)

    return ProcessedDay(date, rxs, txs, Fc, Fs, tg,
                        ew_amp_pT, ns_amp_pT, combined_pT,
                        ew_stitched, ns_stitched, params)
end
"""
    get_processed(cache, source_folder, rx, tx, date, params;
                  recompute=false, baseline=nothing, ref_channel=nothing,
                  dropout_ranges=nothing) -> Union{ProcessedDay,Nothing}

Load a cached [`ProcessedDay`](@ref) or build and cache one. `baseline` is a
pre-cleaned [`ProcessedDay`](@ref) reference and `ref_channel` selects the referencing
mode (see [`build_processed`](@ref)); the reference's full identity is folded into
`params.baseline` via [`baseline_label`](@ref), so a target referenced against a
different baseline — or a differently-built one — is a distinct cache entry. Returns
`nothing` when no raw data exists. With `recompute=false` a cached product whose
provenance differs is returned with a warning (pass `recompute=true` to rebuild).
"""
function get_processed(c::VLFCache, source_folder::AbstractString, rx, tx, date::Date,
                       params::ProcessParams; recompute::Bool = false,
                       baseline::Union{Nothing,ProcessedDay} = nothing,
                       ref_channel::Union{Nothing,Channel} = nothing,
                       dropout_ranges::Union{Nothing,Vector{UnitRange{Int}}} = nothing)
    stamped = _with_baseline_label(params,
                  baseline === nothing ? "" : baseline_label(baseline, ref_channel))
    if !recompute
        cached = _load_entry(processed_path(c, date, rx, tx))
        if cached isa ProcessedDay
            provenance_matches(cached.params, stamped) ||
                @warn "Cached ProcessedDay was built with different parameters; \
                       returning the cached product. Pass recompute=true to rebuild." rx tx date
            return cached
        end
    end
    day = build_processed(c, source_folder, rx, tx, date, stamped;
                          baseline = baseline, ref_channel = ref_channel,
                          dropout_ranges = dropout_ranges, recompute = recompute)
    day === nothing && return nothing
    _save_entry(processed_path(c, date, rx, tx), day)
    _index_add_processed!(c, date, rx, tx)
    return day
end

"""
    get_processed_view(cache, source_folder, rx, tx, date, params;
                       efield=false, db=false, recompute=false,
                       baseline=nothing, ref_channel=nothing,
                       dropout_ranges=nothing) -> Union{ProcessedView,Nothing}

Build/read the canonical linear-pT [`ProcessedDay`](@ref) via [`get_processed`](@ref)
(which caches it), then return a [`ProcessedView`](@ref) in the requested units. Units
are not part of the cache key.
"""
function get_processed_view(c::VLFCache, source_folder::AbstractString, rx, tx, date::Date,
                            params::ProcessParams; efield::Bool = false, db::Bool = false,
                            recompute::Bool = false,
                            baseline::Union{Nothing,ProcessedDay} = nothing,
                            ref_channel::Union{Nothing,Channel} = nothing,
                            dropout_ranges::Union{Nothing,Vector{UnitRange{Int}}} = nothing)
    day = get_processed(c, source_folder, rx, tx, date, params;
                        recompute = recompute, baseline = baseline,
                        ref_channel = ref_channel, dropout_ranges = dropout_ranges)
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
    params_digest(p::ProcessParams) -> String

Deterministic, field-preserving digest of a [`ProcessParams`](@ref): every field is
stringified in declaration order, so a field added later is included automatically
without editing here. Non-cryptographic — a stable unique representation, used to
stamp a target with the full identity of its reference's processing.
"""
params_digest(p::ProcessParams) =
    join((string(f, "=", getfield(p, f)) for f in fieldnames(ProcessParams)), ";")

"""
    baseline_label(ref::ProcessedDay, ref_channel) -> String

Provenance token for a phase reference: the reference receiver/transmitter, the
referencing mode (`RR` for per-channel NS←NS/EW←EW, else the single channel forced
onto both), and a full [`params_digest`](@ref) of how the reference was built. Any
difference makes a distinct cache identity for targets referenced against it.
"""
baseline_label(ref::ProcessedDay, ref_channel::Union{Nothing,Channel}) =
    string(ref.rx, "@", ref.tx, "/",
           ref_channel === nothing ? "RR" : chstr(ref_channel),
           "/{", params_digest(ref.params), "}")

# Copy of `p` with only `baseline` replaced. Field-preserving, mirroring
# `_with_dropout_label`, so it survives future ProcessParams fields without edits.
function _with_baseline_label(p::ProcessParams, label::AbstractString)
    base = (; (f => getfield(p, f) for f in fieldnames(ProcessParams))...)
    return ProcessParams(; merge(base, (; baseline = label))...)
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
    get_processed_network(jobs, tx, date;
                          baseline=nothing, ref_channel=nothing,
                          net_kw=(;), recompute=false) -> Vector{Union{ProcessedDay,Nothing}}

Process a set of receivers of the same `tx`/`date` as a transmitter-dropout
coincidence network. The shared `baseline`/`ref_channel` reference is applied to every
target; its identity is folded into each product's provenance alongside the network
`dropout_label`, so a target cached against one reference does not satisfy a request
against another. Returns one entry per job in input order (`nothing` where a receiver
has no raw data).
"""
function get_processed_network(jobs::AbstractVector{NetworkJob}, tx, date::Date;
                               baseline::Union{Nothing,ProcessedDay} = nothing,
                               ref_channel::Union{Nothing,Channel} = nothing,
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
    blabel    = baseline === nothing ? "" : baseline_label(baseline, ref_channel)

    jparams = [_with_baseline_label(_with_dropout_label(j.params, net_label), blabel)
               for j in jobs]

    cached      = Vector{Union{ProcessedDay,Nothing}}(nothing, length(jobs))
    needs_build = falses(length(jobs))
    for i in eachindex(jobs)
        if recompute
            needs_build[i] = true
        else
            cc = _load_entry(processed_path(jobs[i].cache, date, jobs[i].rx, txs))
            if cc isa ProcessedDay && provenance_matches(cc.params, jparams[i])
                cached[i] = cc
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
                              baseline = baseline, ref_channel = ref_channel,
                              dropout_ranges = net_ranges, recompute = recompute)
        if day !== nothing
            _save_entry(processed_path(jobs[i].cache, date, jobs[i].rx, txs), day)
            _index_add_processed!(jobs[i].cache, date, jobs[i].rx, txs)
        end
        results[i] = day
    end
    return results
end