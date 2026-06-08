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

# --- phase: unwrap ----------------------------------------------------------
"""
    unwrap_phase(phase) -> Vector{Float64}

Remove ±360° wraps so successive samples differ by < 180°. Operates on a copy
(the old in-place version returned its mutated input). NaN samples are left
untouched and do not corrupt the running value.
"""
function unwrap_phase(phase::AbstractVector)
    out = collect(float.(phase))
    @inbounds for i in 2:length(out)
        d = out[i] - out[i-1]
        isnan(d) && continue
        if d > 180
            out[i] -= 360 * ceil((d - 180) / 360)
        elseif d < -180
            out[i] += 360 * ceil((-d - 180) / 360)
        end
    end
    return out
end

# --- phase: n×90° stitch across gaps ----------------------------------------
# Generalizes the original's explicit ±90/180/270 branches (its own TODO asked
# for "a general expression to check n×90 degrees").
function _quarter_correction(del, tol)
    @inbounds for (mult, corr) in ((1, -90), (2, -180), (3, -270),
                                   (-1, 90), (-2, 180), (-3, 270))
        t = 90 * mult
        (t - tol) < del < (t + tol) && return corr
    end
    return 0
end

"""
    stitch_phase(phase; tolerance=10, unwrap=true) -> Vector{Float64}

Correct integer-quarter-cycle (n×90°) jumps introduced across data gaps. In the
gridded model a gap is a run of NaNs, so the jump is measured between the last
valid sample before the gap and the first valid sample after it; if that
difference is within `tolerance` of n×90°, the remainder is shifted to cancel
it. Optionally unwraps first.
"""
function stitch_phase(phase::AbstractVector; tolerance::Real = 10, unwrap::Bool = true)
    out = unwrap ? unwrap_phase(phase) : collect(float.(phase))
    last_valid = nothing       # value of most recent non-NaN sample
    in_gap = false
    @inbounds for i in eachindex(out)
        if isnan(out[i])
            last_valid !== nothing && (in_gap = true)
        else
            if in_gap && last_valid !== nothing
                corr = _quarter_correction(out[i] - last_valid, tolerance)
                corr != 0 && (out[i:end] .+= corr)
            end
            in_gap = false
            last_valid = out[i]
        end
    end
    return out
end

# --- phase: baseline cleaning ----------------------------------------------
"""
    clean_phase(day_phase, baseline_phase) -> Vector{Float64}

Subtract a reference receiver's phase from `day_phase`, elementwise. Both must
be on the same grid (guaranteed for [`RawDay`](@ref)s of equal `Fs`), so this
replaces the old time-alignment loop with a plain subtraction. NaN-aware.
"""
clean_phase(day_phase::AbstractVector, baseline_phase::AbstractVector) =
    day_phase .- baseline_phase

# --- phase: linear-slope detrend (with optional reference baseline) ---------
"""
    detrend_phase(target, baseline, slope, time) -> Vector{Float64}

Resolve a `target` phase against an optional reference `baseline` phase and an
optional linear `slope` (deg/s) over `time` (seconds-from-midnight), elementwise:

- **No baseline** (`baseline === nothing`): subtract `slope · t` from every
  sample. `slope === nothing` is treated as no detrend — the phase is returned
  unchanged (the original behavior).
- **Baseline present**: where the baseline is valid, subtract it
  (`target − baseline`). Where the baseline is `NaN` but the target is valid, the
  reference is *extrapolated forward from its last valid value* at `slope`, i.e.
  `target − (R₀ + slope · (t − t₀))` where `(t₀, R₀)` is the most recent valid
  reference sample. Because the synthesized reference equals `R₀` at `t₀`, the
  result is continuous across the entry into the gap (no step). If no reference
  has been seen yet (a leading gap), it falls back to `target − slope · t`.
  When `slope === nothing`, gaps are `NaN` (the original behavior).

A `NaN` in the target always yields `NaN`. The inputs are not mutated.

!!! note
    Continuity is guaranteed entering a gap. When the real reference resumes the
    output returns to `target − baseline`, which may step if the extrapolation
    drifted away from the resumed reference — real data is trusted over the
    extrapolation there.
"""
function detrend_phase(target::AbstractVector,
                       baseline::Union{Nothing,AbstractVector},
                       slope::Union{Nothing,Real},
                       time::AbstractVector)
    n = length(target)
    baseline === nothing || length(baseline) == n ||
        throw(DimensionMismatch("baseline length $(length(baseline)) ≠ target length $n"))
    length(time) == n ||
        throw(DimensionMismatch("time length $(length(time)) ≠ target length $n"))

    out = Vector{Float64}(undef, n)
    s = slope === nothing ? nothing : float(slope)

    if baseline === nothing
        @inbounds for i in 1:n
            t = float(target[i])
            out[i] = s === nothing ? t : t - s * time[i]
        end
        return out
    end

    have_anchor = false
    t0 = 0.0          # time of the last valid reference sample
    r0 = 0.0          # value of the last valid reference sample
    @inbounds for i in 1:n
        t = float(target[i])
        b = float(baseline[i])
        if !isnan(b)
            out[i] = t - b                       # NaN if target NaN; else cleaned value
            have_anchor = true; t0 = time[i]; r0 = b
        elseif isnan(t) || s === nothing
            out[i] = NaN                         # target gap, or no slope -> drop
        elseif have_anchor
            out[i] = t - (r0 + s * (time[i] - t0))   # extrapolate reference from anchor
        else
            out[i] = t - s * time[i]             # leading gap: detrend from origin
        end
    end
    return out
end

# Full phase pipeline: (optional) unwrap → baseline/slope resolve → n×90° stitch.
# Unwrapping happens here so the slope detrends the UNWRAPPED phase; stitch then
# runs with unwrap already done. With slope === nothing and no baseline this
# reduces exactly to the previous `stitch_phase(target; unwrap=params.unwrap)`.
function _clean_detrend_stitch(target::AbstractVector,
                               baseline::Union{Nothing,AbstractVector},
                               params::ProcessParams, time::AbstractVector)
    tgt  = params.unwrap ? unwrap_phase(target) : collect(float.(target))
    base = baseline === nothing ? nothing :
           (params.unwrap ? unwrap_phase(baseline) : baseline)
    cleaned = detrend_phase(tgt, base, params.slope, time)
    return stitch_phase(cleaned; tolerance = params.tolerance, unwrap = false)
end

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
                         recompute::Bool = false)
    rxs, txs = Symbol(rx), Symbol(tx)
    key(ch, q) = DataKey(date, rxs, txs, ch, q)
    raw(ch, q) = get_raw(c, source_folder, key(ch, q); recompute = recompute)

    ew_amp = raw(EW, AMPLITUDE)
    ns_amp = raw(NS, AMPLITUDE)
    ew_pha = raw(EW, PHASE)
    ns_pha = raw(NS, PHASE)

    # Derive the canonical grid from whichever entry exists (Fs/Fc/time match
    # across channels of one tx-day). Nothing at all -> there is no data here.
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

    # Calibration only runs for channels that are present, so a present amplitude
    # still requires calibration info while an absent one is simply NaN-filled.
    ew_amp_pT   = ew_amp === nothing ? fill(NaN, n) : calibrate(ew_amp, params)
    ns_amp_pT   = ns_amp === nothing ? fill(NaN, n) : calibrate(ns_amp, params)
    combined_pT = combine_quadrature(ns_amp_pT, ew_amp_pT)   # NaN unless both exist

    ew_p = ew_pha === nothing ? fill(NaN, n) : copy(ew_pha.data)
    ns_p = ns_pha === nothing ? fill(NaN, n) : copy(ns_pha.data)

    bew = baseline_ew === nothing ? nothing : baseline_ew.data
    bns = baseline_ns === nothing ? nothing : baseline_ns.data

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
                       baseline_ew::Union{Nothing,RawDay} = nothing)
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
                          recompute = recompute)
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
                            baseline_ew::Union{Nothing,RawDay} = nothing)
    day = get_processed(c, source_folder, rx, tx, date, params;
                        recompute = recompute, baseline_ns = baseline_ns, baseline_ew = baseline_ew)
    return day === nothing ? nothing : view_units(day; efield = efield, db = db)
end