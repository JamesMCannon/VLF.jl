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
function unwrap_phase(phase::AbstractVector;
                      slope::Union{Nothing,Real} = nothing,
                      time::Union{Nothing,AbstractVector} = nothing,
                      max_gap::Real = Inf)
    out = collect(float.(phase))

    # Original datum-resetting behavior, preserved exactly when no anchor given.
    if slope === nothing
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

    # Anchored: fold each valid sample toward the drift-extrapolated prediction
    # from the last valid sample, so a gap no longer resets the datum. Gaps
    # longer than `max_gap` are NOT folded (extrapolation untrustworthy) — the
    # datum resets there, recovering the old behavior for long dropouts.
    time === nothing && throw(ArgumentError("anchored unwrap (slope set) needs `time`"))
    length(time) == length(out) ||
        throw(DimensionMismatch("time length $(length(time)) ≠ phase length $(length(out))"))
    s = float(slope)
    prev_v = nothing; prev_t = 0.0
    @inbounds for i in eachindex(out)
        isnan(out[i]) && continue
        if prev_v === nothing || (float(time[i]) - prev_t) > max_gap
            prev_v = out[i]; prev_t = float(time[i]); continue   # datum origin / reset
        end
        pred = prev_v + s * (float(time[i]) - prev_t)
        out[i] += 360 * round((pred - out[i]) / 360)
        prev_v = out[i]; prev_t = float(time[i])
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
function _clean_detrend_stitch(target, baseline, params, time)
    uw(v) = params.unwrap ? unwrap_phase(v; slope = params.slope, time = time,
                                          max_gap = params.max_gap) :
                            collect(float.(v))
    tgt  = uw(target)
    base = baseline === nothing ? nothing : uw(baseline)

    # Stitch the baseline independently before subtraction. This catches quarter
    # artifacts at baseline-only gap boundaries, which survive into cleaned as
    # steps between valid samples (no NaN in cleaned → post-subtraction stitch
    # never fires there). Same tolerance as the target pass.
    base_s = base === nothing ? nothing :
             stitch_phase(base; tolerance = params.tolerance, unwrap = false)

    cleaned = detrend_phase(tgt, base_s, params.slope, time)
    return stitch_phase(cleaned; tolerance = params.tolerance, unwrap = false)
end

# NaN-aware median of a small buffer (sorts a copy). Order-preserving, so
# median(dB) == dB(median): the dB/linear choice is irrelevant to the reference.
function _nanmedian(vals::AbstractVector{<:Real})
    isempty(vals) && return NaN
    s = sort(vals); n = length(s)
    return isodd(n) ? s[(n+1)÷2] : 0.5 * (s[n÷2] + s[n÷2 + 1])
end

# Per-sample drop flagging, factored out so the single-receiver and network
# detectors share one definition. Returns (flagged, eligible):
#   eligible[i] = sample i is valid AND its trailing window has ≥ min_valid valids
#                 (i.e. i is allowed to "vote"); a NaN sample or a thin window abstains.
#   flagged[i]  = eligible[i] AND adb[i] is > drop_db below the trailing-window median.
function _flag_drops(adb::AbstractVector, Fs::Real;
                     drop_db::Real, window_s::Real, min_valid::Integer)
    n = length(adb)
    w = max(1, round(Int, window_s * Fs))
    flagged  = falses(n)
    eligible = falses(n)
    buf = Float64[]
    @inbounds for i in 1:n
        isnan(adb[i]) && continue
        empty!(buf)
        for j in max(1, i - w + 1):i
            x = adb[j]; isnan(x) || push!(buf, x)
        end
        length(buf) < min_valid && continue
        eligible[i] = true
        adb[i] < _nanmedian(buf) - drop_db && (flagged[i] = true)
    end
    return flagged, eligible
end

# Contiguous flagged runs → dilated, merged 1-based index ranges (pad in samples).
function _ranges_from_flags(flagged::AbstractVector{Bool}, pad::Integer)
    n = length(flagged)
    ranges = UnitRange{Int}[]
    i = 1
    while i <= n
        if flagged[i]
            j = i; while j < n && flagged[j + 1]; j += 1; end
            lo = max(1, i - pad); hi = min(n, j + pad)
            if !isempty(ranges) && lo <= ranges[end].stop + 1
                ranges[end] = ranges[end].start:max(ranges[end].stop, hi)
            else
                push!(ranges, lo:hi)
            end
            i = j + 1
        else
            i += 1
        end
    end
    return ranges
end

"""
    detect_dropouts(baseline_amp::RawDay; drop_db=10.0, window_s=300.0,
                    pad_s=1.0, min_valid=30) -> Vector{UnitRange{Int}}

Single-receiver transmitter-dropout detection on a near-field baseline amplitude
(see module notes). Unchanged in behavior; now expressed over [`_flag_drops`](@ref)
and [`_ranges_from_flags`](@ref), which it shares with [`detect_dropouts_network`](@ref).
"""
function detect_dropouts(baseline_amp::RawDay; drop_db::Real = 10.0,
                         window_s::Real = 300.0, pad_s::Real = 1.0,
                         min_valid::Integer = 30)
    adb = to_db.(baseline_amp.data)
    flagged, _ = _flag_drops(adb, baseline_amp.Fs;
                             drop_db = drop_db, window_s = window_s,
                             min_valid = min_valid)
    return _ranges_from_flags(flagged, round(Int, pad_s * baseline_amp.Fs))
end

"""
    detect_dropouts_network(days::AbstractVector{RawDay}; drop_db=10.0,
                            window_s=300.0, pad_s=1.0, min_valid=30,
                            min_coincident=nothing, min_receivers=2)
        -> Vector{UnitRange{Int}}

Detect transmitter power-dropouts by *coincidence* across a network of receivers
that all observe the same transmitter on the same date. Each receiver's amplitude
is flagged independently with the same trailing-relative-drop test as
[`detect_dropouts`](@ref); a sample is then declared a network dropout when enough
*eligible* receivers (those recording with a sufficient trailing window at that
sample — see [`_flag_drops`](@ref)) flag it simultaneously. Flagged runs are
dilated by `pad_s` and merged. The returned ranges are on the shared grid and
mask every co-registered product of every receiver.

Coincidence rule:
- `min_coincident === nothing` (default): **unanimous** — every eligible receiver
  at that sample must flag, and at least `min_receivers` must be eligible.
  Receivers that are not recording (NaN) or whose window is too thin are excluded
  from the denominator, so a down receiver neither blocks nor triggers detection.
- `min_coincident::Integer`: require at least that many receivers to flag,
  regardless of how many are eligible (`min_receivers` is then unused).

All receivers must share `Fs` and grid length; a mismatch is an error rather than
a silent drop, because dropping a receiver would change the coincidence denominator
without warning. A non-`AMPLITUDE` input warns but is not rejected (the test is
relative per receiver).

!!! warning "Epistemic status"
    This is a *backstop* for when the near-field baseline is unavailable. Unlike
    the ground-wave baseline, far-field amplitudes are subject to ionospheric
    fades, and a solar flare (SID) can drive coherent drops across many paths of a
    transmitter near-simultaneously — so a coincident drop is *consistent with* but
    not *proof of* a transmitter power-down. Discriminators not exploited here:
    depth (a power-down reaches the noise floor, typically ≫ 10 dB), onset sharpness
    (a switch is a sub-second step; an SID ramps over tens of seconds), and cross-TX
    independence (a power-down spares other transmitters' paths; an SID does not).
    Consider a deeper `drop_db` than the near-field default to suppress flare/fade
    false positives. False positives delete a real simultaneous sample across the
    whole network, so bias conservative.
"""
function detect_dropouts_network(days::AbstractVector{RawDay};
                                 drop_db::Real = 10.0, window_s::Real = 300.0,
                                 pad_s::Real = 1.0, min_valid::Integer = 30,
                                 min_coincident::Union{Nothing,Integer} = nothing,
                                 min_receivers::Integer = 2)
    isempty(days) && return UnitRange{Int}[]
    min_coincident === nothing || min_coincident >= 1 ||
        throw(ArgumentError("min_coincident must be ≥ 1 (got $min_coincident)"))
    mr = max(1, min_receivers)

    Fs = days[1].Fs
    n  = length(days[1].data)
    for d in days
        d.quantity == AMPLITUDE ||
            @warn "detect_dropouts_network: expected AMPLITUDE." rx=d.rx quantity=d.quantity
        (isapprox(d.Fs, Fs; rtol = 1e-9) && length(d.data) == n) ||
            error("detect_dropouts_network: receiver $(d.rx) grid (Fs=$(d.Fs), \
                   n=$(length(d.data))) ≠ reference (Fs=$Fs, n=$n). All receivers \
                   must share the canonical grid.")
    end

    R = length(days)
    flags = Vector{BitVector}(undef, R)
    elig  = Vector{BitVector}(undef, R)
    for r in 1:R
        flags[r], elig[r] = _flag_drops(to_db.(days[r].data), Fs;
                                        drop_db = drop_db, window_s = window_s,
                                        min_valid = min_valid)
    end

    coincident = falses(n)
    @inbounds for i in 1:n
        e = 0; f = 0
        for r in 1:R
            elig[r][i]  && (e += 1)
            flags[r][i] && (f += 1)
        end
        coincident[i] = min_coincident === nothing ? (e >= mr && f == e) :
                                                      (f >= min_coincident)
    end
    return _ranges_from_flags(coincident, round(Int, pad_s * Fs))
end

"""
    network_dropout_label(amps; drop_db, window_s, pad_s, min_valid,
                          min_coincident, min_receivers) -> String

Deterministic provenance label for a [`detect_dropouts_network`](@ref) mask: the
sorted contributing receiver codes plus the coincidence settings. Pass the result
as `ProcessParams(...; dropout_label=…)` alongside the same network's ranges so the
cache can tell network-masked products apart. Derive the label and the ranges from
the *same* `amps` and keyword values to keep them consistent.
"""
function network_dropout_label(amps::AbstractVector{RawDay};
                               drop_db, window_s, pad_s, min_valid,
                               min_coincident, min_receivers)
    rxs = join(sort([string(a.rx) for a in amps]), ",")
    mc  = min_coincident === nothing ? "all" : string(min_coincident)
    return "network[$rxs];db=$drop_db,win=$window_s,pad=$pad_s," *
           "mv=$min_valid,mc=$mc,mr=$min_receivers"
end

"""
    mask_dropouts_network(days; kwargs...) -> (; ranges, days)

Run [`detect_dropouts_network`](@ref) over `days`, then return the detected
`ranges` together with a vector of the input days each NaN-masked over those
ranges (via [`mask_dropouts`](@ref); inputs are not mutated). `kwargs` are passed
through to the detector. To also notch phase (or the other channel), reuse
`ranges`: `mask_dropouts(phase_rawday, ranges)`.
"""
function mask_dropouts_network(days::AbstractVector{RawDay}; kwargs...)
    ranges = detect_dropouts_network(days; kwargs...)
    return (ranges = ranges, days = [mask_dropouts(d, ranges) for d in days])
end

"""
    mask_dropouts(day::RawDay, ranges) -> RawDay

Copy of `day` with `data[r] .= NaN` for each range in `ranges`. Does not mutate
`day` — raw stays faithful to the source `.mat`; the NaNs live only in the copy
fed forward into processing.
"""
function mask_dropouts(day::RawDay, ranges)
    d = copy(day.data)
    for r in ranges
        @views d[r] .= NaN
    end
    return RawDay(day.date, day.rx, day.tx, day.channel, day.quantity,
                  day.Fc, day.Fs, day.time, d)
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
    NetworkJob(cache, src, rx, params)

One receiver's inputs for [`get_processed_network`](@ref): its cache, source
folder, receiver code, and per-receiver [`ProcessParams`](@ref). Bundling them
prevents the misalignment possible with parallel vectors. `params.dropout_label`
is overwritten by the orchestrator with the network label, so leave it unset.
"""
struct NetworkJob
    cache::VLFCache
    src::String
    rx::Symbol
    params::ProcessParams
end
NetworkJob(cache::VLFCache, src::AbstractString, rx, params::ProcessParams) =
    NetworkJob(cache, String(src), Symbol(rx), params)

"""
    get_processed_network(jobs::AbstractVector{NetworkJob}, tx, date;
                          baseline_amp=nothing, baseline_ns=nothing, baseline_ew=nothing,
                          net_kw=(;), recompute=false)
        -> Vector{Union{ProcessedDay,Nothing}}

Process a set of receivers of the same `tx`/`date` as a mutual transmitter-dropout
coincidence network. Detection runs once over the jobs' raw EW amplitudes
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
        a = get_raw(j.cache, j.src, DataKey(date, j.rx, txs, EW, AMPLITUDE);
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