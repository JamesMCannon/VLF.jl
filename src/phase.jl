# ============================================================================
# phase.jl  —  Tier 2: phase cleaning pipeline.
#
# Operates on the canonical full-day grid (gaps are explicit NaNs), so the target,
# the opposite channel, and any reference baseline are aligned by construction.
# The pipeline is (optional) unwrap -> baseline/slope resolve -> n*90 deg gap
# stitch, assembled by `_clean_detrend_stitch` and consumed by `build_processed`
# (process.jl).
# ============================================================================

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