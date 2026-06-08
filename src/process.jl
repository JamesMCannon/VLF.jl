# ============================================================================
# process.jl  —  Tier 2: raw counts -> calibrated/cleaned products.
#
# Everything here operates on the canonical full-day grid, so the legacy
# time-matching / deleteat! gymnastics are gone: gaps are explicit NaNs and the
# two channels are aligned by construction.
# ============================================================================

# --- calibration ------------------------------------------------------------
"""
    calibration_factor(amp_day::RawDay, params) -> Float64

Counts→pT scale factor. Uses `params.cal_num` if set (≠ -1.0); otherwise reads
`params.cal_file`, selects the NS/EW curve by `amp_day.channel`, and takes the
response at the frequency nearest `amp_day.Fc` (cal-file frequencies are in kHz).
"""
function calibration_factor(amp_day::RawDay, params::ProcessParams)
    if params.cal_num != -1.0
        return params.cal_num
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

# --- build / get ProcessedDay ----------------------------------------------
"""
    build_processed(cache, source_folder, rx, tx, date, params;
                    baseline_ns=nothing, baseline_ew=nothing) -> ProcessedDay

Assemble a [`ProcessedDay`](@ref) from raw entries (pulled via [`get_raw`](@ref),
so they cache on first read). Requires both amplitude channels; phase channels
are optional (filled with NaN if missing). `baseline_ns`/`baseline_ew` are
optional reference-receiver phase [`RawDay`](@ref)s to subtract before stitching.
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
    (ew_amp === nothing || ns_amp === nothing) &&
        error("Missing amplitude data for $rxs←$txs on $date (need both EW and NS).")

    Fs, Fc, tg = ew_amp.Fs, ew_amp.Fc, ew_amp.time
    ew_amp_pT = calibrate(ew_amp, params)
    ns_amp_pT = calibrate(ns_amp, params)
    combined  = combine_quadrature(ns_amp_pT, ew_amp_pT)

    ew_pha = raw(EW, PHASE)
    ns_pha = raw(NS, PHASE)
    ew_p = ew_pha === nothing ? fill(NaN, length(tg)) : copy(ew_pha.data)
    ns_p = ns_pha === nothing ? fill(NaN, length(tg)) : copy(ns_pha.data)

    baseline_ew !== nothing && (ew_p = clean_phase(ew_p, baseline_ew.data))
    baseline_ns !== nothing && (ns_p = clean_phase(ns_p, baseline_ns.data))

    ew_stitched = stitch_phase(ew_p; tolerance = params.tolerance, unwrap = params.unwrap)
    ns_stitched = stitch_phase(ns_p; tolerance = params.tolerance, unwrap = params.unwrap)

    return ProcessedDay(date, rxs, txs, Fc, Fs, tg,
                        ew_amp_pT, ns_amp_pT, combined,
                        ew_stitched, ns_stitched, params)
end

"""
    get_processed(cache, source_folder, rx, tx, date, params;
                  recompute=false, baseline_ns=nothing, baseline_ew=nothing) -> ProcessedDay

Load a cached [`ProcessedDay`](@ref) or build and cache one.

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
    _save_entry(processed_path(c, date, rx, tx), day)
    _index_add_processed!(c, date, rx, tx)
    return day
end