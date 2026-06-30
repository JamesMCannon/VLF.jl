# ============================================================================
# dropouts.jl  —  Tier 2: transmitter-dropout detection and masking.
#
# Two detectors over the canonical grid, sharing `_flag_drops` /
# `_ranges_from_flags`:
#   * single-receiver    — a relative trailing-median drop on a near-field
#                          baseline amplitude (`detect_dropouts`);
#   * network-coincidence — the same per-receiver test combined across receivers
#                          of one transmitter (`detect_dropouts_network`), with a
#                          deterministic provenance label (`network_dropout_label`).
# Ranges are applied to a RawDay by `mask_dropouts`. The masking is integrated and
# gated per-sample inside `build_processed`, and the network is orchestrated by
# `get_processed_network` (both in process.jl).
# ============================================================================

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
    detect_dropouts(amp::AbstractVector, Fs::Real; drop_db=10.0, window_s=300.0,
                    pad_s=1.0, min_valid=30) -> Vector{UnitRange{Int}}

Trailing-relative-drop detection on an amplitude vector with explicit `Fs`. The
test is dB below a trailing median, so units are irrelevant — calibrated pT and raw
counts give identical ranges. See the [`RawDay`](@ref) method for the gated form.
"""
function detect_dropouts(amp::AbstractVector, Fs::Real; drop_db::Real = 10.0,
                         window_s::Real = 300.0, pad_s::Real = 1.0,
                         min_valid::Integer = 30)
    adb = to_db.(amp)
    flagged, _ = _flag_drops(adb, Fs; drop_db = drop_db, window_s = window_s,
                             min_valid = min_valid)
    return _ranges_from_flags(flagged, round(Int, pad_s * Fs))
end

"""
    detect_dropouts(baseline_amp::RawDay; kwargs...) -> Vector{UnitRange{Int}}

Trailing-relative-drop detection on a single receiver's amplitude. Delegates to the
vector form on `baseline_amp.data`/`baseline_amp.Fs`.
"""
detect_dropouts(baseline_amp::RawDay; kwargs...) =
    detect_dropouts(baseline_amp.data, baseline_amp.Fs; kwargs...)

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
sorted contributing receiver codes, each tagged with the channel that supplied its
evidence (`FSI:EW`), plus the coincidence settings. Pass the result as
`ProcessParams(...; dropout_label=…)` alongside the same network's ranges so the
cache can tell network-masked products apart. Derive the label and the ranges from
the *same* `amps` and keyword values to keep them consistent.
"""
function network_dropout_label(amps::AbstractVector{RawDay};
                               drop_db, window_s, pad_s, min_valid,
                               min_coincident, min_receivers)
    rxs = join(sort([string(a.rx) * ":" * chstr(a.channel) for a in amps]), ",")
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