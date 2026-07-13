"""
    rotate_day(day::ProcessedDay, bearing_deg; polarity=(NS=1, EW=-1), offset_m=0, baseline_label="") -> RotatedDay

Rotate `day`'s NS/EW horizontal phasors into the radial/azimuthal frame of Gross
et al. (2018), Eq. (5). `bearing_deg` is `θ_az` (rx→tx, degrees clockwise from
true north). Two independent per-receiver corrections act on the EW phasor before
the mix: `polarity` is the antenna wiring sign (`(NS, EW)`, each ±1; the Gross
convention is `EW = -1`), and `offset_m ∈ {0,1,2,3}` is the demodulation
quarter-turn `u_rel` (`+90·offset_m°` on EW — `0` for synchronous demodulation,
resolved from the γ ladder for asynchronous). 
Both corrections assume the parent's phases are free of inter-channel sampling
skew (`ProcessParams.ew_dt_s`, removed at build); an uncorrected multiplexed-DAQ
parent leaves a continuous EW phase error that no quarter-turn can absorb,
producing an irreducible B_r/B_azi leakage floor and biased rotated phases. Reconstructs
`v_NS = polarity.NS · A_NS · e^{jψ_NS}` and
`v_EW = polarity.EW · A_EW · e^{j(ψ_EW + 90·offset_m°)}` from the calibrated pT
amplitudes and cleaned phases, then `[B_r; B_azi] = R(θ_az + π/2) · [v_NS; v_EW]`
with `R(α) = [cos α  sin α; −sin α  cos α]`.

Phases are absolute unless `day` was built with a phase baseline (a differential)
or with the slope actually detrended (`subtract_slope = true`); a slope used only
to anchor the unwrap leaves the rotated phase absolute, since whole-turn folds
vanish under rotation. Amplitudes and the ratio ∠(−B_r/B_azi) are unaffected
either way.
"""
function rotate_day(day::ProcessedDay, bearing_deg::Real;
                    polarity::NamedTuple = (NS = 1, EW = -1),
                    offset_m::Integer = 0,
                    baseline_label::AbstractString = "")
    detrended = day.params.slope !== nothing && day.params.subtract_slope
    if !isempty(day.params.baseline) || detrended
        @warn "rotate_day: parent phase is referenced/detrended (baseline=\"$(day.params.baseline)\", \
               slope=$(day.params.slope), subtract_slope=$(day.params.subtract_slope)); \
               rotated phases inherit it and are not absolute." day.rx day.tx day.date
    end

    α      = deg2rad(bearing_deg) + π/2
    cα, sα = cos(α), sin(α)
    pNS, pEW = polarity.NS, polarity.EW
    dϕEW   = 90.0 * offset_m

    n = length(day.time)
    Br_amp   = Vector{Float64}(undef, n);  Bazi_amp = Vector{Float64}(undef, n)
    Br_pha   = Vector{Float64}(undef, n);  Bazi_pha = Vector{Float64}(undef, n)

    @inbounds for i in 1:n
        vNS = pNS * day.NS_amp[i] * cis(deg2rad(day.NS_pha[i]))
        vEW = pEW * day.EW_amp[i] * cis(deg2rad(day.EW_pha[i] + dϕEW))
        Br   =  cα * vNS + sα * vEW
        Bazi = -sα * vNS + cα * vEW
        Br_amp[i]   = abs(Br);   Br_pha[i]   = rad2deg(angle(Br))
        Bazi_amp[i] = abs(Bazi); Bazi_pha[i] = rad2deg(angle(Bazi))
    end

    return RotatedDay(day.date, day.rx, day.tx, day.Fc, day.Fs, day.time,
                      float(bearing_deg), Br_amp, Bazi_amp, Br_pha, Bazi_pha,
                      polarity, Int(offset_m), day.params, String(baseline_label))
end

# Wrap to (-180, 180], the RotatedDay phase convention (matches rad2deg∘angle).
_wrap180(x) = isnan(x) ? x : rem(x, 360, RoundNearest)

"""
    baseline_subtract(target::RotatedDay, reference::RotatedDay; label=string(reference.rx)) -> RotatedDay

Per-component path differential against a MEASURED reference: subtract `reference`'s
rotated phase from `target`'s, component by component, cancelling the common-mode
transmitter source phase (offset, drift, and wander). Amplitudes, `polarity`, and
`offset_m` are the target's — the result carries the target's rotation convention,
not the reference's. `target` and `reference` must share `tx`, `date`, and grid;
`bearing_deg` stays the target's. Result phases are wrapped to (-180, 180].
"""
function baseline_subtract(target::RotatedDay, reference::RotatedDay;
                           label::AbstractString = string(reference.rx))
    (target.tx == reference.tx && target.date == reference.date) ||
        error("baseline_subtract: tx/date mismatch \
               (target $(target.tx)/$(target.date), reference $(reference.tx)/$(reference.date)).")
    (isapprox(target.Fs, reference.Fs; rtol = 1e-9) &&
     length(target.time) == length(reference.time)) ||
        error("baseline_subtract: target and reference grids differ.")

    Br_pha   = _wrap180.(target.Br_pha   .- reference.Br_pha)
    Bazi_pha = _wrap180.(target.Bazi_pha .- reference.Bazi_pha)

    return RotatedDay(target.date, target.rx, target.tx, target.Fc, target.Fs, target.time,
                      target.bearing_deg, copy(target.Br_amp), copy(target.Bazi_amp),
                      Br_pha, Bazi_pha, target.polarity, target.offset_m, target.params,
                      String(label))
end

"""
    baseline_subtract(target::RotatedDay, slope::Real; label="slope=\$slope") -> RotatedDay

Modeled-linear fallback for a transmitter with no measured near-field reference.
Removes the nominal carrier ramp `slope · t` (degrees; `slope` in deg/s, `t` in
seconds-from-midnight) from each rotated component and wraps to (-180, 180].
`polarity` and `offset_m` are the target's, carried through unchanged since this
method only removes a phase ramp. `slope` MUST be the slope the target was
processed with under `subtract_slope = false`, so the ramp it carries is the one
removed.

This removes only the DETERMINISTIC linear drift. Unlike a measured reference it
leaves the constant source-phase offset and any nonlinear clock wander, so the
result is a detrended product for relative phase structure, NOT a source-referenced
path differential. For absolute-phase work use a near-field reference or pin
against a model (LMP).
"""
function baseline_subtract(target::RotatedDay, slope::Real;
                           label::AbstractString = string("slope=", slope))
    ramp     = slope .* target.time
    Br_pha   = _wrap180.(target.Br_pha   .- ramp)
    Bazi_pha = _wrap180.(target.Bazi_pha .- ramp)
    return RotatedDay(target.date, target.rx, target.tx, target.Fc, target.Fs, target.time,
                      target.bearing_deg, copy(target.Br_amp), copy(target.Bazi_amp),
                      Br_pha, Bazi_pha, target.polarity, target.offset_m, target.params,
                      String(label))
end