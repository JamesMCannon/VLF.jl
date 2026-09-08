"""
    rotate_day(day::ProcessedDay, bearing_deg; polarity=(NS=1, EW=-1), offset_m=0, baseline_label="") -> RotatedDay

Rotate `day`'s NS/EW horizontal phasors into the radial/azimuthal frame of Gross
et al. (2018), Eq. (5). `bearing_deg` is `θ_az` (rx→tx, degrees clockwise from
true north). Two independent per-receiver corrections act on the EW phasor before
the mix: `polarity` is the antenna wiring sign (`(NS, EW)`, each ±1; the Gross
convention is `EW = -1`), and `offset_m ∈ {0,1,2,3}` is the demodulation
quarter-turn `u_rel` (`+90·offset_m°` on EW — `0` for synchronous demodulation,
resolved from the γ ladder for asynchronous). 
`rho_deg`/`xi_deg` invert the Wood & Inan (2004) antenna-frame geometry before
the bearing rotation: the measured phasors satisfy `m_NS ∝ cos(θ−ρ)`,
`m_EW ∝ sin(θ−ρ−ξ)` (all azimuths degrees clockwise from true north), i.e.
`m = M·[B_N; B_E]` with `det M = cos ξ`, and the applied map is
`R(θ_az + 90°)·M⁻¹`. Gains are assumed already removed by calibration; only
geometry is inverted here. At ξ = 0 the correction is exactly a bearing shift
`θ_az → θ_az − ρ`; |ξ| ≥ 90° errors (degenerate pair). For a parent built with
`swap_channels`, the slots already carry physical orientation, so ρ/ξ describe
the physical antennas directly.
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
                    rho_deg::Real = 0.0,
                    xi_deg::Real = 0.0,
                    baseline_label::AbstractString = "")
    detrended = day.params.slope !== nothing && day.params.subtract_slope
    if !isempty(day.params.baseline) || detrended
        @warn "rotate_day: parent phase is referenced/detrended (baseline=\"$(day.params.baseline)\", \
               slope=$(day.params.slope), subtract_slope=$(day.params.subtract_slope)); \
               rotated phases inherit it and are not absolute." day.rx day.tx day.date
    end

    # Antenna-frame geometry (Wood & Inan 2004): NS axis at azimuth ρ, EW axis
    # at 90°+ρ+ξ (azimuths CW from true north). The measured phasors are
    # m = M·[B_N; B_E] with det M = cos ξ; the composed map applied per sample is
    # T = R(θ_az + 90°)·M⁻¹, in closed form via a = α−ρ, b = α−ρ−ξ. ξ = ±90°
    # is a degenerate (parallel-axis) pair; the 1/cos ξ factor also sets the
    # amplitude-noise magnification of the de-skew.
    α  = deg2rad(bearing_deg) + π/2
    ρ  = deg2rad(rho_deg)
    ξ  = deg2rad(xi_deg)
    abs(xi_deg) < 90.0 ||
        error("rotate_day: |xi_deg| ≥ 90° is a degenerate antenna pair (det = cos ξ ≤ 0).")
    cξ = cos(ξ)
    a, b = α - ρ, α - ρ - ξ
    t11, t12 =  cos(b) / cξ, -sin(a) / cξ
    t21, t22 = sin(b) / cξ, cos(a) / cξ

    pNS, pEW = polarity.NS, polarity.EW
    dϕEW     = 90.0 * offset_m

    n = length(day.time)
    Br_amp   = Vector{Float64}(undef, n);  Bazi_amp = Vector{Float64}(undef, n)
    Br_pha   = Vector{Float64}(undef, n);  Bazi_pha = Vector{Float64}(undef, n)

    @inbounds for i in 1:n
        vNS = pNS * day.NS_amp[i] * cis(deg2rad(day.NS_pha[i]))
        vEW = pEW * day.EW_amp[i] * cis(deg2rad(day.EW_pha[i] + dϕEW))
        Br   = t11 * vNS + t12 * vEW
        Bazi = t21 * vNS + t22 * vEW
        Br_amp[i]   = abs(Br);   Br_pha[i]   = rad2deg(angle(Br))
        Bazi_amp[i] = abs(Bazi); Bazi_pha[i] = rad2deg(angle(Bazi))
    end

    return RotatedDay(day.date, day.rx, day.tx, day.Fc, day.Fs, day.time,
                      float(bearing_deg), Br_amp, Bazi_amp, Br_pha, Bazi_pha,
                      polarity, Int(offset_m), float(rho_deg), float(xi_deg),
                      day.params, String(baseline_label))
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
                      Br_pha, Bazi_pha, target.polarity, target.offset_m, target.rho_deg, target.xi_deg,
                      target.params, String(label))
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
                      Br_pha, Bazi_pha, target.polarity, target.offset_m, target.rho_deg, target.xi_deg,
                      target.params, String(label))
end