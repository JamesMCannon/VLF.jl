"""
    rotate(day::ProcessedDay, bearing_deg; polarity=(NS=1, EW=-1), baseline_label="") -> RotatedDay

Rotate `day`'s NS/EW horizontal phasors into the radial/azimuthal frame of Gross
et al. (2018), Eq. (5). `bearing_deg` is `θ_az` (rx→tx, degrees clockwise from
true north). Reconstructs `v_NS = polarity.NS · A_NS · e^{jψ_NS}` and
`v_EW = polarity.EW · A_EW · e^{jψ_EW}` from the calibrated pT amplitudes and
cleaned phases, then `[B_r; B_azi] = R(θ_az + π/2) · [v_NS; v_EW]` with
`R(α) = [cos α  sin α; −sin α  cos α]`.

Phases are absolute only when `day` was built with no phase baseline and
`slope === nothing`; a referenced/detrended parent yields referenced rotated
phases (amplitudes and the ratio ∠(−B_r/B_azi) are unaffected either way).
"""
function rotate(day::ProcessedDay, bearing_deg::Real;
                polarity::NamedTuple = (NS = 1, EW = -1),
                baseline_label::AbstractString = "")
    if !isempty(day.params.baseline) || day.params.slope !== nothing
        @warn "rotate: parent phase is referenced/detrended (baseline=\"$(day.params.baseline)\", \
               slope=$(day.params.slope)); rotated phases inherit it and are not absolute." day.rx day.tx day.date
    end

    α      = deg2rad(bearing_deg) + π/2
    cα, sα = cos(α), sin(α)
    pNS, pEW = polarity.NS, polarity.EW

    n = length(day.time)
    Br_amp   = Vector{Float64}(undef, n);  Bazi_amp = Vector{Float64}(undef, n)
    Br_pha   = Vector{Float64}(undef, n);  Bazi_pha = Vector{Float64}(undef, n)

    @inbounds for i in 1:n
        vNS = pNS * day.NS_amp[i] * cis(deg2rad(day.NS_pha[i]))
        vEW = pEW * day.EW_amp[i] * cis(deg2rad(day.EW_pha[i]))
        Br   =  cα * vNS + sα * vEW
        Bazi = -sα * vNS + cα * vEW
        Br_amp[i]   = abs(Br);   Br_pha[i]   = rad2deg(angle(Br))
        Bazi_amp[i] = abs(Bazi); Bazi_pha[i] = rad2deg(angle(Bazi))
    end

    return RotatedDay(day.date, day.rx, day.tx, day.Fc, day.Fs, day.time,
                      float(bearing_deg), Br_amp, Bazi_amp, Br_pha, Bazi_pha,
                      polarity, day.params, String(baseline_label))
end

"""
    baseline_subtract(target::RotatedDay, reference::RotatedDay; label=string(reference.rx)) -> RotatedDay

Per-component path differential: subtract `reference`'s rotated phase from
`target`'s, component by component, cancelling the common-mode transmitter source
phase. Amplitudes are the target's (phase-only reference). `target` and
`reference` must share `tx`, `date`, and grid; `bearing_deg` stays the target's
(the reference's own geometry is recorded in `baseline_label`).
"""
function baseline_subtract(target::RotatedDay, reference::RotatedDay;
                           label::AbstractString = string(reference.rx))
    (target.tx == reference.tx && target.date == reference.date) ||
        error("baseline_subtract: tx/date mismatch \
               (target $(target.tx)/$(target.date), reference $(reference.tx)/$(reference.date)).")
    (isapprox(target.Fs, reference.Fs; rtol = 1e-9) &&
     length(target.time) == length(reference.time)) ||
        error("baseline_subtract: target and reference grids differ.")

    Br_pha   = target.Br_pha   .- reference.Br_pha
    Bazi_pha = target.Bazi_pha .- reference.Bazi_pha

    return RotatedDay(target.date, target.rx, target.tx, target.Fc, target.Fs, target.time,
                      target.bearing_deg, copy(target.Br_amp), copy(target.Bazi_amp),
                      Br_pha, Bazi_pha, target.polarity, target.params, String(label))
end