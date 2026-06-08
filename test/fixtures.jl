# ============================================================================
# fixtures.jl  —  Helpers that synthesize .mat files so the test suite is
# self-contained and needs no real receiver data.
#
# Field set matches what VLF.read_mat_partial reads:
#   start_hour, start_minute, start_second, Fs, Fc, data
# (date / rx / tx / channel / quantity come from the FILENAME, not contents.)
# ============================================================================

using MAT

"Write one source .mat with the given header + data. Returns the path."
function write_mat_file(path; start_hour = 0, start_minute = 0, start_second = 0,
                        Fs = 1.0, Fc = 24.8e3, data = zeros(10))
    matwrite(path, Dict(
        "start_hour"   => Float64(start_hour),
        "start_minute" => Float64(start_minute),
        "start_second" => Float64(start_second),
        "Fs"           => Float64(Fs),
        "Fc"           => Float64(Fc),
        "data"         => Float64.(collect(data)),
    ))
    return path
end

"Create an AVID-format file (e.g. FSI250601000000_NLK_EW_A.mat) in `dir`. Returns its name."
function make_avid(dir; rx = "FSI", date = "250601", time = "000000",
                   tx = "NLK", ch = "EW", q = "A", kwargs...)
    fname = string(rx, date, time, "_", tx, "_", ch, "_", q, ".mat")
    write_mat_file(joinpath(dir, fname); kwargs...)
    return fname
end

"Create an AWESOME-format file (e.g. JU250701000000NLK_100A.mat) in `dir`. Returns its name."
function make_awesome(dir; rx = "JU", date = "250701", time = "000000",
                      tx = "NLK", ch = "100", q = "A", kwargs...)
    fname = string(rx, date, time, tx, "_", ch, q, ".mat")
    write_mat_file(joinpath(dir, fname); kwargs...)
    return fname
end

"""
Write a calibration .mat with NS/EW curves. Each curve is [freq_kHz response].
The response at the row nearest `Fc` (Hz) is `ns_factor` / `ew_factor`.
"""
function write_cal_file(path; ns_factor = 2.0, ew_factor = 3.0)
    freqs = [20.0, 24.8, 30.0]                       # kHz; 24.8 kHz == 24800 Hz
    matwrite(path, Dict(
        "CalibrationNumberNS" => hcat(freqs, [1.0, ns_factor, 9.0]),
        "CalibrationNumberEW" => hcat(freqs, [1.0, ew_factor, 9.0]),
    ))
    return path
end