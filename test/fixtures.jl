# ============================================================================
# fixtures.jl  —  Helpers that synthesize .mat files so the test suite is
# self-contained and needs no real receiver data.
#
# Field set matches what VLF.read_mat_partial reads:
#   start_hour, start_minute, start_second, Fs, Fc, data
# (date / rx / tx / rx_channel / quantity come from the FILENAME, not contents.)
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

"""
Write a minimal MAT Level 4 file with the fields read_mat_partial needs,
plus `pad` trailing zero bytes (mimicking receiver-written files).
"""
function write_v4_file(path; start_hour = 0, start_minute = 0, start_second = 0,
                       Fs = 1.0, Fc = 24.8e3, data = zeros(Float32, 10), pad = 0)
    open(path, "w") do io
        # Level 4 element: MOPT type, mrows, ncols, imagf, namlen, name\0, values
        writevar(name, vals, P) = begin
            write(io, Int32(10P), Int32(length(vals)), Int32(1), Int32(0),
                  Int32(length(name) + 1))
            write(io, codeunits(name), 0x00)
            write(io, vals)
        end
        writevar("start_hour",   [Float64(start_hour)],   0)
        writevar("start_minute", [Float64(start_minute)], 0)
        writevar("start_second", [Float64(start_second)], 0)
        writevar("Fs",           [Float64(Fs)],           0)
        writevar("Fc",           [Float64(Fc)],           0)
        writevar("data",         Float32.(collect(data)), 1)   # P=1 -> Float32
        write(io, zeros(UInt8, pad))
    end
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