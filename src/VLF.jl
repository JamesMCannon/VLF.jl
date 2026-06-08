module VLF

using Dates
using MAT
using JLD2

include("types.jl")
include("mat.jl")
include("cache.jl")
include("process.jl")

# --- data model ---
export Channel, NS, EW
export Quantity, AMPLITUDE, PHASE
export DataKey, RawDay, ProcessedDay, ProcessParams
export TimeGrid, timegrid, SCHEMA_VERSION

# --- source files ---
export parse_filename, scan_keys, source_files_for, build_rawday

# --- cache ---
export VLFCache, get_raw, load_raw, save_raw, load_index
export raw_path, processed_path

# --- processing ---
export calibrate, calibration_factor, combine_quadrature
export pT_to_uVm, to_db, apply_amplitude_units, PT_TO_UVM
export amplitude, amplitudes
export unwrap_phase, stitch_phase, clean_phase, detrend_phase
export build_processed, get_processed

end # module VLF