# ============================================================================
# types.jl  —  Core data model for VLF.jl
#
#   RawDay        : one (date, rx, tx, channel, quantity), gridded to a full
#                   86400 s day, NaN where no sample exists. Mirrors a .mat file.
#   ProcessedDay  : one (date, rx, tx) bundle of calibrated/cleaned products
#                   (EW/NS/combined amplitude in pT, EW/NS cleaned phase).
#
# Time is stored as a derived StepRangeLen (TimeGrid) rather than a Vector:
# smaller on disk, and elements are computed as start + (i-1)*step in twice
# precision, so there is no accumulation drift across a 50 Hz day.
# ============================================================================

"""
    SCHEMA_VERSION

Layout version for serialized [`RawDay`](@ref)/[`ProcessedDay`](@ref) cache
entries. Bump this whenever the meaning or layout of a stored struct changes
(units, field order, sign conventions, ...). On load, a mismatch triggers a
rebuild rather than silently reinterpreting stale bytes. MAT-file immutability
protects the *inputs*; `SCHEMA_VERSION` protects the *cache layout*.

Bumped 1 → 2: [`ProcessParams`](@ref) (serialized inside every cached
[`ProcessedDay`](@ref)) gained a `slope` field and the phase pipeline now
slope-detrends, so v1 entries are rebuilt rather than reinterpreted. Amplitudes
remain stored in linear pT; µV/m and dB are applied at read time.

Bumped 2 → 3: [`ProcessParams`](@ref) `cal_num` widened from `Float64` to
`Union{Float64,NamedTuple}` to allow per-channel calibration (`cal_num =
(EW = a, NS = b)`). Since `ProcessParams` is serialized inside every cached
[`ProcessedDay`](@ref), v2 entries are rebuilt rather than reinterpreted.

Bumped 3 → 4: [`ProcessParams`](@ref) gained a `max_gap` field to control
phase unwrapping across long gaps. 

Bumped 4 → 5: ProcessParams gained transmitter-dropout fields (dropout_db,
dropout_window, dropout_pad, dropout_min_valid). The phase/amplitude pipeline
now NaN-masks dropouts detected on a near-field baseline amplitude, changing the
stored product, so v4 entries are rebuilt rather than reinterpreted.

Bumped 5 → 6: ProcessParams gained `dropout_label`, a provenance string identifying
externally-supplied dropout ranges (e.g. a network-coincidence mask) that are not
otherwise reproducible from the other fields. It changes the provenance identity of
a cached product, so v5 entries are rebuilt rather than reinterpreted.
"""
const SCHEMA_VERSION = 6

"""
    Channel

Physical antenna channel. Only the two physical channels are valid for raw
data; `COMBINED`/`DOMINANT` are products and live as explicit fields on
[`ProcessedDay`](@ref), so they are deliberately absent here.
"""
@enum Channel NS EW

"""
    Quantity

Measured quantity carried by a raw file. `AMPLITUDE` ↔ filename suffix `A`,
`PHASE` ↔ suffix `B`.
"""
@enum Quantity AMPLITUDE PHASE

"Concrete type of `range(0.0; step = 1/Fs, length = n)`."
const TimeGrid = StepRangeLen{Float64,Base.TwicePrecision{Float64},
                              Base.TwicePrecision{Float64},Int}

"""
    timegrid(Fs) -> TimeGrid

Canonical full-day grid for sample rate `Fs` (Hz): `0 : 1/Fs : 86400 - 1/Fs`,
i.e. `round(Int, 86400*Fs)` samples. Exact for the rates worked with:
1 Hz → 86400 samples, 50 Hz → 4_320_000 samples.
"""
timegrid(Fs::Real) = range(0.0; step = 1 / Fs, length = round(Int, 86400 * Fs))

# ----------------------------------------------------------------------------
# DataKey: identifies a raw entry. (date, rx, tx, channel, quantity).
# The filename's start-time field is intentionally NOT part of the key, so
# multiple partial files for the same day (receiver reboots) collapse to one.
# ----------------------------------------------------------------------------
"""
    DataKey(date, rx, tx, channel, quantity)

Identity of a single raw channel-day. Used as the cache key and to match
source `.mat` files. Does not include the file start time, so partial files
from the same day share one key.
"""
struct DataKey
    date::Date
    rx::Symbol
    tx::Symbol
    channel::Channel
    quantity::Quantity
end

Base.:(==)(a::DataKey, b::DataKey) =
    a.date == b.date && a.rx == b.rx && a.tx == b.tx &&
    a.channel == b.channel && a.quantity == b.quantity

Base.hash(k::DataKey, h::UInt) =
    hash(k.date, hash(k.rx, hash(k.tx, hash(k.channel, hash(k.quantity, h)))))

# ----------------------------------------------------------------------------
# ProcessParams: provenance for a ProcessedDay. Recording these lets a load
# warn when the cached product was built with settings other than those the
# caller is now asking for.
# ----------------------------------------------------------------------------
"""
    ProcessParams(; cal_file="", cal_num=-1.0, tolerance=10.0, unwrap=true,
                    baseline="", slope=nothing, max_gap=3600)

Settings that determine a [`ProcessedDay`](@ref). Stored as provenance so a
cached product can be checked against the parameters requested on load. The
defaults are shown in the signature above; pass only the keywords you need.

Amplitude *units* are intentionally NOT here: the cache always stores linear pT,
and pT↔µV/m / linear↔dB conversions are applied at read time (see
[`amplitude`](@ref) / [`amplitudes`](@ref)), so a viewer can switch units ad hoc
without recomputing.

Calibration (counts → pT):
- `cal_file::String = ""` — calibration file path, or `""` when a scalar is used.
- `cal_num::Union{Float64,NamedTuple} = -1.0` — calibration factor. Either a
  single scalar applied to *both* channels, or a per-channel NamedTuple
  `(EW = a, NS = b)` keyed by channel name. The scalar sentinel `-1.0` means
  "unset, read `cal_file` instead" (a NamedTuple is always considered set). When
  a NamedTuple is given it must contain the key for any channel actually being
  calibrated, else an error is raised; `cal_file` is ignored.

Phase:
- `tolerance::Float64 = 10.0` — ± window (deg) for n×90° phase-jump stitching.
- `unwrap::Bool = true` — whether phase is unwrapped before detrend/stitch.
- `baseline::String = ""` — provenance label for the phase baseline, or `""`.
- `max_gap::Float64 = 3600` — maximum gap length (s) for phase unwrapping.
- `slope::Union{Nothing,Float64} = nothing` — linear phase-detrend rate (deg/s).
  See [`detrend_phase`](@ref): with no reference baseline it detrends every
  sample as `slope · t` (`t` = seconds-from-midnight; `nothing` ⇒ no detrend,
  the original behavior). Where a reference baseline drops out but the target is
  valid, the reference is *extrapolated* from its last valid value at this slope
  so the result stays continuous (`nothing` ⇒ `NaN` across the dropout).

- `dropout_label::String = ""` — provenance for dropout ranges supplied to
  [`build_processed`](@ref) via `dropout_ranges` (e.g. from
  [`detect_dropouts_network`](@ref)). Caller-asserted: it labels the mask source
  so the cache can distinguish products, but is not verified against the ranges.

!!! note
    `slope` IS baked into the cached [`ProcessedDay`](@ref) (it changes the phase
    product), so changing it requires `recompute=true` in [`get_processed`](@ref);
    otherwise a provenance mismatch warns and the cached product is served. Unit
    choices (µV/m, dB) are not baked in and never need a recompute.
"""
Base.@kwdef struct ProcessParams
    cal_file::String   = ""
    cal_num::Union{Float64,NamedTuple} = -1.0
    tolerance::Float64 = 10.0
    unwrap::Bool       = true
    baseline::String   = ""
    slope::Union{Nothing,Float64} = nothing
    max_gap::Float64   = 3600      # s; gaps longer than this reset datum (no fold)
    dropout_db::Union{Nothing,Float64} = nothing   # dB below trailing median to flag a drop; nothing ⇒ off
    dropout_window::Float64 = 300.0                # trailing rolling-median window (s)
    dropout_pad::Float64    = 1.0                  # dilation each side of a flagged run (s)
    dropout_min_valid::Int  = 30                # minimum valid samples in a window to compute median
    dropout_label::String   = ""                # provenance for externally-applied dropout ranges (e.g. a network mask); "" ⇒ none
end

"True if two parameter sets would produce the same product."
provenance_matches(a::ProcessParams, b::ProcessParams) =
    a.cal_file == b.cal_file && a.cal_num == b.cal_num &&
    a.tolerance == b.tolerance && a.unwrap == b.unwrap &&
    a.baseline == b.baseline && a.slope == b.slope &&
    a.max_gap == b.max_gap &&
    a.dropout_db == b.dropout_db && a.dropout_window == b.dropout_window &&
    a.dropout_pad == b.dropout_pad && a.dropout_min_valid == b.dropout_min_valid &&
    a.dropout_label == b.dropout_label
# ----------------------------------------------------------------------------
# RawDay
# ----------------------------------------------------------------------------
"""
    RawDay

A single channel-day on the canonical [`timegrid`](@ref). `data` has length
`length(time)` with `NaN` wherever no source sample existed (gaps from reboots
or missing files). Holds *raw counts* for amplitude; calibration happens at the
processed tier.

# Fields
- `date::Date` — UTC day this entry covers.
- `rx::Symbol` — receiver code, e.g. `:FSI`.
- `tx::Symbol` — transmitter code, e.g. `:NLK`.
- `channel::Channel` — `NS` or `EW`.
- `quantity::Quantity` — `AMPLITUDE` or `PHASE`.
- `Fc::Float64` — transmitter center frequency (Hz).
- `Fs::Float64` — sample rate (Hz); 1.0 for narrowband, 50.0 for broadband.
- `time::TimeGrid` — seconds-from-midnight grid, `0 : 1/Fs : 86400 - 1/Fs`.
- `data::Vector{Float64}` — samples aligned to `time`; raw counts (`AMPLITUDE`)
  or degrees (`PHASE`), with `NaN` in any gap.

See also [`ProcessedDay`](@ref), [`build_rawday`](@ref), [`get_raw`](@ref).
"""
struct RawDay
    date::Date
    rx::Symbol
    tx::Symbol
    channel::Channel
    quantity::Quantity
    Fc::Float64
    Fs::Float64
    time::TimeGrid
    data::Vector{Float64}
end

# ----------------------------------------------------------------------------
# ProcessedDay
# ----------------------------------------------------------------------------
"""
    ProcessedDay

Calibrated/cleaned products for one `(date, rx, tx)`, all co-registered on the
same [`timegrid`](@ref). Amplitudes are stored in linear pT (convert to µV/m
and/or dB at read time with [`amplitude`](@ref)/[`amplitudes`](@ref)); phases are
in degrees, cleaned (optionally baseline-subtracted and/or slope-detrended) and
stitched (optionally unwrapped). `dominant` phase is reserved for a future
revision. `params` records provenance.

# Fields
- `date::Date` — UTC day this product covers.
- `rx::Symbol` — receiver code, e.g. `:FSI`.
- `tx::Symbol` — transmitter code, e.g. `:NLK`.
- `Fc::Float64` — transmitter center frequency (Hz).
- `Fs::Float64` — sample rate (Hz).
- `time::TimeGrid` — seconds-from-midnight grid shared by every product below.
- `EW_amp::Vector{Float64}` — east–west amplitude (pT), `NaN` in gaps.
- `NS_amp::Vector{Float64}` — north–south amplitude (pT).
- `combined_amp::Vector{Float64}` — `sqrt(NS² + EW²)` (pT); `NaN` unless both exist.
- `EW_pha::Vector{Float64}` — east–west phase (deg), cleaned + stitched.
- `NS_pha::Vector{Float64}` — north–south phase (deg), cleaned + stitched.
- `params::ProcessParams` — provenance: the settings this product was built with.

See also [`RawDay`](@ref), [`build_processed`](@ref), [`get_processed`](@ref).
"""
struct ProcessedDay
    date::Date
    rx::Symbol
    tx::Symbol
    Fc::Float64
    Fs::Float64
    time::TimeGrid
    EW_amp::Vector{Float64}        # pT
    NS_amp::Vector{Float64}        # pT
    combined_amp::Vector{Float64}  # pT, NS/EW in quadrature
    EW_pha::Vector{Float64}        # deg, cleaned + stitched
    NS_pha::Vector{Float64}        # deg, cleaned + stitched
    # dominant_pha::Vector{Float64}  # TODO: phase of the higher-amplitude channel
    params::ProcessParams
end

"""
    RotatedDay

Horizontal magnetic flux density for one `(date, rx, tx)` expressed in the radial
(`r̂`, away from the source) and azimuthal (`φ̂`, 90° CCW from `r̂`) frame of Gross
et al. (2018), Eq. (5). Built from a [`ProcessedDay`](@ref) by [`rotate`](@ref).
All products share the parent's [`timegrid`](@ref); amplitudes are in pT, phases
in degrees in `(-180, 180]`. A sample is valid only where *both* parent channels
are valid, so `NaN` in either NS/EW amplitude or phase yields `NaN` in both
rotated components.

`bearing_deg` is `θ_az`: the rx→tx forward azimuth, degrees clockwise from true
north. `polarity` is the per-channel ±1 applied to each reconstructed phasor;
the default `(NS = 1, EW = -1)` reproduces Gross's stated winding (the −ŷ EW
surface vector carries the minus sign). `baseline_label` is `""` for an absolute
(un-referenced) product and names the reference receiver after
[`baseline_subtract`](@ref). `params` is carried through from the parent for
provenance.
"""
struct RotatedDay
    date::Date
    rx::Symbol
    tx::Symbol
    Fc::Float64
    Fs::Float64
    time::TimeGrid
    bearing_deg::Float64
    Br_amp::Vector{Float64}     # pT
    Bazi_amp::Vector{Float64}   # pT
    Br_pha::Vector{Float64}     # deg, (-180, 180]
    Bazi_pha::Vector{Float64}   # deg, (-180, 180]
    polarity::NamedTuple        # (NS = +1, EW = -1) default
    params::ProcessParams
    baseline_label::String
end

# ----------------------------------------------------------------------------
# String tokens used in filenames/cache names
# ----------------------------------------------------------------------------
chstr(c::Channel) = c === EW ? "EW" : "NS"
qstr(q::Quantity) = q === AMPLITUDE ? "A" : "B"

# ----------------------------------------------------------------------------
# Pretty printing (handy at the REPL while testing)
# ----------------------------------------------------------------------------
_coverage(v) = isempty(v) ? 0.0 : 100 * count(!isnan, v) / length(v)

function Base.show(io::IO, d::RawDay)
    print(io, "RawDay($(d.rx)←$(d.tx) $(chstr(d.channel))/$(qstr(d.quantity)) ",
          "$(d.date) | $(round(d.Fs; digits=3)) Hz, ",
          "$(round(_coverage(d.data); digits=1))% covered)")
end

function Base.show(io::IO, d::ProcessedDay)
    print(io, "ProcessedDay($(d.rx)←$(d.tx) $(d.date) | ",
          "$(round(d.Fs; digits=3)) Hz, amp $(round(_coverage(d.combined_amp); digits=1))% covered)")
end

# Multi-line display (fires when you type the variable at the REPL). The compact
# one-liners above are kept for when days appear inside a Vector, Dict, etc.
_cov(v) = string(round(_coverage(v); digits = 1), "%")

function Base.show(io::IO, ::MIME"text/plain", d::RawDay)
    println(io, "RawDay: ", d.rx, " ← ", d.tx, "   ",
            chstr(d.channel), "/", qstr(d.quantity), "   ", d.date)
    println(io, "  Fc       = ", d.Fc, " Hz")
    println(io, "  Fs       = ", d.Fs, " Hz")
    println(io, "  time     = ", d.time, "  (", length(d.time), " samples)")
    print(io,   "  data     = ", length(d.data), "-element Vector{Float64}, ",
          _cov(d.data), " covered (NaN in gaps)")
end

function Base.show(io::IO, ::MIME"text/plain", d::ProcessedDay)
    println(io, "ProcessedDay: ", d.rx, " ← ", d.tx, "   ", d.date)
    println(io, "  Fc   = ", d.Fc, " Hz    Fs = ", d.Fs, " Hz    ",
            length(d.time), " samples")
    println(io, "  EW_amp        (pT)   ", _cov(d.EW_amp), " covered")
    println(io, "  NS_amp        (pT)   ", _cov(d.NS_amp), " covered")
    println(io, "  combined_amp  (pT)   ", _cov(d.combined_amp), " covered")
    println(io, "  EW_pha        (deg)  ", _cov(d.EW_pha), " covered")
    println(io, "  NS_pha        (deg)  ", _cov(d.NS_pha), " covered")
    print(io,   "  params: ", d.params)
end