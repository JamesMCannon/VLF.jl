# ============================================================================
# mat.jl  —  Reading source .mat files and gridding them into RawDay.
#
# Two filename formats are supported, detected from the name itself (NOT from a
# hardcoded list of which receivers are which — a site that changes software
# simply starts matching the other pattern):
#
#   AVID (11 receivers):   FSI250601000000_NLK_EW_A.mat
#                          rx date(yymmdd) time(hhmmss) _ tx _ ch(EW|NS) _ q(A|B)
#
#   AWESOME (2 receivers): JU250701000000NLK_100A.mat
#                          rx date(yymmdd) time(hhmmss) tx _ ch(100|101) q(A|B)
#                          with 100 -> EW, 101 -> NS
#
# Both formats have identical .mat internals, so a single reader handles all 13.
# ============================================================================

const RE_AVID = r"^(?<rx>[A-Z]+)(?<date>\d{6})(?<time>\d{6})_(?<tx>[A-Z]+)_(?<ch>EW|NS)_(?<q>[AB])\.mat$"
const RE_AWESOME = r"^(?<rx>[A-Z]+)(?<date>\d{6})(?<time>\d{6})(?<tx>[A-Z]+)_(?<ch>10[01])(?<q>[AB])\.mat$"

_channel_from_token(t) = (t == "EW" || t == "100") ? EW :
                         (t == "NS" || t == "101") ? NS :
                         error("Unrecognized channel token: $t")

_quantity_from_token(t) = t == "A" ? AMPLITUDE :
                         t == "B" ? PHASE :
                         error("Unrecognized quantity token: $t")

"""
    parse_filename(fname) -> Union{DataKey,Nothing}

Parse a source `.mat` filename into a [`DataKey`](@ref), trying the AVID format
then the AWESOME format. Returns `nothing` for names matching neither (so a
directory can be scanned without pre-filtering). The file start time is parsed
but intentionally discarded — it is not part of the key.
"""
function parse_filename(fname::AbstractString)
    m = match(RE_AVID, fname)
    m === nothing && (m = match(RE_AWESOME, fname))
    m === nothing && return nothing

    ds = m[:date]
    yy = parse(Int, ds[1:2]); mm = parse(Int, ds[3:4]); dd = parse(Int, ds[5:6])
    date = Date(2000 + yy, mm, dd)

    return DataKey(date, Symbol(m[:rx]), Symbol(m[:tx]),
                   _channel_from_token(m[:ch]), _quantity_from_token(m[:q]))
end

"""
    source_files_for(folder, key) -> Vector{String}

Filenames in `folder` whose [`DataKey`](@ref) equals `key`. Multiple results
mean partial files (e.g. a receiver reboot mid-day) that will be gridded
together into one [`RawDay`](@ref).
"""
function source_files_for(folder::AbstractString, key::DataKey)
    out = String[]
    for f in readdir(folder)
        k = parse_filename(f)
        k !== nothing && k == key && push!(out, f)
    end
    return sort!(out)
end

"""
    scan_keys(folder) -> Vector{DataKey}

Every distinct [`DataKey`](@ref) for which `folder` holds at least one `.mat`
file. Useful for discovering what is available before requesting it.
"""
function scan_keys(folder::AbstractString)
    keys = Set{DataKey}()
    for f in readdir(folder)
        k = parse_filename(f)
        k !== nothing && push!(keys, k)
    end
    return collect(keys)
end

# A single partial file's payload, in physical units (seconds, Hz).
struct _Partial
    start_seconds::Float64
    Fs::Float64
    Fc::Float64
    data::Vector{Float64}
end

"""
    read_mat_partial(path) -> _Partial

Read one `.mat` file's header + data. `start_seconds` is seconds-from-midnight
(`hour*3600 + minute*60 + second`); contrast with the legacy reader, which
left this in *sample* units by multiplying by `Fs` too early.
"""
function read_mat_partial(path::AbstractString)
    mf = matopen(path)
    try
        get1(name) = read(mf, name)[1]          # MAT scalars come back as 1×1
        hh = Int(get1("start_hour"))
        mn = Int(get1("start_minute"))
        ss = Int(get1("start_second"))
        Fs = Float64(get1("Fs"))
        Fc = Float64(get1("Fc"))
        data = Float64.(vec(read(mf, "data")))
        return _Partial(ss + 60mn + 3600hh, Fs, Fc, data)
    finally
        close(mf)
    end
end

"""
    build_rawday(folder, key) -> Union{RawDay,Nothing}

Read every source file matching `key` in `folder` and scatter each onto a fresh
NaN-filled full-day grid by index (`round(Int, start_seconds*Fs) + 1`). Returns
`nothing` if no matching files exist. Placement is positional, so file order is
irrelevant and there is no fragile same-day append logic. Samples that would
spill past midnight are clipped.
"""
function build_rawday(folder::AbstractString, key::DataKey)
    files = source_files_for(folder, key)
    isempty(files) && return nothing

    partials = [read_mat_partial(joinpath(folder, f)) for f in files]
    Fs = partials[1].Fs
    Fc = partials[1].Fc
    n  = round(Int, 86400 * Fs)
    data = fill(NaN, n)

    for (f, p) in zip(files, partials)
        if !isapprox(p.Fs, Fs; rtol = 1e-9)
            @warn "Skipping file with mismatched Fs" file=f expected=Fs got=p.Fs
            continue
        end
        base = round(Int, p.start_seconds * Fs) + 1
        stop = base + length(p.data) - 1
        lo = max(base, 1)
        hi = min(stop, n)
        hi < lo && continue
        data[lo:hi] .= @view p.data[(lo-base+1):(hi-base+1)]
    end

    return RawDay(key.date, key.rx, key.tx, key.channel, key.quantity,
                  Fc, Fs, timegrid(Fs), data)
end