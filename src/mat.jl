# ============================================================================
# mat.jl  —  Reading source .mat files and gridding them into RawDay.
#
# Two filename formats are supported, detected from the name itself:
#
#   AVID:       FSI250601000000_NLK_EW_A.mat
#               rx date(yymmdd) time(hhmmss) _ tx _ ch(EW|NS) _ q(A|B)
#
#   AWESOME:    JU250701000000NLK_100A.mat
#               rx date(yymmdd) time(hhmmss) tx _ ch(100|101) q(A|B)
#               with 100 -> EW, 101 -> NS
#
# Both formats have identical .mat internals, so a single reader handles all 13.
# ============================================================================

const RE_AVID = r"^(?<rx>[A-Z]+)(?<date>\d{6})(?<time>\d{6})_(?<tx>[A-Z]+)_(?<ch>EW|NS)_(?<q>[AB])\.mat$"
const RE_AWESOME = r"^(?<rx>[A-Z]+)(?<date>\d{6})(?<time>\d{6})(?<tx>[A-Z]+)_(?<ch>10[01])(?<q>[AB])\.mat$"

_rx_channel_from_token(t) = (t == "EW" || t == "100") ? EW :
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
                   _rx_channel_from_token(m[:ch]), _quantity_from_token(m[:q]))
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

# ----------------------------------------------------------------------------
# MAT Level 4 fallback reader.
#
# The receiver acquisition software writes MAT Level 4 (v4) files, some of
# which end in a run of zero bytes after the last variable. MATLAB ignores that
# padding, but MAT.jl's v4 scanner loops `while !eof` and raises EOFError on
# it, making an otherwise-intact file unreadable. This reader scans element by
# element and stops at the first position that is not a valid Level 4 element
# header. A zero tail is benign padding; a nonzero tail (a truly truncated
# variable) is reported with @warn while the intact variables are still
# returned — a missing *required* variable still errors in read_mat_partial,
# so build_rawday's per-file guard skips the file as before.
# ----------------------------------------------------------------------------

# P (precision) digit of the MOPT type code -> element type (Level 4 spec).
const _V4_ELTYPE = (Float64, Float32, Int32, Int16, UInt16, UInt8)

# Parse the five Int32s of a Level 4 element header at the current position.
# Returns nothing when they do not describe a valid element (e.g. padding).
function _v4_element_header(io::IO, swap::Bool)
    rd() = (x = read(io, Int32); swap ? bswap(x) : x)
    dtype = rd(); mrows = rd(); ncols = rd(); imagf = rd(); namlen = rd()
    0 <= dtype <= 9999 || return nothing
    M = div(dtype, 1000) % 10          # numeric format (0 = IEEE little-endian)
    O = div(dtype, 100)  % 10          # reserved, always 0
    P = div(dtype, 10)   % 10          # element precision
    T = dtype % 10                     # 0 numeric, 1 text, 2 sparse
    ok = 0 <= M <= 4 && O == 0 && 0 <= P <= 5 && 0 <= T <= 2 &&
         mrows >= 0 && ncols >= 0 && 0 <= imagf <= 1 && namlen >= 1
    return ok ? (; P, T, mrows, ncols, imagf, namlen) : nothing
end

# Byte order of a Level 4 file, or nothing if the file does not begin with a
# valid Level 4 element header under either byte order (i.e. is not v4).
function _v4_byte_order(io::IO)
    for swap in (false, true)
        seekstart(io)
        try
            _v4_element_header(io, swap) !== nothing && return swap
        catch e
            e isa EOFError || rethrow()
        end
    end
    return nothing
end

"""
    _read_v4_vars(path) -> Union{Dict{String,Any},Nothing}

Tolerant MAT Level 4 reader, used as a fallback when MAT.jl fails on `path`.
Returns `nothing` when the file is not Level 4 (caller rethrows the original
exception). Numeric variables come back as `Matrix`; text (T = 1) as a single
`String` in column-major order (sufficient here — `read_mat_partial` only uses
numeric variables).
"""
function _read_v4_vars(path::AbstractString)
    raw = read(path)                   # source files are small (≤ ~35 MB/day)
    io  = IOBuffer(raw)
    swap = _v4_byte_order(io)
    swap === nothing && return nothing
    seekstart(io)

    vars = Dict{String,Any}()
    while length(raw) - position(io) >= 20
        start = position(io)
        hdr = _v4_element_header(io, swap)
        hdr === nothing && (seek(io, start); break)      # padding / garbage
        elT   = _V4_ELTYPE[hdr.P + 1]
        nvals = hdr.mrows * hdr.ncols
        need  = hdr.namlen + sizeof(elT) * nvals * (hdr.imagf + 1)
        length(raw) - position(io) < need && (seek(io, start); break)  # truncated

        name = String(strip(String(read(io, hdr.namlen)), '\0'))
        readvals() = (v = Vector{elT}(undef, nvals); read!(io, v);
                      swap ? bswap.(v) : v)
        if hdr.T == 1                                    # text: char codes
            vars[name] = String(UInt8.(readvals()))
        else                                             # numeric
            re = readvals()
            v  = hdr.imagf == 1 ? complex.(re, readvals()) : re
            vars[name] = reshape(v, hdr.mrows, hdr.ncols)
        end
    end

    tail = @view raw[position(io)+1:end]
    if !isempty(tail) && !all(iszero, tail)
        @warn "MAT v4 file has an unreadable nonzero tail; returning intact variables" file=path tail_bytes=length(tail)
    end
    return vars
end

"""
    read_mat_partial(path) -> _Partial

Read one `.mat` file's header + data. `start_seconds` is seconds-from-midnight
(`hour*3600 + minute*60 + second`). MAT v5/v7.3 files go through MAT.jl; if
that fails, receiver-written MAT Level 4 files (which may carry trailing zero
padding that MAT.jl's v4 scanner cannot step over) are retried with the
tolerant `_read_v4_vars` reader. Non-v4 failures propagate unchanged.
"""
function read_mat_partial(path::AbstractString)
    needed = ("start_hour", "start_minute", "start_second", "Fs", "Fc", "data")
    vars = try
        mf = matopen(path)
        try
            Dict{String,Any}(n => read(mf, n) for n in needed)
        finally
            close(mf)
        end
    catch
        v4 = _read_v4_vars(path)
        v4 === nothing && rethrow()
        v4
    end
    getvar(n) = haskey(vars, n) ? vars[n] :
                error("required variable '$n' missing from $path")
    get1(n) = getvar(n)[1]                  # MAT scalars come back as 1×1
    hh = Int(get1("start_hour"))
    mn = Int(get1("start_minute"))
    ss = Int(get1("start_second"))
    Fs = Float64(get1("Fs"))
    Fc = Float64(get1("Fc"))
    data = Float64.(vec(getvar("data")))
    return _Partial(ss + 60mn + 3600hh, Fs, Fc, data)
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

    # Per-file guard: a truncated or corrupt .mat costs only its own segment
    # (warn + skip), not the whole channel-day. Fs/Fc template comes from the
    # first READABLE file.
    partials = _Partial[]
    kept     = String[]
    for f in files
        try
            push!(partials, read_mat_partial(joinpath(folder, f)))
            push!(kept, f)
        catch e
            @warn "Skipping unreadable source file" file = joinpath(folder, f) exception = e
        end
    end
    isempty(partials) && return nothing

    Fs = partials[1].Fs
    Fc = partials[1].Fc
    n  = round(Int, 86400 * Fs)
    data = fill(NaN, n)

    for (f, p) in zip(kept, partials)
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

    return RawDay(key.date, key.rx, key.tx, key.rx_channel, key.quantity,
                  Fc, Fs, timegrid(Fs), data)
end
