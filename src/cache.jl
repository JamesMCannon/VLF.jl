# ============================================================================
# cache.jl  —  Per-receiver JLD2 cache.
#
# Layout (per receiver path, self-contained so a site can be synced alone):
#
#   <rx_root>/.vlfcache/
#       index.jld2                          # quick lookup of cached entries
#       raw/<yymmdd>_<RX>_<TX>_<CH>_<Q>.jld2
#       processed/<yymmdd>_<RX>_<TX>.jld2
#
# "index" — deliberately not "manifest": Manifest.toml is Julia's pinned
# dependency graph, and "registry" is also taken by Pkg. Entry paths are
# deterministic from the key, so correctness never depends on the index; it is
# purely an acceleration for "what do I already have cached?" queries.
#
# Every entry stores `schema_version` as a sibling key, checked *before* the
# struct is deserialized, so a layout change can be detected without trying to
# reconstruct an incompatible type.
# ============================================================================

"""
    VLFCache(rx_root) -> VLFCache

Open (creating if needed) the cache rooted at `rx_root/.vlfcache`. `rx_root` is
one receiver's directory in your `data/<rx>/...` tree.
"""
struct VLFCache
    root::String   # the .vlfcache directory

    # Inner constructor: defining it suppresses the default field constructor
    # `VLFCache(::String)`, which would otherwise be MORE specific than an outer
    # `VLFCache(::AbstractString)` and silently win for String args — storing the
    # raw rx_root as `root` and skipping the `.vlfcache` resolution + mkpath.
    function VLFCache(rx_root::AbstractString)
        root = joinpath(rx_root, ".vlfcache")
        mkpath(joinpath(root, "raw"))
        mkpath(joinpath(root, "processed"))
        return new(root)
    end
end

# --- entry paths (deterministic from key) ----------------------------------
_raw_name(k::DataKey) = string(Dates.format(k.date, "yymmdd"), "_", k.rx, "_",
                               k.tx, "_", chstr(k.channel), "_", qstr(k.quantity), ".jld2")
_proc_name(date, rx, tx) = string(Dates.format(date, "yymmdd"), "_", rx, "_", tx, ".jld2")

raw_path(c::VLFCache, k::DataKey) = joinpath(c.root, "raw", _raw_name(k))
processed_path(c::VLFCache, date::Date, rx, tx) =
    joinpath(c.root, "processed", _proc_name(date, Symbol(rx), Symbol(tx)))

# --- cache index ------------------------------------------------------------
const _INDEX_KEYS = (:raw, :processed)

_index_file(c::VLFCache) = joinpath(c.root, "index.jld2")

"""
    load_index(cache) -> Dict

Return the cache index: `Dict(:raw => Set{DataKey}, :processed => Set{Tuple})`.
Rebuilt empty if absent or stale; treated as advisory only.
"""
function load_index(c::VLFCache)
    f = _index_file(c)
    if isfile(f)
        try
            return jldopen(f, "r") do io
                Dict(:raw => io["raw"]::Set{DataKey},
                     :processed => io["processed"]::Set{Tuple{Date,Symbol,Symbol}})
            end
        catch
            # fall through to a fresh index
        end
    end
    return Dict(:raw => Set{DataKey}(),
                :processed => Set{Tuple{Date,Symbol,Symbol}}())
end

function _save_index(c::VLFCache, idx)
    jldsave(_index_file(c); raw = idx[:raw], processed = idx[:processed])
    return nothing
end

function _index_add_raw!(c::VLFCache, k::DataKey)
    idx = load_index(c)
    push!(idx[:raw], k)
    _save_index(c, idx)
end

function _index_add_processed!(c::VLFCache, date, rx, tx)
    idx = load_index(c)
    push!(idx[:processed], (date, Symbol(rx), Symbol(tx)))
    _save_index(c, idx)
end

# --- low-level JLD2 read/write with schema gate -----------------------------
function _save_entry(path, obj)
    jldsave(path; schema_version = SCHEMA_VERSION, payload = obj)
    return path
end

# Returns the stored object if present and schema-current, else `nothing`.
function _load_entry(path)
    isfile(path) || return nothing
    io = jldopen(path, "r")
    try
        (haskey(io, "schema_version") && io["schema_version"] == SCHEMA_VERSION) || return nothing
        return io["payload"]
    catch
        return nothing
    finally
        close(io)
    end
end

# --- raw tier ---------------------------------------------------------------
"""
    save_raw(cache, day) -> String

Write `day` to the cache and record it in the index. Returns the path.
"""
function save_raw(c::VLFCache, day::RawDay)
    k = DataKey(day.date, day.rx, day.tx, day.channel, day.quantity)
    p = _save_entry(raw_path(c, k), day)
    _index_add_raw!(c, k)
    return p
end

"""
    load_raw(cache, key) -> Union{RawDay,Nothing}

Load a cached [`RawDay`](@ref), or `nothing` if absent or schema-stale.
"""
load_raw(c::VLFCache, k::DataKey) =
    let obj = _load_entry(raw_path(c, k)); obj isa RawDay ? obj : nothing end

"""
    get_raw(cache, source_folder, key; recompute=false) -> Union{RawDay,Nothing}

Return the [`RawDay`](@ref) for `key`: from cache if present (and `recompute` is
false), otherwise build it from the `.mat` files in `source_folder`, cache it,
and return it. `nothing` if no source files exist. Since `.mat` files are
immutable, "exists ⇒ trust it" is safe for the raw tier.
"""
function get_raw(c::VLFCache, source_folder::AbstractString, k::DataKey; recompute::Bool = false)
    if !recompute
        cached = load_raw(c, k)
        cached !== nothing && return cached
    end
    day = build_rawday(source_folder, k)
    day === nothing && return nothing
    save_raw(c, day)
    return day
end