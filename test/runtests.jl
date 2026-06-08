using Test
using Dates
using JLD2
using VLF

include("fixtures.jl")

secofday(h, m, s) = s + 60m + 3600h

@testset "VLF.jl" begin

# ---------------------------------------------------------------------------
@testset "types: time grid" begin
    g1 = timegrid(1.0)
    @test g1 isa TimeGrid
    @test length(g1) == 86_400
    @test first(g1) == 0.0
    @test last(g1) == 86_400 - 1/1.0          # 86399.0

    g50 = timegrid(50.0)
    @test length(g50) == 4_320_000
    @test last(g50) == 86_400 - 1/50.0        # 86399.98
    # element computed in twice precision -> no accumulation drift
    @test g50[end] == (length(g50) - 1) * (1/50.0)
end

@testset "types: DataKey & enums" begin
    k1 = DataKey(Date(2025,6,1), :FSI, :NLK, EW, AMPLITUDE)
    k2 = DataKey(Date(2025,6,1), :FSI, :NLK, EW, AMPLITUDE)
    k3 = DataKey(Date(2025,6,1), :FSI, :NLK, NS, AMPLITUDE)
    @test k1 == k2
    @test k1 != k3
    @test hash(k1) == hash(k2)                # usable as a Dict/Set key
    d = Dict(k1 => 1); @test d[k2] == 1
    @test VLF.chstr(EW) == "EW" && VLF.chstr(NS) == "NS"
    @test VLF.qstr(AMPLITUDE) == "A" && VLF.qstr(PHASE) == "B"
end

@testset "types: ProcessParams provenance" begin
    a = ProcessParams(cal_num = 2.0, tolerance = 10.0, unwrap = true)
    b = ProcessParams(cal_num = 2.0, tolerance = 10.0, unwrap = true)
    c = ProcessParams(cal_num = 3.0, tolerance = 10.0, unwrap = true)
    @test VLF.provenance_matches(a, b)
    @test !VLF.provenance_matches(a, c)
end

# ---------------------------------------------------------------------------
@testset "filenames: AVID + AWESOME parsing" begin
    k = parse_filename("FSI250601000000_NLK_EW_A.mat")
    @test k == DataKey(Date(2025,6,1), :FSI, :NLK, EW, AMPLITUDE)

    @test parse_filename("CAL250601000000_NML_NS_B.mat") ==
          DataKey(Date(2025,6,1), :CAL, :NML, NS, PHASE)

    # AWESOME: 100 -> EW, 101 -> NS, quantity glued to channel
    @test parse_filename("JU250701000000NLK_100A.mat") ==
          DataKey(Date(2025,7,1), :JU, :NLK, EW, AMPLITUDE)
    @test parse_filename("JU250701000000NLK_101B.mat") ==
          DataKey(Date(2025,7,1), :JU, :NLK, NS, PHASE)

    # partial files (different start time in name) collapse to one key
    @test parse_filename("FSI250601000000_NLK_EW_A.mat") ==
          parse_filename("FSI250601143015_NLK_EW_A.mat")

    @test parse_filename("not_a_data_file.txt") === nothing
    @test parse_filename("README.md") === nothing
end

@testset "filenames: scan + match in a folder" begin
    mktempdir() do dir
        make_avid(dir; rx="FSI", date="250601", time="000000", ch="EW", q="A")
        make_avid(dir; rx="FSI", date="250601", time="143000", ch="EW", q="A") # partial, same key
        make_avid(dir; rx="FSI", date="250601", ch="NS", q="A")                # different key
        write(joinpath(dir, "ignore.txt"), "x")

        key = DataKey(Date(2025,6,1), :FSI, :NLK, EW, AMPLITUDE)
        @test length(source_files_for(dir, key)) == 2
        @test length(scan_keys(dir)) == 2                                      # EW/A and NS/A
    end
end

# ---------------------------------------------------------------------------
@testset "mat: read_mat_partial seconds-from-midnight" begin
    mktempdir() do dir
        make_avid(dir; start_hour=1, start_minute=2, start_second=3,
                  Fs=1.0, Fc=24.8e3, data=ones(5))
        p = VLF.read_mat_partial(joinpath(dir, readdir(dir)[1]))
        @test p.start_seconds == secofday(1,2,3)        # 3723.0, NOT scaled by Fs
        @test p.Fs == 1.0 && p.Fc == 24.8e3
        @test p.data == ones(5)
    end
end

@testset "mat: build_rawday gridding" begin
    key = DataKey(Date(2025,6,1), :FSI, :NLK, EW, AMPLITUDE)

    # single partial at 14:30:15 lands at index 52216 (1 Hz), rest NaN
    mktempdir() do dir
        make_avid(dir; start_hour=14, start_minute=30, start_second=15,
                  Fs=1.0, data=[1.0,2.0,3.0])
        d = build_rawday(dir, key)
        @test d isa RawDay
        @test length(d.data) == 86_400
        @test d.data[52216:52218] == [1.0,2.0,3.0]
        @test isnan(d.data[52215]) && isnan(d.data[1])
        @test count(!isnan, d.data) == 3
    end

    # two partials (reboot) place by index regardless of file order
    mktempdir() do dir
        make_avid(dir; time="000100", start_minute=1, data=fill(2.0,5))   # 60 s
        make_avid(dir; time="000000", start_second=0, data=fill(1.0,10))  # 0 s
        d = build_rawday(dir, key)
        @test d.data[1:10] == ones(10)
        @test d.data[61:65] == fill(2.0,5)
        @test isnan(d.data[11])
    end

    # samples spilling past midnight are clipped, not an error
    mktempdir() do dir
        make_avid(dir; start_hour=23, start_minute=59, start_second=59, data=[7.0,8.0,9.0])
        d = build_rawday(dir, key)
        @test d.data[end] == 7.0
        @test count(!isnan, d.data) == 1
    end

    # mismatched Fs is warned and skipped
    mktempdir() do dir
        make_avid(dir; time="000000", Fs=1.0, data=ones(3))
        make_avid(dir; time="010000", start_hour=1, Fs=50.0, data=ones(3))
        d = @test_logs (:warn,) match_mode=:any build_rawday(dir, key)
        @test d.Fs == 1.0
        @test count(!isnan, d.data) == 3      # only the 1 Hz file placed
    end

    # no matching files -> nothing
    mktempdir() do dir
        @test build_rawday(dir, key) === nothing
    end
end

# ---------------------------------------------------------------------------
@testset "cache: round-trip (incl. TimeGrid + NaN)" begin
    mktempdir() do dir
        make_avid(dir; data=[1.0, NaN, 3.0])
        key = DataKey(Date(2025,6,1), :FSI, :NLK, EW, AMPLITUDE)
        cache = VLFCache(dir)

        day = build_rawday(dir, key)
        save_raw(cache, day)
        @test isfile(raw_path(cache, key))

        loaded = load_raw(cache, key)
        @test loaded isa RawDay
        @test loaded.time == day.time          # StepRangeLen survives JLD2
        @test loaded.time isa TimeGrid
        @test isequal(loaded.data, day.data)    # isequal so NaN == NaN
        @test loaded.Fc == day.Fc && loaded.Fs == day.Fs
        @test loaded.channel == EW && loaded.quantity == AMPLITUDE
    end
end

@testset "cache: schema gate rejects stale entries" begin
    mktempdir() do dir
        make_avid(dir; data=ones(4))
        key = DataKey(Date(2025,6,1), :FSI, :NLK, EW, AMPLITUDE)
        cache = VLFCache(dir)
        day = build_rawday(dir, key)

        # write with a bumped schema version -> load must decline it
        jldsave(raw_path(cache, key); schema_version = SCHEMA_VERSION + 1, payload = day)
        @test load_raw(cache, key) === nothing
    end
end

@testset "cache: get_raw ingests then serves from cache" begin
    mktempdir() do srcdir
        make_avid(srcdir; data=[5.0,6.0,7.0])
        key = DataKey(Date(2025,6,1), :FSI, :NLK, EW, AMPLITUDE)
        cache = VLFCache(srcdir)

        first = get_raw(cache, srcdir, key)          # builds + caches
        @test first isa RawDay

        # second call against an EMPTY folder must still succeed from cache
        mktempdir() do empty
            second = get_raw(cache, empty, key)
            @test isequal(second.data, first.data)
        end

        # index records the entry
        idx = load_index(cache)
        @test key in idx[:raw]
    end
end

# ---------------------------------------------------------------------------
@testset "process: calibration" begin
    g = timegrid(1.0)
    amp = RawDay(Date(2025,6,1), :FSI, :NLK, EW, AMPLITUDE, 24.8e3, 1.0, g, fill(2.0, length(g)))

    # scalar factor
    @test calibration_factor(amp, ProcessParams(cal_num=2.5)) == 2.5
    @test calibrate(amp, ProcessParams(cal_num=2.5))[1] == 5.0

    # from a calibration file: EW curve, nearest row to Fc=24.8 kHz
    mktempdir() do dir
        cal = write_cal_file(joinpath(dir, "cal.mat"); ns_factor=2.0, ew_factor=3.0)
        @test calibration_factor(amp, ProcessParams(cal_file=cal)) == 3.0
        ns = RawDay(Date(2025,6,1), :FSI, :NLK, NS, AMPLITUDE, 24.8e3, 1.0, g, fill(2.0, length(g)))
        @test calibration_factor(ns, ProcessParams(cal_file=cal)) == 2.0
    end

    @test_throws ErrorException calibration_factor(amp, ProcessParams())  # neither cal supplied
end

@testset "process: combine_quadrature" begin
    @test combine_quadrature([3.0], [4.0]) == [5.0]
    @test isnan(combine_quadrature([NaN], [4.0])[1])
    @test isnan(combine_quadrature([3.0], [NaN])[1])
end

@testset "process: unwrap_phase (and no mutation)" begin
    p = [10.0, 350.0, 20.0]
    u = unwrap_phase(p)
    @test u == [10.0, -10.0, 20.0]
    @test p == [10.0, 350.0, 20.0]            # regression: input untouched
    # NaN is passed through without poisoning
    @test isequal(unwrap_phase([0.0, NaN, 10.0]), [0.0, NaN, 10.0])
end

@testset "process: stitch_phase across NaN gaps" begin
    # +90 deg step across a gap is removed (tolerance window)
    in  = [0.0, 0.0, NaN, NaN, 90.0, 90.0]
    @test isequal(stitch_phase(in; unwrap=false), [0.0, 0.0, NaN, NaN, 0.0, 0.0])

    # a jump outside any n*90 +/- tolerance is left alone
    in2 = [0.0, NaN, 45.0, 45.0]
    @test isequal(stitch_phase(in2; unwrap=false), [0.0, NaN, 45.0, 45.0])

    # no gaps -> unchanged (unwrap off)
    @test stitch_phase([1.0, 2.0, 3.0]; unwrap=false) == [1.0, 2.0, 3.0]
end

@testset "process: clean_phase" begin
    @test isequal(clean_phase([10.0, NaN, 30.0], [5.0, 5.0, 5.0]), [5.0, NaN, 25.0])
end

# ---------------------------------------------------------------------------
@testset "process: build/get ProcessedDay end-to-end" begin
    mktempdir() do dir
        # four channel-files, all 1 Hz partials of length 10 starting at 00:00
        make_avid(dir; ch="EW", q="A", data=fill(3.0, 10))
        make_avid(dir; ch="NS", q="A", data=fill(4.0, 10))
        make_avid(dir; ch="EW", q="B", data=zeros(10))
        make_avid(dir; ch="NS", q="B", data=zeros(10))
        cache  = VLFCache(dir)
        params = ProcessParams(cal_num = 2.0)

        day = get_processed(cache, dir, :FSI, :NLK, Date(2025,6,1), params)
        @test day isa ProcessedDay
        @test length(day.combined_amp) == 86_400
        @test day.EW_amp[5] == 6.0                       # 3 * 2
        @test day.NS_amp[5] == 8.0                       # 4 * 2
        @test day.combined_amp[5] == 10.0                # sqrt(36+64)
        @test isnan(day.combined_amp[20_000])            # outside covered range
        @test VLF.provenance_matches(day.params, params)
        @test isfile(processed_path(cache, Date(2025,6,1), :FSI, :NLK))

        # wrong params + recompute=false -> warn but serve cached product
        day2 = @test_logs (:warn,) match_mode=:any begin
            get_processed(cache, dir, :FSI, :NLK, Date(2025,6,1),
                          ProcessParams(cal_num = 5.0))
        end
        @test day2.EW_amp[5] == 6.0                      # still the cal_num=2 product

        # recompute=true rebuilds with the new params
        day3 = get_processed(cache, dir, :FSI, :NLK, Date(2025,6,1),
                             ProcessParams(cal_num = 5.0); recompute = true)
        @test day3.EW_amp[5] == 15.0                     # 3 * 5
    end
end

@testset "process: missing amplitude errors" begin
    mktempdir() do dir
        make_avid(dir; ch="EW", q="A", data=ones(10))    # no NS amplitude
        cache = VLFCache(dir)
        @test_throws ErrorException get_processed(
            cache, dir, :FSI, :NLK, Date(2025,6,1), ProcessParams(cal_num=1.0))
    end
end

end # @testset "VLF.jl"