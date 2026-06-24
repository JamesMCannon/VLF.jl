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

@testset "process: per-channel cal_num (NamedTuple)" begin
    g  = timegrid(1.0)
    ew = RawDay(Date(2025,6,1), :FSI, :NLK, EW, AMPLITUDE, 24.8e3, 1.0, g, fill(2.0, length(g)))
    ns = RawDay(Date(2025,6,1), :FSI, :NLK, NS, AMPLITUDE, 24.8e3, 1.0, g, fill(2.0, length(g)))

    p = ProcessParams(cal_num = (EW = 3.0, NS = 5.0))
    @test calibration_factor(ew, p) == 3.0
    @test calibration_factor(ns, p) == 5.0
    @test calibrate(ew, p)[1] == 6.0
    @test calibrate(ns, p)[1] == 10.0

    # scalar still applies to both channels
    @test calibration_factor(ew, ProcessParams(cal_num = 4.0)) == 4.0
    @test calibration_factor(ns, ProcessParams(cal_num = 4.0)) == 4.0

    # a NamedTuple missing the channel's key errors
    @test_throws ErrorException calibration_factor(ns, ProcessParams(cal_num = (EW = 3.0,)))
end

@testset "process: per-channel cal_num through build/get" begin
    mktempdir() do dir
        make_avid(dir; ch="EW", q="A", data=fill(3.0, 10))
        make_avid(dir; ch="NS", q="A", data=fill(4.0, 10))
        cache = VLFCache(dir)
        day = get_processed(cache, dir, :FSI, :NLK, Date(2025,6,1),
                            ProcessParams(cal_num = (EW = 2.0, NS = 10.0)))
        @test day.EW_amp[5]       == 6.0                  # 3 * 2
        @test day.NS_amp[5]       == 40.0                 # 4 * 10
        @test day.combined_amp[5] == sqrt(6.0^2 + 40.0^2)
        @test isfile(processed_path(cache, Date(2025,6,1), :FSI, :NLK))   # complete -> cached
    end
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

@testset "process: missing amplitude -> NaN + warn, phase preserved" begin
    mktempdir() do dir
        # NS amplitude absent; EW amplitude + both phase channels present
        make_avid(dir; ch="EW", q="A", data=ones(10))
        make_avid(dir; ch="EW", q="B", data=fill(30.0, 10))
        make_avid(dir; ch="NS", q="B", data=fill(40.0, 10))
        cache = VLFCache(dir)

        day = @test_logs (:warn,) match_mode=:any get_processed(
            cache, dir, :FSI, :NLK, Date(2025,6,1), ProcessParams(cal_num=2.0))

        @test day isa ProcessedDay
        @test day.EW_amp[5]   == 2.0            # present: 1 * 2
        @test isnan(day.NS_amp[5])              # absent -> NaN
        @test isnan(day.combined_amp[5])        # NaN unless both channels
        @test day.EW_pha[5]   == 30.0           # phase still populated
        @test day.NS_pha[5]   == 40.0
        # missing channel usually reflects a permanent instrument gap, so the
        # product is cached like any other.
        @test isfile(processed_path(cache, Date(2025,6,1), :FSI, :NLK))
    end
end

@testset "process: no raw data -> get_processed returns nothing" begin
    mktempdir() do dir
        cache = VLFCache(dir)
        @test get_processed(cache, dir, :FSI, :NLK, Date(2025,6,1),
                            ProcessParams(cal_num=1.0)) === nothing
    end
end

# ---------------------------------------------------------------------------
@testset "process: amplitude units (µV/m, dB) — runtime conversion" begin
    @test VLF.PT_TO_UVM ≈ 299.70 atol = 0.05
    @test pT_to_uVm(2.0) ≈ 2.0 * VLF.PT_TO_UVM
    let v = pT_to_uVm([1.0, NaN])
        @test v[1] ≈ VLF.PT_TO_UVM
        @test isnan(v[2])
    end

    @test to_db(10.0) == 20.0                 # 20*log10(10)
    @test to_db(1.0)  == 0.0
    @test isnan(to_db(NaN))
    @test isnan(to_db(0.0))                   # non-positive -> NaN
    @test isnan(to_db(-3.0))

    # pure vector helper: default pT; efield -> µV/m; db -> dB; preserves NaN
    @test apply_amplitude_units([6.0]) == [6.0]
    @test apply_amplitude_units([6.0]; efield=true)[1] ≈ 6.0 * VLF.PT_TO_UVM
    let v = apply_amplitude_units([10.0, NaN]; db=true)
        @test v[1] == 20.0 && isnan(v[2])
    end
end

@testset "process: cache stays pT, accessors convert at read time" begin
    mktempdir() do dir
        make_avid(dir; ch="EW", q="A", data=fill(3.0, 10))
        make_avid(dir; ch="NS", q="A", data=fill(4.0, 10))
        cache = VLFCache(dir)

        d = get_processed(cache, dir, :FSI, :NLK, Date(2025,6,1),
                          ProcessParams(cal_num=2.0))
        # stored fields are linear pT
        @test d.EW_amp[5]       == 6.0          # 3 * 2
        @test d.combined_amp[5] == 10.0         # sqrt(6²+8²)

        # accessors convert without touching the cache / needing recompute
        @test amplitude(d, :EW)[5]              == 6.0
        @test amplitude(d, :EW; efield=true)[5] ≈ 6.0 * VLF.PT_TO_UVM
        @test amplitude(d, :combined; db=true)[5]            == 20*log10(10.0)
        @test amplitude(d, :combined; efield=true, db=true)[5] ==
              20*log10(10.0 * VLF.PT_TO_UVM)
        @test isnan(amplitude(d, :EW; db=true)[20_000])
        @test_throws ErrorException amplitude(d, :bogus)

        all3 = amplitudes(d; efield=true)
        @test all3.EW[5] ≈ 6.0 * VLF.PT_TO_UVM
        @test all3.combined[5] ≈ 10.0 * VLF.PT_TO_UVM

        # the cached struct is unchanged (still pT) after all those views
        @test d.EW_amp[5] == 6.0
    end
end

# ---------------------------------------------------------------------------
@testset "process: detrend_phase (3a no-ref, 3b ref dropouts w/ anchor)" begin
    t = [0.0, 1.0, 2.0, 3.0]

    # 3a: no reference -> subtract slope*t everywhere
    @test detrend_phase(fill(10.0, 4), nothing, 2.0, t) == [10.0, 8.0, 6.0, 4.0]
    # 3a: slope=nothing -> unchanged (recovers current behavior)
    @test detrend_phase(fill(10.0, 4), nothing, nothing, t) == fill(10.0, 4)
    @test isequal(detrend_phase([10.0, NaN, 30.0], nothing, nothing, [0.0,1.0,2.0]),
                  [10.0, NaN, 30.0])

    # 3b: reference present with dropouts; gap fills extrapolate the reference
    #     forward from its last valid value (continuous, no step at gap entry).
    target = [10.0, 20.0, 30.0, 40.0]
    ref    = [1.0, NaN, 3.0, NaN]
    #  i=1: ref valid -> 10-1=9; anchor (t0=0, r0=1)
    #  i=2: gap -> 20-(1+2*(1-0)) = 17
    #  i=3: ref valid -> 30-3=27; anchor (t0=2, r0=3)
    #  i=4: gap -> 40-(3+2*(3-2)) = 35
    @test detrend_phase(target, ref, 2.0, t) == [9.0, 17.0, 27.0, 35.0]
    # 3b: no slope -> NaN at the dropouts
    @test isequal(detrend_phase(target, ref, nothing, t), [9.0, NaN, 27.0, NaN])
    # target NaN where ref valid -> NaN regardless of slope
    @test isequal(detrend_phase([10.0, NaN, 30.0], [1.0,2.0,3.0], 5.0, [0.0,1.0,2.0]),
                  [9.0, NaN, 27.0])

    # continuity entering a gap: a flat target/ref ramps by exactly -slope/sample
    let out = detrend_phase(fill(100.0,4), [50.0,50.0,NaN,NaN], 2.0, [0.0,1.0,2.0,3.0])
        @test out[1] == 50.0 && out[2] == 50.0
        @test out[3] - out[2] == -2.0       # one step of the slope, no jump
        @test out[4] - out[3] == -2.0
    end

    # leading gap (no anchor yet) -> detrend from origin (== 3a)
    @test detrend_phase([10.0,20.0,30.0], [NaN,NaN,5.0], 2.0, [0.0,1.0,2.0]) ==
          [10.0, 18.0, 25.0]

    @test_throws DimensionMismatch detrend_phase([1.0,2.0], [1.0], nothing, [0.0,1.0])
end

@testset "process: slope detrend through build/get_processed" begin
    mktempdir() do dir
        make_avid(dir; ch="EW", q="A", data=fill(3.0, 10))
        make_avid(dir; ch="NS", q="A", data=fill(4.0, 10))
        make_avid(dir; ch="EW", q="B", data=fill(100.0, 10))
        make_avid(dir; ch="NS", q="B", data=zeros(10))
        cache = VLFCache(dir)

        # 3a: no baseline, slope=10 deg/s detrends the (unwrapped) phase
        day = get_processed(cache, dir, :FSI, :NLK, Date(2025,6,1),
                            ProcessParams(cal_num=1.0, slope=10.0))
        @test day.EW_pha[1]  == 100.0      # t=0
        @test day.EW_pha[2]  == 90.0       # t=1 -> 100-10
        @test day.EW_pha[10] == 10.0       # t=9 -> 100-90
        @test isnan(day.EW_pha[11])

        # default slope=nothing leaves phase untouched (recovers prior behavior)
        day0 = get_processed(cache, dir, :FSI, :NLK, Date(2025,6,1),
                             ProcessParams(cal_num=1.0); recompute=true)
        @test day0.EW_pha[1] == 100.0 && day0.EW_pha[10] == 100.0

        # 3b: baseline valid for samples 1..5 (R=50), dropout 6..10. Reference is
        # extrapolated from R0=50 at t0=4s, so the fill is continuous with the
        # cleaned value 50 (no +40 jump that an absolute detrend would produce).
        g  = timegrid(1.0)
        bl = fill(NaN, length(g)); bl[1:5] .= 50.0
        base_ew = RawDay(Date(2025,6,1), :REF, :NLK, EW, PHASE, 24.8e3, 1.0, g, bl)

        d3b = build_processed(cache, dir, :FSI, :NLK, Date(2025,6,1),
                              ProcessParams(cal_num=1.0, slope=2.0, unwrap=false);
                              baseline_ew = base_ew)
        @test d3b.EW_pha[5]  == 50.0                     # 100 - 50 (baseline valid)
        @test d3b.EW_pha[6]  == 100.0 - (50.0 + 2.0*(5-4))   # = 48, continuous
        @test d3b.EW_pha[10] == 100.0 - (50.0 + 2.0*(9-4))   # = 40
        @test d3b.EW_pha[6] - d3b.EW_pha[5] == -2.0      # no discontinuity at entry

        # 3b with no slope -> NaN across the baseline dropout
        d3b0 = build_processed(cache, dir, :FSI, :NLK, Date(2025,6,1),
                               ProcessParams(cal_num=1.0, unwrap=false);
                               baseline_ew = base_ew)
        @test d3b0.EW_pha[5] == 50.0
        @test isnan(d3b0.EW_pha[6]) && isnan(d3b0.EW_pha[10])
    end
end

@testset "process: provenance includes slope (units excluded)" begin
    base = ProcessParams(cal_num=2.0)
    @test VLF.provenance_matches(base, ProcessParams(cal_num=2.0))
    @test !VLF.provenance_matches(base, ProcessParams(cal_num=2.0, slope=0.0))  # nothing ≠ 0.0
    @test !VLF.provenance_matches(ProcessParams(slope=1.0), ProcessParams(slope=2.0))
    @test VLF.provenance_matches(ProcessParams(slope=1.0), ProcessParams(slope=1.0))

    # per-channel cal_num: scalar vs NamedTuple differ; NamedTuples compare field-wise
    @test !VLF.provenance_matches(ProcessParams(cal_num=2.0),
                                  ProcessParams(cal_num=(EW=2.0, NS=2.0)))
    @test VLF.provenance_matches(ProcessParams(cal_num=(EW=2.0, NS=3.0)),
                                 ProcessParams(cal_num=(EW=2.0, NS=3.0)))
    @test !VLF.provenance_matches(ProcessParams(cal_num=(EW=2.0, NS=3.0)),
                                  ProcessParams(cal_num=(EW=2.0, NS=4.0)))
end

@testset "process: ProcessedView (eager units projection)" begin
    mktempdir() do dir
        make_avid(dir; ch="EW", q="A", data=fill(3.0, 10))
        make_avid(dir; ch="NS", q="A", data=fill(4.0, 10))
        make_avid(dir; ch="EW", q="B", data=fill(30.0, 10))   # phase on EW only
        cache = VLFCache(dir)
        day = get_processed(cache, dir, :FSI, :NLK, Date(2025,6,1), ProcessParams(cal_num=2.0))
        # parent stored in linear pT: EW=6, NS=8, combined=10
        @test day.EW_amp[5] == 6.0 && day.NS_amp[5] == 8.0 && day.combined_amp[5] == 10.0

        # default view is pT and equals the parent's stored amplitudes
        vpt = view_units(day)
        @test vpt isa ProcessedView
        @test vpt.units == AmplitudeUnits(false, false)
        @test unit_label(vpt.units) == "pT"
        @test vpt.EW_amp[5] == 6.0 && vpt.combined_amp[5] == 10.0

        # µV/m view converts each stored pT field (combined is CONVERTED, not recombined)
        vE = view_units(day; efield=true)
        @test unit_label(vE.units) == "µV/m"
        @test vE.EW_amp[5]       ≈ 6.0  * VLF.PT_TO_UVM
        @test vE.NS_amp[5]       ≈ 8.0  * VLF.PT_TO_UVM
        @test vE.combined_amp[5] ≈ 10.0 * VLF.PT_TO_UVM

        # dB view
        vdb = view_units(day; db=true)
        @test unit_label(vdb.units) == "pT (dB)"
        @test vdb.combined_amp[5] == 20*log10(10.0)

        # metadata + phase forward to the parent; phase is borrowed (shared identity)
        @test vE.Fc == day.Fc && vE.Fs == day.Fs && vE.date == day.date
        @test vE.params === day.params
        @test vE.EW_pha === day.EW_pha            # read-through, not copied
        @test vE.EW_pha[5] == 30.0
        @test :EW_pha in propertynames(vE)

        # the view OWNS its amplitude arrays: mutating them cannot reach the parent
        vE.EW_amp[5] = -999.0
        @test day.EW_amp[5] == 6.0

        # parent handle + re-view never stacks conversions
        @test vE.parent === day
        vEE = view_units(vE; efield=true)         # re-view a view -> re-derives from pT parent
        @test vEE.EW_amp[5] ≈ 6.0 * VLF.PT_TO_UVM # single conversion, not squared
        @test vEE.parent === day

        # amplitude(::ProcessedView) returns the already-converted field; asking for
        # further conversion is a hard error (prevents double-conversion)
        @test amplitude(vE, :combined)[5] ≈ 10.0 * VLF.PT_TO_UVM
        @test_throws ErrorException amplitude(vE, :EW; db=true)
        @test_throws ErrorException amplitude(vE, :bogus)

        # views are never cacheable
        @test_throws ErrorException VLF._save_entry(
            processed_path(cache, Date(2025,6,1), :FSI, :NLK), vE)
    end
end

@testset "process: get_processed_view convenience" begin
    mktempdir() do dir
        make_avid(dir; ch="EW", q="A", data=fill(3.0, 10))
        make_avid(dir; ch="NS", q="A", data=fill(4.0, 10))
        cache = VLFCache(dir)

        v = get_processed_view(cache, dir, :FSI, :NLK, Date(2025,6,1),
                               ProcessParams(cal_num=2.0); efield=true)
        @test v isa ProcessedView
        @test v.combined_amp[5] ≈ 10.0 * VLF.PT_TO_UVM
        @test v.parent isa ProcessedDay
        # building the view still cached the canonical pT parent
        @test isfile(processed_path(cache, Date(2025,6,1), :FSI, :NLK))

        # no data -> nothing (mirrors get_processed)
        @test get_processed_view(cache, dir, :ZZZ, :NLK, Date(2025,6,1),
                                 ProcessParams(cal_num=2.0)) === nothing
    end
end

@testset "network-coherence dropouts" begin
    Fs, N = 1.0, 400
    tg = range(0.0; step = 1/Fs, length = N)
    mkamp(v) = RawDay(Date(2025,5,25), :RX, :NLK, EW, AMPLITUDE, 24.8e3, Fs, tg, v)

    base() = fill(100.0, N)                 # 40 dB everywhere
    drop!(v, r) = (v[r] .= 1.0; v)          # 0 dB → 40 dB drop, ≫ 10 dB

    # Three receivers, a coincident drop at 200:209.
    a = drop!(base(), 200:209); b = drop!(base(), 200:209); c = drop!(base(), 200:209)
    # Receiver `a` also has a lone fade at 300:304 the others lack.
    drop!(a, 300:304)
    d1, d2, d3 = mkamp(a), mkamp(b), mkamp(c)

    kw = (; drop_db = 10.0, window_s = 120.0, pad_s = 1.0, min_valid = 30)

    r = VLF.detect_dropouts_network([d1, d2, d3]; kw...)
    @test any(rng -> 200 in rng && 209 in rng, r)         # coincident drop caught
    @test !any(rng -> 300 in rng, r)                       # lone fade rejected (unanimous)

    # Integer override: one path now suffices, so the lone fade is caught.
    r1 = VLF.detect_dropouts_network([d1, d2, d3]; kw..., min_coincident = 1)
    @test any(rng -> 300 in rng, r1)

    # Blind-spot coverage: third receiver entirely down (NaN) must not block the
    # coincident drop seen by the two that are up.
    ddn = mkamp(fill(NaN, N))
    rdn = VLF.detect_dropouts_network([d2, d3, ddn]; kw..., min_receivers = 2)
    @test any(rng -> 200 in rng, rdn)

    # Grid mismatch is an error, not a silent drop.
    short = RawDay(Date(2025,5,25), :RX2, :NLK, EW, AMPLITUDE, 24.8e3, Fs,
                   range(0.0; step = 1/Fs, length = N - 1), base()[1:N-1])
    @test_throws ErrorException VLF.detect_dropouts_network([d2, short]; kw...)

    # Wrapper masks the passed days over the detected ranges; inputs untouched.
    out = VLF.mask_dropouts_network([d1, d2, d3]; kw...)
    @test out.ranges == r
    @test all(isnan, out.days[1].data[200:209])
    @test isequal(d1.data, a)                              # original not mutated
end

end # @testset "VLF.jl"