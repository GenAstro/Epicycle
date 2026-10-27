# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# Frame theory and Earth orientation parameter tables.
#
# Truth: the IERS finals files themselves. Three rows of each series are written below exactly as
# IERS publishes them, loaded through `eop_load`, and read back at a row's date, where the table
# must return the published value unchanged.
#
# The session's theory and tables are shared by every later test, so each testset puts back what
# it found.

using Test
using AstroUniverse
using AstroUniverse: EopIau1980, EopIau2000A

const _EOP_HEADER = "MJD;Year;Month;Day;Type;x_pole;sigma_x_pole;y_pole;sigma_y_pole;x_rate;" *
    "sigma_x_rate;y_rate;sigma_y_rate;Type;UT1-UTC;sigma_UT1-UTC;LOD;sigma_LOD;Type;dPsi;" *
    "sigma_dPsi;dEpsilon;sigma_dEpsilon;dX;sigma_dX;dY;sigma_dY;Type;bulB/x_pole;bulB/y_pole;" *
    "Type;bulB/UT-UTC;Type;bulB/dPsi;bulB/dEpsilon;bulB/dX;bulB/dY"

# finals2000A.all.csv, 1973-01-02 to 1973-01-04: carries the CIP offsets dX, dY.
const _EOP_2000A_ROWS = """
41684;1973;01;02;final;0.120733;0.009786;0.136966;0.015902;;;;;final;0.8084178;0.0002710;0.0000;0.1916;prediction;;;;;-0.766;0.199;-0.720;0.300;final;.143000;.137000;final;.8075000;final;;;-18.637;-3.667
41685;1973;01;03;final;0.118980;0.011039;0.135656;0.013616;;;;;final;0.8056163;0.0002710;3.5563;0.1916;prediction;;;;;-0.751;0.199;-0.701;0.300;final;.141000;.134000;final;.8044000;final;;;-18.636;-3.571
41686;1973;01;04;final;0.117227;0.011039;0.134348;0.013616;;;;;final;0.8027895;0.0002710;2.6599;0.1916;prediction;;;;;-0.738;0.199;-0.662;0.300;final;.139000;.131000;final;.8012000;final;;;-18.669;-3.621
"""

# finals.all.csv, same dates: carries the nutation corrections dPsi, dEpsilon instead.
const _EOP_1980_ROWS = """
41684;1973;01;02;final;0.120733;0.009786;0.136966;0.015902;;;;;final;0.8084178;0.0002710;0.0000;0.1916;prediction;44.969;.500;2.839;.300;;;;;final;.143000;.137000;final;.8075000;final;.000;.000;;
41685;1973;01;03;final;0.118980;0.011039;0.135656;0.013616;;;;;final;0.8056163;0.0002710;3.5563;0.1916;prediction;45.005;.500;2.762;.300;;;;;final;.141000;.134000;final;.8044000;final;.000;.000;;
41686;1973;01;04;final;0.117227;0.011039;0.134348;0.013616;;;;;final;0.8027895;0.0002710;2.6599;0.1916;prediction;45.122;.500;2.851;.300;;;;;final;.139000;.131000;final;.8012000;final;.000;.000;;
"""

const _JD_1973_01_03 = 2441685.5          # MJD 41685

function _eop_file(rows)
    path = tempname() * ".csv"
    write(path, _EOP_HEADER * "\n" * rows)
    return path
end

const AU = AstroUniverse
_url(path) = "file:///" * replace(path, "\\" => "/")

# Every load in this file reads local fixtures: one source per theory served from a file, a cache
# in a temporary directory, an empty directory in place of SatelliteToolboxTransformations' older
# cache, and "today" set to a day the fixtures cover. Nothing reaches the network or a real cache
# unless the live test at the end is asked for. The settings found are put back at the end.
const _EOP_HOOKS_FOUND = (AU._EOP_SOURCES[], AU._EOP_CACHE_DIR[], AU._EOP_LEGACY_DIR[],
                          AU._EOP_AGE_LIMIT[], AU._EOP_TODAY_MJD[], AU._EOP_REPLACE[])
const _FIXTURE_SOURCES = Dict(
    FK5     => [AU._EopSource(_url(_eop_file(_EOP_1980_ROWS)), :csv)],
    IAU2006 => [AU._EopSource(_url(_eop_file(_EOP_2000A_ROWS)), :csv)])

function _use_fixture_sources!()
    AU._EOP_SOURCES[]    = theory -> _FIXTURE_SOURCES[typeof(theory)]
    AU._EOP_CACHE_DIR[]  = mktempdir()
    AU._EOP_LEGACY_DIR[] = mktempdir()
    AU._EOP_AGE_LIMIT[]  = AU._EOP_MAX_AGE
    AU._EOP_TODAY_MJD[]  = 41685.0
    AU._EOP_REPLACE[]    = nothing
end
_use_fixture_sources!()

"""Run `f` and put the session's theory and both table slots, loaded or not, back afterwards."""
function _with_eop_restored(f)
    theory = frame_theory()
    t1980, t2000a = AU._eop_1980[], AU._eop_2000a[]
    try
        f()
    finally
        set_frame_theory!(theory)
        AU._eop_1980[], AU._eop_2000a[] = t1980, t2000a
    end
end

@testset "EOP — the frame theory selects the default table" begin
    _with_eop_restored() do
        @test frame_theory() isa IAU2006                     # the documented default
        @test FK5() isa AbstractFrameTheory && IAU2006() isa AbstractFrameTheory

        @test set_frame_theory!(FK5()) isa FK5
        @test frame_theory() isa FK5
        @test eop() isa EopIau1980
        @test eop() === eop(FK5())

        set_frame_theory!(IAU2006())
        @test eop() isa EopIau2000A
        @test eop() === eop(IAU2006())
    end
end

@testset "EOP — a file loads into its theory's slot and reads back the published values" begin
    _with_eop_restored() do
        p2000a, p1980 = _eop_file(_EOP_2000A_ROWS), _eop_file(_EOP_1980_ROWS)
        before_1980 = eop(FK5())

        t = eop_load(p2000a; theory = IAU2006())
        @test t isa EopIau2000A
        @test eop(IAU2006()) === t
        @test eop(FK5()) === before_1980                     # the other slot is untouched
        @test t.x(_JD_1973_01_03)        ≈ 0.118980  atol = 1e-12   # arcsec
        @test t.y(_JD_1973_01_03)        ≈ 0.135656  atol = 1e-12
        @test t.Δut1_utc(_JD_1973_01_03) ≈ 0.8056163 atol = 1e-12   # s
        @test t.δx(_JD_1973_01_03)       ≈ -0.751    atol = 1e-12   # mas

        # The theory keyword defaults to the session's theory.
        set_frame_theory!(FK5())
        f = eop_load(p1980)
        @test f isa EopIau1980
        @test eop(FK5()) === f
        @test f.δΔψ(_JD_1973_01_03) ≈ 45.005 atol = 1e-12            # mas
        @test f.δΔϵ(_JD_1973_01_03) ≈ 2.762  atol = 1e-12

        rm(p2000a; force = true); rm(p1980; force = true)
        @test_throws ArgumentError eop_load(joinpath(tempdir(), "no_such_eop_file.csv"))
    end
end

@testset "EOP — bulletin = :B reads the final values, and Bulletin A where B has none" begin
    _with_eop_restored() do
        p2000a, p1980 = _eop_file(_EOP_2000A_ROWS), _eop_file(_EOP_1980_ROWS)

        t = eop_load(p2000a; theory = IAU2006(), bulletin = :B)
        @test eop(IAU2006()) === t
        @test t.x(_JD_1973_01_03)        ≈ 0.141     atol = 1e-12   # bulB/x_pole
        @test t.y(_JD_1973_01_03)        ≈ 0.134     atol = 1e-12
        @test t.Δut1_utc(_JD_1973_01_03) ≈ 0.8044    atol = 1e-12
        @test t.δx(_JD_1973_01_03)       ≈ -18.636   atol = 1e-12
        @test t.δy(_JD_1973_01_03)       ≈ -3.571    atol = 1e-12
        @test t.lod(_JD_1973_01_03)      ≈ 3.5563    atol = 1e-12   # Bulletin B has none

        f = eop_load(p1980; theory = FK5(), bulletin = :B)
        @test f.δΔψ(_JD_1973_01_03) ≈ 0.0 atol = 1e-12

        # A line with no Bulletin B keeps its Bulletin A values.
        rows = replace(_EOP_2000A_ROWS, "final;.141000;.134000;final;.8044000;final;;;-18.636;-3.571" => "final;;;final;;final;;;;")
        p = _eop_file(rows)
        a = eop_load(p; theory = IAU2006(), bulletin = :B)
        @test a.x(_JD_1973_01_03)  ≈ 0.118980 atol = 1e-12
        @test a.δx(_JD_1973_01_03) ≈ -0.751   atol = 1e-12

        @test_throws ArgumentError eop_load(p2000a; theory = IAU2006(), bulletin = :C)
        rm(p2000a; force = true); rm(p1980; force = true); rm(p; force = true)
    end
end

@testset "EOP — set_eop! installs a table by its type" begin
    _with_eop_restored() do
        t1980  = eop_load(_eop_file(_EOP_1980_ROWS);  theory = FK5())
        t2000a = eop_load(_eop_file(_EOP_2000A_ROWS); theory = IAU2006())
        live_2000a = eop(IAU2006())

        @test set_eop!(t1980) === t1980
        @test eop(FK5()) === t1980
        @test eop(IAU2006()) === live_2000a
        @test set_eop!(t2000a) === t2000a
        @test eop(IAU2006()) === t2000a
    end
end

@testset "EOP — an empty slot fills once, whichever task asks first" begin
    _with_eop_restored() do
        AstroUniverse._eop_2000a[] = nothing
        tables = fetch.([Threads.@spawn eop(IAU2006()) for _ in 1:8])
        @test all(t -> t isa EopIau2000A, tables)
        @test all(t -> t === tables[1], tables)
    end
end

# --- The fixed-width finals file ------------------------------------------------------------------
#
# Truth: the same three days as the CSV rows above, exactly as the US Naval Observatory publishes
# them in finals2000A.all and finals.all. Read through the fixed-width reader they must give the
# tables the IERS CSV gives, value for value.

const _FINALS_2000A_ROWS = """
73 1 2 41684.00 I  0.120733 0.009786  0.136966 0.015902  I 0.8084178 0.0002710  0.0000 0.1916  P    -0.766    0.199    -0.720    0.300   .143000   .137000   .8075000   -18.637    -3.667
73 1 3 41685.00 I  0.118980 0.011039  0.135656 0.013616  I 0.8056163 0.0002710  3.5563 0.1916  P    -0.751    0.199    -0.701    0.300   .141000   .134000   .8044000   -18.636    -3.571
73 1 4 41686.00 I  0.117227 0.011039  0.134348 0.013616  I 0.8027895 0.0002710  2.6599 0.1916  P    -0.738    0.199    -0.662    0.300   .139000   .131000   .8012000   -18.669    -3.621
"""

const _FINALS_1980_ROWS = """
73 1 2 41684.00 I  0.120733 0.009786  0.136966 0.015902  I 0.8084178 0.0002710  0.0000 0.1916  P    44.969     .500     2.839     .300   .143000   .137000   .8075000      .000      .000
73 1 3 41685.00 I  0.118980 0.011039  0.135656 0.013616  I 0.8056163 0.0002710  3.5563 0.1916  P    45.005     .500     2.762     .300   .141000   .134000   .8044000      .000      .000
73 1 4 41686.00 I  0.117227 0.011039  0.134348 0.013616  I 0.8027895 0.0002710  2.6599 0.1916  P    45.122     .500     2.851     .300   .139000   .131000   .8012000      .000      .000
"""

function _finals_file(rows)
    path = tempname()
    write(path, rows)
    return path
end

const _EOP_FIELDS = (:x, :y, :Δut1_utc, :lod)

@testset "EOP — the fixed-width file reads as the IERS CSV does" begin
    _with_eop_restored() do
        for (theory, csv_rows, fw_rows, nut) in
                ((IAU2006(), _EOP_2000A_ROWS, _FINALS_2000A_ROWS, (:δx, :δy)),
                 (FK5(),     _EOP_1980_ROWS,  _FINALS_1980_ROWS,  (:δΔψ, :δΔϵ)))
            from_csv = eop_load(_eop_file(csv_rows); theory = theory)
            from_fw  = eop_load(_finals_file(fw_rows); theory = theory)
            for jd in (2441684.5, _JD_1973_01_03, 2441686.5), f in (_EOP_FIELDS..., nut...)
                @test getproperty(from_fw, f)(jd) == getproperty(from_csv, f)(jd)
            end
            # Bulletin B, which the fixed-width file carries in its own columns, reads the same too.
            b_csv = eop_load(_eop_file(csv_rows); theory = theory, bulletin = :B)
            b_fw  = eop_load(_finals_file(fw_rows); theory = theory, bulletin = :B)
            for jd in (2441684.5, _JD_1973_01_03, 2441686.5), f in (_EOP_FIELDS..., nut...)
                @test getproperty(b_fw, f)(jd) == getproperty(b_csv, f)(jd)
            end
        end
    end
end

# Earth orientation at an epoch: the table read at the epoch, plus the installed corrections.
# Truth for the table read is the 1973 rows above at a row's date. Truth for the corrections is
# arithmetic: a model that adds a known amount must show up as that amount, in every theory.

using AstroEpochs
using ForwardDiff

const _T_1973_01_03 = Time("1973-01-03T00:00:00", UTC(), ISOT())
const _JD_TT_1973_01_03 = _JD_1973_01_03 + 44.184 / 86_400          # TAI − UTC was 12 s in January 1973

struct _PoleShift <: AbstractEopCorrection
    x::Float64                                                     # arcsec
end
AstroUniverse.eop_correction(m::_PoleShift, jd_utc, jd_tt) = (x = m.x, y = 0.0, ut1_utc = 0.0, lod = 0.0)

# UT1 shifted by TT − UTC in seconds and LOD by half a millisecond, so the test sees both dates arrive.
struct _BothDates <: AbstractEopCorrection end
AstroUniverse.eop_correction(::_BothDates, jd_utc, jd_tt) =
    (x = 0.0, y = 0.0, ut1_utc = (jd_tt - jd_utc) * 86_400, lod = 0.5)

"""Run `f` with the tables, theory and corrections put back afterwards."""
function _with_corrections_restored(f)
    before = eop_corrections()
    try
        _with_eop_restored(f)
    finally
        set_eop_corrections!(before)
    end
end

@testset "EOP — eop(theory, epoch) reads the table at an epoch" begin
    _with_eop_restored() do
        eop_load(_eop_file(_EOP_2000A_ROWS); theory = IAU2006())
        eop_load(_eop_file(_EOP_1980_ROWS);  theory = FK5())

        v = eop(IAU2006(), _T_1973_01_03)
        @test v isa EopValues
        @test v.theory === IAU2006()
        @test v.x                 ≈ 0.118980  atol = 1e-9      # arcsec
        @test v.y                 ≈ 0.135656  atol = 1e-9
        @test v.ut1_utc           ≈ 0.8056163 atol = 1e-9      # s
        @test v.lod               ≈ 3.5563    atol = 1e-9      # ms
        @test v.celestial_pole[1] ≈ -0.751    atol = 1e-9      # δx, mas
        @test v.celestial_pole[2] ≈ -0.701    atol = 1e-9      # δy

        f = eop(FK5(), _T_1973_01_03)
        @test f.theory === FK5()
        @test f.x == v.x && f.y == v.y && f.ut1_utc == v.ut1_utc && f.lod == v.lod
        @test f.celestial_pole[1] ≈ 45.005 atol = 1e-9         # δΔψ, mas
        @test f.celestial_pole[2] ≈ 2.762  atol = 1e-9         # δΔϵ

        # The active theory, and the two-date form the frame chain calls with.
        set_frame_theory!(FK5())
        @test eop(_T_1973_01_03) == f
        w = eop(IAU2006(), _JD_1973_01_03, _JD_TT_1973_01_03)
        @test w.x ≈ v.x atol = 1e-12
        @test w.ut1_utc ≈ v.ut1_utc atol = 1e-12

        @test occursin("δx, δy = -0.751, -0.701 mas", sprint(show, MIME"text/plain"(), v))
        @test occursin("δΔψ, δΔϵ = 45.005, 2.762 mas", sprint(show, MIME"text/plain"(), f))
        @test occursin("x       = 0.11898 arcsec", sprint(show, v))
    end
end

@testset "EOP — corrections add to the table and reach both theories" begin
    _with_corrections_restored() do
        eop_load(_eop_file(_EOP_2000A_ROWS); theory = IAU2006())
        eop_load(_eop_file(_EOP_1980_ROWS);  theory = FK5())
        set_eop_corrections!([])
        plain = eop(IAU2006(), _T_1973_01_03)

        installed = set_eop_corrections!([_PoleShift(0.001), _PoleShift(0.002), _BothDates()])
        @test installed isa Vector{AbstractEopCorrection}
        @test length(installed) == 3
        @test eop_corrections() === installed

        v = eop(IAU2006(), _T_1973_01_03)
        @test v.x       ≈ plain.x + 0.003        atol = 1e-12
        @test v.y       == plain.y
        # A single Float64 Julian date resolves 40 µs, which is what the dates handed to a model carry.
        @test v.ut1_utc ≈ plain.ut1_utc + 44.184 atol = 1e-4
        @test v.lod     ≈ plain.lod + 0.5        atol = 1e-12
        @test v.celestial_pole == plain.celestial_pole
        @test eop(FK5(), _T_1973_01_03).x ≈ plain.x + 0.003 atol = 1e-12

        set_eop_corrections!([])
        @test isempty(eop_corrections())
        @test eop(IAU2006(), _T_1973_01_03) == plain
    end
end

@testset "EOP — a dual date differentiates through the table and the corrections" begin
    _with_corrections_restored() do
        eop_load(_eop_file(_EOP_2000A_ROWS); theory = IAU2006())
        set_eop_corrections!([_BothDates()])
        jd, jd_tt = _JD_1973_01_03 + 0.25, _JD_TT_1973_01_03 + 0.25
        ut1(jd_utc) = eop(IAU2006(), jd_utc, jd_tt).ut1_utc
        d = eop(IAU2006(), ForwardDiff.Dual(jd, 1.0), jd_tt)
        @test d isa EopValues{IAU2006, <:ForwardDiff.Dual}
        @test ForwardDiff.value(d.ut1_utc) ≈ ut1(jd) atol = 1e-12
        # The table is linear between its daily knots, so a wide step loses nothing and keeps the
        # 40 µs resolution of a single Float64 date out of the difference.
        h = 1e-2
        @test ForwardDiff.partials(d.ut1_utc, 1) ≈ (ut1(jd + h) - ut1(jd - h)) / 2h rtol = 1e-6
    end
end

@testset "EOP — fixed-width predictions, signs against the flag, and date-only rows" begin
    # Truth: a predicted row and the last row of finals2000A.all as published, read directly. The
    # UT1-UTC value fills its column up to the flag, and a date-only row carries no values.
    rows = """
271029 61707.00 P  0.179319 0.017516  0.276208 0.028383  P-0.2081431 0.0254096
271218 61757.00
"""
    csv = tempname()
    AstroUniverse._finals_to_csv(_finals_file(rows), csv, IAU2006())
    lines = readlines(csv)
    @test lines[1] == AstroUniverse._IERS_CSV_HEADER
    f = split(lines[2], ';')
    @test length(f) == 37
    @test f[1:5] == ["61707", "2027", "10", "29", "prediction"]
    @test parse(Float64, f[6]) == 0.179319 && parse(Float64, f[8]) == 0.276208
    @test parse(Float64, f[15]) == -0.2081431 && parse(Float64, f[16]) == 0.0254096
    @test all(isempty, f[17:27])                          # no LOD and no nutation offsets
    g = split(lines[3], ';')
    @test g[1:2] == ["61757", "2027"] && all(isempty, g[5:end])
    @test_throws ArgumentError AstroUniverse._finals_to_csv(_finals_file("<html>Request Rejected</html>"),
                                                            tempname(), IAU2006())
end

# --- Sources and the cache --------------------------------------------------------------------------
#
# Truth: what each source holds, served from local files. A source that fails, serves something
# that is not EOP data, or serves data ending before today is passed over for the next; the cache
# changes only when a source succeeds; and when none does, a cached copy is used with a warning and
# no copy is an error naming every source. "Today" is set to a day the fixture covers.

_missing_url() = _url(joinpath(tempdir(), "no_such_eop_source_$(time_ns())"))

# `f(cached)` with `sources` for every theory and a cache of its own, where `cached` is the path the
# IAU 2000A copy is kept at. The file's fixture settings are put back afterwards.
function _with_sources(f, sources; today = 41685.0)
    cache = mktempdir()
    AU._EOP_SOURCES[]   = _ -> sources
    AU._EOP_CACHE_DIR[] = cache
    AU._EOP_TODAY_MJD[] = today
    try
        f(joinpath(cache, AU._eop_cache_name(IAU2006())))
    finally
        _use_fixture_sources!()
    end
end

_good_fw()  = AU._EopSource(_url(_finals_file(_FINALS_2000A_ROWS)), :finals)
_good_csv() = AU._EopSource(_url(_eop_file(_EOP_2000A_ROWS)), :csv)
_html()     = AU._EopSource(_url(_finals_file("<html><body>Request Rejected</body></html>")), :csv)
_reads_as_table(path) = AU.read_iers_eop(path, Val(:IAU2000A)) isa EopIau2000A

@testset "EOP — each source is tried in turn until one serves a usable series" begin
    _with_eop_restored() do
        sources = [AU._EopSource(_missing_url(), :csv), _html(), _good_fw()]
        _with_sources(sources) do cached
            t = eop_refresh!(; theory = IAU2006())
            @test isfile(cached) && AU._is_iers_csv(cached)        # stored in the CSV layout
            @test t.δx(_JD_1973_01_03) == -0.751
            @test eop(IAU2006()) === t
        end
        # The CSV source first, when it answers.
        _with_sources([_good_csv(), _good_fw()]) do cached
            eop_refresh!(; theory = IAU2006())
            @test read(cached, String) == AU._IERS_CSV_HEADER * "\n" * _EOP_2000A_ROWS
        end
    end
end

@testset "EOP — a series that ends before today is passed over" begin
    _with_eop_restored() do
        _with_sources([_good_csv()]; today = 41700.0) do cached
            err = try eop_refresh!(; theory = IAU2006()); nothing catch e; e end
            @test err isa ErrorException
            @test occursin("before today", sprint(showerror, err))
            @test !isfile(cached)
        end
    end
end

@testset "EOP — when every source fails a cached copy is used, and a refresh leaves it alone" begin
    _with_eop_restored() do
        bad = [AU._EopSource(_missing_url(), :csv), AU._EopSource(_missing_url(), :finals)]
        _with_sources([_good_fw()]) do cached
            eop_refresh!(; theory = IAU2006())               # a good copy in the cache
            before = read(cached)
            AU._EOP_SOURCES[] = _ -> bad

            # A fresh copy is read without asking any source, so nothing fails and nothing warns.
            path = @test_logs AU._eop_file(IAU2006())
            @test path == cached

            # A week-old copy: every source fails, so the copy is used, with a warning.
            AU._EOP_AGE_LIMIT[] = 0.0
            path = @test_logs (:warn, r"using the cached copy") AU._eop_file(IAU2006())
            @test path == cached && read(cached) == before

            # A refresh that was asked for throws, and the copy is unchanged.
            err = try eop_refresh!(; theory = IAU2006()); nothing catch e; e end
            @test err isa ErrorException && occursin("is unchanged", sprint(showerror, err))
            @test read(cached) == before
        end

        # No copy at all: an error naming every source and how to load a file by hand.
        _with_sources(bad) do cached
            err = try AU._eop_file(IAU2006()); nothing catch e; e end
            msg = sprint(showerror, err)
            @test err isa ErrorException
            @test all(s -> occursin(s.url, msg), bad)
            @test occursin("eop_load", msg)
            @test !isfile(cached)
        end
    end
end

@testset "EOP — a replacement that fails leaves the cached copy whole" begin
    # Truth: the cached bytes. The rename that installs a download is made to fail after the
    # download has been validated, the last point at which the cache could be lost.
    _with_eop_restored() do
        _with_sources([_good_fw()]) do cached
            eop_refresh!(; theory = IAU2006())
            before = read(cached)
            AU._EOP_REPLACE[] = (src, dst) -> error("simulated rename failure")
            err = try eop_refresh!(; theory = IAU2006()); nothing catch e; e end
            @test err isa ErrorException
            @test occursin("simulated rename failure", sprint(showerror, err))
            @test read(cached) == before
            @test readdir(dirname(cached)) == [basename(cached)]        # no temporary files left
        end
    end
end

@testset "EOP — the rename itself leaves the destination whole when it fails" begin
    # Truth: the destination's bytes. Renaming a file that does not exist fails in the operating
    # system's rename, which must not have removed the file it would have replaced.
    dir = mktempdir()
    dst = joinpath(dir, "finals2000A.all.csv")
    write(dst, "the cached copy")
    @test_throws Base.IOError AU._replace_file(joinpath(dir, "no_such_download"), dst)
    @test read(dst, String) == "the cached copy"

    src = joinpath(dir, "download")
    write(src, "a new copy")
    @test AU._replace_file(src, dst) == dst
    @test read(dst, String) == "a new copy" && !isfile(src)
end

@testset "EOP — refreshes in several tasks at once each install a whole file" begin
    _with_eop_restored() do
        _with_sources([_good_fw()]) do cached
            tables = fetch.([Threads.@spawn eop_refresh!(; theory = IAU2006()) for _ in 1:8])
            @test all(t -> t isa EopIau2000A && t.δx(_JD_1973_01_03) == -0.751, tables)
            @test _reads_as_table(cached)
            @test readdir(dirname(cached)) == [basename(cached)]
        end
    end
end


# The fixture rows with one field of one row replaced.
function _eop_rows_with(row, col, value)
    rows = split(strip(_EOP_2000A_ROWS), '\n')
    f = split(rows[row], ';')
    f[col] = value
    rows[row] = join(f, ';')
    return join(rows, '\n') * "\n"
end

@testset "EOP — every value present must be a finite number, and coverage reads either layout" begin
    # Truth: the fixture's last UT1-UTC date, 41686, in the 37- and 33-column layouts, where
    # UT1-UTC is column 15 and column 11; and fields made infinite or unreadable, which the table
    # reader would accept as a value or as a missing one.
    @test AU._check_eop_csv(_eop_file(_EOP_2000A_ROWS)) == 41686.0
    for (row, col, value) in ((2, 6, "Inf"), (2, 15, "-Inf"), (3, 15, "garbled"), (1, 24, "NaN"))
        err = try AU._check_eop_csv(_eop_file(_eop_rows_with(row, col, value))); nothing catch e; e end
        @test err isa ErrorException && occursin("not a finite number", sprint(showerror, err))
    end

    drop = i -> !(10 <= i <= 13)                         # the polar-motion rate columns
    narrow(line) = join([x for (i, x) in enumerate(split(line, ';')) if drop(i)], ';')
    rows = split(strip(_EOP_2000A_ROWS), '\n')
    path33 = tempname() * ".csv"
    write(path33, narrow(_EOP_HEADER) * "\n" * join(narrow.(rows), '\n') * "\n")
    @test count(==(';'), readline(path33)) + 1 == 33
    @test AU._check_eop_csv(path33) == 41686.0
    @test AU.read_iers_eop(path33, Val(:IAU2000A)).Δut1_utc(_JD_1973_01_03) ≈ 0.8056163 atol = 1e-12
end

@testset "EOP — a download with an infinite value is refused and the cache kept" begin
    # Truth: the cached bytes. The table reader accepts Inf, so this is the check that stops a
    # corrupt polar motion from reaching a frame transformation.
    _with_eop_restored() do
        _with_sources([_good_fw()]) do cached
            eop_refresh!(; theory = IAU2006())
            before = read(cached)
            corrupt = AU._EopSource(_url(_eop_file(_eop_rows_with(2, 6, "Inf"))), :csv)
            AU._EOP_SOURCES[] = _ -> [corrupt]
            err = try eop_refresh!(; theory = IAU2006()); nothing catch e; e end
            @test err isa ErrorException && occursin("not a finite number", sprint(showerror, err))
            @test read(cached) == before
        end
    end
end

@testset "EOP — with no cache of its own, the older SatelliteToolboxTransformations copy is read" begin
    # Truth: the copy placed where SatelliteToolboxTransformations kept it. It is read where it is,
    # with a warning, and never by a refresh that was asked for.
    _with_eop_restored() do
        bad = [AU._EopSource(_missing_url(), :csv)]
        _with_sources(bad) do cached
            legacy = joinpath(AU._EOP_LEGACY_DIR[], AU._eop_cache_name(IAU2006()))
            cp(_eop_file(_EOP_2000A_ROWS), legacy)
            path = @test_logs (:warn, r"SatelliteToolboxTransformations cached earlier") AU._eop_file(IAU2006())
            @test path == legacy
            @test !isfile(cached)                                  # not copied into the cache
            err = try eop_refresh!(; theory = IAU2006()); nothing catch e; e end
            @test err isa ErrorException
        end
    end
end

# --- The published sources ------------------------------------------------------------------------
#
# Only when asked for, with EPICYCLE_NETWORK_TESTS=1. These are third-party servers, so a test that
# needs them fails for reasons that are not Epicycle's; asked for, it says which server is serving.
# Truth: each source, on its own, yields a table that reaches today.

if get(ENV, "EPICYCLE_NETWORK_TESTS", "") == "1"
    @testset "EOP — each published source serves today's series" begin
        for theory in (IAU2006(), FK5()), source in AU._eop_sources(theory)
            @testset "$(source.url)" begin
                _with_sources([source]; today = nothing) do _
                    path = joinpath(AU._EOP_CACHE_DIR[], AU._eop_cache_name(theory))
                    @test AU._download_eop!(path, source, theory) == path
                end
            end
        end
    end
end

# Put back the settings this file found.
AU._EOP_SOURCES[], AU._EOP_CACHE_DIR[], AU._EOP_LEGACY_DIR[], AU._EOP_AGE_LIMIT[],
    AU._EOP_TODAY_MJD[], AU._EOP_REPLACE[] = _EOP_HOOKS_FOUND

nothing
