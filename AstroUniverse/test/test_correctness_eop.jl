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

"""Run `f` and put the session's theory and both tables back afterwards."""
function _with_eop_restored(f)
    theory = frame_theory()
    t1980, t2000a = eop(FK5()), eop(IAU2006())
    try
        f()
    finally
        set_frame_theory!(theory)
        set_eop!(t1980)
        set_eop!(t2000a)
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

# eop_refresh! downloads even when the cache is fresh, which is its only difference from eop().
# Truth: SatelliteToolboxTransformations writes a `_timestamp` file beside each IERS file it
# downloads, in its own scratch space, so a download moves that file's modification time to now.
# The installed table must also be the one eop() returns, and cover today. Needs the network.
using Scratch: get_scratch!

const _STB_UUID = Base.UUID("6b019ec1-7a1e-4f04-96c7-a9db1ca5514d")   # SatelliteToolboxTransformations

@testset "EOP — eop_refresh! downloads the current IERS series and installs it" begin
    _with_eop_restored() do
        for (theory, dir, file) in ((IAU2006(), "eop_iau2000A", "finals2000A.all.csv"),
                                    (FK5(),     "eop_iau1980",  "finals.all.csv"))
            stamp  = joinpath(get_scratch!(_STB_UUID, dir), file * "_timestamp")
            before = isfile(stamp) ? mtime(stamp) : 0.0

            table = eop_refresh!(; theory = theory)

            @test mtime(stamp) > before
            @test mtime(stamp) > time() - 600                  # downloaded now, not earlier
            @test eop(theory) === table
            jd_today = time() / 86_400 + 2440587.5              # Unix epoch as a Julian date
            @test isfinite(table.Δut1_utc(jd_today))
        end
    end
end

nothing
