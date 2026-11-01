# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

using AstroFrames, AstroUniverse, AstroEpochs, AstroStates, ForwardDiff, LinearAlgebra, Test

# AF-12. Truth: the same instant, converted by AstroEpochs to each supported scale.
# The tolerance permits the final round-off in time conversion, well below a microarcsecond.
@testset "Equivalent instants in every time scale [AF-12]" begin
    saved = frame_theory()
    try
        t = Time("2020-03-15T06:30:00", UTC(), ISOT())
        for theory in (FK5(), IAU2006())
            set_frame_theory!(theory)
            for pair in ((ICRF(), ITRF()), (ICRF(), CelestialBodyFixed(mars)))
                truth = axes_rotation(pair..., t)
                for scale in (:tai, :tt, :tdb, :tcg, :tcb, :utc)
                    @test axes_rotation(pair..., getproperty(t, scale)) ≈ truth atol = 2e-15 rtol = 0
                end
            end
            truth = origin_translation(earth, moon, ICRF(), t)
            for scale in (:tai, :tt, :tdb, :tcg, :tcb, :utc)
                @test origin_translation(earth, moon, ICRF(), getproperty(t, scale)) ≈ truth atol = 1e-9 rtol = 0
            end
        end
    finally
        set_frame_theory!(saved)
    end
end

# AF-10, AF-13. Truth: a linear rotation has Jacobian equal to its operator.
@testset "Coordinate state Jacobian [AF-10, AF-13]" begin
    t = Time(2458849., 0.5, TDB(), JD())
    x = [7000., 200., 30., -0.2, 7.5, 0.4]
    f(s) = to_vector(Coordinate(Coordinate(s, EarthICRF, t), EarthFixed).state)
    @test ForwardDiff.jacobian(f, x) ≈ axes_rotation(ICRF(), ITRF(), t) atol = 2e-15 rtol = 0
end

# AF-8, AF-12. Truth: Earth turns continuously through equal SI-second intervals.
# Use TAI arithmetic so a leap-second day's UTC fraction cannot change the interval.
@testset "Earth orientation through a leap second [AF-8, AF-12]" begin
    saved = frame_theory()
    try
        for theory in (FK5(), IAU2006()), date in ("2012-06-30", "2015-06-30", "2016-12-31")
            set_frame_theory!(theory)
            base = Time(date * "T23:59:59", UTC(), ISOT()).tai
            for second in (0., 0.5, 1., 1.5, 2.)
                t1, t2 = base + second/86400, base + (second+0.5)/86400
                A = axes_rotation(ICRF(), ITRF(), t1)[1:3,1:3]
                B = axes_rotation(ICRF(), ITRF(), t2)[1:3,1:3]
                D = B * A'
                angle = asin(norm((D[3,2]-D[2,3], D[1,3]-D[3,1], D[2,1]-D[1,2])) / 2)
                @test angle ≈ 0.5 * 7.292115146706979e-5 rtol = 2e-6
            end
            # Equal instants in different time scales also agree inside the leap second.
            leap = base + 1.5/86400
            truth = axes_rotation(ICRF(), ITRF(), leap)
            for scale in (:tai, :tt, :tdb, :tcg, :tcb, :utc)
                @test axes_rotation(ICRF(), ITRF(), getproperty(leap,scale)) ≈ truth atol=2e-15 rtol=0
            end
        end
    finally
        set_frame_theory!(saved)
    end
end
