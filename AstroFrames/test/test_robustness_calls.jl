# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

using AstroFrames, AstroUniverse, AstroEpochs, AstroStates, ForwardDiff, Test

# AF-16. Truth: the documented argument types and diagnostic contract.
@testset "Malformed public calls [AF-16]" begin
    state = CartesianState([7000., 0, 0, 0, 7.5, 0])
    for (f, field) in ((()->axes_rotation(ICRF(), ITRF(), "2020"), "epoch"),
                      (()->axes_rotation(ICRF(), ITRF(), "2020", (;)), "epoch"),
                      (()->axes_rotation(ICRF(), RIC(), 0., [1.,2.]), "params"),
                      (()->origin_translation(earth, moon, ICRF(), "2020"), "epoch"),
                      (()->origin_translation(earth, moon, ICRF(), "2020", (;)), "epoch"),
                      (()->origin_translation(earth, moon, ICRF(), 0., [1.,2.]), "params"),
                      (()->Coordinate("wrong", EarthICRF, 0.), "state"),
                      (()->Coordinate(state, EarthICRF, "2020"), "epoch"))
        err = try f(); nothing catch e; e end
        @test err isa ArgumentError
        @test occursin(field, sprint(showerror, err))
    end
    for f in (state_of, frame_of, epoch_of)
        err = try f("wrong"); nothing catch e; e end
        @test err isa ArgumentError
        @test occursin("subject interface", sprint(showerror, err))
    end
end

@testset "Lunar data and derivative failures [AF-16.4, AF-16.5]" begin
    for axes in (MoonPA(), MoonME()), reverse in (false, true)
        source, target = reverse ? (axes, ICRF()) : (ICRF(), axes)
        err = try axes_rotation(source, target, ForwardDiff.Dual(2458849.5, 1.)); nothing catch e; e end
        @test err isa ArgumentError
        @test occursin("not differentiable", sprint(showerror, err))
        # Outside the loaded DE440 lunar orientation's 1550–2650 coverage.
        epoch = EpochScales(3000000., 3000000., 3000000.)
        err = try axes_rotation(source, target, epoch); nothing catch e; e end
        @test err isa ArgumentError
        @test occursin("no lunar orientation data covers", sprint(showerror, err))
    end
end

# AF-15. A non-finite date cannot define a time-dependent orientation.
# Static axes do not use the epoch, so their constant rotation remains valid.
@testset "Non-finite epochs [AF-15]" begin
    for x in (NaN, Inf, -Inf), dual in (false, true)
        epoch = dual ? ForwardDiff.Dual(x, 1.) : x
        for call in (() -> axes_rotation(ICRF(), ITRF(), epoch),
                     () -> axes_rotation(ICRF(), ITRF(), epoch, (;)),
                     () -> origin_translation(earth, moon, ICRF(), epoch),
                     () -> body_fixed_rotation(IAU2015(), 499, epoch))
            if isnan(x)
                result = @test_logs call()
            else
                result = @test_logs (:warn, r"epoch must be finite") call()
            end
            @test all(isnan, result)
        end
        @test axes_rotation(ICRF(), MJ2000Eq(), epoch) == axes_rotation(ICRF(), MJ2000Eq(), 0.)
    end
end

@testset "Nonnumeric Cartesian input names the state [AF-16.1]" begin
    for state in (["a", "b", "c", "d", "e", "f"], fill(1+im, 6))
        err = try Coordinate(state, EarthICRF, 0.); nothing catch e; e end
        @test err isa ArgumentError
        @test occursin("state", sprint(showerror, err))
    end
end

@testset "A raw state needs an epoch before conversion parameters [AF-16.1]" begin
    err = try Coordinate([7000.,0,0,0,7.5,0], EarthICRF, (;)); nothing catch e; e end
    @test err isa ArgumentError
    @test occursin("epoch", sprint(showerror, err))
end
