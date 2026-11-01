# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

using AstroFrames, AstroEpochs, ForwardDiff, Test

# AF-15, AF-16. Truth: the reference geometry is undefined when r × v = 0;
# malformed arguments must name the field, while numerical failures return NaNs.
@testset "Undefined orbit-relative geometry [AF-15]" begin
    cases = ([0., 0, 0, 0, 7.5, 0], [7000., 0, 0, 0, 0, 0],
             [7000., 0, 0, 1., 0, 0])
    for axes in (RIC(), LVLH(), VNB()), state in cases, dual in (false, true)
        input = dual ? ForwardDiff.Dual.(state, 1.) : state
        M = @test_logs (:warn, r"zero angular momentum") axes_rotation(
            ICRF(), axes, 0., (; reference_state = input))
        @test all(isnan, M[1:3, 1:3])
        @test all(isnan, M[4:6, 1:3])
        @test all(isnan, M[4:6, 4:6])
        dual && @test eltype(M) <: ForwardDiff.Dual
    end
    for axes in (RIC(), LVLH(), VNB())
        M = @test_logs axes_rotation(ICRF(), axes, 0.,
            (; reference_state = [NaN, 0, 0, 0, 7.5, 0]))
        @test all(isnan, M[1:3, 1:3])
        M = @test_logs (:warn, r"non-finite") axes_rotation(ICRF(), axes, 0.,
            (; reference_state = [Inf, 0, 0, 0, 7.5, 0]))
        @test all(isnan, M[1:3, 1:3])
    end
end

@testset "Malformed reference data [AF-16.1]" begin
    good = [7000., 0, 0, 0, 7.5, 0]
    for axes in (RIC(), LVLH(), VNB()), reverse in (false, true)
        source, target = reverse ? (axes, ICRF()) : (ICRF(), axes)
        for (params, field) in (( (; reference_state = good[1:5]), "reference_state"),
                               ( (; reference_state = ["a", "b", "c", "d", "e", "f"]), "reference_state"),
                               ( (; reference_state = good, reference_accel = [0., 0.]), "reference_accel"),
                               ( (; reference_state = good, reference_accel = ["a", "b", "c"]), "reference_accel"))
            err = try axes_rotation(source, target, 0., params); nothing catch e; e end
            @test err isa ArgumentError
            @test occursin(field, sprint(showerror, err))
        end
    end
end

@testset "Missing and nonnumeric data name their fields [AF-16.1]" begin
    good = [7000., 0, 0, 0, 7.5, 0]
    for axes in (RIC(), LVLH(), VNB()), bad in (nothing, missing, "wrong")
        for (params, field) in (((; reference_state = bad), "reference_state"),
                               ((; reference_state = good, reference_accel = bad), "reference_accel"))
            err = try axes_rotation(ICRF(), axes, 0., params); nothing catch e; e end
            @test err isa ArgumentError
            @test occursin(field, sprint(showerror, err))
        end
    end
end

@testset "Overflow in reference acceleration warns [AF-15]" begin
    for axes in (RIC(), LVLH(), VNB())
        M = @test_logs (:warn, r"overflow") axes_rotation(ICRF(), axes, 0.,
            (; reference_state=[7000.,0,0,0,7.5,0], reference_accel=[0.,0,1e308]))
        @test all(isnan, M[1:3,1:3])
    end
end
