# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# The mass field is private and `total_mass` is the only supported read. These cover the access rules and the notice
# memory; the burn behaviour that uses them is tested in AstroManeuvers.

using Test
using ForwardDiff
using Logging

using EpicycleBase
using AstroStates
using AstroEpochs
using AstroFrames
using AstroUniverse
using AstroModels

_demo_sc(; mass = 1000.0) = Spacecraft(
    state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]),
    time  = Time("2015-09-21T12:23:12", TAI(), ISOT()),
    mass  = mass,
    name  = "SC-mass")

@testset "mass field is private" begin
    sc = _demo_sc()

    # Reading the field would hand back a number whose meaning changes when the mass
    # model lands, so it raises instead.
    @test_throws ErrorException sc.mass
    @test_throws ErrorException sc.mass = 500.0

    # The message has to name the replacement, or the error is just an obstacle.
    err = try; sc.mass; catch e; sprint(showerror, e); end
    @test occursin("total_mass", err)

    # Every other field is untouched by the getproperty hook.
    @test sc.name == "SC-mass"
    @test sc.save_history === true
    @test sc.state isa OrbitState

    @test total_mass(sc) == 1000.0
end

@testset "mass has no setter through the param tag" begin
    sc = _demo_sc()

    # `Mass` still reads. Setting a total would require deciding which tank it came
    # from, so the method is gone rather than guessing.
    @test get_field(sc, Mass()) == 1000.0
    @test_throws MethodError set_field!(sc, Mass(), 500.0)
end

@testset "a notice fires once per spacecraft" begin
    sc = _demo_sc()

    @test_logs (:warn, "first") AstroModels.notify_once!(sc, :demo, "first")
    @test_logs AstroModels.notify_once!(sc, :demo, "first")          # silent second time

    # Distinct tags are independent.
    @test_logs (:warn, "other") AstroModels.notify_once!(sc, :other, "other")

    @test AstroModels.notify_once!(sc, :fresh, "fresh") === true
    @test AstroModels.notify_once!(sc, :fresh, "fresh") === false
end

@testset "notice memory survives promotion and deepcopy" begin
    # This is the invariant that matters under a solver. Every AD pass rebuilds
    # the spacecraft, so a copy that resets `notified` turns "once per spacecraft" into
    # once per iterate — invisible outside an AD or optimizer context.
    sc = _demo_sc()
    dual = Base.promote(sc, ForwardDiff.Dual{Nothing,Float64,1})

    @test_logs (:warn, "once") AstroModels.notify_once!(dual, :shared, "once")
    @test_logs AstroModels.notify_once!(sc, :shared, "once")         # original sees it too

    # And the other way round: notify the original, then promote repeatedly.
    sc2 = _demo_sc()
    @test_logs (:warn, "orig") AstroModels.notify_once!(sc2, :once_only, "orig")
    for _ in 1:3
        iterate_sc = Base.promote(sc2, ForwardDiff.Dual{Nothing,Float64,1})
        @test_logs AstroModels.notify_once!(iterate_sc, :once_only, "orig")
    end

    # deepcopy is an independent object, so it copies rather than shares — but the
    # contents carry over, so it does not warn again.
    sc3 = _demo_sc()
    AstroModels.notify_once!(sc3, :carried, "carried")
    d = deepcopy(sc3)
    @test_logs AstroModels.notify_once!(d, :carried, "carried")
    @test getfield(d, :notified) !== getfield(sc3, :notified)
end

@testset "total_mass is differentiable" begin
    # Drag and SRP accelerations divide by mass, and a solver may differentiate with respect
    # to it, so the accessor has to pass a Dual through unchanged.
    g = ForwardDiff.derivative(m -> total_mass(_demo_sc(mass = m)), 1000.0)
    @test g == 1.0

    sc = _demo_sc()
    dual = Base.promote(sc, ForwardDiff.Dual{Nothing,Float64,1})
    @test total_mass(dual) isa ForwardDiff.Dual
    @test ForwardDiff.value(total_mass(dual)) ≈ total_mass(sc)
end
