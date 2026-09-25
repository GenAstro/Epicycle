# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# =============================================================================
# `history`.
#
# The property that matters is that a walk re-applies each `Calc` to every
# recorded sample rather than evaluating it once — which is why a `Calc` must
# keep its parts, and why a closure cannot be a column.
#
# `report` used to be tested here too. It moved to `EpicycleIO` with the
# function, along with the row-integrity tests that go with it.
# =============================================================================

using Test
using AstroCallbacks
using AstroModels
using AstroStates
using AstroEpochs
using AstroFrames
using AstroUniverse

const _H_EPOCH = Time(2458849.5, 0.0, :tdb, :jd)

"""A spacecraft with a short recorded history, built without propagating."""
function _recorded()
    sat = Spacecraft(CartesianState([7000.0, 0.0, 1300.0, 0.0, 7.35, 1.0]), _H_EPOCH;
                     coord_sys = CoordinateSystem(earth, ICRF()), name = "rec")
    seg = AstroModels.HistorySegment(CoordinateSystem(earth, ICRF()); name = "seg")
    for k in 0:4
        push!(seg.times,  Time(2458849.5 + k * 0.01, 0.0, :tdb, :jd))
        push!(seg.states, CartesianState([7000.0 + 10k, 100.0k, 1300.0, 0.0, 7.35, 1.0]))
    end
    push!(sat.history.segments, seg)
    return sat
end

@testset "a walk re-applies the Calc to every sample" begin
    sat = _recorded()
    eq  = CoordinateSystem(earth, MJ2000Eq())

    x = history(Calc(position_x, sat, eq))
    @test length(x) == 5

    # Each entry is that sample's value, not the spacecraft's current one
    # repeated — which is exactly what a closure would have produced.
    @test length(unique(x)) == 5
    @test !all(==(position_x(sat, eq)), x)
end

@testset "one walk serves several columns, in several frames" begin
    sat = _recorded()
    eq  = CoordinateSystem(earth, MJ2000Eq())
    ec  = CoordinateSystem(earth, MJ2000Ec())

    t, r_eq, r_ec = history(Calc(epoch, sat),
                            Calc(position_vector, sat, eq),
                            Calc(position_vector, sat, ec))

    @test length(t) == length(r_eq) == length(r_ec) == 5
    @test t[1] isa Time

    # The ecliptic column is genuinely a different frame, not a copy. The
    # obliquity rotation moves Y and Z by hundreds of km at this altitude.
    @test maximum(abs.(r_eq[1] .- r_ec[1])) > 100.0

    # X is the rotation axis, so it is untouched.
    @test r_eq[1][1] ≈ r_ec[1][1] rtol = 1e-12
end

@testset "a single Calc returns a bare vector" begin
    sat = _recorded()
    @test history(Calc(position_x, sat, CoordinateSystem(earth, MJ2000Eq()))) isa Vector
end

@testset "different subjects are different series" begin
    # Two recordings with their own sample times, so two walks and two columns
    # that need not be the same length.
    a, b = _recorded(), _recorded()
    push!(b.history.segments[1].times,  Time(2458849.56, 0.0, :tdb, :jd))
    push!(b.history.segments[1].states, CartesianState([7060.0, 500.0, 1300.0, 0.0, 7.35, 1.0]))

    eq = CoordinateSystem(earth, MJ2000Eq())
    xa, xb = history(Calc(position_x, a, eq), Calc(position_x, b, eq))
    @test length(xa) == 5
    @test length(xb) == 6
end

@testset "nothing recorded says so" begin
    bare = Spacecraft(CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]), _H_EPOCH;
                      coord_sys = CoordinateSystem(earth, ICRF()), name = "bare")
    e = try; history(Calc(position_x, bare)); nothing; catch e; e; end
    @test e isa ArgumentError
    @test occursin("no recorded segments", e.msg)
    @test occursin("bare", e.msg)
end
