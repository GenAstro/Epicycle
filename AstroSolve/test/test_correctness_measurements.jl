# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Measurement types — the data layer an estimator reads.
#
# There is no numerical algorithm here, so the truth is structural: a path records participants
# in order and rejects a participant relaying to itself, noise draws are reproducible from a
# seeded generator and have the declared variance, and the type parameters stay concrete so the
# estimator's inner loop is not dispatching dynamically.
#
# Measurements.jl was 0 of 24 coverable lines.

using AstroSolve
using Random
using Statistics
using Test

const _M = AstroSolve.Measurements

# A participant is anything carrying a `name`; the file path accepts bare strings too.
struct _Station
    name::String
end

@testset "Measurements — signal path records participants in order" begin
    gs  = _Station("DSS-14")
    sat = _Station("Explorer")

    # The canonical two-way pattern: out and back, so the same participant appears twice
    # without being adjacent to itself.
    p = _M.SignalPath(gs, sat, gs)
    @test length(p) == 3
    @test _M.participant_names(p) == ("DSS-14", "Explorer", "DSS-14")
    @test p[1] === gs
    @test p[2] === sat
    @test collect(p) == [gs, sat, gs]

    # One-way is a legal path, and so is a single participant.
    @test length(_M.SignalPath(gs, sat)) == 2
    @test length(_M.SignalPath(gs)) == 1

    # Bare strings work, which is what the load-from-file path produces when only the name
    # survived.
    ps = _M.SignalPath("DSS-14", "Explorer")
    @test _M.participant_names(ps) == ("DSS-14", "Explorer")

    # A four-leg relay keeps its order rather than deduplicating.
    relay = _M.SignalPath(gs, sat, _Station("Relay"), gs)
    @test _M.participant_names(relay) == ("DSS-14", "Explorer", "Relay", "DSS-14")

    # The length is a type parameter, so it is known at compile time and the estimator's loop
    # over legs does not allocate.
    @test p isa _M.SignalPath{3}
end

@testset "Measurements — noise draws and variance" begin
    # A zero sigma draws exactly zero, which is what noise-free simulation depends on. Not
    # approximately zero: the branch returns the literal.
    quiet = _M.MeasurementNoise(0.0)
    rng = MersenneTwister(1)
    @test _M.draw(quiet, rng) === 0.0
    @test _M.variance(quiet) == 0.0

    # A non-zero sigma gives variance sigma², which is the R the estimator forms its gain from.
    σ = 0.03
    noisy = _M.MeasurementNoise(σ)
    @test _M.variance(noisy) ≈ σ^2 atol = 1e-18

    # Draws are reproducible from a seeded generator, which is what makes a simulated data set
    # a regression fixture rather than a new problem each run.
    a = [_M.draw(noisy, MersenneTwister(42)) for _ in 1:3]
    @test all(a .== a[1])

    # And the sample statistics match the declaration. Loose bounds: this asserts the draw is
    # scaled by sigma, not that randn is correct.
    rng = MersenneTwister(2024)
    samples = [_M.draw(noisy, rng) for _ in 1:20_000]
    @test abs(mean(samples)) < 5σ / sqrt(length(samples)) * 3
    @test isapprox(std(samples), σ; rtol = 0.05)

    # Integer sigma converts rather than making the caller do it.
    @test _M.MeasurementNoise(1).sigma === 1.0
end

@testset "Measurements — observables carry their own noise" begin
    gs  = _Station("DSS-14")
    sat = _Station("Explorer")
    path = _M.SignalPath(gs, sat, gs)

    # The default is noise-free, which is what geometric simulation wants.
    r = _M.TwoWayRange(path)
    @test _M.variance(r.noise) == 0.0
    @test r.path === path

    d = _M.TwoWayDoppler(path; noise = _M.MeasurementNoise(1e-4))
    @test _M.variance(d.noise) ≈ 1e-8 atol = 1e-20

    # Both are AbstractMeasurement, which is what lets a problem hold a mixed list.
    @test r isa _M.AbstractMeasurement
    @test d isa _M.AbstractMeasurement

    # Type parameters stay concrete, so a list of one observable type is concretely typed and
    # the estimator's evaluation loop is statically dispatched.
    @test isconcretetype(typeof(r))
    @test isconcretetype(typeof(d))
    @test typeof(r).parameters[1] === typeof(path)

    # Two observables on the same path are independent objects, so setting noise on one does
    # not reach the other.
    r2 = _M.TwoWayRange(path; noise = _M.MeasurementNoise(0.01))
    @test _M.variance(r.noise) == 0.0
    @test _M.variance(r2.noise) ≈ 1e-4 atol = 1e-18
end

@testset "Measurements — visibility defaults to unconditional" begin
    # The open default is that anything can see anything.
    #
    # 🔴 The docstring says a participant type narrows this by adding a method, and that is not
    # what happens today. AstroModels defines and exports its own `is_visible` for
    # `GroundStation` (ground_station.jl:165) as a separate generic, not a method on this one,
    # so the elevation cutoff never overrides the catch-all below and a station reads as always
    # visible. Loading both packages also makes the bare name ambiguous. This is the is_visible
    # collision deferred in Beta_Release.md; which package owns the generic is an API decision,
    # so these tests pin current behaviour rather than assert the promise.
    @test _M.is_visible(_Station("DSS-14"), [7000.0, 0.0, 0.0], 0.0) === true
    @test _M.is_visible("bare name", [0.0, 0.0, 0.0], 0.0) === true
    @test _M.is_visible(nothing, nothing, nothing) === true
end

@testset "Measurements — input validation" begin
    # The exception type and its message are both checked.
    @test_throws ArgumentError _M.SignalPath()
    empty_msg = try; _M.SignalPath(); catch e; sprint(showerror, e); end
    @test occursin("at least one participant", empty_msg)

    # A participant cannot relay to itself, but may appear again later — the two-way path is
    # exactly that, so the check has to be on adjacency rather than on repetition.
    gs = _Station("DSS-14")
    @test_throws ArgumentError _M.SignalPath(gs, gs)
    @test_throws ArgumentError _M.SignalPath(gs, _Station("Explorer"), gs, gs)
    adj_msg = try; _M.SignalPath(gs, gs); catch e; sprint(showerror, e); end
    @test occursin("DSS-14", adj_msg)          # names the offender
    @test occursin("relay to itself", adj_msg) # and says why it is wrong

    # Non-adjacent repetition is legal and must not raise, or two-way range is unexpressible.
    @test _M.SignalPath(gs, _Station("Explorer"), gs) isa _M.SignalPath

    # A standard deviation cannot be negative.
    @test_throws DomainError _M.MeasurementNoise(-1.0)
    @test_throws DomainError _M.MeasurementNoise(-1e-12)
    σ_msg = try; _M.MeasurementNoise(-0.5); catch e; sprint(showerror, e); end
    @test occursin("non-negative", σ_msg)
    @test occursin("-0.5", σ_msg)

    # Zero is allowed, and is the documented noise-free case rather than an edge case.
    @test _M.MeasurementNoise(0.0).sigma == 0.0
end
