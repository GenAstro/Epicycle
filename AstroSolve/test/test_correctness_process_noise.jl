# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Process noise models — analytic truth.
#
# The contract is small and entirely closed form: a diagonal SNC accumulates psd·Δt over a
# step, zero entries are pruned, and the (G, q) pair must reconstruct the intended G·diag(q)·Gᵀ.
# The reconstruction is the check that matters, because that product is what the Thornton time
# update actually consumes.
#
# ProcessNoiseModels.jl was 0 of 14 coverable lines.

using AstroSolve
using LinearAlgebra
using Test

const _PN = AstroSolve.ProcessNoiseModels

@testset "ProcessNoise — NoNoise is the deterministic case" begin
    for n in (1, 3, 6)
        G, q = _PN.discretize(_PN.NoNoise(), 60.0, n)
        @test size(G) == (n, 0)
        @test isempty(q)
        # This is the shape Thornton accepts as "no process noise", so the product it forms is
        # an n × n zero and the covariance is propagated by Φ alone.
        @test size(G * Diagonal(q) * G') == (n, n)
        @test all(G * Diagonal(q) * G' .== 0.0)
    end

    # Independent of Δt, including zero and a long step.
    for Δt in (0.0, 1.0, 86400.0)
        G, q = _PN.discretize(_PN.NoNoise(), Δt, 4)
        @test size(G) == (4, 0)
        @test isempty(q)
    end
end

@testset "ProcessNoise — DiagonalSNC accumulates psd·Δt" begin
    psd = [0.0, 0.0, 0.0, 1e-9, 2e-9, 4e-9]
    snc = _PN.DiagonalSNC(psd)

    Δt = 60.0
    G, q = _PN.discretize(snc, Δt, 6)

    # Zero entries are pruned, so the mapping is narrower than the slot.
    @test size(G) == (6, 3)
    @test length(q) == 3
    @test q ≈ [1e-9, 2e-9, 4e-9] .* Δt atol = 1e-20

    # G selects the elements that carry noise, one channel each.
    @test G[4, 1] == 1.0
    @test G[5, 2] == 1.0
    @test G[6, 3] == 1.0
    @test sum(G) == 3.0                       # exactly one entry per channel
    @test all(G[1:3, :] .== 0.0)              # the deterministic elements receive nothing

    # The product Thornton consumes puts the variance on the right diagonal entries and
    # nowhere else. This is the assertion that would catch a transposed or misindexed G.
    Qfull = G * Diagonal(q) * G'
    @test size(Qfull) == (6, 6)
    @test diag(Qfull) ≈ psd .* Δt atol = 1e-20
    @test all(abs.(Qfull - Diagonal(diag(Qfull))) .< 1e-20)   # no off-diagonal coupling

    # Variance accumulates linearly in the step, which is what "per unit time" means. A model
    # quoting variance per step rather than a PSD would fail this.
    _, q1 = _PN.discretize(snc, 1.0, 6)
    _, q2 = _PN.discretize(snc, 2.0, 6)
    @test q2 ≈ 2 .* q1 atol = 1e-20

    # Zero duration gives zero variance rather than an error, which a filter update at the same
    # epoch as the previous one depends on.
    _, q0 = _PN.discretize(snc, 0.0, 6)
    @test all(q0 .== 0.0)
    @test length(q0) == 3                     # still three channels, just no accumulation
end

@testset "ProcessNoise — an all-zero psd equals NoNoise" begin
    G, q   = _PN.discretize(_PN.DiagonalSNC(zeros(5)), 30.0, 5)
    Gn, qn = _PN.discretize(_PN.NoNoise(), 30.0, 5)
    @test size(G) == size(Gn)
    @test q == qn
end

@testset "ProcessNoise — construction accepts any real vector" begin
    # The convenience constructor converts, so an integer or a range is accepted and stored as
    # Float64 rather than making the caller convert.
    @test _PN.DiagonalSNC([1, 2, 3]).psd == [1.0, 2.0, 3.0]
    @test _PN.DiagonalSNC(1:3).psd       == [1.0, 2.0, 3.0]
    @test _PN.DiagonalSNC([1, 2, 3]).psd isa Vector{Float64}
end

@testset "ProcessNoise — input validation" begin
    # The exception type and its message are both asserted.
    snc = _PN.DiagonalSNC([1e-9, 1e-9])

    @test_throws ArgumentError _PN.discretize(snc, 60.0, 3)
    @test_throws ArgumentError _PN.discretize(snc, 60.0, 1)
    msg = try; _PN.discretize(snc, 60.0, 3); catch e; sprint(showerror, e); end
    @test occursin("3", msg)                  # what was expected
    @test occursin("2", msg)                  # and what arrived

    @test_throws DomainError _PN.discretize(snc, -1.0, 2)
    dmsg = try; _PN.discretize(snc, -60.0, 2); catch e; sprint(showerror, e); end
    @test occursin("negative", dmsg)
    @test occursin("-60", dmsg)
end
