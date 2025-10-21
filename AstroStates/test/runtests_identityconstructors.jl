using Test
using AstroStates




# Reuse μ, tol, states from runtests.jl
@assert isdefined(Main, :μ)
@assert isdefined(Main, :tol)
@assert isdefined(Main, :states)

μ   = Main.μ
tol = Main.tol

include("state_truthdata_elliptic_orbits.jl")
@testset "Identity constructors for state types" begin
    for (name, truth) in states
        Tparam = typeof(truth)                       # e.g., AlternateEquinoctialState{Float64}
        ctor   = getfield(AstroStates, nameof(Tparam))  # UnionAll constructor, e.g., AlternateEquinoctialState
        svec   = to_vector(truth)

        @testset "$(name) identity" begin
            tested = 0

            # (state::SameType, μ::Real)
            if applicable(ctor, truth, μ)
                r1 = ctor(truth, μ)
                @test typeof(r1) === Tparam
                @test isapproxvec_percent(to_vector(r1), to_vector(truth); tol=tol)
                tested += 1
            end

            # (state::SameType)
            if applicable(ctor, truth)
                r2 = ctor(truth)
                @test typeof(r2) === Tparam
                @test isapproxvec_percent(to_vector(r2), to_vector(truth); tol=tol)
                tested += 1
            end

            # (s::Vector{<:Real}, μ::Real)
            if applicable(ctor, svec, μ)
                r3 = ctor(svec, μ)
                @test typeof(r3) === Tparam
                @test isapproxvec_percent(to_vector(r3), svec; tol=tol)
                tested += 1
            end

            # (s::Vector{<:Real})
            if applicable(ctor, svec)
                r4 = ctor(svec)
                @test typeof(r4) === Tparam
                @test isapproxvec_percent(to_vector(r4), svec; tol=tol)
                tested += 1
            end

            @test tested ≥ 1  # ensure at least one identity form exists per type
        end
    end
end