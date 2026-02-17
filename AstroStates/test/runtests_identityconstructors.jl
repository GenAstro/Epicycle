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

@testset "OrbitState deepcopy and copy" begin
    # Create an OrbitState from a Cartesian state
    state_vec = [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]
    os_original = OrbitState(state_vec, Cartesian())
    
    @testset "deepcopy_internal" begin
        os_deep = deepcopy(os_original)
        
        # Verify the copy has the same values
        @test os_deep.state == os_original.state
        @test os_deep.statetype === os_original.statetype
        
        # Verify it's not aliased (different objects)
        @test os_deep.state !== os_original.state
        
        # Verify modifying copy doesn't affect original
        original_val = os_original.state[1]
        os_deep.state[1] = 9999.0
        @test os_original.state[1] == original_val  # Original unchanged
        @test os_deep.state[1] == 9999.0
    end
    
    @testset "copy" begin
        os_copy = copy(os_original)
        
        # Verify the copy has the same values
        @test os_copy.state == os_original.state
        @test os_copy.statetype === os_original.statetype
        
        # Verify it's not aliased (different objects)
        @test os_copy.state !== os_original.state
        
        # Verify modifying copy doesn't affect original
        original_val = os_original.state[2]
        os_copy.state[2] = 8888.0
        @test os_original.state[2] == original_val  # Original unchanged
        @test os_copy.state[2] == 8888.0
    end
    
    @testset "deepcopy with different state types" begin
        # Test with Keplerian
        ks = KeplerianState(8000.0, 0.05, 0.1, 0.2, 0.3, 0.4)
        os_kep = OrbitState(ks)
        os_kep_copy = deepcopy(os_kep)
        
        @test os_kep_copy.state == os_kep.state
        @test os_kep_copy.statetype === os_kep.statetype
        @test os_kep_copy.state !== os_kep.state
        
        original_val = os_kep.state[1]
        os_kep_copy.state[1] = 9000.0
        @test os_kep.state[1] == original_val  # Original unchanged
    end
end

