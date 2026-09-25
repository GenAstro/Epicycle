# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

using AstroFrames
using Test
using InteractiveUtils 
using AstroUniverse  # for earth object

cs = CoordinateSystem(earth, ICRF())
@testset "Construction and Field Access" begin
    @test cs.origin === earth           
    @test cs.axes isa ICRF        

    origin_field = cs.origin  
    axes_field = cs.axes       
    
    @test origin_field === earth           
    @test axes_field isa ICRF         
    @test typeof(axes_field) === ICRF
end

@testset "Field Access for inlined field" begin
    cs = CoordinateSystem(earth, ICRF())
    
    # Force dynamic access that can't be optimized away
    field_name = :axes
    axes_value = getfield(cs, field_name)  # Dynamic field access
    @test axes_value isa ICRF
    
    # Or use reflection to force access
    @test hasfield(typeof(cs), :axes)
    @test getfield(cs, :axes) isa AbstractAxes
end

@testset "Display Methods" begin
    # Test the show methods that are currently uncovered
    io = IOBuffer()
    show(io, MIME("text/plain"), cs)   
    output = String(take!(io))
    
    @test contains(output, "CoordinateSystem:")
    @test contains(output, "origin =")
    @test contains(output, "axes   =")
    
    # Test the delegation method
    io2 = IOBuffer()
    show(io2, cs)                      
    output2 = String(take!(io2))
    @test output == output2            
end

@testset "Helper Function Coverage" begin
    # Test _maybe_get with different scenarios
    cs = CoordinateSystem(earth, VNB())
    
    # This should exercise the _maybe_get function 
    io = IOBuffer()
    show(io, MIME("text/plain"), cs) 
    output = String(take!(io))
    @test !isempty(output)
end

@testset "Abstract Types" begin
    @test AbstractCoordinateSystem isa Type
    @test AbstractAxes isa Type
    @test isabstracttype(AbstractCoordinateSystem)
    @test isabstracttype(AbstractAxes)
end


@testset "Type Construction" begin
    @test ICRF() isa ICRF
    @test MJ2000Eq() isa MJ2000Eq
    @test VNB() isa VNB
end

@testset "Type Hierarchy" begin
    @test ICRF <: AbstractAxes
    @test MJ2000Eq <: AbstractAxes
    @test VNB <: AbstractAxes
end

@testset "Type Uniqueness" begin
    # Test that axes types are singletons (empty structs)
    @test ICRF() === ICRF()
    @test MJ2000Eq() === MJ2000Eq()
    @test VNB() === VNB()
end

@testset "Exported Types" begin
    @test @isdefined ICRF
    @test @isdefined MJ2000Eq
    @test @isdefined VNB
    @test @isdefined AbstractAxes
end


@testset "CoordinateSystem" begin        
    @testset "Export" begin
        @test @isdefined CoordinateSystem
        @test @isdefined AbstractCoordinateSystem
    end
    
    @testset "Type Definition" begin
        # Test that CoordinateSystem is properly defined
        @test CoordinateSystem <: AbstractCoordinateSystem
        @test isabstracttype(AbstractCoordinateSystem)
        @test !isabstracttype(CoordinateSystem)
    end
end

@testset "Type Discovery" begin
    @testset "subtypes function" begin
        # Test that users can discover available axes types
        axes_subtypes = subtypes(AbstractAxes)
        
        # Every exported axes type must be discoverable. Asserting a fixed
        # count instead would break on every new axes type, which is what
        # happened here: this read `== 4` from the era when there were four.
        for A in (ICRF, GCRF, CIRS, TIRS, ITRF, TEME,
                  MJ2000Eq, MJ2000Ec, MODEq, TODEq, MODEc, TODEc, PEF,
                  MoonPA, MoonME, CelestialBodyFixed, RIC, LVLH, VNB, Inertial)
            @test A in axes_subtypes
        end
    end
    
    @testset "supertype relationships" begin
        @test supertype(ICRF) === AbstractAxes
        @test supertype(MJ2000Eq) === AbstractAxes
        @test supertype(VNB) === AbstractAxes
        @test supertype(CoordinateSystem) === AbstractCoordinateSystem
    end
end
    
nothing

include("test_correctness_invariants.jl")
include("test_correctness_translation.jl")
include("test_correctness_orbit_relative.jl")
include("test_correctness_orbit_relative_truth.jl")
include("test_correctness_body_fixed.jl")
include("test_robustness_origin_coupling.jl")
include("test_robustness_eop_range.jl")
include("test_correctness_user_body.jl")
include("test_correctness_coordinate.jl")
include("test_correctness_subject_origin.jl")
include("test_correctness_earth_truth.jl")
include("test_correctness_earth_rates.jl")
include("test_correctness_body_truth.jl")
include("test_correctness_translation_truth.jl")
include("test_correctness_moon_truth.jl")
include("test_correctness_routing_paths.jl")
