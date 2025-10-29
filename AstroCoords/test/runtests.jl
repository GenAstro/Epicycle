using AstroCoords
using Test

@testset "AstroCoords.jl" begin
    
    @testset "Abstract Types" begin
        @test AbstractCoordinateSystem isa Type
        @test AbstractAxes isa Type
        @test isabstracttype(AbstractCoordinateSystem)
        @test isabstracttype(AbstractAxes)
    end
    
    @testset "Axes Types" begin
        @testset "Type Construction" begin
            @test ICRFAxes() isa ICRFAxes
            @test MJ2000Axes() isa MJ2000Axes
            @test VNB() isa VNB
            @test Inertial() isa Inertial
        end
        
        @testset "Type Hierarchy" begin
            @test ICRFAxes <: AbstractAxes
            @test MJ2000Axes <: AbstractAxes
            @test VNB <: AbstractAxes
            @test Inertial <: AbstractAxes
        end
        
        @testset "Type Uniqueness" begin
            # Test that axes types are singletons (empty structs)
            @test ICRFAxes() === ICRFAxes()
            @test MJ2000Axes() === MJ2000Axes()
            @test VNB() === VNB()
            @test Inertial() === Inertial()
        end
        
        @testset "Exported Types" begin
            @test @isdefined ICRFAxes
            @test @isdefined MJ2000Axes
            @test @isdefined VNB
            @test @isdefined Inertial
            @test @isdefined AbstractAxes
        end
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
            
            @test ICRFAxes in axes_subtypes
            @test MJ2000Axes in axes_subtypes
            @test VNB in axes_subtypes
            @test Inertial in axes_subtypes
            @test length(axes_subtypes) == 4
        end
        
        @testset "supertype relationships" begin
            @test supertype(ICRFAxes) === AbstractAxes
            @test supertype(MJ2000Axes) === AbstractAxes
            @test supertype(VNB) === AbstractAxes
            @test supertype(Inertial) === AbstractAxes
            @test supertype(CoordinateSystem) === AbstractCoordinateSystem
        end
    end
    
end
nothing