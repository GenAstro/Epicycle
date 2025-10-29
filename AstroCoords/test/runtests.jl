using AstroCoords
using Test
using InteractiveUtils 
using AstroUniverse  # for earth object

cs = CoordinateSystem(earth, ICRFAxes())
@testset "Construction and Field Access" begin
    @test cs.origin === earth           
    @test cs.axes isa ICRFAxes        

    origin_field = cs.origin  
    axes_field = cs.axes       
    
    @test origin_field === earth           
    @test axes_field isa ICRFAxes         
    @test typeof(axes_field) === ICRFAxes
end

@testset "Field Access for inlined field" begin
    cs = CoordinateSystem(earth, ICRFAxes())
    
    # Force dynamic access that can't be optimized away
    field_name = :axes
    axes_value = getfield(cs, field_name)  # Dynamic field access
    @test axes_value isa ICRFAxes
    
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
    
nothing