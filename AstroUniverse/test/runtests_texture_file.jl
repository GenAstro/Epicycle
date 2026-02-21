using Test
using AstroUniverse

@testset "CelestialBody texture_file field" begin
    
    @testset "Positional constructor with texture_file" begin
        b = CelestialBody("TestBody", 398600.0, 6378.0, 0.1, 123, "path/to/texture.jpg")
        @test b.texture_file == "path/to/texture.jpg"
        @test b.name == "TestBody"
        @test b.mu == 398600.0
    end
    
    @testset "Positional constructor with default texture_file" begin
        b = CelestialBody("TestBody", 398600.0, 6378.0, 0.1, 123)
        @test b.texture_file == ""
    end
    
    @testset "Keyword constructor with texture_file" begin
        b = CelestialBody(
            name="TestBody", 
            mu=398600.0, 
            equatorial_radius=6378.0, 
            flattening=0.1, 
            naifid=123,
            texture_file="assets/test.jpg"
        )
        @test b.texture_file == "assets/test.jpg"
    end
    
    @testset "Keyword constructor default texture_file" begin
        b = CelestialBody(name="TestBody", mu=398600.0, equatorial_radius=6378.0, flattening=0.1, naifid=123)
        @test b.texture_file == ""
    end
    
    @testset "String type preservation" begin
        # Verify texture_file is always String, not SubString or other
        b = CelestialBody("Test", 1.0, 2.0, 0.0, 1, "texture.jpg")
        @test b.texture_file isa String
        
        # With empty string
        b2 = CelestialBody("Test", 1.0, 2.0, 0.0, 1, "")
        @test b2.texture_file isa String
        @test b2.texture_file == ""
    end
    
    @testset "Show method with texture" begin
        b = CelestialBody("TestBody", 398600.0, 6378.0, 0.1, 123, "path/to/texture.jpg")
        output = sprint(show, MIME"text/plain"(), b)
        
        @test occursin("CelestialBody:", output)
        @test occursin("Texture File       = path/to/texture.jpg", output)
        @test !occursin("(none)", output)
    end
    
    @testset "Show method without texture" begin
        b = CelestialBody("TestBody", 398600.0, 6378.0, 0.1, 123, "")
        output = sprint(show, MIME"text/plain"(), b)
        
        @test occursin("CelestialBody:", output)
        @test occursin("Texture File       = (none)", output)
    end
    
    @testset "Built-in bodies have expected texture files" begin
        # Bodies with textures
        @test !isempty(earth.texture_file)
        @test occursin("EarthTexture.jpg", earth.texture_file)
        @test occursin("data", earth.texture_file)
        
        @test !isempty(moon.texture_file)
        @test occursin("MoonTexture.jpg", moon.texture_file)
        
        @test !isempty(mars.texture_file)
        @test occursin("MarsTexture.jpg", mars.texture_file)
        
        @test !isempty(sun.texture_file)
        @test occursin("SunTexture.jpg", sun.texture_file)
        
        @test !isempty(mercury.texture_file)
        @test occursin("MercuryTexture.jpg", mercury.texture_file)
        
        @test !isempty(venus.texture_file)
        @test occursin("VenusTexture.jpg", venus.texture_file)
        
        @test !isempty(jupiter.texture_file)
        @test occursin("JupiterTexture.jpg", jupiter.texture_file)
        
        @test !isempty(saturn.texture_file)
        @test occursin("SaturnTexture.jpg", saturn.texture_file)
        
        @test !isempty(uranus.texture_file)
        @test occursin("UranusTexture.jpg", uranus.texture_file)
        
        @test !isempty(neptune.texture_file)
        @test occursin("NeptuneTexture.jpg", neptune.texture_file)
        
        # Pluto has no texture
        @test pluto.texture_file == ""
    end
    
    @testset "No validation - accepts any string" begin
        # Non-existent path accepted
        b1 = CelestialBody("Test", 1.0, 2.0, 0.0, 1, "/totally/fake/path/notreal.jpg")
        @test b1.texture_file == "/totally/fake/path/notreal.jpg"
        
    end
    
end

nothing
