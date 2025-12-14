using Test
using AstroModels

@testset "CADModel default constructor" begin
    model = CADModel()
    
    @test model.file_path == ""
    @test model.scale == 1.0
    @test model.visible == false
end

@testset "CADModel keyword constructor with valid inputs" begin
    # With file path and default scale/visible
    model1 = CADModel(file_path="assets/model.obj")
    @test model1.file_path == "assets/model.obj"
    @test model1.scale == 1.0
    @test model1.visible == false
    
    # With all parameters
    model2 = CADModel(file_path="assets/DeepSpace1.obj", scale=10.0, visible=true)
    @test model2.file_path == "assets/DeepSpace1.obj"
    @test model2.scale == 10.0
    @test model2.visible == true
    
    # With different scale values
    model3 = CADModel(file_path="model.stl", scale=0.5)
    @test model3.scale == 0.5
    
    model4 = CADModel(file_path="model.ply", scale=100.0)
    @test model4.scale == 100.0
    
    # Visible false with file path (valid)
    model5 = CADModel(file_path="hidden.obj", visible=false)
    @test model5.file_path == "hidden.obj"
    @test model5.visible == false
end

@testset "CADModel scale validation" begin
    # Zero scale should throw
    @test_throws ArgumentError CADModel(scale=0.0)
    
    # Negative scale should throw
    @test_throws ArgumentError CADModel(scale=-1.0)
    @test_throws ArgumentError CADModel(scale=-10.5)
    
    # Verify error message
    err = try
        CADModel(scale=-5.0)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("scale must be positive", err.msg)
    @test occursin("-5.0", err.msg)
end

@testset "CADModel visible without file_path validation" begin
    # visible=true with empty file_path should throw
    @test_throws ArgumentError CADModel(visible=true)
    @test_throws ArgumentError CADModel(file_path="", visible=true)
    
    # Verify error message
    err = try
        CADModel(visible=true)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("visible=true requires a non-empty file_path", err.msg)
end

@testset "CADModel Base.show for empty model" begin
    model = CADModel()
    
    output = repr(MIME"text/plain"(), model)
    @test output == "CADModel: (no model)"
end

@testset "CADModel Base.show for configured model" begin
    model = CADModel(file_path="assets/DeepSpace1.obj", scale=10.0, visible=true)
    
    output = repr(MIME"text/plain"(), model)
    
    # Check that output contains expected components
    @test occursin("CADModel:", output)
    @test occursin("file_path = \"assets/DeepSpace1.obj\"", output)
    @test occursin("scale     = 10.0", output)
    @test occursin("visible   = true", output)
    
    # Test with visible=false
    model2 = CADModel(file_path="model.stl", scale=5.0, visible=false)
    output2 = repr(MIME"text/plain"(), model2)
    
    @test occursin("visible   = false", output2)
end

nothing
