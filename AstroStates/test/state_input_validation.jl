

@testset "cart_to_kep Input Validation" begin
    # Invalid length
    err = try
        cart_to_kep([1.0,2.0,3.0,4.0,5.0,6.0,7.0], μ)
        nothing
    catch e
        e
    end
    @test isa(err, ErrorException)
    @test occursin("Input vector must have exactly six elements", err.msg)

    # Degenerate position vector
    @testset "Degenerate r" begin
        io = IOBuffer()
        Logging.with_logger(ConsoleLogger(io, Logging.Warn)) do
            result = cart_to_kep([0.0,0.0,0.0, 1.0,1.0,1.0], μ)
            @test all(isnan.(result))
        end
        log_output = String(take!(io))
        @test occursin("Orbit is singular due to degenerate position", log_output)
    end

    # Degenerate velocity vector
    @testset "Degenerate v" begin
        io = IOBuffer()
        Logging.with_logger(ConsoleLogger(io, Logging.Warn)) do
            result = cart_to_kep([1.0,1.0,1.0, 0.0,0.0,0.0], μ)
            @test all(isnan.(result))
        end
        log_output = String(take!(io))
        @test occursin("Conversion failed: Orbit is singular due to degenerate position or velocity vector", log_output)
    end

    # r ∥ v (zero angular momentum)
    @testset "Parallel r and v" begin
        io = IOBuffer()
        Logging.with_logger(ConsoleLogger(io, Logging.Warn)) do
            result = cart_to_kep([1.0,1.0,1.0, 2.0,2.0,2.0], μ)
            @test all(isnan.(result))
        end
        log_output = String(take!(io))
        @test occursin("Conversion Failed: Orbit is singular due to degenerate angular momentum", log_output)
    end

    # Invalid μ
    @testset "Zero mu" begin
        io = IOBuffer()
        Logging.with_logger(ConsoleLogger(io, Logging.Warn)) do
            result = cart_to_kep([7000.0,0.0,0.0, 0.0,7.5,0.0], 0.0)
            @test all(isnan.(result))
        end
        log_output = String(take!(io))
        @test occursin("Conversion Failed: μ < tolerance", log_output)
    end
end
