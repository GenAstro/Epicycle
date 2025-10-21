
function cart_kep_run_roundtrip_test(kep::Vector{Float64}, μ::Float64; label="")
    cart = kep_to_cart(kep, μ)
    kep2 = cart_to_kep(cart, μ)
    cart2 = kep_to_cart(kep2, μ)
    @testset "Round-trip test: $label" begin
        @test isapprox_orbit(cart, cart2)
    end
end



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

@testset "kep_to_cart Input Validation" begin
    # Too many elements
    err = try
        kep_to_cart([7000.0, 0.0, 0.0, 0.0, pi, pi/2, 3.0], μ)
        nothing
    catch e
        e
    end
    @test isa(err, ErrorException)
    @test occursin("Input vector must have exactly six elements: a, e, i, Ω, ω, ν.", err.msg)

    # Parabolic orbit
    io = IOBuffer()
    Logging.with_logger(ConsoleLogger(io, Logging.Warn)) do
        state = kep_to_cart([7000.0, 1.0, 0.0, 0.0, pi, pi/2], μ)
        @test all(isnan.(state))
    end
    output = String(take!(io))
    @test occursin("Conversion Failed: Orbit is parabolic or singular.", output)

    # μ = 0 case
    io = IOBuffer()
    Logging.with_logger(ConsoleLogger(io, Logging.Warn)) do
        state = kep_to_cart([7000.0, 0.1, 0.0, 0.0, pi, pi/2], 0.0)
        @test all(isnan.(state))
    end
    output = String(take!(io))
    @test occursin("Conversion Failed: μ < tolerance.", output)
end

@testset "cart_to_kep Quadrant Tests" begin
    sma = 42000.0
    ecc = 0.05
    i_tests = [pi/4, pi/2, 3pi/4, pi]
    angle_tests = [pi/4.1, 3.1*pi/4.0, 5.0*pi/4.1, 7.0*pi/4.0]

    for i in i_tests
        kep = [sma, ecc, i, 0.1, 1.1, 0.5]
        cart = kep_to_cart(kep, μ)
        kep2 = cart_to_kep(cart, μ)
        cart2 = kep_to_cart(kep2, μ)
        if !isapproxvec_percent(kep,kep2; tol = 1e-13)
            println("❌ Conversion failed: cart to kep round trip inclination quadrant test)")
            @test false
        else
            @test true
        end
        if !isapproxvec_percent(cart,cart2; tol = 1e-13)
            println("❌ Conversion failed: cart to kep round trip inclination quadrant test)")
            @test false
        else
            @test true
        end
    end

    for Ω in angle_tests, ω in angle_tests, ν in angle_tests
        kep = [sma, ecc, pi/3, Ω, ω, ν]
        cart = kep_to_cart(kep, μ)
        kep2 = cart_to_kep(cart, μ)
        cart2 = kep_to_cart(kep2, μ)
        if !isapproxvec_percent(kep,kep2; tol = 1e-13)
            println(kep)
            println(kep2)
            println("❌ Conversion failed: cart to kep round trip inclination quadrant test)")
            @test false
        else
            @test true
        end
        if !isapproxvec_percent(cart,cart2; tol = 1e-13)
            println(cart)
            println(cart2)
            println("❌ Conversion failed: cart to kep round trip inclination quadrant test)")
            @test false
        else
            @test true
        end
    end
end

@testset "cart to kep Special Case Orbits Tests" begin
    tests = Dict(
        :elliptical_equatorial => [7000.0, 0.08, 0.0, 0.5, 3*pi/2, pi/7],
        :circular_equatorial   => [10000.0, 0.0, 0.0, 2.0, 4.0, pi/7],
        :circular_inclined     => [10000.0, 0.0, 1.0, 0.1, 0.2, 2pi/7],
        :elliptic_inclined     => [10000.0, 0.1, 1.0, 2.0, 3.0, 2pi/7],
        :hyperbolic            => [-6000.0, 2.16, pi/10, pi/20, 3pi/2, 3.5pi/2],
    )

    for (label, kep) in tests
        cart = kep_to_cart(kep, μ)
        kep2 = cart_to_kep(cart, μ)
        cart2 = kep_to_cart(kep2, μ)

        if label == :elliptic_inclined || label == :hyperbolic
            if !isapproxvec_percent(kep,kep2; tol = 1e-12)
                println(kep)
                println(kep2)
                println("❌ Conversion failed: cart to kep special case round trip test)")
                @test false
            else
                @test true
            end
            if !isapproxvec_percent(cart,cart2; tol = 1e-12)
                println(cart)
                println(cart2)
                println("❌ Conversion failed: cart to kep special case round trip test)")
                @test false
            else
                @test true
            end
        end

        if label == :circular_inclined
            keptest = kep
            # lump aop into ta
            keptest[6] = kep[5] + kep[6]
            keptest[5] = 0.0;
            if !isapproxvec_percent(keptest,kep2; tol = 1e-12)
                println(keptest)
                println(kep2)
                println("❌ Conversion failed: cart to kep special case round trip test)")
                @test false
            else
                @test true
            end
            if !isapproxvec_percent(cart,cart2; tol = 1e-12)
                println(cart)
                println(cart2)
                println("❌ Conversion failed: cart to kep special case round trip test)")
                @test false
            else
                @test true
            end
        end

        if label == :circular_equatorial
            keptest = kep
            # lump aop and raan into ta
            keptest[6] = mod(kep[4] + kep[5] + kep[6],2*pi)
            keptest[4] = 0.0
            keptest[5] = 0.0;
            if !isapproxvec_percent(keptest,kep2; tol = 1e-12)
                println(keptest)
                println(kep2)
                println("❌ Conversion failed: cart to kep special case round trip test)")
                @test false
            else
                @test true
            end
            if !isapproxvec_percent(cart,cart2; tol = 1e-12)
                println(cart)
                println(cart2)
                println("❌ Conversion failed: cart to kep special case round trip test)")
                @test false
            else
                @test true
            end
        end

        if label == :elliptic_equatorial
            keptest = kep
            # lump raan into ta
            keptest[6] = mod(kep[4] + kep[6],2*pi)
            keptest[4] = 0.0
            if !isapproxvec_percent(keptest,kep2; tol = 1e-12)
                println(keptest)
                println(kep2)
                println("❌ Conversion failed: cart to kep special case round trip test)")
                @test false
            else
                @test true
            end
            if !isapproxvec_percent(cart,cart2; tol = 1e-12)
                println(cart)
                println(cart2)
                println("❌ Conversion failed: cart to kep special case round trip test)")
                @test false
            else
                @test true
            end
        end
    end
end

