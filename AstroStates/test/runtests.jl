using Test
using AstroStates
using Logging

μ = 398600.4415
tol = 1e-12

function isapproxvec_percent(x, y; tol=1e-12)
    all(percent_error(xi, yi) ≤ tol for (xi, yi) in zip(x, y));
end

# Utilities
function percent_error(x, y)
    if x isa Number && y isa Number
        denom = abs(y) < 1 ? 1.0 : abs(y)
        return abs(x - y) / denom;
    elseif x isa AbstractVector && y isa AbstractVector && eltype(x) <: Number && eltype(y) <: Number
        return maximum(percent_error.(x, y))
    else
        error("Unsupported type combination in percent_error")
    end
end

include("runtests_kep_cart.jl")

include("state_truthdata_elliptic_orbits.jl")
types = collect(keys(states))

# Pairwise State Conversion Tests
@testset "Elliptic state Conversion Permutations" begin
    for from_type in types
        for to_type in types
            try
                from_state = states[from_type]
                expected = states[to_type]

                result = to_type === :Cartesian              ? CartesianState(from_state, μ) :
                         to_type === :Keplerian              ? KeplerianState(from_state, μ) :
                         to_type === :SphericalRADEC         ? SphericalRADECState(from_state, μ) :
                         to_type === :SphericalAZIFPAState   ? SphericalAZIFPAState(from_state, μ) :
                         to_type === :ModifiedEquinoctial    ? ModifiedEquinoctialState(from_state, μ) :
                         to_type === :OutGoingAsymptoteState ? OutGoingAsymptoteState(from_state, μ) :
                         to_type === :IncomingAsymptoteState ? IncomingAsymptoteState(from_state, μ) :
                         to_type === :ModifiedKeplerianState ? ModifiedKeplerianState(from_state, μ) :
                         to_type === :EquinoctialState       ? EquinoctialState(from_state,μ) :
                         to_type === :AlternateEquinoctialState ? AlternateEquinoctialState(from_state,μ) :
                         nothing

                @test typeof(result) == typeof(expected);

                result_vec = to_vector(result)
                expect_vec = to_vector(expected)

                if !isapproxvec_percent(result_vec, expect_vec; tol=tol)
                    println("❌ Conversion failed: $from_type → $to_type")
                    @show result_vec
                    @show expect_vec
                    maxdiff = maximum(abs.(result_vec .- expect_vec))
                    println("Maximum absolute difference: ", maxdiff)
                    @test false
                end

            catch e
                println("⚠️ Exception during conversion: $from_type → $to_type")
                println(e)
                @test false
            end
        end
    end
end

truth_kep = KeplerianState(-146500.0, 1.05, deg2rad(97.0), deg2rad(345.0), deg2rad(120.0), deg2rad(300.0))
truth_modkep = ModifiedKeplerianState(7325.0, -300325.0, deg2rad(97.0), deg2rad(345.0), deg2rad(120.0), deg2rad(300.00))
truth_cart = CartesianState([4486.62555067877, -2278.090328380258, 8463.948027345208, -8.831137842558398 ,2.38254720654504, -0.1278435615976432])
truth_sphazifpa = SphericalAZIFPAState(9846.721311475412, deg2rad(-26.91924567746847), deg2rad(59.26835712578854), 9.147779553636203, deg2rad(-13.79705843709202), deg2rad(120.8067694505813))
truth_sphradec = SphericalRADECState(9846.721311475412, deg2rad(-26.91924567746847), deg2rad(59.26835712578854), 9.147779553636203, deg2rad(164.9016731348416), deg2rad(-0.8007555204065349))
truth_mee = ModifiedEquinoctialState(15016.25,  -0.2717599973576463, 1.014222117603522, 1.091780539096372, -0.2925417137628885, deg2rad(45.0))
truth_outasymptote = OutGoingAsymptoteState(7325.0, 2.720822126279856, deg2rad(14.31107485333013), deg2rad(-75.92005522066884), deg2rad(210.0630007380261),deg2rad(300.00))
truth_inasymptote  = IncomingAsymptoteState(7325.0 , 2.720822126279856, deg2rad(351.3162147414849), deg2rad(-41.86050335870885), deg2rad(189.4178404199244),deg2rad(300.0))

states = Dict(
    :Cartesian => truth_cart,
    :Keplerian => truth_kep,
    :SphericalRADEC => truth_sphradec,
    :SphericalAZIFPAState => truth_sphazifpa,
    :ModifiedEquinoctial => truth_mee,
    :OutGoingAsymptoteState => truth_outasymptote,
    :IncomingAsymptoteState => truth_inasymptote,
    :ModifiedKeplerianState => truth_modkep,
)

types = collect(keys(states))

# Pairwise State Conversion Tests
@testset "Hyperbolic State Conversion Permutations" begin
    for from_type in types
        for to_type in types
            try
                from_state = states[from_type]
                expected = states[to_type]

                result = to_type === :Cartesian              ? CartesianState(from_state, μ) :
                         to_type === :Keplerian              ? KeplerianState(from_state, μ) :
                         to_type === :SphericalRADEC         ? SphericalRADECState(from_state, μ) :
                         to_type === :SphericalAZIFPAState   ? SphericalAZIFPAState(from_state, μ) :
                         to_type === :ModifiedEquinoctial    ? ModifiedEquinoctialState(from_state, μ) :
                         to_type === :OutGoingAsymptoteState ? OutGoingAsymptoteState(from_state, μ) :
                         to_type === :IncomingAsymptoteState ? IncomingAsymptoteState(from_state, μ) :
                         to_type === :ModifiedKeplerianState ? ModifiedKeplerianState(from_state, μ) :
                         nothing

                @test typeof(result) == typeof(expected);

                result_vec = to_vector(result)
                expect_vec = to_vector(expected)

                if !isapproxvec_percent(result_vec, expect_vec; tol=tol)
                    println("❌ Conversion failed: $from_type → $to_type")
                    @show result_vec
                    @show expect_vec
                    maxdiff = maximum(abs.(result_vec .- expect_vec))
                    println("Maximum absolute difference: ", maxdiff)
                    @test false
                end

            catch e
                println("⚠️ Exception during conversion: $from_type → $to_type")
                println(e)
                @test false
            end
        end
    end
end
# Final Clean Test Summary
println("\n=========================================")
println("  AstroStates tests completed.")
println("=========================================")

@testset "OrbitState Round-Trip Tests" begin
    for (name, state) in states
        # 1. Concrete → OrbitState → Concrete
        os = OrbitState(state)
        state2 = to_state(os)
        @test typeof(state2) == typeof(state)
        @test isapproxvec_percent(to_vector(state2), to_vector(state); tol=tol)

        # 2. Concrete → OrbitState → Concrete → OrbitState → Concrete
        os2 = OrbitState(state2)
        state3 = to_state(os2)
        @test typeof(state3) == typeof(state)
        @test isapproxvec_percent(to_vector(state3), to_vector(state); tol=tol)
    end
end

include("state_truthdata_elliptic_orbits.jl")

@testset "to_vector OrbitState tests" begin
    @test !isempty(states)
    for (name, st) in states
        os = OrbitState(st)
        v_original = to_vector(st)
        v_wrapped  = to_vector(os)
        # Value equality
        @test v_wrapped == v_original
        # (Optional) Same object identity (depends on implementation; keep as info if it ever changes)
        @test v_wrapped === os.state
    end
end

@testset "OrbitState show method" begin
    for (name, state) in states
        os = OrbitState(state)
        io_os = IOBuffer()
        show(io_os, os)
        output_os = String(take!(io_os))

        io_state = IOBuffer()
        show(io_state, state)
        output_state = String(take!(io_state))

        # Split the OrbitState output into header and body
        lines = split(output_os, '\n')
        @test occursin("OrbitState:", lines[1])
        @test occursin("statetype:", lines[2])

        # The rest should match the concrete state's show output
        body_os = join(lines[3:end], "\n")
        output_state = strip(output_state)  # Remove leading/trailing whitespace
        body_os = strip(body_os)
        @test body_os == output_state
    end
end

# Dummy state type with no defined conversions
struct FakeState{T<:Real} <: AbstractOrbitState
    data::Vector{T}
end

μ = 398600.4418
fake = FakeState([1.0, 2.0, 3.0])  # contents irrelevant—fallbacks never inspect

@testset "Fallback conversion constructors" begin
    # Ensure the fallback methods exist (dispatchable on FakeState, μ::Float64)
    @test hasmethod(KeplerianState, Tuple{FakeState{Float64}, Float64})
    @test hasmethod(SphericalRADECState, Tuple{FakeState{Float64}, Float64})
    @test hasmethod(ModifiedEquinoctialState, Tuple{FakeState{Float64}, Float64})
    @test hasmethod(CartesianState, Tuple{FakeState{Float64}, Float64})
    @test hasmethod(SphericalAZIFPAState, Tuple{FakeState{Float64}, Float64})
    @test hasmethod(ModifiedKeplerianState, Tuple{FakeState{Float64}, Float64})
    @test hasmethod(EquinoctialState, Tuple{FakeState{Float64}, Float64})
    @test hasmethod(AlternateEquinoctialState, Tuple{FakeState{Float64}, Float64})
    @test hasmethod(OutGoingAsymptoteState, Tuple{FakeState{Float64}, Float64})
    @test hasmethod(IncomingAsymptoteState, Tuple{FakeState{Float64}, Float64})

    # Each should throw the defined fallback error
    @test_throws ErrorException KeplerianState(fake, μ)
    @test_throws ErrorException SphericalRADECState(fake, μ)
    @test_throws ErrorException ModifiedEquinoctialState(fake, μ)
    @test_throws ErrorException CartesianState(fake, μ)
    @test_throws ErrorException SphericalAZIFPAState(fake, μ)
    @test_throws ErrorException ModifiedKeplerianState(fake, μ)
    @test_throws ErrorException EquinoctialState(fake, μ)
    @test_throws ErrorException AlternateEquinoctialState(fake, μ)
    @test_throws ErrorException OutGoingAsymptoteState(fake, μ)
    @test_throws ErrorException IncomingAsymptoteState(fake, μ)

    # Optionally verify an error message pattern for one (spot check)
    msg = try
        KeplerianState(fake, μ)
        nothing
    catch e
        sprint(showerror, e)
    end
    @test occursin("No conversion defined", msg)
end

@testset "to_vector raw Vector fallback" begin
    v = [1.0,2,3,4,5,6]
    msg = try
        AstroStates.to_vector(v)
        nothing
    catch e
        @test e isa ErrorException
        sprint(showerror, e)
    end
    @test msg !== nothing
    @test occursin("Cannot call `to_vector` on raw Vector", msg)
end

include("runtests_types.jl")
include("runtests_identityconstructors.jl")
include("runtests_kinematictypes_nomu.jl")
include("runtests_parabolicorbits.jl")
include("runtests_cart_kep_inputValidation.jl")
include("runtests_cart_mee_inputValidation.jl")
include("runtests_cart_sphazfpa_inputValidation.jl")
include("runtests_cart_sphradec_inputValidation.jl")
include("runtests_outasymp_kep_inputValidation.jl")
include("runtests_kep_modkep_inputValidation.jl")
include("runtests_alt_equinoctial_equinoctial_inputValidation.jl")

println(" ")
println(" Starting state differentation tests... precompile takes several minutes")
println(" ")
#include("state_differentiation.jl")