include("dummy_force.jl")

using Test
using ForwardDiff
#using StaticArrays
using AstroBase
using AstroEpochs
using AstroProp

# Dummy spacecraft (assumes Spacecraft has a default constructor)
dummy_spacecraft = Spacecraft()

# Parameters using empty struct types as keys
dummy_params = Dict(
    DummyParam() => 3e-8,
    DummyThrust() => -1e-8 .* ones(3)
)

# AD wrapper
function the_force_wrapped(x)
    accel_eval(TheForceIsWithYou(), Time("2015-09-21T12:23:12", TDB(), ISOT()), x, dummy_spacecraft, dummy_params)
end

# Test state and time
t = Time("2015-09-21T12:23:12", TDB(), ISOT())
posvel = [200000.0, 20000.0, -400000.0, 5.0, -6.0, -7.0]

# Compute some quantitues needed to test output
r̄ = posvel[1:3]
v̄ = posvel[4:6]
ĥ = cross(r̄, v̄) / norm(cross(r̄, v̄))
force = dummy_params[DummyParam()] * ĥ + dummy_params[DummyThrust()]

@testset "TheForceIsWithYou - Accel" begin
    acc = zeros(6)
    compute_the_force!(t, posvel, acc, dummy_params)
    sol = zeros(6, 1)
    sol[1:3] = v̄
    sol[4:6] = force
    @test isapprox(acc, sol; rtol=1e-12)
end
    
# Define the wrapper for ForwardDiff
compute_the_force_wrapped(x::Vector{T}) where T = begin
    acc = zeros(T, 6)
    t = Time("2015-09-21T12:23:12", TDB(), ISOT())
    dummy_params = Dict(
    DummyParam() => 3e-8,
    DummyThrust() => -1e-8 .* ones(3)
    )
    compute_the_force!(t, x, acc, dummy_params)
    return acc
end

# Prep for jacobian tests
acc = zeros(6)
jacs = Dict(PosVel => zeros(6, 6), DummyThrust => zeros(6, 3), DummyParam => zeros(6, 1))
compute_the_force!(t, posvel, acc, dummy_params; jac = jacs)
J_posvel = ForwardDiff.jacobian(compute_the_force_wrapped, posvel)
sol = zeros(6, 3)
sol[4:6,1:3] = Matrix{Float64}(I, 3, 3)
sol_thrust = sol;
sol = zeros(6, 1)
sol[1:3] = zeros(3)
sol[4:6] = ĥ
sol_param = sol;

@testset "TheForceIsWithYou Jacobians" begin
    @test isapprox(jacs[PosVel], J_posvel; rtol=1e-12)
    @test isapprox(jacs[DummyThrust], sol_thrust; rtol=1e-12)
    @test isapprox(jacs[DummyParam], sol_param; rtol=1e-12)
end

# Compute jacs using propagator interface
sat = Spacecraft(
    state=CartesianState(posvel), 
    time= t
    )

the_force = TheForceIsWithYou();
jacs = Dict(PosVel => zeros(6, 6), DummyThrust => zeros(6, 3), DummyParam => zeros(6, 1))
acc = zeros(6)
accel_eval!(the_force, t, posvel, acc, sat, dummy_params; jac = jacs)

@testset "TheForceIsWithYou Jacobians: Prop Interface" begin
    @test isapprox(jacs[PosVel], J_posvel; rtol=1e-12)
    @test isapprox(jacs[DummyThrust], sol_thrust; rtol=1e-12)
    @test isapprox(jacs[DummyParam], sol_param; rtol=1e-12)
end
