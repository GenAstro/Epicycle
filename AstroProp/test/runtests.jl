
using ForwardDiff
using Test
using OrdinaryDiffEq
using LinearAlgebra

using AstroBase
using AstroStates
using AstroEpochs
using AstroUniverse
using AstroCoords
using AstroModels
using AstroMan
using AstroFun
using AstroProp

earth.mu = 398600.4415
sun.mu =132712440017.99
moon.mu = 4902.8005821478

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

function isapproxvec_percent(x, y; tol=1e-12)
    all(percent_error(xi, yi) ≤ tol for (xi, yi) in zip(x, y));
end

gravity_wrapped(posvel::Vector{T}) where T = begin
    acc = zeros(T, 6)
    t = Time("2015-09-21T12:23:12", TDB(), ISOT())
    compute_point_mass_gravity!(t, posvel, acc, earth, (moon,sun))
    return acc
end

t = Time("2015-09-21T12:23:12", TDB(), ISOT())
posvel = [200000.0, 20000.000, -400000, 5.0, -6.0, -7.0]
acc = zeros(6)

@testset "Point Mass Grav - Total accel (earth)" begin
    acc = zeros(6)
    compute_point_mass_gravity!(t,posvel,acc,earth,(moon,sun))

    # Extract acceleration (skip velocity part): truth date is from GMAT
    acc_expected = [posvel[4], posvel[5], posvel[6], -8.793174990480165e-07 , -7.275450071077514e-08 , 1.812613891537005e-06 ]

    @test isapprox(acc, acc_expected; rtol=1e-12)

end

@testset "Point Mass Grav - Pert accel" begin
    acc = zeros(6)
    compute_point_mass_gravity!(t,posvel,acc,earth,(moon,sun);include_center = false)

    # Extract acceleration (skip velocity part): truth date is from GMAT
    acc_expected = [posvel[4], posvel[5], posvel[6], 9.312960128176688e-09, 1.610854520684418e-08, 3.535297318461839e-08]

    @test isapprox(acc, acc_expected; rtol=1e-12)

end

@testset "Point Mass Grav - PosVel Jac" begin
    # Compute PosVel Jacobian
    acc = zeros(6)
    jac_posvel = Dict(PosVel => zeros(6, 6))
    compute_point_mass_gravity!(t,posvel,acc,earth,(moon,sun);jac = jac_posvel)
    J = ForwardDiff.jacobian(gravity_wrapped, posvel)
    jac_posvel[PosVel] - J

    @test isapprox(jac_posvel[PosVel], J; rtol=1e-14)
end

@testset "Point Mass Grav - error r = 0 " begin
    t = Time(2451545.0, 0.0, TT(), JD())
    posvel = zeros(6)
    acc = zeros(6)
    central = CelestialBody("Earth", 398600.4418, 6378.137, 1/298.257223563, 399)
    pert_bodies = ()

    err = try
        compute_point_mass_gravity!(t, posvel, acc, central, pert_bodies)
        nothing
    catch e
        e
    end

    @test isa(err, ErrorException)
    @test occursin("approaching singularity", err.msg)
end

@testset "Point Mass Grav - error r = moon" begin

    #  Put spacecraft at center of the moon and test error is thrown
    t_tdb = t.tdb
    pos = AstroUniverse.translate(earth,moon,t_tdb.jd)
    posvel = zeros(6)
    posvel[1:3] = pos;
    err = try
        compute_point_mass_gravity!(t, posvel, acc, earth, (moon,))
        nothing
    catch e
        e
    end

    @test isa(err, ErrorException)
    @test occursin("Perturbing body vector is less than tol", err.msg)

end

@testset "Point Mass Grav - duplicate check" begin

    err = try
        PointMassGravity(earth, (moon, earth)) 
        nothing
    catch e
        e
    end

    @test isa(err, ErrorException)
    @test occursin("The CelestialBody Earth is included in force model multiple times.", err.msg)
end

@testset "Point Mass Gravity Propagation -Earth" begin
    
    # Create spacecraft
    sat = Spacecraft(
        state=CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]), 
        time=Time("2015-09-21T12:23:12", TAI(), ISOT())
        )

    # Create force models, integrator, and dynamics system
    pm_grav = PointMassGravity(earth,(moon,sun))
    forces = ForceModel(pm_grav)
    integ = IntegratorConfig(DP8(); abstol = 1e-11, reltol = 1e-11, dt = 200)

    # Define which spacecraft to propagate and which force model to use
    dynsys = DynSys(
        forces = forces, 
        spacecraft = [sat])

    # Propagate for 10 days
    propagate(dynsys, integ, StopAtDays(10.0))
    expect_vec = [ -6767.7365586489, -733.99885811324, -1.6788146940718, 1.1437377582848, -7.6333483270125, -0.0307369480771];
    result_vec = to_posvel(sat);
    dx = abs.(expect_vec - result_vec)
    if any(dx[1:3] .> 1e-4)
        println("❌ Point Mass Propagtion (Earth-centered) Failed")
        @test false
    else
        @test true
    end

    if any(dx[4:6] .> 1e-7)
        println("❌ Point Mass Propagtion (Earth-centered) Failed")
        @test false
    else
        @test true
    end

end

t = Time("2023-09-21T12:23:12", TDB(), ISOT())
posvel = [11.540658048439733, 2331.8764015262946, 385.6799889600088, -1.8701042509951207, 0.6293064348540126, 0.21705285968096613];
acc = zeros(6)

@testset "Point Mass Grav - Total accel (moon)" begin
    acc = zeros(6)
    compute_point_mass_gravity!(t,posvel,acc,moon,(earth,sun))

    # Extract acceleration (skip velocity part): truth date is from GMAT
    acc_expected = [posvel[4], posvel[5], posvel[6], -4.270867892330212e-06 , -0.0008658160451536338, -0.0001431876516428954];

    @test isapprox(acc, acc_expected; rtol=1e-14)

end

@testset "Point Mass Gravity Propagation -Moon" begin
    
    # Create spacecraft
    sat = Spacecraft(
        state = CartesianState([11.540658048439733, 2331.8764015262946, 385.6799889600088, -1.8701042509951207, 0.6293064348540126, 0.21705285968096613]), 
        time = Time("2023-09-21T12:23:12", TDB(), ISOT())
        )

    # Create force models, integrator, and dynamics system
    pm_grav = PointMassGravity(moon,(earth,sun))
    forces = ForceModel(pm_grav)
    integ = IntegratorConfig(DP8(); abstol = 1e-11, reltol = 1e-11, dt = 200)

    # Define which spacecraft to propagate and which force model to use
    dynsys = DynSys(
        forces = forces, 
        spacecraft = [sat])

    # Propagate for 10 day
    propagate(dynsys, integ, StopAtDays(5.0))
    expect_vec = [-15192.169060110,  9323.5690156746,  6396.2477992693,  -0.5245710649352, -0.2561104305432,  0.0426186946048 ];
    result_vec = to_posvel(sat);
    dx = abs.(expect_vec - result_vec)
    if any(dx[1:3] .> 1e-4)
        println("❌ Point Mass Propagtion (Earth-centered) Failed")
        @test false
    else
        @test true
    end

    if any(dx[4:6] .> 1e-8)
        println("❌ Point Mass Propagtion (Earth-centered) Failed")
        @test false
    else
        @test true
    end

end

@testset "Point Mass Grav - PosVel Jac Accel Interface" begin

    # Compute PosVel Jacobian
    time=Time("2015-09-21T12:23:12", TAI(), ISOT())
    posvel = [7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]
    sat = Spacecraft(
        state=CartesianState(posvel), 
        time=Time("2015-09-21T12:23:12", TAI(), ISOT())
        )

    pm_grav = PointMassGravity(earth,())
    params = []
    jacobian = Dict(PosVel => zeros(6, 6))
    acc = zeros(eltype(posvel), 6)
    accel_eval!(pm_grav, time, posvel, acc, sat, params; jac = jacobian)

    J = ForwardDiff.jacobian(gravity_wrapped, posvel)

    @test isapprox(jacobian[PosVel], J; rtol=1e-12)
end

@testset "History Storage and Type Safety" begin
    # Create spacecraft with Float64 state
    sat = Spacecraft(
        state=CartesianState([7000.0, 300.0, 0.0, 0.0, 8.5, 0.03]), 
        time=Time("2015-09-21T12:23:12", TAI(), ISOT())
    )

    # Initial history should be empty
    @test isempty(sat.history)

    # Create force models and propagate for a short time
    gravity = PointMassGravity(earth,(moon,sun))
    forces  = ForceModel(gravity)
    integ   = IntegratorConfig(Tsit5(); dt=10.0, reltol=1e-9, abstol=1e-9)
    prop    = OrbitPropagator(forces, integ)

    # Propagate to apoapsis
    propagate(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=-1))

    # Verify history was populated
    @test !isempty(sat.history)
    @test length(sat.history) == 1  # Should have one segment

    # Check the first segment
    segment = sat.history[1]
    @test !isempty(segment)

    # Verify all entries have correct types
    for (time_entry, state_entry) in segment
        @test time_entry isa AstroEpochs.Time{Float64}
        @test state_entry isa Vector{Float64}
        @test length(state_entry) == 6  # Position and velocity
    end

    # Verify chronological ordering
    times = [entry[1].jd for entry in segment]
    @test issorted(times)

    # Test that we can access the history without type issues
    first_time = segment[1][1]
    last_time = segment[end][1]
    @test last_time.jd > first_time.jd

    # Verify the state values are reasonable (not NaN or Inf)
    for (_, state) in segment
        @test all(isfinite.(state))
    end
end

include("runtests_dummyforce.jl")
#include("runtests_DynamicSystem.jl")
include("runtests_orbitcalcstop.jl")
include("runtests_legacy_stopconditions.jl")