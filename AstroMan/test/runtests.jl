# Copyright 2025 Gen Astro LLC. All Rights Reserved.
#
# This software is licensed under the GNU AGPL v3.0,
# WITHOUT ANY WARRANTY, including implied warranties of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# This file may also be used under a commercial license,
# if one has been purchased from Gen Astro LLC.
#
# By modifying this software, you agree to the terms of the
# Gen Astro LLC Contributor License Agreement.

using Test
using LinearAlgebra
using ForwardDiff
using Zygote


using AstroBase
using AstroStates
using AstroEpochs
using AstroUniverse
using AstroCoords
using AstroModels
using AstroMan



# Create spacecraft
sat1 = Spacecraft(
            state=CartesianState([7000.0, 300.0, 0.0, 0.1, 7.5, 1.0]), 
            )

@testset "Delta V: VNB Vel., Mass" begin

    # Apply VNB maneuver
    deltav1 = ImpulsiveManeuver(
        axes = VNB(),
        element1 = 0.1,
        element2 = 0.2,
        element3 = 0.3
    )
    maneuver(sat1, deltav1)

    # Extract updated velocity
    state = to_posvel(sat1)
    vel = state[4:6];

    # GMAT reference values
    ref_vel = [0.4024239562931, 7.5689619624405, 1.2092462625626]
    ref_mass = 840.8269951084173
    @test isapprox(vel, ref_vel; rtol=1e-10, atol=1e-12)
    @test isapprox(sat1.mass, ref_mass; rtol=1e-10, atol=1e-12)

    # TODO. Test history update
end

@testset "Delta V Inertial Vel., Mass" begin
    # Apply inertial maneuver
    deltav2 = ImpulsiveManeuver(
        axes = Inertial(),
        element1 = 0.4,
        element2 = 0.5,
        element3 = 0.6
    )
    maneuver(sat1, deltav2)

    # Extract updated velocity
    state = to_posvel(sat1)
    vel = state[4:6]

    # GMAT reference values
    ref_vel = [0.8024239562931, 8.0689619624405, 1.8092462625626]
    ref_mass = 559.9227073368959
    @test isapprox(vel, ref_vel; rtol=1e-10, atol=1e-12)
    @test isapprox(sat1.mass, ref_mass; rtol=1e-10, atol=1e-12)

    # TODO. Test history update
end

# Note this test tests sat1.history from first to tests!
@testset "Spacecraft history update after maneuvers" begin
    # Expect two segments in history after the two maneuvers above
    @test length(sat1.history) == 2
    @test length(sat1.history[1]) == 1
    @test length(sat1.history[2]) == 1

    # Common initial position/velocity and applied Δv's
    pos0 = [7000.0, 300.0, 0.0]
    vel0 = [0.1, 7.5, 1.0]
    dv1_vnb = [0.1, 0.2, 0.3]      # first maneuver in VNB
    dv2_inertial = [0.4, 0.5, 0.6] # second maneuver in inertial

    # Expected velocities
    R1 = AstroMan.rot_mat_vnb_to_inertial(vcat(pos0, vel0))
    vel1 = vel0 .+ R1 * dv1_vnb
    vel2 = vel1 .+ dv2_inertial

    # First history entry
    t1, pv1 = sat1.history[1][1]
    @test t1 === sat1.time
    @test pv1[1:3] == pos0
    @test isapprox(pv1[4:6], vel1; rtol=1e-10, atol=1e-12)

    # Second history entry
    t2, pv2 = sat1.history[2][1]
    @test t2 === sat1.time
    @test pv2[1:3] == pos0
    @test isapprox(pv2[4:6], vel2; rtol=1e-10, atol=1e-12)
end

@testset "ImpulsiveManeuver Show" begin
    # Cover Base.show printing of all fields
    m = ImpulsiveManeuver(axes=Inertial(), g0=9.81, Isp=220.0,
                          element1=0.1, element2=0.2, element3=0.3)
    io = IOBuffer()
    show(io, m)
    out = String(take!(io))
    @test occursin("ImpulsiveManeuver(", out)
    @test occursin("axes = Inertial()", out)
    @test occursin("g0 = ", out)
    @test occursin("Isp = ", out)
    @test occursin("element1 = ", out)
    @test occursin("element2 = ", out)
    @test occursin("element3 = ", out)
end

@testset "ImpulsiveManeuver Input Validation Tests" begin

    m = ImpulsiveManeuver(axes=Inertial(), g0=9.81, Isp=220.0,
                          element1=0.1, element2=0.2, element3=0.3)
    # Cover get_deltav_elements one-liner
    @test AstroMan.get_deltav_elements(m) == (0.1, 0.2, 0.3)

    # Cover VNB DCM error on zero-norm r or v
    pv_zero_r = [0.0, 0.0, 0.0, 0.0, 7.5, 0.0]
    @test_throws ErrorException AstroMan.rot_mat_vnb_to_inertial(pv_zero_r)

    import AstroCoords: AbstractAxes
    struct BadAxes <: AbstractAxes end
  
    # 1) Axes must be one of MANEUVER_AXES (positional constructor -> inner validation)
    @test_throws ArgumentError ImpulsiveManeuver(BadAxes(), 9.81, 220.0, 0.0, 0.0, 0.0)

    # 2) g0 must be > 0 (keyword constructor)
    @test_throws ArgumentError ImpulsiveManeuver(axes=VNB(), g0=0.0, Isp=220.0)
    @test_throws ArgumentError ImpulsiveManeuver(axes=VNB(), g0=-9.81, Isp=220.0)

    # 3) Isp must be > 0 (keyword constructor)
    @test_throws ArgumentError ImpulsiveManeuver(axes=VNB(), g0=9.81, Isp=0.0)
    @test_throws ArgumentError ImpulsiveManeuver(axes=VNB(), g0=9.81, Isp=-1.0)

    # 4) Valid constructions succeed for both allowed axes
    m1 = ImpulsiveManeuver(axes=VNB(), g0=9.81, Isp=220.0, element1=0.0, element2=0.0, element3=0.0)
    @test m1 isa ImpulsiveManeuver

    m2 = ImpulsiveManeuver(axes=Inertial(), g0=9.81, Isp=300.0, element1=0.1, element2=0.0, element3=0.0)
    @test m2 isa ImpulsiveManeuver
    @test AstroMan.get_deltav_elements(m2) == (0.1, 0.0, 0.0)
end

@testset "ForwardDiff promotion behavior" begin
    # Common maneuver in inertial frame
    m = ImpulsiveManeuver(axes=Inertial(), g0=9.81, Isp=300.0,
                          element1=0.01, element2=0.0, element3=0.0)

    # Helper to create a Dual vector from Float64 posvel
    to_dual_vec(v::AbstractVector{<:Real}) = map(x -> ForwardDiff.Dual{Nothing}(float(x), 0.0), v)

    # Promote BOTH state and mass to Duals for AD-safe workflow
    sc = Spacecraft(
        state=CartesianState(to_dual_vec([7000.0, 300.0, 0.0, 0.1, 7.5, 1.0])),
        mass=ForwardDiff.Dual{Nothing}(1000.0, 1.0),
    )

    # Apply maneuver: state and mass remain Dual-typed
    maneuver(sc, m)
    @test eltype(to_posvel(sc)) <: ForwardDiff.Dual
    @test sc.mass isa ForwardDiff.Dual

    # compute_mass_used promotes to Dual when initial_mass is Dual
    used = AstroMan.compute_mass_used(m, ForwardDiff.Dual{Nothing}(1000.0, 1.0), m.Isp)
    @test used isa ForwardDiff.Dual
end

@testset "AD Jacobians: Δv_VNB → posvel and mass" begin
    # Spacecraft pos/vel (km, km/s) for VNB frame construction
    pv0 = [7000.0, 300.0, 0.0, 0.0, 7.5, 0.05]    # non-degenerate r, v
    R = AstroMan.rot_mat_vnb_to_inertial(pv0)     # 3x3 DCM VNB→inertial

    # Function mapping Δv in VNB (km/s) to new posvel (length-6)
    f_posvel(dv::AbstractVector) = vcat(pv0[1:3], pv0[4:6] .+ R * dv)

    # Expectation: d(pos)/d(dv)=0, d(vel)/d(dv)=R
    dv0 = [0.01, 0.02, 0.03]
    J_fd_posvel = ForwardDiff.jacobian(f_posvel, dv0)
    J_zyg_posvel = first(Zygote.jacobian(f_posvel, dv0))

    J_expected = [zeros(3,3); R]
    @test size(J_fd_posvel) == (6,3)
    @test size(J_zyg_posvel) == (6,3)
    @test isapprox(J_fd_posvel, J_expected; rtol=1e-12, atol=1e-12)
    @test isapprox(J_zyg_posvel, J_expected; rtol=1e-12, atol=1e-12)

    # Mass derivative wrt Δv in VNB using compute_mass_used
    g0 = 9.81; Isp = 300.0; m0 = 1000.0
    g_mass(dv::AbstractVector) = begin
        m = ImpulsiveManeuver(axes=VNB(), g0=g0, Isp=Isp,
                              element1=dv[1], element2=dv[2], element3=dv[3])
        m0 - AstroMan.compute_mass_used(m, m0, Isp)
    end

    # ForwardDiff/Zygote gradients agree and match closed form
    grad_fd = ForwardDiff.gradient(g_mass, dv0)
    grad_zyg = first(Zygote.gradient(g_mass, dv0))

    denom = (g0/1000) * Isp
    m_final = g_mass(dv0)
    nrm = norm(dv0)
    grad_expected = -(m_final/denom) * (dv0 / nrm)   # ∂m_final/∂dv (km/s)

    @test length(grad_fd) == 3
    @test length(grad_zyg) == 3
    @test isapprox(grad_fd, grad_expected; rtol=1e-10, atol=1e-12)
    @test isapprox(grad_zyg, grad_expected; rtol=1e-10, atol=1e-12)
end

nothing
