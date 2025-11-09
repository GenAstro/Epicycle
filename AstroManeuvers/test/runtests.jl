# Copyright (C) 2025 Gen Astro LLC


using Test
using LinearAlgebra
using ForwardDiff
using Zygote


using AstroBase
using AstroStates
using AstroEpochs
using AstroUniverse
using AstroFrames
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

    # First history entry - now stored as Float64
    t1, pv1 = sat1.history[1][1]
    @test t1 isa Time{Float64}  # History always uses Float64
    @test pv1[1:3] == pos0
    @test isapprox(pv1[4:6], vel1; rtol=1e-10, atol=1e-12)

    # Second history entry - now stored as Float64
    t2, pv2 = sat1.history[2][1]
    @test t2 isa Time{Float64}  # History always uses Float64
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

    import AstroFrames: AbstractAxes
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

@testset "ImpulsiveManeuver promote function" begin
    # Create a baseline maneuver with mixed input types
    m_base = ImpulsiveManeuver(
        axes=VNB(), 
        g0=9.81f0,        # Float32
        Isp=300.0,        # Float64
        element1=0.01,    # Float64
        element2=0.02f0,  # Float32
        element3=0.03     # Float64
    )
    
    # Test 1: Original type promotion in constructor
    @test typeof(m_base.g0) == Float64      # All promoted to Float64
    @test typeof(m_base.Isp) == Float64
    @test typeof(m_base.element1) == Float64
    @test typeof(m_base.element2) == Float64
    @test typeof(m_base.element3) == Float64
    
    # Test 2: Promote to ForwardDiff.Dual
    m_dual = promote(m_base, ForwardDiff.Dual{Float64})
    @test m_dual.axes === m_base.axes  # Axes unchanged
    @test typeof(m_dual.g0) <: ForwardDiff.Dual
    @test typeof(m_dual.Isp) <: ForwardDiff.Dual
    @test typeof(m_dual.element1) <: ForwardDiff.Dual
    @test typeof(m_dual.element2) <: ForwardDiff.Dual
    @test typeof(m_dual.element3) <: ForwardDiff.Dual
    
    # Test 3: Values preserved during promotion
    @test m_dual.g0.value == m_base.g0
    @test m_dual.Isp.value == m_base.Isp
    @test m_dual.element1.value == m_base.element1
    @test m_dual.element2.value == m_base.element2
    @test m_dual.element3.value == m_base.element3
    
    # Test 4: Promote to BigFloat
    m_big = promote(m_base, BigFloat)
    @test typeof(m_big.g0) == BigFloat
    @test typeof(m_big.Isp) == BigFloat
    @test typeof(m_big.element1) == BigFloat
    @test typeof(m_big.element2) == BigFloat
    @test typeof(m_big.element3) == BigFloat
    @test m_big.g0 ≈ m_base.g0
    @test m_big.Isp ≈ m_base.Isp
    @test m_big.element1 ≈ m_base.element1
    
    # Test 5: Promote to Float32 (downcast)
    m_f32 = promote(m_base, Float32)
    @test typeof(m_f32.g0) == Float32
    @test typeof(m_f32.Isp) == Float32
    @test typeof(m_f32.element1) == Float32
    @test m_f32.g0 ≈ Float32(m_base.g0)
    @test m_f32.Isp ≈ Float32(m_base.Isp)
    
    # Test 6: Promote with Inertial axes (verify axes preservation)
    m_inertial = ImpulsiveManeuver(axes=Inertial(), g0=9.81, Isp=220.0, 
                                   element1=0.1, element2=0.2, element3=0.3)
    m_inertial_dual = promote(m_inertial, ForwardDiff.Dual{Float64})
    @test m_inertial_dual.axes isa Inertial
    @test typeof(m_inertial_dual.g0) <: ForwardDiff.Dual
    
    # Test 7: Chain promotions (should work)
    m_chain = promote(promote(m_base, BigFloat), Float64)
    @test typeof(m_chain.g0) == Float64
    @test m_chain.g0 ≈ m_base.g0
end

@testset "promote function with AD workflow" begin
    # Test realistic AD workflow using promoted maneuver
    m_orig = ImpulsiveManeuver(axes=VNB(), g0=9.81, Isp=300.0,
                               element1=0.01, element2=0.02, element3=0.03)
    
    # Function that uses a promoted maneuver for AD
    function test_ad_workflow(dv_scale::Real)
        # Scale the delta-v and promote for AD
        m_scaled = ImpulsiveManeuver(axes=VNB(), g0=9.81, Isp=300.0,
                                     element1=dv_scale * 0.01, 
                                     element2=dv_scale * 0.02, 
                                     element3=dv_scale * 0.03)
        m_dual = promote(m_scaled, typeof(dv_scale))
        
        # Return total delta-v magnitude
        return sqrt(m_dual.element1^2 + m_dual.element2^2 + m_dual.element3^2)
    end
    
    # Test gradient computation
    scale_0 = 1.0
    grad = ForwardDiff.derivative(test_ad_workflow, scale_0)
    
    # Expected: d/ds[sqrt((s*0.01)² + (s*0.02)² + (s*0.03)²)] at s=1
    # = (0.01² + 0.02² + 0.03²) / sqrt(0.01² + 0.02² + 0.03²)
    expected_grad = (0.01^2 + 0.02^2 + 0.03^2) / sqrt(0.01^2 + 0.02^2 + 0.03^2)
    @test isapprox(grad, expected_grad; rtol=1e-12, atol=1e-14)
    
    # Test that regular (non-AD) evaluation gives same result
    regular_result = test_ad_workflow(1.0)
    expected_result = sqrt(0.01^2 + 0.02^2 + 0.03^2)
    @test isapprox(regular_result, expected_result; rtol=1e-12, atol=1e-14)
end

nothing
