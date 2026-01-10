using ForwardDiff
using Test

using EpicycleBase
using AstroStates
using AstroProp
using AstroEpochs
using AstroUniverse
using ForwardDiff
using AstroFrames
using Test
using OrdinaryDiffEq
using LinearAlgebra
using SciMLSensitivity
using Zygote
earth.mu = 398600.4415
sun.mu =132712440017.99
moon.mu = 4902.8005821478

#=
function test_propagate_forwarddiff()

    x0 = ForwardDiff.Dual.(collect(1:6), [7000.0, 300.0, 0.0, 0.0, 7.5, 0.03])

    # Create spacecraft
    sat = Spacecraft(
        state = CartesianState(x0), 
        time = Time("2015-09-21T12:23:12", TAI(), ISOT())
    )

    # Create force models, integrator, and dynamics system
    pm_grav = PointMassGravity(earth, ())
    forces = ForceModel(pm_grav)
    integ = IntegratorConfig(DP8(); abstol = 1e-11, reltol = 1e-11, dt = 200)

    # Define which spacecraft to propagate and which force model to use
    dynsys = DynSys(
        forces = forces, 
        spacecraft = [sat]
    )

    # Propagate for 1 hour in seconds
    result = propagate(dynsys, integ, StopAtSeconds(3600.0))
    @show result
end
=# 

function propagate_final_x(x0_vec)
    # Create spacecraft with Float64 state (Zygote works with Float64 arrays)
    sat = Spacecraft(
        state = CartesianState(x0_vec),
        time = Time("2015-09-21T12:23:12", TAI(), ISOT())
    )
    pm_grav = PointMassGravity(earth, ())
    forces = ForceModel(pm_grav)
    integ = IntegratorConfig(DP8(); abstol = 1e-11, reltol = 1e-11, dt = 200)
    dynsys = DynSys(forces = forces, spacecraft = [sat])
    result = propagate!(dynsys, integ, StopAtSeconds(3600.0))
    # Return the final x position (or any scalar you want to differentiate)
    return result.u[end][1]
end

x0 = [7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]
grad = Zygote.gradient(propagate_final_x, x0)
println("Gradient: ", grad)
#x0 = ForwardDiff.Dual.(7000.0, 300.0, 0.0, 0.0, 7.5, 0.03), ... # etc.

#@testset "Propagate ForwardDiff Test" begin
#    test_propagate_forwarddiff()
#end

#f(x) = propagate(..., state=CartesianState(x), ...)
#J = ForwardDiff.jacobian(f, x0)
