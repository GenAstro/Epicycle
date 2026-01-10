
using ForwardDiff
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

# Wrapper for differentiable propagation
function propagate2_state(posvel, spacecraft, forces, integ, stopcond)
    println("entered prop_state")

    @show typeof(posvel)
    @show typeof(spacecraft.time)
    @show typeof(to_posvel(spacecraft))
    spacecraft.state = CartesianState(posvel)
    println("set state to dual")

    dynsys = DynSys(forces = forces, spacecraft = [sc])
    propagate!(dynsys, integ, stopcond)
    final_state = to_posvel(sc)
    return final_state
end

# Example usage for differentiation
posvel = [7000.0, 300.0, 0.0, 0.0, 7.5, 1.0]
sat1 = Spacecraft(state = CartesianState(posvel), time = Time("2020-09-21T12:23:12", TAI(), ISOT()))
pm_grav = PointMassGravity(earth, (moon, sun))
forces = ForceModel(pm_grav)
integ = IntegratorConfig(DP8(); abstol = 1e-11, reltol = 1e-11, dt = 4000)
stopcond = StopAtDays(1.0)

# Test the closure before using ForwardDiff
#result = propagate2_state(posvel, sat1, forces, integ, stopcond)
#println("Propagated final state:")
#println(result)

# Compute Jacobian of final state w.r.t. initial state
jac = ForwardDiff.jacobian(x -> propagate2_state(x, sat1, forces, integ, stopcond), posvel)
#println("Jacobian of final state w.r.t. initial state:")
#println(jac)