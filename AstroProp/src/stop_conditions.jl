# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier:

"""
    apoapsis_condition(sc::Spacecraft)

Returns a closure that evaluates the apoapsis condition `r ⋅ v = 0` 
using the spacecraft's registered `:posvel` index in the ODE registry.
"""
function apsis_condition(sc::Spacecraft)
    return function(u, t, integrator)
        idxs = integrator.p.odereg[sc][:posvel]
        pos = u[idxs[1:3]]
        vel = u[idxs[4:6]]
        return dot(pos, vel)
    end
end

function ascnode_condition(sc::Spacecraft)
    return function(u, t, integrator)
        idxs = integrator.p.odereg[sc][:posvel]
        # Condtion is z = 0
        return u[idxs[3]]
    end
end

"""
    stop_affect!(integrator)

Callback effect to stop the integration.
"""
function stop_affect!(integrator)
    terminate!(integrator)
end

"""
    StopAtApo(sc::Spacecraft)

Constructs a `ContinuousCallback` for stopping at apoapsis of `sc`.
"""
function StopAtApoapsis(sc::Spacecraft)
    cond = apsis_condition(sc)

    function affect_pos!(integrator)
        if integrator.dt > 0   # forward propagation
            return terminate!(integrator)
        end
    end

    function affect_neg!(integrator)
        if integrator.dt < 0   # backward propagation
            return terminate!(integrator)
        end
    end

    return ContinuousCallback(cond, affect_neg!, affect_pos!)
end

function StopAtAscendingNode(sc::Spacecraft)
    cond = ascnode_condition(sc)

    function affect_pos!(integrator)
        if integrator.dt < 0   # forward propagation
            return terminate!(integrator)
        end
    end

    function affect_neg!(integrator)
        if integrator.dt > 0   # backward propagation
            return terminate!(integrator)
        end
    end
    return ContinuousCallback(cond, affect_neg!, affect_pos!)
end

function StopAtPeriapsis(sc::Spacecraft)
    cond = apsis_condition(sc)

    function affect_pos!(integrator)
        if integrator.dt < 0   # forward propagation
            return terminate!(integrator)
        end
    end

    function affect_neg!(integrator)
        if integrator.dt > 0   # backward propagation
            return terminate!(integrator)
        end
    end

    return ContinuousCallback(cond, affect_neg!, affect_pos!)
end

"""
    StopAtRadius(sc::Spacecraft, target_radius::Float64)

Returns a `ContinuousCallback` that triggers when `norm(pos) - target_radius = 0`.
"""
function StopAtRadius(sc::Spacecraft, target_radius::Float64)
    condition = (u, t, integrator) -> begin
        odereg = integrator.p.odereg
        idxs = odereg[sc][:posvel]
        pos = u[idxs[1:3]]
        return norm(pos) - target_radius
    end
    return ContinuousCallback(condition, nothing, stop_affect!)
end

"""
    StopAtSeconds(seconds::Float64)

Returns a `DiscreteCallback` that stops the integrator at `t = seconds`.
"""
function StopAtSeconds(seconds::Float64)
    condition = (u, t, integrator) -> begin
        Float64(t) - seconds
    end
    return ContinuousCallback(condition, stop_affect!, stop_affect!; abstol=1e-12)
end

"""
    StopAtDays(days::Float64)

Returns a `DiscreteCallback` that stops the integrator after `days` of simulation time.
"""
StopAtDays(days::Float64) = StopAtSeconds(days * 86400.0)
