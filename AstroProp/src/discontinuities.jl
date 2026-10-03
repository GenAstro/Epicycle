# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0
#
# Ending a step where a force's acceleration has a kink.
#
# The dual-cone shadow's lighting factor is continuous, but its derivative jumps where the
# occulting body's disk first touches the Sun's (the penumbra's outer edge) and where one disk first
# lies inside the other (the umbra, or an annular eclipse). The penumbra is a few seconds wide in
# LEO. A step whose stages straddle a kink commits an error its error estimate does not see: in a
# one-day LEO full-force case Vern9 was 4 to 40 cm from its converged answer at tolerances from
# 1e-9 to 1e-12, and 0.1 to 4 mm with the shadow taken out (EpicycleEnterprise/benchmark/full_force).
#
# Locating the kink after the fact does not help. A continuous callback finds the root on the
# interpolant of a step that already straddles it, and that interpolant carries the same error.
# The step has to end at the kink instead. So after each accepted step a predictor looks across the
# step the integrator proposes next, and 5 % beyond, for a zero of each kink function and adds its
# time as a tstop: that step ends at the kink and the next starts there, and no step straddles it.
# A kink beyond the look-ahead is found at a later step, and one within _KINK_MIN_GAP of a stop
# already pending is that stop, so each gets one.
#
# The predictor integrates two-body motion plus the rest of the acceleration, taken as the constant
# the last step measured, with RK4 at sixteen sub-steps, and the Sun moving linearly across the
# look-ahead. Over a LEO step that places a shadow edge to within milliseconds. A kink nearer than
# _KINK_MIN_GAP is taken as reached: a step landing a hair short of one would otherwise get a stop
# a fraction of a millisecond on, and the step size would take several steps to grow back, while
# straddling a kink by 0.05 s in a step of a minute costs about (0.05/60)² of the error it causes.
#
# A stop cuts the step it ends short, and the step-size controller would grow back from the short
# step over several steps. The dynamics on the far side of the kink are as smooth as on the near
# side, so the step after a kink is given the size the controller had proposed before the cut. The acceleration is continuous at
# a kink, so a state transition matrix needs no saltation.

# How many kink functions a force has, and their values at position `r` with the Sun at `r_sun`
# (both from the central body, km) written into `out` after index `i`; returns the last index
# written. Forces without kinks have none.
_n_kinks(::OrbitODE) = 0
_kink_values!(out, i, ::OrbitODE, r, r_sun) = i

# Each shadow model's kink functions of the Sun's apparent radius a, the occulting body's b, and
# their apparent separation c. DualCone: the penumbra's outer edge, c = a + b, and its inner edge,
# c = |b − a|, the umbra when b > a and an annular eclipse when a > b. SmoothedConical is smooth.
_shadow_kinks(::AbstractShadowModel) = ()
_shadow_kinks(::DualCone) = ((a, b, c) -> c - (a + b), (a, b, c) -> c - abs(b - a))

_n_kinks(f::SolarRadiationPressure) = length(_shadow_kinks(f.shadow))

function _kink_values!(out, i, f::SolarRadiationPressure, r, r_sun)
    kinks = _shadow_kinks(f.shadow)
    isempty(kinks) && return i
    R_ss = r .- r_sun
    a = asin(f.R_sun / norm(R_ss))
    b = asin(f.R_occ / norm(r))
    c = _angle_between(R_ss, r)
    for g in kinks
        i += 1
        out[i] = g(a, b, c)
    end
    return i
end

const _KINK_SAMPLES = 16          # predictor sub-steps over the look-ahead
const _KINK_MIN_GAP = 0.05        # s; a kink nearer than this is taken as reached

# The callback that ends a step at every kink of every force, for each spacecraft whose position
# and velocity are at the indices in `posvels`, or `nothing` when no force has a kink. `t` is
# seconds from `start_epoch`, as in the right-hand side.
function _kink_callback(forces::ForceModel, start_epoch, posvels)
    nk = sum(_n_kinks, forces.forces; init = 0)
    nk == 0 && return nothing
    center = forces.center
    center === nothing && return nothing
    μ = center.mu
    g0 = zeros(nk); g1 = zeros(nk)
    stops = Float64[]                 # the kink stops added and not yet reached
    resume_dt = Ref(0.0)              # the step proposed when the last of them was added

    # Values of every kink function at r with the Sun at r_sun.
    kinks!(out, r, r_sun) = foldl((i, f) -> _kink_values!(out, i, f, r, r_sun), forces.forces; init = 0)
    two_body(r) = -μ * r / norm(r)^3

    function predict!(integrator)
        t0 = integrator.t
        dir = sign(integrator.tdir)
        # On a kink stop: resume with the step proposed before the cut.
        i = findfirst(k -> abs(k - t0) ≤ 1e-9 * max(1.0, abs(t0)), stops)
        if i !== nothing
            deleteat!(stops, i)
            abs(resume_dt[]) > abs(get_proposed_dt(integrator)) && set_proposed_dt!(integrator, resume_dt[])
        end
        H = dir * clamp(1.05 * abs(get_proposed_dt(integrator)), 1e-3, 900.0)
        tend = last(integrator.sol.prob.tspan)
        dir * (t0 + H - tend) > 0 && (H = tend - t0)
        abs(H) < 1e-6 && return
        e0 = start_epoch + t0 / 86400.0
        s0 = SVector{3}(force_position(nothing, center, sun, e0))
        s1 = SVector{3}(force_position(nothing, center, sun, start_epoch + (t0 + H) / 86400.0))
        u, uprev = integrator.u, integrator.uprev
        hprev = t0 - integrator.tprev
        h = H / _KINK_SAMPLES
        for p in posvels
            r = SVector{3}(u[p[1]], u[p[2]], u[p[3]])
            v = SVector{3}(u[p[4]], u[p[5]], u[p[6]])
            # The acceleration beyond two-body, held at what the last step measured.
            δa = zero(r)
            if abs(hprev) > 0
                rp = SVector{3}(uprev[p[1]], uprev[p[2]], uprev[p[3]])
                vp = SVector{3}(uprev[p[4]], uprev[p[5]], uprev[p[6]])
                δa = (v - vp) / hprev - (two_body(r) + two_body(rp)) / 2
            end
            acc(x) = two_body(x) + δa
            kinks!(g0, r, s0)
            τ0 = 0.0
            for k in 1:_KINK_SAMPLES
                # One RK4 step of the predictor.
                k1r = v;              k1v = acc(r)
                k2r = v + h/2 * k1v;  k2v = acc(r + h/2 * k1r)
                k3r = v + h/2 * k2v;  k3v = acc(r + h/2 * k2r)
                k4r = v + h * k3v;    k4v = acc(r + h * k3r)
                rn = r + h/6 * (k1r + 2k2r + 2k3r + k4r)
                vn = v + h/6 * (k1v + 2k2v + 2k3v + k4v)
                τ1 = k * h
                kinks!(g1, rn, s0 + (τ1 / H) * (s1 - s0))
                for j in 1:nk
                    if sign(g0[j]) != sign(g1[j]) && g1[j] != 0
                        # Linear interpolation of the root on this sub-step: the kink function is
                        # smooth, and a sub-step is a sixteenth of the look-ahead.
                        τ = τ0 + (τ1 - τ0) * g0[j] / (g0[j] - g1[j])
                        if abs(τ) > _KINK_MIN_GAP &&
                           all(k -> abs(k - (t0 + τ)) > _KINK_MIN_GAP, stops)
                            add_tstop!(integrator, t0 + τ)
                            push!(stops, t0 + τ)
                            resume_dt[] = get_proposed_dt(integrator)
                        end
                    end
                end
                r, v, τ0 = rn, vn, τ1
                g0 .= g1
            end
        end
        u_modified!(integrator, false)
        return nothing
    end
    return DiscreteCallback((u, t, integrator) -> true, predict!;
                            initialize = (c, u, t, integrator) -> predict!(integrator),
                            save_positions = (false, false))
end
_kink_callback(::OrbitODE, start_epoch, posvels) = nothing

# Callbacks combined, either of which may be `nothing`.
_with_callback(cb, ::Nothing) = cb
_with_callback(::Nothing, extra) = extra
_with_callback(::Nothing, ::Nothing) = nothing
_with_callback(cb, extra) = CallbackSet(cb, extra)
