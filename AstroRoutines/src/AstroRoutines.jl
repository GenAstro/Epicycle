# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

"""
Module containing low-level routines for classical astrodynamics: anomaly conversions for elliptic
and hyperbolic orbits, two-body quantities of a single state, and the circular restricted
three-body problem.
"""
module AstroRoutines

export mean_to_eccentric_anomaly, eccentric_to_mean_anomaly
export eccentric_to_true_anomaly, true_to_eccentric_anomaly
export mean_to_true_anomaly, true_to_mean_anomaly
export true_to_hyperbolic_anomaly, hyperbolic_to_true_anomaly
export hyperbolic_to_mean_anomaly, mean_to_hyperbolic_anomaly

# Two-body quantities of a single state. See two_body.jl.
export specific_energy, angular_momentum, angular_momentum_magnitude, orbit_normal
export orbital_angular_velocity, eccentricity_vector, orbital_period, mean_motion
export hyperbolic_excess_speed, apoapsis_speed, periapsis_speed

# Circular restricted three-body problem, in the rotating frame and normalized
# units. See cr3bp.jl.
export cr3bp_mass_ratio, cr3bp_accel, cr3bp_eom!, cr3bp_jacobian
export jacobi_constant, libration_point
export cr3bp_stm_eom!, cr3bp_stm_initial

# Internal: validate eccentricity is in the elliptic domain [0, 1). Fails loudly on bad input
# (FR-ANOM-7; coding standards §9.2, §9.4) rather than returning a plausible wrong value.
function _check_elliptic_eccentricity(e::Real)
    if e < 0 || e >= 1
        throw(ArgumentError("eccentricity e must be in [0, 1) for elliptic anomaly conversions; got e = $(e)"))
    end
    return nothing
end

# Internal: validate eccentricity is in the hyperbolic domain (1, ∞) (FR-ANOM-12).
function _check_hyperbolic_eccentricity(e::Real)
    if !(e > 1)
        throw(ArgumentError("eccentricity e must be greater than 1 for hyperbolic anomaly conversions; got e = $(e)"))
    end
    return nothing
end

include("mean_to_eccentric_anomaly.jl")
include("eccentric_to_mean_anomaly.jl")
include("eccentric_to_true_anomaly.jl")
include("true_to_eccentric_anomaly.jl")
include("mean_to_true_anomaly.jl")
include("true_to_mean_anomaly.jl")
include("true_to_hyperbolic_anomaly.jl")
include("hyperbolic_to_true_anomaly.jl")
include("hyperbolic_to_mean_anomaly.jl")
include("mean_to_hyperbolic_anomaly.jl")
include("two_body.jl")
include("cr3bp.jl")

end # module AstroRoutines
