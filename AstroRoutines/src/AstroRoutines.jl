# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

"""
Module containing low-level routines for classical astrodynamics. 


"""
module AstroRoutines

export mean_to_eccentric_anomaly, eccentric_to_mean_anomaly
export eccentric_to_true_anomaly, true_to_eccentric_anomaly
export mean_to_true_anomaly, true_to_mean_anomaly

# Circular restricted three-body problem, in the rotating frame and normalized
# units. See cr3bp.jl.
export cr3bp_mass_ratio, cr3bp_accel, cr3bp_eom!, cr3bp_jacobian
export jacobi_constant, libration_point
export cr3bp_stm_eom!, cr3bp_stm_initial

# Internal: validate eccentricity is in the elliptic domain [0, 1). Fails loudly on bad input
# (FR-ANOM-7; CodingStandards §9.2, §9.4) rather than returning a plausible wrong value.
function _check_elliptic_eccentricity(e::Real)
    if e < 0 || e >= 1
        throw(ArgumentError("eccentricity e must be in [0, 1) for elliptic anomaly conversions; got e = $(e)"))
    end
    return nothing
end

include("mean_to_eccentric_anomaly.jl")
include("eccentric_to_mean_anomaly.jl")
include("eccentric_to_true_anomaly.jl")
include("true_to_eccentric_anomaly.jl")
include("mean_to_true_anomaly.jl")
include("true_to_mean_anomaly.jl")
include("cr3bp.jl")

end # module AstroRoutines
