# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # Impulsive Maneuvers
#'
#' Apply two impulsive maneuvers to one spacecraft and compare inertial and
#' velocity-normal-binormal (VNB) component definitions. Each maneuver updates the spacecraft
#' state and reduces its mass according to its specific impulse.

using Epicycle

#' ## Configuration
#'
#' The default spacecraft is sufficient here because the focus is the maneuver, not the orbit.

# Create the spacecraft
sat = Spacecraft()

# Define an inertial-frame burn
inertial_burn = ImpulsiveManeuver(axes = Inertial(),
                                  g0 = 9.80665,
                                  Isp = 250.0,
                                  element1 = 0.04,
                                  element2 = -0.3,
                                  element3 = 0.1)

# Define a VNB-frame burn
vnb_burn = ImpulsiveManeuver(axes = VNB(),
                             g0 = 9.80665,
                             Isp = 250.0,
                             element1 = 0.2,
                             element2 = 0.1,
                             element3 = -0.2)

#' ## Apply the maneuvers
#'
#' Apply the burns in sequence. The VNB axes for the second burn are computed from the state left
#' by the first.

# Apply and report the inertial burn
println("Initial mass: ", total_mass(sat))
maneuver!(sat, inertial_burn)
println("Mass after the inertial maneuver: ", total_mass(sat))
println("State after the inertial maneuver: \n", get_state(sat, Cartesian()))

# Apply and report the VNB burn
maneuver!(sat, vnb_burn)
println("Mass after the VNB maneuver: ", total_mass(sat))
println("State after the VNB maneuver: \n", get_state(sat, Cartesian()))
