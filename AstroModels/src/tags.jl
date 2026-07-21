# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

"""
    Mass <: AbstractParamTag

Tag identifying the total mass field on a `Spacecraft` (kg).

Note: `Mass` is currently `AbstractParamTag` because spacecraft mass is treated
as a fixed parameter in the simple force models.  When the mass ODE is added
(finite-thrust / rocket equation), `Mass` will become `AbstractStateTag`.
"""
struct Mass <: AbstractParamTag end

"""Return the total mass of `sc` (kg)."""
get_field(sc::Spacecraft, ::Mass) = sc.mass

"""Set the total mass of `sc` to `v` (kg)."""
function set_field!(sc::Spacecraft, ::Mass, v::Real)
    sc.mass = convert(typeof(sc.mass), v)
    return nothing
end
