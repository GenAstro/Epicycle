# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

"""
    Mass <: AbstractParamTag

Tag identifying the total mass field on a `Spacecraft` (kg).

Note: `Mass` is currently `AbstractParamTag` because spacecraft mass is treated
as a fixed parameter in the simple force models.  When the mass ODE is added
(finite-thrust / rocket equation), `Mass` will become `AbstractStateTag`.

# Example
```jldoctest
Mass() isa Mass

# output
true
```
"""
struct Mass <: AbstractParamTag end

"""Return the total mass of `sc` (kg)."""
get_field(sc::Spacecraft, ::Mass) = total_mass(sc)

# No `set_field!` for `Mass`. Setting a total requires deciding which tank it came
# from, which the caller cannot know, so the operation is not well posed. Mass changes
# by applying a maneuver. The method that used to be here had no caller anywhere in the
# repository. When the mass model lands, the design variable a solver varies has to be
# named explicitly — dry mass or propellant load, not "mass".
