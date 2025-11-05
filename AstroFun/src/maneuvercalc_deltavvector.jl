# Copyright (C) 2025 Gen Astro LLC
SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

"""
    DeltaVVector <: AbstractManeuverVar

Tag struct indicating full delta-V vector [Δvx, Δvy, Δvz] of a maneuver.


Notes
- set_calc!(::ManeuverCalc, …) mutates the maneuver parameters and returns nothing.
- Applying the maneuver to the spacecraft is performed by events/sequence execution, not by setters.

Examples
```julia
# dv::ImpulsiveManeuver — replace with your maneuver
dvvec_calc = ManeuverCalc(dv, DeltaVVector())
Δv = get_calc(dvvec_calc)                 # e.g., [0.1, 0.2, 0.3]
set_calc!(dvvec_calc, [0.2, 0.3, 0.4])    # set the delta-V vector
```
"""

struct DeltaVVector <: AbstractManeuverVar end

calc_is_settable(::DeltaVVector) = true    # COV_EXCL_LINE (inlined)
calc_numvars(::DeltaVVector) = 3           # COV_EXCL_LINE (inlined)
@inline _evaluate(::DeltaVVector, man, _sc) = [man.element1, man.element2, man.element3]

function _set!(c::ManeuverCalc{M,S,DeltaVVector}, newdv::AbstractVector{<:Real}) where {M,S}
    length(newdv) == 3 || error("DeltaVVector set! expects 3 elements (got $(length(newdv))).")
    man = c.man
    setfield!(man, :element1, oftype(man.element1, newdv[1]))
    setfield!(man, :element2, oftype(man.element2, newdv[2]))
    setfield!(man, :element3, oftype(man.element3, newdv[3]))
    return man
end
