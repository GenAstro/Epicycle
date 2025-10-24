# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier:

"""
    PosDotVel <: AbstractOrbitVar

Tag struct indicating dot product of Cartesian position and velocity vectors.

Examples
```julia
sc = Spacecraft()
posdotvel_calc = OrbitCalc(Spacecraft(), PosDotVel())
r = get_calc(posdotvel_calc)
set_calc!(posdotvel_calc, [7000.0, 300.0, 0.0])
```
"""
struct PosDotVel <: AbstractOrbitVar end
calc_input_statetag(::PosDotVel) = Cartesian()
calc_is_settable(::PosDotVel) = false   # COV_EXCL_LINE (inlined)
calc_numvars(::PosDotVel) = 1           # COV_EXCL_LINE (inlined)
_evaluate(::PosDotVel, s::CartesianState) = dot(s.posvel[1:3], s.posvel[4:6])


