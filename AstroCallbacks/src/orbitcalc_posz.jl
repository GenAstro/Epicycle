# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

"""
    PosZ <: AbstractOrbitVar

Tag struct indicating Cartesian position vector z-component (km) of a spacecraft.

Examples
```julia
posz_calc = OrbitCalc(Spacecraft(), PosZ())
r = get_calc(posz_calc)          
set_calc!(posz_calc, 1000.0)    
```
"""
struct PosZ <: AbstractOrbitVar end
calc_input_statetag(::PosZ) = Cartesian()
calc_is_settable(::PosZ) = true        # COV_EXCL_LINE (inlined)
calc_numvars(::PosZ) = 1               # COV_EXCL_LINE (inlined)
_evaluate(::PosZ, s::CartesianState) = s.position[3]

function _set!(::PosZ, s::CartesianState, newpos::Vector{<:Real}) 
    length(newpos) == 1 || error("PosZ requires 1 element.")
    s.position = SVector{3}(s.position[1], s.position[2], newpos[1])
    return s
end