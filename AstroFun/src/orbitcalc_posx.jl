# Copyright (C) 2025 Gen Astro LLC
SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

"""
    PosX <: AbstractOrbitVar

Tag struct indicating Cartesian position vector x-component (km) of a spacecraft.

Examples
```julia
posx_calc = OrbitCalc(Spacecraft(), PosX())
r = get_calc(posx_calc)          
set_calc!(posx_calc, 10000.0)    
```

See also
- PositionVector, VelocityVector, SMA
- `subtypes(AbstractOrbitVar)` for a full list of supported variables
"""
struct PosX <: AbstractOrbitVar end
calc_input_statetag(::PosX) = Cartesian()
calc_is_settable(::PosX) = true        # COV_EXCL_LINE (inlined)
calc_numvars(::PosX) = 1               # COV_EXCL_LINE (inlined)
_evaluate(::PosX, s::CartesianState) = s.posvel[1]

function _set!(::PosX, s::CartesianState, newpos::Vector{<:Real}) 
    length(newpos) == 1 || error("PosX requires 1 element.")
    @inbounds begin
        s.posvel[1] = newpos[1];
    end
    return s
end