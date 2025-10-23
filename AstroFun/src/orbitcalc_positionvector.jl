# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier:

"""
    PositionVector <: AbstractOrbitVar

Tag struct indicating Cartesian position vector (x, y, z) of a spacecraft.

Examples
```julia
# sc::Spacecraft — replace with your Spacecraft instance
posvec_calc = OrbitCalc(Spacecraft(), PositionVector())
r = get_calc(posvec_calc)
set_calc!(posvec_calc, [7000.0, 300.0, 0.0])
```

See also
- PosMag, VelocityVector
- `subtypes(AbstractOrbitVar)` for a full list of supported variables
"""
struct PositionVector <: AbstractOrbitVar end
calc_input_statetag(::PositionVector) = Cartesian()
calc_is_settable(::PositionVector) = true   # COV_EXCL_LINE (inlined)
calc_numvars(::PositionVector) = 3          # COV_EXCL_LINE (inlined)
_evaluate(::PositionVector, s::CartesianState) = s.posvel[1:3]

function _set!(::PositionVector, s::CartesianState, newpos::Vector{<:Real}) 
    length(newpos) == 3 || error("PositionVector requires 3 elements.")
    @inbounds begin
        s.posvel[1:3] = newpos
    end
    return s
end
