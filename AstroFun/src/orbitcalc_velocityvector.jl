# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier:

"""
    VelocityVector <: AbstractOrbitVar

Tag struct indicating Cartesian velocity vector (vx, vy, vz) of a spacecraft.

Examples
```julia
# sc::Spacecraft — replace with your Spacecraft instance
velvec_calc = OrbitCalc(Spacecraft(), VelocityVector())
r = get_calc(velvec_calc)
set_calc!(velvec_calc, [7000.0, 300.0, 0.0])
```

See also
- PosMag, PositionVector
- `subtypes(AbstractOrbitVar)` for a full list of supported variables
"""
struct VelocityVector <: AbstractOrbitVar end
calc_input_statetag(::VelocityVector) = Cartesian()
calc_is_settable(::VelocityVector) = true    # COV_EXCL_LINE (inlined)
calc_numvars(::VelocityVector) = 3           # COV_EXCL_LINE (inlined)
_evaluate(::VelocityVector, s::CartesianState) = s.posvel[4:6]

function _set!(::VelocityVector, s::CartesianState, newvel::Vector{<:Real}) 
    length(newvel) == 3 || error("VelocityVector requires 3 elements.")
    @inbounds begin
        s.posvel[4:6] = newvel
    end
    return s
end