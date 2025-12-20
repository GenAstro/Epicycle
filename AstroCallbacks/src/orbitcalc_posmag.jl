# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

"""
    PosMag <: AbstractOrbitVar

Tag struct indicating Cartesian position vector magnitude (km) of a spacecraft.

Examples
```julia
posmag_calc = OrbitCalc(Spacecraft(), PosMag())
r = get_calc(posmag_calc)          # e.g., 7020.31
set_calc!(posmag_calc, 10000.0)    # set |r| to 10000 km (keeps direction)
```

See also
- PositionVector, VelocityVector, SMA
- `subtypes(AbstractOrbitVar)` for a full list of supported variables
"""
struct PosMag <: AbstractOrbitVar end
calc_input_statetag(::PosMag) = Cartesian()
calc_is_settable(::PosMag) = true       # COV_EXCL_LINE (inlined)
calc_numvars(::PosMag) = 1              # COV_EXCL_LINE (inlined)
_evaluate(::PosMag, s::CartesianState) = norm(s.position)

function _set!(::PosMag, s::CartesianState, newpos::Vector{<:Real}) 
    length(newpos) == 1 || error("PosMag requires 1 element.")
    s.position = s.position / norm(s.position) * newpos[1]
    return s
end