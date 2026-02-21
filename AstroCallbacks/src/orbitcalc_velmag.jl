# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

"""
    VelMagMag <: AbstractOrbitVar

Tag struct indicating Cartesian velocity vector magnitude (km/s) of a spacecraft.

Examples
```julia
velmag_calc = OrbitCalc(Spacecraft(), VelMag())
r = get_calc(velmag_calc)          # e.g., 7.5
set_calc!(velmag_calc, 10.0)       # set |v| to 10 km/s (keeps direction)
```

See also
- PositionVector, VelocityVector, SMA
- `subtypes(AbstractOrbitVar)` for a full list of supported variables
"""
struct VelMag <: AbstractOrbitVar end
calc_input_statetag(::VelMag) = Cartesian()
calc_is_settable(::VelMag) = true             # COV_EXCL_LINE (inlined)
calc_numvars(::VelMag) = 1                    # COV_EXCL_LINE (inlined)
_evaluate(::VelMag, s::CartesianState) = norm(s.velocity)

function _set!(::VelMag, s::CartesianState, new_vel::Vector{<:Real}) 
    length(new_vel) == 1 || error("VelMag requires 1 element.")
    s.velocity = s.velocity / norm(s.velocity) * new_vel[1]
    return s
end