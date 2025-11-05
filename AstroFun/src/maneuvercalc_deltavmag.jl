# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

"""
    DeltaVMag <: AbstractManeuverVar

Tag struct indicating magnitude of the delta-V vector of a maneuver.

Examples
```julia
# dv::ImpulsiveManeuver — replace with your maneuver
dvmag_calc = ManeuverCalc(dv, DeltaVMag())
Δv = get_calc(dvmag_calc)     
```

See also
- DeltaVVector
- `subtypes(AbstractManeuverVar)` for a full list of supported variables
"""
struct DeltaVMag <: AbstractManeuverVar end
calc_numvars(::DeltaVMag) = 1          # COV_EXCL_LINE (inlined)
calc_is_settable(::DeltaVMag) = false  # COV_EXCL_LINE (inlined)
@inline _evaluate(::DeltaVMag, man, _sc) = [sqrt(man.element1^2 + man.element2^2 + man.element3^2)]