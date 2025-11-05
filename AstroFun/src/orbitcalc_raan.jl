# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

"""
    RAAN <: AbstractOrbitVar

Tag struct indicating Keplerian right ascension of ascending node (radians) of a spacecraft.

Examples
```julia
raan_calc = OrbitCalc(Spacecraft(), RAAN())
Ω = get_calc(raan_calc)            # e.g., 1.234
set_calc!(raan_calc, pi/2)         # set RAAN to 90 degrees
```

See also
- SMA, TA, INC
- `subtypes(AbstractOrbitVar)` for a full list of supported variables
"""
struct RAAN <: AbstractOrbitVar end
calc_input_statetag(::RAAN) = Keplerian()
calc_numvars(::RAAN) = 1           # COV_EXCL_LINE (inlined)
calc_is_settable(::RAAN) = true    # COV_EXCL_LINE (inlined)
_evaluate(::RAAN, kep::KeplerianState) = kep.raan

function _set!(::RAAN, s::KeplerianState, newval::Vector{<:Real}) 
    length(newval) == 1 || error("RAAN requires 1 element.")
    # Build a new KeplerianState with updated RAAN
    @inbounds begin
        s = KeplerianState(s.sma, s.ecc, s.inc, newval[1], s.aop, s.ta)
    end
    return s
end