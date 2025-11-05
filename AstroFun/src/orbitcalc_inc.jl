# Copyright (C) 2025 Gen Astro LLC
SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

"""
    Inc <: AbstractOrbitVar

Tag struct indicating Keplerian inclination (radians) of a spacecraft.

Examples
```julia
sc = Spacecraft()
inc_calc = OrbitCalc(sc, Inc())
inc = get_calc(inc_calc)                # e.g., 0.4974 (28.5 degrees)
set_calc!(inc_calc, deg2rad(2.0))     # set inclination to 2 degrees
```
"""
struct Inc <: AbstractOrbitVar end
calc_numvars(::Inc) = 1             # COV_EXCL_LINE (inlined)
calc_is_settable(::Inc) = true      # COV_EXCL_LINE (inlined)
calc_input_statetag(::Inc) = Keplerian()
_evaluate(::Inc, kep::KeplerianState) = kep.inc

function _set!(::Inc, s::KeplerianState, newval::Vector{<:Real}) 
    length(newval) == 1 || error("Inc requires 1 element.")
    # Build a new KeplerianState with updated inclination
    @inbounds begin
        s = KeplerianState(s.sma, s.ecc, newval[1], s.aop, s.raan, s.ta)
    end
    return s
end