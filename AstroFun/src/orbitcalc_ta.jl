# Copyright (C) 2025 Gen Astro LLC
SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

"""
    TA <: AbstractOrbitVar

Tag struct indicating Keplerian true anomaly (radians) of a spacecraft.

Examples
```julia
ta_calc = OrbitCalc(Spacecraft(), TA())
θ = get_calc(ta_calc)            # e.g., 0.5235987756
set_calc!(ta_calc, pi/6)         # set TA to 30 degrees
```
See also
- SMA, RAAN, INC
- `subtypes(AbstractOrbitVar)` for a full list of supported variables
"""
struct TA <: AbstractOrbitVar end
calc_numvars(::TA) = 1                     # COV_EXCL_LINE (inlined)
calc_input_statetag(::TA) = Keplerian()
calc_is_settable(::TA) = true              # COV_EXCL_LINE (inlined)
_evaluate(::TA, kep::KeplerianState) = kep.ta

function _set!(::TA, s::KeplerianState, newval::Vector{<:Real}) 
    length(newval) == 1 || error("TA requires 1 element.")
    # Build a new KeplerianState with updated TA
    @inbounds begin
        s = KeplerianState(s.sma, s.ecc, s.inc, s.aop, s.raan, newval[1])
    end
    return s
end