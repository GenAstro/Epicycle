# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier:

"""
    Ecc <: AbstractOrbitVar

Tag struct indicating Keplerian eccentricity of an orbit.

Examples
```julia
# sc::Spacecraft — replace with your Spacecraft instance
ecc_calc = OrbitCalc(Spacecraft(), Ecc())
a = get_calc(ecc_calc)           
set_calc!(ecc_calc, 0.02)        # set Ecc to 0.02
```
See also
- TA, RAAN, SMA
- `subtypes(AbstractOrbitVar)` for a full list of supported variables
"""
struct Ecc <: AbstractOrbitVar end
calc_numvars(::Ecc) = 1                   # COV_EXCL_LINE
calc_is_settable(::Ecc) = true            # COV_EXCL_LINE
calc_input_statetag(::Ecc) = Keplerian()  
_evaluate(::Ecc, kep::KeplerianState) = kep.ecc

function _set!(::Ecc, s::KeplerianState, newval::Vector{<:Real}) 
    length(newval) == 1 || error("Ecc requires 1 element.")
    # Build a new KeplerianState with updated Ecc
    @inbounds begin
        s = KeplerianState(s.sma, newval[1], s.inc, s.aop, s.raan, s.ta)
    end
    return s
end