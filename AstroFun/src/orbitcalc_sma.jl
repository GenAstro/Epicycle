# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier:

"""
    SMA <: AbstractOrbitVar

Tag struct indicating Keplerian semi-major axis (km) of a spacecraft.

Examples
```julia
sma_calc = OrbitCalc(Spacecraft(), SMA())
a = get_calc(sma_calc)           # e.g., 7000.0
set_calc!(sma_calc, 10000.0)     # set SMA to 10000 km
```
See also
- TA, RAAN, PositionVector
- `subtypes(AbstractOrbitVar)` for a full list of supported variables
"""
struct SMA <: AbstractOrbitVar end
calc_numvars(::SMA) = 1             # COV_EXCL_LINE (inlined)
calc_is_settable(::SMA) = true      # COV_EXCL_LINE (inlined)
calc_input_statetag(::SMA) = Keplerian()
_evaluate(::SMA, kep::KeplerianState) = kep.sma

function _set!(::SMA, s::KeplerianState, newval::Vector{<:Real}) 
    length(newval) == 1 || error("SMA requires 1 element.")
    # Build a new KeplerianState with updated SMA
    @inbounds begin
        s = KeplerianState(newval[1], s.ecc, s.inc, s.aop, s.raan, s.ta)
    end
    return s
end