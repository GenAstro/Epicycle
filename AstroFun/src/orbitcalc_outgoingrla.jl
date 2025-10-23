# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier:

"""
    OutGoingRLA <: AbstractOrbitVar

Tag struct indicating outgoing hyperbolic asymptote right ascension (radians).

Examples
```julia
# sc::Spacecraft — replace with your Spacecraft instance
rla_calc = OrbitCalc(Spacecraft(), OutGoingRLA())
Ω_out = get_calc(rla_calc)       # e.g., 1.047
set_calc!(rla_calc, pi/3)        # set outgoing RLA to 60 degrees
```

See also
- SMA, TA, RAAN
- `subtypes(AbstractOrbitVar)` for a full list of supported variables
"""
struct OutGoingRLA <: AbstractOrbitVar end
calc_input_statetag(::OutGoingRLA) = OutGoingAsymptote()
calc_numvars(::OutGoingRLA) = 1           # COV_EXCL_LINE (inlined)
calc_is_settable(::OutGoingRLA) = true    # COV_EXCL_LINE (inlined)
_evaluate(::OutGoingRLA, state::OutGoingAsymptoteState) = state.rla

function _set!(::OutGoingRLA, s::OutGoingAsymptoteState, newval::Vector{<:Real}) 
    length(newval) == 1 || error("OutGoingRLA requires 1 element.")
    # Build a new OutGoingAsymptoteState with updated RLA
    @inbounds begin
        s = OutGoingAsymptoteState(s.rp, s.c3, newval[1], s.dla, s.bpa, s.ta)
    end
    return s
end