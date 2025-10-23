# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier:

"""
    IncomingAsymptote <: AbstractOrbitVar

Tag struct indicating full incoming hyperbolic asymptote state vector (length 6).

Examples
```julia
# sc::Spacecraft — replace with your Spacecraft instance
ia_calc = OrbitCalc(Spacecraft(), IncomingAsymptote())
x = get_calc(ia_calc)             # [rp, c3, rla, dla, bpa, ta]
set_calc!(ia_calc, [7000.0, 0.01, 1.0, 0.5, 0.2, 0.0])
```

See also
- OutGoingRLA, TA, SMA
- `subtypes(AbstractOrbitVar)` for a full list of supported variables
"""
# struct IncomingAsymptote <: AbstractOrbitVar end
calc_input_statetag(::IncomingAsymptote) = IncomingAsymptote()
calc_numvars(::IncomingAsymptote) = 6             # COV_EXCL_LINE (inlined)
calc_is_settable(::IncomingAsymptote) = true      # COV_EXCL_LINE (inlined)
_evaluate(::IncomingAsymptote, state::IncomingAsymptoteState) =
    to_vector(state)

function _set!(::IncomingAsymptote, s::IncomingAsymptoteState, newval::Vector{<:Real}) 
    length(newval) == 6 || error("IncomingAsymptote requires 6 elements.")
    # Build a new IncomingAsymptoteState with updated values
    @inbounds begin
        s = IncomingAsymptoteState(newval[1], newval[2], newval[3], newval[4], newval[5], newval[6])
    end
    return s
end

"""
    Keplerian <: AbstractOrbitVar

Tag struct indicating the full Keplerian orbital state vector (length 6).

Examples
```julia
sc = Spacecraft() 
kep_calc = OrbitCalc(sc, Keplerian())
x = get_calc(kep_calc)  
set_calc!(kep_calc, [7000.0, 0.01, 0.1, 1.0, 0.5, 0.0])

See also

SMA, ECC, INC, RAAN, AOP, TA
subtypes(AbstractOrbitVar) for a full list of supported variables
"""

calc_input_statetag(::Keplerian) = Keplerian()  
calc_numvars(::Keplerian) = 6                   # COV_EXCL_LINE (inlined)
calc_is_settable(::Keplerian) = true            # COV_EXCL_LINE (inlined)

_evaluate(::Keplerian, state::KeplerianState) = to_vector(state)

function _set!(::Keplerian, s::KeplerianState, newval::Vector{<:Real})
    length(newval) == 6 || error("Keplerian requires 6 elements.")
    @inbounds begin s = KeplerianState(newval[1], newval[2], newval[3], 
    newval[4], newval[5], newval[6]) end 
    return s 
end 