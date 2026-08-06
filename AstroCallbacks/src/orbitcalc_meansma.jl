# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

"""
    MeanSMA <: AbstractOrbitVar

Tag struct indicating Brouwer-Lyddane mean-long semi-major axis (km) of a spacecraft.

Uses `BrouwerMeanLongState` as the underlying representation; a non-Brouwer
spacecraft state is converted through `BrouwerMeanLongState`, requiring μ on the
`CoordinateSystem` origin.

Examples
```julia
sc = Spacecraft(
    state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0]),
    time  = Time("2020-09-21T12:23:12", TAI(), ISOT()),
)
mean_sma_calc = OrbitCalc(sc, MeanSMA())
a̅ = get_calc(mean_sma_calc)
set_calc!(mean_sma_calc, 10000.0)
```

See also
- SMA, Ecc, Inc
- `subtypes(AbstractOrbitVar)` for a full list of supported variables
"""
struct MeanSMA <: AbstractOrbitVar end
calc_numvars(::MeanSMA) = 1             # COV_EXCL_LINE (inlined)
calc_is_settable(::MeanSMA) = true      # COV_EXCL_LINE (inlined)
calc_input_statetag(::MeanSMA) = BrouwerMeanLong()
_evaluate(::MeanSMA, bml::BrouwerMeanLongState) = bml.sma

function _set!(::MeanSMA, s::BrouwerMeanLongState, newval::Vector{<:Real})
    length(newval) == 1 || error("MeanSMA requires 1 element.")
    # Build a new BrouwerMeanLongState with updated mean SMA
    @inbounds begin
        s = BrouwerMeanLongState(newval[1], s.ecc, s.inc, s.raan, s.aop, s.ma)
    end
    return s
end
