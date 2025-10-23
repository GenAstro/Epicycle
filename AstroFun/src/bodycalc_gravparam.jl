# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier:

"""
    GravParam <: AbstractBodyVar

Tag struct indicating gravitational parameter μ of a celestial body.

Examples
```julia
mu_calc = BodyCalc(earth, GravParam())
μ = get_calc(mu_calc)            
set_calc!(mu_calc, 3.986e5)       
```

See also
- BodyCalc
- `subtypes(AbstractBodyVar)` for a full list of supported variables
"""
struct GravParam <: AbstractBodyVar end
calc_is_settable(::GravParam) = true   # COV_EXCL_LINE (inlined)
calc_numvars(::GravParam) = 1          # COV_EXCL_LINE (inlined)
@inline _evaluate(::GravParam, body) = AstroUniverse.get_gravparam(body)

function set_calc!(c::BodyCalc{B,GravParam}, newmu::Real) where {B}
    body = c.body
    AstroUniverse.set_gravparam!(body, newmu)
    return body
end