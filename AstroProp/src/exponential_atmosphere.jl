# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation
#
# Exponential atmosphere — a concrete AbstractDensityModel plugged into AtmosphericDrag.
# Native reimplementation adapted from SatelliteToolboxAtmosphericModels.jl (MIT).

"""
    Exponential()

Analytic exponential atmosphere: density falls off piecewise-exponentially with geodetic altitude,
using tabulated base density and scale height per altitude band.

# Notes
Base altitudes, densities, and scale heights are from Vallado, *Fundamentals of Astrodynamics and
Applications*, 4th ed., Table 8-4. Valid altitude range: ``[0, 1000]`` km (values above 1000 km use
the topmost band). Ignores space weather. Native reimplementation adapted from
`SatelliteToolboxAtmosphericModels.jl` (MIT).

# Examples
```julia
AtmosphericDrag(earth; model = Exponential())
```
"""
struct Exponential <: AbstractDensityModel end

# Exponential atmosphere (Vallado, Fundamentals of Astrodynamics 4th ed., Table 8-4).
# _EXP_H0 base altitude [km], _EXP_RHO0 nominal density [kg/m³], _EXP_H scale height [km].
const _EXP_H0 = (0.0, 25.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0,
                 130.0, 140.0, 150.0, 180.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0,
                 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0)
const _EXP_RHO0 = (1.225, 3.899e-2, 1.774e-2, 3.972e-3, 1.057e-3, 3.206e-4, 8.770e-5,
                   1.905e-5, 3.396e-6, 5.297e-7, 9.661e-8, 2.438e-8, 8.484e-9, 3.845e-9,
                   2.070e-9, 5.464e-10, 2.789e-10, 7.248e-11, 2.418e-11, 9.518e-12,
                   3.725e-12, 1.585e-12, 6.967e-13, 1.454e-13, 3.614e-14, 1.170e-14,
                   5.245e-15, 3.019e-15)
const _EXP_H  = (7.249, 6.349, 6.682, 7.554, 8.382, 7.714, 6.549, 5.799, 5.382, 5.877,
                 7.263, 9.473, 12.636, 16.149, 22.523, 29.740, 37.105, 45.546, 53.628,
                 53.298, 58.515, 60.828, 63.822, 71.835, 88.667, 124.64, 181.05, 268.00)

function _geodetic(jd, x̄, eop)
    R    = r_eci_to_ecef(DCM, J2000(), ITRF(), jd, eop)
    ecef = R * SVector{3}(x̄[1], x̄[2], x̄[3])
    return ecef_to_geodetic(ecef .* 1.0e3)          # (lat, lon, alt) [rad, rad, m]
end

function _exponential_density(h_m::Real)
    h_m ≥ 0 || throw(ArgumentError("altitude must be ≥ 0; got altitude = $h_m m"))
    h  = h_m / 1000                          # → km
    id = 28
    @inbounds for i in 1:28
        if _EXP_H0[i] - h > 0
            id = i - 1
            break
        end
    end
    id = max(id, 1)
    @inbounds return _EXP_RHO0[id] * exp(-(h - _EXP_H0[id]) / _EXP_H[id])
end

density(::Exponential, jd, x̄, eop) = _exponential_density(_geodetic(jd, x̄, eop)[3])
