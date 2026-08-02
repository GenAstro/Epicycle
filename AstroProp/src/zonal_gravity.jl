# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation
#
# Zonal (J2–J5) geopotential — a concrete AbstractGeopotential plugged into HarmonicGravity.

"""
    Zonal()

Earth zonal gravity field — the harmonics ``J_2`` through ``J_5`` (order 0). These capture Earth's
oblateness, the largest departure from a spherical field, and cover most low-Earth-orbit analysis.

# Notes
Coefficients ``J_2 \\ldots J_5`` and the reference constants ``\\mu`` and ``R_\\oplus`` are the
EGM96 values. Earth only; the maximum supported degree is 5 and the maximum supported order is 0 —
[`HarmonicGravity`](@ref) validates your request against these via [`max_degree`](@ref) and
[`max_order`](@ref). Acceleration in the Earth-fixed frame is computed with a stable Legendre
recurrence (Vallado, *Fundamentals of Astrodynamics and Applications*, §8.7). Cross-validated
against GMAT with EGM96 at degree 5, order 0.

# Examples
```julia
HarmonicGravity(earth; degree = 5, order = 0, model = Zonal())
```
"""
struct Zonal <: AbstractGeopotential end

max_degree(::Zonal) = 5
max_order(::Zonal)  = 0

# EGM96 unnormalized zonal coefficients Jₙ = -√(2n+1)·C̄ₙ₀, and the EGM96 scale.
const _ZONAL_GM_KM3 = 398600.4415          # EGM96 GM [km³/s²]
const _ZONAL_RE_KM  = 6378.1363            # EGM96 Rₑ [km]
const _ZONAL_J = (0.0, 0.0,                # J0, J1 (unused)
                  1.0826266835531513e-3,   # J2
                 -2.5326564853322355e-6,   # J3
                 -1.6196215913670311e-6,   # J4
                 -2.2729608286869628e-7)   # J5

"Cached scale for the native zonal model (SI units)."
struct ZonalData
    mu_m::Float64        # GM [m³/s²]
    Re_m::Float64        # Rₑ [m]
end

geopotential_data(::Zonal, body::CelestialBody, degree::Int, order::Int) =
    ZonalData(_ZONAL_GM_KM3 * 1.0e9, _ZONAL_RE_KM * 1.0e3)

# Native zonal acceleration (central + J2..Jₙ) in the Earth-fixed frame, SI. Legendre
# polynomials and their derivatives use stable recurrences:
#   Pₙ  = ((2n-1) u Pₙ₋₁ - (n-1) Pₙ₋₂) / n
#   Pₙ' = u Pₙ₋₁' + n Pₙ₋₁
# with u = z/r = sin(geocentric latitude). Verified against the closed-form J2 acceleration.
function geopotential_accel(::Zonal, data::ZonalData, r_itrf, tsec,
                            degree::Int, order::Int)
    μ  = data.mu_m
    Re = data.Re_m
    x, y, z = r_itrf[1], r_itrf[2], r_itrf[3]
    r  = sqrt(x * x + y * y + z * z)
    u  = z / r
    r3 = r^3
    ax = -μ * x / r3
    ay = -μ * y / r3
    az = -μ * z / r3
    n_max = min(degree, 5)
    if n_max ≥ 2
        Pnm2, Pnm1  = 1.0, u        # P₀, P₁
        dPnm2, dPnm1 = 0.0, 1.0     # P₀', P₁'
        for n in 2:n_max
            Pn  = ((2n - 1) * u * Pnm1 - (n - 1) * Pnm2) / n
            dPn = u * dPnm1 + n * Pnm1
            Jn  = _ZONAL_J[n + 1]
            if Jn != 0.0
                c  = -μ * Jn * (Re / r)^n
                tc = -(n + 1) * Pn / r3
                ax += c * (tc * x + dPn * (-z * x / r^4))
                ay += c * (tc * y + dPn * (-z * y / r^4))
                az += c * (tc * z + dPn * (1.0 / r^2 - z * z / r^4))
            end
            Pnm2, Pnm1   = Pnm1, Pn
            dPnm2, dPnm1 = dPnm1, dPn
        end
    end
    return SVector{3}(ax, ay, az)
end
