# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

"""
    IAU2015Orientation{T<:Real}

Pole orientation and prime-meridian angle for a celestial body per IAU 2015
recommendations, plus their time derivatives.

# Fields
- `ra_pole::T`             — right ascension of the pole (α₀) [rad]
- `dec_pole::T`            — declination of the pole (δ₀) [rad]
- `prime_meridian::T`      — prime meridian angle (W) [rad]
- `ra_pole_rate::T`        — α̇₀ [rad/s]
- `dec_pole_rate::T`       — δ̇₀ [rad/s]
- `prime_meridian_rate::T` — Ẇ [rad/s]

# Notes
- All angles in radians; rates in rad/s.
- The IAU 2015 formulas use TDB Julian date; callers supply that.
- Reference: Archinal et al., "Report of the IAU Working Group on Cartographic
  Coordinates and Rotational Elements: 2015" (Cel. Mech. Dyn. Astron. 130:22, 2018).

# Example
```jldoctest
julia> using AstroUniverse

julia> iau2015_orientation(mars, 2451545.0) isa IAU2015Orientation
true
```
"""
struct IAU2015Orientation{T<:Real}
    ra_pole::T
    dec_pole::T
    prime_meridian::T
    ra_pole_rate::T
    dec_pole_rate::T
    prime_meridian_rate::T
end

const _SECONDS_PER_DAY     = 86400.0
const _SECONDS_PER_CENTURY = _SECONDS_PER_DAY * 36525.0
const _J2000_TDB_JD        = 2451545.0

# =============================================================================
# Per-body implementations
#
# Each function returns an `IAU2015Orientation` computed from the IAU 2015
# recommendations at TDB Julian date `jd`. Time argument variables:
#   d = jd - _J2000_TDB_JD           # days from J2000 TDB
#   T = d / 36525                    # Julian centuries from J2000 TDB
# =============================================================================

function _iau2015_orientation_sun(jd::Real)
    d = jd - _J2000_TDB_JD
    α₀ = deg2rad(286.13)
    δ₀ = deg2rad(63.87)
    W  = deg2rad(84.176 + 14.1844000 * d)
    return IAU2015Orientation(
        α₀, δ₀, W,
        0.0, 0.0,
        deg2rad(14.1844000) / _SECONDS_PER_DAY,
    )
end

function _iau2015_orientation_mercury(jd::Real)
    d = jd - _J2000_TDB_JD
    T = d / 36525.0

    # Mean anomaly arguments (radians).
    M1 = deg2rad(174.7910857 +  4.092335  * d)
    M2 = deg2rad(349.5821714 +  8.184670  * d)
    M3 = deg2rad(164.3732571 + 12.277005  * d)
    M4 = deg2rad(339.1643429 + 16.369340  * d)
    M5 = deg2rad(153.9554286 + 20.461675  * d)

    α₀ = deg2rad(281.0103 - 0.0328 * T)
    δ₀ = deg2rad( 61.4155 - 0.0049 * T)

    W_deg = 329.5988 + 6.1385108 * d +
            0.01067257 * sin(M1) -
            0.00112309 * sin(M2) -
            0.00011040 * sin(M3) -
            0.00002539 * sin(M4) -
            0.00000571 * sin(M5)
    W = deg2rad(W_deg)

    # dW/dt in deg/day, then convert to rad/s.
    dW_deg_per_day = 6.1385108 +
            0.01067257 * cos(M1) * deg2rad( 4.092335) -
            0.00112309 * cos(M2) * deg2rad( 8.184670) -
            0.00011040 * cos(M3) * deg2rad(12.277005) -
            0.00002539 * cos(M4) * deg2rad(16.369340) -
            0.00000571 * cos(M5) * deg2rad(20.461675)

    return IAU2015Orientation(
        α₀, δ₀, W,
        deg2rad(-0.0328) / _SECONDS_PER_CENTURY,
        deg2rad(-0.0049) / _SECONDS_PER_CENTURY,
        deg2rad(dW_deg_per_day) / _SECONDS_PER_DAY,
    )
end

function _iau2015_orientation_venus(jd::Real)
    d = jd - _J2000_TDB_JD
    α₀ = deg2rad(272.76)
    δ₀ = deg2rad(67.16)
    W  = deg2rad(160.20 - 1.4813688 * d)
    return IAU2015Orientation(
        α₀, δ₀, W,
        0.0, 0.0,
        deg2rad(-1.4813688) / _SECONDS_PER_DAY,
    )
end

function _iau2015_orientation_mars(jd::Real)
    d = jd - _J2000_TDB_JD
    T = d / 36525.0

    # α₀ terms
    α_poly = 317.269202 - 0.10927547 * T
    α_pert =  0.000068 * sin(deg2rad(198.991226 + 19139.4819985 * T)) +
              0.000238 * sin(deg2rad(226.292679 + 38280.8511281 * T)) +
              0.000052 * sin(deg2rad(249.663391 + 57420.7251593 * T)) +
              0.000009 * sin(deg2rad(266.183510 + 76560.6367950 * T)) +
              0.419057 * sin(deg2rad( 79.398797 +     0.5042615 * T))
    α₀ = deg2rad(α_poly + α_pert)

    dα_pert_per_century =
              0.000068 * cos(deg2rad(198.991226 + 19139.4819985 * T)) * deg2rad(19139.4819985) +
              0.000238 * cos(deg2rad(226.292679 + 38280.8511281 * T)) * deg2rad(38280.8511281) +
              0.000052 * cos(deg2rad(249.663391 + 57420.7251593 * T)) * deg2rad(57420.7251593) +
              0.000009 * cos(deg2rad(266.183510 + 76560.6367950 * T)) * deg2rad(76560.6367950) +
              0.419057 * cos(deg2rad( 79.398797 +     0.5042615 * T)) * deg2rad(    0.5042615)
    dα₀_dt = deg2rad(-0.10927547 + dα_pert_per_century) / _SECONDS_PER_CENTURY

    # δ₀ terms
    δ_poly = 54.432516 - 0.05827105 * T
    δ_pert =  0.000051 * cos(deg2rad(122.433576 + 19139.9407476 * T)) +
              0.000141 * cos(deg2rad( 43.058401 + 38280.8753272 * T)) +
              0.000031 * cos(deg2rad( 57.663379 + 57420.7517205 * T)) +
              0.000005 * cos(deg2rad( 79.476401 + 76560.6495004 * T)) +
              1.591274 * cos(deg2rad(166.325722 +     0.5042615 * T))
    δ₀ = deg2rad(δ_poly + δ_pert)

    dδ_pert_per_century =
            -(0.000051 * sin(deg2rad(122.433576 + 19139.9407476 * T)) * deg2rad(19139.9407476)) -
             (0.000141 * sin(deg2rad( 43.058401 + 38280.8753272 * T)) * deg2rad(38280.8753272)) -
             (0.000031 * sin(deg2rad( 57.663379 + 57420.7517205 * T)) * deg2rad(57420.7517205)) -
             (0.000005 * sin(deg2rad( 79.476401 + 76560.6495004 * T)) * deg2rad(76560.6495004)) -
             (1.591274 * sin(deg2rad(166.325722 +     0.5042615 * T)) * deg2rad(    0.5042615))
    dδ₀_dt = deg2rad(-0.05827105 + dδ_pert_per_century) / _SECONDS_PER_CENTURY

    # W terms
    W_linear = 176.049863 + 350.891982443297 * d
    W_pert   = 0.000145 * sin(deg2rad(129.071773 + 19140.0328244 * T)) +
               0.000157 * sin(deg2rad( 36.352167 + 38281.0473591 * T)) +
               0.000040 * sin(deg2rad( 56.668646 + 57420.9295360 * T)) +
               0.000001 * sin(deg2rad( 67.364003 + 76560.2552215 * T)) +
               0.000001 * sin(deg2rad(104.792680 + 95700.4387578 * T)) +
               0.584542 * sin(deg2rad( 95.391654 +     0.5042615 * T))
    W = deg2rad(W_linear + W_pert)

    dW_pert_per_century =
               0.000145 * cos(deg2rad(129.071773 + 19140.0328244 * T)) * deg2rad(19140.0328244) +
               0.000157 * cos(deg2rad( 36.352167 + 38281.0473591 * T)) * deg2rad(38281.0473591) +
               0.000040 * cos(deg2rad( 56.668646 + 57420.9295360 * T)) * deg2rad(57420.9295360) +
               0.000001 * cos(deg2rad( 67.364003 + 76560.2552215 * T)) * deg2rad(76560.2552215) +
               0.000001 * cos(deg2rad(104.792680 + 95700.4387578 * T)) * deg2rad(95700.4387578) +
               0.584542 * cos(deg2rad( 95.391654 +     0.5042615 * T)) * deg2rad(    0.5042615)
    dW_dt = deg2rad(350.891982443297) / _SECONDS_PER_DAY +
            deg2rad(dW_pert_per_century) / _SECONDS_PER_CENTURY

    return IAU2015Orientation(α₀, δ₀, W, dα₀_dt, dδ₀_dt, dW_dt)
end

function _iau2015_orientation_jupiter(jd::Real)
    d = jd - _J2000_TDB_JD
    T = d / 36525.0

    Ja = deg2rad( 99.360714  + 4850.4046 * T)
    Jb = deg2rad(175.895369  + 1191.9605 * T)
    Jc = deg2rad(300.323162  +  262.5475 * T)
    Jd = deg2rad(114.012305  + 6070.2476 * T)
    Je = deg2rad( 49.511251  +   64.3000 * T)

    α_poly = 268.056595 - 0.006499 * T
    α_pert = 0.000117 * sin(Ja) +
             0.000938 * sin(Jb) +
             0.001432 * sin(Jc) +
             0.000030 * sin(Jd) +
             0.002150 * sin(Je)
    α₀ = deg2rad(α_poly + α_pert)

    dα_pert_per_century =
             0.000117 * cos(Ja) * deg2rad(4850.4046) +
             0.000938 * cos(Jb) * deg2rad(1191.9605) +
             0.001432 * cos(Jc) * deg2rad( 262.5475) +
             0.000030 * cos(Jd) * deg2rad(6070.2476) +
             0.002150 * cos(Je) * deg2rad(  64.3000)
    dα₀_dt = deg2rad(-0.006499 + dα_pert_per_century) / _SECONDS_PER_CENTURY

    δ_poly = 64.495303 + 0.002413 * T
    δ_pert = 0.000050 * cos(Ja) +
             0.000404 * cos(Jb) +
             0.000617 * cos(Jc) -
             0.000013 * cos(Jd) +
             0.000926 * cos(Je)
    δ₀ = deg2rad(δ_poly + δ_pert)

    dδ_pert_per_century =
            -(0.000050 * sin(Ja) * deg2rad(4850.4046)) -
             (0.000404 * sin(Jb) * deg2rad(1191.9605)) -
             (0.000617 * sin(Jc) * deg2rad( 262.5475)) +
             (0.000013 * sin(Jd) * deg2rad(6070.2476)) -
             (0.000926 * sin(Je) * deg2rad(  64.3000))
    dδ₀_dt = deg2rad(0.002413 + dδ_pert_per_century) / _SECONDS_PER_CENTURY

    W  = deg2rad(284.95 + 870.5360000 * d)
    dW_dt = deg2rad(870.5360000) / _SECONDS_PER_DAY

    return IAU2015Orientation(α₀, δ₀, W, dα₀_dt, dδ₀_dt, dW_dt)
end

function _iau2015_orientation_saturn(jd::Real)
    d = jd - _J2000_TDB_JD
    T = d / 36525.0

    α₀ = deg2rad(40.589 - 0.036 * T)
    δ₀ = deg2rad(83.537 - 0.004 * T)
    W  = deg2rad(38.90 + 810.7939024 * d)

    return IAU2015Orientation(
        α₀, δ₀, W,
        deg2rad(-0.036) / _SECONDS_PER_CENTURY,
        deg2rad(-0.004) / _SECONDS_PER_CENTURY,
        deg2rad(810.7939024) / _SECONDS_PER_DAY,
    )
end

function _iau2015_orientation_uranus(jd::Real)
    d = jd - _J2000_TDB_JD
    α₀ = deg2rad(257.311)
    δ₀ = deg2rad(-15.175)
    W  = deg2rad(203.81 - 501.1600928 * d)
    return IAU2015Orientation(
        α₀, δ₀, W,
        0.0, 0.0,
        deg2rad(-501.1600928) / _SECONDS_PER_DAY,
    )
end

function _iau2015_orientation_neptune(jd::Real)
    d = jd - _J2000_TDB_JD
    T = d / 36525.0

    N        = deg2rad(357.85 + 52.316 * T)         # Neptune's precession angle
    dN_dt    = deg2rad(52.316) / _SECONDS_PER_CENTURY  # [rad/s]

    α₀     = deg2rad(299.36 + 0.70 * sin(N))
    δ₀     = deg2rad(43.46  - 0.51 * cos(N))
    W      = deg2rad(249.978 + 541.1397757 * d - 0.48 * sin(N))

    dα₀_dt = deg2rad(0.70) * cos(N) * dN_dt
    dδ₀_dt = deg2rad(0.51) * sin(N) * dN_dt
    dW_dt  = deg2rad(541.1397757) / _SECONDS_PER_DAY - deg2rad(0.48) * cos(N) * dN_dt

    return IAU2015Orientation(α₀, δ₀, W, dα₀_dt, dδ₀_dt, dW_dt)
end

function _iau2015_orientation_pluto(jd::Real)
    d = jd - _J2000_TDB_JD
    α₀ = deg2rad(132.993)
    δ₀ = deg2rad(-6.163)
    W  = deg2rad(302.695 + 56.3625225 * d)
    return IAU2015Orientation(
        α₀, δ₀, W,
        0.0, 0.0,
        deg2rad(56.3625225) / _SECONDS_PER_DAY,
    )
end

# =============================================================================
# Public dispatch
# =============================================================================

"""
    iau2015_orientation(body::CelestialBody, jd::Real) → IAU2015Orientation

Return the IAU 2015 orientation parameters for `body` at TDB Julian date `jd`.

Supported bodies: Sun, Mercury, Venus, Mars, Jupiter, Saturn, Uranus,
Neptune, Pluto.

# Not supported
- **Earth (naifid 399).** Earth rotation is model-family dependent (FK5 IAU 1980
  vs IAU 2006/2010) and requires EOP corrections; the IAU 2015 planet-rotation
  format is not the right primitive. Use the Earth-frame edges in AstroFrames
  (`ICRS ↔ GCRS ↔ CIRS ↔ TIRS ↔ ITRS`, or the FK5 chain).
- **Moon (naifid 301).** The IAU stopped publishing a simplified lunar rotation
  model; production lunar work uses SPICE PCK files with full physical
  librations (`MOON_PA`, `MOON_ME`). Use AstroFrames' SPICE-backed Moon-frame
  edges.

# Arguments
- `body::CelestialBody`
- `jd::Real` — TDB Julian date. Callers are responsible for converting from
  other time scales.

# Errors
Throws `ArgumentError` for Earth and Moon (see above) and for any other body
without a registered IAU 2015 implementation.

# Reference
Archinal et al., "Report of the IAU Working Group on Cartographic Coordinates
and Rotational Elements: 2015" (Cel. Mech. Dyn. Astron. 130:22, 2018).

# Returns
An [`IAU2015Orientation`](@ref) containing pole, prime-meridian, and rate values.

# Example
```jldoctest
julia> using AstroUniverse

julia> iau2015_orientation(mars, 2451545.0) isa IAU2015Orientation
true
```
"""
iau2015_orientation(body::CelestialBody, jd::Real) = iau2015_orientation(body.naifid, jd)

function iau2015_orientation(n::Integer, jd::Real)
    if     n == 10 ; return _iau2015_orientation_sun(jd)
    elseif n == 199; return _iau2015_orientation_mercury(jd)
    elseif n == 299; return _iau2015_orientation_venus(jd)
    elseif n == 399
        throw(ArgumentError(
            "Earth rotation is not covered by the IAU 2015 planet-rotation format. " *
            "Use IAU 2006/2010 (with EOP) or FK5 (with EOP) Earth-frame edges in AstroFrames."))
    elseif n == 301
        throw(ArgumentError(
            "Moon rotation is not covered by the IAU 2015 planet-rotation format. " *
            "Use SPICE PCK-backed Moon frames (MOON_PA / MOON_ME) via AstroFrames."))
    elseif n == 499; return _iau2015_orientation_mars(jd)
    elseif n == 599; return _iau2015_orientation_jupiter(jd)
    elseif n == 699; return _iau2015_orientation_saturn(jd)
    elseif n == 799; return _iau2015_orientation_uranus(jd)
    elseif n == 899; return _iau2015_orientation_neptune(jd)
    elseif n == 999; return _iau2015_orientation_pluto(jd)
    else
        throw(ArgumentError(
            "IAU 2015 orientation not implemented for NAIF ID $(n)."))
    end
end
