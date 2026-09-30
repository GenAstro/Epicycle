# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# GMAT's Earth-fixed axes, as an orientation model for the `orientation` keyword of HarmonicGravity
# and AtmosphericDrag.
#
# GMAT propagates in EarthMJ2000Eq and reaches its Earth-fixed frame by the IAU-76/FK5 chain with
# polar motion and UT1 but without the IERS celestial-pole corrections δΔψ and δΔε, and treats its
# inertial axes as FK5 mean J2000, with no frame bias to ICRF. SatelliteToolbox's
# `r_eci_to_ecef(J2000(), ITRF(), jd_utc, eop_iau1980)` is that chain.
#
# Epicycle's own Earth axes are the frame theory's ITRF chain from ICRF, which applies the
# corrections and the bias. Measured on the one-day LEO case the GMAT comparisons run, that moves
# the final state 12 cm under IAU2006 and 17 cm under FK5 from where GMAT's axes put it. Those
# are modelling differences, not errors, so the GMAT comparisons evaluate gravity and drag in GMAT's axes
# and hold the tight tolerance they had.

using AstroUniverse: AbstractOrientationModel
import AstroUniverse: body_axes_rotation

struct GmatEarthAxes <: AbstractOrientationModel end

const _GMAT_EOP = AstroProp.fetch_iers_eop()

function body_axes_rotation(::GmatEarthAxes, naifid::Integer, jd_tdb::Real)
    jd_utc = Time(jd_tdb, zero(jd_tdb), TDB(), JD()).utc.jd
    R = AstroProp.r_eci_to_ecef(AstroProp.J2000(), AstroProp.ITRF(), jd_utc, _GMAT_EOP)
    # The rate block gives the atmosphere GMAT's spin, constant, about the inertial z axis:
    # with Ṙ = −R [ω×], drag's Rᵀ(R v + Ṙ r) is v − ω × r. Gravity does not read it.
    ω = AstroProp.EARTH_ANGULAR_SPEED
    W = [0.0 -ω 0.0; ω 0.0 0.0; 0.0 0.0 0.0]
    Z = zero(R)
    return AstroProp.SMatrix{6,6}([R Z; -R * W R])
end
