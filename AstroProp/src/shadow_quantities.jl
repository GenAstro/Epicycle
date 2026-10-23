# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Shadow as a quantity (#11, section 4): how much of the Sun's disk the Earth or the Moon hides
# from a subject, from one of the eclipse models SolarRadiationPressure takes. Given the model the
# force uses, a report, a constraint or a stopping condition reads the same shadow the force
# applies. The quantity needs no force model; the eclipse model alone defines it.

"""
    earth_shadow_fraction(subject[, model])

Fraction of the Sun's disk hidden by the Earth, seen from `subject`: 0 in full sun, 1 in umbra,
in between in penumbra.

# Arguments
- `subject`: a `Spacecraft`, an `AstroFrames.Coordinate`, or anything else implementing the
  subject interface.
- `model`: the eclipse model, [`DualCone`](@ref) by default or [`SmoothedConical`](@ref). Pass
  the one an SRP force uses to read the shadow that force applies.

# Notes
One minus the model's lighting factor, which [`SolarRadiationPressure`](@ref) multiplies its
acceleration by. No force model is needed. The Sun and the Earth are spheres of their equatorial
radii; the Earth's atmosphere and oblateness are not modeled. The Sun's position comes from the
loaded planetary ephemeris at the subject's epoch.

# Returns
The hidden fraction, from 0 to 1.

# Examples
```julia
sat = Spacecraft(state = CartesianState([-7000.0, 0.0, 0.0, 0.0, -7.5, 0.0]))
earth_shadow_fraction(sat)
earth_shadow_fraction(sat, SmoothedConical())
```
"""
earth_shadow_fraction(subject, model::AbstractShadowModel = DualCone()) =
    _hidden_fraction(subject, earth, model)

"""
    moon_shadow_fraction(subject[, model])

Fraction of the Sun's disk hidden by the Moon, seen from `subject`: 0 in full sun, 1 in umbra,
in between in penumbra.

# Arguments
- `subject`: a `Spacecraft`, an `AstroFrames.Coordinate`, or anything else implementing the
  subject interface.
- `model`: the eclipse model, [`DualCone`](@ref) by default or [`SmoothedConical`](@ref).

# Notes
As [`earth_shadow_fraction`](@ref), with the Moon as the occulting sphere of its equatorial
radius.

# Returns
The hidden fraction, from 0 to 1.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
moon_shadow_fraction(sat)
```
"""
moon_shadow_fraction(subject, model::AbstractShadowModel = DualCone()) =
    _hidden_fraction(subject, moon, model)

"One minus the lighting factor of `subject` behind `occulter`, from eclipse `model`."
function _hidden_fraction(subject, occulter::CelestialBody, model::AbstractShadowModel)
    x = to_vector(CartesianState(Coordinate(subject, CoordinateSystem(occulter, ICRF()))))
    r_sat = SVector{3}(x[1], x[2], x[3])
    r_sun = translate(occulter, sun, epoch_of(subject).tdb.jd)
    return 1 - _shadow_factor(model, r_sat, r_sun, sun.equatorial_radius,
                              occulter.equatorial_radius)
end

for (f, l) in ((:earth_shadow_fraction, "Earth shadow fraction"),
               (:moon_shadow_fraction,  "Moon shadow fraction"))
    @eval EpicycleBase.label(::typeof($f)) = $l
    @eval AstroCallbacks.subject_type(::typeof($f)) = Union{Spacecraft, Coordinate}
end
