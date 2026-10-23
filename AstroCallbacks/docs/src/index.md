```@meta
CurrentModule = AstroCallbacks
```

# AstroCallbacks

The AstroCallbacks package provides interfaces for reading and writing
spacecraft, maneuver, and celestial-body quantities among others. Sometimes these are
simple interfaces to simply get and set properties on structs.  Other times, 
they enable computing more complex quantities; for example the right ascension of the outgoing asymptote
of a spacecraft in a desired coordinate system. These quantities are used
throughout Epicycle: AstroProp uses them for stopping conditions, AstroSolve uses them for
optimization/estimation variables and
constraints, and reporting and plotting tools use them for data export and visualization.

## Installation

Versions through 0.4.0 are in Julia's General registry. From the next version AstroCallbacks is
released under the Gen Astro Source Available License, which General does not carry, so later
versions come from the Gen Astro registry. Add both registries once, then install as usual:

```julia
using Pkg
pkg"registry add General https://github.com/GenAstro/GenAstroRegistry.git"
pkg"add AstroCallbacks"
```

General is named in that command for two reasons: the dependencies live there, and Julia
installs it by itself only while no registry is present at all, so adding the Gen Astro
registry alone would leave it out. Installing without the Gen Astro registry resolves to 0.4.0,
the last version General carries, and reports nothing about the newer ones.

## Quick Start

The examples below compute orbital elements in different coordinate systems, compute maneuver and body quantities, set a quantity on a spacecraft and on a maneuver, and define a quantity held unapplied, the way `StopAt` and `Constraint` take it.

```julia
using AstroCallbacks, AstroModels, AstroStates, AstroEpochs, AstroFrames, AstroUniverse, AstroManeuvers

# Create a Spacecraft first
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]),
                 time  = Time("2024-01-01T12:00:00", UTC(), ISOT()))

# Compute quantities of the spacecraft in the spacecraft's coordinates
semi_major_axis(sat)
eccentricity(sat)
inclination(sat)
position_magnitude(sat)

# Compute quantities of the spacecraft in different (in this case built-in) Coordinates
inclination(sat, EarthMJ2000Ec)
raan(sat, EarthMJ2000Ec)
position_z(sat, EarthMJ2000Ec)
velocity_magnitude(sat, EarthFixed) 

# Compute quantities of a maneuver and a celestial body
toi = ImpulsiveManeuver(axes = VNB(), element1 = 0.5)
delta_v(toi)                             # km/s, in the maneuver's axes
gravitational_parameter(earth)           # km³/s²

# Set quantities, which holds the rest of the representation the value belongs to
semi_major_axis!(sat; to = 7200.0)
inclination!(sat, EarthMJ2000Ec; to = deg2rad(30.0))
delta_v!(toi; to = [0.5, 0.0, 0.0])

```

## Quantities Overview

A **quantity** is a value you can read from an object in Epicycle, such as a spacecraft's inclination, a maneuver's delta-v, or a celestial body's gravitational parameter. The object is called the **subject**. Quantity functions provide a common interface for reading these values and, where supported, changing them.

### Reading quantities

Call a quantity function with its subject as the first argument:

```julia
semi_major_axis(sat)                 # Spacecraft semi-major axis, in km
delta_v(toi)                         # Maneuver delta-v vector, in km/s
gravitational_parameter(earth)       # Earth's gravitational parameter, in km³/s²
```

A quantity may return a scalar, a vector, or another value such as a `Time`. The reference table below lists the units for each quantity; angles are in radians.

### Choosing a coordinate system or body

For frame-dependent quantities, pass the coordinate system as the second argument. If you omit it, the quantity uses the subject's coordinate system.

```julia
inclination(sat)                     # In the spacecraft's coordinate system
inclination(sat, EarthMJ2000Ec)       # In EarthMJ2000Ec coordinates
velocity_vector(sat, EarthFixed)      # In Earth-fixed coordinates, in km/s
```

The coordinate system determines the origin and axes used in the calculation. Orbital quantities use the gravitational parameter of that origin. When a coordinate system requires additional data, such as a reference state for an orbit-relative frame, pass a `params` named tuple as the third argument. See the quantity's docstring for the supported arguments.

Some quantities take a celestial body as the second argument instead of a coordinate system. For example, `beta_angle(sat, moon)` uses the Moon; omitting the body uses the central body of the subject's frame.

### Setting quantities

For quantities that support writing, call the function with `!` appended to its name and supply the new value with the `to` keyword:

```julia
semi_major_axis!(sat; to = 7200.0)
inclination!(sat, EarthMJ2000Ec; to = deg2rad(30.0))
delta_v!(toi; to = [0.5, 0.0, 0.0])
```

The `!` follows Julia's convention for functions that modify an argument. Supply values in the same units used by the corresponding reader.

Setting a state element preserves the other elements of the representation used by that setter. This may change values in other representations. Each setter's docstring lists the elements it preserves. For example, `equinoctial_h!` preserves the semi-major axis, `equinoctial_k`, `equinoctial_p`, `equinoctial_q`, and mean longitude.

Use `is_settable(subject, quantity)` to check whether a quantity can be changed on a particular subject:

```julia
is_settable(sat, semi_major_axis)     # true
is_settable(sat, position_dot_velocity) # false
```

## Available quantities

The table summarizes the built-in quantities and their supported subjects. A coordinate subject is an `AstroFrames.Coordinate`. The second argument is optional where shown. The Settable column indicates which subjects support the corresponding `!` function.

Some quantities are defined only for particular orbit types. For example, `orbital_period` is undefined for a hyperbolic orbit and throws an `ArgumentError` explaining the restriction. Check the quantity's docstring for applicability and other requirements.

| Quantity | Subject | Units | Second argument | Settable |
|---|---|---|---|---|
| `position_vector`, `velocity_vector` | spacecraft or coordinate | km, km/s | coordinate system | on a spacecraft |
| `position_x`, `position_y`, `position_z` | spacecraft or coordinate | km | coordinate system | on a spacecraft |
| `position_magnitude`, `velocity_magnitude` | spacecraft or coordinate | km, km/s | coordinate system | on a spacecraft |
| `velocity_x`, `velocity_y`, `velocity_z` | spacecraft or coordinate | km/s | coordinate system | on a spacecraft |
| `right_ascension`, `declination`, `velocity_right_ascension`, `velocity_declination`, `velocity_azimuth`, `flight_path_angle` | spacecraft or coordinate | rad | coordinate system | on a spacecraft |
| `semi_major_axis` | spacecraft or coordinate | km | coordinate system | on a spacecraft |
| `eccentricity` | spacecraft or coordinate | none | coordinate system | on a spacecraft |
| `inclination`, `raan`, `argument_of_periapsis`, `true_anomaly` | spacecraft or coordinate | rad | coordinate system | on a spacecraft |
| `mean_anomaly` | spacecraft or coordinate | rad | coordinate system | on a spacecraft |
| `eccentric_anomaly`, `hyperbolic_anomaly` | spacecraft or coordinate | rad | coordinate system | no |
| `argument_of_latitude`, `longitude_of_periapsis`, `descending_node_right_ascension` | spacecraft or coordinate | rad | coordinate system | no |
| `true_longitude` | spacecraft or coordinate | rad | coordinate system | on a spacecraft |
| `apoapsis_radius`, `periapsis_radius` | spacecraft or coordinate | km | coordinate system | on a spacecraft |
| `semilatus_rectum` | spacecraft or coordinate | km | coordinate system | on a spacecraft |
| `specific_energy`, `c3` | spacecraft or coordinate | km²/s² | coordinate system | no |
| `hyperbolic_excess_speed`, `apoapsis_speed`, `periapsis_speed` | spacecraft or coordinate | km/s | coordinate system | no |
| `angular_momentum`, `angular_momentum_magnitude` | spacecraft or coordinate | km²/s | coordinate system | no |
| `orbit_normal`, `eccentricity_vector` | spacecraft or coordinate | none | coordinate system | no |
| `orbital_angular_velocity`, `mean_motion` | spacecraft or coordinate | rad/s | coordinate system | no |
| `orbital_period` | spacecraft or coordinate | s | coordinate system | no |
| `b_dot_t`, `b_dot_r`, `b_magnitude` | spacecraft or coordinate on a hyperbolic orbit | km | coordinate system | no |
| `b_plane_angle` | spacecraft or coordinate on a hyperbolic orbit | rad | coordinate system | no |
| `equinoctial_h`, `equinoctial_k`, `equinoctial_p`, `equinoctial_q`, `alternate_equinoctial_p`, `alternate_equinoctial_q` | spacecraft or coordinate on a closed orbit | none | coordinate system | on a spacecraft |
| `mean_longitude` | spacecraft or coordinate on a closed orbit | rad | coordinate system | on a spacecraft |
| `modified_equinoctial_f`, `modified_equinoctial_g`, `modified_equinoctial_h`, `modified_equinoctial_k` | spacecraft or coordinate | none | coordinate system | on a spacecraft |
| `delaunay_l`, `delaunay_g`, `delaunay_h` | spacecraft or coordinate on a closed orbit | rad | coordinate system | on a spacecraft |
| `delaunay_L`, `delaunay_G`, `delaunay_H` | spacecraft or coordinate on a closed orbit | km²/s | coordinate system | on a spacecraft |
| `brouwer_long_semi_major_axis`, `brouwer_long_eccentricity`, `brouwer_long_inclination`, `brouwer_long_raan`, `brouwer_long_argument_of_periapsis`, `brouwer_long_mean_anomaly`, and the same six beginning `brouwer_short_` | spacecraft or coordinate | km, none, rad | coordinate system | on a spacecraft |
| `outgoing_asymptote_right_ascension`, `outgoing_asymptote_declination`, `outgoing_b_vector_azimuth`, and the same three beginning `incoming_` | spacecraft or coordinate | rad | coordinate system | on a spacecraft |
| `beta_angle`, `lunar_beta_angle`, `orbit_normal_sun_angle` | spacecraft or coordinate | rad | body | no |
| `sun_earth_spacecraft_angle`, `sun_spacecraft_earth_angle`, `sun_spacecraft_moon_angle`, `moon_earth_sun_angle`, `moon_earth_spacecraft_angle`, `moon_spacecraft_earth_angle` | spacecraft or coordinate | rad | none | no |
| `earth_shadow_fraction`, `moon_shadow_fraction`, from AstroProp | spacecraft or coordinate | none | none | no |
| `state` | spacecraft | km, km/s | none | yes |
| `epoch` | spacecraft or coordinate | `Time` | none | on a spacecraft |
| `delta_v` | maneuver | km/s | none | yes |
| `delta_v_magnitude` | maneuver | km/s | none | yes |
| `gravitational_parameter` | celestial body | km³/s² | none | yes |
Use `quantities` to list the quantity functions available for a subject type:

```julia
quantities(Spacecraft)
quantities(ImpulsiveManeuver)
quantities(CelestialBody)
```

By default, `quantities` searches `AstroCallbacks`. To find quantities defined in your script, pass `Main` as the second argument, for example `quantities(Spacecraft, Main)`.

### Extending the interface

Quantity functions can be extended with methods for your own types. Scripts that call these functions can then use those types through the same interface. You can also define new quantities for existing types without changing their fields. See [Writing your own Quantity](@ref) for details.

To store a quantity and its arguments for evaluation later, use `Calc`, as described in [Defining a quantity for later evaluation](@ref).

## Using Quantities

The examples below illustrate how to evaluate and set various quantities on different structs such as `Spacecraft`, and `ImpulsiveManeuver`.

```julia
using AstroCallbacks, AstroModels, AstroStates, AstroEpochs, AstroFrames, AstroManeuvers

# Create the spacecraft
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]),
                 time  = Time("2024-01-01T12:00:00", UTC(), ISOT()))

# Read the semi-major axis in the spacecraft's own coordinate system
semi_major_axis(sat)                            # km

# Read the inclination in EarthMJ2000Ec
inclination(sat, EarthMJ2000Ec)                 # rad

# Set the semi-major axis, holding the other five orbital elements
semi_major_axis!(sat; to = 7200.0)
semi_major_axis(sat)                            # 7200.0

# Set the inclination as expressed in EarthMJ2000Ec
inclination!(sat, EarthMJ2000Ec; to = deg2rad(30.0))
inclination(sat, EarthMJ2000Ec)                 # 0.5236

# Replace the whole state
state!(sat; to = [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0])

# Set the delta-v on an impulsive maneuver
toi = ImpulsiveManeuver(axes = VNB())
delta_v!(toi; to = [0.5, 0.0, 0.0])             # km/s, in the maneuver's axes
delta_v_magnitude(toi)                          # 0.5

# Ask whether a subject and quantity pair can be written
is_settable(sat, semi_major_axis)               # true
is_settable(sat, position_dot_velocity)         # false
```

## Writing your own Quantity

A quantity function takes its subject as the first argument. Custom quantities
require no registry or subtype and work with `StopAt`, `Constraint`, `Vary`, and
`history`. Five optional declarations provide the same capabilities as the
shipped quantities:

- `label(::typeof(q))`: the name it takes in a report or an error message.
- `AstroCallbacks.subject_type(::typeof(q))`: the type of subject it reads, so `quantities` lists it for that type.
- a method taking a coordinate system after the subject: reads it in another frame. Building it from shipped frame-aware quantities makes it frame-aware too.
- a setter `q!(subject; to)` and one method `set_quantity!(subject, ::typeof(q); to)` that calls it: writes it, so a user can set it and a solver can vary it. `is_settable(subject, q)` then returns `true`. The setter takes the same coordinate system argument as the reader, if the reader takes one.
- `output_partial(subject, ::typeof(q))`: its derivative with respect to the subject's state, so a solver uses it instead of automatic differentiation.

```julia
# Geocentric latitude, in radians, in any coordinate system; the spacecraft's own by default
geocentric_latitude(sc, cs) = asin(position_z(sc, cs) / position_magnitude(sc, cs))
geocentric_latitude(sc)     = geocentric_latitude(sc, frame_of(sc))
AstroCallbacks.label(::typeof(geocentric_latitude)) = "Geocentric latitude"
AstroCallbacks.subject_type(::typeof(geocentric_latitude)) = Spacecraft

geocentric_latitude(sat, EarthMJ2000Ec)
quantities(Spacecraft, Main)             # the quantities this script has defined

# A settable quantity: the drag area on the spacecraft's drag geometry, in m²
drag_area(sc) = sc.drag.drag_area
drag_area!(sc::Spacecraft; to::Real) =
    (sc.drag = SphericalDrag(c_d = sc.drag.c_d, drag_area = to); sc)
AstroCallbacks.set_quantity!(sc::Spacecraft, ::typeof(drag_area); to) = drag_area!(sc; to = to)

dragged = Spacecraft(drag = SphericalDrag(c_d = 2.2, drag_area = 4.0))
is_settable(dragged, drag_area)          # true, so AstroSolve can vary it
drag_area!(dragged; to = 5.0)
drag_area(dragged)
```

## Legacy Calc Interface

!!! note "Deprecation Notice"

    This interface will be deprecated in a future release.

AstroCallbacks also provides an older form, the "Calcs": structs that pair a subject with a calculation tag, such as `OrbitCalc(sc, SMA())`. `StopAt` and the solvers still accept it, and existing scripts keep working. It cannot name a coordinate system, so new scripts use the quantities above. The OrbitCalc interface will be deprecated in a future release.

Each supports reading its quantity and, where defined, writing a target value for
optimization.

- **OrbitCalc**: Semi-major axis, eccentricity, inclination, incoming/outgoing asymptotes, periapsis conditions, and other orbital properties
- **BodyCalc**: Celestial body gravitational parameter
- **ManeuverCalc**: Δv components and magnitude

```julia
using AstroCallbacks, AstroStates, AstroModels, AstroEpochs, AstroUniverse, AstroManeuvers

# Create a spacecraft with orbital state
sc = Spacecraft(state = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]),
                time = Time("2024-01-01T12:00:00", UTC(), ISOT()),
                mass = 1000.0)

# Get semi-major axis from current state
sma_calc = OrbitCalc(sc, SMA())
a = get_calc(sma_calc)
set_calc!(sma_calc, 10000.0)

# Set target incoming asymptote (rp = 6900, C3 = 14.0)
hyp = OrbitCalc(sc, IncomingAsymptote())
set_calc!(hyp, [6900.0, 14.0, 0.0, 0.0, 0.0, 0.0])

# Set and get Earth's mu
mu_calc = BodyCalc(earth, GravParam())
μ = get_calc(mu_calc)
set_calc!(mu_calc, 3.986e5)

# Set and get maneuver elements
toi = ImpulsiveManeuver()
dvvec_calc = ManeuverCalc(toi, sc, DeltaVVector())
Δv = get_calc(dvvec_calc)
set_calc!(dvvec_calc, [0.2, 0.3, 0.4])
```

```@index
```

## API Reference

```@autodocs
Modules = [AstroCallbacks]
Order = [:type, :function, :macro, :constant]
Public = true
Private = false
```
