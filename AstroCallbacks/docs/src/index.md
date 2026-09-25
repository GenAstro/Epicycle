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
versions come from the Gen Astro registry. Add it once, then install as usual:

```julia
using Pkg
Pkg.Registry.add(RegistrySpec(url = "https://github.com/GenAstro/GenAstro.git"))
Pkg.add("AstroCallbacks")
```

General is still required, since these packages depend on packages registered there. Installing
without the Gen Astro registry resolves to 0.4.0, the last version General carries, and
reports nothing about the newer ones.

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

# Define a quantity held unapplied for later execution
# This is used for constraints and where execution should be deferred 
# to the desired point in a simulation.  
q = Calc(inclination, sat, EarthMJ2000Ec)
q()
```

## Quantities Overview

A quantity in Epicycle is a function of the thing it describes rather than a method on it. Where an object-oriented library would offer `sat.getSMA()` and `sat.setSMA(value)`, Epicycle reads a spacecraft's semi-major axis with `semi_major_axis(sat)` and sets it with `semi_major_axis!(sat; to = value)`. The `!` is the Julia convention marking a function that changes one of its arguments. Written this way, a new spacecraft or state type supplies its own `semi_major_axis` and every script that already calls that name keeps working, with no change to Epicycle and no recompilation of the original struct.

In Epicycle, a quantity like `semi_major_axis` is a function of a subject, such as a `Spacecraft` named `sat` (i.e. `semi_major_axis(sat)`).  A subject may be a `Spacecraft`, a maneuver, or a
`CelestialBody` for example. Frame-dependent quantities accept a coordinate system after the
subject, as in `inclination(sat, EarthMJ2000Ec)`; without one, they use the
subject's coordinate system.

A settable quantity is set with its name followed by `!`, with the value as the keyword `to`: `semi_major_axis!(sat; to = 7200.0)`. `is_settable(subject, quantity)` answers whether a pair can be written.

Every quantity takes its subject first. Quantities that depend on a frame take an optional coordinate system second, and an optional `params` third for an orbit-relative coordinate system whose origin does not carry the reference orbit.  The table below illustrates a few built-in quantities. 

| Quantity | Subject | Units | Takes a coordinate system | Settable |
|---|---|---|---|---|
| `position_vector`, `velocity_vector` | spacecraft or coordinate | km, km/s | yes | on a spacecraft |
| `position_x`, `position_y`, `position_z` | spacecraft or coordinate | km | yes | on a spacecraft |
| `position_magnitude`, `velocity_magnitude` | spacecraft or coordinate | km, km/s | yes | on a spacecraft |
| `semi_major_axis` | spacecraft or coordinate | km | yes | on a spacecraft |
| `eccentricity` | spacecraft or coordinate | none | yes | on a spacecraft |
| `inclination`, `raan`, `argument_of_periapsis`, `true_anomaly` | spacecraft or coordinate | rad | yes | on a spacecraft |
| `mean_long_sma` | spacecraft or coordinate | km | yes | on a spacecraft |
| `outgoing_rla` | spacecraft or coordinate on a hyperbolic orbit | rad | yes | no |
| `state` | spacecraft | km, km/s | no | yes |
| `epoch` | spacecraft or coordinate | `Time` | no | on a spacecraft |
| `delta_v` | maneuver | km/s | no | yes |
| `delta_v_magnitude` | maneuver | km/s | no | yes |
| `gravitational_parameter` | celestial body | km³/s² | no | yes |

For a complete list of implemented quantities (you can easily write your own don't forget!) you can use the quantities function as illustrated below.

```julia
# List the quantities implemented for various structs
quantities(Spacecraft)
quantities(ImpulsiveManeuver)
quantities(CelestialBody)
```

!!! note "Extensible Property Access"

    A struct's fields are always reachable directly, as `sat.state`, and a mutable field can be
    assigned with `sat.state = value`. Reaching a value that way is fixed by how the struct was
    defined, so it cannot be extended without editing Epicycle itself. The quantity functions,
    `inclination(sat)` and `inclination!(sat; to = value)`, are the extension point: a third party
    adds a method for their own type and the core structs stay as they are.

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

## Defining a quantity for later evaluation

There are times when you need to define a quantity to be evaluated later. For example, you may want to report a quantity or plot it during or after a simulation. The code `inclination(sat)` computes the quantity at the point where it appears in the code. A `Calc` allows you to define a quantity to be evaluated later for use in stopping conditions, plots, reports, solver variables, and constraints.

`Calc(quantity, subject, deps...)` holds a quantity together with its subject and coordinate system, unapplied. Calling it evaluates the quantity. 

```julia

# Define a calc for a Spacecraft's Z-component of postion in EarthMJ2000Ec coordinates
sat = Spacecraft()
q = Calc(position_z, sat, EarthMJ2000Ec)

# Evaluate the calc
q()                                  

# Apply the the same quantity and frame to a different Spacecraft
other = Spacecraft(state = CartesianState([0.0, 7000.0, 0.0, -5.3, 0.0, 5.3]),
                   time  = Time("2024-01-01T12:00:00", UTC(), ISOT()))
reapply(q, other)                    
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
