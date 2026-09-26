```@meta
CurrentModule = AstroSolve
```

# AstroSolve

AstroSolve provides parameter optimization, optimal control, and orbit estimation capablity. Parameter
optimization adjusts maneuver components, spacecraft states, epochs, model parameters, and other
finite sets of values to meet mission constraints. The optimal-control subsystem supports methods
including Hermite-Simpson and Legendre-Gauss-Lobatto collocation, Sims-Flanagan low-thrust
optimization, EMTG's MGA-nDSM formulation, and zero-order-hold finite-burn multiple shooting.
Legendre-Gauss-Lobatto collocation and zero-order-hold multiple shooting are available in
Enterprise.

The estimation subsystem provides batch least squares, an extended Kalman filter with
UDU-factorized covariance, and Rauch-Tung-Striebel smoothing. AstroSolve models two-way range and
Doppler measurements, simulates tracking data, and reads and writes CCSDS Tracking Data Messages.

All three problem types use the same concepts for variables, contraints, and objectives wherever 
logical to configure and solve problems. Trajectories are directed acyclic graphs of events and
intervals, following the approach used in NASA's Copernicus system. A trajectory may contain
propagation events, impulsive manuevers, optimal-control phases, branches, merges, or several 
transcription methods.

Partial derivatives may be supplied analytically or computed with automatic differentiation, and
both sources may be used within one problem. Parameter-optimization sequences use finite
differences.

## Installation

To install the latest version of AstroSovle, first add the local registry (the app store, for those unfamiliar with Julia), then install as usual:

```julia
using Pkg
Pkg.Registry.add(RegistrySpec(url = "https://github.com/GenAstro/GenAstroRegistry.git"))
Pkg.add("AstroSolve")
```
!!! note
    Some packages originally registered in the Julia General registry, including AstroSolve, have moved to the GenAstro local registry. If you do not add the local registry as shown above, you will install only the first MVP release of AstroSolve.

## Quick Start

This example below solves the Hohmann transfer from a 7,000 km circular orbit to
geostationary radius. `Vary` identifies the burn component the solver may change, `Constraint`
sets the radius that must be reached at apoapsis, and `solve!` runs the optimization.
Examples that use formal transcriptions, and estimation examples, are documented in later sections. 

```julia
using Epicycle

sat  = Spacecraft(state = KeplerianState(7000.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                  time  = Time("2020-09-21T12:23:12", TAI(), ISOT()), name = "Sat")
prop = OrbitPropagator(ForceModel(PointMassGravity(earth, ())),
                       IntegratorConfig(DP8(); abstol = 1e-11, reltol = 1e-11, dt = 300.0))
toi  = ImpulsiveManeuver(axes = VNB(), element1 = 1.0)

# Write the flight sequence, in the order it flies
result = target!(method = Optimize(derivatives = :fd, print_level = 0)) do

    # Vary and apply the maneuver
    Vary(delta_v, toi; lower_bound = [0.0, 0.0, 0.0], upper_bound = [3.0, 0.0, 0.0])
    maneuver!(sat, toi)

    # Coast to apoapsis and constrain the radius there
    propagate!(prop, sat, StopAt(position_dot_velocity, sat; equals = 0.0, direction = -1))
    Constraint(position_magnitude, sat; equals = 42164.0)
end

result.info, delta_v(toi)[1]        # :Solve_Succeeded, 2.336796 km/s
```

In the REPL, `?` enters help mode: `?Vary` gives every way a variable is declared, and
`?Constraint` gives the forms of `at =`.

## Parameter optimization

Parameter optimization adjusts a finite set of variables to meet a set of mission contraints. 
It supports problems such as targeting an apogee, designing a maneuver sequence,
selecting an epoch, or identifying a model parameter.  In Epicycle, paramater optimization currenty 
uses finite differencing for partial derivatives. 

The [Parameter optimization](optimization.md) guide develops the quick-start problem, then writes a
three-burn GEO transfer in both supported forms including an Event graph and a simple Domain Specific Language
simlar to GMAT's Target command. A `target!` defines the event sequence `Event` sequences that defines the
problem structure including branching and merging elements. The guide also covers bounds, scaling, reports, and convergence.

## Optimal control

Optimal control determines the state and control histories over one or more trajectory phases. It
supports boundary and path constraints, terminal and integrated objectives, linked phases, and
problems that combine collocation and shooting methods.

The [Optimal control](optimal_control.md) guide begins with a collocation solution of the
brachistochrone, then introduces shooting with a Sims-Flanagan interplanetary transfer. It covers
states, controls, dynamics, bounds, boundary and path constraints, objectives, derivatives, and
phase links. A mixed-transcription example shows how different methods can be used within one
trajectory. 

The section concludes with a library of solved examples. 

## Estimation

AstroSolve estimates a spacecraft state and other properties from measurement data. It provides
batch processing for a complete tracking arc, sequential updates as observations arrive, and
smoothing after the arc is complete.

The [Estimation](estimation.md) guide builds a tracking problem from simulated range and Doppler
data. It first fits the complete arc with batch least squares, then processes the same kind of data
sequentially with an extended Kalman filter and smoother. The guide covers a priori
covariance, measurement noise, process noise, CCSDS tracking data, residuals, and solution
covariance.

## Where to look next

- [Concepts](concepts.md) explains the graph, what nodes and edges carry, and where the solver's
  variables and equations come from.
- [Parameter optimization](optimization.md) covers variables, constraints, and the two ways to
  write a sequence of maneuvers and propagations.
- [Optimal control](optimal_control.md) covers phases, transcriptions, and control histories.
- [Estimation](estimation.md) covers orbit determination, filtering, and smoothing.
- [API reference](api.md) lists the exported names and presents the most common interfaces first.

Docstrings carry the full API reference.
