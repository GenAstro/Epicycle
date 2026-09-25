```@meta
CurrentModule = AstroSolve
```

# AstroSolve

AstroSolve provides tools for parameter optimization, optimal control, and orbit estimation.
Parameter optimization solves for maneuvers, states, epochs, and other adjustable quantities along a
propagated trajectory. Optimal control supports Hermite-Simpson and Legendre-Gauss-Lobatto
collocation, Sims-Flanagan low-thrust transcription, EMTG's MGA-nDSM transcription, and
zero-order-hold finite-burn multiple shooting. Orbit estimation includes batch least squares, an
extended Kalman filter with UDU-factorized covariance, and a Rauch-Tung-Striebel smoother.
Legendre-Gauss-Lobatto collocation and zero-order-hold multiple shooting are available in
Enterprise.

AstroSolve represents a trajectory as a directed acyclic graph of events and intervals, following
the approach used in NASA's Copernicus system. Optimization and estimation use the same forms for
declaring variables and constraints.

Partial derivatives can be supplied analytically or computed with automatic differentiation, and a
problem can use both sources. Each phase selects its own transcription, so one trajectory can
combine collocation and shooting methods. Parameter-optimization sequences use finite differences.

References:

- Betts, J. T. (2010), *Practical Methods for Optimal Control and Estimation Using Nonlinear
  Programming*, 2nd ed., SIAM.
- Bryson, A. E., and Ho, Y.-C. (1975), *Applied Optimal Control*, Hemisphere.
- Hughes, S. P. (2026), "A Transcription-Agnostic Formulation for Optimal Control and Estimation
  in Astrodynamics," AAS/AIAA Astrodynamics Specialist Conference, Vancouver, British Columbia,
  July 2026.
- Sims, J., and Flanagan, S. (1999), "Preliminary Design of Low-Thrust Interplanetary Missions,"
  AAS/AIAA Astrodynamics Specialist Conference, AAS 99-338.
- Tapley, B. D., Schutz, B. E., and Born, G. H. (2004), *Statistical Orbit Determination*, Elsevier.
- Williams, J., Falck, R., and Beekman, I. (2018), "Application of Modern Fortran to Spacecraft
  Trajectory Design and Optimization," AIAA/AAS Space Flight Mechanics Meeting.
  [Online](https://ntrs.nasa.gov/api/citations/20180000413/downloads/20180000413.pdf)

## Installation

Versions through 0.4.0 are in Julia's General registry. From the next version AstroSolve is
released under the Gen Astro Source Available License, which General does not carry, so later
versions come from the Gen Astro registry. Add it once, then install as usual:

```julia
using Pkg
Pkg.Registry.add(RegistrySpec(url = "https://github.com/GenAstro/GenAstro.git"))
Pkg.add("AstroSolve")
```

General is still required, since these packages depend on packages registered there. Installing
without the Gen Astro registry resolves to 0.4.0, the last version General carries, and
reports nothing about the newer ones.

## Quick Start

A targeting problem uses three functions. `Vary` declares what the solver may change,
`Constraint` declares what must hold, and `solve!` solves the problem. This example varies one burn
and constrains the radius reached after the coast. The result is the first burn of a Hohmann
transfer from 7,000 km to geostationary radius.

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

The unknowns are a finite set of values with propagation between them: the components of a
maneuver, a state at epoch, or a launch date. A sequence is written either as a `target!` block that
reads in flight order, as above, or as `Event`s assembled into a `Sequence` when the shape of the
problem is decided while it is being built. The block is a layer over the events and reaches the
same solver.

[Parameter optimization](optimization.md) covers variables and bounds, where a constraint may be
placed, both ways of writing a sequence, and what to check when one does not converge. Five worked
examples run from a single targeted burn to a three-burn GEO transfer written both ways.

## Optimal control

The unknown is a control history rather than a few numbers, so the state and control at every point
of a discretized arc become variables. A phase holds the dynamics, the state and control types, the
span, and the transcription that discretizes it.

<!-- doc-fragment -->
```julia
phase = CollocationPhase(name = :L1_to_L2, transcription = HermiteSimpson(n_steps = 50),
                         dynamics = cr3bp!, model = μ, state = CRState, control = CRControl,
                         tspan = (0.0, tf))
```

Open Epicycle ships Hermite-Simpson collocation, Sims-Flanagan, and MGA with deep-space
maneuvers; Enterprise adds Legendre-Gauss-Lobatto collocation and zero-order-hold multiple
shooting. Naming a different transcription changes nothing else in the problem, and one trajectory
may use several.

[Optimal control](optimal_control.md) covers the dynamics signature, declared partials and the
automatic-differentiation fallback, each transcription's own keywords, where constraints and
objectives attach, and linking phases. Eleven worked examples run, four of them classical problems
with published answers.

## Estimation

The unknown is a state that was never measured directly, and the data are range and Doppler
observations taken from the ground. The same `Vary` declares what is estimated, with the covariance
that says how well it is known going in and how well the data determined it coming out.

<!-- doc-fragment -->
```julia
y0 = Vary(state, sat; guess = guess,
          covariance = Diagonal([1e2, 1e2, 1e2, 1e-2, 1e-2, 1e-2]))

fit = solve!(problem, records; method = Batch(n_iters = 10, tol = 1e-9))
```

`Batch` fits the whole arc at once. `Sequential` runs an extended Kalman filter forward through the
observations and sweeps a Rauch-Tung-Striebel smoother back over them.

[Estimation](estimation.md) covers measurements and the units that carry them, simulating and
reading tracking data, both estimators, process noise, and what comes back. Three worked examples
run.

## Where to look next

- [Concepts](concepts.md) explains the graph, what nodes and edges carry, and where the solver's
  variables and equations come from.
- [Parameter optimization](optimization.md) covers variables, constraints, and the two ways to
  write a sequence of maneuvers and propagations.
- [Optimal control](optimal_control.md) covers phases, transcriptions, and control histories.
- [Estimation](estimation.md) covers orbit determination, filtering, and smoothing.
- [API reference](api.md) lists the exported names and presents the most common interfaces first.

Docstrings carry the rest.
