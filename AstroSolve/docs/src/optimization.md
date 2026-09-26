```@meta
CurrentModule = AstroSolve
```

# Parameter Optimization

Parameter optimization changes a finite set of values to satisfy mission constraints. The values
may be maneuver components, epochs, spacecraft states, or model parameters. A problem may also
include an objective, but an objective is not required for targeting.

AstroSolve provides two interfaces for these problems. The event-sequence interface exposes the
trajectory graph directly through `Event` and `Sequence`. The Domain Specific Language records a
linear sequence from the operations inside a `target!` block. Both interfaces produce the same
kind of optimization problem.

## GEO Transfer Problem

The example targets a three-burn transfer to geostationary orbit. Transfer-orbit insertion raises
apoapsis to 85,000 km. A correction at the next equatorial crossing lowers the inclination and
places periapsis at geostationary radius. Mission-orbit insertion then targets the final
semi-major axis.

The solver varies the tangential component of the first and final maneuvers and the tangential and
normal components of the correction. Propagation stops at equatorial crossings, apoapsis, and
periapsis so that each constraint is evaluated at the intended point in the trajectory.

## Configure the Transfer

Both interfaces use the same spacecraft, force model, maneuver guesses, and stopping conditions.
The setup is wrapped in a function because a solve changes the spacecraft and maneuver objects.
Calling the function again gives the second formulation a fresh set of inputs.

```julia
using Epicycle

function geo_transfer_inputs()
    sat = Spacecraft(
        state = CartesianState([
            3737.792, -4607.692, -2845.644,
            5.411, 5.367, -1.566,
        ]),
        time = Time("2000-01-01T11:59:28.000", UTC(), ISOT()),
        name = "GeoSat-1",
    )

    gravity = PointMassGravity(earth, ())
    forces = ForceModel(gravity)
    integ = IntegratorConfig(DP8(); abstol = 1e-12, reltol = 1e-12, dt = 60.0)
    prop = OrbitPropagator(forces, integ)

    toi = ImpulsiveManeuver(
        axes = VNB(), element1 = 2.518, element2 = 0.0, element3 = 0.0,
    )
    mcc = ImpulsiveManeuver(
        axes = VNB(), element1 = 0.559, element2 = 0.588, element3 = 0.0,
    )
    moi = ImpulsiveManeuver(
        axes = VNB(), element1 = 0.282, element2 = 0.0, element3 = 0.0,
    )

    z_crossing = StopAt(position_z, sat, EarthMJ2000Eq; equals = 0.0)
    apoapsis = StopAt(position_dot_velocity, sat; equals = 0.0, direction = -1)
    periapsis = StopAt(position_dot_velocity, sat; equals = 0.0, direction = +1)

    return (; sat, prop, toi, mcc, moi, z_crossing, apoapsis, periapsis)
end
```

## Build the Event Sequence

An `Event` keeps an action with the variables and constraints that apply there. Variables are set
before the action runs, and constraints are evaluated afterward. A maneuver event therefore owns
the components varied for that maneuver, while a propagation event owns constraints on the state
at its stopping condition.

Start with fresh inputs and create the events in mission sequence:

```julia
(; sat, prop, toi, mcc, moi, z_crossing, apoapsis, periapsis) = geo_transfer_inputs()

to_equator_1 = Event(
    name = "Propagate to first equatorial crossing",
    event = () -> propagate!(prop, sat, z_crossing),
)

apply_toi = Event(
    name = "Transfer-orbit insertion",
    event = () -> maneuver!(sat, toi),
    vars = [
        Vary(delta_v, toi;
             guess = [2.518, 0.0, 0.0],
             lower_bound = [0.0, 0.0, 0.0],
             upper_bound = [8.0, 0.0, 0.0]),
    ],
)

to_apoapsis = Event(
    name = "Propagate to apoapsis",
    event = () -> propagate!(prop, sat, apoapsis),
    funcs = [Constraint(position_magnitude, sat; equals = 85000.0)],
)

to_periapsis_1 = Event(
    name = "Propagate to first periapsis",
    event = () -> propagate!(prop, sat, periapsis),
)

to_equator_2 = Event(
    name = "Propagate to second equatorial crossing",
    event = () -> propagate!(prop, sat, z_crossing),
)

apply_mcc = Event(
    name = "Mid-course correction",
    event = () -> maneuver!(sat, mcc),
    vars = [
        Vary(delta_v, mcc;
             guess = [0.559, 0.588, 0.0],
             lower_bound = [-1.0, -1.0, -0.001],
             upper_bound = [4.0, 1.0, 0.001]),
    ],
)

to_periapsis_2 = Event(
    name = "Propagate to second periapsis",
    event = () -> propagate!(prop, sat, periapsis),
    funcs = [
        Constraint(inclination, sat, EarthMJ2000Eq; equals = deg2rad(2.0)),
        Constraint(position_magnitude, sat; equals = 42195.0),
    ],
)

apply_moi = Event(
    name = "Mission-orbit insertion",
    event = () -> maneuver!(sat, moi),
    vars = [
        Vary(delta_v, moi;
             guess = [0.282, 0.0, 0.0],
             lower_bound = [-1.0, -0.001, -0.001],
             upper_bound = [4.0, 0.001, 0.001]),
    ],
    funcs = [Constraint(semi_major_axis, sat; equals = 42166.90)],
)
```

Each `Vary` supplies an initial guess and component-wise bounds. Equal lower and upper bounds hold a
component fixed, as they do for the normal and binormal components of transfer-orbit insertion.
The correction and final insertion allow small binormal components instead. The optional `scale`
keyword sets solver scaling when the natural magnitudes of varied or constrained quantities differ.

`add_sequence!` connects the events as a linear path through the trajectory graph:

```julia
seq = Sequence()
add_sequence!(seq,
    to_equator_1,
    apply_toi,
    to_apoapsis,
    to_periapsis_1,
    to_equator_2,
    apply_mcc,
    to_periapsis_2,
    apply_moi,
)
```

For a graph that branches or merges, use `add_events!` to give each event its dependencies instead
of placing every event in one linear sequence.

## Solve and Check the Sequence

`Optimize` configures the nonlinear solve. Event-sequence problems use finite-difference
derivatives because their propagation and maneuver actions do not provide declared partials.

```julia
result = solve!(seq; method = Optimize(
    max_iter = 1000,
    tol = 1e-6,
    derivatives = :fd,
    print_level = 5,
))

report_sequence(seq)
report_solution(seq, result)
```

`result.info` reports why the solver stopped. A returned result may contain the last iterate from
an unsuccessful solve, so the status and targeted quantities must both be checked. After the solve,
the maneuver objects contain their solved delta-v values and the spacecraft contains the state
produced by the final event.

```julia
println("status            : ", result.info)
println("TOI delta-v (km/s): ", round(delta_v(toi)[1], digits = 6))
println("MCC delta-v (km/s): ", round.(delta_v(mcc)[1:2], digits = 6))
println("MOI delta-v (km/s): ", round(delta_v(moi)[1], digits = 6))
println("final SMA (km)    : ", round(semi_major_axis(sat), digits = 3))
```

## A Simplified Interface

`target!` is a Domain Specific Language for linear parameter-optimization sequences. Operations are
written directly inside the block. A `Vary` applies to the operation immediately after it, and a
`Constraint` evaluates the state produced by the operation immediately before it.

Create fresh inputs before solving the same transfer again:

```julia
(; sat, prop, toi, mcc, moi, z_crossing, apoapsis, periapsis) = geo_transfer_inputs()

result = target!(method = Optimize(derivatives = :fd, print_level = 5)) do
    propagate!(prop, sat, z_crossing)

    Vary(delta_v, toi;
         lower_bound = [0.0, 0.0, 0.0],
         upper_bound = [8.0, 0.0, 0.0])
    maneuver!(sat, toi)

    propagate!(prop, sat, apoapsis)
    Constraint(position_magnitude, sat; equals = 85000.0)

    propagate!(prop, sat, periapsis)
    propagate!(prop, sat, z_crossing)

    Vary(delta_v, mcc;
         lower_bound = [-1.0, -1.0, -0.001],
         upper_bound = [4.0, 1.0, 0.001])
    maneuver!(sat, mcc)

    propagate!(prop, sat, periapsis)
    Constraint(inclination, sat, EarthMJ2000Eq; equals = deg2rad(2.0))
    Constraint(position_magnitude, sat; equals = 42195.0)

    Vary(delta_v, moi;
         lower_bound = [-1.0, -0.001, -0.001],
         upper_bound = [4.0, 0.001, 0.001])
    maneuver!(sat, moi)
    Constraint(semi_major_axis, sat; equals = 42166.90)
end
```

The block records the operations and declarations, builds an event sequence, and solves it. The
result and the updated spacecraft and maneuver objects are read in the same way as they are for an
explicit sequence. Use the explicit interface when the graph branches or merges, or when the
assembled sequence must be inspected before it is solved.
