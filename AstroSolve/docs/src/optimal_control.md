```@meta
CurrentModule = AstroSolve
```

# Optimal Control

Optimal control solves for the state and control histories over one or more trajectory phases.
Each phase selects a transcription, declares the quantities the solver may change, and carries its
own constraints and objectives. AstroSolve supports collocation and shooting transcriptions through
the same `Vary`, `Constraint`, `Objective`, and `Sequence` interface.

The first example uses Hermite-Simpson collocation to introduce the complete phase workflow. The
second uses Sims-Flanagan multiple shooting and focuses on the quantities that differ from a
collocation problem.

## Collocation: The Brachistochrone

The brachistochrone is the minimum-time path between two points under constant gravity. Its state
contains horizontal position, downward vertical position, and speed. The control is the path angle.
The final time is also a decision variable.

### Define the State, Control, and Dynamics

A collocation phase uses concrete state and control types. Their fields define the component order
used by guesses, bounds, and derivative matrices.

```julia
using Epicycle

struct BrachModel
    g::Float64
end

struct BrachState{T} <: AbstractState
    x::T
    y::T
    v::T
end

struct BrachControl{T} <: AbstractControl
    theta::T
end

function brach_dynamics!(dy, y::BrachState, u::BrachControl, p, t, model)
    dy[1] = y.v * sin(u.theta)
    dy[2] = y.v * cos(u.theta)
    dy[3] = model.g * cos(u.theta)
end
```

The dynamics function writes one derivative for each state component. The parameter vector `p` and
time `t` are present in the interface even when a model does not use them.

### Build the Collocation Phase

`CollocationPhase` combines the dynamics with a transcription, model, state and control types, and
time span. `HermiteSimpson(n_steps = 20)` supplies the mesh and the defect equations that enforce
the dynamics between nodes.

```julia
phase = CollocationPhase(
    name = :brachistochrone,
    transcription = HermiteSimpson(n_steps = 20),
    dynamics = brach_dynamics!,
    model = BrachModel(9.80665),
    state = BrachState,
    control = BrachControl,
    tspan = (0.0, 1.0),
)

Vary(state, phase;
     guess = [0.0 0.5 1.0 1.5 2.0
              0.0 0.5 1.0 1.5 2.0
              0.0 1.0 1.5 2.0 2.5],
     lower_bound = [-Inf, -Inf, 0.0],
     upper_bound = [Inf, Inf, Inf])

Vary(control, phase;
     guess = reshape([0.3, 0.5, 0.7, 0.9, 1.0], 1, 5),
     lower_bound = [-pi / 2],
     upper_bound = [pi / 2])

Vary(final_time, phase; lower_bound = 0.1, upper_bound = 10.0)
```

State and control guesses may contain fewer columns than the transcription mesh; AstroSolve
interpolates them to the required nodes. Component bounds apply across the complete history. The
final-time guess comes from the end of `tspan` because `Vary(final_time, ...)` does not replace it.

### Apply Constraints and the Objective

Constraint functions read values from the phase context. `Initial()`, `Final()`, and `Path()` place
the same function at the start, end, or mesh points of the phase.

```julia
start_state(c) = [state(c).x, state(c).y, state(c).v]
final_position(c) = [state(c).x, state(c).y]
speed(c) = [state(c).v]
duration(c) = final_time(c)

Constraint(start_state, phase; equals = [0.0, 0.0, 0.0], at = Initial())
Constraint(final_position, phase; equals = [2.0, 2.0], at = Final())
Constraint(speed, phase; upper_bound = 10.0, at = Path())

Objective(duration, phase; sense = Min())
```

The path placement applies the speed limit at every mesh node. With no `at` keyword, the objective
is a terminal term. An objective at `Path()` is integrated over the phase.

### Supply and Check Derivatives

`@partial` associates a derivative with the function and varied quantity it differentiates. A
dynamics partial writes into the supplied matrix; a constraint or objective partial returns its
matrix.

```julia
@partial(brach_dynamics!, state) do dF, y, u, p, t, model
    dF[1, 3] = sin(u.theta)
    dF[2, 3] = cos(u.theta)
end

@partial(brach_dynamics!, control) do dF, y, u, p, t, model
    dF[1, 1] = y.v * cos(u.theta)
    dF[2, 1] = -y.v * sin(u.theta)
    dF[3, 1] = -model.g * sin(u.theta)
end

@partial(start_state, state) do c
    [1.0 0.0 0.0
     0.0 1.0 0.0
     0.0 0.0 1.0]
end

@partial(final_position, state) do c
    [1.0 0.0 0.0
     0.0 1.0 0.0]
end

@partial(speed, state) do c
    [0.0 0.0 1.0]
end

@partial(duration, final_time) do c
    [1.0]
end

check_partials(phase)
```

`check_partials` compares the declared partials with finite differences. It should be run after the
phase is fully configured and before a long solve.

### Solve and Inspect the Trajectory

A phase is added to a `Sequence` and solved with `Optimize`. The solution remains on the phase, where
the state and control histories have one column per transcription node.

```julia
result = solve!(Sequence(phase); method = Optimize())

println("status: ", result.info)
println("time  : ", round(get_final_time(phase), digits = 6), " s")

times = get_node_times(phase)
states = state(phase)
controls = control(phase)
```

Check `result.info` before using the trajectory. The endpoint accessors
`get_initial_state(phase)` and `get_final_state(phase)` return the declared state type; `state(phase)`
and `control(phase)` return the complete histories.

## Shooting: An Earth-to-Mars Transfer

Sims-Flanagan divides a low-thrust transfer into segments, applies one impulse per segment, and
propagates from both endpoints to a match point. The phase varies segment controls and endpoint
quantities rather than a state and control value at every collocation node.

This example connects circular, coplanar Earth and Mars ephemerides over the Hohmann transfer time.
It minimizes propellant use for a spacecraft with a one-newton engine.

### Build the Shooting Phase

The phase needs an ephemeris at each endpoint and a propulsion model for its internal propagation.
`matchpoint_scale` makes position, velocity, and mass errors comparable in the continuity equations.

```julia
using LinearAlgebra

const MU_SUN = 1.32712440018e11
const AU = 1.495978707e8
const R_EARTH = AU
const R_MARS = 1.524 * AU
const VC_EARTH = sqrt(MU_SUN / R_EARTH)
const VC_MARS = sqrt(MU_SUN / R_MARS)
const OM_EARTH = VC_EARTH / R_EARTH
const OM_MARS = VC_MARS / R_MARS
const M0 = 1500.0
const ISP = 3000.0
const TMAX = 1.0e-3
const G0 = 9.80665e-3
const TOF = pi * sqrt(((R_EARTH + R_MARS) / 2)^3 / MU_SUN)
const THETA_MARS_0 = pi - OM_MARS * TOF

function earth_ephemeris(t::Real)
    theta = OM_EARTH * t
    r = R_EARTH * [cos(theta), sin(theta), 0.0]
    v = VC_EARTH * [-sin(theta), cos(theta), 0.0]
    a = -(MU_SUN / R_EARTH^2) * [cos(theta), sin(theta), 0.0]
    return r, v, a
end

function mars_ephemeris(t::Real)
    theta = THETA_MARS_0 + OM_MARS * t
    r = R_MARS * [cos(theta), sin(theta), 0.0]
    v = VC_MARS * [-sin(theta), cos(theta), 0.0]
    a = -(MU_SUN / R_MARS^2) * [cos(theta), sin(theta), 0.0]
    return r, v, a
end

const N_SEGMENTS = 60
const SMOOTHING = 0.015

transfer = SimsFlanaganPhase(
    name = :earth_to_mars,
    transcription = SimsFlanagan(
        n_segments = N_SEGMENTS,
        throttle_smoothing = SMOOTHING,
    ),
    model = PropulsionModel(mu = MU_SUN, Isp = ISP, Tmax = TMAX, g0 = G0),
    ephemeris_left = earth_ephemeris,
    ephemeris_right = mars_ephemeris,
    tspan = (0.0, TOF),
    matchpoint_scale = [
        R_EARTH, R_EARTH, R_EARTH,
        VC_EARTH, VC_EARTH, VC_EARTH,
        M0,
    ],
)
```

### Define the Shooting Variables

The forward and backward control blocks cover the segments propagated from each endpoint. Initial
mass is fixed with equal bounds, while final mass remains free.

```julia
_, departure_velocity, _ = earth_ephemeris(0.0)
_, arrival_velocity, _ = mars_ephemeris(TOF)

Vary(forward_control, transfer;
     guess = 0.5 * departure_velocity / norm(departure_velocity),
     lower_bound = fill(-2.0, 3),
     upper_bound = fill(2.0, 3))

Vary(backward_control, transfer;
     guess = 0.5 * arrival_velocity / norm(arrival_velocity),
     lower_bound = fill(-2.0, 3),
     upper_bound = fill(2.0, 3))

Vary(initial_mass, transfer; lower_bound = M0, upper_bound = M0)

Vary(final_mass, transfer;
     guess = 0.85 * M0,
     lower_bound = 100.0,
     upper_bound = M0,
     scale = M0)
```

The component bounds permit vectors whose magnitude is greater than one. A path constraint is
therefore required to enforce the throttle limit on every segment.

### Apply the Path Constraint and Objective

The objective uses the same smoothed throttle magnitude as the transcription's mass propagation.
Its partials cover both control blocks.

```julia
thrust_ball(c) = [dot(control(c), control(c))]
Constraint(thrust_ball, transfer;
           lower_bound = 0.0,
           upper_bound = 1.0,
           at = Path())

const KG_PER_THROTTLE = TMAX * (TOF / N_SEGMENTS) / (G0 * ISP)

throttle_magnitudes(u) = sqrt.(vec(sum(abs2, u, dims = 1)) .+ SMOOTHING^2)

propellant(c) = KG_PER_THROTTLE * (
    sum(throttle_magnitudes(forward_control(c))) +
    sum(throttle_magnitudes(backward_control(c)))
)

@partial(propellant, forward_control) do c
    u = forward_control(c)
    KG_PER_THROTTLE .* vec(u ./ throttle_magnitudes(u)')
end

@partial(propellant, backward_control) do c
    u = backward_control(c)
    KG_PER_THROTTLE .* vec(u ./ throttle_magnitudes(u)')
end

Objective(propellant, transfer; sense = Min())
```

### Solve and Inspect the Transfer

The shooting phase uses the same sequence and solve interface as the collocation phase. Its results
are read through the quantities declared with `Vary`.

```julia
check_partials(transfer)
result = solve!(Sequence(transfer);
                method = Optimize(max_iter = 500, print_level = 5))

println("status         : ", result.info)
println("delivered mass : ", round(final_mass(transfer), digits = 3), " kg")
println("propellant     : ", round(M0 - final_mass(transfer), digits = 3), " kg")

u_forward = forward_control(transfer)
u_backward = backward_control(transfer)
```

The phase enforces its match-point equations internally. A successful solver status should still be
followed by checks of the endpoint conditions, throttle constraint, and quantities used downstream.

## Link Multiple Phases

A sequence may contain several phases, and each phase keeps its own transcription. `Link` identifies
two adjacent phases; `Constraint(continuity, link)` requires their boundary state and time to agree.
For shooting phases, the continuity vector also includes mass.

The [mixed-transcription example](https://github.com/GenAstro/Epicycle/blob/main/Epicycle/examples/Ex_MixedTranscription.jl)
uses Sims-Flanagan for the first half of an orbit-raising transfer and Hermite-Simpson for the
second. The phases are added to one `Sequence`, linked by continuity, and solved as one nonlinear
program. A custom constraint on the same link can represent a handoff that is not continuous.

## Solved Examples

The Epicycle examples apply the same interface to several problem classes:

- [Brachistochrone](https://github.com/GenAstro/Epicycle/blob/main/Epicycle/examples/Ex_Brachistochrone.jl): Hermite-Simpson collocation with free final time.
- [Hull problem](https://github.com/GenAstro/Epicycle/blob/main/Epicycle/examples/Ex_HullProblem.jl): terminal and integrated objectives in one phase.
- [Goddard rocket](https://github.com/GenAstro/Epicycle/blob/main/Epicycle/examples/Ex_GoddardRocket.jl): linked collocation phases and a singular arc.
- [Low-thrust orbit raising](https://github.com/GenAstro/Epicycle/blob/main/Epicycle/examples/Ex_LowThrustOrbitRaising.jl): continuous-thrust trajectory optimization.
- [Soft lunar landing](https://github.com/GenAstro/Epicycle/blob/main/Epicycle/examples/Ex_MoonLanding.jl): final-mass maximization with endpoint constraints.
- [Earth-to-Mars transfer](https://github.com/GenAstro/Epicycle/blob/main/Epicycle/examples/Ex_MarsTransferSimsFlanagan.jl): Sims-Flanagan multiple shooting.
- [Earth-to-Apophis rendezvous](https://github.com/GenAstro/Epicycle/blob/main/Epicycle/examples/Ex_ApophisRendezvousSimsFlanagan.jl): low-thrust rendezvous with moving endpoints.
- [Earth-Earth-Venus transfer](https://github.com/GenAstro/Epicycle/blob/main/Epicycle/examples/Ex_GravityAssistMGA.jl): linked MGA-nDSM phases with a flyby constraint.
