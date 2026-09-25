# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

__precompile__()
"""
    module AstroSolve

Tools for trajectory targeting, optimal control, and orbit estimation.

AstroSolve represents trajectories as events and intervals, with shared variable,
constraint, objective, derivative, and solver interfaces across its problem types.
"""
module AstroSolve

using ForwardDiff
using NLsolve
using SNOW
using LinearAlgebra
using Printf
import Random        # check_sparsity samples the variable bounds repeatably
using Dates          # the spec layer stamps a solver log with the time

using EpicycleBase
using AstroManeuvers
using AstroCallbacks
# A spacecraft has a state and so does a phase. Importing the function is what
# makes the phase methods below extend it rather than shadow it.
import AstroCallbacks: state
using AstroModels: Spacecraft, to_posvel, set_posvel!
# The spec vocabulary lets a force model be a phase's right-hand side, and `simulate` flies the
# truth spacecraft between observation epochs.
using AstroProp: ForceModel, accel_eval!, OrbitPropagator, propagate!, StopAt,
                 PropDurationSeconds
using AstroFrames: Coordinate
# `using Dates` above makes `Time` mean Dates.Time here, so an epoch is named AstroEpochs.Time.
import AstroEpochs
using AstroStates: CartesianState

# The export list is what a user writes to set up a problem, solve it and read the answer. The
# rest of the package, including the interface a new transcription, phase type or measurement is
# written against, is imported by name: `using AstroSolve: jacobian_chunk`.

# Targeting with events, maneuvers and propagation.
export SolverVariable, Event, Sequence, add_events!, add_sequence!
export solve_trajectory!, target!, report_sequence, report_solution

# Posing and solving a problem: what varies, what must hold, what is minimized, and how.
export Vary, Constraint, Objective, Link, continuity, solve!
export Optimize, Batch, Sequential, RTS
export Initial, Final, Path, Boundary, Min, Max
export SolveFor, Consider, Fixed
export partial, @partial, check_partials
export check_sparsity, dense_jacobian!, dense_jacobian

# Phases and the transcriptions that fly them.
export CollocationPhase, HermiteSimpson, PropulsionModel
export SimsFlanagan, SimsFlanaganPhase, MGAnDSMs, MGAnDSMsPhase

# What a spec can name: the quantities a phase holds. One name serves both directions, so
# `state` is what `Vary(state, phase)` varies and what `@partial(f, state)` differentiates.
export state, control, parameter, segment_durations
export initial_time, final_time, initial_state, final_state
export departure_vinf, arrival_vinf, initial_mass, final_mass
export deep_space_dv, arc_fractions, forward_control, backward_control

# Reading a solved phase.
export get_initial_state, get_final_state, get_initial_time, get_final_time
export get_node_times, get_param_value, n_segments, n_dsm, subject_at

# Estimation: the problem, what is observed, how it is weighted, and how it is read in.
export ODProblem, SignalPath, AbstractMeasurement, TwoWayRange, TwoWayDoppler, MeasurementNoise
export DiagonalSNC, NoNoise
export TrackingDataFile, CCSDS_KVN, read_records, write_records, simulate
export ObservationRecord, TDMHeader, TDMSegmentMeta, build_od_closures

# Driving the filter a step at a time, for a caller that wants to look or intervene between
# updates.
export init_ekf, time_update!, measurement_update!, current_state, current_covariance

"""
    is_astrosolve_stateful(T::Type) -> Bool

Whether objects of type `T` carry state a solve iteration changes, so the
framework must restore them before the next one.

A `Spacecraft` is stateful: an event propagates it and its epoch and state move.
Run the next iteration from where the last one finished and the solver is
differentiating a moving target, which shows up as a Jacobian that will not
converge rather than as an error. A coordinate system or a maneuver's delta-V is
not stateful in this sense.

Declare a method for your own type to opt in. `false` by default, so a type that
holds no such state needs no declaration.

# Arguments
- `T::Type`: Type to test for state restored between solve iterations.

# Returns
`true` when the framework must reset objects of this type between iterations,
`false` otherwise.

# Examples
```jldoctest
using AstroModels
using AstroSolve: is_astrosolve_stateful
is_astrosolve_stateful(Spacecraft), is_astrosolve_stateful(Float64)

# output
(true, false)
```
"""
is_astrosolve_stateful(::Type) = false 
is_astrosolve_stateful(::Type{T}) where {T<:Spacecraft} = true
is_astrosolve_stateful(::Type{T}) where {T<:ImpulsiveManeuver} = true

"""
    Event

An action performed at one point in a trajectory sequence.

# Fields
- `name::String`: Name used in reports.
- `event::Function`: Action executed when the sequence reaches the event.
- `vars::Vector{Any}`: Solver variables associated with the event.
- `funcs::Vector{Any}`: Constraints and objectives evaluated at the event.

# Notes
Variables are applied before the action runs. Constraints and objectives are
evaluated afterward, against the state produced by the event.

# Examples
<!-- doc-fragment -->
```julia
# Simple propagation event
prop_event = Event(
    name = "Propagate to Apoapsis",
    event = () -> propagate!(prop, sat, StopAt(position_dot_velocity, sat; equals = 0.0, direction = -1)),
)

# Event with solver variables and constraints
sc = Spacecraft(); toi = ImpulsiveManeuver();
var_dv = SolverVariable(
    calc = ManeuverCalc(toi, sc, DeltaVVector()),
)
fun_dv = Constraint(delta_v_magnitude, toi; equals = 0.1)
maneuver_event = Event(
    event = () -> maneuver!(sc, toi),
    vars = [var_dv],
    funcs = [fun_dv],
    name = "Departure Maneuver",
)
```
"""
struct Event
    name::String
    event::Function
    vars::Vector{Any}
    funcs::Vector{Any}
end

"""
    Base.show(io::IO, event::Event)

Write the event name and scalar variable and function counts to `io`.
"""
function Base.show(io::IO, event::Event)
    name_str = isempty(event.name) ? "<unnamed>" : "\"$(event.name)\""
    
    # Count actual scalar variables (not just SolverVariable structs)
    total_vars = 0
    for var in event.vars
        if hasfield(typeof(var), :numvars)
            total_vars += var.numvars
        else
            total_vars += 1  # Fallback for unknown types
        end
    end
    
    # Count actual constraint functions (not just Constraint structs)  
    total_funcs = 0
    for func in event.funcs
        if hasfield(typeof(func), :numfuncs)
            total_funcs += func.numfuncs
        elseif hasmethod(length, (typeof(func),))
            # Try to get length if it's a collection
            try
                total_funcs += length(func)
            catch
                total_funcs += 1
            end
        else
            total_funcs += 1  # Fallback for unknown types
        end
    end
    
    print(io, "Event($name_str; $total_vars vars, $total_funcs funcs)")
end 

"""
    Event(; name::String = "", event::Function = () -> nothing, vars = [], funcs = [])

Create an event from its name, action, variables, and functions.
"""
Event(; name::String = "", event::Function = () -> nothing, vars = [], funcs = []) =
    Event(name, event, vars, funcs)

"""
    AbstractRole

What a solver is to do with a variable.

`SolveFor` is the default and the only one an optimizer understands: the value moves. The other two
belong to estimation, and are accepted on the `Vary` form rather than on a transcribed block.
`Fixed` holds the value and ignores its uncertainty. `Consider` is intended to hold the value while
letting its uncertainty propagate into the solution covariance, and is not implemented: both
estimators throw on it. See [`Consider`](@ref).
"""
abstract type AbstractRole end

"""
    SolveFor()

The role of a variable whose value the solver moves. It is the default for `Vary`, and the only
role an optimizer accepts.

See [`AbstractRole`](@ref) for the three roles together.

# Example
<!-- doc-fragment -->
```julia
Vary(state, sat; guess = y_guess, covariance = [1e2, 1e2, 1e2, 1e-2, 1e-2, 1e-2])
```
"""
struct SolveFor <: AbstractRole end

"""
    Consider()

The role of a variable an estimator would hold at its value while letting its uncertainty
propagate into the solution covariance, which is how a parameter the tracking data cannot observe
degrades the answer rather than being ignored.

This role is not implemented. `Vary` accepts it and both estimators reject it, so a `Batch` or
`Sequential` solve throws an `ArgumentError` naming the variable and its role. [`Fixed`](@ref)
holds a value, dropping its uncertainty rather than propagating it, and is the role to use until
this one is built.

See [`AbstractRole`](@ref) for the three roles together.

# Example
```julia
epoch = Time("2020-03-01T00:00:00.000", TT(), ISOT())
sat = Spacecraft(state = CartesianState([6878.137, 0.0, 0.0, 0.0, 4.71754, 5.9982]),
                 time = epoch, name = "Sat")

# Accepted here, and rejected by both estimators when the problem is solved
Vary(state, sat; guess = [6878.137, 0.0, 0.0, 0.0, 4.71754, 5.9982],
     covariance = [1e2, 1e2, 1e2, 1e-2, 1e-2, 1e-2], role = Consider())
```
"""
struct Consider <: AbstractRole end

"""
    Fixed()

The role of a variable an estimator holds at its value and whose uncertainty it ignores, so the
variable contributes nothing to the solution covariance and nothing to the estimated state.

Both estimators accept it, on the `Vary` form rather than on a transcribed block: a bounded and
scaled transcription variable takes no role, and passing one throws.

See [`AbstractRole`](@ref) for the three roles together, and [`Consider`](@ref) for the role that
would keep the uncertainty.

# Example
```julia
epoch = Time("2020-03-01T00:00:00.000", TT(), ISOT())
sat = Spacecraft(state = CartesianState([6878.137, 0.0, 0.0, 0.0, 4.71754, 5.9982]),
                 time = epoch, name = "Sat")

Vary(state, sat; guess = [6878.137, 0.0, 0.0, 0.0, 4.71754, 5.9982],
     covariance = [1e2, 1e2, 1e2, 1e-2, 1e-2, 1e-2], role = Fixed())
```
"""
struct Fixed    <: AbstractRole end

"""
    SolverVariable{C<:AbstractCalc}

Represents a solver-controlled variable defined by a calculation container for 
trajectory optimization.

# Fields
- `calc::C`: Calculation container (e.g., `OrbitCalc`, `ManeuverCalc`, `BodyCalc`)
- `numvars::Int`: Number of scalar variables for this calculation
- `lower_bound::Vector{Float64}`: Lower bounds for optimization (length `numvars`, 
  defaults to `-Inf`)
- `upper_bound::Vector{Float64}`: Upper bounds for optimization (length `numvars`, 
  defaults to `+Inf`)
- `shift::Vector{Float64}`: Variable shifting for numerical conditioning 
  (length `numvars`, defaults to `0.0`)
- `scale::Vector{Float64}`: Variable scaling for numerical conditioning 
  (length `numvars`, defaults to `1.0`)
- `name::String`: Optional human-readable identifier

# Constructor
    SolverVariable(; calc, lower_bound=nothing, upper_bound=nothing, shift=nothing,
    scale=nothing, name="")

# Arguments
- `calc`: Calculation container implementing `AbstractCalc` interface
- `lower_bound`: Scalar or vector of lower bounds (broadcast to `numvars` if scalar)
- `upper_bound`: Scalar or vector of upper bounds (broadcast to `numvars` if scalar)  
- `shift`: Scalar or vector of shift values for numerical conditioning
- `scale`: Scalar or vector of scale factors for numerical conditioning
- `name`: Optional descriptive name for the variable

# Notes
- See `AbstractCalc` for supported calculation types
- `numvars` is automatically determined from the calculation variable type
- Bounds, shift, and scale vectors are stored as `Float64` for optimizer compatibility
- The calculation must support setting values (`calc_is_settable(calc.var) == true`) 

# Examples
```julia
using Epicycle

# Spacecraft and maneuver objects
sc = Spacecraft(); toi = ImpulsiveManeuver()

# 3D delta-V vector variable with individual component bounds
var_toi = SolverVariable(
    calc = ManeuverCalc(toi, sc, DeltaVVector()),
    lower_bound = [-2.0, -2.0, -2.0],
    upper_bound = [2.0, 2.0, 2.0],
    scale = [1.0, 1.0, 1.0],  
    name = "Departure DV Vector"
)

# Orbital element variable
var_sma = SolverVariable(
    calc = OrbitCalc(sc, SMA()),
    lower_bound = 6700.0,   
    upper_bound = 42000.0,  
    name = "Target Semi-Major Axis"
)
```
"""
mutable struct SolverVariable{C}
    calc::C
    numvars::Int
    lower_bound::Vector{Float64}
    upper_bound::Vector{Float64}
    shift::Vector{Float64}
    scale::Vector{Float64}
    name::String

    # An optimizer bounds a variable in a box and scales it. An estimator gives
    # it a prior covariance and, for a filter, process noise. A variable carries
    # whichever its problem needs; `Vary` refuses the combinations that mean
    # nothing, so a box bound on an estimated state is an error rather than a
    # setting no filter reads.
    role::AbstractRole
    covariance::Any
    process_noise::Any

    function SolverVariable{C}(; 
        calc::C,
        lower_bound = nothing, 
        upper_bound = nothing, 
        shift = nothing, 
        scale = nothing, 
        name::String = "",
        role::AbstractRole = SolveFor(),
        covariance = nothing,
        process_noise = nothing
    ) where {C}
        n = calc_numvars(calc)

        to_vec(x, n, default) = x === nothing ? fill(default, n) :
                                (x isa AbstractVector ? Float64.(x) : fill(Float64(x), n))

        lb = to_vec(lower_bound, n, -Inf)
        ub = to_vec(upper_bound, n, +Inf)
        sh = to_vec(shift,       n, 0.0)
        sc = to_vec(scale,       n, 1.0)

        (length(lb) == n && length(ub) == n && length(sh) == n && length(sc) == n) ||
            throw(ArgumentError("SolverVariable: all bound/shift/scale vectors must have length $n."))

        new{C}(calc, n, lb, ub, sh, sc, name, role, covariance, process_noise)
    end
end

# The keyword constructor that infers `C` from `calc` is in direct_variable.jl,
# because a second variable kind wants the same name. Julia dispatches on
# positional arguments, and both forms take none, so only one method can exist
# and it decides from the keywords it was given. Nothing is defined here.


""" 
   Base.show(io::IO, sv::SolverVariable)

Show method for `SolverVariable` that displays key information about the variable.
"""
function Base.show(io::IO, sv::SolverVariable)
    calc_ty = nameof(typeof(sv.calc))
    var_str = hasfield(typeof(sv.calc), :var) ? string(getfield(sv.calc, :var)) : "?"
    println(io, "SolverVariable($calc_ty; var=$var_str)")
    println(io, "  numvars:     ", sv.numvars)
    println(io, "  lower_bound: ", sv.lower_bound)
    println(io, "  upper_bound: ", sv.upper_bound)
    println(io, "  shift:       ", sv.shift)
    println(io, "  scale:       ", sv.scale)
    println(io, "  name:        ", sv.name)
end

"""
    Vary(quantity, subject, deps...; guess, lower_bound, upper_bound, scale, shift, name)

Let a solver change a quantity.

Its first arguments are the same as every other spec. There is no second
vocabulary for writing: the quantity that reads a value also carries the
function that writes it, as the `setter` trait, so naming the quantity once
gets both directions.

```julia
delta_v(toi)                                     # read it
Vary(delta_v, toi; lower_bound = [0.0, 0.0, 0.0],   # let the solver change it
                   upper_bound = [8.0, 0.0, 0.0])
```

The subject is the maneuver, not the spacecraft — a delta-V belongs to the
maneuver, and applying it to a spacecraft is what `maneuver!` does.

`guess` is written immediately, so the starting value lives in the spec rather
than being implied by how the subject was built. It is a guess and not an
`initial` anything: in this domain an initial condition is a boundary condition
of the equations of motion, which is a different idea entirely.

# Returns
The `SolverVariable` it created and attached to the subject. Keeping it is only necessary where the
variable is named again later, as an `ODProblem`'s `solve_for` list does; a phase already holds its
own.

Throws an `ArgumentError` when the quantity cannot be set on that subject, since a solver with no
way to write the value has nothing to vary.
"""
function Vary(quantity::Function, subject, deps...;
              guess = nothing, lower_bound = nothing, upper_bound = nothing,
              scale = nothing, shift = nothing, name::String = "",
              role::AbstractRole = SolveFor(),
              covariance = nothing, process_noise = nothing)
    calc = Calc(quantity, subject, deps...)

    calc_is_settable(calc) || throw(ArgumentError(
        "Vary: $(EpicycleBase.label(quantity)) cannot be set on a " *
        "$(nameof(typeof(subject))) — there is no `set_quantity!` method for " *
        "that pair, so a solver has no way to change it."))

    # A box bound and a covariance are two different claims about one variable,
    # and no solver reads both. Saying so here is what keeps the mistake from
    # becoming a setting that is quietly ignored.
    if covariance !== nothing
        for (kw, v) in (("lower_bound", lower_bound), ("upper_bound", upper_bound),
                        ("scale", scale), ("shift", shift))
            v === nothing || throw(ArgumentError(
                "Vary: $kw and covariance cannot both be given. An estimated " *
                "quantity is bounded by its covariance, and there is no filter " *
                "for a box bound to act on."))
        end
    end
    process_noise === nothing || covariance !== nothing || throw(ArgumentError(
        "Vary: process_noise needs a covariance. Process noise is added to a " *
        "covariance between updates, so there is nothing for it to grow."))

    guess === nothing || set_calc!(calc, guess)

    sv = SolverVariable(; calc = calc, lower_bound = lower_bound,
                        upper_bound = upper_bound, scale = scale,
                        shift = shift, name = name, role = role,
                        covariance = covariance, process_noise = process_noise)
    return _register_vary!(sv)       # attaches to the next step inside a target! block
end

# ── Reading and writing a variable ───────────────────────────────────────────
#
# What a solver needs of a variable, whether it optimises or estimates: how many
# scalars it is, what it reads now, and how to write an iterate back. Everything
# goes through the Calc, so the subject is updated in place and the next
# propagation sees it.

"""
    length_of(v::SolverVariable)

How many scalars the variable contributes to the solver's vector.
"""
length_of(v::SolverVariable) = v.numvars

"""
    current_value(v::SolverVariable)

Read the variable's present value from its subject.
"""
current_value(v::SolverVariable) = get_calc(v.calc)

"""
    assign!(v::SolverVariable, x)

Write a solver iterate back onto the variable's subject.

# Returns
`nothing`. The subject is modified in place.
"""
function assign!(v::SolverVariable, x)
    set_calc!(v.calc, x)
    return nothing
end

"""
    set_sol_var(var::SolverVariable, val::Vector)

Set the value(s) of the solver variable struct
"""
function set_sol_var(var::SolverVariable,val::Vector)
    # Test this calc is settable
    if !calc_is_settable(var.calc)
        throw(ArgumentError("set_sol_var: calc type $(typeof(var.calc)) does not support setting values."))
    end

    # Delegate to AstroCallbacks; accept vectors for both scalar and vector calcs
    n = var.numvars
    length(val) == n || throw(ArgumentError("set_sol_var: expected length $n (got $(length(val)))."))
    if n == 1
        set_calc!(var.calc, val[1])
    else
        set_calc!(var.calc, val)
    end
    return var
end

"""
    get_sol_var(var::SolverVariable)

Get the solver variable values from the struct.
"""
function get_sol_var(var::SolverVariable)
    # Always return a Vector for SequenceManager
    vals = get_calc(var.calc)
    return vals isa AbstractVector ? vals : [vals]
end

"""
    apply_event(event::Event)

Execute the event's function closure.
"""
function apply_event(event::Event)
    event.event()
    return nothing
end

include("constraint.jl")
include("transcription.jl")
include("hermite_simpson.jl")
include("sequence.jl")
include("sequence_report.jl")

# The direct transcription layer. Order matters: each file assumes what the ones
# before it define, and the tags and the phase come before anything that
# dispatches on them.
include("element_interface.jl")
include("autodiff.jl")
include("collocation_phase.jl")
include("direct_variable.jl")
include("path_constraint.jl")
include("collocation.jl")
include("linkage.jl")
include("jacobian_storage.jl")
include("system_function.jl")
include("collocation_sequence.jl")
include("spec_api.jl")
include("diagnostics.jl")
include("collocation_eval.jl")


# Shooting. The manager and the do-block verbs come first, then the two
# transcriptions that subtype AbstractShootingPhase.
include("kepler_time_domain.jl")
include("shooting_sequence.jl")
include("spec_api_shooting.jl")
include("sims_flanagan.jl")
include("mga_ndsms.jl")

# One NLP from any mix of collocation and shooting phases, so it comes after
# every phase type it may be handed.
include("oc_manager.jl")

# The vocabulary a problem is written in. Last, because every verb in it
# dispatches on something defined above.
include("spec_vocabulary.jl")

# A targeting problem written in flight order, GMAT-style. After the vocabulary it records.
include("target_block.jl")

# Estimation. These stay submodules rather than being flattened in, because each
# is a self-contained subject — a measurement, a factorisation, a file format —
# and their names are re-exported below so a user need not say which.
include("UDU.jl")
include("ProcessNoiseModels.jl")
include("Measurements.jl")
include("TrackingDataIO.jl")

# The estimators. Each reaches its sibling submodules above through the parent module,
# with `..`.
include("SpringMassEstimator.jl")
include("BatchLeastSquares.jl")
include("ExtendedKalmanFilter.jl")

using .UDU: udu_from_P, udu_to_P, thornton_time_update!, bierman_measurement_update!
using .ProcessNoiseModels: ProcessNoiseModel, NoNoise, DiagonalSNC, discretize
using .Measurements: SignalPath, participant_names, AbstractMeasurement,
                     TwoWayRange, TwoWayDoppler, AbstractMeasurementNoise,
                     MeasurementNoise, draw, variance, is_visible
using .BatchLeastSquares: ODProblem, solve_batch_ls!, BatchLSResult,
                          build_od_closures
using .ExtendedKalmanFilter: run_ekf!, run_rts, run_iterated_rts!,
                             init_ekf, time_update!, measurement_update!,
                             current_state, current_covariance
using .SpringMassEstimator: solve_spring_mass_batch!

using Random: Random
using .TrackingDataIO: ObservationRecord, TrackingDataFile,
                       AbstractTrackingDataFormat, CCSDS_KVN, TDMHeader,
                       TDMSegmentMeta, write_records, read_records, TDM_VERSION,
                       ALLOWED_TIME_SYSTEMS, ALLOWED_MODES, ALLOWED_RANGE_UNITS,
                       SUPPORTED_OBSERVABLES, REQUIRED_META_KEYS

# Simulated tracking data. It uses the measurement predictors, the noise draw and the record type
# above, so it is included after them.
include("MeasurementSimulation.jl")

include("precompile.jl")

end
