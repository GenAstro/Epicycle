# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# CollocationPhase, BoundaryFunction and the two objective forms.
#

# ─────────────────────────────────────────────────────────────────────────────
# CollocationPhase
# ─────────────────────────────────────────────────────────────────────────────

"""
    CollocationPhase(; name, transcription, dynamics = nothing, model = nothing,
                       state = nothing, control = nothing, tspan = nothing)

One arc of an optimal control problem, transcribed by direct collocation into a nonlinear
program. The phase fixes the dynamics, the state and control dimensions and the time interval;
`Vary`, `Constraint` and `Objective` then say what the solver may change, what must hold and what
it minimizes, and `solve!(Sequence(phase))` solves it.

# Arguments
- `name::Symbol`: the phase's name, used in messages and when the phase is displayed.
- `transcription`: the collocation scheme and its mesh, such as `HermiteSimpson(n_steps = 50)`.
- `dynamics`: the right-hand side `f!(dy, y, u, p, t, model)`, which writes the state derivative
  into `dy`. `y` and `u` arrive as instances of the `state` and `control` types.
- `model`: any value the dynamics needs, such as a mass parameter, passed to it as `model`.
- `state`, `control`: struct types whose fields are the state and control components, subtypes of
  `AbstractState` and `AbstractControl`. Their field counts set the phase's dimensions.
- `tspan`: `(t0, tf)` in the time units of the dynamics. Both ends are held fixed until
  `Vary(initial_time, phase; ...)` or `Vary(final_time, phase; ...)` frees them.

# Notes
Throws `ArgumentError` when `tspan` does not increase or `state` has no fields. A phase without
a `state` type can be built, but varying its state or control then raises an `ArgumentError`
naming the missing keyword.

# Returns
The phase. After `solve!`, `state(phase)` and `control(phase)` return the solution at the mesh
nodes, one column per node, and `get_node_times(phase)` the node times.

# Example
```julia
using EpicycleBase, AstroSolve

struct Pos{T} <: AbstractState;   x::T end
struct Rate{T} <: AbstractControl; u::T end
integrator!(dy, y, u, p, t, model) = (dy[1] = u.u)

# Drive x from 1.5 toward zero with least control effort plus a terminal penalty
phase = CollocationPhase(name = :hull, transcription = HermiteSimpson(n_steps = 20),
                         dynamics = integrator!, state = Pos, control = Rate,
                         tspan = (0.0, 1.0))
Vary(state, phase;   guess = [1.5 0.25], lower_bound = [-10.0], upper_bound = [10.0])
Vary(control, phase; guess = [-1.25 -1.25], lower_bound = [-10.0], upper_bound = [10.0])
Constraint(c -> [state(c).x], phase; equals = [1.5], at = Initial())
Objective(c -> 2.5 * state(c).x^2, phase; sense = Min())
Objective(c -> 0.5 * control(c).u^2, phase; sense = Min(), at = Path())

solve!(Sequence(phase); method = Optimize(print_level = 0))
state(phase)[1, end]    # final x, about 0.25
```
"""
mutable struct CollocationPhase
    name::Symbol
    transcription::Any   # LGL, HermiteSimpson, or future transcription type
    dynamics::Any
    model::Any
    state_var::Any
    control_var::Any
    param_var::Any               # DirectSolverVariable with AbstractParameter var, or nothing
    t0_var::Any
    tf_var::Any
    constraints::Vector{Any}
    path_constraints::Vector{Any}
    objective::Any
    _y0::Vector{Float64}   # populated by framework before boundary closures run
    _yf::Vector{Float64}
    _t0::Float64
    _tf::Float64
    _n_states::Int              # set by set_state!
    _n_controls::Int            # set by set_control!
    _n_params::Int              # set by set_parameter!
    _Y::Matrix{Float64}         # n_states  × N, filled by set_decision_vector!
    _U::Matrix{Float64}         # n_controls × N, filled by set_decision_vector!
    _params::Vector{Float64}    # length _n_params, filled by set_decision_vector!
    dynamics_jac::Dict{UInt64, Function}   # keyed by objectid(var); one fn per variable block
    # A phase is complete when it exists: it knows its dimensions, its physics
    # and its interval. `state` and `control` are types, so the counts derive
    # rather than being stated twice, and the invariants can be checked here
    # (CodingStandards 7.8). Bounds and guesses are not phase definition —
    # they say what a solver may change, and arrive through `Vary`.
    function CollocationPhase(; name, transcription,
                                dynamics = nothing, model = nothing,
                                state    = nothing, control = nothing,
                                tspan    = nothing)
        ns = state   === nothing ? 0 : n_components(state)
        nc = control === nothing ? 0 : n_components(control)
        state === nothing || ns > 0 ||
            throw(ArgumentError("CollocationPhase: state type $(state) has no components."))
        tspan === nothing || tspan[2] > tspan[1] ||
            throw(ArgumentError("CollocationPhase: tspan must increase, got $(tspan)."))
        t0, tf = tspan === nothing ? (0.0, 0.0) : (float(tspan[1]), float(tspan[2]))

        p = new(name, transcription,
                nothing, model, nothing, nothing, nothing, nothing, nothing,
                Any[], Any[], nothing,
                Float64[], Float64[], t0, tf,
                ns, nc, 0,
                Matrix{Float64}(undef,0,0), Matrix{Float64}(undef,0,0), Float64[],
                Dict{UInt64, Function}())

        if state !== nothing
            _phase_registry[objectid(p)] = _PhaseRegistry(state, control, nothing, ns, nc)
        end
        dynamics === nothing || set_dynamics!(p, dynamics; model = model)
        if tspan !== nothing
            set_initial_time!(p; equality = t0)
            set_final_time!(p;   equality = tf)
        end
        return p
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# BoundaryFunction
# ─────────────────────────────────────────────────────────────────────────────

struct BoundaryFunction
    phases  ::Tuple                    # sparsity declarations (which NLP variable blocks this touches)
    fn      ::Function                 # (ctx::BoundaryContext) -> Vector; framework passes current values
    name    ::String                   # display label, e.g. "start_bc"
    jac_fns ::Dict{UInt64, Function}   # analytic Jacobian closures, keyed by objectid(var)
end

# do-block sugar: BoundaryFunction(phases...; name="") do ... end
BoundaryFunction(fn::Function, phases...; name::String = "") =
    BoundaryFunction(phases, fn, name, Dict{UInt64, Function}())

# ─────────────────────────────────────────────────────────────────────────────
# MayerObjective
# ─────────────────────────────────────────────────────────────────────────────

struct MayerObjective
    phases  ::Tuple
    sense   ::Symbol
    fn      ::Function                 # (ctx::BoundaryContext) -> scalar; framework passes current values
    jac_fns ::Dict{UInt64, Function}   # keyed by objectid(var), same pattern as BoundaryFunction
end

# do-block sugar: MayerObjective(phase; sense=:Min) do ... end
MayerObjective(fn::Function, phases...; sense::Symbol = :Min) =
    MayerObjective(phases, sense, fn, Dict{UInt64, Function}())

# ─────────────────────────────────────────────────────────────────────────────
# BolzaObjective
#
# Bolza-form objective:  J = φ(y_f) + ∫[t0,tf] L(y,u) dt
#
#   mayer_fn()            → scalar   (terminal cost; use () -> 0.0 for pure Lagrange)
#   lagrange_fn(y_k, u_k) → scalar   (integrand evaluated at one node)
#
# The Lagrange integral is approximated via LGL quadrature (lgl_quadrature_weights).
# Gradient closures are registered with the same add_objective_jacobian! API as
# MayerObjective. When the phase's Mayer term came from add_mayer!, it is a
# MayerObjective with its own partials, and this table holds the Lagrange half's
# alone; objective_gradient_chunk sums the two.
# ─────────────────────────────────────────────────────────────────────────────

struct BolzaObjective
    phases      ::Tuple
    sense       ::Symbol
    mayer_fn    ::Function          # () -> scalar
    lagrange_fn ::Function          # (y_k::Vector, u_k::Vector) -> scalar
    jac_fns     ::Dict{UInt64, Function}
end

BolzaObjective(mayer::Function, lagrange::Function, phases...;
               sense::Symbol = :Min) =
    BolzaObjective(phases, sense, mayer, lagrange, Dict{UInt64, Function}())
