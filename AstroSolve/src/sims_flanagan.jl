# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0


# SimsFlanagan.jl
#
# Sims-Flanagan low-thrust transcription for use with the Epicycle optimal
# control framework.
#
# Reference: Ellison, Englander, Ozimek, Conway — AAS 14-310
#   "Analytical Partial Derivative Calculation of the Sims-Flanagan
#    Transcription Match Point Constraints"
#
# Phase structure (Figure 1 of the paper):
#
#   [left ctrl pt] ──[ fwd N/2 segments ]──► [◇ match pt] ◄──[ bwd N/2 segments ]── [right ctrl pt]
#
# Control points (blue circles in Fig. 1) are planets/bodies at phase boundaries.
# Their position comes from an ephemeris callback, NOT from the decision vector.
# The spacecraft state at each control point is:
#
#   x_sc(t0) = [ r_planet(t0);  v_planet(t0) + v∞_dep ]
#   x_sc(tf) = [ r_planet(tf);  v_planet(tf) + v∞_arr ]
#
# Decision variables per phase: v∞_dep (3), v∞_arr (3), u_fwd (3×N/2),
#                                 u_bwd (3×N/2), t0 (1), tf (1), m0 (1), mf (1)
#
# The match-point defect (7 constraints) is internal to the phase.
# Linkage between phases connects match points of adjacent phases.
# Flyby constraints (v∞ magnitude equality, safe altitude) are
# add_constraint! calls on the linkage — no special-casing needed.

# ─────────────────────────────────────────────────────────────────────────────
# AbstractShootingPhase
# Defined in OptControlStubs.jl so ShootingManager / add_sequence! can
# dispatch on it before SimsFlanagan.jl is included.
# ─────────────────────────────────────────────────────────────────────────────

# (abstract type AbstractShootingPhase end  ← declared in OptControlStubs.jl)

# ─────────────────────────────────────────────────────────────────────────────
# SimsFlanagan  —  transcription descriptor (parallel role to LGL in
#                  collocation)
# ─────────────────────────────────────────────────────────────────────────────

"""
    SimsFlanagan(; n_segments = 10, duty_cycle = 1.0, n_thrusters = 1,
                   throttle_smoothing = 1.0e-10)

Transcription descriptor for a Sims-Flanagan phase, which divides the phase into `n_segments`
equal segments and models the thrust on each one as an impulse applied at the segment midpoint.

# Fields
- `n_segments::Int`: total number of segments, which must be even. Half propagate forward from the
  left control point and half backward from the right, meeting at the match point.
- `duty_cycle::Float64`: fraction of each segment over which the thruster fires, in `(0, 1]`.
- `n_thrusters::Int`: number of thrusters contributing thrust.
- `throttle_smoothing::Float64`: throttle magnitude below which a segment reads as coasting.

# Notes
The propellant a segment burns is proportional to the magnitude of its throttle vector, and that
magnitude is not differentiable where the throttle is zero. `throttle_smoothing` replaces it with
`sqrt(u·u + throttle_smoothing^2)`, which is smooth through the origin, so a coasting segment has
a well-defined derivative. The mass update and the match-point Jacobian use the same value, so the
propellant the constraint burns is the propellant an objective prices.

Smaller values track the true propellant more closely and make the problem harder to solve, and the
default is small enough to be no regularisation at all, which is why a mass-maximizing solve usually
needs it raised. On an Earth to Apophis rendezvous with 40 segments, 0.015 converges at the solver's
default tolerance and 0.01 does not, and the two differ by 1.7 kg of propellant out of 287.

The cost of raising it is that a coasting segment is charged `throttle_smoothing` of throttle rather
than none, so a trajectory that coasts pays propellant it does not burn: eight segments at 0.015 lose
1.9e-3 of a unit initial mass. That is why the default is not raised for everyone, and why a solve
that is meant to coast wants it left alone.

# Example
```julia
using AstroSolve
SimsFlanagan(n_segments = 40, throttle_smoothing = 0.015)
```
"""
struct SimsFlanagan
    n_segments        ::Int        # total number of segments; must be even (N/2 fwd, N/2 bwd)
    duty_cycle        ::Float64    # thruster duty cycle D ∈ (0, 1]
    n_thrusters       ::Int        # number of thrusters contributing thrust
    throttle_smoothing::Float64    # throttle magnitude below which a segment reads as coasting
end

function SimsFlanagan(; n_segments::Int = 10,
                        duty_cycle::Float64 = 1.0,
                        n_thrusters::Int = 1,
                        throttle_smoothing::Float64 = 1.0e-10)
    iseven(n_segments) || throw(ArgumentError(
        "n_segments must be even (got $n_segments); each half-phase needs N/2 segments"))
    throttle_smoothing > 0 || throw(ArgumentError(
        "throttle_smoothing must be positive (got $throttle_smoothing); it is the throttle " *
        "magnitude below which a segment reads as coasting"))
    SimsFlanagan(n_segments, duty_cycle, n_thrusters, throttle_smoothing)
end

n_fwd(t::SimsFlanagan) = t.n_segments ÷ 2
n_bwd(t::SimsFlanagan) = t.n_segments ÷ 2

# ─────────────────────────────────────────────────────────────────────────────
# SimsFlanaganPhase <: AbstractShootingPhase
# ─────────────────────────────────────────────────────────────────────────────

"""
    SimsFlanaganPhase(; name, transcription, model, ephemeris_left, ephemeris_right,
                        tspan, matchpoint_scale = ones(7))

A low-thrust arc between two control points, transcribed by [`SimsFlanagan`](@ref). Each end follows
a body through an ephemeris callback, the segments propagate from both ends to a match point in the
middle, and the seven-component disagreement there is the phase's built-in constraint.

# Arguments
- `name::Symbol`: the phase's name, which labels it in reports.
- `transcription::SimsFlanagan`: the segment count and the throttle smoothing.
- `model`: a `PropulsionModel` carrying `mu`, `Isp`, `Tmax` and `g0`. The four can be passed
  individually instead, and `mu`, `Isp` and `Tmax` have no default, so one of the two forms is
  required.
- `ephemeris_left`, `ephemeris_right`: callables `eph(t) -> (r, v, a)`, each a 3-vector, giving the
  body's position, velocity and gravitational acceleration at the control point. The acceleration is
  what the time partial of the match point needs.
- `tspan::Tuple`: the phase's start and end times, which must increase.
- `matchpoint_scale::Vector{Float64}`: seven divisors, one per match-point residual, in the order
  three position, three velocity, one mass. They make the two halves' disagreement comparable across
  quantities with different sizes.

# Fields
The struct's fields are the values above, with `model` decomposed into `mu`, `Isp`, `Tmax` and `g0`,
alongside the solver variables and the cached decision values the framework fills during a solve. A
caller sets none of them directly; `Vary`, `Constraint` and `Objective` do.

# Notes
Units are the caller's throughout and have only to agree with each other: `mu` sets the length and
time units, and `Tmax`, `g0` and the ephemeris follow them.

A Sims-Flanagan phase varies throttle blocks and masses rather than a state and a control history,
so its quantities have their own names: `forward_control`, `backward_control`, `initial_mass` and
`final_mass`, along with `initial_time` and `final_time`.

# Example
```julia
using AstroSolve
phase = SimsFlanaganPhase(name = :transfer,
                          transcription = SimsFlanagan(n_segments = 40),
                          model = PropulsionModel(mu = 1.32712440018e11,
                                                  Isp = 3000.0,
                                                  Tmax = 1.0e-3,
                                                  g0 = 9.80665e-3),
                          ephemeris_left = t -> ([1.5e8, 0.0, 0.0], [0.0, 29.8, 0.0], zeros(3)),
                          ephemeris_right = t -> ([2.3e8, 0.0, 0.0], [0.0, 24.1, 0.0], zeros(3)),
                          tspan = (0.0, 3.1e7),
                          matchpoint_scale = [1.5e8, 1.5e8, 1.5e8, 29.8, 29.8, 29.8, 1500.0])
```
"""
mutable struct SimsFlanaganPhase <: AbstractShootingPhase
    name            ::Symbol
    transcription   ::SimsFlanagan

    # ── Spacecraft / propulsion constants ──────────────────────────────────
    mu              ::Float64   # gravitational parameter (km³/s² or m³/s²)
    Isp             ::Float64   # specific impulse (s)
    Tmax            ::Float64   # max thrust per thruster (N, consistent units)
    g0              ::Float64   # reference gravity (m/s², consistent units)

    # ── Ephemeris callbacks ────────────────────────────────────────────────
    # Signature: eph(t::Float64) -> (r::Vector{Float64}(3),
    #                                v::Vector{Float64}(3),
    #                                a::Vector{Float64}(3))
    # r, v = planet position/velocity; a = acceleration (needed for ∂x_planet/∂t
    # in the time-partial of the match-point defect).
    ephemeris_left  ::Any       # ephemeris at left  control point (t0)
    ephemeris_right ::Any       # ephemeris at right control point (tf)

    # ── Decision variable blocks ───────────────────────────────────────────
    x0_var          ::Any       # SolverVariable: 6-vector free state at left  ctrl pt
    xf_var          ::Any       # SolverVariable: 6-vector free state at right ctrl pt
    vinf_dep_var    ::Any       # SolverVariable: 3-vector v∞ at left  control pt
    vinf_arr_var    ::Any       # SolverVariable: 3-vector v∞ at right control pt
    u_fwd_var       ::Any       # SolverVariable: 3 × (N/2) throttle matrix, fwd
    u_bwd_var       ::Any       # SolverVariable: 3 × (N/2) throttle matrix, bwd
    t0_var          ::Any       # SolverVariable: phase start time
    tf_var          ::Any       # SolverVariable: phase end   time
    m0_var          ::Any       # SolverVariable: mass at left  control point
    mf_var          ::Any       # SolverVariable: mass at right control point

    # ── User-supplied boundary constraints and objective ───────────────────
    constraints         ::Vector{Any}
    objective           ::Any

    # ── Match-point constraint scaling ────────────────────────────────────
    # g_nlp = (g_phys - matchpoint_shift) ./ matchpoint_scale  (element-wise, 7-vector)
    matchpoint_scale    ::Vector{Float64}   # default: ones(7)
    matchpoint_shift    ::Vector{Float64}   # default: zeros(7)

    # ── Cached decision variable values (set_decision_vector! fills these) ─
    _vinf_dep       ::Vector{Float64}   # [3]  departure v∞
    _vinf_arr       ::Vector{Float64}   # [3]  arrival   v∞
    _u_fwd          ::Matrix{Float64}   # [3, N/2]  forward  throttle vectors
    _u_bwd          ::Matrix{Float64}   # [3, N/2]  backward throttle vectors
    _t0             ::Float64
    _tf             ::Float64
    _m0             ::Float64
    _mf             ::Float64

    # ── Derived control point states (built from ephemeris + v∞) ──────────
    # These are constructed in evaluate_matchpoint! and cached for use in
    # matchpoint_jacobian so the ephemeris is only called once per evaluation.
    _x0             ::Vector{Float64}   # [6]  left  ctrl pt spacecraft state
    _xf             ::Vector{Float64}   # [6]  right ctrl pt spacecraft state
    _planet_left    ::NTuple{3, Vector{Float64}}  # (r, v, a) at t0
    _planet_right   ::NTuple{3, Vector{Float64}}  # (r, v, a) at tf

    # ── Match-point propagated states (filled by evaluate_matchpoint!) ─────
    _x_match_fwd    ::Vector{Float64}   # [6]  fwd propagation end at match pt
    _x_match_bwd    ::Vector{Float64}   # [6]  bwd propagation end at match pt
    _m_match_fwd    ::Float64
    _m_match_bwd    ::Float64

    # ── Segment boundary history (warm-starting / diagnostics) ────────────
    _states_fwd     ::Matrix{Float64}   # [6, N/2+1]
    _states_bwd     ::Matrix{Float64}   # [6, N/2+1]
    _masses_fwd     ::Vector{Float64}   # [N/2+1]
    _masses_bwd     ::Vector{Float64}   # [N/2+1]

    # mu, Isp, Tmax and g0 are a model. They live behind `model = …` like every
    # other physical constant, which is where collocation already puts them.
    # The explicit keywords still work so existing scripts do not move.
    function SimsFlanaganPhase(; name, transcription, model = nothing,
                                 mu = nothing, Isp = nothing, Tmax = nothing,
                                 g0 = nothing,
                                 ephemeris_left = nothing, ephemeris_right = nothing,
                                 tspan = nothing, matchpoint_scale = nothing)
        mu   = _model_field(model, :mu,   mu,   nothing)
        Isp  = _model_field(model, :Isp,  Isp,  nothing)
        Tmax = _model_field(model, :Tmax, Tmax, nothing)
        g0   = _model_field(model, :g0,   g0,   9.80665)
        (mu === nothing || Isp === nothing || Tmax === nothing) && throw(ArgumentError(
            "SimsFlanaganPhase needs mu, Isp and Tmax — give `model = PropulsionModel(...)`."))
        tspan === nothing || tspan[2] > tspan[1] ||
            throw(ArgumentError("SimsFlanaganPhase: tspan must increase, got $(tspan)."))
        matchpoint_scale === nothing || length(matchpoint_scale) == 7 || throw(ArgumentError(
            "SimsFlanaganPhase: matchpoint_scale has one entry per match-point residual, " *
            "position, velocity and mass, so 7; got $(length(matchpoint_scale))."))
        t0, tf = tspan === nothing ? (0.0, 0.0) : (float(tspan[1]), float(tspan[2]))
        N2    = n_fwd(transcription)
        zero3 = zeros(3)
        new(
            name, transcription, mu, Isp, Tmax, g0,
            ephemeris_left, ephemeris_right,
            nothing, nothing,          # x0_var, xf_var
            # decision variable slots (vinf_dep, vinf_arr, u_fwd, u_bwd, t0, tf, m0, mf)
            nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing,
            # constraints / objective
            Any[], nothing,
            # matchpoint scaling (default: no scaling)
            matchpoint_scale === nothing ? ones(7) : float.(collect(matchpoint_scale)), zeros(7),
            # cached decision values
            zeros(3), zeros(3),
            zeros(3, N2), zeros(3, N2),
            t0, tf, 0.0, 0.0,
            # derived control point states
            zeros(6), zeros(6),
            (copy(zero3), copy(zero3), copy(zero3)),
            (copy(zero3), copy(zero3), copy(zero3)),
            # match-point states
            zeros(6), zeros(6), 0.0, 0.0,
            # segment boundary history
            zeros(6, N2 + 1), zeros(6, N2 + 1),
            zeros(N2 + 1),    zeros(N2 + 1),
        )
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Accessors
# ─────────────────────────────────────────────────────────────────────────────

"""
    n_segments(phase::SimsFlanaganPhase) -> Int

The total number of segments the phase is divided into, half propagated forward from the left
control point and half backward from the right.

# Returns
The segment count as an `Int`, always even.

# Example
<!-- doc-fragment -->
```julia
n_segments(phase)
```
"""
n_segments(p::SimsFlanaganPhase) = p.transcription.n_segments
n_fwd(p::SimsFlanaganPhase)      = n_fwd(p.transcription)
n_bwd(p::SimsFlanaganPhase)      = n_bwd(p.transcription)

# ─────────────────────────────────────────────────────────────────────────────
# set_*! configuration API  (mirrors CollocationPhase pattern)
# ─────────────────────────────────────────────────────────────────────────────

set_initial_state!(p::SimsFlanaganPhase, var)   = (p.x0_var = var)
set_final_state!(p::SimsFlanaganPhase,   var)   = (p.xf_var = var)
set_departure_vinf!(p::SimsFlanaganPhase, var)  = (p.vinf_dep_var = var)
set_arrival_vinf!(p::SimsFlanaganPhase,   var)  = (p.vinf_arr_var = var)
set_ephemeris_left!(p::SimsFlanaganPhase,  eph)  = (p.ephemeris_left  = eph)
set_ephemeris_right!(p::SimsFlanaganPhase, eph)  = (p.ephemeris_right = eph)
set_fwd_control!(p::SimsFlanaganPhase,   var)  = (p.u_fwd_var = var)
set_bwd_control!(p::SimsFlanaganPhase,   var)  = (p.u_bwd_var = var)
set_initial_time!(p::SimsFlanaganPhase,  var)  = (p.t0_var  = var)
set_final_time!(p::SimsFlanaganPhase,    var)  = (p.tf_var  = var)
set_initial_mass!(p::SimsFlanaganPhase,  var)  = (p.m0_var  = var)
set_final_mass!(p::SimsFlanaganPhase,    var)  = (p.mf_var  = var)
set_objective!(p::SimsFlanaganPhase, obj::MayerObjective) = (p.objective = obj)
set_matchpoint_scale!(p::SimsFlanaganPhase, sc) = (p.matchpoint_scale = Float64.(sc); nothing)
set_matchpoint_shift!(p::SimsFlanaganPhase, sh) = (p.matchpoint_shift = Float64.(sh); nothing)

function add_constraint!(p::SimsFlanaganPhase, con)
    push!(p.constraints, con)
    return p
end

# ─────────────────────────────────────────────────────────────────────────────
# evaluate_matchpoint!
#
# Propagates both half-phases from their respective control points to the
# internal match point and caches the results in the phase.  Called by the
# framework's get_functions dispatch before assembling the NLP residual.
#
# Requires: kepler_propagate_time_domain  (kepler_time_domain.jl)
# ─────────────────────────────────────────────────────────────────────────────

function evaluate_matchpoint!(p::SimsFlanaganPhase)
    N2       = n_fwd(p)
    dt       = (p._tf - p._t0) / p.transcription.n_segments
    mdot_max = p.Tmax / (p.Isp * p.g0)
    D        = p.transcription.duty_cycle
    nT       = p.transcription.n_thrusters
    eps2     = p.transcription.throttle_smoothing^2

    # ── Build control-point spacecraft states from ephemeris + v∞ ─────────
    # The planet position is fixed by the ephemeris; only v∞ is a decision var.
    # Cache planet (r, v, a) so matchpoint_jacobian can reuse without re-calling.
    if p.ephemeris_left !== nothing
        r0, v0, a0 = p.ephemeris_left(p._t0)
        p._planet_left = (r0, v0, a0)
        p._x0 = vcat(r0, v0 .+ p._vinf_dep)
    end
    if p.ephemeris_right !== nothing
        rf, vf, af = p.ephemeris_right(p._tf)
        p._planet_right = (rf, vf, af)
        p._xf = vcat(rf, vf .+ p._vinf_arr)
    end

    Fwd = _propagate_half_with_stm(p._x0, p._m0, p._u_fwd, dt, p.mu,
                                    p.Tmax, mdot_max, D, nT, -1.0, +1.0, eps2)
    Bwd = _propagate_half_with_stm(p._xf, p._mf, p._u_bwd, dt, p.mu,
                                    p.Tmax, mdot_max, D, nT, +1.0, -1.0, eps2)

    p._states_fwd   = Fwd.x_seg
    p._masses_fwd   = Fwd.m
    p._x_match_fwd  = Fwd.x_seg[:, N2+1]
    p._m_match_fwd  = Fwd.m[N2+1]

    p._states_bwd   = Bwd.x_seg
    p._masses_bwd   = Bwd.m
    p._x_match_bwd  = Bwd.x_seg[:, N2+1]
    p._m_match_bwd  = Bwd.m[N2+1]

    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# matchpoint_defect
#
# Returns the 7-vector match-point constraint (Eq. 6 from the paper):
#   c_mp = [X_bwd - X_fwd; m_bwd - m_fwd]  ∈ ℝ⁷  (= 0 at feasibility)
# ─────────────────────────────────────────────────────────────────────────────

function matchpoint_defect(p::SimsFlanaganPhase)
    vcat(p._x_match_bwd .- p._x_match_fwd,
         p._m_match_bwd  - p._m_match_fwd)
end

# ─────────────────────────────────────────────────────────────────────────────
# _propagate_half_with_stm
#
# Propagates one half-phase (forward or backward) storing STMs and ancillary
# quantities needed by matchpoint_jacobian.
#
# Arguments:
#   x_cp    : 6-vector, control point state (start of this half-phase)
#   m_cp    : scalar, mass at control point
#   u       : 3 × N2 throttle matrix
#   dt      : scalar, time step per full segment (half-phase uses dt/2 each side)
#   mu      : gravitational parameter
#   Tmax    : max thrust
#   mdot_max: Tmax / (Isp * g0)
#   D       : duty cycle
#   nT      : number of thrusters
#   sign_m  : +1 for forward (mass depletes), -1 for backward (mass accumulates)
#   sign_dt : +1 for forward (propagate +dt/2), -1 for backward (propagate -dt/2)
#
# Returns a NamedTuple with:
#   stm1  :: Vector{Matrix{Float64}}  — STM for first  half-prop of each segment [N2]
#   stm2  :: Vector{Matrix{Float64}}  — STM for second half-prop of each segment [N2]
#   dxdt1 :: Vector{Vector{Float64}}  — ∂x/∂Δt for first  half-prop [N2] (6-vec)
#   dxdt2 :: Vector{Vector{Float64}}  — ∂x/∂Δt for second half-prop [N2] (6-vec)
#   m     :: Vector{Float64}          — mass at start of each segment [N2+1]
#   dvmax :: Vector{Float64}          — Δv_max for each segment [N2]
#   x_seg :: Matrix{Float64}          — 6×(N2+1) states at segment boundaries
#   x_mid :: Matrix{Float64}          — 6×N2 post-impulse midpoint states
# ─────────────────────────────────────────────────────────────────────────────

function _propagate_half_with_stm(x_cp, m_cp, u, dt, mu, Tmax, mdot_max, D, nT,
                                   sign_m, sign_dt, eps2)
    N2   = size(u, 2)
    dt2  = dt / 2.0

    # Promote so ForwardDiff Duals flow through when differentiating wrt x_cp/m_cp/u/dt
    T = promote_type(eltype(x_cp), typeof(m_cp), eltype(u), typeof(dt))

    stm1  = [zeros(T, 6, 6) for _ in 1:N2]
    stm2  = [zeros(T, 6, 6) for _ in 1:N2]
    dxdt1 = [zeros(T, 6)    for _ in 1:N2]
    dxdt2 = [zeros(T, 6)    for _ in 1:N2]
    m     = zeros(T, N2 + 1)
    dvmax = zeros(T, N2)
    x_seg = zeros(T, 6, N2 + 1)
    x_mid = zeros(T, 6, N2)    # post-impulse midpoint states

    x        = convert(Vector{T}, x_cp)
    m[1]     = m_cp
    x_seg[:, 1] = x

    for k in 1:N2
        u_k     = u[:, k]
        m_k     = m[k]
        dv_mx   = Tmax * D * nT * dt / m_k
        dvmax[k] = dv_mx
        dv      = dv_mx * u_k

        # first half-propagation
        r1 = kepler_propagate_time_domain(x, mu, sign_dt * dt2;
                 need_stm = true, need_time_partials = true)
        stm1[k]  = r1.stm
        dxdt1[k] = r1.dstate_dt
        x        = r1.state

        # apply impulse at midpoint
        x[4:6] .+= dv
        x_mid[:, k] = x

        # second half-propagation
        r2 = kepler_propagate_time_domain(x, mu, sign_dt * dt2;
                 need_stm = true, need_time_partials = true)
        stm2[k]  = r2.stm
        dxdt2[k] = r2.dstate_dt
        x        = r2.state

        # mass update — hypot-style regularisation: sqrt(||u||² + ε²)
        # ForwardDiff differentiates this as u/sqrt(||u||²+ε²) = u/norm_uk,
        # exactly matching the analytic Jacobian formula u_k/nu_k. eps2 is the square of the
        # transcription's throttle_smoothing, so the mass chain and the Jacobian share one number.
        norm_uk   = sqrt(dot(u_k, u_k) + eps2)
        m[k+1]   = m_k + sign_m * mdot_max * D * norm_uk * dt
        x_seg[:, k+1] = x
    end

    return (stm1 = stm1, stm2 = stm2, dxdt1 = dxdt1, dxdt2 = dxdt2,
            m = m, dvmax = dvmax, x_seg = x_seg, x_mid = x_mid)
end

# ─────────────────────────────────────────────────────────────────────────────
# matchpoint_jacobian
#
# Computes the analytic Jacobian of the 7-constraint match-point defect
#   c_mp = [X_bwd - X_fwd; m_bwd - m_fwd]  ∈ ℝ⁷
# with respect to all decision variable blocks.
#
# Reference: Ellison et al. AAS 14-310, Equations (11)–(33).
#
# Time partial derivation (ephemeris formulation):
#   x0 = [r_planet(t0); v_planet(t0) + v∞_dep]
#   ∂x0/∂t0 = [v_planet(t0); a_planet(t0)]
#   ∂c_mp/∂t0|_full = ∂c_mp/∂t0|_Kepler + J_x0 · [v_planet(t0); a_planet(t0)]
#   (and similarly for tf / xf)
#
# Returns a NamedTuple of Jacobian blocks (all as dense matrices / vectors):
#   dvinf_dep :: Matrix{Float64}  7 × 3        ∂c_mp / ∂v∞_dep   (= J_x0[:,4:6])
#   dvinf_arr :: Matrix{Float64}  7 × 3        ∂c_mp / ∂v∞_arr   (= J_xf[:,4:6])
#   dx0       :: Matrix{Float64}  7 × 6        ∂c_mp / ∂x0  (internal; useful for debugging)
#   dxf       :: Matrix{Float64}  7 × 6        ∂c_mp / ∂xf  (internal; useful for debugging)
#   du_fwd    :: Matrix{Float64}  7 × 3*(N/2)  ∂c_mp / ∂vec(u_fwd)
#   du_bwd    :: Matrix{Float64}  7 × 3*(N/2)  ∂c_mp / ∂vec(u_bwd)
#   dm0       :: Vector{Float64}  7            ∂c_mp / ∂m0
#   dmf       :: Vector{Float64}  7            ∂c_mp / ∂mf
#   dt0       :: Vector{Float64}  7            ∂c_mp / ∂t0  (includes ephemeris term)
#   dtf       :: Vector{Float64}  7            ∂c_mp / ∂tf  (includes ephemeris term)
#
# Columns of du_fwd/du_bwd are ordered: [u1x, u1y, u1z, u2x, u2y, u2z, ...]
# ─────────────────────────────────────────────────────────────────────────────

function matchpoint_jacobian(p::SimsFlanaganPhase)
    N2       = n_fwd(p)
    dt       = (p._tf - p._t0) / p.transcription.n_segments
    mdot_max = p.Tmax / (p.Isp * p.g0)
    D        = p.transcription.duty_cycle
    nT       = p.transcription.n_thrusters
    ε        = p.transcription.throttle_smoothing^2   # sqrt(||u||²+ε) form, ε = smoothing²
    eps2     = ε

    # ── Propagate both half-phases, collecting STMs ──────────────────────────
    Fwd = _propagate_half_with_stm(p._x0, p._m0, p._u_fwd, dt, p.mu,
                                    p.Tmax, mdot_max, D, nT, -1.0, +1.0, eps2)
    Bwd = _propagate_half_with_stm(p._xf, p._mf, p._u_bwd, dt, p.mu,
                                    p.Tmax, mdot_max, D, nT, +1.0, -1.0, eps2)

    # ── Segment STMs: Φ_seg_k = stm2[k] · stm1[k] ──────────────────────────
    seg_stm_fwd = [Fwd.stm2[k] * Fwd.stm1[k] for k in 1:N2]
    seg_stm_bwd = [Bwd.stm2[k] * Bwd.stm1[k] for k in 1:N2]

    # ── Tail STMs: tail[k] = Φ_{N2} · ... · Φ_{k+1}, tail[N2] = I ──────────
    # tail_fwd[k] propagates a sensitivity at end of segment k to the match pt
    I6 = Matrix{Float64}(I, 6, 6)
    tail_fwd = Vector{Matrix{Float64}}(undef, N2)
    tail_bwd = Vector{Matrix{Float64}}(undef, N2)
    tail_fwd[N2] = I6
    tail_bwd[N2] = I6
    for k in (N2-1):-1:1
        tail_fwd[k] = tail_fwd[k+1] * seg_stm_fwd[k+1]
        tail_bwd[k] = tail_bwd[k+1] * seg_stm_bwd[k+1]
    end

    # Cumulative STM from control point to match point
    cum_fwd = tail_fwd[1] * seg_stm_fwd[1]   # ∂X_fwd / ∂x0  (6×6)
    cum_bwd = tail_bwd[1] * seg_stm_bwd[1]   # ∂X_bwd / ∂xf  (6×6)

    # ── Velocity impulse sensitivity matrix B_k = [0_{3×3}; I_{3×3}] ────────
    # post-impulse velocity sensitivity to an applied Δv: ∂x_after/∂Δv = B
    B = vcat(zeros(3,3), I(3))   # 6×3

    # ─────────────────────────────────────────────────────────────────────────
    # 1. ∂c_mp / ∂x0  (7×6)
    #    = [- ∂X_fwd/∂x0;  - ∂m_fwd/∂x0]
    #    mass at match point is independent of initial state (CSI model, Eq.33)
    # ─────────────────────────────────────────────────────────────────────────
    J_x0_pos = -cum_fwd               # 6×6
    J_x0_mass = zeros(1, 6)           # ∂m_fwd/∂x0 = 0 for CSI
    J_x0 = vcat(J_x0_pos, J_x0_mass) # 7×6

    # ─────────────────────────────────────────────────────────────────────────
    # 2. ∂c_mp / ∂xf  (7×6)
    #    = [+ ∂X_bwd/∂xf;  + ∂m_bwd/∂xf]
    # ─────────────────────────────────────────────────────────────────────────
    J_xf_pos  = +cum_bwd
    J_xf_mass = zeros(1, 6)
    J_xf = vcat(J_xf_pos, J_xf_mass)

    # ─────────────────────────────────────────────────────────────────────────
    # 3. ∂c_mp / ∂u_fwd  (7 × 3*N2)
    #
    # For segment k, column block (7×3):
    #   Position/velocity (6×3):
    #     "direct" term:  tail[k] · stm2[k] · B · Δvmax_k     (Eq. 25)
    #     "mass correction" from downstream k' > k:
    #       sum_{k'>k} tail[k'] · stm2[k'] · B · outer(u_{k'}) · α_{kk'}
    #       where α_{kk'} = -Δvmax_{k'}/m_{k'} · ∂m_{k'}/∂||u_k|| · ∂||u_k||/∂u_{jk}
    #              ∂m_{k'}/∂u_{jk} = -mdot_max·D·Δt · u_{jk}/(||u_k||+ε)  (Eqs. 26,31)
    #   Mass row (1×3):
    #     ∂m_fwd_match/∂u_k = -mdot_max·D·Δt · u_k/(||u_k||+ε)  (Eq. 31)
    #     propagated through chain: for CSI, ∂m_{k'}/∂u_k = const for k'>k,
    #     so ∂m_match/∂u_k = -mdot_max·D·Δt · u_k/(||u_k||+ε)
    # ─────────────────────────────────────────────────────────────────────────
    J_uf = zeros(7, 3 * N2)
    for k in 1:N2
        u_k    = p._u_fwd[:, k]
        nu_k   = sqrt(dot(u_k, u_k) + ε)
        m_k    = Fwd.m[k]
        dv_k   = Fwd.dvmax[k]
        col    = (k-1)*3 + 1 : k*3

        # "direct" term at segment k: tail[k] · stm2[k] · B · Δvmax_k  (6×3)
        sv_k = tail_fwd[k] * Fwd.stm2[k] * B   # 6×3 sensitivity matrix at k
        J_uf[1:6, col] .-= sv_k * dv_k      # c_mp = X_bwd - X_fwd → ∂c/∂u_fwd = -∂X_fwd/∂u_fwd

        # mass correction terms from downstream segments k' > k
        # ∂m_{k'}/∂u_{jk} = -mdot_max·D·Δt · u_{jk}/(||u_k||+ε)
        # So correction to velocity at k' via Δvmax_{k'}: -Δvmax_{k'}/m_{k'} times same
        dm_uk_factor = mdot_max * D * dt / nu_k   # scalar: ∂m_k/∂||u_k|| · ∂||u_k||/∂u_k = u_k/nu_k
        for kp in (k+1):N2
            dv_kp = Fwd.dvmax[kp]
            m_kp  = Fwd.m[kp]
            u_kp  = p._u_fwd[:, kp]
            # ∂X_fwd/∂u_{jk} from Δvmax_{kp} variation: (6×1)·(1×3) outer product
            sv_kp = tail_fwd[kp] * Fwd.stm2[kp] * B   # 6×3
            # ∂Δvmax_{kp}/∂m_{kp} · ∂m_{kp}/∂u_k = (-dv_kp/m_kp) · (-dm_uk_factor · u_k^T)
            # Net sign: positive (two negatives)
            # Contribution to ∂X_fwd/∂u_k (6×3): sv_kp · u_kp · (dv_kp/m_kp * dm_uk_factor) · u_k^T
            coeff = (dv_kp / m_kp) * dm_uk_factor
            J_uf[1:6, col] .-= coeff * (sv_kp * u_kp) * u_k'  # (6×1)·(1×3) outer product → 6×3
        end

        # mass row: ∂m_fwd_match/∂u_k (Eq. 31)
        # The match-point mass = m after all N2 forward segments.
        # Since m_{k+1..N2} each have ∂m/∂u_k the same (all depend on m_k which depends on ||u_k||)
        # For CSI: ∂m_match/∂u_{jk} = -mdot_max·D·Δt · u_{jk}/(||u_k||+ε)
        J_uf[7, col] .+= (mdot_max * D * dt / nu_k) .* u_k   # sign: + because c_mp = m_bwd - m_fwd
    end

    # ─────────────────────────────────────────────────────────────────────────
    # 4. ∂c_mp / ∂u_bwd  (7 × 3*N2)  — same structure, opposite sign in c_mp
    # ─────────────────────────────────────────────────────────────────────────
    J_ub = zeros(7, 3 * N2)
    for k in 1:N2
        u_k    = p._u_bwd[:, k]
        nu_k   = sqrt(dot(u_k, u_k) + ε)
        m_k    = Bwd.m[k]
        dv_k   = Bwd.dvmax[k]
        col    = (k-1)*3 + 1 : k*3

        sv_k = tail_bwd[k] * Bwd.stm2[k] * B
        J_ub[1:6, col] .+= sv_k * dv_k

        dm_uk_factor = mdot_max * D * dt / nu_k
        for kp in (k+1):N2
            dv_kp = Bwd.dvmax[kp]
            m_kp  = Bwd.m[kp]
            u_kp  = p._u_bwd[:, kp]
            sv_kp = tail_bwd[kp] * Bwd.stm2[kp] * B
            coeff = (dv_kp / m_kp) * dm_uk_factor
            J_ub[1:6, col] .-= coeff * (sv_kp * u_kp) * u_k'
        end

        J_ub[7, col] .+= (mdot_max * D * dt / nu_k) .* u_k   # bwd adds mass → ∂m_bwd/∂u_bwd > 0
    end

    # ─────────────────────────────────────────────────────────────────────────
    # 5. ∂c_mp / ∂m0  (7-vector)
    #    Position/velocity: ∂X_fwd/∂m0 via Δvmax_k = Tmax·D·nT·Δt/m_k dependence
    #      ∂X_fwd/∂m0 = sum_k {tail[k] · stm2[k] · B · u_k · (-Δvmax_k/m_k)}
    #      (since ∂m_k/∂m0 = 1 for all k, from Eq. 33)
    #    Mass: ∂m_fwd_match/∂m0 = 1.0  (Eq. 33, chain through all segments)
    # ─────────────────────────────────────────────────────────────────────────
    dm0_pos = zeros(6)
    for k in 1:N2
        u_k   = p._u_fwd[:, k]
        m_k   = Fwd.m[k]
        dv_k  = Fwd.dvmax[k]
        sv_k  = tail_fwd[k] * Fwd.stm2[k] * B   # 6×3
        dm0_pos .+= sv_k * (u_k * (-dv_k / m_k))  # correction per segment
    end
    J_dm0 = vcat(-dm0_pos, [-1.0])   # 7-vector; c_mp = X_bwd - X_fwd → sign = -∂X_fwd

    # ─────────────────────────────────────────────────────────────────────────
    # 6. ∂c_mp / ∂mf  (7-vector) — same for backward half-phase
    # ─────────────────────────────────────────────────────────────────────────
    dmf_pos = zeros(6)
    for k in 1:N2
        u_k   = p._u_bwd[:, k]
        m_k   = Bwd.m[k]
        dv_k  = Bwd.dvmax[k]
        sv_k  = tail_bwd[k] * Bwd.stm2[k] * B
        dmf_pos .+= sv_k * (u_k * (-dv_k / m_k))
    end
    J_dmf = vcat(+dmf_pos, [+1.0])   # c_mp = X_bwd - X_fwd → +∂X_bwd

    # ─────────────────────────────────────────────────────────────────────────
    # 7. ∂c_mp / ∂tf  and  ∂c_mp / ∂t0  (7-vectors each)
    #
    # Δt = (tf - t0)/N for each full segment, dt2 = Δt/2 for each half-prop.
    # ∂Δt/∂tf = +1/N,   ∂Δt/∂t0 = -1/N
    #
    # Two contributions per segment k:
    #   (a) Kepler time partials: propagation duration changes
    #       ∂X_match/∂tf |_via_prop_k = tail[k]·stm2[k]·dxdt1_k·(∂dt2/∂tf)
    #                                  + tail[k]·dxdt2_k·(∂dt2/∂tf)
    #   (b) Δvmax change: Δvmax_k = Tmax·D·nT·Δt/m_k → ∂Δvmax_k/∂tf = Δvmax_k/Δt · (∂Δt/∂tf)
    #       ∂X_match/∂tf |_via_dvmax_k = tail[k]·stm2[k]·B·u_k · (Δvmax_k/Δt) · (1/N)
    #
    # Total ∂X_fwd/∂tf = sum_k [tail[k]·stm2[k]·dxdt1[k]·(1/(2N)) + tail[k]·dxdt2[k]·(1/(2N))
    #                           + tail[k]·stm2[k]·B·u_k·(Δvmax_k/Δt)·(1/N)]
    # ─────────────────────────────────────────────────────────────────────────
    N_total = p.transcription.n_segments
    ddtf = 1.0 / N_total         # ∂Δt/∂tf
    ddt0 = -1.0 / N_total        # ∂Δt/∂t0
    ddt2f = ddtf / 2.0           # ∂(dt/2)/∂tf
    ddt2_0 = ddt0 / 2.0          # ∂(dt/2)/∂t0

    dtf_fwd = zeros(6)
    dt0_fwd = zeros(6)
    for k in 1:N2
        T2k   = tail_fwd[k] * Fwd.stm2[k]   # 6×6
        Tk    = tail_fwd[k]                   # 6×6

        # Kepler time partial contributions
        dtf_fwd .+= T2k * Fwd.dxdt1[k] * ddt2f + Tk * Fwd.dxdt2[k] * ddt2f
        dt0_fwd .+= T2k * Fwd.dxdt1[k] * ddt2_0 + Tk * Fwd.dxdt2[k] * ddt2_0

        # Δvmax time partial contribution
        # dvmax_k = Tmax·D·nT·dt/m_k, so ∂dvmax_k/∂dt = Tmax·D·nT/m_k - dvmax_k/m_k · ∂m_k/∂dt
        # ∂m_k/∂dt = -mdot_max·D · sum_{j<k} norm(u_j)   (sign_m = -1 for forward)
        # ⟹ ∂dvmax_k/∂dt = dvmax_k/dt - dvmax_k/m_k · ∂m_k/∂dt
        #                 = dvmax_k/dt + dvmax_k/m_k · mdot_max·D · sum_{j<k} norm(u_j)
        u_k    = p._u_fwd[:, k]
        dv_k   = Fwd.dvmax[k]
        m_k    = Fwd.m[k]
        dm_k_ddt = -sum(mdot_max * D * sqrt(dot(p._u_fwd[:, j], p._u_fwd[:, j]) + ε) for j in 1:k-1; init=0.0)
        ddvmax_k_ddt = dv_k / dt - dv_k / m_k * dm_k_ddt   # ∂dvmax_k/∂dt
        sv_k  = T2k * B * u_k   # 6-vector: sensitivity to Δv
        dtf_fwd .+= sv_k * (ddvmax_k_ddt * ddtf)
        dt0_fwd .+= sv_k * (ddvmax_k_ddt * ddt0)
    end

    dtf_bwd = zeros(6)
    dt0_bwd = zeros(6)
    for k in 1:N2
        T2k = tail_bwd[k] * Bwd.stm2[k]
        Tk  = tail_bwd[k]

        # dxdt1/dxdt2 from backward propagator = ∂x/∂(-dt2)  (signed).
        # ∂(-dt2)/∂tf = -ddt2f,  ∂(-dt2)/∂t0 = +ddt2f
        dtf_bwd .+= T2k * Bwd.dxdt1[k] * (-ddt2f) + Tk * Bwd.dxdt2[k] * (-ddt2f)
        dt0_bwd .+= T2k * Bwd.dxdt1[k] * (+ddt2f) + Tk * Bwd.dxdt2[k] * (+ddt2f)

        # Backward sign_m = +1, so m_k increases with dt:
        # ∂m_k/∂dt = +mdot_max·D · sum_{j<k} norm(u_j)
        u_k    = p._u_bwd[:, k]
        dv_k   = Bwd.dvmax[k]
        m_k    = Bwd.m[k]
        dm_k_ddt = +sum(mdot_max * D * sqrt(dot(p._u_bwd[:, j], p._u_bwd[:, j]) + ε) for j in 1:k-1; init=0.0)
        ddvmax_k_ddt = dv_k / dt - dv_k / m_k * dm_k_ddt
        sv_k  = T2k * B * u_k
        dtf_bwd .+= sv_k * (ddvmax_k_ddt * ddtf)
        dt0_bwd .+= sv_k * (ddvmax_k_ddt * ddt0)
    end

    # c_mp = X_bwd - X_fwd → ∂c_mp/∂tf = ∂X_bwd/∂tf - ∂X_fwd/∂tf
    # For mass: m_match_fwd and m_match_bwd both change with Δt through mflow
    #   ∂m_fwd/∂tf = -sum_k(mdot_max·D·||u_k||) · ddtf   (Eq. 14 with ∂Δt/∂tf)
    #   ∂m_bwd/∂tf = +sum_k(mdot_max·D·||u_k||) · ddtf
    dm_fwd_dtf = -sum(mdot_max * D * sqrt(dot(p._u_fwd[:, k], p._u_fwd[:, k]) + ε) * ddtf for k in 1:N2)
    dm_bwd_dtf = +sum(mdot_max * D * sqrt(dot(p._u_bwd[:, k], p._u_bwd[:, k]) + ε) * ddtf for k in 1:N2)
    dm_fwd_dt0 = -sum(mdot_max * D * sqrt(dot(p._u_fwd[:, k], p._u_fwd[:, k]) + ε) * ddt0 for k in 1:N2)
    dm_bwd_dt0 = +sum(mdot_max * D * sqrt(dot(p._u_bwd[:, k], p._u_bwd[:, k]) + ε) * ddt0 for k in 1:N2)

    J_dtf = vcat(dtf_bwd .- dtf_fwd, [dm_bwd_dtf - dm_fwd_dtf])
    J_dt0 = vcat(dt0_bwd .- dt0_fwd, [dm_bwd_dt0 - dm_fwd_dt0])

    # ── Ephemeris time-partial contributions ─────────────────────────────────
    # x0 = [r_planet(t0); v_planet(t0) + v∞_dep]
    # ∂x0/∂t0 = [v_planet(t0); a_planet(t0)]   →  6-vector
    # Chain rule: ∂c_mp/∂t0|_full += J_x0 · ∂x0/∂t0
    # (and symmetrically for tf / xf)
    if p.ephemeris_left !== nothing
        v0_p, a0_p = p._planet_left[2], p._planet_left[3]
        J_dt0 .+= J_x0 * vcat(v0_p, a0_p)
    end
    if p.ephemeris_right !== nothing
        vf_p, af_p = p._planet_right[2], p._planet_right[3]
        J_dtf .+= J_xf * vcat(vf_p, af_p)
    end

    # ── v∞ Jacobians: J_x0/xf restricted to velocity columns ────────────────
    # ∂x0/∂v∞_dep = [0_{3×3}; I_{3×3}]  →  ∂c_mp/∂v∞_dep = J_x0[:, 4:6]
    J_vinf_dep = J_x0[:, 4:6]   # 7×3
    J_vinf_arr = J_xf[:, 4:6]   # 7×3

    return (dvinf_dep = J_vinf_dep,
            dvinf_arr = J_vinf_arr,
            dx0    = J_x0,
            dxf    = J_xf,
            du_fwd = J_uf,
            du_bwd = J_ub,
            dm0    = J_dm0,
            dmf    = J_dmf,
            dtf    = J_dtf,
            dt0    = J_dt0)
end

# ═════════════════════════════════════════════════════════════════════════════
# NLP Interface
#
# Implements the same variable_list / function_list / get_decision_vector /
# get_functions / get_constraint_bounds / get_variable_bounds interface as
# CollocationPhase, so SimsFlanaganPhase can be managed by the same
# PhaseManager / SequenceManager infrastructure.
#
# DEPENDENCY: PhaseFunction, BoundaryConstraint, BoundaryFunction, and
# DirectSolverVariable are defined in OptControlStubs.jl.  Those types must
# have been loaded before any of the interface functions below are called.
#
# Decision vector layout  (total length = 8 + 6·N2):
#   vinf_dep  [3]     indices  1 : 3
#   vinf_arr  [3]     indices  4 : 6
#   u_fwd     [3×N2]  indices  7 : 6+3·N2
#   u_bwd     [3×N2]  indices  6+3·N2+1 : 6+6·N2
#   t0        [1]     index    6+6·N2+1
#   tf        [1]     index    6+6·N2+2
#   m0        [1]     index    6+6·N2+3
#   mf        [1]     index    6+6·N2+4
#
# Function vector layout  (total length = 7 + n_bc):
#   matchpoint_defect [7]   equality constraints (lb = ub = 0)
#   boundary constraints    user-supplied (variable bounds)
# ═════════════════════════════════════════════════════════════════════════════

# ─────────────────────────────────────────────────────────────────────────────
# Variable type tags  (used as the `var=` value of a DirectSolverVariable)
#
# These are singleton types; the identity of the decision block lives in the
# DirectSolverVariable wrapper.  They play the same role as AbstractStateArray /
# AbstractControlArray / AbstractTime / AbstractParameter in collocation.
# ─────────────────────────────────────────────────────────────────────────────

struct SFVInfinity3    end   # v∞ departure or arrival:   3 scalars (no tiling)
struct SFStateBlock    end   # x0 or xf: a free endpoint state, 6 scalars.
                             # The alternative to a planetary endpoint — the
                             # transcription is two half-propagations meeting at
                             # a match point, and where the ends come from is a
                             # separate question from how the middle is solved.
struct SFThrustBlock   end   # u_fwd or u_bwd:  per-component bounds × N2 segments
struct SFMassParam     end   # m0 or mf:  scalar mass parameter
struct SFTime          end   # t0 or tf:  scalar time  (alternative to AbstractTime)
struct SFMatchPointBlock end # sentinel in function_list (like CollocationPhase's DefectBlock)

# ─────────────────────────────────────────────────────────────────────────────
# variable_list
#
# Returns the DirectSolverVariables in NLP order, skipping any that are nothing.
# ─────────────────────────────────────────────────────────────────────────────

function variable_list(p::SimsFlanaganPhase)
    filter(!isnothing, Any[
        p.x0_var,       p.xf_var,
        p.vinf_dep_var, p.vinf_arr_var,
        p.u_fwd_var,    p.u_bwd_var,
        p.t0_var,       p.tf_var,
        p.m0_var,       p.mf_var,
    ])
end

# ─────────────────────────────────────────────────────────────────────────────
# nlp_length  —  number of NLP scalars contributed by variable v
# ─────────────────────────────────────────────────────────────────────────────

function nlp_length(p::SimsFlanaganPhase, v)
    N2 = n_fwd(p)
    vt = v.var
    vt isa SFStateBlock  && return 6
    vt isa SFVInfinity3  && return 3
    vt isa SFThrustBlock && return 3 * N2
    vt isa SFMassParam   && return 1
    (vt isa SFTime || vt isa AbstractTime) && return 1
    throw(ArgumentError(
        "a SimsFlanaganPhase variable must be a state, control, mass, time or v-infinity " *
        "variable; got $(typeof(vt))"))
end

nlp_length(p::SimsFlanaganPhase) = sum(nlp_length(p, v) for v in variable_list(p))

# ─────────────────────────────────────────────────────────────────────────────
# nlp_bounds  —  (lb, ub) for the NLP scalars contributed by variable v
#
# For SFThrustBlock, v.lower_bounds / v.upper_bounds are 3-element
# (per-component) and are tiled N2 times to produce 3×N2 scalars.
# For all other types the bounds are used as-is.
# ─────────────────────────────────────────────────────────────────────────────

# _var_nlp_scale / _var_nlp_shift
# Full-length scale/shift vector for variable v in NLP space.
# SFThrustBlock has 3-component scale/shift tiled N2 times.
function _var_nlp_scale(p::SimsFlanaganPhase, v::DirectSolverVariable)
    v.var isa SFThrustBlock && return repeat(v.scale, n_fwd(p))
    return v.scale
end

function _var_nlp_shift(p::SimsFlanaganPhase, v::DirectSolverVariable)
    v.var isa SFThrustBlock && return repeat(v.shift, n_fwd(p))
    return v.shift
end

function nlp_bounds(p::SimsFlanaganPhase, v)
    lb = v.var isa SFThrustBlock ? repeat(v.lower_bounds, n_fwd(p)) : copy(v.lower_bounds)
    ub = v.var isa SFThrustBlock ? repeat(v.upper_bounds, n_fwd(p)) : copy(v.upper_bounds)
    sc = _var_nlp_scale(p, v)
    sh = _var_nlp_shift(p, v)
    return (lb .- sh) ./ sc, (ub .- sh) ./ sc
end

# ─────────────────────────────────────────────────────────────────────────────
# get_variable_bounds  —  flat (lb, ub) for the full decision vector
# ─────────────────────────────────────────────────────────────────────────────

function get_variable_bounds(p::SimsFlanaganPhase)
    lb = Float64[]
    ub = Float64[]
    for v in variable_list(p)
        lo, hi = nlp_bounds(p, v)
        append!(lb, lo)
        append!(ub, hi)
    end
    return lb, ub
end

# ─────────────────────────────────────────────────────────────────────────────
# variable_ranges  —  phase-local column slices matching variable_list order
# ─────────────────────────────────────────────────────────────────────────────

function variable_ranges(p::SimsFlanaganPhase)
    ranges = UnitRange{Int}[]
    offset = 0
    for v in variable_list(p)
        len = nlp_length(p, v)
        push!(ranges, offset+1 : offset+len)
        offset += len
    end
    ranges
end

# ─────────────────────────────────────────────────────────────────────────────
# get_decision_vector  —  pack cached phase values into a flat vector
# ─────────────────────────────────────────────────────────────────────────────

function get_decision_vector(p::SimsFlanaganPhase)
    N2 = n_fwd(p)
    x  = zeros(nlp_length(p))
    off = 0
    if p.x0_var !== nothing
        v = p.x0_var
        x[off+1 : off+6] .= (p._x0 .- v.shift) ./ v.scale;                 off += 6
    end
    if p.xf_var !== nothing
        v = p.xf_var
        x[off+1 : off+6] .= (p._xf .- v.shift) ./ v.scale;                 off += 6
    end
    if p.vinf_dep_var !== nothing
        v = p.vinf_dep_var
        x[off+1 : off+3] .= (p._vinf_dep .- v.shift) ./ v.scale;           off += 3
    end
    if p.vinf_arr_var !== nothing
        v = p.vinf_arr_var
        x[off+1 : off+3] .= (p._vinf_arr .- v.shift) ./ v.scale;           off += 3
    end
    if p.u_fwd_var !== nothing
        v = p.u_fwd_var
        sc = repeat(v.scale, N2); sh = repeat(v.shift, N2)
        x[off+1 : off+3N2] .= (vec(p._u_fwd) .- sh) ./ sc;                 off += 3N2
    end
    if p.u_bwd_var !== nothing
        v = p.u_bwd_var
        sc = repeat(v.scale, N2); sh = repeat(v.shift, N2)
        x[off+1 : off+3N2] .= (vec(p._u_bwd) .- sh) ./ sc;                 off += 3N2
    end
    if p.t0_var !== nothing
        v = p.t0_var
        x[off+1] = (p._t0 - v.shift[1]) / v.scale[1];                      off += 1
    end
    if p.tf_var !== nothing
        v = p.tf_var
        x[off+1] = (p._tf - v.shift[1]) / v.scale[1];                      off += 1
    end
    if p.m0_var !== nothing
        v = p.m0_var
        x[off+1] = (p._m0 - v.shift[1]) / v.scale[1];                      off += 1
    end
    if p.mf_var !== nothing
        v = p.mf_var
        x[off+1] = (p._mf - v.shift[1]) / v.scale[1];                      off += 1
    end
    return x
end

# ─────────────────────────────────────────────────────────────────────────────
# set_decision_vector!  —  unpack a flat vector into the cached phase fields.
# Called by the framework at every NLP iteration before get_functions.
# ─────────────────────────────────────────────────────────────────────────────

function set_decision_vector!(p::SimsFlanaganPhase, x::Vector{Float64})
    N2  = n_fwd(p)
    off = 0
    if p.x0_var !== nothing
        v = p.x0_var
        p._x0 = x[off+1 : off+6] .* v.scale .+ v.shift;                     off += 6
        v.value = copy(p._x0)
    end
    if p.xf_var !== nothing
        v = p.xf_var
        p._xf = x[off+1 : off+6] .* v.scale .+ v.shift;                     off += 6
        v.value = copy(p._xf)
    end
    if p.vinf_dep_var !== nothing
        v = p.vinf_dep_var
        p._vinf_dep = x[off+1 : off+3] .* v.scale .+ v.shift;               off += 3
        v.value = copy(p._vinf_dep)
    end
    if p.vinf_arr_var !== nothing
        v = p.vinf_arr_var
        p._vinf_arr = x[off+1 : off+3] .* v.scale .+ v.shift;               off += 3
        v.value = copy(p._vinf_arr)
    end
    if p.u_fwd_var !== nothing
        v = p.u_fwd_var
        sc = repeat(v.scale, N2); sh = repeat(v.shift, N2)
        p._u_fwd = reshape(x[off+1 : off+3N2] .* sc .+ sh, 3, N2);         off += 3N2
        v.value = vec(p._u_fwd)
    end
    if p.u_bwd_var !== nothing
        v = p.u_bwd_var
        sc = repeat(v.scale, N2); sh = repeat(v.shift, N2)
        p._u_bwd = reshape(x[off+1 : off+3N2] .* sc .+ sh, 3, N2);         off += 3N2
        v.value = vec(p._u_bwd)
    end
    if p.t0_var !== nothing
        v = p.t0_var
        p._t0 = x[off+1] * v.scale[1] + v.shift[1];                         off += 1
        v.value = [p._t0]
    end
    if p.tf_var !== nothing
        v = p.tf_var
        p._tf = x[off+1] * v.scale[1] + v.shift[1];                         off += 1
        v.value = [p._tf]
    end
    if p.m0_var !== nothing
        v = p.m0_var
        p._m0 = x[off+1] * v.scale[1] + v.shift[1];                         off += 1
        v.value = [p._m0]
    end
    if p.mf_var !== nothing
        v = p.mf_var
        p._mf = x[off+1] * v.scale[1] + v.shift[1];                         off += 1
        v.value = [p._mf]
    end
    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# SFBoundaryContext  —  value object passed to boundary constraint closures.
#
# Carries all S-F phase state so closures can read v∞, orbital states,
# times, and masses without capturing the mutable phase directly.
# ─────────────────────────────────────────────────────────────────────────────

struct SFBoundaryContext
    x0       # [6]  left  ctrl pt spacecraft state
    xf       # [6]  right ctrl pt spacecraft state
    vinf_dep # [3]  departure v∞
    vinf_arr # [3]  arrival   v∞
    t0
    tf
    m0
    mf
    u_fwd    # [3, N/2]  forward  throttle vectors
    u_bwd    # [3, N/2]  backward throttle vectors
end

"""
    SFPathContext

One segment of a Sims-Flanagan phase, handed to a path constraint. A path
constraint is evaluated once per segment, so what it sees is that segment's
throttle rather than the whole matrix.

`half` is `:forward` or `:backward`. The two halves are propagated from
opposite ends and meet at the match point, so a constraint that cares which
end it is near can ask.
"""
struct SFPathContext{U}
    control  ::U
    segment  ::Int
    half     ::Symbol
end

_sf_boundary_context(p::SimsFlanaganPhase) =
    SFBoundaryContext(copy(p._x0), copy(p._xf),
                      copy(p._vinf_dep), copy(p._vinf_arr),
                      p._t0, p._tf, p._m0, p._mf,
                      copy(p._u_fwd), copy(p._u_bwd))

# Generic shooting-sequence dispatch — SF↔SF (and SF↔MGA, ZOH↔SF) linkage
# constraints look up the boundary context via this method.
_shoot_boundary_context(p::SimsFlanaganPhase) = _sf_boundary_context(p)

# ─────────────────────────────────────────────────────────────────────────────
# _sf_context_replace / _bc_jacobian_ad_chunk
#
# The counterpart of an analytic Jacobian registered with add_jacobian!. Every
# other phase type falls back to automatic differentiation when none is
# registered; this one returned zeros, which is a wrong derivative rather than a
# missing one.
#
# SFBoundaryContext has untyped fields, so a dual flows into it without the
# promotion MGAnDSMs needs for its parameterised context.
#
# Differentiates with respect to the *physical* value, so the caller applies the
# same scaling it applies to a registered Jacobian.
# ─────────────────────────────────────────────────────────────────────────────

function _sf_context_replace(ctx::SFBoundaryContext, p::SimsFlanaganPhase,
                             var::DirectSolverVariable, phys::AbstractVector)
    oid = objectid(var)
    x0  = ctx.x0;       xf  = ctx.xf
    vd  = ctx.vinf_dep; va  = ctx.vinf_arr
    t0  = ctx.t0;       tf  = ctx.tf
    m0  = ctx.m0;       mf  = ctx.mf
    uf  = ctx.u_fwd;    ub  = ctx.u_bwd
    N2  = n_fwd(p)
    !isnothing(p.x0_var)       && oid == objectid(p.x0_var)       && (x0 = phys)
    !isnothing(p.xf_var)       && oid == objectid(p.xf_var)       && (xf = phys)
    !isnothing(p.vinf_dep_var) && oid == objectid(p.vinf_dep_var) && (vd = phys)
    !isnothing(p.vinf_arr_var) && oid == objectid(p.vinf_arr_var) && (va = phys)
    !isnothing(p.t0_var)       && oid == objectid(p.t0_var)       && (t0 = phys[1])
    !isnothing(p.tf_var)       && oid == objectid(p.tf_var)       && (tf = phys[1])
    !isnothing(p.m0_var)       && oid == objectid(p.m0_var)       && (m0 = phys[1])
    !isnothing(p.mf_var)       && oid == objectid(p.mf_var)       && (mf = phys[1])
    !isnothing(p.u_fwd_var)    && oid == objectid(p.u_fwd_var)    && (uf = reshape(phys, 3, N2))
    !isnothing(p.u_bwd_var)    && oid == objectid(p.u_bwd_var)    && (ub = reshape(phys, 3, N2))
    SFBoundaryContext(x0, xf, vd, va, t0, tf, m0, mf, uf, ub)
end

function _bc_jacobian_ad_chunk(p::SimsFlanaganPhase, bf::BoundaryFunction,
                               var::DirectSolverVariable)
    ctx_f64 = _sf_boundary_context(p)
    return ForwardDiff.jacobian(
        phys -> bf.fn(_sf_context_replace(ctx_f64, p, var, phys)),
        copy(var.value))
end

# ─────────────────────────────────────────────────────────────────────────────
# function_list
#
# Declares the NLP function structure at setup time.  SequenceManager/
# PhaseManager iterates this to assign global row offsets and collect bounds.
#
# Requires PhaseFunction (from OptControlStubs.jl) at call time.
# ─────────────────────────────────────────────────────────────────────────────

function function_list(p::SimsFlanaganPhase)
    # Match-point residuals: 7 equality constraints (lb = ub = 0)
    result = PhaseFunction[
        PhaseFunction(SFMatchPointBlock(), 7, zeros(7), zeros(7), "matchpoint_defect")
    ]
    for c in p.constraints
        if c isa BoundaryConstraint
            n = length(c.lower_bounds)
            push!(result, PhaseFunction(c, n, c.lower_bounds, c.upper_bounds,
                                        c.calc.name))
        end
    end
    return result
end

# ─────────────────────────────────────────────────────────────────────────────
# get_functions
#
# Called by the framework at every NLP iteration (after set_decision_vector!).
# Propagates both half-phases, then returns the concatenated residual vector:
#   [matchpoint_defect(7);  boundary_constraint_vals...]
# ─────────────────────────────────────────────────────────────────────────────

function get_functions(p::SimsFlanaganPhase)
    evaluate_matchpoint!(p)
    F = (matchpoint_defect(p) .- p.matchpoint_shift) ./ p.matchpoint_scale

    for c in p.constraints
        if c isa BoundaryConstraint
            ctx  = _sf_boundary_context(p)
            vals = applicable(c.calc.fn, ctx) ? c.calc.fn(ctx) : c.calc.fn()
            append!(F, (vals .- c.shift) ./ c.scale)
        end
    end
    return F
end

# ─────────────────────────────────────────────────────────────────────────────
# get_constraint_bounds  —  (lb, ub) aligned with get_functions output
# ─────────────────────────────────────────────────────────────────────────────

function get_constraint_bounds(p::SimsFlanaganPhase)
    # match-point defect equality constraints scaled to NLP units
    lb = (zeros(7) .- p.matchpoint_shift) ./ p.matchpoint_scale
    ub = copy(lb)
    for c in p.constraints
        if c isa BoundaryConstraint
            append!(lb, (c.lower_bounds .- c.shift) ./ c.scale)
            append!(ub, (c.upper_bounds .- c.shift) ./ c.scale)
        end
    end
    return lb, ub
end

n_constraints(p::SimsFlanaganPhase) =
    7 + sum(length(c.lower_bounds) for c in p.constraints
            if c isa BoundaryConstraint; init = 0)

# ─────────────────────────────────────────────────────────────────────────────
# set_initial_guess!
#
# Sets phase-cached decision variable values and syncs var.value fields in
# one call.  Any keyword omitted leaves the corresponding field unchanged.
#
#   set_initial_guess!(phase;
#       vinf_dep = [...],   # 3-vector (km/s)
#       vinf_arr = [...],   # 3-vector (km/s)
#       u_fwd    = [...],   # 3×N2 matrix (unit throttle)
#       u_bwd    = [...],   # 3×N2 matrix (unit throttle)
#       t0       = 0.0,     # scalar (s)
#       tf       = ...,     # scalar (s)
#       m0       = ...,     # scalar (kg)
#       mf       = ...,     # scalar (kg)
#   )
# ─────────────────────────────────────────────────────────────────────────────

function set_initial_guess!(p::SimsFlanaganPhase;
                             vinf_dep = nothing,
                             vinf_arr = nothing,
                             u_fwd    = nothing,
                             u_bwd    = nothing,
                             t0       = nothing,
                             tf       = nothing,
                             m0       = nothing,
                             mf       = nothing)
    if vinf_dep !== nothing
        p._vinf_dep = Float64.(vinf_dep)
        p.vinf_dep_var !== nothing && (p.vinf_dep_var.value = copy(p._vinf_dep))
    end
    if vinf_arr !== nothing
        p._vinf_arr = Float64.(vinf_arr)
        p.vinf_arr_var !== nothing && (p.vinf_arr_var.value = copy(p._vinf_arr))
    end
    if u_fwd !== nothing
        p._u_fwd = Float64.(u_fwd)
        p.u_fwd_var !== nothing && (p.u_fwd_var.value = vec(p._u_fwd))
    end
    if u_bwd !== nothing
        p._u_bwd = Float64.(u_bwd)
        p.u_bwd_var !== nothing && (p.u_bwd_var.value = vec(p._u_bwd))
    end
    if t0 !== nothing
        p._t0 = Float64(t0)
        p.t0_var !== nothing && (p.t0_var.value = [p._t0])
    end
    if tf !== nothing
        p._tf = Float64(tf)
        p.tf_var !== nothing && (p.tf_var.value = [p._tf])
    end
    if m0 !== nothing
        p._m0 = Float64(m0)
        p.m0_var !== nothing && (p.m0_var.value = [p._m0])
    end
    if mf !== nothing
        p._mf = Float64(mf)
        p.mf_var !== nothing && (p.mf_var.value = [p._mf])
    end
    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# sparsity_structure
#
# Block-level dependency: does function_i depend on variable_j?
# Returns Vector{Bool} of length length(variable_list(p)) per PhaseFunction,
# or Matrix{Bool} (n_functions × n_variables) for the whole phase.
# ─────────────────────────────────────────────────────────────────────────────

function sparsity_structure(p::SimsFlanaganPhase, pf::PhaseFunction)
    vlist = variable_list(p)
    n     = length(vlist)
    if pf.source isa SFMatchPointBlock
        # Match-point defect depends on all 8 decision variable blocks
        return fill(true, n)
    elseif pf.source isa BoundaryConstraint
        bc = pf.source
        # Default: depends on v∞, times, masses.
        # Also include any variable for which an analytic Jacobian was registered
        # (e.g. thrust blocks for a per-segment control-magnitude path constraint).
        return Bool[
            v === p.x0_var || v === p.xf_var ||
            v === p.vinf_dep_var || v === p.vinf_arr_var ||
            v === p.t0_var || v === p.tf_var ||
            v === p.m0_var || v === p.mf_var ||
            has_jacobian(bc.calc, v)
            for v in vlist
        ]
    end
    return fill(false, n)
end

function sparsity_structure(p::SimsFlanaganPhase)
    flist = function_list(p)
    vlist = variable_list(p)
    Bool[sparsity_structure(p, f)[j] for f in flist, j in eachindex(vlist)]
end

# ─────────────────────────────────────────────────────────────────────────────
# jacobian_chunk
#
# Returns the Jacobian block (Matrix{Float64}) for one (PhaseFunction, variable)
# pair.  Called by ShootingManager when assembling the global sparse Jacobian.
#
# For SFMatchPointBlock: delegates to matchpoint_jacobian (analytic).
# For BoundaryConstraint: uses registered analytic Jacobian if available,
#   otherwise returns a zero block (all non-zero Jacobians must be registered
#   via add_jacobian! before solve).
# ─────────────────────────────────────────────────────────────────────────────

function jacobian_chunk(p::SimsFlanaganPhase, pf::PhaseFunction,
                        var::DirectSolverVariable)
    var_sc = _var_nlp_scale(p, var)   # column scales for this variable block
    if pf.source isa SFMatchPointBlock
        jac = matchpoint_jacobian(p)
        oid = objectid(var)
        chunk_phys =
            !isnothing(p.x0_var)       && oid == objectid(p.x0_var)       ? jac.dx0                 :
            !isnothing(p.xf_var)       && oid == objectid(p.xf_var)       ? jac.dxf                 :
            !isnothing(p.vinf_dep_var) && oid == objectid(p.vinf_dep_var) ? jac.dvinf_dep           :
            !isnothing(p.vinf_arr_var) && oid == objectid(p.vinf_arr_var) ? jac.dvinf_arr           :
            !isnothing(p.u_fwd_var)    && oid == objectid(p.u_fwd_var)    ? jac.du_fwd              :
            !isnothing(p.u_bwd_var)    && oid == objectid(p.u_bwd_var)    ? jac.du_bwd              :
            !isnothing(p.t0_var) && oid == objectid(p.t0_var) ? reshape(jac.dt0, 7, 1)             :
            !isnothing(p.tf_var) && oid == objectid(p.tf_var) ? reshape(jac.dtf, 7, 1)             :
            !isnothing(p.m0_var) && oid == objectid(p.m0_var) ? reshape(jac.dm0, 7, 1)             :
            !isnothing(p.mf_var) && oid == objectid(p.mf_var) ? reshape(jac.dmf, 7, 1)             :
            throw(ArgumentError(
                "jacobian_chunk: $(repr(var.name)) must be a variable declared on phase " *
                ":$(p.name); declare it with Vary before asking for its Jacobian"))
        # Scale: J_nlp[i,j] = J_phys[i,j] * var_scale[j] / con_scale[i]
        return chunk_phys .* transpose(var_sc) ./ p.matchpoint_scale
    elseif pf.source isa BoundaryConstraint
        bc = pf.source
        chunk_phys = has_jacobian(bc.calc, var) ?
            get_jacobian(bc.calc, var) :
            _bc_jacobian_ad_chunk(p, bc.calc, var)
        return chunk_phys .* transpose(var_sc) ./ bc.scale
    end
    throw(ArgumentError("SimsFlanaganPhase.jacobian_chunk: unknown function source " *
        "$(typeof(pf.source)). Returning zeros here would be a wrong " *
        "derivative rather than a missing one."))
end


# ─────────────────────────────────────────────────────────────────────────────
# get_objective
#
# Evaluates the scalar Mayer objective.  fn receives an SFBoundaryContext.
# Returns the NLP-sense value (negated for :Max).
# ─────────────────────────────────────────────────────────────────────────────

function get_objective(p::SimsFlanaganPhase)
    isnothing(p.objective) && return 0.0
    obj = p.objective
    ctx = _sf_boundary_context(p)
    raw = applicable(obj.fn, ctx) ? obj.fn(ctx) : obj.fn()
    val = Float64(raw)
    return obj.sense === :Max ? -val : val
end

# ─────────────────────────────────────────────────────────────────────────────
# objective_gradient_chunk
#
# Returns ∂J/∂var as a Vector{Float64}.  Uses a registered analytic Jacobian
# if available; otherwise returns zeros (the variable does not affect the
# objective).  Register non-zero blocks via add_objective_jacobian!.
# Sign convention: negated for :Max.
# ─────────────────────────────────────────────────────────────────────────────

function objective_gradient_chunk(p::SimsFlanaganPhase, var::DirectSolverVariable)
    isnothing(p.objective) && return zeros(nlp_length(p, var))
    obj = p.objective

    # A registered analytic gradient takes priority.
    if has_objective_jacobian(obj, var)
        g_phys   = get_objective_jacobian(obj, var)
        g_scaled = g_phys .* _var_nlp_scale(p, var)
        return obj.sense === :Max ? .-g_scaled : g_scaled
    end

    # None registered. The objective may still depend on this variable, so
    # differentiate it rather than assume it does not. An objective closure that
    # takes no context cannot read a variable at all, and zero is right there.
    ctx_f64 = _sf_boundary_context(p)
    applicable(obj.fn, ctx_f64) || return zeros(nlp_length(p, var))
    g_phys = ForwardDiff.gradient(
        phys -> obj.fn(_sf_context_replace(ctx_f64, p, var, phys)), copy(var.value))
    g_scaled = g_phys .* _var_nlp_scale(p, var)
    return obj.sense === :Max ? .-g_scaled : g_scaled
end

# Quantities on a segment of a Sims-Flanagan phase. `control` is the one a
# path constraint reads; the thrust ball |u| <= 1 is written with it.
control(c::SFPathContext)              = c.control

# Quantities on a Sims-Flanagan phase.
state(p::SimsFlanaganPhase)            = p._x0
final_state(p::SimsFlanaganPhase)      = p._xf
forward_control(p::SimsFlanaganPhase)  = p._u_fwd
backward_control(p::SimsFlanaganPhase) = p._u_bwd
initial_time(p::SimsFlanaganPhase)     = p._t0
final_time(p::SimsFlanaganPhase)       = p._tf
initial_mass(p::SimsFlanaganPhase)     = p._m0
final_mass(p::SimsFlanaganPhase)       = p._mf
