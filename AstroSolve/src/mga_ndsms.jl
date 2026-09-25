# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0


# MGAnDSMs.jl
#
# Multiple Gravity Assists with n Deep-Space Maneuvers (MGAnDSMs) transcription
# for use with the Epicycle optimal control framework.
#
# Reference: Englander (PhD thesis) Section 3.2 / 3.4
#   "Multiple Gravity-Assist with n Deep-Space Maneuvers Trajectory Model"
#
# Phase structure (Figure 3.3):
#
#   [left ctrl pt] ──[ fwd ceil(n/2) arcs + DSMs ]──► [◇ match pt] ◄──[ bwd floor(n/2)+1 arcs ]── [right ctrl pt]
#
#   With n DSMs there are n+1 Keplerian arcs total.
#   N_fwd = ceil((n+1)/2)   arcs in the forward  half-phase
#   N_bwd = floor((n+1)/2)  arcs in the backward half-phase
#   N_dsm_fwd = N_fwd - 1   DSMs applied within the forward  half-phase
#   N_dsm_bwd = N_bwd - 1   DSMs applied within the backward half-phase
#
# Control variables:
#   Δv_k  (3-vector, km/s)  — impulsive maneuver at each mid-course node   [3n total]
#   α_k   (scalar ∈ [0,1])  — fractional time allocation for each arc      [n+1 total]
#     where  Δt_k = α_k · Δt_p  and  Δt_p = tf - t0
#   Subject to:  sum(α_k) = 1   (built-in equality constraint)
#
# Boundary variables (same as SimsFlanaganPhase):
#   v∞_dep  [3]   departure excess velocity
#   v∞_arr  [3]   arrival   excess velocity
#   t0, tf  [1]   phase start/end times
#   m0, mf  [1]   masses at left/right control points
#
# Match-point defect (7 equalities, lb = ub = 0):
#   c_mp = [X_bwd_mp - X_fwd_mp;  m_bwd_mp - m_fwd_mp]
#
# Alpha-sum constraint (1 equality, built-in):
#   sum(α_k, k=1..n+1) - 1 = 0
#
# Total built-in constraints = 8 (7 match-point + 1 alpha-sum)
#
# Notation conventions (matching the paper):
#   Φ_k  — 6×6 two-body STM for arc k
#   M_k  — 9×9 augmented MTM at DSM k     (paper §3.4, Eq. 3.49)
#   The augmented state is X = [r; v; m]  (9-vector, mass appended)
#   STM-MTM chains propagate 9×9 sensitivity matrices.
#
# TIME CONVENTION (departure from paper):
#   Phase times are ABSOLUTE (seconds from epoch), not relative to other phases.
#   The paper's Δt_previous terms are NOT included; the augmented state X
#   does not carry the extra time slots from Eq. (3.38).
#   Our augmented state is simply X = [r; v; m] (7-vector in physical units).
#
# DEPENDENCY: OptControlStubs.jl, kepler_propagator_time_domain.jl

# ─────────────────────────────────────────────────────────────────────────────
# MGAnDSMs  —  transcription descriptor  (parallel role to SimsFlanagan / LGL)
# ─────────────────────────────────────────────────────────────────────────────

"""
    MGAnDSMs(; n_dsm = 1, Isp = nothing, g0 = 9.80665e-3)

Transcription descriptor for a gravity-assist phase carrying `n_dsm` deep-space maneuvers, each an
impulse between two Keplerian arcs. A phase with `n` maneuvers has `n + 1` arcs, split between a
forward half-phase of `ceil((n+1)/2)` arcs and a backward half-phase of `floor((n+1)/2)`.

# Fields
- `n_dsm::Int`: number of deep-space maneuvers, which must be zero or more. Zero gives a single
  Keplerian arc between the control points.
- `ve::Float64`: exhaust velocity `Isp * g0`, in the phase's velocity units, and `NaN` when `Isp`
  was not given.

# Notes
`Isp` and `g0` belong to the propulsion system, so the usual place for them is the phase, as
`model = PropulsionModel(...)`. Setting them here also works and computes `ve` from the pair.

Follows Englander's thesis sections 3.2 and 3.4, with one departure: phase times are absolute
seconds from the epoch rather than relative to the previous phase, so the augmented state is
`[r; v; m]` without the paper's extra time slots.

# Example
```julia
MGAnDSMs(n_dsm = 2)
```
"""
struct MGAnDSMs
    n_dsm  ::Int    # number of deep-space maneuvers (≥ 0)
    ve     ::Float64  # exhaust velocity  ve = Isp * g0  (km/s)
end

# Isp and g0 belong to the propulsion system; pass them to the phase as
# `model = PropulsionModel(...)`. Setting them here still works.
function MGAnDSMs(; n_dsm::Int = 1, Isp = nothing, g0::Float64 = 9.80665e-3)
    n_dsm >= 0 || throw(ArgumentError("n_dsm must be ≥ 0 (got $n_dsm)"))
    MGAnDSMs(n_dsm, Isp === nothing ? NaN : Isp * g0)
end

# total keplerian arcs = n_dsm + 1, split across fwd/bwd half-phases
n_arcs_fwd(t::MGAnDSMs) = (t.n_dsm + 1 + 1) ÷ 2     # ceil((n+1)/2)
n_arcs_bwd(t::MGAnDSMs) = (t.n_dsm + 1)     ÷ 2     # floor((n+1)/2)
n_dsm_fwd(t::MGAnDSMs)  = n_arcs_fwd(t) - 1          # DSMs in fwd half
n_dsm_bwd(t::MGAnDSMs)  = n_arcs_bwd(t) - 1          # DSMs in bwd half
n_alphas(t::MGAnDSMs)   = t.n_dsm + 1                 # one α per arc

# ─────────────────────────────────────────────────────────────────────────────
# Variable type tags  (parallel to SFVInfinity3, SFThrustBlock, etc.)
# ─────────────────────────────────────────────────────────────────────────────

struct MGAVInfinity3     end   # v∞ departure or arrival:   3 scalars
struct MGADVBlock        end   # Δv maneuvers:  3 × n_dsm vectors (column-major)
struct MGAAlphaBlock     end   # fractional arc durations: (n_dsm+1) scalars
struct MGAMassParam      end   # m0 or mf: scalar
struct MGATime           end   # t0 or tf: scalar
struct MGAMatchPointBlock end  # sentinel in function_list
struct MGAAlphaSumBlock  end   # sentinel for the alpha-sum equality constraint

# ─────────────────────────────────────────────────────────────────────────────
# MGAnDSMsPhase <: AbstractShootingPhase
# ─────────────────────────────────────────────────────────────────────────────

"""
    MGAnDSMsPhase(; name, transcription, model, ephemeris_left, ephemeris_right)

A gravity-assist arc between two control points, transcribed by [`MGAnDSMs`](@ref). Each end follows
a body through an ephemeris callback, the Keplerian arcs and their deep-space maneuvers propagate
from both ends to a match point in the middle, and the phase carries eight built-in constraints: the
seven-component match-point disagreement and one equality holding the arc-duration fractions to one.

# Arguments
- `name::Symbol`: the phase's name, which labels it in reports.
- `transcription::MGAnDSMs`: the maneuver count.
- `model`: a `PropulsionModel` carrying `mu`, and `Isp` with `g0` when the maneuvers are to draw on
  a propulsion system. The three can be passed individually instead. `mu` is required; `Isp` given
  without `g0` throws, because the pair sets the units.
- `ephemeris_left`, `ephemeris_right`: callables `eph(t) -> (r, v, a)`, each a 3-vector, giving the
  body's position, velocity and gravitational acceleration at the control point.

# Fields
The struct's fields are the values above, with `model` decomposed into `mu`, `Isp` and `g0`, alongside
the solver variables and the cached decision values the framework fills during a solve. A caller sets
none of them directly; `Vary`, `Constraint` and `Objective` do.

# Notes
Units are the caller's and have only to agree with each other, with `mu` setting the length and time
units and `g0` following them, so a phase in kilometres and seconds takes `g0 = 9.80665e-3`.

The phase varies departure and arrival v-infinity, the maneuver impulses, the fractional arc
durations, the two times and the two masses, rather than a state and a control history.

# Example
```julia
phase = MGAnDSMsPhase(name = :earth_to_venus,
                      transcription = MGAnDSMs(n_dsm = 1),
                      model = PropulsionModel(mu = 1.32712440018e11,
                                              Isp = 3000.0,
                                              g0 = 9.80665e-3),
                      ephemeris_left = t -> ([1.5e8, 0.0, 0.0], [0.0, 29.8, 0.0], zeros(3)),
                      ephemeris_right = t -> ([1.1e8, 0.0, 0.0], [0.0, 35.0, 0.0], zeros(3)))
```
"""
mutable struct MGAnDSMsPhase <: AbstractShootingPhase
    name            ::Symbol
    transcription   ::MGAnDSMs

    # ── Gravitational parameter ────────────────────────────────────────────
    mu              ::Float64

    # Exhaust velocity is a property of the propulsion system, not of the
    # discretisation. It reached the transcription first and is still read from
    # there when a script sets it that way; `model = …` sets it here.
    Isp             ::Union{Float64,Nothing}
    g0              ::Union{Float64,Nothing}

    # ── Ephemeris callbacks ────────────────────────────────────────────────
    # Signature: eph(t::Float64) -> (r::Vector{Float64}(3), v::Vector{Float64}(3),
    #                                 a::Vector{Float64}(3))
    ephemeris_left  ::Any
    ephemeris_right ::Any

    # ── Decision variable blocks ───────────────────────────────────────────
    vinf_dep_var    ::Any   # DirectSolverVariable (MGAVInfinity3):  3
    vinf_arr_var    ::Any   # DirectSolverVariable (MGAVInfinity3):  3
    dv_var          ::Any   # DirectSolverVariable (MGADVBlock):     3 × n_dsm
    alpha_var       ::Any   # DirectSolverVariable (MGAAlphaBlock):  n_dsm+1
    t0_var          ::Any   # DirectSolverVariable (MGATime):        1
    tf_var          ::Any   # DirectSolverVariable (MGATime):        1
    m0_var          ::Any   # DirectSolverVariable (MGAMassParam):   1
    mf_var          ::Any   # DirectSolverVariable (MGAMassParam):   1

    # ── User-supplied boundary constraints and objective ───────────────────
    constraints         ::Vector{Any}
    objective           ::Any

    # ── Match-point constraint scaling (7-vector) ─────────────────────────
    matchpoint_scale    ::Vector{Float64}
    matchpoint_shift    ::Vector{Float64}

    # ── Cached decision variable values ───────────────────────────────────
    _vinf_dep   ::Vector{Float64}   # [3]
    _vinf_arr   ::Vector{Float64}   # [3]
    _dv         ::Matrix{Float64}   # [3, n_dsm]  columns are Δv_1 … Δv_n
    _alpha      ::Vector{Float64}   # [n_dsm+1]   α_1 … α_{n+1}, sum ≈ 1
    _t0         ::Float64
    _tf         ::Float64
    _m0         ::Float64
    _mf         ::Float64

    # ── Derived control-point states (built from ephemeris + v∞) ──────────
    _x0             ::Vector{Float64}   # [6]  left  ctrl pt spacecraft state
    _xf             ::Vector{Float64}   # [6]  right ctrl pt spacecraft state
    _planet_left    ::NTuple{3, Vector{Float64}}   # (r, v, a) at t0
    _planet_right   ::NTuple{3, Vector{Float64}}   # (r, v, a) at tf

    # ── Match-point propagated states (filled by evaluate_matchpoint!) ─────
    _x_match_fwd    ::Vector{Float64}   # [6]
    _x_match_bwd    ::Vector{Float64}   # [6]
    _m_match_fwd    ::Float64
    _m_match_bwd    ::Float64

    # ── Propagation history (diagnostics / warm-starting) ─────────────────
    # _states_fwd[:, k] = state at start of arc k (k = 1..N_fwd+1),  [6, N_fwd+1]
    # _masses_fwd[k]    = mass   at start of arc k
    _states_fwd     ::Matrix{Float64}   # [6, N_fwd+1]
    _masses_fwd     ::Vector{Float64}   # [N_fwd+1]
    _states_bwd     ::Matrix{Float64}   # [6, N_bwd+1]
    _masses_bwd     ::Vector{Float64}   # [N_bwd+1]

    # mu is a model, same as for Sims-Flanagan and collocation. The explicit
    # keyword still works so existing scripts do not move.
    function MGAnDSMsPhase(; name, transcription, model = nothing, mu = nothing,
                             Isp = nothing, g0 = nothing,
                             ephemeris_left = nothing, ephemeris_right = nothing)
        mu  = _model_field(model, :mu,  mu,  nothing)
        Isp = _model_field(model, :Isp, Isp, nothing)
        g0  = _model_field(model, :g0,  g0,  nothing)
        mu === nothing && throw(ArgumentError(
            "MGAnDSMsPhase needs mu — give `model = PropulsionModel(mu = ...)`."))
        Isp !== nothing && g0 === nothing && throw(ArgumentError(
            "MGAnDSMsPhase: Isp given without g0. Set both — their units must " *
            "match the units the phase flies in (km/s here, so g0 = 9.80665e-3)."))
        n  = transcription.n_dsm
        Nf = n_arcs_fwd(transcription)
        Nb = n_arcs_bwd(transcription)
        z3 = zeros(3)
        new(
            name, transcription, mu, Isp, g0,
            ephemeris_left, ephemeris_right,
            # decision variable slots (all nothing until set by caller)
            nothing, nothing, nothing, nothing,
            nothing, nothing, nothing, nothing,
            # constraints / objective
            Any[], nothing,
            # matchpoint scaling (default: no scaling)
            ones(7), zeros(7),
            # cached decision values
            zeros(3), zeros(3),
            zeros(3, max(n, 0)),          # _dv  [3 × n_dsm]
            ones(n + 1) ./ (n + 1),       # _alpha  uniform initial guess
            0.0, 0.0, 0.0, 0.0,
            # derived control point states
            zeros(6), zeros(6),
            (copy(z3), copy(z3), copy(z3)),
            (copy(z3), copy(z3), copy(z3)),
            # match-point states
            zeros(6), zeros(6), 0.0, 0.0,
            # propagation history
            zeros(6, Nf + 1), zeros(Nf + 1),
            zeros(6, Nb + 1), zeros(Nb + 1),
        )
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Accessors
# ─────────────────────────────────────────────────────────────────────────────

"""
    n_dsm(phase::MGAnDSMsPhase) -> Int

The number of deep-space maneuvers in the phase. The phase carries one more Keplerian arc than it
has maneuvers, since a maneuver separates two arcs.

# Returns
The maneuver count as an `Int`, zero or more.

# Example
<!-- doc-fragment -->
```julia
n_dsm(phase)
```
"""
n_dsm(p::MGAnDSMsPhase)       = p.transcription.n_dsm
n_arcs_fwd(p::MGAnDSMsPhase)  = n_arcs_fwd(p.transcription)
n_arcs_bwd(p::MGAnDSMsPhase)  = n_arcs_bwd(p.transcription)
n_dsm_fwd(p::MGAnDSMsPhase)   = n_dsm_fwd(p.transcription)
n_dsm_bwd(p::MGAnDSMsPhase)   = n_dsm_bwd(p.transcription)
n_alphas(p::MGAnDSMsPhase)     = n_alphas(p.transcription)
ve(p::MGAnDSMsPhase)           = p.Isp === nothing ? p.transcription.ve : p.Isp * p.g0

# ─────────────────────────────────────────────────────────────────────────────
# Configuration API  (mirrors SimsFlanaganPhase)
# ─────────────────────────────────────────────────────────────────────────────

set_departure_vinf!(p::MGAnDSMsPhase, var)  = (p.vinf_dep_var = var)
set_arrival_vinf!(p::MGAnDSMsPhase,   var)  = (p.vinf_arr_var = var)
set_ephemeris_left!(p::MGAnDSMsPhase,  eph)  = (p.ephemeris_left  = eph)
set_ephemeris_right!(p::MGAnDSMsPhase, eph)  = (p.ephemeris_right = eph)
set_dsm_control!(p::MGAnDSMsPhase,    var)  = (p.dv_var    = var)
set_alpha!(p::MGAnDSMsPhase,          var)  = (p.alpha_var = var)
set_initial_time!(p::MGAnDSMsPhase,   var)  = (p.t0_var    = var)
set_final_time!(p::MGAnDSMsPhase,     var)  = (p.tf_var    = var)
set_initial_mass!(p::MGAnDSMsPhase,   var)  = (p.m0_var    = var)
set_final_mass!(p::MGAnDSMsPhase,     var)  = (p.mf_var    = var)
set_objective!(p::MGAnDSMsPhase, obj::MayerObjective) = (p.objective = obj)
set_matchpoint_scale!(p::MGAnDSMsPhase, sc) = (p.matchpoint_scale = Float64.(sc); nothing)
set_matchpoint_shift!(p::MGAnDSMsPhase, sh) = (p.matchpoint_shift = Float64.(sh); nothing)

function add_constraint!(p::MGAnDSMsPhase, con)
    push!(p.constraints, con)
    return p
end

# ─────────────────────────────────────────────────────────────────────────────
# _propagate_half_mgandsms
#
# Propagates one half-phase of an MGAnDSMs transcription, applying impulsive
# Δv maneuvers between Keplerian arcs and updating mass via Tsiolkovsky.
#
# Arguments:
#   x_cp     : 6-vector, control point Cartesian state [r; v]  (km, km/s)
#   m_cp     : scalar, mass at control point (kg)
#   dv       : 3 × n_dsm_half  matrix of Δv vectors (km/s) for this half
#   alpha    : (n_arcs_half)-vector of fractional arc durations (sum ≤ 1)
#   dt_phase : scalar, total phase flight time = tf - t0  (s)
#   mu       : gravitational parameter (km³/s²)
#   ve       : exhaust velocity = Isp * g0  (km/s)
#   sign_dt  : +1 for forward propagation,  -1 for backward
#   sign_m   : -1 for forward (mass depletes), +1 for backward (mass accumulates)
#
# Returns a NamedTuple:
#   Phi      :: Vector{Matrix{Float64}}  — 6×6 STMs,  one per arc  [n_arcs_half]
#   dxdt     :: Vector{Vector{Float64}}  — 6-vec ∂x/∂(arc_time),  one per arc
#   m        :: Vector{Float64}          — mass at start  of each arc  [n_arcs_half+1]
#   x_nodes  :: Matrix{Float64}          — 6×(n_arcs_half+1) states at arc endpoints
# ─────────────────────────────────────────────────────────────────────────────

function _propagate_half_mgandsms(x_cp, m_cp, dv, alpha, dt_phase, mu, ve,
                                   sign_dt, sign_m)
    n_arcs = length(alpha)
    n_dv   = size(dv, 2)   # = n_arcs - 1

    T = promote_type(eltype(x_cp), typeof(m_cp), eltype(dv), eltype(alpha),
                     typeof(dt_phase))

    Phi    = [zeros(T, 6, 6) for _ in 1:n_arcs]
    dxdt   = [zeros(T, 6)    for _ in 1:n_arcs]
    m      = zeros(T, n_arcs + 1)
    x_nodes = zeros(T, 6, n_arcs + 1)

    x       = convert(Vector{T}, x_cp)
    m[1]    = m_cp
    x_nodes[:, 1] = x

    for k in 1:n_arcs
        dt_k = alpha[k] * dt_phase   # arc duration (signed via sign_dt below)

        # Propagate Keplerian arc k
        res = kepler_propagate_time_domain(x, mu, sign_dt * dt_k;
                  need_stm = true, need_time_partials = true)
        Phi[k]  = res.stm
        dxdt[k] = res.dstate_dt   # ∂x_end/∂(arc_time), i.e. ∂x/∂(sign_dt·dt_k)
        x       = res.state

        # Apply DSM (if this is not the last arc in the half-phase)
        if k < n_arcs
            dv_k      = dv[:, k]
            dv_mag    = norm(dv_k)
            # Tsiolkovsky — forward (sign_m=-1): depletes  m_{k+1} = m_k·exp(-|Δv|/ve)
            #              backward (sign_m=+1): accumulates m_{k+1} = m_k·exp(+|Δv|/ve)
            m[k+1]    = m[k] * exp(sign_m * dv_mag / ve)
            x[4:6] .+= sign_dt * dv_k            # Δv applied in direction of propagation
        else
            m[k+1]    = m[k]   # no DSM after last arc
        end
        x_nodes[:, k+1] = x
    end

    return (Phi = Phi, dxdt = dxdt, m = m, x_nodes = x_nodes)
end

# ─────────────────────────────────────────────────────────────────────────────
# evaluate_matchpoint!
#
# Propagates both half-phases from their control points to the match point
# and caches results.  Called before get_functions and matchpoint_jacobian.
# ─────────────────────────────────────────────────────────────────────────────

function evaluate_matchpoint!(p::MGAnDSMsPhase)
    n    = n_dsm(p)
    Nf   = n_arcs_fwd(p)
    Nb   = n_arcs_bwd(p)
    ndf  = n_dsm_fwd(p)
    ndb  = n_dsm_bwd(p)
    dt_p = p._tf - p._t0

    # Build control-point spacecraft states from ephemeris + v∞
    if p.ephemeris_left !== nothing
        r0, v0, a0    = p.ephemeris_left(p._t0)
        p._planet_left = (r0, v0, a0)
        p._x0 = vcat(r0, v0 .+ p._vinf_dep)
    end
    if p.ephemeris_right !== nothing
        rf, vf, af    = p.ephemeris_right(p._tf)
        p._planet_right = (rf, vf, af)
        p._xf = vcat(rf, vf .+ p._vinf_arr)
    end

    # Split alphas and Δvs between the two half-phases
    alpha_fwd = p._alpha[1:Nf]            # arcs 1 … N_fwd
    alpha_bwd = p._alpha[Nf+1:end]        # arcs N_fwd+1 … N_fwd+N_bwd

    dv_fwd = ndf > 0 ? p._dv[:, 1:ndf]                   : zeros(3, 0)
    dv_bwd = ndb > 0 ? p._dv[:, ndf+1:ndf+ndb]           : zeros(3, 0)

    Fwd = _propagate_half_mgandsms(p._x0, p._m0, dv_fwd, alpha_fwd, dt_p,
                                    p.mu, ve(p), +1.0, -1.0)
    Bwd = _propagate_half_mgandsms(p._xf, p._mf, dv_bwd, alpha_bwd, dt_p,
                                    p.mu, ve(p), -1.0, +1.0)

    # The match point is the end of the forward propagation (= end of backward)
    p._x_match_fwd = Fwd.x_nodes[:, Nf+1]
    p._m_match_fwd = Fwd.m[Nf+1]
    p._x_match_bwd = Bwd.x_nodes[:, Nb+1]
    p._m_match_bwd = Bwd.m[Nb+1]

    p._states_fwd  = Fwd.x_nodes
    p._masses_fwd  = Fwd.m
    p._states_bwd  = Bwd.x_nodes
    p._masses_bwd  = Bwd.m

    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# matchpoint_defect  →  7-vector  [X_bwd_mp - X_fwd_mp; m_bwd_mp - m_fwd_mp]
# ─────────────────────────────────────────────────────────────────────────────

function matchpoint_defect(p::MGAnDSMsPhase)
    vcat(p._x_match_bwd .- p._x_match_fwd,
         p._m_match_bwd  - p._m_match_fwd)
end

# ─────────────────────────────────────────────────────────────────────────────
# alpha_sum_defect  →  scalar  sum(α) - 1  (= 0 at feasibility)
# ─────────────────────────────────────────────────────────────────────────────

alpha_sum_defect(p::MGAnDSMsPhase) = sum(p._alpha) - 1.0

# ─────────────────────────────────────────────────────────────────────────────
# matchpoint_jacobian
#
# Analytic Jacobian of the 7-constraint match-point defect with respect to all
# decision variable blocks.  Implements the STM-MTM chain method from
# Englander thesis §3.4 specialised to the MGAnDSMs transcription.
#
# STM (Φ_k): 6×6 two-body state transition matrix for arc k (from propagator).
# MTM (M_k): maps 7-vector [r; v; m]⁻ before DSM k to [r; v; m]⁺ after DSM k.
#   From Eq. (3.49):
#     M_k = [I₃  0₃  0₃₁]   r is unchanged: r⁺ = r⁻
#           [0₃  I₃  0₃₁]   v⁺ = v⁻ + Δv_k  (identity for v, Δv handled separately)
#           [0₁₃ 0₁₃ M33]   m⁺ = m⁻ · exp(-|Δv|/ve)
#   where M33 = exp(-|Δv_k|/ve) = m_{k+1}/m_k  (Eq. 3.50)
#
# The match-point Jacobian ∂c_mp/∂x = ∂X_bwd/∂x - ∂X_fwd/∂x  (Eq. 3.36).
# The chain propagates 7×7 augmented sensitivity matrices.
#
# Returns a NamedTuple with fields:
#   dvinf_dep  :: 7×3  matrix
#   dvinf_arr  :: 7×3  matrix
#   ddv        :: 7 × 3·n_dsm  matrix  (columns: dv1x,dv1y,dv1z, dv2x,...)
#   dalpha     :: 7 × (n_dsm+1)  matrix
#   dm0        :: 7-vector
#   dmf        :: 7-vector
#   dt0        :: 7-vector
#   dtf        :: 7-vector
#   dx0        :: 7×6  (internal; for testing)
#   dxf        :: 7×6  (internal; for testing)
# ─────────────────────────────────────────────────────────────────────────────

function matchpoint_jacobian(p::MGAnDSMsPhase)
    n    = n_dsm(p)
    Nf   = n_arcs_fwd(p)
    Nb   = n_arcs_bwd(p)
    ndf  = n_dsm_fwd(p)
    ndb  = n_dsm_bwd(p)
    dt_p = p._tf - p._t0
    ve_p = ve(p)

    # ── Re-propagate to collect STMs ─────────────────────────────────────────
    alpha_fwd = p._alpha[1:Nf]
    alpha_bwd = p._alpha[Nf+1:end]
    dv_fwd = ndf > 0 ? p._dv[:, 1:ndf]         : zeros(3, 0)
    dv_bwd = ndb > 0 ? p._dv[:, ndf+1:ndf+ndb] : zeros(3, 0)

    Fwd = _propagate_half_mgandsms(p._x0, p._m0, dv_fwd, alpha_fwd, dt_p,
                                    p.mu, ve_p, +1.0, -1.0)
    Bwd = _propagate_half_mgandsms(p._xf, p._mf, dv_bwd, alpha_bwd, dt_p,
                                    p.mu, ve_p, -1.0, +1.0)

    I7 = Matrix{Float64}(I, 7, 7)
    I6 = Matrix{Float64}(I, 6, 6)

    # ── Build augmented 7×7 STMs: Φ̄_k = blockdiag(Φ_k, 1) ──────────────────
    # Mass is unchanged across a Keplerian arc.
    function aug_stm(Phi6)
        A = zeros(7, 7)
        A[1:6, 1:6] = Phi6
        A[7, 7]     = 1.0
        return A
    end

    # ── Build MTMs at each DSM ─────────────────────────────────────────────
    # M_k (7×7) for a single impulsive DSM:  r unchanged, v+Δv, m scaled.
    # ∂m⁺/∂m⁻ = exp(sign_m·|Δv|/ve) = M33
    #   sign_m = -1  (forward):  M33 = exp(-|Δv|/ve) < 1  (depletion)
    #   sign_m = +1  (backward): M33 = exp(+|Δv|/ve) > 1  (accumulation)
    function build_mtm(dv_k, sign_m)
        dv_mag = norm(dv_k)
        M33    = exp(sign_m * dv_mag / ve_p)   # = m_{k+1}/m_k
        M      = Matrix{Float64}(I, 7, 7)
        M[7, 7] = M33
        return M
    end

    # Forward half: n_arcs_fwd Φ̄s and n_dsm_fwd MTMs
    Phi_aug_fwd = [aug_stm(Fwd.Phi[k]) for k in 1:Nf]
    Mtm_fwd     = [build_mtm(dv_fwd[:, k], -1.0) for k in 1:ndf]

    # Backward half: n_arcs_bwd Φ̄s and n_dsm_bwd MTMs
    Phi_aug_bwd = [aug_stm(Bwd.Phi[k]) for k in 1:Nb]
    Mtm_bwd     = [build_mtm(dv_bwd[:, k], +1.0) for k in 1:ndb]

    # ── STM-MTM chains (Eq. 3.51): right-to-left sweep ───────────────────────
    #
    # For the forward half with N_fwd arcs and N_dsm_fwd = N_fwd-1 DSMs:
    #   Segment sequence: Φ̄₁, M₁, Φ̄₂, M₂, …, M_{Nf-1}, Φ̄_{Nf}
    #
    # chain_fwd[k] = product from end of arc k to match point
    #   chain_fwd[Nf] = I  (at match point)
    #   chain_fwd[k]  = chain_fwd[k+1] · Φ̄_{k+1}  if k < Nf and arc k+1 is the next
    #   But MTMs sit between arcs: after arc k comes DSM k, then arc k+1.
    #   So the product from just-before-DSM-k to match point is:
    #     chain_before_dsm_fwd[k] = chain_fwd[k+1] · Φ̄_{k+1} · ... (rest of fwd chain)
    #
    # We track two chains:
    #   tail_arc_fwd[k]:  7×7  sensitivity of match-pt state w.r.t. state at END of arc k
    #   tail_arc_fwd[Nf] = I
    #   tail_arc_fwd[k]  = (DSM k present?) M_k · Φ̄_{k+1} · tail_arc_fwd[k+1]  for k < Nf

    tail_arc_fwd = Vector{Matrix{Float64}}(undef, Nf)
    tail_arc_fwd[Nf] = I7
    for k in (Nf-1):-1:1
        # after arc k: apply DSM k, then arc k+1, then rest
        tail_arc_fwd[k] = tail_arc_fwd[k+1] * Phi_aug_fwd[k+1] * Mtm_fwd[k]
    end

    tail_arc_bwd = Vector{Matrix{Float64}}(undef, Nb)
    if Nb > 0
        tail_arc_bwd[Nb] = I7
        for k in (Nb-1):-1:1
            tail_arc_bwd[k] = tail_arc_bwd[k+1] * Phi_aug_bwd[k+1] * Mtm_bwd[k]
        end
    end

    # Cumulative STM-MTM chain from control point to match point
    # cum_fwd = tail_arc_fwd[1] · Φ̄_1  (7×7)
    cum_fwd = (Nf > 0) ? tail_arc_fwd[1] * Phi_aug_fwd[1] : I7
    cum_bwd = (Nb > 0) ? tail_arc_bwd[1] * Phi_aug_bwd[1] : I7

    # ─────────────────────────────────────────────────────────────────────────
    # 1. ∂c_mp / ∂x0  (7×6) and ∂c_mp / ∂xf  (7×6)
    #    Mass is the 7th component; ∂m_fwd_mp/∂x0 = 0 (mass depends only on m0 and Δvs)
    #    ∂m_bwd_mp/∂xf = 0 similarly.
    #    c_mp = X_bwd - X_fwd → ∂c/∂x0 = -∂X_fwd/∂x0,  ∂c/∂xf = +∂X_bwd/∂xf
    # ─────────────────────────────────────────────────────────────────────────
    J_x0 = -cum_fwd[:, 1:6]   # 7×6:  extract columns for r,v (mass col is zero)
    J_xf = +cum_bwd[:, 1:6]   # 7×6

    # ─────────────────────────────────────────────────────────────────────────
    # 2. ∂c_mp / ∂m0  (7-vector) and ∂c_mp / ∂mf  (7-vector)
    #    cum_fwd[:, 7] is ∂X_fwd_mp/∂m0  (via 7th column of the chain)
    # ─────────────────────────────────────────────────────────────────────────
    J_dm0 = -cum_fwd[:, 7]    # 7-vec; c_mp = X_bwd - X_fwd → minus sign
    J_dmf = +cum_bwd[:, 7]    # 7-vec

    # ─────────────────────────────────────────────────────────────────────────
    # 3. ∂c_mp / ∂Δv_k  (7×3 per DSM, assembled into 7 × 3·n_dsm)
    #
    # For a DSM at position k in the forward half (k = 1 … n_dsm_fwd):
    #   The impulse adds [0₃; Δv_k; 0] to the augmented state immediately after
    #   arc k.  So the sensitivity connection matrix (Ξ_k, Eq. 3.62) is:
    #     Ξ_k = ∂X_k⁺/∂Δv_k = [0₃ₓ₃; I₃ₓ₃; ∂m_{k+1}/∂Δv_k]  (7×3)
    #   where ∂m_{k+1}/∂Δv_k = -m_k/ve · Δv_k/|Δv_k| · exp(-|Δv_k|/ve)
    #                         = -(m_{k+1}/ve) · Δv_k/|Δv_k|  (Eq. 3.61)
    #
    #   The match-point sensitivity:
    #     ∂X_fwd_mp/∂Δv_k = tail_arc_fwd[k+1] · Φ̄_{k+1} · Ξ_k   (Eq. 3.52)
    #     (tail chain starts after the DSM, so it sees the post-DSM arc k+1)
    # ─────────────────────────────────────────────────────────────────────────
    J_ddv = zeros(7, 3 * n)

    B = vcat(zeros(3, 3), I(3))   # 6×3: velocity rows

    # Forward DSMs
    for k in 1:ndf
        dv_k   = dv_fwd[:, k]
        dv_mag = norm(dv_k)
        m_k    = Fwd.m[k+1]   # post-burn mass
        # ∂m_{k+1}/∂Δv_k (1×3)
        dm_dv  = -(m_k / ve_p) .* dv_k ./ (dv_mag + 1e-20)   # 3-vec
        # Ξ_k: 7×3
        Xi_k   = vcat(zeros(3, 3), I(3), dm_dv')
        # chain after the DSM: tail_arc_fwd[k+1] · Φ̄_{k+1}
        chain_after = tail_arc_fwd[k+1] * Phi_aug_fwd[k+1]
        dXfwd_dv = chain_after * Xi_k   # 7×3
        col = (k-1)*3+1 : k*3
        J_ddv[:, col] .= -dXfwd_dv   # c_mp = X_bwd - X_fwd
    end

    # Backward DSMs (DSM indices ndf+1 … n in the global _dv matrix)
    # Backward propagator applies sign_dt * dv_k = -dv_k, so ∂v_post/∂dv_k = -I.
    # Mass: sign_m = +1 → ∂m_after/∂dv_k = +(m_after/ve)·dv/|dv|
    for k in 1:ndb
        dv_k   = dv_bwd[:, k]
        dv_mag = norm(dv_k)
        m_k    = Bwd.m[k+1]   # post-burn mass (in backward propagation)
        dm_dv  = +(m_k / ve_p) .* dv_k ./ (dv_mag + 1e-20)
        Xi_k   = vcat(zeros(3, 3), -I(3), dm_dv')   # -I: backward sign_dt = -1
        chain_after = tail_arc_bwd[k+1] * Phi_aug_bwd[k+1]
        dXbwd_dv = chain_after * Xi_k
        col = (ndf + k - 1)*3+1 : (ndf + k)*3
        J_ddv[:, col] .= +dXbwd_dv   # c_mp = X_bwd - X_fwd → + sign for bwd
    end

    # ─────────────────────────────────────────────────────────────────────────
    # 4. ∂c_mp / ∂α_k  (7 per alpha, assembled into 7 × (n+1))
    #
    # α_k controls the duration of arc k: Δt_k = α_k · Δt_p.
    # The propagation time partial from kepler_propagate_time_domain is
    #   dxdt[k] = ∂x_end/∂(arc_time)  where arc_time = sign_dt · Δt_k.
    # So ∂x_end/∂α_k = dxdt[k] · (sign_dt · Δt_p).
    #
    # For α_k in the forward half (k = 1 … N_fwd):
    #   Ξ_k = [dxdt_fwd[k] · Δt_p; 0]  (7-vector: no mass change from time)
    #   ∂X_fwd_mp/∂α_k = (chain from end of arc k to mp) · Ξ_k
    #                  = tail_arc_fwd[k] · Ξ_k
    #
    # For α_k in the backward half (k = N_fwd+1 … N_fwd+N_bwd):
    #   (backward arc index is kb = k - N_fwd)
    #   Ξ_k = [dxdt_bwd[kb] · (-1) · Δt_p; 0]  (sign_dt = -1 for backward)
    #   ∂X_bwd_mp/∂α_k = tail_arc_bwd[kb] · Ξ_k
    # ─────────────────────────────────────────────────────────────────────────
    J_dal = zeros(7, n + 1)

    # Forward alphas (α_1 … α_{N_fwd})
    for k in 1:Nf
        dxdt_k = vcat(Fwd.dxdt[k], 0.0)   # augment: mass row zero
        Xi_k   = dxdt_k .* dt_p            # 7-vec: ∂X_arc_end/∂α_k (sign_dt=+1)
        dXfwd_dal = tail_arc_fwd[k] * Xi_k  # 7-vec
        J_dal[:, k] .= -dXfwd_dal          # c_mp = X_bwd - X_fwd
    end

    # Backward alphas (α_{N_fwd+1} … α_{N_fwd+N_bwd})
    for kb in 1:Nb
        k      = Nf + kb
        dxdt_k = vcat(Bwd.dxdt[kb], 0.0)
        Xi_k   = dxdt_k .* (-dt_p)         # sign_dt = -1 for backward
        dXbwd_dal = tail_arc_bwd[kb] * Xi_k
        J_dal[:, k] .= +dXbwd_dal          # c_mp uses +∂X_bwd
    end

    # ─────────────────────────────────────────────────────────────────────────
    # 5. ∂c_mp / ∂t0  and  ∂c_mp / ∂tf  (7-vectors each)
    #
    # Two contributions:
    #   (a) Ephemeris: ∂x0/∂t0 = [v_planet(t0); a_planet(t0)]
    #       chain-rules through cum_fwd (same as SF)
    #   (b) Phase duration: Δt_p = tf - t0, so ∂Δt_k/∂tf = +α_k, ∂Δt_k/∂t0 = -α_k
    #       This means every arc time changes → sum over forward arcs:
    #         ∂X_fwd_mp/∂tf|duration = Σ_k tail_arc_fwd[k] · [dxdt_fwd[k]·α_k; 0]
    #       and backward arcs (sign_dt = -1):
    #         ∂X_bwd_mp/∂tf|duration = Σ_k tail_arc_bwd[kb] · [dxdt_bwd[kb]·(-α_{Nf+kb}); 0]
    # ─────────────────────────────────────────────────────────────────────────
    dtf_fwd_dur = zeros(7)
    dt0_fwd_dur = zeros(7)
    for k in 1:Nf
        dxdt_k = vcat(Fwd.dxdt[k], 0.0)
        contrib = tail_arc_fwd[k] * (dxdt_k .* alpha_fwd[k])
        dtf_fwd_dur .+= contrib     # ∂Δt_k/∂tf = +α_k
        dt0_fwd_dur .-= contrib     # ∂Δt_k/∂t0 = -α_k
    end

    dtf_bwd_dur = zeros(7)
    dt0_bwd_dur = zeros(7)
    for kb in 1:Nb
        dxdt_k = vcat(Bwd.dxdt[kb], 0.0)
        # backward sign_dt = -1, so ∂x_end/∂dt_p = dxdt[kb] · (-1) · α_{Nf+kb}
        contrib = tail_arc_bwd[kb] * (dxdt_k .* (-alpha_bwd[kb]))
        dtf_bwd_dur .+= contrib     # ∂Δt_p/∂tf = +1
        dt0_bwd_dur .-= contrib     # ∂Δt_p/∂t0 = -1
    end

    # c_mp = X_bwd - X_fwd → combine signs
    J_dtf_kepler = dtf_bwd_dur .- dtf_fwd_dur
    J_dt0_kepler = dt0_bwd_dur .- dt0_fwd_dur

    # Ephemeris contributions
    J_dt0 = copy(J_dt0_kepler)
    J_dtf = copy(J_dtf_kepler)
    if p.ephemeris_left !== nothing
        v0_p, a0_p = p._planet_left[2], p._planet_left[3]
        dX0_dt0 = vcat(v0_p, a0_p, 0.0)      # ∂x0/∂t0 augmented (mass row zero)
        J_dt0 .+= J_x0 * vcat(v0_p, a0_p)    # J_x0 is 7×6, vcat gives 6-vec
    end
    if p.ephemeris_right !== nothing
        vf_p, af_p = p._planet_right[2], p._planet_right[3]
        J_dtf .+= J_xf * vcat(vf_p, af_p)
    end

    # ── v∞ Jacobians: chain-rule ∂x0/∂v∞ = [0; I; 0] restriction ────────────
    J_vinf_dep = J_x0[:, 4:6]   # 7×3: only velocity rows matter
    J_vinf_arr = J_xf[:, 4:6]   # 7×3

    return (dvinf_dep = J_vinf_dep,
            dvinf_arr = J_vinf_arr,
            ddv       = J_ddv,
            dalpha    = J_dal,
            dm0       = J_dm0,
            dmf       = J_dmf,
            dt0       = J_dt0,
            dtf       = J_dtf,
            dx0       = J_x0,
            dxf       = J_xf)
end

# ─────────────────────────────────────────────────────────────────────────────
# alpha_sum_jacobian
#
# Jacobian of the alpha-sum constraint (scalar) w.r.t. decision variables.
# Only the alpha_var block is non-zero: ∂(Σα - 1)/∂α_k = 1 for all k.
# Returns a NamedTuple parallel to matchpoint_jacobian.
# ─────────────────────────────────────────────────────────────────────────────

function alpha_sum_jacobian(p::MGAnDSMsPhase)
    na = n_alphas(p)
    # All other blocks are zero; only dalpha is non-zero.
    # Returns 1×(n_dsm+1) matrix (1 constraint row).
    dal = ones(1, na)
    return (dalpha = dal,)
end

# ═════════════════════════════════════════════════════════════════════════════
# NLP Interface
# ═════════════════════════════════════════════════════════════════════════════

# ─────────────────────────────────────────────────────────────────────────────
# variable_list
# ─────────────────────────────────────────────────────────────────────────────

function variable_list(p::MGAnDSMsPhase)
    filter(!isnothing, Any[
        p.vinf_dep_var, p.vinf_arr_var,
        p.dv_var, p.alpha_var,
        p.t0_var, p.tf_var,
        p.m0_var, p.mf_var,
    ])
end

# ─────────────────────────────────────────────────────────────────────────────
# nlp_length
# ─────────────────────────────────────────────────────────────────────────────

function nlp_length(p::MGAnDSMsPhase, v::DirectSolverVariable)
    vt = v.var
    vt isa MGAVInfinity3   && return 3
    vt isa MGADVBlock      && return 3 * n_dsm(p)
    vt isa MGAAlphaBlock   && return n_alphas(p)
    vt isa MGAMassParam    && return 1
    (vt isa MGATime || vt isa AbstractTime) && return 1
    throw(ArgumentError(
        "an MGAnDSMsPhase variable must be a state, burn, mass, time or v-infinity " *
        "variable; got $(typeof(vt))"))
end

nlp_length(p::MGAnDSMsPhase) = sum(nlp_length(p, v) for v in variable_list(p))

# ─────────────────────────────────────────────────────────────────────────────
# _var_nlp_scale / _var_nlp_shift  (tile for block variables)
# ─────────────────────────────────────────────────────────────────────────────

function _var_nlp_scale(p::MGAnDSMsPhase, v::DirectSolverVariable)
    v.var isa MGADVBlock    && return repeat(v.scale, n_dsm(p))
    return v.scale
end

function _var_nlp_shift(p::MGAnDSMsPhase, v::DirectSolverVariable)
    v.var isa MGADVBlock    && return repeat(v.shift, n_dsm(p))
    return v.shift
end

# ─────────────────────────────────────────────────────────────────────────────
# nlp_bounds
# ─────────────────────────────────────────────────────────────────────────────

function nlp_bounds(p::MGAnDSMsPhase, v::DirectSolverVariable)
    lb = v.var isa MGADVBlock ? repeat(v.lower_bounds, n_dsm(p)) : copy(v.lower_bounds)
    ub = v.var isa MGADVBlock ? repeat(v.upper_bounds, n_dsm(p)) : copy(v.upper_bounds)
    sc = _var_nlp_scale(p, v)
    sh = _var_nlp_shift(p, v)
    return (lb .- sh) ./ sc, (ub .- sh) ./ sc
end

# ─────────────────────────────────────────────────────────────────────────────
# get_variable_bounds
# ─────────────────────────────────────────────────────────────────────────────

function get_variable_bounds(p::MGAnDSMsPhase)
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
# variable_ranges
# ─────────────────────────────────────────────────────────────────────────────

function variable_ranges(p::MGAnDSMsPhase)
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
# get_decision_vector
# ─────────────────────────────────────────────────────────────────────────────

function get_decision_vector(p::MGAnDSMsPhase)
    n  = n_dsm(p)
    na = n_alphas(p)
    x  = zeros(nlp_length(p))
    off = 0
    if p.vinf_dep_var !== nothing
        v = p.vinf_dep_var
        x[off+1:off+3] .= (v.value .- v.shift) ./ v.scale;         off += 3
    end
    if p.vinf_arr_var !== nothing
        v = p.vinf_arr_var
        x[off+1:off+3] .= (v.value .- v.shift) ./ v.scale;         off += 3
    end
    if p.dv_var !== nothing && n > 0
        v  = p.dv_var
        sc = repeat(v.scale, n); sh = repeat(v.shift, n)
        x[off+1:off+3n] .= (v.value .- sh) ./ sc;                  off += 3n
    end
    if p.alpha_var !== nothing
        v  = p.alpha_var
        x[off+1:off+na] .= (v.value .- v.shift) ./ v.scale;        off += na
    end
    if p.t0_var !== nothing
        v = p.t0_var
        x[off+1] = (v.value[1] - v.shift[1]) / v.scale[1];         off += 1
    end
    if p.tf_var !== nothing
        v = p.tf_var
        x[off+1] = (v.value[1] - v.shift[1]) / v.scale[1];         off += 1
    end
    if p.m0_var !== nothing
        v = p.m0_var
        x[off+1] = (v.value[1] - v.shift[1]) / v.scale[1];         off += 1
    end
    if p.mf_var !== nothing
        v = p.mf_var
        x[off+1] = (v.value[1] - v.shift[1]) / v.scale[1];         off += 1
    end
    return x
end

# ─────────────────────────────────────────────────────────────────────────────
# set_decision_vector!
# ─────────────────────────────────────────────────────────────────────────────

function set_decision_vector!(p::MGAnDSMsPhase, x::Vector{Float64})
    n  = n_dsm(p)
    na = n_alphas(p)
    off = 0
    if p.vinf_dep_var !== nothing
        v = p.vinf_dep_var
        v.value = x[off+1:off+3] .* v.scale .+ v.shift
        p._vinf_dep = v.value;                                       off += 3
    end
    if p.vinf_arr_var !== nothing
        v = p.vinf_arr_var
        v.value = x[off+1:off+3] .* v.scale .+ v.shift
        p._vinf_arr = v.value;                                       off += 3
    end
    if p.dv_var !== nothing && n > 0
        v  = p.dv_var
        sc = repeat(v.scale, n); sh = repeat(v.shift, n)
        v.value = x[off+1:off+3n] .* sc .+ sh
        p._dv = reshape(v.value, 3, n);                              off += 3n
    end
    if p.alpha_var !== nothing
        v = p.alpha_var
        v.value = x[off+1:off+na] .* v.scale .+ v.shift
        p._alpha = v.value;                                          off += na
    end
    if p.t0_var !== nothing
        v = p.t0_var
        v.value = [x[off+1] * v.scale[1] + v.shift[1]]
        p._t0 = v.value[1];                                          off += 1
    end
    if p.tf_var !== nothing
        v = p.tf_var
        v.value = [x[off+1] * v.scale[1] + v.shift[1]]
        p._tf = v.value[1];                                          off += 1
    end
    if p.m0_var !== nothing
        v = p.m0_var
        v.value = [x[off+1] * v.scale[1] + v.shift[1]]
        p._m0 = v.value[1];                                          off += 1
    end
    if p.mf_var !== nothing
        v = p.mf_var
        v.value = [x[off+1] * v.scale[1] + v.shift[1]]
        p._mf = v.value[1];                                          off += 1
    end
    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# MGABoundaryContext  —  value object passed to boundary constraint closures
# ─────────────────────────────────────────────────────────────────────────────

struct MGABoundaryContext{T<:Real}
    x0       ::Vector{T}   # [6]  left  ctrl pt spacecraft state
    xf       ::Vector{T}   # [6]  right ctrl pt spacecraft state
    vinf_dep ::Vector{T}   # [3]  departure v∞
    vinf_arr ::Vector{T}   # [3]  arrival   v∞
    t0       ::T
    tf       ::T
    m0       ::T
    mf       ::T
    dv       ::Matrix{T}   # [3, n_dsm]  all DSM vectors
    alpha    ::Vector{T}   # [n_dsm+1]
end

# Upcast constructor: MGABoundaryContext{T}(ctx::MGABoundaryContext) promotes all
# fields to type T.  Required so ForwardDiff Dual numbers flow through sequence
# constraint AD Jacobians.
function MGABoundaryContext{T}(c::MGABoundaryContext) where {T<:Real}
    MGABoundaryContext{T}(
        T.(c.x0), T.(c.xf),
        T.(c.vinf_dep), T.(c.vinf_arr),
        T(c.t0), T(c.tf), T(c.m0), T(c.mf),
        T.(c.dv), T.(c.alpha),
    )
end

# _mga_context_replace: returns MGABoundaryContext{T} with one variable's physical
# value replaced by phys (a Vector{T}).  Used by the AD fallback for
# SequenceConstraint Jacobians.
function _mga_context_replace(ctx::MGABoundaryContext,
                               p::MGAnDSMsPhase,
                               var::DirectSolverVariable,
                               phys::AbstractVector{T}) where {T<:Real}
    oid = objectid(var)
    x0  = T.(ctx.x0);   xf  = T.(ctx.xf)
    vd  = T.(ctx.vinf_dep);  va  = T.(ctx.vinf_arr)
    t0  = T(ctx.t0);    tf  = T(ctx.tf)
    m0  = T(ctx.m0);    mf  = T(ctx.mf)
    dv  = T.(ctx.dv);   al  = T.(ctx.alpha)
    !isnothing(p.vinf_dep_var) && oid == objectid(p.vinf_dep_var) && (vd = phys)
    !isnothing(p.vinf_arr_var) && oid == objectid(p.vinf_arr_var) && (va = phys)
    !isnothing(p.t0_var)       && oid == objectid(p.t0_var)       && (t0 = phys[1])
    !isnothing(p.tf_var)       && oid == objectid(p.tf_var)       && (tf = phys[1])
    !isnothing(p.m0_var)       && oid == objectid(p.m0_var)       && (m0 = phys[1])
    !isnothing(p.mf_var)       && oid == objectid(p.mf_var)       && (mf = phys[1])
    !isnothing(p.dv_var)       && oid == objectid(p.dv_var)       && (dv = reshape(phys, 3, :))
    !isnothing(p.alpha_var)    && oid == objectid(p.alpha_var)    && (al = phys)
    MGABoundaryContext{T}(x0, xf, vd, va, t0, tf, m0, mf, dv, al)
end

# Generic sequence-constraint dispatch entry (called from ShootingManager).
# Delegates to the MGA-specific builder.
_shoot_boundary_context(p::MGAnDSMsPhase) = _mga_boundary_context(p)

function _mga_boundary_context(p::MGAnDSMsPhase)
    vinf_dep = p.vinf_dep_var !== nothing ? p.vinf_dep_var.value : p._vinf_dep
    vinf_arr = p.vinf_arr_var !== nothing ? p.vinf_arr_var.value : p._vinf_arr
    t0  = p.t0_var    !== nothing ? p.t0_var.value[1]    : p._t0
    tf  = p.tf_var    !== nothing ? p.tf_var.value[1]    : p._tf
    m0  = p.m0_var    !== nothing ? p.m0_var.value[1]    : p._m0
    mf  = p.mf_var    !== nothing ? p.mf_var.value[1]    : p._mf
    dv  = p.dv_var    !== nothing ? reshape(p.dv_var.value, 3, n_dsm(p)) : p._dv
    al  = p.alpha_var !== nothing ? p.alpha_var.value     : p._alpha
    MGABoundaryContext(copy(p._x0), copy(p._xf),
                       copy(vinf_dep), copy(vinf_arr),
                       t0, tf, m0, mf, copy(dv), copy(al))
end

# ─────────────────────────────────────────────────────────────────────────────
# function_list
# ─────────────────────────────────────────────────────────────────────────────

function function_list(p::MGAnDSMsPhase)
    result = PhaseFunction[
        PhaseFunction(MGAMatchPointBlock(), 7, zeros(7), zeros(7), "matchpoint_defect"),
    ]
    # Emit the α-sum equality only when α is an actual decision variable;
    # otherwise the row is identically zero with a zero Jacobian, which
    # makes the constraint Jacobian rank-deficient and can trip IPOPT's
    # "local infeasibility" detection.
    if p.alpha_var !== nothing
        push!(result, PhaseFunction(MGAAlphaSumBlock(), 1, [0.0], [0.0], "alpha_sum"))
    end
    for c in p.constraints
        if c isa BoundaryConstraint
            n_c = length(c.lower_bounds)
            push!(result, PhaseFunction(c, n_c, c.lower_bounds, c.upper_bounds,
                                        c.calc.name))
        end
    end
    return result
end

# ─────────────────────────────────────────────────────────────────────────────
# get_functions
# ─────────────────────────────────────────────────────────────────────────────

function get_functions(p::MGAnDSMsPhase)
    evaluate_matchpoint!(p)
    F = (matchpoint_defect(p) .- p.matchpoint_shift) ./ p.matchpoint_scale

    # Alpha-sum built-in constraint (no scaling needed) — only when α is
    # a registered decision variable; otherwise the row is omitted from
    # function_list and must not appear here either.
    if p.alpha_var !== nothing
        push!(F, alpha_sum_defect(p))
    end

    for c in p.constraints
        if c isa BoundaryConstraint
            ctx  = _mga_boundary_context(p)
            vals = applicable(c.calc.fn, ctx) ? c.calc.fn(ctx) : c.calc.fn()
            append!(F, (vals .- c.shift) ./ c.scale)
        end
    end
    return F
end

# ─────────────────────────────────────────────────────────────────────────────
# get_constraint_bounds
# ─────────────────────────────────────────────────────────────────────────────

function get_constraint_bounds(p::MGAnDSMsPhase)
    lb = (zeros(7) .- p.matchpoint_shift) ./ p.matchpoint_scale
    ub = copy(lb)
    # Alpha-sum equality — only when α is a registered decision variable.
    if p.alpha_var !== nothing
        push!(lb, 0.0)
        push!(ub, 0.0)
    end
    for c in p.constraints
        if c isa BoundaryConstraint
            append!(lb, (c.lower_bounds .- c.shift) ./ c.scale)
            append!(ub, (c.upper_bounds .- c.shift) ./ c.scale)
        end
    end
    return lb, ub
end

n_constraints(p::MGAnDSMsPhase) =
    7 + (p.alpha_var !== nothing ? 1 : 0) +
    sum(length(c.lower_bounds) for c in p.constraints
        if c isa BoundaryConstraint; init = 0)

# ─────────────────────────────────────────────────────────────────────────────
# set_initial_guess!
#
# Sets phase-cached decision variable values and syncs var.value fields in
# one call.  Any keyword omitted leaves the corresponding field unchanged.
#
#   set_initial_guess!(phase;
#       vinf_dep = [...],   # 3-vector (km/s or m/s)
#       vinf_arr = [...],   # 3-vector (km/s or m/s)
#       dv       = [...],   # 3×n_dsm matrix (km/s or m/s)
#       alpha    = [...],   # (n_dsm+1)-vector (arc fractions, sum = 1)
#       t0       = 0.0,     # scalar (s)
#       tf       = ...,     # scalar (s)
#       m0       = ...,     # scalar (kg)
#       mf       = ...,     # scalar (kg)
#   )
# ─────────────────────────────────────────────────────────────────────────────

function set_initial_guess!(p::MGAnDSMsPhase;
                             vinf_dep = nothing,
                             vinf_arr = nothing,
                             dv       = nothing,
                             alpha    = nothing,
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
    if dv !== nothing
        p._dv = Float64.(dv)
        p.dv_var !== nothing && (p.dv_var.value = vec(p._dv))
    end
    if alpha !== nothing
        p._alpha = Float64.(alpha)
        p.alpha_var !== nothing && (p.alpha_var.value = copy(p._alpha))
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
# ─────────────────────────────────────────────────────────────────────────────

# Derived from the evaluation context, the same rule the collocation phase uses.
# MGABoundaryContext carries x0, xf, vinf_dep, vinf_arr, t0, tf, m0, mf, dv and
# alpha, which is every variable `variable_list` returns, so a boundary function
# here can reach all of them.
#
# What this replaces enumerated those same eight variables with `===` and then
# OR'd in `has_jacobian`. It gave the same answer by a longer route, and the
# `has_jacobian` term could never change it: a variable reaches this loop only
# by being in `variable_list`, and every one of those was already named.
function sparsity_structure(p::MGAnDSMsPhase, pf::PhaseFunction)
    vlist = variable_list(p)
    if pf.source isa MGAMatchPointBlock || pf.source isa BoundaryConstraint
        return fill(true, length(vlist))
    elseif pf.source isa MGAAlphaSumBlock
        # sum(alpha) == 1 reads the alphas and nothing else. Exact, not a bound.
        return Bool[v === p.alpha_var for v in vlist]
    end
    return fill(false, length(vlist))
end

function sparsity_structure(p::MGAnDSMsPhase)
    flist = function_list(p)
    vlist = variable_list(p)
    Bool[sparsity_structure(p, f)[j] for f in flist, j in eachindex(vlist)]
end

# ─────────────────────────────────────────────────────────────────────────────
# jacobian_chunk
# ─────────────────────────────────────────────────────────────────────────────

function jacobian_chunk(p::MGAnDSMsPhase, pf::PhaseFunction,
                        var::DirectSolverVariable)
    var_sc = _var_nlp_scale(p, var)
    if pf.source isa MGAMatchPointBlock
        jac = matchpoint_jacobian(p)
        oid = objectid(var)
        chunk_phys =
            !isnothing(p.vinf_dep_var) && oid == objectid(p.vinf_dep_var) ? jac.dvinf_dep                   :
            !isnothing(p.vinf_arr_var) && oid == objectid(p.vinf_arr_var) ? jac.dvinf_arr                   :
            !isnothing(p.dv_var)       && oid == objectid(p.dv_var)       ? jac.ddv                         :
            !isnothing(p.alpha_var)    && oid == objectid(p.alpha_var)    ? jac.dalpha                      :
            !isnothing(p.t0_var)       && oid == objectid(p.t0_var)       ? reshape(jac.dt0, 7, 1)          :
            !isnothing(p.tf_var)       && oid == objectid(p.tf_var)       ? reshape(jac.dtf, 7, 1)          :
            !isnothing(p.m0_var)       && oid == objectid(p.m0_var)       ? reshape(jac.dm0, 7, 1)          :
            !isnothing(p.mf_var)       && oid == objectid(p.mf_var)       ? reshape(jac.dmf, 7, 1)          :
            throw(ArgumentError(
                "jacobian_chunk: $(repr(var.name)) must be a variable declared on phase " *
                ":$(p.name); declare it with Vary before asking for its Jacobian"))
        return chunk_phys .* transpose(var_sc) ./ p.matchpoint_scale

    elseif pf.source isa MGAAlphaSumBlock
        oid = objectid(var)
        if !isnothing(p.alpha_var) && oid == objectid(p.alpha_var)
            # ∂(Σα)/∂α = ones; no scaling on this constraint
            jac_al = alpha_sum_jacobian(p)
            return jac_al.dalpha .* transpose(var_sc)
        end
        return zeros(1, nlp_length(p, var))

    elseif pf.source isa BoundaryConstraint
        bc = pf.source
        chunk_phys = has_jacobian(bc.calc, var) ?
            get_jacobian(bc.calc, var) :
            _bc_jacobian_ad_chunk(p, bc.calc, var)
        return chunk_phys .* transpose(var_sc) ./ bc.scale
    end
    throw(ArgumentError("MGAnDSMsPhase.jacobian_chunk: unknown function source " *
        "$(typeof(pf.source)). Returning zeros here would be a wrong " *
        "derivative rather than a missing one."))
end

# ─────────────────────────────────────────────────────────────────────────────
# Boundary-constraint Jacobian by AD
#
# The counterpart of an analytic Jacobian registered with add_jacobian!. Every
# other phase type already falls back to AD when none is registered; this phase
# returned zeros instead, which is a wrong answer rather than a slow one, and it
# is what made uc7's gradient vanish.
#
# Differentiates with respect to the *physical* value, so the caller applies the
# same scaling it applies to a registered Jacobian.
# ─────────────────────────────────────────────────────────────────────────────

function _bc_jacobian_ad_chunk(p::MGAnDSMsPhase, bf::BoundaryFunction,
                               var::DirectSolverVariable)
    ctx_f64 = _mga_boundary_context(p)
    return ForwardDiff.jacobian(
        phys -> bf.fn(_mga_context_replace(ctx_f64, p, var, phys)),
        copy(var.value))
end

# ─────────────────────────────────────────────────────────────────────────────
# get_objective / objective_gradient_chunk
# ─────────────────────────────────────────────────────────────────────────────

function get_objective(p::MGAnDSMsPhase)
    isnothing(p.objective) && return 0.0
    obj = p.objective
    ctx = _mga_boundary_context(p)
    raw = applicable(obj.fn, ctx) ? obj.fn(ctx) : obj.fn()
    val = Float64(raw)
    return obj.sense === :Max ? -val : val
end

function objective_gradient_chunk(p::MGAnDSMsPhase, var::DirectSolverVariable)
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
    ctx_f64 = _mga_boundary_context(p)
    applicable(obj.fn, ctx_f64) || return zeros(nlp_length(p, var))
    g_phys = ForwardDiff.gradient(
        phys -> obj.fn(_mga_context_replace(ctx_f64, p, var, phys)), copy(var.value))
    g_scaled = g_phys .* _var_nlp_scale(p, var)
    return obj.sense === :Max ? .-g_scaled : g_scaled
end

# Quantities on an MGAnDSMs phase.
departure_vinf(p::MGAnDSMsPhase) = p._vinf_dep
arrival_vinf(p::MGAnDSMsPhase)   = p._vinf_arr
initial_time(p::MGAnDSMsPhase)   = p._t0
final_time(p::MGAnDSMsPhase)     = p._tf
initial_mass(p::MGAnDSMsPhase)   = p._m0
final_mass(p::MGAnDSMsPhase)     = p._mf
deep_space_dv(p::MGAnDSMsPhase)  = p._dv
arc_fractions(p::MGAnDSMsPhase)  = p._alpha
