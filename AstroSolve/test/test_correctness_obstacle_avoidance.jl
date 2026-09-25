# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Steering between two points while staying clear of two obstacles.
#
# Truth: **analytic**, for the objective. The vehicle cruises at a fixed speed, so the integrand
# `(V cos θ)² + (V sin θ)²` is identically `V²` whatever the heading does, and the cost is `V² tf`
# exactly. That makes the objective a closed-form number that no amount of steering can change —
# which is precisely what makes it a clean check on the quadrature, since the path constraint is
# active and the trajectory is not analytic at all.
#
#   minimise   ∫₀¹ (ẋ² + ẏ²) dt
#   subject to ẋ = V cos θ,  ẏ = V sin θ,  V = 2.138
#   from (0, 0) to (1.2, 1.6), outside two circles of radius π/10
#
# What this covers that the other collocation cases do not:
#
#   · a **vector path constraint** — two components, one per obstacle, evaluated at every node,
#     where the other cases constrain a single scalar along the path;
#   · a **lower** bound on a path constraint rather than an upper one;
#   · an **active** path constraint, where the optimum is pressed against the boundary instead of
#     sitting clear of it as in the tracking and brachistochrone cases;
#   · `scale =` on `Vary`, the variable scaling nothing else in the open suite reaches.
#
# The problem is the regression suite's uc19, recorded at 4.571044 under LGL. This file transcribes
# it with Hermite-Simpson.

using SNOW
using LinearAlgebra
using EpicycleBase
using AstroSolve
using Test

const _OA_V  = 2.138
const _OA_R  = π / 10
const _OA_R2 = _OA_R^2
const _OA_C1 = (0.4, 0.5)
const _OA_C2 = (0.8, 1.5)
const _OA_P0 = (0.0, 0.0)
const _OA_PF = (1.2, 1.6)
const _OA_TF = 1.0

struct _OaState{T}   <: AbstractState;   x::T; y::T end
struct _OaControl{T} <: AbstractControl; θ::T       end

function _oa_dynamics!(dy, y::_OaState, u::_OaControl, p, t, model)
    dy[1] = _OA_V * cos(u.θ)
    dy[2] = _OA_V * sin(u.θ)
end

# ∂f/∂y is structurally zero and left undeclared, so the state Jacobian falls back to AD.
@partial(_oa_dynamics!, control) do dF, y, u, p, t, model
    dF[1, 1] = -_OA_V * sin(u.θ)
    dF[2, 1] =  _OA_V * cos(u.θ)
end

_oa_position(c) = [state(c).x, state(c).y]
@partial(_oa_position, state) do c
    Matrix{Float64}(I, 2, 2)
end

# Two components, one per obstacle, bounded below: the squared distance to each centre must stay
# at or above the squared exclusion radius.
_oa_clearance(c) = [(state(c).x - _OA_C1[1])^2 + (state(c).y - _OA_C1[2])^2,
                    (state(c).x - _OA_C2[1])^2 + (state(c).y - _OA_C2[2])^2]
@partial(_oa_clearance, state) do c
    [2*(state(c).x - _OA_C1[1])  2*(state(c).y - _OA_C1[2])
     2*(state(c).x - _OA_C2[1])  2*(state(c).y - _OA_C2[2])]
end

_oa_speed_squared(c) = (_OA_V * cos(control(c).θ))^2 + (_OA_V * sin(control(c).θ))^2
@partial(_oa_speed_squared, control) do c
    θ = control(c).θ
    [-2*_OA_V^2 * cos(θ) * sin(θ) + 2*_OA_V^2 * sin(θ) * cos(θ)]
end

const _OA_WX = [_OA_P0[1], 0.3, 0.6, 0.9, _OA_PF[1]]
const _OA_WY = [_OA_P0[2], 0.4, 0.8, 1.2, _OA_PF[2]]

# Heading taken from consecutive waypoints, so the solver starts pointing the way it has to go.
const _OA_THETA = let θ = [atan(_OA_WY[i+1] - _OA_WY[i], _OA_WX[i+1] - _OA_WX[i]) for i in 1:4]
    reshape(push!(θ, θ[end]), 1, 5)
end

"""Build the avoidance phase on a Hermite-Simpson mesh of `n_steps` intervals."""
function _oa_phase(; n_steps = 35)
    phase = CollocationPhase(
        name          = :obstacle_avoidance,
        transcription = HermiteSimpson(n_steps = n_steps),
        dynamics      = _oa_dynamics!,
        state         = _OaState,      control = _OaControl,
        tspan         = (0.0, _OA_TF))

    Vary(state, phase;
         guess       = [_OA_WX'; _OA_WY'],
         lower_bound = [0.0, 0.0],
         upper_bound = [_OA_PF[1], _OA_PF[2]],
         scale       = [_OA_PF[1], _OA_PF[2]])

    Vary(control, phase; guess = _OA_THETA,
                         lower_bound = [-10.0], upper_bound = [10.0])

    Constraint(_oa_position,  phase; equals = [_OA_P0...], at = Initial())
    Constraint(_oa_position,  phase; equals = [_OA_PF...], at = Final())
    Constraint(_oa_clearance, phase; lower_bound = [_OA_R2, _OA_R2], at = Path())

    Objective(_oa_speed_squared, phase; sense = Min(), at = Path())
    return phase
end

@testset "obstacle avoidance — the cost is V² tf whatever the path does" begin
    # The integrand is V²(cos²θ + sin²θ) = V² for every heading, so the cost cannot depend on the
    # trajectory. That makes it a check on the quadrature alone: the path is being bent around two
    # obstacles at the same time, and the number must not move.
    phase = _oa_phase()
    result = solve!(Sequence(phase))

    @test result.info === :Solve_Succeeded
    @test result.objective ≈ _OA_V^2 * _OA_TF atol = 1e-6
    @test result.objective ≈ 4.571044 atol = 1e-5      # the recorded LGL baseline

    # The endpoints the problem is posed between.
    yf = get_final_state(phase)
    @test yf.x ≈ _OA_PF[1] atol = 1e-6
    @test yf.y ≈ _OA_PF[2] atol = 1e-6
    y0 = get_initial_state(phase)
    @test y0.x ≈ _OA_P0[1] atol = 1e-6
    @test y0.y ≈ _OA_P0[2] atol = 1e-6
end

@testset "obstacle avoidance — the path clears both obstacles at every node" begin
    # A vector path constraint with two components. Both must hold everywhere, and a transcription
    # that applied only the first would still reach the target and still report the same cost,
    # because the cost cannot see the path.
    phase = _oa_phase()
    solve!(Sequence(phase))

    xs = phase._Y
    @test size(xs, 1) == 2

    d1 = [(xs[1, k] - _OA_C1[1])^2 + (xs[2, k] - _OA_C1[2])^2 for k in axes(xs, 2)]
    d2 = [(xs[1, k] - _OA_C2[1])^2 + (xs[2, k] - _OA_C2[2])^2 for k in axes(xs, 2)]

    @test minimum(d1) >= _OA_R2 - 1e-6
    @test minimum(d2) >= _OA_R2 - 1e-6

    # And the straight line between the endpoints does not clear them, so the constraint is doing
    # something rather than being satisfied by the initial guess.
    straight = [( _OA_P0[1] + s*(_OA_PF[1]-_OA_P0[1]),
                  _OA_P0[2] + s*(_OA_PF[2]-_OA_P0[2]) ) for s in range(0, 1; length = 201)]
    worst1 = minimum((p[1]-_OA_C1[1])^2 + (p[2]-_OA_C1[2])^2 for p in straight)
    @test worst1 < _OA_R2
end

@testset "obstacle avoidance — the constraint is active, not slack" begin
    # This is the case the path machinery exists for, and the one the other collocation files do
    # not have: the optimum is pressed against at least one boundary rather than sitting clear of
    # it. If both obstacles were slack the test above would pass on a straight line.
    phase = _oa_phase()
    solve!(Sequence(phase))
    xs = phase._Y

    d1 = minimum((xs[1, k] - _OA_C1[1])^2 + (xs[2, k] - _OA_C1[2])^2 for k in axes(xs, 2))
    d2 = minimum((xs[1, k] - _OA_C2[1])^2 + (xs[2, k] - _OA_C2[2])^2 for k in axes(xs, 2))

    # At least one obstacle is grazed, to within the mesh's ability to resolve the tangency.
    @test min(d1 - _OA_R2, d2 - _OA_R2) < 5e-3
end

@testset "obstacle avoidance — the trajectory moves at the cruise speed throughout" begin
    # The dynamics fix the speed, so the distance between consecutive nodes must be V times the
    # time between them. This checks what the transcription integrated, which the cost cannot,
    # since the cost is constant by construction.
    phase = _oa_phase()
    solve!(Sequence(phase))

    ts = get_node_times(phase)
    xs = phase._Y

    # Chord versus arc. On the true trajectory the chord can only be shorter than V dt, since the
    # vehicle travels an arc of that length between the two nodes. On the converged mesh the node
    # positions carry Hermite-Simpson's discretization error, so the chord can exceed V dt by that
    # much — measured at 2e-4 relative here — and an exact one-sided bound is the wrong claim.
    # The tolerance is the mesh error; the interesting half is the lower bound, which says the
    # vehicle really is moving at the cruise speed rather than loitering.
    for k in 1:(length(ts) - 1)
        dt = ts[k+1] - ts[k]
        ds = hypot(xs[1, k+1] - xs[1, k], xs[2, k+1] - xs[2, k])
        @test ds <= _OA_V * dt * (1 + 1e-3)
        @test ds >= 0.95 * _OA_V * dt
    end

    # The whole arc, where the per-step errors cannot hide: total path length is V tf to within
    # the accumulated chord-versus-arc shortfall.
    total = sum(hypot(xs[1, k+1] - xs[1, k], xs[2, k+1] - xs[2, k])
                for k in 1:(size(xs, 2) - 1))
    @test total ≈ _OA_V * _OA_TF rtol = 1e-3
end
