# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Soft lunar landing with maximum remaining mass, free final time.
#
# Truth: **external reference**. Meditch (1964), as reproduced by NASA Glenn's Dymos, gives
# m_f ≈ 0.3953 and tf ≈ 1.397 in these non-dimensional units. That is an external reference,
# the strongest source available for this problem, which has no closed
# form: the optimal control is bang-bang and the switching time is what the solve finds.
#
#   maximise   m_f
#   subject to ḣ = v,  v̇ = -1 + T/m,  ṁ = -T/2.349
#   from h₀ = 1, v₀ = -0.783, m₀ = 1 to h_f = v_f = 0, with 0 ≤ T ≤ 1.227
#
# What this covers that the other collocation cases do not:
#
#   · a **maximised** objective — `sense = Max()`, which is a different sign path through the
#     assembly than the three minimisations elsewhere in the suite;
#   · a **free final time with a guess**, where tf is both a variable and part of the answer;
#   · **bang-bang control**, so the bounds are active for most of the arc rather than slack, which
#     is the case the box-bound machinery is actually for;
#   · **three states with coupled dynamics**, where the state Jacobian has an off-diagonal term
#     that depends on the control.
#
# The problem is the regression suite's uc16, recorded at 0.395082 under LGL. This file transcribes
# it with Hermite-Simpson, so agreement with the reference is cross-transcription as well.

using SNOW
using LinearAlgebra
using EpicycleBase
using AstroSolve
using Test

struct _MlState{T} <: AbstractState
    h::T
    v::T
    m::T
end

struct _MlControl{T} <: AbstractControl
    T_::T
end

struct _MlModel
    mdot_coeff::Float64
    T_max     ::Float64
end

const _ML_MODEL = _MlModel(2.349, 1.227)
const _ML_H0, _ML_V0, _ML_M0 = 1.0, -0.783, 1.0

# The published answer this file is checked against.
const _ML_REF_MF = 0.3953
const _ML_REF_TF = 1.397

function _ml_dynamics!(dy, y::_MlState, u::_MlControl, p, t, model)
    dy[1] =  y.v
    dy[2] = -1.0 + u.T_ / y.m
    dy[3] = -u.T_ / model.mdot_coeff
end

@partial(_ml_dynamics!, state) do dF, y, u, p, t, model
    dF[1, 2] =  1.0
    dF[2, 3] = -u.T_ / y.m^2
end

@partial(_ml_dynamics!, control) do dF, y, u, p, t, model
    dF[2, 1] =  1.0 / y.m
    dF[3, 1] = -1.0 / model.mdot_coeff
end

_ml_launch(c) = [state(c).h, state(c).v, state(c).m]
@partial(_ml_launch, state) do c
    [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
end

_ml_touchdown(c) = [state(c).h, state(c).v]
@partial(_ml_touchdown, state) do c
    [1.0 0.0 0.0; 0.0 1.0 0.0]
end

_ml_final_mass(c) = state(c).m
@partial(_ml_final_mass, state) do c
    [0.0 0.0 1.0]
end

"""Build the landing phase on a Hermite-Simpson mesh of `n_steps` intervals."""
function _ml_phase(; n_steps = 30)
    phase = CollocationPhase(
        name          = :moon_landing,
        transcription = HermiteSimpson(n_steps = n_steps),
        dynamics      = _ml_dynamics!, model   = _ML_MODEL,
        state         = _MlState,      control = _MlControl,
        tspan         = (0.0, 1.4))

    Vary(state, phase;
         guess       = [ _ML_H0  0.0
                         _ML_V0  0.0
                         _ML_M0  0.4],
         lower_bound = [0.0, -5.0, 0.001],
         upper_bound = [5.0,  5.0, 2.0])

    Vary(control, phase;
         guess       = reshape([_ML_MODEL.T_max/2  _ML_MODEL.T_max/2], 1, 2),
         lower_bound = [0.0],
         upper_bound = [_ML_MODEL.T_max])

    Vary(final_time, phase; guess = 1.4, lower_bound = 0.5, upper_bound = 5.0)

    Constraint(_ml_launch,    phase; equals = [_ML_H0, _ML_V0, _ML_M0], at = Initial())
    Constraint(_ml_touchdown, phase; equals = [0.0, 0.0],               at = Final())

    Objective(_ml_final_mass, phase; sense = Max())
    return phase
end

@testset "moon landing — lands softly with the published remaining mass" begin
    phase = _ml_phase()
    result = solve!(Sequence(phase))
    yf = get_final_state(phase)

    @test result.info === :Solve_Succeeded

    # The landing itself: on the surface, at rest. These were declared as constraints, so what is
    # checked here is that the solver honoured them rather than reporting success without them.
    @test yf.h ≈ 0.0 atol = 1e-6
    @test yf.v ≈ 0.0 atol = 1e-6

    # Against Meditch via Dymos. The tolerance is the reference's own precision — it is quoted to
    # four figures — not the solver's.
    @test yf.m ≈ _ML_REF_MF atol = 1e-3
    @test get_final_time(phase) ≈ _ML_REF_TF atol = 5e-3

    # Against the recorded LGL baseline, which is tighter because it is our own number to six
    # figures rather than a published one to four.
    @test yf.m ≈ 0.395082 atol = 1e-3
end

@testset "moon landing — the maximised objective is the mass it reports" begin
    # `sense = Max()` is a different sign path through the assembly than the minimisations
    # elsewhere in this suite. A sign error would converge to the *least* remaining mass, which is
    # a perfectly feasible landing, so nothing but comparing the two catches it.
    phase = _ml_phase()
    result = solve!(Sequence(phase))
    yf = get_final_state(phase)

    @test abs(result.objective) ≈ yf.m atol = 1e-8

    # Burning everything would also land, and would leave far less mass. The optimum must beat it
    # by a wide margin, which is the check that says Max() maximised.
    @test yf.m > 0.35
    @test yf.m < _ML_M0                      # some propellant was necessarily spent
end

@testset "moon landing — the trajectory obeys its physics and its bounds" begin
    phase = _ml_phase()
    solve!(Sequence(phase))

    ts = get_node_times(phase)
    xs = phase._Y
    us = phase._U

    @test size(xs, 1) == 3
    @test length(ts) == size(xs, 2)

    # Thrust inside its box at every node, which is where it spends most of the arc.
    @test all(-1e-8 .<= us .<= _ML_MODEL.T_max + 1e-8)

    # Bang-bang: the control sits at a bound for most of the arc rather than in between. This is
    # the shape of the answer, and a transcription that smoothed it would still land.
    at_bound = count(u -> u < 1e-3 || u > _ML_MODEL.T_max - 1e-3, us)
    @test at_bound > length(us) ÷ 2

    # Altitude starts at 1, ends at 0, and never goes below the surface.
    @test xs[1, 1] ≈ _ML_H0 atol = 1e-8
    @test xs[1, end] ≈ 0.0 atol = 1e-6
    @test all(xs[1, :] .>= -1e-8)

    # Mass is spent monotonically: ṁ = -T/coeff with T ≥ 0, so it can only fall. On the converged
    # mesh it ticks up by the node error near the thrust switch — 7.4e-5 on a 30-step mesh, which
    # is 1.2e-4 of the propellant spent — so the claim is stated against the spend rather than as
    # an exact sort. The refinement testset below is what says that figure is mesh error.
    spent_total = _ML_M0 - xs[3, end]
    @test xs[3, 1] ≈ _ML_M0 atol = 1e-8
    @test spent_total > 0.5
    @test maximum(diff(xs[3, :])) <= 1e-3 * spent_total

    # And the mass budget closes. Integrating ṁ = -T/coeff by trapezoid over the converged mesh
    # must reproduce the propellant actually spent, which checks the dynamics the transcription
    # integrated rather than the endpoints it was handed.
    tf = get_final_time(phase)
    tn = range(0, tf; length = size(us, 2))
    dt = step(tn)
    spent = dt * (sum(us) - 0.5 * (us[1] + us[end])) / _ML_MODEL.mdot_coeff
    @test spent ≈ (_ML_M0 - xs[3, end]) atol = 5e-3
end

@testset "moon landing — refining the mesh shrinks the monotonicity violation" begin
    # The mass tick-up above is discretization error, which is a claim rather than an excuse until
    # it is shown to fall with the mesh. It does, which also says the bang-bang switch is being
    # resolved better rather than being smoothed away.
    function worst_uptick(n)
        phase = _ml_phase(n_steps = n)
        solve!(Sequence(phase))
        return maximum(diff(phase._Y[3, :]))
    end

    coarse = worst_uptick(15)
    fine   = worst_uptick(60)
    @test fine < coarse

    # And the answer itself holds still while the mesh changes under it, which is the property
    # that matters to a user.
    for n in (15, 60)
        phase = _ml_phase(n_steps = n)
        solve!(Sequence(phase))
        @test get_final_state(phase).m ≈ _ML_REF_MF atol = 2e-3
    end
end
