# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# A Sims-Flanagan phase handing off to a collocation phase, written in the spec vocabulary.
#
# `Constraint(continuity, Link(a, b))` joined two collocation phases only. A shooting phase joined
# a collocation phase through `add_continuity!`, which is the development interface, so the
# examples that fly one transcription into another could not be written from the exports. Two
# things had to change: the link had to route a mixed pair to the unified manager, and the pending
# links had to be registered before the sequence chose how to initialise itself. Registered after,
# a sequence holding a Sims-Flanagan phase was initialised as a shooting sequence, its collocation
# phase never was, and IPOPT stopped at iteration 0 on an invalid number.
#
# Truth: **two implementations of one condition**. The link written in the vocabulary must give the
# rows `add_continuity!` gives, at the same point, and the sequence must assemble with no NaN.
# The Sims-Flanagan setup that used to need field assignments (`tspan`, `matchpoint_scale`, a
# `Vary` guess, an objective) is checked against the values it was given.

using LinearAlgebra
using EpicycleBase
using AstroModels
using AstroCallbacks
using AstroSolve
using AstroSolve: add_continuity!, get_functions, get_jacobian, initialize_sequence!,
    objective_gradient_chunk, sequence_oc_links
using Test

const _ML_MU, _ML_M0, _ML_TMAX, _ML_MDOT = 1.0, 1.0, 0.1405, 0.0749
const _ML_TMID, _ML_TF = 1.66, 3.32

struct _MlState{T}   <: AbstractState;   x::T; y::T; z::T; vx::T; vy::T; vz::T; m::T end
struct _MlControl{T} <: AbstractControl; ux::T; uy::T; uz::T end
struct _MlModel; mu::Float64; thrust::Float64; mdot::Float64 end
struct _MlSubjState{T} <: AbstractState; x::T; y::T; z::T; vx::T; vy::T; vz::T; m::T end

function _ml_dynamics!(dy, y, u, p, t, m)
    r = sqrt(y.x^2 + y.y^2 + y.z^2)
    a = m.thrust / y.m
    dy[1] = y.vx;  dy[2] = y.vy;  dy[3] = y.vz
    dy[4] = -m.mu * y.x / r^3 + a * u.ux
    dy[5] = -m.mu * y.y / r^3 + a * u.uy
    dy[6] = -m.mu * y.z / r^3 + a * u.uz
    dy[7] = -m.mdot * sqrt(u.ux^2 + u.uy^2 + u.uz^2 + 1e-12)
end

# The guess is an arc flown from a circular orbit at r = 1 with thrust along the velocity, so the
# two legs start as one trajectory and every defect starts small.
function _ml_arc(t_end; steps = 1000, y0 = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, _ML_M0])
    f(y) = (r = sqrt(y[1]^2 + y[2]^2 + y[3]^2); v = sqrt(y[4]^2 + y[5]^2 + y[6]^2);
            vcat(y[4:6], -_ML_MU .* y[1:3] ./ r^3 .+ (_ML_TMAX / y[7]) .* y[4:6] ./ v, -_ML_MDOT))
    y, h = copy(y0), t_end / steps
    for _ in 1:steps
        k1 = f(y); k2 = f(y .+ h/2 .* k1); k3 = f(y .+ h/2 .* k2); k4 = f(y .+ h .* k3)
        y = y .+ (h/6) .* (k1 .+ 2k2 .+ 2k3 .+ k4)
    end
    return y
end

const _ML_Y0   = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
const _ML_YMID = _ml_arc(_ML_TMID)
const _ML_YEND = _ml_arc(_ML_TF)

"""The two legs, each fully declared. `coast` asks for the most final mass instead of radius."""
function _ml_phases(; coast = false)
    sf = SimsFlanaganPhase(name = :first_leg, transcription = SimsFlanagan(n_segments = 8),
                           model = PropulsionModel(mu = _ML_MU, Isp = _ML_TMAX / _ML_MDOT,
                                                   Tmax = _ML_TMAX, g0 = 1.0),
                           tspan = (0.0, _ML_TMID))
    Vary(state,       sf; guess = _ML_Y0, lower_bound = _ML_Y0, upper_bound = _ML_Y0)
    Vary(final_state, sf; guess = _ML_YMID[1:6],
                          lower_bound = fill(-5.0, 6), upper_bound = fill(5.0, 6))
    Vary(forward_control,  sf; guess = [0.0, 0.5, 0.0],
                               lower_bound = fill(-1.0, 3), upper_bound = fill(1.0, 3))
    Vary(backward_control, sf; lower_bound = fill(-1.0, 3), upper_bound = fill(1.0, 3))
    Vary(initial_mass, sf; lower_bound = _ML_M0, upper_bound = _ML_M0)
    Vary(final_mass,   sf; guess = _ML_YMID[7], lower_bound = 0.5, upper_bound = _ML_M0)

    hs = CollocationPhase(name = :second_leg, transcription = HermiteSimpson(n_steps = 6),
                          dynamics = _ml_dynamics!, model = _MlModel(_ML_MU, _ML_TMAX, _ML_MDOT),
                          state = _MlState, control = _MlControl, tspan = (_ML_TMID, _ML_TF))
    Vary(state, hs; guess = hcat(_ML_YMID, _ML_YEND),
                    lower_bound = vcat(fill(-5.0, 6), 0.5), upper_bound = vcat(fill(5.0, 6), 1.0))
    Vary(control, hs; guess = repeat([0.0, 1.0, 0.0], 1, 2),
                      lower_bound = fill(-1.0, 3), upper_bound = fill(1.0, 3))
    final_radius(c) = (y = state(c); sqrt(y.x^2 + y.y^2 + y.z^2))
    final_mass_c(c) = state(c).m
    Objective(coast ? final_mass_c : final_radius, hs; sense = Max())
    return sf, hs
end

function _ml_sequence(join; coast = false)
    sf, hs = _ml_phases(; coast = coast)
    seq = Sequence()
    add_sequence!(seq, sf)
    add_sequence!(seq, hs)
    join(seq, sf, hs)
    initialize_sequence!(seq)
    return seq, sf, hs
end

@testset "Mixed links — Constraint(continuity, Link) joins a shooting and a collocation phase" begin
    seq, sf, hs = _ml_sequence((seq, a, b) -> Constraint(continuity, Link(a, b)))

    # Registered before the sequence chose how to initialise, so it went to the unified manager.
    links = sequence_oc_links(seq)
    @test length(links) == 1
    @test length(links[1].lower_bounds) == 8        # six states, the mass, and the time

    om = AstroSolve._get_or_build_oc(seq)
    F = get_functions(om)
    J = get_jacobian(om)
    @test !any(isnan, F)                            # the collocation leg was initialised
    @test !any(isnan, J)

    # The rows are what add_continuity! writes for the same pair, at the same point.
    ref, _, _ = _ml_sequence((seq, a, b) -> add_continuity!(seq, a, b))
    omr = AstroSolve._get_or_build_oc(ref)
    @test F ≈ get_functions(omr) atol = 1e-12
    @test J ≈ get_jacobian(omr) atol = 1e-8

    # And the link residual is the handoff: the Sims-Flanagan end, mass last, minus the start of
    # the collocation leg, then the time.
    c1 = AstroSolve._unified_boundary_context(sf)
    c2 = AstroSolve._unified_boundary_context(hs)
    @test F[end-7:end] ≈ vcat(c1.yf .- c2.y0, c1.tf - c2.t0) atol = 1e-12
end

@testset "Mixed links — a function of its own on a Link reads both ends by name" begin
    handoff(c1, c2) = initial_state(c2) .- final_state(c1)
    seq, sf, hs = _ml_sequence((seq, a, b) ->
        Constraint(handoff, Link(a, b); equals = zeros(7), name = :handoff))
    F = get_functions(AstroSolve._get_or_build_oc(seq))
    @test !any(isnan, F)
    c1 = AstroSolve._unified_boundary_context(sf)
    c2 = AstroSolve._unified_boundary_context(hs)
    @test F[end-6:end] ≈ c2.y0 .- c1.yf atol = 1e-12
end

@testset "Mixed links — a solve joins the two legs" begin
    # The most final mass is a coast: no thrust on either leg, so the mass stays at its initial
    # value and the orbit stays circular at r = 1. That answer is known, and it is reached only if
    # the link carries the state and the mass across the handoff.
    #
    # IPOPT's exit status is not the claim. Mass flow goes as |u|, which has a kink at zero thrust,
    # so the coast optimum is degenerate and IPOPT reaches it without meeting its termination test.
    # The answer it reaches is checked instead, to 1e-5 in mass and 1e-3 in radius.
    seq, sf, hs = _ml_sequence((seq, a, b) -> Constraint(continuity, Link(a, b)); coast = true)
    solve!(seq; method = Optimize(max_iter = 300, tol = 1e-6, print_level = 0))
    Y = state(hs)
    @test maximum(abs.(Y[:, 1] .- vcat(final_state(sf), final_mass(sf)))) < 1e-8
    @test get_initial_time(hs) ≈ final_time(sf)
    @test Y[7, end] ≈ _ML_M0 atol = 1e-5
    @test norm(Y[1:3, end]) ≈ 1.0 atol = 1e-3
end

@testset "Sims-Flanagan — setup the vocabulary now carries" begin
    sf, _ = _ml_phases()

    # tspan sets the times a phase keeps when they are not varied.
    @test initial_time(sf) == 0.0
    @test final_time(sf) == _ML_TMID

    # A guess reaches the phase, one direction repeated over every segment; equal bounds pin.
    @test forward_control(sf) == repeat([0.0, 0.5, 0.0], 1, 4)
    @test backward_control(sf) == zeros(3, 4)       # no guess: left as it was
    @test initial_mass(sf) == _ML_M0
    @test final_mass(sf) == _ML_YMID[7]

    # The match-point scaling is a constructor keyword, and a wrong length or a backwards span is
    # refused by name rather than failing inside the NLP.
    kw = (name = :s, transcription = SimsFlanagan(n_segments = 4),
          model = PropulsionModel(mu = 1.0, Isp = 1.0, Tmax = 0.1, g0 = 1.0))
    scaled = SimsFlanaganPhase(; kw..., matchpoint_scale = fill(2.0, 7))
    @test scaled.matchpoint_scale == fill(2.0, 7)
    @test_throws ArgumentError SimsFlanaganPhase(; kw..., matchpoint_scale = ones(6))
    @test_throws ArgumentError SimsFlanaganPhase(; kw..., tspan = (1.0, 0.0))
    @test_throws ArgumentError Vary(forward_control, sf; guess = ones(5),
                                    lower_bound = fill(-1.0, 3), upper_bound = fill(1.0, 3))
end

@testset "Sims-Flanagan — an Objective reads the boundary context, declared or differentiated" begin
    # The regularised final mass of the Mars example. The declared partials and the framework's
    # differentiation of the same function must agree block by block.
    effort(c) = sum(abs2, forward_control(c)) + sum(abs2, backward_control(c))
    declared(c) = final_mass(c) - 1e-3 * effort(c)
    @partial(declared, final_mass)       do c; [1.0]; end
    @partial(declared, forward_control)  do c; -2e-3 .* vec(forward_control(c));  end
    @partial(declared, backward_control) do c; -2e-3 .* vec(backward_control(c)); end
    bare(c) = final_mass(c) - 1e-3 * effort(c)

    a, _ = _ml_phases(); Objective(declared, a; sense = Max())
    b, _ = _ml_phases(); Objective(bare,     b; sense = Max())
    for (va, vb) in ((a.mf_var, b.mf_var), (a.u_fwd_var, b.u_fwd_var), (a.u_bwd_var, b.u_bwd_var))
        @test objective_gradient_chunk(a, va) ≈ objective_gradient_chunk(b, vb) atol = 1e-12
    end
    @test_throws ArgumentError Objective(declared, SimsFlanaganPhase(
        name = :bare, transcription = SimsFlanagan(n_segments = 4),
        model = PropulsionModel(mu = 1.0, Isp = 1.0, Tmax = 0.1, g0 = 1.0)); sense = Max())
end

@testset "subject_at — a solved phase's end as a subject quantities read" begin
    # A phase flown by a force model keeps the spacecraft as its model; its time is in seconds.
    sc = Spacecraft()
    p = CollocationPhase(name = :subject, transcription = HermiteSimpson(n_steps = 2),
                         model = sc, state = _MlSubjState, tspan = (0.0, 600.0))
    p._Y  = hcat([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0, 1000.0],
                 [0.0, 7100.0, 0.0, -7.4, 0.0, 0.0, 990.0])
    p._tf = 600.0

    fin = subject_at(p, Final())
    @test position_magnitude(fin) ≈ 7100.0
    @test fin.coord_sys === sc.coord_sys
    @test (fin.time - sc.time) * 86400 ≈ 600.0 atol = 1e-4     # a Julian date's resolution
    @test position_magnitude(subject_at(p, Initial())) ≈ 7000.0
end
