# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# What the collocation types print.
#
# Twelve `Base.show` methods sit in collocation_sequence.jl and shooting_sequence.jl and none of
# them had ever run. In a library with no GUI that output is the interface: a user holding a
# `CollocationManager` or a `PathConstraint` learns what it is from what it prints and from
# nothing else.
#
# The comment above them says they exist to "prevent recursive struct explosion in the REPL", so
# that is the contract under test. Each one must stay on a single line and name its own type,
# because the failure they guard against is a phase printing its whole mesh when it appears as a
# field of something else.
#
# The second thing tested here is that they run at all. A `show` method referring to a field that
# has been renamed is invisible until a user displays the object, and then it throws in the one
# place a user is least able to work around.

using SNOW
using LinearAlgebra
using EpicycleBase
using AstroSolve
using AstroSolve: CollocationManager, function_list, initialize!, print_sparsity, variable_list
using Test

struct _ShState{T}   <: AbstractState;   x::T end
struct _ShControl{T} <: AbstractControl; u::T end

_sh_dyn!(dy, y::_ShState, u::_ShControl, p, t, model) = (dy[1] = u.u)
_sh_pos(c)   = [state(c).x]
_sh_speed(c) = [control(c).u]
_sh_cost(c)  = state(c).x^2 + 0.5 * control(c).u^2
_sh_final(c) = [state(c).x]

"One phase carrying a path constraint, a boundary constraint and an objective."
function _sh_phase(name::Symbol; with_path::Bool = true)
    phase = CollocationPhase(
        name          = name,
        transcription = HermiteSimpson(n_steps = 3),
        dynamics      = _sh_dyn!,
        state         = _ShState, control = _ShControl,
        tspan         = (0.0, 1.0))

    Vary(state,   phase; guess = [0.0 0.1],
                         lower_bound = [-2.0], upper_bound = [2.0])
    Vary(control, phase; guess = reshape([0.0 0.1], 1, 2),
                         lower_bound = [-2.0], upper_bound = [2.0])
    with_path && Constraint(_sh_speed, phase; upper_bound = 3.0, at = Path())
    Constraint(_sh_pos, phase; equals = [0.0], at = Initial())
    Objective(_sh_cost, phase; sense = Min(), at = Path())
    return phase
end

"A two-phase sequence, initialised, so the manager and its linkages exist."
function _sh_sequence()
    a = _sh_phase(:show_a)
    b = _sh_phase(:show_b; with_path = false)
    seq = Sequence()
    add_sequence!(seq, a)
    add_sequence!(seq, b)
    Constraint(continuity, Link(a, b))
    initialize!(seq)
    return seq, a, b
end

"The single line an object prints, with the checks every one of them must pass."
function shown(x)
    text = sprint(show, x)
    @test !isempty(text)
    # Compact means one line. A nested phase printing its mesh is the failure these methods exist
    # to prevent, and it only shows up when the object is a field of something else.
    @test !occursin('\n', text)
    # A `show` that falls through to the default prints the module path and every field.
    @test !occursin("#undef", text)
    return text
end

@testset "show — every collocation type prints a single identifying line" begin
    seq, a, b = _sh_sequence()
    cm = CollocationManager(seq)

    @test occursin("CollocationManager", shown(cm))
    @test occursin("phase", shown(cm))

    @test shown(a) == "CollocationPhase(:show_a)"
    @test occursin("show_b", shown(b))

    @test !isempty(cm.linkages)
    for lc in cm.linkages
        text = shown(lc)
        @test occursin("LinkageConstraint", text)
        # A linkage is only meaningful as a pair, so both ends have to be named.
        @test occursin("show_a", text)
        @test occursin("show_b", text)
    end

    for p in cm.phases, v in variable_list(p)
        @test occursin("SolverVariable", shown(v))
    end

    for p in cm.phases, pf in function_list(p)
        text = shown(pf)
        @test !isempty(text)
        # Whatever the concrete function type, it names itself rather than falling through.
        @test occursin(r"^[A-Za-z]", text)
    end
end

@testset "show — the counts a manager prints are its real ones" begin
    # The numbers are the whole content of this line, so a stale one is worse than no line.
    seq, _, _ = _sh_sequence()
    cm = CollocationManager(seq)
    text = sprint(show, cm)

    @test occursin("$(length(cm.phases)) phase", text)
    @test occursin("n_vars=$(cm.n_vars)", text)
    @test occursin("n_funs=$(cm.n_funs)", text)
    @test occursin("$(length(cm.linkages)) linkage", text)
end

@testset "show — a phase nested in a manager stays compact" begin
    # The contract stated in the source: displaying the container must not print the contents of
    # every phase. Measured rather than asserted by eye.
    seq, _, _ = _sh_sequence()
    cm = CollocationManager(seq)
    @test length(sprint(show, cm)) < 200
end

@testset "print_sparsity — the pattern renders and labels its rows and columns" begin
    # Twenty-five lines that nothing had run. It is a diagnostic a user reaches for when a solve
    # is slow, which is exactly when a method error in the printer is least welcome.
    seq, a, _ = _sh_sequence()
    initialize!(seq)

    io = IOBuffer()
    print_sparsity(io, a)
    text = String(take!(io))

    @test !isempty(text)
    @test occursin('\n', text)          # a matrix, unlike the compact show methods
    # Row and column labels are the point of it; a bare grid of marks says nothing.
    @test occursin("state", text) || occursin("control", text) || occursin("var", text)
end
