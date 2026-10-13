# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# The live solver status box.
#
# While a solve runs, a box in the terminal shows each variable and goal by name, with its unit,
# value and error, and a ✓ or ✗ per goal, and is redrawn in place every iteration. The solver
# backend (ipopt_backend.jl) calls `show_iteration!` from IPOPT's per-iteration callback; nothing
# here knows about IPOPT. A solver hands over one solver-neutral `IterationRecord` per iteration,
# built by a `MonitorRows` adapter that knows the problem's variables and goals.
#
# Rules the display keeps:
#   - It never erases output it did not draw. A redraw moves up over exactly the lines the last
#     frame drew.
#   - It is at most 80 columns, narrows to fit a smaller terminal, and becomes one status line
#     below 64 columns.
#   - The width is checked before every frame. Wider: redraw in place. Narrower: the terminal has
#     re-wrapped the old frame, so its line count can't be trusted; leave it and draw a fresh frame
#     below.
#   - Messages logged during a solve go above the box (`_MonitorLogger`).
#   - When stdout is not an interactive terminal, only the final frame is printed.

import Base.CoreLogging: AbstractLogger, handle_message, shouldlog, min_enabled_level,
                         catch_exceptions, with_logger, current_logger

# ── On and off ────────────────────────────────────────────────────────────────────────────────────

const _MONITOR_DEFAULT = Ref(true)

"""
    solver_monitor!(on::Bool) -> Bool

Turn the live solver status box on or off for every solve that does not set `monitor` itself.
On by default.

# Arguments
- `on::Bool`: `true` to show the box, `false` for IPOPT's own log as before.

# Notes
A solve's own `Optimize(monitor = true | false)` takes precedence. With the box on, it is redrawn
in place when output goes to an interactive terminal; anywhere else (a file, CI, a docs build)
only its final frame is printed.

# Returns
`on`.

# Example
```julia
solver_monitor!(false)   # IPOPT's log for every solve from here on
solver_monitor!(true)    # back to the status box
```
"""
solver_monitor!(on::Bool) = (_MONITOR_DEFAULT[] = on)

# Test hooks: a fixed terminal width, and whether output counts as an interactive terminal.
const _MONITOR_COLS = Ref{Union{Nothing, Int}}(nothing)
const _MONITOR_LIVE = Ref{Union{Nothing, Bool}}(nothing)

_term_cols() = something(_MONITOR_COLS[], displaysize(stdout)[2])
_term_live() = something(_MONITOR_LIVE[], stdout isa Base.TTY)

# ── The solver-neutral iteration record ──────────────────────────────────────────────────────────

"""One variable row: a name, the unit shown and the value in that unit."""
struct MonitorVar
    name::String
    unit::String
    value::Float64
end

"""
One goal row. `desired` and `achieved` are display text, so a bound or a summary fits the same
columns as an equality. `error` is the violation in the displayed unit, 0 within bounds; `met`
is the solver's own verdict.
"""
struct MonitorGoal
    name::String
    desired::String
    achieved::String
    unit::String
    error::Float64
    met::Bool
end

"""
    IterationRecord

What one solver iteration reports to the status box. Every solver adapter produces this and the
display reads nothing else, so a new solver needs an adapter and no change to the display.
"""
struct IterationRecord
    solver::String
    problem::String
    iter::Int
    vars::Vector{MonitorVar}
    notes::Vector{String}          # summary lines under Variables, e.g. a decision-variable count
    goals::Vector{MonitorGoal}
    objective::Float64             # NaN for a feasibility problem
    inf_pr::Float64
    inf_du::Float64
    mu::Float64
    worst::Float64                 # the progress bar's measure: worst constraint violation
    tol::Float64                   # where the bar ends
    elapsed::Float64
    status::Symbol                 # :running, or the solver's exit status
end

_converged(s::Symbol) = s === :Solve_Succeeded || s === :Solved_To_Acceptable_Level ||
                        s === :Feasible_Point_Found

# ── Rows: what a problem's variables and goals are called ────────────────────────────────────────

"""
    MonitorRows

Turns a solver's current iterate into display rows. `monitor_rows(rows, x, g)` returns
`(vars, notes, goals, worst)`; `monitor_title(rows)` names the problem in the box's header.
"""
abstract type MonitorRows end

monitor_title(::MonitorRows) = "Optimizer"

# Radians read better as degrees in a status display; everything else shows in its own unit.
_shown(v::Real, u::AbstractString) = u == "rad" ? (rad2deg(v), "deg") : (Float64(v), String(u))

function _violation(g, lo, hi)
    g < lo && return lo - g
    g > hi && return g - hi
    return 0.0
end

# The desired column: the value for an equality, the bound for one side, both for a range.
function _desired_text(lo, hi, u)
    lo_s, _ = _shown(lo, u)
    hi_s, _ = _shown(hi, u)
    lo == hi                    && return _fmt(lo_s)
    isfinite(lo) && !isfinite(hi) && return "≥ " * strip(_fmt(lo_s))
    !isfinite(lo) && isfinite(hi) && return "≤ " * strip(_fmt(hi_s))
    isfinite(lo) && isfinite(hi)  && return @sprintf("%.5g…%.5g", lo_s, hi_s)
    return "free"
end

# A problem with no named structure: counts and worst violations only. Used by every solve path
# that has not supplied rows of its own, so every solver gets a box.
struct GenericRows <: MonitorRows
    title::String
    nx::Int
    lg::Vector{Float64}
    ug::Vector{Float64}
    ctol::Float64
end

monitor_title(r::GenericRows) = r.title

function monitor_rows(r::GenericRows, x, g)
    notes = ["$(r.nx) decision variables"]
    goals = MonitorGoal[]
    worst = 0.0
    for (label, pick) in (("Equalities", (lo, hi) -> lo == hi),
                          ("Inequalities", (lo, hi) -> lo != hi))
        idx = [i for i in eachindex(r.lg) if pick(r.lg[i], r.ug[i])]
        isempty(idx) && continue
        v   = [_violation(g[i], r.lg[i], r.ug[i]) for i in idx]
        met = count(<=(r.ctol), v)
        w   = maximum(v)
        worst = max(worst, w)
        push!(goals, MonitorGoal("$label ($(length(idx)))", "", "$met met", "", w,
                                 met == length(idx)))
    end
    return MonitorVar[], notes, goals, worst
end

# ── Rows for a targeting sequence (the event graph) ──────────────────────────────────────────────

_quantity(calc) = calc isa Calc ? calc.f : nothing
_quantity_label(calc) = (q = _quantity(calc); q === nothing ? string(nameof(typeof(calc))) :
                                                               EpicycleBase.label(q))
_quantity_unit(calc)  = (q = _quantity(calc); q === nothing ? "" : EpicycleBase.unit(q))

# Component names: a burn's axes (VNB → V, N, B), X Y Z for a Cartesian vector, else 1, 2, 3.
function _component_names(calc, n::Int)
    n == 1 && return [""]
    q = _quantity(calc)
    if q !== nothing && !isempty(calc.args) && hasproperty(calc.args[1], :axes) && n == 3
        a = string(nameof(typeof(getproperty(calc.args[1], :axes))))
        a == "Inertial" && return ["X", "Y", "Z"]
        (length(a) == 3 && all(isuppercase, a)) && return [string(c) for c in a]
    end
    if q !== nothing && n == 3 && EpicycleBase.unit(q) in ("km", "km/s")
        return ["X", "Y", "Z"]
    end
    return [string(i) for i in 1:n]
end

# Names that repeat get a running number, so two "Position magnitude" goals read #1 and #2.
function _number_repeats(names::Vector{String})
    seen = Dict{String, Int}()
    total = Dict{String, Int}()
    for n in names
        total[n] = get(total, n, 0) + 1
    end
    return map(names) do n
        total[n] == 1 && return n
        seen[n] = get(seen, n, 0) + 1
        return "$n #$(seen[n])"
    end
end

struct _VarSlot
    name::String
    unit::String
    index::Int          # position in the solver's variable vector
end

struct _GoalSlot
    name::String
    unit::String
    index::Int          # position in the solver's constraint vector
    lo::Float64
    hi::Float64
end

"""Rows for an event-graph targeting sequence: every `Vary` component and every `Constraint`."""
struct EventGraphRows <: MonitorRows
    vars::Vector{_VarSlot}
    goals::Vector{_GoalSlot}
    ctol::Float64
end

monitor_title(::EventGraphRows) = "Targeter"

function EventGraphRows(sm, ctol::Real)
    bases = _number_repeats([isempty(v.name) ? _quantity_label(v.calc) : v.name
                             for v in sm.ordered_vars])
    vars = _VarSlot[]
    k = 0
    for (v, base) in zip(sm.ordered_vars, bases)
        comps = _component_names(v.calc, v.numvars)
        u = _quantity_unit(v.calc)
        for i in 1:v.numvars
            k += 1
            v.lower_bound[i] == v.upper_bound[i] && continue      # held fixed: nothing to watch
            push!(vars, _VarSlot(isempty(comps[i]) ? base : "$base[$(comps[i])]", u, k))
        end
    end

    cons = [c for e in sm.sorted_events for c in e.funcs]
    bases = _number_repeats([isempty(c.name) ? _quantity_label(c.calc) : c.name for c in cons])
    goals = _GoalSlot[]
    k = 0
    for (c, base) in zip(cons, bases)
        comps = _component_names(c.calc, c.numvars)
        u = _quantity_unit(c.calc)
        for i in 1:c.numvars
            k += 1
            push!(goals, _GoalSlot(isempty(comps[i]) ? base : "$base[$(comps[i])]", u, k,
                                   c.lower_bound[i], c.upper_bound[i]))
        end
    end
    return EventGraphRows(vars, goals, Float64(ctol))
end

function monitor_rows(r::EventGraphRows, x, g)
    vars = map(r.vars) do s
        val, u = _shown(x[s.index], s.unit)
        MonitorVar(s.name, u, val)
    end
    worst = 0.0
    goals = map(r.goals) do s
        v = _violation(g[s.index], s.lo, s.hi)
        worst = max(worst, v)
        achieved, u = _shown(g[s.index], s.unit)
        err, _ = _shown(v, s.unit)
        MonitorGoal(s.name, _desired_text(s.lo, s.hi, s.unit), _fmt(achieved), u, err,
                    v <= r.ctol)
    end
    return vars, String[], goals, worst
end

# ── The box ──────────────────────────────────────────────────────────────────────────────────────

const _MAX_INNER = 78                       # inner width at most 78, so the box is at most 80
const _MIN_COLS  = 64                       # below this, one status line instead of the box
const _SPIN      = collect("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏")
const _MAX_VARS  = 8                        # rows beyond these collapse into "+ N more"
const _MAX_GOALS = 8

_inner_width(cols) = cols >= _MIN_COLS ? min(_MAX_INNER, cols - 2) : 0

mutable struct StatusBox
    W::Int                                  # inner width; 0 means the one-line form
    live::Bool                              # interactive terminal: redraw in place
    nlines::Int                             # lines the last frame drew
    top::Float64                            # the bar's left end: the first frame's measure
    last::Union{Nothing, IterationRecord}
end

StatusBox() = StatusBox(_inner_width(_term_cols()), _term_live(), 0, NaN, nothing)

_name_width(box) = clamp(box.W - 50, 12, 28)
_bar_width(box)  = clamp(box.W - 46, 10, 32)
_box_row(box, s) = "│ " * rpad(_clip(s, box.W - 2), box.W - 2) * " │"

_clip(s, n) = textwidth(s) <= n ? s : first(s, n - 1) * "…"
_fit(s, n)  = rpad(_clip(s, n), n)

# A fixed-width number that switches to exponent form when plain decimals would not fit.
function _fmt(v::Real; w::Int = 13)
    a = abs(v)
    s = (a == 0 || 1e-3 <= a < 1e7) ? @sprintf("%.6f", v) : @sprintf("%.5e", v)
    textwidth(s) > w && (s = @sprintf("%.4e", v))
    return lpad(s, w)
end

_rcol(s, w) = lpad(_clip(s, w), w)

# Log-scale progress from the first frame's measure down to the tolerance.
function _bar(err, top, tol, n)
    tol = max(tol, floatmin())
    top = max(top, 10tol)
    f = (log10(top) - log10(max(err, tol))) / (log10(top) - log10(tol))
    k = round(Int, n * clamp(f, 0, 1))
    return "█"^k * "░"^(n - k)
end

function _state_text(rec)
    rec.status === :running && return "$(_SPIN[mod1(rec.iter + 1, length(_SPIN))]) running"
    _converged(rec.status) && return "✓ converged"
    return "✗ " * replace(string(rec.status), "_" => " ")
end

function _frame(box::StatusBox, rec::IterationRecord)
    nw   = _name_width(box)
    head = "─ $(rec.solver) · $(rec.problem) ─── iter $(lpad(rec.iter, 3)) ─── $(_state_text(rec)) ─── " *
           @sprintf("%.1f s ", rec.elapsed)
    lines = ["┌" * rpad(_clip(head, box.W), box.W, '─') * "┐"]

    push!(lines, _box_row(box, @sprintf("%-*s %13s  %s", nw, "Variables", "value", "unit")))
    for v in rec.vars[1:min(end, _MAX_VARS)]
        push!(lines, _box_row(box, "  " * _fit(v.name, nw - 2) * " " * _fmt(v.value) * "  " * v.unit))
    end
    length(rec.vars) > _MAX_VARS &&
        push!(lines, _box_row(box, "  + $(length(rec.vars) - _MAX_VARS) more"))
    for n in rec.notes
        push!(lines, _box_row(box, "  " * n))
    end
    push!(lines, _box_row(box, ""))

    if !isempty(rec.goals)
        push!(lines, _box_row(box, @sprintf("%-*s %13s %13s  %-6s %9s", nw, "Goals", "desired",
                                            "achieved", "unit", "error")))
        shown = length(rec.goals) <= _MAX_GOALS ? rec.goals :
                sort(rec.goals; by = g -> g.met ? -Inf : g.error, rev = true)[1:_MAX_GOALS]
        for g in shown
            push!(lines, _box_row(box, "  " * _fit(g.name, nw - 2) * " " * _rcol(g.desired, 13) *
                                       " " * _rcol(g.achieved, 13) * "  " * _fit(g.unit, 6) *
                                       @sprintf(" %9.1e ", g.error) * (g.met ? "✓" : "✗")))
        end
        length(rec.goals) > _MAX_GOALS &&
            push!(lines, _box_row(box, "  + $(length(rec.goals) - _MAX_GOALS) more"))
        push!(lines, _box_row(box, ""))
    end

    if !isnan(rec.objective)
        push!(lines, _box_row(box, @sprintf("%-*s %s   inf_du %.1e   μ %.1e", nw, "Objective",
                                            _fmt(rec.objective), rec.inf_du, rec.mu)))
    end
    met = count(g -> g.met, rec.goals)
    lead = isempty(rec.goals) ? "Optimality      " : @sprintf("Goals met %d/%d    ", met, length(rec.goals))
    push!(lines, _box_row(box, lead * _bar(rec.worst, box.top, rec.tol, _bar_width(box)) *
                               @sprintf("  %.1e", rec.worst)))
    push!(lines, "└" * "─"^box.W * "┘")
    return lines
end

function _status_line(rec::IterationRecord, cols)
    met = count(g -> g.met, rec.goals)
    s = isempty(rec.goals) ?
        @sprintf("%s · iter %d · optimality %.1e", _state_text(rec), rec.iter, rec.worst) :
        @sprintf("%s · iter %d · goals %d/%d · err %.1e", _state_text(rec), rec.iter, met,
                 length(rec.goals), rec.worst)
    return _clip(s, cols - 1)
end

# Remove the last frame: move up over exactly the lines it drew and clear below.
function _erase!(box::StatusBox)
    box.nlines > 0 && print(stdout, "\e[$(box.nlines)A\e[J")
    box.nlines = 0
    return nothing
end

# Before a redraw. Wider terminal: the old frame still fits, so clear it and redraw in place.
# Narrower: the terminal re-wrapped the old frame and its line count can't be trusted, so leave it
# and start a fresh frame below.
function _clear_for_redraw!(box::StatusBox)
    W = _inner_width(_term_cols())
    W >= box.W ? _erase!(box) : (box.nlines = 0)
    box.W = W
    return nothing
end

function _draw!(box::StatusBox, rec::IterationRecord)
    lines = box.W == 0 ? [_status_line(rec, _term_cols())] : _frame(box, rec)
    foreach(l -> println(stdout, l), lines)
    flush(stdout)
    box.nlines = length(lines)
    return nothing
end

"""
    show_iteration!(box::StatusBox, rec::IterationRecord) -> Bool

Show one iteration. In an interactive terminal the box is redrawn in place; elsewhere nothing is
printed until the final record. Returns `true`; a solver that reads `false` would stop.
"""
function show_iteration!(box::StatusBox, rec::IterationRecord)
    isnan(box.top) && (box.top = rec.worst)
    box.last = rec
    if box.live
        _clear_for_redraw!(box)
        _draw!(box, rec)
    end
    return true
end

"""
    finish!(box::StatusBox, status::Symbol)

Draw the final frame with the solver's exit status, and a one-line summary below it.
"""
function finish!(box::StatusBox, status::Symbol)
    rec = box.last
    rec === nothing && return nothing
    final = IterationRecord(rec.solver, rec.problem, rec.iter, rec.vars, rec.notes, rec.goals,
                            rec.objective, rec.inf_pr, rec.inf_du, rec.mu, rec.worst, rec.tol,
                            rec.elapsed, status)
    box.last = final
    box.live ? _clear_for_redraw!(box) : (box.W = _inner_width(_term_cols()))
    _draw!(box, final)
    verb = _converged(status) ? "Converged in" : "Stopped after"
    @printf(stdout, "%s %d iterations, %.1f s (%s)\n", verb, final.iter, final.elapsed, status)
    box.nlines = 0                          # the final frame stays; nothing redraws over it
    return nothing
end

# A message during a solve: clear the box, print the message, redraw the box below it.
function _note!(f, box::StatusBox)
    live = box.live && box.nlines > 0
    live && _clear_for_redraw!(box)
    f()
    flush(stdout); flush(stderr)
    live && box.last !== nothing && _draw!(box, box.last)
    return nothing
end

# Logging during a solve goes through the box, so a warning lands above it rather than being
# drawn over.
struct _MonitorLogger{L<:AbstractLogger} <: AbstractLogger
    box::StatusBox
    inner::L
end

min_enabled_level(l::_MonitorLogger) = min_enabled_level(l.inner)
shouldlog(l::_MonitorLogger, args...) = shouldlog(l.inner, args...)
catch_exceptions(l::_MonitorLogger) = catch_exceptions(l.inner)
handle_message(l::_MonitorLogger, args...; kwargs...) =
    _note!(() -> handle_message(l.inner, args...; kwargs...), l.box)

"""Run `f` with log messages routed above `box`."""
with_box_logging(f, box::StatusBox) = with_logger(f, _MonitorLogger(box, current_logger()))
