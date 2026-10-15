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
One row standing for a group of constraint rows, for a problem too large to show a row per goal:
one kind of constraint in one phase, such as a phase's defects. `met` of `total` rows are within
tolerance; `error` is the worst violation among them, and `where` says which row that is.
"""
struct MonitorGroup
    name::String
    met::Int
    total::Int
    where::String
    error::Float64
end

_is_met(g::MonitorGoal)  = g.met
_is_met(g::MonitorGroup) = g.met == g.total

"""
    IterationRecord

What one solver iteration reports to the status box. Every solver adapter produces this and the
display reads nothing else, so a new solver needs an adapter and no change to the display. The
goals are either one row per goal or one row per group of goals.
"""
struct IterationRecord
    solver::String
    problem::String
    iter::Int
    vars::Vector{MonitorVar}
    notes::Vector{String}          # summary lines under Variables, e.g. a decision-variable count
    goals::Union{Vector{MonitorGoal}, Vector{MonitorGroup}}
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
`(vars, notes, goals, worst)`, where `goals` holds `MonitorGoal`s or `MonitorGroup`s;
`monitor_title(rows)` names the problem in the box's header.
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
    scale::Float64      # the solver sees (value - shift) / scale
    shift::Float64
end

_VarSlot(name, unit, index) = _VarSlot(name, unit, index, 1.0, 0.0)

# The value the user wrote, from the solver's scaled one.
_physical(s::_VarSlot, x) = x[s.index] * s.scale + s.shift

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

# ── Rows for a phased problem: collocation and shooting ──────────────────────────────────────────
#
# Hundreds or thousands of decision variables and constraint rows, so the box summarizes.
# Variables: the free scalars (times, masses, parameters, v-infinities), a line per phase and the
# total count. Goals: one group per phase and kind of constraint (defects or the match point,
# path, boundary) and one per pair of linked phases, each with its met count and where its worst
# row is. `CollocationRows` and `ShootingRows` build the same `GroupedRows` from their managers.

# A run of constraint rows from one constraint, and what is needed to say which row is which.
struct _RowBlock
    rows::UnitRange{Int}                    # positions in the solver's constraint vector
    kind::Symbol                            # :defect, :path, :segment_path, :match_point,
                                            # :boundary or :linkage
    label::String                           # the constraint's name; empty for defects
    comps::Vector{String}                   # component names; "" for a scalar constraint
    phase::Int                              # the phase whose node times place a row, 0 for none
    steps::Vector{Tuple{UnitRange{Int}, Int}}   # defects: each step's rows in the block, first node
    segments::NTuple{2, Int}                # :segment_path: forward and backward segment counts
end

_RowBlock(rows, kind, label, comps, phase, steps) =
    _RowBlock(rows, kind, label, comps, phase, steps, (0, 0))

struct _RowGroup
    name::String
    blocks::Vector{_RowBlock}
end

# Where a phase's initial and final times sit in the decision vector, and its mesh.
struct _PhaseClock
    t0::Int
    tf::Int
    mesh::Any
end

"""Rows for a phased problem: free scalars, a line per phase, grouped goals."""
struct GroupedRows <: MonitorRows
    vars::Vector{_VarSlot}
    notes::Vector{String}
    groups::Vector{_RowGroup}
    clocks::Vector{_PhaseClock}             # by phase index; empty where rows carry no times
    lg::Vector{Float64}
    ug::Vector{Float64}
    ctol::Float64
end

# A `Link` written with a name (`Link(a, b; name = :earth_flyby)`), keyed by the two phases, so a
# link's group is called what the script calls it. Filled when the links are registered.
const _LINK_NAMES = Dict{Tuple{UInt64, UInt64}, String}()

_link_key(p1, p2) = (objectid(p1), objectid(p2))

# Add a link's rows to the group for its pair of phases, after the phases' own groups.
function _add_link_block!(groups, index, p1, p2, block)
    key = _link_key(p1, p2)
    if haskey(index, key)
        push!(groups[index[key]].blocks, block)
    else
        name = get(_LINK_NAMES, key, "$(p1.name)→$(p2.name)")
        push!(groups, _RowGroup("$name · link", [block]))
        index[key] = length(groups)
    end
    return groups
end

# One line per phase, the first few when there are many, then the decision-variable count.
function _phase_notes(lines::Vector{String}, nvars::Int)
    if length(lines) > _MAX_PHASE_LINES
        lines = [lines[1:_MAX_PHASE_LINES-1]; "+ $(length(lines) - _MAX_PHASE_LINES + 1) more phases"]
    end
    return [lines; "$nvars decision variables"]
end

const _MAX_PHASE_LINES = 4                  # more phases than this: the first few and a count

_comp_names(n::Int) = n == 1 ? [""] : [string(i) for i in 1:n]
_with_comp(label, comp) = isempty(comp) ? label : "$label[$comp]"
_generated(name) = isempty(name) || startswith(name, "##")

# The name a function was written with, or `nothing` for an anonymous one. A function defined
# inside another carries a suffix (`accel#12`), which is dropped.
function _readable_fname(f)
    s = first(split(lstrip(string(nameof(typeof(f))), '#'), '#'))
    return isempty(s) || isdigit(first(s)) ? nothing : String(s)
end

# A state's components by the field names of the phase's state type, or by index.
function _state_names(p)
    reg = get(_phase_registry, objectid(p), nothing)
    names = try
        reg === nothing || reg.state_type === nothing ? Symbol[] : collect(fieldnames(reg.state_type))
    catch
        Symbol[]
    end
    return length(names) == p._n_states ? string.(names) : _comp_names(p._n_states)
end

# The names of a phase's boundary or path constraints, in the order the phase holds them. A
# constraint written without `name` carries a generated one (`##con#280`), so it is named by the
# function the user wrote instead, plus the end it applies at. Constraints and the specs recorded
# as they were written come in the same order; when the counts differ the position can't be
# trusted, and the constraint is named by its place.
function _constraint_labels(p, path::Bool)
    cons  = path ? [c for c in p.path_constraints if c isa PathConstraint] :
                   [c for c in p.constraints if c isa BoundaryConstraint]
    specs = [sp for sp in phase_specs(p) if sp.kind === :constraint && (sp.at isa Path) == path]
    return map(enumerate(cons)) do (k, c)
        name = path ? c.name : c.calc.name
        _generated(name) || return name
        length(specs) == length(cons) || return "constraint $k"
        f = _readable_fname(specs[k].f)
        f === nothing && return "constraint $k"
        at = specs[k].at
        return at isa Initial ? "$f (initial)" : at isa Final ? "$f (final)" : f
    end
end

"""Rows for a collocation sequence, from its `CollocationManager`."""
function CollocationRows(cm, ctol::Real)
    multi = length(cm.phases) > 1
    tag(p, what) = multi ? "$(p.name) · $what" : what
    vars, notes = _VarSlot[], String[]
    groups, clocks = _RowGroup[], _PhaseClock[]

    for (i, p) in enumerate(cm.phases)
        col = cm.phase_var_offsets[i]
        pre = multi ? "$(p.name) · " : ""
        t0 = tf = 0
        for (v, r) in zip(variable_list(p), variable_ranges(p))
            v === p.t0_var && (t0 = col + first(r))
            v === p.tf_var && (tf = col + first(r))
            base = v === p.t0_var ? "initial_time" : v === p.tf_var ? "final_time" :
                   v === p.param_var ? "parameter" : continue
            # "t0", "tf" and "param" are the framework's own names; a `Vary(...; name)` replaces them.
            base = v.name in ("", "t0", "tf", "param") ? pre * base : v.name
            for (k, j) in enumerate(r)
                v.lower_bounds[k] == v.upper_bounds[k] && continue      # held fixed
                push!(vars, _VarSlot(length(r) == 1 ? base : "$base[$k]", "", col + j))
            end
        end
        push!(clocks, _PhaseClock(t0, tf, build_mesh(p.transcription)))

        ns     = p._n_states
        steps  = [(rows, first(nodes)) for (rows, nodes) in defect_row_groups(p.transcription, ns)]
        blabel = _constraint_labels(p, false)
        plabel = _constraint_labels(p, true)
        defect, path, bound = _RowBlock[], _RowBlock[], _RowBlock[]
        row = cm.phase_fun_offsets[i]
        for pf in function_list(p)
            r = row+1 : row+pf.n_nlp
            row += pf.n_nlp
            isempty(r) && continue
            if pf.source isa DefectBlock
                push!(defect, _RowBlock(r, :defect, "", _state_names(p), i, steps))
            elseif pf.source isa PathConstraintBlock
                push!(path, _RowBlock(r, :path, plabel[length(path) + 1],
                                      _comp_names(pf.source.pc.n_pc), i, []))
            else
                push!(bound, _RowBlock(r, :boundary, blabel[length(bound) + 1],
                                       _comp_names(pf.n_nlp), i, []))
            end
        end
        isempty(defect) || push!(groups, _RowGroup(tag(p, "defects"), defect))
        isempty(path)   || push!(groups, _RowGroup(tag(p, "path"), path))
        isempty(bound)  || push!(groups, _RowGroup(tag(p, "boundary"), bound))
    end

    # One group per pair of linked phases, after the phases.
    links = Dict{Tuple{UInt64, UInt64}, Int}()
    for (k, lc) in enumerate(cm.linkages)
        n = length(lc.lower_bounds)
        n == 0 && continue
        r = cm.linkage_fun_offsets[k]+1 : cm.linkage_fun_offsets[k]+n
        label = _generated(lc.name) ? "link" : lc.name
        comps = label == "continuity" && n == lc.phase1._n_states + 1 ?
                [_state_names(lc.phase1); "t"] : _comp_names(n)
        _add_link_block!(groups, links, lc.phase1, lc.phase2,
                         _RowBlock(r, :linkage, label, comps, 0, []))
    end

    lines = ["$(p.name) · $(nameof(typeof(p.transcription))), $(_count(n_control_nodes(p), "node"))"
             for p in cm.phases]
    notes = _phase_notes(lines, cm.n_vars)

    lg, ug = get_constraint_bounds(cm)
    return GroupedRows(vars, notes, groups, clocks, lg, ug, Float64(ctol))
end

# ── Shooting: Sims-Flanagan and MGA phases ───────────────────────────────────────────────────────

# The components of a shooting phase's match-point defect: position, velocity, then mass.
const _MATCH_POINT_COMPS = ["x", "y", "z", "vx", "vy", "vz", "m"]

# The quantities a shooting phase's scalar and three-component variables stand for, by the field
# that holds each. The large blocks (throttles, deep-space burns, arc fractions, endpoint states)
# are counted in the total rather than listed.
const _SHOOTING_SCALARS = (:t0_var => "initial_time", :tf_var => "final_time",
                           :m0_var => "initial_mass", :mf_var => "final_mass",
                           :vinf_dep_var => "departure_vinf", :vinf_arr_var => "arrival_vinf")

# A shooting phase's boundary constraints, which also hold its path constraints lowered to one
# row per segment: each with its name and whether it was written `at = Path()`. Named as for
# collocation when the user gave no name.
function _shooting_constraint_labels(p)
    cons  = [c for c in _constraints_of(p) if c isa BoundaryConstraint]
    specs = [sp for sp in phase_specs(p) if sp.kind === :constraint]
    trusted = length(specs) == length(cons)
    return map(enumerate(cons)) do (k, c)
        path = trusted && specs[k].at isa Path
        _generated(c.calc.name) || return (c.calc.name, path)
        trusted || return ("constraint $k", false)
        f = _readable_fname(specs[k].f)
        f === nothing && return ("constraint $k", path)
        at = specs[k].at
        return (at isa Initial ? "$f (initial)" : at isa Final ? "$f (final)" : f, path)
    end
end

# The phase line under Variables: its transcription and how finely it is divided. The phase types
# are defined after this file, so they are tested here rather than dispatched on.
function _shooting_phase_line(p)
    p isa SimsFlanaganPhase && return "$(p.name) · SimsFlanagan, $(_count(n_fwd(p) + n_bwd(p), "segment"))"
    p isa MGAnDSMsPhase     && return "$(p.name) · MGAnDSMs, $(_count(n_dsm(p), "deep-space maneuver"))"
    return "$(p.name) · $(nameof(typeof(p)))"
end

# A Sims-Flanagan phase's forward and backward segment counts; (0, 0) for any other phase.
_segments(p) = p isa SimsFlanaganPhase ? (n_fwd(p), n_bwd(p)) : (0, 0)

"""Rows for a shooting sequence (Sims-Flanagan and MGA phases), from its `ShootingManager`."""
function ShootingRows(sm, ctol::Real)
    multi = length(sm.phases) > 1
    tag(p, what) = multi ? "$(p.name) · $what" : what
    vars, groups = _VarSlot[], _RowGroup[]

    # A variable shared by two phases is one column; it is listed once, under the first.
    listed = Set{UInt64}()
    for p in sm.phases, (field, quantity) in _SHOOTING_SCALARS
        hasproperty(p, field) || continue
        v = getproperty(p, field)
        (v === nothing || objectid(v) in listed) && continue
        push!(listed, objectid(v))
        cols = sm.global_var_ranges[objectid(v)]
        base = v.name == quantity || _generated(v.name) ?
               (multi ? "$(p.name) · $quantity" : quantity) : v.name
        for (k, j) in enumerate(cols)
            v.lower_bounds[k] == v.upper_bounds[k] && continue      # held fixed
            push!(vars, _VarSlot(length(cols) == 1 ? base : "$base[$k]", "", j,
                                 v.scale[k], v.shift[k]))
        end
    end

    for (i, p) in enumerate(sm.phases)
        labels = _shooting_constraint_labels(p)
        match, path, bound = _RowBlock[], _RowBlock[], _RowBlock[]
        row = sm.phase_fun_offsets[i] - 1                         # offsets here are 1-based
        nb  = 0
        for pf in function_list(p)
            r = row+1 : row+pf.n_nlp
            row += pf.n_nlp
            isempty(r) && continue
            if pf.source isa BoundaryConstraint
                nb += 1
                label, at_path = labels[nb]
                segs = _segments(p)
                if at_path && sum(segs) > 0 && pf.n_nlp % sum(segs) == 0
                    push!(path, _RowBlock(r, :segment_path, label,
                                          _comp_names(pf.n_nlp ÷ sum(segs)), 0, [], segs))
                else
                    push!(bound, _RowBlock(r, :boundary, label, _comp_names(pf.n_nlp), 0, []))
                end
            elseif pf.n_nlp == length(_MATCH_POINT_COMPS) && pf.name == "matchpoint_defect"
                push!(match, _RowBlock(r, :match_point, "", _MATCH_POINT_COMPS, 0, []))
            else                                                  # e.g. MGA's arc-fraction sum
                push!(bound, _RowBlock(r, :boundary, replace(pf.name, "_" => " "),
                                       _comp_names(pf.n_nlp), 0, []))
            end
        end
        row == sm.phase_fun_offsets[i] - 1 + n_constraints(p) || throw(ArgumentError(
            "phase $(p.name): its function list covers $(row - sm.phase_fun_offsets[i] + 1) " *
            "rows, but it reports $(n_constraints(p))"))
        isempty(match) || push!(groups, _RowGroup(tag(p, "match point"), match))
        isempty(path)  || push!(groups, _RowGroup(tag(p, "path"), path))
        isempty(bound) || push!(groups, _RowGroup(tag(p, "boundary"), bound))
    end

    links = Dict{Tuple{UInt64, UInt64}, Int}()
    for (k, sc) in enumerate(sm.seq_constraints)
        n = length(sc.lower_bounds)
        n == 0 && continue
        r = sm.seq_con_offsets[k] : sm.seq_con_offsets[k] + n - 1
        label = _generated(sc.name) ? "link" : sc.name
        _add_link_block!(groups, links, sc.phase1, sc.phase2,
                         _RowBlock(r, :linkage, label, _comp_names(n), 0, []))
    end

    notes = _phase_notes([_shooting_phase_line(p) for p in sm.phases], sm.n_vars)
    lg, ug = get_constraint_bounds(sm)
    return GroupedRows(vars, notes, groups, _PhaseClock[], lg, ug, Float64(ctol))
end

_tfmt(t) = @sprintf("%.4g", t)

# The time of node `k` of phase `i` at the current iterate.
function _node_time(r::GroupedRows, i::Int, k::Int, x)
    c = r.clocks[i]
    (c.t0 == 0 || c.tf == 0) && return nothing
    ts = node_times(c.mesh, Float64(x[c.t0]), Float64(x[c.tf]))
    return 1 <= k <= length(ts) ? ts[k] : nothing
end

# The first of these that fits the "worst at" column, or the last one cut to fit.
function _fit_where(texts...)
    for s in texts
        textwidth(s) <= _WHERE_W && return s
    end
    return _clip(last(texts), _WHERE_W)
end

# Which row of `b` the solver-vector row `i` is: a component, and a step or node with its time.
function _where(r::GroupedRows, b::_RowBlock, i::Int, x)
    l = i - first(b.rows) + 1
    if b.kind === :defect
        ns   = length(b.comps)
        comp = ns > 0 && length(b.rows) % ns == 0 ? b.comps[mod1(l, ns)] : ""
        k    = findfirst(s -> l in s[1], b.steps)
        k === nothing && return isempty(comp) ? "row $l" : comp
        t    = _node_time(r, b.phase, b.steps[k][2], x)
        at   = "$(isempty(comp) ? "" : comp * " ")@ step $k"
        return t === nothing ? at : _fit_where("$at (t $(_tfmt(t)))", at)
    elseif b.kind === :path
        n    = length(b.comps)
        node = (l - 1) ÷ n + 1
        name = _with_comp(b.label, b.comps[mod1(l, n)])
        t    = _node_time(r, b.phase, node, x)
        at   = "$name @ node $node"
        return t === nothing ? _fit_where(at) :
               _fit_where("$at (t $(_tfmt(t)))", at, "node $node (t $(_tfmt(t)))")
    elseif b.kind === :segment_path
        # One row per component per segment: the forward half's segments, then the backward's.
        n    = length(b.comps)
        seg  = (l - 1) ÷ n + 1
        nf   = b.segments[1]
        name = _with_comp(b.label, b.comps[mod1(l, n)])
        half, k = seg <= nf ? ("fwd", seg) : ("bwd", seg - nf)
        return _fit_where("$name @ $half seg $k", "$half seg $k")
    elseif b.kind === :match_point
        return b.comps[mod1(l, length(b.comps))]
    else
        return _fit_where(_with_comp(b.label, b.comps[min(l, end)]))
    end
end

# A NaN constraint is as unmet as it gets.
_violation_or_inf(g, lo, hi) = isnan(g) ? Inf : _violation(g, lo, hi)

function monitor_rows(r::GroupedRows, x, g)
    vars  = [MonitorVar(s.name, s.unit, _physical(s, x)) for s in r.vars]
    worst = 0.0
    groups = map(r.groups) do grp
        met, total, werr = 0, 0, 0.0
        wb, wi = first(grp.blocks), 0
        for b in grp.blocks, i in b.rows
            v = _violation_or_inf(g[i], r.lg[i], r.ug[i])
            total += 1
            v <= r.ctol && (met += 1)
            if v > werr
                werr, wb, wi = v, b, i
            end
        end
        worst = max(worst, werr)
        where = met < total && wi > 0 ? _where(r, wb, wi, x) : ""
        MonitorGroup(grp.name, met, total, where, werr)
    end
    return vars, r.notes, groups, worst
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

# A row name of the form "phase · what" that is too long loses characters from the phase, not
# from what follows it: "earth_to_fl… · match point", since the kind or quantity is what tells
# rows apart. Without the separator, or when even that will not fit, the end is cut as usual.
const _SEP = " · "

function _fit_name(s, n)
    textwidth(s) <= n && return rpad(s, n)
    i = findlast(_SEP, s)
    if i !== nothing
        head, tail = s[1:prevind(s, first(i))], s[nextind(s, last(i)):end]
        room = n - textwidth(_SEP) - textwidth(tail)
        room >= 4 && return rpad(_clip(head, room) * _SEP * tail, n)
    end
    return _fit(s, n)
end

_count(n, noun) = "$n $noun" * (n == 1 ? "" : "s")

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

# One row per goal: the desired and achieved values, unit and error. More than fit show the unmet
# goals with the largest errors first.
function _goal_lines!(lines, box, nw, goals::Vector{MonitorGoal})
    push!(lines, _box_row(box, @sprintf("%-*s %13s %13s  %-6s %9s", nw, "Goals", "desired",
                                        "achieved", "unit", "error")))
    shown = length(goals) <= _MAX_GOALS ? goals :
            sort(goals; by = g -> g.met ? -Inf : g.error, rev = true)[1:_MAX_GOALS]
    for g in shown
        push!(lines, _box_row(box, "  " * _fit_name(g.name, nw - 2) * " " * _rcol(g.desired, 13) *
                                   " " * _rcol(g.achieved, 13) * "  " * _fit(g.unit, 6) *
                                   @sprintf(" %9.1e ", g.error) * (g.met ? "✓" : "✗")))
    end
    length(goals) > _MAX_GOALS &&
        push!(lines, _box_row(box, "  + $(length(goals) - _MAX_GOALS) more"))
    return lines
end

const _MET_W   = 11                         # "4850/4872" and wider
const _WHERE_W = 22                         # the "worst at" column; adapters keep within it

# One row per group: how many of its rows are met, where the worst one is and its error. Up to
# `_MAX_GOALS` groups show in the order the adapter gives, which is trajectory order. More than
# that show the unmet groups first, then the largest errors, and a last line counts the rest.
function _goal_lines!(lines, box, nw, groups::Vector{MonitorGroup})
    push!(lines, _box_row(box, @sprintf("%-*s %*s  %-*s %9s", nw, "Goals", _MET_W, "met",
                                        _WHERE_W, "worst at", "error")))
    shown, hidden = groups, MonitorGroup[]
    if length(groups) > _MAX_GOALS
        order  = sort(groups; by = g -> (_is_met(g), -g.error))
        shown  = order[1:_MAX_GOALS]
        hidden = order[_MAX_GOALS+1:end]
    end
    for g in shown
        ok = _is_met(g)
        push!(lines, _box_row(box, "  " * _fit_name(g.name, nw - 2) * " " *
                                   _rcol("$(g.met)/$(g.total)", _MET_W) * "  " *
                                   _fit(ok ? "" : g.where, _WHERE_W) *
                                   @sprintf(" %9.1e ", g.error) * (ok ? "✓" : "✗")))
    end
    if !isempty(hidden)
        unmet = count(!_is_met, hidden)
        push!(lines, _box_row(box, "  + $(length(hidden)) more · " *
                                   (unmet == 0 ? "all met" : "$unmet unmet")))
    end
    return lines
end

function _frame(box::StatusBox, rec::IterationRecord)
    nw   = _name_width(box)
    head = "─ $(rec.solver) · $(rec.problem) ─── iter $(lpad(rec.iter, 3)) ─── $(_state_text(rec)) ─── " *
           @sprintf("%.1f s ", rec.elapsed)
    lines = ["┌" * rpad(_clip(head, box.W), box.W, '─') * "┐"]

    push!(lines, _box_row(box, @sprintf("%-*s %13s  %s", nw, "Variables", "value", "unit")))
    for v in rec.vars[1:min(end, _MAX_VARS)]
        push!(lines, _box_row(box, "  " * _fit_name(v.name, nw - 2) * " " * _fmt(v.value) * "  " * v.unit))
    end
    length(rec.vars) > _MAX_VARS &&
        push!(lines, _box_row(box, "  + $(length(rec.vars) - _MAX_VARS) more"))
    for n in rec.notes
        push!(lines, _box_row(box, "  " * n))
    end
    push!(lines, _box_row(box, ""))

    if !isempty(rec.goals)
        _goal_lines!(lines, box, nw, rec.goals)
        push!(lines, _box_row(box, ""))
    end

    if !isnan(rec.objective)
        push!(lines, _box_row(box, @sprintf("%-*s %s   inf_du %.1e   μ %.1e", nw, "Objective",
                                            _fmt(rec.objective), rec.inf_du, rec.mu)))
    end
    met = count(_is_met, rec.goals)
    lead = isempty(rec.goals) ? "Optimality      " : @sprintf("Goals met %d/%d    ", met, length(rec.goals))
    push!(lines, _box_row(box, lead * _bar(rec.worst, box.top, rec.tol, _bar_width(box)) *
                               @sprintf("  %.1e", rec.worst)))
    push!(lines, "└" * "─"^box.W * "┘")
    return lines
end

function _status_line(rec::IterationRecord, cols)
    met = count(_is_met, rec.goals)
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
