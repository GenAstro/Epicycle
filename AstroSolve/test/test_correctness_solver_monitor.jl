# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# The live solver status box and the IPOPT driver behind it.
#
# The box is checked through what it prints: widths at several terminal sizes, the in-place
# redraw (the cursor-up count must equal the lines the last frame drew, so nothing earlier is ever
# erased), the resize rules, messages landing above the box, and final-frame-only output when not
# a terminal. The driver is checked against SNOW's own IPOPT driver on a small problem, then on a
# targeting sequence end to end. Terminal width and "is a terminal" are injected, since a test's
# stdout is neither.

using Test
using Printf
using AstroProp: Vern9

const _SM = AstroSolve

# Run f with stdout captured, a fixed terminal width, and live redraw on or off.
function _sm_capture(f; cols = 100, live = true)
    _SM._MONITOR_COLS[] = cols
    _SM._MONITOR_LIVE[] = live
    path, io = mktemp()
    try
        try
            redirect_stdout(f, io)
        finally
            close(io)                    # before reading, and before rm on Windows
        end
        return read(path, String)
    finally
        _SM._MONITOR_COLS[] = nothing
        _SM._MONITOR_LIVE[] = nothing
        rm(path; force = true)
    end
end

_sm_strip(s) = replace(s, r"\e\[\d+A" => "", "\e[J" => "")
_sm_box_lines(s) = [l for l in split(_sm_strip(s), '\n') if startswith(l, r"[┌│└]")]
_sm_ups(s) = [parse(Int, m.captures[1]) for m in eachmatch(r"\e\[(\d+)A", s)]

function _sm_record(iter; nvars = 4, ngoals = 4, status = :running, worst = 10.0^(-iter),
                    objective = NaN, longname = false)
    vars  = [_SM.MonitorVar(longname && i == 1 ? "a variable name far too long for any column" :
                            "toi.delta_v[$i]", "km/s", 2.5 + i) for i in 1:nvars]
    goals = [_SM.MonitorGoal("Goal $i", _SM._fmt(100.0i), _SM._fmt(100.0i + worst), "km",
                             worst, worst <= 1e-6) for i in 1:ngoals]
    return _SM.IterationRecord("IPOPT", "Targeter", iter, vars, String[], goals, objective,
                               worst, 1e-3, 1e-5, worst, 1e-6, 0.1iter, status)
end

@testset "Solver status box" begin

    @testset "fits the terminal: at most 80 columns, narrower when the terminal is" begin
        for (cols, width) in ((200, 80), (80, 80), (70, 70), (64, 64))
            out = _sm_capture(cols = cols, live = false) do
                box = _SM.StatusBox()
                _SM.show_iteration!(box, _sm_record(0; longname = true))
                _SM.finish!(box, :Solve_Succeeded)
            end
            lines = _sm_box_lines(out)
            @test !isempty(lines)
            @test all(textwidth.(lines) .== width)
        end
        # Below 64 columns: one status line that fits
        out = _sm_capture(cols = 50, live = false) do
            box = _SM.StatusBox()
            _SM.show_iteration!(box, _sm_record(3))
            _SM.finish!(box, :Solve_Succeeded)
        end
        @test isempty(_sm_box_lines(out))
        first_line = first(split(out, '\n'))
        @test occursin("iter 3", first_line) && textwidth(first_line) <= 49
    end

    @testset "names, units, the ✓/✗ marks and fixed height" begin
        out = _sm_capture(live = false) do
            box = _SM.StatusBox()
            _SM.show_iteration!(box, _sm_record(0; longname = true))
            _SM.finish!(box, :Solve_Succeeded)
        end
        @test occursin("…", out)                         # a long name is cut, not wrapped
        @test occursin("km/s", out) && occursin("km", out)
        @test occursin("✗", out)                          # worst = 1.0 at iteration 0
        @test occursin("✓ converged", out)
        @test occursin("Converged in 0 iterations", out)

        # More rows than fit collapse into "+ N more"; every frame has the same height.
        heights = map(0:3) do k
            out = _sm_capture(live = false) do
                box = _SM.StatusBox()
                _SM.show_iteration!(box, _sm_record(k; nvars = 11, ngoals = 12))
                _SM.finish!(box, :Solve_Succeeded)
            end
            @test occursin("+ 3 more", out) && occursin("+ 4 more", out)
            length(_sm_box_lines(out))
        end
        @test allequal(heights)
    end

    @testset "bounds read as bounds in the desired column" begin
        @test strip(_SM._desired_text(42166.9, 42166.9, "km")) == "42166.900000"
        @test startswith(_SM._desired_text(6700.0, Inf, "km"), "≥")
        @test startswith(_SM._desired_text(-Inf, 0.1, ""), "≤")
        @test occursin("…", _SM._desired_text(1.0, 2.0, ""))
        @test _SM._shown(π, "rad") == (180.0, "deg")
    end

    @testset "redraws in place, never over earlier output" begin
        out = _sm_capture(live = true) do
            println("earlier output")
            box = _SM.StatusBox()
            for k in 0:4
                _SM.show_iteration!(box, _sm_record(k))
            end
            _SM.finish!(box, :Solve_Succeeded)
        end
        @test startswith(out, "earlier output\n")
        frame_height = length(_sm_box_lines(_sm_capture(live = false) do
            box = _SM.StatusBox(); _SM.show_iteration!(box, _sm_record(0))
            _SM.finish!(box, :Solve_Succeeded)
        end))
        ups = _sm_ups(out)
        @test length(ups) == 5                 # every frame after the first, plus the final one
        @test all(==(frame_height), ups)       # exactly the lines the last frame drew
    end

    @testset "a narrower terminal starts a fresh frame; a wider one redraws in place" begin
        out = _sm_capture(live = true, cols = 100) do
            box = _SM.StatusBox()
            _SM.show_iteration!(box, _sm_record(0))
            _SM.show_iteration!(box, _sm_record(1))
            _SM._MONITOR_COLS[] = 66                    # narrowed
            _SM.show_iteration!(box, _sm_record(2))
            _SM.show_iteration!(box, _sm_record(3))
            _SM._MONITOR_COLS[] = 100                   # widened back
            _SM.show_iteration!(box, _sm_record(4))
        end
        frames = split(out, "┌")[2:end]
        before = [split(out, "┌")[i] for i in 1:length(frames)]
        erased = [occursin(r"\e\[\d+A\e\[J$", b) for b in before]
        widths = [textwidth("┌" * first(split(f, '\n'))) for f in frames]
        @test erased == [false, true, false, true, true]
        @test widths == [80, 80, 66, 66, 80]
    end

    @testset "messages land above the box, and the box is redrawn below them" begin
        out = _sm_capture(live = true) do
            box = _SM.StatusBox()
            _SM.show_iteration!(box, _sm_record(0))
            _SM._note!(() -> println("a message mid-solve"), box)
            _SM.show_iteration!(box, _sm_record(1))
        end
        clean = _sm_strip(out)
        @test count("a message mid-solve", clean) == 1
        msg_at = findfirst("a message mid-solve", clean)[1]
        @test occursin("┌", clean[msg_at:end])               # a box follows the message
        @test endswith(rstrip(clean), "┘")
    end

    @testset "not a terminal: only the final frame" begin
        out = _sm_capture(live = false) do
            box = _SM.StatusBox()
            for k in 0:5
                _SM.show_iteration!(box, _sm_record(k))
            end
            _SM.finish!(box, :Solve_Succeeded)
        end
        @test count("┌", out) == 1
        @test isempty(_sm_ups(out))
    end
end

# ── The IPOPT driver ─────────────────────────────────────────────────────────────────────────────

# minimize (x1 - 1)² + (x2 - 2)²  subject to  x1 + x2 = 2,  x1 ≥ 0.2
function _sm_nlp!(g, x)
    g[1] = x[1] + x[2]
    g[2] = x[1]
    return (x[1] - 1)^2 + (x[2] - 2)^2
end

function _sm_solve(solver; rows = nothing)
    opts = AstroSolve.SNOW.Options(derivatives = AstroSolve.SNOW.ForwardFD(), solver = solver)
    rows === nothing || (opts = _SM._with_rows(opts, rows))
    return AstroSolve.SNOW.minimize(_sm_nlp!, [0.0, 0.0], 2, [-10.0, -10.0], [10.0, 10.0],
                                    [2.0, 0.2], [2.0, Inf], opts)
end

_sm_opts() = Dict{String,Any}("print_level" => 0, "tol" => 1e-8,
                              "output_file" => tempname() * "_ipopt.out")

struct _SMThrowingRows <: _SM.MonitorRows
    err::Exception
end
_SM.monitor_rows(r::_SMThrowingRows, x, g) = throw(r.err)

@testset "IPOPT driver" begin

    @testset "solves as SNOW's IPOPT driver does, with the monitor on or off" begin
        xs, fs, is, _ = _sm_solve(AstroSolve.SNOW.IPOPT(_sm_opts()))
        local xo, fo, io
        out_off = _sm_capture(live = true) do
            xo, fo, io, _ = _sm_solve(_SM._ipopt(_sm_opts(); monitor = false))
        end
        local xm, fm, im
        out_on = _sm_capture(live = true) do
            xm, fm, im, _ = _sm_solve(_SM._ipopt(_sm_opts(); monitor = true))
        end
        @test is == io == im == :Solve_Succeeded
        @test isapprox(xo, xs; atol = 1e-8) && isapprox(xm, xs; atol = 1e-8)
        @test isapprox(xs, [0.5, 1.5]; atol = 1e-6)
        @test !occursin("┌", out_off)                       # monitor off: no box
        @test occursin("✓ converged", out_on)
        @test occursin("Equalities (1)", out_on)            # generic rows for an unnamed problem

        # One frame per IPOPT iteration (0 through N), plus the final frame.
        n = parse(Int, match(r"Converged in (\d+) iterations", out_on).captures[1])
        @test count("┌", out_on) == n + 2
    end

    @testset "a display error turns the box off; the solve is unaffected" begin
        local x, info
        out = _sm_capture(live = true) do
            @test_logs (:warn, r"status display hit an error") match_mode = :any begin
                x, _, info, _ = _sm_solve(_SM._ipopt(_sm_opts(); monitor = true);
                                          rows = _SMThrowingRows(ErrorException("boom")))
            end
        end
        @test info == :Solve_Succeeded
        @test isapprox(x, [0.5, 1.5]; atol = 1e-6)
    end

    @testset "an interrupt stops the solve cleanly" begin
        local info
        _sm_capture(live = true) do
            _, _, info, _ = _sm_solve(_SM._ipopt(_sm_opts(); monitor = true);
                                      rows = _SMThrowingRows(InterruptException()))
        end
        @test info == :User_Requested_Stop
    end
end

# ── End to end: a targeting sequence ─────────────────────────────────────────────────────────────

@testset "Targeting sequence with the status box" begin
    function _sm_hohmann()
        sat = Spacecraft(state = KeplerianState(7000.0, 0.001, 0.0, 0.0, 7.5, 1.0),
                         time = Time("2020-09-21T12:23:12", TAI(), ISOT()), name = "Sat")
        prop = OrbitPropagator(ForceModel(PointMassGravity(earth, ())),
                               IntegratorConfig(Vern9(); abstol = 1e-11, reltol = 1e-11, dt = 300.0))
        toi = ImpulsiveManeuver(axes = VNB(), element1 = 0.1, element2 = 0.0, element3 = 0.0)
        moi = ImpulsiveManeuver(axes = VNB(), element1 = 0.4, element2 = 0.0, element3 = 0.0)
        seq = Sequence()
        add_sequence!(seq,
            Event(name = "TOI", event = () -> maneuver!(sat, toi),
                  vars = [Vary(delta_v, toi; lower_bound = [0.0, 0.0, 0.0],
                               upper_bound = [2.5, 0.0, 0.0], name = "toi")]),
            Event(name = "Coast", event = () -> propagate!(prop, sat,
                      StopAt(position_dot_velocity, sat; equals = 0.0, direction = -1))),
            Event(name = "MOI", event = () -> maneuver!(sat, moi),
                  vars = [Vary(delta_v, moi; lower_bound = [0.0, 0.0, 0.0],
                               upper_bound = [3.0, 0.0, 0.0], name = "moi")],
                  funcs = [Constraint(position_magnitude, sat; equals = 45000.0),
                           Constraint(eccentricity, sat; equals = 0.0)]))
        return seq, toi, moi
    end

    seq_off, toi_off, moi_off = _sm_hohmann()
    local r_off
    out_off = _sm_capture(live = true) do
        r_off = solve!(seq_off; method = Optimize(derivatives = :fd, monitor = false, print_level = 0))
    end

    seq_on, toi_on, moi_on = _sm_hohmann()
    local r_on
    out_on = _sm_capture(live = true) do
        r_on = solve!(seq_on; method = Optimize(derivatives = :fd, monitor = true))
    end

    @test r_off.info == r_on.info == :Solve_Succeeded
    @test isapprox(r_on.variables, r_off.variables; rtol = 1e-8)
    @test !occursin("┌", out_off)

    # Named rows: the burns by their Vary names and axes, the goals by their quantity labels, with
    # units. The fixed N and B components are not shown.
    final = last(split(_sm_strip(out_on), "┌"))
    @test occursin("toi[V]", final) && occursin("moi[V]", final)
    @test !occursin("toi[N]", final)
    @test occursin("Position magnitude", final) && occursin("Eccentricity", final)
    @test occursin("km/s", final)
    @test occursin("Goals met 2/2", final)
    @test occursin("Targeter", final)

    # The session default turns it off for a solve that does not say.
    seq3, _, _ = _sm_hohmann()
    out3 = _sm_capture(live = true) do
        solver_monitor!(false)
        try
            solve!(seq3; method = Optimize(derivatives = :fd, print_level = 0))
        finally
            solver_monitor!(true)
        end
    end
    @test !occursin("┌", out3)
end

# ── Grouped goals: a collocation problem's constraints, one row per group ────────────────────────

function _sm_group_record(iter; ngroups = 4, status = :running, unmet = Int[], err = 10.0^(-iter),
                          longname = false)
    groups = map(1:ngroups) do i
        bad = i in unmet
        name = longname && i == 1 ? "a_phase_name_far_too_long_for_any_column · defects" :
               "phase_$i · defects"
        _SM.MonitorGroup(name, bad ? 590 : 600, 600,
                         bad ? "vz @ step 1234 (t 12345.6789)" : "", bad ? 1e-2 * i : err * i)
    end
    notes = ["phase_1 · HermiteSimpson, 101 nodes", "911 decision variables"]
    vars  = [_SM.MonitorVar("final_time", "", 7.1)]
    worst = maximum(g.error for g in groups)
    return _SM.IterationRecord("IPOPT", "Optimizer", iter, vars, notes, groups, 3e-4,
                               worst, 1e-3, 1e-5, worst, 1e-4, 0.1iter, status)
end

_sm_final_frame(out) = last(split(_sm_strip(out), "┌"))
_sm_goal_names(frame) = [m.captures[1] for m in eachmatch(r"│   (phase_\d+ · defects)", frame)]

@testset "Grouped goals" begin

    @testset "fit the terminal at every width, with long names and locations" begin
        for (cols, width) in ((200, 80), (80, 80), (70, 70), (64, 64))
            out = _sm_capture(cols = cols, live = false) do
                box = _SM.StatusBox()
                _SM.show_iteration!(box, _sm_group_record(0; unmet = [1, 2], longname = true))
                _SM.finish!(box, :Maximum_Iterations_Exceeded)
            end
            lines = _sm_box_lines(out)
            @test !isempty(lines)
            @test all(textwidth.(lines) .== width)
        end
        out = _sm_capture(cols = 50, live = false) do
            box = _SM.StatusBox()
            _SM.show_iteration!(box, _sm_group_record(2; unmet = [3]))
            _SM.finish!(box, :Solve_Succeeded)
        end
        first_line = first(split(out, '\n'))
        @test occursin("goals 3/4", first_line) && textwidth(first_line) <= 49
    end

    @testset "met count, location only where unmet, groups counted in Goals met" begin
        out = _sm_capture(live = false) do
            box = _SM.StatusBox()
            _SM.show_iteration!(box, _sm_group_record(3; unmet = [2]))
            _SM.finish!(box, :Maximum_Iterations_Exceeded)
        end
        frame = _sm_final_frame(out)
        @test occursin(r"Goals\s+met\s+worst at\s+error", frame)
        @test occursin(r"phase_2 · defects\s+590/600\s+vz @ step 1234", frame)
        @test occursin(r"phase_1 · defects\s+600/600\s+\d", frame)    # no location when met
        @test occursin("…", frame)                                     # a long location is cut
        @test occursin("Goals met 3/4", frame)
        @test occursin("phase_1 · HermiteSimpson, 101 nodes", frame)
        @test occursin("911 decision variables", frame)
    end

    @testset "up to 8 in the order given; more show the worst and count the rest" begin
        out = _sm_capture(live = false) do
            box = _SM.StatusBox()
            _SM.show_iteration!(box, _sm_group_record(3; ngroups = 8, unmet = [7]))
            _SM.finish!(box, :Maximum_Iterations_Exceeded)
        end
        @test _sm_goal_names(_sm_final_frame(out)) == ["phase_$i · defects" for i in 1:8]
        @test !occursin("more ·", out)

        # 12 groups: the unmet ones first, largest error first, then the met ones by error.
        out = _sm_capture(live = false) do
            box = _SM.StatusBox()
            _SM.show_iteration!(box, _sm_group_record(3; ngroups = 12, unmet = [2, 9]))
            _SM.finish!(box, :Maximum_Iterations_Exceeded)
        end
        frame = _sm_final_frame(out)
        @test _sm_goal_names(frame) == ["phase_$i · defects" for i in (9, 2, 12, 11, 10, 8, 7, 6)]
        @test occursin("+ 4 more · all met", frame)
        out = _sm_capture(live = false) do
            box = _SM.StatusBox()
            _SM.show_iteration!(box, _sm_group_record(3; ngroups = 12, unmet = 1:10))
            _SM.finish!(box, :Maximum_Iterations_Exceeded)
        end
        @test occursin("+ 4 more · 2 unmet", out)

        # The height stays the same as groups go from unmet to met.
        heights = map(0:3) do k
            out = _sm_capture(live = false) do
                box = _SM.StatusBox()
                _SM.show_iteration!(box, _sm_group_record(k; ngroups = 12, unmet = 1:(3 - k)))
                _SM.finish!(box, :Solve_Succeeded)
            end
            length(_sm_box_lines(out))
        end
        @test allequal(heights)
    end

    @testset "a live solve redraws grouped frames in place" begin
        out = _sm_capture(live = true) do
            println("earlier output")
            box = _SM.StatusBox()
            for k in 0:3
                _SM.show_iteration!(box, _sm_group_record(k; ngroups = 12, unmet = 1:(3 - k)))
            end
            _SM.finish!(box, :Solve_Succeeded)
        end
        @test startswith(out, "earlier output\n")
        ups = _sm_ups(out)
        @test length(ups) == 4 && allequal(ups)
    end
end

# ── End to end: a two-phase collocation problem ──────────────────────────────────────────────────

# A double integrator in two linked arcs, started from a guess that violates every kind of
# constraint: the defects (position moves with zero velocity), the path bound on acceleration and
# the continuity link (the arcs' guesses don't meet).
struct _SMDIState{T}   <: AbstractState;   x::T; v::T end
struct _SMDIControl{T} <: AbstractControl; a::T end
_sm_di!(dy, y::_SMDIState, u::_SMDIControl, p, t, model) = (dy[1] = y.v; dy[2] = u.a)

function _sm_di_arc(name, t0, tf, xs, a)
    ph = CollocationPhase(name = name, transcription = HermiteSimpson(n_steps = 10),
                          dynamics = _sm_di!, state = _SMDIState, control = _SMDIControl,
                          tspan = (t0, tf))
    Vary(state, ph; guess = [xs[1] xs[2]; 0.0 0.0], lower_bound = [-10.0, -10.0],
         upper_bound = [10.0, 10.0])
    Vary(control, ph; guess = fill(a, 1, 2), lower_bound = [-10.0], upper_bound = [10.0])
    return ph
end

_sm_di_rest(c)   = [state(c).x, state(c).v]
_sm_di_accel(c)  = [control(c).a]
_sm_di_effort(c) = control(c).a^2

function _sm_di_problem(; linked = true)
    first_arc = _sm_di_arc(:first_arc, 0.0, 1.0, (0.0, 0.4), 0.0)
    seq = Sequence()
    add_sequence!(seq, first_arc)
    Constraint(_sm_di_rest, first_arc; equals = [0.0, 0.0], at = Initial())
    if !linked
        Vary(final_time, first_arc; guess = 1.0, lower_bound = 0.5, upper_bound = 3.0)
        Constraint(c -> [state(c).x], first_arc; equals = [1.0], at = Final(), name = :arrival)
        Objective(_sm_di_effort, first_arc; sense = Min(), at = Path())
        return seq
    end
    second_arc = _sm_di_arc(:second_arc, 1.0, 2.0, (0.6, 1.0), 5.0)
    add_sequence!(seq, second_arc)
    Vary(final_time, second_arc; guess = 2.0, lower_bound = 1.5, upper_bound = 3.0)
    Constraint(c -> [state(c).x, state(c).v], second_arc; equals = [1.0, 0.0], at = Final(),
               name = :arrival)
    Constraint(_sm_di_accel, second_arc; lower_bound = [-2.0], upper_bound = [2.0], at = Path())
    Objective(_sm_di_effort, second_arc; sense = Min(), at = Path())
    Constraint(continuity, Link(first_arc, second_arc))
    return seq
end

@testset "Collocation sequence with the status box" begin

    @testset "groups name the phase, the kind and where the worst row is" begin
        local r
        out = _sm_capture(live = true) do
            r = solve!(_sm_di_problem(); method = Optimize(monitor = true, max_iter = 1))
        end
        @test r.info == :Maximum_Iterations_Exceeded
        frame = _sm_final_frame(out)
        @test occursin("second_arc.final_time", frame)
        @test !occursin("initial_time", frame)                  # fixed by tspan: not shown
        @test occursin("first_arc · HermiteSimpson, 21 nodes", frame)
        @test occursin("130 decision variables", frame)
        @test occursin(r"first_arc · defects\s+\d+/40\s+[xv] @ step \d+ \(t ", frame)
        @test occursin(r"second_arc · path\s+\d+/21\s+_sm_di_accel @ node \d+", frame)
        @test occursin(r"first_arc→second_arc · li…\s+\d/3\s+continuity\[[xvt]\]", frame)
        @test occursin(r"first_arc · boundary\s+2/2", frame)
        @test occursin(r"Goals met \d/6", frame)

        # Every frame of the solve has the same height.
        frames = split(_sm_strip(out), "┌")[2:end]
        @test allequal(count(==('\n'), first(split(f, "┘"))) for f in frames)
    end

    @testset "the solution is the same with the box on or off" begin
        local on, off
        out = _sm_capture(live = true) do
            on = solve!(_sm_di_problem(); method = Optimize(monitor = true))
        end
        _sm_capture(live = true) do
            off = solve!(_sm_di_problem(); method = Optimize(monitor = false, print_level = 0))
        end
        @test on.info == off.info == :Solve_Succeeded
        @test on.variables == off.variables
        @test occursin("Goals met 6/6", _sm_final_frame(out))
    end

    @testset "one phase: no phase prefix, names from the script" begin
        out = _sm_capture(live = true) do
            solve!(_sm_di_problem(linked = false); method = Optimize(monitor = true, max_iter = 2))
        end
        frame = _sm_final_frame(out)
        @test occursin(r"│   final_time\s", frame)
        @test occursin(r"│   defects\s", frame) && occursin(r"│   boundary\s", frame)
        @test !occursin("first_arc · defects", frame)
        phase = first(AstroSolve.sequence_phases(_sm_di_problem(linked = false)))
        @test _SM._constraint_labels(phase, false) == ["_sm_di_rest (initial)", "arrival"]
    end

    @testset "rows that cannot be named fall back to counts, with one warning" begin
        opts = AstroSolve.SNOW.Options(derivatives = AstroSolve.SNOW.UserDeriv(),
                                       solver = _SM._ipopt(_sm_opts(); monitor = true))
        local same
        @test_logs (:warn, r"could not name this problem's rows") begin
            same = AstroSolve._with_collocation_rows(opts, nothing)
        end
        @test same === opts
        off = AstroSolve.SNOW.Options(derivatives = AstroSolve.SNOW.UserDeriv(),
                                      solver = _SM._ipopt(_sm_opts(); monitor = false))
        @test AstroSolve._with_collocation_rows(off, nothing) === off     # box off: not built
    end
end
