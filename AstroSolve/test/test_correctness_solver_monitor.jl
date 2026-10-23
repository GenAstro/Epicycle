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

# Run f with stdout captured, a fixed terminal size, and live redraw on or off. Every iteration is
# drawn unless `interval` says otherwise.
function _sm_capture(f; cols = 100, rows = 200, live = true, interval = 0.0)
    _SM._MONITOR_COLS[]     = cols
    _SM._MONITOR_ROWS[]     = rows
    _SM._MONITOR_LIVE[]     = live
    _SM._MONITOR_INTERVAL[] = interval
    path, io = mktemp()
    try
        try
            redirect_stdout(f, io)
        finally
            close(io)                    # before reading, and before rm on Windows
        end
        return read(path, String)
    finally
        _SM._MONITOR_COLS[]     = nothing
        _SM._MONITOR_ROWS[]     = nothing
        _SM._MONITOR_LIVE[]     = nothing
        _SM._MONITOR_INTERVAL[] = 0.1
        rm(path; force = true)
    end
end

_sm_strip(s) = replace(s, r"\e\[\d+A\r?" => "", "\e[J" => "", "\e[K" => "")
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
        over   = [occursin(r"\e\[\d+A\r$", b) for b in before]
        widths = [textwidth("┌" * _sm_strip(first(split(f, '\n')))) for f in frames]
        @test over == [false, true, false, true, true]
        @test widths == [80, 80, 66, 66, 80]
    end

    @testset "a redraw overwrites in place and never blanks the box first" begin
        out = _sm_capture(live = true) do
            box = _SM.StatusBox()
            for k in 0:3
                _SM.show_iteration!(box, _sm_record(k))
            end
        end
        @test length(_sm_ups(out)) == 3
        @test !occursin(r"\e\[\d+A\e\[J", out)        # no erase between moving up and drawing
        @test occursin(r"\e\[\d+A\r┌", out)            # the new frame starts on the old one's line

        # A shorter frame clears only the old lines left below it.
        out = _sm_capture(live = true) do
            box = _SM.StatusBox()
            _SM.show_iteration!(box, _sm_record(0; nvars = 6))
            _SM.show_iteration!(box, _sm_record(1; nvars = 2))
        end
        @test endswith(out, "\e[J")
    end

    @testset "a live frame fits the terminal's height; the final frame is whole" begin
        full = length(_sm_box_lines(_sm_capture(live = false) do
            box = _SM.StatusBox(); _SM.show_iteration!(box, _sm_record(0; nvars = 8, ngoals = 8))
            _SM.finish!(box, :Solve_Succeeded)
        end))
        rows = 16
        @test full > rows - 1                           # the case this guards
        out = _sm_capture(live = true, rows = rows) do
            box = _SM.StatusBox()
            for k in 0:3
                _SM.show_iteration!(box, _sm_record(k; nvars = 8, ngoals = 8))
            end
            _SM.finish!(box, :Solve_Succeeded)
        end
        frames = split(_sm_strip(out), "┌")[2:end]
        heights = [count(==('\n'), first(split(f, "┘"))) + 1 for f in frames]
        @test all(<=(rows - 1), heights[1:end-1])
        @test occursin("more", frames[1])
        @test heights[end] == full
        @test all(<=(rows - 1), _sm_ups(out))           # never asks to move above the window

        # Too short for any box: one status line.
        out = _sm_capture(live = true, rows = 6) do
            box = _SM.StatusBox()
            _SM.show_iteration!(box, _sm_record(0))
        end
        @test !occursin("┌", out) && occursin("iter 0", out)
    end

    @testset "live redraws are spaced, and the final frame is always drawn" begin
        out = _sm_capture(live = true, interval = 60.0) do
            box = _SM.StatusBox()
            for k in 0:9
                _SM.show_iteration!(box, _sm_record(k))
            end
            _SM.finish!(box, :Solve_Succeeded)
        end
        @test count("┌", out) == 2                      # the first frame, then the final one
        @test occursin("iter   9", out)
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

function _sm_di_problem(; linked = true, link_name = :link)
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
    Constraint(continuity, Link(first_arc, second_arc; name = link_name))
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
        @test occursin("second_arc · final_time", frame)
        @test !occursin("initial_time", frame)                  # fixed by tspan: not shown
        @test occursin("first_arc · HermiteSimpson, 21 nodes", frame)
        @test occursin("130 decision variables", frame)
        @test occursin(r"first_arc · defects\s+\d+/40\s+[xv] @ step \d+ \(t ", frame)
        @test occursin(r"second_arc · path\s+\d+/21\s+_sm_di_accel @ node \d+", frame)
        # Too long for the column: the phase names are cut, not the kind after them.
        @test occursin(r"first_arc→second_a… · link\s+\d/3\s+continuity\[[xvt]\]", frame)
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

# ── Names that fit: the phase is cut, the kind after it is kept ──────────────────────────────────

@testset "Row names and counts" begin
    @test _SM._fit_name("defects", 12) == "defects     "
    @test _SM._fit_name("earth_to_flyby · match point", 26) == "earth_to_fl… · match point"
    @test textwidth(_SM._fit_name("earth_to_flyby · match point", 26)) == 26
    # Variables use the same form, so their quantity survives too.
    @test endswith(rstrip(_SM._fit_name("earth_to_flyby · departure_vinf[1]", 26)),
                   "departure_vinf[1]")
    # No separator, or no room for one: the end is cut as before.
    @test _SM._fit_name("a name far too long for any column", 12) == "a name far …"
    @test _SM._fit_name("phase · a kind far too long to keep", 12) == "phase · a k…"
    @test _SM._count(1, "node") == "1 node" && _SM._count(21, "node") == "21 nodes"
end

# ── Named links ──────────────────────────────────────────────────────────────────────────────────

@testset "A named Link names its group" begin
    out = _sm_capture(live = true) do
        solve!(_sm_di_problem(link_name = :handover); method = Optimize(monitor = true, max_iter = 1))
    end
    frame = _sm_final_frame(out)
    @test occursin(r"handover · link\s+\d/3\s+continuity\[[xvt]\]", frame)
    @test !occursin("first_arc→second_a", frame)
end

# ── Shooting: a Sims-Flanagan phase ──────────────────────────────────────────────────────────────

# An Earth-Mars rendezvous on circular, coplanar orbits, as in the Sims-Flanagan tests: 20
# segments, so 10 forward and 10 backward, a throttle ball as a path constraint and the final
# mass varied with a scale.
const _SMSF_MU  = 1.32712440018e11
const _SMSF_AU  = 1.495978707e8
const _SMSF_RE  = 1.0 * _SMSF_AU
const _SMSF_RM  = 1.524 * _SMSF_AU
const _SMSF_VE  = sqrt(_SMSF_MU / _SMSF_RE)
const _SMSF_VM  = sqrt(_SMSF_MU / _SMSF_RM)
const _SMSF_TOF = π * sqrt(((_SMSF_RE + _SMSF_RM) / 2)^3 / _SMSF_MU)
const _SMSF_TH0 = π - (_SMSF_VM / _SMSF_RM) * _SMSF_TOF
const _SMSF_M0  = 1500.0

_smsf_body(r, v, θ0, ω) = t -> (r * [cos(θ0 + ω * t), sin(θ0 + ω * t), 0.0],
                                v * [-sin(θ0 + ω * t), cos(θ0 + ω * t), 0.0],
                                -(_SMSF_MU / r^2) * [cos(θ0 + ω * t), sin(θ0 + ω * t), 0.0])
const _smsf_earth = _smsf_body(_SMSF_RE, _SMSF_VE, 0.0, _SMSF_VE / _SMSF_RE)
const _smsf_mars  = _smsf_body(_SMSF_RM, _SMSF_VM, _SMSF_TH0, _SMSF_VM / _SMSF_RM)

_smsf_thrust_ball(c) = [dot(control(c), control(c))]

function _smsf_problem()
    phase = SimsFlanaganPhase(name = :earth_mars, transcription = SimsFlanagan(n_segments = 20),
                              model = PropulsionModel(mu = _SMSF_MU, Isp = 3000.0,
                                                      Tmax = 1.0e-3, g0 = 9.80665e-3),
                              ephemeris_left = _smsf_earth, ephemeris_right = _smsf_mars)
    uf = Vary(forward_control,  phase; lower_bound = fill(-2.0, 3), upper_bound = fill(2.0, 3))
    ub = Vary(backward_control, phase; lower_bound = fill(-2.0, 3), upper_bound = fill(2.0, 3))
    mf = Vary(final_mass, phase; lower_bound = 100.0, upper_bound = _SMSF_M0, scale = _SMSF_M0)
    AstroSolve.set_matchpoint_scale!(phase, [_SMSF_RE, _SMSF_RE, _SMSF_RE,
                                  _SMSF_VE, _SMSF_VE, _SMSF_VE, _SMSF_M0])
    Constraint(_smsf_thrust_ball, phase; lower_bound = 0.0, upper_bound = 1.0, at = Path())
    obj = AstroSolve.MayerObjective(ctx -> ctx.mf, phase; sense = :Max)
    AstroSolve.set_objective!(phase, obj)
    AstroSolve.add_objective_jacobian!(() -> [1.0], obj, mf)
    phase._t0, phase._tf, phase._m0 = 0.0, _SMSF_TOF, _SMSF_M0
    _, v0, _ = _smsf_earth(0.0)
    _, vf, _ = _smsf_mars(_SMSF_TOF)
    for k in 1:10
        phase._u_fwd[:, k] = 0.5 * v0 / norm(v0)
        phase._u_bwd[:, k] = 0.5 * vf / norm(vf)
    end
    phase._mf = 0.85 * _SMSF_M0
    seq = Sequence()
    add_sequence!(seq, phase)
    return seq, phase
end

# Constraint values that meet every bound, to violate one row at a time.
_sm_feasible(r) = [isfinite(lo) ? lo : (isfinite(hi) ? hi : 0.0) for (lo, hi) in zip(r.lg, r.ug)]

_sm_group(groups, name) = only(filter(g -> g.name == name, groups))

@testset "Sims-Flanagan rows" begin
    seq, phase = _smsf_problem()
    AstroSolve.initialize_sequence!(seq)
    sm   = AstroSolve._get_or_build_sm(seq)
    rows = _SM.ShootingRows(sm, 1e-4)
    x    = AstroSolve.get_decision_vector(sm)

    vars, notes, groups, worst = _SM.monitor_rows(rows, x, _sm_feasible(rows))
    @test [g.name for g in groups] == ["match point", "path"]
    @test (_sm_group(groups, "match point").total, _sm_group(groups, "path").total) == (7, 20)
    @test all(_SM._is_met, groups) && worst == 0.0
    @test notes == ["earth_mars · SimsFlanagan, 20 segments", "$(sm.n_vars) decision variables"]
    # The solver sees final_mass divided by its scale; the box shows kilograms.
    @test only(v.name for v in vars) == "final_mass"
    @test only(v.value for v in vars) ≈ 0.85 * _SMSF_M0

    # One row violated at a time: each group reports it, by component and segment.
    mp, pa = _sm_group(groups, "match point"), _sm_group(groups, "path")
    # The match point is an equality, so +2 violates it by 2; the throttle ball sits at its lower
    # bound of 0, so +2 passes its upper bound of 1 by 1. "_smsf_thrust_ball @ fwd seg 3" is
    # wider than the column, so the location keeps the segment and drops the name.
    for (row, group, where, viol) in ((4, "match point", "vx", 2.0), (7, "match point", "m", 2.0),
                                      (7 + 3,  "path", "fwd seg 3", 1.0),
                                      (7 + 12, "path", "bwd seg 2", 1.0))
        g = _sm_feasible(rows)
        g[row] += 2.0
        _, _, gs, w = _SM.monitor_rows(rows, x, g)
        hit = _sm_group(gs, group)
        @test hit.met == hit.total - 1
        @test hit.where == where
        @test w ≈ viol
        @test all(_SM._is_met, filter(gg -> gg.name != group, gs))
    end
end

@testset "Sims-Flanagan solve with the status box" begin
    local on, off
    out = _sm_capture(live = true) do
        on = solve!(first(_smsf_problem()); method = Optimize(monitor = true, max_iter = 3))
    end
    _sm_capture(live = true) do
        off = solve!(first(_smsf_problem()); method = Optimize(monitor = false, max_iter = 3,
                                                                print_level = 0))
    end
    @test on.info == off.info
    @test on.variables == off.variables
    frame = _sm_final_frame(out)
    @test occursin(r"│   final_mass\s+\d", frame)
    @test occursin(r"match point\s+\d/7", frame) && occursin(r"path\s+\d+/20", frame)
    frames = split(_sm_strip(out), "┌")[2:end]
    @test allequal(count(==('\n'), first(split(f, "┘"))) for f in frames)
end

# ── Shooting: two MGA legs and a named flyby link ────────────────────────────────────────────────

const _SMMG_RV   = 0.723 * _SMSF_AU
const _SMMG_TDEP = 12222.0 * 86400.0
const _SMMG_TFLY = _SMMG_TDEP + 365.0 * 86400.0
const _SMMG_TARR = _SMMG_TFLY + 175.0 * 86400.0
const _smmg_earth = _smsf_body(_SMSF_RE, _SMSF_VE, 0.0, _SMSF_VE / _SMSF_RE)
const _smmg_venus = _smsf_body(_SMMG_RV, sqrt(_SMSF_MU / _SMMG_RV), 1.1,
                               sqrt(_SMSF_MU / _SMMG_RV) / _SMMG_RV)

_smmg_epoch(c1, c2)  = [final_time(c1) - initial_time(c2)]
_smmg_mass(c1, c2)   = [final_mass(c1) - initial_mass(c2)]
_smmg_launch_c3(c)   = [dot(departure_vinf(c), departure_vinf(c))]
@partial(_smmg_launch_c3, departure_vinf) do c
    reshape(2 .* departure_vinf(c), 1, 3)
end
_smmg_speed(c1, c2)  = [dot(arrival_vinf(c1), arrival_vinf(c1)) -
                        dot(departure_vinf(c2), departure_vinf(c2))]

function _smmg_problem()
    trans = MGAnDSMs(n_dsm = 1)
    model = PropulsionModel(mu = _SMSF_MU, Isp = 320.0, g0 = 9.80665e-3)
    p1 = MGAnDSMsPhase(name = :earth_to_flyby, transcription = trans, model = model,
                       ephemeris_left = _smmg_earth, ephemeris_right = _smmg_earth)
    p2 = MGAnDSMsPhase(name = :flyby_to_venus, transcription = trans, model = model,
                       ephemeris_left = _smmg_earth, ephemeris_right = _smmg_venus)
    for (p, t0, tf) in ((p1, _SMMG_TDEP, _SMMG_TFLY), (p2, _SMMG_TFLY, _SMMG_TARR))
        p._alpha = [0.5, 0.5]
        Vary(departure_vinf, p; guess = [3.0, 1.0, 0.0],
             lower_bound = fill(-6.0, 3), upper_bound = fill(6.0, 3))
        Vary(arrival_vinf, p; guess = [-3.0, 1.0, 0.0],
             lower_bound = fill(-6.0, 3), upper_bound = fill(6.0, 3))
        Vary(initial_time, p; guess = t0,
             lower_bound = t0 - 30 * 86400.0, upper_bound = t0 + 30 * 86400.0)
        Vary(final_time, p; guess = tf,
             lower_bound = tf - 60 * 86400.0, upper_bound = tf + 60 * 86400.0)
        Vary(initial_mass, p; guess = 2000.0, lower_bound = 10.0, upper_bound = 4000.0,
             scale = 2000.0)
    end
    flyby = Link(p1, p2; name = :earth_flyby)
    Constraint(_smmg_epoch, flyby; equals = 0.0, name = :epoch)
    Constraint(_smmg_mass,  flyby; equals = 0.0, name = :mass)
    Constraint(_smmg_speed, flyby; equals = 0.0)                 # unnamed: its function names it
    Constraint(_smmg_launch_c3, p1; lower_bound = 0.0, upper_bound = 50.0)   # a leg's own bound
    seq = Sequence()
    add_sequence!(seq, p1)
    add_sequence!(seq, p2)
    return seq, p1, p2
end

@testset "MGA rows" begin
    seq, p1, p2 = _smmg_problem()
    AstroSolve.initialize_sequence!(seq)
    sm   = AstroSolve._get_or_build_sm(seq)
    rows = _SM.ShootingRows(sm, 1e-4)
    x    = AstroSolve.get_decision_vector(sm)

    vars, notes, groups, _ = _SM.monitor_rows(rows, x, _sm_feasible(rows))
    @test [g.name for g in groups] ==
          ["earth_to_flyby · match point", "earth_to_flyby · boundary",
           "flyby_to_venus · match point", "earth_flyby · link"]
    # A leg's own constraint is a boundary row, named by its function.
    g = _sm_feasible(rows)
    g[first(rows.groups[2].blocks[1].rows)] += 100.0              # C3 past its upper bound
    _, _, gs, _ = _SM.monitor_rows(rows, x, g)
    @test _sm_group(gs, "earth_to_flyby · boundary").where == "_smmg_launch_c3"
    @test _sm_group(groups, "earth_flyby · link").total == 3
    @test notes[1:2] == ["earth_to_flyby · MGAnDSMs, 1 deep-space maneuver",
                         "flyby_to_venus · MGAnDSMs, 1 deep-space maneuver"]
    names = [v.name for v in vars]
    @test "earth_to_flyby · initial_time" in names && "flyby_to_venus · final_time" in names
    @test "earth_to_flyby · departure_vinf[2]" in names
    # Scaled for the solver, shown as written: seconds and kilograms.
    @test vars[findfirst(==("earth_to_flyby · initial_time"), names)].value ≈ _SMMG_TDEP
    @test vars[findfirst(==("flyby_to_venus · initial_mass"), names)].value ≈ 2000.0

    # The link's rows by their names; the unnamed one by its function.
    link = rows.groups[end].blocks
    @test [b.label for b in link] == ["epoch", "mass", "_smmg_speed"]
    for (k, label) in enumerate(("epoch", "mass", "_smmg_speed"))
        g = _sm_feasible(rows)
        g[first(link[k].rows)] += 1.0
        _, _, gs, _ = _SM.monitor_rows(rows, x, g)
        @test _sm_group(gs, "earth_flyby · link").where == label
    end
    g = _sm_feasible(rows)
    g[first(rows.groups[3].blocks[1].rows) + 2] += 1.0                # leg 2, match point z
    _, _, gs, _ = _SM.monitor_rows(rows, x, g)
    @test _sm_group(gs, "flyby_to_venus · match point").where == "z"
end

@testset "A shooting problem whose rows cannot be named still solves" begin
    opts = AstroSolve.SNOW.Options(derivatives = AstroSolve.SNOW.UserDeriv(),
                                   solver = _SM._ipopt(_sm_opts(); monitor = true))
    local same
    @test_logs (:warn, r"could not name this problem's rows") begin
        same = AstroSolve._with_grouped_rows(opts, _SM.ShootingRows, nothing)
    end
    @test same === opts
end

# ── Mixed transcriptions: a Sims-Flanagan half joined to a collocation half ──────────────────────

# The orbit raise of Ex_MixedTranscription in canonical units, on a coarse mesh: the first half
# Sims-Flanagan, the second Hermite-Simpson, joined by continuity in state, mass and time.
const _SMMX_MU, _SMMX_TMAX, _SMMX_MDOT, _SMMX_M0 = 1.0, 0.1, 0.05, 1.0
const _SMMX_TMID, _SMMX_TF = 1.5, 3.0
const _SMMX_Y0 = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
const _SMMX_LO = [-5.0, -5.0, -5.0, -3.0, -3.0, -3.0]
const _SMMX_HI = [5.0, 5.0, 5.0, 3.0, 3.0, 3.0]

struct _SMMXState{T} <: AbstractState
    x::T; y::T; z::T; vx::T; vy::T; vz::T; m::T
end
struct _SMMXControl{T} <: AbstractControl
    ux::T; uy::T; uz::T
end
struct _SMMXModel
    mu::Float64; thrust::Float64; mdot::Float64
end

function _smmx_raise!(dy, y::_SMMXState, u::_SMMXControl, p, t, m::_SMMXModel)
    r = sqrt(y.x^2 + y.y^2 + y.z^2)
    a = m.thrust / y.m
    dy[1], dy[2], dy[3] = y.vx, y.vy, y.vz
    dy[4] = -m.mu * y.x / r^3 + a * u.ux
    dy[5] = -m.mu * y.y / r^3 + a * u.uy
    dy[6] = -m.mu * y.z / r^3 + a * u.uz
    dy[7] = -m.mdot * sqrt(u.ux^2 + u.uy^2 + u.uz^2 + 1e-12)
end

# Thrust along the velocity, integrated with RK4: a guess both halves agree with.
function _smmx_flown(t_end; steps = 400)
    f(y) = (r = norm(y[1:3]); v = norm(y[4:6]);
            vcat(y[4:6], -_SMMX_MU .* y[1:3] ./ r^3 .+ (_SMMX_TMAX / y[7]) .* y[4:6] ./ v,
                 -_SMMX_MDOT))
    y, h = vcat(_SMMX_Y0, _SMMX_M0), t_end / steps
    for _ in 1:steps
        k1 = f(y); k2 = f(y .+ h / 2 .* k1); k3 = f(y .+ h / 2 .* k2); k4 = f(y .+ h .* k3)
        y = y .+ (h / 6) .* (k1 .+ 2k2 .+ 2k3 .+ k4)
    end
    return y
end

_smmx_throttle(c)    = [dot(control(c), control(c))]
_smmx_unit_thrust(c) = [control(c).ux^2 + control(c).uy^2 + control(c).uz^2]
_smmx_radius(c)      = (y = state(c); sqrt(y.x^2 + y.y^2 + y.z^2))
_smmx_circular(c)    = (y = state(c); [y.x * y.vx + y.y * y.vy + y.z * y.vz])

function _smmx_problem(; link_name = :link)
    y_mid, y_end = _smmx_flown(_SMMX_TMID), _smmx_flown(_SMMX_TF)
    model = _SMMXModel(_SMMX_MU, _SMMX_TMAX, _SMMX_MDOT)
    first_half = SimsFlanaganPhase(name = :first_half, transcription = SimsFlanagan(n_segments = 6),
                                   model = PropulsionModel(mu = _SMMX_MU,
                                                           Isp = _SMMX_TMAX / _SMMX_MDOT,
                                                           Tmax = _SMMX_TMAX, g0 = 1.0),
                                   tspan = (0.0, _SMMX_TMID))
    Vary(state, first_half; guess = _SMMX_Y0, lower_bound = _SMMX_Y0, upper_bound = _SMMX_Y0)
    Vary(final_state, first_half; guess = y_mid[1:6], lower_bound = _SMMX_LO,
         upper_bound = _SMMX_HI)
    Vary(forward_control,  first_half; lower_bound = fill(-1.0, 3), upper_bound = fill(1.0, 3))
    Vary(backward_control, first_half; lower_bound = fill(-1.0, 3), upper_bound = fill(1.0, 3))
    Vary(initial_mass, first_half; lower_bound = _SMMX_M0, upper_bound = _SMMX_M0)
    Vary(final_mass, first_half; guess = y_mid[7], lower_bound = 0.5, upper_bound = _SMMX_M0)
    Constraint(_smmx_throttle, first_half; lower_bound = 0.0, upper_bound = 1.0, at = Path())

    second_half = CollocationPhase(name = :second_half, transcription = HermiteSimpson(n_steps = 6),
                                   dynamics = _smmx_raise!, model = model,
                                   state = _SMMXState, control = _SMMXControl,
                                   tspan = (_SMMX_TMID, _SMMX_TF))
    Vary(state, second_half; guess = hcat(y_mid, 0.5 .* (y_mid .+ y_end), y_end),
         lower_bound = vcat(_SMMX_LO, 0.5), upper_bound = vcat(_SMMX_HI, _SMMX_M0))
    Vary(control, second_half; guess = repeat([0.0, 1.0, 0.0], 1, 3),
         lower_bound = fill(-1.0, 3), upper_bound = fill(1.0, 3))
    Constraint(_smmx_unit_thrust, second_half; equals = 1.0, at = Path())

    seq = Sequence()
    add_sequence!(seq, first_half)
    add_sequence!(seq, second_half)
    Constraint(continuity, Link(first_half, second_half; name = link_name))
    Constraint(_smmx_circular, second_half; equals = 0.0, at = Final())
    Objective(_smmx_radius, second_half; sense = Max())
    return seq
end

@testset "Mixed-transcription rows" begin
    seq = _smmx_problem()
    AstroSolve.initialize_sequence!(seq)
    om   = AstroSolve._get_or_build_oc(seq)
    rows = _SM.OCRows(om, 1e-4)
    x    = AstroSolve.get_decision_vector(om)

    vars, notes, groups, _ = _SM.monitor_rows(rows, x, _sm_feasible(rows))
    # Each phase reports as its own transcription does; the link comes after the phases.
    @test [g.name for g in groups] ==
          ["first_half · match point", "first_half · path", "second_half · defects",
           "second_half · path", "second_half · boundary", "first_half→second_half · link"]
    @test notes[1:2] == ["first_half · SimsFlanagan, 6 segments",
                         "second_half · HermiteSimpson, 13 nodes"]
    @test "first_half · final_mass" in [v.name for v in vars]
    @test only(v.value for v in vars if v.name == "first_half · final_mass") ≈
          _smmx_flown(_SMMX_TMID)[7]

    # The collocation half's defects carry its own node times, though it is the second phase.
    defects = _sm_group(groups, "second_half · defects")
    block   = only(rows.groups[3].blocks)
    g = _sm_feasible(rows)
    g[first(block.rows) + 7 * 2 + 3] += 1.0                      # step 2, first equation, vx
    _, _, gs, _ = _SM.monitor_rows(rows, x, g)
    @test startswith(_sm_group(gs, "second_half · defects").where, "vx @ step 2 (t 1.")

    # Continuity runs over the Sims-Flanagan half's state, mass and time.
    link = only(rows.groups[end].blocks)
    @test link.label == "continuity"
    @test link.comps == ["x", "y", "z", "vx", "vy", "vz", "m", "t"]
    g = _sm_feasible(rows)
    g[first(link.rows) + 6] += 1.0
    _, _, gs, _ = _SM.monitor_rows(rows, x, g)
    @test _sm_group(gs, "first_half→second_half · link").where == "continuity[m]"

    # The Sims-Flanagan half's throttle, by segment: 3 forward, then 3 backward. With its name the
    # location is wider than the column, so the segment is kept and the name dropped.
    g = _sm_feasible(rows)
    g[first(rows.groups[2].blocks[1].rows) + 4] += 2.0
    _, _, gs, _ = _SM.monitor_rows(rows, x, g)
    @test _sm_group(gs, "first_half · path").where == "bwd seg 2"
end

@testset "Mixed-transcription solve with the status box" begin
    local on, off
    out = _sm_capture(live = true) do
        on = solve!(_smmx_problem(link_name = :handoff); method = Optimize(monitor = true, max_iter = 5))
    end
    _sm_capture(live = true) do
        off = solve!(_smmx_problem(); method = Optimize(monitor = false, max_iter = 5, print_level = 0))
    end
    @test on.info == off.info
    @test on.variables == off.variables
    frame = _sm_final_frame(out)
    @test occursin("handoff · link", frame)
    @test occursin("first_half · match point", frame) && occursin("second_half · defects", frame)
    @test occursin(r"Goals met \d/6", frame)
    frames = split(_sm_strip(out), "┌")[2:end]
    @test allequal(count(==('\n'), first(split(f, "┘"))) for f in frames)
end

@testset "Phases listed in the order they fly" begin
    # Stand-ins: only identity matters to the ordering.
    a, b, c, d = Ref(:a), Ref(:b), Ref(:c), Ref(:d)
    @test _SM._flight_order([a, b, c], Tuple{Any,Any}[]) == [1, 2, 3]       # no links: as given
    @test _SM._flight_order([c, a, b], [(a, b), (b, c)]) == [2, 3, 1]       # a → b → c
    @test _SM._flight_order([b, a, d], [(a, b)]) == [2, 1, 3]               # d is unlinked
    @test sort(_SM._flight_order([a, b], [(a, b), (b, a)])) == [1, 2]       # a cycle still ends
end

# ── Where IPOPT's log is ─────────────────────────────────────────────────────────────────────────

@testset "A failed solve names IPOPT's log" begin
    finished(status; log = raw"C:\tmp\ipopt_run.out", record = true) =
        _sm_capture(live = false) do
            box = _SM.StatusBox(log = log)
            record && _SM.show_iteration!(box, _sm_record(2))
            _SM.finish!(box, status)
        end

    out = finished(:Maximum_Iterations_Exceeded)
    @test occursin("Stopped after 2 iterations", out)
    @test occursin(raw"IPOPT log: C:\tmp\ipopt_run.out", out)
    @test endswith(rstrip(out), "ipopt_run.out")                      # the last line, in full

    @test !occursin("IPOPT log", finished(:Solve_Succeeded))           # success: nothing to look up
    @test !occursin("IPOPT log", finished(:Maximum_Iterations_Exceeded; log = ""))   # empty log

    # Stopped before any iteration: no frame, but the status and the log are still given.
    out = finished(:Not_Enough_Degrees_Of_Freedom; record = false)
    @test !occursin("┌", out)
    @test occursin("Stopped before the first iteration (Not_Enough_Degrees_Of_Freedom)", out)
    @test occursin("IPOPT log: ", out)

    # Through the IPOPT driver: the named file is the one IPOPT wrote. Level 0 writes nothing, so
    # it is not named.
    logfile = tempname() * "_ipopt.out"
    opts = Dict{String,Any}("print_level" => 5, "max_iter" => 1, "output_file" => logfile)
    out = _sm_capture(live = false) do
        _sm_solve(_SM._ipopt(opts; monitor = true))
    end
    @test occursin("IPOPT log: " * abspath(logfile), out)
    @test isfile(logfile) && filesize(logfile) > 0
    out0 = _sm_capture(live = false) do
        _sm_solve(_SM._ipopt(merge(opts, Dict("print_level" => 0,
                                              "output_file" => tempname() * "_ipopt.out"));
                             monitor = true))
    end
    @test occursin("Stopped after", out0) && !occursin("IPOPT log", out0)
end

# ── Goal names: a Constraint's name, else the quantity, with its event where that repeats ───────

@testset "Targeting goals named apart" begin
    geo = [("", "Position magnitude", "Prop to apoapsis"), ("", "Inclination", "Prop to perigee 2"),
           ("", "Position magnitude", "Prop to perigee 2"), ("", "Semi-major axis", "MOI")]
    @test _SM._goal_names(geo) == ["Position magnitude @ Prop to apoapsis", "Inclination",
                                   "Position magnitude @ Prop to perigee 2", "Semi-major axis"]
    # A name the script gives is used as written, and leaves the other goal unique.
    @test _SM._goal_names([("Radius at apoapsis", "Position magnitude", "Prop to apoapsis"),
                           ("", "Position magnitude", "Prop to perigee 2")]) ==
          ["Radius at apoapsis", "Position magnitude"]
    # The same quantity twice in one event still needs a number.
    @test _SM._goal_names([("", "Eccentricity", "A"), ("", "Eccentricity", "A")]) ==
          ["Eccentricity @ A #1", "Eccentricity @ A #2"]
    # An event with no name adds nothing to the name.
    @test _SM._goal_names([("", "Eccentricity", ""), ("", "Eccentricity", "")]) ==
          ["Eccentricity #1", "Eccentricity #2"]

    # Too long for the column: both the quantity and the event are cut, so both still show.
    fitted = _SM._fit_name("Position magnitude @ Prop to perigee 2", 26)
    @test fitted == "Position m… @ Prop to per…" && textwidth(fitted) == 26
    @test _SM._fit_name("Radius @ Prop to perigee 2 long event", 26) == "Radius @ Prop to perigee …"
end

# ── One phase solved on its own ──────────────────────────────────────────────────────────────────

@testset "Rows for a single phase" begin
    # A shooting phase: the same groups as inside a manager, with no phase prefix.
    _, sf = _smsf_problem()
    rows  = _SM.PhaseRows(sf, 1e-4)
    x     = AstroSolve.get_decision_vector(sf)
    vars, notes, groups, _ = _SM.monitor_rows(rows, x, _sm_feasible(rows))
    @test [g.name for g in groups] == ["match point", "path"]
    @test notes == ["earth_mars · SimsFlanagan, 20 segments",
                    "$(AstroSolve.nlp_length(sf)) decision variables"]
    @test only(v.value for v in vars) ≈ 0.85 * _SMSF_M0

    # A collocation phase: defects with their node times.
    seq = _sm_di_problem(linked = false)
    AstroSolve.initialize_sequence!(seq)
    ph   = first(AstroSolve.sequence_phases(seq))
    rows = _SM.PhaseRows(ph, 1e-4)
    x    = AstroSolve.get_decision_vector(ph)
    _, _, groups, _ = _SM.monitor_rows(rows, x, _sm_feasible(rows))
    @test [g.name for g in groups] == ["defects", "boundary"]
    g = _sm_feasible(rows)
    g[2 * 2 * 3 + 2] += 1.0                                   # step 4, first equation, v
    _, _, gs, _ = _SM.monitor_rows(rows, x, g)
    @test startswith(_sm_group(gs, "defects").where, "v @ step 4 (t ")
end

# ── A transcription the box has never seen ───────────────────────────────────────────────────────

# A stand-in for a transcription with per-interval defect blocks and a row of its own, as an
# extension package's might have: the box knows only the names its function list gives.
struct _SMIntervalPhase
    name::Symbol
end
AstroSolve.function_list(::_SMIntervalPhase) = [
    AstroSolve.PhaseFunction(nothing, 2, zeros(2), zeros(2), "defect"),       # one interval
    AstroSolve.PhaseFunction(nothing, 2, zeros(2), zeros(2), "defect_2"),     # numbered
    AstroSolve.PhaseFunction(nothing, 1, [1.0], [1.0], "alpha_sum")]
AstroSolve.n_constraints(::_SMIntervalPhase)         = 5
AstroSolve.nlp_length(::_SMIntervalPhase)            = 0
AstroSolve.variable_list(::_SMIntervalPhase)         = []
AstroSolve.variable_ranges(::_SMIntervalPhase)       = UnitRange{Int}[]
AstroSolve.get_constraint_bounds(::_SMIntervalPhase) = (vcat(zeros(4), 1.0), vcat(zeros(4), 1.0))

@testset "A transcription known only by its function list" begin
    rows = _SM.PhaseRows(_SMIntervalPhase(:mock), 1e-4)
    _, notes, groups, _ = _SM.monitor_rows(rows, Float64[], _sm_feasible(rows))
    @test [g.name for g in groups] == ["defects", "boundary"]
    @test notes[1] == "mock · _SMIntervalPhase"                 # no transcription line it knows
    @test _sm_group(groups, "defects").total == 4

    for (row, group, where) in ((2, "defects", "state[2]"),            # the unnumbered block
                                (3, "defects", "state[1] @ interval 2"),
                                (5, "boundary", "alpha sum"))
        g = _sm_feasible(rows)
        g[row] += 1.0
        _, _, gs, _ = _SM.monitor_rows(rows, Float64[], g)
        @test _sm_group(gs, group).where == where
    end
end

@testset "Many phases: the first few listed, the rest counted" begin
    lines = ["p$i · HermiteSimpson, 21 nodes" for i in 1:6]
    @test _SM._phase_notes(lines, 900) ==
          ["p1 · HermiteSimpson, 21 nodes", "p2 · HermiteSimpson, 21 nodes",
           "p3 · HermiteSimpson, 21 nodes", "+ 3 more phases", "900 decision variables"]
    @test length(_SM._phase_notes(lines[1:4], 900)) == 5          # four fit: all listed
end
