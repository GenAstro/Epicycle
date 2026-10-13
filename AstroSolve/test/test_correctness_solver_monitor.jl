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
