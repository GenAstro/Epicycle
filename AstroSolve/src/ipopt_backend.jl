# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# The IPOPT driver every mission-design solve goes through.
#
# SNOW.minimize builds the derivative cache (finite differences, ForwardDiff or the user's
# partials, dense or sparse) and then calls `SNOW.optimize(solver, cache, ...)`, dispatching on the
# solver type. `AstroIPOPT` is AstroSolve's own solver type, so SNOW keeps doing the derivatives and
# this file owns the IPOPT call. What it adds over SNOW's IPOPT driver is the per-iteration
# callback, which feeds the live status box (solver_monitor.jl).
#
# Options follow SNOW's IPOPT driver, so a solve with the monitor off runs as it did through
# SNOW.IPOPT: every option is passed to IPOPT, a limited-memory Hessian is used, and the log file
# is opened at the given print level.

const Ipopt = SNOW.Ipopt

"""
    AstroIPOPT

AstroSolve's IPOPT solver for `SNOW.minimize`. It solves exactly as `SNOW.IPOPT` does and, with the
monitor on, shows the live status box while it runs.
"""
struct AstroIPOPT <: SNOW.AbstractSolver
    options::Dict{String, Any}      # passed to IPOPT as they are
    monitor::Bool
    file_level::Int                 # print level of IPOPT's log file
    rows::Union{Nothing, MonitorRows}
end

"""
    _ipopt(options; monitor = <session default>) -> AstroIPOPT

The IPOPT solver for a set of IPOPT options. With the monitor on, IPOPT's console output is turned
off so it cannot draw over the box; its log file keeps the print level the options asked for.
"""
function _ipopt(options::AbstractDict; monitor::Bool = _MONITOR_DEFAULT[])
    o = Dict{String, Any}(String(k) => v for (k, v) in options)
    file_level = Int(get(o, "print_level", 5))
    if monitor
        o["print_level"] = 0
        o["sb"] = "yes"                     # IPOPT's start-up banner
    end
    return AstroIPOPT(o, monitor, file_level, nothing)
end

"""
    _with_rows(options::SNOW.Options, rows) -> SNOW.Options

The same options with the status box's rows named by the caller, which knows the problem. Options
whose solver is not `AstroIPOPT` are returned unchanged.
"""
function _with_rows(options::SNOW.Options, rows::MonitorRows)
    s = options.solver
    s isa AstroIPOPT || return options
    return SNOW.Options(options.sparsity, options.derivatives,
                        AstroIPOPT(s.options, s.monitor, s.file_level, rows))
end

"""The constraint-violation tolerance IPOPT will use, which is what a goal's ✓ is judged by."""
_constr_viol_tol(options::SNOW.Options) =
    options.solver isa AstroIPOPT ?
        Float64(get(options.solver.options, "constr_viol_tol", 1e-4)) : 1e-4

_ipopt_option!(prob, k, v::AbstractString) = Ipopt.AddIpoptStrOption(prob, k, String(v))
_ipopt_option!(prob, k, v::Integer)        = Ipopt.AddIpoptIntOption(prob, k, Int(v))
_ipopt_option!(prob, k, v::Real)           = Ipopt.AddIpoptNumOption(prob, k, Float64(v))

function SNOW.optimize(solver::AstroIPOPT, cache, x0, lx, ux, lg, ug, rows, cols)
    nx, ng = length(x0), length(lg)
    lx, ux = Vector{Float64}(lx), Vector{Float64}(ux)
    lg, ug = Vector{Float64}(lg), Vector{Float64}(ug)

    # IPOPT asks for the objective, constraints and derivatives separately; one evaluation at a
    # point fills all of them, and is reused until the point changes.
    #
    # The point is kept as a copy. SNOW's driver keeps `xlast = x`, which is IPOPT's own buffer:
    # when IPOPT reuses that memory for a new point, the new point compares equal to the "last" one
    # and a stale objective or gradient is returned. On Ex_ForceModelCustomForceModel that stopped
    # IPOPT at objective -1.5633 after 23 iterations; with the copy it reaches -1.5796 in 134
    # (checked against SNOW's driver with only that line changed, which matches this one exactly).
    xlast = fill(NaN, nx)
    have  = Ref(false)
    f     = Ref(0.0)
    gbuf  = zeros(ng)
    df    = zeros(nx)
    dg    = zeros(length(rows))

    function at!(x)
        if !have[] || !isequal(x, xlast)
            f[] = SNOW.evaluate!(gbuf, df, dg, x, cache)
            copyto!(xlast, x)
            have[] = true
        end
        return nothing
    end
    ipobj(x)          = (at!(x); f[])
    ipcon(x, g)       = (at!(x); g .= gbuf; nothing)
    ipgrad(x, grad)   = (at!(x); grad .= df; nothing)
    function ipjac(x, r, c, values)
        if values === nothing
            r .= rows
            c .= cols
        else
            at!(x)
            values .= dg
        end
        return nothing
    end
    eval_h(args...) = false                 # a limited-memory Hessian is used instead

    prob = Ipopt.CreateIpoptProblem(nx, lx, ux, ng, lg, ug, length(rows), 1,
                                    ipobj, ipcon, ipgrad, ipjac, eval_h)
    _ipopt_option!(prob, "hessian_approximation", "limited-memory")
    for (k, v) in solver.options
        _ipopt_option!(prob, k, v)
    end
    logfile = String(get(solver.options, "output_file", "ipopt.out"))
    Ipopt.OpenIpoptOutputFile(prob, logfile, solver.file_level)

    # The log's path is named when a solve fails, unless the file was opened at level 0 and so
    # holds nothing.
    box = solver.monitor ?
          StatusBox(log = solver.file_level > 0 ? abspath(logfile) : "") : nothing
    if box !== nothing
        ctol   = Float64(get(solver.options, "constr_viol_tol", 1e-4))
        tol    = Float64(get(solver.options, "tol", 1e-8))
        mrows  = something(solver.rows, GenericRows("Optimizer", nx, lg, ug, ctol))
        title  = monitor_title(mrows)
        feas   = mrows isa EventGraphRows      # targeting: no objective to show
        xcur   = zeros(nx)
        gcur   = zeros(ng)
        t0     = time()
        on     = Ref(true)

        function intermediate(alg_mod, iter, obj, inf_pr, inf_du, mu, d_norm, reg,
                              alpha_du, alpha_pr, ls_trials)
            on[] || return true
            try
                _current_iterate!(xcur, gcur, prob, xlast, gbuf)
                vars, notes, goals, worst = monitor_rows(mrows, xcur, gcur)
                isempty(goals) && (worst = inf_du)
                rec = IterationRecord("IPOPT", title, Int(iter), vars, notes, goals,
                                      feas ? NaN : obj, inf_pr, inf_du, mu, worst,
                                      isempty(goals) ? tol : ctol, time() - t0, :running)
                show_iteration!(box, rec)
            catch err
                err isa InterruptException && rethrow()
                on[] = false
                @warn "The solver status display hit an error and was turned off; the solve " *
                      "continues." exception = (err, catch_backtrace())
            end
            return true
        end
        Ipopt.SetIntermediateCallback(prob, intermediate)
    end

    prob.x .= x0
    status = box === nothing ? Ipopt.IpoptSolve(prob) :
             with_box_logging(() -> Ipopt.IpoptSolve(prob), box)
    info = SNOW.ApplicationReturnStatus[Int(status)]
    if box !== nothing
        try
            finish!(box, info)
        catch err
            @warn "The solver status display could not draw its final frame." exception = err
        end
    end
    return prob.x, prob.obj_val, info, nothing
end

# The accepted iterate. IPOPT reports it directly (3.14 and later); otherwise the most recent
# evaluation is used, which is the accepted point or a trial step beside it.
function _current_iterate!(x, g, prob, xlast, gbuf)
    try
        Ipopt.GetIpoptCurrentIterate(prob, false, length(x), x, C_NULL, C_NULL,
                                     length(g), g, C_NULL)
    catch
        copyto!(x, xlast)
        copyto!(g, gbuf)
    end
    return nothing
end
