# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# The nonlinear program solver every mission-design solve goes through.
#
# A problem hands over one function that evaluates the objective and constraints, and either their
# derivatives too or not; `solve_nlp` fills the derivatives a problem does not supply by forward
# finite differences, and calls IPOPT. The per-iteration callback feeds the live status box
# (solver_monitor.jl).
#
# Every IPOPT option given is passed to IPOPT, a limited-memory Hessian is used, and IPOPT's log
# file is opened at the given print level.

"""
    AstroIPOPT

AstroSolve's IPOPT solver: IPOPT's options, and whether the live status box shows while it runs.
"""
struct AstroIPOPT
    options::Dict{String, Any}      # passed to IPOPT as they are
    monitor::Bool
    file_level::Int                 # print level of IPOPT's log file
    rows::Union{Nothing, MonitorRows}
end

"""
    JacobianPattern(rows, cols)

The positions of the constraint Jacobian's declared entries, in the order the problem gives their
values. A position left out is never asked for, so a pattern that omits a real derivative converges
to the wrong answer.
"""
struct JacobianPattern
    rows::Vector{Int}
    cols::Vector{Int}
end

"""
    NLPOptions(; derivatives = :fd, sparsity = nothing, solver = _ipopt(Dict()))

How one nonlinear program is solved.

# Arguments
- `derivatives`: `:user` when the problem's function supplies its derivatives,
  `f = fun!(g, df, dg, x)`, or `:fd` to take them by forward finite differences of
  `f = fun!(g, x)`.
- `sparsity`: a [`JacobianPattern`](@ref) for `:user` derivatives, in which case `dg` is a vector
  of values in the pattern's order; `nothing` for a dense Jacobian, in which case `dg` is the
  constraints-by-variables matrix. Finite differences are always dense.
- `solver`: the IPOPT solver, from `_ipopt`.
"""
struct NLPOptions
    derivatives::Symbol
    sparsity::Union{Nothing, JacobianPattern}
    solver::AstroIPOPT

    function NLPOptions(derivatives::Symbol, sparsity, solver::AstroIPOPT)
        derivatives in (:user, :fd) || throw(ArgumentError(
            "derivatives must be :user or :fd; got :$derivatives."))
        derivatives === :fd && sparsity !== nothing && throw(ArgumentError(
            "finite-difference derivatives are dense; a sparsity pattern needs derivatives = :user."))
        return new(derivatives, sparsity, solver)
    end
end

NLPOptions(; derivatives::Symbol = :fd, sparsity = nothing,
           solver::AstroIPOPT = _ipopt(Dict{String, Any}())) =
    NLPOptions(derivatives, sparsity, solver)

"""
    solve_nlp(fun!, x0, ng, lx, ux, lg, ug, options::NLPOptions) -> (x, f, status)

Minimize f(x) subject to `lg ≤ g(x) ≤ ug` and `lx ≤ x ≤ ux`, with IPOPT. `fun!` fills the `ng`
constraints `g` and returns f; how it is called depends on `options.derivatives`, see
[`NLPOptions`](@ref). An equality constraint has `lg[i] == ug[i]`.

# Returns
The solution, its objective, and IPOPT's exit status as a symbol such as `:Solve_Succeeded`.
"""
function solve_nlp(fun!, x0, ng::Integer, lx, ux, lg, ug, options::NLPOptions)
    nx = length(x0)
    rows, cols = options.sparsity === nothing ? _dense_positions(ng, nx) :
                 (options.sparsity.rows, options.sparsity.cols)
    return _ipopt_solve(options.solver, _evaluator(fun!, options, nx, ng),
                        x0, lx, ux, lg, ug, rows, cols)
end

"""Every position of an `ng` by `nx` Jacobian, column by column, the order `vec` reads a matrix."""
function _dense_positions(ng::Integer, nx::Integer)
    rows = repeat(1:ng, nx)
    cols = repeat(1:nx; inner = ng)
    return collect(Int, rows), collect(Int, cols)
end

# One evaluation at x: the objective returned, and the constraints, objective gradient and Jacobian
# values written into g, df and dg, the values in the order of the positions IPOPT was given.
function _evaluator(fun!, options::NLPOptions, nx::Integer, ng::Integer)
    if options.derivatives === :user
        options.sparsity === nothing || return (g, df, dg, x) -> fun!(g, df, dg, x)
        J = zeros(ng, nx)
        return function (g, df, dg, x)
            f = fun!(g, df, J, x)
            copyto!(dg, J)                  # column-major, the dense positions' order
            return f
        end
    end
    # Forward differences, each step √eps relative to the variable and never below √eps.
    xp = zeros(nx)
    gp = zeros(ng)
    return function (g, df, dg, x)
        f = fun!(g, x)
        copyto!(xp, x)
        for j in 1:nx
            h = sqrt(eps(Float64)) * max(1.0, abs(x[j]))
            xp[j] = x[j] + h
            fp = fun!(gp, xp)
            df[j] = (fp - f) / h
            base = (j - 1) * ng
            for i in 1:ng
                dg[base + i] = (gp[i] - g[i]) / h
            end
            xp[j] = x[j]
        end
        return f
    end
end

# IPOPT's exit codes, as the symbols a solve reports.
const _IPOPT_STATUS = Dict(
    0 => :Solve_Succeeded, 1 => :Solved_To_Acceptable_Level, 2 => :Infeasible_Problem_Detected,
    3 => :Search_Direction_Becomes_Too_Small, 4 => :Diverging_Iterates, 5 => :User_Requested_Stop,
    6 => :Feasible_Point_Found, -1 => :Maximum_Iterations_Exceeded, -2 => :Restoration_Failed,
    -3 => :Error_In_Step_Computation, -4 => :Maximum_CpuTime_Exceeded,
    -5 => :Maximum_WallTime_Exceeded, -10 => :Not_Enough_Degrees_Of_Freedom,
    -11 => :Invalid_Problem_Definition, -12 => :Invalid_Option, -13 => :Invalid_Number_Detected,
    -100 => :Unrecoverable_Exception, -101 => :NonIpopt_Exception_Thrown,
    -102 => :Insufficient_Memory, -199 => :Internal_Error)

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
    _with_rows(options::NLPOptions, rows) -> NLPOptions

The same options with the status box's rows named by the caller, which knows the problem.
"""
function _with_rows(options::NLPOptions, rows::MonitorRows)
    s = options.solver
    return NLPOptions(options.derivatives, options.sparsity,
                      AstroIPOPT(s.options, s.monitor, s.file_level, rows))
end

"""The constraint-violation tolerance IPOPT will use, which is what a goal's ✓ is judged by."""
_constr_viol_tol(options::NLPOptions) =
    Float64(get(options.solver.options, "constr_viol_tol", 1e-4))

_ipopt_option!(prob, k, v::AbstractString) = Ipopt.AddIpoptStrOption(prob, k, String(v))
_ipopt_option!(prob, k, v::Integer)        = Ipopt.AddIpoptIntOption(prob, k, Int(v))
_ipopt_option!(prob, k, v::Real)           = Ipopt.AddIpoptNumOption(prob, k, Float64(v))

function _ipopt_solve(solver::AstroIPOPT, evaluate!, x0, lx, ux, lg, ug, rows, cols)
    nx, ng = length(x0), length(lg)
    lx, ux = Vector{Float64}(lx), Vector{Float64}(ux)
    lg, ug = Vector{Float64}(lg), Vector{Float64}(ug)

    # IPOPT asks for the objective, constraints and derivatives separately; one evaluation at a
    # point fills all of them, and is reused until the point changes.
    #
    # The point is kept as a copy, not as IPOPT's own buffer: when IPOPT reuses that memory for a
    # new point, the new point would compare equal to the "last" one and a stale objective or
    # gradient would be returned. Keeping the buffer stopped Ex_ForceModelCustomForceModel at
    # objective -1.5633 after 23 iterations; with the copy it reaches -1.5796 in 134.
    xlast = fill(NaN, nx)
    have  = Ref(false)
    f     = Ref(0.0)
    gbuf  = zeros(ng)
    df    = zeros(nx)
    dg    = zeros(length(rows))

    function at!(x)
        if !have[] || !isequal(x, xlast)
            f[] = evaluate!(gbuf, df, dg, x)
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
    info = _IPOPT_STATUS[Int(status)]
    if box !== nothing
        try
            finish!(box, info)
        catch err
            @warn "The solver status display could not draw its final frame." exception = err
        end
    end
    return prob.x, prob.obj_val, info
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
