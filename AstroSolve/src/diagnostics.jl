# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# report_variables, report_functions, and the finite-difference partial checks.
#

# ─────────────────────────────────────────────────────────────────────────────
# report_variables / report_functions — diagnostic tables
#
#   report_variables(phase)   prints each NLP variable with lb / value / ub
#                             in scaled (NLP) space; flags active bounds
#
#   report_functions(phase)   prints each constraint row with lb / value / ub
#                             in scaled (NLP) space; flags violated rows
# ─────────────────────────────────────────────────────────────────────────────

# Values and bounds come from the decision vector, in the NLP scaling report_functions also uses.
# A variable's own `value` is what it was declared with: on a Sims-Flanagan phase that is one
# three-component throttle, while the phase holds one per segment, so reading it printed a few rows
# of stale numbers for a vector of dozens.
function report_variables(p::AbstractShootingPhase)
    println()
    println("── Variables (nlp-scaled) ────────────────────────────────────────────────")
    @printf "  %-20s  %16s  %16s  %16s\n" "name" "lower" "value" "upper"
    x_nlp  = get_decision_vector(p)
    for (v, r) in zip(variable_list(p), variable_ranges(p))
        lb_nlp, ub_nlp = nlp_bounds(p, v)
        xn = x_nlp[r]
        n  = length(r)
        for i in 1:n
            atlo  = xn[i] <= lb_nlp[i] + 1e-6*(1 + abs(lb_nlp[i]))
            athi  = xn[i] >= ub_nlp[i] - 1e-6*(1 + abs(ub_nlp[i]))
            flag  = atlo ? "  ← LB" : athi ? "  ← UB" : ""
            label = n == 1 ? v.name : "$(v.name)[$i]"
            @printf "  %-20s  %16.8g  %16.8g  %16.8g%s\n" label lb_nlp[i] xn[i] ub_nlp[i] flag
        end
    end
    println()
end

function report_functions(p::AbstractShootingPhase)
    get_functions(p)          # ensure cached trajectory state is fresh
    F        = get_functions(p)
    lb, ub   = get_constraint_bounds(p)
    println()
    println("── Functions (nlp-scaled) ────────────────────────────────────────────────")
    @printf "  %-20s  %16s  %16s  %16s\n" "name" "lower" "value" "upper"
    row = 0
    for pf in function_list(p)
        for i in 1:pf.n_nlp
            row  += 1
            viol  = max(0.0, lb[row] - F[row], F[row] - ub[row])
            flag  = viol > 1e-6 ? @sprintf("  ← viol %.2e", viol) : ""
            label = pf.n_nlp == 1 ? pf.name : "$(pf.name)[$i]"
            @printf "  %-20s  %16.8g  %16.8g  %16.8g%s\n" label lb[row] F[row] ub[row] flag
        end
    end
    println()
end

# ─────────────────────────────────────────────────────────────────────────────
# test_partials_fd
#
# Central-difference check of every analytic Jacobian block and the objective
# gradient against finite differences.  Reports max absolute and relative
# error per (constraint_group, variable) block, and flags failures.
#
#   test_partials_fd(phase)
#   test_partials_fd(phase; h = 1e-6, tol = 1e-4, verbose = true)
#
# Returns true if all blocks pass, false otherwise.
#
# NOTE: AD version blocked until phase fields are made generic (Vector{T}).
#       Restructure evaluate_matchpoint! to accept a plain x-vector for that.
# ─────────────────────────────────────────────────────────────────────────────

function test_partials_fd(p::AbstractShootingPhase;
                           h       = 1e-6,
                           tol     = 1e-4,
                           verbose = true)

    x0    = get_decision_vector(p)
    F0    = get_functions(p)
    n_x   = length(x0)
    n_f   = length(F0)
    vlist = variable_list(p)
    flist = function_list(p)
    rngs  = variable_ranges(p)

    # ── Full central-difference Jacobian ──────────────────────────────────────
    J_fd = zeros(n_f, n_x)
    for j in 1:n_x
        xp = copy(x0);  xp[j] += h
        xm = copy(x0);  xm[j] -= h
        set_decision_vector!(p, xp);  Fp = get_functions(p)
        set_decision_vector!(p, xm);  Fm = get_functions(p)
        J_fd[:, j] = (Fp .- Fm) ./ (2h)
    end
    set_decision_vector!(p, x0)

    # ── Assemble analytic Jacobian from jacobian_chunk ─────────────────────
    J_an = zeros(n_f, n_x)
    f_off = 0
    for pf in flist
        for (v, r) in zip(vlist, rngs)
            chunk = jacobian_chunk(p, pf, v)
            J_an[f_off+1:f_off+pf.n_nlp, r] .+= chunk
        end
        f_off += pf.n_nlp
    end

    # ── Report by (constraint_group × variable) block ─────────────────────
    println()
    println("── Constraint Jacobian check  (central FD, h=$h) ─────────────────")
    @printf "  %-24s  %-16s  %12s  %12s\n" "function" "variable" "max_abs" "max_rel"

    any_fail = false
    f_off = 0
    for pf in flist
        frows = f_off+1 : f_off+pf.n_nlp
        for (v, r) in zip(vlist, rngs)
            an      = J_an[frows, r]
            fd      = J_fd[frows, r]
            err     = abs.(an .- fd)
            denom   = max.(abs.(fd), 1e-8)
            max_abs = maximum(err)
            max_rel = maximum(err ./ denom)
            fail    = max_abs > tol
            any_fail |= fail
            if verbose || fail
                flag = fail ? "  ✗" : "  ✓"
                @printf "  %-24s  %-16s  %12.3e  %12.3e%s\n" pf.name v.name max_abs max_rel flag
            end
        end
        f_off += pf.n_nlp
    end

    # ── Objective gradient ────────────────────────────────────────────────────
    if !isnothing(p.objective)
        g_an = zeros(n_x)
        for (v, r) in zip(vlist, rngs)
            g_an[r] .= objective_gradient_chunk(p, v)
        end

        g_fd = zeros(n_x)
        for j in 1:n_x
            xp = copy(x0);  xp[j] += h
            xm = copy(x0);  xm[j] -= h
            set_decision_vector!(p, xp);  Jp = get_objective(p)
            set_decision_vector!(p, xm);  Jm = get_objective(p)
            g_fd[j] = (Jp - Jm) / (2h)
        end
        set_decision_vector!(p, x0)

        println()
        println("── Objective gradient check ──────────────────────────────────────")
        @printf "  %-16s  %12s  %12s\n" "variable" "max_abs" "max_rel"
        for (v, r) in zip(vlist, rngs)
            an      = g_an[r]
            fd      = g_fd[r]
            err     = abs.(an .- fd)
            denom   = max.(abs.(fd), 1e-8)
            max_abs = maximum(err)
            max_rel = maximum(err ./ denom)
            fail    = max_abs > tol
            any_fail |= fail
            if verbose || fail
                flag = fail ? "  ✗" : "  ✓"
                @printf "  %-16s  %12.3e  %12.3e%s\n" v.name max_abs max_rel flag
            end
        end
    end

    println()
    if any_fail
        println("FAIL — one or more blocks exceed tol=$tol")
    else
        println("PASS — all blocks within tol=$tol")
    end
    return !any_fail
end

# ─────────────────────────────────────────────────────────────────────────────
# dump_partials_fd
#
# Full element-wise dump of every analytic and FD partial, plus errors.
# One row per (function_row, variable_col) scalar entry.
#
#   dump_partials_fd(phase)
#   dump_partials_fd(phase; h = 1e-6)
#
# Columns: function, variable, analytic, FD, abs_err, rel_err
# ─────────────────────────────────────────────────────────────────────────────

function dump_partials_fd(p::AbstractShootingPhase; h = 1e-6)

    x0    = get_decision_vector(p)
    n_x   = length(x0)
    n_f   = length(get_functions(p))
    vlist = variable_list(p)
    flist = function_list(p)
    rngs  = variable_ranges(p)

    # ── Central-difference Jacobian ───────────────────────────────────────────
    J_fd = zeros(n_f, n_x)
    for j in 1:n_x
        xp = copy(x0);  xp[j] += h
        xm = copy(x0);  xm[j] -= h
        set_decision_vector!(p, xp);  Fp = get_functions(p)
        set_decision_vector!(p, xm);  Fm = get_functions(p)
        J_fd[:, j] = (Fp .- Fm) ./ (2h)
    end
    set_decision_vector!(p, x0)

    # ── Analytic Jacobian ─────────────────────────────────────────────────────
    J_an = zeros(n_f, n_x)
    f_off = 0
    for pf in flist
        for (v, r) in zip(vlist, rngs)
            chunk = jacobian_chunk(p, pf, v)
            J_an[f_off+1:f_off+pf.n_nlp, r] .+= chunk
        end
        f_off += pf.n_nlp
    end

    # ── Dump constraint Jacobian ──────────────────────────────────────────────
    println()
    println("── Partial derivative dump  (central FD, h=$h) ──────────────────────────────────")
    @printf "  %-22s  %-16s  %13s  %13s  %13s  %13s\n" "function" "variable" "analytic" "fd" "abs_err" "rel_err"

    f_off = 0
    for pf in flist
        for i in 1:pf.n_nlp
            fname = pf.n_nlp == 1 ? pf.name : "$(pf.name)[$i]"
            grow  = f_off + i
            for (v, r) in zip(vlist, rngs)
                for (k, gcol) in enumerate(r)
                    an      = J_an[grow, gcol]
                    fd      = J_fd[grow, gcol]
                    abs_err = abs(an - fd)
                    rel_err = abs_err / max(abs(fd), 1e-8)
                    vname   = length(r) == 1 ? v.name : "$(v.name)[$k]"
                    @printf "  %-22s  %-16s  %13.5g  %13.5g  %13.5g  %13.5g\n" fname vname an fd abs_err rel_err
                end
            end
        end
        f_off += pf.n_nlp
    end

    # ── Objective gradient dump ───────────────────────────────────────────────
    if !isnothing(p.objective)
        g_an = zeros(n_x)
        for (v, r) in zip(vlist, rngs)
            g_an[r] .= objective_gradient_chunk(p, v)
        end

        g_fd = zeros(n_x)
        for j in 1:n_x
            xp = copy(x0);  xp[j] += h
            xm = copy(x0);  xm[j] -= h
            set_decision_vector!(p, xp);  Jp = get_objective(p)
            set_decision_vector!(p, xm);  Jm = get_objective(p)
            g_fd[j] = (Jp - Jm) / (2h)
        end
        set_decision_vector!(p, x0)

        println()
        println("── Objective gradient dump ───────────────────────────────────────────────────────")
        @printf "  %-16s  %13s  %13s  %13s  %13s\n" "variable" "analytic" "fd" "abs_err" "rel_err"
        for (v, r) in zip(vlist, rngs)
            for (k, gcol) in enumerate(r)
                an      = g_an[gcol]
                fd      = g_fd[gcol]
                abs_err = abs(an - fd)
                rel_err = abs_err / max(abs(fd), 1e-8)
                vname   = length(r) == 1 ? v.name : "$(v.name)[$k]"
                @printf "  %-16s  %13.5g  %13.5g  %13.5g  %13.5g\n" vname an fd abs_err rel_err
            end
        end
    end
    println()
end
