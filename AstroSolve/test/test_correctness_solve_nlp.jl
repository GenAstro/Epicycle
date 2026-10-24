# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# The nonlinear program solver every mission-design solve goes through: `solve_nlp` and the
# derivatives it hands IPOPT.
#
# Truth, for the solves: **closed form**. Minimizing (x1 - 1)² + (x2 - 2)² on x1 + x2 = 2 gives
# x = (0.5, 1.5) and f = 0.5, with the bound x1 ≥ 0.2 inactive. Each derivative route reaches it.
#
# Truth, for the finite differences: the analytic derivatives of a function chosen so no entry is
# constant, compared to forward-difference accuracy.

using Test
using LinearAlgebra
using AstroSolve
using AstroSolve: NLPOptions, JacobianPattern, solve_nlp, _ipopt, _evaluator

_nlp_ipopt(; kw...) = _ipopt(merge(Dict{String,Any}("print_level" => 0, "tol" => 1e-10,
                                                    "output_file" => tempname() * "_ipopt.out"),
                                   Dict{String,Any}(String(k) => v for (k, v) in kw));
                             monitor = false)

_nlp_f(x) = (x[1] - 1)^2 + (x[2] - 2)^2
_nlp_g!(g, x) = (g[1] = x[1] + x[2]; g[2] = x[1]; g)

function _nlp_fd!(g, x)
    _nlp_g!(g, x)
    return _nlp_f(x)
end

function _nlp_dense!(g, df, dg, x)
    _nlp_g!(g, x)
    df .= [2(x[1] - 1), 2(x[2] - 2)]
    dg .= [1.0 1.0; 1.0 0.0]                # constraints by variables
    return _nlp_f(x)
end

# Only the three nonzeros, as (row, column): (1, 1), (2, 1), (1, 2).
const _NLP_PATTERN = JacobianPattern([1, 2, 1], [1, 1, 2])

function _nlp_sparse!(g, df, dg, x)
    _nlp_g!(g, x)
    df .= [2(x[1] - 1), 2(x[2] - 2)]
    dg .= [1.0, 1.0, 1.0]
    return _nlp_f(x)
end

_nlp_solve(fun!, options) =
    solve_nlp(fun!, [0.0, 0.0], 2, [-10.0, -10.0], [10.0, 10.0], [2.0, 0.2], [2.0, Inf], options)

@testset "solve_nlp" begin

    @testset "every derivative route reaches the known optimum" begin
        # Each route at a tolerance its derivatives can meet. Forward differences give gradients to
        # about √eps, 1e-8, so IPOPT cannot drive the dual infeasibility below that and reports an
        # acceptable solve at tol = 1e-10; analytic derivatives meet 1e-10.
        cases = (("forward differences", _nlp_fd!,     1e-6,
                  NLPOptions(derivatives = :fd, solver = _nlp_ipopt(tol = 1e-6))),
                 ("user, dense",         _nlp_dense!,  1e-8,
                  NLPOptions(derivatives = :user, solver = _nlp_ipopt())),
                 ("user, sparse",        _nlp_sparse!, 1e-8,
                  NLPOptions(derivatives = :user, sparsity = _NLP_PATTERN, solver = _nlp_ipopt())))
        for (name, fun!, atol, options) in cases
            @testset "$name" begin
                x, f, status = _nlp_solve(fun!, options)
                @test status === :Solve_Succeeded
                @test x ≈ [0.5, 1.5] atol = atol
                @test f ≈ 0.5 atol = atol
            end
        end
    end

    @testset "forward differences match the analytic derivatives" begin
        f(x) = sin(x[1]) * exp(x[2])
        fun!(g, x) = (g[1] = x[1]^2 * x[2]; g[2] = cos(x[2]) + x[1]; f(x))
        x  = [0.3, -1.2]
        df_true = [cos(x[1]) * exp(x[2]), sin(x[1]) * exp(x[2])]
        J_true  = [2x[1]*x[2]  x[1]^2;
                   1.0         -sin(x[2])]

        eval! = _evaluator(fun!, NLPOptions(derivatives = :fd, solver = _nlp_ipopt()), 2, 2)
        g, df, dg = zeros(2), zeros(2), zeros(4)
        @test eval!(g, df, dg, x) == f(x)
        @test g == [x[1]^2 * x[2], cos(x[2]) + x[1]]
        # A forward difference is accurate to about √eps times the second derivative's size.
        @test df ≈ df_true atol = 1e-7
        @test dg ≈ vec(J_true) atol = 1e-7              # column-major, the dense positions' order
    end

    @testset "a dense Jacobian is handed over column by column" begin
        J = [1.0 2.0 3.0; 4.0 5.0 6.0]
        fun!(g, df, dg, x) = (g .= 0; df .= 0; dg .= J; 0.0)
        eval! = _evaluator(fun!, NLPOptions(derivatives = :user, solver = _nlp_ipopt()), 3, 2)
        dg = zeros(6)
        eval!(zeros(2), zeros(3), dg, zeros(3))
        @test dg == [1.0, 4.0, 2.0, 5.0, 3.0, 6.0]
        rows, cols = AstroSolve._dense_positions(2, 3)
        @test [J[r, c] for (r, c) in zip(rows, cols)] == dg
    end

    @testset "options refuse what the solver cannot do" begin
        @test_throws ArgumentError NLPOptions(derivatives = :central, solver = _nlp_ipopt())
        @test_throws ArgumentError NLPOptions(derivatives = :fd, sparsity = _NLP_PATTERN,
                                              solver = _nlp_ipopt())
    end

    @testset "IPOPT's exit status is reported by name" begin
        _, _, status = _nlp_solve(_nlp_fd!, NLPOptions(derivatives = :fd,
                                                       solver = _nlp_ipopt(max_iter = 1)))
        @test status === :Maximum_Iterations_Exceeded
    end
end
