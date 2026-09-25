# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# =============================================================================
# Output partials, against a finite difference.
#
# The property under test is that a declared d(quantity)/d(state) is the
# derivative of the reader beside it. A wrong analytic partial is the worst
# failure mode in the framework: it does not raise, it does not return zero,
# it moves the solver in a direction that is merely wrong, and the run either
# converges somewhere else or does not converge at all.
#
# Central differences, so the error is O(h^2) and a 1e-6 step leaves room
# between truncation and round-off.
# =============================================================================

using Test
using AstroCallbacks
using AstroModels
using AstroStates
using AstroEpochs
using AstroFrames
using AstroUniverse
using EpicycleBase: output_partial, has_output_partial, tag
using LinearAlgebra: norm

_op_epoch() = Time("2020-09-21T12:23:12", TAI(), ISOT())

# A state with every component non-zero and non-degenerate, so no term drops out
# of the comparison by accident.
_op_state() = [7100.0, 320.0, -410.0, 0.35, 7.42, 1.15]
_op_sat()   = Spacecraft(state = CartesianState(_op_state()), time = _op_epoch())

"Central-difference d(quantity)/d(state) for a subject built from a state vector."
function _fd_partial(quantity, x0; h = 1e-6)
    cols = map(1:6) do j
        step = h * max(abs(x0[j]), 1.0)
        xp = copy(x0); xp[j] += step
        xm = copy(x0); xm[j] -= step
        fp = quantity(Spacecraft(state = CartesianState(xp), time = _op_epoch()))
        fm = quantity(Spacecraft(state = CartesianState(xm), time = _op_epoch()))
        (collect(fp) .- collect(fm)) ./ (2 * step)
    end
    return hcat(cols...)
end

const _OP_QUANTITIES = (
    (position_vector,       "position_vector",       3),
    (velocity_vector,       "velocity_vector",       3),
    (position_x,            "position_x",            1),
    (position_y,            "position_y",            1),
    (position_z,            "position_z",            1),
    (position_magnitude,    "position_magnitude",    1),
    (velocity_magnitude,    "velocity_magnitude",    1),
    (position_dot_velocity, "position_dot_velocity", 1),
    (semi_major_axis,       "semi_major_axis",       1),
)

@testset "output partials match a finite difference" begin
    x0  = _op_state()
    sat = _op_sat()
    for (q, name, nrows) in _OP_QUANTITIES
        @testset "$name" begin
            @test has_output_partial(sat, q)

            analytic = output_partial(sat, q)
            @test size(analytic) == (nrows, 6)

            numeric = _fd_partial(q, x0)
            scale   = max(maximum(abs, numeric), 1.0)
            @test maximum(abs, analytic .- numeric) / scale < 1e-6
        end
    end
end

@testset "a quantity without an output partial says so" begin
    sat = _op_sat()
    # Declared for none of these; they fall to automatic differentiation.
    @test !has_output_partial(sat, eccentricity)
    @test !has_output_partial(sat, inclination)
    @test !has_output_partial(sat, true_anomaly)

    # A scalar field on a model takes the other path entirely: it has a tag and
    # its derivative comes from param_jac!, not from here.
    @test !has_output_partial(earth, gravitational_parameter)
end
