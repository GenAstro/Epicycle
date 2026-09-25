# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# =============================================================================
# Constraint — the spec, and what it evaluates to.
#
# `Constraint` and `Vary` take the same first arguments, so a script states a
# goal and a freedom the same way. The property under test is that the quantity
# form builds the same object the calc form does, that `func_eval` always hands
# the solver a Vector, and that a Constraint missing its goal says so.
#
# These tests came from AstroCallbacks with the type itself.
# =============================================================================

using Test
using AstroCallbacks
using AstroManeuvers
using AstroModels
using AstroStates
using AstroEpochs
using AstroFrames
using AstroUniverse
using AstroSolve
using AstroSolve: func_eval

const _C_EPOCH = Time(2458849.5, 0.0, :tdb, :jd)
_c_sat() = Spacecraft(CartesianState([7000.0, 1000.0, 2000.0, 1.5, 6.5, 1.0]), _C_EPOCH;
                      coord_sys = CoordinateSystem(earth, ICRF()))

@testset "Constraint from a quantity" begin
    sat = _c_sat()

    con = Constraint(position_magnitude, sat; equals = 8000.0)
    @test con.lower_bound == [8000.0]
    @test con.upper_bound == [8000.0]
    @test con.scale        == [1.0]
    @test con.numvars      == 1
    @test func_eval(con)[1] ≈ position_magnitude(sat) rtol = 1e-12

    # The dependency slot works here exactly as it does in StopAt.
    ec  = CoordinateSystem(earth, MJ2000Ec())
    con_ec = Constraint(inclination, sat, ec; equals = 0.1)
    @test func_eval(con_ec)[1] ≈ inclination(sat, ec) rtol = 1e-12
    @test func_eval(con_ec)[1] != inclination(sat)     # genuinely a different frame

    # An inequality is the other form.
    con_ineq = Constraint(position_magnitude, sat; lower_bound = 7000.0, upper_bound = 9000.0)
    @test con_ineq.lower_bound == [7000.0]
    @test con_ineq.upper_bound == [9000.0]

    # A vector quantity gets one bound per component.
    con_vec = Constraint(position_vector, sat; equals = 0.0)
    @test con_vec.numvars == 3
    @test length(con_vec.lower_bound) == 3
end

@testset "a Constraint says what is missing" begin
    sat = _c_sat()
    e1 = try; Constraint(position_magnitude, sat); nothing; catch e; e; end
    @test e1 isa ArgumentError
    @test occursin("equals", e1.msg)

    e2 = try; Constraint(position_magnitude, sat; equals = 1.0, lower_bound = 2.0); nothing; catch e; e; end
    @test e2 isa ArgumentError
    @test occursin("not both", e2.msg)
end

# Rationale: func_eval always returns a Vector, for scalar and vector-valued calcs.
@testset "func_eval output normalization" begin
    sc = Spacecraft(
        state=CartesianState([7000.0,300.0,0.0, 0.0,7.5,1.0]),
        time=Time("2020-01-01T00:00:00", TAI(), ISOT()),
    )
    c_vec = Constraint(position_vector, sc; lower_bound = [-1.0,-1.0,-1.0], upper_bound = [1.0,1.0,1.0], scale = [1.0,1.0,1.0])
    v = func_eval(c_vec)
    @test v == [7000.0, 300.0, 0.0]

    c_sca = Constraint(gravitational_parameter, earth; lower_bound = 0.0, upper_bound = 1e7)
    s = func_eval(c_sca)
    @test s isa Vector
    @test length(s) == 1
end

# Rationale: func_eval errors on a quantity that is neither a number nor a vector.
# A user writes a quantity as a function, so that is what the bad case is here.
_matrix_quantity(_subject) = [1.0 2.0; 3.0 4.0]

@testset "func_eval rejects a quantity that is not a number or a vector" begin
    c = Constraint(_matrix_quantity, _c_sat(); equals = 0.0)
    try
        func_eval(c)
        @test false
    catch e
        msg = sprint(showerror, e)
        # §9.4: the message names the constraint the quantity failed, why it exists, and the
        # type that was actually returned. The wording was tightened in the 2026-09-15 scrub;
        # what is pinned here is the content, not the old phrasing.
        @test occursin("Number or an AbstractVector", msg)
        @test occursin("Matrix", msg)
    end
end
