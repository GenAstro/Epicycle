# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

using AstroSolve: func_eval

# Three quantities, one Constraint each. Moved here from AstroCallbacks with the
# type and rewritten for the quantity form when the calc constructors were
# retired. Names are file-local so this does not shadow the `sat` the rest of
# the AstroSolve suite shares.

_con_posvel = [7000.0, 300.0, -412.0, 0.0, 7.5, 0.03]
_con_sat = Spacecraft(
    state=CartesianState(_con_posvel),
)

@testset "Constraint on position magnitude" begin
    c = Constraint(position_magnitude, _con_sat; equals = 15000.0, scale = 6378.0)

    @test c.lower_bound == [15000.0]
    @test c.upper_bound == [15000.0]
    @test c.scale == [6378.0]
    @test c.numvars == 1

    @test isapprox(func_eval(c), [position_magnitude(_con_sat)]; atol=1e-14)
end

@testset "Constraint on velocity magnitude" begin
    c = Constraint(velocity_magnitude, _con_sat; lower_bound = 7.0, upper_bound = 8.0)

    @test c.lower_bound == [7.0]
    @test c.upper_bound == [8.0]
    @test c.scale == [1.0]          # scale defaults to one per component
    @test c.numvars == 1

    @test isapprox(func_eval(c), [velocity_magnitude(_con_sat)]; atol=1e-14)
end

@testset "Constraint on eccentricity" begin
    c = Constraint(eccentricity, _con_sat; lower_bound = 0.0, upper_bound = 0.9)

    @test c.lower_bound == [0.0]
    @test c.upper_bound == [0.9]
    @test c.numvars == 1

    @test isapprox(func_eval(c), [eccentricity(_con_sat)]; atol=1e-14)
end
