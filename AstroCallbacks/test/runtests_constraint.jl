posvel = [7000.0, 300.0, -412.0, 0.0, 7.5, 0.03]
sat = Spacecraft(
    state=CartesianState(posvel),
)

@testset "Constraint constructor with PosMag OrbitCalc" begin
    calc = OrbitCalc(sat, PosMag())
    lower_bounds = [15000.0]
    upper_bounds = [15000.0]
    scale = [6378.0]

    c = Constraint(
        calc = calc,
        lower_bounds = lower_bounds,
        upper_bounds = upper_bounds,
        scale = scale,
    )

    @test c.calc === calc
    @test c.lower_bounds == [15000.0]
    @test c.upper_bounds == [15000.0]
    @test c.scale == [6378.0]
    @test c.numvars == 1

    val = func_eval(c)
    @test isapprox(val, [get_calc(calc)]; atol=1e-14)
end

@testset "Constraint constructor with VelMag OrbitCalc" begin
    calc = OrbitCalc(sat, VelMag())
    lower_bounds = [7.0]
    upper_bounds = [8.0]
    scale = [1.0]

    c = Constraint(
        calc = calc,
        lower_bounds = lower_bounds,
        upper_bounds = upper_bounds,
        scale = scale,
    )

    @test c.calc === calc
    @test c.lower_bounds == [7.0]
    @test c.upper_bounds == [8.0]
    @test c.scale == [1.0]
    @test c.numvars == 1

    val = func_eval(c)
    @test isapprox(val, [get_calc(calc)]; atol=1e-14)
end

@testset "Constraint constructor with Ecc OrbitCalc" begin
    calc = OrbitCalc(sat, Ecc())
    lower_bounds = [0.0]
    upper_bounds = [0.9]
    scale = [1.0]

    c = Constraint(
        calc = calc,
        lower_bounds = lower_bounds,
        upper_bounds = upper_bounds,
        scale = scale,
    )

    @test c.calc === calc
    @test c.lower_bounds == lower_bounds
    @test c.upper_bounds == upper_bounds
    @test c.scale == scale
    @test c.numvars == 1

    val = func_eval(c)
    @test isapprox(val, [get_calc(calc)]; atol=1e-14)
end
