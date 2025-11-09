using Test
using LinearAlgebra
using AstroFun
using AstroUniverse

@testset "GravParam (BodyCalc)" begin

    testbody = CelestialBody("TestBody",     1.32712440018e11, 696342.0, 0.0,         10)

    c = BodyCalc(testbody, GravParam())

    # get_calc returns scalar mu
    @test isapprox(get_calc(c), 1.32712440018e11; atol=1e-8)

    # set updates body.mu and get_calc reflects the change
    set_calc!(c, 4.0e14)
    @test testbody.mu == 4.0e14
    @test get_calc(c) == 4.0e14

    # Optional: API flag
    @test calc_is_settable(GravParam())

    # Test trait methods
    @test calc_numvars(c.var) == 1
    @test calc_is_settable(c.var) == true
end