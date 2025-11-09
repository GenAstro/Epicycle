using Test
using AstroSolve
using AstroCallbacks
using AstroManeuvers
using AstroEpochs
using AstroStates
using AstroUniverse
using AstroModels: Spacecraft

# Helpers
make_sat() = Spacecraft(
    state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]),
    time  = Time("2020-01-01T00:00:00", TAI(), ISOT()),
)

make_impulse() = ImpulsiveManeuver(axes=Inertial(), element1=0.1, element2=0.2, element3=0.3)

@testset "is_astrosolve_stateful trait" begin
    @test AstroSolve.is_astrosolve_stateful(Int) == false
    @test AstroSolve.is_astrosolve_stateful(Spacecraft) == true
    @test AstroSolve.is_astrosolve_stateful(ImpulsiveManeuver) == true
end

@testset "SolverVariable show prints fields" begin
    sat = make_sat()
    man = make_impulse()
    calc = ManeuverCalc(man, sat, DeltaVVector())
    sv = SolverVariable(calc=calc, lower_bound=[-1.0,-1.0,-1.0], upper_bound=[1.0,1.0,1.0], shift=[0,0,0], scale=[1,1,1], name="dv")
    io = IOBuffer()
    show(io, sv)
    s = String(take!(io))
    @test occursin("SolverVariable", s)
    @test occursin("numvars", s)
    @test occursin("name", s)
end

@testset "set_sol_var error for non-settable calc" begin
    sat = make_sat()
    man = make_impulse()
    calc = ManeuverCalc(sat, man, DeltaVMag())
    sv = SolverVariable(calc=calc, name="dv_mag")
    @test_throws ArgumentError set_sol_var(sv, [42.0])
end

@testset "topo_sort detects cycle" begin
    e1 = Event(name="A", event=()->nothing, vars=[], funcs=[])
    e2 = Event(name="B", event=()->nothing, vars=[], funcs=[])
    seq = Sequence()
    add_events!(seq, e1, [e2])  # A depends on B
    add_events!(seq, e2, [e1])  # B depends on A -> cycle
    @test_throws ErrorException topo_sort(seq)
end


nothing