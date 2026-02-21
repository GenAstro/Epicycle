
using Test
using LinearAlgebra
using Epicycle

# Create spacecraft
sat1 = Spacecraft(
            state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0]), 
            time = Time("2020-09-21T12:23:12", TAI(), ISOT())
            )

# Create force models, integrator, and dynamics system
pm_grav = PointMassGravity(earth,(moon,sun))
forces = ForceModel(pm_grav)
integ = IntegratorConfig(DP8(); abstol = 1e-11, reltol = 1e-11, dt = 4000)

# Define which spacecraft to propagate and which force model to use
dynsys = DynSys(
          forces = forces, 
          spacecraft = [sat1]
          )

sat2 = deepcopy(sat1)
dynsys2 = DynSys(
          forces = forces, 
          spacecraft = [sat2]
          )

# Create maneuver models for the hohmann transfer
toi = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 0.4,
    element2 = 0.0,
    element3 = 0.3
)

moi = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 0.1,
    element2 = 0.0,
    element3 = 0.0
)

# Define toi as a solver variable
var_toi = SolverVariable(
    calc = ManeuverCalc(toi, sat1, DeltaVVector()),
    name = "toi",
    lower_bound = [-1.0, 0.0, 0.0],
    upper_bound = [1.0, 0.0, 0.0]
)

# Define moi as a solver variable
var_moi = SolverVariable(
    calc = ManeuverCalc(moi, sat1, DeltaVVector()),
    name = "moi",
    lower_bound = [-1.0, 0.0, 0.0],
    upper_bound = [1.0, 0.0, 0.0]
)

@testset "Event Propagation Test" begin

    # === Define propagation function ===
    prop_to_moi() = propagate!(dynsys, integ, StopAtApoapsis(sat1))

    propagate!(dynsys2, integ, StopAtApoapsis(sat2))

    # === Create the Event ===
    prop_moi = Event(event = prop_to_moi)
    apply_event(prop_moi)

    @test length(prop_moi.vars) == 0
    @test length(prop_moi.funcs) == 0

    # === Run the event, just confirm it does not throw ===
    @test begin
        #event_prop = CartesianState(sat1.state).posvel  
        event_prop = to_posvel(sat1)
        #truth_prop = CartesianState(sat2.state).posvel 
        truth_prop = to_posvel(sat2)      
        isapprox(event_prop, truth_prop; rtol=1e-14, atol=1e-14)
    end

end
