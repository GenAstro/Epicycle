using LinearAlgebra
using AstroFrames
using AstroProp
using AstroEpochs
using AstroStates
using AstroManeuvers
using AstroSolve
using AstroUniverse

# Create spacecraft
sat1 = Spacecraft(
            state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0]), 
            time = Time("2020-09-21T12:23:12", TAI(), ISOT())
            )

sat2 = Spacecraft(
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

dynsys2 = DynSys(
          forces = forces, 
          spacecraft = [sat2]
          )
          
# Create maneuver models for the hohmann transfer
toi = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 0.1,
    element2 = 0.2,
    element3 = 0.3
)

moi = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 0.4,
    element2 = 0.5,
    element3 = 0.6
)

# Define toi as a solver variable
var_toi = SolverVariable(
    calc = ManeuverCalc(toi, sat1, DeltaVVector()),
    name = "toi",
    lower_bound = [-1.0, -2.0, -3.0],
    upper_bound = [1.0, 2.0, 3.0],
    shift = [0.01, 0.02, 0.03],
    scale = [1.1, 1.2, 1.3]
)

# Define moi as a solver variable
var_moi = SolverVariable(
    calc = ManeuverCalc(moi, sat1, DeltaVVector()),
    name = "moi",
    lower_bound = [-4.0, -5.0, -6.0],
    upper_bound = [4.0, 5.0, 6.0],
    shift = [0.04, 0.05, 0.06],
    scale = [1.4, 1.5, 1.6]
)

pos_target = 45000.0
pos_con = Constraint(
    calc = OrbitCalc(sat1, PosMag()),
    lower_bounds = [pos_target],
    upper_bounds =[pos_target],
    scale = [1.0],
)

vel_target = sqrt(earth.mu / pos_target)
vel_con = Constraint(
    calc = OrbitCalc(sat1, VelMag()),
    lower_bounds = [vel_target],
    upper_bounds = [vel_target],
    scale = [1.0],
)

# Create the MOI Event
toi_fun() = maneuver!(sat1, toi) 
toi_event = Event(name = "toi", 
                  event = toi_fun, 
                  vars = [var_toi],
                  funcs = []
                  )

# Create the prop to apopasis event
prop_apo_fun() = propagate!(dynsys, integ, StopAtApoapsis(sat1))
prop_event = Event(name = "prop_apo", event = prop_apo_fun)

# Create the TOI event. 
moi_fun() = maneuver!(sat1, moi)

moi_event = Event(event = moi_fun, vars = [var_moi],
                  funcs = [pos_con, vel_con])

# Build sequence and solve
seq = Sequence()
add_events!(seq, prop_event, [toi_event]) 
add_events!(seq, moi_event, [prop_event])

sm = SequenceManager(seq)
