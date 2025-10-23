


using SNOW
using OrdinaryDiffEq
using LinearAlgebra

using AstroEpochs
using AstroStates
using AstroBase
using AstroUniverse
using AstroCoords
using AstroProp
using AstroMan
using AstroSolve
using AstroFun

# Create spacecraft
sat = Spacecraft(
            state = KeplerianState(6578.0, 0.001, 28.5, 66.999, 355.0, 250.0), 
            time = Time("2030-07-04T12:23:12", UTC(), ISOT())
            )

# Create force models, integrator, and dynamics system
gravity = PointMassGravity(earth,(moon,sun))
forces  = ForceModel(gravity)
integ   = IntegratorConfig(DP8(); abstol = 1e-11, reltol = 1e-11, dt = 4000)

# Define which spacecraft to propagate and which force model to use
dynsys = DynSys(
          forces = forces, 
          spacecraft = [sat]
          )

# Dfine maneuver models
toi = ImpulsiveManeuver(
    axes = VNB(),      # Local VNB axes about the spacecraft
    element1 = 0.0,        # dv_V
    element2 = 0.0,        # dv_N
    element3 = 0.0,        # dv_B
)

mcc = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 0.0,
    element2 = 0.0,
    element3 = 0.0,
)

moi = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 0.0,
    element2 = 0.0,
    element3 = 0.0,
)

# Define toi as a solver variable
var_toi = SolverVariable(
    calc = ManeuverCalc(toi, sat, DeltaVVector()),
    name = "toi",
    lower_bound = [-10.0, 0.0, 0.0],
    upper_bound = [10.0, 0.0, 0.0],
)

var_mcc = SolverVariable(
    calc = ManeuverCalc(mcc, sat, DeltaVVector()),
    name = "mcc",
    lower_bound = [-10.0, 0.0, 0.0],
    upper_bound = [10.0, 0.0, 0.0],
)

# Define moi as a solver variable
var_moi = SolverVariable(
    calc = ManeuverCalc(moi, sat, DeltaVVector()),
    name = "moi",
    lower_bound = [-10.0, 0.0, 0.0],
    upper_bound = [10.0, 0.0, 0.0]
)

pos_target = 45000.0
pos_con = Constraint(
    calc = OrbitCalc(sat, PosMag()),
    lower_bounds = [pos_target],
    upper_bounds = [pos_target],
    scale = [1.0],
    numvars=1
)

ecc_con = Constraint(
    calc = OrbitCalc(sat, Ecc()),
    lower_bounds = [0.0],
    upper_bounds = [0.0], 
    scale = [1.0],
    numvars = 1
)

vel_target = sqrt(earth.mu / pos_target)
vel_con = Constraint(
    calc = OrbitCalc(sat, VelMag()),
    lower_bounds = [vel_target],
    upper_bounds = [vel_target],
    scale = [1.0],
    numvars = 1
)

# Create the TOI Event
toi_fun() = maneuver(sat, toi) 
toi_event = Event(name = "toi", 
                  event = toi_fun, 
                  vars = [var_toi],
                  funcs = [])

# Create the prop to apopasis event
prop_apo_fun() = propagate(dynsys, integ, StopAtApoapsis(sat))
prop_event = Event(name = "prop_apo", event = prop_apo_fun)

# Create the TOI event. 
moi_fun() = maneuver(sat, moi)
moi_event = Event(event = moi_fun, 
                  vars = [var_moi],
                  funcs = [pos_con, ecc_con])

# Build sequence and solve
seq = Sequence()
add_events!(seq, prop_event, [toi_event]) 
add_events!(seq, moi_event, [prop_event])

sm = SequenceManager(seq)

f = get_fun_values(sm)

x0 = get_var_values(sm)

set_var_values(sm, x0)
