
using SNOW
using OrdinaryDiffEq
using LinearAlgebra

using AstroBase
using AstroStates
using AstroEpochs
using AstroUniverse
using AstroFrames
using AstroModels
using AstroManeuvers
using AstroCallbacks
using AstroProp
using AstroSolve

# Create spacecraft
sat = Spacecraft(
            state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0]), 
            time = Time("2020-09-21T12:23:12", TAI(), ISOT())
            )

# Create force models, integrator, and dynamics system
pm_grav = PointMassGravity(earth,(moon,sun))
#pm_grav = PointMassGravity(earth,())
forces = ForceModel(pm_grav)
integ = IntegratorConfig(DP8(); abstol = 1e-11, reltol = 1e-11, dt = 4000)

# Define which spacecraft to propagate and which force model to use
dynsys = DynSys(
          forces = forces, 
          spacecraft = [sat]
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
    calc = ManeuverCalc(toi, sat, DeltaVVector()),
    name = "toi",
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
)

ecc_con = Constraint(
    calc = OrbitCalc(sat, Ecc()),
    lower_bounds = [0.0],
    upper_bounds = [0.0], 
    scale = [1.0],
)

vel_target = sqrt(earth.mu / pos_target)
vel_con = Constraint(
    calc = OrbitCalc(sat, VelMag()),
    lower_bounds = [vel_target],
    upper_bounds = [vel_target],
    scale = [1.0],
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

x0 = get_var_values(sm)
lx = get_var_lower_bounds(sm)
ux = get_var_upper_bounds(sm)
lg = get_fun_lower_bounds(sm)
ug = get_fun_upper_bounds(sm)
ng = length(lg)

ip_options = Dict(
        "max_iter" => 1000,
        "tol" => 1e-6,
        "output_file" => "ipopt_$(rand(UInt)).out",
        "file_print_level" => 0,
        )

options = Options(derivatives= ForwardFD(), solver = IPOPT(ip_options))

# Define a closure that matches SNOW's expected signature
snow_solver_fun!(F, x) = solver_fun!(F, x, sm)
xopt, fopt, info = minimize(snow_solver_fun!, x0, ng, lx, ux, lg, ug, options)

println("Optimal Variables: ", xopt)