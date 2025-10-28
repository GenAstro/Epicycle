
using SNOW
using OrdinaryDiffEq
using LinearAlgebra

using Epicycle

# Create spacecraft
sat = Spacecraft(
            state = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 1.0]), 
            time = Time("2020-09-21T12:23:12", TAI(), ISOT()),
            )

# Create the propagator
gravity = PointMassGravity(earth,(moon,sun))
forces  = ForceModel(gravity)
integ   = IntegratorConfig(Tsit5(); dt=10.0, reltol=1e-9, abstol=1e-9)
prop    = OrbitPropagator(forces, integ)

# Create two impulsive maneuvers 
toi = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 1.86,
    element2 = 0.0,
    element3 = 0.0
)

moi = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 0.5,
    element2 = 0.0,
    element3 = 0.0
)

# Define maneuvers as solver variables
var_toi = SolverVariable(
    calc = ManeuverCalc(toi, sat, DeltaVVector()),
    name = "toi",
    lower_bound = [-10.0, 0.0, 0.0],
    upper_bound = [10.0, 0.0, 0.0],
)

var_moi = SolverVariable(
    calc = ManeuverCalc(moi, sat, DeltaVVector()),
    name = "moi", 
    lower_bound = [-1.0, 0.0, 0.0],
    upper_bound = [1.0, 0.0, 0.0],
)

# Define constraints on position magnitude at two different points
# These should have DIFFERENT values if constraints are evaluated correctly
pos_target_1 = 25000.0  # Target after first propagation
pos_target_2 = 55000.0  # Target after second propagation

pos_con_1 = Constraint(
    calc = OrbitCalc(sat, PosMag()),
    lower_bounds = [pos_target_1],
    upper_bounds = [pos_target_1],
    scale = [1.0],
)

pos_con_2 = Constraint(
    calc = OrbitCalc(sat, PosMag()),
    lower_bounds = [pos_target_2], 
    upper_bounds = [pos_target_2],
    scale = [1.0],
)

# Create Events
# Event 1: Apply TOI maneuver
fun_toi() = maneuver(sat, toi) 
toi_event = Event(name = "toi", 
                  event = fun_toi,
                  vars = [var_toi],
                  funcs = [])

# Event 2: Propagate to apoapsis and check first constraint
fun_prop_apo_1() = propagate(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=-1))
prop_event_1 = Event(name = "prop_apo_1", 
                     event = fun_prop_apo_1,
                     funcs = [pos_con_1])

# Event 3: Apply MOI maneuver  
fun_moi() = maneuver(sat, moi)
moi_event = Event(name = "moi",
                  event = fun_moi,
                  vars = [var_moi],
                  funcs = [pos_con_2])

# Event 4: Propagate to apoapsis again and check second constraint
fun_prop_apo_2() = propagate(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=-1))
prop_event_2 = Event(name = "prop_apo_2",
                     event = fun_prop_apo_2, 
                     funcs = [pos_con_2])

# Build sequence with dependencies
seq = Sequence()
add_events!(seq, prop_event_1, [toi_event])      # prop_1 after toi
add_events!(seq, moi_event, [prop_event_1])      # moi after prop_1  
add_events!(seq, prop_event_2, [moi_event])      # prop_2 after moi 

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
println("Final Function Values: ", fopt)
