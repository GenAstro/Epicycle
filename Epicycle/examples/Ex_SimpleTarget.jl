using Epicycle

# Create spacecraft
sat = Spacecraft(
            state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0]), 
            time = Time("2020-09-21T12:23:12", TAI(), ISOT()),
            )

# Create the propagator
gravity = PointMassGravity(earth,(moon,sun))
forces  = ForceModel(gravity)
integ   = IntegratorConfig(Tsit5(); dt=10.0, reltol=1e-9, abstol=1e-9)
prop    = OrbitPropagator(forces, integ)

# Create an impulsive maneuver 
toi = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 0.1,
    element2 = 0.2,
    element3 = 0.3
)

# Define toi as a solver variable
var_toi = SolverVariable(
    calc = ManeuverCalc(toi, sat, DeltaVVector()),
    name = "toi",
    lower_bound = [-10.0, 0.0, 0.0],
    upper_bound = [10.0, 0.0, 0.0],
)

# Define a constraint on magnitude of orbit position
pos_target = 55000.0
pos_con = Constraint(
    calc = OrbitCalc(sat, PosMag()),
    lower_bounds = [pos_target],
    upper_bounds = [pos_target],
    scale = [1.0],
)

# Create the TOI Event that applies the maneuver.
# The event is a maneuver and the variable is the toi variable
fun_toi() = maneuver(sat, toi) 
toi_event = Event(name = "toi", 
                  event = fun_toi,
                  vars = [var_toi],
                  funcs = [])

# Create the prop to apopasis event
fun_prop_apo() = propagate(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=-1))
prop_event = Event(name = "prop_apo", 
                   event = fun_prop_apo,
                   funcs = [pos_con])

# Build sequence and solve
seq = Sequence()
add_events!(seq, prop_event, [toi_event]) 

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