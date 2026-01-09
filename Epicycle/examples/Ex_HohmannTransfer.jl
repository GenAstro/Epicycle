
using Epicycle

# ============================================================================
# Shared Resources
# ============================================================================

# Create spacecraft
sat = Spacecraft(
            state = KeplerianState(7000.0, 0.001, 0.0, 0.0, 7.5, 1.0), 
            time = Time("2020-09-21T12:23:12", TAI(), ISOT()),
            name = "Sat",
            )

# Create force models, integrator, and dynamics system
gravity = PointMassGravity(earth, (moon,sun))  
forces  = ForceModel(gravity)
integ   = IntegratorConfig(DP8(); abstol=1e-9, reltol=1e-9, dt=300.0)
prop    = OrbitPropagator(forces, integ)

# ============================================================================
# TOI Event - Transfer Orbit Insertion
# ============================================================================

# Define TOI maneuver
toi = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 0.1,
    element2 = 0.2,
    element3 = 0.3
)

# Define solver variable for TOI delta-V
toi_var = SolverVariable(
    calc = ManeuverCalc(toi, sat, DeltaVVector()),
    name = "toi",
    lower_bound = [0.0, 0.0, 0.0],
    upper_bound = [2.5, 0.0, 0.0],
)

# Define TOI event struct with event function, solver variables, and constraints
toi_fun() = maneuver(sat, toi) 
toi_event = Event(
    name = "TOI", 
    event = toi_fun, 
    vars = [toi_var],
    funcs = []
)

# ============================================================================
# Propagation Event - Coast to Apoapsis
# ============================================================================

prop_fun() = propagate(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=-1))
prop_event = Event(
    name = "Prop to Apoapsis", 
    event = prop_fun
)

# ============================================================================
# MOI Event - Mars Orbit Insertion
# ============================================================================

# Define MOI maneuver 
moi = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 0.4,
    element2 = 0.5,
    element3 = 0.6
)

# Define solver variable for MOI delta-V
moi_var = SolverVariable(
    calc = ManeuverCalc(moi, sat, DeltaVVector()),
    name = "moi",
    lower_bound = [0.0, 0.0, 0.0],
    upper_bound = [3.0, 0.0, 0.0]
)

# Define constraints for MOI event
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

# Define MOI event struct with event function, solver variables, and constraints
moi_fun() = maneuver(sat, moi)
moi_event = Event(
    name = "MOI", 
    event = moi_fun,
    vars = [moi_var],
    funcs = [pos_con, ecc_con]
)

# ============================================================================
# Trajectory Optimization
# ============================================================================

# Build sequence and solve
seq = Sequence()
add_sequence!(seq, toi_event, prop_event, moi_event)

# Solve the sequence and report the solution
result = trajectory_solve(seq; record_iterations=true)
sequence_report(seq)
solution_report(seq, result)

# Plot the trajectory 3D
view = View3D()
add_spacecraft!(view, sat; show_iterations=true)
display_view(view)
