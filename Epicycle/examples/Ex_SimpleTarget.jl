using Epicycle

# ============================================================================
# Shared Resources
# ============================================================================

# Create spacecraft with default orbital state
sat = Spacecraft()

# Create force models, integrator, and propagator
gravity = PointMassGravity(earth, (moon, sun))
forces  = ForceModel(gravity)
integ   = IntegratorConfig(Tsit5(); dt=10.0, reltol=1e-9, abstol=1e-9)
prop    = OrbitPropagator(forces, integ)

# ============================================================================
# Event 1: TOI - Transfer Orbit Insertion
# ============================================================================

# Define TOI maneuver
toi = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 0.1,
)

# Define solver variable for TOI delta-V
toi_var = SolverVariable(
    calc = ManeuverCalc(toi, sat, DeltaVVector()),
    name = "toi",
    lower_bound = [-10.0, 0.0, 0.0],
    upper_bound = [10.0, 0.0, 0.0],
)

# Define TOI event struct with event function, solver variables, and constraints
toi_fun() = maneuver!(sat, toi) 
toi_event = Event(
    name = "TOI Maneuver", 
    event = toi_fun,
    vars = [toi_var],
    funcs = []
)

# ============================================================================
# Event 2: Propagate to Apoapsis
# ============================================================================

# Define constraint on position magnitude at apoapsis
pos_target = 55000.0
pos_con = Constraint(
    calc = OrbitCalc(sat, PosMag()),
    lower_bounds = [pos_target],
    upper_bounds = [pos_target],
    scale = [1.0],
)

# Define propagation event to apoapsis with position constraint
prop_fun() = propagate!(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=-1))
prop_event = Event(
    name = "Propagate to Apoapsis", 
    event = prop_fun,
    funcs = [pos_con]
)

# ============================================================================
# Trajectory Optimization
# ============================================================================

# Create sequence and add events
seq = Sequence()
add_sequence!(seq, toi_event, prop_event) 

# Solve trajectory optimization using default settings (finite differences, IPOPT)
result = solve_trajectory!(seq)

# Write a report documenting sequence and solution
report_sequence(seq)
report_solution(seq, result)