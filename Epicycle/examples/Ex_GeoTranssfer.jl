
using Epicycle

# ========== Create Models =================================================================
# Create spacecraft
sat = Spacecraft(
    state = CartesianState([3737.792, -4607.692, -2845.644, 5.411, 5.367, -1.566]),
    time = Time("2000-01-01T11:59:28.000", UTC(), ISOT())  # J2000 epoch from GMAT
)

# Create simple Earth point mass dynamics (no third bodies for now)
gravity = PointMassGravity(earth, ())  # Only Earth gravity
forces  = ForceModel(gravity)
integ   = IntegratorConfig(DP8(); abstol=1e-12, reltol=1e-12, dt=60.0)
prop    = OrbitPropagator(forces, integ)

# Define maneuver models
toi = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 1.518,
    element2 = 0.0,
    element3 = 0.0,
)

mcc = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 0.559,
    element2 = 0.588,
    element3 = 0.0,
)

moi = ImpulsiveManeuver(
    axes = VNB(),
    element1 = -0.282,
    element2 = 0.0,
    element3 = 0.0,
)

# ========== Define Solver Variables ======================================================
# TOI: Only vary V component (Element1)
var_toi_v = SolverVariable(
    calc = ManeuverCalc(toi, sat, DeltaVVector()),
    name = "toi_v",
    lower_bound = [-5.0, 0.0, 0.0],  # Only V component varies significantly
    upper_bound = [5.0, 0.0, 0.0],
)

# MCC: Vary V and N components (Element1 and Element2)
var_mcc_vn = SolverVariable(
    calc = ManeuverCalc(mcc, sat, DeltaVVector()),
    name = "mcc_vn",
    lower_bound = [-2.0, -2.0, -0.001],   # V and N components vary
    upper_bound = [2.0, 2.0, 0.001],
)

# MOI: Only vary V component (Element1)
var_moi_v = SolverVariable(
    calc = ManeuverCalc(moi, sat, DeltaVVector()),
    name = "moi_v",
    lower_bound = [-2.0, -0.001, -0.001],  # Only V component varies significantly
    upper_bound = [2.0, 0.001, 0.001],
)

# ========== Define Constraints ===========================================================
# Constraint: Apoapsis radius = 85,000 km after TOI
apogee_radius_con = Constraint(
    calc = OrbitCalc(sat, PosMag()),
    lower_bounds = [85000.0],
    upper_bounds = [85000.0],
    scale = [1.0],
)

# Constraint: Inclination = 2° after MCC
inclination_con = Constraint(
    calc = OrbitCalc(sat, Inc()),
    lower_bounds = [deg2rad(2.0)],
    upper_bounds = [deg2rad(2.0)],
    scale = [1.0],
)

# Constraint: Perigee radius = 42,195 km after MCC
perigee_radius_con = Constraint(
    calc = OrbitCalc(sat, PosMag()),
    lower_bounds = [42195.0],
    upper_bounds = [42195.0],
    scale = [1.0],
)

# Constraint: Final SMA = 42,166.90 km (GEO) after MOI
final_sma_con = Constraint(
    calc = OrbitCalc(sat, SMA()),
    lower_bounds = [42166.90],
    upper_bounds = [42166.90],
    scale = [1.0],
)

# ========== Define events (variables, constraints, actions) ==============================
# Event 1: Propagate to Z=0 crossing (equatorial plane)
prop_to_z_crossing_1_fun() = propagate(prop, sat, StopAt(sat, PosZ(), 0.0))
prop_to_z_crossing_1_event = Event(
    name = "Prop to Z 1",
    event = prop_to_z_crossing_1_fun,
)

# Event 2: Apply TOI maneuver
toi_fun() = maneuver(sat, toi)
toi_event = Event(
    name = "TOI",
    event = toi_fun,
    vars = [var_toi_v],
)

# Event 3: Propagate to apoapsis and check radius constraint
prop_to_apogee_fun() = propagate(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=-1))
prop_to_apogee_event = Event(
    name = "Prop to Apoapsis",
    event = prop_to_apogee_fun,
    funcs = [apogee_radius_con],
)

# Event 4: Propagate to perigee
prop_to_perigee_1_fun() = propagate(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=1))
prop_to_perigee_1_event = Event(
    name = "Prop to Perigee 1",
    event = prop_to_perigee_1_fun,
)

# Event 5: Propagate to Z=0 crossing again
prop_to_z_crossing_2_fun() = propagate(prop, sat, StopAt(sat, PosZ(), 0.0))
prop_to_z_crossing_2_event = Event(
    name = "Prop to Z 2",
    event = prop_to_z_crossing_2_fun,
)

# Event 6: Apply MCC maneuver
mcc_fun() = maneuver(sat, mcc)
mcc_event = Event(
    name = "MCC",
    event = mcc_fun,
    vars = [var_mcc_vn],
)

# Event 7: Propagate to perigee and check constraints
prop_to_perigee_2_fun() = propagate(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=1))
prop_to_perigee_2_event = Event(
    name = "Prop to Perigee 2",
    event = prop_to_perigee_2_fun,
    funcs = [inclination_con, perigee_radius_con],
)

# Event 8: Apply MOI maneuver and check final SMA
moi_fun() = maneuver(sat, moi)
moi_event = Event(
    name = "MOI",
    event = moi_fun,
    vars = [var_moi_v],
    funcs = [final_sma_con],
)

# ========== Build the Sequence and Solve =================================================
seq = Sequence()
add_events!(seq, toi_event, [prop_to_z_crossing_1_event])
add_events!(seq, prop_to_apogee_event, [toi_event])
add_events!(seq, prop_to_perigee_1_event, [prop_to_apogee_event])
add_events!(seq, prop_to_z_crossing_2_event, [prop_to_perigee_1_event])
add_events!(seq, mcc_event, [prop_to_z_crossing_2_event])
add_events!(seq, prop_to_perigee_2_event, [mcc_event])
add_events!(seq, moi_event, [prop_to_perigee_2_event])

# Set up IPOPT options
ipopt_options = Dict(
    "max_iter" => 1000,
    "tol" => 1e-6,
    "output_file" => "ipopt_geo_transfer$(rand(UInt)).out",
    "file_print_level" => 5,
    "print_level" => 5,
)
snow_options = Options(derivatives=ForwardFD(), solver=IPOPT(ipopt_options))

result = trajectory_solve(seq, snow_options)
sequence_report(seq)
solution_report(seq, result)