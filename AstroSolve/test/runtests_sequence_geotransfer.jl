


using Test
using SNOW
using OrdinaryDiffEq
using LinearAlgebra

@testset "GEO Transfer Sequence vs GMAT" begin

# Create spacecraft matching GMAT initial state
sat = Spacecraft(
    state = CartesianState([3737.792, -4607.692, -2845.644, 5.411, 5.367, -1.566]),
    time = Time("2000-01-01T11:59:28.000", UTC(), ISOT())  # J2000 epoch from GMAT
)

# Create simple Earth point mass dynamics (no third bodies for now)
gravity = PointMassGravity(earth,())  # Only Earth gravity
forces  = ForceModel(gravity)
integ   = IntegratorConfig(DP8(); abstol=1e-12, reltol=1e-12, dt=60.0)
prop    = OrbitPropagator(forces, integ)

# Define maneuver models matching GMAT
toi = ImpulsiveManeuver(
    axes = VNB(),        # VNB coordinates like GMAT
    element1 = 2.818,      # V component (to be varied)
    element2 = 0.0,      # N component  
    element3 = 0.0,      # B component
)

mcc = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 0.759,      # V component (to be varied)
    element2 = 0.788,      # N component (to be varied)
    element3 = 0.0,      # B component
)

moi = ImpulsiveManeuver(
    axes = VNB(), 
    element1 = -0.482,      # V component (to be varied)
    element2 = 0.0,      # N component
    element3 = 0.0,      # B component
)

# Define solver variables matching GMAT targeting
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

# Define constraints matching GMAT targets

# Target 1: Apoapsis radius = 85,000 km after TOI
apogee_radius_con = Constraint(
    calc = OrbitCalc(sat, PosMag()),
    lower_bounds = [85000.0],
    upper_bounds = [85000.0],
    scale = [1.0],
)

# Target 2a: Inclination = 2° after MCC
inclination_con = Constraint(
    calc = OrbitCalc(sat, Inc()),   
    lower_bounds = [deg2rad(2.0)],
    upper_bounds = [deg2rad(2.0)], 
    scale = [1.0],
)

# Target 2b: Perigee radius = 42,195 km after MCC
perigee_radius_con = Constraint(
    calc = OrbitCalc(sat, PosMag()),
    lower_bounds = [42195.0],
    upper_bounds = [42195.0],
    scale = [1.0],
)

# Target 3: Final SMA = 42,166.90 km (GEO) after MOI
final_sma_con = Constraint(
    calc = OrbitCalc(sat, SMA()),
    lower_bounds = [42166.90],
    upper_bounds = [42166.90],
    scale = [1.0],
)

# Event 1: Propagate to Z=0 crossing (equatorial plane)
prop_to_z_crossing_1_fun() = propagate(prop, sat, StopAt(sat, PosZ(), 0.0))
prop_to_z_crossing_1_event = Event(
    name = "prop_to_z_crossing_1",
    event = prop_to_z_crossing_1_fun,
)

# Event 2: Apply TOI maneuver 
toi_fun() = maneuver(sat, toi)
toi_event = Event(
    name = "toi",
    event = toi_fun,
    vars = [var_toi_v],
)

# Event 3: Propagate to apoapsis and check radius constraint
prop_to_apogee_fun() = propagate(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=-1))
prop_to_apogee_event = Event(
    name = "prop_to_apogee", 
    event = prop_to_apogee_fun,
    funcs = [apogee_radius_con],
)

# Event 4: Propagate to perigee
prop_to_perigee_1_fun() = propagate(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=1))
prop_to_perigee_1_event = Event(
    name = "prop_to_perigee_1",
    event = prop_to_perigee_1_fun,
)

# Event 5: Propagate to Z=0 crossing again
prop_to_z_crossing_2_fun() = propagate(prop, sat, StopAt(sat, PosZ(), 0.0))
prop_to_z_crossing_2_event = Event(
    name = "prop_to_z_crossing_2", 
    event = prop_to_z_crossing_2_fun,
)

# Event 6: Apply MCC maneuver
mcc_fun() = maneuver(sat, mcc)
mcc_event = Event(
    name = "mcc",
    event = mcc_fun,
    vars = [var_mcc_vn],
)

# Event 7: Propagate to perigee and check constraints
prop_to_perigee_2_fun() = propagate(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=1))
prop_to_perigee_2_event = Event(
    name = "prop_to_perigee_2",
    event = prop_to_perigee_2_fun,
    funcs = [inclination_con, perigee_radius_con],
)

# Event 8: Apply MOI maneuver and check final SMA
moi_fun() = maneuver(sat, moi)
moi_event = Event(
    name = "moi",
    event = moi_fun,
    vars = [var_moi_v],
    funcs = [final_sma_con],
)

# Build the complete sequence with dependencies
seq = Sequence()
add_events!(seq, toi_event, [prop_to_z_crossing_1_event])
add_events!(seq, prop_to_apogee_event, [toi_event])
add_events!(seq, prop_to_perigee_1_event, [prop_to_apogee_event])
add_events!(seq, prop_to_z_crossing_2_event, [prop_to_perigee_1_event])
add_events!(seq, mcc_event, [prop_to_z_crossing_2_event])
add_events!(seq, prop_to_perigee_2_event, [mcc_event])
add_events!(seq, moi_event, [prop_to_perigee_2_event])

# Create sequence manager and execute
sm = SequenceManager(seq)

# Get initial variable values and test function evaluation
x0 = get_var_values(sm)
F = zeros(length(sm.ordered_funcs))

# Execute sequence with GMAT solution values
status = solver_fun!(F, x0, sm)

# GMAT truth data for regression testing
F_GMAT = [84593.426030930; deg2rad(2.0374382158683); 42006.900021350; 41969.105396615]

# Test: Compare our results to GMAT with relative tolerance
rel_tol = 1e-7
@test isapprox(F, F_GMAT; rtol=rel_tol)

end  # testset
nothing


    