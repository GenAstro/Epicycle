


using Epicycle

using SNOW
using OrdinaryDiffEq
using LinearAlgebra
using ForwardDiff

# Create Dual types for the initial values with proper tag ordering
struct EpicycleADTag end

# Define tag ordering to resolve conflicts with OrdinaryDiffEq
import ForwardDiff: ≺
≺(::Type{EpicycleADTag}, ::Type{T}) where {T<:ForwardDiff.Tag} = true
≺(::Type{T}, ::Type{EpicycleADTag}) where {T<:ForwardDiff.Tag} = false
≺(::Type{EpicycleADTag}, ::Type{ForwardDiff.Tag{F,V}}) where {F,V} = true

SolverDataType = ForwardDiff.Dual{EpicycleADTag, Float64, 3}
#SolverDataType = Float64

sat = Spacecraft(
    state = CartesianState(SolverDataType.([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0])), 
    time = Time(SolverDataType(2.459114e6), SolverDataType(0.5), TAI(), JD()),
    mass = SolverDataType(1000.0),
)

# Create an impulsive maneuver with Dual types
toiman = ImpulsiveManeuver(
    axes = VNB(),
    element1 = SolverDataType(0.1),
    element2 = SolverDataType(0.2),
    element3 = SolverDataType(0.3),
    Isp = SolverDataType(300.0),
    g0 = SolverDataType(9.81),
)

# Create the propagator
gravity = PointMassGravity(earth, ())
forces  = ForceModel(gravity)
integ   = IntegratorConfig(Tsit5(); dt=10.0, reltol=1e-9, abstol=1e-9)
prop    = OrbitPropagator(forces, integ)

# Promote spacecraft and maneuver to Dual numbers for AD
# ForwardDiff extension will handle the Dual type conversion automatically
# sat = promote(sat, ForwardDiff.Dual)
#toi = promote(toi, ForwardDiff.Dual)

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

F = SolverDataType.([0.0])
x = SolverDataType.([0.2, -1.0, 0.4])

status = solver_fun!(F, x, sm)
println("Status: ", status)
println("F values: ", F)

# New function

function solver_jac!(F, x, sm)
    J = zeros(length(F), length(x))
    ForwardDiff.jacobian!(J, (y) -> solver_fun!(similar(F), y, sm), x)
    return J
end

# Test the Jacobian computation
x_test = [0.2, -1.0, 0.4]  # Float64 inputs for ForwardDiff
F_test = [0.0]  # Float64 output vector
J = solver_jac!(F_test, x_test, sm)
println("Jacobian: ")
println(J)
println("Jacobian size: ", size(J))