
using SNOW
using DifferentialEquations
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
using ForwardDiff

# Create spacecraft
sat1 = Spacecraft(
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
    obj = toi, 
    variable = DeltaV(), 
    name = "toi",
    lower_bound = [-10.0, 0.0, 0.0],
    upper_bound = [10.0, 0.0, 0.0],
)

# Define moi as a solver variable
var_moi = SolverVariable(
    obj = moi, 
    variable = DeltaV(), 
    name = "moi",
    lower_bound = [-10.0, 0.0, 0.0],
    upper_bound = [10.0, 0.0, 0.0]
)

pos_target = 45000
pos_con_fun() = pos_mag(sat1)
pos_con = Constraint(
    func = pos_con_fun,
    lower_bounds = [pos_target],
    upper_bounds = [pos_target],
    scale = [1.0],
    numvars=1
)

ecc_con_fun() = eccentricity(sat1, earth.mu)
ecc_con = Constraint(
    func = ecc_con_fun,
    lower_bounds = [0],
    upper_bounds = [0], 
    scale = [1.0],
    numvars = 1
)

vel_target = sqrt(earth.mu / pos_target)
vel_con_fun() = vel_mag(sat1)
vel_con = Constraint(
    func = vel_con_fun,
    lower_bounds = [vel_target],
    upper_bounds = [vel_target],
    scale = [1.0],
    numvars = 1
)

# Create the MOI Event
toi_fun() = maneuver(sat1, toi) 
toi_event = Event(name = "toi", 
                  event = toi_fun, 
                  vars = [var_toi],
                  funcs = [])

# Create the prop to apopasis event
prop_apo_fun() = propagate(dynsys, integ, StopAtApoapsis(sat1))
prop_event = Event(name = "prop_apo", event = prop_apo_fun)

# Create the TOI event. 
moi_fun() = maneuver(sat1, moi)

moi_event = Event(event = moi_fun, 
                  vars = [var_moi],
                  funcs = [pos_con, ecc_con])

# Build sequence and solve
seq = Sequence()
add_events!(seq, prop_event, [toi_event]) 
add_events!(seq, moi_event, [prop_event])

sm = SequenceManager(seq)

ordered_vars = sm.ordered_vars

#function set_all_vars_to_dual(ordered_vars::Vector{SolverVariable})
for v in ordered_vars
    val = get_sol_var(v.obj, v.variable)
    dual_val = ForwardDiff.Dual.(val, ones(length(val)))
    set_sol_var(v.obj, v.variable, dual_val)
    # Check and print types for debugging
    println("Variable: ", v.name)
    println("Type of obj: ", typeof(v.obj))
    println("Type of field after set: ", typeof(get_sol_var(v.obj, v.variable)))
    println("Field values: ", get_sol_var(v.obj, v.variable))
end

#end
#=
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
#options = Options(derivatives=ForwardFD(), ipopt_options=Dict("print_level"=>0))

# Define a closure that matches SNOW's expected signature
snow_solver_fun!(F, x) = solver_fun!(F, x, sm)

xopt, fopt, info = minimize(snow_solver_fun!, x0, ng, lx, ux, lg, ug, options)

=#
