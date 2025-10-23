


using SNOW
using DifferentialEquations
using LinearAlgebra
using ForwardDiff

using AstroEpochs
using AstroStates
using AstroBase
using AstroUniverse
using AstroCoords
using AstroProp
using AstroMan
using AstroSolve
using AstroFun

# Promote numeric leaves of a struct via its Functors.functor to eltype T
_promote_via_functor(x, ::Type{T}) where {T} = begin
    ch, re = Functors.functor(x)
    chT = Functors.fmap(ch) do v
        if v isa AbstractArray{<:Real}
            T.(v)                 # broadcast preserves container shape (StaticArrays too)
        elseif v isa Real
            convert(T, v)
        else
            v                     # leave non-numeric leaves (e.g., axes, symbols) untouched
        end
    end
    re(chT)
end

function _promote_spacecraft(sc::Spacecraft, ::Type{T}) where {T}
    st_T   = _promote_via_functor(sc.state, T)                 # positions/velocities -> T
    mass_T = sc.mass isa Real ? convert(T, sc.mass) : sc.mass  # optional: promote mass
    hist_F64 = Vector{Vector{Tuple{Time, Vector{Float64}}}}()
    return Spacecraft(
        state     = st_T,
        time      = sc.time,        # keep Time intact
        mass      = mass_T,
        name      = sc.name,
        history   = hist_F64,
        coord_sys = sc.coord_sys,
    )
end

# Prepare an AD-ready model by promoting Spacecraft and Maneuvers
function prepare_ad_model(sat::Spacecraft, mans::AbstractVector{<:ImpulsiveManeuver}; T::Type=Float64)
    sat_T  = _promote_spacecraft(sat, T)
    mans_T = [_promote_via_functor(m, T) for m in mans]
    return (sat = sat_T, maneuvers = mans_T)
end

# Promote model to ForwardDiff duals (independent of problem size)
function promote_model_to_dual(sat::Spacecraft, mans::AbstractVector{<:ImpulsiveManeuver})
    Tdual = ForwardDiff.Dual{Nothing, Float64, 1}
    return prepare_ad_model(sat, mans; T=Tdual)
end

# Optional convenience for a single maneuver
function promote_model_to_dual(sat::Spacecraft, man::ImpulsiveManeuver)
    prom = promote_model_to_dual(sat, [man])
    return (sat = prom.sat, maneuvers = prom.maneuvers[1:1])
end

# Create spacecraft
sat = Spacecraft(
            state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0]), 
            time = Time("2020-09-21T12:23:12", TAI(), ISOT())
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

# Promote model here (choose T as needed; Float64 by default, AD will pass a Dual type)
prom = promote_model_to_dual(sat, [toi, moi])
sat = prom.sat; toi, moi = prom.maneuvers 

# Create force models, integrator, and dynamics system
pm_grav = PointMassGravity(earth,(moon,sun))
forces = ForceModel(pm_grav)
integ = IntegratorConfig(DP8(); abstol = 1e-11, reltol = 1e-11, dt = 4000)

# Define which spacecraft to propagate and which force model to use
dynsys = DynSys(
          forces = forces, 
          spacecraft = [sat]
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
    scale = [1.0]
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
#= 
f = get_fun_values(sm)

x0 = get_var_values(sm)

set_var_values(sm, x0)
=#

nothing