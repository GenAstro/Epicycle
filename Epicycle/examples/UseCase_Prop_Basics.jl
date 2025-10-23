# Create spacecraft
using Epicycle

sat = Spacecraft(
    state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]), 
    time = Time("2015-09-21T12:23:12", TAI(), ISOT())
    )

sat = Spacecraft(
    state = OrbitState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03],Cartesian()), 
    time = Time("2015-09-21T12:23:12", TAI(), ISOT())
    )

# Create force models, integrator, and dynamics system
gravity = PointMassGravity(earth,())
forces = ForceModel(gravity)
integ = IntegratorConfig(DP8(); abstol = 1e-11, reltol = 1e-11, dt = 200)

# Define which spacecraft to propagate and which force model to use
dynsys = DynSys(
     forces = forces, 
     spacecraft = [sat]
     )

# Propagate for 1 hour in seconds
propagate(dynsys, integ, StopAtSeconds(3600.0))
println("prop 10 days: ", sat.time)

# Propagate for 10 days
propagate(dynsys, integ, StopAtDays(10.0))
println("prop 10 days: ", sat.time)

# Propagate to apoapsis
propagate(dynsys, integ, StopAtApoapsis(sat))
println("prop to apoapsis: ", get_state(sat, ModifiedKeplerian()))

# Propagte to periapsis, print Modified Keplerian state
propagate(dynsys, integ, StopAtPeriapsis(sat))
println("prop to periapsis: ", get_state(sat, Keplerian()))

# Propagate to ascending node
propagate(dynsys, integ, StopAtAscendingNode(sat))
println("prop to node crossing: ", get_state(sat, Cartesian()))

# Propagate to a radial distance from central body
propagate(dynsys, integ, StopAtRadius(sat,7100.0))
println("prop to radius = 7100: ", get_state(sat, SphericalRADEC()))

# Use an OrbitCalc in the stopping condition
#propagate(dynsys, integ, StopAt(OrbitCalc(sat, PosMag()),7100.0))
#println("prop to radius = 7100: ", get_state(sat, SphericalRADEC()))