using Epicycle

# Spacecraft
sat = Spacecraft(
    state=CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]),
    time=Time("2015-09-21T12:23:12", TAI(), ISOT()),
    #name="SC-StopAt",
    coord_sys=CoordinateSystem(earth, ICRFAxes()),
)

# Forces + integrator
gravity = PointMassGravity(earth,(moon,sun))
forces  = ForceModel(gravity)
integ   = IntegratorConfig(Tsit5(); dt=10.0, reltol=1e-9, abstol=1e-9)
prop    = OrbitPropagator(forces, integ)

# Propagate to periapsis
propagate(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=+1))
println(get_state(sat, Keplerian()))

# Propagate to apoapsis
propagate(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=-1))
println(get_state(sat, Keplerian()))

# Stop when |r| reaches 7000 km 
propagate(prop, sat, StopAt(sat, PosMag(), 7000.0))
println(get_state(sat, SphericalRADEC()))       

# Propagate to x-position crossing (increasing)
sol = propagate(prop, sat, StopAt(sat, PosX(), 7.5; direction=+1))
println(get_state(sat, Cartesian()))

nothing
