using Epicycle

# Spacecraft
sat = Spacecraft(
    state = CartesianState([5000.0, 5000.0, 0.0, -3.8, 3.8, 5.4]),
    time = Time("2015-09-21T12:23:12", TAI(), ISOT()),
    coord_sys = CoordinateSystem(earth, ICRFAxes()),
    name = "Sat",
)

# Forces + integrator
gravity = PointMassGravity(earth,(moon,sun))
forces  = ForceModel(gravity)
integ   = IntegratorConfig(Tsit5(); dt=10.0, reltol=1e-9, abstol=1e-9)
prop    = OrbitPropagator(forces, integ)

# Propagate for 1 hour
propagate(prop, sat, StopAt(sat, PropDurationSeconds(), 3600.0))

# Propagate for 0.5 days
propagate(prop, sat, StopAt(sat, PropDurationDays(), 0.5))

# Propagate to an epoch
stop_time = Time("2015-09-22T12:00:00", TDB(), ISOT())
propagate(prop, sat, StopAt(sat, stop_time))

# Propagate to periapsis 
propagate(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=+1))
println(get_state(sat, Keplerian()))

# Propagate to ascending node crossing (increasing)
# Note: in StopAt(), direction=+1 means locate root when function is increasing
sol = propagate(prop, sat, StopAt(sat, PosZ(), 0.0; direction=+1))
println(get_state(sat, Cartesian()))

# Propagate to apoapsis
# Note: in StopAt(), direction=-1 means locate root when function is decreasing
propagate(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=-1))
println(get_state(sat, Keplerian()))

# Propagate to |r| = 7250 km 
propagate(prop, sat, StopAt(sat, PosMag(), 7250.0))
println(get_state(sat, SphericalRADEC()))       

# Propagate backwards for 2 hours
propagate(prop, sat, StopAt(sat, PropDurationSeconds(), -7200.0); direction=:infer)   

# Propagate backwards to an epoch
stop_time_back = sat.time - 0.05
propagate(prop, sat, StopAt(sat, stop_time_back); direction=:infer)   

# Plot the trajectory 3D
view = View3D()
add_spacecraft!(view, sat)
display_view(view)

nothing
