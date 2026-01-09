using Epicycle

# Spacecraft
mars_icrf = CoordinateSystem(mars, ICRFAxes())  
k =  KeplerianState(6500.0,0.3,deg2rad(30),deg2rad(145),deg2rad(180),deg2rad(0.01))
sat = Spacecraft(
    state = CartesianState(k,mars.mu),
    time = Time("2015-09-21T12:23:12", TAI(), ISOT()),
    coord_sys = mars_icrf,
    name = "MarsSat"
)

# Forces + integrator
gravity = PointMassGravity(mars, ())
forces  = ForceModel(gravity)
integ   = IntegratorConfig(Tsit5(); dt=10.0, reltol=1e-9, abstol=1e-9)
prop    = OrbitPropagator(forces, integ)

# Propagate to periapsis
propagate(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=+1))

# Create 3D view and display
view = View3D(coord_sys = mars_icrf)
add_spacecraft!(view, sat)
display_view(view)
