using Epicycle

# Create spacecraft and define time and state
sat = Spacecraft(
    state=KeplerianState(8000.0,0.15,pi/4,pi/2,0.0,pi/2),
    time=Time("2015-09-21T12:23:12", TAI(), ISOT()),
    name = "sat",
)

# Define a propagator
gravity = PointMassGravity(earth,(moon,sun))
forces  = ForceModel(gravity)
integ   = IntegratorConfig(Tsit5(); dt=10.0, reltol=1e-9, abstol=1e-9)
prop    = OrbitPropagator(forces, integ)

# Propagate for 5000 seconds
propagate!(prop, sat, StopAt(sat, PropDurationSeconds(), 5000.0))
println(get_state(sat, Keplerian()))

# Visualize the orbit
view = View3D()
add_spacecraft!(view,sat)
display_view(view)