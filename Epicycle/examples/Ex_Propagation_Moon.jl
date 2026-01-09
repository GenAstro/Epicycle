using Epicycle

# Spacecraft
moon_icrf = CoordinateSystem(moon, ICRFAxes())  
k =  KeplerianState(5500.0,0.38,deg2rad(80),deg2rad(145),deg2rad(180),deg2rad(0.01))
sat = Spacecraft(
    state = CartesianState(k,moon.mu),
    time = Time("2015-09-21T12:23:12", TAI(), ISOT()),
    coord_sys = moon_icrf,
    name = "MoonSat",
    cad_model = CADModel(
        file_path = joinpath(@__DIR__, "data", "DeepSpace1.obj"),
        scale = 100.0,
        visible = true
    )
)

# Forces + integrator
gravity = PointMassGravity(moon,(earth,sun))
forces  = ForceModel(gravity)
integ   = IntegratorConfig(Tsit5(); dt=10.0, reltol=1e-9, abstol=1e-9)
prop    = OrbitPropagator(forces, integ)

# Propagate to periapsis
propagate(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=+1))

# Create 3D view and display
view = View3D(coord_sys = moon_icrf)
add_spacecraft!(view, sat)
display_view(view)
