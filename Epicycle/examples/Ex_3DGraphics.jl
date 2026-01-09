using Epicycle

# Create a spacecraft and use a CAD model to visualize it
sat = Spacecraft(
    state = CartesianState([3737.792, -4607.692, -2845.644, 5.411, 5.367, -1.566]),
    name = "Deep Space 1", 
    cad_model = CADModel(
        file_path = joinpath(pkgdir(Epicycle), "assets", "DeepSpace1.obj"),
        scale = 100.0,
        visible = true
    )
)

# Propagate the spacecraft to ascending node
gravity = PointMassGravity(earth, (sun,moon)) 
forces  = ForceModel(gravity)
integ   = IntegratorConfig(DP8(); abstol=1e-12, reltol=1e-12, dt=60.0)
prop    = OrbitPropagator(forces, integ)
propagate(prop, sat, StopAt(sat, PosZ(), 0.0; direction = +1))

# Create a view, add the spacecraft, and display it
view = View3D()
add_spacecraft!(view, sat)
display_view(view)
