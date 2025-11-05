using Epicycle
using ForwardDiff
using FiniteDiff


# ==== FULL APPROACH: Include state and mass as input variables ====

# Function to simulate maneuver with all variables as inputs
# Input vector: [x, y, z, vx, vy, vz, mass, dv1, dv2, dv3]
function sim_maneuver_full(input_vec)
    # Unpack input vector
    x, y, z = input_vec[1], input_vec[2], input_vec[3]
    vx, vy, vz = input_vec[4], input_vec[5], input_vec[6]
    mass = input_vec[7]
    dv1, dv2, dv3 = input_vec[8], input_vec[9], input_vec[10]
    
    # Create spacecraft with Dual-typed initial state and mass
    sat = Spacecraft(
        state = CartesianState([x, y, z, vx, vy, vz]),
        mass = mass
    )
    
    # Create maneuver with Dual-typed delta-v components
    deltav = ImpulsiveManeuver(
        axes = Inertial(),
        g0 = 9.80665,
        Isp = 250.0,
        element1 = dv1,
        element2 = dv2,
        element3 = dv3
    )

    gravity = PointMassGravity(earth,())
    forces  = ForceModel(gravity)
    integ   = IntegratorConfig(Tsit5(); dt=10.0, reltol=1e-9, abstol=1e-9)
    prop    = OrbitPropagator(forces, integ)

    # Apply the maneuver to the spacecraft
    maneuver(sat, deltav)

    # Propagate to periapsis
    propagate(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=+1))
    
    # Return final position
    return sat.state.state[1:3]
end

# Set up nominal input values
# Default spacecraft: [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0] for state, 1000.0 for mass
nominal_input = [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0, 1000.0, 0.04, -0.3, 0.1]

# Test the full function
println("\n" * "="^50)
println("TESTING FULL APPROACH")
println("="^50)
println("Testing nominal function call...")
pos_final_full = sim_maneuver_full(nominal_input)
println("Final position: ", pos_final_full)

# Compute full Jacobian
println("\nComputing full Jacobian with ForwardDiff...")

J_full = ForwardDiff.jacobian(sim_maneuver_full, nominal_input)
println("Full Jacobian size: ", size(J_full))

# Extract only the ∂pos/∂dv portion (columns 8, 9, 10)
J_pos_dv = J_full[:, 8:10]
println("\n∂(position)/∂(delta-v) Jacobian:")
println(J_pos_dv)
println("Size: ", size(J_pos_dv))

println("\nPartial derivatives breakdown:")
println("∂x/∂(dv1, dv2, dv3) = ", J_pos_dv[1, :])
println("∂y/∂(dv1, dv2, dv3) = ", J_pos_dv[2, :])
println("∂z/∂(dv1, dv2, dv3) = ", J_pos_dv[3, :])
    
nothing # suppress output from last command
=#