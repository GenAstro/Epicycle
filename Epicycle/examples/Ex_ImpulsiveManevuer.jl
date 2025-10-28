using Epicycle

# Create a default spacecraft
sat1 = Spacecraft()

# Create an impulsive maneuver in the Inertial frame
deltav2 = ImpulsiveManeuver(
      axes = Inertial(),
      g0 = 9.80665,
      Isp = 250.0,
      element1 = 0.04,
      element2 = -0.3,
      element3 = 0.1
     )

# Apply the maneuver to the spacecraft
println("Initial mass: ", sat1.mass)
maneuver(sat1, deltav2)
println("Mass after Inertial maneuver: ", sat1.mass)
println("State after Inertial maneuver: \n", get_state(sat1, Cartesian()))

# Create an impulsive maneuver in the VNB frame
deltav1 = ImpulsiveManeuver(
      axes = VNB(),
      g0 = 9.80665,
      Isp = 250.0,
      element1 = 0.2,
      element2 = 0.1,
      element3 = -0.2
      )

# Apply the maneuver to the spacecraft
maneuver(sat1, deltav1)
println("Mass after VNB maneuver: ", sat1.mass)
println("State after VNB maneuver: \n", get_state(sat1, Cartesian()))

nothing # suppress output from last command