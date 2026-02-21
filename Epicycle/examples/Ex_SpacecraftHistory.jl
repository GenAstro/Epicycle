using Epicycle
using LinearAlgebra
    
# Spacecraft
sat = Spacecraft(
    state = CartesianState([3737.792, -4607.692, -2845.644, 5.411, 5.367, -1.566]),
    time = Time("2000-01-01T11:59:28.000", UTC(), ISOT()), 
    name = "GeoSat-1"
)

# Dynamics, Integator, Propagator
gravity = PointMassGravity(earth, ())  # Only Earth gravity
forces  = ForceModel(gravity)
integ   = IntegratorConfig(DP8(); abstol=1e-12, reltol=1e-12, dt=60.0)
prop    = OrbitPropagator(forces, integ)

# Maneuver model
toi = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 2.518,
    element2 = 0.0,
    element3 = 0.0,
)

#==============================================================================
# Multi-Phase Mission: Propagate → Maneuver → Propagate
# (History segments are created automatically)
==============================================================================#

propagate!(prop, sat, StopAt(sat, PropDurationDays(), 0.5))

maneuver!(sat, toi)

propagate!(prop, sat, StopAt(sat, PropDurationDays(), 0.5))

#==============================================================================
# History Container Basics
==============================================================================#

println("\n=== History Container Interface ===")

# Display entire history
println(sat.history)

# Number of segments
println("\nNumber of segments: ", length(sat.history))

# Check if empty
println("Is history empty? ", isempty(sat.history))

# Indexing
println("\nFirst segment: ", sat.history[1].name)
println("Last segment: ", sat.history[end].name)
println("Second segment: ", sat.history[2].name)

#==============================================================================
# Segment Metadata and Information
==============================================================================#

println("\n=== Segment Metadata ===")

# Access first segment
seg1 = sat.history[1]

println("Name: ", seg1.name)
println("Number of points: ", length(seg1.times))
println("Origin: ", seg1.coordinate_system.origin.name)
println("Axes: ", split(string(typeof(seg1.coordinate_system.axes)), ".")[end])

# Check if segment is empty
println("Is segment empty? ", isempty(seg1))

# Time span of segment
if !isempty(seg1)
    start_time = seg1.times[1]
    end_time = seg1.times[end]
    duration = end_time - start_time
    println("\nTime span:")
    println("  Start: ", start_time)
    println("  End: ", end_time)
    println("  Duration: ", duration, " days (", duration*24, " hours)")
end

#==============================================================================
# Extracting Data from History
==============================================================================#

println("\n=== Data Extraction ===")

# Access data from a single segment
seg1 = sat.history[1]
times_seg1 = seg1.times        # Vector{Time{Float64}}
states_seg1 = seg1.states      # Vector{CartesianState{Float64}}

println("First segment: $(length(times_seg1)) time-state pairs")
println("  First time: ", times_seg1[1])
println("  First position: ", states_seg1[1].position, " km")

# Extract data from all segments
all_positions = [state.position for seg in sat.history for state in seg.states]
all_times = [time for seg in sat.history for time in seg.times]

println("\nAll segments combined: $(length(all_positions)) total points")

# Convert to plotting-friendly formats
x = [p[1] for p in all_positions]
y = [p[2] for p in all_positions]
z = [p[3] for p in all_positions]
times_mjd = [t.mjd for t in all_times]

println("  Extracted X, Y, Z coordinates and MJD times")

# To visualize the data, uncomment and run:
# using Plots
# plot(x, y, z, label="Trajectory", xlabel="X (km)", ylabel="Y (km)", zlabel="Z (km)")
# plot(times_mjd, x, label="X position", xlabel="Time (MJD)", ylabel="X (km)")

