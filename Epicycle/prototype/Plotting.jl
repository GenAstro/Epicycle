using GLMakie
using FileIO  # For loading images
using Random  # For reproducible star field
using MeshIO  # For loading 3D model files
using GeometryBasics  # For mesh handling
using LinearAlgebra  # For cross product and normalize

# Screen management - reuse windows instead of creating new ones
const SCREEN_3D = Ref{Union{Nothing,GLMakie.Screen}}(nothing)
const SCREEN_2D = Ref{Union{Nothing,GLMakie.Screen}}(nothing)

# Load texture image
earth_img = load(raw"C:\Users\steve\Dev\Epicycle\Epicycle\prototype\earthTexture.jpg")
milkyway_img_raw = load(raw"C:\Users\steve\Dev\Epicycle\Epicycle\prototype\MilkyWay.jpg")

# Filter out dim stars - keep only bright features (galactic cloud, bright stars, galaxies)
# Convert to brightness and threshold
brightness_threshold = 0.2 # Adjust this: higher = fewer stars, lower = more stars
milkyway_filtered = map(milkyway_img_raw) do pixel
    # Calculate brightness (luminance)
    brightness = 0.299 * pixel.r + 0.587 * pixel.g + 0.114 * pixel.b
    
    if brightness < brightness_threshold
        # Replace dim pixels with black
        typeof(pixel)(0, 0, 0)
    else
        # Keep bright pixels, slightly dimmed
        pixel * 0.5
    end
end

milkyway_img = milkyway_filtered 

"""
    extract_trajectory(sc::Spacecraft) -> (Vector{Float64}, Vector{Float64}, Vector{Float64})

Extract position trajectory from spacecraft history.
Returns separate x, y, z coordinate vectors in km.

# Arguments
- sc: Spacecraft with populated history field

# Returns
- Tuple of (x, y, z) vectors representing the trajectory

# Notes
- History structure: Vector{Vector{Tuple{Time{Float64}, Vector{Float64}}}}
- Each segment contains time-state pairs: (time, [x,y,z,vx,vy,vz])
- All segments are concatenated into a single trajectory
"""
function extract_trajectory(sc)
    x_traj = Float64[]
    y_traj = Float64[]
    z_traj = Float64[]
    
    # Iterate through all history segments
    for segment in sc.history
        for (time, posvel) in segment
            push!(x_traj, posvel[1])
            push!(y_traj, posvel[2])
            push!(z_traj, posvel[3])
        end
    end
    
    return x_traj, y_traj, z_traj
end

"""
    extract_trajectory_segments(sc::Spacecraft) -> Vector{Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}}

Extract position trajectory from spacecraft history, keeping segments separate.
Returns a vector of (x, y, z) tuples, one per history segment.

# Arguments
- sc: Spacecraft with populated history field

# Returns
- Vector of tuples, each containing (x_vec, y_vec, z_vec) for one segment

# Notes
- Useful for plotting each segment in a different color
- Each segment represents a continuous propagation arc
"""
function extract_trajectory_segments(sc)
    segments = []
    
    # Iterate through all history segments
    for segment in sc.history
        x_seg = Float64[]
        y_seg = Float64[]
        z_seg = Float64[]
        
        for (time, posvel) in segment
            push!(x_seg, posvel[1])
            push!(y_seg, posvel[2])
            push!(z_seg, posvel[3])
        end
        
        push!(segments, (x_seg, y_seg, z_seg))
    end
    
    return segments
end

segments = extract_trajectory_segments(sat)

# Define colors for different segments
segment_colors = [:cyan, :magenta, :yellow, :green, :orange, :red, :lightblue]

# Create figure with space-like theme
fig = Figure(size=(800, 600), backgroundcolor=:black)
lscene = LScene(fig[1, 1], 
    show_axis=false,
    scenekw=(
        backgroundcolor=:black,
        clear=true
    ))

# Fix rotation center - disable zoom from shifting the lookat point
cc = cameracontrols(lscene.scene)
cc.settings.zoom_shift_lookat = false

# Plot each segment in a different color (do this first to set camera scale)
for (i, (x, y, z)) in enumerate(segments)
    color = segment_colors[mod1(i, length(segment_colors))]
    lines!(lscene, x, y, z, 
        color=color, 
        linewidth=1.5)
end

# Add textured Earth sphere
n_lon = 100
n_lat = 100
θ_sphere = range(0, 2π, length=n_lon)  # Longitude
φ_sphere = range(0, π, length=n_lat)   # Latitude

R_earth = 6378.0  # km

# Create meshgrid for sphere
x_earth = [R_earth * sin(φ) * cos(θ) for φ in φ_sphere, θ in θ_sphere]
y_earth = [R_earth * sin(φ) * sin(θ) for φ in φ_sphere, θ in θ_sphere]
z_earth = [R_earth * cos(φ) for φ in φ_sphere, θ in θ_sphere]

surface!(lscene, x_earth, y_earth, z_earth, 
    color=earth_img,  # Use texture image
    shading=FastShading)

# Add dome over Earth location
dome_lat = 30.0  # degrees South
dome_lon = 80.0  # degrees East
dome_radius = 2500.0  # km - dome radius
dome_height = R_earth + 100.0  # km - altitude of dome base above Earth center

# Convert lat/lon to radians
lat_rad = deg2rad(dome_lat)
lon_rad = deg2rad(dome_lon)

# Calculate dome center position on Earth surface
dome_center_x = dome_height * cos(lat_rad) * cos(lon_rad)
dome_center_y = dome_height * cos(lat_rad) * sin(lon_rad)
dome_center_z = dome_height * sin(lat_rad)
dome_center = Vec3f(dome_center_x, dome_center_y, dome_center_z)

# Create local coordinate system for dome
# Radial direction (outward from Earth)
radial = normalize(dome_center)

# Find tangent vectors (local east and north)
north_pole = Vec3f(0, 0, 1)
if abs(dot(radial, north_pole)) > 0.99
    east = normalize(cross(radial, Vec3f(1, 0, 0)))
else
    east = normalize(cross(radial, north_pole))
end
north = normalize(cross(east, radial))

# Generate dome mesh (hemisphere)
n_dome_lat = 30
n_dome_lon = 30
φ_dome = range(0, π/2, length=n_dome_lat)  # 0 to 90 degrees (hemisphere)
θ_dome = range(0, 2π, length=n_dome_lon)

# Create dome vertices in local coordinates, then transform to world
dome_vertices = Point3f[]
for φ in φ_dome
    for θ in θ_dome
        # Spherical dome point in local frame
        local_x = dome_radius * sin(φ) * cos(θ)
        local_y = dome_radius * sin(φ) * sin(θ)
        local_z = dome_radius * cos(φ)
        
        # Transform to world coordinates
        world_pos = dome_center + local_x * east + local_y * north + local_z * radial
        push!(dome_vertices, Point3f(world_pos))
    end
end

# Create dome mesh faces
dome_faces = GeometryBasics.TriangleFace[]
for i in 1:(n_dome_lat-1)
    for j in 1:(n_dome_lon-1)
        idx1 = (i-1) * n_dome_lon + j
        idx2 = (i-1) * n_dome_lon + j + 1
        idx3 = i * n_dome_lon + j
        idx4 = i * n_dome_lon + j + 1
        
        push!(dome_faces, GeometryBasics.TriangleFace(idx1, idx2, idx3))
        push!(dome_faces, GeometryBasics.TriangleFace(idx2, idx4, idx3))
    end
    # Wrap around
    idx1 = (i-1) * n_dome_lon + n_dome_lon
    idx2 = (i-1) * n_dome_lon + 1
    idx3 = i * n_dome_lon + n_dome_lon
    idx4 = i * n_dome_lon + 1
    
    push!(dome_faces, GeometryBasics.TriangleFace(idx1, idx2, idx3))
    push!(dome_faces, GeometryBasics.TriangleFace(idx2, idx4, idx3))
end

dome_mesh = GeometryBasics.Mesh(dome_vertices, dome_faces)
mesh!(lscene, dome_mesh,
    color=RGBAf(1.0, 1.0, 0.0, 0.4),  # More vibrant yellow with transparency
    transparency=true,
    shading=FastShading)

# Add circle at the base of the dome
base_circle_points = Point3f[]
θ_base = range(0, 2π, length=n_dome_lon+1)
for θ in θ_base
    # Base is at φ = π/2 (equator of dome hemisphere)
    local_x = dome_radius * cos(θ)
    local_y = dome_radius * sin(θ)
    local_z = 0.0
    
    # Transform to world coordinates
    world_pos = dome_center + local_x * east + local_y * north + local_z * radial
    push!(base_circle_points, Point3f(world_pos))
end

lines!(lscene, base_circle_points,
    color=RGBAf(1.0, 1.0, 0.0, 0.9),  # Bright yellow, mostly opaque
    linewidth=1.0)

# OPTION: Star field
star_radius = 10_000_000.0  # km - very large radius
n_stars = 5000
Random.seed!(42)  # For reproducible star positions
θ_stars = rand(n_stars) .* 2π
φ_stars = acos.(2 .* rand(n_stars) .- 1)  # Uniform distribution on sphere
x_stars = star_radius .* sin.(φ_stars) .* cos.(θ_stars)
y_stars = star_radius .* sin.(φ_stars) .* sin.(θ_stars)
z_stars = star_radius .* cos.(φ_stars)
scatter!(lscene, x_stars, y_stars, z_stars, 
    color=:white, 
    markersize=1.5,
    glowwidth=0.3,
    glowcolor=:white)

# Load and display spacecraft 3D model
# Supported formats: .obj, .stl, .ply, .off
try
    spacecraft_model = FileIO.load(raw"C:\Users\steve\Dev\Epicycle\Epicycle\prototype\DeepSpace1.obj")
    
    # Get a position along the first trajectory segment (e.g., midpoint)
    seg_idx = 10
    pos_idx = length(segments[seg_idx][1]) ÷ 2
    
    x_pos = segments[seg_idx][1][end]
    y_pos = segments[seg_idx][2][end]
    z_pos = segments[seg_idx][3][end]
    
    # Display the spacecraft model at the trajectory position
    # Use Makie.qrotation for quaternion rotation (axis, angle)
    mesh!(lscene, spacecraft_model,
        color=:yellow,
        transformation = (
            translation = Vec3f(x_pos, y_pos, z_pos),
            rotation = qrotation(Vec3f(1, 0, 0), 3*π/2),  # 90° around x-axis
            scale = Vec3f(750.0)
        ))
    
    # Add body x
    vec_length = 5000.0  # km - visual length of arrow
    vec_dir = Vec3f(0, 0, 1)  # Placeholder direction (along x-axis)
    
    arrows!(lscene, 
        [Point3f(x_pos, y_pos, z_pos)],  # Start point at spacecraft
        [vec_dir * vec_length],           # Direction and length
        color=:cyan,
        linewidth=200,           # Thicker shaft line
        arrowsize=500,        # Smaller cone relative to length
        lengthscale=1.0)

    vec_dir = Vec3f(0, 1, 0)  # Placeholder direction (along x-axis)
    
    arrows!(lscene, 
        [Point3f(x_pos, y_pos, z_pos)],  # Start point at spacecraft
        [vec_dir * vec_length],           # Direction and length
        color=:cyan,
        linewidth=200,           # Thicker shaft line
        arrowsize=500,        # Smaller cone relative to length
        lengthscale=1.0)

    vec_dir = Vec3f(-1, 0, 0)  # Placeholder direction (along x-axis)
    
    arrows!(lscene, 
        [Point3f(x_pos, y_pos, z_pos)],  # Start point at spacecraft
        [vec_dir * vec_length],           # Direction and length
        color=:cyan,
        linewidth=200,           # Thicker shaft line
        arrowsize=500,        # Smaller cone relative to length
        lengthscale=1.0)
    
    # Add sensor cone
    cone_height = 6500.0  # km - how far the cone extends
    cone_angle = deg2rad(20)  # Half-angle of the cone in degrees
    cone_direction = normalize(Vec3f(0, 0, -1))  # Direction cone points (adjust to match sensor pointing)
    
    # Create cone geometry
    n_segments = 32  # Number of segments around the cone
    cone_radius = cone_height * tan(cone_angle)
    
    # Generate cone vertices
    θ_cone = range(0, 2π, length=n_segments+1)
    
    # Cone tip at spacecraft
    tip = Point3f(x_pos, y_pos, z_pos)
    
    # Create two perpendicular vectors to cone_direction for the base circle
    # Find a vector not parallel to cone_direction
    if abs(cone_direction[3]) < 0.9
        perp1 = normalize(cross(cone_direction, Vec3f(0, 0, 1)))
    else
        perp1 = normalize(cross(cone_direction, Vec3f(0, 1, 0)))
    end
    perp2 = normalize(cross(cone_direction, perp1))
    
    # Calculate base center point
    base_center = Vec3f(x_pos, y_pos, z_pos) + cone_direction * cone_height
    
    # Base circle points using perpendicular vectors
    base_points = [Point3f(
        base_center + cone_radius * (cos(θ) * perp1 + sin(θ) * perp2)
    ) for θ in θ_cone]
    
    # # Draw cone edges (lines from tip to base circle)
    # for i in 1:8  # Draw 8 edges for visibility
    #     idx = 1 + (i-1) * (n_segments ÷ 8)
    #     lines!(lscene, 
    #         [tip, base_points[idx]],
    #         color=RGBAf(1, 1, 0, 0.4),  # Yellow with transparency
    #         linewidth=1.5)
    # end
    
    # Draw base circle
    lines!(lscene, base_points,
        color=RGBAf(1, 1, 0, 0.5),  # Yellow with transparency
        linewidth=2.0)
    
    # Fill cone with semi-transparent mesh
    # Create triangular mesh for cone surface
    cone_vertices = vcat([tip], base_points[1:end-1])
    cone_faces = [GeometryBasics.TriangleFace(1, i+1, i+2) for i in 1:n_segments-1]
    push!(cone_faces, GeometryBasics.TriangleFace(1, n_segments, 2))  # Close the cone
    
    cone_mesh = GeometryBasics.Mesh(cone_vertices, cone_faces)
    mesh!(lscene, cone_mesh,
        color=RGBAf(1, 1, 0, 0.15),  # Very transparent yellow fill
        transparency=true)
    
    println("Spacecraft model loaded at position ($x_pos, $y_pos, $z_pos) km")
catch e
    println("Could not load spacecraft model: ", e)
end

# OPTION: Milky Way skybox (commented out in favor of star field)
# skybox_radius = 10_000_000.0  # km - very large radius
# n_lon_sky = 100
# n_lat_sky = 100
# θ_sky = range(0, 2π, length=n_lon_sky)
# φ_sky = range(0, π, length=n_lat_sky)
# x_sky = [-skybox_radius * sin(φ) * cos(θ) for φ in φ_sky, θ in θ_sky]
# y_sky = [-skybox_radius * sin(φ) * sin(θ) for φ in φ_sky, θ in θ_sky]
# z_sky = [-skybox_radius * cos(φ) for φ in φ_sky, θ in θ_sky]
# surface!(lscene, x_sky, y_sky, z_sky,
#     color=milkyway_img,
#     shading=NoShading)
n_lon = 100
n_lat = 100
θ_sphere = range(0, 2π, length=n_lon)  # Longitude
φ_sphere = range(0, π, length=n_lat)   # Latitude

R_earth = 6378.0  # km

# Create meshgrid for sphere
x_earth = [R_earth * sin(φ) * cos(θ) for φ in φ_sphere, θ in θ_sphere]
y_earth = [R_earth * sin(φ) * sin(θ) for φ in φ_sphere, θ in θ_sphere]
z_earth = [R_earth * cos(φ) for φ in φ_sphere, θ in θ_sphere]

surface!(lscene, x_earth, y_earth, z_earth, 
    color=earth_img,  # Use texture image
    shading=FastShading)

# Add equatorial plane (xy plane at z=0) with two-level grid
plane_size = 1250000  # km, extends beyond typical orbits
fine_spacing = 10000  # km between fine grid lines
coarse_spacing = 50000  # km between coarse grid lines (5x fine spacing)

n_fine = Int(plane_size / fine_spacing)
n_coarse = Int(plane_size / coarse_spacing)

# Draw fine grid lines (skip positions where coarse lines will be drawn)
for i in -n_fine:n_fine
    # Skip if this position will have a coarse line
    if i * fine_spacing % coarse_spacing == 0
        continue
    end
    
    y_grid = i * fine_spacing
    x_line = [-plane_size, plane_size]
    y_line = [y_grid, y_grid]
    z_line = [0.0, 0.0]
    lines!(lscene, x_line, y_line, z_line, 
        color=RGBAf(0, 1, 0, 0.15),  # Green, more transparent
        linewidth=0.4)
end

for i in -n_fine:n_fine
    # Skip if this position will have a coarse line
    if i * fine_spacing % coarse_spacing == 0
        continue
    end
    
    x_grid = i * fine_spacing
    x_line = [x_grid, x_grid]
    y_line = [-plane_size, plane_size]
    z_line = [0.0, 0.0]
    lines!(lscene, x_line, y_line, z_line, 
        color=RGBAf(0, 1, 0, 0.15),  # Green, more transparent
        linewidth=0.4)
end

# Draw coarse grid lines (slightly brighter, same z-level as fine grid)
for i in -n_coarse:n_coarse
    y_grid = i * coarse_spacing
    x_line = [-plane_size, plane_size]
    y_line = [y_grid, y_grid]
    z_line = [0.0, 0.0]
    lines!(lscene, x_line, y_line, z_line, 
        color=RGBAf(0, 1, 0, 0.32),  # Green, just slightly brighter than fine lines
        linewidth=0.5)
end

for i in -n_coarse:n_coarse
    x_grid = i * coarse_spacing
    x_line = [x_grid, x_grid]
    y_line = [-plane_size, plane_size]
    z_line = [0.0, 0.0]
    lines!(lscene, x_line, y_line, z_line, 
        color=RGBAf(0, 1, 0, 0.22),  # Green, just slightly brighter than fine lines
        linewidth=0.5)
end

println("Displaying 3D trajectory figure...")
if SCREEN_3D[] === nothing || !GLMakie.is_displayed(SCREEN_3D[])
    SCREEN_3D[] = display(GLMakie.Screen(), fig)
else
    empty!(SCREEN_3D[])
    display(SCREEN_3D[], fig)
end
sleep(0.1)

# Velocity magnitude plot for first 8 segments
println("Creating velocity plot...")
fig2 = Figure(size=(1000, 600), backgroundcolor=:gray95)
ax = Axis(fig2[1, 1], 
    xlabel="Time (MJD)", 
    ylabel="Velocity (km/s)",
    title="Velocity Magnitude per Segment",
    backgroundcolor=:white)

# Better colors for white background
plot_colors = [:blue, :red, :darkgreen, :purple, :darkorange, :brown, :deeppink, :navy]

# Extract velocity data from first 8 segments
for (i, segment) in enumerate(sat.history[1:min(8, length(sat.history))])
    times = Float64[]
    vx_vals = Float64[]
    vy_vals = Float64[]
    vz_vals = Float64[]
    vmag_vals = Float64[]
    
    for entry in segment
        time = entry[1]
        state = entry[2]
        push!(times, time.mjd)  # Extract MJD from Time struct
        push!(vx_vals, state[4])
        push!(vy_vals, state[5])
        push!(vz_vals, state[6])
        push!(vmag_vals, sqrt(state[4]^2 + state[5]^2 + state[6]^2))
    end
    
    # Plot velocity magnitude with better colors for white background
    color = plot_colors[mod1(i, length(plot_colors))]
    lines!(ax, times, vmag_vals, 
        label="Segment $i",
        color=color,
        linewidth=2.5)
end

axislegend(ax, position=:rt)
if SCREEN_2D[] === nothing || !GLMakie.is_displayed(SCREEN_2D[])
    SCREEN_2D[] = display(GLMakie.Screen(), fig2)
else
    empty!(SCREEN_2D[])
    display(SCREEN_2D[], fig2)
end