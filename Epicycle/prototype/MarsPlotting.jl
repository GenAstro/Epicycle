using GLMakie
using FileIO  # For loading images
using Random  # For reproducible star field
using MeshIO  # For loading 3D model filesex
using GeometryBasics  # For mesh handling

# Screen management - reuse windows instead of creating new ones
const SCREEN_3D = Ref{Union{Nothing,GLMakie.Screen}}(nothing)
const SCREEN_2D = Ref{Union{Nothing,GLMakie.Screen}}(nothing)

# Load texture image
mars_img = load(raw"C:\Users\steve\Dev\Epicycle\Epicycle\prototype\Mars_JPLCaltechUSGS.jpg")
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

# Add textured Mars sphere
n_lon = 100
n_lat = 100
θ_sphere = range(0, 2π, length=n_lon)  # Longitude
φ_sphere = range(0, π, length=n_lat)   # Latitude

R_mars = 3389.5  # km

# Create meshgrid for sphere
x_mars = [R_mars * sin(φ) * cos(θ) for φ in φ_sphere, θ in θ_sphere]
y_mars = [R_mars * sin(φ) * sin(θ) for φ in φ_sphere, θ in θ_sphere]
z_mars = [R_mars * cos(φ) for φ in φ_sphere, θ in θ_sphere]

surface!(lscene, x_mars, y_mars, z_mars, 
    color=mars_img,  # Use texture image
    shading=FastShading)

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

surface!(lscene, x_mars, y_mars, z_mars, 
    color=mars_img,  # Use texture image
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
