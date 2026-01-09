# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

"""
    View3D

Container for 3D visualization of spacecraft orbits and celestial bodies.

# Fields
- `coord_sys::CoordinateSystem` - Coordinate system for the view (default: Earth GCRF)
- `spacecraft::Vector{Spacecraft}` - Spacecraft to visualize
- `_scene::Union{Nothing, Scene}` - Internal Makie scene (private)

# V1.0 Limitations
- Single spacecraft only
- All visual elements always on (no customization)

# Examples
```julia
# Use default Earth GCRF coordinate system
view = View3D()
add_spacecraft!(view, my_spacecraft)
display(view)

# Specify different coordinate system
view = View3D(coord_sys=CoordinateSystem(moon, ICRFAxes()))
add_spacecraft!(view, lunar_spacecraft)
display(view)
```
"""
mutable struct View3D
    coord_sys::CoordinateSystem
    spacecraft::Vector{Spacecraft}
    options::Dict{Spacecraft, Dict{Symbol, Any}}
    _scene::Union{Nothing, Scene}
end

"""
    View3D(; coord_sys::CoordinateSystem = CoordinateSystem(earth, ICRFAxes()))

Create a new 3D visualization view with a default or specified coordinate system.

# Keyword Arguments
- `coord_sys::CoordinateSystem` - Coordinate system for the view (default: Earth GCRF)

# Returns
- `View3D` - Empty view ready to add spacecraft

# Examples
```julia
# Default Earth GCRF
view = View3D()

# Moon-centered inertial
view = View3D(coord_sys=CoordinateSystem(moon, ICRFAxes()))
```
"""
function View3D(; coord_sys::CoordinateSystem = CoordinateSystem(earth, ICRFAxes()))
    View3D(coord_sys, Spacecraft[], Dict{Spacecraft, Dict{Symbol, Any}}(), nothing)
end

"""
    Base.show(io::IO, ::MIME"text/plain", view::View3D)

Display View3D showing coordinate system and number of spacecraft.
"""
function Base.show(io::IO, ::MIME"text/plain", view::View3D)
    println(io, "View3D:")
    println(io, "  Coordinate System = ", view.coord_sys.origin.name, " ", typeof(view.coord_sys.axes))
    println(io, "  Spacecraft        = ", length(view.spacecraft), " spacecraft")
end

"""
    Base.show(io::IO, view::View3D)

Delegate generic show to text/plain variant for print/println.
"""
function Base.show(io::IO, view::View3D)
    show(io, MIME"text/plain"(), view)
end

"""
    add_spacecraft!(view::View3D, sc::Spacecraft)

Add a spacecraft to the 3D view.

# Arguments
- `view::View3D` - View to add spacecraft to
- `sc::Spacecraft` - Spacecraft with history to visualize

# V1.0 Limitations
- Only supports single spacecraft (replacing existing if called multiple times)
- Spacecraft must have history populated

# Throws
- `ArgumentError` if spacecraft history is empty

# Field-Level Validation
This function validates that the spacecraft has data to visualize. Coordinate system
matching is validated at display time.

# Examples
```julia
view = View3D()
add_spacecraft!(view, mission_sat)
```
"""
function add_spacecraft!(view::View3D, sc::Spacecraft; show_iterations::Bool=false)
    # Field-level validation: check spacecraft has data
    if isempty(sc.history)
        error("Spacecraft history is empty. Cannot visualize spacecraft without trajectory data.")
    end
    
    # V1.0: Single spacecraft only (append to vector, but display will use first)
    push!(view.spacecraft, sc)
    
    # Store options for this spacecraft
    view.options[sc] = Dict(:show_iterations => show_iterations)
    
    return view
end

"""
    Base.display(view::View3D; size=(800, 600))

Display the 3D visualization.

# Arguments
- `view::View3D` - View to display

# Keyword Arguments
- `size::Tuple{Int,Int}` - Window size in pixels (default: (800, 600))

# V1.0 Features
- Black background with space theme
- Camera controls (rotation, zoom)
- All elements always visible:
  - Spacecraft trajectory
  - Central body with texture (Earth, Moon, or Mars based on coord_sys)
  - Star field
  - Equatorial plane
  - Spacecraft model (if available)
  - Spacecraft label

# Coupling Validation
Validates that spacecraft coordinate system matches the view coordinate system.

# Examples
```julia
view = View3D()
add_spacecraft!(view, my_spacecraft)
display_view(view)
```
"""
function display_view(view::View3D; size=(800, 600))
    # Coupling validation: check spacecraft added
    if isempty(view.spacecraft)
        error("No spacecraft added to view. Use add_spacecraft!(view, spacecraft) first.")
    end
    
    # Validate all spacecraft coordinate systems match view
    for sc in view.spacecraft
        sc_origin_name = sc.coord_sys.origin.name
        view_origin_name = view.coord_sys.origin.name
        sc_axes_type = typeof(sc.coord_sys.axes)
        view_axes_type = typeof(view.coord_sys.axes)
        
        if sc_origin_name != view_origin_name || sc_axes_type != view_axes_type
            error("Spacecraft '$(sc.name)' coordinate system ($sc_origin_name $sc_axes_type) " *
                  "does not match View3D coordinate system ($view_origin_name $view_axes_type)")
        end
    end
    
    # Create figure with space theme
    fig = Figure(size=size, backgroundcolor=:black)
    lscene = LScene(fig[1, 1],
        show_axis=false,
        scenekw=(
            backgroundcolor=:black,
            clear=true,
            camera=cam3d!
        ))
    
    # Disable auto-centering so our camera positioning works
    cc = cameracontrols(lscene.scene)
    cc.settings.center = false
    cc.settings.zoom_shift_lookat = false
    
    # Render all spacecraft trajectories, models, and labels
    for sc in view.spacecraft
        render_trajectory!(lscene, view, sc)
        render_spacecraft_model!(lscene, sc)
        render_spacecraft_label!(lscene, sc)
    end
    
    # Render scene elements (once)
    render_body!(lscene, view.coord_sys)
    render_stars!(lscene)
    render_equatorial_plane!(lscene)
    
    # Set camera to reasonable default distance
    # Calculate max position magnitude from ALL spacecraft trajectories
    max_r = 0.0
    for sc in view.spacecraft
        for segment in sc.history
            for state in segment.states
                r = sqrt(state.position[1]^2 + state.position[2]^2 + state.position[3]^2)
                max_r = max(max_r, r)
            end
        end
    end
    
    # Set camera distance to 2x max orbital radius
    cam_distance = 2.0 * max_r
    
    # If origin is a celestial body, ensure camera is at least 3x body radius away
    if view.coord_sys.origin isa CelestialBody
        body_radius = view.coord_sys.origin.equatorial_radius
        min_distance = 5.0 * body_radius
        cam_distance = max(cam_distance, min_distance)
    end
    
    # Set camera position using update_cam! (from Makie docs)
    # Position camera at 45 degree angle at calculated distance
    d = cam_distance / sqrt(3)  # Distance along each axis for isometric view
    update_cam!(lscene.scene, Vec3f(d, d, d), Vec3f(0, 0, 0), Vec3f(0, 0, 1))
    
    # Store scene reference
    view._scene = lscene.scene
    
    # Display the figure
    Base.display(fig)
    
    return nothing
end
