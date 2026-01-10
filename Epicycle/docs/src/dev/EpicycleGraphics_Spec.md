
# EpicycleGraphics Design Specification

## Table of Contents
1. [Features & Roadmap](#features--roadmap)
2. [Architecture Principles](#architecture-principles)
3. [Data Ownership](#data-ownership)
4. [Core API Design](#core-api-design)
5. [Attributes & Customization](#attributes--customization)
6. [Default Behaviors](#default-behaviors)
7. [Use Cases & Examples](#use-cases--examples)

---

## Features & Roadmap

### Objects
- Spacecraft (with 3D models)
- Stars (from AstroUniverse)
- Celestial Bodies (Earth, Moon, Sun, planets)
- Ground Stations
- Labels (for spacecraft and bodies)

### Overlays
- Trajectory (multi-segment support)
- Reference Planes (equatorial, ecliptic, orbital)
- Sensor Cones (field of view)
- Visibility Masks
- Reference Axes (inertial, body-fixed, LVLH)
- Vectors (velocity, sun direction, custom)
- Milky Way skybox

### Attributes
- Colors (RGBA with alpha channel)
- Per-segment trajectory colors (for multi-phase missions, maneuvers)
- Coordinate systems (from AstroFrames: CoordinateSystem with origin and axes types)
- Coordinate frames (for axes display)
- Line widths, sizes, scales
- Transparency/alpha values
- Labels (text, color, size)
- Iteration trajectory display (for solver diagnostics)

### Animation
- Real-time playback
- Playback controls: play/pause/stop/step
- Variable speed (faster/slower than real-time)
- Time display and scrubbing
- FPS control

### Version 1.0 MVP - Minimal Viable Product
**Core Feature: Single Spacecraft Visualization (Use Case 0)**
- Single spacecraft with trajectory from history
- Earth with texture (always visible)
- Star field background (always visible)
- Equatorial plane with grid (always visible)
- Spacecraft 3D model (if present in spacecraft.model)
- Spacecraft marker (simple, always visible)
- Spacecraft label showing spacecraft.name (always visible)
- Basic camera controls (GLMakie defaults)
- GCRF coordinate system only (validated, errors on non-GCRF)

**Explicitly OUT of V1.0:**
- Multiple spacecraft
- Customization (colors, line widths, sizes, alpha, etc.)
- User control of stars/plane visibility
- Animation
- Multiple celestial bodies
- Overlays (vectors, axes, custom planes)
- Solver/optimization integration
- 2D plotting (Use Case 7)
- Multiple coordinate systems (V1.1+)
- Ground stations, sensor cones, custom vectors

---

## Architecture Principles

### 1. Separation of Concerns
- **Domain data lives in domain modules** (AstroModels, AstroUniverse)
- **Visualization metadata on domain objects** where sensible (textures, models)
- **EpicycleGraphics handles scene composition and rendering only**

### 2. Default-Driven Design
- **Zero configuration produces professional results** - smart defaults
- **Customization is explicit and clear** - setter functions
- **Convention over configuration** - sensible choices for common cases

### 3. Composability
- Builder pattern for scene construction
- Immutable domain objects
- Mutable view configuration via setters

### 4. Epicycle Consistency
- Kwargs for flexibility
- Time handling via Time structs
- History-based data extraction

### 5. Appropriate Abstraction Levels
- **Custom APIs for domain-specific features** - spacecraft data extraction, trajectory rendering, color handling
- **Standard APIs for standard features** - leverage Makie directly for axis labels, legends, grids, fonts, styling
- **ColorLike everywhere** - unified color system across 3D and 2D
- **Less code to maintain, better documentation leverage** - don't reinvent wheels that Makie already provides

---

## Data Ownership

### AstroFrames
**Coordinate systems used by View3D:**
```julia
# CoordinateSystem from AstroFrames defines origin and axes
struct CoordinateSystem
    origin::CelestialBody  # e.g., Earth, Sun, Moon
    axes::AxesType         # e.g., GCRF(), ITRF(), etc.
end

# Common axes types:
# - GCRF() - Geocentric Celestial Reference Frame (inertial)
# - ITRF() - International Terrestrial Reference Frame (Earth-fixed)
# - Other axes types as defined in AstroFrames
```

### AstroModels
**Additions to existing Spacecraft struct:**
```julia
struct Spacecraft
    # ... existing fields ...
    history::Vector{Vector{Tuple{Time{Float64}, Vector{Float64}}}}
    model::Union{SpacecraftModel, Nothing}  # NEW - 3D model metadata
end

struct SpacecraftModel
    mesh_path::String      # Path to .obj/.stl/.ply file
    scale::Float64         # Default size scaling
    offset::Vec3f          # Model origin offset (optional)
    rotation::Quaternion   # Default orientation (optional)
end
```

### AstroUniverse
**Additions to existing CelestialBody:**
```julia
struct CelestialBody
    # ... existing fields (name, radius, μ, etc.) ...
    texture_path::Union{String, Nothing}  # NEW - texture image path
    rotation_elements::Union{RotationElements, Nothing}  # For future animation
end
```

**Star field data:**
```julia
struct StarField
    positions::Vector{Vec3f}     # Pre-computed positions on celestial sphere
    magnitudes::Vector{Float64}  # Brightness values
    colors::Vector{RGB}          # Star colors (white/yellow/blue based on type)
end

struct Universe
    # ... existing fields ...
    stars::Union{StarField, Nothing}  # NEW - star catalog
end
```

### EpicycleGraphics
**Visualization-only structures:**
```julia
struct View3D
    # Object references (domain data)
    spacecraft::Vector{Spacecraft}
    celestial_bodies::Vector{CelestialBody}
    
    # Visualization primitives (no domain data)
    planes::Vector{Plane3D}
    vectors::Vector{Vector3D}
    axes::Vector{Axes3D}
    
    # Universe reference
    universe::Universe
    
    # Display preferences (per-object customization)
    options::Dict{Any, Dict{Symbol, Any}}
    
    # Scene-level settings
    show_stars::Bool
    background::Symbol  # :black, :space
    coordinate_system::CoordinateSystem  # From AstroFrames (default: GCRF origin with GCRF axes)
    
    # Internal rendering state (private)
    _scene::Union{Nothing, Scene}
end

struct Plane3D
    type::Symbol  # :equatorial, :ecliptic, :orbital, :custom
    reference_body::Union{CelestialBody, Nothing}
    color::RGBA
    grid_spacing::Union{Tuple{Float64, Float64}, Nothing}  # (fine, coarse)
end

struct Vector3D
    name::String
    source::Spacecraft  # Which object it's attached to
    vector_type::Symbol  # :velocity, :sun, :nadir, :custom
    color::RGBA
    scale::Float64
    vector_func::Union{Function, Nothing}  # (time) -> Vec3 for custom
end

struct Axes3D
    frame::Symbol  # :inertial, :body, :lvlh, :orbital
    reference::Union{Spacecraft, CelestialBody}
    colors::Tuple{RGBA, RGBA, RGBA}  # X, Y, Z axis colors
    length::Float64
    linewidth::Float64
end
```

---

## Core API Design

### Scene Construction (Builder Pattern)

```julia
# Create view with universe reference
view = View3D(universe)  # Default: GCRF origin with GCRF() axes
# view = View3D(universe; coordinate_system=CoordinateSystem(Earth, ITRF()))  # ECEF

# Add objects (plural functions accept multiple objects)
add_spacecraft!(view, mysat, yoursat)
add_bodies!(view, earth, moon)
add_planes!(view, (:equatorial, earth), (:ecliptic, sun))
add_vectors!(view, (mysat, :velocity), (yoursat, :sun))
add_axes!(view, (mysat, :body), (earth, :inertial))

# Can also add single objects
add_spacecraft!(view, another_sat)
add_bodies!(view, mars)

# Remove objects
remove_spacecraft!(view, mysat)
remove_bodies!(view, moon)

# Update display (refresh after modifying object data)
update_display!(view)

# Display with rendering options
display(view; 
    size = (800, 600),
    camera_position = :auto,  # or Vec3f(x,y,z)
    background = :black,
    time = nothing  # Optional: specific time to display, or nothing for latest/current
)

# Animate with playback options
animate(view;
    start_time = nothing,  # Auto from earliest in history
    end_time = nothing,    # Auto from latest in history
    speed = 10.0,          # 10x real-time
    show_time = true,      # Display time overlay
    fps = 30
)
# If time extends beyond history, spacecraft simply won't be shown for those times
```

### Separation of "What" from "How"
- **Scene content (WHAT)**: Objects added to View3D
- **Rendering options (HOW)**: Kwargs to `display()` and `animate()`
- **Time range**: Auto-detected from spacecraft history

---

## Attributes & Customization

### Philosophy
1. **Smart defaults make common cases beautiful**
2. **Kwargs-first for initial setup** - declare attributes when adding objects
3. **Setters for later changes** - update attributes after adding (solver iterations, interactive changes)
4. **Flexible color input, consistent internal representation** (ColorLike → RGBA{Float32})
5. **First-class multi-segment support** - scalar for uniform, vector for per-segment

### Kwargs Pattern (Primary)

```julia
# Single spacecraft - simple attributes
add_spacecraft!(view, sat1;
    color = :cyan,
    linewidth = 2.5,
    show_model = true,
    model_scale = 1000.0,
    show_label = true,  # Shows sat1.name
    show_iterations = false)  # Default: hide solver iterations

# Single spacecraft - per-segment attributes (multi-phase mission)
add_spacecraft!(view, transfer_sat;
    color = [:red, :orange, :yellow, :green],      # Per-segment colors
    linewidth = [3.0, 1.5, 1.5, 3.0],              # Thick for burns, thin for coasts
    alpha = [1.0, 0.8, 0.6, 0.4],                  # Progressive fade
    show_model = true,
    show_label = true)  # Shows transfer_sat.name

# Multiple spacecraft - different configurations
add_spacecraft!(view,
    sat1,  # Uses defaults (auto color from palette)
    sat2;  # Also uses defaults
    linewidth = 2.0)  # Applied to both

add_spacecraft!(view, sat3;
    color = :magenta,
    linewidth = 3.5,
    show_model = false,
    label = "Debris")

# Full example with all trajectory attributes
add_spacecraft!(view, mission_sat;
    # Trajectory appearance
    color = [:red, :orange, :green],     # Or single: :cyan
    linewidth = [3.0, 2.0, 2.0],         # Or single: 2.5
    alpha = [1.0, 0.8, 0.6],             # Or single: 1.0
    
    # Model display
    show_model = true,
    model_scale = 2000.0,
    
    # Marker (for distant views)
    show_marker = true,                  # Or :auto for automatic
    marker = :circle,
    marker_size = 10.0,
    marker_color = :cyan,                # Or :auto to match trajectory
    
    # Label
    show_label = true,                   # Shows mission_sat.name
    label_color = :white,
    label_size = 14.0,
    
    # Advanced
    show_attitude = false,
    show_trajectory = true,
    show_iterations = false)             # Solver iteration diagnostics
```

### Iteration Display

**Purpose**: Visualize solver convergence by showing iteration history from `solve_trajectory!(...; record_iterations=true)`

```julia
# Default: show only final trajectory (from spacecraft.history.segments)
add_spacecraft!(view, sat;
    color = :cyan,
    linewidth = 2.5)

# Enable iteration display (from spacecraft.history.iterations)
add_spacecraft!(view, sat;
    color = :cyan,
    linewidth = 2.5,
    show_iterations = true)  # Shows solver iterations behind final trajectory

# Typical diagnostic workflow:
# 1. Run solve_trajectory! with record_iterations=true
result = solve_trajectory!(seq, options; record_iterations=true)

# 2. Visualize convergence
view = View3D(universe)
add_spacecraft!(view, sat; show_iterations=true)
display(view)  # See how solver converged to solution
```

**Implementation details:**
- Iterations rendered first (background layer) with default styling
- Final trajectory rendered on top (foreground layer)
- Default iteration style: gray, thin, semi-transparent (hard-coded for V1.0)
- Only shown if `spacecraft.history.iterations` contains data
- If `show_iterations=true` but no iteration data exists, silently ignored

**Future enhancements (V1.1+):**
- `iteration_color`: Custom color for iteration trajectories
- `iteration_linewidth`: Custom line width
- `iteration_alpha`: Custom transparency
- `show_progression`: Progressive styling showing convergence pattern

### Setter Functions (Secondary - for changes after adding)

```julia
# Generic attribute setter
set_attribute!(view, object, attribute::Symbol, value)

# Color convenience setter  
set_color!(view, object, color)  # For trajectory color

# Use cases for setters:
# 1. Changing attributes after initial setup
set_attribute!(view, mysat, :linewidth, 5.0)
set_color!(view, mysat, :red)  # Change trajectory color

# 2. Solver iterations - update colors during optimization
for iteration in solver
    set_color!(view, current_solution, :gray)  # Dim old solution
    # ... compute new solution ...
    set_color!(view, new_solution, :cyan)      # Highlight new
    update_display!(view)
end

# 3. Per-segment colors for complex missions
set_color!(view, mysat, [
    RGBA(1, 0, 0, 1),     # Phase 1: red
    RGBA(1, 0.5, 0, 1),   # Phase 2: orange  
    RGBA(0, 1, 0, 1)      # Phase 3: green
])
```

### Kwargs Defaults

Smart defaults ensure zero-configuration produces professional results:

**Spacecraft:**
```julia
DEFAULT_SPACECRAFT_KWARGS = (
    color = :auto,              # Auto-cycle from palette
    linewidth = 2.0,
    alpha = 1.0,
    show_trajectory = :auto,    # true if history exists
    show_model = :auto,         # true if spacecraft.model !== nothing
    show_attitude = false,
    show_marker = :auto,        # true when camera far from spacecraft (model too small)
    marker = :circle,           # Shape: :circle, :rect, :diamond, :cross, :xcross, etc.
    marker_size = 8.0,          # Marker size in pixels
    marker_color = :auto,       # Same as trajectory color if not specified
    show_label = false,         # Display spacecraft.name as label
    label_color = :white,
    label_size = 12.0
)
```

**Celestial Bodies:**
```julia
DEFAULT_BODY_KWARGS = (
    show_texture = :auto,       # true if body.texture_path !== nothing
    show_axes = false,
    show_grid = false,
    grid_spacing = :auto,       # Based on body radius
    show_label = false,         # Display body.name as label
    label_color = :white,
    label_size = 14.0
)
```

**Planes:**
```julia
DEFAULT_PLANE_KWARGS = (
    color = :auto,              # Depends on plane type (equatorial=green, etc.)
    alpha = 0.2,
    grid_spacing = :auto,       # Based on orbit size
    extent = :auto              # Auto to encompass trajectories
)
```

**Vectors:**
```julia
DEFAULT_VECTOR_KWARGS = (
    color = :cyan,
    scale = 1.0,
    linewidth = 2.0,
    arrowsize = :auto           # Based on vector length
)
```

**Auto-cycling color palette for spacecraft:**
```julia
DEFAULT_TRAJECTORY_COLORS = [
    :cyan,      # RGBA(0.0, 1.0, 1.0, 1.0)
    :magenta,   # RGBA(1.0, 0.0, 1.0, 1.0)
    :yellow,    # RGBA(1.0, 1.0, 0.0, 1.0)
    :green,     # RGBA(0.0, 1.0, 0.0, 1.0)
    :orange,    # RGBA(1.0, 0.5, 0.0, 1.0)
    :red        # RGBA(1.0, 0.0, 0.0, 1.0)
]
# Cycles: sat1=cyan, sat2=magenta, sat3=yellow, sat4=green, sat5=orange, sat6=red, sat7=cyan, ...
```

### Available Attributes (Kwargs)

#### Spacecraft
**Trajectory appearance:**
- `color` → ColorLike or Vector{ColorLike} (default: auto from palette, cycles through DEFAULT_TRAJECTORY_COLORS)
- `linewidth` → Float64 or Vector{Float64} (default: 2.0)
- `alpha` → Float64 or Vector{Float64} (default: 1.0)
- `show_trajectory` → Bool (default: true if history exists)

**Model display:**
- `show_model` → Bool (default: true if spacecraft.model !== nothing)
- `model_scale` → Float64 (default: 1.0, multiplies spacecraft.model.scale)
- `show_attitude` → Bool (default: false)

**Marker (for distant views):**
- `show_marker` → Bool (default: :auto - true when model too small to see)
- `marker` → Symbol (default: :circle, options: :rect, :diamond, :cross, :xcross, :utriangle, :dtriangle, etc.)
- `marker_size` → Float64 (default: 8.0 pixels)
- `marker_color` → ColorLike (default: :auto - same as trajectory color)

**Label:**
- `show_label` → Bool (default: false - displays spacecraft.name from domain object)
- `label_color` → ColorLike (default: :white)
- `label_size` → Float64 (default: 12.0)

#### Celestial Body
**Appearance:**
- `show_texture` → Bool (default: true if body.texture_path !== nothing)
- `show_axes` → Bool (default: false)
- `show_grid` → Bool (default: false)
- `grid_spacing` → Tuple{Float64, Float64} (default: auto based on body radius)

**Label:**
- `show_label` → Bool (default: false - displays body.name from domain object)
- `label_color` → ColorLike (default: :white)
- `label_size` → Float64 (default: 14.0)

#### Plane
- `color` → ColorLike (default: depends on type - see Reference Planes defaults)
- `alpha` → Float64 (default: 0.2)
- `grid_spacing` → Tuple{Float64, Float64} (default: auto based on orbit size)
- `extent` → Float64 (default: auto to encompass trajectories)

#### Vector
- `color` → ColorLike (default: :cyan)
- `scale` → Float64 (default: 1.0)
- `linewidth` → Float64 (default: 2.0)
- `arrowsize` → Float64 (default: auto based on vector length)

---

## Default Behaviors

### Kwargs-First Philosophy
- **Primary configuration**: kwargs at `add_*!()` time (declarative)
- **Secondary updates**: setters after adding (imperative, for solver iterations)
- **Smart defaults**: zero configuration produces professional results
- **Per-segment support**: scalar for uniform, vector for per-segment

### Color Palette (Auto-Cycling)
When spacecraft are added without explicit `color` kwarg, colors auto-cycle from:
```julia
DEFAULT_TRAJECTORY_COLORS = [:cyan, :magenta, :yellow, :green, :orange, :red]
```
First spacecraft gets cyan, second gets magenta, etc. Cycles back to cyan after red.

### Visibility Defaults
- Show trajectory if spacecraft has history (`show_trajectory = true`)
- Show model if `spacecraft.model !== nothing` (`show_model = true`)
- Show marker when `show_marker = :auto` and model appears too small (camera far from spacecraft)
  - Threshold: model projected size < ~10 pixels → show marker instead
  - Can force with `show_marker = true` or disable with `show_marker = false`
- Show texture if `body.texture_path !== nothing` (`show_texture = true`)
- Show stars by default (`view.show_stars = true`)
- Don't show attitude/axes by default (reduces clutter: `show_attitude = false`)
- Don't show labels by default (add via `label` kwarg)

### Attribute Polymorphism (Scalar vs Vector)
- **Scalar value**: applies to all segments uniformly
  ```julia
  add_spacecraft!(view, sat; color = :cyan, linewidth = 2.0)
  ```
- **Vector value**: per-segment (must match number of history segments)
  ```julia
  add_spacecraft!(view, sat; 
      color = [:red, :orange, :yellow],
      linewidth = [3.0, 2.0, 2.0])
  ```

### Coordinate System Defaults
- Default coordinate system: `CoordinateSystem(Earth, GCRF())` (Earth-centered inertial)
- From AstroFrames: CoordinateSystem has origin (e.g., Earth, Sun) and axes type (e.g., GCRF(), ITRF())
- All position/velocity data in spacecraft history assumed to be in this frame
- Alternative: `CoordinateSystem(Earth, ITRF())` - Earth-fixed, rotating with Earth
- Coordinate system affects:
  - How trajectories are displayed
  - Reference plane orientations
  - Star field rotation (for ITRF)
  - Ground station positions

### Scale Defaults
- Camera auto-scales to fit the trajectories (not all object like far way planets)
- Trajectory linewidth: 2.0
- Model scale: 1.0 (assumes models are pre-scaled)
- Plane grid spacing: auto based on orbit size
- Star field radius: 10,000,000 km (large enough to appear infinite)

### Reference Planes
- Equatorial: Green, RGBA(0, 1, 0, 0.2)
- Ecliptic: Blue, RGBA(0, 0.5, 1, 0.2)
- Orbital: Gray, RGBA(0.5, 0.5, 0.5, 0.2)

### Time Handling
- Time range specified in `animate()` kwargs (start_time, end_time)
- Default: use union of all spacecraft history time ranges
- If display time is outside spacecraft history, spacecraft not rendered
- Interpolation between history points as needed
- Bodies and stars always rendered (independent of time)

---

## 2D Plotting Integration

### Time-Series API with AstroCallbacks

EpicycleGraphics integrates with AstroCallbacks to provide type-safe extraction of time-series data from spacecraft history and celestial bodies.

**Primary API - Full Calc Struct:**
```julia
plot_timeseries!(ax::Axis, calc::AbstractCalc; color, linewidth, label, kwargs...)
```
Accepts any Calc type from AstroCallbacks. This is the most flexible form and supports:
- `OrbitCalc(spacecraft, var_tag; dependency=nothing, coordinate_system=nothing)`
- `BodyCalc(body, var_tag)`
- `ManeuverCalc(maneuver, spacecraft, var_tag)`

**Convenience Dispatch - Simple OrbitCalc:**
```julia
plot_timeseries!(ax::Axis, spacecraft::Spacecraft, var::AbstractOrbitVar; kwargs...)
```
Automatically constructs `OrbitCalc(spacecraft, var)` for simple cases without optional kwargs.

**How it works:**
1. Extracts time values from spacecraft history (or appropriate source)
2. Evaluates `get_calc(calc)` at each time point
3. Plots using Makie's `lines!()` or `scatter!()` with provided attributes
4. Returns plot object for legend integration

**Supported kwargs:**
- `color` → ColorLike (default: :blue)
- `linewidth` → Float64 (default: 2.0)
- `linestyle` → Symbol (default: :solid, options: :dash, :dot, :dashdot)
- `label` → String (default: nothing, for legend)
- `alpha` → Float64 (default: 1.0)
- Plus any standard Makie line attributes

**Time extraction:**
- For `OrbitCalc`: Extract times from `spacecraft.history`
- For `BodyCalc`: Use provided time range or match other plots
- For `ManeuverCalc`: Extract from maneuver and spacecraft context

---

## Use Cases & Examples

### Use Case 0: Minimal MVP - Single Spacecraft (V1.0)
```julia
# V1.0 MVP: Minimal viable visualization
# Shows single spacecraft orbit around Earth with stars
# All visual elements always on, no customization

using GLMakie
using EpicycleGraphics

# Spacecraft must be in GCRF (validated, errors otherwise)
# Spacecraft must have history for trajectory display
view = View3D(universe)
add_spacecraft!(view, mission_sat)  # Adds Earth automatically

# Simple display - uses all defaults
# - Stars: always visible
# - Equatorial plane: always visible
# - Earth texture: always visible (if available)
# - Spacecraft model: visible if spacecraft.model exists
# - Spacecraft marker: always visible
# - Spacecraft label: always visible (shows mission_sat.name)
display(view)

# That's it for V1.0 MVP!
# No customization options in this release
# Coordinate system assumed/validated as Earth GCRF
```

### Use Case 1: Simple Mission Visualization (V1.1+)
```julia
# Just works - zero configuration, smart defaults
view = View3D(universe)
add_spacecraft!(view, mission_sat)  # Auto: cyan color, linewidth=2.0, shows model if present
add_bodies!(view, earth)             # Auto: shows texture if present
display(view)

# Or with minimal customization
view = View3D(universe)
add_spacecraft!(view, mission_sat; show_label = true)  # Shows mission_sat.name, other defaults apply
add_bodies!(view, earth)
display(view)
```

### Use Case 2: Scene Attributes (Stars, Planes, Camera) (V1.1+)
```julia
# Control scene-level display options
view = View3D(universe)
add_spacecraft!(view, sat)
add_bodies!(view, earth, moon)
add_planes!(view, (:equatorial, earth))

# Scene attributes set at construction or via setters
view.show_stars = true  # or set in View3D constructor
view.coordinate_system = earth_icrf  # Earth-Centered Inertial (default)
view.show_plane = true

# Display with camera and rendering options
display(view; 
    background = :black,
    size = (1200, 800),
    camera_position = Vec3f(50000, 0, 20000),  # Custom view position
    camera_lookat = Vec3f(0, 0, 0),            # Look at origin
    camera_up = Vec3f(0, 0, 1)                 # Z-up orientation
)

# Or use predefined camera views
display(view; camera_position = :earth_centered)
display(view; camera_position = :spacecraft_following)
```

### Use Case 3: Overlays (Vectors, Axes, Planes) (V1.1+)
```julia
# Add various overlay elements
view = View3D(universe)
add_spacecraft!(view, sat)
add_bodies!(view, earth)

# Add reference planes
add_planes!(view, (:equatorial, earth), (:ecliptic, sun))

# Add vectors
add_vectors!(view, (sat, :velocity), (sat, :sun))

# Add coordinate axes
add_axes!(view, (sat, :body), (earth, :inertial))

display(view)
```

### Use Case 4: Models and Textures (V1.1+)
```julia
# Spacecraft models are defined on the Spacecraft struct in AstroModels
# Model is automatically displayed if present
sc_with_model = Spacecraft(
    # ... spacecraft parameters ...
    model = SpacecraftModel(
        mesh_path = "path/to/spacecraft.obj",
        scale = 1000.0,  # Scale factor for model size
        offset = Vec3f(0, 0, 0),
        rotation = Quaternion(1, 0, 0, 0)  # Identity rotation
    )
)

view = View3D(universe)
add_spacecraft!(view, sc_with_model)
add_bodies!(view, earth)

# Model is shown by default, can be controlled via attributes
set_attribute!(view, sc_with_model, :show_model, true)
set_attribute!(view, sc_with_model, :model_scale, 2000.0)  # Override scale
set_attribute!(view, sc_with_model, :show_label, true)     # Show spacecraft name

# Celestial body textures are defined in AstroUniverse
# Earth texture is automatically applied if texture_path is set
earth = CelestialBody(
    # ... body parameters ...
    texture_path = "path/to/earth_texture.jpg"
)

# Texture is shown by default
add_bodies!(view, earth)

# Can control texture display
set_attribute!(view, earth, :show_texture, true)
set_attribute!(view, earth, :show_axes, false)

display(view)
```

### Use Case 5: Custom Colors and Attributes (V1.1+)
```julia
# Customize trajectory colors and visual attributes
view = View3D(universe)
add_bodies!(view, earth)
add_planes!(view, (:equatorial, earth))

# Single spacecraft - simple uniform attributes
add_spacecraft!(view, sat1;
    color = :white,
    linewidth = 3.0,
    alpha = 0.8,
    show_label = true)  # Shows sat1.name

# Single spacecraft - per-segment attributes (multi-phase mission)
add_spacecraft!(view, sat2;
    color = [:red, :orange, :yellow, :green],    # Launch, transfer, approach, orbit
    linewidth = [3.0, 2.0, 2.0, 2.5],            # Thicker during burns
    alpha = [1.0, 0.9, 0.8, 0.7],                # Slight fade over time
    show_model = true,
    model_scale = 1500.0,
    show_label = true,                           # Shows sat2.name
    label_color = :cyan,
    label_size = 12.0)

# Multiple spacecraft - mix of default and custom
add_spacecraft!(view, 
    sat3,  # Default: auto color from palette, standard linewidth
    sat4;  # Default: next color from palette
    show_model = false)  # Both don't show models

add_spacecraft!(view, sat5;
    color = RGBColor(255, 128, 0),  # Custom RGB
    linewidth = 4.0,
    marker = :diamond,              # Different marker shape
    marker_size = 12.0,             # Larger marker
    show_label = true)              # Shows sat5.name

display(view)

# Later: change attributes using setters (e.g., during solver iteration)
set_color!(view, sat2, :gray)  # Dim previous solution
set_attribute!(view, sat2, :alpha, 0.3)
update_display!(view)
```

### Use Case 6: Solver/Optimization Workflow (V1.2+)
```julia
# TBD - Showing solver iterations in real-time
# - Display current iteration with subdued color
# - Update display as solver progresses
# - Replace with final solution in bright color
# Details to be refined based on solver integration needs
```

### Use Case 7: 2D Data Plots (V1.1+)
```julia
# Plot time-series data from spacecraft history using AstroCallbacks
# Use Makie's Figure/Axis directly for standard plotting features
# Use Epicycle's plot_timeseries! for domain-specific data extraction with Calc framework
# ColorLike system works seamlessly with Makie
# Note: GLMakie handles both 3D and 2D plotting

using GLMakie
using AstroCallbacks

view = View3D(universe)
add_spacecraft!(view, sat1, sat2, sat3)
add_bodies!(view, earth, moon)

# Create Makie figure and axis directly - no wrappers needed
fig = Figure(size=(1000, 600))
ax = Axis(fig[1, 1],
    xlabel = "Time (MJD)",
    ylabel = "Velocity (km/s)",
    title = "Velocity Comparison")

# PRIMARY API: Full Calc struct (supports all Calc types and optional kwargs)
# OrbitCalc with simple variable tag
plot_timeseries!(ax, OrbitCalc(sat1, VelMag()); 
    color = :cyan,
    label = "Sat 1")

# OrbitCalc with dependency (for relative calculations)
plot_timeseries!(ax, OrbitCalc(sat2, OutgoingRLA(); dependency=sat1);
    color = RGBColor(255, 0, 255),   # RGB integer
    label = "Sat 2 relative to Sat 1")

# OrbitCalc with coordinate system (future feature)
# plot_timeseries!(ax, OrbitCalc(sat3, SMA(); coordinate_system=earth_itrf);
#     color = HexColor("#FFFF00"),
#     label = "Sat 3 SMA in ITRF")

# BodyCalc for celestial body properties
plot_timeseries!(ax, BodyCalc(earth, GravParam());
    color = :green,
    label = "Earth μ")

# CONVENIENCE DISPATCH: Simple OrbitCalc cases (no optional kwargs)
# Automatically constructs OrbitCalc(spacecraft, var_tag)
plot_timeseries!(ax, sat1, VelMag();
    color = :cyan,
    label = "Sat 1 (convenience)")

plot_timeseries!(ax, sat2, SMA();
    color = :magenta,
    label = "Sat 2 SMA")

# ManeuverCalc requires full struct (no convenience dispatch)
plot_timeseries!(ax, ManeuverCalc(toi_maneuver, sat3, DeltaVMag());
    color = :orange,
    label = "TOI ΔV magnitude")

# Makie API: standard plotting features
axislegend(ax; position=:rt)
ax.xgridvisible = true
ax.ygridvisible = true

display(fig)

# Or display 3D view and 2D plot side-by-side using Makie's layout
fig2 = Figure(size=(1600, 600))
# 3D view in left panel (Epicycle renders into Makie axis)
display(view; figure=fig2[1, 1])
# 2D plot in right panel
ax2 = Axis(fig2[1, 2], 
    xlabel = "Time (MJD)", 
    ylabel = "Semi-Major Axis (km)")

# Mix of full Calc and convenience dispatch
plot_timeseries!(ax2, OrbitCalc(sat1, SMA()); color=:cyan, label="Sat 1")
plot_timeseries!(ax2, sat2, SMA(); color=:magenta, label="Sat 2")  # Convenience

axislegend(ax2)
display(fig2)

# ColorLike system is consistent across 3D and 2D
# All accept: Symbol, RGBA, RGBColor(), HexColor(), NamedColor()
# Makie handles these natively or via simple conversion
```

### Use Case 8: Animation (V1.2+)
```julia
# Animate trajectory over time
view = View3D(universe)
add_spacecraft!(view, probe)
add_bodies!(view, earth, mars)
add_planes!(view, (:ecliptic, sun))

# Playback controls
animate(view;
    speed = 10.0,        # 10x real-time
    show_time = true,    # Display time overlay
    fps = 30,
    start_time = nothing,  # Auto from history
    end_time = nothing
)

# Future: Interactive controls
# - play/pause/stop/step
# - time scrubbing
# - speed adjustment during playback
```

---

## Implementation Notes

### Rendering Pipeline
1. **Scene setup**: Create Makie scene with background
2. **Add celestial bodies**: Textured spheres from Universe (always rendered)
3. **Add stars**: Scatter plot from Universe.stars (always rendered)
4. **Add trajectories**: Multi-segment lines with per-segment colors
5. **Add spacecraft models**: Load and transform meshes at specified time
6. **Add overlays**: Planes, vectors, axes
7. **Camera setup**: Auto-scale or user-specified
8. **For animation**: Update time, interpolate spacecraft states, refresh
   - If time outside history range, spacecraft not rendered for that frame
   - Celestial bodies and static elements remain visible

### Multi-Segment Trajectory Rendering
```julia
function render_trajectory!(scene, view, sc)
    segments = extract_trajectory_segments(sc)
    segment_colors = get_attribute(view, sc, :trajectory_colors)
    
    if segment_colors !== nothing
        # Per-segment colors
        for (i, (x, y, z)) in enumerate(segments)
            color = i <= length(segment_colors) ? segment_colors[i] : segment_colors[end]
            lines!(scene, x, y, z, color=color, linewidth=...)
        end
    else
        # Single color
        color = get_attribute(view, sc, :trajectory_color)
        for (x, y, z) in segments
            lines!(scene, x, y, z, color=color, linewidth=...)
        end
    end
end
```

### Camera Control
- Free rotation (default from prototype)
- Predefined views: `:earth_centered`, `:spacecraft_following`, `:inertial`
- Manual position: `Vec3f(x, y, z)`
- Zoom control without shifting lookat point

---

## Open Questions for Review

1. **Animation in V1?** Currently stretch goal - should it be core?
2. **Ground stations in V1?** Or defer to V1.1?
3. **Color specification**: RGBA only, or allow named colors that convert?
4. **Plane grid auto-sizing**: Algorithm for determining good grid spacing?
5. **Multiple viewports**: Single scene only for V1?
6. **Export capabilities**: Save images/videos in V1 or later?
7. **Interactive controls**: Mouse picking objects, interactive time scrubbing?

---

## Version Roadmap

### V1.0 - Minimal MVP (Current Focus - Business Launch)
**Goal: Ship fast, get market feedback, create "wow factor"**
- Single spacecraft visualization (Use Case 0)
- Earth with texture (always on)
- Star field (always on)
- Equatorial plane (always on)
- Spacecraft model, marker, label (always on)
- GCRF coordinate system only
- No customization options
- No animation
- Static display only

**Implementation Priority:**
1. Basic View3D struct
2. Single spacecraft trajectory rendering
3. Earth sphere with texture
4. Star field scatter plot
5. Equatorial plane grid
6. Simple display() function

### V1.1 - Customization & Multi-Object
**Goal: Production-ready with flexibility**
- Multiple spacecraft (Use Case 1)
- Attribute customization (colors, sizes, visibility) (Use Cases 2, 5)
- Multiple celestial bodies
- Multiple reference planes
- User control of stars/plane visibility
- Builder API with kwargs
- Setter functions for updates
- 2D plotting integration (Use Case 7)
- Multi-segment trajectory colors
- Coordinate system selection (GCRF, ITRF, etc.)

### V1.2 - Overlays & Interaction
**Goal: Advanced visualization**
- Vectors and axes overlays (Use Case 3)
- Ground stations
- Sensor cones (FOV)
- Animation controls (Use Case 8)
- Time scrubbing
- Interactive camera presets
- Solver/optimization integration (Use Case 6)
- Object picking

### V2.0 - Advanced Features
**Goal: Professional grade**
- Multiple viewports
- Video export
- Real-time updates
- Custom shaders
- VR support (?)
- Performance optimization


# Uitlities

## Color

```julia
module EpicycleColor

using ColorTypes
using Colors

export ColorLike, RGBColor, HexColor, NamedColor

"Types accepted as colors across Epicycle graphics."
const ColorLike = Union{Symbol, Colorant}

"0–1 float RGB(A)."
RGBColor(r::Real, g::Real, b::Real; a::Real = 1.0) =
    RGBA{Float32}(r, g, b, a)

"0–255 integer RGB(A)."
RGBColor(r::Integer, g::Integer, b::Integer; a::Integer = 255) =
    RGBA{Float32}(r/255, g/255, b/255, a/255)

"Hex string '#RRGGBB' or '#RRGGBBAA' (with or without leading '#')."
function HexColor(s::AbstractString)::Colorant
    s = strip(s)
    s = startswith(s, "#") ? s : "#" * s
    parse(Colorant, s)   # Colors.jl
end

"Semantic or named color symbol (Makie-compatible)."
NamedColor(name::Symbol) = name
NamedColor(name::AbstractString) = Symbol(name)

end # module
```
