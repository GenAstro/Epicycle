# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

# Module-level texture cache to avoid reloading
const TEXTURE_CACHE = Dict{String, Any}()

"""
    load_texture_cached(texture_path::String)

Load a texture from the specified path, using cache to avoid repeated loads.
Converts to RGB{N0f8} format for efficient GPU upload.

# Arguments
- `texture_path::String` - Full path to texture file (provided by CelestialBody.texture_file)
"""
function load_texture_cached(texture_path::String)
    # Check cache first
    if haskey(TEXTURE_CACHE, texture_path)
        return TEXTURE_CACHE[texture_path]
    end
    
    # Load from disk (path provided by CelestialBody)
    img = load(texture_path)
    
    # Convert to efficient format (RGB with 8-bit channels)
    # This prevents Float32 conversion and reduces memory by ~4x
    texture = RGB{N0f8}.(img)
    
    # Cache it
    TEXTURE_CACHE[texture_path] = texture
    
    return texture
end

"""
    render_body!(lscene::LScene, coord_sys::CoordinateSystem)

Render central celestial body with texture mapping from CelestialBody.texture_file field.

# Arguments
- `lscene::LScene` - Makie scene to render into
- `coord_sys::CoordinateSystem` - Coordinate system determining which body to render

# Texture Handling
- Texture path from CelestialBody.texture_file field
- Empty texture_file → solid color fallback
- Texture loading cached for performance

# V1.0 Rendering
- Uses GeometryBasics.Sphere for proper UV mapping at poles
- FastShading for performance
- Sphere mesh with equatorial radius from body
"""
function render_body!(lscene::LScene, coord_sys::CoordinateSystem)
    # Get body from coordinate system origin
    body = coord_sys.origin
    body_name = body.name
    
    # Get radius from the celestial body
    R = body.equatorial_radius  # km
    
    # Use GeometryBasics.Sphere with proper UV coordinates
    # This avoids pole singularities that cause texture artifacts
    sphere = GeometryBasics.Sphere(GeometryBasics.Point3f(0, 0, 0), Float32(R))
    
    # Texture path from body field (empty string if none)
    texture_path = body.texture_file
    
    # Try to load texture, fall back to solid gray if unavailable or fails
    color = :gray
    if !isempty(texture_path)
        try
            color = load_texture_cached(texture_path)
        catch e
            @warn "Failed to load texture from $texture_path: $e. Using solid gray."
        end
    end
    
    mesh!(lscene, sphere,
        color=color,
        shading=FastShading,
        interpolate=!isa(color, Symbol))  # Only interpolate for textures
    
    return nothing
end
