# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

"""
    render_spacecraft_model!(lscene::LScene, sc::Spacecraft)

Render spacecraft 3D model if available.

# Arguments
- `lscene::LScene` - Makie scene to render into
- `sc::Spacecraft` - Spacecraft with cad_model field

# V1.0 Rendering
- Loads .obj file from sc.cad_model.file_path if visible
- Positions at current spacecraft state location
- Yellow color (fixed)
- Applies scale from sc.cad_model.scale
- Falls back to simple marker if model unavailable or fails to load

# Supported Formats
- .obj (primary)
- .stl, .ply, .off (supported by FileIO/MeshIO)
"""
function render_spacecraft_model!(lscene::LScene, sc::Spacecraft)
    # Check if spacecraft has visible model
    if isempty(sc.cad_model.file_path) || !sc.cad_model.visible
        # V1.0: Render simple marker as fallback
        render_spacecraft_marker!(lscene, sc)
        return nothing
    end
    
    # Get current spacecraft position (last point in history)
    last_segment = sc.history[end]
    last_state = last_segment.states[end]
    
    x_pos = last_state.position[1]
    y_pos = last_state.position[2]
    z_pos = last_state.position[3]
    
    # Try to load and render the model
    try
        model_mesh = load(sc.cad_model.file_path)
        
        # Render mesh at spacecraft position
        mesh!(lscene, model_mesh,
            color=:yellow,
            transformation=(
                translation=Vec3f(x_pos, y_pos, z_pos),
                scale=Vec3f(sc.cad_model.scale)
            ))
    catch e
        @warn "Failed to load spacecraft model, using marker instead" exception=e file_path=sc.cad_model.file_path
        render_spacecraft_marker!(lscene, sc)
    end
    
    return nothing
end

"""
    render_spacecraft_marker!(lscene::LScene, sc::Spacecraft)

Render simple marker at spacecraft position (fallback when model unavailable).

# Arguments
- `lscene::LScene` - Makie scene to render into
- `sc::Spacecraft` - Spacecraft to mark

# V1.0 Rendering
- Yellow scatter point
- Markersize 8
- Positioned at current spacecraft state
"""
function render_spacecraft_marker!(lscene::LScene, sc::Spacecraft)
    # Get current spacecraft position
    last_segment = sc.history[end]
    last_state = last_segment.states[end]
    
    x_pos = last_state.position[1]
    y_pos = last_state.position[2]
    z_pos = last_state.position[3]
    
    # Render marker
    scatter!(lscene, [x_pos], [y_pos], [z_pos],
        color=:yellow,
        markersize=8)
    
    return nothing
end
