# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

"""
    render_spacecraft_label!(lscene::LScene, sc::Spacecraft)

Render text label for spacecraft showing its name.

# Arguments
- `lscene::LScene` - Makie scene to render into
- `sc::Spacecraft` - Spacecraft to label

# V1.0 Rendering
- Label shows spacecraft.name
- Positioned at current spacecraft state location
- White text, size 12
- Always visible
"""
function render_spacecraft_label!(lscene::LScene, sc::Spacecraft)
    # Get current spacecraft position (last point in history)
    last_segment = sc.history[end]
    last_state = last_segment.states[end]
    
    x_pos = last_state.position[1]
    y_pos = last_state.position[2]
    z_pos = last_state.position[3]
    
    # Render label at spacecraft position
    text!(lscene,
        sc.name,
        position=Point3f(x_pos, y_pos, z_pos),
        color=:white,
        fontsize=12)
    
    return nothing
end
