# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

"""
    render_equatorial_plane!(lscene::LScene)

Render equatorial reference plane (xy plane at z=0) with two-level grid.

# Arguments
- `lscene::LScene` - Makie scene to render into

# V1.0 Rendering
- Green grid lines (Earth equatorial plane)
- Two-level grid: fine (10,000 km) and coarse (50,000 km)
- Extends ±1,250,000 km (beyond typical orbits)
- Transparent lines (alpha 0.15 for fine, 0.22-0.32 for coarse)
"""
function render_equatorial_plane!(lscene::LScene)
    # Grid parameters
    plane_size = 1_250_000.0  # km, extends beyond typical orbits
    fine_spacing = 10_000.0    # km between fine grid lines
    coarse_spacing = 50_000.0  # km between coarse grid lines (5x fine)
    
    n_fine = Int(plane_size / fine_spacing)
    n_coarse = Int(plane_size / coarse_spacing)
    
    # Draw fine grid lines (skip positions where coarse lines will be drawn)
    for i in -n_fine:n_fine
        # Skip if this position will have a coarse line
        if i * fine_spacing % coarse_spacing == 0
            continue
        end
        
        # Horizontal line
        y_grid = i * fine_spacing
        x_line = [-plane_size, plane_size]
        y_line = [y_grid, y_grid]
        z_line = [0.0, 0.0]
        lines!(lscene, x_line, y_line, z_line,
            color=RGBAf(0, 1, 0, 0.15),  # Green, more transparent
            linewidth=0.4)
        
        # Vertical line
        x_grid = i * fine_spacing
        x_line = [x_grid, x_grid]
        y_line = [-plane_size, plane_size]
        z_line = [0.0, 0.0]
        lines!(lscene, x_line, y_line, z_line,
            color=RGBAf(0, 1, 0, 0.15),
            linewidth=0.4)
    end
    
    # Draw coarse grid lines (brighter)
    for i in -n_coarse:n_coarse
        # Horizontal line
        y_grid = i * coarse_spacing
        x_line = [-plane_size, plane_size]
        y_line = [y_grid, y_grid]
        z_line = [0.0, 0.0]
        lines!(lscene, x_line, y_line, z_line,
            color=RGBAf(0, 1, 0, 0.32),  # Brighter green
            linewidth=0.5)
        
        # Vertical line
        x_grid = i * coarse_spacing
        x_line = [x_grid, x_grid]
        y_line = [-plane_size, plane_size]
        z_line = [0.0, 0.0]
        lines!(lscene, x_line, y_line, z_line,
            color=RGBAf(0, 1, 0, 0.22),
            linewidth=0.5)
    end
    
    return nothing
end
