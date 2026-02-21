# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

"""
    render_trajectory!(lscene::LScene, view::View3D, sc::Spacecraft)

Render spacecraft trajectory from history with different colors per segment.
Optionally renders solver iterations if show_iterations=true.

# Arguments
- `lscene::LScene` - Makie scene to render into
- `view::View3D` - View containing options
- `sc::Spacecraft` - Spacecraft with history

# V1.0 Rendering
- Cycles through distinct colors per segment
- Linewidth 1.5 (fixed)
- Each history segment rendered separately with unique color
- Iterations (if enabled): gray, linewidth 0.5, alpha 0.3
"""
function render_trajectory!(lscene::LScene, view::View3D, sc::Spacecraft)
    # Render iterations first (background layer) if requested
    show_iters = get(get(view.options, sc, Dict()), :show_iterations, false)
    
    # Warn if user requested iterations but none exist
    if show_iters && isempty(sc.history.iterations)
        @warn "show_iterations=true for spacecraft '$(sc.name)', but spacecraft history does not contain iterations. " *
              "Use solve_trajectory!(...; record_iterations=true) to capture solver iterations."
    end
    
    if show_iters && !isempty(sc.history.iterations)
        for iteration in sc.history.iterations
            x_iter = Float64[]
            y_iter = Float64[]
            z_iter = Float64[]
            
            for state in iteration.states
                push!(x_iter, state.position[1])
                push!(y_iter, state.position[2])
                push!(z_iter, state.position[3])
            end
            
            if !isempty(x_iter)
                lines!(lscene, x_iter, y_iter, z_iter,
                    color=(:gray, 0.3),
                    linewidth=0.5)
            end
        end
    end
    
    # Color palette - cycle through these colors for segments
    # Optimized for dark background (black space theme)
    colors = [:cyan, :magenta, :yellow, :green, :orange, :red, :lightblue, :pink]
    
    # Render each segment with a different color
    for (seg_idx, segment) in enumerate(sc.history)
        # Extract trajectory for this segment
        x_traj = Float64[]
        y_traj = Float64[]
        z_traj = Float64[]
        
        for state in segment.states
            push!(x_traj, state.position[1])
            push!(y_traj, state.position[2])
            push!(z_traj, state.position[3])
        end
        
        # Skip empty segments
        if isempty(x_traj)
            continue
        end
        
        # Cycle through colors
        color = colors[mod1(seg_idx, length(colors))]
        
        # Render this segment
        lines!(lscene, x_traj, y_traj, z_traj,
            color=color,
            linewidth=1.5)
    end
    
    return nothing
end
