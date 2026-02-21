# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

"""
    render_stars!(lscene::LScene)

Render star field background.

# Arguments
- `lscene::LScene` - Makie scene to render into

# V1.0 Rendering
- Generated star field (reproducible with seed 42)
- White stars with subtle glow
- 5000 stars on large sphere (10,000,000 km radius)
"""
function render_stars!(lscene::LScene)
    # V1.0: Generate star field
    star_radius = 100_000_000_000.0  # km - very large radius
    n_stars = 5000
    
    Random.seed!(42)  # Reproducible positions
    θ_stars = rand(n_stars) .* 2π
    φ_stars = acos.(2 .* rand(n_stars) .- 1)  # Uniform distribution on sphere
    
    x_stars = star_radius .* sin.(φ_stars) .* cos.(θ_stars)
    y_stars = star_radius .* sin.(φ_stars) .* sin.(θ_stars)
    z_stars = star_radius .* cos.(φ_stars)
    
    # Render stars
    scatter!(lscene, x_stars, y_stars, z_stars,
        color=:white,
        markersize=1.5,
        glowwidth=0.3,
        glowcolor=:white)
    
    return nothing
end
