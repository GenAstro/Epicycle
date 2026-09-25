# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# A position column is a vector of 3-element vectors. It should draw as three lines
# without you reshaping anything.
#
# Review questions:
#   - Is `position[1]` the right automatic name, or would you rather it were something else?
#   - Is the nesting rule for per-point attributes at the bottom acceptable, or too clever?

using EpicycleIO
include(joinpath(@__DIR__, "synthetic.jl"))

t, _, r = fake_orbit()          # `r` is what `Calc(position_vector, sat, frame)` returns

# ── Three lines from one call ─────────────────────────────────────────────────
# A scalar name is suffixed, because three legend entries all saying "position"
# would be useless.

xyplot(t, r; name = "position")
# legend reads: position[1], position[2], position[3]

# ── Name them yourself ────────────────────────────────────────────────────────

xyplot(t, r; name = ["x", "y", "z"])

# ── Style each component ──────────────────────────────────────────────────────
# One value per line. This is the ordinary case, not an advanced one.

xyplot(t, r;
     name       = ["x", "y", "z"],
     line_color = ["red", "green", "blue"],
     line_width = [1, 2, 3])

# ── A scalar applies to all three ─────────────────────────────────────────────

xyplot(t, r; name = ["x", "y", "z"], line_width = 1.5, line_dash = "dot")

# ── Where a vector means something else ───────────────────────────────────────
# Plotly already uses arrays on some attributes to mean "one value per data point".
# `marker_size` is one of those, so a flat vector keeps Plotly's meaning and
# nesting is how you say per-line instead.

n = length(t)

xyplot(t, r;
     mode        = "markers",
     marker_size = 3 .+ 3 .* sin.(range(0, 4π; length = n)))   # per data point

xyplot(t, r;
     mode        = "markers",
     marker_size = [[3], [6], [9]])                            # per line

# ── Independent pairs in one call ─────────────────────────────────────────────
# Nothing requires the pairs to relate to each other.

t_b, alt_b = fake_orbit_b()
_, alt, _  = fake_orbit()

xyplot(t, alt, t_b, alt_b; name = ["Sat A", "Sat B"])
