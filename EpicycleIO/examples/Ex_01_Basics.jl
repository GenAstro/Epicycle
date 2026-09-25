# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# The common case: draw a time history, style it, add a second one.
#
# Review question: is `xyplot(t, alt)` the call you want to type after a propagation?

using EpicycleIO
include(joinpath(@__DIR__, "synthetic.jl"))

t, alt, _ = fake_orbit()

# ── The whole thing ───────────────────────────────────────────────────────────
# Opens the dashboard if it is not open, creates the default panel, draws a line.

xyplot(t, alt)

# ── With attributes ───────────────────────────────────────────────────────────
# These are Plotly's names. Any other Plotly scatter attribute works the same way,
# because the value is handed over untouched.

xyplot(t, alt;
     name       = "altitude",
     mode       = "lines",
     line_color = "cyan",
     line_width = 2)

# ── Markers instead of a line ─────────────────────────────────────────────────

xyplot(t[1:20:end], alt[1:20:end];
     mode        = "markers",
     name        = "sampled",
     marker_size = 6)

# ── A second spacecraft, sampled differently ──────────────────────────────────
# `xyplot` replaces the panel, `xyplot!` adds to it. The two series have different
# lengths and different spans, which needs no special handling — they are just
# two pairs of arrays.

t_b, alt_b = fake_orbit_b()

xyplot(t,   alt;   name = "Sat A")
xyplot!(t_b, alt_b; name = "Sat B")

# ── Index on the horizontal axis ──────────────────────────────────────────────
# One argument means "plot this against its own index".

xyplot(alt; name = "altitude by sample")
