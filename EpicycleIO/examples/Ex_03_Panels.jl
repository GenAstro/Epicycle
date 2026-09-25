# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Panels are named with a leading string. They all live in one browser window.
#
# Review questions:
#   - Is a leading positional string the right way to name a panel, or would you
#     rather it were a keyword?
#   - Is "re-running the script replaces the panel" the behaviour you expect?

using EpicycleIO
include(joinpath(@__DIR__, "synthetic.jl"))

t, alt, r = fake_orbit()
t_b, alt_b = fake_orbit_b()

# ── Two panels, side by side in one window ────────────────────────────────────

xyplot("Altitude", t, alt; name = "Sat A")
xyplot("Position", t, r;   name = ["x", "y", "z"])

# ── Adding to a named panel ───────────────────────────────────────────────────

xyplot!("Altitude", t_b, alt_b; name = "Sat B")

# ── Replacing one ─────────────────────────────────────────────────────────────
# Same name, so this wipes the panel rather than adding a third line. Run this
# file twice and you still have two panels, not four.

xyplot("Altitude", t, alt; name = "Sat A, second look")

# ── Layout ────────────────────────────────────────────────────────────────────
# Keywords go to Plotly's layout, so Plotly's documentation applies directly.

panel!("Altitude";
       title       = "Altitude above the ellipsoid",
       xaxis_title = "hours from epoch",
       yaxis_title = "km")

panel!("Position";
       xaxis_title = "hours from epoch",
       yaxis_title = "km",
       legend_orientation = "h")

# ── A log axis, and a reversed one ────────────────────────────────────────────
# Limits are an ordered pair, so reversing an axis is just giving them backwards:
# time runs right to left here, which is how a decay margin is often read.

xyplot("Decay", t, alt .- minimum(alt) .+ 1e-3; name = "margin")
panel!("Decay"; yaxis_type = "log", xaxis_range = [6, 0],
                xaxis_title = "hours from epoch", yaxis_title = "km above minimum")

# ── Emptying without removing ─────────────────────────────────────────────────
# `clear!` takes the traces out and leaves the panel, layout and all. The panel goes
# blank for a moment, then the next plot fills it again — the axes below are still
# log and still reversed, because clearing data is not the same as forgetting how the
# panel is meant to look.

clear!("Decay")

xyplot("Decay", t, alt .- minimum(alt) .+ 1e-3 .+ 20; name = "margin, second run")
