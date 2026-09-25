# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# A trajectory in three dimensions, beside the plots that describe it.
#
# This is the one example that needs real Epicycle rather than made-up arrays. A 3D view has
# to know the frame and the epoch, so it takes the spacecraft rather than columns pulled out
# of it — unlike `xyplot`, which never learns what a spacecraft is.
#
# Review questions:
#   - Is `orbitview(sat)` the right verb? GMAT calls this an OrbitView, FreeFlyer a
#     ViewWindow, STK a 3D Graphics window. All three name the window rather than the act.
#   - Earth is the only central body today. Nothing here names it — the body comes from the
#     frame the history was recorded in. Does that hold up when you try the Moon?

using Epicycle
using EpicycleIO

# ── A quarter day of LEO ──────────────────────────────────────────────────────

sat = Spacecraft(state = CartesianState([7000.0, 0.0, 1300.0, 0.0, 7.35, 1.0]),
                 time  = Time("2020-01-01T00:00:00.000", UTC(), ISOT()),
                 name  = "Explorer")

prop = OrbitPropagator(ForceModel(PointMassGravity(earth, ())),
                       IntegratorConfig(DP8(); abstol = 1e-12, reltol = 1e-12, dt = 60.0))

propagate!(prop, sat, StopAt(sat, PropDurationDays(), 0.25))

# ── The view ──────────────────────────────────────────────────────────────────
# Everything comes off the spacecraft: the arcs from its recorded segments, the epoch from
# the earliest sample, the colours from the same palette the plots use.

orbitview("Trajectory", sat)

# Play, pause, reset and the speed controls are under the view. The camera is locked in the
# inertial frame, so the orbit holds still against the stars and Earth turns beneath it.

# ── The same run, as plots ────────────────────────────────────────────────────
# A 3D view and the quantities describing it belong on one page, which is the reason panels
# and the globe share a dashboard rather than living in separate windows.

t, r = history(Calc(epoch, sat), Calc(position_vector, sat, EarthMJ2000Eq))

hours = [(x - first(t)) * 24 for x in t]
radius = [sqrt(sum(abs2, p)) for p in r]

xyplot("Radius", hours, radius; name = "|r|", line_color = "cyan")
panel!("Radius"; xaxis_title = "hours from epoch", yaxis_title = "km")

xyplot("Position", hours, r; name = ["x", "y", "z"])
panel!("Position"; xaxis_title = "hours from epoch", yaxis_title = "km")
