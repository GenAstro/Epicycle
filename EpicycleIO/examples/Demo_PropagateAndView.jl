# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Propagate, burn, propagate — then look at it three ways on one page.
#
#   include(joinpath(pkgdir(EpicycleIO), "examples", "Demo_PropagateAndView.jl"))
#
# The burn splits the history into two segments, which the 3D view draws as two coloured
# arcs. The plots below it are the same run read as quantities.

using Epicycle
using EpicycleIO

EpicycleIO.clear_all!()          # start from an empty dashboard

# ── A spacecraft, and something to propagate it with ──────────────────────────

sat = Spacecraft(state = CartesianState([7000.0, 0.0, 1300.0, 0.0, 7.35, 1.0]),
                 time  = Time("2020-01-01T00:00:00.000", UTC(), ISOT()),
                 name  = "Explorer")

prop = OrbitPropagator(ForceModel(PointMassGravity(earth, ())),
                       IntegratorConfig(DP8(); abstol = 1e-12, reltol = 1e-12, dt = 60.0))

# ── Coast, burn, coast ────────────────────────────────────────────────────────
# Two propagation segments, so the view has two arcs to colour differently.

propagate!(prop, sat, StopAt(sat, PropDurationDays(), 0.08))

maneuver!(sat, ImpulsiveManeuver(axes = VNB(), element1 = 0.20,
                                 element2 = 0.0, element3 = 0.0))

propagate!(prop, sat, StopAt(sat, PropDurationDays(), 0.25))

# ── The trajectory ────────────────────────────────────────────────────────────
# Takes the spacecraft, because a trajectory on a globe needs the frame and the epoch.
# Play, pause, reset and the speed controls are under the view.

orbitview("Trajectory", sat)

# ── The same run, as quantities ───────────────────────────────────────────────
# One walk of the history, several columns out of it.

t, r, rmag, vmag, sma = history(Calc(epoch,              sat),
                                Calc(position_vector,    sat, EarthMJ2000Eq),
                                Calc(position_magnitude, sat),
                                Calc(velocity_magnitude, sat),
                                Calc(semi_major_axis,    sat))

hours = [(x - first(t)) * 24 for x in t]        # Time is not a number; render it first

xyplot("Radius and speed", hours, rmag; name = "|r|  km", line_color = "cyan")
xyplot!("Radius and speed", hours, vmag .* 1000; name = "|v|  m/s", line_color = "orange")
panel!("Radius and speed"; xaxis_title = "hours from epoch")

# A position column arrives as a vector of 3-vectors and draws as three lines.
xyplot("Position", hours, r; name = ["x", "y", "z"])
panel!("Position"; xaxis_title = "hours from epoch", yaxis_title = "km")

# The burn is visible here as a step, which is the point of plotting it beside the view.
xyplot("Semi-major axis", hours, sma; name = "a", line_color = "rgb(55,255,55)")
panel!("Semi-major axis"; xaxis_title = "hours from epoch", yaxis_title = "km")

# ── And the same numbers as a file ────────────────────────────────────────────

report(joinpath(tempdir(), "explorer.txt");
       hours = hours, radius = rmag, speed = vmag, position = r)

println("dashboard : ", EpicycleIO.dashboard_url())
println("report    : ", joinpath(tempdir(), "explorer.txt"))
