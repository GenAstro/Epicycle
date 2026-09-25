# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# A ground station pass drawn the way an operator expects to see it:
# north at the top, bearings increasing clockwise, zenith at the centre.
#
# ── The trap in this example ──────────────────────────────────────────────────
# Plotly has `thetaunit`, so the ANGLE can stay in radians and it converts. There is
# no equivalent for the radial axis. Elevation therefore has to be converted, and if
# you forget, the plot does not fail — every point lands within a degree or so of the
# rim and you get a plausible-looking ring instead of a pass. That is what the first
# version of this example did.
#
# Epicycle computes in radians. Plotly draws in axis units. The conversion belongs to
# whoever knows which is which, and that is you.
#
# Review question: should there be an `runit` equivalent, given Plotly has none? It
# would be an invention, and §1 says we do not make those — but this is the one place
# the asymmetry actually bites.

using EpicycleIO
include(joinpath(@__DIR__, "synthetic.jl"))

az, el = fake_pass()            # both radians, as Epicycle computes them

# ── The pass ──────────────────────────────────────────────────────────────────
# `thetaunit` is a Plotly attribute and handles the angle. The elevation is converted
# here, because nothing will do it for us.

scatterpolar("Sky view", az, rad2deg.(el);
             thetaunit = "radians",
             mode      = "lines",
             name      = "Pass 3")

# ── The convention ────────────────────────────────────────────────────────────
# Two of these are correctness, not taste. Mathematics puts zero to the right and
# increases counterclockwise; a compass bearing puts zero at the top and increases
# clockwise. An azimuth drawn the mathematical way is MIRRORED, not merely rotated,
# and it looks entirely plausible while being wrong.
#
# The third says elevation runs 90 degrees at the centre to 0 at the rim, which falls
# out of the range being an ordered pair.

panel!("Sky view";
       polar_angularaxis_direction = "clockwise",
       polar_angularaxis_rotation  = 90,
       polar_radialaxis_range      = [90, 0],
       polar_radialaxis_title      = "elevation, deg")

# ── A second pass on the same sky ─────────────────────────────────────────────

az2, el2 = fake_pass(90)
scatterpolar!("Sky view", az2 .+ 1.9, rad2deg.(el2 .* 0.6);
              thetaunit = "radians",
              mode      = "lines",
              name      = "Pass 4",
              line_dash = "dash")

# ── The horizon mask ──────────────────────────────────────────────────────────
# A station will not track below some elevation. Drawn as a ring of constant
# elevation all the way round — in degrees, like everything else on this axis.

θ = range(0, 2π; length = 181)
scatterpolar!("Sky view", collect(θ), fill(10.0, length(θ));
              thetaunit  = "radians",
              mode       = "lines",
              name       = "10 deg mask",
              line_color = "gray",
              line_dash  = "dot")
