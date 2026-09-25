# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# The standard orbit determination picture: residuals inside a covariance band.
#
# This is the example that justifies `band` existing at all. Plotly draws a filled
# region as one trace whose x runs forward then backward, with `fill = "toself"`.
# Nobody should have to know that.
#
# Review question: is `band(x, lo, hi)` the right shape, or would you rather pass
# a centre and a half-width — `band(x, centre, ±3σ)`?

using EpicycleIO
include(joinpath(@__DIR__, "synthetic.jl"))

t, resid, σ = fake_residuals()

# ── The band goes down first, so the points draw on top of it ─────────────────

band("Range residuals", t, -3 .* σ, 3 .* σ;
     name      = "3σ",
     fillcolor = "rgba(120,120,120,0.25)",
     line_width = 0)

xyplot!("Range residuals", t, resid;
      mode        = "markers",
      name        = "residual",
      marker_size = 4)

panel!("Range residuals";
       title       = "Range residuals with 3σ covariance",
       xaxis_title = "hours",
       yaxis_title = "km")

# ── A second measurement type on its own panel ────────────────────────────────

t2, resid2, σ2 = fake_residuals(180)

band("Range-rate residuals", t2, -3 .* σ2 ./ 40, 3 .* σ2 ./ 40;
     name      = "3σ",
     fillcolor = "rgba(120,120,120,0.25)",
     line_width = 0)

xyplot!("Range-rate residuals", t2, resid2 ./ 40;
      mode        = "markers",
      name        = "residual",
      marker_size = 4)

panel!("Range-rate residuals"; xaxis_title = "hours", yaxis_title = "km/s")

# ── Colouring points by how many sigma they are ───────────────────────────────
# `marker_color` accepts one value per data point, which is Plotly's own meaning
# and passes straight through. Useful for spotting the ones about to be edited out.

nσ = abs.(resid) ./ σ

xyplot("Residual outliers", t, resid;
     mode         = "markers",
     name         = "residual",
     marker_size  = 5,
     marker_color = nσ,
     marker_colorscale = "Viridis",
     marker_colorbar_title = "sigma")
