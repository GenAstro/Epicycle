# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Anything EpicycleIO does not wrap is still reachable.
#
# There are wrappers for about ten trace types. Plotly has around forty, and adds
# more. Rather than chase that list, `xyplot` accepts a PlotlyBase trace directly, so
# every trace type is available whether or not anyone wrote a convenience for it.
#
# This is what makes the interface future-proof without a list to maintain: a
# wrapper is a shortcut for a common case, never a gate.
#
# Review question: this is the one place a user's script names PlotlyBase. Should
# EpicycleIO re-export the trace constructors so the dependency stays invisible and
# could be swapped later without breaking scripts?

using EpicycleIO
using PlotlyBase          # for trace types with no wrapper

include(joinpath(@__DIR__, "synthetic.jl"))

t, resid, σ = fake_residuals()

# ── A violin plot, which has no wrapper ───────────────────────────────────────
# Comparing how two measurement types are behaving.
#
# The residuals are NORMALISED — each divided by its own sigma — because range is in km
# and range-rate in km/s, and raw values in different units share no axis. Normalised,
# both should look like a unit normal centred on zero. Range does. Range-rate is fat and
# offset, which is what a wrong noise model or an unmodelled bias looks like.
#
# Plotly positions a violin numerically when only `y` is given, so the tick would read 0
# and 1. The legend already names them, so the ticks are turned off rather than faked with
# a category column — an earlier version passed one and the plots stopped drawing entirely.

range_n, rangerate_n = fake_normalised_residuals()

xyplot("Residual distribution",
     violin(y = range_n, box_visible = true, name = "range"))

xyplot!("Residual distribution",
      violin(y = rangerate_n, box_visible = true, name = "range-rate"))

panel!("Residual distribution";
       yaxis_title = "residual / sigma",
       xaxis_showticklabels = false)

# ── A box plot, likewise ──────────────────────────────────────────────────────
# Insertion altitude over 500 Monte Carlo runs. A box plot answers "how spread out is
# the outcome", so it wants a set of outcomes — not a quantity sampled over time, which
# is what an earlier version of this example wrongly handed it.

mc = fake_monte_carlo()

xyplot("Monte Carlo dispersion",
     box(y = mc, name = "insertion altitude", boxmean = "sd"))

panel!("Monte Carlo dispersion";
       yaxis_title = "km",
       xaxis_showticklabels = false,
       showlegend = true)

# ── Mixing a raw trace with a wrapped one ─────────────────────────────────────
# Both end up as traces on the same cartesian panel, so they compose.

xyplot("Mixed", t, resid; mode = "markers", name = "residual", marker_size = 4)
xyplot!("Mixed", scatter(x = t, y = 3 .* σ, mode = "lines", name = "3σ",
                       line_dash = "dash", line_color = "gray"))

# ── Full control of the layout ────────────────────────────────────────────────
# `panel!` takes layout keywords, but a Layout object works where that is easier.

panel!("Mixed", Layout(title = "Residuals against covariance",
                       xaxis_title = "hours",
                       yaxis_title = "km"))
