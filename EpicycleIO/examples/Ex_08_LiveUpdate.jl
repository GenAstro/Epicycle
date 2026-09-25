# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Watching a solver converge while it runs.
#
# This is the example the whole publishing design exists for. `xyplot!` is safe to
# call from inside a callback: writes to the same panel coalesce, the file write is
# atomic so the browser never sees half a figure, and the call returns before
# anything is drawn so the solver is never waiting on a picture.
#
# Review questions:
#   - Is appending one point per call the right pattern, or should there be a
#     `push!`-style call that does not resend the whole series?
#   - Panels are named strings precisely so a callback has nothing to hold on to.
#     Does that read as a simplification or as a missing handle?

using EpicycleIO
include(joinpath(@__DIR__, "synthetic.jl"))

iters, cost, viol = fake_convergence()

# ── What a solver callback would do ───────────────────────────────────────────
# Stands in for `solve_trajectory!(seq, options; on_iteration = ...)`.

function on_iteration(k, c, v)
    xyplot!("Convergence", [k], [c];
          mode         = "markers",
          name         = "cost",
          marker_color = "cyan")

    xyplot!("Feasibility", [k], [v];
          mode         = "markers",
          name         = "max violation",
          marker_color = "orange")
end

# Start both panels empty so a re-run does not accumulate.
clear!("Convergence")
clear!("Feasibility")

for k in eachindex(iters)
    on_iteration(iters[k], cost[k], viol[k])
end

panel!("Convergence";  xaxis_title = "iteration", yaxis_title = "cost")
panel!("Feasibility";  xaxis_title = "iteration", yaxis_title = "max violation",
                       yaxis_type  = "log")

# ── The same thing by resending the series ────────────────────────────────────
# Cheaper on the browser, and what you would do if you already hold the history.
# Both work; the first is what a callback can do without keeping state.

xyplot("Convergence, resent", iters, cost;
     mode = "lines+markers", name = "cost")
panel!("Convergence, resent"; xaxis_title = "iteration", yaxis_title = "cost")
