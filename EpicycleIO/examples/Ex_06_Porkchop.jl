# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# A porkchop plot. Departure date against arrival date, coloured by launch energy.
#
# This one is here to prove the call shape can carry a grid — `(x, y, Z)` — rather
# than only paired series. If it could not, adding contours later would be a
# redesign instead of an addition.
#
# Review question: is `Z` indexed [arrival, departure] the orientation you expect,
# or should the first index follow the first axis argument?

using EpicycleIO
include(joinpath(@__DIR__, "synthetic.jl"))

dep, arr, C3 = fake_porkchop()          # C3[i, j] is arrival arr[i], departure dep[j]

contour("Porkchop", dep, arr, C3;
        colorscale     = "Viridis",
        contours_start = 10,
        contours_end   = 60,
        contours_size  = 5,
        contours_showlabels = true,
        colorbar_title = "C3 km²/s²")

panel!("Porkchop";
       title       = "Earth to Mars, 2031 opportunity",
       xaxis_title = "departure, days from 2031-01-01",
       yaxis_title = "arrival, days from 2031-01-01")

# ── The same grid as a heatmap ────────────────────────────────────────────────
# Same data, same call shape, different trace name.

heatmap("Porkchop, filled", dep, arr, C3;
        colorscale     = "Viridis",
        colorbar_title = "C3 km²/s²")

# ── Marking a selected trajectory ─────────────────────────────────────────────
# A single point drawn over the contours. Cartesian trace on a cartesian panel,
# so this is allowed; a `scatterpolar` here would be an error naming both kinds.

xyplot!("Porkchop", [55.0], [300.0];
      mode          = "markers",
      name          = "selected",
      marker_size   = 12,
      marker_symbol = "x",
      marker_color  = "white")
