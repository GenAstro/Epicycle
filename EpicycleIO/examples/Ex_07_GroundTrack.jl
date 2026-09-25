# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# A ground track on a map, and the stations it passes over.
#
# Review question: Plotly wants degrees. Epicycle computes in radians. Here the
# conversion is explicit with `rad2deg`, whereas the polar example pushed it into
# Plotly's own `thetaunit` attribute. Is explicit conversion right here, or should
# there be an equivalent for geographic traces?

using EpicycleIO
include(joinpath(@__DIR__, "synthetic.jl"))

lon, lat = fake_groundtrack()

scattergeo("Ground track", lon, lat;
           mode       = "lines",
           name       = "Sat A",
           line_color = "cyan",
           line_width = 2)

panel!("Ground track";
       geo_projection_type = "equirectangular",
       geo_showcoastlines  = true,
       geo_coastlinecolor  = "gray",
       geo_showland        = true,
       geo_landcolor       = "rgb(30,30,30)")

# ── Ground stations ───────────────────────────────────────────────────────────
# Added to the same geographic panel.

station_lon = [-77.0, 11.9, 148.98]
station_lat = [38.9, 57.0, -35.4]
station_name = ["Wallops", "Kiruna", "Canberra"]

scattergeo!("Ground track", station_lon, station_lat;
            mode         = "markers+text",
            name         = "stations",
            text         = station_name,
            textposition = "top center",
            marker_size  = 9,
            marker_color = "orange")

# ── A polar view of the same track ────────────────────────────────────────────
# Different projection, so the high-latitude passes stop being misleading.

scattergeo("Ground track, polar", lon, lat;
           mode       = "lines",
           name       = "Sat A",
           line_color = "cyan")

panel!("Ground track, polar";
       geo_projection_type     = "orthographic",
       geo_projection_rotation_lon = 0,
       geo_projection_rotation_lat = 60,
       geo_showcoastlines      = true)
