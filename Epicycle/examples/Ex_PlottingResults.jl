# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # Plotting Results
#'
#' Plot a propagated orbit, draw each plot type provided by EpicycleIO, and prepare a panel for
#' print. The first panel uses a six-hour propagation; the remaining panels use generated data to
#' keep the example focused on plotting.
#'
#' Plots are drawn by Plotly in a browser page served from the local machine, and the page updates
#' while a script is still running, so a long propagation or a solver iteration can be watched as it
#' proceeds. Attribute names are Plotly's own and pass through unchanged, so Plotly's reference
#' documentation applies directly.
#'
#' Only the first panel is real. The rest draw from the generators below, each named for the plot
#' type it feeds, because a porkchop needs a Lambert sweep and a sky view needs a station pass, and
#' generating either here would teach a second interface to explain the first.

using Epicycle
using EpicycleIO

#' ## Configuration
#'
#' Configure a spacecraft and propagator for the orbit panel.

# Configure the propagator
forces = ForceModel(PointMassGravity(earth, ()))
integ = IntegratorConfig(DP8();
                        dt = 60.0,
                        reltol = 1e-12,
                        abstol = 1e-12)
propagator = OrbitPropagator(forces, integ)

# Place the spacecraft in low Earth orbit
start_epoch = Time("2024-03-01T00:00:00.000", UTC(), ISOT())
sat = Spacecraft(state = CartesianState([6878.137, 0.0, 0.0, 0.0, 4.71754, 5.99820]),
                 time = start_epoch,
                 name = "Explorer")

#' ## Plot a propagated orbit
#'
#' Propagate the orbit and plot altitude and speed from the recorded history.

# Propagate for six hours
propagate!(propagator, sat, StopAt(sat, PropDurationSeconds(), 6 * 3600.0))

# Read the recorded trajectory back as quantities rather than as raw state vectors
t, r_mag, v_mag = history(Calc(epoch, sat),
                          Calc(position_magnitude, sat),
                          Calc(velocity_magnitude, sat))

hours = [24 * (ti - t[1]) for ti in t]
altitude = r_mag .- earth.equatorial_radius

# Plot altitude and speed against elapsed time
xyplot("Orbit", hours, altitude;
     name = "altitude  km",
     line_color = "rgb(30,90,200)",
     line_width = 3)

xyplot!("Orbit", hours, v_mag .* 100;
      name = "speed  cm/s × 10⁵",
      line_color = "rgb(230,90,30)",
      line_width = 2)

panel!("Orbit";
       xaxis_title = "hours from epoch",
       legend_orientation = "h")

#' ## Data for the rest of the page
#'
#' One generator per plot shape, each named for the plot type it feeds. `grid_data` serves the
#' contour, the heatmap and the surface, because all three take two axes and a matrix.

# Two axes and a matrix, indexed [row, column] with rows following y and columns following x
grid_data() = (collect(range(0, 120; length = 60)),
               collect(range(180, 400; length = 70)),
               [12 + 0.004(d - 55)^2 + 0.0025(a - 300)^2 + 6sin(d / 23) * cos(a / 37)
                for a in range(180, 400; length = 70), d in range(0, 120; length = 60)])

# A series and the one-sigma envelope around it
xy_data(n = 240) = (collect(range(0, 12; length = n)),
                    [0.004 + 0.010 * abs(sin(2π * x / 9)) for x in range(0, 12; length = n)])

# Azimuth and elevation through one station pass, in radians
polar_data(n = 120) = (deg2rad.(35 .+ 210 .* range(0, 1; length = n)),
                       deg2rad.(78 .* sin.(π .* range(0, 1; length = n)) .+ 2))

# Values whose distribution is the point
histogram_data(n = 600) = [0.014 * sin(11k) + 0.009 * cos(3k) for k in 1:n]

# Points in three axes
scatter_data(n = 250) = ([90 * sin(0.7k) + 12 * cos(5k) for k in 1:n],
                         [70 * cos(0.7k) + 12 * sin(5k) for k in 1:n],
                         [0.4 * sin(0.3k) for k in 1:n])

# Longitude and latitude in degrees
ground_track_data(n = 500) = (collect(range(-180, 180; length = n)),
                              [42 * sin(deg2rad(2.2x)) for x in range(-180, 180; length = n)])

# A value per category
bar_data() = (["TOI", "MCC-1", "MCC-2", "MOI"], [2.82, 0.14, 0.06, 1.48])

#' ## A panel for each plot type
#'
#' A contour for a porkchop, a heatmap for the same grid filled, a surface for it in three axes, a
#' band for a covariance envelope, a polar plot for a station pass, a bar chart for a delta-v
#' budget, a histogram for a residual distribution, a geographic plot for a ground track, and a 3D
#' scatter for a B-plane.

# Draw a porkchop plot with labeled contours
departure, arrival, c3 = grid_data()

contour("Porkchop", departure, arrival, c3;
        colorscale = "Viridis",
        contours_start = 10,
        contours_end = 60,
        contours_size = 5,
        contours_showlabels = true,
        colorbar_title = "C3  km²/s²")

# Draw the grid as a heatmap
heatmap("Porkchop, filled", departure, arrival, c3;
        colorscale = "Viridis",
        colorbar_title = "C3  km²/s²")

panel!("Porkchop, filled";
       xaxis_title = "departure, days",
       yaxis_title = "arrival, days")

# Draw the grid as a surface
surface("Porkchop, in relief", departure, arrival, c3;
        colorscale = "Plasma",
        colorbar_title = "C3  km²/s²")

panel!("Porkchop, in relief";
       scene_xaxis_title = "departure, days",
       scene_yaxis_title = "arrival, days",
       scene_zaxis_title = "C3  km²/s²")

# Overlay residuals on a covariance envelope
t_resid, sigma = xy_data()
residual = [0.9 * s * sin(7x) for (x, s) in zip(t_resid, sigma)]

band("Range residuals", t_resid, -3 .* sigma, 3 .* sigma;
     name = "3σ",
     fillcolor = "rgba(120,130,150,0.25)",
     line_width = 0)

xyplot!("Range residuals", t_resid, residual;
      mode = "markers",
      name = "residual",
      marker_size = 4,
      marker_color = "rgb(30,90,200)")

panel!("Range residuals";
       title = "Range residuals with 3σ covariance",
       xaxis_title = "hours",
       yaxis_title = "km")

# Draw a sky view with north at the top and zenith at the center
az, el = polar_data()

scatterpolar("Sky view", az, rad2deg.(el);
             thetaunit = "radians",
             mode = "lines",
             name = "Pass 3",
             line_color = "rgb(230,90,30)",
             line_width = 3)

panel!("Sky view";
       polar_angularaxis_direction = "clockwise",
       polar_angularaxis_rotation = 90,
       polar_radialaxis_range = [90, 0],
       polar_radialaxis_title = "elevation, deg")

# Draw a delta-v budget, one bar per maneuver
maneuvers, delta_vs = bar_data()

bar("Delta-v budget", maneuvers, delta_vs;
    marker_color = "rgb(30,90,200)")

panel!("Delta-v budget"; yaxis_title = "delta-v  km/s")

# Draw the residual distribution
histogram("Residual distribution", histogram_data();
          nbinsx = 40,
          marker_color = "rgb(120,140,190)")

panel!("Residual distribution";
       xaxis_title = "postfit residual  km",
       yaxis_title = "count")

# Draw a ground track and tracking stations
track_lon, track_lat = ground_track_data()

scattergeo("Ground track", track_lon, track_lat;
           mode = "lines",
           name = "Explorer",
           line_color = "rgb(0,170,200)",
           line_width = 2)

scattergeo!("Ground track", [-77.0, 11.9, 148.98], [38.9, 57.0, -35.4];
            mode = "markers+text",
            name = "stations",
            text = ["Wallops", "Kiruna", "Canberra"],
            textposition = "top center",
            marker_size = 9,
            marker_color = "rgb(255,170,40)")

panel!("Ground track";
       geo_projection_type = "equirectangular",
       geo_showcoastlines = true,
       geo_showland = true)

# Draw one B-plane point per Monte Carlo case
b_dot_r, b_dot_t, arrival_offset = scatter_data()

scatter3d("B-plane", b_dot_r, b_dot_t, arrival_offset;
          mode = "markers",
          name = "dispersion",
          marker_size = 3,
          marker_color = arrival_offset,
          marker_colorscale = "Turbo")

panel!("B-plane";
       scene_xaxis_title = "B·R  km",
       scene_yaxis_title = "B·T  km",
       scene_zaxis_title = "arrival offset  hours")

#' ## Adding to a panel, and replacing one
#'
#' A leading string names a panel. Plotting into a name that already exists replaces what is in it,
#' so running a script twice leaves one figure rather than two, and the `!` form adds instead.

# Add a second series to a panel that already holds one
xyplot!("Delta-v budget", maneuvers, delta_vs .* 1.08;
      name = "with margin",
      mode = "markers",
      marker_size = 10,
      marker_color = "rgb(230,90,30)")

#' ## A trace with no wrapper
#'
#' Plotly has around forty trace types and nine have wrappers here. The rest are reached by building
#' the trace with `PlotlyBase`, which EpicycleIO re-exports, and handing it to `xyplot`. A wrapper is
#' a convenience for a common case rather than a gate.

# Build an unwrapped Plotly violin trace
xyplot("Residual spread",
     PlotlyBase.violin(y = histogram_data(),
                       box_visible = true,
                       name = "range"))

panel!("Residual spread";
       yaxis_title = "residual  km",
       xaxis_showticklabels = false)

#' ## Styling a panel for print
#'
#' Everything above takes the dashboard's default appearance, which is made for a screen. A figure
#' headed for a document or a slide wants near-black text on white, gridlines light enough to sit
#' behind the data, and a font size that survives being shrunk into a column.
#'
#' Metadata and appearance are both `panel!` keywords, and settings accumulate, so both blocks below
#' apply to a panel that is already drawn without redrawing it.

# Set the title, axis labels, and legend placement
panel!("Porkchop";
       title = "Earth to Mars, 2031 opportunity",
       xaxis_title = "departure, days from 2031-01-01",
       yaxis_title = "arrival, days from 2031-01-01",
       legend_orientation = "h")

# Use print-friendly colors, gridlines, and type
panel!("Porkchop";
       font_size = 14,
       font_color = "rgb(20,20,25)",
       paper_bgcolor = "white",
       plot_bgcolor = "white",
       xaxis_gridcolor = "rgb(226,230,236)",
       yaxis_gridcolor = "rgb(226,230,236)",
       xaxis_linecolor = "rgb(20,20,25)",
       yaxis_linecolor = "rgb(20,20,25)")

# Report what was drawn
println("panels drawn : 11")
println("real data    : the Orbit panel, from a six hour propagation")
println("invented data: every other panel, to show the call rather than the analysis")
