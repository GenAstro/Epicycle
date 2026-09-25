# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# A thrust-direction figure for the Earth to Apophis rendezvous.
#
# The example is the data. This file includes it, reads the solved phase, and draws the transfer
# in the ecliptic plane with an arrow on each segment along its thrust direction, scaled by
# throttle magnitude. Segments the optimizer chose to coast have no arrow, which is the point of
# the picture: the solution has structure rather than a constant burn.
#
# Nothing here writes to the example.

include(joinpath(@__DIR__, "..", "examples", "Ex_ApophisRendezvousSimsFlanagan.jl"))

using EpicycleIO

const PANEL = "Earth to Apophis, low thrust"
const COAST_THRESHOLD = 0.01        # throttle magnitude below which a segment reads as coasting

# Read the solved trajectory off the phase. States are stored at segment boundaries, forward from
# Earth and backward from Apophis, meeting at the match point.
fwd_states = phase._states_fwd                      # 6 x (N/2 + 1), Earth to the match point
bwd_states = phase._states_bwd                      # 6 x (N/2 + 1), Apophis to the match point
u_fwd = forward_control(phase)
u_bwd = backward_control(phase)

# Assemble the whole arc: Earth, through the match point, to Apophis
arc = hcat(fwd_states, bwd_states[:, end:-1:1])
arc_x = arc[1, :] ./ AU
arc_y = arc[2, :] ./ AU

# Place each impulse at the midpoint of the segment it acts on, and pair it with its throttle
function segment_arrows(states, u)
    n = size(u, 2)
    x = [0.5 * (states[1, k] + states[1, k + 1]) / AU for k in 1:n]
    y = [0.5 * (states[2, k] + states[2, k + 1]) / AU for k in 1:n]
    mag = sqrt.(vec(sum(abs2, u, dims = 1)))
    return x, y, u, mag
end

fx, fy, fu, fmag = segment_arrows(fwd_states, u_fwd)
bx, by, bu, bmag = segment_arrows(bwd_states, u_bwd)

seg_x = vcat(fx, bx)
seg_y = vcat(fy, by)
seg_u = hcat(fu, bu)
seg_mag = vcat(fmag, bmag)

# Scale the longest arrow to a readable fraction of the figure
const ARROW_SPAN = 0.05
span = max(maximum(arc_x) - minimum(arc_x), maximum(arc_y) - minimum(arc_y))
scale = ARROW_SPAN * span / maximum(seg_mag)

# Build the arrows as one trace, separating each by a gap
function arrow_trace(x, y, u, scale)
    ax = Float64[]
    ay = Float64[]
    for k in eachindex(x)
        push!(ax, x[k], x[k] + scale * u[1, k], NaN)
        push!(ay, y[k], y[k] + scale * u[2, k], NaN)
    end
    return ax, ay
end

thrusting = seg_mag .>= COAST_THRESHOLD
arrow_x, arrow_y = arrow_trace(seg_x[thrusting], seg_y[thrusting],
                               seg_u[:, thrusting], scale)

# Sample each body's orbit over one revolution for context
function orbit_track(ephemeris, period)
    ts = range(0.0, period; length = 400)
    r = [ephemeris(t)[1] for t in ts]
    return [p[1] / AU for p in r], [p[2] / AU for p in r]
end

earth_x, earth_y = orbit_track(earth_ephemeris, 2π * sqrt(EARTH_A^3 / MU_SUN))
apophis_x, apophis_y = orbit_track(apophis_ephemeris, 2π * sqrt(APOPHIS_A^3 / MU_SUN))

r_depart, _, _ = earth_ephemeris(0.0)
r_arrive, _, _ = apophis_ephemeris(TOF)

# Draw the figure, back to front
xyplot(PANEL, earth_x, earth_y;
     name = "Earth orbit", mode = "lines",
     line_color = "rgb(150,160,175)", line_width = 1)

xyplot!(PANEL, apophis_x, apophis_y;
      name = "Apophis orbit", mode = "lines",
      line_color = "rgb(190,140,110)", line_width = 1)

xyplot!(PANEL, arc_x, arc_y;
      name = "transfer", mode = "lines",
      line_color = "rgb(90,150,255)", line_width = 3)

xyplot!(PANEL, arrow_x, arrow_y;
      name = "thrust direction", mode = "lines",
      line_color = "rgb(255,125,55)", line_width = 2)

xyplot!(PANEL, seg_x[.!thrusting], seg_y[.!thrusting];
      name = "coasting", mode = "markers",
      marker_size = 9, marker_color = "black",
      marker_line_color = "rgb(90,150,255)", marker_line_width = 2)

xyplot!(PANEL, [0.0], [0.0];
      name = "Sun", mode = "markers",
      marker_size = 14, marker_color = "rgb(245,190,60)")

xyplot!(PANEL, [r_depart[1] / AU, r_arrive[1] / AU], [r_depart[2] / AU, r_arrive[2] / AU];
      name = "departure and arrival", mode = "markers+text",
      text = ["Earth departure", "Apophis rendezvous"], textposition = ["bottom center", "bottom right"],
      textfont_size = 21,
      marker_size = 10, marker_color = "rgb(225,230,240)")

# Style for print on a dark page: light text on black,
# and gridlines dark enough to sit behind the data rather than compete with it
panel!(PANEL;
       font_size = 14,
       font_color = "rgb(232,236,244)",
       paper_bgcolor = "black",
       plot_bgcolor = "black",
       xaxis_title = "x  AU",
       yaxis_title = "y  AU",
       xaxis_gridcolor = "rgb(42,48,60)",
       yaxis_gridcolor = "rgb(42,48,60)",
       xaxis_zerolinecolor = "rgb(66,74,90)",
       yaxis_zerolinecolor = "rgb(66,74,90)",
       xaxis_linecolor = "rgb(232,236,244)",
       yaxis_linecolor = "rgb(232,236,244)",
       yaxis_scaleanchor = "x",
       yaxis_scaleratio = 1.0,
       showlegend = false)

# Report what the figure shows
@printf("\nsegments thrusting : %d of %d\n", count(thrusting), length(seg_mag))
@printf("segments coasting  : %d\n", count(.!thrusting))
@printf("throttle magnitude : min %.4f, max %.4f\n", minimum(seg_mag), maximum(seg_mag))
@printf("propellant         : %.3f kg of %.1f\n", M0 - final_mass(phase), M0)
