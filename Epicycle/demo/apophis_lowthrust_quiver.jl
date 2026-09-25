# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# The Earth to Apophis rendezvous flown with a tenth of the thrust.
#
# The same spacecraft and the same target as Ex_ApophisRendezvousSimsFlanagan, with a 0.1 N engine
# instead of 1 N and 1200 days instead of 360. The example is fixed at the one newton case, so this
# script carries its own setup rather than editing it.
#
# The trade is the point: a tenth the thrust delivers more mass, because the slower trajectory is
# the cheaper one. It also stops coasting — with time to spare a propellant-optimal solution
# throttles down everywhere rather than shutting off anywhere.

using Epicycle
using EpicycleIO
using LinearAlgebra
using Printf

const PANEL = "Earth to Apophis, 0.1 N"

const MU_SUN = 1.32712440018e11           # km^3/s^2
const AU = 1.495978707e8                  # km

const M0 = 1500.0                         # kg
const ISP = 3000.0                        # s
const TMAX = 1.0e-4                       # kN, a tenth of a newton
const G0 = 9.80665e-3                     # km/s^2

const TOF = 1200.0 * 86400.0              # s
const N_SEGMENTS = 60
const SMOOTHING = 0.015

# Earth and Apophis, with mean anomaly given at departure
const EARTH_A, EARTH_E, EARTH_I = 1.00000011 * AU, 0.01671, 0.0
const EARTH_RAAN, EARTH_AOP, EARTH_M = 0.0, deg2rad(102.9372), deg2rad(147.5124)

const APOPHIS_A, APOPHIS_E, APOPHIS_I = 0.9224 * AU, 0.1914, deg2rad(3.3393)
const APOPHIS_RAAN, APOPHIS_AOP, APOPHIS_M = deg2rad(203.9609), deg2rad(126.7213), deg2rad(178.8933)

# Create an ephemeris that advances mean anomaly and converts to a Cartesian state
function keplerian_ephemeris(a, e, i, raan, aop, m_departure)
    mean_motion = sqrt(MU_SUN / a^3)
    return function (t::Real)
        true_anomaly = mean_to_true_anomaly(m_departure + mean_motion * t, e)
        state = CartesianState(KeplerianState(a, e, i, raan, aop, true_anomaly), MU_SUN)
        r = collect(state.position)
        v = collect(state.velocity)
        return r, v, -(MU_SUN / norm(r)^3) .* r
    end
end

earth_ephemeris = keplerian_ephemeris(EARTH_A, EARTH_E, EARTH_I,
                                      EARTH_RAAN, EARTH_AOP, EARTH_M)
apophis_ephemeris = keplerian_ephemeris(APOPHIS_A, APOPHIS_E, APOPHIS_I,
                                        APOPHIS_RAAN, APOPHIS_AOP, APOPHIS_M)

_, v_departure, _ = earth_ephemeris(0.0)
_, v_arrival, _ = apophis_ephemeris(TOF)

# Create the phase
phase = SimsFlanaganPhase(name = :apophis_low_thrust,
                          transcription = SimsFlanagan(n_segments = N_SEGMENTS,
                                                       throttle_smoothing = SMOOTHING),
                          model = PropulsionModel(mu = MU_SUN,
                                                  Isp = ISP,
                                                  Tmax = TMAX,
                                                  g0 = G0),
                          ephemeris_left = earth_ephemeris,
                          ephemeris_right = apophis_ephemeris,
                          tspan = (0.0, TOF),
                          matchpoint_scale = [AU, AU, AU,
                                              norm(v_departure), norm(v_departure),
                                              norm(v_departure), M0])

# Declare the throttle blocks and the masses as optimization variables
Vary(forward_control, phase;
     guess = 0.5 * v_departure / norm(v_departure),
     lower_bound = fill(-2.0, 3),
     upper_bound = fill(2.0, 3))

Vary(backward_control, phase;
     guess = 0.5 * v_arrival / norm(v_arrival),
     lower_bound = fill(-2.0, 3),
     upper_bound = fill(2.0, 3))

Vary(initial_mass, phase; lower_bound = M0, upper_bound = M0)

Vary(final_mass, phase;
     guess = 0.80 * M0,
     lower_bound = 100.0,
     upper_bound = M0,
     scale = M0)

# Hold the throttle inside the unit ball at every segment
thrust_ball(c) = [dot(control(c), control(c))]
Constraint(thrust_ball, phase; lower_bound = 0.0, upper_bound = 1.0, at = Path())

# Minimize the propellant, which is the segment throttle magnitudes summed and scaled to kilograms
const KG_PER_THROTTLE = TMAX * (TOF / N_SEGMENTS) / (G0 * ISP)

throttle_magnitudes(u) = sqrt.(vec(sum(abs2, u, dims = 1)) .+ SMOOTHING^2)

propellant(c) = KG_PER_THROTTLE * (sum(throttle_magnitudes(forward_control(c))) +
                                   sum(throttle_magnitudes(backward_control(c))))

@partial(propellant, forward_control) do c
    u = forward_control(c)
    KG_PER_THROTTLE .* vec(u ./ throttle_magnitudes(u)')
end

@partial(propellant, backward_control) do c
    u = backward_control(c)
    KG_PER_THROTTLE .* vec(u ./ throttle_magnitudes(u)')
end

Objective(propellant, phase; sense = Min())

# Solve
result = solve!(Sequence(phase); method = Optimize(max_iter = 1000, print_level = 0))

# Read the solved trajectory off the phase
fwd_states = phase._states_fwd
bwd_states = phase._states_bwd
u_fwd = forward_control(phase)
u_bwd = backward_control(phase)

arc = hcat(fwd_states, bwd_states[:, end:-1:1])
arc_x = arc[1, :] ./ AU
arc_y = arc[2, :] ./ AU

# Place each impulse at the midpoint of the segment it acts on
function segment_arrows(states, u)
    n = size(u, 2)
    x = [0.5 * (states[1, k] + states[1, k + 1]) / AU for k in 1:n]
    y = [0.5 * (states[2, k] + states[2, k + 1]) / AU for k in 1:n]
    return x, y, sqrt.(vec(sum(abs2, u, dims = 1)))
end

fx, fy, fmag = segment_arrows(fwd_states, u_fwd)
bx, by, bmag = segment_arrows(bwd_states, u_bwd)

seg_x = vcat(fx, bx)
seg_y = vcat(fy, by)
seg_u = hcat(u_fwd, u_bwd)
seg_mag = vcat(fmag, bmag)

# Scale the longest arrow to a readable fraction of the figure
const ARROW_SPAN = 0.04
span = max(maximum(arc_x) - minimum(arc_x), maximum(arc_y) - minimum(arc_y))
scale = ARROW_SPAN * span / maximum(seg_mag)

arrow_x = Float64[]
arrow_y = Float64[]
for k in eachindex(seg_x)
    push!(arrow_x, seg_x[k], seg_x[k] + scale * seg_u[1, k], NaN)
    push!(arrow_y, seg_y[k], seg_y[k] + scale * seg_u[2, k], NaN)
end

# Sample each body's orbit for context
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
     line_color = "rgb(120,132,150)", line_width = 1, line_dash = "dot")

xyplot!(PANEL, apophis_x, apophis_y;
      name = "Apophis orbit", mode = "lines",
      line_color = "rgb(175,120,90)", line_width = 1, line_dash = "dot")

xyplot!(PANEL, arc_x, arc_y;
      name = "transfer", mode = "lines",
      line_color = "rgb(30,90,200)", line_width = 2)

xyplot!(PANEL, arrow_x, arrow_y;
      name = "thrust direction", mode = "lines",
      line_color = "rgb(230,90,30)", line_width = 2)

xyplot!(PANEL, [0.0], [0.0];
      name = "Sun", mode = "markers",
      marker_size = 14, marker_color = "rgb(245,190,60)")

xyplot!(PANEL, [r_depart[1] / AU, r_arrive[1] / AU], [r_depart[2] / AU, r_arrive[2] / AU];
      name = "departure and arrival", mode = "markers+text",
      text = ["Earth departure", "Apophis rendezvous"], textposition = "bottom center",
      marker_size = 10, marker_color = "rgb(40,40,45)")

# Style for print: near-black text on white, gridlines light enough to sit behind the data
panel!(PANEL;
       font_size = 14,
       font_color = "rgb(20,20,25)",
       paper_bgcolor = "white",
       plot_bgcolor = "white",
       xaxis_title = "x  AU",
       yaxis_title = "y  AU",
       xaxis_gridcolor = "rgb(226,230,236)",
       yaxis_gridcolor = "rgb(226,230,236)",
       xaxis_linecolor = "rgb(20,20,25)",
       yaxis_linecolor = "rgb(20,20,25)",
       yaxis_scaleanchor = "x",
       yaxis_scaleratio = 1.0,
       legend_orientation = "h")

# Report the case
@printf("\nstatus           : %s\n", result.info)
@printf("thrust           : %.1f N over %.0f days\n", TMAX * 1000, TOF / 86400)
@printf("delivered mass   : %.3f kg\n", final_mass(phase))
@printf("propellant       : %.3f kg\n", M0 - final_mass(phase))
@printf("throttle         : min %.4f, max %.4f\n", minimum(seg_mag), maximum(seg_mag))
@printf("coasting         : %d of %d segments\n", count(<(0.01), seg_mag), N_SEGMENTS)
