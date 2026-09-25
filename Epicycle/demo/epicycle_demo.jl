# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# A guided tour of Epicycle, narrated in the terminal and on the dashboard.
#
#   include("Epicycle/demo/epicycle_demo.jl")
#   EpicycleDemo.epicycle_demo()                      # timed pauses
#   EpicycleDemo.epicycle_demo(wait_for_enter = true) # presenter mode: Enter advances
#   EpicycleDemo.epicycle_demo(pace = 0.0)            # no pauses, for a quick check
#
# Seven parts: a spacecraft in orbit with its 3D view, one orbit read in several frames with its
# ground track, a targeted Hohmann transfer, a GEO transfer written as a GMAT-style target block,
# orbit determination by batch least squares and by Kalman filter, a low-cost transfer between
# Earth-Moon L1 and L2 Lyapunov orbits, and a minimum-fuel lunar landing. Every number printed is
# computed during the run.
#
# The narration goes to the terminal and to a caption strip across the top of the dashboard, so an
# audience watching the browser reads the same story. The code the strip shows is read from this
# file at the step that runs it, so it cannot drift from what actually ran.

module EpicycleDemo

using LinearAlgebra
using Printf
using Random

using EpicycleBase
using AstroEpochs
using AstroStates
using AstroUniverse
using AstroFrames
using AstroModels
using AstroManeuvers
using AstroCallbacks
using AstroProp
using AstroSolve
using AstroRoutines
using EpicycleIO

# Vern9 comes from the solver package AstroProp loads, so the demo needs nothing further.
using AstroProp: Vern9

export epicycle_demo

# The narration helpers, the caption strip and the @step macro.
include(joinpath(@__DIR__, "narration.jl"))


# ═══════════════════════════════════════════════════════════════════════════════
# Part 1 — a spacecraft in orbit
# ═══════════════════════════════════════════════════════════════════════════════

function _act_orbit()
    _act(1, "A spacecraft in orbit",
         "State, epoch and frame in one object; a propagator; a trajectory you can fly through";
         panels = ["1 · Explorer in 3D", "1 · Altitude and speed"])

    _say("  A spacecraft carries its state, its epoch and the coordinate system they are in.")
    _say("  This one flies an orbit like the International Space Station's, from 1 January 2020.")
    _pause(1.5)

    @step "Building the spacecraft" begin
        kep = KeplerianState(6778.0, 0.0005, deg2rad(51.6), deg2rad(40.0), 0.0, 0.0)
        sat = Spacecraft(state = kep,
                         time  = Time("2020-01-01T00:00:00.000", UTC(), ISOT()),
                         name  = "Explorer")
    end
    _fact("semi-major axis", @sprintf("%.1f", kep.sma); unit = "km")
    _fact("eccentricity", @sprintf("%.4f", kep.ecc))
    _fact("inclination", @sprintf("%.2f", rad2deg(kep.inc)); unit = "deg")
    _fact("period", @sprintf("%.1f", 2π * sqrt(kep.sma^3 / earth.mu) / 60); unit = "min")
    _pause(2.0)

    println()
    _say("  The force model is Earth's gravity with the Moon and Sun as third bodies,")
    _say("  their positions read from the JPL DE440 ephemeris. Vern9 integrates it.")
    _pause(1.0)

    @step "Propagating a quarter of a day" begin
        forces = ForceModel(PointMassGravity(earth, (moon, sun)))
        integ  = IntegratorConfig(Vern9(); abstol = 1e-11, reltol = 1e-11, dt = 60.0)
        prop   = OrbitPropagator(forces, integ)
        propagate!(prop, sat, StopAt(sat, PropDurationDays(), 0.25))
    end

    @step "Drawing the trajectory in 3D" begin
        orbitview("1 · Explorer in 3D", sat)
    end

    @step "Plotting altitude and speed" begin
        t, rmag, vmag = history(Calc(epoch, sat), Calc(position_magnitude, sat),
                                Calc(velocity_magnitude, sat))
        altitude = rmag .- earth.equatorial_radius

        xyplot("1 · Altitude and speed", _hours(t), altitude;
             name = "altitude  km", line_color = "cyan", line_width = 3)
        xyplot!("1 · Altitude and speed", _hours(t), vmag .* 100;
              name = "speed  ×100 km/s", line_color = "orange", line_width = 2)
        panel!("1 · Altitude and speed"; xaxis_title = "hours from epoch")
    end

    println()
    _fact("samples recorded", string(length(t)))
    _fact("altitude range", @sprintf("%.1f to %.1f", minimum(altitude), maximum(altitude));
          unit = "km")
    println()
    _say("  The 3D view shows the orbit on a turning Earth. Press play under it; the orbit")
    _say("  holds still against the stars while the planet rotates beneath it.", color = _DIM)
    _pause(4.0)
    return sat
end

# ═══════════════════════════════════════════════════════════════════════════════
# Part 2 — one orbit, many frames
# ═══════════════════════════════════════════════════════════════════════════════

function _act_frames(sat)
    _act(2, "One orbit, many frames",
         "The same state read in inertial, ecliptic and Earth-fixed axes";
         panels = ["2 · Ground track"])

    _say("  A number like RAAN means nothing until you say which frame it was measured in.")
    _say("  Every quantity takes the coordinate system as an argument, so it is always said.")
    _pause(1.5)
    println()

    @step "Reading RAAN and inclination in three frames" begin
        frames = ("Earth mean equator J2000"  => EarthMJ2000Eq,
                  "Earth mean ecliptic J2000" => EarthMJ2000Ec,
                  "Earth-fixed ITRF"          => EarthFixed)
        angles = [name => (rad2deg(raan(sat, cs)), rad2deg(inclination(sat, cs)))
                  for (name, cs) in frames]
    end

    println()
    printstyled(@sprintf("    %-30s %14s %14s\n", "frame", "RAAN  deg", "incl  deg");
                color = _ACCENT, bold = true)
    for (name, (Ω, i)) in angles
        push!(_CAPTION.facts, name => @sprintf("RAAN %.3f°   incl %.3f°", Ω, i))
        _publish_caption()
        printstyled(@sprintf("    %-30s", name); color = :normal)
        printstyled(@sprintf(" %14.3f %14.3f\n", Ω, i); color = _NUMBER, bold = true)
        _wait(0.6)
    end
    println()
    _say("  Each value is correct in the frame it was measured in; the three differ because the",
         color = _DIM)
    _say("  reference planes differ.", color = _DIM)
    _pause(2.5)

    println()
    _say("  Reading the whole recorded trajectory in Earth-fixed axes gives the ground track.")
    @step "Converting every sample to ITRF and drawing the ground track" begin
        t, r_fixed = history(Calc(epoch, sat), Calc(position_vector, sat, EarthFixed))
        lat = [rad2deg(asin(p[3] / norm(p))) for p in r_fixed]
        lon = [rad2deg(atan(p[2], p[1])) for p in r_fixed]

        scattergeo("2 · Ground track", lon, lat;
                   mode = "lines", name = "Explorer", line_color = "cyan", line_width = 2)
        scattergeo!("2 · Ground track", [-116.89, 148.98, -4.25], [35.43, -35.40, 40.43];
                    mode = "markers+text", name = "Deep Space Network",
                    text = ["Goldstone", "Canberra", "Madrid"], textposition = "top center",
                    marker_size = 10, marker_color = "orange")
        panel!("2 · Ground track";
               geo_projection_type = "natural earth", geo_showcoastlines = true,
               geo_coastlinecolor = "gray", geo_showland = true,
               geo_landcolor = "rgb(35,40,48)", geo_showocean = true,
               geo_oceancolor = "rgb(12,18,30)", geo_bgcolor = "rgba(0,0,0,0)")
    end

    println()
    _fact("ground-track samples", string(length(lat)))
    _fact("latitude reached", @sprintf("±%.1f", maximum(abs, lat)); unit = "deg")
    _pause(3.5)
    return nothing
end

# ═══════════════════════════════════════════════════════════════════════════════
# Part 3 — targeting a Hohmann transfer
# ═══════════════════════════════════════════════════════════════════════════════

function _act_targeting()
    _act(3, "Targeting a transfer",
         "Two burns, one coast, two goals: the solver finds the burns";
         panels = ["3 · Hohmann transfer", "3 · Radius through the transfer"])

    _say("  The mission: raise a 7,000 km circular orbit to a circular orbit at 45,000 km.")
    _say("  The script declares what may vary and what must hold; the solver picks the values.")
    _pause(1.5)

    @step "Describing the mission" begin
        sat  = Spacecraft(state = KeplerianState(7000.0, 0.001, 0.0, 0.0, 0.0, 1.0),
                          time  = Time("2020-09-21T12:23:12", TAI(), ISOT()),
                          name  = "Transfer")
        prop = OrbitPropagator(ForceModel(PointMassGravity(earth, (moon, sun))),
                               IntegratorConfig(Vern9(); abstol = 1e-11, reltol = 1e-11, dt = 300.0))

        toi = ImpulsiveManeuver(axes = VNB(), element1 = 0.1)
        moi = ImpulsiveManeuver(axes = VNB(), element1 = 0.4)

        toi_event = Event(name = "Transfer orbit insertion",
                          event = () -> maneuver!(sat, toi),
                          vars = [SolverVariable(calc = ManeuverCalc(toi, sat, DeltaVVector()),
                                                 name = "toi", lower_bound = [0.0, 0.0, 0.0],
                                                 upper_bound = [2.5, 0.0, 0.0])],
                          funcs = [])
        coast_event = Event(name = "Coast to apoapsis",
                            event = () -> propagate!(prop, sat,
                                                     StopAt(sat, PosDotVel(), 0.0; direction = -1)))
        moi_event = Event(name = "Mission orbit insertion",
                          event = () -> maneuver!(sat, moi),
                          vars = [SolverVariable(calc = ManeuverCalc(moi, sat, DeltaVVector()),
                                                 name = "moi", lower_bound = [0.0, 0.0, 0.0],
                                                 upper_bound = [3.0, 0.0, 0.0])],
                          funcs = [Constraint(position_magnitude, sat; equals = 45000.0),
                                   Constraint(eccentricity, sat; equals = 0.0)])

        seq = Sequence()
        add_sequence!(seq, toi_event, coast_event, moi_event)
    end

    println()
    printstyled("    burn   →   coast to apoapsis   →   burn   ⇒   |r| = 45,000 km,  e = 0\n";
                color = _ACCENT, bold = true)
    _pause(2.0)

    result = @step "Solving with IPOPT, propagating the coast at every iteration" begin
        solve!(seq; method = Optimize(print_level = 0, derivatives = :fd,
                                      extra = Dict{String,Any}("sb" => "yes")))
    end

    dv1, dv2 = toi.element1, moi.element1
    r1, r2 = 7000.0, 45000.0
    h1 = sqrt(earth.mu / r1) * (sqrt(2r2 / (r1 + r2)) - 1)
    h2 = sqrt(earth.mu / r2) * (1 - sqrt(2r1 / (r1 + r2)))

    println()
    _fact("solver status", string(result.info))
    _fact("first burn", @sprintf("%.4f", dv1); unit = "km/s")
    _fact("second burn", @sprintf("%.4f", dv2); unit = "km/s")
    _fact("total ΔV", @sprintf("%.4f", dv1 + dv2); unit = "km/s")
    _fact("textbook two-body Hohmann", @sprintf("%.4f", h1 + h2); unit = "km/s")
    println()
    _say("  The textbook figure ignores the Moon and Sun; the solve includes them, so the two")
    _say("  agree closely without matching exactly.", color = _DIM)
    _pause(2.5)

    @step "Flying the final orbit for a day and drawing it" begin
        propagate!(prop, sat, StopAt(sat, PropDurationDays(), 1.0))
        orbitview("3 · Hohmann transfer", sat)

        t, rmag, sma = history(Calc(epoch, sat), Calc(position_magnitude, sat),
                               Calc(semi_major_axis, sat))
        xyplot("3 · Radius through the transfer", _hours(t), rmag;
             name = "|r|  km", line_color = "cyan", line_width = 3)
        xyplot!("3 · Radius through the transfer", _hours(t), sma;
              name = "semi-major axis  km", line_color = "rgb(55,255,55)", line_width = 2)
        panel!("3 · Radius through the transfer"; xaxis_title = "hours from first burn",
               yaxis_title = "km")
    end
    println()
    _say("  In the 3D view each burn is marked with the ΔV measured across it.", color = _DIM)
    _pause(3.5)
    return nothing
end

# ═══════════════════════════════════════════════════════════════════════════════
# Part 4 — a GEO transfer written as a target block
# ═══════════════════════════════════════════════════════════════════════════════

# The GEO transfer of AstroSolve's runtests_sequence_geotransfer.jl, whose initial state and targets
# come from a GMAT script. There it is assembled as an event graph; here it is written in flight
# order with `target!`.

function _act_target_block()
    _act(4, "A GMAT-style target block",
         "Three burns and four goals, written in the order the spacecraft flies them";
         panels = ["4 · GEO transfer", "4 · Radius and inclination"])

    _say("  GMAT users write a Target sequence: propagate, vary a burn, apply it, achieve a goal.")
    _say("  target! reads the same way. The block is recorded, then solved as one problem.")
    _pause(1.5)

    @step "Setting up the spacecraft, burns and stopping conditions" begin
        sat = Spacecraft(
            state = CartesianState([3737.792, -4607.692, -2845.644, 5.411, 5.367, -1.566]),
            time  = Time("2000-01-01T11:59:28.000", UTC(), ISOT()),
            name  = "GeoSat-1")
        prop = OrbitPropagator(ForceModel(PointMassGravity(earth, ())),
                               IntegratorConfig(Vern9(); abstol = 1e-12, reltol = 1e-12, dt = 60.0))

        toi = ImpulsiveManeuver(axes = VNB(), element1 = 2.518, element2 = 0.0,   element3 = 0.0)
        mcc = ImpulsiveManeuver(axes = VNB(), element1 = 0.559, element2 = 0.588, element3 = 0.0)
        moi = ImpulsiveManeuver(axes = VNB(), element1 = 0.282, element2 = 0.0,   element3 = 0.0)

        z_crossing = StopAt(position_z, sat, EarthMJ2000Eq; equals = 0.0)
        apoapsis   = StopAt(position_dot_velocity, sat; equals = 0.0, direction = -1)
        perigee    = StopAt(position_dot_velocity, sat; equals = 0.0, direction =  1)
    end
    _pause(2.0)

    println()
    _say("  Now the mission itself. Each Vary belongs to the burn after it, and each Constraint")
    _say("  checks the state the step before it produced.")
    _pause(1.5)

    result = @step "Solving the target block with IPOPT" begin
        target!(method = Optimize(derivatives = :fd, print_level = 0,
                                  extra = Dict{String,Any}("sb" => "yes"))) do
            propagate!(prop, sat, z_crossing)

            Vary(delta_v, toi; lower_bound = [0.0, 0.0, 0.0], upper_bound = [8.0, 0.0, 0.0])
            maneuver!(sat, toi)

            propagate!(prop, sat, apoapsis)
            Constraint(position_magnitude, sat; equals = 85000.0)

            propagate!(prop, sat, perigee)
            propagate!(prop, sat, z_crossing)

            Vary(delta_v, mcc; lower_bound = [-1.0, -1.0, -0.001], upper_bound = [4.0, 1.0, 0.001])
            maneuver!(sat, mcc)

            propagate!(prop, sat, perigee)
            Constraint(inclination, sat, EarthMJ2000Eq; equals = deg2rad(2.0))
            Constraint(position_magnitude, sat; equals = 42195.0)

            Vary(delta_v, moi; lower_bound = [-1.0, -0.001, -0.001], upper_bound = [4.0, 0.001, 0.001])
            maneuver!(sat, moi)
            Constraint(semi_major_axis, sat; equals = 42166.90)
        end
    end

    c = result.constraints
    println()
    _fact("solver status", string(result.info))
    _fact("transfer orbit insertion", @sprintf("%.4f", norm(delta_v(toi))); unit = "km/s")
    _fact("mid-course plane change", @sprintf("%.4f", norm(delta_v(mcc))); unit = "km/s")
    _fact("mission orbit insertion", @sprintf("%.4f", norm(delta_v(moi))); unit = "km/s")
    _fact("apoapsis radius", @sprintf("%.3f", c[1]); unit = "km")
    _fact("inclination at perigee", @sprintf("%.4f", rad2deg(c[2])); unit = "deg")
    _fact("perigee radius", @sprintf("%.3f", c[3]); unit = "km")
    _fact("final semi-major axis", @sprintf("%.3f", c[4]); unit = "km")
    println()
    _say("  AstroSolve's tests check this problem's event values against GMAT to 1e-7. Here it is")
    _say("  the same problem with no Event or Sequence assembled.", color = _DIM)
    _pause(2.5)

    @step "Flying a day in the new orbit and drawing the transfer" begin
        propagate!(prop, sat, StopAt(sat, PropDurationDays(), 1.0))
        orbitview("4 · GEO transfer", sat)

        t, rmag, inc = history(Calc(epoch, sat), Calc(position_magnitude, sat),
                               Calc(inclination, sat, EarthMJ2000Eq))
        xyplot("4 · Radius and inclination", _hours(t), rmag ./ 1000;
             name = "|r|  1000 km", line_color = "cyan", line_width = 3)
        xyplot!("4 · Radius and inclination", _hours(t), rad2deg.(inc);
              name = "inclination  deg", line_color = "orange", line_width = 2)
        panel!("4 · Radius and inclination"; xaxis_title = "hours from epoch")
    end
    println()
    _fact("epoch to final burn", @sprintf("%.1f", _hours(t)[end] - 24.0); unit = "h")
    _say("  The transfer climbs to 85,000 km, changes plane and raises perigee at an equator")
    _say("  crossing, and circularises at the 42,195 km perigee.", color = _DIM)
    _pause(3.5)
    return nothing
end

# ═══════════════════════════════════════════════════════════════════════════════
# Part 5 — orbit determination
# ═══════════════════════════════════════════════════════════════════════════════

# The scenario of the regression cases uc8, uc10 and uc11: a low Earth orbit tracked by the three
# Deep Space Network complexes with two-way range and Doppler.
const OD_TRUTH       = [4382.084365606881, -2755.139441815347, -5071.185336852633,
                        4.470494673424, 5.860946905163, 0.987952901499]    # km, km/s
const OD_GUESS_ERROR = [1.0, -1.0, 0.5, 0.005, -0.005, 0.002]            # km, km/s
const OD_PRIOR_SIGMA = [5.0, 5.0, 5.0, 0.05, 0.05, 0.05]                 # km, km/s
const OD_SIGMA_RANGE   = 15.0e-3      # km, 15 m
const OD_SIGMA_DOPPLER = 2.0e-5       # km/s, 2 cm/s
const OD_ARC_HOURS   = 6.0
const OD_CADENCE_S   = 60.0

const _STATION_COLORS = Dict("Goldstone" => "cyan", "Madrid" => "orange", "Canberra" => "rgb(55,255,55)")

function _act_estimation()
    _act(5, "Orbit determination",
         "Three DSN stations, six hours of tracking: batch least squares, a Kalman filter, and the filter in a loop";
         panels = ["5 · Tracking data", "5 · Batch residuals",
                   "5 · Filter and smoother", "5 · The filter, live"])

    _say("  Goldstone, Madrid and Canberra track a spacecraft in low Earth orbit for six hours,")
    _say("  measuring two-way range and Doppler. The truth is known here, so every estimate can")
    _say("  be scored against it.")
    _pause(1.5)

    @step "Simulating the tracking and writing it as a CCSDS tracking data file" begin
        t0 = Time("2010-06-10T11:40:00", TT(), ISOT())
        stations = [GroundStation(name = name, body = earth, latitude = lat, longitude = lon,
                                  altitude = alt, min_elevation = 5.0)
                    for (name, lat, lon, alt) in (("Goldstone",  35.4267, -116.89, 1.0014),
                                                  ("Madrid",     40.43,     -4.25, 0.8),
                                                  ("Canberra",  -35.40,    148.98, 0.7))]
        prop = OrbitPropagator(ForceModel(PointMassGravity(earth, ())),
                               IntegratorConfig(Vern9(); dt = 10.0, reltol = 1e-12, abstol = 1e-12))

        truth = Spacecraft(state = CartesianState(copy(OD_TRUTH)), time = t0,
                           coord_sys = CoordinateSystem(earth, ICRF()), name = "Sat")
        rng = MersenneTwister(2010)
        observations = ObservationRecord[]
        truth_at = Dict{Int, Vector{Float64}}()          # seconds from epoch => true state
        for k in 1:round(Int, OD_ARC_HOURS * 3600 / OD_CADENCE_S)
            propagate!(prop, truth, StopAt(truth, PropDurationSeconds(), OD_CADENCE_S))
            y = to_posvel(truth)
            truth_at[round(Int, k * OD_CADENCE_S)] = y
            t = Time(truth.time.jd1, truth.time.jd2, truth.time.scale, :jd)   # the file writes ISO
            for gs in stations
                AstroModels.is_visible(gs, y[1:3], truth.time) || continue
                r_gs, v_gs = get_state(gs, truth.time)
                ρ = y[1:3] .- r_gs
                push!(observations, ObservationRecord(:RANGE, t,
                    2norm(ρ) + OD_SIGMA_RANGE * randn(rng), gs.name))
                push!(observations, ObservationRecord(:DOPPLER, t,
                    2dot(ρ, y[4:6] .- v_gs) / norm(ρ) + OD_SIGMA_DOPPLER * randn(rng), gs.name))
            end
        end

        tdm = TrackingDataFile(joinpath(tempdir(), "epicycle_demo_dsn.tdm"))
        write_records(tdm, observations, TDMHeader("2010-06-10T00:00:00", "EPICYCLE DEMO"),
                      [TDMSegmentMeta(time_system = "TT", participant_1 = gs.name,
                                      participant_2 = "Sat") for gs in stations])
        records, _, _ = read_records(tdm)
    end

    _fact("observations", string(length(records)))
    _fact("range noise, 1σ", @sprintf("%.0f", OD_SIGMA_RANGE * 1e3); unit = "m")
    _fact("Doppler noise, 1σ", @sprintf("%.0f", OD_SIGMA_DOPPLER * 1e5); unit = "cm/s")
    _fact("tracking file", basename(tdm.path))

    hours_of(r) = ((r.t_receive.jd1 - t0.jd1) + (r.t_receive.jd2 - t0.jd2)) * 24
    for gs in stations
        rs = filter(r -> r.participant_1 == gs.name && r.measurement_type === :RANGE, records)
        isempty(rs) && continue
        f = gs.name == stations[1].name ? plot : plot!
        f("5 · Tracking data", hours_of.(rs), [r.observed / 2 for r in rs];
          mode = "markers", name = "$(gs.name) range", marker_size = 5,
          marker_color = _STATION_COLORS[gs.name])
    end
    panel!("5 · Tracking data"; xaxis_title = "hours", yaxis_title = "one-way range  km")
    _pause(2.5)

    println()
    _say("  First, batch least squares: every observation at once, from a guess a kilometre off.")
    @step "Estimating the epoch state by batch least squares" begin
        y_guess = OD_TRUTH .+ OD_GUESS_ERROR
        sc = Spacecraft(state = CartesianState(copy(y_guess)), time = t0,
                        coord_sys = CoordinateSystem(earth, ICRF()), name = "Sat")
        measurements = vcat(
            [TwoWayRange(SignalPath(gs, sc, gs); noise = MeasurementNoise(OD_SIGMA_RANGE))
             for gs in stations],
            [TwoWayDoppler(SignalPath(gs, sc, gs); noise = MeasurementNoise(OD_SIGMA_DOPPLER))
             for gs in stations])

        x0 = Vary(state, sc; guess = y_guess, covariance = Diagonal(OD_PRIOR_SIGMA .^ 2))
        problem = ODProblem(spacecraft = sc, propagator = prop, measurements = measurements,
                            solve_for = [x0])
        batch = solve!(problem, records; method = Batch(n_iters = 15, tol = 1e-9))
    end

    err = batch.X_hat .- OD_TRUTH
    _fact("iterations", string(batch.iters))
    _fact("position error", @sprintf("%.2f", norm(err[1:3]) * 1e3); unit = "m")
    _fact("position formal σ", @sprintf("%.2f", norm(batch.sigma[1:3]) * 1e3); unit = "m")
    _fact("velocity error", @sprintf("%.2f", norm(err[4:6]) * 1e6); unit = "mm/s")

    is_range = [r.measurement_type === :RANGE for r in records]
    res = [first(v) for v in batch.residuals]
    xyplot("5 · Batch residuals", hours_of.(records[is_range]), res[is_range] ./ OD_SIGMA_RANGE;
         mode = "markers", name = "range / σ", marker_size = 5, marker_color = "cyan")
    xyplot!("5 · Batch residuals", hours_of.(records[.!is_range]),
          res[.!is_range] ./ OD_SIGMA_DOPPLER;
          mode = "markers", name = "Doppler / σ", marker_size = 5, marker_color = "orange")
    panel!("5 · Batch residuals"; xaxis_title = "hours",
           yaxis_title = "post-fit residual in sigmas", yaxis_range = [-5, 5])
    rms(v) = sqrt(sum(abs2, v) / length(v))
    _fact("range residual RMS", @sprintf("%.2f", rms(res[is_range]) / OD_SIGMA_RANGE); unit = "σ")
    _fact("Doppler residual RMS", @sprintf("%.2f", rms(res[.!is_range]) / OD_SIGMA_DOPPLER); unit = "σ")
    println()
    _say("  An RMS near one sigma means the fit is consistent with the noise that was simulated.",
         color = _DIM)
    _pause(3.0)

    println()
    _say("  A filter processes the same data one observation at a time, with process noise so it")
    _say("  never becomes overconfident. The smoother then runs back over the arc.")
    @step "Kalman filter with a smoother, iterated over the arc" begin
        sc = Spacecraft(state = CartesianState(copy(y_guess)), time = t0,
                        coord_sys = CoordinateSystem(earth, ICRF()), name = "Sat")
        x0 = Vary(state, sc; guess = y_guess, covariance = OD_PRIOR_SIGMA .^ 2,
                  process_noise = DiagonalSNC([0.0, 0.0, 0.0, 1e-12, 1e-12, 1e-12]))
        problem = ODProblem(spacecraft = sc, propagator = prop, measurements = measurements,
                            solve_for = [x0])
        filtered = solve!(problem, records; method = Sequential(iterations = 3, smoother = RTS()))
    end

    ekf, rts = filtered.ekf, filtered.rts
    filter_hours = [r.t / 3600 for r in ekf.records]
    filter_err   = [norm(r.y_post[1:3] .- truth_at[round(Int, r.t)][1:3]) * 1e3 for r in ekf.records]
    smooth_err   = [norm(rts.y_smooth[k][1:3] .- truth_at[round(Int, ekf.records[k].t)][1:3]) * 1e3
                    for k in eachindex(ekf.records)]
    xyplot("5 · Filter and smoother", filter_hours, filter_err;
         name = "filter", line_color = "orange", line_width = 2)
    xyplot!("5 · Filter and smoother", filter_hours, smooth_err;
          name = "smoother", line_color = "cyan", line_width = 3)
    panel!("5 · Filter and smoother"; xaxis_title = "hours",
           yaxis_title = "position error  m", yaxis_type = "log")
    _fact("filter iterations", string(filtered.iters_run))
    _fact("smoothed position error, median",
          @sprintf("%.2f", sort(smooth_err)[cld(length(smooth_err), 2)]); unit = "m")
    _pause(3.0)

    println()
    _say("  The same filter, stepped by hand. Between any two updates a script can look at the")
    _say("  estimate, act on it, or change what comes next. The plot fills in as each update")
    _say("  lands: the 3σ bound, then the error against it.")
    @step "The filter in a loop, one observation at a time" begin
        sc = Spacecraft(state = CartesianState(copy(y_guess)), time = t0,
                        coord_sys = CoordinateSystem(earth, ICRF()), name = "Sat")
        x0 = Vary(state, sc; guess = y_guess, covariance = OD_PRIOR_SIGMA .^ 2,
                  process_noise = DiagonalSNC([0.0, 0.0, 0.0, 1e-12, 1e-12, 1e-12]))
        problem = ODProblem(spacecraft = sc, propagator = prop, measurements = measurements)
        dyn, meas, obs_times, obs_data, R = build_od_closures(records, problem)

        ekf = init_ekf([x0], dyn, meas; model = problem, R = R[1], t0 = 0.0)
        hours, error_m, sigma_m = Float64[], Float64[], Float64[]
        for k in eachindex(obs_times)
            time_update!(ekf, obs_times[k])
            measurement_update!(ekf, obs_data[k]; R = R[k])

            # Everything a controller would need is in hand here: the estimate, its covariance,
            # and the truth this demo happens to know.
            push!(hours,   obs_times[k] / 3600)
            push!(error_m, norm(current_state(ekf)[1:3] .- truth_at[round(Int, obs_times[k])][1:3]) * 1e3)
            push!(sigma_m, 3 * sqrt(sum(diag(current_covariance(ekf))[1:3])) * 1e3)

            xyplot("5 · The filter, live", hours, sigma_m;
                 name = "3σ", line_color = "cyan", line_width = 2)
            xyplot!("5 · The filter, live", hours, error_m;
                  name = "position error", mode = "markers", marker_color = "orange", marker_size = 6)
            _wait(0.03)
        end
    end
    panel!("5 · The filter, live"; xaxis_title = "hours", yaxis_title = "position  m",
           yaxis_type = "log")
    _fact("updates", string(length(obs_times)))
    _fact("final position error",
          @sprintf("%.2f", norm(current_state(ekf)[1:3] .- truth_at[round(Int, obs_times[end])][1:3]) * 1e3);
          unit = "m")
    _fact("final position 3σ",
          @sprintf("%.2f", 3 * sqrt(sum(diag(current_covariance(ekf))[1:3])) * 1e3); unit = "m")
    inside = count(error_m .<= sigma_m)
    _fact("updates with error inside 3σ", "$inside of $(length(error_m))")
    println()
    _say(inside == length(error_m) ?
         "  The error stays inside the 3σ bound at every update, so the covariance the filter reports is consistent with the error it makes." :
         "  Orange is the true error and cyan the 3σ the filter reports; the count above says how often the two agree.",
         color = _DIM)
    _pause(3.5)
    return nothing
end

# ═══════════════════════════════════════════════════════════════════════════════
# Part 6 — optimal control: an L1 to L2 Lyapunov transfer
# ═══════════════════════════════════════════════════════════════════════════════

# The regression case uc23. The Earth-Moon mass ratio and both Lyapunov orbits are JPL's (SSD
# periodic orbits API, family lyapunov, near C = 3.15), and the connection between equal-energy
# L1 and L2 orbits is Koon, Lo, Marsden and Ross, Chaos 10(2), 2000.

const LYAP_MU = 1.215058560962404e-02
const LYAP_L1 = [8.1596252146384562e-01, 0.0, 0.0, 0.0, 2.0722124749217649e-01, 0.0]
const LYAP_L2 = [1.1182825695532028e+00, 0.0, 0.0, 0.0, 1.8601928389638619e-01, 0.0]
const LYAP_L1_C, LYAP_L1_PERIOD = 3.15001683280912, 2.8447547942812985
const LYAP_L2_C, LYAP_L2_PERIOD = 3.15000013081292, 3.4205694071950448
const LYAP_DAYS_PER_TU = 4.342          # Earth-Moon time unit, days
const LYAP_KMS_PER_VU  = 1.02322        # Earth-Moon velocity unit, km/s

const LYAP_MODEL_SOURCE = raw"""
struct CRState{T}   <: AbstractState;   x::T; y::T; z::T; vx::T; vy::T; vz::T end
struct CRControl{T} <: AbstractControl; ux::T; uy::T; uz::T                 end
struct CRModel;     mu::Float64                                              end

# The ballistic three-body acceleration is AstroRoutines'; the control adds to it.
function cr_dynamics!(dy, y::CRState, u::CRControl, p, t, model)
    a = cr3bp_accel([y.x, y.y, y.z, y.vx, y.vy, y.vz], model.mu)
    dy[1], dy[2], dy[3] = y.vx, y.vy, y.vz
    dy[4], dy[5], dy[6] = a[1] + u.ux, a[2] + u.uy, a[3] + u.uz
end

@partial(cr_dynamics!, state) do dF, y, u, p, t, model
    dF .= cr3bp_jacobian([y.x, y.y, y.z, y.vx, y.vy, y.vz], model.mu)
end

@partial(cr_dynamics!, control) do dF, y, u, p, t, model
    dF[4, 1] = dF[5, 2] = dF[6, 3] = 1.0
end

endpoint(c) = [state(c).x, state(c).y, state(c).z, state(c).vx, state(c).vy, state(c).vz]
@partial(endpoint, state) do c
    Matrix(1.0I, 6, 6)
end

effort(c) = control(c).ux^2 + control(c).uy^2 + control(c).uz^2
@partial(effort, control) do c
    [2control(c).ux  2control(c).uy  2control(c).uz]
end
"""
include_string(@__MODULE__, LYAP_MODEL_SOURCE)

"""The x and y of a state flown ballistically for `tf` time units, for drawing an orbit."""
function _cr3bp_xy(s0, tf; n = 400)
    V    = AstroProp.OrdinaryDiffEqVerner
    prob = V.ODEProblem((ds, s, p, t) -> cr3bp_eom!(ds, s, LYAP_MU), collect(s0), (0.0, tf))
    sol  = V.solve(prob, Vern9(); abstol = 1e-13, reltol = 1e-13,
                   saveat = range(0.0, tf; length = n))
    return [u[1] for u in sol.u], [u[2] for u in sol.u]
end

function _act_lyapunov()
    _act(6, "Optimal control · between libration orbits",
         "Earth-Moon L1 to L2 at equal energy, where a nearly free path exists";
         panels = ["6 · L1 to L2 transfer", "6 · Control and energy"])

    _say("  Beyond the Moon's libration points sit families of periodic orbits. At the same")
    _say("  Jacobi energy, an L1 orbit and an L2 orbit are joined by their invariant manifolds,")
    _say("  so a transfer between them should cost almost nothing.")
    _pause(1.5)

    @step "Checking the model against JPL's three-body data" begin
        L1x = libration_point(LYAP_MU, :L1)[1]
        L2x = libration_point(LYAP_MU, :L2)[1]
        C1  = jacobi_constant(LYAP_L1, LYAP_MU)
        C2  = jacobi_constant(LYAP_L2, LYAP_MU)
    end
    println()
    _fact("L1 point, x", @sprintf("%.15f   JPL 0.836915125772357", L1x))
    _fact("L2 point, x", @sprintf("%.14f    JPL 1.15568216544488", L2x))
    _fact("Jacobi constant, L1 orbit", @sprintf("%.14f    JPL %s", C1, LYAP_L1_C))
    _fact("Jacobi constant, L2 orbit", @sprintf("%.14f    JPL %s", C2, LYAP_L2_C))
    _pause(2.5)

    _CAPTION.code, _CAPTION.status = strip(LYAP_MODEL_SOURCE), "the model a user writes"
    _say("  A user writes the state, the control, the dynamics and its partials. The three-body")
    _say("  acceleration and its Jacobian come from AstroRoutines.")
    _pause(3.5)

    @step "Posing the problem: pinned ends, free time, minimum control effort" begin
        tf_guess = (LYAP_L1_PERIOD + LYAP_L2_PERIOD) / 2
        phase = CollocationPhase(name = :lyapunov_transfer,
                                 transcription = HermiteSimpson(n_steps = 50),
                                 dynamics = cr_dynamics!, model = CRModel(LYAP_MU),
                                 state = CRState, control = CRControl, tspan = (0.0, tf_guess))
        Vary(state, phase; guess = hcat(LYAP_L1, LYAP_L2),
             lower_bound = [0.70, -0.45, -0.05, -3.0, -3.0, -0.5],
             upper_bound = [1.30,  0.45,  0.05,  3.0,  3.0,  0.5])
        Vary(control, phase; guess = zeros(3, 2),
             lower_bound = [-0.5, -0.5, -0.1], upper_bound = [0.5, 0.5, 0.1])
        Vary(final_time, phase; guess = tf_guess,
             lower_bound = 0.5tf_guess, upper_bound = 4.0tf_guess)
        Constraint(endpoint, phase; equals = LYAP_L1, at = Initial())
        Constraint(endpoint, phase; equals = LYAP_L2, at = Final())
        Objective(effort, phase; sense = Min(), at = Path())
    end

    result = @step "Transcribing and solving the nonlinear program" begin
        solve!(Sequence(phase); method = Optimize(print_level = 0,
                                                  extra = Dict{String,Any}("sb" => "yes")))
    end

    tf = get_final_time(phase)
    Y, U = state(phase), control(phase)
    ts = get_node_times(phase)
    Cs = [jacobi_constant(Y[:, k], LYAP_MU) for k in axes(Y, 2)]

    println()
    _fact("solver status", string(result.info))
    _fact("transfer time", @sprintf("%.3f TU, %.2f", tf, tf * LYAP_DAYS_PER_TU); unit = "days")
    _fact("control effort ∫|u|² dt", @sprintf("%.7e", result.objective))
    # ∫|u| dt ≤ sqrt(tf ∫|u|² dt), so this bounds the delta-V from above.
    _fact("delta-V, upper bound",
          @sprintf("%.1f", 1000 * sqrt(result.objective * tf) * LYAP_KMS_PER_VU); unit = "m/s")
    _fact("energy change ΔC", @sprintf("%.2e", Cs[end] - Cs[1]))
    println()
    _say("  A ballistic connection is asymptotic at both ends, so a finite arc pinned to the two")
    _say("  orbits cannot reach zero cost.", color = _DIM)
    _pause(2.5)

    @step "Drawing both orbits and the transfer in the rotating frame" begin
        l1x, l1y = _cr3bp_xy(LYAP_L1, LYAP_L1_PERIOD)
        l2x, l2y = _cr3bp_xy(LYAP_L2, LYAP_L2_PERIOD)
        frame = "6 · L1 to L2 transfer"
        xyplot(frame, l1x, l1y; name = "L1 Lyapunov orbit", line_color = "royalblue", line_width = 4)
        xyplot!(frame, l2x, l2y; name = "L2 Lyapunov orbit", line_color = "rgb(55,255,55)",
              line_width = 4)
        xyplot!(frame, Y[1, :], Y[2, :]; name = "transfer", line_color = "orange", line_width = 3)
        xyplot!(frame, [1 - LYAP_MU, L1x, L2x], [0.0, 0.0, 0.0]; mode = "markers+text",
              name = "Moon and libration points", text = ["Moon", "L1", "L2"],
              textposition = "bottom center", marker_size = [12, 8, 8],
              marker_color = ["gray", "white", "white"])
        panel!(frame; xaxis_title = "x  Earth-Moon distances", yaxis_title = "y",
               yaxis_scaleanchor = "x", yaxis_scaleratio = 1)

        umag = [norm(U[:, k]) for k in axes(U, 2)]
        xyplot("6 · Control and energy", ts, umag;
             name = "|u|  control acceleration", line_color = "orange", line_width = 3)
        xyplot!("6 · Control and energy", ts, (Cs .- LYAP_L2_C) .* 1e5;
              name = "(C − C of L2 orbit) × 1e5", line_color = "cyan", line_width = 2)
        panel!("6 · Control and energy"; xaxis_title = "time  TU")
    end
    println()
    _say("  The orange arc leaves the blue L1 orbit, passes the Moon, and joins the green L2 orbit.",
         color = _DIM)
    _pause(3.5)
    return nothing
end

# ═══════════════════════════════════════════════════════════════════════════════
# Part 7 — optimal control: landing on the Moon
# ═══════════════════════════════════════════════════════════════════════════════

const LANDER_MODEL_SOURCE = raw"""
struct LanderState{T}   <: AbstractState;   h::T; v::T; m::T end
struct LanderControl{T} <: AbstractControl; thrust::T      end
struct LanderModel;     mdot_coeff::Float64; thrust_max::Float64 end

function lander_dynamics!(dy, y::LanderState, u::LanderControl, p, t, model)
    dy[1] =  y.v                            # altitude rate
    dy[2] = -1.0 + u.thrust / y.m           # gravity against thrust
    dy[3] = -u.thrust / model.mdot_coeff    # propellant burned
end

@partial(lander_dynamics!, state) do dF, y, u, p, t, model
    dF[1, 2] =  1.0
    dF[2, 3] = -u.thrust / y.m^2
end

@partial(lander_dynamics!, control) do dF, y, u, p, t, model
    dF[2, 1] =  1.0 / y.m
    dF[3, 1] = -1.0 / model.mdot_coeff
end

lander_start(c)     = [state(c).h, state(c).v, state(c).m]
lander_touchdown(c) = [state(c).h, state(c).v]
lander_mass(c)      = state(c).m
"""
include_string(@__MODULE__, LANDER_MODEL_SOURCE)

function _act_landing()
    _act(7, "Optimal control · landing on the Moon",
         "Maximise the fuel left at touchdown; the answer is a bang-bang burn";
         panels = ["7 · Lunar descent", "7 · Thrust"])

    _say("  A lander descends at 0.78 units per unit time. It must touch down at rest and keep")
    _say("  as much propellant as it can. Meditch solved this in 1964; the published figures")
    _say("  are a remaining mass of 0.3953 and a landing time of 1.397.")
    _pause(1.5)

    _CAPTION.code, _CAPTION.status = strip(LANDER_MODEL_SOURCE), "the model a user writes"
    _say("  The model is three states, one control and three lines of dynamics.")
    _pause(3.5)

    @step "Posing the landing" begin
        model = LanderModel(2.349, 1.227)
        phase = CollocationPhase(name = :moon_landing,
                                 transcription = HermiteSimpson(n_steps = 30),
                                 dynamics = lander_dynamics!, model = model,
                                 state = LanderState, control = LanderControl, tspan = (0.0, 1.4))
        Vary(state, phase; guess = [1.0 0.0; -0.783 0.0; 1.0 0.4],
             lower_bound = [0.0, -5.0, 0.001], upper_bound = [5.0, 5.0, 2.0])
        Vary(control, phase; guess = reshape([model.thrust_max / 2, model.thrust_max / 2], 1, 2),
             lower_bound = [0.0], upper_bound = [model.thrust_max])
        Vary(final_time, phase; guess = 1.4, lower_bound = 0.5, upper_bound = 5.0)
        Constraint(lander_start,     phase; equals = [1.0, -0.783, 1.0], at = Initial())
        Constraint(lander_touchdown, phase; equals = [0.0, 0.0],         at = Final())
        Objective(lander_mass, phase; sense = Max())
    end

    result = @step "Solving for the thrust history" begin
        solve!(Sequence(phase); method = Optimize(print_level = 0,
                                                  extra = Dict{String,Any}("sb" => "yes")))
    end

    tf = get_final_time(phase)
    Y  = state(phase)
    U  = control(phase)
    tt = collect(range(0.0, tf; length = size(Y, 2)))

    println()
    _fact("solver status", string(result.info))
    _fact("remaining mass", @sprintf("%.4f", Y[3, end]); unit = "(published 0.3953)")
    _fact("landing time", @sprintf("%.3f", tf); unit = "(published 1.397)")
    _fact("touchdown speed", @sprintf("%.1e", abs(Y[2, end])))
    _pause(2.0)

    @step "Plotting the descent and the thrust profile" begin
        xyplot("7 · Lunar descent", tt, Y[1, :]; name = "altitude",
             line_color = "cyan", line_width = 3)
        xyplot!("7 · Lunar descent", tt, Y[2, :]; name = "vertical speed",
              line_color = "orange", line_width = 2)
        xyplot!("7 · Lunar descent", tt, Y[3, :]; name = "mass",
              line_color = "rgb(55,255,55)", line_width = 2)
        panel!("7 · Lunar descent"; xaxis_title = "time")

        xyplot("7 · Thrust", tt, U[1, :]; name = "thrust", line_shape = "hv",
             line_color = "rgb(255,80,80)", line_width = 3, fill = "tozeroy")
        panel!("7 · Thrust"; xaxis_title = "time", yaxis_title = "thrust")
    end
    println()
    _say("  Coast, then full throttle to the ground. The switch time is what the solve found;")
    _say("  nothing in the script set it.", color = _DIM)
    _pause(3.5)
    return nothing
end

# ═══════════════════════════════════════════════════════════════════════════════
# The tour
# ═══════════════════════════════════════════════════════════════════════════════

"""
    epicycle_demo(; pace = 1.0, wait_for_enter = false, open_browser = true)

Run a narrated tour of Epicycle: propagation with a 3D view, coordinate frames and a ground
track, a targeted Hohmann transfer, a GEO transfer written as a target block, orbit determination,
an Earth-Moon L1 to L2 Lyapunov transfer, and a minimum-fuel lunar landing.

# Arguments
- `pace`: scales every pause and the typing speed. `0` removes them, `2` doubles them.
- `wait_for_enter`: pause between steps until Enter is pressed, for presenting.
- `open_browser`: open the dashboard in a browser tab before the first plot.

# Returns
`nothing`. The narration is printed to the terminal and shown in a caption strip on the
EpicycleIO dashboard, beside the graphics.

# Notes
The first run downloads SPICE kernels and Earth orientation data if they are not already cached.
"""
function epicycle_demo(; pace::Real = 1.0, wait_for_enter::Bool = false,
                         open_browser::Bool = true)
    _TABS[] = open_browser
    pace >= 0 || throw(ArgumentError("pace must be zero or positive; got $pace"))
    _PACE[]  = Float64(pace)
    _ENTER[] = wait_for_enter

    t_start = time()
    EpicycleIO.clear_all!()
    if open_browser
        printstyled("  Each act opens a browser tab of its own; the caption strip tells the same story.\n";
                    color = _DIM)
    end

    _title_card()
    _pause(2.0)

    # Library notices, such as the one-time note that a burn draws from total mass, are right
    # for a script and noise in a narrated tour. Errors still show.
    Base.CoreLogging.with_logger(Base.CoreLogging.SimpleLogger(stderr, Base.CoreLogging.Error)) do
        sat = _act_orbit()
        _act_frames(sat)
        _act_targeting()
        _act_target_block()
        _act_estimation()
        _act_lyapunov()
        _act_landing()
    end

    c = _CAPTION
    c.title, c.subtitle, c.code, c.status = "That was Epicycle", "", "", ""
    empty!(c.panels)
    empty!(c.lines); empty!(c.facts); c.fresh = false

    println("\n")
    printstyled("  ", "━"^78, "\n"; color = _ACCENT)
    printstyled("   THAT WAS EPICYCLE\n"; color = _ACCENT, bold = true)
    printstyled("  ", "━"^78, "\n\n"; color = _ACCENT)
    _fact("parts", "7")
    _fact("panels on the dashboard", string(length(EpicycleIO.panels())))
    _fact("wall time", @sprintf("%.0f", time() - t_start); unit = "s")
    println()
    _say("  Every panel stays live. Scroll the dashboard, fly the 3D views, zoom the plots.")
    _say("  The source of this tour is a single file you can read and change.", color = _DIM)
    printstyled("    dashboard  ", EpicycleIO.dashboard_url(), "\n"; color = _DIM)
    println()
    return nothing
end

end # module
