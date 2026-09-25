# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # Stepped Extended Kalman Filter Orbit Determination
#'
#' Estimate an orbit from six hours of two-way range and Doppler data, advancing an extended Kalman
#' filter one observation at a time. The stepped interface exposes the state, covariance, and
#' residual between updates.
#'
#' This is the same scenario as the extended Kalman filter example, which hands the whole arc to
#' `solve!` and reads the answer at the end. Stepping it is what a navigation filter running
#' alongside a spacecraft does, and what an example needs when something has to happen between one
#' observation and the next.

using Epicycle
using LinearAlgebra

#' ## Configuration
#'
#' Configure two ground stations, a truth spacecraft, and the orbit propagator.

# Configure the tracking stations
goldstone = GroundStation(name = "DSS-14",
                          body = earth,
                          latitude = 35.4267,
                          longitude = -116.89,
                          altitude = 1.0,
                          min_elevation = 5.0)

canberra = GroundStation(name = "DSS-43",
                         body = earth,
                         latitude = -35.4,
                         longitude = 148.98,
                         altitude = 0.7,
                         min_elevation = 5.0)

stations = (goldstone, canberra)

# Configure the orbit propagator
forces = ForceModel(PointMassGravity(earth, ()))
integ = IntegratorConfig(DP8();
                        dt = 60.0,
                        reltol = 1e-12,
                        abstol = 1e-12)
propagator() = OrbitPropagator(forces, integ)

# Set the truth state
epoch   = Time("2020-03-01T00:00:00.000", TT(), ISOT())
y_truth = [6878.137, 0.0, 0.0, 0.0, 4.71754, 5.99820]
truth = Spacecraft(state = CartesianState(copy(y_truth)),
                   time = epoch,
                   name = "Sat (truth)")

#' ## Simulate the tracking data
#'
#' Simulate noisy two-way range and Doppler measurements from both stations.

# Set the measurement noise to 15 m for range and 2 cm/s for Doppler
range_data(sc) = [TwoWayRange(SignalPath(gs, sc, gs);
                              noise = MeasurementNoise(15.0e-3)) for gs in stations]

doppler_data(sc) = [TwoWayDoppler(SignalPath(gs, sc, gs);
                                  noise = MeasurementNoise(2.0e-5)) for gs in stations]

tracking(sc) = AbstractMeasurement[range_data(sc); doppler_data(sc)]

# Simulate six hours of tracking data
# TODO - interface cleanup
records = simulate(tracking(truth), truth, propagator(), 30.0:30.0:21600.0; seed = 42)

#' ## Build the estimation problem
#'
#' Perturb the truth state and assign its a priori covariance and process noise.

# Form the a priori estimate with position and velocity errors
guess = y_truth .+ [1.0, -1.0, 0.5, 1e-3, -1e-3, 5e-4]
sat   = Spacecraft(state = CartesianState(copy(guess)),
                   time = epoch,
                   name = "Sat (estimate)")

# Vary the spacecraft state with its covariance and process noise
y0 = Vary(state, sat;
          guess = guess,
          covariance = [1e2, 1e2, 1e2, 1e-2, 1e-2, 1e-2],
          process_noise = DiagonalSNC([0.0, 0.0, 0.0, 1e-12, 1e-12, 1e-12]))

# Build the estimation problem
problem = ODProblem(spacecraft = sat,
                    propagator = propagator(),
                    measurements = tracking(sat),
                    solve_for = [y0])

#' ## Open the filter
#'
#' `build_od_closures` returns the pieces `solve!` would have used: the dynamics, the measurement
#' model, the observation times and data, and the measurement covariance for each observation.
#'
#' The filter starts at the spacecraft epoch rather than at the first observation, because
#' observation times are seconds from that epoch and the a priori state belongs to it. Starting at
#' the first observation hands the filter a state that is already hours out of date, and the
#' covariance converges anyway, so the residuals are what show it and the sigma is not.

# Extract the dynamics, measurements, times, data, and measurement covariances
dyn, meas, obs_times, obs_data, R_per_obs = build_od_closures(records, problem)

# Initialize the filter at the spacecraft epoch
ekf = init_ekf([y0], dyn, meas;
               model = problem,
               R = R_per_obs[1],
               t0 = 0.0)

#' ## Step through the observations
#'
#' Each pass propagates to the observation time, then updates on the observation. The estimate, its
#' covariance and the residual that just arrived are all readable inside the loop.

# Process each observation and record the postfit residual and position sigma
n_obs = length(obs_times)
postfit = Vector{Float64}(undef, n_obs)
sigma_position = Vector{Float64}(undef, n_obs)

for k in 1:n_obs
    time_update!(ekf, obs_times[k])
    info = measurement_update!(ekf, obs_data[k]; R = R_per_obs[k])

    postfit[k] = info.postfit[1]
    sigma_position[k] = sqrt(diag(current_covariance(ekf))[1])
end

#' ## Report the residuals
#'

# Separate the range and Doppler residuals
is_range = [r.measurement_type === :RANGE for r in records]
range_post = postfit[is_range]
doppler_post = postfit[.!is_range]

# Compute the RMS over the last third of the arc, past the filter transient
settled(v) = (tail = v[(2 * length(v) ÷ 3 + 1):end]; sqrt(sum(abs2, tail) / length(tail)))

# Report the residuals against the noise each measurement was simulated with
println("records             : ", length(records))
println("postfit RMS, range  : ", round(settled(range_post) * 1e3, digits = 2),
        " m    (noise 15 m)")
println("postfit RMS, Doppler: ", round(settled(doppler_post) * 1e6, digits = 2),
        " mm/s (noise 20 mm/s)")
println("final sigma, position: ", round(sigma_position[end] * 1e3, digits = 2), " m")
