# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # Extended Kalman Filter Orbit Determination
#'
#' Estimate an orbit from two-way range and Doppler measurements. The example
#' simulates six hours of tracking from two ground stations, processes the
#' observations with an extended Kalman filter, and applies a Rauch-Tung-Striebel
#' smoother to the filtered solution.

using Epicycle
using LinearAlgebra

#' ## Configuration
#'
#' Define the tracking stations, force model, and truth trajectory used to
#' generate the measurements.

# Define the ground stations
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

# Configure the truth and estimation propagators
forces = ForceModel(PointMassGravity(earth, ()))
integ = IntegratorConfig(DP8();
                        dt = 60.0,
                        reltol = 1e-12,
                        abstol = 1e-12)
propagator() = OrbitPropagator(forces, integ)

# Set the truth state at the estimation epoch
epoch   = Time("2020-03-01T00:00:00.000", TT(), ISOT())
y_truth = [6878.137, 0.0, 0.0, 0.0, 4.71754, 5.99820]
truth = Spacecraft(state = CartesianState(copy(y_truth)),
                   time = epoch,
                   name = "Sat (truth)")

#' ## Simulate the Tracking Data
#'
#' Define two-way range and Doppler measurements between each station and the
#' spacecraft. The simulated one-sigma noise is 15 m for range and 2 cm/s for
#' Doppler.

# Build the range and Doppler measurement models
range_data(sc) = [TwoWayRange(SignalPath(gs, sc, gs);
                              noise = MeasurementNoise(15.0e-3)) for gs in stations]

doppler_data(sc) = [TwoWayDoppler(SignalPath(gs, sc, gs);
                                  noise = MeasurementNoise(2.0e-5)) for gs in stations]

tracking(sc) = AbstractMeasurement[range_data(sc); doppler_data(sc)]

# Simulate measurements every 30 seconds for six hours
records = simulate(tracking(truth), truth, propagator(), 30.0:30.0:21600.0; seed = 42)

#' ## Run the Extended Kalman Filter
#'
#' Begin from a perturbed Cartesian state and specify the a priori covariance and
#' process noise. `Sequential` processes the observations forward in time, then
#' the RTS smoother updates the estimates backward through the arc.

# Perturb the truth state to form the a priori estimate
guess = y_truth .+ [1.0, -1.0, 0.5, 1e-3, -1e-3, 5e-4]
sat   = Spacecraft(state = CartesianState(copy(guess)),
                   time = epoch,
                   name = "Sat (estimate)")

# Vary the spacecraft state and assign its covariance and process noise
y0 = Vary(state, sat;
          guess = guess,
          covariance = [1e2, 1e2, 1e2, 1e-2, 1e-2, 1e-2],
          process_noise = DiagonalSNC([0.0, 0.0, 0.0, 1e-12, 1e-12, 1e-12]))

# Assemble the orbit-determination problem
problem = ODProblem(spacecraft = sat,
                    propagator = propagator(),
                    measurements = tracking(sat),
                    solve_for = [y0])

# Run the EKF and RTS smoother
fit = solve!(problem, records; method = Sequential(iterations = 1, smoother = RTS()))

#' ## Report the Residuals
#'
#' Separate the postfit range and Doppler residuals, then compare their settled
#' RMS values with the noise used to generate the measurements.

# Collect postfit residuals by measurement type
range_post = [r.postfit[1] for (r, obs) in zip(fit.ekf.records, records)
              if obs.measurement_type === :RANGE]

doppler_post = [r.postfit[1] for (r, obs) in zip(fit.ekf.records, records)
                if obs.measurement_type === :DOPPLER]

# Compute RMS after the filter has settled
settled(v) = (tail = v[(2 * length(v) ÷ 3 + 1):end]; sqrt(sum(abs2, tail) / length(tail)))

# Compare postfit RMS with the simulated measurement noise
println("records             : ", length(records))
println("postfit RMS, range  : ", round(settled(range_post) * 1e3, digits = 2),
        " m    (noise 15 m)")
println("postfit RMS, Doppler: ", round(settled(doppler_post) * 1e6, digits = 2),
        " mm/s (noise 20 mm/s)")
