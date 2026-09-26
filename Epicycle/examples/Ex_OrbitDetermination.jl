# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # Batch Least Squares Orbit Determination
#'
#' Estimate a spacecraft orbit from six hours of simulated two-way range and Doppler data. The
#' batch least squares solution starts from an a priori state with position and velocity errors and
#' reports the recovered state error and posterior uncertainty.

using Epicycle
using EpicycleIO
using LinearAlgebra

#' ## Configuration
#'
#' Configure the ground station, force model, and truth spacecraft.

# Configure the tracking station
station = GroundStation(name = "DSS-14", 
                        body = earth, 
                        latitude = 35.4267,
                        longitude = -116.89, 
                        altitude = 1.0, 
                        min_elevation = 5.0)

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
#' Simulate noisy two-way range and Doppler measurements between the station and spacecraft.

# Set the measurement noise to 15 m for range and 2 cm/s for Doppler
tracking(sc) = AbstractMeasurement[
    TwoWayRange(SignalPath(station, sc, station);   noise = MeasurementNoise(15.0e-3)),
    TwoWayDoppler(SignalPath(station, sc, station); noise = MeasurementNoise(2.0e-5))]

# Simulate six hours of tracking data
# TODO - interface cleanup
records = simulate(tracking(truth), truth, propagator(), 30.0:30.0:21600.0; seed = 42)

#' ## Run the batch least squares estimator
#'

# Perturb the truth state to form the a priori estimate
guess = y_truth .+ [1.0, -1.0, 0.5, 1e-3, -1e-3, 5e-4]
sat   = Spacecraft(state = CartesianState(copy(guess)), 
                   time = epoch,
                   name = "Sat (estimate)")

# Vary the spacecraft state with its a priori covariance
y0 = Vary(state, sat; 
          guess = guess, 
          covariance = Diagonal([1e2, 1e2, 1e2, 1e-2, 1e-2, 1e-2]))

# Build the estimation problem
problem = ODProblem(spacecraft = sat, 
                    propagator = propagator(),
                    measurements = tracking(sat), 
                    solve_for = [y0])

# Solve with batch least squares
fit = solve!(problem, records; method = Batch(n_iters = 10, tol = 1e-9, verbose = true))

# Report the solution, errors, and computed posterior covariance
err = fit.X_hat .- y_truth
println("records        : ", length(records))
println("position error : ", round.(err[1:3] .* 1e3, digits = 2), " m")
println("velocity error : ", round.(err[4:6] .* 1e6, digits = 2), " mm/s")
println("formal sigma   : ", round.(fit.sigma[1:3] .* 1e3, digits = 2), " m")
println("formal sigma   : ", round.(fit.sigma[4:6] .* 1e6, digits = 2), " mm/s")

#' ## Plot the Residuals
#'
#' `fit.residuals` holds one vector per observation, in the order the records were supplied, so
#' the residuals split by measurement type the same way the records do.

# Split the converged residuals and their times by measurement type
hours   = [r.t_receive - epoch for r in records] ./ 3600
resid   = [r[1] for r in fit.residuals]
is_rng  = [r.measurement_type === :RANGE for r in records]

# Plot the range residuals in metres
xyplot("Batch Range Residuals", hours[is_rng], resid[is_rng] .* 1e3;
       mode = "markers", name = "residual")
panel!("Batch Range Residuals";
       xaxis_title = "hours from epoch", yaxis_title = "range residual (m)",
       yaxis_range = [-60.0, 60.0])     # four sigma of the 15 m simulated noise

# Plot the Doppler residuals in millimetres per second
xyplot("Batch Doppler Residuals", hours[.!is_rng], resid[.!is_rng] .* 1e6;
       mode = "markers", name = "residual")
panel!("Batch Doppler Residuals";
       xaxis_title = "hours from epoch", yaxis_title = "Doppler residual (mm/s)",
       yaxis_range = [-80.0, 80.0])     # four sigma of the 20 mm/s simulated noise
