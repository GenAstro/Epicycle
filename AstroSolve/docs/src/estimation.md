```@meta
CurrentModule = AstroSolve
```

# Estimation

AstroSolve estimates spacecraft states and other properties from measurement data. `Batch` fits a
complete observation arc with batch least squares. `Sequential` processes observations with an
extended Kalman filter and can apply a Rauch-Tung-Striebel smoother after the forward pass. Both
methods solve an `ODProblem` built from the same propagator, measurement models, and solve-for
variables.

The examples below use simulated data so the estimation error can be checked against truth. The
same solve interface accepts observations read from a CCSDS Tracking Data Message.

## Tracking Scenario

The example estimates a Cartesian state from six hours of two-way range and Doppler measurements.
Two ground stations track a spacecraft propagated with point-mass gravity.

```julia
using Epicycle
using LinearAlgebra

goldstone = GroundStation(
    name = "DSS-14",
    body = earth,
    latitude = 35.4267,
    longitude = -116.89,
    altitude = 1.0,
    min_elevation = 5.0,
)

canberra = GroundStation(
    name = "DSS-43",
    body = earth,
    latitude = -35.4,
    longitude = 148.98,
    altitude = 0.7,
    min_elevation = 5.0,
)

stations = (goldstone, canberra)

forces = ForceModel(PointMassGravity(earth, ()))
integ = IntegratorConfig(DP8(); dt = 60.0, reltol = 1e-12, abstol = 1e-12)
propagator() = OrbitPropagator(forces, integ)

epoch = Time("2020-03-01T00:00:00.000", TT(), ISOT())
y_truth = [6878.137, 0.0, 0.0, 0.0, 4.71754, 5.99820]

truth = Spacecraft(
    state = CartesianState(copy(y_truth)),
    time = epoch,
    name = "Sat (truth)",
)
```

## Measurement Models and Tracking Data

`SignalPath` identifies the participants in a measurement. The range and Doppler models below use
the same station for transmission and reception. `MeasurementNoise` takes a one-sigma value in the
measurement's units: kilometers for range and kilometers per second for Doppler.

```julia
range_models(sc) = [
    TwoWayRange(
        SignalPath(station, sc, station);
        noise = MeasurementNoise(15.0e-3),
    ) for station in stations
]

doppler_models(sc) = [
    TwoWayDoppler(
        SignalPath(station, sc, station);
        noise = MeasurementNoise(2.0e-5),
    ) for station in stations
]

tracking(sc) = AbstractMeasurement[range_models(sc); doppler_models(sc)]

records = simulate(
    tracking(truth),
    truth,
    propagator(),
    30.0:30.0:21600.0;
    seed = 42,
)
```

Each observation record retains its epoch, measurement type, and participants. The estimator uses
those fields to select the matching measurement model and orders records by epoch before processing
them.

### CCSDS Tracking Data

`TrackingDataFile` identifies a tracking file and its format. AstroSolve reads CCSDS TDM KVN files
into the same observation records produced by `simulate`.

```julia
file = TrackingDataFile(tempname() * ".tdm", CCSDS_KVN())

write_records(file, records,
              TDMHeader("2020-03-01T12:00:00", "GEN ASTRO"),
              [TDMSegmentMeta(time_system = "TT", participant_1 = "DSS-14",
                              participant_2 = "Sat"),
               TDMSegmentMeta(time_system = "TT", participant_1 = "DSS-43",
                              participant_2 = "Sat")])

tdm_records, header, segments = read_records(file)
```

A TDM may contain segments from several stations. Segment metadata preserves the participant names
used to match records to measurement models. `solve!` also accepts the `TrackingDataFile` directly
and reads it before estimation. `write_records` writes observation records with a `TDMHeader` and
the corresponding `TDMSegmentMeta` entries.

## Build the Estimation Problem

`Vary` declares the quantities to estimate and their a priori uncertainty. A scalar covariance is a
variance applied to every component, a vector supplies diagonal variances, and a matrix supplies the
complete covariance. These values are variances, unlike the standard deviation passed to
`MeasurementNoise`.

The setup function below creates a fresh estimate for each method. This matters because a solve
updates the varied spacecraft state.

```julia
prior_variance = [1e2, 1e2, 1e2, 1e-2, 1e-2, 1e-2]
initial_offset = [1.0, -1.0, 0.5, 1e-3, -1e-3, 5e-4]

function estimation_problem(; process_noise = nothing)
    guess = y_truth .+ initial_offset
    sat = Spacecraft(
        state = CartesianState(copy(guess)),
        time = epoch,
        name = "Sat (estimate)",
    )

    state_variable = Vary(
        state,
        sat;
        guess = guess,
        covariance = prior_variance,
        process_noise = process_noise,
    )

    problem = ODProblem(
        spacecraft = sat,
        propagator = propagator(),
        measurements = tracking(sat),
        solve_for = [state_variable],
    )

    return (; problem, sat, state_variable)
end
```

The `ODProblem` holds the estimated spacecraft, its propagation model, the available measurement
models, and the solve-for declarations. Additional estimated properties are added to `solve_for`
with their own `Vary` declarations.

## Batch Least Squares

Batch least squares processes every observation on each iteration. `n_iters` limits the number of
linearization and correction cycles, and `tol` tests the relative correction between iterations.

```julia
batch_case = estimation_problem()

batch_fit = solve!(
    batch_case.problem,
    records;
    method = Batch(n_iters = 10, tol = 1e-9),
)
```

### Check the Batch Solution

The result contains the estimated state and parameters, posterior covariance, formal standard
deviations, correlation matrix, postfit residuals, and convergence status.

```julia
batch_error = batch_fit.X_hat .- y_truth

println("converged      : ", batch_fit.converged)
println("iterations     : ", batch_fit.iters)
println("position error : ", round.(batch_error[1:3] .* 1e3, digits = 2), " m")
println("velocity error : ", round.(batch_error[4:6] .* 1e6, digits = 2), " mm/s")
println("position sigma : ", round.(batch_fit.sigma[1:3] .* 1e3, digits = 2), " m")
println("velocity sigma : ", round.(batch_fit.sigma[4:6] .* 1e6, digits = 2), " mm/s")
```

`batch_fit.P_hat` and `batch_fit.sigma` cover the solve-for components in
`batch_fit.solve_for_idx`. `batch_fit.X_hat` is the complete estimated state and parameter vector.

## Sequential Estimation

The extended Kalman filter propagates the state and covariance to each observation and then applies
a measurement update. Process noise is attached to the varied quantity because it describes
uncertainty added while that quantity is propagated, not uncertainty in an observation.

`DiagonalSNC` takes one white-noise power spectral density per component. Each value has units of
the squared component unit per unit time; the filter accumulates it over each propagation step.

```julia
process_noise = DiagonalSNC([
    0.0, 0.0, 0.0,
    1e-12, 1e-12, 1e-12,
])

sequential_case = estimation_problem(; process_noise = process_noise)

sequential_fit = solve!(
    sequential_case.problem,
    records;
    method = Sequential(iterations = 1, smoother = RTS()),
)
```

### Rauch-Tung-Striebel Smoothing

`RTS()` adds a backward pass after the forward filter. The forward result is available as
`sequential_fit.ekf`; the smoothed states and covariances at each observation epoch are available as
`sequential_fit.rts`. With `iterations = 1`, the method performs one forward and backward pass.
Larger values repeat the process from the revised initial estimate until `tol` is met or the
iteration limit is reached.

### Check the Sequential Solution

Each forward-filter record contains the predicted and updated state and covariance, state
transition matrix, Kalman gain, and prefit and postfit residuals for one observation.

```julia
range_postfit = [
    update.postfit[1]
    for (update, observation) in zip(sequential_fit.ekf.records, records)
    if observation.measurement_type === :RANGE
]

doppler_postfit = [
    update.postfit[1]
    for (update, observation) in zip(sequential_fit.ekf.records, records)
    if observation.measurement_type === :DOPPLER
]

settled_rms(values) = begin
    tail = values[(div(2 * length(values), 3) + 1):end]
    sqrt(sum(abs2, tail) / length(tail))
end

println("range RMS  : ", round(settled_rms(range_postfit) * 1e3, digits = 2), " m")
println("Doppler RMS: ", round(settled_rms(doppler_postfit) * 1e6, digits = 2), " mm/s")
```

The residual records belong to the forward filter. `sequential_fit.ekf.P_hat` is the covariance at
the final observation, while `sequential_fit.rts.P_smooth` contains the smoothed covariance across
the arc.

## Step the Filter

The stateful interface exposes the estimate between observations. Use it when an application must
record intermediate values, change settings, or perform another action between updates.

```julia
stepped_case = estimation_problem(; process_noise = process_noise)

dyn, meas, obs_times, obs_data, R_per_obs =
    build_od_closures(records, stepped_case.problem)

ekf = init_ekf(
    [stepped_case.state_variable],
    dyn,
    meas;
    model = stepped_case.problem,
    R = R_per_obs[1],
    t0 = 0.0,
)

postfit = Vector{Float64}(undef, length(obs_times))
position_sigma = similar(postfit)

for k in eachindex(obs_times)
    time_update!(ekf, obs_times[k])
    update = measurement_update!(ekf, obs_data[k]; R = R_per_obs[k])

    postfit[k] = update.postfit[1]
    position_sigma[k] = sqrt(current_covariance(ekf)[1, 1])
end
```

The filter starts at the spacecraft epoch, `t0 = 0.0`, and propagates to the first observation.
Starting it at the first observation would associate the a priori state with the wrong epoch.

## Evaluating an Estimate

A completed solve is not enough to establish that an estimate is useful. Check convergence where
the method reports it, inspect prefit and postfit residuals for bias and outliers, compare settled
residual statistics with the measurement noise, and confirm that the covariance behaves as
expected. Results should also be tested against reasonable changes to the a priori state,
covariance, and process noise.

## Solved Examples

- [Batch least squares](https://github.com/GenAstro/Epicycle/blob/main/Epicycle/examples/Ex_OrbitDetermination.jl) estimates an initial state from a complete tracking arc.
- [Extended Kalman filter](https://github.com/GenAstro/Epicycle/blob/main/Epicycle/examples/Ex_ExtendedKalmanFilter.jl) processes range and Doppler data and applies RTS smoothing.
- [Stepped Kalman filter](https://github.com/GenAstro/Epicycle/blob/main/Epicycle/examples/Ex_SteppedKalmanFilter.jl) exposes the state, covariance, and residual after every observation.
