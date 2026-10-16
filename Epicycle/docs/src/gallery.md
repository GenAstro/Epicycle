# Gallery

Output from the example scripts that ship in `Epicycle/examples/` and from work in development.

## A Trajectory on the Globe

A view of Earth and an orbit.

![Earth view with trajectory](assets/gallery/Ex_EarthViewWithTrajectory.png)

## Working in the Editor

This is `Ex_HohmannTransfer.jl` open in VS Code with the solver report in the terminal.

![Epicycle in VS Code](assets/gallery/Ex_IDE_View.png)

## OD Solutions from Real Tracking Data

Two days of the Lunar Reconnaissance Orbiter from the S-band Doppler and range of the NASA ground
network, differenced with the LOLA team's GEODYN reconstruction and the Flight Dynamics Facility's
daily reconstruction in the radial, along-track and cross-track axes. Against GEODYN the RMS
difference is 0.3 m radial, 4.4 m along-track and 2.8 m cross-track.

![LRO against the GEODYN and FDF orbits](assets/gallery/LRO_OD_DefinitiveSolutionCompare.png)

A month of Juno's 2014 cruise from DSN two-way and three-way Doppler and sequential range, with
the uplink ramps, the DSN's troposphere and ionosphere calibrations and the spacecraft's spin in the
measurement models. Eleven attitude slews from the small-forces file enter the fit as impulsive
maneuvers, the same events used in mission design, and their ΔV is estimated. Post-fit residuals of
1.0 mHz and 0.23 m.

![Juno DSN residuals](assets/gallery/Juno_Slews_Residuals.png)

The same fit differenced with JPL's reconstruction, with the slews estimated, held at the file's
values and left out. Estimating them brings the RMS difference to 1.4 km radial, 7.1 km along-track
and 20 km cross-track.

![Juno against JPL's reconstruction](assets/gallery/Juno_Slews_Compare.png)

A week of LAGEOS-1 laser ranging from the ILRS normal points, with the IERS station displacements,
sub-daily Earth orientation and the Mendes-Pavlis troposphere. Post-fit residuals of 6.5 mm over
twenty-one stations; the hollow points are the normal points the editing left out.

![LAGEOS range residuals](assets/gallery/LAGEOS_Residuals.png)

Thirty-five days of MESSENGER's approach to Mercury in 2011 from DSN Doppler, range and ΔDOR, with
an eight-plate sunshade and solar panels oriented by the reconstructed attitude, and Moyer's solar
corona model. Post-fit residuals in Hz, meters and seconds.

![MESSENGER residuals](assets/gallery/MESSENGER_Residuals.png)

The same fit differenced with the navigation team's od431 reconstruction: 1.0 km radial, 2.1 km
along-track, 1.2 km cross-track and 13 nrad on the plane of the sky. Without ΔDOR the cross-track
difference grows to 13.5 km.

![MESSENGER against od431, with and without ΔDOR](assets/gallery/MESSENGER_DDOR_Compare.png)

## Propagation and Targeting

A Hohmann transfer between two circular orbits.

![Hohmann transfer](assets/gallery/Ex_HohmannTransfer.png)

A geosynchronous orbit transfer.

![GEO transfer](assets/gallery/Ex_GeoTransfer.png)

A close-proximity transfer onto a safety ellipse, in the chief's rotating frame. The arrows are
the thrust direction.

![Safety ellipse transfer](assets/gallery/Ex_SafetyEllipseTransfer.png)

## Optimal Control

Low-thrust orbit raising.

![Low thrust orbit raising](assets/gallery/Ex_LowThrustOrbitRaising.png)

Low-thrust orbit raising with three phases, each with a different transcription — zero-order hold, Hermite-Simpson and
Legendre-Gauss-Lobatto.

![Hybrid transcription](assets/gallery/Ex_LowThrustHybridTranscription.png)

The Goddard rocket problem.

![Goddard rocket](assets/gallery/Ex_GoddardRocket.png)

The linear tangent steering problem.

![Linear tangent steering](assets/gallery/Ex_LinearTangentSteering.png)

A minimum-time route between two keep-out zones, where the clearance is a path constraint.

![Obstacle avoidance](assets/gallery/Ex_ObstacleAvoidance.png)

A transfer between L1 and L2 Lyapunov orbits.

![Lyapunov transfer](assets/gallery/Ex_LypunovTransfer.png)

## Interplanetary

An Earth to Apophis rendezvous with Sims-Flanagan transcription.

![Apophis rendezvous](assets/gallery/Ex_ApophisRendezvousSimsFlanagan.png)

An Earth-Earth-Venus gravity assist.

![Gravity assist](assets/gallery/Ex_GravityAssistMGA.avif)

## Dynamical Systems

A preview of work in progress; this example is not yet in `Epicycle/examples/`.

A Sun-Earth L2 halo orbit with an out-of-plane amplitude of about 200,000 km, started from JPL's
three-body periodic orbit catalog, with its stable manifold in the circular restricted three-body
problem. The red trajectories are the transfers from near Earth with the lowest perigee.

![Halo orbit and its stable manifold](assets/gallery/HaloManifolds.png)

The Sun-Earth L2 northern halo family, continued from that orbit with BifurcationKit.jl.

![Halo orbit family](assets/gallery/HaloOrbitFamily.png)

## Orbit Determination

A batch least-squares fit to six hours of range and Doppler tracking from two ground stations.

![Batch orbit determination](assets/gallery/Ex_OrbitDetermination.png)

A filter set up to process measurements as they arrive.

![Filter evolution](assets/gallery/Ex_FilterEvolution.png)

A filter-smoother solution.

![Filter and smoother](assets/gallery/Ex_FilterSmoother.png)

Stepping the filter one observation at a time, with the postfit residuals and the formal sigma
recorded at each update.

![Stepped Kalman filter](assets/gallery/Ex_SteppedKalmanFilter.png)


## Plotting

A gallery of supported data visualization styles.

![Plot examples](assets/gallery/Ex_PlotExamples.png)
