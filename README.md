# Welcome to Epicycle

[![CI](https://github.com/GenAstro/Epicycle/workflows/CI/badge.svg)](https://github.com/GenAstro/Epicycle/actions)
[![codecov](https://codecov.io/gh/GenAstro/Epicycle/branch/main/graph/badge.svg?token=FNHOVC5O5N)](https://codecov.io/gh/GenAstro/Epicycle)
[![Documentation (stable)](https://img.shields.io/badge/docs-stable-blue.svg)](https://genastro.github.io/Epicycle/Epicycle/stable/)
[![Documentation (dev)](https://img.shields.io/badge/docs-dev-blue.svg)](https://genastro.github.io/Epicycle/Epicycle/dev/)

Epicycle is an application for space systems with a nod to the giants before us and to how far we've come. 

## Documentation

### The Epicycle Application

- **[Epicycle Application](https://genastro.github.io/Epicycle/Epicycle/dev/)** - Complete application guide
- **[Video Overview](https://www.youtube.com/watch?v=Ph51mj0SVP4)** - Video walkthrough of the application

### Component Documentation

- **[EpicycleBase Documentation](https://genastro.github.io/Epicycle/EpicycleBase/dev/)** - Foundation types and abstractions
- **[AstroStates Documentation](https://genastro.github.io/Epicycle/AstroStates/dev/)** - State vector representations and conversions
- **[AstroEpochs Documentation](https://genastro.github.io/Epicycle/AstroEpochs/dev/)** - Time systems and epoch handling
- **[AstroUniverse Documentation](https://genastro.github.io/Epicycle/AstroUniverse/dev/)** - Celestial body models and ephemeris
- **[AstroFrames Documentation](https://genastro.github.io/Epicycle/AstroFrames/dev/)** - Coordinate systems and transformations
- **[AstroModels Documentation](https://genastro.github.io/Epicycle/AstroModels/dev/)** - Spacecraft and physical models
- **[AstroManeuvers Documentation](https://genastro.github.io/Epicycle/AstroManeuvers/dev/)** - Maneuver models and algorithms
- **[AstroCallbacks Documentation](https://genastro.github.io/Epicycle/AstroCallbacks/dev/)** - I/O, stopping conditions, cost, and constraints
- **[AstroProp Documentation](https://genastro.github.io/Epicycle/AstroProp/dev/)** - Trajectory propagation algorithms
- **[AstroSolve Documentation](https://genastro.github.io/Epicycle/AstroSolve/dev/)** - Optimization and constraint solving

## What's New

**AstroStates**
- Brouwer mean-element state types (`BrouwerMeanShortState`, `BrouwerMeanLongState`).
- `MeanSMA` calc (in `AstroCallbacks`) for targeting and reporting the Brouwer long-period
  mean semi-major axis.

**AstroProp**
- Zonal Earth gravity (J₂–J₅) — `HarmonicGravity(earth; model = Zonal())`. Cross-validated
  against GMAT EGM96 at degree 5.
- Exponential Earth atmosphere — `Exponential()` density model for `AtmosphericDrag`.
- Spherical atmospheric drag — `AtmosphericDrag` force, reading `SphericalDrag` geometry
  from the spacecraft (see AstroModels).
- Spherical solar radiation pressure — `SolarRadiationPressure` force with a dual-cone
  eclipse shadow, reading `SphericalSRP` geometry from the spacecraft (see AstroModels).
- Faster and more flexible stopping conditions — `StopAt` gains `detection` (`:discrete`
  default, `:continuous` escape hatch) and `rootfind_tol` kwargs. Default `:discrete` polls
  once per accepted step and bisects on the Vern9 interpolant when a sign change appears —
  same root precision, ~2× faster on typical LEO drag workloads. Drop-in for existing calls.

**Epicycle** (umbrella)
- Earth station-keeping example — `Epicycle/examples/Ex_StationKeeping.jl`. LEO satellite
  maintained above a mean-SMA trigger with periodic Hohmann re-boosts; solver refines
  analytic ΔV guesses against a mean-SMA constraint at MOI. Demonstrates the full stack:
  states, epochs, forces, propagation, maneuvers, and sequences composed via `using Epicycle`.

**AstroModels**
- `SphericalDrag` spacecraft geometry (`SphericalDrag(; c_d, drag_area)`) — isotropic
  (spherical-body) drag geometry attached to `Spacecraft.drag`.
- `SphericalSRP` spacecraft geometry (`SphericalSRP(; c_r, srp_area)`) — isotropic
  (spherical-body) SRP geometry attached to `Spacecraft.srp`.
- `Spacecraft.save_history` field — opt out of trajectory-segment storage during
  propagation. `sc.state` and `sc.time` still update; only the segment push is skipped.

**Enterprise: EpicycleEnterprise**
- Full-field spherical-harmonic gravity — `EGM96` and `EGM2008` via `HarmonicGravity`.
- NRLMSISE-00 atmospheric density — `MSISE00` via `AtmosphericDrag`.

## Contributing to Epicycle 

Contributing is easy.

1. Fork the project
2. Create a new feature branch
3. Make your changes
4. Submit a pull request

We use the Linux Kernel's Developer's Certificate of Origin (DCO) as detailed in CONTRIBUTING.txt.

## License

We believe in the power of open source to foster innovation and community-driven 
development and also recognize the need for a sustainable business model and a model
that can handle export-controlled aerospace content. 

For these reasons, Epicycle is offered under a tri-licensing model. The license allows
users to choose between the following three options:

1) LGPL V3.0
2) Evaluation and Education use Only
3) Commercial License

See LICENSE.txt for terms each license option.  For commercial licensing, 
email licensing at genastro.org.

## Notes
Claude Sonnet and ChatGPT are used in the development of Epicycle.