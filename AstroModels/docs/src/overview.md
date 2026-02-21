# Overview

## Spacecraft Model

The `Spacecraft` type represents a spacecraft with:
- Orbital state (position and velocity)
- Epoch (time)
- Physical properties (mass, reflectivity)
- Trajectory history
- Coordinate system
- Optional 3D model for visualization

## History System

Trajectory history is automatically recorded during propagation and organized into segments. Segments delineate mission phases (coast, maneuver) and enable multi-phase mission analysis.

See [History Guide](history.md) for details on accessing and working with trajectory data.
