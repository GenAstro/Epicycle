```@meta
CurrentModule = AstroProp
```

# Mission Analysis Functions

This page covers the orbit propagation and trajectory analysis functions you can use for spacecraft mission planning. Each function is designed for practical mission analysis without requiring deep knowledge of numerical integration.

## Setting Up Propagation

Start your trajectory analysis by defining the simulation environment:

```@docs
DynSys
ForceModel
IntegratorConfig
```

## Running Simulations

Core functions for propagating spacecraft trajectories:

```@docs
propagate
OrbitPropagator
```

## Force Models

Physical models that affect spacecraft motion:

### Gravitational Forces
```@docs
PointMassGravity
TwoBodyGravity
```

### Environmental Forces
```@docs
ExponentialAtmosphere
```

## Stop Conditions

Define when to end propagation based on mission events:

### Orbital Events
```@docs
StopAtApoapsis
StopAtPeriapsis
StopAtAscendingNode
```

### Time-Based Events
```@docs
StopAtDays
StopAtSeconds
```

### Position-Based Events
```@docs
StopAtRadius
```

### Custom Events
```@docs
StopAt
```

## State Management

Handle spacecraft states during propagation:

```@docs
PosVel
```

## Advanced Features

Specialized functions for complex mission analysis:

```@docs
nbody_perts
compute_point_mass_gravity!
evaluate
accel_eval!
```