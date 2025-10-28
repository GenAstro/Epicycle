```@meta
CurrentModule = AstroProp
```

# Developer API Reference

This page documents internal implementation details for developers contributing to AstroProp.

!!! warning "Internal API"
    These functions are implementation details and may change without notice. 
    Use the public API documented in the main reference guide for stable interfaces.

## Internal Types

```@docs
OrbitODE
CartesianODE
```

## ODE System Construction

Internal functions for building the differential equation system:

```@docs
build_odereg
build_odes!
build_state
update_structs!
```

## Force Model Internals

Low-level force calculation functions:

```@docs
find_center
```

## Propagator Internals

Implementation details for the orbit propagator:

```@docs
_calc_subject
_posvel_from_u
_build_callback
```

## Stop Condition Internals

Internal functions for event detection:

```@docs
apsis_condition
ascnode_condition
stop_affect!
```

## State Management Internals

Low-level state handling functions:

```@docs
push_history_segment!
```