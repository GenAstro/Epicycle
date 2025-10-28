```@meta
CurrentModule = AstroFun
```

# Mission Analysis Functions

This page covers the calculations you can perform for spacecraft mission analysis. Each calculation is designed to be intuitive - you don't need to understand the underlying Julia mechanics.

## Setting Up Calculations

Start your mission analysis by setting up calculation objects:

```@docs
OrbitCalc
ManeuverCalc
BodyCalc
Constraint
```

## Getting and Setting Values

These are the main functions you'll use to get current orbital parameters and set target values for mission planning:

```@docs
func_eval
set_calc!
```

## Orbital Element Calculations

Calculate the fundamental parameters that describe spacecraft orbits:

### Orbit Size and Shape
```@docs
SMA
Ecc
```

### Orbit Orientation
```@docs
Inc
RAAN
```

### Spacecraft Position in Orbit
```@docs
TA
```

## Position and Velocity Analysis

Determine where your spacecraft is and how fast it's moving:

### Position Information
```@docs
PositionVector
PosMag
PosX
PosY
PosZ
```

### Velocity Information
```@docs
VelocityVector
VelMag
PosDotVel
```

## Maneuver Planning

Calculate fuel costs and requirements for orbital maneuvers:

```@docs
DeltaVMag
DeltaVVector
```

## Gravitational Environment

Work with different celestial bodies for mission analysis:

```@docs
GravParam
```