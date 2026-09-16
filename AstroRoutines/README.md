# AstroRoutines

The AstroRoutines module provides low-level routines for classical astrodynamics. It currently 
covers conversions between the mean, eccentric, and true anomalies of an elliptic orbit, and the
circular restricted three-body problem in the rotating frame. For the three-body problem it 
provides the equations of motion, their Jacobian, the Jacobi constant, the five libration points,
and the variational equations for the state transition matrix. The routines are primarily
function-based rather than struct-based. The anomaly conversions follow Vallado (2013), and the 
three-body formulation follows Koon, Lo, Marsden, and Ross (2000). The three-body routines are 
verified against the NASA/JPL SSD three-body periodic orbit database.

## Installation

```julia
using Pkg
Pkg.add("AstroRoutines")
```

## Usage

Full documentation is published for this package.
