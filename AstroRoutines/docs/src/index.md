```@meta
CurrentModule = AstroRoutines
```

# AstroRoutines

The AstroRoutines module provides low-level routines for classical astrodynamics. It currently 
covers conversions between the mean, eccentric, and true anomalies of an elliptic orbit, and the
circular restricted three-body problem in the rotating frame. For the three-body problem it 
provides the equations of motion, their Jacobian, the Jacobi constant, the five libration points,
and the variational equations for the state transition matrix. The routines are primarily
function-based rather than struct-based. The anomaly conversions follow Vallado (2013), and the 
three-body formulation follows Koon, Lo, Marsden, and Ross (2000). The three-body routines are 
verified against the NASA/JPL SSD three-body periodic orbit database.

References:

- Vallado, D. A. (2013), *Fundamentals of Astrodynamics and Applications*, 4th ed., Microcosm
  Press / Springer.
- Koon, W. S., Lo, M. W., Marsden, J. E., and Ross, S. D. (2000), Heteroclinic connections between
  periodic orbits and resonance transitions in celestial mechanics, *Chaos* 10(2), 427–469.

## Quick Start

Convert between the three anomalies of an elliptic orbit. Angles are in radians.

```@example quickstart
using AstroRoutines

e = 0.2                                   # eccentricity, 0 ≤ e < 1
M = 1.0                                   # mean anomaly

E = mean_to_eccentric_anomaly(M, e)       # solves Kepler's equation
ν = eccentric_to_true_anomaly(E, e)
```

Evaluate the three-body dynamics at a state. The mass ratio is the only parameter, and the state
is `[x, y, z, vx, vy, vz]` in the rotating frame.

```@example quickstart
using AstroRoutines

mu = 1.215058560962404e-02                            # Earth–Moon
s  = [0.81596252146384562, 0.0, 0.0, 0.0, 0.20722124749217649, 0.0]

C  = jacobi_constant(s, mu)                           # the integral of motion
ds = zeros(6)
cr3bp_eom!(ds, s, mu)                                 # derivative, for a solver
C, ds[4]
```

Locate the libration points of that system.

```@example quickstart
using AstroRoutines

mu = 1.215058560962404e-02
libration_point(mu, :L1), libration_point(mu, :L4)
```

## Anomaly conversions

AstroRoutines provides all six conversions among mean anomaly `M`, eccentric anomaly `E`, and true
anomaly `ν`, so any pair converts directly rather than through a chain.

```@example anomalies
using AstroRoutines

e = 0.2
M = 1.0

E = mean_to_eccentric_anomaly(M, e)
ν = eccentric_to_true_anomaly(E, e)

# Every pair converts directly, so the round trip returns the angle you started with
true_to_mean_anomaly(ν, e) - M
```

Angles are in radians. Outputs are single-revolution and unwrapped: `M` and `E` return as the
smooth continuation of the input, and `ν` returns in the principal branch `(−π, π]`. The routines
do not apply `[0, 2π)` normalization, which would introduce a point where the derivative does not
exist.

`mean_to_eccentric_anomaly` and `mean_to_true_anomaly` iterate, and accept `tol` (default `1e-12`)
and `maxiter` (default `50`). The other four are closed form. Eccentricity outside `[0, 1)` throws
`ArgumentError`.

All six are differentiable with ForwardDiff with respect to both arguments, including the two that
iterate.

```@example anomalies
using AstroRoutines, ForwardDiff

ForwardDiff.derivative(m -> mean_to_true_anomaly(m, 0.2), 1.0)
```

## The circular restricted three-body problem

These routines use the rotating frame and normalized units: the primaries sit at `(−μ, 0, 0)` and
`(1−μ, 0, 0)`, their separation is 1, their orbital rate is 1, and the total mass is 1. The only
parameter is the mass ratio `μ`, and every state argument is the six-element vector
`[x, y, z, vx, vy, vz]`.

`cr3bp_mass_ratio` computes `μ` from the two primary masses, so the value carries whichever mass
constants you supply. The Earth–Moon databases at JPL are built on `1.215058560962404e-02`, which
is the value used below.

```@example cr3bp
using AstroRoutines

mu = cr3bp_mass_ratio(5.9722e24, 7.342e22)
```

`libration_point` takes `:L1` through `:L5` and returns the three-vector `[x, y, 0]`. The
collinear points are computed by Newton iteration rather than tabulated, since their positions move
with `μ`; `L4` and `L5` are closed form. A mass ratio outside `(0, 1)`, or a symbol other than the
five, throws `ArgumentError`.

```@example cr3bp
using AstroRoutines

mu = 1.215058560962404e-02
libration_point(mu, :L1)
```

Four functions evaluate the dynamics at a state. `cr3bp_accel` returns the acceleration,
`cr3bp_jacobian` the 6×6 `∂f/∂s`, and `jacobi_constant` the integral of motion. `cr3bp_eom!` writes
the full derivative into a six-element vector, in the signature an ODE solver expects; it is
ballistic motion only, so a propagator that applies thrust adds the control acceleration to
`ds[4:6]` after calling it.

```@example cr3bp
using AstroRoutines

mu = 1.215058560962404e-02

# An Earth–Moon L1 Lyapunov orbit from the JPL SSD database
s = [0.81596252146384562, 0.0, 0.0, 0.0, 0.20722124749217649, 0.0]

a = cr3bp_accel(s, mu)
A = cr3bp_jacobian(s, mu)
C = jacobi_constant(s, mu)

ds = zeros(6)
cr3bp_eom!(ds, s, mu)
C, ds
```

### State transition matrix

`cr3bp_stm_eom!` evaluates the equations of motion and the variational terms together, writing the
derivative of a 42-element vector: the six-element state, then the state transition matrix `Φ`
flattened column-major. The variational block is `Φ̇ = A(s) Φ`, with `A` the analytic Jacobian
`cr3bp_jacobian` returns. Integration is the caller's, so the same integrator advances the state
and `Φ` on the same steps and to the same tolerance. `cr3bp_stm_initial` packs a six-element state
into that layout with `Φ(0) = I`. Either function throws `ArgumentError` on a vector of the wrong
length.

```@example cr3bp
using AstroRoutines

mu = 1.215058560962404e-02
s0 = [0.81596252146384562, 0.0, 0.0, 0.0, 0.20722124749217649, 0.0]

z0 = cr3bp_stm_initial(s0)
dz = zeros(42)
cr3bp_stm_eom!(dz, z0, mu)
length(z0)
```

## Reference

The entry points of each area. The full list follows below.

```@docs
mean_to_eccentric_anomaly
cr3bp_mass_ratio
libration_point
cr3bp_eom!
jacobi_constant
```

## API

```@index
```

```@autodocs
Modules = [AstroRoutines]
Public  = true
Private = false
Order   = [:function, :type, :constant]
Filter  = t -> !(t in (mean_to_eccentric_anomaly, cr3bp_mass_ratio, libration_point,
                       cr3bp_eom!, jacobi_constant))
```
