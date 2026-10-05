```@meta
CurrentModule = AstroProp
```

# API Reference

Every public name in AstroProp, with its full signature, grouped by what it is for. How to choose
and combine the forces is explained on the [Force Models](force_models.md) page, and how to set up
and run a propagation on the [AstroProp](index.md) page.

## Propagation

The propagator, the integrator settings, the stopping conditions, and the problem and
variational-equation types under them.

```@autodocs
Modules = [AstroProp]
Order   = [:type, :function, :macro, :constant]
Public  = true
Private = false
Pages   = ["AstroProp.jl", "orbit_propagator.jl", "orbit_ode_problem.jl", "variational.jl",
           "jacobian_config.jl", "discontinuities.jl"]
```

## Forces

The forces, the gravity fields and atmospheres open AstroProp provides, and the interfaces a new
force, field or atmosphere implements.

```@autodocs
Modules = [AstroProp]
Order   = [:type, :function, :macro, :constant]
Public  = true
Private = false
Pages   = ["point_mass_gravity.jl", "harmonic_gravity.jl", "zonal_gravity.jl",
           "atmospheric_drag.jl", "exponential_atmosphere.jl", "spherical_srp.jl",
           "external_force.jl", "force_context.jl"]
```

A new gravity field declares the axes its coefficients are defined in by extending this function,
which is not exported:

```@docs
AstroProp.field_orientation
```

## Spacecraft geometry

The drag and SRP properties a spacecraft carries, defined in AstroModels.

```@docs
SphericalDrag
SphericalSRP
```

## Index

```@index
```
