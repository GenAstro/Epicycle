```@meta
CurrentModule = AstroManeuvers
```

# AstroManeuvers

The AstroManeuvers module provides utilities and functions for orbital maneuver calculations in astrodynamics applications. The module includes impulsive maneuver models and functions for applying maneuvers to spacecraft objects.

## Installation

Versions through 0.2.0 are in Julia's General registry. From the next version AstroManeuvers is
released under the Gen Astro Source Available License, which General does not carry, so later
versions come from the Gen Astro registry. Add it once, then install as usual:

```julia
using Pkg
Pkg.Registry.add(RegistrySpec(url = "https://github.com/GenAstro/GenAstro.git"))
Pkg.add("AstroManeuvers")
```

General is still required, since these packages depend on packages registered there. Installing
without the Gen Astro registry resolves to 0.2.0, the last version General carries, and
reports nothing about the newer ones.

## Quick Start

Apply an impulsive orbital maneuver:

```julia
using AstroModels, AstroManeuvers
m = ImpulsiveManeuver(axes=Inertial(), 
                      Isp=300.0, 
                      element1=0.01, 
                      element2=0.0, 
                      element3=0.0)

sc = Spacecraft()
maneuver!(sc, m)
```

## Table of Contents

```@index
```

## API Reference

```@autodocs
Modules = [AstroManeuvers]
Public  = true
Private = false
Order = [:type, :function, :macro, :constant]
```





