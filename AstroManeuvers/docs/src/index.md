```@meta
CurrentModule = AstroManeuvers
```

# AstroManeuvers

The AstroManeuvers module provides utilities and functions for orbital maneuver calculations in astrodynamics applications. The module includes impulsive maneuver models and functions for applying maneuvers to spacecraft objects.

## Installation

Versions through 0.2.0 are in Julia's General registry. From the next version AstroManeuvers is
released under the Gen Astro Source Available License, which General does not carry, so later
versions come from the Gen Astro registry. Add both registries once, then install as usual:

```julia
using Pkg
pkg"registry add General https://github.com/GenAstro/GenAstroRegistry.git"
pkg"add AstroManeuvers"
```

General is named in that command for two reasons: the dependencies live there, and Julia
installs it by itself only while no registry is present at all, so adding the Gen Astro
registry alone would leave it out. Installing without the Gen Astro registry resolves to 0.2.0,
the last version General carries, and reports nothing about the newer ones.

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





