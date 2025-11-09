```@meta
CurrentModule = AstroManeuvers
```

# AstroManeuvers

The AstroManeuvers module provides utilities and functions for orbital maneuver calculations in astrodynamics applications. The module includes impulsive maneuver models and functions for applying maneuvers to spacecraft objects.

## Quick Start

Apply an impulsive orbital maneuver:

```julia
using Epicycle
m = ImpulsiveManeuver(axes=Inertial(), 
                      Isp=300.0, 
                      element1=0.01, 
                      element2=0.0, 
                      element3=0.0)

sc = Spacecraft()
maneuver(sc, m)
```

## Table of Contents

```@index
```

## API Reference

```@autodocs
Modules = [AstroManeuvers]
Order = [:type, :function, :macro, :constant]
```





