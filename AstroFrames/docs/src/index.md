```@meta
CurrentModule = AstroFrames
```

# AstroFrames

The AstroFrames package provides coordinate systems and reference frames for astrodynamics applications. AstroFrames provides types for defining coordinate systems with customizable origins and axes and conversions between different coordinates.

# Quick Start

A coordinate system consists of an origin and axis system.  

```julia
using Epicycle

# Create an Earth-centered ICRF coordinate system
earth_icrf = CoordinateSystem(earth, ICRFAxes())

# Create a spacecraft-relative VNB coordinate system
sc = Spacecraft(name = "sat")
sat_vnb = CoordinateSystem(sc, VNB())
```

## Available Axes Types

To see all supported coordinate axes types, use Julia's type system:

```julia
# Find all concrete subtypes of AbstractAxes
subtypes(AbstractAxes)
4-element Vector{Any}:
 ICRFAxes
 Inertial
 MJ2000Axes
 VNB
```

