# AstroBase.jl

*Foundation package for space systems analysis and astrodynamics*

## Overview

AstroBase provides the core abstract type hierarchy that forms the foundation for the entire Epicycle space systems application. It defines the fundamental types for variables, states, controls, time systems, functions, and geometric points that are used throughout the higher-level packages.

## Installation

AstroBase is part of the Epicycle monorepo. To use it:

```julia
using Pkg
Pkg.develop(url="https://github.com/GenAstro/Epicycle", subdir="AstroBase")
using AstroBase
```

Or from within the Epicycle monorepo:

```julia
using AstroBase
```

## Type Hierarchy

AstroBase defines several key abstract type hierarchies:

### Variables and Parameters
- `AbstractVar` - Base type for all variable kinds
- `AbstractState` - State variables (position, velocity, etc.)
- `AbstractControl` - Control variables (thrust, torque, etc.)
- `AbstractTime` - Time-related variables
- `AbstractParam` - Parameters and constants

### Functions
- `AbstractFun` - Base type for function objects
- `AlgebraicFun` - Algebraic function representations

### Calculation Variables
- `AbstractCalcVariable` - Base for calculated quantities
- `AbstractOrbitVar` - Orbital element calculations
- `AbstractBodyVar` - Celestial body properties
- `AbstractManeuverVar` - Maneuver-related calculations

### Geometric Types
- `AbstractPoint` - Geometric points in space
- `AbstractOrbitStateType` - Orbit state representations

## Usage Example

```julia
using AstroBase

# Define custom types using AstroBase foundations
struct MyOrbitState <: AbstractState
    position::Vector{Float64}
    velocity::Vector{Float64}
end

struct MySpacecraftPoint <: AbstractPoint
    coordinates::Vector{Float64}
    frame::String
end

# Type checking works as expected
@assert MyOrbitState <: AbstractVar
@assert MySpacecraftPoint <: AbstractPoint
```

## Philosophy

AstroBase follows the principle that "good abstractions enable powerful compositions." By providing a solid foundation of abstract types, the higher-level packages in Epicycle can:

- Share common interfaces
- Ensure type safety across the system
- Enable generic algorithms that work with multiple concrete types
- Maintain consistency in the overall architecture

## Integration with Epicycle

AstroBase is designed to be the foundation layer for:

- **AstroStates** - Concrete state vector implementations
- **AstroEpochs** - Time system implementations
- **AstroFrames** - Coordinate system types
- **AstroModels** - Mathematical model types
- **AstroProp** - Propagator implementations
- **AstroSolve** - Solver algorithm types

All higher-level packages in Epicycle build upon the abstract types defined here.

## Contributing

AstroBase is part of the larger Epicycle project. See the main repository for contribution guidelines.

## License

Copyright (C) 2025 Gen Astro LLC  
*(License to be determined)*