```@meta
CurrentModule = AstroFun
```

# AstroFun

AstroFun provides fundamental calculations used throughout the Epicycle ecosystem. These calculations can be used in targeting and optimization variables and constraints, propagator stopping conditions, I/O operations, and eventally solver-for and consider parameters for orbit determination.

The module serves as a foundation for higher-level astrodynamics operations, providing core mathematical functions and utilities needed across multiple Epicycle packages.

## Calculation Framework

AstroFun contains the "Calcs" - a collection of structs that provide a unified interface for setting and getting quantities throughout Epicycle. These calculation objects can be used in stopping conditions, solver variables and constraints, and estimator solve-for and consider parameters.

The Calc framework provides type-stable access to commonly needed astrodynamics quantities while maintaining automatic differentiation compatibility. Each Calc struct implements standardized interfaces for both retrieving values from spacecraft states, celestial bodies, and setting target values for optimization.

### Calc Types

The framework includes several categories of calculations, and is designed for extensibility. A few types of Calcs include:

- **OrbitCalc**: Semi-major axis, eccentricity, inclination, incoming/outgoing asymptotes, periapsis conditions, and other orbital properties
- **BodyCalc**: Celestial body gravitational parameters, physical properties, and other body-specific quantities  
- **ManeuverCalc**: Δv components and magnitude, thrust direction, and other maneuver-related properties

Each Calc type supports both getter operations (extracting values) and setter operations (defining target values for optimization or constraints).

## Quick Start

Examples of common orbital calculations: (see AstroProp and AstroSolve docs for interation into those packages in stopping conditions, optimization variables, and constraints.)

```julia
using AstroFun, AstroStates, AstroModels

# Create a spacecraft with orbital state
sc = Spacecraft(CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]), 
                Time("2024-01-01T12:00:00"), 1000.0)

# Get semi-major axis from current state
ecc_calc = OrbitCalc(sc, SMA())
a = get_calc(ecc_calc)           
set_calc!(ecc_calc, 10000.0)  

# Set target incoming asymptote (rp = 6900, C3 = 14.0)
hyp = OrbitCalc(sc, IncomingAsymptote())
set_calc!(hyp, [6900.0, 14.0, 0.0, 0.0, 0.0, 0.0])  
    
# Set and get Earth's mu
mu_calc = BodyCalc(earth, GravParam())
μ = get_calc(mu_calc)            
set_calc!(mu_calc, 3.986e5)      

# Set and get maneuver elements
toi = ImpulsiveManeuver()
dvvec_calc = ManeuverCalc(toi, sc, DeltaVVector())
Δv = get_calc(dvvec_calc)   
set_calc!(dvvec_calc, [0.2, 0.3, 0.4])
```

## Table of Contents

```@index
```

## API Reference

```@autodocs
Modules = [AstroFun]
Order = [:type, :function, :macro, :constant]
```
