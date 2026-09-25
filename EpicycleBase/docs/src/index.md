```@meta
CurrentModule = EpicycleBase
```

# EpicycleBase

## Overview

The EpicycleBase package provides the abstract types and shared contracts used across Epicycle. It defines common categories for model variables, functions, calculated quantities, geometric points, and scalar variable identifiers so that higher-level packages can work with the same types without depending on one another. EpicycleBase contains no concrete astrodynamics models or user workflows; packages such as AstroStates, AstroEpochs, AstroFrames, AstroProp, and AstroSolve provide those implementations.

## Model Variables and Functions

`AbstractVar` is the common type for values that may vary in a model. Its subtypes separate state variables, control variables, time variables, and parameters as `AbstractState`, `AbstractControl`, `AbstractTime`, and `AbstractParam`.

`AbstractFun` is the common type for function objects used by Epicycle models. `AlgebraicFun` identifies functions that do not define differential equations.

## Calculated Quantities

`AbstractCalcVariable` is the common type for quantities calculated from model data. Its subtypes group orbital, celestial-body, and maneuver quantities as `AbstractOrbitVar`, `AbstractBodyVar`, and `AbstractManeuverVar`. `AbstractOrbitStateType` identifies orbital state representations within the orbital group.

## Geometric Points

`AbstractPoint` is the common type for objects that can serve as geometric points or coordinate-system origins. Celestial bodies and spacecraft are concrete examples defined by other Epicycle packages.

## Variable Identifiers

Variable tags identify scalar fields on model objects. `AbstractVarTag` is their common type, with `AbstractStateTag`, `AbstractParamTag`, `AbstractControlTag`, and `AbstractTimeTag` separating the corresponding variable categories.

`ModelVariable` pairs a model object with a variable tag, which identifies one scalar field on one model instance. `DirectVariable` represents a scalar decision variable that does not belong to a model object, such as a time of flight, scale factor, or maneuver component.

Higher-level packages define how tagged fields are read and written, how model and parameter Jacobians are evaluated, and which behaviors a calculated quantity supports. EpicycleBase provides the common function contracts for those definitions.

```@index
```

