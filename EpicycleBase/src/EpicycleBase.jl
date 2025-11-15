# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

__precompile__()

"""
    module EpicycleBase

Core abstract types shared across Epicycle (variables, states, controls, time, functions, points).
These form the public type hierarchy used by higher-level packages (AstroStates, AstroEpochs, AstroFrames,
AstroProp, AstroSolve, …).
"""
module EpicycleBase

export AbstractVar, AbstractState, AbstractControl, AbstractTime, AbstractParam
export AbstractFun, AlgebraicFun
export AbstractCalcVariable, AbstractOrbitVar, AbstractBodyVar, AbstractManeuverVar
export AbstractOrbitStateType
export AbstractPoint

"""
    AbstractVar

Base tag for all variable kinds (states, controls, time, parameters).

# Notes:
- Serves as the common supertype for variable categories.
"""
abstract type AbstractVar end

"""
    AbstractState <: AbstractVar

Base type for all state variables.

"""
abstract type AbstractState <: AbstractVar end

"""
    AbstractControl <: AbstractVar

Base type for all control variables.
"""
abstract type AbstractControl <: AbstractVar end

"""
    AbstractTime <: AbstractVar

Base type for time variables.
"""
abstract type AbstractTime <: AbstractVar end

"""
    AbstractParam <: AbstractVar

Base type for parameter variables.
"""
abstract type AbstractParam <: AbstractVar end

"""
    AbstractFun

Base type for function objects (e.g., dynamics, outputs).
"""
abstract type AbstractFun end

"""
    AlgebraicFun <: AbstractFun

Base type for algebraic (non-differential) function objects.
"""
abstract type AlgebraicFun <: AbstractFun end

"""
    AbstractCalcVariable

Base type for calculation variables (orbit, body, maneuver).
"""
abstract type AbstractCalcVariable end

"""
    AbstractOrbitVar <: AbstractCalcVariable

Base type for orbit calculation variables.
"""
abstract type AbstractOrbitVar    <: AbstractCalcVariable end

"""
    AbstractBodyVar <: AbstractCalcVariable 

Base type for body calculation variables.
"""
abstract type AbstractBodyVar     <: AbstractCalcVariable end

"""
    AbstractManeuverVar <: AbstractCalcVariable 

Base type for maneuver calculation variables.
"""
abstract type AbstractManeuverVar <: AbstractCalcVariable end

"""
    AbstractOrbitStateType <: AbstractOrbitVar

Base type for orbit state representation types (e.g., `Cartesian()`, `Keplerian()`, etc.).
"""
abstract type AbstractOrbitStateType <: AbstractOrbitVar end  

"""
    AbstractPoint

Base type for geometric points (e.g., Spacecraft, CelestialBody).
"""
abstract type AbstractPoint end

"""
    no_op()

A no-op function that returns nothing. Useful as a default callback placeholder.

# Notes:
- Unexported utility; reference as EpicycleBase.no_op.
"""
function no_op()
    return nothing
end

end