```@meta
CurrentModule = AstroSolve
```

# API reference

The names a user writes to set up a problem, solve it and read the answer. The interfaces a new
transcription, phase type or measurement is written against are not exported and are imported by
name, as in `using AstroSolve: jacobian_chunk`.

## The interfaces most problems use

```@docs
Vary
Constraint
Objective
Link
continuity
Sequence
solve!
Optimize
target!
Event
CollocationPhase
HermiteSimpson
ODProblem
simulate
Batch
Sequential
check_partials
```

## All exported names

```@index
```

```@autodocs
Modules = [AstroSolve, AstroSolve.Measurements, AstroSolve.TrackingDataIO,
           AstroSolve.ProcessNoiseModels, AstroSolve.BatchLeastSquares,
           AstroSolve.ExtendedKalmanFilter]
Order   = [:type, :function, :macro, :constant]
Public  = true
Private = false
Filter  = t -> !(t in (Vary, Constraint, Objective, Link, continuity, Sequence, solve!,
                       Optimize, target!, Event, CollocationPhase, HermiteSimpson,
                       ODProblem, simulate, Batch, Sequential, check_partials))
```
