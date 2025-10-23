using LinearAlgebra
#using OrdinaryDiffEq

using AstroBase
using AstroStates
using AstroEpochs
using AstroUniverse
using AstroCoords
using AstroModels
using AstroMan
using AstroFun
using AstroProp
using AstroSolve

# Create a spacecraft posvel SolverVariable
time=Time("2015-09-21T12:23:12", TAI(), ISOT())
posvel = [7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]
sat = Spacecraft(
    state=CartesianState(posvel), 
    time=Time("2015-09-21T12:23:12", TAI(), ISOT())
    )

using Test

include("runtests_construction_solvervariable.jl")
include("runtests_setget_solvervariable.jl")
include("runtests_events.jl")
include("runtests_events_deltav.jl")
include("runtests_sequence.jl")
include("runtests_sequence_manager.jl")
include("runtests_gaps.jl")
include("runtests_optimize.jl")



