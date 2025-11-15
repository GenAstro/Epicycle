using Test
using LinearAlgebra

using AstroStates
using AstroEpochs
using EpicycleBase
using AstroUniverse
using AstroFrames
using AstroModels: Spacecraft, get_state 
using AstroManeuvers
using AstroCallbacks 

include("runtests_orbitcalcs.jl")
include("runtests_constraint.jl")
include("runtests_maneuvercalcs.jl")
include("runtests_bodycalcs.jl")
include("runtests_inputvalidation.jl")
include("runtests_infrastructure.jl")

nothing