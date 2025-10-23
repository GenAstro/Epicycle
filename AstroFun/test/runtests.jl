using Test
using LinearAlgebra

using AstroStates
using AstroEpochs
using AstroBase
using AstroUniverse
using AstroCoords
using AstroModels: Spacecraft, get_state 
using AstroMan
using AstroFun 

include("runtests_orbitcalcs.jl")
include("runtests_constraint.jl")
include("runtests_maneuvercalcs.jl")
include("runtests_bodycalcs.jl")
include("runtests_inputvalidation.jl")
include("runtests_infrastructure.jl")

nothing