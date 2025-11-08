


using LinearAlgebra
using AstroFrames
using AstroProp
using AstroEpochs
using AstroStates
using AstroMan
using AstroSolve
using AstroUniverse

include("config_hohmann_transfer.jl")

println("States immediately after initialization of sm")
println(sm.initial_stateful_structs[3].state.state)
println(sat1.state.state)

# Set up x to match toi2 and moi2
x = [ 0.15, 0.25, 0.35, 0.45, 0.55, 0.65]

F = similar(get_fun_values(sm))

# Call solver_fun! (this sets vars and applies events to sat1)
solver_fun!(F, x, sm)
println("States immediately after call to solver_fun!")
println(sm.initial_stateful_structs[3].state)

