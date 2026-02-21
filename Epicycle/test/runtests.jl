using Epicycle
using Test

@testset "Epicycle.jl" begin
    include("graphics/runtests_graphics.jl")
    # Visual regression tests excluded from CI - run manually before releases
    # include("graphics/visual_regression.jl")
end
