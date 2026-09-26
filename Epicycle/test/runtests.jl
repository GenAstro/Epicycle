# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

using Epicycle
using Test

# The umbrella re-exports; it implements almost nothing, so what is worth testing here is that
# the re-export works and that a name from each layer arrives. The graphics tests that used to
# be the whole of this file were removed with the native graphics code on 2026-09-13.
#
# Two files in this directory are not run by anything: UseCase_Propagation_Advanced.jl and
# Benchmark_Propagation.jl. They predate this file and were never included.

include("test_correctness_examples.jl")
include("test_robustness_examples.jl")

@testset "Epicycle.jl" begin
    @testset "the umbrella re-exports each layer" begin
        for name in (:Spacecraft,          # AstroModels
                     :CartesianState,      # AstroStates
                     :Time,                # AstroEpochs
                     :earth,               # AstroUniverse
                     :propagate!,          # AstroProp
                     :Sequence,            # AstroSolve
                     :solve!)              # AstroSolve, resolved here deliberately
            @test isdefined(Epicycle, name)
        end
    end

    @testset "graphics is not here" begin
        # The native GLMakie window is out of the release. Re-exporting these would put GLMakie
        # back in the umbrella's dependency tree, which is the cost removing it bought.
        for name in (:View3D, :add_spacecraft!, :display_view)
            @test !isdefined(Epicycle, name)
        end
    end
end
