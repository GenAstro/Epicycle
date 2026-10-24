# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# The examples are the interface's specification, so the suite runs each one: it must finish
# without a warning, and every plot panel it draws must hold data. What the pictures look like
# stays a visual check (manual_tests.md).

using Test
using EpicycleIO
using JSON
using Logging

const EXAMPLES_DIR = joinpath(dirname(@__DIR__), "examples")

# Whether a published trace carries any data, on whichever keys its type uses.
has_data(trace) = any(k -> haskey(trace, k) && !isempty(trace[k]), ("x", "y", "z", "r", "lon"))

@testset "examples run without a warning and draw data" begin
    EpicycleIO.auto_open!(false)
    files = sort(filter(f -> startswith(f, "Ex_") && endswith(f, ".jl"), readdir(EXAMPLES_DIR)))
    @test !isempty(files)
    for f in files
        @testset "$f" begin
            EpicycleIO.reset_panels!()
            logger = Test.TestLogger(min_level = Logging.Warn)
            with_logger(logger) do
                # Each in a module of its own, as a user's script runs: the names the other
                # test files leave in Main would otherwise collide with the example's.
                m = Module(Symbol(splitext(f)[1]))
                Core.eval(m, :(include(path) = Base.include($m, path)))   # as any module has
                Base.include(m, joinpath(EXAMPLES_DIR, f))
            end
            EpicycleIO.flush_publishes()
            @test isempty(logger.logs) || (@info "warnings in $f" [l.message for l in logger.logs]; false)
            @test !isempty(EpicycleIO.panels())
            for p in EpicycleIO.panels()
                p.family === :scene && continue                     # a Cesium view, checked in test_orbitview.jl
                data = JSON.parse(EpicycleIO.payload(p))["data"]
                @test !isempty(data) && any(has_data, data)
            end
        end
    end
    EpicycleIO.reset_panels!()
end
