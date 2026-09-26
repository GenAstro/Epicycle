# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# What `Epicycle.Examples` does with input it cannot use.
#
# This is among the first things a new user calls, and the most likely mistake is the name of an
# example. The message has to name what was asked for and say where the list is, because a bare
# `no method matching` at that moment reads as the package being broken.

using Epicycle
using Test

const EXR = Epicycle.Examples

@testset "Examples — robustness" begin

    @testset "an unknown name is refused by all three entry points" begin
        for f in (EXR.run_example, EXR.example_source, EXR.example_path)
            err = try
                f("Ex_NotAnExample")
                nothing
            catch e
                e
            end
            @test err isa ArgumentError
            @test occursin("Ex_NotAnExample", err.msg)
            @test occursin("example_names()", err.msg)
        end
    end

    @testset "an empty file runs and produces nothing" begin
        # Not an error. An empty script is a script; refusing it would be a rule with no purpose.
        dir = mktempdir()
        try
            path = joinpath(dir, "Ex_Empty.jl")
            write(path, "")
            mod = Module(:EmptyProbe)
            @test EXR.run_example(path; mod = mod, echo = false) === nothing
            @test EXR.example_source(path) == ""
        finally
            rm(dir; recursive = true, force = true)
        end
    end

    @testset "an error inside the script is not swallowed" begin
        # A run that fails should fail where it failed, rather than reporting success over it.
        dir = mktempdir()
        try
            path = joinpath(dir, "Ex_Throws.jl")
            write(path, "error(\"deliberate\")\n")
            mod = Module(:ThrowProbe)
            @test_throws Exception EXR.run_example(path; mod = mod, echo = false)
        finally
            rm(dir; recursive = true, force = true)
        end
    end

    @testset "a directory is not mistaken for an example" begin
        dir = mktempdir()
        try
            mkdir(joinpath(dir, "Ex_Directory.jl"))
            @test_throws ArgumentError EXR.example_path(joinpath(dir, "Ex_Directory.jl"))
        finally
            rm(dir; recursive = true, force = true)
        end
    end
end
