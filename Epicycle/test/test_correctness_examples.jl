# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Running a shipped example, `Epicycle.Examples`.
#
# The tests run a synthetic script written to a temporary file rather than a shipped one. A
# shipped example propagates or solves, which would make this suite minutes long and would test
# the physics again rather than the runner.
#
# Everything is reached through `Epicycle.Examples`. The four public names are deliberately not
# exported — `Epicycle.jl` says so where it imports them — so a test that called them unqualified
# would be testing a namespace the package does not offer.

using Epicycle
using Test

const EX = Epicycle.Examples

# A licence header to be dropped, Literate markers to be flattened, an ordinary comment to be left
# alone, and two assignments so a run leaves something behind.
const _SYNTHETIC = """
# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # Synthetic Example
#'
#' Prose that Literate would render on the page.

#' ## A Section

# Add two numbers
a = 2 + 3

# Scale the result
b = a * 4
"""

"Write `text` to a temporary `.jl` file, hand the path to `f`, and remove it afterwards."
function with_example(f, text = _SYNTHETIC)
    dir = mktempdir()
    try
        path = joinpath(dir, "Ex_Synthetic.jl")
        write(path, text)
        f(path)
    finally
        rm(dir; recursive = true, force = true)
    end
end

"""
Run `f` with stdout captured, and return `(value, text)`.

`redirect_stdout` needs a real stream rather than an `IOBuffer`, so this goes through a
temporary file.
"""
function captured(f)
    sink, io = mktemp()
    try
        value = redirect_stdout(io) do
            f()
        end
        close(io)
        return value, read(sink, String)
    finally
        isopen(io) && close(io)
        rm(sink; force = true)
    end
end

@testset "Examples" begin

    @testset "examples_dir and example_names" begin
        dir = EX.examples_dir()
        @test isdir(dir)
        @test basename(dir) == "examples"

        names = EX.example_names()
        @test !isempty(names)
        @test names == sort(names)
        @test all(!endswith(n, ".jl") for n in names)
        # Every name listed resolves to a file, which is what the other three rely on.
        @test all(isfile(joinpath(dir, n * ".jl")) for n in names)
    end

    @testset "list_examples prints every name and returns nothing" begin
        # The reason this exists: the REPL abbreviates a 32-element vector, so `example_names()`
        # shows a reader perhaps a dozen of them and a vertical ellipsis. Printing has to be
        # complete, and it has to return nothing, or the abbreviated vector comes straight back.
        value, text = captured() do
            EX.list_examples()
        end
        @test value === nothing

        names = EX.example_names()
        for n in names
            @test occursin(n, text)
        end
        # One per line, and no more lines than there are names.
        @test length(split(strip(text), '\n')) == length(names)
        @test !occursin("⋮", text)      # the ellipsis the REPL would have inserted
    end

    @testset "example_path resolves a name, a name with .jl, and a path" begin
        expected = joinpath(EX.examples_dir(), "Ex_GettingStarted.jl")
        @test EX.example_path("Ex_GettingStarted") == abspath(expected)
        @test EX.example_path("Ex_GettingStarted.jl") == abspath(expected)
        @test isabspath(EX.example_path("Ex_GettingStarted"))

        with_example() do path
            @test EX.example_path(path) == abspath(path)
        end
    end

    @testset "example_source drops the licence and flattens the Literate markers" begin
        with_example() do path
            src = EX.example_source(path)

            @test !occursin("SPDX-License-Identifier", src)
            @test !occursin("Copyright", src)

            # `#'` is a build marker, so the echo shows the prose as an ordinary comment.
            @test !occursin("#'", src)
            @test occursin("# Synthetic Example", src)
            @test occursin("# Prose that Literate would render on the page.", src)
            @test occursin("# A Section", src)
            # A `#' ##` heading flattens to one `#`, not two.
            @test !occursin("# ##", src)
            @test !occursin("# # ", src)

            # An ordinary comment and the code are untouched.
            @test occursin("# Add two numbers", src)
            @test occursin("a = 2 + 3", src)
            @test occursin("b = a * 4", src)
        end
    end

    @testset "what example_source prints still runs" begin
        with_example() do path
            mod = Module(:SourceRunsProbe)
            Base.include_string(mod, EX.example_source(path))
            @test Base.invokelatest(getfield, mod, :b) == 20
        end
    end

    @testset "run_example runs the script and leaves its results behind" begin
        with_example() do path
            mod = Module(:RunProbe)
            result, _ = captured() do
                EX.run_example(path; mod = mod)
            end
            @test result === nothing
            @test Base.invokelatest(getfield, mod, :a) == 5
            @test Base.invokelatest(getfield, mod, :b) == 20
        end
    end

    @testset "run_example prints the script, then the path" begin
        with_example() do path
            mod = Module(:EchoProbe)
            _, text = captured() do
                EX.run_example(path; mod = mod)
            end
            @test occursin("a = 2 + 3", text)
            @test occursin("# Synthetic Example", text)
            @test occursin(path, text)
            # The script comes before the path, so a reader sees the code, then the results, then
            # where to find the file.
            @test findfirst("a = 2 + 3", text)[1] < findfirst(path, text)[1]
        end
    end

    @testset "echo = false prints neither" begin
        with_example() do path
            mod = Module(:QuietProbe)
            _, text = captured() do
                EX.run_example(path; mod = mod, echo = false)
            end
            @test !occursin("a = 2 + 3", text)
            @test !occursin("Script:", text)
            @test Base.invokelatest(getfield, mod, :b) == 20
        end
    end

    @testset "a name with or without .jl runs the same script" begin
        with_example() do path
            for given in (path, first(path, length(path) - 3))
                mod = Module(Symbol("SuffixProbe", hash(given)))
                captured() do
                    EX.run_example(given; mod = mod, echo = false)
                end
                @test Base.invokelatest(getfield, mod, :a) == 5
            end
        end
    end
end
