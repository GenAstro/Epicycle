# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Running a shipped example from the REPL.
#
# `run_example` prints the script, runs it, and says where the file is. Output on its own reports
# what Epicycle computed and nothing about how it was asked, and the how is the part worth seeing;
# the path at the end is what turns a run into something the reader can open and change.
#
# Code is evaluated in `Main` by default, so when an example ends the spacecraft, the propagator
# and whatever it solved for are at the prompt.
#
# The files are never modified, and none of this is needed to run one: they are ordinary scripts
# in `Epicycle/examples/`, and `include` works as well as this does.

module Examples

"""The directory the shipped examples live in."""
examples_dir() = joinpath(dirname(@__DIR__), "examples")

# Two lines that say nothing about the example, dropped from the echo.
_is_licence(line) = occursin(r"^#\s*(Copyright \(C\)|SPDX-License-Identifier)", lstrip(line))

# `#'` marks the prose Literate renders on the documentation page, and `#' ##` its headings.
# Both are build markers rather than anything a reader of the script needs, so the echo shows
# them as the ordinary comments they read as.
_plain_comment(line) = replace(line, r"^(\s*)#'[ \t]?#*[ \t]?" => s"\1# ")

"""
    example_names() -> Vector{String}

The names of the examples that ship with Epicycle, in alphabetical order. Each one can be handed
to [`run_example`](@ref), [`example_source`](@ref) or [`example_path`](@ref).

# Notes
This is the programmatic form, for filtering or iterating. The REPL abbreviates a long vector,
so [`list_examples`](@ref) is what to call to read the whole set.

# Returns
A vector of names without the `.jl` extension.

# Example
```julia
Epicycle.example_names()
```
"""
example_names() = sort([first(f, length(f) - 3) for f in readdir(examples_dir())
                        if endswith(f, ".jl")])

"""
    list_examples()

Print the name of every shipped example, one per line.

# Notes
The REPL abbreviates a vector longer than the window with a `⋮`, which hides most of the set from
a reader who does not already know how to page through an array. This prints all of them and
returns nothing, so nothing is abbreviated on the way back out either.

# Returns
`nothing`. [`example_names`](@ref) returns the same names as a vector.

# Example
```julia
Epicycle.list_examples()
```
"""
function list_examples()
    for name in example_names()
        println("  ", name)
    end
    return nothing
end

"""
    example_path(name) -> String

Where an example's script is on disk, so it can be opened and changed. `name` comes from
[`example_names`](@ref).

# Arguments
- `name::AbstractString`: the example, with or without `.jl`. A path to an existing file is taken
  as given, so a script of your own is accepted wherever a shipped name is.

# Notes
Raises `ArgumentError` for a name that resolves to no file, naming it and pointing at
`example_names()`.

# Returns
An absolute path to a `.jl` file.

# Example
```julia
Epicycle.example_path("Ex_GettingStarted")
```
"""
function example_path(name::AbstractString)
    file = endswith(name, ".jl") ? name : name * ".jl"
    isfile(file) && return abspath(file)
    path = joinpath(examples_dir(), file)
    isfile(path) || throw(ArgumentError(
        "no example called $(repr(name)). `Epicycle.example_names()` lists the " *
        "$(length(example_names())) that ship with Epicycle."))
    return abspath(path)
end

"""
    example_source(name) -> String

The text of an example as a reader wants it: the licence header dropped, and the `#'` markers
Literate uses for the documentation page shown as ordinary comments. `name` comes from
[`example_names`](@ref).

# Arguments
- `name::AbstractString`: the example, with or without `.jl`, or a path to a script.

# Notes
The code is untouched. Only comment markers change, so what is printed still runs.

# Returns
The script as a `String`.

# Example
```julia
print(Epicycle.example_source("Ex_GettingStarted"))
```
"""
example_source(name::AbstractString) =
    join((_plain_comment(l) for l in eachline(example_path(name)) if !_is_licence(l)), "\n")

"""
    run_example(name; mod = Main, echo = true)

Print a shipped example, run it, and say where the script is. `name` comes from
[`example_names`](@ref).

# Arguments
- `name::AbstractString`: the example, with or without `.jl`, or a path to a script.
- `mod::Module`: where the script is evaluated, `Main` by default so what it built is at the
  prompt when it ends.
- `echo::Bool`: print the script and its path around the run. On by default, since results alone
  show what was computed and not what was asked.

# Notes
Raises `ArgumentError` for a name that resolves to no file. An error inside the example
propagates rather than being caught, so a run stops where it broke.

# Returns
`nothing`. Everything the example built is left in `mod`.

# Example
```julia
Epicycle.run_example("Ex_GettingStarted")
```
"""
function run_example(name::AbstractString; mod::Module = Main, echo::Bool = true)
    path = example_path(name)
    echo && println("\n", example_source(name), "\n")
    Base.include(mod, path)
    echo && println("\nScript: ", path)
    return nothing
end

end # module Examples
