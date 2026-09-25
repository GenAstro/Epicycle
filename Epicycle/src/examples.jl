# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Running a shipped example as a guided REPL session.
#
# An example already carries what this needs. `#'` marks the prose Literate puts on the
# documentation page, an ordinary `#` labels a step inside the code, and a blank line separates one
# step from the next. That is a walkthrough: prose, then a labelled block, then what it printed.
#
# Code is evaluated in `Main`, so when a walkthrough ends the spacecraft, the phase and the result
# are at the prompt. Julia 1.12 redefines structs and constants, so a second example in the same
# session is fine; an object built before a redefinition prints with a `@world(...)` tag, which is
# Julia saying it belongs to the superseded definition.
#
# The example files are never modified, and none of this is needed to run one normally.

module Examples

const _ACCENT = :cyan
const _PROMPT = :green
const _VALUE  = :yellow
const _DIM    = :light_black

# What a block's code and its result line up under, which is the width of "  julia> ".
const _GUTTER = " "^9

"""The directory the shipped examples live in."""
examples_dir() = joinpath(dirname(@__DIR__), "examples")

"""One step of a walkthrough: the prose above a block, and the block."""
struct Step
    heading::String
    prose::String
    code::String
end

# ── Parsing ──────────────────────────────────────────────────────────────────────────────────

_is_prose(line) = startswith(lstrip(line), "#'")
_is_licence(line) = occursin(r"^#\s*(Copyright \(C\)|SPDX-License-Identifier)", lstrip(line))
_prose_text(line) = strip(replace(lstrip(line), r"^#'\s?" => ""))

"""Split a script into steps: a run of `#'` lines opens one, the code up to the next run is its."""
function _parse(path::AbstractString)
    steps = Step[]
    heading, prose, code = "", String[], String[]
    in_prose = false

    function flush!()
        if !isempty(prose) || !isempty(code)
            push!(steps, Step(heading, join(prose, "\n"), rstrip(join(code, "\n"))))
        end
        heading, prose, code = "", String[], String[]
    end

    for line in readlines(path)
        if _is_licence(line)
            continue
        elseif _is_prose(line)
            in_prose || isempty(code) || flush!()
            in_prose = true
            text = _prose_text(line)
            if startswith(text, "## ")
                heading = strip(text[4:end])
            elseif startswith(text, "# ")
                heading = strip(text[3:end])
            elseif !isempty(text)
                push!(prose, text)
            end
        else
            in_prose = false
            push!(code, line)
        end
    end
    flush!()
    return steps
end

"""Whether `str` is one or more complete top-level expressions."""
function _complete(str::AbstractString)
    try
        pos = 1
        while pos <= lastindex(str)
            expr, pos = Meta.parse(str, pos; greedy = true, raise = false)
            expr === nothing && break
            Meta.isexpr(expr, :incomplete) && return false
        end
        return true
    catch
        return false
    end
end

"""Split a step's code into the blocks a reader would paste one at a time."""
function _blocks(code::AbstractString)
    raw, current = String[], String[]
    for line in split(code, '\n')
        if isempty(strip(line)) && !isempty(current)
            push!(raw, join(current, "\n")); empty!(current)
        elseif !isempty(strip(line))
            push!(current, line)
        end
    end
    isempty(current) || push!(raw, join(current, "\n"))

    merged, pending = String[], ""
    for b in raw
        candidate = isempty(pending) ? b : pending * "\n\n" * b
        if _complete(candidate)
            push!(merged, candidate); pending = ""
        else
            pending = candidate
        end
    end
    isempty(pending) || push!(merged, pending)
    return merged
end

# ── Display ──────────────────────────────────────────────────────────────────────────────────

function _type(text; color = :normal, pace = 1.0)
    if pace <= 0
        printstyled(text; color = color); return
    end
    for c in text
        printstyled(c; color = color)
        c == ' ' || sleep(0.005 * pace)
    end
end

function _wait(seconds, pace, enter)
    if enter
        printstyled("\n    ⏎  Enter to continue"; color = _DIM)
        readline()
    elseif pace > 0
        sleep(seconds * pace)
    end
end

"""Show a block the way a user would type it, run it, and print what the REPL would print."""
function _run_block(src, mod, pace, max_lines)
    printstyled("\n  julia> "; color = _PROMPT, bold = true)
    _type(replace(src, "\n" => "\n" * _GUTTER); color = :white, pace = pace)
    println()
    value = Base.include_string(mod, src)
    if value !== nothing
        # The REPL limits long output but does not compact it. Compacting rounds a Julian Date to
        # six figures, which makes epochs in different time scales print identically.
        shown = rstrip(sprint(show, MIME"text/plain"(), value;
                              context = (:limit => true, :color => false)))
        lines = split(shown, '\n')
        # A blank line and the prompt's own indent separate the result from the call, so a long
        # block and what it returned do not read as one wall of text.
        println()
        for line in lines[1:min(end, max_lines)]
            printstyled(_GUTTER, line, "\n"; color = _VALUE)
        end
        length(lines) > max_lines &&
            printstyled(_GUTTER, "… ", length(lines) - max_lines, " more lines\n"; color = _DIM)
    end
    return value
end

"""Names this example defines that `mod` already has, so a reader is not surprised by a @world tag."""
function _shadowed(steps, mod)
    found = String[]
    for s in steps
        for m in eachmatch(r"(?m)^(?:mutable struct|struct|const)\s+([A-Za-z_][A-Za-z0-9_!]*)", s.code)
            name = Symbol(m.captures[1])
            isdefined(mod, name) && push!(found, String(name))
        end
    end
    return unique(found)
end

# ── The interface ────────────────────────────────────────────────────────────────────────────

"""The packages that carry a tutorial, in the order a reader would meet them."""
const TUTORIAL_PACKAGES = ("AstroEpochs", "AstroStates", "AstroUniverse", "AstroFrames",
                           "AstroModels", "AstroManeuvers", "AstroProp")

function _tutorial_path(pkg::AbstractString)
    src = Base.find_package(pkg)
    src === nothing && return nothing
    path = joinpath(dirname(dirname(src)), "examples", "Tutorial.jl")
    return isfile(path) ? path : nothing
end

"""
    tutorial_names() -> Vector{String}

The packages that have a tutorial written so far. Each one can be handed to
[`run_tutorial`](@ref).

# Returns
A vector of package names.

# Example
```julia
tutorial_names()
```
"""
tutorial_names() = [p for p in TUTORIAL_PACKAGES if _tutorial_path(p) !== nothing]

"""
    run_tutorial(package; pace = 1.0, enter = true, mod = Main)

Walk a package's tutorial a step at a time, the way [`run_example`](@ref) walks an example. A
tutorial covers one package on its own, where an example shows several working together.

# Returns
`nothing`. Everything the tutorial built is left in `mod`.

# Example
```julia
run_tutorial("AstroEpochs"; pace = 0, enter = false)
```
"""
function run_tutorial(pkg::AbstractString; kwargs...)
    path = _tutorial_path(pkg)
    path === nothing && throw(ArgumentError(
        "no tutorial for $(repr(pkg)). `tutorial_names()` lists the packages that have one: " *
        join(tutorial_names(), ", ")))
    return run_example(path; kwargs...)
end

"""
    example_names() -> Vector{String}

The names of the examples that ship with Epicycle, in alphabetical order. Each one can be handed
to [`run_example`](@ref).

# Returns
A vector of names without the `.jl` extension.

# Example
```julia
example_names()
```
"""
example_names() = sort([first(f, length(f) - 3) for f in readdir(examples_dir())
                        if endswith(f, ".jl")])

"""
    run_example(name; pace = 1.0, enter = true, mod = Main)

Walk a shipped example a step at a time: the prose, then each block of code typed and run, then
what it printed. `name` comes from [`example_names`](@ref).

# Arguments
- `name::AbstractString`: the example to run, with or without `.jl`.
- `pace::Real`: scales the typing and the pauses. `0` removes both.
- `enter::Bool`: wait for a keypress between steps rather than timing them.
- `mod::Module`: where the code is evaluated, `Main` by default so the results are at the prompt
  when the walkthrough ends.

# Returns
`nothing`. Everything the example built is left in `mod`.

# Example
```julia
run_example("Ex_HohmannTransfer"; pace = 0, enter = false)
```
"""
function run_example(name::AbstractString; pace::Real = 1.0, enter::Bool = true,
                     mod::Module = Main, max_result_lines::Int = 12)
    file = endswith(name, ".jl") ? name : name * ".jl"
    path = isfile(file) ? file : joinpath(examples_dir(), file)
    isfile(path) || throw(ArgumentError(
        "no example called $(repr(name)). `example_names()` lists the $(length(example_names())) " *
        "that ship with Epicycle."))

    steps = _parse(path)
    isempty(steps) && throw(ArgumentError("$path carries no `#'` prose, so there is nothing to walk"))

    title = isempty(first(steps).heading) ? basename(path) : first(steps).heading
    println()
    printstyled("  ", title, "\n"; color = _ACCENT, bold = true)
    printstyled("  ", "─"^max(length(title), 20), "\n"; color = _DIM)
    printstyled("  ", basename(path), "\n"; color = _DIM)

    shadowed = _shadowed(steps, mod)
    if !isempty(shadowed)
        printstyled("\n  Redefines ", join(shadowed, ", "), ", which this session already has.\n";
                    color = _DIM)
        printstyled("  Julia allows that. Objects made earlier print with a @world tag.\n";
                    color = _DIM)
    end

    for (i, s) in enumerate(steps)
        if !isempty(s.heading) && i > 1
            printstyled("\n  ", s.heading, "\n"; color = _ACCENT, bold = true)
        end
        if !isempty(s.prose)
            println()
            _type("  " * replace(s.prose, "\n" => "\n  "); color = :white, pace = pace)
            println()
        end
        for block in _blocks(s.code)
            _run_block(block, mod, pace, max_result_lines)
        end
        isempty(strip(s.code)) || _wait(1.2, pace, enter)
    end

    printstyled("\n  Done. Everything it built is at the prompt.\n"; color = _PROMPT)
    return nothing
end

end # module Examples
