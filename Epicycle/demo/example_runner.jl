# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Run an example script as a guided REPL session.
#
# An example already carries everything this needs. `#'` marks the prose Literate puts on the
# documentation page, an ordinary `#` labels a step inside the code, and a blank line separates one
# step from the next. That is a walkthrough: prose, then a labelled block, then what it printed.
#
#   include(joinpath(@__DIR__, "example_runner.jl"))
#   ExampleRunner.run_example("Ex_OrbitDetermination")
#
# Code is evaluated into `Main`, so when the walkthrough ends the spacecraft, the phase and the
# result are sitting at the prompt to be inspected. Julia 1.12 redefines structs and constants, so
# a second example in the same session is fine; objects built before a redefinition print with a
# `@world(...)` tag, which is Julia saying they belong to the superseded definition.
#
# The example files are never modified, and nothing here is needed to run one normally.

module ExampleRunner

# narration.jl expects its including module to have these already
using Printf
using EpicycleIO: narrate, clear_narration

include(joinpath(@__DIR__, "narration.jl"))

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")

"""One step of a walkthrough: the prose above a block, and the block."""
struct Step
    heading::String        # "" unless the prose opened a section
    prose::String
    code::String
end

# ═══════════════════════════════════════════════════════════════════════════════
# Parsing
# ═══════════════════════════════════════════════════════════════════════════════

is_prose(line) = startswith(lstrip(line), "#'")
is_licence(line) = occursin(r"^#\s*(Copyright \(C\)|SPDX-License-Identifier)", lstrip(line))

"""Strip the `#'` marker from a prose line."""
prose_text(line) = strip(replace(lstrip(line), r"^#'\s?" => ""))

"""
    example_title(path) -> (title, intro)

The `#' # Title` line and the prose under it, which is the example's own introduction.
"""
function example_title(steps)
    isempty(steps) && return ("", "")
    first(steps).heading, first(steps).prose
end

"""
    parse_example(path) -> Vector{Step}

Split a script into steps. A run of `#'` lines opens a step and everything up to the next run of
`#'` lines is its code. A blank line inside the code separates one labelled block from the next,
and blocks are merged until they parse, so a multi-line call stays whole.
"""
function parse_example(path::AbstractString)
    lines = readlines(path)
    steps = Step[]
    heading, prose, code = "", String[], String[]
    in_prose = false

    flush_step!() = begin
        if !isempty(prose) || !isempty(code)
            push!(steps, Step(heading, join(prose, "\n"), rstrip(join(code, "\n"))))
        end
        heading, prose, code = "", String[], String[]
    end

    for line in lines
        if is_licence(line)
            continue
        elseif is_prose(line)
            # A new run of prose after code closes the step before it
            in_prose || isempty(code) || flush_step!()
            in_prose = true
            text = prose_text(line)
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
    flush_step!()
    return steps
end

"""
    code_blocks(code) -> Vector{String}

Split a step's code into the blocks a reader would paste one at a time. A blank line ends a block,
and a block that does not parse on its own is joined to the next, so a call spanning several lines
or a `do` block stays in one piece.
"""
function code_blocks(code::AbstractString)
    blocks, current = String[], String[]
    for line in split(code, '\n')
        if isempty(strip(line)) && !isempty(current)
            push!(blocks, join(current, "\n"))
            empty!(current)
        elseif !isempty(strip(line))
            push!(current, line)
        end
    end
    isempty(current) || push!(blocks, join(current, "\n"))

    # Join anything that is not a complete expression to the block after it
    merged, pending = String[], ""
    for b in blocks
        candidate = isempty(pending) ? b : pending * "\n\n" * b
        if parses_completely(candidate)
            push!(merged, candidate)
            pending = ""
        else
            pending = candidate
        end
    end
    isempty(pending) || push!(merged, pending)
    return merged
end

"""Whether `str` is one or more complete top-level expressions."""
function parses_completely(str::AbstractString)
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

# ═══════════════════════════════════════════════════════════════════════════════
# Running
# ═══════════════════════════════════════════════════════════════════════════════

"""Names this example defines that are already bound in `mod`, so a reader is not surprised."""
function shadowed_names(steps, mod::Module)
    found = String[]
    for s in steps, m in eachmatch(r"(?m)^(?:mutable struct|struct|const)\s+([A-Za-z_][A-Za-z0-9_!]*)",
                                   s.code)
        name = Symbol(m.captures[1])
        isdefined(mod, name) && push!(found, String(name))
    end
    return unique(found)
end

resolve(name::AbstractString) =
    endswith(name, ".jl") ? name : joinpath(EXAMPLES_DIR, name * ".jl")

"""
    run_example(name; pace = 1.0, enter = true, mod = Main)

Walk an example a step at a time. `name` is the file's name with or without `.jl`.

`pace` scales the typing and the pauses, and `pace = 0` removes both. `enter = true` waits for a
keypress between steps rather than timing them. Code is evaluated in `mod`, which defaults to
`Main` so the results are at the prompt afterwards.
"""
function run_example(name::AbstractString; pace::Real = 1.0, enter::Bool = true,
                     mod::Module = Main)
    path = resolve(name)
    isfile(path) || error("no example at $path")

    _PACE[] = pace
    _ENTER[] = enter
    steps = parse_example(path)
    isempty(steps) && error("$path has no `#'` prose, so there is nothing to walk")

    title, intro = example_title(steps)
    _title_card_for(title, basename(path))

    shadowed = shadowed_names(steps, mod)
    if !isempty(shadowed)
        _say("  This example redefines " * join(shadowed, ", ") *
             ", which this session already has."; color = _DIM)
        _say("  Julia allows that. Objects made before the change print with a @world tag.";
             color = _DIM)
        println()
    end

    isempty(intro) || (_say(intro; color = :white); println())
    _pause(1.5)

    act = 0
    for s in steps
        if !isempty(s.heading) && s.heading != title
            act += 1
            _act(act, s.heading, "")
            isempty(s.prose) || _say(s.prose; color = :white)
            println()
        elseif !isempty(s.prose) && s.heading == title
            # the introduction, already shown
        elseif !isempty(s.prose)
            _say(s.prose; color = :white)
            println()
        end

        for block in code_blocks(s.code)
            _repl(block; mod = mod, max_result_lines = 12)
        end
        isempty(strip(s.code)) || _pause()
    end

    printstyled("\n  Done. Everything the example built is at the prompt.\n"; color = _GOOD)
    return nothing
end

"""A heading for the walkthrough, in the demo's visual language."""
function _title_card_for(title, file)
    println()
    printstyled("  ", isempty(title) ? file : title, "\n"; color = _ACCENT, bold = true)
    printstyled("  ", "─"^max(length(isempty(title) ? file : title), 20), "\n"; color = _DIM)
    printstyled("  ", file, "\n\n"; color = _DIM)
end

end # module
