# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Narration shared by the demos: what the terminal prints, what the dashboard's caption strip
# shows, and how a step is timed and its source displayed.
#
# Each demo includes this file into its own module, so the settings and the caption belong to that
# demo. A demo supplies its own parts.
#
#   include(joinpath(@__DIR__, "narration.jl"))
#
# `_repl` is for a demo that shows what a user types: it types the line, evaluates it in the
# demo's module, and prints the value the way the REPL would.

# ═══════════════════════════════════════════════════════════════════════════════
# Narration
# ═══════════════════════════════════════════════════════════════════════════════

const _ACCENT = :cyan
const _GOOD   = :green
const _NUMBER = :yellow
const _DIM    = :light_black

# Settings for one run, set by `epicycle_demo` and read by the narration helpers.
const _PACE  = Ref(1.0)
const _ENTER = Ref(false)
const _TABS  = Ref(true)      # open a browser tab per part, scoped to that part's panels

# What the dashboard strip shows. A pause ends a paragraph, so the next line starts a new one.
mutable struct _Caption
    title::String
    subtitle::String
    lines::Vector{String}
    facts::Vector{Pair{String,String}}
    code::String
    status::String
    fresh::Bool
    panels::Vector{String}      # this part's panels, so its own tab keeps this caption
end
const _CAPTION = _Caption("", "", String[], Pair{String,String}[], "", "", false, String[])

function _publish_caption()
    c = _CAPTION
    EpicycleIO.narrate(join(c.lines, " "); title = c.title, subtitle = c.subtitle, code = c.code,
            facts = c.facts, status = c.status, panels = c.panels)
end

function _start_paragraph_if_fresh()
    if _CAPTION.fresh
        empty!(_CAPTION.lines)
        empty!(_CAPTION.facts)
        _CAPTION.fresh = false
    end
end

"""Sleep for `seconds` scaled by the pace. Pace zero skips every pause."""
_wait(seconds) = _PACE[] > 0 && sleep(seconds * _PACE[])

"""Say a line: typed into the terminal a character at a time, and added to the dashboard strip."""
function _say(text::AbstractString; color = :normal, delay = 0.012)
    _start_paragraph_if_fresh()
    push!(_CAPTION.lines, strip(text))
    _publish_caption()
    if _PACE[] <= 0
        printstyled(text, "\n"; color = color)
        return
    end
    for c in text
        printstyled(c; color = color)
        c == ' ' || sleep(delay * _PACE[])
    end
    println()
end

"""A labelled value, aligned in the terminal and listed in the strip."""
function _fact(label, value; unit = "")
    _start_paragraph_if_fresh()
    push!(_CAPTION.facts, String(label) => (isempty(unit) ? String(value) : "$value $unit"))
    _publish_caption()
    printstyled(@sprintf("    %-34s", label); color = _DIM)
    printstyled(value; color = _NUMBER, bold = true)
    isempty(unit) || printstyled(" ", unit; color = _DIM)
    println()
    _wait(0.35)
end

"""Hold the screen until the reader is ready: Enter in presenter mode, a timed pause otherwise."""
function _pause(seconds = 2.5)
    if _ENTER[]
        printstyled("\n    ⏎  press Enter to continue"; color = _DIM)
        readline()
    else
        _wait(seconds)
    end
    _CAPTION.fresh = true
end

function _act(number::Int, title::AbstractString, subtitle::AbstractString;
              panels::AbstractVector{<:AbstractString} = String[])
    c = _CAPTION
    c.title, c.subtitle, c.code, c.status = "$number · $title", subtitle, "", ""
    empty!(c.lines); empty!(c.facts); c.fresh = false
    c.panels = collect(String, panels)
    _publish_caption()

    println("\n")
    bar = "━"^78
    printstyled("  ", bar, "\n"; color = _ACCENT)
    printstyled(@sprintf("   %d  ", number); color = :black, bold = true, reverse = true)
    printstyled("  ", uppercase(title), "\n"; color = _ACCENT, bold = true)
    printstyled("   ", subtitle, "\n"; color = _DIM)
    printstyled("  ", bar, "\n\n"; color = _ACCENT)

    # A tab of its own keeps each part's graphics readable. The caption strip shows on every tab,
    # and the panels appear in it as the part draws them.
    if _TABS[] && !isempty(panels)
        open_dashboard(panels...)
        _wait(2.0)
    end
    _wait(1.2)
end

function _step_begin(label::AbstractString, code::AbstractString)
    _CAPTION.code   = code
    _CAPTION.status = "▸ $label …"
    _publish_caption()
    printstyled("  ▸ ", color = _ACCENT)
    printstyled(label, " … "; color = :normal)
    # Give the audience a moment to read the code before the result replaces what they look at.
    _wait(1.5)
    return time()
end

function _step_end(label::AbstractString, t0::Float64)
    dt = time() - t0
    _CAPTION.status = @sprintf("✔ %s · %.1f s", label, dt)
    _publish_caption()
    printstyled("done", color = _GOOD, bold = true)
    printstyled(@sprintf(" (%.1f s)\n", dt); color = _DIM)
end

"""
The lines of `file` between line `line`, which opens a block, and the `end` that closes it at the
same indentation, with the common indentation removed.
"""
function _source_block(file::AbstractString, line::Integer)
    lines  = readlines(file)
    indent = l -> length(l) - length(lstrip(l))
    open_indent = indent(lines[line])
    body = String[]
    for l in lines[line+1:end]
        strip(l) == "end" && indent(l) == open_indent && break
        push!(body, l)
    end
    nonblank = filter(l -> !isempty(strip(l)), body)
    d = isempty(nonblank) ? 0 : minimum(indent, nonblank)
    return join((isempty(strip(l)) ? "" : l[d+1:end] for l in body), "\n")
end

"""
    @step "label" begin ... end

Run a block as a narrated step: announce it, show its source on the dashboard, time it. The
block runs in the enclosing scope, so what it assigns is visible after it, and the step's value
is the block's value.
"""
macro step(label, block)
    code = _source_block(String(__source__.file), __source__.line)
    return quote
        local t0 = _step_begin($(esc(label)), $code)
        local value = $(esc(block))
        _step_end($(esc(label)), t0)
        value
    end
end

"""Hours from the first sample, for a plot axis.

Each epoch is put in TT first: a history can hold segments recorded in different time scales,
and `Time - Time` refuses a difference across them.
"""
function _hours(times)
    tt = [x.scale === :tt ? x : x.tt for x in times]
    return [((x.jd1 - tt[1].jd1) + (x.jd2 - tt[1].jd2)) * 24 for x in tt]
end

"""The Epicycle banner, then a subtitle and a tagline the demo chooses."""
function _title_card(subtitle::AbstractString = "A guided tour",
                     tagline::AbstractString =
                         "   Propagation · frames · targeting · estimation · optimal control · 3D graphics")
    art = raw"""
     ███████╗██████╗ ██╗ ██████╗██╗   ██╗ ██████╗██╗     ███████╗
     ██╔════╝██╔══██╗██║██╔════╝╚██╗ ██╔╝██╔════╝██║     ██╔════╝
     █████╗  ██████╔╝██║██║      ╚████╔╝ ██║     ██║     █████╗
     ██╔══╝  ██╔═══╝ ██║██║       ╚██╔╝  ██║     ██║     ██╔══╝
     ███████╗██║     ██║╚██████╗   ██║   ╚██████╗███████╗███████╗
     ╚══════╝╚═╝     ╚═╝ ╚═════╝   ╚═╝    ╚═════╝╚══════╝╚══════╝
    """
    colors = (:blue, :light_blue, :cyan, :light_cyan, :cyan, :light_blue)
    println()
    for (line, color) in zip(split(chomp(art), '\n'), colors)
        printstyled(line, "\n"; color = color, bold = true)
        _wait(0.08)
    end
    println()
    c = _CAPTION
    c.title, c.subtitle, c.code, c.status = "Epicycle", subtitle, "", ""
    empty!(c.lines); empty!(c.facts)
    _say("   Astrodynamics for mission design and navigation, written in Julia.";
         color = :light_cyan)
    _say(tagline; color = _DIM, delay = 0.006)
    println()
end

# ═══════════════════════════════════════════════════════════════════════════════
# What a user types
#
# A demo that teaches an interface shows the script, not a description of it. `_repl` types one
# line the way a user would, evaluates it in the demo's own module, and prints what the REPL would
# print. Later lines see what earlier ones defined, because each is evaluated at module scope.
# ═══════════════════════════════════════════════════════════════════════════════

"""Type `text` into the terminal a character at a time, without a newline."""
function _type(text::AbstractString; color = :normal, delay = 0.012)
    if _PACE[] <= 0
        printstyled(text; color = color)
        return
    end
    for c in text
        printstyled(c; color = color)
        c == ' ' || sleep(delay * _PACE[])
    end
end

"""
    @repl "kep = KeplerianState(cart, mu)"

Show a line of script and run it: the prompt and the line typed out, the value printed as the REPL
prints it, and both on the dashboard strip. Returns the value.

A comment-only line prints as a comment and returns `nothing`; `show_result = false` runs a line
whose value is noise.
"""
function _repl(src::AbstractString; show_result::Bool = true, mod::Module = @__MODULE__,
               max_result_lines::Int = typemax(Int))
    _CAPTION.code = isempty(_CAPTION.code) ? String(src) : _CAPTION.code * "
" * src
    _CAPTION.status = ""
    _publish_caption()
    printstyled("  julia> "; color = _GOOD, bold = true)
    _type(src; color = :white)
    println()

    value = Base.include_string(mod, src)
    if show_result && value !== nothing
        shown = rstrip(sprint(show, MIME"text/plain"(), value;
                              context = (:limit => true, :compact => true, :color => false)))
        lines = split(shown, '
')
        # A solver result prints pages of nested structs, which buries the walkthrough it belongs to
        clipped = length(lines) > max_result_lines
        for line in lines[1:min(end, max_result_lines)]
            printstyled("  ", line, "
"; color = _NUMBER)
        end
        clipped && printstyled("  … ", length(lines) - max_result_lines,
                               " more lines
"; color = _DIM)
        _CAPTION.status = "⇒ " * first(replace(shown, '
' => "  "), 120)
        _publish_caption()
    end
    _wait(0.9)
    return value
end

"""Start a fresh transcript: the strip shows only what follows."""
_clear_transcript() = (_CAPTION.code = ""; _CAPTION.status = ""; _publish_caption())
