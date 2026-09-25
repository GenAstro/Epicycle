# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# A caption strip across the top of the dashboard, for words that belong beside the plots.
#
# A presenter talking over a live run, or a tutorial that explains each step as it draws, needs the
# explanation where the audience is already looking. The strip holds one moment at a time: a
# heading, the text, a few labelled values, and the code that produced what is on screen. It goes
# through `_publish` like a panel, so it is atomic and throttled the same way.
#
# A caption can belong to a set of panels. It is then also written under a name derived from
# those panels, and a tab scoped to them reads that copy instead of the latest caption, so the tab
# keeps its own words after the script has moved on to other panels.

const NARRATION_FILE = "narration.json"

const _NARRATION_SEQ = Ref(0)

"""The file a caption for `panels` is kept in. The dashboard derives the same name from its scope."""
narration_file(panels) =
    isempty(panels) ? NARRATION_FILE :
    "narration_" * join(sort(panel_id.(collect(panels))), "~") * ".json"

"""
    narrate(text = ""; title = "", subtitle = "", code = "", facts = Pair[], status = "",
            panels = String[])

Show a caption across the top of the dashboard.

Each call replaces what the strip shows, so a presenter's script calls it once per moment. The
strip appears when there is something to show and disappears after [`clear_narration`](@ref).

# Arguments
- `text`: The narration, as plain text. Line breaks are kept.
- `title`: A heading, such as the name of a section of a tour.
- `subtitle`: A line under the heading.
- `code`: Julia source to show beside the text, highlighted. Show the code that actually ran, so
  what a reader copies is what produced the plots.
- `facts`: Labelled values to set out as a small table, as `"label" => "value"` pairs. Values are
  shown as given, so format numbers and units before passing them.
- `status`: A short note on what is happening now, such as "solving…" or "done in 1.2 s".
- `panels`: The panel titles this caption describes. A tab opened with
  `open_dashboard(panels...)` for exactly those panels shows this caption and keeps it after
  later captions go to other panels. The full dashboard always shows the latest caption.

# Returns
`nothing`. The caption is published to the dashboard.

# Examples
```julia
narrate("Propagating a quarter of a day with the Moon and Sun as third bodies.";
        title  = "A spacecraft in orbit",
        code   = "propagate!(prop, sat, StopAt(sat, PropDurationDays(), 0.25))",
        facts  = ["inclination" => "51.6 deg", "period" => "92.6 min"],
        panels = ["Orbit", "Altitude"])

clear_narration()
```
"""
function narrate(text::AbstractString = ""; title::AbstractString = "",
                 subtitle::AbstractString = "", code::AbstractString = "",
                 facts = Pair[], status::AbstractString = "",
                 panels::AbstractVector{<:AbstractString} = String[])
    for f in facts
        f isa Pair || throw(ArgumentError(
            "narrate: facts must be \"label\" => \"value\" pairs; got $(repr(f))"))
    end
    ensure_server()
    _NARRATION_SEQ[] += 1
    payload = JSON.json(Dict(:seq      => _NARRATION_SEQ[],
                             :title    => String(title),
                             :subtitle => String(subtitle),
                             :text     => String(text),
                             :code     => String(code),
                             :facts    => [[string(first(f)), string(last(f))] for f in facts],
                             :status   => String(status)))
    _publish(NARRATION_FILE, payload)
    isempty(panels) || _publish(narration_file(panels), payload)
    return nothing
end

"""
    clear_narration()

Remove the caption strip from the dashboard, including captions kept for scoped tabs.

# Returns
`nothing`.

# Examples
```julia
clear_narration()
```
"""
function clear_narration()
    _SERVER[] === nothing && return nothing
    sweep_narration_files()
    _NARRATION_SEQ[] += 1
    _publish(NARRATION_FILE, JSON.json(Dict(:seq => _NARRATION_SEQ[])))
    return nothing
end

"""Delete the captions kept for scoped tabs. Like panel files, they outlive the process."""
function sweep_narration_files()
    dir = assets_dir()
    for f in readdir(dir)
        startswith(f, "narration_") && endswith(f, ".json") && rm(joinpath(dir, f); force = true)
    end
    return nothing
end
