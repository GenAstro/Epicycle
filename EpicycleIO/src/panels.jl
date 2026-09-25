# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Panels are the other thing EpicycleIO adds to Plotly. Plotly has a figure; it has
# no notion of a long-lived named panel that a running solver keeps updating.
#
# A panel is named by a string. Plotting into a name that already exists replaces its contents,
# which is what makes a re-run script leave one figure rather than twenty, and is the same rule
# live update needs.

const DEFAULT_PANEL = "Plot"

# `layout` and `czml` hold whatever Plotly and Cesium accept, which is arbitrary nested JSON, so
# their value types stay open (§14.4). A panel is touched once per publish rather than per
# sample, so nothing here is specialization-sensitive.
mutable struct Panel
    id::String
    title::String
    traces::Vector{PlotlyBase.GenericTrace}
    layout::Dict{Symbol, Any}
    family::Symbol          # :cartesian, :polar, :geo, :scene — Plotly cannot mix these
    order::Int
    revision::Int           # bumped when the contents are replaced; see `figure`
    czml::Union{Nothing, Vector}   # a scene panel carries CZML instead of traces
end

const _PANELS = Dict{String, Panel}()
const _ORDER  = Ref(0)

"""Panel id: the title with anything awkward replaced, matching the CZML entity id rule."""
panel_id(title::AbstractString) = replace(String(title), r"[^A-Za-z0-9_\-]" => "_")

# Which coordinate system a trace type draws in. Plotly draws polar, geographic and
# three-dimensional traces on their own subplot objects, so they cannot share a panel with a
# cartesian one. This is Plotly's constraint reported, not a limit we impose.
function trace_family(t::Symbol)
    t in (:scatterpolar, :scatterpolargl, :barpolar) && return :polar
    t in (:scattergeo, :choropleth)                  && return :geo
    t in (:scatter3d, :surface, :mesh3d, :cone)      && return :scene
    return :cartesian
end

trace_type_of(tr::PlotlyBase.GenericTrace) = Symbol(get(getfield(tr, :fields), :type, "scatter"))

"""Fetch a panel, creating it if this is the first time the name has been used."""
function get_panel(title::AbstractString)
    id = panel_id(title)
    get!(_PANELS, id) do
        _ORDER[] += 1
        Panel(id, String(title), PlotlyBase.GenericTrace[], Dict{Symbol, Any}(),
              :cartesian, _ORDER[], 1, nothing)
    end
end

"""
Warn when a name that is meant to reach an existing panel does not.

Configuring a panel before plotting into it is legitimate, so this cannot be an error. A
mistyped name is the far more common case, and without this it silently makes a second, empty
panel and leaves the real one unchanged, so the settings appear to have been ignored.
"""
function warn_unknown_panel(title::AbstractString, verb::AbstractString)
    haskey(_PANELS, panel_id(title)) && return nothing
    isempty(_PANELS) && return nothing          # nothing to have meant instead
    known = join(("\"" * p.title * "\"" for p in panels()), ", ")
    @warn "EpicycleIO: $verb was given \"$title\", which is not a panel yet. If you meant an " *
          "existing one, check the spelling — there is $known. Otherwise this creates it, and " *
          "it stays empty until something is plotted into it."
    return nothing
end

"""Every panel, in the order they were first created."""
panels() = sort(collect(values(_PANELS)); by = p -> p.order)

function _check_family!(p::Panel, traces)
    isempty(traces) && return nothing
    fam = trace_family(trace_type_of(first(traces)))

    if isempty(p.traces)
        p.family = fam
        return nothing
    end

    fam == p.family && return nothing
    throw(ArgumentError(
        "panel \"$(p.title)\" holds $(p.family) traces and a $(fam) trace cannot join them. " *
        "Plotly draws these on different subplot types, so they cannot share a panel. " *
        "Give this one a panel of its own."))
end

"""
Put a 3D scene in a panel, replacing whatever was there.

A scene panel carries CZML rather than Plotly traces, and the dashboard hands it to Cesium
rather than to Plotly. Everything else about a panel works the same way: its name, its place in
the manifest, and replace-versus-add.
"""
function set_scene!(title::AbstractString, czml)
    p = get_panel(title)
    empty!(p.traces)
    p.family = :scene
    p.czml   = czml
    publish_panel(p)
    return p
end

"""Replace a panel's contents."""
function set_traces!(title::AbstractString, traces)
    p = get_panel(title)
    empty!(p.traces)
    p.revision += 1              # different data: the old zoom should not survive it
    p.czml   = nothing           # a panel that held a scene is a plot again now
    p.family = :cartesian
    _check_family!(p, traces)
    append!(p.traces, traces)
    publish_panel(p)
    return p
end

"""Add to a panel, keeping what is already there."""
function add_traces!(title::AbstractString, traces)
    p = get_panel(title)
    _check_family!(p, traces)
    append!(p.traces, traces)
    publish_panel(p)
    return p
end

"""
Merge layout settings into a panel.

Deeply, because the settings arrive already nested: `panel!(…; xaxis_title=…)` followed by
`panel!(…; xaxis_type=…)` are two `:xaxis` dictionaries, and a shallow merge would let the
second erase the first.
"""
function set_layout!(title::AbstractString, layout::Dict{Symbol, Any})
    warn_unknown_panel(title, "panel!")
    p = get_panel(title)
    _deepmerge!(p.layout, layout)
    publish_panel(p)
    return p
end

function _deepmerge!(into::AbstractDict, from::AbstractDict)
    for (k, v) in from
        if v isa AbstractDict && get(into, k, nothing) isa AbstractDict
            _deepmerge!(into[k], v)
        else
            into[k] = v
        end
    end
    return into
end

"""
    clear!(panel = "Plot")

Empty a panel but leave it on the dashboard.

# Arguments
- `panel`: the panel to empty. Without one you get the default panel.

# Notes
The layout survives: a panel set to a log axis is still a log axis after clearing, so refilling
it does not mean configuring it again.

Clearing a panel that does not exist does nothing, silently. Emptying a script's panels before
filling them, so a re-run starts fresh rather than accumulating, is an ordinary thing to write —
and the first time it runs, those panels are not there yet.

# Returns
`nothing`. The emptied panel is republished, so the browser follows.

# Examples
```julia
clear!("Altitude")
```
"""
function clear!(title::AbstractString = DEFAULT_PANEL)
    haskey(_PANELS, panel_id(title)) || return nothing
    p = get_panel(title)
    empty!(p.traces)
    p.czml   = nothing
    p.family = :cartesian
    publish_panel(p)
    return nothing
end

"""
    clear_all!()

Remove every panel from the dashboard.

[`clear!`](@ref) empties one panel and leaves it in place; this takes them all away, layout and
all, and the next plot starts from an empty page.

# Notes
The browser follows. Forgetting the panels here without saying so would leave the page showing
panels that no longer exist, and the next plot would appear among the ghosts.

Panel files left in the served directory are deleted too, all of them rather than only this
session's. The directory outlives the process, so otherwise it fills with orphans from every
earlier run.

# Returns
`nothing`.

# Examples
```julia
clear_all!()          # a script that starts from a clean dashboard
```
"""
clear_all!() = reset_panels!()

"""
Forget every panel, and tell the browser.

The manifest is republished and the panel files are deleted, because the dashboard decides what
to draw from the manifest and reads the data from those files. Clearing only the Julia side left
the page exactly as it was, which looked like the call had done nothing at all.
"""
function reset_panels!()
    stale = collect(values(_PANELS))
    empty!(_PANELS)
    _ORDER[] = 0

    # Nothing is serving, so there is nobody to tell.
    _SERVER[] === nothing && return nothing

    sweep_panel_files()
    sweep_narration_files()          # captions kept for scoped tabs describe panels now gone
    _LAST_MANIFEST[] = ""            # the empty manifest must go out, not be skipped as unchanged
    publish_manifest()
    return nothing
end

"""
Delete every panel file in the served directory.

All of them, not only the ones this session knows about. The directory outlives the process, so
without this it fills with orphans from every earlier run, and a panel file whose panel no
longer exists is a file the dashboard may still be asked for.
"""
function sweep_panel_files()
    dir = assets_dir()
    for f in readdir(dir)
        startswith(f, "panel_") && (endswith(f, ".json") || endswith(f, ".czml")) &&
            rm(joinpath(dir, f); force = true)
    end
    return nothing
end

# ─── What reaches the browser ─────────────────────────────────────────────────────────────────

"""
A panel as the Plotly figure the dashboard renders.

`uirevision` is what decides whether a zoom survives an update, and it is Plotly's own mechanism
rather than anything of ours. Holding it steady tells Plotly to keep what the user did to the
view; changing it tells Plotly to start again.

That is exactly the difference between the two verbs. `xyplot!` adds to a panel, so the view is
still about the same picture and a zoom must survive, which is the whole point of updating with
`react` while a solver runs. `xyplot` replaces the contents, and keeping the old axis range then
leaves the new data drawn somewhere off screen, in a panel that looks empty until autoscale is
pressed.
"""
function figure(p::Panel)
    layout = copy(p.layout)
    haskey(layout, :title) || (layout[:title] = p.title)
    haskey(layout, :uirevision) || (layout[:uirevision] = p.revision)
    return Dict(:data   => [getfield(t, :fields) for t in p.traces],
                :layout => layout)
end

"""
What a panel publishes: CZML for a 3D view, a Plotly figure for anything else.

The dashboard reads the panel's kind from the manifest and hands the file to Cesium or to
Plotly accordingly, so this is the only place the two diverge.
"""
payload(p::Panel) = p.family === :scene ? JSON.json(p.czml) : JSON.json(figure(p))

"""The manifest: which panels exist, in what order, and what kind each is."""
function manifest()
    return Dict(:panels => [Dict(:id     => p.id,
                                 :title  => p.title,
                                 :family => String(p.family),
                                 :file   => panel_file(p),
                                 :kind   => p.family === :cartesian ? "xy" : String(p.family))
                            for p in panels()])
end
