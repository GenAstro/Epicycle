# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# What the user types. Everything here is either shape handling or panel management; attribute
# values pass through to Plotly untouched.

# A leading string names a panel. Peel it off before the arrays.
_split_panel(args) = (!isempty(args) && first(args) isa AbstractString) ?
                     (String(first(args)), Base.tail(args)) : (DEFAULT_PANEL, args)

# The palette, cycled. Shared with the Cesium side so a trajectory and its altitude plot come
# out the same colour without anyone asking.
const PALETTE = ["rgb(0,220,255)", "rgb(255,165,0)", "rgb(55,255,55)", "rgb(255,80,80)",
                 "rgb(140,100,255)", "rgb(255,220,0)", "rgb(255,120,200)", "rgb(55,200,255)"]

palette_color(i::Integer) = PALETTE[mod1(i, length(PALETTE))]

# Above this many points, a scatter is drawn by WebGL instead of SVG.
#
# SVG gives one DOM node per point, and a browser starts to stutter on pan and zoom somewhere in
# the low thousands; WebGL does not, but it renders slightly differently and is worth avoiding
# for the small plots that are most of them. A day of propagation at a 60 s step is about 1400
# points, so the ordinary case stays on SVG and a long or dense trajectory moves across on its
# own. Chosen for that boundary rather than measured — if a plot feels slow, this is the number
# to revisit, and the WebGL path has no polar equivalent, so only cartesian scatter is promoted.
const SCATTERGL_POINTS = 4_000

"""
Build the traces one call describes: expand the arrays into series, then give each series its
share of the attributes.
"""
function build_traces(trace_type::Symbol, args, kwargs; offset::Int = 0)
    series = expand_series(args...)
    n = length(series)

    effective = (trace_type === :scatter &&
                 maximum(length(s.x) for s in series) > SCATTERGL_POINTS) ?
                :scattergl : trace_type

    traces = PlotlyBase.GenericTrace[]
    for (i, s) in enumerate(series)
        fields = Dict{Symbol, Any}()
        for (k, v) in pairs(kwargs)
            k === :name && continue                       # handled below, it rewrites
            warn_unknown_attribute(effective, k)
            fields[k] = distribute(k, v, i, n, effective)
        end

        nm = series_name(get(kwargs, :name, nothing), s)
        nm === nothing || (fields[:name] = nm isa AbstractVector ?
                                           distribute(:name, nm, i, n, effective) : nm)

        haskey(fields, :line) || haskey(fields, :line_color) ||
            (fields[:line_color] = palette_color(offset + i))

        push!(traces, _trace(effective, s, fields))
    end
    return traces
end

# Which keys carry the data differs by trace type: cartesian uses x/y, polar theta/r,
# geographic lon/lat.
function _trace(trace_type::Symbol, s::Series, fields::Dict{Symbol, Any})
    fam = trace_family(trace_type)
    x, y = plotdata(s.x), plotdata(s.y)      # NaN and Inf become gaps rather than an error
    if fam === :polar
        return PlotlyBase.GenericTrace(trace_type; theta = x, r = y, fields...)
    elseif fam === :geo
        return PlotlyBase.GenericTrace(trace_type; lon = x, lat = y, fields...)
    else
        return PlotlyBase.GenericTrace(trace_type; x = x, y = y, fields...)
    end
end

# ─── plot and plot! ───────────────────────────────────────────────────────────────────────────

"""
    xyplot([panel], x, y; attributes...)
    xyplot([panel], y; attributes...)

Draw a line, replacing whatever was in the panel.

The leading string names a panel; without one you get the default panel. Plotting into a name
that already exists replaces it, so re-running a script gives you one figure rather than twenty.
Use [`xyplot!`](@ref) to add instead.

Attributes are Plotly's, passed through untouched:

```julia
xyplot(t, altitude; name = "altitude", mode = "lines", line_color = "cyan", line_width = 2)
```

A matrix or a vector of vectors draws one line per component, which is what a position column
from `history` already is:

```julia
xyplot(t, r; name = ["x", "y", "z"], line_color = ["red", "green", "blue"])
```

# Arguments
- `panel`: optional leading string naming the panel.
- `x`, `y`: paired arrays. Several pairs may follow one another, each with its own length.
  Given `y` alone, the horizontal axis is its index.

# Notes
`x` and `y` are plain arrays and carry no units, frame or time system of their own — whatever
you computed them in is what gets drawn, so label the axes with `panel!`. A `Time` is not a
number: render it first, or subtract the epoch to get hours.

A `NaN` or an `Inf` draws as a break in the line rather than raising.

# Returns
`nothing`. The panel is drawn on the dashboard, and a browser tab opens if none is watching.
"""
function xyplot(args...; kwargs...)
    title, rest = _split_panel(args)
    isempty(rest) && throw(ArgumentError(
        "plot needs data. Give it y, or x and y — for example xyplot(t, altitude). A leading " *
        "string names the panel: xyplot(\"Altitude\", t, altitude)."))
    traces = build_traces(:scatter, rest, kwargs)
    set_traces!(title, traces)
    ensure_browser()
    return nothing
end

"""
    xyplot!([panel], x, y; attributes...)

Add a line to a panel instead of replacing what is there. Takes everything [`xyplot`](@ref) does.

Lines need not be the same length, so two spacecraft sampled differently share a panel without
any fuss.

Safe to call from inside a propagation or a solver callback: writes are throttled and happen on
their own task, so publishing cannot slow the run feeding it. A zoom you set survives the
update, where [`xyplot`](@ref) resets the view because the data underneath it changed.

# Returns
`nothing`. The panel is redrawn on the dashboard.

# Examples
```julia
xyplot("Altitude",  t_a, alt_a; name = "Sat A")
xyplot!("Altitude", t_b, alt_b; name = "Sat B")
```
"""
function xyplot!(args...; kwargs...)
    title, rest = _split_panel(args)
    isempty(rest) && throw(ArgumentError(
        "plot! needs data. Give it y, or x and y — for example xyplot!(t, altitude). A leading " *
        "string names the panel: xyplot!(\"Altitude\", t, altitude)."))
    p = get_panel(title)
    traces = build_traces(:scatter, rest, kwargs; offset = length(p.traces))
    add_traces!(title, traces)
    ensure_browser()
    return nothing
end

# ─── The trace wrappers ───────────────────────────────────────────────────────────────────────

# Each kind takes its data differently. A contour takes two axes and a grid, a histogram takes
# one array, a scatter takes paired series — so each wrapper is routed to the builder that
# matches its shape.
const GRID_WRAPPERS = (:contour, :heatmap, :surface)

# A 3D scatter takes three coordinate arrays rather than two axes and a grid, so it needs its own
# builder: the paired-series rule below reads three arrays as an x/y pair with one left over.
const SPATIAL_WRAPPERS = (:scatter3d,)

# The first argument is a panel title or the first data array. Constraining it keeps these
# methods clear of PlotlyBase's own, which take no positional argument or a single Dict.
const WrapperArg = Union{AbstractString, AbstractArray, Tuple}

for f in TRACE_WRAPPERS
    bang  = Symbol(f, "!")
    build = f === :band            ? :(build_band(rest, kwargs)) :
            f in GRID_WRAPPERS     ? :(build_grid($(QuoteNode(f)), rest, kwargs)) :
            f in SPATIAL_WRAPPERS  ? :(build_spatial($(QuoteNode(f)), rest, kwargs)) :
                                     :(build_traces($(QuoteNode(f)), rest, kwargs))

    # Every `!` form differs from its partner in exactly one way, so it says that and points at
    # the partner for the rest (§11.1). Writing nine near-identical docstrings by hand would
    # leave nine places to drift.
    bang_doc = "    $bang([panel], data...; attributes...)\n\n" *
               "Add to a panel instead of replacing what is in it. Takes the same arguments " *
               "and attributes as [`$f`](@ref).\n\n" *
               "# Returns\n\n" *
               "`nothing`. The panel is redrawn on the dashboard.\n"

    @eval begin
        function $f(a::WrapperArg, more...; kwargs...)
            title, rest = _split_panel((a, more...))
            traces = $build
            set_traces!(title, traces)
            ensure_browser()
            return nothing
        end
        function $bang(a::WrapperArg, more...; kwargs...)
            title, rest = _split_panel((a, more...))
            get_panel(title)
            traces = $build
            add_traces!(title, traces)
            ensure_browser()
            return nothing
        end
        @doc $bang_doc $bang
    end
end

# ─── What the wrappers document ───────────────────────────────────────────────────────────────
# The functions above are generated, so their docstrings are attached here by name. Without
# these, `?contour` reaches PlotlyBase's own trace constructor and says nothing about panels,
# which is the whole difference.

"""
    scatterpolar([panel], θ, r; attributes...)

Draw against angle and radius — a ground station pass, an antenna pattern.

# Arguments
- `panel`: optional leading string naming the panel. Without one you get the default panel.
- `θ`: angle, in degrees unless you pass `thetaunit = "radians"`.
- `r`: radius, in whatever units the radial axis is set to.

# Notes
`thetaunit` converts the angle for you. There is no equivalent for the radial axis, so an
elevation in radians plotted against a degree range puts every point on the rim — convert it
yourself. A polar trace cannot share a panel with a cartesian one; Plotly draws them on
different subplot types.

# Returns
`nothing`. The panel is drawn on the dashboard.

# Examples
```julia
scatterpolar("Sky view", azimuth, rad2deg.(elevation); mode = "lines", name = "Pass 3")

panel!("Sky view";
       polar_angularaxis_direction = "clockwise",   # compass sense
       polar_angularaxis_rotation  = 90,            # north at the top
       polar_radialaxis_range      = [90, 0])       # zenith at the centre
```
"""
scatterpolar

"""
    bar([panel], categories, values; attributes...)

Draw one bar per category — ΔV by maneuver, passes by station.

# Arguments
- `panel`: optional leading string naming the panel.
- `categories`: the label for each bar. Strings or numbers.
- `values`: bar heights, one per category.

# Returns
`nothing`. The panel is drawn on the dashboard.

# Examples
```julia
bar("ΔV budget", ["TOI", "MCC", "MOI"], [2.82, 1.14, 0.48])
panel!("ΔV budget"; yaxis_title = "km/s")
```
"""
bar

"""
    histogram([panel], values; attributes...)

Draw how a set of values is distributed — residuals, Monte Carlo misses.

Plotly chooses the bins. Set `nbinsx` for an upper bound on their number, or `xbins_size` to fix
the width.

# Arguments
- `panel`: optional leading string naming the panel.
- `values`: the sample. One array, not a pair.

# Returns
`nothing`. The panel is drawn on the dashboard.

# Examples
```julia
histogram("Range residuals", residuals; nbinsx = 40)
```
"""
histogram

"""
    scattergeo([panel], longitude, latitude; attributes...)

Draw a track on a map of the body — a ground track, a station network.

# Arguments
- `panel`: optional leading string naming the panel.
- `longitude`: degrees east, in the body-fixed frame you computed them in.
- `latitude`: degrees north, same frame.

# Notes
Degrees, not radians, and geodetic rather than geocentric latitude if you want the track to sit
on the coastlines the way a map expects. The projection and the coastlines are Plotly's, set
through `panel!` with `geo_` attributes.

# Returns
`nothing`. The panel is drawn on the dashboard.

# Examples
```julia
scattergeo("Ground track", lon_deg, lat_deg; mode = "lines", name = "Explorer")
panel!("Ground track"; geo_projection_type = "equirectangular", geo_showland = true)
```
"""
scattergeo

"""
    scatter3d([panel], x, y, z; attributes...)

Draw a line or points in three axes of data — a B-plane, a parameter sweep.

This is a data plot with three axes, not a view of a body. For a trajectory around a central
body, with a clock and a globe, use [`orbitview`](@ref).

# Arguments
- `panel`: optional leading string naming the panel.
- `x`, `y`, `z`: one value each per point, all the same length.

# Returns
`nothing`. The panel is drawn on the dashboard.

# Examples
```julia
scatter3d("B-plane", b_dot_r, b_dot_t, epoch_offset; mode = "markers")
```
"""
scatter3d

"""
    band([panel], x, lower, upper; attributes...)

Draw a filled region between two bounds — a covariance envelope around residuals.

Draw the band before the data that sits inside it, so the points land on top rather than under.

# Arguments
- `panel`: optional leading string naming the panel.
- `x`: the horizontal axis, shared by both bounds.
- `lower`, `upper`: the two edges, each the same length as `x`.

# Notes
This is the one wrapper with no Plotly trace behind it. Plotly has no band, so this builds a
closed polygon — `x` forward along the upper bound and back along the lower — and fills it.

# Returns
`nothing`. The panel is drawn on the dashboard.

# Examples
```julia
band("Residuals", t, -3 .* σ, 3 .* σ; name = "3σ", fillcolor = "rgba(120,120,120,0.25)")
xyplot!("Residuals", t, residual; mode = "markers", name = "residual")
```
"""
band

for grid in GRID_WRAPPERS
    drawn = grid === :contour ? "as labelled contour lines" :
            grid === :heatmap ? "as coloured cells"         : "as a surface in three axes"

    @eval @doc """
        $($grid)([panel], x, y, Z; attributes...)
        $($grid)([panel], Z; attributes...)

    Draw a scalar field over a grid, $($drawn) — a porkchop plot, a coverage map.

    # Arguments
    - `panel`: optional leading string naming the panel.
    - `x`: the horizontal axis, one value per column of `Z`.
    - `y`: the vertical axis, one value per row of `Z`.
    - `Z`: the field, indexed `[row, column]`, so `size(Z) == (length(y), length(x))`.

    Given `Z` alone, the axes are its indices.

    # Notes
    Rows follow `y` and columns follow `x`. That is Plotly's convention, and it is why the
    horizontal axis is the second index rather than the first — the shape is checked, so a
    transposed grid is refused rather than drawn sideways.

    # Returns
    `nothing`. The panel is drawn on the dashboard.

    # Examples
    ```julia
    $($grid)("Porkchop", departure_days, arrival_days, C3;
             colorscale = "Viridis", colorbar_title = "C3 km²/s²")
    ```
    """ $grid
end

"""
Build the traces for one grid call: two axes and a matrix, or a matrix alone.

`Z` arrives indexed `[row, column]` and leaves as an array of rows, because handing Plotly a
Julia matrix leaves the orientation to the serializer.
"""
function build_grid(trace_type::Symbol, args, kwargs)
    if length(args) == 1
        Z = args[1]
        x, y = nothing, nothing
    elseif length(args) == 3
        x, y, Z = collect(args[1]), collect(args[2]), args[3]
    else
        throw(ArgumentError(
            "$trace_type takes x, y and Z — two axes and a grid — or Z alone. " *
            "Got $(length(args)) arrays."))
    end

    Z isa AbstractMatrix || throw(ArgumentError(
        "$trace_type needs a matrix for Z; got a $(typeof(Z))."))

    if x !== nothing
        size(Z) == (length(y), length(x)) || throw(ArgumentError(
            "Z is $(size(Z)) but x has $(length(x)) points and y has $(length(y)). " *
            "Z is indexed [row, column] with rows following y and columns following x, so " *
            "it should be $((length(y), length(x)))."))
    end

    # Plotly wants z as an array of rows. Handing it a Julia matrix leaves the orientation to
    # whatever the serializer decides, which is not a thing to leave to chance.
    rows = [plotdata(collect(view(Z, i, :))) for i in 1:size(Z, 1)]

    fields = Dict{Symbol, Any}(:z => rows)
    x === nothing || (fields[:x] = x)
    y === nothing || (fields[:y] = y)
    for (k, v) in pairs(kwargs)
        warn_unknown_attribute(trace_type, k)
        fields[k] = v
    end
    return [PlotlyBase.GenericTrace(trace_type; fields...)]
end

"""
    build_spatial(trace_type, args, kwargs) -> Vector{GenericTrace}

Build a trace from three coordinate arrays, for a kind that draws in three axes rather than over a
grid. Every array is a coordinate of the same points, so all three must be the same length.
"""
function build_spatial(trace_type::Symbol, args, kwargs)
    length(args) == 3 || throw(ArgumentError(
        "$trace_type takes x, y and z — three coordinate arrays of the same length. " *
        "Got $(length(args)) arrays."))

    x, y, z = collect(args[1]), collect(args[2]), collect(args[3])
    length(x) == length(y) == length(z) || throw(ArgumentError(
        "$trace_type needs x, y and z the same length; got $(length(x)), $(length(y)) and " *
        "$(length(z))."))

    fields = Dict{Symbol, Any}(:x => plotdata(x), :y => plotdata(y), :z => plotdata(z))
    for (k, v) in pairs(kwargs)
        warn_unknown_attribute(trace_type, k)
        fields[k] = v
    end
    return [PlotlyBase.GenericTrace(trace_type; fields...)]
end

"""
    band([panel], x, lo, hi; attributes...)

Draw a filled region between two bounds — a covariance band around a residual series.

Plotly draws this as one trace whose x runs forward then backward and whose y is the upper bound
followed by the reversed lower bound, filled to itself. Nobody should have to know that.
"""
function build_band(args, kwargs)
    length(args) == 3 || throw(ArgumentError(
        "band takes x, lo and hi — three arrays. Got $(length(args))."))
    x, lo, hi = collect(args[1]), collect(args[2]), collect(args[3])
    (length(x) == length(lo) == length(hi)) || throw(ArgumentError(
        "band got x of $(length(x)), lo of $(length(lo)) and hi of $(length(hi)). " *
        "All three must agree."))

    fields = Dict{Symbol, Any}(:fill => "toself", :mode => "lines", :hoverinfo => "skip")
    haskey(kwargs, :fillcolor) || (fields[:fillcolor] = "rgba(140,140,140,0.25)")
    haskey(kwargs, :line_width) || (fields[:line_width] = 0)
    for (k, v) in pairs(kwargs)
        warn_unknown_attribute(:scatter, k)
        fields[k] = v
    end
    return [PlotlyBase.GenericTrace(:scatter; x = plotdata(vcat(x, reverse(x))),
                                    y = plotdata(vcat(hi, reverse(lo))), fields...)]
end

# ─── panel! ───────────────────────────────────────────────────────────────────────────────────

"""
    panel!([panel]; layout...)
    panel!([panel], layout::PlotlyBase.Layout)

Set a panel's titles, axis labels, scales and limits.

# Arguments
- `panel`: the panel to configure. Without one you get the default panel.
- `layout...`: Plotly layout attributes. Nested settings are spelled with underscores, so
  `yaxis_type` reaches Plotly as `yaxis.type`.

# Notes
Settings accumulate rather than replace, so calling this twice for the same axis keeps both.
Limits are an ordered pair, so an axis runs backwards if you give them backwards. A `Dict` is
not a way in — write `legend_orientation = "h"`, not `legend = Dict("orientation" => "h")`.

Configuring a panel that does not exist yet is allowed and creates it empty, which is what lets
you set a panel up before plotting into it. A name that was meant to reach an existing panel and
does not earns a warning naming the panels there are.

# Returns
`nothing`. The panel is republished with the new layout.

# Examples
```julia
panel!("Altitude"; xaxis_title = "hours from epoch", yaxis_title = "km", yaxis_type = "log")

panel!("Decay"; xaxis_range = [6, 0])            # runs backwards
```
"""
function panel!(title::AbstractString = DEFAULT_PANEL; kwargs...)
    for (k, _) in pairs(kwargs)
        warn_unknown_layout(k)
    end
    set_layout!(title, layout_fields(; kwargs...))
    return nothing
end

# What an empty Layout already carries: a default margin and the whole light-theme template,
# some 15 kB of it.
const _LAYOUT_DEFAULTS = getfield(PlotlyBase.Layout(), :fields)

"""
Expand layout keywords the way Plotly needs them, keeping only what the caller actually set.

Going through `PlotlyBase.Layout` is what turns `yaxis_type` into `yaxis.type`. plotly.js does
not know the flattened spelling and ignores it in silence, so building the dictionary directly
made every `panel!` call do nothing at all.

Its defaults are then dropped. Carrying the default template would put 15 kB of light theme in
every panel on every poll, and it would argue with the dark styling the dashboard applies.
A value the caller set to something else is kept, defaults being compared by value rather than
by name.
"""
function layout_fields(; kwargs...)
    lay = try
        PlotlyBase.Layout(; kwargs...)
    catch e
        e isa MethodError || rethrow()
        # A Dict with String keys is the way into this: PlotlyBase expands nested layout
        # settings from underscored keywords, and cannot walk into one.
        throw(ArgumentError(
            "EpicycleIO: panel! could not build this layout. If a value is a Dict, write it " *
            "with underscores instead — `legend_orientation = \"h\"` rather than " *
            "`legend = Dict(\"orientation\" => \"h\")`. Underscores are how Plotly's nested " *
            "settings are spelled here, and they are what the documentation shows. " *
            "(Underlying error: $(sprint(showerror, e)))"))
    end
    return strip_layout_defaults(getfield(lay, :fields))
end

"""
Drop what an empty `Layout` already carries, comparing by value so a caller who
sets something to its default still gets it.

Both ways of setting a layout go through here. They did not always, and the
`Layout` form used to publish the whole 15 kB light-theme template on every poll
while the keyword form did not.
"""
strip_layout_defaults(f) =
    Dict{Symbol, Any}(k => v for (k, v) in f
                      if !(haskey(_LAYOUT_DEFAULTS, k) && _LAYOUT_DEFAULTS[k] == v))

function panel!(title::AbstractString, layout::PlotlyBase.Layout)
    set_layout!(title, strip_layout_defaults(getfield(layout, :fields)))
    return nothing
end

# ─── Handing over a trace directly ────────────────────────────────────────────────────────────
# Anything unwrapped is still reachable, which is what makes the interface future-proof without
# a list of trace types to maintain.

xyplot(title::AbstractString, tr::PlotlyBase.GenericTrace, more::PlotlyBase.GenericTrace...) =
    (set_traces!(title, [tr, more...]); ensure_browser(); nothing)

xyplot(tr::PlotlyBase.GenericTrace, more::PlotlyBase.GenericTrace...) =
    xyplot(DEFAULT_PANEL, tr, more...)

xyplot!(title::AbstractString, tr::PlotlyBase.GenericTrace, more::PlotlyBase.GenericTrace...) =
    (add_traces!(title, [tr, more...]); ensure_browser(); nothing)

xyplot!(tr::PlotlyBase.GenericTrace, more::PlotlyBase.GenericTrace...) =
    xyplot!(DEFAULT_PANEL, tr, more...)
