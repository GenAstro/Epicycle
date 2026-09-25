# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

using Test
using EpicycleIO
using PlotlyBase
using JSON

const IO_ = EpicycleIO

# Do not pop a browser tab while testing. The server still starts, on localhost, which is what
# the publish path needs.
IO_.auto_open!(false)

fields(tr) = getfield(tr, :fields)

@testset "EpicycleIO" begin

# ── Shape expansion ──────────────────────────────────────────────────────────────
@testset "shape expansion" begin
    x = collect(1.0:5.0)
    y = collect(6.0:10.0)

    @test length(IO_.expand_series(y)) == 1
    @test IO_.expand_series(y)[1].x == collect(1:5)          # index becomes the axis

    s = IO_.expand_series(x, y)
    @test length(s) == 1 && s[1].x == x && s[1].y == y

    # A matrix: columns are series.
    Y = hcat(y, y .* 2, y .* 3)
    s = IO_.expand_series(x, Y)
    @test length(s) == 3
    @test s[2].y == y .* 2
    @test s[3].component == 3 && s[3].ncomponents == 3

    # Sample-major vector of vectors, which is what `history` returns for a position column.
    r = [[i, 2i, 3i] for i in 1.0:5.0]
    s = IO_.expand_series(x, r)
    @test length(s) == 3
    @test s[1].y == collect(1.0:5.0)
    @test s[3].y == 3 .* collect(1.0:5.0)

    # Series-major: one full-length vector per series.
    s = IO_.expand_series(x, [y, y .* 2])
    @test length(s) == 2 && s[2].y == y .* 2

    # Independent pairs, different lengths.
    s = IO_.expand_series(x, y, collect(1.0:3.0), collect(4.0:6.0))
    @test length(s) == 2
    @test length(s[1].y) == 5 && length(s[2].y) == 3

    @test_throws ArgumentError IO_.expand_series(x, collect(1.0:4.0))   # mismatched pair
    @test_throws ArgumentError IO_.expand_series(x, y, x)              # odd count
    @test_throws ArgumentError IO_.expand_series()
end

# ── arrayOk drives what a vector means ───────────────────────────────────
@testset "arrayOk from the shipped schema" begin
    @test IO_.is_array_ok(:scatter, :marker_size)          # per data point
    @test IO_.is_array_ok(:scatter, :marker_color)
    @test !IO_.is_array_ok(:scatter, :name)                # so a vector can only be per-trace
    @test !IO_.is_array_ok(:scatter, :line_width)
    @test !IO_.is_array_ok(:scatter, :line_color)

    @test IO_.attr_exists(:scatter, :line_width)
    @test IO_.attr_exists(:scatter, :mode)
    @test !IO_.attr_exists(:scatter, :nonsense_attribute)
    @test IO_.layout_attr_exists(:xaxis_title)
    @test !IO_.layout_attr_exists(:not_a_layout_thing)
end

@testset "attribute distribution" begin
    # Not arrayOk: one per trace.
    @test IO_.distribute(:line_color, ["red", "green", "blue"], 2, 3, :scatter) == "green"
    @test IO_.distribute(:line_width, [1, 2, 3], 3, 3, :scatter) == 3
    # Scalar: every trace.
    @test IO_.distribute(:line_width, 1.5, 2, 3, :scatter) == 1.5
    # arrayOk: Plotly's per-point meaning is left alone.
    @test IO_.distribute(:marker_size, [4, 5, 6], 2, 3, :scatter) == [4, 5, 6]
    # Nesting always means per-trace.
    @test IO_.distribute(:marker_size, [[4], [8], [12]], 3, 3, :scatter) == [12]
    # A wrong count is an error, not a recycle.
    @test_throws ArgumentError IO_.distribute(:line_color, ["red", "green"], 1, 3, :scatter)
end

@testset "name is the one attribute that rewrites" begin
    s3 = IO_.Series([1.0], [1.0], 2, 3)
    s1 = IO_.Series([1.0], [1.0], 0, 1)
    @test IO_.series_name("position", s3) == "position[2]"   # never three identical entries
    @test IO_.series_name("altitude", s1) == "altitude"      # nothing to disambiguate
    @test IO_.series_name(nothing, s3) === nothing
end

# ── Traces ───────────────────────────────────────────────────────────────────────
@testset "trace building" begin
    IO_.reset_panels!()
    t = collect(0.0:0.1:1.0)
    r = [[i, 2i, 3i] for i in t]

    traces = IO_.build_traces(:scatter, (t, r), (; name = "position", line_width = 2))
    @test length(traces) == 3
    @test fields(traces[1])[:name] == "position[1]"
    @test fields(traces[2])[:line][:width] == 2                # underscore flattening
    @test haskey(fields(traces[1])[:line], :color)             # palette applied

    traces = IO_.build_traces(:scatter, (t, r),
                              (; name = ["x", "y", "z"], line_color = ["red", "green", "blue"]))
    @test [fields(tr)[:name] for tr in traces] == ["x", "y", "z"]
    @test fields(traces[3])[:line][:color] == "blue"

    # Polar carries theta/r, geographic lon/lat.
    pol = IO_.build_traces(:scatterpolar, (t, t), (; thetaunit = "radians"))
    @test haskey(fields(pol[1]), :theta) && haskey(fields(pol[1]), :r)
    geo = IO_.build_traces(:scattergeo, (t, t), NamedTuple())
    @test haskey(fields(geo[1]), :lon) && haskey(fields(geo[1]), :lat)

    # A long series is promoted to the WebGL path, which is what keeps a trajectory
    # responsive. PlotlyBase stores the type as a Symbol.
    long  = collect(1.0:5000.0)
    short = collect(1.0:10.0)
    @test fields(IO_.build_traces(:scatter, (long, long), NamedTuple())[1])[:type] === :scattergl
    @test fields(IO_.build_traces(:scatter, (short, short), NamedTuple())[1])[:type] === :scatter
    # No scatterpolargl exists, so a polar trace is never promoted.
    @test fields(IO_.build_traces(:scatterpolar, (long, long), NamedTuple())[1])[:type] ===
          :scatterpolar
end

@testset "band is a shape helper, not a kind" begin
    x  = collect(1.0:4.0)
    lo = fill(-1.0, 4)
    hi = fill(1.0, 4)
    tr = IO_.build_band((x, lo, hi), (; name = "3sigma"))[1]
    f  = fields(tr)
    @test f[:fill] == "toself"
    @test length(f[:x]) == 8                      # forward then back
    @test f[:y] == vcat(hi, reverse(lo))
    @test_throws ArgumentError IO_.build_band((x, lo), NamedTuple())
end

@testset "grids: contour, heatmap, surface" begin
    IO_.reset_panels!()
    x = collect(1.0:4.0)          # columns
    y = collect(1.0:3.0)          # rows
    Z = [10i + j for i in 1:3, j in 1:4]

    tr = IO_.build_grid(:contour, (x, y, Z), NamedTuple())[1]
    f  = fields(tr)
    @test length(f[:z]) == 3 && length(f[:z][1]) == 4   # array of rows, rows follow y
    @test f[:z][2] == [21.0, 22.0, 23.0, 24.0]
    @test f[:x] == x && f[:y] == y

    @test length(fields(IO_.build_grid(:surface, (Z,), NamedTuple())[1])[:z]) == 3

    # A transposed grid is caught rather than drawn sideways.
    @test_throws ArgumentError IO_.build_grid(:contour, (y, x, Z), NamedTuple())
    @test_throws ArgumentError IO_.build_grid(:contour, (x, y), NamedTuple())
    @test_throws ArgumentError IO_.build_grid(:contour, (x, y, x), NamedTuple())
end

@testset "wrappers extend PlotlyBase rather than shadow it" begin
    # Our wrappers are named after the traces they draw, which is deliberate — and it means
    # PlotlyBase exports every one of them. Rival functions would make each name ambiguous the
    # moment a user loaded both, which the escape hatch asks them to do. Adding methods to the
    # same binding is what keeps `using EpicycleIO, PlotlyBase` working.
    for f in IO_.TRACE_WRAPPERS
        f === :band && continue                      # ours alone, PlotlyBase has no band
        @test getfield(IO_, f) === getfield(PlotlyBase, f)
    end

    # PlotlyBase's own way of building a trace still works and is untouched.
    @test contour(z = [1 2; 3 4]) isa PlotlyBase.GenericTrace
    @test bar(x = [1, 2], y = [3, 4]) isa PlotlyBase.GenericTrace

    # And ours draws.
    IO_.reset_panels!()
    contour("Grid", collect(1.0:3.0), collect(1.0:2.0), [1.0 2.0 3.0; 4.0 5.0 6.0])
    @test length(IO_.get_panel("Grid").traces) == 1

    # Attributes reach a grid trace too, not just a scatter.
    contour("Grid", collect(1.0:3.0), collect(1.0:2.0), [1.0 2.0 3.0; 4.0 5.0 6.0];
            colorscale = "Viridis")
    @test fields(IO_.get_panel("Grid").traces[1])[:colorscale] == "Viridis"

    # Every wrapper has a `!` form that adds rather than replaces. The suite used
    # to reach these only through the error path, so the succeeding one was never
    # actually run.
    IO_.reset_panels!()
    az = collect(range(0, 3.0; length = 5))
    scatterpolar("Sky", az, az)
    scatterpolar!("Sky", az, 2 .* az; name = "second pass")
    @test length(IO_.get_panel("Sky").traces) == 2

    IO_.reset_panels!()
    bar("Budget", ["TOI", "MCC"], [120.0, 4.0])
    bar!("Budget", ["TOI", "MCC"], [10.0, 1.0]; name = "margin")
    @test length(IO_.get_panel("Budget").traces) == 2
end

@testset "the escape hatch the docs promise" begin
    # `data_plotting.md` tells the reader that a wrapper is a shortcut, never a
    # gate: build any of Plotly's forty-odd traces with PlotlyBase and hand it
    # over. That promise is only worth making if it is checked.
    IO_.reset_panels!()

    xyplot("Residual distribution", PlotlyBase.violin(y = randn(50), name = "range"))
    p = IO_.get_panel("Residual distribution")
    @test length(p.traces) == 1
    # PlotlyBase's own constructors record the type as a String, where the traces
    # we build carry a Symbol. Both serialise the same, so compare as text.
    @test string(fields(p.traces[1])[:type]) == "violin"
    @test fields(p.traces[1])[:name] == "range"

    # Several at once, and a PlotlyBase.Layout rather than keywords.
    IO_.reset_panels!()
    xyplot("Two", PlotlyBase.violin(y = randn(10), name = "a"),
                PlotlyBase.violin(y = randn(10), name = "b"))
    @test length(IO_.get_panel("Two").traces) == 2

    panel!("Two", PlotlyBase.Layout(title = "By station", yaxis_title = "km"))
    lay = IO_.figure(IO_.get_panel("Two"))[:layout]
    @test lay[:title] == "By station"
    @test lay[:yaxis][:title][:text] == "km"    # underscores still nest

    # This form used to publish PlotlyBase's whole light-theme template — 15 kB
    # on every poll, arguing with the dashboard's dark styling — because only the
    # keyword path stripped the defaults.
    @test !haskey(lay, :template)
    @test length(JSON.json(lay)) < 500
end

# ── Panels ───────────────────────────────────────────────────────────────────────
@testset "panels" begin
    IO_.reset_panels!()
    t = collect(0.0:0.5:2.0)

    xyplot("Altitude", t, t; name = "A")
    @test length(IO_.get_panel("Altitude").traces) == 1

    xyplot!("Altitude", t, 2 .* t; name = "B")
    @test length(IO_.get_panel("Altitude").traces) == 2

    # Same name replaces, so re-running a script leaves one figure rather than twenty.
    xyplot("Altitude", t, t; name = "A again")
    @test length(IO_.get_panel("Altitude").traces) == 1

    clear!("Altitude")
    @test isempty(IO_.get_panel("Altitude").traces)

    @test IO_.panel_id("Ground track, polar") == "Ground_track__polar"

    # Plotly cannot draw a polar trace on a cartesian panel, and says so.
    xyplot("Mixed", t, t)
    @test_throws ArgumentError scatterpolar!("Mixed", t, t)
end

@testset "layout and figure" begin
    IO_.reset_panels!()
    t = collect(0.0:0.5:2.0)
    xyplot("Alt", t, t)
    panel!("Alt"; xaxis_title = "hours", yaxis_type = "log")

    fig = IO_.figure(IO_.get_panel("Alt"))
    @test fig[:layout][:xaxis][:title][:text] == "hours"   # nested, as Plotly needs
    @test length(fig[:data]) == 1
    @test !isempty(JSON.json(fig))

    m = IO_.manifest()
    @test any(p -> p[:title] == "Alt" && p[:kind] == "xy", m[:panels])
end

@testset "layout reaches the browser in the shape Plotly wants" begin
    IO_.reset_panels!()
    t = collect(0.0:0.5:2.0)
    xyplot("Alt", t, t)

    # yaxis_type has to arrive as yaxis.type. plotly.js does not know the flattened spelling
    # and ignores it in silence, which made every panel! call a no-op that looked fine.
    panel!("Alt"; xaxis_title = "hours", yaxis_type = "log")
    lay = IO_.figure(IO_.get_panel("Alt"))[:layout]
    @test lay[:xaxis][:title][:text] == "hours"
    @test lay[:yaxis][:type] == "log"
    @test !haskey(lay, :xaxis_title)

    # A second call must not erase the first: both live under :xaxis.
    panel!("Alt"; xaxis_range = [6, 0])
    lay = IO_.figure(IO_.get_panel("Alt"))[:layout]
    @test lay[:xaxis][:range] == [6, 0]
    @test lay[:xaxis][:title][:text] == "hours"

    # PlotlyBase's default light template is 15 kB and would be sent on every poll.
    @test !haskey(lay, :template)
    @test length(JSON.json(lay)) < 500

    # An explicitly set default is still kept: dropped by value, not by name.
    @test haskey(IO_.layout_fields(margin_l = 10), :margin)

    # Polar nests three deep.
    IO_.reset_panels!()
    az = collect(range(0, 3.0; length = 10))
    scatterpolar("Sky", az, az; thetaunit = "radians")
    panel!("Sky"; polar_angularaxis_direction = "clockwise", polar_radialaxis_range = [90, 0])
    lay = IO_.figure(IO_.get_panel("Sky"))[:layout]
    @test lay[:polar][:angularaxis][:direction] == "clockwise"
    @test lay[:polar][:radialaxis][:range] == [90, 0]

    # A Dict with String keys cannot be walked into. Say so, and say what to write instead.
    IO_.reset_panels!()
    xyplot("Legend", t, t)
    err = try
        panel!("Legend"; legend = Dict("orientation" => "h")); nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("legend_orientation", sprint(showerror, err))
    @test IO_.layout_fields(legend_orientation = "h")[:legend][:orientation] == "h"

    # The radial axis has no unit attribute, so a pass given in radians against a degree
    # range lands on the rim rather than failing. Guard the example against regressing.
    IO_.reset_panels!()
    el_deg = [2.0, 45.0, 80.0]
    scatterpolar("Pass", [0.1, 0.2, 0.3], el_deg; thetaunit = "radians")
    r = fields(IO_.get_panel("Pass").traces[1])[:r]
    @test maximum(r) > 45           # elevation in degrees, not radians

    # Clearing takes the data out and leaves the panel looking how it was told to look.
    IO_.reset_panels!()
    xyplot("Keep", t, t)
    panel!("Keep"; yaxis_type = "log")
    clear!("Keep")
    @test IO_.figure(IO_.get_panel("Keep"))[:layout][:yaxis][:type] == "log"
end

# ── Publishing ──────────────────────────────────────────────────────────────────
@testset "atomic write" begin
    dir  = mktempdir()
    path = joinpath(dir, "x.json")
    IO_.atomic_write(path, "{\"a\":1}")
    @test read(path, String) == "{\"a\":1}"
    IO_.atomic_write(path, "{\"a\":2}")
    @test read(path, String) == "{\"a\":2}"
    @test isempty(filter(startswith(".tmp_"), readdir(dir)))    # no debris left behind
end

@testset "throttling coalesces and never loses the last value" begin
    IO_.reset_panels!()
    t = collect(0.0:1.0:3.0)
    for k in 1:40                    # a solver publishing every iteration
        xyplot!("Convergence", [Float64(k)], [1.0 / k]; mode = "markers")
    end
    IO_.flush_publishes()
    sleep(0.4)
    written = joinpath(IO_.assets_dir(), "panel_Convergence.json")
    @test isfile(written)
    fig = JSON.parse(read(written, String))
    @test length(fig["data"]) == 40           # every point is there, however few writes ran
end

@testset "a closed window is reopened" begin
    IO_.reset_panels!()
    IO_.auto_open!(false)                    # count decisions, do not launch anything
    t = collect(0.0:0.5:2.0)

    # Cold: nothing has ever polled, so a viewer is needed.
    IO_._LAST_REQUEST[] = 0.0
    IO_._LAST_OPEN[]    = 0.0
    @test IO_.wants_browser()

    # A tab is polling: leave it alone rather than spawning another on every call.
    IO_._LAST_REQUEST[] = time()
    @test !IO_.wants_browser()

    # The window was closed, so polling stopped. This is the case the old latch got wrong:
    # it remembered opening one and never offered another.
    IO_._LAST_REQUEST[] = time() - 2 * IO_.VIEWER_TIMEOUT[]
    IO_._LAST_OPEN[]    = time() - 2 * IO_.OPEN_GRACE[]
    @test IO_.wants_browser()

    # Just opened one; it has not started polling yet. Do not open a second.
    IO_._LAST_OPEN[] = time()
    @test !IO_.wants_browser()

    # A cold browser can take longer to appear than a script takes to run. Treating that
    # silence as "no viewer" opened a second tab onto the same plots while the first was
    # still starting.
    IO_._LAST_REQUEST[] = 0.0
    IO_._LAST_OPEN[]    = time() - 2 * IO_.OPEN_GRACE[]
    @test !IO_.wants_browser()

    # But not forever: if nothing ever connects, offer another eventually.
    IO_._LAST_OPEN[] = time() - 2 * IO_.COLD_GRACE[]
    @test IO_.wants_browser()
end

# ── report ─────────────────────────────────────────────────────────────────────
@testset "report" begin
    dir  = mktempdir()
    path = joinpath(dir, "flight.txt")
    t = collect(1.0:4.0)
    r = [[i, 2i, 3i] for i in t]

    report(path; time = t, position = r)
    lines = readlines(path)
    @test occursin("position_1", lines[1]) && occursin("position_3", lines[1])
    @test length(lines) == 5                      # header plus four rows

    # Every row carries as many fields as the header. A column whose printed form
    # broke would misalign the file rather than fail, and nothing else would notice.
    ncol = length(split(lines[1]))
    @test all(length(split(l)) == ncol for l in lines[2:end])

    @test_throws ArgumentError report(path; a = [1.0, 2.0], b = [1.0])
    @test_throws ArgumentError report(path)

    # A value whose printed form spans lines would silently turn one row into several.
    struct _Tall end
    Base.show(io::IO, ::_Tall) = print(io, "line one\nline two")
    @test_throws ArgumentError report(path; bad = [_Tall(), _Tall()])
end

end # testset

EpicycleIO.close_dashboard()

include("test_misuse.jl")
include("test_orbitview.jl")
include("test_narration.jl")
include("test_server.jl")     # last: it closes the server on its way out

nothing
