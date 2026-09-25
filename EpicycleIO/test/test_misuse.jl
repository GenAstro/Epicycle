# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# What happens when the interface is used wrongly.
#
# Asserting that something throws is the easy half and the less useful one. These tests check
# what the message SAYS, because a user meets these while already confused. The bar every one
# of them is held to:
#
#   1. name the thing that is wrong — which attribute, which panel, which argument
#   2. give the actual numbers, not "lengths must match"
#   3. say what to do instead, not only what went wrong
#
# The third is the one that gets dropped, so it is checked explicitly.

using Test
using EpicycleIO
using PlotlyBase
using JSON

const IO_ = EpicycleIO
IO_.auto_open!(false)

"""Run `f` and return the error message it produced, or `nothing` if it did not throw."""
function message(f)
    try
        f()
        return nothing
    catch e
        return sprint(showerror, e)
    end
end

"""Every fragment must appear in the message. Reports which one was missing rather than false."""
function says(msg, fragments...)
    msg === nothing && return "did not throw at all"
    missing_ = [f for f in fragments if !occursin(f, msg)]
    isempty(missing_) && return true
    return "message is missing $(missing_): " * msg
end

@testset "misuse is trapped, and says something useful" begin

@testset "plot, no data" begin
    m = message(() -> xyplot())
    @test says(m, "plot needs data") === true
    @test says(m, "xyplot(t, altitude)") === true            # shows the right shape
end

@testset "mismatched x and y" begin
    m = message(() -> xyplot(collect(1.0:5.0), collect(1.0:3.0)))
    @test says(m, "5", "3") === true                       # both actual counts
    @test says(m, "agree") === true
end

@testset "an array without a partner" begin
    m = message(() -> xyplot(collect(1.0:3.0), collect(1.0:3.0), collect(1.0:3.0)))
    @test says(m, "3 positional arguments") === true
    @test says(m, "x/y pairs") === true                    # what to do instead
end

@testset "a matrix the wrong way round" begin
    x = collect(1.0:4.0)
    Y = rand(3, 2)                                          # 3 rows against 4 points
    m = message(() -> xyplot(x, Y))
    @test says(m, "4", "3 rows") === true
    @test says(m, "Columns are series") === true           # states the convention
end

@testset "ragged components" begin
    x = collect(1.0:3.0)
    y = [[1.0, 2.0, 3.0], [1.0, 2.0], [1.0, 2.0, 3.0]]     # sample 2 is short
    m = message(() -> xyplot(x, y))
    @test says(m, "ragged", "sample 2") === true           # names WHICH sample
end

@testset "a vector of vectors that is neither shape" begin
    x = collect(1.0:5.0)
    y = [rand(4), rand(4)]                                  # 4-long, against 5 points
    m = message(() -> xyplot(x, y))
    @test says(m, "5 points", "2 vectors") === true
    @test says(m, "one entry per sample", "one full-length vector per series") === true
end

@testset "an attribute with the wrong number of values" begin
    t = collect(0.0:0.1:1.0)
    r = [[i, 2i, 3i] for i in t]                            # draws three traces
    m = message(() -> xyplot(t, r; line_color = ["red", "green"]))
    @test says(m, "line_color", "2 entries", "3 traces") === true
    # And why it cannot mean anything else, which is the non-obvious part.
    @test says(m, "one value per point") === true
end

@testset "nesting with the wrong number of groups" begin
    t = collect(0.0:0.1:1.0)
    r = [[i, 2i, 3i] for i in t]
    m = message(() -> xyplot(t, r; marker_size = [[4], [8]]))
    @test says(m, "marker_size", "2 entries", "3 traces") === true
end

@testset "band with the wrong arguments" begin
    x = collect(1.0:4.0)
    @test says(message(() -> band(x, fill(-1.0, 4))),
               "band takes x, lo and hi", "Got 2") === true
    @test says(message(() -> band(x, fill(-1.0, 4), fill(1.0, 3))),
               "4", "3", "must agree") === true
end

@testset "a grid that is transposed" begin
    x = collect(1.0:4.0)          # columns
    y = collect(1.0:3.0)          # rows
    Z = rand(3, 4)
    m = message(() -> contour(y, x, Z))                     # axes swapped
    @test says(m, "[row, column]") === true                 # states the convention
    @test says(m, "should be") === true                     # and the shape wanted
end

@testset "a grid that is not a grid" begin
    x = collect(1.0:3.0)
    @test says(message(() -> contour(x, x, x)), "needs a matrix for Z") === true
    @test says(message(() -> contour(x, x)),
               "two axes and a grid", "Got 2") === true
end

@testset "mixing coordinate families" begin
    IO_.reset_panels!()
    t = collect(0.0:0.5:2.0)
    xyplot("Mixed", t, t)
    m = message(() -> scatterpolar!("Mixed", t, t))
    @test says(m, "Mixed", "cartesian", "polar") === true   # names the panel and both kinds
    @test says(m, "panel of its own") === true              # what to do about it
    # Plotly's constraint, not ours — the message should say so rather than sound arbitrary.
    @test says(m, "different subplot types") === true
end

@testset "a layout Dict that cannot be expanded" begin
    IO_.reset_panels!()
    t = collect(0.0:0.5:2.0)
    xyplot("L", t, t)
    m = message(() -> panel!("L"; legend = Dict("orientation" => "h")))
    @test says(m, "legend_orientation") === true            # the exact thing to write
    @test says(m, "underscores") === true
end

@testset "report misuse" begin
    dir  = mktempdir()
    path = joinpath(dir, "f.txt")

    @test says(message(() -> report(path)),
               "at least one column", "time = t") === true

    @test says(message(() -> report(path; a = [1.0, 2.0], b = [1.0])),
               "a 2", "b 1") === true                       # names both and their lengths

    # An epoch is the value that turns up in practice, so the message names the way out.
    struct Multiline end
    Base.show(io::IO, ::Multiline) = print(io, "one\ntwo")
    m = message(() -> report(path; bad = [Multiline(), Multiline()]))
    @test says(m, "spans", "t.jd", "t.isot") === true
end

# ─── Warnings, for the misuse that cannot be an error ─────────────────────────────────────────

@testset "an attribute Plotly does not have" begin
    IO_.reset_panels!()
    t = collect(0.0:0.5:2.0)
    # PlotlyBase drops what it does not recognise, so without this the panel is simply wrong
    # with nothing said anywhere.
    msg = @test_logs (:warn,) match_mode = :any xyplot(t, t; line_colour = "cyan")
    @test true
end

@testset "the unknown-attribute warning names the fix" begin
    IO_.reset_panels!()
    t = collect(0.0:0.5:2.0)
    out = IOBuffer()
    logger = Base.CoreLogging.SimpleLogger(out, Base.CoreLogging.Warn)
    Base.CoreLogging.with_logger(logger) do
        xyplot(t, t; line_colour = "cyan")          # British spelling: not a Plotly attribute
    end
    text = String(take!(out))
    @test occursin("line_colour", text)           # names the offender
    @test occursin("dropped", text)               # says what will happen
    @test occursin("plotly.com", text)            # where to look it up

    # The schema ships with PlotlyBase and describes plotly.js 2.0 or 2.1, while the dashboard
    # renders with 2.35.2. An attribute newer than the schema is real and warns anyway, so the
    # message must not claim the attribute will be dropped as though that were certain.
    @test occursin("schema shipped with PlotlyBase", text)
    @test occursin("behind the renderer", text)
end

@testset "the schema is older than the renderer, and says so" begin
    # Attributes plotly.js 2.35.2 has and the shipped schema does not. If a PlotlyBase bump ever
    # makes these true, the caveat in the warning can be softened and this test says so.
    @test !IO_.attr_exists(:scatter, :zorder)
    @test !IO_.attr_exists(:scatter, :legend_grouptitlefont)

    # And what it does still catch, which is the reason to keep it.
    @test !IO_.attr_exists(:scatter, :line_colour)
    @test IO_.attr_exists(:scatter, :line_color)
end

@testset "a layout attribute Plotly does not have" begin
    IO_.reset_panels!()
    t = collect(0.0:0.5:2.0)
    xyplot("A", t, t)
    out = IOBuffer()
    logger = Base.CoreLogging.SimpleLogger(out, Base.CoreLogging.Warn)
    Base.CoreLogging.with_logger(logger) do
        panel!("A"; xaxis_headline = "hours")
    end
    text = String(take!(out))
    @test occursin("xaxis_headline", text)
    @test occursin("layout", text)
end

@testset "clear! before the panel exists is a silent no-op" begin
    IO_.reset_panels!()
    t = collect(0.0:0.5:2.0)
    xyplot("Real", t, t)

    # Emptying a script's panels before filling them, so a re-run does not accumulate, runs
    # before those panels exist the first time. It is ordinary, not a mistake.
    out = IOBuffer()
    Base.CoreLogging.with_logger(
        Base.CoreLogging.SimpleLogger(out, Base.CoreLogging.Warn)) do
        clear!("Not yet plotted")
    end
    @test isempty(String(take!(out)))                            # says nothing
    @test !haskey(IO_._PANELS, IO_.panel_id("Not yet plotted"))  # and creates nothing
    @test length(IO_.panels()) == 1
end

@testset "a mistyped panel name" begin
    IO_.reset_panels!()
    t = collect(0.0:0.5:2.0)
    xyplot("Altitude", t, t)

    out = IOBuffer()
    logger = Base.CoreLogging.SimpleLogger(out, Base.CoreLogging.Warn)
    Base.CoreLogging.with_logger(logger) do
        panel!("Altitde"; yaxis_title = "km")     # typo: would silently make a second panel
    end
    text = String(take!(out))
    @test occursin("Altitde", text)               # what was given
    @test occursin("Altitude", text)              # what probably was meant
    @test occursin("stays empty", text)           # what happens otherwise

    # Configuring a panel before plotting into it is legitimate, so the first one is silent.
    IO_.reset_panels!()
    out2 = IOBuffer()
    logger2 = Base.CoreLogging.SimpleLogger(out2, Base.CoreLogging.Warn)
    Base.CoreLogging.with_logger(logger2) do
        panel!("Fresh"; yaxis_title = "km")
    end
    @test !occursin("Fresh", String(take!(out2)))
end

@testset "gaps in the data are drawn as gaps, not refused" begin
    IO_.reset_panels!()
    t = [1.0, 2.0, 3.0, 4.0]

    # JSON has no NaN, and JSON.json refuses it outright, so one gap used to take the whole
    # call down. Gaps are ordinary here: between passes, a circular orbit's argument of
    # periapsis, a solver step that did not converge.
    @test xyplot("Gappy", t, [1.0, NaN, Inf, 4.0]) === nothing

    j = JSON.json(IO_.figure(IO_.get_panel("Gappy")))
    parsed = JSON.parse(j)                     # would throw if we emitted a bare NaN
    @test parsed["data"][1]["y"] == [1.0, nothing, nothing, 4.0]

    # A clean series is not copied.
    clean = [1.0, 2.0, 3.0]
    @test IO_.plotdata(clean) === clean

    # And through the other builders.
    @test band("B", t, [0.0, NaN, 0.0, 0.0], [1.0, 1.0, 1.0, 1.0]) === nothing
    @test contour("C", [1.0, 2.0], [1.0, 2.0], [1.0 NaN; 2.0 3.0]) === nothing
    @test JSON.parse(JSON.json(IO_.figure(IO_.get_panel("C"))))["data"][1]["z"][1] == [1.0, nothing]
end

@testset "clearing reaches the browser, not just Julia" begin
    t = collect(0.0:0.5:2.0)
    xyplot("Gone1", t, t)
    xyplot("Gone2", t, 2t)
    IO_.flush_publishes(); sleep(0.3)
    dir = IO_.assets_dir()
    panel_files() = filter(f -> startswith(f, "panel_") && endswith(f, ".json"), readdir(dir))
    titles() = [p["title"] for p in
                JSON.parse(read(joinpath(dir, "dashboard.json"), String))["panels"]]

    @test "Gone1" in titles()
    @test !isempty(panel_files())

    # Clearing only the Julia side left the page showing everything, which looked like the
    # call had done nothing at all.
    clear_all!()
    IO_.flush_publishes(); sleep(0.3)
    @test isempty(titles())
    @test isempty(panel_files())          # no orphans for the dashboard to be asked for

    # And the board is usable again afterwards.
    xyplot("After", t, t)
    IO_.flush_publishes(); sleep(0.3)
    @test titles() == ["After"]
end

# ─── Misuse that should NOT be trapped ────────────────────────────────────────────────────────
# Refusing something legitimate is its own failure, so the boundary is worth pinning.

@testset "these are fine and must stay fine" begin
    IO_.reset_panels!()
    t = collect(0.0:0.5:2.0)

    @test xyplot(t, t) === nothing                                  # the ordinary case
    @test xyplot(t) === nothing                                     # index as the axis
    @test xyplot(t, t, collect(1.0:3.0), collect(1.0:3.0)) === nothing   # unequal pairs
    @test xyplot(Float64[], Float64[]) === nothing                  # empty is not an error
    @test xyplot(t, t; marker_size = collect(1.0:5.0)) === nothing   # arrayOk, per point
    @test panel!("Plot"; xaxis_title = "hours") === nothing
end

end # testset
