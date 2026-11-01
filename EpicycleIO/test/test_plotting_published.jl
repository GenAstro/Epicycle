# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# What each plotting function publishes: called the way a user calls it, then read back as the
# JSON Plotly receives. The data keys and values for the trace type, the color attribute it
# draws with, nothing that changes how Plotly reads the data, and the panel's kind and file.
# Documented attributes are checked to pass without a warning, as well as misspelled ones to warn.

using Test
using EpicycleIO
using JSON

const PUB = EpicycleIO

published(id) = JSON.parse(PUB.payload(PUB._PANELS[PUB.panel_id(id)]))
traces(id) = published(id)["data"]
entry(id) = only(e for e in PUB.manifest()[:panels] if e[:id] == PUB.panel_id(id))
quiet(f) = @test_logs min_level = Base.CoreLogging.Warn f()          # passes only with no warning

@testset "plotting: what each function publishes" begin
    PUB.auto_open!(false)

    @testset "xy, polar and geographic data reach the keys Plotly reads" begin
        quiet(() -> xyplot("Values xy", [1.0, 2, 3], [4.0, 5, 6]))
        t = only(traces("Values xy"))
        @test t["type"] == "scatter" && t["x"] == [1, 2, 3] && t["y"] == [4, 5, 6]
        quiet(() -> scatterpolar("Values polar", [0.0, 90, 180], [1.0, 2, 3]))
        t = only(traces("Values polar"))
        @test t["theta"] == [0, 90, 180] && t["r"] == [1, 2, 3] && !haskey(t, "x")
        quiet(() -> scattergeo("Values geo", [-77.0, 2.3], [38.9, 48.9]))
        t = only(traces("Values geo"))
        @test t["lon"] == [-77, 2.3] && t["lat"] == [38.9, 48.9] && !haskey(t, "x")
        quiet(() -> bar("Values bar", ["TOI", "MCC"], [2.82, 1.14]))
        t = only(traces("Values bar"))
        @test t["x"] == ["TOI", "MCC"] && t["y"] == [2.82, 1.14]
    end

    @testset "histogram bins the values" begin
        v = [0.1, -0.2, 0.3, 0.25]
        quiet(() -> histogram("Hist", v; nbinsx = 4))
        t = only(traces("Hist"))
        @test t["type"] == "histogram" && t["x"] == v && !haskey(t, "y")
        @test_throws "takes one array of values" histogram("Hist", v, v)
    end

    @testset "contour and heatmap keep the grid's orientation" begin
        x, y = [1.0, 2.0, 3.0], [10.0, 20.0]
        Z = [1.0 2.0 3.0; 4.0 5.0 6.0]                                     # rows follow y
        for (f, id) in ((contour, "Grid contour"), (heatmap, "Grid heatmap"))
            quiet(() -> f(id, x, y, Z))
            t = only(traces(id))
            @test t["transpose"] == false && t["x"] == x && t["y"] == y
            @test length(t["z"]) == length(y) && t["z"][2][3] == Z[2, 3]
        end
        # A gap in an axis becomes null, as in a data series, rather than failing to publish.
        quiet(() -> contour("Grid gap", [1.0, NaN, 3.0], y, Z))
        @test traces("Grid gap")[1]["x"] == [1.0, nothing, 3.0]
    end

    @testset "3D traces publish a Plotly figure, not a Cesium scene" begin
        quiet(() -> surface("3D surface", [1.0 2.0; 3.0 4.0]))
        quiet(() -> scatter3d("3D points", [1.0, 2], [3.0, 4], [5.0, 6]))
        for id in ("3D surface", "3D points")
            e = entry(id)
            @test e[:kind] == "plot3d" && endswith(e[:file], ".json")
        end
        @test only(traces("3D surface"))["z"] == [[1, 2], [3, 4]]
        t = only(traces("3D points"))
        @test t["type"] == "scatter3d" && t["x"] == [1, 2] && t["y"] == [3, 4] && t["z"] == [5, 6]
        # The dashboard draws every kind Julia can publish: Cesium for a scene, Plotly otherwise.
        html = read(joinpath(PUB.assets_dir(), "dashboard.html"), String)
        @test occursin("p.kind === 'scene'", html)
    end

    @testset "error bars arrive under error_x, error_y and error_z" begin
        quiet(() -> xyplot("Errors", [1.0, 2, 3], [1.0, 2, 3]; error_y_type = "data", error_y_array = [0.1, 0.2, 0.3]))
        quiet(() -> bar("Errors bar", ["a", "b"], [1.0, 2]; error_y_type = "data", error_y_array = [0.1, 0.2]))
        quiet(() -> scatter3d("Errors 3D", [1.0, 2], [1.0, 2], [1.0, 2]; error_z_type = "data", error_z_array = [0.1, 0.2]))
        for (id, group, values) in (("Errors", "error_y", [0.1, 0.2, 0.3]), ("Errors bar", "error_y", [0.1, 0.2]),
                                    ("Errors 3D", "error_z", [0.1, 0.2]))
            t = only(traces(id))
            @test t[group] == Dict("type" => "data", "array" => values) && !haskey(t, "error")
        end
        # Above 4,000 points the scatter is WebGL, and its error bars arrive the same way.
        n = 5000
        quiet(() -> xyplot("Errors gl", collect(1.0:n), ones(n); error_y_type = "constant", error_y_value = 0.1))
        t = only(traces("Errors gl"))
        @test t["type"] == "scattergl" && t["error_y"] == Dict("type" => "constant", "value" => 0.1)
        # One array per trace when a call draws several.
        quiet(() -> xyplot("Errors two", [1.0, 2], [[1.0, 2], [3.0, 4]]; error_y_array = [[0.1, 0.2], [0.3, 0.4]]))
        @test [t["error_y"]["array"] for t in traces("Errors two")] == [[0.1, 0.2], [0.3, 0.4]]
    end

    @testset "documented attributes pass without a warning" begin
        quiet(() -> bar("Quiet bar", ["a", "b"], [1.0, 2]))
        quiet(() -> panel!("Quiet bar"; barmode = "group", plot_bgcolor = "white", yaxis_title = "km/s"))
        quiet(() -> histogram("Quiet hist", [1.0, 2, 3]; xbins_size = 0.5, histnorm = "probability"))
        quiet(() -> contour("Quiet contour", [1.0, 2], [1.0, 2], [1.0 2; 3 4]; contours_showlabels = true, colorbar_title = "C3"))
        quiet(() -> scattergeo("Quiet geo", [0.0], [0.0]))
        quiet(() -> panel!("Quiet geo"; geo_projection_type = "orthographic"))
        # and a misspelled one still warns
        @test_logs (:warn, r"barmod") panel!("Quiet bar"; barmod = "group")
    end

    @testset "colors: the plot palette, on the attribute each trace draws with" begin
        P = PUB.PALETTE
        for (f, f!, id, args) in ((bar, bar!, "Color bar", (["a", "b"], [1.0, 2])),
                                  (histogram, histogram!, "Color hist", ([1.0, 2, 3],)))
            f(id, args...); f!(id, args...)
            @test [t["marker"]["color"] for t in traces(id)] == P[1:2]
        end
        for (f, f!, id, args) in ((xyplot, xyplot!, "Color xy", ([1.0, 2], [1.0, 2])),
                                  (scatterpolar, scatterpolar!, "Color polar", ([0.0, 90], [1.0, 2])),
                                  (scattergeo, scattergeo!, "Color geo", ([0.0, 1], [0.0, 1])))
            f(id, args...); f!(id, args...)
            @test [t["line"]["color"] for t in traces(id)] == P[1:2]            # the `!` form continues the cycle
        end
    end

    @testset "band and its ! form" begin
        x, lo, hi = [1.0, 2, 3], [0.0, 0, 0], [1.0, 2, 3]
        band("Band", x, lo, hi)
        t = only(traces("Band"))
        @test t["x"] == vcat(x, reverse(x)) && t["y"] == vcat(hi, reverse(lo))
        band!("Band", x, lo .- 1, hi .+ 1)
        @test length(traces("Band")) == 2
    end

    @testset "a one-entry name vector names the one trace" begin
        xyplot("Named", [1.0, 2], [1.0, 2]; name = ["altitude"])
        @test only(traces("Named"))["name"] == "altitude"
        @test_throws "`name` has 2 entries" xyplot("Named", [1.0, 2], [1.0, 2]; name = ["a", "b"])
    end

    @testset "the empty dashboard names a function that exists" begin
        html = read(joinpath(PUB.assets_dir(), "dashboard.html"), String)
        m = match(r"Nothing plotted yet\. Call <code>(\w+)\(", html)
        @test m !== nothing && isdefined(EpicycleIO, Symbol(m[1])) && Symbol(m[1]) in names(EpicycleIO)
    end

    @testset "panels survive close_dashboard and a restart" begin
        xyplot("Restart A", [1.0, 2], [1.0, 2]); xyplot("Restart B", [1.0, 2], [2.0, 1])
        PUB.flush_publishes()
        PUB.close_dashboard()
        xyplot("Restart C", [1.0, 2], [3.0, 3])
        PUB.flush_publishes()
        files = readdir(PUB.assets_dir())
        for id in ("Restart A", "Restart B", "Restart C")
            @test any(e -> e[:id] == PUB.panel_id(id), PUB.manifest()[:panels])
            @test "panel_$(PUB.panel_id(id)).json" in files
        end
    end
end
