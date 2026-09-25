# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Checks the PlotlyBase behaviours EpicycleIO's design assumes.
#
# The spec was written from Plotly's published schema, not from the installed package,
# and carried a red flag saying so. This is that flag being discharged. Run it after
# bumping PlotlyBase.
#
#   julia --project=<environment> EpicycleIO/test/plotlybase_assumptions.jl

using PlotlyBase

pass = fail = 0

function check(f::Function, label)     # f first, so `check(label) do … end` works
    global pass, fail
    try
        result = f()
        if result === true || result isa AbstractString
            pass += 1
            println("  ok    ", label, result isa AbstractString ? "  ($result)" : "")
        else
            fail += 1
            println("  FAIL  ", label, "  got ", result)
        end
    catch e
        fail += 1
        println("  ERROR ", label, "  ", first(sprint(showerror, e), 120))
    end
end

# Reach the raw attribute dictionary a trace or layout carries.
raw(x) = getfield(x, :fields)

nested(d, path...) = begin
    cur = d
    for k in path
        cur isa AbstractDict || return nothing
        haskey(cur, k) || return nothing
        cur = cur[k]
    end
    cur
end

println("PlotlyBase v", pkgversion(PlotlyBase))
println("\n── Underscore flattening ──")

check("line_width nests to line.width") do
    nested(raw(scatter(x = [1, 2], y = [3, 4], line_width = 2)), :line, :width) == 2
end

check("line_color nests to line.color") do
    nested(raw(scatter(x = [1, 2], y = [3, 4], line_color = "cyan")), :line, :color) == "cyan"
end

check("marker_size nests to marker.size") do
    nested(raw(scatter(x = [1, 2], y = [3, 4], marker_size = 8)), :marker, :size) == 8
end

check("layout polar attributes nest three deep") do
    l = raw(Layout(polar_angularaxis_direction = "clockwise",
                   polar_angularaxis_rotation  = 90,
                   polar_radialaxis_range      = [90, 0]))
    nested(l, :polar, :angularaxis, :direction) == "clockwise" &&
        nested(l, :polar, :radialaxis, :range) == [90, 0]
end

println("\n── Trace constructors the wrappers rely on ──")

for (name, f) in ("scatter"      => () -> scatter(x = [1, 2], y = [3, 4]),
                  "scatterpolar" => () -> scatterpolar(r = [1, 2], theta = [0.1, 0.2],
                                                       thetaunit = "radians"),
                  "bar"          => () -> bar(x = ["a", "b"], y = [1, 2]),
                  "histogram"    => () -> histogram(x = [1, 2, 3]),
                  "contour"      => () -> contour(x = [1, 2], y = [3, 4], z = [1 2; 3 4]),
                  "heatmap"      => () -> heatmap(x = [1, 2], y = [3, 4], z = [1 2; 3 4]),
                  "scattergeo"   => () -> scattergeo(lon = [1, 2], lat = [3, 4]),
                  "scatter3d"    => () -> scatter3d(x = [1], y = [2], z = [3]),
                  "surface"      => () -> surface(z = [1 2; 3 4]),
                  "scattergl"    => () -> scattergl(x = [1, 2], y = [3, 4]),
                  "violin"       => () -> violin(y = [1, 2, 3]))
    check(() -> (f(); true), "$name constructs")
end

println("\n── The shipped schema ──")

check("get_plotschema is callable") do
    s = PlotlyBase.get_plotschema()
    string(typeof(s))
end

check("schema top-level keys") do
    d = getfield(PlotlyBase.get_plotschema(), :fields)
    string(sort(collect(keys(d))))
end

check("arrayOk is reachable for scatter attributes") do
    d = getfield(PlotlyBase.get_plotschema(), :fields)
    attrs = d[:traces][:scatter][:attributes]
    nm = get(attrs[:name], :arrayOk, missing)
    ms = get(attrs[:marker][:size], :arrayOk, missing)
    "name arrayOk=$nm  |  marker.size arrayOk=$ms"
end

check("unknown attributes: dropped or carried?") do
    f = raw(scatter(x = [1], y = [2], nonsense_attribute = 1))
    haskey(f, :nonsense_attribute) ? "carried through" : "DROPPED silently by PlotlyBase"
end

println("\n", fail == 0 ? "all $pass checks passed" : "$pass passed, $fail FAILED")
exit(fail == 0 ? 0 : 1)
