# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# Draws time_scales.svg for the AstroEpochs guide from AstroEpochs itself. It needs nothing else, and
# the numbers are the package's own. Rerun it when the conversions change:
#
#   julia --project=<env> AstroEpochs/docs/src/assets/make_time_scales.jl

using AstroEpochs

const INK, MUTED, GRID, MARK = "#1f2933", "#52606d", "#d9dee3", "#1f6f8b"
const PANEL_W, PANEL_H, LEFT, TOP, GAP_X, GAP_Y = 400, 200, 70, 44, 90, 84

year_jd(y) = 2451545.0 + (y - 2000.0) * 365.25
date_at(y, scale) = Time(year_jd(y), scale, JD())
parts_minus(a, b) = ((a.jd1 - b.jd1) + (a.jd2 - b.jd2)) * 86400.0   # seconds between the numbers

# Each panel's numbers.
const LONG = range(1972.01, 2030.0; length = 240)

function leap_series()
    ys = [tai_minus_utc(date_at(y, UTC())) for y in range(1972.01, 2030.0; length = 2320)]
    return collect(range(1972.01, 2030.0; length = 2320)), ys
end

tcg_series() = (collect(LONG), [parts_minus(date_at(y, TT()).tcg, date_at(y, TT())) for y in LONG])
tcb_series() = (collect(LONG), [parts_minus(date_at(y, TDB()).tcb, date_at(y, TDB())) for y in LONG])

function tdb_series()
    xs = collect(range(2020.0, 2024.0; length = 480))
    return xs, [tdb_minus_tt(date_at(y, TT())) * 1e3 for y in xs]
end

function panel(io, col, row, title, xs, ys, xticks, yticks)
    x0 = LEFT + col * (PANEL_W + GAP_X)
    y0 = TOP + row * (PANEL_H + GAP_Y)
    px(x) = x0 + (x - first(xticks)) / (last(xticks) - first(xticks)) * PANEL_W
    py(y) = y0 + PANEL_H - (y - first(yticks)) / (last(yticks) - first(yticks)) * PANEL_H
    println(io, """<text x="$x0" y="$(y0 - 14)" font-size="16" font-weight="600" fill="$INK">$title</text>""")
    for t in yticks
        println(io, """<line x1="$x0" x2="$(x0 + PANEL_W)" y1="$(py(t))" y2="$(py(t))" stroke="$GRID" stroke-width="1"/>""")
        println(io, """<text x="$(x0 - 10)" y="$(py(t) + 4)" font-size="13" text-anchor="end" fill="$MUTED">$(t)</text>""")
    end
    for t in xticks
        println(io, """<text x="$(px(t))" y="$(y0 + PANEL_H + 22)" font-size="13" text-anchor="middle" fill="$MUTED">$(Int(t))</text>""")
    end
    points = join(("$(round(px(x); digits = 2)),$(round(py(y); digits = 2))" for (x, y) in zip(xs, ys)), " ")
    println(io, """<polyline points="$points" fill="none" stroke="$MARK" stroke-width="2" stroke-linejoin="round"/>""")
end

function draw(path)
    width  = LEFT + 2PANEL_W + GAP_X + 30
    height = TOP + 2PANEL_H + GAP_Y + 44
    open(path, "w") do io
        println(io, """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 $width $height" font-family="Helvetica, Arial, sans-serif" role="img" aria-label="Four panels against time: TAI minus UTC steps from 10 to 37 seconds, TDB minus TT oscillates within about 1.7 milliseconds each year, and TCG minus TT and TCB minus TDB each drift steadily.">""")
        println(io, """<rect width="$width" height="$height" fill="#ffffff"/>""")
        decades = 1970.0:20.0:2030.0
        panel(io, 0, 0, "TAI − UTC, s: a step at each leap second", leap_series()..., decades, 0.0:10.0:40.0)
        panel(io, 1, 0, "TDB − TT, ms: a yearly cycle (2020 to 2024)", tdb_series()..., 2020.0:1.0:2024.0, -2.0:1.0:2.0)
        panel(io, 0, 1, "TCG − TT, s: a steady drift", tcg_series()..., decades, -0.5:0.5:1.5)
        panel(io, 1, 1, "TCB − TDB, s: a faster steady drift", tcb_series()..., decades, -10.0:10.0:30.0)
        println(io, "</svg>")
    end
end

draw(joinpath(@__DIR__, "time_scales.svg"))
