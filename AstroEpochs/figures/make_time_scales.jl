# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# Draws time_scales.svg for the AstroEpochs guide from AstroEpochs itself. It needs nothing else, and
# the numbers are the package's own. Rerun it when the conversions change:
#
#   julia --project=<env> AstroEpochs/figures/make_time_scales.jl
#
# The form follows Figure 1 of Seidelmann and Fukushima, "Why new time scales?", Astronomy and
# Astrophysics 265, 833 (1992): every scale against TAI on one axis, with the periodic terms of TDB
# and TCB magnified by 100. UTC starts at 1972, where AstroEpochs starts modeling it, and stops at the
# last year this figure was drawn for, since later leap seconds are unknown. The names on the lines
# stand in for a legend, since TT, TDB and TCG run close together and color alone would not tell them.

using AstroEpochs

const INK, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#dcdbd5", "#fcfcfb"
const COLORS = (TAI = "#2a78d6", TT = "#eb6834", TDB = "#1baf7a", TCG = "#eda100", TCB = "#e87ba4", UTC = "#008300")
const WIDTH, HEIGHT = 860, 510
const LEFT, RIGHT, TOP, BOTTOM = 70, 24, 24, 52

const FIRST_YEAR, LAST_YEAR, LAST_UTC_YEAR = 1950.0, 2030.0, 2026.0
const MAGNIFY_PERIODIC = 100.0
const TT_MINUS_TAI = 32.184

year_jd(y) = 2451545.0 + (y - 2000.0) * 365.25
parts_minus(a, b) = ((a.jd1 - b.jd1) + (a.jd2 - b.jd2)) * 86400.0     # seconds between the two dates

tai_at(y) = Time(year_jd(y), TAI(), JD())
periodic(y) = tdb_minus_tt(tai_at(y))

tt_offset(y)  = parts_minus(tai_at(y).tt, tai_at(y))
tcg_offset(y) = parts_minus(tai_at(y).tcg, tai_at(y))
tdb_offset(y) = TT_MINUS_TAI + MAGNIFY_PERIODIC * periodic(y)
tcb_offset(y) = parts_minus(tai_at(y).tcb, tai_at(y).tdb) + TT_MINUS_TAI + MAGNIFY_PERIODIC * periodic(y)
utc_offset(y) = -tai_minus_utc(Time(year_jd(y), UTC(), JD()))

function draw(path)
    plot_w, plot_h = WIDTH - LEFT - RIGHT, HEIGHT - TOP - BOTTOM
    ymin, ymax = -40.0, 60.0
    px(x) = round(LEFT + (x - FIRST_YEAR) / (LAST_YEAR - FIRST_YEAR) * plot_w; digits = 1)
    py(y) = round(TOP + plot_h - (y - ymin) / (ymax - ymin) * plot_h; digits = 1)

    fine = collect(range(FIRST_YEAR, LAST_YEAR; length = 1600))
    coarse = [FIRST_YEAR, LAST_YEAR]
    leap = collect(range(1972.01, LAST_UTC_YEAR; length = 2160))
    series = ((:TAI, coarse, y -> 0.0, ""),
              (:TDB, fine, tdb_offset, ""),
              (:TT, coarse, tt_offset, ""),      # over the TDB curve, whose middle it runs through
              (:TCG, coarse, tcg_offset, "7 4"),
              (:TCB, fine, tcb_offset, ""),
              (:UTC, leap, utc_offset, ""))

    open(path, "w") do io
        println(io, """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 $WIDTH $HEIGHT" font-family="Helvetica, Arial, sans-serif" role="img" aria-label="Each time scale's offset from TAI from 1950 to 2030, with the periodic terms of TDB and TCB magnified by 100. TT is a constant 32.184 seconds. TDB and TCG stay close to TT. TCB drifts from about 19 to about 58 seconds. UTC steps down from minus 10 seconds at 1972 to minus 37 seconds, and ends at 2026.">""")
        println(io, """<rect width="$WIDTH" height="$HEIGHT" fill="$SURFACE"/>""")

        for y in ymin:20.0:ymax
            println(io, """<line x1="$LEFT" x2="$(LEFT + plot_w)" y1="$(py(y))" y2="$(py(y))" stroke="$GRID" stroke-width="1"/>""")
            println(io, """<text x="$(LEFT - 10)" y="$(py(y) + 4)" font-size="13" text-anchor="end" fill="$MUTED">$(Int(y))</text>""")
        end
        for x in FIRST_YEAR:10.0:LAST_YEAR
            println(io, """<text x="$(px(x))" y="$(TOP + plot_h + 22)" font-size="13" text-anchor="middle" fill="$MUTED">$(Int(x))</text>""")
        end
        println(io, """<text x="$(LEFT + plot_w / 2)" y="$(HEIGHT - 8)" font-size="13" text-anchor="middle" fill="$MUTED">Year</text>""")
        println(io, """<text transform="translate(18 $(TOP + plot_h / 2)) rotate(-90)" font-size="13" text-anchor="middle" fill="$MUTED">Offset from TAI, s</text>""")

        for (name, xs, f, dash) in series
            points = join(("$(px(x)),$(py(f(x)))" for x in xs), " ")
            dashes = isempty(dash) ? "" : """ stroke-dasharray="$dash\""""
            println(io, """<polyline points="$points" fill="none" stroke="$(COLORS[name])" stroke-width="2" stroke-linejoin="round"$dashes/>""")
        end

        println(io, """<text x="$(px(1952.0))" y="$(py(13.0))" font-size="13" fill="$MUTED">Periodic terms of TDB and TCB × 100</text>""")
        println(io, """<text x="$(px(LAST_UTC_YEAR) + 6)" y="$(py(utc_offset(LAST_UTC_YEAR)) + 4)" font-size="13" fill="$MUTED">$(Int(LAST_UTC_YEAR))</text>""")
        for (name, year, value, label_year, label_value) in (("TCB", 2006.0, tcb_offset(2006.0), 1999.0, 52.0),
                                                             ("TDB", 1962.0, tdb_offset(1962.0), 1956.0, 43.0),
                                                             ("TCG", 2018.0, tcg_offset(2018.0), 2022.0, 42.0),
                                                             ("TT", 1993.0, TT_MINUS_TAI, 1999.0, 21.0),
                                                             ("TAI", 2010.0, 0.0, 2014.0, 7.0),
                                                             ("UTC", 1990.0, utc_offset(1990.0), 1982.0, -33.0))
            println(io, """<line x1="$(px(label_year))" y1="$(py(label_value))" x2="$(px(year))" y2="$(py(value))" stroke="$MUTED" stroke-width="1"/>""")
            println(io, """<text x="$(px(label_year))" y="$(py(label_value) + (label_value > value ? -5 : 15))" font-size="14" text-anchor="middle" fill="$INK">$name</text>""")
        end
        println(io, "</svg>")
    end
end

draw(joinpath(@__DIR__, "..", "docs", "src", "assets", "time_scales.svg"))
