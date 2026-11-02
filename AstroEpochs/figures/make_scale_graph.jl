# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# Draws scale_graph.svg for the AstroEpochs guide: which time scales convert directly, and what
# separates them. It stops with an error if the links drawn here are not the links AstroEpochs has, so
# the figure cannot drift from the code. Rerun it when the conversion graph changes:
#
#   julia --project=<env> AstroEpochs/figures/make_scale_graph.jl

using AstroEpochs

const INK, MUTED, SURFACE, LINK = "#0b0b0b", "#52514e", "#fcfcfb", "#52514e"
const WIDTH, HEIGHT = 860, 310
const NODE_W, NODE_H = 92, 58

# The colors match time_scales.svg, so a scale looks the same in both figures.
const NODES = (tai = (name = "TAI", about = "atomic", x = 248, y = 90, color = "#2a78d6"),
               utc = (name = "UTC", about = "civil", x = 70, y = 90, color = "#008300"),
               tt  = (name = "TT", about = "terrestrial", x = 426, y = 90, color = "#eb6834"),
               tdb = (name = "TDB", about = "solar system", x = 604, y = 90, color = "#1baf7a"),
               tcb = (name = "TCB", about = "barycentric", x = 782, y = 90, color = "#e87ba4"),
               tcg = (name = "TCG", about = "geocentric", x = 426, y = 256, color = "#eda100"))

# Each link, and the two lines that say what separates its scales.
const LINKS = ((:utc, :tai, "leap seconds", "10 s to 37 s"),
               (:tai, :tt, "constant", "32.184 s"),
               (:tt, :tdb, "yearly cycle", "about 1.7 ms"),
               (:tdb, :tcb, "steady rate", "0.49 s a year"),
               (:tt, :tcg, "steady rate", "0.022 s a year"))

function check_against_code()
    drawn = Set(Iterators.flatten(((a, b), (b, a)) for (a, b, _, _) in LINKS))
    actual = Set(keys(AstroEpochs.SCALE_TRANSFORMS))
    drawn == actual || error("the figure's links are not AstroEpochs' links: $(symdiff(drawn, actual))")
    return nothing
end

function draw(path)
    check_against_code()
    open(path, "w") do io
        println(io, """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 $WIDTH $HEIGHT" font-family="Helvetica, Arial, sans-serif" role="img" aria-label="The conversion graph. UTC joins TAI by leap seconds. TAI joins TT by a constant 32.184 seconds. TT joins TDB by a yearly cycle of about 1.7 milliseconds. TDB joins TCB by a steady rate. TT also joins TCG by a steady rate.">""")
        println(io, """<defs><marker id="end" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0 0 L10 5 L0 10 z" fill="$LINK"/></marker></defs>""")
        println(io, """<rect width="$WIDTH" height="$HEIGHT" fill="$SURFACE"/>""")
        for (a, b, line1, line2) in LINKS
            from, to = NODES[a], NODES[b]
            if from.y == to.y
                x1, x2 = from.x + NODE_W / 2 + 3, to.x - NODE_W / 2 - 3
                mid = (x1 + x2) / 2
                println(io, """<line x1="$x1" x2="$x2" y1="$(from.y)" y2="$(to.y)" stroke="$LINK" stroke-width="1.5" marker-start="url(#end)" marker-end="url(#end)"/>""")
                println(io, """<text x="$mid" y="$(from.y - 24)" font-size="12" text-anchor="middle" fill="$INK">$line1</text>""")
                println(io, """<text x="$mid" y="$(from.y - 9)" font-size="12" text-anchor="middle" fill="$MUTED">$line2</text>""")
            else
                y1, y2 = from.y + NODE_H / 2 + 3, to.y - NODE_H / 2 - 3
                mid = (y1 + y2) / 2
                println(io, """<line x1="$(from.x)" x2="$(to.x)" y1="$y1" y2="$y2" stroke="$LINK" stroke-width="1.5" marker-start="url(#end)" marker-end="url(#end)"/>""")
                println(io, """<text x="$(from.x + 12)" y="$(mid - 3)" font-size="12" fill="$INK">$line1</text>""")
                println(io, """<text x="$(from.x + 12)" y="$(mid + 12)" font-size="12" fill="$MUTED">$line2</text>""")
            end
        end
        for node in NODES
            println(io, """<rect x="$(node.x - NODE_W / 2)" y="$(node.y - NODE_H / 2)" width="$NODE_W" height="$NODE_H" rx="4" fill="#ffffff" stroke="$(node.color)" stroke-width="2.5"/>""")
            println(io, """<text x="$(node.x)" y="$(node.y - 3)" font-size="16" font-weight="600" text-anchor="middle" fill="$INK">$(node.name)</text>""")
            println(io, """<text x="$(node.x)" y="$(node.y + 15)" font-size="12" text-anchor="middle" fill="$MUTED">$(node.about)</text>""")
        end
        println(io, "</svg>")
    end
end

draw(joinpath(@__DIR__, "..", "docs", "src", "assets", "scale_graph.svg"))
