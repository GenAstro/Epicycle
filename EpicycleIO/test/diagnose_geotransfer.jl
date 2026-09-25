# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Does a real multi-burn trajectory survive the trip into `orbitview`?
#
# Ex_GeoTransfer is the case that matters: three impulsive maneuvers (TOI, MCC,
# MOI) recorded by `solve_trajectory!` rather than pushed in by hand. Two things
# are being asked.
#
#   1. What does the solver actually write into the history? Segment count,
#      names, and which of them are 1-sample maneuver markers. Everything the
#      3D path does is downstream of that, and it has never been looked at.
#
#   2. Would the maneuver-id collision bite this case? `czml.jl` builds a
#      maneuver packet id from its name alone, so two burns sharing a name would
#      collide into one Cesium entity. Detection is suppressed today, so nothing
#      collides now — this asks whether it would if it were switched back on.
#
#   julia --project=<environment> EpicycleIO/test/diagnose_geotransfer.jl

using Epicycle
using EpicycleIO

const IO_ = EpicycleIO
IO_.auto_open!(false)

const EXAMPLE = normpath(joinpath(@__DIR__, "..", "..", "Epicycle", "examples", "Ex_GeoTransfer.jl"))

# Run the example, but stop before its Makie view — we want the spacecraft, not
# a window. Everything above that line is the solve.
#
# This runs at top level rather than inside a function: `include_string` creates
# the `sat` binding, and a function already running cannot see a binding that new.
println("Solving Ex_GeoTransfer …")
let src = read(EXAMPLE, String)
    cut = findfirst("# Visualize with iterations", src)
    cut === nothing && error("the example no longer has the visualise block; check the tail")
    include_string(Main, src[1:first(cut)-1], EXAMPLE)
end

function main(sat)
    hsegs = sat.history.segments
    println("\n── what the solver recorded ──")
    println("spacecraft: ", sat.name)
    println("segments:   ", length(hsegs))
    for (i, s) in enumerate(hsegs)
        kind = length(s.times) == 1 ? "MARKER" : "arc"
        println("  ", lpad(i, 2), "  ", rpad(kind, 7),
                lpad(length(s.times), 5), " samples   name=\"", s.name, "\"")
    end

    markers = [s for s in hsegs if length(s.times) == 1]
    names   = [s.name for s in markers]
    println("\n── the id-collision question ──")
    println("1-sample markers: ", length(markers))
    println("their names:      ", names)
    if isempty(markers)
        println("VERDICT: the solver records no 1-sample markers, so detection")
        println("         method 1 would find nothing here at all.")
    elseif length(unique(names)) == length(names)
        println("VERDICT: names are distinct here, so the name alone would have been")
        println("         enough for this case. It is not in general.")
    else
        dupes = [n for n in unique(names) if count(==(n), names) > 1]
        println("VERDICT: names REPEAT ", dupes, ". Building the packet id from the")
        println("         name alone gave every burn the same id, and Cesium merges")
        println("         packets by id — three burns, one point. The id now carries")
        println("         the index too, so this case draws all three.")
    end

    println("\n── what orbitview draws today ──")
    epoch_iso, segments, maneuvers = IO_._extract(sat)
    println("epoch (UTC):    ", epoch_iso)
    println("arcs drawn:     ", length(segments))
    println("maneuvers drawn:", length(maneuvers))
    for s in segments
        println("  arc  \"", rpad(s.name, 22), "\"  ", lpad(length(s.times_s), 5),
                " samples  colour=", s.color)
    end

    # The ΔV is read from the velocity jump across each marker, never handed in.
    # Comparing it with what the solver printed above is the check that the
    # extraction is looking at the right pair of arcs.
    println("\n── the burns, differenced from the history ──")
    println("(compare with the solver's Total ΔV for toi_v, mcc_vn, moi_v above)")
    for (i, m) in enumerate(maneuvers)
        mag = sqrt(sum(abs2, m.dv_kms))
        println("  ", i, "  t=", lpad(round(m.time_s, digits = 1), 10), " s",
                "   |ΔV|=", round(mag, digits = 6), " km/s")
    end

    packets = IO_.build_orbitview(sat)
    ids = [p["id"] for p in packets]
    println("\nCZML packets:   ", length(packets))
    println("ids:            ", ids)
    println("all ids unique: ", length(unique(ids)) == length(ids))

    # And it must survive being written out, which is how it reaches the page.
    orbitview("GEO transfer", sat)
    p = IO_.get_panel("GEO transfer")
    println("\npanel family:   ", p.family)
    println("payload bytes:  ", length(IO_.payload(p)))
    println("\nDone.")
    return nothing
end

main(sat)
