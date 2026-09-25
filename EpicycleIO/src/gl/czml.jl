# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Ported from the Cesium prototype unchanged except for what the panel system now owns.
# The CZML itself — sampled positions, path lead and trail, the clock packet, the
# inertial reference frame — took time to get right and is deliberately untouched.

# CZML serialization. Internal — the public interface is `orbitview`.
#
# Emits packets that Cesium's CzmlDataSource can load:
#   1. document packet: clock
#   2. one path entity per trajectory segment (its own color arc)
#   3. one marker entity carrying an interval-tagged position property so
#      Lagrange interpolation never crosses a maneuver boundary
#   4. one point entity per maneuver, positioned at the discontinuity

"Internal representation of one continuous propagation arc."
struct _Segment
    times_s::Vector{Float64}                   # seconds since trajectory epoch
    positions_km::Vector{NTuple{3, Float64}}
    name::String
    color::NTuple{4, Int}
    linewidth::Float64
end

"Internal representation of an impulsive maneuver event."
struct _Maneuver
    name::String
    time_s::Float64
    dv_kms::NTuple{3, Float64}
end

# Convert (segments, maneuvers, epoch, ...) into a CZML packet array.
function _build_czml(
    segments::Vector{_Segment}, maneuvers::Vector{_Maneuver},
    epoch_iso::AbstractString;
    id::AbstractString,
    name::AbstractString,
    marker_color::NTuple{4, Int},
    marker_pixel_size::Real,
    multiplier::Real,
    inertial::Bool,
    pos_scale::Real = 1000.0,   # km → m
)
    isempty(segments) && throw(ArgumentError(
        "orbitview: the trajectory has no arcs to draw. An arc needs at least 2 samples, and " *
        "every recorded segment holds fewer. Propagate the spacecraft before drawing it."))

    t_min = minimum(minimum(s.times_s) for s in segments)
    t_max = maximum(maximum(s.times_s) for s in segments)
    epoch_start_iso = _iso_shift(epoch_iso, t_min)
    epoch_end_iso   = _iso_shift(epoch_iso, t_max)
    availability    = string(epoch_start_iso, "/", epoch_end_iso)
    total_duration  = t_max - t_min

    packets = Any[
        Dict(
            "id"      => "document",
            "name"    => "Epicycle trajectory",
            "version" => "1.0",
            "clock"   => Dict(
                "interval"    => availability,
                "currentTime" => epoch_start_iso,
                "multiplier"  => multiplier,
                "range"       => "LOOP_STOP",
                "step"        => "SYSTEM_CLOCK_MULTIPLIER",
            ),
        ),
    ]

    for (i, seg) in enumerate(segments)
        length(seg.times_s) >= 2 || throw(ArgumentError(
            "orbitview: arc $i must have at least 2 samples to interpolate between; " *
            "got $(length(seg.times_s))."))

        cart = Float64[]
        for (t, p) in zip(seg.times_s, seg.positions_km)
            push!(cart, t, p[1] * pos_scale, p[2] * pos_scale, p[3] * pos_scale)
        end

        seg_position = Dict{String, Any}(
            "epoch"                     => epoch_iso,
            "cartesian"                 => cart,
            "interpolationAlgorithm"    => "LAGRANGE",
            "interpolationDegree"       => min(5, length(seg.times_s) - 1),
            "forwardExtrapolationType"  => "NONE",
            "backwardExtrapolationType" => "NONE",
        )
        inertial && (seg_position["referenceFrame"] = "INERTIAL")

        push!(packets, Dict(
            "id"           => "$id-seg-$i",
            "name"         => isempty(seg.name) ? "segment $i" : seg.name,
            "availability" => availability,
            "position"     => seg_position,
            "path"         => Dict(
                "show"       => true,
                "leadTime"   => total_duration,
                "trailTime"  => total_duration,
                "width"      => seg.linewidth,
                "resolution" => 60.0,
                "material"   => Dict("solidColor" =>
                    Dict("color" => Dict("rgba" => collect(seg.color)))),
            ),
        ))
    end

    # Marker position as an ARRAY of interval-tagged position properties —
    # Lagrange interpolation stays inside each segment, so the marker never
    # smooths across a velocity discontinuity at a maneuver.
    marker_position = Any[]
    all_pairs = Tuple{Float64, NTuple{3, Float64}}[]
    for seg in segments, (t, p) in zip(seg.times_s, seg.positions_km)
        push!(all_pairs, (t, p))
    end
    sort!(all_pairs, by = p -> p[1])

    for seg in segments
        seg_cart = Float64[]
        for (t, p) in zip(seg.times_s, seg.positions_km)
            push!(seg_cart, t, p[1] * pos_scale, p[2] * pos_scale, p[3] * pos_scale)
        end
        entry = Dict{String, Any}(
            "interval"                  => string(_iso_shift(epoch_iso, first(seg.times_s)),
                                                  "/",
                                                  _iso_shift(epoch_iso, last(seg.times_s))),
            "epoch"                     => epoch_iso,
            "cartesian"                 => seg_cart,
            "interpolationAlgorithm"    => "LAGRANGE",
            "interpolationDegree"       => min(5, length(seg.times_s) - 1),
            "forwardExtrapolationType"  => "NONE",
            "backwardExtrapolationType" => "NONE",
        )
        inertial && (entry["referenceFrame"] = "INERTIAL")
        push!(marker_position, entry)
    end

    push!(packets, Dict(
        "id"           => id,
        "name"         => name,
        "availability" => availability,
        "position"     => marker_position,
        "point"        => Dict(
            "pixelSize"    => marker_pixel_size,
            "color"        => Dict("rgba" => collect(marker_color)),
            "outlineColor" => Dict("rgba" => [0, 0, 0, 255]),
            "outlineWidth" => 2,
        ),
        "label" => Dict(
            "text"             => name,
            "font"             => "14px sans-serif",
            "horizontalOrigin" => "CENTER",
            "verticalOrigin"   => "BOTTOM",
            "pixelOffset"      => Dict("cartesian2" => [0, -18]),
            "fillColor"        => Dict("rgba" => [255, 255, 255, 255]),
            "showBackground"   => true,
            "backgroundColor"  => Dict("rgba" => [0, 0, 0, 128]),
        ),
    ))

    # The index is part of the id, as it is for the arcs. Naming alone is not
    # enough: solve_trajectory! names every marker it records "maneuver", so
    # three burns all produced the same id and Cesium — which merges packets by
    # id — drew one point where there should have been three.
    for (i, m) in enumerate(maneuvers)
        pos = _interpolate(all_pairs, m.time_s)
        pos === nothing && continue
        maneuver_pos = Dict{String, Any}(
            "cartesian" => [pos[1] * pos_scale, pos[2] * pos_scale, pos[3] * pos_scale],
        )
        inertial && (maneuver_pos["referenceFrame"] = "INERTIAL")
        dv_mag = sqrt(m.dv_kms[1]^2 + m.dv_kms[2]^2 + m.dv_kms[3]^2)
        push!(packets, Dict(
            "id"       => "$id-mnv-$i-$(m.name)",
            "name"     => m.name,
            "position" => maneuver_pos,
            "point"    => Dict(
                "pixelSize"    => 12,
                "color"        => Dict("rgba" => [255, 0, 128, 255]),
                "outlineColor" => Dict("rgba" => [255, 255, 255, 255]),
                "outlineWidth" => 2,
            ),
            "label" => Dict(
                "text"            => "$(m.name)  ΔV=$(round(dv_mag, digits=3)) km/s",
                "font"            => "11px sans-serif",
                "pixelOffset"     => Dict("cartesian2" => [14, -12]),
                "fillColor"       => Dict("rgba" => [255, 200, 200, 255]),
                "showBackground"  => true,
                "backgroundColor" => Dict("rgba" => [40, 0, 20, 200]),
            ),
        ))
    end

    return packets
end


# Linearly interpolate the state at `t_query` from a time-sorted list of pairs.
function _interpolate(pairs::Vector{Tuple{Float64, NTuple{3, Float64}}}, t_query::Real)
    isempty(pairs) && return nothing
    t_query < pairs[1][1]   && return nothing
    t_query > pairs[end][1] && return nothing
    for i in 2:length(pairs)
        t1, p1 = pairs[i-1]
        t2, p2 = pairs[i]
        if t_query <= t2
            α = t2 == t1 ? 0.0 : (t_query - t1) / (t2 - t1)
            return (p1[1] + α * (p2[1] - p1[1]),
                    p1[2] + α * (p2[2] - p1[2]),
                    p1[3] + α * (p2[3] - p1[3]))
        end
    end
    return pairs[end][2]
end

# Shift an ISO 8601 UTC string by `delta_s` seconds (returns ISO with trailing Z).
function _iso_shift(iso::AbstractString, delta_s::Real)
    clean = replace(replace(iso, "Z" => ""), "z" => "")
    dt = DateTime(clean)
    dt2 = dt + Millisecond(round(Int, delta_s * 1000))
    return string(dt2) * "Z"
end

# User-passed panel specs coerced to a JSON-friendly Vector{Dict}.
