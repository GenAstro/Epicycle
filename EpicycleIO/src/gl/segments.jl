# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Extract trajectory segments and impulsive-maneuver events from a Spacecraft's
# history. This is the "no user boilerplate" layer — turns
# `sat.history.segments` (any mixture of multi-sample propagation arcs and
# 1-sample maneuver markers) into the internal representation the CZML writer
# consumes.
#
# An impulsive maneuver is found as a 1-sample HistorySegment, which is how
# solve_trajectory! records TOI, MCC and the rest. The maneuver's time is that
# sample's time and its name is the segment's; the ΔV is the velocity jump
# between the propagation segments either side, so nothing has to be passed in.
#
# A burn recorded some other way — two arcs meeting at a shared time with a
# velocity discontinuity and no marker between them — is not found. Nothing
# writes a history in that shape today.

# Below this, a velocity jump at a segment boundary is not a maneuver (km/s).
#
# 1e-6 km/s is 1 mm/s. It sits well above the disagreement you get between the last state of one
# propagation arc and the first of the next when both were integrated to a tight tolerance, and
# well below any burn anyone models — a cold-gas attitude pulse is still tens of mm/s. Anything
# between those two is not a real case, which is why a single fixed number does the job here.
const _MANEUVER_DV_TOL_KMS = 1e-6

"""
    _extract(sat::Spacecraft; colors) -> (epoch_iso, segments, maneuvers)

Walk `sat.history.segments` and return everything the CZML writer needs.
`epoch_iso` is the ISO 8601 UTC string of the earliest sample.
`segments` are colored propagation arcs; `maneuvers` are ΔV points.
"""
function _extract(sat::Spacecraft; colors=nothing, linewidth::Real=1.5)
    hist = sat.history
    (hist === nothing || isempty(hist)) && throw(ArgumentError(
        "orbitview: spacecraft \"$(sat.name)\" must have a propagated history to draw; it has " *
        "none. Call propagate! before orbitview."))

    hsegs = hist.segments

    # Trajectory epoch: the earliest sample across all segments with data.
    t_ref = nothing
    for seg in hsegs
        isempty(seg.times) && continue
        candidate = seg.times[1]
        if t_ref === nothing || _time_lt(candidate, t_ref)
            t_ref = candidate
        end
    end
    t_ref === nothing && throw(ArgumentError(
        "orbitview: spacecraft \"$(sat.name)\" has a history, but every segment in it is empty, " *
        "so there is no epoch to start the clock from."))

    epoch_iso = _iso_from_time(t_ref)

    # Collect segments (≥2 samples) with colors; note where each came from
    # so we can find neighbors during maneuver detection.
    segments = _Segment[]
    prop_seg_indices = Int[]              # into hsegs, aligned with `segments`
    color_counter    = 0
    for (idx, seg) in enumerate(hsegs)
        length(seg.times) < 2 && continue
        color_counter += 1
        color = _pick_color(colors, color_counter)
        times_s = Float64[(_time_sub_seconds(t, t_ref)) for t in seg.times]
        positions_km = NTuple{3, Float64}[
            (Float64(s.position[1]), Float64(s.position[2]), Float64(s.position[3]))
            for s in seg.states]
        name = isempty(seg.name) ? "segment $color_counter" : seg.name
        push!(segments, _Segment(times_s, positions_km, name, color, Float64(linewidth)))
        push!(prop_seg_indices, idx)
    end

    # Impulsive maneuvers. A 1-sample segment is a marker rather than an arc, and
    # the burn it stands for is the velocity jump across it: the last velocity of
    # the arc before, against the first velocity of the arc after. Neither is the
    # marker's own velocity, which is why the neighbours have to be found.
    maneuvers = _Maneuver[]
    for (idx, seg) in enumerate(hsegs)
        length(seg.times) == 1 || continue

        prev_prop = _nearest_prop(hsegs, idx, -1)
        next_prop = _nearest_prop(hsegs, idx, +1)

        # Nothing before it means nothing to difference against, so the size of
        # the burn is unknown. Mark it at zero rather than dropping it.
        v_pre  = prev_prop === nothing ? nothing : _vec3(hsegs[prev_prop].states[end].velocity)
        v_post = next_prop === nothing ? _vec3(seg.states[1].velocity) :
                                         _vec3(hsegs[next_prop].states[1].velocity)

        dv = v_pre === nothing ? (0.0, 0.0, 0.0) :
             (v_post[1] - v_pre[1], v_post[2] - v_pre[2], v_post[3] - v_pre[3])

        # A boundary that happens to share a time without a real velocity change
        # is not a maneuver.
        sqrt(dv[1]^2 + dv[2]^2 + dv[3]^2) < _MANEUVER_DV_TOL_KMS && continue

        name = isempty(seg.name) ? "maneuver $(length(maneuvers) + 1)" : seg.name
        push!(maneuvers, _Maneuver(name, _time_sub_seconds(seg.times[1], t_ref), dv))
    end

    return epoch_iso, segments, maneuvers
end

# ————————————————————————————————————————————————————————————————
# Small helpers — kept private, no need to leak on user
# ————————————————————————————————————————————————————————————————

"A velocity as a plain 3-tuple, so differencing two of them needs no array package."
_vec3(v) = (Float64(v[1]), Float64(v[2]), Float64(v[3]))

# Find the nearest multi-sample segment to `from_idx` in the given `dir` (±1).
function _nearest_prop(hsegs, from_idx::Integer, dir::Integer)
    j = from_idx + dir
    while 1 <= j <= length(hsegs)
        length(hsegs[j].times) >= 2 && return j
        j += dir
    end
    return nothing
end

# Time comparison / subtraction, in days, after putting both times in TT.
#
# A history can hold segments recorded in different time scales: a spacecraft built in TAI and
# then maneuvered and propagated by a solver carries TAI and TT segments side by side, and
# AstroEpochs refuses `Time - Time` across scales or formats. Comparing the Julian dates in one
# scale is what the view needs, since it only ever orders samples and offsets them from the first.
function _days_between(a::Time, b::Time)
    a_tt = a.scale === :tt ? a : a.tt
    b_tt = b.scale === :tt ? b : b.tt
    return (a_tt.jd1 - b_tt.jd1) + (a_tt.jd2 - b_tt.jd2)
end

_time_lt(a::Time, b::Time) = _days_between(a, b) < 0
_time_sub_seconds(a::Time, b::Time) = _days_between(a, b) * 86400.0

# Convert an AstroEpochs Time to an ISO 8601 UTC string ending in "Z".
# `t.utc` returns the equivalent Time in UTC scale; `.isot` gives the ISO string.
function _iso_from_time(t::Time)
    utc_time = t.scale === :utc ? t : getproperty(t, :utc)
    iso = String(getproperty(utc_time, :isot))
    return endswith(iso, 'Z') || endswith(iso, 'z') ? iso : iso * "Z"
end
