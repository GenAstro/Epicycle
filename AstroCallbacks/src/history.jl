# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# =============================================================================
# `history` — a quantity over every recorded sample.
#
# `report` used to live here too. It moved to `EpicycleIO`, which is where getting results
# out of Epicycle belongs; this package is about quantities.
#
# The subject does not vary during a walk; the *sample* does. Every recorded
# state is a snapshot of the same spacecraft, and the walk evaluates the same
# quantity against each in turn. That is why a `Calc` carries its subject and
# why it must keep its parts: `reapply` swaps the subject for a sample and
# leaves the dependencies alone.
#
# Several columns come from one walk. Each `Calc` carries its own coordinate
# system, so columns in *different* frames still share the traversal — only the
# per-sample conversion differs, and that is the irreducible cost. Frame
# conversion is ~111 µs into ITRF against ~3.9 µs same-frame, so asking column
# by column would repeat the dominant cost.
#
# Calcs naming different subjects are walked separately, because they are
# different recordings with their own sample times. Reporting altitude for two
# spacecraft is two series, not one.
# =============================================================================

"""
    history(calcs::Calc...) -> values or (values...)

Evaluate each `Calc` over every recorded sample of its subject.

Returns one vector per `Calc`, in the order given. A single `Calc` returns a
bare vector. Columns naming the same subject share one walk of its history.

Time is a quantity like any other: ask for it with `Calc(epoch, sat)`.

# Examples
```julia
t, r, v = history(Calc(epoch,           sat),
                  Calc(position_vector, sat, EarthMJ2000Eq),
                  Calc(velocity_vector, sat, EarthMJ2000Eq))
```

Columns in different frames still cost one traversal:

```julia
t, r_eq, r_ec = history(Calc(epoch,           sat),
                        Calc(position_vector, sat, EarthMJ2000Eq),
                        Calc(position_vector, sat, EarthMJ2000Ec))
```

# Notes
Reads recorded segments. A spacecraft records by default; if
`sat.history.record_segments` was turned off there is nothing to walk and this
says so rather than returning empty columns.
"""
function history(calcs::Calc...)
    isempty(calcs) && throw(ArgumentError("`history` needs at least one Calc."))

    columns = Vector{Any}(undef, length(calcs))

    # One walk per distinct subject. `===` rather than `==`: two spacecraft
    # with identical states are still two recordings.
    subjects = Any[]
    for c in calcs
        s = _subject_of(c)
        any(x -> x === s, subjects) || push!(subjects, s)
    end

    for s in subjects
        which = [i for i in eachindex(calcs) if _subject_of(calcs[i]) === s]
        walked = _walk(s, [calcs[i] for i in which])
        for (k, i) in enumerate(which)
            columns[i] = walked[k]
        end
    end

    return length(columns) == 1 ? columns[1] : Tuple(columns)
end

"""The subject of a `Calc`, stored as its first argument."""
_subject_of(c::Calc) = first(c.args)

"""
Evaluate every `Calc` against each recorded sample of `subject`, in one pass.

A sample is rebuilt as a `Coordinate` so it carries its own frame. Segments may
have been recorded in different frames, and because each sample states its own,
a column asking for `EarthMJ2000Ec` converts correctly from whichever that was
with no bookkeeping here.
"""
function _walk(subject, calcs::Vector{<:Calc})
    h = _history_of(subject)
    isempty(h.segments) && throw(ArgumentError(
        "$(_name_of(subject)) has no recorded segments to walk. Propagate first, " *
        "or check that `history.record_segments` is on."))

    columns = [Any[] for _ in calcs]
    for seg in h.segments, i in eachindex(seg.times)
        sample = Coordinate(seg.states[i], seg.coordinate_system, seg.times[i])
        for (k, c) in enumerate(calcs)
            push!(columns[k], reapply(c, sample))
        end
    end
    return [identity.(col) for col in columns]     # narrow Any[] to a real type
end

_history_of(sc) = sc.history
_name_of(sc) = hasproperty(sc, :name) && !isempty(sc.name) ? sc.name : string(typeof(sc))
