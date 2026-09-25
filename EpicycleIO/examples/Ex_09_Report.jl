# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Writing results to a file, and the one rule that changed when `report` moved here.
#
# `report` no longer knows what an epoch is. It used to special-case `Time` and
# silently write a bare Julian date — which says nothing about whether it is TT or
# UTC, and made a text writer carry a domain decision it had no business making.
#
# You render the epoch. `t.jd`, `t.mjd` and `t.isot` are all there, and whichever
# you pick, the file says what you meant.
#
# Review question: is requiring the conversion the right call, or too blunt? The
# alternative is a numeric epoch quantity so you ask `history` for a number in the
# first place — which is the open question in Quantities_Core §10.

using EpicycleIO
include(joinpath(@__DIR__, "synthetic.jl"))

t, alt, r = fake_orbit()

# Standing in for `history(Calc(epoch, sat), ...)`, whose first column would be
# `Time` values rather than these numbers.
jd = 2_459_580.5 .+ t ./ 24

# ── Columns become headings ───────────────────────────────────────────────────
# The keywords are your words, not the quantity's, because someone who did not
# write the run has to read the file. A vector column expands into components, so
# `position` becomes position_1, position_2, position_3.

# Written to a temporary directory rather than the working directory. A relative path
# lands wherever the example happened to be run from, which for a script people run
# out of curiosity means a stray file in a repository.
const REPORT = joinpath(tempdir(), "flight.txt")

report(REPORT;
       jd_utc   = jd,
       altitude = alt,
       position = r)

# ── What is no longer allowed ─────────────────────────────────────────────────
# Passing the `Time` values straight through used to work by accident. It now
# fails, and the error names the fix. An epoch does not render on one line, and a
# value that spans lines would quietly turn one row into four.
#
#   t_epochs = history(Calc(epoch, sat))
#   report("flight.txt"; time = t_epochs, altitude = alt)      # errors
#   report("flight.txt"; time = [x.isot for x in t_epochs], altitude = alt)   # say what you mean

# ── The same data, reported and plotted ───────────────────────────────────────
# The two take the same arrays from the same place. They are not alike in syntax,
# because a file has column headings and a plot has series, and those are named
# differently for good reason.

xyplot("Position", t, r; name = ["x", "y", "z"])
panel!("Position"; xaxis_title = "hours", yaxis_title = "km")
