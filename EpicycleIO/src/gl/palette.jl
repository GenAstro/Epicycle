# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# CZML wants colours as RGBA integer tuples; the plotting side wants CSS strings.
# Same eight colours either way, so a trajectory and the plots of it match.

const DEFAULT_PALETTE = NTuple{4, Int}[
    (255,  80,  80, 255), (255, 165,   0, 255), (255, 220,   0, 255), ( 55, 255,  55, 255),
    ( 55, 200, 255, 255), (140, 100, 255, 255), (255, 120, 200, 255), (  0, 220, 255, 255),
]

_pick_color(colors, i::Integer) =
    (colors === nothing ? DEFAULT_PALETTE : _normalize_palette(colors))[
        mod1(i, length(colors === nothing ? DEFAULT_PALETTE : _normalize_palette(colors)))]

_normalize_palette(c::NTuple{4, <:Integer}) = NTuple{4, Int}[Int.(c)]
_normalize_palette(v::AbstractVector) = NTuple{4, Int}[NTuple{4, Int}(Int.(c)) for c in v]
