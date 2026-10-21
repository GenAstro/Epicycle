# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Frame, centre and time-system names in ephemeris files, and the Epicycle objects they mean.
#
# Sources: CCSDS 502.0-B-3 §3.2.3.2 and §3.2.3.3, and the SANA registries the standard makes
# normative (Annex B): Celestial Body Reference Frames, Orbit Centers and Time Systems. The two
# disagree on some spellings (the standard's TOD and TEME are TOD_EARTH and TEMEOFDATE in the
# registry), so a reader accepts both and a writer uses the spelling found in both where one
# exists, and the standard's otherwise.

# --- frames -------------------------------------------------------------------------------

# Writing: Epicycle axes → the name written. Axes not listed cannot be written.
const _FRAME_NAME = Dict{DataType, String}(
    ICRF => "ICRF", MJ2000Eq => "EME2000", GCRF => "GCRF", ITRF => "ITRF", TEME => "TEME",
    TODEq => "TOD", MODEq => "MOD_EARTH", MJ2000Ec => "J2000_ECLIPTIC", CIRS => "CIRS",
    TIRS => "TIRS", PEF => "GTOD", MoonME => "MOON_ME", MoonPA => "MOON_PA440",
)

# Reading: name → (axes, note). A non-empty note says what reading it this way assumes.
const _FRAME_AXES = Dict{String, Tuple{Any, String}}(
    "EME2000" => (MJ2000Eq(), ""), "J2000" => (MJ2000Eq(), ""),
    "ICRF" => (ICRF(), ""), "ICRF1" => (ICRF(), ""), "ICRF2" => (ICRF(), ""), "ICRF3" => (ICRF(), ""),
    "GCRF" => (GCRF(), ""), "GCRF1" => (GCRF(), ""), "GCRF2" => (GCRF(), ""), "GCRF3" => (GCRF(), ""),
    "ITRF" => (ITRF(), ""),
    "FIXED_EARTH" => (ITRF(), "FIXED_EARTH read as ITRF; the registry allows IAU 1976/1980 or 2000A realizations too"),
    "TEME" => (TEME(), ""), "TEMEOFDATE" => (TEME(), ""),
    "TOD" => (TODEq(), ""), "TOD_EARTH" => (TODEq(), ""),
    "MOD_EARTH" => (MODEq(), ""),
    "J2000_ECLIPTIC" => (MJ2000Ec(), ""),
    "CIRS" => (CIRS(), ""), "TIRS" => (TIRS(), ""),
    "GTOD" => (PEF(), ""), "GRC" => (PEF(), ""), "TDR" => (PEF(), ""),
    "MOON_ME" => (MoonME(), ""), "MOON_PA440" => (MoonPA(), ""),
)

# Registered names Epicycle has no axes for. Reading one is an error that says so, rather than
# "unknown frame".
const _FRAMES_REGISTERED_UNSUPPORTED = (
    "ALIGN_CB", "ALIGN_EARTH", "B1950", "EFG", "INERTIAL_CB", "J2000A", "MCI", "MOD_CB",
    "MOD_MOON", "MOE_CB", "MOE_EARTH", "MOON_MEIAUE", "TEMEOFEPOCH", "TOD_CB", "TOD_MOON",
    "TOE_CB", "TOE_EARTH", "TOE_MOON", "TRUE_ECLIPTIC", "UVW_GO_INERTIAL", "WGS84",
)

"""
    _frame_axes(name, center) -> (axes, note, problem)

The Epicycle axes a file's REF_FRAME names. `axes` is `nothing` when the name cannot be read,
with `problem` saying why; `note` is non-empty when the reading assumes something.
"""
function _frame_axes(name::AbstractString, center::AbstractString)
    n = uppercase(strip(name))
    haskey(_FRAME_AXES, n) && return (_FRAME_AXES[n]..., "")
    if occursin(r"^(ITRF-?\d{2,4}|DTRF\d{4})$", n)
        return ITRF(), "$n read as ITRF; realizations differ at the centimetre level", ""
    end
    if occursin(r"^MOON_PA\d{3}$", n)
        return nothing, "", "$n is the Moon principal-axes frame of a different planetary " *
            "ephemeris; Epicycle's is MOON_PA440 (DE440). To accept the difference, override " *
            "it with ref_frame = \"MOON_PA440\""
    end
    if n == "FIXED_CB"
        center in ("", "EARTH") &&
            return ITRF(), "FIXED_CB with an Earth centre read as ITRF", ""
        return CelestialBodyFixed(), "", ""
    end
    n in _FRAMES_REGISTERED_UNSUPPORTED &&
        return nothing, "", "$n is a registered frame that Epicycle does not provide"
    known = vcat(collect(keys(_FRAME_AXES)), "FIXED_CB", "ITRF2000", collect(_FRAMES_REGISTERED_UNSUPPORTED))
    return nothing, "", "$(repr(String(strip(name)))) is not a known frame.$(_did_you_mean(n, known))"
end

"""
    _frame_name(axes, center) -> String

The REF_FRAME name a writer gives `axes`. Throws `ArgumentError` for axes no file format names.
"""
function _frame_name(axes, center_name::AbstractString)
    if axes isa CelestialBodyFixed
        center_name == "EARTH" && return "ITRF"
        return "FIXED_CB"
    end
    name = get(_FRAME_NAME, typeof(axes), nothing)
    name === nothing && throw(ArgumentError(
        "$(nameof(typeof(axes))) axes have no name in ephemeris files; supported axes are " *
        join(sort!(string.(nameof.(collect(keys(_FRAME_NAME))))), ", ") * ", and a body's CelestialBodyFixed"))
    return name
end

# --- centres -------------------------------------------------------------------------------

const _BODIES = Dict("SUN" => sun, "MERCURY" => mercury, "VENUS" => venus, "EARTH" => earth,
                     "MOON" => moon, "MARS" => mars, "JUPITER" => jupiter, "SATURN" => saturn,
                     "URANUS" => uranus, "NEPTUNE" => neptune, "PLUTO" => pluto)
const _BARYCENTERS = ("SOLAR SYSTEM BARYCENTER", "MERCURY BARYCENTER", "VENUS BARYCENTER",
                      "EARTH BARYCENTER", "MARS BARYCENTER", "JUPITER BARYCENTER",
                      "SATURN BARYCENTER", "URANUS BARYCENTER", "NEPTUNE BARYCENTER",
                      "PLUTO BARYCENTER")

"""
    _center_body(name) -> (body, problem)

The Epicycle body a file's CENTER_NAME names, or `nothing` and why not.
"""
function _center_body(name::AbstractString)
    n = uppercase(strip(name))
    haskey(_BODIES, n) && return _BODIES[n], ""
    n in _BARYCENTERS && return nothing,
        "$n is a barycentre; Epicycle's origins are body centres. If the file's states are " *
        "relative to the body itself, override with center = \"$(replace(n, " BARYCENTER" => ""))\""
    occursin(r"L[1-5]$", n) && return nothing, "$n is a libration point, which is not an Epicycle origin"
    return nothing, "$(repr(String(strip(name)))) is not a centre Epicycle provides.$(_did_you_mean(n, keys(_BODIES)))"
end

"""
    _center_name(body) -> String

The CENTER_NAME a writer gives `body`. Throws `ArgumentError` for a body with no registered name.
"""
function _center_name(body)
    for (k, v) in _BODIES
        v === body && return k
    end
    n = uppercase(body.name)
    haskey(_BODIES, n) && return n
    throw(ArgumentError("no CENTER_NAME for body $(body.name); supported centres are " *
                        join(sort!(collect(keys(_BODIES))), ", ")))
end

# --- time systems --------------------------------------------------------------------------

"""
    _time_system(name) -> (system, note, problem)

The supported time system a file's TIME_SYSTEM names: one of the keys of `_TIME_SYSTEMS`.
"""
function _time_system(name::AbstractString)
    n = uppercase(strip(name))
    haskey(_TIME_SYSTEMS, n) && return n, "", ""
    n == "ET" && return "TDB", "ET read as TDB, as SPICE defines it", ""
    n in (_UNSUPPORTED_TIME_SYSTEMS..., "BEIDOU", "GALILEO", "GLONASS", "NAVIC") &&
        return "", "", "$n is a registered time system that Epicycle does not convert"
    return "", "", "$(repr(String(strip(name)))) is not a known time system.$(_did_you_mean(n, keys(_TIME_SYSTEMS)))"
end

"""
    _time_system_name(scale) -> String

The TIME_SYSTEM name for an AstroEpochs scale tag (`UTC()`, `TDB()`, …) or the string `"GPS"`.
"""
_time_system_name(s::AbstractString) = (n = uppercase(s); haskey(_TIME_SYSTEMS, n) ? n :
    throw(ArgumentError("time system $(repr(s)) cannot be written; supported: " *
                        join(sort!(collect(keys(_TIME_SYSTEMS))), ", "))))
_time_system_name(s) = _time_system_name(uppercase(String(nameof(typeof(s)))))
