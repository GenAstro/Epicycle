# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# =============================================================================
# Frame theory and Earth Orientation Parameters (EOP)
#
# The frame theory — which precession-nutation model the Earth chain uses —
# is a property of the universe, and it lives here. The EOP tables live with
# it, because the theory is what selects them.
#
# **Both tables are held at once.** They are not interchangeable: `EopIau1980`
# carries the nutation corrections `δΔψ`/`δΔϵ`, `EopIau2000A` carries the CIP
# offsets `δx`/`δy`. An FK5 edge needs fields that do not exist on a 2000A
# table — a missing field, not a precision difference. And which table an edge
# needs is fixed by the *edge*, not by this setting: `TODEq` and `PEF` are FK5
# nodes, `CIRS` and `TIRS` are IAU-2006 nodes, and a caller may name either at
# any time. A single active table would make every edge of the other theory
# unserviceable.
#
# The theory setting therefore does one thing: choose the default chain when
# the endpoints do not name a branch. It does not decide which data exists.
#
# Each slot fills lazily and independently; asking for one never fetches the
# other. The tables come from a cache of the IERS files, refreshed from the
# network when it is a week old; a reproducible run should `set_eop!` both
# tables up front rather than rely on the network.
#
# **No single server is relied on.** The IERS data centre has refused every
# request for days at a time, so each series has an ordered list of sources:
# the IERS CSV, then the US Naval Observatory's fixed-width file of the same
# series, which USNO produces as the IERS Rapid Service/Prediction Centre. A
# download replaces the cached copy only once it has been converted to the CSV
# layout, read back as a table, and found to reach today, and then in one
# rename that leaves the old copy in place if it fails. When no source supplies
# usable data the cached copy is used, with a warning saying how old it is; only
# with no copy at all is it an error.
# =============================================================================

using SatelliteToolboxTransformations: EopIau1980, EopIau2000A, read_iers_eop

# --- Frame theory -----------------------------------------------------------

"""
    AbstractFrameTheory

Supertype for the Earth precession-nutation theory in force.

The theory names the *model*, not the data product: which EOP series is
required follows from it and is not something a caller selects.

See also: [`FK5`](@ref), [`IAU2006`](@ref), [`frame_theory`](@ref).

# Example
```jldoctest
julia> using AstroUniverse

julia> FK5() isa AbstractFrameTheory
true
```
"""
abstract type AbstractFrameTheory <: AbstractOrientationModel end

# A frame theory is Earth's orientation model, but its rotation is the ITRF
# chain, which reads UT1 and TT and lives in AstroFrames. A TDB Julian date alone
# cannot give it, so this points there rather than approximating.
function body_axes_rotation(th::AbstractFrameTheory, naifid::Integer, jd_tdb::Real)
    throw(ArgumentError(
        "$(th) is Earth's orientation, the ITRF chain, which needs UT1 and TT as well " *
        "as TDB. Evaluate it with AstroFrames: `body_fixed_rotation($(th), 399, epoch)` " *
        "or `axes_rotation(ICRF(), ITRF(), epoch)` with an `AstroEpochs.Time` epoch."))
end

orientation_parameters(::AbstractFrameTheory) = NamedTuple()
_orients(::AbstractFrameTheory, naifid) = naifid == 399

"""
    FK5()

IAU-76 precession with IAU-80 nutation, the classical system, and the one heritage
tools and GMAT use. Uses the IERS IAU-1980 EOP series.

# Example
```jldoctest
julia> using AstroUniverse

julia> FK5() isa AbstractFrameTheory
true
```
"""
struct FK5 <: AbstractFrameTheory end

"""
    IAU2006()

IAU-2006/2010 precession-nutation, CIO-based, which is the modern standard. Uses
the IERS IAU-2000A EOP series.

# Example
```jldoctest
julia> using AstroUniverse

julia> IAU2006() isa AbstractFrameTheory
true
```
"""
struct IAU2006 <: AbstractFrameTheory end

const _eop_1980     = Ref{Union{Nothing, EopIau1980}}(nothing)
const _eop_2000a    = Ref{Union{Nothing, EopIau2000A}}(nothing)
const _frame_theory = Ref{AbstractFrameTheory}(IAU2006())   # modern default
const _eop_lock     = ReentrantLock()

# The `Val{:IAU1980}` / `Val{:IAU2000A}` names are SatelliteToolboxTransformations'
# argument spelling for the EOP *series*. They stay here, at the boundary, and
# are not part of our vocabulary.
_stb_theory_symbol(::FK5)     = Val(:IAU1980)
_stb_theory_symbol(::IAU2006) = Val(:IAU2000A)

_eop_slot(::FK5)     = _eop_1980
_eop_slot(::IAU2006) = _eop_2000a

# --- Public API -------------------------------------------------------------

export AbstractFrameTheory, FK5, IAU2006
export eop, frame_theory, set_frame_theory!, eop_load, eop_refresh!, set_eop!

"""
    frame_theory() -> AbstractFrameTheory

Earth precession-nutation theory in force for this session. Defaults to
[`IAU2006`](@ref).

This selects the default chain when a transform's endpoints do not name a
branch. It does not restrict which EOP data is available, since both series can
be loaded at once (see [`eop`](@ref)).

# Returns
The active [`AbstractFrameTheory`](@ref).

# Example
```jldoctest
julia> using AstroUniverse

julia> frame_theory() isa AbstractFrameTheory
true
```
"""
frame_theory() = _frame_theory[]

"""
    set_frame_theory!(theory::AbstractFrameTheory) -> theory

Set the Earth precession-nutation theory for this session.

Changing it changes computed results, because the two theories are different
physical models. A result that has to be reproducible is recorded with the theory
that produced it.

The theory is Earth's orientation model, so `set_orientation!(earth, theory)`
is the same call. It governs Earth frames, the default EOP series, and Earth
gravity built after it is set.

# Returns
The installed `theory`.

# Example
```jldoctest
julia> using AstroUniverse

julia> set_frame_theory!(IAU2006())
IAU2006()
```
"""
function set_frame_theory!(theory::AbstractFrameTheory)
    lock(_eop_lock) do
        _frame_theory[] = theory
    end
    return theory
end

"""
    eop(theory::AbstractFrameTheory) -> EopIau1980 | EopIau2000A
    eop() -> EopIau1980 | EopIau2000A

EOP table for `theory`, or for the active theory if none is given.

Each theory's table occupies its own slot and fills on first use, so both can
be live simultaneously and asking for one never fetches the other. The table is
read from the cached IERS file, which is downloaded when it is missing or a week
old. Each series has more than one source, IERS first and then the US Naval
Observatory; when none can be reached a cached copy is used, with a warning
giving its age.

Thread-safe.

# Returns
An `EopIau1980` table for `FK5()` or an `EopIau2000A` table for `IAU2006()`.

# Example
```julia
table = eop(IAU2006())
```
"""
function eop(theory::AbstractFrameTheory)
    slot = _eop_slot(theory)
    slot[] === nothing && _load_eop!(theory)
    return slot[]
end

eop() = eop(frame_theory())

"""
    set_eop!(table) -> table

Install a pre-built EOP table (`EopIau1980` or `EopIau2000A`).

The table's type selects its slot, so installing one never disturbs the other.
Bypasses disk and network entirely.

A reproducible run is one whose numbers do not move between sessions: a
Monte-Carlo, a regression, a delivered analysis. Install both tables up front for
one, which is the only way to guarantee that no evaluation reaches the network.

# Returns
The installed table.

# Example
```julia
set_eop!(table)
```
"""
function set_eop!(table::EopIau1980)
    lock(_eop_lock) do
        _eop_1980[] = table
    end
    return table
end

function set_eop!(table::EopIau2000A)
    lock(_eop_lock) do
        _eop_2000a[] = table
    end
    return table
end

"""
    eop_load(path::AbstractString; theory = frame_theory()) -> table

Parse an IERS EOP file from `path` for `theory` and install it in that
theory's slot. Does not touch the network.

# Arguments
- `path`: an IERS finals file, either the CSV the IERS data centre publishes
  (`finals2000A.all.csv`, `finals.all.csv`) or the fixed-width file the US Naval
  Observatory publishes (`finals2000A.all`, `finals.all`). The format is read
  from the file.
- `theory`: the theory the file is for. A `finals2000A` file carries the IAU
  2000A offsets [`IAU2006`](@ref) needs, a `finals` file the IAU 1980
  corrections [`FK5`](@ref) needs.

# Returns
The parsed and installed EOP table.

# Example
```julia
table = eop_load("finals2000A.all"; theory = IAU2006())
```
"""
function eop_load(path::AbstractString;
                  theory::AbstractFrameTheory = frame_theory())
    isfile(path) || throw(ArgumentError(
        "EOP file not found: $(path)"))
    _is_iers_csv(path) && return set_eop!(read_iers_eop(path, _stb_theory_symbol(theory)))
    csv = tempname() * ".csv"
    try
        _finals_to_csv(path, csv, theory)
        return set_eop!(read_iers_eop(csv, _stb_theory_symbol(theory)))
    finally
        rm(csv; force = true)
    end
end

"""
    eop_refresh!(; theory = frame_theory()) -> table

Download the latest EOP series for `theory` now, bypassing the freshness check,
and install it in that theory's slot.

The sources are tried in order, IERS first and then the US Naval Observatory.
Unlike the refresh `eop` makes when the cache is a week old, this throws when
no source supplies usable data, because a caller who asks for a refresh wants to
know it did not happen. The cached copy is left as it was.

# Returns
The downloaded and installed EOP table.

# Example
```julia
table = eop_refresh!(; theory = IAU2006())
```
"""
function eop_refresh!(; theory::AbstractFrameTheory = frame_theory())
    # Under the same lock as a first load, so the two never write the cache at once.
    lock(_eop_lock) do
        set_eop!(read_iers_eop(_eop_file(theory; refresh = true), _stb_theory_symbol(theory)))
    end
end

# --- Internals --------------------------------------------------------------

function _load_eop!(theory::AbstractFrameTheory)
    lock(_eop_lock) do
        slot = _eop_slot(theory)
        # Re-check under the lock: another task may have filled it while we
        # waited.
        slot[] === nothing &&
            (slot[] = read_iers_eop(_eop_file(theory), _stb_theory_symbol(theory)))
        return slot[]
    end
end

# --- Sources and the cache ---------------------------------------------------

# One place an EOP series can be downloaded from, and the layout it comes in: the IERS CSV
# (`:csv`) or the fixed-width finals file (`:finals`).
struct _EopSource
    url::String
    format::Symbol
end

_eop_sources(::FK5) = [
    _EopSource("https://datacenter.iers.org/data/csv/finals.all.csv", :csv),
    _EopSource("https://maia.usno.navy.mil/ser7/finals.all", :finals)]
_eop_sources(::IAU2006) = [
    _EopSource("https://datacenter.iers.org/data/csv/finals2000A.all.csv", :csv),
    _EopSource("https://maia.usno.navy.mil/ser7/finals2000A.all", :finals)]

_eop_cache_name(::FK5)     = "finals.all.csv"
_eop_cache_name(::IAU2006) = "finals2000A.all.csv"

# Where SatelliteToolboxTransformations kept the same files before AstroUniverse had a cache of its
# own. Read in place, never written, as the copy of last resort.
const _STB_UUID = "6b019ec1-7a1e-4f04-96c7-a9db1ca5514d"
_legacy_eop_key(::FK5)     = "eop_iau1980"
_legacy_eop_key(::IAU2006) = "eop_iau2000A"

const _EOP_MAX_AGE = 7 * 86_400.0          # seconds before a cached file is refreshed

# Test hooks: the sources tried, the cache directory, the directory of the older cache, how old a
# cached file may get, the Modified Julian Date a download must reach (today by default), and how
# a new file replaces the cached one.
const _EOP_SOURCES    = Ref{Any}(nothing)            # theory -> sources
const _EOP_CACHE_DIR  = Ref{Union{Nothing, String}}(nothing)
const _EOP_LEGACY_DIR = Ref{Union{Nothing, String}}(nothing)
const _EOP_AGE_LIMIT  = Ref(_EOP_MAX_AGE)
const _EOP_TODAY_MJD  = Ref{Union{Nothing, Float64}}(nothing)
const _EOP_REPLACE    = Ref{Any}(nothing)

_eop_source_list(theory) = (f = _EOP_SOURCES[]; f === nothing ? _eop_sources(theory) : f(theory))
_today_mjd()             = something(_EOP_TODAY_MJD[], time() / 86_400 + 40_587.0)

# The package's Scratch space is created only when no directory is set.
function _eop_cache_dir()
    dir = _EOP_CACHE_DIR[]
    return dir === nothing ? @get_scratch!("eop") : dir
end

function _legacy_eop_file(theory)
    dir = _EOP_LEGACY_DIR[]
    dir === nothing &&
        (dir = joinpath(first(DEPOT_PATH), "scratchspaces", _STB_UUID, _legacy_eop_key(theory)))
    path = joinpath(dir, _eop_cache_name(theory))
    return isfile(path) ? path : nothing
end

_age_days(path) = round((time() - mtime(path)) / 86_400; digits = 1)

"""
    _eop_file(theory; refresh = false) -> path

The cached IERS CSV for `theory`, downloaded first when it is missing, a week old, or `refresh`
is asked for. Each source is tried in turn. When none supplies usable data, the cached copy is
used with a warning, or, with no copy of its own, the file SatelliteToolboxTransformations cached
earlier, read where it is. A `refresh` that was asked for throws instead and leaves the copy as it
was. Called under `_eop_lock`.
"""
function _eop_file(theory::AbstractFrameTheory; refresh::Bool = false)
    path = joinpath(_eop_cache_dir(), _eop_cache_name(theory))
    have = isfile(path)
    have && !refresh && time() - mtime(path) < _EOP_AGE_LIMIT[] && return path

    failures = String[]
    for source in _eop_source_list(theory)
        try
            _download_eop!(path, source, theory)
            return path
        catch err
            err isa InterruptException && rethrow()
            push!(failures, "  $(source.url)\n    $(_one_line(err))")
        end
    end

    tried = join(failures, "\n")
    if have
        refresh && error("Could not refresh the Earth orientation data for $(theory); the cached " *
                         "copy, $(_age_days(path)) days old, is unchanged. No source supplied " *
                         "usable data:\n$(tried)")
        @warn "Could not refresh the Earth orientation data for $(theory); using the cached " *
              "copy, $(_age_days(path)) days old. No source supplied usable data:\n$(tried)"
        return path
    end
    legacy = refresh ? nothing : _legacy_eop_file(theory)
    if legacy !== nothing
        @warn "No Earth orientation data for $(theory) in AstroUniverse's cache, and no source " *
              "supplied usable data; using the copy SatelliteToolboxTransformations cached " *
              "earlier, $(_age_days(legacy)) days old, at $(legacy).\n$(tried)"
        return legacy
    end
    error("No Earth orientation data for $(theory): there is no cached copy and no source " *
          "supplied usable data:\n$(tried)\nDownload one of these files by other means and load " *
          "it with `eop_load(path; theory = $(theory))`.")
end

_one_line(err) = first(split(sprint(showerror, err), '\n'))

# Replace `dst` with `src` in one step, with the operating system's rename: on Windows and POSIX it
# replaces an existing file atomically and leaves it as it was when it fails. `mv(; force = true)`
# deletes the destination first, so a failure there would lose the cache. This is the call
# `Base.Filesystem.rename` makes, without its fallback to copying.
function _replace_file(src, dst)
    err = ccall(:jl_fs_rename, Int32, (Cstring, Cstring), src, dst)
    err < 0 && Base.uv_error("rename($(repr(src)), $(repr(dst)))", err)
    return dst
end

# Download one source into the cache. Nothing replaces the cached file until the download has
# been converted to the CSV layout, read back as a table, found to hold only finite values, and
# found to reach today, so a failed, truncated, corrupt, out-of-date or wrong file leaves the
# cache as it was. Each attempt has its own
# temporary files beside the cache, so attempts in other tasks or processes never share one and
# the final rename stays on one volume.
function _download_eop!(path, source::_EopSource, theory)
    dir = mkpath(dirname(path))
    raw = tempname(dir; cleanup = false)
    csv = tempname(dir; cleanup = false)
    try
        Downloads.download(source.url, raw; timeout = 120)
        if source.format === :finals
            _finals_to_csv(raw, csv, theory)
        else
            _replace_file(raw, csv)
        end
        read_iers_eop(csv, _stb_theory_symbol(theory))
        last = _check_eop_csv(csv)
        last >= _today_mjd() || error(
            "its UT1-UTC series ends at MJD $(last), before today (MJD $(floor(Int, _today_mjd())))")
        something(_EOP_REPLACE[], _replace_file)(csv, path)
    finally
        rm(raw; force = true)
        rm(csv; force = true)
    end
    return path
end

# The columns of an IERS CSV a table is built from, by layout: polar motion and its errors,
# UT1-UTC, LOD and their errors, and the nutation corrections of both theories with their errors.
# The 33-column layout lacks the four polar-motion rate columns, so the later columns move up four.
function _eop_value_columns(ncol)
    ncol == 37 && return (ut1 = 15, all = [6:9; 15:18; 20:27])
    ncol == 33 && return (ut1 = 11, all = [6:9; 11:14; 16:23])
    error("an IERS CSV has 33 or 37 columns; this one has $(ncol)")
end

"""
    _check_eop_csv(csv) -> mjd

Check the values an IERS CSV gives a table, and return the last date with a UT1-UTC value, which
is how far the table can be used. Throws when a value that is present is not a finite number. An
empty field is a missing value, which the table reader allows only at the end of a series and
rejects elsewhere; the reader takes `Inf` as a value, so it is refused here.
"""
function _check_eop_csv(csv)
    open(csv) do io
        cols = _eop_value_columns(count(==(';'), readline(io)) + 1)
        last = -Inf
        for line in eachline(io)
            f   = split(line, ';')
            mjd = tryparse(Float64, f[1])
            for c in cols.all
                c <= length(f) && !isempty(f[c]) || continue
                v = tryparse(Float64, f[c])
                v !== nothing && isfinite(v) || error(
                    "the row for MJD $(f[1]) has $(repr(String(f[c]))) in column $(c), which is " *
                    "not a finite number")
                c == cols.ut1 && mjd !== nothing && (last = max(last, mjd))
            end
        end
        return last
    end
end

# --- The fixed-width finals file ---------------------------------------------------------------
#
# The US Naval Observatory publishes each series as a fixed-width text file, documented in
# https://maia.usno.navy.mil/ser7/readme.finals2000A, and the IERS data centre republishes the
# same values as a semicolon-separated file. The fixed-width file is rewritten here in the CSV
# layout, field for field and in the same units, so both are read by one reader.
#
# The two series differ only in the nutation columns: `finals2000A` carries the IAU 2000A offsets
# dX, dY and `finals` the IAU 1980 corrections dPsi, dEpsilon, in the same columns of the
# fixed-width file and in different columns of the CSV.

const _IERS_CSV_HEADER =
    "MJD;Year;Month;Day;Type;x_pole;sigma_x_pole;y_pole;sigma_y_pole;x_rate;sigma_x_rate;" *
    "y_rate;sigma_y_rate;Type;UT1-UTC;sigma_UT1-UTC;LOD;sigma_LOD;Type;dPsi;sigma_dPsi;" *
    "dEpsilon;sigma_dEpsilon;dX;sigma_dX;dY;sigma_dY;Type;bulB/x_pole;bulB/y_pole;Type;" *
    "bulB/UT-UTC;Type;bulB/dPsi;bulB/dEpsilon;bulB/dX;bulB/dY"

_is_iers_csv(path) = open(io -> startswith(readline(io), "MJD;"), path)

# Columns of the fixed-width file, from the readme.
_field(line, r) = strip(SubString(line, first(r), min(last(r), lastindex(line))))
_flag(c) = c == "I" ? "final" : c == "P" ? "prediction" : ""

"""
    _finals_to_csv(src, dst, theory)

Rewrite the fixed-width finals file `src` as an IERS CSV at `dst`, for the series `theory` uses.
Throws when `src` holds no line in the fixed-width layout.
"""
function _finals_to_csv(src, dst, theory::AbstractFrameTheory)
    iau1980 = theory isa FK5
    rows = 0
    open(dst, "w") do io
        println(io, _IERS_CSV_HEADER)
        for line in eachline(src)
            length(line) >= 15 || continue
            mjd = tryparse(Float64, _field(line, 8:15))
            yy  = tryparse(Int, _field(line, 1:2))
            (mjd === nothing || yy === nothing) && continue
            year = yy + (mjd <= 51_543 ? 1900 : 2000)
            f(r) = length(line) >= first(r) ? _field(line, r) : ""
            nut  = (f(98:106), f(107:115), f(117:125), f(126:134))
            empty4 = ("", "", "", "")
            dpsi, dxy = iau1980 ? (nut, empty4) : (empty4, nut)
            bnut = (f(166:175), f(176:185))
            bdpsi, bdxy = iau1980 ? (bnut, ("", "")) : (("", ""), bnut)
            fields = (string(Int(mjd)), string(year), lpad(_field(line, 3:4), 2, '0'),
                      lpad(_field(line, 5:6), 2, '0'), _flag(f(17:17)),
                      f(19:27), f(28:36), f(38:46), f(47:55), "", "", "", "",
                      _flag(f(58:58)), f(59:68), f(69:78), f(80:86), f(87:93),
                      _flag(f(96:96)), dpsi..., dxy...,
                      isempty(f(135:144)) ? "" : "final", f(135:144), f(145:154),
                      isempty(f(155:165)) ? "" : "final", f(155:165),
                      isempty(bnut[1]) ? "" : "final", bdpsi..., bdxy...)
            println(io, join(fields, ';'))
            rows += 1
        end
    end
    rows > 0 || throw(ArgumentError(
        "$(src) is not an IERS finals file: no line has the fixed-width layout."))
    return dst
end
