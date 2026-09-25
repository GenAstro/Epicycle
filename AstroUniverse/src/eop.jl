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
# other. On a cold cache `fetch_iers_eop` downloads, so a reproducible run
# should `set_eop!` both tables up front rather than rely on the network.
# =============================================================================

using SatelliteToolboxTransformations:
    EopIau1980, EopIau2000A, fetch_iers_eop, read_iers_eop

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
abstract type AbstractFrameTheory end

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
be live simultaneously and asking for one never fetches the other. On a cold
cache this calls `SatelliteToolboxTransformations.fetch_iers_eop`, which owns
its own on-disk cache and downloads only when that is missing or stale.

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
    return set_eop!(read_iers_eop(path, _stb_theory_symbol(theory)))
end

"""
    eop_refresh!(; theory = frame_theory()) -> table

Force-download the latest IERS EOP series for `theory`, bypassing the
freshness check, and install it in that theory's slot.

# Returns
The downloaded and installed EOP table.

# Example
```julia
table = eop_refresh!(; theory = IAU2006())
```
"""
function eop_refresh!(; theory::AbstractFrameTheory = frame_theory())
    return set_eop!(fetch_iers_eop(_stb_theory_symbol(theory);
                                   force_download = true))
end

# --- Internals --------------------------------------------------------------

function _load_eop!(theory::AbstractFrameTheory)
    lock(_eop_lock) do
        slot = _eop_slot(theory)
        # Re-check under the lock: another task may have filled it while we
        # waited.
        slot[] === nothing && (slot[] = fetch_iers_eop(_stb_theory_symbol(theory)))
        return slot[]
    end
end
