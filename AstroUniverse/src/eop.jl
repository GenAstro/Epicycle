# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

# =============================================================================
# Earth Orientation Parameters (EOP)
#
# Single global service for the active EOP table used across the Epicycle
# stack (frame transforms, force-model adapters, etc.). On a cold cache,
# `eop()` triggers `SatelliteToolboxTransformations.fetch_iers_eop`, which
# owns its own on-disk cache; we just hold the parsed table.
#
# Theory selection is limited to what `SatelliteToolboxTransformations`
# supports: FK5 / IAU-1980 and IERS Conventions 2010 / IAU-2000A. Default is
# `IAU2000A()`. Legacy callers (e.g., the GMAT FK5-mode regression suite)
# opt in to `IAU1980()` per-test.
# =============================================================================

using SatelliteToolboxTransformations:
    EopIau1980, EopIau2000A, fetch_iers_eop, read_iers_eop

# --- Theory selection -----------------------------------------------------

"""Abstract supertype for an EOP/precession-nutation theory choice."""
abstract type AbstractEopTheory end

"""FK5: IAU-1976 precession + IAU-1980 nutation. GMAT's legacy default."""
struct IAU1980  <: AbstractEopTheory end

"""IERS Conventions 2010, IAU-2000A precession-nutation."""
struct IAU2000A <: AbstractEopTheory end

const _EopTable     = Union{Nothing, EopIau1980, EopIau2000A}
const _eop_theory   = Ref{AbstractEopTheory}(IAU2000A())   # modern default
const _eop_table    = Ref{_EopTable}(nothing)
const _eop_lock     = ReentrantLock()

# Mapping to the `Val{:IAU1980}` / `Val{:IAU2000A}` symbols consumed by
# SatelliteToolboxTransformations.
_stb_theory_symbol(::IAU1980)  = Val(:IAU1980)
_stb_theory_symbol(::IAU2000A) = Val(:IAU2000A)

# --- Public API -----------------------------------------------------------

export AbstractEopTheory, IAU1980, IAU2000A
export eop, eop_theory, set_eop_theory!, eop_load, eop_refresh!, set_eop!

"""
    eop_theory() -> AbstractEopTheory

Active EOP/precession-nutation theory for the current session. Default
`IAU2000A()` (IERS Conventions 2010).
"""
eop_theory() = _eop_theory[]

"""
    set_eop_theory!(theory::AbstractEopTheory) -> AbstractEopTheory

Switch the active EOP theory. Clears the cached table; the next call to
[`eop`](@ref) reloads (or downloads) the matching IERS file.
"""
function set_eop_theory!(theory::AbstractEopTheory)
    lock(_eop_lock) do
        _eop_theory[] = theory
        _eop_table[]  = nothing
    end
    return theory
end

"""
    eop() -> Union{EopIau1980, EopIau2000A}

Return the active EOP table. On first access this calls
`SatelliteToolboxTransformations.fetch_iers_eop`, which loads from its own
on-disk cache (downloading from IERS only if the cache is missing or stale).
Thread-safe.
"""
function eop()
    _eop_table[] === nothing && _load_active_eop!()
    return _eop_table[]
end

"""
    set_eop!(table) -> table

Install a pre-built EOP table (`EopIau1980` or `EopIau2000A`) as the active
table. Bypasses disk and network entirely. Useful for tests and Monte-Carlo
where the EOP series must not vary between runs.
"""
function set_eop!(table::Union{EopIau1980, EopIau2000A})
    lock(_eop_lock) do
        _eop_table[] = table
    end
    return table
end

"""
    eop_load(path::AbstractString; theory = eop_theory()) -> table

Parse an IERS EOP file from `path` using the given theory and install it as
the active table. Does not touch the network.
"""
function eop_load(path::AbstractString;
                  theory::AbstractEopTheory = eop_theory())
    isfile(path) || throw(ArgumentError("EOP file not found: $path"))
    table = read_iers_eop(path, _stb_theory_symbol(theory))
    return set_eop!(table)
end

"""
    eop_refresh!(; theory = eop_theory()) -> table

Force-download the latest IERS EOP series (bypassing STB's freshness check)
and install it as the active table.
"""
function eop_refresh!(; theory::AbstractEopTheory = eop_theory())
    lock(_eop_lock) do
        _eop_table[] = fetch_iers_eop(_stb_theory_symbol(theory);
                                       force_download = true)
    end
    return _eop_table[]
end

# --- Internals ------------------------------------------------------------

function _load_active_eop!()
    lock(_eop_lock) do
        theory = _eop_theory[]
        _eop_table[] = fetch_iers_eop(_stb_theory_symbol(theory))
    end
    return _eop_table[]
end
