# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# State-transform test harness. See the rotation test matrix, an internal
# design note that is not shipped with this package.
#
# A truth vector is a Vector{AxesRotationTruthRow}. Each row carries the
# names of source and target axes, an epoch, and a hand-generated 6×6 truth
# matrix from SPICE, GMAT, Astropy, or analytic derivation. `run_truth_row`
# resolves the axes-name strings to instances, invokes `axes_rotation`,
# and Frobenius-diffs against the truth matrix.

using AstroFrames
using AstroUniverse
using StaticArrays
using LinearAlgebra: norm
using Test

struct AxesRotationTruthRow
    id::String
    source_frame::String
    target_frame::String
    epoch_jd_tdb::Float64
    M::SMatrix{6,6,Float64,36}
    tolerance::Float64
    source::Symbol            # :SPICE | :GMAT | :Astropy | :Analytic
    source_note::String
end

# Axes-name resolver. Extend as new axes are added.
const _AXES_BY_NAME = Dict{String,Any}(
    "ICRF"     => ICRF(),
    "GCRF"     => GCRF(),
    "CIRS"     => CIRS(),
    "TIRS"     => TIRS(),
    "ITRF"     => ITRF(),
    "MJ2000Eq"    => MJ2000Eq(),
    "MJ2000Ec" => MJ2000Ec(),
    "MODEq"    => MODEq(),
    "TODEq"    => TODEq(),
    "MODEc"    => MODEc(),
    "TODEc"    => TODEc(),
    "PEF"      => PEF(),
    "MoonPA"   => MoonPA(),
    "MoonME"   => MoonME(),
    "IAU_SUN"     => CelestialBodyFixed(AstroUniverse.sun),
    "IAU_MERCURY" => CelestialBodyFixed(AstroUniverse.mercury),
    "IAU_VENUS"   => CelestialBodyFixed(AstroUniverse.venus),
    "IAU_MARS"    => CelestialBodyFixed(AstroUniverse.mars),
    "IAU_JUPITER" => CelestialBodyFixed(AstroUniverse.jupiter),
    "IAU_SATURN"  => CelestialBodyFixed(AstroUniverse.saturn),
    "IAU_URANUS"  => CelestialBodyFixed(AstroUniverse.uranus),
    "IAU_NEPTUNE" => CelestialBodyFixed(AstroUniverse.neptune),
    "IAU_PLUTO"   => CelestialBodyFixed(AstroUniverse.pluto),
)

_resolve_axes(name::AbstractString) = get(_AXES_BY_NAME, String(name)) do
    error("Unknown axes name '$name' in state-transform truth row.")
end

"""
    run_truth_row(row::AxesRotationTruthRow)

Compute the state transform for `row` and Frobenius-diff against `row.M`.
Records a `@test` assertion. Returns the numeric residual for logging.
"""
function run_truth_row(row::AxesRotationTruthRow)
    src = _resolve_axes(row.source_frame)
    tgt = _resolve_axes(row.target_frame)
    M_ours = axes_rotation(src, tgt, row.epoch_jd_tdb)
    resid = norm(M_ours - row.M)
    @test resid < row.tolerance
    return resid
end

"""
    run_truth_rows(rows::AbstractVector{AxesRotationTruthRow}; verbose=false)

Iterate a truth vector, run each row inside a nested `@testset` named by
`row.id`, and (optionally) print each residual. Returns a vector of residuals
in row order.
"""
function run_truth_rows(rows::AbstractVector{AxesRotationTruthRow}; verbose::Bool=false)
    residuals = Float64[]
    for row in rows
        @testset "$(row.id)" begin
            r = run_truth_row(row)
            push!(residuals, r)
            verbose && @info "state-transform truth row" id=row.id resid=r tol=row.tolerance
        end
    end
    return residuals
end
