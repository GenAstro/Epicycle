# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0
# Process-noise models map local solve-for noise into discrete filter covariance.

module ProcessNoiseModels

using LinearAlgebra

export ProcessNoiseModel, NoNoise, DiagonalSNC, discretize

"""
    ProcessNoiseModel

Extension interface (abstract type) for filter process noise.

# Notes
A subtype implements `discretize(model, Δt, n_local) -> (G, q)`, where `n_local` is the
length of the solved-for quantity's slot in the filter state, `G` is `n_local × nq` selecting
which slots each noise channel reaches, and `q` holds the `nq` discrete variances over `Δt`.

The filter stacks the `(G, q)` pairs from every solved-for quantity into one `(G, Q)` and
hands them to the Thornton time update, so a model never sees the full state and cannot
couple two quantities. A model that needs that coupling is a different interface.

`NoNoise` and `DiagonalSNC` ship today.
"""
abstract type ProcessNoiseModel end

# ── No noise ────────────────────────────────────────────────────────────────
"""
    NoNoise()

Process noise model that adds none, for a quantity treated as deterministic.

# Notes
`discretize` returns an `n_local × 0` mapping and an empty variance vector, which the Thornton
time update accepts as the deterministic case. This is the default for a solved-for quantity
that is a constant of the problem rather than a state that wanders.

# Example
```julia
using AstroSolve
G, q = AstroSolve.discretize(NoNoise(), 60.0, 6)
```
"""
struct NoNoise <: ProcessNoiseModel end

discretize(::NoNoise, Δt::Real, n_local::Integer) =
    (zeros(Float64, n_local, 0), Float64[])

# ── Diagonal SNC ────────────────────────────────────────────────────────────
"""
    DiagonalSNC(psd)

State noise compensation applied independently to each element of a solved-for quantity.

# Fields
- `psd::Vector{Float64}`: White-noise power spectral density per element, in (element units)² per unit time.
  One entry per element of the quantity's slot.

# Notes
The discrete variance accumulated over a step is `psd[i] * Δt`, which is the continuous white
noise assumption integrated over the step rather than a variance quoted per step. Doubling the
step doubles the variance.

Zero entries are pruned rather than carried as zero-variance channels, so a quantity with
noise on some elements and not others produces a `G` narrower than its slot. An all-zero `psd`
is therefore equivalent to [`NoNoise`](@ref).

Independent per element: this model cannot express correlation between two elements of the
same quantity, let alone between quantities.

# Example
```julia
using AstroSolve
# Position deterministic, velocity wandering at 1e-9 km²/s² per second.
snc = DiagonalSNC([0.0, 0.0, 0.0, 1e-9, 1e-9, 1e-9])
G, q = AstroSolve.discretize(snc, 60.0, 6)
```
"""
struct DiagonalSNC <: ProcessNoiseModel
    psd :: Vector{Float64}
end
DiagonalSNC(psd::AbstractVector) = DiagonalSNC(Float64.(collect(psd)))

"""
    discretize(model, Δt, n_local) -> (G, q)

Discretize a process noise model over one step.

# Arguments
- `model`: The process noise model.
- `Δt`: Step duration, in the filter's time unit. Must not be negative.
- `n_local`: Length of the solved-for quantity's slot in the filter state.

# Returns
`(G, q)` where `G` is `n_local × nq` selecting which slots each noise channel reaches and `q`
holds the `nq` discrete variances accumulated over `Δt`. Both are empty in the second
dimension for a deterministic model.

# Notes
Throws `ArgumentError` if the model's size disagrees with `n_local`, and `DomainError` for a
negative `Δt`. A zero `Δt` is allowed and yields zero variances, which is what a filter update
at the same epoch as the previous one needs.

# Example
```jldoctest
using AstroSolve
G, q = AstroSolve.discretize(DiagonalSNC([0.0, 2.0]), 3.0, 2)
(size(G), q)

# output
((2, 1), [6.0])
```
"""
function discretize(m::DiagonalSNC, Δt::Real, n_local::Integer)
    length(m.psd) == n_local || throw(ArgumentError(
        "discretize(DiagonalSNC, ...): psd must have one entry per element of the " *
        "quantity's slot, so length $(n_local); got $(length(m.psd))"))
    Δt >= 0 || throw(DomainError(Δt,
        "discretize(DiagonalSNC, ...): Δt must not be negative; got $(Δt)"))
    nz_idx = findall(>(0.0), m.psd)
    nq     = length(nz_idx)
    G      = zeros(Float64, n_local, nq)
    q      = zeros(Float64, nq)
    for (j, i) in enumerate(nz_idx)
        G[i, j] = 1.0
        q[j]    = m.psd[i] * Δt
    end
    return G, q
end

end # module ProcessNoiseModels
