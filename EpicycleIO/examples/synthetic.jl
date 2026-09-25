# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Made-up data for the examples.
#
# None of this is Epicycle output and none of it is physically right — it is shaped like real
# results so the examples exercise the interface without needing a propagation. When the
# interface is settled these arrays get replaced by `history(Calc(...), ...)` and nothing else
# in the examples changes. That is the point: the plotting side never knew the difference.

using Random

const _RNG = MersenneTwister(20260901)

"""Hours, altitude in km, and position as a vector of 3-element vectors."""
function fake_orbit(n = 400; hours = 6.0)
    t   = range(0, hours; length = n)
    ω   = 2π / 1.55                      # about a 93 minute period
    alt = 500 .+ 60 .* sin.(ω .* t)
    r   = [[(6878 + 60sin(ω*τ)) * cos(ω*τ),
            (6878 + 60sin(ω*τ)) * sin(ω*τ),
            900 * sin(0.5ω*τ)] for τ in t]
    return collect(t), alt, r
end

"""A second spacecraft, sampled differently and for a different span."""
function fake_orbit_b(n = 173; hours = 4.5)
    t   = range(0, hours; length = n)
    alt = 780 .+ 45 .* cos.(2π / 1.72 .* t .+ 0.7)
    return collect(t), alt
end

"""One ground station pass. Azimuth and elevation in radians."""
function fake_pass(n = 120)
    s  = range(0, 1; length = n)
    az = deg2rad.(35 .+ 210 .* s)                    # sweeping roughly N through S
    el = deg2rad.(78 .* sin.(π .* s) .+ 2)           # rises, culminates, sets
    return collect(az), collect(el)
end

"""Range residuals with a 1σ that grows between measurement passes."""
function fake_residuals(n = 240)
    t = collect(range(0, 12; length = n))
    σ = 0.004 .+ 0.010 .* abs.(sin.(2π .* t ./ 9))
    resid = σ .* randn(_RNG, n) .* 0.8
    return t, resid, σ
end

"""A porkchop-shaped grid. Departure days, arrival days, and a C3-like surface."""
function fake_porkchop(nx = 60, ny = 70)
    dep = collect(range(0, 120; length = nx))
    arr = collect(range(180, 400; length = ny))
    C3  = [12 + 0.004(d - 55)^2 + 0.0025(a - 300)^2 + 6sin(d/23)*cos(a/37)
           for a in arr, d in dep]                    # rows = arrival, cols = departure
    return dep, arr, C3
end

"""A ground track in degrees, wrapped the way a real one would be."""
function fake_groundtrack(n = 500; revs = 2.5)
    s   = range(0, revs; length = n)
    lat = 51.6 .* sin.(2π .* s)
    lon = mod.(-180 .+ 360 .* (0.85 .* s .% 1.0) .- 22 .* s, 360) .- 180
    return collect(lon), collect(lat)
end

"""Cost and constraint violation for a solver that converges."""
function fake_convergence(n = 25)
    it   = 1:n
    cost = 4.2 .* exp.(-0.28 .* it) .+ 0.61 .+ 0.02 .* randn(_RNG, n)
    viol = 3.0 .* exp.(-0.42 .* it) .+ 1e-6
    return collect(it), cost, viol
end

"""
Insertion altitude across many Monte Carlo runs, in km.

A dispersion is what a box plot is for: a spread of outcomes from repeated runs, not a
quantity sampled over time. Slightly skewed, because a real one usually is.
"""
function fake_monte_carlo(n = 500; target = 500.0)
    return target .+ 8.0 .* randn(_RNG, n) .+ 2.5 .* abs.(randn(_RNG, n))
end

"""
Normalised residuals for two measurement types: residual divided by its own sigma.

Dividing by sigma is what makes them comparable — range is in km and range-rate in km/s, so
the raw values share no axis. Normalised, both should look like a unit normal, and one that
comes out fat is a measurement whose noise model is wrong. That is the thing an analyst is
actually looking for.
"""
function fake_normalised_residuals(n = 400)
    range_n     = randn(_RNG, n)                              # behaving
    rangerate_n = 1.7 .* randn(_RNG, n) .+ 0.35               # too fat, and biased
    return range_n, rangerate_n
end
