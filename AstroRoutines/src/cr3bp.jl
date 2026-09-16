# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# Circular restricted three-body problem, in the rotating frame and normalized
# units: the two primaries sit at (-mu, 0, 0) and (1-mu, 0, 0), their separation
# is 1, their orbital rate is 1, and the total mass is 1. Everything here is a
# function of the state and the mass ratio alone, which is what keeps this
# package free of dependencies.

"""
    cr3bp_mass_ratio(m1::Real, m2::Real)

Mass ratio of a circular restricted three-body system.

# Arguments
- `m1`: Mass of the larger primary, in any unit.
- `m2`: Mass of the smaller primary, in the same unit.

# Returns
`μ = m2 / (m1 + m2)`, the dimensionless mass ratio the other routines here take.
Gravitational parameters may be passed instead of masses; the ratio is the same.

# Notes
Throws `ArgumentError` unless both masses are positive.

# Example
```jldoctest
cr3bp_mass_ratio(3.986004418e5, 4.9028e3)   # Earth-Moon, from the GMs

# output

0.012150583916324809
```
"""
function cr3bp_mass_ratio(m1::Real, m2::Real)
    (m1 > 0 && m2 > 0) || throw(ArgumentError(
        "cr3bp_mass_ratio: both masses must be positive; got m1 = $(m1), m2 = $(m2)"))
    return m2 / (m1 + m2)
end

# Distances to the two primaries. Everything below is written in terms of these.
@inline function _cr3bp_radii(x, y, z, mu)
    r1 = sqrt((x + mu)^2     + y^2 + z^2)
    r2 = sqrt((x - 1 + mu)^2 + y^2 + z^2)
    return r1, r2
end

"""
    cr3bp_accel(s::AbstractVector, mu::Real)

Acceleration of the third body in the rotating frame of a circular restricted
three-body system.

# Arguments
- `s`: State `[x, y, z, ẋ, ẏ, ż]` in normalized units.
- `mu`: Mass ratio, `0 < μ < 1`. See [`cr3bp_mass_ratio`](@ref).

# Returns
The three-vector `[ẍ, ÿ, z̈]` in normalized units, carrying the Coriolis,
centrifugal and gravitational terms and no thrust. A control acceleration is
added by the caller.

# Notes
- Differentiable with respect to the state and `μ` (ForwardDiff).
- Singular at the primaries themselves, where `r1` or `r2` is zero.

# Example
```julia
cr3bp_accel([0.8159625214638456, 0.0, 0.0, 0.0, 0.2072212474921765, 0.0],
            0.01215058560962404)
```
"""
function cr3bp_accel(s::AbstractVector, mu::Real)
    x, y, z, vx, vy, _ = s[1], s[2], s[3], s[4], s[5], s[6]
    r1, r2 = _cr3bp_radii(x, y, z, mu)
    c1 = (1 - mu) / r1^3
    c2 = mu / r2^3
    return [ 2vy + x - c1 * (x + mu)     - c2 * (x - 1 + mu),
            -2vx + y - c1 * y            - c2 * y,
                     - c1 * z            - c2 * z ]
end

"""
    cr3bp_eom!(ds, s::AbstractVector, mu::Real)

Write the time derivative of a circular restricted three-body state into `ds`.

# Arguments
- `ds`: Six-element output, overwritten.
- `s`: State `[x, y, z, ẋ, ẏ, ż]` in normalized units.
- `mu`: Mass ratio, `0 < μ < 1`.

# Returns
`nothing`. The derivative is written into `ds`.

# Notes
Ballistic motion only. A propagator that applies thrust adds the control
acceleration to `ds[4:6]` after calling this.

# Example
```julia
ds = zeros(6)
cr3bp_eom!(ds, [0.8159625214638456, 0.0, 0.0, 0.0, 0.2072212474921765, 0.0],
           0.01215058560962404)
ds
```
"""
function cr3bp_eom!(ds, s::AbstractVector, mu::Real)
    ds[1] = s[4]
    ds[2] = s[5]
    ds[3] = s[6]
    a = cr3bp_accel(s, mu)
    ds[4] = a[1]
    ds[5] = a[2]
    ds[6] = a[3]
    return nothing
end

"""
    cr3bp_jacobian(s::AbstractVector, mu::Real)

Jacobian of the circular restricted three-body equations of motion with respect
to the state.

# Arguments
- `s`: State `[x, y, z, ẋ, ẏ, ż]` in normalized units.
- `mu`: Mass ratio, `0 < μ < 1`.

# Returns
The 6x6 matrix `∂ṡ/∂s`. The upper right block is the identity, the lower left is
the Hessian of the pseudo-potential, and the lower right carries the Coriolis
terms.

# Notes
- Analytic, not a finite difference, so it carries no step-size error.
- Differentiable with respect to the state and `μ` (ForwardDiff).
- Singular at the primaries themselves, where `r1` or `r2` is zero.
- This is the `A(s)` of the variational equation `Φ̇ = A(s) Φ` that
  [`cr3bp_stm_eom!`](@ref) evaluates.

# Example
```julia
cr3bp_jacobian([0.8159625214638456, 0.0, 0.0, 0.0, 0.2072212474921765, 0.0],
               0.01215058560962404)
```
"""
function cr3bp_jacobian(s::AbstractVector, mu::Real)
    x, y, z = s[1], s[2], s[3]
    r1, r2  = _cr3bp_radii(x, y, z, mu)
    a1, a2  = (1 - mu) / r1^3, mu / r2^3
    b1, b2  = 3 * (1 - mu) / r1^5, 3 * mu / r2^5
    dx1, dx2 = x + mu, x - 1 + mu

    # Second derivatives of the pseudo-potential. The +1 terms on xx and yy are
    # the centrifugal contribution; z has none, which is why out-of-plane motion
    # is purely oscillatory near the primaries.
    Uxx = 1 - a1 - a2 + b1 * dx1^2 + b2 * dx2^2
    Uyy = 1 - a1 - a2 + b1 * y^2   + b2 * y^2
    Uzz =   - a1 - a2 + b1 * z^2   + b2 * z^2
    Uxy = b1 * dx1 * y + b2 * dx2 * y
    Uxz = b1 * dx1 * z + b2 * dx2 * z
    Uyz = b1 * y   * z + b2 * y   * z

    T = promote_type(eltype(s), typeof(mu))
    A = zeros(T, 6, 6)
    A[1, 4] = one(T); A[2, 5] = one(T); A[3, 6] = one(T)
    A[4, 1] = Uxx; A[4, 2] = Uxy; A[4, 3] = Uxz; A[4, 5] =  2one(T)
    A[5, 1] = Uxy; A[5, 2] = Uyy; A[5, 3] = Uyz; A[5, 4] = -2one(T)
    A[6, 1] = Uxz; A[6, 2] = Uyz; A[6, 3] = Uzz
    return A
end

"""
    jacobi_constant(s::AbstractVector, mu::Real)

Jacobi constant of a state in the circular restricted three-body problem.

# Arguments
- `s`: State `[x, y, z, ẋ, ẏ, ż]` in normalized units.
- `mu`: Mass ratio, `0 < μ < 1`.

# Returns
`C = 2U − v²`, where `U` is the pseudo-potential.

# Notes
- Constant along a ballistic arc, and changes only when something does work on
  the spacecraft. Holding it along an integrated arc checks the integration
  tolerance; along a low-thrust arc its change measures the energy the thrust
  supplied.
- Sign convention: `C` decreases as speed increases, and the forbidden region
  for a given `C` is where `2U < C`.
- Differentiable with respect to the state and `μ` (ForwardDiff).

# Example
```jldoctest
jacobi_constant([0.8159625214638456, 0.0, 0.0, 0.0, 0.2072212474921765, 0.0],
                0.01215058560962404)

# output

3.150016832809122
```
"""
function jacobi_constant(s::AbstractVector, mu::Real)
    x, y, z, vx, vy, vz = s[1], s[2], s[3], s[4], s[5], s[6]
    r1, r2 = _cr3bp_radii(x, y, z, mu)
    U = (x^2 + y^2) / 2 + (1 - mu) / r1 + mu / r2
    return 2U - (vx^2 + vy^2 + vz^2)
end

"""
    libration_point(mu::Real, which::Symbol)

Position of one of the five libration points of a circular restricted
three-body system.

# Arguments
- `mu`: Mass ratio, `0 < μ < 1`.
- `which`: `:L1`, `:L2`, `:L3`, `:L4` or `:L5`.

# Returns
The three-vector `[x, y, 0]` in normalized rotating-frame units.

# Notes
- `L1`, `L2` and `L3` lie on the x-axis and are iterated to a residual below
  `1e-14`, so they carry that much error rather than being exact. `L4` and `L5`
  are closed form, at `x = 1/2 − μ` and `y = ±√3/2`.
- Throws `ArgumentError` for `μ ∉ (0, 1)` or a symbol outside the five.

# Example
```julia
libration_point(0.01215058560962404, :L1)   # Earth-Moon L1
```
"""
function libration_point(mu::Real, which::Symbol)
    (mu > 0 && mu < 1) || throw(ArgumentError(
        "libration_point: mass ratio must be in (0, 1); got mu = $(mu)"))

    if which === :L4 || which === :L5
        y = which === :L4 ? sqrt(3) / 2 : -sqrt(3) / 2
        return [1 / 2 - mu, y, zero(mu)]
    end
    which in (:L1, :L2, :L3) || throw(ArgumentError(
        "libration_point: which must be one of :L1, :L2, :L3, :L4, :L5; got :$(which)"))

    # Collinear points are the roots of dU/dx = 0 on the x-axis. Each branch
    # carries the sign of (x + mu) and (x - 1 + mu), which is what distinguishes
    # the three, so the residual is written per branch rather than generically.
    rh = cbrt(mu / 3)                       # hill radius, the classical start
    x  = which === :L1 ? 1 - mu - rh :
         which === :L2 ? 1 - mu + rh :
                        -1 - 5mu / 12       # L3, on the far side of the larger body

    # f(x) = x - (1-mu)(x+mu)/|x+mu|^3 - mu(x-1+mu)/|x-1+mu|^3, the x-component
    # of the pseudo-potential gradient on the axis. Since d/du (u/|u|^3) is
    # -2/|u|^3, f' carries no sign factor and is positive everywhere, which is
    # what keeps each branch converging to its own root.
    for _ in 1:100
        d1 = x + mu
        d2 = x - 1 + mu
        f  = x - (1 - mu) * sign(d1) / d1^2 - mu * sign(d2) / d2^2
        df = 1 + 2 * (1 - mu) / abs(d1)^3 + 2 * mu / abs(d2)^3
        dx = f / df
        x -= dx
        abs(dx) < 1e-14 && break
    end
    return [x, zero(x), zero(x)]
end

"""
    cr3bp_stm_eom!(dz, z::AbstractVector, mu::Real)

Write the derivative of a state augmented with its state transition matrix.

# Arguments
- `dz`: 42-element output, overwritten.
- `z`: `[s; vec(Φ)]` — the six-element state followed by the 6x6 state
  transition matrix in column-major order, 42 elements in all.
- `mu`: Mass ratio, `0 < μ < 1`.

# Returns
`nothing`. The derivative is written into `dz`.

# Notes
- The variational block is `Φ̇ = A(s) Φ`, where `A` is
  [`cr3bp_jacobian`](@ref). This function evaluates that derivative; integration
  is the caller's, so the same integrator advances the state and `Φ` on the same
  steps and to the same tolerance.
- Start from [`cr3bp_stm_initial`](@ref), which sets `Φ(0) = I`.
- The flow is Hamiltonian, so `Φ` stays symplectic and `det Φ = 1` along an exact
  arc. That determinant is the well-conditioned check on an integration
  tolerance; the symplectic residual scales with the square of the matrix norm,
  which is large on an unstable orbit.
- Throws `ArgumentError` unless `z` has 42 elements.

# Example
```julia
z0 = cr3bp_stm_initial([0.8159625214638456, 0.0, 0.0, 0.0, 0.2072212474921765, 0.0])
dz = zeros(42)
cr3bp_stm_eom!(dz, z0, 0.01215058560962404)
dz[1:6]
```
"""
function cr3bp_stm_eom!(dz, z::AbstractVector, mu::Real)
    length(z) == 42 || throw(ArgumentError(
        "cr3bp_stm_eom!: z must be the 6 state elements followed by the 36 of " *
        "the transition matrix, 42 in all; got $(length(z))"))

    s = @view z[1:6]
    cr3bp_eom!(view(dz, 1:6), s, mu)

    A = cr3bp_jacobian(s, mu)
    # Phi is stored column-major, so element (i, j) is z[6 + i + 6(j-1)]. The
    # product is written out rather than reshaped, to keep this allocation-free
    # inside an integrator's inner loop.
    @inbounds for j in 1:6, i in 1:6
        acc = zero(eltype(dz))
        for k in 1:6
            acc += A[i, k] * z[6 + k + 6 * (j - 1)]
        end
        dz[6 + i + 6 * (j - 1)] = acc
    end
    return nothing
end

"""
    cr3bp_stm_initial(s::AbstractVector)

Build the augmented initial condition [`cr3bp_stm_eom!`](@ref) expects.

# Arguments
- `s`: The six-element state the arc starts from.

# Returns
A 42-element vector: the state followed by the identity matrix flattened
column-major, which is the state transition matrix at the start of any arc.

# Notes
Throws `ArgumentError` unless `s` has six elements.

# Example
```julia
cr3bp_stm_initial([0.8159625214638456, 0.0, 0.0, 0.0, 0.2072212474921765, 0.0])
```
"""
function cr3bp_stm_initial(s::AbstractVector)
    length(s) == 6 || throw(ArgumentError(
        "cr3bp_stm_initial: state must have six elements; got $(length(s))"))
    z = zeros(eltype(s), 42)
    z[1:6] = s
    for i in 1:6
        z[6 + i + 6 * (i - 1)] = one(eltype(s))
    end
    return z
end
