# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# The state types: named elements in the documented order [AS-1, AS-11], no assumed central
# body [AS-3], and a display with labeled elements and angles in degrees [AS-12].

using Test
using AstroStates

# Each type's elements in the order `to_vector` documents, and which of them are angles.
const _ELEMENTS = [
    (CartesianState,            (:x, :y, :z, :vx, :vy, :vz),           ()),
    (KeplerianState,            (:sma, :ecc, :inc, :raan, :aop, :ta),  (3, 4, 5, 6)),
    (ModifiedKeplerianState,    (:rp, :ra, :inc, :raan, :aop, :ta),    (3, 4, 5, 6)),
    (EquinoctialState,          (:a, :h, :k, :p, :q, :mlong),          (6,)),
    (AlternateEquinoctialState, (:a, :h, :k, :altp, :altq, :mlong),    (6,)),
    (ModifiedEquinoctialState,  (:p, :f, :g, :h, :k, :L),              (6,)),
    (SphericalRADECState,       (:r, :ra, :dec, :v, :rav, :decv),      (2, 3, 5, 6)),
    (SphericalAZIFPAState,      (:r, :ra, :dec, :v, :vazi, :fpa),      (2, 3, 5, 6)),
    (IncomingAsymptoteState,    (:rp, :c3, :rla, :dla, :bpa, :ta),     (3, 4, 5, 6)),
    (OutGoingAsymptoteState,    (:rp, :c3, :rla, :dla, :bpa, :ta),     (3, 4, 5, 6)),
    (BrouwerMeanShortState,     (:sma, :ecc, :inc, :raan, :aop, :ma),  (3, 4, 5, 6)),
    (BrouwerMeanLongState,      (:sma, :ecc, :inc, :raan, :aop, :ma),  (3, 4, 5, 6)),
]

_element(s::CartesianState, i) = i <= 3 ? s.position[i] : s.velocity[i - 3]
_element(s, i) = getfield(s, fieldnames(typeof(s))[i])

@testset "Each type's elements are named and in the documented order [AS-1, AS-11]" begin
    @testset "$(nameof(T))" for (T, names, _) in _ELEMENTS
        s = T(collect(1.0:6.0))
        @test to_vector(s) == collect(1.0:6.0)
        T === CartesianState || @test fieldnames(T) == names
        @test [_element(s, i) for i in 1:6] == collect(1.0:6.0)
    end
end

@testset "No conversion assumes a central body [AS-3]" begin
    kinematic = (CartesianState, SphericalRADECState, SphericalAZIFPAState)
    states = [T(collect(1.0:6.0)) for (T, _, _) in _ELEMENTS]
    for (T, _, _) in _ELEMENTS, s in states
        typeof(s) <: T && continue                       # the identity constructor
        T in kinematic && typeof(s).name.wrapper in kinematic && continue
        @test !hasmethod(T, Tuple{typeof(s)})
    end
end

@testset "Each state prints its elements labeled, angles in degrees [AS-12]" begin
    @testset "$(nameof(T))" for (T, names, angles) in _ELEMENTS
        text = sprint(show, T(ones(6)))
        lines = split(chomp(text), '\n')
        @test lines[1] == "$(nameof(T)):"
        @test length(lines) == 7
        for i in 1:6
            label, value = strip.(split(lines[i + 1], '='))
            @test startswith(label, string(names[i]))
            if i in angles
                @test occursin("deg", label)
                @test isapprox(parse(Float64, value), rad2deg(1.0); atol = 1e-5)   # one radian, in degrees
            else
                @test !occursin("deg", label)
                @test parse(Float64, value) == 1.0
            end
        end
    end
end

@testset "CartesianState equality, hashing and copies" begin
    a = CartesianState([7000.0, 0.0, -0.0, 0.0, 7.5, 0.0])
    b = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0])
    @test a == b                                   # numerically equal
    @test !isequal(a, b)                           # but -0.0 and 0.0 differ under isequal, as for Float64
    @test isequal(a, CartesianState(to_vector(a)))
    @test hash(a) == hash(CartesianState(to_vector(a)))   # isequal states hash alike
    @test length(Set([b, CartesianState(to_vector(b))])) == 1

    # the identity conversion of a mutable state is a copy
    c = CartesianState(b, 398600.4415)
    @test c == b && c !== b
    c.position = c.position .+ 1.0
    @test b.position[1] == 7000.0
end

@testset "State types take any vector, any state and an OrbitState [AS-1, AS-2]" begin
    μ = 398600.4415
    using StaticArrays
    v = [7000.0, 300.0, 0.0, 0.0, 7.5, 1.0]
    @test CartesianState(SVector{6}(v), μ) == CartesianState(v)
    @test CartesianState(view(v, 1:6), μ) == CartesianState(v)
    for (T, _, _) in _ELEMENTS
        s = T(collect(1.0:6.0))
        @test T(s) == s || T(s) === s            # the identity form without μ, for every type
    end
    cart = CartesianState(v)
    kep = KeplerianState(cart, μ)
    os = OrbitState(to_vector(cart), Cartesian())
    @test KeplerianState(os, μ) ≈ kep
    @test SphericalRADECState(os) ≈ SphericalRADECState(cart)
    @test CartesianState(OrbitState(to_vector(kep), Keplerian()), μ) ≈ cart
end

@testset "Float32 stays Float32 through the alternate equinoctial set" begin
    alt = Float32[7000, 0.01, 0.02, 0.1, 0.05, 1.0]
    @test eltype(alt_equinoctial_to_equinoctial(alt)) == Float32
    @test eltype(equinoctial_to_alt_equinoctial(alt)) == Float32
end
