# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# One quantity per element of a state representation.
#
# Truth, for the readers: **closed form**. Each element is computed here from the Keplerian
# elements or the Cartesian state by its textbook definition, independently of the AstroStates
# conversion the reader goes through. The Brouwer mean elements have no closed form, so they are
# checked against the AstroStates conversion they read, which has its own tests against GMAT.
#
# Truth, for the setters: **round trip**. What is written reads back, in the coordinate system it
# was written in, and the rest of its representation reads back unchanged. Each is checked in the
# spacecraft's own frame and in another one, and through `set_quantity!`.

using Test
using Logging
using LinearAlgebra
using AstroCallbacks
using AstroCallbacks: _STATE_ELEMENT_QUANTITIES, _BROUWER_QUANTITIES, subject_type
using AstroStates
using AstroEpochs
using AstroFrames
using AstroModels
using AstroUniverse
using EpicycleBase: label, unit, is_cyclic, cycle, is_settable, set_quantity!

const _SE_EPOCH = Time(2458849.5, 0.0, TDB(), JD())
const _SE_μ     = earth.mu

_se_closed()     = Spacecraft(CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]), _SE_EPOCH;
                              coord_sys = EarthMJ2000Eq)
_se_hyperbolic() = Spacecraft(CartesianState([7000.0, 300.0, 0.0, 0.0, 12.5, 1.0]), _SE_EPOCH;
                              coord_sys = EarthMJ2000Eq)

const _SE_ALL = (_STATE_ELEMENT_QUANTITIES..., _BROUWER_QUANTITIES...)

_angle_diff(a, b) = abs(rem2pi(a - b, RoundNearest))

"`a` and `b` agree, as angles modulo 2π when `cyclic`."
_agree(a, b, cyclic; atol) = cyclic ? _angle_diff(a, b) <= atol : isapprox(a, b; atol = atol)

_rv(sat, cs) = (position_vector(sat, cs), velocity_vector(sat, cs))

# --- the readers, against closed forms ---------------------------------------------------------

@testset "the spherical elements are the directions of r and v" begin
    for cs in (EarthMJ2000Eq, EarthMJ2000Ec)
        sat  = _se_closed()
        r, v = _rv(sat, cs)
        λ, δ = mod(atan(r[2], r[1]), 2π), asin(r[3] / norm(r))
        @test right_ascension(sat, cs) ≈ λ atol = 1e-12
        @test declination(sat, cs) ≈ δ atol = 1e-12
        @test velocity_right_ascension(sat, cs) ≈ mod(atan(v[2], v[1]), 2π) atol = 1e-12
        @test velocity_declination(sat, cs) ≈ asin(v[3] / norm(v)) atol = 1e-12
        @test flight_path_angle(sat, cs) ≈ acos(dot(r, v) / (norm(r) * norm(v))) atol = 1e-12

        east  = [-sin(λ), cos(λ), 0.0]
        north = [-sin(δ) * cos(λ), -sin(δ) * sin(λ), cos(δ)]
        @test velocity_azimuth(sat, cs) ≈ mod(atan(dot(v, east), dot(v, north)), 2π) atol = 1e-12
    end
end

@testset "the equinoctial sets follow from the Keplerian elements" begin
    for cs in (EarthMJ2000Eq, EarthMJ2000Ec)
        sat = _se_closed()
        a, e, i = semi_major_axis(sat, cs), eccentricity(sat, cs), inclination(sat, cs)
        Ω, ω, M = raan(sat, cs), argument_of_periapsis(sat, cs), mean_anomaly(sat, cs)
        ν = true_anomaly(sat, cs)

        @test equinoctial_h(sat, cs) ≈ e * sin(ω + Ω) atol = 1e-12
        @test equinoctial_k(sat, cs) ≈ e * cos(ω + Ω) atol = 1e-12
        @test equinoctial_p(sat, cs) ≈ tan(i / 2) * sin(Ω) atol = 1e-12
        @test equinoctial_q(sat, cs) ≈ tan(i / 2) * cos(Ω) atol = 1e-12
        @test alternate_equinoctial_p(sat, cs) ≈ sin(i / 2) * sin(Ω) atol = 1e-12
        @test alternate_equinoctial_q(sat, cs) ≈ sin(i / 2) * cos(Ω) atol = 1e-12
        @test _angle_diff(mean_longitude(sat, cs), Ω + ω + M) < 1e-10

        @test modified_equinoctial_f(sat, cs) ≈ e * cos(ω + Ω) atol = 1e-12
        @test modified_equinoctial_g(sat, cs) ≈ e * sin(ω + Ω) atol = 1e-12
        @test modified_equinoctial_h(sat, cs) ≈ tan(i / 2) * cos(Ω) atol = 1e-12
        @test modified_equinoctial_k(sat, cs) ≈ tan(i / 2) * sin(Ω) atol = 1e-12
        @test semilatus_rectum(sat, cs) ≈ a * (1 - e^2) rtol = 1e-12
        @test _angle_diff(true_longitude(sat, cs), Ω + ω + ν) < 1e-10
    end
end

@testset "the Delaunay elements follow from the Keplerian elements" begin
    for cs in (EarthMJ2000Eq, EarthMJ2000Ec)
        sat = _se_closed()
        a, e, i = semi_major_axis(sat, cs), eccentricity(sat, cs), inclination(sat, cs)
        L, G    = sqrt(_SE_μ * a), sqrt(_SE_μ * a * (1 - e^2))
        @test _angle_diff(delaunay_l(sat, cs), mean_anomaly(sat, cs)) < 1e-10
        @test _angle_diff(delaunay_g(sat, cs), argument_of_periapsis(sat, cs)) < 1e-10
        @test _angle_diff(delaunay_h(sat, cs), raan(sat, cs)) < 1e-10
        @test delaunay_L(sat, cs) ≈ L rtol = 1e-12
        @test delaunay_G(sat, cs) ≈ G rtol = 1e-12
        @test delaunay_G(sat, cs) ≈ angular_momentum_magnitude(sat, cs) rtol = 1e-12
        @test delaunay_H(sat, cs) ≈ G * cos(i) rtol = 1e-12
    end
end

@testset "the asymptotes are the directions at true anomaly ±ν∞" begin
    for cs in (EarthMJ2000Eq, EarthMJ2000Ec)
        sat  = _se_hyperbolic()
        e    = eccentricity(sat, cs)
        ν∞   = acos(-1 / e)
        P    = eccentricity_vector(sat, cs) ./ e
        ĥ    = orbit_normal(sat, cs)
        Q    = cross(ĥ, P)
        S_out = cos(ν∞) .* P .+ sin(ν∞) .* Q
        @test outgoing_asymptote_right_ascension(sat, cs) ≈ mod(atan(S_out[2], S_out[1]), 2π) atol = 1e-10
        @test outgoing_asymptote_declination(sat, cs) ≈ asin(S_out[3]) atol = 1e-10

        # The incoming asymptote points from the body back along the approach, at true anomaly −ν∞.
        S_in = cos(ν∞) .* P .- sin(ν∞) .* Q
        @test incoming_asymptote_right_ascension(sat, cs) ≈ mod(atan(S_in[2], S_in[1]), 2π) atol = 1e-10
        @test incoming_asymptote_declination(sat, cs) ≈ asin(S_in[3]) atol = 1e-10

        # The B-plane is built on the incoming velocity at infinity, Ŝ = −S_in. B points from the
        # centre to where that line pierces the plane, B̂ = Ŝ × ĥ, and the azimuth runs from R̂
        # toward T̂.
        Ŝ = -S_in
        T̂ = normalize(cross(Ŝ, [0.0, 0.0, 1.0]))
        R̂ = cross(Ŝ, T̂)
        B̂ = cross(Ŝ, ĥ)
        @test _angle_diff(incoming_b_vector_azimuth(sat, cs), atan(dot(B̂, T̂), dot(B̂, R̂))) < 1e-10
    end
end

@testset "the Brouwer mean elements are the AstroStates conversion's" begin
    sat  = _se_closed()
    cart = CartesianState(Coordinate(sat, EarthMJ2000Ec))
    for (fn, S, field, _...) in _BROUWER_QUANTITIES
        truth = getfield(getfield(AstroStates, S)(cart, _SE_μ), field)
        @test getfield(AstroCallbacks, fn)(sat, EarthMJ2000Ec) ≈ truth rtol = 1e-12
    end
    # Mean and osculating differ by the short-period J2 terms, kilometres in a, not more.
    @test abs(brouwer_long_semi_major_axis(sat) - semi_major_axis(sat)) < 100.0
    @test abs(brouwer_short_semi_major_axis(sat) - semi_major_axis(sat)) < 100.0
end

@testset "a Coordinate reads the same as the spacecraft it was taken from" begin
    sat   = _se_closed()
    coord = Coordinate(sat, EarthMJ2000Eq)
    for (fn, _, _, _, _, _, _, _, _, _, hyperbolic) in _SE_ALL
        hyperbolic && continue
        f = getfield(AstroCallbacks, fn)
        @test f(coord) == f(sat)
    end
end

@testset "a reader refuses an orbit its representation does not describe" begin
    hyp = _se_hyperbolic()
    with_logger(NullLogger()) do
        for f in (equinoctial_h, mean_longitude, alternate_equinoctial_p, delaunay_L, delaunay_l)
            @test_throws ArgumentError f(hyp)
        end
        circular = Spacecraft(CartesianState([7000.0, 0.0, 0.0, 0.0, sqrt(_SE_μ / 7000.0), 0.0]),
                              _SE_EPOCH; coord_sys = EarthMJ2000Eq)
        @test_throws ArgumentError outgoing_asymptote_right_ascension(circular)
    end
end

# --- the setters --------------------------------------------------------------------------------

"Fields of representation `S` of `sat` in `cs`."
_fields(sat, S, cs) =
    (s = getfield(AstroStates, S)(CartesianState(Coordinate(sat, cs)), _SE_μ);
     [getfield(s, f) for f in fieldnames(typeof(s))])

@testset "each setter writes its element and holds the rest of its representation" begin
    for (fn, S, field, _, _, cyclic, _, _, _, example, hyperbolic) in _SE_ALL,
            cs in (EarthMJ2000Eq, EarthMJ2000Ec)
        @testset "$(fn)!, $(cs.axes)" begin
            brouwer = S in (:BrouwerMeanLongState, :BrouwerMeanShortState)
            tol     = brouwer ? 1e-6 : 1e-9
            f, f!   = getfield(AstroCallbacks, fn), getfield(AstroCallbacks, Symbol(fn, "!"))
            sat     = hyperbolic ? _se_hyperbolic() : _se_closed()
            names   = fieldnames(getfield(AstroStates, S))
            k       = findfirst(==(field), names)
            to      = parse(Float64, example)
            before  = _fields(sat, S, cs)

            @test f!(sat, cs; to = to) === sat
            @test _agree(f(sat, cs), to, cyclic; atol = tol * max(1, abs(to)))
            after = _fields(sat, S, cs)
            for j in eachindex(names)
                j == k && continue
                scale = max(1, abs(before[j]))
                @test _agree(after[j], before[j], true; atol = tol * scale) ||
                      isapprox(after[j], before[j]; atol = tol * scale)
            end
            @test frame_of(sat) === EarthMJ2000Eq
        end
    end
end

@testset "set_quantity! is the setter, and only a spacecraft is settable" begin
    for (fn, _, _, _, _, _, _, _, _, example, hyperbolic) in _SE_ALL
        f, f! = getfield(AstroCallbacks, fn), getfield(AstroCallbacks, Symbol(fn, "!"))
        make  = hyperbolic ? _se_hyperbolic : _se_closed
        a, b  = make(), make()
        to    = parse(Float64, example)
        f!(a, EarthMJ2000Ec; to = to)
        set_quantity!(b, f, EarthMJ2000Ec; to = to)
        @test state(a) == state(b)
        @test is_settable(a, f)
        @test is_settable(a, f, EarthMJ2000Ec)
        @test !is_settable(Coordinate(a, EarthMJ2000Eq), f)
    end
end

@testset "a setter refuses elements that describe no orbit, and leaves the state" begin
    cases = (("declination past the pole",      _se_closed,     declination!,                     2.0),
             ("velocity declination past the pole", _se_closed, velocity_declination!,           -2.0),
             ("flight path angle beyond π",     _se_closed,     flight_path_angle!,               4.0),
             ("equinoctial eccentricity of 1.2", _se_closed,    equinoctial_h!,                   1.2),
             ("Delaunay G above L",             _se_closed,     delaunay_G!,                      1e6),
             ("Delaunay H above G",             _se_closed,     delaunay_H!,                      1e6),
             ("negative semi-latus rectum",     _se_closed,     semilatus_rectum!,               -1.0),
             ("asymptote declination past the pole", _se_hyperbolic, outgoing_asymptote_declination!, 2.0),
             ("Brouwer eccentricity of 0.995",  _se_closed,     brouwer_long_eccentricity!,       0.995),
             ("Brouwer periapsis below 3000 km", _se_closed,    brouwer_short_semi_major_axis!,   2000.0),
             ("a right ascension that is NaN",  _se_closed,     right_ascension!,                 NaN))
    for (name, make, set!, value) in cases
        @testset "$name" begin
            sat = make()
            x0  = state(sat)
            @test_throws ArgumentError set!(sat; to = value)
            @test state(sat) == x0
        end
    end
end

@testset "semilatus_rectum! and true_longitude! hold the rest of the modified equinoctial set" begin
    for (set!, k, to) in ((semilatus_rectum!, 1, 7100.0), (true_longitude!, 6, 2.0))
        sat    = _se_closed()
        before = _fields(sat, :ModifiedEquinoctialState, EarthMJ2000Ec)
        set!(sat, EarthMJ2000Ec; to = to)
        after = _fields(sat, :ModifiedEquinoctialState, EarthMJ2000Ec)
        @test after[k] ≈ to rtol = 1e-10
        @test after[setdiff(1:6, k)] ≈ before[setdiff(1:6, k)] rtol = 1e-9 atol = 1e-12
    end
    @test is_settable(_se_closed(), semilatus_rectum)
    @test is_settable(_se_closed(), true_longitude)
end

# --- the velocity components --------------------------------------------------------------------

@testset "velocity_x, velocity_y and velocity_z read and set one component" begin
    for (k, (f, f!)) in enumerate(((velocity_x, velocity_x!), (velocity_y, velocity_y!),
                                   (velocity_z, velocity_z!))),
            cs in (EarthMJ2000Eq, EarthMJ2000Ec)
        sat = _se_closed()
        @test f(sat, cs) == velocity_vector(sat, cs)[k]
        x0 = vcat(_rv(sat, cs)...)
        f!(sat, cs; to = 1.25)
        x1 = vcat(_rv(sat, cs)...)
        @test x1[3 + k] ≈ 1.25 rtol = 1e-12
        @test x1[setdiff(1:6, 3 + k)] ≈ x0[setdiff(1:6, 3 + k)] rtol = 1e-12
        @test label(f) == "V" * "XYZ"[k:k]
        @test unit(f) == "km/s"
        @test is_settable(sat, f)
    end
    sat = _se_closed()
    @test output_partial(sat, velocity_y) == reshape(Float64[0, 0, 0, 0, 1, 0], 1, 6)
end

# --- traits and discovery -----------------------------------------------------------------------

@testset "each quantity carries its label, unit and wrap, and is found" begin
    found = Set(quantities(Spacecraft))
    for (fn, _, _, lbl, u, cyclic, _...) in _SE_ALL
        f = getfield(AstroCallbacks, fn)
        @test label(f) == lbl
        @test unit(f) == u
        @test is_cyclic(f) == cyclic
        cyclic && @test cycle(f) == 2π
        @test subject_type(f) == Union{Spacecraft, Coordinate}
        @test f in found
    end
    for f in (velocity_x, velocity_y, velocity_z)
        @test f in found
    end
end
