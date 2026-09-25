# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# =============================================================================
# Setting one Keplerian element changes that element and nothing else.
#
# This is the property the shipped setters violated. Each rebuilt a
# `KeplerianState` positionally, and four of the five passed the arguments in
# the wrong order:
#
#     KeplerianState(a, e, i, raan, aop, ta)          # the constructor
#     KeplerianState(newval, s.ecc, s.inc, s.aop, s.raan, s.ta)   # sma, as shipped
#
# so `set_calc!(OrbitCalc(sat, SMA()), 8000.0)` also swapped RAAN and AOP —
# silently, because the resulting orbit is perfectly valid, just rotated in its
# plane and about its pole. Measured before the fix:
#
#     before : sma 6328.468  raan 5.034140  aop 5.346374
#     after  : sma 8000.000  raan 5.346374  aop 5.034140
#
# Only `raan` was correct, and only because its new value happens to land in
# slot four. There is no `aop` setter, which is why nothing caught it.
#
# The test is written as the invariant rather than against the four known
# cases, so any future permutation error fails here — including in a setter
# that does not exist yet.
# =============================================================================

using Test
using AstroCallbacks
using AstroModels
using AstroStates
using AstroEpochs
using AstroFrames
using AstroUniverse

const _ELEMENTS = (:sma, :ecc, :inc, :raan, :aop, :ta)

"""The subject's Keplerian elements, read straight from its state."""
function _elements(sat)
    cart = CartesianState(to_vector(AstroFrames.state_of(sat)))
    return KeplerianState(cart, get_gravparam(AstroFrames.frame_of(sat).origin))
end

"""A spacecraft whose six elements are all distinct, so a swap cannot hide."""
function _subject()
    return Spacecraft(CartesianState([7000.0, 1000.0, 2000.0, 1.5, 6.5, 1.0]),
                      Time(2458849.5, 0.0, :tdb, :jd);
                      coord_sys = CoordinateSystem(earth, ICRF()))
end

@testset "setting one element leaves the others alone" begin
    # `aop` is deliberately absent — there is no setter for it. Every element
    # that *can* be set is covered.
    settable = ((SMA(),  :sma,  8000.0),
                (Ecc(),  :ecc,  0.20),
                (Inc(),  :inc,  deg2rad(35.0)),
                (RAAN(), :raan, deg2rad(200.0)),
                (TA(),   :ta,   deg2rad(120.0)))

    for (tag, name, target) in settable
        @testset "$(name)" begin
            sat = _subject()
            before = _elements(sat)

            set_calc!(OrbitCalc(sat, tag), target)
            after = _elements(sat)

            # The one asked for moved to the value asked for.
            @test getfield(after, name) ≈ target atol = 1e-9

            # Every other element is untouched. This is the assertion the
            # shipped code failed: RAAN and AOP traded places.
            for other in _ELEMENTS
                other === name && continue
                @test getfield(after, other) ≈ getfield(before, other) atol = 1e-9
            end
        end
    end
end

@testset "the elements used are genuinely distinct" begin
    # The invariant above is only meaningful if no two elements start equal —
    # otherwise a swap between them would pass. Guards the fixture.
    k = _elements(_subject())
    values = [getfield(k, e) for e in _ELEMENTS]
    for i in eachindex(values), j in eachindex(values)
        i < j || continue
        @test !isapprox(values[i], values[j]; atol = 1e-6)
    end
end

@testset "reading an element round-trips through set" begin
    # Setting an element to the value it already has must be a no-op. A
    # permutation error shows up here too, and it needs no target value.
    for (tag, name) in ((SMA(), :sma), (Ecc(), :ecc), (Inc(), :inc),
                        (RAAN(), :raan), (TA(), :ta))
        sat = _subject()
        before = _elements(sat)

        set_calc!(OrbitCalc(sat, tag), getfield(before, name))
        after = _elements(sat)

        for e in _ELEMENTS
            @test getfield(after, e) ≈ getfield(before, e) atol = 1e-9
        end
    end
end
