# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Automatic differentiation through the quantities a stopping condition or a constraint reads.
#
# A solver differentiates through `StopAt(position_z, sat, EarthMJ2000Ec; ...)` and
# `Constraint(inclination, sat, EarthMJ2000Ec; ...)` with dual numbers, and a quantity read in a
# named frame goes through a frame transform first. Nothing else in the suite differentiates
# through one.
#
# Truth: central finite differences of the same quantity, step 1e-3 km in position and 1e-6 km/s
# in velocity. The frames are chosen so the transform is not the identity: ecliptic axes and the
# rotating Earth-fixed axes.

using Test
using ForwardDiff
using AstroEpochs, AstroStates, AstroFrames, AstroUniverse
using AstroCallbacks: position_z, inclination, velocity_magnitude, position_dot_velocity

const _QAD_EPOCH = Time(2458849.5, 0.0, TDB(), JD())
const _QAD_X     = [7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]
const _QAD_H     = [1e-3, 1e-3, 1e-3, 1e-6, 1e-6, 1e-6]

_qad_subject(x) = Coordinate(x, CoordinateSystem(earth, MJ2000Eq()), _QAD_EPOCH)

function _qad_central(f, x)
    g = zeros(length(x))
    for i in eachindex(x)
        e = zeros(length(x)); e[i] = _QAD_H[i]
        g[i] = (f(x .+ e) - f(x .- e)) / (2 * _QAD_H[i])
    end
    return g
end

@testset "AD through a quantity read in a named frame matches finite differences" begin
    cases = (("z in ecliptic axes",          x -> position_z(_qad_subject(x), EarthMJ2000Ec)),
             ("inclination in ecliptic axes", x -> inclination(_qad_subject(x), EarthMJ2000Ec)),
             ("speed in Earth-fixed axes",    x -> velocity_magnitude(_qad_subject(x), EarthFixed)),
             ("r·v in the subject's frame",   x -> position_dot_velocity(_qad_subject(x))))
    for (label, f) in cases
        @testset "$label" begin
            g_ad = ForwardDiff.gradient(f, _QAD_X)
            @test g_ad ≈ _qad_central(f, _QAD_X) rtol = 1e-6
        end
    end
end

nothing
