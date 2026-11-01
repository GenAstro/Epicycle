# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# An Earth orientation correction installed in AstroUniverse reaches the Earth frame chain.
#
# Truth is the geometry of each edge. A polar motion correction δ turns the polar motion edge by
# exactly δ; a UT1 correction turns the spin edge by δ times the rotation rate per UT1 second;
# a length-of-day correction scales the spin edge's rotation rate. Removing the corrections gives
# back the rotation that was there before, bit for bit.

using AstroFrames
using AstroUniverse
using AstroEpochs
using LinearAlgebra
using Test

struct _EopShift <: AbstractEopCorrection
    x::Float64            # arcsec
    ut1_utc::Float64      # s
    lod::Float64          # ms
end
AstroUniverse.eop_correction(m::_EopShift, jd_utc, jd_tt) =
    (x = m.x, y = 0.0, ut1_utc = m.ut1_utc, lod = m.lod)

# The angle of a small rotation, from its skew part.
_angle(M) = norm((M[3, 2] - M[2, 3], M[1, 3] - M[3, 1], M[2, 1] - M[1, 2])) / 2
_R(M) = M[1:3, 1:3]
_Rdot(M) = M[4:6, 1:3]

const _ARCSEC = deg2rad(1 / 3600)
const _ERA_RATE = 2π * 1.00273781191135448 / 86_400      # rad per second of UT1

@testset "an EOP correction reaches the Earth frame chain" begin
    before = eop_corrections()
    epoch = Time("2020-03-15T06:30:00", UTC(), ISOT())
    try
        set_eop_corrections!([])
        pairs = ((TIRS(), ITRF()), (PEF(), ITRF()), (CIRS(), TIRS()), (TODEq(), PEF()),
                 (ICRF(), ITRF()))
        plain = Dict(pair => axes_rotation(pair..., epoch) for pair in pairs)

        # Polar motion: the edge turns by the correction.
        set_eop_corrections!([_EopShift(1.0, 0.0, 0.0)])
        for pair in ((TIRS(), ITRF()), (PEF(), ITRF()))
            M = _R(axes_rotation(pair..., epoch)) * _R(plain[pair])'
            @test _angle(M) ≈ _ARCSEC rtol = 1e-9
        end

        # UT1: the spin edge turns by the correction times the rotation rate.
        set_eop_corrections!([_EopShift(0.0, 0.5, 0.0)])
        for pair in ((CIRS(), TIRS()), (TODEq(), PEF()))
            M = _R(axes_rotation(pair..., epoch)) * _R(plain[pair])'
            @test _angle(M) ≈ 0.5 * _ERA_RATE rtol = 1e-6
        end

        # Length of day: the spin rate, and so Ṙ, scales by (1 − lod / 86400 s).
        set_eop_corrections!([_EopShift(0.0, 0.0, 86.4)])
        for pair in ((CIRS(), TIRS()), (TODEq(), PEF()))
            ratio = norm(_Rdot(axes_rotation(pair..., epoch))) / norm(_Rdot(plain[pair]))
            @test ratio ≈ 1 - 1e-6 rtol = 1e-9
        end

        # The whole chain moves, and comes back.
        set_eop_corrections!([_EopShift(0.1, 0.1, 1.0)])
        @test axes_rotation(ICRF(), ITRF(), epoch) != plain[(ICRF(), ITRF())]
        set_eop_corrections!([])
        @test axes_rotation(ICRF(), ITRF(), epoch) == plain[(ICRF(), ITRF())]
    finally
        set_eop_corrections!(before)
    end
end

nothing
