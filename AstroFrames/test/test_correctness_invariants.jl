# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# =============================================================================
# Structural invariants of `axes_rotation`.
#
# These need no truth data and no SPICE kernels: they check properties every
# rotation must have regardless of which model produced it. A broken rotation
# fails orthonormality or the round trip immediately, which makes this the
# cheapest regression net available and the first thing to run after any
# refactor.
#
# Invariants checked: I-1 through I-3, I-10, I-11.
# Truth-data comparison against SPICE/GMAT is a separate concern — see
# `axes_rotation_harness.jl` and the internal rotation test matrix.
# =============================================================================

using AstroFrames
using AstroUniverse
using AstroEpochs
using LinearAlgebra
using Test

# Every body with an analytic IAU 2015 orientation, at three epochs spanning
# the model's useful range (J2000, recent, 50-year projection).
const _BODIES = (AstroUniverse.sun, AstroUniverse.mercury, AstroUniverse.venus,
                 AstroUniverse.mars, AstroUniverse.jupiter, AstroUniverse.saturn,
                 AstroUniverse.uranus, AstroUniverse.neptune, AstroUniverse.pluto)

const _EPOCHS = (2451545.0, 2458849.5, 2469854.0)

# Tolerances are named, not bare literals.
const _TOL_ORTHONORMAL = 1e-15   # I-2
const _TOL_ROUNDTRIP   = 1e-14   # I-3

"""
    _rotation_arcsec(R) -> Float64

Rotation angle of `R`, in arcseconds.

Uses the skew part rather than `acos((tr(R)-1)/2)`. Every angle in this file
is tiny — 23 mas for the frame bias, 0.3″ for polar motion — and near θ = 0
the trace formula evaluates `acos` at `1 − θ²/2`, where machine epsilon
propagates to a relative error of about `eps/θ²`. At θ ≈ 1.4e-6 rad that is
~1e-4, which is large enough to fail a legitimate check. The skew part gives
`sin θ` directly and stays well conditioned.
"""
function _rotation_arcsec(R::AbstractMatrix)
    v = 0.5 .* (R[3,2] - R[2,3], R[1,3] - R[3,1], R[2,1] - R[1,2])
    return rad2deg(asin(sqrt(sum(abs2, v)))) * 3600
end

@testset "axes_rotation invariants" begin

    @testset "I-1 block form" begin
        for b in _BODIES, jd in _EPOCHS
            M = axes_rotation(ICRF(), CelestialBodyFixed(b), jd)
            # Position must not depend on source velocity.
            @test all(iszero, M[1:3, 4:6])
            # Velocity rotates by the same R as position.
            @test M[1:3, 1:3] == M[4:6, 4:6]
        end
    end

    @testset "I-2 orthonormal, proper rotation" begin
        for b in _BODIES, jd in _EPOCHS
            R = axes_rotation(ICRF(), CelestialBodyFixed(b), jd)[1:3, 1:3]
            @test norm(R' * R - I) < _TOL_ORTHONORMAL
            # Proper (right-handed): a reflection would also be orthonormal.
            @test abs(det(R) - 1) < _TOL_ORTHONORMAL
        end
    end

    @testset "I-3 round trip" begin
        for b in _BODIES, jd in _EPOCHS
            ax = CelestialBodyFixed(b)
            @test norm(axes_rotation(ICRF(), ax, jd) *
                       axes_rotation(ax, ICRF(), jd) - I) < _TOL_ROUNDTRIP
        end
        # Static edges, epoch-independent.
        for (a, b) in ((MJ2000Eq(), MJ2000Ec()), (ICRF(), MJ2000Eq()))
            @test norm(axes_rotation(a, b, _EPOCHS[1]) * axes_rotation(b, a, _EPOCHS[1]) - I) < _TOL_ROUNDTRIP
        end
    end

    @testset "identity" begin
        for A in (ICRF(), GCRF(), ITRF(), MJ2000Eq(), MJ2000Ec(), MODEq(), TODEq())
            @test axes_rotation(A, A, 2451545.0) == I
        end
    end

    @testset "MJ2000Ec is a rotation about X by ε₀" begin
        # Analytic truth: the IAU 1976 / FK5 mean obliquity of
        # J2000, 23°26'21.448″. This is the one value that fixes the frame, so
        # it is checked directly rather than inferred from a truth row.
        M = axes_rotation(MJ2000Eq(), MJ2000Ec(), _EPOCHS[1])
        @test rad2deg(atan(M[2, 3], M[2, 2])) ≈ 23.4392911111 atol = 1e-10
        @test M[1, 1] == 1.0                 # X axis unchanged
        @test all(iszero, M[4:6, 1:3])       # static: Ṙ is exactly zero
    end

    @testset "ICRF ↔ MJ2000Eq is a 23 mas static frame bias" begin
        M = axes_rotation(ICRF(), MJ2000Eq(), 2451545.0)
        R = M[1:3, 1:3]

        # Epoch-independent by construction — a bias, not a precession.
        @test M == axes_rotation(ICRF(), MJ2000Eq(), 2469854.0)
        @test all(iszero, M[4:6, 1:3])          # Ṙ exactly zero
        @test norm(R' * R - I) < _TOL_ORTHONORMAL

        # Magnitude, checked against theory rather than against ourselves.
        # The IERS frame-bias angles are dα₀ = −14.6 mas, ξ₀ = −16.617 mas,
        # η₀ = −6.819 mas, so the total rotation is
        #     √(dα₀² + ξ₀² + η₀²) = 23.15 mas.
        # Worth checking as a value: a bias of the wrong size is invisible by
        # eye and reaches ≈16 km over 1 AU.
        θ_mas = _rotation_arcsec(R) * 1000
        @test θ_mas ≈ hypot(14.6, 16.617, 6.819) atol = 0.01

        # Orientation, not just magnitude: composing the bias with the static
        # obliquity must shift the recovered ε by η₀ ≈ −6.82 mas, one of the
        # three standard frame-bias angles. A sign-flipped or mis-axed bias
        # has the right magnitude and fails this.
        M_ec    = axes_rotation(MJ2000Eq(), MJ2000Ec(), _EPOCHS[1]) * M
        ε_mas   = rad2deg(atan(M_ec[2, 3], M_ec[2, 2])) * 3600 * 1000
        ε₀_mas  = 23.4392911111 * 3600 * 1000
        @test (ε_mas - ε₀_mas) ≈ 6.82 atol = 0.05
    end

    @testset "MJ2000Eq ↔ MODEq reproduces the precession rate" begin
        # Physics, not just structure. General precession is ~50.3"/yr and is
        # zero at the reference epoch by definition, so the rotation angle is
        # a direct check on the model rather than on the plumbing.
        R2000 = axes_rotation(MJ2000Eq(), MODEq(), 2451545.0)[1:3, 1:3]
        @test _rotation_arcsec(R2000) < 1e-6

        # 20 years on: 50.3"/yr × 20 ≈ 1006". 1% is loose enough that the
        # rate's own nonlinearity is not being tested, tight enough that a
        # wrong time scale or a factor error fails.
        R2020 = axes_rotation(MJ2000Eq(), MODEq(), 2458849.5)[1:3, 1:3]
        @test _rotation_arcsec(R2020) ≈ 1006.0 rtol = 0.01

        for jd in _EPOCHS
            M = axes_rotation(MJ2000Eq(), MODEq(), jd)
            @test norm(M[1:3,1:3]' * M[1:3,1:3] - I) < _TOL_ORTHONORMAL
            @test norm(M * axes_rotation(MODEq(), MJ2000Eq(), jd) - I) < _TOL_ROUNDTRIP
            @test all(iszero, M[4:6, 1:3])   # Ṙ neglected for precession
        end
    end

    @testset "MODEq ↔ TODEq nutation is bounded, not secular" begin
        # Nutation is periodic — longitude up to ~17.2", obliquity ~9.2" — so
        # the total rotation stays in a narrow band forever. Precession, by
        # contrast, grows ~50"/yr. Checking the band therefore catches the
        # error that matters here: wiring a secular model in where a periodic
        # one belongs, which orthonormality would happily accept.
        angles = Float64[]
        for jd in _EPOCHS
            M = axes_rotation(MODEq(), TODEq(), jd)
            R = M[1:3, 1:3]
            push!(angles, _rotation_arcsec(R))
            @test norm(R' * R - I) < _TOL_ORTHONORMAL
            @test norm(M * axes_rotation(TODEq(), MODEq(), jd) - I) < _TOL_ROUNDTRIP
            @test all(iszero, M[4:6, 1:3])   # Ṙ neglected for nutation
        end
        @test all(a -> 10.0 < a < 20.0, angles)
        # Over 50 years the spread stays small; precession would move ~2500".
        @test maximum(angles) - minimum(angles) < 10.0
    end

    @testset "TODEq ↔ PEF carries Earth's spin rate" begin
        # The first edge with a real Ṙ. Every check here is on the rate,
        # because that is what distinguishes this edge from the ones whose
        # rate is neglected — and getting it wrong leaves position exactly
        # right while velocity is wrong by up to ω·r.
        for jd in (2451545.0, 2458849.5)
            M = axes_rotation(TODEq(), PEF(), jd)
            R, Ṙ = M[1:3, 1:3], M[4:6, 1:3]

            @test norm(R' * R - I) < _TOL_ORTHONORMAL
            @test norm(M * axes_rotation(PEF(), TODEq(), jd) - I) < _TOL_ROUNDTRIP

            # Ṙ must NOT be zero here. A regression to _rotation_no_rate would
            # pass every structural check above and fail only this one.
            @test !all(iszero, Ṙ)

            # Ṙ = −[ω]ₓR, so −ṘRᵀ must be the skew-symmetric [ω]ₓ. This is a
            # structural proof of the rate's form, independent of its size.
            Ω = -Ṙ * R'
            @test norm(Ω + Ω') < 1e-18

            # And its size is Earth's spin rate. The LOD correction moves it
            # by ~1e-8 relative, so the tolerance sits above that and far
            # below anything a modelling error would produce.
            ω = Ω[2, 1]
            @test ω ≈ 7.292115146706979e-5 rtol = 1e-7

            # Same fact stated as a period: the sidereal day.
            @test 2π / ω / 3600 ≈ 23.9345 atol = 1e-3
        end
    end

    @testset "PEF ↔ ITRF is polar motion" begin
        for jd in (2451545.0, 2458849.5)
            M = axes_rotation(PEF(), ITRF(), jd)
            R = M[1:3, 1:3]
            @test norm(R' * R - I) < _TOL_ORTHONORMAL
            @test norm(M * axes_rotation(ITRF(), PEF(), jd) - I) < _TOL_ROUNDTRIP
            @test all(iszero, M[4:6, 1:3])   # polar-motion rate neglected

            # The rotation angle is the magnitude of the pole offset, so it
            # can be checked against the table directly rather than trusted.
            # Also bounds it: the pole wanders a few tenths of an arcsecond,
            # never degrees, so a units error (arcsec taken as radians) fails
            # here by five orders of magnitude.
            tbl   = AstroUniverse.eop(AstroUniverse.FK5())
            jdutc = epoch_utc(AstroFrames._scales(jd))
            expected = hypot(tbl.x(jdutc), tbl.y(jdutc))
            @test _rotation_arcsec(R) ≈ expected rtol = 1e-6
            @test expected < 1.0
        end
    end

    @testset "the FK5 chain composes: ICRF → ITRF" begin
        # Five edges end to end — bias, precession, nutation, spin, polar
        # motion. The point is that composition preserves the properties each
        # edge has individually; a sign or ordering error in any one of them
        # shows up here even when that edge passes in isolation.
        jd = 2458849.5
        C = axes_rotation(PEF(), ITRF(), jd) *
            axes_rotation(TODEq(), PEF(), jd) *
            axes_rotation(MODEq(), TODEq(), jd) *
            axes_rotation(MJ2000Eq(), MODEq(), jd) *
            axes_rotation(ICRF(), MJ2000Eq(), jd)
        R, Ṙ = C[1:3, 1:3], C[4:6, 1:3]

        @test norm(R' * R - I) < 1e-14
        @test C[1:3, 1:3] ≈ C[4:6, 4:6]

        # Earth's spin must survive composition: the only rotating edge in the
        # chain is the GAST one, so the composed rate is still Earth's.
        Ω = -Ṙ * R'
        @test norm(Ω + Ω') < 1e-18
        @test Ω[2, 1] ≈ 7.292115146706979e-5 rtol = 1e-7
    end

    @testset "ICRF ↔ GCRF is exactly identity" begin
        # By definition, not by approximation: the GCRS is kinematically
        # non-rotating with respect to the ICRS. Tested as exact equality —
        # a tolerance here would invite someone to slip a small rotation in.
        for jd in _EPOCHS
            @test axes_rotation(ICRF(), GCRF(), jd) == I
            @test axes_rotation(GCRF(), ICRF(), jd) == I
        end
    end

    @testset "the IAU-2006 chain: GCRF → CIRS → TIRS → ITRF" begin
        jd = 2458849.5
        for (a, b, spins) in ((GCRF(), CIRS(), false),
                              (CIRS(), TIRS(), true),
                              (TIRS(), ITRF(), false))
            M = axes_rotation(a, b, jd)
            R = M[1:3, 1:3]
            @test norm(R' * R - I) < _TOL_ORTHONORMAL
            @test norm(M * axes_rotation(b, a, jd) - I) < _TOL_ROUNDTRIP
            # Only the ERA edge rotates. Getting this backwards is the
            # velocity-only failure the two rotation helpers exist to prevent.
            @test all(iszero, M[4:6, 1:3]) == !spins
        end

        # The ERA edge must carry Earth's rate, exactly as its FK5 counterpart
        # does — same physics, different formulation.
        M = axes_rotation(CIRS(), TIRS(), jd)
        Ω = -M[4:6, 1:3] * M[1:3, 1:3]'
        @test norm(Ω + Ω') < 1e-18
        @test Ω[2, 1] ≈ 7.292115146706979e-5 rtol = 1e-7

        # Polar motion is the same physical quantity as on the FK5 side, so
        # the two polar-motion edges must agree on its magnitude.
        tbl   = AstroUniverse.eop(AstroUniverse.IAU2006())
        jdutc = epoch_utc(AstroFrames._scales(jd))
        @test _rotation_arcsec(axes_rotation(TIRS(), ITRF(), jd)[1:3, 1:3]) ≈
              hypot(tbl.x(jdutc), tbl.y(jdutc)) rtol = 1e-6
    end

    @testset "the two theories agree to the frame-bias scale" begin
        # Cross-family: FK5 and IAU-2006 are different models, so their ITRF
        # orientations differ — but by how much is a meaningful check. The
        # FK5 chain refers its precession to MJ2000Eq, which is itself the
        # frame bias away from ICRF, so the two chains should land about one
        # frame bias apart and not further.
        jd = 2458849.5
        C6 = axes_rotation(TIRS(), ITRF(), jd) * axes_rotation(CIRS(), TIRS(), jd) *
             axes_rotation(GCRF(), CIRS(), jd) * axes_rotation(ICRF(), GCRF(), jd)
        C5 = axes_rotation(PEF(), ITRF(), jd) * axes_rotation(TODEq(), PEF(), jd) *
             axes_rotation(MODEq(), TODEq(), jd) * axes_rotation(MJ2000Eq(), MODEq(), jd) *
             axes_rotation(ICRF(), MJ2000Eq(), jd)
        Δ_mas = _rotation_arcsec(C6[1:3,1:3] * C5[1:3,1:3]') * 1000
        @test 20.0 < Δ_mas < 30.0
    end

    @testset "obliquity edges reach the ecliptic of date" begin
        # At J2000 the mean obliquity must be the same constant MJ2000Ec uses.
        # That ties this edge to an independently-fixed value rather than to
        # its own output, and would catch a wrong series or a wrong epoch.
        ε_2000 = _rotation_arcsec(axes_rotation(MODEq(), MODEc(), 2451545.0)[1:3,1:3]) / 3600
        @test ε_2000 ≈ 23.4392911111 atol = 1e-6

        # And it drifts: the obliquity decreases ~46.8"/century, so 20 years
        # from J2000 is about −9.4". A static value would pass the check above
        # and fail this one.
        ε_2020 = _rotation_arcsec(axes_rotation(MODEq(), MODEc(), 2458849.5)[1:3,1:3]) / 3600
        @test (ε_2020 - ε_2000) * 3600 ≈ -9.36 atol = 0.2

        for jd in (2451545.0, 2458849.5)
            for (a, b) in ((MODEq(), MODEc()), (TODEq(), TODEc()))
                M = axes_rotation(a, b, jd)
                R = M[1:3, 1:3]
                @test norm(R' * R - I) < _TOL_ORTHONORMAL
                @test norm(M * axes_rotation(b, a, jd) - I) < _TOL_ROUNDTRIP
                @test all(iszero, M[4:6, 1:3])     # obliquity drift neglected
                @test R[1, 1] ≈ 1.0                # rotation about X only
            end

            # True obliquity is the mean plus the nutation in obliquity, which
            # is bounded by ~9.2". Larger would mean the wrong term; zero would
            # mean TODEc had silently reused the mean value.
            mean_ε = _rotation_arcsec(axes_rotation(MODEq(), MODEc(), jd)[1:3,1:3])
            true_ε = _rotation_arcsec(axes_rotation(TODEq(), TODEc(), jd)[1:3,1:3])
            @test 0.0 < abs(true_ε - mean_ε) < 9.5
        end
    end

    @testset "TEME and the equation of the equinoxes" begin
        jd = 2458849.5
        M  = axes_rotation(TODEq(), TEME(), jd)
        R  = M[1:3, 1:3]

        @test norm(R' * R - I) < _TOL_ORTHONORMAL
        @test norm(M * axes_rotation(TEME(), TODEq(), jd) - I) < _TOL_ROUNDTRIP
        @test all(iszero, M[4:6, 1:3])      # no spin between them

        # The equation of the equinoxes is Δψ·cos ε, bounded near 16″.
        @test 0.0 < _rotation_arcsec(R) < 20.0

        # Independent cross-check, and the point of attaching TEME to TODEq.
        # STB also offers a direct TEME → PEF rotation built on GMST. Since
        # GAST = GMST + the equation of the equinoxes, our routed path — which
        # composes the equation of the equinoxes with GAST — must agree with
        # it. Two different formulations, so this catches an error in either
        # the TEME edge or the GAST edge, which self-consistency cannot.
        STB = Base.require(Base.PkgId(
            Base.UUID("6b019ec1-7a1e-4f04-96c7-a9db1ca5514d"),
            "SatelliteToolboxTransformations"))

        original = AstroUniverse.frame_theory()
        try
            AstroUniverse.set_frame_theory!(AstroUniverse.FK5())
            e      = AstroFrames._scales(jd)
            jd_utc = epoch_utc(e)
            jd_ut1 = jd_utc + AstroUniverse.eop(AstroUniverse.FK5()).Δut1_utc(jd_utc) / 86_400

            ours   = axes_rotation(TEME(), PEF(), jd)[1:3, 1:3]
            theirs = STB.r_teme_to_pef(STB.DCM, jd_ut1)
            @test maximum(abs.(ours .- theirs)) < 1e-14
        finally
            AstroUniverse.set_frame_theory!(original)
        end
    end

    @testset "edge_theory declares each edge's theory" begin
        # The mechanism that makes passing the wrong EOP table unreachable:
        # an edge names its theory, and the theory selects the table.
        @test edge_theory(MJ2000Eq(), MODEq()) === AstroUniverse.FK5()
        @test edge_theory(MODEq(), TODEq())    === AstroUniverse.FK5()
        @test edge_theory(TODEq(), MODEq())    === AstroUniverse.FK5()
        @test edge_theory(TODEq(), PEF())      === AstroUniverse.FK5()
        @test edge_theory(PEF(), ITRF())       === AstroUniverse.FK5()
        @test edge_theory(MODEq(), MODEc())    === AstroUniverse.FK5()
        @test edge_theory(TODEq(), TODEc())    === AstroUniverse.FK5()
        @test edge_theory(TODEq(), TEME())     === AstroUniverse.FK5()
        @test edge_theory(GCRF(), CIRS())      === AstroUniverse.IAU2006()
        @test edge_theory(CIRS(), TIRS())      === AstroUniverse.IAU2006()
        @test edge_theory(TIRS(), ITRF())      === AstroUniverse.IAU2006()
        @test edge_theory(ICRF(), GCRF())      === nothing

        # Edges needing no Earth orientation data say so explicitly.
        @test edge_theory(ICRF(), MJ2000Eq())  === nothing
        @test edge_theory(MJ2000Eq(), MJ2000Ec()) === nothing
        @test edge_theory(ICRF(), MoonPA())    === nothing
    end

    @testset "both EOP series are available at once" begin
        # The defect this replaced: one active table meant every edge of the
        # other theory was unserviceable. They carry different fields, so one
        # cannot stand in for the other.
        fk5  = AstroUniverse.eop(AstroUniverse.FK5())
        iau  = AstroUniverse.eop(AstroUniverse.IAU2006())
        @test typeof(fk5) != typeof(iau)
        @test :δΔψ in propertynames(fk5)   # nutation corrections
        @test :δx  in propertynames(iau)   # CIP offsets
        # Switching the theory must not evict either.
        AstroUniverse.set_frame_theory!(AstroUniverse.FK5())
        @test AstroUniverse.eop(AstroUniverse.IAU2006()) === iau
        AstroUniverse.set_frame_theory!(AstroUniverse.IAU2006())
        @test AstroUniverse.eop(AstroUniverse.FK5()) === fk5
    end

    @testset "both epoch forms are supported and agree" begin
        # Two public forms, deliberately: `Time` carries its scale, a bare
        # Julian date is terser at a REPL. Both are supported, so both are
        # tested — most of this suite calls the scalar form, which would
        # otherwise leave the `Time` path almost unexercised.
        jd = 2458849.5
        t  = Time(jd, 0.0, :tdb, :jd)
        for (a, b) in ((ICRF(), MJ2000Eq()), (MJ2000Eq(), MODEq()),
                       (MODEq(), TODEq()),   (TODEq(), PEF()),
                       (PEF(), ITRF()),      (GCRF(), CIRS()),
                       (CIRS(), TIRS()),     (ICRF(), ITRF()))
            @test axes_rotation(a, b, jd) == axes_rotation(a, b, t)
        end
        # The scalar form means TDB, and nothing but the docstring says so —
        # which is why the Time form exists.
        @test axes_rotation(ICRF(), ITRF(), jd) ==
              axes_rotation(ICRF(), ITRF(), Time(jd, 0.0, :tdb, :jd))
    end

    @testset "Correctness: edge_theory is declared both ways round" begin
        # A route may traverse an edge in either direction, and the theory
        # filter reads `edge_theory` for the direction it is walking. A
        # forward-only declaration would silently drop reverse edges out of
        # their own theory's subgraph.
        for (a, b) in ((MJ2000Eq(), MODEq()), (MODEq(), TODEq()), (TODEq(), PEF()),
                       (PEF(), ITRF()), (MODEq(), MODEc()), (TODEq(), TODEc()),
                       (TODEq(), TEME()))
            @test edge_theory(a, b) === FK5()
            @test edge_theory(b, a) === FK5()
        end
        for (a, b) in ((GCRF(), CIRS()), (CIRS(), TIRS()), (TIRS(), ITRF()))
            @test edge_theory(a, b) === IAU2006()
            @test edge_theory(b, a) === IAU2006()
        end
    end

    @testset "Correctness: the epoch accessors read the scales they name" begin
        # The public way a user's own frame reads the epoch. Reading the wrong
        # scale is a ~70 s error, invisible in a rotation matrix.
        t = Time(_EPOCHS[2], 0.0, :tdb, :jd)
        e = AstroFrames._scales(t)

        @test epoch_tdb(e) == t.tdb.jd
        @test epoch_tt(e)  == t.tt.jd
        @test epoch_utc(e) == t.utc.jd

        # They are genuinely different scales, not three names for one number.
        @test epoch_tdb(e) != epoch_tt(e)
        @test epoch_tt(e)  != epoch_utc(e)
        # TT − TAI is 32.184 s and TAI − UTC was 37 s in 2020, so TT − UTC is
        # 69.184 s. Bounded rather than exact: leap seconds come from a table.
        @test 60 < (epoch_tt(e) - epoch_utc(e)) * 86_400 < 80
    end

    @testset "Correctness: the Earth Rotation Angle keeps its low bits" begin
        # Found by comparison against ERFA: forming ERA as
        #   2π(0.7790572732640 + 1.00273781191135448 d)
        # lets the bracket reach ~9200 by 2025, where a Float64's spacing is
        # already 2.4 µas of angle, and the whole turns are then discarded by
        # the reduction after they have cost the precision. That lost 0.7 µas
        # by 2025 and grew with epoch.
        #
        # The property, stated without needing ERFA: the split form and the
        # direct form must agree, and where they do not it is the direct form
        # that is wrong. Checked against a high-precision evaluation of the
        # same expression, which is the arbiter neither float form can be.
        for jd_ut1 in (2451545.0, 2455197.5, 2458849.5, 2460676.5, 2469854.0)
            θ = AstroFrames._earth_rotation_angle(jd_ut1)
            @test 0 <= θ < 2π

            d_big = big(jd_ut1) - big(2451545)
            θ_big = mod(2 * big(π) * (parse(BigFloat, "0.7790572732640") +
                                      parse(BigFloat, "1.00273781191135448") * d_big),
                        2 * big(π))
            err_μas = abs(rad2deg(Float64(θ - θ_big))) * 3600 * 1e6
            @test err_μas < 0.05        # measured: under 0.01 across this span

            # And the direct form really is worse, so this is not a no-op.
            θ_direct = mod2pi(2π * (0.7790572732640 + 1.00273781191135448 *
                                    (jd_ut1 - 2451545.0)))
            direct_μas = abs(rad2deg(Float64(θ_direct - θ_big))) * 3600 * 1e6
            @test direct_μas >= err_μas
        end
    end

    @testset "Correctness: mean obliquity is the IERS arcsecond expression" begin
        # Found by comparison against ERFA: taking ε_A from `nutation_fk5` gave
        # Vallado's decimal-degree rounding, 23.439291° = 84381.4476″, which is
        # 0.4 mas below the 84381.448″ the IERS conventions publish — and the
        # linear terms differ by another 0.12 mas/century. It showed up only in
        # the ecliptic edges, where ε appears as a single rotation with nothing
        # to cancel against.
        e0 = AstroFrames._scales(Time(2451545.0, 0.0, :tt, :jd))
        @test rad2deg(AstroFrames._mean_obliquity(e0)) * 3600 ≈ 84_381.448 atol = 1e-6

        # The rotation actually applied by the edge carries the same value.
        R = axes_rotation(MODEq(), MODEc(), Time(2451545.0, 0.0, :tt, :jd))[1:3, 1:3]
        @test rad2deg(atan(R[2, 3], R[2, 2])) * 3600 ≈ 84_381.448 atol = 1e-6

        # And the drift matches the published rate over a century.
        e1 = AstroFrames._scales(Time(2451545.0 + 36_525.0, 0.0, :tt, :jd))
        drift = (AstroFrames._mean_obliquity(e1) - AstroFrames._mean_obliquity(e0))
        @test rad2deg(drift) * 3600 ≈ -46.8150 - 0.00059 + 0.001813 atol = 1e-6
    end

    @testset "Robustness: an integer Julian date is accepted" begin
        # The signature says `Real`, and an Int is a Real. It used to die as
        # `InexactError: Int64(0.5)` several calls down, because `Time` splits
        # the date in two and the split is fractional. Found by a generated
        # test file emitting `2451545` without a decimal point.
        for jd in (2451545, 2458849)
            @test axes_rotation(ICRF(), MJ2000Eq(), jd) ≈
                  axes_rotation(ICRF(), MJ2000Eq(), float(jd))
        end
        @test axes_rotation(ICRF(), ITRF(), 2458849) isa AbstractMatrix
    end

    @testset "epoch carries its scale" begin
        # The point of taking a `Time` rather than a bare number: the scale
        # travels with the value. Passing the same numeric JD tagged UTC
        # instead of TDB must give a different answer, because it is a
        # different instant — 64 s apart at J2000.
        jd    = 2458849.5
        M_tdb = axes_rotation(MJ2000Eq(), MODEq(), Time(jd, 0.0, :tdb, :jd))
        M_utc = axes_rotation(MJ2000Eq(), MODEq(), Time(jd, 0.0, :utc, :jd))
        @test M_tdb != M_utc

        # And the bare-number entry point means TDB.
        @test axes_rotation(MJ2000Eq(), MODEq(), jd) == M_tdb
    end

    @testset "I-10 unsupported pair fails loudly" begin
        # Naming both axes is the contract; a bare MethodError would not.
        # `Inertial` has no edges at all, so no route exists under any theory.
        # (MODEq → CIRS used to sit here; it now warns and routes instead —
        # see the warn-once test below.)
        err = try
            axes_rotation(Inertial(), ITRF(), 2451545.0)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("Inertial", err.msg) && occursin("ITRF", err.msg)
    end

    @testset "FR-FRAME-4: routing outside the active theory warns once" begin
        # MODEq is an FK5 node; under IAU2006 it is unreachable within the
        # active theory. Per the ruling this proceeds and warns rather than
        # failing — a mixed chain is a modelling choice, not an error.
        original = AstroUniverse.frame_theory()
        try
            AstroUniverse.set_frame_theory!(AstroUniverse.IAU2006())
            empty!(AstroFrames._OUT_OF_THEORY_WARNED)

            @test_logs (:warn,) axes_rotation(MODEq(), CIRS(), 2451545.0)

            # Once per distinct pair: a transform inside a propagation loop
            # must not warn on every step.
            @test_logs axes_rotation(MODEq(), CIRS(), 2451545.0)
            @test_logs axes_rotation(MODEq(), CIRS(), 2458849.5)

            # And it produced a usable answer, not just a warning.
            M = axes_rotation(MODEq(), CIRS(), 2451545.0)
            @test norm(M[1:3,1:3]' * M[1:3,1:3] - I) < _TOL_ORTHONORMAL
        finally
            AstroUniverse.set_frame_theory!(original)
            empty!(AstroFrames._OUT_OF_THEORY_WARNED)
        end
    end

    @testset "I-11 origin coupling" begin
        @test CoordinateSystem(AstroUniverse.earth, ITRF()) isa CoordinateSystem
        @test_throws ArgumentError CoordinateSystem(AstroUniverse.sun, ITRF())
        @test_throws ArgumentError CoordinateSystem(AstroUniverse.sun, TEME())
        @test_throws ArgumentError CoordinateSystem(AstroUniverse.earth, MoonPA())
        # Body-fixed axes must match their origin.
        @test_throws ArgumentError CoordinateSystem(AstroUniverse.venus,
                                                    CelestialBodyFixed(AstroUniverse.mars))
    end
end
