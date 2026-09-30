# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# =============================================================================
# Body orientation models.
#
# The interface a user writes their own model against. The integration case —
# a user-defined body reached from every other frame — lives in AstroFrames;
# this file covers the pieces in isolation.
# =============================================================================

using AstroUniverse
using ForwardDiff
using LinearAlgebra
using Test

const _JD_OR = 2458849.5

@testset "pole_axes_rotation" begin

    @testset "block form and orthonormality" begin
        M = pole_axes_rotation(0.3, 0.7, 1.1, 0.0, 0.0, 1e-4)
        @test size(M) == (6, 6)

        R, Ṙ = M[1:3, 1:3], M[4:6, 1:3]
        @test M[1:3, 4:6] ≈ zeros(3, 3)        # upper-right block is zero
        @test M[4:6, 4:6] ≈ R                  # lower-right repeats R
        @test norm(R' * R - I) < 1e-15
        @test det(R) ≈ 1
    end

    @testset "a still body has no rate" begin
        M = pole_axes_rotation(0.3, 0.7, 1.1, 0.0, 0.0, 0.0)
        @test M[4:6, 1:3] ≈ zeros(3, 3)
    end

    @testset "spin about the pole turns at the rate given" begin
        # The rate vector must lie along the pole and have the magnitude asked
        # for. A sign or axis error here is invisible in R alone.
        Ẇ = 7.0e-5
        α₀, δ₀ = 0.3, 0.7
        M = pole_axes_rotation(α₀, δ₀, 1.1, 0.0, 0.0, Ẇ)
        R, Ṙ = M[1:3, 1:3], M[4:6, 1:3]

        Ω = Ṙ * R'
        @test norm(Ω + Ω') < 1e-18                       # skew, so it is a rate
        ω = [Ω[3,2], Ω[1,3], Ω[2,1]]
        @test norm(ω) ≈ Ẇ rtol = 1e-12

        # `ṘRᵀ` is the rate in BODY coordinates, so it lies along the body's
        # own z. Checking it against the ICRF pole instead gives sin(δ₀) — the
        # first draft of this test did exactly that.
        @test abs(dot(ω / norm(ω), [0, 0, 1])) ≈ 1 rtol = 1e-12

        # Rotated back out, it points at the pole, which is the statement that
        # actually matters: the body spins about the axis it was given.
        pole = [cos(δ₀)cos(α₀), cos(δ₀)sin(α₀), sin(δ₀)]
        @test abs(dot(R' * (ω / norm(ω)), pole)) ≈ 1 rtol = 1e-12
    end

    @testset "the third row is the pole" begin
        α₀, δ₀ = 1.9, -0.4
        R = pole_axes_rotation(α₀, δ₀, 2.2, 0.0, 0.0, 0.0)[1:3, 1:3]
        @test R[3, :] ≈ [cos(δ₀)cos(α₀), cos(δ₀)sin(α₀), sin(δ₀)] atol = 1e-14
    end

    @testset "mixed argument types promote" begin
        @test pole_axes_rotation(0, 0.7, 1.1, 0, 0, 1e-4) isa AbstractMatrix
    end
end

@testset "the shipped polynomial model" begin

    @testset "every body it covers produces a proper rotation" begin
        for n in (10, 199, 299, 499, 599, 699, 799, 899, 999)
            M = body_axes_rotation(IAU2015(), n, _JD_OR)
            R = M[1:3, 1:3]
            @test norm(R' * R - I) < 1e-14
            @test det(R) ≈ 1
        end
    end

    @testset "reading orientation by NAIF ID matches reading it by body" begin
        # The NAIF-ID method is what the model interface calls; the body method
        # is what existing code calls. They must not drift apart.
        for (b, n) in ((sun, 10), (mars, 499), (jupiter, 599), (pluto, 999))
            @test iau2015_orientation(b, _JD_OR) === iau2015_orientation(n, _JD_OR)
        end
    end

    @testset "it has nothing to estimate, which is the honest answer" begin
        # Its coefficients are fixed. Solving for a pole means registering a
        # model that carries the pole as a parameter.
        @test isempty(orientation_parameters(IAU2015()))
    end

    @testset "the old name still works, as a type" begin
        # It was an exported type in 0.3, so code may test or dispatch on it.
        @test IauPolynomialOrientation === IAU2015
        @test IauPolynomialOrientation() === IAU2015()
        @test IAU2015() isa IauPolynomialOrientation
        _old(::IauPolynomialOrientation) = :dispatched
        @test _old(IAU2015()) === :dispatched
    end
end

@testset "the IAU 1991 model" begin
    @testset "a proper rotation, with the rate block its derivative" begin
        M = body_axes_rotation(IAU1991(), 499, _JD_OR)
        R, Ṙ = M[1:3, 1:3], M[4:6, 1:3]
        @test norm(R' * R - I) < 1e-14
        @test det(R) ≈ 1
        # ForwardDiff through the date: a finite difference on a date near
        # 2.46e6 resolves only about 40 µs.
        dR = ForwardDiff.derivative(d -> Matrix(body_axes_rotation(IAU1991(), 499, d)[1:3, 1:3]),
                                    _JD_OR)
        @test norm(Ṙ - dR / 86_400) / norm(Ṙ) < 1e-12
    end

    @testset "Mars only" begin
        e = try; body_axes_rotation(IAU1991(), 599, _JD_OR); nothing; catch e; e; end
        @test e isa ArgumentError
        @test occursin("Mars", e.msg)
    end

    @testset "0.236° from IAU 2015 at J2000, nearly all of it prime meridian" begin
        A = body_axes_rotation(IAU1991(), 499, 2451545.0)[1:3, 1:3]
        B = body_axes_rotation(IAU2015(), 499, 2451545.0)[1:3, 1:3]
        @test rad2deg(acos((tr(A * B') - 1) / 2)) ≈ 0.236 atol = 0.001
        @test rad2deg(acos(clamp(dot(A[3, :], B[3, :]), -1, 1))) < 0.001      # poles
    end
end

@testset "which model a body uses" begin

    struct _TestSpin{T<:Real} <: AbstractOrientationModel
        pole_ra::T
        pole_dec::T
        spin_rate::T
    end
    AstroUniverse.body_axes_rotation(m::_TestSpin, naifid, jd) =
        pole_axes_rotation(m.pole_ra, m.pole_dec, m.spin_rate * (jd - 2451545.0),
                           zero(m.pole_ra), zero(m.pole_ra), m.spin_rate / 86_400)

    @testset "defaults are the published models, unasked" begin
        for b in (sun, mercury, venus, mars, jupiter, saturn, uranus, neptune, pluto)
            @test orientation_model(b) isa IAU2015
            @test has_orientation_model(b)
        end
        @test orientation_model(moon) === LunarPA()
        @test has_orientation_model(moon) && has_orientation_model(earth)
    end

    @testset "Earth's orientation is the frame theory, one setting under two names" begin
        original = frame_theory()
        try
            set_frame_theory!(IAU2006())
            @test orientation_model(earth) === IAU2006()
            @test set_orientation!(earth, FK5()) === FK5()
            @test frame_theory() === FK5()
            @test orientation_model(earth) === FK5()
            @test orientation_model(399) === FK5()
        finally
            set_frame_theory!(original)
        end

        # Anything else in Earth's place would silently replace the ITRF chain.
        e = try; set_orientation!(earth, IAU2015()); nothing; catch e; e; end
        @test e isa ArgumentError
        @test occursin("frame theory", e.msg)
        @test frame_theory() === original

        # The chain needs UT1 and TT, so a TDB date alone points to AstroFrames.
        e = try; body_axes_rotation(IAU2006(), 399, _JD_OR); nothing; catch e; e; end
        @test e isa ArgumentError
        @test occursin("body_fixed_rotation", e.msg)
        @test FK5() isa AbstractOrientationModel
    end

    @testset "a body nobody has heard of names the fix" begin
        rock = CelestialBody("Test Rock", 1e-9, 1.0, 0.0, 2999999)
        @test !has_orientation_model(rock)
        e = try; orientation_model(rock); nothing; catch e; e; end
        @test e isa ArgumentError
        @test occursin("set_orientation!", e.msg)
        @test occursin("2999999", e.msg)
    end

    @testset "registering one, and taking it back" begin
        rock  = CelestialBody("Test Rock 2", 1e-9, 1.0, 0.0, 2999998)
        model = _TestSpin(0.3, 0.7, 1000.0)

        set_orientation!(rock, model)
        @test orientation_model(rock) === model
        @test has_orientation_model(rock)

        # Keyed on NAIF ID, so a second reference to the same body agrees.
        same = CelestialBody("Test Rock 2 again", 1e-9, 1.0, 0.0, 2999998)
        @test orientation_model(same) === model

        # And it overrides a shipped default rather than being ignored.
        original = orientation_model(mars)
        try
            set_orientation!(mars, model)
            @test orientation_model(mars) === model
        finally
            set_orientation!(mars, original)
        end
        @test orientation_model(mars) isa IAU2015
    end
end

@testset "orientation parameters" begin

    struct _P3{T<:Real} <: AbstractOrientationModel
        pole_ra::T
        pole_dec::T
        spin_rate::T
    end

    @testset "read off the model's fields, with no extra code" begin
        m = _P3(1.0, 2.0, 3.0)
        p = orientation_parameters(m)
        @test keys(p) == (:pole_ra, :pole_dec, :spin_rate)
        @test values(p) == (1.0, 2.0, 3.0)
    end

    @testset "writing replaces only what is named" begin
        m = _P3(1.0, 2.0, 3.0)
        m2 = set_orientation_parameters(m, (; pole_dec = 20.0))
        @test orientation_parameters(m2) == (pole_ra = 1.0, pole_dec = 20.0, spin_rate = 3.0)
        @test orientation_parameters(m) == (pole_ra = 1.0, pole_dec = 2.0, spin_rate = 3.0)
    end

    @testset "an unknown name is refused, and the message names the real ones" begin
        e = try
            set_orientation_parameters(_P3(1.0, 2.0, 3.0), (; albedo = 0.1))
            nothing
        catch e; e; end
        @test e isa ArgumentError
        @test occursin("albedo", e.msg)
        @test occursin("pole_ra", e.msg)
    end

    @testset "writing a wider type widens the model" begin
        # An estimated parameter arrives carrying derivative information, so
        # the model has to accept a type it was not built with. Stood in for
        # here by Float64 -> Rational, which has the same shape.
        m2 = set_orientation_parameters(_P3(1.0, 2.0, 3.0), (; pole_ra = 1//2))
        @test m2 isa _P3
        @test orientation_parameters(m2).pole_ra == 1//2
    end

    @testset "through the body" begin
        rock = CelestialBody("Test Rock 3", 1e-9, 1.0, 0.0, 2999997)
        set_orientation!(rock, _P3(1.0, 2.0, 3.0))

        @test orientation_parameters(rock).pole_ra == 1.0
        set_orientation_parameters!(rock, (; pole_ra = 9.0))
        @test orientation_parameters(rock).pole_ra == 9.0
        @test orientation_parameters(rock).spin_rate == 3.0     # untouched
    end
end

@testset "the SPICE-backed model" begin

    # `IAU_<BODY>` frames come from a text PCK, which is in the kernel manifest
    # but is not one of the defaults — nothing in the shipped models needs it,
    # since Epicycle carries its own IAU polynomials. So fetch it here, and put
    # the kernel pool back as it was afterwards.
    #
    # A failure to fetch skips the checks rather than failing them: what would
    # be broken then is the network, not the model.
    _pck = "pck00011.tpc"

    _pck_loaded = try
        body_axes_rotation(SpiceOrientation("IAU_MARS"), 499, _JD_OR)
        true
    catch
        false
    end

    _pck_fetched = false
    if !_pck_loaded
        try
            download_spice_kernel(_pck)
            load_spice_kernel(_pck)
            _pck_loaded = _pck_fetched = true
        catch e
            @info "Could not fetch $(_pck); skipping the two SpiceOrientation " *
                  "numerical checks. ($(sprint(showerror, e)))"
        end
    end

    if _pck_loaded
        @testset "it works where a kernel defines the frame" begin
            R = body_axes_rotation(SpiceOrientation("IAU_MARS"), 499, _JD_OR)[1:3, 1:3]
            @test norm(R' * R - I) < 1e-12
            @test det(R) ≈ 1
        end

        @testset "it agrees with the polynomial model to microarcseconds" begin
            # Two independent implementations of the same IAU constants: our
            # polynomials, and NAIF's in pck00011.tpc. Agreement is at the
            # level of floating-point roundoff, so the tolerance is set to
            # catch a changed constant or a dropped term rather than to
            # accommodate a real difference between the models.
            #
            # Measured at this epoch, in arcseconds:
            #   Mercury 6.4e-8   Venus 5.7e-12   Mars 9.0e-7   Jupiter 5.2e-7
            #   Saturn  1.5e-6   Uranus 1.8e-7   Neptune 1.7e-6  Pluto 2.3e-8
            #
            # This compares against pck00011 specifically. An older text PCK
            # carries the 2009 constants and would legitimately differ by more.
            for (naifid, frame) in ((199, "IAU_MERCURY"), (299, "IAU_VENUS"),
                                    (499, "IAU_MARS"),    (599, "IAU_JUPITER"),
                                    (699, "IAU_SATURN"),  (799, "IAU_URANUS"),
                                    (899, "IAU_NEPTUNE"), (999, "IAU_PLUTO"))
                A = body_axes_rotation(SpiceOrientation(frame), naifid, _JD_OR)[1:3, 1:3]
                B = body_axes_rotation(IAU2015(), naifid, _JD_OR)[1:3, 1:3]
                D = A * B'
                v = 0.5 .* (D[3,2] - D[2,3], D[1,3] - D[3,1], D[2,1] - D[1,2])
                @test rad2deg(asin(sqrt(sum(abs2, v)))) * 3600 < 1e-3   # arcseconds
            end
        end
    end

    @testset "it carries nothing to estimate, and says why" begin
        @test isempty(orientation_parameters(SpiceOrientation("IAU_MARS")))
        e = try
            set_orientation_parameters(SpiceOrientation("IAU_MARS"), (; pole_ra = 1.0))
            nothing
        catch e; e; end
        @test e isa ArgumentError
        @test occursin("kernel", e.msg)
    end

    _pck_fetched && unload_spice_kernel(_pck)
end

@testset "IAU1991 against NAIF's PCK of the 1991 report" begin
    # preliminary_pck00004.tpc carries the 1991 constants, so SPICE's IAU_MARS
    # with it loaded is an independent evaluation of the same model. Loaded
    # alone and unloaded after: its Mars constants would otherwise replace
    # pck00011's. A failure to fetch skips, as above.
    name = "preliminary_pck00004.tpc"
    url  = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/a_old_versions/" * name
    fetched = try
        download_spice_kernel(name, url)
        load_spice_kernel(name)
        true
    catch e
        @info "Could not fetch $(name); skipping the IAU1991 truth check. ($(sprint(showerror, e)))"
        false
    end
    # Recorded as skipped, so an offline run shows it in the summary rather than passing quietly.
    fetched || @test_skip "IAU1991 against preliminary_pck00004.tpc"
    if fetched
        try
            for jd in (2451545.0, _JD_OR, 2462502.5)          # 2000, the test epoch, 2030
                A = body_axes_rotation(SpiceOrientation("IAU_MARS"), 499, jd)
                B = body_axes_rotation(IAU1991(), 499, jd)
                D = A[1:3, 1:3] * B[1:3, 1:3]'
                v = 0.5 .* (D[3,2] - D[2,3], D[1,3] - D[3,1], D[2,1] - D[1,2])
                @test rad2deg(asin(sqrt(sum(abs2, v)))) * 3600 < 1e-3   # arcseconds
                @test norm(A[4:6, 1:3] - B[4:6, 1:3]) / norm(B[4:6, 1:3]) < 1e-9
            end
        finally
            unload_spice_kernel(name)
        end
    end
end

@testset "the lunar models" begin
    M = body_axes_rotation(LunarPA(), 301, _JD_OR)
    @test norm(M[1:3, 1:3]' * M[1:3, 1:3] - I) < 1e-12
    # They read the versioned DE440 frames, which is what the aliases name in
    # the frame kernel loaded at startup.
    @test M == body_axes_rotation(SpiceOrientation("MOON_PA_DE440"), 301, _JD_OR)
    @test body_axes_rotation(LunarME(), 301, _JD_OR) ==
          body_axes_rotation(SpiceOrientation("MOON_ME_DE440_ME421"), 301, _JD_OR)
    # Principal axes and mean-Earth axes: about 875 m apart on the surface.
    A, B = M[1:3, 1:3], body_axes_rotation(LunarME(), 301, _JD_OR)[1:3, 1:3]
    @test 800 < acos((tr(A * B') - 1) / 2) * 1737.4e3 < 950
    @test isempty(orientation_parameters(LunarPA()))
    e = try; body_axes_rotation(LunarPA(), 499, _JD_OR); nothing; catch e; e; end
    @test e isa ArgumentError && occursin("Moon", e.msg)

    # A dual date cannot pass into SPICE, and the error says that rather than blaming the kernels.
    for m in (LunarPA(), LunarME(), SpiceOrientation("MOON_PA_DE440"))
        e = try
            ForwardDiff.derivative(d -> body_axes_rotation(m, 301, d)[1, 1], _JD_OR)
            nothing
        catch e; e; end
        @test e isa ArgumentError
        @test occursin("not differentiable", e.msg) && !occursin("load_spice_kernel", e.msg)
    end

    # Outside the kernel's coverage, SPICE's failure is reported with what to do about it.
    e = try; body_axes_rotation(LunarPA(), 301, 1.0e5); nothing; catch e; e; end
    @test e isa ArgumentError
    @test occursin("1550 to 2650", e.msg) && occursin("SPICE reported", e.msg)

    # The mean-Earth axes can be the Moon's default instead, for maps and surface work.
    try
        @test set_orientation!(moon, LunarME()) === LunarME()
        @test orientation_model(moon) === LunarME()
    finally
        set_orientation!(moon, LunarPA())
    end
    @test orientation_model(moon) === LunarPA()
end

@testset "a shipped model is refused for a body it cannot orient" begin
    for (body, model, says) in ((moon, IAU2015(), "planets"), (mars, LunarPA(), "Moon only"),
                                (jupiter, IAU1991(), "Mars only"), (mars, FK5(), "Earth"))
        before = orientation_model(body)
        e = try; set_orientation!(body, model); nothing; catch e; e; end
        @test e isa ArgumentError && occursin(says, e.msg)
        @test orientation_model(body) === before                  # nothing was registered
    end
    try
        @test set_orientation!(mars, IAU1991()) === IAU1991()
    finally
        set_orientation!(mars, IAU2015())
    end
end

@testset "setting no parameters writes nothing" begin
    # For the Earth that would re-set the frame theory it had just read.
    @test set_orientation_parameters!(earth, (;)) === frame_theory()
    @test set_orientation_parameters!(mars, (;)) === orientation_model(mars)
end
