# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# ForceContext: what the forces in one evaluation share.
#
# The context exists for speed and must not change a result. The truth is the direct computation
# each accessor stands in for: the time scales from the Time, body_fixed_rotation, translate and
# translate_state. The tests check that a value from the context is that value exactly, that a
# second request returns the first value rather than recomputing it, that a context for another
# epoch is ignored, that a Dual epoch leaves the context unused, and that a propagation gives the
# same state with and without one.

using Test
using AstroProp
using AstroModels, AstroStates, AstroEpochs
using AstroUniverse: earth, moon, sun, translate, translate_state, orientation_model
using AstroFrames: body_fixed_rotation
using OrdinaryDiffEqVerner: Vern9
using ForwardDiff

const _APC = AstroProp

@testset "ForceContext" begin
    t  = Time("2020-10-20T12:00:00", UTC(), ISOT()).tt
    t2 = t + 30.0 / 86400.0
    om = orientation_model(earth)
    ctx = _APC.ForceContext()
    params = (context = ctx,)
    _APC._reset!(ctx, t)

    @testset "values are the direct computation" begin
        e = _APC.force_epoch(params, t)
        @test e.tdb == t.tdb.jd && e.tt == t.tt.jd && e.utc == t.utc.jd
        @test _APC.force_rotation(params, om, 399, t) == body_fixed_rotation(om, 399, t)
        @test _APC.force_position(params, earth, sun, t) == translate(earth, sun, t.tdb.jd)
        @test _APC.force_position(params, earth, moon, t) == translate(earth, moon, t.tdb.jd)
        @test _APC.force_state(params, sun, earth, t) == translate_state(sun, earth, t.tdb.jd)
    end

    @testset "computed once per evaluation" begin
        @test length(ctx.rotations) == 1 && length(ctx.positions) == 2 && length(ctx.states) == 1
        # A second request is served from the context: change the stored value and it comes back.
        ctx.positions[1] = ctx.positions[1] .+ 1.0
        @test _APC.force_position(params, earth, sun, t) == translate(earth, sun, t.tdb.jd) .+ 1.0
        @test length(ctx.positions) == 2
        # A reset empties it.
        _APC._reset!(ctx, t)
        @test isempty(ctx.positions) && ctx.epoch === nothing
        @test _APC.force_position(params, earth, sun, t) == translate(earth, sun, t.tdb.jd)
    end

    @testset "another epoch, or none, is computed directly" begin
        _APC._reset!(ctx, t)
        @test _APC.force_position(params, earth, sun, t2) == translate(earth, sun, t2.tdb.jd)
        @test isempty(ctx.positions)
        @test _APC.force_rotation(nothing, om, 399, t2) == body_fixed_rotation(om, 399, t2)
        @test _APC.force_position([], earth, moon, t2) == translate(earth, moon, t2.tdb.jd)
        # A Dual epoch, as when the time itself is differentiated, leaves the context unused.
        td = Time(ForwardDiff.Dual(t.tt.jd, 1.0), 0.0, TT(), JD())
        _APC._reset!(ctx, td)
        @test ctx.time === nothing
    end

    @testset "a propagation is unchanged" begin
        sc() = Spacecraft(; state = CartesianState([6878.137, 0.0, 0.0, 0.0, 4.71754, 5.99820]),
                          time = Time("2020-10-20T12:00:00", UTC(), ISOT()), mass = 1000.0,
                          drag = SphericalDrag(c_d = 2.2, drag_area = 10.0),
                          srp = SphericalSRP(c_r = 1.8, srp_area = 10.0))
        forces = ForceModel(HarmonicGravity(earth; degree = 4, order = 0),
                            PointMassGravity(earth, (moon, sun); include_center = false),
                            AtmosphericDrag(earth; model = Exponential()),
                            SolarRadiationPressure(earth))
        integ = IntegratorConfig(Vern9(); reltol = 1e-12, abstol = 1e-12, dt = 60.0)
        s1 = sc()
        y1 = propagate!(OrbitPropagator(forces, integ), s1,
                        StopAt(s1, PropDurationSeconds(), 7200.0)).u[end]
        # The same right-hand side with no context: each force computes its own values. It ends
        # steps at the shadow boundaries as propagate! does, so the steps are the same.
        s2 = sc()
        y0 = collect(to_posvel(s2))
        start = s2.time.tt
        f!(dy, y, _p, τ) = _APC._eval_all!(forces, start + τ / 86400.0, y, dy, s2)
        y2 = _APC.solve(_APC.ODEProblem(f!, y0, (0.0, 7200.0)), Vern9();
                        callback = _APC._kink_callback(forces, start, (1:6,)),
                        reltol = 1e-12, abstol = 1e-12, dt = 60.0).u[end]
        @test y1 == y2
    end
end
