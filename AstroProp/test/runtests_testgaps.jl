using Test
using LinearAlgebra
using SciMLBase

using AstroEpochs
using AstroStates
using AstroUniverse
using AstroModels: Spacecraft, to_posvel, set_posvel!
using AstroCallbacks: OrbitCalc, PosMag, PosDotVel
using AstroProp

# Helpers
make_sat() = Spacecraft(
    state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]),
    time  = Time("2015-09-21T12:23:12", TAI(), ISOT()),
)

forces_earth_only() = ForceModel(PointMassGravity(earth, ()))
integ_fast() = IntegratorConfig(Tsit5(); dt=10.0, reltol=1e-9, abstol=1e-9)

# 1) find_center branches: return nothing and conflict error
@testset "find_center branches" begin
    # A tiny dummy OrbitODE for exercising the 'no centers' branch
    struct _DummyForce <: AstroProp.OrbitODE end
    fm_none = ForceModel((_DummyForce(),))
    @test fm_none.center === nothing  # covers: return nothing

    # Conflicting centers error
    pm1 = PointMassGravity(earth, ())
    pm2 = PointMassGravity(moon,  ())
    @test_throws ErrorException ForceModel((pm1, pm2))  # covers: error("Multiple conflicting...")
end

# FR-FORCE-16: at most one force may add the central (monopole) term for a given body.
@testset "double-count central-term guard (FR-FORCE-16)" begin
    # Two point-mass forces both adding Earth's central term → double-count → error.
    @test_throws ErrorException ForceModel((PointMassGravity(earth, ()),
                                            PointMassGravity(earth, ())))

    # The real case: point-mass Earth + spherical-harmonic Earth → central term twice → error.
    sh = HarmonicGravity(earth; degree = 2, order = 0, model = Zonal())
    @test_throws ErrorException ForceModel((PointMassGravity(earth, ()), sh))

    # Valid: Sun & Moon as perturbers only (include_center = false) alongside the harmonic
    # Earth field — the central term is provided once, so this composes without error.
    third = PointMassGravity(earth, (moon, sun); include_center = false)
    fm = ForceModel((sh, third))
    @test fm.center === earth
end

# 2) propagate! direction handling: :infer default (no time condition) and the StopAt-direction /
#    propagation-direction mismatch guards in _compute_tf.
@testset "direction inference and mismatch guards" begin
    forces = forces_earth_only()
    integ  = integ_fast()

    # _infer_direction: state-only stop under :infer → defaults to :forward (covers that branch).
    sat = make_sat()
    sol = propagate!(DynSys(spacecraft=[sat], forces=forces), integ,
                     StopAt(sat, PosDotVel(), 0.0; direction=+1); direction=:infer)
    @test sol.retcode in (ReturnCode.Success, ReturnCode.Terminated)

    # _compute_tf: StopAt direction contradicts the inferred propagation direction → error.
    #   negative duration infers :backward, but stop_dir = +1 (increasing) contradicts it
    @test_throws ErrorException propagate!(DynSys(spacecraft=[make_sat()], forces=forces), integ,
        StopAt(make_sat(), PropDurationSeconds(), -3600.0; direction=+1); direction=:infer)
    #   positive duration infers :forward, but stop_dir = -1 (decreasing) contradicts it
    @test_throws ErrorException propagate!(DynSys(spacecraft=[make_sat()], forces=forces), integ,
        StopAt(make_sat(), PropDurationSeconds(),  3600.0; direction=-1); direction=:infer)
end

# 4) Legacy StopAt* backward-prop terminate branches (integrator.dt < 0)
@testset "Legacy StopAt* terminate on backward propagation" begin
    sat1 = make_sat()
    sat2 = make_sat()
    sat3 = make_sat()
    prop_forces = forces_earth_only()
    integ = integ_fast()

    sol1 = propagate!(DynSys(spacecraft=[sat1], forces=prop_forces), integ, StopAtSeconds(-10.0); direction = :backward)
    @test sol1.retcode in (ReturnCode.Success, ReturnCode.Terminated)

    sol2 = propagate!(DynSys(spacecraft=[sat2], forces=prop_forces), integ, StopAtDays(-1/86400); direction = :backward)
    @test sol2.retcode in (ReturnCode.Success, ReturnCode.Terminated)

    sol3 = propagate!(DynSys(spacecraft=[sat3], forces=prop_forces), integ, StopAtRadius(sat3, 7100.0); direction = :backward)
    @test sol3.retcode in (ReturnCode.Success, ReturnCode.Terminated)
end

nothing