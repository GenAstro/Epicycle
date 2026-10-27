# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

using AstroFrames, AstroUniverse, AstroEpochs, AstroStates, ForwardDiff, LinearAlgebra, Test

# AF-13. Truth: central differences of the complete operator, not its approximate
# velocity-rate block. Steps keep date rounding below the derivative being checked.
@testset "Earth epoch derivatives [AF-13]" begin
    saved_theory, saved_cip = frame_theory(), cip_interpolation()
    try
        for theory in (FK5(), IAU2006()), interpolate in (false, true)
            set_frame_theory!(theory)
            set_cip_interpolation!(interpolate)
            pairs = theory isa FK5 ?
                ((ICRF(), MODEq()), (MODEq(), TODEq()), (TODEq(), PEF()),
                 (PEF(), ITRF()), (MODEq(), MODEc()), (TODEq(), TODEc()),
                 (TODEq(), TEME()), (ICRF(), ITRF())) :
                ((GCRF(), CIRS()), (CIRS(), TIRS()), (TIRS(), ITRF()), (ICRF(), ITRF()))
            for pair in pairs, reverse in (false, true)
                source, target = reverse ? (pair[2], pair[1]) : pair
                f(x) = vec(Matrix(axes_rotation(source, target,
                    Time(2458849., 0.371 + x[1], TDB(), JD()))))
                ad = ForwardDiff.jacobian(f, [0.])[:, 1]
                # A spin edge changes over a day; other edges change much more slowly.
                # Away from daily EOP knots, a slow edge uses a 0.001-day step. Mean
                # obliquity needs 0.01 day to resolve its smaller secular rate.
                obliquity = (source isa MODEq && target isa MODEc) ||
                            (source isa MODEc && target isa MODEq)
                h = norm(ad) > 1 ? 1e-5 : obliquity ? 0.01 : 0.001
                fd = (f([h]) - f([-h])) / (2h)
                @test all(isfinite, ad)
                @test norm(ad - fd) <= 3e-6 * max(norm(fd), 1e-12)
            end
        end
    finally
        set_frame_theory!(saved_theory)
        set_cip_interpolation!(saved_cip)
    end
end

@testset "Orbit-relative reference Jacobians [AF-13]" begin
    for axes in (RIC(), LVLH(), VNB()), reverse in (false, true),
        state in ([7000., 400., -200., -0.2, 7.5, 0.8],
                  [0., -9000., 2000., 5., 0.4, -1.])
        source, target = reverse ? (axes, ICRF()) : (ICRF(), axes)
        q = vcat(state, [-0.007, 0.0002, 0.0005])
        f(x) = vec(Matrix(axes_rotation(source, target, 0.,
            (; reference_state = x[1:6], reference_accel = x[7:9]))))
        ad = ForwardDiff.jacobian(f, q)
        @test all(isfinite, ad)
        for i in eachindex(q)
            # Scale a position perturbation to the orbital radius, including a zero
            # component; a micrometre step there loses the derivative to subtraction.
            scale = i <= 3 ? norm(q[1:3]) : i <= 6 ? norm(q[4:6]) : 1.
            h = cbrt(eps(Float64)) * max(abs(q[i]), scale)
            hi, lo = copy(q), copy(q)
            hi[i] += h
            lo[i] -= h
            fd = (f(hi) - f(lo)) / (2h)
            @test norm(ad[:, i] - fd) <= 3e-7 * max(norm(fd), 1e-12)
        end
    end
end

# The conversion docstring says state derivatives survive an origin change at a fixed epoch: the
# ephemeris offset is constant there. Truth: central differences of the same conversion.
@testset "State Jacobian through an origin change [AF-13]" begin
    t = Time(2458849.5, 0.0, TDB(), JD())
    f(x) = to_vector(CartesianState(Coordinate(x, CoordinateSystem(earth, ICRF()), t),
                                    CoordinateSystem(moon, MoonME())))
    x = [7000., 300., 100., 0.1, 7.5, 0.03]
    ad = ForwardDiff.jacobian(f, x)
    h = 1e-3
    fd = hcat([(f(x .+ h .* (1:6 .== i)) - f(x .- h .* (1:6 .== i))) / 2h for i in 1:6]...)
    @test all(isfinite, ad)
    @test norm(ad - fd) <= 1e-7 * norm(fd)
end
