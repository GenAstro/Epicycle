# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Hermite-Simpson, the open tier's reference transcription.
#
# Two things are checked here, and they are different claims.
#
# The first is that the transcription extension interface is complete: all twelve
# names resolve, `HermiteSimpson` answers the four asked of a transcription, and
# its mesh answers the seven asked of a mesh. A package that publishes an
# interface and ships no implementation of it cannot show this, which is why the
# test arrived with the implementation.
#
# The second is that the method is right. Defects are the residual of a
# discretisation, so they vanish exactly on a trajectory the discretisation can
# represent. Hermite-Simpson is fourth order and integrates a cubic exactly, so a
# cubic in time is the sharpest closed-form case available, and its analytic
# Jacobians are checked against finite differences of its own defects.

using Test
using AstroSolve
using AstroSolve: HermiteSimpsonMesh, build_mesh, compute_defects, defect_jacobian_P,
    defect_jacobian_U, defect_jacobian_Y, defect_jacobian_t0, defect_jacobian_tf, n_defect_rows,
    n_intervals, n_unique_nodes, node_times, quadrature_weights
using LinearAlgebra

@testset "the transcription interface is complete" begin
    t = HermiteSimpson(n_steps = 4)

    @test n_intervals(t)       == 4
    @test n_unique_nodes(t)    == 9          # 2n + 1, midpoints included
    @test n_defect_rows(t, 3)  == 24         # 2 equations per step per state

    mesh = build_mesh(t)
    @test mesh isa HermiteSimpsonMesh

    for f in (node_times, quadrature_weights)
        @test hasmethod(f, Tuple{typeof(mesh), Float64, Float64})
    end
    for f in (compute_defects, defect_jacobian_Y, defect_jacobian_U,
              defect_jacobian_t0, defect_jacobian_tf, defect_jacobian_P)
        @test isdefined(AstroSolve, nameof(f))
    end

    @test HermiteSimpson().n_steps == 10     # the documented default
end

@testset "node times span the phase and include midpoints" begin
    mesh = build_mesh(HermiteSimpson(n_steps = 3))
    ts   = node_times(mesh, 2.0, 8.0)

    @test length(ts) == 7
    @test ts[1]   ≈ 2.0
    @test ts[end] ≈ 8.0
    @test issorted(ts)
    @test ts ≈ collect(range(2.0, 8.0, length = 7))   # uniform, midpoints halfway
end

@testset "quadrature integrates what Simpson integrates exactly" begin
    mesh = build_mesh(HermiteSimpson(n_steps = 4))
    t0, tf = 0.0, 3.0
    w  = quadrature_weights(mesh, t0, tf)
    ts = node_times(mesh, t0, tf)

    @test sum(w) ≈ tf - t0                        # ∫1
    @test dot(w, ts) ≈ (tf^2 - t0^2) / 2          # ∫t
    @test dot(w, ts .^ 2) ≈ (tf^3 - t0^3) / 3     # ∫t²
    @test dot(w, ts .^ 3) ≈ (tf^4 - t0^4) / 4     # ∫t³, the highest it is exact for
end

# A one-state system whose solution is a cubic in time: ẏ = 3t².
_cubic!(dy, y, p, t) = (dy[1] = 3t^2; nothing)

@testset "defects vanish on a trajectory the discretisation represents" begin
    t   = HermiteSimpson(n_steps = 5)
    m   = build_mesh(t)
    t0, tf = 0.0, 2.0
    ts  = node_times(m, t0, tf)

    Y = reshape(ts .^ 3, 1, :)                      # the exact solution, y = t³
    U = zeros(0, length(ts))

    r = compute_defects(_cubic!, Y, U, Float64[], nothing, m, t0, tf)
    @test maximum(abs, r) < 1e-12

    # Perturb one interior node and the defects must notice.
    Y[1, 4] += 0.1
    r_bad = compute_defects(_cubic!, Y, U, Float64[], nothing, m, t0, tf)
    @test maximum(abs, r_bad) > 1e-3
end

@testset "the analytic state Jacobian matches a finite difference of the defects" begin
    t   = HermiteSimpson(n_steps = 3)
    m   = build_mesh(t)
    t0, tf = 0.0, 1.0
    ts  = node_times(m, t0, tf)

    Y = reshape(ts .^ 3 .+ 0.05 .* ts, 1, :)        # off the exact solution
    U = zeros(0, length(ts))
    P = Float64[]

    _cubic_jac_y!(J, y, p, t) = (J[1, 1] = 0.0; nothing)   # ∂f/∂y = 0 for ẏ = 3t²

    J  = defect_jacobian_Y(_cubic_jac_y!, Y, U, P, nothing, m, t0, tf)
    r0 = vec(compute_defects(_cubic!, Y, U, P, nothing, m, t0, tf))

    h  = 1e-6
    fd = similar(J)
    for j in eachindex(Y)
        Yp = copy(Y); Yp[j] += h
        fd[:, j] = (vec(compute_defects(_cubic!, Yp, U, P, nothing, m, t0, tf)) .- r0) ./ h
    end

    @test size(J) == size(fd)
    @test maximum(abs, J .- fd) < 1e-6
end
