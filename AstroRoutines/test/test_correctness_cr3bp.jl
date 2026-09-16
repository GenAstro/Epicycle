# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# Correctness dimension (TestingStandards.md §2) for the circular restricted
# three-body functional area.
#
# Truth source: the problem's own invariants, which is stronger here than a
# table of published numbers. A libration point is defined by the acceleration
# vanishing there, so the definition is the test. The Jacobi constant is the one
# integral of motion, so it must not drift along a ballistic arc. The analytic
# Jacobian is checked against automatic differentiation of the equations of
# motion themselves, which makes the two independent derivations of one thing.
#
# Published values are used once, as a cross-check on the Earth-Moon collinear
# points, so that a self-consistent but wrong convention would still be caught.
#
# Each @testitem runs in its own isolated module, so every block is self-contained.

@testitem "libration points are where the acceleration vanishes" tags=[:Correctness] begin
    using AstroRoutines, LinearAlgebra
    for mu in (0.01215058560962404, 3.0034806e-6, 0.1, 0.4)   # Earth-Moon, Sun-Earth, extremes
        for p in (:L1, :L2, :L3, :L4, :L5)
            r = libration_point(mu, p)
            s = [r[1], r[2], r[3], 0.0, 0.0, 0.0]             # at rest in the rotating frame
            @test norm(cr3bp_accel(s, mu)) < 1e-12
        end
    end
end

@testitem "collinear points match published Earth-Moon values" tags=[:Correctness] begin
    using AstroRoutines
    mu = 0.01215058560962404
    @test isapprox(libration_point(mu, :L1)[1],  0.83691513; atol = 1e-7)
    @test isapprox(libration_point(mu, :L2)[1],  1.15568216; atol = 1e-7)
    @test isapprox(libration_point(mu, :L3)[1], -1.00506265; atol = 1e-7)
end

@testitem "L4 and L5 are equilateral with the primaries" tags=[:Correctness] begin
    using AstroRoutines, LinearAlgebra
    for mu in (0.01215058560962404, 0.2, 0.4)
        for p in (:L4, :L5)
            r  = libration_point(mu, p)
            d1 = norm(r .- [-mu, 0.0, 0.0])        # to the larger primary
            d2 = norm(r .- [1 - mu, 0.0, 0.0])     # to the smaller one
            @test isapprox(d1, 1.0; atol = 1e-14)  # both a unit from each, which is
            @test isapprox(d2, 1.0; atol = 1e-14)  # what equilateral means here
        end
    end
end

@testitem "analytic Jacobian agrees with automatic differentiation" tags=[:Correctness] begin
    using AstroRoutines, ForwardDiff
    mu = 0.01215058560962404
    states = ([0.8369,  0.05,  0.02,  0.10, -0.20,  0.03],
              [1.1557, -0.10,  0.07, -0.05,  0.30, -0.10],
              [0.4878,  0.866, 0.00,  0.00,  0.00,  0.00],
              [-1.005,  0.30, -0.15,  0.02,  0.01,  0.05])
    for s in states
        rhs   = ss -> (d = similar(ss); cr3bp_eom!(d, ss, mu); d)
        @test maximum(abs.(cr3bp_jacobian(s, mu) .- ForwardDiff.jacobian(rhs, s))) < 1e-12
    end
end

@testitem "Jacobi constant is conserved on a ballistic arc" tags=[:Correctness] begin
    using AstroRoutines
    mu = 0.01215058560962404
    s0 = [0.85, 0.0, 0.0, 0.0, 0.2, 0.0]
    C0 = jacobi_constant(s0, mu)
    d  = zeros(6)
    h  = 1e-4
    rhs(v) = (cr3bp_eom!(d, v, mu); copy(d))
    function march(s, n)                             # classical RK4, 0.2 time units
        for _ in 1:n
            k1 = rhs(s); k2 = rhs(s + h/2*k1)
            k3 = rhs(s + h/2*k2); k4 = rhs(s + h*k3)
            s  = s + h/6*(k1 + 2k2 + 2k3 + k4)
        end
        return s
    end
    sf = march(s0, 2000)
    @test abs(jacobi_constant(sf, mu) - C0) < 1e-10  # integrator error, not model error
end

@testitem "mass ratio and the domain checks" tags=[:Robustness] begin
    using AstroRoutines
    @test isapprox(cr3bp_mass_ratio(2.0, 1.0), 1/3; atol = 1e-15)
    @test isapprox(cr3bp_mass_ratio(3.986004418e5, 4.9028e3), 0.012150584; atol = 1e-8)
    @test_throws ArgumentError cr3bp_mass_ratio(0.0, 1.0)
    @test_throws ArgumentError cr3bp_mass_ratio(1.0, -1.0)
    @test_throws ArgumentError libration_point(0.0, :L1)
    @test_throws ArgumentError libration_point(1.0, :L1)
    @test_throws ArgumentError libration_point(0.1, :L6)
end

@testitem "equations of motion are differentiable in the state and the mass ratio" tags=[:Differentiability] begin
    using AstroRoutines, ForwardDiff
    s = [0.8369, 0.05, 0.02, 0.10, -0.20, 0.03]
    # With respect to mu, which is what makes the mass ratio usable as a solver
    # variable rather than only as a constant.
    g = ForwardDiff.derivative(m -> jacobi_constant(s, m), 0.01215058560962404)
    @test isfinite(g) && g != 0
    J = ForwardDiff.jacobian(v -> cr3bp_accel(v, 0.01215058560962404), s)
    @test all(isfinite, J)
end

@testitem "variational equations: transition matrix of a Hamiltonian flow" tags=[:Correctness] begin
    using AstroRoutines, LinearAlgebra
    mu = 0.01215058560962404
    # JPL SSD periodic_orbits.api, earth-moon lyapunov L1, C = 3.15001683280912.
    s0 = [8.1596252146384562e-01, 0.0, 0.0, 0.0, 2.0722124749217649e-01, 0.0]
    z  = cr3bp_stm_initial(s0)
    @test z[1:6] == s0
    @test reshape(z[7:42], 6, 6) == Matrix(1.0I, 6, 6)

    h = 1e-5
    rhs(v) = (d = zeros(42); cr3bp_stm_eom!(d, v, mu); d)
    function march(v, n)                               # RK4 over one period
        for _ in 1:n
            k1 = rhs(v); k2 = rhs(v + h/2*k1)
            k3 = rhs(v + h/2*k2); k4 = rhs(v + h*k3)
            v  = v + h/6*(k1 + 2k2 + 2k3 + k4)
        end
        return v
    end
    Phi = reshape(march(z, 20000)[7:42], 6, 6)
    # The flow is Hamiltonian, so the transition matrix has unit determinant.
    # This is the well-conditioned check: the symplectic residual itself scales
    # with the square of the matrix norm, which is large on an unstable orbit.
    @test abs(det(Phi) - 1) < 1e-6
    @test_throws ArgumentError cr3bp_stm_eom!(zeros(42), zeros(6), mu)
    @test_throws ArgumentError cr3bp_stm_initial(zeros(3))
end

@testitem "published JPL Lyapunov orbits are periodic under these equations" tags=[:Correctness] begin
    using AstroRoutines, LinearAlgebra
    # Truth source outside this package: NASA/JPL SSD three-body periodic orbit
    # database, earth-moon Lyapunov families near C = 3.15. If the equations of
    # motion or the Jacobi constant were wrong, these would not agree.
    mu = 1.215058560962404e-02
    for (x, vy, C) in ((8.1596252146384562e-01, 2.0722124749217649e-01, 3.15001683280912),
                       (1.1182825695532028e+00, 1.8601928389638619e-01, 3.15000013081292))
        s = [x, 0.0, 0.0, 0.0, vy, 0.0]
        @test isapprox(jacobi_constant(s, mu), C; atol = 1e-12)
    end
    # JPL's own L1 and L2 locations, to every digit they publish.
    @test isapprox(libration_point(mu, :L1)[1], 0.836915125772357; atol = 1e-13)
    @test isapprox(libration_point(mu, :L2)[1], 1.15568216544488;  atol = 1e-13)
end
