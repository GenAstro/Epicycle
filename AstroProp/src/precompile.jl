# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# A short propagation run while the package precompiles, so the code it compiles is cached with
# the package rather than compiled on a user's first call.
#
# Measured on 2026-09-16 before this existed: the first `propagate!` in a session took 37 s, all
# of it compilation, and the second took a millisecond. The ODE solvers are precompiled by their
# own packages, but only for their own types; the right-hand side here closes over a spacecraft,
# a force model and an epoch, and that combination is compiled nowhere until someone propagates.
#
# The workload runs only while generating the package image, and only exercises the two-body
# force model, which needs no ephemeris kernels. It covers the two integrators AstroProp loads,
# both state representations a spacecraft is commonly built with, and the stopping conditions a
# first script reaches for: a fixed duration and an apsis.
#
# This uses the compiler's own "am I generating output" flag rather than PrecompileTools, which
# would add a dependency for the same effect.

if ccall(:jl_generating_output, Cint, ()) == 1
    let
        epoch = Time("2020-01-01T00:00:00.000", UTC(), ISOT())
        states = (CartesianState([7000.0, 0.0, 1300.0, 0.0, 7.35, 1.0]),
                  KeplerianState(7000.0, 0.001, 0.9, 0.5, 0.0, 0.0))
        for integrator in (Vern9(), Tsit5()), state in states
            sc   = Spacecraft(state = state, time = epoch, name = "precompile")
            prop = OrbitPropagator(ForceModel(PointMassGravity(earth, ())),
                                   IntegratorConfig(integrator; dt = 60.0,
                                                    reltol = 1e-10, abstol = 1e-10))
            propagate!(prop, sc, StopAt(sc, PropDurationSeconds(), 120.0))
            propagate!(prop, sc, StopAt(sc, PropDurationDays(), 0.001))
            propagate!(prop, sc, StopAt(sc, AstroCallbacks.PosDotVel(), 0.0; direction = 1))
        end
    end
end
