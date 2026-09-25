# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Propagation, compiled again here, while the package precompiles.
#
# AstroProp caches the same workload, and on its own that brings a first `propagate!` from 37 s to
# about 1 s. Loading AstroSolve undid it: SNOW brings in ReverseDiff, whose methods invalidate the
# code AstroProp cached, and the first propagation went back to 9-16 s (measured 2026-09-16 by
# loading each dependency on its own). Code compiled here is compiled with ReverseDiff already
# loaded, so it stays valid for a script that loads the whole stack, which is every script that
# says `using Epicycle`.

import AstroEpochs, AstroStates, AstroUniverse, AstroModels, AstroProp, AstroCallbacks

if ccall(:jl_generating_output, Cint, ()) == 1
    let
        epoch = AstroEpochs.Time("2020-01-01T00:00:00.000", AstroEpochs.UTC(), AstroEpochs.ISOT())
        states = (AstroStates.CartesianState([7000.0, 0.0, 1300.0, 0.0, 7.35, 1.0]),
                  AstroStates.KeplerianState(7000.0, 0.001, 0.9, 0.5, 0.0, 0.0))
        for integrator in (AstroProp.Vern9(), AstroProp.Tsit5()), state in states
            sc   = AstroModels.Spacecraft(state = state, time = epoch, name = "precompile")
            prop = AstroProp.OrbitPropagator(
                       AstroProp.ForceModel(AstroProp.PointMassGravity(AstroUniverse.earth, ())),
                       AstroProp.IntegratorConfig(integrator; dt = 60.0,
                                                  reltol = 1e-10, abstol = 1e-10))
            AstroProp.propagate!(prop, sc,
                AstroProp.StopAt(sc, AstroProp.PropDurationSeconds(), 120.0))
            AstroProp.propagate!(prop, sc,
                AstroProp.StopAt(sc, AstroProp.PropDurationDays(), 0.001))
            AstroProp.propagate!(prop, sc,
                AstroProp.StopAt(sc, AstroCallbacks.PosDotVel(), 0.0; direction = 1))
        end
    end
end
