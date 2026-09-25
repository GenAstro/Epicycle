# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0
#
# Benchmark: the cost of a stopping condition written in the quantity form against the same stop
# written as a calculation tag. The stop is evaluated at every accepted step and again while the
# root is bracketed, so its cost is paid on every propagation that uses one.
#
# Three stops over the same ten-day point-mass propagation, each run to the tenth apoapsis:
#   tag             StopAt(sat, PosDotVel(), 0.0; direction = -1)
#   quantity        StopAt(position_dot_velocity, sat; equals = 0.0, direction = -1)
#   quantity+frame  StopAt(position_dot_velocity, sat, EarthMJ2000Ec; equals = 0.0, direction = -1)
# The last reads the quantity through a frame transform, which is the case the quantity form
# exists for.
#
# Run manually (not wired into runtests.jl):
#   using Pkg; Pkg.activate("<environment>")
#   include(joinpath(pkgdir(AstroProp), "test", "bench_quantity_stop.jl"))

using AstroProp
using AstroModels, AstroStates, AstroEpochs, AstroFrames
using AstroCallbacks: PosDotVel, position_dot_velocity
using AstroUniverse: earth
using BenchmarkTools

make_sat() = Spacecraft(state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]),
                        time  = Time("2020-01-01T00:00:00", UTC(), ISOT()),
                        save_history = false)

const _PROP = OrbitPropagator(ForceModel(PointMassGravity(earth, ())),
                              IntegratorConfig(Vern9(); dt = 60.0, reltol = 1e-10, abstol = 1e-10))

function ten_apoapses(stop_for)
    sat = make_sat()
    for _ in 1:10
        propagate!(_PROP, sat, stop_for(sat))
    end
    return sat
end

stops = (
    "tag"            => sat -> StopAt(sat, PosDotVel(), 0.0; direction = -1),
    "quantity"       => sat -> StopAt(position_dot_velocity, sat; equals = 0.0, direction = -1),
    "quantity+frame" => sat -> StopAt(position_dot_velocity, sat, EarthMJ2000Ec;
                                      equals = 0.0, direction = -1),
)

for (label, stop_for) in stops
    ten_apoapses(stop_for)                                  # compile
    t = @belapsed ten_apoapses($stop_for) samples = 5 evals = 1
    println(rpad(label, 16), lpad(round(t * 1e3; digits = 1), 8), " ms")
end
