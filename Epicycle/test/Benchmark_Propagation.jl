using BenchmarkTools
using Profile
using ProfileView
using StatProfilerHTML

# ==========================
# Execution
# ==========================

# Create spacecraft
sat = Spacecraft(
    state=CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]), 
    time=Time("2015-09-21T12:23:12", TAI(), ISOT())
    )

# Create force models, integrator, and dynamics system
pm_grav = PointMassGravity(earth,(moon,sun))
forces = ForceModel(pm_grav)
integ = IntegratorConfig(DP8(); abstol = 1e-11, reltol = 1e-11, dt = 4000)

# Define which spacecraft to propagate and which force model to use
dynsys = DynSys(
           forces = forces,
           spacecraft =  [sat]
            )

# Loop with point mass gravity for 10 days in LEO
for i in 1:10
    local val, elapsed_time, bytes, gctime2
    val, elapsed_time, bytes, gctime2 = @timed propagate!(dynsys, integ, StopAtSeconds(864000.0); prop_stm = true)
    println("Elapsed: $elapsed_time s, GC time: $gctime2 s, Allocated: $bytes bytes")
end

