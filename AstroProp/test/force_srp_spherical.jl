# Solar radiation pressure — open, Epicycle interface vs GMAT.
# Two-body (point-mass Earth) + cannonball SRP, near-GEO, 1 day. Dependency-light (@assert + println).
#
# NOTE — assumes GMAT used POINT-MASS Earth gravity in this force model. At GEO a gravity field
# (e.g. the default 4×4) would differ by ~km over a day; if that's what the run used, switch the
# gravity here to match (and it becomes an EGM test).
#
# GMAT config: SRP On, Flux = 1367, SRPModel = Spherical (cannonball), Nominal_Sun = 149597870.691;
#              Cr = 1.8, SRPArea = 10, DryMass = 1000. Epoch/duration as the other scripts.

using Test
using AstroProp
using AstroModels, AstroStates, AstroEpochs
using AstroUniverse: earth
using OrdinaryDiffEq: Vern9
using LinearAlgebra: norm
earth.mu = 398600.4415

@testset "SRP vs GMAT" begin

sc = Spacecraft(;
    state = CartesianState([41795.780, 0.0, 0.0, 0.0, 3.0838, 0.26980]),
    time  = Time("2020-10-20T12:00:00", UTC(), ISOT()),
    mass  = 1000.0,
    name  = "GEO",
    srp   = SphericalSRP(c_r = 1.8, srp_area = 10.0),
)

forces = ForceModel(
    PointMassGravity(earth, ()),
    SolarRadiationPressure(earth; shadow = DualCone(), solar_flux = 1367.0, nominal_sun = 149597870.691),
)

integ = IntegratorConfig(Vern9(); reltol = 1e-12, abstol = 1e-12, dt = 2700.0)
prop  = OrbitPropagator(forces, integ)
sol   = propagate!(prop, sc, StopAt(sc, PropDurationSeconds(), 5*86400.0))
yf    = sol.u[end]

# GMAT truth: two-body + cannonball SRP, 5 days.
gmat = [ 40213.595533343,  11374.748122627,  995.17095544591,
           -0.8414481463335, 2.9671200023382, 0.2595917089176]

Δr = norm(gmat[1:3] .- yf[1:3]) * 1e3      # m
Δv = norm(gmat[4:6] .- yf[4:6]) * 1e6      # mm/s

@test Δr < 50.0
@test Δv < 50.0
end
