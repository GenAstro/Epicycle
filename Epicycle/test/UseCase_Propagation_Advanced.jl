include("UseCase_Init.jl")

# Create the four mms spacecraft
mms1 = Spacecraft(
    state = KeplerianState(42256.691, 0.80229, deg2rad(10.0481), deg2rad(90.2462), deg2rad(-0.06204), deg2rad(180.1627)), 
    time = Time("2015-09-21T12:23:12", TAI(), ISOT())
    )

mms2 = Spacecraft(
    state = KeplerianState(42256.691, 0.80221, deg2rad(10.0493), deg2rad(90.2623), deg2rad(-0.06873),deg2rad(180.1585)),
    time = Time("2015-09-21T12:23:12", TAI(), ISOT())
    )

mms3 = Spacecraft(
    state = KeplerianState(42256.691, 0.80214, deg2rad(10.0486), deg2rad(90.2540), deg2rad(-0.07096), deg2rad(180.1634)), 
    time = Time("2015-09-21T12:23:12", TAI(), ISOT())
    )

mms4 = Spacecraft(
    state = KeplerianState(42256.691, 0.80220, deg2rad(10.0426), deg2rad(90.2591), deg2rad(-0.07084), deg2rad(180.1608)), 
    time = Time("2015-09-21T12:23:12", TAI(), ISOT())
    )

# Create force models, integrator, and dynamics system
pm_grav = PointMassGravity(earth,(moon,sun))
forces = ForceModel(pm_grav)
integ = IntegratorConfig(DP8(); abstol = 1e-11, reltol = 1e-11, dt = 4000)

# Create a model that includes all four spacecraft
dynsys = DynSys(
                forces = forces, 
                spacecraft = [mms1,mms2,mms3,mms4]
                )

# Propagate to apoapsis of mms1
propagate(dynsys, integ, StopAtApoapsis(mms1))
println(KeplerianState(mms1.state,earth.mu))

# Propagate backwards for an hour
println((mms1.time.tai).isot)
t1 = mms1.time.jd;
propagate(dynsys, integ, StopAtSeconds(-3600.0); direction = :backward);
t2 = mms2.time.jd;
(t2 - t1)*86400
println((mms1.time.tai).isot)

