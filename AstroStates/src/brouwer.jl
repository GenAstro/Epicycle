# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation
#
# Brouwer-Lyddane mean element conversions (BROU), reimplemented in idiomatic Julia
# from GMAT's StateConversionUtil.cpp (NASA, public domain). Radians throughout;
# distances in km. Earth-only in this slice (constants hard-coded; guarded).
#
# GMAT structure: Cartesian -> osculating Keplerian -> (iterate) -> mean. Per project
# direction these public routines take/return OSCULATING KEPLERIAN directly:
#   kep_to_brouwer_mean_{long,short}(osc[a,e,i,Ω,ω,ν], μ) -> mean[a,e,i,Ω,ω,M]
#   brouwer_mean_{long,short}_to_kep(mean[a,e,i,Ω,ω,M], μ) -> osc[a,e,i,Ω,ω,ν]

# ---- GMAT-consistent Earth constants (transcribed from StateConversionUtil.cpp) ----
const _BROUWER_RE       = 6378.1363
const _BROUWER_J2       = 1.082626925638815e-3
const _BROUWER_J3       = -0.2532307818191774e-5
const _BROUWER_J4       = -0.1620429990000000e-5
const _BROUWER_J5       = -0.2270711043920343e-6
const _BROUWER_MU_EARTH = 398600.4415

"""
    _brouwer_guard_mu(μ, name="Brouwer")

Assert `μ` is Earth's gravitational parameter (within 1 km³/s²); the Brouwer element theory is
Earth-only in this release. Errors with a diagnostic naming `name` otherwise. Internal.
"""
@inline function _brouwer_guard_mu(μ::Real, name::AbstractString="Brouwer")
    if abs(μ - _BROUWER_MU_EARTH) > 1.0
        error("$name is applicable only to the Earth (require |μ - $(_BROUWER_MU_EARTH)| ≤ 1; got μ=$μ).")
    end
    return nothing
end

# ---- local anomaly helpers (AstroStates must not depend on AstroRoutines; O-1) ----
"""
    _true_to_mean_anom(ν, e)

Convert true anomaly `ν` to mean anomaly `M` for an elliptic orbit via the eccentric anomaly:
`E = atan2(√(1-e²)·sin ν, e + cos ν)`, `M = E - e·sin E`. Angles in radians. Internal helper so
`AstroStates` need not depend on the anomaly-conversion package across the public/dev seam.
"""
@inline function _true_to_mean_anom(ν::Real, e::Real)
    s = sqrt(max(1.0 - e^2, 0.0))
    E = atan(s*sin(ν), e + cos(ν))          # eccentric anomaly, (-π, π]
    return E - e*sin(E)
end

"""
    _mean_to_true_anom(M, e; tol=1e-13, maxiter=100)

Convert mean anomaly `M` to true anomaly `ν` in `[0, 2π)` for an elliptic orbit by Newton iteration
on Kepler's equation `E - e·sin E = M` (Danby starter), then `ν = atan2(√(1-e²)·sin E, cos E - e)`.
Angles in radians. Internal helper (see `_true_to_mean_anom`).
"""
function _mean_to_true_anom(M::Real, e::Real; tol::Real=1e-13, maxiter::Int=100)
    Mr = mod(M, 2π)
    E  = Mr + (sin(Mr) < 0 ? -1.0 : 1.0)*0.85*e     # Danby starter (robust to high e)
    for _ in 1:maxiter
        f = E - e*sin(E) - Mr
        abs(f) < tol && break
        E -= f/(1.0 - e*cos(E))
    end
    ν = atan(sqrt(max(1.0 - e^2, 0.0))*sin(E), cos(E) - e)
    ν < 0 && (ν += 2π)
    return ν
end

# ---- equinoctial-like helpers used by the osculating->mean fixed point (GMAT aeq) ----
"""
    _brouwer_kep_to_eq6(k)

Map Keplerian-with-mean-anomaly `[a, e, i, Ω, ω, M]` to the equinoctial-like 6-tuple GMAT differences
in during the osculating→mean iteration: `(a, e·sin(ω+Ω), e·cos(ω+Ω), sin(i/2)·sin Ω, sin(i/2)·cos Ω, Ω+ω+M)`.
Internal.
"""
@inline function _brouwer_kep_to_eq6(k)
    a, e, i, Ω, ω, M = k
    return SVector{6}(a,
                      e*sin(ω+Ω),
                      e*cos(ω+Ω),
                      sin(i/2)*sin(Ω),
                      sin(i/2)*cos(Ω),
                      Ω+ω+M)
end

"""
    _brouwer_eq6_to_kep(q)

Inverse of `_brouwer_kep_to_eq6`: recover `[a, e, i, Ω, ω, M]` (Ω, ω wrapped to `[0, 2π)`) from the
equinoctial-like 6-vector. Internal.
"""
@inline function _brouwer_eq6_to_kep(q)
    a = q[1]
    e = sqrt(q[2]^2 + q[3]^2)
    s = q[4]^2 + q[5]^2
    i = acos(1.0 - 2.0*min(s, 1.0))
    Ω = mod(atan(q[4], q[5]), 2π)
    ω = mod(atan(q[2], q[3]) - Ω, 2π)
    M = q[6] - atan(q[2], q[3])
    return (a, e, i, Ω, ω, M)
end

# ============================================================================
# CORE (direct): Brouwer mean -> osculating Keplerian (osc carries MEAN anomaly)
# Faithful port of BrouwerMeanShortToOsculatingElements (radians, km).
# ============================================================================
"""
    _brouwer_short_mean_to_osc(blms, μ)

Direct Brouwer short-period (J2) transform from mean elements `blms = [a, e, i, Ω, ω, M]` to
osculating Keplerian elements `[a, e, i, Ω, ω, M]` (osculating **mean anomaly** in slot 6).
Radians, km, Earth-only. Faithful port of GMAT `BrouwerMeanShortToOsculatingElements`. Internal.
"""
function _brouwer_short_mean_to_osc(blms::AbstractVector, μ::Real)
    _brouwer_guard_mu(μ, "BrouwerMeanShort")
    re = _BROUWER_RE
    j2 = _BROUWER_J2
    smap     = blms[1]/re
    eccp     = float(blms[2])
    incp     = float(blms[3])
    raanp    = float(blms[4])
    aopp     = float(blms[5])
    meanAnom = float(blms[6])

    (incp < 0.0 || incp > π) && error("BrouwerMeanShort: mean INC must be in [0, 180°] (got $(rad2deg(incp))°).")
    radper = blms[1]*(1.0 - blms[2])
    radper < 3000.0 && error("BrouwerMeanShort: mean RadPer must be > 3000 km (got $radper).")
    radper < 6378.0 && @warn "BrouwerMeanShort: mean RadPer < 6378 km — possible inside-Earth singularity inaccuracy."
    if eccp < 0.0                                            # COV_EXCL_START
        eccp = -eccp; meanAnom -= π; aopp += π
    end                                                      # COV_EXCL_STOP
    eccp > 0.99 && error("BrouwerMeanShort: mean ECC must be < 0.99 (got $eccp).")

    pseudostate = 0
    if incp > deg2rad(175.0)
        incp = π - incp; raanp = -raanp; pseudostate = 1
    end
    raanp = mod(raanp, 2π); aopp = mod(aopp, 2π); meanAnom = mod(meanAnom, 2π)

    eta   = sqrt(1.0 - eccp^2)
    theta = cos(incp)
    p     = smap*eta^2
    k2    = 0.5*j2
    gm2   = k2/smap^2
    gm2p  = gm2/eta^4

    tap = _mean_to_true_anom(meanAnom, eccp; tol=1e-13)
    rp  = p/(1.0 + eccp*cos(tap))
    adr = smap/rp

    sma1 = smap + smap*gm2*((adr^3 - 1.0/eta^3)*(-1.0 + 3.0*theta^2) +
           3.0*(1.0 - theta^2)*adr^3*cos(2.0*aopp + 2.0*tap))

    decc = eta^2/2.0*(3.0*(1.0/eta^6)*gm2*(1.0 - theta^2)*cos(2.0*aopp + 2.0*tap)*
             (3.0*eccp*cos(tap)^2 + 3.0*cos(tap) + eccp^2*cos(tap)^3 + eccp)
           - gm2p*(1.0 - theta^2)*(3.0*cos(2.0*aopp + tap) + cos(3.0*tap + 2.0*aopp))
           + (3.0*theta^2 - 1.0)*gm2/eta^6*(eccp*eta + eccp/(1.0 + eta) +
             3.0*eccp*cos(tap)^2 + 3.0*cos(tap) + eccp^2*cos(tap)^3))

    dinc = gm2p/2.0*theta*sin(incp)*(3.0*cos(2.0*aopp + 2.0*tap) +
             3.0*eccp*cos(2.0*aopp + tap) + eccp*cos(2.0*aopp + 3.0*tap))

    draan = -gm2p/2.0*theta*(6.0*(tap - meanAnom + eccp*sin(tap)) -
              3.0*sin(2.0*aopp + 2.0*tap) - 3.0*eccp*sin(2.0*aopp + tap) -
              eccp*sin(2.0*aopp + 3.0*tap))

    # NOTE: GMAT computes a short-period aop1 and ma1 here, but both are overwritten
    # below (dead code in the C++). Omitted; results identical.

    lgh = raanp + aopp + meanAnom + gm2p/4.0*(6.0*(-1.0 - 2.0*theta + 5.0*theta^2)*
             (tap - meanAnom + eccp*sin(tap)) +
             (3.0 + 2.0*theta - 5.0*theta^2)*(3.0*sin(2.0*aopp + 2.0*tap) +
               3.0*eccp*sin(2.0*aopp + tap) + eccp*sin(2.0*aopp + 3.0*tap))) +
          gm2p/4.0*eta^2/(eta + 1.0)*eccp*(3.0*(1.0 - theta^2)*
             (sin(3.0*tap + 2.0*aopp)*(1.0/3.0 + adr^2*eta^2 + adr) +
              sin(2.0*aopp + tap)*(1.0 - adr^2*eta^2 - adr)) +
             2.0*sin(tap)*(3.0*theta^2 - 1.0)*(1.0 + adr^2*eta^2 + adr))

    eccpdl = -eta^3/4.0*gm2p*(2.0*(-1.0 + 3.0*theta^2)*(adr^2*eta^2 + adr + 1.0)*sin(tap) +
               3.0*(1.0 - theta^2)*((-adr^2*eta^2 - adr + 1.0)*sin(2.0*aopp + tap) +
                 (adr^2*eta^2 + adr + 1.0/3.0)*sin(2.0*aopp + 3.0*tap)))

    ecosl = (eccp + decc)*cos(meanAnom) - eccpdl*sin(meanAnom)
    esinl = (eccp + decc)*sin(meanAnom) + eccpdl*cos(meanAnom)
    ecc1  = sqrt(ecosl^2 + esinl^2)
    ma1   = 0.0
    if ecc1 >= 1e-11
        ma1 = atan(esinl, ecosl)
        ma1 < 0.0 && (ma1 += 2π)
    end

    sinhalfisinh = (sin(0.5*incp) + cos(0.5*incp)*0.5*dinc)*sin(raanp) +
                   0.5*sin(incp)/cos(incp/2.0)*draan*cos(raanp)
    sinhalficosh = (sin(0.5*incp) + cos(0.5*incp)*0.5*dinc)*cos(raanp) -
                   0.5*sin(incp)/cos(incp/2.0)*draan*sin(raanp)
    mag  = sqrt(sinhalfisinh^2 + sinhalficosh^2)
    inc1 = mag > 1.0 ? 2.0*asin(1.0) : (mag < -1.0 ? 2.0*asin(-1.0) : 2.0*asin(mag))

    local raan1, aop1
    if inc1 == 0.0     # GMAT also tests inc1==180.0, dead in radians
        raan1 = 0.0
        aop1  = lgh - ma1 - raan1
    else
        raan1 = atan(sinhalfisinh, sinhalficosh)
        raan1 < 0.0 && (raan1 += 2π)
        aop1  = lgh - ma1 - raan1
    end
    aop1 = mod(aop1, 2π)

    inc1o, raan1o = inc1, raan1
    if pseudostate != 0
        inc1o  = π - inc1
        raan1o = 2π - raan1
    end
    return [sma1*re, ecc1, inc1o, raan1o, aop1, ma1]
end

# ============================================================================
# CORE (direct): Brouwer-Lyddane mean(long) -> osculating Keplerian (MEAN anomaly)
# Faithful port of BrouwerMeanLongToOsculatingElements (radians, km).
# ============================================================================
"""
    _brouwer_long_mean_to_osc(blml, μ)

Direct Brouwer-Lyddane long-period (J2–J5) transform from mean elements `blml = [a, e, i, Ω, ω, M]`
to osculating Keplerian elements `[a, e, i, Ω, ω, M]` (osculating **mean anomaly** in slot 6).
Zeroes the long-period terms near the critical inclination (GMAT `bisubc ≥ 0.001`). Radians, km,
Earth-only. Faithful port of GMAT `BrouwerMeanLongToOsculatingElements`. Internal.
"""
function _brouwer_long_mean_to_osc(blml::AbstractVector, μ::Real)
    _brouwer_guard_mu(μ, "BrouwerMeanLong")
    re = _BROUWER_RE
    j2 = _BROUWER_J2; j3 = _BROUWER_J3; j4 = _BROUWER_J4; j5 = _BROUWER_J5
    smadp    = blml[1]/re
    eccdp    = float(blml[2])
    incdp    = float(blml[3])
    raandp   = float(blml[4])
    aopdp    = float(blml[5])
    meanAnom = float(blml[6])

    pseudostate = 0
    if incdp > deg2rad(175.0)
        incdp = π - incdp; raandp = -raandp; pseudostate = 1
    end
    eccdp > 0.99 && error("BrouwerMeanLong: mean ECC must be < 0.99 (got $eccdp).")
    radper = blml[1]*(1.0 - blml[2])
    radper < 3000.0 && error("BrouwerMeanLong: mean RadPer must be > 3000 km (got $radper).")
    radper < 6378.0 && @warn "BrouwerMeanLong: mean RadPer < 6378 km — possible inside-Earth singularity inaccuracy."
    blml[3] > π && error("BrouwerMeanLong: INC must be < 180° (got $(rad2deg(blml[3]))°).")
    raandp = mod(raandp, 2π); aopdp = mod(aopdp, 2π); meanAnom = mod(meanAnom, 2π)

    bk2 = 0.5*(j2)
    bk3 = -j3
    bk4 = -(3.0/8.0)*j4
    bk5 = -j5
    eccdp2 = eccdp^2
    cn2 = 1.0 - eccdp2
    cn  = sqrt(cn2)
    gm2  = bk2/smadp^2
    gmp2 = gm2/(cn2*cn2)
    gm4  = bk4/smadp^4
    gmp4 = gm4/cn^8
    theta = cos(incdp); theta2 = theta^2; theta4 = theta2^2
    gm3  = bk3/smadp^3
    gmp3 = gm3/(cn2*cn2*cn2)
    gm5  = bk5/smadp^5
    gmp5 = gm5/cn^10
    g3dg2 = gmp3/gmp2
    g4dg2 = gmp4/gmp2
    g5dg2 = gmp5/gmp2
    sinMADP = sin(meanAnom); cosMADP = cos(meanAnom)
    sinraandp = sin(raandp); cosraandp = cos(raandp)

    tadp = _mean_to_true_anom(meanAnom, eccdp; tol=1e-13)
    rp  = smadp*(1.0 - eccdp2)/(1.0 + eccdp*cos(tadp))
    adr = smadp/rp
    sinta = sin(tadp); costa = cos(tadp)
    cs2gta = cos(2.0*aopdp + 2.0*tadp)
    adr2 = adr^2; adr3 = adr2*adr; costa2 = costa^2

    a1 = (0.125*gmp2*cn2)*(1.0 - 11.0*theta2 - 40.0*theta4/(1.0 - 5.0*theta2))
    a2 = ((5.0/12.0)*g4dg2*cn2)*(1.0 - 8.0*theta4/(1.0 - 5.0*theta2) - 3.0*theta2)
    a3 = g5dg2*(3.0*eccdp2 + 4.0)
    a4 = g5dg2*(1.0 - 24.0*theta4/(1.0 - 5.0*theta2) - 9.0*theta2)
    a5 = (g5dg2*(3.0*eccdp2 + 4.0))*(1.0 - 24.0*theta4/(1.0 - 5.0*theta2) - 9.0*theta2)
    a6 = g3dg2*0.25
    sinI = sin(incdp)
    a10 = cn2*sinI
    a7  = a6*a10
    a8p = g5dg2*eccdp*(1.0 - 16.0*theta4/(1.0 - 5.0*theta2) - 5.0*theta2)
    a8  = a8p*eccdp
    b13 = eccdp*(a1 - a2)
    b14 = a7 + (5.0/64.0)*a5*a10
    b15 = a8*a10*(35.0/384.0)
    a11 = 2.0 + eccdp2
    a12 = 3.0*eccdp2 + 2.0
    a13 = theta2*a12
    a14 = (5.0*eccdp2 + 2.0)*(theta4/(1.0 - 5.0*theta2))
    a17 = theta4/((1.0 - 5.0*theta2)^2)
    a15 = (eccdp2*theta4*theta2)/((1.0 - 5.0*theta2)^2)
    a16 = theta2/(1.0 - 5.0*theta2)
    a18 = eccdp*sinI
    a19 = a18/(1.0 + cn)
    a21 = eccdp*theta
    a22 = eccdp2*theta
    sinI2 = sin(incdp/2.0); cosI2 = cos(incdp/2.0); tanI2 = tan(incdp/2.0)
    a26 = 16.0*a16 + 40.0*a17 + 3.0
    a27 = a22*0.125*(11.0 + 200.0*a17 + 80.0*a16)

    b1 = cn*(a1 - a2) - ((a11 - 400.0*a15 - 40.0*a14 - 11.0*a13)*(1.0/16.0) +
           (11.0 + 200.0*a17 + 80.0*a16)*a22*0.125)*gmp2 +
         ((-80.0*a15 - 8.0*a14 - 3.0*a13 + a11)*(5.0/24.0) + (5.0/12.0)*a26*a22)*g4dg2
    b2 = a6*a19*(2.0 + cn - eccdp2) + (5.0/64.0)*a5*a19*cn2 - (15.0/32.0)*a4*a18*cn*cn2 +
         ((5.0/64.0)*a5 + a6)*a21*tanI2 + (9.0*eccdp2 + 26.0)*(5.0/64.0)*a4*a18 +
         (15.0/32.0)*a3*a21*a26*sinI*(1.0 - theta)
    b3 = ((80.0*a17 + 5.0 + 32.0*a16)*a22*sinI*(theta - 1.0)*(35.0/576.0)*g5dg2*eccdp) -
         ((a22*tanI2 + (2.0*eccdp2 + 3.0*(1.0 - cn2*cn))*sinI)*(35.0/1152.0)*a8p)
    b4 = cn*eccdp*(a1 - a2)
    b5 = ((9.0*eccdp2 + 4.0)*a10*a4*(5.0/64.0) + a7)*cn
    b6 = (35.0/384.0)*a8*cn2*cn*sinI
    b7 = ((cn2*a18)/(1.0 - 5.0*theta2))*(0.125*gmp2*(1.0 - 15.0*theta2) +
           (1.0 - 7.0*theta2)*g4dg2*(-(5.0/12.0)))
    b8 = (5.0/64.0)*(a3*cn2*(1.0 - 9.0*theta2 - 24.0*theta4/(1.0 - 5.0*theta2))) + a6*cn2
    b9 = a8*(35.0/384.0)*cn2
    b10 = sinI*(a22*a26*g4dg2*(5.0/12.0) - a27*gmp2)
    b11 = a21*(a5*(5.0/64.0) + a6 + a3*a26*(15.0/32.0)*sinI^2)
    b12 = -((80.0*a17 + 32.0*a16 + 5.0)*(a22*eccdp*sinI^2*(35.0/576.0)*g5dg2) +
            (a8*a21*(35.0/1152.0)))

    sma = smadp*(1.0 + gm2*((3.0*theta2 - 1.0)*(eccdp2/(cn2*cn2*cn2))*(cn + 1.0/(1.0 + cn)) +
            ((3.0*theta2 - 1.0)/(cn2*cn2*cn2))*(eccdp*costa)*(3.0 + 3.0*eccdp*costa + eccdp2*costa2) +
            3.0*(1.0 - theta2)*adr3*cs2gta))
    sn2gta = sin(2.0*aopdp + 2.0*tadp)
    snf2gd = sin(2.0*aopdp + tadp)
    csf2gd = cos(2.0*aopdp + tadp)
    sn2gd  = sin(2.0*aopdp)
    cs2gd  = cos(2.0*aopdp)
    sin3gd = sin(3.0*aopdp)
    cs3gd  = cos(3.0*aopdp)
    sn3fgd = sin(3.0*tadp + 2.0*aopdp)
    cs3fgd = cos(3.0*tadp + 2.0*aopdp)
    sinGD  = sin(aopdp)
    cosGD  = cos(aopdp)

    bisubc = ((1.0 - 5.0*theta2)^(-2))*((25.0*theta4*theta)*(gmp2*eccdp2))
    local blghp, eccdpdl, dltI, sinDH, dlt1e
    if bisubc >= 0.001
        @warn "BrouwerMeanLong: mean INC near critical (63°/117°) — possible inaccuracy."
        dlt1e = 0.0; blghp = 0.0; eccdpdl = 0.0; dltI = 0.0; sinDH = 0.0
    else
        blghp = raandp + aopdp + meanAnom + b3*cs3gd + b1*sn2gd + b2*cosGD
        blghp = mod(blghp, 2π)
        dlt1e = b14*sinGD + b13*cs2gd - b15*sin3gd
        eccdpdl = b4*sn2gd - b5*cosGD + b6*cs3gd -
                  0.25*cn2*cn*gmp2*(2.0*(3.0*theta2 - 1.0)*(adr2*cn2 + adr + 1.0)*sinta +
                    3.0*(1.0 - theta2)*((-adr2*cn2 - adr + 1.0)*snf2gd +
                      (adr2*cn2 + adr + (1.0/3.0))*sn3fgd))
        dltI = 0.5*theta*gmp2*sinI*(eccdp*cs3fgd + 3.0*(eccdp*csf2gd + cs2gta)) -
               (a21/cn2)*(b8*sinGD + b7*cs2gd - b9*sin3gd)
        sinDH = (1.0/cosI2)*(0.5*(b12*cs3gd + b11*cosGD + b10*sn2gd -
                  (0.5*gmp2*theta*sinI*(6.0*(eccdp*sinta - meanAnom + tadp) -
                    (3.0*(sn2gta + eccdp*snf2gd) + eccdp*sn3fgd)))))
    end

    blgh = blghp + ((1.0/(cn + 1.0))*0.25*eccdp*gmp2*cn2*(3.0*(1.0 - theta2)*
             (sn3fgd*((1.0/3.0) + adr2*cn2 + adr) + snf2gd*(1.0 - (adr2*cn2 + adr))) +
             2.0*sinta*(3.0*theta2 - 1.0)*(adr2*cn2 + adr + 1.0))) +
           gmp2*1.5*((-2.0*theta - 1.0 + 5.0*theta2)*(eccdp*sinta + tadp - meanAnom)) +
           (3.0 + 2.0*theta - 5.0*theta2)*(gmp2*0.25*(eccdp*sn3fgd +
             3.0*(sn2gta + eccdp*snf2gd)))
    blgh = mod(blgh, 2π)

    dlte = dlt1e + (0.5*cn2*((3.0*(1.0/(cn2*cn2*cn2))*gm2*(1.0 - theta2)*cs2gta*
             (3.0*eccdp*costa2 + 3.0*costa + eccdp2*costa*costa2 + eccdp)) -
             (gmp2*(1.0 - theta2)*(3.0*csf2gd + cs3fgd)) +
             (3.0*theta2 - 1.0)*gm2*(1.0/(cn2*cn2*cn2))*(eccdp*cn + (eccdp/(1.0 + cn)) +
               3.0*eccdp*costa2 + 3.0*costa + eccdp2*costa*costa2)))
    ecc  = sqrt(eccdpdl^2 + (eccdp + dlte)^2)
    squar = (dltI*cosI2*0.5 + sinI2)^2
    sqrI  = sqrt(sinDH^2 + squar)
    inc  = 2.0*asin(sqrI)
    inc  = mod(inc, 2π)

    local ma, raan, aop
    if ecc <= 1e-11                                          # COV_EXCL_START
        aop = 0.0
        if inc <= 1e-7
            raan = 0.0; ma = blgh
        else
            arg1 = sinDH*cosraandp + sinraandp*(0.5*dltI*cosI2 + sinI2)
            arg2 = cosraandp*(0.5*dltI*cosI2 + sinI2) - sinDH*sinraandp
            raan = atan(arg1, arg2); ma = blgh - aop - raan
        end                                                  # COV_EXCL_STOP
    else
        arg1 = eccdpdl*cosMADP + (eccdp + dlte)*sinMADP
        arg2 = (eccdp + dlte)*cosMADP - eccdpdl*sinMADP
        ma = atan(arg1, arg2); ma = mod(ma, 2π)
        if inc <= 1e-7
            raan = 0.0; aop = blgh - raan - ma
        else
            arg1b = sinDH*cosraandp + sinraandp*(0.5*dltI*cosI2 + sinI2)
            arg2b = cosraandp*(0.5*dltI*cosI2 + sinI2) - sinDH*sinraandp
            raan = atan(arg1b, arg2b); aop = blgh - ma - raan
        end
    end
    ma < 0.0 && (ma += 2π)
    raan = mod(raan, 2π); aop = mod(aop, 2π)

    inc_o, raan_o = inc, raan
    if pseudostate != 0
        inc_o  = π - inc
        raan_o = 2π - raan
    end
    return [sma*re, ecc, inc_o, raan_o, aop, ma]
end

# ============================================================================
# Iterative osculating Keplerian -> Brouwer mean (GMAT CartesianTo... engine,
# driven from osculating Keplerian per project direction).
# core = mean->osc function; use_proposal selects GMAT's final vector (short uses
# the proposal aeqmean2; long uses the accepted aeqmean).
# ============================================================================
"""
    _brouwer_osc_to_mean(osc, μ, core; use_proposal)

Iterative osculating→mean solve driven from osculating Keplerian `osc = [a, e, i, Ω, ω, ν]`, using
`core` (a mean→osculating transform) inside a Cartesian-error fixed point. `use_proposal` selects
GMAT's final vector (short uses the proposal, long the accepted iterate). Returns `(mean, pseudostate)`
with mean `[a, e, i, Ω, ω, M]`. On near-singular divergence it warns and returns the best estimate.
Faithful port of GMAT's `CartesianToBrouwerMean*` iteration. Internal.
"""
function _brouwer_osc_to_mean(osc::AbstractVector, μ::Real, core; use_proposal::Bool)
    tol = 1e-8
    maxiter = 75
    a, e, i, Ω, ω, ν = float.(osc)
    cart = kep_to_cart(collect(float.(osc)), μ)
    M = _true_to_mean_anom(ν, e)
    kep = SVector{6}(a, e, i, Ω, ω, M)

    pseudostate = 0
    if i > deg2rad(175.0)
        i2 = π - i; Ω2 = -Ω
        kep  = SVector{6}(a, e, i2, Ω2, ω, M)
        cart = kep_to_cart([a, e, i2, Ω2, ω, _mean_to_true_anom(M, e)], μ)
        pseudostate = 1
    end

    kep2     = SVector{6}(core(kep, μ)...)
    aeq      = _brouwer_kep_to_eq6(kep)
    aeq2     = _brouwer_kep_to_eq6(kep2)
    aeqmean  = _brouwer_kep_to_eq6(kep)
    aeqmean2 = aeqmean + (aeq - aeq2)

    emag = 0.9; emag_old = 1.0; ii = 0
    diverged = false
    while emag > tol
        blmean2 = _brouwer_eq6_to_kep(aeqmean2)
        # A near-singular mean estimate (e.g. near critical inclination at high
        # eccentricity) can overshoot the core's domain. Match GMAT: interrupt,
        # warn, and return the best estimate so far rather than throwing.
        kep2v = try
            core(collect(blmean2), μ)
        catch
            diverged = true
            nothing
        end
        diverged && break
        kep2    = SVector{6}(kep2v...)
        kep2TA  = [kep2[1], kep2[2], kep2[3], kep2[4], kep2[5], _mean_to_true_anom(kep2[6], kep2[2])]
        cart2   = kep_to_cart(kep2TA, μ)
        emag    = norm(cart .- cart2)/norm(cart)
        if emag_old > emag
            emag_old = emag
            aeq2     = _brouwer_kep_to_eq6(kep2)
            aeqmean  = aeqmean2
            aeqmean2 = aeqmean + (aeq - aeq2)
        else
            break
        end
        ii > maxiter && break
        ii += 1
    end
    diverged && @warn "Brouwer osculating→mean iteration did not converge (near-singular regime); returning best estimate."

    final = diverged ? aeqmean : (use_proposal ? aeqmean2 : aeqmean)
    b = collect(float.(_brouwer_eq6_to_kep(final)))
    return b, pseudostate
end

# ============================================================================
# Public math API (radian vectors)
# ============================================================================
"""
    kep_to_brouwer_mean_short(osc::AbstractVector, μ::Real)

Convert osculating Keplerian elements to Brouwer mean (short-period-only, J2) elements.

# Arguments
- `osc::AbstractVector`: osculating Keplerian elements `[a, e, i, Ω, ω, ν]` (see `kep_to_brouwer_mean_long`)
- `μ`: Gravitational parameter (Earth only)

# Returns
A 6-element vector `[a, e, i, Ω, ω, M]` of mean elements, where `M` is the **mean anomaly** (rad).

# Example
mean = kep_to_brouwer_mean_short([7000.0, 0.01, pi/6, 0.0, 0.0, pi/3], 398600.4415)

# Notes
- As `kep_to_brouwer_mean_long`, but removes only the first-order J2 short-period variation
  (long-period and secular terms are retained); has no critical-inclination degradation.
- Inverse of `brouwer_mean_short_to_kep`. See `kep_to_brouwer_mean_long` for the shared units,
  domain, and differentiability contract.
"""
function kep_to_brouwer_mean_short(osc::AbstractVector, μ::Real)
    _brouwer_guard_mu(μ, "BrouwerMeanShort")
    a, e, i = osc[1], osc[2], osc[3]
    i > π && error("BrouwerMeanShort: osculating INC must be ≤ 180° (got $(rad2deg(i))°).")
    (e > 0.99 || e < 0.0) && error("BrouwerMeanShort: applicable only for 0 ≤ ECC < 0.99 (got e=$e).")
    radper = a*(1.0 - e)
    radper < 3000.0 && error("BrouwerMeanShort: osculating RadPer must be > 3000 km (got $radper).")
    radper < 6378.0 && @warn "BrouwerMeanShort: RadPer < 6378 km — possible inside-Earth singularity inaccuracy."
    b, pseudostate = _brouwer_osc_to_mean(osc, μ, _brouwer_short_mean_to_osc; use_proposal=true)
    a1, e1, i1, Ω1, ω1, M1 = b
    if e1 < 0.0                                              # COV_EXCL_START
        e1 = -e1; ω1 += π; M1 -= π
    end                                                      # COV_EXCL_STOP
    if pseudostate != 0
        i1 = π - i1; Ω1 = -Ω1
    end
    return [a1, e1, i1, mod(Ω1, 2π), mod(ω1, 2π), mod(M1, 2π)]
end

"""
    kep_to_brouwer_mean_long(osc::AbstractVector, μ::Real)

Convert osculating Keplerian elements to Brouwer-Lyddane mean (long-period) elements.

# Arguments
- `osc::AbstractVector`: osculating Keplerian elements `[a, e, i, Ω, ω, ν]`
- `μ`: Gravitational parameter (Earth only)
- `a`: semi-major axis (km)
- `e`: eccentricity
- `i`: inclination (rad)
- `Ω`: right ascension of ascending node (rad)
- `ω`: argument of periapsis (rad)
- `ν`: true anomaly (rad)

# Returns
A 6-element vector `[a, e, i, Ω, ω, M]` of mean elements, where `M` is the **mean anomaly** (rad).

# Example
mean = kep_to_brouwer_mean_long([7000.0, 0.01, pi/6, 0.0, 0.0, pi/3], 398600.4415)

# Notes
- Angles in radians; distances in km consistent with the Earth `μ`.
- Earth-only: errors if `μ` is not Earth's gravitational parameter.
- Domain: `0 ≤ e < 0.99` and mean periapsis radius `a(1-e) > 3000 km`. Warns when periapsis is
  below 6378 km (inside-Earth) and near the critical inclination (≈ 63.4°/116.6°), where the
  long-period theory is degraded.
- Removes secular, long-period, and short-period zonal (J2–J5) variations. Inverse of
  `brouwer_mean_long_to_kep`.
- Differentiable (ForwardDiff) with respect to `osc`.
- On a near-singular, non-convergent input the osculating→mean iteration warns and returns its
  best estimate rather than throwing.
- Brouwer (1959) / Lyddane (1963) theory; transcribed from GMAT `StateConversionUtil.cpp`.
"""
function kep_to_brouwer_mean_long(osc::AbstractVector, μ::Real)
    _brouwer_guard_mu(μ, "BrouwerMeanLong")
    a, e, i = osc[1], osc[2], osc[3]
    (e > 0.99 || e < 0.0) && error("BrouwerMeanLong: applicable only for 0 ≤ ECC < 0.99 (got e=$e).")
    radper = a*(1.0 - e)
    radper < 3000.0 && error("BrouwerMeanLong: osculating RadPer must be > 3000 km (got $radper).")
    radper < 6378.0 && @warn "BrouwerMeanLong: RadPer < 6378 km — possible inside-Earth singularity inaccuracy."
    i > π && error("BrouwerMeanLong: osculating INC must be ≤ 180° (got $(rad2deg(i))°).")
    idg = rad2deg(i)
    if (58.80 < idg < 65.78) || (114.22 < idg < 121.2)
        @warn "BrouwerMeanLong: INC near critical (63°/117°) — possible inaccuracy."
    end
    b, pseudostate = _brouwer_osc_to_mean(osc, μ, _brouwer_long_mean_to_osc; use_proposal=false)
    a1, e1, i1, Ω1, ω1, M1 = b
    if pseudostate != 0
        i1 = π - i1; Ω1 = -Ω1
    end
    return [a1, e1, i1, mod(Ω1, 2π), mod(ω1, 2π), mod(M1, 2π)]
end

"""
    brouwer_mean_short_to_kep(mean::AbstractVector, μ::Real)

Convert Brouwer mean (short-period-only, J2) elements to osculating Keplerian elements.

# Arguments
- `mean::AbstractVector`: Brouwer mean (short) elements `[a, e, i, Ω, ω, M]` (see `brouwer_mean_long_to_kep`)
- `μ`: Gravitational parameter (Earth only)

# Returns
A 6-element vector `[a, e, i, Ω, ω, ν]` of osculating Keplerian elements, where `ν` is the **true anomaly** (rad).

# Example
osc = brouwer_mean_short_to_kep([7000.0, 0.01, pi/6, 0.0, 0.0, pi/3], 398600.4415)

# Notes
- As `brouwer_mean_long_to_kep`, but applies only the first-order J2 short-period correction.
- Inverse of `kep_to_brouwer_mean_short`. See `brouwer_mean_long_to_kep` for the shared contract.
"""
brouwer_mean_short_to_kep(mean::AbstractVector, μ::Real) =
    (o = _brouwer_short_mean_to_osc(mean, μ); [o[1], o[2], o[3], o[4], o[5], _mean_to_true_anom(o[6], o[2])])
"""
    brouwer_mean_long_to_kep(mean::AbstractVector, μ::Real)

Convert Brouwer-Lyddane mean (long-period) elements to osculating Keplerian elements.

# Arguments
- `mean::AbstractVector`: Brouwer mean (long) elements `[a, e, i, Ω, ω, M]`
- `μ`: Gravitational parameter (Earth only)
- `a`: semi-major axis (km)
- `e`: eccentricity
- `i`: inclination (rad)
- `Ω`: right ascension of ascending node (rad)
- `ω`: argument of periapsis (rad)
- `M`: mean anomaly (rad)

# Returns
A 6-element vector `[a, e, i, Ω, ω, ν]` of osculating Keplerian elements, where `ν` is the **true anomaly** (rad).

# Example
osc = brouwer_mean_long_to_kep([7000.0, 0.01, pi/6, 0.0, 0.0, pi/3], 398600.4415)

# Notes
- Angles in radians; distances in km consistent with the Earth `μ`.
- Earth-only: errors if `μ` is not Earth's gravitational parameter.
- Domain: `0 ≤ e < 0.99` and mean periapsis radius `a(1-e) > 3000 km`; warns below 6378 km.
- Applies the secular + long-period + short-period zonal (J2–J5) corrections. Inverse of
  `kep_to_brouwer_mean_long`.
- Differentiable (ForwardDiff) with respect to `mean`.
- Brouwer (1959) / Lyddane (1963) theory; transcribed from GMAT `StateConversionUtil.cpp`.
"""
brouwer_mean_long_to_kep(mean::AbstractVector, μ::Real) =
    (o = _brouwer_long_mean_to_osc(mean, μ); [o[1], o[2], o[3], o[4], o[5], _mean_to_true_anom(o[6], o[2])])

# ============================================================================
# State-type conversion constructors (pivot through KeplerianState)
# ============================================================================
BrouwerMeanLongState(s::AbstractVector{T}, μ::Real) where {T<:Real} =
    (@assert length(s) == 6; BrouwerMeanLongState{T}(s[1], s[2], s[3], s[4], s[5], s[6]))
BrouwerMeanLongState(s::AbstractVector{T}) where {T<:Real} = BrouwerMeanLongState(s, zero(T))
BrouwerMeanShortState(s::AbstractVector{T}, μ::Real) where {T<:Real} =
    (@assert length(s) == 6; BrouwerMeanShortState{T}(s[1], s[2], s[3], s[4], s[5], s[6]))
BrouwerMeanShortState(s::AbstractVector{T}) where {T<:Real} = BrouwerMeanShortState(s, zero(T))

BrouwerMeanLongState(s::BrouwerMeanLongState, μ::Real)  = s
BrouwerMeanShortState(s::BrouwerMeanShortState, μ::Real) = s

BrouwerMeanLongState(kep::KeplerianState, μ::Real) =
    BrouwerMeanLongState(kep_to_brouwer_mean_long(to_vector(kep), μ))
KeplerianState(b::BrouwerMeanLongState, μ::Real) =
    KeplerianState(brouwer_mean_long_to_kep(to_vector(b), μ), μ)
BrouwerMeanShortState(kep::KeplerianState, μ::Real) =
    BrouwerMeanShortState(kep_to_brouwer_mean_short(to_vector(kep), μ))
KeplerianState(b::BrouwerMeanShortState, μ::Real) =
    KeplerianState(brouwer_mean_short_to_kep(to_vector(b), μ), μ)

BrouwerMeanLongState(s::AbstractOrbitState, μ::Real)  = BrouwerMeanLongState(KeplerianState(s, μ), μ)
BrouwerMeanShortState(s::AbstractOrbitState, μ::Real) = BrouwerMeanShortState(KeplerianState(s, μ), μ)
CartesianState(b::BrouwerMeanLongState, μ::Real)  = CartesianState(KeplerianState(b, μ), μ)
CartesianState(b::BrouwerMeanShortState, μ::Real) = CartesianState(KeplerianState(b, μ), μ)
