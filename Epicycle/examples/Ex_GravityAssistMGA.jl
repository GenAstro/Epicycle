# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # Earth-Earth-Venus Gravity Assist
#'
#' Maximize the mass delivered to Venus using an Earth gravity assist. Two
#' `MGAnDSMsPhase`s represent the Earth-to-Earth and Earth-to-Venus legs, with
#' one optional deep-space maneuver on each leg.
#'
#' The patched-conic model joins Keplerian legs at encounters. At the Earth
#' flyby, the incoming excess velocity turns while retaining its magnitude. The
#' solver varies encounter epochs, excess velocities, and boundary masses.
#'
#' The reference solution delivers 91.55 kg after 512.2 days. The paper's Table 4.5 reports a
#' different local optimum at 540 days.

using Epicycle

#' ## Mission Models
#'
#' Define simplified planet, ephemeris, flyby, launch, and arrival models so the
#' example is self-contained. The circular coplanar ephemerides are suitable for
#' demonstrating the transcription, not for mission analysis.

#' ### Planet data


const MU_SUN = 1.32712440018e11   # km³/s²  heliocentric gravitational parameter
const AU     = 1.495978707e8      # km       1 astronomical unit
const G0     = 9.80665e-3         # km/s²    standard gravity

# Store the constants for a body on a circular heliocentric orbit

struct PlanetData
    name      ::Symbol
    mu        ::Float64   # km³/s²  planetary gravitational parameter
    r_orbit   ::Float64   # km      mean heliocentric orbit radius (circular)
    r_surface ::Float64   # km      mean equatorial surface radius
    theta0    ::Float64   # rad     ecliptic longitude at J2000
end

"""Circular heliocentric speed for a planet (km/s)."""
vc_orbit(p::PlanetData) = sqrt(MU_SUN / p.r_orbit)

"""Mean motion for a planet (rad/s)."""
om_orbit(p::PlanetData) = vc_orbit(p) / p.r_orbit

# Define the Earth and Venus approximations
const EARTH   = PlanetData(:earth,   3.986004418e5, 1.000 * AU,  6371.0, 0.0           )
const VENUS   = PlanetData(:venus,   3.24859e5,     0.723 * AU,  6051.8, deg2rad(130.0))

#' ### The ephemeris


"""
    circular_ephemeris(planet::PlanetData) → Function

Returns a closure `eph(t::Float64) → (r, v, a)` representing a body on a
circular heliocentric orbit in the ecliptic plane.

  r : [3] km     heliocentric position
  v : [3] km/s   heliocentric velocity
  a : [3] km/s²  gravitational acceleration (for STM-based propagators)

The planet's ecliptic longitude at J2000 (t = 0) is `planet.theta0`.
"""
function circular_ephemeris(planet::PlanetData)
    vc = vc_orbit(planet)
    om = om_orbit(planet)
    θ0 = planet.theta0
    r0 = planet.r_orbit
    μs = MU_SUN
    return function(t::Float64)
        θ = θ0 + om * t
        r = r0 * [cos(θ),  sin(θ), 0.0]
        v = vc * [-sin(θ), cos(θ), 0.0]
        a = -(μs / r0^2) * [cos(θ), sin(θ), 0.0]
        return r, v, a
    end
end

#' ### Flyby geometry


function flyby_turn_angle(vinf_in::AbstractVector, vinf_out::AbstractVector)
    mag_in  = norm(vinf_in)
    mag_out = norm(vinf_out)
    (mag_in < 1e-12 || mag_out < 1e-12) && return 0.0
    # Keep the acos derivative finite near parallel and antiparallel vectors
    cos_δ = clamp(dot(vinf_in, vinf_out) / (mag_in * mag_out), -1.0 + 1e-10, 1.0 - 1e-10)
    return acos(cos_δ)
end

# Convert flyby turn angle to hyperbolic eccentricity
function flyby_eccentricity(delta::Real)
    sin_half = sin(delta / 2)
    sin_half < 1e-12 && return Inf   # no deflection → straight-line (e → ∞)
    return 1.0 / sin_half
end

# Compute flyby periapsis radius from the excess-velocity turn
function flyby_rp(vinf_in::AbstractVector, vinf_out::AbstractVector,
                  planet::PlanetData)
    vinf_mag = 0.5 * (norm(vinf_in) + norm(vinf_out))
    δ        = flyby_turn_angle(vinf_in, vinf_out)
    return flyby_rp(vinf_mag, δ, planet)
end

function flyby_rp(vinf_mag::Real, delta::Real, planet::PlanetData)
    vinf_mag < 1e-12 && return Inf
    e  = flyby_eccentricity(delta)
    isinf(e) && return Inf
    return planet.mu * (e - 1.0) / vinf_mag^2
end

#' ### Launch vehicle and orbit insertion


# Define the departure and capture burn relations
"""Hyperbolic periapsis speed: sqrt(v∞² + 2μ/r)  [km/s]"""
_v_hyp(vinf_mag, mu, r) = sqrt(vinf_mag^2 + 2*mu/r)

"""Circular orbit speed at radius r: sqrt(μ/r)  [km/s]"""
_v_circ(mu, r) = sqrt(mu / r)

"""ΔV for a hyperbolic departure or capture burn at periapsis  [km/s]"""
_dv_burn(vinf_mag, mu, r) = _v_hyp(vinf_mag, mu, r) - _v_circ(mu, r)

# Model launch mass as a function of departure excess velocity
"""
    GALLOPModel

Parameters for the GALLOP launch vehicle model.

Fields:
  zeta1  — total launch mass budget [kg]
  zeta2  — upper stage structural mass [kg]  (0 if no upper stage)
  Isp    — propulsion system Isp [s]
  g0     — standard gravity [km/s²]
  r_park — departure parking orbit radius [km]
  planet — departure planet (supplies μ for the TMI hyperbolic burn)
"""
struct GALLOPModel
    zeta1  ::Float64
    zeta2  ::Float64
    Isp    ::Float64
    g0     ::Float64
    r_park ::Float64
    planet ::PlanetData
end

"""
    lv_m0(model::GALLOPModel, vinf_dep) → m0 [kg]

Heliocentric departure mass delivered to escape trajectory by the launch vehicle.
"""
function lv_m0(model::GALLOPModel, vinf_dep::AbstractVector)
    ve     = model.Isp * model.g0
    vinf   = norm(vinf_dep)
    dv_tmi = _dv_burn(vinf, model.planet.mu, model.r_park)
    m0star = model.zeta1 * exp(-dv_tmi / ve) - model.zeta2
    return m0star   # simplified: k_lv = k_adp = 0
end

"""
    lv_dm0_dvinf(model::GALLOPModel, vinf_dep) → 1×3 matrix  [kg/(km/s)]

∂m0/∂v∞_dep_i  (row Jacobian for NLP)
"""
function lv_dm0_dvinf(model::GALLOPModel, vinf_dep::AbstractVector)
    ve    = model.Isp * model.g0
    vinf  = norm(vinf_dep)
    vinf < 1e-12 && return zeros(1, 3)
    v_hyp = _v_hyp(vinf, model.planet.mu, model.r_park)
    m0    = lv_m0(model, vinf_dep)
    coeff = -(m0 + model.zeta2) / (ve * v_hyp)
    return reshape(coeff .* vinf_dep, 1, 3)
end

# Model the propulsive capture at Venus
"""
    OrbitInsertionModel

Parameters for a propulsive orbit insertion burn at arrival.

Fields:
  Isp    — propulsion system Isp [s]
  g0     — standard gravity [km/s²]
  r_park — capture parking orbit radius [km]
  planet — arrival planet (supplies μ for the capture hyperbolic burn)
"""
struct OrbitInsertionModel
    Isp    ::Float64
    g0     ::Float64
    r_park ::Float64
    planet ::PlanetData
end

"""
    orbit_insertion_delivered(model, mf_helio, vinf_arr) → m_delivered [kg]

Mass delivered to the capture parking orbit after the insertion burn (Tsiolkovsky).
"""
function orbit_insertion_delivered(model::OrbitInsertionModel,
                                    mf_helio::Real,
                                    vinf_arr::AbstractVector)
    ve     = model.Isp * model.g0
    vinf   = norm(vinf_arr)
    dv_moi = _dv_burn(vinf, model.planet.mu, model.r_park)
    return mf_helio * exp(-dv_moi / ve)
end

"""
    orbit_insertion_dm_dmf(model, vinf_arr) → scalar

∂m_delivered/∂mf_helio
"""
function orbit_insertion_dm_dmf(model::OrbitInsertionModel,
                                 vinf_arr::AbstractVector)
    ve     = model.Isp * model.g0
    vinf   = norm(vinf_arr)
    dv_moi = _dv_burn(vinf, model.planet.mu, model.r_park)
    return exp(-dv_moi / ve)
end

"""
    orbit_insertion_dm_dvinf(model, mf_helio, vinf_arr) → 3-vector  [kg/(km/s)]

∂m_delivered/∂v∞_arr_i
"""
function orbit_insertion_dm_dvinf(model::OrbitInsertionModel,
                                   mf_helio::Real,
                                   vinf_arr::AbstractVector)
    ve    = model.Isp * model.g0
    vinf  = norm(vinf_arr)
    vinf < 1e-12 && return zeros(3)
    v_hyp = _v_hyp(vinf, model.planet.mu, model.r_park)
    m_del = orbit_insertion_delivered(model, mf_helio, vinf_arr)
    return (-m_del / (ve * v_hyp)) .* vinf_arr
end
using LinearAlgebra, Printf, SNOW

#' ## Configuration
#'
#' Define the parking orbits, flyby altitude limit, launch-energy limit, and
#' encounter epochs used to initialize the problem.


const M0_KG    = 1000.0
const ISP      = 320.0
const R_PARK_E = EARTH.r_surface + 300.0
const R_PARK_V = VENUS.r_surface + 300.0
const H_P_MIN  = 200.0
const R_P_MIN  = EARTH.r_surface + H_P_MIN
const C3_MAX   = 20.0

const T_DEP = 12222.0 * 86400.0                 # 2033-Jun-19
const T_FLY = T_DEP + 365.0 * 86400.0           # 2034-Jun-19
const T_ARR = T_FLY + 175.0 * 86400.0           # 2034-Dec-11

const VC_E = vc_orbit(EARTH);  const OM_E = om_orbit(EARTH)
const VC_V = vc_orbit(VENUS);  const OM_V = om_orbit(VENUS)

earth_eph = circular_ephemeris(EARTH)
venus_eph = circular_ephemeris(VENUS)

const lv = GALLOPModel(M0_KG, 0.0, ISP, G0, R_PARK_E, EARTH)
const oi = OrbitInsertionModel(ISP, G0, R_PARK_V, VENUS)

#' ## Build the Two Legs
#'
#' Create one phase from Earth to the Earth flyby and one from the flyby to
#' Venus. Each leg permits one deep-space maneuver.


trans = MGAnDSMs(n_dsm = 1)
model = PropulsionModel(mu = MU_SUN, Isp = ISP, g0 = G0)

phase1 = MGAnDSMsPhase(name = :earth_to_flyby, transcription = trans, model = model,
                       ephemeris_left = earth_eph, ephemeris_right = earth_eph)
phase2 = MGAnDSMsPhase(name = :flyby_to_venus, transcription = trans, model = model,
                       ephemeris_left = earth_eph, ephemeris_right = venus_eph)
phase1._alpha = [0.5, 0.5]          # split the leg at its midpoint
phase2._alpha = [0.5, 0.5]

#' ## Set the Initial Guess
#'
#' Initialize the encounter velocities from the reference solution. Use a radial
#' incoming velocity at the flyby to avoid the singular derivative associated
#' with a 180-degree turn.


r_dep, v_dep, _ = earth_eph(T_DEP)
r_fly, v_fly, _ = earth_eph(T_FLY)
r_arr, v_arr, _ = venus_eph(T_ARR)
vhat(v) = v ./ norm(v)

# Set the excess-velocity guesses at departure, flyby, and arrival
vinf_dep1 =  3.944 .* vhat(v_dep)
vinf_arr1 =  3.962 .* vhat(r_fly)
vinf_dep2 = -3.962 .* vhat(v_fly)
vinf_arr2 = -3.141 .* vhat(v_arr)
m0_guess  = lv_m0(lv, vinf_dep1)

#' ## Define the Optimization Variables
#'
#' Vary the departure and arrival excess velocities, boundary epochs, and
#' boundary masses for each leg. Scale each quantity by its characteristic
#' orbital value.

Vary(departure_vinf, phase1; guess = vinf_dep1, lower_bound = fill(-6.0, 3),
                             upper_bound = fill(6.0, 3), scale = fill(VC_E, 3))
Vary(arrival_vinf,   phase1; guess = vinf_arr1, lower_bound = fill(-6.0, 3),
                             upper_bound = fill(6.0, 3), scale = fill(VC_E, 3))
Vary(initial_time,   phase1; guess = T_DEP, lower_bound = T_DEP - 30 * 86400.0,
                             upper_bound = T_DEP + 30 * 86400.0, scale = 1 / OM_E)
Vary(final_time,     phase1; guess = T_FLY, lower_bound = T_DEP + 200 * 86400.0,
                             upper_bound = T_DEP + 500 * 86400.0, scale = 1 / OM_E)
Vary(initial_mass,   phase1; guess = m0_guess, lower_bound = 10.0,
                             upper_bound = M0_KG, scale = M0_KG)
Vary(final_mass,     phase1; guess = m0_guess, lower_bound = 10.0,
                             upper_bound = M0_KG, scale = M0_KG)

Vary(departure_vinf, phase2; guess = vinf_dep2, lower_bound = fill(-6.0, 3),
                             upper_bound = fill(6.0, 3), scale = fill(VC_E, 3))
Vary(arrival_vinf,   phase2; guess = vinf_arr2, lower_bound = fill(-5.0, 3),
                             upper_bound = fill(5.0, 3), scale = fill(VC_V, 3))
Vary(initial_time,   phase2; guess = T_FLY, lower_bound = T_DEP + 200 * 86400.0,
                             upper_bound = T_DEP + 500 * 86400.0, scale = 1 / OM_E)
Vary(final_time,     phase2; guess = T_ARR, lower_bound = T_DEP + 300 * 86400.0,
                             upper_bound = T_DEP + 800 * 86400.0, scale = 1 / OM_V)
Vary(initial_mass,   phase2; guess = m0_guess, lower_bound = 10.0,
                             upper_bound = M0_KG, scale = M0_KG)
Vary(final_mass,     phase2; guess = m0_guess, lower_bound = 10.0,
                             upper_bound = M0_KG, scale = M0_KG)

#' ## Constrain the Flyby
#'
#' Join the phases at Earth. Epoch and mass remain continuous, incoming and
#' outgoing excess speeds are equal, and periapsis remains above the minimum
#' flyby altitude.


flyby = Link(phase1, phase2; body = EARTH, name = :earth_flyby)

link_epoch(c1, c2)   = [final_time(c1) - initial_time(c2)]
link_mass(c1, c2)    = [final_mass(c1) - initial_mass(c2)]
flyby_speed(c1, c2)  = [dot(arrival_vinf(c1), arrival_vinf(c1)) -
                        dot(departure_vinf(c2), departure_vinf(c2))]
flyby_height(c1, c2) = [flyby_rp(arrival_vinf(c1), departure_vinf(c2), EARTH)]

Constraint(link_epoch,   flyby; equals = 0.0,            name = :epoch)
Constraint(link_mass,    flyby; equals = 0.0,            name = :mass)
Constraint(flyby_speed,  flyby; equals = 0.0,            name = :speed)
Constraint(flyby_height, flyby; lower_bound = R_P_MIN,   name = :height)

#' ## Apply the Leg Constraints
#'
#' Limit launch energy, tie departure mass to launch-vehicle performance, and
#' bound the mass delivered after Venus orbit insertion.


launch_c3(c)     = [dot(departure_vinf(c), departure_vinf(c))]
launch_mass(c)   = [initial_mass(c) - lv_m0(lv, departure_vinf(c))]
delivered_mass(c) = [orbit_insertion_delivered(oi, final_mass(c), arrival_vinf(c))]

# Supply launch-energy and launch-mass partials
@partial(launch_c3, departure_vinf) do c
    reshape(2 .* departure_vinf(c), 1, 3)
end
@partial(launch_mass, initial_mass)   do c; reshape([1.0], 1, 1)                end
@partial(launch_mass, departure_vinf) do c; -lv_dm0_dvinf(lv, departure_vinf(c)) end

# Supply delivered-mass partials for arrival mass and excess velocity
@partial(delivered_mass, final_mass) do c
    reshape([orbit_insertion_dm_dmf(oi, arrival_vinf(c))], 1, 1)
end
@partial(delivered_mass, arrival_vinf) do c
    reshape(orbit_insertion_dm_dvinf(oi, final_mass(c), arrival_vinf(c)), 1, 3)
end

Constraint(launch_c3,      phase1; lower_bound = 0.0, upper_bound = C3_MAX,
                                   scale = C3_MAX)
Constraint(launch_mass,    phase1; equals = 0.0, scale = M0_KG)
Constraint(delivered_mass, phase2; lower_bound = 50.0, upper_bound = M0_KG,
                                   scale = M0_KG)

#' ## Define the Objective
#'
#' Maximize the mass remaining in the Venus parking orbit after the capture burn.


mayer_delivered(c) = orbit_insertion_delivered(oi, final_mass(c), arrival_vinf(c))

# Supply the Mayer-objective gradient as vectors
@partial(mayer_delivered, final_mass) do c
    [orbit_insertion_dm_dmf(oi, arrival_vinf(c))]
end
@partial(mayer_delivered, arrival_vinf) do c
    orbit_insertion_dm_dvinf(oi, final_mass(c), arrival_vinf(c))
end

Objective(mayer_delivered, phase2; sense = Max())

seq = Sequence()
add_sequence!(seq, phase1)
add_sequence!(seq, phase2)

# Check all supplied partials before solving
println("check_partials, phase 1:"); check_partials(phase1)
println("check_partials, phase 2:"); check_partials(phase2)

#' ## Solve the Problem
#'
#' Solve the two-leg trajectory with IPOPT and compare the encounter speeds and
#' duration with the reference solution.

result = solve!(seq; method = Optimize(max_iter = 500, tol = 1e-6, print_level = 5))

@printf("status         : %s
", result.info)
@printf("delivered mass : %.3f kg
",
        orbit_insertion_delivered(oi, phase2.mf_var.value[1], phase2.vinf_arr_var.value))
@printf("departure Vinf : %.4f km/s   (truth 3.944)
", norm(phase1.vinf_dep_var.value))
@printf("flyby Vinf     : %.4f km/s   (truth 3.962)
", norm(phase1.vinf_arr_var.value))
@printf("arrival Vinf   : %.4f km/s   (truth 3.141)
", norm(phase2.vinf_arr_var.value))
@printf("duration       : %.1f days   (truth 540)
",
        (phase2.tf_var.value[1] - phase1.t0_var.value[1]) / 86400.0)
