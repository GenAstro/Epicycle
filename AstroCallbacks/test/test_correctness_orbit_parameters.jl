# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# =============================================================================
# Orbit parameters of the current state (#11, sections 1, 2 and 6): the quantities that read a
# state representation or call an AstroRoutines formula, tested here end to end through a
# Spacecraft, in two frames.
#
# Truth, strongest first:
# - External reference: GMAT R2026a, from truth/orbit_parameters_truth.script, for the mean,
#   eccentric and hyperbolic anomalies, specific energy, the angular momentum vector and
#   magnitude, period, apoapsis and periapsis speed, and C3 (for the excess speed). Six orbits:
#   near-circular LEO, Molniya, polar, retrograde, and a hyperbola outbound and inbound.
# - Analytic: Kepler's equation, vis-viva, orthogonality and the identities between the
#   parameters, for the ones the reference does not report.
# - Differentiability: ForwardDiff against central differences, through a Coordinate subject.
# =============================================================================

using Test
using ForwardDiff
using AstroCallbacks
using AstroModels
using AstroStates
using AstroEpochs
using AstroFrames
using AstroUniverse
using EpicycleBase: label, unit, is_cyclic, is_settable, set_quantity!
using LinearAlgebra: cross, dot, norm
using Logging: with_logger, NullLogger

const _OP_MU    = 398600.4418                       # Earth, as EarthMJ2000Eq carries it
const _OP_EPOCH = Time("2024-01-01T12:00:00", UTC(), ISOT())

const _OP_STATES = Dict(
    :leo        => [7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9],
    :molniya    => [4637.031328726552, 178.53697947901998, -5679.055240387162,
                    6.252424682730313, 6.928411997008257, 2.5730558589825416],
    :polar      => [1731.6965430023222, -2999.386395771405, -7427.270441381914,
                    -3.0611710106575205, 5.302103721115798, -2.5783940176845728],
    :retrograde => [-5967.240072871029, 10335.562987173626, 6.230821664531257e-13,
                    3.6791756283823123, 3.4619800616646623, -2.8389753744567603],
    :hyp_out    => [-3701.281373968369, 7485.70648182266, 4011.578013161604,
                    -9.268115532432553, -1.7246408988602304, 4.13619682261008],
    :hyp_in     => [6729.261014450535, 6167.260026326277, -1584.29588904432,
                    -8.0185180311041, 3.909714301621651, 5.137896317799944])

_op_sat(name) = Spacecraft(state = CartesianState(_OP_STATES[name]), time = _OP_EPOCH,
                           coord_sys = EarthMJ2000Eq)

# GMAT output, angles in degrees as GMAT reports them. `anomaly` is EA on a closed orbit and HA on
# a hyperbola; `period` is C3 on a hyperbola.
const _OP_TRUTH = Dict(
    :leo        => (MA = 122.5442363250099, anomaly = 123.0840766021676, energy = -28.20920520414685,
                    h = [-8610.0, -6780.0, 51920.0], hmag = 53063.98873058828,
                    period = 5909.984957722748, sma = 7065.077497139203,
                    vapo = 7.42722369879468, vperi = 7.59616415181483),
    :molniya    => (MA = 3.137014995390691, anomaly = 11.8269083043028, energy = -7.492489507518812,
                    h = [39806.22008026456, -47439.20578825107, 31010.97447117744],
                    hmag = 69258.16876405635, period = 43175.10828214536,
                    sma = 26599.99999999995, vapo = 1.496373882206735, vperi = 10.01419444246043),
    :polar      => (MA = 204.2134756141651, anomaly = 202.0614423370398, energy = -26.57336278666666,
                    h = [47113.7581845659, 27201.14097039404, 5.456968210637569e-12],
                    hmag = 54402.2819407881, period = 6464.02273990878, sma = 7500.0,
                    vapo = 6.5942159928228, vperi = 8.059597324561201),
    :retrograde => (MA = 76.9807494122247, anomaly = 88.43563434076363, energy = -16.60835174166667,
                    h = [-29342.40880173268, -16940.84762035241, -58684.81760346539],
                    hmag = 67763.39048140965, period = 13082.26221134972, sma = 12000.0,
                    vapo = 4.705791005653448, vperi = 7.058686508480172),
    :hyp_out    => (MA = 7.222077693204831, anomaly = 17.15418710553493, energy = 9.96501104499999,
                    h = [37880.88687557353, -21870.54023475415, 75761.77375114706],
                    hmag = 87482.16093901657, period = 19.93002208999998, sma = -20000.00000000002,
                    vapo = NaN, vperi = 10.93527011737707),
    :hyp_in     => (MA = -7.22207769320485, anomaly = -17.15418710553496, energy = 9.965011044999997,
                    h = [37880.88687557353, -21870.54023475414, 75761.77375114708],
                    hmag = 87482.16093901657, period = 19.93002208999999, sma = -20000.0,
                    vapo = NaN, vperi = 10.93527011737707))

const _OP_CLOSED = (:leo, :molniya, :polar, :retrograde)
const _OP_OPEN   = (:hyp_out, :hyp_in)

# GMAT's MJ2000Eq and MJ2000Ec belong to the IAU 1976/1980 chain, so the comparisons with it run
# with Epicycle's frame theory set to FK5, and restored after. The planets come from the same
# DE440 kernel in both tools (see the truth script).
function _op_fk5(f)
    original = frame_theory()
    set_frame_theory!(FK5())
    try
        f()
    finally
        set_frame_theory!(original)
    end
end

_op_fk5() do
@testset "agree with an external reference (GMAT)" begin
    for name in (_OP_CLOSED..., _OP_OPEN...)
        @testset "$name" begin
            sat, t = _op_sat(name), _OP_TRUTH[name]
            @test rad2deg(mean_anomaly(sat)) ≈ t.MA rtol = 1e-11
            @test specific_energy(sat) ≈ t.energy rtol = 1e-12
            @test angular_momentum(sat) ≈ t.h rtol = 1e-12 atol = 1e-8
            @test angular_momentum_magnitude(sat) ≈ t.hmag rtol = 1e-12
            @test mean_motion(sat) ≈ sqrt(_OP_MU / abs(t.sma)^3) rtol = 1e-11
            @test periapsis_speed(sat) ≈ t.vperi rtol = 1e-11
            if name in _OP_CLOSED
                @test rad2deg(eccentric_anomaly(sat)) ≈ t.anomaly rtol = 1e-11
                @test orbital_period(sat) ≈ t.period rtol = 1e-11
                @test apoapsis_speed(sat) ≈ t.vapo rtol = 1e-11
            else
                @test rad2deg(hyperbolic_anomaly(sat)) ≈ t.anomaly rtol = 1e-11
                @test hyperbolic_excess_speed(sat)^2 ≈ t.period rtol = 1e-11   # C3
            end
        end
    end
end
end

@testset "satisfy the two-body identities" begin
    for name in (_OP_CLOSED..., _OP_OPEN...)
        @testset "$name" begin
            sat  = _op_sat(name)
            x    = _OP_STATES[name]
            r, v = x[1:3], x[4:6]
            e    = eccentricity(sat)
            h    = angular_momentum(sat)

            # Kepler's equation, in both forms
            M = mean_anomaly(sat)
            if e < 1
                E = eccentric_anomaly(sat)
                @test mod(E - e * sin(E), 2π) ≈ M atol = 1e-12
            else
                H = hyperbolic_anomaly(sat)
                @test e * sinh(H) - H ≈ M atol = 1e-12
                @test sign(H) == sign(dot(r, v))     # negative inbound, positive outbound
            end

            # Energy is vis-viva, and the speeds at the apsides follow from it
            a = semi_major_axis(sat)
            @test specific_energy(sat) ≈ -_OP_MU / (2a) rtol = 1e-12
            rp = a * (1 - e)
            @test periapsis_speed(sat) ≈ sqrt(_OP_MU * (2 / rp - 1 / a)) rtol = 1e-12
            if e < 1
                ra = a * (1 + e)
                @test apoapsis_speed(sat) ≈ sqrt(_OP_MU * (2 / ra - 1 / a)) rtol = 1e-12
                @test orbital_period(sat) * mean_motion(sat) ≈ 2π rtol = 1e-14
            end

            # The orbit plane
            n̂ = orbit_normal(sat)
            @test norm(n̂) ≈ 1 atol = 1e-15
            @test n̂ ≈ h ./ norm(h)
            @test abs(dot(n̂, r)) < 1e-12 * norm(r) && abs(dot(n̂, v)) < 1e-12 * norm(v)
            @test orbital_angular_velocity(sat) ≈ h ./ dot(r, r)

            # The eccentricity vector: magnitude e, perpendicular to h, pointing to periapsis
            ev = eccentricity_vector(sat)
            @test norm(ev) ≈ e rtol = 1e-12
            @test abs(dot(ev, h)) < 1e-12 * norm(h)
            @test acos(clamp(dot(ev, r) / (e * norm(r)), -1, 1)) ≈
                  abs(rem2pi(true_anomaly(sat), RoundNearest)) atol = 1e-9

            # Angles built from the elements
            Ω, ω, ν = raan(sat), argument_of_periapsis(sat), true_anomaly(sat)
            @test argument_of_latitude(sat) ≈ mod(ω + ν, 2π) atol = 1e-12
            @test longitude_of_periapsis(sat) ≈ mod(Ω + ω, 2π) atol = 1e-12
            @test descending_node_right_ascension(sat) ≈ mod(Ω + π, 2π) atol = 1e-12
        end
    end
end

@testset "argument of latitude stays defined on a circular orbit" begin
    # Circular and inclined: ω and ν are separately arbitrary, their sum is not. Truth: the angle
    # from the ascending node vector to the position, in the orbit plane.
    vc  = sqrt(_OP_MU / 7000.0)
    i   = deg2rad(40.0)
    x   = [7000.0 * cos(1.0), 7000.0 * sin(1.0) * cos(i), 7000.0 * sin(1.0) * sin(i),
           -vc * sin(1.0), vc * cos(1.0) * cos(i), vc * cos(1.0) * sin(i)]
    sat = Spacecraft(state = CartesianState(x), time = _OP_EPOCH, coord_sys = EarthMJ2000Eq)
    @test eccentricity(sat) < 1e-12
    @test argument_of_latitude(sat) ≈ 1.0 atol = 1e-9
end

@testset "throw where undefined" begin
    ellipse, hyperbola = _op_sat(:leo), _op_sat(:hyp_out)
    @test_throws ArgumentError eccentric_anomaly(hyperbola)
    @test_throws ArgumentError apoapsis_speed(hyperbola)
    @test_throws ArgumentError orbital_period(hyperbola)
    @test_throws ArgumentError hyperbolic_anomaly(ellipse)
    @test_throws ArgumentError hyperbolic_excess_speed(ellipse)

    msg = try orbital_period(hyperbola); "" catch err; sprint(showerror, err) end
    @test occursin("closed orbit", msg) && occursin("specific energy", msg)
end

@testset "read in any coordinate system" begin
    sat = _op_sat(:molniya)
    # A vector quantity is expressed in the frame asked for; its length is not.
    h_ec = angular_momentum(sat, EarthMJ2000Ec)
    @test h_ec ≈ cross(position_vector(sat, EarthMJ2000Ec), velocity_vector(sat, EarthMJ2000Ec))
    @test !(h_ec ≈ angular_momentum(sat))
    @test norm(h_ec) ≈ angular_momentum_magnitude(sat) rtol = 1e-12
    # An angle measured from the frame's equator changes with it.
    @test !(argument_of_latitude(sat, EarthMJ2000Ec) ≈ argument_of_latitude(sat))
    # So does a speed in a rotating frame: the frame's own rotation enters.
    @test !(specific_energy(sat, EarthFixed) ≈ specific_energy(sat))
end

@testset "mean_anomaly! sets M and holds the other elements" begin
    for name in (:molniya, :retrograde, :hyp_out)
        @testset "$name" begin
            sat  = _op_sat(name)
            held = (semi_major_axis(sat), eccentricity(sat), inclination(sat), raan(sat),
                    argument_of_periapsis(sat))
            to   = name === :hyp_out ? -0.8 : 5.0
            mean_anomaly!(sat; to = to)
            @test mean_anomaly(sat) ≈ to atol = 1e-10
            @test all(isapprox.((semi_major_axis(sat), eccentricity(sat), inclination(sat),
                                 raan(sat), argument_of_periapsis(sat)), held; rtol = 1e-10))
        end
    end

    sat = _op_sat(:leo)
    @test is_settable(sat, mean_anomaly)
    set_quantity!(sat, mean_anomaly; to = 2π + 0.25)     # wraps on a closed orbit
    @test mean_anomaly(sat) ≈ 0.25 atol = 1e-10
    mean_anomaly!(sat, EarthMJ2000Ec; to = 1.5)
    @test mean_anomaly(sat, EarthMJ2000Ec) ≈ 1.5 atol = 1e-10

    before = copy(_OP_STATES[:leo])
    fresh  = _op_sat(:leo)
    @test_throws ArgumentError mean_anomaly!(fresh; to = NaN)
    @test position_vector(fresh) ≈ before[1:3]
    @test !is_settable(sat, orbital_period)
end

# --- B-plane (#11, section 6) -------------------------------------------------

# GMAT B-plane, in EarthMJ2000Eq and EarthMJ2000Ec: B·T, B·R, |B| in km, angle in degrees.
const _OP_BPLANE_TRUTH = Dict(
    :hyp_out => (eq = (19013.84571804537, 4740.640358682364, 19595.91794226543, 13.99987568053744),
                 ec = (17625.14972566177, 8564.700645557023, 25.91680101329081)),
    :hyp_in  => (eq = (19013.84571804536, 4740.64035868236, 19595.91794226542, 13.99987568053744),
                 ec = (17625.14972566177, 8564.700645557023, 25.91680101329081)))

_op_fk5() do
@testset "B-plane agrees with an external reference (GMAT)" begin
    for name in _OP_OPEN
        @testset "$name" begin
            sat, t = _op_sat(name), _OP_BPLANE_TRUTH[name]
            @test b_dot_t(sat)     ≈ t.eq[1] rtol = 1e-11
            @test b_dot_r(sat)     ≈ t.eq[2] rtol = 1e-11
            @test b_magnitude(sat) ≈ t.eq[3] rtol = 1e-12
            @test rad2deg(b_plane_angle(sat)) ≈ t.eq[4] atol = 1e-9
            # The components and the angle depend on the frame's Z axis
            @test b_dot_t(sat, EarthMJ2000Ec) ≈ t.ec[1] rtol = 1e-10
            @test b_dot_r(sat, EarthMJ2000Ec) ≈ t.ec[2] rtol = 1e-10
            @test rad2deg(b_plane_angle(sat, EarthMJ2000Ec)) ≈ t.ec[3] atol = 1e-8
        end
    end
end
end

# GMAT about the Moon and Mars, in Luna- and Mars-centred MJ2000Eq frames with the GMs Epicycle
# uses. Angles in degrees; GMAT reports a hyperbola's mean anomaly in degrees too.
const _OP_LUNA = CoordinateSystem(moon, MJ2000Eq())
const _OP_MARS = CoordinateSystem(mars, MJ2000Eq())
const _OP_BODY_TRUTH = (
    (name = "lunar orbit", cs = _OP_LUNA, x = [1900.0, 300.0, -200.0, 0.2, 1.55, 0.4],
     MA = 69.7726056615227, energy = -1.23392752650758, h = [430.0, -800.0, 2885.0],
     period = 7945.927046849652, ra = 2476.483212765156, rp = 1496.84576699797,
     p = 1865.898058252427, tlong = 8.113027814548889),
    (name = "Mars orbit", cs = _OP_MARS, x = [4000.0, -1500.0, 900.0, 1.0, 2.9, 1.1],
     MA = 25.93171952074533, energy = -4.500025825036821,
     h = [-4259.999999999998, -3499.999999999998, 13099.99999999999],
     period = 9966.529681395987, ra = 5205.766311992815, rp = 4311.595783067953,
     p = 4716.676712357893, tlong = 341.5084555550559))
const _OP_FLYBY_TRUTH = (x = [-8000.0, 3000.0, 1500.0, 2.5, -4.0, 1.0], MA = -99.28898203902548,
                         HA = -49.93835373824788, energy = 6.687826984161758,
                         h = [9000.0, 11749.99999999999, 24500.00000000001], rp = 5254.192206971382,
                         p = 19130.13267270009, bt = 6932.778552291419, br = -3631.885995759158,
                         angle = -27.64876754629629)

_op_fk5() do
@testset "agree with an external reference (GMAT) about the Moon and Mars" begin
    for t in _OP_BODY_TRUTH
        @testset "$(t.name)" begin
            sat = Spacecraft(state = CartesianState(t.x), time = _OP_EPOCH, coord_sys = t.cs)
            @test rad2deg(mean_anomaly(sat, t.cs)) ≈ t.MA rtol = 1e-11
            @test specific_energy(sat, t.cs) ≈ t.energy rtol = 1e-12
            @test angular_momentum(sat, t.cs) ≈ t.h rtol = 1e-12
            @test orbital_period(sat, t.cs) ≈ t.period rtol = 1e-12
            @test apoapsis_radius(sat, t.cs) ≈ t.ra rtol = 1e-12
            @test periapsis_radius(sat, t.cs) ≈ t.rp rtol = 1e-12
            @test semilatus_rectum(sat, t.cs) ≈ t.p rtol = 1e-12
            @test rad2deg(true_longitude(sat, t.cs)) ≈ t.tlong rtol = 1e-11
        end
    end
    @testset "Mars flyby" begin
        t = _OP_FLYBY_TRUTH
        sat = Spacecraft(state = CartesianState(t.x), time = _OP_EPOCH, coord_sys = _OP_MARS)
        @test rad2deg(mean_anomaly(sat, _OP_MARS)) ≈ t.MA rtol = 1e-11
        @test rad2deg(hyperbolic_anomaly(sat, _OP_MARS)) ≈ t.HA rtol = 1e-11
        @test specific_energy(sat, _OP_MARS) ≈ t.energy rtol = 1e-12
        @test angular_momentum(sat, _OP_MARS) ≈ t.h rtol = 1e-12
        @test periapsis_radius(sat, _OP_MARS) ≈ t.rp rtol = 1e-12
        @test semilatus_rectum(sat, _OP_MARS) ≈ t.p rtol = 1e-12
        @test b_dot_t(sat, _OP_MARS) ≈ t.bt rtol = 1e-11
        @test b_dot_r(sat, _OP_MARS) ≈ t.br rtol = 1e-11
        @test rem(rad2deg(b_plane_angle(sat, _OP_MARS)) - t.angle, 360, RoundNearest) ≈ 0 atol = 1e-9
    end
end
end

@testset "B-plane satisfies its geometry" begin
    sat = _op_sat(:hyp_out)
    x   = _OP_STATES[:hyp_out]
    h   = norm(cross(x[1:3], x[4:6]))
    # The impact parameter is h/v∞, in any frame
    @test b_magnitude(sat) ≈ h / hyperbolic_excess_speed(sat) rtol = 1e-12
    @test b_magnitude(sat, EarthMJ2000Ec) ≈ b_magnitude(sat) rtol = 1e-12
    @test b_dot_t(sat)^2 + b_dot_r(sat)^2 ≈ b_magnitude(sat)^2 rtol = 1e-12
    @test b_plane_angle(sat) ≈ mod(atan(b_dot_r(sat), b_dot_t(sat)), 2π)
    @test is_cyclic(b_plane_angle) && unit(b_dot_t) == "km" && label(b_plane_angle) == "B-plane angle"
    # GMAT's BVectorAngle is 90° less the B-vector azimuth the incoming asymptote state carries
    azimuth = IncomingAsymptoteState(CartesianState(x), _OP_MU).bpa
    @test rem2pi(b_plane_angle(sat) - (π / 2 - azimuth), RoundNearest) ≈ 0 atol = 1e-12
    @test_throws ArgumentError b_dot_t(_op_sat(:leo))
    # A circular orbit has no asymptote; the error says so instead of passing a NaN on
    circular = Spacecraft(state = CartesianState([7000.0, 0.0, 0.0, 0.0, sqrt(_OP_MU / 7000.0), 0.0]),
                          time = _OP_EPOCH, coord_sys = EarthMJ2000Eq)
    message = try
        with_logger(NullLogger()) do
            b_plane_angle(circular)
        end
        ""
    catch err
        sprint(showerror, err)
    end
    @test occursin("incoming asymptote could not be formed", message)
end

# --- differentiability --------------------------------------------------------

const _OP_H = [1e-3, 1e-3, 1e-3, 1e-6, 1e-6, 1e-6]

_op_subject(x) = Coordinate(x, EarthMJ2000Eq, _OP_EPOCH)

function _op_central(f, x)
    cols = map(eachindex(x)) do j
        δ = zeros(length(x)); δ[j] = _OP_H[j]
        (f(x .+ δ) .- f(x .- δ)) ./ (2 * _OP_H[j])
    end
    return hcat(cols...)
end

@testset "differentiable: ForwardDiff matches central differences" begin
    scalars = (mean_anomaly, eccentric_anomaly, specific_energy, apoapsis_speed, periapsis_speed,
               angular_momentum_magnitude, orbital_period, mean_motion, argument_of_latitude,
               longitude_of_periapsis, descending_node_right_ascension, apoapsis_radius,
               periapsis_radius, semilatus_rectum, c3, true_longitude)
    vectors = (angular_momentum, orbit_normal, orbital_angular_velocity, eccentricity_vector)
    x = _OP_STATES[:polar]
    for q in scalars
        @testset "$(nameof(q))" begin
            f = x -> q(_op_subject(x))
            @test ForwardDiff.gradient(f, x) ≈ vec(_op_central(f, x)) rtol = 1e-6
        end
    end
    for q in vectors
        @testset "$(nameof(q))" begin
            f = x -> q(_op_subject(x))
            @test ForwardDiff.jacobian(f, x) ≈ _op_central(f, x) rtol = 1e-6
        end
    end
    xh = _OP_STATES[:hyp_out]
    for q in (hyperbolic_anomaly, hyperbolic_excess_speed, mean_anomaly, b_dot_t, b_dot_r,
              b_magnitude, b_plane_angle)
        @testset "$(nameof(q)), hyperbola" begin
            f = x -> q(_op_subject(x))
            @test ForwardDiff.gradient(f, xh) ≈ vec(_op_central(f, xh)) rtol = 1e-6
        end
    end
end

@testset "labels, units and wrapping" begin
    @test label(mean_anomaly) == "Mean anomaly"
    @test label(descending_node_right_ascension) == "RA of descending node"
    @test unit(orbital_period) == "s"
    @test unit(specific_energy) == "km²/s²"
    @test unit(angular_momentum) == "km²/s"
    @test unit(mean_motion) == "rad/s"
    @test is_cyclic(argument_of_latitude) && is_cyclic(mean_anomaly)
    @test !is_cyclic(hyperbolic_anomaly) && !is_cyclic(orbital_period)
    @test orbital_period in quantities(Spacecraft) && mean_anomaly in quantities(Coordinate)
end

nothing
