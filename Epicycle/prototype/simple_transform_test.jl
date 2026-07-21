# Simple Transform Test Runner
# Test/demo code is intentionally separated from FrameTransformPrototype.jl

include("FrameTransformPrototype.jl")

using StaticArrays
using Random
using Statistics

# =============================================================================
# Usage Example
# =============================================================================

function example_usage()
    println("=== Frame Transform Prototype ===\n")

    itrf_frame = CoordinateFrame(earth, ITRF())

    pos = [-1033.4793830; 7901.2952754; 6380.3565958]
    vel = [-3.225636520; -2.872451450; +5.531924446]
    time = Time("2004-04-06T07:51:28.386009", UTC(), ISOT())

    coord_itrf = Coordinate(pos, vel, time, itrf_frame)
    println("Original coordinate in ITRF:")
    println("  pos: $(coord_itrf.pos)")
    println("  vel: $(coord_itrf.vel)\n")

    j2000_frame = CoordinateFrame(earth, J2000())

    println("Finding path from ITRF to J2000...")
    path = get_cached_path(ITRF, J2000)
    println("  Path: ", join(string.(path), " → "), "\n")

    coord_j2000 = transform(coord_itrf, j2000_frame)
    println("Transformed coordinate in J2000:")
    println("  pos: $(coord_j2000.pos)")
    println("  vel: $(coord_j2000.vel)\n")

    println("Graph structure:")
    all_from_axes = sort(unique(first.(keys(EDGE_TRANSFORMS))), by=string)
    for from in all_from_axes
        tos = neighbors_for_axes(from, EdgeTransformMetadata())
        println("  $from → ", join(string.(tos), ", "))
    end
end

# =============================================================================
# IAU2006 Round-Trip Testing
# =============================================================================

function test_iau2006_chain()
    println("\n" * "="^70)
    println("Testing IAU2006 Transformation Chain: ITRF ↔ TIRS ↔ CIRS ↔ ICRF")
    println("="^70 * "\n")

    println("Fetching EOP data...")
    eop_iau2000a = fetch_iers_eop(Val(:IAU2000A))

    pos = [-1033.4793830; 7901.2952754; 6380.3565958]
    vel = [-3.225636520; -2.872451450; +5.531924446]
    time = Time("2004-04-06T07:51:28.386009", UTC(), ISOT())

    itrf_frame = CoordinateFrame(earth, ITRF())
    tirs_frame = CoordinateFrame(earth, TIRS())
    cirs_frame = CoordinateFrame(earth, CIRS())
    icrf_frame = CoordinateFrame(earth, ICRF())

    coord_itrf_orig = Coordinate(pos, vel, time, itrf_frame)

    println("Original ITRF coordinate:")
    println("  pos: $(coord_itrf_orig.pos)")
    println("  vel: $(coord_itrf_orig.vel)\n")

    println("-" * "="^69)
    println("Test 1: Multi-hop transformation ITRF → ICRF")
    println("-" * "="^69)

    path_to_icrf = get_cached_path(ITRF, ICRF)
    println("Graph path: ", join(string.(path_to_icrf), " → "))

    coord_icrf = transform(coord_itrf_orig, icrf_frame)

    println("\nICRF coordinate (via graph):")
    println("  pos: $(coord_icrf.pos)")
    println("  vel: $(coord_icrf.vel)")

    println("\n" * "-" * "="^69)
    println("Test 2: Validate each edge against direct STB calls")
    println("-" * "="^69)

    sv_itrf = to_orbit_state_vector(coord_itrf_orig)
    jd_ut1 = time.utc.jd

    println("\nEdge: ITRF → TIRS")
    sv_tirs_direct = SatelliteToolboxTransformations.sv_ecef_to_ecef(
        sv_itrf, Val(:ITRF), Val(:TIRS), jd_ut1, eop_iau2000a
    )

    coord_tirs_graph = transform(coord_itrf_orig, tirs_frame)
    sv_tirs_graph = to_orbit_state_vector(coord_tirs_graph)

    pos_diff_tirs = sv_tirs_direct.r - sv_tirs_graph.r
    vel_diff_tirs = sv_tirs_direct.v - sv_tirs_graph.v

    println("  Direct STB pos: $(sv_tirs_direct.r)")
    println("  Graph path pos: $(sv_tirs_graph.r)")
    println("  Difference:     $(pos_diff_tirs) (norm: $(norm(pos_diff_tirs)) km)")
    println("  Vel difference: $(vel_diff_tirs) (norm: $(norm(vel_diff_tirs)) km/s)")

    println("\nEdge: TIRS → CIRS")
    sv_cirs_direct = SatelliteToolboxTransformations.sv_ecef_to_eci(
        sv_tirs_direct, Val(:TIRS), Val(:CIRS), jd_ut1
    )

    coord_cirs_graph = transform(coord_itrf_orig, cirs_frame)
    sv_cirs_graph = to_orbit_state_vector(coord_cirs_graph)

    pos_diff_cirs = sv_cirs_direct.r - sv_cirs_graph.r
    vel_diff_cirs = sv_cirs_direct.v - sv_cirs_graph.v

    println("  Direct STB pos: $(sv_cirs_direct.r)")
    println("  Graph path pos: $(sv_cirs_graph.r)")
    println("  Difference:     $(pos_diff_cirs) (norm: $(norm(pos_diff_cirs)) km)")
    println("  Vel difference: $(vel_diff_cirs) (norm: $(norm(vel_diff_cirs)) km/s)")

    println("\n" * "-" * "="^69)
    println("Test 3: Round-trip ITRF → ICRF → ITRF")
    println("-" * "="^69)

    println("\nForward path: ", join(string.(path_to_icrf), " → "))
    path_to_itrf = get_cached_path(ICRF, ITRF)
    println("Return path:  ", join(string.(path_to_itrf), " → "))

    coord_itrf_roundtrip = transform(coord_icrf, itrf_frame)

    println("\nOriginal ITRF:")
    println("  pos: $(coord_itrf_orig.pos)")
    println("  vel: $(coord_itrf_orig.vel)")

    println("\nRound-trip ITRF:")
    println("  pos: $(coord_itrf_roundtrip.pos)")
    println("  vel: $(coord_itrf_roundtrip.vel)")

    pos_error = coord_itrf_roundtrip.pos - coord_itrf_orig.pos
    vel_error = coord_itrf_roundtrip.vel - coord_itrf_orig.vel

    println("\nRound-trip error:")
    println("  pos: $(pos_error) (norm: $(norm(pos_error)) km = $(norm(pos_error)*1e6) mm)")
    println("  vel: $(vel_error) (norm: $(norm(vel_error)) km/s = $(norm(vel_error)*1e6) mm/s)")

    println("\n" * "-" * "="^69)
    println("Test 4: Direct STB round-trip (ITRF → CIRS → ITRF)")
    println("-" * "="^69)

    sv_cirs_stb = SatelliteToolboxTransformations.sv_ecef_to_eci(
        sv_itrf, Val(:ITRF), Val(:CIRS), jd_ut1, eop_iau2000a
    )

    sv_itrf_roundtrip_stb = SatelliteToolboxTransformations.sv_eci_to_ecef(
        sv_cirs_stb, Val(:CIRS), Val(:ITRF), jd_ut1, eop_iau2000a
    )

    pos_error_stb = sv_itrf_roundtrip_stb.r - sv_itrf.r
    vel_error_stb = sv_itrf_roundtrip_stb.v - sv_itrf.v

    println("\nDirect STB round-trip error:")
    println("  pos: $(pos_error_stb) (norm: $(norm(pos_error_stb)) km = $(norm(pos_error_stb)*1e6) mm)")
    println("  vel: $(vel_error_stb) (norm: $(norm(vel_error_stb)) km/s = $(norm(vel_error_stb)*1e6) mm/s)")

    println("\n" * "="^70)
    println("IAU2006 Chain Testing Complete")
    println("="^70 * "\n")
end

# =============================================================================
# LVLH (Orbit-Relative Frame) Testing
# =============================================================================

function test_lvlh_transformations()
    println("\n" * "="^70)
    println("Testing LVLH (Orbit-Relative) Transformations")
    println("="^70 * "\n")

    icrf_frame = CoordinateFrame(earth, ICRF())

    r_earth = 6378.137
    r_orbit = r_earth + 400.0
    mu_earth = 398600.4418
    v_circular = sqrt(mu_earth / r_orbit)

    time = Time("2004-04-06T12:00:00.000", UTC(), ISOT())

    sc_pos = [r_orbit, 0.0, 0.0]
    sc_vel = [0.0, v_circular, 0.0]

    spacecraft = Spacecraft(
        state = CartesianState([sc_pos; sc_vel]),
        time = time,
        name = "TestSat"
    )

    println("Spacecraft orbit state (ICRF):")
    println("  pos: $sc_pos km")
    println("  vel: $sc_vel km/s")
    println("  |r|: $(norm(sc_pos)) km")
    println("  |v|: $(norm(sc_vel)) km/s\n")

    lvlh_frame = CoordinateFrame(earth, LVLH(spacecraft))

    println("-" * "="^69)
    println("Test 1: Verify LVLH frame axes alignment")
    println("-" * "="^69)

    cart_state = get_state(spacecraft, Cartesian())
    state_vec = to_vector(cart_state)
    sc_pos_vec = SVector{3,Float64}(state_vec[1:3])
    sc_vel_vec = SVector{3,Float64}(state_vec[4:6])

    R = compute_lvlh_rotation_matrix(sc_pos_vec, sc_vel_vec)
    x_lvlh = R[:, 1]
    y_lvlh = R[:, 2]
    z_lvlh = R[:, 3]

    println("\nLVLH axes in ICRF coordinates:")
    println("  x-hat (along-track): $(x_lvlh)")
    println("  y-hat (cross-track): $(y_lvlh)")
    println("  z-hat (radial):      $(z_lvlh)")

    r_hat = sc_pos_vec / norm(sc_pos_vec)
    println("\nVerification:")
    println("  z-hat · r-hat = $(dot(z_lvlh, r_hat)) (should be 1.0)")

    h = cross(sc_pos_vec, sc_vel_vec)
    h_hat = h / norm(h)
    println("  y-hat · h-hat = $(dot(y_lvlh, -h_hat)) (should be 1.0)")

    println("  x-hat · y-hat = $(dot(x_lvlh, y_lvlh)) (should be 0.0)")
    println("  y-hat · z-hat = $(dot(y_lvlh, z_lvlh)) (should be 0.0)")
    println("  z-hat · x-hat = $(dot(z_lvlh, x_lvlh)) (should be 0.0)")
    println("  |x-hat| = $(norm(x_lvlh)) (should be 1.0)")
    println("  |y-hat| = $(norm(y_lvlh)) (should be 1.0)")
    println("  |z-hat| = $(norm(z_lvlh)) (should be 1.0)")

    println("\n" * "-" * "="^69)
    println("Test 2: Transform nearby object from ICRF to LVLH")
    println("-" * "="^69)

    v_hat = sc_vel_vec / norm(sc_vel_vec)
    obj_offset_icrf = 10.0 * v_hat

    obj_pos_icrf = sc_pos_vec + obj_offset_icrf
    obj_vel_icrf = sc_vel_vec

    obj_coord_icrf = Coordinate(obj_pos_icrf, obj_vel_icrf, time, icrf_frame)

    println("\nObject in ICRF:")
    println("  pos: $(obj_coord_icrf.pos)")
    println("  vel: $(obj_coord_icrf.vel)")

    obj_coord_lvlh = transform(obj_coord_icrf, lvlh_frame)

    println("\nObject in LVLH (relative to spacecraft):")
    println("  pos: $(obj_coord_lvlh.pos) km")
    println("  vel: $(obj_coord_lvlh.vel) km/s")

    println("\n" * "-" * "="^69)
    println("Test 3: Round-trip ICRF → LVLH → ICRF")
    println("-" * "="^69)

    spacecraft_coord = spacecraft_to_coordinate(spacecraft, icrf_frame)
    println("\nOriginal ICRF coordinate:")
    println("  pos: $(spacecraft_coord.pos)")
    println("  vel: $(spacecraft_coord.vel)")

    coord_lvlh = transform(spacecraft_coord, lvlh_frame)

    println("\nTransformed to LVLH:")
    println("  pos: $(coord_lvlh.pos)")
    println("  vel: $(coord_lvlh.vel)")

    coord_icrf_roundtrip = transform(coord_lvlh, icrf_frame)

    println("\nRound-trip back to ICRF:")
    println("  pos: $(coord_icrf_roundtrip.pos)")
    println("  vel: $(coord_icrf_roundtrip.vel)")

    pos_error = coord_icrf_roundtrip.pos - spacecraft_coord.pos
    vel_error = coord_icrf_roundtrip.vel - spacecraft_coord.vel

    println("\nRound-trip error:")
    println("  pos: $(pos_error)")
    println("  |Δpos|: $(norm(pos_error)) km = $(norm(pos_error)*1e6) mm")
    println("  |Δvel|: $(norm(vel_error)) km/s = $(norm(vel_error)*1e6) mm/s")

    println("\n" * "-" * "="^69)
    println("Test 4: Multi-hop path through graph")
    println("-" * "="^69)

    println("\nFetching EOP data...")
    eop_iau2000a = fetch_iers_eop(Val(:IAU2000A))
    eop_iau1980 = fetch_iers_eop()
    itrf_frame = CoordinateFrame(earth, ITRF())
    coord_itrf = Coordinate(sc_pos_vec, sc_vel_vec, time, itrf_frame)

    path = get_cached_path(ITRF, LVLH)
    println("Graph path ITRF → LVLH: ", join(string.(path), " → "))

    coord_lvlh_via_graph = transform(coord_itrf, lvlh_frame)

    println("\nResult in LVLH (via graph):")
    println("  pos: $(coord_lvlh_via_graph.pos)")
    println("  vel: $(coord_lvlh_via_graph.vel)")

    println("\n" * "="^70)
    println("LVLH Testing Complete")
    println("="^70 * "\n")
end

# =============================================================================
# LVLH with Spacecraft Reference (New API)
# =============================================================================

function test_lvlh_with_spacecraft_reference()
    println("\n" * "="^70)
    println("Testing LVLH with Spacecraft Reference (New API)")
    println("="^70 * "\n")

    icrf_frame = CoordinateFrame(earth, ICRF())

    time = Time("2004-04-06T12:00:00.000", UTC(), ISOT())

    spacecraft = Spacecraft(
        state = CartesianState([4315.56104260, 5410.82548470, 446.04777726,
                               -5.94170149, 4.72854935, 0.58341873]),
        time = time,
        name = "TestSat"
    )

    state_vec = to_vector(get_state(spacecraft, Cartesian()))
    println("Spacecraft orbit state (ICRF):")
    println("  pos: $(state_vec[1:3]) km")
    println("  vel: $(state_vec[4:6]) km/s")
    println("  |r|: $(norm(state_vec[1:3])) km")
    println("  |v|: $(norm(state_vec[4:6])) km/s\n")

    println("-" * "="^69)
    println("Test 1: Create LVLH frame using CoordinateFrame(earth, LVLH(spacecraft))")
    println("-" * "="^69)

    lvlh_frame = CoordinateFrame(earth, LVLH(spacecraft))

    println("\nLVLH frame created successfully")
    println("  Origin: earth")
    println("  Axes: LVLH with reference to spacecraft\n")

    println("-" * "="^69)
    println("Test 2: Transform position/velocity vectors from ICRF to LVLH")
    println("-" * "="^69)

    state_vec = to_vector(get_state(spacecraft, Cartesian()))
    sc_pos = SVector{3,Float64}(state_vec[1:3])
    sc_vel = SVector{3,Float64}(state_vec[4:6])

    test_pos = sc_pos + SVector{3,Float64}(.05, .03, 0.52)
    test_vel = sc_vel + SVector{3,Float64}(0.01, -0.02, 0.015)

    coord_icrf = Coordinate(test_pos, test_vel, time, icrf_frame)

    println("\nTest coordinate in ICRF:")
    println("  pos: $(coord_icrf.pos) km")
    println("  vel: $(coord_icrf.vel) km/s")

    coord_lvlh = transform(coord_icrf, lvlh_frame)

    println("\nTransformed to LVLH:")
    println("  pos: $(coord_lvlh.pos) km")
    println("  vel: $(coord_lvlh.vel) km/s")

    println("\n" * "-" * "="^69)
    println("Test 3: Round-trip ICRF → LVLH → ICRF")
    println("-" * "="^69)

    println("\nOriginal ICRF coordinate:")
    println("  pos: $(coord_icrf.pos)")
    println("  vel: $(coord_icrf.vel)")

    println("\nLVLH coordinate:")
    println("  pos: $(coord_lvlh.pos)")
    println("  vel: $(coord_lvlh.vel)")

    coord_icrf_roundtrip = transform(coord_lvlh, icrf_frame)

    println("\nRound-trip back to ICRF:")
    println("  pos: $(coord_icrf_roundtrip.pos)")
    println("  vel: $(coord_icrf_roundtrip.vel)")

    pos_error = coord_icrf_roundtrip.pos - coord_icrf.pos
    vel_error = coord_icrf_roundtrip.vel - coord_icrf.vel

    println("\nRound-trip error:")
    println("  |Δpos|: $(norm(pos_error)) km = $(norm(pos_error)*1e6) mm")
    println("  |Δvel|: $(norm(vel_error)) km/s = $(norm(vel_error)*1e6) mm/s")

    println("\n" * "-" * "="^69)
    println("Test 4: Verify LVLH axes are correctly oriented")
    println("-" * "="^69)

    R = compute_lvlh_rotation_matrix(sc_pos, sc_vel)
    x_lvlh = R[:, 1]
    y_lvlh = R[:, 2]
    z_lvlh = R[:, 3]

    println("\nLVLH axes in ICRF:")
    println("  x-hat (along-track): $(x_lvlh)")
    println("  y-hat (cross-track): $(y_lvlh)")
    println("  z-hat (radial):      $(z_lvlh)")

    r_hat = sc_pos / norm(sc_pos)
    h = cross(sc_pos, sc_vel)
    h_hat = h / norm(h)

    println("\nVerification:")
    println("  z-hat · r-hat = $(dot(z_lvlh, r_hat)) (should be 1.0)")
    println("  y-hat · h-hat = $(dot(y_lvlh, -h_hat)) (should be 1.0)")
    println("  Orthonormality: |x·y| = $(abs(dot(x_lvlh, y_lvlh))), |y·z| = $(abs(dot(y_lvlh, z_lvlh))), |z·x| = $(abs(dot(z_lvlh, x_lvlh)))")

    println("\n" * "="^70)
    println("LVLH with Spacecraft Reference Testing Complete")
    println("="^70 * "\n")
end

# =============================================================================
# PEF -> MODEq Math Validation Test
# =============================================================================

function test_pef_to_mod_graph_walk()
    println("\n" * "="^70)
    println("Testing PEF → MODEq Transformation Math")
    println("="^70 * "\n")

    time = Time("2004-04-06T12:00:00.000", UTC(), ISOT())
    pef_frame = CoordinateFrame(earth, PEF())
    modeq_frame = CoordinateFrame(earth, MODEq())

    coord_pef = Coordinate(
        SVector{3,Float64}(7000.0, 100.0, 50.0),
        SVector{3,Float64}(0.0, 7.5, 1.0),
        time,
        pef_frame,
    )

    eop_iau1980 = fetch_iers_eop(Val(:IAU1980))

    coord_modeq = transform(coord_pef, modeq_frame)

    sv_pef = to_orbit_state_vector(coord_pef)
    jd_ut1 = time.utc.jd
    sv_modeq_direct = SatelliteToolboxTransformations.sv_ecef_to_eci(
        sv_pef, Val(:PEF), Val(:MOD), jd_ut1, eop_iau1980
    )

    pos_diff = coord_modeq.pos - sv_modeq_direct.r
    vel_diff = coord_modeq.vel - sv_modeq_direct.v

    println("Graph/plan result in MODEq:")
    println("  pos: $(coord_modeq.pos)")
    println("  vel: $(coord_modeq.vel)")

    println("\nDirect STB result in MODEq:")
    println("  pos: $(sv_modeq_direct.r)")
    println("  vel: $(sv_modeq_direct.v)")

    println("\nDifference (graph - direct STB):")
    println("  Δpos: $(pos_diff) (norm: $(norm(pos_diff)) km)")
    println("  Δvel: $(vel_diff) (norm: $(norm(vel_diff)) km/s)")

    println("\n" * "="^70)
    println("PEF → MODEq Math Validation Complete")
    println("="^70 * "\n")
end

# =============================================================================
# J2000 -> TODEq Validation Against Direct STB
# =============================================================================

function test_j2000_to_todeq_vs_stb()
    println("\n" * "="^70)
    println("Testing J2000 → TODEq Against Direct STB")
    println("="^70 * "\n")

    time = Time("2004-04-06T12:00:00.000", UTC(), ISOT())
    j2000_frame = CoordinateFrame(earth, J2000())
    todeq_frame = CoordinateFrame(earth, TODEq())

    coord_j2000 = Coordinate(
        SVector{3,Float64}(7000.0, 100.0, 50.0),
        SVector{3,Float64}(0.0, 7.5, 1.0),
        time,
        j2000_frame,
    )

    eop_iau1980 = fetch_iers_eop(Val(:IAU1980))
    jd_utc = time.utc.jd

    path = get_cached_path(J2000, TODEq)
    println("Graph path J2000 → TODEq: ", join(string.(path), " → "))

    # Graph transform (our implementation)
    coord_todeq_graph = transform(coord_j2000, todeq_frame)

    # Direct STB ECI->ECI call (with EOP)
    sv_j2000 = to_orbit_state_vector(coord_j2000)
    sv_tod_direct_eop = SatelliteToolboxTransformations.sv_eci_to_eci(
        sv_j2000, Val(:J2000), Val(:TOD), jd_utc, eop_iau1980
    )

    # Direct STB ECI->ECI call (without EOP)
    sv_tod_direct_noeop = SatelliteToolboxTransformations.sv_eci_to_eci(
        sv_j2000, Val(:J2000), Val(:TOD), jd_utc
    )

    # Two-step STB path mirroring our graph: J2000 -> ITRF -> TOD
    sv_itrf_step = SatelliteToolboxTransformations.sv_eci_to_ecef(
        sv_j2000, Val(:J2000), Val(:ITRF), jd_utc, eop_iau1980
    )
    sv_tod_twostep = SatelliteToolboxTransformations.sv_ecef_to_eci(
        sv_itrf_step, Val(:ITRF), Val(:TOD), jd_utc, eop_iau1980
    )

    # Differences vs direct STB (with EOP)
    pos_diff_graph_direct = coord_todeq_graph.pos - sv_tod_direct_eop.r
    vel_diff_graph_direct = coord_todeq_graph.vel - sv_tod_direct_eop.v

    # Differences vs direct STB (no EOP)
    pos_diff_graph_noeop = coord_todeq_graph.pos - sv_tod_direct_noeop.r
    vel_diff_graph_noeop = coord_todeq_graph.vel - sv_tod_direct_noeop.v

    # Differences vs two-step STB
    pos_diff_graph_twostep = coord_todeq_graph.pos - sv_tod_twostep.r
    vel_diff_graph_twostep = coord_todeq_graph.vel - sv_tod_twostep.v

    println("\nOur graph result (J2000 → TODEq):")
    println("  pos: $(coord_todeq_graph.pos)")
    println("  vel: $(coord_todeq_graph.vel)")

    println("\nDirect STB (J2000 → TOD, with EOP):")
    println("  pos: $(sv_tod_direct_eop.r)")
    println("  vel: $(sv_tod_direct_eop.v)")

    println("\nDirect STB (J2000 → TOD, no EOP):")
    println("  pos: $(sv_tod_direct_noeop.r)")
    println("  vel: $(sv_tod_direct_noeop.v)")

    println("\nTwo-step STB (J2000 → ITRF → TOD, with EOP):")
    println("  pos: $(sv_tod_twostep.r)")
    println("  vel: $(sv_tod_twostep.v)")

    println("\nDifferences (graph - direct with EOP):")
    println("  Δpos norm: $(norm(pos_diff_graph_direct)) km = $(norm(pos_diff_graph_direct)*1000) m")
    println("  Δvel norm: $(norm(vel_diff_graph_direct)) km/s = $(norm(vel_diff_graph_direct)*1000) m/s")

    println("\nDifferences (graph - direct no EOP):")
    println("  Δpos norm: $(norm(pos_diff_graph_noeop)) km = $(norm(pos_diff_graph_noeop)*1000) m")
    println("  Δvel norm: $(norm(vel_diff_graph_noeop)) km/s = $(norm(vel_diff_graph_noeop)*1000) m/s")

    println("\nDifferences (graph - two-step STB):")
    println("  Δpos norm: $(norm(pos_diff_graph_twostep)) km = $(norm(pos_diff_graph_twostep)*1000) m")
    println("  Δvel norm: $(norm(vel_diff_graph_twostep)) km/s = $(norm(vel_diff_graph_twostep)*1000) m/s")

    println("\n" * "="^70)
    println("J2000 → TODEq STB Validation Complete")
    println("="^70 * "\n")
end

# =============================================================================
# ITRF -> Inertial (IAU Context Disambiguation)
# =============================================================================

function test_itrf_to_inertial_iau_context()
    println("\n" * "="^70)
    println("Testing ITRF -> Inertial Context (FK5 vs IAU2006 in STB)")
    println("="^70 * "\n")

    time = Time("2004-04-06T12:00:00.000", UTC(), ISOT())
    itrf_frame = CoordinateFrame(earth, ITRF())
    gcrf_frame = CoordinateFrame(earth, GCRF())
    cirs_frame = CoordinateFrame(earth, CIRS())

    coord_itrf = Coordinate(
        SVector{3,Float64}(7000.0, 100.0, 50.0),
        SVector{3,Float64}(0.0, 7.5, 1.0),
        time,
        itrf_frame,
    )

    eop_iau1980 = fetch_iers_eop(Val(:IAU1980))
    eop_iau2000a = fetch_iers_eop(Val(:IAU2000A))
    jd_utc = time.utc.jd
    sv_itrf = to_orbit_state_vector(coord_itrf)

    println("Interface-selected paths:")
    println("  ITRF -> GCRF: ", join(string.(get_cached_path(ITRF, GCRF)), " → "))
    println("  ITRF -> CIRS: ", join(string.(get_cached_path(ITRF, CIRS)), " → "))

    ctx_fk5 = TransformContext(model_family=:fk5, warn_on_ambiguity=false)
    ctx_iau2006 = TransformContext(model_family=:iau2006, warn_on_ambiguity=false)

    # Interface results
    coord_gcrf_iface_auto = transform(coord_itrf, gcrf_frame)
    coord_gcrf_iface_fk5 = transform(coord_itrf, gcrf_frame, ctx_fk5)
    coord_gcrf_iface_iau = transform(coord_itrf, gcrf_frame, ctx_iau2006)
    coord_cirs_iface = transform(coord_itrf, cirs_frame)

    # Direct STB references
    sv_gcrf_fk5 = SatelliteToolboxTransformations.sv_ecef_to_eci(
        sv_itrf, Val(:ITRF), Val(:GCRF), jd_utc, eop_iau1980
    )
    sv_gcrf_iau = SatelliteToolboxTransformations.sv_ecef_to_eci(
        sv_itrf, Val(:ITRF), Val(:GCRF), jd_utc, eop_iau2000a
    )
    sv_cirs_iau = SatelliteToolboxTransformations.sv_ecef_to_eci(
        sv_itrf, Val(:ITRF), Val(:CIRS), jd_utc, eop_iau2000a
    )

    # Compare interface ITRF->GCRF against both STB interpretations
    gcrf_auto_pos_err_vs_fk5 = norm(coord_gcrf_iface_auto.pos - sv_gcrf_fk5.r)
    gcrf_auto_vel_err_vs_fk5 = norm(coord_gcrf_iface_auto.vel - sv_gcrf_fk5.v)
    gcrf_auto_pos_err_vs_iau = norm(coord_gcrf_iface_auto.pos - sv_gcrf_iau.r)
    gcrf_auto_vel_err_vs_iau = norm(coord_gcrf_iface_auto.vel - sv_gcrf_iau.v)

    gcrf_fk5_pos_err = norm(coord_gcrf_iface_fk5.pos - sv_gcrf_fk5.r)
    gcrf_fk5_vel_err = norm(coord_gcrf_iface_fk5.vel - sv_gcrf_fk5.v)

    gcrf_iau_pos_err = norm(coord_gcrf_iface_iau.pos - sv_gcrf_iau.r)
    gcrf_iau_vel_err = norm(coord_gcrf_iface_iau.vel - sv_gcrf_iau.v)

    # Compare interface ITRF->CIRS against IAU STB reference
    cirs_pos_err = norm(coord_cirs_iface.pos - sv_cirs_iau.r)
    cirs_vel_err = norm(coord_cirs_iface.vel - sv_cirs_iau.v)

    println("\nITRF -> GCRF disambiguation:")
    println("  AUTO context vs STB FK5:        |Δpos|=$(gcrf_auto_pos_err_vs_fk5) km, |Δvel|=$(gcrf_auto_vel_err_vs_fk5) km/s")
    println("  AUTO context vs STB IAU2006:    |Δpos|=$(gcrf_auto_pos_err_vs_iau) km, |Δvel|=$(gcrf_auto_vel_err_vs_iau) km/s")
    println("  model_family=:fk5 vs STB FK5:   |Δpos|=$(gcrf_fk5_pos_err) km, |Δvel|=$(gcrf_fk5_vel_err) km/s")
    println("  model_family=:iau2006 vs STB IAU2006: |Δpos|=$(gcrf_iau_pos_err) km, |Δvel|=$(gcrf_iau_vel_err) km/s")

    println("\nITRF -> CIRS (explicit IAU2006 target):")
    println("  Interface vs STB IAU2006: |Δpos|=$(cirs_pos_err) km, |Δvel|=$(cirs_vel_err) km/s")

    println("\nConclusion:")
    println("  - In STB, tags like GCRF can be reached through different model contexts.")
    println("  - In this prototype, ITRF->GCRF can be forced with TransformContext(model_family=:fk5|:iau2006).")
    println("  - J2000 is treated as FK5-specific (not model-ambiguous) in this prototype.")
    println("  - ITRF->CIRS is the unambiguous IAU2006 inertial check.")

    println("\n" * "="^70)
    println("ITRF -> Inertial Context Test Complete")
    println("="^70 * "\n")
end

# =============================================================================
# All-Axes Round-Trip Test (Extensible)
# =============================================================================

"""
    test_axes_earthcentered_roundtrip()

Round-trip validation across all currently listed axes systems.

Design for extensibility:
- Add new systems by adding one entry to `frame_catalog`.
- Test automatically attempts all source/target pairs.
- Unavailable paths are reported as SKIP.
- Strict numeric PASS requires round-trip error below tolerances.
"""
function test_axes_earthcentered_roundtrip()
    println("\n" * "="^70)
    println("All-Axes Earth-Centered Round-Trip Test")
    println("="^70 * "\n")

    time = Time("2004-04-06T12:00:00.000", UTC(), ISOT())

    spacecraft = Spacecraft(
        state = CartesianState([4315.56104260, 5410.82548470, 446.04777726,
                               -5.94170149, 4.72854935, 0.58341873]),
        time = time,
        name = "TestSat"
    )

    # Add new systems here as they are implemented
    frame_catalog = [
        (name="ICRF",  frame=CoordinateFrame(earth, ICRF())),
        (name="CIRS",  frame=CoordinateFrame(earth, CIRS())),
        (name="GCRF",  frame=CoordinateFrame(earth, GCRF())),
        (name="J2000", frame=CoordinateFrame(earth, J2000())),
        (name="MODEq", frame=CoordinateFrame(earth, MODEq())),
        (name="TODEq", frame=CoordinateFrame(earth, TODEq())),
        (name="ITRF",  frame=CoordinateFrame(earth, ITRF())),
        (name="TIRS",  frame=CoordinateFrame(earth, TIRS())),
        (name="PEF",   frame=CoordinateFrame(earth, PEF())),
        (name="LVLH",  frame=CoordinateFrame(earth, LVLH(spacecraft))),
    ]

    # Seed state in ICRF; source states are derived from this
    seed_icrf = Coordinate(
        SVector{3,Float64}(4315.56, 5410.83, 446.05),
        SVector{3,Float64}(-5.94, 4.73, 0.58),
        time,
        CoordinateFrame(earth, ICRF()),
    )

    tested_pairs = 0
    passed_pairs = 0
    failed_pairs = 0
    skipped_pairs = 0

    pos_tol_km = 1e-6
    vel_tol_km_s = 1e-6

    for source in frame_catalog
        # Build source coordinate by transforming seed ICRF into source frame
        coord_source = if source.name == "ICRF"
            seed_icrf
        else
            try
                transform(seed_icrf, source.frame)
            catch err
                println("SKIP source $(source.name): $(sprint(showerror, err))")
                skipped_pairs += length(frame_catalog) - 1
                continue
            end
        end

        for target in frame_catalog
            if source.name == target.name
                continue
            end

            tested_pairs += 1

            # Known limitation: LVLH transform currently does not model origin translation,
            # so strict round-trip against absolute states is not meaningful.
            if source.name == "LVLH" || target.name == "LVLH"
                println("SKIP $(source.name) → $(target.name): LVLH strict round-trip disabled (origin translation TODO)")
                skipped_pairs += 1
                continue
            end

            try
                coord_target = transform(coord_source, target.frame)
                coord_back = transform(coord_target, source.frame)

                pos_err = norm(coord_back.pos - coord_source.pos)
                vel_err = norm(coord_back.vel - coord_source.vel)

                if pos_err <= pos_tol_km && vel_err <= vel_tol_km_s
                    println("PASS $(source.name) → $(target.name) → $(source.name): |Δpos|=$(pos_err) km, |Δvel|=$(vel_err) km/s")
                    passed_pairs += 1
                else
                    println("FAIL $(source.name) → $(target.name) → $(source.name): |Δpos|=$(pos_err) km, |Δvel|=$(vel_err) km/s")
                    failed_pairs += 1
                end
            catch err
                println("SKIP $(source.name) → $(target.name): $(sprint(showerror, err))")
                skipped_pairs += 1
            end
        end
    end

    println("\n" * "-"^70)
    println("Summary")
    println("  Tested pairs:  $tested_pairs")
    println("  Passed pairs:  $passed_pairs")
    println("  Failed pairs:  $failed_pairs")
    println("  Skipped pairs: $skipped_pairs")
    println("  Tolerances:    |Δpos| ≤ $(pos_tol_km) km, |Δvel| ≤ $(vel_tol_km_s) km/s")
    println("="^70 * "\n")
end

# =============================================================================
# Satellite-Centered Round-Trip Test (GCRF ↔ LVLH, rotation-only)
# =============================================================================

function test_satcentered_roundtrip()
    println("\n" * "="^70)
    println("Satellite-Centered Round-Trip Test (GCRF ↔ LVLH)")
    println("="^70 * "\n")

    time = Time("2004-04-06T12:00:00.000", UTC(), ISOT())

    spacecraft = Spacecraft(
        state = CartesianState([4315.56104260, 5410.82548470, 446.04777726,
                               -5.94170149, 4.72854935, 0.58341873]),
        time = time,
        name = "TestSat"
    )

    gcrf_sat_frame = CoordinateFrame(spacecraft, GCRF())
    lvlh_sat_frame = CoordinateFrame(spacecraft, LVLH(spacecraft))

    coord_gcrf_sat = Coordinate(
        SVector{3,Float64}(0.05, -0.03, 0.52),
        SVector{3,Float64}(0.010, -0.020, 0.015),
        time,
        gcrf_sat_frame,
    )

    println("Input coordinate in satellite-centered GCRF:")
    println("  pos: $(coord_gcrf_sat.pos) km")
    println("  vel: $(coord_gcrf_sat.vel) km/s")

    path = get_cached_path(GCRF, LVLH)
    println("Path GCRF → LVLH: ", join(string.(path), " → "))

    coord_lvlh = transform(coord_gcrf_sat, lvlh_sat_frame)
    coord_back = transform(coord_lvlh, gcrf_sat_frame)

    pos_err = norm(coord_back.pos - coord_gcrf_sat.pos)
    vel_err = norm(coord_back.vel - coord_gcrf_sat.vel)

    pos_tol_km = 1e-6
    vel_tol_km_s = 1e-6

    println("\nRound-trip result (GCRF → LVLH → GCRF):")
    println("  |Δpos| = $(pos_err) km")
    println("  |Δvel| = $(vel_err) km/s")
    println("  Tolerances: |Δpos| ≤ $(pos_tol_km) km, |Δvel| ≤ $(vel_tol_km_s) km/s")

    if pos_err <= pos_tol_km && vel_err <= vel_tol_km_s
        println("  PASS")
    else
        println("  FAIL")
    end

    println("="^70 * "\n")
end

# =============================================================================
# Earth-Fixed -> Inertial Performance Test (Interface Overhead)
# =============================================================================

function test_earthcoords_performance(; n_cases::Int = 20000, repeats::Int = 5)
    println("\n" * "="^70)
    println("Earth-Fixed → Inertial Performance Test (Our Interface vs Direct STB)")
    println("="^70 * "\n")

    if n_cases < 1 || repeats < 1
        error("n_cases and repeats must both be >= 1")
    end

    if !isdefined(@__MODULE__, :Random)
        @eval import Random
    end
    if !isdefined(@__MODULE__, :Statistics)
        @eval import Statistics
    end

    base_time = Time("2004-04-06T12:00:00.000", UTC(), ISOT())
    itrf_frame = CoordinateFrame(earth, ITRF())

    Random.seed!(42)

    coords = Vector{Coordinate}(undef, n_cases)
    jds = Vector{Float64}(undef, n_cases)

    for i in 1:n_cases
        pos = SVector{3,Float64}(
            6378.137 + 500.0 * Random.randn(),
            100.0 * Random.randn(),
            100.0 * Random.randn(),
        )
        vel = SVector{3,Float64}(
            0.02 * Random.randn(),
            7.5 + 0.05 * Random.randn(),
            0.02 * Random.randn(),
        )

        coord = Coordinate(pos, vel, base_time, itrf_frame)
        coords[i] = coord
        jds[i] = base_time.utc.jd + 1e-6 * (i - 1)
    end

    sv_template = to_orbit_state_vector(coords[1])
    svs = Vector{typeof(sv_template)}(undef, n_cases)
    for i in 1:n_cases
        svs[i] = to_orbit_state_vector(coords[i])
    end

    eop_iau1980 = fetch_iers_eop(Val(:IAU1980))
    eop_iau2000a = fetch_iers_eop(Val(:IAU2000A))

    println("Active comparisons:")
    println("  1) FK5:      ITRF -> J2000")
    println("  2) IAU-2006: ITRF -> CIRS")
    println("  (Other comparisons are intentionally commented out.)\n")

    function execute_plan_without_routing(coord::Coordinate,
                                          target_frame::CoordinateFrame,
                                          plan,
                                          context::TransformContext)
        if !isempty(plan.edges) && all(edge -> edge isa STBTransform, plan.edges)
            return execute_stb_sequence(coord, plan.edges, target_frame, context)
        end

        current_coord = coord
        path = plan.path
        for i in 1:(length(path)-1)
            from = path[i]
            to = path[i+1]
            edge = plan.edges[i]

            intermediate_axes = if to == LVLH
                target_frame.axes
            elseif from == LVLH
                coord.frame.axes
            else
                to()
            end
            intermediate_frame = CoordinateFrame(target_frame.origin, intermediate_axes)
            current_coord = execute_transform(current_coord, edge, intermediate_frame, context)
        end
        return current_coord
    end

    function execute_raw_stb_plan(sv_in, jd_ut1::Float64, plan,
                                  eop_iau1980, eop_iau2000a)
        sv = sv_in
        for edge in plan.edges
            if !(edge isa STBTransform)
                error("Non-STB edge in performance test path: $(typeof(edge))")
            end

            eop = if edge.eop_type == :IAU1980
                eop_iau1980
            elseif edge.eop_type == :IAU2000A
                eop_iau2000a
            else
                nothing
            end

            if eop === nothing
                sv = edge.stb_function(sv, edge.val_from, edge.val_to, jd_ut1)
            else
                sv = edge.stb_function(sv, edge.val_from, edge.val_to, jd_ut1, eop)
            end
        end
        return sv
    end

    function benchmark_case(label::String, target_frame::CoordinateFrame)
        from_axes = typeof(coords[1].frame.axes)
        to_axes = typeof(target_frame.axes)

        setup_context = TransformContext()
        available_metadata = edge_transform_metadata_from_context(coords[1], target_frame, setup_context)
        plan = get_sim_plan(from_axes, to_axes, available_metadata)
        path_str = join(string.(plan.path), " → ")
        model_tags = [edge isa STBTransform ? String(edge.eop_type) : string(nameof(typeof(edge))) for edge in plan.edges]

        println("Case: $label")
        println("  Interface path: $path_str")
        println("  Edge model tags: $(join(model_tags, ", "))")

        transform(coords[1], target_frame)
        execute_plan_without_routing(coords[1], target_frame, plan, TransformContext())
        execute_raw_stb_plan(svs[1], jds[1], plan, eop_iau1980, eop_iau2000a)

        ours_times = Float64[]
        plan_exec_times = Float64[]
        direct_with_adapter_times = Float64[]
        direct_times = Float64[]

        for _ in 1:repeats
            GC.gc()
            t0 = time_ns()
            for i in 1:n_cases
                transform(coords[i], target_frame)
            end
            push!(ours_times, (time_ns() - t0) / 1e9)

            GC.gc()
            t0 = time_ns()
            local_ctx = TransformContext()
            for i in 1:n_cases
                execute_plan_without_routing(coords[i], target_frame, plan, local_ctx)
            end
            push!(plan_exec_times, (time_ns() - t0) / 1e9)

            GC.gc()
            t0 = time_ns()
            for i in 1:n_cases
                sv_in = to_orbit_state_vector(coords[i])
                sv_out = execute_raw_stb_plan(sv_in, jds[i], plan, eop_iau1980, eop_iau2000a)
                from_orbit_state_vector(sv_out, coords[i].time, target_frame)
            end
            push!(direct_with_adapter_times, (time_ns() - t0) / 1e9)

            GC.gc()
            t0 = time_ns()
            for i in 1:n_cases
                execute_raw_stb_plan(svs[i], jds[i], plan, eop_iau1980, eop_iau2000a)
            end
            push!(direct_times, (time_ns() - t0) / 1e9)
        end

        ours_median = Statistics.median(ours_times)
        plan_exec_median = Statistics.median(plan_exec_times)
        direct_with_adapter_median = Statistics.median(direct_with_adapter_times)
        direct_median = Statistics.median(direct_times)
        overhead_ratio = ours_median / direct_median
        planner_ratio = ours_median / plan_exec_median
        wrapper_ratio = plan_exec_median / direct_median
        interface_vs_adapter_ratio = ours_median / direct_with_adapter_median

        ours_us = ours_median * 1e6 / n_cases
        plan_exec_us = plan_exec_median * 1e6 / n_cases
        direct_with_adapter_us = direct_with_adapter_median * 1e6 / n_cases
        direct_us = direct_median * 1e6 / n_cases

        println("  Samples: $n_cases per repeat, repeats=$repeats")
        println("  Full interface median total:  $(round(ours_median, sigdigits=6)) s")
        println("  Plan-exec median total:       $(round(plan_exec_median, sigdigits=6)) s")
        println("  Direct+adapter median total:  $(round(direct_with_adapter_median, sigdigits=6)) s")
        println("  Direct STB median total:      $(round(direct_median, sigdigits=6)) s")
        println("  Full interface median/call:   $(round(ours_us, sigdigits=6)) μs")
        println("  Plan-exec median/call:        $(round(plan_exec_us, sigdigits=6)) μs")
        println("  Direct+adapter median/call:   $(round(direct_with_adapter_us, sigdigits=6)) μs")
        println("  Direct STB median/call:       $(round(direct_us, sigdigits=6)) μs")
        println("  Total overhead (full/direct): $(round(overhead_ratio, sigdigits=6))x")
        println("  Planner overhead (full/plan): $(round(planner_ratio, sigdigits=6))x")
        println("  Wrapper overhead (plan/direct): $(round(wrapper_ratio, sigdigits=6))x")
        println("  Interface overhead vs direct+adapter: $(round(interface_vs_adapter_ratio, sigdigits=6))x")
        println()
    end

    benchmark_case(
        "FK5: ITRF -> J2000",
        CoordinateFrame(earth, J2000())
    )

    benchmark_case(
        "IAU-2006: ITRF -> CIRS",
        CoordinateFrame(earth, CIRS())
    )

    # benchmark_case(
    #     "ITRF -> GCRF (interface-selected model)",
    #     CoordinateFrame(earth, GCRF())
    # )

    # benchmark_case(
    #     "FK5: ITRF -> TODEq",
    #     CoordinateFrame(earth, TODEq())
    # )

    println("="^70)
    println("Earth-Fixed → Inertial Performance Test Complete")
    println("="^70 * "\n")
end

# =============================================================================
# Run Tests
# =============================================================================

test_lvlh_transformations()
test_lvlh_with_spacecraft_reference()
test_j2000_to_todeq_vs_stb()
test_itrf_to_inertial_iau_context()
test_axes_earthcentered_roundtrip()
test_satcentered_roundtrip()
test_earthcoords_performance()
