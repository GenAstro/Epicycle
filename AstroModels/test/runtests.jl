using Test
using ForwardDiff
using Logging

using EpicycleBase
using AstroStates
using AstroEpochs
using AstroFrames
using AstroUniverse
using AstroModels

@testset "AstroModels.Spacecraft Constructor Tests" begin
    state0 = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03])
    t0     = Time("2015-09-21T12:23:12", TAI(), ISOT())
    name0  = "SC-001"
    mass0  = 1000.0

    sc = Spacecraft(state=state0, time=t0, mass=mass0, name=name0)

    @test sc.state isa OrbitState
    @test to_vector(sc.state) == to_vector(state0)
    @test sc.time == t0
    @test sc.mass == mass0
    @test sc.name == name0

        # Coordinate system on constructor (Moon, ICRF)
    cs_moon_icrf = CoordinateSystem(moon, ICRFAxes())
    sc_cs = Spacecraft(state=state0, time=t0, mass=mass0, name=name0, coord_sys=cs_moon_icrf)

    @test sc_cs.coord_sys isa CoordinateSystem{<:CelestialBody, ICRFAxes}
    @test sc_cs.coord_sys.origin === moon
    @test sc_cs.coord_sys.axes isa ICRFAxes
end

@testset "set_posvel! and to_posvel tests" begin
    state0 = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03])
    t0     = Time("2015-09-21T12:23:12", TAI(), ISOT())
    name0  = "SC-001"
    mass0  = 1000.0

    sc = Spacecraft(state=state0, time=t0, mass=mass0, name=name0)

    new_pv = [7050.0, 0.0, 0.0, 0.0, 7.6, 0.0]
   set_posvel!(sc, new_pv)          # mutates in place; no return value
   # Verify pos/vel updated
   @test to_posvel(sc) == new_pv
   # Verify conversion to Keplerian matches (allowing FP tolerance)
   kepstate = KeplerianState(CartesianState(new_pv), earth.mu)
   @test isapprox(
       to_vector(get_state(sc, Keplerian())),
       to_vector(kepstate);
       rtol=1e-12, atol=1e-12,
   )

   # Test get_state asking for same type uses shortcut
   sc = Spacecraft(state=CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]))
   st = get_state(sc, Cartesian() )
   @test sc.state.state == st.posvel
end

@testset "AstroModels.Spacecraft show" begin
    # Construct the spacecraft matching the expected pretty-print
    state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03])
    t     = Time("2015-09-21T12:23:12", TAI(), ISOT())
    cs    = CoordinateSystem(earth, ICRFAxes())
    sc    = Spacecraft(state=state, time=t, mass=1000.0, name="SC-001", coord_sys=cs)

    expected = """
Spacecraft: SC-001
  AstroEpochs.Time
    value  = 2015-09-21T12:23:12.000
    scale  = TAI()
    format = ISOT()
  OrbitState:
    statetype: Cartesian
  CartesianState:
    x   =  7000.00000000
    y   =   300.00000000
    z   =     0.00000000
    vx  =     0.00000000
    vy  =     7.50000000
    vz  =     0.03000000
  CoordinateSystem:
    origin = Earth
    axes   = ICRFAxes
  Total Mass = 1000.0 kg
  CADModel: (no model)
"""

    # Capture pretty text/plain representation and compare (ignore trailing newline)
    shown = repr(MIME"text/plain"(), sc)
    @test chomp(shown) == chomp(expected)
end

@testset "AstroModels.Spacecraft posvel accessors" begin
    # Construct a Cartesian OrbitState spacecraft (fast path for posvel getters/setters)
    pv0 = [7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]
    state0 = CartesianState(pv0)
    t0     = Time("2015-09-21T12:23:12", TAI(), ISOT())
    cs     = CoordinateSystem(earth, ICRFAxes())
    sc     = Spacecraft(state=state0, time=t0, mass=1000.0, name="SC-001", coord_sys=cs)

    # Getter returns the current pos/vel vector
    @test to_posvel(sc) == pv0

    # Setter updates the underlying OrbitState to new values
    pv1 = [7050.0, 0.0, 1.0, 0.0, 7.6, 0.0]
    set_posvel!(sc, pv1)
    @test to_posvel(sc) == pv1

    # Setter enforces length==6
    @test_throws ArgumentError set_posvel!(sc, [1.0, 2.0, 3.0])
end

@testset "AstroModels.Spacecraft deepcopy" begin
    pv0    = [7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]
    state0 = CartesianState(pv0)
    t0     = Time("2015-09-21T12:23:12", TAI(), ISOT())
    cs0    = CoordinateSystem(earth, ICRFAxes())
    sc0    = Spacecraft(state=state0, time=t0, mass=1000.0, name="SC-001", coord_sys=cs0)

    sc1 = deepcopy(sc0)

    # Different objects; equal field values
    @test sc1 !== sc0
    @test sc1.name == sc0.name == "SC-001"
    @test sc1.mass == sc0.mass == 1000.0
    @test sc1.time == sc0.time
    @test sc1.coord_sys.origin === earth
    @test sc1.coord_sys.axes isa ICRFAxes
    @test to_posvel(sc1) == to_posvel(sc0) == pv0

    # Test time mutation independence
    sc0.time = Time("2024-10-03T12:00:00", UTC(), ISOT())
    @test sc1.time == t0

    # Test state mutation independence
    set_posvel!(sc0, [7050.0, 0.0, 1.0, 0.0, 7.6, 0.0])
    @test to_posvel(sc0) == [7050.0, 0.0, 1.0, 0.0, 7.6, 0.0]
    @test to_posvel(sc1) == pv0
    to_posvel(sc1) == to_posvel(sc0) == pv0

    # Test name mutation independence
    sc0.name = "SC-CHANGED"
    @test sc1.name == "SC-001"

    # Test coord_sys mutation independence
    sc0.coord_sys = CoordinateSystem(moon, ICRFAxes())
    @test sc1.coord_sys.origin == earth
    @test sc1.coord_sys.axes == ICRFAxes()

    # History independence
    @test length(sc1.history) == 0
    seg = HistorySegment([t0], [CartesianState([1.0,2.0,3.0,4.0,5.0,6.0])], sc0.coord_sys)
    push_segment!(sc0.history, seg)
    @test length(sc0.history) == 1
    @test length(sc1.history) == 0
end

@testset "AstroModels.Spacecraft to_posvel copy semantics" begin
    pv0 = [7000.0, 300.0, 0.0, 0.0, 7.5, 1.0]
    sc  = Spacecraft(state=CartesianState(pv0), time=Time("2015-09-21T12:23:12", TAI(), ISOT()))

    # 1) to_posvel returns a copy; mutating it does not affect sc
    v = to_posvel(sc)
    v[6] += 0.26566252414
    @test sc.state.statetype == Cartesian()
    @test sc.state.state[6] == pv0[6]
    @test to_posvel(sc)[6] == pv0[6]

    # Snapshot stays stable even if later calls are mutated
    snapshot = to_posvel(sc)
    tmp = to_posvel(sc); tmp[1] += 1.0
    @test snapshot[1] == pv0[1]

    # 2) set_posvel! applies updates only when called
    newpv = copy(v)
    set_posvel!(sc, newpv)
    @test to_posvel(sc) == newpv
end

@testset "Spacecraft numeric promotion with Dual mass" begin
    # Build inputs: Float64 state, concrete Time, Dual mass
    state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03])
    t     = Time("2015-09-21T12:23:12", TAI(), ISOT())
    m     = ForwardDiff.Dual(1000.0, 1.0)  # 1 active partial

    sc = Spacecraft(state=state, time=t, mass=m, coord_sys=CoordinateSystem(earth, ICRFAxes()))

    # Mass promoted to Dual
    @test sc.mass isa ForwardDiff.Dual

    # State elements promoted to Dual and statetype preserved
    @test eltype(sc.state.state) <: ForwardDiff.Dual
    @test sc.state.statetype == Cartesian()

    # Time value preserved (do not assert on internal numeric type)
    @test sc.time == t

    # History always stores Float64 regardless of spacecraft promotion
    seg = HistorySegment([sc.time], [CartesianState(zeros(6))], sc.coord_sys)
    push_segment!(sc.history, seg)
    @test length(sc.history) == 1
    @test length(sc.history.segments[1].times) == 1
    @test sc.history.segments[1].times[1] isa Time{Float64}
    @test sc.history.segments[1].states[1] isa CartesianState{Float64}
end

@testset "Spacecraft with Dual mass" begin
    state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03])
    t  = Time("2015-09-21T12:23:12", TAI(), ISOT())
    m  = ForwardDiff.Dual(1000.0, 1.0)  # 1 active partial
    sc = Spacecraft(state=state, time=t, mass=m)

    @test sc.state.state isa Vector{<:ForwardDiff.Dual}
    @test sc.mass isa ForwardDiff.Dual
    @test sc.time.jd1 isa ForwardDiff.Dual
    @test sc.time.jd2 isa ForwardDiff.Dual
    @test sc.history isa SpacecraftHistory
end

@testset "Spacecraft with Dual time" begin
    state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03])
    m = 1500.0
    jd1 = ForwardDiff.Dual(215450.0, 1.0)
    jd2 = ForwardDiff.Dual(0.5, 1.0)
    t = Time(jd1, jd2, TDB(), JD())
    sc = Spacecraft(state=state, time=t, mass=m)

    @test sc.state.state isa Vector{<:ForwardDiff.Dual}
    @test sc.mass isa ForwardDiff.Dual
    @test sc.time.jd1 isa ForwardDiff.Dual
    @test sc.time.jd2 isa ForwardDiff.Dual
    @test sc.history isa SpacecraftHistory
end

@testset "Spacecraft with Dual state" begin
    state = CartesianState(ForwardDiff.Dual.([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03], 6.0))
    t  = Time("2015-09-21T12:23:12", TAI(), ISOT())
    m = 1500.0
    sc = Spacecraft(state=state, time=t, mass=m)

    @test sc.state.state isa Vector{<:ForwardDiff.Dual}
    @test sc.mass isa ForwardDiff.Dual
    @test sc.time.jd1 isa ForwardDiff.Dual
    @test sc.time.jd2 isa ForwardDiff.Dual
    @test sc.history isa SpacecraftHistory
end

@testset "gaps" begin
  # Prepare a baseline spacecraft in Cartesian
  state_cart = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.546, 0.0])
  t_tt       = Time("2015-09-21T12:23:12", TT(), ISOT())
  cs_icrf    = CoordinateSystem(earth, ICRFAxes())
  sc_cart    = Spacecraft(state=state_cart, time=t_tt, mass=500.0, coord_sys=cs_icrf)

  # Test get pos_vel function (fast path for Cartesian)
  @test AstroModels.to_posvel(sc_cart) == to_vector(state_cart)

  # Test that to_posvel now works with non-Cartesian states (previously threw error)
  sc_kep_for_to = Spacecraft(state=KeplerianState(7000.0, 0.01, 0.1, 0.2, 0.3, 0.4),
                            time=t_tt, mass=1.0, coord_sys=cs_icrf)
  posvel = AstroModels.to_posvel(sc_kep_for_to)
  @test length(posvel) == 6
  @test posvel isa Vector

  # Test that set_posvel! now works with non-Cartesian states (previously threw error)
  sc_kep_for_set = Spacecraft(state=KeplerianState(7000.0, 0.01, 0.1, 0.2, 0.3, 0.4),
                              time=t_tt, mass=1.0, coord_sys=cs_icrf)
  AstroModels.set_posvel!(sc_kep_for_set, [7100.0, 100.0, 200.0, 0.5, 7.5, 0.1])
  @test sc_kep_for_set.state.statetype == Keplerian()  # State type preserved
  @test length(AstroModels.to_posvel(sc_kep_for_set)) == 6

end

# Construct a CoordinateSystem with an origin that has no μ field
struct FakeBody <: EpicycleBase.AbstractPoint
    name::String
end

@testset "get_state error when μ is missing in coord system origin" begin
    fake_coords = CoordinateSystem(FakeBody("FakePoint"), ICRFAxes())
    sc = Spacecraft(coord_sys=fake_coords)  # defaults to Cartesian state

    # Capture the exact error and assert on message substrings
    err = try
        get_state(sc, Keplerian())
        nothing
    catch e
        e
    end

    @test err isa ErrorException
    msg = sprint(showerror, err)

    @test occursin("get_state: μ is required to convert from Cartesian() to Keplerian()", msg)
    @test occursin("system does not have a celestial body (with a μ) as its origin", msg)
    @test occursin("Change to state types", msg)
end

@testset "Spacecraft promotion for AD" begin
    using ForwardDiff
    
    # Create a Float64 spacecraft with some history
    sc = Spacecraft(
        state=CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]),
        time=Time("2015-09-21T12:23:12", TAI(), ISOT()),
        mass=1000.0,
        name="test_sat"
    )
    
    # Add some history to verify it's preserved
    times = [Time("2015-09-21T12:23:12", TAI(), ISOT()), Time("2015-09-21T12:24:12", TAI(), ISOT())]
    states = [CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]), CartesianState([7001.0, 301.0, 1.0, 0.1, 7.4, 0.04])]
    seg = HistorySegment(times, states, sc.coord_sys)
    push_segment!(sc.history, seg)
    
    # Promote to ForwardDiff.Dual type
    DualType = ForwardDiff.Dual{Nothing,Float64,3}
    sc_dual = promote(sc, DualType)
    
    # Verify state was promoted
    @test eltype(sc_dual.state.state) === DualType
    @test sc_dual.state.statetype == sc.state.statetype
    
    # Verify time was promoted
    @test typeof(sc_dual.time.jd1) === DualType
    @test typeof(sc_dual.time.jd2) === DualType
    @test sc_dual.time.scale == sc.time.scale
    @test sc_dual.time.format == sc.time.format
    
    # Verify mass was promoted  
    @test typeof(sc_dual.mass) === DualType
    @test ForwardDiff.value(sc_dual.mass) ≈ sc.mass
    
    # Verify history remains Float64 (efficient storage)
    @test sc_dual.history isa SpacecraftHistory
    @test length(sc_dual.history) == length(sc.history)
    @test sc_dual.history.segments[1].times[1] isa Time{Float64}
    @test sc_dual.history.segments[1].states[1] isa CartesianState{Float64}
    
    # Verify non-numeric fields preserved
    @test sc_dual.name == sc.name
    @test sc_dual.coord_sys == sc.coord_sys
    
    # Test that original spacecraft is unchanged
    @test eltype(sc.state.state) === Float64
    @test typeof(sc.mass) === Float64
end

include("runtests_cadmodel.jl")
include("runtests_spacecraft_historysegment.jl")
include("runtests_spacecraft_history.jl")
include("runtests_morestatetests.jl")

nothing
