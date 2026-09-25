# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# The astronomy core, as a user meets it: epochs, states, bodies and frames.
#
#   include("Epicycle/demo/astronomy_core_demo.jl")
#   AstronomyCoreDemo.astronomy_core_demo()                      # timed pauses
#   AstronomyCoreDemo.astronomy_core_demo(wait_for_enter = true)  # presenter mode
#   AstronomyCoreDemo.astronomy_core_demo(pace = 0.0)             # no pauses, for a quick check
#
# Four parts, one per package, each following that package's Quick Start on the documentation site:
# AstroEpochs, AstroStates, AstroUniverse, AstroFrames. Every line is typed as a user would type
# it and then evaluated, so what the screen shows is a session rather than a description of one.
# The figures between the transcripts are computed from what the lines just built.

module AstronomyCoreDemo

using LinearAlgebra
using Printf
using InteractiveUtils          # subtypes, as the AstroStates Quick Start uses it

using EpicycleBase
using AstroEpochs
using AstroStates
using AstroUniverse
using AstroFrames
using AstroModels
using AstroCallbacks
using EpicycleIO

export astronomy_core_demo

# The narration helpers, the caption strip, the @step macro and the REPL transcript.
include(joinpath(@__DIR__, "narration.jl"))

# ═══════════════════════════════════════════════════════════════════════════════
# Part 1 — epochs
# ═══════════════════════════════════════════════════════════════════════════════

function _act_epochs()
    _act(1, "Epochs", "AstroEpochs: one instant, every scale and format")

    _say("  An epoch is stored as two Float64 Julian dates, which holds microseconds across")
    _say("  decades, and it carries its own time scale.")
    _pause(1.5)
    _clear_transcript()

    t1 = _repl("t1 = Time(2451545.0, TT(), JD())")
    _repl("t2 = Time(51544.5, UTC(), MJD())")
    _repl("""t3 = Time("2000-01-01T12:00:00.000", TAI(), ISOT())""")
    _pause(1.5)

    _say("  The same instant, read back in whichever format the next tool wants.")
    _repl("t1.jd")
    _repl("t1.mjd")
    _repl("t1.isot")
    _pause(1.5)

    _say("  A scale is a property, so a conversion is an access rather than a function call.")
    _repl("t_utc = t1.utc")
    _repl("t_tdb = t1.tdb")
    _pause(2.0)

    # The offsets, computed from the epochs above rather than quoted.
    @step "Reading the offsets between the scales" begin
        tai_utc_2000 = (t1.tai.jd - t1.utc.jd) * 86400
        now_utc      = Time("2026-09-17T00:00:00.000", UTC(), ISOT())
        tai_utc_now  = (now_utc.tai.jd - now_utc.utc.jd) * 86400
        tt_tai       = (t1.tt.jd - t1.tai.jd) * 86400
        tdb_tt_2000  = (t1.tdb.jd - t1.tt.jd) * 86400
        tdb_tt_july  = let t = Time("2000-07-01T12:00:00.000", TT(), ISOT())
            (t.tdb.jd - t.tt.jd) * 86400
        end
    end
    println()
    _fact("TAI − UTC at J2000", @sprintf("%.0f", tai_utc_2000); unit = "s")
    _fact("TAI − UTC in 2026", @sprintf("%.0f", tai_utc_now); unit = "s")
    _fact("TT − TAI, fixed by definition", @sprintf("%.3f", tt_tai); unit = "s")
    _fact("TDB − TT in January", @sprintf("%+.5f", tdb_tt_2000); unit = "s")
    _fact("TDB − TT in July", @sprintf("%+.5f", tdb_tt_july); unit = "s")
    println()
    _say("  UTC needs a table because of the leap seconds. TDB − TT is the relativistic term,",
         color = _DIM)
    _say("  which passes through zero twice a year: the two dates above sit either side of it.",
         color = _DIM)
    _pause(3.0)
    return nothing
end

# ═══════════════════════════════════════════════════════════════════════════════
# Part 2 — states
# ═══════════════════════════════════════════════════════════════════════════════

function _act_states()
    _act(2, "States", "AstroStates: one orbit, whichever representation suits the task")

    _say("  A state is a struct with named elements, not six numbers in an agreed order. A")
    _say("  conversion is a constructor: the type you want, from the state you have.")
    _pause(1.5)
    _clear_transcript()

    cart = _repl("cart = CartesianState([7000.0, 0.0, 100.0, 0.0, 7.5, 2.5])")
    mu = _repl("mu = 398600.4418")
    kep = _repl("kep = KeplerianState(cart, mu)")
    _pause(1.5)

    _say("  Elements are read by name, and the whole state comes out as a vector when a solver")
    _say("  wants one.")
    _repl("kep.sma")
    _repl("rad2deg(kep.inc)")
    _repl("to_vector(kep)")
    _pause(1.5)

    _say("  Convert back, and the round trip has to return what it started with.")
    cart2 = _repl("cart2 = CartesianState(kep, mu)")
    _pause(1.0)

    @step "Checking the round trip and the orbit it describes" begin
        round_trip = maximum(abs, to_vector(cart2) .- to_vector(cart))
        period_min = 2π * sqrt(kep.sma^3 / mu) / 60
        altitude   = kep.sma * (1 - kep.ecc) - earth.equatorial_radius
    end
    println()
    _fact("largest round-trip difference", @sprintf("%.2e", round_trip); unit = "km, km/s")
    _fact("period", @sprintf("%.2f", period_min); unit = "min")
    _fact("perigee altitude", @sprintf("%.1f", altitude); unit = "km")
    _pause(2.0)

    _say("  Twelve representations ship, including the ones an interplanetary arrival needs.")
    _repl("subtypes(AbstractOrbitState)")
    println()
    _say("  Each is a type, so a signature says which representation a function takes.",
         color = _DIM)
    _pause(3.0)
    return nothing
end

# ═══════════════════════════════════════════════════════════════════════════════
# Part 3 — bodies
# ═══════════════════════════════════════════════════════════════════════════════

function _act_universe()
    _act(3, "Bodies", "AstroUniverse: the solar system, and room for a body of your own")

    _say("  Eleven bodies ship with their constants. A body is an object you pass around, so a")
    _say("  force model or a frame takes the body rather than a number copied out of a table.")
    _pause(1.5)
    _clear_transcript()

    _repl("earth.mu")
    _repl("venus.naifid")
    _repl("moon.equatorial_radius")
    _pause(1.5)

    _say("  Your own body is the same struct, so everything downstream accepts it.")
    phobos = _repl("""phobos = CelestialBody(name = "Phobos", naifid = 401, mu = 7.0875e-4,
                      equatorial_radius = 11.1)""")
    _pause(1.5)

    _say("  Positions come from JPL ephemerides through SPICE, read at an epoch you name.")
    epoch = _repl("epoch = Time(\"2026-09-17T00:00:00.000\", UTC(), ISOT())")
    _repl("r_moon = translate(earth, moon, epoch.tdb.jd)")
    _pause(1.0)

    @step "Measuring the solar system at that epoch" begin
        d_moon  = norm(translate(earth, moon, epoch.tdb.jd))
        d_sun   = norm(translate(earth, sun, epoch.tdb.jd))
        d_mars  = norm(translate(earth, mars, epoch.tdb.jd))
        au      = 1.495978707e8
        g_phobos = phobos.mu / phobos.equatorial_radius^2 * 1e3      # m/s² at the surface
    end
    println()
    _fact("Earth to Moon", @sprintf("%.0f", d_moon); unit = "km")
    _fact("Earth to Sun", @sprintf("%.6f", d_sun / au); unit = "AU")
    _fact("Earth to Mars", @sprintf("%.3f", d_mars / au); unit = "AU")
    _fact("surface gravity on Phobos", @sprintf("%.4f", g_phobos); unit = "m/s²")
    println()
    _say("  Six millimetres per second squared at the surface: a visiting spacecraft hovers")
    _say("  rather than lands.", color = _DIM)
    _pause(3.0)
    return nothing
end

# ═══════════════════════════════════════════════════════════════════════════════
# Part 4 — frames
# ═══════════════════════════════════════════════════════════════════════════════

function _act_frames()
    _act(4, "Frames", "AstroFrames: one state, read against different origins and axes")

    _say("  A coordinate system is an origin and a set of axes. A coordinate is a state, a")
    _say("  coordinate system and an epoch, and a conversion is again a constructor.")
    _pause(1.5)
    _clear_transcript()

    cs_j2000 = _repl("cs_j2000 = CoordinateSystem(earth, MJ2000Eq())")
    cs_fixed = _repl("cs_fixed = CoordinateSystem(earth, ITRF())")
    epoch = _repl("epoch = Time(2458849.5, 0.0, TDB(), JD())")
    c = _repl("c = Coordinate([7000.0, 0.0, 0.0, 0.0, 7.546, 0.0], cs_j2000, epoch)")
    _pause(1.5)

    _say("  Earth-fixed axes turn with the planet, so the same instant reads differently there.")
    c_fixed = _repl("c_fixed = Coordinate(c, cs_fixed)")
    _pause(1.0)

    @step "Reading the Earth-fixed state as a place on the ground" begin
        x_fixed = to_vector(CartesianState(c_fixed))
        r_fixed = x_fixed[1:3]
        lat     = rad2deg(asin(r_fixed[3] / norm(r_fixed)))
        lon     = rad2deg(atan(r_fixed[2], r_fixed[1]))
        kept    = norm(r_fixed) - norm(to_vector(CartesianState(c))[1:3])
        speed_fixed   = norm(x_fixed[4:6])
        speed_inertial = norm(to_vector(CartesianState(c))[4:6])
    end
    println()
    _fact("sub-satellite latitude", @sprintf("%+.3f", lat); unit = "deg")
    _fact("sub-satellite longitude", @sprintf("%+.3f", lon); unit = "deg")
    _fact("change in |r| across the rotation", @sprintf("%.2e", kept); unit = "km")
    _fact("speed, inertial", @sprintf("%.4f", speed_inertial); unit = "km/s")
    _fact("speed, Earth-fixed", @sprintf("%.4f", speed_fixed); unit = "km/s")
    println()
    _say("  The rotation leaves |r| unchanged and subtracts the surface speed of the turning",
         color = _DIM)
    _say("  planet from the velocity.", color = _DIM)
    _pause(2.5)

    _say("  A spacecraft is a subject too, and a quantity takes the frame it is measured in.")
    sc = _repl("sc = Spacecraft(state = CartesianState([7000.0, 1000.0, 2000.0, -1.0, 7.0, 2.0]),
                 time = epoch)")
    _repl("rad2deg(raan(sc, cs_j2000))")
    _repl("rad2deg(raan(sc, CoordinateSystem(earth, MJ2000Ec())))")
    println()

    @step "The same orbit, measured against two planes" begin
        raan_eq = rad2deg(raan(sc, cs_j2000))
        raan_ec = rad2deg(raan(sc, CoordinateSystem(earth, MJ2000Ec())))
        inc_eq  = rad2deg(inclination(sc, cs_j2000))
        inc_ec  = rad2deg(inclination(sc, CoordinateSystem(earth, MJ2000Ec())))
    end
    println()
    _fact("RAAN, mean equator of J2000", @sprintf("%.3f", raan_eq); unit = "deg")
    _fact("RAAN, mean ecliptic of J2000", @sprintf("%.3f", raan_ec); unit = "deg")
    _fact("inclination, equator", @sprintf("%.3f", inc_eq); unit = "deg")
    _fact("inclination, ecliptic", @sprintf("%.3f", inc_ec); unit = "deg")
    println()
    _say("  Both values are correct, measured against different reference planes, which is why",
         color = _DIM)
    _say("  every quantity takes the coordinate system as an argument.", color = _DIM)
    _pause(3.0)
    return nothing
end

# ═══════════════════════════════════════════════════════════════════════════════
# The tour
# ═══════════════════════════════════════════════════════════════════════════════

"""
    astronomy_core_demo(; pace = 1.0, wait_for_enter = false, open_browser = true)

Run a narrated tour of Epicycle's astronomy core: epochs, orbital states, celestial bodies and
coordinate frames, following each package's Quick Start.

Every line is typed as a user would type it and evaluated live, so the terminal shows a session.
The figures between the transcripts are computed from the lines that ran.

# Arguments
- `pace`: scales every pause and the typing speed. `0` removes them, `2` doubles them.
- `wait_for_enter`: pause between steps until Enter is pressed, for presenting.
- `open_browser`: open the dashboard, whose caption strip carries the same script and figures.

# Returns
`nothing`.

# Notes
The first run downloads SPICE kernels and Earth orientation data if they are not already cached.
"""
function astronomy_core_demo(; pace::Real = 1.0, wait_for_enter::Bool = false,
                               open_browser::Bool = true)
    pace >= 0 || throw(ArgumentError("pace must be zero or positive; got $pace"))
    _PACE[]  = Float64(pace)
    _ENTER[] = wait_for_enter
    _TABS[]  = false                 # this tour draws nothing, so one tab carries all of it

    t_start = time()
    EpicycleIO.clear_all!()
    open_browser && open_dashboard()

    _title_card("The astronomy core",
                "   Epochs · orbital states · celestial bodies · coordinate frames")
    _pause(2.0)

    # Library notices are right for a script and noise in a narrated tour. Errors still show.
    Base.CoreLogging.with_logger(Base.CoreLogging.SimpleLogger(stderr,
                                                               Base.CoreLogging.Error)) do
        _act_epochs()
        _act_states()
        _act_universe()
        _act_frames()
    end

    c = _CAPTION
    c.title, c.subtitle, c.code, c.status = "That is the core", "", "", ""
    empty!(c.lines); empty!(c.facts); c.fresh = false

    println("\n")
    printstyled("  ", "━"^78, "\n"; color = _ACCENT)
    printstyled("   THAT IS THE CORE\n"; color = _ACCENT, bold = true)
    printstyled("  ", "━"^78, "\n\n"; color = _ACCENT)
    _fact("packages", "4")
    _fact("wall time", @sprintf("%.0f", time() - t_start); unit = "s")
    println()
    _say("  Four packages, and every line above is one you can paste into your own session.")
    _say("  The documentation for each is the page these transcripts follow.", color = _DIM)
    println()
    return nothing
end

end # module
