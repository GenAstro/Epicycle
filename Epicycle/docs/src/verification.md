# Verification

Epicycle's numerical models are verified against established software, published standards and
mathematical definitions. AstroEpochs matches Astropy
bit for bit, apart from its deliberate treatment of UTC before 1972, and AstroProp generally
agrees with GMAT and Orekit to within 2 cm in the individual Earth-orbit force-model comparisons,
with exceptions for drag and GMAT's Earth orientation. Full force-model comparisons about the
Moon, Mars, Phobos and Bennu agree to millimetres when time conventions match.

For each package, this page describes the models compared, the reference data and how they were
produced, the test cases, the measured agreement and every known difference from the reference.
Comparisons included in a package's automated test suite run in continuous integration on every
change, without the reference software installed.

| Package | Reference | Cases | Measured agreement | Test tolerance |
|:---|:---|:---|:---|:---|
| [AstroEpochs](@ref verification-astroepochs) | Astropy 8.0.1 (ERFA 2.0.1) | 1,790 epochs in six scales; 94 ERFA routine cases | Bit identical | 2e-11 s |
| [AstroProp](@ref verification-astroprop) | GMAT R2026A, Orekit 13.1 | 5 Earth orbits, 10 force models, 80 comparisons; full models about the Sun, Moon, Mars, Phobos and Bennu | Within 2 cm, except for drag and the Earth's field compared with GMAT, whose pole differs; millimetres about the Moon, Mars, Phobos and Bennu with the reference's time conventions | Not yet in the test suite |

Further packages are added to this page as their verification is written up.

## Approach

Each comparison uses one of three kinds of reference. A definition is a value fixed by a standard
or by the mathematics, such as the 32.184 s offset between TT and TAI. A reference implementation
is established software whose output is taken as correct, such as Astropy for time and GMAT for
orbit states. An internal consistency check, such as a round trip or a comparison of automatic and
finite-difference derivatives, needs no external value but bounds the error a calculation can
introduce.

When Epicycle ports a reference algorithm, agreement shows that the port is faithful; the model
has the reference's accuracy. When Epicycle implements its own model, agreement with an
independent reference provides evidence of correctness. Each section identifies which case
applies.

Tools can differ in convention as well as in model, including the axes used to evaluate a field
and the time scales used for integration and ephemeris evaluation. In these cases, the comparison
is run twice: first with Epicycle's conventions, at the same physical instants, then with the
reference tool's conventions. Matching conventions isolates the model comparison. If the models
agree, the original difference is due to convention alone. The section then explains which
tool's convention is more correct.

The test suite checks each comparison against a tolerance chosen for the quantity's precision.
Passing a tolerance establishes only an upper bound on error, so this page also reports the
measured agreement. Scripts kept with the tests generate the reference data once and save them
as data files. Rerunning those scripts against a newer reference release produces a new benchmark.

## [AstroEpochs: time scales and formats](@id verification-astroepochs)

AstroEpochs is tested against Astropy, and its ported ERFA routines against pyerfa. The tests
compare scale and format conversions across 1,790 epochs, including edge cases such as instants
inside leap seconds and conversions that round across a boundary. AstroEpochs agrees with Astropy
to the last bit at every epoch tested, in all six time scales and in every format. The one
deliberate difference is UTC before 1972, described at the end of this section.

### Models

AstroEpochs represents an epoch as a two-part Julian date in one of six scales: TAI, TT, UTC, TDB,
TCB and TCG. It converts between these scales and between Julian date, modified Julian date and
ISO 8601 formats.

The scale relations follow the IAU resolutions as implemented in ERFA, the open derivative of the
IAU SOFA library. TT is TAI plus 32.184 s; TCG and TCB follow from TT and TDB by their defining
rates; TDB − TT is the Fairhead and Bretagnon series evaluated at the geocentre; and UTC is TAI
minus the leap-second count from the IERS list.

The ERFA routines are ported to Julia line by line. Calendar conversion, normalisation of the
two-part date and rounding of formatted strings follow Astropy. The comparison therefore tests
the faithfulness of the port.

The two-part date represents an instant to about 2e-11 s, two units in the last place of the day
fraction.

### Reference data

The reference is Astropy 8.0.1 with pyerfa 2.0.1.5 (ERFA 2.0.1) and astropy-iers-data
0.2026.9.28. Astropy's automatic IERS download is disabled so that each run uses the bundled
leap-second table and is reproducible. Each case records an input epoch in one scale and format,
along with Astropy's two-part Julian date and ISO string for that instant in all six scales.

The data and the scripts that produced them are in
[AstroEpochs/test/astropy](https://github.com/GenAstro/Epicycle/tree/main/AstroEpochs/test/astropy):
`reference.csv` holds the 1,793 epochs, written by `make_reference.py`, and `make_erfa_parity.py`
produces the per-routine values.

### Test cases and results

The test cases cover conditions where time conversions are most likely to fail: leap seconds,
rounding boundaries and dates far from the present. Each epoch is converted into all six scales and formatted
as an ISO string in each, and every result is compared with Astropy's.

| Case kind | Epochs | What it exercises | Result |
|:---|:---|:---|:---|
| Random, 1972–2100 | 1,200 | Every scale and format as input, across the leap-second era | Bit identical |
| Either side of a leap second | 270 | UTC and TAI just before and after every leap second | Bit identical |
| Inside a leap second | 162 | UTC instants at 23:59:60.x | Bit identical |
| Rounding across a boundary | 15 | ISO strings that round across a second, minute, hour or day | Bit identical |
| Fixed epochs | 18 | J2000, the 1977 TCG/TCB epoch, MJD 0, and dates from 1600 to 2500 | Bit identical |
| Two-part MJD | 125 | Modified Julian dates whose first part is fractional | Bit identical |
| UTC before 1972 | 3 | Dates before the leap-second system | Differs by design, below |

Bit identical means every one of the 10,740 converted dates (1,790 epochs in six scales) equals
Astropy's in both parts of the two-part Julian date, bit for bit, and every one of the 10,740 ISO
strings matches character for character. The test tolerance is 2e-11 s.

The ported ERFA routines are also compared directly with pyerfa in 94 cases across 13 routines.
These cover the TAI, TT, TCG, TDB and TCB transforms in both directions, the TDB − TT series,
UTC-to-TAI conversions and their reverse across a leap second, and the calendar and time-of-day
routines. The inputs cover each routine's boundary cases: half days, where ERFA rounds away from
zero and Julia rounds to even; both orders of a two-part date; leap-second days; and times that
round across midnight. All 94 cases agree exactly.

Converting each epoch to every other scale and back gives 53,700 round trips. Each returns the
starting date to within 2.9e-11 s. The longest route, TCB to UTC, takes four conversions each way,
with one rounding per conversion. The resulting error is the expected accumulation of a few
units in the last place. The tolerance is 1e-10 s.

The leap-second table is read from the IERS list published by IANA. The test suite checks the
latest entry, 37 s from 1 January 2017, verifies that the list has not expired, and checks the
transforms on both sides of every change.

Epochs accept dual numbers, allowing derivatives through a scale conversion. Forward and reverse
mode agree and match the analytic rate of one day per 86,400 s.

### Differences from the reference

UTC before 1972 is the one deliberate difference. Before the present leap-second system, UTC
drifted against TAI at rates that changed several times. Astropy reproduces those offsets;
AstroEpochs takes TAI − UTC as zero and warns once per session. The difference at the three
reference epochs is 0.943 s on 1 January 1960, 3.855 s on 15 June 1965 and 9.892 s on 31 December
1971. The test suite records these as known failures, so they stay visible, and the test fails if
the behaviour changes without the test being updated. An epoch in any other scale is unaffected
unless it is converted to or from UTC.

ISO strings are formatted to the millisecond, matching Astropy's default. The underlying
two-part date retains full precision.

## [AstroProp: orbit propagation](@id verification-astroprop)

AstroProp is tested against two independent astrodynamics tools, GMAT and Orekit. When configured
identically, Epicycle agrees with them to millimetres in almost all cases, across the supported
force models and on orbits about the Earth, the Moon, Mars, Phobos and Bennu, and to 1.4 cm over a
month about the Sun. The results below include both comparisons as Epicycle runs and
comparisons with Epicycle configured as the reference tool runs. The key differences between the
tools found during testing are:

- Time scale. Epicycle integrates in Terrestrial Time (TT) about the Earth and in Barycentric
  Dynamical Time (TDB) about every other body. GMAT and Orekit integrate in atomic time, which runs
  at the rate of TT, about every body, and read the ephemeris at a two-term approximation of TDB
  that differs from the series Epicycle uses by 57 µs at the test epoch.
- Earth-fixed frame. GMAT's Earth-fixed frame leaves out the IERS celestial-pole offsets, which
  Epicycle applies, so GMAT's pole is 46.6 mas from the IERS pole. This moves orbits under the
  Earth's gravity field by up to 19 cm in a day of low Earth orbit and 4.1 m over three days on
  HEO. Its effect on drag is below a millimetre.
- Atmospheric density above each model's altitude limit. Epicycle returns zero density above
  1,000 km for NRLMSISE-00 and Harris-Priester, and above 3,000 km for JB2008. Orekit continues
  to evaluate NRLMSISE-00 and JB2008 above those limits, where its density is small but not zero.

With Epicycle configured as each reference tool runs, a configuration used only in testing that
the user interface does not offer, agreement is 4.9 mm with GMAT for EGM96 70×70 over a day of low
Earth orbit, between 0.9 and 2.6 mm about the Moon, Mars, Phobos and Bennu, and 1.4 cm over a month
about the Sun.

With individual gravity, third-body, SRP, tide and relativity models, Epicycle agrees with Orekit
to within 2 cm on every tested orbit and 1.5 mm in low Earth orbit. Agreement with GMAT is as
close except where Earth orientation affects the comparison; in the one case run in GMAT's Earth
axes, EGM96 70×70 for a day at 500 km, the gravity field agrees to 4.9 mm.

Drag agreement depends on the atmosphere model. In low Earth orbit, the exponential atmosphere
agrees with GMAT to 15 cm out of a 40 km effect, and NRLMSISE-00 agrees with Orekit to 5 mm.
JB2008 and Harris-Priester differ from Orekit by up to 1.1 m and 11 m. With all forces combined,
the Orekit differences are 0.7 to 9.3 mm in low Earth orbit, 0.1 mm at GEO, 2.1 mm at MEO and 3.6 m
on HEO.
Each combined difference also appears in the individual-model comparisons. The 3.6 m on HEO comes
from NRLMSISE-00 above Epicycle's 1,000 km density cutoff, where the spacecraft spends most of its
time and Orekit's density is not zero.

With full force models and matching time conventions, agreement with Orekit and GMAT is at the
millimetre level about the Moon, Mars, Phobos and Bennu, and within 1.4 cm over a month about the
Sun. Epicycle integrates about these bodies in TDB; Orekit and GMAT use TT. The results below
separate these convention differences from model differences.

### Models

The Earth-orbit comparisons test the following AstroProp and EpicycleEnterprise force models
individually:

- Harmonic gravity: EGM96 to degree and order 70, evaluated by Pines' formulation.
- Third bodies: the Sun and Moon as point masses, from the DE440 ephemeris.
- Solar radiation pressure: a cannonball model with the Earth's conical shadow, umbra and
  penumbra. The integrator's steps end on the shadow boundaries, where the force switches on and
  off.
- Drag: a cannonball model with four atmospheres. The exponential atmosphere is Vallado's
  table; NRLMSISE-00, JB2008 and Harris-Priester are SatelliteToolbox's implementations, reading
  the published space-weather indices.
- Solid tides: IERS Conventions (2010) step 1, with the Love numbers k₂ = 0.30190 and
  k₃ = 0.093. The tides are carried by the gravity field, whose reference radius and tide system
  they use. The comparisons use EGM96 4×4.
- Relativity: the post-Newtonian correction of IERS Conventions (2010) eq. 10.12: the
  Schwarzschild, Lense-Thirring and de Sitter terms.

The gravity, tide and relativity models are Epicycle's own implementations, so GMAT and Orekit are
independent checks of them. The density models are SatelliteToolbox's, so for them the comparison
tests both the models and how AstroProp applies them. The integrator is Vern9 from
OrdinaryDiffEq.

### Reference data

The references are GMAT R2026A and Orekit 13.1. GMAT uses the IAU 1976/FK5 chain for Earth
orientation, while Orekit uses the IAU 2006 theory. Epicycle runs every case twice, once with
each tool's frame theory.

GMAT's Earth-fixed axes differ from FK5 as defined by the IERS conventions. GMAT reads polar
motion and UT1 from the daily Earth orientation parameters but leaves out the observed
celestial-pole offsets δΔψ and δΔε, which correct the 1976 precession and 1980 nutation theory to
the measured pole. On 2020-10-20, GMAT's Earth-fixed rotation is within 0.35 mas of the
IAU 1976/FK5 chain without these offsets. It differs by 46.6 mas from both of Epicycle's frame
theories, which apply the offsets and agree with each other to 0.25 mas. GMAT's pole is therefore
about 1.4 m from the IERS pole at the Earth's surface. This affects low-orbit propagation under
a full gravity field:

| EGM96 70×70, one day at 500 km | Difference from GMAT |
|:---|:---|
| Epicycle in GMAT's Earth axes | 4.9 mm |
| Epicycle in its own FK5 axes | 11.5 cm |

The case above uses GMAT's axes. All other GMAT comparisons use Epicycle's FK5 axes, so their
gravity-field and tide results include this orientation difference. Its effect on drag is below a
millimetre.

Constants outside the model under test are matched between the tools. GMAT is given Epicycle's
gravitational parameters of the Earth, Sun and Moon. Orekit is given Epicycle's for the Earth, but
its third bodies take theirs from the DE440 header with no way to give them others, so Epicycle's
Orekit-matched runs use Orekit's Sun and Moon values; the Moon's differs from Epicycle's by
2.4×10⁻⁸, which unmatched was 1.9 mm at GEO. Both tools are given Epicycle's Earth and Sun radii
for the shadow, and its solar flux. All tools use DE440. Orekit's NRLMSISE-00 and JB2008 use the same
indices as Epicycle.

Every run in each tool is converged: tolerances are tight enough and steps short enough that
reducing the step does not change the result. Differences therefore reflect the force models.
GMAT does not end its steps at shadow boundaries and needs a shorter maximum step than the other
tools. With a 2 s step, its LEO runs with SRP were 2 cm from their converged result; with a 0.5 s
step, they were within 0.1 mm of Epicycle.

Epicycle measures the integration floor for each run by repeating it with a tenfold tighter
tolerance and half the step. The floor is below 0.3 mm in every run without drag and below 1 cm
in every run.

The cases, the scripts that ran GMAT and Orekit, and their ephemerides are kept in Gen Astro's
internal test suite and are available on request.

### Test cases

The orbits are those of the GMAT verification paper[^gmatvv], at the epoch
2020-10-20T12:00:00 UTC, so the Earth orientation and space-weather data are final values. Each
orbit is propagated with the point-mass Earth alone, with one force model added at a time, and
with all forces together, and the state is compared at every output step.

| Orbit | Altitude (km) | Inclination | Duration | Output step |
|:---|:---|:---|:---|:---|
| LEO, Sun-synchronous | 400, circular | 97.0° | 1 day | 60 s |
| LEO, ISS | 358 × 380 | 51.7° | 1 day | 60 s |
| GEO | 35,786, circular | 0.0° | 7 days | 600 s |
| MEO | 19,758 × 20,604 | 56.3° | 2 days | 120 s |
| HEO, Molniya | 500 × 39,850 | 63.4° | 3 days | 300 s |

The spacecraft has a mass of 1000 kg, a drag coefficient of 2.2 and a reflectivity coefficient of
1.8. Drag and SRP each use an area of 10 m². The comparisons use the drag models available in each
reference tool: the exponential atmosphere in GMAT, and NRLMSISE-00, JB2008 and Harris-Priester
in Orekit.

[^gmatvv]: S. P. Hughes, R. H. Qureshi, D. S. Cooley, J. J. K. Parker and T. G. Grubb, "Verification and Validation of the General Mission Analysis Tool (GMAT)," AIAA/AAS Astrodynamics Specialist Conference, AIAA 2014-4151, 2014.

### Results

Each cell is the largest distance between Epicycle's position and the reference tool's over the
whole propagation. A dash marks a model with no effect on that orbit, such as drag at GEO.

The GMAT harmonic-gravity and solid-tide results include the Earth-orientation difference
described under Reference data. GMAT's pole differs from the IERS pole by 46.6 mas, or
2.3×10⁻⁷ rad. The position difference is that fraction of the field's effect to within a factor of
two: 17 to 19 cm out of 510 to 1,100 km in low Earth orbit, and 4.1 m out of 17,000 km on HEO. In GMAT's own axes,
the same field agrees to 4.9 mm.

#### Against GMAT

| Orbit | Point mass | Harmonic gravity | Solid tides | Sun and Moon | SRP | Drag, exponential | Relativity |
|:---|---:|---:|---:|---:|---:|---:|---:|
| LEO, Sun-synchronous | 24 µm | 19 cm | 19 cm | 85 µm | 29 µm | 13 cm | 46 µm |
| LEO, ISS | 110 µm | 17 cm | 16 cm | 47 µm | 63 µm | 15 cm | 66 µm |
| GEO | 120 µm | 1.5 cm | 1.5 cm | 1.2 cm | 280 µm | 180 µm | 26 µm |
| MEO | 39 µm | 2.0 cm | 2.0 cm | 2.9 mm | 81 µm | 48 µm | 67 µm |
| HEO | 1.2 mm | 4.1 m | 4.1 m | 1.9 cm | 400 µm | 2.0 cm | 1.7 mm |

#### Against Orekit

| Orbit | Point mass | Harmonic gravity | Solid tides | Sun and Moon | SRP | Drag, NRLMSISE-00 | Drag, JB2008 | Drag, Harris-Priester | Relativity |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| LEO, Sun-synchronous | 4.6 µm | 730 µm | 770 µm | 6.5 µm | 6.7 µm | 350 µm | 61 cm | 11 m | 13 µm |
| LEO, ISS | 3.4 µm | 1.5 mm | 1.5 mm | 26 µm | 5.9 µm | 5.1 mm | 1.0 m | 11 m | 16 µm |
| GEO | 22 µm | 100 µm | 94 µm | 16 µm | 14 µm | – | – | – | 58 µm |
| MEO | 5.6 µm | 200 µm | 220 µm | 1.4 µm | 5.4 µm | 2.0 mm | 2.8 cm | – | 3.2 µm |
| HEO | 110 µm | 1.6 cm | 1.6 cm | 25 µm | 48 µm | 4.4 m | 1.1 m | 10 m | 300 µm |

#### Size of each model's effect

| Orbit | Harmonic gravity | Solid tides | Sun and Moon | SRP | Drag, exponential | Drag, NRLMSISE-00 | Drag, JB2008 | Drag, Harris-Priester | Relativity |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| LEO, Sun-synchronous | 510 km | 34 m | 71 m | 18 m | 23 km | 7.4 km | 4.6 km | 27 km | 2.6 m |
| LEO, ISS | 1,100 km | 38 m | 29 m | 22 m | 40 km | 14 km | 9.5 km | 44 km | 2.7 m |
| GEO | 100 km | 58 cm | 100 km | 2.7 km | – | – | – | – | 1.2 m |
| MEO | 59 km | 99 cm | 2.3 km | 310 m | – | – | – | – | 63 cm |
| HEO | 17,000 km | 120 m | 32 km | 1.4 km | 2.3 km | 700 m | 450 m | 5.7 km | 83 m |

#### Every force at once

| Orbit | GMAT | GMAT, without tides | Orekit | Orekit, without tides | Effect of all forces |
|:---|---:|---:|---:|---:|---:|
| LEO, Sun-synchronous | 28 cm | 28 cm | 710 µm | 700 µm | 520 km |
| LEO, ISS | 32 cm | 32 cm | 9.3 mm | 9.3 mm | 1,100 km |
| GEO | 2.3 cm | 2.3 cm | 87 µm | 100 µm | 18 km |
| MEO | 1.8 cm | 1.8 cm | 2.1 mm | 2.1 mm | 62 km |
| HEO | 4.1 m | 4.1 m | 3.6 m | 3.6 m | 17,000 km |

Combined force models agree as closely as the individual models. Each combined difference is
also present in an individual-model comparison.

For Orekit, the differences come from the gravity field in low Earth orbit, NRLMSISE-00 on the
ISS orbit, and NRLMSISE-00 above 1,000 km at MEO and on HEO, where
Epicycle returns zero density and Orekit does not. For GMAT, the differences come from its pole
on all orbits, together with the exponential atmosphere in low Earth orbit and the third bodies
at GEO. Omitting tides changes no result by more than 0.3 mm.

Velocity differences follow the position differences: the largest ratio of the two is
1.2e-3 per second, so a 1 cm position difference comes with at most 12 µm/s.

### Full force models about other bodies

These cases compare the full force model supported by each body, with all forces applied together.
The Sun, Moon and Mars cases use Orekit as the reference; the Phobos and Bennu cases use GMAT.
Epicycle does not ship body definitions for Phobos or Bennu.

Each orbit about a body that casts a shadow is circular, with its plane containing the Sun
direction, so it passes through the body's shadow on every revolution. The spacecraft has an
area-to-mass ratio of 1 m² per kg and a reflectivity coefficient of 1.8. The area-to-mass ratio is
a hundred times that of a typical spacecraft, so SRP and the shadow boundaries weigh heavily in the
comparison. The Bennu case uses 0.02 m² per kg because the higher ratio would push the spacecraft
away from Bennu's weak gravity.

| Case | Orbit | Duration | Forces | Reference |
|:---|:---|:---|:---|:---|
| Sun | Near 1 AU, 0.02 AU outside the Earth's orbit | 30 days | The Sun; the eight planets as third bodies; SRP; relativity | Orekit |
| Moon | 100 km, circular | 1 day | GL0660B 50×50 with solid tides (k₂ = 0.02405, k₃ = 0.0089, raised by the Earth and Sun); the Earth and Sun; SRP in the Moon's shadow; relativity | Orekit |
| Mars | 400 km, circular | 1 day | JGM85F01 20×20 with solid tides (k₂ = 0.169, raised by the Sun); the Sun, Earth and Jupiter; SRP in Mars's shadow; relativity | Orekit |
| Phobos | 20 km | 6 hours | Phobos, Mars and the Sun as point masses; SRP in Phobos's shadow | GMAT |
| Bennu | 1 km, circular | 1 day | Bennu, the Sun, Earth and Jupiter as point masses; SRP in Bennu's shadow | GMAT |

Every orbit with a shadow spends between 12 % and 40 % of its time in umbra. The SRP effect ranges
from 0.56 km about Bennu to 2.8 km about the Moon. About the Sun, it is 28,000 km.

Phobos's sphere of influence barely extends beyond its surface, so no orbit about it is bound.
The case covers a quarter-day arc, about one revolution. Its body definition uses a gravitational
parameter of 7.087e-4 km³/s², a radius of 11.08 km and its NAIF ID. Its ephemeris comes from NAIF's
Mars satellite kernel, built on DE440.

Bennu is defined in the same way, with a gravitational parameter of 4.892e-9 km³/s² and a radius
of 0.2825 km. Its ephemeris comes from the OSIRIS-REx orbit solution (Farnocchia et al. 2021),
built on DE424. Each tool loads the same kernels.

In the Orekit cases, the relativity model includes only the Schwarzschild term: Orekit's
Lense-Thirring model is for the Earth, and its de Sitter term is set up differently from
Epicycle's. The GMAT
cases, Phobos and Bennu, have no relativity.

The Orekit cases use Orekit's gravitational parameters, taken from the DE440 header. The GMAT
cases use Epicycle's gravitational parameters, supplied to GMAT.

For bodies other than the Earth, the tools use different time conventions. Epicycle integrates
in TDB, the time scale used by the barycentric dynamics and ephemeris. Orekit and GMAT integrate
in SI seconds counted in TT, regardless of the central body. They evaluate the ephemeris using
their own two-term approximations of TDB, which differ from Epicycle's IAU series by 57 µs at
this epoch.

TDB and TT have the same mean rate but differ by a periodic term of 1.66 ms over a year. The
clocks therefore drift apart by up to 0.44 ms over 30 days. Each case is compared first using
Epicycle's conventions, at the same physical instants, then using the reference tool's
conventions: integrating in TT and evaluating the ephemeris and body axes at the tool's TDB.

| Case | Difference with Epicycle's time conventions, at the same instants | Difference with the reference's time conventions | Epicycle's integration floor |
|:---|:---|:---|:---|
| Sun | 13.2 m | 1.4 cm | 5.3 mm |
| Moon | 4.9 cm | 2.5 mm | 1.9 mm |
| Mars | 10.4 cm | 1.1 mm | 0.34 mm |
| Phobos | 4.1 mm | 2.6 mm | 0.03 mm |
| Bennu | 0.9 mm | 0.9 mm | 0.7 µm |

With the reference's time conventions, Epicycle agrees to millimetres about the Moon, Mars,
Phobos and Bennu through every shadow crossing. About the Sun, it agrees to 1.4 cm over a month,
10⁻¹³ of the orbit's radius. The difference between the two results columns is therefore due to
the time conventions; about the Sun, the Moon and Mars it is nearly all of the first column, and
about Phobos and Bennu, whose arcs are short and slow, little of it. The Sun case separates the two
time-convention effects:

| Epicycle, Sun case, 30 days | Difference from Orekit |
|:---|:---|
| Epicycle's conventions: integrated in TDB | 13.2 m |
| Integrated in TT, as Orekit does | 5.8 cm |
| Also reading the ephemeris at Orekit's approximation of TDB | 1.4 cm |

The 13.2 m difference is the spacecraft's speed of 29.6 km/s multiplied by the 0.44 ms clock
drift. Integrating in TT uses a different clock from the heliocentric dynamics. Epicycle's
convention is correct; Orekit's approximation produces a difference of about 10⁻¹⁰ of the
orbit's radius over a month.

### Differences from the reference

Epicycle returns zero density above each model's altitude limit: 1,000 km for
NRLMSISE-00 and Harris-Priester, and 3,000 km for JB2008. Orekit evaluates NRLMSISE-00 and JB2008
above those altitudes. This accounts for the full MEO drag differences of 2.0 mm and 2.8 cm, where
Epicycle's drag is zero. The 4.4 m NRLMSISE-00 difference on the HEO orbit grows through each
apogee, where Epicycle's drag is also zero, as well as at perigee.

JB2008 differs from Orekit's by 0.6 to 1.1 m and Harris-Priester by 10 to 11 m, with the same
space-weather indices. In low Earth orbit that is about one part in 10,000 of JB2008's effect and
one part in 4,000 of Harris-Priester's; on HEO, where the effects are smaller, about one part in
400 and 600.
Neither difference has yet been attributed to a specific model term.

GMAT omits the IERS celestial-pole offsets, while both of Epicycle's frame theories apply them,
as described under Reference data. The resulting 46.6 mas pole difference at the test epoch
accounts for the harmonic-gravity and tide differences, from 1.5 cm at GEO to 4.1 m on HEO.
Evaluating the field in GMAT's axes reduces the difference to 4.9 mm.

With the Sun and Moon as third bodies, GMAT differs from Epicycle by 1.2 cm at GEO and 1.9 cm on
the HEO orbit, where Orekit agrees with Epicycle to 16 µm and 25 µm. With the exponential
atmosphere, GMAT differs by 13 to 15 cm in low Earth orbit, about four parts in a million of the
effect. Neither the third-body nor the exponential-atmosphere difference has yet been attributed.

Two comparison-setup choices are separate from differences in Epicycle's models. The Orekit
comparison omits the de Sitter relativity term from both tools, because Orekit 13.1 sets it up
differently from Epicycle. The GMAT comparison includes it.

GMAT's NRLMSISE-00 and Jacchia-Roberts use GMAT's own space-weather file. These runs do not yet
match that file to Epicycle's indices, so the GMAT drag comparison uses only the exponential
atmosphere.
