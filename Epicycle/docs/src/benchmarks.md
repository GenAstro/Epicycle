# Benchmarks

Epicycle's propagator is optimized for performance. This page compares it with GMAT R2026A and
Orekit 13.1 on one day of low Earth orbit under a full force model, in two ways. The first runs each
tool with matched integration methods and tolerances and measures each run's error from that
tool's own high-accuracy reference. The second runs all three with the same fixed step, set up as a
user of each tool would set it.

The results describe this case on one machine. Each tool makes its own design choices, in its
defaults, features and conventions, and timings in other cases will differ. Comparisons of
individual force models are on the [Verification](@ref verification-astroprop) page.

## Benchmark case

The benchmark includes gravity, third-body attraction, drag, solar radiation pressure (SRP), solid
tides and relativity. Epicycle uses a separate configuration for each comparison, matched to each
tool's force model configuration.

| Parameter | Value |
|:---|:---|
| Epoch | 2020-10-20T12:00:00 UTC |
| Initial state | [6878.137, 0, 0] km, [0, 4.71754, 5.99820] km/s, ICRF |
| Duration | 86,400 s |
| Spacecraft | 1000 kg; drag Cd 2.2, area 10 m²; SRP Cr 1.8, area 10 m² |
| Gravity | EGM96 70×70, with IERS 2010 solid tides (step 1) |
| Third bodies | Sun and Moon, DE440 |
| Drag | Exponential (GMAT comparison), NRLMSISE-00 (Orekit comparison) |
| SRP | Cannonball, 1367 W/m² at 1 AU, Earth's dual-cone shadow |
| Relativity | Schwarzschild, Lense-Thirring, and de Sitter in the GMAT comparison only |

The GMAT comparison uses an exponential atmosphere and GMAT's value for Earth's angular momentum.
The Orekit comparison uses NRLMSISE-00, with both tools reading the same space weather data, and
omits the de Sitter term in the relativistic correction.

## Comparison method

Each tool makes a reference run with a tight tolerance and a 2 s maximum step, so that its
integration error is small next to the errors measured here. Epicycle's runs at 1e-12 and 1e-13
without the step limit land within 0.6 mm of its reference. A timed run's error is its distance
from the same tool's reference after one day, and the tools are compared with each other through
their reference runs, so that difference reflects the force models rather than the integrators.

The step limit matters because the spacecraft enters and leaves Earth's shadow on every orbit, and
SRP switches on and off within a few seconds at each boundary. A step across a boundary can
introduce error that the integrator's error estimate misses, so tightening the tolerance does not
always reduce the propagation error. Epicycle ends its steps at the shadow boundaries.

The variable-step runs use corresponding integration methods and error-control settings where
available. Orekit's absolute tolerance in metres is 1000 times Epicycle's tolerance in kilometres,
with the same error norm. GMAT uses RSS error control relative to the state. Every variable-step
run starts with a 60 s step and has no maximum step; tool defaults are not used.

Timings were measured on a Windows 11 desktop with an Intel Core Ultra 7 265, the High Performance
power plan, and each process pinned to the eight performance cores. The tools ran in turn over
three rounds, with 10 timings per setting in each round. Reported times are the median of the
three round medians, which agreed within a few percent. Epicycle and Orekit were timed after
warm-up runs. GMAT's propagation time was calculated from the difference between runs with 1 and 11
propagations to exclude its start-up cost.

## Agreement between tools

The reference runs differ by the following amounts after one day.

| Pair | Position difference | Velocity difference |
|:---|---:|---:|
| Epicycle and Orekit | 0.6 mm | 0.7 µm/s |
| Epicycle and GMAT | 11.0 cm | 99.5 µm/s |

Epicycle and Orekit agree to 0.6 mm with all forces enabled. In comparisons that add the forces
one at a time, their position difference is at most 0.7 mm.

The 11.0 cm difference from GMAT comes from a difference in Earth orientation conventions:
Epicycle applies the IERS celestial-pole offsets to the Earth's axes, and GMAT's Earth-fixed frame
does not. Over one day with the 70×70 gravity field, that moves the orbit by about 11 cm. An earlier
run with Epicycle using GMAT's convention agreed with GMAT to 2.3 cm.

## Variable-step runs

Each entry gives the position error relative to that tool's own reference run, followed by the
propagation time. The methods are matched where available: GMAT uses PrinceDormand78 in place of
Dormand-Prince 8(5,3), and Orekit has no Verner method.

| Method and tolerance | Epicycle | Orekit | GMAT |
|:---|:---|:---|:---|
| Dormand-Prince 5(4), 1e-10 | 1.05 m, 0.33 s | 1.24 m, 0.71 s | 19 cm, 0.65 s |
| Dormand-Prince 8(5,3), 1e-10 | 9.2 m, 0.14 s | 25 m, 0.34 s | 16 cm, 0.44 s |
| Dormand-Prince 8(5,3), 1e-12 | 12 cm, 0.26 s | 4.1 cm, 0.62 s | 2.5 cm, 0.78 s |
| Verner 9 (GMAT RungeKutta89), 1e-12 | 0.1 mm, 0.35 s | — | 2.9 mm, 0.89 s |

Epicycle's run times at these settings are 2.0 to 3.1 times shorter. Ending its steps at the shadow
boundaries lets Verner 9 reach 0.1 mm at a tolerance of 1e-12. With Dormand-Prince 8(5,3),
Epicycle's error does not fall consistently as the tolerance tightens, even with a force model that
has no shadow, and this is under investigation.

The following settings are the fastest measured for each tool to reach errors of a few
millimetres. Epicycle has one entry for each comparison configuration.

| Tool | Setting | Position error | Time |
|:---|:---|---:|---:|
| Epicycle | Verner 7, 1e-11 | 0.2 mm | 0.25 s |
| Epicycle | Verner 9, 1e-11 | 2.0 mm | 0.29 s |
| GMAT | RungeKutta89, 1e-9 | 3.1 mm | 0.49 s |
| Orekit | Dormand-Prince 8(5,3), 1e-8 m | 8.0 mm | 1.20 s |

The Verner 7 result uses the GMAT comparison configuration and Verner 9 the Orekit comparison
configuration; for a given setting, both configurations take the same time. With shadow boundaries
in the orbit, the measured error is a better basis for choosing a setting than the tolerance alone.

## Fixed-step runs

Each tool propagates with Dormand-Prince 5(4) and a 10 s step, set up as a user of that tool would
set a fixed step:

| Tool | Setting |
|:---|:---|
| Epicycle | `DP5()`, `dt = 10`, `adaptive = false` |
| GMAT | PrinceDormand45; initial, minimum and maximum step 10 s; `ErrorControl = None` |
| Orekit | DormandPrince54; minimum and maximum step 10 s, with tolerances loose enough that no step is rejected |

| Tool | Error from own reference | Time |
|:---|---:|---:|
| Epicycle, GMAT comparison | 4.5 mm | 0.93 s |
| Epicycle, Orekit comparison | 7.3 mm | 0.95 s |
| GMAT | 2.6 mm | 2.11 s |
| Orekit | 4.8 mm | 2.11 s |

Every tool ends within a few millimetres of its reference, and Epicycle's run time is about 2.2
times shorter. Comparing final states directly, Epicycle and GMAT differ by 10.9 cm, the same Earth
orientation difference as their reference runs, and Epicycle and Orekit by 1.1 cm. The steps are not
identical even with a fixed step size, because each tool treats the shadow boundaries in its own
way: Epicycle ends 18 of its steps early at a boundary, for 8,658 steps in all, and Orekit stops at
its eclipse events. The remaining difference is each run's few millimetres of integration error, so
the reference runs above are the closer measure of how well the force models agree.
