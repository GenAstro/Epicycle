```@meta
CurrentModule = EpicycleIO
```

# Ephemeris Files

EpicycleIO reads and writes spacecraft ephemerides in the formats other tools exchange: CCSDS
Orbit Ephemeris Messages (`.oem`, versions 1.0 to 3.0 read, 3.0 and 2.0 written), STK ephemeris
files (`.e`) and JPL SPK files (`.bsp`). A file reads into a `SpacecraftHistory`, the same type a
propagation records, so a file and a propagation answer the same calls. Writing is strict: every
rule of the format is checked before a file is opened. Reading is tolerant but conservative: what
can be corrected without changing the data is corrected and reported, and the rest stops the read
with every problem listed. The readers and writers are tested against files GMAT writes, against
CSPICE's evaluation of SPK files, and by round trips through every format.

## Reading

`read_ephemeris` chooses the format from the extension. Each segment of the file becomes a
segment of the history, in the coordinate system the file names, and interpolation never crosses
from one segment to the next.

```@example eph
using Epicycle

# Write a short ephemeris to read
epoch = Time("2026-10-21T00:00:00", UTC(), ISOT())
sat = Spacecraft(state = CartesianState([7000.0, 0.0, 1300.0, 0.0, 7.35, 1.0]),
                 time = epoch,
                 name = "LEO")
forces = ForceModel(PointMassGravity(earth, ()))
integ = IntegratorConfig(DP8(); reltol = 1e-12, abstol = 1e-12)
propagate!(OrbitPropagator(forces, integ), sat, StopAt(sat, PropDurationSeconds(), 3600.0))
path = joinpath(mktempdir(), "leo.oem")
write_ephemeris(path, sat; step = 60.0)

# Read it
eph = read_ephemeris(path)
length(epochs(eph)), span(eph)
```

Beyond times and states, each segment's `metadata` holds what the file declares: the
interpolation method (`"interpolation"`), accelerations and covariance where the format carries
them (`"accelerations"`, `"covariances"`, in km, km/s, km/s² and their products), the object's
name and ID, the usable span, and the file's details (`"source"`).

`ref_frame`, `center` and `time_system` read a file as if it named different ones, for a vendor
file with a non-standard name. An SPK file holding more than one object needs `naif_id`, and an
SPK segment type that stores no samples, such as Chebyshev, is sampled through SPICE at a `step`
the caller gives, with the interpolation error between samples measured and reported.

## Evaluating

`state_at` interpolates the segment that covers a time with its declared method, or Lagrange
degree 7, and returns a `Coordinate`: the state with its coordinate system and epoch.

```@example eph
c = state_at(eph, Time("2026-10-21T00:30:30", UTC(), ISOT()))
to_vector(CartesianState(c))
```

At an epoch two segments share, such as a maneuver, the later segment applies, which is the state
after the event. A time in a gap between segments, or outside them, is an error naming the
nearest covered times; there is no extrapolation. `coverage` lists the covered intervals, since
`span` alone does not show gaps. Where two segments' usable spans overlap, `state_at` asks which
to use through `overlap = :earlier` or `:later`.

`position_difference` interpolates a second trajectory at the first one's epochs and returns the
distance at each, in km:

```@example eph
maximum(position_difference(eph, sat))
```

## Writing

`write_ephemeris` takes a `Spacecraft`, whose history it writes, or any `SpacecraftHistory`, such
as one just read. `step` resamples each segment at a fixed interval, keeping its first and last
epochs; `coordinate_system` writes in another frame, converting every state.

```@example eph
stk = joinpath(mktempdir(), "leo.e")
write_ephemeris(stk, eph; interpolation = Hermite(7))

spk = joinpath(mktempdir(), "leo.bsp")
write_ephemeris(spk, eph; naif_id = -900001)

maximum(position_difference(read_ephemeris(spk), sat))
```

Each format has rules a writer must not get wrong, and each writer states every value the format
would otherwise default:

- OEM: values are written to 16 significant digits, the most CCSDS 502.0-B-3 allows, and epochs to
  the nanosecond. A file written in version 2.0 for older tools cannot carry the 3.0-only
  `MESSAGE_ID` or `CLASSIFICATION`.
- STK: distances default to metres in the format, so `DistanceUnit` is always written.
  `InterpolationSamplesM1` is the number of points less one, so Hermite degree 7 is written as 3.
  Segments are joined by `SegmentBoundaryTimes` with the boundary epoch repeated, and the format
  has no way to show a gap, so segments must share their endpoints.
- SPK: Hermite segments are written as type 13 and Lagrange as type 9, in ICRF axes, which SPICE
  calls J2000. A source in EME2000 is rotated, not relabelled, since the two differ by the 23 mas
  frame bias. SPK holds no accelerations or covariance; a source carrying them is refused unless
  `discard_unsupported = true`. Its 64-bit ephemeris time resolves about 0.1 µs.

## Problems in a File

A read reports every problem in a file at once, each with its line. A departure from the format
that can be corrected without changing the data, such as a keyword in lower case or an exact
duplicate record, is corrected and listed in one warning. An assumption made reading a valid file,
such as an override or ITRF2000 read as Epicycle's ITRF, is recorded as a note without a warning.
`diagnostics` returns both. Anything else stops the read with an `EphemerisReadError`:

```
ERROR: EphemerisReadError: vendor.oem has 2 problems that need fixing
  line 14     REF_FRAME       "EME200" is not a known frame. Did you mean EME2000?
  line 102    data            data line has 5 values; expected 7 (epoch and state) or 10 (with accelerations)
  Fix the file, or override a field: read_ephemeris(path; ref_frame = "…", center = "…", time_system = "…")
```

`strict = true` treats every correction as an error, which is how the tests read files Epicycle
writes. `skip_bad_lines = true` drops an unreadable data line instead, splits the segment there,
and checks each piece still has enough points for its interpolation.

## Frames, Centres and Time Systems

Names follow CCSDS 502.0-B-3 and the SANA registries it makes normative. Frames read include
EME2000 (also J2000), ICRF, GCRF, ITRF and its realizations, TEME, TOD, MOD_EARTH,
J2000_ECLIPTIC, CIRS, TIRS, GTOD, MOON_ME, MOON_PA440 and a body's FIXED_CB. Centres are the Sun,
the planets, Pluto and the Moon; a barycentre or a libration point is not an Epicycle origin and
is reported as such. Time systems are UTC, TAI, TT, TDB, TCB, TCG and GPS, with ET read as TDB.

STK coordinate systems map to the same frames: J2000, ICRF, Fixed, TrueOfDate, MeanOfDate and
TEMEOfDate about the Earth. STK's Inertial depends on the central body and the STK version, and
the Moon's Fixed may be either lunar frame, so for those a read asks for `ref_frame` rather than
guess.
