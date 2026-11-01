```@meta
CurrentModule = EpicycleIO
```

# Ephemeris Files

An ephemeris describes an object's position and velocity over time. Use
[`read_ephemeris`](@ref) to load an ephemeris into Epicycle, and
[`write_ephemeris`](@ref) to export a propagated trajectory or convert an existing
file to another supported format.

Ephemerides are represented as `SpacecraftHistory` objects, the same history type
used to record a spacecraft's propagation. A history contains one or more segments,
each with its own epochs, states, coordinate system, and metadata. Segment boundaries
preserve discontinuities such as impulsive maneuvers.

## Supported formats

The filename extension selects the format. Use the `format` keyword when a filename
has a different extension or no extension.

| Format | Extension | Explicit selection | Support |
|:--|:--|:--|:--|
| CCSDS Orbit Ephemeris Message (OEM) | `.oem` | `OEM()` | Read versions 1.0–3.0 in keyword/value form; write 3.0 or 2.0 |
| STK ephemeris | `.e` | `STKEphemeris()` | Read Cartesian position, position/velocity, and position/velocity/acceleration data; write Cartesian ephemerides |
| JPL SPK | `.bsp` | `SPK()` | Read stored samples or sample a segment through SPICE; write type 9 (Lagrange) or type 13 (Hermite) |

## Read an ephemeris

```julia
using Epicycle

eph = read_ephemeris("trajectory.oem")

# Specify the format if the extension does not identify it.
eph = read_ephemeris("trajectory.txt"; format = OEM())
```

Each file segment becomes a history segment in the coordinate system specified by
the file. Positions and velocities are represented in km and km/s, regardless of
the file's distance units. Accelerations, when present, use km/s²; covariance
entries use the corresponding products of these units.

### Inspect the history

```julia
sample_times = epochs(eph)
t_start, t_stop = span(eph)
intervals = coverage(eph)
read_notes = diagnostics(eph)
```

`epochs` returns the stored sample epochs. `span` returns the beginning and end of
the overall coverage, which may contain gaps. Use `coverage` to find the intervals
where the history can be evaluated: it respects each segment's usable span and
merges intervals that touch or overlap.

Each segment's `metadata` retains the object's name and ID, the declared
interpolation method, and any accelerations, covariance records, or usable span,
under the keys `"object_name"`, `"object_id"`, `"interpolation"`, `"accelerations"`,
`"covariances"`, `"useable_start"` and `"useable_stop"`. The `"source"` entry
records file details and read diagnostics.

### Read an SPK file

Select an object with `naif_id` if the SPK contains more than one object:

```julia
eph = read_ephemeris("mission.bsp"; naif_id = -900001)
```

SPK types 9 and 13 are read using their stored samples, epochs, and interpolation
degree. For other supported types, including Chebyshev segments, provide a sampling
interval in seconds:

```julia
eph = read_ephemeris("mission.bsp"; naif_id = -900001, step = 60.0)
diagnostics(eph)
```

The reader evaluates these segments through SPICE and builds a sampled history.
It also compares the history's interpolation with SPICE at the midpoints between
samples and records the largest position difference in the diagnostics. Use this
estimate to assess the sampling interval; it is not a bound on the error at every
time. The read `step` keyword does not resample type 9 or type 13 segments.

## Evaluate and compare trajectories

[`state_at`](@ref) returns a `Coordinate` containing the interpolated state, its
coordinate system, and its epoch:

```julia
t = Time("2026-10-21T00:30:30", UTC(), ISOT())
c = state_at(eph, t)
x = to_vector(CartesianState(c))  # [x, y, z, vx, vy, vz], in km and km/s
```

The query time must fall within the history's coverage. Interpolation uses the
segment's declared method, or `Lagrange(7)` when no method is specified (lower for a
segment of fewer than eight samples), and never
crosses a segment boundary. Times in a gap or outside the coverage raise an error;
Epicycle does not extrapolate.

At a shared segment endpoint, the later segment supplies the state, so a query at
an impulsive maneuver returns the state after the event. If usable spans overlap
beyond a shared endpoint, select a segment explicitly:

```julia
c = state_at(eph, t; overlap = :later)  # or :earlier
```

Use [`position_difference`](@ref) to compare trajectories. It evaluates the second
trajectory at each of the first trajectory's sample epochs that the second covers,
converts it to the first's coordinate system, and returns the distances in km.

```julia
errors_km = position_difference(eph, sat)
maximum(errors_km)
```

## Write an ephemeris

`write_ephemeris` accepts a `Spacecraft` or a `SpacecraftHistory` and returns the
output path. For a spacecraft, it writes the recorded history; propagate the
spacecraft before exporting it. An existing output file is replaced.

The following example propagates a spacecraft for one hour and exports an OEM
with samples at 60-second intervals:

```@example eph
using Epicycle

epoch = Time("2026-10-21T00:00:00", UTC(), ISOT())
sat = Spacecraft(state = CartesianState([7000.0, 0.0, 1300.0, 0.0, 7.35, 1.0]),
                 time = epoch,
                 name = "LEO")
forces = ForceModel(PointMassGravity(earth, ()))
integ = IntegratorConfig(DP8(); reltol = 1e-12, abstol = 1e-12)
propagate!(OrbitPropagator(forces, integ), sat,
           StopAt(sat, PropDurationSeconds(), 3600.0))

output_dir = mktempdir()
path = joinpath(output_dir, "leo.oem")
write_ephemeris(path, sat; step = 60.0)

eph = read_ephemeris(path)
length(epochs(eph))
```

### Sampling, interpolation, and coordinates

By default, the writer keeps the source epochs. Set `step` to resample each segment
at a fixed interval in seconds, including its first and last epochs. The last
interval may be shorter than `step`.

The `interpolation` keyword specifies the interpolation method declared in the
output file. It defaults to the source's method, or `Lagrange(7)` if none is
recorded. Resampling evaluates the source using its existing interpolation method;
setting `interpolation` alone does not create new samples. Each segment must
contain enough samples for the selected method.

For OEM and STK, the writer uses each segment's coordinate system unless you
supply `coordinate_system`, in which case it converts the states before writing.
SPK is always written in ICRF: leave `coordinate_system` out and each segment is
rotated to ICRF; any other choice is refused. OEM epochs default to UTC; use
`time_system` to select another supported time system. STK times are UTC-based, and SPK uses TDB ephemeris time.

```julia
write_ephemeris("trajectory.oem", sat;
                step = 60.0,
                interpolation = Hermite(7),
                coordinate_system = EarthMJ2000Eq,
                time_system = UTC())
```

Segment boundaries are preserved. A redundant one-state segment recorded by a
maneuver is omitted when its state matches a neighboring segment's endpoint.
Other isolated one-state segments cannot be exported because they cannot be
interpolated.

### Convert between formats

Read a file, then write the resulting history with the desired extension:

```@example eph
stk_path = joinpath(output_dir, "leo.e")
write_ephemeris(stk_path, eph; interpolation = Hermite(7))

spk_path = joinpath(output_dir, "leo.bsp")
write_ephemeris(spk_path, eph; naif_id = -900001)

maximum(position_difference(read_ephemeris(spk_path), sat))
```

Conversion preserves supported trajectory data, but formats differ in what they
can represent. In particular:

- **OEM:** The writer defaults to version 3.0. Use `version = "2.0"` for tools that
  require it. `message_id` and `classification` require version 3.0. A new file
  does not reuse the source message ID. The source classification marking is
  retained unless explicitly overridden.
- **STK:** Segments must share their endpoints; this format cannot represent gaps
  between segments. The writer explicitly records the distance unit and converts
  the interpolation degree to STK's `InterpolationSamplesM1` convention. STK has no
  usable span, so a segment usable over less than its data is refused unless you
  pass `discard_unsupported = true`; the file then reads back with its whole data
  span as coverage. All segments must share one central body, frame, time system
  and interpolation method. STK output supports EME2000, ICRF, and Earth's ITRF,
  TOD, MOD and TEME; use `coordinate_system` to convert other frames.
- **SPK:** Supply the object's `naif_id` when writing. Hermite interpolation produces
  type 13 segments; Lagrange produces type 9. SPK output is always in ICRF axes
  (called J2000 by SPICE). States in any other frame, EME2000 and J2000_ECLIPTIC
  included, are rotated to ICRF rather than relabeled, and an ECLIPJ2000 segment
  read from a file is rotated to ICRF the same way. SPK's 64-bit ephemeris time has
  a resolution of approximately 0.1 µs in this century. The writer refuses
  position-only data, a `naif_id` that belongs to a natural body, and segments
  whose usable spans overlap; the object name becomes the SPK segment name, at
  most 40 printable ASCII characters. The reader takes J2000 and ECLIPJ2000
  segments, and a segment sampled through SPICE needs at least eight samples.

### Accelerations, covariance, and metadata

Accelerations and covariance are written where the output format supports them.
The writer rejects resampling with `step` when the source contains these data.
It also rejects a change of coordinate system for such a source, because that
conversion transforms states only.

SPK cannot carry accelerations or covariance, and STK cannot carry a usable span.
To explicitly omit those data, use `discard_unsupported = true`. For SPK they are
dropped before the states are resampled or rotated to ICRF; without the option,
writing such a source to SPK is refused.

STK and SPK cannot carry a classification marking. Exporting a marked history to
either format requires an explicit `classification = ""`.

Set `object_name` and `object_id` to identify the object (in SPK, the name becomes
the segment name), and, for OEM, `originator` to identify the file producer. Object names and IDs otherwise come from source metadata where
available; a spacecraft's name is used as a fallback object name, and missing
identifiers are written as `"UNKNOWN"`. The default originator is `"EPICYCLE"`.

For OEM, `digits` sets the decimal places of seconds in every epoch, from 0 to 9,
default 9. For STK it sets only the ScenarioEpoch; data times are measured from
it, so nothing shifts. It does not control state precision; OEM state values use
16 significant digits. A rounded OEM epoch keeps its state, so rounding would shift the
trajectory in time; the writer refuses a `digits` setting that moves any epoch by
more than 1 ns. Epochs on whole seconds can be written with `digits = 0`, while
those that are not, such as a UTC grid written in TDB, need the default.

The writer validates the output before opening the file. Validation failures raise
an `ArgumentError` and leave an existing file unchanged.

## Handle read problems

By default, readers correct supported, recoverable deviations, such as lowercase
keywords or exact duplicate records, and issue a combined warning. Assumptions,
such as an explicit metadata override, are recorded as notes. Inspect both with
`diagnostics(eph)` after a successful read.

Problems that prevent a read are collected in an `EphemerisReadError`. For text
formats, the error identifies the affected lines and fields so you can locate the
problem in the file.

| Option | Use |
|:--|:--|
| `strict = true` | Reject recoverable deviations as errors instead of correcting them. |
| `skip_bad_lines = true` | For OEM and STK, drop unreadable data rows with a warning and split the segment at each omission. Each remaining piece must have enough samples for interpolation. |
| `ref_frame`, `center` | For OEM and STK, override the declared frame or central body. |
| `time_system` | For OEM only, override the declared time system. |

Overrides change how the reader interprets the numbers; they do not transform the
states. Use them only when you know the file's declaration is incorrect or uses a
name Epicycle does not recognize. To transform correctly labeled data, use the
writer's `coordinate_system` option instead.

```julia
eph = read_ephemeris("vendor.oem"; ref_frame = "EME2000", center = "EARTH")
diagnostics(eph)
```

SPK identifies frames and centers by NAIF codes and does not accept these overrides.
STK does not accept a time-system override because its times are UTC-based.

## Frames, centers, and time systems

Supported file frame names include EME2000 (also J2000 in OEM), ICRF, GCRF, ITRF and
its realizations, TEME, TOD, MOD_EARTH, J2000_ECLIPTIC, CIRS, TIRS, GTOD, MOON_ME,
MOON_PA440, and body-fixed FIXED_CB. ITRF realizations map to Epicycle's ITRF, with
a diagnostic note describing the assumption.

Supported centers are the Sun, planets, Pluto, and Moon. Barycenters and libration
points are not supported as Epicycle origins. OEM time systems include UTC, TAI,
TT, TDB, TCB, TCG, and GPS; ET is interpreted as TDB.

STK frame names include J2000, ICRF, Fixed, TrueOfDate, MeanOfDate, and TEMEOfDate
about Earth. STK's Inertial frame depends on the central body and STK version, and
the Moon's Fixed frame is ambiguous. For those cases, the reader requires an
explicit `ref_frame` instead of choosing a frame automatically.

See the [API reference](api.md) for complete keyword descriptions.
