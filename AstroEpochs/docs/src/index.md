```@meta
CurrentModule = AstroEpochs
```

# AstroEpochs

The AstroEpochs module provides time system implementations for astronomical applications. AstroEpochs supports high-precision time representations using dual-float Julian Date storage, parameterized for differentiability, and conversions between time scales and formats.

**Key Features:**
- **High-precision storage** using dual Float64 values (`jd1`, `jd2`) to represent Julian Dates
- **Automatic scale conversion** via property access (e.g., `t.tt`, `t.utc`, `t.tdb`)
- **Multiple input formats** including Julian Date, Modified Julian Date, and ISO 8601 strings
- **Time arithmetic** supporting addition and subtraction of time intervals
- **Type stability** preserving numeric types through operations
- **Differentiability** using standard packages such as FiniteDiff and Zygote

## Acknowledgements

The API for AstroEpochs is inspired by Astropy.Time. The numerics are built on Julia Space Mission Design's Tempo.jl library. AstroEpochs.jl is tested against Astropy.Time. 

## Comparison with Other Julia Time-Keeping Libraries

Tempo.jl and AstroTime.jl also handle astronomical time in Julia. AstroTime.jl, from the JuliaAstro community, supports six time scales (TAI, TT, TCG, TCB, TDB and UT1) with a separate type for each scale, so a conversion changes the type. Tempo.jl supports UTC, TAI, TT, TDB, TCG and TCB with allocation-free conversions and changes scale without changing the type, which Epicycle's propagation and optimization rely on for performance. AstroEpochs builds on Tempo.jl, keeps the IERS leap-second list current itself, and follows the interface of Astropy's `Time`. 

## Quick Start

```julia
using AstroEpochs

# Create from Julian Date
t1 = Time(2451545.0, TT(), JD())

# Create from Modified Julian Date  
t2 = Time(51544.5, UTC(), MJD())

# Create from ISO string
t3 = Time("2000-01-01T12:00:00.000", TAI(), ISOT())

# Access different representations
t1.jd        # Julian Date value
t1.mjd       # Modified Julian Date value  
t1.isot      # ISO 8601 string

# Convert between scales (creates new Time object)
t_utc = t1.utc
t_tdb = t1.tdb
```

Scale properties return new `Time` values. Format properties return a number
for `jd` and `mjd`, or a string for `isot`.

## Time Struct

The `Time` struct is the core type for representing astronomical epochs with high precision. It uses a split Julian Date representation to maintain numerical accuracy over long time spans and supports automatic conversions between different time scales and formats.

**Fields:**
- `jd1` — Primary component of the split Julian Date (typically the integer part)
- `jd2` — Secondary component of the split Julian Date (typically the fractional part)  
- `scale` — the time scale as a Symbol (`:tt`, `:tai`, `:utc`, `:tdb`, `:tcb`, `:tcg`); constructors take the tags `TT()`, `TAI()`, `UTC()`, `TDB()`, `TCB()`, `TCG()`
- `format` — the time format as a Symbol (`:jd`, `:mjd`, `:isot`); constructors take the tags `JD()`, `MJD()`, `ISOT()`

The split representation maintains precision by keeping `jd2` small (∈ [-0.5, 0.5)) while `jd1` carries the large offset. The complete Julian Date is `jd1 + jd2`.

**Precision Guidelines:**
For maximum precision, follow these best practices:
- Keep `jd2` magnitude small (< 1.0 day) to preserve floating-point precision
- Use `jd1` for large epoch offsets (e.g., set `jd1` to the integer Julian Date)
- Avoid fractional parts in `jd1` that have more than a few decimal significant figures
- The internal `_rebalance()` function automatically maintains these constraints 

## Time Scales

AstroEpochs supports various astronomical time scales.

| Scale | Description |
|:------|:------------|
| **TAI** | International Atomic Time - Uniform atomic time scale based on cesium atomic clocks |
| **TT** | Terrestrial Time - Theoretical uniform time scale for Earth-based observations (TT = TAI + 32.184s) |
| **TDB** | Barycentric Dynamical Time - Time scale for solar system dynamics, corrected for relativistic effects |
| **UTC** | Coordinated Universal Time - Civil time standard with leap seconds to maintain alignment with Earth rotation |
| **TCB** | Barycentric Coordinate Time - Coordinate time in the barycentric reference system |
| **TCG** | Geocentric Coordinate Time - Coordinate time in the geocentric reference system |

## Time Formats

AstroEpochs supports multiple time formats for input and output:

| Format | Description | Example |
|:-------|:------------|:--------|
| **JD** | Julian Date - Days since January 1, 4713 BCE at noon | 2451545.0 |
| **JD (precision)** | Julian Date with split representation for high precision | jd1=2451545.0, jd2=0.378264 |
| **MJD** | Modified Julian Date - JD minus 2400000.5 | 51544.5 |
| **ISOT** | ISO 8601 timestamp string | "2000-01-01T12:00:00.000" |

## Usage Examples

The examples below illustrate how to create a time struct in various time scales.

```julia
using AstroEpochs

# Time using TAI
t_tai = Time(51545.0, TAI(), MJD())

# Time using TDB
t_tdb = Time(51545.0, TDB(), MJD())

# Time using UTC
t_utc = Time(51545.0, UTC(), MJD())

# Time using TCB
t_tcb = Time(51545.0, TCB(), MJD())

# Time using TCG
t_tcg = Time(51545.0, TCG(), MJD())

# View all supported scales
subtypes(AstroEpochs.AbstractTimeScale)
```
Converting between time scales creates a new Time object with the converted epoch:

```julia
using AstroEpochs

# Convert from TAI to TT
t_tai = Time(51545.0, TAI(), MJD())
t_tt = t_tai.tt

# Convert from UTC to TDB
t_utc = Time(51545.0, UTC(), MJD())
t_tdb = t_utc.tdb

# Chain conversions while preserving format
t_final = t_utc.tai.tt.tdb
```

## Leap Seconds

UTC differs from TAI by a whole number of seconds, TAI − UTC, which changes when the IERS adds a
leap second. The IERS announces each one in Bulletin C, about six months ahead; TAI − UTC has
been 37 s since 2017-01-01.

AstroEpochs reads TAI − UTC from `leap-seconds.list`, the IERS list as IANA publishes it at
<https://data.iana.org/time-zones/tzdb/leap-seconds.list>. Each data line gives the date of a
change, in seconds since 1900-01-01, and the value of TAI − UTC from that date. The line that
begins `#@` gives the date the list expires, which the IERS extends each time it confirms that
no leap second is coming. The change takes effect at 0h UTC on the date given.

The file is stored in AstroEpochs' Scratch space:

```text
<depot>/scratchspaces/241dcde3-d7a4-450d-948f-15f2ea2ba1fa/leap_seconds/leap-seconds.list
```

where `<depot>` is the first entry of `DEPOT_PATH`, usually `~/.julia`. AstroEpochs downloads it
the first time a session converts to or from UTC and no stored copy exists, and again once the
stored copy has passed its expiry date. `refresh_leap_seconds!()` downloads it immediately.

Without a network connection, an expired copy is used with a warning. With no copy at all, the
table built into Tempo.jl is used, which ends at the 2017-01-01 leap second, also with a warning.
A UTC date before 1972-01-01 has no leap-second value; AstroEpochs uses 0 and warns.

```@raw html
<!-- doc-fragment -->
```
```julia
using AstroEpochs

# Where the list is stored, and its contents
path = joinpath(DEPOT_PATH[1], "scratchspaces", "241dcde3-d7a4-450d-948f-15f2ea2ba1fa",
                "leap_seconds", "leap-seconds.list")
isfile(path)                # true once a UTC conversion has run
print(read(path, String))

# Download the list now, rather than when the stored copy expires
refresh_leap_seconds!()
```

## Format Conversion

Lowercase format properties expose the same epoch in another representation.

```julia
t = Time(2451545.25, UTC(), JD())

jd = t.jd
mjd = t.mjd
timestamp = t.isot
```

Construct a new `Time` when a different stored format is required:

```julia
t_mjd = Time(t.mjd, UTC(), MJD())
```

## Precision And Arithmetic

`Time` stores a Julian Date as `jd1 + jd2`. The internal representation keeps
`jd2` near zero to retain precision while `jd1` carries the large epoch offset.
Supplying split values avoids losing a small offset when constructing a distant
epoch.

Arithmetic uses days. Adding a real number advances an epoch by that many days,
and subtracting two epochs in the same scale returns their separation in days.

```julia
t0 = Time(2451545.0, TT(), JD())
t1 = t0 + 0.5

elapsed_days = t1 - t0
```

The numeric type is preserved through construction, arithmetic, and scale
conversion so time-dependent calculations can participate in differentiation.

## API Reference

```@index
Pages = ["index.md"]
```

```@autodocs
Modules = [AstroEpochs]
Public  = true
Private = false
Order = [:type, :function, :macro, :constant]
```
