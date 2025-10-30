```@meta
CurrentModule = AstroEpochs
```

# AstroEpochs

The AstroEpochs module provides time system implementations for astronomical applications. AstroEpochs supports high-precision time representations using dual-float  Julian Date storage -  parameterized for differentiability - and conversions between time scales and formats.

**Key Features:**
- **High-precision storage** using dual Float64 values (`jd1`, `jd2`) to represent Julian Dates
- **Automatic scale conversion** via property access (e.g., `t.tt`, `t.utc`, `t.tdb`)
- **Multiple input formats** including Julian Date, Modified Julian Date, and ISO 8601 strings
- **Time arithmetic** supporting addition and subtraction of time intervals
- **Type stability** preserving numeric types through operations
- **Differentiability** using standard packages such as FiniteDiff and Zygote

## Acknowledgements

The API for AstroEpochs is inspired by Astropy.Time. The numerics are built on Julia Space Mission Design's Tempo.jl library. AstroEpochs.jl is tested against Astropy.Time. 

## Quickstart

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

## Time Struct

The `Time` struct is the core type for representing astronomical epochs with high precision. It uses a split Julian Date representation to maintain numerical accuracy over long time spans and supports automatic conversions between different time scales and formats.

**Fields:**
- `jd1` — Primary component of the split Julian Date (typically the integer part)
- `jd2` — Secondary component of the split Julian Date (typically the fractional part)  
- `scale` — Time scale tag struct (`TT()`, `TAI()`, `UTC()`, `TDB()`, `TCB()`, `TCG()`)
- `format` — Time format tag struct (`JD()`, `MJD()`, `ISOT()`)

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
# Convert from TAI to TT
t_tai = Time(51545.0, TAI(), MJD())
t_tt = t_tai.tt

# Convert from UTC to TDB
t_utc = Time(51545.0, UTC(), MJD())
t_tdb = t_utc.tdb

# Chain conversions while preserving format
t_final = t_utc.tai.tt.tdb
``` 
## Time Formats

AstroEpochs supports multiple time formats for input and output:

| Format | Description | Example |
|:-------|:------------|:--------|
| **JD** | Julian Date - Days since January 1, 4713 BCE at noon UTC | 2451545.0 |
| **JD (precision)** | Julian Date with split representation for high precision | jd1=2451545.0, jd2=0.378264 |
| **MJD** | Modified Julian Date - JD minus 2400000.5 | 51544.5 |
| **ISOT** | ISO 8601 timestamp string | "2000-01-01T12:00:00.000" |

The examples below illustrate how to create time objects using different time formats.

```julia
using AstroEpochs

# Julian Date format
t_jd = Time(2451545.0, TT(), JD())

# Modified Julian Date format  
t_mjd = Time(51544.5, TT(), MJD())

# ISO 8601 string format
t_iso = Time("2000-01-01T12:00:00.000", TT(), ISOT())

# High-precision Julian Date using split representation
t_precise = Time(2451545.0, 0.37826388888889, TT(), JD())  

# View all supported formats
subtypes(AstroEpochs.AbstractTimeFormat)
```

Converting between formats (returns numeric values or strings, not new Time objects):

```julia
# Start with a time in JD format
t = Time(2451545.25, UTC(), JD())

# Access different format representations
jd_value = t.jd      # 2451545.25 (Julian Date)
mjd_value = t.mjd    # 51544.75 (Modified Julian Date)  
iso_string = t.isot  # "2000-01-01T18:00:00.000" (ISO string)

# Note: Format conversions return values, not new Time structs
# To create a new Time with different format, use the constructor
t_mjd_format = Time(t.mjd, UTC(), MJD())
```

## Time Differentiation

AstroEpochs supports automatic differentiation for time-dependent calculations using standard Julia AD packages. The `Time` struct preserves numeric types through operations, enabling differentiation of functions that depend on time.

**Using FiniteDiff.jl:**

```julia
using AstroEpochs, FiniteDiff

# Define a function that depends on time
function time_dependent_function(jd_offset)
    t = Time(2451545.0 + jd_offset, TT(), JD())
    # Convert to TDB and extract Julian Date
    return t.tdb.jd
end

# Compute derivative with respect to Julian Date offset
jd_offset = 0.5  # half day offset
derivative = FiniteDiff.finite_difference_derivative(time_dependent_function, jd_offset)
println("d(TDB)/d(JD) ≈ $derivative")
```

**Using Zygote.jl:**

```julia
using AstroEpochs, Zygote

# Create base time object outside the differentiated function
t_base = Time(2451545.0, TAI(), JD())

# Define a function that adds time in a specific scale
function add_time_in_scale(seconds_offset, scale_sym)
    t_in_scale = getproperty(t_base, scale_sym)  # Convert to desired scale
    dt_days = seconds_offset / 86400.0           # Convert seconds to days
    t_future = t_in_scale + dt_days              # Add time offset
    t_back = getproperty(t_future, t_base.scale) # Convert back to original scale
    return t_back.jd                             # Return Julian Date
end

# Compute gradient with respect to seconds added in TDB scale
seconds = 3600.0  # 1 hour in seconds
grad = Zygote.gradient(s -> add_time_in_scale(s, :tdb), seconds)
println("∇(JD)/∇(TDB_seconds) = $(grad[1])")
```
## Index

```@index
Pages = ["index.md"]
```
