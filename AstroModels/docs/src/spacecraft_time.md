# Time

The `time` field stores the spacecraft's epoch using the `Time` type from AstroEpochs.

## Basic Usage

```julia
using AstroModels, AstroEpochs

# Create spacecraft with specific epoch
sc = Spacecraft(
    time = Time("2015-09-21T12:23:12", TAI(), ISOT())
)

# Access current epoch
current_time = sc.time
```

## Time Scales

Time can be specified in different time scales (TAI, UTC, TT, TDB):

```julia
# TAI (International Atomic Time)
sc = Spacecraft(
    time = Time("2015-09-21T12:23:12", TAI(), ISOT())
)

# UTC (Coordinated Universal Time)  
sc = Spacecraft(
    time = Time("2015-09-21T12:23:12", UTC(), ISOT())
)

# TT (Terrestrial Time) - for Earth-centered dynamics
sc = Spacecraft(
    time = Time("2015-09-21T12:23:12", TT(), ISOT())
)
```

## Type Promotion

When using automatic differentiation, the time's Julian day components promote to Dual numbers:

```julia
using ForwardDiff

sc = Spacecraft(
    time = Time("2015-09-21T12:23:12", TAI(), ISOT()),
    mass = ForwardDiff.Dual(1000.0, 1.0)
)

# Time components promoted to Dual
sc.time.jd1  # Dual number
sc.time.jd2  # Dual number
```
