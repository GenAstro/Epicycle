
# AstroEpochs

## Epicycle Overview

Epicycle is an application and package ecosystem for space mission design and navigation. It contains packages that handle astrodynamics models and algorithms that integrate seamlessly to allow users to setup and solve hard problems, fast.

## AstroEpochs Overview

The AstroEpochs module provides time system implementations for astronomical applications. AstroEpochs supports high-precision time representations using dual-float Julian Date storage - parameterized for differentiability - and conversions between time scales and formats.  

The API for AstroEpochs is inspired by AstroPy.time which allows for an interface that works across Epicycle components and provides type stability when time systems must change during a simulation. 

Acknowledgement: AstroEpochs uses Julia Space Mission Design's Tempo.jl for time conversion algorithms. 

**Key Features:**
- **High-precision storage** using dual Float64 values (`jd1`, `jd2`) to represent Julian Dates
- **Automatic scale conversion** via property access (e.g., `t.tt`, `t.utc`, `t.tdb`)
- **Multiple input formats** including Julian Date, Modified Julian Date, and ISO 8601 strings
- **Time arithmetic** supporting addition and subtraction of time intervals
- **Type stability** preserving numeric types through operations
- **Differentiability** using standard packages such as FiniteDiff and Zygote

## Installation

```julia
using Pkg
Pkg.add("AstroEpochs")
```

## Documentation

Full documentation is available at: [AstroEpochs Documentation](https://genastro.github.io/Epicycle/AstroEpochs/dev/)

## Comparison with Other Julia Time-Keeping Libraries

Tempo.jl and AstroTime.jl also handle astronomical time in Julia. AstroTime.jl, from the JuliaAstro community, supports six time scales (TAI, TT, TCG, TCB, TDB and UT1) with a separate type for each scale, so a conversion changes the type. Tempo.jl supports UTC, TAI, TT, TDB, TCG and TCB with allocation-free conversions and changes scale without changing the type, which Epicycle's propagation and optimization rely on for performance. AstroEpochs builds on Tempo.jl, keeps the IERS leap-second list current itself, and follows the interface of Astropy's `Time`.

## Contributing

The terms for contributions are in [LICENSE.md](LICENSE.md).

## License

AstroEpochs is open source under the [MIT License](LICENSE.md).

## Notes
Claude Sonnet and ChatGPT are used in the development of Epicycle.
