
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

Tempo.jl and AstroTime.jl are other high-quality Julia packages for astronomical time handling with distinct design philosophies. AstroTime.jl, developed by the JuliaAstro community, supports six astronomical time scales (TAI, TT, TCG, TCB, TDB, and UT1) using scale-specific types that change with each conversion. Tempo.jl supports UTC, TAI, TT, TDB, TCG, and TCB with efficient, allocation-free conversions and a type-stable architecture that allows time scale changes without changing the struct type—critical for performance in Epicycle's propagation and optimization algorithms. AstroEpochs.jl is an API built on Tempo.jl that provides type stability and seamless integration with the Epicycle ecosystem while using an interface styled after AstroPy.Time, a widely adopted standard in the astronomical community. 

## Contributing to Epicycle 

Contributing is easy.

1. Fork the project
2. Create a new feature branch
3. Make your changes
4. Submit a pull request

We use the Linux Kernel's Developer's Certificate of Origin (DCO) as detailed in CONTRIBUTING.txt.

## License

We believe in the power of open source to foster innovation and community-driven 
development and also recognize the need for a sustainable business model and a model
that can handle export-controlled aerospace content. 

For these reasons, Epicycle is offered under a tri-licensing model. The license allows
users to choose between the following three options:

1) LGPL V3.0
2) Evaluation and Education use Only
3) Commercial License

See LICENSE.txt for terms each license option.  For commercial licensing, 
email licensing at genastro.org.

## Notes
Claude Sonnet and ChatGPT are used in the development of Epicycle.