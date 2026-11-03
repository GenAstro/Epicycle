# AstroEpochs

AstroEpochs represents astronomical epochs and converts between TAI, TT, UTC, TDB, TCB, and TCG. Its `Time` type stores a Julian date in two parts and supports Julian dates, Modified Julian Dates, and ISOT timestamps. The package manages a cached IERS leap-second table. Its API is inspired by Astropy, with algorithms adapted from ERFA and Astropy and results tested against Astropy reference data.

## Installation

```julia
using Pkg
Pkg.add("AstroEpochs")
```

## Usage

Full documentation is available on the AstroEpochs documentation site.
