# EpicycleIO

The EpicycleIO module provides plotting, reporting, and three-dimensional trajectory visualization
for Epicycle results. Output is drawn in a browser page served from the local machine.

- Data plots through Plotly, using Plotly's own attribute names and passing values through
  unchanged, so Plotly's reference documentation applies directly.
- Three-dimensional trajectory views through Cesium, built from a spacecraft's recorded history.
- A dashboard that holds several named panels on one page and updates as a script plots into it,
  including from inside a running solve.
- Formatted text reports of tabular results.

## Installation

EpicycleIO is registered in the Gen Astro registry rather than Julia's General registry, because it
is released under the Gen Astro Source Available License. Add the registry once, then install
as usual:

```julia
using Pkg
Pkg.Registry.add(RegistrySpec(url = "https://github.com/GenAstro/GenAstro.git"))
Pkg.add("EpicycleIO")
```

General is still required, since these packages depend on packages registered there.

## Usage

Full documentation is published for this package.

## Contributing

The terms for contributions are in [LICENSE.md](LICENSE.md).

## License

EpicycleIO is source available under the [Gen Astro Source Available License 1.0](LICENSE.md).
