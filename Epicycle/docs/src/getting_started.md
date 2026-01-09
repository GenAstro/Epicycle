# Getting Started

## Quick Start

Get up and running with Epicycle in just a few lines of code:

```julia
using Pkg
Pkg.add("Epicycle")

using Epicycle

# Spacecraft
sat = Spacecraft(
    state=KeplerianState(8000.0,0.15,pi/4,pi/2,0.0,pi/2),
    time=Time("2015-09-21T12:23:12", TAI(), ISOT()),
    coord_sys=CoordinateSystem(earth, ICRFAxes()),
    name = "sat",
)

# Forces + integrator
gravity = PointMassGravity(earth,(moon,sun))
forces  = ForceModel(gravity)
integ   = IntegratorConfig(Tsit5(); dt=10.0, reltol=1e-9, abstol=1e-9)
prop    = OrbitPropagator(forces, integ)

# Propagate to periapsis
propagate(prop, sat, StopAt(sat, PropDurationSeconds(), 5000.0))
println(get_state(sat, Keplerian()))

view = View3D()
add_spacecraft!(view,sat)
display_view(view)
```
This example creates an orbit and propagates to periapis

## Installing Julia

Epicycle requires Julia 1.10 or later. 

For installation instructions, see the [Julia Downloads page](https://julialang.org/install/). Platform-specific guides are available for Windows, macOS, and Linux.

## Installing VS Code

Visual Studio Code is the recommended editor for using Epicycle. It provides excellent support for Julia through the Julia Language Server, including syntax highlighting, intelligent code completion, debugging, and integrated REPL.

For complete installation and setup instructions, see the [VS Code Julia Tutorial](https://code.visualstudio.com/docs/languages/julia). This guide covers:

- Installing VS Code
- Installing the Julia extension
- Configuring the Julia Language Server
- Using the integrated REPL
- Debugging Julia code  

## Installing Epicycle

### From the Julia Package Registry

The easiest way to install Epicycle is through Julia's built-in package manager:

```julia
using Pkg
Pkg.add("Epicycle")
```

This will automatically install Epicycle and all its dependencies.

### Development Installation

If you want to contribute to Epicycle or need the latest development version:

```julia
using Pkg
Pkg.develop(url="https://github.com/GenAstro/Epicycle.jl")
```

### Verification

Test your installation by running:

```julia
using Epicycle

# Basic functionality test
sat = Spacecraft(
    state=CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]),
    time=Time("2015-09-21T12:23:12", TAI(), ISOT()),
    coord_sys=CoordinateSystem(earth, ICRFAxes()),
)
```

### Common Installation Issues

**Package not found:**
- Ensure you're using Julia 1.9 or later: `julia --version`
- Update your package registry: `Pkg.Registry.update()`

**Dependency conflicts:**
- Start with a fresh environment: `Pkg.activate(temp=true)`
- Try installing in isolated environment first

**Network issues:**
- If behind a corporate firewall, configure Julia's package server
- Check proxy settings in your Julia startup file

### Getting Help

If you encounter installation issues:

1. Check the [GitHub Issues](https://github.com/GenAstro/Epicycle.jl/issues) for known problems
2. Search [Julia Discourse](https://discourse.julialang.org/) for installation help
3. Open a new issue with your Julia version and error message

### Next Steps

Once installed, explore the documentation:
- [Unit Examples](unit_examples.md) - Learn specific concepts
- [Complete Examples](complete_examples.md) - See full mission simulations
- [Components](components.md) - Understand the package structure

