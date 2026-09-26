# Getting Started


## Installing Julia

Epicycle requires Julia 1.10 or later. For installation instructions, see the [Julia Downloads page](https://julialang.org/install/). Platform-specific guides are available for Windows, macOS, and Linux.

## Installing VS Code

Visual Studio Code is the recommended editor for using Epicycle. It provides excellent support for Julia through the Julia Language Server, including syntax highlighting, intelligent code completion, debugging, and integrated REPL.

For complete installation and setup instructions, see the [VS Code Julia Tutorial](https://code.visualstudio.com/docs/languages/julia). This guide covers:

- Installing VS Code
- Installing the Julia extension
- Configuring the Julia Language Server
- Using the integrated REPL
- Debugging Julia code  

## Installing Epicycle

### From the Gen Astro registry

To install the latest version of Epicycle, first add the local registry (the app store, for those unfamiliar with Julia), then install as usual:

```julia
using Pkg
Pkg.Registry.add(
    RegistrySpec(url = "https://github.com/GenAstro/GenAstro.git")
)
Pkg.add("Epicycle")
```

!!! note
    Some packages originally registered in the Julia General registry, including AstroModels, AstroProp, and AstroSolve, have moved to the GenAstro local registry. If you do not add the local registry as shown above, you will install only the first MVP release of Epicycle.


The first `using Epicycle` also downloads about 110 MB of SPICE kernels, most of it the DE440
ephemeris that planet and Moon positions come from. That first load needs a network connection and
waits for the download to finish. Later sessions use the stored copy and work offline.

### From source

The repository is a monorepo. The umbrella `Epicycle` package and each component package sit in
their own subdirectory, so a development install names the subdirectory it wants:

```julia
using Pkg
Pkg.develop(url = "https://github.com/GenAstro/Epicycle.git", subdir = "Epicycle")
```

To work on more than one package, clone the repository yourself and develop each one by path, so
that all of them resolve to your clone rather than to the registry:

```julia
Pkg.develop(path = "Epicycle/AstroProp")
```

### Verification

Confirm the installation by building a spacecraft:

```julia
using Epicycle

sat = Spacecraft(
    state=CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]),
    time=Time("2015-09-21T12:23:12", TAI(), ISOT()),
    coord_sys=CoordinateSystem(earth, ICRF()),
)
```

### Common Installation Issues

`Pkg.add("Epicycle")` reporting that the package is not found almost always means the Gen Astro
registry has not been added. Run the `Pkg.Registry.add` line above, then `Pkg.Registry.update()`,
and try again.

Epicycle requires Julia 1.10 or later. On an older release Pkg reports that the julia version
requirement is not satisfied rather than that the package is missing; check with `julia --version`.

A dependency conflict in an environment that already holds other packages is usually a version
bound Epicycle cannot meet. Installing into a fresh environment with `Pkg.activate(temp=true)`
separates that case from a broken installation.

Behind a corporate firewall, both the registry clone and the SPICE kernel download go out over
HTTPS. Configure Julia's package server and proxy settings in your startup file if either times
out.

### Getting Help

If you encounter installation issues:

1. Check the [GitHub Issues](https://github.com/GenAstro/Epicycle/issues) for known problems
2. Search [Julia Discourse](https://discourse.julialang.org/) for installation help
3. Open a new issue with your Julia version and error message

## Running Examples

Epicycle provides a library runnable examples covering propagation, targeting, optimal control and orbit
determination.  The code below shows how to run the "getting started" example, and how to get the names of
all examples to run others in the suite.

```julia
# Run the example named "Ex_GettingStarted"
using Epicycle
Epicycle.run_example("Ex_GettingStarted")

# Print the names of all examples
Epicycle.list_examples()
```

## Next Steps

Once installed, explore the documentation:
- [Examples](examples/Ex_PropagationBasics.md) - Runnable scripts, from propagation to optimal control
- [Packages](packages.md) - Understand the package structure

