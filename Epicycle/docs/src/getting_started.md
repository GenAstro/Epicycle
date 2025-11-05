# Getting Started

## Quick Start

Get up and running with Epicycle in just a few lines of code:

```julia
using Pkg
Pkg.add("Epicycle")

using Epicycle

# Spacecraft
sat = Spacecraft(
    state=CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]),
    time=Time("2015-09-21T12:23:12", TAI(), ISOT()),
    coord_sys=CoordinateSystem(earth, ICRFAxes()),
)

# Forces + integrator
gravity = PointMassGravity(earth,(moon,sun))
forces  = ForceModel(gravity)
integ   = IntegratorConfig(Tsit5(); dt=10.0, reltol=1e-9, abstol=1e-9)
prop    = OrbitPropagator(forces, integ)

# Propagate to periapsis
propagate(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=+1))
println(get_state(sat, Keplerian()))
```
This example creates an orbit and propagates to periapis

## Installing Julia

Epicycle requires Julia 1.9 or later. Here's how to install Julia on different platforms:

### Windows
1. Download the installer from [julialang.org](https://julialang.org/downloads/)
2. Run the `.exe` installer and follow the setup wizard
3. Add Julia to your PATH when prompted
4. Verify installation by opening Command Prompt and typing `julia --version`

### macOS
**Option 1: Official Installer**
1. Download the `.dmg` file from [julialang.org](https://julialang.org/downloads/)
2. Mount the disk image and drag Julia to Applications
3. Add to PATH: `sudo ln -s /Applications/Julia-1.12.app/Contents/Resources/julia/bin/julia /usr/local/bin/julia`

**Option 2: Homebrew**
```bash
brew install julia
```

### Linux
**Ubuntu/Debian:**
```bash
# Add Julia repository
curl -fsSL https://install.julialang.org | sh
# Follow the installation script prompts
```

**From tarball (all Linux distributions):**
```bash
# Download and extract
wget https://julialang-s3.julialang.org/bin/linux/x64/1.12/julia-1.12.1-linux-x86_64.tar.gz
tar zxvf julia-1.12.1-linux-x86_64.tar.gz
sudo mv julia-1.12.1 /opt/
sudo ln -s /opt/julia-1.12.1/bin/julia /usr/local/bin/julia
```

**Verify Installation:**
```bash
julia --version
```

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

## Installing VS Code

Visual Studio Code is the recommended editor for using Epicycle. It provides excellent support for Julia through the Julia Language Server, including syntax highlighting, intelligent code completion, debugging, and integrated REPL (Like the MATLAB workspace).

### Installing VS Code

Download and install Visual Studio Code from the [official website](https://code.visualstudio.com/). Follow the platform-specific installation guides:

- **Windows**: [VS Code on Windows](https://code.visualstudio.com/docs/setup/windows)
- **macOS**: [VS Code on macOS](https://code.visualstudio.com/docs/setup/mac)
- **Linux**: [VS Code on Linux](https://code.visualstudio.com/docs/setup/linux)

### Julia Extension

After installing VS Code, add Julia support:

1. Open VS Code
2. Go to Extensions (Ctrl+Shift+X / Cmd+Shift+X)
3. Search for "Julia" and install the official Julia extension by Julia Team
4. Restart VS Code

### Configuring Julia in VS Code

The Julia extension needs to know where your Julia installation is located:

1. Open VS Code settings (Ctrl+, / Cmd+,)
2. Search for "julia executable"
3. Set the path to your Julia binary:
   - **Windows**: `C:\Users\YourName\AppData\Local\Programs\Julia-1.12.1\bin\julia.exe`
   - **macOS**: `/Applications/Julia-1.12.app/Contents/Resources/julia/bin/julia`
   - **Linux**: `/usr/local/bin/julia` (or wherever you installed Julia)

### Julia Language Server

The Julia Language Server provides advanced IDE features and is automatically installed with the Julia extension. For detailed configuration and troubleshooting:

- [Julia VS Code Documentation](https://www.julia-vscode.org/)
- [Julia Language Server Features](https://www.julia-vscode.org/docs/stable/userguide/overview/)

### Recommended VS Code Settings for Julia

Add these settings to your VS Code configuration for the best Julia experience:

```json
{
    "julia.enableTelemetry": false,
    "julia.execution.resultDisplay.repl": true,
    "julia.execution.codeInREPL": true,
    "julia.lint.run": true,
    "julia.format.indent": 4
}
```

### Getting Started with Julia in VS Code

1. Create a new file with `.jl` extension
2. Start the Julia REPL: Ctrl+Shift+P (Cmd+Shift+P) → "Julia: Start REPL"
3. Execute code: Ctrl+Enter (Cmd+Enter) to run current line/selection
4. Use F5 to run the entire file

For comprehensive guides on using Julia with VS Code, see:
- [Julia VS Code User Guide](https://www.julia-vscode.org/docs/stable/userguide/getting_started/)
- [VS Code Julia Tutorial](https://code.visualstudio.com/docs/languages/julia)  