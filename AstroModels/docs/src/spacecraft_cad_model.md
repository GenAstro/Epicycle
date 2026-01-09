# CAD Model

The `cad_model` field provides 3D visualization support for spacecraft by linking to external CAD model files.

## Basic Usage

```julia
using AstroModels

# Add a 3D model to spacecraft
sc = Spacecraft(
    cad_model = CADModel(
        file_path = "path/to/spacecraft.obj",
        scale = 100.0,
        visible = true
    )
)
```

See `CADModel` documentation for field details and default values.

The `cad_model` field stores only the file path and display settings. The visualization package (e.g., Epicycle) must support the referenced file format. See [Epicycle documentation](https://genastro.github.io/Epicycle/Epicycle/dev/) for supported formats and visualization details.

## Best Practices

### File Paths

Use relative paths from your project directory:

```julia
# Good - relative path
cad_model = CADModel(file_path = "models/satellite.obj")

# Avoid - absolute paths are not portable
cad_model = CADModel(file_path = "C:/Users/Me/models/satellite.obj")
```

## See Also

- [Spacecraft Overview](spacecraft.md) - Main spacecraft documentation
- [Reference](reference.md) - Complete API documentation
