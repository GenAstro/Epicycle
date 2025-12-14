"""
    CADModel(; file_path="", scale=1.0, visible=false)

Represents a 3D CAD model for visualization of spacecraft or other objects.

# Keyword Arguments
- `file_path::String = ""` — Path to CAD model file (empty string means no model). If `visible=true`, must not be empty.
- `scale::Float64 = 1.0` — Scaling factor for the model. Must be positive (> 0.0).
- `visible::Bool = false` — Whether the model should be displayed (default hidden until model is set)

# Examples
```julia
model = CADModel(file_path="data/SpacecraftCADModel.obj", scale=10.0, visible=true)
```
"""
struct CADModel
    file_path::String
    scale::Float64
    visible::Bool
end

"""
    CADModel(; file_path="", scale=1.0, visible=false)

Kwarg constructor with defaults for all fields.
"""
function CADModel(; file_path::String="", scale::Real=1.0, visible::Bool=false)
    # Validate scale is positive
    if scale <= 0.0
        throw(ArgumentError("CADModel: scale must be positive, got $scale"))
    end
    
    # Validate that visible models have a file path
    if visible && isempty(file_path)
        throw(ArgumentError("CADModel: visible=true requires a non-empty file_path"))
    end
    
    return CADModel(file_path, Float64(scale), visible)
end

"""
    Base.show(io::IO, model::CADModel)

Pretty-print CADModel configuration.
"""
function Base.show(io::IO, model::CADModel)
    if isempty(model.file_path)
        print(io, "CADModel: (no model)")
    else
        println(io, "CADModel:")
        println(io, "  file_path = \"", model.file_path, "\"")
        println(io, "  scale     = ", model.scale)
        print(io, "  visible   = ", model.visible)
    end
end
