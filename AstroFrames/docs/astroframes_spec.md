# PackageName – Spec (v0.1)

## Development Plan

- High level scope and design 
- Prototype and test the design
- Update spec and finalize design
- Migrate prototype with refactoring to production code and implement remaining functionality.
- Test/Doc

## Todo (Acceptance Checklist)

### Phase 1 — AstroUniverse model policy (single source of truth)
- [ ] Add a global Earth frame-theory policy in AstroUniverse (`FK5()`, `IAU2006_2010()`).
- [ ] Document default policy and override behavior.

**Done when:** identical calls resolve to the same model without hidden assumptions.

### Phase 2 — EOP ownership and lifecycle
- [ ] Define AstroFrames as owner/provider of EOP data for frame transforms.
- [ ] Specify cache lifecycle (load, refresh, invalid/missing data behavior).
- [ ] Expose a stable interface for downstream consumers (AstroFrames/Epicycle internals and user API where needed).

**Done when:** transform callers do not need ad-hoc EOP handling logic.

### Phase 3 — Ambiguity spec and naming contract
- [ ] Mark which tags are model-ambiguous (e.g., `GCRF`) vs model-specific (e.g., `J2000`).
- [ ] Specify ambiguity behavior: warn/error/default + explicit override path.
- [ ] Keep intermediate systems (`CIRS`, `TIRS`) documented as internal chain elements unless explicitly requested.

**Done when:** a user can infer model semantics from API docs without reading STB internals.

### Phase 4 — Implementation alignment
- [ ] Wire AstroFrames/Epicycle model selection to AstroUniverse policy.
- [ ] Remove prototype-only ambiguity hacks once policy integration is complete.
- [ ] Add focused regression tests for Earth-fixed → Earth-inertial under both model families.

**Done when:** Earth-fixed → Earth-inertial calls are unambiguous, reproducible, and test-covered.

### Follow-up cleanup
- [ ] Fix context handling in AstroFrames and remove raw/prototype-only context pathways where no longer needed.
## 1. Overview and Scope (What are we building?)

**Purpose:**  
Provide coordinate frame transformations for astrodynamics, enabling conversion of states (position, velocity, acceleration) between different reference frames with support for both rotations and translations.

**In-scope:**
- Transform position, velocity, and acceleration vectors between coordinate frames
- Support both pure rotations and rotation+translation transformations
- Support historical (FK5) and modern (IAU 2006/2000) reference system conventions
- Provide high-performance simplified models using interpolation where appropriate
- Support common axes types:
  - Inertial frames (GCRF, J2000, etc.)
  - Earth-fixed frames (ITRF, PEF, etc.)
  - Orbit-relative frames (LVLH, TNW, etc.)
  - Ephemeris-based frames (centered on solar system bodies)
  - IAU Cartograhic coordinates for celestial bodies
- Integration with SatelliteToolboxTransformations.jl for established transformations
- Extensible design allowing users to add custom axes types
- Graph-based transformation routing (find transformation path between arbitrary frames)
- Warn users when mixing incompatible fundamental reference systems
- Fast, and differentiable

**Out-of-scope:**
- Gravitational or atmospheric models
- Numerical integration (handled by AstroProp)
  
## 2. Requirements (What we are building.)

List of abstract, testable requirements describing capabilities:
- R1: Must provide [capability]...
- R2: Must support [feature]...
- Focus on WHAT, not HOW
- No implementation details here

## 3. Design (How we are building it.)

### Core Components

The design centers around several key types, AbstarctPoint, which is in AstroBase, AbstractAxes, defined in AstroFrames, CoordinateSystem, defined in AstroFrames, and a new struct called Coordinate.

Here are the core types:

```julia
"""
    AbstractPoint

Base type for geometric points (e.g., Spacecraft, CelestialBody).
"""
abstract type AbstractPoint end
```

```julia
"""
    AbstractAxes

Abstract type for coordinate system axes definitions.

Concrete types specify the orientation of coordinate axes.
See also: [`ICRFAxes`](@ref), [`VNB`](@ref), [`Inertial`](@ref).
"""
abstract type AbstractAxes end
```

```julia
mutable struct CoordinateSystem{O<:AbstractPoint, A<:AbstractAxes} <: AbstractCoordinateSystem
    origin::O
    axes::A        
end
```


```julia

"""
    Coordinate{T<:Real}

Represents a state in space-time with position, velocity, acceleration,
time, and coordinate frame. Parametric type T supports Float64, ForwardDiff.Dual, etc.
for automatic differentiation compatibility.
"""
struct Coordinate{T<:Real}
    pos::SVector{3,T}
    vel::SVector{3,T}
    acc::SVector{3,T}
    time::Time
    frame::CoordinateFrame
end
```

**SatelliteToolboxTransformations.jl Integration**

STB provides comprehensive implementations of Earth-fixed ↔ Inertial transformations for both IAU-76/FK5 and IAU-2006/2010 models, including critical velocity transforms with angular velocity corrections (ω × r terms).

**Integration Strategy:**
- **STB provides:** ~73 Earth-fixed/Inertial frame transformations with proper velocity/acceleration handling
- **Epicycle implements:** Orbit-relative frames (LVLH, TNW, VNB), ephemeris-based frames, and custom axes types
- **Type mapping:** `ICRFAxes` ↔ `Val(:GCRF)`, `ITRFAxes` ↔ `Val(:ITRF)`, etc.
- **State conversion:** `Coordinate{T}` ↔ `OrbitStateVector` adapter layer

**AstroFrames-Owned EOP Policy and Loading API:**

EOP loading is owned by AstroFrames (not AstroUniverse) to avoid STB coupling in AstroUniverse.

Conceptual API:

```julia
struct EopPolicy
    source::Symbol          # e.g., :iers_auto, :local_file
    allow_stale::Bool
    max_age_days::Float64
end

set_eop_policy!(policy::EopPolicy)
get_eop_policy()::EopPolicy

get_eop_iau1980(; epoch=nothing)
get_eop_iau2000a(; epoch=nothing)
refresh_eop_cache!()
```

Resolution contract:
- Frame theory default (`FK5()` or `IAU2006_2010()`) comes from AstroUniverse policy.
- Required EOP data is loaded/cached by AstroFrames according to `EopPolicy`.
- Missing/stale-data behavior is enforced by AstroFrames EOP policy.

**STB Model Selection Semantics (Important):**

STB frame tags are primarily frame identifiers, not standalone declarations of a single global "fundamental system".
The effective precession-nutation/Earth-rotation model is selected by the specific transform call signature and the
source/target frame pair that STB supports.

- Do **not** assume a tag like `GCRF` always means one unique model chain in isolation.
- The same tag can participate in FK5-era or IAU-2006/2010-compatible transforms depending on the selected pair/path.
- In practice, model family is an edge/path property (transform operation), not a node-label property.
- STB documents mixed-family calls as unsupported; AstroFrames must not route through mixed-family chains implicitly.

Concrete interpretation for AstroFrames:
- `GCRF` is treated as an inertial label that can appear in different valid STB transform families.
- `J2000`, `MOD`, `TOD`, `PEF` are FK5-associated in common STB tables.
- `CIRS`, `TIRS` and CIO-style paths are IAU-2006/2010-associated in common STB tables.
- If multiple candidate paths exist, prefer a same-family direct edge over a mixed-family multi-hop.
- `J2000` is treated as FK5-specific in AstroFrames/Epicycle APIs (not model-ambiguous).

Ambiguity policy (required behavior):
- For ambiguous tags (notably `GCRF` in Earth-fixed/inertial transforms), expose or honor explicit model-theory selection.
- Default model-theory should come from AstroUniverse (`FK5()` or `IAU2006_2010()`).
- If caller behavior is ambiguous under current defaults, emit a user-visible warning with explicit override options.

Design requirement derived from this:
- Attach explicit model metadata to transform edges (e.g., `FK5`, `IAU2006_2010`) and use it in graph pruning.
- Keep public axes tags readable (`MODEq`, `TODEq`, etc.) but do not infer model-family correctness from tag names alone.

**Axes Tag Convention (Naming Table):**

To avoid ambiguity between equatorial and ecliptic frame families, use explicit suffix tags in AstroFrames names.

| Tag | Meaning | Plane Family | Typical STB/Base Mapping | Status | Notes |
|:----|:--------|:-------------|:--------------------------|:-------|:------|
| `ICRF` | International Celestial Reference Frame | Equatorial inertial | `Val(:GCRF)` (STB-compatible inertial mapping) | Prototyped | Bias-level differences may be modeled separately |
| `CIRS` | Celestial Intermediate Reference System | Equatorial inertial | `Val(:CIRS)` | Prototyped | CIO-based intermediate frame |
| `GCRF` | Geocentric Celestial Reference Frame | Equatorial inertial | `Val(:GCRF)` | Prototyped | FK5/modern compatibility depends on transform chain |
| `J2000` | Mean equator/equinox of J2000 | Equatorial inertial | `Val(:J2000)` | Prototyped | Current prototype tag; may be renamed `J2000Eq` later |
| `MOD` | Mean of Date (equatorial) | Equatorial inertial | `Val(:MOD)` | Prototyped | Current prototype tag for equatorial MOD |
| `MODEq` | Mean of Date (equatorial) | Equatorial inertial | `Val(:MOD)` | Planned naming | Preferred public-equatorial naming style |
| `TODEq` | True of Date (equatorial) | Equatorial inertial | `Val(:TOD)` | Planned | Preferred public-equatorial naming style |
| `ITRF` | International Terrestrial Reference Frame | Earth-fixed | `Val(:ITRF)` | Prototyped | Earth-fixed hub frame |
| `TIRS` | Terrestrial Intermediate Reference System | Earth-fixed | `Val(:TIRS)` | Prototyped | CIO-based Earth-fixed intermediate |
| `PEF` | Pseudo Earth Fixed | Earth-fixed | `Val(:PEF)` | Prototyped | FK5-era Earth-fixed intermediate |
| `LVLH` | Local Vertical Local Horizontal | Orbit-relative | Custom/orbit-relative transform | Prototyped | Requires spacecraft reference |
| `MODEc` | Mean of Date (ecliptic) | Ecliptic | Not directly available in STB today | Planned | Distinct from `MODEq`; do not alias |
| `TODEc` | True of Date (ecliptic) | Ecliptic | Not directly available in STB today | Planned | Distinct from `TODEq`; do not alias |

**Rule:** when a frame family has both equatorial and ecliptic variants, always include `Eq` or `Ec` in the public AstroFrames tag (`MODEq`, `TODEq`, `MODEc`, ...).

**Current STB status:** STB provides direct equatorial tags such as `MOD` and `TOD`; ecliptic tags (`MODEc`, `TODEc`) are AstroFrames-defined and currently require package-side/custom transforms.

**Design rationale:** Reimplementing STB's transformations would be error-prone and wasteful. STB handles EOP data, precession, nutation correctly. Epicycle focuses on transformations STB doesn't provide.

**Transformation Graph Design**

Coordinate frame transformations form a directed graph where:
- **Nodes** represent axes types (ICRFAxes, ITRFAxes, LVLHAxes, etc.)
- **Edges** represent known transformations between axes types
- **Paths** through the graph enable indirect transformations

**Graph Structure:**

The graph is defined statically at package initialization, specifying which direct transformations exist:

```julia
# Conceptual example - actual implementation details TBD
const AXES_GRAPH = Dict{DataType, Vector{DataType}}(
    ICRFAxes => [CIRSAxes, IAUCartographicAxes],
    CIRSAxes => [ICRFAxes, TIRSAxes],
    TIRSAxes => [CIRSAxes, ITRFAxes],
    ITRFAxes => [TIRSAxes, PlanetocentricAxes, PlanetographicAxes],
    # ... etc
)
```

**Path Finding:**

To transform from `SourceAxes` → `TargetAxes`:
1. Check if direct transformation exists (edge in graph)
2. If not, use BFS/Dijkstra to find shortest path
3. Compose transformations along the path

**Hub Nodes:**

Certain axes types act as "hubs" to reduce graph complexity:
- **ICRF**: Central hub for modern IAU transformations
- **ITRF**: Hub for Earth-fixed and surface coordinates

This allows adding new axes types by connecting to a hub, rather than defining transformations to every existing axes type.

**Extensibility:**

Users can register custom axes types:
```julia
# Add custom axes to graph
register_axes_transformation(MyCustomAxes, ICRFAxes, my_transform_function)
```

The graph automatically becomes available for path finding to any other registered axes type.

### Proposed High-Level Architecture (For Review)

This section captures the intended production design before implementation.

#### 1) Single Source of Truth: Edge Registry

- Maintain one mutable registry of direct transforms keyed by `(FromAxesType, ToAxesType)`.
- Do **not** maintain a separate hard-coded adjacency map.
- Build graph neighbors from registered edges at runtime.

Result: adding a new axes type only requires registering at least one edge.

#### 2) Dynamic Graph + Walker

- Nodes are discovered from registered edge endpoints.
- Path search uses BFS (Breadth-First Search; unweighted shortest-hop) by default.
- Future option: Dijkstra/A* if weighted routing is needed.
- Edges can expose capability requirements (e.g., EOP required, LVLH reference required).
- Walker should prune edges that are invalid for the current call context.

Result: routing remains extensible while avoiding invalid paths.

#### 3) Two Cache Layers

- **L1 `SIM_PLAN_CACHE` (simulation/session cache):**
    - Stores fully resolved transform plans for frame pairs used repeatedly in the current simulation.
    - Intended lifetime: one simulation/run context.
    - Key: `(from_axes, to_axes, edge_transform_metadata, graph_version)`.
    - Value: compiled execution plan (ordered edges + edge transform metadata snapshot).
- **L2 `PRECOMPUTED_PATH_CACHE` (precomputed and/or short path cache):**
    - Stores resolved graph paths reusable across simulations in the same Julia process.
    - Intended lifetime: process/module lifetime.
    - Key: `(from_axes, to_axes, edge_transform_metadata, graph_version)`.
    - Value: path definition (axes/edge sequence) + edge transform metadata snapshot.

Both caches are lazy-populated and read-mostly.

#### 4) Version Counter for Safe Invalidation

- Maintain `graph_version::UInt64`.
- Increment on every `register_axes_transformation!` call.
- Edge metadata is immutable after registration.
- Cache entries store the version they were built against.
- On lookup, stale entries (`entry.version != graph_version`) are ignored and rebuilt.

Result: O(1) stale-check and deterministic behavior after graph edits.

#### 5) Resolution Order (Runtime)

1. Check explicit fast-path implementation (STB direct mapping) when available.
2. Check `SIM_PLAN_CACHE` for a valid plan.
3. Check `PRECOMPUTED_PATH_CACHE` for a valid path.
4. Run graph walker with edge transform metadata filtering.
5. Store result in `PRECOMPUTED_PATH_CACHE`, then materialize/store plan in `SIM_PLAN_CACHE`.
6. Execute composed transform plan.

#### 6) API Shape (Conceptual)

```julia
register_axes_transformation!(from_axes, to_axes, edge)

transform(coord, target_frame; context=default_context())

clear_transform_caches!()            # optional admin API
graph_version()::UInt64              # diagnostic API
```

#### 7) Thread-Safety Strategy

- Registry mutations are rare; guard with a write lock.
- Read path (transform calls) should be lock-light:
    - immutable snapshots or read locks for registry view,
    - concurrent-safe cache lookups,
    - recompute on miss/stale.

This keeps hot-path performance while preserving correctness.

#### 8) Non-Goals for First Implementation

- No weighted cost model beyond shortest-hop BFS.
- No persistent on-disk transform cache.
- No automatic edge quality ranking.

These can be added incrementally after baseline dynamic graph is stable.

**Fast Path Implementations**

For Earth-fixed ↔ Inertial transformations, STB provides optimized multi-edge fast paths that bypass graph decomposition. The graph walker checks these first before attempting path finding.

**FAST_PATH_IMPLEMENTATIONS Registry:**

Maps axes type pairs directly to STB function calls:

```julia
const FAST_PATH_IMPLEMENTATIONS = Dict{Tuple{DataType, DataType}, Function}(
    (ICRFAxes, ITRFAxes) => stb_gcrf_to_itrf,
    (ITRFAxes, ICRFAxes) => stb_itrf_to_gcrf,
    (J2000Axes, ITRFAxes) => stb_j2000_to_itrf,
    # ... ~70 more STB-provided permutations
)
```

**STB Coverage Example:**

Sample from `sv_eci_to_ecef` function (see STB docstrings for complete tables):

|   Model                     |   ECI    |  ECEF  |    EOP Data     |
|:----------------------------|:---------|:-------|:----------------|
| IAU-76/FK5                  | `GCRF`   | `ITRF` | EOP IAU1980     |
| IAU-76/FK5                  | `J2000`  | `ITRF` | EOP IAU1980     |
| IAU-76/FK5                  | `MOD`    | `ITRF` | EOP IAU1980     |
| IAU-76/FK5                  | `GCRF`   | `PEF`  | EOP IAU1980     |
| IAU-76/FK5                  | `J2000`  | `PEF`  | Not required    |
| IAU-2006/2010 CIO-based     | `CIRS`   | `ITRF` | EOP IAU2000A    |
| IAU-2006/2010 CIO-based     | `GCRF`   | `ITRF` | EOP IAU2000A    |
| ...                         | ...      | ...    | ...             |

Full transformation tables available in STB docstrings for:
- `sv_eci_to_ecef` (ECI → ECEF transformations)
- `sv_ecef_to_eci` (ECEF → ECI transformations)  
- `sv_eci_to_eci` (ECI → ECI transformations)
- `sv_ecef_to_ecef` (ECEF → ECEF transformations)

**Transform Execution Order:**

1. Check `FAST_PATH_IMPLEMENTATIONS` for direct STB function
2. If found, convert `Coordinate{T}` → `OrbitStateVector`, call STB, convert back
3. If not found, use graph path finding and compose transformations
4. Cache the path for future use

**Coverage Requirement:** STB must provide all Earth-fixed/Inertial permutations shown in tables. Epicycle will NOT implement intermediate transforms for these frame types—we rely entirely on STB for Earth rotation transformations.

**Graph Walking and Performance**

Coordinate transformations happen frequently in simulations (thousands of times), but users typically use only a handful of frame pairs.

**Design:**
- **Three-level approach:** Fast paths → Cached graph paths → Fresh graph walk
- **Type-based keys:** Cache uses axes types `(FromAxes, ToAxes)` as keys
- **Lazy population:** Cache starts empty and populates via `get!`
- **What's cached:** Transformation *paths* (sequence of intermediate axes types)
- **What's computed fresh:** Transformation matrices (time-dependent: precession, nutation, EOP)

**Performance characteristics:**
- Fast path hit: Direct STB function call (~microseconds)
- Cached graph path: Dictionary lookup + compose transforms (~microseconds)
- Fresh graph walk: BFS + cache store (once per frame pair, ~milliseconds)
- Typical simulation: ~5-10 unique frame pairs, 99% fast path or cache hits

```julia
const TRANSFORM_PATH_CACHE = Dict{Tuple{DataType, DataType}, Vector{DataType}}()

function get_axes_path(from_axes::Type, to_axes::Type)
    # Check fast path first
    key = (from_axes, to_axes)
    if haskey(FAST_PATH_IMPLEMENTATIONS, key)
        return [from_axes, to_axes]  # Direct transformation
    end
    
    # Check cache
    get!(TRANSFORM_PATH_CACHE, key) do
        # Only called on cache miss - walk graph
        find_graph_path(from_axes, to_axes)
    end
end
```

**Thread safety:** Not currently addressed; add locking if multi-threaded transforms are needed.

**Types:**
- List key structs/types with brief purpose
- Note mutability where relevant
- Example: `Spacecraft` - mutable, holds state/time/properties

**Functions:**
- List key functions/methods with brief purpose
- Group by category if helpful
- Example: `propagate()` - integrates equations of motion

**Patterns:**
- Common usage patterns (constructor patterns, access patterns, etc.)
- Example: `OrbitCalc(subject, variable_type)` - bind subject to quantity

### Usage Examples

Concrete code examples demonstrating key design patterns and how they satisfy requirements. Each example should:
- Show realistic usage
- Reference which requirement(s) it addresses
- Highlight key design decisions

Example format:
```
**UC-1: [Brief Title] (Addresses R1, R3)**

Description: What this demonstrates

Code:
[actual code example]

Key Points:
- Why this pattern
- What requirement it satisfies
- Important details
```

### Design Decisions and Rationale

- Why we chose this approach
- Trade-offs made
- Alternatives considered and rejected
- Example: "Use closures for event functions because..."

### Conventions & Constraints

**Naming:**
- Any package-specific naming patterns
- Example: "Calc types named by subject: OrbitCalc, ManeuverCalc"

**Numeric Types:**
- Float64 default
- AD compatibility requirements
- Any special handling

**Mutability:**
- What's mutable, what's not
- Why (e.g., "Spacecraft mutable for in-place propagation")

**Error Handling:**
- When to throw errors vs return sentinels
- Validation points
- Error message conventions

## 4. Package Interactions (Where it fits)

- Package Dependecies:

- Used by:

- Example interactions and use cases

## 5. Testing (How we verify it)

### Test categories:

**Input Validation**
- Invalid types rejected
- Out-of-range values caught
- Inconsistent combinations detected
- Clear error messages provided

**Correctness**
- Core functionality produces expected results
- API contracts honored (e.g., mutability guarantees)
- Type stability maintained
- Return values have correct dimensions/units

**Numeric Accuracy**
- Results match reference values within specified tolerance
- Typical tolerance: 1e-12 for dimensionless, appropriate for units
- Reference sources documented (analytical, published data, validated tools)
- Known numeric edge cases handled (singularities, near-zero, etc.)

**Integration**
- Works correctly with dependent packages
- Used correctly by consuming packages
- Cross-package data flows validated
- Conversion/compatibility tested

**Regression**
- Known bugs have test coverage
- Issue number referenced in test
- Prevents reintroduction of fixed bugs

### Notes:
- Feature-specific edge cases documented in feature tests, not here
- Performance benchmarks optional, documented separately if needed