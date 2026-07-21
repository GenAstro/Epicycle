# AstroUniverse – Spec (v0.1)

## 1. Overview and Scope

Purpose:
- Provide canonical universe/body metadata and policy defaults used by other packages.
- Provide Earth frame-theory defaults that resolve reference-system ambiguity in frame transforms.

In-scope:
- Core body metadata (Earth, Moon, planets) and retrieval APIs.
- Theory-tag policy objects for Earth rotation/precession-nutation conventions.
- Stable integration contract for AstroFrames/Epicycle.

Out-of-scope:
- Implementing frame transforms directly (AstroFrames/Epicycle responsibility).
- Orbit propagation, force models, or numerical integration.
- Replacing STB algorithms.

## 2. Requirements

R1: Must define model policy with strongly typed tags (not Symbols).

R2: Must provide deterministic Earth model selection default for downstream packages.

R3: Must expose frame-theory policy/query API that downstream packages can consume without internal coupling.

R4: Must make ambiguous model choices observable (explicit config and inspectable active policy).

R5: Must remain dependency-light and avoid direct transform-library coupling.

R6: Must be backward-compatible at package boundary where feasible, with clear deprecation path.

## 3. Design

### 3.1 Core Types

Frame theory types (type-safe, dispatch-friendly):

```julia
abstract type AbstractFrameTheory end
struct FK5 <: AbstractFrameTheory end
struct IAU2006_2010 <: AbstractFrameTheory end
```

Earth frame policy object:

```julia
struct EarthFramePolicy{T<:AbstractFrameTheory}
	theory::T
end
```

Universe-level runtime settings:

```julia
mutable struct UniverseSettings
	earth_frame_policy::EarthFramePolicy
end
```

Notes:
- No auto theory mode in v0.1. Theory is explicit (`FK5()` or `IAU2006_2010()`).
- If an override is not provided by caller, consumers use UniverseSettings.earth_frame_policy.

### 3.2 Public API (v0.1)

User-facing configuration/query:

```julia
set_frame_theory!(theory::AbstractFrameTheory)
get_frame_theory()::AbstractFrameTheory
```

Advanced/operational API (optional exposure):

```julia
set_earth_frame_policy!(policy::EarthFramePolicy)
get_earth_frame_policy()::EarthFramePolicy
```

Integration helpers for downstream packages:

```julia
resolve_earth_theory()::AbstractFrameTheory
```

Behavior contract:
- If theory is `FK5()`, downstream transforms requiring Earth-fixed/inertial ambiguity resolution default to FK5 branch.
- If theory is `IAU2006_2010()`, downstream defaults use the STB modern CIO-based branch.
- EOP loading and cache policy are owned by AstroFrames (or transform integration layer), not AstroUniverse.

STB mapping note (must be verified precisely during implementation):
- FK5 branch corresponds to STB calls documented as IAU-76/FK5 with IAU1980 EOP inputs.
- Modern branch is currently documented in STB as IAU-2006/2010 CIO-based with IAU2000A EOP inputs.
- Final AstroUniverse naming must match verified STB behavior in code/docstrings for each transform family.

### 3.3 Usage Pattern

Startup configuration:

```julia
set_frame_theory!(IAU2006_2010())
```

Consumer resolution (AstroFrames/Epicycle conceptual flow):
1. Query `resolve_earth_theory()`.
2. Map theory to model family branch.
3. Request required EOP from AstroFrames-owned EOP APIs.
4. Execute transform via STB-backed path.

### 3.4 Design Decisions and Rationale

Decision: Frame theory is represented by concrete types, not Symbols.
- Why: type safety, dispatch clarity, typo resistance, explicit extension path.

Decision: No “auto” theory in v0.1.
- Why: avoids hidden behavior while ambiguity semantics are being stabilized.

Decision: Centralize defaults in AstroUniverse.
- Why: one source of truth for all consuming packages, deterministic results.

Decision: EOP ownership is outside AstroUniverse.
- Why: avoids direct STB coupling and keeps AstroUniverse dependency-light.

### 3.5 Conventions and Constraints

Naming:
- Frame theory uses concrete types (`FK5`, `IAU2006_2010`) under `AbstractFrameTheory`.
- Policy types end with `Policy`.

Mutability:
- Settings container is mutable for runtime configuration updates.
- Policy payload types are immutable.

Error handling:
- Invalid policy values throw ArgumentError.

## 4. Package Interactions

Depends on:
- AstroBase (shared point/type abstractions)
- AstroEpochs (epoch/time types)

Used by:
- AstroFrames: resolves default model family.
- Epicycle prototype/production transform layer: consumes frame-theory defaults.

Interaction contract:
- AstroUniverse does not execute transforms.
- AstroUniverse provides policy + data; AstroFrames/Epicycle perform transform execution.

## 5. Testing

Input validation:
- Reject unsupported theory tag types in policy setters.

Correctness:
- `set_*` then `get_*` roundtrip equality.
- `resolve_earth_theory()` returns active configured theory type.
- `set_frame_theory!` then `get_frame_theory()` roundtrip equality.

Integration:
- Mock consumer reads policy and selects expected model branch.

Regression:
- Coverage for policy persistence and migration/backward-compat behavior.

## 6. Acceptance Checklist (Step 1)

- [ ] Theory tag structs defined and exported.
- [ ] `AbstractFrameTheory`, `FK5`, and `IAU2006_2010` defined and exported.
- [ ] EarthFramePolicy implemented and query/set APIs available.
- [ ] User API `set_frame_theory!` and `get_frame_theory()` implemented.
- [ ] EarthFramePolicy implemented (internal/advanced API).
- [ ] EOP ownership delegated to AstroFrames/integration layer (documented boundary).
- [ ] AstroFrames integration contract documented and reviewed.
- [ ] Tests cover policy resolution and boundary behavior.
- [ ] STB model-family mapping verified and documented precisely (FK5 vs modern CIO-based branch).