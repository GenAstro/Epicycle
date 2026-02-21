# PackageName – Spec (v0.1)

## 1. Overview and Scope (What are we building?)

Purpose:

- In-scope:
  
- Out-of-scope:
  
## 2. Requirements (What we are building.)

List of abstract, testable requirements describing capabilities:
- R1: Must provide [capability]...
- R2: Must support [feature]...
- Focus on WHAT, not HOW
- No implementation details here

## 3. Design (How we are building it.)

### Core Components

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