# AstroProp Stopping Conditions – Spec (v0.1)

## 1. Overview and Scope (What are we building?)

**Purpose:**

Stopping conditions enable termination of orbital propagation when specific events occur (state-based) or time thresholds are reached (time-based). They integrate with the propagation loop to detect events and halt integration at precise conditions.

**In-scope:**
- State-based stopping conditions (position, velocity, orbital elements)
- Time-based stopping conditions (elapsed duration, absolute epoch)
- Event detection and root-finding during propagation
- Direction-aware crossing detection (increasing, decreasing, any)
- Time scale handling (TT for Earth-centered, TDB for other centers)
- Optimization-compatible duration variables (accept positive/negative values)

**Out-of-scope:**
- Continuous event logging (only stop at first occurrence)
- Multi-event chains (stop conditions are single-target)
- Multiple time-based stopping conditions (only one time-based stop allowed per propagation)
- Event sequencing logic (that's AstroSolve)
- Custom event functions beyond StopAt pattern

## 2. Requirements (What we are building.)

**R1: State-Based Stopping**
- Must support stopping when any calc variable (position, velocity, orbital elements) crosses a target value
- Must detect crossings with specified direction (increasing, decreasing, or any)
- Must achieve target value within specified tolerance

**R2: Time-Based Duration Stopping**
- Must support stopping after elapsed time duration (seconds, days)
- Must accept positive durations (forward propagation) and negative durations (backward propagation)
- Must infer propagation direction from duration sign when direction=:infer
- Must work as optimization variables (no rejection of negative values)
- Must use appropriate dynamical time scale (TT for Earth-centered, TDB otherwise)

**R3: Absolute Time Stopping**
- Must support stopping at a specific absolute epoch (past or future)
- Must handle time scale conversions automatically
- Must infer propagation direction from time comparison when direction=:infer

**R4: Type Safety**
- Must use AbstractCalcVariable type hierarchy for compile-time validation
- Must reject invalid variable types at construction
- Must provide clear error messages for type mismatches

**R5: Integration with Propagation**
- Must work with OrdinaryDiffEq.jl callback system for state-based conditions
- Must set integration time span directly for time-based conditions (no callbacks needed)
- Must preserve propagation history when stopping

## 3. Design (How we are building it.)

### Core Components

**Types:**
- `StopAt{S,V<:AbstractCalcVariable,T}` - Generic stopping condition wrapper
  - Fields: `subject`, `var`, `target`, `direction`
  - Immutable struct holding stopping criteria
  - Type parameters ensure compile-time safety
  
- `IntegratorTimeCalc <: AbstractCalcVariable` - Abstract base for time-based stops
  - Marker type distinguishing time-based from state-based
  - Handled specially by propagate() (no callback, sets time span)
  - **Not general-purpose calcs** - propagation-specific only
  
- `PropDurationSeconds <: IntegratorTimeCalc` - Elapsed time in seconds
  - Positive value = forward propagation
  - Negative value = backward propagation
  - Zero not allowed
  
- `PropDurationDays <: IntegratorTimeCalc` - Elapsed time in days
  - Same sign semantics as PropDurationSeconds
  - Converted to seconds internally (×86400)

**Functions:**
- `StopAt(subject, var, target; direction=0)` - Main constructor
  - `direction`: -1 (decreasing), 0 (any), +1 (increasing)
  - Validates `var <: AbstractCalcVariable`
  
- `StopAt(subject::Spacecraft, target_time::Time; direction=0)` - Convenience for absolute time
  - Converts absolute time to elapsed seconds
  - Uses TT for Earth-centered, TDB otherwise (based on force model central body)
  - Supports past times with direction=:infer, errors with direction=:forward
  
- `propagate(prop, spacecraft, stop_condition; direction=:forward)` - Integration driver
  - `direction` keyword: `:forward` (default), `:backward`, or `:infer` (analyze stop condition)
  - `:infer` determines direction from stop condition (duration sign, time comparison)
  - Explicit `:forward`/`:backward` can validate agreement with stop condition

**Direction Inference:**
- For `PropDuration*`: Direction inferred from duration sign (positive=forward, negative=backward)
- For `StopAt(Time)`: Direction inferred from time comparison (target>current=forward, target<current=backward)
- For state-based: Direction defaults to `:forward` when `:infer` specified
- Explicit `:forward`/`:backward` validates against inferred direction

**Patterns:**
- **State-based**: `StopAt(subject, CalcVariable(), target_value; direction=...)`
- **Duration-based**: `StopAt(subject, PropDurationSeconds|Days(), duration)` with `direction=:infer`
- **Absolute time**: `StopAt(subject, Time(...))` with `direction=:infer`

### Usage Examples

**UC-1: State-Based Stopping (Addresses R1)**

Stop when spacecraft crosses specific position/velocity threshold.

```julia
# Stop at periapsis (r·v = 0, increasing)
stop_peri = StopAt(sat, PosDotVel(), 0.0; direction=+1)
propagate(prop, sat, stop_peri)

# Stop at radius = 7000 km (any crossing)
stop_radius = StopAt(sat, PosMag(), 7000.0)
propagate(prop, sat, stop_radius)

# Stop at x = 500 km (decreasing)
stop_x = StopAt(sat, PosX(), 500.0; direction=-1)
propagate(prop, sat, stop_x)
```

**Key Points:**
- Uses AstroCallbacks calc variables (PosDotVel, PosMag, PosX)
- Direction control for precise event detection
- Integrator finds exact crossing via root-finding

---

**UC-2: Forward Duration Stopping (Addresses R2)**

Propagate forward for specified duration.

```julia
# Propagate for 3600 seconds (1 hour) - default :forward works
stop_1hr = StopAt(sat, PropDurationSeconds(), 3600.0)
propagate(prop, sat, stop_1hr)  # direction=:forward (default)

# Propagate for 2.5 days
stop_days = StopAt(sat, PropDurationDays(), 2.5)
propagate(prop, sat, stop_days)  # Default :forward matches positive duration

# Can explicitly validate if desired
propagate(prop, sat, stop_days; direction=:forward)
```

**Key Points:**
- Positive duration with default :forward direction works naturally
- No need for :infer in common forward case
- Uses TT for Earth-centered, TDB otherwise
- No callback overhead (sets integration time span directly)

---

**UC-3: Backward Duration Stopping (Addresses R2)**

Propagate backward for specified duration.

```julia
# Propagate backward for 7200 seconds (2 hours) - use :infer or :backward
stop_back = StopAt(sat, PropDurationSeconds(), -7200.0)
propagate(prop, sat, stop_back; direction=:infer)  # Infers :backward from negative
# OR
propagate(prop, sat, stop_back; direction=:backward)  # Explicit backward

# Propagate backward for 1.5 days
stop_back_days = StopAt(sat, PropDurationDays(), -1.5)
propagate(prop, sat, stop_back_days; direction=:infer)  # Infers from sign
```

**Key Points:**
- Negative duration requires explicit direction (default :forward would conflict)
- Use `direction=:infer` to infer from sign automatically
- Or use explicit `direction=:backward`
- Spacecraft time decreases after propagation

---

**UC-4: Duration in Optimization (Addresses R2)**

Use duration as optimization variable.

```julia
# Duration is optimizable parameter
duration_calc = VariableCalc(1.5)  # Initial guess: 1.5 days
var = SolverVariable(calc=duration_calc, lower=-10.0, upper=10.0)

# In objective/constraint function - use :infer for automatic direction
function propagate_and_evaluate(duration_calc)
    dur_val = get_calc(duration_calc)
    stop = StopAt(sat, PropDurationDays(), dur_val)
    
    # direction=:infer handles positive or negative duration automatically
    propagate(prop, sat, stop; direction=:infer)
    
    # Evaluate cost/constraint...
end
```

**Key Points:**
- Optimizer can explore negative durations (backward propagation)
- `direction=:infer` automatically handles sign changes (one line, no if-test!)
- Much cleaner than manual direction checking
- No error thrown for negative values

---

**UC-5: Absolute Time Stopping (Addresses R3)**

Stop at specific epoch (past or future).

```julia
# Target time 1 day in the future - default :forward works
target_future = Time("2015-09-22T12:00:00", UTC(), ISOT())
propagate(prop, sat, StopAt(sat, target_future))  # direction=:forward (default)

# Target time in the past - use :infer or explicit :backward
target_past = Time("2015-09-20T12:00:00", UTC(), ISOT())
propagate(prop, sat, StopAt(sat, target_past); direction=:infer)  # Infers :backward
# OR
propagate(prop, sat, StopAt(sat, target_past); direction=:backward)  # Explicit

# Internally converts to elapsed seconds in correct time scale
# For Earth: uses TT
# For Mars/others: uses TDB
```

**Key Points:**
- Time scale conversion automatic (TT for Earth, TDB for others based on force model central body)
- Future times work with default :forward
- Past times work with :infer (infers :backward from comparison) or explicit :backward
- Past times error with default :forward (conflict detection)
- Direction inferred from time comparison when :infer used
- Implemented as PropDurationSeconds internally with appropriate sign

---

**UC-6: Error Cases (Addresses R4)**

Invalid configurations are rejected.

```julia
# Error: invalid direction value
stop = StopAt(sat, PropDurationDays(), 1.5)
propagate(prop, sat, stop; direction=:invalid)
# ERROR: Invalid direction :invalid, must be :infer, :forward, or :backward

# Error: explicit direction conflicts with inferred
stop = StopAt(sat, PropDurationDays(), -1.5)  # Negative → infers :backward
propagate(prop, sat, stop; direction=:forward)  # But user says :forward
# ERROR: Inferred direction is :backward (negative duration) but explicit direction is :forward

# Error: zero duration
stop = StopAt(sat, PropDurationDays(), 0.0)
propagate(prop, sat, stop)
# ERROR: Duration must be non-zero

# Error: multiple time-based stops
stop1 = StopAt(sat, PropDurationSeconds(), 3600.0)
stop2 = StopAt(sat, PropDurationDays(), 1.0)
propagate(prop, sat, stop1, stop2)
# ERROR: Multiple time-based stopping conditions not allowed

# Error: positive duration with :backward
stop = StopAt(sat, PropDurationDays(), 1.5)  # Positive → infers :forward
propagate(prop, sat, stop; direction=:backward)  # But user says :backward
# ERROR: Duration is positive (forward) but explicit direction is :backward
```

**Key Points:**
- Explicit direction validated against inferred direction (both positive→:backward and negative→:forward conflict)
- Zero duration rejected
- Multiple time-based stops rejected (ambiguous which to use)
- Clear error messages with context

---

### Design Decisions and Rationale

**Why PropDuration* instead of general elapsed time calcs?**
- Problem: Elapsed time requires anchor point, unclear who owns it
- Tried: RelativeTimeCalc(sat, ElapsedDays()) with anchor-at-construction
- Issue: Inline construction in loops creates new anchor each iteration
- Solution: PropDuration* are propagation-specific, anchor managed by propagate()
- Future: May add general RelativeTimeCalc later with explicit anchor management

**Why duration sign encodes direction with :infer?**
- Alternative 1: Always positive, direction keyword required → errors in optimization
- Alternative 2: Signed duration, manual direction keyword → boilerplate in optimization
- Solution: Signed duration + direction=:infer (default) infers from sign
- Benefits: Clean optimization (no if-test), explicit override available, works for absolute time too
- Trade-off: Slightly "magical" but ergonomics win justifies it

**Why separate StopAt(Time) constructor?**
- Common use case: "propagate to epoch"
- Convenience: Handles time scale conversion automatically
- Implementation: Converts to PropDurationSeconds internally
- Alternative considered: User manually computes elapsed time - rejected for ergonomics

**Why IntegratorTimeCalc subtypes instead of TimeCalc?**
- TimeCalc is for general time queries (scale conversions, formatting)
- PropDuration* are propagation-specific, not general-purpose
- Inheritance from AbstractCalcVariable for type safety in StopAt
- No make_calc() implementations (can't create OrbitCalc with these)

**Why TT for Earth, TDB for others?**
- Earth-based missions: TT (Terrestrial Time) is standard dynamical time
- Solar system missions: TDB (Barycentric Dynamical Time) is standard
- Consistency: Matches how integrator time is managed internally
- Detection: Based on force model central body (forces.center), not coordinate system
- Reason: Propagation dynamics are governed by the central gravitational body

### Conventions & Constraints

**Naming:**
- PropDuration* prefix signals propagation-specific scope
- Leaves ElapsedDays/Seconds/Hours available for future general calcs
- StopAt named for action (stop at event/time)

**Numeric Types:**
- Duration: Float64 (can be negative for backward)
- Target values: Float64 or Vector{Float64} (matches calc output)
- Time: Julian Date (Float64) internally
- Tolerance: Integrator-dependent (typically 1e-9 relative)

**Mutability:**
- StopAt struct is immutable (just configuration)
- Spacecraft is mutated in-place by propagate()
- Propagation history updated during integration

**Error Handling:**
- Throw errors for invalid types at construction
- Throw errors for past target times
- Throw errors for direction keyword conflicts
- Validation before propagation starts (fail fast)
- Clear error messages: "Duration is negative (-1.5) but direction is :forward"

**Duration Sign Semantics:**
- Positive: Forward propagation (time increases)
- Negative: Backward propagation (time decreases)
- Zero: Not allowed (would cause division by zero in validation)
- Validation: If direction keyword provided, must match sign

**Direction Keyword:**
- Default: `:forward` (most common case - forward propagation)
- `:infer`: Opt-in analysis of stop condition to determine direction
  - Time-based (PropDuration*): Inferred from duration sign (positive=forward, negative=backward)
  - Time-based (absolute): Inferred from time comparison (target>current=forward, target<current=backward)
  - State-based: Defaults to `:forward` when `:infer` specified
- Explicit `:backward`: Required for backward propagation with default behavior
- Explicit `:forward`/`:backward` can validate against inferred direction when using `:infer`
- Values: `:forward` (default), `:backward`, `:infer`

## 4. Package Interactions (Where it fits)

**Package Dependencies:**
- AstroCallbacks - OrbitCalc, AbstractCalcVariable hierarchy
- AstroEpochs - Time type, scale conversions
- AstroModels - Spacecraft type
- OrdinaryDiffEq - Integration and callback system

**Used by:**
- User propagation scripts
- AstroSolve - optimization with time-based variables
- Example/tutorial code

**Example Interactions:**

```julia
# AstroCallbacks provides calc variables
using AstroCallbacks: PosMag, PosDotVel, PosX

# AstroProp provides StopAt and propagate
stop_condition = StopAt(sat, PosMag(), 7000.0)
propagate(prop, sat, stop_condition)

# AstroSolve uses duration as optimization variable
duration_var = VariableCalc(2.0)
solver_var = SolverVariable(calc=duration_var, lower=-5.0, upper=5.0)
```

## 5. Testing (How we verify it)

### Test categories:

**Input Validation**
- Invalid variable types rejected at StopAt construction
- Past target time throws error
- Direction-duration sign conflicts throw error
- Zero duration throws error
- Clear error messages for all validation failures

**Correctness**
- State-based stops achieve target within tolerance
- Direction constraints honored (increasing/decreasing crossings)
- Duration-based stops advance time by correct amount
- Absolute time stops arrive at target epoch
- Backward propagation decreases time correctly

**Numeric Accuracy**
- Duration-based: elapsed time matches target ±1e-6 seconds
- State-based: target value achieved ±1e-3 km or ±1e-6 rad
- Time scale conversions accurate ±1e-9 days
- Reference: Analytical solutions, published ephemerides

**Integration**
- Works with OrbitCalc variables from AstroCallbacks
- Time conversions via AstroEpochs
- Propagation updates Spacecraft state correctly
- History segments preserved after stopping

**Regression**
- Duration sign semantics (positive=forward, negative=backward)
- Direction validation catches conflicts
- TT/TDB selection based on central body
- StopAt(Time) past-time detection

### Test Coverage:

**PropDuration Tests:**
- Forward propagation (seconds, days)
- Backward propagation with :infer (seconds, days)
- Backward propagation with explicit :backward
- Zero duration error
- Direction conflict validation (negative+:forward, positive+:backward)
- Multiple time-based stops error

**Absolute Time Tests:**
- Future time stop
- Past time error
- Time scale conversion (UTC→TT, various scales)

**Time Scale Tests:**
- Earth-centered uses TT
- Mars-centered uses TDB
- Other bodies use TDB

**State-Based Tests:**
- Radius crossing (any direction)
- Position component (increasing/decreasing)
- Orbital element crossings

### Notes:
- Tolerances reflect integrator accuracy, not event detection
- Event detection accuracy depends on callback root-finding
- PropDuration* are propagation-specific, not tested as general calcs