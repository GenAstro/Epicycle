# AstroSolve Design Specification

## Iteration Recording Integration

### Requirements

1. **Multiple Spacecraft Support**
   - A Sequence can contain multiple spacecraft
   - Events may reference different spacecraft, maneuvers, etc.
   - History flag management must handle all spacecraft in the sequence

2. **Flag Management Responsibility**
   - SequenceManager should handle history recording flag management
   - Do NOT clear existing history (user may have prior data)
   - Save original flag values and restore them after optimization
   - Apply flag settings to ALL spacecraft in the sequence

3. **Optimization Workflow**
   - **During optimization iterations:**
     - `record_segments = false` (don't accumulate solution segments)
     - `record_iterations = record_iterations` (based on user kwarg)
   - **After optimization (final solution):**
     - Restore original flag values
     - Run ONE final propagation with optimized variables
     - This populates the solution trajectory in spacecraft.history.segments

4. **API Design**
   ```julia
   solve_trajectory!(seq::Sequence; record_iterations::Bool=false)
   ```
   - Kwarg default: `false` (most users want solution only)
   - Non-breaking change (existing code works)
   - Opt-in for diagnostic iteration recording

### Design Discovery

**Key Insight:** Spacecraft are already in `sm.stateful_structs`!

From [AstroSolve.jl](c:\Users\steve\Dev\Epicycle\AstroSolve\src\AstroSolve.jl#L36-L38):
```julia
is_astrosolve_stateful(::Type{T}) where {T<:Spacecraft} = true
```

This means:
1. `find_all_stateful_structs()` already finds all spacecraft in the sequence
2. They're stored in `sm.stateful_structs` (and `sm.initial_stateful_structs`)
3. `sm` is available throughout `solve_trajectory!` - no closure problem!
4. We can filter `sm.stateful_structs` to get all `Spacecraft` instances

**Solution:** Access spacecraft through `sm.stateful_structs`, manipulate their flags directly

### Proposed Solution

```julia
function solve_trajectory!(seq::Sequence, options::SNOW.Options; record_iterations::Bool=false)
    # Get all spacecraft from the sequence (filter stateful structs)
    sm = SequenceManager(seq)
    spacecraft = [obj for obj in sm.stateful_structs if obj isa Spacecraft]
    
    # Save original flags for ALL spacecraft
    original_flags = [(sc.history.record_segments, sc.history.record_iterations) 
                      for sc in spacecraft]
    
    # Set flags for optimization iterations
    for sc in spacecraft
        sc.history.record_segments = false
        sc.history.record_iterations = record_iterations
    end
    
    # IMPORTANT: Recreate SequenceManager after flag changes
    # This ensures initial_stateful_structs captures the new flag settings
    # (reset_stateful_structs! copies from initial_stateful_structs on each iteration)
    sm = SequenceManager(seq)
    
    # Run optimization (standard flow)
    x0 = get_var_values(sm)
    lx = get_var_lower_bounds(sm)
    ux = get_var_upper_bounds(sm)
    lg = get_fun_lower_bounds(sm)
    ug = get_fun_upper_bounds(sm)
    ng = length(lg)
    
    snow_solver_fun!(F, x) = solver_fun!(F, x, sm)
    xopt, fopt, info = minimize(snow_solver_fun!, x0, ng, lx, ux, lg, ug, options)
    
    # Restore original flags before final evaluation
    for (sc, (seg_flag, iter_flag)) in zip(spacecraft, original_flags)
        sc.history.record_segments = seg_flag
        sc.history.record_iterations = iter_flag
    end
    
    # Run final evaluation with optimal variables
    # This populates spacecraft.history.segments with the solution trajectory
    constraint_values = Vector{Float64}(undef, ng)
    snow_solver_fun!(constraint_values, xopt)
    
    return (variables=xopt, objective=fopt, constraints=constraint_values, info=info)
end
```

**Key Points:**
1. Filter `sm.stateful_structs` to get `Spacecraft` instances
2. Save original flags before modification
3. **Recreate SequenceManager** after setting flags (critical!)
4. Optimization runs with iterations recording (if enabled)
5. Restore flags, run final eval to get solution trajectory

### Implementation Checklist

- [ ] Update both `solve_trajectory!` signatures with `record_iterations::Bool=false` kwarg
- [ ] Implement spacecraft filtering: `[obj for obj in sm.stateful_structs if obj isa Spacecraft]`
- [ ] Implement flag save/restore logic
- [ ] Recreate SequenceManager after flag changes (ensures initial_stateful_structs has correct flags)
- [ ] Test with single spacecraft sequence
- [ ] Test with multi-spacecraft sequence  
- [ ] Test with `record_iterations=true` (verify iteration data captured)
- [ ] Test with `record_iterations=false` (verify clean solution only)
- [ ] Verify original flags restored after optimization
- [ ] Update docstrings
- [ ] Test with Ex_GeoTransfer.jl example (should show clean trajectory, not 2035 segments)

### Open Questions & Decisions

1. **How do we know optimization is complete?**
   - **RESOLVED:** `minimize()` is a blocking call - when it returns, optimization is done
   - Everything after `minimize()` executes post-optimization
   - The final `snow_solver_fun!(constraint_values, xopt)` call happens AFTER convergence
   - This final call evaluates constraints AND triggers full propagation with optimal variables
   - Timing is perfect: restore flags, then final propagation populates solution trajectory

2. **Does the final `snow_solver_fun!(constraint_values, xopt)` call run a full propagation?**
   - **YES** - `solver_fun!` calls `apply_event` for all events in sequence
   - Events include propagation, maneuvers, etc.
   - This populates `spacecraft.history.segments` with the solution (flags restored before this call)
   - Pattern is correct: evaluate constraints + record solution trajectory in one call

3. **History clearing:**
   - Spec says "do NOT clear existing history"
   - If user runs solve multiple times, history accumulates
   - **Decision:** Leave this behavior as-is - user has reason for multiple runs
   - User can manually clear history between runs if desired
   - Consider adding `clear_history!()` helper in future if requested

4. **Flag restoration timing:**
   - Restore flags BEFORE final `snow_solver_fun!` call ✓
   - This ensures solution trajectory goes to segments (not iterations)
   - Correct as shown in Proposed Solution
