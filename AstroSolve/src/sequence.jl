# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

"""
    Sequence

A directed acyclic graph (DAG) representation of trajectory events for optimization.

In Epicycle, trajectories are modeled as sequences of events (maneuvers, propagations, etc.) 
that must occur in a specific order. The `Sequence` struct represents this as a DAG where
events are nodes and dependencies are edges.

# Fields
- `events::Vector{Event}`: Collection of all events in the sequence
- `adj_map::Dict{Event, Vector{Event}}`: Adjacency map defining event dependencies

# Examples
```julia
seq = Sequence()
add_events!(seq, event1, Event[])          # Add root event
add_events!(seq, event2, [event1])         # event2 depends on event1
add_events!(seq, event3, [event1, event2]) # event3 depends on both
```
"""
struct Sequence
    events::Vector{Event}
    adj_map::Dict{Event, Vector{Event}}
end

"""
    Sequence()

Create an empty trajectory sequence with no events or dependencies.
"""
function Sequence()
    Sequence(Event[], Dict{Event, Vector{Event}}())
end

"""
    SequenceManager

Orchestrator for trajectory sequence execution and optimization.

The `SequenceManager` prepares and manages a `Sequence` for iterative optimization.
It topologically sorts events, orders variables and constraints, extracts bounds and
scaling parameters, and manages stateful object restoration.

# Fields
- `sequence::Sequence`: Original sequence definition
- `sorted_events::Vector{Event}`: Topologically sorted events for execution
- `ordered_vars::Vector{SolverVariable}`: Unique variables in dependency order
- `ordered_funcs::Vector{Constraint}`: All constraints in execution order
- `fun_sizes::Vector{Int}`: Number of scalar constraints per constraint function
- `var_shift::Vector`: Concatenated variable shift parameters for scaling
- `var_scale::Vector`: Concatenated variable scale parameters for scaling  
- `var_lower_bounds::Vector`: Concatenated lower bounds for all variables
- `var_upper_bounds::Vector`: Concatenated upper bounds for all variables
- `stateful_structs::Vector{Any}`: All stateful objects referenced in sequence
- `initial_stateful_structs::Vector{Any}`: Deep copies for state restoration

# Examples
```julia
seq = Sequence()
add_events!(seq, event1, Event[])
add_events!(seq, event2, [event1])
sm = SequenceManager(seq)
```
"""
struct SequenceManager
    sequence::Sequence
    sorted_events::Vector{Event}
    ordered_vars::Vector{SolverVariable}
    ordered_funcs::Vector{Constraint}
    fun_sizes::Vector{Int}
    var_shift::Vector
    var_scale::Vector
    var_lower_bounds::Vector
    var_upper_bounds::Vector
    stateful_structs::Vector{Any}
    initial_stateful_structs::Vector{Any}
end

"""
    SequenceManager(seq::Sequence)

Construct a sequence manager from a trajectory sequence.

This constructor performs several optimization steps:
1. Topologically sorts events to ensure dependency ordering
2. Extracts unique variables while preserving dependency order
3. Collects all constraints with their sizes
4. Identifies and snapshots all stateful objects for reset capability
5. Pre-computes variable bounds and scaling parameters

# Notes
- `ErrorException`: If sequence contains cycles (not a valid DAG)
"""
function SequenceManager(seq::Sequence)
    sorted_events = topo_sort(seq)
    ordered_vars = order_unique_vars(sorted_events)
    ordered_funcs = Constraint[]
    fun_sizes = Int[]
    for event in sorted_events
        for c in event.funcs
            push!(ordered_funcs, c)
            push!(fun_sizes, c.numvars)
        end
    end
    stateful_structs = find_all_stateful_structs(ordered_vars, sorted_events)
    initial_stateful_structs = deepcopy(stateful_structs)
    var_shift = get_var_shifts(ordered_vars)
    var_scale = get_var_scales(ordered_vars)
    var_lower_bounds = get_var_lower_bounds(ordered_vars)
    var_upper_bounds = get_var_upper_bounds(ordered_vars)
    SequenceManager(
        seq, sorted_events, ordered_vars, ordered_funcs, fun_sizes,
        var_shift, var_scale, var_lower_bounds, var_upper_bounds,
        stateful_structs, initial_stateful_structs
    )
end

"""
    _push_stateful!(out::Vector{Any}, seen::Set{UInt}, objs)

Push unique stateful objects to output vector, preserving discovery order.

Only objects that are marked as stateful (via `is_astrosolve_stateful`) and not
already seen are added to the output collection.

# Arguments
- `out::Vector{Any}`: Output collection to append to
- `seen::Set{UInt}`: Set tracking object IDs of already discovered objects
- `objs`: Collection of objects to check and potentially add

# Returns
- `Vector{Any}`: The modified output vector (for chaining)
"""
@inline function _push_stateful!(out::Vector{Any}, seen::Set{UInt}, objs)
    for obj in objs
        if is_astrosolve_stateful(typeof(obj)) && !(objectid(obj) in seen)
            push!(out, obj)
            push!(seen, objectid(obj))
        end
    end
    return out
end

"""
    find_all_stateful_structs(ordered_vars, sorted_events)

Find all unique stateful structs referenced by SolverVariables and events.

Stateful structs are objects that maintain state during trajectory propagation and
need to be reset between optimization iterations. This function discovers all such
objects referenced throughout the sequence.

# Arguments
- `ordered_vars::Vector{SolverVariable}`: Ordered optimization variables
- `sorted_events::Vector{Event}`: Topologically sorted events

# Returns
- `Vector{Any}`: Collection of unique stateful objects in discovery order
"""
function find_all_stateful_structs(ordered_vars, sorted_events)
    seen = Set{UInt}()
    out  = Any[]

    # 1) From SolverVariables (via calc containers)
    for sv in ordered_vars
        _push_stateful!(out, seen, _subjects_from_calc(sv.calc))
    end

    # 2) From Events: vars (SolverVariables) and funcs (Constraints with calc)
    for event in sorted_events
        for sv in event.vars
            _push_stateful!(out, seen, _subjects_from_calc(sv.calc))
        end
        for con in event.funcs
            _push_stateful!(out, seen, _subjects_from_calc(con.calc))
        end
    end

    return out
end

"""
    copy_struct_fields!(dest, src)

Recursively copy all fields from `src` to `dest` struct, preserving mutable references.

This function performs deep copying for stateful struct restoration. For arrays,
it overwrites contents rather than replacing references. For nested mutable structs,
it recursively copies fields.

# Arguments
- `dest`: Destination struct to copy into
- `src`: Source struct to copy from

# Returns
- Modified `dest` struct
"""
function copy_struct_fields!(dest, src)
    @assert typeof(dest) == typeof(src)
    for field in fieldnames(typeof(dest))
        # Skip history field for Spacecraft (preserves recording flags across resets)
        if dest isa Spacecraft && field == :history
            continue
        end
        
        val = getfield(src, field)
        if typeof(val) <: AbstractArray
            # Overwrite array contents
            dest_array = getfield(dest, field)
            if length(dest_array) == length(val)
                dest_array .= val
            else
                setfield!(dest, field, copy(val))
            end
        elseif ismutable(val) && typeof(val) == typeof(getfield(dest, field))
            copy_struct_fields!(getfield(dest, field), val)
        else
            setfield!(dest, field, val)
        end
    end
    return dest
end
# TODO: Should deep copy be handled on the struct?  Seems like the better place.  

"""
    reset_stateful_structs!(sm::SequenceManager)

Reset all stateful structs in the sequence manager to their initial state.

This function is called before each optimization iteration to ensure a clean
starting state for trajectory propagation. It restores all spacecraft, maneuvers,
and other stateful objects to their initial conditions.

# Arguments
- `sm::SequenceManager`: Sequence manager containing stateful objects to reset
"""
function reset_stateful_structs!(sm::SequenceManager)
    for (obj, init) in zip(sm.stateful_structs, sm.initial_stateful_structs)
        copy_struct_fields!(obj, init)
    end
    return nothing
end

"""
    add_events!(seq::Sequence, event::Event, dependencies::Vector{Event})

Add an event to the sequence with specified dependencies.

This function adds `event` to the sequence DAG, ensuring that all events in 
`dependencies` must occur before `event` during execution. The internal adjacency 
map is updated to reflect these dependencies.

# Arguments
- `seq::Sequence`: Sequence to modify
- `event::Event`: Event to add to the sequence
- `dependencies::Vector{Event}`: Events directly linked to and preceding `event`.

# Examples
```julia
seq = Sequence()
add_events!(seq, maneuver_event, [prop_event])  # maneuver after propagation
add_events!(seq, final_event, [maneuver_event, other_event])  # final after both
```

See also: [`Sequence`](@ref), [`topo_sort`](@ref)
"""
function add_events!(seq::Sequence, event::Event, dependencies::Vector{Event})
    # Ensure all dependencies are keys in the map
    for dep in dependencies
        if !haskey(seq.adj_map, dep)
            seq.adj_map[dep] = Event[]
        end
        push!(seq.adj_map[dep], event)
    end
    # Ensure event is also a key (even if it has no dependents yet)
    if !haskey(seq.adj_map, event)
        seq.adj_map[event] = Event[]
    end
end

"""
    add_sequence!(seq::Sequence, events::Event...)

Add a linear sequence of events where each event depends on the previous one.

This is a convenience function for the common case of a linear event chain.
Each event is added only with dependencies from the event immediately before it.

# Arguments
- `seq::Sequence`: Sequence to add events to
- `events::Event...`: Events in execution order (first executes first)

# Examples
```julia
seq = Sequence()
# These two are equivalent:
add_sequence!(seq, toi_event, prop_event, moi_event)

# Equivalent to:
add_events!(seq, toi_event, Event[])
add_events!(seq, prop_event, [toi_event])
add_events!(seq, moi_event, [prop_event])
```
"""
function add_sequence!(seq::Sequence, events::Event...)
    if isempty(events)
        return
    end
    
    # Each event depends on the previous (first event implicitly has no dependencies)
    for i in 2:length(events)
        add_events!(seq, events[i], [events[i-1]])
    end
end

"""
    topo_sort(seq::Sequence)

Perform topological sorting of events using Kahn's algorithm.

This function sorts the events in the sequence. If event A depends on event B, 
then B will appear before A in the sorted output. 

# Arguments  
- `seq::Sequence`: Sequence containing events and their dependencies

# Returns
- `Vector{Event}`: Events sorted in dependency order (dependencies first)

# Throws
- `ErrorException`: If the sequence contains cycles (not a valid DAG)

# Algorithm
Uses Kahn's algorithm:
1. Calculate in-degree for all events
2. Start with events having no dependencies (in-degree 0)  
3. Remove events and update in-degrees until all processed
4. Detect cycles if any events remain unprocessed
"""
function topo_sort(seq::Sequence)

    adj_map = seq.adj_map
    # Collect all nodes (keys and all dependents)
    nodes = Set(keys(adj_map))
    for deps in values(adj_map)
        for dep in deps
            push!(nodes, dep)
        end
    end

    # Compute in-degrees
    in_degree = Dict(node => 0 for node in nodes)
    for deps in values(adj_map)
        for dep in deps
            in_degree[dep] += 1
        end
    end

    # Initialize queue with nodes of in-degree zero
    queue = [node for node in nodes if in_degree[node] == 0]
    sorted = Event[]

    while !isempty(queue)
        node = popfirst!(queue)
        push!(sorted, node)
        for dep in get(adj_map, node, Event[])
            in_degree[dep] -= 1
            if in_degree[dep] == 0
                push!(queue, dep)
            end
        end
    end

    # Check for cycles
    if length(sorted) != length(nodes)
        error("Graph has at least one cycle, topological sort not possible.")
    end

    return sorted
end

"""
    order_unique_vars(sorted_events::Vector{Event})

Extract unique solver variables from events in dependency order.

This function processes topologically sorted events and extracts all unique
SolverVariables while preserving the order of first appearance. This ensures
variables are ordered according to their dependency relationships.

# Arguments
- `sorted_events::Vector{Event}`: Events in topological order

# Returns  
- `Vector{SolverVariable}`: Unique variables in dependency order
"""
function order_unique_vars(sorted_events::Vector{Event})
    seen = Set{Any}()
    ordered_vars = SolverVariable[]
    for event in sorted_events
        for v in event.vars
            if !(v in seen)
                push!(ordered_vars, v)
                push!(seen, v)
            end
        end
    end
    return ordered_vars
end


"""
    set_var_values(sm::SequenceManager, x::AbstractVector)

Set all solver variables in the sequence manager to values from vector `x`.

# Arguments
- `sm::SequenceManager`: Sequence manager containing ordered variables
- `x::AbstractVector`: Flat vector of variable values
"""
function set_var_values(sm::SequenceManager, x::AbstractVector)
   set_var_values(x, sm.ordered_vars)
end

"""
    get_var_values(ordered_vars)

Extract current values from all SolverVariables into a flat vector.

# Arguments
- `ordered_vars::Vector{SolverVariable}`: Variables to extract values from

# Returns
- `Vector`: Concatenated values from all variables
"""
function get_var_values(ordered_vars)
    vcat((get_sol_var(v) for v in ordered_vars)...)
end

"""
    get_var_values(sm::SequenceManager)

Extract current values from all variables in the sequence manager.

# Arguments
- `sm::SequenceManager`: Sequence manager containing variables

# Returns
- `Vector`: Concatenated values from all variables
"""
function get_var_values(sm::SequenceManager)
    get_var_values(sm.ordered_vars)
end

"""
    get_var_shifts(ordered_vars)

Extract shift parameters from all solver variables into a flat vector.

# Arguments
- `ordered_vars::Vector{SolverVariable}`: Variables to extract shifts from

# Returns
- `Vector`: Concatenated shift values for optimization scaling
"""
function get_var_shifts(ordered_vars)
    vcat([v.shift for v in ordered_vars]...)
end

"""
    get_var_scales(ordered_vars)

Extract scale parameters from all solver variables into a flat vector.

# Arguments
- `ordered_vars::Vector{SolverVariable}`: Variables to extract scales from

# Returns
- `Vector`: Concatenated scale values for optimization scaling
"""
function get_var_scales(ordered_vars)
    vcat([v.scale for v in ordered_vars]...)
end

"""
    get_var_lower_bounds(ordered_vars)

Extract lower bounds from all solver variables into a flat vector.

# Arguments
- `ordered_vars::Vector{SolverVariable}`: Variables to extract bounds from

# Returns
- `Vector`: Concatenated lower bounds for optimization constraints
"""
function get_var_lower_bounds(ordered_vars)
    vcat([v.lower_bound for v in ordered_vars]...)
end

"""
    get_var_lower_bounds(sm::SequenceManager)

Extract lower bounds from all variables in the sequence manager.

# Arguments
- `sm::SequenceManager`: Sequence manager containing variables

# Returns  
- `Vector`: Concatenated lower bounds for optimization constraints
"""
function get_var_lower_bounds(sm::SequenceManager)
    get_var_lower_bounds(sm.ordered_vars)
end

"""
    get_var_upper_bounds(ordered_vars)

Extract upper bounds from all solver variables into a flat vector.

# Arguments
- `ordered_vars::Vector{SolverVariable}`: Variables to extract bounds from

# Returns
- `Vector`: Concatenated upper bounds for optimization constraints
"""
function get_var_upper_bounds(ordered_vars)
    vcat([v.upper_bound for v in ordered_vars]...)
end

"""
    get_var_upper_bounds(sm::SequenceManager)

Extract upper bounds from all variables in the sequence manager.

# Arguments
- `sm::SequenceManager`: Sequence manager containing variables

# Returns
- `Vector`: Concatenated upper bounds for optimization constraints
"""
function get_var_upper_bounds(sm::SequenceManager )
    get_var_upper_bounds(sm.ordered_vars)
end

"""
    get_fun_upper_bounds(sm::SequenceManager)

Extract upper bounds from all constraint functions in the sequence manager.

# Arguments
- `sm::SequenceManager`: Sequence manager containing constraint functions

# Returns
- `Vector`: Concatenated upper bounds for all constraints
"""
function get_fun_upper_bounds(sm::SequenceManager)
    vcat([c.upper_bounds for c in sm.ordered_funcs]...)
end 

"""
    get_fun_lower_bounds(sm::SequenceManager)

Extract lower bounds from all constraint functions in the sequence manager.

# Arguments
- `sm::SequenceManager`: Sequence manager containing constraint functions

# Returns
- `Vector`: Concatenated lower bounds for all constraints
"""
function get_fun_lower_bounds(sm::SequenceManager)
    vcat([c.lower_bounds for c in sm.ordered_funcs]...)
end 

"""
    set_var_values(vec::Vector, ordered_vars::Vector{SolverVariable})

Set values for an ordered collection of solver variables from a flat vector.

This function distributes values from a flat optimization vector to the appropriate
SolverVariables, handling the different sizes of vector vs scalar variables correctly.

# Arguments
- `vec::Vector`: Flat vector of values to distribute
- `ordered_vars::Vector{SolverVariable}`: Variables to set in order

# Examples
```julia
# Set 3 variables: 2 scalars + 1 3-vector
x = [1.0, 2.0, 0.1, 0.2, 0.3]  # 5 total values
set_var_values(x, [scalar1, scalar2, vector3])
```
"""
function set_var_values(vec::Vector, ordered_vars::Vector{SolverVariable})
    idx = 1
    for v in ordered_vars
        n = v.numvars
        subvec = vec[idx:idx+n-1]
        set_sol_var(v, subvec)
        idx += n
    end
    return nothing
end

"""
    get_fun_values(sm::SequenceManager)

Evaluate all constraint functions in the sequence manager and return 
concatenated function values.

This function evaluates all constraints in the current state and returns
their values as a flat vector. The ordering matches the constraint order
established during SequenceManager construction.

# Arguments
- `sm::SequenceManager`: Sequence manager containing constraints to evaluate

# Returns
- `Vector{Float64}`: Concatenated constraint function values

# Notes
This function assumes all stateful objects are in the correct state for
constraint evaluation (typically after event sequence execution).
"""
function get_fun_values(sm::SequenceManager)
    total_size = sum(sm.fun_sizes)
    out = Vector{Float64}(undef, total_size)
    idx = 1
    for (c, n) in zip(sm.ordered_funcs, sm.fun_sizes)
        vals = func_eval(c)
        out[idx:idx+n-1] = vals
        idx += n
    end
    return out
end

"""
    solver_fun!(F::AbstractVector, x::AbstractVector, sm::SequenceManager)

Core optimization function: set variables, execute sequence, and evaluate constraints.

This is the primary interface between optimization solvers and Epicycle trajectory
sequences. It performs a complete trajectory execution cycle:

1. **Reset**: Restore all stateful objects to initial conditions
2. **Set Variables**: Apply optimization variables to events 
3. **Execute**: Run events in topological order, applying maneuvers and propagations, etc.
4. **Evaluate**: Collect function values at appropriate times during execution

# Arguments
- `F::AbstractVector`: Output vector to fill with constraint values
- `x::AbstractVector`: Input optimization variables (flat vector)
- `sm::SequenceManager`: Configured sequence manager

# Returns
- `Int`: Status code (0 for success)

# Notes
- Function values are evaluated immediately after their associated events to capture
  the correct spacecraft state at that point in the trajectory
"""
function solver_fun!(F::AbstractVector, x::AbstractVector, sm::SequenceManager)
   
    # Reset all stateful structs to their initial state
    reset_stateful_structs!(sm)

    # Set all variables to the values passed in from solver
    set_var_values(x, sm.ordered_vars)

    # Collect constraint values as we execute events (preserves timing)
    collected_values = eltype(F)[]
    
    # Execute events and collect their constraints at the right time
    for event in sm.sorted_events
        apply_event(event)  # Modify spacecraft state
        
        # Immediately evaluate constraints for this event at current state
        for constraint in event.funcs
            vals = func_eval(constraint)  # Evaluate at current spacecraft state
            append!(collected_values, vals)
        end
    end

    # Fill F with the correctly-timed constraint values
    F[:] = collected_values
    return 0
end

"""
    solve_trajectory!(seq::Sequence; record_iterations::Bool=false)

Solve trajectory sequence using default SNOW/IPOPT configuration.

# Keyword Arguments
- `record_iterations::Bool=false`: Record solver iteration history for diagnostics.
  When `true`, iteration segments are stored in `spacecraft.history.iterations`.
  When `false` (default), only the final solution is recorded in `spacecraft.history.segments`.

# Notes
- Spacecraft history flags are automatically managed during optimization
- Original flag values are restored after optimization completes
- Final solution respects original recording settings (no forced recording)
"""
function solve_trajectory!(seq::Sequence; record_iterations::Bool=false)
    solve_trajectory!(seq, default_snow_options(); record_iterations=record_iterations)
end

"""
    solve_trajectory!(seq::Sequence, options::SNOW.Options; record_iterations::Bool=false)

Solve trajectory sequence using specified SNOW optimization options.

# Arguments
- `seq::Sequence`: Event sequence defining the trajectory optimization problem
- `options::SNOW.Options`: SNOW/IPOPT solver options

# Keyword Arguments
- `record_iterations::Bool=false`: Record solver iteration history for diagnostics.
  When `true`, iteration segments are stored in `spacecraft.history.iterations`.
  When `false` (default), only the final solution is recorded in `spacecraft.history.segments`.

# Returns
Named tuple with:
- `variables`: Optimal variable values
- `objective`: Optimal objective function value
- `constraints`: Constraint values at optimal solution
- `info`: Solver convergence information

# Notes
- Automatically manages spacecraft history recording flags during optimization
- During iterations: records to `iterations` vector if `record_iterations=true`
- After convergence: restores original flags and records final solution to `segments`
- Original flag values are preserved and restored after optimization
"""
function solve_trajectory!(seq::Sequence, options::SNOW.Options; record_iterations::Bool=false)
    # Get all spacecraft from the sequence
    sm = SequenceManager(seq)
    spacecraft = [obj for obj in sm.stateful_structs if obj isa Spacecraft]
    
    # Save original history recording flags for all spacecraft
    original_flags = [(sc.history.record_segments, sc.history.record_iterations) 
                      for sc in spacecraft]
    
    # Configure flags for optimization iterations
    for sc in spacecraft
        sc.history.record_segments = false
        sc.history.record_iterations = record_iterations
    end
    
    # Recreate SequenceManager to capture new flag settings in initial_stateful_structs
    # (reset_stateful_structs! copies from initial_stateful_structs on each iteration)
    sm = SequenceManager(seq)
    
    # Setup optimization problem
    x0 = get_var_values(sm)
    lx = get_var_lower_bounds(sm)
    ux = get_var_upper_bounds(sm)
    lg = get_fun_lower_bounds(sm)
    ug = get_fun_upper_bounds(sm)
    ng = length(lg)
    
    # Create closure for SNOW interface
    snow_solver_fun!(F, x) = solver_fun!(F, x, sm)
    
    # Run optimization
    xopt, fopt, info = minimize(snow_solver_fun!, x0, ng, lx, ux, lg, ug, options)
    
    # Restore original flags before final evaluation
    for (sc, (seg_flag, iter_flag)) in zip(spacecraft, original_flags)
        sc.history.record_segments = seg_flag
        sc.history.record_iterations = iter_flag
    end
    
    # Update initial_stateful_structs with restored flags
    # (Cannot recreate SequenceManager as it would capture post-optimization state)
    for (obj, init) in zip(sm.stateful_structs, sm.initial_stateful_structs)
        if obj isa Spacecraft
            init.history.record_segments = obj.history.record_segments
            init.history.record_iterations = obj.history.record_iterations
        end
    end
    
    # Evaluate constraints at optimal solution and record final trajectory
    # (solver_fun! will reset to initial state with correct flags via reset_stateful_structs!)
    constraint_values = Vector{Float64}(undef, ng)
    snow_solver_fun!(constraint_values, xopt)
    
    return (variables=xopt, objective=fopt, constraints=constraint_values, info=info)
end

"""
    default_snow_options()

Create default SNOW optimization options.
"""
function default_snow_options()
    # Sensible aerospace defaults
    ip_options = Dict(
        "max_iter" => 1000,
        "tol" => 1e-6,
        "file_print_level" => 0,
        "output_file" => "ipopt_$(time_ns()).out"
    )
    return Options(derivatives=ForwardFD(), solver=IPOPT(ip_options))
end

