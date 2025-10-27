# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier:

struct Sequence
    events::Vector{Event}
    adj_map::Dict{Event, Vector{Event}}
end

function Sequence()
    Sequence(Event[], Dict{Event, Vector{Event}}())
end

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

function SequenceManager(seq::Sequence)
    sorted_events = topo_sort(seq)
    ordered_vars = order_unique_vars(sorted_events)
    ordered_funcs = Constraint[]
    fun_sizes = Int[]
    found = IdSet()
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

# Utility to push unique, stateful objects preserving discovery order
@inline function _push_stateful!(out::Vector{Any}, seen::IdSet, objs)
    for obj in objs
        if is_astrosolve_stateful(typeof(obj)) && !(obj in seen)
            push!(out, obj)
            push!(seen, obj)
        end
    end
    return out
end

# Utility to find all unique stateful structs referenced by SolverVariables and events
function find_all_stateful_structs(ordered_vars, sorted_events)
    seen = IdSet()
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
    add_events(seq::Sequence, event::Event, dependencies::Vector{Event})

Add `event` to the sequence, specifying that all `dependencies` must occur before `event`.
Updates the internal adjacency map (adj_map) of the sequence.
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

function topo_sort(seq::Sequence)
    #Kahn's algorithm for topological sorting
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

function get_var_values(ordered_vars)
    vcat((get_sol_var(v) for v in ordered_vars)...)
end

function get_var_values(sm::SequenceManager)
    get_var_values(sm.ordered_vars)
end

function get_var_shifts(ordered_vars)
    vcat([v.shift for v in ordered_vars]...)
end

function get_var_scales(ordered_vars)
    vcat([v.scale for v in ordered_vars]...)
end

function get_var_lower_bounds(ordered_vars)
    vcat([v.lower_bound for v in ordered_vars]...)
end

function get_var_lower_bounds(sm::SequenceManager)
    get_var_lower_bounds(sm.ordered_vars)
end

function get_var_upper_bounds(ordered_vars)
    vcat([v.upper_bound for v in ordered_vars]...)
end

function get_var_upper_bounds(sm::SequenceManager )
    get_var_upper_bounds(sm.ordered_vars)
end

function get_fun_upper_bounds(sm::SequenceManager)
    vcat([c.upper_bounds for c in sm.ordered_funcs]...)
end 

function get_fun_lower_bounds(sm::SequenceManager)
    vcat([c.lower_bounds for c in sm.ordered_funcs]...)
end 

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
    solver_fun!(F, x, sm)

Set variables, execute the event sequence, and fill F with all constraint/objective values.
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

function copy_struct_fields!(dest, src)
    @assert typeof(dest) == typeof(src)
    for field in fieldnames(typeof(dest))
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

function reset_stateful_structs!(sm::SequenceManager)
    for (obj, init) in zip(sm.stateful_structs, sm.initial_stateful_structs)
        copy_struct_fields!(obj, init)
    end
    return nothing
end

function set_var_values(sm::SequenceManager, x::AbstractVector)
   set_var_values(x, sm.ordered_vars)
end
