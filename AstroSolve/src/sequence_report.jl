"""
    report_sequence(seq::Sequence)

Generate a comprehensive report summarizing the trajectory sequence configuration.

This function creates a human-readable summary of the sequence structure showing:
- Sequence overview (total events, variables, constraints)
- Event details with dependencies and execution order
- Variable information (names, sizes, bounds)
- Constraint summary (types, sizes, bounds)

# Arguments
- `seq::Sequence`: Trajectory sequence to analyze

# Returns
- Nothing (prints directly to stdout for better formatting)
"""
function report_sequence(seq::Sequence)
    # Create sequence manager to get ordered information
    sm = SequenceManager(seq)
    
    println()
    println("TRAJECTORY SEQUENCE SUMMARY")
    println("="^50)
    println()
    
    # Overview
    println("Sequence Overview:")
    println("- Total Events: $(length(sm.sorted_events))")
    total_vars = sum(v.numvars for v in sm.ordered_vars; init=0)
    println("- Variable Objects: $(length(sm.ordered_vars)) ($(total_vars) optimization variables)")
    total_constraints = sum(sm.fun_sizes; init=0)
    println("- Constraint Objects: $(length(sm.ordered_funcs)) ($(total_constraints) constraint functions)")
    
    # Execution order summary
    event_names = ["\"$(event.name)\"" for event in sm.sorted_events]
    if length(event_names) <= 3
        order_str = join(event_names, " → ")
    else
        order_str = join(event_names[1:2], " → ") * " → ... → " * event_names[end]
    end
    println("- Execution Order: [$order_str]")
    println()
    
    # Event Details
    println("EVENT DETAILS:")
    println("-"^20)
    println()
    
    for (i, event) in enumerate(sm.sorted_events)
        println("Event $i: \"$(event.name)\"")
        
        # Variables section
        if !isempty(event.vars)
            total_event_vars = sum(v.numvars for v in event.vars; init=0)
            println("├─ Variable Objects ($(length(event.vars))): $(total_event_vars) optimization variables")
            for (j, var) in enumerate(event.vars)
                is_last_var = (j == length(event.vars))
                prefix = is_last_var && isempty(event.funcs) ? "│  └─ " : "│  ├─ "
                
                # Extract variable information
                var_desc = get_enhanced_variable_description(var)
                println("$prefix$var_desc")
                
                # Show component breakdown for multi-component variables
                if var.numvars > 1
                    show_variable_components(var, is_last_var && isempty(event.funcs))
                end
            end
        else
            println("├─ Variable Objects (0): None")
        end
        
        # Constraints section
        if !isempty(event.funcs)
            total_event_constraints = sum(f.numvars for f in event.funcs; init=0)
            println("└─ Constraint Objects ($(length(event.funcs))): $(total_event_constraints) constraint functions")
            for (j, func) in enumerate(event.funcs)
                is_last = (j == length(event.funcs))
                prefix = is_last ? "   └─ " : "   ├─ "
                
                # Extract constraint information
                constraint_desc = get_enhanced_constraint_description(func)
                println("$prefix$constraint_desc")
                
                # Show component breakdown for multi-component constraints
                if func.numvars > 1
                    show_constraint_components(func, is_last)
                end
            end
        else
            if isempty(event.vars)
                println("└─ Constraint Objects (0): None")
            else
                println("└─ Constraint Objects (0): None")
            end
        end
        
        println()
    end
    
    # Stateful Objects Summary
    if !isempty(sm.stateful_structs)
        println("STATEFUL OBJECTS:")
        println("-"^20)
        
        # Group by type for cleaner display
        type_counts = Dict{String, Int}()
        for obj in sm.stateful_structs
            type_name = get_simplified_type_name(obj)
            type_counts[type_name] = get(type_counts, type_name, 0) + 1
        end
        
        for (type_name, count) in sort(collect(type_counts))
            if count == 1
                println("- $type_name")
            else
                println("- $type_name (×$count)")
            end
        end
        println()
    end
    
    println("="^50)
    return nothing
end

"""
    format_bounds(lower, upper, size=1)

Format numeric bounds for display with consistent precision.
"""
function format_bounds(lower, upper, size=1)
    if size == 1
        return "[$(round(lower, digits=3)), $(round(upper, digits=3))]"
    else
        return "$(size)×[$(round(lower, digits=3)), $(round(upper, digits=3))]"
    end
end

"""
    get_enhanced_variable_description(var::SolverVariable)

Get enhanced variable description with calc info for display.
"""
function get_enhanced_variable_description(var::SolverVariable)
    # Get variable name
    var_name = hasfield(typeof(var), :name) ? var.name : "unnamed"
    
    # Extract calc information
    calc_desc = get_calc_description(var.calc)
    
    # Get bounds description for single component
    if var.numvars == 1
        bounds_desc = get_bounds_description(var.lower_bound[1], var.upper_bound[1])
        return "$var_name: $calc_desc $bounds_desc"
    else
        return "$var_name: $calc_desc ($(var.numvars) components)"
    end
end

"""
    get_enhanced_constraint_description(func)

Get enhanced constraint description with calc info for display.
"""
function get_enhanced_constraint_description(func)
    # Extract calc information (same as variables)
    calc_desc = get_calc_description(func.calc)
    
    # Get bounds description for single component
    if func.numvars == 1
        bounds_desc = get_bounds_description(func.lower_bounds[1], func.upper_bounds[1])
        return "$calc_desc $bounds_desc"
    else
        return "$calc_desc ($(func.numvars) components)"
    end
end

"""
    get_calc_description(calc)

Extract calculation description from calc objects.
"""
function get_calc_description(calc)
    try
        calc_type = string(typeof(calc))
        calc_type = replace(calc_type, r"^.*\.([^.]+)$" => s"\1")  # Remove module prefix
        calc_type = replace(calc_type, r"\{.*\}" => "")  # Remove type parameters
        
        # Try to extract variable type from calc.var field (more reliable)
        calc_fields = fieldnames(typeof(calc))
        
        if :var in calc_fields
            var_obj = getfield(calc, :var)
            var_name = string(typeof(var_obj))
            var_name = replace(var_name, r"^.*\.([^.]+)$" => s"\1")
            var_name = replace(var_name, r"\{.*\}" => "")
            var_name = replace(var_name, r"\(\)" => "")  # Remove empty parentheses
            return "$var_name() ($calc_type)"
        else
            return calc_type
        end
    catch e
        # If anything fails, fall back to basic type name
        calc_type = string(typeof(calc))
        calc_type = replace(calc_type, r"^.*\.([^.]+)$" => s"\1")
        calc_type = replace(calc_type, r"\{.*\}" => "")
        return calc_type
    end
end

"""
    get_subject_description(subject)

Get subject description with type and name information.
"""
function get_subject_description(subject)
    try
        subject_type = string(typeof(subject))
        subject_type = replace(subject_type, r"^.*\.([^.]+)$" => s"\1")
        subject_type = replace(subject_type, r"\{.*\}" => "")
        
        # Try to get subject name/ID if available
        if hasfield(typeof(subject), :name)
            return "$subject_type(\"$(subject.name)\")"
        elseif hasfield(typeof(subject), :id)
            return "$subject_type(\"$(subject.id)\")"
        else
            return subject_type
        end
    catch
        return "UnknownSubject"
    end
end

"""
    get_bounds_description(lower, upper)

Describe bounds as equality, inequality, or range constraints.
"""
function get_bounds_description(lower, upper)
    tol = 1e-12
    if abs(lower - upper) < tol
        return "= $(round(upper, digits=6))"
    elseif lower == -Inf && upper != Inf
        return "≤ $(round(upper, digits=6))"
    elseif lower != -Inf && upper == Inf
        return "≥ $(round(lower, digits=6))"
    else
        return "∈ [$(round(lower, digits=6)), $(round(upper, digits=6))]"
    end
end

"""
    show_variable_components(var::SolverVariable, is_last_var::Bool)

Show variable component breakdown with proper indentation.
"""
function show_variable_components(var::SolverVariable, is_last_var::Bool)
    base_prefix = is_last_var ? "│     " : "│  │  "
    
    for i in 1:var.numvars
        is_last_component = (i == var.numvars)
        comp_prefix = is_last_component ? "└─ " : "├─ "
        bounds_desc = get_bounds_description(var.lower_bound[i], var.upper_bound[i])
        println("$base_prefix$comp_prefix Component $i: $bounds_desc")
    end
end

"""
    show_constraint_components(func, is_last_constraint::Bool)

Show constraint component breakdown with proper indentation.
"""
function show_constraint_components(func, is_last_constraint::Bool)
    base_prefix = is_last_constraint ? "      " : "   │  "
    
    for i in 1:func.numvars
        is_last_component = (i == func.numvars)
        comp_prefix = is_last_component ? "└─ " : "├─ "
        bounds_desc = get_bounds_description(func.lower_bounds[i], func.upper_bounds[i])
        println("$base_prefix$comp_prefix Component $i: $bounds_desc")
    end
end

"""
    get_simplified_type_name(obj)

Simplify type names by removing module prefixes and parameters.
"""
function get_simplified_type_name(obj)
    type_name = string(typeof(obj))
    # Remove module prefixes and type parameters
    type_name = replace(type_name, r"^.*\.([^.]+)$" => s"\1")
    type_name = replace(type_name, r"\{.*\}" => "")
    return type_name
end

"""
    report_solution(seq::Sequence, result)

Generate a comprehensive report showing the optimized trajectory solution.

This function displays the optimization results including variable values,
constraint satisfaction, and solution summary. It works with the result
from `solve_trajectory!()`.

# Arguments
- `seq::Sequence`: Original trajectory sequence
- `result`: Result tuple from solve_trajectory! with (variables, objective, info)

# Returns
- Nothing (prints directly to stdout for better formatting)
"""
function report_solution(seq::Sequence, result)
    # Create sequence manager to get ordered information
    sm = SequenceManager(seq)
    
    println()
    println("TRAJECTORY SOLUTION REPORT")
    println("="^50)
    println()
    
    # Optimization Status
    println("OPTIMIZATION STATUS:")
    println("- Converged: $(result.info)")
    
    # Calculate actual counts for consistency
    total_var_objects = length(sm.ordered_vars)
    total_optimization_vars = sum(v.numvars for v in sm.ordered_vars; init=0)
    total_constraint_objects = length(sm.ordered_funcs)
    total_constraint_functions = sum(sm.fun_sizes; init=0)
    
    println("- Variable Objects: $total_var_objects ($total_optimization_vars optimization variables)")
    
    # Use the constraint values from the result
    if haskey(result, :constraints)
        constraint_values = result.constraints
        actual_constraints_returned = length(constraint_values)
    else
        # Fall back to objective if constraints field doesn't exist (backwards compatibility)
        constraint_values = isa(result.objective, AbstractVector) ? result.objective : [result.objective]
        actual_constraints_returned = length(constraint_values)
    end
    
    println("- Constraint Objects: $total_constraint_objects ($total_constraint_functions constraint functions)")
    
    if actual_constraints_returned != total_constraint_functions
        println("- ⚠️  WARNING: Expected $total_constraint_functions constraints but solver returned $actual_constraints_returned")
    end
    
    println()
    
    # Optimization Variables Section
    println("OPTIMIZATION VARIABLES:")
    println("-"^25)
    
    # Parse the flattened variable vector back to individual SolverVariables
    var_idx = 1
    for var in sm.ordered_vars
        calc_desc = get_calc_description(var.calc)
        println("$(var.name) ($calc_desc):")
        
        # Store component values for potential ΔV calculation
        component_values = Float64[]
        
        for i in 1:var.numvars
            opt_value = result.variables[var_idx]
            lower_bound = var.lower_bound[i]
            upper_bound = var.upper_bound[i]
            
            # Store value for ΔV calculation
            push!(component_values, opt_value)
            
            # Format bounds for context
            if abs(lower_bound - upper_bound) < 1e-12
                bounds_str = "(fixed at $(round(upper_bound, digits=6)))"
            else
                bounds_str = "(bounds: [$(round(lower_bound, digits=3)), $(round(upper_bound, digits=3))])"
            end
            
            println("  Component $i: $(round(opt_value, digits=6)) $bounds_str")
            var_idx += 1
        end
        
        # Add Total ΔV for maneuver variables
        if occursin("ManeuverCalc", calc_desc) && occursin("DeltaVVector", calc_desc)
            delta_v_magnitude = norm(component_values)
            println("  Total ΔV: $(round(delta_v_magnitude, digits=6))")
        end
        
        println()
    end
    
    # Constraint Satisfaction Section
    println("CONSTRAINT SATISFACTION:")
    println("-"^25)
    
    # Use constraint values from result
    if haskey(result, :constraints)
        constraint_values = result.constraints
    else
        # Fall back for backwards compatibility
        constraint_values = isa(result.objective, AbstractVector) ? result.objective : [result.objective]
    end
    
    # Parse the constraint vector back to individual constraints
    constraint_idx = 1
    event_constraint_map = get_event_constraint_mapping(sm)
    
    for (event_name, constraints) in event_constraint_map
        if !isempty(constraints)
            println("Event \"$event_name\":")
            
            for (constraint, size) in constraints
                # Extract values for this constraint from the constraint vector
                if constraint_idx + size - 1 <= length(constraint_values)
                    constraint_vals = constraint_values[constraint_idx:constraint_idx + size - 1]
                else
                    # Handle case where we don't have enough constraint values
                    println("  $(get_calc_description(constraint.calc)): No constraint data available")
                    constraint_idx += size
                    continue
                end
                
                # Get constraint description and targets
                calc_desc = get_calc_description(constraint.calc)
                
                if size == 1
                    # Single constraint
                    achieved = constraint_vals[1]
                    target_lower = constraint.lower_bounds[1]
                    target_upper = constraint.upper_bounds[1]
                    
                    # Show achieved value and target bounds
                    if abs(target_lower - target_upper) < 1e-12
                        target_str = "$(round(target_upper, digits=6))"
                    else
                        target_str = "[$(round(target_lower, digits=6)), $(round(target_upper, digits=6))]"
                    end
                    
                    println("  $calc_desc: $(round(achieved, digits=6)) (target: $target_str)")
                else
                    # Multi-component constraint
                    println("  $calc_desc ($size components):")
                    for i in 1:size
                        achieved = constraint_vals[i]
                        target_lower = constraint.lower_bounds[i]
                        target_upper = constraint.upper_bounds[i]
                        
                        if abs(target_lower - target_upper) < 1e-12
                            target_str = "$(round(target_upper, digits=6))"
                        else
                            target_str = "[$(round(target_lower, digits=6)), $(round(target_upper, digits=6))]"
                        end
                        
                        println("    Component $i: $(round(achieved, digits=6)) (target: $target_str)")
                    end
                end
                
                constraint_idx += size
            end
            println()
        end
    end
    
    println("="^50)
    return nothing
end

"""
    get_event_constraint_mapping(sm::SequenceManager)

Create mapping from events to their constraints for reporting.
"""
function get_event_constraint_mapping(sm::SequenceManager)
    event_constraints = Dict{String, Vector{Tuple{Any, Int}}}()
    
    # Initialize all events with empty constraint lists
    for event in sm.sorted_events
        event_constraints[event.name] = Tuple{Any, Int}[]
    end
    
    # Map constraints to their events
    constraint_idx = 1
    for event in sm.sorted_events
        for constraint in event.funcs
            push!(event_constraints[event.name], (constraint, constraint.numvars))
        end
    end
    
    return event_constraints
end




