# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Reading Plotly's own schema, which PlotlyBase ships.
#
# Two questions are answered here and they come from the same place:
#
#   is_array_ok  — does an array on this attribute mean one value per data point? That decides
#                  whether a vector argument is per-trace or per-point, and it happens on every
#                  call, so it is parsing rather than checking.
#   attr_exists  — is this attribute real? PlotlyBase silently drops one it does not know, so
#                  without this a misspelling is a panel that quietly does nothing.

const _SCHEMA = Ref{Union{Nothing, Dict{Symbol, Any}}}(nothing)

"""Plotly's schema as a plain nested `Dict`, loaded once per session."""
function schema()
    if _SCHEMA[] === nothing
        _SCHEMA[] = getfield(PlotlyBase.get_plotschema(), :fields)
    end
    return _SCHEMA[]
end

"""
Attribute table for one trace type, e.g. `:scatter`. Returns `nothing` for a trace the schema
does not describe rather than throwing, so an unknown trace degrades to "no opinion".
"""
function trace_attributes(trace_type::Symbol)
    traces = get(schema(), :traces, nothing)
    traces === nothing && return nothing
    entry = get(traces, trace_type, nothing)
    entry === nothing && return nothing
    return get(entry, :attributes, nothing)
end

# Walk an attribute key to its schema node. `line_width` is stored nested as line → width, but
# some attributes contain an underscore in their own name (`error_x`, `plot_bgcolor`), so the
# whole key is tried before splitting it.
function _attr_node(attrs, key::Symbol)
    attrs isa AbstractDict || return nothing
    haskey(attrs, key) && return attrs[key]

    parts = Symbol.(split(String(key), '_'))
    length(parts) == 1 && return nothing

    node = attrs
    for p in parts
        node isa AbstractDict || return nothing
        # Nested groups sometimes hold their children under a further :attributes key.
        if !haskey(node, p) && haskey(node, :attributes) && node[:attributes] isa AbstractDict
            node = node[:attributes]
        end
        haskey(node, p) || return nothing
        node = node[p]
    end
    return node
end

"""
    is_array_ok(trace_type, key) -> Bool

Whether Plotly reads an array on this attribute as one value per data point.

`false` for anything the schema does not know, which is the conservative answer: a vector on an
unknown attribute is then treated as per-trace, and per-trace never displaces a meaning Plotly
assigned because there was no known meaning to displace.
"""
function is_array_ok(trace_type::Symbol, key::Symbol)
    attrs = trace_attributes(trace_type)
    attrs === nothing && return false
    node = _attr_node(attrs, key)
    node isa AbstractDict || return false
    return get(node, :arrayOk, false) === true
end

"""
    attr_exists(trace_type, key) -> Bool

Whether the schema knows this attribute. Used to catch typos, not to gate: an attribute the
schema does not know is still passed to PlotlyBase (§3.2), it just earns a warning first.
"""
function attr_exists(trace_type::Symbol, key::Symbol)
    attrs = trace_attributes(trace_type)
    attrs === nothing && return true          # no opinion about an undescribed trace
    return _attr_node(attrs, key) !== nothing
end

"""Whether the schema knows this layout attribute, e.g. `:polar_angularaxis_direction`."""
function layout_attr_exists(key::Symbol)
    lay = get(schema(), :layout, nothing)
    lay isa AbstractDict || return true
    attrs = get(lay, :layoutAttributes, lay)
    return _attr_node(attrs, key) !== nothing
end

# Checking is on by default. It costs a dictionary lookup per attribute, which is nothing beside
# serializing a figure, but a solver callback publishing per iteration can turn it off.
const CHECK_ATTRIBUTES = Ref(true)

"""
    check_attributes!(on::Bool)

Turn attribute-name checking on or off for the session. On by default.
"""
check_attributes!(on::Bool) = (CHECK_ATTRIBUTES[] = on)

# The schema PlotlyBase ships describes plotly.js 2.0 or 2.1, and the dashboard renders with
# 2.35.2, so the schema is the older of the two by about three years. An attribute added to
# plotly.js since then is real, works, and is unknown here — `zorder`, `legend_grouptitlefont`
# and `marker_cornerradius` among them. The check still earns its place, because a misspelling
# is far more common than a newer attribute, but the message has to admit it can be wrong.
# Saying "will be dropped" flatly would be a lie in exactly the case the user is right.
const _SCHEMA_CAVEAT = "If it is a recent Plotly attribute the check may be behind the " *
                       "renderer, in which case it will work; if it is a typo it is dropped " *
                       "in silence and the plot will not show it."

function warn_unknown_attribute(trace_type::Symbol, key::Symbol)
    CHECK_ATTRIBUTES[] || return nothing
    attr_exists(trace_type, key) && return nothing
    @warn "EpicycleIO: `$key` is not an attribute of a Plotly `$trace_type` trace in the " *
          "schema shipped with PlotlyBase. $_SCHEMA_CAVEAT Check it against " *
          "https://plotly.com/javascript/reference/$trace_type/"
    return nothing
end

function warn_unknown_layout(key::Symbol)
    CHECK_ATTRIBUTES[] || return nothing
    layout_attr_exists(key) && return nothing
    @warn "EpicycleIO: `$key` is not a Plotly layout attribute in the schema shipped with " *
          "PlotlyBase. $_SCHEMA_CAVEAT Check it against " *
          "https://plotly.com/javascript/reference/layout/"
    return nothing
end
