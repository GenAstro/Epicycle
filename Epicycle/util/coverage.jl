# Convenience script for running coverage analysis on monorepo packages
# Usage: julia coverage.jl AstroManeuvers

include("Epicycle/util/coverage_util.jl")

# Get module name from command line or use default
module_name = length(ARGS) > 0 ? ARGS[1] : "AstroManeuvers"

println("Running coverage analysis for $module_name...")
run_module_coverage(module_name; clean=true)   
report_uncovered(module_name; context=1)