#!/usr/bin/env julia
"""
CI Generate Coverage Script

Generates and processes code coverage information.
This runs after tests have completed successfully.
"""

println("📈 Generating coverage...")

using Pkg
Pkg.activate(".")

# Add coverage packages if not already present
try
    using Coverage
catch
    println("📥 Installing Coverage.jl...")
    Pkg.add("Coverage")
    using Coverage
end

# also load CoverageTools and LCOV helpers
try
    using CoverageTools
catch
    println("📥 Installing CoverageTools.jl...")
    Pkg.add("CoverageTools")
    using CoverageTools
end

println("🔍 Processing coverage files...")

# Look for coverage files in all package directories
# Process each package and collect FileCoverage objects
all_coverage = Coverage.FileCoverage[]

packages = [
    "AstroBase", "AstroStates", "AstroEpochs", "AstroUniverse",
    "AstroCoords", "AstroModels", "AstroMan", "AstroFun", 
    "AstroProp", "AstroSolve", "Epicycle"
]

for pkg in packages
    pkg_coverage = Coverage.process_folder("$pkg/src")
    append!(all_coverage, pkg_coverage)
    println("  → Found $(length(pkg_coverage)) coverage files in $pkg/src")
end

if !isempty(all_coverage)
    println("📊 Coverage summary:")
    covered_lines, total_lines = Coverage.get_summary(all_coverage)
    percentage = round(covered_lines / total_lines * 100, digits=2)
    println("  → Coverage: $covered_lines/$total_lines lines ($percentage%)")
    
    # Generate LCOV file for upload
    println("📄 Generating LCOV file...")
    Coverage.LCOV.writefile("coverage.info", all_coverage)
    println("  → Coverage data written to coverage.info")
    
    println("✅ Coverage generation complete!")
else
    println("⚠️  No coverage files found")
end