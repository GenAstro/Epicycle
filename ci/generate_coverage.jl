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

# Load Glob for file pattern matching
try
    using Glob
catch
    println("📥 Installing Glob.jl...")
    Pkg.add("Glob")
    using Glob
end

println("🔍 Processing coverage files...")

# Look for .cov files directly in package src directories
all_coverage = Coverage.FileCoverage[]

packages = [
    "AstroBase", "AstroStates", "AstroEpochs", "AstroUniverse",
    "AstroCoords", "AstroModels", "AstroMan", "AstroFun", 
    "AstroProp", "AstroSolve", "Epicycle"
]

for pkg in packages
    src_path = "$pkg/src"
    if isdir(src_path)
        println("🔍 Processing coverage for $pkg...")
        
        # Use the standard Coverage.process_folder function
        # This should find .cov files automatically
        try
            pkg_coverage = Coverage.process_folder(src_path)
            if !isempty(pkg_coverage)
                append!(all_coverage, pkg_coverage)
                println("  ✅ Processed $(length(pkg_coverage)) coverage files for $pkg")
            else
                println("  ⚠️  No coverage processed for $pkg")
            end
        catch e
            println("  ❌ Error processing coverage for $pkg: $e")
        end
    else
        println("  → Directory $src_path does not exist")
    end
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