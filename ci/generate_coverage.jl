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

# DEBUG: Show current working directory and environment
println("🐛 DEBUG: Current working directory: $(pwd())")
println("🐛 DEBUG: JULIA_CODE_COVERAGE env var: $(get(ENV, "JULIA_CODE_COVERAGE", "NOT SET"))")
println("🐛 DEBUG: Contents of current directory:")
for item in readdir(".")
    if isdir(item)
        println("  📁 $item/")
    else
        println("  📄 $item")
    end
end

# Look for .cov files directly in package src directories
all_coverage = Coverage.FileCoverage[]

packages = [
    "AstroBase", "AstroStates", "AstroEpochs", "AstroUniverse",
    "AstroCoords", "AstroModels", "AstroMan", "AstroFun", 
    "AstroProp", "AstroSolve", "Epicycle"
]

# DEBUG: Search for .cov files everywhere first
println("🐛 DEBUG: Searching for ALL .cov files in entire directory tree...")
function find_cov_files(dir=".")
    cov_files = []
    try
        for (root, dirs, files) in walkdir(dir)
            for file in files
                if endswith(file, ".cov")
                    push!(cov_files, joinpath(root, file))
                end
            end
        end
    catch e
        println("  ⚠️  Error walking directory $dir: $e")
    end
    return cov_files
end

all_cov_files = find_cov_files()
println("🐛 DEBUG: Found $(length(all_cov_files)) .cov files total:")
for cov_file in all_cov_files[1:min(20, length(all_cov_files))]  # Limit output
    println("  📊 $cov_file")
end
if length(all_cov_files) > 20
    println("  ... and $(length(all_cov_files) - 20) more")
end

for pkg in packages
    src_path = "$pkg/src"
    println("🔍 Processing coverage for $pkg...")
    println("🐛 DEBUG: Checking directory: $src_path")
    
    if isdir(src_path)
        println("  ✅ Directory exists")
        
        # DEBUG: Show contents of src directory
        println("🐛 DEBUG: Contents of $src_path:")
        try
            for item in readdir(src_path)
                full_path = joinpath(src_path, item)
                if isdir(full_path)
                    println("    📁 $item/")
                else
                    println("    📄 $item")
                    if endswith(item, ".cov")
                        println("      🎯 FOUND .cov file!")
                    end
                end
            end
        catch e
            println("    ⚠️  Error reading directory: $e")
        end
        
        # Use the standard Coverage.process_folder function
        try
            println("🐛 DEBUG: Calling Coverage.process_folder(\"$src_path\")")
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
        println("  ❌ Directory $src_path does not exist")
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