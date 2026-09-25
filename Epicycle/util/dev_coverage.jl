# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Fast development coverage for individual modules
# Usage: include("dev_coverage.jl"); dev_coverage("AstroManeuvers")

using Coverage, Printf, Test

"""
    dev_coverage(module_name::String; run_tests::Bool=true)

Fast coverage analysis for development - runs in current Julia session.
Much faster than full recompilation approach.
"""
function dev_coverage(module_name::String; run_tests::Bool=true)
    println("🚀 Fast coverage analysis for $module_name...")
    
    # Paths
    srcdir = joinpath(module_name, "src")
    testjl = joinpath(module_name, "test", "runtests.jl")
    
    if !isdir(srcdir)
        error("Source directory not found: $srcdir")
    end
    
    # Clean existing .cov files
    Coverage.clean_folder(srcdir)
    if isdir(dirname(testjl))
        Coverage.clean_folder(dirname(testjl))
    end
    
    if run_tests && isfile(testjl)
        println("📋 Running tests with coverage tracking...")
        
        # Enable coverage for this session
        # Note: This won't track as precisely as --code-coverage=user but it's much faster
        try
            # Run the tests (assumes module is already loaded)
            include(testjl)
            println("✅ Tests completed")
        catch e
            println("❌ Test error: ", e)
            return nothing
        end
    end
    
    # Process coverage
    println("📊 Processing coverage data...")
    cov = Coverage.process_folder(srcdir)
    
    if isempty(cov)
        println("⚠️  No coverage data found. Try running with:")
        println("   julia --code-coverage=user")
        println("   Then: include(\"dev_coverage.jl\"); dev_coverage(\"$module_name\", run_tests=false)")
        return nothing
    end
    
    # Quick summary
    total_cov = 0; total_lines = 0
    println("\n📈 Coverage Summary:")
    println("="^50)
    
    for f in sort(cov; by = x -> x.filename)
        lines = f.coverage
        n_total = count(!isnothing, lines)
        n_covered = count(x -> x !== nothing && x > 0, lines)
        pct = n_total == 0 ? 100.0 : 100 * n_covered / n_total
        
        status = pct >= 90 ? "🟢" : pct >= 70 ? "🟡" : "🔴"
        filename = basename(f.filename)
        
        println(@sprintf("%s %-30s %5.1f%% (%3d/%3d)", status, filename, pct, n_covered, n_total))
        total_cov += n_covered; total_lines += n_total
    end
    
    println("="^50)
    tpct = total_lines == 0 ? 100.0 : 100 * total_cov / total_lines
    status = tpct >= 90 ? "🟢" : tpct >= 70 ? "🟡" : "🔴"
    println(@sprintf("%s %-30s %5.1f%% (%3d/%3d)", status, "TOTAL", tpct, total_cov, total_lines))
    
    return cov
end

"""
    show_uncovered(module_name::String; context::Int=2)

Show uncovered lines with context for quick gap analysis.
"""
function show_uncovered(module_name::String; context::Int=2)
    srcdir = joinpath(module_name, "src")
    cov = Coverage.process_folder(srcdir)
    
    if isempty(cov)
        println("⚠️  No coverage data found for $module_name")
        return
    end
    
    total_uncovered = 0
    
    for f in sort(cov; by = x -> x.filename)
        lines = f.coverage
        zeros = Int[]
        
        for (i, v) in enumerate(lines)
            v === 0 && push!(zeros, i)
        end
        
        if !isempty(zeros)
            total_uncovered += length(zeros)
            println("\n🔍 ", basename(f.filename))
            println("   Uncovered lines: ", join(string.(zeros), ", "))
            
            if context > 0
                src = readlines(f.filename)
                for line_num in zeros[1:min(5, length(zeros))]  # Show first 5 uncovered lines
                    sctx = max(1, line_num - context)
                    ectx = min(length(src), line_num + context)
                    
                    println("   Lines $(sctx)-$(ectx):")
                    for i in sctx:ectx
                        mark = i == line_num ? "❌" : "  "
                        println(@sprintf("   %s %3d | %s", mark, i, src[i]))
                    end
                    println()
                    
                    length(zeros) > 5 && line_num == zeros[5] && println("   ... ($(length(zeros)-5) more uncovered lines)")
                end
            end
        end
    end
    
    if total_uncovered == 0
        println("🎉 No uncovered lines found!")
    else
        println("📋 Total uncovered lines: $total_uncovered")
    end
end

"""
    quick_test_coverage(module_name::String, testset_name::String)

Run a specific testset and show coverage for just that.
"""
function quick_test_coverage(module_name::String, testset_name::String)
    println("🏃 Quick test: $testset_name in $module_name")
    
    # Clean coverage
    srcdir = joinpath(module_name, "src")
    Coverage.clean_folder(srcdir)
    
    # Run specific testset (this is a bit hacky but works for development)
    testjl = joinpath(module_name, "test", "runtests.jl")
    if isfile(testjl)
        # You'd need to modify this based on your test structure
        # For now, just run all tests
        include(testjl)
    end
    
    return dev_coverage(module_name; run_tests=false)
end

println("📚 Dev coverage tools loaded!")
println("Usage:")
println("  dev_coverage(\"AstroManeuvers\")              # Full coverage")
println("  show_uncovered(\"AstroManeuvers\")            # Show gaps")
println("  dev_coverage(\"AstroManeuvers\", run_tests=false) # Just analyze existing .cov files")

nothing
