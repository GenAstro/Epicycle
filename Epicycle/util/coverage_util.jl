#= 

# ===== To run from monorepo root (c:\Users\steve\Dev\Epicycle). 
module_name = "AstroBase"      # Done
module_name = "AstroStates"
module_name = "AstroEpochs"    # Done
module_name = "AstroUniverse"
module_name = "AstroCoords"
module_name = "AstroModels"
module_name = "AstroMan"
module_name = "AstroFun"

module_name = "AstroProp"
module_name = "AstroSolve"

# Usage:
module_name = "AstroMan"
include("Epicycle/util/coverage_util.jl")                
run_module_coverage(module_name; clean=true)   
report_uncovered(module_name; context=1)

report_uncovered(module_name; files=["AstroStates.jl"], context = 1)

report_uncovered(module_name) # list ranges
report_uncovered(module_name; context=1) # show 1 line of context
report_uncovered(["AstroBase/src","AstroStates/src"])# multiple packages
report_uncovered(module_name; exclude=["AstroStates/src/AstroStates.jl"]) # skip loader

=#

using Coverage, Printf

# Collect src coverage

@info "Collecting coverage for $module_name"

# Simple per-file + total report
function simple_report(cov)
    total_cov = 0; total_lines = 0
    for f in sort(cov; by = x -> x.filename)
        lines = f.coverage
        n_total   = count(!isnothing, lines)
        n_covered = count(x -> x !== nothing && x > 0, lines)
        pct = n_total == 0 ? 100.0 : 100 * n_covered / n_total
        println(rpad(f.filename, 70), @sprintf("%5.1f%% (%d/%d)", pct, n_covered, n_total))
        total_cov += n_covered; total_lines += n_total
    end
    println("---")
    tpct = total_lines == 0 ? 100.0 : 100 * total_cov / total_lines
    println(@sprintf("Total coverage: %5.1f%% (%d/%d)", tpct, total_cov, total_lines))
end

"""
function run_module_coverage(module_name::AbstractString; project::AbstractString=pwd(), run_tests::Bool=true, clean::Bool=true)

Run `module_name/test/runtests.jl` under coverage in a fresh Julia process, then
summarize coverage for `module_name/src` and write `cov/<module>_lcov.info`.

If you already started Julia with --code-coverage=user and ran tests, you can call
this with `run_tests=false` to only summarize.
"""
function run_module_coverage(module_name::AbstractString; project::AbstractString=pwd(), run_tests::Bool=true, clean::Bool=true)
    root   = project
    srcdir = joinpath(root, module_name, "src")
    testjl = joinpath(root, module_name, "test", "runtests.jl")
    modsrc = joinpath(root, module_name, "src", "$(module_name).jl")
    
    # Create cov directory if it doesn't exist
    covdir = joinpath(root, "cov")
    mkpath(covdir)
    out = joinpath(covdir, "$(module_name)_lcov.info")

    @info "Coverage paths" root srcdir testjl modsrc out

    if clean
        rm_cov!([srcdir, dirname(testjl)])
    end

    if run_tests && isfile(testjl)
        # Ensure we include the exact src file path we will later process for coverage.
        script = """
            cd($(repr(root)))
            # Generate .cov for $(modsrc) without defining Main.$(module_name)
            module __COV__$(module_name)
                include($(repr(modsrc)))
            end
            # Now run the package tests (they may include the module again into Main)
            try
                include($(repr(testjl)))
            catch e
                showerror(stderr, e, catch_backtrace()); println(stderr)
                rethrow()
            end
        """
        mktemp() do path, io
            write(io, script); close(io)
            # Use user or all; if user misses, try :all
            run(`$(Base.julia_cmd()) --project=$(root) --code-coverage=user -O0 --inline=no $(path)`)
        end
    end

    @info "Processing coverage" srcdir
    if !isdir(srcdir)
        error("Source directory does not exist: $srcdir")
    end
    cov = Coverage.process_folder(srcdir)
    simple_report(cov)
    Coverage.LCOV.writefile(out, cov)
    @info "LCOV written" out
    return cov
end

# Collapse a sorted list of line numbers into ranges like ["10-14","17","21-22"]
function _collapse_ranges(lines::Vector{Int})
    isempty(lines) && return String[]
    ranges = String[]
    s = lines[1]; e = s
    for i in Iterators.drop(lines, 1)
        if i == e + 1
            e = i
        else
            push!(ranges, s == e ? string(s) : string(s, "-", e))
            s = e = i
        end
    end
    push!(ranges, s == e ? string(s) : string(s, "-", e))
    return ranges
end

# Print uncovered (zero-hit) executable lines per file, with optional context
function report_uncovered(paths::Vector{<:AbstractString};
                          files::Vector{<:AbstractString} = String[],
                          exclude::Vector{<:AbstractString} = String[],
                          context::Int = 0)
    cov = reduce(vcat, map(Coverage.process_folder, paths))
    # Optional: restrict to specific files (match by basename)
    if !isempty(files)
        wanted = Set(files)
        cov = filter(f -> basename(f.filename) in wanted, cov)
    end
    total_exec = 0
    total_uncovered = 0

    for f in sort(cov; by = x -> x.filename)
        any(occursin(p, f.filename) for p in exclude) && continue

        lines = f.coverage
        # Executable lines are those where coverage[i] !== nothing
        exec_idxs = findall(!isnothing, lines)
        zeros = Int[]
        for i in exec_idxs
            v = lines[i]
            # v is Int or nothing; we already filtered nothing above
            v === 0 && push!(zeros, i)
        end

        total_exec += length(exec_idxs)
        total_uncovered += length(zeros)

        if !isempty(zeros)
            ranges = _collapse_ranges(sort(zeros))
            println("\n", f.filename)
            println("  Uncovered ranges: ", join(ranges, ", "))
            if context > 0
                src = readlines(f.filename)
                for r in ranges
                    s, e = occursin("-", r) ? parse.(Int, split(r, "-")) : (parse(Int, r), parse(Int, r))
                    sctx = max(1, s - context)
                    ectx = min(length(src), e + context)
                    for i in sctx:ectx
                        mark = (i >= s && i <= e) ? ">" : " "
                        covv = lines[i]
                        hit = covv === nothing ? "·" : string(covv)
                        @printf("  %s %6d | %3s | %s\n", mark, i, hit, src[i])
                    end
                    println()
                end
            end
        end
    end

    covered = total_exec - total_uncovered
    pct = total_exec == 0 ? 100.0 : 100 * covered / total_exec
    println("\n---")
    @printf("Total executable lines: %d, covered: %d, uncovered: %d (%.1f%%)\n",
            total_exec, covered, total_uncovered, pct)
    return nothing
end

function rm_cov!(roots::Vector{<:AbstractString})
    n = 0
    for r in roots
        isdir(r) || continue
        for (dir, _, files) in walkdir(r)
            for f in files
                endswith(f, ".cov") || continue
                rm(joinpath(dir, f); force=true)
                n += 1
            end
        end
    end
    @info "Removed $n .cov files" roots
    return n
end

# Convenience for a single module directory, e.g. "AstroBase" or "AstroStates"
report_uncovered(mod::AbstractString; kwargs...) =
    report_uncovered([joinpath(mod, "src")]; kwargs...)

nothing