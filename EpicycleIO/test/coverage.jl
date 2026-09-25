# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Line coverage for EpicycleIO.
#
#   julia --project=<environment> EpicycleIO/test/coverage.jl --run
#
# Or in two steps, instrumenting this package only:
#
#   julia --code-coverage=@<path to EpicycleIO> --project=… EpicycleIO/test/runtests.jl
#   julia --project=… EpicycleIO/test/coverage.jl
#
# Use `@<dir>` rather than `user`. `user` instruments every package the tests touch, which wrote
# some seven thousand .cov files into the shared depot alongside PlotlyBase and the Astro
# packages. They are harmless and invisible until someone greps the depot, and they are not ours
# to leave there.
#
# Coverage says which lines ran, not whether the behaviour is right. Everything the browser
# does is invisible to it, so a high number here is not evidence the plots are correct.

using Coverage
using Printf

const PKG = dirname(@__DIR__)
const SRC = joinpath(PKG, "src")

"Delete .cov files under a directory, including the gl subdirectory."
function clean_cov(dir)
    for (root, _, files) in walkdir(dir), f in files
        endswith(f, ".cov") && rm(joinpath(root, f); force = true)
    end
    return nothing
end

function run_tests_with_coverage()
    clean_cov(PKG)
    project = Base.active_project()
    run(`$(Base.julia_cmd()) --startup-file=no --code-coverage=@$PKG --project=$project
         $(joinpath(@__DIR__, "runtests.jl"))`)
    return nothing
end

function report()
    cov = process_folder(SRC)
    isempty(cov) && error("no coverage files in $SRC — run the tests with --code-coverage=user first")

    rows = NTuple{4, Any}[]
    covered_total = tested_total = 0
    for f in cov
        c, t = get_summary(f)
        covered_total += c
        tested_total  += t
        push!(rows, (basename(f.filename), c, t, t == 0 ? 100.0 : 100c / t))
    end
    sort!(rows; by = r -> r[4])

    println("\nEpicycleIO line coverage\n")
    for (name, c, t, pct) in rows
        @printf("  %-22s %6.1f%%   %4d/%-4d\n", name, pct, c, t)
    end
    @printf("  %-22s %6.1f%%   %4d/%-4d\n", "TOTAL",
            100covered_total / tested_total, covered_total, tested_total)

    println("\nUncovered lines\n")
    for f in cov
        misses = [i for (i, n) in enumerate(f.coverage) if n !== nothing && n == 0]
        isempty(misses) && continue
        println("  ", basename(f.filename), ": ", join(misses, ", "))
    end
    println()
    return 100covered_total / tested_total
end

if abspath(PROGRAM_FILE) == @__FILE__
    "--run" in ARGS && run_tests_with_coverage()
    report()
end
