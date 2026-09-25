# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Runs every example headlessly and reports which ones fail.
#
# The examples are the interface's specification, so "do they still run" is a real check rather
# than a nicety. It does not look at the pictures — that stays a visual check, honestly named.
#
#   julia --project=<environment> EpicycleIO/test/run_examples.jl

using EpicycleIO

EpicycleIO.auto_open!(false)

const EXAMPLES = joinpath(dirname(@__DIR__), "examples")

function main()
    files = sort(filter(f -> startswith(f, "Ex_") && endswith(f, ".jl"), readdir(EXAMPLES)))
    isempty(files) && error("no examples found in $EXAMPLES")

    ok, failed = 0, String[]
    for f in files
        EpicycleIO.reset_panels!()
        try
            Main.include(joinpath(EXAMPLES, f))
            ok += 1
            println("  ok    ", f)
        catch e
            push!(failed, f)
            println("  FAIL  ", f)
            println("          ", first(sprint(showerror, e), 300))
        end
    end

    println()
    if isempty(failed)
        println("all $ok examples ran")
    else
        println("$ok ran, $(length(failed)) failed: ", join(failed, ", "))
    end
    EpicycleIO.close_dashboard()
    return isempty(failed)
end

success = main()
exit(success ? 0 : 1)
