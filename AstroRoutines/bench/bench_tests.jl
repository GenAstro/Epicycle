# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# In-process benchmark for AstroRoutines tests.
# When include()'d multiple times in the same Julia session, each pass
# reports its own timing. First pass = cold JIT; later passes = warm.

if !@isdefined(__bench_state)
    global __bench_state = Ref(0)
end
__bench_state[] += 1
run_id = __bench_state[]

t0 = time_ns()
using Pkg
t1 = time_ns()
using TestItemRunner
t2 = time_ns()
using AstroRoutines
t3 = time_ns()

t4 = time_ns()
Base.eval(Main, :(@run_package_tests))
t5 = time_ns()

ms(a, b) = round((b - a) / 1e6; digits = 1)
s(a, b)  = round((b - a) / 1e9; digits = 3)

println()
println("======== AstroRoutines timing (in-process run #$run_id) ========")
println(rpad("using Pkg",            22), lpad(string(ms(t0, t1), " ms"), 12))
println(rpad("using TestItemRunner", 22), lpad(string(ms(t1, t2), " ms"), 12))
println(rpad("using AstroRoutines",  22), lpad(string(ms(t2, t3), " ms"), 12))
println(rpad("@run_package_tests",   22), lpad(string(s(t4,  t5),  " s"), 12))
println(rpad("--- run subtotal ---", 22), lpad(string(s(t0,  t5),  " s"), 12))
println("=================================================================")
