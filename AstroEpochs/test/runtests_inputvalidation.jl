# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

using Test
using AstroEpochs

@testset "Time Input Validation" begin

# Scales and formats are tags. A symbol or a string fails with a message that shows the tag form.
for bad in (() -> Time(2451545.0, 0.0, :tt, JD()), () -> Time(2451545.0, :tt, JD()),
            () -> Time("2024-01-01T00:00:00", :tt, ISOT()), () -> Time(2451545.0, 0.0, "tt", JD()))
    err = try bad(); nothing catch e; e end
    @test err isa ArgumentError
    @test occursin("scale must be a time scale tag such as TT()", sprint(showerror, err))
end
for bad in (() -> Time(2451545.0, 0.0, TT(), :jd), () -> Time(2451545.0, TT(), :jd))
    err = try bad(); nothing catch e; e end
    @test err isa ArgumentError
    @test occursin("format must be a time format tag such as JD()", sprint(showerror, err))
end
# Tests error in input format coupling 
@test begin
    err = try
        AstroEpochs._validate_inputcoupling(1.0, :foo)
        nothing
    catch e
        e
    end
    err isa ArgumentError && occursin("Time: unsupported format :foo in this constructor", sprint(showerror, err))
end

# Test error is thrown when Real is provided for ISOT format. 
@test begin
    err = try
        Time(59000.0, TT(), ISOT())  # Real with :isot triggers Unsupported time format
        nothing
    catch e
        e
    end
    err isa ArgumentError && occursin("Time: format ISOT() requires AbstractString input", sprint(showerror, err))
end



# get_conversion_path (scales/graph.jl): the reverse MULTI_HOPS path
@test AstroEpochs.get_conversion_path(:tcb, :tai) == [:tcb, :tdb, :tt, :tai]


# Test missing path is trapped and thrown
@test begin
    err = try
        AstroEpochs.get_conversion_path(:tt, :notascale)
        nothing
    catch e
        e
    end
    err !== nothing && occursin("No known time scale conversion path", sprint(showerror, err))
end

# Invalid input types path
@test_throws ArgumentError Time(2451545.0, 0.25, 42, 99)
   
end
nothing
