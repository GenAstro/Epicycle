using Test
using AstroEpochs

@testset "Time Input Validation" begin

# Invalid scale 
#@test begin
#    err = try
#        Time(2451545.0, 0.0, "tt", :jd)
#        nothing
#    catch e
#        e
#    end
#    err isa ArgumentError && occursin("invalid time input types. scale must be Symbol; got String.", sprint(showerror, err))
#end
#=
# Invalid scale (scale passed as String)
err_scale = let e = try
        Time(2451545.0, 0.0, "tt", :jd)
        nothing
    catch err
        err
    end
    e
end
@test err_scale isa ArgumentError
@test sprint(showerror, err_scale) == "ArgumentError: Time: invalid time input types. scale must be Symbol; got String."

# Invalid format (format passed as String)
err_format = let e = try
        Time(2451545.0, 0.0, :tt, "jd")
        nothing
    catch err
        err
    end
    e
end
@test err_format isa ArgumentError
@test sprint(showerror, err_format) == "ArgumentError: Time: invalid time input types. format must be Symbol; got String."
=#
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
        Time(59000.0, :tt, :isot)  # Real with :isot triggers Unsupported time format
        nothing
    catch e
        e
    end
    err isa ArgumentError && occursin("Time: format ISOT() requires AbstractString input", sprint(showerror, err))
end



# Lines 439–440: reverse MULTI_HOPS path
#@test begin
#    p = AstroEpochs.get_conversion_path(:tcb, :tai)  
#    p == [:tcb, :tdb, :tt, :tai]
#end


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

# Invalid input types path (non-Symbol scale/format)
@test_throws ArgumentError Time(2451545.0, 0.25, 42, 99)
   
end
nothing