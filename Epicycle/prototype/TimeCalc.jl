using AstroModels 
using AstroEpochs  
import AstroCallbacks: AbstractCalcVariable, AbstractCalc, set_calc!
import AstroEpochs: AbstractTimeScale, AbstractTimeFormat

# New abstract type for time variables  
abstract type AbstractTimeVar <: AbstractCalcVariable end

# Main TimeCalc struct - reuse existing AstroEpochs types
struct TimeCalc{S<:Spacecraft, TS<:AbstractTimeScale, TF<:AbstractTimeFormat} <: AbstractCalc
    sc::S
    scale::TS
    format::TF
end

# Traits
calc_is_settable(::TimeCalc) = true
calc_numvars(::TimeCalc) = 1

# Evaluation
function get_calc(c::TimeCalc)
    time = c.sc.time
    # Convert to requested scale and format using AstroEpochs functionality
    converted_time = Time(time.jd1, time.jd2, c.scale, c.format)  # or whatever the conversion API is
    return converted_time.mjd  # Return the numeric value
end

# Setting  
function _set_calc_type!(c::TimeCalc, newval::Real)
    # Create new time with the specified scale/format
    new_time = Time(newval, c.scale, c.format)
    # Convert back to spacecraft's native scale/format and assign
    c.sc.time = Time(new_time, c.sc.time.scale, c.sc.time.format)
    return c.sc.time
end

# Usage examples using existing AstroEpochs tags:
sat = Spacecraft()
tt_mjd   = TimeCalc(sat, TT(), MJD())
utc_jd   = TimeCalc(sat, UTC(), JD())
tai_isot = TimeCalc(sat, TAI(), ISOT())

get_calc(tt_mjd)  # Get time in TT MJD
set_calc!(tt_mjd, 59000.0)  # Set time in TT