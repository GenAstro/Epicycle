```@meta
CurrentModule = AstroFun
```

# Developer API Reference

This page documents internal implementation details for developers contributing to AstroFun.

!!! warning "Internal API"
    These functions are implementation details and may change without notice. 
    Use the public API documented in the main reference guide for stable interfaces.

## Internal Functions

```@docs
_evaluate
_set!
_extract_mu
_state_for_calc
_set_calc_type!
_infer_numvars
_subjects_from_calc
```

## Type System Functions

```@docs
calc_numvars
calc_is_settable
calc_input_statetag
to_concrete_state
convert_orbitcalc_state
```

## Abstract Types

```@docs
AbstractCalc
AbstractCalcVariable
AbstractOrbitVar
AbstractManeuverVar
AbstractBodyVar
```