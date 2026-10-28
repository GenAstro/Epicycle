# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# Loaded only when ForwardDiff is: gives interpolation the value of a dual number, so that range
# checks and window choice ignore the partials. ForwardDiff 1.x orders duals with equal values by
# their partials, so Dual(6.0, 1.0) > 6.0, which would put a dual at the last tabulated time
# outside the table.
module AstroRoutinesForwardDiffExt

using AstroRoutines, ForwardDiff

AstroRoutines._primal(x::ForwardDiff.Dual) = AstroRoutines._primal(ForwardDiff.value(x))

end
