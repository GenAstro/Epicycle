# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT


using Test
using AstroEpochs

include("test_construction.jl")
include("test_arithmetic.jl")
include("test_formatconversions.jl")
include("test_scaleconversions.jl")
include("test_io.jl")
include("runtests_types.jl")
include("runtests_inputvalidation.jl")
include("test_differentiation.jl")
include("test_correctness_date_parts.jl")
include("test_correctness_leap_seconds.jl")



