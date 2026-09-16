# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# Requirements-traceability pilot (TestingStandards.md §1.1 + "Requirements Traceability").
# Test files follow test_<axis>_<feature>.jl — test_correctness_anom.jl, test_differentiability_anom.jl,
# test_robustness_anom.jl — and each @testitem is tagged with its requirement ID(s) and dimension.
# @run_package_tests discovers every @testitem in the package automatically (no includes needed).
using TestItemRunner
@run_package_tests
