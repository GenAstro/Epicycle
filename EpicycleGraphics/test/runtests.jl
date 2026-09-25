# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

using EpicycleGraphics
using Test

@testset "EpicycleGraphics" begin
    include("runtests_graphics.jl")

    # Visual regression compares rendered frames against the images in test/reference/. It is
    # excluded from CI because the comparison is sensitive to the GPU and driver that drew the
    # frame, so a fresh runner fails it for reasons that are not defects. Run it by hand before
    # a release, on the machine whose references are checked in.
    # include("visual_regression.jl")
end
