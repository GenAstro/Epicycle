# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

using AstroFrames, Test

# Documentation contract: every plain Julia docstring block must run after loading
# AstroFrames, with the remaining imports and inputs supplied by the example itself.
@testset "Plain docstring examples" begin
    source = joinpath(pkgdir(AstroFrames), "src")
    for file in sort(filter(f -> endswith(f, ".jl"), readdir(source; join = true)))
        for (i, block) in enumerate(eachmatch(r"```julia\r?\n(.*?)```"s, read(file, String)))
            @testset "$(basename(file)) example $(i)" begin
                example = Module(gensym(:AstroFramesExample))
                Core.eval(example, :(using AstroFrames))
                @test begin
                    Base.include_string(example, block.captures[1], file)
                    true
                end
            end
        end
    end
end
