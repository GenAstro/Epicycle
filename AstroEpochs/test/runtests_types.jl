using Test
using AstroEpochs

@testset "Time scale tags <-> symbols" begin
    scale_pairs = [
        (TT(),  :tt),
        (TDB(), :tdb),
        (UTC(), :utc),
        (TCB(), :tcb),
        (TCG(), :tcg),
        (TAI(), :tai),
    ]

    for (tag, sym) in scale_pairs
        # tag -> symbol
        @test AstroEpochs._scale_symbol(tag) == sym
        @test Core.invoke(AstroEpochs._scale_symbol, Tuple{typeof(tag)}, tag) == sym

        # symbol -> tag
        t2 = AstroEpochs._scale_tag(Val(sym))
        @test t2 isa typeof(tag)
        @test Core.invoke(AstroEpochs._scale_tag, Tuple{Val{sym}}, Val(sym)) isa typeof(tag)

        # round-trip consistency
        @test AstroEpochs._scale_symbol(t2) == sym
    end

end

@testset "Time format tags <-> symbols" begin
    # Note: PrecJD maps to :jd
    format_pairs = [
        (JD(),     :jd),
        (MJD(),    :mjd),
        (ISOT(),   :isot),
    ]

    for (tag, sym) in format_pairs
        # tag -> symbol
        @test AstroEpochs._format_symbol(tag) == sym
        @test Core.invoke(AstroEpochs._format_symbol, Tuple{typeof(tag)}, tag) == sym

        # symbol -> tag
        t2 = AstroEpochs._format_tag(Val(sym))
        @test t2 isa typeof(tag)
        @test Core.invoke(AstroEpochs._format_tag, Tuple{Val{sym}}, Val(sym)) isa typeof(tag)

        # round-trip consistency
        @test AstroEpochs._format_symbol(t2) == sym
    end

end