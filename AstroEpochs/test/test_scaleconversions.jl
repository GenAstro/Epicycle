

using Test
using AstroEpochs

function isapproxrel(a::Float64, b::Float64; rtol::Float64=1e-14, atol::Float64=1e-14)::Bool
    abs(a - b) ≤ atol || abs(a - b) / max(abs(a), abs(b), eps()) ≤ rtol
end

@testset "Scale Conversion — (AstroPy Truth)" begin
    
    # Define the time
    t = Time("2017-09-21T12:23:12", TAI(), ISOT())
    @test isapproxrel(t.jd1, 2458018.0)
    @test isapproxrel(t.jd2, 0.016111111111111076)

    # Symbol based construction from UTC
    t_utc = t.utc
    @test isapproxrel(t_utc.jd1, 2458018.0)
    @test isapproxrel(t_utc.jd2, 0.015682870370370305; atol = 1e-9, rtol = 1e-9)  # Due to using D2000 Ref in conversion

    # TT
    t_tt = t.tt
    @test isapproxrel(t_tt.jd1, 2458018.0)
    @test isapproxrel(t_tt.jd2, 0.016483611111111077)

    # TDB
    t_tdb = t.tdb
    @test isapproxrel(t_tdb.jd1, 2458018.0)
    @test isapproxrel(t_tdb.jd2, 0.016483592271460668; atol = 1e-9, rtol = 1e-9) # Due to lower precision tdb conversion and different model

    # TCB
    t_tcb = t.tcb
    @test isapproxrel(t_tcb.jd1, 2458018.0)
    @test isapproxrel(t_tcb.jd2, 0.016714209840637442; atol = 1e-9, rtol = 1e-9) # Due to lower precision tdb conversion and different model

    # TCG
    t_tcg = t.tcg
    @test isapproxrel(t_tcg.jd1, 2458018.0)
    @test isapproxrel(t_tcg.jd2, 0.01649397689602741)
end

# TODO: Round trip tests probably 

# Identity test for conversion paths 
@test begin
    p = AstroEpochs.get_conversion_path(:tt, :tt)
    p == [:tt]
end
@test begin
    p = AstroEpochs.get_conversion_path(:tdb, :tdb)
    p == [:tdb]
end
@test begin
    p = AstroEpochs.get_conversion_path(:utc, :utc)
    p == [:utc]
end
@test begin
    p = AstroEpochs.get_conversion_path(:tai, :tai)
    p == [:tai]
end
@test begin
    p = AstroEpochs.get_conversion_path(:tcg, :tcg)
    p == [:tcg]
end
@test begin
    p = AstroEpochs.get_conversion_path(:tcb, :tcb)
    p == [:tcb]
end

# Lines 443–444: reverse single-edge path (temporarily remove forward edge)
@test begin
    key_fwd = (:tt, :tai)
    key_rev = (:tai, :tt)
    fwd = haskey(AstroEpochs.OFFSET_TABLE, key_fwd) ? AstroEpochs.OFFSET_TABLE[key_fwd] : nothing
    try
        if fwd !== nothing
            delete!(AstroEpochs.OFFSET_TABLE, key_fwd)
        end
        p = AstroEpochs.get_conversion_path(:tt, :tai)
        p == [:tt, :tai]
    finally
        if fwd !== nothing && !haskey(AstroEpochs.OFFSET_TABLE, key_fwd)
            AstroEpochs.OFFSET_TABLE[key_fwd] = fwd
        end
    end
end

# Canonical JD truth data (AstroPy) and named Time instances
const JD_TRUTH = 2458018.0
const JD2_TAI = 0.016111111111111076
const JD2_UTC = 0.015682870370370305
const JD2_TT  = 0.016483611111111077
const JD2_TDB = 0.016483592271460668
const JD2_TCB = 0.016714209840637442
const JD2_TCG = 0.01649397689602741

t_tai = Time(JD_TRUTH, JD2_TAI, TAI(), JD())
t_utc = Time(JD_TRUTH, JD2_UTC, UTC(), JD())
t_tt  = Time(JD_TRUTH, JD2_TT,  TT(),  JD())
t_tdb = Time(JD_TRUTH, JD2_TDB, TDB(), JD())
t_tcb = Time(JD_TRUTH, JD2_TCB, TCB(), JD())
t_tcg = Time(JD_TRUTH, JD2_TCG, TCG(), JD())

# Tolerances (looser for scales involving relativistic offsets/model differences)
JD2_ATOL = Dict(
    :tai => 1e-9,
    :utc => 1e-9,
    :tt  => 1e-9,
    :tdb => 1e-9,
    :tcb => 1e-9,
    :tcg => 1e-9,
)
# Map scale symbol -> canonical jd2 truth
JD2_MAP = Dict(
    :tai => JD2_TAI,
    :utc => JD2_UTC,
    :tt  => JD2_TT,
    :tdb => JD2_TDB,
    :tcb => JD2_TCB,
    :tcg => JD2_TCG,
)

# Diagnostic helper: prints from/to scales (and diff) if comparison fails
function assert_isapprox_scale(val, ref; from::Symbol, to::Symbol, phase::Symbol, atol, rtol)
    ok = isapproxrel(val, ref; atol=atol, rtol=rtol)
    if !ok
        @info "Scale conversion mismatch" phase=phase from=from to=to got=val expected=ref diff=(val-ref) atol=atol rtol=rtol
    end
    @test ok
end

@testset "Scale forward conversions (AstroPy no round-trip)" begin
    refs = [
        (:tai, t_tai),
        (:utc, t_utc),
        (:tt,  t_tt),
        (:tdb, t_tdb),
        (:tcb, t_tcb),
        (:tcg, t_tcg),
    ]

    for (from_sym, t_from_ref) in refs
        @testset "from $(from_sym)" begin
            # Sanity of source reference vs its canonical jd2
            assert_isapprox_scale(t_from_ref.jd1, JD_TRUTH;
                from=from_sym, to=from_sym, phase=:source_jd1, atol=0.0, rtol=0.0)
            assert_isapprox_scale(t_from_ref.jd2, JD2_MAP[from_sym];
                from=from_sym, to=from_sym, phase=:source_jd2, atol=JD2_ATOL[from_sym], rtol=1e-12)

            for (to_sym, _) in refs
                from_sym == to_sym && continue
                t_to = getproperty(t_from_ref, to_sym)
                # Forward path expectation: jd1 stays partitioned same integer, jd2 matches canonical of target
                assert_isapprox_scale(t_to.jd1, JD_TRUTH;
                    from=from_sym, to=to_sym, phase=:forward_jd1, atol=0.0, rtol=0.0)
                assert_isapprox_scale(t_to.jd2, JD2_MAP[to_sym];
                    from=from_sym, to=to_sym, phase=:forward_jd2, atol=JD2_ATOL[to_sym], rtol=1e-12)
            end
        end
    end
end