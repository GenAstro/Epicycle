include("state_truthdata_elliptic_orbits.jl")

percent_error(x, y) = begin
    denom = abs(y) < 1 ? 1.0 : abs(y)
    abs(x - y) / denom
end
isapproxvec_percent(x::AbstractVector, y::AbstractVector; tol=1e-12) =
    maximum(percent_error.(x, y)) ≤ tol

state_pairs = [
    (Cartesian(),            CartesianState),
    (Keplerian(),            KeplerianState),
    (ModifiedKeplerian(),    ModifiedKeplerianState),
    (Equinoctial(),          EquinoctialState),
    (AlternateEquinoctial(), AlternateEquinoctialState),
    (ModifiedEquinoctial(),  ModifiedEquinoctialState),
    (SphericalRADEC(),       SphericalRADECState),
    (SphericalAZIFPA(),      SphericalAZIFPAState),
    (IncomingAsymptote(),    IncomingAsymptoteState),
    (OutGoingAsymptote(),    OutGoingAsymptoteState),
]

@testset "state_tag_to_type (tag → type) hits" begin
    for (tag, Tstate) in state_pairs
        # Direct call
        @test state_tag_to_type(tag) === Tstate
        # Also via invokelatest to ensure a dynamic call path gets counted
        @test Base.invokelatest(AstroStates.state_tag_to_type, tag) === Tstate
    end
end

@testset "state_tag_to_type and state_type_to_tag mapping" begin

    # Helper to build a minimal instance of a given state type (Float64)
    mk(::Type{CartesianState})            = CartesianState([1.0,2.0,3.0,4.0,5.0,6.0])
    mk(::Type{KeplerianState})            = KeplerianState([7000.0, 0.01, 0.2, 0.3, 0.4, 0.5])
    mk(::Type{ModifiedKeplerianState})    = ModifiedKeplerianState([300.0, 8000.0, 0.2, 0.3, 0.4, 0.5])
    mk(::Type{EquinoctialState})          = EquinoctialState([7000.0, 0.01, 0.02, 0.1, 0.2, 0.3])
    mk(::Type{AlternateEquinoctialState}) = AlternateEquinoctialState([7000.0, 0.01, 0.02, 0.1, 0.2, 0.3])
    mk(::Type{ModifiedEquinoctialState})  = ModifiedEquinoctialState([7000.0, 0.01, 0.02, 0.1, 0.2, 0.3])
    mk(::Type{SphericalRADECState})       = SphericalRADECState([7000.0, 0.1, 0.2, 7.5, 0.1, 0.2])
    mk(::Type{SphericalAZIFPAState})      = SphericalAZIFPAState([7000.0, 0.2, 0.1, 7.5, 0.3, 0.05])
    mk(::Type{IncomingAsymptoteState})    = IncomingAsymptoteState([300.0, 10.0, 0.2, 0.1, 0.05, 2.5])
    mk(::Type{OutGoingAsymptoteState})    = OutGoingAsymptoteState([300.0, 10.0, 0.2, 0.1, 0.05, 2.5])

    for (tag, Tstate) in state_pairs
        # Tag → Type
        @test state_tag_to_type(tag) === Tstate

        # Type → Tag (UnionAll or concrete)
        @test state_type_to_tag(Tstate) == tag

        # Parametric instantiation only when Tstate is parametric
        if Tstate isa UnionAll
            @test state_type_to_tag(Tstate{Float64}) == tag
        end

        # Instance → Tag (exercise instance overload once, avoid round-trip recursion)
        inst = mk(Tstate)
        @test state_type_to_tag(inst) == tag

        # Round-trip via non-instance path (avoids instance→typeof recursion)
        @test state_tag_to_type(state_type_to_tag(Tstate)) === Tstate
    end
end

@testset "force-hit state_tag_to_type one-liners" begin
    @test Core.invoke(AstroStates.state_tag_to_type, Tuple{AstroStates.Cartesian},            AstroStates.Cartesian())            === AstroStates.CartesianState
    @test Core.invoke(AstroStates.state_tag_to_type, Tuple{AstroStates.Keplerian},            AstroStates.Keplerian())            === AstroStates.KeplerianState
    @test Core.invoke(AstroStates.state_tag_to_type, Tuple{AstroStates.ModifiedKeplerian},    AstroStates.ModifiedKeplerian())    === AstroStates.ModifiedKeplerianState
    @test Core.invoke(AstroStates.state_tag_to_type, Tuple{AstroStates.Equinoctial},          AstroStates.Equinoctial())          === AstroStates.EquinoctialState
    @test Core.invoke(AstroStates.state_tag_to_type, Tuple{AstroStates.AlternateEquinoctial}, AstroStates.AlternateEquinoctial()) === AstroStates.AlternateEquinoctialState
    @test Core.invoke(AstroStates.state_tag_to_type, Tuple{AstroStates.ModifiedEquinoctial},  AstroStates.ModifiedEquinoctial())  === AstroStates.ModifiedEquinoctialState
    @test Core.invoke(AstroStates.state_tag_to_type, Tuple{AstroStates.SphericalRADEC},       AstroStates.SphericalRADEC())       === AstroStates.SphericalRADECState
    @test Core.invoke(AstroStates.state_tag_to_type, Tuple{AstroStates.SphericalAZIFPA},      AstroStates.SphericalAZIFPA())      === AstroStates.SphericalAZIFPAState
    @test Core.invoke(AstroStates.state_tag_to_type, Tuple{AstroStates.IncomingAsymptote},    AstroStates.IncomingAsymptote())    === AstroStates.IncomingAsymptoteState
    @test Core.invoke(AstroStates.state_tag_to_type, Tuple{AstroStates.OutGoingAsymptote},    AstroStates.OutGoingAsymptote())    === AstroStates.OutGoingAsymptoteState
end

@testset "force-hit state_type_to_tag instance overload" begin
    # Use one concrete instance to exercise st::AbstractOrbitState method
    inst = AstroStates.CartesianState([1.0,2.0,3.0,4.0,5.0,6.0])
    @test Core.invoke(AstroStates.state_type_to_tag, Tuple{AstroStates.AbstractOrbitState}, inst) == AstroStates.Cartesian()
end
nothing