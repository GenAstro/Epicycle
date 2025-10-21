

include("state_truthdata_elliptic_orbits.jl")

@testset "Kinematic mu-less round-trip conversions" begin
    kinematic_truth = Dict(
        :CartesianState          => truth_cart,
        :SphericalAZIFPAState    => truth_sphazifpa,
        :SphericalRADECState     => truth_sphradec,
    )

    approx(v1, v2) = isapproxvec_percent(to_vector(v1), to_vector(v2); tol=1e-12)

    types_seq = [
        (CartesianState, SphericalAZIFPAState),
        (CartesianState, SphericalRADECState),
        (SphericalAZIFPAState, CartesianState),
        (SphericalAZIFPAState, SphericalRADECState),
        (SphericalRADECState, CartesianState),
        (SphericalRADECState, SphericalAZIFPAState),
    ]

    for (FromT, ToT) in types_seq
        from_sym = Symbol(nameof(FromT))
        haskey(kinematic_truth, from_sym) || continue
        src = kinematic_truth[from_sym]

        @testset "$(nameof(FromT)) -> $(nameof(ToT)) (no μ)" begin
            if hasmethod(ToT, Tuple{FromT})
                tgt = ToT(src)
                @test tgt isa ToT                      # was: typeof(tgt) === ToT
                if hasmethod(FromT, Tuple{ToT})
                    back = FromT(tgt)
                    @test back isa FromT               # was: typeof(back) === FromT
                    @test approx(back, src)
                end
            else
                @info "Skipping $(nameof(FromT)) -> $(nameof(ToT)) (no μ constructor not defined)"
            end
        end
    end
end


nothing