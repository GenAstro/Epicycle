
# marker_ad_permutations.jl (add/include in your test suite)

using Test
using ForwardDiff
using Zygote
using AstroStates

truth_cart = CartesianState([6759.747343616322723, 1115.043329211011041, 1344.722777534846955,
                             -2.660243619064134, 7.541202154282467, 0.640887592324028])

truth_kep = KeplerianState(8000.0, 0.2, deg2rad(12.85), deg2rad(310.0), deg2rad(120.0), deg2rad(300.0))

truth_sphradec = SphericalRADECState(6981.818181818182893, deg2rad(9.366786036853014), deg2rad(11.10476195323573), 
                                     8.022304092374013, deg2rad(109.4308653184559), deg2rad(4.582140635361989))

truth_sphazifpa = SphericalAZIFPAState(6981.818181818182893, deg2rad(9.366786036853014), deg2rad(11.10476195323573), 
                                     8.022304092374013, deg2rad(83.49318158128649), deg2rad(98.94827556462707))

truth_mee = ModifiedEquinoctialState(7680.000000000004547, 0.068404028665134, 0.187938524157182,
                                     0.072384194340782, -0.086264123652664, deg2rad(10.0))

truth_outasymptote = OutGoingAsymptoteState(6400.000000000003, -49.8250551875, deg2rad(250.633213963147),
                                            deg2rad(-11.10476195323574), deg2rad(96.50681841871349), deg2rad(300.0))

# Test case is elliptic so in and out asymptote state are the same (there is no asymptote)
truth_inasymptote = IncomingAsymptoteState(6400.000000000003, -49.8250551875, deg2rad(250.633213963147),
                                            deg2rad(-11.10476195323574), deg2rad(96.50681841871349), deg2rad(300.0));  

truth_modkep = ModifiedKeplerianState(6400.000000000003, 9600.0, deg2rad(12.85),
                                           deg2rad(310.0), deg2rad(120.0), deg2rad(300.0));  

truth_equinoct = EquinoctialState(8000.0 , 0.1879385241571815, 0.0684040286651338, -0.08626412365266437, 0.07238419434078193, deg2rad(28.36066564454829))

truth_altequinoct = AlternateEquinoctialState(8000.0 , 0.1879385241571815, 0.0684040286651338, -0.08572231482903742, 0.07192956275670796, deg2rad(28.36066564454829))


# Assumes μ and the truth_* state instances already exist in scope.
# If not, define them before including this file.
# μ = 398600.4415

# 1. Seed vectors keyed by marker instances
seed_vectors = Dict(
    Cartesian()            => to_vector(truth_cart),
    Keplerian()            => to_vector(truth_kep),
    ModifiedKeplerian()    => to_vector(truth_modkep),
    ModifiedEquinoctial()  => to_vector(truth_mee),
    Equinoctial()          => to_vector(truth_equinoct),
    AlternateEquinoctial() => to_vector(truth_altequinoct),
    SphericalRADEC()       => to_vector(truth_sphradec),
    SphericalAZIFPA()      => to_vector(truth_sphazifpa),
    #Delaunay()             => to_vector(truth_delaunay),
    IncomingAsymptote()    => to_vector(truth_inasymptote),
    OutGoingAsymptote()    => to_vector(truth_outasymptote),
)

# 2. Marker -> state constructor
state_ctor = Dict(
    Cartesian()            => CartesianState,
    Keplerian()            => KeplerianState,
    ModifiedKeplerian()    => ModifiedKeplerianState,
    ModifiedEquinoctial()  => ModifiedEquinoctialState,
    Equinoctial()          => EquinoctialState,
    AlternateEquinoctial() => AlternateEquinoctialState,
    SphericalRADEC()       => SphericalRADECState,
    SphericalAZIFPA()      => SphericalAZIFPAState,
    #Delaunay()             => DelaunayState,
    IncomingAsymptote()    => IncomingAsymptoteState,
    OutGoingAsymptote()    => OutGoingAsymptoteState,
)

# 3. Builders (vector -> source state)
function build_state(marker, v::AbstractVector, μ)
    ctor = state_ctor[marker]
    if hasmethod(ctor, Tuple{Vector{<:Real}, Real})
        return ctor(v, μ)
    elseif hasmethod(ctor, Tuple{Vector{<:Real}})
        return ctor(v)
    else
        return ctor(v...)
    end
end

# 4. Conversion (source state -> dest state)
function convert_state(dst_marker, src_state, μ)
    ctor = state_ctor[dst_marker]

    # Prefer (state, μ)
    if hasmethod(ctor, Tuple{typeof(src_state), Real})
        return ctor(src_state, μ)
    elseif hasmethod(ctor, Tuple{typeof(src_state)})
        return ctor(src_state)
    end

    v = to_vector(src_state)
    if hasmethod(ctor, Tuple{Vector{<:Real}, Real})
        return ctor(v, μ)
    elseif hasmethod(ctor, Tuple{Vector{<:Real}})
        return ctor(v)
    else
        return ctor(v...)
    end
end


# ...existing code...

# TEMP: use manual reverse Jacobian (Zygote.jacobian is choking on tuple pullbacks)
# Reason: even after vector to_vector methods, Zygote still sees the KeplerianState
# vector constructor as six scalar arguments (due to field extraction) and its
# per‑output pullbacks return NTuple{6,T}. Zygote.jacobian then tries _gradcopy!
# row <- NTuple and errors. Manual loop avoids that internal path.

zygote_rev_jac(f, x) = begin
    y = f(x)
    m = length(y); n = length(x)
    J = zeros(eltype(x), m, n)
    for i in 1:m
        g = Zygote.gradient(z -> f(z)[i], x)[1]   # returns a Vector (alloc OK for tests)
        @inbounds J[i, :] = g
    end
    J
end

tol = 1e-10
@testset "All marker instance permutations (Forward vs Reverse AD)" begin;
    for (src_marker, src_vec) in seed_vectors
        for (dst_marker, _) in seed_vectors
            @testset "$(src_marker) -> $(dst_marker)" begin
                f(x) = begin
                    s = build_state(src_marker, x, μ)
                    d = convert_state(dst_marker, s, μ)
                    to_vector(d)
                end
                J_fwd = ForwardDiff.jacobian(f, src_vec);
                J_rev = zygote_rev_jac(f, src_vec);

                abs_err = abs.(J_fwd .- J_rev);

                # Symmetric scale for relative part
                scale   = max.(abs.(J_fwd), abs.(J_rev));

                # Tunable tolerances (tight for AD vs AD; loosen for FD)
                atol = 1e-10
                rtol = 1e-10

                # Normalized error: <= 1 passes
                norm_err = abs_err ./ (atol .+ rtol .* scale);

                @test size(J_fwd) == (6,6);
                @test maximum(norm_err) ≤ 1 || begin
                    val_norm, idx_norm = findmax(norm_err)
                    val_abs,  idx_abs  = findmax(abs_err)
                    @info "Jacobian mismatch" max_norm=val_norm idx_norm=idx_norm max_abs=val_abs idx_abs=idx_abs
                    false
                end;
            end;
        end;
    end;
end;

nothing