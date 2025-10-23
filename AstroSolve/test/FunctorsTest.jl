using Functors
using LinearAlgebra

# Minimal nested types resembling Spacecraft and Maneuver
struct SCState{T}
    state::Vector{T}
end
struct Spacecraft{T}
    state::SCState{T}
    mass::T
    id::Int
end
struct ImpulsiveManeuver{T}
    element1::T
    element2::T
    element3::T
end

# Only expose differentiable leaves via custom functors:
# - SCState: the inner state vector
# - Spacecraft: the SCState and mass (exclude id)
# - Maneuver: elements 1..3
Functors.functor(s::SCState) = ((state=s.state,), c -> SCState(c.state))
Functors.functor(sc::Spacecraft) = ((state=sc.state, mass=sc.mass), c -> Spacecraft(c.state, c.mass, sc.id))
Functors.functor(m::ImpulsiveManeuver) = ((e1=m.element1, e2=m.element2, e3=m.element3),
                                          c -> ImpulsiveManeuver(c.e1, c.e2, c.e3))

# A small context that holds subjects (this mimics where closures capture from)
struct Ctx{SC,M}
    sc::SC
    man::M
end
Functors.functor(ctx::Ctx) = ((sc=ctx.sc, man=ctx.man), c -> Ctx(c.sc, c.man))

# A closure that reads from ctx; reassignment of ctx will be observed by this closure
function cost(ctx::Ctx)
    # simple expression that touches all promoted leaves
    s = ctx.sc.state.state
    m = ctx.sc.mass
    dv = (ctx.man.element1, ctx.man.element2, ctx.man.element3)
    return norm(s) + m + dv[1]^2 + 2*dv[2] + sin(dv[3])
end

# Build initial objects
ctx = Ctx(
    Spacecraft(SCState([1.0, 2.0, 3.0]), 500.0, 42),
    ImpulsiveManeuver(0.1, -0.2, 0.3),
)

stored = () -> cost(ctx)

println("Before promotion:")
val_f64 = stored()
println("  stored() value: ", val_f64, " :: ", typeof(val_f64))

# Promote leaves to BigFloat using Functors.functor and rebuild
children, rebuild = Functors.functor(ctx)

promote_leaf(x) = x
promote_leaf(x::Real) = BigFloat(x)
promote_leaf(x::AbstractArray) = BigFloat.(x)

promoted_children = Functors.fmap(promote_leaf, children)
ctx = rebuild(promoted_children)  # reassign binding captured by stored

println("\nAfter promotion to BigFloat:")
val_big = stored()
println("  stored() value: ", val_big, " :: ", typeof(val_big))