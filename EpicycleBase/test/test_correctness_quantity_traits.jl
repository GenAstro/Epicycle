# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# =============================================================================
# Traits on quantity functions.
#
# The claim being tested is that a trait on `typeof(f)` behaves exactly as one
# on a struct would — which is what lets a quantity be a plain function and
# still carry a label and a wrap period. If that fails, quantities need wrapper
# types instead.

# Settability is the exception, and deliberately so: it is not a trait on the
# quantity alone. A quantity is settable on a particular kind of subject, so
# `set_quantity!` dispatches on the pair and `is_settable` asks whether that
# pair has a method.
#
# Two properties matter beyond the values coming back right. Declaring a trait
# must be the whole opt-in, with no registry, or a user's own quantity is a
# second-class citizen. And the traits must resolve at compile time, or every
# `is_settable` check in a solver loop costs a lookup.
# =============================================================================

using Test
using InteractiveUtils
using EpicycleBase

# A quantity that declares nothing, and one that declares everything.
undeclared(subject) = subject
declared(subject) = subject

# Two subjects, so settability can differ between them for the same quantity.
mutable struct _Writable; value::Float64 end
struct _ReadOnly end

EpicycleBase.label(::typeof(declared))     = "Declared quantity"
EpicycleBase.is_cyclic(::typeof(declared)) = true
EpicycleBase.cycle(::typeof(declared))     = 2π

EpicycleBase.set_quantity!(s::_Writable, ::typeof(declared); to) = (s.value = to)

# A quantity read with an optional frame, and its writer declared with the frame defaulted: one
# method definition, two positional shapes, and a typed value.
struct _Frame2 end
framed(s, f::_Frame2 = _Frame2()) = s
EpicycleBase.set_quantity!(s::_Writable, ::typeof(framed), f::_Frame2 = _Frame2(); to::Real) =
    (s.value = to)

struct _ATag <: AbstractParamTag end
tagged(subject) = subject
EpicycleBase.tag(::typeof(tagged)) = _ATag()

@testset "an undeclared quantity still works everywhere it is read" begin
    # No registration, no subtyping — a bare function has every trait, at its
    # default. Nothing has to be declared for a quantity to be readable.
    @test label(undeclared)       == "quantity"
    @test is_cyclic(undeclared)   == false
    @test cycle(undeclared)       === nothing
    @test tag(undeclared)         === nothing
    @test is_settable(_Writable(0.0), undeclared) == false
end

@testset "declaring a method is the whole opt-in" begin
    @test label(declared)     == "Declared quantity"
    @test is_cyclic(declared) == true
    @test cycle(declared)     ≈ 2π
    @test tag(tagged)         isa _ATag
end

@testset "is_settable is derived, so the two cannot disagree" begin
    # The method is the primitive. There is no separate is_settable to fall out
    # of step with it — which is the failure a boolean field would invite.
    w = _Writable(0.0)
    @test is_settable(w, declared) == hasmethod(EpicycleBase.set_quantity!,
                                                Tuple{_Writable, typeof(declared)})
    @test is_settable(w, declared) == true

    set_quantity!(w, declared; to = 3.5)
    @test w.value == 3.5
end

@testset "settability belongs to the pair, not to the quantity" begin
    # The same quantity, two subjects, two answers. A trait on the quantity
    # alone could not express this, which is why `setter` was retired.
    @test is_settable(_Writable(0.0), declared) == true
    @test is_settable(_ReadOnly(),    declared) == false

    @test_throws MethodError set_quantity!(_ReadOnly(), declared; to = 1.0)
end

@testset "is_settable answers only for the argument shapes a writer accepts" begin
    # A default argument declares each arity it covers, so the frame can be given or left out.
    # Another type in the frame's place has no method, and is_settable says so rather than
    # answering yes to a call the writer would refuse, which is what a `deps...` signature does.
    w = _Writable(0.0)
    @test is_settable(w, framed)
    @test is_settable(w, framed, _Frame2())
    @test !is_settable(w, framed, 1.0)

    set_quantity!(w, framed, _Frame2(); to = 2.0)
    @test w.value == 2.0

    # The value is a keyword, outside the probed signature, so it can be typed.
    @test_throws TypeError set_quantity!(w, framed; to = "two")
end

@testset "traits resolve at compile time" begin
    # A solver checks settability once per variable per iteration; a runtime
    # lookup there is a real cost. These must fold to a constant.
    settable()   = is_settable(_Writable(0.0), declared)
    unsettable() = is_settable(_Writable(0.0), undeclared)
    cyclic()     = is_cyclic(declared)

    folds_to(f, needle) = occursin(needle,
        sprint(io -> code_llvm(io, f, Tuple{}; debuginfo = :none)))

    @test folds_to(settable,   "ret i8 1")
    @test folds_to(unsettable, "ret i8 0")
    @test folds_to(cyclic,     "ret i8 1")
end

@testset "has_output_partial is derived from output_partial, per subject and deps" begin
    # The same rule as is_settable: the method is the declaration, and a quantity read with an
    # extra argument declares a separate method.
    struct _Frame end
    EpicycleBase.output_partial(::_ReadOnly, ::typeof(declared))          = zeros(1, 6)
    EpicycleBase.output_partial(::_ReadOnly, ::typeof(declared), ::_Frame) = ones(1, 6)

    @test has_output_partial(_ReadOnly(), declared)
    @test has_output_partial(_ReadOnly(), declared, _Frame())
    @test !has_output_partial(_ReadOnly(), undeclared)
    @test !has_output_partial(_Writable(0.0), declared)             # another subject
    @test !has_output_partial(_ReadOnly(), declared, 1.0)           # another dep type

    @test output_partial(_ReadOnly(), declared)           == zeros(1, 6)
    @test output_partial(_ReadOnly(), declared, _Frame()) == ones(1, 6)
end

@testset "traits dispatch on the function, not on a value" begin
    # `typeof(f)` is a singleton type, so two distinct functions never collide
    # and a trait declared for one says nothing about the other.
    @test typeof(declared) !== typeof(undeclared)
    @test label(declared)  != label(undeclared)
end
