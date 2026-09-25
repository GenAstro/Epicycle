# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

using Test

using EpicycleBase

@testset "no_op function call" begin
    @test EpicycleBase.no_op() === nothing
end 

@testset "EpicycleBase exports" begin
    for sym in (:AbstractVar, :AbstractState, :AbstractControl, :AbstractTime, :AbstractParam,
                :AbstractFun, :AlgebraicFun, :AbstractPoint)
        @test Base.isexported(EpicycleBase, sym)
    end
end

@testset "EpicycleBase type hierarchy" begin
    @test isabstracttype(AbstractVar)
    @test isabstracttype(AbstractState)
    @test isabstracttype(AbstractControl)
    @test isabstracttype(AbstractTime)
    @test isabstracttype(AbstractParam)
    @test isabstracttype(AbstractFun)
    @test isabstracttype(AlgebraicFun)
    @test isabstracttype(AbstractPoint)

    @test AbstractState   <: AbstractVar
    @test AbstractControl <: AbstractVar
    @test AbstractTime    <: AbstractVar
    @test AbstractParam   <: AbstractVar

    @test AlgebraicFun <: AbstractFun
    @test AbstractPoint <: Any
end

@testset "EpicycleBase abstractness (non-instantiable)" begin
    @test_throws MethodError AbstractVar()
    @test_throws MethodError AbstractState()
    @test_throws MethodError AbstractControl()
    @test_throws MethodError AbstractTime()
    @test_throws MethodError AbstractParam()
    @test_throws MethodError AbstractFun()
    @test_throws MethodError AlgebraicFun()
    @test_throws MethodError AbstractPoint()
end

@testset "EpicycleBase subtyping works for user types" begin
    struct MyState    <: EpicycleBase.AbstractState   end
    struct MyControl  <: EpicycleBase.AbstractControl end
    struct MyTime     <: EpicycleBase.AbstractTime    end
    struct MyParam    <: EpicycleBase.AbstractParam   end
    struct MyAlgFun   <: EpicycleBase.AlgebraicFun    end
    struct MyPoint    <: EpicycleBase.AbstractPoint   end

    # Construct trivial instances to ensure no conflicts
    @test MyState()    isa MyState
    @test MyControl()  isa MyControl
    @test MyTime()     isa MyTime
    @test MyParam()    isa MyParam
    @test MyAlgFun()   isa MyAlgFun
    @test MyPoint()    isa MyPoint

    # And confirm subtyping
    @test MyState    <: EpicycleBase.AbstractState
    @test MyControl  <: EpicycleBase.AbstractControl
    @test MyTime     <: EpicycleBase.AbstractTime
    @test MyParam    <: EpicycleBase.AbstractParam
    @test MyAlgFun   <: EpicycleBase.AlgebraicFun
    @test MyPoint    <: EpicycleBase.AbstractPoint
end

# =============================================================================
# Tag / Variable System
# =============================================================================

@testset "AbstractVarTag exports" begin
    for sym in (:AbstractVarTag, :AbstractStateTag, :AbstractParamTag,
                :AbstractControlTag, :AbstractTimeTag,
                :ModelVariable, :DirectVariable,
                :get_field, :set_field!)
        @test Base.isexported(EpicycleBase, sym)
    end
end

@testset "AbstractVarTag hierarchy" begin
    @test isabstracttype(AbstractVarTag)
    @test isabstracttype(AbstractStateTag)
    @test isabstracttype(AbstractParamTag)
    @test isabstracttype(AbstractControlTag)
    @test isabstracttype(AbstractTimeTag)

    @test AbstractStateTag   <: AbstractVarTag
    @test AbstractParamTag   <: AbstractVarTag
    @test AbstractControlTag <: AbstractVarTag
    @test AbstractTimeTag    <: AbstractVarTag
end

@testset "AbstractVarTag concrete subtypes" begin
    struct MyStateTag   <: AbstractStateTag   end
    struct MyParamTag   <: AbstractParamTag   end
    struct MyControlTag <: AbstractControlTag end

    @test MyStateTag()   isa AbstractStateTag
    @test MyParamTag()   isa AbstractParamTag
    @test MyControlTag() isa AbstractControlTag
    @test MyStateTag()   isa AbstractVarTag
    @test MyParamTag()   isa AbstractVarTag
end

@testset "ModelVariable construction and identity" begin
    struct DummyTag <: AbstractParamTag end
    mutable struct DummyModel; val::Float64 end

    m1 = DummyModel(1.0)
    m2 = DummyModel(2.0)
    tag = DummyTag()

    v1a = ModelVariable(m1, tag)
    v1b = ModelVariable(m1, tag)
    v2  = ModelVariable(m2, tag)

    # Same model instance + same tag type -> identical key
    @test (objectid(v1a.model), typeof(v1a.tag)) == (objectid(v1b.model), typeof(v1b.tag))
    # Different model instance -> different key
    @test (objectid(v1a.model), typeof(v1a.tag)) != (objectid(v2.model),  typeof(v2.tag))

    @test v1a.model === m1
    @test v1a.tag   isa DummyTag
end

@testset "DirectVariable" begin
    dv = DirectVariable(value=10.0, lower_bound=-100.0, upper_bound=100.0, name="dvx")
    @test dv.value        ≈  10.0
    @test dv.lower_bound ≈ -100.0
    @test dv.upper_bound ≈  100.0
    @test dv.name         == "dvx"
end

@testset "get_field / set_field! convenience via ModelVariable" begin
    struct DummyTag2 <: AbstractParamTag end
    mutable struct DummyModel2; val::Float64 end

    EpicycleBase.get_field(m::DummyModel2, ::DummyTag2)           = m.val
    EpicycleBase.set_field!(m::DummyModel2, ::DummyTag2, v::Real) = (m.val = v; nothing)

    m   = DummyModel2(3.14)
    var = ModelVariable(m, DummyTag2())

    @test get_field(var)   ≈ 3.14
    set_field!(var, 2.71)
    @test get_field(var)   ≈ 2.71
    @test m.val            ≈ 2.71
end

@testset "ModelVariable equality follows the model's identity" begin
    # The rule the docstring states: same tag type, and an identical model. For a mutable model
    # that is the same instance; for an immutable one it is equal field values.
    struct EqTag      <: AbstractParamTag end
    struct OtherEqTag <: AbstractParamTag end
    mutable struct EqModel; val::Float64 end

    a = EqModel(1.0)
    b = EqModel(1.0)                                   # same type, same value, another instance
    @test ModelVariable(a, EqTag()) == ModelVariable(a, EqTag())
    @test ModelVariable(a, EqTag()) != ModelVariable(b, EqTag())
    @test ModelVariable(a, EqTag()) != ModelVariable(a, OtherEqTag())

    keys_mutable = Dict(ModelVariable(a, EqTag()) => 1, ModelVariable(b, EqTag()) => 2)
    @test length(keys_mutable) == 2

    # Two separately built NamedTuples with the same values are one variable.
    n1 = ModelVariable((val = 1.0,), EqTag())
    n2 = ModelVariable((val = 1.0,), EqTag())
    @test n1 == n2
    @test length(Dict(n1 => 1, n2 => 2)) == 1
    @test ModelVariable((val = 1.0,), EqTag()) != ModelVariable((val = 2.0,), EqTag())
end

@testset "DirectVariable converts its inputs and requires every keyword" begin
    dv = DirectVariable(value = 3, lower_bound = 1, upper_bound = 10, name = SubString("tof x", 1, 3))
    @test dv.value === 3.0 && dv.lower_bound === 1.0 && dv.upper_bound === 10.0
    @test dv.name === "tof"

    @test_throws UndefKeywordError DirectVariable(value = 1.0, lower_bound = 0.0, upper_bound = 2.0)
    @test_throws UndefKeywordError DirectVariable(lower_bound = 0.0, upper_bound = 2.0, name = "x")
    @test_throws TypeError DirectVariable(value = "1", lower_bound = 0.0, upper_bound = 2.0,
                                          name = "x")
end

@testset "get_field / set_field! carry a vector field through a ModelVariable" begin
    # A state tag's value is a vector, as AstroProp's PosVel is.
    struct VecStateTag <: AbstractStateTag end
    mutable struct VecModel; posvel::Vector{Float64} end

    EpicycleBase.get_field(m::VecModel, ::VecStateTag) = m.posvel
    EpicycleBase.set_field!(m::VecModel, ::VecStateTag, v::AbstractVector) = (m.posvel = v; nothing)

    m   = VecModel([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0])
    var = ModelVariable(m, VecStateTag())

    @test get_field(var) == [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]
    set_field!(var, [7100.0, 1.0, 2.0, 0.1, 7.4, 0.2])
    @test m.posvel == [7100.0, 1.0, 2.0, 0.1, 7.4, 0.2]
end

include("test_correctness_quantity_traits.jl")

nothing
