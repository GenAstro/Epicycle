# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

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
                :get_field, :set_field!, :differentiate_wrt, :fd_differentiate_wrt)
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
    dv = DirectVariable(value=10.0, lower_bounds=-100.0, upper_bounds=100.0, name="dvx")
    @test dv.value        ≈  10.0
    @test dv.lower_bounds ≈ -100.0
    @test dv.upper_bounds ≈  100.0
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

@testset "differentiate_wrt default and fd_differentiate_wrt oracle" begin
    # Scalar model: f(p) = p^2, so df/dp = 2p
    struct ScalarTag <: AbstractParamTag end
    mutable struct ScalarModel; p::Float64 end

    EpicycleBase.get_field(m::ScalarModel, ::ScalarTag)           = m.p
    EpicycleBase.set_field!(m::ScalarModel, ::ScalarTag, v::Real) = (m.p = float(v); nothing)

    m   = ScalarModel(3.0)
    rhs = () -> [m.p^2]       # returns a 1-element vector

    deriv_default = differentiate_wrt(rhs, m, ScalarTag())
    deriv_fd      = fd_differentiate_wrt(rhs, m, ScalarTag())

    @test length(deriv_default) == 1
    @test isapprox(deriv_default[1], 2 * 3.0; rtol=1e-8)   # 2p = 6.0
    @test isapprox(deriv_fd[1],      2 * 3.0; rtol=1e-8)
    @test isapprox(deriv_default,    deriv_fd; rtol=1e-10)

    # model field is restored after differentiation
    @test m.p ≈ 3.0
end

nothing