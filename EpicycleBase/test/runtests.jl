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
nothing