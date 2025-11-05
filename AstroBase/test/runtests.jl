# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

using Test

using AstroBase

@testset "no_op function call" begin
    @test AstroBase.no_op() === nothing
end 

@testset "AstroBase exports" begin
    for sym in (:AbstractVar, :AbstractState, :AbstractControl, :AbstractTime, :AbstractParam,
                :AbstractFun, :AlgebraicFun, :AbstractPoint)
        @test Base.isexported(AstroBase, sym)
    end
end

@testset "AstroBase type hierarchy" begin
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

@testset "AstroBase abstractness (non-instantiable)" begin
    @test_throws MethodError AbstractVar()
    @test_throws MethodError AbstractState()
    @test_throws MethodError AbstractControl()
    @test_throws MethodError AbstractTime()
    @test_throws MethodError AbstractParam()
    @test_throws MethodError AbstractFun()
    @test_throws MethodError AlgebraicFun()
    @test_throws MethodError AbstractPoint()
end

@testset "AstroBase subtyping works for user types" begin
    struct MyState    <: AstroBase.AbstractState   end
    struct MyControl  <: AstroBase.AbstractControl end
    struct MyTime     <: AstroBase.AbstractTime    end
    struct MyParam    <: AstroBase.AbstractParam   end
    struct MyAlgFun   <: AstroBase.AlgebraicFun    end
    struct MyPoint    <: AstroBase.AbstractPoint   end

    # Construct trivial instances to ensure no conflicts
    @test MyState()    isa MyState
    @test MyControl()  isa MyControl
    @test MyTime()     isa MyTime
    @test MyParam()    isa MyParam
    @test MyAlgFun()   isa MyAlgFun
    @test MyPoint()    isa MyPoint

    # And confirm subtyping
    @test MyState    <: AstroBase.AbstractState
    @test MyControl  <: AstroBase.AbstractControl
    @test MyTime     <: AstroBase.AbstractTime
    @test MyParam    <: AstroBase.AbstractParam
    @test MyAlgFun   <: AstroBase.AlgebraicFun
    @test MyPoint    <: AstroBase.AbstractPoint
end
nothing