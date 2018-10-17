using ApproxFun, SingularIntegralEquations, SpecialFunctions, LinearAlgebra, Test
    import SingularIntegralEquations: x̄sqrtx2real, sqrtx2, joukowskyinverse,
            joukowskyinversereal, joukowskyinverseabs, ⁺, ⁻, logabslegendremoment,
            stieltjeslegendremoment, stieltjesjacobimoment, stieltjesmoment, Directed,
            HypergeometricFunctions
    import SingularIntegralEquations.HypergeometricFunctions: _₂F₁general,_₂F₁Inf


## Special functions

@testset "Special function tests" begin
    for (a,b,c) in ((1,1,2),(2,2,4)), x in (1.1,10.1,1.5)
        @test _₂F₁(a,b,c,x+eps()im) ≈ _₂F₁(a,b,c,(x)⁺)
        @test _₂F₁(a,b,c,x-eps()im) ≈ _₂F₁(a,b,c,(x)⁻)
    end

    @test x̄sqrtx2real(2.0+3.0im) ≈ real(sqrtx2(2.0+3.0im)*(2.0-3.0im))

    for s in (true,false)
        for z in (2.0+3.0im,2.0,0.1)
            @test real(joukowskyinverse(Val{s},z+0im)) ≈ joukowskyinversereal(Val{s},z)
            @test abs(joukowskyinverse(Val{s},z+0im)) ≈ joukowskyinverseabs(Val{s},z)
        end

        p=0.1
        @test joukowskyinverse(Val{s},(p)⁺) ≈ joukowskyinverse(Val{s},p+0im)
        @test joukowskyinverse(Val{s},(p)⁻) ≈ joukowskyinverse(Val{s},p-0im)
    end


    x=Fun()
    @test sum(logabs(x-2.0)) ≈ logabslegendremoment(2.0)
    @test sum(logabs(x-(2.0+im))) ≈ logabslegendremoment(2.0+im)

    @test sqrt(Directed{true}(-0.1)) ≈ sqrt(-0.1-0.0im)
    @test sqrt(Directed{false}(-0.1)) ≈ sqrt(-0.1+0.0im)

    @test stieltjesmoment(Legendre(),0,0.1+0im) ≈ stieltjesmoment(Legendre(),0,Directed{true}(0.1))
    @test stieltjesmoment(Legendre(),0,0.1-0im) ≈ stieltjesmoment(Legendre(),0,Directed{false}(0.1))
    @test stieltjesjacobimoment(0.,0.,1,0.1+0im) ≈ stieltjesjacobimoment(0.,0.,1,0.1⁺)
    @test stieltjesjacobimoment(0.,0.,1,0.1-0im) ≈ stieltjesjacobimoment(0.,0.,1,0.1⁻)

    @test (stieltjesjacobimoment(0,0,0,0.1+0im)-stieltjesjacobimoment(0,0,0,0.1-0im))/(-2π*im) ≈ 1.0
    @test (stieltjesjacobimoment(0,0,1,0.1+0im)-stieltjesjacobimoment(0,0,1,0.1-0im))/(-2π*im) ≈ 0.1

    @test stieltjesjacobimoment(0.5,0,0,Directed{true}(0.1)) ≈ stieltjesjacobimoment(0.5,0,0,0.1+0im)
    @test stieltjesjacobimoment(0.5,0,0,Directed{false}(0.1)) ≈ stieltjesjacobimoment(0.5,0,0,0.1-0.0im)
end

include("HilbertTest.jl")
include("stieltjesmomentTest.jl")
include("stieltjesintegraltest.jl")
include("logkerneltest.jl")
include("CurveTest.jl")
include("FundamentalSolutionsTest.jl")
include("convolutionProductFunTest.jl")
include("GreensFunTest.jl")
include("NonlocalOperatorsTest.jl")


# Testing Travis CI abilities. Currently, Examples tests are slow.
#println("\nExamples test\n")
#
#include("ExamplesTest.jl")
