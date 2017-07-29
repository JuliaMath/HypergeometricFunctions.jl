using ApproxFun, SingularIntegralEquations, Base.Test
    import SingularIntegralEquations: x̄sqrtx2real, sqrtx2, joukowskyinverse,
            joukowskyinversereal, joukowskyinverseabs, ⁺, ⁻, logabslegendremoment,
            stieltjeslegendremoment, stieltjesjacobimoment, stieltjesmoment, Directed,
            HypergeometricFunctions
    import HypergeometricFunctions: _₂F₁general,_₂F₁Inf


## Special functions

println("Special function tests")


for (a,b,c) in ((1,1,2),(2,2,4)), x in (1.1,10.1,1.5)
    @test _₂F₁(a,b,c,x+eps()im) ≈ _₂F₁(a,b,c,x*⁺)
    @test _₂F₁(a,b,c,x-eps()im) ≈ _₂F₁(a,b,c,x*⁻)
end


@test x̄sqrtx2real(2.0+3.0im) ≈ real(sqrtx2(2.0+3.0im)*(2.0-3.0im))

for s in (true,false)
    for z in (2.0+3.0im,2.0,0.1)
        @test real(joukowskyinverse(Val{s},z+0im)) ≈ joukowskyinversereal(Val{s},z)
        @test abs(joukowskyinverse(Val{s},z+0im)) ≈ joukowskyinverseabs(Val{s},z)
    end

    p=0.1
    @test joukowskyinverse(Val{s},p*⁺) ≈ joukowskyinverse(Val{s},p+0im)
    @test joukowskyinverse(Val{s},p*⁻) ≈ joukowskyinverse(Val{s},p-0im)
end


x=Fun()
@test sum(logabs(x-2.0)) ≈ logabslegendremoment(2.0)
@test sum(logabs(x-(2.0+im))) ≈ logabslegendremoment(2.0+im)

@test stieltjesmoment(Legendre(),0,0.1+0im) ≈ stieltjesmoment(Legendre(),0,Directed{true}(0.1))
@test stieltjesmoment(Legendre(),0,0.1-0im) ≈ stieltjesmoment(Legendre(),0,Directed{false}(0.1))
@test stieltjesjacobimoment(0.,0.,1,0.1+0im) ≈ stieltjesjacobimoment(0.,0.,1,0.1*⁺)
@test stieltjesjacobimoment(0.,0.,1,0.1-0im) ≈ stieltjesjacobimoment(0.,0.,1,0.1*⁻)

@test (stieltjesjacobimoment(0,0,0,0.1+0im)-stieltjesjacobimoment(0,0,0,0.1-0im))/(-2π*im) ≈ 1.0
@test (stieltjesjacobimoment(0,0,1,0.1+0im)-stieltjesjacobimoment(0,0,1,0.1-0im))/(-2π*im) ≈ 0.1

println("Hilbert test")
include("HilbertTest.jl")

println("Stieltjes moment test")
include("stieltjesmomentTest.jl")

println("Stieltjes integral test")
include("stieltjesintegraltest.jl")

println("Log kernel test")
include("logkerneltest.jl")

println("Curve Test")
include("CurveTest.jl")

println("FundamentalSolutions test")
include("FundamentalSolutionsTest.jl")

println("Convolution ProductFun test")
include("convolutionProductFunTest.jl")

## TODO: Move these into ApproxFun tests (since the code is there now).
println("LowRankMatrix test")
include("LowRankMatrixTest.jl")

println("GreensFun test")
include("GreensFunTest.jl")

println("Fractal test")
include("FractalTest.jl")

println("HierarchicalVector test")
include("HierarchicalVectorTest.jl")

println("Hierarchical solve test")
include("hierarchicalsolveTest.jl")

# Testing Travis CI abilities. Currently, Examples tests are slow.
#println("\nExamples test\n")
#
#include("ExamplesTest.jl")
