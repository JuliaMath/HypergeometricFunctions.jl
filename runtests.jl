using ApproxFun, SingularIntegralEquations, Base.Test
    import SingularIntegralEquations: x̄sqrtx2real, sqrtx2, joukowskyinverse,
            joukowskyinversereal, joukowskyinverseabs

## Special functions

println("Special function tests")

@test_approx_eq x̄sqrtx2real(2.0+3.0im) real(sqrtx2(2.0+3.0im)*(2.0-3.0im))

for s in (true,false), z in (2.0+3.0im,2.0,0.1)
    @test_approx_eq real(joukowskyinverse(Val{s},z+0im)) joukowskyinversereal(Val{s},z)
    @test_approx_eq abs(joukowskyinverse(Val{s},z+0im)) joukowskyinverseabs(Val{s},z)
end


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

println("LowRankMatrix test")
include("LowRankMatrixTest.jl")

println("HierarchicalVector test")
include("HierarchicalVectorTest.jl")

println("Hierarchical solve test")
include("hierarchicalsolveTest.jl")

# Testing Travis CI abilities. Currently, Examples tests are slow.
#println("\nExamples test\n")
#
#include("ExamplesTest.jl")

println("WienerHopfTest")
include("WienerHopfTest.jl")


println("Ideal Fluid Flow tests")
include("IdealFluidFlowTest.jl")
