using ApproxFun, SingularIntegralEquations, Base.Test

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
