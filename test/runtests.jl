using HypergeometricFunctions, SpecialFunctions, Test
import LinearAlgebra: norm
import HypergeometricFunctions: iswellpoised, isalmostwellpoised, M, U,
                                pochhammer, _₂F₁general, _₂F₁general2,
                                pFqdrummond, pFqweniger, pFq2string

const rtol = 1.0e-3
const NumberType = Float64

include("specialfunctions_tests.jl")
include("algorithm_tests.jl")
include("confluent_tests.jl")
include("regression_tests.jl")
