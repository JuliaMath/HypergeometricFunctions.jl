module HypergeometricFunctions

using DualNumbers, LinearAlgebra, SpecialFunctions

export _₁F₁, _₂F₁, _₃F₂, pFq

const KMAX = 1_048_576

include("specialfunctions.jl")
include("gauss.jl")
include("confluent.jl")
include("generalized.jl")
include("drummond.jl")
include("weniger.jl")

end #module
