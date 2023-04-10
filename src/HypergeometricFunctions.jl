"""
This module calculates (generalized) hypergeometric functions:

    pFq(α, β; z) = Σ_{k=0}^∞ (α_1)ₖ ⋯ (α_p)ₖ / (β_1)ₖ ⋯ (β_q)ₖ zᵏ/k!
"""
module HypergeometricFunctions

using DualNumbers, LinearAlgebra, SpecialFunctions

export _₁F₁, _₂F₁, _₃F₂, pFq

const KMAX = 1_000_000

include("specialfunctions.jl")
include("gauss.jl")
include("confluent.jl")
include("generalized.jl")
include("drummond.jl")
include("weniger.jl")

end #module
