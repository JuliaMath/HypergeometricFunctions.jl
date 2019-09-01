"""
This module calculates (generalized) hypergeometric functions:

    mFn(a;b;z) = Σ_{k=0}^∞ (a_1,k) ⋯ (a_m,k) / (b_1,k) ⋯ (b_n,k) zᵏ/k!
"""
module HypergeometricFunctions

using DualNumbers, SpecialFunctions

export _₂F₁, _₃F₂, mFn

include("Gauss.jl")
include("specialfunctions.jl")

end #module


