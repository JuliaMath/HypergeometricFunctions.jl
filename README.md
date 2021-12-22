# HypergeometricFunctions.jl

[![Build Status](https://github.com/JuliaMath/HypergeometricFunctions.jl/workflows/CI/badge.svg)](https://github.com/JuliaMath/HypergeometricFunctions.jl/actions?query=workflow%3ACI) [![codecov](https://codecov.io/gh/JuliaMath/HypergeometricFunctions.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaMath/HypergeometricFunctions.jl)

A Julia package for calculating hypergeometric functions

This package implements the generalized hypergeometric function `pFq([a1,…,am], [b1,…,bn], z)`. In particular, the Gauss hypergeometric function is available as `_₂F₁(a, b, c, z)`, confluent hypergeometric function is available as `_₁F₁(a, b, z) ≡ HypergeometricFunctions.M(a, b, z)` and `HypergeometricFunctions.U(a, b, z)`, as well as `_₃F₂([a1, a2, a3], [b1, b2], z)`.
