# HypergeometricFunctions.jl

[![Build Status](https://github.com/JuliaMath/HypergeometricFunctions.jl/workflows/CI/badge.svg)](https://github.com/JuliaMath/HypergeometricFunctions.jl/actions?query=workflow%3ACI) [![codecov](https://codecov.io/gh/JuliaMath/HypergeometricFunctions.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaMath/HypergeometricFunctions.jl) [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaMath.github.io/HypergeometricFunctions.jl/stable) [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaMath.github.io/HypergeometricFunctions.jl/dev)

This package provides an implementation of the generalized hypergeometric function `pFq(α, β, z)`.

```julia
julia> using HypergeometricFunctions

julia> pFq((), (), 0.1) # ≡ exp(0.1)
1.1051709180756477

julia> pFq((0.5, ), (), 1.0+0.001im) # ≡ exp(-0.5*log1p(-1.0-0.001im))
22.360679774997894 + 22.36067977499789im

julia> pFq((), (1.5, ), -π^2/4) # A root of a spherical Bessel function
4.042865030283967e-17

julia> pFq((1/3, ), (2/3, ), -1000) # A confluent hypergeometric with large argument
0.050558053946448855

julia> pFq((1, 2+im), (3.5, ), exp(im*π/3)) # ₂F₁ at that special point in ℂ
0.6786952632946592 + 0.4523504929285015im

```

# References

[1] N. Michel and M. V. Stoitsov, [Fast computation of the Gauss hypergeometric function with all its parameters complex with application to the Pöschl–Teller–Ginocchio potential wave functions](https://doi.org/10.1016/j.cpc.2007.11.007), *Comp. Phys. Commun.*, **178**:535–551, 2008.

[2] J. W. Pearson, S. Olver and M. A. Porter, [Numerical methods for the computation of the confluent and Gauss hypergeometric functions](https://doi.org/10.1007/s11075-016-0173-0), *Numer. Algor.*, **74**:821–866, 2017.

[3] R. M. Slevinsky, [Fast and stable rational approximation of generalized hypergeometric functions](https://arxiv.org/abs/2307.06221), arXiv:2307.06221, 2023.
