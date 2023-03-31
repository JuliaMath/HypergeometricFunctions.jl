# HypergeometricFunctions.jl

[![Build Status](https://github.com/JuliaMath/HypergeometricFunctions.jl/workflows/CI/badge.svg)](https://github.com/JuliaMath/HypergeometricFunctions.jl/actions?query=workflow%3ACI) [![codecov](https://codecov.io/gh/JuliaMath/HypergeometricFunctions.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaMath/HypergeometricFunctions.jl)

A Julia package for calculating hypergeometric functions

This package implements the generalized hypergeometric function `pFq((α1,…,αp), (β1,…,βq), z)`. In particular, the Gauss hypergeometric function is available as `_₂F₁(a, b, c, z)`, confluent hypergeometric function is available as `_₁F₁(a, b, z) ≡ HypergeometricFunctions.M(a, b, z)` and `HypergeometricFunctions.U(a, b, z)`, as well as `_₃F₂(a1, a2, a3, b1, b2, z)`.

```julia
julia> using HypergeometricFunctions

julia> pFq((), (), 0.1) # ≡ exp(0.1)
1.1051709180756477

julia> pFq((0.5, ), (), 1.0+0.001im) # ≡ exp(-0.5*log1p(-1.0-0.001im))
22.360679774997894 + 22.36067977499789im

julia> pFq((1, ), (2, ), 0.01) # ≡ expm1(0.01)/0.01
1.0050167084168058

julia> pFq((1/3, ), (2/3, ), -1000) # ₁F₁
0.05055805394644901

julia> pFq((1, 2), (4, ), 1) # a well-poised ₂F₁
2.9999999999999996

julia> pFq((1, 2+im), (3, ), exp(im*π/3)) # ₂F₁ at that special point in ℂ
0.6036052434504436 + 0.47936537129685325im

julia> pFq((1, 2+im), (3, ), exp(im*big(π)/3)) # More digits, you say?
0.6036052434504420675589582633940696514609205969239602403078725550900852966247618 + 0.4793653712968461935342247941087557738444578105073821284003193221963949116720305im

julia> pFq((1, 2+im, 2.5), (3, 4), exp(im*π/3)) # ₃F₂ because why not
0.7994439061487193 + 0.3796162420050217im

julia> pFq((1, 1), (), -1) # A divergent series
0.5963473623232072

julia> pFq((1, 1), (), -big(1))
0.5963473623231940743410784993692793760741778601525487815734849104823272190869544

```
