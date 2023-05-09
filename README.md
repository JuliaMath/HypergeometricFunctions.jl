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
0.05055805394644902

julia> pFq((1, 2), (4, ), 1) # a well-poised ₂F₁
2.9999999999999996

julia> pFq((1, 2+im), (3.5, ), exp(im*π/3)) # ₂F₁ at that special point in ℂ
0.6786952632946589 + 0.45235049292850116im

julia> pFq((1, 2+im), (3.5, ), exp(im*big(π)/3)) # More digits, you say?
0.6786952632946589823300834090168381068073515492901393549193461972311801512528478 + 0.4523504929285013648194489713901658143893464679689810112119412310631860619948458im

julia> pFq((1, 2+im, 2.5), (3.5, 4), exp(im*π/3)) # ₃F₂ because why not
0.843443403161569 + 0.3417550761546328im

julia> pFq((1, 2+im, 2.5), (3.5, 4), exp(im*big(π)/3)) # Also in extended precision
0.8434434031615690763389963048175253868863156451003855955719081209861492349265671 + 0.3417550761546319732614495656712509723030350666571102474299311122586948108410529im

julia> pFq((1, 1), (), -1) # A divergent series
0.5963473623231935

julia> pFq((1, 1), (), -big(1))
0.5963473623231940743410784993692793760741778601525487815734849104823272191158165

```
