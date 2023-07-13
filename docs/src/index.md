# HypergeometricFunctions.jl Documentation

## Introduction

[`HypergeometricFunctions.jl`](https://github.com/JuliaMath/HypergeometricFunctions.jl) provides a numerical computation of generalized hypergeometric functions. The main exported function and recommended interface is [`pFq`](@ref), but there are a few others for specialists' convenience.

```@docs
pFq
```

## Complex phase portraits

Broadly speaking, there are three classes of generalized hypergeometric functions: when $p\le q$ they are entire functions of the complex variable $z$; when $p = q+1$, they are analytic functions in the cut plane $\mathbb{C}\setminus[1,\infty)$; and, when $p > q+1$, they are analytic functions in the cut plane $\mathbb{C}\setminus[0,\infty)$.

Examples of each of these classes are illustrated over $\left\{z\in\mathbb{C} : -10<\Re z<10, -10<\Im z<10\right\}$ with [complex phase portraits](https://en.wikipedia.org/wiki/Domain_coloring), a beautiful tool in computational complex analysis.

```@example
using ComplexPhasePortrait, HypergeometricFunctions, Images
x = range(-10, stop=10, length=300)
y = range(-10, stop=10, length=300)
z = x' .+ im*y

import Logging # To avoid printing warnings
Logging.with_logger(Logging.SimpleLogger(Logging.Error)) do
    img = portrait(map(z->pFq((), (), z), z), ctype = "nist")
    save("0F0.png", img)
    img = portrait(map(z->pFq((), (1.0, ), z), z), ctype = "nist")
    save("0F1.png", img)
    img = portrait(map(z->pFq((0.5, ), (0.75, ), z), z), ctype = "nist")
    save("1F1.png", img)
    img = portrait(map(z->pFq((3.5+7.5im, ), (), z), z), ctype = "nist")
    save("1F0.png", img)
    img = portrait(map(z->pFq((1.0, 3.5+7.5im), (0.75, ), z), z), ctype = "nist")
    save("2F1.png", img)
    img = portrait(map(z->pFq((1.0, 1.5+7.5im), (), z), z), ctype = "nist")
    save("2F0.png", img)
end
nothing # hide
```

|p\q | 0 | 1 |
| :---: | :---: | :---: |
| 0 | ![₀F₀](0F0.png) | ![₀F₁](0F1.png) |
| 1 | ![₁F₀](1F0.png) | ![₁F₁](1F1.png) |
| 2 | ![₂F₀](2F0.png) | ![₂F₁](2F1.png) |

## Library

```@docs
_₁F₁
_₂F₁
_₃F₂
```

## Internals

```@docs
HypergeometricFunctions.M
HypergeometricFunctions.U
HypergeometricFunctions._₂F₁positive
HypergeometricFunctions._₂F₁general
HypergeometricFunctions._₂F₁general2
HypergeometricFunctions.pFqdrummond
HypergeometricFunctions.pFqweniger
HypergeometricFunctions.pFqcontinuedfraction
HypergeometricFunctions.pochhammer
HypergeometricFunctions.@clenshaw
HypergeometricFunctions.@lanczosratio
HypergeometricFunctions.G
HypergeometricFunctions.P
```
