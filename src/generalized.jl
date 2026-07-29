@doc raw"""
    pFq(α, β, z)
Compute the generalized hypergeometric function, defined by
```math
{}_pF_q(α, β, z) = \sum_{k=0}^\infty \dfrac{(\alpha_1)_k\cdots(\alpha_p)_k}{(\beta_1)_k\cdots(\beta_q)_k}\dfrac{z^k}{k!},
```
where the series converges and elsewhere by analytic continuation.

External links: [DLMF](https://dlmf.nist.gov/16.2.1), [Wikipedia](https://en.wikipedia.org/wiki/Generalized_hypergeometric_function).

# Examples
```jldoctest
julia> pFq((), (), 0.1) # ≡ exp(0.1)
1.1051709180756477

julia> pFq((0.5, ), (), 1.0+0.001im) # ≡ exp(-0.5*log1p(-1.0-0.001im))
22.360679774997894 + 22.36067977499789im

julia> pFq((), (1.5, ), -π^2/4) # A root of a spherical Bessel function
4.042865030283967e-17

julia> pFq((), (1.5, ), -big(π)^2/4) # In extended precision
8.674364372518869408017614476652675180406967418943475242812160199356160822272727e-78

julia> pFq((1, ), (2, ), 0.01) # ≡ expm1(0.01)/0.01
1.0050167084168058

julia> pFq((1/3, ), (2/3, ), -1000) # A confluent hypergeometric with large argument
0.050558053946448855

julia> pFq((1, 2), (4, ), 1) # a well-poised ₂F₁
2.9999999999999996

julia> pFq((1, 2+im), (3.5, ), exp(im*π/3)) # ₂F₁ at that special point in ℂ
0.6786952632946592 + 0.4523504929285015im

julia> pFq((1, 2+im), (3.5, ), exp(im*big(π)/3)) # More digits, you say?
0.6786952632946589823300834090168381068073515492901393549193461972311801512528996 + 0.4523504929285013648194489713901658143893464679689810112119412310631860619947939im

julia> pFq((1, 2+im, 2.5), (3.5, 4), exp(im*π/3)) # ₃F₂ because why not
0.8434434031615691 + 0.34175507615463174im

julia> pFq((1, 2+im, 2.5), (3.5, 4), exp(im*big(π)/3)) # Also in extended precision
0.8434434031615690763389963048175253868863156451003855955719081209861492349268002 + 0.3417550761546319732614495656712509723030350666571102474299311122586948108413206im

julia> pFq((1, 1), (), -1) # A divergent series
0.5963473623231942

julia> pFq((1, 1), (), -big(1))
0.5963473623231940743410784993692793760741778601525487815734849104823272191142015
```
"""
pFq

pFq(::Tuple{}, ::Tuple{}, z; kwds...) = exp(z)
pFq(α::NTuple{1}, ::Tuple{}, z; kwds...) = isinteger(α[1]) ? (1-z)^-α[1] : exp(-α[1]*log1p(-z))
pFq(α::NTuple{1, <: Integer}, ::Tuple{}, z; kwds...) = (1-z)^-α[1]
pFq(α::NTuple{1}, β::NTuple{1}, z; kwds...) = _₁F₁(α[1], β[1], z; kwds...)
pFq(α::NTuple{2, Any}, β::NTuple{1}, z; kwds...) = _₂F₁(α[1], α[2], β[1], z; kwds...)

function pFq(α::NTuple{p, Any}, β::NTuple{q, Any}, z; kwds...) where {p, q}
    T1 = isempty(α) ? Any : mapreduce(typeof, promote_type, α)
    T2 = isempty(β) ? Any : mapreduce(typeof, promote_type, β)
    pFq(T1.(α), T2.(β), z; kwds...)
end
function pFq(α::NTuple{p, T1}, β::Tuple{}, z; kwds...) where {p, T1}
    z = float(z)
    T = promote_type(T1, typeof(z))
    if isterminating(α)
        return pFqmaclaurin(α, β, z, getterminatinglength(α); kwds...)::T
    else
        return pFqweniger(α, β, z; kwds...)::T
    end
end
function pFq(α::Tuple{}, β::NTuple{q, T2}, z; kwds...) where {q, T2}
    z = float(z)
    T = promote_type(T2, typeof(z))
    if real(z) > 0
        return pFqmaclaurin(α, β, z; kwds...)::T
    else
        return pFqweniger(α, β, z; kwds...)::T
    end
end
function pFq(α::NTuple{p, T1}, β::NTuple{q, T2}, z; kwds...) where {p, q, T1, T2}
    z = float(z)
    T = promote_type(T1, T2, typeof(z))
    if isacommonparameter(α, β)
        return pFq(cancelacommonparameter(α, β)..., z; kwds...)::T
    end
    if isterminating(α)
        return pFqmaclaurin(α, β, z, getterminatinglength(α); kwds...)::T
    end
    if p ≤ q
        if real(z) > 0
            return pFqmaclaurin(α, β, z; kwds...)::T
        else
            return pFqweniger(α, β, z; kwds...)::T
        end
    elseif p == q + 1
        if ispositivelybalanced(α, β)
            #return pFqconformalrationalgeneralizedeuler(α, β, z; kwds...)::T
            return pFqconformalrational(α, β, z; kwds...)::T
        else
            return pFqconformalrational(α, β, z; kwds...)::T
        end
    else
        return pFqweniger(α, β, z; kwds...)::T
    end
end

pFq(α::AbstractVector, β::AbstractVector, z; kwds...) = pFq(Tuple(α), Tuple(β), z; kwds...)

"""
Compute the generalized hypergeometric function `₃F₂(a₁, 1, 1, b₁, 2, z)`.
"""
_₃F₂(a₁, b₁, z; kwds...) = _₃F₂(a₁, 1, 1, b₁, 2, z; kwds...)
"""
Compute the generalized hypergeometric function `₃F₂(a₁, a₂, a₃, b₁, b₂; z)`.
"""
_₃F₂(a₁, a₂, a₃, b₁, b₂, z; kwds...) = pFq((a₁, a₂, a₃), (b₁, b₂), z; kwds...)
