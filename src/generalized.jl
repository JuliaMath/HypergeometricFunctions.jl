@doc raw"""
    pFq(α, β; z)
Compute the generalized hypergeometric function, defined by
```math
\operatorname{pFq}(α, β; z) = \sum_{k=0}^\infty \dfrac{(\alpha_1)_k\cdots(\alpha_p)_k}{(\beta_1)_k\cdots(\beta_q)_k}\dfrac{z^k}{k!}.
```
"""
pFq

pFq(::Tuple{}, ::Tuple{}, z; kwds...) = exp(z)
pFq(α::NTuple{1}, ::Tuple{}, z; kwds...) = exp(-α[1]*log1p(-z))
pFq(α::NTuple{1}, β::NTuple{1}, z; kwds...) = _₁F₁(α[1], β[1], float(z); kwds...)
pFq(α::NTuple{2, Any}, β::NTuple{1}, z; kwds...) = _₂F₁(α[1], α[2], β[1], float(z); kwds...)

function pFq(α::NTuple{p, Any}, β::NTuple{q, Any}, z; kwds...) where {p, q}
    if p ≤ q
        if real(z) ≥ 0
            return pFqmaclaurin(α, β, float(z); kwds...)
        else
            return pFqweniger(α, β, float(z); kwds...)
        end
    elseif p == q + 1
        if abs(z) ≤ ρ
            return pFqmaclaurin(α, β, float(z); kwds...)
        else
            return pFqweniger(α, β, float(z); kwds...)
        end
    else
        return pFqweniger(α, β, float(z); kwds...)
    end
end

pFq(α::AbstractVector, β::AbstractVector, z; kwds...) = pFq(Tuple(α), Tuple(β), z; kwds...)

"""
Compute the generalized hypergeometric function `pFq(α, β; z)` by continued fraction.
"""
function pFqcontinuedfraction(α::AbstractVector{S}, β::AbstractVector{U}, z::V) where {S, U, V}
    T = promote_type(S, U, V)
    numerator(i) = - z * prod(i .+ α) / prod(i .+ β) / (i + 1)
    denominator(i) = 1 - numerator(i)
    K = continuedfraction(denominator, numerator, 10eps(real(T))) - denominator(0)
    # iszero(K + 1) leads to Inf, when e.g. Float64s run out of digits
    return 1 + z * prod(α) / prod(β) / (1 + K)
end

"""
Compute the generalized hypergeometric function `₃F₂(a₁, a₂, a₃, b₁, b₂; z)`.
"""
_₃F₂(a₁, a₂, a₃, b₁, b₂, z; kwds...) = pFq((a₁, a₂, a₃), (b₁, b₂), z; kwds...)
_₃F₂(a₁, b₁, z; kwds...) = _₃F₂(a₁, 1, 1, b₁, 2, z; kwds...)
