@doc raw"""
    pFq(α, β; z)
Compute the generalized hypergeometric function, defined by
```math
\operatorname{pFq}(α, β; z) = \sum_{k=0}^\infty \dfrac{(\alpha_1)_k\cdots(\alpha_p)_k}{(\beta_1)_k\cdots(\beta_q)_k}\dfrac{z^k}{k!}.
```
"""
function pFq(α::AbstractVector, β::AbstractVector, z; kwds...)
    if length(α) == length(β) == 0
        return exp(z)
    elseif length(α) == 1 && length(β) == 0
        return exp(-α[1]*log1p(-z))
    elseif length(α) == 1 && length(β) == 1
        return _₁F₁(α[1], β[1], float(z); kwds...)
    elseif length(α) == 2 && length(β) == 1
        return _₂F₁(α[1], α[2], β[1], float(z); kwds...)
    elseif length(α) ≤ length(β)
        if real(z) ≥ 0
            return pFqmaclaurin(α, β, float(z); kwds...)
        else
            return pFqweniger(α, β, float(z); kwds...)
        end
    elseif length(α) == length(β) + 1
        if abs(z) ≤ ρ
            return pFqmaclaurin(α, β, float(z); kwds...)
        else
            return pFqweniger(α, β, float(z); kwds...)
        end
    else
        return pFqweniger(α, β, float(z); kwds...)
    end
end

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
function _₃F₂(a₁, a₂, a₃, b₁, b₂, z; kwds...)
    if abs(z) ≤ ρ
        _₃F₂maclaurin(a₁, a₂, a₃, b₁, b₂, float(z); kwds...)
    else
        pFqweniger([a₁, a₂, a₃], [b₁, b₂], float(z); kwds...)
    end
end
_₃F₂(a₁, b₁, z; kwds...) = _₃F₂(1, 1, a₁, 2, b₁, z; kwds...)
