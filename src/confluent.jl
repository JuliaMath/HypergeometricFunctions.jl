# The references to special cases are to NIST's DLMF.

"""
Compute Kummer's confluent hypergeometric function `₁F₁(a, b, z)`.
"""
function _₁F₁(a, b, z; kwds...)
    z = float(z)
    if isequal(a, b) # 13.6.1
        return exp(z)
    elseif -b ∈ ℕ
        if -a ∈ ℕ₀ && real(a) ≥ real(b)
            return _₁F₁maclaurin(a, b, z; kwds...)
        else
            return throw(DomainError(b, "M(a, b, z) = ₁F₁(a, b, z) is not defined for negative integer b unless a is a nonpositive integer and a ≥ b."))
        end
    elseif -a ∈ ℕ₀
        return _₁F₁maclaurin(a, b, z; kwds...)
    elseif isequal(a, 1)
        if isequal(b, 2) # 13.6.2
            return iszero(z) ? one(z) : expm1(z)/z
        end
    end
    return _₁F₁general(a, b, z; kwds...)
end

function _₁F₁general(a, b, z; kwds...)
    if real(z) > 0
        return _₁F₁maclaurin(a, b, z; kwds...)
    else
        return pFqweniger((a, ), (b, ), z; kwds...)
    end
end

"""
Compute Kummer's confluent hypergeometric function `M(a, b, z) = ₁F₁(a, b, z)`.
"""
const M = _₁F₁

"""
Compute Tricomi's confluent hypergeometric function `U(a, b, z) ∼ z⁻ᵃ ₂F₀((a, a-b+1), (), -z⁻¹)`.
"""
function U(a, b, z; kwds...)
    return z^-a*pFq((a, a-b+1), (), -inv(z); kwds...)
end
