"""
Compute Kummer's confluent hypergeometric function `M(a, b, z) = ₁F₁(a; b; z)`.
"""
function _₁F₁(a, b, z)
    if real(z) ≥ 0
        return _₁F₁maclaurin(a, b, z)
    else
        return exp(z)*_₁F₁(b-a, b, -z)
    end
end

"""
Compute Kummer's confluent hypergeometric function `M(a, b, z) = ₁F₁(a; b; z)`.
"""
const M = _₁F₁

"""
Compute Tricomi's confluent hypergeometric function `U(a, b, z) ∼ z⁻ᵃ ₂F₀([a, a-b+1]; []; -z⁻¹)`.
"""
function U(a, b, z)
    return z^-a*pFq([a, a-b+1], typeof(z)[], -inv(z))
end
