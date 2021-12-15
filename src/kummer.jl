"""
Compute Kummer's confluent hypergeometric function `₁F₁(a; b; z)`.
"""
function _₁F₁(a, b, z)
    if real(z) ≥ 0
        return _₁F₁maclaurin(a, b, z)
    else
        return exp(z)*_₁F₁(b-a, b, -z)
    end
end

const M = _₁F₁
