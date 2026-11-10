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
    z = float(z)
    if z isa Real && z < 0
        z = complex(z)
    end

    if isequal(a, 1) && isequal(b, 1)
        return _U11(z)
    elseif -a ∈ ℕ₀
        return _U_polynomial(-round(Int, real(a)), b, z)
    elseif b ∈ ℤ
        n = round(Int, real(b))
        if n > 0
            return _U_integer_b(a, n - 1, z)
        else
            return exp((1 - b) * log(z)) * _U_integer_b(a - b + 1, 1 - n, z)
        end
    else
        return gamma(1 - b) * ogamma(a - b + 1) * _₁F₁(a, b, z; kwds...) +
               gamma(b - 1) * ogamma(a) * exp((1 - b) * log(z)) * _₁F₁(a - b + 1, 2 - b, z; kwds...)
    end
end

_U11(z::Float32) = Float32(_U11(Float64(z)))
_U11(z::ComplexF32) = ComplexF32(_U11(ComplexF64(z)))
_U11(z) = exp(z) * expint(z)

function _U_polynomial(m::Int, b, z)
    T = promote_type(typeof(b), typeof(z))
    ret = zero(T)
    for s = 0:m
        ret += (-1)^s * binomial(m, s) * pochhammer(b + s, m - s) * z^s
    end
    return ret
end

function _U_integer_b(a, n::Int, z)
    T = promote_type(typeof(a), typeof(z))
    logz = log(z)
    sum1 = zero(T)
    term = one(T)
    tol = 10eps(real(T))
    for k = 0:KMAX
        coeff = logz + digamma(a + k) - digamma(k + 1) - digamma(n + k + 1)
        addend = term * coeff
        sum1 += addend
        if k > 0 && !errcheck(sum1, sum1 - addend, tol)
            break
        end
        term *= (a + k) * z / ((n + k + 1) * (k + 1))
    end

    sum2 = zero(T)
    for k = 1:n
        sum2 += factorial(k - 1) * pochhammer(1 - a + k, n - k) * z^(-k) / factorial(n - k)
    end

    return (-1)^(n + 1) * ogamma(a - n) / factorial(n) * sum1 + ogamma(a) * sum2
end
