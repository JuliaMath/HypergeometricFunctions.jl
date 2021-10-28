# The references to special cases are to Table of Integrals, Series, and Products, § 9.121, followed by NIST's DLMF.

"""
Compute the Gauss hypergeometric function `₂F₁(a, b; c; z)`.
"""
function _₂F₁(a, b, c, z)
    if real(b) < real(a)
        return _₂F₁(b, a, c, z) # ensure a ≤ b
    elseif isequal(a, c) # 1. 15.4.6
        return exp(-b*log1p(-z))
    elseif isequal(b, c) # 1. 15.4.6
        return exp(-a*log1p(-z))
    elseif isequal(c, 0.5)
        if isequal(a+b, 0) # 31. 15.4.11 & 15.4.12
            return cosnasinsqrt(2b, z)
        elseif isequal(a+b, 1) # 32. 15.4.13 & 15.4.14
            return cosnasinsqrt(1-2b, z)*exp(-log1p(-z)/2)
        elseif isequal(b-a, 0.5) # 15.4.7 & 15.4.8
            return expnlog1pcoshatanhsqrt(-2a, z)
        end
    elseif isequal(c, 1.5)
        if abeqcd(a, b, 0.5) # 13. 15.4.4 & 15.4.5
            return sqrtasinsqrt(z)
        elseif abeqcd(a, b, 1) # 14.
            return sqrtasinsqrt(z)*exp(-log1p(-z)/2)
        elseif abeqcd(a, b, 0.5, 1) # 15. 15.4.2 & 15.4.3
            return sqrtatanhsqrt(z)
        elseif isequal(a+b, 1) # 29. 15.4.15 & 15.4.16
            return sinnasinsqrt(1-2b, z)
        elseif isequal(a+b, 2) # 30.
            return sinnasinsqrt(2-2b, z)*exp(-log1p(-z)/2)
        elseif isequal(b-a, 0.5) # 4. 15.4.9 & 15.4.10
            return expnlog1psinhatanhsqrt(1-2a, z)
        end
    elseif isequal(c, 2)
        if abeqcd(a, b, 1) # 6. 15.4.1
            return log1pover(-z)
        elseif a ∈ ℤ && b == 1 # 5.
            return expm1nlog1p(1-a, -z)
        elseif a == 1 && b ∈ ℤ # 5.
            return expm1nlog1p(1-b, -z)
        end
    elseif isequal(c, 4)
        if abeqcd(a, b, 2)
            return logandpoly(z)
        end
    elseif isequal(c, 2.5) && abeqcd(a, b, 1, 1.5)
        return speciallog(z)
    end
    return _₂F₁general(a, b, c, z) # catch-all
end

# Special case of (-x)^a*_₂F₁ to handle LogNumber correctly in RiemannHilbert.jl
function mxa_₂F₁(a, b, c, z)
    if isequal(c, 2)
        if abeqcd(a, b, 1) # 6. 15.4.1
            return log1p(-z)
        end
    elseif isequal(c, 4)
        if abeqcd(a, b, 2)
            return 6*(-2 + (1-2/z)*log1p(-z))
        end
    end
    return (-z)^a*_₂F₁(a, b, c, z)
end

"""
Compute the Gauss hypergeometric function `₂F₁(a, b; c; z)` with general parameters a, b, and c.
This polyalgorithm is designed based on the paper

N. Michel and M. V. Stoitsov, Fast computation of the Gauss hypergeometric function with all its parameters complex with application to the Pöschl–Teller–Ginocchio potential wave functions, Comp. Phys. Commun., 178:535–551, 2008.
"""
function _₂F₁general(a, b, c, z)
    T = promote_type(typeof(a), typeof(b), typeof(c), typeof(z))

    real(b) < real(a) && return _₂F₁general(b, a, c, z)
    real(c) < real(a) + real(b) && return exp((c-a-b)*log1p(-z))*_₂F₁general(c-a, c-b, c, z)

    if abs(z) ≤ ρ || -a ∈ ℕ₀ || -b ∈ ℕ₀
        return _₂F₁maclaurin(a, b, c, z)
    elseif abs(z/(z-1)) ≤ ρ
        return exp(-a*log1p(-z))_₂F₁maclaurin(a, c-b, c, z/(z-1))
    elseif abs(inv(z)) ≤ ρ
        return _₂F₁Inf(a, b, c, z)
    elseif abs(1-inv(z)) ≤ ρ
        return exp(-a*log1p(-z))*_₂F₁Inf(a, c-b, c, z/(z-1))
    elseif abs(1-z) ≤ ρ
        return _₂F₁one(a, b, c, z)
    elseif abs(inv(1-z)) ≤ ρ
        return exp(-a*log1p(-z))*_₂F₁one(a, c-b, c, z/(z-1))
    else
        return pFqweniger([a, b], [c], z) # _₂F₁taylor(a, b, c, z)
    end
end

"""
Compute the Gauss hypergeometric function `₂F₁(a, b; c; z)` with general parameters a, b, and c.
This polyalgorithm is designed based on the review

J. W. Pearson, S. Olver and M. A. Porter, Numerical methos for the computation of the confluent and Gauss hypergeometric functions, arXiv:1407.7786, 2014.
"""
function _₂F₁general2(a, b, c, z)
    T = promote_type(typeof(a), typeof(b), typeof(c), typeof(z))
    if abs(z) ≤ ρ || -a ∈ ℕ₀ || -b ∈ ℕ₀
        return _₂F₁maclaurin(a, b, c, z)
    elseif abs(z / (z - 1)) ≤ ρ && absarg(1 - z) < convert(real(T), π) # 15.8.1
        w = z/(z-1)
        return _₂F₁maclaurin(a, c-b, c, w)*exp(-a*log1p(-z))
    elseif abs(inv(z)) ≤ ρ && absarg(-z) < convert(real(T), π)
        w = inv(z)
        if isapprox(a, b) # 15.8.8
            return gamma(c)/gamma(a)/gamma(c-a)*(-w)^a*_₂F₁logsumalt(a, c-a, z, w)
        elseif a-b ∉ ℤ # 15.8.2
            return gamma(c)*((-w)^a*gamma(b-a)/gamma(b)/gamma(c-a)*_₂F₁maclaurin(a, a-c+1, a-b+1, w)+(-w)^b*gamma(a-b)/gamma(a)/gamma(c-b)*_₂F₁maclaurin(b, b-c+1, b-a+1, w))
        end # TODO: full 15.8.8
    elseif abs(inv(1-z)) ≤ ρ && absarg(-z) < convert(real(T), π)
        w = inv(1-z)
        if isapprox(a, b) # 15.8.9
            return gamma(c)*exp(-a*log1p(-z))/gamma(a)/gamma(c-b)*_₂F₁logsum(a, c-a, z, w, 1)
        elseif a-b ∉ ℤ # 15.8.3
            return gamma(c)*(exp(-a*log1p(-z))*gamma(b-a)/gamma(b)/gamma(c-a)*_₂F₁maclaurin(a, c-b, a-b+1, w)+exp(-b*log1p(-z))*gamma(a-b)/gamma(a)/gamma(c-b)*_₂F₁maclaurin(b, c-a, b-a+1, w))
        end # TODO: full 15.8.9
    elseif abs(1-z) ≤ ρ && absarg(z) < convert(real(T), π) && absarg(1-z) < convert(real(T), π)
        w = 1-z
        if isapprox(c, a + b) # 15.8.10
            return gamma(c)/gamma(a)/gamma(b)*_₂F₁logsum(a, b, z, w, -1)
        elseif c - a - b ∉ ℤ # 15.8.4
            return gamma(c)*(gamma(c-a-b)/gamma(c-a)/gamma(c-b)*_₂F₁maclaurin(a, b, a+b-c+1, w)+exp((c-a-b)*log1p(-z))*gamma(a+b-c)/gamma(a)/gamma(b)*_₂F₁maclaurin(c-a, c-b, c-a-b+1, w))
        end # TODO: full 15.8.10
    elseif abs(1-inv(z)) ≤ ρ && absarg(z) < convert(real(T), π) && absarg(1-z) < convert(real(T), π)
        w = 1-inv(z)
        if isapprox(c, a + b) # 15.8.11
            return gamma(c)*z^(-a)/gamma(a)*_₂F₁logsumalt(a, b, z, w)
        elseif c - a - b ∉ ℤ # 15.8.5
            return gamma(c)*(z^(-a)*gamma(c-a-b)/gamma(c-a)/gamma(c-b)*_₂F₁maclaurin(a, a-c+1, a+b-c+1, w)+z^(a-c)*(1-z)^(c-a-b)*gamma(a+b-c)/gamma(a)/gamma(b)*_₂F₁maclaurin(c-a, 1-a, c-a-b+1, w))
        end # TODO: full 15.8.11
    elseif abs(z-0.5) > 0.5
        if isapprox(a, b) && !isapprox(c, a+0.5)
            return gamma(c)/gamma(a)/gamma(c-a)*(0.5-z)^(-a)*_₂F₁continuationalt(a, c, 0.5, z)
        elseif a-b ∉ ℤ
            return gamma(c)*(gamma(b-a)/gamma(b)/gamma(c-a)*(0.5-z)^(-a)*_₂F₁continuation(a, a+b, c, 0.5, z) + gamma(a-b)/gamma(a)/gamma(c-b)*(0.5-z)^(-b)*_₂F₁continuation(b, a+b, c, 0.5, z))
        end
    end
    return pFqweniger([a, b], [c], z)
end
