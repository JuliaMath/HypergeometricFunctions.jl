@inline errcheck(x, y, tol) = isfinite(x) && isfinite(y) && (norm(x-y) > max(norm(x), norm(y))*tol)

kind2string(::Val{0}) = "₀"
kind2string(::Val{1}) = "₁"
kind2string(::Val{2}) = "₂"
kind2string(::Val{3}) = "₃"
kind2string(::Val{4}) = "₄"
kind2string(::Val{5}) = "₅"
kind2string(::Val{6}) = "₆"
kind2string(::Val{7}) = "₇"
kind2string(::Val{8}) = "₈"
kind2string(::Val{9}) = "₉"
kind2string(::Val{p}) where p = kind2string(Val(div(p, 10)))*kind2string(Val(rem(p, 10)))
pFq2string(::Val{p}, ::Val{q}) where {p, q} = kind2string(Val(p))*"F"*kind2string(Val(q))

# Same as in FastTransforms.jl
"""
Pochhammer symbol ``(x)_n = \\frac{\\Gamma(x+n)}{\\Gamma(x)}`` for the rising factorial.
"""
function pochhammer(x::T, n::Integer) where T
    S = typeof(zero(T)/one(T))
    ret = one(S)
    if n ≥ 0
        for i = 0:n-1
            ret *= x+i
        end
    else
        ret /= pochhammer(x+n, -n)
    end
    ret
end

pochhammer(x::T, n::Number) where T = isinteger(n) ? pochhammer(x, Int(n)) : ogamma(x)/ogamma(x+n)

function pochhammer(x::Number, n::UnitRange{T}) where T <: Real
    ret = Vector{promote_type(typeof(x),T)}(undef, length(n))
    ret[1] = pochhammer(x, first(n))
    for i = 2:length(n)
        ret[i] = (x+n[i]-1)*ret[i-1]
    end
    ret
end

ogamma(x::Number) = (isinteger(x) && x<0) ? zero(float(x)) : inv(gamma(x))

unsafe_gamma(z) = gamma(z)
function unsafe_gamma(x::Real)
    try
        gamma(x)
    catch err
        err isa DomainError || rethrow()
        x == 0 && return inv(zero(float(x)))
        oftype(float(x), NaN)
    end
end

"""
    @clenshaw(x, c...)

Evaluate the Chebyshev polynomial series ``\\displaystyle \\sum_{k=1}^N c[k] T_{k-1}(x)`` by the Clenshaw algorithm.

External links: [DLMF](https://dlmf.nist.gov/3.11#ii), [Wikipedia](https://en.wikipedia.org/wiki/Clenshaw_algorithm).

# Examples
```jldoctest
julia> HypergeometricFunctions.@clenshaw(1, 1, 2, 3)
6

julia> HypergeometricFunctions.@clenshaw(0.5, 1, 2, 3)
0.5
```
"""
macro clenshaw(x, c...)
    a, b = :(zero(t)), :(zero(t))
    as = []
    for k = length(c):-1:2
        ak = Symbol("a", k)
        push!(as, :($ak = $a))
        a = :(muladd(t, $a, $(esc(c[k]))-$b))
        b = :($ak)
    end
    ex = Expr(:block, as..., :(muladd($(esc(x)), $a, $(esc(c[1]))-$b)))
    Expr(:block, :(t = $(esc(2))*$(esc(x))), ex)
end


const ρ = 0.72
const ρϵ = 0.71

struct ℕ end

Base.in(n::Integer, ::Type{ℕ}) = n > 0
Base.in(n::Number, ::Type{ℕ}) = isinteger(n) && round(Int, n) ∈ ℕ

struct ℕ₀ end

Base.in(n::Integer, ::Type{ℕ₀}) = n ≥ 0
Base.in(n::Number, ::Type{ℕ₀}) = isinteger(n) && round(Int, n) ∈ ℕ₀

struct ℤ end

Base.in(n::Number, ::Type{ℤ}) = isinteger(n)

abeqcd(a, b, cd) = isequal(a, b) && isequal(b, cd)
abeqcd(a, b, c, d) = isequal(a, c) && isequal(b, d)

iswellpoised(a, b, c) = real(c - a - b) > 0
isalmostwellpoised(a, b, c) = real(c - a - b) ≥ 0

absarg(z) = abs(angle(z))

sqrtatanhsqrt(x) = x == 0 ? one(x) : (s = sqrt(-x); atan(s)/s)
sqrtasinsqrt(x) = x == 0 ? one(x) : (s = sqrt(x); asin(s)/s)
sinnasinsqrt(n, x) = x == 0 ? one(x) : (s = sqrt(x); sin(n*asin(s))/(n*s))
cosnasinsqrt(n, x) = cos(n*asin(sqrt(x)))
expnlog1pcoshatanhsqrt(n, x) = x == 0 ? one(x) : (s = sqrt(x); (exp(n*log1p(s))+exp(n*log1p(-s)))/2)
expnlog1psinhatanhsqrt(n, x) = x == 0 ? one(x) : (s = sqrt(x); (exp(n*log1p(s))-exp(n*log1p(-s)))/(2n*s))

sqrtatanhsqrt(x::Real) = x == 0 ? one(x) : x > 0 ? (s = sqrt(x); atanh(s)/s) : (s = sqrt(-x); atan(s)/s)
sqrtasinsqrt(x::Real) = x == 0 ? one(x) : x > 0 ? (s = sqrt(x); asin(s)/s) : (s = sqrt(-x); asinh(s)/s)
sinnasinsqrt(n, x::Real) = x == 0 ? one(x) : x > 0 ? (s = sqrt(x); sin(n*asin(s))/(n*s)) : (s = sqrt(-x); sinh(n*asinh(s))/(n*s))
cosnasinsqrt(n, x::Real) = x > 0 ? cos(n*asin(sqrt(x))) : cosh(n*asinh(sqrt(-x)))
expnlog1pcoshatanhsqrt(n, x::Real) = x == 0 ? one(x) : x > 0 ? exp(n/2*log1p(-x))*cosh(n*atanh(sqrt(x))) : exp(n/2*log1p(-x))*cos(n*atan(sqrt(-x)))
expnlog1psinhatanhsqrt(n, x::Real) = x == 0 ? one(x) : x > 0 ? (s = sqrt(x); exp(n/2*log1p(-x))*sinh(n*atanh(s))/(n*s)) : (s = sqrt(-x); exp(n/2*log1p(-x))*sin(n*atan(s))/(n*s))

expm1nlog1p(n, x) = x == 0 ? one(x) : expm1(n*log1p(x))/(n*x)

log1pover(x) = iszero(x) ? one(x) : log1p(x)/x
logandpoly(x) = iszero(x) ? one(x) : 6*(-2x+(x-2)*log1p(-x))/x^3
function logandpoly(x::Union{Float64, ComplexF64})
    if abs(x) > 0.2
        6*(-2x+(x-2)*log1p(-x))/x^3
    else
        logandpolyseries(x)
    end
end

logandpolyseries(x::Union{Float64, ComplexF64}) = @evalpoly(x, 1.0, 1.0, 0.9, 0.8, 0.7142857142857143, 0.6428571428571429, 0.5833333333333334, 0.5333333333333333, 0.4909090909090909, 0.45454545454545453, 0.4230769230769231, 0.3956043956043956, 0.37142857142857144, 0.35, 0.33088235294117646, 0.3137254901960784, 0.2982456140350877, 0.28421052631578947, 0.2714285714285714, 0.2597402597402597)

speciallog(x::Real) = iszero(x) ? one(x) : (x > 0 ? (s = sqrt(x); 3(atanh(s)-s)/s^3) : (s = sqrt(-x); 3(s-atan(s))/s^3))
speciallog(x) = iszero(x) ? one(x) : (s = sqrt(-x); 3(s-atan(s))/s^3)
function speciallog(x::Float64)
    if x > 0.2
        s = sqrt(x)
        3(atanh(s)-s)/s^3
    elseif x < -0.2
        s = sqrt(-x)
        3(s-atan(s))/s^3
    else
        speciallogseries(x)
    end
end
function speciallog(x::ComplexF64)
    if abs(x) > 0.2
        s = sqrt(-x)
        3(s-atan(s))/s^3
    else
        speciallogseries(x)
    end
end
# The Maclaurin series fails to be accurate to 1e-15 near x ≈ ±0.2. So we use a highly accurate Chebyshev expansion.
speciallogseries(x::Float64) = @clenshaw(5.0x, 1.0087391788544393911192, 1.220474262857857637288e-01, 8.7957928919918696061703e-03, 6.9050958578444820505037e-04, 5.7037120050065804396306e-05, 4.8731405131379353370205e-06, 4.2648797509486828820613e-07, 3.800372208946157617901e-08, 3.434168059359993493634e-09, 3.1381484326392473547608e-10, 2.8939845618385022798906e-11, 2.6892186934806386106143e-12, 2.5150879096374730760324e-13, 2.3652490233687788117887e-14, 2.2349973917002118259929e-15, 2.120769988408948118084e-16)
speciallogseries(x::ComplexF64) = @evalpoly(x, 1.0000000000000000000000, 5.9999999999999999999966e-01, 4.2857142857142857142869e-01, 3.3333333333333333333347e-01, 2.7272727272727272727292e-01, 2.3076923076923076923072e-01, 1.9999999999999999999996e-01, 1.7647058823529411764702e-01, 1.5789473684210526315786e-01, 1.4285714285714285714283e-01, 1.3043478260869565217384e-01, 1.2000000000000000000000e-01, 1.1111111111111111111109e-01, 1.0344827586206896551722e-01, 9.6774193548387096774217e-02, 9.0909090909090909090938e-02, 8.5714285714285714285696e-02, 8.1081081081081081081064e-02, 7.6923076923076923076907e-02, 7.3170731707317073170688e-02)


"""
    @lanczosratio(z, ϵ, c₀, c...)

Evaluate ``\\dfrac{\\displaystyle \\sum_{k=0}^{N-1} \\frac{c[k+1]}{(z+k)(z+k+\\epsilon)}}{\\displaystyle c_0 + \\sum_{k=0}^{N-1} \\frac{c[k+1]}{z+k}}``.

This ratio is used in the Lanczos approximation of ``\\log\\frac{\\Gamma(z+\\epsilon)}{\\Gamma(z)}`` in

> N. Michel and M. V. Stoitsov, [Fast computation of the Gauss hypergeometric function with all its parameters complex with application to the Pöschl–Teller–Ginocchio potential wave functions](https://doi.org/10.1016/j.cpc.2007.11.007), *Comp. Phys. Commun.*, **178**:535–551, 2008.
"""
macro lanczosratio(z, ϵ, c₀, c...)
    ex_num = :(zero(z))
    ex_den = esc(c₀)
    for k = 0:length(c)-1
        temp = :(inv(z+$k))
        ex_num = :(muladd($(esc(c[k+1])), $temp/(zpϵ+$k), $ex_num))
        ex_den = :(muladd($(esc(c[k+1])), $temp, $ex_den))
    end
    ex = :($ex_num/$ex_den)
    Expr(:block, :(z = $(esc(z))), :(zpϵ = $(esc(z))+$(esc(ϵ))), ex)
end

lanczosratio(z::Union{Float64, ComplexF64}, ϵ::Union{Float64, ComplexF64}) = @lanczosratio(z, ϵ, 0.99999999999999709182, 57.156235665862923517, -59.597960355475491248, 14.136097974741747174, -0.49191381609762019978, 0.33994649984811888699E-4, 0.46523628927048575665E-4, -0.98374475304879564677E-4, 0.15808870322491248884E-3, -0.21026444172410488319E-3, 0.21743961811521264320E-3, -0.16431810653676389022E-3, 0.84418223983852743293E-4, -0.26190838401581408670E-4, 0.36899182659531622704E-5)

function lanczosapprox(z::Union{Float64, ComplexF64}, ϵ::Union{Float64, ComplexF64})
    zm0p5 = z-0.5
    zpgm0p5 = zm0p5+4.7421875
    return zm0p5*log1p(ϵ/zpgm0p5) + ϵ*log(zpgm0p5+ϵ) - ϵ + log1p(-ϵ*lanczosratio(z, ϵ))
end

function H(z::Union{Float64, ComplexF64}, ϵ::Union{Float64, ComplexF64})
    zm0p5 = z-0.5
    zpgm0p5 = zm0p5+4.7421875
    if real(z) ≥ 1/2
        if z == z+ϵ # ϵ is numerical 0
            zm0p5/zpgm0p5 + log(zpgm0p5) - 1 - lanczosratio(z, ϵ)
        else
            expm1( zm0p5*log1p(ϵ/zpgm0p5) + ϵ*log(zpgm0p5+ϵ) - ϵ + log1p(-ϵ*lanczosratio(z, ϵ)) )/ϵ
        end
    else
        tpz = sinpi(z)/cospi(z)
        if z == z+ϵ # ϵ is numerical 0
            H(1-z, ϵ) - π/tpz
        else
            temp = (cospi(ϵ) + sinpi(ϵ)/tpz)*H(1-z, -ϵ) + .5ϵ*(π*sinc(.5ϵ))^2 - π*sinc(ϵ)/tpz
            temp/(1-ϵ*temp)
        end
    end
end

"""
Compute the function ``\\dfrac{\\frac{1}{\\Gamma(z)}-\\frac{1}{\\Gamma(z+\\epsilon)}}{\\epsilon}`` by the method dscribed in

> N. Michel and M. V. Stoitsov, [Fast computation of the Gauss hypergeometric function with all its parameters complex with application to the Pöschl–Teller–Ginocchio potential wave functions](https://doi.org/10.1016/j.cpc.2007.11.007), *Comp. Phys. Commun.*, **178**:535–551, 2008.
"""
function G(z::Union{Float64, ComplexF64}, ϵ::Union{Float64, ComplexF64})
    n, zpϵ = round(Int, real(z)), z+ϵ
    if abs(ϵ) > 0.1
        (inv(unsafe_gamma(z))-inv(unsafe_gamma(zpϵ)))/ϵ
    elseif z ≠ zpϵ
        m = round(Int, real(zpϵ))
        if z == n && n ≤ 0
            -inv(ϵ*unsafe_gamma(zpϵ))
        elseif zpϵ == m && m ≤ 0
            inv(ϵ*unsafe_gamma(z))
        elseif abs(z+abs(n)) < abs(zpϵ+abs(m))
            H(z, ϵ)/unsafe_gamma(zpϵ)
        else
            H(zpϵ, -ϵ)/unsafe_gamma(z)
        end
    else # ϵ is numerical 0
        if z == n && n ≤ 0
            (-1)^(n+1)*unsafe_gamma(1-n)
        else
            digamma(z)/unsafe_gamma(z)
        end
    end
end

G(z::T, ϵ::T) where {T<:Number} = ϵ == 0 ? digamma(z)/unsafe_gamma(z) : (inv(unsafe_gamma(z))-inv(unsafe_gamma(z+ϵ)))/ϵ
G(z::Number, ϵ::Number) = G(promote(z, ϵ)...)

"""
Compute the function ``\\dfrac{(z+\\epsilon)_m-(z)_m}{\\epsilon}`` by the method dscribed in

> N. Michel and M. V. Stoitsov, [Fast computation of the Gauss hypergeometric function with all its parameters complex with application to the Pöschl–Teller–Ginocchio potential wave functions](https://doi.org/10.1016/j.cpc.2007.11.007), *Comp. Phys. Commun.*, **178**:535–551, 2008.
"""
function P(z::Number, ϵ::Number, m::Int)
    n₀ = -round(Int, real(z))
    if ϵ == 0
        if 0 ≤ n₀ < m
            ret1, ret2, n = one(z), zero(z), 0
            while n < m
                n == n₀ && (n+=1; continue)
                ret1 *= z+n
                ret2 += inv(z+n)
                n += 1
            end
            ret1 + pochhammer(z, m)*ret2
        else
            ret = zero(z)
            for n=0:m-1
                ret += inv(z+n)
            end
            pochhammer(z, m)*ret
        end
    else
        if 0 ≤ n₀ < m
            zpϵ, ret1, ret2, n = z+ϵ, one(z), zero(z), 0
            while n < m
                n == n₀ && (n+=1; continue)
                ret1 *= zpϵ+n
                ret2 += log1p(ϵ/(z+n))
                n += 1
            end
            ret1 + pochhammer(z, m)*expm1(ret2)/ϵ
        else
            ret = zero(z)
            for n=0:m-1
                ret += log1p(ϵ/(z+n))
            end
            pochhammer(z, m)*expm1(ret)/ϵ
        end
    end
end

E(z::Number, ϵ::Number) = iszero(ϵ) ? z : expm1(ϵ*z)/ϵ

G(z::AbstractVector{BigFloat}, ϵ::BigFloat) = BigFloat[G(zi, ϵ) for zi in z]

# Transformation formula w = 1-z

function reconeα₀(a, b, c, m::Int, ϵ)
    _a, _b, _c, _ϵ = promote(a, b, c, ϵ)
    return _reconeα₀(_a, _b, _c, m, _ϵ)
end
function _reconeα₀(a::T, b::T, c::T, m::Int, ϵ::T) where {T}
    if ϵ == 0
        return (-1)^m*gamma(real(T)(m))*gamma(c)/(gamma(a+m)*gamma(b+m))
    else
        return gamma(c)/(ϵ*gamma(1-m-ϵ)*gamma(a+m+ϵ)*gamma(b+m+ϵ))
    end
end
function reconeβ₀(a, b, c, w, m::Int, ϵ)
    _a, _b, _c, _, _ϵ = promote(a, b, c, real(w), ϵ)
    _w, _ = promote(w, zero(_a))
    return _reconeβ₀(_a, _b, _c, _w, m, _ϵ)
end
function _reconeβ₀(a::T, b::T, c::T, w::Number, m::Int, ϵ::T) where {T}
    if abs(ϵ) > 0.1
        return ( pochhammer(a, m)*pochhammer(b, m)/(gamma(1-ϵ)*gamma(a+m+ϵ)*gamma(b+m+ϵ)*gamma(real(T)(m)+1)) - w^ϵ/(gamma(a)*gamma(b)*gamma(m+1+ϵ)) )*gamma(c)*w^m/ϵ
    else
        return ( (G(1, -ϵ)/gamma(real(T)(m)+1)+G(m+1, ϵ))/(gamma(a+m+ϵ)*gamma(b+m+ϵ)) - (G(a+m, ϵ)/gamma(b+m+ϵ)+G(float(b)+m, ϵ)/gamma(a+m))/gamma(m+1+ϵ) - E(log(w), ϵ)/(gamma(a+m)*gamma(b+m)*gamma(m+1+ϵ)) )*gamma(c)*pochhammer(a, m)*pochhammer(b, m)*w^m
    end
end
function reconeγ₀(a, b, c, w, m::Int, ϵ)
    _a, _b, _c, _, _ϵ = promote(a, b, c, real(w), ϵ)
    _w, _ = promote(w, zero(_a))
    return _reconeγ₀(_a, _b, _c, _w, m, _ϵ)
end
_reconeγ₀(a::T, b::T, c::T, w::Number, m::Int, ϵ::T) where {T} = gamma(c)*pochhammer(a, m)*pochhammer(b, m)*w^m/(gamma(a+m+ϵ)*gamma(b+m+ϵ)*gamma(real(T)(m)+1)*gamma(1-ϵ))

function Aone(a, b, c, w, m::Int, ϵ)
    αₙ = reconeα₀(a, b, c, m, ϵ)*one(w)
    ret = m ≤ 0 ? zero(w) : αₙ
    for n = 0:m-2
        αₙ *= (a+n)*(b+n)/((n+1)*(1-m-ϵ+n))*w
        ret += αₙ
    end
    ret
end

function Bone(a, b, c, w, m::Int, ϵ)
    T = promote_type(typeof(a), typeof(b), typeof(c), typeof(w), typeof(m), typeof(ϵ))
    βₙ, γₙ = reconeβ₀(a, b, c, w, m, ϵ)*one(w), reconeγ₀(a, b, c, w, m, ϵ)*w
    ret, n = βₙ, 0
    while abs(βₙ) > 8abs(ret)*eps(real(T)) || n ≤ 0
        βₙ = (a+m+n+ϵ)*(b+m+n+ϵ)/((m+n+1+ϵ)*(n+1))*w*βₙ + ( (a+m+n)*(b+m+n)/(m+n+1) - (a+m+n) - (b+m+n) - ϵ + (a+m+n+ϵ)*(b+m+n+ϵ)/(n+1) )*γₙ/((m+n+1+ϵ)*(n+1-ϵ))
        ret += βₙ
        γₙ *= (a+m+n)*(b+m+n)/((m+n+1)*(n+1-ϵ))*w
        n += 1
    end
    ret
end

function _₂F₁one(a, b, c, z)
    m = round(Int, real(c-(a+b)))
    ϵ = c-(a+b)-m
    w = 1-z
    (-1)^m/sinc(ϵ)*(Aone(a, b, c, w, m, ϵ) + Bone(a, b, c, w, m, ϵ))
end

# Transformation formula w = 1/z

function recInfα₀(a, b, c, m::Int, ϵ)
    _a, _b, _c, _ϵ = promote(a, b, c, ϵ)
    return _recInfα₀(_a, _b, _c, m, _ϵ)
end
function _recInfα₀(a::T, b::T, c::T, m::Int, ϵ::T) where {T}
    if ϵ == 0
        return (-1)^m*gamma(real(T)(m))*gamma(c)/(gamma(a+m)*gamma(c-a))
    else
        return gamma(c)/(ϵ*gamma(1-m-ϵ)*gamma(a+m+ϵ)*gamma(c-a))
    end
end
function recInfβ₀(a, b, c, w, m::Int, ϵ)
    _a, _b, _c, _, _ϵ = promote(a, b, c, real(w), ϵ)
    _w, _ = promote(w, zero(_a))
    return _recInfβ₀(_a, _b, _c, _w, m, _ϵ)
end
function _recInfβ₀(a::T, b::T, c::T, w::Number, m::Int, ϵ::T) where {T}
    if abs(ϵ) > 0.1
        return ( pochhammer(a, m)*pochhammer(1-c+a, m)/(gamma(1-ϵ)*gamma(a+m+ϵ)*gamma(c-a)*gamma(real(T)(m)+1)) -
            (-w)^ϵ*pochhammer(1-c+a+ϵ, m)/(gamma(a)*gamma(c-a-ϵ)*gamma(m+1+ϵ)) )*gamma(c)*w^m/ϵ
    else
        return ( (pochhammer(1-c+a+ϵ, m)*G(1, -ϵ)-P(1-c+a, ϵ, m)/gamma(1-ϵ))/(gamma(c-a)*gamma(a+m+ϵ)*gamma(real(T)(m)+1)) +
            pochhammer(1-c+a+ϵ, m)*( (G(m+1, ϵ)/gamma(a+m+ϵ) - G(a+m, ϵ)/gamma(m+1+ϵ))/gamma(c-a) -
            (G(c-a, -ϵ) - E(-log(-w), -ϵ)/gamma(c-a-ϵ))/(gamma(m+1+ϵ)*gamma(a+m)) ) )*gamma(c)*pochhammer(a, m)*w^m
    end
end
function recInfγ₀(a, b, c, w, m::Int, ϵ)
    _a, _b, _c, _, _ϵ = promote(a, b, c, real(w), ϵ)
    _w, _ = promote(w, zero(_a))
    return _recInfγ₀(_a, _b, _c, _w, m, _ϵ)
end
_recInfγ₀(a::T, b::T, c::T, w::Number, m::Int, ϵ::T) where {T} = gamma(c)*pochhammer(a, m)*pochhammer(1-c+a, m)*w^m/(gamma(a+m+ϵ)*gamma(c-a)*gamma(real(T)(m)+1)*gamma(1-ϵ))

function AInf(a, b, c, w, m::Int, ϵ)
    αₙ = recInfα₀(a, b, c, m, ϵ)*one(w)
    ret = m ≤ 0 ? zero(w) : αₙ
    for n = 0:m-2
        αₙ *= (a+n)*(1-c+a+n)/((n+1)*(1-m-ϵ+n))*w
        ret += αₙ
    end
    ret
end

function BInf(a, b, c, win, m::Int, ϵ)
    T = promote_type(typeof(a), typeof(b), typeof(c), typeof(win), typeof(m), typeof(ϵ))
    w = win
    βₙ, γₙ = recInfβ₀(a, b, c, win, m, ϵ)*one(w), recInfγ₀(a, b, c, win, m, ϵ)*w
    ret, n = βₙ, 0
    while abs(βₙ) > 8abs(ret)*eps(real(T)) || n ≤ 0
        βₙ = (a+m+n+ϵ)*(1-c+a+m+n+ϵ)/((m+n+1+ϵ)*(n+1))*w*βₙ + ( (a+m+n)*(1-c+a+m+n)/(m+n+1) - (a+m+n) - (1-c+a+m+n) - ϵ + (a+m+n+ϵ)*(1-c+a+m+n+ϵ)/(n+1) )*γₙ/((m+n+1+ϵ)*(n+1-ϵ))
        ret += βₙ
        γₙ *= (a+m+n)*(1-c+a+m+n)/((m+n+1)*(n+1-ϵ))*w
        n += 1
    end
    ret
end

function _₂F₁Inf(a, b, c, z)
    m = round(Int, real(b-a))
    ϵ = b-a-m
    w = inv(z)
    (-1)^m*(-w)^a/sinc(ϵ)*(AInf(a, b, c, w, m, ϵ) + BInf(a, b, c, w, m, ϵ))
end


function _₂F₁maclaurin(a::Number, b::Number, c::Number, z::Number)
    T = promote_type(typeof(a), typeof(b), typeof(c), typeof(z))
    S₀, S₁, j = one(T), one(T)+a*b*z/c, 1
    while errcheck(S₀, S₁, 8eps(real(T)))
        rⱼ = (a+j)*z/(j+1)*(b+j)/(c+j)
        S₀, S₁ = S₁, S₁+(S₁-S₀)*rⱼ
        j += 1
    end
    return S₁
end

function _₂F₁maclaurinalt(a::Number, b::Number, c::Number, z::Number)
    T = promote_type(typeof(a), typeof(b), typeof(c), typeof(z))
    C, S, j = one(T), one(T), 0
    while abs(C) > 8abs(S)*eps(real(T)) || j ≤ 1
        C *= (a+j)/(j+1)*(b+j)/(c+j)*z
        S += C
        j += 1
    end
    return S
end

function _₂F₁continuation(s::Number, t::Number, c::Number, z₀::Number, z::Number)
    T = promote_type(typeof(s), typeof(t), typeof(c), typeof(z₀), typeof(z))
    izz₀, d0, d1 = inv(z-z₀), one(T), s/(2s-t+one(T))*((s+1)*(1-2z₀)+(t+1)*z₀-c)
    S₀, S₁, izz₀j, j = one(T), one(T)+d1*izz₀, izz₀, 2
    while errcheck(S₀, S₁, 8eps(real(T))) || j ≤ 2
        d0, d1, izz₀j = d1, (j+s-one(T))/j/(j+2s-t)*(((j+s)*(1-2z₀)+(t+1)*z₀-c)*d1 + z₀*(1-z₀)*(j+s-2)*d0), izz₀j*izz₀
        S₀, S₁ = S₁, S₁+d1*izz₀j
        j += 1
    end
    return S₁
end

function _₂F₁continuationalt(a::Number, c::Number, z₀::Number, z::Number)
    T = promote_type(typeof(a), typeof(c), typeof(z₀), typeof(z))
    izz₀ = inv(z-z₀)
    e0, e1 = one(T), (a+one(T))*(one(T)-2z₀)+(2a+one(T))*z₀-c
    f0, f1 = zero(T), one(T)-2z₀
    cⱼ = log(z₀-z)+2digamma(one(real(T)))-digamma(a)-digamma(c-a)
    S₀ = cⱼ
    cⱼ += 2/one(T)-one(T)/a
    C = a*izz₀
    S₁, j = S₀+(e1*cⱼ-f1)*C, 2
    while errcheck(S₀, S₁, 8eps(real(T))) || j ≤ 2
        f0, f1 = f1, (((j+a)*(1-2z₀)+(2a+1)*z₀-c)*f1+z₀*(1-z₀)*(j-1)*f0+(1-2z₀)*e1+2z₀*(1-z₀)*e0)/j
        e0, e1 = e1, (((j+a)*(1-2z₀)+(2a+1)*z₀-c)*e1+z₀*(1-z₀)*(j-1)*e0)/j
        C *= (a+j-1)*izz₀/j
        cⱼ += 2/T(j)-one(T)/(a+j-one(T))
        S₀, S₁ = S₁, S₁+(e1*cⱼ-f1)*C
        j += 1
    end
    return S₁
end

function _₂F₁logsum(a::Number, b::Number, z::Number, w::Number, s::Int)
    T = promote_type(typeof(a), typeof(b), typeof(z), typeof(w))
    cⱼ = 2digamma(one(real(T)))-digamma(a)-digamma(b)+s*log1p(-z)
    C, S, j = one(T), cⱼ, 0
    while abs(C) > 8abs(S)*eps(real(T)) || j ≤ 1
        C *= (a+j)/(j+1)^2*(b+j)*w
        cⱼ += 2/(j+one(T))-one(T)/(a+j)-one(T)/(b+j)
        S += C*cⱼ
        j += 1
    end
    return S
end

function _₂F₁logsumalt(a::Number, b::Number, z::Number, w::Number)
    T = promote_type(typeof(a), typeof(b), typeof(z), typeof(w))
    d, cⱼ = one(T)-b, 2digamma(one(real(T)))-digamma(a)-digamma(b)-log(-w)
    C, S, j = one(T), cⱼ, 0
    while abs(C) > 8abs(S)*eps(real(T)) || j ≤ 1
        C *= (a+j)/(j+1)^2*(d+j)*w
        cⱼ += 2/(j+one(T))-one(T)/(a+j)+one(T)/(b-(j+one(T)))
        S += C*cⱼ
        j += 1
    end
    return S
end

function _₂F₁taylor(a::Number, b::Number, c::Number, z::Number)
    T = promote_type(typeof(a), typeof(b), typeof(c), typeof(z))
    z₀ = abs(z) < 1 ? ρϵ*sign(z) : sign(z)/ρϵ
    q₀, q₁ = _₂F₁(a, b, c, z₀), a*b/c*_₂F₁(a+1, b+1, c+1, z₀)
    S₀, zz₀ = q₀, z-z₀
    S₁, zz₀j, j = S₀+q₁*zz₀, zz₀, 0
    while errcheck(S₀, S₁, 8eps(real(T))) || j ≤ 1
        q₀, q₁ = q₁, ((j*(2z₀-one(T))-c+(a+b+one(T))*z₀)*q₁ + (a+j)*(b+j)/(j+one(T))*q₀)/(z₀*(one(T)-z₀)*(j+2))
        zz₀j *= zz₀
        S₀, S₁ = S₁, S₁+q₁*zz₀j
        j += 1
    end
    return S₁
end

function _₁F₁maclaurin(a::Number, b::Number, z::Number)
    T = promote_type(typeof(a), typeof(b), typeof(z))
    S₀, S₁, j = one(T), one(T)+a*z/b, 1
    while errcheck(S₀, S₁, 8eps(real(T)))
        rⱼ = (a+j)*z/((b+j)*(j+1))
        S₀, S₁ = S₁, S₁+(S₁-S₀)*rⱼ
        j += 1
    end
    return S₁
end

function pFqmaclaurin(α::NTuple{p, Any}, β::NTuple{q, Any}, z; kwds...) where {p, q}
    T1 = isempty(α) ? Any : mapreduce(typeof, promote_type, α)
    T2 = isempty(β) ? Any : mapreduce(typeof, promote_type, β)
    pFqmaclaurin(T1.(α), T2.(β), z; kwds...)
end

function pFqmaclaurin(a::NTuple{p, S}, b::NTuple{q, U}, z::V; kmax::Int = KMAX) where {p, q, S, U, V}
    T = promote_type(eltype(a), eltype(b), V)
    S₀, S₁, k = one(T), one(T)+prod(a)*z/prod(b), 1
    while k < kmax && errcheck(S₀, S₁, 8eps(real(T)))
        rₖ = z/(k+one(T))
        for i=1:p rₖ *= a[i]+k end
        for i=1:q rₖ /= b[i]+k end
        S₀, S₁ = S₁, S₁+(S₁-S₀)*rₖ
        k += 1
    end
    k < kmax || @warn "Maclaurin approximation to "*pFq2string(Val(p), Val(q))*" reached the maximum degree of "*string(kmax)*"."
    return S₁
end

function continuedfraction(v::V, u::U, rtol::T) where {V<:Function, U<:Function, T<:Number}
    n = 4
    a, b = continuedfraction(v, u, n, 2n)
    n *= 2
    while n <= 2^16
        isapprox(a, b, rtol=rtol) && return b
        n *= 2
        a = deepcopy(b)
        b = continuedfraction(v, u, n)
    end
    isapprox(a, b, rtol=rtol) && return b
    error("Convergence failure for continued fraction for hypergeometric function")
end

function continuedfraction(v::V, u::U, n::Int, m::Int) where {V<:Function, U<:Function}
    @assert m > n
    output_m = u(m) / v(m)
    for i ∈ reverse(n+1:m-1)
        output_m = u(i) / (v(i) + output_m)
    end
    ui, vi = u(n), v(n)
    output_m = ui / (vi + output_m)
    output_n = ui / vi
    for i ∈ reverse(1:n-1)
        ui, vi = u(i), v(i)
        output_m = ui / (vi + output_m)
        output_n = ui / (vi + output_n)
    end
    vi = v(0)
    return vi + output_n, vi + output_m
end

function continuedfraction(v::V, u::U, n::Int) where {V<:Function, U<:Function}
    output = u(n) / v(n)
    for i ∈ reverse(1:n-1)
        output = u(i) / (v(i) + output)
    end
    return v(0) + output
end
