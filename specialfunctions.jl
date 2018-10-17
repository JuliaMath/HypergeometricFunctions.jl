const ρ = 0.72
const ρϵ = 0.71

struct ℕ end

Base.in(n::Integer,::Type{ℕ}) = n > 0
Base.in(n::Real,::Type{ℕ}) = (ν = round(Int,n); n == ν && ν ∈ ℕ)
Base.in(n::Complex,::Type{ℕ}) = imag(n) == 0 && real(n) ∈ ℕ
Base.in(n::Dual,::Type{ℕ}) = dualpart(n) == 0 && realpart(n) ∈ ℕ

struct ℕ₀ end

Base.in(n::Integer,::Type{ℕ₀}) = n ≥ 0
Base.in(n::Real,::Type{ℕ₀}) = (ν = round(Int,n); n == ν && ν ∈ ℕ₀)
Base.in(n::Complex,::Type{ℕ₀}) = imag(n) == 0 && real(n) ∈ ℕ₀
Base.in(n::Dual,::Type{ℕ₀}) = dualpart(n) == 0 && realpart(n) ∈ ℕ₀

struct ℤ end

Base.in(n::Integer,::Type{ℤ}) = true
Base.in(n::Real,::Type{ℤ}) = n == round(Int,n)
Base.in(n::Complex,::Type{ℤ}) = imag(n) == 0 && real(n) ∈ ℤ
Base.in(n::Dual,::Type{ℤ}) = dualpart(n) == 0 && realpart(n) ∈ ℤ

abeqcd(a,b,cd) = isequal(a,b) && isequal(b,cd)
abeqcd(a,b,c,d) = isequal(a,c) && isequal(b,d)

absarg(z) = abs(angle(z))

sqrtatanhsqrt(x) = x == 0 ? one(x) : (s = sqrt(-x); atan(s)/s)
sqrtasinsqrt(x) = x == 0 ? one(x) : (s = sqrt(x); asin(s)/s)
sinnasinsqrt(n,x) = x == 0 ? one(x) : (s = sqrt(x); sin(n*asin(s))/(n*s))
cosnasinsqrt(n,x) = cos(n*asin(sqrt(x)))
expnlog1pcoshatanhsqrt(n,x) = x == 0 ? one(x) : (s = sqrt(x); (exp(n*log1p(s))+exp(n*log1p(-s)))/2)
expnlog1psinhatanhsqrt(n,x) = x == 0 ? one(x) : (s = sqrt(x); (exp(n*log1p(s))-exp(n*log1p(-s)))/(2n*s))

sqrtatanhsqrt(x::Union{T,Dual{T}}) where {T<:Real} = x == 0 ? one(x) : x > 0 ? (s = sqrt(x); atanh(s)/s) : (s = sqrt(-x); atan(s)/s)
sqrtasinsqrt(x::Union{T,Dual{T}}) where {T<:Real} = x == 0 ? one(x) : x > 0 ? (s = sqrt(x); asin(s)/s) : (s = sqrt(-x); asinh(s)/s)
sinnasinsqrt(n,x::Union{T,Dual{T}}) where {T<:Real} = x == 0 ? one(x) : x > 0 ? (s = sqrt(x); sin(n*asin(s))/(n*s)) : (s = sqrt(-x); sinh(n*asinh(s))/(n*s))
cosnasinsqrt(n,x::Union{T,Dual{T}}) where {T<:Real} = x > 0 ? cos(n*asin(sqrt(x))) : cosh(n*asinh(sqrt(-x)))
expnlog1pcoshatanhsqrt(n,x::Union{T,Dual{T}}) where {T<:Real} = x == 0 ? one(x) : x > 0 ? exp(n/2*log1p(-x))*cosh(n*atanh(sqrt(x))) : exp(n/2*log1p(-x))*cos(n*atan(sqrt(-x)))
expnlog1psinhatanhsqrt(n,x::Union{T,Dual{T}}) where {T<:Real} = x == 0 ? one(x) : x > 0 ? (s = sqrt(x); exp(n/2*log1p(-x))*sinh(n*atanh(s))/(n*s)) : (s = sqrt(-x); exp(n/2*log1p(-x))*sin(n*atan(s))/(n*s))

expm1nlog1p(n,x) = x == 0 ? one(x) : expm1(n*log1p(x))/(n*x)

logandpoly(x) = undirected(x) == 0 ? one(x) : 6*(-2undirected(x)+(undirected(x)-2)*log1p(-x))/undirected(x)^3
function logandpoly(x::Union{Float64, ComplexF64})
    if abs(x) > 0.2
        6*(-2x+(x-2)*log1p(-x))/x^3
    else
        logandpolyseries(x)
    end
end

logandpolyseries(x::Union{Float64, Dual128, ComplexF64, DualComplex256}) = @evalpoly(x, 1.0, 1.0, 0.9, 0.8, 0.7142857142857143, 0.6428571428571429, 0.5833333333333334, 0.5333333333333333, 0.4909090909090909, 0.45454545454545453, 0.4230769230769231, 0.3956043956043956, 0.37142857142857144, 0.35, 0.33088235294117646, 0.3137254901960784, 0.2982456140350877, 0.28421052631578947, 0.2714285714285714, 0.2597402597402597)

speciallog(x) = x == 0 ? one(x) : (x > 0 ? (s = sqrt(x); 3(atanh(s)-s)/s^3) : (s = sqrt(-x); 3(s-atan(s))/s^3))
speciallog(x::Directed) = (s = sqrt(-x); 3(s-atan(s))/s^3)
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
speciallogseries(x::Union{Float64,Dual128}) = @clenshaw(5.0x,1.0087391788544393911192,1.220474262857857637288e-01,8.7957928919918696061703e-03,6.9050958578444820505037e-04,5.7037120050065804396306e-05,4.8731405131379353370205e-06,4.2648797509486828820613e-07,3.800372208946157617901e-08,3.434168059359993493634e-09,3.1381484326392473547608e-10,2.8939845618385022798906e-11,2.6892186934806386106143e-12,2.5150879096374730760324e-13,2.3652490233687788117887e-14,2.2349973917002118259929e-15,2.120769988408948118084e-16)
speciallogseries(x::Union{ComplexF64,DualComplex256}) = @evalpoly(x,1.0000000000000000000000,5.9999999999999999999966e-01,4.2857142857142857142869e-01,3.3333333333333333333347e-01,2.7272727272727272727292e-01,2.3076923076923076923072e-01,1.9999999999999999999996e-01,1.7647058823529411764702e-01,1.5789473684210526315786e-01,1.4285714285714285714283e-01,1.3043478260869565217384e-01,1.2000000000000000000000e-01,1.1111111111111111111109e-01,1.0344827586206896551722e-01,9.6774193548387096774217e-02,9.0909090909090909090938e-02,8.5714285714285714285696e-02,8.1081081081081081081064e-02,7.6923076923076923076907e-02,7.3170731707317073170688e-02)

tanpi(z) = sinpi(z)/cospi(z)

const libm = Base.libm_name

unsafe_gamma(x::Float64) = ccall((:tgamma,libm),  Float64, (Float64,), x)
unsafe_gamma(x::Float32) = ccall((:tgammaf,libm),  Float32, (Float32,), x)
unsafe_gamma(x::Real) = unsafe_gamma(float(x))
function unsafe_gamma(x::BigFloat)
    z = BigFloat()
    ccall((:mpfr_gamma, :libmpfr), Int32, (Ref{BigFloat}, Ref{BigFloat}, Int32), z, x, Base.MPFR.ROUNDING_MODE[])
    return z
end
unsafe_gamma(z::Complex) = gamma(z)
unsafe_gamma(z::Dual) = (r = realpart(z);w = unsafe_gamma(r); dual(w, w*digamma(r)*dualpart(z)))

"""
`@lanczosratio(z,ϵ,c₀,c...)`

Compute `∑_{i=1}^N cᵢ/(z-1+i)/(z-1+i+ϵ) / ( c₀ + ∑_{i=1}^N cᵢ/(z-1+i) )`
"""
macro lanczosratio(z, ϵ, c₀, c...)
    ex_num = :(zero(zm1))
    ex_den = esc(c₀)
    for i = 0:length(c)-1
        temp = :(inv(zm1+$i))
        ex_num = :(muladd($(esc(c[i+1])),$temp*inv(zm1pϵ+$i),$ex_num))
        ex_den = :(muladd($(esc(c[i+1])),$temp,$ex_den))
    end
    ex = :($ex_num/$ex_den)
    Expr(:block, :(zm1 = $(esc(z))), :(zm1pϵ = $(esc(z))+$(esc(ϵ))), ex)
end

lanczosratio(z::Union{Float64,ComplexF64,Dual128,DualComplex256},ϵ::Union{Float64,ComplexF64,Dual128,DualComplex256}) = @lanczosratio(z,ϵ,0.99999999999999709182,57.156235665862923517,-59.597960355475491248,14.136097974741747174,-0.49191381609762019978,0.33994649984811888699E-4,0.46523628927048575665E-4,-0.98374475304879564677E-4,0.15808870322491248884E-3,-0.21026444172410488319E-3,0.21743961811521264320E-3,-0.16431810653676389022E-3,0.84418223983852743293E-4,-0.26190838401581408670E-4,0.36899182659531622704E-5)

function lanczosapprox(z::Union{Float64,ComplexF64,Dual128,DualComplex256}, ϵ::Union{Float64,ComplexF64,Dual128,DualComplex256})
    zm0p5 = z-0.5
    zpgm0p5 = zm0p5+4.7421875
    return zm0p5*log1p(ϵ/zpgm0p5) + ϵ*log(zpgm0p5+ϵ) - ϵ + log1p(-ϵ*lanczosratio(z,ϵ))
end

function H(z::Union{Float64,ComplexF64,Dual128,DualComplex256},ϵ::Union{Float64,ComplexF64,Dual128,DualComplex256})
    zm0p5 = z-0.5
    zpgm0p5 = zm0p5+4.7421875
    if real(z) ≥ 1/2
        if z == z+ϵ # ϵ is numerical 0
            zm0p5/zpgm0p5 + log(zpgm0p5) - 1 - lanczosratio(z,ϵ)
        else
            expm1( zm0p5*log1p(ϵ/zpgm0p5) + ϵ*log(zpgm0p5+ϵ) - ϵ + log1p(-ϵ*lanczosratio(z,ϵ)) )/ϵ
        end
    else
        tpz = tanpi(z)
        if z == z+ϵ # ϵ is numerical 0
            H(1-z,ϵ) - π/tpz
        else
            temp = (cospi(ϵ) + sinpi(ϵ)/tpz)*H(1-z,-ϵ) + .5ϵ*(π*sinc(.5ϵ))^2 - π*sinc(ϵ)/tpz
            temp/(1-ϵ*temp)
        end
    end
end

"""
Compute the function (1/Γ(z)-1/Γ(z+ϵ))/ϵ
"""
function G(z::Union{Float64,ComplexF64,Dual128,DualComplex256},ϵ::Union{Float64,ComplexF64,Dual128,DualComplex256})
    n,zpϵ = round(Int,real(z)),z+ϵ
    if abs(ϵ) > 0.1
        (inv(unsafe_gamma(z))-inv(unsafe_gamma(zpϵ)))/ϵ
    elseif z ≠ zpϵ
        m = round(Int,real(zpϵ))
        if z == n && n ≤ 0
            -inv(ϵ*unsafe_gamma(zpϵ))
        elseif zpϵ == m && m ≤ 0
            inv(ϵ*unsafe_gamma(z))
        elseif abs(z+abs(n)) < abs(zpϵ+abs(m))
            H(z,ϵ)/unsafe_gamma(zpϵ)
        else
            H(zpϵ,-ϵ)/unsafe_gamma(z)
        end
    else # ϵ is numerical 0
        if z == n && n ≤ 0
            (-1)^(n+1)*unsafe_gamma(1-n)
        else
            digamma(z)/unsafe_gamma(z)
        end
    end
end

G(z::Number,ϵ::Number) = ϵ == 0 ? digamma(z)/unsafe_gamma(z) : (inv(unsafe_gamma(z))-inv(unsafe_gamma(z+ϵ)))/ϵ


"""
Compute the function ((z+ϵ)ₘ-(z)ₘ)/ϵ
"""
function P(z::Number,ϵ::Number,m::Int)
    n₀ = -round(Int,real(z))
    if ϵ == 0
        if 0 ≤ n₀ < m
            ret1,ret2,n = one(z),zero(z),0
            while n < m
                n == n₀ && (n+=1; continue)
                ret1 *= z+n
                ret2 += inv(z+n)
                n+=1
            end
            ret1 + pochhammer(z,m)*ret2
        else
            ret = zero(z)
            for n=0:m-1
                ret += inv(z+n)
            end
            pochhammer(z,m)*ret
        end
    else
        if 0 ≤ n₀ < m
            zpϵ,ret1,ret2,n = z+ϵ,one(z),zero(z),0
            while n < m
                n == n₀ && (n+=1; continue)
                ret1 *= zpϵ+n
                ret2 += log1p(ϵ/(z+n))
                n+=1
            end
            ret1 + pochhammer(z,m)*expm1(ret2)/ϵ
        else
            ret = zero(z)
            for n=0:m-1
                ret += log1p(ϵ/(z+n))
            end
            pochhammer(z,m)*expm1(ret)/ϵ
        end
    end
end

E(z::Number,ϵ::Number) = ϵ == 0 ? z : expm1(ϵ*z)/ϵ

G(z::AbstractVector{BigFloat},ϵ::BigFloat) = BigFloat[G(zi,ϵ) for zi in z]

# Transformation formula w = 1-z

reconeα₀(a,b,c,m::Int,ϵ) = ϵ == 0 ? (-1)^m*gamma(m)*gamma(c)/(gamma(a+m)*gamma(b+m)) : gamma(c)/(ϵ*gamma(1-m-ϵ)*gamma(a+m+ϵ)*gamma(b+m+ϵ))
reconeβ₀(a,b,c,w,m::Int,ϵ) = abs(ϵ) > 0.1 ? ( pochhammer(float(a),m)*pochhammer(b,m)/(gamma(1-ϵ)*gamma(a+m+ϵ)*gamma(b+m+ϵ)*gamma(m+1)) - w^ϵ/(gamma(a)*gamma(b)*gamma(m+1+ϵ)) )*gamma(c)*w^m/ϵ : ( (G(1.0,-ϵ)/gamma(m+1)+G(m+1.0,ϵ))/(gamma(a+m+ϵ)*gamma(b+m+ϵ)) - (G(float(a)+m,ϵ)/gamma(b+m+ϵ)+G(float(b)+m,ϵ)/gamma(a+m))/gamma(m+1+ϵ) - E(log(w),ϵ)/(gamma(a+m)*gamma(b+m)*gamma(m+1+ϵ)) )*gamma(c)*pochhammer(float(a),m)*pochhammer(b,m)*w^m
reconeγ₀(a,b,c,w,m::Int,ϵ) = gamma(c)*pochhammer(float(a),m)*pochhammer(b,m)*w^m/(gamma(a+m+ϵ)*gamma(b+m+ϵ)*gamma(m+1)*gamma(1-ϵ))

function Aone(a,b,c,w,m::Int,ϵ)
    αₙ = reconeα₀(a,b,c,m,ϵ)*one(w)
    ret = m ≤ 0 ? zero(w) : αₙ
    for n = 0:m-2
        αₙ *= (a+n)*(b+n)/((n+1)*(1-m-ϵ+n))*w
        ret += αₙ
    end
    ret
end

function Bone(a,b,c,w,m::Int,ϵ)
    βₙ,γₙ = reconeβ₀(a,b,c,w,m,ϵ)*one(w),reconeγ₀(a,b,c,w,m,ϵ)*w
    ret,err,n = βₙ,1.0,0
    while err > 10eps()
        βₙ = (a+m+n+ϵ)*(b+m+n+ϵ)/((m+n+1+ϵ)*(n+1))*w*βₙ + ( (a+m+n)*(b+m+n)/(m+n+1) - (a+m+n) - (b+m+n) - ϵ + (a+m+n+ϵ)*(b+m+n+ϵ)/(n+1) )*γₙ/((m+n+1+ϵ)*(n+1-ϵ))
        ret += βₙ
        err = errcheck(βₙ,ret)
        γₙ *= (a+m+n)*(b+m+n)/((m+n+1)*(n+1-ϵ))*w
        n+=1
    end
    ret
end

function _₂F₁one(a,b,c,z)
    m = round(Int,real(c-(a+b)))
    ϵ = c-(a+b)-m
    w = 1-z
    (-1)^m/sinc(ϵ)*(Aone(a,b,c,undirected(w),m,ϵ) + Bone(a,b,c,w,m,ϵ))
end

# Transformation formula w = 1/z

recInfα₀(a,b,c,m::Int,ϵ) = ϵ == 0 ? (-1)^m*gamma(m)*gamma(c)/(gamma(a+m)*gamma(c-a)) : gamma(c)/(ϵ*gamma(1-m-ϵ)*gamma(a+m+ϵ)*gamma(c-a))
recInfβ₀(a,b,c,w,m::Int,ϵ) = abs(ϵ) > 0.1 ?
                        ( pochhammer(float(a),m)*pochhammer(float(1-c+a),m)/(gamma(1-ϵ)*gamma(a+m+ϵ)*gamma(c-a)*gamma(m+1)) -
                                (-w)^ϵ*pochhammer(float(1-c+a)+ϵ,m)/(gamma(a)*gamma(c-a-ϵ)*gamma(m+1+ϵ)) )*gamma(c)*w^m/ϵ :
                        ( (pochhammer(float(1-c+a)+ϵ,m)*G(1.0,-ϵ)-P(1-c+a,ϵ,m)/gamma(1-ϵ))/(gamma(c-a)*gamma(a+m+ϵ)*gamma(m+1)) +
                                pochhammer(float(1-c+a)+ϵ,m)*( (G(m+1.0,ϵ)/gamma(a+m+ϵ) - G(float(a)+m,ϵ)/gamma(m+1+ϵ))/gamma(c-a) -
                                    (G(float(c-a),-ϵ) - E(-log(-w),-ϵ)/gamma(c-a-ϵ))/(gamma(m+1+ϵ)*gamma(a+m)) ) )*gamma(c)*pochhammer(float(a),m)*w^m
recInfγ₀(a,b,c,w,m::Int,ϵ) = gamma(c)*pochhammer(float(a),m)*pochhammer(float(1-c+a),m)*w^m/(gamma(a+m+ϵ)*gamma(c-a)*gamma(m+1)*gamma(1-ϵ))

function AInf(a,b,c,w,m::Int,ϵ)
    αₙ = recInfα₀(a,b,c,m,ϵ)*one(w)
    ret = m ≤ 0 ? zero(w) : αₙ
    for n = 0:m-2
        αₙ *= (a+n)*(1-c+a+n)/((n+1)*(1-m-ϵ+n))*w
        ret += αₙ
    end
    ret
end

function BInf(a,b,c,win,m::Int,ϵ)
    w=undirected(win)
    βₙ,γₙ = recInfβ₀(a,b,c,win,m,ϵ)*one(w),recInfγ₀(a,b,c,win,m,ϵ)*w
    ret,err,n = βₙ,1.0,0
    while err > 10eps()
        βₙ = (a+m+n+ϵ)*(1-c+a+m+n+ϵ)/((m+n+1+ϵ)*(n+1))*w*βₙ + ( (a+m+n)*(1-c+a+m+n)/(m+n+1) - (a+m+n) - (1-c+a+m+n) - ϵ + (a+m+n+ϵ)*(1-c+a+m+n+ϵ)/(n+1) )*γₙ/((m+n+1+ϵ)*(n+1-ϵ))
        ret += βₙ
        err = errcheck(βₙ,ret)
        γₙ *= (a+m+n)*(1-c+a+m+n)/((m+n+1)*(n+1-ϵ))*w
        n+=1
    end
    ret
end

function _₂F₁Inf(a,b,c,z)
    m = round(Int,real(b-a))
    ϵ = b-a-m
    w = reverseorientation(inv(z))  # we've swapped the branch cut
    (-1)^m*(-w)^a/sinc(ϵ)*(AInf(a,b,c,undirected(w),m,ϵ) + BInf(a,b,c,w,m,ϵ))
end


function _₂F₁maclaurin(a::Number,b::Number,c::Number,z::Number)
    T = promote_type(typeof(a),typeof(b),typeof(c),typeof(z))
    S₀,S₁,err,j = one(T),one(T)+a*b*z/c,one(real(T)),1
    while err > 10eps(T)
        rⱼ = (a+j)/(j+1)*(b+j)/(c+j)
        S₀,S₁ = S₁,S₁+(S₁-S₀)*rⱼ*z
        err = errcheck(S₁-S₀,S₀)
        j+=1
    end
    return S₁
end

function _₂F₁maclaurinalt(a::Number,b::Number,c::Number,z::Number)
    T = promote_type(typeof(a),typeof(b),typeof(c),typeof(z))
    C,S,err,j = one(T),one(T),one(real(T)),0
    while err > 10eps(T)
        C *= (a+j)/(j+1)*(b+j)/(c+j)*z
        S += C
        err = errcheck(C,S)
        j+=1
    end
    return S
end

function _₂F₁continuation(s::Number,t::Number,c::Number,z₀::Number,z::Number)
    T = promote_type(typeof(s),typeof(t),typeof(c),typeof(z₀),typeof(z))
    izz₀,d0,d1 = inv(z-z₀),one(T),s/(2s-t+one(T))*((s+1)*(1-2z₀)+(t+1)*z₀-c)
    S₀,S₁,izz₀j,err,j = one(T),one(T)+d1*izz₀,izz₀,one(real(T)),2
    while err > 10eps(T)
        d0,d1,izz₀j = d1,(j+s-one(T))/j/(j+2s-t)*(((j+s)*(1-2z₀)+(t+1)*z₀-c)*d1 + z₀*(1-z₀)*(j+s-2)*d0),izz₀j*izz₀
        S₀,S₁ = S₁,S₁+d1*izz₀j
        err = errcheck(S₁-S₀,S₀)
        j+=1
    end
    return S₁
end

function _₂F₁continuationalt(a::Number,c::Number,z₀::Number,z::Number)
    T = promote_type(typeof(a),typeof(c),typeof(z₀),typeof(z))
    izz₀ = inv(z-z₀)
    e0,e1 = one(T),(a+one(T))*(one(T)-2z₀)+(2a+one(T))*z₀-c
    f0,f1 = zero(T),one(T)-2z₀
    cⱼ = log(z₀-z)+2digamma(one(T))-digamma(a)-digamma(c-a)
    S₀ = cⱼ
    cⱼ += 2/one(T)-one(T)/a
    C = a*izz₀
    S₁,err,j = S₀+(e1*cⱼ-f1)*C,one(real(T)),2
    while err > 10eps(T)
        f0,f1 = f1,(((j+a)*(1-2z₀)+(2a+1)*z₀-c)*f1+z₀*(1-z₀)*(j-1)*f0+(1-2z₀)*e1+2z₀*(1-z₀)*e0)/j
        e0,e1 = e1,(((j+a)*(1-2z₀)+(2a+1)*z₀-c)*e1+z₀*(1-z₀)*(j-1)*e0)/j
        C *= (a+j-1)*izz₀/j
        cⱼ += 2/T(j)-one(T)/(a+j-one(T))
        S₀,S₁ = S₁,S₁+(e1*cⱼ-f1)*C
        err = errcheck(S₁-S₀,S₀)
        j+=1
    end
    return S₁
end

function _₂F₁logsum(a::Number,b::Number,z::Number,w::Number,s::Int)
    T = promote_type(typeof(a),typeof(b),typeof(z),typeof(w))
    cⱼ = 2digamma(one(T))-digamma(a)-digamma(b)+s*log1p(-z)
    C,S,err,j = one(T),cⱼ,one(real(T)),0
    while err > 10eps(T)
        C *= (a+j)/(j+1)^2*(b+j)*w
        cⱼ += 2/(j+one(T))-one(T)/(a+j)-one(T)/(b+j)
        S += C*cⱼ
        err = errcheck(C,S)
        j+=1
    end
    return S
end

function _₂F₁logsumalt(a::Number,b::Number,z::Number,w::Number)
    T = promote_type(typeof(a),typeof(b),typeof(z),typeof(w))
    d,cⱼ = one(T)-b,2digamma(one(T))-digamma(a)-digamma(b)-log(-w)
    C,S,err,j = one(T),cⱼ,one(real(T)),0
    while err > 10eps(T)
        C *= (a+j)/(j+1)^2*(d+j)*w
        cⱼ += 2/(j+one(T))-one(T)/(a+j)+one(T)/(b-(j+one(T)))
        S += C*cⱼ
        err = errcheck(C,S)
        j+=1
    end
    return S
end

function _₂F₁taylor(a::Number,b::Number,c::Number,z::Number)
    T = promote_type(typeof(a),typeof(b),typeof(c),typeof(z))
    z₀ = abs(z) < 1 ? ρϵ*sign(z) : sign(z)/ρϵ
    q₀,q₁ = _₂F₁(a,b,c,z₀),a*b/c*_₂F₁(a+1,b+1,c+1,z₀)
    S₀,zz₀ = q₀,z-z₀
    S₁,err,zz₀j,j = S₀+q₁*zz₀,one(real(T)),zz₀,0
    while err > 10eps(T)
        q₀,q₁ = q₁,((j*(2z₀-one(T))-c+(a+b+one(T))*z₀)*q₁ + (a+j)*(b+j)/(j+one(T))*q₀)/(z₀*(one(T)-z₀)*(j+2))
        zz₀j *= zz₀
        S₀,S₁ = S₁,S₁+q₁*zz₀j
        err = errcheck(S₀-S₁,S₀)
        j+=1
    end
    return S₁
end

function _₃F₂maclaurin(a₁::Number,a₂::Number,a₃::Number,b₁::Number,b₂::Number,z::Number)
    T = promote_type(typeof(a₁),typeof(a₂),typeof(a₃),typeof(b₁),typeof(b₂),typeof(z))
    S₀,S₁,err,j = one(T),one(T)+(a₁*a₂*a₃*z)/(b₁*b₂),one(real(T)),1
    while err > 100eps(T)
        rⱼ = ((a₁+j)*(a₂+j)*(a₃+j))/((b₁+j)*(b₂+j)*(j+1))
        S₀,S₁ = S₁,S₁+(S₁-S₀)*rⱼ*z
        err = errcheck(S₁-S₀,S₀)
        j+=1
    end
    return S₁
end

function mFnmaclaurin(a::AbstractVector{S},b::AbstractVector{V},z::Number) where {S<:Number,V<:Number}
    T = promote_type(S,V,typeof(z))
    S₀,S₁,err,j = one(T),one(T)+prod(a)*z/prod(b),one(real(T)),1
    while err > 100eps(T)
        rⱼ = inv(j+one(T))
        for i=1:length(a) rⱼ *= a[i]+j end
        for i=1:length(b) rⱼ /= b[i]+j end
        S₀,S₁ = S₁,S₁+(S₁-S₀)*rⱼ*z
        err = errcheck(S₁-S₀,S₀)
        j+=1
    end
    return S₁
end

hypot2(a,b) = hypot(a,b)
hypot2(a::Complex{T},b::Complex{T}) where {T<:AbstractFloat} = hypot(hypot(a.re,b.re),hypot(a.im,b.im))
hypot2(a::T,b::Complex{T}) where {T<:AbstractFloat} = hypot(hypot(a,b.re),b.im)
hypot2(a::Complex{T},b::T) where {T<:AbstractFloat} = hypot(hypot(a.re,b),a.im)

errcheck(x,y) = abs(x/y)
errcheck(x::Dual,y::Dual) = hypot2(realpart(x),dualpart(x))/hypot2(realpart(y),dualpart(y))

# This is from Base because real(::Type{Dual{BigFloat}}) doesn't exist in the scop of rtoldefault, required by isapprox.

# isapprox: approximate equality of numbers
function isapprox(x::Number, y::Number; rtol::Real=rtoldefault(x,y), atol::Real=0)
    x == y || (isfinite(x) && isfinite(y) && abs(x-y) <= atol + rtol*max(abs(x), abs(y)))
end

const ≈ = isapprox
≉(x,y) = !(x ≈ y)

# default tolerance arguments
rtoldefault(::Type{T}) where {T<:AbstractFloat} = sqrt(eps(T))
rtoldefault(::Type{T}) where {T<:Real} = 0
rtoldefault(::Type{Dual{T}}) where {T<:Real} = rtoldefault(T)
rtoldefault(x::Union{T,Type{T}}, y::Union{S,Type{S}}) where {T<:Number,S<:Number} = rtoldefault(promote_type(real(T),real(S)))
