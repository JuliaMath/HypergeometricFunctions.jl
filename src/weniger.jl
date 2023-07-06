# Rational approximations to generalized hypergeometric functions
# using Weniger's sequence transformation

# ₀F₀(;z), γ = 2.
function pFqweniger(::Tuple{}, ::Tuple{}, z::T; kmax::Int = KMAX) where T
    if norm(z) < eps(real(T))
        return one(T)
    end
    μlo = one(T)
    μhi = inv(2-z)
    Rlo = one(T)
    Rhi = Rlo + 2z*μhi
    k = 1
    z2 = z*z
    while k < kmax && errcheck(Rlo, Rhi, 8eps(real(T)))
        μhi, μlo = inv(4k+2 + z2*μhi), μhi
        Rhi, Rlo = ((4k+2)*Rhi + z2*Rlo*μlo)*μhi, Rhi
        k += 1
    end
    k < kmax || @warn "Rational approximation to "*pFq2string(Val{0}(), Val{0}())*" reached the maximum type of ("*string(kmax, ", ", kmax)*")."
    return isfinite(Rhi) ? Rhi : Rlo
end

# ₁F₀(α;z), γ = 2.
function pFqweniger(α::Tuple{T1}, ::Tuple{}, z::T2; kmax::Int = KMAX) where {T1, T2}
    α = α[1]
    T = promote_type(T1, T2)
    absα = abs(T(α))
    if norm(z) < eps(real(T)) || norm(α) < eps(absα)
        return one(T)
    end
    μlo = α
    μhi = inv(2-(α+1)*z)
    Rlo = one(T)
    Rhi = Rlo + 2z*μhi*μlo
    if norm(α+1) < eps(absα+1)
        return Rhi
    end
    μhi *= α+1
    k = 1
    z2 = z*z
    while k < kmax && errcheck(Rlo, Rhi, 8eps(real(T)))
        μhi, μlo = inv((2k+1)*(2-z) - (k-α)*z2*μhi), μhi
        Rhi, Rlo = ((2k+1)*(2-z)*Rhi - (k-α)*z2*Rlo*μlo)*μhi, Rhi
        if norm(α+k+1) < eps(absα+k+1)
            return Rhi
        end
        μhi *= α+k+1
        k += 1
    end
    k < kmax || @warn "Rational approximation to "*pFq2string(Val{1}(), Val{0}())*" reached the maximum type of ("*string(kmax, ", ", kmax)*")."
    return isfinite(Rhi) ? Rhi : Rlo
end

# ₂F₀(α,β;z), algorithm γ = 2.
function pFqweniger(α::Tuple{T1, T1}, ::Tuple{}, z::T2; kmax::Int = KMAX) where {T1, T2}
    (α, β) = α
    T = promote_type(T1, T2)
    absα = abs(T(α))
    absβ = abs(T(β))
    if norm(z) < eps(real(T)) || norm(α*β) < eps(absα*absβ)
        return one(T)
    end
    μlo = T(α*β)
    Rlo = one(T)
    a0 = T((α+1)*(β+1))
    μmid = inv(2-a0*z)
    Rmid = Rlo + 2z*μmid*μlo
    if norm(a0) < eps((absα+1)*(absβ+1))
        return Rmid
    end
    μmid *= a0
    k = 1
    a0 = T((α+2)*(β+2))
    t0 = 6-(6-3*α*β)*z
    t1 = 6-2*T(2*α*β+α+β-1)*z
    μhi = inv(t0 - t1*z*μmid)
    Rhi = (t0*Rmid - (t1*Rlo + 6*z*μlo)*z*μmid)*μhi
    if norm(a0) < eps((absα+2)*(absβ+2))
        return Rhi
    end
    μhi *= a0
    k = 2
    z2 = z*z
    while k < 3 || (k < kmax && errcheck(Rmid, Rhi, 8eps(real(T))))
        a0 = T((α+k+1)*(β+k+1))
        t0 = (4k+2)-T((k*(α+β+3k)-(α+1)*(β+1)))*T(2k+1)/T(2k-1)*z
        t1 = (4k+2)+T(k*(3k-α-β)-(α+1)*(β+1))*z
        a2 = T((α+1-k)*(β+1-k))*T(2k+1)/T(2k-1)*z2
        μhi, μmid, μlo = inv(t0 - (t1 + a2*μmid)*z*μhi), μhi, μmid
        Rhi, Rmid, Rlo = (t0*Rhi - (t1*Rmid + a2*Rlo*μlo)*z*μmid)*μhi, Rhi, Rmid
        if norm(a0) < eps((absα+k+1)*(absβ+k+1))
            return Rhi
        end
        μhi *= a0
        k += 1
    end
    k < kmax || @warn "Rational approximation to "*pFq2string(Val{2}(), Val{0}())*" reached the maximum type of ("*string(kmax, ", ", kmax)*")."
    return isfinite(Rhi) ? Rhi : isfinite(Rmid) ? Rmid : Rlo
end

# ₁F₁(α,β;z), algorithm γ = 2.
function pFqweniger(α::Tuple{T1}, β::Tuple{T2}, z::T3; kmax::Int = KMAX) where {T1, T2, T3}
    α = α[1]
    β = β[1]
    T = promote_type(T1, T2, T3)
    absα = abs(T(α))
    if norm(z) < eps(real(T)) || norm(α) < eps(absα)
        return one(T)
    end
    ζ = inv(z)
    Nlo = β*ζ/α
    Dlo = β*ζ/α
    Tlo = Nlo/Dlo
    a0 = T(α+1)
    b0 = T(2*(β+1))
    Nmid = (b0*ζ-a0)*Nlo + b0*ζ
    Dmid = (b0*ζ-a0)*Dlo
    Tmid = Nmid/Dmid
    if norm(a0) < eps(absα+1)
        return Tmid
    end
    Nmid /= a0
    Dmid /= a0
    k = 1
    a0 = T(α+2)
    #a1 = -T(4α+2)
    b0 = T(6*(β+2))
    b1 = T(-6*β)
    t0 = b0*ζ+3α
    t1 = b1*ζ+4α+2
    Nhi = t0*Nmid + t1*Nlo + b1*ζ
    Dhi = t0*Dmid + t1*Dlo
    #Nhi = -a0*Nmid - a1*(Nmid+Nlo) + ζ*(b0*Nmid + b1*Nlo) + b1*ζ
    #Dhi = -a0*Dmid - a1*(Dmid+Dlo) + ζ*(b0*Dmid + b1*Dlo)
    Thi = Nhi/Dhi
    if norm(a0) < eps(absα+2)
        return Thi
    end
    Nhi /= a0
    Dhi /= a0
    k = 2
    while k < 5 || (k < kmax && errcheck(Tmid, Thi, 8eps(real(T))))
        a0 = T(α+k+1)
        #a1 = -T(2α+1)*T(2k)/T(2k-1)
        a2 = T(α+1-k)*T(2k+1)/T(2k-1)
        b0 = T(4k+2)*T(β+k+1)
        b1 = T(4k+2)*T(k-β-1)
        t0 = b0*ζ+a2
        t1 = b1*ζ+a0
        Nhi, Nmid, Nlo = t0*Nhi + t1*Nmid - a2*Nlo, Nhi, Nmid
        Dhi, Dmid, Dlo = t0*Dhi + t1*Dmid - a2*Dlo, Dhi, Dmid
        #Nhi, Nmid, Nlo = -a0*Nhi - a1*(Nhi+Nmid) - a2*(Nmid+Nlo) + ζ*(b0*Nhi + b1*Nmid), Nhi, Nmid
        #Dhi, Dmid, Dlo = -a0*Dhi - a1*(Dhi+Dmid) - a2*(Dmid+Dlo) + ζ*(b0*Dhi + b1*Dmid), Dhi, Dmid
        Thi, Tmid, Tlo = Nhi/Dhi, Thi, Tmid
        if norm(a0) < eps(absα+k+1)
            return Thi
        end
        Nhi /= a0
        Dhi /= a0
        k += 1
    end
    k < kmax || @warn "Rational approximation to "*pFq2string(Val{1}(), Val{1}())*" reached the maximum type of ("*string(kmax, ", ", kmax)*")."
    return isfinite(Thi) ? Thi : isfinite(Tmid) ? Tmid : Tlo
end

# ₂F₁(α,β,γ;z), algorithm γ = 2.
function pFqweniger(α::Tuple{T1, T1}, β::Tuple{T2}, z::T3; kmax::Int = KMAX) where {T1, T2, T3}
    γ = β[1]
    (α, β) = α
    T = promote_type(T1, T2, T3)
    absα = abs(T(α))
    absβ = abs(T(β))
    if norm(z) < eps(real(T)) || norm(α*β) < eps(absα*absβ)
        return one(T)
    end
    μlo = T(α*β)/T(γ)
    Rlo = one(T)
    a0 = T((α+1)*(β+1))
    b0 = T(2*(γ+1))
    μmid = inv(b0-a0*z)
    Rmid = Rlo + b0*z*μmid*μlo
    if norm(a0) < eps((absα+1)*(absβ+1))
        return Rmid
    end
    μmid *= a0
    k = 1
    a0 = T((α+2)*(β+2))
    b0 = T(6*(γ+2))
    b1 = T(-6*γ)
    t0 = b0-(6-3*α*β)*z
    t1 = b1+2*T(2*α*β+α+β-1)*z
    μhi = inv(t0 + t1*z*μmid)
    Rhi = (t0*Rmid + (t1*Rlo + b1*z*μlo)*z*μmid)*μhi
    if norm(a0) < eps((absα+2)*(absβ+2))
        return Rhi
    end
    μhi *= a0
    k = 2
    z2 = z*z
    while k < 5 || (k < kmax && errcheck(Rmid, Rhi, 8eps(real(T))))
        a0 = T((α+k+1)*(β+k+1))
        a2 = T((α+1-k)*(β+1-k))*T(2k+1)/T(2k-1)*z2
        b0 = T(4k+2)*T(γ+k+1)
        b1 = T(4k+2)*T(k-γ-1)
        t0 = b0-T((k*(α+β+3k)-(α+1)*(β+1)))*T(2k+1)/T(2k-1)*z
        t1 = b1-T(k*(3k-α-β)-(α+1)*(β+1))*z
        μhi, μmid, μlo = inv(t0 + (t1 - a2*μmid)*z*μhi), μhi, μmid
        Rhi, Rmid, Rlo = (t0*Rhi + (t1*Rmid - a2*Rlo*μlo)*z*μmid)*μhi, Rhi, Rmid
        if norm(a0) < eps((absα+k+1)*(absβ+k+1))
            return Rhi
        end
        μhi *= a0
        k += 1
    end
    k < kmax || @warn "Rational approximation to "*pFq2string(Val{2}(), Val{1}())*" reached the maximum type of ("*string(kmax, ", ", kmax)*")."
    return isfinite(Rhi) ? Rhi : isfinite(Rmid) ? Rmid : Rlo
end

# ₘFₙ(α;β;z)
# γ ∉ ℕ
function pFqweniger(α::AbstractVector{T1}, β::AbstractVector{T2}, z::T3, args...; kwds...) where {T1, T2, T3}
    pFqweniger(Tuple(α), Tuple(β), z, args...; kwds...)
end
function pFqweniger(α::NTuple{p, Any}, β::NTuple{q, Any}, z, args...; kwds...) where {p, q}
    T1 = isempty(α) ? Any : mapreduce(typeof, promote_type, α)
    T2 = isempty(β) ? Any : mapreduce(typeof, promote_type, β)
    pFqweniger(T1.(α), T2.(β), z, args...; kwds...)
end
function pFqweniger(α::NTuple{p, T1}, β::NTuple{q, T2}, z::T3; kmax::Int = KMAX, γ::T4 = 2) where {p, q, T1, T2, T3, T4}
    T = promote_type(eltype(α), eltype(β), T3, T4)
    absα = abs.(T.(α))
    if norm(z) < eps(real(T)) || norm(prod(α)) < eps(real(T)(prod(absα)))
        return one(T)
    end
    γ = T(γ)
    ζ = inv(z)
    r = max(p+1, q+1)
    N = zeros(T, r+3)
    D = zeros(T, r+3)
    R = zeros(T, r+3)
    N[r+3] = prod(β)*ζ/prod(α)/(γ-1)
    D[r+3] = prod(β)*ζ/prod(α)/(γ-1)
    R[r+3] = N[r+3]/D[r+3]
    err = real(T)(1)
    @inbounds for j in 1:p
        err *= absα[j]+1
    end
    P̂ = zeros(T, r+2)
    t = T(1)
    @inbounds for j in 1:p
        t *= α[j]+1
    end
    P̂[1] = t
    Q = zeros(T, r+1)
    t = T(2)
    @inbounds for j in 1:q
        t *= β[j]+1
    end
    Q[1] = t
    P̂R = γ+2
    QR = T(1)
    k = 0
    @inbounds while k ≤ r+2 || (k < kmax && errcheck(R[r+2], R[r+3], 8eps(real(T))))
        for j in 1:r+2
            N[j] = N[j+1]
            D[j] = D[j+1]
            R[j] = R[j+1]
        end
        t1 = zero(T)
        for j in 0:r
            t1 += Q[j+1]*(γ+2k-2j-1)*N[r-j+2]
        end
        if k ≤ r
            for j in 0:k
                t2 = one(T)
                for i in 1:q
                    t2 *= β[i]+j+1
                end
                t2 *= j+2
                t1 += binomial(k, j)*(-one(T))^(k-j)*t2
            end
        end
        t2 = zero(T)
        t2 += P̂[1]*(γ+k-1)*N[r+2]
        for j in 1:r+1
            t2 += P̂[j+1]*(N[r-j+3]+(γ+k-j-1)*N[r-j+2])
        end
        N[r+3] = ζ*t1-t2
        t1 = zero(T)
        for j in 0:r
            t1 += Q[j+1]*(γ+2k-2j-1)*D[r-j+2]
        end
        t2 = zero(T)
        t2 += P̂[1]*(γ+k-1)*D[r+2]
        for j in 1:r+1
            t2 += P̂[j+1]*(D[r-j+3]+(γ+k-j-1)*D[r-j+2])
        end
        D[r+3] = ζ*t1-t2
        R[r+3] = N[r+3]/D[r+3]
        if norm(P̂[1]) < eps(err)
            return R[r+3]
        end
        N[r+3] /= P̂[1]
        D[r+3] /= P̂[1]
        k += 1
        err = real(T)(1)
        for j in 1:p
            err *= absα[j]+k+1
        end
        if k ≤ r+1
            t2 = T(1)
            for j in 1:p
                t2 *= α[j]+k+1
            end
            t2 /= P̂R
            P̂R *= ((γ+2k+1)*(γ+2k+2))/(γ+k+1)
            t1 = k*((γ+2k)*t2 - P̂[1])
            for j in 2:k
                s = ((k-j+1)*((γ+2k-j+1)*t1+k*P̂[j-1]) - k*P̂[j])/j
                P̂[j-1] = t2
                t2 = t1
                t1 = s
            end
            P̂[k] = t2
            P̂[k+1] = t1
            t2 = T(k+2)
            for j in 1:q
                t2 *= β[j]+k+1
            end
            QR *= ((γ+2k-2)*(γ+2k-1))/(γ+k-1)
            t2 /= QR
            t1 = k*((γ+2k-1)*t2 - Q[1])
            for j in 2:min(k, r)
                s = ((k-j+1)*((γ+2k-j)*t1+k*Q[j-1]) - k*Q[j])/j
                Q[j-1] = t2
                t2 = t1
                t1 = s
            end
            Q[min(k, r)] = t2
            Q[min(k, r)+1] = t1
        else
            t2 = γ+k
            for j in 1:p
                t2 *= α[j]+k+1
            end
            t2 /= P̂R
            P̂R *= ((γ+2k+1)*(γ+2k+2))/((γ+2k-r-1)*(γ+2k-r))
            t1 = k*((γ+2k)*t2 - (γ+2k-1-r-2)*P̂[1])
            for j in 2:r+1
                s = ((k-j+1)*((γ+2k-j+1)*t1-(r-j+3)*k*P̂[j-1]) - (γ+2k-j-r-2)*k*P̂[j])/j
                P̂[j-1] = t2
                t2 = t1
                t1 = s
            end
            P̂[r+1] = t2
            P̂[r+2] = t1
            t2 = T(k+2)
            for j in 1:q
                t2 *= β[j]+k+1
            end
            QR *= ((γ+2k-2)*(γ+2k-1))/((γ+2k-r-3)*(γ+2k-r-2))
            t2 /= QR
            t1 = k*((γ+2k-1)*t2 - (γ+2k-1-r-2)*Q[1])
            for j in 2:r
                s = ((k-j+1)*((γ+2k-j)*t1-(r-j+2)*k*Q[j-1]) - (γ+2k-j-r-2)*k*Q[j])/j
                Q[j-1] = t2
                t2 = t1
                t1 = s
            end
            Q[r] = t2
            Q[r+1] = t1
        end
    end
    k < kmax || @warn "Rational approximation to "*pFq2string(Val{p}(), Val{q}())*" reached the maximum type of ("*string(kmax, ", ", kmax)*")."
    return isfinite(R[r+3]) ? R[r+3] : R[r+2]
end
