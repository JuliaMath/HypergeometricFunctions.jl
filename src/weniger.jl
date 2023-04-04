# Rational approximations to generalized hypergeometric functions
# using Weniger's sequence transformation

# ₀F₀(;z), γ = 2.
function pFqweniger(::Tuple{}, ::Tuple{}, z::T; kmax::Int = 10_000) where T
    if norm(z) < eps(real(T))
        return one(T)
    end
    ζ = inv(z)
    Nlo = ζ
    Dlo = ζ
    Tlo = Nlo/Dlo
    Nhi = (2ζ - 1)*Nlo + 2ζ
    Dhi = (2ζ - 1)*Dlo
    Thi = Nhi/Dhi
    k = 1
    while k < kmax && errcheck(Tlo, Thi, 10eps(real(T)))
        Nhi, Nlo = (4k+2)*ζ*Nhi + Nlo, Nhi
        Dhi, Dlo = (4k+2)*ζ*Dhi + Dlo, Dhi
        Thi, Tlo = Nhi/Dhi, Thi
        k += 1
    end
    return isfinite(Thi) ? Thi : Tlo
end

# ₁F₀(α;z), γ = 2.
function pFqweniger(α::Tuple{T1}, ::Tuple{}, z::T2; kmax::Int = 10_000) where {T1, T2}
    α = α[1]
    T = promote_type(T1, T2)
    absα = abs(T(α))
    if norm(z) < eps(real(T)) || norm(α) < eps(absα)
        return one(T)
    end
    ζ = inv(z)
    Nlo = ζ/α
    Dlo = ζ/α
    Tlo = Nlo/Dlo
    Nhi = (2ζ - (α+1))*Nlo + 2ζ
    Dhi = (2ζ - (α+1))*Dlo
    Thi = Nhi/Dhi
    if norm(α+1) < eps(absα+1)
        return Thi
    end
    Nhi /= α+1
    Dhi /= α+1
    k = 1
    while k < kmax && errcheck(Tlo, Thi, 10eps(real(T)))
        Nhi, Nlo = (2k+1)*(2ζ-1)*Nhi - (k-α)*Nlo, Nhi
        Dhi, Dlo = (2k+1)*(2ζ-1)*Dhi - (k-α)*Dlo, Dhi
        Thi, Tlo = Nhi/Dhi, Thi
        if norm(α+k+1) < eps(absα+k+1)
            return Thi
        end
        Nhi /= α+k+1
        Dhi /= α+k+1
        k += 1
    end
    return isfinite(Thi) ? Thi : Tlo
end

# ₂F₀(α,β;z), algorithm γ = 2.
function pFqweniger(α::Tuple{T1, T1}, ::Tuple{}, z::T2; kmax::Int = 10_000) where {T1, T2}
    (α, β) = α
    T = promote_type(T1, T2)
    absα = abs(T(α))
    absβ = abs(T(β))
    if norm(z) < eps(real(T)) || norm(α*β) < eps(absα*absβ)
        return one(T)
    end
    ζ = inv(z)
    Nlo = ζ/(α*β)
    Dlo = ζ/(α*β)
    Tlo = Nlo/Dlo
    a0 = (α+1)*(β+1)
    Nmid = (2ζ-a0)*Nlo + 2ζ
    Dmid = (2ζ-a0)*Dlo
    Tmid = Nmid/Dmid
    if norm(a0) < eps((absα+1)*(absβ+1))
        return Tmid
    end
    Nmid /= a0
    Dmid /= a0
    k = 1
    a0 = (α+2)*(β+2)
    a1 = 2*(2-(2*α*β+α+β+1))
    Nhi = -a0*Nmid - a1*(Nmid+Nlo) + 6ζ*(Nmid - Nlo) - 6ζ
    Dhi = -a0*Dmid - a1*(Dmid+Dlo) + 6ζ*(Dmid - Dlo)
    Thi = Nhi/Dhi
    if norm(a0) < eps((absα+2)*(absβ+2))
        return Thi
    end
    Nhi /= a0
    Dhi /= a0
    k = 2
    while k < 3 || (k < kmax && errcheck(Tmid, Thi, 10eps(real(T))))
        a0 = (α+k+1)*(β+k+1)
        a1 = T(2*k*(2*k*k-(2*α*β+α+β+1)))/(2*k-1)
        a2 = T((α+1-k)*(β+1-k)*(2*k+1))/(2*k-1)
        Nhi, Nmid, Nlo = -a0*Nhi - a1*(Nhi+Nmid) - a2*(Nmid+Nlo) + (4k+2)*ζ*(Nhi-Nmid), Nhi, Nmid
        Dhi, Dmid, Dlo = -a0*Dhi - a1*(Dhi+Dmid) - a2*(Dmid+Dlo) + (4k+2)*ζ*(Dhi-Dmid), Dhi, Dmid
        Thi, Tmid, Tlo = Nhi/Dhi, Thi, Tmid
        if norm(a0) < eps((absα+k+1)*(absβ+k+1))
            return Thi
        end
        Nhi /= a0
        Dhi /= a0
        k += 1
    end
    return isfinite(Thi) ? Thi : isfinite(Tmid) ? Tmid : Tlo
end

# ₂F₁(α,β,γ;z), algorithm γ = 2.
function pFqweniger(α::Tuple{T1, T1}, β::Tuple{T2}, z::T3; kmax::Int = 10_000) where {T1, T2, T3}
    γ = β[1]
    (α, β) = α
    T = promote_type(T1, T2, T3)
    absα = abs(T(α))
    absβ = abs(T(β))
    if norm(z) < eps(real(T)) || norm(α*β) < eps(absα*absβ)
        return one(T)
    end
    ζ = inv(z)
    Nlo = γ*ζ/(α*β)
    Dlo = γ*ζ/(α*β)
    Tlo = Nlo/Dlo
    a0 = (α+1)*(β+1)
    b0 = 2*(γ+1)
    Nmid = (b0*ζ-a0)*Nlo + b0*ζ
    Dmid = (b0*ζ-a0)*Dlo
    Tmid = Nmid/Dmid
    if norm(a0) < eps((absα+1)*(absβ+1))
        return Tmid
    end
    Nmid /= a0
    Dmid /= a0
    k = 1
    a0 = (α+2)*(β+2)
    a1 = 2*(2-(2*α*β+α+β+1))
    b0 = 6*(γ+2)
    b1 = -6*γ
    Nhi = -a0*Nmid - a1*(Nmid+Nlo) + ζ*(b0*Nmid + b1*Nlo) + b1*ζ
    Dhi = -a0*Dmid - a1*(Dmid+Dlo) + ζ*(b0*Dmid + b1*Dlo)
    Thi = Nhi/Dhi
    if norm(a0) < eps((absα+2)*(absβ+2))
        return Thi
    end
    Nhi /= a0
    Dhi /= a0
    k = 2
    while k < 3 || (k < kmax && errcheck(Tmid, Thi, 10eps(real(T))))
        a0 = (α+k+1)*(β+k+1)
        a1 = T(2*k*(2*k*k-(2*α*β+α+β+1)))/(2*k-1)
        a2 = T((α+1-k)*(β+1-k)*(2*k+1))/(2*k-1)
        b0 = (4k+2)*(γ+k+1)
        b1 = (4k+2)*(k-γ-1)
        Nhi, Nmid, Nlo = -a0*Nhi - a1*(Nhi+Nmid) - a2*(Nmid+Nlo) + ζ*(b0*Nhi + b1*Nmid), Nhi, Nmid
        Dhi, Dmid, Dlo = -a0*Dhi - a1*(Dhi+Dmid) - a2*(Dmid+Dlo) + ζ*(b0*Dhi + b1*Dmid), Dhi, Dmid
        Thi, Tmid, Tlo = Nhi/Dhi, Thi, Tmid
        if norm(a0) < eps((absα+k+1)*(absβ+k+1))
            return Thi
        end
        Nhi /= a0
        Dhi /= a0
        k += 1
    end
    return isfinite(Thi) ? Thi : isfinite(Tmid) ? Tmid : Tlo
end

# ₘFₙ(α;β;z)
# γ ∉ ℕ
function pFqweniger(α::AbstractVector{T1}, β::AbstractVector{T2}, z::T3; kwds...) where {T1, T2, T3}
    pFqweniger(Tuple(α), Tuple(β), z; kwds...)
end
function pFqweniger(α::NTuple{p, Any}, β::NTuple{q, Any}, z; kwds...) where {p, q}
    T1 = isempty(α) ? Any : mapreduce(typeof, promote_type, α)
    T2 = isempty(β) ? Any : mapreduce(typeof, promote_type, β)
    pFqweniger(T1.(α), T2.(β), z; kwds...)
end
function pFqweniger(α::NTuple{p, T1}, β::NTuple{q, T2}, z::T3; kmax::Int = 10_000) where {p, q, T1, T2, T3}
    T = promote_type(eltype(α), eltype(β), T3)
    absα = abs.(T.(α))
    if norm(z) < eps(real(T)) || norm(prod(α)) < eps(real(T)(prod(absα)))
        return one(T)
    end
    γ = T(3)/2
    ζ = inv(z)
    r = max(p, q)+3
    ρ = max(p, q)+1
    C = zeros(T, r)
    C[1] = one(T)
    Cρ = zeros(T, ρ+2)
    Cρ[ρ+2] = one(T)
    @inbounds for s in ρ:-1:0
        Cρ[s+1] = -(s+1)*Cρ[s+2]/(ρ+1-s)
    end
    C1 = zeros(T, ρ+1)
    C2 = zeros(T, ρ+2)
    C2[ρ+2] = one(T)/pochhammer(γ-ρ-2, ρ+2)
    C3 = zeros(T, ρ+2)
    C3[ρ+1] = one(T)/pochhammer(γ-ρ-1, ρ+2)
    C3[ρ+2] = one(T)/pochhammer(γ-ρ, ρ+2)
    P = zeros(T, p+2)
    t = γ
    @inbounds for j in 1:p
        t *= α[j]+1
    end
    P[1] = t
    err = abs(γ)
    @inbounds for j in 1:p
        err *= absα[j]+1
    end
    Q = zeros(T, q+2)
    t = one(T)
    @inbounds for j in 1:q
        t *= β[j]+1
    end
    Q[1] = 2t
    N = zeros(T, r+1)
    ΔN = zeros(T, r+1)
    ΔNold = zeros(T, r+1)
    D = zeros(T, r+1)
    ΔD = zeros(T, r+1)
    ΔDold = zeros(T, r+1)
    R = zeros(T, r+1)
    N[r+1] = prod(β)*ζ/prod(α)/(γ-1)
    ΔN[r] = N[r+1]/pochhammer(γ-ρ-1, ρ)
    D[r+1] = prod(β)*ζ/prod(α)/(γ-1)
    ΔD[r] = D[r+1]/pochhammer(γ-ρ-1, ρ)
    R[r+1] = N[r+1]/D[r+1]
    k = 0
    @inbounds while k < r || (k < kmax && errcheck(R[r], R[r+1], 10eps(real(T))))
        for j in 1:r
            N[j] = N[j+1]
            D[j] = D[j+1]
            R[j] = R[j+1]
        end
        t1 = zero(T)
        for j in 0:min(k, q+1)
            t1 += C[j+1]*Q[j+1]*ΔN[r-j]
        end
        if k ≤ ρ
            for j in 0:k
                t2 = one(T)
                for i in 1:q
                    t2 *= β[i]+j+1
                end
                t2 *= j+2
                t1 += C[j+1]*(-one(T))^(k-j)*pochhammer(j+γ, k-ρ-1)*t2
            end
        end
        t2 = zero(T)
        s2 = zero(T)
        for s in max(0, ρ+1-k):ρ
            s2 += Cρ[s+1]*C1[s+1]*(N[r-ρ+s]+(γ+k-ρ+s-2)*N[r-ρ+s-1])
        end
        s2 += (γ+k-1)*N[r]/pochhammer(γ+2k-ρ-1, ρ+2)
        t2 += P[1]*s2
        s2 = zero(T)
        for s in max(0, ρ+1-k):ρ+1
            s2 += Cρ[s+1]*C2[s+1]*(γ+2k-2ρ+2s-3)*N[r-ρ+s-1]
        end
        ΔNold[r+1] = s2
        for j in 1:min(k, p+1)
            t2 += C[j+1]*P[j+1]*(ΔNold[r+2-j] + ΔNold[r+1-j])
        end
        N[r+1] = ζ*t1-t2
        t1 = zero(T)
        for j in 0:min(k, q+1)
            t1 += C[j+1]*Q[j+1]*ΔD[r-j]
        end
        t2 = zero(T)
        s2 = zero(T)
        for s in max(0,ρ+1-k):ρ
            s2 += Cρ[s+1]*C1[s+1]*(D[r-ρ+s]+(γ+k-ρ+s-2)*D[r-ρ+s-1])
        end
        s2 += (γ+k-1)*D[r]/pochhammer(γ+2k-ρ-1, ρ+2)
        t2 += P[1]*s2
        s2 = zero(T)
        for s in max(0,ρ+1-k):ρ+1
            s2 += Cρ[s+1]*C2[s+1]*(γ+2k-2ρ+2s-3)*D[r-ρ+s-1]
        end
        ΔDold[r+1] = s2
        for j in 1:min(k, p+1)
            t2 += C[j+1]*P[j+1]*(ΔDold[r+2-j] + ΔDold[r+1-j])
        end
        D[r+1] = ζ*t1-t2
        if norm(P[1]) < eps(err)
            return N[r+1]/D[r+1]
        end
        N[r+1] /= P[1]/pochhammer(γ+2k-ρ-1, ρ+2)
        D[r+1] /= P[1]/pochhammer(γ+2k-ρ-1, ρ+2)
        R[r+1] = N[r+1]/D[r+1]
        s1 = zero(T)
        for s in max(0, ρ-k):ρ+1
            s1 += Cρ[s+1]*C3[s+1]*(γ+2k-2ρ+2s-1)*N[r-ρ+s]
        end
        ΔN[r+1] = s1
        s1 = zero(T)
        for s in max(0, ρ-k):ρ+1
            s1 += Cρ[s+1]*C3[s+1]*(γ+2k-2ρ+2s-1)*D[r-ρ+s]
        end
        ΔD[r+1] = s1
        k += 1
        for j in min(k, p+1):-1:0
            ΔNold[r-j] = (γ+2k-ρ-j-4)*ΔNold[r-j+1] + (k-j-1)*ΔNold[r-j]
        end
        for j in min(k, q+1):-1:0
            ΔN[r-j] = (γ+2k-ρ-j-2)*ΔN[r-j+1] + (k-j)*ΔN[r-j]
        end
        for j in min(k, p+1):-1:0
            ΔDold[r-j] = (γ+2k-ρ-j-4)*ΔDold[r-j+1] + (k-j-1)*ΔDold[r-j]
        end
        for j in min(k, q+1):-1:0
            ΔD[r-j] = (γ+2k-ρ-j-2)*ΔD[r-j+1] + (k-j)*ΔD[r-j]
        end
        for j in min(k, ρ):-1:1
            C[j+1] += C[j]
        end
        if k ≤ ρ+1
            for s in max(0, ρ+1-k):ρ
                C1[s+1] = pochhammer(k-ρ+s, ρ+1-s)/pochhammer(γ+2k-2ρ+s-2, ρ+2)
            end
            for s in max(0, ρ+1-k):ρ+1
                C2[s+1] = pochhammer(k-ρ+s, ρ+1-s)/pochhammer(γ+2k-2ρ+s-3, ρ+2)
            end
        else
            for s in 0:ρ
                C1[s+1] *= k/(k-one(T)-ρ+s)*(γ+2k-2ρ+s-4)/(γ+2k-ρ+s-2)*(γ+2k-2ρ+s-3)/(γ+2k-ρ+s-1)
            end
            for s in 0:ρ+1
                C2[s+1] *= k/(k-one(T)-ρ+s)*(γ+2k-2ρ+s-5)/(γ+2k-ρ+s-3)*(γ+2k-2ρ+s-4)/(γ+2k-ρ+s-2)
            end
        end
        if k ≤ ρ
            for s in max(0, ρ-k):ρ+1
                C3[s+1] = pochhammer(k-ρ+s+1, ρ+1-s)/pochhammer(γ+2k-2ρ+s-1, ρ+2)
            end
        else
            for s in max(0, ρ-k):ρ+1
                C3[s+1] *= (k+one(T))/(k-ρ+s)*(γ+2k-2ρ+s-3)/(γ+2k-ρ+s-1)*(γ+2k-2ρ+s-2)/(γ+2k-ρ+s)
            end
        end
        t = γ+k
        for j in 1:p
            t *= α[j]+k+1
        end
        err = abs(γ)+k
        for j in 1:p
            err *= absα[j]+k+1
        end
        for j in 2:p+2
            s = t - P[j-1]
            P[j-1] = t
            t = s
        end
        P[p+2] = t
        t = one(T)
        for j in 1:q
            t *= β[j]+k+1
        end
        t *= k+2
        for j in 2:q+2
            s = t - Q[j-1]
            Q[j-1] = t
            t = s
        end
        Q[q+2] = t
    end
    return isfinite(R[r+1]) ? R[r+1] : R[r]
end
