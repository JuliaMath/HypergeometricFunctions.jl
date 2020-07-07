# Rational approximations to generalized hypergeometric functions
# using Drummond's sequence transformation

# ₀F₀(;z)
function drummond0F0(z::T; kmax::Int = 10_000) where T
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
    Nhi, Nlo = ((k+2)*ζ-1)*Nhi + k*ζ*Nlo + ζ, Nhi
    Dhi, Dlo = ((k+2)*ζ-1)*Dhi + k*ζ*Dlo, Dhi
    Thi, Tlo = Nhi/Dhi, Thi
    k += 1
    while k < kmax && errcheck(Tlo, Thi, 10eps(real(T)))
        Nhi, Nlo = ((k+2)*ζ-1)*Nhi + k*ζ*Nlo, Nhi
        Dhi, Dlo = ((k+2)*ζ-1)*Dhi + k*ζ*Dlo, Dhi
        Thi, Tlo = Nhi/Dhi, Thi
        k += 1
    end
    return isfinite(Thi) ? Thi : Tlo
end

# ₁F₀(α;z)
function drummond1F0(α::T1, z::T2; kmax::Int = 10_000) where {T1, T2}
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
    Nhi, Nlo = ((k+2)*ζ-(α+2k+1))*Nhi + k*(ζ-1)*Nlo + ζ, Nhi
    Dhi, Dlo = ((k+2)*ζ-(α+2k+1))*Dhi + k*(ζ-1)*Dlo, Dhi
    Thi, Tlo = Nhi/Dhi, Thi
    if norm(α+k+1) < eps(absα+k+1)
        return Thi
    end
    Nhi /= α+k+1
    Dhi /= α+k+1
    k += 1
    while k < kmax && errcheck(Tlo, Thi, 10eps(real(T)))
        Nhi, Nlo = ((k+2)*ζ-(α+2k+1))*Nhi + k*(ζ-1)*Nlo, Nhi
        Dhi, Dlo = ((k+2)*ζ-(α+2k+1))*Dhi + k*(ζ-1)*Dlo, Dhi
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

# ₀F₁(β;z)
function drummond0F1(β::T1, z::T2; kmax::Int = 10_000) where {T1, T2}
    T = promote_type(T1, T2)
    if norm(z) < eps(real(T))
        return one(T)
    end
    ζ = inv(z)
    Nlo = β*ζ
    Dlo = β*ζ
    Tlo = Nlo/Dlo
    Nmid = ((β+1)*(2)*ζ - 1)*Nlo + (β+1)*(2)*ζ
    Dmid = ((β+1)*(2)*ζ - 1)*Dlo
    Tmid = Nmid/Dmid
    Nhi = ((β+2)*(3)*ζ - 1)*Nmid + (β+4)*ζ*Nlo + (β+4)*ζ
    Dhi = ((β+2)*(3)*ζ - 1)*Dmid + (β+4)*ζ*Dlo
    Thi = Nhi/Dhi
    k = 2
    Nhi, Nmid, Nlo = ((β+k+1)*(k+2)*ζ-1)*Nhi + k*(β+2k+2)*ζ*Nmid + k*(k-1)*ζ*Nlo + 2ζ, Nhi, Nmid
    Dhi, Dmid, Dlo = ((β+k+1)*(k+2)*ζ-1)*Dhi + k*(β+2k+2)*ζ*Dmid + k*(k-1)*ζ*Dlo, Dhi, Dmid
    Thi, Tmid, Tlo = Nhi/Dhi, Thi, Tmid
    k += 1
    while k < kmax && errcheck(Tmid, Thi, 10eps(real(T)))
        Nhi, Nmid, Nlo = ((β+k+1)*(k+2)*ζ-1)*Nhi + k*(β+2k+2)*ζ*Nmid + k*(k-1)*ζ*Nlo, Nhi, Nmid
        Dhi, Dmid, Dlo = ((β+k+1)*(k+2)*ζ-1)*Dhi + k*(β+2k+2)*ζ*Dmid + k*(k-1)*ζ*Dlo, Dhi, Dmid
        Thi, Tmid, Tlo = Nhi/Dhi, Thi, Tmid
        k += 1
    end
    return isfinite(Thi) ? Thi : isfinite(Tmid) ? Tmid : Tlo
end

# ₂F₀(α,β;z)
function drummond2F0(α::T1, β::T2, z::T3; kmax::Int = 10_000) where {T1, T2, T3}
    T = promote_type(T1, T2, T3)
    absα = abs(T(α))
    absβ = abs(T(β))
    if norm(z) < eps(real(T)) || norm(α*β) < eps(absα*absβ)
        return one(T)
    end
    ζ = inv(z)
    Nlo = ζ/(α*β)
    Dlo = ζ/(α*β)
    Tlo = Nlo/Dlo
    Nmid = (2ζ-(α+1)*(β+1))*Nlo + 2ζ
    Dmid = (2ζ-(α+1)*(β+1))*Dlo
    Tmid = Nmid/Dmid
    if norm((α+1)*(β+1)) < eps((absα+1)*(absβ+1))
        return Tmid
    end
    Nmid /= (α+1)*(β+1)
    Dmid /= (α+1)*(β+1)
    Nhi = (3ζ-(α+2)*(β+2)-(α+β+3))*Nmid - (α+β+3-ζ)*Nlo + ζ
    Dhi = (3ζ-(α+2)*(β+2)-(α+β+3))*Dmid - (α+β+3-ζ)*Dlo
    Thi = Nhi/Dhi
    if norm((α+2)*(β+2)) < eps((absα+2)*(absβ+2))
        return Thi
    end
    Nhi /= (α+2)*(β+2)
    Dhi /= (α+2)*(β+2)
    k = 2
    while k < kmax && errcheck(Tmid, Thi, 10eps(real(T)))
        Nhi, Nmid, Nlo = ((k+2)*ζ-(α+k+1)*(β+k+1)-k*(α+β+2k+1))*Nhi - k*(α+β+3k-ζ)*Nmid - k*(k-1)*Nlo, Nhi, Nmid
        Dhi, Dmid, Dlo = ((k+2)*ζ-(α+k+1)*(β+k+1)-k*(α+β+2k+1))*Dhi - k*(α+β+3k-ζ)*Dmid - k*(k-1)*Dlo, Dhi, Dmid
        Thi, Tmid, Tlo = Nhi/Dhi, Thi, Tmid
        if norm((α+k+1)*(β+k+1)) < eps((absα+k+1)*(absβ+k+1))
            return Thi
        end
        Nhi /= (α+k+1)*(β+k+1)
        Dhi /= (α+k+1)*(β+k+1)
        k += 1
    end
    return isfinite(Thi) ? Thi : isfinite(Tmid) ? Tmid : Tlo
end

# ₁F₁(α,β;z)
function drummond1F1(α::T1, β::T2, z::T3; kmax::Int = 10_000) where {T1, T2, T3}
    T = promote_type(T1, T2, T3)
    absα = abs(T(α))
    if norm(z) < eps(real(T)) || norm(α) < eps(absα)
        return one(T)
    end
    ζ = inv(z)
    Nlo = β*ζ/α
    Dlo = β*ζ/α
    Tlo = Nlo/Dlo
    Nmid = ((β+1)*(2)*ζ - (α+1))*Nlo + (β+1)*(2)*ζ
    Dmid = ((β+1)*(2)*ζ - (α+1))*Dlo
    Tmid = Nmid/Dmid
    if norm(α+1) < eps(absα+1)
        return Tmid
    end
    Nmid /= α+1
    Dmid /= α+1
    Nhi = ((β+2)*(3)*ζ - (α+3))*Nmid + ((β+4)*ζ-1)*Nlo + (β+4)*ζ
    Dhi = ((β+2)*(3)*ζ - (α+3))*Dmid + ((β+4)*ζ-1)*Dlo
    Thi = Nhi/Dhi
    if norm(α+2) < eps(absα+2)
        return Thi
    end
    Nhi /= α+2
    Dhi /= α+2
    k = 2
    Nhi, Nmid, Nlo = ((β+k+1)*(k+2)*ζ-(α+2k+1))*Nhi + k*((β+2k+2)*ζ-1)*Nmid + k*(k-1)*ζ*Nlo + 2ζ, Nhi, Nmid
    Dhi, Dmid, Dlo = ((β+k+1)*(k+2)*ζ-(α+2k+1))*Dhi + k*((β+2k+2)*ζ-1)*Dmid + k*(k-1)*ζ*Dlo, Dhi, Dmid
    Thi, Tmid, Tlo = Nhi/Dhi, Thi, Tmid
    if norm(α+k+1) < eps(absα+k+1)
        return Thi
    end
    Nhi /= α+k+1
    Dhi /= α+k+1
    k += 1
    while k < kmax && errcheck(Tmid, Thi, 10eps(real(T)))
        Nhi, Nmid, Nlo = ((β+k+1)*(k+2)*ζ-(α+2k+1))*Nhi + k*((β+2k+2)*ζ-1)*Nmid + k*(k-1)*ζ*Nlo, Nhi, Nmid
        Dhi, Dmid, Dlo = ((β+k+1)*(k+2)*ζ-(α+2k+1))*Dhi + k*((β+2k+2)*ζ-1)*Dmid + k*(k-1)*ζ*Dlo, Dhi, Dmid
        Thi, Tmid, Tlo = Nhi/Dhi, Thi, Tmid
        if norm(α+k+1) < eps(absα+k+1)
            return Thi
        end
        Nhi /= α+k+1
        Dhi /= α+k+1
        k += 1
    end
    return isfinite(Thi) ? Thi : isfinite(Tmid) ? Tmid : Tlo
end

# ₀F₂(α,β;z)
function drummond0F2(α::T1, β::T2, z::T3; kmax::Int = 10_000) where {T1, T2, T3}
    T = promote_type(T1, T2, T3)
    if norm(z) < eps(real(T)) || norm(α) < eps(real(T)) || norm(β) < eps(real(T))
        return one(T)
    end
    ζ = inv(z)
    Nlo = α*β*ζ
    Dlo = α*β*ζ
    Tlo = Nlo/Dlo
    k = 0
    Nmid2 = (ζ*(k+2)*(α+k+1)*(β+k+1)-1)*Nlo + (2*(α+1)*(β+1))*ζ
    Dmid2 = (ζ*(k+2)*(α+k+1)*(β+k+1)-1)*Dlo
    Tmid2 = Nmid2/Dmid2
    k = 1
    Nmid1 = (ζ*(k+2)*(α+k+1)*(β+k+1)-1)*Nmid2 + ζ*k*((k+1)*(α+β+2k)+(α+k)*(β+k)+α+β+3k+2)*Nlo + (3*(α+β+3)+(α+1)*(β+1))*ζ
    Dmid1 = (ζ*(k+2)*(α+k+1)*(β+k+1)-1)*Dmid2 + ζ*k*((k+1)*(α+β+2k)+(α+k)*(β+k)+α+β+3k+2)*Dlo
    Tmid1 = Nmid1/Dmid1
    k = 2
    Nhi = (ζ*(k+2)*(α+k+1)*(β+k+1)-1)*Nmid1 + ζ*k*((k+1)*(α+β+2k)+(α+k)*(β+k)+α+β+3k+2)*Nmid2 + ζ*k*(k-1)*(3k+α+β+1)*Nlo + (14+2*(α+β))*ζ
    Dhi = (ζ*(k+2)*(α+k+1)*(β+k+1)-1)*Dmid1 + ζ*k*((k+1)*(α+β+2k)+(α+k)*(β+k)+α+β+3k+2)*Dmid2 + ζ*k*(k-1)*(3k+α+β+1)*Dlo
    Thi = Nhi/Dhi
    k = 3
    Nhi, Nmid1, Nmid2, Nlo = (ζ*(k+2)*(α+k+1)*(β+k+1)-1)*Nhi + ζ*k*((k+1)*(α+β+2k)+(α+k)*(β+k)+α+β+3k+2)*Nmid1 + ζ*k*(k-1)*(3k+α+β+1)*Nmid2 + ζ*k*(k-1)*(k-2)*Nlo + 6ζ, Nhi, Nmid1, Nmid2
    Dhi, Dmid1, Dmid2, Dlo = (ζ*(k+2)*(α+k+1)*(β+k+1)-1)*Dhi + ζ*k*((k+1)*(α+β+2k)+(α+k)*(β+k)+α+β+3k+2)*Dmid1 + ζ*k*(k-1)*(3k+α+β+1)*Dmid2 + ζ*k*(k-1)*(k-2)*Dlo, Dhi, Dmid1, Dmid2
    Thi, Tmid1, Tmid2, Tlo = Nhi/Dhi, Thi, Tmid1, Tmid2
    k += 1
    while k < kmax && errcheck(Tmid1, Thi, 10eps(real(T)))
        Nhi, Nmid1, Nmid2, Nlo = (ζ*(k+2)*(α+k+1)*(β+k+1)-1)*Nhi + ζ*k*((k+1)*(α+β+2k)+(α+k)*(β+k)+α+β+3k+2)*Nmid1 + ζ*k*(k-1)*(3k+α+β+1)*Nmid2 + ζ*k*(k-1)*(k-2)*Nlo, Nhi, Nmid1, Nmid2
        Dhi, Dmid1, Dmid2, Dlo = (ζ*(k+2)*(α+k+1)*(β+k+1)-1)*Dhi + ζ*k*((k+1)*(α+β+2k)+(α+k)*(β+k)+α+β+3k+2)*Dmid1 + ζ*k*(k-1)*(3k+α+β+1)*Dmid2 + ζ*k*(k-1)*(k-2)*Dlo, Dhi, Dmid1, Dmid2
        Thi, Tmid1, Tmid2, Tlo = Nhi/Dhi, Thi, Tmid1, Tmid2
        k += 1
    end
    return isfinite(Thi) ? Thi : isfinite(Tmid1) ? Tmid1 : isfinite(Tmid2) ? Tmid2 : Tlo
end

# ₂F₁(α,β,γ;z)
function drummond2F1(α::T1, β::T2, γ::T3, z::T4; kmax::Int = 10_000) where {T1, T2, T3, T4}
    T = promote_type(T1, T2, T3, T4)
    absα = abs(T(α))
    absβ = abs(T(β))
    if norm(z) < eps(real(T)) || norm(α*β) < eps(absα*absβ)
        return one(T)
    end
    ζ = inv(z)
    Nlo = γ*ζ/(α*β)
    Dlo = γ*ζ/(α*β)
    Tlo = Nlo/Dlo
    Nmid = ((2)*(γ+1)*ζ-(α+1)*(β+1))*Nlo + 2*(γ+1)*ζ
    Dmid = ((2)*(γ+1)*ζ-(α+1)*(β+1))*Dlo
    Tmid = Nmid/Dmid
    if norm((α+1)*(β+1)) < eps((absα+1)*(absβ+1))
        return Tmid
    end
    Nmid /= (α+1)*(β+1)
    Dmid /= (α+1)*(β+1)
    Nhi = ((3)*(γ+2)*ζ-(α+2)*(β+2)-(α+β+3))*Nmid + ((γ+4)*ζ-(α+β+3))*Nlo + (γ+4)*ζ
    Dhi = ((3)*(γ+2)*ζ-(α+2)*(β+2)-(α+β+3))*Dmid + ((γ+4)*ζ-(α+β+3))*Dlo
    Thi = Nhi/Dhi
    if norm((α+2)*(β+2)) < eps((absα+2)*(absβ+2))
        return Thi
    end
    Nhi /= (α+2)*(β+2)
    Dhi /= (α+2)*(β+2)
    k = 2
    Nhi, Nmid, Nlo = ((k+2)*(γ+k+1)*ζ-(α+k+1)*(β+k+1)-k*(α+β+2k+1))*Nhi + k*((γ+2k+2)*ζ-(α+β+3k))*Nmid + k*(k-1)*(ζ-1)*Nlo + 2ζ, Nhi, Nmid
    Dhi, Dmid, Dlo = ((k+2)*(γ+k+1)*ζ-(α+k+1)*(β+k+1)-k*(α+β+2k+1))*Dhi + k*((γ+2k+2)*ζ-(α+β+3k))*Dmid + k*(k-1)*(ζ-1)*Dlo, Dhi, Dmid
    Thi, Tmid, Tlo = Nhi/Dhi, Thi, Tmid
    if norm((α+k+1)*(β+k+1)) < eps((absα+k+1)*(absβ+k+1))
        return Thi
    end
    Nhi /= (α+k+1)*(β+k+1)
    Dhi /= (α+k+1)*(β+k+1)
    k += 1
    while k < kmax && errcheck(Tmid, Thi, 10eps(real(T)))
        Nhi, Nmid, Nlo = ((k+2)*(γ+k+1)*ζ-(α+k+1)*(β+k+1)-k*(α+β+2k+1))*Nhi + k*((γ+2k+2)*ζ-(α+β+3k))*Nmid + k*(k-1)*(ζ-1)*Nlo, Nhi, Nmid
        Dhi, Dmid, Dlo = ((k+2)*(γ+k+1)*ζ-(α+k+1)*(β+k+1)-k*(α+β+2k+1))*Dhi + k*((γ+2k+2)*ζ-(α+β+3k))*Dmid + k*(k-1)*(ζ-1)*Dlo, Dhi, Dmid
        Thi, Tmid, Tlo = Nhi/Dhi, Thi, Tmid
        if norm((α+k+1)*(β+k+1)) < eps((absα+k+1)*(absβ+k+1))
            return Thi
        end
        Nhi /= (α+k+1)*(β+k+1)
        Dhi /= (α+k+1)*(β+k+1)
        k += 1
    end
    return isfinite(Thi) ? Thi : isfinite(Tmid) ? Tmid : Tlo
end

# ₘFₙ(α;β;z)
function pFqdrummond(α::AbstractVector{T1}, β::AbstractVector{T2}, z::T3; kmax::Int = 10_000) where {T1, T2, T3}
    T = promote_type(T1, T2, T3)
    absα = abs.(T.(α))
    if norm(z) < eps(real(T)) || norm(prod(α)) < eps(prod(absα))
        return one(T)
    end
    ζ = inv(z)
    p = length(α)
    q = length(β)
    r = max(p+1, q+2)
    C = zeros(T, r)
    C[1] = one(T)
    Ĉ = zeros(T, r)
    Ĉ[1] = one(T)
    P = zeros(T, p+1)
    t = one(T)
    for j in 1:p
        t *= α[j]+1
    end
    P[1] = t
    err = one(real(T))
    for j in 1:p
        err *= absα[j]+1
    end
    Q = zeros(T, q+2)
    t = one(T)
    for j in 1:q
        t *= β[j]+1
    end
    Q[1] = 2t
    N = zeros(T, r+1)
    D = zeros(T, r+1)
    R = zeros(T, r+1)
    N[r+1] = prod(β)*ζ/prod(α)
    D[r+1] = prod(β)*ζ/prod(α)
    R[r+1] = N[r+1]/D[r+1]
    k = 0
    while k < r || (k < kmax && errcheck(R[r], R[r+1], 10eps(real(T))))
        for j in 1:r
            N[j] = N[j+1]
            D[j] = D[j+1]
            R[j] = R[j+1]
        end
        t1 = zero(T)
        for j in 0:min(k, q+1)
            t1 += Ĉ[j+1]*Q[j+1]*N[r-j]
        end
        if k ≤ q+1
            if p > q
                t1 += Q[k+1]
            else
                t1 += Q[k+1] / T(factorial(k+1))
            end
        end
        t2 = zero(T)
        t2 += Ĉ[1]*P[1]*N[r]
        for j in 1:min(k, p)
            t2 += P[j+1]*(C[j+1]*N[r-j+1] + Ĉ[j+1]*N[r-j])
        end
        N[r+1] = ζ*t1-t2
        t1 = zero(T)
        for j in 0:min(k, q+1)
            t1 += Ĉ[j+1]*Q[j+1]*D[r-j]
        end
        t2 = zero(T)
        t2 += Ĉ[1]*P[1]*D[r]
        for j in 1:min(k, p)
            t2 += P[j+1]*(C[j+1]*D[r-j+1] + Ĉ[j+1]*D[r-j])
        end
        D[r+1] = ζ*t1-t2
        R[r+1] = N[r+1]/D[r+1]
        if norm(P[1]) < eps(err)
            return R[r+1]
        end
        N[r+1] /= P[1]
        D[r+1] /= P[1]
        k += 1
        if p > q
            for j in min(k, max(p, q+1)):-1:1
                C[j+1] += C[j]
                Ĉ[j+1] += Ĉ[j]
            end
        else
            for j in min(k, max(p, q+1)):-1:1
                C[j+1] = (C[j+1]*(k+1-j) + C[j])/(k+1)
                Ĉ[j+1] = (Ĉ[j+1]*(k-j) + Ĉ[j])/(k+1)
            end
            Ĉ[1] = Ĉ[1]*k/(k+1)
        end
        t = one(T)
        for j in 1:p
            t *= α[j]+k+1
        end
        err = one(real(T))
        for j in 1:p
            err *= absα[j]+k+1
        end
        for j in 2:p+1
            s = t - P[j-1]
            P[j-1] = t
            t = s
        end
        P[p+1] = t
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
