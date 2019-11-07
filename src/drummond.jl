# Rational approximations to generalized hypergeometric functions
# using Drummond's sequence transformation

# ₀F₀(;z)
function drummond0F0(z::T) where T
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
    while abs(Thi-Tlo) > 10*abs(Thi)*eps(real(T)) && k < 10_000
        Nhi, Nlo = ((k+2)*ζ-1)*Nhi + k*ζ*Nlo, Nhi
        Dhi, Dlo = ((k+2)*ζ-1)*Dhi + k*ζ*Dlo, Dhi
        Thi, Tlo = Nhi/Dhi, Thi
        k += 1
    end
    return isnan(Thi) ? Tlo : Thi
end

# ₁F₀(α;z)
function drummond1F0(α::T1, z::T2) where {T1, T2}
    T = promote_type(T1, T2)
    if norm(z) < eps(real(T)) || norm(α) < eps(real(T))
        return one(T)
    end
    ζ = inv(z)
    Nlo = ζ/α
    Dlo = ζ/α
    Tlo = Nlo/Dlo
    if norm(α+1) < eps(real(T))
        return Tlo
    end
    Nhi = ((2ζ - (α+1))*Nlo + 2ζ)/(α+1)
    Dhi = ((2ζ - (α+1))*Dlo)/(α+1)
    Thi = Nhi/Dhi
    k = 1
    Nhi, Nlo = (((k+2)*ζ-(α+2k+1))*Nhi + k*(ζ-1)*Nlo + ζ)/(α+k+1), Nhi
    Dhi, Dlo = (((k+2)*ζ-(α+2k+1))*Dhi + k*(ζ-1)*Dlo)/(α+k+1), Dhi
    Thi, Tlo = Nhi/Dhi, Thi
    k += 1
    while abs(Thi-Tlo) > 10*abs(Thi)*eps(real(T)) && k < 10_000
        Nhi, Nlo = (((k+2)*ζ-(α+2k+1))*Nhi + k*(ζ-1)*Nlo)/(α+k+1), Nhi
        Dhi, Dlo = (((k+2)*ζ-(α+2k+1))*Dhi + k*(ζ-1)*Dlo)/(α+k+1), Dhi
        Thi, Tlo = Nhi/Dhi, Thi
        k += 1
    end
    return isnan(Thi) ? Tlo : Thi
end

# ₀F₁(β;z)
function drummond0F1(β::T1, z::T2) where {T1, T2}
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
    while abs(Thi-Tmid) > 10*abs(Thi)*eps(real(T)) && abs(Tmid-Tlo) > 10*abs(Tmid)*eps(real(T)) && k < 10_000
        Nhi, Nmid, Nlo = ((β+k+1)*(k+2)*ζ-1)*Nhi + k*(β+2k+2)*ζ*Nmid + k*(k-1)*ζ*Nlo, Nhi, Nmid
        Dhi, Dmid, Dlo = ((β+k+1)*(k+2)*ζ-1)*Dhi + k*(β+2k+2)*ζ*Dmid + k*(k-1)*ζ*Dlo, Dhi, Dmid
        Thi, Tmid, Tlo = Nhi/Dhi, Thi, Tmid
        k += 1
    end
    return isnan(Thi) ? isnan(Tmid) ? Tlo : Tmid : Thi
end

# ₂F₀(α,β;z)
function drummond2F0(α::T1, β::T2, z::T3) where {T1, T2, T3}
    T = promote_type(T1, T2, T3)
    if norm(z) < eps(real(T)) || norm(α*β) < eps(real(T))
        return one(T)
    end
    ζ = inv(z)
    Nlo = ζ/(α*β)
    Dlo = ζ/(α*β)
    Tlo = Nlo/Dlo
    if norm((α+1)*(β+1)) < eps(real(T))
        return Tlo
    end
    Nmid = ((2ζ-(α+1)*(β+1))*Nlo + 2ζ)/((α+1)*(β+1))
    Dmid = (2ζ-(α+1)*(β+1))*Dlo/((α+1)*(β+1))
    Tmid = Nmid/Dmid
    if norm((α+2)*(β+2)) < eps(real(T))
        return Tmid
    end
    Nhi = ((3ζ-(α+2)*(β+2)-(α+β+3))*Nmid - (α+β+3-ζ)*Nlo + ζ)/((α+2)*(β+2))
    Dhi = ((3ζ-(α+2)*(β+2)-(α+β+3))*Dmid - (α+β+3-ζ)*Dlo)/((α+2)*(β+2))
    Thi = Nhi/Dhi
    if norm((α+3)*(β+3)) < eps(real(T))
        return Thi
    end
    k = 2
    while abs(Thi-Tmid) > 10*abs(Thi)*eps(real(T)) && abs(Tmid-Tlo) > 10*abs(Tmid)*eps(real(T)) && k < 10_000
        Nhi, Nmid, Nlo = (((k+2)*ζ-(α+k+1)*(β+k+1)-k*(α+β+2k+1))*Nhi - k*(α+β+3k-ζ)*Nmid - k*(k-1)*Nlo)/((α+k+1)*(β+k+1)), Nhi, Nmid
        Dhi, Dmid, Dlo = (((k+2)*ζ-(α+k+1)*(β+k+1)-k*(α+β+2k+1))*Dhi - k*(α+β+3k-ζ)*Dmid - k*(k-1)*Dlo)/((α+k+1)*(β+k+1)), Dhi, Dmid
        Thi, Tmid, Tlo = Nhi/Dhi, Thi, Tmid
        k += 1
    end
    return isnan(Thi) ? isnan(Tmid) ? Tlo : Tmid : Thi
end

# ₁F₁(α,β;z)
function drummond1F1(α::T1, β::T2, z::T3) where {T1, T2, T3}
    T = promote_type(T1, T2, T3)
    if norm(z) < eps(real(T)) || norm(α) < eps(real(T))
        return one(T)
    end
    ζ = inv(z)
    Nlo = β*ζ/α
    Dlo = β*ζ/α
    Tlo = Nlo/Dlo
    if norm(α+1) < eps(real(T))
        return Tlo
    end
    Nmid = (((β+1)*(2)*ζ - (α+1))*Nlo + (β+1)*(2)*ζ)/(α+1)
    Dmid = ((β+1)*(2)*ζ - (α+1))*Dlo/(α+1)
    Tmid = Nmid/Dmid
    if norm(α+2) < eps(real(T))
        return Tmid
    end
    Nhi = (((β+2)*(3)*ζ - (α+3))*Nmid + ((β+4)*ζ-1)*Nlo + (β+4)*ζ)/(α+2)
    Dhi = (((β+2)*(3)*ζ - (α+3))*Dmid + ((β+4)*ζ-1)*Dlo)/(α+2)
    Thi = Nhi/Dhi
    if norm(α+3) < eps(real(T))
        return Thi
    end
    k = 2
    Nhi, Nmid, Nlo = (((β+k+1)*(k+2)*ζ-(α+2k+1))*Nhi + k*((β+2k+2)*ζ-1)*Nmid + k*(k-1)*ζ*Nlo + 2ζ)/(α+k+1), Nhi, Nmid
    Dhi, Dmid, Dlo = (((β+k+1)*(k+2)*ζ-(α+2k+1))*Dhi + k*((β+2k+2)*ζ-1)*Dmid + k*(k-1)*ζ*Dlo)/(α+k+1), Dhi, Dmid
    Thi, Tmid, Tlo = Nhi/Dhi, Thi, Tmid
    k += 1
    while abs(Thi-Tmid) > 10*abs(Thi)*eps(real(T)) && abs(Tmid-Tlo) > 10*abs(Tmid)*eps(real(T)) && k < 10_000
        Nhi, Nmid, Nlo = (((β+k+1)*(k+2)*ζ-(α+2k+1))*Nhi + k*((β+2k+2)*ζ-1)*Nmid + k*(k-1)*ζ*Nlo)/(α+k+1), Nhi, Nmid
        Dhi, Dmid, Dlo = (((β+k+1)*(k+2)*ζ-(α+2k+1))*Dhi + k*((β+2k+2)*ζ-1)*Dmid + k*(k-1)*ζ*Dlo)/(α+k+1), Dhi, Dmid
        Thi, Tmid, Tlo = Nhi/Dhi, Thi, Tmid
        k += 1
    end
    return isnan(Thi) ? isnan(Tmid) ? Tlo : Tmid : Thi
end

# ₀F₂(α,β;z)
function drummond0F2(α::T1, β::T2, z::T3) where {T1, T2, T3}
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
    while abs(Thi-Tmid1) > 10*abs(Thi)*eps(real(T)) && abs(Tmid1-Tmid2) > 10*abs(Tmid1)*eps(real(T)) && abs(Tmid2-Tlo) > 10*abs(Tmid2)*eps(real(T)) && k < 10_000
        Nhi, Nmid1, Nmid2, Nlo = (ζ*(k+2)*(α+k+1)*(β+k+1)-1)*Nhi + ζ*k*((k+1)*(α+β+2k)+(α+k)*(β+k)+α+β+3k+2)*Nmid1 + ζ*k*(k-1)*(3k+α+β+1)*Nmid2 + ζ*k*(k-1)*(k-2)*Nlo, Nhi, Nmid1, Nmid2
        Dhi, Dmid1, Dmid2, Dlo = (ζ*(k+2)*(α+k+1)*(β+k+1)-1)*Dhi + ζ*k*((k+1)*(α+β+2k)+(α+k)*(β+k)+α+β+3k+2)*Dmid1 + ζ*k*(k-1)*(3k+α+β+1)*Dmid2 + ζ*k*(k-1)*(k-2)*Dlo, Dhi, Dmid1, Dmid2
        Thi, Tmid1, Tmid2, Tlo = Nhi/Dhi, Thi, Tmid1, Tmid2
        k += 1
    end
    return isnan(Thi) ? isnan(Tmid1) ? isnan(Tmid2) ? Tlo : Tmid1 : Tmid2 : Thi
end

# ₂F₁(α,β,γ;z)
function drummond2F1(α::T1, β::T2, γ::T3, z::T4) where {T1, T2, T3, T4}
    T = promote_type(T1, T2, T3, T4)
    if norm(z) < eps(real(T)) || norm(α*β) < eps(real(T))
        return one(T)
    end
    ζ = inv(z)
    Nlo = γ*ζ/(α*β)
    Dlo = γ*ζ/(α*β)
    Tlo = Nlo/Dlo
    if norm((α+1)*(β+1)) < eps(real(T))
        return Tlo
    end
    Nmid = (((2)*(γ+1)*ζ-(α+1)*(β+1))*Nlo + 2*(γ+1)*ζ)/((α+1)*(β+1))
    Dmid = ((2)*(γ+1)*ζ-(α+1)*(β+1))*Dlo/((α+1)*(β+1))
    Tmid = Nmid/Dmid
    if norm((α+2)*(β+2)) < eps(real(T))
        return Tmid
    end
    Nhi = (((3)*(γ+2)*ζ-(α+2)*(β+2)-(α+β+3))*Nmid + ((γ+4)*ζ-(α+β+3))*Nlo + (γ+4)*ζ)/((α+2)*(β+2))
    Dhi = (((3)*(γ+2)*ζ-(α+2)*(β+2)-(α+β+3))*Dmid + ((γ+4)*ζ-(α+β+3))*Dlo)/((α+2)*(β+2))
    Thi = Nhi/Dhi
    if norm((α+3)*(β+3)) < eps(real(T))
        return Thi
    end
    k = 2
    Nhi, Nmid, Nlo = (((k+2)*(γ+k+1)*ζ-(α+k+1)*(β+k+1)-k*(α+β+2k+1))*Nhi + k*((γ+2k+2)*ζ-(α+β+3k))*Nmid + k*(k-1)*(ζ-1)*Nlo + 2ζ)/((α+k+1)*(β+k+1)), Nhi, Nmid
    Dhi, Dmid, Dlo = (((k+2)*(γ+k+1)*ζ-(α+k+1)*(β+k+1)-k*(α+β+2k+1))*Dhi + k*((γ+2k+2)*ζ-(α+β+3k))*Dmid + k*(k-1)*(ζ-1)*Dlo)/((α+k+1)*(β+k+1)), Dhi, Dmid
    Thi, Tmid, Tlo = Nhi/Dhi, Thi, Tmid
    k += 1
    while abs(Thi-Tmid) > 10*abs(Thi)*eps(real(T)) && abs(Tmid-Tlo) > 10*abs(Tmid)*eps(real(T)) && k < 10_000
        Nhi, Nmid, Nlo = (((k+2)*(γ+k+1)*ζ-(α+k+1)*(β+k+1)-k*(α+β+2k+1))*Nhi + k*((γ+2k+2)*ζ-(α+β+3k))*Nmid + k*(k-1)*(ζ-1)*Nlo)/((α+k+1)*(β+k+1)), Nhi, Nmid
        Dhi, Dmid, Dlo = (((k+2)*(γ+k+1)*ζ-(α+k+1)*(β+k+1)-k*(α+β+2k+1))*Dhi + k*((γ+2k+2)*ζ-(α+β+3k))*Dmid + k*(k-1)*(ζ-1)*Dlo)/((α+k+1)*(β+k+1)), Dhi, Dmid
        Thi, Tmid, Tlo = Nhi/Dhi, Thi, Tmid
        k += 1
    end
    return isnan(Thi) ? isnan(Tmid) ? Tlo : Tmid : Thi
end
