"""
Implement `p(z, ::Val{m}) = 1-(1-z)^(1/m)` without subtractive cancellation near the origin.
"""
function p(z)
    return p(z, Val(4))
end
function p(z, ::Val{1})
    return z
end
function p(z, ::Val{2})
    t = sqrt(1-z)
    return z/(1+t)
end
function p(z, ::Val{3})
    t = (one(z)-z)^(one(z)/3)
    u = one(z)+t
    return z/muladd(t, u, one(z))
end
function p(z, ::Val{4})
    t = (one(z)-z)^(one(z)/4)
    u = one(z)+t
    return z/muladd(t, muladd(t, u, one(z)), one(z))
end
@generated function p(z, ::Val{m}) where m
    str = """
    begin
        t = (one(z) - z)^(one(z) / m)
        u = one(t) + t
        s = u
        for k in 3:m
            s = muladd(t, s, one(z))
        end
        return z / s
    end"""
    Meta.parse(str)
end

"""
Implement `r(z, ::Val{m}) = (1-(1-z)^(1/m))/(1+(1-z)^(1/m))` without subtractive cancellation near the origin.
"""
function r(z)
    return r(z, Val(4))
end
function r(z, ::Val{1})
    return z/(2-z)
end
function r(z, ::Val{2})
    t = sqrt(1-z)
    return z/(1+t)^2
end
function r(z, ::Val{3})
    t = (one(z)-z)^(one(z)/3)
    u = one(z)+t
    return z/(u*muladd(t, u, one(z)))
end
function r(z, ::Val{4})
    t = (one(z)-z)^(one(z)/4)
    u = one(z)+t
    return z/(u*muladd(t, muladd(t, u, one(z)), one(z)))
end
@generated function r(z, ::Val{m}) where m
    str = """
    begin
        t = (one(z) - z)^(one(z) / m)
        u = one(t) + t
        s = u
        for k in 3:m
            s = muladd(t, s, one(z))
        end
        return z / (u * s)
    end"""
    Meta.parse(str)
end

"""
Maclaurin series of ₚ₊₁Fₚ(α, β; z = 1 - (1-w)^m)
"""
pFqconformalpolynomial(α, β, z; kwds...) = pFqconformalpolynomial(α, β, z, Val(4); kwds...)
function pFqconformalpolynomial(α::Tuple{Any, Vararg{Any, q}}, β::NTuple{q, Any}, z, m; kwds...) where q
    T1 = isempty(α) ? Any : mapreduce(typeof, promote_type, α)
    T2 = isempty(β) ? Any : mapreduce(typeof, promote_type, β)
    pFqconformalpolynomial(T1.(α), T2.(β), z, m; kwds...)
end


function pFqconformalpolynomial(α::Tuple{T1, T1}, β::Tuple{T2}, z::T3, ::Val{2}; kmax::Int = KMAX) where {T1, T2, T3}
    T = promote_type(T1, T2, T3)
    γ = T(β[1])
    (α, β) = T.(α)
    if abs(z) < eps(real(T))
        return one(T)
    end
    w = p(T(z), Val(2))
    B = 1 + 2*(α+β)
    C = 4*α*β
    k = 0
    ulo = zero(T)
    uhi = one(T)
    Slo = one(T)
    uhi, ulo = (w/(2*γ))*C*uhi, uhi
    Shi = Slo + uhi
    k = 1
    while (k < kmax && errcheck(Shi, Slo, 8eps(real(T))))
        alo = -(k+B-2)*(k-1) - C
        ahi = (2*B+3*k-3)*k + C
        uhi, ulo = (w/(2*(k+1)*(k+γ)))*(ahi*uhi+w*alo*ulo), uhi
        Shi, Slo = Shi + uhi, Shi
        k += 1
    end
    k < kmax || @warn "Iteration limit reached before convergence."
    return Shi
end

function pFqconformalpolynomial(α::Tuple{T1, T1}, β::Tuple{T2}, z::T3, ::Val{3}; kmax::Int = KMAX) where {T1, T2, T3}
    T = promote_type(T1, T2, T3)
    γ = T(β[1])
    (α, β) = T.(α)
    if abs(z) < eps(real(T))
        return one(T)
    end
    w = p(T(z), Val(3))
    B = 1 + 3*(α+β)
    C = 9*α*β
    k = 0
    ulo = zero(T)
    umid = zero(T)
    uhi = one(T)
    Slo = one(T)
    uhi, umid, ulo = (w/(3*γ))*C*uhi, uhi, umid
    Shi = Slo + uhi
    k = 1
    while (k < kmax && errcheck(Shi, Slo, 8eps(real(T))))
        alo = (k+B-3)*(k-2)+C
        amid = -(4*k+3*B-8)*(k-1) - 2*C
        ahi = (6*k+3*B-6)*k + C
        uhi, umid, ulo = (w/(3*(k+1)*(k+γ)))*(ahi*uhi+w*(amid*umid+w*alo*ulo)), uhi, umid
        Shi, Slo = Shi + uhi, Shi
        k += 1
    end
    k < kmax || @warn "Iteration limit reached before convergence."
    return Shi
end

function pFqconformalpolynomial(α::Tuple{T1, T1}, β::Tuple{T2}, z::T3, ::Val{4}; kmax::Int = KMAX) where {T1, T2, T3}
    T = promote_type(T1, T2, T3)
    γ = T(β[1])
    (α, β) = T.(α)
    if abs(z) < eps(real(T))
        return one(T)
    end
    w = p(T(z), Val(4))
    B = 1 + 4*(α+β)
    C = 16*α*β
    k = 0
    ulo = zero(T)
    umid2 = zero(T)
    umid1 = zero(T)
    uhi = one(T)
    Slo = one(T)
    uhi, umid1, umid2, ulo = (w/(4*γ))*C*uhi, uhi, umid1, umid2
    Shi = Slo + uhi
    k = 1
    while (k < kmax && errcheck(Shi, Slo, 8eps(real(T))))
        alo = -(k+B-4)*(k-3) - C
        amid2 = (5*k+4*B-15)*(k-2) + 3*C
        amid1 = -(10*k+6*B-20)*(k-1) - 3*C
        ahi = (10*k+4*B-10)*k + C
        uhi, umid1, umid2, ulo = (w/(4*(k+1)*(k+γ)))*(ahi*uhi+w*(amid1*umid1+w*(amid2*umid2+w*alo*ulo))), uhi, umid1, umid2
        Shi, Slo = Shi + uhi, Shi
        k += 1
    end
    k < kmax || @warn "Iteration limit reached before convergence."
    return Shi
end

"""
Maclaurin series of ₚ₊₁Fₚ(α, β; z = 1 - ((1-w)/(1+w))^m)
"""
pFqconformalrational(α, β, z; kwds...) = pFqconformalrational(α, β, z, Val(4); kwds...)
function pFqconformalrational(α::Tuple{Any, Vararg{Any, q}}, β::NTuple{q, Any}, z, m; kwds...) where q
    T1 = isempty(α) ? Any : mapreduce(typeof, promote_type, α)
    T2 = isempty(β) ? Any : mapreduce(typeof, promote_type, β)
    pFqconformalrational(T1.(α), T2.(β), z, m; kwds...)
end

function pFqconformalrational(α::Tuple{T1, T1}, β::Tuple{T2}, z::T3, ::Val{2}; kmax::Int = KMAX) where {T1, T2, T3}
    T = promote_type(T1, T2, T3)
    γ = T(β[1])
    (α, β) = T.(α)
    if abs(z) < eps(real(T))
        return one(T)
    end
    w = r(T(z), Val(2))
    k = 0
    ulo = zero(T)
    umid = zero(T)
    uhi = one(T)
    Slo = one(T)
    uhi, umid, ulo = (w/γ)*4*α*β*uhi, uhi, umid
    Shi = Slo + uhi
    k = 1
    while (k < kmax && errcheck(Shi, Slo, 8eps(real(T))))
        alo = T(-(k-2)*(γ-k+1))
        amid = T(-((k-1)*(3γ-4α-4β-k)+4α*β))
        ahi = T(-(k*(k-1+3γ-4α-4β)-4α*β))
        uhi, umid, ulo = (w/((k+1)*(k+γ)))*(ahi*uhi+w*(amid*umid+w*alo*ulo)), uhi, umid
        Shi, Slo = Shi + uhi, Shi
        k += 1
    end
    k < kmax || @warn "Iteration limit reached before convergence."
    return Shi
end

function pFqconformalrational(α::Tuple{T1, T1}, β::Tuple{T2}, z::T3, ::Val{3}; kmax::Int = KMAX) where {T1, T2, T3}
    T = promote_type(T1, T2, T3)
    γ = T(β[1])
    (α, β) = T.(α)
    if abs(z) < eps(real(T))
        return one(T)
    end
    w = r(T(z), Val(3))
    k = 0
    ulo = zero(T)
    umedlo = zero(T)
    umid = zero(T)
    umedhi = zero(T)
    uhi = one(T)
    Slo = one(T)
    uhi, umedhi, umid, umedlo, ulo = (w/γ)*6*α*β*uhi, uhi, umedhi, umid, umedlo
    Shi = Slo + uhi
    k = 1
    while (k < kmax && errcheck(Shi, Slo, 8eps(real(T))))
        alo = T(2*(k-5)*(k-4) + 4*(k-4))
        amedlo = T(2*(k-4)*(k-3) - 6*γ*(k-3) + 12*(α+β)*(k-3) + 4*(k-3))
        amid = T(4*(k-3)*(k-2) - 24*γ*(k-2) + 12*(α+β)*(k-2) + 12*(k-2) + 36*α*β)
        amedhi = T(4*(k-2)*(k-1) - 36*γ*(k-1) + 36*(α+β)*(k-1) + 12*(k-1) - 72*α*β)
        ahi = T(-6*(k-1)*k - 24*γ*k + 36*(α+β)*k + 36*α*β)
        uhi, umedhi, umid, umedlo, ulo = (w/(6*(k+1)*(k+γ)))*(ahi*uhi+w*(amedhi*umedhi+w*(amid*umid+w*(amedlo*umedlo+w*alo*ulo)))), uhi, umedhi, umid, umedlo
        Shi, Slo = Shi + uhi, Shi
        k += 1
    end
    k < kmax || @warn "Iteration limit reached before convergence."
    return Shi
end

function pFqconformalrational(α::Tuple{T1, T1}, β::Tuple{T2}, z::T3, ::Val{4}; kmax::Int = KMAX) where {T1, T2, T3}
    T = promote_type(T1, T2, T3)
    γ = T(β[1])
    (α, β) = T.(α)
    if abs(z) < eps(real(T))
        return one(T)
    end
    w = r(T(z), Val(4))
    k = 0
    ulo = zero(T)
    umedlo = zero(T)
    umid = zero(T)
    umedhi = zero(T)
    uhi = one(T)
    Slo = one(T)
    uhi, umedhi, umid, umedlo, ulo = (w/γ)*8*α*β*uhi, uhi, umedhi, umid, umedlo
    Shi = Slo + uhi
    k = 1
    while (k < kmax && errcheck(Shi, Slo, 8eps(real(T))))
        #alo = T(8*(k-5)*(k-4) - 8*γ*(k-4) + 16*(k-4))
        alo = (8*(k-γ) - 24)*(k-4)
        #amedlo = T(8*(k-4)*(k-3) - 40*γ*(k-3) + 64*(α+β)*(k-3) + 16*(k-3) - 64*α*β)
        amedlo = (8*k - 40*γ + 64*(α+β) - 16)*(k-3) - 64*α*β
        #amid = T(-80*γ*(k-2) + 64*(α+β)*(k-2) + 16*(k-2) + 192*α*β)
        amid = (-80*γ + 64*(α+β) + 16)*(k-2) + 192*α*β
        #amedhi = T(-80*γ*(k-1) + 64*(α+β)*(k-1) + 16*(k-1) - 192*α*β)
        amedhi = (-80*γ + 64*(α+β) + 16)*(k-1) - 192*α*β
        #ahi = T(-8*(k-1)*k - 40*γ*k + 64*(α+β)*k + 64*α*β)
        ahi = (-8*(k-1) - 40*γ + 64*(α+β))*k + 64*α*β
        uhi, umedhi, umid, umedlo, ulo = (w/(8*(k+1)*(k+γ)))*(ahi*uhi+w*(amedhi*umedhi+w*(amid*umid+w*(amedlo*umedlo+w*alo*ulo)))), uhi, umedhi, umid, umedlo
        Shi, Slo = Shi + uhi, Shi
        k += 1
    end
    k < kmax || @warn "Iteration limit reached before convergence."
    return Shi
end

function pFqconformalrational(α::Tuple{T1, Vararg{T1, q}}, β::NTuple{q, T2}, z::T3, ::Val{m}; kmax::Int = KMAX) where {q, m, T1, T2, T3}
    T = promote_type(T1, T2, T3)
    if abs(z) < eps(real(T))
        return one(T)
    end
    # Transform alpha and beta
    α = Tuple(2*m .* T.(α))
    β = Tuple(2*m .* (T.(β) .- 1))

    # Defining factors L_k and R_k as operators     # L_0⋯L_q u = R_1⋯R_{q+1} u
    L = zeros(T, m+2, q)
    for j in 1:q-1
        for k in 0:m+1
            L[k+1, j] = (binomial(m, k+1) * (1 - (-1)^(k+1)) - binomial(m, k-1) * (1 - (-1)^(k-1))) * (0 - k) + ((q - j) * (m - 1) + q - j - 1) * (binomial(m, k) * (1 - (-1)^k) + binomial(m, k-1) * (1 - (-1)^(k-1))) + (-1)^k * β[j] * binomial(m, k)
        end
    end
    for k in 0:m+1
        L[k+1, q] = 2*binomial(m, 2*floor(Int, k/2) + 1) * (0 - k) + (-1)^k * β[q] * binomial(m-1, k)
    end
    ΔL = zeros(T, m+2)
    for k in 0:2:m+1
        ΔL[k+1] = 2*(binomial(m, k+1) - binomial(m, k-1))
    end
    ΔLq = zeros(T, m+2)
    for k in 0:m+1
        ΔLq[k+1] = 2*binomial(m, 2*floor(Int, k/2) + 1)
    end

    R = zeros(T, m+2, q+1)
    for j in 1:q
        for k in 0:m+1
            R[k+1, j] = (binomial(m, k+1) * (1 - (-1)^(k+1)) - binomial(m, k-1) * (1 - (-1)^(k-1))) * (0 - k) + ((q - j + 1) * (m - 1) + q - j) * (binomial(m, k) * (1 - (-1)^k) + binomial(m, k-1) * (1 - (-1)^(k-1))) + (-1)^k * α[j] * binomial(m, k)
        end
    end
    for k in 0:m+1
        R[k+1, q+1] = 2*binomial(m, 2*floor(Int, k/2) + 1) * (0 - k) + (-1)^k * α[q+1] * binomial(m-1, k)
    end
    ΔR = zeros(T, m+2)
    for k in 0:2:m+1
        ΔR[k+1] = 2*(binomial(m, k+1) - binomial(m, k-1))
    end
    ΔRqp1 = zeros(T, m+2)
    for k in 0:m+1
        ΔRqp1[k+1] = 2*binomial(m, 2*floor(Int, k/2) + 1)
    end

    Lt0 = zeros(T, m+2)
    for k in 0:m+1
        Lt0[k+1] = (binomial(m, k+1) - binomial(m, k-1)) * (0 - k) + (q*m - 1) * binomial(m+1, k)
    end
    ΔLt0 = zeros(T, m+2)
    for k in 0:m+1
        ΔLt0[k+1] = binomial(m, k+1) - binomial(m, k-1)
    end

    # Initial values for the sum
    w   = r(T(z), Val(m))
    Slo = zero(T)
    Shi = one(T)

    ##############################################################################################################
    # Matrices of auxiliary sequences, Z and Zt
    Z  = zeros(T, m+2, q+1)
    Zt = zeros(T, m+2, q)
    # Initialize Z (n=0)
    Z[m+2, 1] = one(T)
    for i in 0:q-1
        Z[m+2, i+2] = sum(L[k+1, q-i]*Z[m+3-(k+1), i+1] for k in 0:m+1)/w
    end

    # Initialize Zt (n=0)
    Zt[m+2, 1] = sum(R[k+1, q+1]*Z[m+3-(k+1), 1] for k in 0:m+1)
    for i in 1:q-1
        Zt[m+2, i+1] = sum(R[k+1, q+1-i]*Zt[m+3-(k+1), i] for k in 0:m+1)/w
    end

    #############   Main loop updating matrices and summing

    n = 0
    @inbounds while (n < kmax && errcheck(Shi, Slo, 8eps(real(T))))
        # Increment L
        for j in 1:q-1
            for k in 0:2:m+1
                L[k+1, j] += ΔL[k+1]
            end
        end
        for k in 0:m+1
            L[k+1, q] += ΔLq[k+1]
        end

        # Increment Z
        vn_new = (sum(R[k+1, 1]*Zt[m+3-(k+1), q] for k in 0:m+1) - w*sum(Lt0[k+1]*Z[m+3-(k+1), q+1] for k in 0:m+1)) / (n+1)
        for j in 1:m+1
            Z[j, q+1] = w*Z[j+1, q+1]
        end
        Z[m+2, q+1] = vn_new

        for i in 1:q
            zni_new = w*(Z[m+2, q+1-i+1] - sum(L[k+1, i]*Z[m+4-(k+1), q+1-i] for k in 1:m+1)) / (2*(n+1)*m + β[i])
            for j in 1:m+1
                Z[j, q+1-i] = w*Z[j+1, q+1-i]
            end
            Z[m+2, q+1-i] = zni_new
        end

        # Increment Lt0
        for k in 0:m+1
            Lt0[k+1] += ΔLt0[k+1]
        end

        # Increment R
        for j in 1:q
            for k in 0:2:m+1
                R[k+1, j] += ΔR[k+1]
            end
        end
        for k in 0:m+1
            R[k+1, q+1] += ΔRqp1[k+1]
        end

        # Increment Zt
        zt1_new = sum(R[k+1, q+1]*Z[m+3-(k+1), 1] for k in 0:m+1)
        for j in 1:m+1
            Zt[j, 1] = w*Zt[j+1, 1]
        end
        Zt[m+2, 1] = zt1_new

        for i in 1:q-1
            ztni_new = sum(R[k+1, q+1-i]*Zt[m+3-(k+1), i] for k in 0:m+1)/w
            for j in 1:m+1
                Zt[j, i+1] = w*Zt[j+1, i+1]
            end
            Zt[m+2, i+1] = ztni_new
        end

        Shi, Slo = Shi + Z[m+2, 1], Shi
        n += 1
    end

    return Shi
end

"""
Maclaurin series of ₚ₊₁Fₚ(α, β; z) via the generalized Euler transformation with δ = ∑β-∑α:

ₚ₊₁Fₚ(α, β; z) = (1-z)^δ v(z = 1 - ((1-w)/(1+w))^m)

"""
pFqconformalrationalgeneralizedeuler(α, β, z; kwds...) = pFqconformalrationalgeneralizedeuler(α, β, z, Val(4); kwds...)
function pFqconformalrationalgeneralizedeuler(α::Tuple{Any, Vararg{Any, q}}, β::NTuple{q, Any}, z, m; kwds...) where q
    T1 = isempty(α) ? Any : mapreduce(typeof, promote_type, α)
    T2 = isempty(β) ? Any : mapreduce(typeof, promote_type, β)
    pFqconformalrationalgeneralizedeuler(T1.(α), T2.(β), z, m; kwds...)
end


function pFqconformalrationalgeneralizedeuler(α::Tuple{T1, Vararg{T1, q}}, β::NTuple{q, T2}, z::T3, ::Val{m}; kmax::Int = KMAX) where {q, m, T1, T2, T3}
    T = promote_type(T1, T2, T3)
    if abs(z) < eps(real(T))
        return one(T)
    end
    δ = sum(β)-sum(α)
    # Transform alpha and beta
    β = Tuple(2*m .* (T.(β) .- 1))
    α = Tuple(2*m .* T.(α))
    # Defining factors L_k and R_k as operators     # L_0⋯L_q u = R_1⋯R_{q+1} u

    L = zeros(T, m+2, q)
    for j in 1:q
        for k in 0:m+1
            L[k+1, j] = (binomial(m, k+1) * (1 - (-1)^(k+1)) - binomial(m, k-1) * (1 - (-1)^(k-1))) * (0 - k)- (1 - (-1)^k)* binomial(m, k) *((q - j) * m + β[j]+2m*(δ-q+j))+ (1 - (-1)^(k-1))*binomial(m, k-1)*(q-j)*m+ β[j] * binomial(m, k)
        end
    end
    ΔL = zeros(T, m+2)
    for k in 0:2:m+1
        ΔL[k+1] = 2*(binomial(m, k+1) - binomial(m, k-1))
    end

    R = zeros(T, m+2, q+1)
    for j in 1:q+1
        for k in 0:m+1
            R[k+1, j] = (binomial(m, k+1) * (1 - (-1)^(k+1)) - binomial(m, k-1) * (1 - (-1)^(k-1))) * (0 - k) - (1 - (-1)^k)*binomial(m, k) *((q - j + 1) *m + α[j] + 2m*(δ-q+j-1)) +  (1 - (-1)^(k-1))*binomial(m, k-1)*(q-j+1)*m+ α[j] * binomial(m, k)
        end
    end
    ΔR = zeros(T, m+2)
    for k in 0:2:m+1
        ΔR[k+1] = 2*(binomial(m, k+1) - binomial(m, k-1))
    end

    # Truncated L_0
    Lt0 = zeros(T, m+3)
    for k in 0:m+2
        Lt0[k+1] = (binomial(m, k+1) - binomial(m, k-1)) * (0 - k) + q*m*binomial(m, k-1)+m*(q-2δ)*binomial(m,k)
    end
    ΔLt0 = zeros(T, m+2)
    for k in 0:m+1
        ΔLt0[k+1] = binomial(m, k+1) - binomial(m, k-1)
    end

    # Initial values for the sum
    w   = r(T(z), Val(m))
    Slo = zero(T)
    Shi = one(T)

    ##############################################################################################################
    # Matrices of auxiliary sequences, Z and Zt
    Z  = zeros(T, m+2, q+1)
    Zt = zeros(T, m+2, q)
    # Initialize Z (n=0)
    Z[m+2, 1] = one(T)
    for i in 0:q-1
        Z[m+2, i+2] = sum(L[k+1, q-i]*Z[m+3-(k+1), i+1] for k in 0:m+1)/w
    end

    # Initialize Zt (n=0)
    Zt[m+2, 1] = sum(R[k+1, q+1]*Z[m+3-(k+1), 1] for k in 0:m+1)
    for i in 1:q-1
        Zt[m+2, i+1] = sum(R[k+1, q+1-i]*Zt[m+3-(k+1), i] for k in 0:m+1)/w
    end

    #############   Main loop updating matrices and summing

    n = 0
    @inbounds while (n < kmax && errcheck(Shi, Slo, 8eps(real(T))))
        # Increment L
        for j in 1:q
            for k in 0:2:m+1
                L[k+1, j] += ΔL[k+1]
            end
        end

        # Increment Z
        vn_new = (sum(R[k+1, 1]*Zt[m+3-(k+1), q] for k in 0:m+1) - w*sum(Lt0[k+1]*Z[m+3-(k+1), q+1] for k in 0:m+1)) / (n+1)
        for j in 1:m+1
            Z[j, q+1] = w*Z[j+1, q+1]
        end
        Z[m+2, q+1] = vn_new

        for i in 1:q
            zni_new = w*(Z[m+2, q+1-i+1] - sum(L[k+1, i]*Z[m+4-(k+1), q+1-i] for k in 1:m+1)) / (2*(n+1)*m + β[i])
            for j in 1:m+1
                Z[j, q+1-i] = w*Z[j+1, q+1-i]
            end
            Z[m+2, q+1-i] = zni_new
        end

        # Increment Lt0
        for k in 0:m+1
            Lt0[k+1] += ΔLt0[k+1]
        end

        # Increment R
        for j in 1:q+1
            for k in 0:2:m+1
                R[k+1, j] += ΔR[k+1]
            end
        end

        # Increment Zt
        zt1_new = sum(R[k+1, q+1]*Z[m+3-(k+1), 1] for k in 0:m+1)
        for j in 1:m+1
            Zt[j, 1] = w*Zt[j+1, 1]
        end
        Zt[m+2, 1] = zt1_new

        for i in 1:q-1
            ztni_new = sum(R[k+1, q+1-i]*Zt[m+3-(k+1), i] for k in 0:m+1)/w
            for j in 1:m+1
                Zt[j, i+1] = w*Zt[j+1, i+1]
            end
            Zt[m+2, i+1] = ztni_new
        end

        Shi, Slo = Shi + Z[m+2, 1], Shi
        n += 1
    end

    return Shi*exp(δ*log1p(-z))
end
