@testset "M" begin
    @test M(-3, -3, 0.5) ≡ exp(0.5)
    @test M(0, -1, 10) ≡ 1.0
    @test M(1.2,  0.0, 0.1) ==  Inf # Mimick gamma( 0.0) =  Inf
    @test M(1.2, -0.0, 0.1) == -Inf # and    gamma(-0.0) = -Inf
    @test_throws DomainError M(1, -2, 0.5)
    @test_throws DomainError M(-3, -2, 0.5)
    @test M(-2, -3, 0.5) ≡ 1.375
    @test M(0.5, 1.5, -1000) ≈ 0.028024956081989644 # From #46
    @test M(1, 2, 0) == 1
    @test M(1, 2, 0.25) == expm1(0.25)/0.25
    for (S, T) in ((Float64, BigFloat),)
        a = T(8.9)
        b = T(0.5)
        for x in T(-36):T(2):T(70)
            @test M(S(a), S(b), S(x)) ≈ S(M(a, b, x)) # From #45
        end
        a = T(5)/6
        b = T(1)/2
        for x in T(-5):T(0.25):T(5)
            @test M(S(a), S(b), -S(x)^2) ≈ S(M(a, b, -x^2)) # From #66
        end
        b = 1
        x = T(1)/3
        for a in S(1):S(0.5):S(7)
            @test M(a, b, S(x)) ≈ S(M(a, b, x))
        end
    end
end

@testset "Confluent near-singular parameters" begin
    for δ in (1e-8, 1e-10)
        @test M(1.2, 1 + δ, 0.3) ≈ Float64(M(big(1.2), big(1 + δ), big(0.3)))
        @test M(-2.5, -1 + δ, 0.2) ≈ Float64(M(big(-2.5), big(-1 + δ), big(0.2)))
        @test _₂F₁(0.8, 0.7, 1.5 + δ, 0.9) ≈ Float64(_₂F₁(big(0.8), big(0.7), big(1.5 + δ), big(0.9)))
    end
end

@testset "U" begin
    @test U(1, 1, 1.f0) ≈ 0.5963473623231942 # the Euler series
    @test U(1, 1, 1) == 0.5963473623231942
    for (S, T) in ((Float64, BigFloat),)
        b = 0
        x = T(1)/3
        for a in S(1):S(0.5):S(7)
            @test U(a, b, S(x)) ≈ S(U(a, b, x)) # From #55
        end
    end
end
