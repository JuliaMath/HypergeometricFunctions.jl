@testset "Exact identities" begin
    z = 0.3
    @test M(2.5, 2.5, z) ≈ exp(z)
    @test pFq((3,), (), z) ≈ (1 - z)^(-3)
    @test pFq([1.0, 2.0], [4.0], 0.25) == pFq((1.0, 2.0), (4.0,), 0.25)
    @test _₂F₁(0, 2, 3, 0.7) == 1.0
    @test _₂F₁(2, 3, 2, 0.25) == pFq((3,), (), 0.25)
    @test _₂F₁(2, 3, 3, 0.25) == pFq((2,), (), 0.25)
end

@testset "Singular parameters" begin
    @test_throws DomainError M(1.2, -1, 0.1)
    @test isinf(_₂F₁(1, 2, 2, 1.0))
    @test isnan(_₂F₁(1, 2 + im, 3, 1))
end

@testset "Polyalgorithm boundaries" begin
    α, β = (1.0, 2.0), (3.0,)
    @test pFq(α, β, 0.72) ≈ pFq(Tuple(α), Tuple(β), 0.72)
    @test pFq(α, β, 0.73) ≈ pFq(Tuple(α), Tuple(β), 0.73)
    @test pFq((1.0,), (2.0,), 0.0) == 1.0
    @test pFq((1.0, 2.0), (3.0,), 0.0) == 1.0
end

@testset "Branch cut consistency" begin
    for ϵ in (1e-8, 1e-10)
        z₊ = -0.3 + im * ϵ
        z₋ = -0.3 - im * ϵ
        @test pFq((0.5,), (), z₊) == conj(pFq((0.5,), (), z₋))

        z₊ = 1.2 + im * ϵ
        z₋ = 1.2 - im * ϵ
        @test M(0.75, 1.25, z₊) == conj(M(0.75, 1.25, z₋))

        z₊ = 0.9 + im * ϵ
        z₋ = 0.9 - im * ϵ
        @test _₂F₁(0.8, 0.7, 1.5, z₊) == conj(_₂F₁(0.8, 0.7, 1.5, z₋))
    end
end

@testset "Warning symbols" begin
    @test pFq2string(Val(86420), Val(97531)) == "₈₆₄₂₀F₉₇₅₃₁"
end
