@testset "Special function" begin
    @test pochhammer(2,3) == 24
    @test pochhammer(0.5,3) == 0.5*1.5*2.5
    @test pochhammer(0.5,0.5) == 1/sqrt(pi)
    @test pochhammer(0,1) == 0
    @test pochhammer(-1,2) == 0
    @test pochhammer(-5,3) == -60
    @test pochhammer(-1,-0.5) == 0
    @test 1.0/pochhammer(-0.5,-0.5) == 0
    @test pochhammer(-1+0im,-1) == -0.5
    @test pochhammer(2,1) == pochhammer(2,1.0) == pochhammer(2.0,1) == 2
    @test pochhammer(1.1,2.2) ≈ gamma(3.3)/gamma(1.1)
    @test pochhammer(-2,1) == pochhammer(-2,1.0) == pochhammer(-2.0,1) == -2
    @test pochhammer(3, 1:5) == [3, 12, 60, 360, 2520]
end

@testset "Helper functions" begin
    @test HypergeometricFunctions.log1pover(-0.25) ≈ log1p(-0.25) / -0.25
    @test HypergeometricFunctions.log1pover(0.0) == 1.0

    @test HypergeometricFunctions.speciallog(0.0) ≈ 1.0
    @test HypergeometricFunctions.speciallog(0.1) ≈ HypergeometricFunctions.speciallog(big(0.1)) rtol = 1e-14
    @test HypergeometricFunctions.speciallog(-0.1) ≈ HypergeometricFunctions.speciallog(big(-0.1)) rtol = 1e-14

    @test HypergeometricFunctions.logandpoly(0.1) ≈ HypergeometricFunctions.logandpoly(big(0.1)) rtol = 1e-14
    @test HypergeometricFunctions.logandpoly(-0.1) ≈ HypergeometricFunctions.logandpoly(big(-0.1)) rtol = 1e-14

    @test HypergeometricFunctions.sqrtasinsqrt(0.0) == 1.0
    @test HypergeometricFunctions.sqrtasinsqrt(0.2) ≈ asin(sqrt(0.2)) / sqrt(0.2)
    @test HypergeometricFunctions.sqrtatanhsqrt(0.2) ≈ atanh(sqrt(0.2)) / sqrt(0.2)
    @test HypergeometricFunctions.sinnasinsqrt(3, 0.2) ≈ sin(3 * asin(sqrt(0.2))) / (3 * sqrt(0.2))
    @test HypergeometricFunctions.cosnasinsqrt(3, 0.2) ≈ cos(3 * asin(sqrt(0.2)))
end
