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
