using HypergeometricFunctions, SpecialFunctions, Test
import LinearAlgebra: norm
import HypergeometricFunctions: iswellpoised, isalmostwellpoised, M, U,
                                pochhammer, _₂F₁general, _₂F₁general2,
                                pFqdrummond, pFqweniger, pFq2string

const rtol = 1.0e-3
const NumberType = Float64

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

@testset "Hypergeometric Functions" begin
    @testset "_₂F₁ vs _₂F₁general vs _₂F₁general2" begin
        e = exp(1.0)
        regression_max_accumulated_error = 2^12*eps() # discovered by running the test
        for z in (.9rand(Float64, 10), 10rand(ComplexF64, 10))
            j = 1
            for (a,b,c) in ((√2/e, 1.3, 1.3), (1.2, √3, 1.2), (-0.4, 0.4, 0.5),
                            (-0.3, 1.3, 0.5), (0.35, 0.85, 0.5), (0.5, 0.5, 1.5),
                            (1.0, 1.0, 1.5), (0.5, 1.0, 1.5), (0.3, 0.7, 1.5),
                            (0.7, 1.3, 1.5), (0.35, 0.85, 1.5), (1.0, 1.0, 2.0),
                            (3.0, 1.0, 2.0), (-2.0, 1.0, 2.0), (-3.0, 1.0, 2.0),
                            (1.0, -4.0, 2.0), (2.0, 2.0, 4.0), (1.0, 1.5, 2.5))
                error_accum = 0.0
                for zi in z
                    twoFone = _₂F₁(a,b,c,zi)
                    aa,bb,cc = big(a),big(b),big(c)
                    twoFonegeneral = convert(Complex{Float64},_₂F₁general(aa, bb, cc, big(zi)))
                    norm(twoFone / twoFonegeneral - 1) > sqrt(eps()) && println("This is ₂F₁($a,$b;$c;zi) - ₂F₁general($a,$b;$c;zi): ", norm(twoFone / twoFonegeneral - 1), "   ", twoFone, "   ", twoFonegeneral, "   ", isfinite(twoFone), "   ", isfinite(twoFonegeneral), " this is zi: ", zi)
                    error_accum += Float64(norm(twoFone / twoFonegeneral - 1))
                end
                @test error_accum < regression_max_accumulated_error
                j += 1
            end
            for (a,b,c) in ((√2/e, 1.3, 1.3), (1.2, √3, 1.2), (-0.4, 0.4, 0.5),
                            (-0.3, 1.3, 0.5), (0.35, 0.85, 0.5), (0.5, 1.0, 1.5),
                            (0.3, 0.7, 1.5), (0.7, 1.3, 1.5), (0.35, 0.85, 1.5),
                            (3.0, 1.0, 2.0), (-2.0, 1.0, 2.0), (-3.0, 1.0, 2.0),
                            (1.0, -4.0, 2.0), (1.0, 1.5, 2.5))
                error_accum = 0.0
                for zi in z
                    twoFone = _₂F₁(a,b,c,zi)
                    aa,bb,cc = big(a),big(b),big(c)
                    twoFonegeneral2 = convert(Complex{Float64},_₂F₁general2(aa, bb, cc, big(zi)))
                    norm(twoFone / twoFonegeneral2 - 1) > sqrt(eps()) && println("This is ₂F₁($a,$b;$c;zi) - ₂F₁general2($a,$b;$c;zi): ", norm(twoFone / twoFonegeneral2 - 1), "   ", twoFone, "   ", twoFonegeneral2, "   ", isfinite(twoFone), "   ", isfinite(twoFonegeneral2), " this is zi: ", zi)
                    error_accum += Float64(norm(twoFone / twoFonegeneral2 - 1))
                end
                @test error_accum < regression_max_accumulated_error
                j += 1
            end
        end
        @test _₂F₁general(-0.4, 0.4, 0.5, 0.75+0.75im) ≈ _₂F₁general2(-0.4, 0.4, 0.5, 0.75+0.75im)
    end

    @testset "Test that _₂F₁ is inferred for Float32 arguments" begin
        @test @inferred(_₂F₁(0.3f0, 0.7f0, 1.3f0, 0.1f0)) ≈ Float32(_₂F₁(0.3, 0.7, 1.3, 0.1))
    end

    @testset "method = positive" begin
        for (a, b, c, z) in ((1, 2, 3, 0.5), (3, 5, 7, 0.75), (1, 8537, 6042, 0.25))
            positivetwoFone = _₂F₁(a, b, c, z; method = :positive)
            twoFone = Float64(pFqweniger((BigFloat(a), BigFloat(b)), (BigFloat(c), ), big(z)))
            @test positivetwoFone > 0
            @test positivetwoFone ≈ twoFone
        end
    end

    @testset "_₂F₁ argument unity" begin
        for (a, b, c, z) in ((1, 2, 3, 1), (1, 2, 4, 1), (3, 5, 7, 1.0), (1, 8537, 6042, 1.0))
            twoFone = _₂F₁(a, b, c, z)
            twoFonep = _₂F₁(a, b, c, z; method = :positive)
            if iswellpoised(a, b, c)
                @test isfinite(twoFone)
                @test isfinite(twoFonep)
            else
                @test isinf(twoFone)
                @test isinf(twoFonep)
            end
        end
        for (a, b, c, z) in ((1, 2+im, 3, 1), (1, 2+im, 3+im, 1), (1, 2 + 1im, 2, 1))
            if isalmostwellpoised(a, b, c)
                if a+b==c
                    @test isinf(_₂F₁(a, b, c, z))
                else
                    @test isnan(_₂F₁(a, b, c, z))
                end
            else
                @test isinf(_₂F₁(a, b, c, z))
            end
        end
    end

    @testset "pFq vs mpmath" begin
        (a, b, c, result) = NumberType.([-1.138]), NumberType.([-1.17865]), NumberType(0.29524203533943405), 1.3211939223293958
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-0.00209496]), NumberType.([1.55444]), NumberType(1.9964304489951332), 0.9956895798837077
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-0.550861]), NumberType.([-2.76028]), NumberType(1.686164802648714), 2.26529646755787
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([0.522481]), NumberType.([2.43767]), NumberType(0.07397564666157663), 1.0161190845825492
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([1.14629]), NumberType.([-1.71237]), NumberType(0.865999806211927), 6.09627259912178
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-0.977563]), NumberType.([-0.943489]), NumberType(1.3295519796473885), 2.9650972591473237
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([2.19487]), NumberType.([-2.84673]), NumberType(-2.946765152780697), -10.159918912696533
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-1.40648]), NumberType.([-2.09642]), NumberType(0.003973775403141477), 1.002667944757093
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([2.78115]), NumberType.([-0.974578]), NumberType(-0.289876068484225), -9.436220915244991
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-2.0195]), NumberType.([2.97415]), NumberType(1.9088632413219395), 0.020342030150391748
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-1.29058]), NumberType.([1.12475, 2.20863]), NumberType(2.5233777332943026), -0.23709345988925976
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([1.93422]), NumberType.([-0.403018, -1.68818]), NumberType(1.8841553851312423), -368.4667336251361
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([1.38965]), NumberType.([1.34519, 2.18274]), NumberType(-2.967668498991821), 0.12382284309201524
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([0.317445]), NumberType.([2.83319, 2.08577]), NumberType(-0.8292545561139764), 0.9574454648148745
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-1.53034]), NumberType.([-2.00218, 2.71893]), NumberType(-0.055950986916451395), 0.984386941287138
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([1.8068]), NumberType.([2.785, -2.27362]), NumberType(-2.965224079697172), 3.343788228467509
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([2.01277]), NumberType.([0.0297853, 0.638178]), NumberType(2.87029742430684), 1862.664910896795
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-1.332]), NumberType.([-0.571085, 2.40734]), NumberType(2.5191976052414615), 2.6748631921200974
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([0.585852]), NumberType.([-0.840541, -0.914398]), NumberType(-0.7114454462840016), 13.397306612348904
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-1.30148]), NumberType.([-0.0397773, 0.452281]), NumberType(2.9085746526416285), 134.70529132344416
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([1/2]), NumberType.([1, 3/2]), NumberType(-10.0), 0.1215046687753564
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([1/2]), NumberType.([1, 3/2]), NumberType(-100.0), 0.05291894107105639
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([1/2]), NumberType.([1, 3/2]), NumberType(-1000.0), 0.015220788412455796
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([1/2]), NumberType.([1, 3/2]), NumberType(-10000.0), 0.004728870002692929
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([1.79664]), NumberType.([2.5454, -2.0446, -1.84776]), NumberType(-1.5647341841975129), 10.224086634371103
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-1.08954]), NumberType.([-1.09992, -0.0898055, -1.21452]), NumberType(2.911146030454528), -281.8950930467646
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([2.50123]), NumberType.([1.481, 1.04535, -1.1706]), NumberType(1.2093267908720815), 4.297922628682247
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-0.918156]), NumberType.([-1.66879, -1.26548, -2.71329]), NumberType(-1.4539176485483885), 0.7245506668664914
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-2.25871]), NumberType.([0.734527, 2.69253, 1.52627]), NumberType(1.412857816927049), 0.0007367729643306418
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([1.75347]), NumberType.([-1.19881, 1.46644, 0.377671]), NumberType(-2.894612064681616), 31.986877557865302
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-2.55356]), NumberType.([2.31273, -1.79583, 2.71976]), NumberType(-1.8211352681511421), 0.6524733226915893
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-2.96471]), NumberType.([0.10731, 0.150247, -0.749495]), NumberType(-0.5456496601437784), -364.7310910695041
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([2.28595]), NumberType.([1.24272, -1.3733, -1.97931]), NumberType(2.694315046446391), 1427.85053117312
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-1.02891]), NumberType.([-0.203779, 0.414445, -2.72581]), NumberType(2.7264202122777323), -11.325191332268156
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([1.21643, 2.56177]), NumberType.([-2.85355]), NumberType(-1.173143685351698), -2.6927101990763944
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-0.884016, -0.278806]), NumberType.([-2.95573]), NumberType(-0.5131107181888153), 1.0423537583061855
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-0.489218, 1.02043]), NumberType.([-1.8805]), NumberType(-2.054787166536263), 1.3985297910285344
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-2.61829, -0.0405846]), NumberType.([-2.52321]), NumberType(-2.8547909742932127), 1.111640235796153
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-1.78806, 0.244808]), NumberType.([-2.03175]), NumberType(0.5005127629084343), 1.0071139633716848
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([2.79663, 2.95841]), NumberType.([-2.23659]), NumberType(-2.1658329294772676), -1.298805755547899
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([0.473781, -0.844289]), NumberType.([-0.228894]), NumberType(0.8762669945656332), 3.026596343521232
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-1.99653, -2.15667]), NumberType.([1.45492]), NumberType(0.5669144675224369), 2.901074111158336
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([1.89145, -2.17305]), NumberType.([-1.97603]), NumberType(-2.2566276457677343), 145.35756041830106
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([0.998819, -0.225309]), NumberType.([0.447173]), NumberType(-0.10835817353459554), 1.051599100560819
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-0.383828, 1.77538]), NumberType.([-1.0875, -2.84315]), NumberType(1.4065086054194986), 715.336897395762
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([0.733652, -1.6307]), NumberType.([-0.0608787, 0.896927]), NumberType(0.10210464529041374), 3.1665794151940694
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-0.202523, -1.68826]), NumberType.([-1.92541, 2.88428]), NumberType(0.23168614635513318), 0.9854522687722247
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-0.946942, 1.58094]), NumberType.([-0.150083, 2.30757]), NumberType(-2.109785441324441), -7.76796376651872
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([1.55118, -2.75387]), NumberType.([1.46214, 2.18634]), NumberType(2.6348816640393258), -0.31994633361734
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-2.06745, -0.28521]), NumberType.([-1.99816, -0.766774]), NumberType(-0.7965581191453945), 5.472490531872681
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-2.20859, -1.41744]), NumberType.([1.85829, 1.92196]), NumberType(-1.2487286295352975), -0.053122867013287486
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([0.260072, 0.317936]), NumberType.([-1.4299, 1.60909]), NumberType(-1.554201646172852), 1.0698093358492726
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-1.81949, 2.65844]), NumberType.([2.84994, -1.67876]), NumberType(1.5573520464968), 4.557008143231845
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-2.09665, 1.27147]), NumberType.([-0.682723, 0.192544]), NumberType(-1.175497640118377), -118.5195721586634
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([0.472347, 2.38973]), NumberType.([1.23953, -2.45145, -0.0816949]), NumberType(2.9793981957869673), 562.1327877798817
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([2.70886, -2.26978]), NumberType.([-1.83008, -1.6277, -1.33282]), NumberType(2.7790611310114053), -12231.7908950255
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([1.48786, -2.73106]), NumberType.([-1.73601, -1.676, -2.17108]), NumberType(-0.6847891664895385), -36.681955083411964
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-1.72604, -1.83087]), NumberType.([-2.68909, -0.904331, 1.93759]), NumberType(1.651593934432523), 0.9572670911067562
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([0.27953, 2.57139]), NumberType.([-1.72806, 2.13894, -0.674948]), NumberType(1.3654741919794335), -8.061422671050726
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-0.363687, 1.0458]), NumberType.([-1.13771, -2.02746, -2.86789]), NumberType(2.026179003725624), -5919.391531381284
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([2.49209, -1.06563]), NumberType.([1.47807, 1.01142, -1.54966]), NumberType(-0.9883650134564901), -0.0985713724396932
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([0.416307, 1.24566]), NumberType.([1.56196, 1.32576, -1.82953]), NumberType(-2.4172898932733053), 1.0404646069316064
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([2.52485, 0.109503]), NumberType.([-1.3165, 0.177332, 2.58557]), NumberType(-2.84982197723157), 2.9369917768818667
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-2.45525, -1.61598]), NumberType.([-1.75953, -1.59772, 2.05578]), NumberType(-1.1643127384783538), 0.5520841179227342
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-1.7029, 1.4496, -2.82896]), NumberType.([-2.77533]), NumberType(-1.5171018909198049), 8.437330705993508
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-1.34783, -0.324078, 2.44238]), NumberType.([-0.777003]), NumberType(-1.059855902026896), 3.9754790978546266
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-1.87289, 2.99671, 2.99742]), NumberType.([-1.4312]), NumberType(-2.433424057786612), 638.863801591021
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-0.963631, -1.36196, -2.90998]), NumberType.([-0.950805]), NumberType(-1.1275151312466227), -1.9506846998048282
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([2.20585, 1.21729, 0.419755]), NumberType.([-2.95176]), NumberType(-0.8474785819006732), -0.9897576725972367
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-2.50248, -2.01281, -1.10468]), NumberType.([-1.47035]), NumberType(-2.161852515937607), -4.205644305802463
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-1.373, 0.0281506, 1.31409]), NumberType.([-0.775952]), NumberType(-2.0586457735217407), 0.6176889810366223
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([1.4476, 0.286052, 1.71704]), NumberType.([1.29826]), NumberType(-0.9902327227931127), 0.7640485766337997
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-2.49284, 1.34398, 1.38886]), NumberType.([-2.04354]), NumberType(-0.9516109943455211), -137.34441283528875
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([2.19033, -0.446436, -2.73578]), NumberType.([0.758367]), NumberType(-1.115415216535085), -8.890490669826645
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([1.21084, -1.43893, 1.56842]), NumberType.([-0.598469, 0.52221]), NumberType(-2.2391561626302376), -64.66480492684083
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([0.607277, -0.376803, -2.08424]), NumberType.([1.92106, 2.80899]), NumberType(-1.5589765316931263), 0.8516398920401931
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-2.63828, 2.59771, 2.54461]), NumberType.([2.62021, -0.0157026]), NumberType(0.5526474385388149), -29.571181883739737
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([0.314356, 2.59205, -1.28807]), NumberType.([-1.22202, -1.66475]), NumberType(-0.7121680260550485), 0.856714022601912
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([0.694337, -2.06707, 0.214575]), NumberType.([-1.31576, 0.266519]), NumberType(-2.9461185849069706), 21.80248076992691
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([1.90394, 2.56441, 0.07523]), NumberType.([0.0390543, -1.65175]), NumberType(-0.7780439739130411), -3.184516491992518
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-2.80752, 2.7876, -1.37287]), NumberType.([-2.9165, 0.49188]), NumberType(0.39686964583141693), -1.3014023309227971
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([2.36231, -0.901808, 1.00893]), NumberType.([1.39985, -2.34968]), NumberType(0.970806084927347), 42731.87930657499
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        @test HypergeometricFunctions.pFqcontinuedfraction(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-0.715302, 2.94101, -0.498188]), NumberType.([-0.393808, -0.0265878]), NumberType(0.4375358657735875), 61.68955383799459
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-1.14105, 1.9464, -2.44787]), NumberType.([0.54636, -0.654353]), NumberType(-2.9505198429994763), -53.189709409599516
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([0.228056, -1.24092, 0.843438]), NumberType.([1.25348, 2.42447, -1.59626]), NumberType(-2.32373843841519), 0.8942755488371047
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-1.78019, -2.49326, -1.62843]), NumberType.([-1.85186, -2.72869, -2.36236]), NumberType(-1.6400906330598044), 0.2981496993840609
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([1.98727, 2.31523, -0.0735596]), NumberType.([-0.0644426, 2.54497, -1.08036]), NumberType(-1.0705461999772998), 10.575143699481718
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-1.45036, -2.25443, -2.74298]), NumberType.([-1.39353, -1.20841, 0.225826]), NumberType(-2.950872722868197), 982.2060223576477
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-0.36497, 1.79508, 2.36814]), NumberType.([1.55155, 0.704328, -0.686943]), NumberType(0.2186808938362801), 1.7064306875530209
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-1.51876, -1.37851, -1.12183]), NumberType.([1.3349, -2.95129, 1.10139]), NumberType(0.7674847040296018), 1.4158192973160308
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-1.99362, -2.09819, -2.16471]), NumberType.([-1.76093, 2.0683, -0.947339]), NumberType(-1.8729057277220629), -41.660320326270636
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([1.82466, 1.28428, 2.67447]), NumberType.([2.10754, 0.871646, 1.39154]), NumberType(2.658591557834844), 81.52661900335721
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-1.79631, -0.0269564, 0.915417]), NumberType.([-2.17853, 1.0933, -0.000243005]), NumberType(2.661306310651646), -83.69268822607185
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([1.36163, -0.758606, -0.716056]), NumberType.([-1.04856, 0.266574, -2.47814]), NumberType(1.8461964795702661), -113.75411958323524
        @test pFq(a, b, c) ≈ result atol=eps() rtol=rtol
    end
    @testset "_₂F₁ vs mpmath" begin
        (a, b, c, result) = NumberType.([-0.509718, 2.74761]), NumberType.([-0.650532]), NumberType(-2.9555898601070316), 1.2790428246176058
        @test _₂F₁(a..., b..., c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-2.5215, -1.50451]), NumberType.([-0.759254]), NumberType(-0.5749286896537393), 1.1388556009105488
        @test _₂F₁(a..., b..., c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-2.86391, 1.25904]), NumberType.([0.517432]), NumberType(-2.940843956490598), 189.9291251487297
        @test _₂F₁(a..., b..., c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([0.0877358, 0.682746]), NumberType.([2.84558]), NumberType(-1.0261339665475955), 0.9822590779186833
        @test _₂F₁(a..., b..., c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([2.3191, 1.94423]), NumberType.([1.88828]), NumberType(-1.3591292051020116), 0.12818348728367582
        @test _₂F₁(a..., b..., c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([2.56295, -0.753746]), NumberType.([-0.159836]), NumberType(-0.8246738184622906), -6.594728994176734
        @test _₂F₁(a..., b..., c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([1.67125, 0.101391]), NumberType.([0.107343]), NumberType(-0.7965765135904737), 0.40781958014265957
        @test _₂F₁(a..., b..., c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-1.99678, -0.201406]), NumberType.([-1.95704]), NumberType(0.05697188547835097), 0.9880143477892287
        @test _₂F₁(a..., b..., c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-1.45054, 0.439005]), NumberType.([0.908234]), NumberType(0.8723075146005943), 0.49786449277738715
        @test _₂F₁(a..., b..., c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([1.01798, -2.7207]), NumberType.([-2.9141]), NumberType(-2.164369297479041), 8.70430095368102
        @test _₂F₁(a..., b..., c) ≈ result atol=eps() rtol=rtol
    end
    @testset "_₃F₂ vs mpmath" begin
        (a, b, c, result) = NumberType.([-2.76481, 0.00672772, 0.741886]), NumberType.([-2.53237, 2.79933]), NumberType(-0.21967338526219615), 0.9995944744632735
        @test _₃F₂(a..., b..., c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([1.00388, -2.50276, -1.20821]), NumberType.([-1.552, 1.0715]), NumberType(0.6716163785419433), -0.050120516406253396
        @test _₃F₂(a..., b..., c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([0.483392, -1.42038, -1.82681]), NumberType.([2.98456, 1.87387]), NumberType(-1.241922052077597), 0.7292186256180657
        @test _₃F₂(a..., b..., c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-2.80985, 1.57763, 0.344136]), NumberType.([-2.00075, -0.973491]), NumberType(-0.3899943196975246), 14466.732269418837
        @test _₃F₂(a..., b..., c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-2.3729, 1.09037, -2.76876]), NumberType.([-1.13014, -1.81025]), NumberType(-2.3731909909995714), -1351.450695124059
        @test _₃F₂(a..., b..., c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([2.60036, 1.54435, -2.72411]), NumberType.([-0.141692, -2.06915]), NumberType(0.813329638552919), -79135.94437145953
        @test _₃F₂(a..., b..., c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-1.32094, 1.1992, -0.270438]), NumberType.([-0.483108, -2.41227]), NumberType(-2.1529259273055366), 0.6103979954183846
        @test _₃F₂(a..., b..., c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-1.39151, 0.699785, -0.457295]), NumberType.([-2.53724, -0.45027]), NumberType(0.4811844516239163), 1.1524122102606527
        @test _₃F₂(a..., b..., c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-2.22368, 2.92094, -1.8566]), NumberType.([-2.35922, 2.78358]), NumberType(-0.06437623035744533), 1.121248505469995
        @test _₃F₂(a..., b..., c) ≈ result atol=eps() rtol=rtol
        (a, b, c, result) = NumberType.([-2.95746, -1.43804, -0.363603]), NumberType.([-1.59028, 1.87905]), NumberType(0.6408597427285083), 1.3015755922356282
        @test _₃F₂(a..., b..., c) ≈ result atol=eps() rtol=rtol
        @test _₃F₂(2.5, 3.25, 0.125) == _₃F₂(2.5, 1, 1, 3.25, 2, 0.125)
    end

end

@testset "₀F₀" begin
    for (S, T) in ((Float16, Float32),
                   (Float32, Float64),
                   (Float64, BigFloat),
                   (BigFloat, BigFloat))
        CS = Complex{S}
        CT = Complex{T}
        atol = rtol = 1000eps(S)
        for z in S(-4):S(0.25):S(1)
            @test pFqdrummond((), (), z) ≈ S(pFq((), (), T(z))) atol=atol rtol=rtol
            @test pFqweniger((), (), z) ≈ S(pFq((), (), T(z))) atol=atol rtol=rtol
        end
        atol *= 2
        rtol *= 2
        for x in S(-4):S(0.25):S(1), y in S(-4):S(0.25):S(1)
            z = complex(x, y)
            @test pFqdrummond((), (), z) ≈ CS(pFq((), (), CT(z))) atol=atol rtol=rtol
            @test pFqweniger((), (), z) ≈ CS(pFq((), (), CT(z))) atol=atol rtol=rtol
        end
    end
end

@testset "₁F₀" begin
    for (S, T) in ((Float32, Float64),
                   (Float64, BigFloat),
                   (BigFloat, BigFloat))
        CS = Complex{S}
        CT = Complex{T}
        atol = rtol = 1000eps(S)
        for α in S(-1.5):S(0.5):S(1.5), z in S(-0.75):S(0.25):S(0.25)
            @test pFqdrummond((α, ), (), z) ≈ S(pFq((T(α), ), (), T(z))) atol=atol rtol=rtol
            @test pFqweniger((α, ), (), z) ≈ S(pFq((T(α), ), (), T(z))) atol=atol rtol=rtol
        end
        atol *= 2
        rtol *= 2
        for α in S(-1.5):S(0.5):S(1.5), x in S(-0.75):S(0.25):S(0.25), y in S(-0.75):S(0.25):S(0.25)
            z = complex(x, y)
            @test pFqdrummond((α, ), (), z) ≈ CS(pFq((T(α), ), (), CT(z))) atol=atol rtol=rtol
            @test pFqweniger((α, ), (), z) ≈ CS(pFq((T(α), ), (), CT(z))) atol=atol rtol=rtol
        end
    end
end

@testset "₀F₁" begin
    for (S, T) in ((Float32, Float64),
                   (Float64, BigFloat),
                   (BigFloat, BigFloat))
        CS = Complex{S}
        CT = Complex{T}
        atol = rtol = 1000eps(S)
        for α in S(-1.5):S(1.0):S(1.5), z in S(-0.75):S(0.25):S(0.75)
            @test pFqdrummond((), (α, ), z) ≈ S(pFq((), (T(α), ), T(z))) atol=atol rtol=rtol
            @test pFqweniger((), (α, ), z) ≈ S(pFq((), (T(α), ), T(z))) atol=atol rtol=rtol
        end
        atol *= 2
        rtol *= 2
        for α in S(-1.5):S(1.0):S(1.5), x in S(-0.75):S(0.25):S(0.75), y in S(-0.75):S(0.25):S(0.75)
            z = complex(x, y)
            @test pFqdrummond((), (α, ), z) ≈ CS(pFq((), (T(α), ), CT(z))) atol=atol rtol=rtol
            @test pFqweniger((), (α, ), z) ≈ CS(pFq((), (T(α), ), CT(z))) atol=atol rtol=rtol
        end
    end
end

@testset "₂F₀" begin
    for (S, T) in ((Float32, Float64),
                   (Float64, BigFloat),
                   (BigFloat, BigFloat))
        CS = Complex{S}
        CT = Complex{T}
        atol = rtol = sqrt(eps(S))
        for α in S(-1.5):S(1.0):S(1.5), β in S(-1.5):S(1.0):S(1.5), z in S(-1.0):S(0.25):S(0.0)
            @test pFqdrummond((α, β), (), z) ≈ S(pFq((T(α), T(β)), (), T(z))) atol=atol rtol=rtol
            @test pFqweniger((α, β), (), z) ≈ S(pFq((T(α), T(β)), (), T(z))) atol=atol rtol=rtol
        end
        for α in S(-0.5):S(1.0):S(0.5), β in S(-0.5):S(1.0):S(0.5), x in S(-0.5):S(0.25):S(0.0), y in S(-0.5):S(0.25):S(0.0)
            z = complex(x, y)
            @test pFqdrummond((α, β), (), z) ≈ CS(pFq((T(α), T(β)), (), CT(z))) atol=atol rtol=rtol
            @test pFqweniger((α, β), (), z) ≈ CS(pFq((T(α), T(β)), (), CT(z))) atol=atol rtol=rtol
        end
    end
end

@testset "₁F₁" begin
    for (S, T) in ((Float32, Float64),
                   (Float64, BigFloat),
                   (BigFloat, BigFloat))
        CS = Complex{S}
        CT = Complex{T}
        atol = rtol = 1000eps(S)
        for α in S(-1.5):S(0.5):S(1.5), β in S(-1.5):S(1.0):S(1.5), z in S(-0.75):S(0.25):S(0.75)
            @test pFqdrummond((α, ), (β, ), z) ≈ S(pFq((T(α), ), (T(β), ), T(z))) atol=atol rtol=rtol
            @test pFqweniger((α, ), (β, ), z) ≈ S(pFq((T(α), ), (T(β), ), T(z))) atol=atol rtol=rtol
        end
        atol *= 2
        rtol *= 2
        for α in S(-1.5):S(0.5):S(1.5), β in S(-1.5):S(1.0):S(1.5), x in S(-0.75):S(0.25):S(0.75), y in S(-0.75):S(0.25):S(0.75)
            z = complex(x, y)
            @test pFqdrummond((α, ), (β, ), z) ≈ CS(pFq((T(α), ), (T(β), ), CT(z))) atol=atol rtol=rtol
            @test pFqweniger((α, ), (β, ), z) ≈ CS(pFq((T(α), ), (T(β), ), CT(z))) atol=atol rtol=rtol
        end
    end
end

@testset "₀F₂" begin
    for (S, T) in ((Float32, Float64),
                   (Float64, BigFloat),
                   (BigFloat, BigFloat))
        CS = Complex{S}
        CT = Complex{T}
        atol = rtol = 1000eps(S)
        for α in S(-1.5):S(1.0):S(1.5), β in S(-1.5):S(1.0):S(1.5), z in S(-0.75):S(0.25):S(0.75)
            @test pFqdrummond((), (α, β), z) ≈ S(pFq((), (T(α), T(β)), T(z))) atol=atol rtol=rtol
            @test pFqweniger((), (α, β), z) ≈ S(pFq((), (T(α), T(β)), T(z))) atol=atol rtol=rtol
        end
        atol *= 2
        rtol *= 2
        for α in S(-1.5):S(1.0):S(1.5), β in S(-1.5):S(1.0):S(1.5), x in S(-0.75):S(0.25):S(0.75), y in S(-0.75):S(0.25):S(0.75)
            z = complex(x, y)
            @test pFqdrummond((), (α, β), z) ≈ CS(pFq((), (T(α), T(β)), CT(z))) atol=atol rtol=rtol
            @test pFqweniger((), (α, β), z) ≈ CS(pFq((), (T(α), T(β)), CT(z))) atol=atol rtol=rtol
        end
    end
end

@testset "₂F₁" begin
    for (S, T) in ((Float32, Float64),
                   (Float64, BigFloat),
                   (BigFloat, BigFloat))
        CS = Complex{S}
        CT = Complex{T}
        atol = rtol = 1000eps(S)
        for α in S(-1.5):S(0.5):S(1.5), β in S(-1.5):S(0.5):S(1.5), γ in S(-1.5):S(1.0):S(1.5), z in S(-0.625):S(0.25):S(0.125)
            @test pFqdrummond((α, β), (γ, ), z) ≈ S(pFq((T(α), T(β)), (T(γ), ), T(z))) atol=atol rtol=rtol
            @test pFqweniger((α, β), (γ, ), z) ≈ S(pFq((T(α), T(β)), (T(γ), ), T(z))) atol=atol rtol=rtol
        end
        atol *= 2
        rtol *= 2
        for α in S(-1.5):S(0.5):S(1.5), β in S(-1.5):S(0.5):S(1.5), γ in S(-1.5):S(1.0):S(1.5), x in S(-0.375):S(0.25):S(0.125), y in S(-0.375):S(0.25):S(0.125)
            z = complex(x, y)
            @test pFqdrummond((α, β), (γ, ), z) ≈ CS(pFq((T(α), T(β)), (T(γ), ), CT(z))) atol=atol rtol=rtol
            @test pFqweniger((α, β), (γ, ), z) ≈ CS(pFq((T(α), T(β)), (T(γ), ), CT(z))) atol=atol rtol=rtol
        end
    end

    @testset "The points exp(±im*π/3)" begin
        for (S, T) in ((Float32, Float64),
                       (Float64, BigFloat))
            zS = exp(im*S(π)/3)
            zT = exp(im*T(π)/3)
            @test _₂F₁(1, 2, 3, zS) ≈ _₂F₁(1, 2, 3, zT)
            @test _₂F₁(1, 2, 3, conj(zS)) ≈ _₂F₁(1, 2, 3, conj(zT))
            @test _₂F₁(1, 2, 3, zS) == conj(_₂F₁(1, 2, 3, conj(zS)))
            @test pFq((1, 2+im), (S(3.5), ), zS) ≈ pFq((1, 2+im), (T(3.5), ), zT)
            @test pFq((1, 2+im), (S(3.5), ), conj(zS)) ≈ pFq((1, 2+im), (T(3.5), ), conj(zT))
        end
    end

    @testset "Integer arguments" begin
        @test _₂F₁(1, 0, 3, -1) ≡ _₂F₁(1.0, 0, 3, -1) ≡ 1.0
    end
end

@testset "₅F₄" begin # Tests the generic algorithms for pFqdrummond and pFqweniger
    for (S, T) in ((Float32, Float64),
                   (Float64, BigFloat),
                   (BigFloat, BigFloat))
        atol = rtol = sqrt(eps(S))
        @testset "Small argument" begin
            α = (1, 2, 3, 4, 5)
            β = (6, 7, 8, 9)
            z = eps(T)^2
            @test pFqdrummond(α, β, S(z)) ≈ S(pFq(α, β, z)) atol=atol rtol=rtol
            @test pFqweniger(α, β, S(z)) ≈ S(pFq(α, β, z)) atol=atol rtol=rtol
        end
        @testset "Terminating hypergeometrics, multiple dispatch" begin
            α = Real[1//1, 2, 3, 4, -5]
            β = Real[6//1, 7, 8, 9]
            for z in T(-5):T(0.5):T(5)
                @test pFqdrummond(α, β, S(z)) ≈ S(pFq(α, β, z)) atol=atol rtol=rtol
                @test pFqweniger(α, β, S(z)) ≈ S(pFq(α, β, z)) atol=atol rtol=rtol
            end
        end
        @testset "Generic algorithms" begin
            α = (T(1)/T(6), T(2)/T(6), T(3)/T(6), T(4)/T(6), T(5)/T(6))
            β = (T(1)/T(7), T(2)/T(7), T(3)/T(7), T(4)/T(7))
            for z in T(-5):T(0.5):T(0)
                @test pFqdrummond(S.(α), S.(β), S(z)) ≈ S(pFq(α, β, z)) atol=atol rtol=rtol
                @test pFqweniger(S.(α), S.(β), S(z)) ≈ S(pFq(α, β, z)) atol=atol rtol=rtol
            end
        end
    end
end

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

@testset "Warning symbols" begin
    @test pFq2string(Val(86420), Val(97531)) == "₈₆₄₂₀F₉₇₅₃₁"
end
