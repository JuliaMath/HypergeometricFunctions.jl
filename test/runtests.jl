using HypergeometricFunctions, Test
import LinearAlgebra: norm
import HypergeometricFunctions: _₂F₁general # not exported 

@testset "Hypergeometric Function tests" begin
    e = exp(1.0)
    regression_max_accumulated_error = 6.0e-14 # discovered by running the test
    for z in (.9rand(Float64,10), 10rand(ComplexF64, 10))
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
            println("This is the cumulative error for Case $j: ", error_accum)
            @test error_accum < regression_max_accumulated_error
            j += 1
        end
    end
end

