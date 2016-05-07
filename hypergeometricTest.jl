using ApproxFun, SingularIntegralEquations, Base.Test

import SingularIntegralEquations: _₂F₁general

println()
println("Testing Hypergeometric series ₂F₁(a,b;c;z) against known transcendental cases")
println()

for z in (.9rand(Float64,10),10rand(Complex128,10))
    j=1
    for (a,b,c) in ((√2/e,1.3,1.3),(1.2,√3,1.2),
                    (-.4,.4,.5),(-.3,1.3,.5),(.35,.85,.5),
                    (.5,.5,1.5),(1.,1.,1.5),(.5,1.,1.5),(.3,.7,1.5),(.7,1.3,1.5),(.35,.85,1.5),
                    (1.,1.,2.),(3.,1.,2.),(-2.,1.,2.),(-3.,1.,2.),(1.,-4.,2.))
        normj = 0.0
        for zi in z
            twoFone = _₂F₁(a,b,c,zi)
            aa,bb,cc = big(a)+400eps(Float64),big(b)+700eps(Float64),big(c)+900eps(Float64)
            twoFonegeneral = convert(Complex{Float64},_₂F₁general(aa,bb,cc,big(zi)))
            norm(twoFone/twoFonegeneral-1) > sqrt(eps()) && println("This is ₂F₁($a,$b;$c;zi) - ₂F₁general($a,$b;$c;zi): ",norm(twoFone/twoFonegeneral-1),"   ",twoFone,"   ",twoFonegeneral,"   ",isfinite(twoFone),"   ",isfinite(twoFonegeneral)," this is zi: ",zi)
            normj += Float64(norm(twoFone/twoFonegeneral-1))
        end
        println("This is the cumulative error for Case $j: ",normj)
        j+=1
    end
end
