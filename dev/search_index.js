var documenterSearchIndex = {"docs":
[{"location":"#HypergeometricFunctions.jl-Documentation","page":"HypergeometricFunctions.jl Documentation","title":"HypergeometricFunctions.jl Documentation","text":"","category":"section"},{"location":"#Introduction","page":"HypergeometricFunctions.jl Documentation","title":"Introduction","text":"","category":"section"},{"location":"","page":"HypergeometricFunctions.jl Documentation","title":"HypergeometricFunctions.jl Documentation","text":"HypergeometricFunctions.jl provides a numerical computation of generalized hypergeometric functions. The main exported function and recommended interface is pFq, but there are a few others for specialists' convenience.","category":"page"},{"location":"","page":"HypergeometricFunctions.jl Documentation","title":"HypergeometricFunctions.jl Documentation","text":"pFq","category":"page"},{"location":"#HypergeometricFunctions.pFq","page":"HypergeometricFunctions.jl Documentation","title":"HypergeometricFunctions.pFq","text":"pFq(α, β, z)\n\nCompute the generalized hypergeometric function, defined by\n\n_pF_q(α β z) = sum_k=0^infty dfrac(alpha_1)_kcdots(alpha_p)_k(beta_1)_kcdots(beta_q)_kdfracz^kk\n\nwhere the series converges and elsewhere by analytic continuation.\n\nExternal links: DLMF, Wikipedia.\n\nExamples\n\njulia> pFq((), (), 0.1) # ≡ exp(0.1)\n1.1051709180756477\n\njulia> pFq((0.5, ), (), 1.0+0.001im) # ≡ exp(-0.5*log1p(-1.0-0.001im))\n22.360679774997894 + 22.36067977499789im\n\njulia> pFq((), (1.5, ), -π^2/4) # A root of a spherical Bessel function\n4.042865030283967e-17\n\njulia> pFq((), (1.5, ), -big(π)^2/4) # In extended precision\n8.674364372518869408017614476652675180406967418943475242812160199356160822272727e-78\n\njulia> pFq((1, ), (2, ), 0.01) # ≡ expm1(0.01)/0.01\n1.0050167084168058\n\njulia> pFq((1/3, ), (2/3, ), -1000) # A confluent hypergeometric with large argument\n0.050558053946448855\n\njulia> pFq((1, 2), (4, ), 1) # a well-poised ₂F₁\n2.9999999999999996\n\njulia> pFq((1, 2+im), (3.5, ), exp(im*π/3)) # ₂F₁ at that special point in ℂ\n0.6786952632946592 + 0.4523504929285015im\n\njulia> pFq((1, 2+im), (3.5, ), exp(im*big(π)/3)) # More digits, you say?\n0.6786952632946589823300834090168381068073515492901393549193461972311801512528996 + 0.4523504929285013648194489713901658143893464679689810112119412310631860619947939im\n\njulia> pFq((1, 2+im, 2.5), (3.5, 4), exp(im*π/3)) # ₃F₂ because why not\n0.8434434031615691 + 0.34175507615463174im\n\njulia> pFq((1, 2+im, 2.5), (3.5, 4), exp(im*big(π)/3)) # Also in extended precision\n0.8434434031615690763389963048175253868863156451003855955719081209861492349268002 + 0.3417550761546319732614495656712509723030350666571102474299311122586948108413206im\n\njulia> pFq((1, 1), (), -1) # A divergent series\n0.5963473623231942\n\njulia> pFq((1, 1), (), -big(1))\n0.5963473623231940743410784993692793760741778601525487815734849104823272191142015\n\n\n\n\n\n","category":"function"},{"location":"#Complex-phase-portraits","page":"HypergeometricFunctions.jl Documentation","title":"Complex phase portraits","text":"","category":"section"},{"location":"","page":"HypergeometricFunctions.jl Documentation","title":"HypergeometricFunctions.jl Documentation","text":"Broadly speaking, there are three classes of generalized hypergeometric functions: when ple q they are entire functions of the complex variable z; when p = q+1, they are analytic functions in the cut plane mathbbCsetminus1infty); and, when p  q+1, they are analytic functions in the cut plane mathbbCsetminus0infty).","category":"page"},{"location":"","page":"HypergeometricFunctions.jl Documentation","title":"HypergeometricFunctions.jl Documentation","text":"Examples of each of these classes are illustrated over leftzinmathbbC  -10Re z10 -10Im z10right with complex phase portraits, a beautiful tool in computational complex analysis.","category":"page"},{"location":"","page":"HypergeometricFunctions.jl Documentation","title":"HypergeometricFunctions.jl Documentation","text":"using ComplexPhasePortrait, HypergeometricFunctions, Images\nx = range(-10, stop=10, length=300)\ny = range(-10, stop=10, length=300)\nz = x' .+ im*y\n\nimport Logging # To avoid printing warnings\nLogging.with_logger(Logging.SimpleLogger(Logging.Error)) do\n    img = portrait(map(z->pFq((), (), z), z), ctype = \"nist\")\n    save(\"0F0.png\", img)\n    img = portrait(map(z->pFq((), (1.0, ), z), z), ctype = \"nist\")\n    save(\"0F1.png\", img)\n    img = portrait(map(z->pFq((0.5, ), (0.75, ), z), z), ctype = \"nist\")\n    save(\"1F1.png\", img)\n    img = portrait(map(z->pFq((3.5+7.5im, ), (), z), z), ctype = \"nist\")\n    save(\"1F0.png\", img)\n    img = portrait(map(z->pFq((1.0, 3.5+7.5im), (0.75, ), z), z), ctype = \"nist\")\n    save(\"2F1.png\", img)\n    img = portrait(map(z->pFq((1.0, 1.5+7.5im), (), z), z), ctype = \"nist\")\n    save(\"2F0.png\", img)\nend\nnothing # hide","category":"page"},{"location":"","page":"HypergeometricFunctions.jl Documentation","title":"HypergeometricFunctions.jl Documentation","text":"p\\q 0 1\n0 (Image: ₀F₀) (Image: ₀F₁)\n1 (Image: ₁F₀) (Image: ₁F₁)\n2 (Image: ₂F₀) (Image: ₂F₁)","category":"page"},{"location":"#Library","page":"HypergeometricFunctions.jl Documentation","title":"Library","text":"","category":"section"},{"location":"","page":"HypergeometricFunctions.jl Documentation","title":"HypergeometricFunctions.jl Documentation","text":"_₁F₁\n_₂F₁\n_₃F₂","category":"page"},{"location":"#HypergeometricFunctions._₁F₁","page":"HypergeometricFunctions.jl Documentation","title":"HypergeometricFunctions._₁F₁","text":"Compute Kummer's confluent hypergeometric function ₁F₁(a, b, z).\n\n\n\n\n\n","category":"function"},{"location":"#HypergeometricFunctions._₂F₁","page":"HypergeometricFunctions.jl Documentation","title":"HypergeometricFunctions._₂F₁","text":"Compute the Gauss hypergeometric function ₂F₁(a, b, c, z).\n\n\n\n\n\n","category":"function"},{"location":"#HypergeometricFunctions._₃F₂","page":"HypergeometricFunctions.jl Documentation","title":"HypergeometricFunctions._₃F₂","text":"Compute the generalized hypergeometric function ₃F₂(a₁, 1, 1, b₁, 2, z).\n\n\n\n\n\nCompute the generalized hypergeometric function ₃F₂(a₁, a₂, a₃, b₁, b₂; z).\n\n\n\n\n\n","category":"function"},{"location":"#Internals","page":"HypergeometricFunctions.jl Documentation","title":"Internals","text":"","category":"section"},{"location":"","page":"HypergeometricFunctions.jl Documentation","title":"HypergeometricFunctions.jl Documentation","text":"HypergeometricFunctions.M\nHypergeometricFunctions.U\nHypergeometricFunctions._₂F₁positive\nHypergeometricFunctions._₂F₁general\nHypergeometricFunctions._₂F₁general2\nHypergeometricFunctions.pFqdrummond\nHypergeometricFunctions.pFqweniger\nHypergeometricFunctions.pFqcontinuedfraction\nHypergeometricFunctions.pochhammer\nHypergeometricFunctions.@clenshaw\nHypergeometricFunctions.@lanczosratio\nHypergeometricFunctions.G\nHypergeometricFunctions.P","category":"page"},{"location":"#HypergeometricFunctions.M","page":"HypergeometricFunctions.jl Documentation","title":"HypergeometricFunctions.M","text":"Compute Kummer's confluent hypergeometric function M(a, b, z) = ₁F₁(a, b, z).\n\n\n\n\n\n","category":"function"},{"location":"#HypergeometricFunctions.U","page":"HypergeometricFunctions.jl Documentation","title":"HypergeometricFunctions.U","text":"Compute Tricomi's confluent hypergeometric function U(a, b, z) ∼ z⁻ᵃ ₂F₀((a, a-b+1), (), -z⁻¹).\n\n\n\n\n\n","category":"function"},{"location":"#HypergeometricFunctions._₂F₁positive","page":"HypergeometricFunctions.jl Documentation","title":"HypergeometricFunctions._₂F₁positive","text":"Compute the Gauss hypergeometric function ₂F₁(a, b, c, z) with positive parameters a, b, and c and argument 0 ≤ z ≤ 1. Useful for statisticians.\n\n\n\n\n\n","category":"function"},{"location":"#HypergeometricFunctions._₂F₁general","page":"HypergeometricFunctions.jl Documentation","title":"HypergeometricFunctions._₂F₁general","text":"Compute the Gauss hypergeometric function ₂F₁(a, b, c, z) with general parameters a, b, and c. This polyalgorithm is designed based on the paper\n\nN. Michel and M. V. Stoitsov, Fast computation of the Gauss hypergeometric function with all its parameters complex with application to the Pöschl–Teller–Ginocchio potential wave functions, Comp. Phys. Commun., 178:535–551, 2008.\n\n\n\n\n\n","category":"function"},{"location":"#HypergeometricFunctions._₂F₁general2","page":"HypergeometricFunctions.jl Documentation","title":"HypergeometricFunctions._₂F₁general2","text":"Compute the Gauss hypergeometric function ₂F₁(a, b, c, z) with general parameters a, b, and c. This polyalgorithm is designed based on the review\n\nJ. W. Pearson, S. Olver and M. A. Porter, Numerical methods for the computation of the confluent and Gauss hypergeometric functions, Numer. Algor., 74:821–866, 2017.\n\n\n\n\n\n","category":"function"},{"location":"#HypergeometricFunctions.pFqdrummond","page":"HypergeometricFunctions.jl Documentation","title":"HypergeometricFunctions.pFqdrummond","text":"pFqdrummond(α, β, z; kmax)\n\nCompute the generalized hypergeometric function pFq by rational approximations of type (k, k) generated by Drummond's sequence transformation described in\n\nR. M. Slevinsky, Fast and stable rational approximation of generalized hypergeometric functions, Numer. Algor., 98:587–624, 2025.\n\n\n\n\n\n","category":"function"},{"location":"#HypergeometricFunctions.pFqweniger","page":"HypergeometricFunctions.jl Documentation","title":"HypergeometricFunctions.pFqweniger","text":"pFqweniger(α, β, z; kmax, γ = 2)\n\nCompute the generalized hypergeometric function pFq by rational approximations of type (k, k) generated by a factorial Levin-type sequence transformation described in\n\nR. M. Slevinsky, Fast and stable rational approximation of generalized hypergeometric functions, Numer. Algor., 98:587–624, 2025.\n\n\n\n\n\n","category":"function"},{"location":"#HypergeometricFunctions.pFqcontinuedfraction","page":"HypergeometricFunctions.jl Documentation","title":"HypergeometricFunctions.pFqcontinuedfraction","text":"Compute the generalized hypergeometric function pFq(α, β, z) by continued fraction.\n\n\n\n\n\n","category":"function"},{"location":"#HypergeometricFunctions.pochhammer","page":"HypergeometricFunctions.jl Documentation","title":"HypergeometricFunctions.pochhammer","text":"Pochhammer symbol (x)_n = fracGamma(x+n)Gamma(x) for the rising factorial.\n\n\n\n\n\n","category":"function"},{"location":"#HypergeometricFunctions.@clenshaw","page":"HypergeometricFunctions.jl Documentation","title":"HypergeometricFunctions.@clenshaw","text":"@clenshaw(x, c...)\n\nEvaluate the Chebyshev polynomial series displaystyle sum_k=1^N ck T_k-1(x) by the Clenshaw algorithm.\n\nExternal links: DLMF, Wikipedia.\n\nExamples\n\njulia> HypergeometricFunctions.@clenshaw(1, 1, 2, 3)\n6\n\njulia> HypergeometricFunctions.@clenshaw(0.5, 1, 2, 3)\n0.5\n\n\n\n\n\n","category":"macro"},{"location":"#HypergeometricFunctions.@lanczosratio","page":"HypergeometricFunctions.jl Documentation","title":"HypergeometricFunctions.@lanczosratio","text":"@lanczosratio(z, ϵ, c₀, c...)\n\nEvaluate dfracdisplaystyle sum_k=0^N-1 fracck+1(z+k)(z+k+epsilon)displaystyle c_0 + sum_k=0^N-1 fracck+1z+k.\n\nThis ratio is used in the Lanczos approximation of logfracGamma(z+epsilon)Gamma(z) in\n\nN. Michel and M. V. Stoitsov, Fast computation of the Gauss hypergeometric function with all its parameters complex with application to the Pöschl–Teller–Ginocchio potential wave functions, Comp. Phys. Commun., 178:535–551, 2008.\n\n\n\n\n\n","category":"macro"},{"location":"#HypergeometricFunctions.G","page":"HypergeometricFunctions.jl Documentation","title":"HypergeometricFunctions.G","text":"Compute the function dfracfrac1Gamma(z)-frac1Gamma(z+epsilon)epsilon by the method dscribed in\n\nN. Michel and M. V. Stoitsov, Fast computation of the Gauss hypergeometric function with all its parameters complex with application to the Pöschl–Teller–Ginocchio potential wave functions, Comp. Phys. Commun., 178:535–551, 2008.\n\n\n\n\n\n","category":"function"},{"location":"#HypergeometricFunctions.P","page":"HypergeometricFunctions.jl Documentation","title":"HypergeometricFunctions.P","text":"Compute the function dfrac(z+epsilon)_m-(z)_mepsilon by the method dscribed in\n\nN. Michel and M. V. Stoitsov, Fast computation of the Gauss hypergeometric function with all its parameters complex with application to the Pöschl–Teller–Ginocchio potential wave functions, Comp. Phys. Commun., 178:535–551, 2008.\n\n\n\n\n\n","category":"function"}]
}
