using HypergeometricFunctions, Test
import LinearAlgebra: norm
import HypergeometricFunctions: _₂F₁general # not exported 
#import SingularIntegralEquations: mFn

const rtol = 1.0e-3

@testset "Hypergeometric Function tests" begin
  @testset "_₂F₁ vs _₂F₁general" begin
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

  @testset "_₂F₁ vs mpmath" begin
(a, b, c, result) = ([2.1692, 1.93762], [-0.0476108], -2.2931638546847504, -0.6599275400463818)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-0.400129, 1.91538], [-2.81784], 0.6811778271570272, 115.20573579649673)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-1.195, -2.34769], [-2.6397], 0.2615452629117505, 0.7281348361596178)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([1.13309, -0.367523], [-1.9268], -2.5986485146899376, 2.460363672158983)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-2.35708, 1.61642], [1.30338], -2.327580427222807, 22.835862384814842)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([0.954224, -2.44124], [0.732855], -0.9291840329978731, 6.477958868848867)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([2.66308, -1.57175], [2.61258], 0.951911371268817, -0.0014906095077265441)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([1.70608, 1.98557], [-2.47355], 0.7288912607897422, -172822.00668870765)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
  end

  @testset "_₃F₂ vs mpmath" begin
(a, b, c, result) = ([-1.82239, -0.520611, -0.351354], [0.120754, -0.927537], -0.22258191520734938, 0.10847716887590778)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-2.13579, 0.882367, 1.36492], [-0.698264, -2.4247], 0.9241926789225885, -12420.846337330078)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([0.630819, -1.0692, -2.17273], [-1.94888, 0.592168], -2.8764584839306466, 5.859296099858778)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-0.506505, 1.81628, 2.79278], [-2.20183, -2.83999], -0.45417732691014967, 3.654026483672525)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-0.680207, 1.2615, -1.06195], [-0.333544, -1.04341], -1.5760912049164775, -0.8718948204728744)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-0.0969323, -0.0938326, 2.55666], [-2.6199, 1.4186], 0.16573730123455022, 0.9989470589931153)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([0.8332, -0.851457, -0.564525], [0.579897, -1.93074], 0.8983750100832251, 1.9662567799536188)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
  end

  @testset "mFn vs mpmath" begin
(m, n) = (1, 1)
(a, b, c, result) = ([0.440755], [-2.51575], -2.7636058217262978, 1.3196127822630042)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-0.736412], [-2.46369], 1.1942629308763637, 1.4510642408138872)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([0.16622], [1.36032], 2.363692185412819, 1.5993489928773104)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([1.09305], [0.929982], 1.0663092391866908, 3.345547702290961)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-0.458657], [0.533445], 1.8289308422934223, -1.35371849169052)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([1.81437], [-0.648648], 1.1397965558600038, -46.443556364437065)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-0.204083], [-0.742601], 0.6665081386549785, 1.4463673236315062)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([0.360231], [1.67732], 2.7330356479094924, 2.4248370700541133)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([0.917765], [-1.05831], 0.4051850210490664, 4.199890925860252)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([1.55751], [0.713481], -0.4877473498368956, 0.25127947921492916)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(m, n) = (1, 2)
(a, b, c, result) = ([2.14012], [-2.44988, -1.92993], -1.1496121094103158, -19.210706713662876)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([0.10468], [2.66697, 2.41039], -1.9245285721545295, 0.9711585423906209)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-2.86873], [1.24933, -1.71116], -0.562614304959415, 0.5360567296241636)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-2.01274], [-2.74831, 1.93418], 1.2387159052950691, 1.5265045306760903)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([2.11422], [0.827965, 1.00587], 0.4844712020883666, 2.5035188859994624)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-1.4574], [-0.523284, 0.576072], -2.3589893153505828, -17.76330988948321)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-2.12378], [0.154787, -0.857248], 1.2260648133937222, -59.6305970582522)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([2.18923], [-2.86515, 0.304549], -0.3028639307184604, 1.9165619789597454)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([0.0917145], [-1.04493, 2.45165], 0.9383303605005469, 1.0956259127995325)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-0.374778], [-1.91622, -2.20963], 2.5837768019903433, 56.892149246264204)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(m, n) = (1, 3)
(a, b, c, result) = ([1.6433], [2.92957, -0.937415, -0.885855], -2.2561003004794333, 92.76061227189554)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-0.695529], [1.797, 0.516484, 1.36271], 1.0025792514731209, 0.440142960696577)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([1.83793], [-0.387405, -0.138869, -2.07554], 2.643824501428736, -6272.931477720918)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-1.12319], [-2.32979, 1.77287, -1.11033], 0.013172167352968778, 0.996780400598354)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-0.827685], [1.20337, -0.749147, 0.915409], -0.0005311620162200548, 0.9994672894815547)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-1.23053], [-1.141, -2.80502, -2.84556], -0.3606922506897261, 0.9567442665474625)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-1.45151], [2.7884, -0.0818559, 2.35224], 0.4974569685983159, 2.331917932076122)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-2.9264], [-1.52897, 1.33085, 2.27067], 2.5642144802930846, 3.5010763774333102)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-0.578959], [1.59537, 1.51563, 0.478805], -1.568855909062452, 1.7583353402171065)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([0.163702], [1.54918, 2.72876, 0.644636], -1.1945369617189758, 0.9313706492155626)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(m, n) = (2, 1)
(a, b, c, result) = ([-1.0564, 2.98783], [-0.303294], -2.8376924774009495, -32.95087177615518)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([1.52436, -0.0081617], [-0.283169], -0.5400346592629188, 0.9888954821464697)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-1.04946, 0.662651], [1.34829], -2.049532638789996, 2.083811851316413)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-0.580822, -2.277], [1.67567], -1.8571855227709224, -0.7537893805700949)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-0.410205, -1.31748], [1.82558], -1.8954126596178278, 0.40824104205380396)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([2.62166, 0.548989], [0.97727], -2.9234420984459444, 0.21452815708847717)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([0.496529, -0.861281], [2.93887], -0.1489573992576565, 1.0215928543032071)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(m, n) = (2, 2)
(a, b, c, result) = ([-0.623775, -1.5497], [2.3707, -1.76367], 2.2390692630191973, 0.40681977671392205)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-0.68464, -0.504461], [1.61667, -0.754288], 1.8496053892836195, 0.31609870171743226)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([0.950995, 1.753], [-0.309264, 1.64383], -0.07159849615385117, 1.211100536633986)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-2.64391, 0.593623], [-2.05458, 2.72044], 0.8204769103667688, 1.4112589566509304)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-0.919535, 0.665718], [-2.73075, 0.637213], 2.9937645610269086, 5.6315228169161875)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-0.151991, 1.27694], [-2.27906, -1.28123], -1.7697303173446897, 2.2974001581723296)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-2.38042, -0.0438634], [-2.36484, -0.0119172], 2.382381594304804, 36.876818209636276)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([1.95593, 0.400897], [0.782395, 0.609921], -0.40970974194028464, 0.49433380474362487)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([2.22909, 1.19391], [0.252987, -1.80782], -1.0380919197535223, -25.47692703711477)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([1.36612, 2.38078], [1.48043, -0.468035], 0.5694089934240094, -9.590144992510018)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(m, n) = (2, 3)
(a, b, c, result) = ([-0.339652, 1.68456], [2.48436, -2.7906, -2.90924], -0.5146019395545789, 1.0120705321211192)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([2.02812, 0.044901], [-2.57825, 1.35943, -2.7112], -0.934470150577039, 1.0049424096707185)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([1.03967, -0.506077], [0.712296, -0.572497, 1.0473], 2.2403798465641955, 6.5379635620966035)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-0.176333, 1.09293], [-0.498712, 2.26034, -1.73552], 2.7056329045109218, 3.712058614623512)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-0.0118928, -1.14223], [2.87683, -2.41144, -2.56972], 0.6326529220667685, 1.0004783006209943)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([1.65295, 1.79902], [-0.157296, -0.556574, 0.0678397], -2.951396835947348, -6712.636146441764)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-1.96683, -2.94163], [-1.98571, 1.71498, -1.52858], 2.5062568694008402, 3.5615591275296152)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-0.404805, 2.97669], [0.297451, -0.0231223, -0.578194], 0.1841449285549288, -79.31181066534661)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([0.67252, 2.33849], [1.82024, 0.121729, 1.31276], 2.9144896952488657, 43.492619145312425)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([1.08546, -0.293141], [-2.31911, -1.85707, -0.46119], 2.9061825898499647, -1200.8399960407003)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(m, n) = (3, 1)
(a, b, c, result) = ([-2.29252, 2.02566, -0.142973], [2.40018], 0.014343129555648648, 1.003939671478494)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-1.37482, 1.15773, 2.95406], [1.03428], -1.957076153703674, 16.64958394810834)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-1.62896, -1.33134, 0.887817], [2.23886], -0.6425686175372713, 0.4683206310399534)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-0.616014, 1.59686, -2.44411], [-2.59645], -0.9411688106771599, 1.5821016066720122)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(m, n) = (3, 2)
(a, b, c, result) = ([-1.98721, 0.0503148, 0.980717], [2.85255, -1.16388], 0.6486161733912801, 1.0394739225656187)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([1.73836, 2.44785, -1.34564], [1.41883, 1.91207], -0.588936842379483, 2.3913670710288426)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-2.89573, 1.38051, 0.102979], [1.68537, 2.17109], -1.339289847468002, 1.2203798316914103)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-2.4812, -2.98054, -1.78085], [1.30398, -1.75694], -2.0427920669019644, 3.6314918617117207)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([0.321834, -2.30457, -0.537795], [1.12969, 2.23388], -0.9848834358334244, 0.8352454691630404)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([2.89183, -2.58316, -0.134141], [-2.47993, 1.31591], -0.9897486110094627, 1.222684196322197)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(m, n) = (3, 3)
(a, b, c, result) = ([-0.921838, -0.19101, 2.91581], [-2.22718, -0.166386, 0.0240815], -1.5115132677881888, -109.8511998267201)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-1.48543, 0.625028, -0.151949], [1.74899, 2.43635, 2.93188], 1.8569592122119403, 1.0206114407906073)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-0.0544257, -1.62261, -1.03717], [0.402305, 2.00103, -2.05284], -0.3432596240775254, 0.9809562074458447)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([-2.55021, 0.874857, -2.64777], [1.54557, -1.72532, -0.719343], 1.7182953664316063, -55.92605023366547)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([0.144753, 1.1023, -0.812289], [-2.40491, -2.85321, 1.65655], 1.371042050287076, 0.7268084170518301)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([2.03127, -1.73788, 0.846854], [-0.769023, -1.55802, 0.0718722], 0.4358169961903138, -163.33548768303692)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([2.51962, -0.215367, -2.49022], [-2.51441, -1.91864, 0.37046], 1.8607306671835948, -543.1105622122955)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([1.80472, -2.47465, -1.52207], [-0.520596, -2.89333, -0.996009], -2.748332570461503, 16654.433680736995)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([2.4027, 0.291704, 2.94719], [-2.57265, 0.988027, -1.38968], -1.81109489555055, 80.90898175755316)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
(a, b, c, result) = ([2.25047, 2.76235, -2.59358], [-2.56857, -1.26227, 1.60436], -0.6346325745613317, 3.059595029498082)
@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
  end

end

