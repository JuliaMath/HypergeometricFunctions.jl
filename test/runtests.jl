using HypergeometricFunctions, Test
import LinearAlgebra: norm
import HypergeometricFunctions: _₂F₁general # not exported 
#import SingularIntegralEquations: mFn

const rtol = 1.0e-3
const NumberType = Float128

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
    (a, b, c, result) = NumberType[-2.15894, -1.41515], NumberType[-2.21784], NumberType(-0.9036697220161636), 2.4460583867912056
    @show a, b, c
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-1.65804, 2.30939], NumberType[-2.60282], NumberType(-0.3348500117337361), 0.6244788846687951
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[1.16832, -1.47593], NumberType[1.43303], NumberType(-0.5743816667545976), 1.7688091506986645
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[2.27829, 1.06092], NumberType[1.52085], NumberType(-2.146696498266405), 0.17945084553186508
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[1.39759, 2.1962], NumberType[1.2911], NumberType(0.3761207731832421), 3.047461664797781
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-1.53199, -2.61604], NumberType[2.7967], NumberType(0.7044850300486005), 2.0889142563689536
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[1.37039, 1.88608], NumberType[-4.61418e-02], NumberType(-2.877787131547951), 0.832725184088605
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-2.51538, -1.1265], NumberType[2.19666], NumberType(-0.5339175425492884), 0.32250554406579257
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[2.14893, 2.82534], NumberType[-5.01687e-01], NumberType(-2.939989251628186), -0.1054783365038995
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-3.56049e-01, 1.99076], NumberType[-2.94152], NumberType(0.40126941541382477), 4.660598174324774
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
  end

  @testset "_₃F₂ vs mpmath" begin
    (a, b, c, result) = NumberType[2.29443, -2.85876, 1.30521e-01], NumberType[2.18536, -1.1185e-01], NumberType(-2.2819453099409284), -45.15742677026391
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-2.93285, -2.96773, 1.12393], NumberType[-2.70198, 4.99942e-01], NumberType(-1.5378834160530541), 61.499801802257956
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[1.00525, -1.89826, 2.09855], NumberType[-6.97571e-01, -2.96297], NumberType(-2.550623073647634), 81.50691123919889
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-1.90897, -2.20448, -2.89402], NumberType[1.16712, 2.49403], NumberType(-2.4051617409405024), 14.37421350842237
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[2.61764, 2.38183, 7.19534e-01], NumberType[-6.94957e-03, -1.82867], NumberType(-1.194685621921007), -239.6084221295957
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-1.21187e-01, 2.64157, -2.25808], NumberType[-1.88609, -9.88632e-02], NumberType(-1.5454565225149501), 74.8140699891001
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-2.4202, 1.17118e-02, -2.70763], NumberType[-2.61018, -2.17417], NumberType(-2.107168165586075), 0.9140988500722098
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[6.26809e-01, 2.84814, -8.41167e-01], NumberType[-1.17618, 2.91723], NumberType(0.9941976483625536), -58.46114677039432
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-2.68225, 1.79677e-01, -2.03469], NumberType[1.02078, 1.20166], NumberType(-2.0541594107981287), 0.13350092077423076
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[1.00647, 5.25768e-01, 1.98715], NumberType[-9.63588e-01, -1.50051], NumberType(0.18393761428600364), -79.01874030197901
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
  end

  @testset "mFn vs mpmath" begin
    (m, n) = (1, 1)
    (a, b, c, result) = NumberType[1.70676], NumberType[1.89631], NumberType(-2.372374943987298), 0.13366001302730574
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-9.75991e-01], NumberType[-1.9617], NumberType(-1.5809434631994352), 0.3058096082553889
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-2.9531e-01], NumberType[-6.80749e-01], NumberType(1.1879408611195883), 2.6964843972328683
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[7.57579e-01], NumberType[1.6613], NumberType(-0.44500246962720436), 0.8238374413997872
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-2.68464], NumberType[5.84702e-01], NumberType(-1.683705273991436), 16.642877426610667
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[5.84799e-01], NumberType[-4.95155e-01], NumberType(0.69327050708339), -1.1695997248677032
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[2.3905], NumberType[2.22256], NumberType(-2.7298953220689492), 0.0441348741520008
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[1.2061], NumberType[1.22404], NumberType(-1.3694262069972267), 0.2615515763584455
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[2.30457], NumberType[6.57249e-01], NumberType(-1.495274392445984), -0.3517753320280511
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-1.14305], NumberType[1.8597], NumberType(0.7184172341572581), 0.5668200640707556
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (m, n) = (1, 2)
    (a, b, c, result) = NumberType[2.55657], NumberType[-2.78535, 2.73316], NumberType(2.374678426129718), -1.6941431709676191
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-1.40918], NumberType[-2.29006, 2.53045], NumberType(-2.325332487247478), 0.5098238857677238
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[1.86397], NumberType[-8.47673e-01, -2.93765], NumberType(-1.79085084929), 166.97974182727535
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[2.12172], NumberType[-1.98193, -8.86197e-01], NumberType(-0.06190360006129492), 1.1238742598512195
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-8.38017e-01], NumberType[-8.67518e-01, 1.38063], NumberType(-2.8062844431181593), 0.10860169312981847
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[9.94448e-01], NumberType[1.81165e-01, -1.17758], NumberType(0.016871288847900345), 0.9277236964018801
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[2.34928], NumberType[-8.10193e-01, 2.44249], NumberType(-1.8807793344726926), -3.023026855461745
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-2.00517], NumberType[-8.11401e-01, 2.98206], NumberType(1.6286957423368738), 0.8793244176558546
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-2.05028], NumberType[-2.15705, 9.00933e-01], NumberType(-1.7913314984294453), -0.1284823844186217
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[9.38309e-01], NumberType[2.01844, 1.03747], NumberType(-2.323332718985186), 0.276500302035838
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (m, n) = (1, 3)
    (a, b, c, result) = NumberType[-6.38817e-01], NumberType[1.1624, 1.75917, 2.45639], NumberType(-2.23623310602093), 1.2789362918406582
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[2.16614], NumberType[1.9022, 2.13265, -1.72594], NumberType(1.3734235655770184), 0.7788923045024129
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-2.74802], NumberType[4.30747e-01, -2.98604, 1.04665], NumberType(-0.9523498839275031), -0.6720509434647786
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-1.5942], NumberType[1.40658, -3.59324e-01, -2.69924], NumberType(1.610478224991608), -1.2032846564251822
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-1.34483], NumberType[-1.8497, -1.54209, 1.71444], NumberType(-1.3463891389166909), 1.3653721688721143
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-6.01003e-01], NumberType[-1.24284, 1.50623, -1.20916], NumberType(0.4405501422143554), 0.7942159209442613
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[2.64209], NumberType[-2.93569, -2.23575, -2.17681], NumberType(-2.7525480448608497), -1295.5127683549124
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-7.42972e-02], NumberType[-3.65231e-01, 2.9594, -1.89975], NumberType(1.750864068928636), 0.9931339381882528
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[4.4676e-01], NumberType[2.45553, 9.04256e-01, -2.94047], NumberType(-2.0411500839758077), 1.1551253034591646
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-2.73282], NumberType[-3.15481e-01, 1.81549e-01, -2.49563e-01], NumberType(2.0478579099425716), 664.494252671576
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (m, n) = (2, 1)
    (a, b, c, result) = NumberType[-1.78572, -2.1286], NumberType[-2.01773], NumberType(-1.4273191047800493), 4.432126882174454
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[6.61096e-01, 1.73908], NumberType[-2.40653e-02], NumberType(-1.1400005772799728), 10.41567183575274
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[2.45092, 1.02137], NumberType[-2.55196], NumberType(-0.5879698287201567), 1.3186506624438314
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[2.47837, -7.76862e-01], NumberType[-1.14999e-01], NumberType(0.41363823046832016), 10.080254420522024
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-2.64483, 8.40198e-01], NumberType[-6.65097e-02], NumberType(0.1718613580645023), 5.2308694058837295
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[9.43251e-01, 3.15523e-01], NumberType[2.92612], NumberType(-2.9748172138827806), 0.8302227587306925
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[1.28669, 6.5402e-01], NumberType[1.39542], NumberType(-2.31198875289835), 0.48076901781464365
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[1.85218, -2.8234], NumberType[-1.7853], NumberType(0.3448525521981125), 1.050407573846761
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[1.33061e-02, -1.98165], NumberType[1.60944], NumberType(-0.8255873480681988), 1.015648940440939
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[1.85238, 2.9345], NumberType[-5.2404e-01], NumberType(-2.3566581483108653), -0.28281694060671875
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (m, n) = (2, 2)
    (a, b, c, result) = NumberType[1.58121, 1.34768], NumberType[-8.79313e-01, -1.01424], NumberType(-2.270564255898779), 503.96940253864847
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-1.09043, 2.05], NumberType[1.0318, -2.27541], NumberType(2.544841297846053), -6.375629681812616
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-1.53558, -1.2336], NumberType[-1.42632, 2.09199], NumberType(-1.3750931765964216), 1.9263836655397244
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[1.82387, -2.96938], NumberType[1.03103, -1.66024], NumberType(-2.6717149592913105), 187.3211548956439
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-2.53547, 1.75006], NumberType[2.30864, -8.37104e-01], NumberType(-0.637179724256359), -4.4164591555407675
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[2.80981, 1.77346], NumberType[1.92165, 4.78668e-01], NumberType(1.5128235377941257), 53.20227469320951
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-1.28487, 2.26301], NumberType[9.48739e-01, 1.80404], NumberType(-0.9290540577429596), 2.6935759802046193
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[9.7987e-02, -2.58837], NumberType[-2.28646, -2.50159], NumberType(1.733756200030439), -0.4403727735119794
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[1.02025, 1.87127], NumberType[5.9682e-01, 1.77169], NumberType(-2.744375165569857), -0.19376051419417767
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-3.74486e-01, 7.62436e-01], NumberType[1.509, 1.86942], NumberType(1.6921522441697139), 0.801442160309535
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (m, n) = (2, 3)
    (a, b, c, result) = NumberType[2.70055, -6.50649e-01], NumberType[-8.78809e-02, 1.30714, -4.1776e-01], NumberType(1.9972707551580977), -189.29495098808098
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[1.91062, 1.11843], NumberType[2.06884e-01, -6.81032e-01, -5.99952e-01], NumberType(-0.14317207727168269), 6.331890865967105
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-2.30765, -1.78008], NumberType[-1.96917, -2.91966, -2.88729], NumberType(-0.5853524274698896), 1.092522249691883
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-1.39832, 7.42913e-01], NumberType[2.83979e-01, -2.9874, 2.97792], NumberType(-0.2717073564702357), 0.8893041773648375
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-1.38971, -1.74202e-01], NumberType[2.48285, -2.05977, -1.64466], NumberType(-2.205873022548032), 0.870814513352312
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-2.08571, -1.87755], NumberType[1.24674e-01, 2.72268, 2.83482], NumberType(2.6424658817019266), 12.596595185645997
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[2.20494, -2.74201e-01], NumberType[2.80949, -8.54445e-02, 2.63613], NumberType(-1.7631861567205362), -0.4374002687867758
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-1.38148, -5.4878e-01], NumberType[-1.03515, 1.25495, -1.81551], NumberType(2.9204917090691676), -9.10570431958492
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[3.61147e-01, 2.13107], NumberType[1.29714, 5.97651e-01, 4.52751e-02], NumberType(1.5453036548237513), 73.59876298822317
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-4.782e-01, -1.03659], NumberType[-6.90553e-01, -1.31267, 8.79193e-01], NumberType(-0.9355716977049893), 0.4420008917645215
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (m, n) = (3, 1)
    (a, b, c, result) = NumberType[-2.19918, 1.03617, 1.13665], NumberType[1.45005], NumberType(-0.9372659588697774), 4.539098436417541
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[2.50212, -2.43696e-01, 1.18427], NumberType[-1.53833], NumberType(-2.2341736065923996), 0.9552685731642808
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-5.49382e-01, -8.72184e-01, 2.83068], NumberType[8.20955e-01], NumberType(-1.8946549586973322), -1.9522027471390941
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[9.77175e-01, -1.96155, 2.47999], NumberType[6.94459e-01], NumberType(-2.867599864667291), 124.42232304126127
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-8.21767e-01, 2.91009, 1.1674], NumberType[-4.1248e-01], NumberType(-0.16052226332117892), 0.0588408557345908
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[2.85341, -5.85145e-01, -1.98007], NumberType[-1.78095], NumberType(-1.0306099186554394), 1.0704674243962757
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-1.67321, 8.84656e-01, 7.46252e-01], NumberType[-1.76112], NumberType(-0.29023315657615933), 0.8547273022303031
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-2.53615, -3.72469e-01, -2.22535e-01], NumberType[1.89038], NumberType(-2.7406079107304047), 1.4438768242845073
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-2.31352, 5.62164e-01, 2.52853], NumberType[1.32734], NumberType(-1.715537936406721), 20.658201472924887
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[1.59228, -2.98202, 2.68859], NumberType[1.59924], NumberType(-1.6860828667753758), 306.9296410280478
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (m, n) = (3, 2)
    (a, b, c, result) = NumberType[2.68772, 1.06722, 2.74101], NumberType[2.15777, 1.86408], NumberType(0.5740222352798141), 6.00503559065529
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[1.27397, -2.93904, -2.8797], NumberType[-3.96577e-02, 1.85523], NumberType(-1.9417537853168843), -356.291462713288
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[3.26026e-01, -2.34176, 5.89925e-01], NumberType[-2.27193e-01, -1.78432], NumberType(-1.6453518987131943), -15.914036829015224
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-2.00464e-01, 1.71386, 2.95797], NumberType[3.51295e-01, 2.53807], NumberType(0.7617252851953022), -2.8218075896818884
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-1.53999, -1.8412, 2.39342], NumberType[-2.70016, 1.82632], NumberType(-2.99673925055985), 7.061478869649158
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-2.63468, 1.84885, -1.90141e-01], NumberType[1.54065, 1.65157], NumberType(-2.848436109585445), -1.1011889132740955
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-9.1968e-01, -1.80394, -2.6924], NumberType[1.02727, 1.35637], NumberType(-1.6546515047853756), 6.203285968923722
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-1.96263, 4.87085e-01, 2.15902e-02], NumberType[2.26626, 2.34871], NumberType(-0.3587246821544561), 1.0014243050836682
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-9.46665e-01, 2.58564, -1.76849], NumberType[1.83186, 2.30154], NumberType(0.682805982528933), 1.697228401658023
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[1.38037e-01, -2.9377, -1.6562], NumberType[-2.77071, 2.08648], NumberType(-2.0042763557118337), 1.286279026974602
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (m, n) = (3, 3)
    (a, b, c, result) = NumberType[2.5713, 2.28672, 1.12178], NumberType[-1.22312, 2.90957e-01, -9.79813e-01], NumberType(0.5491289379969375), -297384.03887138446
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-1.39551, -2.79561, 2.41967], NumberType[6.26155e-01, -2.38208, -1.71469], NumberType(-0.8943871721336043), -4.4491341838866045
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[9.67741e-01, 8.42717e-01, 2.5836], NumberType[-2.82689, 1.67769, 2.0795], NumberType(-2.2875369855870535), 0.3616687899341795
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[1.62908, 1.7015, -2.6736e-01], NumberType[-1.95855, 1.54934, 2.55148], NumberType(-2.958396660915621), 2.450601906472814
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[-2.34528, 1.9183, 1.17128], NumberType[1.62052, 1.55863, -2.80072], NumberType(1.8826051923616278), 4.227776320275677
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[2.68063, -2.4882, -3.8737e-01], NumberType[-2.92169, -5.70777e-01, -4.85122e-02], NumberType(-2.794781130629891), -1926.403901528597
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[6.38698e-01, 2.39734e-01, -2.52872], NumberType[-9.43027e-01, -8.19748e-01, 7.65437e-01], NumberType(-1.1012032003207244), 90.48000665987023
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[5.59714e-01, 6.53033e-01, 2.94006], NumberType[-1.8711, -1.42919, 6.86401e-02], NumberType(0.2753544371339096), 284.03115461294715
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[2.92603, 2.2942, 1.77154], NumberType[5.82292e-01, -2.59065, -2.7667], NumberType(-0.10310736306314716), 2.0875854065105055
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
    (a, b, c, result) = NumberType[5.01323e-01, -9.76258e-01, -1.60155], NumberType[2.75106, -8.0503e-01, 1.85558], NumberType(-1.6287125316657916), 1.3131884579930222
    @test mFn(a, b, c) ≈ result atol=eps() rtol=rtol
  end

end

