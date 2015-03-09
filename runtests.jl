using ApproxFun, SIE, Base.Test

println("Hilbert test")

x=Fun(identity)
f=(exp(x)/(sqrt(1-x)*sqrt(x+1)))
@test_approx_eq (Hilbert(f|>space,0)*f)[.1] (-0.8545003781055088)
@test_approx_eq (Hilbert(0)*f)[.1] (-0.8545003781055088)
@test_approx_eq (Hilbert()*f)[.1] 1.1404096104609646386

x=Fun(identity,[-1,2])
f=(exp(x)/(sqrt(2-x)*sqrt(x+1)))
@test_approx_eq (Hilbert(f|>space,0)*f)[.1] 0.49127801561694168644
@test_approx_eq (Hilbert(0)*f)[.1] 0.49127801561694168644
@test_approx_eq (Hilbert()*f)[.1] 1.6649936695644078289



x=Fun(identity)
f=(exp(x)*(sqrt(1-x)*sqrt(x+1)))
@test_approx_eq (Hilbert()*f)[.1] 0.43723982258866913063

x=Fun(identity,[-1,2])
f=(exp(x)*(sqrt(2-x)*sqrt(x+1)))
@test_approx_eq (Hilbert()*f)[.1] 2.1380903070701673244

x=Fun(identity)
d=domain(x)
w=1/sqrt(1-x^2)
H=Hilbert(d)
B=ldirichlet(d)

for a in [sqrt(sqrt(5)-2)/2,1.,10.]
    L=H[w]+1/a/sqrt(1+a^2)*x
    u=[B,L]\[1.]
    usol = (1+a^2)/(x^2+a^2)
    @test norm(u-usol) <= eps(100/a)
end


x = Fun(identity)
w = 1/sqrt(1-x^2)
H = Hilbert(space(w))
@test_approx_eq  (H[w]*exp(x))[.1] hilbert(w*exp(x))[.1]


x = Fun(identity)
w = sqrt(1-x^2)
H = Hilbert(space(w))
@test_approx_eq (H[w]*exp(x))[.1] hilbert(w*exp(x))[.1]


println("Stieltjes test")

ds1 = JacobiWeight(-.5,-.5,ApproxFun.ChebyshevDirichlet{1,1}())
ds2 = JacobiWeight(-.5,-.5,Chebyshev())
rs = Chebyshev([2.,4.+3im])
f1 = Fun(x->exp(x)/sqrt(1-x^2),ds1)
f2 = Fun(x->exp(x)/sqrt(1-x^2),ds2)
S = Stieltjes(ds1,rs)

z = 3.+1.5im
@test_approx_eq (S*f1)[z] stieltjes(f2,z) #val,err = quadgk(x->f1[x]./(z-x),-1.,1.)
# Operator 1.1589646343327578 - 0.7273679005911196im
# Function 1.1589646343327455 - 0.7273679005911283im


ds1 = JacobiWeight(.5,.5,Ultraspherical{1}())
ds2 = JacobiWeight(.5,.5,Chebyshev())
rs = Chebyshev([2.,4.])
f1 = Fun(x->exp(x)*sqrt(1-x^2),ds1)
f2 = Fun(x->exp(x)*sqrt(1-x^2),ds2)
S = Stieltjes(ds1,rs)

z = 3.
@test_approx_eq (S*f1)[z] stieltjes(f2,z) #val,err = quadgk(x->f1[x]./(z-x),-1.,1.;reltol=eps())
# Operator 0.6616422557285478 + 0.0im
# Function 0.661642255728541 - 0.0im

println("Stieltjes integral test")

ds1 = JacobiWeight(-.5,-.5,ApproxFun.ChebyshevDirichlet{1,1}())
ds2 = JacobiWeight(-.5,-.5,Chebyshev())
rs = Chebyshev([2.,4.])
f1 = Fun(x->exp(x)/sqrt(1-x^2),ds1)
f2 = Fun(x->exp(x)/sqrt(1-x^2),ds2)
S = Stieltjes(ds1,rs,0)

z = 3.
@test_approx_eq (S*f1)[z] SIE.stieltjesintegral(f2,z)
# Operator 3.6322473044237698 + 0.0im
# Function 3.6322473044237515

ds1 = JacobiWeight(.5,.5,Ultraspherical{1}())
ds2 = JacobiWeight(.5,.5,Chebyshev())
rs = Chebyshev([2.0,4.0])
f1 = Fun(x->exp(x)*sqrt(1-x^2),ds1)
f2 = Fun(x->exp(x)*sqrt(1-x^2),ds2)
S = Stieltjes(ds1,rs,0)

z = 3.0
@test_approx_eq (S*f1)[z] SIE.stieltjesintegral(f2,z)
# Operator 1.7772163062194861 + 0.0im
# Function 1.7772163062194637

println("Cauchy test")

Γ=Circle()∪Circle(0.5)
f=depiece([Fun(z->z^(-1),Γ[1]),Fun(z->z,Γ[2])])
A=I-(f-1)*Cauchy(-1)
u=A\(f-1)
@test_approx_eq 1+cauchy(u,.1) 1
@test_approx_eq 1+cauchy(u,.8) 1/0.8
@test_approx_eq 1+cauchy(u,2.) 1

c1=0.5+0.1;r1=3.;
c2=-0.1+.2im;r2=0.3;
d1=Circle(c1,r1)
d2=Circle(c2,r2)
z=Fun(identity,d2);
C=Cauchy(Space(d1),Space(d2))
@test norm((C*Fun(exp,d1)-Fun(exp,d2)).coefficients)<100eps()

C2=Cauchy(Space(d2),Space(d1))
@test norm((C2*Fun(z->exp(1/z)-1,d2)+Fun(z->exp(1/z)-1,d1)).coefficients)<10000eps()

c1=0.1+.1im;r1=.4;
c2=-2.+.2im;r2=0.3;
d1=Circle(c1,r1)
d2=Circle(c2,r2)
@test norm((Cauchy(d1,d2)*Fun(z->exp(1/z)-1,d1)+Fun(z->exp(1/z)-1,d2)).coefficients)<200eps()


println("Arc test")

a=Arc(0.,1.,0.,π/2)
ζ=Fun(identity,a)
f=Fun(exp,a)*sqrt(abs((ζ-1)*(ζ-im)))
z=.1+.2im
@test_approx_eq cauchy(f,z) sum(f/(ζ-z))/(2π*im)
z=exp(.1im)
@test_approx_eq hilbert(f,z) im*(cauchy(+1,f,z)+cauchy(-1,f,z))


println("Functional test")
z=.1+.2im
x=Fun(identity)
f=exp(x)*sqrt(1-x^2)
@test_approx_eq Stieltjes(space(f),z)*f stieltjes(f,z)

a=Arc(0.,1.,0.,π/2)
ζ=Fun(identity,a)
f=Fun(exp,a)*sqrt(abs((ζ-1)*(ζ-im)))
H=Hilbert()
z=exp(.1im)
@test_approx_eq (H*f)[z] hilbert(f,z)


println("Logkernel test")

a=1.0;b=2.0
d=Interval(a,b)
z=Fun(d)
f=real(exp(z)/(sqrt(z-a)*sqrt(b-z)))
S=space(f)
x=4.0+2im;
@test_approx_eq linesum(f*log(abs(x-z))) logkernel(f,x)

a=1.0;b=2.0+im
d=Interval(a,b)
z=Fun(d)
f=real(exp(z)/(sqrt(z-a)*sqrt(b-z)))
S=space(f)
x=4.0+2im;
@test_approx_eq linesum(f*log(abs(x-z))) logkernel(f,x)

a=1.0;b=2.0+im
d=Interval(a,b)
z=Fun(d)
f=real(exp(z)*(sqrt(z-a)*sqrt(b-z)))
S=space(f)
x=4.0+2im;
@test_approx_eq linesum(f*log(abs(x-z))) logkernel(f,x)


a=1.0;b=2.0
d=Interval(a,b)
z=Fun(d)
f=real(exp(z)/(sqrt(z-a)*sqrt(b-z)))
x=1.5
@test_approx_eq (Hilbert(space(f),0)*f)[x] logkernel(f,x)/π

f=real(exp(z)*(sqrt(z-a)*sqrt(b-z)))
x=1.5
@test_approx_eq (Hilbert(space(f),0)*f)[x] logkernel(f,x)/π

a=1.0;b=2.0+im
d=Interval(a,b)
z=Fun(d)
f=real(exp(z)/(sqrt(z-a)*sqrt(b-z)))
x=1.5+0.5im
H=Hilbert(space(f),0)
@test_approx_eq (H*f)[x] logkernel(f,x)/π


f=real(exp(z)*(sqrt(z-a)*sqrt(b-z)))
@test_approx_eq (Hilbert(space(f),0)*f)[x] logkernel(f,x)/π

a=1.0;b=2.0+im
d=Interval(a,b)
z=Fun(d)
f=real(exp(z)/(sqrt(z-a)*sqrt(b-z)))
S=JacobiWeight(-0.5,-0.5,ChebyshevDirichlet{1,1}(d))
H=OffHilbert(S,Chebyshev([3,4]),0)
@test_approx_eq (H*f)[3.5] logkernel(f,3.5)/π

H=OffHilbert(S,Chebyshev([3,4.0+im]),0)
@test_approx_eq (H*f)[3.5+0.5im] logkernel(f,3.5+0.5im)/π


## Circle

d=Circle(0.2,3.0)
S=Fourier(d)
ζ=Fun(d)
f=real(ζ+1/(ζ-0.1))
z=0.1+0.1im;@test_approx_eq linesum(log(abs(ζ-z))*f) logkernel(f,z)
z=5.0+0.1im;@test_approx_eq linesum(log(abs(ζ-z))*f) logkernel(f,z)

d=Circle(0.2,3.0)
S=Fourier(d)
H=Hilbert(S,0)
ζ=Fun(d)
f=real(ζ+1/(ζ-0.1))
z=0.2+3im;@test_approx_eq (H*f)[z] logkernel(f,z)/π


println("KernelFun test")

include("KernelFunTest.jl")

println("Example test")

include("ExamplesTest.jl")
