using ApproxFun, SIE, Base.Test


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
w=1/sqrt(1-x^2)
H=Hilbert()
B=ldirichlet(space(x))

for a in [sqrt(sqrt(5)-2)/2,1.,10.]
    L=H[w]+1/a/sqrt(1+a^2)*x
    u=[B,L]\[1.]
    usol = (1+a^2)/(x^2+a^2)
    @test norm(u-usol) <= eps(100/a)
end



Γ=Circle()∪Circle(0.5)
f=devec([Fun(z->z^(-1),Γ[1]),Fun(z->z,Γ[2])])
A=I-(f-1)*Cauchy(-1)
u=A\(f-1)
@test_approx_eq 1+cauchy(u,.1) 1
@test_approx_eq 1+cauchy(u,.8) 1/0.8
@test_approx_eq 1+cauchy(u,2.) 1





x = Fun(identity)
w = 1/sqrt(1-x^2)
H = Hilbert(space(w))
@test_approx_eq  (H[w]*exp(x))[.1] hilbert(w*exp(x))[.1]


x = Fun(identity)
w = sqrt(1-x^2)
H = Hilbert(space(w))
@test_approx_eq (H[w]*exp(x))[.1] hilbert(w*exp(x))[.1]

