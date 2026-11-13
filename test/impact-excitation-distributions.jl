# Standalone: julia --project=. test/impact-excitation-distributions.jl
using Test, JenaAtomicCalculator, QuadGK, SpecialFunctions, Logging
const NMI = JenaAtomicCalculator.ImpactExcitation
const nmempty = Level()
nmlevel(index,E) = Level(AngularJ64(0),nmempty.M,Basics.plus,index,E,nmempty.relativeOcc,nmempty.hasStateRep,nmempty.basis,nmempty.mc)
const nminitial = nmlevel(1,0.)
nmline(f,E,omega;shift=0.) = NMI.Line(nminitial,f,E+f.energy-shift,E,0.,omega,Float64[],Float64[],NMI.Channel[],0.)
const nmT = Defaults.convertUnits("temperature: from atomic to Kelvin",1.)
nmsettings(T=nmT) = NMI.Settings(NMI.Settings();temperatures=[T])
function nmrate(lines;model=:Kappa,kap=3.,T=nmT)
    with_logger(NullLogger()) do
        only(NMI.computeRateCoefficients(lines,nmsettings(T);distribution=model,kappa=kap))
    end
end

@testset "Non-Maxwellian excitation rates" begin
    @testset "PDF normalization and mean-energy convention" begin
        for model in [:Maxwellian,:Kappa,:Druyvesteyn], kap in (model==:Kappa ? [1.6,3.,1.e6] : [3.]), T in [1.e3,1.e5]
            theta=Defaults.convertUnits("temperature: from Kelvin to (Hartree) units",T)
            function moment(t,n)
                logE=log(theta)+t
                logE >= log(floatmax(Float64)) && return 0.
                E=exp(logE);f=NMI.electronEnergyDistribution(E,T;distribution=model,kappa=kap)
                f==0. && return 0.
                exp(log(f)+(n+1)*logE)
            end
            mass=first(quadgk(t->moment(t,0),-Inf,0.,Inf;rtol=1.e-9))
            mean=first(quadgk(t->moment(t,1),-Inf,0.,Inf;rtol=1.e-9))
            @test mass ≈ 1. rtol=2.e-7
            @test mean ≈ 1.5theta rtol=2.e-7
            @test NMI.electronEnergyDistribution(0.,T;distribution=model,kappa=kap)==0.
        end
    end
    @testset "Analytic constant-strength rates" begin
        prefactor=8.6291269e-6/sqrt(nmT)
        for kap in [1.6,2.,3.,10.,1.e4]
            lines=[nmline(nmlevel(2,.5),E,2.) for E in [0.,1.,1000.]]
            r=nmrate(lines;kap=kap)
            # Independent high-precision gamma ratio, not the PDF's log-beta implementation.
            factor=setprecision(128) do
                k=BigFloat(kap);a=k-BigFloat(1.5)
                Float64(gamma(k+1)/gamma(k-BigFloat(.5))/sqrt(a)/k*(1+BigFloat(.5)/a)^(-k))
            end
            @test r.alphas[1] ≈ prefactor*2factor rtol=2.e-8
            @test r.distribution==:kappa && r.kappa==kap
            @test !hasproperty(r,:effOmegas)
        end
        b=1.5gamma(.75)/gamma(1.25)
        for delta in [.01,.5,5.,100.]
            lines=[nmline(nmlevel(2,delta),E,2.) for E in [0.,1.,1000.]]
            r=nmrate(lines;model=:Druyvesteyn)
            expected=prefactor*2*pi/(2sqrt(b)*gamma(.75))*erfc(delta/b)
            @test r.alphas[1] ≈ expected rtol=2.e-8
            @test isnothing(r.kappa) && r.distribution==:druyvesteyn
        end
    end
    @testset "Maxwellian limit, shift and arbitrary grouping" begin
        lines=[nmline(nmlevel(2,.5),E,2+.1*E) for E in [0.,1.,10.,100.]]
        m=nmrate(lines;model=:Maxwellian)
        k=nmrate(lines;kap=1.e6)
        @test k.alphas[1] ≈ m.alphas[1] rtol=1.e-5
        shifted=[nmline(nmlevel(3,.8),E,2+.1*E;shift=.3) for E in [0.,1.,10.,100.]]
        for model in [:Kappa,:Druyvesteyn]
            @test nmrate(shifted;model=model).alphas[1] ≈ nmrate(lines;model=model).alphas[1] rtol=1.e-10
            shuffled=deepcopy(vcat(lines[[3,1]],shifted[[4,1,2]],lines[[4,2]]))
            r=with_logger(NullLogger()) do
                NMI.computeRateCoefficients(shuffled,nmsettings();distribution=model)
            end
            @test length(r)==2
            @test r[1].alphas[1] ≈ r[2].alphas[1] rtol=1.e-10
        end
    end
    @testset "Extrapolation diagnostics and labelled output" begin
        lines=[nmline(nmlevel(2,.5),E,2.) for E in [.1,.2]]
        @test_logs (:warn,r"extrapolated collision strengths") NMI.computeRateCoefficients(lines,nmsettings();distribution=:Kappa,kappa=3.)
        r=nmrate(lines;kap=3.);a=1.5
        @test r.lowerExtrapolationFractions[1] ≈ 1-((a+.5+.1)/(a+.5))^(-3) rtol=1.e-8
        @test r.upperExtrapolationFractions[1] ≈ ((a+.5+.2)/(a+.5))^(-3) rtol=1.e-8
        io=IOBuffer();NMI.displayResults(io,[r]);text=String(take!(io))
        @test occursin("kappa",text) && occursin("3kBT/2",text)
        @test !occursin("Eff Coll Strength",text)
        for model in [:Kappa,:Druyvesteyn]
            z=nmrate([nmline(nmlevel(2,.5),E,0.) for E in [0.,1.]];model=model)
            @test z.alphas==[0.] && z.lowerExtrapolationFractions==[0.] && z.upperExtrapolationFractions==[0.]
        end
    end
    @testset "Validation" begin
        for kap in [1.5,1.,Inf,NaN]
            @test_throws ErrorException NMI.electronEnergyDistribution(.1,nmT;distribution=:Kappa,kappa=kap)
            @test_throws ErrorException NMI.computeRateCoefficients(NMI.Line[],nmsettings();distribution=:Kappa,kappa=kap)
        end
        @test_throws ErrorException NMI.electronEnergyDistribution(-1.,nmT;distribution=:Druyvesteyn)
        @test_throws ErrorException NMI.electronEnergyDistribution(.1,0.;distribution=:Druyvesteyn)
        @test_throws ErrorException NMI.computeRateCoefficients(NMI.Line[],nmsettings();distribution=:unknown)
        @test_throws ErrorException NMI.computeRateCoefficients([nmline(nmlevel(2,.5),0.,1.)],nmsettings();distribution=:Kappa)
    end
    @testset "Dipole extrapolation with independent energy-space integration" begin
        f=Level(AngularJ64(1),nmempty.M,Basics.minus,2,.5,nmempty.relativeOcc,nmempty.hasStateRep,nmempty.basis,nmempty.mc)
        xs=[0.,.1,.2];ys=[1.,1.5,2.]
        lines=[nmline(f,x,y) for (x,y) in zip(xs,ys)]
        for model in [:Kappa,:Druyvesteyn]
            r=nmrate(lines;model=model)
            function reference(u)
                omega=NMI.interpolateCS(u-.5,xs,ys,true)
                weight = if model==:Kappa
                    gamma(4.)/gamma(2.5)/1.5^1.5*(1+u/1.5)^(-4)
                else
                    b=1.5gamma(.75)/gamma(1.25)
                    sqrt(pi)/(b^1.5*gamma(.75))*exp(-(u/b)^2)
                end
                omega*weight
            end
            direct=first(quadgk(reference,.5,.6,.7,Inf;rtol=1.e-10))*8.6291269e-6/sqrt(nmT)
            @test r.alphas[1] ≈ direct rtol=1.e-8
            @test r.upperExtrapolationFractions[1] > .01
        end
    end
    @testset "Independent SI cross-section integrals" begin
        a0=5.29177210903e-11;Eh=4.3597447222071e-18;me=9.1093837015e-31;kB=1.380649e-23
        T=5.e4;delta=.5*Eh
        lines=[nmline(nmlevel(2,.5),E,2.) for E in [0.,1.,1000.]]
        for model in [:Kappa,:Druyvesteyn]
            function reference(u)
                E=u*kB*T
                f = if model==:Kappa
                    a=1.5*kB*T
                    2/sqrt(pi)*gamma(4.)/gamma(2.5)/a^1.5*sqrt(E)*(1+E/a)^(-4)
                else
                    b=1.5*kB*T*gamma(.75)/gamma(1.25)
                    2/(b^1.5*gamma(.75))*sqrt(E)*exp(-(E/b)^2)
                end
                sigma=pi*a0^2*2/(2*E/Eh)
                sigma*sqrt(2*E/me)*f*kB*T*1.e6
            end
            direct=first(quadgk(reference,delta/(kB*T),delta/(kB*T)+1,Inf;rtol=1.e-10))
            @test nmrate(lines;model=model,T=T).alphas[1] ≈ direct rtol=2.e-5
        end
    end
end
