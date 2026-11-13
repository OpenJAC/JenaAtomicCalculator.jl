# Standalone: julia --project=. test/impact-excitation-rates.jl
using Test, JenaAtomicCalculator, QuadGK, Logging
using JenaAtomicCalculator: ImpactExcitation
const IE = ImpactExcitation
setDefaults("unit: energy", "Hartree")
const emptyLevel = Level()
makeLevel(index, E, J=0, parity=Basics.plus) = Level(AngularJ64(J), emptyLevel.M, parity, index, E,
    emptyLevel.relativeOcc, emptyLevel.hasStateRep, emptyLevel.basis, emptyLevel.mc)
const initial = makeLevel(1, 0.)
const final = makeLevel(2, .5)
makeLine(f, E, omega; shift=0., initialState=initial) = IE.Line(initialState, f,
    E+f.energy-initialState.energy-shift, E, 0., omega, Float64[], Float64[], IE.Channel[], 0.)
const temperature = Defaults.convertUnits("temperature: from atomic to Kelvin", 1.)
settings(n; temps=[temperature]) = IE.Settings(IE.Settings(); numElectronEnergies=n, temperatures=temps)

@testset "ImpactExcitation thermal rates" begin
    @testset "Analytic integrals" begin
        for kT in [1.e-5, .1, 1., 100., 1.e8]
            @test IE.maxwellianStrength([0., 1., 10.], [2., 2., 2.], kT, false) ≈ 2. rtol=2.e-13
            # Linear table with a constant continuation after L: Upsilon=2+3kT(1-exp(-L/kT)).
            @test IE.maxwellianStrength([0., 1., 10.], [2., 5., 32.], kT, false) ≈
                2+3kT*(-expm1(-10/kT)) rtol=2.e-12
        end
        energies = collect(range(0., 10.; length=2001))
        triangle(E) = max(0., 1-abs(E-.75)/.25)
        lines = [makeLine(final, E, triangle(E)) for E in energies]
        exact = (exp(-.5)-2exp(-.75)+exp(-1.))/.25
        @test only(IE.computeEffStrengths(lines, settings(length(lines)))).effOmegas[1] ≈ exact rtol=2.e-12
        @test IE.maxwellianStrength([0., .5, .75, 1., 10.], [0., 0., 1., 0., 0.], 1., false) ≈ exact rtol=2.e-12
    end
    @testset "Kinematic shift and generated grid" begin
        for shift in [0., .2, -.1]
            lines = [makeLine(final, E, 2.; shift=shift) for E in [0., 1., 100.]]
            r = only(IE.computeRateCoefficients(lines, settings(999)))
            expected = 8.6291269e-6*2exp(-(.5-shift))/sqrt(temperature)
            @test r.alphas[1] ≈ expected rtol=1.e-13
            s = IE.Settings(settings(10); calcRateCoefficient=true, electronEnergies=Float64[], energyShift=shift)
            generated = IE.determineLines(Multiplet("f", [final]), Multiplet("i", [initial]), s)
            @test minimum(l.finalElectronEnergy for l in generated) ≈ .0003 atol=1.e-14
            @test IE.kineticThreshold(generated) ≈ .5-shift rtol=1.e-13
        end
        s = IE.Settings(settings(10); calcRateCoefficient=true, electronEnergies=Float64[], energyShift=.6)
        @test_throws ErrorException IE.determineLines(Multiplet("f", [final]), Multiplet("i", [initial]), s)
        ji2 = makeLevel(4, 0., 2)
        r = only(IE.computeRateCoefficients([makeLine(final,E,2.;initialState=ji2) for E in [0.,100.]],settings(2)))
        @test r.alphas[1] ≈ 8.6291269e-6*2exp(-.5)/(5sqrt(temperature)) rtol=1.e-13
        crossSettings = IE.Settings(settings(1);calcRateCoefficient=false,electronEnergies=[1.],energyShift=.2)
        crossLines = IE.determineLines(Multiplet("f",[final]),Multiplet("i",[initial]),crossSettings)
        @test only(crossLines).initialElectronEnergy == 1.
        @test only(crossLines).finalElectronEnergy ≈ .7
    end
    @testset "Actual transition grouping" begin
        other = makeLevel(3, .8)
        lines = [makeLine(final, 2., 2.), makeLine(other, 1., 3.), makeLine(final, 0., 2.),
                 makeLine(other, 0., 3.), makeLine(final, 1., 2.)]
        groups = IE.groupLines(lines, settings(999))
        @test length.(groups) == [3, 2]
        @test all(g -> length(unique(l.finalLevel.index for l in g)) == 1, groups)
        @test all(g -> issorted([l.finalElectronEnergy for l in g]), groups)
        rates = with_logger(NullLogger()) do
            IE.computeRateCoefficients(lines, settings(999))
        end
        @test length(rates) == 2
        @test rates[1].effOmegas[1] ≈ 2.
        @test rates[2].effOmegas[1] ≈ 3.
        copied = [deepcopy(makeLine(final,E,2.)) for E in [2.,0.,1.]]
        @test length(IE.groupLines(copied,settings(99))) == 1
        @test only(IE.computeEffStrengths(copied,settings(99))).effOmegas[1] ≈ 2.
    end
    @testset "Finite, positive anchored E1 tails" begin
        for endpoint in [.5, 1., 2.]
            xs = [0., endpoint/2, endpoint]; ys = [0.5, 1., 2.]
            @test IE.interpolateCS(endpoint, xs, ys, true) == 2.
            @test isfinite(IE.interpolateCS(10endpoint, xs, ys, true))
            @test IE.interpolateCS(10endpoint, xs, ys, true) >= 2.
            # Independently integrate each declared interval and extrapolation with adaptive quadrature.
            ref = sum(first(quadgk(t -> IE.interpolateCS(t, xs, ys, true)*exp(-t), a, b; rtol=1.e-11))
                      for (a,b) in [(0., endpoint/2), (endpoint/2, endpoint), (endpoint, Inf)])
            @test IE.maxwellianStrength(xs, ys, 1., true) ≈ ref rtol=2.e-10
        end
        @test IE.interpolateCS(2., [.1, .5], [2., 1.], true) == 1.
        @test IE.interpolateCS(2., [.1, 1.], [0., 0.], true) == 0.
    end
    @testset "Validation and coverage" begin
        @test_throws ErrorException IE.maxwellianStrength([0., 0.], [1., 1.], 1., false)
        @test_throws ErrorException IE.maxwellianStrength([0., 1.], [1., -1.], 1., false)
        @test_throws ErrorException IE.maxwellianStrength([0., 1.], [1., NaN], 1., false)
        @test_throws ErrorException IE.maxwellianStrength([0., 1.], [1., 1.], 0., false)
        @test_throws ErrorException IE.interpolateCS(2.,[0.,1.],[1.,1.])
        @test_throws ErrorException IE.interpolateCS(.5,[1.,0.],[1.,1.])
        @test_throws ErrorException IE.computeEffStrengths([makeLine(final,0.,1.),makeLine(final,0.,1.)], settings(2))
        @test_throws ErrorException IE.computeEffStrengths([makeLine(final,0.,1.),makeLine(final,1.,1.)], settings(2;temps=[0.]))
        @test_throws ErrorException IE.computeEffStrengths([makeLine(final,0.,1.),makeLine(final,1.,1.;shift=.1)], settings(2))
        @test_logs (:warn,r"outside the sampled range") IE.computeEffStrengths([makeLine(final,.1,1.),makeLine(final,1.,1.)], settings(2))
    end
    @testset "Independent SI Maxwellian cross-section integral" begin
        # Modern SI constants are independent of JAC's rounded Kelvin/Hartree conversion and rate prefactor.
        a0 = 5.29177210903e-11; Eh = 4.3597447222071e-18
        me = 9.1093837015e-31; kB = 1.380649e-23
        for T in [5.e4, 1.e5, 1.e6]
            kT = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", T)
            lines = [makeLine(final,E,2+.1*E;shift=.2) for E in [0., 1., 1000.]]
            rate = only(IE.computeRateCoefficients(lines,settings(3;temps=[T]))).alphas[1]
            delta = .3*Eh
            function integrand(u)
                u > 745. && return 0.
                E = delta+kB*T*u
                omega = 2+.1*(E-delta)/Eh
                sigma = pi*a0^2*omega/(2*E/Eh) # gi=1; nonrelativistic incoming momentum.
                v = sqrt(2E/me)
                pdf = 2sqrt(E)/(sqrt(pi)*(kB*T)^1.5)*exp(-E/(kB*T))
                sigma*v*pdf*kB*T*1.e6
            end
            direct = first(quadgk(integrand,0.,Inf;rtol=1.e-10))
            @test rate ≈ direct rtol=1.e-5
        end
    end
end
