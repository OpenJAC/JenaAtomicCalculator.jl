# Bounded correctness checks, including a phase-sensitive comparison of both amplitude implementations.
# Run from any writable working directory with --project=<JAC-root> -t 2.
using JenaAtomicCalculator, Test, FastGaussQuadrature, Printf
using JenaAtomicCalculator: ImpactExcitation, Continuum, RadialIntegrals
const IE = ImpactExcitation

@testset "Consecutive, weighted shell convergence" begin
    @test IE.updateShellConvergence(1e-7, 1., 4) == (1e-7, 5)
    @test IE.updateShellConvergence(1e-3, 1., 5) == (1e-3, 0)
    @test IE.updateShellConvergence(0., 0., 5) == (Inf, 0)
    @test_throws ErrorException IE.updateShellConvergence(NaN, 1., 0)
    @test_throws ErrorException IE.updateShellConvergence(-1., 1., 0)
    @test_throws ErrorException IE.collisionStrengthNormalization(1., 0.)
end

setDefaults("unit: energy", "eV")
setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, Alok")
grid = Radial.Grid(Radial.Grid(true), rnt=4e-6, h=0.05, hp=0.008, rbox=30.)
@testset "Relativistic grid-resolution guard" begin
    # At this energy 15hp lies between the nonrelativistic and relativistic wavelengths.
    fastgrid=Radial.Grid(grid;hp=0.0028,rbox=10.)
    @test_throws ErrorException Continuum.gridConsistency(10000.,fastgrid)
end
@testset "Projection with nonorthogonal bound functions and finite continuum extent" begin
    r=grid.r; z=zeros(length(r)); s1=Subshell(1,-1); s2=Subshell(2,-1)
    b1=Orbital(s1,true,true,-1.,r.*exp.(-r),z,z,z,grid)
    b2=Orbital(s2,true,true,-.5,r.^2 .* exp.(-r),z,z,z,grid)
    basis=Basis(true,2,[s1,s2],CsfR[],Subshell[],Dict(s1=>b1,s2=>b2))
    l=Level()
    target=Level(l.J,l.M,l.parity,l.index,l.energy,l.relativeOcc,l.hasStateRep,basis,Float64[])
    n=findlast(x->x<=3.,r); zz=zeros(n)
    wave=Orbital(Subshell(101,-1),false,true,1.,sin.(r[1:n]),zz,zz,zz,grid)
    projected=Continuum.orthogonalizeToBasis(wave,target,grid)
    @test abs(RadialIntegrals.overlap(b1,projected,grid)) < 1e-12
    @test abs(RadialIntegrals.overlap(b2,projected,grid)) < 1e-12
    @test Continuum.orthogonalizeToBasis(projected,target,grid).P ≈ projected.P atol=1e-12
end
@testset "Coulomb tail validity" begin
    E = 1.; c = Defaults.getDefaults("speed of light: c")
    q = sqrt(E*(E+2c^2))/c
    N = sqrt((E+2c^2)/(pi*c^2*q))
    P = N .* sin.(q .* grid.r); z = zeros(length(P))
    orb = Orbital(Subshell(101,-1), false, true, E, P, z, N*q .* cos.(q .* grid.r), z, grid)
    tail = IE.continuumTail(orb, grid, 0.)
    @test tail.A ≈ N rtol=1e-9
    @test abs(tail.B) < 1e-10
    @test tail.fitError < 1e-9
    # A zero wave in an allowed region used to acquire a fabricated, perfectly fitted Coulomb tail.
    bad = Orbital(Subshell(101,-1), false, true, E, z, z, z, z, grid)
    @test_throws ErrorException IE.continuumTail(bad, grid, 0.)
    @test all(isfinite, IE.coulombFG(tail, 2tail.r0))
    regular,_ = IE.coulombContinuumOrbital(E,Subshell(101,100),grid,0.,grid.NoPoints)
    analytic = IE.continuumTail(regular,grid,0.;regularCoulomb=true)
    @test analytic.A ≈ N
    @test analytic.B == 0. && isnan(analytic.fitError)
    @test all(isfinite,IE.coulombFG(analytic,200.))
end

x, weights = gausslegendre(32)
angles = rad2deg.(acos.(x))
function calculation(method; topup=false)
    settings = IE.Settings(IE.Settings(); electronEnergies=[50.], maxKappa=6,
        calcRateCoefficient=false, continuumExchange=:none, amplitudeMethod=method,
        calcDCS=true, angles=angles, asymptoticTail=true, topUp=topup, printBefore=false)
    comp = Atomic.Computation(Atomic.Computation(); name="IE reliability: C3+",
        grid=grid, nuclearModel=Nuclear.Model(6.),
        initialConfigs=[Configuration("1s^2 2s")], finalConfigs=[Configuration("1s^2 2p")],
        processSettings=settings)
    redirect_stdout(devnull) do
        perform(comp; output=true)[IE.Settings]
    end
end

@testset "Complex amplitudes, DCS and its angular integral" begin
    a = calculation(:csfTable); b = calculation(:density)
    @test length(a) == length(b) == 2
    for (table, density) in zip(a,b)
        @test [(c.initialKappa,c.finalKappa,c.symmetry) for c in table.channels] ==
              [(c.initialKappa,c.finalKappa,c.symmetry) for c in density.channels]
        @test maximum(max(abs(c.initialKappa),abs(c.finalKappa)) for c in density.channels) == 6
        ta = [c.amplitude for c in table.channels]; da = [c.amplitude for c in density.channels]
        @test maximum(abs.(ta-da)) < 1e-10 * maximum(abs.(ta))
        @test table.crossSection ≈ density.crossSection rtol=1e-10
        @test table.dcs ≈ density.dcs rtol=1e-10
        @test IE.computeDCS(density, angles) ≈ density.dcs rtol=1e-12
        integral = 2pi * sum(weights .* density.dcs)
        @printf("Jf=%s: sigma=%.10g, integral(DCS)/sigma=%.10g, max relative amplitude difference=%.3e\n",
            string(density.finalLevel.J), density.crossSection, integral/density.crossSection,
            maximum(abs.(ta-da))/maximum(abs.(ta)))
        @test integral ≈ density.crossSection rtol=1e-11
    end
end
