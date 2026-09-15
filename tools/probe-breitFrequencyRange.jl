#
# probe-breitFrequencyRange.jl   --   WHERE THE FREQUENCY-DEPENDENT BREIT STRENGTH STOPS BEING TRUSTWORTHY.
#
# PRIORITY ITEM 19.  `CoulombBreit(factor)` sets omega = factor |E_a - E_c| / c from the ORBITAL ENERGIES, and
# `InteractionStrength.checkFrequencyIsMeaningful` guards only the omega -> 0 end (four orbitals at exactly zero,
# an unset field).  NOTHING GUARDS THE OTHER END.  In an Auger or dielectronic-capture amplitude one orbital is a
# CONTINUUM electron, so |E_a - E_c| is the continuum energy itself -- 727 a.u. at Z = 53 rising to 2338 at Z = 92,
# i.e. omega = 5.3 to 17 a.u. -- and the amplitudes measured on 09-Sep-2026 degrade from a healthy 0.2 % correction
# at Z = 53 to values wrong by 10^6-10^7 at Z = 92, with |b/a| jumping 0.25, 0.048, 39.9, 3.2, 6.1, 6.2 from
# structure input that varies perfectly smoothly.  That is the signature of an oscillatory integrand the quadrature
# cannot follow, not of a wrong kernel.
#
# THIS FILE LOCATES THE BOUNDARY DIRECTLY, and it does so on ONE orbital quadruple with omega as the only variable,
# which the Z scan could not:  there, Z changed the orbitals, the energies and the mixing all at once.  Here the
# orbitals are FIXED and only the energy difference is varied, so any departure belongs to omega alone.
#
# HOW TO READ IT.  The frequency-dependent strength must approach the omega -> 0 value smoothly and quadratically:
# the leading retardation correction is O(omega^2), which is exactly the law `example-Ad.jl` branch 4 validates for
# BOUND-BOUND transitions at small omega.  So the ratio XL(omega)/XL(0) should start at 1 and bend away gently.
# The breakdown is where it stops doing that -- where the ratio moves by orders of magnitude, or changes sign, for a
# smooth change in omega.  `omega*rbox` is printed beside it because that is the quantity a guard can test cheaply.
#
using JenaAtomicCalculator, Printf
const IS = JenaAtomicCalculator.InteractionStrength

# TWO SYSTEMS WITH VERY DIFFERENT BOXES, because that is what discriminates the criterion.  If the breakdown sits
# at the same omega in both, the guard must test OMEGA;  if it sits at the same omega*rbox, it must test that
# product instead.  The item proposed omega*rbox, but the Auger measurements of 09-Sep were HEALTHY at
# omega*rbox = 21 (Z = 53, a 4 a.u. box) while the first run of this probe already turned over at 8.8 in a 0.8 a.u.
# box -- so the two cannot both be on the same curve in that variable.
for (Z, conf) in [(54.0, [Configuration("1s^2 2s^2 2p^6")]), (6.0, [Configuration("1s^2 2s^2 2p^2")])]
Z    = Z
conf = conf
nm   = Nuclear.Model(Z)
grid = Basics.recommendedGrid(conf, nm; printout=false)
set  = AsfSettings(AsfSettings(); scField = Basics.DFSField(), eeInteraction = CoulombInteraction(),
                                  eeInteractionCI = CoulombInteraction(), gridStopper = false)
mp   = redirect_stdout(devnull) do
           SelfConsistent.performSCF(conf, nm, grid, set; printout=false)
       end
orbs = mp.levels[1].basis.orbitals
rbox = grid.r[end]
# a fine scan through the collapse region, expressed in omega*rbox so both systems are sampled at the same products
dEs = sort(unique(vcat([0.0], [ p * 137.035999 / grid.r[end] for p in
            [0.5,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,18.0,20.0] ])))
a0 = orbs[Subshell("2p_1/2")];   b0 = orbs[Subshell("2s_1/2")]
c0 = orbs[Subshell("2p_3/2")];   d0 = orbs[Subshell("2s_1/2")]

"the same orbital with its energy replaced, so omega is the ONLY thing that varies"
withEnergy(o::Orbital, e::Float64) = Orbital(o.subshell, o.isBound, o.useStandardGrid, e, o.P, o.Q, o.Pprime, o.Qprime, o.grid)

ref = IS.XL_Breit(1, a0, b0, c0, d0, grid, CoulombBreit(0.))
@printf("\nZ = %.0f,  rbox = %.1f a.u.,  XL_Breit at omega = 0:  %.6e\n\n", Z, rbox, ref)
@printf("%10s %12s %12s %16s %14s\n", "dE [a.u.]", "omega", "omega*rbox", "XL_Breit", "ratio to omega=0")
println("-"^70)
for dE in dEs
    a = withEnergy(a0, dE);    c = withEnergy(c0, 0.0)
    v = IS.XL_Breit(1, a, b0, c, d0, grid, CoulombBreit(1.))
    om = dE / 137.035999
    @printf("%10.1f %12.4f %12.2f %16.6e %14.4f\n", dE, om, om*rbox, v, v/ref)
    flush(stdout)
end
end
println("\nTHE RATIO SHOULD LEAVE 1 SMOOTHLY AND QUADRATICALLY (the O(omega^2) retardation law).  Where it instead")
println("moves by orders of magnitude, or changes sign, for a smooth change in omega, the quadrature has lost the")
println("integrand -- and `omega*rbox` at that point is what a guard should test.")
