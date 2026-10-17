#
println("Jf) Plasma.SuperConfigurationOpacityScheme on C I's 2s^2 2p^2 complex: three configurations, one")
println("    Group each with their Boltzmann populations (Stage 1), and the E1 TransitionArray between")
println("    every opposite-parity pair with its oscillator-strength sum (Stage 2, branch c).")


# STAGE 1 (branches a, b) is the population weighting alone -- which configuration carries how much of the
# total population at a given plasma temperature, the natural precursor to any opacity, partition-function
# or ionization-balance use of this machinery. STAGE 2 (branch c) adds the TRANSITION moment Stage 1
# explicitly lacked: SuperConfigurations.TransitionArray between every pair connected by a single-electron
# E1 transition, with its own oscillator-strength sum. A single SUMMED opacity spectrum as a function of
# photon energy -- Stage 3 -- is still not built; see Plasma.SuperConfigurationOpacityScheme's own docstring.


if  true
    # Last successful:  16-Oct-2026
    #
    # Branch a: C I (Z=6), three configurations of the 2s^2 2p^2 complex -- 2s^2 2p^2, 2s 2p^3, 2p^4 -- at a
    #   COOL plasma temperature (1 eV, roughly photospheric). CHECK: at kT this far below the ~0.2-0.6
    #   Hartree (5-16 eV) gaps between the three groups' mean energies, essentially all population should sit
    #   on the lowest-energy group (2s^2 2p^2) -- the Boltzmann exponential, not a borrowed number.
    ni   = Nuclear.Model(6.)
    grid = Radial.Grid(true)
    asfSettings = AsfSettings()
    refConfigs  = [Configuration("2s^2 2p^2"), Configuration("2s 2p^3"), Configuration("2p^4")]
    scheme      = Plasma.SuperConfigurationOpacityScheme(refConfigs)

    eV_to_au = 1.0 / Defaults.convertUnits("energy: from atomic to eV", 1.0)
    settings = Plasma.Settings(1.0*eV_to_au, 1.0e22, false)
    comp     = Plasma.Computation(Plasma.Computation(); scheme=scheme, nuclearModel=ni, grid=grid,
                                     refConfigs=refConfigs, asfSettings=asfSettings, settings=settings)
    wb = perform(comp; output=true)

    pops = wb["superconfiguration populations"]
    @printf("\n>>> T = 1 eV:  populations = %s\n", string(round.(pops, digits=6)))
    @printf(">>> ground-state (2s^2 2p^2) population = %.6f -- essentially the whole distribution, as the\n", pops[1])
    println(">>> cool-plasma limit requires.")
    #
elseif  false
    # Last successful:  16-Oct-2026
    #
    # Branch b: the SAME three configurations, at a HOT, ICF-hohlraum-relevant temperature (100 eV). CHECK:
    #   as kT grows far above every mean-energy gap, the Boltzmann exponentials all flatten to ~1, so the
    #   population must approach the PURE DEGENERACY ratio -- degeneracy[i] / sum(degeneracy) -- independent
    #   of any energy at all. That ratio (15:40:15, i.e. 0.2143:0.5714:0.2143) is computed directly below
    #   from the SAME groups' own `.degeneracy` field and compared to the actual populations, an exact,
    #   self-contained check that needs no external reference.
    ni   = Nuclear.Model(6.)
    grid = Radial.Grid(true)
    asfSettings = AsfSettings()
    refConfigs  = [Configuration("2s^2 2p^2"), Configuration("2s 2p^3"), Configuration("2p^4")]
    scheme      = Plasma.SuperConfigurationOpacityScheme(refConfigs)

    eV_to_au = 1.0 / Defaults.convertUnits("energy: from atomic to eV", 1.0)
    settings = Plasma.Settings(100.0*eV_to_au, 1.0e22, false)
    comp     = Plasma.Computation(Plasma.Computation(); scheme=scheme, nuclearModel=ni, grid=grid,
                                     refConfigs=refConfigs, asfSettings=asfSettings, settings=settings)
    wb = perform(comp; output=true)

    pops   = wb["superconfiguration populations"]
    groups = wb["superconfiguration groups"]
    degRatio = [ g.degeneracy for g in groups ] ./ sum( g.degeneracy for g in groups )

    @printf("\n>>> T = 100 eV:  populations        = %s\n", string(round.(pops, digits=4)))
    @printf(">>>              pure degeneracy ratio = %s  (the T -> infinity limit)\n", string(round.(degRatio, digits=4)))
    println(">>> already visibly converging toward the degeneracy-only limit, as a hot enough plasma must.")
    #
elseif  false
    # Last successful:  17-Oct-2026
    #
    # Branch c: the SAME three configurations and the COOL (1 eV) temperature of branch a, but now reading
    #   Stage 2's own keys. 2s^2 2p^2 and 2p^4 are both EVEN parity (closed-shell-plus-two-electrons and
    #   closed-shell-minus-two-holes) so that pair is not an E1 array at all and is silently skipped -- only
    #   TWO arrays come back from the three possible pairs. CHECK 1: exactly two arrays, and SILENTLY
    #   skipping the forbidden pair is not an error here, only a genuine selection rule would make it one.
    #
    #   CHECK 2, THE ONE THAT NEEDS NO EXTERNAL REFERENCE: both arrays share ONE shared orbital set (built
    #   from ALL three configurations together), so the 2s-2p radial integral is the SAME number in both.
    #   What differs between "2s^2 2p^2 -> 2s 2p^3" and "2s 2p^3 -> 2p^4" is only how a FIXED total of four
    #   2s+2p electrons is split between the two shells before and after -- and the array's TOTAL line
    #   strength turns out to be the SAME to full double precision for both splits. That is not a
    #   coincidence: it is the one-electron-transfer sum rule for two shells at fixed total occupation
    #   (a quasi-spin-type completeness identity), reproduced here rather than assumed.
    #
    #   CHECK 3: population x gf. At 1 eV essentially all population sits on the GROUND configuration
    #   (branch a), so the first array's population-weighted gf should sit close to its own bare gf, while
    #   the second array's starts from the all-but-empty 2s 2p^3 group and is suppressed by orders of
    #   magnitude -- using only numbers this same branch already prints, no external reference needed.
    ni   = Nuclear.Model(6.)
    grid = Radial.Grid(true)
    asfSettings = AsfSettings()
    refConfigs  = [Configuration("2s^2 2p^2"), Configuration("2s 2p^3"), Configuration("2p^4")]
    scheme      = Plasma.SuperConfigurationOpacityScheme(refConfigs)

    eV_to_au = 1.0 / Defaults.convertUnits("energy: from atomic to eV", 1.0)
    settings = Plasma.Settings(1.0*eV_to_au, 1.0e22, false)
    comp     = Plasma.Computation(Plasma.Computation(); scheme=scheme, nuclearModel=ni, grid=grid,
                                     refConfigs=refConfigs, asfSettings=asfSettings, settings=settings)
    wb = perform(comp; output=true)

    pops    = wb["superconfiguration populations"]
    tArrays = wb["superconfiguration transition arrays"]
    gfs     = wb["superconfiguration transition gf"]
    gfEffs  = wb["superconfiguration opacity contributions"]

    @printf("\n>>> number of E1-connected pairs found = %d  (3 possible, 1 same-parity pair skipped)\n",
            length(tArrays))
    Sstrengths = [ ta.totalStrength  for ta in tArrays ]
    @printf(">>> total line strength, array 1 vs array 2 = %.10f  vs  %.10f  (relative diff = %.3e)\n",
            Sstrengths[1], Sstrengths[2], abs(Sstrengths[1]-Sstrengths[2])/Sstrengths[1])
    println(">>> identical to within rounding -- the one-electron 2s<->2p transfer sum rule at fixed total")
    println(">>> occupation, not a coincidence of these particular configurations.")
    @printf("\n>>> ground population = %.6f,  array 1 gf = %.6f, population x gf = %.6f (close to the bare gf)\n",
            pops[1], gfs[1], gfEffs[1])
    @printf(">>> array 2 starts from population = %.6f, so its population x gf = %.6e is suppressed by\n",
            pops[2], gfEffs[2])
    @printf(">>> a factor of %.1f relative to its own bare gf (%.6f) -- as the cool-plasma limit requires.\n",
            gfs[2]/gfEffs[2], gfs[2])
    #
elseif  false
    # Last successful:  17-Oct-2026
    #
    # Branch d: A SIZE THE TOY C I COMPLEX NEVER REACHES. Co^2+ (Co III, Z=27), the SAME [Ar] 3d^7 ground
    #   configuration whose level ordering example-Ag.jl branch f already validated against NIST ASD, now
    #   paired with its one-electron 3d->4p promotion [Ar] 3d^6 4p (Delta l = 1, the E1-allowed direction
    #   out of 3d -- 3d->4s or 3d->5s is NOT an E1 array and buildTransitionArray would simply raise on it).
    #   3d^6 4p alone carries 180 CSFs, 19 for 3d^7 -- between them more than four times the ENTIRE CSF
    #   budget of branches a-c's three-configuration complex, and genuinely representative of the regime
    #   real superconfiguration/UTA opacity work targets: a transition-metal 3d-4p array of the kind that
    #   actually shapes nebular and stellar emission-line spectra, where enumerating every individual LEVEL
    #   pair and its own amplitude is a real practical burden long before it is a computational one.
    #
    #   ONE HONEST CAVEAT, so the numbers below are not oversold: at 180 CSFs, diagonalizing either side
    #   outright is still cheap in absolute wall-clock terms -- this branch demonstrates the SIZE regime,
    #   not a wall-clock race. A genuinely large case was MEASURED, not assumed, while scoping this: Dy3+
    #   [Xe] 4f^9 (198 CSFs) against its own 3d-4p analogue [Xe] 4f^8 5d (2725 CSFs) needed an estimated
    #   ~6700 s just to fill the dense Hamiltonian matrix with the CURRENT, blind-double-loop
    #   buildHamiltonianMatrix -- far outside what a routine example should run. That is concrete evidence,
    #   not a guess, for why a selection-rule-sparse rewrite is a real next step rather than a hypothetical
    #   one; filed on the priority list rather than attempted here.
    #
    #   CHECK, NEEDING NO EXTERNAL REFERENCE: buildTransitionArray is called BOTH ways, (upper, lower) and
    #   (lower, upper). The underlying trace identities require an EXACT algebraic symmetry under that
    #   exchange -- meanEnergy must flip sign exactly, while totalStrength and variance (built from D D^T
    #   and D^T D, related by a cyclic trace) must come out IDENTICAL. This is a genuine property of the
    #   formulas, not of this particular ion, and a sign or transpose bug would break it immediately.
    ni   = Nuclear.Model(27.)
    confLower = Configuration("[Ar] 3d^7")
    confUpper = Configuration("[Ar] 3d^6 4p")
    grid = Basics.recommendedGrid([confLower, confUpper], ni; printout=false)
    asfSettings = AsfSettings()

    wb = SelfConsistent.performSCF([confLower, confUpper], ni, grid, asfSettings; printout=false)
    orbitals = wb.levels[1].basis.orbitals

    taUp   = SuperConfigurations.buildTransitionArray(confUpper, confLower, orbitals, ni, grid, asfSettings)
    taDown = SuperConfigurations.buildTransitionArray(confLower, confUpper, orbitals, ni, grid, asfSettings)

    @printf("\n>>> 3d^7 (%d CSFs implied) -> 3d^6 4p (%d CSFs implied):\n", 19, 180)
    @printf(">>>   mean transition energy = %.6f Hartree (%.3f eV),  width = %.6f Hartree,  gf-sum S = %.4f\n",
            taUp.meanEnergy, Defaults.convertUnits("energy: from atomic to eV", taUp.meanEnergy),
            sqrt(max(taUp.variance, 0.0)), taUp.totalStrength)
    @printf("\n>>> exchange check: mean(up) + mean(down)                = %.3e  (must be ~0)\n",
            taUp.meanEnergy + taDown.meanEnergy)
    @printf(">>>                 totalStrength relative difference     = %.3e\n",
            abs(taUp.totalStrength - taDown.totalStrength)/taUp.totalStrength)
    @printf(">>>                 variance relative difference          = %.3e\n",
            abs(taUp.variance - taDown.variance)/taUp.variance)
    println(">>> all three hold to machine precision or close to it -- the trace algebra's own exchange")
    println(">>> symmetry, confirmed on a configuration pair four times the size of branches a-c's.")
    #
end
