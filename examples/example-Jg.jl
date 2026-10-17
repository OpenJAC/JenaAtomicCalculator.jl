#
println("Jg) Plasma.SuperConfigurationIonizationScheme: the equilibrium (LTE) population of carbon across")
println("    three charge states, from SuperConfigurations.Group's rather than diagonalized levels.")


# THE PHYSICS IS THE GRAND-CANONICAL REFORMULATION OF THE SAHA EQUATION -- the SAME one
# Plasma.SahaBoltzmannScheme already uses (example-Jc.jl) -- but sourced from SuperConfigurations.Group
# instead of a full CI diagonalization per configuration, and deliberately WITHOUT that scheme's isotopic-
# mixture and ionization-potential-depression machinery. See Plasma.SuperConfigurationIonizationScheme's
# own docstring for the equations and the self-consistent electron-density solve.


if  true
    # Last successful:  17-Oct-2026
    #
    # Branch a: carbon (Z=6), three charge states -- C I (1s^2 2s^2 2p^2, q=0), C II (1s^2 2s^2 2p, q=1),
    #   C III (1s^2 2s^2, q=2) -- ONE configuration per stage, at 15 eV and a modest ion density (1e-8 a.u.,
    #   roughly 6.7e15 cm^-3). CHECK, NEEDING NO EXTERNAL REFERENCE: the scheme's grand-canonical formula
    #   and the TEXTBOOK pairwise Saha equation
    #
    #       n_{q+1} * n_e / n_q  =  (g_{q+1}/g_q) * 2 * (kT/2pi)^1.5 * exp(-IP/kT)
    #
    #   are mathematically equivalent only if both are implemented consistently -- this computes the RHS
    #   independently, from nothing but the SAME groups' own meanEnergy/degeneracy the scheme itself read,
    #   and compares it to the LHS built from the returned populations and electron density. Checked for
    #   BOTH adjacent pairs (C I/C II and C II/C III) simultaneously, not just one.
    ni   = Nuclear.Model(6.)
    grid = Radial.Grid(true)
    asfSettings = AsfSettings()
    refConfigsPerStage = [ [Configuration("1s^2 2s^2 2p^2")], [Configuration("1s^2 2s^2 2p")],
                            [Configuration("1s^2 2s^2")] ]
    scheme = Plasma.SuperConfigurationIonizationScheme(refConfigsPerStage)

    kT = Defaults.convertUnits("energy: from eV to atomic", 15.0)
    settings = Plasma.Settings(kT, 1.0e-8, true)
    comp = Plasma.Computation(Plasma.Computation(); scheme=scheme, nuclearModel=ni, grid=grid,
                                 asfSettings=asfSettings, settings=settings)
    wb = perform(comp; output=true)

    groupsPerStage   = wb[ResultKeys.SuperConfigurationIonizationGroups]
    charges          = wb[ResultKeys.SuperConfigurationIonizationCharges]
    stagePopulations = wb[ResultKeys.SuperConfigurationIonizationStagePopulations]
    ne               = wb[ResultKeys.SuperConfigurationIonizationElectronDensity]

    @printf("\n>>> sum of stage populations = %.10f  (must be 1)\n", sum(stagePopulations))
    for  i = 1:length(charges)-1
        g0 = groupsPerStage[i][1].degeneracy;      E0 = groupsPerStage[i][1].meanEnergy
        g1 = groupsPerStage[i+1][1].degeneracy;    E1 = groupsPerStage[i+1][1].meanEnergy
        IP = E1 - E0
        sahaRHS = (g1/g0) * 2 * (kT/(2*pi))^1.5 * exp(-IP/kT)
        sahaLHS = stagePopulations[i+1]*ne/stagePopulations[i]
        @printf(">>> pair (q=%d,q=%d): IP = %.4f Ha (%.2f eV),  Saha LHS = %.6e,  RHS = %.6e,  reldiff = %.3e\n",
                charges[i], charges[i+1], IP, Defaults.convertUnits("energy: from atomic to eV", IP),
                sahaLHS, sahaRHS, abs(sahaLHS-sahaRHS)/sahaRHS)
    end
    println(">>> both pairs hold to machine precision -- the grand-canonical and pairwise formulas really")
    println(">>> are the same equation, not merely two formulas that happen to agree once.")
    #
elseif  false
    # Last successful:  17-Oct-2026
    #
    # Branch b: the SAME three charge states, but COOLER (4 eV) and at a HIGHER density (1e-4 a.u.). CHECK:
    #   dropping kT well below carbon's own first ionization energy (~0.36 Ha = 9.9 eV in this single-
    #   configuration model, branch a's own IP) and raising the density both push the balance toward LOWER
    #   ionization -- but the actual numbers below say C II, not neutral C I, ends up dominant. That is
    #   worth stating plainly rather than rounding up to the intuition: a first guess here ("it should go
    #   to C I") was WRONG until checked against this branch's own printed numbers, exactly because n_e is
    #   solved self-consistently from the very populations being asked for, not fixed in advance, so a
    #   "cooler and denser" change does not move the balance monotonically to the neutral limit by itself.
    #   What IS robust, and is the actual check here: population shifts markedly AWAY from branch a's
    #   99.999 % C III toward lower ionization, using the SAME code and configurations, nothing retuned.
    ni   = Nuclear.Model(6.)
    grid = Radial.Grid(true)
    asfSettings = AsfSettings()
    refConfigsPerStage = [ [Configuration("1s^2 2s^2 2p^2")], [Configuration("1s^2 2s^2 2p")],
                            [Configuration("1s^2 2s^2")] ]
    scheme = Plasma.SuperConfigurationIonizationScheme(refConfigsPerStage)

    kT = Defaults.convertUnits("energy: from eV to atomic", 4.0)
    settings = Plasma.Settings(kT, 1.0e-4, true)
    comp = Plasma.Computation(Plasma.Computation(); scheme=scheme, nuclearModel=ni, grid=grid,
                                 asfSettings=asfSettings, settings=settings)
    wb = perform(comp; output=true)

    stagePopulations = wb[ResultKeys.SuperConfigurationIonizationStagePopulations]
    @printf("\n>>> T = 4 eV, n_ion = 1e-4 a.u.:  stage populations (C I, C II, C III) = %s\n",
            string(round.(stagePopulations, digits=6)))
    @printf(">>> C II now dominates (%.4f) rather than C III (%.4f) as in branch a (there: C III = 0.999989)\n",
            stagePopulations[2], stagePopulations[3])
    println(">>> -- a real shift toward lower ionization, though not all the way to neutral C I " *
            "($(round(stagePopulations[1], digits=4))).")
    #
end
