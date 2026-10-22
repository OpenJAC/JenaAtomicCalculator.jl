## Functions in this file cover: plasma shift and plasma environment effects.
## Alphabetical order within this file.


"""
`TestFrames.testModule_AverageAtom(; short::Bool=true)`
    ... tests the average-atom plasma model of `Plasma.AverageAtomScheme` and the self-consistent field behind it;
        `true` is returned if every check below passes.

        **NOTHING HERE IS COMPARED AGAINST A STORED NUMBER OF OUR OWN.**  Five of the seven checks are exact
        identities that the model must satisfy whatever it computes, and the remaining two are values PUBLISHED by
        W. R. Johnson, C. Guet and G. F. Bertsch, J. Quant. Spectrosc. Radiat. Transfer 99, 327 (2006), Table 1,
        for the case this test runs:  aluminium at 0.27 g/cm^3 and T = 5 eV.  So the test cannot pass on a stale
        reference, and it cannot pass by agreeing with yesterday's version of itself.

        The seven are chosen so that each one fails for a different reason:

        1. the relativistic angular factor summed over the two spin-orbit partners of each l is exactly `2 l_>`,
           which is pure algebra and needs no field.  It fixes the one constant that multiplies every conductivity
           the module returns;
        2. NEUTRALITY -- the cell holds Z electrons;
        3. the RADIAL DENSITY integrates over the cell to the same number the STATE SUM counts.  These are two
           different code paths, `averageAtomDensity` against `averageAtomElectronCount`, and nothing forces them
           to agree;
        4. the FREE box states, counted inside the cell, reproduce the ideal electron gas there.  This is
           completeness, and it is what Johnson's free-electron background term rests on.  It fails by a factor of
           3.3 if the B-spline box is taken to BE the Wigner-Seitz cell, which is the error this model is most
           prone to;
        5. the density MERGES onto that background outside the cell, which a cell-sized box cannot do;
        6. the THOMAS-REICHE-KUHN sum rule for one bound state, with the occupation factors left OUT:
           `sum_b (2/3) w_ab |<kappa_a||C^(1)||kappa_b>|^2 R_ab^2 / (2j_a+1) = 1` exactly, for any state of any
           one-electron Hamiltonian.  It tests the dipole radial integrals, the angular factors and the
           completeness of the basis in a single number, and it is the only check here that would catch a wrong
           matrix element;
        7. the published bound-electron count and chemical potential.

        **WHY THE OCCUPATIONS ARE LEFT OUT OF CHECK 6.**  The conductivity sum rule that `Plasma` reports carries
        the factor `f_a - f_b` and is therefore NOT equal to Z:  for this case it comes to 11.63 against Z = 13,
        the difference being Pauli blocking of transitions into occupied shells, and the 2p even exceeds one per
        electron because downward transitions carry negative oscillator strength.  Asserting Z there would be
        asserting something false.  Dropping the occupations restores an exact identity.

        Runs in about 15 s.  lMax = 3 is enough for every check:  measured 11-Oct-2026, all seven are unchanged to
        the quoted tolerances at lMax = 3, 4 and 5.
"""
function testModule_AverageAtom(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-AverageAtom-new.sum")
    printstyled("\n\nTest the average-atom model of  Plasma  ... \n", color=:cyan)
    success = true;    printTest, iostream = Defaults.getDefaults("test flag/stream")

    # (1) The angular factor, which needs no field at all.  Summed over the two spin-orbit partners of each l,
    #     JAC's relativistic reduced matrix element must give exactly twice the non-relativistic line-strength
    #     factor l_>;  the factor of two is the electron spin.
    for  (la, lb)  in  ((0,1), (1,0), (1,2), (2,1), (2,3), (3,2), (3,4))
        ksa = la == 0 ? [-1] : [-(la+1), la];    ksb = lb == 0 ? [-1] : [-(lb+1), lb]
        wa  = 0.
        for  ka in ksa,  kb in ksb
            wa = wa + AngularMomentum.CL_reduced_me(Subshell(9,ka), 1, Subshell(9,kb))^2
        end
        if  abs(wa - 2*max(la,lb)) > 1.0e-10
            success = false
            if printTest   info(iostream, "sum of |<kappa_a||C^(1)||kappa_b>|^2 over the spin-orbit partners of " *
                                          "l_a = $la, l_b = $lb is $wa, not 2 l_> = $(2*max(la,lb))")   end
        end
    end

    # The one self-consistent field every remaining check is read from:  Johnson, Guet and Bertsch's own example.
    nm    = Nuclear.Model(13.0, 26.98);    temp = Defaults.convertUnits("energy: from eV to atomic", 5.0)
    rWS   = redirect_stdout(devnull) do;   Plasma.determineWignerSeitzRadius(0.27, nm)         end
    grid  = redirect_stdout(devnull) do;   Radial.generateGrid(Radial.Grid(false), boxSize = 3*rWS)   end
    prim  = redirect_stdout(devnull) do;   Bsplines.generatePrimitives(grid)                   end
    subsh = Basics.generateSubshellList( Basics.generateShellList(1, 6, 3) )
    wScf  = redirect_stdout(devnull) do
        orbitals = Bsplines.generateOrbitalsHydrogenic(subsh, nm, prim; printout=false)
        SelfConsistent.solveAverageAtomSpectra(orbitals, nm, Basics.AaDFSField(), temp, rWS, prim; printout=false)
    end
    rho0   = SelfConsistent.freeElectronDensity(wScf.chemMu, temp)
    volume = 4pi/3 * rWS^3

    # (2) Neutrality:  the cell holds Z electrons.
    if  abs(wScf.nBound + wScf.nCont - nm.Z) > 1.0e-8
        success = false
        if printTest   info(iostream, "the cell holds $(wScf.nBound + wScf.nCont) electrons, not Z = $(nm.Z)")   end
    end

    # (3) The radial density and the state sum are two different code paths and must give the same number.
    wInt = RadialIntegrals.V0(wScf.rhot, wScf.iCell, grid)
    if  abs(wInt - (wScf.nBound + wScf.nCont)) > 1.0e-8
        success = false
        if printTest   info(iostream, "the radial density integrates over the cell to $wInt, while the state sum " *
                                      "counts $(wScf.nBound + wScf.nCont);  the two paths disagree")   end
    end

    # (4) Completeness:  the FREE box states, counted inside the cell, are the ideal gas there.
    function freeInCell(spectra, norms)
        wa = 0.
        for  kappa  in  spectra.kappasFree
            occ = 2 * abs(kappa);    wc = spectra.freeEigen[kappa];    wn = norms.freeNorm[kappa]
            for  i = spectra.freeStart[kappa]:length(wc.values)
                wa = wa + occ * Basics.FermiDirac(wc.values[i], wScf.chemMu, temp) * wn[i]
            end
        end
        return( wa )
    end
    wFree = freeInCell(wScf.spectra, wScf.norms)
    if  abs(wFree/(rho0*volume) - 1.0) > 5.0e-3
        success = false
        if printTest   info(iostream, "the free box states hold $wFree electrons inside the cell against an ideal " *
                                      "gas of $(rho0*volume) -- a ratio of $(wFree/(rho0*volume)), which should " *
                                      "be 1;  a box equal to the cell gives about 0.3 here")   end
    end

    # (5) ... and the density merges onto that background outside the cell.
    i15   = something( findlast(r -> r < 1.5*rWS, grid.r), grid.NoPoints )
    wMerge = wScf.rhot[i15] / (4pi * grid.r[i15]^2) / rho0
    if  abs(wMerge - 1.0) > 1.0e-3
        success = false
        if printTest   info(iostream, "at 1.5 R^(WS) the density is $wMerge times the free-electron background, " *
                                      "not 1;  the density does not merge into the background")   end
    end

    # (6) Thomas-Reiche-Kuhn for the 1s, with the occupations left out;  see the docstring.
    spectra = wScf.spectra;    kappaA = -1;    iA = spectra.trueStart[kappaA]
    epsA    = spectra.trueEigen[kappaA].values[iA]
    orbA    = Bsplines.generateOrbitalFromPrimitives(Subshell(100,kappaA), epsA, grid.NoPoints,
                                                     spectra.trueEigen[kappaA].vectors[iA], prim)
    jA2     = Basics.subshell_2j( Subshell(100,kappaA) );    wTrk = 0.
    for  kappaB  in  spectra.kappas
        if  abs(Basics.subshell_l(Subshell(100,kappaA)) - Basics.subshell_l(Subshell(100,kappaB))) != 1   continue  end
        cl = AngularMomentum.CL_reduced_me(Subshell(100,kappaA), 1, Subshell(100,kappaB))
        if  cl == 0.    continue    end
        wc = spectra.trueEigen[kappaB]
        for  ib = spectra.trueStart[kappaB]:length(wc.values)
            orbB = Bsplines.generateOrbitalFromPrimitives(Subshell(100,kappaB), wc.values[ib], grid.NoPoints,
                                                          wc.vectors[ib], prim)
            wR   = RadialIntegrals.rkDiagonal(1, orbA, orbB, grid)
            wTrk = wTrk + (2/3) * (wc.values[ib] - epsA) * cl^2 * wR^2 / (jA2 + 1)
        end
    end
    if  abs(wTrk - 1.0) > 2.0e-2
        success = false
        if printTest   info(iostream, "the Thomas-Reiche-Kuhn sum of the 1s comes to $wTrk and must be 1;  the " *
                                      "dipole radial integrals, the angular factors or the completeness of the " *
                                      "basis are wrong")   end
    end

    # (7) The two published numbers of Johnson, Guet and Bertsch, Table 1.
    for  (name, ours, theirs, tol)  in  ( ("bound electrons",    wScf.nBound, 11.5059, 5.0e-3),
                                          ("chemical potential", wScf.chemMu, -0.3823, 1.5e-2) )
        if  abs(ours/theirs - 1.0) > tol
            success = false
            if printTest   info(iostream, "the $name comes to $ours against the published $theirs of Johnson, " *
                                          "Guet and Bertsch (2006), Table 1 -- a relative difference of " *
                                          "$(abs(ours/theirs - 1.0)), beyond the tolerance $tol")   end
        end
    end

    # (8) THE DC CONDUCTIVITY FROM THE SCATTERING PHASE SHIFTS, which is the one quantity of the optical chain
    #     that used to be supplied from outside.  Three things are asserted, and the first two are identities
    #     that cannot be satisfied by a wrong number:
    #       (a) the transport cross section must be positive and must FALL with energy.  A fast electron is
    #           scattered less, and the 4 pi/k^2 prefactor makes that quantitative:  measured 127.2, 32.0 and
    #           15.6 a.u.^2 at 0.2, 0.8 and 2.0 Ha.  A wrong k-dependence -- the commonest way to get a
    #           transport cross section wrong -- cannot survive this, while a constant factor would.
    #       (b) the sum rule relating gamma to sigma(0) must hold EXACTLY, not approximately.
    #       (c) and the published value, which is what pins the constant factor (a) cannot see:  Johnson, Guet and Bertsch quote sigma(0) = 0.0262 a.u. and check it
    #           against 0.125 (mOhm m)^-1 from the Ziman formula;  this gives 0.0276, i.e. 2 % from the Ziman
    #           number and 5 % from theirs.  The tolerance is set at 12 %, loose enough that the mesh and lMax
    #           defaults may change without re-approving anything and tight enough that a factor-of-two error in
    #           the phase convention, the prefactor or the density of scatterers cannot pass.
    wDC = Plasma.dcConductivity(wScf, temp, rWS, grid)
    wSt = [ Plasma.transportCrossSection(wScf, en, grid)  for en in (0.2, 0.8, 2.0) ]
    if  !(wSt[1] > wSt[2] > wSt[3] > 0.)
        success = false
        if printTest   info(iostream, "the transport cross sections at 0.2, 0.8 and 2.0 Ha come to $wSt and must " *
                                      "be positive and decreasing;  a faster electron is scattered less")   end
    end
    if  abs(wDC.gamma - wScf.nCont/(volume*wDC.sigma0)) > 1.0e-12 * abs(wDC.gamma)
        success = false
        if printTest   info(iostream, "the collision rate $(wDC.gamma) does not satisfy its own sum rule " *
                                      "N_free/(Omega sigma(0)) = $(wScf.nCont/(volume*wDC.sigma0))")   end
    end
    if  abs(wDC.sigma0/0.0262 - 1.0) > 0.12
        success = false
        if printTest   info(iostream, "sigma(0) comes to $(wDC.sigma0) a.u. against the 0.0262 of Johnson, Guet " *
                                      "and Bertsch (2006) and the 0.125 (mOhm m)^-1 Ziman value they check it "  *
                                      "against -- a relative difference of $(abs(wDC.sigma0/0.0262 - 1.0))")   end
    end

    println(iostream, "AverageAtom: the relativistic angular factor against 2 l_> as pure algebra; neutrality of " *
                      "the Wigner-Seitz cell; the radial density and the state sum as two paths to one number; "   *
                      "the free box states reproducing the ideal gas inside the cell, and the density merging "    *
                      "onto that background outside it -- both of which fail if the B-spline box is taken to BE "  *
                      "the cell; the Thomas-Reiche-Kuhn sum of the 1s with the occupations left out, which is an " *
                      "exact identity where the conductivity sum rule is NOT; and the bound-electron count and "   *
                      "chemical potential against Johnson, Guet and Bertsch (2006), Table 1; and the DC "         *
                      "conductivity from the scattering phase shifts, whose transport cross section must fall "   *
                      "with energy, whose sum rule must hold exactly, and whose sigma(0) must agree with the "    *
                      "published value. No approved data and no stored number of our own are used.")
    Defaults.setDefaults("print summary: close", "")
    testPrint("testModule_AverageAtom()::", success)
    return( success )
end


"""
`TestFrames.testModule_SuperConfigurationOpacity(; short::Bool=true)`
    ... tests `Plasma.SuperConfigurationOpacityScheme` (STAGE 1: population weighting only, see that
        struct's own docstring) on C I (Z=6), three configurations of `2s^2 2p^2`'s own complex:
        `2s^2 2p^2`, `2s 2p^3`, `2p^4`. Checks: each group's degeneracy against the EXACT, independently
        known CSF/microstate count (C(2,2)*C(6,2)=15, 2*C(6,3)=40, C(6,4)=15); the returned populations sum
        to 1; and, the real physics check, that in the T -> 0 limit the population concentrates entirely on
        whichever group has the LOWEST mean energy -- the defining property of a Boltzmann distribution, not
        a number borrowed from elsewhere. A success::Bool is returned.
"""
function testModule_SuperConfigurationOpacity(; short::Bool=true)
    printstyled("\n\nTest the module  Plasma (SuperConfigurationOpacityScheme)  ... \n", color=:cyan)
    success = true

    ni   = Nuclear.Model(6.)
    grid = Radial.Grid(true)
    asfSettings = AsfSettings()
    refConfigs  = [Configuration("2s^2 2p^2"), Configuration("2s 2p^3"), Configuration("2p^4")]

    scheme   = Plasma.SuperConfigurationOpacityScheme(refConfigs)
    settings = Plasma.Settings(1.0, 1.0e22, false)
    wa = Plasma.Computation(Plasma.Computation(); scheme=scheme, nuclearModel=ni, grid=grid,
                              refConfigs=refConfigs, asfSettings=asfSettings, settings=settings)
    wb = perform(wa; output=true)

    groups      = wb[ResultKeys.SuperConfigurationGroups]
    populations = wb[ResultKeys.SuperConfigurationPopulations]

    expectedDegeneracies = [15.0, 40.0, 15.0]
    for  (group, expected) in zip(groups, expectedDegeneracies)
        if  group.degeneracy != expected
            success = false
            println("  ** degeneracy for $(group.conf) is $(group.degeneracy), expected $expected")
        end
    end
    if  abs(sum(populations) - 1.0) > 1.0e-10
        success = false;   println("  ** populations do not sum to 1: $(sum(populations))")
    end

    # T -> 0: population must concentrate on the lowest-mean-energy group.
    coldSettings = Plasma.Settings(1.0e-6, 1.0e22, false)
    waCold = Plasma.Computation(wa; settings=coldSettings)
    wbCold = perform(waCold; output=true)
    coldPopulations = wbCold[ResultKeys.SuperConfigurationPopulations]
    groundIdx = argmin( [group.meanEnergy for group in groups] )
    if  coldPopulations[groundIdx] < 0.999999
        success = false
        println("  ** at T -> 0 the lowest-energy group should carry essentially all the population, got " *
                 "$(coldPopulations[groundIdx])")
    end

    testPrint("testModule_SuperConfigurationOpacity()::", success)
    return( success )
end


"""
`TestFrames.testModule_SuperConfigurationIonization(; short::Bool=true)`
    ... tests `Plasma.SuperConfigurationIonizationScheme` on three charge states of carbon (Z=6):
        `1s^2 2s^2 2p^2` (q=0), `1s^2 2s^2 2p` (q=1), `1s^2 2s^2` (q=2). Checks: the stage populations sum
        to 1; electron density equals the charge-neutrality sum built from those very populations (a
        self-consistency floor, not an external reference); and -- the decisive test -- BOTH adjacent-pair
        population ratios independently satisfy the TEXTBOOK two-level Saha equation,
        `n_{q+1}*n_e/n_q = (g_{q+1}/g_q) * 2*(kT/2pi)^1.5 * exp(-IP/kT)`, built here from nothing but the
        SAME groups' own `meanEnergy`/`degeneracy` fields the scheme itself read. The scheme's
        grand-canonical formula and the textbook pairwise formula are mathematically equivalent only if
        implemented consistently; this is what actually proves that identity in code, not merely asserts
        it. A success::Bool is returned.
"""
function testModule_SuperConfigurationIonization(; short::Bool=true)
    printstyled("\n\nTest the module  Plasma (SuperConfigurationIonizationScheme)  ... \n", color=:cyan)
    success = true

    ni   = Nuclear.Model(6.)
    grid = Radial.Grid(true)
    asfSettings = AsfSettings()
    refConfigsPerStage = [ [Configuration("1s^2 2s^2 2p^2")], [Configuration("1s^2 2s^2 2p")],
                            [Configuration("1s^2 2s^2")] ]

    scheme   = Plasma.SuperConfigurationIonizationScheme(refConfigsPerStage)
    kT       = 15.0 / 27.211386   # 15 eV in Hartree, close to carbon's own ionization energies
    settings = Plasma.Settings(kT, 1.0e-8, true)
    wa = Plasma.Computation(Plasma.Computation(); scheme=scheme, nuclearModel=ni, grid=grid, asfSettings=asfSettings,
                              settings=settings)
    wb = perform(wa; output=true)

    groupsPerStage   = wb[ResultKeys.SuperConfigurationIonizationGroups]
    charges          = wb[ResultKeys.SuperConfigurationIonizationCharges]
    stagePopulations = wb[ResultKeys.SuperConfigurationIonizationStagePopulations]
    ne               = wb[ResultKeys.SuperConfigurationIonizationElectronDensity]

    if  abs(sum(stagePopulations) - 1.0) > 1.0e-10
        success = false;   println("  ** stage populations do not sum to 1: $(sum(stagePopulations))")
    end

    neFromNeutrality = 1.0e-8 * sum( charges[i]*stagePopulations[i]   for i = 1:length(charges) )
    if  abs(ne - neFromNeutrality) > 1.0e-10 * ne
        success = false
        println("  ** electron density $ne is not the charge-neutrality sum $neFromNeutrality built from the " *
                 "SAME populations the scheme itself returned")
    end

    for  i = 1:length(charges)-1
        g0 = groupsPerStage[i][1].degeneracy;      E0 = groupsPerStage[i][1].meanEnergy
        g1 = groupsPerStage[i+1][1].degeneracy;    E1 = groupsPerStage[i+1][1].meanEnergy
        sahaRHS = (g1/g0) * 2 * (kT/(2*pi))^1.5 * exp( -(E1-E0)/kT )
        sahaLHS = stagePopulations[i+1] * ne / stagePopulations[i]
        if  abs(sahaLHS - sahaRHS) > 1.0e-6 * sahaRHS
            success = false
            println("  ** pair (q=$(charges[i]), q=$(charges[i+1])): Saha check failed, " *
                     "LHS = $sahaLHS, RHS = $sahaRHS")
        end
    end

    testPrint("testModule_SuperConfigurationIonization()::", success)
    return( success )
end


"""
`TestFrames.testModule_SuperConfigurationMixtureOpacity(; short::Bool=true)`
    ... tests `Plasma.SuperConfigurationMixtureOpacityScheme` on a two-species mixture, carbon (Z=6, the
        SAME three-configuration complex `testModule_SuperConfigurationOpacity` already validates) and neon
        (Z=10, an analogous `2s^2 2p^4`/`2s 2p^5`/`2p^6` complex). Checks (1)-(4) are the ORIGINAL
        (17-Oct-2026) species-axis checks, unaffected by the stage axis added 18-Oct-2026 and run here
        UNEDITED as the backward-compatibility proof: (1) the mixing rule itself, by an independent
        re-derivation -- `kappa_mix` must equal `numberFractions[1]*kappa_C + numberFractions[2]*kappa_Ne`
        to machine precision; (2) the species REDUCTION check -- a one-species-only mixture
        (`numberFractions=[1.0]`) must reproduce `SuperConfigurationOpacityScheme`'s own `kappa(omega)`
        exactly; (3) a parameter-free check on `Plasma.planckMeanOpacity` alone -- a FREQUENCY-INDEPENDENT
        `kappa` must return exactly that constant; (4) that no Rosseland-mean key is present, documenting
        the deliberate finding that it has no stable value for a line-only spectrum.

        CHECKS (5)-(7) are NEW 18-Oct-2026, for the STAGE axis: (5) the STAGE reduction check -- one
        species alone, carried through THREE charge states (the SAME `1s^2 2s^2 2p^2`/`...2p`/`1s^2 2s^2`
        case `testModule_SuperConfigurationIonization` already validates against the textbook Saha
        equation) -- must reproduce that scheme's OWN stage populations and electron density to its own
        1e-10 convergence tolerance, proving the coupled N-species solve reduces correctly to the
        already-validated single-species one; (6) a HOMOGENEITY check -- splitting that SAME species into
        TWO IDENTICAL species at 0.5/0.5 abundance must reproduce the single full-abundance species
        EXACTLY (electron density, both species' own populations, and the mixture spectrum all to machine
        precision) -- the genuinely new coupled-across-species physics, exercised directly; (7) a real
        TWO-STAGE opacity-mixing check (carbon's own `2s^2 2p^2`/`2s 2p^3` complex as one stage, its
        `2s^2 2p`/`2s 2p^2` ionized partner as the next) -- the mixture's own `kappa` must equal the
        STAGE-population-weighted sum of each stage's OWN independently-computed `kappa`, built by an
        independent re-derivation exactly as check (1) does for species, and must differ from EITHER
        stage's own `kappa` alone (proving a genuine mix, not a degenerate one-stage result).

        CHECK (8) IS NEW 18-Oct-2026, for BOUND-FREE: on the SAME stageA/stageB pair as check (7), now with
        `boundFreeApprox=Empirical.ScaledHydrogenic()`, it independently re-derives the scheme's own
        bound-free spectrum by calling `Empirical.photoionizationCrossSection` directly for every
        (initialConf,finalConf) pair differing by exactly one electron in exactly one shell, weighted by
        stageA's own within-stage populations (from `SuperConfigurationOpacityScheme`'s own result for
        stageA alone) times stage A's balance population -- built by a DIFFERENT route than
        `Plasma.boundFreeOpacityContribution` itself uses, to machine precision. Also checks that turning
        bound-free on leaves the bound-bound-only spectrum UNCHANGED from check (7), that the TOTAL spectrum
        is exactly bound-bound+bound-free, and that the ORIGINAL (no `boundFreeApprox`) mixture from this
        test's own opening carries NONE of the new keys -- the backward-compatibility proof for this
        addition, mirroring check (4)'s "no Rosseland key" proof.

        CHECK (9) IS NEW 19-Oct-2026, for FREE-FREE: on the SAME stageA/stageB pair, now with
        `freeFreeIncluded=true` (`boundFreeApprox` left at `nothing`, isolating free-free alone), it
        independently re-derives the scheme's own free-free spectrum by calling
        `Empirical.freeFreeAbsorptionCoefficient` directly for EACH stage's own ion charge -- no adjacency
        needed, unlike bound-free -- weighted by that stage's own balance population, built by a DIFFERENT
        route than `Plasma.freeFreeOpacityContribution` itself uses, to machine precision. Also checks (a)
        turning free-free on leaves the bound-bound-only spectrum UNCHANGED, (b) the TOTAL spectrum is
        exactly bound-bound+free-free, (c) the ORIGINAL (no `freeFreeIncluded`) mixture carries NONE of the
        new keys, and (d) two literal values of the embedded Sutherland (1998) Gaunt-factor table against
        the downloaded file itself -- a data-transcription check independent of the physics above. A
        success::Bool is returned.
"""
function testModule_SuperConfigurationMixtureOpacity(; short::Bool=true)
    printstyled("\n\nTest the module  Plasma (SuperConfigurationMixtureOpacityScheme)  ... \n", color=:cyan)
    success = true

    nmC   = Nuclear.Model(6.)
    nmNe  = Nuclear.Model(10.)
    refsC  = [Configuration("2s^2 2p^2"), Configuration("2s 2p^3"), Configuration("2p^4")]
    refsNe = [Configuration("2s^2 2p^4"), Configuration("2s 2p^5"), Configuration("2p^6")]
    omegas = collect(0.05:0.05:1.2)

    scheme   = Plasma.SuperConfigurationMixtureOpacityScheme([nmC, nmNe], [refsC, refsNe], [0.6, 0.4], omegas)
    settings = Plasma.Settings(0.2, 1.0e22, false)
    wa = Plasma.Computation(Plasma.Computation(); scheme=scheme, nuclearModel=nmC, grid=Radial.Grid(true),
                              asfSettings=AsfSettings(), settings=settings)
    wb = perform(wa; output=true)

    speciesKappas = wb[ResultKeys.SuperConfigurationMixtureSpeciesOpacities]
    kappaMix      = wb[ResultKeys.SuperConfigurationMixtureOpacitySpectrum]
    kappaP        = wb[ResultKeys.SuperConfigurationMixturePlanckOpacity]

    maxMixDiff = maximum( abs.(kappaMix .- (0.6 .* speciesKappas[1] .+ 0.4 .* speciesKappas[2])) )
    if  maxMixDiff > 1.0e-12
        success = false
        println("  ** mixing rule failed: max |kappa_mix - (0.6*kappa_C+0.4*kappa_Ne)| = $maxMixDiff")
    end

    schemeOnlyC = Plasma.SuperConfigurationMixtureOpacityScheme([nmC], [refsC], [1.0], omegas)
    waOnlyC     = Plasma.Computation(wa; scheme=schemeOnlyC)
    wbOnlyC     = perform(waOnlyC; output=true)
    kappaOnlyC  = wbOnlyC[ResultKeys.SuperConfigurationMixtureOpacitySpectrum]

    schemeSingle = Plasma.SuperConfigurationOpacityScheme(refsC, omegas)
    gridC        = Basics.recommendedGrid(refsC, nmC; printout=false)
    waSingle     = Plasma.Computation(wa; scheme=schemeSingle, nuclearModel=nmC, grid=gridC, refConfigs=refsC)
    wbSingle     = perform(waSingle; output=true)
    kappaSingle  = wbSingle[ResultKeys.SuperConfigurationOpacitySpectrum]

    maxReductionDiff = maximum(abs.(kappaOnlyC .- kappaSingle))
    if  maxReductionDiff > 1.0e-12
        success = false
        println("  ** reduction check failed: max |mixture(1 species) - single-species scheme| = $maxReductionDiff")
    end

    # a DEDICATED, wide synthetic x-grid for this check -- the physical `omegas`/`kT` above cover only the
    # lines themselves (x up to ~6) and are not meant to resolve planckWeight's own tail out past x~20-30.
    xGridWide  = collect(0.01:0.01:30.0)
    kappaConst = fill(2.3, length(xGridWide))
    planckOfConst = Plasma.planckMeanOpacity(xGridWide, kappaConst)
    if  abs(planckOfConst - 2.3) > 1.0e-4
        success = false
        println("  ** Planck mean of a constant kappa=2.3 returned $planckOfConst, not 2.3")
    end

    if  haskey(wb, "superconfiguration mixture rosseland opacity")
        success = false
        println("  ** a Rosseland-mean key is present -- it should not be, see this scheme's own docstring")
    end

    # --- CHECK (5): STAGE reduction -- one species, three charge states, against
    #     SuperConfigurationIonizationScheme's OWN already-validated numbers for the identical case.
    refConfigsPerStage = [ [Configuration("1s^2 2s^2 2p^2")], [Configuration("1s^2 2s^2 2p")],
                            [Configuration("1s^2 2s^2")] ]
    ionScheme   = Plasma.SuperConfigurationIonizationScheme(refConfigsPerStage)
    kTIon       = 15.0 / 27.211386
    settingsIon = Plasma.Settings(kTIon, 1.0e-8, true)
    waIon       = Plasma.Computation(wa; scheme=ionScheme, settings=settingsIon)
    wbIon       = perform(waIon; output=true)
    popIon      = wbIon[ResultKeys.SuperConfigurationIonizationStagePopulations]
    neIon       = wbIon[ResultKeys.SuperConfigurationIonizationElectronDensity]

    mixSchemeOneSpecies = Plasma.SuperConfigurationMixtureOpacityScheme([nmC], [refConfigsPerStage], [1.0], omegas)
    waMixOne  = Plasma.Computation(wa; scheme=mixSchemeOneSpecies, settings=settingsIon)
    wbMixOne  = perform(waMixOne; output=true)
    popMixOne = wbMixOne[ResultKeys.SuperConfigurationMixtureStagePopulations][1]
    neMixOne  = wbMixOne[ResultKeys.SuperConfigurationMixtureElectronDensity]

    if  maximum(abs.(popIon .- popMixOne)) > 1.0e-8   ||   abs(neIon-neMixOne) > 1.0e-8*neIon
        success = false
        println("  ** stage-reduction check failed: populations $popIon vs $popMixOne, ne $neIon vs $neMixOne")
    end

    # --- CHECK (6): HOMOGENEITY -- the SAME species split into two identical halves (0.5/0.5) must
    #     reproduce the single full-abundance species exactly: electron density, both halves' own
    #     populations (equal to each other AND to the single-species case), and the mixture spectrum.
    mixSchemeSplit = Plasma.SuperConfigurationMixtureOpacityScheme([nmC, nmC],
                                                                       [refConfigsPerStage, refConfigsPerStage],
                                                                       [0.5, 0.5], omegas)
    waMixSplit = Plasma.Computation(wa; scheme=mixSchemeSplit, settings=settingsIon)
    wbMixSplit = perform(waMixSplit; output=true)
    neMixSplit = wbMixSplit[ResultKeys.SuperConfigurationMixtureElectronDensity]
    popSplitA, popSplitB = wbMixSplit[ResultKeys.SuperConfigurationMixtureStagePopulations]
    kappaMixSplit = wbMixSplit[ResultKeys.SuperConfigurationMixtureOpacitySpectrum]
    kappaMixOne   = wbMixOne[ResultKeys.SuperConfigurationMixtureOpacitySpectrum]

    if  abs(neMixOne-neMixSplit) > 1.0e-10*neMixOne   ||   maximum(abs.(popSplitA.-popSplitB)) > 1.0e-10   ||
        maximum(abs.(popMixOne.-popSplitA)) > 1.0e-10   ||   maximum(abs.(kappaMixOne.-kappaMixSplit)) > 1.0e-10
        success = false
        println("  ** homogeneity check failed: splitting one species into two identical halves moved the result")
    end

    # --- CHECK (7): a REAL two-stage opacity mix -- each stage now has its OWN E1-connected pair, so
    #     kappa is genuinely non-zero per stage, and the across-stage population weighting does something.
    stageA = [Configuration("2s^2 2p^2"), Configuration("2s 2p^3")]   # 4 electrons (C2+)
    stageB = [Configuration("2s^2 2p"),   Configuration("2s 2p^2")]   # 3 electrons (C3+)
    twoStageScheme = Plasma.SuperConfigurationMixtureOpacityScheme([nmC], [[stageA, stageB]], [1.0], omegas)
    waTwoStage = Plasma.Computation(wa; scheme=twoStageScheme, settings=settingsIon)
    wbTwoStage = perform(waTwoStage; output=true)
    popTwoStage   = wbTwoStage[ResultKeys.SuperConfigurationMixtureStagePopulations][1]
    kappaTwoStage = wbTwoStage[ResultKeys.SuperConfigurationMixtureOpacitySpectrum]

    gridA = Basics.recommendedGrid(stageA, nmC; printout=false)
    schemeA = Plasma.SuperConfigurationOpacityScheme(stageA, omegas)
    waA = Plasma.Computation(wa; scheme=schemeA, nuclearModel=nmC, grid=gridA, refConfigs=stageA, settings=settingsIon)
    kappaA = perform(waA; output=true)[ResultKeys.SuperConfigurationOpacitySpectrum]

    gridB = Basics.recommendedGrid(stageB, nmC; printout=false)
    schemeB = Plasma.SuperConfigurationOpacityScheme(stageB, omegas)
    waB = Plasma.Computation(wa; scheme=schemeB, nuclearModel=nmC, grid=gridB, refConfigs=stageB, settings=settingsIon)
    kappaB = perform(waB; output=true)[ResultKeys.SuperConfigurationOpacitySpectrum]

    handMixed = popTwoStage[1].*kappaA .+ popTwoStage[2].*kappaB
    if  maximum(abs.(kappaTwoStage .- handMixed)) > 1.0e-10
        success = false
        println("  ** two-stage mixing rule failed: max |kappa - (p1*kappaA+p2*kappaB)| = " *
                 "$(maximum(abs.(kappaTwoStage .- handMixed)))")
    end
    if  maximum(abs.(kappaTwoStage.-kappaA)) < 1.0e-6   ||   maximum(abs.(kappaTwoStage.-kappaB)) < 1.0e-6
        success = false
        println("  ** two-stage mix is degenerate -- it equals one stage's own kappa alone")
    end

    # --- CHECK (8), NEW 18-Oct-2026: BOUND-FREE, on the SAME stageA/stageB pair as check (7), now with
    #     boundFreeApprox turned on. Independently re-derives the scheme's own bound-free spectrum by
    #     calling Empirical.photoionizationCrossSection directly for every (iConf,fConf) pair that differs
    #     by exactly one electron in exactly one shell, weighted by stageA's own WITHIN-stage populations
    #     (from SuperConfigurationOpacityScheme's own result for stageA alone) times the balance population
    #     of stage A itself -- EXACTLY the quantity Plasma.boundFreeOpacityContribution computes, built here
    #     by a different route. Also checks that (a) the bound-bound-only keys are UNCHANGED from check (7)'s
    #     own numbers, (b) total = bound-bound + bound-free exactly, and (c) the ORIGINAL (no-boundFreeApprox)
    #     mixture from the top of this test carries NONE of the new keys -- the backward-compatibility proof
    #     for this addition, mirroring check (4)'s "no Rosseland key" proof for the Rosseland mean itself.
    popsA = perform(waA; output=true)[ResultKeys.SuperConfigurationPopulations]

    bfScheme = Plasma.SuperConfigurationMixtureOpacityScheme([nmC], [[stageA, stageB]], [1.0], omegas,
                                                                 Empirical.ScaledHydrogenic())
    waBF = Plasma.Computation(wa; scheme=bfScheme, settings=settingsIon)
    wbBF = perform(waBF; output=true)
    kappaBound     = wbBF[ResultKeys.SuperConfigurationMixtureOpacitySpectrum]
    kappaBF        = wbBF[ResultKeys.SuperConfigurationMixtureBoundFreeOpacity]
    kappaTotalBF   = wbBF[ResultKeys.SuperConfigurationMixtureTotalOpacitySpectrum]
    popTwoStageBF  = wbBF[ResultKeys.SuperConfigurationMixtureStagePopulations][1]

    if  maximum(abs.(kappaBound .- kappaTwoStage)) > 1.0e-10
        success = false
        println("  ** turning boundFreeApprox ON changed the BOUND-BOUND-only spectrum -- it must not")
    end

    zOld = Defaults.getDefaults("nuclear: charge");   Defaults.setDefaults("nuclear: charge", nmC.Z)
    handBF = zeros(Float64, length(omegas))
    try
        for  (g, iConf)  in  enumerate(stageA)
            popAbs = popTwoStageBF[1] * popsA[g]
            popAbs <= 0.0   &&   continue
            for  fConf  in  stageB
                occDiff = Basics.extractFromConfigurations(Basics.OccupationDifference(), iConf, fConf)
                ( length(occDiff) == 1  &&  occDiff[1][2] == 1 )   ||   continue
                css = Empirical.photoionizationCrossSection(omegas, iConf, fConf, Empirical.ScaledHydrogenic())
                handBF .+= popAbs .* css
            end
        end
    finally
        Defaults.setDefaults("nuclear: charge", zOld)
    end

    if  maximum(handBF) <= 0.0
        success = false
        println("  ** bound-free check is vacuous -- no compatible channel found between stageA and stageB")
    end
    if  maximum(abs.(kappaBF .- handBF)) > 1.0e-10
        success = false
        println("  ** bound-free check failed: max |kappaBF - hand-derived| = $(maximum(abs.(kappaBF.-handBF)))")
    end
    if  maximum(abs.(kappaTotalBF .- (kappaBound.+kappaBF))) > 1.0e-12
        success = false
        println("  ** total spectrum is not exactly bound-bound + bound-free")
    end
    if  haskey(wb, "superconfiguration mixture bound-free opacity")
        success = false
        println("  ** the ORIGINAL (no boundFreeApprox) mixture carries a bound-free key -- it must not")
    end

    # --- CHECK (9), NEW 19-Oct-2026: FREE-FREE, on the SAME stageA/stageB pair, now with freeFreeIncluded
    #     turned on (boundFreeApprox left at nothing, to isolate free-free alone). Independently re-derives
    #     the scheme's own free-free spectrum by calling Empirical.freeFreeAbsorptionCoefficient directly
    #     for EACH stage's own charge -- no adjacency needed, unlike bound-free -- weighted by that stage's
    #     own balance population, EXACTLY the quantity Plasma.freeFreeOpacityContribution computes, built
    #     here by a different route. Also checks (a) turning it on leaves the bound-bound-only spectrum
    #     UNCHANGED, (b) total = bound-bound + free-free exactly, (c) the ORIGINAL (freeFreeIncluded=false)
    #     mixture carries NONE of the new keys, and (d) two literal values of the embedded Sutherland (1998)
    #     Gaunt-factor table -- a transcription check independent of the physics above, since a wrong digit
    #     anywhere in 3321 embedded numbers would not necessarily show up as an obviously wrong cross section.
    ffScheme = Plasma.SuperConfigurationMixtureOpacityScheme([nmC], [[stageA, stageB]], [1.0], omegas,
                                                                 nothing, true)
    waFF = Plasma.Computation(wa; scheme=ffScheme, settings=settingsIon)
    wbFF = perform(waFF; output=true)
    kappaBoundFF    = wbFF[ResultKeys.SuperConfigurationMixtureOpacitySpectrum]
    kappaFF         = wbFF[ResultKeys.SuperConfigurationMixtureFreeFreeOpacity]
    kappaTotalFF    = wbFF[ResultKeys.SuperConfigurationMixtureTotalOpacitySpectrum]
    popTwoStageFF   = wbFF[ResultKeys.SuperConfigurationMixtureStagePopulations][1]
    neFF            = wbFF[ResultKeys.SuperConfigurationMixtureElectronDensity]

    if  maximum(abs.(kappaBoundFF .- kappaTwoStage)) > 1.0e-10
        success = false
        println("  ** turning freeFreeIncluded ON changed the BOUND-BOUND-only spectrum -- it must not")
    end

    chargesAB = [6 - stageA[1].NoElectrons, 6 - stageB[1].NoElectrons]
    handFF    = zeros(Float64, length(omegas))
    for  (i, q)  in  enumerate(chargesAB)
        handFF .+= popTwoStageFF[i] .* Empirical.freeFreeAbsorptionCoefficient(omegas, q, neFF, kTIon)
    end

    if  maximum(handFF) <= 0.0
        success = false
        println("  ** free-free check is vacuous -- both stage charges gave zero")
    end
    if  maximum(abs.(kappaFF .- handFF)) > 1.0e-10
        success = false
        println("  ** free-free check failed: max |kappaFF - hand-derived| = $(maximum(abs.(kappaFF.-handFF)))")
    end
    if  maximum(abs.(kappaTotalFF .- (kappaBoundFF.+kappaFF))) > 1.0e-12
        success = false
        println("  ** total spectrum is not exactly bound-bound + free-free")
    end
    if  haskey(wb, "superconfiguration mixture free-free opacity")
        success = false
        println("  ** the ORIGINAL (no freeFreeIncluded) mixture carries a free-free key -- it must not")
    end

    # literal spot values from Sutherland's (1998) gffgu.dat itself, at its own grid nodes (u,gamma^2) =
    # (1e-4,1e-4) and (1.0,1.0), read directly off the downloaded, md5-verified file -- NOT recomputed from
    # any formula, so this check is independent of everything else in this test.
    gffCorner = Empirical.freeFreeGauntFactor(1.0e-4, 1.0e-4)
    gffCenter = Empirical.freeFreeGauntFactor(1.0, 1.0)
    if  abs(gffCorner - 5.524282) > 1.0e-5   ||   abs(gffCenter - 1.293904) > 1.0e-5
        success = false
        println("  ** embedded Sutherland table spot check failed: gff(1e-4,1e-4) = $gffCorner " *
                 "(expect 5.524282), gff(1,1) = $gffCenter (expect 1.293904)")
    end

    testPrint("testModule_SuperConfigurationMixtureOpacity()::", success)
    return( success )
end


"""
`TestFrames.testModule_KineticChargeBalance(; short::Bool=true)`
    ... tests `Plasma.KineticChargeBalanceScheme` on a 3-stage carbon chain (C^3+/C^4+/C^5+), restricted to
        `Plasma.EmpiricalTreatment` for every leg so the test stays fast (no SCF/CI). Checks: the scheme runs
        through the generic `perform()` route; the returned (Coulomb, Babushkin) population vectors are
        non-negative and sum to 1; at T_e = 1e4 K the lowest-charge stage's ground slot holds essentially all
        the population (carbon's ionization potentials vastly exceed kT there); the `nStages < 2` guard raises;
        and -- the regression this test exists for -- the global "nuclear: charge" default is restored to
        whatever it was BEFORE the call, not left at the scheme's own nuclearModel.Z (Empirical's plasma-rate
        functions read that global default rather than taking Z as an argument; a run that leaves it changed
        would silently corrupt whatever the caller computes next in the same session). A SECOND, minimal case
        then exercises the part of the scheme the carbon chain above never touches: a GENERAL (non-classic)
        `Plasma.ChargeStage` -- one carrying a single configuration with 2 retained CI levels -- paired with
        `Plasma.StructureBasedTreatment` on every leg, checking it completes with non-negative, sum-to-1
        populations; and `Plasma.validateKineticChargeBalanceSupport`'s own loud guard, checking it RAISES when
        that same general stage is instead paired with `Plasma.EmpiricalTreatment` on `peTreatment`. A
        success::Bool is returned.
"""
function testModule_KineticChargeBalance(; short::Bool=true)
    printstyled("\n\nTest the module  Plasma (KineticChargeBalanceScheme)  ... \n", color=:cyan)
    success = true

    ni   = Nuclear.Model(6.)
    grid = Radial.Grid(true)

    stages = [ Plasma.ChargeStage(Configuration("1s^2 2s^1"), Configuration("1s^2 2p^1")),
               Plasma.ChargeStage(Configuration("1s^2"),      Configuration("1s^1 2p^1")),
               Plasma.ChargeStage(Configuration("1s^1"),      Configuration("2p^1"))      ]
    scheme = Plasma.KineticChargeBalanceScheme(Plasma.KineticChargeBalanceScheme(); chargeStages=stages, electronTemperatures=[1.0e4])
    comp   = Plasma.Computation(Plasma.Computation(); scheme=scheme, nuclearModel=ni, grid=grid,
                                settings=Plasma.Settings(0., 1.0e10, true))

    oldZ = Defaults.getDefaults("nuclear: charge");   Defaults.setDefaults("nuclear: charge", 99.0)
    wb = redirect_stdout(devnull) do
             perform(comp; output=true)
         end
    if  Defaults.getDefaults("nuclear: charge") != 99.0
        success = false
        println("  ** the global \"nuclear: charge\" default was left at $(Defaults.getDefaults("nuclear: charge")) " *
                 "instead of being restored to 99.0")
    end
    Defaults.setDefaults("nuclear: charge", oldZ)

    popVectors = wb["populations:"]
    pv = popVectors[1.0e4]
    for  (label, pop)  in  (("Coulomb", pv.Coulomb), ("Babushkin", pv.Babushkin))
        if  any(pop .< 0.)
            success = false;   println("  ** $label populations contain a negative entry: $pop")
        end
        if  abs(sum(pop) - 1.0) > 1.0e-8
            success = false;   println("  ** $label populations do not sum to 1: $(sum(pop))")
        end
        if  pop[1] < 0.999999
            success = false
            println("  ** at T_e = 1e4 K the lowest-charge stage's ground slot should hold essentially all the " *
                     "population, got $label = $(pop[1])")
        end
    end

    oneStageScheme = Plasma.KineticChargeBalanceScheme(scheme; chargeStages=[stages[1]])
    oneStageComp   = Plasma.Computation(comp; scheme=oneStageScheme)
    try
        redirect_stdout(devnull) do
            perform(oneStageComp; output=true)
        end
        success = false;   println("  ** a 1-stage KineticChargeBalanceScheme should raise, but did not")
    catch  ex
        if  ex isa UndefVarError || ex isa MethodError
            success = false;   println("  ** the 1-stage guard raised the wrong kind of error: $(typeof(ex))")
        end
    end

    # A GENERAL (non-classic) stage pair, minimal: the higher-charge stage carries TWO ionization-product
    # configurations -- "1s^2" (remove the 2s electron) and "1s^1 2s^1" (remove a 1s electron instead), the
    # second with 2 retained CI levels (J=0/J=1, same occupation) -- against a single-configuration lower-charge
    # neighbour, every leg StructureBasedTreatment. This is the one shape the carbon-chain case above never
    # exercises, cf. Plasma.ChargeStage, and the two product configurations are chosen so every (level, level)
    # pair against the neighbour is a genuine single-electron ionization (cf. ThreeBodyRecombination.capturedSubshell):
    # a general stage whose extra levels instead came from splitting ONE shell's own j (e.g. "1s^0 2p^1" with
    # nLevels=2 giving 2p_1/2 and 2p_3/2) was tried first and correctly REJECTED by capturedSubshell -- pairing its
    # 2p_3/2 level against a neighbour fixed at 2p_1/2 is not a single capture, it also silently changes the
    # spectator electron's own j, and this is exactly the "user's own responsibility" caveat in Plasma.ChargeStage's
    # own docstring, not a bug to route around here.
    generalStages = [ Plasma.ChargeStage([Configuration("1s^2 2s^1")]; nLevels=[1]),
                      Plasma.ChargeStage([Configuration("1s^2"), Configuration("1s^1 2s^1")]; nLevels=[1,2]) ]
    generalScheme = Plasma.KineticChargeBalanceScheme(Plasma.KineticChargeBalanceScheme();
                        chargeStages=generalStages, electronTemperatures=[1.0e4],
                        peTreatment =Plasma.StructureBasedTreatment(), ieTreatment=Plasma.StructureBasedTreatment(),
                        eiiTreatment=Plasma.StructureBasedTreatment(), rrTreatment=Plasma.StructureBasedTreatment(),
                        tbrTreatment=Plasma.StructureBasedTreatment())
    # ImpactExcitation/PhotoRecombination need a continuum, which needs a grid with an explicit hp (cf. Rule 12);
    # the carbon-chain grid above has none, since the EmpiricalTreatment legs it tests never touch a continuum.
    # This box is deliberately NOT tuned against the orbitals (that is a job for /test Mx on a real example, not a
    # mechanism-only scaffold test): it may print ImpactExcitation's own "slowest continuum electron does not fit
    # inside the box" warning at the lowest mesh energy, which is expected here and does not affect what this test
    # actually checks (non-negativity, sum-to-1, and the guard), cf. testModule_Liouville's own "NO physics is
    # compared against approved data" disclaimer for the same kind of scaffold test.
    generalGrid   = Radial.Grid(Radial.Grid(true), rnt=4.0e-6, h=5.0e-2, hp=1.0e-2, rbox=30.0)
    generalComp   = Plasma.Computation(Plasma.Computation(); scheme=generalScheme, nuclearModel=ni, grid=generalGrid,
                                       settings=Plasma.Settings(0., 1.0e10, true))
    wg = redirect_stdout(devnull) do
             perform(generalComp; output=true)
         end
    popG = wg["populations:"][1.0e4]
    for  (label, pop)  in  (("Coulomb", popG.Coulomb), ("Babushkin", popG.Babushkin))
        if  any(pop .< 0.)
            success = false;   println("  ** general-stage $label populations contain a negative entry: $pop")
        end
        if  abs(sum(pop) - 1.0) > 1.0e-8
            success = false;   println("  ** general-stage $label populations do not sum to 1: $(sum(pop))")
        end
    end

    # The SAME general stages, but with peTreatment left at Plasma.EmpiricalTreatment() -- a stage/treatment
    # combination Plasma.validateKineticChargeBalanceSupport must refuse rather than guess at.
    badScheme = Plasma.KineticChargeBalanceScheme(generalScheme; peTreatment=Plasma.EmpiricalTreatment(Empirical.ScaledHydrogenic()))
    try
        Plasma.validateKineticChargeBalanceSupport(badScheme)
        success = false;   println("  ** a general stage with peTreatment = EmpiricalTreatment should raise, but did not")
    catch  ex
        if  !(ex isa ErrorException)
            success = false;   println("  ** the general-stage guard raised the wrong kind of error: $(typeof(ex))")
        end
    end

    testPrint("testModule_KineticChargeBalance()::", success)
    return( success )
end


## PlasmaShift: no test exists and none is planned here. The former testModule_PlasmaShift was deleted
## on 09-Aug-2026 as obsolete: PlasmaShift is no longer an Atomic property but a Plasma.Computation
## scheme, so a test for it belongs with the Plasma module and has to be written afresh.
