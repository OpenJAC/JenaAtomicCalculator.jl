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


## PlasmaShift: no test exists and none is planned here. The former testModule_PlasmaShift was deleted
## on 09-Aug-2026 as obsolete: PlasmaShift is no longer an Atomic property but a Plasma.Computation
## scheme, so a test for it belongs with the Plasma module and has to be written afresh.
