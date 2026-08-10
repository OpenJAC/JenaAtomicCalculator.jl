#
println("Hb) Tests of (empirical) impact-ionization cross sections.")
setDefaults("unit: energy", "eV")
#
setDefaults("print summary: open", "S-plus.sum")
grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)

if  false
    #
    # Last successful:  14-Jul-2026
    # Calculate partial EII cross section for K-subshell of Ne IX in different models
    # BEBmodel(), BEDmodel(), RelativisticBEBmodel(), RelativisticBEDmodel(), FittedBEDmodel()
    approx      = ImpactIonization.RelativisticBEBmodel()
    multipleN   = 1
    iEnergies   = [1200., 1400., 1600., 1800., 2000., 3000., 4000., 5000., 6000., 8000., 
                   10000., 20000., 40000., 60000., 80000., 100000.] ## unit: eV. The incident energies should be > epsilon_subshell.
    shells      = Basics.generateShellList(1,1, [0])
    selection   = ShellSelection(true, shells, Int64[])
    configs     = [Configuration("1s^2")]
    name        = "EII cross section for K-subshell of Ne IX."
    nucModel    = Nuclear.Model(10.0)
    eiiSettings = ImpactIonization.Settings(approx, multipleN, iEnergies, true, true, selection)
    comp        = Empirical.Computation(name, nucModel, grid, configs, eiiSettings)
    #
    perform(comp; output=true)
    #
    setDefaults("print summary: close", "")
    #
elseif false
    #
    # Last successful:  14-Jul-2026
    # Calculate partial EII cross section for M1-M3 subshells (3s, 3p) of Kr XIX in different models.
    # Kr XIX = Kr^18+ with [Ar] ground config; Z=36 gives M-shell binding energies 759-831 eV.
    # BEBmodel(), BEDmodel(), RelativisticBEBmodel(), RelativisticBEDmodel(), FittedBEDmodel()
    approx      = ImpactIonization.RelativisticBEBmodel()
    multipleN   = 1
    iEnergies   = [800, 1000, 1200., 1400., 1600., 1800., 2000., 3000., 4000., 5000., 6000., 8000., 
                   10000., 20000., 40000., 60000., 80000., 100000.] ## unit: eV. The incident energies should be > epsilon_subshell.
    shells      = Basics.generateShellList(3,3, [0,1])
    selection   = ShellSelection(true, shells, Int64[])
    configs     = [Configuration("[Ar]")]
    name        = "EII cross section for M-subshell of Kr XIX."
    nucModel    = Nuclear.Model(36.0)
    eiiSettings = ImpactIonization.Settings(approx, multipleN, iEnergies, true, true, selection)
    comp        = Empirical.Computation(name, nucModel, grid, configs, eiiSettings)
    #
    perform(comp; output=true)
    #
    setDefaults("print summary: close", "")
    #
elseif false
    #
    # Last successful:  14-Jul-2026
    # Calculate partial EII cross section for M-subshells of U in different models
    approx      = ImpactIonization.FittedBEDmodel()
    multipleN   = 1
    iEnergies   = [6000., 8000., 10000., 12000., 14000., 16000., 18000., 20000., 30000., 40000., 60000., 70000., 
                   80000., 90000., 100000., 120000., 140000., 160000., 180000.,  200000., 500000., 1000000., 2000000., 
                   3000000., 5000000., 10000000.] ## unit: eV. The incident energies should be > epsilon_subshell.
    shells      = Basics.generateShellList(3,3, [0,1,2])
    selection   = ShellSelection(true, shells, Int64[])
    configs     = [Configuration("[Rn] 5f^3 6d 7s^2")]
    name        = "EII cross section for M-subshells of U."
    nucModel    = Nuclear.Model(92.0)
    eiiSettings = ImpactIonization.Settings(approx, multipleN, iEnergies, false, true, selection)
    comp        = Empirical.Computation(name, nucModel, grid, configs, eiiSettings)
    #
    perform(comp; output=true)
    #
    setDefaults("print summary: close", "")
    #
elseif false
    #
    # Last successful:  14-Jul-2026
    # Calculate partial EII cross section for for L-subshells of Bi in different models
    approx      = ImpactIonization.BEBmodel() 
    multipleN   = 1
    iEnergies   = [10000.0, 12000., 14000., 16000., 18000., 20000., 30000., 40000., 50000., 60000., 70000., 80000., 
                   100000., 200000., 400000., 600000., 800000., 1000000., 2000000., 4000000., 6000000., 8000000., 
                   10000000., 100000000, 1000000000] ## unit: eV. The incident energies should be > epsilon_subshell.                
    shells      = Basics.generateShellList(2,2, [0,1])  
    selection   = ShellSelection(true, shells=shells)
    configs     = [Configuration("[Xe] 4f^14 5d^10 6s^2 6p^3")] 
    name        = "EII cross section for L-subshells of Bi."
    nucModel    = Nuclear.Model(83.0)
    eiiSettings = ImpactIonization.Settings( approx, multipleN, iEnergies, false, true, selection )
    comp        = Empirical.Computation(name, nucModel, grid, configs, eiiSettings)
    #
    perform(comp; output=true)
    #
    setDefaults("print summary: close", "")
    #
elseif false
    #
    # Last successful:  14-Jul-2026
    # Calculate total EIMI cross section for Be-like Boron (B+) via DirectMultipleModel.
    # Bug fix (14-Jul-2026): Belenger formula outputs Mb, so multiply by 1e6 before /Barn2Au.
    # multipleN=1 now errors correctly; multipleN=2..10 uses param table; >10 uses fit.
    approx      = ImpactIonization.DirectMultipleModel()
    multipleN   = 2
    iEnergies   = [1., 1200., 1400., 1600., 1800., 2000., 3000., 4000., 5000., 6000., 8000.,
                   10000., 20000., 40000., 60000., 80000., 100000.]
    shells      = Basics.generateShellList(1,1, [0])
    selection   = ShellSelection(true, shells, Int64[])
    configs     = [Configuration("1s^2 2s^2")]
    name        = "EIMI double-ionization cross section for B+ via DirectMultipleModel."
    nucModel    = Nuclear.Model(5.0)
    eiiSettings = ImpactIonization.Settings(approx, multipleN, iEnergies, false, true, selection)
    comp        = Empirical.Computation(name, nucModel, grid, configs, eiiSettings)
    #
    perform(comp; output=true)
    #
    setDefaults("print summary: close", "")
    #
elseif false
    #
    # Last successful:  14-Jul-2026
    # Test DoubleExperimentModel on He I — verifies the NaN fix (c_gamma=0 guard)
    # and the (x3+phi_gamma) correction replacing the hardcoded (x3+5).
    # Before fix: sigma_indir = NaN for He (phi_gamma=0 → 0.3/0/log(...)=Inf, 0*Inf=NaN).
    # After fix:  sigma_indir = 0 for He (c_gamma=0 → continue), result = sigma_dir only.
    approx      = ImpactIonization.DoubleExperimentModel()
    multipleN   = 2
    iEnergies   = [80., 100., 150., 200., 300., 500., 1000., 2000., 5000.]
    shells      = Basics.generateShellList(1,1, [0])
    selection   = ShellSelection(true, shells, Int64[])
    configs     = [Configuration("1s^2")]
    name        = "EIMI double-ionization cross section for He I via DoubleExperimentModel."
    nucModel    = Nuclear.Model(2.0)
    eiiSettings = ImpactIonization.Settings(approx, multipleN, iEnergies, false, true, selection)
    comp        = Empirical.Computation(name, nucModel, grid, configs, eiiSettings)
    #
    perform(comp; output=true)
    #
    setDefaults("print summary: close", "")
    #
elseif false
    #
    # Last successful:  14-Jul-2026
    # Test LotzMultipleModel on neutral Boron — verifies unit fix (Lotz2Au) and t_inner fix.
    # Before fix: t = epsilon/totalEnergy (wrong); coefficient 4.5*1e4/Barn2Au (wrong units).
    # After fix:  t_inner = epsilon/I_gamma; Lotz2Au = 4.5e-14/(Barn2Au*1e-24) a.u./eV².
    # Inner 1s threshold for B = 200.33 eV; energies below give sigma=0.
    approx      = ImpactIonization.LotzMultipleModel()
    multipleN   = 2
    iEnergies   = [100., 200., 300., 500., 1000., 2000., 5000., 10000.]
    shells      = Basics.generateShellList(1,1, [0])
    selection   = ShellSelection(true, shells, Int64[])
    configs     = [Configuration("1s^2 2s^2 2p^1")]
    name        = "EIMI double-ionization cross section for B I via LotzMultipleModel."
    nucModel    = Nuclear.Model(5.0)
    eiiSettings = ImpactIonization.Settings(approx, multipleN, iEnergies, false, true, selection)
    comp        = Empirical.Computation(name, nucModel, grid, configs, eiiSettings)
    #
    perform(comp; output=true)
    #
    setDefaults("print summary: close", "")
    #
elseif true
    #
    # Last visit:  10-Aug-2026
    # Last successful:  10-Aug-2026
    # EII of an OPEN-SHELL target: C I (1s^2 2s^2 2p^2) with the BEB model.
    #
    # WHY THIS BRANCH EXISTS. Every other EII test and example here uses a CLOSED shell -- He 1s^2, Ne IX
    # 1s^2, Kr XIX [Ar] -- whose basis holds a SINGLE CSF. TestFrames.testModule_ImpactIonization uses
    # He 1s^2 for the same reason. That is a blind spot: the BEB cross section is LINEAR in the mean
    # subshell occupation, and Basics.computeMeanSubshellOccupation(sh, basis) is a per-CSF AVERAGE, so a
    # defect in that averaging is exactly inert when there is only one CSF to average over. One was present
    # until 10-Aug-2026 (commit cef7c63): a spurious outer loop made the function return the SUM over CSFs
    # instead of their mean, i.e. a result too large by the number of CSFs. The C I basis below holds 5
    # CSFs, and its 2p cross sections were correspondingly 5.2-7.3x too large -- ENERGY-DEPENDENTLY so, so
    # the shape of the curve was wrong too, not merely its scale. No closed-shell case could ever have shown
    # this.
    #
    # WHAT TO CHECK, and it is an EXACT invariant rather than a tolerance: the mean subshell occupations
    # must sum to the number of electrons. That equality is what the defect violated (it gave 30.0 for the
    # 6 electrons of C I), and it holds independently of any model, grid or measurement. The branch prints
    # the occupations and their sum for exactly this reason. Observed 10-Aug-2026:
    #     1s_1/2 = 2.0,  2s_1/2 = 2.0,  2p_1/2 = 0.8,  2p_3/2 = 1.2,  summing to 6.0 for 6 electrons.
    # Note that 0.8 : 1.2 is NOT the degeneracy-weighted statistical ratio 2 : 4 (which would give
    # 0.667 : 1.333). It should not be: Basics.computeMeanSubshellOccupation averages over the CSFs of the
    # basis with EQUAL weight, not with the weight of each subshell's degeneracy. The 2p pair still sums to
    # 2.0 electrons, which is the invariant that matters. Do not "correct" 0.8 : 1.2 towards 0.667 : 1.333.
    #
    # For the MAGNITUDE: the total cross section peaks here at 2.53e-16 cm^2 (2.527e+08 barn) near 70 eV,
    # which is the right order for neutral carbon. Compare against Kim & Desclaux, Phys. Rev. A 66, 012708
    # (2002), who applied BEB to exactly C, N and O, and against the measurements of Brook, Harrison &
    # Smith, J. Phys. B 11, 3115 (1978), still the recommended data for neutral carbon. No literature number
    # is asserted here because none was verified at first hand this session; the comparison is left explicit
    # rather than baked into a tolerance that would then be believed.
    approx      = ImpactIonization.BEBmodel()
    multipleN   = 1
    iEnergies   = [20., 30., 50., 70., 100., 150., 200., 300., 500., 1000.]   ## C I ionization potential 11.26 eV
    shells      = [Shell("1s"), Shell("2s"), Shell("2p")]
    selection   = ShellSelection(true, shells, Int64[])
    configs     = [Configuration("1s^2 2s^2 2p^2")]
    name        = "EII cross section for the open-shell target C I (BEBmodel)."
    nucModel    = Nuclear.Model(6.0)
    eiiSettings = ImpactIonization.Settings(approx, multipleN, iEnergies, false, true, selection)
    comp        = Empirical.Computation(name, nucModel, grid, configs, eiiSettings)
    #
    ## The occupation sum rule, printed so that it is checked rather than assumed
    wa    = Atomic.Computation(Atomic.Computation(), name="C I basis", grid=grid, nuclearModel=nucModel, configs=configs)
    basis = perform(wa; output=true)["multiplet:"].levels[1].basis
    println("\n>> C I basis holds $(length(basis.csfs)) CSFs and $(basis.NoElectrons) electrons.")
    qList = [Basics.computeMeanSubshellOccupation(sh, basis)  for sh in basis.subshells]
    for  (k, sh) in enumerate(basis.subshells)    println("     mean occupation of $sh  = $(qList[k])")    end
    println(">> The occupations sum to $(sum(qList)) and must equal $(basis.NoElectrons) electrons exactly.")
    #
    perform(comp; output=true)
    #
    setDefaults("print summary: close", "")
    #
end
 
