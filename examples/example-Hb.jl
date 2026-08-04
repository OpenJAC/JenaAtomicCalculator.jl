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
elseif true
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
end
 
