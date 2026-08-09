
"""
`module  ForPedestrians`  
	... a submodel of JAC that comprises all functions to support a "pedestrian approach" to atomic computations.
	    The idea of this approach is to provide a set of very simple functions, i.e. of functions with simplified argument 
	    lists, in order to encourage useful computations. Other information, which is needed for such computation,
        are provided as defaults. This approach also provides a number of explanations and hints, how more advanced 
        computations can be carried out by means of the JAC toolbox.
        
        We hope and expect to add further functionality to this "pedestrian approach" but without giving up the basic
        idea to KISS: Keep Input Short and Simple.
"""
module ForPedestrians


using  Printf, ..AngularMomentum, ..Atomic, ..AutoIonization, ..Basics, ..Cascade, ..Defaults, ..DielectronicRecombination,
               ..Empirical, ..ImpactIonization, ..ManyElectron,  ..Nuclear, ..PhotoEmission, ..PhotoIonization,
               ..PhotoRecombination, ..Radial

export computeBranchingFractions,  computeChargeStateDistribution,  computeCrossSections,  computeForPedestrians,  computeLevelEnergies,
       computeLifetimes,  computeResonanceStrength,  computeTransitionRates,
       displayCouplings,  displaySpectrum,  estimateCrossSections


       
"""
`ForPedestrians.computeCrossSections(theme::Basics.ForPhotoIonization, initialConfigs::Array{Configuration,1},
                                     finalConfigs::Array{Configuration,1};
                                     grid::Radial.Grid=Radial.Grid(), asfSettings::AsfSettings=AsfSettings(),
                                     printout::Bool=false)` 
    ... computes the photoionization cross sections for all levels that are defined by the given initial and final 
        configurations. The default settings of the grid and asfSettings are used but can be overwritten on demand.
        The results are printed to screen but nothing is returned otherwise.
        
        Simplified call:   setDefaults("nuclear: charge", 10.0)
                           initialConfigs = [Configuration("[Ne]")]
                           finalConfigs   = [Configuration("[He] 2s 2p^6"), Configuration("[He] 2s^2 2p^5")]
                           computeCrossSections(Basics.ForPhotoIonization(), initialConfigs, finalConfigs)
"""
function computeCrossSections(theme::Basics.ForPhotoIonization, initialConfigs::Array{Configuration,1},
                              finalConfigs::Array{Configuration,1};
                              grid::Radial.Grid=Radial.Grid(), asfSettings::AsfSettings=AsfSettings(),
                              printout::Bool=false)
    Basics.checkConfigurations(Basics.NumberOfElectrons(), initialConfigs, finalConfigs)
    # Collect explanations
    configs = copy(initialConfigs);    append!(configs, finalConfigs)
    Basics.displayConfigurations(stdout, configs, details = "photoionization computations")
        
    sa =    "\n* Compute the photoionization cross sections between all levels from the given initial and final configurations above; " *
            "the following assumptions/simplifications are made: " *
            "\n    + All photoionization cross sections are based on the electric-dipole (E1) approximation only. " *
            "\n    + Use the optional argument  printout = true  to generate intermediate printout." *
            "\n    + Use the optional argument  grid = Radial.Grid(...)  to refine the radial grid, if needed." *
            "\n    + For more elaborate computations, make use of perform(comp::Atomic.Computation) " *
            "\n    + Call ? Atomic.Computation for further details. " *
            "\n    + Call ? setDefaults ... to define user-specified units for the computations. \n"
    println(sa)
    
    # Determine a useful grid
    if      grid.NoPoints == 0    currentGrid = Radial.Grid(Radial.Grid(false), rnt = 1.0e-5, h = 5.0e-2, hp = 2.0e-2, rbox = 15.0)
    else                          currentGrid = grid
    end
    
    photoSettings = PhotoIonization.Settings(PhotoIonization.Settings(), printBefore =true, photonEnergies = [2.0^(i-1) for i=1:12],
                                             lValues = [0, 1, 2, 3])
     
    # Specify the atomic computations
    function atomic_code()    
        Defaults.setDefaults("standard grid", currentGrid)
        Z = Defaults.getDefaults("nuclear: charge")
        
        comp = Atomic.Computation(Atomic.Computation(), name="Photoionization cross sections",  
                                  grid=currentGrid, nuclearModel=Nuclear.Model(Z);
                                  initialConfigs = initialConfigs, finalConfigs = finalConfigs, 
                                  processSettings = photoSettings ); 
        results = perform(comp, output=true)
        return( results )
    end
    
    # Print or suppress the standard output
    if    printout  atomic_code()
    else  
          results = redirect_stdout(devnull) do   
                       atomic_code()  end
          photoLines = results["photoionization lines:"]
          PhotoIonization.displayResults(stdout, photoLines, photoSettings)
    end
        
    return( nothing )
end


#################################################################################################################################
#################################################################################################################################


"""
`ForPedestrians.computeCrossSections(theme::Basics.ForPhotoRecombination, initialConfigs::Array{Configuration,1};
                                     grid::Radial.Grid=Radial.Grid(), asfSettings::AsfSettings=AsfSettings(),
                                     printout::Bool=false)`
    ... computes the photorecombination (radiative capture) cross sections for all levels defined by the given
        initial-ion configurations. The shells into which the electron is captured are taken from theme.intoShells.
        Final-state configurations are generated automatically. The default settings of the grid and asfSettings
        are used but can be overwritten on demand.  Results are printed to screen but nothing is returned otherwise.

        Simplified call:   setDefaults("nuclear: charge", 10.0)
                           setDefaults("unit: energy",         "eV")
                           setDefaults("unit: cross section",  "barn")
                           initialConfigs = [Configuration("1s")]
                           intoShells     = [Shell("1s")]
                           computeCrossSections(Basics.ForPhotoRecombination(intoShells), initialConfigs)
"""
function computeCrossSections(theme::Basics.ForPhotoRecombination, initialConfigs::Array{Configuration,1};
                              grid::Radial.Grid=Radial.Grid(), asfSettings::AsfSettings=AsfSettings(),
                              printout::Bool=false)
    Basics.checkConfigurations(Basics.NumberOfElectrons(), initialConfigs)
    Basics.displayConfigurations(stdout, initialConfigs, details = "photorecombination computations")

    sa =    "\n* Compute the photorecombination cross sections between all levels from the given initial " *
            "and the auto-generated final configurations above; " *
            "the following assumptions/simplifications are made: " *
            "\n    + All photorecombination cross sections are based on the electric-dipole (E1) approximation only. " *
            "\n    + Continuum orbitals are generated with the B-spline Galerkin method and pure-sine normalization. " *
            "\n    + The electron is captured into shells: $(theme.intoShells). " *
            "\n    + Use the optional argument  printout = true  to generate intermediate printout." *
            "\n    + Use the optional argument  grid = Radial.Grid(...)  to refine the radial grid, if needed." *
            "\n    + For more elaborate computations, make use of perform(comp::Atomic.Computation) " *
            "\n    + Call ? Atomic.Computation for further details. " *
            "\n    + Call ? setDefaults ... to define user-specified units for the computations. \n"
    println(sa)

    # Determine a useful grid
    if      grid.NoPoints == 0    currentGrid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)
    else                          currentGrid = grid
    end

    phSettings = PhotoRecombination.Settings(PhotoRecombination.Settings(), multipoles = [E1],
                                             gauges = [UseCoulomb, UseBabushkin],
                                             electronEnergies = [1., 5., 10., 50., 100., 500.],
                                             calcTotalCs = true, printBefore = true, maxKappa = 2)

    # Specify the atomic computations
    function atomic_code()
        Defaults.setDefaults("standard grid",         currentGrid)
        Defaults.setDefaults("method: continuum, Galerkin")
        Defaults.setDefaults("method: normalization, pure sine")
        Z            = Defaults.getDefaults("nuclear: charge")
        finalConfigs = Basics.generateConfigurations(theme, initialConfigs)

        comp = Atomic.Computation(Atomic.Computation(), name="Photorecombination cross sections",
                                  grid=currentGrid, nuclearModel=Nuclear.Model(Z);
                                  initialConfigs = initialConfigs, finalConfigs = finalConfigs,
                                  processSettings = phSettings );
        results = perform(comp, output=true)
        return( results )
    end

    # Print or suppress the standard output
    if    printout  atomic_code()
    else
          results  = redirect_stdout(devnull) do
                         atomic_code()  end
          phLines  = results["photo recombination lines:"]
          PhotoRecombination.displayResults(stdout, phLines, phSettings)
    end

    return( nothing )
end


#################################################################################################################################
#################################################################################################################################


"""
`ForPedestrians.computeForPedestrians()`
    ... prints an overview of all pedestrian functions available in this module, together with brief usage hints.
"""
function computeForPedestrians()

    sa =    "\n* ForPedestrians — a pedestrian approach to atomic computations with JAC." *
            "\n  Set the nuclear charge first:   setDefaults(\"nuclear: charge\", Z::Float64)" *
            "\n  Set output units as needed:     setDefaults(\"unit: energy\", \"eV\")   etc." *
            "\n" *
            "\n  Level energies and fine-structure:" *
            "\n    computeLevelEnergies(ForGivenConfigs(), configs)" *
            "\n      ... level energies and J^P for all levels of the given configurations (Dirac-Coulomb)." *
            "\n    computeLevelEnergies(ForIsoelectronicSequence(), Zvalues, configs)" *
            "\n      ... configuration-averaged energies for each Z in Zvalues; one row per Z, one column per config." *
            "\n    displayCouplings(FineStructure(), configs)" *
            "\n      ... jj-coupling fine-structure levels (J values and multiplicities)." *
            "\n    displayCouplings(FineStructureLS(), configs)" *
            "\n      ... LS-coupling term symbols (^{2S+1}L) for the given configurations." *
            "\n" *
            "\n  Transition rates and lifetimes:" *
            "\n    computeTransitionRates(ForPhotoEmission(), initialConfigs, finalConfigs)" *
            "\n      ... E1 radiative rates and oscillator strengths; initial and final configs have the same N." *
            "\n    computeTransitionRates(ForAutoIonization(), initialConfigs, finalConfigs)" *
            "\n      ... Auger (autoionization) rates; finalConfigs must have N-1 electrons." *
            "\n    computeBranchingFractions(ForPhotoEmission(), initialConfigs, finalConfigs)" *
            "\n      ... E1 radiative branching fractions BF(i->f) [%]; Coulomb and Babushkin gauge." *
            "\n    computeBranchingFractions(ForAutoIonization(), initialConfigs, finalConfigs)" *
            "\n      ... Auger branching fractions BF(i->f) [%]; Coulomb interaction only." *
            "\n    displaySpectrum(ForPhotoEmission(), initialConfigs, finalConfigs; plotfile=...)" *
            "\n      ... photon emission spectrum as ASCII bar chart; optional Plots.jl figure." *
            "\n    displaySpectrum(ForAutoIonization(), initialConfigs, finalConfigs; plotfile=...)" *
            "\n      ... Auger electron spectrum as ASCII bar chart; optional Plots.jl figure." *
            "\n    computeLifetimes(ForPhotoEmission(), configs)" *
            "\n      ... total radiative lifetime of levels with inner-shell holes (E1 only)." *
            "\n    computeLifetimes(ForAutoIonization(), configs)" *
            "\n      ... total Auger lifetime; final configs (N-1 electrons) are generated automatically." *
            "\n" *
            "\n  Photoionization and photorecombination cross sections:" *
            "\n    computeCrossSections(ForPhotoIonization(), initialConfigs, finalConfigs)" *
            "\n      ... photoionization cross sections (E1); finalConfigs have N-1 electrons." *
            "\n    computeCrossSections(ForPhotoRecombination(intoShells), initialConfigs)" *
            "\n      ... photorecombination cross sections (E1); final configs generated from intoShells." *
            "\n" *
            "\n  Dielectronic recombination resonance strengths:" *
            "\n    computeResonanceStrength(ForDielectronicRecombination(fromShells,toShells,intoShells,decayShells), initialConfigs)" *
            "\n      ... DR resonance strengths [cm^2 eV]; intermediate and final configs auto-generated." *
            "\n      ... Use  setDefaults(\"unit: strength\", \"cm^2 eV\")  for the output unit." *
            "\n" *
            "\n  Empirical cross section estimates:" *
            "\n    estimateCrossSections(ForImpactIonization(), initialConfigs)" *
            "\n      ... relativistic BEB electron-impact ionization cross sections for all shells." *
            "\n      ... Optional:  electronEnergies = [...]  to set the impact-energy grid [eV]." *
            "\n" *
            "\n  Decay cascades and charge state distributions:" *
            "\n    computeChargeStateDistribution(ForStepwiseDecay(N), initialConfigs)" *
            "\n      ... charge state distribution after stepwise radiative+Auger cascade." *
            "\n      ... N = max. electron-loss steps; use N=3 for K-hole in Ne, N=5 for Ar." *
            "\n" *
            "\n  For more elaborate computations:  perform(comp::Atomic.Computation)" *
            "\n    Call ? Atomic.Computation   for the full interface." *
            "\n    Call ? Cascade.Computation  for multi-step radiative and Auger cascades." *
            "\n    Call ? Plasma.Computation   for plasma-shift and average-atom computations." *
            "\n    Call ? RacahAlgebra.RacahExpression  for symbolic Racah-algebra evaluations." *
            "\n" *
            "\n  Reference:  S. Fritzsche, Comp. Phys. Commun. 240, 1-14 (2019)."
    println(sa)
    return( nothing )
end



#################################################################################################################################
#################################################################################################################################


"""
`ForPedestrians.computeLevelEnergies(theme::Basics.ForGivenConfigs, configs::Array{Configuration,1};
                                     grid::Radial.Grid=Radial.Grid(true), asfSettings::AsfSettings=AsfSettings(),
                                     printout::Bool=false)` 
    ... computes the level energies and leading configurations for all levels that are defined by the given configurations.
        The default values of the grid and asfSettings are used but can be overwritten on demand.
        The results are printed to screen but nothing is returned otherwise.
        
        Simplified call:   setDefaults("nuclear: charge", 10.0)
                           configs = [Configuration("[Ne]"), Configuration("[He] 2s^2 2p^5 3s")]
                           computeLevelEnergies(Basics.ForGivenConfigs(), configs)
"""
function computeLevelEnergies(theme::Basics.ForGivenConfigs, configs::Array{Configuration,1};
                              grid::Radial.Grid=Radial.Grid(true), asfSettings::AsfSettings=AsfSettings(),
                              printout::Bool=false)
    # Collect explanations
    Basics.displayConfigurations(stdout, configs, details = "energy level computations")

    sa =    "\n* Compute the level energies for all levels of the configurations above; " *
            "the following assumptions/simplifications are made: " *
            "\n    + All level energies are based on the Dirac-Coulomb Hamiltonian only. " *
            "\n    + Configurations with different numbers of electrons are handled in separate computations." *
            "\n    + Use the optional argument  printout = true  to generate intermediate printout." *
            "\n    + Use the optional argument  grid = Radial.Grid(...)  to refine the radial grid, if needed." *
            "\n    + For more elaborate computations, make use of perform(comp::Atomic.Computation) " *
            "\n    + Call ? Atomic.Computation for further details. " *
            "\n    + Call ? setDefaults ... to define user-specified units for the computations. \n"
    println(sa)

    # Run a separate Atomic.Computation for each group of configurations with the same number of electrons
    function atomic_code(confs)
        Defaults.setDefaults("standard grid", grid)
        Z    = Defaults.getDefaults("nuclear: charge")
        comp = Atomic.Computation(Atomic.Computation(), name="Level energies",
                                  grid=grid, nuclearModel=Nuclear.Model(Z), configs=confs);
        results = perform(comp, output=true)
        return( results )
    end

    Ns         = sort( unique( Basics.extractFromConfigurations(Basics.NumberOfElectrons(), configs) ) )
    multiplets = ManyElectron.Multiplet[]
    for  N  in  Ns
        confs = Basics.extractConfigurations(Basics.ByNumber([N]), configs)
        if    printout  atomic_code(confs)
        else
              results = redirect_stdout(devnull) do
                            atomic_code(confs)  end
              push!(multiplets, results["multiplet:"])
        end
    end

    if  !printout   Basics.displayLevels(stdout, multiplets; N=200)   end

    return( nothing )
end


"""
`ForPedestrians.computeLevelEnergies(theme::Basics.ForIsoelectronicSequence,
                                     Zvalues::Array{Float64,1},
                                     configs::Array{Configuration,1};
                                     grid::Radial.Grid=Radial.Grid(true),
                                     asfSettings::AsfSettings=AsfSettings(),
                                     printout::Bool=false)`
    ... computes the configuration-averaged level energy for each configuration in configs
        and for each nuclear charge in Zvalues, tracing the given configurations along an
        isoelectronic sequence.  For each (Z, config) pair a separate SCF+CI computation
        is performed and the (2J+1)-weighted mean energy is returned.
        Results are printed as a compact table; nothing is returned otherwise.

        Simplified call:   Zvalues = collect(4.0:1.0:10.0)
                           configs = [Configuration("1s^2"), Configuration("1s^2 2s")]
                           computeLevelEnergies(Basics.ForIsoelectronicSequence(), Zvalues, configs)
"""
function computeLevelEnergies(theme::Basics.ForIsoelectronicSequence,
                               Zvalues::Array{Float64,1},
                               configs::Array{Configuration,1};
                               grid::Radial.Grid=Radial.Grid(true),
                               asfSettings::AsfSettings=AsfSettings(),
                               printout::Bool=false,
                               plotfile::String="")
    Basics.checkConfigurations(Basics.NumberOfElectrons(), configs)
    sa =    "\n* Compute configuration-averaged level energies along an isoelectronic sequence; " *
            "the following assumptions/simplifications are made: " *
            "\n    + One separate SCF+CI computation per (Z, configuration) pair." *
            "\n    + Configuration-averaged energy = (2J+1)-weighted mean over all levels." *
            "\n    + All level energies are based on the Dirac-Coulomb Hamiltonian only. " *
            "\n    + Use the optional argument  printout = true  to generate intermediate printout." *
            "\n    + For more elaborate computations, make use of perform(comp::Atomic.Computation) \n"
    println(sa)

    # Compute the configuration-averaged energy for a single (Z, config) pair
    function avg_energy(Z::Float64, conf::Configuration)
        function atomic_code()
            Defaults.setDefaults("standard grid", grid)
            comp = Atomic.Computation(Atomic.Computation(), name="Isoelectronic level energies",
                                      grid=grid, nuclearModel=Nuclear.Model(Z), configs=[conf])
            return( perform(comp, output=true) )
        end
        results  = redirect_stdout(devnull) do
                       atomic_code()  end
        multiplet = results["multiplet:"]
        # (2J+1)-weighted average energy
        weight = 0.;   Etotal = 0.
        for  lev  in  multiplet.levels
            w = Float64(Basics.twice(lev.J) + 1)
            weight  += w
            Etotal  += w * lev.energy
        end
        return( weight > 0. ? Etotal/weight : 0. )
    end

    # Build result matrix: rows = Z values, cols = configurations
    nZ   = length(Zvalues)
    nC   = length(configs)
    Eavg = zeros(Float64, nZ, nC)
    for  (iz, Z)    in  enumerate(Zvalues)
        for  (ic, conf)  in  enumerate(configs)
            Eavg[iz, ic] = avg_energy(Z, conf)
        end
    end

    # Display as a compact table using plain string formatting (TableStrings not imported)
    eunit  = Defaults.getDefaults("unit: energy")
    col_w  = 22
    nx     = 12 + col_w * nC
    hline  = "  " * "-"^nx
    # Helper: center a string in a field of given width
    ctr = (s, w) -> begin n = length(s); n >= w && return s;
                    l = div(w-n,2); " "^l * s * " "^(w-n-l) end
    println("  Configuration-averaged energies along an isoelectronic sequence:\n")
    println(hline)
    sa = "  " * @sprintf("%10s  ", "Z")
    for  conf  in  configs
        sa = sa * ctr(string(conf) * " [" * eunit * "]", col_w)
    end
    println(sa);   println(hline)
    for  (iz, Z)  in  enumerate(Zvalues)
        sb = "  " * @sprintf("%10.1f  ", Z)
        for  ic  in  1:nC
            sb = sb * ctr(@sprintf("%.6e", Defaults.convertUnits("energy: from atomic", Eavg[iz,ic])), col_w)
        end
        println(sb)
    end
    println(hline)

    # Optional: generate and save a line plot if Plots.jl is loaded and a filename is given
    if  !isempty(plotfile)
        if  isdefined(Main, :Plots)
            Plots = Main.Plots
            labels = reshape([string(c) for c in configs], 1, :)
            p = Plots.plot(Zvalues, Eavg,
                           xlabel = "Nuclear charge Z",
                           ylabel = "Config.-averaged energy [" * eunit * "]",
                           label  = labels,
                           marker = :circle, linewidth = 2,
                           title  = "Isoelectronic sequence")
            Plots.savefig(p, plotfile)
            println("  Plot saved to: " * plotfile)
        else
            println("  Plotting skipped — load Plots.jl first:  using Plots")
        end
    end

    return( nothing )
end


#################################################################################################################################
#################################################################################################################################


"""
`ForPedestrians.computeLifetimes(theme::Basics.ForPhotoEmission, configs::Array{Configuration,1};
                                 grid::Radial.Grid=Radial.Grid(true), asfSettings::AsfSettings=AsfSettings(),
                                 printout::Bool=false)` 
    ... computes the radiative lifetimes of all levels that are defined by the given configurations.
        The default values of the grid and asfSettings are used but can be overwritten on demand.
        The results are printed to screen but nothing is returned otherwise.
        
        Simplified call:   setDefaults("nuclear: charge", 10.0)
                           configs = [Configuration("1s 2s^2 2p^6")]
                           computeLifetimes(Basics.ForPhotoEmission(), configs)
"""
function computeLifetimes(theme::Basics.ForPhotoEmission, configs::Array{Configuration,1};
                          grid::Radial.Grid=Radial.Grid(true), asfSettings::AsfSettings=AsfSettings(),
                          printout::Bool=false)
    Basics.checkConfigurations(Basics.NumberOfElectrons(), configs)
    # Collect explanations
    Basics.displayConfigurations(stdout, configs, details = "radiative lifetime computations")
        
    sa =    "\n* Compute the radiative lifetimes for all levels of the configurations above; " *
            "the following assumptions/simplifications are made: " *
            "\n    + Lifetimes are calculated for levels with inner-shell holes (at least, one sub-valence hole) only." *
            "\n    + All radiative lifetimes are based on the electric-dipole (E1) approximation only. " *
            "\n    + Use the optional argument  printout = true  to generate intermediate printout." *
            "\n    + Use the optional argument  grid = Radial.Grid(...)  to refine the radial grid, if needed." *
            "\n    + For more elaborate computations, make use of perform(comp::Atomic.Computation) " *
            "\n    + Call ? Atomic.Computation for further details. " *
            "\n    + Call ? setDefaults ... to define user-specified units for the computations. \n"
    println(sa)
    
    photoSettings = PhotoEmission.Settings(PhotoEmission.Settings(), multipoles=[E1], gauges=[UseCoulomb, UseBabushkin], 
                                           printBefore=true)
     
    # Specify the atomic computations
    function atomic_code()        
        Defaults.setDefaults("standard grid", grid)
        Z             = Defaults.getDefaults("nuclear: charge")
        finalConfigs  = Basics.generateConfigurations(Basics.ForPhotoEmission(), configs)
       
        comp = Atomic.Computation(Atomic.Computation(), name="PhotoEmission (radiative) lifetimes",  
                                  grid=grid, nuclearModel=Nuclear.Model(Z);
                                  initialConfigs = configs, finalConfigs = finalConfigs, 
                                  processSettings = photoSettings ); 
        results = perform(comp, output=true)
        return( results )
    end
    
    # Print or suppress the standard output
    if    printout  atomic_code()
    else  
          results = redirect_stdout(devnull) do   
                       atomic_code()  end
          lines = results["radiative lines:"]
          PhotoEmission.displayLifetimes(stdout, lines, photoSettings)
    end
    
    return( nothing )
end 


"""
`ForPedestrians.computeLifetimes(theme::Basics.ForAutoIonization, configs::Array{Configuration,1};
                                 grid::Radial.Grid=Radial.Grid(), asfSettings::AsfSettings=AsfSettings(),
                                 printout::Bool=false)` 
    ... computes the non-radiative (Auger) lifetimes of all levels that are defined by the given configurations.
        The default values of the grid and asfSettings are used but can be overwritten on demand.
        The results are printed to screen but nothing is returned otherwise.
        
        Simplified call:   setDefaults("nuclear: charge", 10.0)
                           configs = [Configuration("1s 2s^2 2p^6")]
                           computeLifetimes(Basics.ForAutoIonization(), configs)
"""
function computeLifetimes(theme::Basics.ForAutoIonization, configs::Array{Configuration,1};
                          grid::Radial.Grid=Radial.Grid(), asfSettings::AsfSettings=AsfSettings(),
                          printout::Bool=false)
    Basics.checkConfigurations(Basics.NumberOfElectrons(), configs)
    # Collect explanations
    Basics.displayConfigurations(stdout, configs, details = "Auger lifetime computations")
        
    sa =    "\n* Compute the non-radiative (Auger) lifetimes for all levels of the configurations above; " *
            "the following assumptions/simplifications are made: " *
            "\n    + All Auger lifetimes are based on the instantaneous Coulomb interaction. " *
            "\n    + Continuum orbitals are generated with the B-spline Galerkin method and pure-sine normalization. " *
            "\n    + Use the optional argument  printout = true  to generate intermediate printout." *
            "\n    + Use the optional argument  grid = Radial.Grid(...)  to refine the radial grid, if needed." *
            "\n    + For more elaborate computations, make use of perform(comp::Atomic.Computation) " *
            "\n    + Call ? Atomic.Computation for further details. " *
            "\n    + Call ? setDefaults ... to define user-specified units for the computations. \n"
    println(sa)
    
    # Determine a useful grid
    if      grid.NoPoints == 0    currentGrid = Radial.Grid(Radial.Grid(false), rnt = 1.0e-5, h = 5.0e-2, hp = 2.0e-2, rbox = 15.0)
    else                          currentGrid = grid
    end
    
    # Specify the atomic computations
    function atomic_code()
        Defaults.setDefaults("standard grid",         currentGrid)
        Defaults.setDefaults("method: continuum, Galerkin")
        Defaults.setDefaults("method: normalization, pure sine")
        Z             = Defaults.getDefaults("nuclear: charge")
        finalConfigs  = Basics.generateConfigurations(Basics.ForAutoIonization(), configs)
        augerSettings = AutoIonization.Settings()

        comp = Atomic.Computation(Atomic.Computation(), name="Non-radiative (Auger) lifetimes",
                                  grid=currentGrid, nuclearModel=Nuclear.Model(Z);
                                  initialConfigs = configs, finalConfigs = finalConfigs,
                                  processSettings = augerSettings );
        results = perform(comp, output=true)
        return( results )
    end

    # Print or suppress the standard output
    if    printout  atomic_code()
    else
          results = redirect_stdout(devnull) do
                        atomic_code()  end
          lines = results["AutoIonization lines:"]
          AutoIonization.displayLifetimes(stdout, lines)
    end
    
    return( nothing )
end 

#################################################################################################################################
#################################################################################################################################


"""
`ForPedestrians.computeResonanceStrength(theme::Basics.ForDielectronicRecombination, initialConfigs::Array{Configuration,1};
                                         grid::Radial.Grid=Radial.Grid(), asfSettings::AsfSettings=AsfSettings(),
                                         printout::Bool=false)` 
    ... computes the dielectronic recombination resonance strength of all levels of the given configurations.
        The default values of the grid and asfSettings are used but can be overwritten on demand.
        The results are printed to screen but nothing is returned otherwise.
        
        Simplified call:   setDefaults("nuclear: charge", 10.0)
                           setDefaults("unit: strength", "cm^2 eV")   
                           initialConfigs = [Configuration("1s^2 2s")]
                           fromShells     = [Shell("2s")]
                           toShells       = [Shell("2p")]
                           intoShells     = Basics.generateShellList( 7,  7, 3)
                           decayShells    = Basics.generateShellList( 2,  4, 3)
                           theme          = Basics.ForDielectronicRecombination(fromShells, toShells, intoShells, decayShells)
                           computeResonanceStrength(theme, initialConfigs, printout=false)
"""
function computeResonanceStrength(theme::Basics.ForDielectronicRecombination, initialConfigs::Array{Configuration,1};
                                  grid::Radial.Grid=Radial.Grid(), asfSettings::AsfSettings=AsfSettings(), printout::Bool=false)
    Basics.checkConfigurations(Basics.NumberOfElectrons(), initialConfigs)
    Basics.displayConfigurations(stdout, initialConfigs, details = "DR resonance strength (initial configurations)")
        
    # Collect explanations
    sa =    "\n* computes the dielectronic recombination resonance strength of all levels from the configurations above; " *
            "the following assumptions/simplifications are made: " *
            "\n    + The intermediate (doubly-excited) and final-state configurations are generated automatically." *
            "\n    + The intermediate configurations include excitations fromShells --> toShells + the capture intoShells." *
            "\n    + The final-state configurations include the de-excitation toShells, intoshells --> decayShells." *
            "\n    + These two lists of configurations can be controlled by generating proper shell lists." *
            "\n    + Call ? Basics.generateShellList  to understand how useful shell lists can be generated." *
            "\n    + Use the optional argument  printout = true/false  to generate intermediate printout." *
            "\n    + Use the optional argument  grid = Radial.Grid(...)  to refine the radial grid, if needed." *
            "\n    + For more elaborate computations, make use of perform(comp::Atomic.Computation) " *
            "\n    + Call ? Atomic.Computation for further details. " *
            "\n    + Call ? setDefaults ... to define user-specified units for the computations. \n"
    println(sa)
    
    # Generate the intermediate and final-state configurations
    (intermediateConfs, finalConfs) = Basics.generateConfigurations(theme, initialConfigs)
    Basics.displayConfigurations(stdout, intermediateConfs, details = "dielectronic capture")
    Basics.displayConfigurations(stdout, finalConfs,        details = "(radiative) stabilization")
    
    # Determine a useful grid
    if      grid.NoPoints == 0    currentGrid = Radial.Grid(Radial.Grid(false), rnt = 1.0e-5, h = 5.0e-2, hp = 2.0e-2, rbox = 20.0)
    else                          currentGrid = grid
    end
    
    # Specify physical data
    Z           = Defaults.getDefaults("nuclear: charge")
    drSettings  = DielectronicRecombination.Settings(DielectronicRecombination.Settings(), multipoles = [E1], gauges = [UseCoulomb, UseBabushkin],
                                                     printBefore = true, electronEnergyShift = 0.)
                                            
    # Specify the atomic computations
    function atomic_code()
        comp        = Atomic.Computation(Atomic.Computation(), name="Dielectronic recombination resonance strength computations", 
                                         grid=currentGrid, nuclearModel=Nuclear.Model(Z), 
                                         initialConfigs = initialConfigs, intermediateConfigs = intermediateConfs,  
                                         finalConfigs = finalConfs, processSettings = drSettings )

        results     = perform(comp, output=true)
        return( results )
    end
    
    # Print or suppress the standard output
    if    printout  atomic_code()
    else  
          results = redirect_stdout(devnull) do   
                       atomic_code()  end
          captureLines, photonLines = results["dielectronic recombination lines:"]
          ## The CaptureLine's already carry the resonance strength, so no separate aggregation step is needed;
          ## the per-final-level satellite table is printed only if settings.calcPhotonSpectrum is set.
          DielectronicRecombination.displayResults(stdout, captureLines, photonLines, drSettings)
    end
    
    return( nothing )
end 


#################################################################################################################################
#################################################################################################################################


"""
`ForPedestrians.computeTransitionRates(theme::Basics.ForAutoIonization, initialConfigs::Array{Configuration,1},
                                       finalConfigs::Array{Configuration,1};
                                       grid::Radial.Grid=Radial.Grid(), asfSettings::AsfSettings=AsfSettings(),
                                       printout::Bool=false)` 
    ... computes the autoionization (Auger) rates for all levels that are defined by the given initial and final 
        configurations. The final-state configurations must have one electron less than the initial-state configurations.
        The default settings of the grid and asfSettings are used but can be overwritten on demand.
        The results are printed to screen but nothing is returned otherwise.
        
        Simplified call:   setDefaults("nuclear: charge", 10.0)
                           initialConfigs = [Configuration("1s 2s^2 2p^6")]
                           finalConfigs   = [Configuration("[He] 2s^0 2p^6"), Configuration("[He] 2s 2p^5"),
                                             Configuration("[He] 2s^2 2p^4")]
                           computeTransitionRates(Basics.ForAutoIonization(), initialConfigs, finalConfigs)
"""
function computeTransitionRates(theme::Basics.ForAutoIonization, initialConfigs::Array{Configuration,1},
                                finalConfigs::Array{Configuration,1};
                                grid::Radial.Grid=Radial.Grid(), asfSettings::AsfSettings=AsfSettings(),
                                printout::Bool=false)
    Basics.checkConfigurations(Basics.NumberOfElectrons(), initialConfigs, finalConfigs)
    # Collect explanations
    configs = copy(initialConfigs);    append!(configs, finalConfigs)
    Basics.displayConfigurations(stdout, configs, details = "autoionization computations")
        
    sa =    "\n* Compute the autoionization (Auger) rates between all levels from the given initial and final configurations above; " *
            "the following assumptions/simplifications are made: " *
            "\n    + All autoionization (Auger) rates are based on the instantaneous Coulomb interaction only. " *
            "\n    + Continuum orbitals are generated with the B-spline Galerkin method and pure-sine normalization. " *
            "\n    + Use the optional argument  printout = true  to generate intermediate printout." *
            "\n    + Use the optional argument  grid = Radial.Grid(...)  to refine the radial grid, if needed." *
            "\n    + For more elaborate computations, make use of perform(comp::Atomic.Computation) " *
            "\n    + Call ? Atomic.Computation for further details. " *
            "\n    + Call ? setDefaults ... to define user-specified units for the computations. \n"
    println(sa)
    
    # Determine a useful grid
    if      grid.NoPoints == 0    currentGrid = Radial.Grid(Radial.Grid(false), rnt = 1.0e-5, h = 5.0e-2, hp = 2.0e-2, rbox = 15.0)
    else                          currentGrid = grid
    end
        
    augerSettings = AutoIonization.Settings()    ## AutoIonization.Settings(), printout=true, operator=CoulombInteraction())
    
    # Specify the atomic computations
    function atomic_code()
        Defaults.setDefaults("standard grid",         currentGrid)
        Defaults.setDefaults("method: continuum, Galerkin")
        Defaults.setDefaults("method: normalization, pure sine")
        Z = Defaults.getDefaults("nuclear: charge")

        comp = Atomic.Computation(Atomic.Computation(), name="Autoionization (Auger) rates",
                                grid=currentGrid, nuclearModel=Nuclear.Model(Z);
                                initialConfigs = initialConfigs, finalConfigs = finalConfigs, 
                                processSettings = augerSettings ); 
        results = perform(comp, output=true)
        return( results )
    end
    
    # Print or suppress the standard output
    if    printout  atomic_code()
    else
          results = redirect_stdout(devnull) do
                        atomic_code()  end
          lines = results["AutoIonization lines:"]
          AutoIonization.displayRates(stdout, lines, augerSettings)
    end

    return( nothing )
end


"""
`ForPedestrians.computeTransitionRates(theme::Basics.ForPhotoEmission, initialConfigs::Array{Configuration,1},
                                       finalConfigs::Array{Configuration,1};
                                       grid::Radial.Grid=Radial.Grid(true), asfSettings::AsfSettings=AsfSettings(),
                                       printout::Bool=false)` 
    ... computes the photoemission rates and oscillator strengths for all levels that are defined by the given 
        initial and final configurations. Obviously, the initial- and final-state configurations must share the same 
        number of electrons. The default settings of the grid and asfSettings are used but can be overwritten on demand.
        The results are printed to screen but nothing is returned otherwise.
        
        Simplified call:   setDefaults("nuclear: charge", 10.0)
                           initialConfigs = [Configuration("1s 2s^2 2p^6")]
                           finalConfigs   = [Configuration("[He] 2s 2p^6"), Configuration("[He] 2s^2 2p^5")]
                           computeTransitionRates(Basics.ForPhotoEmission(), initialConfigs, finalConfigs)
"""
function computeTransitionRates(theme::Basics.ForPhotoEmission, initialConfigs::Array{Configuration,1},
                                finalConfigs::Array{Configuration,1};
                                grid::Radial.Grid=Radial.Grid(true), asfSettings::AsfSettings=AsfSettings(),
                                printout::Bool=false)
    configs = copy(initialConfigs);    append!(configs, finalConfigs)
    Basics.checkConfigurations(Basics.NumberOfElectrons(), configs)
    # Collect explanations
    Basics.displayConfigurations(stdout, configs, details = "photoemission computations")
        
    sa =    "\n* Compute the photoemission (radiative) rates and oscillator strengths between all levels " *
            "from the given initial and final configurations above; " *
            "\n  the following assumptions/simplifications are made: " *
            "\n    + All photoemission rates are based on the electric-dipole (E1) approximation only. " *
            "\n    + Use the optional argument  printout = true  to generate intermediate printout." *
            "\n    + Use the optional argument  grid = Radial.Grid(...)  to refine the radial grid, if needed." *
            "\n    + For more elaborate computations, make use of perform(comp::Atomic.Computation) " *
            "\n    + Call ? Atomic.Computation for further details. " *
            "\n    + ... \n"
    println(sa)
    
    photoSettings = PhotoEmission.Settings(PhotoEmission.Settings(), multipoles=[E1], gauges=[UseCoulomb, UseBabushkin], 
                                           printBefore=true)
    
    # Specify the atomic computations
    function atomic_code()    
        Defaults.setDefaults("standard grid", grid)
        Z = Defaults.getDefaults("nuclear: charge")
        
        comp = Atomic.Computation(Atomic.Computation(), name="Photoemission (radiative) rates and oscillator strengths",  
                                  grid=grid, nuclearModel=Nuclear.Model(Z);
                                  initialConfigs = initialConfigs, finalConfigs = finalConfigs, 
                                  processSettings = photoSettings ); 
        results = perform(comp, output=true)
        return( results )
    end
    
    # Print or suppress the standard output
    if    printout  atomic_code()
    else  
          results = redirect_stdout(devnull) do   
                       atomic_code()  end
          lines = results["radiative lines:"]
          PhotoEmission.displayRates(stdout, lines, photoSettings)
    end

    return( nothing )
end


#################################################################################################################################
#################################################################################################################################


"""
`ForPedestrians.computeBranchingFractions(theme::Basics.ForPhotoEmission, initialConfigs::Array{Configuration,1}, finalConfigs::Array{Configuration,1})`
    ... computes branching fractions for E1 radiative decay from all levels of initialConfigs to
        finalConfigs.  Both Coulomb and Babushkin gauge fractions are shown; all other rate quantities
        are suppressed.  The results are printed to screen but nothing is returned otherwise.

        Simplified call:   setDefaults("nuclear: charge", 10.0)
                           initialConfigs = [Configuration("1s 2s^2 2p^6")]
                           finalConfigs   = [Configuration("[He] 2s 2p^6"), Configuration("[He] 2s^2 2p^5")]
                           computeBranchingFractions(Basics.ForPhotoEmission(), initialConfigs, finalConfigs)
"""
function computeBranchingFractions(theme::Basics.ForPhotoEmission, initialConfigs::Array{Configuration,1},
                                   finalConfigs::Array{Configuration,1};
                                   grid::Radial.Grid=Radial.Grid(true), asfSettings::AsfSettings=AsfSettings(),
                                   printout::Bool=false)
    configs = copy(initialConfigs);   append!(configs, finalConfigs)
    Basics.checkConfigurations(Basics.NumberOfElectrons(), configs)
    Basics.displayConfigurations(stdout, configs, details = "radiative branching fraction computations")

    sa =    "\n* Compute E1 radiative branching fractions between all levels from the given configurations; " *
            "the following assumptions/simplifications are made: " *
            "\n    + BF(i->f) = A(i->f) / sum_f A(i->f);  Coulomb and Babushkin gauge shown separately." *
            "\n    + All rates are based on the electric-dipole (E1) approximation only. " *
            "\n    + Use the optional argument  printout = true  to generate intermediate printout." *
            "\n    + For more elaborate computations, make use of perform(comp::Atomic.Computation) \n"
    println(sa)

    photoSettings = PhotoEmission.Settings(PhotoEmission.Settings(), multipoles=[E1],
                                           gauges=[UseCoulomb, UseBabushkin], printBefore=true)
    function atomic_code()
        Defaults.setDefaults("standard grid", grid)
        Z = Defaults.getDefaults("nuclear: charge")
        comp = Atomic.Computation(Atomic.Computation(), name="Radiative branching fractions",
                                  grid=grid, nuclearModel=Nuclear.Model(Z);
                                  initialConfigs=initialConfigs, finalConfigs=finalConfigs,
                                  processSettings=photoSettings)
        return( perform(comp, output=true) )
    end

    if    printout  atomic_code()
    else
          results = redirect_stdout(devnull) do
                        atomic_code()  end
          lines = results["radiative lines:"]
          # Display branching fractions grouped by initial level
          nx    = 74
          hline = "  " * "-"^nx
          println(" ")
          println("  Radiative branching fractions (E1, Coulomb and Babushkin gauge):")
          println(" ")
          println(hline)
          println(@sprintf("  %-14s  %-18s  %12s  %10s  %10s",
                            "i-level-f", "i--J^P--f", "Energy [eV]", "BF Cou[%]", "BF Bab[%]"))
          println(hline)
          iIndices = unique(l.initialLevel.index for l in lines)
          for  iIdx  in  iIndices
              iLines = [l for l in lines if l.initialLevel.index == iIdx]
              sumCou = sum(l.photonRate.Coulomb   for l in iLines)
              sumBab = sum(l.photonRate.Babushkin for l in iLines)
              il  = iLines[1].initialLevel
              sym = LevelSymmetry(il.J, il.parity)
              println("  Level $(il.index)  ($(string(sym)))")
              for  l  in  iLines
                  bfC  = sumCou > 0. ? 100.0 * l.photonRate.Coulomb   / sumCou : 0.
                  bfB  = sumBab > 0. ? 100.0 * l.photonRate.Babushkin / sumBab : 0.
                  Eph  = Defaults.convertUnits("energy: from atomic", l.omega)
                  symI = string(LevelSymmetry(l.initialLevel.J, l.initialLevel.parity))
                  symF = string(LevelSymmetry(l.finalLevel.J,   l.finalLevel.parity))
                  lf   = @sprintf("%3i --> %3i", il.index, l.finalLevel.index)
                  jp   = @sprintf("%-5s --> %-7s", symI, symF)
                  println(@sprintf("  %-14s  %-18s  %12.4e  %10.4f  %10.4f", lf, jp, Eph, bfC, bfB))
              end
              println(@sprintf("  %-14s  %-18s  %12s  %10.2f  %10.2f", "", "", "Sum:", 100., 100.))
              println(hline)
          end
    end

    return( nothing )
end


"""
`ForPedestrians.computeBranchingFractions(theme::Basics.ForAutoIonization, initialConfigs::Array{Configuration,1}, finalConfigs::Array{Configuration,1})`
    ... computes branching fractions for Auger decay from all levels of initialConfigs to
        finalConfigs (N-1 electrons).  A single rate column is shown (Coulomb interaction).
        The results are printed to screen but nothing is returned otherwise.

        Simplified call:   setDefaults("nuclear: charge", 10.0)
                           initialConfigs = [Configuration("1s 2s^2 2p^6")]
                           finalConfigs   = [Configuration("[He] 2p^6"), Configuration("[He] 2s 2p^5"),
                                             Configuration("[He] 2s^2 2p^4")]
                           computeBranchingFractions(Basics.ForAutoIonization(), initialConfigs, finalConfigs)
"""
function computeBranchingFractions(theme::Basics.ForAutoIonization, initialConfigs::Array{Configuration,1},
                                   finalConfigs::Array{Configuration,1};
                                   grid::Radial.Grid=Radial.Grid(), asfSettings::AsfSettings=AsfSettings(),
                                   printout::Bool=false)
    Basics.checkConfigurations(Basics.NumberOfElectrons(), initialConfigs, finalConfigs)
    configs = copy(initialConfigs);   append!(configs, finalConfigs)
    Basics.displayConfigurations(stdout, configs, details = "Auger branching fraction computations")

    sa =    "\n* Compute Auger branching fractions between all levels from the given configurations; " *
            "the following assumptions/simplifications are made: " *
            "\n    + BF(i->f) = Gamma_A(i->f) / sum_f Gamma_A(i->f);  instantaneous Coulomb interaction." *
            "\n    + Continuum orbitals are generated with the B-spline Galerkin method and pure-sine normalization." *
            "\n    + Use the optional argument  printout = true  to generate intermediate printout." *
            "\n    + For more elaborate computations, make use of perform(comp::Atomic.Computation) \n"
    println(sa)

    if      grid.NoPoints == 0    currentGrid = Radial.Grid(Radial.Grid(false), rnt=1.0e-5, h=5.0e-2, hp=2.0e-2, rbox=15.0)
    else                          currentGrid = grid
    end
    augerSettings = AutoIonization.Settings()

    function atomic_code()
        Defaults.setDefaults("standard grid",         currentGrid)
        Defaults.setDefaults("method: continuum, Galerkin")
        Defaults.setDefaults("method: normalization, pure sine")
        Z = Defaults.getDefaults("nuclear: charge")
        comp = Atomic.Computation(Atomic.Computation(), name="Auger branching fractions",
                                  grid=currentGrid, nuclearModel=Nuclear.Model(Z);
                                  initialConfigs=initialConfigs, finalConfigs=finalConfigs,
                                  processSettings=augerSettings)
        return( perform(comp, output=true) )
    end

    if    printout  atomic_code()
    else
          results = redirect_stdout(devnull) do
                        atomic_code()  end
          lines = results["AutoIonization lines:"]
          # Display branching fractions grouped by initial level
          nx    = 64
          hline = "  " * "-"^nx
          println(" ")
          println("  Auger branching fractions (Coulomb interaction):")
          println(" ")
          println(hline)
          println(@sprintf("  %-14s  %-18s  %12s  %10s",
                            "i-level-f", "i--J^P--f", "e_kin [eV]", "BF [%]"))
          println(hline)
          iIndices = unique(l.initialLevel.index for l in lines)
          for  iIdx  in  iIndices
              iLines = [l for l in lines if l.initialLevel.index == iIdx]
              sumRate = sum(l.totalRate for l in iLines)
              il  = iLines[1].initialLevel
              sym = LevelSymmetry(il.J, il.parity)
              println("  Level $(il.index)  ($(string(sym)))")
              for  l  in  iLines
                  bf   = sumRate > 0. ? 100.0 * l.totalRate / sumRate : 0.
                  Ekin = Defaults.convertUnits("energy: from atomic", l.electronEnergy)
                  symI = string(LevelSymmetry(l.initialLevel.J, l.initialLevel.parity))
                  symF = string(LevelSymmetry(l.finalLevel.J,   l.finalLevel.parity))
                  lf   = @sprintf("%3i --> %3i", il.index, l.finalLevel.index)
                  jp   = @sprintf("%-5s --> %-7s", symI, symF)
                  println(@sprintf("  %-14s  %-18s  %12.4e  %10.4f", lf, jp, Ekin, bf))
              end
              println(@sprintf("  %-14s  %-18s  %12s  %10.2f", "", "", "Sum:", 100.))
              println(hline)
          end
    end

    return( nothing )
end


#################################################################################################################################
#################################################################################################################################


"""
`ForPedestrians.displaySpectrum(theme::Basics.ForPhotoEmission, initialConfigs::Array{Configuration,1}, finalConfigs::Array{Configuration,1})`
    ... computes E1 radiative transition rates and displays the photon emission spectrum as an
        ASCII bar chart (lines sorted by energy, bar height = relative intensity).
        If plotfile is given and Plots.jl is loaded, a figure is also saved.

        Simplified call:   setDefaults("nuclear: charge", 3.0)
                           initialConfigs = [Configuration("1s^2 2p"), Configuration("1s^2 3s")]
                           finalConfigs   = [Configuration("1s^2 2s"), Configuration("1s^2 2p")]
                           displaySpectrum(Basics.ForPhotoEmission(), initialConfigs, finalConfigs,
                                           plotfile="spectrum-photon.pdf")
"""
function displaySpectrum(theme::Basics.ForPhotoEmission, initialConfigs::Array{Configuration,1},
                         finalConfigs::Array{Configuration,1};
                         grid::Radial.Grid=Radial.Grid(true), asfSettings::AsfSettings=AsfSettings(),
                         plotfile::String="", printout::Bool=false)
    configs = copy(initialConfigs);   append!(configs, finalConfigs)
    Basics.checkConfigurations(Basics.NumberOfElectrons(), configs)
    Basics.displayConfigurations(stdout, configs, details = "photon emission spectrum computations")

    sa =    "\n* Display the E1 photon emission spectrum (bar chart, relative intensities); " *
            "the following assumptions/simplifications are made: " *
            "\n    + All rates are based on the electric-dipole (E1) approximation only. " *
            "\n    + BF[%] shown for both Coulomb and Babushkin gauge." *
            "\n    + ASCII bar height = A_Bab(line) / max_line(A_Bab) x 100  (Babushkin gauge)." *
            "\n    + Use the optional argument  plotfile = \"filename.pdf\"  to save a figure. " *
            "\n    + For more elaborate computations, make use of perform(comp::Atomic.Computation) \n"
    println(sa)

    photoSettings = PhotoEmission.Settings(PhotoEmission.Settings(), multipoles=[E1],
                                           gauges=[UseCoulomb, UseBabushkin], printBefore=true)
    function atomic_code()
        Defaults.setDefaults("standard grid", grid)
        Z = Defaults.getDefaults("nuclear: charge")
        comp = Atomic.Computation(Atomic.Computation(), name="Photon emission spectrum",
                                  grid=grid, nuclearModel=Nuclear.Model(Z);
                                  initialConfigs=initialConfigs, finalConfigs=finalConfigs,
                                  processSettings=photoSettings)
        return( perform(comp, output=true) )
    end

    if    printout  atomic_code()
    else
          results = redirect_stdout(devnull) do
                        atomic_code()  end
          lines = results["radiative lines:"]
          # Sort by photon energy; normalise each gauge to its own strongest line
          slines   = sort(lines, by = l -> l.omega, rev=true)
          maxBab   = maximum(l.photonRate.Babushkin for l in slines)
          sumCou   = sum(l.photonRate.Coulomb   for l in slines)
          sumBab   = sum(l.photonRate.Babushkin for l in slines)
          eunit    = Defaults.getDefaults("unit: energy")
          # ASCII bar chart: two BF columns; bar shows Babushkin relative intensity
          println("  Photon emission spectrum  (bar = Babushkin relative intensity):\n")
          println("  " * "-"^102)
          println(@sprintf("  %-14s  %-16s  %10s  %10s  %-42s", "Energy [" * eunit * "]",
                            "i--J^P--f", "BF Cou[%]", "BF Bab[%]", "Rel. intensity (Bab)"))
          println("  " * "-"^102)
          for  l  in  slines
              Eph  = Defaults.convertUnits("energy: from atomic", l.omega)
              symI = string(LevelSymmetry(l.initialLevel.J,  l.initialLevel.parity))
              symF = string(LevelSymmetry(l.finalLevel.J,    l.finalLevel.parity))
              bfC  = sumCou > 0. ? 100. * l.photonRate.Coulomb   / sumCou : 0.
              bfB  = sumBab > 0. ? 100. * l.photonRate.Babushkin / sumBab : 0.
              relB = maxBab > 0. ? l.photonRate.Babushkin / maxBab : 0.
              bars = "█" ^ max(0, round(Int, relB * 40))
              jp   = @sprintf("%-5s->%-5s", symI, symF)
              println(@sprintf("  %-14.4e  %-16s  %10.4f  %10.4f  %-42s", Eph, jp, bfC, bfB, bars))
          end
          println("  " * "-"^102)
          # Optional Plots.jl figure: Coulomb and Babushkin side-by-side coloured bars
          if  !isempty(plotfile)
              if  isdefined(Main, :Plots)
                  Plots   = Main.Plots
                  Eph_vec = [Defaults.convertUnits("energy: from atomic", l.omega) for l in slines]
                  BF_Cou  = [sumCou > 0. ? 100. * l.photonRate.Coulomb   / sumCou : 0. for l in slines]
                  BF_Bab  = [sumBab > 0. ? 100. * l.photonRate.Babushkin / sumBab : 0. for l in slines]
                  Erange  = length(Eph_vec) > 1 ? maximum(Eph_vec) - minimum(Eph_vec) : 1.0
                  bar_w   = max(Erange * 0.04, 0.05) / 3   # visible bar width per gauge
                  p = Plots.bar(Eph_vec .- bar_w/2, BF_Cou, bar_width=bar_w,
                                color=:steelblue, label="Coulomb",
                                xlabel="Photon energy [" * eunit * "]",
                                ylabel="BF [%]", title="Photon emission spectrum")
                  Plots.bar!(Eph_vec .+ bar_w/2, BF_Bab, bar_width=bar_w,
                             color=:darkorange, label="Babushkin")
                  Plots.savefig(p, plotfile)
                  println("  Plot saved to: " * plotfile)
              else
                  println("  Plotting skipped — load Plots.jl first:  using Plots")
              end
          end
    end

    return( nothing )
end


"""
`ForPedestrians.displaySpectrum(theme::Basics.ForAutoIonization, initialConfigs::Array{Configuration,1}, finalConfigs::Array{Configuration,1})`
    ... computes Auger rates and displays the electron (Auger) emission spectrum as an
        ASCII bar chart (lines sorted by kinetic energy, bar height = relative intensity).
        If plotfile is given and Plots.jl is loaded, a figure is also saved.

        Simplified call:   setDefaults("nuclear: charge", 10.0)
                           initialConfigs = [Configuration("1s 2s^2 2p^6")]
                           finalConfigs   = [Configuration("[He] 2p^6"), Configuration("[He] 2s 2p^5"),
                                             Configuration("[He] 2s^2 2p^4")]
                           displaySpectrum(Basics.ForAutoIonization(), initialConfigs, finalConfigs,
                                           plotfile="spectrum-auger.pdf")
"""
function displaySpectrum(theme::Basics.ForAutoIonization, initialConfigs::Array{Configuration,1},
                         finalConfigs::Array{Configuration,1};
                         grid::Radial.Grid=Radial.Grid(), asfSettings::AsfSettings=AsfSettings(),
                         plotfile::String="", printout::Bool=false)
    Basics.checkConfigurations(Basics.NumberOfElectrons(), initialConfigs, finalConfigs)
    configs = copy(initialConfigs);   append!(configs, finalConfigs)
    Basics.displayConfigurations(stdout, configs, details = "Auger electron spectrum computations")

    sa =    "\n* Display the Auger electron emission spectrum (bar chart, relative intensities); " *
            "the following assumptions/simplifications are made: " *
            "\n    + All Auger rates are based on the instantaneous Coulomb interaction only." *
            "\n    + Continuum orbitals: B-spline Galerkin method, pure-sine normalization." *
            "\n    + Bar height = Gamma_A(line) / max_line(Gamma_A) x 100  (relative intensity)." *
            "\n    + Use the optional argument  plotfile = \"filename.pdf\"  to save a figure. " *
            "\n    + For more elaborate computations, make use of perform(comp::Atomic.Computation) \n"
    println(sa)

    if      grid.NoPoints == 0    currentGrid = Radial.Grid(Radial.Grid(false), rnt=1.0e-5, h=5.0e-2, hp=2.0e-2, rbox=15.0)
    else                          currentGrid = grid
    end
    augerSettings = AutoIonization.Settings()

    function atomic_code()
        Defaults.setDefaults("standard grid",         currentGrid)
        Defaults.setDefaults("method: continuum, Galerkin")
        Defaults.setDefaults("method: normalization, pure sine")
        Z = Defaults.getDefaults("nuclear: charge")
        comp = Atomic.Computation(Atomic.Computation(), name="Auger electron spectrum",
                                  grid=currentGrid, nuclearModel=Nuclear.Model(Z);
                                  initialConfigs=initialConfigs, finalConfigs=finalConfigs,
                                  processSettings=augerSettings)
        return( perform(comp, output=true) )
    end

    if    printout  atomic_code()
    else
          results = redirect_stdout(devnull) do
                        atomic_code()  end
          lines   = results["AutoIonization lines:"]
          # Sort by electron kinetic energy (descending) and normalise
          slines  = sort(lines, by = l -> l.electronEnergy, rev=true)
          maxRate = maximum(l.totalRate for l in slines)
          sumRate = sum(l.totalRate for l in slines)
          eunit   = Defaults.getDefaults("unit: energy")
          # ASCII bar chart
          println("  Auger electron spectrum  (Coulomb interaction, relative to strongest line):\n")
          println("  " * "-"^84)
          println(@sprintf("  %-14s  %-16s  %6s  %-42s", "e_kin [" * eunit * "]",
                            "i--J^P--f", "BF[%]", "Relative intensity"))
          println("  " * "-"^84)
          for  l  in  slines
              Ekin = Defaults.convertUnits("energy: from atomic", l.electronEnergy)
              symI = string(LevelSymmetry(l.initialLevel.J, l.initialLevel.parity))
              symF = string(LevelSymmetry(l.finalLevel.J,   l.finalLevel.parity))
              rel  = maxRate > 0. ? l.totalRate / maxRate : 0.
              bfC  = sumRate > 0. ? 100. * l.totalRate / sumRate : 0.
              bars = "█" ^ max(0, round(Int, rel * 40))
              jp   = @sprintf("%-5s->%-5s", symI, symF)
              println(@sprintf("  %-14.4e  %-16s  %6.2f  %-42s", Ekin, jp, bfC, bars))
          end
          println("  " * "-"^84)
          # Optional Plots.jl figure
          if  !isempty(plotfile)
              if  isdefined(Main, :Plots)
                  Plots = Main.Plots
                  Ek_vec  = [Defaults.convertUnits("energy: from atomic", l.electronEnergy) for l in slines]
                  BF_vec  = [sumRate > 0. ? 100. * l.totalRate / sumRate : 0. for l in slines]
                  dE      = length(Ek_vec) > 1 ? minimum(diff(sort(Ek_vec))) * 0.3 : 1.0
                  p = Plots.bar(Ek_vec, BF_vec, bar_width=dE, legend=false,
                                xlabel="Electron kinetic energy [" * eunit * "]",
                                ylabel="Relative intensity [%]",
                                title="Auger electron spectrum")
                  Plots.savefig(p, plotfile)
                  println("  Plot saved to: " * plotfile)
              else
                  println("  Plotting skipped — load Plots.jl first:  using Plots")
              end
          end
    end

    return( nothing )
end


#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

"""
`ForPedestrians.displayCouplings(theme::Basics.FineStructure, configs::Array{Configuration,1})` 
    ... displays the (open-shell) configurations along with the total angular momenta J and multiplicities of the associated 
        fine-structure levels. Each configurations is treated separately and can have a different number of electrons. 
        The coupling information is printed to screen but nothing is returned otherwise.
        
        Simplified call:   configs   = [Configuration("[He] 2p^6"), Configuration("[He] 2s 2p^5"), 
                                        Configuration("[He] 2s^2 2p^4"), Configuration("1s 2s 2p^3")]
                           displayCouplings(Basics.FineStructure(), configs)
"""
function displayCouplings(theme::Basics.FineStructure, configs::Array{Configuration,1})
    # Collect explanations
    sa =    "\n* Selected configurations along with the total angular momenta J and multiplicities of the associated " *
            "fine-structure levels: " *
            "\n    + The total J are derived from the subsequent coupling of the (open-subshell) states. " *
            "\n    + Call ? displayConfigurationFineStructure(...) for further details. \n"
    println(sa)

    was = Basics.extractFromConfigurations(NumberOfElectrons(), configs)
    was = sort( unique(was) )
    for  wa in was
        confs = Basics.extractConfigurations(Basics.ByNumber([wa]), configs)
        for  conf in confs
            Basics.displayConfiguration(stdout, theme, conf, header=false)
        end
    end 
    
    return( nothing )
end 

            
"""
`ForPedestrians.displayCouplings(theme::Basics.FineStructureLS, configs::Array{Configuration,1})` 
    ... displays the (open-shell) configurations along with the total angular momenta J and multiplicities of the associated 
        fine-structure levels. Each configurations is treated separately and can have a different number of electrons. 
        The coupling information is printed to screen but nothing is returned otherwise.
        
        Simplified call:   configs   = [Configuration("[He] 2p^6"), Configuration("[He] 2s 2p^5"), 
                                        Configuration("[He] 2s^2 2p^4"), Configuration("1s 2s 2p^3")]
                           displayCouplings(Basics.FineStructureLS(), configs)
"""
function displayCouplings(theme::Basics.FineStructureLS, configs::Array{Configuration,1})
    # Collect explanations
    sa =    "\n* Selected configurations along with the total angular momenta J and multiplicities of the associated " *
            "fine-structure levels: " *
            "\n    + The total J are derived from the subsequent coupling of the (open-subshell) states. " *
            "\n    + Call ? displayCouplings(FineStructureLS(), ...) for further details. \n"
    println(sa)

    was = Basics.extractFromConfigurations(NumberOfElectrons(), configs)
    was = sort( unique(was) )
    for  wa in was
        confs = Basics.extractConfigurations(Basics.ByNumber([wa]), configs)
        for  conf in confs
            Basics.displayConfiguration(stdout, theme, conf, header=false)
        end
    end 
    
    return( nothing )
end 



#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

"""
`ForPedestrians.estimateCrossSections(theme::Basics.ForImpactIonization, initialConfigs::Array{Configuration,1};
                                      grid::Radial.Grid=Radial.Grid(true), printout::Bool=false)` 
    ... estimates the electron impact-ionization cross sections of all shells in the given configurations.
        The results are printed to screen but nothing is returned otherwise.
        
        Simplified call:   setDefaults("nuclear: charge", 10.0)
                           initialConfigs = [Configuration("1s^2 2s^2 2p^6")]
                           estimateCrossSections(Basics.ForImpactIonization(), initialConfigs)
                           estimateCrossSections(Basics.ForImpactIonization(), initialConfigs,
                                                 electronEnergies = [20., 50., 100., 500., 1000.])
"""
function estimateCrossSections(theme::Basics.ForImpactIonization, initialConfigs::Array{Configuration,1};
                               electronEnergies::Array{Float64,1} = [2.0^(i-1) for i=1:18],
                               grid::Radial.Grid=Radial.Grid(true), printout::Bool=false)
    Basics.checkConfigurations(Basics.NumberOfElectrons(), initialConfigs)
    Basics.displayConfigurations(stdout, initialConfigs, details = "electron impact-ionization cross sections")

    sa =    "\n* Estimate the electron impact-ionization cross sections for the shells of the configurations above; " *
            "the following assumptions/simplifications are made: " *
            "\n    + The relativistic binary-encounter Bethe (BEB) model is applied." *
            "\n    + Impact energies [eV]: $(electronEnergies[1]) ... $(electronEnergies[end])  ($(length(electronEnergies)) points)." *
            "\n    + Use the optional argument  electronEnergies = [...]  to set a custom energy grid [eV]." *
            "\n    + Use the optional argument  printout = true  to generate intermediate printout." *
            "\n    + For more elaborate computations, make use of perform(comp::Empirical.Computation) " *
            "\n    + Call ? Empirical.Computation for further details. " *
            "\n    + Call ? setDefaults ... to define user-specified units for the computations. \n"
    println(sa)

    # Assign physics parameters
    Z           = Defaults.getDefaults("nuclear: charge")
    approx      = ImpactIonization.RelativisticBEBmodel()
    multipleN   = 1
    shells      = Basics.extractFromConfigurations(Basics.AllShells(), initialConfigs)
    selection   = ShellSelection(true, shells, Int64[])
    name        = "EII cross section estimates."
    nucModel    = Nuclear.Model(Z)
    eiiSettings = ImpactIonization.Settings(approx, multipleN, electronEnergies, true, true, selection)

    # Specify the atomic computations
    function atomic_code()
        comp    = Empirical.Computation(name, nucModel, grid, initialConfigs, eiiSettings)
        results = perform(comp, output=true)
        return( results )
    end

    # Print or suppress the standard output
    if    printout  atomic_code()
    else
          results = redirect_stdout(devnull) do
                        atomic_code()  end
          cs = results["EII cross sections:"]
          ImpactIonization.displayCrossSections(stdout, cs, eiiSettings)
    end

    return( nothing )
end


#################################################################################################################################
#################################################################################################################################


"""
`ForPedestrians.computeChargeStateDistribution(theme::Basics.ForStepwiseDecay,
                                               initialConfigs::Array{Configuration,1};
                                               grid::Radial.Grid=Radial.Grid(false),
                                               printout::Bool=false)`
    ... computes the charge state distribution that arises from the stepwise radiative and Auger
        decay of all levels from the given initial configurations (typically a single configuration
        with one inner-shell hole).  Both radiative (E1) and Auger channels are included at each
        step; the cascade is followed until at most theme.maximallyReleased electrons have been
        emitted.  Mean-field (AverageSCA) orbitals are used throughout.  The initial population
        is placed entirely on the lowest cascade level (weight 1.0).
        The final ion distribution is printed to screen; nothing is returned otherwise.

        Simplified call:   setDefaults("nuclear: charge", 10.0)
                           initialConfigs = [Configuration("1s 2s^2 2p^6")]
                           computeChargeStateDistribution(Basics.ForStepwiseDecay(3), initialConfigs)
                           computeChargeStateDistribution(Basics.ForStepwiseDecay(5), initialConfigs,
                                                          grid = Radial.Grid(Radial.Grid(false), rbox=25.0))
"""
function computeChargeStateDistribution(theme::Basics.ForStepwiseDecay,
                                        initialConfigs::Array{Configuration,1};
                                        grid::Radial.Grid=Radial.Grid(),
                                        printout::Bool=false)
    Basics.displayConfigurations(stdout, initialConfigs, details = "charge state distribution computations")

    sa =    "\n* Compute the charge state distribution from stepwise (radiative + Auger) decay; " *
            "the following assumptions/simplifications are made: " *
            "\n    + Both radiative (E1) and Auger channels are included at each decay step." *
            "\n    + Cascade is followed for at most $(theme.maximallyReleased) electron-loss steps." *
            "\n    + Mean-field (AverageSCA) orbitals are used for all configurations." *
            "\n    + Continuum orbitals are generated with the Galerkin method, pure-sine normalization." *
            "\n    + The lowest cascade level carries the full initial population (weight 1.0)." *
            "\n    + Use the optional argument  printout = true  to generate intermediate printout." *
            "\n    + For more elaborate computations, make use of perform(comp::Cascade.Computation) " *
            "\n    + Call ? Cascade.Computation for further details. \n"
    println(sa)

    # Determine a grid suitable for high-energy continuum orbitals in the cascade
    if      grid.NoPoints == 0    currentGrid = Radial.Grid(Radial.Grid(false), rnt=2.0e-6, h=5.0e-2, hp=0.5e-2, rbox=20.0)
    else                          currentGrid = grid
    end

    Z      = Defaults.getDefaults("nuclear: charge")
    scheme = Cascade.StepwiseDecayScheme([Auger(), Radiative()], theme.maximallyReleased,
                                          Dict{Int64,Float64}(), 0, Shell[], Shell[], Shell[])

    function cascade_code()
        Defaults.setDefaults("method: continuum, Galerkin")
        Defaults.setDefaults("method: normalization, pure sine")
        comp = Cascade.Computation(Cascade.Computation(), name="Charge state distribution",
                                   nuclearModel=Nuclear.Model(Z), grid=currentGrid,
                                   scheme=scheme, approach=Cascade.AverageSCA(),
                                   initialConfigs=initialConfigs)
        return( perform(comp; output=true) )
    end

    if    printout   wb = cascade_code()
    else
          wb = redirect_stdout(devnull) do
                   cascade_code()  end
    end

    sim = Cascade.Simulation(Cascade.Simulation(), name="Charge state distribution",
                              property=Cascade.IonDistribution(),
                              settings=Cascade.SimulationSettings(),
                              computationData=[Dict{String,Any}("results" => wb)])
    perform(sim; output=true)

    return( nothing )
end


end  ## module
