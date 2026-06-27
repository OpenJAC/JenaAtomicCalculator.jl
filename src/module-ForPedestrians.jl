
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


using  Printf, ..AngularMomentum, ..Atomic, ..AutoIonization, ..Basics, ..Defaults, ..DielectronicRecombination,
               ..Empirical, ..ImpactIonization, ..ManyElectron,  ..Nuclear, ..PhotoEmission, ..PhotoIonization,
               ..PhotoRecombination, ..Radial

export computeCrossSections,  computeForPedestrians,  computeLevelEnergies,  computeLifetimes,  computeResonanceStrength,
       computeTransitionRates,  displayCouplings,  estimateCrossSections


       
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
    ... computes with minimal (very simplified) input different excitation energies, rates and ionization cross sections.
        Please, use the functions computeLevelEnergies(...),  computeCrossSections(...),  computeTransitionRates(...),
        displayCoulings(...), estimateCrossSections(...), ...
"""
function computeForPedestrians()
        
    sa =    "\n* We here provide a simple-man's (pedestrian) approach to the computation of atomic level energies, transition rates " *
            "and cross sections of various kind: \n" *
            "\n    + Useful functions are:  computeLevelEnergies(...),  computeCrossSections(...),  computeTransitionRates(...),  " *
            "\n                             displayCoulings(...), estimateCrossSections(...)." *
            "\n    + Call ? computeLevelEnergies  ... for an example, the argument and further details. " *
            "\n    " *
            "\n    + Call ? setDefaults ... to define user-specified units for the computations. " *
            "\n    + Call setDefaults(``nuclear: charge'', Z::Float64) ... to define the nuclear charge Z. " *
            "\n    " *
            "\n    + For more elaborate atomic computations, make use of perform(comp::Atomic.Computation) " *
            "\n    + Call ? Atomic.Computation for further details. " *
            "\n    + For more elaborate cascade or plasma computations, make use of perform(comp::Atomic.Computation) or " *
            "\n    + perform(comp::Plasma.Computation), respectively. " *
            "\n    + Call ? Cascade.Computation  or  ? Plasma.Computation  for further details. " *
            "\n    " *
            "\n    + For the algebraic simplification of formal expressions from Racah's algebra, call ? RacahAlgebra.RacahExpression  and " *
            "\n      ? RacahAlgebra.evaluate for further details. " *
            "\n    " *
            "\n    + Donate properly with a reference to: S. Fritzsche, Comp. Phys. Commun. 240, 1–14 (2019).   :)"
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
          pathways   = results["dielectronic recombination pathways:"]
          resonances = DielectronicRecombination.computeResonances(pathways, drSettings)
          DielectronicRecombination.displayResults(stdout, pathways, drSettings)
          DielectronicRecombination.displayResults(stdout, resonances, drSettings)
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

end  ## module
