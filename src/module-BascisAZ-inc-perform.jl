
export  perform

# 7Apr25

"""
`Basics.perform(computation::Atomic.Computation)`  
    ... to perform the computation as prescribed by comp. All relevant intermediate and final results are printed to screen (stdout). 
        Nothing is returned.

`Basics.perform(computation::Atomic.Computation; output=true)`  
    ... to perform the same but to return the complete output in a dictionary; the particular output depends on the type and 
        specifications of the computations but can easily accessed by the keys of this dictionary.
"""
function Basics.perform(computation::Atomic.Computation; output::Bool=false)
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    nModel = computation.nuclearModel

    # Distinguish between the computation of level energies and properties and the simulation of atomic processes
    if   length(computation.configs) != 0
        multiplet = SelfConsistent.performSCF(computation.configs, nModel, computation.grid, computation.asfSettings)
        LSjj.expandLevelsIntoLS(multiplet, computation.asfSettings.jjLS)
        #
        if output    results = Base.merge( results, Dict("multiplet:" => multiplet) ) 
                        results = Base.merge( results, Dict("grid:"      => computation.grid) )  end
        
        # Now compute all requested properties
        for settings  in computation.propertySettings
            if      typeof(settings) == Einstein.Settings
                outcome = Einstein.computeLines(multiplet,        computation.grid, settings)    
                if output    results = Base.merge( results, Dict("Einstein lines:" => outcome) )                  end
                #
            ## elseif  typeof(settings) == Hfs.Settings    && settings.calcIJFexpansion  
            ##     outcome = Hfs.computeHyperfineMultiplet(multiplet, nModel, computation.grid, settings)         
            ##     if output    results = Base.merge( results, Dict("IJF multiplet:" => outcome) )                   end
            ##     #
            elseif  typeof(settings) == Hfs.Settings
                outcome = Hfs.computeOutcomes(multiplet, nModel,  computation.grid, settings)         
                if output    results = Base.merge( results, Dict("HFS outcomes:" => outcome) )                    end
                #
            elseif  typeof(settings) == LandeZeeman.Settings 
                outcome = LandeZeeman.computeOutcomes(multiplet, nModel,  computation.grid, settings)      
                if output    results = Base.merge( results, Dict("Zeeman parameter outcomes:" => outcome) )       end
                #
            elseif  typeof(settings) == StarkShift.Settings 
                outcome = StarkShift.computeOutcomes(multiplet, nModel,  computation.grid, settings)      
                if output    results = Base.merge( results, Dict("Stark-shift outcomes:" => outcome) )            end
                #
            elseif  typeof(settings) == IsotopeShift.Settings 
                outcome = IsotopeShift.computeOutcomes(multiplet, nModel, computation.grid, settings)         
                if output    results = Base.merge( results, Dict("Isotope parameter outcomes:" => outcome) )      end
                #
            elseif  typeof(settings) == AlphaVariation.Settings 
                outcome = AlphaVariation.computeOutcomes(multiplet, nModel, computation.grid, settings)         
                if output    results = Base.merge( results, Dict("alpha variation parameter outcomes:" => outcome) )      end
                #
            elseif  typeof(settings) == FormFactor.Settings 
                outcome = FormFactor.computeOutcomes(multiplet, nModel, computation.grid, settings)         
                if output    results = Base.merge( results, Dict("Form factor outcomes:" => outcome) )            end
                #
            elseif  typeof(settings) == DecayYield.Settings 
                outcome = DecayYield.computeOutcomes(computation.configs, computation.asfSettings, 
                                                        multiplet, nModel, computation.grid, settings)     
                if output    results = Base.merge( results, Dict("Fluorescence and AutoIonization yield outcomes:" => outcome) )   end
                #
            elseif  typeof(settings) == MultipolePolarizibility.Settings 
                outcome = MultipolePolarizibility.computeOutcomes(multiplet, nModel, computation.grid, settings)         
                if output    results = Base.merge( results, Dict("Polarizibility outcomes:" => outcome) )         end
                #
            elseif  typeof(settings) == ReducedDensityMatrix.Settings 
                outcome = ReducedDensityMatrix.computeOutcomes(multiplet, nModel, computation.grid, settings)         
                if output    results = Base.merge( results, Dict("RDM outcomes:" => outcome) )         end
                #
            end
        end
        
    else
        initialMultiplet = SelfConsistent.performSCF(computation.initialConfigs, nModel, computation.grid, computation.initialAsfSettings)
        LSjj.expandLevelsIntoLS(initialMultiplet, computation.initialAsfSettings.jjLS)
        finalMultiplet   = SelfConsistent.performSCF(computation.finalConfigs, nModel, computation.grid, computation.finalAsfSettings)
        LSjj.expandLevelsIntoLS(finalMultiplet, computation.finalAsfSettings.jjLS)
        if  output   results["initialMultiplet"] = initialMultiplet;   results["finalMultiplet"] = finalMultiplet    end 
        #
        if typeof(computation.processSettings) in [PhotoExcitationFluores.Settings, PhotoExcitationAutoion.Settings, PhotoIonizationFluores.Settings, 
                                                   PhotoIonizationAutoion.Settings, ImpactExcitationAutoion.Settings, 
                                                   DielectronicRecombination.Settings, ResonantInelastic.Settings]
            intermediateMultiplet = SelfConsistent.performSCF(computation.intermediateConfigs, nModel, computation.grid, computation.intermediateAsfSettings)
            if  output   results["intermediateMultiplet"] = intermediateMultiplet    end 
        end
        #
        if      typeof(computation.processSettings) == AutoIonization.Settings 
            outcome = AutoIonization.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("AutoIonization lines:" => outcome) )                  end
        elseif  typeof(computation.processSettings) == RayleighCompton.Settings 
            outcome = RayleighCompton.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("Rayleigh-Compton lines:" => outcome) )                end
        elseif  typeof(computation.processSettings) == DoubleAutoIonization.Settings   
            outcome = DoubleAutoIonization.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("Double-Auger lines:" => outcome) )                    end
        elseif  typeof(computation.processSettings) == DielectronicRecombination.Settings 
            outcome = DielectronicRecombination.computePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, nModel, 
                                                                computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("dielectronic recombination pathways:" => outcome) )   end
        elseif  typeof(computation.processSettings) == MultiPhotonDeExcitation.Settings
            outcome = MultiPhotonDeExcitation.computeLines(computation.processSettings.process,
                                                            finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("multi-photon-de-excitation lines:" => outcome) )      end
        elseif  typeof(computation.processSettings) == PhotoIonization.Settings   
            outcome = PhotoIonization.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("photoionization lines:" => outcome) )                 end
        elseif  typeof(computation.processSettings) == PhotoDoubleIonization.Settings   
            outcome = PhotoDoubleIonization.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("Single-photon double-ionization lines:" => outcome) )      end
        elseif  typeof(computation.processSettings) == PhotoExcitation.Settings
            outcome = PhotoExcitation.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("photo-excitation lines:" => outcome) )                     end
        elseif  typeof(computation.processSettings) == PhotoExcitationAutoion.Settings  
            outcome = PhotoExcitationAutoion.computePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, nModel, 
                                                                computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("photo-excitation-autoionization pathways:" => outcome) )   end
        elseif  typeof(computation.processSettings) == PhotoExcitationFluores.Settings
            outcome = PhotoExcitationFluores.computePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, 
                                                                computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("photo-excitation-fluorescence pathways:" => outcome) )     end
        elseif  typeof(computation.processSettings) == PhotoEmission.Settings
            outcome = PhotoEmission.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("radiative lines:" => outcome) )                            end
        elseif  typeof(computation.processSettings) == ResonantInelastic.Settings 
            outcome = ResonantInelastic.computePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, 
                                                        computation.grid, computation.processSettings) 
        elseif  typeof(computation.processSettings) == CoulombExcitation.Settings
            outcome = CoulombExcitation.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("Coulomb excitation lines:" => outcome) )                   end
        elseif  typeof(computation.processSettings) == RadiativeAuger.Settings
            outcome = RadiativeAuger.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("radiative Auger sharings:" => outcome) )                   end
        elseif  typeof(computation.processSettings) == PhotoRecombination.Settings 
            outcome = PhotoRecombination.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("photo recombination lines:" => outcome) )                  end
        elseif  typeof(computation.processSettings) == ImpactExcitation.Settings 
            outcome = ImpactExcitation.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("impact-excitation lines:" => outcome) )                    end
        elseif  typeof(computation.processSettings) == InternalRecombination.Settings 
            outcome = InternalRecombination.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("internal-recombination lines:" => outcome) )          end
        elseif  typeof(computation.processSettings) == TwoElectronOnePhoton.Settings 
            outcome = TwoElectronOnePhoton.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("two-electron-one-photon lines:" => outcome) )         end
        elseif  typeof(computation.processSettings) == ParticleScattering.Settings 
            outcome = ParticleScattering.computeEvents(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("particle-scattering events:" => outcome) )         end
        elseif  typeof(computation.processSettings) == BeamPhotoExcitation.Settings 
            outcome = BeamPhotoExcitation.computeOutcomes(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("beam-assisted photo-excitation:" => outcome) )         end
        elseif  typeof(computation.processSettings) == HyperfineInduced.Settings 
            outcome = HyperfineInduced.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("hyperfine-induced transitions:" => outcome) )         end
            #
            #
        elseif  typeof(computation.processSettings) == PairA1P()  
            outcome = PairAnnihilation1Photon.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("pair-annihilation-1-photon lines:" => outcome) )       end
        elseif  typeof(computation.processSettings) == Coulex()
            outcome = CoulombExcitation.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("Coulomb excitation lines:" => outcome) )               end
        elseif  typeof(computation.processSettings) == Coulion()
            outcome = CoulombIonization.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("Coulomb ionization lines:" => outcome) )               end
        elseif  typeof(computation.processSettings) == PhotoIonAuto()   
            outcome = PhotoIonizationAutoion.computePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, 
                                                                computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("photo-ionization-autoionization pathways:" => outcome) )      end
        elseif  typeof(computation.processSettings) == PhotoIonFluor()  
            outcome = PhotoIonizationFluores.computePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, 
                                                                computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("photo-ionizatiton-fluorescence pathways:" => outcome) )        end
        elseif  typeof(computation.processSettings) == ImpactExcAuto()  
            outcome = ImpactExcitationAutoion.computePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, 
                                                                computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("impact-excitation-autoionization pathways:" => outcome) )     end
        elseif  typeof(computation.processSettings) == MultiPI()   
            outcome = MultiPhotonIonization.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("multi-photon single ionization:" => outcome) )        end
        elseif  typeof(computation.processSettings) == MultiPDI()   
            outcome = MultiPhotonDoubleIon.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("multi-photon double ionization:" => outcome) )        end
        elseif  typeof(computation.processSettings) == InternalConv()   
            outcome = InternalConversion.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("internal conversion lines:" => outcome) )        end
        else
            error("stop b")
        end
    end
    
    Defaults.warn(PrintWarnings())
    Defaults.warn(ResetWarnings())
    return( results )
end 



#== 6Apr2025
"""
`Basics.perform(comp::Empirical.Computation)`  
    ... to set-up and perform an empirical computation that starts from a given nuclear model and set of configurations,
        and which is mainly controlled by its settings. The results are printed to screen but nothing is returned otherwise.

`Basics.perform(comp::Empirical.Computation; output=true)`   
    ... to perform the same but to return the complete output in a dictionary;  the particular output depends on the kind
        and specifications of the empirical computation but can easily accessed by the keys of this dictionary.
"""
function Basics.perform(comp::Empirical.Computation; output::Bool=false)
    Empirical.perform(comp::Empirical.Computation, output=output)
end ==#



#== 6Apr2025
"""
`Basics.perform("computation: mutiplet from orbitals, no CI, CSF diagonal", configs::Array{Configuration,1}, 
                initalOrbitals::Dict{Subshell, Orbital}, nuclearModel::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings; printout::Bool=true)` 
    ... to generate from the given initial orbitals a multiplet of single-CSF levels by just using the diagonal 
        part of the Hamiltonian matrix; a multiplet::Multiplet is returned.  
"""
function Basics.perform(sa::String, configs::Array{Configuration,1}, initalOrbitals::Dict{Subshell, Orbital}, nuclearModel::Nuclear.Model, 
                grid::Radial.Grid, settings::AsfSettings; printout::Bool=true)
    !(sa == "computation: mutiplet from orbitals, no CI, CSF diagonal")   &&   error("Unsupported keystring = $sa")
    return( Basics.performCI(configs::Array{Configuration,1}, initalOrbitals::Dict{Subshell, Orbital}, nuclearModel::Nuclear.Model, 
                                grid::Radial.Grid, settings::AsfSettings; printout=printout) )
end ==#                


#== 7Apr25
"""
`Basics.performSCF(configs::Array{Configuration,1}, nuclearModel::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings;
                    printout::Bool=true)`  
    ... to generate an atomic basis and to compute the self-consistent field (SCF) for this basis due to the given settings; 
        a basis::Basis is returned.   
"""
function Basics.performSCF(configs::Array{Configuration,1}, nuclearModel::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings;
                            printout::Bool=true)
    if  printout    println("\n... in performSCF ...")    end
    
    # Generate a list of relativistic configurations and determine an ordered list of subshells for these configurations
    relconfList = ConfigurationR[]
    for  conf in configs
        wa = Basics.generateConfigurationRs(conf)
        append!( relconfList, wa)
    end
    if  printout    for  i = 1:length(relconfList)    println(">> include ", relconfList[i])    end   end
    subshellList = Basics.generateSubshellList(relconfList)
    Defaults.setDefaults("relativistic subshell list", subshellList; printout=printout)

    # Generate the relativistic CSF's for the given subshell list
    csfList = CsfR[]
    for  relconf in relconfList
        newCsfs = Basics.generateCsfRs(relconf, subshellList)
        append!( csfList, newCsfs)
    end

    # Determine the number of electrons and the list of coreSubshells
    NoElectrons      = sum( csfList[1].occupation )
    coreSubshellList = Subshell[]
    for  k in 1:length(subshellList)
        mocc = Basics.subshell_2j(subshellList[k]) + 1;    is_filled = true
        for  csf in csfList
            if  csf.occupation[k] != mocc    is_filled = false;    break   end
        end
        if   is_filled    push!( coreSubshellList, subshellList[k])    end
    end
        
    wa = Bsplines.generatePrimitives(grid)

    # Generate start orbitals
    if  typeof(settings.startScfFrom) == StartFromHydrogenic
        if  printout   println("Start SCF process with hydrogenic orbitals.")   end
        # Generate start orbitals for the SCF field by using B-splines
        orbitals  = Bsplines.generateOrbitalsHydrogenic(wa, nuclearModel, subshellList; printout=printout)
    elseif  typeof(settings.startScfFrom) == StartFromPrevious
        if  printout   println("Start SCF process from given list of orbitals.energy")    end
        # Taking starting orbitals for the given dictionary; non-relativistic orbitals with a proper nuclear charge
        # are adapted if no orbital is found
        orbitals = Dict{Subshell, Orbital}()
        for  subsh in subshellList
            if  haskey(settings.startScfFrom.orbitals, subsh)  
                orb = settings.startScfFrom.orbitals[subsh]
                orbitals = Base.merge( orbitals, Dict( subsh => orb) )
            else
                println("Start orbitals do not contain an Orbital for subshell $subsh ")
                orb      = HydrogenicIon.radialOrbital(subsh, nuclearModel.Z, grid)
                orb      = Orbital(orb.subshell, orb.isBound, true, orb.energy, orb.P, orb.Q, orb.Pprime, orb.Qprime, Radial.Grid())
                orbitals = Base.merge( orbitals, Dict( subsh => orb) ) 
            end
        end
    else  error("stop a")
    end

    # Perform the SCF computations for the orbitals due to the given settings
    basis    = Basis(true, NoElectrons, subshellList, csfList, coreSubshellList, orbitals)
    if   typeof(settings.scField)  in  [Basics.DFSField, Basics.DFSwCPField, Basics.HSField]
        basis = Bsplines.solveSelfConsistentMeanField(wa, nuclearModel, basis, settings; printout=printout) 
    elseif   settings.scField in [Basics.NuclearField()]  && settings.startScfFrom == StartFromHydrogenic() 
        # Return the basis as already generated.
    elseif   settings.scField in [Basics.ALField()] 
        basis = Bsplines.solveSelfConsistentALField(wa, nuclearModel, basis, settings; printout=printout)
    elseif   typeof(settings.scField) == Basics.EOLField 
        basis = Bsplines.solveSelfConsistentEOLField(wa, nuclearModel, basis, settings; printout=printout)
    else  error("stop b")
    end
    
    return( basis )  
end  ==#
                


#== 7Apr25
"""
`Basics.performSCF(basis::Basis, nuclearModel::Nuclear.Model, grid::Radial.Grid, frozenShells::Array{Shell,1}, 
                    settings::AtomicState.RasSettings; printout::Bool=true)`  
    ... to generate an atomic basis and to compute the self-consistent field (SCF) for this basis due to the given settings;
        all subshells due to the list of frozenShells are kept fixed. A newBasis::Basis is returned.  
"""
function Basics.performSCF(basis::Basis, nuclearModel::Nuclear.Model, grid::Radial.Grid, frozenShells::Array{Shell,1}, 
                            settings::AtomicState.RasSettings; printout::Bool=true)
    if  printout    println("\n... in performSCF['for RAS step'] ...")    end
    
    # Determine the list of frozen subshells
    frozenSubshells = Subshell[]
    for  sh in frozenShells
        if  sh.l == 0    push!( frozenSubshells, Subshell(sh.n, -1))
        else             push!( frozenSubshells, Subshell(sh.n, sh.l));    push!( frozenSubshells, Subshell(sh.n, -sh.l -1))
        end
    end
    
    wa = Bsplines.generatePrimitives(grid)
    # Specify the AsfSettings for the given RAS step
    asfSettings = AsfSettings(true, CoulombInteraction(), Basics.DFSField(), StartFromHydrogenic(),    
                                settings.maxIterationsScf, settings.accuracyScf, Subshell[], frozenSubshells, 
                                CoulombInteraction(), NoneQed(), LSjjSettings(false), settings.levelSelectionCI ) 

    basis = Bsplines.solveSelfConsistentMeanField(wa, nuclearModel, basis, asfSettings; printout=printout) 
    
    return( basis )  
end  ==#


#== 7Apr25
"""
`Basics.performCI(configs::Array{Configuration,1}, initalOrbitals::Dict{Subshell, Orbital}, nuclearModel::Nuclear.Model, 
                    grid::Radial.Grid, settings::AsfSettings; printout::Bool=true)` 
    ... to generate from the given initial orbitals a multiplet of single-CSF levels by just using the diagonal 
        part of the Hamiltonian matrix; a multiplet::Multiplet is returned.  
"""
function Basics.performCI(configs::Array{Configuration,1}, initalOrbitals::Dict{Subshell, Orbital}, nuclearModel::Nuclear.Model, 
                            grid::Radial.Grid, settings::AsfSettings; printout::Bool=true)
    if  printout    println("\n... in perform['computation: mutiplet from orbitals, no CI, CSF diagonal'] ...")    end
    
    if  settings.eeInteractionCI == DiagonalCoulomb()   error("Not yet implemented.")   end
    
    # Generate a list of relativistic configurations and determine an ordered list of subshells for these configurations
    relconfList = ConfigurationR[]
    for  conf in configs
        wa = Basics.generateConfigurationRs(conf)
        append!( relconfList, wa)
    end
    if  printout    for  i = 1:length(relconfList)    println(">> include ", relconfList[i])    end   end
    subshellList = Basics.generateSubshellList(relconfList)
    Defaults.setDefaults("relativistic subshell list", subshellList; printout=printout)

    # Generate the relativistic CSF's for the given subshell list
    csfList = CsfR[]
    for  relconf in relconfList
        newCsfs = Basics.generateCsfRs(relconf, subshellList)
        append!( csfList, newCsfs)
    end

    # Determine the number of electrons and the list of coreSubshells
    NoElectrons      = sum( csfList[1].occupation )
    coreSubshellList = Subshell[]
    for  k in 1:length(subshellList)
        mocc = Basics.subshell_2j(subshellList[k]) + 1;    is_filled = true
        for  csf in csfList
            if  csf.occupation[k] != mocc    is_filled = false;    break   end
        end
        if   is_filled    push!( coreSubshellList, subshellList[k])    end
    end
    
    # Set-up a basis for calculating the Hamiltonian matrix
    basis = Basis(true, NoElectrons, subshellList, csfList, coreSubshellList, initalOrbitals)

    # Calculate the diagonal part of the Hamiltonian matrix and define a multiplet for these approximate ASf;
    # Generate first an effective nuclear charge Z(r) on the given grid
    potential = Nuclear.nuclearPotential(nuclearModel, grid)
    
    keep = true
    InteractionStrength.XL_Coulomb_reset_storage(keep)
    #
    n = length(csfList);    matrix = zeros(Float64, n, n)
    for  r = 1:n
        # Calculate the spin-angular coefficients
        if  Defaults.saRatip()
            waR = compute("angular coefficients: e-e, Ratip2013", basis.csfs[r], basis.csfs[r])
            wa  = waR       
        end
        if  Defaults.saGG()
            subshellList = basis.subshells
            opa  = SpinAngular.OneParticleOperator(0, plus, true)
            waG1 = SpinAngular.computeCoefficients(opa, basis.csfs[r], basis.csfs[r], subshellList) 
            opa  = SpinAngular.TwoParticleOperator(0, plus, true)
            waG2 = SpinAngular.computeCoefficients(opa, basis.csfs[r], basis.csfs[r], subshellList)
            wa   = [waG1, waG2]
        end
        if  Defaults.saRatip() && Defaults.saGG() && true
            if  length(waR[1]) != 0     println(  ">> Angular coeffients from Ratip2013   = $(waR[1]) ")    end
            if  length(waG1)   != 0     println("\n>> Angular coeffients from SpinAngular = $waG1 ")        end
            if  length(waR[2]) != 0     println(  ">> Angular coeffients from Ratip2013   = $(waR[2]) ")    end
            if  length(waG2)   != 0     println("\n>> Angular coeffients from SpinAngular = $waG2 ")        end
        end
        #
        me = 0.
        for  coeff in wa[1]
            jj = Basics.subshell_2j(basis.orbitals[coeff.a].subshell)
            me = me + coeff.T * sqrt( jj + 1) * RadialIntegrals.GrantIab(basis.orbitals[coeff.a], basis.orbitals[coeff.b], grid, potential)
        end

        for  coeff in wa[2]
            if  typeof(settings.eeInteractionCI) in [CoulombInteraction, CoulombBreit, CoulombGaunt]
                me = me + coeff.V * InteractionStrength.XL_Coulomb(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                                basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid, keep=keep)  end
            if  typeof(settings.eeInteractionCI) in [BreitInteraction, CoulombBreit]
                me = me + coeff.V * InteractionStrength.XL_Breit(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                           basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid,
                                                                           settings.eeInteractionCI)                                                end
        end
        matrix[r,r] = me
    end 
    #
    eigen  = Basics.diagonalize("matrix: LinearAlgebra", matrix)
    levels = Level[]
    for  ev = 1:length(eigen.values)
        vector = eigen.vectors[ev];   ns = 0
        for  r = 1:length(vector)   
            if      (vector[r])^2 > 0.99    ns = r;   break    
            elseif  (vector[r])^2 > 0.01    error("stop a; unexpected mixing coefficieint; vector = $vector ")    
            end
        end
        J = csfList[ns].J
        newlevel = Level( J, AngularM64(J.num//J.den), csfList[ns].parity, 0, eigen.values[ev], 0., true, basis, vector ) 
        push!( levels, newlevel)
    end
    mp = Multiplet("noName", levels)
    mp = Basics.sortByEnergy(mp)
    
    levelNos = Int64[]
    for (ilev, level) in  enumerate(mp.levels)       push!(levelNos, ilev)   end
    
    # Display all level energies and energy splittings
    if  printout
        Basics.tabulate(stdout, "multiplet: energies", mp, levelNos)
        Basics.tabulate(stdout, "multiplet: energy relative to immediately lower level",    mp, levelNos)
        Basics.tabulate(stdout, "multiplet: energy of each level relative to lowest level", mp, levelNos)
    end
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary  &&  printout    
        Basics.tabulate(iostream, "multiplet: energies", mp, levelNos)
        Basics.tabulate(iostream, "multiplet: energy relative to immediately lower level",    mp, levelNos)
        Basics.tabulate(iostream, "multiplet: energy of each level relative to lowest level", mp, levelNos)
    end

    return( mp )
end  ==#


#== 7Apr25
"""
`Basics.performCI(basis::Basis, nuclearModel::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings; printout::Bool=true)`  
    ... to  set-up and diagonalize from the (SCF) basis the configuration-interaction matrix and to derive and display the 
        level structure of the corresponding multiplet due to the given settings; a multiplet::Multiplet is returned.   
"""
function Basics.performCI(basis::Basis, nuclearModel::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings; printout::Bool=true)
    
    # Determine the J^P symmetry blocks
    symmetries = Dict{LevelSymmetry,Int64}()
    for  csf in basis.csfs
        sym = LevelSymmetry(csf.J, csf.parity)
        if     haskey(symmetries, sym)    symmetries[sym] = symmetries[sym] + 1
        else   symmetries = Base.merge( symmetries, Dict( sym => 1 ) )
        end
    end

    # Test the total number of CSF
    NoCsf = 0
    for (sym,v) in symmetries   NoCsf = NoCsf + v   end

    # Calculate for each symmetry block the corresponding CI matrix, diagonalize it and append a Multiplet for this block
    multiplets = Multiplet[]
    for  (sym,v) in  symmetries
        # Skip the symmetry block if it not selected
            if  !Basics.selectSymmetry(sym, settings.levelSelectionCI)     continue    end
        matrix = compute("matrix: CI, J^P symmetry", sym, basis, nuclearModel, grid, settings; printout=printout)
        #
        eigen  = Basics.diagonalize("matrix: LinearAlgebra", matrix)
        levels = Level[]
        for  ev = 1:length(eigen.values)
            # Construct the eigenvector with regard to the given basis (not w.r.t the symmetry block)
            evSym = eigen.vectors[ev];    vector = zeros( length(basis.csfs) );   ns = 0
            for  r = 1:length(basis.csfs) 
                if LevelSymmetry(basis.csfs[r].J, basis.csfs[r].parity) == sym    ns = ns + 1;   vector[r] = evSym[ns]   end
            end
            newlevel = Level( sym.J, AngularM64(sym.J.num//sym.J.den), sym.parity, 0, eigen.values[ev], 0., true, basis, vector ) 
            push!( levels, newlevel)
        end
        wa = Multiplet("noName", levels)
        push!( multiplets, wa)
    end
    
    # Merge all multiplets into a single one
    mp = Basics.merge(multiplets)
    mp = Basics.sortByEnergy(mp)
    
    # Determine the level list to be printed out
    levelNos = Int64[]
    for (ilev, level) in  enumerate(mp.levels)
        if  Basics.selectLevel(level, settings.levelSelectionCI)    push!(levelNos, ilev)    end
    end
    
    # Display all level energies and energy splittings
    if  printout
        Basics.tabulate(stdout, "multiplet: energies", mp, levelNos)
        Basics.tabulate(stdout, "multiplet: energy relative to immediately lower level",    mp, levelNos)
        Basics.tabulate(stdout, "multiplet: energy of each level relative to lowest level", mp, levelNos)
    end
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary     
        Basics.tabulate(iostream, "multiplet: energies", mp, levelNos)
        Basics.tabulate(iostream, "multiplet: energy relative to immediately lower level",    mp, levelNos)
        Basics.tabulate(iostream, "multiplet: energy of each level relative to lowest level", mp, levelNos)
    end

    return( mp )
end  ==#


#== 7Apr25
"""
`Basics.performCI(basis::Basis, nuclearModel::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings, 
                  plasmaModel::Basics.AbstractPlasmaModel; printout::Bool=true)`  
    ... to  set-up and diagonalize from the given (SCF) basis the configuration-interaction matrix and to derive and
        display the level structure of the corresponding multiplet due to the given settings. Here, the CI matrix
        includes the modifications of the Hamiltonian due to the given plasmaSettings; a multiplet::Multiplet is returned.   
"""
function Basics.performCI(basis::Basis, nuclearModel::Nuclear.Model, grid::Radial.Grid, 
                          settings::AsfSettings, plasmaModel::Basics.AbstractPlasmaModel; printout::Bool=true)
    
    # Determine the J^P symmetry blocks
    symmetries = Dict{LevelSymmetry,Int64}()
    for  csf in basis.csfs
        sym = LevelSymmetry(csf.J, csf.parity)
        if     haskey(symmetries, sym)    symmetries[sym] = symmetries[sym] + 1
        else   symmetries = Base.merge( symmetries, Dict( sym => 1 ) )
        end
    end

    # Test the total number of CSF
    NoCsf = 0
    for (sym,v) in symmetries   NoCsf = NoCsf + v   end

    # Calculate for each symmetry block the corresponding CI matrix for the given plasma model and parameters as specified by
    # plasmaSettings, diagonalize it and append a Multiplet for this block
    multiplets = Multiplet[]
    for  (sym,v) in  symmetries
        ## matrix = compute("matrix: CI for plasma, J^P symmetry", sym, basis, nuclearModel, grid, settings, plasmaSettings; printout=printout)
        matrix = compute(sym, basis, nuclearModel, grid, settings, plasmaModel; printout=printout)
        eigen  = Basics.diagonalize("matrix: LinearAlgebra", matrix)
        levels = Level[]
        for  ev = 1:length(eigen.values)
            # Construct the eigenvector with regard to the given basis (not w.r.t the symmetry block)
            evSym = eigen.vectors[ev];    vector = zeros( length(basis.csfs) );   ns = 0
            for  r = 1:length(basis.csfs) 
                if LevelSymmetry(basis.csfs[r].J, basis.csfs[r].parity) == sym    ns = ns + 1;   vector[r] = evSym[ns]   end
            end
            newlevel = Level( sym.J, AngularM64(sym.J.num//sym.J.den), sym.parity, 0, eigen.values[ev], 0., true, basis, vector ) 
            push!( levels, newlevel)
        end
        wa = Multiplet("noName", levels)
        push!( multiplets, wa)
    end
    
    # Merge all multiplets into a single one
    mp = Basics.merge(multiplets)
    mp = Basics.sortByEnergy(mp)
    
    levelNos = Int64[]
    for (ilev, level) in  enumerate(mp.levels)       push!(levelNos, ilev)   end
    
    # Display all level energies and energy splittings
    if  printout
        Basics.tabulate(stdout, "multiplet: energies", mp, levelNos)
        Basics.tabulate(stdout, "multiplet: energy relative to immediately lower level",    mp, levelNos)
        Basics.tabulate(stdout, "multiplet: energy of each level relative to lowest level", mp, levelNos)
    end
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary     
        Basics.tabulate(iostream, "multiplet: energies", mp, levelNos)
        Basics.tabulate(iostream, "multiplet: energy relative to immediately lower level",    mp, levelNos)
        Basics.tabulate(iostream, "multiplet: energy of each level relative to lowest level", mp, levelNos)
    end

    return( mp )
end  ==#
