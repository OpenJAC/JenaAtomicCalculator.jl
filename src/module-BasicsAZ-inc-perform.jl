
export  perform

# 7Apr25

"""
`Basics.perform(computation::Atomic.Computation)`  
    ... to perform the computation as prescribed by comp. All relevant intermediate and final results are printed to screen (stdout). 
        Nothing is returned.

`Basics.perform(computation::Atomic.Computation; output::Bool=false)`  
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
            elseif  typeof(settings) == StarkZeeman.Settings
                outcome = StarkZeeman.computeOutcomes(multiplet, nModel, computation.grid, settings)
                if output    results = Base.merge( results, Dict("Stark-Zeeman outcomes:" => outcome) )           end
                #
            elseif  typeof(settings) == IsotopeShift.Settings 
                outcome = IsotopeShift.computeOutcomes(multiplet, nModel, computation.grid, settings)         
                if output    results = Base.merge( results, Dict("Isotope parameter outcomes:" => outcome) )      end
                #
            elseif  typeof(settings) == AlphaVariation.Settings
                outcome = AlphaVariation.computeOutcomes(multiplet, nModel, computation.grid, computation.configs,
                                                          computation.asfSettings, settings)
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
            elseif  typeof(settings) == WeakInteractionEnhancement.Settings
                outcome = WeakInteractionEnhancement.computeOutcomes(multiplet, nModel, computation.grid, settings)
                if output    results = Base.merge( results, Dict("Weak-interaction enhancement outcomes:" => outcome) )   end
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
        if typeof(computation.processSettings) in [PhotoExcitationFluores.Settings, PhotoExcitationAutoion.Settings,
                                                   ImpactExcitationAutoion.Settings,
                                                   DielectronicRecombination.Settings, ResonantInelastic.Settings,
                                                   PhotoRecombinationInterference.Settings]
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
            ## settings.calcHyperfineResolved selects the hyperfine-resolved route, which builds the hyperfine
            ## multiplets itself from these same three ELECTRONIC multiplets and recouples the electronic
            ## amplitudes; it is a strict analogue of the fine-structure one and needs no separate entry.
            if  computation.processSettings.calcHyperfineResolved
                outcome = DielectronicRecombination.computeHfCaptureLines(finalMultiplet, intermediateMultiplet, initialMultiplet,
                                                                          nModel, computation.grid, computation.processSettings)
                if output    results = Base.merge( results, Dict("hyperfine-resolved dielectronic recombination lines:" => outcome) )   end
            else
                outcome = DielectronicRecombination.computeCaptureLines(finalMultiplet, intermediateMultiplet, initialMultiplet, nModel,
                                                                        computation.grid, computation.processSettings)
                if output    results = Base.merge( results, Dict("dielectronic recombination lines:" => outcome) )   end
            end
        elseif  typeof(computation.processSettings) == MultiPhotonTransition.Settings
            outcome = MultiPhotonTransition.computeLines(computation.processSettings.scheme,
                                                            finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("multi-photon transition lines:" => outcome) )         end
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
            ## The Green-function (intermediate) multiplet is NOT generated here; it is prepared by the user in a
            ## separate, standard Atomic.Computation -- whose configurations are chosen to include just those levels
            ## that are expected to contribute strongly to the second-order amplitude -- and handed over explicitly
            ## as settings.gMultiplet.
            outcome = TwoElectronOnePhoton.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings)
            if output    results = Base.merge( results, Dict("two-electron-one-photon lines:" => outcome) )         end
        elseif  typeof(computation.processSettings) == ParticleScattering.Settings 
            outcome = ParticleScattering.computeEvents(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("particle-scattering events:" => outcome) )         end
        elseif  typeof(computation.processSettings) == PhotonScattering.Settings
            outcome = PhotonScattering.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings)
            if output    results = Base.merge( results, Dict("photon-scattering lines:" => outcome) )               end
        elseif  typeof(computation.processSettings) == BeamPhotoExcitation.Settings 
            outcome = BeamPhotoExcitation.computeOutcomes(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("beam-assisted photo-excitation:" => outcome) )         end
        elseif  typeof(computation.processSettings) == HyperfineInduced.Settings 
            outcome = HyperfineInduced.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("hyperfine-induced transitions:" => outcome) )         end
            #
            #
        ## THE SEVEN TESTS BELOW USED TO COMPARE A TYPE WITH AN INSTANCE -- `typeof(x) == Coulex()` and its six
        ## siblings -- which is false for every x, so all seven processes fell through to error("stop b") and were
        ## unreachable.  They also tested the wrong THING: computation.processSettings is declared
        ## ::Basics.AbstractProcessSettings, so it is a Settings object and never one of the AbstractProcess
        ## singletons.  Each now tests the Settings type of the module whose function the branch calls, which is
        ## the form the rest of this chain already uses.  Fixed 24-Aug-2026.
        elseif  typeof(computation.processSettings) == CoulombExcitation.Settings
            outcome = CoulombExcitation.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("Coulomb excitation lines:" => outcome) )               end
        elseif  typeof(computation.processSettings) == CoulombIonization.Settings
            ## Reserved and empty until 25-Aug-2026, when the module was written; this branch used to raise and say
            ## that CoulombIonization held no computation to reach.  It now does.
            outcome = CoulombIonization.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid,
                                                     computation.processSettings)
            if output    results = Base.merge( results, Dict("Coulomb ionization lines:" => outcome) )           end
        elseif  typeof(computation.processSettings) == ImpactExcitationAutoion.Settings
            outcome = ImpactExcitationAutoion.computePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, 
                                                                computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("impact-excitation-autoionization pathways:" => outcome) )     end
        elseif  typeof(computation.processSettings) == MultiPhotonIonization.Settings
            outcome = MultiPhotonIonization.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings) 
            if output    results = Base.merge( results, Dict("multi-photon single ionization:" => outcome) )        end
        elseif  typeof(computation.processSettings) == InternalConversion.Settings
            outcome = InternalConversion.computeLines(finalMultiplet, initialMultiplet, nModel, computation.grid, computation.processSettings)
            if output    results = Base.merge( results, Dict("internal conversion lines:" => outcome) )        end
        elseif  typeof(computation.processSettings) == CrystalFieldEmission.Settings
            outcome = CrystalFieldEmission.computeLines(finalMultiplet, initialMultiplet, computation.grid, computation.processSettings)
            if output    results = Base.merge( results, Dict("crystal-field-resolved emission lines:" => outcome) )    end
        elseif  typeof(computation.processSettings) == PhotoRecombinationInterference.Settings
            outcome = PhotoRecombinationInterference.computePathways(finalMultiplet, intermediateMultiplet, initialMultiplet, nModel,
                                                                     computation.grid, computation.processSettings)
            if output    results = Base.merge( results, Dict("photorecombination-interference pathways:" => outcome) )    end
        else
            error("stop b")
        end
    end
    
    Defaults.warn(PrintWarnings())
    Defaults.warn(ResetWarnings())
    return( results )
end
