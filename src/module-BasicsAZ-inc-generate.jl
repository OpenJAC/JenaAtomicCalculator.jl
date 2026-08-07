
using   QuadGK, ..Hamiltonian
export  generate


"""
`Basics.generate(representation::AtomicState.Representation)`  
    ... to generate an atomic representation as specified by the representation.repType::AbstractRepresentationType.
        All relevant intermediate and final results are printed to screen (stdout). Nothing is returned.

`Basics.perform(representation::AtomicState.Representation; output=true)`  
    ... to generate the same but to return the complete output in a dictionary; the particular output depends on the type and 
        specifications of the representation but can easily accessed by the keys of this dictionary.
"""
function Basics.generate(representation::AtomicState.Representation; output::Bool=false)
    results = Basics.generate(representation.repType, representation; output=output)
    
    Defaults.warn(PrintWarnings())
    Defaults.warn(ResetWarnings())
    println(" ")
    return( results )
end



"""
`Basics.generate(repType::AtomicState.MeanFieldBasis, representation::AtomicState.Representation)`  
    ... to generate a mean-field basis (representation) for a set of reference configurations; all relevant intermediate 
        and final results are printed to screen (stdout). Nothing is returned.

`Basics.generate(repType::AtomicState.MeanFieldBasis, representation::AtomicState.Representation; output=true)`  
    ... to generate the same but to return the complete output in a dictionary; the particular output depends on the type and 
        specifications of the representation but can easily accessed by the keys of this dictionary.
"""
function Basics.generate(repType::AtomicState.MeanFieldBasis, rep::AtomicState.Representation; output::Bool=false)
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    nModel    = rep.nuclearModel

    # The asfSettings only define the SCF part and are partly derived from the MeanFieldSettings
    asfSettings = AsfSettings(AsfSettings(); scField=repType.settings.scField) 
    
    multiplet      = SelfConsistent.performSCF(rep.refConfigs, nModel, rep.grid, asfSettings; printout=true)
    basis          = multiplet.levels[1].basis
    if output    results = Base.merge( results, Dict("mean-field basis" => basis) )          end
    
    return( results )
end


"""
`Basics.generate(repType::AtomicState.MeanFieldMultiplet, representation::AtomicState.Representation)`  
    ... to generate a mean-field basis (representation) for a set of reference configurations; all relevant intermediate 
        and final results are printed to screen (stdout). Nothing is returned.

`Basics.generate(repType::AtomicState.MeanFieldMultiplet, representation::AtomicState.Representation; output=true)`  
    ... to generate the same but to return the complete output in a dictionary; the particular output depends on the type and 
        specifications of the representation but can easily accessed by the keys of this dictionary.
"""
function Basics.generate(repType::AtomicState.MeanFieldMultiplet, rep::AtomicState.Representation; output::Bool=false)
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    nModel    = rep.nuclearModel

    # The asfSettings only define the SCF part and are partly derived from the MeanFieldSettings
    asfSettings = AsfSettings(AsfSettings(); scField=repType.settings.scField) 
    
    multiplet  = SelfConsistent.performSCF(rep.refConfigs, nModel, rep.grid, asfSettings; printout=true)
    if output    
        ## `levels`, not `level` (fixed 06-Aug-2026): this typo made generate(MeanFieldMultiplet, ...) raise a
        ## FieldError for every call with output=true, i.e. the only useful call. It survived because the one
        ## example that exercised this path -- the old branch 1 of example-Dh.jl -- had never been run.
        results = Base.merge( results, Dict("mean-field basis" => multiplet.levels[1].basis) )
        results = Base.merge( results, Dict("mean-field multiplet" => multiplet) )          
    end
    
    return( results )
end




"""
`Basics.generate(repType::AtomicState.OneElectronSpectrum, representation::AtomicState.Representation)`  
    ... to generate a one-electron spectrum for the atomic potential from the (given) levels, based on a set of reference 
        configurations as well as for given settings. Relevant intermediate and final results are printed to screen (stdout). 
        Nothing is returned in this case.

    + `(repType::AtomicState.OneElectronSpectrum, representation::AtomicState.Representation; output=true)`  
    ... to generate the same but to return the complete output in a orbitals::Dict{Subshell, Orbital}.
"""
function Basics.generate(repType::AtomicState.OneElectronSpectrum, rep::AtomicState.Representation; output::Bool=true)
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    nModel    = rep.nuclearModel
    settings  = repType.settings
    # First perform a SCF+CI computations for the reference configurations below to generate a spectrum of start orbitals
    asfSettings   = AsfSettings()   ## Use default settings to define a first multiplet from the reference configurations;
                                    ## all further details are specified for each step
    refMultiplet  = SelfConsistent.performSCF(rep.refConfigs, nModel, rep.grid, asfSettings; printout=true)
    refBasis      = refMultiplet.levels[1].basis
    Basics.display(stdout, refBasis.orbitals, rep.grid; longTable=false)
    nuclearPot    = Nuclear.nuclearPotential(nModel, rep.grid)
    @warn("The potential is generated for the mean-field basis but not (yetc hosen for the selected levels.")
    ## A NuclearField carries no electronic part at all; adding a DFS screening term on top of it would silently
    ## reintroduce the very mismatch this field exists to remove.
    if  typeof(settings.scField) == Basics.NuclearField
        electronicPot = Radial.Potential("zero", zeros(length(rep.grid.r)), rep.grid)
    else
        electronicPot = Basics.computePotential(Basics.DFSField(1.0), rep.grid, refBasis)
    end
    meanPot       = Basics.add(nuclearPot, electronicPot)
    
    println("")
    printstyled("Compute an one-electron spectrum for the selected partial waves ... \n", color=:light_green)
    printstyled("------------------------------------------------------------------- \n", color=:light_green)
    
    # Generate all non-relativistic and relativistic subshells and the associated single-electron spectrum for this potential
    shellList = Basics.generateShellList(1, settings.nMax, settings.lValues)
    subshellList = Subshell[]
    for  shell in shellList     append!(subshellList, Basics.shellSplitIntoSubshells(shell))    end
    primitives = Bsplines.generatePrimitives(rep.grid)
    orbitals   = Bsplines.generateOrbitals(subshellList, meanPot, nModel, primitives; printout=true)
    
    # Print all results to screen
    Basics.display(stdout, orbitals, rep.grid; longTable=true)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    Basics.display(iostream, orbitals, rep.grid; longTable=false)         end
    
    if output    
        results = Base.merge( results, Dict("mean potential" => meanPot) )              
        results = Base.merge( results, Dict("orbitals" => orbitals) )              
    end
    
    return( results )
end



"""
`Basics.generate(repType::AtomicState.CiExpansion, representation::AtomicState.Representation)`  
    ... to generate a configuration-interaction expansion for a single level symmetry and based on a set of reference configurations
        and a number of pre-specified steps. All relevant intermediate and final results are printed to screen (stdout). 
        Nothing is returned.

`Basics.generate(repType::AtomicState.CiExpansion, representation::AtomicState.Representation; output=true)`  
    ... to generate the same but to return the complete output in a dictionary; the particular output depends on the type and 
        specifications of the computations but can easily accessed by the keys of this dictionary.
"""
function Basics.generate(repType::AtomicState.CiExpansion, rep::AtomicState.Representation; output::Bool=false)
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    nModel    = rep.nuclearModel
    orbitals  = repType.applyOrbitals
    
    # Generate a list of relativistic configurations and  CSF's for the given subshell list
    relconfList = ConfigurationR[]
    for  conf in rep.refConfigs
        wa = Basics.generateConfigurations(Basics.RelativisticConfigurations(), conf)
        append!( relconfList, wa)
    end
    subshellList = Basics.generateSubshellList(relconfList)
    csfList = CsfR[]
    for  relconf in relconfList
        newCsfs = Basics.generateCsfRs(relconf, subshellList)
        append!( csfList, newCsfs)
    end
    
    # Generate all level symmetries that need to be taken into account
    if     length(repType.settings.levelSelectionCI.symmetries) != 0    symmetries = repType.settings.levelSelectionCI.symmetries
    else   symmetries = LevelSymmetry[]
        for csf in csfList
            if      LevelSymmetry(csf.J, csf.parity)  in  symmetries
            else    push!( symmetries, LevelSymmetry(csf.J, csf.parity) )
            end
        end
    end
    println("*** Level symmetries = $symmetries ")

    # The asfSettings only define the CI part and are partly derived from the CiSettings
    asfSettings = AsfSettings(true, CoulombInteraction(), Basics.DFSField(), StartFromHydrogenic(),    0, 0., Subshell[], Subshell[], 
                                repType.settings.eeInteractionCI, NoneQed(), LSjjSettings(false), 
                                repType.settings.levelSelectionCI) 
    
    basis      = Basics.generateBasis(rep.refConfigs, symmetries, repType.excitations)
    basis      = Basis( true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, orbitals )  
    # Test that all orbitals are available 
    for  subshell in  basis.subshells   
        if   !(haskey(orbitals, subshell))  error("No orbital available for subshell $subshell ")     end
    end
    multiplet  = Hamiltonian.performCI(basis,  nModel, rep.grid, asfSettings; printout=true) 
    if output    results = Base.merge( results, Dict("CI multiplet" => Multiplet("CI multiplet:", multiplet.levels)) )              end
    
    return( results )
end



"""
`Basics.generate(repType::AtomicState.RasExpansion, representation::AtomicState.Representation)`  
    ... to generate a restricted active-space expansion for a single level symmetry and based on a set of reference configurations
        and a number of pre-specified steps. All relevant intermediate and final results are printed to screen (stdout). 
        Nothing is returned.

`Basics.generate(repType::AtomicState.RasExpansion, representation::AtomicState.Representation; output=true)`  
    ... to generate the same but to return the complete output in a dictionary; the particular output depends on the type and 
        specifications of the computations but can easily accessed by the keys of this dictionary.
"""
function Basics.generate(repType::AtomicState.RasExpansion, rep::AtomicState.Representation; output::Bool=false)
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    nModel   = rep.nuclearModel
    # First perform a SCF+CI computations for the reference configurations below to generate a spectrum of start orbitals
    asfSettings    = AsfSettings()  ## Use default settings to define a first multiplet from the reference configurations;
                                    ## all further details are specified for each step
    priorMultiplet = SelfConsistent.performSCF(rep.refConfigs, nModel, rep.grid, asfSettings; printout=true)
    nuclearPot     = Nuclear.nuclearPotential(nModel, rep.grid)
    ## electronicPot  = Basics.compute("radial potential: Dirac-Fock-Slater", rep.grid, priorMultiplet.levels[1].basis)
    ## meanPot        = Basics.add(nuclearPot, electronicPot)
    subshellList   = Basics.extractRelativisticSubshellList(rep)             ## extract all subshells that occur in the RAS computation
    primitives     = Bsplines.generatePrimitives(rep.grid)
    startOrbitals  = Bsplines.generateOrbitals(subshellList, nuclearPot, nModel, primitives, printout=true)  ## generate a spectrum of sufficient size
    if output    results = Base.merge( results, Dict("reference multiplet" => Multiplet("Reference multiplet:", priorMultiplet.levels) ) )  end

    # Now, cycle over all steps of the RasExpansion. Each step re-optimizes ONLY the shells newly introduced
    # at this layer (EOL, so the SAME target level(s) specified once in repType.settings.levelSelectionCI stay
    # the optimization target as the CSF space grows -- an AL/configuration-average functional would let
    # higher-lying states entering the growing space dilute the correlation for the levels actually wanted),
    # freezing every shell already present/optimized in an earlier step -- GRASP's own default behavior
    # (LFIX(:NW)=.TRUE. unless explicitly varied).
    for (istep, step)  in  enumerate(repType.steps)
        println("")
        printstyled("++ Compute the orbitals, orbitals and multiplet for step $istep ... \n", color=:light_green)
        printstyled("--------------------------------------------------------------      \n", color=:light_green)
        basis      = Basics.generateBasis(rep.refConfigs, repType.symmetries, step)
        orbitals   = Basics.generateOrbitalsForBasis(basis, step.frozenShells, priorMultiplet.levels[1].basis, startOrbitals)
        basis      = Basis( true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, orbitals )

        # step.frozenShells is a list of non-relativistic Shell(n,l); translate to the concrete, relativistic
        # Subshell(n,kappa) instances actually present in THIS step's basis, as AsfSettings.frozenSubshells needs.
        frozenSubshellsThisStep = [ sh  for  shell in step.frozenShells  for sh in basis.subshells
                                        if  sh.n == shell.n  &&  Basics.subshell_l(sh) == shell.l ]

        stepSettings = AsfSettings( AsfSettings();  scField=Basics.EOLField(),  frozenSubshells=frozenSubshellsThisStep,
                                     eeInteractionCI=repType.settings.eeInteractionCI,  levelSelectionCI=repType.settings.levelSelectionCI,
                                     maxIterationsScf=repType.settings.maxIterationsScf,  accuracyScf=repType.settings.accuracyScf )

        multiplet  = SelfConsistent.performSCF(basis, nModel, rep.grid, stepSettings; printout=true)
        if output    results = Base.merge( results, Dict("step"*string(istep) => Multiplet("Multiplet:", multiplet.levels)) )              end
        priorMultiplet = multiplet
    end
    
    return( results )
end



"""
`Basics.generate(repType::AtomicState.GreenExpansion, representation::AtomicState.Representation)`  
    ... to generate a Green (function) expansion for a given approach and excitation scheme of the electron,
        based on a set of reference configurations, a list of level symmetries as well as for given settings.
        All relevant intermediate and final results are printed to screen (stdout). Nothing is returned.

`Basics.generate(repType::AtomicState.GreenExpansion, representation::AtomicState.Representation; output=true)`  
    ... to generate the same but to return the complete output in a dictionary; the particular output depends on the type and 
        specifications of the representation but can easily accessed by the keys of this dictionary.
"""
function Basics.generate(repType::AtomicState.GreenExpansion, rep::AtomicState.Representation; output::Bool=false)
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    nModel    = rep.nuclearModel
    settings  = repType.settings
    # First perform a SCF+CI computations for the reference configurations below to generate a spectrum of start orbitals
    ## THE POTENTIAL IS NOW TAKEN FROM settings.scField (07-Aug-2026), not hardcoded. Both places below used to
    ## be fixed to DFS -- the reference AsfSettings and the mean potential -- so a Green expansion could never be
    ## matched to the potential of the computation that consumes it. That matters because a Green expansion is
    ## normally used as the INTERMEDIATE spectrum of a second-order calculation, where gauge invariance requires
    ## the same one-body Hamiltonian throughout.
    asfSettings   = AsfSettings(AsfSettings(); scField=settings.scField)
    refMultiplet  = SelfConsistent.performSCF(rep.refConfigs, nModel, rep.grid, asfSettings; printout=true)
    refBasis      = refMultiplet.levels[1].basis
    nuclearPot    = Nuclear.nuclearPotential(nModel, rep.grid)
    ## A NuclearField carries no electronic part at all; adding a DFS screening term on top of it would silently
    ## reintroduce the very mismatch this field exists to remove.
    if  typeof(settings.scField) == Basics.NuclearField
        electronicPot = Radial.Potential("zero", zeros(length(rep.grid.r)), rep.grid)
    else
        electronicPot = Basics.computePotential(Basics.DFSField(1.0), rep.grid, refBasis)
    end
    meanPot       = Basics.add(nuclearPot, electronicPot)
    
    println("")
    printstyled("Compute an approximate Green function expansion ... \n", color=:light_green)
    printstyled("--------------------------------------------------- \n", color=:light_green)
    
    channels = AtomicState.GreenChannel[]
    
    # Generate all (non-relativistic) configurations from the bound configurations due to the given excitation scheme 
    confList = Basics.generateConfigurationsForExcitationScheme(rep.refConfigs, repType.excitationScheme, settings.nMax, settings.lValues)
    # Print (if required) information about the generated configuration list
    if  settings.printBefore    Basics.displayConfigurations(stdout, confList)    end
    
    # Generate shell list abd a full single-electron spectrum for this potential
    subshellList = Basics.extractRelativisticSubshellList(confList)                      ## extract all subshells that occur in confList
    primitives   = Bsplines.generatePrimitives(rep.grid)
    orbitals     = Bsplines.generateOrbitals(subshellList, meanPot, nModel, primitives, printout=true) ## generate a spectrum of sufficient size
    Defaults.setDefaults("relativistic subshell list", subshellList; printout=true)
    Basics.display(stdout, orbitals, rep.grid)

    # The asfSettings only define the CI part of the Green channels and are partly derived from the GreenSettings
    asfSettings = AsfSettings(true, CoulombInteraction(), Basics.DFSField(), StartFromHydrogenic(),    0, 0., Subshell[], Subshell[], 
                                CoulombInteraction(), NoneQed(), LSjjSettings(false), settings.levelSelection ) 
    
    # Cycle over all selected level symmetries to generate the requested channels
    for  levelSymmetry  in  repType.levelSymmetries
        basis      = Basics.generateBasis(confList, [levelSymmetry])
        basis      = Basis( true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, orbitals )  
        multiplet  = Basics.computeMultipletForGreenApproach(repType.approach, basis, nModel, rep.grid, asfSettings, settings; printout=true) 
        push!( channels, GreenChannel( levelSymmetry, multiplet) )
    end        
    
    # Print all results to screen
    Basics.display(stdout, channels)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    Basics.display(iostream, channels)         end
    
    if output    results = Base.merge( results, Dict("Green channels" => channels) )   end
    return( results )
end



"""
`Basics.generate(::CondensedMultiplet, multiplet::Multiplet)`
    ... to condense/reduce the number of CSF in the basis of the given multiplet due to a single 'weight';
        a multiplet::Multiplet is returned.  **Not yet implemented !**
"""
function Basics.generate(::CondensedMultiplet, multiplet::Multiplet)
    error("Not yet implemented !")
end



"""
`Basics.generate(::ConfigurationListNRFromBasis, basis::Basis)`
    ... to (re-) generate the list of NR configurations from the given basis; a confList::Array{Configuration,1} is returned.
"""
function Basics.generate(::ConfigurationListNRFromBasis, basis::Basis)
    confList    = Configuration[]
    NoElectrons = sum( basis.csfs[1].occupation )
    for  csf in basis.csfs
        shellList = Dict{Shell,Int64}[]
        for  i = 1:length(csf.occupation)
            n = csf.subshell[i].n;    l = Basics.subshell_l( csf.subshell[i] );    sh = Shell(n, l)
            if     haskey(shellList, sh)    shellList(sh) = shellList(sh) + csf.occupation[i]
            else   shellList = merge( shellList, Dict(sh => csf.occupation[i]) )
            end
        end
        confNew = Configuration( shellList, NoElectrons )
        push!( confList, confNew )
    end
    confListNew = ConfigurationsExcludeDoubles(confList)
    
    return( confListNew )
end


"""
`Basics.generate(::ConfigurationListNRFromConfiguration, refConf::Configuration, NoExcitations::Int64, fromShells::Array{Shell,1},`
                    toShells::Array{Shell,1})
    ... to generate a non-relativistic configuration list, including the given reference configuration (refConf) and with
        all configurations that differ by NoExcitations from the fromShells into the toShells; an Array{Configuration,1}
        is returned.
"""
function Basics.generate(::ConfigurationListNRFromConfiguration, refConf::Configuration, NoExcitations::Int64, fromShells::Array{Shell,1}, toShells::Array{Shell,1})
    confList = [refConf]
    # First prepare a proper reference configuration that also includes all fromShells and toShells with zero occupation
    shellDict = deepcopy(refConf.shells)
    allShells = deepcopy(fromShells)
    for  sh in toShells
        if     sh in allShells
        else   append!(allShells, [sh] )
        end
    end
    
    for  sh in allShells
    if  !haskey(shellDict, sh )    shellDict = merge( shellDict, Dict( sh => 0 ) )    end
    end

    rconf = Configuration(shellDict, refConf.NoElectrons)
    println("aa", rconf)

    # Add single excitations
    if      NoExcitations == 0
    elseif  NoExcitations == 1
        for  fromsh in fromShells
            for  tosh in toShells
                confnew = deepcopy(rconf)
                confnew.shells[fromsh] = confnew.shells[fromsh] - 1
                confnew.shells[tosh]   = confnew.shells[tosh]   + 1
                println("bb", confnew)

                # Check that this is valid configuration and append, if appropriate
                add = true
                for (k,v) in  confnew.shells
                    if  v <= -1  ||  Parity(rconf) != Parity(confnew)    add = false    end
                end
                if  add 
                    for (k,v) in confnew.shells   
                    if  v == 0   delete!(confnew.shells, k)   end
                    end
                    push!(confList, confnew);   println("fromsh = $fromsh  tosh = $tosh")   
                end
            end
        end

    # Add double excitations
    elseif  NoExcitations == 2
        for  fromsha in fromShells
        for  fromshb in fromShells
            for  tosha in toShells
            for  toshb in toShells
                confnew = deepcopy(rconf)
                confnew.shells[fromsha] = confnew.shells[fromsha] - 1
                confnew.shells[fromshb] = confnew.shells[fromshb] - 1
                confnew.shells[tosha]   = confnew.shells[tosha]   + 1
                confnew.shells[toshb]   = confnew.shells[toshb]   + 1
                println("bb", confnew)

                # Check that this is valid configuration and append, if appropriate
                add = true
                for (k,v) in  confnew.shells
                    if  v <= -1  ||  Parity(rconf) != Parity(confnew)    add = false    end
                end
                if  add    
                    for (k,v) in confnew.shells   
                    if  v == 0   delete!(confnew.shells, k)   end
                    end
                    push!(confList, confnew)    
                end
            end
            end
        end
        end

    else
        error("stop a")
    end

    # Exclude all configurations that appear twice (and more) in the list
    confList = Basics.ConfigurationsExcludeDoubles(confList)

    return( confList )
end


"""
`Basics.generateCsfRs(conf::ConfigurationR, subshellList::Array{Subshell,1})` 
    ... to construct from a given (relativistic) configuration all possible CSF with regard to the subshell order as specified 
        by subshellList; a list::Array{CsfR,1} is returned.
"""
function Basics.generateCsfRs(conf::ConfigurationR, subshellList::Array{Subshell,1})
    parity  = Basics.extractFromConfiguration(Basics.GetParity(), conf)
    csfList = CsfR[];   useStandardSubshells = true;    first = true;    previousCsfs = CsfR[]
    # subhshellList = Subshell[];   
    for  subsh in subshellList
        if   subsh in keys(conf.subshells)    occ = conf.subshells[subsh]    else    occ = 0    end
        if   first
            stateList   = ManyElectron.provideSubshellStates(subsh, occ)
            currentCsfs = CsfR[]
            for  state in stateList
                push!( currentCsfs, CsfR( true, AngularJ64(state.Jsub2//2), parity, [state.occ], [state.nu],
                                            [AngularJ64(state.Jsub2//2)], [AngularJ64(state.Jsub2//2)], Subshell[]) )
            end
            previousCsfs = copy(currentCsfs)
            first        = false
        else
            # Now support also all couplings of the subshell states with the CSFs that were built-up so far
            stateList   = ManyElectron.provideSubshellStates(subsh, occ)
            currentCsfs = CsfR[]
            for  csf in  previousCsfs
                for  state in stateList
                    occupation = deepcopy(csf.occupation);    seniorityNr = deepcopy(csf.seniorityNr);    
                    subshellJ  = deepcopy(csf.subshellJ);     subshells = deepcopy(csf.subshells)
                    push!(occupation, state.occ);   push!(seniorityNr, state.nu);   push!(subshellJ, AngularJ64(state.Jsub2//2) ) 
                    push!(subshells, subsh)
                    newXList = oplus( csf.subshellX[end], AngularJ64(state.Jsub2//2) )
                    for  newX in newXList
                        subshellX = deepcopy(csf.subshellX);   push!(subshellX, newX) 
                        push!( currentCsfs, CsfR( true, subshellX[end], parity, occupation, seniorityNr, subshellJ, subshellX, Subshell[]) ) 
                    end
                end
            end
            previousCsfs = copy(currentCsfs)
        end
    end
    
    return( previousCsfs )
end
    


"""
`Basics.generate(::OrderedShellList, confs::Array{Configuration,1})`
    ... to generate for confs, i.e. all the given (non-relativistic) configurations, a common and ordered shell list;
        a list::Array{Shell,1} is returned.
"""
function Basics.generate(::OrderedShellList, confs::Array{Configuration,1})
    shells = Shell[]

    wa = Defaults.getDefaults("ordered shell list: non-relativistic", 11)
    ## wa = Defaults.getDefaults("ordered shell list: non-relativistic", 19)
    for  a in wa
        for  cf in 1:length(confs)
            ks = keys(confs[cf].shells)
            if  a in ks   push!(shells, a);    break    end
        end 
    end

    return( shells )
end


"""
`Basics.generateSubshellList(confs::Array{ConfigurationR,1})`  
    ... to generate for confs, i.e. all the given (relativistic) configurations, common and ordered subshell list; 
        a list::Array{Subshell,1} is returned.
"""
function Basics.generateSubshellList(confs::Array{ConfigurationR,1})
    subshells = Subshell[]   

    for  conf in confs
        for  subsh in keys(conf.subshells)      push!(subshells, subsh)     end
    end
    subshells = Base.unique(subshells)
    subshells = Base.sort( subshells, lt=Base.isless)

    return( subshells )
end



"""
`Basics.generate(::OrderedSubshellList, basisA::Basis, basisB::Basis)`
    ... to generate common and ordered subshell list for the two basis A and B; a list::Array{Subshell,1} is returned.
"""
function Basics.generate(::OrderedSubshellList, basisA::Basis, basisB::Basis)
    function areEqual(nn::Int64, sha::Array{Subshell,1}, shb::Array{Subshell,1})
        # Determines whether the first nn subshells are equal in sha and shb (true) or not (false)
        for  i = 1:nx
            if    sha[i] != shb[i]    return( false )   end
        end
        return( true )
    end

    subshells = Subshell[]

    nx = min(length(basisA.subshells), length(basisB.subshells))
    if  areEqual(nx, basisA.subshells, basisB.subshells)
        # If subshell order is equal, any subshell order is accepted for those subshells that occur in both basis
        for  i = 1:nx   
            if    basisA.subshells[i] != basisB.subshells[i]   error("Inconsistent subshells of two bases.")
            else  push!( subshells, basisA.subshells[i])
            end
        end
        #
        if       length(basisA.subshells) > nx   
            for  i = nx+1:length(basisA.subshells)    push!(subshells, basisA.subshells[i])    end
        elseif   length(basisB.subshells) > nx   
            for  i = nx+1:length(basisB.subshells)    push!(subshells, basisB.subshells[i])    end
        end
    else
        ## println("basisA.subshells = $(basisA.subshells)")
        ## println("basisB.subshells = $(basisB.subshells)")
        # If subshell order is NOT equal, all subshell in basisA and basisB must follow standard order
        standardList = Defaults.getDefaults("ordered subshell list: relativistic", 7)
        na = 0;  nb = 0
        for  sh in standardList   
            if  sh in basisA.subshells  ||   sh in basisB.subshells   push!( subshells, sh)    end    
        end
        nn = 0;   
        for sh in basisA.subshells   
            if    !(sh in subshells) push!( subshells, sh)    
            else  wb = findall(x->x==sh, subshells);    if wb[1] <= nn   error("stop a")   else   nn = wb[1]  end    
            end
        end
        nn = 0;   
        for sh in basisB.subshells   
            if    !(sh in subshells) push!( subshells, sh)    
            else  wb = findall(x->x==sh, subshells);    if wb[1] <= nn   error("stop a")   else   nn = wb[1]  end    
            end
        end
        println(">>> Extended subshells from two basis = $subshells")
    end

    return( subshells )
end


"""
`Basics.generate(::SlaterTypeSpectrum, N::Int64, potential::Radial.Potential, grid::Radial.Grid; N_0::Int64=30, alpha_0::Float64=1.0,`
                    beta_0::Float64=1.1)
    ... to generate a complete one-electron spectrum with N positive and N negative states, and by using even-tempered Slater-type
        orbitals (STO) with parameters ``\\alpha_i = \\alpha_0 \\beta_0^i``; a spectrum::SingleElecSpectrum is returned where just
        N0 positive and N_0 negative are kept for later use.  **Not yet implemented !**
"""
function Basics.generate(::SlaterTypeSpectrum, N::Int64, potential::Radial.Potential, grid::Radial.Grid; N_0::Int64=30, alpha_0::Float64=1.0, beta_0::Float64=1.1)
    error("Not yet implemented !")
    return( nothing )
end


"""
`Basics.generate(::SlaterTypeSpectrumPositive, N::Int64, potential::Radial.Potential, grid::Radial.Grid; N_0::Int64=30, alpha_0::Float64=1.0,`
                    beta_0::Float64=1.1)
    ... to generate the same but to return only the N_0 positive states.  **Not yet implemented !**
"""
function Basics.generate(::SlaterTypeSpectrumPositive, N::Int64, potential::Radial.Potential, grid::Radial.Grid; N_0::Int64=30, alpha_0::Float64=1.0, beta_0::Float64=1.1)
    error("Not yet implemented !")
    return( nothing )
end



"""
`Basics.generateBasis(confList::Array{Configuration,1}, symmetries::Array{LevelSymmetry,1})`  
    ... generates the CSF basis for the given configuration list and level symmetries; a basis::Basis is returned but 
        without a valid representation of the radial orbitals.
"""
function Basics.generateBasis(confList::Array{Configuration,1}, symmetries::Array{LevelSymmetry,1})
    #
    relconfList = ConfigurationR[]
    for  conf in confList
        wa = Basics.generateConfigurations(Basics.RelativisticConfigurations(), conf)
        append!( relconfList, wa)
    end
    subshellList = Basics.generateSubshellList(relconfList)
    Defaults.setDefaults("relativistic subshell list", subshellList; printout=true)

    # Generate the relativistic CSF's for the given subshell list
    csfList = CsfR[]
    for  relconf in relconfList
        newCsfs = Basics.generateCsfRs(relconf, subshellList)
        for  csf in newCsfs     if  LevelSymmetry(csf.J, csf.parity) in symmetries   push!(csfList, csf)   end   end
    end
    #
    if length(csfList) == 0     error("There are no CSF with $symmetries in the given configuration list.")      end

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
    
    basis = Basis(true, NoElectrons, subshellList, csfList, coreSubshellList, Dict{Subshell, Orbital}())
    
    println("Construct a basis with $(length(basis.csfs)) CSF for J^P = $symmetries with $(length(basis.subshells)) subshells: " *
            "$(basis.subshells[1])  $(basis.subshells[2]) ...  $(basis.subshells[end-1])  $(basis.subshells[end])")
    return( basis )
end


"""
`Basics.generateBasis(refConfigs::Array{Configuration,1}, symmetries::Array{LevelSymmetry,1}, step::AtomicState.RasStep)`  
    ... generates the CSF basis for the given reference configurations and single, double, ... excitations of
        electrons from the corresponding fromShells --> toShells. A basis::Basis is returned but without a valid
        representation of the radial orbitals
"""
function Basics.generateBasis(refConfigs::Array{Configuration,1}, symmetries::Array{LevelSymmetry,1}, step::AtomicState.RasStep)
    # Single excitations
    if      step.seFrom == Shell[]  ||  step.seTo == Shell[]    confSingles = Configuration[]
    else    confSingles = Basics.generateConfigurations(refConfigs, step.seFrom, step.seTo)
    end
    # Double excitations
    if      step.deFrom == Shell[]  ||  step.deTo == Shell[]    confDoubles = Configuration[]
    else    confDoubles = Basics.generateConfigurations(refConfigs, step.deFrom, step.deTo)
            confDoubles = Basics.generateConfigurations(confDoubles,     step.deFrom, step.deTo)
    end
    # Triple excitations
    if      step.teFrom == Shell[]  ||  step.teTo == Shell[]    confTriples = Configuration[]
    else    confTriples = Basics.generateConfigurations(refConfigs,      step.teFrom, step.teTo)
            confTriples = Basics.generateConfigurations(confTriples,     step.teFrom, step.teTo)
            confTriples = Basics.generateConfigurations(confTriples,     step.teFrom, step.teTo)
    end
    # Quadruple excitations
    if      step.teFrom == Shell[]  ||  step.teTo == Shell[]    confQuadruples = Configuration[]
    else    confQuadruples = Basics.generateConfigurations(refConfigs, step.qeFrom, step.qeTo)
            confQuadruples = Basics.generateConfigurations(confQuadruples,  step.qeFrom, step.qeTo)
            confQuadruples = Basics.generateConfigurations(confQuadruples,  step.qeFrom, step.qeTo)
            confQuadruples = Basics.generateConfigurations(confQuadruples,  step.qeFrom, step.qeTo)
    end
    
    # Now get a unique set of configurations and generate the relativistic configurations
    # configurations = Base.unique((==), configurations) does not work here
    confList       = Configuration[]  
    configurations = deepcopy(refConfigs);      append!(configurations, confSingles);    append!(configurations, confDoubles)
                                                append!(configurations, confTriples);    append!(configurations, confQuadruples)
    for  confa in configurations
        addTo = true
        for confb in confList    if   confa == confb    addTo = false;    break     end     end
        if  addTo    push!(confList, confa)     end
    end
    for  i = 1:length(confList)    println(">> include ", confList[i])    end
    #
    relconfList = ConfigurationR[]
    for  conf in confList
        wa = Basics.generateConfigurations(Basics.RelativisticConfigurations(), conf)
        append!( relconfList, wa)
    end
    subshellList = Basics.generateSubshellList(relconfList)
    Defaults.setDefaults("relativistic subshell list", subshellList; printout=true)

    # Generate the relativistic CSF's for the given subshell list
    csfList = CsfR[]
    for  relconf in relconfList
        newCsfs = Basics.generateCsfRs(relconf, subshellList)
        for  csf in newCsfs     if  LevelSymmetry(csf.J, csf.parity) in symmetries   push!(csfList, csf)   end   end
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
    
    basis = Basis(true, NoElectrons, subshellList, csfList, coreSubshellList, Dict{Subshell, Orbital}())
    
    println("Construct a basis with $(length(basis.csfs)) CSF for J^P = $symmetries with $(length(basis.subshells)) subshells: " *
            "$(basis.subshells[1])  $(basis.subshells[2]) ...  $(basis.subshells[end-1])  $(basis.subshells[end])")
    return( basis )
end


"""
`Basics.generateConfigurations(refConfigs::Array{Configuration,1}, fromShells::Array{Shell,1}, toShells::Array{Shell,1})`  
    ... generates all nonrelativistic configurations due to the (single) excitation of an electron fromShells --> toShells
        for all given reference configurations. A confList::Array{Configuration,1} is returned.
"""
function Basics.generateConfigurations(refConfigs::Array{Configuration,1}, fromShells::Array{Shell,1}, toShells::Array{Shell,1})
    confList = Configuration[]
    for  config in refConfigs
        for  fromShell in fromShells
            for  toShell in toShells
                # Cycle through all shells of config to built-up a new configuration with a single excitation if fromShell occurs
                if  haskey(config.shells, fromShell )
                    newShells = Dict{Shell,Int64}();   addShell = true;   addConfiguration = true
                    for (k,v) in config.shells
                        if        k == fromShell == toShell       newShells = Base.merge( newShells, Dict( k => v))
                                    addShell = false   
                        elseif    k == fromShell    if  v-1 > 0   newShells = Base.merge( newShells, Dict( fromShell => v-1))  end
                        elseif    k == toShell      
                            if    v+1 <= 2*(2*toShell.l + 1)      newShells = Base.merge( newShells, Dict( toShell => v+1))
                                    addShell = false   
                            else  addConfiguration = false
                            end
                        else      newShells = Base.merge( newShells, Dict( k => v))
                        end
                        #
                        if  addShell                              newShells = Base.merge( newShells, Dict( toShell => 1))      end
                    end
                    if  addConfiguration   push!(confList, Configuration(newShells, config.NoElectrons ))          end
                else    # No configuration to add if fromShell does not occur in the reference configuration
                end
            end
        end
    end
    
    confList = Base.unique(confList)
    
    return( confList )
end



"""
`Basics.generateConfigurations(refConfigs::Array{Configuration,1}, fromShells::Array{Shell,1}, toShells::Array{Shell,1}, 
                                noex::Int64; restrictions::Array{AbstractConfigurationRestriction,1}=AbstractConfigurationRestriction[])`  
    ... generates all nonrelativistic configurations with excitation of up to noex electrons fromShells --> toShells
        for all given reference configurations. Moreover, the list of restrictions (if any is given) is finally applied 
        to restrict the configurations to a given set of limitations. It remains the reponsibility of the user to make sure that
        the given restrictions do not contradict each other and are consistent with what is to be achieved. The given set of 
        restrictions can be easily extended if this need arises by the users. A confList::Array{Configuration,1} is returned.
"""
function Basics.generateConfigurations(refConfigs::Array{Configuration,1}, fromShells::Array{Shell,1}, toShells::Array{Shell,1}, 
                                        noex::Int64; restrictions::Array{AbstractConfigurationRestriction,1}=AbstractConfigurationRestriction[])
    if     noex == 0
        newConfigs = refConfigs
    elseif noex > 0
        newConfigs = Basics.generateConfigurations(refConfigs, fromShells, toShells)
        newConfigs = Base.unique(newConfigs)
        for ne = 2:noex    
            newConfigs = Basics.generateConfigurations(newConfigs, fromShells, toShells)    
            newConfigs = Base.unique(newConfigs)
        end
    end
    
    # Now apply in turn all given restrictions, if any, and append if no restriction is violated
    confList = Configuration[]
    for  conf  in  newConfigs
        addConf = true
        for res in restrictions
            if  Basics.isViolated(conf,res)     addConf = false;    break   end
        end
        if  addConf     push!(confList, conf)   end
    end
    
    return( confList )
end



"""
`Basics.generateConfigurationsForExcitationScheme(confs::Array{Configuration,1}, exScheme::Basics.NoExcitationScheme, 
                                                    nMax::Int64, lValues::Array{Int64,1})`  
    ... generates a list of non-relativistic configurations for the given (reference) confs and the excitation scheme. 
        All orbitals in standard order are considered up to maxShell = (n_max, l_max).
"""
function Basics.generateConfigurationsForExcitationScheme(confs::Array{Configuration,1}, exScheme::Basics.NoExcitationScheme, 
                                                        nMax::Int64, lValues::Array{Int64,1})
    @warn(  "No excitations are included if exScheme::NoExcitationScheme. ")
    println("No excitations are included if exScheme::NoExcitationScheme. ")
    #
    newConfList = confs
    return( newConfList )
end



"""
`Basics.generateConfigurationsForExcitationScheme(confs::Array{Configuration,1}, exScheme::Basics.DeExciteSingleElectron, 
                                                    nMax::Int64, lValues::Array{Int64,1})`  
    ... generates a list of non-relativistic configurations for the given (reference) confs and the excitation scheme. 
        All orbitals in standard order are considered up to maxShell = (n_max, l_max).
"""
function Basics.generateConfigurationsForExcitationScheme(confs::Array{Configuration,1}, exScheme::Basics.DeExciteSingleElectron, 
                                                            nMax::Int64, lValues::Array{Int64,1})
    confList  = Configuration[];    NoElectrons = confs[1].NoElectrons
    shellList = Basics.generateShellList(confs, nMax, lValues)
    
    # Create all configurations with single excitations from the given list of configurations
    for  conf in confs
        # Take one electron fromShell and add toShell
        for  (fromShell,occ)  in  conf.shells
            for  toShell  in  shellList
                newShells = deepcopy( conf.shells )
                if      fromShell == toShell    
                elseif  haskey(conf.shells, toShell )   fromOcc = occ;   toOcc = conf.shells[toShell]
                        if  fromOcc - 1 < 0                        println("..");    continue    end
                        if  toOcc   + 1 > 2*(2*toShell.l + 1)      println("..");    continue    end
                        newShells[fromShell] = newShells[fromShell] - 1
                        newShells[toShell]   = newShells[toShell] + 1
                else    fromOcc = occ
                        if  fromOcc - 1 < 0                        println("..");    continue    end
                        newShells[fromShell] = newShells[fromShell] - 1
                        newShells = Base.merge( newShells, Dict( toShell => 1))
                end
                # Add a new configuration
                push!( confList, Configuration( newShells, NoElectrons))
                if  true  println(">> Generate $(Configuration( newShells, NoElectrons))")   end
            end
        end
    end

    nbefore     = length(confList)
    confList = unique(confList)
    ## newConfList = Configuration[]
    ## for  conf  in  confList
    ##     if  conf in newConfList    continue;    else    push!( newConfList,  conf)    end
    ## end
    nafter      = length(confList)
    println(">> Number of generated configurations for $exScheme is: $nbefore (before) and $nafter (after).")

    return( confList )
end



"""
`Basics.generateConfigurationsForExcitationScheme(confs::Array{Configuration,1}, exScheme::Basics.DeExciteTwoElectrons, 
                                                    nMax::Int64, lValues::Array{Int64,1})`  
    ... generates a list of non-relativistic configurations for the given (reference) confs and the excitation scheme. 
        All orbitals in standard order are considered up to maxShell = (n_max, l_max).
"""
function Basics.generateConfigurationsForExcitationScheme(confs::Array{Configuration,1}, exScheme::Basics.DeExciteTwoElectrons, 
                                                            nMax::Int64, lValues::Array{Int64,1})
    newConfList = Basics.generateConfigurationsForExcitationScheme(confs, Basics.DeExciteSingleElectron(), nMax, lValues)
    newConfList = Basics.generateConfigurationsForExcitationScheme(newConfList, Basics.DeExciteSingleElectron(), nMax, lValues)
    newConfList = unique(newConfList)
    return( newConfList )
end



"""
`Basics.generateConfigurationsForExcitationScheme(confs::Array{Configuration,1}, exScheme::Basics.ExciteByCapture, 
        fromShells::Array{Shell,1}, toShells::Array{Shell,1}, intoShells::Array{Shell,1}, noex::Int64)`  
    ... generates a list of non-relativistic configurations for the given (reference) confs with the excitation
        of upto noex electrons fromShell --> toShells and the additional capture of one electron --> intoShells.
        All configuration will therefore contain N+1 electrons.
"""
function Basics.generateConfigurationsForExcitationScheme(confs::Array{Configuration,1}, exScheme::Basics.ExciteByCapture, 
                fromShells::Array{Shell,1}, toShells::Array{Shell,1}, intoShells::Array{Shell,1}, noex::Int64)
    NoElectrons = confs[1].NoElectrons
    confList    = Basics.generateConfigurations(confs, fromShells, toShells, noex)
    shellList   = Basics.extractNonrelativisticShellList(confList)
    #
    # Now add one electron to the intoShells
    newConfList = Configuration[]
    for  conf in confList
        for intoShell in intoShells
            newShells = deepcopy( conf.shells )
            if   haskey(conf.shells, intoShell )   occ = conf.shells[intoShell]
                if  occ   + 1 > 2*(2*intoShell.l + 1)      println("..");    continue    end
                newShells[intoShell] = newShells[intoShell] + 1
            else  
                newShells = Base.merge( newShells, Dict( intoShell => 1))
            end
            # Add a new configuration
            push!( newConfList, Configuration( newShells, NoElectrons+1))
            if  true  println(">> Generate $(Configuration( newShells, NoElectrons+1))")   end
        end
    end
    newConfList = unique(newConfList)
    
    nbefore     = length(confList)
    nafter      = length(newConfList)
    println(">> Number of generated configurations for $exScheme is: $nbefore (before) and $nafter (after).")

    return( newConfList )
end


"""
`Basics.generateFieldCoordinates(mesh::Basics.AbstractMesh)`  
    ... generates a list of field values of proper type due to the specification by the given mesh; 
        this specification determines both, the number and kind of coordinates as well as the type of the array.
        A list::Array{...FieldValue{Type},1} is returned.
"""
function Basics.generateFieldCoordinates(mesh::Basics.AbstractMesh)
    #
    if      typeof(mesh) == Basics.Cartesian2DMesh
        fValues = Basics.Cartesian2DFieldValue{Float64}[]
        xs      = Basics.generateMeshCoordinates(mesh.xMesh)
        ys      = Basics.generateMeshCoordinates(mesh.yMesh)
        for  x in xs, y in ys             push!(fValues, Basics.Cartesian2DFieldValue{Float64}(x,y, 0.))     end
    elseif  typeof(mesh) == Basics.PolarMesh 
        fValues = Basics.PolarFieldValue{Float64}[]
        rhos    = Basics.generateMeshCoordinates(mesh.rhoMesh)
        phis    = Basics.generateMeshCoordinates(mesh.phiMesh)
        for  rho in rhos, phi in phis     push!(fValues, Basics.PolarFieldValue{Float64}(rho,phi, 0.))       end
    elseif  typeof(mesh) == Basics.SphericalMesh
        fValues = Basics.SphericalFieldValue{Float64}[]
        rs      = Basics.generateMeshCoordinates(mesh.rMesh)
        thetas  = Basics.generateMeshCoordinates(mesh.thetaMesh)
        phis    = Basics.generateMeshCoordinates(mesh.phiMesh)
        for  r in rs, theta in thetas, phi in phis     push!(fValues, Basics.SphericalFieldValue{Float64}(r,theta,phi, 0.))    end
    else    error("Unknown mesh type: typeof(mesh) = $(typeof(mesh))")
    end 
    
    return( fValues )
end
        


"""
`Basics.generateLevelWithExtraElectron(newOrbital, symt::LevelSymmetry, level::Level)`  
    ... generates a (new) level with one extra electron in subshell newOrbital.subshell and with overall symmetry symt. The function 
        assumes that all CSF in the basis of level have the same symmetry as the level itself, that the new subshell is not yet part 
        of the basis and that CSF with the total symmetry can be constructed; it terminates with an error if one of this assumptions
        is violated. A newLevel::Level with the same number of CSF and the same representation (eigenvector) as the given level is 
        returned. From a physics viewpoint, newLevel refers to a antisymmetrized product state of [level x orbital(sh, energy)] J^P 
        (symt) and with an total energy = level.energy + newOrbital.energy.
"""
function Basics.generateLevelWithExtraElectron(newOrbital, symt::LevelSymmetry, level::Level)
    basis = level.basis;    newSubshells = copy(basis.subshells);    newCsfs = CsfR[];   J = level.J;   parity = level.parity
    
    push!(newSubshells, newOrbital.subshell)
    Defaults.setDefaults("relativistic subshell list", newSubshells; printout=false)
    newOrbitals = Base.merge( copy(basis.orbitals), Dict( newOrbital.subshell => newOrbital ))

    for  i = 1:length(basis.csfs)
        stateList   = ManyElectron.provideSubshellStates(newOrbital.subshell, 1)
        if  basis.csfs[i].J != J   ||   basis.csfs[i].parity != parity    error("Improper symmetry of CSF.")             end
        if  length(stateList) != 1                                        error("Improper number of subshell states.")   end
        if  AngularMomentum.triangularDelta(J, Basics.subshell_j(newOrbital.subshell), symt.J) != 1    
                                                                        error("Improper coupling of (new) subshell & total J.")   end
        substate   = stateList[1]
        occupation = copy(basis.csfs[i].occupation);    push!(occupation, substate.occ)
        seniorityNr  = copy(basis.csfs[i].seniorityNr);     push!(seniorityNr,  substate.nu)
        subshellJ  = copy(basis.csfs[i].subshellJ);     push!(subshellJ,  AngularJ64(substate.Jsub2//2) )
        subshellX  = copy(basis.csfs[i].subshellX);     push!(subshellX,  symt.J)
        push!(newCsfs, CsfR(true, symt.J, symt.parity, occupation, seniorityNr, subshellJ, subshellX, Subshell[] ) )
    end

    newBasis = Basis(true, basis.NoElectrons+1, newSubshells, newCsfs, copy(basis.coreSubshells), newOrbitals)
    newLevel = Level(symt.J, AngularM64(symt.J), symt.parity, level.index, level.energy+newOrbital.energy, level.relativeOcc, 
                    true, newBasis, level.mc)
    return( newLevel )
end


"""
`Basics.generateLevelWithExtraTwoElectrons(orb1::Orbital, symx::LevelSymmetry, orb2::Orbital, symt::LevelSymmetry, level::Level)`  
    ... generates a (new) level with two extra electrons in subshells orb1.subshell and orb2.subshell and with intermediate
        symmetry symx and overall symmetry symt. The function assumes that all CSF in the basis of level have the 
        same symmetry as the level itself, that the new subshells are not yet part of the basis and that CSF with the 
        intermediate and total symmetries can be constructed; it terminates with an error if one of this assumptions
        is violated. A newLevel::Level with the same number of CSF and the same representation (eigenvector) as the given level is 
        returned. From a physics viewpoint, newLevel refers to a antisymmetrized product state of 
        [ [level x orbital1(sh, energy)] Jx^Px  x orbital1(sh, energy)] J^P  (symt) and with an 
        total energy = level.energy + orb1.energy + orb2.energy.
"""
function Basics.generateLevelWithExtraTwoElectrons(orb1::Orbital, symx::LevelSymmetry, orb2::Orbital, symt::LevelSymmetry, level::Level)
    basis = level.basis;    newSubshells = copy(basis.subshells);    newCsfs = CsfR[];   J = level.J;   parity = level.parity
    
    push!(newSubshells, orb1.subshell);     push!(newSubshells, orb2.subshell)
    Defaults.setDefaults("relativistic subshell list", newSubshells; printout=false)
    newOrbitals = Base.merge( copy(basis.orbitals), Dict( orb1.subshell => orb1, orb2.subshell => orb2 ))

    for  i = 1:length(basis.csfs)
        stateList1   = ManyElectron.provideSubshellStates(orb1.subshell, 1)
        stateList2   = ManyElectron.provideSubshellStates(orb2.subshell, 1)
        if  basis.csfs[i].J != J     ||   basis.csfs[i].parity != parity    error("Improper symmetry of CSF.")             end
        if  length(stateList1) != length(stateList2) != 1                   error("Improper number of subshell states.")   end
        if  AngularMomentum.triangularDelta(J,      Basics.subshell_j(orb1.subshell), symx.J) != 1    
                                                            error("Improper coupling of (new) subshell & intermediate J.")  end
        if  AngularMomentum.triangularDelta(symx.J, Basics.subshell_j(orb2.subshell), symt.J) != 1    
                                                            error("Improper coupling of (new) subshell & total J.")         end
        substate1   = stateList1[1];    substate2  = stateList2[1]
        occupation  = copy(basis.csfs[i].occupation);    push!(occupation, substate1.occ);    push!(occupation, substate2.occ)
        seniorityNr = copy(basis.csfs[i].seniorityNr);   push!(seniorityNr,substate1.nu);     push!(seniorityNr,substate2.nu)
        subshellJ   = copy(basis.csfs[i].subshellJ);     push!(subshellJ,  AngularJ64(substate1.Jsub2//2) )
                                                            push!(subshellJ,  AngularJ64(substate2.Jsub2//2) )
        subshellX   = copy(basis.csfs[i].subshellX);     push!(subshellX,  symx.J);           push!(subshellX,  symt.J)
        push!(newCsfs, CsfR(true, symt.J, symt.parity, occupation, seniorityNr, subshellJ, subshellX, Subshell[] ) )
    end

    newBasis = Basis(true, basis.NoElectrons+2, newSubshells, newCsfs, copy(basis.coreSubshells), newOrbitals)
    newLevel = Level(symt.J, AngularM64(symt.J), symt.parity, level.index, level.energy+orb1.energy+orb2.energy, level.relativeOcc, 
                    true, newBasis, level.mc)
    return( newLevel )
end


"""
`Basics.generateLevelWithExtraSubshell(sh::Subshell, level::Level)`  
    ... generates a (new) level with one extra subshell sh but with the same overall symmetry as before. The function assumes that the 
        new subshell is not yet part of the basis; it terminates with an error if one of this assumptions is violated. 
        A newLevel::Level with the same symmetry, number of CSF, energy and the same representation (eigenvector) as the given 
        level is returned. 
"""
function Basics.generateLevelWithExtraSubshell(sh::Subshell, level::Level)
    basis = level.basis;    newSubshells = copy(basis.subshells);    newCsfs = CsfR[];   J = level.J;   parity = level.parity
    newOrbitals = copy(basis.orbitals)

    push!(newSubshells, sh)
    newOrbitals = Base.merge( newOrbitals, Dict( sh => Orbital(sh, -10000.) ))
    Defaults.setDefaults("relativistic subshell list", newSubshells; printout=false)

    for  i = 1:length(basis.csfs)
        stateList   = ManyElectron.provideSubshellStates(sh, 0)
        if  length(stateList) != 1                                        error("Improper number of subshell states.")   end
        substate   = stateList[1]
        occupation = copy(basis.csfs[i].occupation);    push!(occupation, substate.occ)
        seniorityNr  = copy(basis.csfs[i].seniorityNr);     push!(seniorityNr,  substate.nu)
        subshellJ  = copy(basis.csfs[i].subshellJ);     push!(subshellJ,  AngularJ64(substate.Jsub2//2) )
        subshellX  = copy(basis.csfs[i].subshellX);     push!(subshellX,  subshellX[end])
        push!(newCsfs, CsfR(true, basis.csfs[i].J, basis.csfs[i].parity, occupation, seniorityNr, subshellJ, subshellX, Subshell[] ) )
    end

    newBasis = Basis(true, basis.NoElectrons, newSubshells, newCsfs, copy(basis.coreSubshells), newOrbitals)
    newLevel = Level(level.J, level.M, level.parity, level.index, level.energy, level.relativeOcc, true, newBasis, level.mc)
    return( newLevel )
end


"""
`Basics.generateLevelWithExtraSubshells(subshells::Array{Subshell,1}, level::Level)`  
    ... generates a (new) level with one or several extra subshells but with the same overall symmetry as before. 
        The function assumes that the new subshell is not yet part of the basis. A newLevel::Level with the same symmetry, 
        number of CSF, energy and the same representation (eigenvector) as the given level is returned. 
"""
function Basics.generateLevelWithExtraSubshells(subshells::Array{Subshell,1}, level::Level)
    newLevel = deepcopy(level)
    for  sh in subshells    newLevel = Basics.generateLevelWithExtraSubshell(sh::Subshell, newLevel::Level)    end
    return( newLevel )
end


"""
`Basics.generateLevelWithExtraSubshell(sh::Subshell, level::Level, ns::Int64)`  
    ... generates a (new) level with one extra subshell sh but with the same overall symmetry as before. The function assumes that the 
        new subshell is not yet part of the basis and terminates if this assumptions is not fulfilled. 
        It builts the subshell sh at position ns and moves all other subshells further.
        A newLevel::Level with the same symmetry, number of CSF, energy and the same representation (eigenvector) as the given 
        level is returned. 
"""
function Basics.generateLevelWithExtraSubshell(sh::Subshell, level::Level, ns::Int64)
    basis = level.basis;    newCsfs = CsfR[];   J = level.J;   parity = level.parity
    newOrbitals = copy(basis.orbitals)
    
    # Build the new subshell list and add the associated orbital
    newSubshells = basis.subshells[1:ns-1];   push!(newSubshells, sh);    append!(newSubshells, basis.subshells[ns:end])
    newOrbitals = Base.merge( newOrbitals, Dict( sh => Orbital(sh, -10000.) ))
    Defaults.setDefaults("relativistic subshell list", newSubshells; printout=false)

    for  i = 1:length(basis.csfs)
        stateList   = ManyElectron.provideSubshellStates(sh, 0)
        if    length(stateList) != 1        error("Improper number of subshell states.")   
        else                                substate = stateList[1]
        end
        
        occupation  = copy(basis.csfs[i].occupation[1:ns-1]);  push!(occupation, substate.occ);    append!(occupation, basis.csfs[i].occupation[ns:end]) 
        seniorityNr = copy(basis.csfs[i].seniorityNr[1:ns-1]); push!(seniorityNr,  substate.nu);   append!(seniorityNr, basis.csfs[i].seniorityNr[ns:end]) 
        subshellJ   = copy(basis.csfs[i].subshellJ[1:ns-1]);   push!(subshellJ,  AngularJ64(substate.Jsub2//2) );  append!(subshellJ, basis.csfs[i].subshellJ[ns:end]) 
        if  ns  > 1    subshellX   = copy(basis.csfs[i].subshellX[1:ns-1]);    push!(subshellX,  subshellX[ns-1])
        else           subshellX   = AngularJ64[AngularJ64(0)]                           
        end
        append!(subshellX, basis.csfs[i].subshellX[ns:end])
        push!(newCsfs, CsfR(true, basis.csfs[i].J, basis.csfs[i].parity, occupation, seniorityNr, subshellJ, subshellX, Subshell[] ) )
    end

    newBasis = Basis(true, basis.NoElectrons, newSubshells, newCsfs, copy(basis.coreSubshells), newOrbitals)
    newLevel = Level(level.J, level.M, level.parity, level.index, level.energy, level.relativeOcc, true, newBasis, level.mc)
    return( newLevel )
end


"""
`Basics.generateLevelWithSymmetryReducedBasis(level::Level, subshells::Array{Subshell,1})`  
    ... generates a level with a new basis and representation that only includes the CSF with the same symmetry as the level 
        itself. It also ensures that the new basis is based on the given subshells so that it can be later used with levels
        which have one or several subshells more (at it often appears in ionization processes). Both, the underlying CSF basis 
        and the eigenvectors are adopted accordingly. The procedures only supports the 'addition' of subshells but it stops 
        if (1) the individual sequence differs or if the basis of level contains more subshells than the given subshell-list.
        A (new) level::levelR is returned.
"""
function Basics.generateLevelWithSymmetryReducedBasis(level::Level, subshells::Array{Subshell,1})
    testBasis = level.basis;    newCsfs = CsfR[];    newMc = Float64[]

    # Decide of whether the basis need to be extended by some subshells
    nt = length(testBasis.subshells)
    if      testBasis.subshells == subshells      basis = testBasis;    nLevel = level
    elseif  nt >= length(subshells)               error("stop a")
    elseif  nt + 1 == length(subshells)  &&  testBasis.subshells == subshells[1:nt]
            nLevel = Basics.generateLevelWithExtraSubshell(subshells[end], level)
            basis  = nLevel.basis
    elseif  testBasis.subshells == subshells[1:nt]
        nLevel  = level
        for  ns = nt+1:length(subshells)
            nLevel = Basics.generateLevelWithExtraSubshell(subshells[ns], nLevel)
        end
        basis  = nLevel.basis
    else
        nLevel  = level
        for  (ns, sh) in enumerate(subshells)
            if  sh in testBasis.subshells    continue   
            else      nLevel = Basics.generateLevelWithExtraSubshell(sh, nLevel, ns)
            end
        end
        basis  = nLevel.basis
    end
        
    
    for  i = 1:length(basis.csfs)
        if  basis.csfs[i].J == nLevel.J  &&   basis.csfs[i].parity == nLevel.parity        &&    abs(nLevel.mc[i]) > 1.0e-7
            push!(newCsfs, deepcopy(basis.csfs[i]));    push!(newMc, deepcopy(nLevel.mc[i]))
        end
    end
    newBasis = Basis(true, basis.NoElectrons, deepcopy(basis.subshells), newCsfs, deepcopy(basis.coreSubshells), deepcopy(basis.orbitals))
    newLevel = Level(nLevel.J, nLevel.M, nLevel.parity, nLevel.index, nLevel.energy, nLevel.relativeOcc, true, newBasis, newMc)
    return( newLevel )
end


"""
`Basics.generateMeshCoordinates(mesh::Basics.AbstractMesh)`  
    ... generates a list of coordinates for the given  (1D)mesh; a list::Array{Float64,1} is returned.
        Note: No 'weights' for integration are returned in the present form ... though this could be done readily.
"""
## function Basics.generateMeshCoordinates(mesh::Basics.GLegenreMesh)
function Basics.generateMeshCoordinates(mesh::Basics.AbstractMesh)
    #
    ## using  QuadGK
    if      typeof(mesh) == Basics.GLegenreMesh
        # The weights of the GL coordinates are determined here ... but not returned in the present version.
        wax  = QuadGK.gauss(mesh.NoZeros);    t = wax[1];     wt = wax[2]        
        fac1 = 0.5 * ( mesh.b - mesh.a );   fac2 = 0.5 * ( mesh.b + mesh.a )
        for  j = 1:mesh.NoZeros
            t[j] = fac1 * t[j] + fac2;    wt[j] = fac1 * wt[j]
        end
        coords = t
    elseif  typeof(mesh) == Basics.LinearMesh 
        coords = Float64[];    x0 = (mesh.b - mesh.a) / (mesh.NoPoints-1);   x = mesh.a
        for  i = 1:mesh.NoPoints    x = x + x0;    push!(coords, x)     end
    else    error("Unknown mesh type for generating coordinates: typeof(mesh) = $(typeof(mesh))")
    end 
    
    return( coords )
end


"""
`Basics.generateOrbitalsForBasis(basis::Basis, frozenShells::Array{Shell,1}, priorBasis::Basis, startOrbitals::Dict{Subshell, Orbital})`  
    ... generates a dict of (relativistic) orbitals as specificed by subshells of the basis. These orbitals are taken from priorBasis
        if contained in frozen-shells, and from spectrum (start orbitals) otherwises. An error message is issued if some
        requested orbital is not found. A orbList::Dict{Subshell, Orbital} is returned.
"""
function Basics.generateOrbitalsForBasis(basis::Basis, frozenShells::Array{Shell,1}, priorBasis::Basis, 
                                            startOrbitals::Dict{Subshell, Orbital})
    orbitals = Dict{Subshell, Orbital}()
    for  subsh in basis.subshells
        if  subsh in frozenShells    &&   haskey(priorBasis.orbitals, subsh)
            orbitals = Base.merge( orbitals, Dict( subsh => priorBasis.orbitals[subsh]) )
        elseif  subsh in frozenShells     error("Frozen orbital $subsh not found in prior basis.")
        elseif  haskey(priorBasis.orbitals, subsh)
            println(">> Start orbital $subsh is taken from prior basis")
            orbitals = Base.merge( orbitals, Dict( subsh => priorBasis.orbitals[subsh]) )
        else
            println(">> Start orbital $subsh is taken from hydrogenic orbitals")
            orbitals = Base.merge( orbitals, Dict( subsh => startOrbitals[subsh]) )
        end
    end
    return( orbitals )
end



"""
`Basics.generateOrbitalSuperposition(a::Orbital, b::Orbital, cx::Float64, grid::Radial.Grid)`  
    ... generates a superposition  a + cx * b of two given orbitals; the function just takes the linear combination
        of the large and small components as well as the energy. The function assumes that both orbitals are defined on the same
        grid. A re-normalized newOrbital::Orbital is returned for which the derivatives are not defined. 
        This function is used to accelerate the convergence (hopefully).
"""
function Basics.generateOrbitalSuperposition(a::Orbital, b::Orbital, cx::Float64, grid::Radial.Grid)
    if  a.subshell != b.subshell  ||  !(a.isBound)  ||  !(b.isBound)   error("stop a")     end
    mtp  = max(size(a.P, 1), size(b.P, 1));     newP = zeros(mtp);  newQ = zeros(mtp);   newEnergy = (a.energy + cx * b.energy) / (1 + cx)
    mtpa = size(a.P, 1);    newP[1:mtpa] = a.P[1:mtpa];    newQ[1:mtpa] = a.Q[1:mtpa]
    mtpb = size(b.P, 1);    newP[1:mtpb] = newP[1:mtpb] + cx * b.P[1:mtpb]    
                            newQ[1:mtpb] = newQ[1:mtpb] + cx * b.Q[1:mtpb]
    newOrbital = Orbital(a.subshell, a.isBound, a.useStandardGrid, newEnergy, newP, newQ, Float64[], Float64[], a.grid)
    norm       = RadialIntegrals.overlap(newOrbital, newOrbital, grid)
    newP       = newP / sqrt(norm);      newQ = newQ / sqrt(norm)
    newOrbital = Orbital(a.subshell, a.isBound, a.useStandardGrid, newEnergy, newP, newQ, Float64[], Float64[], a.grid)

    return( newOrbital )
end


"""
`Basics.generateShellList(confs::Array{Configuration,1}, nMax::Int64, lValues::Array{Int64,1})`  
    ... generates a list of (non-relativistic) shells in standard order that contains all shells from the given list of
        configurations as well as the shells up to the principal quantum number nMax and the list of orbital quantum
        numbers lValues. A shellList::Array{Shell,1} is returned.
"""
function Basics.generateShellList(confs::Array{Configuration,1}, nMax::Int64, lValues::Array{Int64,1})
    shellList = Shell[]
    # Add all shells from the given configurations
    for  conf in confs
        for (k,v) in conf.shells    push!( shellList, k)    end
    end
    # Add all shells for the given quantum numbers
    for  n = 1:nMax
        for  ll in  lValues  
            if  n >= ll + 1     push!( shellList, Shell(n,ll))      end    
        end
    end
    
    # Now bring the shells in standard order
    newShellList = Shell[]
    for  n = 1:nMax
        for  l = 0:30  
            if  n >= l + 1  && Shell(n,l) in shellList    push!( newShellList, Shell(n,l))     end    
        end
    end
    println(">> From configurations generated shell list $newShellList ")
    
    return( newShellList )
end


"""
`Basics.generateShellList(nMin::Int64, nMax::Int64, lValues::Array{Int64,1})`  
    ... generates a list of (non-relativistic) shells in standard order that contains all shells with principal quantum
        numbers from nMin .. nMaxthe and orbital angular momenta from lValues. A shellList::Array{Shell,1} is returned.
"""
function Basics.generateShellList(nMin::Int64, nMax::Int64, lValues::Array{Int64,1})
    shellList = Shell[]
    # Add all shells for the given quantum numbers
    for  n = nMin:nMax
        for  ll in  lValues  
            if  n >= ll + 1     push!( shellList, Shell(n,ll))      end    
        end
    end
    
    # Now bring the shells in standard order
    newShellList = Shell[]
    for  n = 1:nMax
        for  l = 0:30  
            if  n >= l + 1  && Shell(n,l) in shellList    push!( newShellList, Shell(n,l))     end    
        end
    end
    println(">> Generated shell list $newShellList ")
    
    return( newShellList )
end


"""
`Basics.generateShellList(nMin::Int64, nMax::Int64, lMax::Int64)`  
    ... generates a list of (non-relativistic) shells in standard order that contains all shells with principal quantum
        numbers from nMin .. nMaxthe and orbital angular momenta l <= lMax. A shellList::Array{Shell,1} is returned.
"""
function Basics.generateShellList(nMin::Int64, nMax::Int64, lMax::Int64)
    shellList = Shell[]
    # Add all shells for the given quantum numbers
    for  n = nMin:nMax
        for  ll = 0:lMax 
            if  n >= ll + 1     push!( shellList, Shell(n,ll))      end    
        end
    end
    
    # Now bring the shells in standard order
    newShellList = Shell[]
    for  n = 1:nMax
        for  l = 0:30  
            if  n >= l + 1  && Shell(n,l) in shellList    push!( newShellList, Shell(n,l))     end    
        end
    end
    ## println(">> Generated shell list $newShellList ")
    
    return( newShellList )
end


"""
`Basics.generateShellList(nMin::Int64, nMax::Int64, symMax::String)`  
    ... generates a list of (non-relativistic) shells in standard order that contains all shells with principal quantum
        numbers from nMin .. nMaxthe and orbital angular momenta l <= l(symMax), and where symMax is the symmetry character
        string ("s", "p", ...,"z") of the corresponding orbital angular momentum. A shellList::Array{Shell,1} is returned.
"""
function Basics.generateShellList(nMin::Int64, nMax::Int64, symMax::String)
    lMax = Basics.shellNotation(symMax)
    shellList = Shell[]
    # Add all shells for the given quantum numbers
    for  n = nMin:nMax
        for  ll = 0:lMax 
            if  n >= ll + 1     push!( shellList, Shell(n,ll))      end    
        end
    end
    
    # Now bring the shells in standard order
    newShellList = Shell[]
    for  n = 1:nMax
        for  l = 0:30  
            if  n >= l + 1  && Shell(n,l) in shellList    push!( newShellList, Shell(n,l))     end    
        end
    end
    println(">> Generated shell list $newShellList ")
    
    return( newShellList )
end


"""
`Basics.generateSubshellList(shells::Array{Shell,1})`  
    ... generates a list of relativistic subshells from the given shell list, and by keeping the same order.
        A subshellList::Array{Subshell,1} is returned.
"""
function Basics.generateSubshellList(shells::Array{Shell,1})
    subshellList = Subshell[]
    for  sh in shells
        if  sh.l == 0   push!(subshellList, Subshell(sh.n, -1))
        else            push!(subshellList, Subshell(sh.n, sh.l))
                        push!(subshellList, Subshell(sh.n, -sh.l - 1))
        end
    end
        
    return( subshellList )
end


"""
`Basics.generateSpectrumLorentzian(xIntensities::Array{Tuple{Float64,Float64},1}, widths::Float64; 
                                    energyShift::Float64=0., resolution::Int64=500)`  
    ... to generate the Lorentzian spectrum for the given intensities (x-position, y-intensity) and the (constant) width.
        A tuple of an xList::Array{Float64,1} (positions) and yList::Array{Float64,1} (valiues) is returned that can immediately
        be plotted. The energy shift is provided already in 'current' units and is simply added to all energies.
"""
function Basics.generateSpectrumLorentzian(xIntensities::Array{Tuple{Float64,Float64},1}, widths::Float64; 
                                            energyShift::Float64=0., resolution::Int64=500)
    function normalLorentzian(en::Float64, ga::Float64)
        gaover2 = ga /2.0;  wa = gaover2 / (en^2 + gaover2^2);      wa = wa / pi
    end
    #
    # Convert energies and determine their range 
    wc = Defaults.convertUnits("energy: from atomic", 1.0);  xwidths = widths * wc
    energies = Float64[];     intensities = Float64[]
    for xInt in xIntensities    push!(energies, xInt[1]*wc + energyShift);    push!(intensities, xInt[2])   end
    xMin = minimum(energies) - 3.0xwidths
    xMax = maximum(energies) + 3.0xwidths
    @show xMin, xMax
    res = 1 / resolution;   nn = resolution * Int64( floor(xMax-xMin) + 1.0 )
    x = zeros(nn);    y = zeros(nn)
    for  n = 1:nn     x[n] = xMin +n*res    end
    for  n = 1:nn 
        for ie = 1:length(energies)
            en = energies[ie];   int = intensities[ie]
            y[n] = y[n] + int * normalLorentzian(x[n] - en, xwidths)
        end
    end
    
    return(x,y)
end


