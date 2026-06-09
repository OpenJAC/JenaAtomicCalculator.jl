
##x export  checkConfigurations,  displayConfiguration,  displayConfigurations,  extractConfiguration,  extractConfigurations,  
##x         extractFromConfiguration,  extractFromConfigurations,  generateConfiguration,  generateConfigurations


"""
`Basics.checkConfigurations(theme::Basics.NumberOfElectrons, confs::Array{Configuration,1})`  
    ... checks for and terminates with a comment, if not all the given configurations have the same number of 
        electrons > 0.
"""
function  Basics.checkConfigurations(theme::Basics.NumberOfElectrons, confs::Array{Configuration,1})
    was   = Basics.extractFromConfigurations(theme, confs)
    wcs   = unique(was)
    freqs = Dict(x => count(==(x), was) for x in unique(was))
    wbs   = Int64[];   for (k,v) in freqs   push!(wbs, v)   end;    mbs = maximum(wbs) 
    wa    = 0;         for (k,v) in freqs   if  freqs[k] == mbs     wa  = k;   break   end   end
    # wa ... is the most frequent number of electrons
    if  length(wcs) != 1
        println("Configurations found with a different number of electron; problems may arise for:")
        for  i = 1:length(was)
            if  was[i] != wa    println("> $(confs[i])")   end 
        end 
        error("\n ... the executation terminates. ")
    end
    
    return( nothing )
end   

    
"""
`Basics.checkConfigurations(theme::Basics.NumberOfElectrons, confs::Array{Configuration,1}, 
                            confsMinus1::Array{Configuration,1})`  
    ... checks for and terminates with a comment, if not all the given configurations have the same number of 
        electrons > 0.
"""
function  Basics.checkConfigurations(theme::Basics.NumberOfElectrons, confs::Array{Configuration,1}, 
                                     confsMinus1::Array{Configuration,1})
    Basics.checkConfigurations(Basics.NumberOfElectrons(), confs)
    Basics.checkConfigurations(Basics.NumberOfElectrons(), confsMinus1)
    #
    was = Basics.extractFromConfigurations(theme, confs);         was = unique(was)
    wbs = Basics.extractFromConfigurations(theme, confsMinus1);   wbs = unique(wbs)

    if  was[1] != wbs[1] + 1
        println("Configurations number of electrons $(was[1]) != $(wbs[1]) + 1 are given; " *
                "check for consistency with the given process.")
        error("\n ... the executation terminates. ")
    end
    
    return( nothing )
end   


#################################################################################################################################
#################################################################################################################################

"""
`Basics.displayConfiguration(stream::IO, theme::Basics.FineStructure, config::Configuration; 
                             details::String="", header::Bool=true)`  
    ... displays the given configuration along with its fine-structure, i.e. the total J's and their number of levels. 
        Nothing is returned in this case.
"""
function  Basics.displayConfiguration(stream::IO, theme::Basics.FineStructure, config::Configuration; 
                                      details::String="", header::Bool=true)
    if  header
        if  details == ""   println(stream,"  Fine-structure of configuration:")
        else                println(stream,"  Fine-structure of configuration for $details:") 
        end
    end
    
    totalJs = Basics.extractFromConfiguration(TotalAM(true, AngularJ64[]), config)
    tJs     = Basics.extractFromConfiguration(TotalAM(false, AngularJ64[]), config);   tJs = sort(tJs)
    counts  = Dict{AngularJ64, Int}()
    for totalJ in totalJs   counts[totalJ] = get(counts, totalJ, 0) + 1   end
    
    sa = "    " * string(config) * "  J's = "
    for  tJ in tJs   sa = sa * string(tJ) * " ($(counts[tJ]))  "   end
    println(stream, sa)

    return( nothing )
end


"""
`Basics.displayConfiguration(stream::IO, theme::Basics.FineStructureLS, config::Configuration; 
                             details::String="", header::Bool=true)`  
    ... displays the given configuration along with its LSJ-fine-structure, i.e. the total ^(2S+1) L, J's and 
        the number of these LSJ-levels. Nothing is returned in this case.
"""
function  Basics.displayConfiguration(stream::IO, theme::Basics.FineStructureLS, config::Configuration; 
                                      details::String="", header::Bool=true)
    if  header
        if  details == ""   println(stream,"  LSJ-coupled fine-structure of configuration:")
        else                println(stream,"  LSJ-coupled fine-structure of configuration for $details:") 
        end
    end 
    
    newTerms = Tuple{AngularJ64, AngularJ64}[];   currentTerms = [(AngularJ64(0), AngularJ64(0))]
    for  (shell, occ) in config.shells
        # Collect all LS-terms from the open shells
        newTerms = Tuple{AngularJ64, AngularJ64}[]
        was   = LSjj.provideShellStates(shell, occ)
        for  wa in was,   (L,S) in currentTerms
            newLs = Basics.oplus( L, AngularJ64(wa.LL//2));    newSs = Basics.oplus( S, AngularJ64(wa.SS//2))
            for  LL in newLs, SS in newSs   push!(newTerms, (LL, SS))   end
        end 
        currentTerms = copy(newTerms)
    end 
    
    newTerms = sort(newTerms)
    nTerms   = unique(newTerms)
    counts   = Dict{Tuple{AngularJ64, AngularJ64}, Int}()
    for newTerm in newTerms     counts[newTerm] = get(counts, newTerm, 0) + 1   end
    
    sa = "  " * string(config) * "  ^(2S+1) L's = "
    for  nTerm in nTerms
        multiplicity = Basics.twice(nTerm[2]) + 1
        if       nTerm[1] ==  AngularJ64(0)   sb = "S"
        elseif   nTerm[1] ==  AngularJ64(1)   sb = "P"
        elseif   nTerm[1] ==  AngularJ64(2)   sb = "D"
        elseif   nTerm[1] ==  AngularJ64(3)   sb = "F"
        elseif   nTerm[1] ==  AngularJ64(4)   sb = "G"
        elseif   nTerm[1] ==  AngularJ64(5)   sb = "H"
        elseif   nTerm[1] ==  AngularJ64(6)   sb = "I"
        elseif   nTerm[1] ==  AngularJ64(7)   sb = "J"
        elseif   nTerm[1] ==  AngularJ64(8)   sb = "K"
        elseif   nTerm[1] ==  AngularJ64(9)   sb = "L"
        elseif   nTerm[1] ==  AngularJ64(10)  sb = "M"
        elseif   nTerm[1] ==  AngularJ64(11)  sb = "N"
        elseif   nTerm[1] ==  AngularJ64(12)  sb = "O"
        elseif   nTerm[1] ==  AngularJ64(13)  sb = "Q"
        elseif   nTerm[1] ==  AngularJ64(14)  sb = "R"
        else                                  sb = "XX"
        end
        sa = sa * string("^$multiplicity$sb") * " ($(counts[nTerm]))  "   
    end
    println(stream, sa)

    return( nothing )
end

"""
`Basics.displayConfiguration(stream::IO, theme::Basics.HyperfineStructure, config::Configuration; 
                             details::String="", header::Bool=true)`  
    ... displays the given configuration along with its hyperfine-structure, i.e. the total F's and their number of 
        levels for the nuclear spin I. Nothing is returned in this case.
"""
function  Basics.displayConfiguration(stream::IO, theme::Basics.HyperfineStructure, config::Configuration; 
                                      details::String="", header::Bool=true)
    if  header
        if  details == ""   println(stream,"  Hyperfine-structure of configuration for nuclear spin $theme.spinI:")
        else                println(stream,"  Hyperfine-structure of configuration for nuclear spin $theme.spinI and $details:") 
        end
    end
    
    totalJs = Basics.extractFromConfiguration(TotalAM(true, AngularJ64[]), config)
    tJs     = unique(totalJs);   tJs = sort(tJs)
    totalFs = AngularJ64[];      for totalJ in totalJs    append!(totalFs, Basics.oplus(theme.spinI, totalJ) )   end
    tFs     = unique(totalFs);   tFs = sort(tFs)

    counts  = Dict{AngularJ64, Int}()
    for totalF in totalFs     counts[totalF] = get(counts, totalF, 0) + 1   end
    
    sa = "    " * string(config) * "  F's = "
    for  tF in tFs   sa = sa * string(tF) * " ($(counts[tF]))  "   end
    println(stream, sa)

    return( nothing )
end


"""
`Basics.displayConfigurations(stream::IO, theme::Basics.HundsRules, config::Configuration; string::String="..")`  
    ... displays ... the generated list of configurations in a compact form; nothing is returned in this case.
    ... with the  occupation numbers
"""
function  Basics.displayConfigurations(stream::IO, theme::Basics.HundsRules, config::Configuration; string::String="..")
    error("Not yet implemented")
    
    return( nothing )
end


"""
`Basics.displayConfiguration(stream::IO, theme::Basics.MeanConfiguration, confs::Array{Configuration,1}; details::String="")`  
    ... displays a mean configuration for all given configurations confs, i.e. an electron configuration with mean occuation 
        numbers. The procedure terminates if configurations with a different number of electrons are given. Nothing is returned.
"""
function  Basics.displayConfiguration(stream::IO, theme::Basics.MeanConfiguration, confs::Array{Configuration,1}; details::String="")
    shellList  = Basics.generateShellList(1, 10, "k");   sa = "\n   +  "
    meanShells = Dict{Shell,Float64}();    
    #
    if  details == ""   println(stream,"  Mean configuration:")
    else                println(stream,"  Mean configurations for $details:") 
    end
    
    if      length(confs) == 0
        println(stream, "\n   + No configurations are given.") 
    else
        for conf in confs
            for  (shell, occ)  in  conf.shells
                if  haskey(meanShells, shell)    meanShells[shell] = meanShells[shell] + occ
                else                             meanShells[shell] = occ
                end 
            end
        end
        #  Divide occupations be number of configurations
        for   (shell, occ)  in  meanShells    meanShells[shell] = meanShells[shell] / length(confs)   end
        # Now printout the configuration with mean occupation numbers
        for  shell in  shellList
            if  haskey(meanShells, shell)
                occ = meanShells[shell]
                sa  = sa * " " * string(shell) * "^(" * @sprintf("%.3f", occ) * ")"
            end
        end
        println(stream, sa)
    end
    
    return( nothing )
end   


#################################################################################################################################
#################################################################################################################################

"""
`Basics.displayConfigurations(stream::IO, configs::Array{Configuration,1}; details::String="", longForm::Bool=false)`  
    ... displays the generated list of configurations in a compact form; nothing is returned in this case.
        If details are given, they are printed as "... configurations for " * details * ":"
        If longForm == true,  the full configurations are displayed with all shells; otherwise, the filled core-shells 
        of all configurations are displayed separately and only the other (valence) shells are displayed explicitly.
        longForm=false is suggested for medium and heavy elements.
"""
function  Basics.displayConfigurations(stream::IO, configs::Array{Configuration,1}; details::String="", longForm::Bool=false)
    # Replace in code: Basics.displayConfigurations(Z::Float64, confs::Array{Configuration,1}; sa::String="")
    nx = 85
    println(stream," ")
    if  details == ""   println(stream,"  Generated configurations:")
    else                println(stream,"  Generated configurations for $details:") 
    end
    
    if      length(configs) == 0
        println(stream, "\n   + No configurations are given.") 
    elseif  longForm  ||  configs[1].NoElectrons < 21   
        # Display the full configurations with all shells
        println(stream, "  ", TableStrings.hLine(nx))
        for  (i, conf)  in  enumerate(configs)
            sa = string("(", i, ")" );    sa = TableStrings.flushright(9, sa, na=3);   sa = sa * string(conf)
            println(stream, sa)
        end 
        println(stream, "  ", TableStrings.hLine(nx))
    else 
        # Display the filled core separately [Closed], common to all configurations, and only the other shells explicitly.
        closedConf = Basics.extractFromConfigurations(Basics.ClosedCore(), configs)
        println(stream, "\n   + [Core] = " * string(closedConf) ) 
        println(stream, "  ", TableStrings.hLine(nx))
        for  (i, conf)  in  enumerate(configs)
            valenceConf = Basics.extractFromConfiguration(Basics.ValenceOccupation(), conf, closedConf)
            sa = string("(", i, ")" );    sa = TableStrings.flushright(9, sa, na=3)
            sa = sa * "[Core] " * string(valenceConf)
            println(stream, sa)
        end 
        println(stream, "  ", TableStrings.hLine(nx))
    end
    
    return( nothing )
end

    
"""
`Basics.displayConfigurations(stream::IO, theme::Basics.ByNumber, configs::Array{Configuration,1}; 
                              details::String="", longForm::Bool=false)`  
    ... displays and summarizes all configurations from the given list, ordered by the numbers of electrons.
        The theme.NoElectrons remains unconsidered in this case. If not specified otherwise, these configurations 
        are displayed in a compact form and nothing is returned in this case.
        If details are given, they are printed as "... configurations for " * details * " and with ne electrons:"
        If longForm,  the full configurations are displayed with all shells; otherwise, the filled core-shells 
        of all configurations are displayed separately and only the other shells are displayed explicitly.
        longForm=false is suggested for medium and heavy elements.
"""
function  Basics.displayConfigurations(stream::IO, theme::Basics.ByNumber, configs::Array{Configuration,1}; 
                                       details::String="", longForm::Bool=false)
    for  no in theme.NoElectrons
        newDetails = details * " and with $no electrons"
        selConfigs = Basics.extractConfigurations(Basics.ByNumber([no]), configs)
        Basics.displayConfigurations(stream, selConfigs, details=newDetails, longForm=longForm)
    end
    
    return( nothing )
end

#################################################################################################################################
#################################################################################################################################


"""
`Basics.displayConfigurationThemes()`  
    ... displays and explains the various  configurationTheme's  <: AbstractConfigurationTheme that are defined to generate,
        deal with and manipulate individual as well as sets of (electron) configurations. These themes provide a simple means 
        in order to extract information and to manipulate such configuration and are, hence, deep within the use of the atomic-shell 
        model. Whereas further details can be extracted for each of these themes, this (dummy) routine aims to summarize the 
        different themes along with its input parameters. Nothing is done (and returned) by this function. The following topics 
        are distinguished:
            
    Add, excite or remove electrons from some reference configurations
    ------------------------------------------------------------------
    + AddElectrons(ne::Int64, intoShells::Array{Shell,1})  
        ... to add to the given reference configurations one or several (ne) electrons into the specified intoShells. 
    + ExciteElectrons(ne::Int64, fromShells::Array{Shell,1}, intoShells::Array{Shell,1})         
        ... to excite for the given configurations one or several (ne) electrons fromShells --> intoShells.
    + RemoveElectrons(ne::Int64, fromShells::Array{Shell,1})          
        ... to remove from the given configurations one or several (ne) electrons from the specified fromShells.
        
    Generate configurations due to different atomic processes
    ---------------------------------------------------------
    + ForAutoIonization()
        ... to generate configurations that are related by autoionization (deexcitation + single remove).
    + ForDielectronicCapture(fromShells::Array{Shell,1}, toShells::Array{Shell,1}, intoShells::Array{Shell,1})
        ... to generate configurations for the dielectronic capture into a given configuration, and by including
            the excitation of an electron fromShells --> toShells as well as the capture of one electron into 
            the intoShells.
    + ForDielectronicRecombination(fromShells::Array{Shell,1}, toShells::Array{Shell,1}, intoShells::Array{Shell,1},
                                   decayShells::Array{Shell,1})
        ... to generate configurations for the dielectronic recombination of ions in given initial configurations:
            It includes the excitation of an electron fromShells --> toShells, the capture of one electron into 
            the intoShells, and the decay of the generated (intermediate) configurations into the decayShells.
    + ForElectronCapture()
        ... to generate configurations for the capture of an electron into a given configuration, i.e. the capture
            of one electron into the intoShells.
    + ForHollowIons()
        ... to generate configurations that are related by multiple capture into high-n shells.
    + ForPhotoEmission()
        ... to generate configurations that are related by photoemission (single de-excitation of an electron).
    + ForPhotoIonization()
        ... to generate configurations that are related by photoionization (single removement of electrons).
    + ForPhotoRecombination(intoShells::Array{Shell,1})
        ... to generate configurations for the capture of an electron into a given configuration, i.e. the 
            radiative recombination of one electron into the intoShells.
    + ForStepwiseDecay(maximallyReleased::Int64)
        ... to generate all configurations that may arise from a given configuration by the release of 
            ne <= maximallyReleased electrons due to photoemission and autoionization.
       
    Generate specific configuration
    -------------------------------    
    + GroundConfiguration(Z::Float64, NoElectrons::Int64)
        ... to generate the ground configuration for an ion with nuclear charge Z and the given number of electrons.
    + MeanConfiguration()
        ... to generate the mean configuration, i.e. a configuration with mean occupation numbers.
    + RelativisticConfigurations()
        ... to generate/deal with relativistic configurations.
    + SuperConfiguration()
        ... to generate configurations that are described by a given super-configuration.
       
    Extract selected configurations from given lists or information about them
    --------------------------------------------------------------------------    
    + ByMultipoles()
        ... to extract the configurations due to certain multipole selection themes (not yet).
    + ByParity(P::Basics.Parity)
        ... to extract the configurations due to the given parity.
    + ClosedCore()
        ... to extract the closed core from one or several given configurations.
    + ClosedShells()
        ... to extract the closed shells from one or several given configurations.
    + ClosedSubshells()
        ... to extract the closed subshells from one or several given configurations.
    + ContractShells()
        ... to extract/contract the shells to the occupied ones.
    + ExcitationLevel()
        ... to extract an excitation level for a given configuration (presently not used).
    + ExpandShells()
        ... to extract/expand the shells to a given shell list.
    + FromBasis()
        ... to extract configuration(s) from a given many-electron basis::Basis.
    + GeneralizedConfigurations()
        .. to extract the generalized configuration from given configurations (not yet).
    + GetParity()
        ... to extract the parity of a -- relativistic or non-relativistic -- configuration.
    + IsOccupied()
        ... to extract of whether a shell or subshell is occupied in some given configuration(s).
    + LeadingConfiguration()
        ... to extract the leading configuration from a given level::Level.
    + LeadingConfigurationR()
        ... to extract the leading relativistic configuration from a given level::Level.
    + Multiplicity()
        ... to extract the multiplicity of a configuration.
    + NonrelativisticBasis()
        ... to extract the configurations from a non-relativistic basis.
    + NumberOfElectrons()
        ... to extract the numbers of electrons from a given set of configurations.
    + OccupationDifference()
        ... to extract differences in the occupation numbers between two given -- relativistic  or
            non-relativistic -- configurations.
    + OpenShells()
        ... to extract the open shells from a given configurations.
    + OpenSubshells()
        ... to extract the open subshells from a given -- relativistic  or non-relativistic -- configurations.
    + TotalAM()
        ... to extract all total angular momenta J, to which a configuration gives rise to.
    + ValenceOccupation()
        ... to extract the occupation of valence-shell electrons from a given configurations with regard to 
            a closed core configuration.
       
    Display configurations with additional information 
    --------------------------------------------------
    + FineStructure()
        ... to display the total J fine-structure levels a configuration (without energies).
    + FineStructureLS()
        ... to display the total LSJ fine-structure levels a configuration (without energies).
    + HundsRules()
        ... to display the total LSJ fine-structure levels, ordered by Hund's rules (not yet).
    + HyperfineStructure()
        ... to display the total F hyperfine-structure levels a configuration (without energies).
        
"""
function  Basics.displayConfigurationThemes()
    return( nothing )
end


#################################################################################################################################
#################################################################################################################################

"""
`Basics.extractConfiguration(theme::Basics.ContractShells, conf::Configuration)`
    ... to extract/contract a configuration by removing all empty shells from the shell-list; a (nonrelativistic)
        conf::Configuration is returned.
"""
function Basics.extractConfiguration(theme::Basics.ContractShells, conf::Configuration)
    shells = Dict{Shell,Int64}()
    for  (shell, occ)  in  conf.shells
        if  occ > 0    shells[shell] = occ   end
    end

    return( Configuration(shells, conf.NoElectrons) )
end


"""
`Basics.extractConfiguration(theme::Basics.ExpandShells, conf::Configuration)`
    ... to extract/expand an (electron) configuration by adding also 'empty' shells; the procedure terminates 
        if shells are occupied in the configuration, which do not occur in theme.shells. A (nonrelativistic)
        conf::Configuration is returned.
"""
function Basics.extractConfiguration(theme::Basics.ExpandShells, conf::Configuration)
    shells = Dict{Shell,Int64}();   noe = 0
    for  shell in theme.shells
        if  haskey(conf.shells, shell)    shells[shell] = conf.shells[shell];   noe = noe + conf.shells[shell]
        else                              shells[shell] = 0
        end
    end
    #
    if  noe != conf.NoElectrons
        error("No all shells $(theme.shells) occur in $conf")
    end

    return( Configuration(shells, conf.NoElectrons) )
end


"""
`Basics.extractConfiguration(theme::Basics.FromBasis, basis::Basis, csf::CsfR)`  
    ... extract the (nonrelativistic) configuration of the given CSF as it is defined in basis.csfs. 
        A conf::Configuration is returned.
"""
function  Basics.extractConfiguration(theme::Basics.FromBasis, basis::Basis, csf::CsfR)
    subshells = basis.subshells
    newShells = Dict{Shell,Int64}();    NoElectrons = 0
    for  i = 1:length(subshells)
        n = subshells[i].n;    l = Basics.subshell_l(subshells[i]);    occ = csf.occupation[i]
        shell = Shell(n,l)
        if    haskey(newShells, shell)    newShells[shell] = newShells[shell] + occ;   NoElectrons = NoElectrons + occ
        else  
            if   occ > 0  newShells = Base.merge( newShells, Dict( shell => occ));     NoElectrons = NoElectrons + occ   end
        end
    end
    
    if      basis.NoElectrons != NoElectrons    error("stop a")
    else    conf = Configuration(newShells, NoElectrons)
    end 
    
    return( conf )
end


"""
`Basics.extractConfiguration(theme::Basics.GroundConfiguration)`  
    ... determines the ground (reference) configuration for an ion with nuclear charge Z and NoElectrons. 
        A refConfig::Configuration is returned.
"""
function Basics.extractConfiguration(theme::Basics.GroundConfiguration)
    # Replace in code:  Plasma.determineReferenceConfiguration(ne::Int64)
    ne = theme.NoElectrons;   refConfig = Configuration("1s")
    
    if theme.Z - theme.NoElectrons >= 4.0  ||  ne <= 18
        # These reference configurations are independent of Z
        if      ne == 1         refConfig = Configuration("1s")
        elseif  ne == 2         refConfig = Configuration("1s^2")
        elseif  ne == 3         refConfig = Configuration("1s^2 2s")
        elseif  ne == 4         refConfig = Configuration("1s^2 2s^2")
        elseif  3 <= ne < 11    refConfig = Configuration("1s^2 2s^2 2p^$(ne-4)")
        elseif  ne == 11        refConfig = Configuration("1s^2 2s^2 2p^6 3s")
        elseif  ne == 12        refConfig = Configuration("1s^2 2s^2 2p^6 3s^2")
        elseif  13 <= ne < 19   refConfig = Configuration("1s^2 2s^2 2p^6 3s^2 3p^$(ne-12)")
        elseif  19 <= ne < 29   refConfig = Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^$(ne-18)")
        elseif  29 <= ne < 31   refConfig = Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^$(ne-28)")
        elseif  31 <= ne < 37   refConfig = Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^$(ne-30)")
        elseif  37 <= ne < 47   refConfig = Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^$(ne-36)")
        elseif  47 <= ne < 61   refConfig = Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 4f^$(ne-46)")
        else    error("stop a")
        end
        
    elseif  0. <= theme.Z - ne <= 0.2
        if      ne == 19        refConfig = Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 4s")
        elseif  ne == 20        refConfig = Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 4s^2")
        elseif  21 <= ne < 27   refConfig = Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^$(ne-20) 4s^2")
        elseif  ne == 28        refConfig = Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^8 4s^2")
        elseif  ne == 29        refConfig = Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s")
        elseif  ne == 29        refConfig = Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2")
        elseif  31 <= ne < 37   refConfig = Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^$(ne-20)")
        elseif  ne == 37        refConfig = Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 5s")
        elseif  ne == 38        refConfig = Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 5s^2")
        elseif  ne == 48        refConfig = Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2")
        elseif  ne == 49        refConfig = Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2")
        elseif  49 <= ne < 54   refConfig = Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^$(ne-48)")
        else    error("stop b")
        end
            
    else        error("stop c")
    end             
    
    return ( refConfig )
end

    
"""
`Basics.extractConfiguration(theme::Basics.LeadingConfiguration, level::ManyElectron.Level)`  
    ... extract the leading configuration of the given level; a conf::Configuration is returned.
"""
function Basics.extractConfiguration(theme::Basics.LeadingConfiguration, level::ManyElectron.Level)
    allConfs = Basics.extractConfigurations(Basics.FromBasis(), level.basis)
    weights  = zeros(length(allConfs))
    for  (ia, allConf) in  enumerate(allConfs)
        for (ic, csf) in enumerate(level.basis.csfs)
            if  allConf == Basics.extractConfiguration(Basics.FromBasis(), level.basis, csf)     
                weights[ia] = weights[ia] + level.mc[ic]^2
            end
        end
    end
    # Determine index of maximum and return the corresponding configuration
    wx   = findmax(weights)
    conf = allConfs[ wx[2] ]
    
    return( conf )
end

    
"""
`Basics.extractConfiguration(theme::Basics.LeadingConfigurationR, level::ManyElectron.Level)`  
    ... extract the leading relativistic configuration of the given level; a conf::ConfigurationR is returned.
"""
function Basics.extractConfiguration(theme::Basics.LeadingConfigurationR, level::ManyElectron.Level)
    allConfs = Basics.extractConfigurations(Basics.RelativisticConfigurations(), level.basis)
    weights  = zeros(length(allConfs))
    for  (ia, allConf) in  enumerate(allConfs)
        for (ic, csf) in enumerate(level.basis.csfs)
            if  allConf == Basics.extractRelativisticConfigurationFromCsfR(csf, level.basis)     
                weights[ia] = weights[ia] + level.mc[ic]^2
            end
        end
    end
    # Determine index of maximum and return the corresponding configuration
    wx   = findmax(weights)
    conf = allConfs[ wx[2] ]
    
    return( conf )
end


"""
`Basics.extractConfiguration(theme::Basics.NonrelativisticBasis, csfNR::LSjj.CsfNR, basisNR::LSjj.BasisNR)`
    ... to extract the nonrelativistic configuration from the given csfNR, if this is part of basisNR.
        A nonrelativistic conf::Configuration is returned.
"""
function Basics.extractConfiguration(theme::Basics.NonrelativisticBasis, csfNR::LSjj.CsfNR, basisNR::LSjj.BasisNR)
    ##x  Replace in code: LSjj.extractConfigurationFromCsfNR(csfNR::LSjj.CsfNR, basisNR::LSjj.BasisNR)
    
    shells = Dict{Shell,Int64}()
    for  s = 1:length(basisNR.shells)
        shell = basisNR.shells[s];    occ = csfNR.occupation[s]
        if  occ > 0     shells = Base.merge( shells, Dict( shell => occ))   end
    end
    conf   = Configuration( shells, basisNR.NoElectrons)
    return( conf )
end


#################################################################################################################################
#################################################################################################################################

    
"""
`Basics.extractConfigurations(theme::Basics.ByNumber, confs::Array{Configuration,1})`  
    ... extract from confs all those configurations with NoElectrons in theme.NoElectrons; a list of 
        newConfs::Array{Configuration,1} is returned. 
"""
function Basics.extractConfigurations(theme::Basics.ByNumber, confs::Array{Configuration,1})
    ##x Replace in code:  Basics.excludeConfigurations
    newConfs = Configuration[]
    for  conf in confs   if conf.NoElectrons in theme.NoElectrons    push!(newConfs, conf)   end   end
    
    return( newConfs )        
end


"""
`Basics.extractConfigurations(theme::Basics.ByParity, confs::Array{Configuration,1})`  
    ... to extract all configurations with the given parity; a list of possible newConfs::Array{Configuration,1} is returned.
"""
function Basics.extractConfigurations(theme::Basics.ByParity, confs::Array{Configuration,1})
    newConfList = Configuration[]
    for conf  in  confs
        parity = Basics.extractFromConfiguration(GetParity(), conf)
        if  parity == theme.P         push!(newConfList, conf)   end 
    end 
    
    return( newConfList )
end


"""
`Basics.extractConfigurations(theme::Basics.ContractShells, confs::Array{Configuration,1})`
    ... to extract/contract all configuration by removing all empty shells; a list of possible newConfs::Array{Configuration,1}
        is returned.
"""
function Basics.extractConfigurations(theme::Basics.ContractShells, confs::Array{Configuration,1})
    newConfs = Configuration[]
    for conf in confs   push!(newConfs, Basics.extractConfiguration(theme::Basics.ContractShells, conf) )   end 

    return( newConfs )
end


"""
`Basics.extractConfigurations(theme::Basics.FromBasis, basis::Basis)`  
    ... extract all (non-relativistic) configurations that contribute to the given set of CSF in basis.csfs. 
        A confList::Array{Configuration,1} is returned.
"""
function Basics.extractConfigurations(theme::Basics.FromBasis, basis::Basis)
    confList  = Configuration[]
    subshells = basis.subshells
    
    for  csf  in basis.csfs
        newShells = Dict{Shell,Int64}();    NoElectrons = 0
        for  i = 1:length(subshells)
            n = subshells[i].n;    l = Basics.subshell_l(subshells[i]);    occ = csf.occupation[i]
            shell = Shell(n,l)
            if    haskey(newShells, shell)    newShells[shell] = newShells[shell] + occ;   NoElectrons = NoElectrons + occ
            else  
                if   occ > 0  newShells = Base.merge( newShells, Dict( shell => occ));     NoElectrons = NoElectrons + occ   end
            end
        end
        if  basis.NoElectrons != NoElectrons    error("stop a")
        else
            conf = Configuration(newShells, NoElectrons)
            if  conf in confList    continue;    else    push!( confList,  conf)    end
        end
    end 
    
    return( confList )
end


"""
`Basics.extractConfigurations(theme::Basics.FromMultiplet, multiplet::Multiplet)`  
    ... extract all (non-relativistic) configurations that contribute to the given multiplet. 
        A confList::Array{Configuration,1} is returned.
"""
function Basics.extractConfigurations(theme::Basics.FromMultiplet, multiplet::Multiplet)  
    confList = Basics.extractConfigurations(Basics.FromBasis(), multiplet.levels[1].basis)
    
    return( confList )
end


"""
`Basics.extractConfigurations(theme::Basics.FromMultiplet, multiplets::Array{Multiplet,1})`  
    ... extract all (non-relativistic) configurations that contribute to the given multiplets. 
        A confList::Array{Configuration,1} is returned.
"""
function Basics.extractConfigurations(theme::Basics.FromMultiplet, multiplets::Array{Multiplet,1}) 
    confList = Configuration[]
    for multiplet in multiplets    append!(confList, Basics.extractConfigurations(Basics.FromBasis(), multiplet.levels[1].basis))   end 
    confList = unique(confList)
    
    return( confList )
end


"""
`Basics.extractConfigurations(theme::Basics.RelativisticConfigurations, basis::Basis)`  
    ... to extract all relativistic configurations from the given basis. 
        A confList::Array{ConfigurationR,1} is returned.
"""
function Basics.extractConfigurations(theme::Basics.RelativisticConfigurations, basis::Basis)
    confList  = ConfigurationR[]
    subshells = basis.subshells
    
    for  csf  in basis.csfs
        newSubshells = Dict{Subshell,Int64}();    NoElectrons = 0
        for  (i, subshell)  in   enumerate(subshells)
            occ = csf.occupation[i]
            if    haskey(newSubshells, subshell)    
                newSubshells[subshell] = newSubshells[subshell] + occ;   NoElectrons = NoElectrons + occ
            else  
                if   occ > 0             newSubshells[subshell] = occ;   NoElectrons = NoElectrons + occ   end
            end
        end
        if  basis.NoElectrons != NoElectrons    error("stop a")
        else
            confR = ConfigurationR(newSubshells, NoElectrons)
            if  confR in confList    continue;    else    push!( confList, confR)    end
        end
    end 
    confList = unique(confList)
    
    return( confList )
end


"""
`Basics.extractConfigurations(theme::Basics.RelativisticConfigurations, basis::Basis, totalJ::AngularJ64)`  
    ... to extract all relativistic configurations from the given basis. Only the CSF with total angular momentum
        totalJ are considered. A confList::Array{ConfigurationR,1} is returned.
"""
function Basics.extractConfigurations(theme::Basics.RelativisticConfigurations, basis::Basis, totalJ::AngularJ64)
    confList  = ConfigurationR[]
    subshells = basis.subshells
    
    for  csf  in basis.csfs
        if  csf.J !=  totalJ    continue   end
        newSubshells = Dict{Subshell,Int64}();    NoElectrons = 0
        for  (i, subshell)  in  enumerate(subshells)
            occ = csf.occupation[i]
            if    haskey(newSubshells, subshell)    
                newSubshells[subshell] = newSubshells[subshell] + occ;   NoElectrons = NoElectrons + occ
            else  
                if   occ > 0             newSubshells[subshell] = occ;   NoElectrons = NoElectrons + occ   end
            end
        end
        if  basis.NoElectrons != NoElectrons    error("stop a")
        else
            confR = ConfigurationR(newSubshells, NoElectrons)
            if  confR in confList    continue;    else    push!( confList, confR)    end
        end
    end 
    confList = unique(confList)
    
    return( confList )
end


"""
`Basics.extractConfigurations(theme::Basics.RestrictExcitations, confs::Array{Configuration,1})`  
    ... extracts a (reduced) list of non-relativistic configurations in which the theme.restrictions are applied to 
        remove unwanted configurations with regard to the given confs. A list of possible newConfs::Array{Configuration,1} 
        is returned.
"""
function Basics.extractConfigurations(theme::Basics.RestrictExcitations, confs::Array{Configuration,1})
    newConfs = Configuration[];    
    for  conf  in  confs
        addConf = true
        for  restriction  in  theme.restrictions
            if  Basics.isViolated(conf, restriction)    addConf = false;    break    end
        end
        
        if  addConf   push!(newConfs, conf)    end 
    end
    
    return( newConfs )
end


"""
`Basics.extractConfigurations(theme::Basics.TotalAM, confs::Array{Configuration,1})`  
    ... to extract all configurations that contribute to CSFR with totalJs;
        a list of possible newConfs::Array{Configuration,1} is returned.
"""
function  Basics.extractConfigurations(theme::Basics.TotalAM, confs::Array{Configuration,1})
    newConfList = Configuration[]   
    for conf  in  confs
        totalJs = Basics.extractFromConfiguration(TotalAM(false, AngularJ64[]), conf)
        if  length( intersect(totalJs, theme.totalJs) ) > 0    push!(newConfList, conf)   end 
    end 
    
    return( newConfList )
end


#################################################################################################################################
#################################################################################################################################

    
"""
`Basics.extractFromConfiguration(theme::Basics.AllShells, conf::Configuration)`  
    ... extract all (occupied) shells from the given conf; a list::Array{Shell,1} is returned.
"""
function Basics.extractFromConfiguration(theme::Basics.AllShells, conf::Configuration)
    shells    = Shell[]
    for (shell, occ)  in  conf.shells   push!(shells, shell)   end
    shells = sort(shells)
    
    return( shells )        
end

    
"""
`Basics.extractFromConfiguration(theme::Basics.ClosedCore, conf::Configuration)`  
    ... extract the (nonrelativistic) closed-cire shells in conf; a closedConf::Configuration is returned.
"""
function Basics.extractFromConfiguration(theme::Basics.ClosedCore, conf::Configuration)
    closedShells = Basics.extractFromConfiguration(Basics.ClosedShells(), conf::Configuration)
    shellList    = Basics.generateShellList(1, 6, 5)
    coreShellD   = Dict{Shell,Int64}();   noElectrons = 0
    #
    for shell in shellList 
        if      !(haskey(conf.shells, shell))           # do nothing, jump over
        elseif  conf.shells[shell] == 2*(2*shell.l + 1)    
            coreShellD[shell] = 2*(2*shell.l + 1);      noElectrons = noElectrons + 2*(2*shell.l + 1) 
        else    break
        end 
    end
    
    return( Configuration(coreShellD, noElectrons) )        
end

    
"""
`Basics.extractFromConfiguration(theme::Basics.ClosedShells, conf::Configuration)`  
    ... extract the (nonrelativistic) closed shells in conf; a list::Array{Shell,1} is returned.
"""
function Basics.extractFromConfiguration(theme::Basics.ClosedShells, conf::Configuration)
    shells = Shell[]
    for  (shell, occ)  in conf.shells
        if  occ == 2*(2*shell.l + 1)    push!(shells, shell)    end
    end
    
    return( shells )        
end

    
"""
`Basics.extractFromConfiguration(theme::Basics.ClosedSubshells, conf::ConfigurationR)`  
    ... extract the (relativistic) closed subshells in conf; a list::Array{Subshell,1} is returned.
"""
function Basics.extractFromConfiguration(theme::Basics.ClosedSubshells, conf::ConfigurationR)
    subshells = Subshell[]
    for  (subshell, occ)  in conf.subshells
        if  occ == Basics.twice(Basics.subshell_j(subshell)) + 1    push!(subshells, subshell)    end
    end
    
    return( subshells )        
end



"""
`Basics.extractFromConfiguration(theme::Basics.ClosedSubshells, "[Ne]")`  
    ... to extract the list of (relativistic) subshells for the given closed-shell configuration. Allowed closed-shell strings 
        are [He], [Ne], [Mg], [Ar], [Kr] and [Xe].
"""
function Basics.extractFromConfiguration(theme::Basics.ClosedSubshells, sa::String)
    if       sa == "[He]"    wa = [ Subshell("1s_1/2")]    
    elseif   sa == "[Ne]"    wa = [ Subshell("1s_1/2"), Subshell("2s_1/2"), Subshell("2p_1/2"), Subshell("2p_3/2")]   
    elseif   sa == "[Mg]"    wa = [ Subshell("1s_1/2"), Subshell("2s_1/2"), Subshell("2p_1/2"), Subshell("2p_3/2"), Subshell("3s_1/2")]   
    elseif   sa == "[Ar]"    wa = [ Subshell("1s_1/2"), Subshell("2s_1/2"), Subshell("2p_1/2"), Subshell("2p_3/2"), 
                                                        Subshell("3s_1/2"), Subshell("3p_1/2"), Subshell("3p_3/2")]   
    elseif   sa == "[Kr]"    wa = [ Subshell("1s_1/2"), Subshell("2s_1/2"), Subshell("2p_1/2"), Subshell("2p_3/2"), 
                                                        Subshell("3s_1/2"), Subshell("3p_1/2"), Subshell("3p_3/2"),   
                                                        Subshell("3d_3/2"), Subshell("3d_5/2"), 
                                    Subshell("4s_1/2"), Subshell("4p_1/2"), Subshell("4p_3/2")]   
    elseif   sa == "[Xe]"    wa = [ Subshell("1s_1/2"), Subshell("2s_1/2"), Subshell("2p_1/2"), Subshell("2p_3/2"), 
                                                        Subshell("3s_1/2"), Subshell("3p_1/2"), Subshell("3p_3/2"),   
                                                        Subshell("3d_3/2"), Subshell("3d_5/2"), 
                                    Subshell("4s_1/2"), Subshell("4p_1/2"), Subshell("4p_3/2"), 
                                                        Subshell("4d_3/2"), Subshell("4d_5/2"), 
                                    Subshell("5s_1/2"), Subshell("5p_1/2"), Subshell("5p_3/2")]   
    else
        error("Unsupported keystring = $sa ")
    end
    
    return( wa )
end

    
"""
`Basics.extractFromConfiguration(theme::Basics.ExcitationLevel, conf::Configuration)`  
    ... determines the excitation level of a given configuration, which can be used to compare the degree of excitation between 
        two configurations:  exLevelA <= exLevelB.  An exLevel::Int64 is returned. This degree of excitation works nicely for multiply 
        and highly-charged ions with the same number of electrons as well as for the comparison of inner-shell holes but 
        (may) fail for pure valence-shell excitations, whose relative energies may overlap and are difficult to predict without 
        explicit computations. This procedure can be further advanced, if the need arises.
"""
function Basics.extractFromConfiguration(theme::Basics.ExcitationLevel, conf::Configuration)
    exLevel = 0
    for  (shell, occ)  in  conf.shells
        exLevel = exLevel + occ * (-1000 + 10*shell.n + shell.l)
    end
    
    return( exLevel )        
end


"""
`Basics.extractFromConfiguration(theme::Basics.GetParity, conf::Configuration)`  
    ... to determine the parity::Parity (Enum Parity) of a given non-relativistic configuration.
"""
function Basics.extractFromConfiguration(theme::Basics.GetParity , conf::Configuration)
    par = 1
    for  (k,v) in conf.shells
    if   iseven(k.l)    p = 1   else   p = -1    end    
    par = par * (p^v)
    end

    if       par == 1    return( Basics.plus )  
    elseif   par == -1   return( Basics.minus )
    else     error("stop a")
    end  
end


"""
`Basics.extractFromConfiguration(theme::Basics.GetParity, conf::ConfigurationR)`  
    ... to determine the parity::Parity (Enum Parity) of a given relativistic configuration.
"""
function Basics.extractFromConfiguration(theme::Basics.GetParity, conf::ConfigurationR)
    par = 1
    for  (k,v) in conf.subshells
    l = Basics.subshell_l(k)
    if   iseven(l)    p = 1   else   p = -1    end    
    par = par * (p^v)
    end

    if       par == 1    return( Basics.plus )  
    elseif   par == -1   return( Basics.minus )
    else     error("stop a")
    end  
end

    
"""
`Basics.extractFromConfiguration(theme::Basics.IsOccupied, conf::Configuration, shell::Shell)`  
    ... determines of whether the shell is occupied (> 0) in the given conf.
"""
function Basics.extractFromConfiguration(theme::Basics.IsOccupied, conf::Configuration, shell::Shell)
    if  !(haskey(conf.shells, shell))     isOccupied = false
    elseif  0 < conf.shells[shell]        isOccupied = true 
    else                                  isOccupied = false 
    end
    
    return( isOccupied )        
end

    
"""
`Basics.extractFromConfiguration(theme::Basics.Multiplicity, conf::Configuration)`  
    ... extracts the multiplicity, i.e. the g-factor, of the configuration conf if all of its levels are considered to be 
        degenerate; this multiplicity is frequently used for empirical computations that are based on configuration-averaged 
        plasma levels. A g::Int64 is returned.
"""
function Basics.extractFromConfiguration(theme::Basics.Multiplicity, conf::Configuration)
    g = 0
    totalAMs = Basics.extractFromConfiguration(Basics.TotalAM(true, AngularJ64[]), conf::Configuration)
    for  totalAM   in  totalAMs    g = g + Basics.twice(totalAM) + 1   end
    
    return( g )
end

    
"""
`Basics.extractFromConfiguration(theme::Basics.OpenShellNumber, conf::Configuration)`  
    ... determine the number of open (nonrelativistic) shells in conf; a singleton of type LSjj.AbstractOpenShell 
        is returned. This method still need to be combined with the module LSjj.
"""
function Basics.extractFromConfiguration(theme::Basics.OpenShellNumber, conf::Configuration)
    ##x Replace in code: Basics.extractNoOpenShells(conf::Configuration)
    ns = 0
    for  (shell, occ)  in conf.shells
        if      occ == 0  ||  occ == 2*(2*shell.l + 1)
        else    ns = ns + 1
        end
    end
    
    if      ns == 0   return ( LSjj.ZeroOpenShell() )
    elseif  ns == 1   return ( LSjj.OneOpenShell() )
    elseif  ns == 2   return ( LSjj.TwoOpenShells() )
    elseif  ns == 3   return ( LSjj.ThreeOpenShells() )
    else    error("stop a")
    end
    
end

    
"""
`Basics.extractFromConfiguration(theme::Basics.OpenShells, conf::Configuration)`  
    ... extract the open (nonrelativistic) shells in conf; a list::Array{Shell,1} is returned.
"""
function Basics.extractFromConfiguration(theme::Basics.OpenShells, conf::Configuration)
    shells = Shell[]
    for  (shell, occ)  in conf.shells
        if      occ == 0  ||  occ == 2*(2*shell.l + 1)
        else    push!(shells, shell)
        end
    end
    
    return( shells )        
end

    
"""
`Basics.extractFromConfiguration(theme::Basics.OpenSubshells, conf::Configuration)`  
    ... extract the open (relativistic) subshells in conf; a list::Array{Subshell,1} is returned.
"""
function Basics.extractFromConfiguration(theme::Basics.OpenSubshells, conf::Configuration)
    shells    = Basics.extractFromConfiguration(Basics.OpenShells(), conf::Configuration)
    subshells = Subshell[]
    
    for shell  in  shells
        subshellOccupations = Basics.shellSplitOccupation(shell, conf.shells[shell]);   @show subshellOccupations
        #
        for subshellOcc  in  subshellOccupations
            for (subsh, occ)  in  subshellOcc
                if      occ == 0  ||  occ == 2*(Basics.twice( Basics.subshell_j(subsh) ) + 1)
                else    push!(subshells, subsh)
                end
            end
        end
    end
    subshells = unique(subshells)
    
    return( subshells )        
end

    
"""
`Basics.extractFromConfiguration(theme::Basics.OpenSubshells, conf::ConfigurationR)`  
    ... extract the open (relativistic) subshells in conf; a list::Array{Subshell,1} is returned.
"""
function Basics.extractFromConfiguration(theme::Basics.OpenSubshells, conf::ConfigurationR)
    subshells = Subshell[]
    for  (subsh, occ)  in conf.subshells
        if      0 < occ < Basics.twice(Basics.subshell_j(subsh)) + 1    push!(subshells, subsh)    end
    end
    
    return( subshells )        
end

    
"""
`Basics.extractFromConfiguration(theme::Basics.TotalAM, conf::Configuration)`  
    ... extract the total angular momenta that are associated with the fine-structure of the given configuration; 
        a totalJ::Array{AngularJ64,1} is returned.
"""
function Basics.extractFromConfiguration(theme::Basics.TotalAM, conf::Configuration)
    function  extractJsFromSubshellOccupations(subshellOcc::Dict{Basics.Subshell,Int64})
        currentJs  = [ AngularJ64(0) ]
        for (subsh, occ) in  subshellOcc
            if  occ == 0    continue    end
            stateList = ManyElectron.provideSubshellStates(subsh, occ);   newJs = AngularJ64[]
            for  J  in  currentJs, state  in  stateList    
                # Couple angular momenta from currentJs with new state-Js and add the result to newJs
                append!(newJs, Basics.oplus( J, AngularJ64(state.Jsub2//2)) )
            end
            currentJs  = deepcopy(newJs)
        end
        return( currentJs )
    end
    
    # Couple all angular momenta due to a cycle through all open shells, associated subshellOccupations
    currentJs  = [ AngularJ64(0) ];  
    openShells = Basics.extractFromConfiguration(Basics.OpenShells(), conf::Configuration)
    
    for shell  in  openShells
        subshellOccupations = Basics.shellSplitOccupation(shell, conf.shells[shell]);   newJs = AngularJ64[]
        for  J  in  currentJs,  subshellOcc  in  subshellOccupations
            subshellOccJs = extractJsFromSubshellOccupations(subshellOcc)
            for subshJ in subshellOccJs
                # Couple angular momenta from currentJs with subshellOccJs and add the result to newJs
                append!(newJs, Basics.oplus( J, subshJ) )
            end
        end
        currentJs  = deepcopy(newJs)
    end
    
    if     theme.allJ    finalJs = sort(currentJs)
    else                finalJs = sort(  unique(currentJs) ) 
    end
    
    return( finalJs )        
end

    
"""
`Basics.extractFromConfiguration(theme::Basics.ValenceOccupation, conf::Configuration, coreConf::Configuration)`  
    ... extracts the valence-part of the given configuration conf without those shells that are given in the core 
        configuration coreConf already. The procedure terminates if the occupation of coreConf does not coincide with 
        those of conf. A valenceConf::Configuration is returned.
"""
function Basics.extractFromConfiguration(theme::Basics.ValenceOccupation, conf::Configuration, coreConf::Configuration)
    coreShells = coreConf.shells;    valenceShells = Dict{Shell,Int64}()
    for  (shell, occ)  in  conf.shells
        if       haskey(coreShells, shell)  &&   coreShells[shell] == occ
        elseif   haskey(coreShells, shell)       error("stop a") 
        else     valenceShells[shell] = occ
        end 
    end
    
    return( Configuration(valenceShells, conf.NoElectrons - coreConf.NoElectrons) )        
end

    
"""
`Basics.extractFromConfiguration(theme::Basics.ValenceShells, conf::Configuration)`  
    ... determines the valence shells of a configuration, i.e. all shells including and beyond those with the 
        most inner-shell hole in conf. Hence, the first shell always indicates the inner-shell hole. 
        A valenceShells::Array{Shell,1} is returned
"""
function Basics.extractFromConfiguration(theme::Basics.ValenceShells, conf::Configuration)
    valenceShells = Shell[];   firstShellFound = false
    shellList     = Basics.generateShellList(1, 8, 7)
    for  shell  in  shellList
        if      !firstShellFound  &&  haskey(conf.shells, shell)   &&   conf.shells[shell] < 2*(2*shell.l + 1)
                firstShellFound = true;                            push!(valenceShells, shell)
        elseif  firstShellFound   &&  haskey(conf.shells, shell)   push!(valenceShells, shell) 
        end 
    end
    
    return( valenceShells )        
end


#################################################################################################################################
#################################################################################################################################

    
"""
`Basics.extractFromConfigurations(theme::Basics.AllShells, confs::Array{Configuration,1})`  
    ... extract all (occupied) shells from the given configurations; a list::Array{Shell,1} is returned.
"""
function Basics.extractFromConfigurations(theme::Basics.AllShells, confs::Array{Configuration,1})
    shells    = Shell[]
    for conf in confs
       for (shell, occ)  in  conf.shells   if !(shell in shells)    push!(shells, shell)   end  end 
    end 
    shells = unique(shells)
    shells = sort(shells)
    
    return( shells )        
end

    
"""
`Basics.extractFromConfigurations(theme::Basics.ClosedCore, confs::Array{Configuration,1})`  
    ... extract the (nonrelativistic) closed-core configuration that is common to all the given configurations confs; 
        a closedConf::Configuration is returned.
"""
function Basics.extractFromConfigurations(theme::Basics.ClosedCore, confs::Array{Configuration,1})
    coreShells = Dict{Shell,Int64}();   newCoreShells = Dict{Shell,Int64}();   newNoElectrons = 0
    if  length(confs) == 0   
    else
        for  (i, conf) in enumerate(confs)
            if  i == 1   coreShells = deepcopy(conf.shells)
            else
                # Reduce the shells in coreShellD if the occupation does not agree
                for  (shell, occ) in coreShells
                    if    haskey(conf.shells, shell)  &&  conf.shells[shell] == coreShells[shell]  # do nothing
                    else  coreShells[shell] = 0
                    end
                end
            end
        end
    end
    
    # Now determine which shells are still equally filled
    for  (shell, occ) in coreShells
        if  occ > 0    newCoreShells[shell] = occ;    newNoElectrons = newNoElectrons + occ  end
    end
    
    return( Configuration(newCoreShells, newNoElectrons) )        
end

    
"""
`Basics.extractFromConfigurations(theme::Basics.ClosedShells, confs::Array{Configuration,1})`  
    ... extract the (nonrelativistic) shells that are closed in all configurations confs; a list::Array{Shell,1} is returned.
"""
function Basics.extractFromConfigurations(theme::Basics.ClosedShells, confs::Array{Configuration,1})
    closedShells = Shell[];    shells = Shell[]
    
    if    length(confs) == 0      # trivial case
    else  # Select all shells from first configuration and test which of them are closed in all confs
        for  (shell, occ)  in confs[1].shells    push!(shells, shell)    end  
        for  shell  in  shells
            for  conf  in  confs
                if      !(haskey(conf.shells, shell))                  break
                elseif  conf.shells[shell] != 2*(2*shell.l + 1)        break
                else    # do nothing
                end 
            end
            push!(closedShells, shell) 
        end
    end 
    
    return( closedShells )        
end

    
"""
`Basics.extractFromConfigurations(theme::Basics.MeanOccupation, confs::Array{Configuration,1})`  
    ... extract the (mean) occupation numbers from the given configurations. A dictionary meanShells::Dict{Shell,Float64} 
        is returned.
"""
function Basics.extractFromConfigurations(theme::Basics.MeanOccupation, confs::Array{Configuration,1})
    shellList  = Basics.generateShellList(1, 10, "k");   sa = "\n   +  "
    meanShells = Dict{Shell,Float64}();    
    
    if      length(confs) == 0
        println(stream, "\n   + No configurations are given.") 
    else
        for conf in confs
            for  (shell, occ)  in  conf.shells
                if  haskey(meanShells, shell)    meanShells[shell] = meanShells[shell] + occ
                else                             meanShells[shell] = occ
                end 
            end
        end
        #  Divide occupations be number of configurations
        for   (shell, occ)  in  meanShells    meanShells[shell] = meanShells[shell] / length(confs)   end
    end 
    
    return( meanShells )
end

    
"""
`Basics.extractFromConfigurations(theme::Basics.NumberOfElectrons, confs::Array{Configuration,1})`  
    ... extracts the number of electrons of all given configurations in a list::Array{Int64,1}.
        This list can readily be checked for different requests, for instance, that all configurations possess
        the same number of electrons. A intList::Array{Int64,1} is returned.
"""
function Basics.extractFromConfigurations(theme::Basics.NumberOfElectrons, confs::Array{Configuration,1})
    NoElectrons = Int64[]
    for  conf in confs    push!(NoElectrons, conf.NoElectrons)   end
    
    return( NoElectrons )        
end


"""
`Basics.extractFromConfigurations(theme::Basics.OccupationDifference, confa::Configuration, confb::Configuration)`  
    ... extract the differences in the occupation of shells: occupation(confa) - occupation(confb).
        Only those shells are returned for which qa - qb != 0.  A list::Array{Tuple(Shell, Int64),1} is returned. 
"""
function Basics.extractFromConfigurations(theme::Basics.OccupationDifference, confa::Configuration, confb::Configuration)
    occList   = Tuple{Shell, Int64}[]
    shellList = Basics.extractShellList([confa, confb])
    
    for shell in shellList
        if haskey(confa.shells, shell)  qa = confa.shells[shell]   else   qa = 0  end
        if haskey(confb.shells, shell)  qb = confb.shells[shell]   else   qb = 0  end
        if qa - qb != 0     push!(occList, (shell, qa - qb) )                     end
    end
    
    return( occList )
end


"""
`Basics.extractFromConfigurations(theme::Basics.OccupationDifference, confa::ConfigurationR, confb::ConfigurationR)`  
    ... extract the differences in the occupation of subshells: occupation(confa) - occupation(confb).
        Only those subshells are returned for which qa - qb != 0.  A list::Array{Tuple(Subshell, Int64),1} is returned. 
"""
function Basics.extractFromConfigurations(theme::Basics.OccupationDifference, confa::ConfigurationR, confb::ConfigurationR)
    occList   = Tuple{Subshell, Int64}[]
    shList    = Basics.extractSubshellList([confa, confb])
    
    for subsh in shList
        if haskey(confa.subshells, subsh)  qa = confa.subshells[subsh]   else   qa = 0  end
        if haskey(confb.subshells, subsh)  qb = confb.subshells[subsh]   else   qb = 0  end
        if qa - qb != 0     push!(occList, (subsh, qa - qb) )                           end
    end
    
    return( occList )
end


#################################################################################################################################
#################################################################################################################################


"""
`Basics.generateConfigurations(theme::Basics.AddElectrons, confs::Array{Configuration,1})`  
    ... generates a list of non-relativistic configurations in which theme.ne electrons are added into the theme.intoshells.
        The routine updates stepwise the list of configurations until all the required electrons are added, and by
        taking the "filling themes" of the shells into account. A list of possible newConfs::Array{Configuration,1} 
        is returned.
"""
function Basics.generateConfigurations(theme::Basics.AddElectrons, confs::Array{Configuration,1})
    currentConfs = confs;    newConfs = Configuration[]
    #
    for na = 1:theme.ne
        newConfs = Configuration[]
        for  shell in theme.intoShells
            for  conf in currentConfs
                newShells = deepcopy(conf.shells);   NoElectrons = conf.NoElectrons
                if  haskey(newShells, shell)  
                    if  newShells[shell] + 1 > 4*shell.l + 2     continue     end    
                        newShells[shell] = newShells[shell] + 1
                else    newShells[shell] = 1
                end
                push!(newConfs, Configuration(newShells, NoElectrons+1))
            end
        end 
        newConfs     = unique(newConfs);    currentConfs = deepcopy(newConfs)
    end        
    
    return( newConfs )
end


"""
`Basics.generateConfigurations(theme::Basics.ExciteElectrons, confs::Array{Configuration,1})`  
    ... generates a list of non-relativistic configurations in which theme.ne electrons are excited from the theme.fromShells
        into the theme.intoShells. The routine updates stepwise the list of configurations until all the required electrons 
        are excited, and by taking the (de-)excitation of the shells into account. A list of possible 
        newConfs::Array{Configuration,1} is returned.
"""
function Basics.generateConfigurations(theme::Basics.ExciteElectrons, confs::Array{Configuration,1})
    currentConfs = confs;   newConfs = Configuration[] 
    #
    for na = 1:theme.ne
        newConfs = Configuration[]
        for  fromShell in theme.fromShells,   intoShell in theme.intoShells
            for  conf in currentConfs
                newShells = deepcopy(conf.shells);   NoElectrons = conf.NoElectrons
                if  fromShell == intoShell                               continue     end
                if  haskey(newShells, fromShell)  
                    if  newShells[fromShell] == 0                        continue     end    
                        newShells[fromShell] = newShells[fromShell] - 1
                else                                                     continue
                end
                if  haskey(newShells, intoShell)  
                    if  newShells[intoShell] + 1 > 4*intoShell.l + 2     continue     end    
                        newShells[intoShell] = newShells[intoShell] + 1
                else    newShells[intoShell] = 1
                end
                push!(newConfs, Configuration(newShells, NoElectrons))
            end
        end 
        newConfs = unique(newConfs);    currentConfs = deepcopy(newConfs)
    end        
    
    return( newConfs )
end


"""
`Basics.generateConfigurations(theme::Basics.RemoveElectrons, confs::Array{Configuration,1})`  
    ... generates a list of non-relativistic configurations in which theme.ne electrons are removed from the theme.fromShells.
        The routine updates stepwise the list of configurations until all the required electrons are removed, and by
        taking the removal of electrons from  the shells into account. A list of possible newConfs::Array{Configuration,1} 
        is returned.
"""
function Basics.generateConfigurations(theme::Basics.RemoveElectrons, confs::Array{Configuration,1})
    currentConfs = confs;   newConfs = Configuration[]    
    #
    for na = 1:theme.ne
        newConfs = Configuration[]
        for  shell in theme.fromShells
            for  conf in currentConfs
                newShells = deepcopy(conf.shells);   NoElectrons = conf.NoElectrons
                if  haskey(newShells, shell)  
                    if  newShells[shell] == 0     continue     end    
                        newShells[shell] = newShells[shell] - 1
                else    continue
                end
                push!(newConfs, Configuration(newShells, NoElectrons-1))
            end
        end 
        newConfs = unique(newConfs);    currentConfs = deepcopy(newConfs)
    end        
    
    return( newConfs )
end


"""
`Basics.generateConfigurations(theme::Basics.RestrictExcitations, refConfigs::Array{Configuration,1})`  
    ... generates a (reduced) list of non-relativistic configurations in which theme.ne electrons are excited from the 
        theme.fromShells into the theme.toShells and where, in addition, the set of theme.restrictions is taken into 
        account. This is a simplified procedure to generate RAS sets of configurations. A list of possible 
        confList::Array{Configuration,1} is returned.
"""
function Basics.generateConfigurations(theme::Basics.RestrictExcitations, refConfigs::Array{Configuration,1})
    ## Replace in code: Basics.generateConfigurations(refConfigs::Array{Configuration,1}, fromShells::Array{Shell,1}, toShells::Array{Shell,1}, 
    ##                         noex::Int64; restrictions::Array{AbstractConfigurationRestriction,1}=AbstractConfigurationRestriction[])
    newConfigs = Configuration[];    
    #
    if     theme.ne == 0
        newConfigs = refConfigs
    elseif theme.ne > 0
        exciteE    = Basics.ExciteElectrons(theme.ne, theme.fromShells, theme.toShells)
        newConfigs = Basics.generateConfigurations(exciteE, refConfigs)
        newConfigs = Base.unique(newConfigs)
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
`Basics.generateConfigurations(theme::Basics.ForAutoIonization, confs::Array{Configuration,1})`  
    ... generates a list of non-relativistic configurations in which an inner-shell hole is filled in each given configuration 
        and two new holes appear in some outer shells. All combinations of two outer shells are taken into account, while the most 
        inner-shell hole is always filled. The procedure terminates if 
            (i)   no inner-shell hole can be identified, 
            (ii)  configurations with a different number of electrons are given or 
            (iii) if the applied algorithm does result in configurations not appropriate for a single Auger emission. 
        The routine excludes the generation of configurations with intra-shell autoionization paths, for which the free-energy 
        of the emitted electron arises purely from the re-coupling of the valence-shell electrons. A list of possible 
        newConfs::Array{Configuration,1} is returned.
"""
function Basics.generateConfigurations(theme::Basics.ForAutoIonization, confs::Array{Configuration,1})
    newConfs = Configuration[]
    # Test for equal number of electrons
    wa = Basics.extractFromConfigurations(NumberOfElectrons(), confs)
    if  length(wa) != 1   error("Configurations with different NoElectrons are given.")   else   noe = wa[1]    end 
    
    for  conf in confs
        valenceShells  = Basics.extractFromConfiguration(Basics.ValenceShells(),   conf)
        # Determine configuration due to an de-excitation of an electron within valenceShells
        if     length(valenceShells) < 2   @warn("Less than two valence shells are found for autoionization: $valenceShells ")
               ## newConfigs = Configuration[]
               return( Configuration[] )
        else 
               addConfigs = Basics.generateConfigurations(AddElectrons(1, valenceShells[1:1]), [conf])
               newConfigs = Basics.generateConfigurations(RemoveElectrons(2, valenceShells[2:end]), addConfigs)
        end
        append!(newConfs, newConfigs)
    end 
    newConfs = unique(newConfs)
    
    # Check that all generated configurations have one electron less than given initially
    wa = unique( Basics.extractFromConfigurations(NumberOfElectrons(), newConfs) )
    if  length(wa) != 1   ||   wa[1] != noe - 1        
        Basics.displayConfigurations(stdout, newConfs)
        error("Configurations with unexpected NoElectrons are generated.") 
    end 
    
    return( newConfs )
end


"""
`Basics.generateConfigurations(theme::Basics.ForDielectronicCapture, confs::Array{Configuration,1})`  
    ... generates a list of non-relativistic configurations in which an electron from the theme.fromShells is excited 
        to the theme.toShells and an additional electron is capture into the theme.intoShells. These excitation-with-capture 
        refers to the first steps of the dielectronic recombination process and is realized for each given configuration.
        The procedure terminates if 
            (i)  configurations with a different number of electrons are given   or 
            (ii) if the applied algorithm does not results in configurations just one electron `more' than given originally.    
        A list of possible newConfs::Array{Configuration,1} is returned.
"""
function Basics.generateConfigurations(theme::Basics.ForDielectronicCapture, confs::Array{Configuration,1})
    # Check that all configurations have the same number of electrons
    wa = Basics.extractFromConfigurations(NumberOfElectrons(), confs)
    if  length(wa) != 1   error("Configurations with different NoElectrons are given.")   else   noe = wa[1]    end 
    
    newConfs = Basics.generateConfigurations(ExciteElectrons(1, theme.fromShells, theme.toShells), confs)
    newConfs = Basics.generateConfigurations(AddElectrons(1, theme.intoShells), newConfs)
    
    # Check that all generated configurations have one electron less than given initially
    wa = Basics.extractFromConfigurations(NumberOfElectrons(), newConfs)
    if       length(wa) != 1     error("Configurations with different NoElectrons are generated.")   
    elseif   wa[1] != noe + 1    error("Configurations with unexpected NoElectrons are generated.") 
             Basics.displayConfigurations(stdout, newConfs)
    end 
    
    return( newConfs )
end


"""
`Basics.generateConfigurations(theme::Basics.ForDielectronicRecombination, confs::Array{Configuration,1})`  
    ... generates a list of (non-relativistic) intermediate and final configurations, based on theme.fromShells, theme.toShells, 
        theme.intoShells as well as theme.decayShells. The intermediate configurations are formed by exciting one electron 
        fromShells --> toShells and, in addition, by adding one electron into the into shells. All generated configurations must 
        have N+1 electrons of course. For the final configurations, we start from these intermediate configurations and allow 
        de-excitations from toShells   --> fromShells (+ the low-lying decayShells)   as well as
        de-excitations from intoShells --> decayShells. The procedure terminates if 
            (i)  configurations with a different number of electrons are given   or 
            (ii) if the applied algorithm does not results in configurations with a proper number of electrons `relative' to the 
                 given ones. 
        A tuple of two lists (intermediateConfs::Array{Configuration,1}, finalConfs::Array{Configuration,1}) is returned.
"""
function Basics.generateConfigurations(theme::Basics.ForDielectronicRecombination, confs::Array{Configuration,1})
    intermediateConfs = Configuration[];    finalConfs = Configuration[]
    wa = Basics.extractFromConfigurations(NumberOfElectrons(), confs);     wa = unique(wa)
    if  length(wa) != 1   error("Configurations with different NoElectrons are given.")   else   noe = wa[1]    end 
    
    # Generate the intermediate configurations
    newConfs = Basics.generateConfigurations(ExciteElectrons(1, theme.fromShells, theme.toShells), confs)
    newConfs = Basics.generateConfigurations(AddElectrons(1, theme.intoShells), newConfs)
    newConfs = Basics.extractConfigurations(ContractShells(), newConfs)
    
    # Check that all generated configurations have one electron more than given initially 
    wa = Basics.extractFromConfigurations(NumberOfElectrons(), newConfs);     wa = unique(wa)
    if       length(wa) != 1     error("Intermediate configurations with different NoElectrons are generated.")   
    elseif   wa[1] != noe + 1    error("Intermediate configurations with unexpected NoElectrons are generated.") 
             Basics.displayConfigurations(stdout, newConfs)
    end
    intermediateConfs = deepcopy(newConfs)
    
    # Generate the final configurations
    toShells = deepcopy(theme.fromShells)
    for shell in theme.decayShells    if shell < theme.toShells[end]    push!(toShells, shell)   end    end
    aConfs   = Basics.generateConfigurations(ExciteElectrons(1, theme.toShells, toShells), intermediateConfs); 
    aConfs   = Basics.extractConfigurations(ContractShells(), aConfs)
    bConfs   = Basics.generateConfigurations(ExciteElectrons(1, theme.intoShells, theme.decayShells), intermediateConfs);      
    bConfs   = Basics.extractConfigurations(ContractShells(), bConfs)
    finalConfs = unique( append!(aConfs, bConfs) );    
    
    # Check that all generated configurations have one electron more than given initially 
    wa = Basics.extractFromConfigurations(NumberOfElectrons(), finalConfs);     wa = unique(wa)
    if       length(wa) != 1     error("Final configurations with different NoElectrons are generated.")   
    elseif   wa[1] != noe + 1    error("Final configurations with unexpected NoElectrons are generated.") 
             Basics.displayConfigurations(stdout, finalConfs)
    end
    
    # Intermediate and final configurations need to be defined with the (sequence of) subshells; check and expand configurations
    nShells = sort( Basics.extractFromConfigurations(Basics.AllShells(), intermediateConfs) )
    fShells = sort( Basics.extractFromConfigurations(Basics.AllShells(), finalConfs) )
    if  nShells != fShells
        println(">> Expand shells for (intermediate) configurations")
        allShells = sort( unique(append!(nShells, fShells)) )
        nConf  = Basics.extractConfiguration(Basics.ExpandShells(allShells), intermediateConfs[1])
        nConfs = [nConf];   append!(nConfs, intermediateConfs[2:end])
        fConf  = Basics.extractConfiguration(Basics.ExpandShells(allShells), finalConfs[1])
        fConfs = [fConf];   append!(fConfs, finalConfs[2:end])
        (intermediateConfs, finalConfs) = (nConfs, fConfs)
    end
    
    return( (intermediateConfs, finalConfs) )
end


"""
`Basics.generateConfigurations(theme::Basics.ForElectronCapture, confs::Array{Configuration,1})`  
    ... generates a list of non-relativistic configurations in which an additional electron is captured into one of the 
        theme.intoShells, so that each configuration has now one electron more. This is done for all given configurations.
        This is formally the same as the theme::Basics.ForRadiativeRecombination but may later allow to fullfill some 
        additional requirements for non-radiative electron capturForElectronCapturee processes. A list of possible 
        newConfs::Array{Configuration,1} is returned.
"""
function Basics.generateConfigurations(theme::Basics.ForElectronCapture, confs::Array{Configuration,1})
    # Add one electron in each of the theme.intoShells    
    newConfs = Basics.generateConfigurations(AddElectrons(1, theme.intoShells), confs)
    
    return( newConfs )
end


"""
`Basics.generateConfigurations(theme::Basics.ForHollowIons, confs::Array{Configuration,1})`  
    ... generates a list of non-relativistic configurations in which one or several additional electron are captured into 
        the theme.intoShells, so that each configuration has a theme.ne electrons more. This is done for all given configurations.
        A list of possible newConfs::Array{Configuration,1} is returned.
"""
function Basics.generateConfigurations(theme::Basics.ForHollowIons, confs::Array{Configuration,1})
    # Test for equal number of electrons
    wa = Basics.extractFromConfigurations(NumberOfElectrons(), confs)
    if  length(wa) != 1   error("Configurations with different NoElectrons are given.")   else   noe = wa[1]    end
    
    # First capture theme.ne electrons into theme.intoShells
    newConfs = Basics.generateConfigurations(AddElectrons(theme.ne, theme.intoShells), confs)
    #
    # Build configurations with all decayShells 'in between', even if zero occupation
    currentConfs = Configuration[]
    for  conf  in  newConfs 
        nshells = copy(conf.shells)
        for shell in theme.decayShells    if  haskey(nshells, shell)   else    nshells[shell] = 0    end   end
        push!(currentConfs, Configuration(nshells, conf.NoElectrons))
    end
    
    # Now add all decay configurations
    while  true
        newConfigs = Configuration[]
        for  conf in currentConfs
            valenceShells  = Basics.extractFromConfiguration(Basics.ValenceShells(),   conf)
            # Determine configuration due to an de-excitation of an electron within valenceShells
            if     length(valenceShells) < 2
            else 
                addConfigs = Basics.generateConfigurations(AddElectrons(1, valenceShells[1:1]), [conf])
                newConfsPE = Basics.generateConfigurations(RemoveElectrons(1, valenceShells[2:end]), addConfigs)
                newConfsAI = Basics.generateConfigurations(RemoveElectrons(2, valenceShells[2:end]), addConfigs)
                append!(newConfigs, newConfsPE, newConfsAI)
            end
        end
        newConfigs   = unique(newConfigs)
        currentConfs = deepcopy(newConfigs)
        append!(newConfs, newConfigs)
        if  length(newConfigs) == 0   break   end
    end 
    #
    newConfs = unique(newConfs)
    wb = Basics.extractFromConfigurations(NumberOfElectrons(), newConfs)
    
    return( newConfs )
end


"""
`Basics.generateConfigurations(theme::Basics.ForPhotoIonization, confs::Array{Configuration,1})`  
    ... generates a list of non-relativistic configurations in which an additional inner-shell hole is created in all occupied 
        shells and for each given configuration. This creation of an additional hole is independent of whether the photon energy is 
        sufficient to release such an photo-electron. A list of possible newConfs::Array{Configuration,1} is returned.
"""
function Basics.generateConfigurations(theme::Basics.ForPhotoIonization, confs::Array{Configuration,1})
    occShells = Shell[]
    for  conf in confs
        for (shell, occ) in conf.shells 
            if shell in occShells    else    push!(occShells, shell)    end 
        end 
    end

    # Remove one electron from all occupied shells 
    newConfs = Basics.generateConfigurations(RemoveElectrons(1, occShells), confs)
    
    return( newConfs )
end


"""
`Basics.generateConfigurations(theme::ForPhotoEmission, confs::Array{Configuration,1})`  
    ... generates a list of non-relativistic configurations in which an inner-shell hole is filled due to the creation/move to 
        an outer-shell hole. All given configurations are excluded from the generated list, even if they can be formed by 
        filling a lower-lying shell. The routine terminates if 
            (i)   no inner-shell hole can be identified, 
            (ii)  configurations with a different number of electrons are given   or 
            (iii) if the applied algorithm does result in configurations with a different number of electrons. 
        The routine also excludes the generation of configurations due to intra-shell photoemission, for which the emitted photon 
        energy purely arises from the re-coupling of valence-shell electrons. A list of possible newConfs::Array{Configuration,1} 
        is returned.
"""
function Basics.generateConfigurations(theme::Basics.ForPhotoEmission, confs::Array{Configuration,1})
    newConfs = Configuration[]
    # Test for equal number of electrons
    wa = unique( Basics.extractFromConfigurations(NumberOfElectrons(), confs) )
    if  length(wa) != 1   error("Configurations with different NoElectrons are given.")   else   noe = wa[1]    end 
    
    for  conf in confs
        valenceShells  = Basics.extractFromConfiguration(Basics.ValenceShells(),   conf)
        # Determine configuration due to an de-excitation of an electron within valenceShells
        if     length(valenceShells) < 2   @warn("Less than two valence shells are found for photoemission: $valenceShells ")
               newConfigs = Configuration[]
               ##x return( Configuration[] )
        else 
               addConfigs = Basics.generateConfigurations(AddElectrons(1, valenceShells[1:1]), [conf])
               newConfigs = Basics.generateConfigurations(RemoveElectrons(1, valenceShells[2:end]), addConfigs)
        end
        append!(newConfs, newConfigs)
    end 
    newConfs = unique(newConfs)
    
    # Check that all generated configurations have one electron less than given initially
    wa = Basics.extractFromConfigurations(NumberOfElectrons(), newConfs);    wa = unique(wa)
    if  length(wa) != 1   ||   wa[1] != noe        
        Basics.displayConfigurations(stdout, newConfs)
        error("Configurations with unexpected NoElectrons are generated.") 
    end 
    
    return( newConfs )
end


"""
`Basics.generateConfigurations(theme::Basics.ForPhotoRecombination, confs::Array{Configuration,1})`  
    ... generates a list of non-relativistic configurations in which an additional electron is captured into one of the 
        theme.intoShells, so that each configuration has now one electron more. This is done for all given configurations.
        Cf. also Basics.ForElectronCapture; a list of possible newConfs::Array{Configuration,1} is returned.
"""
function Basics.generateConfigurations(theme::Basics.ForPhotoRecombination, confs::Array{Configuration,1}) 
    # Add one electron in each of the theme.intoShells    
    newConfs = Basics.generateConfigurations(AddElectrons(1, theme.intoShells), confs)
    
    return( newConfs )
end


"""
`Basics.generateConfigurations(theme::Basics.ForRasExcitations, confs::Array{Configuration,1})`  
    ... generates a list of non-relativistic configurations in which single, double, triple, etc. excitations are considered 
        for each given configuration due to the boolean values theme.singles, themes.doubles, ... . Excitation of electrons
        are considered from theme.fromShells to the theme.intoShells. The procedure terminates if 
            (i) configurations with a different number of electrons are given. 
        A list of possible newConfs::Array{Configuration,1} is returned.
"""
function Basics.generateConfigurations(theme::Basics.ForRasExcitations, confs::Array{Configuration,1})
    # Check that all configurations have the same number of electrons
    wa = Basics.extractFromConfigurations(NumberOfElectrons(), confs)
    if  length(wa) != 1   error("Configurations with different NoElectrons are given.")   else   noe = wa[1]    end 
    
    newConfs = deepcopy(confs)
    if  theme.se   seConfs  = Basics.generateConfigurations(ExciteElectrons(1, theme.fromShells, theme.intoShells), confs)
        newConfs = Basics.merge(newConfs, seConfs)   end
    
    if  theme.de   deConfs  = Basics.generateConfigurations(ExciteElectrons(2, theme.fromShells, theme.intoShells), confs) 
        newConfs = Basics.merge(newConfs, deConfs)   end
    
    if  theme.te   teConfs  = Basics.generateConfigurations(ExciteElectrons(3, theme.fromShells, theme.intoShells), confs)
        newConfs = Basics.merge(newConfs, teConfs)   end
    
    if  theme.qe   qeConfs  = Basics.generateConfigurations(ExciteElectrons(4, theme.fromShells, theme.intoShells), confs) 
        newConfs = Basics.merge(newConfs, qeConfs)   end
    newConfs = unique(newConfs)
  
    return( newConfs )
end


"""
`Basics.generateConfigurations(theme::Basics.ForStepwiseDecay, confs::Array{Configuration,1})`  
    ... generates a list of non-relativistic configurations that model the stepwise decay of an inner-shell hole by steptwise 
        photoemission and autoionization processes until a certain number of electrons is released. The procedure first determines 
        the most inner-shell hole and continues the stepwise decay until theme.maximallyReleased is reached or no further decay 
        can occur. Intra-shell processes can be considered but are not reflected in the generated configurations. A list of possible newConfs::Array{Configuration,1} is returned.
"""
function Basics.generateConfigurations(theme::Basics.ForStepwiseDecay, confs::Array{Configuration,1})
    newConfs = Configuration[];   currentConfs = deepcopy(confs)
    # Test for equal number of electrons
    wa = Basics.extractFromConfigurations(NumberOfElectrons(), confs)
    if  length(wa) != 1   error("Configurations with different NoElectrons are given.")   else   noe = wa[1]    end 
    
    neMax = max(1, theme.maximallyReleased)
    for  ne = 1:neMax
        newConfigs = Configuration[]
        for  conf in currentConfs
            valenceShells  = Basics.extractFromConfiguration(Basics.ValenceShells(),   conf)
            # Determine configuration due to an de-excitation of an electron within valenceShells
            if     length(valenceShells) < 2
            else 
                addConfigs = Basics.generateConfigurations(AddElectrons(1, valenceShells[1:1]), [conf])
                newConfsPE = Basics.generateConfigurations(RemoveElectrons(1, valenceShells[2:end]), addConfigs)
                newConfsAI = Basics.generateConfigurations(RemoveElectrons(2, valenceShells[2:end]), addConfigs)
                # Special treatment is necessary if no electron is released
                if   theme.maximallyReleased == 0   append!(newConfigs, newConfsPE)   
                else                                append!(newConfigs, newConfsPE, newConfsAI)
                end
            end
        end
        newConfigs   = unique(newConfigs)
        currentConfs = deepcopy(newConfigs)
        append!(newConfs, newConfigs)
    end 
    #
    newConfs = unique(newConfs)
    wb = Basics.extractFromConfigurations(NumberOfElectrons(), newConfs);   @show wb
    
    return( newConfs )
end

    
"""
`Basics.generateConfigurations(theme::Basics.GeneralizedConfigurations, confs::Array{Configuration,1})`  
    ... to generate all generalized configurations that are related to the given set of configurations.
        At present, neither the detailed definition nor their unique generation from a set of non-relativistic configurations has 
        been worked out in good detail. This this therefore a dummy routine whose code still need to be worked out. 
        An error is returned in all cases.
"""
function Basics.generateConfigurations(theme::Basics.GeneralizedConfigurations, confs::Array{Configuration,1})
    error("Not yet implemented")
    
    return( nothing )        
end


"""
`Basics.generateConfigurations(theme::Basics.SuperConfiguration, NoElectrons::Int64, fixedShells::Dict{Shell,Int64}, 
                               mutableShells::Array{Shell,1})`  
    ... to generate all configurations that are described by the superconfiguration as encoded by the NoElectrons, fixedShells 
        as well as mutableShells. The fixedShells are tyically taken from a configuration conf.shells, while 
        k = NoElectrons - conf.NoElectrons > 0 are the number of electrons, which are distributed among the mutable shells.  
        The routine terminates if k <= 0  or if any shell appears both in fixedShells and mutableShells simultaneously. 
        A confList::Array{Configuration,1}.
"""
function Basics.generateConfigurations(theme::Basics.SuperConfiguration, NoElectrons::Int64, fixedShells::Dict{Shell,Int64}, 
                                       mutableShells::Array{Shell,1})
    function distribute_items(m::Int, k::Int)
        # Function to generate all distributions of m identical items into k boxes
        results = Int[] |> x -> Vector{Vector{Int}}()
        
        function backtrack(remaining::Int, box::Int, current::Vector{Int})
            if box == k    # Last box gets whatever is left
                push!(results, vcat(current, [remaining]))
            else
                for i in 0:remaining   backtrack(remaining - i, box + 1, vcat(current, [i]))   end
            end
        end
        
        backtrack(m, 1, Int[])
        return( results )
    end
    
    confList = Configuration[];         nc = 0;   newShells = deepcopy(fixedShells)
    for  (shell, occ) in fixedShells    nc = nc + fixedShells[shell]   end
    k = NoElectrons - nc;               if  k <= 0    error("Improper number of electrons.")    end
    
    # Now distribute k items over length(fixedShells) items
    itemsList = distribute_items(k, length(mutableShells))
    for  items  in  itemsList
        newShells = deepcopy(fixedShells);    newNoElectrons = nc;   breakOut = false
        for  (it, item)  in  enumerate(items)
            mShell = mutableShells[it];   maxNoe = 2*(2 * mShell.l + 1)
            if   item > maxNoe    breakOut = true;    break
            else                  newShells[mShell] = item;   newNoElectrons = newNoElectrons + item
            end 
        end
        if      breakOut
        elseif  newNoElectrons != NoElectrons    error("stop a")    
        else                                     push!(confList, Configuration(newShells, NoElectrons))
        end
    end
    confList = unique(confList)
    
    return( confList )        
end
    

"""
`Basics.generateConfigurations(theme::Basics.RelativisticConfigurations, conf::Configuration)`  
    ... to split/decompose a (single) non-relativistic configuration into an list of relativistic ConfigurationR[]. The proper 
        occupuation of the relativistic subshells is taken into account. A confList::Array{ConfigurationR,1} is returned.
"""
function Basics.generateConfigurations(theme::Basics.RelativisticConfigurations, conf::Configuration)
    subshellList = Subshell[]
    confList     = ConfigurationR[]

    initialize = true;    NoElectrons = 0
    for (k, occ)  in conf.shells
        NoElectrons     = NoElectrons + occ
        wa              = Basics.shellSplitOccupation(k, occ)
        subshellListNew = Dict{Subshell,Int64}[]

        if  initialize
            subshellListNew = wa;    initialize = false
        else
            for  s in 1:length(subshellList)
                for  a in 1:length(wa)
                    wb = Base.merge( subshellList[s], wa[a] )
                    push!(subshellListNew, wb)
                end
            end
        end
        subshellList = deepcopy(subshellListNew)
    end
    
    for subsh in subshellList
        wa = ConfigurationR(subsh, NoElectrons)
        push!(confList, wa)
    end

    return( confList )
end


#################################################################################################################################
#################################################################################################################################


"""
`Basics.merge(aList::Array{Configuration,1}, bList::Array{Configuration,1}, ...)`  
    ... to merge two (or more) configuration list into a single list and to unify them. 
        A cList::Array{Configuration,1} is returned.
"""
function Basics.merge(lists...)
    cList = Configuration[]
    for list in lists
        append!(cList, list);      cList = Base.unique(cList)
    end
    return( cList )
end


#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

# These routines should be later "removed" from the code since they can be readily replaced by other, more transparent
# routines.

"""
`Basics.displayConfigurations(Z::Float64, confs::Array{Configuration,1}; sa::String="")` 
    ... group & display the configuration list into sublists with the same No. of electrons; this lists are displayed together 
        with an estimated total energy. An ordered confList::Array{Configuration,1} is returned with configurations of decreasing
        number of electrons.
"""
function Basics.displayConfigurations(Z::Float64, confs::Array{Configuration,1}; sa::String="")
    #  Replace in code: Cascade.groupDisplayConfigurationList(Z::Float64, confs::Array{Configuration,1}; sa::String="") ... just renamed in the code
    
    minNoElectrons = 1000;   maxNoElectrons = 0  
    for  conf in confs
        minNoElectrons = min(minNoElectrons, conf.NoElectrons)
        maxNoElectrons = max(maxNoElectrons, conf.NoElectrons)
    end
    #
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    #
    println("\n* Electron configuration(s) used:    " * sa)
    ## @warn "*** Limit to just 4 configurations for each No. of electrons. ***"                       ## delete nxx
    if  printSummary   println(iostream, "\n* Electron configuration(s) used:    " * sa)    end
    confList = Configuration[];   nc = 0
    for  n = maxNoElectrons:-1:minNoElectrons
        nxx = 0                                                                                        ## delete nxx
        println("\n  Configuration(s) with $n electrons:")
        if  printSummary   println(iostream, "\n    Configuration(s) with $n electrons:")      end
        nd = 0
        for  conf in confs  nd = max(nd, length("      " * string(conf)))   end
        for  conf in confs
            if n == conf.NoElectrons  
                ## nxx = nxx + 1;    if nxx > 4   break    end                                         ## delete nxx
                nc = nc + 1
                push!(confList, conf ) 
                ##x wa = Semiempirical.estimate("binding energy: XrayDataBooklet", round(Int64, Z), conf);  
                wa = -Empirical.totalEnergy(round(Int64, Z), conf, data = PeriodicTable.XrayDataBooklet() )
                wa = Defaults.convertUnits("energy: from atomic", wa)
                ##x if  Z > 36.0    wa = 0.
                ##x else            wa = Semiempirical.estimate("binding energy: XrayDataBooklet", round(Int64, Z), conf);    
                ##x                 wa = Defaults.convertUnits("energy: from atomic", wa)
                ##x end
                sb = "   av. BE = "  * string( round(-wa) ) * "  " * TableStrings.inUnits("energy")
                sd = "      " * string(conf) * "                                "
                println(sd[1:nd+3] * sb * "      ($nc)" )
                if  printSummary   println(iostream, sd[1:nd+3] * sb * "      ($nc)")      end
            end  
        end
    end
    
    println("\n  A total of $nc configuration have been defined for this " * sa * "cascade, and selected configurations could be " *
            "removed here:  [currently not supported]")
    if  printSummary   println(iostream, "\n* A total of $nc configuration have been defined for this cascade, and selected " *
                                            "configurations could be removed here:  [currently not supported]")      end
    return( confList )
end
