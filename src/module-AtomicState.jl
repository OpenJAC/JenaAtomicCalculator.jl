
"""
`module  JAC.AtomicState`  
	... a submodel of JAC that contains all methods to set-up and process atomic representations.
"""
module AtomicState


## using Interact
using  ..Basics, ..ManyElectron, ..Nuclear, ..Radial

export  MeanFieldSettings, MeanFieldBasis, MeanFieldMultiplet, OneElectronSettings, OneElectronSpectrum, CiSettings, CiExpansion, 
        RasSettings, RasStep, RasExpansion, GreenSettings, GreenChannel, GreenExpansion, Representation


"""
`abstract type AtomicState.AbstractRepresentationType` 
    ... defines an abstract type and a number of data types to work with and distinguish different atomic representations; see also:
    
    + struct MeanFieldBasis       ... to represent (and generate) a mean-field basis and, especially, a set of 
                                        (mean-field) orbitals.
    + struct MeanFieldMultiplet   ... to represent (and generate) a mean-field multiplet, based on (mean-field) orbitals.
    + struct OneElectronSpectrum  ... to represent (and generate) a one-electron spectrum for the mean-field 
                                        potential of refConfigs.
    + struct CiExpansion          ... to represent (and generate) a configuration-interaction representation.
    + struct RasExpansion         ... to represent (and generate) a restricted active-space representation.
    + struct GreenExpansion       ... to represent (and generate) an approximate (many-electron) Green functions.
"""
abstract type  AbstractRepresentationType                           end


### Mean-field basis ############################################################################################
"""
`struct  AtomicState.MeanFieldSettings`  
    ... a struct for defining the settings for a mean-field basis (orbital) representation.

    + scField           ::AbstractScField   
        ... Specify the (mean) self-consistent field as DFSField() or HSField(); note that not all AbstractScField's 
            allowed here.
"""
struct  MeanFieldSettings
    scField              ::AbstractScField
end

"""
`AtomicState.MeanFieldSettings()`  ... constructor for setting the default values.
"""
function MeanFieldSettings()
    MeanFieldSettings(DFSField())
end


# `Base.show(io::IO, settings::MeanFieldSettings)`  ... prepares a proper printout of the settings::MeanFieldSettings.
function Base.show(io::IO, settings::MeanFieldSettings)
        println(io, "scField:                  $(settings.scField)  ")
end


"""
`struct  AtomicState.MeanFieldBasis  <:  AbstractRepresentationType`  
    ... a struct to represent (and generate) a mean-field orbital basis.

    + settings         ::AtomicState.MeanFieldSettings      ... Settings for the given mean-field orbital basis
"""
struct   MeanFieldBasis  <:  AbstractRepresentationType
    settings           ::AtomicState.MeanFieldSettings
end


# `Base.string(basis::MeanFieldBasis)`  ... provides a String notation for the variable basis::MeanFieldBasi.
function Base.string(basis::MeanFieldBasis)
    if !(basis.settings.scField  in [DFSField(), HSField()])
        error("A MeanFieldBasis presently supports only a DFSField() or HSField() but received:  $(basis.settings.scField)")
    end
    #
    sa = "Mean-field orbital basis for a $(basis.settings.scField) SCF field:"
    return( sa )
end


# `Base.show(io::IO, basis::MeanFieldBasis)`  ... prepares a proper printout of the basis::MeanFieldBasis.
function Base.show(io::IO, basis::MeanFieldBasis)
    sa = Base.string(basis);       print(io, sa, "\n")
end


"""
`struct  AtomicState.MeanFieldMultiplet  <:  AbstractRepresentationType`  
    ... a struct to represent (and generate) a mean-field orbital basis and multiplet.

    + settings         ::AtomicState.MeanFieldSettings      ... Settings for the given mean-field orbital basis and multiplet.
"""
struct   MeanFieldMultiplet  <:  AbstractRepresentationType
    settings           ::AtomicState.MeanFieldSettings
end


# `Base.string(basis::MeanFieldMultiplet)`  ... provides a String notation for the variable basis::MeanFieldMultiplet.
function Base.string(basis::MeanFieldMultiplet)
    #
    sa = "Mean-field multiplet for a $(basis.settings.scField) SCF field:"
    return( sa )
end


# `Base.show(io::IO, basis::MeanFieldMultiplet)`  ... prepares a proper printout of the basis::MeanFieldMultiplet.
function Base.show(io::IO, basis::MeanFieldMultiplet)
    sa = Base.string(basis);       print(io, sa, "\n")
end



### Spectrum: One-electron spectrum for given (reference) configurations ####################################################
"""
`struct  AtomicState.OneElectronSettings`  
    ... a struct for defining the settings for a one-electron spectrum.

    + nMax                 ::Int64            ... maximum principal quantum number.
    + lValues              ::Array{Int64,1}   ... l-values (partial waves) for which orbitals are to be generated.
    + levelSelectionMean   ::LevelSelection   ... Level(s) of the mean-field to generate the underlying atomic potential
"""
struct  OneElectronSettings
    nMax                   ::Int64
    lValues                ::Array{Int64,1}
    levelSelectionMean     ::LevelSelection
end

"""
`AtomicState.OneElectronSettings()`  ... constructor for setting the default values.
"""
function OneElectronSettings()
    OneElectronSettings(10, [0, 1], LevelSelection(true, indices=[1]))
end


"""
`AtomicState.OneElectronSettings(settings::AtomicState.OneElectronSettings;`

        nMax::Union{Nothing,Int64}=nothing,                             lValues::Union{Nothing,Array{Int64,1}}=nothing,
        levelSelectionMean::Union{Nothing,LevelSelection}=nothing)
                    
    ... constructor for modifying the given OneElectronSettings by 'overwriting' the explicitly selected parameters.
"""
function OneElectronSettings(settings::AtomicState.OneElectronSettings;
    nMax::Union{Nothing,Int64}=nothing,                             lValues::Union{Nothing,Array{Int64,1}}=nothing,
    levelSelectionMean::Union{Nothing,LevelSelection}=nothing)
    
    if  nMax               == nothing   nMaxx               = settings.nMax                else   nMaxx               = nMax               end 
    if  lValues            == nothing   lValuesx            = settings.lValues             else   lValuesx            = lValues            end 
    if  levelSelectionMean == nothing   levelSelectionMeanx = settings.levelSelectionMean  else   levelSelectionMeanx = levelSelectionMean end 
    
    OneElectronSettings( nMaxx, lValuesx, levelSelectionMeanx)
end


# `Base.show(io::IO, settings::OneElectronSettings)`  ... prepares a proper printout of the settings::OneElectronSettings.
function Base.show(io::IO, settings::OneElectronSettings)
        println(io, "nMax:                 $(settings.nMax)  ")
        println(io, "lValues:              $(settings.lValues)  ")
        println(io, "levelSelectionMean:   $(settings.levelSelectionMean)  ")
end


"""
`struct  AtomicState.OneElectronSpectrum  <:  AbstractRepresentationType`  
    ... a struct to represent (and generate) a one-electron spectrum for the mean-field potential of given levels 
        from a set of reference configurations.

    + settings         ::AtomicState.OneElectronSettings         ... Settings for the given OneElectronSpectrum
"""
struct  OneElectronSpectrum     <:  AbstractRepresentationType
    settings           ::AtomicState.OneElectronSettings
end


# `Base.string(spectrum::OneElectronSpectrum)`  ... provides a String notation for the variable spectrum::OneElectronSpectrum.
function Base.string(spectrum::OneElectronSpectrum)
    sa = "One-electron spectrum for given atomic potential:"
    return( sa )
end


# `Base.show(io::IO, spectrum::OneElectronSpectrum)`  
#       ... prepares a proper printout of the (individual step of computations) spectrum::OneElectronSpectrum.
function Base.show(io::IO, spectrum::OneElectronSpectrum)
    sa = Base.string(spectrum);       print(io, sa, "\n")
    println(io, "... and the current settings:")
    println(io, "$(spectrum.settings)  ")
end



### RAS: Restricted-Active-Space expansions ##################################################################################
"""
`struct  AtomicState.RasSettings`  
    ... a struct for defining the settings for a restricted active-space computations.

    + levelsScf            ::Array{Int64,1}         ... Levels on which the optimization need to be carried out.
    + maxIterationsScf     ::Int64                  ... maximum number of SCF iterations in each RAS step.
    + accuracyScf          ::Float64                ... convergence criterion for the SCF field.
    + eeInteractionCI      ::AbstractEeInteraction  ... logical flag to include Breit interactions.
    + levelSelectionCI     ::LevelSelection         ... Specifies the selected levels, if any.
"""
struct  RasSettings
    levelsScf              ::Array{Int64,1}
    maxIterationsScf       ::Int64  
    accuracyScf            ::Float64 
    eeInteractionCI        ::AbstractEeInteraction 
    levelSelectionCI       ::LevelSelection
end

"""
`AtomicState.RasSettings()`  ... constructor for setting the default values.
"""
function RasSettings()
    RasSettings(Int64[1], 24, 1.0e-6, CoulombInteraction(), LevelSelection() )
end


# `Base.show(io::IO, settings::RasSettings)`  ... prepares a proper printout of the settings::RasSettings.
function Base.show(io::IO, settings::RasSettings)
        println(io, "levelsScf:            $(settings.levelsScf)  ")
        println(io, "maxIterationsScf:     $(settings.maxIterationsScf)  ")
        println(io, "accuracyScf:          $(settings.accuracyScf)  ")
        println(io, "eeInteractionCI:      $(settings.eeInteractionCI)  ")
        println(io, "levelSelectionCI:     $(settings.levelSelectionCI)  ")
end


"""
`struct  AtomicState.RasStep`  
    ... specifies an individual step of a (relativistic) restricted active space computation for a set of levels. This struct 
        comprises all information to generate the orbital basis and to perform the associated SCF and multiplet computations for a 
        selected number of levels.

    + seFrom            ::Array{Shell,1}        ... Single-excitations from shells   [sh_1, sh_2, ...]
    + seTo              ::Array{Shell,1}        ... Single-excitations to shells  [sh_1, sh_2, ...]
    + deFrom            ::Array{Shell,1}        ... Double-excitations from shells   [sh_1, sh_2, ...]
    + deTo              ::Array{Shell,1}        ... Double-excitations to shells  [sh_1, sh_2, ...]
    + teFrom            ::Array{Shell,1}        ... Triple-excitations from shells   [sh_1, sh_2, ...]
    + teTo              ::Array{Shell,1}        ... Triple-excitations to shells  [sh_1, sh_2, ...]
    + qeFrom            ::Array{Shell,1}        ... Quadrupole-excitations from shells   [sh_1, sh_2, ...]
    + qeTo              ::Array{Shell,1}        ... Quadrupole-excitations to shells  [sh_1, sh_2, ...]
    + frozenShells      ::Array{Shell,1}        ... List of shells that are kept 'frozen' in this step.
    + constraints       ::Array{String,1}       ... List of Strings to define 'constraints/restrictions' 
                                                    to the generated CSF basis.
"""
struct  RasStep
    seFrom              ::Array{Shell,1}
    seTo                ::Array{Shell,1}
    deFrom              ::Array{Shell,1}
    deTo                ::Array{Shell,1}
    teFrom              ::Array{Shell,1}
    teTo                ::Array{Shell,1}
    qeFrom              ::Array{Shell,1}
    qeTo                ::Array{Shell,1}
    frozenShells        ::Array{Shell,1}
    constraints         ::Array{String,1}
end

"""
`AtomicState.RasStep()`  ... constructor for an 'empty' instance of a variable::AtomicState.RasStep
"""
function RasStep()
    RasStep(Shell[], Shell[], Shell[], Shell[], Shell[], Shell[], Shell[], Shell[],    Shell[], String[])
end


"""
`AtomicState.RasStep(rasStep::AtomicState.RasStep;`

                    seFrom::Array{Shell,1}=Shell[], seTo::Array{Shell,1}=Shell[], 
                    deFrom::Array{Shell,1}=Shell[], deTo::Array{Shell,1}=Shell[], 
                    teFrom::Array{Shell,1}=Shell[], teTo::Array{Shell,1}=Shell[], 
                    qeFrom::Array{Shell,1}=Shell[], qeTo::Array{Shell,1}=Shell[], 
                    frozen::Array{Shell,1}=Shell[], constraints::Array{String,1}=String[]  
                    
    ... constructor for modifying the given rasStep by specifying all excitations, frozen shells and 
        constraints optionally.
"""
function RasStep(rasStep::AtomicState.RasStep;
                    seFrom::Array{Shell,1}=Shell[], seTo::Array{Shell,1}=Shell[], 
                    deFrom::Array{Shell,1}=Shell[], deTo::Array{Shell,1}=Shell[], 
                    teFrom::Array{Shell,1}=Shell[], teTo::Array{Shell,1}=Shell[], 
                    qeFrom::Array{Shell,1}=Shell[], qeTo::Array{Shell,1}=Shell[], 
                    frozen::Array{Shell,1}=Shell[], constraints::Array{String,1}=String[])
    if  seFrom == Shell[]   sxFrom = rasStep.seFrom   else      sxFrom = seFrom      end
    if  seTo   == Shell[]   sxTo   = rasStep.seTo     else      sxTo   = seTo        end
    if  deFrom == Shell[]   dxFrom = rasStep.deFrom   else      dxFrom = deFrom      end
    if  deTo   == Shell[]   dxTo   = rasStep.deTo     else      dxTo   = deTo        end
    if  teFrom == Shell[]   txFrom = rasStep.teFrom   else      txFrom = teFrom      end
    if  teTo   == Shell[]   txTo   = rasStep.teTo     else      txTo   = teTo        end
    if  qeFrom == Shell[]   qxFrom = rasStep.qeFrom   else      qxFrom = qeFrom      end
    if  qeTo   == Shell[]   qxTo   = rasStep.qeTo     else      qxTo   = qeTo        end
    if  frozen == Shell[]        frozx  = rasStep.frozenShells   else      frozx  = frozen             end
    if  constraints ==String[]   consx  = rasStep.constraints    else      consx  = constraints        end
    
    RasStep( sxFrom, sxTo, dxFrom, dxTo, txFrom, txTo, qxFrom, qxTo, frozx, consx)
end


# `Base.string(step::AtomicState.RasStep)`  ... provides a String notation for the variable step::AtomicState.RasStep.
function Base.string(step::AtomicState.RasStep)
    sa = "\nCI or RAS step with $(length(step.frozenShells)) (explicitly) frozen shell(s): $(step.frozenShells)  ... and virtual excitations"
    return( sa )
end


# `Base.show(io::IO, step::AtomicState.RasStep)`  ... prepares a proper printout of the (individual step of computations) step::AtomicState.RasStep.
function Base.show(io::IO, step::AtomicState.RasStep)
    sa = Base.string(step);                 print(io, sa, "\n")
    if  length(step.seFrom) > 0
        sa = "   Singles from:          { ";      for  sh in step.seFrom   sa = sa * string(sh) * ", "  end
        sa = sa[1:end-2] * " }   ... to { ";      for  sh in step.seTo     sa = sa * string(sh) * ", "  end;   
        sa = sa[1:end-2] * " }";            print(io, sa, "\n")
    end   
    if  length(step.deFrom) > 0
        sa = "   Doubles from:          { ";      for  sh in step.deFrom   sa = sa * string(sh) * ", "  end
        sa = sa[1:end-2] * " }   ... to { ";      for  sh in step.deTo     sa = sa * string(sh) * ", "  end;   
        sa = sa[1:end-2] * " }";            print(io, sa, "\n")
    end   
    if  length(step.teFrom) > 0
        sa = "   Triples from:          { ";      for  sh in step.teFrom   sa = sa * string(sh) * ", "  end
        sa = sa[1:end-2] * " }   ... to { ";      for  sh in step.teTo     sa = sa * string(sh) * ", "  end;   
        sa = sa[1:end-2] * " }";            print(io, sa, "\n")
    end   
    if  length(step.qeFrom) > 0
        sa = "   Quadruples from:       { ";      for  sh in step.qeFrom   sa = sa * string(sh) * ", "  end
        sa = sa[1:end-2] * " }   ... to { ";      for  sh in step.qeTo     sa = sa * string(sh) * ", "  end;   
        sa = sa[1:end-2] * " }";            print(io, sa, "\n")
    end   
end


"""
`struct  AtomicState.RasExpansion    <:  AbstractRepresentationType`  
    ... a struct to represent (and generate) a restricted active-space representation.

    + symmetries       ::Array{LevelSymmetry,1}         ... Symmetries of the levels/CSF in the many-electron basis.
    + NoElectrons      ::Int64                          ... Number of electrons.
    + steps            ::Array{AtomicState.RasStep,1}   ... List of SCF steps that are to be done in this model 
                                                            computation.
    + settings         ::AtomicState.RasSettings        ... Settings for the given RAS computation
"""
struct         RasExpansion    <:  AbstractRepresentationType
    symmetries         ::Array{LevelSymmetry,1}
    NoElectrons        ::Int64
    steps              ::Array{AtomicState.RasStep,1}
    settings           ::AtomicState.RasSettings 
    end


"""
`AtomicState.RasExpansion()`  ... constructor for an 'empty' instance of the a variable::AtomicState.RasExpansion
"""
function RasExpansion()
    RasExpansion([Basics.LevelSymmetry(0, Basics.plus)], 0, 0, AtomicState.RasStep[], AtomicState.RasSettings())
end


# `Base.string(expansion::RasExpansion)`  ... provides a String notation for the variable expansion::RasExpansion.
function Base.string(expansion::RasExpansion)
    sa = "RAS expansion for symmétry $(expansion.symmetries) and with $(length(expansion.steps)) steps:"
    return( sa )
end


# `Base.show(io::IO, expansion::RasExpansion)`  ... prepares a proper printout of the (individual step of computations) expansion::RasExpansion.
function Base.show(io::IO, expansion::RasExpansion)
    sa = Base.string(expansion);       print(io, sa, "\n")
    println(io, "$(expansion.steps)  ")
    println(io, "... and the current settings:")
    println(io, "$(expansion.settings)  ")
end



### CI: Configuration Interaction expansions ##################################################################################
"""
`struct  AtomicState.CiSettings`  
    ... a struct for defining the settings for a configuration-interaction (CI) expansion.

    + eeInteractionCI      ::AbstractEeInteraction   ... Specifies the treatment of the e-e interaction.
    + levelSelectionCI     ::LevelSelection          ... Specifies the selected levels, if any.
"""
struct  CiSettings
    eeInteractionCI        ::AbstractEeInteraction 
    levelSelectionCI       ::LevelSelection
end

"""
`AtomicState.CiSettings()`  ... constructor for setting the default values.
"""
function CiSettings()
    CiSettings(CoulombInteraction(), LevelSelection() )
end


"""
`AtomicState.CiSettings(settings::AtomicState.CiSettings;`

        eeInteractionCI::Union{Nothing,AbstractEeInteraction}=nothing,    levelSelectionCI::Union{Nothing,LevelSelection}=nothing)
                    
    ... constructor for modifying the given CiSettings by 'overwriting' the explicitly selected parameters.
"""
function CiSettings(settings::AtomicState.CiSettings;
    eeInteractionCI::Union{Nothing,AbstractEeInteraction}=nothing,        levelSelectionCI::Union{Nothing,LevelSelection}=nothing)
    
    if  eeInteractionCI     == nothing   eeInteractionCIx      = settings.eeInteractionCI       else   eeInteractionCIx      = eeInteractionCI      end 
    if  levelSelectionCI    == nothing   levelSelectionCIx     = settings.levelSelectionCI      else   levelSelectionCIx     = levelSelectionCI     end 
    
    CiSettings( eeInteractionCIx, levelSelectionCIx)
end


# `Base.show(io::IO, settings::CiSettings)`  ... prepares a proper printout of the settings::CiSettings.
function Base.show(io::IO, settings::CiSettings)
        println(io, "eeInteractionCI:          $(settings.eeInteractionCI)  ")
        println(io, "levelSelectionCI:         $(settings.levelSelectionCI)  ")
end


"""
`struct  AtomicState.CiExpansion  <:  AbstractRepresentationType`  
    ... a struct to represent (and generate) a configuration-interaction representation.

    + applyOrbitals    ::Dict{Subshell, Orbital}
    + excitations      ::AtomicState.RasStep            ... Excitations beyond refConfigs.
    + settings         ::AtomicState.CiSettings         ... Settings for the given CI expansion
"""
struct         CiExpansion     <:  AbstractRepresentationType
    applyOrbitals      ::Dict{Subshell, Orbital}
    excitations        ::AtomicState.RasStep
    settings           ::AtomicState.CiSettings
end


# `Base.string(expansion::CiExpansion)`  ... provides a String notation for the variable expansion::CiExpansion.
function Base.string(expansion::CiExpansion)
    sa = "CI expansion with (additional) excitations:"
    return( sa )
end


# `Base.show(io::IO, expansion::CiExpansion)`  ... prepares a proper printout of the (individual step of computations) expansion::CiExpansion.
function Base.show(io::IO, expansion::CiExpansion)
    sa = Base.string(expansion);       print(io, sa, "\n")
    println(io, "$(expansion.excitations)  ")
    println(io, "... and the current settings:")
    println(io, "$(expansion.settings)  ")
end


    
### Green function representation ##################################################################################

"""
`abstract type AtomicState.AbstractGreenApproach` 
    ... defines an abstract and a number of singleton types for approximating a many-electron Green
        function expansion for calculating second-order processes.

    + struct SingleCSFwithoutCI        
    ... to approximate the many-electron multiplets (gMultiplet) for every chosen level symmetry by dealing with each CSF 
        independently, and without any configuration interaction. This is a fast but also very rough approximation.
        
    + struct CoreSpaceCI                
    ... to approximate the many-electron multiplets (gMultiplet) by taking the electron-electron interaction between the 
        bound-state orbitals into account.
        
    + struct DampedSpaceCI                
    ... to approximate the many-electron multiplets (gMultiplet) by taking the electron-electron interaction for all, the 
        bound and free-electron, orbitals into account but by including a damping factor e^{tau*r}, tau > 0 into the 
        electron densities rho_ab (r) --> rho_ab (r) * e^{tau*r}
"""
abstract type  AbstractGreenApproach                                  end
struct         SingleCSFwithoutCI  <:  AtomicState.AbstractGreenApproach   end
struct         CoreSpaceCI         <:  AtomicState.AbstractGreenApproach   end
struct         DampedSpaceCI       <:  AtomicState.AbstractGreenApproach   end


"""
`struct  AtomicState.GreenSettings`  
    ... defines a type for defining the details and parameters of the approximate Green (function) expansion.

    + nMax                     ::Int64            ... maximum principal quantum numbers of (single-electron) 
                                                        excitations that are to be included into the representation.
    + lValues                  ::Array{Int64,1}   ... List of (non-relativistic) orbital angular momenta for which
                                                        (single-electron) excitations are to be included.
    + dampingTau               ::Float64          ... factor tau (> 0.) that is used to 'damp' the one- and two-electron
                                                        interactions strength: exp( - tau * r)
    + printBefore              ::Bool             ... True if a short overview is to be printed before. 
    + levelSelection           ::LevelSelection   ... Specifies the selected levels, if any.
"""
struct GreenSettings 
    nMax                       ::Int64
    lValues                    ::Array{Int64,1}
    dampingTau                 ::Float64
    printBefore                ::Bool 
    levelSelection             ::LevelSelection
end 


"""
`AtomicState.GreenSettings()`  ... constructor for an `empty` instance of AtomicState.GreenSettings.
"""
function GreenSettings()
    Settings( 0, Int64[], 0., false, LevelSelection() )
end


# `Base.show(io::IO, settings::AtomicState.GreenSettings)`  ... prepares a proper printout of the variable settings::AtomicState.GreenSettings.
function Base.show(io::IO, settings::AtomicState.GreenSettings) 
    println(io, "nMax:                     $(settings.nMax)  ")
    println(io, "lValues:                  $(settings.lValues)  ")
    println(io, "dampingTau:               $(settings.dampingTau)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "levelSelection:           $(settings.levelSelection)  ")
end


"""
`struct  AtomicState.GreenChannel`  
    ... defines a type for a single symmetry channel of an (approximate) Green (function) expansion.

    + symmetry          ::LevelSymmetry    ... Level symmetry of this part of the representation.
    + gMultiplet        ::Multiplet        ... Multiplet of (scattering) levels of this symmetry.
"""
struct GreenChannel 
    symmetry            ::LevelSymmetry
    gMultiplet          ::Multiplet
end   


"""
`AtomicState.GreenChannel()`  ... constructor for an `empty` instance of AtomicState.GreenChannel.
"""
function GreenChannel()
    GreenChannel( LevelSymmetry(0, Basics.plus), ManyElectron.Multiplet() )
end


# `Base.show(io::IO, channel::AtomicState.GreenChannel)`  ... prepares a proper printout of the variable channel::AtomicState.GreenChannel.
function Base.show(io::IO, channel::AtomicState.GreenChannel) 
    println(io, "symmetry:                $(channel.symmetry)  ")
    println(io, "gMultiplet:              $(channel.gMultiplet)  ")
end


"""
`struct  AtomicState.GreenExpansion  <:  AbstractRepresentationType`  
    ... defines a type to keep an (approximate) Green (function) expansion that is associated with a given set of reference
        configurations.

    + approach          ::AtomicState.AbstractGreenApproach  ... Approach used to approximate the representation.
    + excitationScheme  ::Basics.AbstractExcitationScheme    ... Applied excitation scheme w.r.t. refConfigs. 
    + levelSymmetries   ::Array{LevelSymmetry,1}             ... Total symmetries J^P to be included into Green expansion.
    + NoElectrons       ::Int64                              ... Number of electrons.
    + settings          ::AtomicState.GreenSettings          ... settings for the Green (function) expansion.
"""
struct GreenExpansion  <:  AbstractRepresentationType
    approach            ::AtomicState.AbstractGreenApproach
    excitationScheme    ::Basics.AbstractExcitationScheme 
    levelSymmetries     ::Array{LevelSymmetry,1}
    NoElectrons         ::Int64 
    settings            ::AtomicState.GreenSettings
end   


"""
`AtomicState.GreenExpansion()`  ... constructor for an `empty` instance of AtomicState.GreenExpansion.
"""
function GreenExpansion()
    GreenExpansion( AtomicState.SingleCSFwithoutCI(), Basics.NoExcitationScheme(), LevelSymmetry[], 0, AtomicState.GreenSettings())
end


# `Base.string(expansion::GreenExpansion)`  ... provides a String notation for the variable expansion::GreenExpansion.
function Base.string(expansion::GreenExpansion)
    sa = "Green (function) expansion in $(expansion.approach) approach and for excitation scheme  $(expansion.excitationScheme)," *
            "\nincluding (Green function channels with) symmetries $(expansion.levelSymmetries):"
    return( sa )
end


# `Base.show(io::IO, expansion::GreenExpansion)`  ... prepares a proper printout of the (individual step of computations) expansion::GreenExpansion.
function Base.show(io::IO, expansion::GreenExpansion)
    sa = Base.string(expansion);       print(io, sa, "\n")
    println(io, "... and the current settings:")
    println(io, "$(expansion.settings)  ")
end



### Representation ##################################################################################
"""
`struct  AtomicState.Representation`  
    ... a struct for defining an atomic state representation. Such representations often refer to approximate wave function approximations of
        one or several levels but may concern also a mean-field basis (for some multiplet of some given configurations) or Green functions,
        etc.

    + name             ::String                      ... to assign a name to the given model.
    + nuclearModel     ::Nuclear.Model               ... Model, charge and parameters of the nucleus.
    + grid             ::Radial.Grid                 ... The radial grid to be used for the computation.
    + refConfigs       ::Array{Configuration,1}      ... List of references configurations, at least 1.
    + repType          ::AbstractRepresentationType  ... Specifies the particular representation.
"""
struct  Representation
    name               ::String  
    nuclearModel       ::Nuclear.Model
    grid               ::Radial.Grid
    refConfigs         ::Array{Configuration,1} 
    repType            ::AbstractRepresentationType
end


"""
`AtomicState.Representation()`  ... constructor for an 'empty' instance of the a variable::AtomicState.Representation
"""
function Representation()
    Representation("", Nuclear.Model(1.0), Radial.Grid(), ManyElectron.Configuration[], CiExpansion())
end


"""
`AtomicState.Representation( ... example for the generation of a mean-field basis)`  

        name        = "Oxygen 1s^2 2s^2 2p^4 ground configuration"
        grid        = Radial.Grid(true)
        nuclearM    = Nuclear.Model(8.)
        refConfigs  = [Configuration("[He] 2s^2 2p^4")]
        mfSettings  = MeanFieldSettings()
        Representation(name, nuclearM, grid, refConfigs, MeanFieldBasis(mfSettings) )
    
`AtomicState.Representation( ... example for the computation of a configuration-interaction (CI) expansion)`  

        name        = "Oxygen 1s^2 2s^2 2p^4 ground configuration"
        grid        = Radial.Grid(true)
        nuclearM    = Nuclear.Model(8.)
        refConfigs  = [Configuration("[He] 2s^2 2p^4")]
        orbitals    = wb["mean-field basis"].orbitals #   get a proper set of orbitals
        ciSettings  = CiSettings(true, false, Int64[], false, LevelSymmetry[] )
        from        = [Shell("2s")]
        to          = [Shell("2s"), Shell("2p")]
        excitations = RasStep(RasStep(), seFrom=from, seTo=to, deFrom=from, deTo=to, frozen=[Shell("1s")])
        Representation(name, nuclearM, grid, refConfigs, CiExpansion(orbitals, excitations, ciSettings) )
    
`AtomicState.Representation( ... example for the computation of a restricted-active-space (RAS) expansion)`  

        name        = "Beryllium 1s^2 2s^2 ^1S_0 ground state"
        refConfigs  = [Configuration("[He] 2s^2")]
        rasSettings = RasSettings([1], 24, 1.0e-6, CoulombInteraction(), true, [1,2,3] )
        from        = [Shell("2s")]
        
        frozen      = [Shell("1s")]
        to          = [Shell("2s"), Shell("2p")]
        step1       = RasStep(RasStep(), seFrom=from, seTo=deepcopy(to), deFrom=from, deTo=deepcopy(to), 
                                frozen=deepcopy(frozen))

        append!(frozen, [Shell("2s"), Shell("2p")])
        append!(to,     [Shell("3s"), Shell("3p"), Shell("3d")])
        step2       = RasStep(step1; seTo=deepcopy(to), deTo=deepcopy(to), frozen=deepcopy(frozen))

        append!(frozen, [Shell("3s"), Shell("3p"), Shell("3d")])
        append!(to,     [Shell("4s"), Shell("4p"), Shell("4d"), Shell("4f")])
        step3       = RasStep(step2, seTo=deepcopy(to), deTo=deepcopy(to), frozen=deepcopy(frozen))

        Representation(name, Nuclear.Model(4.), Radial.Grid(true), refConfigs, 
                        RasExpansion(LevelSymmetry(0, Basics.plus), 4, [step1, step2, step3], rasSettings) )
    
`AtomicState.Representation( ... example for the computation of Green(function) expansion)`  

        name            = "Lithium 1s^2 2s ground configuration"
        refConfigs      = [Configuration("[He] 2s")]
        levelSymmetries = [LevelSymmetry(1//2, Basics.plus), LevelSymmetry(3//2, Basics.plus)]
        greenSettings   = GreenSettings(5, [0, 1, 2], 0.01, true, false, Int64[])
        Representation(name, Nuclear.Model(8.), Radial.Grid(true), refConfigs, 
                        GreenExpansion( AtomicState.DampedSpaceCI(), Basics.DeExciteSingleElectron(), 
                        levelSymmetries, 3, greenSettings) ) 
                        
    ... These simple examples can be further improved by overwriting the corresponding parameters.
"""
function Representation(wa::Bool)    
    AtomicState.Representation()    
end


# `Base.string(rep::Representation)`  ... provides a String notation for the variable rep::AtomicState.Representation
function Base.string(rep::Representation)
    sa = "Atomic representation:   $(rep.name) for Z = $(rep.nuclearModel.Z) and with reference configurations: \n   "
    for  refConfig  in  rep.refConfigs     sa = sa * string(refConfig) * ",  "     end
    return( sa )
end


# `Base.show(io::IO, rep::Representation)`  ... prepares a printout of rep::Representation.
function Base.show(io::IO, rep::Representation)
    sa = Base.string(rep);            print(io, sa, "\n")
    println(io, "representation type:   $(rep.repType)  ")
    println(io, "nuclearModel:          $(rep.nuclearModel)  ")
    println(io, "grid:                  $(rep.grid)  ")
end

end # module
