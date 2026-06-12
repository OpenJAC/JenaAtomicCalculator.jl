
"""
`module  JAC.ManyElectron`  
... a submodel of JAC that contains all structs and methods for many-electron CSF, basis and levels.
"""
module ManyElectron


using   ..Basics, ..Defaults,  ..Radial
export  AsfSettings, CsfR, Basis, Level, Multiplet, 
        AbstractQedModel, NoneQed, QedPetersburg, QedSydney, LSjjSettings



"""
`abstract type ManyElectron.AbstractQedModel` 
    ... defines an abstract and a number of singleton types for dealing with QED estimates in the many-electron
        computations.

    + struct QedPetersburg         
    ... to estimate the QED corrections due to the model Hamiltonian by Shabaev and coworkers (2013).
        
    + struct QedSydney                
    ... to estimate the QED corrections by means of the radiative potential by Flambaum and Ginges (2004).
        
    + struct NoneQed               
    ... No QED estimates are included into the computations.
"""
abstract type  AbstractQedModel                      end
struct         QedPetersburg  <:  AbstractQedModel   end
struct         QedSydney      <:  AbstractQedModel   end
struct         NoneQed        <:  AbstractQedModel   end


"""
`abstract type ManyElectron.AbstractStartOrbitals` 
    ... defines an abstract and a number of singleton types to determine the orbitals from which the SCF computations start.

    + struct StartFromHydrogenic   ... to start the SCF procedure from hydrogenic orbitals.        
    + struct StartFromPrevious     ... to start the SCF procedure from previously calculated orbitals.        
"""
abstract type  AbstractStartOrbitals                          end
struct     StartFromHydrogenic    <:  AbstractStartOrbitals   end

export  AbstractStartOrbitals, StartFromHydrogenic, StartFromPrevious


"""
`struct  ManyElectron.StartFromPrevious  <: ManyElectron.AbstractStartOrbitals`   ... to start from a set of previously calculated orbitals.

    + orbitals      ::Dict{Subshell, Orbital}  ... set of previously-calculated orbitals. 
"""
struct   StartFromPrevious  <: ManyElectron.AbstractStartOrbitals
    orbitals        ::Dict{Subshell, Orbital}
end


function Base.string(from::ManyElectron.StartFromPrevious)
    sa = " "
    for  (sh,v) in from.orbitals  sa = sa * "  " * string(sh)   end
    return( sa )
end

function Base.show(io::IO, from::ManyElectron.StartFromPrevious)
    sa = "Start from previously-calculated orbitals for shells:\n" * string(from);       print(io, sa, "\n")
end


"""
`struct  ManyElectron.LSjjSettings`  ... defines a type for the details and parameters for performing jj-LS expansions.

    + makeIt         ::Bool            ... True, if the jj-LS expansion is to be made and false otherwise.
    + minWeight      ::Float64         ... minimum weight with which a (relativistic) CSF must contribute 
                                            to (at least one) level.
    + printWeight    ::Float64         ... minimum weight of a nonrelativistic CSF to be printed out in the 
                                            final expansion.
    + levelSelection ::LevelSelection  ... Specifies the selected levels, if any.
"""
struct LSjjSettings 
    makeIt          ::Bool
    minWeight       ::Float64 
    printWeight     ::Float64  
    levelSelection  ::LevelSelection
end 


"""
`LSjjSettings(makeIt::Bool)`  ... constructor for the default values of jj-LS transformations.
"""
function LSjjSettings(makeIt::Bool)
    LSjjSettings(makeIt, 0.05, 0.1, LevelSelection() )
end


# `Base.show(io::IO, settings::LSjjSettings)`  ... prepares a proper printout of the variable settings::LSjjSettings.
function Base.show(io::IO, settings::LSjjSettings) 
    println(io, "makeIt:            $(settings.makeIt)  ")
    println(io, "minWeight:         $(settings.minWeight)  ")
    println(io, "printWeight:       $(settings.printWeight)  ")
    println(io, "levelSelection:    $(settings.levelSelection)  ")
end


"""
`struct  ManyElectron.AsfSettings`  
    ... a struct for defining the settings for the atomic state functions, i.e. the self-consistent-field (SCF) 
        and CI computations

    + generateScf          ::Bool                   ... True, if a SCF need to be generated, and false otherwise 
                                                        (frozen orbitals).
    + eeInteraction        ::AbstractEeInteraction  ... Specify the e-e interaction to be included into the SCF 
                                                        computations.
    + scField              ::AbstractScField        ... Specify the self-consistent field, for instance, 
                                                        Basics.ALField(), etc.
    + startScfFrom         ::AbstractStartOrbitals  ... Specify the orbitals to start the SCF computations
    + maxIterationsScf     ::Int64                  ... maximum number of SCF iterations
    + accuracyScf          ::Float64                ... convergence criterion for the SCF field.
    + shellSequenceScf     ::Array{Subshell,1}      ... Sequence of subshells to be optimized.
    + frozenSubshells      ::Array{Subshell,1}      ... Sequence of subshells to be kept frozen.
    
    + eeInteractionCI      ::AbstractEeInteraction  ... Specify the e-e interaction to be included into the 
                                                        CI computations.
    + qedModel             ::AbstractQedModel       ... model for estimating QED corrections {NoneQed(), 
                                                        QedPetersburg(), QedSydney()}.
    + jjLS                 ::LSjjSettings           ... settings to control a jj-LS transformation of atomic 
                                                        level, if requested.
    + levelSelectionCI     ::LevelSelection         ... Specifies the selected levels, if any.
"""
struct  AsfSettings
    generateScf            ::Bool 
    eeInteraction          ::AbstractEeInteraction 
    scField                ::AbstractScField 
    startScfFrom           ::AbstractStartOrbitals
    maxIterationsScf       ::Int64 
    accuracyScf            ::Float64   
    shellSequenceScf       ::Array{Subshell,1}
    frozenSubshells        ::Array{Subshell,1}
    #
    eeInteractionCI        ::AbstractEeInteraction
    qedModel               ::AbstractQedModel 	
    jjLS                   ::LSjjSettings 
    levelSelectionCI       ::LevelSelection
    end


"""
`ManyElectron.AsfSettings()`  ... constructor for setting the default values.
"""
function AsfSettings()
    AsfSettings(true, CoulombInteraction(), Basics.DFSField(), StartFromHydrogenic(), 24, 1.0e-6, Subshell[], Subshell[],  
                CoulombInteraction(), NoneQed(), LSjjSettings(false), LevelSelection() )
end


"""
`ManyElectron.AsfSettings(settings::AsfSettings;`
    
            generateScf=..,       eeInteraction=..,       scField=..,            startScfFrom=..,           
            maxIterationsScf=..,  accuracyScf=..,         shellSequenceScf=..,   frozenSubshells=..,    
            eeInteractionCI=..,   qedModel=..,            jjLS=..,               levelSelectionCI=..,     
            printout::Bool=false)
    ... constructor for re-defining a settings::AsfSettings.
"""
function AsfSettings(settings::AsfSettings; 
    generateScf::Union{Nothing,Bool}=nothing,                       eeInteraction::Union{Nothing,AbstractEeInteraction}=nothing,         
    scField::Union{Nothing,AbstractScField}=nothing,                startScfFrom::Union{Nothing,AbstractStartOrbitals}=nothing,
    maxIterationsScf::Union{Nothing,Int64}=nothing,                 accuracyScf::Union{Nothing,Float64}=nothing,     
    shellSequenceScf::Union{Nothing,Array{Subshell,1}}=nothing,     frozenSubshells::Union{Nothing,Array{Subshell,1}}=nothing, 
    eeInteractionCI::Union{Nothing,AbstractEeInteraction}=nothing,  qedModel::Union{Nothing,AbstractQedModel}=nothing,              
    jjLS::Union{Nothing,LSjjSettings}=nothing,  
    levelSelectionCI::Union{Nothing,LevelSelection}=nothing,        printout::Bool=false)

    if  isnothing(generateScf)           generateScfx          = settings.generateScf           else   generateScfx          = generateScf          end 
    if  isnothing(eeInteraction)         eeInteractionx        = settings.eeInteraction         else   eeInteractionx        = eeInteraction        end 
    if  isnothing(scField)               scFieldx              = settings.scField               else   scFieldx              = scField              end 
    if  isnothing(startScfFrom)          startScfFromx         = settings.startScfFrom          else   startScfFromx         = startScfFrom         end 
    if  isnothing(maxIterationsScf)      maxIterationsScfx     = settings.maxIterationsScf      else   maxIterationsScfx     = maxIterationsScf     end 
    if  isnothing(accuracyScf)           accuracyScfx          = settings.accuracyScf           else   accuracyScfx          = accuracyScf          end 
    if  isnothing(shellSequenceScf)      shellSequenceScfx     = settings.shellSequenceScf      else   shellSequenceScfx     = shellSequenceScf     end 
    if  isnothing(frozenSubshells)       frozenSubshellsx      = settings.frozenSubshells       else   frozenSubshellsx      = frozenSubshells      end 
    if  isnothing(eeInteractionCI)       eeInteractionCIx      = settings.eeInteractionCI       else   eeInteractionCIx      = eeInteractionCI      end 
    if  isnothing(qedModel)              qedModelx             = settings.qedModel              else   qedModelx             = qedModel             end 
    if  isnothing(jjLS)                  jjLSx                 = settings.jjLS                  else   jjLSx                 = jjLS                 end 
    if  isnothing(levelSelectionCI)      levelSelectionCIx     = settings.levelSelectionCI      else   levelSelectionCIx     = levelSelectionCI     end 
    
    AsfSettings(generateScfx, eeInteractionx, scFieldx, startScfFromx, maxIterationsScfx, accuracyScfx, 
                shellSequenceScfx, frozenSubshellsx, eeInteractionCIx, qedModelx, jjLSx, levelSelectionCIx)
end


# `Base.show(io::IO, settings::AsfSettings)`  ... prepares a proper printout of the settings::AsfSettings.
function Base.show(io::IO, settings::AsfSettings)
        println(io, "generateScf:          $(settings.generateScf)  ")
        println(io, "eeInteraction:        $(settings.eeInteraction)  ")
        println(io, "scField:              $(settings.scField)  ")
        println(io, "startScfFrom:         $(settings.startScfFrom)  ")
        println(io, "maxIterationsScf:     $(settings.maxIterationsScf)  ")
        println(io, "accuracyScf:          $(settings.accuracyScf)  ")
        println(io, "shellSequenceScf:     $(settings.shellSequenceScf)  ")
        println(io, "frozenSubshells:      $(settings.frozenSubshells)  ")
        #
        println(io, "eeInteractionCI:      $(settings.eeInteractionCI)  ")
        println(io, "qedModel :            $(settings.qedModel)  ")
        println(io, "jjLS:                 $(settings.jjLS.makeIt)  ")
        println(io, "levelSelectionCI:     $(settings.levelSelectionCI)  ")
end


# `Base.string(settings::AsfSettings)`  ... provides a String notation for the variable settings::AsfSettings.
function Base.string(settings::AsfSettings)
        error("Not yet implemented.")
        sa = "Asf settings: maximum No. of iterations = $(settings.maxIterationsScf), accuracy = (settings.accuracyScf)"
        return( sa )
end


"""
`Base.:(==)(seta::AsfSettings, setb::AsfSettings)`  
    ... compares (recursively) two instances of AsfSettings and return true if all subfields are equal, 
        and false otherwise
"""
function  Base.:(==)(seta::AsfSettings, setb::AsfSettings)

    if  seta.generateScf      !=  setb.generateScf                return( false )    end
    if  seta.eeInteraction    !=  setb.eeInteraction              return( false )    end
    if  seta.scField          !=  setb.scField                    return( false )    end
    if  seta.startScfFrom     !=  setb.startScfFrom               return( false )    end
    #   startOrbitals
    if  seta.maxIterationsScf !=  setb.maxIterationsScf           return( false )    end
    if  seta.accuracyScf      !=  setb.accuracyScf                return( false )    end
    if  seta.shellSequenceScf !=  setb.shellSequenceScf           return( false )    end
    if  seta.frozenSubshells  !=  setb.frozenSubshells            return( false )    end
    if  seta.eeInteractionCI  !=  setb.eeInteractionCI            return( false )    end
    #   qedModel
    #   jjLS
    
        return( true )
end



"""
`struct  ManyElectron.CsfR`  
    ... defines a type for a relativistic configuration state function (CSF) in terms of its orbitals (sequence of orbitals), 
        their occupation as well as the angular momenta, seniority and the coupling of the corresponding antisymmetric subshell states.

    + useStandardSubshells ::Bool                 ... Determines whether the subshell list has to be obtained from some outer structure
                                                        given by the subshells field below. 
    + J                    ::AngularJ64           ... Total angular momentum J.
    + parity               ::Parity               ... Total parity.
    + occupation           ::Array{Int64,1}       ... occupation of the orbitals with regard to the specified orbital list.
    + seniorityNr          ::Array{Int64,1}       ... list of seniority or Nr values of the antisymmetric subshell states
    + subshellJ            ::Array{AngularJ64,1}  ... list of J-values of the antisymmetric subshell states.
    + subshellX            ::Array{AngularJ64,1}  ... intermediate X-values in the coupling of the antisymmetric subshell states.
    + subshells            ::Array{Subshell,1}    ... Explicitly given subshell list if useStandardSubshells == false.
"""
struct  CsfR
    useStandardSubshells   ::Bool 
    J		               ::AngularJ64
    parity  	           ::Parity
    occupation	           ::Array{Int64,1}
    seniorityNr	           ::Array{Int64,1}
    subshellJ	           ::Array{AngularJ64,1}
    subshellX	           ::Array{AngularJ64,1}
    subshells	           ::Array{Subshell,1}
end


"""
`ManyElectron.CsfR(J::AngularJ64, parity::Parity)`  
    ... simple constructor for given J and parity, and where standardOrbitals is set to false.
"""
function CsfR(J::AngularJ64, parity::Parity)
        CsfR(false, J, parity, Int64[], Int64[], AngularJ64[], AngularJ64[], Subshell[] )
end


"""
`ManyElectron.CsfR(sa::String)`  
    ... simple constructor for given closed-shell occupation (such as "[Ne]"), and where standardOrbitals is set to true.
"""
function CsfR(sa::String)
        wa = Basics.extractFromConfiguration(Basics.ClosedSubshells(), sa)

        occupation = Int64[];	  subshellJ = AngularJ64[];	subshellX  = AngularJ64[]
        subshells  = Subshell[];  seniorityNr = Int64[]
    
    for  sb  in  wa
            #x wb = Basics.SubshellQuantumNumbers(string(sb));   kappa = wb[2];   occ = wb[4] + 1
            occ = Basics.subshell_2j(sb) + 1
            occupation   = push!(occupation, occ)
            seniorityNr  = push!(seniorityNr,    0) 
            subshellJ    = push!(subshellJ, AngularJ64(0) ) 
            subshellX    = push!(subshellX, AngularJ64(0) ) 
            subshells    = push!(subshells, sb ) 
        end
        J = AngularJ64(0);    parity = plus

        CsfR(false, J, parity, occupation, seniorityNr, subshellJ, subshellX, subshells)
end


# `Base.show(io::IO, csf::CsfR)`  ... prepares a proper printout of the variable orbital::CsfR.
function Base.show(io::IO, csf::CsfR) 
        print(io, string(csf) ) 
end


# `Base.string(csf::CsfR)`  ... provides a String notation for csf::CsfR.
function Base.string(csf::CsfR) 
    sa = "\n   CSF: "
    if    csf.useStandardSubshells
        for  i in 1:length(csf.occupation)
            sa = sa * Basics.subshellStateString(string(Defaults.GBL_STANDARD_SUBSHELL_LIST[i]), csf.occupation[i], csf.seniorityNr[i], 
                                                    csf.subshellJ[i], csf.subshellX[i])
            sa = sa * ", "
        end
        sa = sa * ": J=" * string(csf.J) * string(csf.parity)
    else
        for  i in 1:length(csf.subshells)
            sa = sa * Basics.subshellStateString(string(csf.subshells[i]), csf.occupation[i], csf.seniorityNr[i], 
                                                    csf.subshellJ[i], csf.subshellX[i])
            sa = sa * ", "
        end
        sa = sa * ": J=" * string(csf.J) * string(csf.parity)
    end
    return( sa )
end


"""
`Base.:(==)(csfa::CsfR, csfb::CsfR)`  
    ... compares (recursively) two relativistic CSFs and return true if all subfields are equal, 
        and false otherwise
"""
function  Base.:(==)(csfa::CsfR, csfb::CsfR)
    # Stop the comparison with an error message if the CSF are not defined w.r.t a common subshell list.
    if   !(csfa.useStandardSubshells)   ||   !(csfb.useStandardSubshells)	## ||   length(csfa.occupation) != length(csfb.occupation)
        error("To determine the equivalence of two CSF requires that both are defined on a common subshell list.")
    end

    if  length(csfa.occupation) != length(csfb.occupation)      return( false )    end
    if  csfa.J  !=  csfb.J  ||  csfa.parity  !=  csfb.parity	return( false )    end
    if  csfa.occupation      !=  csfb.occupation			    return( false )    end
    if  csfa.seniorityNr       !=  csfb.seniorityNr			        return( false )    end
    if  csfa.subshellJ       !=  csfb.subshellJ			        return( false )    end
    if  csfa.subshellX       !=  csfb.subshellX			        return( false )    end
    
        return( true )
end




"""
`ManyElectron.CsfRGrasp92(subshells::Array{Subshell,1}, coreSubshells::Array{Subshell,1}, sa::String, sb::String, sc::String)`  
    ... to construct a CsfR from 3 Grasp-Strings and the given list of subshells and core subshells.
"""
function CsfRGrasp92(subshells::Array{Subshell,1}, coreSubshells::Array{Subshell,1}, sa::String, sb::String, sc::String)
    occupation = Int64[];	 seniorityNr = Int64[];	 subshellJ = AngularJ64[];    subshellX = AngularJ64[];    subshellsx = Subshell[]

    # Define the occupation and quantum numbers of the closed shells
    for i in 1:length(coreSubshells)
        push!(occupation, Basics.subshell_2j( coreSubshells[i] ) +1);    
        push!(seniorityNr, 0);     push!(subshellJ, AngularJ64(0));	 push!(subshellX, AngularJ64(0)) 
    end

    scc = sc * "	  ";   scx = " "
    i  = -1;   ishell = length(coreSubshells)
    while true
        i = i + 1;    
        sax = sa[9i+1:9i+9];   
        #println("sax : ", sax)  
        if length(sb) < 9i+9
            sbx = ""
        else
            sbx = sb[9i+1:9i+9]
        end 
        #println("sbx : ", sbx)  
        scx = strip( scc[9i+6:9i+14] )
        #println("scx : ", scx)  

        sh = strip( sax[1:5] );    sh = Basics.subshellGrasp(sh);    occ = parse( Int64, sax[7:8] )
        #search(sbx, ',') > 0    &&    error("stop a: missing seniorityNr; sb = $sb")   # Include seniorityNr more properly if it occurs
        #subJx = parse( sbx )
        subJx = parseA( sbx )
        #println(subJx)
        if      typeof(subJx) == Nothing     subJ = AngularJ64(0)
        elseif  typeof(subJx) == Int64    subJ = AngularJ64(subJx)
        elseif  typeof(subJx) == Expr     subJ = AngularJ64( subJx.args[2]// subJx.args[3] )
        else    error("stop b; type = $(typeof(subJx))")
        end

        if   length(scx) > 0  &&  scx[end:end] in ["+", "-"]    scx = scx[1:end-1]   end
        subXx = parseA( scx )
        if      ishell == 0  &&  typeof(subXx)     == Nothing                                       subX = subJ
        elseif  ishell >= 1  &&  subshellX[ishell] == AngularJ64(0)  && typeof(subXx) == Nothing    subX = subJ
        elseif  ishell >= 1  &&  typeof(subXx)     == Nothing                                       subX = subshellX[ishell]
        elseif  typeof(subXx) == Nothing                     subX = subshellX[ishell-1]
        elseif  typeof(subXx) == Int64                    subX = AngularJ64(subXx)
        elseif  typeof(subXx) == Expr                     subX = AngularJ64( subXx.args[2]// subXx.args[3] )
        else    error("stop c")
        end

        # Fill the arrays due to the pre-defined sequence of subshells
        while  true
        ishell = ishell + 1
        if  subshells[ishell] == sh
            wa = ManyElectron.provideSubshellStates(sh, occ);    nu = -1
            for  a in wa
            if   AngularJ64( a.Jsub2//2 ) == subJ	 nu = a.nu;    break   end
            end
            nu < 0    &&    error("stop d")
            push!(occupation, occ);    push!(seniorityNr, nu);	 push!(subshellJ, subJ);    push!(subshellX, subX) 
            break
        else
            push!(occupation, 0);      
            push!(seniorityNr, 0);	 
            push!(subshellJ, AngularJ64(0));    
            push!(subshellX, subshellX[ishell-1]) 
            #println(subshellX, ishell)
            #push!(subshellX, AngularJ64(0)) 
        end
        end

        ishell > length(subshells)    &&	error("stop e")

        # Terminate if all subshells are read
        if  9i + 10 >  length(sa)	 break    end
    end

    while  ishell < length(subshells)
        ishell = ishell + 1
        push!(occupation, 0);    
        push!(seniorityNr, 0);	 
        push!(subshellJ, AngularJ64(0));    
        push!(subshellX, subshellX[ishell-1]) 
    end

    # Fill the remaining subshells

    J = subshellX[end];    scx = strip(sc)
    if	  scx[end:end] == "+"  parity = Parity("+")
    elseif  scx[end:end] == "-"  parity = Parity("-")
    else	  error("stop f")
    end

    wa = CsfR(true, J, parity, occupation, seniorityNr, subshellJ, subshellX, subshellsx)
end


"""
`function parseA(x::AbstractString)`
    ... helper function for `ManyElectron.CsfRGrasp92`
"""
function parseA(x::AbstractString)
    y = strip(x)
    if length(y) == 0 
        return nothing
    elseif length(y) <= 2
        return parse(Int64, y)
    elseif length(y) == 3
        a = parse(Int64, y[1]); b = parse(Int64, y[3])
        return :( $a//$b )
    end
end


"""
`ManyElectron.CsfRTransformToStandardSubshells(csf::CsfR, subshells::Array{Subshell,1})`  
    ... constructor to re-define a CSF with regard to the given list of subshells. It incorporates 'empty' subshell state 
        and their coupling (if necessary) but does not support a re-arrangement of the subshell order. This constructor 
        is used to 'combine' individual CsfR into a common list (such as in a Basis) and, thus, use StandardSubshells = true 
        is set. This method terminates with an error message if the given subshell list requires a re-coupling of the CSF.
"""
function CsfRTransformToStandardSubshells(csf::CsfR, subshells::Array{Subshell,1})
    J = csf.J;    parity = csf.parity; 
    occupation = Int64[];	 seniorityNr = Int64[];	  subshellJ = AngularJ64[];	subshellX = AngularJ64[]
    
    ish = 0
    for  subsh in subshells
        push!(subshells, subsh)
        found = false
        for  i in 1:length(csf.subshells)
        if     subsh == csf.subshells[i]  &&   i >   ish    i = ish;  found = true;   break    
        elseif subsh == csf.subshells[i]  &&   i <=  ish    error("CSF requires recoupling; this is not supported.")
        end
        end
        if  found
        push!(occupation, csf.occupation[ish]);    push!(seniorityNr, csf.seniorityNr[ish])   
        push!(subshellJ,  csf.subshellJ[ish]);     push!(subshellX, csf.subshellX[ish])   
        else
        push!(occupation, 0);	 push!(seniorityNr, 0);	 push!(subshellJ,  0)
        if    ish == 0    push!(subshellX, AngularJ64(0) )   
        else		  push!(subshellX, subshellX[end])  
        end 
        end
    end

    CsfR(true, J, parity, occupation, seniorityNr, subshellJ, subshellX, Subshell[] )
end


"""
`struct  ManyElectron.Basis`  
    ... defines a type for a relativistic atomic basis, including the (full) specification of the configuration space 
        and radial orbitals.

    + isDefined      ::Bool                      ... Determines whether this atomic basis 'points' to some a well-defined basis
                                                        or to an empty instance, otherwise.
    + NoElectrons    ::Int64                     ... No. of electrons.
    + subshells      ::Array{Subshell,1}         ... Explicitly given subshell list for this basis.
    + csfs           ::Array{CsfR,1}             ... List of CSF.	    
    + coreSubshells  ::Array{Subshell,1}         ... List of (one-electron) core orbitals.
    + orbitals       ::Dict{Subshell, Orbital}   ... Dictionary of (one-electron) orbitals.
"""
struct  Basis
    isDefined	       ::Bool
        NoElectrons	 ::Int64 
    subshells	       ::Array{Subshell,1}
    csfs             ::Array{CsfR,1}	      
        coreSubshells	 ::Array{Subshell,1}
        orbitals	       ::Dict{Subshell, Orbital}
end 


"""
`ManyElectron.Basis()`  ... constructor for an 'empty' instance::Basis with isDefined == false.
"""
function Basis()
        Basis(false, 0, Subshell[], CsfR[], Subshell[], Dict{Subshell, Orbital}() )
end


"""
`ManyElectron.Basis(basis::Basis, subshells::Array{Subshell,1})`  
    ... constructor to re-expand a given basis with regard to an extended subshell list; this expansion does not change
        the CSF basis but makes sure that the occupation and other quantum numbers are properly expanded. 
        An error message is issued if the sequence of the subshells does not allow such a simple re-expansion.
        A newBasis::Basis is returned.
"""
function Basis(basis::Basis, subshells::Array{Subshell,1})
    if basis.subshells == subshells   return( basis )   end
    # Check that the subshells allow such an re-expansion
    oldSubshells = basis.subshells
    newSubshells = Subshell[]
    for  subsh in subshells
        if subsh in oldSubshells   push!(newSubshells, subsh)   end
    end
    if  oldSubshells != newSubshells         error("Improper subshell order: $oldSubshells != $newSubshells")        end
    if  !basis.csfs[1].useStandardSubshells  error("Re-expansion of basis only if useStandardSubshells == true.")    end
    
    # Define the newCoreSubshells
    nx = length(basis.coreSubshells);   newCoreSubshells = Subshell[];   ny = 0
    for  ny = nx:-1:1
        if basis.coreSubshells[1:ny] == subshells[1:ny]     newCoreSubshells = copy(subshells[1:ny]);   break   end
    end
    ## x@show basis.coreSubshells, newCoreSubshells
    # Define the new CSF basis
    newCsfs = CsfR[]
    for  csf in basis.csfs
        noccupation = Int64[];   nseniorityNr = Int64[];   nsubshellJ = AngularJ64[];   nsubshellX = AngularJ64[]
        for  subsh  in  subshells
            if  subsh in oldSubshells  
                # Determine the index in oldSubshells
                ii = 0;     for i = 1:length(oldSubshells)   if  subsh == oldSubshells[i]   ii = i; break   end     end
                push!(noccupation,  csf.occupation[ii])
                push!(nseniorityNr, csf.seniorityNr[ii])
                push!(nsubshellJ,   csf.subshellJ[ii]) 
                push!(nsubshellX,   csf.subshellX[ii])
            elseif   length(nsubshellX) == 0
                # Extend quantum numbers properly
                push!(noccupation,  0)
                push!(nseniorityNr, 0)
                push!(nsubshellJ,   AngularJ64(0)) 
                push!(nsubshellX,   AngularJ64(0))
            else
                # Extend quantum numbers properly
                push!(noccupation,  0)
                push!(nseniorityNr, 0)
                push!(nsubshellJ,   AngularJ64(0)) 
                push!(nsubshellX,   nsubshellX[end])
            end
        end
        newcsf = CsfR(csf.useStandardSubshells, csf.J, csf.parity, noccupation, nseniorityNr, nsubshellJ, nsubshellX, csf.subshells)
        push!(newCsfs, newcsf)
    end
    
    Basis(basis.isDefined, basis.NoElectrons, subshells, newCsfs, newCoreSubshells, basis.orbitals)
end


"""
`ManyElectron.Basis("from Grasp2013", cslFilename::String, rwfFilename::String)`  
    ... to construct an instance::Basis from the Grasp92/Grasp2013 .csl and .rwf files.
"""
function Basis(sa::String, cslFilename::String, rwfFilename::String)
    if  sa == "from Grasp2013"
        isDefined = true 
        #
        println("ManyElectron.Basis-aa: warning ... The standard grid is set; make an interactive request if not appropriate.") 
        grid     = Radial.Grid(true)
        basis    = Basics.read(ReadCslFileGrasp92(), cslFilename)
        orbitals = Basics.readOrbitalFileGrasp92(rwfFilename, grid)
    else  error("stop a")
    end

        Basis( isDefined, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, orbitals) 
end


"""
`Base.show(io::IO, basis::Basis; full=false)`  ... prepares a proper printout of the variable basis::Basis.
"""
function Base.show(io::IO, basis::Basis) 
    print(io, string(basis), "\n")
        if  basis.isDefined
        println(io, "isDefined:          $(basis.isDefined)  ")
        println(io, "NoElectrons:        $(basis.NoElectrons)  ")
        println(io, "subshells:          $(basis.subshells)  ")
        println(io, "csfs:               $(basis.csfs)  ")
        println(io, "coreSubshells:      $(basis.coreSubshells)  ")
        println(io, "orbitals:           $(basis.orbitals)  ")
        end 
end


# `Base.string(basis::Basis)`  ... provides a String notation for variable basis::Basis.
function Base.string(basis::Basis) 
    if  basis.isDefined
        sa = "Atomic basis:  $(length(basis.csfs)) CSF defined for $(basis.NoElectrons) elecrons."
    else
        sa = "Atomic basis:  Undefined."
    end
    return(sa)
end


"""
`struct  ManyElectron.Level`  
    ... defines a type for an atomic level in terms of its quantum number, energy and with regard to a relativistic basis.

    + J            ::AngularJ64       ... Total angular momentum J.
    + M            ::AngularM64       ... Total projection M, only defined if a particular sublevel is referred to.
    + parity       ::Parity           ... Parity of the level.
    + index        ::Int64            ... Index of this level in its original multiplet, or 0.
    + energy       ::Float64          ... energy
    + relativeOcc  ::Float64          ... Relative occupation of this level, if involved in the evolution of a cascade.
    + hasStateRep  ::Bool             ... Determines whether this level 'point' to a physical representation of an state
                                            (i.e. an basis and corresponding mixing coefficient vector), or just to empty instances, otherwise.
    + basis        ::Basis            ... relativistic atomic basis
    + mc           ::Vector{Float64}  ... Vector of mixing coefficients w.r.t basis.
"""
struct  Level
    J		       ::AngularJ64
    M		       ::AngularM64
    parity  	   ::Parity
    index  	       ::Int64
    energy  	   ::Float64
    relativeOcc    ::Float64
    hasStateRep    ::Bool
    basis	       ::Basis
    mc		       ::Vector{Float64}
end 


"""
`ManyElectron.Level()`  ... constructor an empty level without useful data.
"""
function Level()
    Level( AngularJ64(0),  AngularM64(0), Parity("+"), 0, 0., 0., false, Basis(), Vector{Float64}[] )
end


"""
`ManyElectron.Level(J::AngularJ64, parity::Parity, energy::Float64, relativeOcc::Float64)`  
    ... constructor an atomic level without a state representation: hasStateRep == false.
"""
function Level(J::AngularJ64, parity::Parity, energy::Float64, relativeOcc::Float64)
    Level(J, J, parity, 0, energy, relativeOcc, false, Basis(), Vector{Float64}[] )
end


"""
`ManyElectron.Level(level::Level, subshells::Array{Subshell,1})`  
    ... constructor to re-expand a given level with regard to an extended subshell list; this expansion does not change
        the CSF basis nor the mixing coeffients but makes sure that the occupation and other quantum numbers are properly
        expanded. An error message is issued if the sequence of the subshells does not allow such a simple re-expansion.
        A newLevel::Level is returned.
"""
function Level(level::Level, subshells::Array{Subshell,1})
    basis = Basis(level.basis, subshells)
    Level(level.J, level.M, level.parity, level.index, level.energy, level.relativeOcc, level.hasStateRep, basis, level.mc)
end


# `Base.show(io::IO, level::Level)`  ... prepares a proper printout of the variable level::Level.
function Base.show(io::IO, level::Level) 
    println(io, "Level: J = $(level.J), M = $(level.M), parity = $(level.parity), index = $(level.index) ")
    println(io, "energy:         $(level.energy)  ")
    println(io, "relativeOcc:    $(level.relativeOcc)  ")
    println(io, "hasStateRep:    $(level.hasStateRep)  ")
    println(io, "basis:           (level.basis)  ")
    println(io, "mc:             $(level.mc)  ")
end


"""
`struct  ManyElectron.Multiplet`  
    ... defines a type for an ordered list of atomic levels; has only the default constructor.

    + name  	 ::String	        ... A name associated to the multiplet.
    + levels	 ::Array{Level,1}   ... A list of levels (pointers).
"""
struct  Multiplet
    name         ::String
    levels  	 ::Array{Level,1}
end 


"""
    `ManyElectron.Multiplet()`  ... constructor for providing an 'empty' instance of this struct.
"""
function Multiplet()
    Multiplet("", Level[] )
end


"""
`ManyElectron.Multiplet("from Ratip2012", cslFilename::String, rwfFilename::String, mixFilename::String)`  
    ... to construct an instance::Multiplet from the Grasp92/Ratip2012 .csl, .rwf and .mix files.  Some consistency 
        checks are made and the method terminates with an error message if the files don't fit together. 
    
+ `("from Grasp2018", cslFilename::String, rwfFilename::String, mixFilename::String)` 
    ... to do the same but for files from Grasp2018.
"""
function Multiplet(sa::String, cslFilename::String, rwfFilename::String, mixFilename::String)
    if  sa == "from Ratip2012"
        name      = "from "*cslFilename[1:end-4] 
        basis     = ManyElectron.Basis("from Grasp2013", cslFilename, rwfFilename)
        multiplet = Basics.readMixFileRelci(mixFilename, basis)
    else  error("Unrecognized keystrings:  $sa")
    end

    return( multiplet ) 
end


"""
`ManyElectron.Multiplet(multiplet::Multiplet, subshells::Array{Subshell,1})`  
    ... constructor to re-expand a given Multiplet with regard to an extended subshell list; this expansion does not change
        the CSF basis nor the mixing coeffients but makes sure that the occupation and other quantum numbers are properly
        expanded. An error message is issued if the sequence of the subshells does not allow such a simple re-expansion.
        A newMultiplet::Multiplet is returned.
"""
function Multiplet(multiplet::Multiplet, subshells::Array{Subshell,1})
    basis = Basis(multiplet.levels[1].basis, subshells)
    newLevels = Level[]
    for  level  in  multiplet.levels
        push!( newLevels, Level(level.J, level.M, level.parity, level.index, level.energy, 
                                level.relativeOcc, level.hasStateRep, basis, level.mc) )
    end

    Multiplet(multiplet.name, newLevels) 
end


# `Base.show(io::IO, multiplet::Multiplet)`  ... prepares a proper printout of the variable multiplet::Multiplet.
function Base.show(io::IO, multiplet::Multiplet) 
    println(io, "name:        $(multiplet.name)  ")
    println(io, "levels:      $(multiplet.levels)  ")
end


# Functions/methods that are later added to the module ManyElectron
function matrixElement_Mab                                      end
function matrixElement_Vee                                      end
function provideSubshellStates                                  end

end # module


