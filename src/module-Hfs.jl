
"""
`module  JAC.Hfs`  
    ... a submodel of JAC that contains all methods for computing HFS A, B and C coefficients, hyperfine-level representations, etc.
        We apply IJF-coupling throughout in this module and within JAC. While the hyperfine-interaction usually acts between 
        electronic (CSF) basis states or, more accurately, between IJF-coupled basis vectors, all hyperfine representations
        usually refer to the product basis of nuclear + electronic states/levels in order to keep the number of hyperfine 
        levels moderate.
        
        This (clear) distinctions between a purely IJF-coupled basis and the use of hyperfine levels and multiplets in 
        different applications makes it necessary to implement the following data types:
        
        + IJF_Vector = IJF_Vector(F, parity, isomer::Nuclear.Isomer, csf::CsfR, basisJ::Basis)
          ... such a vector describes a tensor product state nuclear + electronic (state) and is appropriate to form 
              a basis for evaluating hyperfine-interaction amplitudes. The basisJ need to be provided to have full
              access to the definition of the electronic csfR but does not belong to the definition.
                  
        + IJF_Basis  = Vector{IJF_Vector}  ... Such a basis occurs in the evaluation of the hyperfine (Hamiltonian) interaction
              matrix but does not need an own data structure.
        
        + HfBasisVector = HfBasisVector(F, parity, isomer, LevelJ)  ... This data structure enables one to think
              and deal with hyperfine levels as associated with a well-define isomeric state and electronic level J.
              
        + HfBasis = Vector{HfBasisVector}  ... Such a basis occurs in the evaluation of physical interaction 
              matrix elements between hyperfine levels but does not need an own data structure.
              
        + HfLevel = HfLevel(F, parity, hfBasis::Vector{HfBasisVector}, mc::Vector{Float64})
          ... A hyperfine level is the pendant to a (electronic) Level/state, if hyperfine-resolved transitions are considered.
              Each hyperfine level has a representation mc that refers to the hfBasis, and which contains all information
              about the representation of the underlying nuclear and electronic basis states. The electronic basis is formed
              by a selected set of ASF, typically taken from some (electronic) multiplet. In contrast to a pure (electronic)
              CSF basis, the use of ASF simplifies the interpretation of physical findings but cannot reduce the computational 
              effort (perhaps, even slightly increase the computational effort. This would need to analyzed deeper within
              the Julia machinery). Formally, the computational effort should be equivalent.
              
        + HfMultiplet = HfMultiplet(name::String, hfLevels::Vector{HfLevel}) ... this is the pendant to an 
              (electronic) multiplet and enables one to proceed with hyperfine-induced/hyperfine-resolved transitions
              completely analogue as for the electronic transitions themselves. All what is needed to do is to account for 
              the representation (mixing coefficients) of the associated hyperfine levels as well as the matrix elements
              for the nuclear (isomeric) states and those for the associated ASF.
"""
module Hfs


using Printf, ..AngularMomentum, ..Basics,  ..Defaults, ..InteractionStrength, ..ManyElectron, ..Radial, ..Nuclear,
              ..SpinAngular, ..TableStrings, ..PhotoEmission

## Tabulate-theme sentinel types for Basics.tabulate dispatch on Hfs.HfMultiplet.
struct  HfEnergies          end   ## tabulate all hyperfine level energies
struct  HfEnergiesRelative  end   ## tabulate hyperfine level energies relative to the lowest level




#################################################################################################################################
#################################################################################################################################

"""
`struct  Hfs.IJF_Vector`  ... defines a type for a IJF-coupled basis vector, here based on an ASF. Following IJF-coupling,
    an IJF basis vector is always the (tensor) product state of an isomeric state x electronic state. IJF-coupled basis 
    vector are needed and obtained from the diagonalization of the (electronic) HFS Hamiltonian within a given CSF basis.

    + F         ::AngularJ64        ... Total angular momentum F
    + parity    ::Parity            ... Total parity of the basis vector = nuclear x electronic parity.
    + isomer    ::Nuclear.Isomer    ... Isomeric state of the nucleus.
    + csf       ::ManyElectron.CsfR ... Electronic CSF state.
    + basisJ    ::Basis             ... Electronic CSF basis to which the csf refer to.
    
    There is no need to introduce a type IJF_Basis since such a ijfBasis = Hfs.IJF_Vector[...] can be readily formed at all 
    occurrences. 

"""
struct IJF_Vector
    F           ::AngularJ64
    parity      ::Parity
    isomer      ::Nuclear.Isomer
    csf         ::ManyElectron.CsfR 
    basisJ      ::Basis 
end 


"""
`Hfs.IJF_Vector()`  ... constructor for an `empty` instance of IJF_Vector.
"""
function IJF_Vector()
    IJF_Vector(AngularJ64(0), Basics.plus, Nuclear.Isomer(), CsfR("[He]"), Basis())
end


# `Base.show(io::IO, ijfVector::Hfs.IJF_Vector)`  ... prepares a proper printout of the variable ijfVector.
function Base.show(io::IO, ijfVector::Hfs.IJF_Vector) 
    println(io, "F:           $(ijfVector.F)  ")
    println(io, "parity:      $(ijfVector.parity)  ")
    println(io, "isomer:      $(ijfVector.isomer)  ")
    println(io, "csf:         $(ijfVector.csf)  ")
    println(io, "basisJ:      $(ijfVector.basisJ)  ")
end


"""
`struct  Hfs.HfBasisVector`  ... defines a type for a hyperfine basis vector that enables one to think about and deal with 
    individual hyperfine levels. These hyperfine levels have representations with regard to (tensor) product states, which 
    are formed from a set of isomeric states as well as a set of electronic ASF level J. Hyperfine basis vectors are obtained
    from diagonalizing the hyperfine Hamiltonian in a tensor basis of isomeric + ASF states. 

    + F         ::AngularJ64        ... Total angular momentum F
    + parity    ::Parity            ... Total parity of the basis vector = nuclear x electronic parity.
    + isomer    ::Nuclear.Isomer    ... Isomeric state of the nucleus.
    + LevelJ    ::Basis             ... Electronic level that is part of the electronic basis.
    
    There is no need to introduce a type HfBasis since such a hfBasis = Hfs.HfBasisVector[...] can be readily formed at all 
    occurrences. The IJF_Vector's and     HfBasisVector's have different advantages.

"""
struct HfBasisVector
    F           ::AngularJ64
    parity      ::Parity
    isomer      ::Nuclear.Isomer
    levelJ      ::Level
end 
       

"""
`Hfs.HfBasisVector()`  ... constructor for an `empty` instance of HfBasisVector`.
"""
function HfBasisVector()
    HfBasisVector(AngularJ64(0), Basics.plus, Nuclear.Isomer(), Level())
end


# `Base.show(io::IO, hfBasisVector::Hfs.HfBasisVector)`  ... prepares a proper printout of the variable HfBasisVector.
function Base.show(io::IO, hfBasisVector::Hfs.HfBasisVector) 
    println(io, "F:           $(hfBasisVector.F)  ")
    println(io, "parity:      $(hfBasisVector.parity)  ")
    println(io, "isomer:      $(hfBasisVector.isomer)  ")
    println(io, "levelJ:      $(hfBasisVector.levelJ)  ")
end


"""
`struct  Hfs.HfLevel`  ... defines a type for HfLevel with a representation that refers to a product basis of isomeric and 
    ASF states; the HfLevel is the pendant to a (electronic) Level/state, if hyperfine-resolved transitions are considered.
    Each hyperfine level has a representation mc that refers to the hfBasis, and which contains all information about the 
    representation of the underlying nuclear and electronic basis states. The electronic basis is formed by a selected set 
    of ASF, typically taken from some (electronic) multiplet. In contrast to a pure (electronic) IJF_Basis, the use of 
    ASF simplifies the interpretation of physical findings but cannot reduce the computational effort (perhaps, even results
    in slightly increase the computational effort). A HfLevel is defined by:

    + F              ::AngularJ64               ... Total angular momentum F.
    + M              ::AngularM64               ... Total projection M, only defined if a particular magnetic sublevel is referred to.
    + parity         ::Parity                   ... Parity of the level which corresponds to the electronic system.
    + energy         ::Float64                  ... energy
    + hfBbasis       ::Array{HfBasisVector,1}   ... the product basis nuclear (isomeric) state x selected ASF.
    + mc             ::Vector{Float64}          ... Vector of mixing coefficients w.r.t hfBasis.
    
    
"""
struct HfLevel
    F                ::AngularJ64
    M                ::AngularM64
    parity           ::Parity
    energy           ::Float64
    hfBbasis         ::Array{HfBasisVector,1}
    mc               ::Vector{Float64}
end 


"""
`Hfs.HfLevel()`  ... constructor for an `empty` instance of HfLevel.
"""
function HfLevel()
    HfLevel(AngularJ64(0), AngularM64(0), Basics.plus, 0., HfBasisVector[], Float64[])
end


# `Base.show(io::IO, hfLevel::Hfs.HfLevel)`  ... prepares a proper printout of the variable hfLevel::Hfs.HfLevel.
function Base.show(io::IO, hfLevel::Hfs.HfLevel) 
    println(io, "F:              $(hfLevel.F)  ")
    println(io, "M:              $(hfLevel.M)  ")
    println(io, "parity:         $(hfLevel.parity)  ")
    println(io, "energy:         $(hfLevel.energy)  ")
    println(io, "basis:          $(hfLevel.basis)  ")
    println(io, "mc:             $(hfLevel.mc)  ")
end


"""
`struct  Hfs.HfMultiplet`  ... defines a type for a multiplet of hyperfine levels (HfLevel's) which are based on a (tensor)
    product basis of nuclear (isomeric) x ASF states.

    + name     ::String                ... A name associated to the multiplet.
    + hfLevels ::Array{HfLevel,1}      ... List of hyperfine levels (HfLevel's)

"""
struct HfMultiplet
    name       ::String
    hfLevels   ::Array{HfLevel,1}
end 


"""
`Hfs.HfMultiplet()`  ... constructor for an `empty` instance of Hfs.HfMultiplet.
"""
function HfMultiplet()
    HfMultiplet("", HfLevel[])
end

# `Base.show(io::IO, hfMultiplet::Hfs.HfMultiplet)`  ... prepares a proper printout of the variable hfMultiplet::Hfs.HfMultiplet.
function Base.show(io::IO, hfMultiplet::Hfs.HfMultiplet) 
    println(io, "name:           $(hfMultiplet.name)  ")
    println(io, "hfLevels:       $(hfMultiplet.hfLevels)  ")
end


#################################################################################################################################
#################################################################################################################################

"""
`struct  Hfs.InteractionMatrix`  ... defines a type for storing the T^1 and T^2 interaction matrices for a given basis.

    + calcM1   ::Bool               ... true, if the matrixM1 has been calculated and false otherwise.
    + calcE2   ::Bool               ... true, if the matrixE2 has been calculated and false otherwise.
    + calcM3   ::Bool               ... true, if the matrixM3 has been calculated and false otherwise.
    + matrixM1 ::Array{Float64,2}   ... T^M1 interaction matrix
    + matrixE2 ::Array{Float64,2}   ... T^E2 interaction matrix
    + matrixM3 ::Array{Float64,2}   ... T^M3 interaction matrix

"""
struct InteractionMatrix
    calcM1     ::Bool
    calcE2     ::Bool
    calcM3     ::Bool
    matrixM1   ::Array{Float64,2}
    matrixE2   ::Array{Float64,2}
    matrixM3   ::Array{Float64,2}
end 


"""
`Hfs.InteractionMatrix()`  ... constructor for an `empty` instance of InteractionMatrix.
"""
function InteractionMatrix()
    InteractionMatrix(false, false, false, zeros(2,2), zeros(2,2), zeros(2,2))
end


# `Base.show(io::IO, im::Hfs.InteractionMatrix)`  ... prepares a proper printout of the variable InteractionMatrix.
function Base.show(io::IO, im::Hfs.InteractionMatrix) 
    println(io, "calcM1:           $(im.calcM1)  ")
    println(io, "calcE2:           $(im.calcE2)  ")
    println(io, "calcM3:           $(im.calcM3)  ")
    println(io, "matrixM1:         $(im.matrixM1)  ")
    println(io, "matrixE2:         $(im.matrixE2)  ")
    println(io, "matrixM3:         $(im.matrixM3)  ")
end


"""
`struct  Hfs.Outcome`  
    ... defines a type to keep the outcome of a HFS computation, such as the HFS A and B coefficients as well 
        other results.

    + Jlevel                    ::Level            ... Atomic level to which the outcome refers to.
    + gJ                        ::Float64          ... Lande's g_J factor of the level.
    + AIoverMu                  ::Float64          ... HFS A * I / mu value.
    + BoverQ                    ::Float64          ... HFS B / Q value
    + CoverOmega                ::Float64          ... HFS C / Omega value
    + amplitudeM1               ::Complex{Float64} ... M1 amplitude
    + amplitudeE2               ::Complex{Float64} ... E2 amplitude
    + amplitudeM3               ::Complex{Float64} ... M3 amplitude
    + nuclearI                  ::AngularJ64       ... nuclear spin
    + hfsMultiplet              ::HfMultiplet      ... Multiplet of HfLevel's as associated with the JLevel.
"""
struct Outcome 
    Jlevel                      ::Level 
    gJ                          ::Float64 
    AIoverMu                    ::Float64
    BoverQ                      ::Float64
    CoverOmega                  ::Float64
    amplitudeM1                 ::Complex{Float64}
    amplitudeE2                 ::Complex{Float64}
    amplitudeM3                 ::Complex{Float64}
    nuclearI                    ::AngularJ64
    hfsMultiplet                ::HfMultiplet
end 


"""
`Hfs.Outcome()`  ... constructor for an `empty` instance of Hfs.Outcome for the computation of HFS properties.
"""
function Outcome()
    Outcome(Level(), 0., 0., 0., 0., 0., 0., 0., AngularJ64(0), HfMultiplet() )
end


# `Base.show(io::IO, outcome::Hfs.Outcome)`  ... prepares a proper printout of the variable outcome::Hfs.Outcome.
function Base.show(io::IO, outcome::Hfs.Outcome) 
    println(io, "Jlevel:                    $(outcome.Jlevel)  ")
    println(io, "gJ:                        $(outcome.gJ)  ")
    println(io, "AIoverMu:                  $(outcome.AIoverMu)  ")
    println(io, "BoverQ:                    $(outcome.BoverQ)  ")
    println(io, "CoverOmega:                $(outcome.CoverOmega)  ")
    println(io, "amplitudeM1:               $(outcome.amplitudeM1)  ")
    println(io, "amplitudeE2:               $(outcome.amplitudeE2)  ")
    println(io, "amplitudeM3:               $(outcome.amplitudeM3)  ")
    println(io, "nuclearI:                  $(outcome.nuclearI)  ")
    println(io, "hfMultiplet:                (outcome.hfMultiplet)  ")
end


"""
`struct  Settings  <:  AbstractPropertySettings`  ... defines a type for the details and parameters of computing HFS A and B coefficients.

    + calcM1                    ::Bool             ... True if T^M1-amplitudes (HFS A values) need to be calculated, and false otherwise.
    + calcE2                    ::Bool             ... True if T^E2-amplitudes (HFS B values) need to be calculated, and false otherwise.
    + calcM3                    ::Bool             ... True if T^M3-amplitudes (HFS C values) need to be calculated, and false otherwise.
    + calcNondiagonal           ::Bool             
        ... True if also (non-)diagonal hyperfine amplitudes are to be calculated and printed, and false otherwise.
    + calcHfMultiplet           ::Bool             
        ... True if representation of all hyperfine levels need to be generated in the IJF-coupled ASF basis;
            only the nuclear spin and the selected atomic levels are considered in the basis. false otherwise.
    + printBefore               ::Bool             ... True if a list of selected levels is printed before the actual computations start. 
    + levelSelection            ::LevelSelection   ... Specifies the selected levels, if any.
"""
struct Settings  <:  AbstractPropertySettings
    calcM1                      ::Bool
    calcE2                      ::Bool
    calcM3                      ::Bool
    calcNondiagonal             ::Bool 
    calcHfMultiplet             ::Bool 
    printBefore                 ::Bool 
    levelSelection              ::LevelSelection
end 


"""
`Hfs.Settings(; calcM1::Bool=true,` calcE2::Bool=false, calcM3::Bool=false, calcNondiagonal::Bool=false, 
                calcHfMultiplet::Bool=false, printBefore::Bool=false, levelSelection::LevelSelection=LevelSelection()) 
    ... keyword constructor to overwrite selected value of Hfs.Settings. 
"""
function Settings(; calcM1::Bool=true, calcE2::Bool=false, calcM3::Bool=false, calcNondiagonal::Bool=false, 
                    calcHfMultiplet::Bool=false, printBefore::Bool=false, levelSelection::LevelSelection=LevelSelection())
    Settings(calcM1, calcE2, calcM3, calcNondiagonal, calcHfMultiplet, printBefore, levelSelection)
end


# `Base.show(io::IO, settings::Hfs.Settings)`  ... prepares a proper printout of the variable settings::Hfs.Settings.
function Base.show(io::IO, settings::Hfs.Settings) 
    println(io, "calcM1:                   $(settings.calcM1)  ")
    println(io, "calcE2:                   $(settings.calcE2)  ")
    println(io, "calcM3:                   $(settings.calcM3)  ")
    println(io, "calcNondiagonal:          $(settings.calcNondiagonal)  ")
    println(io, "calcHfMultiplet:         $(settings.calcHfMultiplet)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "levelSelection:           $(settings.levelSelection)  ")
end

#######################################################################################################################
#######################################################################################################################


function Base.isless(x::Hfs.HfLevel, y::Hfs.HfLevel)
    return x.energy < y.energy
end


#######################################################################################################################
#######################################################################################################################

    
"""
`Basics.sortByEnergy(multiplet::Hfs.HfMultiplet)`  
    ... to sort all hyperfine levels in the multiplet into a sequence of increasing energy; a (new) multiplet::Hfs.HfMultiplet 
        is returned.
"""
function Basics.sortByEnergy(multiplet::Hfs.HfMultiplet)
    sortedLevels = Base.sort( multiplet.levelFs , lt=Base.isless)
    newLevels = Hfs.HfLevel[];   index = 0
    for lev in sortedLevels
        index = index + 1
        push!(newLevels, Hfs.IJF_Level(lev.I, lev.F, lev.M, lev.parity, lev.energy, lev.basis, lev.mc) )
    end
    
    newMultiplet = Hfs.HfMultiplet(multiplet.name, newLevels)
    
    return( newMultiplet )  
end


"""
`Basics.tabulate(stream::IO, ::HfEnergies, multiplet::Hfs.HfMultiplet)`
    ... tabulates the energies of all hyperfine levels of the given multiplet; nothing is returned.
"""
function Basics.tabulate(stream::IO, ::HfEnergies, multiplet::Hfs.HfMultiplet)
    println(stream, "\n  Eigenenergies for nuclear spin I = $(multiplet.levelFs[1].I):")
    sb = "  Level  F Parity          Hartrees       " * "             eV                   " * TableStrings.inUnits("energy")
    println(stream, "\n", sb, "\n")
    for  i = 1:length(multiplet.levelFs)
        lev          = multiplet.levelFs[i]
        en           = lev.energy
        en_eV        = Defaults.convertUnits("energy: from atomic to eV", en)
        en_requested = Defaults.convertUnits("energy: from atomic",       en)
        sc  = " " * TableStrings.level(i) * "    " * string(LevelSymmetry(lev.F, lev.parity)) * "    "
        @printf(stream, "%s %.15e %s %.15e %s %.15e %s", sc, en, "  ", en_eV, "  ", en_requested, "\n")
    end

    return( nothing )
end


"""
`Basics.tabulate(stream::IO, ::HfEnergiesRelative, multiplet::Hfs.HfMultiplet)`
    ... tabulates the hyperfine level energies relative to the lowest level of the multiplet; nothing is returned.
"""
function Basics.tabulate(stream::IO, ::HfEnergiesRelative, multiplet::Hfs.HfMultiplet)
    println(stream, "\n  Energy of each level relative to lowest level for nuclear spin I = $(multiplet.levelFs[1].I):")
    sb = "  Level  F Parity          Hartrees       " * "             eV                   " * TableStrings.inUnits("energy")
    println(stream, "\n", sb, "\n")
    for  i = 2:length(multiplet.levelFs)
        lev          = multiplet.levelFs[i]
        en           = lev.energy - multiplet.levelFs[1].energy
        en_eV        = Defaults.convertUnits("energy: from atomic to eV", en)
        en_requested = Defaults.convertUnits("energy: from atomic",       en)
        sc    = " " * TableStrings.level(i) * "    " * string(LevelSymmetry(lev.F, lev.parity)) * "    "
        @printf(stream, "%s %.15e %s %.15e %s %.15e %s", sc, en, "  ", en_eV, "  ", en_requested, "\n")
    end

    return( nothing )
end


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################


"""
`Hfs.amplitude(mp::EmMultipole, rLevel::Level, sLevel::Level, grid::Radial.Grid; printout::Bool=true)`
    ... to compute the T^(M1), T^(E2) or T^(M3) hyperfine amplitude <alpha_r J_r || T^(mp) || alpha_s J_s>
        for a given pair of levels. A value::ComplexF64 is returned.
"""
function  amplitude(mp::EmMultipole, rLevel::Level, sLevel::Level, grid::Radial.Grid; printout::Bool=true)
    if  rLevel.parity != sLevel.parity   return( ComplexF64(0.) )   end
    nr = length(rLevel.basis.csfs);    ns = length(sLevel.basis.csfs);    matrix = zeros(ComplexF64, nr, ns)
    if  printout   printstyled("Compute hyperfine $(string(mp)) matrix of dimension $nr x $ns ... \n", color=:light_green)   end
    #
    if        mp == M1    opa = SpinAngular.OneParticleOperator(1, plus, true);  hfsFunc = InteractionStrength.hfs_tM1
    elseif    mp == E2    opa = SpinAngular.OneParticleOperator(2, plus, true);  hfsFunc = InteractionStrength.hfs_tE2
    elseif    mp == M3    opa = SpinAngular.OneParticleOperator(3, plus, true);  hfsFunc = InteractionStrength.hfs_tM3
    else      error("Hfs.amplitude: unsupported nuclear multipole $mp.")
    end
    for  r = 1:nr
        for  s = 1:ns
            me = 0.
            if  rLevel.basis.csfs[r].parity  != rLevel.parity    ||  sLevel.basis.csfs[s].parity  != sLevel.parity  ||
                rLevel.parity != sLevel.parity    continue
            end
            subshellList = sLevel.basis.subshells
            wa           = SpinAngular.computeCoefficients(opa, rLevel.basis.csfs[r], sLevel.basis.csfs[s], subshellList)
            for  coeff in wa
                tamp  = hfsFunc(rLevel.basis.orbitals[coeff.a], sLevel.basis.orbitals[coeff.b], grid)
                me = me + coeff.T * tamp
            end
            matrix[r,s] = me
        end
    end
    if  printout   printstyled("done.\n", color=:light_green)   end
    amplitude = transpose(rLevel.mc) * matrix * sLevel.mc
    return( amplitude )
end


"""
`Hfs.computeInteractionAmplitudeM(mp::EmMultipole, leftIsomer::Nuclear.Isomer, rightIsomer::Nuclear.Isomer)` 
    ... to compute the hyperfine interaction amplitude (<leftIsomer || M^(mp)) || rightIsomer>) for the interaction of two
        nuclear levels; this ME is geometrically fixed if the left and right isomer are the same, and it depends
        on the nuclear ME otherwise. An amplitude::ComplexF64 is returned.
"""
function  computeInteractionAmplitudeM(mp::EmMultipole, leftIsomer::Nuclear.Isomer, rightIsomer::Nuclear.Isomer)
    amplitude = 1.
    # Calculate the geometrical factor if the left- and right-hand isomer is the same
    if  leftIsomer == rightIsomer
        floatI = Basics.twice(leftIsomer.spinI) / 2.
        if       mp == M1       amplitude = leftIsomer.mu * sqrt( (floatI + 1)*(2*floatI+1) / floatI)
        elseif   mp == E2       amplitude = leftIsomer.Q / 2 * sqrt( (floatI + 1)*(2*floatI+1) * (2*floatI + 3)/ (floatI * (2*floatI -1)) )
        else   error("stop a; mp = $mp")
        end
    else
        if mp in leftIsomer.multipoleM && mp in rightIsomer.multipoleM
            lidx = findall(==(mp), leftIsomer.multipoleM)
            ridx = findall(==(mp), rightIsomer.multipoleM)
            if lidx != ridx
                error("stop a; leftIsomer.multipoleM != rightIsomer.multipoleM ") 
            else
                if length(lidx) == 1 && length(ridx) == 1 
                    amplitude = (leftIsomer.elementM[lidx[1]] + rightIsomer.elementM[ridx[1]]) / 2
                    if rightIsomer.energy < leftIsomer.energy
                        amplitude =amplitude *(-1)^(Basics.twice(rightIsomer.spinI)/2-Basics.twice(leftIsomer.spinI)/2)
                    end
                else error("stop b; leftIsomer.multipoleM setting error")
                end
            end
        else
            amplitude = 0.
        end  
    end  
    return( amplitude )
end


"""
`Hfs.computeInteractionAmplitudeT(mp::EmMultipole, aLevel::Level, bLevel, grid::Radial.Grid)` 
    ... to compute the T^(mp) interaction amplitude for two levels of the same basis, i.e. (<aLevel || T^(mp) || bLevel>).
        A me::ComplexF64 is returned.
"""
function  computeInteractionAmplitudeT(mp::EmMultipole, aLevel::Level, bLevel, grid::Radial.Grid)
    #
    ncsf = length(aLevel.basis.csfs);  me = ComplexF64(0.)
    if  ncsf != length(bLevel.basis.csfs)  ||  aLevel.basis.subshells != bLevel.basis.subshells
        error("stop a: both levels must refer to the same electronic basis.")
    end 
    
    # Compute the  T^(mp) matrix element
    for  (ia, csfa)  in  enumerate(aLevel.basis.csfs)
        for  (ib, csfb)  in  enumerate(bLevel.basis.csfs)
            wb  = ComplexF64(0.)
            if  abs(aLevel.mc[ia] * bLevel.mc[ib]) > 1.0e-10
                subshellList = aLevel.basis.subshells
                orbitals     = aLevel.basis.orbitals
                opa = SpinAngular.OneParticleOperator(mp.L, plus, true)
                wa  = SpinAngular.computeCoefficients(opa, aLevel.basis.csfs[ia], bLevel.basis.csfs[ib], subshellList)
                    for  coeff in wa
                        ja   = Basics.subshell_2j(orbitals[coeff.a].subshell)
                        jb   = Basics.subshell_2j(orbitals[coeff.b].subshell)
                        if     mp == M1   
                            if  aLevel.basis.csfs[ia].parity  != bLevel.basis.csfs[ib].parity        tamp = 0.
                            else       tamp = InteractionStrength.hfs_tM1(orbitals[coeff.a], orbitals[coeff.b], grid)   end                            
                        elseif  mp == E2
                            if  aLevel.basis.csfs[ia].parity  != bLevel.basis.csfs[ib].parity        tamp = 0.
                            else       tamp = InteractionStrength.hfs_tE2(orbitals[coeff.a], orbitals[coeff.b], grid)   end                       
                        elseif  mp == E1 
                            if  aLevel.basis.csfs[ia].parity  != bLevel.basis.csfs[ib].parity 
                                       tamp = InteractionStrength.hfs_tE1(orbitals[coeff.a], orbitals[coeff.b], grid)  
                            else       tamp = 0.                                                                        end
                        elseif  mp == E3   
                            if  aLevel.basis.csfs[ia].parity  != bLevel.basis.csfs[ib].parity 
                                       tamp = InteractionStrength.hfs_tE3(orbitals[coeff.a], orbitals[coeff.b], grid)
                            else       tamp = 0.                                                                        end
                        elseif  mp == M2    
                            if  aLevel.basis.csfs[ia].parity  != bLevel.basis.csfs[ib].parity
                                       tamp = InteractionStrength.hfs_tM2(orbitals[coeff.a], orbitals[coeff.b], grid)   
                            else       tamp = 0.                                                                        end
                        elseif  mp == M3  
                            if  aLevel.basis.csfs[ia].parity  != bLevel.basis.csfs[ib].parity        tamp = 0.
                            else       tamp = InteractionStrength.hfs_tM3(orbitals[coeff.a], orbitals[coeff.b], grid)   end                              
                        else    error("stop b")    
                        end 
                        # wb = wb + coeff.T * tamp   #Stephan
                        wb = wb + coeff.T * tamp/ sqrt( ja + 1) * sqrt( (Basics.twice(aLevel.J) + 1))    #Wu
                    end
            end 
            me = me + aLevel.mc[ia] * bLevel.mc[ib] * wb    
        end 
    end 

    return( me )
end 


"""
`Hfs.computeAmplitudesProperties(outcome::Hfs.Outcome, nm::Nuclear.Model, grid::Radial.Grid, settings::Hfs.Settings, 
                                 im::Hfs.InteractionMatrix) 
    ... to compute all amplitudes and properties of for a given level; an outcome::Hfs.Outcome is returned for which the 
        amplitudes and properties are evaluated explicitly.
"""
function  computeAmplitudesProperties(outcome::Hfs.Outcome, nm::Nuclear.Model, grid::Radial.Grid, settings::Hfs.Settings, 
                                      im::Hfs.InteractionMatrix)
    AIoverMu = BoverQ = CoverOmega = amplitudeM1 = amplitudeE2 = amplitudeM3 = 0.;    
    J = AngularMomentum.oneJ(outcome.Jlevel.J)
    #
    if  settings.calcM1  &&  outcome.Jlevel.J != AngularJ64(0)
        if  im.calcM1   amplitudeM1 = transpose(outcome.Jlevel.mc) * im.matrixM1 * outcome.Jlevel.mc
        else            amplitudeM1 = Hfs.amplitude(Basics.M1, outcome.Jlevel, outcome.Jlevel, grid)
        end
        wx       = Defaults.convertUnits("moment: from nuclear magneton to atomic", 1.0)
        AIoverMu = amplitudeM1 / sqrt(J * (J+1)) * wx 
    end
    #
    if  settings.calcE2  &&  outcome.Jlevel.J != AngularJ64(0)
        if  im.calcE2   amplitudeE2 = transpose(outcome.Jlevel.mc) * im.matrixE2 * outcome.Jlevel.mc
        else            amplitudeE2 = Hfs.amplitude(Basics.E2, outcome.Jlevel, outcome.Jlevel, grid)
        end
        wx       = Defaults.convertUnits("cross section: from barn to atomic unit", 1.0)
        BoverQ   = 2 * amplitudeE2 * sqrt( (2J-1) / ((J+1)*(2J+3)) ) * wx  ## * sqrt(J)
    end
    #
    if  settings.calcM3  &&  outcome.Jlevel.J != AngularJ64(0)
        if  im.calcM3   amplitudeM3 = transpose(outcome.Jlevel.mc) * im.matrixM3 * outcome.Jlevel.mc
        else            amplitudeM3 = Hfs.amplitude(Basics.M3, outcome.Jlevel, outcome.Jlevel, grid)
        end
        wx         = Defaults.convertUnits("moment: from nuclear magneton x fm^2 to atomic", 1.0)
        CoverOmega =   - amplitudeM3 * sqrt( J *(J-1)*(2J-1) / ((J+1)*(J+2)*(2J+3)) ) * wx   
    end
    #
    hfMultiplet = Hfs.HfMultiplet()
    if  settings.calcHfMultiplet
        # Determine a HfMultiplet for the given Jlevel/outcome
        error("... still to be done for a single nuclear spin/isomer")
        hfsMultiplet = Hfs.computeHyperfineMultiplet(outcome.Jlevel, nm, grid)
    end
    newOutcome = Hfs.Outcome( outcome.Jlevel, 1., AIoverMu, BoverQ, CoverOmega, amplitudeM1, amplitudeE2, amplitudeM3, 
                              nm.spinI, hfMultiplet)
    return( newOutcome )
end


"""
`Hfs.computeHyperfineMultiplet(level::Level, nm::Nuclear.Model, grid::Radial.Grid)`  
    ... to compute a hyperfine multiplet, i.e. a representation of hyperfine levels within a hyperfine-coupled basis as defined by the
        given (electronic) level; a hfMultiplet::hfMultiplet is returned.
"""
function computeHyperfineMultiplet(level::Level, nm::Nuclear.Model, grid::Radial.Grid)
    #
    hfBasis     = Hfs.defineHyperfineBasis(level, nm)
    hfMultiplet = Hfs.computeHyperfineRepresentation(hfBasis, nm, grid)

    return( hfMultiplet )
end


"""
`Hfs.computeHyperfineMultiplet(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::Hfs.Settings; output=true)`  
    ... to compute a hyperfine multiplet, i.e. a representation of hyperfine levels within a hyperfine-coupled basis as defined by the
        given (electronic) multiplet; a hfsMultiplet::IJF_Multiplet is returned.
"""
function computeHyperfineMultiplet(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::Hfs.Settings; output=true)
    println("")
    printstyled("Hfs.computeHyperfineMultiplet(): The computation of the hyperfine multiplet starts now ... \n", color=:light_green)
    printstyled("------------------------------------------------------------------------------------------ \n", color=:light_green)
    #
    hfBasis     = Hfs.defineHyperfineBasis(multiplet, nm)
    hfMultiplet = Hfs.computeHyperfineRepresentation(hfBasis, nm, grid)
    # Print all results to screen
    hfMultiplet = Basics.sortByEnergy(hfMultiplet)
    Basics.tabulate(stdout,   HfEnergies(),         hfMultiplet)
    Basics.tabulate(stdout,   HfEnergiesRelative(), hfMultiplet)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary
        Basics.tabulate(iostream, HfEnergies(),         hfMultiplet)
        Basics.tabulate(iostream, HfEnergiesRelative(), hfMultiplet)
    end
    #
    if    output    return( hfMultiplet )
    else            return( nothing )
    end
end


"""
`Hfs.computeHyperfineRepresentation(hfBasisVectors::Array{HfBasisVector,1}, nm::Nuclear.Model, grid::Radial.Grid)`  
    ... to set-up and diagonalized the Hamiltonian matrix of H^(DFB) + H^(hfs) within the atomic hyperfine (IJF-coupled) basis;
        a hfsMultiplet::IJF_Multiplet is returned.
"""
function computeHyperfineRepresentation(hfBasisVectors::Array{HfBasisVector,1}, nm::Nuclear.Model, grid::Radial.Grid)
    n = length(hfBasisVectors);   matrix = zeros(n,n)
    for  r = 1:n
        for  s = 1:n
            matrix[r,s] = 0.
            if  r == s    matrix[r,s] = matrix[r,s] +  hfBasis.vectors[r].levelJ.energy                 end
            if  hfBasis.vectors[r].F               !=  hfBasis.vectors[s].F                continue     end
            if  hfBasis.vectors[r].levelJ.parity   !=  hfBasis.vectors[s].levelJ.parity    continue     end
            # Now add the hyperfine matrix element for K = 1 (magnetic) and K = 2 (electric) contributions
            spinI = nm.spinI;   Jr = hfBasis.vectors[r].levelJ.J;   Js = hfBasis.vectors[s].levelJ.J;   F = hfBasis.vectors[r].F
            Ix    = AngularMomentum.oneJ(spinI)
            #
            wa = AngularMomentum.phaseFactor([spinI, +1, Jr, +1, F])  * 
                    AngularMomentum.Wigner_6j(spinI, Jr, F, Js, spinI, AngularJ64(1))
            wb = Hfs.amplitude(Basics.M1, hfBasis.vectors[r].levelJ, hfBasis.vectors[s].levelJ, grid::Radial.Grid)
            wc = nm.mu * sqrt( (Ix+1)/Ix )
            matrix[r,s] = matrix[r,s] + wa * wb * wc     * 1.0e-6 ## fudge-factor to keep HFS interaction small
            #
            if  spinI  in [AngularJ64(0), AngularJ64(1//2)]                                 continue     end
            wa = AngularMomentum.phaseFactor([spinI, +1, Jr, +1, F]) * 
                    AngularMomentum.Wigner_6j(spinI, Jr, F, Js, spinI, AngularJ64(2))
            wb = Hfs.amplitude(Basics.E2, hfBasis.vectors[r].levelJ, hfBasis.vectors[s].levelJ, grid::Radial.Grid)
            wc = nm.Q / 2. * sqrt( (Ix+1)*(2Ix+3)/ (Ix*(2Ix-1)) )
            matrix[r,s] = matrix[r,s] + wa * wb * wc     * 1.0e-6 ## fudge-factor to keep HFS interaction small
        end
    end
    #
    # Diagonalize the matrix and set-up the representation
    eigen    = Basics.diagonalize("matrix: LinearAlgebra", matrix)
    levelFs  = Hfs.HfLevel[]
    for  ev = 1:length(eigen.values)
        # Construct the eigenvector with regard to the given basis (not w.r.t the symmetry block)
        evector   = eigen.vectors[ev];    en = eigen.values[ev]
        parity    = Basics.plus;    F = AngularJ64(0);     MF = AngularM64(0)
        for  r = 1:length(hfBasis.vectors)
            if  abs(evector[r]) > 1.0e-6    
                parity = hfBasis.vectors[r].levelJ.parity
                F      = hfBasis.vectors[r].F;     MF = AngularM64(hfBasis.vectors[r].F);   break    end
        end
        newlevelF = Hfs.HfLevel(nm.spinI, F, MF, parity, en, hfBasis, evector) 
        push!( levelFs, newlevelF)
    end
    hfMultiplet = Hfs.HfMultiplet("hyperfine", levelFs)
    
    return( hfMultiplet )
end


"""
`Hfs.computeInteractionMatrix(basis::Basis, grid::Radial.Grid, settings::Hfs.Settings)` 
    ... to compute the T^M1 and/or T^E2 interaction matrices for the given basis, i.e. (<csf_r || T^(n)) || csf_s>).
        An im::Hfs.InteractionMatrix is returned.
"""
function  computeInteractionMatrix(basis::Basis, grid::Radial.Grid, settings::Hfs.Settings)
    #
    ncsf = length(basis.csfs);    matrixM1 = zeros(ncsf,ncsf);    matrixE2 = zeros(ncsf,ncsf);    matrixM3 = zeros(ncsf,ncsf)
    #
    if  settings.calcM1
        calcM1 = true;    matrixM1 = zeros(ncsf,ncsf)
        for  r = 1:ncsf
            for  s = 1:ncsf
                if  basis.csfs[r].parity  != basis.csfs[s].parity   continue    end 
                # Calculate the spin-angular coefficients
                subshellList = basis.subshells
                opa = SpinAngular.OneParticleOperator(1, plus, true)
                wa  = SpinAngular.computeCoefficients(opa, basis.csfs[r], basis.csfs[s], subshellList) 
                #
                for  coeff in wa
                    ja   = Basics.subshell_2j(basis.orbitals[coeff.a].subshell)
                    jb   = Basics.subshell_2j(basis.orbitals[coeff.b].subshell)
                    tamp = InteractionStrength.hfs_tM1(basis.orbitals[coeff.a], basis.orbitals[coeff.b], grid)
                    matrixM1[r,s] = matrixM1[r,s] + coeff.T * tamp  
                end
            end
        end
    else   
        calcM1 = false;    matrixM1 = zeros(2,2)
    end
    #
    if  settings.calcE2
        calcE2 = true;    matrixE2 = zeros(ncsf,ncsf)
        for  r = 1:ncsf
            for  s = 1:ncsf
                if  basis.csfs[r].parity  != basis.csfs[s].parity   continue    end 
                 # Calculate the spin-angular coefficients
                subshellList = basis.subshells
                opa = SpinAngular.OneParticleOperator(2, plus, true)
                wa  = SpinAngular.computeCoefficients(opa, basis.csfs[r], basis.csfs[s], subshellList) 
                #
                for  coeff in wa
                    ja   = Basics.subshell_2j(basis.orbitals[coeff.a].subshell)
                    jb   = Basics.subshell_2j(basis.orbitals[coeff.b].subshell)
                    tamp  = InteractionStrength.hfs_tE2(basis.orbitals[coeff.a], basis.orbitals[coeff.b], grid)
                    matrixE2[r,s] = matrixE2[r,s] + coeff.T * tamp  
                end
            end
        end
    else   
        calcE2 = false;    matrixE2 = zeros(2,2)
    end
    #
    if  settings.calcM3
        calcM3 = true;    matrixM3 = zeros(ncsf,ncsf)
        for  r = 1:ncsf
            for  s = 1:ncsf
                if  basis.csfs[r].parity  != basis.csfs[s].parity   continue    end 
                 # Calculate the spin-angular coefficients
                subshellList = basis.subshells
                opa = SpinAngular.OneParticleOperator(3, plus, true)
                wa  = SpinAngular.computeCoefficients(opa, basis.csfs[r], basis.csfs[s], subshellList) 
                #
                for  coeff in wa
                    ja   = Basics.subshell_2j(basis.orbitals[coeff.a].subshell)
                    jb   = Basics.subshell_2j(basis.orbitals[coeff.b].subshell)
                    tamp  = InteractionStrength.hfs_tE2(basis.orbitals[coeff.a], basis.orbitals[coeff.b], grid)
                    matrixM3[r,s] = matrixM3[r,s] + coeff.T * tamp  
                end
            end
        end
    else   
        calcM3 = false;    matrixM3 = zeros(2,2)
    end
    #
    im = Hfs.InteractionMatrix(calcM1, calcE2, calcM3, matrixM1, matrixE2, matrixM3)
    #
    return( im )
end


"""
`Hfs.computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::Hfs.Settings; output=true)`  
    ... to compute (as selected) the HFS A, B and C parameters as well as hyperfine energy splittings for the levels 
        of the given multiplet and as specified by the given settings. The results are printed in neat tables to 
        screen and, if requested, an arrays{Hfs.Outcome,1} with all the results are returned.
"""
function computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::Hfs.Settings; output=true)
    println("")
    printstyled("Hfs.computeOutcomes(): The computation of the Hyperfine amplitudes and parameters starts now ... \n", color=:light_green)
    printstyled("------------------------------------------------------------------------------------------------ \n", color=:light_green)
    println("")
    outcomes = Hfs.determineOutcomes(multiplet, settings)
    # Display all selected levels before the computations start
    if  settings.printBefore    Hfs.displayOutcomes(outcomes)    end
    # Calculate all amplitudes and requested properties
    im = Hfs.computeInteractionMatrix(multiplet.levels[1].basis, grid, settings)
    newOutcomes = Hfs.Outcome[]
    for  outcome in outcomes
        newOutcome = Hfs.computeAmplitudesProperties(outcome, nm, grid, settings, im) 
        push!( newOutcomes, newOutcome)
    end
    # Print all results to screen
    Hfs.displayResults(stdout, newOutcomes, nm, settings)
    # Compute and display the non-diagonal hyperfine amplitudes, if requested
    if  settings.calcNondiagonal    Hfs.displayNondiagonal(stdout, multiplet, grid, settings)   end
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    
        Hfs.displayResults(iostream, newOutcomes, nm, settings) 
        if  settings.calcNondiagonal    Hfs.displayNondiagonal(iostream, multiplet, grid, settings)   end
    end
    #
    if    output    return( newOutcomes )
    else            return( nothing )
    end
end


"""
`Hfs.defineHyperfineBasis(level::Level, nm::Nuclear.Model)`  
    ... to define/set-up an atomic hyperfine (IJF-coupled) basis for the given electronic level; 
        a hfsBasis::IJF_Basis is returned.
"""
function  defineHyperfineBasis(level::Level, nm::Nuclear.Model)
    function  display_ijfVector(i::Int64, vector::Hfs.IJF_Vector) 
        si = string(i);   ni = length(si);    sa = repeat(" ", 5);    sa = sa[1:5-ni] * si * ")  "
        sa = sa * "[" * string( LevelSymmetry(vector.levelJ.J, vector.levelJ.parity) ) * "] " * string(vector.F) * repeat(" ", 4)
        return( sa )
    end

    vectors = HfVector[]
    Flist = Basics.oplus(nm.spinI, level.J)
    for  F in Flist     push!(vectors,  HfVector(nm.spinI, F, level) )      end
    #
    hfBasis = HfBasis(vectors, level.basis)
    
    return( hfBasis )
end


"""
`Hfs.defineHyperfineBasis(multiplet::Multiplet, nm::Nuclear.Model)`  
    ... to define/set-up an atomic hyperfine (IJF-coupled) basis for the given electronic multipet; 
        a hfBasis::HfBasis is returned.
"""
function  defineHyperfineBasis(multiplet::Multiplet, nm::Nuclear.Model)
    function  display_ijfVector(i::Int64, vector::Hfs.IJF_Vector) 
        si = string(i);   ni = length(si);    sa = repeat(" ", 5);    sa = sa[1:5-ni] * si * ")  "
        sa = sa * "[" * string( LevelSymmetry(vector.levelJ.J, vector.levelJ.parity) ) * "] " * string(vector.F) * repeat(" ", 4)
        return( sa )
    end


    vectors = HfVector[]
    #
    for i = 1:length(multiplet.levels)
        Flist = Basics.oplus(nm.spinI, multiplet.levels[i].J)
        for  F in Flist
            push!(vectors,  IJF_Vector(nm.spinI, F, multiplet.levels[i]) )
        end
    end
    #
    hfBasis = HfBasis(vectors, multiplet.levels[1].basis)
    #
    println(" ")
    println("  Construction of a atomic hyperfine (IJF-coupled) basis of dimension $(length(hfBasis.vectors)), based on " *
            "$(length(multiplet.levels)) electronic levels and nuclear spin $(hfBasis.vectors[1].I).")
    println("  Basis vectors are defined as [J^P] F: \n")
    sa = "  "
    for  k = 1:length(hfBasis.vectors)
        if  length(sa) > 100    println(sa);   sa = "  "   end
        sa = sa * display_ijfVector(k, hfBasis.vectors[k])
    end
    if   length(sa) > 5    println(sa)   end
    println(" ")
    
    return( hfBasis )
end


"""
`Hfs.computeModifiedEinsteinRates(upperOutcome::Outcome, lowerOutcome::Outcome, multipoles::Array{EmMultipole,1}, gauge::EmGauge,
                                  grid::Radial.Grid)`  
    ... to compute and tabulate the modified Einstein amplitudes and rates for the hyperfine-resolved 
        transitions between the upper and lower outcome. The procedures assumes that the two outcomes provide
        a proper IJF expansion (multiplet) of the hyperfine levels of interest.
        A neat table is printed but nothing is returned otherwise
"""
function  computeModifiedEinsteinRates(upperOutcome::Outcome, lowerOutcome::Outcome, multipoles::Array{EmMultipole,1}, gauge::EmGauge,
                                       grid::Radial.Grid)
    if  upperOutcome.nuclearI != lowerOutcome.nuclearI   
            error("Inconsistent nuclear spins; upper-I=$(upperOutcome.nuclearI)  !=  lower-I=$(lowerOutcome.nuclearI)")
    end
    
    stream      = stdout
    amplitudesJ = ComplexF64[]
    iJsym = LevelSymmetry(upperOutcome.Jlevel.J, upperOutcome.Jlevel.parity)
    fJsym = LevelSymmetry(lowerOutcome.Jlevel.J, lowerOutcome.Jlevel.parity)
    omega = upperOutcome.Jlevel.energy - lowerOutcome.Jlevel.energy
    
    # First compute the multipole transition amplitudes between the J-levels of the upper and lower outcome
    println(stream, " ")
    println(stream, "  Multipole amplitudes for J-levels with symmetry $iJsym --> $fJsym:")
    println(stream, " ")
    for  multipole in multipoles
        if  multipole  in  [M1, M2, M3]    gaugex = Basics.Magnetic    else    gaugex = gauge  end
        ampJ = PhotoEmission.amplitude(Emission(), multipole, gaugex, omega, upperOutcome.Jlevel, lowerOutcome.Jlevel, grid)
        push!(amplitudesJ, ampJ)
        println(stream, "  J-level amplitude for $multipole  = $ampJ  ")
    end
    
    # Print a table header
    nx = 150
    println(stream, " ")
    println(stream, "  HFS modified Einstein amplitudes and rates:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = " ";   sb = "  "
    sa = sa * TableStrings.center(14, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(16)
    sa = sa * TableStrings.center(14, "i--J^P--f"; na=2);                         sb = sb * TableStrings.hBlank(16)
    sa = sa * TableStrings.center(12, "i--F--f";   na=2);                         sb = sb * TableStrings.hBlank(15)
    sa = sa * TableStrings.center(12, "Energy"   ; na=3);               
    sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=3)
    sa = sa * TableStrings.center( 9, "Multipole"; na=1);                         sb = sb * TableStrings.hBlank(10)
    sa = sa * TableStrings.center(11, "Gauge"    ; na=4);                         sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(26, "A--Einstein--B"; na=3);       
    sb = sb * TableStrings.center(26, TableStrings.inUnits("rate")*"           "*TableStrings.inUnits("rate"); na=2)
    sa = sa * TableStrings.center(26, "re-- <Ff |amplitude L| Fi> --im"; na=2);       
    
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #  
    for  upperLevF in  upperOutcome.hfMultiplet.levelFs
        for  lowerLevF in  lowerOutcome.hfMultiplet.levelFs
            for  (m, mp)  in  enumerate(multipoles)
                Ji = upperOutcome.Jlevel.J;   Fi = upperLevF.F
                Jf = lowerOutcome.Jlevel.J;   Ff = lowerLevF.F
                
                ampF = ComplexF64(0.)
                for  umc in upperLevF.mc,   lmc in lowerLevF.mc
                    ampF = ampF  +  umc * lmc * sqrt( AngularMomentum.bracket([Fi, Ff]) ) *
                           AngularMomentum.phaseFactor([Ji, +1, upperOutcome.nuclearI, +1, Ff, +1, AngularJ64(mp.L)]) * 
                           AngularMomentum.AngularMomentum.Wigner_6j(Fi, Ff, mp.L, Jf, Ji, upperOutcome.nuclearI) * amplitudesJ[m]
                end
                sa = ""
                sa = sa * TableStrings.center(14, TableStrings.levels_if(upperOutcome.Jlevel.index, lowerOutcome.Jlevel.index); na=2)
                sa = sa * TableStrings.center(14, TableStrings.symmetries_if(iJsym, iJsym); na=0)
                sc = "         " * string(Fi) * "    " * string(Ff)
                sa = sa * TableStrings.center(12, sc[end-9:end]; na=3)
                en = upperOutcome.Jlevel.energy - lowerOutcome.Jlevel.energy
                sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", en)) * "    "
                sa = sa * TableStrings.center(9,  string(mp); na=4)
                if  mp  in  [M1, M2, M3]    gaugex = Basics.Magnetic    else    gaugex = gauge  end
                sa = sa * TableStrings.flushleft(11, string(gaugex);  na=0)
                chRate =  8pi * Defaults.getDefaults("alpha") * en / (Basics.twice(Ji) + 1) * (abs(ampF)^2) * (Basics.twice(Jf) + 1)
                ## sa = sa * @sprintf("%.6e", Basics.recast("rate: radiative, to Einstein A",  line, chRate)) * "  "
                ## sa = sa * @sprintf("%.6e", Basics.recast("rate: radiative, to Einstein B",  line, chRate)) * "    "
                sa = sa * @sprintf("% .6e", chRate)  * "  " 
                sa = sa * @sprintf("% .6e", chRate)  * "    " 
                sa = sa * @sprintf("% .6e", ampF.re) * "  " 
                sa = sa * @sprintf("% .6e", ampF.im) * "  " 
                println(stream, sa)
            end
            println(stream, " ")
        end
    end
    println(stream, "  ", TableStrings.hLine(nx)) 
    
    return( nothing )
end


"""
`Hfs.determineOutcomes(multiplet::Multiplet, settings::Hfs.Settings)`  
    ... to determine a list of Outcomes's for the computation of HFS A- and B-parameters for the given multiplet. 
        It takes into account the particular selections and settings. An Array{Hfs.Outcome,1} is returned. Apart from the 
        level specification, all physical properties are set to zero during the initialization process.
"""
function  determineOutcomes(multiplet::Multiplet, settings::Hfs.Settings) 
    outcomes = Hfs.Outcome[]
    for  level  in  multiplet.levels
        if  Basics.selectLevel(level, settings.levelSelection)
            push!( outcomes, Hfs.Outcome(level, 0., 0., 0., 0., 0., 0., 0., AngularJ64(0), Hfs.HfMultiplet() ) )
        end
    end
    return( outcomes )
end


"""
`Hfs.displayNondiagonal(stream::IO, multiplet::Multiplet, grid::Radial.Grid, settings::Hfs.Settings)`  
    ... to compute and display all non-diagonal hyperfine amplitudes for the selected levels. A small neat table of 
        all (pairwise) hyperfine amplitudes is printed but nothing is returned otherwise.
"""
function  displayNondiagonal(stream::IO, multiplet::Multiplet, grid::Radial.Grid, settings::Hfs.Settings)
    # Determine pairs to be calculated
    pairs = Tuple{Int64,Int64}[]
    for  (f, fLevel)  in  enumerate(multiplet.levels)
        for  (i, iLevel)  in  enumerate(multiplet.levels)
            if  Basics.selectLevel(fLevel, settings.levelSelection)  &&   Basics.selectLevel(iLevel, settings.levelSelection)
                push!( pairs, (f,i) )
            end
        end
    end
    #
    nx = 107
    println(stream, " ")
    println(stream, "  Selected (non-) diagonal hyperfine amplitudes:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  "
    sa = sa * TableStrings.center(10, "Level_f"; na=2)
    sa = sa * TableStrings.center(10, "Level_i"; na=2)
    sa = sa * TableStrings.center(10, "J^P_f";   na=3)
    sb = repeat(" ", length(sa))
    sa = sa * TableStrings.center(70, "Amplitudes";                                     na=4);              
    sb = sb * TableStrings.center(70, "M1        ----        E2        ----        M3"; na=4);              
    sa = sa * TableStrings.center(10, "J^P_f";   na=4)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #  
    for  (f,i) in pairs
        sa   = "  ";    
        sa   = sa * TableStrings.center(10, string(f); na=2)
        sa   = sa * TableStrings.center(10, string(i); na=2)
        symf = LevelSymmetry( multiplet.levels[f].J, multiplet.levels[f].parity)
        symi = LevelSymmetry( multiplet.levels[i].J, multiplet.levels[i].parity)
        sa   = sa * TableStrings.center(10, string(symf); na=4)
        M1   = Hfs.amplitude(Basics.M1, multiplet.levels[f], multiplet.levels[i], grid; printout=false)
        E2   = Hfs.amplitude(Basics.E2, multiplet.levels[f], multiplet.levels[i], grid; printout=false)
        M3   = Hfs.amplitude(Basics.M3, multiplet.levels[f], multiplet.levels[i], grid; printout=false)
        sa   = sa * @sprintf("%.5e %s %.5e", M1.re, "  ", M1.im) * "    "
        sa   = sa * @sprintf("%.5e %s %.5e", E2.re, "  ", E2.im) * "    "
        sa   = sa * @sprintf("%.5e %s %.5e", M3.re, "  ", M3.im) * "    "
        sa   = sa * TableStrings.center(10, string(symi); na=4)
        println(stream, sa )
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`Hfs.displayOutcomes(outcomes::Array{Hfs.Outcome,1})`  
    ... to display a list of levels that have been selected for the computations A small neat table of all 
        selected levels and their energies is printed but nothing is returned otherwise.
"""
function  displayOutcomes(outcomes::Array{Hfs.Outcome,1})
    nx = 43
    println(" ")
    println("  Selected HFS levels:")
    println(" ")
    println("  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "Energy"; na=4);              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
    #  
    for  outcome in outcomes
        sa  = "  ";    sym = LevelSymmetry( outcome.Jlevel.J, outcome.Jlevel.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(outcome.Jlevel.index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", outcome.Jlevel.energy)) * "    "
        println( sa )
    end
    println("  ", TableStrings.hLine(nx))
    println(" ")
    #
    return( nothing )
end


"""
`Hfs.displayResults(stream::IO, outcomes::Array{Hfs.Outcome,1}, nm::Nuclear.Model, settings::Hfs.Settings)`  
    ... to display the energies, A- and B-values, Delta E_F energy shifts, etc. for the selected levels. All nuclear 
        parameters are taken from the nuclear model. A neat table is printed but nothing is returned otherwise.
"""
function  displayResults(stream::IO, outcomes::Array{Hfs.Outcome,1}, nm::Nuclear.Model, settings::Hfs.Settings)
    nx = 128
    println(stream, " ")
    println(stream, "  HFS amplitudes and g_J factors:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "Energy"; na=10);              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=10)
    sa = sa * TableStrings.center(54, "M1  -- Amplitudes --   E2  -- Amplitudes --   M3"    ; na=8);        
    sb = sb * TableStrings.hBlank(66)
    sa = sa * TableStrings.center(12, "g_J"; na=2);              
    sb = sb * TableStrings.center(12, "   "; na=2); 
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #  
    for  outcome in outcomes
        sa  = "  ";    sym = LevelSymmetry( outcome.Jlevel.J, outcome.Jlevel.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(outcome.Jlevel.index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=4)
        energy = outcome.Jlevel.energy
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", energy))                        * "      "
        sa = sa * @sprintf("% .8e %s % .8e %s % .8e", outcome.amplitudeM1.re, "      ", outcome.amplitudeE2.re, 
                                                                              "      ", outcome.amplitudeM3.re) * "      "
        sa = sa * @sprintf("%.8e", outcome.gJ)                                                                  * "    " 
        println(stream, sa )
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    nx = 140
    println(stream, " ")
    println(stream, "  HFS parameters:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "Energy"; na=10);              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=10)
    sa = sa * TableStrings.center(85, "  A     ---    A/mu     ---     B    ---     B/Q     ---      C    ---     C/Omega     "; na=2);              
    sb = sb * TableStrings.center(85, "[MHz]       [MHz/mu_nuc]      [MHz]       [MHz/barn]        [MHz]     [MHz/mu_nuc fm^2]"; na=2); 
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #  
    for  outcome in outcomes
        sa  = "  ";    sym = LevelSymmetry( outcome.Jlevel.J, outcome.Jlevel.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(outcome.Jlevel.index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=4)
        energy = outcome.Jlevel.energy
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", energy))           * "      "
        we = Defaults.convertUnits("energy: from atomic to Hz", 1.0) / 1.0e9      # Energy factor into MHz
        wa = outcome.AIoverMu / nm.spinI.num * nm.spinI.den * nm.mu * we          # Prepare A
        sa = sa * @sprintf("% .6e", wa)       * "  " 
        wa = outcome.AIoverMu / nm.spinI.num * nm.spinI.den * we                  # Prepare A/mu
        sa = sa * @sprintf("% .6e", wa)       * "  " 
        wa = outcome.BoverQ * nm.Q * we                                           # Prepare B
        sa = sa * @sprintf("% .6e", wa)       * "  " 
        wa = outcome.BoverQ * we                                                  # Prepare B/Q
        sa = sa * @sprintf("% .6e", wa)       * "  " 
        wa = outcome.CoverOmega * nm.Omega * we                                   # Prepare C
        sa = sa * @sprintf("% .6e", wa)       * "  " 
        wa = outcome.CoverOmega * we                                              # Prepare C/Omega
        sa = sa * @sprintf("% .6e", wa)       * "  " 
        println(stream, sa )
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    #
    # Printout the Delta E_F energy shifts of the hyperfine levels |alpha F> with regard to the (electronic) levels |alpha J>
    nx = 90
    println(stream, " ")
    println(stream, "  HFS Delta E_F energy shifts with regard to the (electronic) level energies E_J:")
    println(stream, " ")
    println(stream, "    Nuclear spin I:                             $(nm.spinI) ")
    println(stream, "    Nuclear magnetic-dipole moment      mu    = $(nm.mu)    ")
    println(stream, "    Nuclear electric-quadrupole moment  Q     = $(nm.Q)     ")
    println(stream, "    Nuclear magnetic-octupole moment    Omega = $(nm.Omega) ")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "Energy"; na=4);                            sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(10, "F^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "Delta E_F"; na=4);                         
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(14, "K factor"; na=4);                         
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #  
    for  outcome in outcomes
        sa  = "  ";    sym = LevelSymmetry( outcome.Jlevel.J, outcome.Jlevel.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(outcome.Jlevel.index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=4)
        energy = outcome.Jlevel.energy
        #
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", energy))           * "    "
        J     = AngularMomentum.oneJ(outcome.Jlevel.J);   spinI = AngularMomentum.oneJ(nm.spinI)
        Flist = Basics.oplus(nm.spinI, outcome.Jlevel.J)
        first = true
        for  Fang in Flist
            Fsym    = LevelSymmetry(Fang, outcome.Jlevel.parity)
            F       = AngularMomentum.oneJ(Fang)
            Kfactor = F*(F+1) - J*(J+1) - spinI*(spinI+1)
            energy  = outcome.AIoverMu * nm.mu / spinI * Kfactor / 2.
            if  abs(outcome.BoverQ) > 1.0e-10
                energy  = energy +  outcome.BoverQ * nm.Q * 3/4 * (Kfactor*(Kfactor+1) - spinI*(spinI+1)*J*(J+1) ) /
                                    ( 2spinI*(2spinI-1)*J*(2J-1) )
            end
            sb = TableStrings.center(10, string(Fsym); na=2)
            sb = sb * TableStrings.flushright(16, @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", energy))) * "    "
            sb = sb * TableStrings.flushright(12, @sprintf("%.5e", Kfactor))
            #
            if   first    println(stream,  sa*sb );   first = false
            else          println(stream,  TableStrings.hBlank( length(sa) ) * sb )
            end
        end
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    if  settings.calcHfMultiplet
        nx = 90
        println(stream, " ")
        println(stream, "  IJF-coupled hyperfine levels:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  "
        sa = sa * TableStrings.center(10, "Level"; na=2)
        sa = sa * TableStrings.center(10, "J^P";   na=4)
        sa = sa * TableStrings.center( 6, "F^P";   na=4);   na = length(sa)   
        sa = sa * TableStrings.flushleft(34, "IJF basis (1:)"; na=4)
        println(stream, sa);    sa = repeat(" ", na);   
        sa = sa * TableStrings.flushleft(34, "Mixing coefficients (1:)"; na=4);                         
        println(stream, sa);    println(stream, "  ", TableStrings.hLine(nx)) 
        #  
        for  outcome in outcomes
            sa = "  ";    sym = LevelSymmetry( outcome.Jlevel.J, outcome.Jlevel.parity)
            sa = sa * TableStrings.center(10, TableStrings.level(outcome.Jlevel.index); na=2)
            sa = sa * TableStrings.center(10, string(sym); na=6);    na = length(sa)
            sa = sa * "        "
            if   length(outcome.hfsMultiplet.levelFs) == 0  sa = sa * "No IJF levels";   println(stream, sa);    continue    end
            nvecs = min( length(outcome.hfsMultiplet.levelFs[1].basis.vectors), 6)
            for  nvec = 1:nvecs
                sa = sa * "[" * string(outcome.hfsMultiplet.levelFs[1].basis.vectors[nvec].levelJ.J) * "] " *
                                string(outcome.hfsMultiplet.levelFs[1].basis.vectors[nvec].F) * "       "
            end
            println(stream, sa)
            for  levelF in outcome.hfsMultiplet.levelFs
                sa = repeat(" ", na)
                sa = sa * string(levelF.F) * "       "
                for  nvec = 1:nvecs
                    sa = sa * @sprintf("% .5e", levelF.mc[nvec] ) * "  "
                end
                println(stream, sa)
            end
            println(stream, " ")
        end
        println(stream, "  ", TableStrings.hLine(nx))
    end
    #
    return( nothing )
end

end # module


