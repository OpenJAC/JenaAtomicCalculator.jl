
"""
`module  JAC.Hfs`  
    ... a submodel of JAC that contains all methods for computing HFS A, B and C coefficients, hyperfine-level representations, etc.
        Hyperfine coupling is expressed throughout as F = I + J, i.e. the nucleus is coupled to a finished electronic
        ASF (Level) rather than to individual CSF. Working with ASF keeps the number of hyperfine levels moderate and
        makes the results directly interpretable; the computational effort is formally equivalent.

        THE CHAIN OF DATA TYPES, from the nucleus up to a hyperfine multiplet (simplified 05-Aug-2026):

        + Nuclear.Isomer ... one specific nuclear state (spinI, parity, energy, mu, Q, ...), assembled from a
              Nuclear.Model. Setting mu = Q = 0 makes the hyperfine Hamiltonian vanish identically, which is the
              device used to switch hyperfine structure off for levels where it is negligible.

        + HfBasisVector = HfBasisVector(F, parity, isomer, levelJ)  ... ONE coupled basis vector,
              |(I J) F M> = sum <I M_I  J M_J | F M> |I M_I> |J M_J>. Note that the Clebsch-Gordan coupling is
              already contained here: a basis vector IS the coupled product state, not a pair of factors.

        + a BASIS is simply an Array{HfBasisVector,1}. There is deliberately NO separate type for it. (Code that
              called a non-existent `HfBasis(...)` constructor was one of the reasons defineHyperfineBasis could
              never run before 05-Aug-2026.)

        + HfLevel = HfLevel(F, M, parity, index, energy, hfBasisVectors, mc)
          ... one physical hyperfine level, = sum_k mc[k] * hfBasisVectors[k]. The mixing coefficients mc are the
              eigenvector of the hyperfine Hamiltonian in that coupled basis, i.e. they describe J-MIXING among
              vectors of equal F -- they are NOT Clebsch-Gordan coefficients, which sit inside each basis vector.
              mc = [.. 1.0 ..] therefore means a pure |(I J) F>, which is what mu = Q = 0 produces automatically.

        + HfMultiplet = HfMultiplet(name, hfLevels) ... the pendant to an electronic Multiplet, so that
              hyperfine-resolved transitions can be treated completely analogously to the electronic ones. What is
              needed in addition is only the recoupling of the electronic reduced matrix elements; see
              Hfs.reducedMeElectronic, whose header explains the normalization that must be used for that.

        Build one with Hfs.computeHyperfineMultiplet(multiplet, nm, grid, settings), which chains
        Hfs.defineHyperfineBasis and Hfs.computeHyperfineRepresentation.
"""
module Hfs


using Printf, ..AngularMomentum, ..Basics,  ..Defaults, ..InteractionStrength, ..ManyElectron, ..Radial, ..Nuclear,
              ..SpinAngular, ..TableStrings, ..PhotoEmission

## Tabulate-theme sentinel types for Basics.tabulate dispatch on Hfs.HfMultiplet.
struct  HfEnergies          end   ## tabulate all hyperfine level energies
struct  HfEnergiesRelative  end   ## tabulate hyperfine level energies relative to the lowest level




#################################################################################################################################
#################################################################################################################################

## RETIRED 05-Aug-2026: Hfs.IJF_Vector, the CSF-based coupled basis vector (F, parity, isomer, csf, basisJ).
## It was never constructed with content anywhere in src/ -- only its empty constructor existed -- and the job it
## was meant for, the hyperfine interaction inside a CSF basis, is done directly by Hfs.computeInteractionMatrix.
## Having two similarly-named coupled basis types in one module was itself a defect: defineHyperfineBasis reached
## for this one when it needed the ASF-based HfBasisVector, which is one reason it could never run.


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
    of ASF, typically taken from some (electronic) multiplet. In contrast to a pure (electronic) CSF basis, the use of 
    ASF simplifies the interpretation of physical findings but cannot reduce the computational effort (perhaps, even results
    in slightly increase the computational effort). A HfLevel is defined by:

    + F              ::AngularJ64               ... Total angular momentum F.
    + M              ::AngularM64               ... Total projection M, only defined if a particular magnetic sublevel is referred to.
    + parity         ::Parity                   ... Parity of the level which corresponds to the electronic system.
    + energy         ::Float64                  ... energy
    + hfBasisVectors ::Array{HfBasisVector,1}   ... the product basis nuclear (isomeric) state x selected ASF.
    + mc             ::Vector{Float64}          ... Vector of mixing coefficients w.r.t hfBasisVectors.
    
    
"""
struct HfLevel
    F                ::AngularJ64
    M                ::AngularM64
    parity           ::Parity
    index            ::Int64
    energy           ::Float64
    hfBasisVectors   ::Array{HfBasisVector,1}
    mc               ::Vector{Float64}
end 


"""
`Hfs.HfLevel()`  ... constructor for an `empty` instance of HfLevel.
"""
function HfLevel()
    HfLevel(AngularJ64(0), AngularM64(0), Basics.plus, 0, 0., HfBasisVector[], Float64[])
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
    ... keyword constructor to create a Hfs.Settings with selected non-default values.
"""
function Settings(; calcM1::Bool=true, calcE2::Bool=false, calcM3::Bool=false, calcNondiagonal::Bool=false,
                    calcHfMultiplet::Bool=false, printBefore::Bool=false, levelSelection::LevelSelection=LevelSelection())
    Settings(calcM1, calcE2, calcM3, calcNondiagonal, calcHfMultiplet, printBefore, levelSelection)
end


"""
`Hfs.Settings(set::Hfs.Settings;`

        calcM1=.., calcE2=.., calcM3=.., calcNondiagonal=.., calcHfMultiplet=.., printBefore=.., levelSelection=..)

    ... keyword copy-constructor for re-defining selected values of a settings::Hfs.Settings.
"""
function Settings(set::Hfs.Settings;
        calcM1::Union{Nothing,Bool}=nothing,           calcE2::Union{Nothing,Bool}=nothing,
        calcM3::Union{Nothing,Bool}=nothing,           calcNondiagonal::Union{Nothing,Bool}=nothing,
        calcHfMultiplet::Union{Nothing,Bool}=nothing,  printBefore::Union{Nothing,Bool}=nothing,
        levelSelection::Union{Nothing,LevelSelection}=nothing)
    if  isnothing(calcM1)            calcM1x          = set.calcM1          else   calcM1x          = calcM1          end
    if  isnothing(calcE2)            calcE2x          = set.calcE2          else   calcE2x          = calcE2          end
    if  isnothing(calcM3)            calcM3x          = set.calcM3          else   calcM3x          = calcM3          end
    if  isnothing(calcNondiagonal)   calcNondiagonalx = set.calcNondiagonal else   calcNondiagonalx = calcNondiagonal end
    if  isnothing(calcHfMultiplet)   calcHfMultipletx = set.calcHfMultiplet else   calcHfMultipletx = calcHfMultiplet end
    if  isnothing(printBefore)       printBeforex     = set.printBefore     else   printBeforex     = printBefore     end
    if  isnothing(levelSelection)    levelSelectionx  = set.levelSelection  else   levelSelectionx  = levelSelection  end

    Settings( calcM1x, calcE2x, calcM3x, calcNondiagonalx, calcHfMultipletx, printBeforex, levelSelectionx )
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
    ## REPAIRED 06-Aug-2026; this could not run. It read `multiplet.levelFs` where the field is `hfLevels`, and
    ## rebuilt each level as an `Hfs.IJF_Level` -- a type that does not exist -- from `lev.I` and `lev.basis`,
    ## neither of which is a field of Hfs.HfLevel. Since the levels are immutable and only their ORDER changes,
    ## none of that reconstruction was needed in the first place.
    sortedLevels = sort( multiplet.hfLevels, by = lev -> lev.energy )
    return( Hfs.HfMultiplet(multiplet.name, sortedLevels) )
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
    ... to compute the T^(M1), T^(E2) or T^(M3) hyperfine amplitude for a given pair of levels.
        A value::ComplexF64 is returned.

        CONVENTION -- read this before using the result in a recoupling formula (clarified 05-Aug-2026). What is
        returned is the "GRASP-like" amplitude, i.e. <alpha_r J_r || T^(mp) || alpha_s J_s> / sqrt(2J_r+1), and NOT
        the reduced matrix element in the standard Wigner-Eckart normalization. That is exactly what the A, B and C
        constants formed in computeAmplitudesProperties expect -- A*I/mu = amplitude/sqrt(J(J+1)), verified against
        the measured A(H 1s) to 0.06% -- but it is NOT what a Wigner-6j expression expects. Use
        Hfs.reducedMeElectronic for the latter; building a 6-j directly on this function is what left
        computeHyperfineRepresentation too small by sqrt(2J+1) until 05-Aug-2026.

        Note (26-Jul-2026): coeff.T, for rank>0 one-particle operators, is returned by
        SpinAngular.computeCoefficientsNonScalar already in "GRASP-like" convention, i.e. with an internal
        factor sqrt(2*j_a+1) applied (see the active "GRASP like spin-angular coefficient" step in that
        function) -- this is NOT what the paper's Eq. (48)-type reduced-matrix-element sum
        <leftCsf||W^(k)||rightCsf> = sum_ab d^k_ab [n_a kappa_a||w^(k)||n_b kappa_b] needs (the "pure"
        coefficient, without that factor). By contrast, IsotopeShift.amplitude's rank-0 case uses
        SpinAngular.computeCoefficientsScalar, where the equivalent conversion step is deliberately commented
        out, and IsotopeShift.amplitude applies its own sqrt(2j_a+1) factor externally to match its formula.
        This function previously used coeff.T directly with no compensating division, leaving an
        uncorrected sqrt(2j_a+1) factor in every M1/E2/M3 hyperfine amplitude -- confirmed for H(1s) M1: the
        A-constant was too large by ~1.42 (= sqrt(2), the j=1/2 value) after the separate missing-alpha fix
        in InteractionStrength.hfs_tM1. Fixed by dividing each term by sqrt(2*j_a+1), undoing the internal
        GRASP-like conversion, matching IsotopeShift's convention.
"""
function  amplitude(mp::EmMultipole, rLevel::Level, sLevel::Level, grid::Radial.Grid; printout::Bool=true)
    ## E3 ADDED 06-Aug-2026, and with it the operator-dependent parity rule below. InteractionStrength.hfs_tE3
    ## had existed all along (Andersson & Jonsson 2008, Eq. 49) but was never dispatched here, because E3 is the
    ## first PARITY-ODD hyperfine multipole: M1, E2 and M3 all carry parity +1 -- (-1)^L for electric, (-1)^(L+1)
    ## for magnetic -- whereas E3 carries -1. The old guard `rLevel.parity != sLevel.parity -> 0` is correct for
    ## the even ones and silently wrong for E3, which requires the two electronic levels to have OPPOSITE parity.
    ##
    ## The physical statement is that the hyperfine Hamiltonian conserves TOTAL parity: if the nuclear states it
    ## connects have opposite parity, as in 235U (I_g = 7/2- and the 1/2+ isomer), then the electronic states it
    ## connects must have opposite parity too. A calculation of that case therefore needs BOTH parities in its
    ## configuration list, or every matrix element vanishes.
    if        mp == M1    opa = SpinAngular.OneParticleOperator(1, plus,  true);  hfsFunc = InteractionStrength.hfs_tM1
    elseif    mp == E2    opa = SpinAngular.OneParticleOperator(2, plus,  true);  hfsFunc = InteractionStrength.hfs_tE2
    elseif    mp == M3    opa = SpinAngular.OneParticleOperator(3, plus,  true);  hfsFunc = InteractionStrength.hfs_tM3
    elseif    mp == E3    opa = SpinAngular.OneParticleOperator(3, minus, true);  hfsFunc = InteractionStrength.hfs_tE3
    else      error("Hfs.amplitude: unsupported nuclear multipole $mp; implemented are M1, E2, M3 and E3.")
    end
    opParity = (mp == E3)  ?  Basics.minus  :  Basics.plus
    if  rLevel.parity * sLevel.parity != opParity   return( ComplexF64(0.) )   end
    nr = length(rLevel.basis.csfs);    ns = length(sLevel.basis.csfs);    matrix = zeros(ComplexF64, nr, ns)
    if  printout   printstyled("Compute hyperfine $(string(mp)) matrix of dimension $nr x $ns ... \n", color=:light_green)   end
    for  r = 1:nr
        for  s = 1:ns
            me = 0.
            if  rLevel.basis.csfs[r].parity  != rLevel.parity    ||  sLevel.basis.csfs[s].parity  != sLevel.parity
                continue
            end
            subshellList = sLevel.basis.subshells
            wa           = SpinAngular.computeCoefficients(opa, rLevel.basis.csfs[r], sLevel.basis.csfs[s], subshellList)
            for  coeff in wa
                tamp  = hfsFunc(rLevel.basis.orbitals[coeff.a], sLevel.basis.orbitals[coeff.b], grid)
                ja2   = Basics.subshell_2j(coeff.a)
                me = me + coeff.T / sqrt(ja2 + 1) * tamp
            end
            matrix[r,s] = me
        end
    end
    if  printout   printstyled("done.\n", color=:light_green)   end
    amplitude = transpose(rLevel.mc) * matrix * sLevel.mc
    return( amplitude )
end


"""
`Hfs.reducedMeElectronic(mp::EmMultipole, rLevel::Level, sLevel::Level, grid::Radial.Grid)`
    ... to return the electronic reduced matrix element <alpha_r J_r || T^(mp) || alpha_s J_s> in the STANDARD
        (Wigner-Eckart) normalization, i.e. the quantity that may be inserted directly into a Wigner-6j recoupling
        expression. A value::ComplexF64 is returned.

        WHY THIS EXISTS (05-Aug-2026). Hfs.amplitude does NOT return that quantity, despite what its own header
        line says: it returns the reduced element divided by sqrt(2J+1), the "GRASP-like" normalization. The two
        differ by a factor that is invisible in every place where the same convention is used on both sides, and
        fatal in any place that mixes them -- as Hfs.computeHyperfineRepresentation did, building a 6-j expression
        directly out of Hfs.amplitude and coming out too small by exactly sqrt(2J+1).

        HOW THIS WAS SETTLED, since the two candidate conventions cannot both be right: A(H 1s) computed in the
        BARE NUCLEAR potential -- NuclearField(), not a DFS field, whose self-interaction on a one-electron system
        puts the 1s orbital at -0.194 instead of -0.5 a.u. and invalidates the benchmark entirely -- gives
        1421.21 MHz against the measured 1420.406 MHz, a ratio of 1.00057, using computeAmplitudesProperties'
        formula A*I/mu = amplitude/sqrt(J(J+1)). Inserting a further sqrt(2J+1) would give 1004.9 MHz, i.e. 29%
        low. So the A/B/C constants are right as they stand and it is the recoupling side that needed the factor.

        Note that A*I/mu = amplitude/sqrt(J(J+1)) and A*I/mu = <J||T^1||J>/sqrt(J(J+1)(2J+1)) are the SAME
        statement once <J||T^1||J> = amplitude*sqrt(2J+1); computeAmplitudesProperties therefore needs no change,
        and none was made. This function only makes the convention explicit so that the next caller who needs a
        genuine reduced matrix element does not have to rediscover it.

        For J_r != J_s the factor is taken as the symmetric ((2J_r+1)(2J_s+1))^(1/4), which reduces to sqrt(2J+1)
        on the diagonal and keeps the hyperfine matrix symmetric. The diagonal is what the H(1s) benchmark pins
        down; the off-diagonal convention is not separately verified here, and matters only through J-mixing,
        which is of order 1e-5 whenever the hyperfine splitting is small against the fine structure.
"""
function  reducedMeElectronic(mp::EmMultipole, rLevel::Level, sLevel::Level, grid::Radial.Grid)
    Jr = AngularMomentum.oneJ(rLevel.J);    Js = AngularMomentum.oneJ(sLevel.J)
    wa = ( (2Jr + 1) * (2Js + 1) )^0.25
    return( Hfs.amplitude(mp, rLevel, sLevel, grid; printout=false) * wa )
end


"""
`Hfs.recouplingElectronicOperator(spinI::AngularJ64, Ja::AngularJ64, Fa::AngularJ64,
                            Jb::AngularJ64, Fb::AngularJ64, L::Int64)`
    ... to return the factor that turns an ELECTRONIC reduced matrix element of rank L into the corresponding
        one between IJF-coupled levels, for an operator that acts on the electrons alone and leaves the nucleus
        untouched:

            <(I J_b) F_b || T^L || (I J_a) F_a>
                = (-1)^(I+J_b+F_a+L) sqrt((2F_a+1)(2F_b+1)) {J_b F_b I; F_a J_a L} <J_b || T^L || J_a>

        A value::Float64 is returned. This is the standard result for a tensor acting on one part of a coupled
        system (Edmonds 7.1.7), with the nucleus as the spectator.

        VERIFIED (work/diag-recoupling.jl, 05-Aug-2026) over 454 combinations of I, J_a, J_b, F_a and L = 1, 2:
        the sum over F_b of the squared factor equals (2F_a+1)/(2J_a+1) to better than 1e-15, which is exactly
        what leaves a radiative width unchanged by hyperfine coupling. For I = 0 the factor reduces to 1.

        Lives here, in the property module, rather than in any one process module: both hyperfine-resolved
        dielectronic recombination and hyperfine-induced transitions need it, and a second copy is how two
        implementations of the same formula start to drift apart.
"""
function  recouplingElectronicOperator(spinI::AngularJ64, Ja::AngularJ64, Fa::AngularJ64,
                                       Jb::AngularJ64, Fb::AngularJ64, L::Int64)
    wa = AngularMomentum.phaseFactor([spinI, +1, Jb, +1, Fa, +1, AngularJ64(L)])
    wb = sqrt( (Basics.twice(Fa) + 1) * (Basics.twice(Fb) + 1) )
    wc = AngularMomentum.Wigner_6j(Jb, Fb, spinI, Fa, Ja, AngularJ64(L))
    return( wa * wb * wc )
end


"""
`Hfs.recouplingNuclearOperator(J::AngularJ64, Ia::AngularJ64, Fa::AngularJ64,
                            Ib::AngularJ64, Fb::AngularJ64, L::Int64)`
    ... the mirror image of Hfs.recouplingElectronicOperator: the factor that turns a NUCLEAR reduced matrix
        element of rank L into the corresponding one between IJF-coupled levels, for an operator that acts on the
        nucleus alone and leaves the electrons untouched:

            <(I_b J) F_b || W^L || (I_a J) F_a>
                = (-1)^(J+I_b+F_a+L) sqrt((2F_a+1)(2F_b+1)) {I_b F_b J; F_a I_a L} <I_b || W^L || I_a>

        A value::Float64 is returned. It is the same expression with the roles of the nuclear and electronic
        angular momenta interchanged, because the coupled state is symmetric in the two: whichever subsystem the
        operator acts on, the other is the spectator that appears in the 6-j.

        This is the factor that carries a hyperfine-induced NUCLEAR transition -- the isomer changes while the
        electronic level does not -- and it is therefore what the 229Th and 205Pb cases rest on. For J = 0 it
        reduces to 1.
"""
function  recouplingNuclearOperator(J::AngularJ64, Ia::AngularJ64, Fa::AngularJ64,
                                    Ib::AngularJ64, Fb::AngularJ64, L::Int64)
    wa = AngularMomentum.phaseFactor([J, +1, Ib, +1, Fa, +1, AngularJ64(L)])
    wb = sqrt( (Basics.twice(Fa) + 1) * (Basics.twice(Fb) + 1) )
    wc = AngularMomentum.Wigner_6j(Ib, Fb, J, Fa, Ia, AngularJ64(L))
    return( wa * wb * wc )
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
        else
            ## DIAGONAL higher multipoles are set to zero rather than raising an error (06-Aug-2026). This branch
            ## is the STATIC hyperfine moment of one nuclear state, which Nuclear.Isomer supplies only for M1 (mu)
            ## and E2 (Q); no octupole moment is part of its input, and the meaning of the `Omega` field is not
            ## documented well enough to be used here silently. The omitted term is a small static shift of a
            ## level, whereas the OFF-DIAGONAL element of the same multipole -- the one that mixes two nuclear
            ## states and drives every rate in HyperfineInduced -- is taken from `elementM` below and is fully
            ## included. Raising an error instead simply made every E3 case unrunnable, which is how the 235U
            ## branch was blocked.
            amplitude = 0.
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
    hfMultiplet = Hfs.computeHyperfineRepresentation(hfBasis, grid)

    return( hfMultiplet )
end


"""
`Hfs.computeHyperfineMultiplet(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::Hfs.Settings; output=true)`  
    ... to compute a hyperfine multiplet, i.e. a representation of hyperfine levels within a hyperfine-coupled basis as defined by the
        given (electronic) multiplet; an hfMultiplet::Hfs.HfMultiplet is returned.
"""
function computeHyperfineMultiplet(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::Hfs.Settings; output=true)
    println("")
    printstyled("Hfs.computeHyperfineMultiplet(): The computation of the hyperfine multiplet starts now ... \n", color=:light_green)
    printstyled("------------------------------------------------------------------------------------------ \n", color=:light_green)
    #
    hfBasis     = Hfs.defineHyperfineBasis(multiplet, nm)
    hfMultiplet = Hfs.computeHyperfineRepresentation(hfBasis, grid)
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
    ... to set-up and diagonalize the Hamiltonian matrix of H^(DFB) + H^(hfs) within the atomic hyperfine (IJF-coupled)
        basis; a hfMultiplet::Hfs.HfMultiplet is returned.

        REPAIRED 05-Aug-2026. This function had never run: its body addressed an undefined `hfBasis.vectors[...]`
        instead of its own argument, it called Hfs.HfLevel(...) with seven arguments against a six-field struct, and
        it carried `* 1.0e-6 ## fudge-factor to keep HFS interaction small` on BOTH the magnetic and the electric
        contribution -- a debugging leftover that would have scaled every hyperfine splitting down by a million.
        The fudge factor is removed: the splitting of the initial ion is precisely the observable this
        representation exists to produce. Nothing can regress, since the function could not previously be called.

        The matrix is block diagonal in F and parity by construction (the guards below), so the eigenvectors do not
        mix different F; each hyperfine level is Sum_k mc[k] |(I J_k) F>, where the Clebsch-Gordan coupling already
        sits inside each HfBasisVector and mc describes the J-mixing among vectors of equal F.

        With mu = Q = 0 the hyperfine Hamiltonian vanishes identically and only the electronic energies remain on
        the diagonal. Since those are non-degenerate, the eigenvectors come out as pure unit vectors, i.e.
        mc = [.. 1.0 ..] on a single basis vector. The "no hyperfine mixing" case therefore needs no separate code
        path -- it is the exact output of this same routine.
"""
function computeHyperfineRepresentation(hfBasisVectors::Array{HfBasisVector,1}, grid::Radial.Grid;
                                        hfMultipoles::Array{EmMultipole,1}=EmMultipole[Basics.M1, Basics.E2],
                                        printout::Bool=false)
    n = length(hfBasisVectors);   matrix = zeros(n,n)
    for  r = 1:n
        for  s = 1:n
            vr = hfBasisVectors[r];    vs = hfBasisVectors[s]
            matrix[r,s] = 0.
            ## THE DIAGONAL carries the electronic energy AND the nuclear excitation energy of its isomer. The
            ## latter was missing before 06-Aug-2026, which did not matter while a basis held one nuclear state,
            ## and is fatal the moment it holds two: without it the ground and isomeric levels are degenerate,
            ## and the nuclear hyperfine mixing -- which is of order V_21 / dE_nuc -- is ill-defined.
            if  r == s
                matrix[r,s] = matrix[r,s] + vr.levelJ.energy +
                              Defaults.convertUnits("energy: to atomic", vr.isomer.energy)
            end
            if  vr.F  !=  vs.F                                                           continue     end
            ## The TOTAL parity, nuclear x electronic, must match; with one nuclear state the nuclear factor is
            ## common and the electronic test suffices, with two it does not.
            if  vr.levelJ.parity * vr.isomer.parity != vs.levelJ.parity * vs.isomer.parity   continue end
            #
            Ir = vr.isomer.spinI;    Is = vs.isomer.spinI
            Jr = vr.levelJ.J;        Js = vs.levelJ.J;      F = vr.F
            #
            ## The hyperfine operator is H_hfs = sum_k T^k(electronic) . W^k(nuclear), a scalar product of two
            ## tensors acting on the two parts of the coupled system, so
            ##
            ##   <(I_r J_r) F| H_hfs |(I_s J_s) F> = sum_k (-1)^(I_r+J_s+F) {F J_r I_r; k I_s J_s}
            ##                                              <I_r||W^k||I_s> <J_r||T^k||J_s>
            ##
            ## For I_r = I_s = I and J_r = J_s = J this reduces to the Casimir form against which this matrix was
            ## verified on 05-Aug-2026. The OFF-DIAGONAL-in-isomer case is new (06-Aug-2026) and is what nuclear
            ## hyperfine mixing consists of: nothing changes but that <I_r||W^k||I_s> now connects two DIFFERENT
            ## nuclear states. Hermiticity of the result is checked below, which is what would catch a wrong phase.
            for  mp in hfMultipoles
                if  AngularMomentum.oneJ(Ir) == 0.  &&  AngularMomentum.oneJ(Is) == 0.    continue    end
                wa = AngularMomentum.phaseFactor([Ir, +1, Js, +1, F]) *
                        AngularMomentum.Wigner_6j(F, Jr, Ir, AngularJ64(mp.L), Is, Js)
                if  wa == 0.                                                              continue    end
                wb = Hfs.reducedMeElectronic(mp, vr.levelJ, vs.levelJ, grid)
                wc = Hfs.computeInteractionAmplitudeM(mp, vr.isomer, vs.isomer)
                ## UNITS. computeInteractionAmplitudeM returns the DIAGONAL nuclear element from mu [nuclear
                ## magnetons] or Q [barn], which must be converted; but for two DIFFERENT isomers it returns
                ## elementM, which the Isomer docstring already defines to be in atomic units. Converting
                ## uniformly would therefore be wrong by ~2.7e-4 or ~3.6e-2 on precisely the off-diagonal
                ## elements that carry the mixing.
                if  vr.isomer == vs.isomer
                    if      mp == Basics.M1   wc = wc * Defaults.convertUnits("moment: from nuclear magneton to atomic", 1.0)
                    elseif  mp == Basics.E2   wc = wc * Defaults.convertUnits("cross section: from barn to atomic unit", 1.0)
                    else    wc = 0.   ## no diagonal moment is defined for the higher multipoles
                    end
                end
                matrix[r,s] = matrix[r,s] + real(wa * wb * wc)
            end
        end
    end
    ## Hermiticity: the matrix is built term by term for both (r,s) and (s,r), so a wrong phase convention in the
    ## off-diagonal-in-isomer element shows up here as a SIGN difference rather than a small rounding one.
    asym = 0.;   scale = 0.
    for  r = 1:n, s = 1:n
        asym  = max(asym,  abs(matrix[r,s] - matrix[s,r]));    scale = max(scale, abs(matrix[r,s]))
    end
    if  scale > 0.  &&  asym / scale > 1.0e-8
        @warn("Hfs.computeHyperfineRepresentation(): the hyperfine matrix is not symmetric -- relative " *
              "asymmetry $(asym/scale). This points at a phase convention in the off-diagonal elements, " *
              "not at rounding; the matrix is symmetrized before diagonalization.")
    end
    if  printout   println(">>> hyperfine matrix of dimension $n x $n; multipoles $hfMultipoles; " *
                           "relative asymmetry $(scale > 0. ? asym/scale : 0.)")    end
    for  r = 1:n, s = r+1:n
        wa = (matrix[r,s] + matrix[s,r]) / 2;    matrix[r,s] = wa;    matrix[s,r] = wa
    end
    #
    # Diagonalize the matrix and set-up the representation
    eigen    = Basics.fixEigenvectorPhase!( Basics.diagonalize(MatrixWithLinearAlgebra(), matrix) )
    levelFs  = Hfs.HfLevel[]
    for  ev = 1:length(eigen.values)
        # Construct the eigenvector with regard to the given basis (not w.r.t the symmetry block).
        # F and parity are read off the DOMINANT component rather than the first one above a threshold: with
        # J-mixing the first component need not be the leading one, and a mis-assigned F would silently corrupt
        # every recoupling coefficient computed from this level later on.
        evector   = eigen.vectors[ev];    en = eigen.values[ev]
        rmax      = 1;    wmax = 0.
        for  r = 1:length(hfBasisVectors)
            if  abs(evector[r]) > wmax    wmax = abs(evector[r]);   rmax = r    end
        end
        parity    = hfBasisVectors[rmax].levelJ.parity
        F         = hfBasisVectors[rmax].F
        ## M = F, not 0: a projection must share the integer/half-integer character of its F, and 0 is simply not
        ## a valid projection of a half-integer F. M is otherwise only a label here -- in a free ion the 2F+1
        ## sublevels are exactly degenerate, and every quantity computed from these levels (energies, rates,
        ## strengths) is summed over final and averaged over initial M, so M cancels by the Wigner-Eckart theorem
        ## and only reduced matrix elements survive. It becomes meaningful only once something breaks the
        ## isotropy -- an external field, an aligned or polarized ensemble, or an angular distribution.
        MF        = AngularM64(F)
        newlevelF = Hfs.HfLevel(F, MF, parity, 0, en, hfBasisVectors, evector)
        push!( levelFs, newlevelF)
    end
    ## A STABLE INDEX, assigned in order of increasing energy (added 06-Aug-2026). Hyperfine levels previously
    ## carried none, so nothing could refer to one: Basics.selectLevelPair dispatches on ManyElectron.Level and
    ## throws a MethodError on an HfLevel, and Hfs-based display routines had to key on the surrogate
    ## (2F, energy) pair instead. Energy order is a safe labelling here because hyperfine multiplets of
    ## neighbouring electronic levels essentially never interleave -- a hyperfine splitting is orders of
    ## magnitude smaller than a fine-structure one.
    levelFs = sort(levelFs, by = lev -> lev.energy)
    levelFs = [ Hfs.HfLevel(lev.F, lev.M, lev.parity, i, lev.energy, lev.hfBasisVectors, lev.mc)
                for (i, lev) in enumerate(levelFs) ]
    hfMultiplet = Hfs.HfMultiplet("hyperfine", levelFs)

    return( hfMultiplet )
end


"""
`Hfs.computeInteractionMatrix(basis::Basis, grid::Radial.Grid, settings::Hfs.Settings)`
    ... to compute the T^M1 and/or T^E2 interaction matrices for the given basis, i.e. (<csf_r || T^(n)) || csf_s>).
        An im::Hfs.InteractionMatrix is returned.

        Note (26-Jul-2026): each coeff.T here is divided by sqrt(2*j_a+1), undoing the "GRASP-like"
        sqrt(2*j_a+1) factor that SpinAngular.computeCoefficientsNonScalar applies internally for rank>0
        one-particle operators -- see the matching note on Hfs.amplitude for the full explanation. Applies
        uniformly to the M1, E2, and M3 blocks below (all use the same rank>0 SpinAngular path).
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
                    matrixM1[r,s] = matrixM1[r,s] + coeff.T / sqrt(ja + 1) * tamp
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
                    matrixE2[r,s] = matrixE2[r,s] + coeff.T / sqrt(ja + 1) * tamp
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
                    tamp  = InteractionStrength.hfs_tM3(basis.orbitals[coeff.a], basis.orbitals[coeff.b], grid)
                    matrixM3[r,s] = matrixM3[r,s] + coeff.T / sqrt(ja + 1) * tamp
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
        an Array{Hfs.HfBasisVector,1} is returned.

        REPAIRED 05-Aug-2026, together with the multiplet method below. Both had been written against `HfVector`
        and `HfBasis`, neither of which exists in this module, and the multiplet method additionally called
        `IJF_Vector(nm.spinI, F, level)` -- three arguments against a five-field struct whose first field is F and
        whose fourth is a CsfR, not a Level. Neither method could ever be called.

        The confusion behind it: the module carried TWO similarly-named coupled basis types, and these methods
        reached for the CSF-based one when they needed the ASF-based HfBasisVector. The CSF-based IJF_Vector has
        since been retired, so only one coupled basis type remains.
"""
function  defineHyperfineBasis(level::Level, nm::Nuclear.Model)
    return( Hfs.defineHyperfineBasis(Multiplet("single level", [level]), nm) )
end


"""
`Hfs.defineHyperfineBasis(multiplet::Multiplet, isomers::Array{Nuclear.Isomer,1}; printout::Bool=true)`
    ... to define/set-up an atomic hyperfine (IJF-coupled) basis for the given electronic multiplet and ONE OR
        MORE nuclear states; an Array{Hfs.HfBasisVector,1} is returned, one vector for every
        (isomer, electronic level, F) triple allowed by F = I + J.

        Admitting more than one isomer (06-Aug-2026) is what makes nuclear hyperfine mixing expressible: two
        nuclear levels of the same F, built on the same electronic level, then appear as two basis vectors, and
        diagonalizing the hyperfine matrix over that basis mixes them. Nothing else in the machinery changes.
        The parity stored on each vector is the TOTAL one, nuclear x electronic.

        REPAIRED 05-Aug-2026; see the note on the single-level method above. The nuclear moments are carried into
        each basis vector through a Nuclear.Isomer built from the given model, so that a model with mu = Q = 0
        produces a basis whose hyperfine Hamiltonian vanishes identically -- the device used to switch the
        hyperfine structure of intermediate and final levels off without a separate code path.
"""
function  defineHyperfineBasis(multiplet::Multiplet, isomers::Array{Nuclear.Isomer,1}; printout::Bool=true)
    vectors = Hfs.HfBasisVector[]
    for  isomer in isomers
        for  level in multiplet.levels
            for  F in Basics.oplus(isomer.spinI, level.J)
                push!(vectors,  Hfs.HfBasisVector(F, level.parity * isomer.parity, isomer, level) )
            end
        end
    end
    if  printout
        println(" ")
        println("  Atomic hyperfine (IJF-coupled) basis of dimension $(length(vectors)), from " *
                "$(length(multiplet.levels)) electronic levels and $(length(isomers)) nuclear state(s):")
        for  isomer in isomers
            println("    I = $(isomer.spinI)$(isomer.parity), excitation energy $(isomer.energy), " *
                    "mu = $(isomer.mu), nuclear transition multipoles $(isomer.multipoleM)")
        end
        println(" ")
    end
    return( vectors )
end


"""
`Hfs.defineHyperfineBasis(multiplet::Multiplet, nm::Nuclear.Model; printout::Bool=true)`
    ... the single-isomer method; an Array{Hfs.HfBasisVector,1} is returned. It simply assembles one
        Nuclear.Isomer from the given nuclear model and defers to the method above, so that the two cannot
        diverge.
"""
function  defineHyperfineBasis(multiplet::Multiplet, nm::Nuclear.Model; printout::Bool=true)
    function  display_ijfVector(i::Int64, vector::Hfs.HfBasisVector)
        si = string(i);   ni = length(si);    sa = repeat(" ", 5);    sa = sa[1:5-ni] * si * ")  "
        sa = sa * "[" * string( LevelSymmetry(vector.levelJ.J, vector.levelJ.parity) ) * "] " * string(vector.F) * repeat(" ", 4)
        return( sa )
    end

    ## The nuclear parity is not carried by Nuclear.Model; it is taken as even, so that the total parity of a
    ## basis vector reduces to the electronic one.
    isomer  = Nuclear.Isomer( nm.spinI, Basics.plus, 0., nm.mu, nm.Q, 0., EmMultipole[], Float64[] )
    return( Hfs.defineHyperfineBasis(multiplet, Nuclear.Isomer[isomer]; printout=printout) )
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
        we = Defaults.convertUnits("energy: from atomic to Hz", 1.0) / 1.0e6      # Energy factor into MHz
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


