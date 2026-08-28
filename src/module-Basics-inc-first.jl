

#################################################################################################################################
#################################################################################################################################


"""
`@enum   Basics.EmGauge`  
    ... defines a enumeration for the (allowed) gauges to deal with the radiation field; this differs from the UseGauge which can 
        be selected for explicit computations.
                    
        + NoGauge      ... No gauge define (yet).
        + Coulomb      ... Coulomb gauge (= velocity gauge in the non-relativistic limit; Grant, 1974).
        + Babushkin    ... Babushkin gauge (= length gauge in the non-relativistic limit; Grant, 1974). 
        + Magnetic     ... Magnetic gauge (transverse velocity gauge).
        + Velocity     ... Transverse velocity gauge (Johnson, 2007).
        + Length       ... Length gauge (Johnson, 2007).
"""
@enum   EmGauge    NoGauge    Coulomb   Babushkin   Magnetic   Velocity   Length

@doc "... No gauge define (yet)."                                                                                                              NoGauge
@doc "... Coulomb gauge (= velocity gauge in the non-relativistic limit; Grant, 1974)."                                                        Coulomb
@doc "... Babushkin gauge (= length gauge in the non-relativistic limit; Grant, 1974)."                                                      Babushkin
@doc "... Magnetic gauge (transverse velocity gauge)."                                                                                        Magnetic
@doc "... Transverse velocity gauge (Johnson, 2007)."                                                                                         Velocity
@doc "... Length gauge (Johnson, 2007)."                                                                                                        Length
export  EmGauge,   NoGauge,   Coulomb,  Babushkin,  Magnetic,  Velocity,  Length


"""
`Basics.EmGauge(sa::String)`  ... constructor for a given String.
"""
function EmGauge(sa::String)
    if       sa in ["none", "NoGauge"]         wa = NoGauge
    elseif   sa == "Coulomb"                   wa = Coulomb
    elseif   sa == "Babushkin"                 wa = Babushkin
    elseif   sa in ["Magnetic", "magnetic"]    wa = Magnetic
    else     error("stop a")
    end

    return( wa )
end


# `Base.string(gauge::EmGauge)`  ... provides a proper printout of the variable gauge::EmGauge.
function Base.string(gauge::EmGauge) 
    if      gauge == NoGauge     return( "no gauge" )
    elseif  gauge == Coulomb     return( "Coulomb" )   
    elseif  gauge == Babushkin   return( "Babushkin" )   
    elseif  gauge == Magnetic    return( "Magnetic" )   
    else    error("stop a")
    end
end

#################################################################################################################################
#################################################################################################################################


"""
`struct  Basics.EmMultipole`  ... defines a struct for the em multipoles, E1, M1, E2, ...

    + L           ::Int64    ... multipolarity
    + electric    ::Bool     ... True if the multipole is an electric one, and false otherwise.
"""  
struct  EmMultipole
    L             ::Int64
    electric      ::Bool   
end

export  EmMultipole, E1, M1, E2, M2, E3, M3, E4, M4, multipoleParity


const  E1 = EmMultipole(1, true);     const  M1 = EmMultipole(1, false)
const  E2 = EmMultipole(2, true);     const  M2 = EmMultipole(2, false)
const  E3 = EmMultipole(3, true);     const  M3 = EmMultipole(3, false)
const  E4 = EmMultipole(4, true);     const  M4 = EmMultipole(4, false)

@doc "... the electric dipole (E1) multipole of the radiation field."                                                                               E1
@doc "... the magnetic dipole (M1) multipole of the radiation field."                                                                               M1
@doc "... the electric quadrupole (E2) multipole of the radiation field."                                                                           E2
@doc "... the magnetic quadrupole (M2) multipole of the radiation field."                                                                           M2
@doc "... the electric octupole (E3) multipole of the radiation field."                                                                             E3
@doc "... the magnetic octupole (M3) multipole of the radiation field."                                                                             M3
@doc "... the electric hexadecapole (E4) multipole of the radiation field."                                                                         E4
@doc "... the magnetic hexadecapole (M4) multipole of the radiation field."                                                                         M4


"""
`Basics.EmMultipole(sa::String)`  ... constructor for a given String.
"""
function EmMultipole(sa::String)
    if       sa == "E1"    L = 1;   electric = true
    elseif   sa == "M1"    L = 1;   electric = false
    elseif   sa == "E2"    L = 2;   electric = true
    elseif   sa == "M2"    L = 2;   electric = false
    elseif   sa == "E3"    L = 3;   electric = true
    elseif   sa == "M3"    L = 3;   electric = false
    elseif   sa == "E4"    L = 4;   electric = true
    elseif   sa == "M4"    L = 4;   electric = false
    else     error("stop a")
    end

    EmMultipole(L, electric)
end


# `Base.show(io::IO, mp::EmMultipole)`  ... prepares a proper printout of the variable mp::EmMultipole.
function Base.show(io::IO, mp::EmMultipole) 
    if  mp.electric   sa = "E"*string(mp.L)    else    sa = "M"*string(mp.L)   end
    print(io, string(mp) )
end


# `Base.string(mp::EmMultipole)`  ... provides a proper printout of the variable mp::EmMultipole.
function Base.string(mp::EmMultipole)
    mp.L < 1    &&    error("Improper multipolarity L = $(mp.L)") 
    if  mp.electric   sa = "E"*string(mp.L)    else    sa = "M"*string(mp.L)   end
    return( sa )
end


"""
`Basics.multipoleParity(mp::EmMultipole)`  
    ... returns the PARITY of the one-particle tensor operator that belongs to the given multipole: an electric
        2^L-pole carries parity (-1)^L and a magnetic one (-1)^(L+1), so E1 and E3 are odd while M1, E2 and M3 are
        even. This is the parity that a spin-angular operator must be constructed with, and it decides which
        coefficients survive the parity selection rule -- an even operator connects orbitals with l_a + l_b even,
        an odd one those with l_a + l_b odd. A value::Basics.Parity is returned.

        IT EXISTS BECAUSE HARD-CODING `plus` IS WRONG AND USED TO BE HARMLESS. `SpinAngularGaigalas.computeCoefficients`
        ignores the parity it is given and returns both sets, leaving the radial integral to kill the wrong half;
        `SpinAngular` applies the rule, so a caller that passes `plus` for an E1 operator now silently loses
        exactly the coefficients E1 needs. Measured on 2s <-> 2p at rank 1: SpinAngularGaigalas returns 9 coefficients
        for either parity, SpinAngular returns 5 for plus and 4 for minus.
"""
function multipoleParity(mp::EmMultipole)
    return( iseven(mp.L) == mp.electric  ?  Basics.plus  :  Basics.minus )
end


"""
`Basics.multipole_p(mp::EmMultipole)`  ... returns the kind integer p=0 (magnetic) or p=1 (electric) of the given  mp::EmMultipole.
"""
function multipole_p(mp::EmMultipole)
    if  mp.electric   return(1)    else    return(0)   end
end

#################################################################################################################################
#################################################################################################################################


"""
`struct  EmProperty`  ... defines a type to maintain two gauge forms of a computed result that depends on the radiation field.

    + Coulomb         ::Float64    ... Value for the Coulomb gauge of the radiation field.
    + Babushkin       ::Float64    ... Value for the Coulomb gauge of the radiation field.
"""
struct  EmProperty 
    Coulomb           ::Float64
    Babushkin         ::Float64
end

export EmProperty


"""
`Basics.EmProperty(wa::Float64)`  ... constructor for an `constant` instance of EmProperty.
"""
function EmProperty(wa::Float64)
    EmProperty(wa, wa)
end


Base.:+(a::EmProperty, b::EmProperty) = EmProperty(a.Coulomb + b.Coulomb, a.Babushkin + b.Babushkin)
Base.:+(a::EmProperty, b) = EmProperty(a.Coulomb + b, a.Babushkin + b)
Base.:+(a, b::EmProperty) = b + a
Base.:-(a::EmProperty, b::EmProperty) = EmProperty(a.Coulomb - b.Coulomb, a.Babushkin - b.Babushkin)
Base.:*(a::EmProperty, b::EmProperty) = EmProperty(a.Coulomb * b.Coulomb, a.Babushkin * b.Babushkin)
Base.:*(a, b::EmProperty) = EmProperty(a * b.Coulomb, a * b.Babushkin)
Base.:*(a::EmProperty, b) = EmProperty(a.Coulomb * b, a.Babushkin * b)
Base.:/(a::EmProperty, b::EmProperty) = EmProperty(a.Coulomb / b.Coulomb, a.Babushkin / b.Babushkin)
Base.:/(a, b::EmProperty) = EmProperty(a / b.Coulomb, a / b.Babushkin)


# `Base.show(io::IO, property::EmProperty)`  ... prepares a proper printout of the variable property::EmProperty.
function Base.show(io::IO, property::EmProperty) 
    sa = Base.string(property);                print(io, sa)
end


# `Base.string(property::EmProperty)`  ... provides a String notation for the variable property::EmProperty.
function Base.string(property::EmProperty)
    sa = "$(property.Coulomb) [Coulomb],  $(property.Babushkin) [Babushkin]"
    return( sa )
end


"""
`struct  EmPropertyC`  ... defines a type to maintain two gauge forms of a computed result that depends on the radiation field.

    + Coulomb         ::Complex{Float64}    ... (Complex) Value for the Coulomb gauge of the radiation field.
    + Babushkin       ::Complex{Float64}    ... (Complex) Value for the Coulomb gauge of the radiation field.
"""
struct  EmPropertyC 
    Coulomb           ::Complex{Float64}
    Babushkin         ::Complex{Float64}
end

export EmPropertyC


"""
`Basics.EmPropertyC(wa::Float64)`  ... constructor for an `constant` instance of EmPropertyC.
"""
function EmPropertyC(wa::Float64)
    EmPropertyC(wa, wa)
end


"""
`Basics.EmPropertyC(wa::Complex{Float64})`  ... constructor for a gauge-independent (complex) instance.

    Added 13-Aug-2026 for amplitudes: a MAGNETIC multipole has no gauge dependence, so its amplitude is the
    same number in both components. Saying that with this constructor is what removes the special case in
    which a magnetic contribution has to be added to the Coulomb AND the Babushkin sum by hand.
"""
function EmPropertyC(wa::Complex{Float64})
    EmPropertyC(wa, wa)
end


Base.:+(a::EmPropertyC, b::EmPropertyC) = EmPropertyC(a.Coulomb + b.Coulomb, a.Babushkin + b.Babushkin)
Base.:+(a::EmPropertyC, b) = EmPropertyC(a.Coulomb + b, a.Babushkin + b)
Base.:+(a, b::EmPropertyC) = b + a
Base.:-(a::EmPropertyC, b::EmPropertyC) = EmPropertyC(a.Coulomb - b.Coulomb, a.Babushkin - b.Babushkin)
Base.:*(a::EmPropertyC, b::EmPropertyC) = EmPropertyC(a.Coulomb * b.Coulomb, a.Babushkin * b.Babushkin)
Base.:*(a, b::EmPropertyC) = EmPropertyC(a * b.Coulomb, a * b.Babushkin)
Base.:*(a::EmPropertyC, b) = EmPropertyC(a.Coulomb * b, a.Babushkin * b)
Base.:/(a::EmPropertyC, b::EmPropertyC) = EmPropertyC(a.Coulomb / b.Coulomb, a.Babushkin / b.Babushkin)
Base.:/(a, b::EmPropertyC) = EmPropertyC(a / b.Coulomb, a / b.Babushkin)

## Completed 13-Aug-2026 so that EmPropertyC can carry an AMPLITUDE, of which EmProperty is the observable.
## Both act componentwise, which is the point: two gauges are two parallel calculations of the same quantity
## and must never mix.  A product amplitude * conj(amplitude') therefore keeps Coulomb with Coulomb and
## Babushkin with Babushkin BY CONSTRUCTION -- where a code that carries the gauge as a label has to
## remember to compare it (see PhotoIonization.computeAngularBeta, which does exactly that by hand).
Base.conj(a::EmPropertyC) = EmPropertyC(conj(a.Coulomb), conj(a.Babushkin))
Base.abs2(a::EmPropertyC) = EmProperty(abs2(a.Coulomb), abs2(a.Babushkin))
Base.abs(a::EmPropertyC)  = EmProperty(abs(a.Coulomb),  abs(a.Babushkin))


# `Base.show(io::IO, property::EmPropertyC)`  ... prepares a proper printout of the variable property::EmPropertyC.
function Base.show(io::IO, property::EmPropertyC) 
    sa = Base.string(property);                print(io, sa)
end


# `Base.string(property::EmPropertyC)`  ... provides a String notation for the variable property::EmPropertyC.
function Base.string(property::EmPropertyC)
    sa = "$(property.Coulomb) [Coulomb],  $(property.Babushkin) [Babushkin]"
    return( sa )
end


"""
`struct  Basics.MultipoleAmplitude`
    ... the contribution of ONE multipole of the electron-photon interaction operator to an amplitude.

    + multipole      ::EmMultipole      ... Multipole of the photon emission/absorption.
    + amplitude      ::EmPropertyC      ... Amplitude in BOTH gauges; the gauge is not a label here.

        This lives in Basics because it is the SAME type in every process that expands the electron-photon
        interaction into multipoles -- PhotoIonization, PhotoRecombination, PhotoEmission and, through the last
        of those, DielectronicRecombination. It was defined identically in three modules until 14-Aug-2026.

        What differs BETWEEN those modules is not this container but two other things: the physical amplitude
        itself, which each module computes in its own `amplitude(...)` function; and the structure AROUND this
        type -- a PartialWave with kappa, energy and phase where there is a free electron, a total symmetry
        where one kappa serves several, and nothing at all for a bound-bound transition. Those stay per-module.

        A MAGNETIC multipole simply has equal components, which is what a gauge-labelled form has to express by
        adding Basics.Magnetic to both gauge sums; and a product conj(amplitude) * amplitude' is componentwise,
        so two gauges cannot mix however the caller is written.
"""
struct  MultipoleAmplitude
    multipole        ::EmMultipole
    amplitude        ::EmPropertyC
end

export MultipoleAmplitude


# `Base.show(io::IO, ma::MultipoleAmplitude)`  ... prepares a proper printout of the variable ma.
function Base.show(io::IO, ma::MultipoleAmplitude)
    print(io, "MultipoleAmplitude($(ma.multipole), amp = $(ma.amplitude))")
end


#################################################################################################################################
#################################################################################################################################


"""
`struct  Basics.IsotopicFraction`  
    ... defines a type to deal with different -- individual or averaged -- isotopic fractions in a computation.

    + Z        ::Float64    ... nuclear charge
    + A        ::Float64    ... (mean) mass number of the isotope/isotopic mixture.
    + x        ::Float64    ... fraction 0 <= x <=1 of the given isotope/isotopic mixture.
"""
struct  IsotopicFraction 
    Z          ::Float64
    A          ::Float64
    x          ::Float64 
end

export IsotopicFraction


"""
`Basics.IsotopicFraction()`  ... constructor for an `empty` instance of IsotopicFraction.
"""
function IsotopicFraction()
    IsotopicFraction(0., 0., 0.)
end


# `Base.show(io::IO, frac::IsotopicFraction)`  ... prepares a proper printout of the variable frac::IsotopicFraction.
function Base.show(io::IO, frac::IsotopicFraction) 
    sa = Base.string(frac);                print(io, sa)
end


# `Base.string(frac::IsotopicFraction)`  ... provides a String notation for the variable frac::IsotopicFraction.
function Base.string(frac::IsotopicFraction)
    sa = "Isotopic fraction (Z=$(frac.Z), A=$(frac.A), x=$(frac.x))"
    return( sa )
end

#################################################################################################################################
#################################################################################################################################


"""
`@enum   Parity`  ... defines a enumeration for the allowed values of parity.

    + plus, minus               ... with obvious meaning
"""
@enum   Parity    plus   minus

@doc "... even parity of an atomic level or configuration."                                                                                       plus
@doc "... odd parity of an atomic level or configuration."                                                                                       minus
export  Parity,   plus,  minus


"""
`Basics.Parity(sa::String)`  ... constructor for a given String.
"""
function Parity(sa::String)
    if       sa == "+"    wa = plus
    elseif   sa == "-"    wa = minus
    else     error("stop a")
    end

    return( wa )
end


"""
`Basics.invertParity(p::Parity)`  ... inverts the given parity plus <--> minus.
"""
function invertParity(p::Parity)
    if       p == plus    return( minus )
    elseif   p == minus   return( plus )
    else     error("stop a")
    end
end


# `Base.string(p::Parity)`  ... provides a proper printout of the variable p::Parity.
function Base.string(p::Parity) 
    if      p == plus     return( "+" )
    elseif  p == minus    return( "-" )  
    else    error("stop a")
    end
end


function  Base.:*(p1::Parity, p2::Parity)
    if      p1 == plus   &&  p2 == plus    return( plus ) 
    elseif  p1 == minus  &&  p2 == minus   return( plus ) 
    else                                   return( minus ) 
    end
end

#################################################################################################################################
#################################################################################################################################


"""
`struct  Basics.Shell`  ... defines a enumeration for the allowed values of a non-relativistic shell.  

    + n     ::Int64  ... principal quantum number 
    + l     ::Int64  ... orbital angular quantum number
"""
struct  Shell
    n       ::Int64 
    l       ::Int64 
end

export  Shell


"""
`Basics.Shell(sa::Union{String,SubString{String}})`  ... constructor for a given String, such as 1s, 2s, 2p, ... .
"""
function Shell(sa::Union{String,SubString{String}}) 
    wa = strip( Base.convert(String, sa) )
    n  = parse(Int64, wa[1:end-1]);    l = shellNotation( string(wa[end]) )
    Shell(n, l)   
end


# `Base.show(io::IO, sh::Shell)`  ... prepares a proper printout of the variable sh::Shell.
function Base.show(io::IO, sh::Shell) 
    print(io, string(sh) )
end


# `Base.string(sh::Shell)`  ... provides a proper printout of the variable sh::Shell.
function Base.string(sh::Shell) 
    return( string(sh.n) * shellNotation(sh.l) )
end


"""
`Basics.shellNotation(l::Int64)`  ... returns the corresponding spectroscopy letter for a given orbital angular quantum number l.
"""
function shellNotation(l::Int64)
    !(0 <= l <= 21)      &&   return (string(l))
    wa = [ "s", "p", "d", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "q", "r", "t", "u", "v", "w", "x", "y", "z"]
    return( wa[l+1] )  
    #
end


"""
`Basics.shellNotation(sa::String)`  ... returns for a given spectroscopy letter the orbital angular quantum number l.
"""
function shellNotation(sa::String)
    wa = [ "s", "p", "d", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "q", "r", "t", "u", "v", "w", "x", "y", "z"]
    wb = findall(in([sa]), wa)
    !(length(wb)==1)      &&   error("Unrecognized orbital letter = $sa")
    return( wb[1]-1)   
end


"""
`Basics.shellSplitIntoSubshells(sh::Shell)'  ... to convert a non-relativistic shell into a list of Subshell[..]
"""
function shellSplitIntoSubshells(sh::Shell)
    subshells = Subshell[]
    if     sh.l == 0   push!(subshells, Subshell(sh.n, -1) )
    else   j2 = 2sh.l - 1;    kappa = Int64(  (j2+1)/2 );   push!(subshells, Subshell(sh.n, kappa) )
        j2 = 2sh.l + 1;    kappa = Int64( -(j2+1)/2 );   push!(subshells, Subshell(sh.n, kappa) )
    end
    return( subshells )
end


"""
`Basics.shellSplitOccupation(sh::Shell, occ::Int64)'  
    ... to split the occupation of a shell into proper subshell occupations; an  Array{Dict{Basics.Subshell,Int64},1} is returned.
"""
function shellSplitOccupation(sh::Shell, occ::Int64)

    subshells = shellSplitIntoSubshells(sh)
    wa        = Dict{Subshell,Int64}[]

    for  k1 in 0:occ
        for  k2 in 0:occ
            if  k1 + k2 == occ
                if       length( subshells ) == 1  
                    !(0 <= occ <= Basics.subshell_2j( subshells[1] ) + 1)   &&    error("Wrong pair of subshell = $(subshells[1]) and occupation = $occ") 
                    wb = Dict( subshells[1] => occ )
                    push!(wa, wb)
                    return( wa )

                    elseif  length( subshells ) == 2
                        j2_1 = Basics.subshell_2j(subshells[1])
                    j2_2 = Basics.subshell_2j(subshells[2])
                    if  k1 <= j2_1 + 1   &&    k2 <= j2_2 + 1    
                    wb = Dict( subshells[1] => k1, subshells[2] => k2 )
                    # Delete contributions with zero occupation
                    ## if  k1 == 0    delete!(wb, subshells[1])   end
                    ## if  k2 == 0    delete!(wb, subshells[2])   end
                    push!(wa, wb)
                    end
                else    error("stop a") 
                end
            end
        end
    end

    return( wa )
end

#################################################################################################################################
#################################################################################################################################


"""
`struct  Basics.Subshell`  ... defines a type for the allowed values of a relativistic subhell.  

    + n        ::Int64  ... principal quantum number 
    + kappa    ::Int64  ... relativistic angular quantum number
"""
struct  Subshell
    n          ::Int64 
    kappa      ::Int64 
end

export  Subshell, subshell_2j

"""
`Base.hash(sh::Basics.Shell, h::UInt)` and `Base.hash(sh::Basics.Subshell, h::UInt)`
    ... hash a shell or subshell from its quantum numbers alone, with a fixed seed. A value::UInt is returned.

        THESE EXIST TO MAKE Dict ITERATION ORDER REPRODUCIBLE, and that is a correctness matter rather than a
        tidiness one. Julia's default hash for a composite type folds in the TYPE's identity, which differs between
        precompiled images: measured on 28-Aug-2026, hash(Shell("1s")) was 17510954341873548245 in one build and
        15512695893170165806 in the next, for the same source. Every Dict or Set keyed by one of these types then
        iterated in a different order from build to build -- so `Basics.generateConfigurations` emitted the
        relativistic configurations of 1s^2 2s^2 2p^5 3p^1 in one order or another, the CSF basis followed, and the
        SCF total energy moved in the 9th digit for no reason a reader could see. That was priority item 57, and it
        cost a day of control runs across two migrations before it was tracked down.

        Hashing the quantum numbers with a literal seed removes the type identity from the calculation and makes the
        order a property of the physics. Consistency with equality is preserved: neither type defines its own `==`,
        so equality is field-wise and these are exactly the fields.
"""
Base.hash(sh::Shell,    h::UInt) = hash(sh.l,     hash(sh.n, h ⊻ 0x5368656c6c000001))
Base.hash(sh::Subshell, h::UInt) = hash(sh.kappa, hash(sh.n, h ⊻ 0x5375627368000002))





"""
`Basics.Subshell(sa::String)`  ... constructor for a given String, such as 1s_1/2, 2s_1/2, 2p_3/2, ... .
"""
function Subshell(sa::String) 
    wa = strip( sa );   
    wb = findnext("_", wa, 1);    if  isnothing(wb)   error("No underscore in sa = $sa")           else   wb = wb[1]   end
    wc = findnext("/", wa, wb);   if  isnothing(wc)   error("No / in wa[] = $(wa[wb[1]+1:end])")   else   wc = wc[1]   end
    wd = string( wa[wb+1:wc-1] )
    n  = parse(Int64, wa[1:wb-2]);                  l  = shellNotation( string(wa[ wb-1 ]) )
    j2 = parse(Int64, wd );                         !(string(wa[end]) == "2")    &&  error("Unrecognized subshell string sa = $sa") 
    if  l - j2/2 > 0     kappa = Int64( (j2+1)/2 )    else    kappa = Int64( -(j2+1)/2 )   end
    Subshell(n, kappa)   
end


# `Base.show(io::IO, sh::Subshell)`  ... prepares a proper printout of the variable sh::Subshell.
function Base.show(io::IO, sh::Subshell) 
    print(io, string(sh) )
end


# `Base.string(sh::Subshell)`  ... provides a proper printout of the variable sh::Subshell.
function Base.string(sh::Subshell) 
    sa = string(sh.n) * shellNotation( Basics.subshell_l(sh) ) * "_" * string( Basics.subshell_j(sh) )
    return( sa )
end


"""
`subshellGrasp(sa::AbstractString)`  ... returns a sh::Subshell if sa denotes a (relativistic) Subshell in Grasp-notation.
"""
function subshellGrasp(sa::AbstractString)
        if   sa == "1s"   return( Subshell("1s_1/2")  )   
    elseif   sa == "2s"   return( Subshell("2s_1/2")  )   elseif   sa == "2p-"   return( Subshell("2p_1/2") )
    elseif   sa == "2p"   return( Subshell("2p_3/2")  ) 
    elseif   sa == "3s"   return( Subshell("3s_1/2")  )   elseif   sa == "3p-"   return( Subshell("3p_1/2") )
    elseif   sa == "3p"   return( Subshell("3p_3/2")  )   elseif   sa == "3d-"   return( Subshell("3d_3/2") )
    elseif   sa == "3d"   return( Subshell("3d_5/2")  ) 
    elseif   sa == "4s"   return( Subshell("4s_1/2")  )   elseif   sa == "4p-"   return( Subshell("4p_1/2") )
    elseif   sa == "4p"   return( Subshell("4p_3/2")  )   elseif   sa == "4d-"   return( Subshell("4d_3/2") )
    elseif   sa == "4d"   return( Subshell("4d_5/2")  )   elseif   sa == "4f-"   return( Subshell("4f_5/2") )
    elseif   sa == "4f"   return( Subshell("4f_7/2")  ) 
    elseif   sa == "5s"   return( Subshell("5s_1/2")  )   elseif   sa == "5p-"   return( Subshell("5p_1/2") )
    elseif   sa == "5p"   return( Subshell("5p_3/2")  )   elseif   sa == "5d-"   return( Subshell("5d_3/2") )
    elseif   sa == "5d"   return( Subshell("5d_5/2")  )   elseif   sa == "5f-"   return( Subshell("5f_5/2") )
    elseif   sa == "5f"   return( Subshell("5f_7/2")  )   elseif   sa == "5g-"   return( Subshell("5g_7/2") )
    elseif   sa == "5g"   return( Subshell("5g_9/2")  ) 
    elseif   sa == "6s"   return( Subshell("6s_1/2")  )   elseif   sa == "6p-"   return( Subshell("6p_1/2") )
    elseif   sa == "6p"   return( Subshell("6p_3/2")  )   elseif   sa == "6d-"   return( Subshell("6d_3/2") )
    elseif   sa == "6d"   return( Subshell("6d_5/2")  )   elseif   sa == "6f-"   return( Subshell("6f_5/2") )
    elseif   sa == "6f"   return( Subshell("6f_7/2")  )   elseif   sa == "6g-"   return( Subshell("6g_7/2") )
    elseif   sa == "6g"   return( Subshell("6g_9/2")  ) 
    elseif   sa == "7s"   return( Subshell("7s_1/2")  )   elseif   sa == "7p-"   return( Subshell("7p_1/2") )
    elseif   sa == "7p"   return( Subshell("7p_3/2")  )   elseif   sa == "7d-"   return( Subshell("7d_3/2") )
    else     error("Grasp notation ... not yet implemented; sa = $sa")
    end
end


"""
`Basics.subshell_l(sh::Subshell)`  ... returns the orbital angular quantum number l for a given sh::Subshell; an Int64 is returned
"""
function subshell_l(sh::Subshell)
    if   sh.kappa < 0   l = abs(sh.kappa) -1   else   l = sh.kappa   end
    return( l )   
end


"""
`Basics.subshell_j(sh::Subshell)`  ... returns the total angular quantum number j for a given sh::Subshell; an AngularJ64 is returned.
"""
function subshell_j(sh::Subshell)
    return( AngularJ64(abs(sh.kappa) - 1//2) )   
end


"""
`Basics.subshell_2j(sh::Subshell)`  ... returns the total angular quantum number j for a given sh::Subshell; an AngularJ64 is returned.
"""
function subshell_2j(sh::Subshell)
    return( 2* abs(sh.kappa) - 1 )   
end

#################################################################################################################################
#################################################################################################################################


"""
`struct   Basics.SubshellStateR`  
    ... defines a struct for the relativistic antisymmetric subshell states in the seniority scheme; this struct can be accessed 
        only internally and, therefore, the only the standard constructor is supported.

    + subshell     ::Subshell        ... to refer to a particular subshell. 
    + occ          ::Int64           ... occupation of the subshell
    + nu           ::Int64           ... seniority of the subshell state, often called nu
    + Jsub2        ::Int64           ... 2*Jsub
"""
struct  SubshellStateR
    subshell       ::Subshell
    occ            ::Int64
    nu             ::Int64
    Jsub2          ::Int64
end     

export  SubshellStateR


# `Base.show(io::IO, state::SubshellStateR)`  ... prepares a proper printout of the variable state::SubshellStateR.
function Base.show(io::IO, state::SubshellStateR) 
    print(io, string(state) )
end


# `Base.string(state::SubshellStateR)`  ... provides a proper printout of the variable state::SubshellStateR.
function Base.string(state::SubshellStateR)
    sa = "SubshellState: [" * string(state.subshell) * "^$(state.occ), nu=$(state.nu), 2*J_sub=$(state.Jsub2)]"
end

#################################################################################################################################
#################################################################################################################################


"""
`@enum   Basics.UseGauge`  ... defines a enumeration for the (allowed) gauges to be selected in explicit computations
"""
@enum   UseGauge    UseCoulomb   UseBabushkin

@doc "... to select the Coulomb gauge in an explicit computation."                                                                          UseCoulomb
@doc "... to select the Babushkin gauge in an explicit computation."                                                                      UseBabushkin
export  UseGauge,   UseCoulomb,  UseBabushkin


"""
`Basics.UseGauge(sa::String)`  ... constructor for a given String.
"""
function UseGauge(sa::String)
    if       sa == "Coulomb"                   wa = UseCoulomb
    elseif   sa == "Babushkin"                 wa = UseBabushkin
    else     error("stop a")
    end

    UseGauge(wa)
end


# `Base.string(gauge::UseGauge)`  ... provides a proper printout of the variable gauge::EmGauge.
function Base.string(gauge::UseGauge) 
    if      gauge == UseCoulomb     return( "Coulomb" )   
    elseif  gauge == UseBabushkin   return( "Babushkin" )   
    else    error("stop a")
    end
end


#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################


"""
`struct  Basics.CorePolarization`  
    ... defines a type to describe the influence of core-polarization upon the potential or interaction operator.

    + doApply          ::Bool                ... True if core-polarization is to be applied, and false otherwise.
    + coreAlpha        ::Float64             ... Polarizibility of the ionic core.
    + coreRadius       ::Float64             ... Radius of the ionic core.
    + coreShells       ::Array{Shell,1}      ... Shells of the ionic core
"""
struct  CorePolarization
    doApply            ::Bool     
    coreAlpha          ::Float64 
    coreRadius         ::Float64
    coreShells         ::Array{Shell,1} 
end 

export  CorePolarization


"""
`Basics.CorePolarization()`  ... constructor for the default values of core-polarization contributions.
"""
function CorePolarization()
    CorePolarization(false, 0., 0., Shell[])
end


# `Base.show(io::IO, cp::CorePolarization)`  ... prepares a proper printout of the variable cp::CorePolarization.
function Base.show(io::IO, cp::CorePolarization) 
    sa = Base.string(cp);                print(io, sa)
end


# `Base.string(cp::CorePolarization)`  ... provides a String notation for the variable cp::CorePolarization.
function Base.string(cp::CorePolarization)
    if  cp.doApply  sa = "Apply "   else    sa = "Do not apply "     end
    sa = sa * "core-polarization with alpha_c = $(cp.coreAlpha) a.u., r_c = $(cp.coreRadius) a.u. " *
                "core shells = $(cp.coreShells)."
    return( sa )
end


#################################################################################################################################
#################################################################################################################################

"""
`struct  Basics.Configuration`  
    ... defines a struct for a non-relativistic electron configuration that is fully speficied by its shell notations 
        (such as "1s", "2s", "2p", ....) and the corresponding occupation numbers (>= 0). An electron configuration 
        is independent of any order of the shells, and a zero occupation is assumed for all shells that do not appear as a 
        key in the shell (dictionary). Therefore, the number of keys do not allow any conclusion about the underlying orbital 
        space of any considered computation that include more than just a single configuration.

    + shells         ::Dict{Shell,Int64}   ... Dictionary that maps shells to their occupation.
    + NoElectrons    ::Int64               ... No. of electrons.
"""
struct  Configuration
    shells           ::Dict{Shell,Int64}  
    NoElectrons      ::Int64               
end

export  Configuration 


"""
`Basics.Configuration(NoElectrons::Int64)`
    ... constructor for a given number of electrons which are filled into the shells in the MADELUNG order, i.e. by
        increasing n+l and, at equal n+l, by increasing n. A conf::Configuration is returned.

        This is the (n+l, n) rule and not the strict (n, l) order that was used until 18-Aug-2026, which placed the
        electrons of calcium into 3d instead of 4s and gave `1s^2 2s^2 2p^6 3s^2 3p^6 3d^2`; strontium came out as
        4d^2 and xenon as 4f^8. Madelung reproduces the ground configuration of most elements, but NOT of the roughly
        twenty exceptions -- Cr, Cu, Nb, Mo, Ru, Rh, Pd, Ag, La, Ce, Gd, Pt, Au, Ac, Th, Pa, U, Np, Cm -- where a
        d or f shell is filled at the expense of the outer s. A caller that needs the true ground configuration of a
        specific element must therefore not rely on this constructor alone.
"""
function Configuration(NoElectrons::Int64)
    shellDict = Dict{Shell,Int64}();  N = NoElectrons
    shellList = Basics.generateShellList(1, 8, 7)
    sort!(shellList, by = sh -> (sh.n + sh.l, sh.n))
    for  sh in shellList
        maxocc = 2 * (2*sh.l + 1);    occ    = min(maxocc, N);    shellDict[sh] = occ
        N      = N - occ;             if  N == 0   break   end
    end

    return ( Configuration(shellDict, NoElectrons) )
end


"""
`Basics.Configuration(sa::String)`  
    ... constructor for a given configuration string, such as "[He]", "[Ne]", "[Ne] 3s 3p^6"  or "1s 2p^6 3s^2 3p".
        One can also use "1s^0" to represent a (bare-ion) empty configuration.
"""
function Configuration(sa::String)
    sax = strip(sa)
    shellOccList = split(sax, " ")
    NoElectrons = 0;    shellDict    = Dict{Shell,Int64}()
    for shocc in shellOccList
        if       shocc == " "    continue
        elseif   shocc  in [ "[He]", "[Ne]", "[Ar]", "[Kr]", "[Pd]", "[Xe]", "[Er]", "[Pt]", "[Hg]", "[Rn]", "[Fm]" ]
            shellDict = Base.merge( shellDict, Dict( Shell("1s") => 2));                                       NoElectrons = NoElectrons + 2
            shocc == "[He]"  &&  continue
            shellDict = Base.merge( shellDict, Dict(Shell("2s") => 2,  Shell("2p") => 6));                     NoElectrons = NoElectrons + 8  
            shocc == "[Ne]"  &&  continue
            shellDict = Base.merge( shellDict, Dict(Shell("3s") => 2,  Shell("3p") => 6));                     NoElectrons = NoElectrons + 8  
            shocc == "[Ar]"  &&  continue 
            shellDict = Base.merge( shellDict, Dict(Shell("3d") => 10, Shell("4s") => 2, Shell("4p") => 6));   NoElectrons = NoElectrons + 18   
            shocc == "[Kr]"  &&  continue   
            shellDict = Base.merge( shellDict, Dict(Shell("4d") => 10));                                       NoElectrons = NoElectrons + 10  
            shocc == "[Pd]"  &&  continue 
            shellDict = Base.merge( shellDict, Dict(Shell("5s") => 2, Shell("5p") => 6));                      NoElectrons = NoElectrons + 8  
            shocc == "[Xe]"  &&  continue 
            shellDict = Base.merge( shellDict, Dict(Shell("4f") => 14));                                       NoElectrons = NoElectrons + 14  
            shocc == "[Er]"  &&  continue 
            shellDict = Base.merge( shellDict, Dict(Shell("5d") => 10));                                       NoElectrons = NoElectrons + 10  
            shocc == "[Pt]"  &&  continue 
            shellDict = Base.merge( shellDict, Dict(Shell("6s") => 2));                                        NoElectrons = NoElectrons + 2  
            shocc == "[Hg]"  &&  continue 
            shellDict = Base.merge( shellDict, Dict(Shell("6p") => 6));                                        NoElectrons = NoElectrons + 6  
            shocc == "[Rn]"  &&  continue 
            shellDict = Base.merge( shellDict, Dict(Shell("5f") => 14));                                       NoElectrons = NoElectrons + 14  
            shocc == "[Fm]"  &&  continue 
        else
            ia = findall(in("^"), shocc);    
            if     length(ia) == 0     shellDict = Base.merge( shellDict, Dict( Shell(shocc) => 1));           NoElectrons = NoElectrons + 1
            else   occ = parse(Int64, shocc[ia[1]+1:end] );    sh = Shell( shocc[1:ia[1]-1] )
                    shellDict = Base.merge( shellDict, Dict( sh => occ));                                       NoElectrons = NoElectrons + occ
            end
        end
    end
    return( Configuration(shellDict, NoElectrons) )
end


# `Base.show(io::IO, conf::Configuration)`  ... prepares a proper printout of the non-relativistic configuration conf::Configuration.
function Base.show(io::IO, conf::Configuration)
    sa = "Configuration: " * string(conf)
    print(io, sa)
end


# `Base.unique(confs::Array{Configuration,1})`  ... return a unique list of configurations.
function Base.unique(confs::Array{Configuration,1})
    confList       = Configuration[]  
    configurations = deepcopy(confs)
    for  confa in configurations
        addTo = true
        for confb in confList    if   confa == confb    addTo = false;    break     end     end
        if  addTo    push!(confList, confa)     end
    end

    return( confList )
end


# `Base.string(conf::Configuration)`  ... provides a String notation for the variable conf::Configuration.
function Base.string(conf::Configuration)
    wa = keys(conf.shells);   va = values(conf.shells)
    shells = Shell[];   for  k in wa    push!(shells, k)        end
    sortedShells = Base.sort( shells , lt=Base.isless)

    sa = ""
    for  k  in sortedShells 
        occ = conf.shells[k]
        sa  = sa * string(k) * "^$occ "
    end

    return( sa )
end


"""
`Base.string(conf::Configuration, open::Bool)`  
    ... provides, if open=trueM a String notation for the open-shell part of variable conf::Configuration,
        and the standard printout for open=false.
"""
function Base.string(conf::Configuration, open::Bool)
    if    false   sa = Base.string(conf)
    else
        # Define a new configurations with shells "beyond" the closed core
        openShell = Shell(100,99)
        for  (shell,occ)  in conf.shells
            if  occ <  2*(2*shell.l + 1)  &&  shell < openShell    openShell = shell   end 
        end
        #
        NoElectrons = 0;       newShells  = Dict{Shell,Int64}()
        CoreElectrons = 0;     coreShells = Dict{Shell,Int64}()
        for  (shell,occ)  in conf.shells
            if    shell >= openShell    newShells[shell]  = occ;    NoElectrons = NoElectrons + occ   
            else                        coreShells[shell] = occ;    CoreElectrons = CoreElectrons + occ 
            end 
        end
        newConf  = Configuration(newShells, NoElectrons)
        coreConf = Configuration(coreShells, CoreElectrons)
        if      coreConf == Configuration("[He]")                   sb = "[He] "
        elseif  coreConf == Configuration("[He] 2s^2")              sb = "[Be] "
        elseif  coreConf == Configuration("[Ne]")                   sb = "[Ne] "
        elseif  coreConf == Configuration("[Ne] 3s^2")              sb = "[Mg] "
        elseif  coreConf == Configuration("[Ar]")                   sb = "[Ar] "
        elseif  coreConf == Configuration("[Ar] 4s^2")              sb = "[Ca] "
        elseif  coreConf == Configuration("[Ar] 3d^10")             sb = "[Ni-like] "
        elseif  coreConf == Configuration("[Ar] 3d^10 4s^2")        sb = "[Zn-like] "
        elseif  coreConf == Configuration("[Kr]")                   sb = "[Kr] "
        elseif  coreConf == Configuration("[Kr] 5s^2")              sb = "[Sr] "
        elseif  coreConf == Configuration("[Kr] 4d^10")             sb = "[Pd-like] "
        elseif  coreConf == Configuration("[Kr] 4d^10 5s^2")        sb = "[Cd-like] "
        elseif  coreConf == Configuration("[Kr] 4d^10 4f^14")       sb = "[Nd-like] "
        elseif  coreConf == Configuration("[Kr] 4d^10 4f^14 5s^2")  sb = "[Sm-like] "
        elseif  coreConf == Configuration("[Xe]")                   sb = "[Xe] "
        elseif  coreConf == Configuration("[Xe] 6s^2")              sb = "[Ba] "
        elseif  coreConf == Configuration("[Xe] 5d^10")             sb = "[Gd-like] "
        elseif  coreConf == Configuration("[Xe] 4f^14")             sb = "[Er-like] "
        elseif  coreConf == Configuration("[Xe] 5d^10 6s^2")        sb = "[Dy-like] "
        elseif  coreConf == Configuration("[Xe] 4f^14 5d^10")       sb = "[Pt-like] "
        elseif  coreConf == Configuration("[Hg]")                   sb = "[Hg] "
        else                                                        sb = "[Core] "
        end 
        #
        sa = sb * Base.string(newConf)
    end

    return( sa )
end


"""
`Base.:(==)(confa::Configuration, confb::Configuration)`  
    ... compares (recursively) two non-relativistic configurations and return true if all subfields are equal, 
        and false otherwise
"""
function  Base.:(==)(confa::Configuration, confb::Configuration)
    if   confa.NoElectrons  != confb.NoElectrons    return( false )    end
    allShells = Base.merge(confa.shells, confb.shells)
    wk        = keys(allShells)
    for  k in wk
        if  haskey(confa.shells, k)   &&   haskey(confb.shells, k)  &&   
                                            confa.shells[k] != confb.shells[k]     return( false )    end
        if  !haskey(confa.shells, k)  &&   confb.shells[k] != 0                   return( false )    end
        if  !haskey(confb.shells, k)  &&   confa.shells[k] != 0                   return( false )    end
    end
    return( true )
end


#################################################################################################################################
#################################################################################################################################


"""
`struct  Basics.ConfigurationR`  
    ... defines a type for a relativistic electron configuration that is fully speficied by its subshell notations 
        (such as "1s_1/2", "2s_1/2", "2p_3/2", ....) and the corresponding occupation numbers (>= 0). An electron 
        configuration is independent of any order of the subshells, and a zero occupation is assumed for all shells that 
        do not appears as a key. Therefore, the number of keys do not allow any conclusion about the underlying orbital 
        space of any considered computation that include more than a single configuration.

    + subshells      ::Dict{Subshell,Int64}  ... Dictionary that maps subshells to their occupation.
    + NoElectrons    ::Int64                 ... No. of electrons.
"""
struct  ConfigurationR
    subshells        ::Dict{Subshell,Int64}  
    NoElectrons      ::Int64               
end

export  ConfigurationR 

# `Base.show(io::IO, conf::ConfigurationR)`  ... prepares a proper printout of the relativistic configuration conf::ConfigurationR.
function Base.show(io::IO, conf::ConfigurationR)
    sa = "Configuration: " * string(conf)
    print(io, sa)
end


# `Base.string(conf::ConfigurationR)`  ... provides a String notation for the variable conf::ConfigurationR.
function Base.string(conf::ConfigurationR)
    wa = keys(conf.subshells);   va = values(conf.subshells)
    subshells = Subshell[];   for  k in wa    push!(subshells, k)    end
    sortedSubshells = Base.sort( subshells , lt=Base.isless)

    sa = ""
    for  k  in sortedSubshells 
        occ = conf.subshells[k]
        sa  = sa * string(k) * "^$occ "
    end

    return( sa )
end


"""
`Base.:(==)(confa::ConfigurationR, confb::ConfigurationR)`  
    ... compares (recursively) two relativistic configurations and return true if all subfields are equal, 
        and false otherwise
"""
function  Base.:(==)(confa::ConfigurationR, confb::ConfigurationR)
    if   confa.NoElectrons  != confb.NoElectrons    return( false )    end
    wk = keys(confa.subshells)
    for  k in wk
        if  !haskey(confb.subshells, k)  ||   confa.subshells[k] != confb.subshells[k]    return( false )    end
    end
    return( true )
end

