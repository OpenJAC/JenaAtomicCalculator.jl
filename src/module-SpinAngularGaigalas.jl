#
# This module has mainly been contributed by Gediminas Gaigalas (gediminas.gaigalas@tfai.vu.lt)
"""
`module  JAC.SpinAngularGaigalas`  
... a submodel of JAC that contains all methods for computing the spin-angular coefficients for one- and two-particle operators 
    of given rank and for symmetry-adapted CSF. At the present, these spin-angular can be obtained for all scalar (and symmetric) 
    one and two-particle interactions and they are sometimes refered to as 'pure angular coefficients' in order to distinguish 
    them from those used in GRASP (and its derivatives) where part of the `physical interaction' were included orginal 
    into the angular integrals. The spin-angular coefficients usally appear in the form:
    
        <CSF_l | operator | CSF_r> = sum_{t} coeff_t (L,abcd) * X^L (abcd)
        
    where CSF is a standard (jj-coupled) configuration state function, t a summation index, and where X^L(abcd) denotes the 
    effective interaction strength. 
        STATUS, 28-Aug-2026: SUPERSEDED BY `SpinAngular`, AND KEPT DELIBERATELY. Every caller in `src/` was moved
        to the newer module between 27 and 28-Aug; nothing here is called by JAC any more. It is retained for two
        reasons, and NOT to be picked up again without asking:

        1. IT IS THE CROSS-CHECK ORACLE. `tools/diag-spinangular-compare.jl` runs the two implementations against
           each other coefficient by coefficient, and that comparison found three real defects in two days -- a
           rank-0 substitution silently dropped, two-particle coefficients lost above j = 9/2, and a parity
           argument that eight call sites had wrong. Deleting this module would leave GRASP2018, which needs a
           Fortran build, as the only independent check.
        2. It is twenty years of results, and the newer module has no direct test of its own.

        WHAT IS **NOT** TRUE OF IT, so that a reader does not mistake age for correctness: it IGNORES the parity
        argument of its own operators and returns both parity classes; it emits rank-0 coefficients between CSFs
        of OPPOSITE parity, which a scalar operator cannot connect; and its quasispin tables -- now shared, in
        `SpinAngularTables` -- stop at j = 9/2. The newer module is right on all three.
"""
module SpinAngularGaigalas

using  Printf, ..AngularMomentum, ..Basics,  ..Defaults, ..ManyElectron, ..Radial, ..SpinAngularTables

export  Xronecker

"""
`abstract type SpinAngularGaigalas.AbstractAngularType` 
    ... defines an abstract type and a number of data types to work with one- and two-particle operators of given rank, see also:
    
    + struct OneParticleOperator    ... to represent a single-particle operator with well-defined spherical tensor rank.
    + struct TwoParticleOperator    ... to represent a two-particle operator with well-defined spherical tensor rank.
"""
abstract type  AbstractAngularType                           end

"""
`struct  SpinAngularGaigalas.OneParticleOperator  <:  AbstractAngularType`  
    ... a struct for defining the spherial tensor structure of a single-particle operator.

    + rank            ::Int64             ... Rank of the operator.
    + parity          ::Basics.Parity     ... Parity of the operator (if needed ??)
    + sameOrbitalSet  ::Bool              ... True if orbitals for both CSF are taken from the same orbital set (if needed ??) 
"""
struct  OneParticleOperator  <:  AbstractAngularType 
    rank              ::Int64         
    parity            ::Basics.Parity 
    sameOrbitalSet    ::Bool
end

"""
`SpinAngularGaigalas.OneParticleOperator()`  ... constructor for setting the default values.
"""
function OneParticleOperator()
    OneParticleOperator( 0, Basics.plus, false)
end

# `Base.show(io::IO, op::OneParticleOperator)`  ... prepares a proper printout of the op::OneParticleOperator.
function Base.show(io::IO, op::OneParticleOperator)
        sa = "One-particle operator O^($(op.rank) $(op.parity))"
        if op.sameOrbitalSet  sa = sa * ", defined for CSF from the same orbital set."
        else                  sa = sa * ", defined for CSF from different orbital sets."
        end
        println(io, sa)
end

"""
`struct  SpinAngularGaigalas.TwoParticleOperator`  
    ... a struct for defining the spherial tensor structure of a two-particle operator.

    + rank            ::Int64             ... Rank of the operator (if needed ??).
    + parity          ::Basics.Parity     ... Parity of the operator (if needed ??)
    + sameOrbitalSet  ::Bool              ... True if orbitals for both CSF are taken from the same orbital set (if needed ??) 
"""
struct  TwoParticleOperator
    rank              ::Int64         
    parity            ::Basics.Parity 
    sameOrbitalSet    ::Bool
end

"""
`SpinAngularGaigalas.TwoParticleOperator()`  ... constructor for setting the default values.
"""
function TwoParticleOperator()
    TwoParticleOperator( 0, Basics.plus, false)
end

# `Base.show(io::IO, op::TwoParticleOperator)`  ... prepares a proper printout of the op::TwoParticleOperator.
function Base.show(io::IO, op::TwoParticleOperator)
        sa = "Two-particle operator O^($(op.rank) $(op.parity))"
        if op.sameOrbitalSet  sa = sa * ", defined for CSF from the same orbital set."
        else                  sa = sa * ", defined for CSF from different orbital sets."
        end
        println(io, sa)
end

"""
`struct  SpinAngularGaigalas.Coefficient1p`  
    ... a struct for defining a single spin-angular coefficient for a reduced one-particle matrix element <a || o^(L) || b>.

    + nu       ::Int64             ... Rank (L or nu) of the single-particle interaction strength.
    + a        ::Subshell          ... Left-hand subshell (orbital).
    + b        ::Subshell          ... Right-hand subshell (orbital).
    + T        ::Float64           ... (Value of) spin-angular coefficient.
"""
struct  Coefficient1p
    nu         ::Int64    
    a          ::Subshell 
    b          ::Subshell
    T          ::Float64
end

# `Base.show(io::IO, coeff::Coefficient1p)`  ... prepares a proper printout of the coeff::Coefficient1p.
function Base.show(io::IO, coeff::Coefficient1p)
    # sa = "\n V^($(coeff.L)) [$(coeff.a), $(coeff.b)] = $(coeff.v)"
    sa = "   T^$(coeff.nu) [$(coeff.a), $(coeff.b)] = $(coeff.T)"
    print(io, sa)
end

"""
`struct  SpinAngularGaigalas.Coefficient2p`  
    ... a struct for defining a single spin-angular coefficient for a reduced two-particle matrix element <ab || o^(L) || cd>,
        such as the Coulomb interaction strength.

    + nu       ::Int64             ... Rank (L or nu) of the single-particle interaction strength.
    + a        ::Subshell          ... Left-hand subshell (orbital).
    + b        ::Subshell          ... Left-hand subshell (orbital).
    + c        ::Subshell          ... Right-hand subshell (orbital).
    + d        ::Subshell          ... Right-hand subshell (orbital).
    + V        ::Float64           ... (Value of) spin-angular coefficient.
"""
struct  Coefficient2p
    nu         ::Int64    
    a          ::Subshell 
    b          ::Subshell
    c          ::Subshell 
    d          ::Subshell  
    V          ::Float64
end


# `Base.show(io::IO, term::QspaceTerm)`  ... prepares a proper printout of term::QspaceTerm.
function Base.show(io::IO, term::QspaceTerm)
    if  term.Nr == 0   
        sa = "|$(term.j) $(term.Q) $(term.J) min_odd=$(term.min_odd) max_odd=$(term.max_odd) " * 
                "min_even=$(term.min_even) max_even=$(term.max_even)>"
    else               
        sa = "|$(term.j) $(term.Q) $(term.J); Nr=$(term.Nr) min_odd=$(term.min_odd) " *
                "max_odd=$(term.max_odd) min_even=$(term.min_even) max_even=$(term.max_even) >"
    end
    println(io, sa)
end

# `Base.show(io::IO, coeff::Coefficient2p)`  ... prepares a proper printout of the coeff::Coefficient2p.
function Base.show(io::IO, coeff::Coefficient2p)
    sa = "\n V^($(coeff.nu)) [$(coeff.a), $(coeff.b)| $(coeff.c), $(coeff.d)] = $(coeff.V)"
    print(io, sa)
end

"""
`struct  SpinAngularGaigalas.SchemeEta`  
    ... to define various singleton structs in order to distinguish between different irreducible tensors 
        and matrix elements.

    + Eta_a         ... a^(qj)_mq
    + Eta_W         ... W^(k_12) = [a x a]
    + Eta_aW        ... W^(k_12, k_2) = [a1 x [a2 x a3]]
    + Eta_Wa        ... W^(k_12, k_2) = [[a1 x a2] x a3]
    + Eta_WW        ... W^(kk0) = [[a1 x a2]^k x [a3 x a4]^k]^0
"""
struct  SchemeEta_a    end
struct  SchemeEta_W    end
struct  SchemeEta_aW   end
struct  SchemeEta_Wa   end
struct  SchemeEta_WW   end

"""
`struct  SpinAngularGaigalas.Diagram`  
    ... to defines various singleton() structs in order to distinguish between different coupling schemes of the 
        matrix elements.
"""
struct  DiagramC01   end
struct  DiagramC02   end
struct  DiagramC03   end
struct  DiagramC04   end
struct  DiagramC05   end

## Sentinel type for the ordered-shells overload of recoupling2p; distinguishes it from the unordered overload.
struct  Ordered2p    end


#######################################################################################################################
#######################################################################################################################

"""
`SpinAngularGaigalas.computeCoefficients(op::SpinAngularGaigalas.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1})`  
    ... computes the spin-angular coefficients for the reduced one-particle matrix element

        <leftCsf || op^(k) || rightCsf >  
        if both CSF refer to the same list of subshells; a list of one-particle coefficients coeffs::Array{Coefficient1p,1} is returned.

        leftCsf   - is the bra confuguration state function in relativistic notation;
        rightCsf  - is the ket confuguration state function in relativistic notation;
        subshells - is explicitly given relativistic subshell list if useStandardSubshells == false.
"""
function  computeCoefficients(op::SpinAngularGaigalas.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1})
    if      op.rank == 0   coeffs = computeCoefficientsScalar(op,    leftCsf, rightCsf, subshells)
    elseif  op.rank >  0   coeffs = computeCoefficientsNonScalar(op, leftCsf, rightCsf, subshells)
    else    error("stop a")    
    end
    
    return (coeffs)
end

"""
`SpinAngularGaigalas.computeCoefficients(op::SpinAngularGaigalas.TwoParticleOperator, leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1})`  
    ... computes the spin-angular coefficients for the reduced two-particle matrix element

        <leftCsf || op^(k) || rightCsf >
        if both CSF refer to the same list of subshells; a Tuple of two lists with one- and two-particle coefficients  tpl::Tuple{coeffs1p::Array{Coefficient1p,1}, coeffs2p::Array{Coefficient2p,1}}  is returned. 

        leftCsf   - is the bra confuguration state function in relativistic notation;
        rightCsf  - is the ket confuguration state function in relativistic notation;
        subshells - is explicitly given relativistic subshell list if useStandardSubshells == false.
"""
function  computeCoefficients(op::SpinAngularGaigalas.TwoParticleOperator, leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1})
    if      op.rank == 0   coeffs = computeCoefficientsScalar(op, leftCsf, rightCsf, subshells)
    else    error("stop a")    
    end
    
    return (coeffs)
end

"""
`SpinAngularGaigalas.computeCoefficientsNonScalar(op::SpinAngularGaigalas.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1})`  
    ... computes the spin-angular coefficients for the reduced (scalar) one-particle matrix element

        <leftCsf || op^(0) || rightCsf >
        if both CSF refer to the same list of subshells; a list of one-particle coefficients coeffs::Array{Coefficient1p,1} is returned.

        leftCsf   - is the bra confuguration state function in relativistic notation;
        rightCsf  - is the ket confuguration state function in relativistic notation;
        subshells - is explicitly given relativistic subshell list if useStandardSubshells == false.
"""
function  computeCoefficientsNonScalar(op::SpinAngularGaigalas.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1})
    coeffs1p = Coefficient1p[];   kj = 2 * op.rank;    kj = AngularJ64(kj//2)
    if AngularMomentum.triangularDelta(leftCsf.J, kj, rightCsf.J) == 0 return( coeffs1p ) end
    wa = 0.0;   diff_occ = 0;   count = 0;   creation = 0;   annihilation = 0
    # Cycle through all subshells of the bra- and ket-CSF
    for  (i,shi)  in  enumerate(subshells)
        diff = leftCsf.occupation[i]-rightCsf.occupation[i]
        if abs(diff) > 1   return( coeffs1p )
        elseif diff != 0
            diff_occ = diff_occ + abs(diff);   count = count + 1
        end
        if diff == 1
            if creation == 0;        creation = i;        else return( coeffs1p ) end
        elseif diff == -1
            if annihilation == 0;    annihilation = i;    else return( coeffs1p ) end
        end 
    end
    if count == 0 
        for  (i,shi)  in  enumerate(subshells)
            Ni  = leftCsf.occupation[i];               Nj  = rightCsf.occupation[i]
            if Nj != 0
                ji  = Basics.subshell_j(shi)
                Ji  = leftCsf.subshellJ[i];            Jj  = rightCsf.subshellJ[i]
                si  = leftCsf.seniorityNr[i];          sj  = rightCsf.seniorityNr[i]
                Qi  = qshellTermQ(ji,si);              Qj  = qshellTermQ(ji,sj)
                aT  = SpinAngularGaigalas.getTermNumber(ji, Ni, Qi, Ji)
                aTp = SpinAngularGaigalas.getTermNumber(ji, Nj, Qj, Jj)
                if recouplingCheck(leftCsf, rightCsf, kj, i, i, length(subshells)) != 0.0 
                    wa = recoupling1p(leftCsf, rightCsf, kj, i, length(subshells))
                    if abs(wa) >= 2.0e-10
                        wa = wa * irreducibleTensor(SchemeEta_W(),aT, Ni, AngularM64(1//2), AngularM64(-1//2), op.rank, aTp, Nj)
                        if abs(wa) >= 2.0e-10
                            #  Pure one-particle spin-angular coefficient
                            wa = -wa/sqrt(2. * op.rank + 1.)
                            #  GRASP like spin-angular coefficient
                            wa = wa * sqrt(Basics.twice(ji) + 1.)
                            push!(coeffs1p, Coefficient1p(op.rank, shi, shi, wa))
                        end
                    end
                end
            end
        end
    elseif count == 2  &&  diff_occ == 2
        if creation != annihilation
            shii = subshells[creation];               ja  = Basics.subshell_j(shii)
            shjj = subshells[annihilation];           jb  = Basics.subshell_j(shjj)
            no_one = min(creation, annihilation);     no_two = max(creation, annihilation)
            if recouplingCheck(leftCsf, rightCsf, kj, no_one, no_two, length(subshells)) != 0.0 
                if creation == no_one && annihilation == no_two
                    j1 = ja;      j2 = jb;   wa = 1.0
                elseif creation == no_two && annihilation == no_one
                    j1 = jb;      j2 = ja
                    wa = (-1)^Int64((Basics.twice(ja)+Basics.twice(jb)-2 * op.rank + 2)/2)
                else    error("SpinAngularGaigalas.computeCoefficientsNonScalar: stop a")    
                end
                #
                wa = wa * recoupling1p(leftCsf, rightCsf, j1, j2, kj, no_one, no_two, length(subshells))
                if abs(wa) >= 2.0e-10
                    aN  = leftCsf.occupation[creation];         aNp = rightCsf.occupation[creation]
                    Ji  = leftCsf.subshellJ[creation];          Jj  = rightCsf.subshellJ[creation]
                    si  = leftCsf.seniorityNr[creation];        sj  = rightCsf.seniorityNr[creation]
                    Qi  = qshellTermQ(ja,si) ;                  Qj  = qshellTermQ(ja,sj)
                    aT  = SpinAngularGaigalas.getTermNumber(ja, aN, Qi, Ji)
                    aTp = SpinAngularGaigalas.getTermNumber(ja, aNp, Qj, Jj)
        
                    bN  = leftCsf.occupation[annihilation];     bNp = rightCsf.occupation[annihilation]
                    Ji  = leftCsf.subshellJ[annihilation];      Jj  = rightCsf.subshellJ[annihilation]
                    si  = leftCsf.seniorityNr[annihilation];    sj  = rightCsf.seniorityNr[annihilation]
                    Qi  = qshellTermQ(jb,si);                   Qj  = qshellTermQ(jb,sj)
                    bT  = SpinAngularGaigalas.getTermNumber(jb, bN, Qi, Ji)
                    bTp = SpinAngularGaigalas.getTermNumber(jb, bNp, Qj, Jj)
                    wa  = wa * SpinAngularGaigalas.irreducibleTensor(SchemeEta_a(),aT, aN, AngularM64(1//2), aTp, aNp) *
                                SpinAngularGaigalas.irreducibleTensor(SchemeEta_a(),bT, bN,AngularM64(-1//2), bTp, bNp)
                    if abs(wa) >= 2.0e-10
                        occup = 0
                        for i = no_one:no_two-1    occup = occup + leftCsf.occupation[i]    end 
                        wa = (-1)^(occup+1) * wa 
                        # Purer_one-particle spin-angular coefficient
                        wa = -wa/sqrt(2. * op.rank + 1.)
                        # GRASP like spin-angular coefficient
                        wa = wa * sqrt(Basics.twice(ja) + 1.)
                        push!(coeffs1p, Coefficient1p(op.rank, shii, shjj, wa)) 
                    end
                end
            end
        end
    end
    return( coeffs1p )
end

"""
`SpinAngularGaigalas.computeCoefficientsScalar(op::SpinAngularGaigalas.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1})`  
    ... computes the spin-angular coefficients for the reduced (scalar) one-particle matrix element

        <leftCsf || op^(0) || rightCsf >
        if both CSF refer to the same list of subshells; a list of one-particle coefficients coeffs::Array{Coefficient1p,1} is returned.

        leftCsf   - is the bra confuguration state function in relativistic notation;
        rightCsf  - is the ket confuguration state function in relativistic notation;
        subshells - is explicitly given relativistic subshell list if useStandardSubshells == false.
"""
function  computeCoefficientsScalar(op::SpinAngularGaigalas.OneParticleOperator, leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1})
    coeffs1p = Coefficient1p[]
    if leftCsf.J != rightCsf.J return( coeffs1p ) end
    wa = 0.0; diff_occ = 0; count = 0; creation = 0; annihilation = 0
    # Cycle through all subshells of the bra- and ket-CSF
    for  (i,shi)  in  enumerate(subshells)
        diff = leftCsf.occupation[i]-rightCsf.occupation[i]
        if abs(diff) > 1   return( coeffs1p )
        elseif diff != 0
            diff_occ = diff_occ + abs(diff);   count = count + 1
        end
        if diff == 1
            if creation == 0;      creation = i;        else return( coeffs1p ) end
        elseif diff == -1
            if annihilation == 0;  annihilation = i;    else return( coeffs1p ) end
        end 
    end
    #
    if count == 0 
        for  (i,shi)  in  enumerate(subshells)
            if recouplingCheck(leftCsf, rightCsf, i, i, i, i, length(subshells), 0) != 0.0
                Ni  = leftCsf.occupation[i];               Nj  = rightCsf.occupation[i]
                if Nj != 0
                    ji  = Basics.subshell_j(shi)
                    Ji  = leftCsf.subshellJ[i];            Jj  = rightCsf.subshellJ[i]
                    si  = leftCsf.seniorityNr[i];          sj  = rightCsf.seniorityNr[i]
                    Qi  = qshellTermQ(ji,si);              Qj  = qshellTermQ(ji,sj)
                    aT  = SpinAngularGaigalas.getTermNumber(ji, Ni, Qi, Ji)
                    aTp = SpinAngularGaigalas.getTermNumber(ji, Nj, Qj, Jj)
                    if  leftCsf == rightCsf 
                        wa = irreducibleTensor(SchemeEta_W(),aT, Ni, AngularM64(1//2), AngularM64(-1//2), 0, aTp, Nj)
                    else
                        if recouplingCheck(leftCsf, rightCsf, i, i, i, i, length(subshells), 0) != 0.0
                            wa = irreducibleTensor(SchemeEta_W(),aT, Ni, AngularM64(1//2), AngularM64(-1//2), 0, aTp, Nj)
                        end
                    end
                    if abs(wa) >= 2.0e-11 
                        # Pure one-particle spin-angular coefficient
                        wa = -wa/sqrt(Basics.twice(Ji) + 1.)
                        # GRASP like spin-angular coefficient
                        wa = wa * sqrt(Basics.twice(ji) + 1.)
                        push!(coeffs1p, Coefficient1p(op.rank, shi, shi, wa))
                    end
                end
            end
        end
    elseif count == 2  &&  diff_occ == 2
        if creation != annihilation 
            shii = subshells[creation];               ja  = Basics.subshell_j(shii)
            shjj = subshells[annihilation];           jb  = Basics.subshell_j(shjj)
            if ja == jb  
                no_one = min(creation, annihilation);     no_two = max(creation, annihilation)
                if recouplingCheck(leftCsf, rightCsf, no_one, no_two, no_two, no_two, length(subshells), 1) != 0.0 
                    wa = recoupling1p(leftCsf, rightCsf, ja, jb, AngularJ64(0//2), no_one, no_two, length(subshells))
                    if abs(wa) >= 2.0e-11
                        aN  = leftCsf.occupation[creation];         aNp = rightCsf.occupation[creation]
                        Ji  = leftCsf.subshellJ[creation];          Jj  = rightCsf.subshellJ[creation]
                        si  = leftCsf.seniorityNr[creation];        sj  = rightCsf.seniorityNr[creation]
                        Qi  = qshellTermQ(ja,si);                   Qj  = qshellTermQ(ja,sj)
                        aT  = SpinAngularGaigalas.getTermNumber(ja, aN, Qi, Ji)
                        aTp = SpinAngularGaigalas.getTermNumber(ja, aNp, Qj, Jj)

                        bN  = leftCsf.occupation[annihilation];     bNp = rightCsf.occupation[annihilation]
                        Ji  = leftCsf.subshellJ[annihilation];      Jj  = rightCsf.subshellJ[annihilation]
                        si  = leftCsf.seniorityNr[annihilation];    sj  = rightCsf.seniorityNr[annihilation]
                        Qi  = qshellTermQ(jb,si);                   Qj  = qshellTermQ(jb,sj)
                        bT  = SpinAngularGaigalas.getTermNumber(jb, bN, Qi, Ji)
                        bTp = SpinAngularGaigalas.getTermNumber(jb, bNp, Qj, Jj)
                        wa = wa * SpinAngularGaigalas.irreducibleTensor(SchemeEta_a(),aT, aN, AngularM64(1//2), aTp, aNp) *
                                    SpinAngularGaigalas.irreducibleTensor(SchemeEta_a(),bT, bN,AngularM64(-1//2), bTp, bNp)
                        if abs(wa) >= 2.0e-11
                            occup = 0
                            for i = no_one:no_two-1    occup = occup + leftCsf.occupation[i]    end 
                            wa = (-1)^(occup+1) * wa 
                            #  Purer_one-particle spin-angular coefficient
                            wa = -wa
                            #  GRASP like spin-angular coefficient
                            wa = wa * sqrt(Basics.twice(ja) + 1.)
                            push!(coeffs1p, Coefficient1p(op.rank, shii, shjj, wa))
                        end
                    end
                end
            end
        end
    end
    return( coeffs1p )
end

"""
`SpinAngularGaigalas.computeCoefficientsScalar(op::SpinAngularGaigalas.TwoParticleOperator, leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1})`  
    ... computes the spin-angular coefficients for the reduced (scalar) two-particle matrix element

        <leftCsf || op^(0) || rightCsf >
        a Tuple of two lists with one- and two-particle coefficients tpl::Tuple{coeffs1p::Array{Coefficient1p,1}, coeffs2p::Array{Coefficient2p,1}}  is returned.

        leftCsf   - is the bra confuguration state function in relativistic notation;
        rightCsf  - is the ket confuguration state function in relativistic notation;
        subshells - is explicitly given relativistic subshell list if useStandardSubshells == false.
"""
function  computeCoefficientsScalar(op::SpinAngularGaigalas.TwoParticleOperator, leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1})
    coeffs2p = Coefficient2p[]
    if leftCsf.J != rightCsf.J return( coeffs2p ) end
    diff_occ = 0;   count = 0
    creation_one     = 0;   creation_two     = 0
    annihilation_one = 0;   annihilation_two = 0
    for  (i,shi)  in  enumerate(subshells)
        diff = leftCsf.occupation[i] - rightCsf.occupation[i]
        if abs(diff) > 2   return( coeffs2p )
        elseif diff != 0
            diff_occ = diff_occ + abs(diff);   count = count + 1
        end
        if diff == 1
            if creation_one == 0;            creation_one = i 
            elseif creation_two == 0;        creation_two = i
            else return( coeffs2p ) end
        elseif diff == 2
            if creation_one == 0;            creation_one = i;        creation_two = i
            else return( coeffs2p ) end
        elseif diff == -1
            if annihilation_one == 0;        annihilation_one = i
            elseif annihilation_two == 0;    annihilation_two = i
            else return( coeffs2p ) end
        elseif diff == -2
            if annihilation_one == 0;        annihilation_one = i;    annihilation_two = i
            else return( coeffs2p ) end
        end 
    end
    if count == 0 
        coeffs2p = twoParticleDiagonal(leftCsf, rightCsf, subshells)
    elseif count == 2 
        if diff_occ == 2
            coeffs2p = twoParticleDiffOcc2(leftCsf, rightCsf, creation_one, annihilation_one, subshells)
        elseif diff_occ == 4
            coeffs2p = twoParticle6(leftCsf, rightCsf, creation_one, annihilation_one, subshells)
        end
    elseif count == 3 
        coeffs2p = twoParticle15to18(leftCsf, rightCsf, creation_one,creation_two,annihilation_one,annihilation_two, subshells)
    elseif count == 4
        coeffs2p = twoParticle19to42(leftCsf, rightCsf, creation_one,creation_two,annihilation_one,annihilation_two, subshells)
    end
    return( coeffs2p )
end

"""
`SpinAngularGaigalas.irreducibleTensor(eta::SchemeEta_a, aIndex::Int64, aN::Int64, mq::AngularM64, bIndex::Int64, bN::Int64)`  
    ... computes the submatrix elements

        (j^N alpha Q J || a^(qj)_{mq} || j'^N' alpha' Q' J')
        as defined by Eq. (7) in Gaigalas et al. (CPC, 2001); a value::Float64 is returned.

        eta    - defines the irreducible tensoral form a^(qj)_mq;
        aIndex - is the state number of the bra subshell in quasispin representation;
        aN     - is number of electrons in the bra subshell;
        mq     - is quasispin projection mq;
        bIndex - is the state number of the ket subshell in quasispin representation;
        aN     - is number of electrons in the ket subshell.
"""
function  irreducibleTensor(eta::SchemeEta_a, aIndex::Int64, aN::Int64, mq::AngularM64, bIndex::Int64, bN::Int64)
    wa = 0.0
    if aIndex < 64  &&  bIndex < 64
        aT = qspaceTerms(aIndex);       bT = qspaceTerms(bIndex)
        if aT.min_odd !=  bT.min_even return(wa)  end
        if aN - bN - Basics.twice(mq) != 0 return(wa)  end
        if AngularMomentum.triangularDelta(aT.Q, AngularJ64(1//2), bT.Q) == 0 return(wa) end
        if AngularMomentum.triangularDelta(aT.J, aT.j, bT.J) == 0 return(wa) end
        aMQ = qshellTermM(aT.j, aN);              bMQ = qshellTermM(bT.j, bN)
        if qspacedelta(aT.Q, aMQ) == 0 return(wa)  end
        if qspacedelta(bT.Q, bMQ) == 0 return(wa)  end
        wa = - AngularMomentum.ClebschGordan(bT.Q, bMQ, AngularJ64(1//2), mq, aT.Q, aMQ)
        wa = wa * completlyReducedCfpByIndices(aIndex, bIndex)/sqrt(Basics.twice(aT.Q)+1.)
    elseif aIndex > 64  &&  bIndex > 64
        aJ = Int64(round(mod(aIndex,1000)))
        if mq.num == -1
            j  = Int64(round(mod(aIndex/(1000*1000),1000))) 
#             aQ = Int64(round(mod(aIndex/(1000),1000)))
            bJ = Int64(round(mod(bIndex,1000)))
            wa = (-1)^Int64((aJ-bJ-j+2*bN)/2) * sqrt(bN * (bJ + 1.0))
        elseif mq.num == 1    wa = (-1)^aN * sqrt(aN * (aJ + 1.0))
        end
    end
    return( wa )
end

"""
`SpinAngularGaigalas.irreducibleTensor(eta::SchemeEta_W, aIndex::Int64, aN::Int64, mLeft::AngularM64, mRight::AngularM64, kj::Int64, 
                                                    bIndex::Int64, bN::Int64)`  
    ... computes the submatrix elements

        (j^N alpha Q J || [a^(qj)_{mq1} x a^(qj)_{mq2}]^kj || j'^N' alpha' Q' J')
        as defined by Eq. (8) in Gaigalas et al. (CPC, 2001); a value::Float64 is returned.

        eta    - defines the irreducible tensoral form [a^(qj)_{mq1} x a^(qj)_{mq2}]^kj;
        aIndex - is the state number of the bra subshell in quasispin representation;
        aN     - is number of electrons in the bra subshell;
        mLeft  - is quasispin projection mq1;
        mRight - is quasispin projection mq2;
        bIndex - is the state number of the ket subshell in quasispin representation;
        aN     - is number of electrons in the ket subshell.
"""
function  irreducibleTensor(eta::SchemeEta_W, aIndex::Int64, aN::Int64, mLeft::AngularM64, mRight::AngularM64, kj::Int64, 
                            bIndex::Int64, bN::Int64)
    wa = 0.0
    if aIndex < 64  && bIndex < 64
        aT = qspaceTerms(aIndex);      bT = qspaceTerms(bIndex)
        if aT.min_even !=  bT.min_even  return( wa ) end
        if AngularMomentum.triangularDelta(aT.J, AngularJ64(2*kj//2), bT.J) == 0 return(wa) end
        aMQ = qshellTermM(aT.j, aN);              bMQ = qshellTermM(bT.j, bN)
        if qspacedelta(aT.Q, aMQ) == 0   return(wa)     end
        if qspacedelta(bT.Q, bMQ) == 0   return(wa)     end
        if aN - bN - Basics.twice(mLeft) - Basics.twice(mRight) != 0 return(wa)  end
        if mLeft ==  mRight
            #  cases    a * a    and    a^+ * a^+
            if (-1)^(kj) == 1
                if AngularMomentum.triangularDelta(aT.Q, AngularJ64(2//2), bT.Q) == 0 return(wa) end
                mtot_num = mLeft.num + mRight.num
                mtot = AngularM64(mtot_num//mLeft.den)
                wa = AngularMomentum.ClebschGordan(bT.Q, bMQ, AngularJ64(2//2), mtot, aT.Q, aMQ)
                if abs(wa) < 1.0e-11 return( wa )  end
                wa = wa * completelyReducedWkk(aIndex, bIndex, 1, kj)
                if abs(wa) < 1.0e-11 return( wa )  end
                wa = wa / sqrt(Basics.twice(aT.Q) + 1.)
            end
        else
            # cases    a * a^+    and    a^+ * a
            if kj == 0 
                if aIndex == bIndex
                    if mLeft  == AngularM64(1//2)
                        wa = - aN
                    else
                        wa = Basics.twice(aT.j)+1. - aN
                    end
                    wa = wa * sqrt((Basics.twice(aT.J)+1.)/(Basics.twice(aT.j)+1.))
                end
            else
                if (-1)^(kj) == 1  kq = 1
                else               kq = 0
                end
                if AngularMomentum.triangularDelta(aT.Q, AngularJ64((2*kq)//2), bT.Q) == 0    return(wa)    end
                mtot_num = mLeft.num + mRight.num
                mtot = AngularM64(mtot_num//mLeft.den)
                wa   = AngularMomentum.ClebschGordan(bT.Q, bMQ, AngularJ64((2*kq)//2), mtot, aT.Q, aMQ)
                if  abs(wa) < 1.0e-11 return( wa )  end
                wa = wa * completelyReducedWkk(aIndex, bIndex, kq, kj)
                if abs(wa) < 1.0e-11  return( wa )  end
                wa = wa / sqrt((Basics.twice(aT.Q) + 1.) * 2.)
                if mLeft == AngularM64(-1//2) && (-1)^(kj) == -1  wa = -wa  end
            end
        end
    elseif aIndex > 64  &&  bIndex > 64
        if aN > 2        error("SpinAngularGaigalas.irreducibleTensor SchemeEta_W: stop a")
        elseif bN > 2    error("SpinAngularGaigalas.irreducibleTensor SchemeEta_W: stop b")
        elseif mLeft.num + mRight.num + bN != aN    error("SpinAngularGaigalas.irreducibleTensor SchemeEta_W: stop c")
        end
        j  = Int64(round(mod(aIndex/(1000*1000),1000))) 
        aJ = Int64(round(mod(aIndex,1000)));    bJ = Int64(round(mod(bIndex,1000)))
        rN = bN + mRight.num
        if rN == 0        Q = Int64((j+1)/2);        minRun = 0;    maxRun = 0
        elseif rN == 1    Q = Int64(((j+1)/2)-1);    minRun = j;    maxRun = j
        elseif rN == 2                               minRun = 0;    maxRun = 2*j - 2
        elseif     error("SpinAngularGaigalas.irreducibleTensor SchemeEta_W: stop d")    
        end
        for runI = minRun:4:maxRun
            rJ = AngularJ64(runI//2)
            if rN == 2
                if rJ == 0    Q = Int64((j+1)/2)
                else          Q = Int64(((j+1)/2)-2)
                end
            end
            rIndex = ((j*1000) + Q)*1000 + runI
            wa =  wa + irreducibleTensor(SchemeEta_a(), aIndex, aN, mLeft, rIndex, rN) * 
                        irreducibleTensor(SchemeEta_a(), rIndex, rN, mRight, bIndex, bN) *
                        AngularMomentum.Wigner_6j(AngularJ64(j//2), AngularJ64(j//2), kj, AngularJ64(bJ//2), AngularJ64(aJ//2), rJ)
        end
        wa = (-1)^Int64((aJ+bJ+2*kj)/2) * wa * sqrt(2.0*kj + 1.0)
    end
    return( wa )
end

"""
`SpinAngularGaigalas.irreducibleTensor(eta::SchemeEta_WW, aIndex::Int64, aN::Int64, mLeft1::AngularM64, mRight1::AngularM64, 
                                mLeft2::AngularM64, mRight2::AngularM64, kj::Int64, bIndex::Int64, bN::Int64)`  
    ... computes the submatrix elements

        ( j^N QJ || [ W^(k_j) * W^(k_j)]^(0) || j^N' QJ)
        for subshells with j = 1/2, 3/2, 5/2, 7/2, 9/2 and for j > 9/2 with N = 0, 1, 2; a value::Float64 is returned.

        eta     - defines the irreducible tensoral form [ W^(k_j) * W^(k_j)]^(0);
        aIndex  - is the state number of the bra subshell in quasispin representation;
        aN      - is number of electrons in the bra subshell;
        mLeft1  - is quasispin projection mq1;
        mRight1 - is quasispin projection mq2;
        mLeft2  - is quasispin projection mq3;
        mRight2 - is quasispin projection mq4;
        bIndex  - is the state number of the ket subshell in quasispin representation;
        aN      - is number of electrons in the ket subshell.
"""
function  irreducibleTensor(eta::SchemeEta_WW, aIndex::Int64, aN::Int64, mLeft1::AngularM64, mRight1::AngularM64, 
                            mLeft2::AngularM64, mRight2::AngularM64, kj::Int64, bIndex::Int64, bN::Int64)
    wa = 0.0
    if aIndex < 64  && bIndex < 64
        aT = qspaceTerms(aIndex);          bT = qspaceTerms(bIndex)
        if aT.J != aT.J                    return( wa )  end
        if aT.min_even !=  bT.min_even     return( wa )  end
        if aT.max_even !=  bT.max_even     return( wa )  end
        aMQ = qshellTermM(aT.j, aN);       bMQ = qshellTermM(bT.j, bN)
        if qspacedelta(aT.Q, aMQ) == 0     return( wa )  end
        if qspacedelta(bT.Q, bMQ) == 0     return( wa )  end
        if aN-bN-Basics.twice(mLeft1)-Basics.twice(mRight1)-Basics.twice(mLeft2)-Basics.twice(mRight2) != 0 return( wa )  end
        for rIndex = aT.min_even:aT.max_even
            rT  = qspaceTerms(rIndex)
            rN  = Int64(bN+Basics.twice(mLeft2)+Basics.twice(mRight2))
            rMQ = qshellTermM(rT.j, rN)
            if qspacedelta(rT.Q, rMQ) != 0
                wa6j = AngularMomentum.Wigner_6j(kj, kj, 0, bT.J, aT.J, rT.J)
                if wa6j != 0.0
                    wairr = irreducibleTensor(SchemeEta_W(),aIndex, aN, mLeft1, mRight1, kj, rIndex, rN) *
                            irreducibleTensor(SchemeEta_W(),rIndex, rN, mLeft2, mRight2, kj, bIndex, bN)
                    wa = wa + (-1)^Int64((2*kj-Basics.twice(aT.J)+Basics.twice(rT.J))/2)* wairr
                end
            end
        end 
        wa = wa/sqrt((2*kj + 1.)*(Basics.twice(aT.J) + 1.))
    elseif aIndex > 64  &&  bIndex > 64
        if aIndex != bIndex    return( wa )  end
        if aN > 2        error("SpinAngularGaigalas.irreducibleTensor SchemeEta_WW: stop a")
        elseif bN > 2    error("SpinAngularGaigalas.irreducibleTensor SchemeEta_WW: stop b")
        elseif mLeft1.num + mRight1.num + mLeft2.num + mRight2.num + bN != aN    error("SpinAngularGaigalas.irreducibleTensori SchemeEta_WW: stop c")    
        end
        j  = Int64(round(mod(aIndex/(1000*1000),1000))) 
        aJ = Int64(round(mod(aIndex,1000)));         bJ = Int64(round(mod(bIndex,1000)))
        rN = bN + mLeft2.num + mRight2.num
        if rN == 0        Q = Int64((j+1)/2);        minRun = 0;    maxRun = 0
        elseif rN == 1    Q = Int64(((j+1)/2)-1);    minRun = j;    maxRun = j
        elseif rN == 2                               minRun = 0;    maxRun = 2*j - 2
        elseif     error("SpinAngularGaigalas.irreducibleTensor  SchemeEta_WW: stop d")
        end
        for runI = minRun:4:maxRun
            rJ = AngularJ64(runI//2)
            if rN == 2
                if rJ == 0    Q = Int64((j+1)/2)
                else          Q = Int64(((j+1)/2)-2)
                end
            end
            rIndex = ((j*1000) + Q)*1000 + runI
            wa =  wa + irreducibleTensor(SchemeEta_W(),aIndex, aN, mLeft1, mRight1, kj, rIndex, rN) *
                        irreducibleTensor(SchemeEta_W(),rIndex, rN, mLeft2, mRight2, kj, bIndex, bN) *
                        AngularMomentum.Wigner_6j(kj, kj, 0, AngularJ64(bJ//2), AngularJ64(aJ//2), rJ)
        end
        wa = (-1)^Int64((aJ+bJ)/2) * wa
    end
    return( wa )
end

"""
`SpinAngularGaigalas.irreducibleTensor(eta::SchemeEta_aW,aIndex::Int64, aN::Int64, mLeft1::AngularM64, mLeft2::AngularM64, 
                                mRight2::AngularM64, kj1::Int64, kj::AngularJ64, bIndex::Int64, bN::Int64)`  
    ... computes the submatrix elements

        ( j^N QJ || [ a^(q j)_mq * W^(kj1)]^(kj) || j^N' QJ)
        for subshells with j = 1/2, 3/2, 5/2, 7/2, 9/2 and for j > 9/2 with N = 0, 1, 2; a value::Float64 is returned.

        eta     - defines the irreducible tensoral form [ a^(q j)_mq * W^(kj1)]^(kj);
        aIndex  - is the state number of the bra subshell in quasispin representation;
        aN      - is number of electrons in the bra subshell;
        mLeft1  - is quasispin projection mq;
        mLeft2  - is quasispin projection mq3;
        mRight2 - is quasispin projection mq4;
        bIndex  - is the state number of the ket subshell in quasispin representation;
        aN      - is number of electrons in the ket subshell.
"""
function  irreducibleTensor(eta::SchemeEta_aW,aIndex::Int64, aN::Int64, mLeft1::AngularM64, mLeft2::AngularM64, 
                            mRight2::AngularM64, kj1::Int64, kj::AngularJ64, bIndex::Int64, bN::Int64)
    wa = 0.0;    kkj1 = 2* kj1;   kj1a =AngularJ64(kkj1//2)
    if aIndex < 64  &&  bIndex < 64
            aT = qspaceTerms(aIndex);          bT = qspaceTerms(bIndex)
            if AngularMomentum.triangularDelta(aT.J, kj, bT.J) == 0 return(wa) end
            if aT.min_odd !=  bT.min_even      return( wa )    end
            if aT.max_odd !=  bT.max_even      return( wa )    end
            aMQ = qshellTermM(aT.j, aN);       bMQ = qshellTermM(bT.j, bN)
            if qspacedelta(aT.Q, aMQ) == 0     return( wa )    end
            if qspacedelta(bT.Q, bMQ) == 0     return( wa )    end
            if aN-bN-Basics.twice(mLeft1)-Basics.twice(mLeft2)-Basics.twice(mRight2) != 0 return( wa )  end
            for rIndex=aT.min_odd:aT.max_odd
                    rT = qspaceTerms(rIndex)
                    rN = Int64(bN+Basics.twice(mLeft2)+Basics.twice(mRight2))
                    rMQ = qshellTermM(rT.j, rN)
                    if qspacedelta(rT.Q, rMQ) != 0
                    wa6j = AngularMomentum.Wigner_6j(aT.j, kj1a, kj, bT.J, aT.J, rT.J)
                    if wa6j != 0.0
                        wb =  irreducibleTensor(SchemeEta_a(), aIndex, aN, mLeft1, rIndex, rN)
                        if abs(wb) >= 2.0e-10
                            wb = wb * wa6j * irreducibleTensor(SchemeEta_W(), rIndex, rN, mLeft2, mRight2, kj1, bIndex, bN)
                            wa = wa + wb
                        end
                    end
                end
            end
            wa = wa * sqrt(Basics.twice(kj) + 1.0)
            wa = (-1)^Int64((Basics.twice(kj) + Basics.twice(aT.J) + Basics.twice(bT.J))/2) * wa
    elseif aIndex > 64  &&  bIndex > 64
        if aN > 2        error("SpinAngularGaigalas.irreducibleTensor SchemeEta_aW: stop a")
        elseif bN > 2    error("SpinAngularGaigalas.irreducibleTensor SchemeEta_aW: stop b")
        elseif mLeft1.num + mLeft2.num + mRight2.num + bN != aN    error("SpinAngularGaigalas.irreducibleTensori SchemeEta_aW: stop c")    
        end
        j  = Int64(round(mod(aIndex/(1000*1000),1000))) 
        aJ = Int64(round(mod(aIndex,1000)));         bJ = Int64(round(mod(bIndex,1000)))
        rN = bN + mLeft2.num + mRight2.num
        if rN == 0        Q = Int64((j+1)/2);        minRun = 0;    maxRun = 0
        elseif rN == 1    Q = Int64(((j+1)/2)-1);    minRun = j;    maxRun = j
        elseif rN == 2                               minRun = 0;    maxRun = 2*j - 2
        elseif     error("SpinAngularGaigalas.irreducibleTensor SchemeEta_aW: stop d")
        end
        for runI = minRun:4:maxRun
            rJ = AngularJ64(runI//2)
            if rN == 2
                if rJ == 0    Q = Int64((j+1)/2)
                else          Q = Int64(((j+1)/2)-2)
                end
            end
            rIndex = ((j*1000) + Q)*1000 + runI
            wa =  wa + irreducibleTensor(SchemeEta_a(), aIndex, aN, mLeft1, rIndex, rN) *
                        irreducibleTensor(SchemeEta_W(),rIndex, rN, mLeft2, mRight2, kj1, bIndex, bN) *
                        AngularMomentum.Wigner_6j(AngularJ64(j//2), kj1a, kj, AngularJ64(bJ//2), AngularJ64(aJ//2), rJ)
        end
        wa = (-1)^Int64((aJ+bJ+Basics.twice(kj))/2) * sqrt(Basics.twice(kj) + 1.0) * wa
    end
    return( wa )
end

"""
`SpinAngularGaigalas.irreducibleTensor(eta::SchemeEta_Wa, aIndex::Int64, aN::Int64, mLeft1::AngularM64, mRight1::AngularM64,  
                                mLeft2::AngularM64, kj1::Int64, kj::AngularJ64, bIndex::Int64, bN::Int64)`  
    ... computes the submatrix elements

        ( j^N QJ || [ W^(kj1) * a^(q j)_mq ]^(kj) || j^N' QJ)
        for subshells with j = 1/2, 3/2, 5/2, 7/2, 9/2 and for j > 9/2 with N = 0, 1, 2; a value::Float64 is returned.

        eta     - defines the irreducible tensoral form [ W^(kj1) * a^(q j)_mq ]^(kj);
        aIndex  - is the state number of the bra subshell in quasispin representation;
        aN      - is number of electrons in the bra subshell;
        mLeft1  - is quasispin projection mq1;
        mRight1 - is quasispin projection mq2;
        mRight2 - is quasispin projection mq;
        bIndex  - is the state number of the ket subshell in quasispin representation;
        aN      - is number of electrons in the ket subshell.
"""
function  irreducibleTensor(eta::SchemeEta_Wa, aIndex::Int64, aN::Int64, mLeft1::AngularM64, mRight1::AngularM64, 
                            mLeft2::AngularM64, kj1::Int64, kj::AngularJ64, bIndex::Int64, bN::Int64)
    wa = 0.0;    kkj1 = 2* kj1;   kj1a =AngularJ64(kkj1//2)
    if aIndex < 64  && bIndex < 64
        aT = qspaceTerms(aIndex);          bT = qspaceTerms(bIndex)
        if AngularMomentum.triangularDelta(aT.J, kj, bT.J) == 0 return(wa) end
        if aT.min_even !=  bT.min_odd      return( wa )    end
        if aT.max_even !=  bT.max_odd      return( wa )    end
        aMQ = qshellTermM(aT.j, aN);       bMQ = qshellTermM(bT.j, bN)
        if qspacedelta(aT.Q, aMQ) == 0     return( wa )    end
        if qspacedelta(bT.Q, bMQ) == 0     return( wa )    end
        if aN-bN-Basics.twice(mLeft1)-Basics.twice(mRight1)-Basics.twice(mLeft2) != 0 return( wa )  end
        for rIndex = aT.min_even:aT.max_even
            rT = qspaceTerms(rIndex)
            rN = Int64(bN+Basics.twice(mLeft2))
            rMQ = qshellTermM(rT.j, rN)
            if qspacedelta(rT.Q, rMQ) != 0
                wa6j = AngularMomentum.Wigner_6j(kj1a, aT.j, kj, bT.J, aT.J, rT.J)
                if wa6j != 0.0
                    wb =  irreducibleTensor(SchemeEta_a(), rIndex, rN, mLeft2, bIndex, bN)
                    if abs(wb) >= 2.0e-10
                        wb = wb * wa6j * irreducibleTensor(SchemeEta_W(), aIndex, aN, mLeft1, mRight1, kj1, rIndex, rN)
                        wa = wa + wb
                    end
                end
            end
        end
        wa = wa * sqrt(Basics.twice(kj) + 1.0)
        wa = (-1)^Int64((Basics.twice(kj) + Basics.twice(aT.J) + Basics.twice(bT.J))/2) * wa
    elseif aIndex > 64  &&  bIndex > 64
        if aN > 2        error("SpinAngularGaigalas.irreducibleTensor SchemeEta_Wa: stop a")
        elseif bN > 2    error("SpinAngularGaigalas.irreducibleTensor SchemeEta_Wa: stop b")
        elseif mLeft1.num + mRight1.num + mLeft2.num + bN != aN    error("SpinAngularGaigalas.irreducibleTensor SchemeEta_Wa: stop c")
        end
        j  = Int64(round(mod(aIndex/(1000*1000),1000))) 
        aJ = Int64(round(mod(aIndex,1000)));         bJ = Int64(round(mod(bIndex,1000)))
        rN = bN + mLeft2.num
        if rN == 0        Q = Int64((j+1)/2);        minRun = 0;    maxRun = 0
        elseif rN == 1    Q = Int64(((j+1)/2)-1);    minRun = j;    maxRun = j
        elseif rN == 2                               minRun = 0;    maxRun = 2*j - 2
        elseif     error("SpinAngularGaigalas.irreducibleTensor SchemeEta_Wa: stop d")
        end
        for runI = minRun:4:maxRun
            rJ = AngularJ64(runI//2)
            if rN == 2
                if rJ == 0    Q = Int64((j+1)/2)
                else          Q = Int64(((j+1)/2)-2)
                end
            end
            rIndex = ((j*1000) + Q)*1000 + runI
            wa =  wa + irreducibleTensor(SchemeEta_W(),aIndex, aN, mLeft1, mRight1, kj1, rIndex, rN) *
                        irreducibleTensor(SchemeEta_a(), rIndex, rN, mLeft2, bIndex, bN) *
                        AngularMomentum.Wigner_6j(kj1a, AngularJ64(j//2), kj, AngularJ64(bJ//2), AngularJ64(aJ//2), rJ)
        end
        wa = (-1)^Int64((aJ+bJ+Basics.twice(kj))/2) * sqrt(Basics.twice(kj) + 1.0) * wa
    end
    return( wa )
end

"""
`SpinAngularGaigalas.normalForm(ia::Int64, ib::Int64, ic::Int64)`
    ... odering operators of second quantisation in norma form; a value::Int64 is returned.
    
    ia, ib, ic - is the place of operators in the set of operators of second quantisatio.
"""
function normalForm(ia::Int64, ib::Int64, ic::Int64)
    wa = Array{Int64}(undef,3);    for i=1:3  wa[i] = (0)   end
    wa[1] = ia;     wa[3] = ia
    if wa[1] > ib   wa[1] = ib  end
    if wa[3] < ib   wa[3] = ib  end
    if wa[1] > ic   wa[1] = ic  end
    if wa[3] < ic   wa[3] = ic  end
    if ia  > wa[1]  &&  ia  < wa[3]   wa[2] = ia  end
    if ib  > wa[1]  &&  ib  < wa[3]   wa[2] = ib  end
    if ic  > wa[1]  &&  ic  < wa[3]   wa[2] = ic  end
    return( wa )
end

"""
`SpinAngularGaigalas.normalPhase(ia::Int64, ib::Int64, ic::Int64, id::Int64)`
    ... compute the phase factor from odering operators of second quantisation in norma form; a wa::Int64 is returned.

    ia, ib, ic, id - is the place of operators in the set of operators of second quantisatio.
"""
function normalPhase(ia::Int64, ib::Int64, ic::Int64, id::Int64)
    wa = 1
    if ia > ib  wa = -wa   end
    if ia > ic  wa = -wa   end
    if ia > id  wa = -wa   end
    if ib > ic  wa = -wa   end
    if ib > id  wa = -wa   end
    if ic > id  wa = -wa   end
    return( wa )
end




"""
`SpinAngularGaigalas.recouplingDiagram(diagram::DiagramC01,leftCsf::CsfR, rightCsf::CsfR, rank::AngularJ64, ia::Int64)`  
    ... computes the coefficient C_1 from the paper of G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, 
        Vol 30 3747, Eq. (15).  The phase factor is given in the table 2 in this reference; a value::Float64 is returned.

        leftCsf   - is the bra confuguration state function in relativistic notation;
        rightCsf  - is the ket confuguration state function in relativistic notation;
        rank      - is the rank of operator;
        ia        - is the index of interacting subshell.
"""
function recouplingDiagram(diagram::DiagramC01,leftCsf::CsfR, rightCsf::CsfR, rank::AngularJ64, ia::Int64)
    wa = 0.0
    csf_r_J = leftCsf.subshellJ[ia];         csf_s_J = rightCsf.subshellJ[ia]
    if ia == 1
        csf_r_X = leftCsf.subshellX[2];      csf_s_X = rightCsf.subshellX[2];    J = leftCsf.subshellJ[2]
    elseif ia == 2
        csf_r_X = leftCsf.subshellX[2];      csf_s_X = rightCsf.subshellX[2];    J = leftCsf.subshellJ[1]
    else
        csf_r_X = leftCsf.subshellX[ia];     csf_s_X = rightCsf.subshellX[ia];   J = leftCsf.subshellX[ia-1]
    end
    wa = AngularMomentum.Wigner_6j(rank, csf_s_J, csf_r_J, J, csf_r_X, csf_s_X)
    wa = wa * sqrt((Basics.twice(csf_r_J) + 1.) * (Basics.twice(csf_s_X) + 1.))

    wa = (-1)^Int64((Basics.twice(J)+Basics.twice(csf_r_X)+Basics.twice(csf_s_J)+Basics.twice(rank))/2) * wa
    if ia == 1
        wa = (-1)^Int64((Basics.twice(csf_r_J)+Basics.twice(csf_s_J)+2*Basics.twice(J)-Basics.twice(csf_r_X)-
                            Basics.twice(csf_s_X))/2) * wa
    end
    return( wa )
end

"""
`SpinAngularGaigalas.recouplingDiagram(diagram::DiagramC02, leftCsf::CsfR, rightCsf::CsfR, rank::AngularJ64, ia::Int64, ib::Int64)`  
    ... computes the coefficient C_2 from the paper of G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, 
        Vol 30 3747, Eq. (16).  The phase factor is given in the table 2 in this reference; a value::Float64 is returned.

        leftCsf   - is the bra confuguration state function in relativistic notation;
        rightCsf  - is the ket confuguration state function in relativistic notation;
        rank      - is the rank of operator;
        ia        - is the index of first interacting subshell;
        ib        - is the index of second interacting subshell.
"""
function recouplingDiagram(diagram::DiagramC02, leftCsf::CsfR, rightCsf::CsfR, rank::AngularJ64, ia::Int64, ib::Int64)
    wa = 1.0
    if ia == 1 run = 3
    else  run = ia + 1
    end
    if run < ib
        for i = run:ib-1
            csf_J = leftCsf.subshellJ[i]
            csf_r_X = leftCsf.subshellX[i-1];     csf_s_X = rightCsf.subshellX[i-1]
            csf_r_X_1 = leftCsf.subshellX[i];     csf_s_X_1 = rightCsf.subshellX[i]
            wa_run  = AngularMomentum.Wigner_6j(rank, csf_s_X, csf_r_X, csf_J, csf_r_X_1, csf_s_X_1)
            wa_run = wa_run * sqrt((Basics.twice(csf_r_X) + 1.) * (Basics.twice(csf_s_X_1) + 1.))
            wa_run = (-1)^Int64((Basics.twice(rank)+Basics.twice(csf_J)+Basics.twice(csf_r_X)+Basics.twice(csf_s_X_1))/2) * 
                        wa_run
            wa = wa * wa_run
        end
    end
    return( wa )
end

"""
`SpinAngularGaigalas.recouplingDiagram(diagram::DiagramC03, leftCsf::CsfR, rightCsf::CsfR, rank::AngularJ64, ia::Int64, nwshells::Int64)`
    ... computes the coefficient C_3 from the paper of G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, 
        Vol 30 3747, Eq. (17).  The phase factor is given in the table 3 in this reference; a value::Float64 is returned.

        leftCsf  - is the bra confuguration state function in relativistic notation;
        rightCsf - is the ket confuguration state function in relativistic notation;
        rank     - is the rank of operator;
        ia       - is the index of interacting subshell;
        nwshells - is number of subshells in confuguration state function.
"""
function recouplingDiagram(diagram::DiagramC03, leftCsf::CsfR, rightCsf::CsfR, rank::AngularJ64, ia::Int64, nwshells::Int64)
    wa = 0.0
    T1_r = leftCsf.subshellX[nwshells];        T1_s = rightCsf.subshellX[nwshells]
    if ia == nwshells
        T_r = leftCsf.subshellJ[nwshells];     T_s = rightCsf.subshellJ[nwshells];    JI = leftCsf.subshellX[nwshells-1]
    else
        T_r = leftCsf.subshellX[nwshells-1];   T_s = rightCsf.subshellX[nwshells-1];  JI = leftCsf.subshellJ[nwshells]
    end
    wa = AngularMomentum.Wigner_6j(rank, T_s, T_r, JI, T1_r, T1_s)
    wa = wa * sqrt((Basics.twice(T_r) + 1.) * (Basics.twice(T1_s) + 1.))
    wa = (-1)^Int64((Basics.twice(rank)+Basics.twice(JI)+Basics.twice(T_s)+Basics.twice(T1_r))/2) * wa
    if ia == nwshells return( wa ) end
    wa = (-1)^Int64((Basics.twice(T_r)+Basics.twice(T_s)-Basics.twice(T1_s)-Basics.twice(T1_r)+Basics.twice(JI)*2)/2) * wa
    return( wa )
end

"""
`SpinAngularGaigalas.recouplingDiagram(diagram::DiagramC04, leftCsf::CsfR, rightCsf::CsfR, rank_1::AngularJ64, rank_2::AngularJ64, 
                                rank::AngularJ64, ia::Int64, ib::Int64)`  
    ... computes the coefficient C_4 from the paper of G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, 
        Vol 30 3747, Eq. (23).  The phase factor is given in the table 4 in this reference; a value::Float64 is returned.

        leftCsf  - is the bra confuguration state function in relativistic notation;
        rightCsf - is the ket confuguration state function in relativistic notation;
        rank_1   - is the first rank of operator;
        rank_2   - is the second rank of operator;
        ia       - is the index of first interacting subshell;
        ib       - is the index of second interacting subshell.
"""
function recouplingDiagram(diagram::DiagramC04, leftCsf::CsfR, rightCsf::CsfR, rank_1::AngularJ64, rank_2::AngularJ64, 
                            rank::AngularJ64, ia::Int64, ib::Int64)
    wa = 0.0
    csf_r_J_2 = leftCsf.subshellJ[ib];        csf_s_J_2 = rightCsf.subshellJ[ib]
    csf_r_X   = leftCsf.subshellX[ib];        csf_s_X   = rightCsf.subshellX[ib]
    if ia == 1  &&  ib == 2     csf_r_X_1 = leftCsf.subshellJ[ia];    csf_s_X_1 = rightCsf.subshellJ[ia]
    else                        csf_r_X_1 = leftCsf.subshellX[ib-1];  csf_s_X_1 = rightCsf.subshellX[ib-1]
    end
    wa = AngularMomentum.Wigner_9j(csf_s_X_1, rank_1, csf_r_X_1, csf_s_J_2, rank_2, csf_r_J_2, csf_s_X, rank, csf_r_X)
    wa = wa * sqrt((Basics.twice(csf_r_X_1) + 1.)*(Basics.twice(rank) + 1.)*(Basics.twice(csf_r_J_2) + 1.) * 
                    (Basics.twice(csf_s_X) + 1.))
    return( wa )
end

"""
`SpinAngularGaigalas.recouplingDiagram(diagram::DiagramC05, leftCsf::CsfR, rightCsf::CsfR, rank::AngularJ64, ia::Int64, ib::Int64)`  
    ... computes the coefficient C_5 from the paper of G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, 
        Vol 30 3747, Eq. (23).  The phase factor is given in the table 5 in this reference; a value::Float64 is returned.

        leftCsf  - is the bra confuguration state function in relativistic notation;
        rightCsf - is the ket confuguration state function in relativistic notation;
        rank     - is the rank of operator;
        ia       - is the index of first interacting subshell;
        ib       - is the index of second interacting subshell.
"""
function recouplingDiagram(diagram::DiagramC05, leftCsf::CsfR, rightCsf::CsfR, rank::AngularJ64, ia::Int64, ib::Int64)
    wa = 0.0
    csf_r_J_2 = leftCsf.subshellJ[ib];      csf_s_J_2 = rightCsf.subshellJ[ib]
    csf_X_2   = leftCsf.subshellX[ib]
    if ia == 1  &&  ib == 2     csf_r_J = leftCsf.subshellJ[ia];    csf_s_J = rightCsf.subshellJ[ia]
    else                        csf_r_J = leftCsf.subshellX[ib-1];  csf_s_J = rightCsf.subshellX[ib-1]
    end
    wa = AngularMomentum.Wigner_6j(rank, csf_s_J_2, csf_r_J_2, csf_X_2, csf_r_J, csf_s_J)
    wa = wa * sqrt((Basics.twice(csf_s_J_2) + 1.) * (Basics.twice(csf_r_J) + 1.))
    wa = (-1)^Int64((Basics.twice(csf_X_2)+Basics.twice(csf_s_J)+Basics.twice(csf_r_J_2)+Basics.twice(rank))/2) * wa
    return( wa )
end

"""
`SpinAngularGaigalas.recouplingCheck(leftCsf::CsfR, rightCsf::CsfR, ia::Int64, ib::Int64, ic::Int64, id::Int64, 
                                nwshells::Int64, switch::Int64)`  
    ... checks the angular momentum selection rules for the recoupling coefficients in the cases of one, two, 
        three or four interacting shells.; a value::Float64 is returned.

        leftCsf  - is the bra confuguration state function in relativistic notation;
        rightCsf - is the ket confuguration state function in relativistic notation;
        ia       - is the index of first interacting subshell;
        ib       - is the index of second interacting subshell;
        ic       - is the index of third interacting subshell;
        id       - is the index of fourth interacting subshell;
        nwshells - is number of subshells in confuguration state function;
        switch   - parameter for selecting the checks.
"""
function recouplingCheck(leftCsf::CsfR, rightCsf::CsfR, ia::Int64, ib::Int64, ic::Int64, id::Int64, nwshells::Int64, switch::Int64)
    wa = 1.0
    if nwshells == 1 return( wa ) end
    if switch == 0 
        for i = 1:nwshells
            if i > 1
                if leftCsf.subshellX[i]   != rightCsf.subshellX[i]    wa = 0.0      end
            end
            if leftCsf.subshellJ[i]   != rightCsf.subshellJ[i]        wa = 0.0      end
            if leftCsf.seniorityNr[i] != rightCsf.seniorityNr[i]      wa = 0.0      end
            if wa == 0.0     return( wa)  end
            end
    elseif ia == 1  &&  ib == 2
        for i = 1:nwshells
            if i > 1
                if leftCsf.subshellX[i] != rightCsf.subshellX[i]        wa = 0.0  end
            end
            if i > 2
                if leftCsf.subshellJ[i] != rightCsf.subshellJ[i]        wa = 0.0  end
                if leftCsf.seniorityNr[i] != rightCsf.seniorityNr[i]    wa = 0.0  end
            end
            if wa == 0.0 return( wa)  end
        end
    else
        if ia == 2  ix = 1    else        ix = ia    end
        for i = 1:nwshells
            if i < ix  ||  i >= ib 
                if i > 1
                    if leftCsf.subshellX[i] != rightCsf.subshellX[i]    wa = 0.0  end
                end
            end
            if      i == ia
            elseif  i == ib 
            elseif  switch == 2 && i == ic
            elseif  switch == 3 && i == ic
            elseif  switch == 3 && i == id
            else
                if leftCsf.subshellJ[i]   != rightCsf.subshellJ[i]      wa = 0.0  end
                if leftCsf.seniorityNr[i] != rightCsf.seniorityNr[i]    wa = 0.0  end
            end
            if wa == 0.0 return( wa )  end
        end
    end 
    return( wa )
end

"""
`SpinAngularGaigalas.recouplingCheck(leftCsf::CsfR, rightCsf::CsfR, rank::AngularJ64, ia::Int64, ib::Int64, nwshells::Int64)`  
    ... checks the angular momentum selection rules for the recoupling coefficients in the cases of one, 
        two, three or four interacting shells.; a value::Float64 is returned.

        leftCsf  - is the bra confuguration state function in relativistic notation;
        rightCsf - is the ket confuguration state function in relativistic notation;
        rank     - is the rank of operator;
        ia       - is the index of first interacting subshell;
        ib       - is the index of second interacting subshell;
        nwshells - is number of subshells in confuguration state function.
"""
function recouplingCheck(leftCsf::CsfR, rightCsf::CsfR, rank::AngularJ64, ia::Int64, ib ::Int64, nwshells::Int64)
    wa = 0.0
    if AngularMomentum.triangularDelta(leftCsf.J, rank, rightCsf.J) == 0    return(wa)  end
    wa = 1.0
    if nwshells == 1 return( wa ) end
    for i = 1:nwshells
        if i != ia
            if i != ib
                if leftCsf.subshellJ[i] != rightCsf.subshellJ[i]        wa = 0.0    end
                if leftCsf.seniorityNr[i] != rightCsf.seniorityNr[i]    wa = 0.0    end
            end
        end
        if wa == 0.0 return( wa )  end
    end
    if nwshells <= 2 return( wa ) end
    if nwshells-ib > 0
        for i = ib:nwshells
            if AngularMomentum.triangularDelta(leftCsf.subshellX[i], rank, rightCsf.subshellX[i]) == 0 wa = 0.0 end
        end
    end
    if wa == 0.0 return( wa ) end
    for i = 1:ia-1   if leftCsf.subshellX[i] != rightCsf.subshellX[i]    wa = 0.0  end   end
    return( wa )
end

"""
`SpinAngularGaigalas.recoupling1p(leftCsf::CsfR, rightCsf::CsfR, rank::AngularJ64, ia::Int64, nwshells::Int64)`
    ... computes the recoupling matrix R(ja, jb, Lambda^bra, Lambda^ket) for a non scalar operator in the case 
        of one interacting shells. See G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747, 
        Eq. (14) (p. 3756); a value::Float64 is returned.

        leftCsf  - is the bra confuguration state function in relativistic notation;
        rightCsf - is the ket confuguration state function in relativistic notation;
        rank     - is the rank of operator;
        ia       - is the index of interacting subshell;
        nwshells - is number of subshells in confuguration state function.
"""
function  recoupling1p(leftCsf::CsfR, rightCsf::CsfR, rank::AngularJ64, ia::Int64, nwshells::Int64)
    wa = 1. / sqrt(Basics.twice(leftCsf.subshellJ[ia])+1.)
    if nwshells == 1    return( wa )    end
    if nwshells != 2 
        wa  = wa * recouplingDiagram(DiagramC03(), leftCsf, rightCsf, rank, ia, nwshells)
        if ia == nwshells return( wa ) end
    end
    wa = wa * recouplingDiagram(DiagramC01(), leftCsf, rightCsf, rank, ia)
    run_correction = nwshells - ia
    if ia == 1      run_correction = run_correction - 1 end
    if run_correction <= 1      return( wa )    end
    wa = wa * recouplingDiagram(DiagramC02(), leftCsf, rightCsf, rank, ia, nwshells)
    return( wa )
end

"""
`SpinAngularGaigalas.recoupling1p(leftCsf::CsfR, rightCsf::CsfR, rank_1::AngularJ64, rank_2::AngularJ64, 
                            rank::AngularJ64, ia::Int64, ib::Int64, nwshells::Int64)`
    ... computes the recoupling matrix R(ja, jb, Lambda^bra, Lambda^ket) for a non scalar operator in the case of two 
        interacting shells. See G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747, Eq. (19) (p. 3756);
        a value::Float64 is returned.

        leftCsf  - is the bra confuguration state function in relativistic notation;
        rightCsf - is the ket confuguration state function in relativistic notation;
        rank_1   - is the first rank of operator;
        rank_2   - is the second rank of operator;
        rank     - is the final rank of operator;
        ia       - is the index of the first interacting subshell;
        ib       - is the index of the second interacting subshell;
        nwshells - is number of subshells in confuguration state function.
"""
function  recoupling1p(leftCsf::CsfR, rightCsf::CsfR, rank_1::AngularJ64, rank_2::AngularJ64, 
                        rank::AngularJ64, ia::Int64, ib::Int64, nwshells::Int64)
    wa = 1. / sqrt((Basics.twice(leftCsf.subshellJ[ia])+1.) * (Basics.twice(leftCsf.subshellJ[ib])+1.))
    if nwshells - ib > 1
        wa = wa * recouplingDiagram(DiagramC02(), leftCsf,rightCsf,rank,ib,nwshells)
        if abs(wa) < 1.0e-11   return ( wa )   end
    end 
    if ib != nwshells
        wa = wa *  recouplingDiagram(DiagramC03(), leftCsf, rightCsf, rank, ib, nwshells)
        if abs(wa) < 1.0e-11   return ( wa )   end
    end
    wa = wa * recouplingDiagram(DiagramC04(), leftCsf, rightCsf, rank_1, rank_2, rank, ia, ib)
    if abs(wa) < 1.0e-11   return ( wa )   end
    if ia == 1 && ib == 2 return( wa) end
    wa = wa * recouplingDiagram(DiagramC01(), leftCsf, rightCsf, rank_1, ia)
    if abs(wa) < 1.0e-11   return ( wa )   end
    run_correction = ib - ia
    if ia == 1  run_correction = run_correction - 1  end
    if run_correction <= 1  return( wa ) end
    wa = wa * recouplingDiagram(DiagramC02(), leftCsf, rightCsf, rank_1, ia, ib)
    return( wa )
end

"""
`SpinAngularGaigalas.recoupling2p(leftCsf::CsfR, rightCsf::CsfR, rank::AngularJ64, ia::Int64, ib::Int64)`  
    ... computes the recoupling matrix R(ja, jb, Lambda^bra, Lambda^ket) for a scalar operator in the case of two 
        interacting shells. See G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747, Eq. (22) (p. 3756);
        a value::Float64 is returned.

        leftCsf  - is the bra confuguration state function in relativistic notation;
        rightCsf - is the ket confuguration state function in relativistic notation;
        rank     - is the rank of operator;
        ia       - is the index of the first interacting subshell.
        ib       - is the index of the second interacting subshell.
"""
function  recoupling2p(leftCsf::CsfR, rightCsf::CsfR, rank::AngularJ64, ia::Int64, ib::Int64) 
    wa = 1. / sqrt((Basics.twice(leftCsf.subshellJ[ia])+1.) * (Basics.twice(leftCsf.subshellJ[ib])+1.))
    if rank != 0 
        wa = wa * recouplingDiagram(DiagramC05(), leftCsf, rightCsf, rank, ia, ib)
        if abs(wa) < 1.0e-11   return ( wa )   end
        wa = wa * sqrt((Basics.twice(leftCsf.subshellJ[ib]) + 1.) / 
                    ((Basics.twice(rank) + 1.) * (Basics.twice(rightCsf.subshellJ[ib]) + 1.)))
        if ia == 1 && ib == 2  return( wa )  end
        wa = wa * recouplingDiagram(DiagramC01(), leftCsf, rightCsf, rank, ia)
        if abs(wa) < 1.0e-11   return ( wa )   end
        if      ia == 1;         if ib - 1 - ia <= 1    return( wa )    end
        else;   if ib - ia <= 1                         return( wa )    end
        end 
        wa = wa * recouplingDiagram(DiagramC02(), leftCsf, rightCsf, rank, ia, ib)
    end
    return( wa )
end

"""
`SpinAngularGaigalas.recoupling2p(leftCsf::CsfR, rightCsf::CsfR, rank1::AngularJ64, rank2::AngularJ64, rank3::AngularJ64, 
                            ia::Int64, ib::Int64, ic::Int64)`  
    ... computes the recoupling matrix R(ja, jb, Lambda^bra, Lambda^ket) for a scalar operator in the case of three 
        interacting shells  (shells are not ordered). See G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys,
        Vol 30 3747, Eq. (26) (p. 3756); a value::Float64 is returned.

        leftCsf  - is the bra confuguration state function in relativistic notation;
        rightCsf - is the ket confuguration state function in relativistic notation;
        rank_1   - is the first rank of operator;
        rank_2   - is the second rank of operator;
        rank_3   - is the third rank of operator;
        ia       - is the index of the first interacting subshell;
        ib       - is the index of the second interacting subshell;
        ic       - is the index of the third interacting subshell.
"""
function  recoupling2p(leftCsf::CsfR, rightCsf::CsfR, rank1::AngularJ64,  rank2::AngularJ64, rank3::AngularJ64, 
                        ia::Int64, ib::Int64, ic::Int64) 
    wa = 0.0
    if ic > ia  && ic > ib
        if      ia < ib     wa = recoupling2p(Ordered2p(), leftCsf, rightCsf, rank1, rank2, rank3, ia, ib, ic)
        elseif  ia > ib     wa = recoupling2p(Ordered2p(), leftCsf, rightCsf, rank2, rank1, rank3, ib, ia, ic)
                            wa = (-1)^Int64((Basics.twice(rank1)+Basics.twice(rank2)-Basics.twice(rank3))/2) * wa
        else                error("Recoupling_matrix_3_shells(): program stop A.")    
        end
    elseif ic < ia  && ic < ib
        if      ia < ib     wa = recoupling2p(Ordered2p(), leftCsf, rightCsf, rank3, rank1, rank2, ic, ia, ib)
                            wa = (-1)^Int64(Basics.twice(rank3)) * wa
        elseif  ia > ib     wa = recoupling2p(Ordered2p(), leftCsf, rightCsf, rank3, rank2, rank1, ic, ib, ia)
                            wa = (-1)^Int64((Basics.twice(rank1)+Basics.twice(rank2)+Basics.twice(rank3))/2) * wa
        else                error("Recoupling_matrix_3_shells(): program stop B.")
        end
    else
        if      ia < ib     wa = recoupling2p(Ordered2p(), leftCsf, rightCsf, rank1, rank3, rank2, ia, ic, ib)
                            wa = (-1)^Int64((Basics.twice(rank1)-Basics.twice(rank2)-Basics.twice(rank3))/2) * wa
        elseif  ia > ib     wa = recoupling2p(Ordered2p(), leftCsf, rightCsf, rank2, rank3, rank1, ib, ic, ia)
                            wa = (-1)^Int64(Basics.twice(rank1)) * wa
        else                error("Recoupling_matrix_3_shells(): program stop C.")
        end
    end
    return( wa )
end

"""
`SpinAngularGaigalas.recoupling2p(::Ordered2p, leftCsf::CsfR, rightCsf::CsfR, rank1::AngularJ64, rank2::AngularJ64,
                            rank::AngularJ64, ia::Int64, ib::Int64, ic::Int64)`
    ... computes the recoupling matrix R(ja, jb, Lambda^bra, Lambda^ket) for a scalar operator in the case of three
        interacting shells (shells are ordered). See G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys,
        Vol 30 3747, Eq. (26) (p. 3756); a value::Float64 is returned.

        leftCsf  - is the bra confuguration state function in relativistic notation;
        rightCsf - is the ket confuguration state function in relativistic notation;
        rank_1   - is the first rank of operator;
        rank_2   - is the second rank of operator;
        rank     - is the final rank of operator;
        ia       - is the index of the first interacting subshell;
        ib       - is the index of the second interacting subshell;
        ic       - is the index of the third interacting subshell.
"""
function  recoupling2p(::Ordered2p, leftCsf::CsfR, rightCsf::CsfR, rank1::AngularJ64, rank2::AngularJ64,
                        rank::AngularJ64, ia::Int64, ib::Int64, ic::Int64)
    wa = 1.0 / sqrt((Basics.twice(leftCsf.subshellJ[ia])+1.) * (Basics.twice(leftCsf.subshellJ[ib])+1.) *
                    (Basics.twice(rightCsf.subshellJ[ic])+1.)*(Basics.twice(rank)+1.))
    if ic - ib > 1      wa = wa * recouplingDiagram(DiagramC02(), leftCsf, rightCsf, rank, ib, ic)
                        if abs(wa) < 1.0e-11   return ( wa )   end
    end 
    wa = wa * recouplingDiagram(DiagramC05(), leftCsf, rightCsf, rank, ia, ic)
    if abs(wa) < 1.0e-11   return ( wa )   end
    wa = wa * recouplingDiagram(DiagramC04(), leftCsf, rightCsf, rank1, rank2, rank, ia, ib)
    if abs(wa) < 1.0e-11   return ( wa )   end
    if ia == 1  &&  ib == 2   return( wa)  end
    wa = wa * recouplingDiagram(DiagramC01(), leftCsf, rightCsf, rank1, ia)
    if abs(wa) < 12.0e-11   return ( wa )   end
    if   ia == 1    if ib - 1 - ia <= 1     return( wa )    end
    else;           if ib - ia     <= 1     return( wa )    end
    end
    wa = wa * recouplingDiagram(DiagramC02(), leftCsf, rightCsf, rank1, ia, ib)
    return( wa )
end

"""
`SpinAngularGaigalas.recoupling2p(leftCsf::CsfR, rightCsf::CsfR, rank1::AngularJ64, rank2::AngularJ64, rank3::AngularJ64, 
                            rank4::AngularJ64, rank::AngularJ64, ia::Int64, ib::Int64, ic::Int64, id::Int64)`  
    ... computes the recoupling matrix R(ja, jb, Lambda^bra, Lambda^ket) for a scalar operator in the case of four 
        interacting shells. See G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747, Eq. (33) (p. 3756);
        a value::Float64 is returned.

        leftCsf  - is the bra confuguration state function in relativistic notation;
        rightCsf - is the ket confuguration state function in relativistic notation;
        rank_1   - is the first rank of operator;
        rank_2   - is the second rank of operator;
        rank_3   - is the third rank of operator;
        rank_4   - is the fourth rank of operator;
        rank     - is the final rank of operator;
        ia       - is the index of the first interacting subshell;
        ib       - is the index of the second interacting subshell;
        ic       - is the index of the third interacting subshell;
        id       - is the index of the fourth interacting subshell.
"""
function  recoupling2p(leftCsf::CsfR, rightCsf::CsfR, rank1::AngularJ64, rank2::AngularJ64, rank3::AngularJ64, 
                        rank4::AngularJ64, rank::AngularJ64, ia::Int64, ib::Int64, ic::Int64, id::Int64) 
    wa = 1.0 / sqrt((Basics.twice(leftCsf.subshellJ[ia])+1.) * (Basics.twice(leftCsf.subshellJ[ib])+1.) * 
                    (Basics.twice(leftCsf.subshellJ[ic])+1.)*(Basics.twice(rightCsf.subshellJ[id])+1.) *
                    (Basics.twice(rank4)+1.))
    if ic - ib > 1      wa = wa * recouplingDiagram(DiagramC02(), leftCsf, rightCsf, rank, ib, ic)
                        if  abs(wa) < 1.0e-11   return ( wa )   end
    end
    if id - ic > 1      wa = wa * recouplingDiagram(DiagramC02(), leftCsf, rightCsf, rank4, ic, id)
                        if abs(wa) < 1.0e-11    return ( wa )   end
    end
    wa = wa * recouplingDiagram(DiagramC05(), leftCsf, rightCsf, rank4, ia, id)
    if abs(wa) < 1.0e-11   return ( wa )   end
    wa = wa * recouplingDiagram(DiagramC04(), leftCsf, rightCsf, rank1, rank2, rank, ia, ib)
    if abs(wa) < 1.0e-11   return ( wa )   end
    wa = wa * recouplingDiagram(DiagramC04(), leftCsf, rightCsf, rank, rank3, rank4, ib, ic)
    if abs(wa) < 1.0e-11   return ( wa )   end
    if ia == 1 && ib == 2    return( wa)    end
    wa = wa * recouplingDiagram(DiagramC01(), leftCsf, rightCsf, rank1, ia)
    if abs(wa) < 1.0e-11   return ( wa )   end
    if ia == 1      if ib - 1 - ia <= 1     return( wa)    end
    else;           if ib - ia     <= 1     return( wa)    end
    end
    wa = wa * recouplingDiagram(DiagramC02(), leftCsf, rightCsf, rank1, ia, ib)
    return( wa )
end

"""
`SpinAngularGaigalas.twoParticleDiffOcc2(leftCsf::CsfR, rightCsf::CsfR, creation::Int64, annihilation::Int64, 
                                    subshells::Array{Subshell,1})`
    ... calculates the spin-angular coefficients for a given pair leftCsf, rightCsf of CSF for diagonal matrix elements 
        with respect to configurations. A coeff::Float64 is returned

        leftCsf      - is the bra confuguration state function in relativistic notation;
        rightCsf     - is the ket confuguration state function in relativistic notation;
        creation     - is index of creation operator;
        annihilation - is index of annihilation operator;
        subshells    - is explicitly given relativistic subshell list if useStandardSubshells == false.
"""
function twoParticleDiffOcc2(leftCsf::CsfR, rightCsf::CsfR, creation::Int64, annihilation::Int64, subshells::Array{Subshell,1})
    coeffs2p = Coefficient2p[];    coeffs2pwa = Coefficient2p[]
    for i = 1:length(subshells)
        check = true
        if creation > i         creation_one = i; creation_two = creation
        elseif creation == i    creation_one = creation;    creation_two = creation
        else                    creation_one = creation;    creation_two = i
        end
        if annihilation > i      annihilation_one = i;             annihilation_two = annihilation
            if rightCsf.occupation[i] == 0      check = false    end
        elseif annihilation == i annihilation_one = annihilation;  annihilation_two = annihilation
            if rightCsf.occupation[i] <= 1      check = false    end
        else                     annihilation_one =  annihilation; annihilation_two =  i
            if rightCsf.occupation[i] == 0      check = false    end
        end 
        if check
            coeffs2pwa = twoParticle7to14(leftCsf, rightCsf, creation_one, creation_two, annihilation_one, 
                                                annihilation_two, subshells)
            append!(coeffs2p, coeffs2pwa)
        end
    end
    return( coeffs2p )
end

"""
`SpinAngularGaigalas.twoParticleDiagonal(leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1})`
    ... calculates the spin-angular coefficients for a given pair leftCsf, rightCsf of CSF for diagonal matrix elements 
        with respect to configurations. A coeff::Float64 is returned

        leftCsf   - is the bra confuguration state function in relativistic notation;
        rightCsf  - is the ket confuguration state function in relativistic notation;
        subshells - is explicitly given relativistic subshell list if useStandardSubshells == false.
"""
function twoParticleDiagonal(leftCsf::CsfR, rightCsf::CsfR, subshells::Array{Subshell,1})
    coeffs2p = Coefficient2p[];    coeffs2pwa = Coefficient2p[]
    for  (ib,shib)  in  enumerate(subshells)
        if rightCsf.occupation[ib] != 0 
            if rightCsf.occupation[ib] <= abs(subshells[ib].kappa)*2 
                coeffs2pwa = twoParticle1(leftCsf, rightCsf, ib, subshells)
                append!(coeffs2p, coeffs2pwa)
                if 1 < ib
                    for  (ia,shia)  in  enumerate(subshells)
                        if ia <= ib-1
                            if rightCsf.occupation[ia] != 0 
                                if rightCsf.occupation[ia] <= abs(subshells[ia].kappa)*2 
                                    coeffs2pwa = twoParticle2to5(leftCsf, rightCsf, ia, ib, subshells)
                                    append!(coeffs2p, coeffs2pwa)
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return( coeffs2p )
end

"""
`SpinAngularGaigalas.twoParticle1(leftCsf::CsfR, rightCsf::CsfR, ia::Int64, subshells::Array{Subshell,1})`
    ... calculates the spin-angular coefficients for a given pair leftCsf, rightCsf of CSF for 
        distribution  1  (of Table 1 in G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747): 1.         
        Alpha   Alpha   Alpha   Alpha.  A coeff::Float64 is returned

        leftCsf   - is the bra confuguration state function in relativistic notation;
        rightCsf  - is the ket confuguration state function in relativistic notation;
        ia        - is the index of the interacting subshell;
        subshells - is explicitly given relativistic subshell list if useStandardSubshells == false.
"""
function twoParticle1(leftCsf::CsfR, rightCsf::CsfR, ia::Int64, subshells::Array{Subshell,1})
    coeffs2p = Coefficient2p[]
    wa = 0.0;   waO = 0.0
    if recouplingCheck(leftCsf, rightCsf, ia, ia, ia, ia, length(subshells), 0) != 0.0
        Ni  = leftCsf.occupation[ia];          Nj  = rightCsf.occupation[ia]
        shia = subshells[ia];                  ji  = Basics.subshell_j(shia)
        Ji  = leftCsf.subshellJ[ia];           Jj  = rightCsf.subshellJ[ia]
        si  = leftCsf.seniorityNr[ia];         sj  = rightCsf.seniorityNr[ia]
        Qi  = qshellTermQ(ji,si);              Qj  = qshellTermQ(ji,sj)
        aT  = SpinAngularGaigalas.getTermNumber(ji, Ni, Qi, Ji)
        aTp = SpinAngularGaigalas.getTermNumber(ji, Nj, Qj, Jj)
        waO = irreducibleTensor(SchemeEta_W(), aT, Ni, AngularM64(1//2), AngularM64(-1//2), 0, aTp, Nj)
        if abs(waO) >= 2.0e-10
            waO = waO / sqrt((Basics.twice(ji) + 1.0))
        end
        for kj = 0:Basics.twice(ji)
            wa = (-1)^Int64(Basics.twice(ji) + kj) * waO
            wa = 0.5*((irreducibleTensor(SchemeEta_WW(),aT, Ni, AngularM64(1//2), AngularM64(-1//2), 
                                            AngularM64(1//2), AngularM64(-1//2), kj, aTp, Nj)/sqrt(2.0*kj + 1.0)) - wa) /
                    sqrt(Basics.twice(Ji) + 1.0)
            if  abs(wa) >= 2.0e-10   push!(coeffs2p, Coefficient2p(kj, shia, shia, shia, shia, wa))     end
        end
    end
    return( coeffs2p )
end

"""
`SpinAngularGaigalas.twoParticle2to5(leftCsf::CsfR, rightCsf::CsfR, ia::Int64, ib::Int64, subshells::Array{Subshell,1})`
    ... calculates the spin-angular coefficients for a given pair leftCsf, rightCsf of CSF for distribution 2-5 
        (of Table 1 in G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747). A coeff::Float64 is 
        returned

        leftCsf   - is the bra confuguration state function in relativistic notation;
        rightCsf  - is the ket confuguration state function in relativistic notation;
        ia        - is the index of the first inetracting subshell;
        ib        - is the index of the second inetracting subshell;
        subshells - is explicitly given relativistic subshell list if useStandardSubshells == false.
"""
function twoParticle2to5(leftCsf::CsfR, rightCsf::CsfR, ia::Int64, ib::Int64, subshells::Array{Subshell,1})
    coeffs2p = Coefficient2p[]
    wa = 0.0;       waO = 0.0
    waTemp::Array{Float64,1} = [];     wbTemp::Array{Float64,1} = []
    if ia == ib   return( coeffs2p )   end
    if recouplingCheck(leftCsf, rightCsf, ia, ib, ib, ib, length(subshells), 1) == 0.0   return( coeffs2p )   end
        
    # cases 1212  + + - -     transform to 1122  + - + -
    #       2121                           1122
    aN   = leftCsf.occupation[ia];         aNp  = rightCsf.occupation[ia]
    shia = subshells[ia];                  ja   = Basics.subshell_j(shia)
    Ji   = leftCsf.subshellJ[ia];          Jj   = rightCsf.subshellJ[ia]
    si   = leftCsf.seniorityNr[ia];        sj   = rightCsf.seniorityNr[ia]
    Qi   = qshellTermQ(ja,si);             Qj   = qshellTermQ(ja,sj)
    aT   = SpinAngularGaigalas.getTermNumber(ja, aN, Qi, Ji)
    aTp  = SpinAngularGaigalas.getTermNumber(ja, aNp, Qj, Jj)

    bN   = leftCsf.occupation[ib];         bNp  = rightCsf.occupation[ib]
    shib = subshells[ib];                  jb   = Basics.subshell_j(shib)
    Ji   = leftCsf.subshellJ[ib];          Jj   = rightCsf.subshellJ[ib]
    si   = leftCsf.seniorityNr[ib];        sj   = rightCsf.seniorityNr[ib]
    Qi   = qshellTermQ(jb,si);             Qj   = qshellTermQ(jb,sj)
    bT   = SpinAngularGaigalas.getTermNumber(jb, bN, Qi, Ji)
    bTp  = SpinAngularGaigalas.getTermNumber(jb, bNp, Qj, Jj)

    i_max = min(Basics.twice(ja),Basics.twice(jb))
    for rank = 0:i_max
        kkj  = 2*rank;    kj = AngularJ64(kkj//2)
        push!(waTemp, recoupling2p(leftCsf, rightCsf, kj, ia, ib))
        wb = irreducibleTensor(SchemeEta_W(), aT, aN, AngularM64(1//2), AngularM64(-1//2), rank, aTp, aNp) * 
                irreducibleTensor(SchemeEta_W(), bT, bN, AngularM64(1//2), AngularM64(-1//2), rank, bTp, bNp)
        push!(wbTemp, wb)
        if  abs(wbTemp[rank + 1])  >= 2.0e-10
            wa = waTemp[rank + 1] * wbTemp[rank + 1] / sqrt(2.0*rank + 1.0)
            if abs(wa)  >= 2.0e-10    push!(coeffs2p, Coefficient2p(rank, shia, shib, shia, shib, wa))    end
        end 
    end
    
    # cases 1221  + + - -     transform to 1122  + - + -
    #       2112                           1122
    i_min_2 = Int64(abs(Basics.twice(ja)-Basics.twice(jb))/2)
    i_max_2 = Int64((Basics.twice(ja)+Basics.twice(jb))/2)
    for rank = i_min_2:i_max_2
        wa = 0.0 
        kkj = 2*rank;    kj = AngularJ64(kkj//2)
        for i = 0:i_max
            wb = waTemp[i + 1] * wbTemp[i + 1]
            if abs(wb)  >= 2.0e-10
                wa = wa + wb * sqrt(2.0*i + 1.0) * AngularMomentum.Wigner_6j(ja, jb, kj, jb, ja, i)
            end
        end
        if abs(wa)  >= 2.0e-10    push!(coeffs2p, Coefficient2p(rank, shia, shib, shib, shia, wa))       end
    end
    return( coeffs2p )
end

"""
`SpinAngularGaigalas.twoParticle6(leftCsf::CsfR, rightCsf::CsfR, ia::Int64, ib::Int64, subshells::Array{Subshell,1})`
    ... calculates the spin-angular coefficients for a given pair leftCsf, rightCsf of CSF for 
        distribution 6 (of Table 1 in G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747):  6.         
        Alpha   Alpha   Beta    Beta. A coeff::Float64 is returned.

        leftCsf   - is the bra confuguration state function in relativistic notation;
        rightCsf  - is the ket confuguration state function in relativistic notation;
        ia        - is the index of the first inetracting subshell;
        ib        - is the index of the second inetracting subshell;
        subshells - is explicitly given relativistic subshell list if useStandardSubshells == false.
"""
function twoParticle6(leftCsf::CsfR, rightCsf::CsfR, ia::Int64, ib::Int64, subshells::Array{Subshell,1})
    coeffs2p = Coefficient2p[]
    wa = 0.0;       waO = 0.0
    if  length(subshells) == 1  return( coeffs2p ) end
    no_for_ia = min(ia, ib);               no_for_ib = max(ia, ib)
    if recouplingCheck(leftCsf, rightCsf, no_for_ia, no_for_ib, no_for_ib, no_for_ib, length(subshells), 1) == 0.0
        return( coeffs2p )   
    end
    waTemp::Array{Float64,1} = [];         wbTemp::Array{Float64,1} = []
    aN   = leftCsf.occupation[ia];         aNp  = rightCsf.occupation[ia]
    shia = subshells[ia];                  ja   = Basics.subshell_j(shia)
    Ji   = leftCsf.subshellJ[ia];          Jj   = rightCsf.subshellJ[ia]
    si   = leftCsf.seniorityNr[ia];        sj   = rightCsf.seniorityNr[ia]
    Qi   = qshellTermQ(ja,si);             Qj   = qshellTermQ(ja,sj)
    aT   = SpinAngularGaigalas.getTermNumber(ja, aN, Qi, Ji)
    aTp  = SpinAngularGaigalas.getTermNumber(ja, aNp, Qj, Jj)

    bN   = leftCsf.occupation[ib];         bNp  = rightCsf.occupation[ib]
    shib = subshells[ib];                  jb   = Basics.subshell_j(shib)
    Ji   = leftCsf.subshellJ[ib];          Jj   = rightCsf.subshellJ[ib]
    si   = leftCsf.seniorityNr[ib];        sj   = rightCsf.seniorityNr[ib]
    Qi   = qshellTermQ(jb,si);             Qj   = qshellTermQ(jb,sj)
    bT   = SpinAngularGaigalas.getTermNumber(jb, bN, Qi, Ji)
    bTp  = SpinAngularGaigalas.getTermNumber(jb, bNp, Qj, Jj)

    i_max = min(Basics.twice(ja),Basics.twice(jb))
    for rank = 0:i_max
        kkj = 2*rank;    kj = AngularJ64(kkj//2)
        push!(waTemp, recoupling2p(leftCsf, rightCsf, kj, no_for_ia, no_for_ib))
        wb = irreducibleTensor(SchemeEta_W(), aT, aN, AngularM64(1//2), AngularM64(1//2), rank, aTp, aNp) * 
                irreducibleTensor(SchemeEta_W(), bT, bN, AngularM64(-1//2), AngularM64(-1//2), rank, bTp, bNp)
        push!(wbTemp, wb)
    end
    i_min_2 = Int64(abs(Basics.twice(ja)-Basics.twice(jb))/2)
    i_max_2 = Int64((Basics.twice(ja)+Basics.twice(jb))/2)
    for rank = i_min_2:i_max_2
        wa = 0.0 
        kkj = 2*rank;    kj = AngularJ64(kkj//2)
        for i = 0:i_max
            wb = waTemp[i + 1] * wbTemp[i + 1]
            if abs(wb)  >= 2.0e-10
                wb = (-1)^Int64((Basics.twice(ja)+Basics.twice(jb)+2*rank+2*i)/2) * wb
                wa = wa + sqrt(2.0 * i + 1.0) * AngularMomentum.Wigner_6j(ja, jb, kj, jb, ja, i) * wb
            end
        end
        if abs(wa)  >= 2.0e-10
            wa = -0.5 * wa
            push!(coeffs2p, Coefficient2p(rank, shia, shia, shib, shib, wa))
        end
    end
    return( coeffs2p )
end

"""
`SpinAngularGaigalas.twoParticle7to14(leftCsf::CsfR, rightCsf::CsfR, ia::Int64, ib::Int64, ic::Int64, id::Int64, 
                                    subshells::Array{Subshell,1})`
    ... calculates the spin-angular coefficients for a given pair leftCsf, rightCsf of CSF for distribution 
        7-14 (of Table 1 in G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747). 
        A coeff::Float64 is returned

        leftCsf   - is the bra confuguration state function in relativistic notation;
        rightCsf  - is the ket confuguration state function in relativistic notation;
        ia        - is the index of the first inetracting subshell;
        ib        - is the index of the second inetracting subshell;
        ic        - is the index of the third inetracting subshell;
        id        - is the index of the fourth inetracting subshell;
        subshells - is explicitly given relativistic subshell list if useStandardSubshells == false.
"""
function twoParticle7to14(leftCsf::CsfR, rightCsf::CsfR, ia::Int64, ib::Int64, ic::Int64, id::Int64, 
                                subshells::Array{Subshell,1})
    coeffs2p = Coefficient2p[]
    if  length(subshells) < 2  return( coeffs2p ) end
    if ib == id
        if ia == ib || ic == ib
            if ia == ic    error("SpinAngularGaigalas.twoParticle7to14: stop a")
            elseif ic == ib     coeffs2p = twoParticle7to8(leftCsf, rightCsf, ic, ia, ia, ib, ic, id, subshells)
            else                coeffs2p = twoParticle9to10(leftCsf, rightCsf, ic, ia, ia, ib, ic, id, subshells)
            end
        else                    coeffs2p = twoParticle11to14(leftCsf, rightCsf, ic, ia, ib, ia, ib, ic, id, subshells)
        end
    elseif ia == ic
        if ib  == ia  || id == ia
            if      ib == id    error("SpinAngularGaigalas.twoParticle7to14: stop b")
            elseif  id == ia    coeffs2p = twoParticle7to8(leftCsf, rightCsf, id, ib, ia, ib, ic, id, subshells)
            else                coeffs2p = twoParticle9to10(leftCsf, rightCsf, id, ib, ia, ib, ic, id, subshells)
            end
        else                    coeffs2p = twoParticle11to14(leftCsf, rightCsf, id, ib, ia, ia, ib, ic, id, subshells)
        end
    elseif ia == id
        if ib == ia  || ic == ia
            if      ib == ic    error("SpinAngularGaigalas.twoParticle7to14: stop c")
            elseif  ic == id    coeffs2p = twoParticle7to8(leftCsf, rightCsf, ic, ib, ia, ib, ic, id, subshells)
            else                coeffs2p = twoParticle9to10(leftCsf, rightCsf, ic, ib, ia, ib, ic, id, subshells)
            end
        else                    coeffs2p = twoParticle11to14(leftCsf, rightCsf, ic, ib, ia, ia, ib, ic, id, subshells)
        end
    elseif ib == ic
        if ia == ib || id == ib
            if      ia == id    error("SpinAngularGaigalas.twoParticle7to14: stop d")
            elseif  id == ib    coeffs2p = twoParticle7to8(leftCsf, rightCsf, id, ia, ia, ib, ic, id, subshells)
            else                coeffs2p = twoParticle9to10(leftCsf, rightCsf, id, ia, ia, ib, ic, id, subshells)
            end
        else                    coeffs2p = twoParticle11to14(leftCsf, rightCsf, id, ia, ib, ia, ib, ic, id, subshells)
        end
    else    error("SpinAngularGaigalas.twoParticle7to14: stop e")
    end
    return( coeffs2p )
end

"""
`SpinAngularGaigalas.twoParticle7to8(leftCsf::CsfR, rightCsf::CsfR, ia::Int64, ib::Int64, 
                                iaw::Int64, ibw::Int64, icw::Int64, idw::Int64, subshells::Array{Subshell,1})`
    ... calculates the spin-angular coefficients for a given pair leftCsf, rightCsf of CSF for distribution 7-8 
        (of Table 1 in G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747). A coeff::Float64 
        is returned

        leftCsf   - is the bra confuguration state function in relativistic notation;
        rightCsf  - is the ket confuguration state function in relativistic notation;
        ia        - is the index of the first inetracting subshell;
        ib        - is the index of the second inetracting subshell;
        iaw       - is the first index of the radial integral;
        ibw       - is the second index of the radial integral;
        icw       - is the third index of the  radial integral;
        idw       - is the fourth index of the radial integral;
        subshells - is explicitly given relativistic subshell list if useStandardSubshells == false.
"""
function twoParticle7to8(leftCsf::CsfR, rightCsf::CsfR, ia::Int64, ib::Int64, 
                            iaw::Int64, ibw::Int64, icw::Int64, idw::Int64, subshells::Array{Subshell,1})
    coeffs2p = Coefficient2p[]
    if  length(subshells) == 1  return( coeffs2p ) end
    no_for_ia = min(ia, ib);               no_for_ib = max(ia, ib)
    if recouplingCheck(leftCsf, rightCsf, no_for_ia, no_for_ib, no_for_ib, no_for_ib, length(subshells), 1) == 0.0   
        return( coeffs2p )   
    end
    wbTemp::Array{Float64,1} = [];         wa = 0.0
    aN   = leftCsf.occupation[ia];         aNp  = rightCsf.occupation[ia]
    shia = subshells[ia];                  ja   = Basics.subshell_j(shia)
    Ji   = leftCsf.subshellJ[ia];          Jj   = rightCsf.subshellJ[ia]
    si   = leftCsf.seniorityNr[ia];        sj   = rightCsf.seniorityNr[ia]
    Qi   = qshellTermQ(ja,si);             Qj   = qshellTermQ(ja,sj)
    aT   = SpinAngularGaigalas.getTermNumber(ja, aN, Qi, Ji)
    aTp  = SpinAngularGaigalas.getTermNumber(ja, aNp, Qj, Jj)

    bN   = leftCsf.occupation[ib];         bNp  = rightCsf.occupation[ib]
    shib = subshells[ib];                  jb   = Basics.subshell_j(shib)
    Ji   = leftCsf.subshellJ[ib];          Jj   = rightCsf.subshellJ[ib]
    si   = leftCsf.seniorityNr[ib];        sj   = rightCsf.seniorityNr[ib]
    Qi   = qshellTermQ(jb,si);             Qj   = qshellTermQ(jb,sj)
    bT   = SpinAngularGaigalas.getTermNumber(jb, bN, Qi, Ji)
    bTp  = SpinAngularGaigalas.getTermNumber(jb, bNp, Qj, Jj)

    wc = recoupling2p(leftCsf, rightCsf, jb, no_for_ia, no_for_ib)
    wc = wc * SpinAngularGaigalas.irreducibleTensor(SchemeEta_a(),bT, bN, AngularM64(1//2), bTp, bNp)
    if abs(wc)  <= 2.0e-10   return( coeffs2p )   end
    i_min = Int64(abs(Basics.twice(ja)-Basics.twice(jb))/2)
    i_max = min(Int64(Basics.twice(ja)+Basics.twice(ja)),Int64(Basics.twice(ja)+Basics.twice(jb)))/2
    i_max = Int64(i_max)
    for i = i_min:i_max
        kkj = 2*i;    kj = AngularJ64(kkj//2)
        push!(wbTemp, irreducibleTensor(SchemeEta_aW(), aT, aN, AngularM64(1//2), AngularM64(-1//2), AngularM64(-1//2),
                                i::Int64, jb, aTp, aNp))
    end
    for rank = i_min:i_max
        kkj = 2*rank;    kj = AngularJ64(kkj//2)
        wa = 0.0 
        for i = i_min:i_max
            if (-1)^Int64(Basics.twice(jb)-i + 1) == 1
                ii = i - i_min + 1 
                wb = wbTemp[ii]
                if abs(wb)  >= 2.0e-10
                    wb = sqrt(2.0 * i + 1.0) * AngularMomentum.Wigner_6j(jb, ja, rank, ja, ja, i) * wb
                    wb = (-1)^Int64(Basics.twice(ja)+rank+i) * wb
                    wa = wa + wb
                end
            end
        end
        wa = -wa * wc
        if abs(wa)  >= 2.0e-10
            occup = 0
            if no_for_ia <= no_for_ib-1 
                for i = no_for_ia:no_for_ib-1    occup = occup + leftCsf.occupation[i]    end
                wa = (-1)^(occup+1) * wa
            end
            shiaw = subshells[iaw];    shibw = subshells[ibw]
            shicw = subshells[icw];    shidw = subshells[idw]
            push!(coeffs2p, Coefficient2p(rank, shiaw, shibw, shicw, shidw, wa))
        end
    end
    return( coeffs2p )
end

"""
`SpinAngularGaigalas.twoParticle9to10(leftCsf::CsfR, rightCsf::CsfR, ia::Int64, ib::Int64, 
                                    iaw::Int64, ibw::Int64, icw::Int64, idw::Int64, subshells::Array{Subshell,1})`
    ... calculates the spin-angular coefficients for a given pair leftCsf, rightCsf of CSF for distribution 9-10 
        (of Table 1 in G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747). A coeff::Float64 
        is returned

        leftCsf   - is the bra confuguration state function in relativistic notation;
        rightCsf  - is the ket confuguration state function in relativistic notation;
        ia        - is the index of the first inetracting subshell;
        ib        - is the index of the second inetracting subshell;
        iaw       - is the first index of the radial integral;
        ibw       - is the second index of the radial integral;
        icw       - is the third index of the radial integral;
        idw       - is the fourth index of the radial integral;
        subshells - is explicitly given relativistic subshell list if useStandardSubshells == false.
"""
function twoParticle9to10(leftCsf::CsfR, rightCsf::CsfR, ia::Int64, ib::Int64, 
                                iaw::Int64, ibw::Int64, icw::Int64, idw::Int64, subshells::Array{Subshell,1})
    coeffs2p = Coefficient2p[]
    if  length(subshells) == 1  return( coeffs2p ) end
    no_for_ia = min(ia, ib);               no_for_ib = max(ia, ib)
    if recouplingCheck(leftCsf, rightCsf, no_for_ia, no_for_ib, no_for_ib, no_for_ib, length(subshells), 1) == 0.0   
        return( coeffs2p )   
    end
    wbTemp::Array{Float64,1} = [];         wa = 0.0
    aN   = leftCsf.occupation[ia];         aNp  = rightCsf.occupation[ia]
    shia = subshells[ia];                  ja   = Basics.subshell_j(shia)
    Ji   = leftCsf.subshellJ[ia];          Jj   = rightCsf.subshellJ[ia]
    si   = leftCsf.seniorityNr[ia];        sj   = rightCsf.seniorityNr[ia]
    Qi   = qshellTermQ(ja,si);             Qj   = qshellTermQ(ja,sj)
    aT   = SpinAngularGaigalas.getTermNumber(ja, aN, Qi, Ji)
    aTp  = SpinAngularGaigalas.getTermNumber(ja, aNp, Qj, Jj)

    bN   = leftCsf.occupation[ib];         bNp  = rightCsf.occupation[ib]
    shib = subshells[ib];                  jb   = Basics.subshell_j(shib)
    Ji   = leftCsf.subshellJ[ib];          Jj   = rightCsf.subshellJ[ib]
    si   = leftCsf.seniorityNr[ib];        sj   = rightCsf.seniorityNr[ib]
    Qi   = qshellTermQ(jb,si);             Qj   = qshellTermQ(jb,sj)
    bT   = SpinAngularGaigalas.getTermNumber(jb, bN, Qi, Ji)
    bTp  = SpinAngularGaigalas.getTermNumber(jb, bNp, Qj, Jj)

    wc = recoupling2p(leftCsf, rightCsf, ja, no_for_ia, no_for_ib)
    wc = wc * SpinAngularGaigalas.irreducibleTensor(SchemeEta_a(), aT, aN, AngularM64(-1//2), aTp, aNp)
    if abs(wc)  <= 2.0e-10   return( coeffs2p )   end
    i_min = Int64(abs(Basics.twice(ja)-Basics.twice(jb))/2)
    i_max = min(Int64(Basics.twice(jb)+Basics.twice(jb)),Int64(Basics.twice(ja)+Basics.twice(jb)))/2
    i_max = Int64(i_max)
    for i = i_min:i_max
        kkj = 2*i;    kj = AngularJ64(kkj//2)
        push!(wbTemp, irreducibleTensor(SchemeEta_Wa(),bT, bN, AngularM64(1//2), AngularM64(1//2), AngularM64(-1//2), i, ja, bTp, bNp))
    end
    for rank = i_min:i_max
        kkj = 2*rank;    kj = AngularJ64(kkj//2)
        wa = 0.0 
        for i = i_min:i_max
            if (-1)^Int64(Basics.twice(ja)-i + 1) == 1
                ii = i - i_min + 1 
                wb = wbTemp[ii]
                if abs(wb)  >= 2.0e-10
                    wb = sqrt(2.0 * i + 1.0) * AngularMomentum.Wigner_6j(jb, ja, rank, jb, jb, i) * wb
                    wb = (-1)^Int64(Basics.twice(jb)+rank+i) * wb
                    wa = wa + wb
                end
            end
        end
        wa = -wa * wc
        if abs(wa)  >= 2.0e-10
            occup = 0
            if no_for_ia <= no_for_ib-1 
                for i = no_for_ia:no_for_ib-1    occup = occup + leftCsf.occupation[i]    end
                wa = (-1)^(occup+1) * wa
            end
            shiaw = subshells[iaw];    shibw = subshells[ibw]
            shicw = subshells[icw];    shidw = subshells[idw]
            push!(coeffs2p, Coefficient2p(rank, shiaw, shibw, shicw, shidw, wa))
        end
    end
    return( coeffs2p )
end

"""
`SpinAngularGaigalas.twoParticle11to14(leftCsf::CsfR, rightCsf::CsfR, ia::Int64, ib::Int64, ic::Int64, 
                                    iaw::Int64, ibw::Int64, icw::Int64, idw::Int64, subshells::Array{Subshell,1})`
    ... calculates the spin-angular coefficients for a given pair leftCsf, rightCsf of CSF for distribution 11-14 
        (of Table 1 in G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747). 
        A coeff::Float64 is returned

        leftCsf   - is the bra confuguration state function in relativistic notation;
        rightCsf  - is the ket confuguration state function in relativistic notation;
        ia        - is the index of the first inetracting subshell;
        ib        - is the index of the second inetracting subshell;
        ic        - is the index of the third inetracting subshell;
        iaw       - is the first index of the radial integral;
        ibw       - is the second index of the radial integral;
        icw       - is the third index of the radial integral;
        idw       - is the fourth index of the radial integral;
        subshells - is explicitly given relativistic subshell list if useStandardSubshells == false.
"""
function twoParticle11to14(leftCsf::CsfR, rightCsf::CsfR, ia::Int64, ib::Int64, ic::Int64, 
                                iaw::Int64, ibw::Int64, icw::Int64, idw::Int64, subshells::Array{Subshell,1})
    coeffs2p = Coefficient2p[]
    if  length(subshells) < 3   return( coeffs2p )      end
    iabc = normalForm(ia, ib, ic)
    if recouplingCheck(leftCsf, rightCsf, iabc[1], iabc[3], iabc[2], iabc[2], length(subshells), 2) == 0.0   
        return( coeffs2p )   
    end
    wbTemp::Array{Float64,1} = [];   wcTemp::Array{Float64,1} = [];    wa = 0.0
    aN   = leftCsf.occupation[ia];         aNp  = rightCsf.occupation[ia]
    shia = subshells[ia];                  ja   = Basics.subshell_j(shia)
    Ji   = leftCsf.subshellJ[ia];          Jj   = rightCsf.subshellJ[ia]
    si   = leftCsf.seniorityNr[ia];        sj   = rightCsf.seniorityNr[ia]
    Qi   = qshellTermQ(ja,si);             Qj   = qshellTermQ(ja,sj)
    aT   = SpinAngularGaigalas.getTermNumber(ja, aN, Qi, Ji)
    aTp  = SpinAngularGaigalas.getTermNumber(ja, aNp, Qj, Jj)

    bN   = leftCsf.occupation[ib];         bNp  = rightCsf.occupation[ib]
    shib = subshells[ib];                  jb   = Basics.subshell_j(shib)
    Ji   = leftCsf.subshellJ[ib];          Jj   = rightCsf.subshellJ[ib]
    si   = leftCsf.seniorityNr[ib];        sj   = rightCsf.seniorityNr[ib]
    Qi   = qshellTermQ(jb,si);             Qj   = qshellTermQ(jb,sj)
    bT   = SpinAngularGaigalas.getTermNumber(jb, bN, Qi, Ji)
    bTp  = SpinAngularGaigalas.getTermNumber(jb, bNp, Qj, Jj)

    wd = SpinAngularGaigalas.irreducibleTensor(SchemeEta_a(), aT, aN, AngularM64(-1//2), aTp, aNp) * 
            SpinAngularGaigalas.irreducibleTensor(SchemeEta_a(), bT, bN, AngularM64(1//2), bTp, bNp)
    if abs(wd)  <= 2.0e-10   return( coeffs2p )   end

    cN   = leftCsf.occupation[ic];         cNp  = rightCsf.occupation[ic]
    shic = subshells[ic];                  jc   = Basics.subshell_j(shic)
    Ji   = leftCsf.subshellJ[ic];          Jj   = rightCsf.subshellJ[ic]
    si   = leftCsf.seniorityNr[ic];        sj   = rightCsf.seniorityNr[ic]
    Qi   = qshellTermQ(jc,si);             Qj   = qshellTermQ(jc,sj)
    cT   = SpinAngularGaigalas.getTermNumber(jc, cN, Qi, Ji)
    cTp  = SpinAngularGaigalas.getTermNumber(jc, cNp, Qj, Jj)

    i_min = Int64(abs(Basics.twice(ja)-Basics.twice(jb))/2)
    i_max = min(Int64(Basics.twice(ja)+Basics.twice(jb)),Int64(Basics.twice(jc)+Basics.twice(jc)))/2
    i_max = Int64(i_max)
    for rank = i_min:i_max
        k3 = 2*rank;   rank3 = AngularJ64(k3//2)
        push!(wbTemp, irreducibleTensor(SchemeEta_W(),cT, cN, AngularM64(1//2), AngularM64(-1//2), rank, cTp, cNp))
        push!(wcTemp, recoupling2p(leftCsf, rightCsf, jb, ja, rank3, ib, ia, ic))
    end
    phase = normalPhase(ib, ia, ic, ic)
    occup = 0
    if min(ia,ib) <= max(ia,ib)-1 
        for i = min(ia,ib):max(ia,ib)-1    occup = occup + leftCsf.occupation[i]    end
        phase = (-1)^(occup+1) * phase
    end

    # cases 2313  + + - -     transform to 2133  + - + -
    #       3231                           2133
    for rank = i_min:i_max
        if (-1)^Int64(2*rank) == 1
            ii = rank - i_min + 1 
            wa = phase * wbTemp[ii] *  wcTemp[ii] * wd / sqrt(2.0 * rank + 1.0)
            if abs(wa)  >= 2.0e-10
                shiaw = subshells[iaw];    shibw = subshells[ibw]
                shicw = subshells[icw];    shidw = subshells[idw]
                push!(coeffs2p, Coefficient2p(rank, shiaw, shibw, shicw, shidw, wa))
            end
        end
    end

    # cases 3213  + + - -     transform to 2133  + - + -
    #       2331                           2133
    i_1 = Int64(abs(Basics.twice(ja)-Basics.twice(jc))/2)
    i_2 = Int64(abs(Basics.twice(jb)-Basics.twice(jc))/2)
    i_min1 = max(i_1,i_2)
    i_max1 = min(Int64(Basics.twice(ja)+Basics.twice(jc)),Int64(Basics.twice(jb)+Basics.twice(jc)))/2
    i_max1 = Int64(i_max1)
    if i_min1 <= i_max1
        for rank = i_min1:i_max1
            wa = 0.0
            for i = i_min:i_max
                ii = i - i_min + 1 
                we = phase * wbTemp[ii] *  wcTemp[ii] * wd
                if abs(we)  >= 2.0e-10
                    we = sqrt(2.0 * i + 1.0) * AngularMomentum.Wigner_6j(ja, jc, rank, jc, jb, i) * we
                    wa = wa + we
                end
            end
            if abs(wa)  >= 2.0e-10
                shiaw = subshells[iaw];    shibw = subshells[ibw]
                shicw = subshells[idw];    shidw = subshells[icw]
                push!(coeffs2p, Coefficient2p(rank, shiaw, shibw, shicw, shidw, wa))
            end
        end
    end
    return( coeffs2p )
end

"""
`SpinAngularGaigalas.twoParticle15to18(leftCsf::CsfR, rightCsf::CsfR, ia::Int64, ib::Int64, ic::Int64, id::Int64, 
                                    subshells::Array{Subshell,1})`
    ... calculates the spin-angular coefficients for a given pair leftCsf, rightCsf of CSF for distribution 15-18 
        (of Table 1 in G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747). 
        A coeff::Float64 is returned

        leftCsf   - is the bra confuguration state function in relativistic notation;
        rightCsf  - is the ket confuguration state function in relativistic notation;
        ia        - is the index of the first inetracting subshell;
        ib        - is the index of the second inetracting subshell;
        ic        - is the index of the third inetracting subshell;
        id        - is the index of the fourth inetracting subshell;
        subshells - is explicitly given relativistic subshell list if useStandardSubshells == false.
"""
function twoParticle15to18(leftCsf::CsfR, rightCsf::CsfR, ia::Int64, ib::Int64, ic::Int64, id::Int64, 
                                subshells::Array{Subshell,1})
    coeffs2p = Coefficient2p[]
    if  length(subshells) < 3  return( coeffs2p )   end
    if      ia == ib    coeffs2p = twoParticle15to18order(leftCsf, rightCsf, ic, id, ia, ia, ib, ic, id, 1, subshells)
    elseif  ic == id    coeffs2p = twoParticle15to18order(leftCsf, rightCsf, ia, ib, ic, ia, ib, ic, id, 2, subshells)
    else                error("SpinAngularGaigalas.twoParticle15to18: stop a")
    end
    return( coeffs2p )
end

"""
`SpinAngularGaigalas.twoParticle15to18order(leftCsf::CsfR, rightCsf::CsfR, ia::Int64, ib::Int64, ic::Int64, iaw::Int64, ibw::Int64, 
                                        icw::Int64, idw::Int64, irez::Int64, subshells::Array{Subshell,1})`
    ... calculates the spin-angular coefficients for a given pair leftCsf, rightCsf of CSF for distribution 15-18 
        (of Table 1 in G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747). 
        A coeff::Float64 is returned

        leftCsf   - is the bra confuguration state function in relativistic notation;
        rightCsf  - is the ket confuguration state function in relativistic notation;
        ia        - is the index of the first inetracting subshell;
        ib        - is the index of the second inetracting subshell;
        ic        - is the index of the third inetracting subshell;
        iaw       - is the first index of the radial integral;
        ibw       - is the second index of the radial integral;
        icw       - is the third index of the radial integral;
        idw       - is the fourth index of the radial integral;
        irez      - parameter for selecting calculations;
        subshells - is explicitly given relativistic subshell list if useStandardSubshells == false.
"""
function twoParticle15to18order(leftCsf::CsfR, rightCsf::CsfR, ia::Int64, ib::Int64, ic::Int64, iaw::Int64, ibw::Int64, 
                                    icw::Int64, idw::Int64, irez::Int64, subshells::Array{Subshell,1})
    coeffs2p = Coefficient2p[]
    if  length(subshells) < 3  return( coeffs2p ) end
    iabc = normalForm(ia, ib, ic)
    if recouplingCheck(leftCsf, rightCsf, iabc[1], iabc[3], iabc[2], iabc[2], length(subshells), 2) == 0.0   
        return( coeffs2p )   
    end
    wbTemp::Array{Float64,1} = [];   wcTemp::Array{Float64,1} = [];    wa = 0.0
    if      irez == 1   mLeft1 = AngularM64(-1//2);  mRight1 = AngularM64(-1//2)
                        mLeft2 = AngularM64(1//2);   mRight2 = AngularM64(1//2)
    elseif  irez == 2   mLeft1 = AngularM64(1//2);   mRight1 = AngularM64(1//2)
                        mLeft2 = AngularM64(-1//2);  mRight2 = AngularM64(-1//2)
    else                error("SpinAngularGaigalas.twoParticle15to18order: stop a")    
    end
    aN   = leftCsf.occupation[ia];         aNp  = rightCsf.occupation[ia]
    shia = subshells[ia];                  ja   = Basics.subshell_j(shia)
    Ji   = leftCsf.subshellJ[ia];          Jj   = rightCsf.subshellJ[ia]
    si   = leftCsf.seniorityNr[ia];        sj   = rightCsf.seniorityNr[ia]
    Qi   = qshellTermQ(ja,si);             Qj   = qshellTermQ(ja,sj)
    aT   = SpinAngularGaigalas.getTermNumber(ja, aN, Qi, Ji)
    aTp  = SpinAngularGaigalas.getTermNumber(ja, aNp, Qj, Jj)

    bN   = leftCsf.occupation[ib];         bNp  = rightCsf.occupation[ib]
    shib = subshells[ib];                  jb   = Basics.subshell_j(shib)
    Ji   = leftCsf.subshellJ[ib];          Jj   = rightCsf.subshellJ[ib]
    si   = leftCsf.seniorityNr[ib];        sj   = rightCsf.seniorityNr[ib]
    Qi   = qshellTermQ(jb,si);             Qj   = qshellTermQ(jb,sj)
    bT   = SpinAngularGaigalas.getTermNumber(jb, bN, Qi, Ji)
    bTp  = SpinAngularGaigalas.getTermNumber(jb, bNp, Qj, Jj)

    wd = SpinAngularGaigalas.irreducibleTensor(SchemeEta_a(), aT, aN, mLeft1, aTp, aNp) * 
            SpinAngularGaigalas.irreducibleTensor(SchemeEta_a(), bT, bN, mRight1, bTp, bNp)
    if abs(wd)  <= 2.0e-10   return( coeffs2p )   end

    cN   = leftCsf.occupation[ic];         cNp  = rightCsf.occupation[ic]
    shic = subshells[ic];                  jc   = Basics.subshell_j(shic)
    Ji   = leftCsf.subshellJ[ic];          Jj   = rightCsf.subshellJ[ic]
    si   = leftCsf.seniorityNr[ic];        sj   = rightCsf.seniorityNr[ic]
    Qi   = qshellTermQ(jc,si);             Qj   = qshellTermQ(jc,sj)
    cT   = SpinAngularGaigalas.getTermNumber(jc, cN, Qi, Ji)
    cTp  = SpinAngularGaigalas.getTermNumber(jc, cNp, Qj, Jj)

    i_min = Int64(abs(Basics.twice(ja)-Basics.twice(jb))/2)
    i_max = min(Int64(Basics.twice(ja)+Basics.twice(jb)),Int64(Basics.twice(jc)+Basics.twice(jc)))/2
    i_max = Int64(i_max)
    for rank = i_min:i_max
        k3 = 2*rank;    rank3 = AngularJ64(k3//2)
        push!(wbTemp, irreducibleTensor(SchemeEta_W(),cT, cN, mLeft2, mRight2, rank, cTp, cNp))
        push!(wcTemp, recoupling2p(leftCsf, rightCsf, ja, jb, rank3, ia, ib, ic))
    end
    phase =  normalPhase(ic, ia, ib, ic)
    occup = 0
    if min(ia,ib) <= max(ia,ib)-1 
        for i = min(ia,ib):max(ia,ib)-1    occup = occup + leftCsf.occupation[i]    end
        phase = (-1)^(occup+1) * phase
    end

    # cases 3312  + + - -     transform to 1233  - - + +
    #       3321                           1233

    # cases 1233  + + - -     transform to 1233  + + - -
    #       2133                           1233
    i_1 = Int64(abs(Basics.twice(ja)-Basics.twice(jc))/2)
    i_2 = Int64(abs(Basics.twice(jb)-Basics.twice(jc))/2)
    i_min1 = max(i_1,i_2)
    i_max1 = min(Int64(Basics.twice(ja)+Basics.twice(jc)),Int64(Basics.twice(jb)+Basics.twice(jc)))/2
    i_max1 = Int64(i_max1)
    if i_min1 <= i_max1
        for rank = i_min1:i_max1
            wa = 0.0
            for i = i_min:i_max
                ii = i - i_min + 1 
                if irez == 1
                    phase_b = Int64(Basics.twice(jb) - i + 1)
                elseif irez == 2
                    phase_b = Int64(Basics.twice(ja) - i + 1)
                end
                if (-1)^Int64(phase_b) == 1
                    we = phase * wbTemp[ii] *  wcTemp[ii] * wd
                    if abs(we)  >= 2.0e-10
                        we = -sqrt(2.0 * i + 1.0) * AngularMomentum.Wigner_6j(jc, ja, rank, jb, jc, i) * we
                        if      irez == 1   we = (-1)^Int64((Basics.twice(ja)+Basics.twice(jc)+2*i+2*rank)/2) * we
                        elseif  irez == 2   we = (-1)^Int64((Basics.twice(jb)+Basics.twice(jc)+2*i+2*rank)/2) * we
                        end
                        wa = wa + we
                    end
                end
            end
            if abs(wa)  >= 2.0e-10
                shiaw = subshells[iaw];    shibw = subshells[ibw]
                shicw = subshells[icw];    shidw = subshells[idw]
                push!(coeffs2p, Coefficient2p(rank, shiaw, shibw, shicw, shidw, wa))
            end
        end
    end
    return( coeffs2p )
end

"""
`SpinAngularGaigalas.twoParticle19to42(leftCsf::CsfR, rightCsf::CsfR, ia::Int64, ib::Int64, ic::Int64, id::Int64, 
                                    subshells::Array{Subshell,1})`
    ... calculates the spin-angular coefficients for a given pair leftCsf, rightCsf of CSF for distribution 19-42 
        (of Table 1 in G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747). 
        A coeff::Float64 is returned

        leftCsf   - is the bra confuguration state function in relativistic notation;
        rightCsf  - is the ket confuguration state function in relativistic notation;
        ia        - is the index of the first inetracting subshell;
        ib        - is the index of the second inetracting subshell;
        ic        - is the index of the third inetracting subshell;
        id        - is the index of the fourth inetracting subshell;
        subshells - is explicitly given relativistic subshell list if useStandardSubshells == false.
"""
function twoParticle19to42(leftCsf::CsfR, rightCsf::CsfR, ia::Int64, ib::Int64, ic::Int64, id::Int64, subshells::Array{Subshell,1})
    coeffs2p = Coefficient2p[]
    if  length(subshells) < 4  return( coeffs2p ) end
    if ib < ic                                   coeffs2p = twoParticle19to26(leftCsf, rightCsf, ia, ib, ic, id, 1, subshells)
    elseif ia > id  &&  ib >  id                 coeffs2p = twoParticle19to26(leftCsf, rightCsf, ic, id, ia, ib, 2, subshells)
    elseif ib > ic  &&  ib <  id  &&  ia < ic    coeffs2p = twoParticle27to34(leftCsf, rightCsf, ia, ic, ib, id, 1, subshells)
    elseif ib > ic  &&  ib >  id  &&  ia > ic    coeffs2p = twoParticle27to34(leftCsf, rightCsf, ic, ia, id, ib, 2, subshells)
    elseif ib > ic  &&  ib >  id  &&  ia < ic    coeffs2p = twoParticle35to42(leftCsf, rightCsf, ia, ic, id, ib, 1, subshells)
    elseif ib > ic  &&  ib <  id  &&  ia > ic    coeffs2p = twoParticle35to42(leftCsf, rightCsf, ic, ia, ib, id, 2, subshells)
    else                                         error("SpinAngularGaigalas.twoParticle19to42: stop a")    
    end
    return( coeffs2p )
end

"""
`SpinAngularGaigalas.twoParticle19to26(leftCsf::CsfR, rightCsf::CsfR, ia::Int64, ib::Int64, ic::Int64, id::Int64, 
                                    irez::Int64, subshells::Array{Subshell,1})`
    ... calculates the spin-angular coefficients for a given pair leftCsf, rightCsf of CSF for distribution 
        19-26 (of Table 1 in G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747). 
        A coeff::Float64 is returned

        leftCsf   - is the bra confuguration state function in relativistic notation;
        rightCsf  - is the ket confuguration state function in relativistic notation;
        ia        - is the index of the first inetracting subshell;
        ib        - is the index of the second inetracting subshell;
        ic        - is the index of the third inetracting subshell;
        id        - is the index of the fourth inetracting subshell;
        irez      - parameter for selecting calculations;
        subshells - is explicitly given relativistic subshell list if useStandardSubshells == false.
"""
function twoParticle19to26(leftCsf::CsfR, rightCsf::CsfR, ia::Int64, ib::Int64, ic::Int64, id::Int64, 
                                irez::Int64, subshells::Array{Subshell,1})
    coeffs2p = Coefficient2p[]
    if  length(subshells) < 4  return( coeffs2p ) end
    if recouplingCheck(leftCsf, rightCsf, ia, id, ic, ib, length(subshells), 3) == 0.0   return( coeffs2p )   end
    wcTemp::Array{Float64,1} = [];    wa = 0.0
    if      irez == 1   mLeft1 = AngularM64(1//2);    mRight1 = AngularM64(1//2)
                        mLeft2 = AngularM64(-1//2);   mRight2 = AngularM64(-1//2)
    elseif  irez == 2   mLeft1 = AngularM64(-1//2);   mRight1 = AngularM64(-1//2)
                        mLeft2 = AngularM64(1//2);    mRight2 = AngularM64(1//2)
    else                error("SpinAngularGaigalas.twoParticle19to26: stop a")    
    end
    aN   = leftCsf.occupation[ia];         aNp  = rightCsf.occupation[ia]
    shia = subshells[ia];                  ja   = Basics.subshell_j(shia)
    Ji   = leftCsf.subshellJ[ia];          Jj   = rightCsf.subshellJ[ia]
    si   = leftCsf.seniorityNr[ia];        sj   = rightCsf.seniorityNr[ia]
    Qi   = qshellTermQ(ja,si);             Qj   = qshellTermQ(ja,sj)
    aT   = SpinAngularGaigalas.getTermNumber(ja, aN, Qi, Ji)
    aTp  = SpinAngularGaigalas.getTermNumber(ja, aNp, Qj, Jj)

    bN   = leftCsf.occupation[ib];         bNp  = rightCsf.occupation[ib]
    shib = subshells[ib];                  jb   = Basics.subshell_j(shib)
    Ji   = leftCsf.subshellJ[ib];          Jj   = rightCsf.subshellJ[ib]
    si   = leftCsf.seniorityNr[ib];        sj   = rightCsf.seniorityNr[ib]
    Qi   = qshellTermQ(jb,si);             Qj   = qshellTermQ(jb,sj)
    bT   = SpinAngularGaigalas.getTermNumber(jb, bN, Qi, Ji)
    bTp  = SpinAngularGaigalas.getTermNumber(jb, bNp, Qj, Jj)

    cN   = leftCsf.occupation[ic];         cNp  = rightCsf.occupation[ic]
    shic = subshells[ic];                  jc   = Basics.subshell_j(shic)
    Ji   = leftCsf.subshellJ[ic];          Jj   = rightCsf.subshellJ[ic]
    si   = leftCsf.seniorityNr[ic];        sj   = rightCsf.seniorityNr[ic]
    Qi   = qshellTermQ(jc,si);             Qj   = qshellTermQ(jc,sj)
    cT   = SpinAngularGaigalas.getTermNumber(jc, cN, Qi, Ji)
    cTp  = SpinAngularGaigalas.getTermNumber(jc, cNp, Qj, Jj)

    dN   = leftCsf.occupation[id];         dNp  = rightCsf.occupation[id]
    shid = subshells[id];                  jd   = Basics.subshell_j(shid)
    Ji   = leftCsf.subshellJ[id];          Jj   = rightCsf.subshellJ[id]
    si   = leftCsf.seniorityNr[id];        sj   = rightCsf.seniorityNr[id]
    Qi   = qshellTermQ(jd,si);             Qj   = qshellTermQ(jd,sj)
    dT   = SpinAngularGaigalas.getTermNumber(jd, dN, Qi, Ji)
    dTp  = SpinAngularGaigalas.getTermNumber(jd, dNp, Qj, Jj)

    wd = SpinAngularGaigalas.irreducibleTensor(SchemeEta_a(), aT, aN, mLeft1, aTp, aNp)  * 
            SpinAngularGaigalas.irreducibleTensor(SchemeEta_a(), bT, bN, mRight1, bTp, bNp) * 
            SpinAngularGaigalas.irreducibleTensor(SchemeEta_a(), cT, cN, mLeft2, cTp, cNp)  * 
            SpinAngularGaigalas.irreducibleTensor(SchemeEta_a(), dT, dN, mRight2, dTp, dNp)
    if abs(wd)  <= 2.0e-10   return( coeffs2p )   end
    i_1 = Int64(abs(Basics.twice(ja)-Basics.twice(jb))/2)
    i_2 = Int64(abs(Basics.twice(jc)-Basics.twice(jd))/2)
    i_min = max(i_1,i_2)
    i_max = min(Int64(Basics.twice(ja)+Basics.twice(jb)),Int64(Basics.twice(jc)+Basics.twice(jd)))/2
    i_max = Int64(i_max)
    if i_min > i_max   return( coeffs2p )   end
    for rank = i_min:i_max
        k3 = 2*rank;   rank3 = AngularJ64(k3//2)
        push!(wcTemp, recoupling2p(leftCsf, rightCsf, ja, jb, jc, jd, rank3, ia, ib, ic, id))
    end
    phase = 1
    occup = 0
    for i = ia: ib-1    occup = occup + leftCsf.occupation[i]    end
    phase = (-1)^(occup+1) * phase
    occup = 0
    for i = ic: id-1    occup = occup + leftCsf.occupation[i]    end
    phase = (-1)^(occup+1) * phase

    # cases 1234  + + - -     transform to 1234  + + - -
    #       2134                           1234
    #                                                   (irez = 1)

    # cases 3412  + + - -     transform to 1234  - - + +
    #       3421                           1234
    #                                                   (irez = 2)
    i_1 = Int64(abs(Basics.twice(ja)-Basics.twice(jc))/2)
    i_2 = Int64(abs(Basics.twice(jb)-Basics.twice(jd))/2)
    i_min1 = max(i_1,i_2)
    i_max1 = min(Int64(Basics.twice(ja)+Basics.twice(jc)),Int64(Basics.twice(jb)+Basics.twice(jd)))/2
    i_max1 = Int64(i_max1)
    if i_min1 <= i_max1
        for rank = i_min1:i_max1
            wa = 0.0
            for i = i_min:i_max
                ii = i - i_min + 1 
                if      irez == 1   phase_b = Int64(Basics.twice(ja) + Basics.twice(jd) - 2*i)
                elseif  irez == 2   phase_b = Int64(Basics.twice(jb) + Basics.twice(jc) - 2*i)
                end
                if (-1)^Int64(phase_b) == 1
                    we = phase * wcTemp[ii] * wd
                    if abs(we)  >= 2.0e-10
                        we = -sqrt(2.0 * i + 1.0) * AngularMomentum.Wigner_6j(ja, jc, rank, jd, jb, i) * we
                        if      irez == 1   we = (-1)^Int64((Basics.twice(jb)+Basics.twice(jc)+2*i+2*rank)/2) * we
                        elseif  irez == 2   we = (-1)^Int64((Basics.twice(ja)+Basics.twice(jd)-2*i-2*rank)/2) * we
                        end
                        wa = wa + we
                    end
                end
            end
            if abs(wa)  >= 2.0e-10
                if irez == 1
                    shiaw = subshells[ia];    shibw = subshells[ib]
                    shicw = subshells[ic];    shidw = subshells[id]
                    push!(coeffs2p, Coefficient2p(rank, shiaw, shibw, shicw, shidw, wa))
                elseif irez == 2
                    shiaw = subshells[ic];    shibw = subshells[id]
                    shicw = subshells[ia];    shidw = subshells[ib]
                    push!(coeffs2p, Coefficient2p(rank, shiaw, shibw, shicw, shidw, wa))
                end
            end
        end 
    end

    # cases 1243  + + - -     transform to 1234  + + - -
    #       2134                           1234
    #                                                   (irez = 1)

    # cases 3421  + + - -     transform to 1234  - - + +
    #       4321                           1234
    #                                                   (irez = 2)
    i_1    = Int64(abs(Basics.twice(ja)-Basics.twice(jd))/2)
    i_2    = Int64(abs(Basics.twice(jb)-Basics.twice(jc))/2)
    i_min1 = max(i_1,i_2)
    i_max1 = min(Int64(Basics.twice(ja)+Basics.twice(jd)),Int64(Basics.twice(jb)+Basics.twice(jc)))/2
    i_max1 = Int64(i_max1)
    if i_min1 <= i_max1
        for rank = i_min1:i_max1
            wa = 0.0
            for i = i_min:i_max
                ii = i - i_min + 1 
                if      irez == 1   phase_b = Int64(Basics.twice(ja) + Basics.twice(jd))
                elseif  irez == 2   phase_b = Int64(Basics.twice(jc) + Basics.twice(jb))
                end
                if (-1)^Int64(phase_b) == 1
                    we = phase * wcTemp[ii] * wd
                    if abs(we)  >= 2.0e-10
                        we = sqrt(2.0 * i + 1.0) * AngularMomentum.Wigner_6j(ja, jd, rank, jc, jb, i) * we
                        if      irez == 1   we = (-1)^Int64((Basics.twice(jb)+Basics.twice(jc)+2*Basics.twice(jd)-2*rank)/2)*we
                        elseif  irez == 2   we = (-1)^Int64((Basics.twice(ja)+2*Basics.twice(jb)+Basics.twice(jd)+4*i+2*rank)/2)*we
                        end
                        wa = wa + we
                    end
                end
            end
            if abs(wa)  >= 2.0e-10
                if      irez == 1
                    shiaw = subshells[ia];    shibw = subshells[ib]
                    shicw = subshells[id];    shidw = subshells[ic]
                    push!(coeffs2p, Coefficient2p(rank, shiaw, shibw, shicw, shidw, wa))
                elseif  irez == 2
                    shiaw = subshells[ic];    shibw = subshells[id]
                    shicw = subshells[ib];    shidw = subshells[ia]
                    push!(coeffs2p, Coefficient2p(rank, shiaw, shibw, shicw, shidw, wa))
                end
            end
        end 
    end
    return( coeffs2p )
end

"""
`SpinAngularGaigalas.twoParticle27to34(leftCsf::CsfR, rightCsf::CsfR, ia::Int64, ib::Int64, ic::Int64, id::Int64, 
                                    irez::Int64, subshells::Array{Subshell,1})`
    ... calculates the spin-angular coefficients for a given pair leftCsf, rightCsf of CSF for distribution 27-34 
        (of Table 1 in G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747). 
        A coeff::Float64 is returned

        leftCsf   - is the bra confuguration state function in relativistic notation;
        rightCsf  - is the ket confuguration state function in relativistic notation;
        ia        - is the index of the first inetracting subshell;
        ib        - is the index of the second inetracting subshell;
        ic        - is the index of the third inetracting subshell;
        id        - is the index of the fourth inetracting subshell;
        irez      - parameter for selecting calculations;
        subshells - is explicitly given relativistic subshell list if useStandardSubshells == false.
"""
function twoParticle27to34(leftCsf::CsfR, rightCsf::CsfR, ia::Int64, ib::Int64, ic::Int64, id::Int64, 
                                irez::Int64, subshells::Array{Subshell,1})
    coeffs2p = Coefficient2p[]
    if  length(subshells) < 4  return( coeffs2p ) end
    if recouplingCheck(leftCsf, rightCsf, ia, id, ic, ib, length(subshells), 3) == 0.0   return( coeffs2p )   end
    wcTemp::Array{Float64,1} = [];    wa = 0.0
    if      irez == 1   mLeft1 = AngularM64(1//2);    mRight1 = AngularM64(-1//2)
                        mLeft2 = AngularM64(1//2);    mRight2 = AngularM64(-1//2)
    elseif  irez == 2   mLeft1 = AngularM64(-1//2);   mRight1 = AngularM64(1//2)
                        mLeft2 = AngularM64(-1//2);   mRight2 = AngularM64(1//2)
    else                error("SpinAngularGaigalas.twoParticle27to34: stop a")
    end
    
    aN   = leftCsf.occupation[ia];         aNp  = rightCsf.occupation[ia]
    shia = subshells[ia];                  ja   = Basics.subshell_j(shia)
    Ji   = leftCsf.subshellJ[ia];          Jj   = rightCsf.subshellJ[ia]
    si   = leftCsf.seniorityNr[ia];        sj   = rightCsf.seniorityNr[ia]
    Qi   = qshellTermQ(ja,si);             Qj   = qshellTermQ(ja,sj)
    aT   = SpinAngularGaigalas.getTermNumber(ja, aN, Qi, Ji)
    aTp  = SpinAngularGaigalas.getTermNumber(ja, aNp, Qj, Jj)

    bN   = leftCsf.occupation[ib];         bNp  = rightCsf.occupation[ib]
    shib = subshells[ib];                  jb   = Basics.subshell_j(shib)
    Ji   = leftCsf.subshellJ[ib];          Jj   = rightCsf.subshellJ[ib]
    si   = leftCsf.seniorityNr[ib];        sj   = rightCsf.seniorityNr[ib]
    Qi   = qshellTermQ(jb,si);             Qj   = qshellTermQ(jb,sj)
    bT   = SpinAngularGaigalas.getTermNumber(jb, bN, Qi, Ji)
    bTp  = SpinAngularGaigalas.getTermNumber(jb, bNp, Qj, Jj)

    cN   = leftCsf.occupation[ic];         cNp  = rightCsf.occupation[ic]
    shic = subshells[ic];                  jc   = Basics.subshell_j(shic)
    Ji   = leftCsf.subshellJ[ic];          Jj   = rightCsf.subshellJ[ic]
    si   = leftCsf.seniorityNr[ic];        sj   = rightCsf.seniorityNr[ic]
    Qi   = qshellTermQ(jc,si);             Qj   = qshellTermQ(jc,sj)
    cT   = SpinAngularGaigalas.getTermNumber(jc, cN, Qi, Ji)
    cTp  = SpinAngularGaigalas.getTermNumber(jc, cNp, Qj, Jj)

    dN   = leftCsf.occupation[id];         dNp  = rightCsf.occupation[id]
    shid = subshells[id];                  jd   = Basics.subshell_j(shid)
    Ji   = leftCsf.subshellJ[id];          Jj   = rightCsf.subshellJ[id]
    si   = leftCsf.seniorityNr[id];        sj   = rightCsf.seniorityNr[id]
    Qi   = qshellTermQ(jd,si);             Qj   = qshellTermQ(jd,sj)
    dT   = SpinAngularGaigalas.getTermNumber(jd, dN, Qi, Ji)
    dTp  = SpinAngularGaigalas.getTermNumber(jd, dNp, Qj, Jj)

    wd = SpinAngularGaigalas.irreducibleTensor(SchemeEta_a(), aT, aN, mLeft1, aTp, aNp) * 
            SpinAngularGaigalas.irreducibleTensor(SchemeEta_a(), bT, bN, mRight1, bTp, bNp) * 
            SpinAngularGaigalas.irreducibleTensor(SchemeEta_a(), cT, cN, mLeft2, cTp, cNp) * 
            SpinAngularGaigalas.irreducibleTensor(SchemeEta_a(), dT, dN, mRight2, dTp, dNp)
    if abs(wd)  <= 2.0e-10   return( coeffs2p )   end
    i_1   = Int64(abs(Basics.twice(ja)-Basics.twice(jb))/2)
    i_2   = Int64(abs(Basics.twice(jc)-Basics.twice(jd))/2)
    i_min = max(i_1,i_2)
    i_max = min(Int64(Basics.twice(ja)+Basics.twice(jb)),Int64(Basics.twice(jc)+Basics.twice(jd)))/2
    i_max = Int64(i_max)
    if i_min > i_max   return( coeffs2p )   end
    for rank = i_min:i_max
        k3 = 2*rank;   rank3 = AngularJ64(k3//2)
        push!(wcTemp, recoupling2p(leftCsf, rightCsf, ja, jb, jc, jd, rank3, ia, ib, ic, id))
    end
    phase = 1
    occup = 0
    for i = ia: ib-1    occup = occup + leftCsf.occupation[i]    end
    phase = (-1)^(occup+1) * phase
    occup = 0
    for i = ic: id-1    occup = occup + leftCsf.occupation[i]    end
    phase = (-1)^(occup+1) * phase

    # cases 1324  + + - -     transform to 1234  - + - +
    #       1342                           1234
    #                                                   (irez = 1)

    # cases 2413  + + - -     transform to 1234  + - + -
    #       4231                           1234
    #                                                   (irez = 2)
    for rank = i_min:i_max
        ii = rank - i_min + 1
        wa = phase * wcTemp[ii] * wd /sqrt(2.0 * rank + 1.0)
        if abs(wa)  >= 2.0e-10
            if      irez == 1
                shiaw = subshells[ia];    shibw = subshells[ic]
                shicw = subshells[ib];    shidw = subshells[id]
                push!(coeffs2p, Coefficient2p(rank, shiaw, shibw, shicw, shidw, wa))
            elseif  irez == 2
                wa = (-1)^Int64((Basics.twice(ja)+Basics.twice(jb)+Basics.twice(jc)+Basics.twice(jd)+4*rank)/2) * wa
                shiaw = subshells[ib];    shibw = subshells[id]
                shicw = subshells[ia];    shidw = subshells[ic]
                push!(coeffs2p, Coefficient2p(rank, shiaw, shibw, shicw, shidw, wa))
            end
        end
    end

    # cases 1342  + + - -     transform to 1234  + - + -
    #       3124                           1234
    #                                                   (irez = 1)

    # cases 2431  + + - -     transform to 1234  - + - +
    #       4213                           1234
    #                                                   (irez = 2)
    i_1    = Int64(abs(Basics.twice(ja)-Basics.twice(jd))/2)
    i_2    = Int64(abs(Basics.twice(jb)-Basics.twice(jc))/2)
    i_min1 = max(i_1,i_2)
    i_max1 = min(Int64(Basics.twice(ja)+Basics.twice(jd)),Int64(Basics.twice(jb)+Basics.twice(jc)))/2
    i_max1 = Int64(i_max1)
    if i_min1 <= i_max1
        for rank = i_min1:i_max1
            wa = 0.0
            for i = i_min:i_max
                ii = i - i_min + 1 
                if      irez == 1   phase_b = Int64(Basics.twice(ja) - Basics.twice(jc)) + 2*i
                elseif  irez == 2   phase_b = Int64(Basics.twice(jd) - Basics.twice(jb)) - 2*i
                end
                if (-1)^Int64(phase_b) == 1
                    we = phase * wcTemp[ii] * wd
                    if abs(we)  >= 2.0e-10
                        we =  -sqrt(2.0 * i + 1.0) * AngularMomentum.Wigner_6j(ja, jd, rank, jc, jb, i) * we
                        if      irez == 1   we = (-1)^Int64((2*Basics.twice(jc)-4*i-4*rank)/2) * we
                        elseif  irez == 2   we = (-1)^Int64((Basics.twice(ja)+Basics.twice(jb)+Basics.twice(jc) -
                                                                Basics.twice(jd)-4*rank)/2) * we
                        end
                        wa = wa + we
                    end
                end
            end
            if abs(wa)  >= 2.0e-10
                if      irez == 1   
                    shiaw = subshells[ia];    shibw = subshells[ic]
                    shicw = subshells[id];    shidw = subshells[ib]
                    push!(coeffs2p, Coefficient2p(rank, shiaw, shibw, shicw, shidw, wa))
                elseif  irez == 2
                    shiaw = subshells[ib];    shibw = subshells[id]
                    shicw = subshells[ic];    shidw = subshells[ia]
                    push!(coeffs2p, Coefficient2p(rank, shiaw, shibw, shicw, shidw, wa))
                end
            end
        end
    end
    return( coeffs2p )
end

"""
`SpinAngularGaigalas.twoParticle35to42(leftCsf::CsfR, rightCsf::CsfR, ia::Int64, ib::Int64, ic::Int64, id::Int64, 
                                    irez::Int64, subshells::Array{Subshell,1})`
    ... calculates the spin-angular coefficients for a given pair leftCsf, rightCsf of CSF for distribution 35-42 
        (of Table 1 in G. Gaigalas et al., 1997 J. Phys. B: At. Mol. Opt. Phys, Vol 30 3747). 
        A coeff::Float64 is returned

        leftCsf   - is the bra confuguration state function in relativistic notation;
        rightCsf  - is the ket confuguration state function in relativistic notation;
        ia        - is the index of the first inetracting subshell;
        ib        - is the index of the second inetracting subshell;
        ic        - is the index of the third inetracting subshell;
        id        - is the index of the fourth inetracting subshell;
        irez      - parameter for selecting calculations;
        subshells - is explicitly given relativistic subshell list if useStandardSubshells == false.
"""
function twoParticle35to42(leftCsf::CsfR, rightCsf::CsfR, ia::Int64, ib::Int64, ic::Int64, id::Int64, 
                                irez::Int64, subshells::Array{Subshell,1})
    coeffs2p = Coefficient2p[]
    if  length(subshells) < 4  return( coeffs2p ) end
    if recouplingCheck(leftCsf, rightCsf, ia, id, ic, ib, length(subshells), 3) == 0.0   return( coeffs2p )   end
    wcTemp::Array{Float64,1} = [];    wa = 0.0
    if      irez == 1   mLeft1 = AngularM64(1//2);    mRight1 = AngularM64(-1//2)
                        mLeft2 = AngularM64(-1//2);   mRight2 = AngularM64(1//2)
    elseif  irez == 2   mLeft1 = AngularM64(-1//2);   mRight1 = AngularM64(1//2)
                        mLeft2 = AngularM64(1//2);    mRight2 = AngularM64(-1//2)
    else                error("SpinAngularGaigalas.twoParticle35to42: stop a")
    end
    
    aN   = leftCsf.occupation[ia];         aNp  = rightCsf.occupation[ia]
    shia = subshells[ia];                  ja   = Basics.subshell_j(shia)
    Ji   = leftCsf.subshellJ[ia];          Jj   = rightCsf.subshellJ[ia]
    si   = leftCsf.seniorityNr[ia];        sj   = rightCsf.seniorityNr[ia]
    Qi   = qshellTermQ(ja,si);             Qj   = qshellTermQ(ja,sj)
    aT   = SpinAngularGaigalas.getTermNumber(ja, aN, Qi, Ji)
    aTp  = SpinAngularGaigalas.getTermNumber(ja, aNp, Qj, Jj)

    bN   = leftCsf.occupation[ib];         bNp  = rightCsf.occupation[ib]
    shib = subshells[ib];                  jb   = Basics.subshell_j(shib)
    Ji   = leftCsf.subshellJ[ib];          Jj   = rightCsf.subshellJ[ib]
    si   = leftCsf.seniorityNr[ib];        sj   = rightCsf.seniorityNr[ib]
    Qi   = qshellTermQ(jb,si);             Qj   = qshellTermQ(jb,sj)
    bT   = SpinAngularGaigalas.getTermNumber(jb, bN, Qi, Ji)
    bTp  = SpinAngularGaigalas.getTermNumber(jb, bNp, Qj, Jj)

    cN   = leftCsf.occupation[ic];         cNp  = rightCsf.occupation[ic]
    shic = subshells[ic];                  jc   = Basics.subshell_j(shic)
    Ji   = leftCsf.subshellJ[ic];          Jj   = rightCsf.subshellJ[ic]
    si   = leftCsf.seniorityNr[ic];        sj   = rightCsf.seniorityNr[ic]
    Qi   = qshellTermQ(jc,si);             Qj   = qshellTermQ(jc,sj)
    cT   = SpinAngularGaigalas.getTermNumber(jc, cN, Qi, Ji)
    cTp  = SpinAngularGaigalas.getTermNumber(jc, cNp, Qj, Jj)

    dN   = leftCsf.occupation[id];         dNp  = rightCsf.occupation[id]
    shid = subshells[id];                  jd   = Basics.subshell_j(shid)
    Ji   = leftCsf.subshellJ[id];          Jj   = rightCsf.subshellJ[id]
    si   = leftCsf.seniorityNr[id];        sj   = rightCsf.seniorityNr[id]
    Qi   = qshellTermQ(jd,si);             Qj   = qshellTermQ(jd,sj)
    dT   = SpinAngularGaigalas.getTermNumber(jd, dN, Qi, Ji)
    dTp  = SpinAngularGaigalas.getTermNumber(jd, dNp, Qj, Jj)

    wd = SpinAngularGaigalas.irreducibleTensor(SchemeEta_a(), aT, aN, mLeft1, aTp, aNp)  * 
            SpinAngularGaigalas.irreducibleTensor(SchemeEta_a(), bT, bN, mRight1, bTp, bNp) * 
            SpinAngularGaigalas.irreducibleTensor(SchemeEta_a(), cT, cN, mLeft2, cTp, cNp)  * 
            SpinAngularGaigalas.irreducibleTensor(SchemeEta_a(), dT, dN, mRight2, dTp, dNp)
    if abs(wd)  <= 2.0e-10   return( coeffs2p )   end
    i_1   = Int64(abs(Basics.twice(ja)-Basics.twice(jb))/2)
    i_2   = Int64(abs(Basics.twice(jc)-Basics.twice(jd))/2)
    i_min = max(i_1,i_2)
    i_max = min(Int64(Basics.twice(ja)+Basics.twice(jb)),Int64(Basics.twice(jc)+Basics.twice(jd)))/2
    i_max = Int64(i_max)
    if i_min > i_max   return( coeffs2p )   end
    for rank = i_min:i_max
        k3 = 2*rank;   rank3 = AngularJ64(k3//2)
        push!(wcTemp, recoupling2p(leftCsf, rightCsf, ja, jb, jc, jd, rank3, ia, ib, ic, id))
    end
    phase = 1
    occup = 0
    for i = ia: ib-1    occup = occup + leftCsf.occupation[i]    end
    phase = (-1)^(occup+1) * phase
    occup = 0
    for i = ic: id-1    occup = occup + leftCsf.occupation[i]    end
    phase = (-1)^(occup+1) * phase

    # cases 1423  + + - -     transform to 1234  + - - +
    #       4132                           1234
    #                                                   (irez = 1)

    # cases 2314  + + - -     transform to 1234  - + + -
    #       3241                           1234
    #                                                   (irez = 2)
    for rank = i_min:i_max
        ii = rank - i_min + 1
        wa = phase * wcTemp[ii] * wd /sqrt(2.0 * rank + 1.0)
        if abs(wa)  >= 2.0e-10
            if      irez == 1
                wa = (-1)^Int64((Basics.twice(jc)+Basics.twice(jd)-2*rank+2)/2) * wa
                shiaw = subshells[ia];    shibw = subshells[id]
                shicw = subshells[ib];    shidw = subshells[ic]
                push!(coeffs2p, Coefficient2p(rank, shiaw, shibw, shicw, shidw, wa))
            elseif  irez == 2
                wa = (-1)^Int64((Basics.twice(ja)+Basics.twice(jb)-2*rank+2)/2) * wa
                shiaw = subshells[ib];    shibw = subshells[ic]
                shicw = subshells[ia];    shidw = subshells[id]
                push!(coeffs2p, Coefficient2p(rank, shiaw, shibw, shicw, shidw, wa))
            end
        end
    end

    # cases 1432  + + - -     transform to 1234  + - - +
    #       4132                           1234
    #                                                   (irez = 1)

    # cases 2341  + + - -     transform to 1234  - + + -
    #       3214                           1234
    #                                                   (irez = 2)
    i_1    = Int64(abs(Basics.twice(ja)-Basics.twice(jc))/2)
    i_2    = Int64(abs(Basics.twice(jb)-Basics.twice(jd))/2)
    i_min1 = max(i_1,i_2)
    i_max1 = min(Int64(Basics.twice(ja)+Basics.twice(jc)),Int64(Basics.twice(jb)+Basics.twice(jd)))/2
    i_max1 = Int64(i_max1)
    if i_min1 <= i_max1
        for rank = i_min1:i_max1
            wa = 0.0
            for i = i_min:i_max
                ii = i - i_min + 1 
                if      irez == 1   phase_b = Int64(Basics.twice(ja) + Basics.twice(jd)) + 2*i
                elseif  irez == 2   phase_b = Int64(Basics.twice(jb) - Basics.twice(jc)) + 2*i
                end
                if (-1)^Int64(phase_b) == 1
                    we = phase * wcTemp[ii] * wd
                    if abs(we)  >= 2.0e-10
                        we =  sqrt(2.0 * i + 1.0) * AngularMomentum.Wigner_6j(ja, jc, rank, jd, jb, i) * we
                        if      irez == 1   we = (-1)^Int64((Basics.twice(jc)-Basics.twice(jd)-4*rank+2*i)/2) * we
                        elseif  irez == 2   we = (-1)^Int64((Basics.twice(ja)-Basics.twice(jb)-4*rank+2*i)/2) * we
                        end
                        wa = wa + we
                    end
                end
            end
            if abs(wa)  >= 2.0e-10
                if      irez == 1
                    shiaw = subshells[ia];    shibw = subshells[id]
                    shicw = subshells[ic];    shidw = subshells[ib]
                    push!(coeffs2p, Coefficient2p(rank, shiaw, shibw, shicw, shidw, wa))
                elseif  irez == 2
                    shiaw = subshells[ib];    shibw = subshells[ic]
                    shicw = subshells[id];    shidw = subshells[ia]
                    push!(coeffs2p, Coefficient2p(rank, shiaw, shibw, shicw, shidw, wa))
                end
            end
        end
    end
    return( coeffs2p )
end

end # module
