
"""
`module JAC.AngularMomentum`  
	... a submodel of JAC that contains various methods to calculate (standard) symbols, coefficients and functions from 
	    the theory of angular momentum and Racah's algebra.
"""
module AngularMomentum

using  SpecialFunctions, WignerSymbols, ..Basics
using  GSL: sf_coupling_3j, sf_coupling_6j, sf_coupling_9j, sf_legendre_sphPlm

"""
`AngularMomentum.allowedDoubleKappaCouplingSequence(syma::LevelSymmetry, symb::LevelSymmetry, maxKappa::Int64)`  
    ... to determine all allowed coupling sequences that fulfill
        syma + (kappa1, symx, kappa2) --> symb  ==  syma + kappa1 --> symx + kappa2 --> symb,
        and where kappa1|, |kappa2| <= maxKappa. A list::Array{Tuple{Int64,LevelSymmetry,Int64},1}) of (kappa1, symx, kappa2)
        is returned.
"""
function  allowedDoubleKappaCouplingSequence(syma::LevelSymmetry, symb::LevelSymmetry, maxKappa::Int64)
    couplings = Tuple{Int64,LevelSymmetry,Int64}[];     kx = abs(maxKappa)
    for  kappa1 = -kx:kx
        if  kappa1 == 0          continue    end
        for  kappa2 = -kx:kx
            if  kappa2 == 0      continue    end
            symx = allowedDoubleKappaSymmetries(syma, kappa1, kappa2, symb)
            for  sx in  symx    push!( couplings, (kappa1, sx, kappa2))      end
        end
    end

    return( couplings )         
end

"""
`AngularMomentum.allowedDoubleKappas(syma::LevelSymmetry, symb::LevelSymmetry, maxKappa::Int64)`  
    ... to determine all allowed pairs of kappa that fulfill  syma + (kappa1, kappa2) --> symb,
        and where |kappa1|, |kappa2| <= maxKappa. A list::Array{Tuple{Int64,Int64},1}) of (kappa1, kappa2) is returned.
"""
function  allowedDoubleKappas(syma::LevelSymmetry, symb::LevelSymmetry, maxKappa::Int64)
    kappaPairs = Tuple{Int64,Int64}[];     kx = abs(maxKappa)
    for  kappa1 = -kx:kx
        if  kappa1 == 0          continue    end
        for  kappa2 = -kx:kx
            if  kappa2 == 0      continue    end
            symx = allowedDoubleKappaSymmetries(syma, kappa1, kappa2, symb)
            if  length(symx) > 0    push!( kappaPairs, (kappa1, kappa2))      end
        end
    end

    return( kappaPairs )         
end

"""
`AngularMomentum.allowedDoubleKappaSymmetries(syma::LevelSymmetry, kappa1::Int64, kappa2::Int64, symb::LevelSymmetry)`  
    ... to determine all allowed level symmetries symx that can be coupled to the sequence 
        syma + kappa1 --> {symx} + kappa2 --> symb. A list::Array{LevelSymmetry,1} of symx is returned.
"""
function  allowedDoubleKappaSymmetries(syma::LevelSymmetry, kappa1::Int64, kappa2::Int64, symb::LevelSymmetry)
    symx1 = AngularMomentum.allowedTotalSymmetries(syma, kappa1)
    symx2 = AngularMomentum.allowedTotalSymmetries(symb, kappa2)
    symx  = intersect(symx1, symx2)
    return( symx )         
end

"""
`AngularMomentum.allowedKappaSymmetries(syma::LevelSymmetry, symb::LevelSymmetry)`  
    ... to determine all allowed single-electron symmetries/partial waves kappa (l,j) that can be coupled to the given 
        level symmetries. A list::Array{Int64,1} of kappa-values is returned.
"""
function  allowedKappaSymmetries(syma::LevelSymmetry, symb::LevelSymmetry)
    kappaList = Int64[];    JList = Basics.oplus(syma.J, symb.J);    kappaMax = syma.J.num + symb.J.num + 2
    for  kappa = -kappaMax:kappaMax
            kappa == 0  &&  continue
            if  syma.parity == Basics.plus   la = 0               else   la = 1       end
            if  symb.parity == Basics.plus   lb = 0               else   lb = 1       end
            if  kappa < 0                 l  = abs(kappa) -1   else   l  = kappa   end
            isodd( la + lb + l )          && continue
            for  j in JList
                if  j == AngularJ64(abs(kappa) - 1//2)   push!( kappaList, kappa)  end
            end
    end
    return( kappaList )         
end

"""
`AngularMomentum.allowedMultipoleSymmetries(syma::LevelSymmetry, multipole::EmMultipole)`  
    ... to determine all allowed level symmetries for which the given multipole can give rise to a non-zero (transition) 
        amplitude; a symList::Array{LevelSymmetry,1} is returned.
"""
function  allowedMultipoleSymmetries(syma::LevelSymmetry, multipole::EmMultipole)
    symList = LevelSymmetry[];    JList = Basics.oplus(syma.J, AngularJ64(multipole.L) )
    for  J in JList
        if      parityEmMultipolePi(syma.parity, multipole, Basics.plus)     push!( symList, LevelSymmetry(J, Basics.plus) )
        elseif  parityEmMultipolePi(syma.parity, multipole, Basics.minus)    push!( symList, LevelSymmetry(J, Basics.minus) )
        else    error("stop a")
        end
    end
    return( symList )           
end

"""
`AngularMomentum.allowedTotalSymmetries(syma::LevelSymmetry, kappa::Int64)`  
    ... to determine all allowed total symmetries J^P that can be constructed by coupling a partial wave kappa (l,j) 
        to the given level symmetry syma. A list::Array{LevelSymmetry,1} of total symmetries is returned.
"""
function allowedTotalSymmetries(syma::LevelSymmetry, kappa::Int64) 
    symtList = LevelSymmetry[]
    if  kappa < 0    l  = abs(kappa) -1   else   l  = kappa   end
    j = AngularJ64( abs(kappa) - 1//2 )
    # 
    JList = Basics.oplus(syma.J, j)
    for J in JList
        if    iseven(l)   push!( symtList, LevelSymmetry(J, syma.parity) )
        else              push!( symtList, LevelSymmetry(J, Basics.invertParity(syma.parity)) )
        end
    end
    return( symtList )         
end

"""
`AngularMomentum.allowedTotalSymmetries(symf::LevelSymmetry, mp2::EmMultipole, mp1::EmMultipole, symi::LevelSymmetry)`  
    ... to determine all allowed total symmetries J^P that can be constructed by coupling a multipole wave mp1 to the 
        initial symmetry symi, and which can be further coupled with mp2 to the final symmetry symf. 
        A list::Array{LevelSymmetry,1} of total symmetries is returned.
"""
function allowedTotalSymmetries(symf::LevelSymmetry, mp2::EmMultipole, mp1::EmMultipole, symi::LevelSymmetry) 
    symtList = LevelSymmetry[]
    waList   = oplus(symi.J, mp1.L);    wbList   = oplus(symf.J, mp2.L)
    for  wa in waList
        if  wa in wbList
            if      AngularMomentum.parityEmMultipolePi(symi.parity, mp1, Basics.plus)    &&
                    AngularMomentum.parityEmMultipolePi(symf.parity, mp2, Basics.plus)    push!( symtList, LevelSymmetry(wa, Basics.plus) )
            elseif  AngularMomentum.parityEmMultipolePi(symi.parity, mp1, Basics.minus)   &&
                    AngularMomentum.parityEmMultipolePi(symf.parity, mp2, Basics.minus)   push!( symtList, LevelSymmetry(wa, Basics.minus) )
            end
        end
    end

    return( symtList )         
end

"""
`AngularMomentum.bracket(jList::Array{AngularJ64,1})`  
    ... to compute the bracket [a, b, c, ... ] = (2a+1) * (2b+1) * (2b+1) * ... of the given angular momenta. 
        A value::Int64 is returned.
"""
function bracket(jList::Array{AngularJ64,1})
    value = 1;    for  j in jList    value = value * (Basics.twice(j) + 1)    end
    return( value )         
end

"""
`AngularMomentum.ChengI

+ (kapa::Int64, ma::AngularM64, kapb::Int64, mb::AngularM64, L::AngularJ64, M::AngularM64)` 
    ... evaluates the angular I (kappa m, kappa' m', LM) integral as defined by Cheng, NATO summerschool (198x), 
        Eq. (A4.5), and including the full magnetic (orientational) dependence. A value::Float64 is returned.
"""
function ChengI(kapa::Int64, ma::AngularM64, kapb::Int64, mb::AngularM64, L::AngularJ64, M::AngularM64)
    ja = Basics.subshell_j( Subshell(9, kapa) );    jb = Basics.subshell_j( Subshell(9, kapb) )   # Use principal QN n=9 arbitrarely here 
    la = Basics.subshell_l( Subshell(9, kapa) );    lb = Basics.subshell_l( Subshell(9, kapb) ) 
    #
    # Test for parity
    if  L.den != 1   error("stop a")                end 
    if  isodd( la + lb + L.num )    return( 0. )    end
    # 
    wa = AngularMomentum.phaseFactor([jb, +1, L, -1, ja]) * sqrt( (Basics.twice(jb)+1)*(Basics.twice(L)+1) / (4*pi*(Basics.twice(ja)+1) ) ) *
            AngularMomentum.ClebschGordan(jb, AngularM64(1//2), L, AngularM64(0), ja, AngularM64(1//2)) *
            AngularMomentum.ClebschGordan(jb, mb, L, M, ja, ma) 
    return( wa )
end

"""
+ (kapa::Int64, kapb::Int64, L::AngularJ64)` 
    ... evaluates the same angular I (kappa m, kappa' m', LM) integral but without the magnetic (orientational) dependence. 
        A value::Float64 is returned.
"""
function ChengI(kapa::Int64, kapb::Int64, L::AngularJ64)
    ja = Basics.subshell_j( Subshell(9, kapa) );    jb = Basics.subshell_j( Subshell(9, kapb) )   # Use principal QN n=9 arbitrarely here 
    la = Basics.subshell_l( Subshell(9, kapa) );    lb = Basics.subshell_l( Subshell(9, kapb) ) 
    #
    # Test for parity
    if  L.den != 1   error("stop a")                end 
    if  isodd( la + lb + L.num )    return( 0. )    end
    # 
    # There occurs here two changes with regard to Cheng formulas: 
    # i)  We now divide by sqrt(twoJ(jb)+1) instead of sqrt(twoJ(ja)+1) ... which is likely related to change emission - absorption
    # ii) The phase (-1)^(ja + L - jb) is replaced by (-1)^L  to get a proper phase between 1/2 --> 3/2 ME (compared to 3/2 --> 3/2) ...
    #     but which already comes from the Clebsch-Gordan
    wa = AngularMomentum.phaseFactor([ja, +1, L, -1, jb]) *  (-1)^L.num  * 
    ## wa = (-1)^L.num  * sqrt( (Basics.twice(jb)+1)*(Basics.twice(L)+1) / (4*pi) ) / sqrt(Basics.twice(jb)+1) *  
            AngularMomentum.ClebschGordan(jb, AngularM64(1//2), L, AngularM64(0), ja, AngularM64(1//2)) 
    return( wa )
end

"""
`AngularMomentum.ClebschGordan(ja, ma, jb, mb, Jab, Mab)`
    ... calculates the Clebsch-Gordan coefficient  <ja, ma, jb, mb; Jab, Mab> for given quantum numbers by
        a proper call to a Wigner 3-j symbol. A value::Float64 is returned.

        THIS IS THE ONE DEFINITION (14-Aug-2026). A second one, `ClebschGordan_old`, stood beside it and was
        called from EIGHTEEN places in four modules -- StrongField (12), StrongField-inc-hydrogenic (2),
        SpinAngular (3) and LandeZeeman (1). The two were the SAME FUNCTION written twice:
        both form  (-1)^(ja-jb+Mab) sqrt(2 Jab + 1) * Wigner_3j(ja, jb, Jab, ma, mb, -Mab), differing only in
        whether the phase is spelled out here or delegated to AngularMomentum.phaseFactor. Measured over
        35728 combinations with ja, jb <= 3 and Jab <= 6 across all m: the worst difference was EXACTLY zero,
        2348 of them non-vanishing, and both raised on precisely the same improper combinations. The twelve
        call sites were therefore migrated with NO phase factor and no other change.

        One difference did not survive, and it was the retired version's only advantage: phaseFactor RAISES
        on an improper combination (ja - jb + Mab not an integer), whereas the exponent is formed here as a
        Float64 and (-1)^x of a half-integer x raises a DomainError instead. Both stop; only the message
        differs.
"""
function ClebschGordan(ja, ma, jb, mb, Jab, Mab)
    mab = - Basics.twice(Mab) / 2
    pp  = (Basics.twice(ja) - Basics.twice(jb) + Basics.twice(Mab))/2
    cg  = (-1)^pp * sqrt(Basics.twice(Jab) + 1) * AngularMomentum.Wigner_3j(ja, jb, Jab, ma, mb, mab)
    return( cg )
end

"""
`AngularMomentum.CL_reduced_me(suba::Subshell, L::Int64, subb::Subshell)`  
    ... calculates the reduced matrix element of the C^L spherical tensor <suba || C^(L) || subb> in Grant's
        normalisation, sqrt((2j_a+1)(2j_b+1)); a value::Float64 is returned.

        THIS IS THE ONE DEFINITION (10-Aug-2026). It previously lacked the PARITY SELECTION RULE: the
        reduced matrix element of C^(L) vanishes unless l_a + l_b + L is even, and the 3-j symbol alone does
        NOT enforce that -- it only imposes the triangle and m-sum conditions. Without the rule the function
        returned nonzero values where the physical quantity must be zero, e.g. <2s||C^(1)||2s> = -0.816497.
        Where a different normalisation is wanted, divide EXPLICITLY at the call site: the convention that
        was carried by the former CL_reduced_me_rb is this quantity divided by sqrt(2j_a+1), with j_a taken
        from the FIRST argument as passed.
"""
function  CL_reduced_me(suba::Subshell, L::Int64, subb::Subshell)   
    la = Basics.subshell_l(suba);    ja2 = Basics.subshell_2j(suba);    
    lb = Basics.subshell_l(subb);    jb2 = Basics.subshell_2j(subb)
    rem(ja2+1, 2) != 0    &&    error("stop a")
    if  rem(la + lb + L, 2) != 0    return( 0. )    end

    redme = ((-1)^((ja2+1)/2)) * sqrt( (ja2+1)*(jb2+1) ) * 
            Wigner_3j(AngularJ64(ja2//2), AngularJ64(L), AngularJ64(jb2//2),  AngularM64(1//2), AngularM64(0), AngularM64(-1//2) )

    return( redme )
end



"""
`AngularMomentum.sigma_TtL_reduced_me(kapa::Int64, L::Int64, t::Int64, kapb::Int64)`
    ... calculates the reduced matrix element of the vector spherical harmonic tensor <kapa || sigma . T^(tL) || kapb>,
        following the RATIP convention (Grant's formalism); a value::Float64 is returned. L must equal t-1, t or t+1.
"""
function  sigma_TtL_reduced_me(kapa::Int64, L::Int64, t::Int64, kapb::Int64)
    if      L == t + 1
        redme = sqrt( (t+1.0)/(4pi) ) * (1.0 + (kapa+kapb)/(t+1.0)) * CL_reduced_me(Subshell(1,-kapa), t, Subshell(1,kapb))
    elseif  L == t
        redme = sqrt( (2t+1.0)/(4pi*t*(t+1.0)) ) * (kapb-kapa) * CL_reduced_me(Subshell(1,kapa), t, Subshell(1,kapb))
    elseif  L == t - 1
        redme = sqrt( t/(4pi) ) * (-1.0 + (kapa+kapb)/t) * CL_reduced_me(Subshell(1,-kapa), t, Subshell(1,kapb))
    else
        error("sigma_TtL_reduced_me(): L must equal t-1, t, or t+1.")
    end
    return( redme )
end


"""
`AngularMomentum.isAllowedMultipole(syma::LevelSymmetry, multipole::EmMultipole, symb::LevelSymmetry)`
    ... evaluates to true if the given multipole may connect the two level symmetries, and false otherwise.
"""
function  isAllowedMultipole(syma::LevelSymmetry, multipole::EmMultipole, symb::LevelSymmetry)
    if  isTriangle(syma.J, AngularJ64(multipole.L), symb.J)  &&
        parityEmMultipolePi(syma.parity, multipole, symb.parity)    return( true )
    else                                                            return( false )
    end
end

"""
`AngularMomentum.isTriangle(ja::AngularJ64, jb::AngularJ64, jc::AngularJ64)`  
    ... evaluates to true if Delta(ja,jb,jc) = 1, ie. if the angular momenta ja, jb and jc can couple to each other, 
        and false otherwise.
"""
function isTriangle(ja::AngularJ64, jb::AngularJ64, jc::AngularJ64) 
    if  ja.den == 1   ja2 = 2ja.num   else   ja2 = ja.num   end
    if  jb.den == 1   jb2 = 2jb.num   else   jb2 = jb.num   end
    if  jc.den == 1   jc2 = 2jc.num   else   jc2 = jc.num   end
    isodd(ja2 + jb2 + jc2)   &&    error("Angular momenta do no fullfill proper coupling rules; 2ja = $ja2, 2jb = $jb2, 2jc = $jc2") 
    if  ja2 + jb2 >= jc2     &&    jc2 + ja2 >= jb2   &&   jb2 + jc2 >= ja2    return( true )
    else                                                                       return( false )
    end
end

"""
`AngularMomentum.isTriangle(ja::Int64, jb::Int64, jc::Int64)`  
    ... evaluates to true if Delta(ja,jb,jc) = 1, ie. if the three integer (length) ja, jb and jc can form a triangle, 
        and false otherwise.
"""
function isTriangle(ja::Int64, jb::Int64, jc::Int64) 
    if  ja + jb >= jc     &&    jc + ja >= jb   &&   jb + jc >= ja    return( true )
    else                                                              return( false )
    end
end

"""
`AngularMomentum.JohnsonI(kapa::Int64, kapb::Int64, L::AngularJ64)` 
    ... evaluates the angular CL (kappa m, kappa' m', L M) integral as defined in his book. A value::Float64 is returned.
"""
function JohnsonI(kapa::Int64, kapb::Int64, L::AngularJ64)
    ja = Basics.subshell_j( Subshell(9, kapa) );    jb = Basics.subshell_j( Subshell(9, kapb) )   # Use principal QN n=9 arbitrarely here 
    la = Basics.subshell_l( Subshell(9, kapa) );    lb = Basics.subshell_l( Subshell(9, kapb) ) 
    #
    # Test for parity
    if  L.den != 1   error("stop a")                end 
    if  isodd( la + lb + L.num )    return( 0. )    end
    ## THE PHASE IS REAL, and was complex only through a slip (corrected 09-Aug-2026). The Racah phase here is
    ## (-1)^(j_b + 1/2); j_b is half-integer, so that exponent is an INTEGER and the phase is +-1, ALTERNATING
    ## with j_b. The previous expression raised (-1+0im) to `jb.num + 1/2` -- the NUMERATOR of j_b, i.e. 1, 3,
    ## 5, ... -- making the exponent half-integer, so it returned a CONSTANT -i for every j_b:
    ##      j_b      1/2     3/2     5/2     7/2
    ##      was       -i      -i      -i      -i
    ##      correct   -1      +1      -1      +1
    ## Two things followed. Every multipole matrix element came out purely imaginary, which is why these
    ## routines were documented as returning Float64 while in fact returning ComplexF64; and, more seriously,
    ## the j_b-dependent SIGN was lost, so contributions from orbitals of different j_b carried the same phase
    ## instead of opposite ones. For a single contributing orbital that is an overall phase and rates are
    ## unaffected -- which is why the hydrogen benchmarks never showed it -- but wherever such contributions
    ## interfere the relative signs were wrong.
    wa =  (-1)^Int((Basics.twice(jb)+1)/2) *sqrt((Basics.twice(jb)+1)*(Basics.twice(L)+1)*(Basics.twice(ja)+1)*(L.num+1) / 
            (4*pi*L.num ) ) *AngularMomentum.Wigner_3j(ja, L, jb, AngularM64(1//2), AngularM64(0),  AngularM64(-1//2))
    return( wa )
end

"""
`AngularMomentum.kappa_j(kappa::Int64)`  ... calculates the j::AngularJ64 value of a given kappa.
"""
function  kappa_j(kappa::Int64)  
    j = AngularJ64( abs(kappa) - 1//2 )
    return( j )
end

"""
`AngularMomentum.kappa_l(kappa::Int64)`  ... calculates the l::AngularJ64 value of a given kappa.
"""
function  kappa_l(kappa::Int64)  
    if  kappa < 0    l  = abs(kappa) -1   else   l  = kappa   end
    return( AngularJ64(l) )
end

"""
`AngularMomentum.j_values(j1::AngularJ64, j2::AngularJ64)`  
    ... returns a list of j-values that all fulfill the triangular condition delta(j1, j2, j) == 1;
        jList::Array{AngularJ64,1} is returned.
"""
function  j_values(j1::AngularJ64, j2::AngularJ64) 
    jList = AngularJ64[]
    
    if       j1.den == 1   &&   j2.den == 1    for  jx = abs(j1.num - j2.num):j1.num + j2.num         push!( jList, AngularJ64(jx))       end
    elseif   j1.den == 1   &&   j2.den == 2    for  jx = abs(2*j1.num - j2.num):2:2*j1.num + j2.num   push!( jList, AngularJ64(jx//2))    end
    elseif   j1.den == 2   &&   j2.den == 1    for  jx = abs(j1.num - 2*j2.num):2:j1.num + 2*j2.num   push!( jList, AngularJ64(jx//2))    end
    elseif   j1.den == 2   &&   j2.den == 2    for  jx = abs(j1.num - j2.num):2:j1.num + j2.num       push!( jList, AngularJ64(jx//2))    end
    else     error("stop a")
    end

    return( jList )
end

"""
`AngularMomentum.m_values(j::AngularJ64)`  ... returns a list of m-values for given j::AngularJ64.
"""
function  m_values(j::AngularJ64) 
    mList = AngularM64[]
    if       j.den == 1   for  jx = -j.num:j.num      push!( mList, AngularM64(jx))       end
    elseif   j.den == 2   for  jx = -j.num:2:j.num    push!( mList, AngularM64(jx//2))    end
    else     error("stop a")
    end

    return( mList )
end

"""
`AngularMomentum.oneJ(ja::AngularJ64)`  ... calculates ja; a (positive) value::Float64 is returned.
"""
function  oneJ(ja::AngularJ64)  
    if  ja.den  == 1    ja1 = 1.0 * ja.num   else   ja1 = ja.num / 2.   end
    return( ja1 )
end

"""
`AngularMomentum.oneM(ma::AngularM64)`  ... calculates ma; a value::Float64 is returned.
"""
function  oneM(ma::AngularM64)  
    if  ma.den  == 1    ma1 = 1.0 * ma.num   else   ma1 = ma.num / 2.   end
    return( ma1 )
end

"""
`AngularMomentum.parityEmMultipolePi(pa::Parity, multipole::EmMultipole, pb::Parity)`  
    ... evaluates to true if the given multipole fullfills the parity selection rule pi(a, multipole, b) = 1, 
        and false otherwise. This includes a proper test for both, electric and magnetic multipoles, based on 
        multipole.electric.
"""
function  parityEmMultipolePi(pa::Parity, multipole::EmMultipole, pb::Parity)
    # The the multipolarity into account to 'interprete' electric multipoles
    if       isodd(multipole.L)    Le = multipole.electric    else   Le = !(multipole.electric)   end

    if        pa == Basics.plus    &&   Le   &&   pb == Basics.minus        return( true )
    elseif    pa == Basics.minus   &&   Le   &&   pb == Basics.plus         return( true )
    elseif    pa == pb          &&   !(Le)                                  return( true )
    else                                                                    return( false )
    end
end

"""
`AngularMomentum.phaseFactor(list::Array{Any,1})` 
    ... checks and calculates the phase factor (-1)^(ja + mb -jc ...) that occur frequently in angular momentum theory; 
        a value +1. or -1. is returned. Use phaseFactor([ja::Union{AngularJ64,AngularM64), -1, mb::Union{AngularJ64,AngularM64), 
        ..., 1, jc::Union{AngularJ64,AngularM64)]) to specify the phase.
"""
function phaseFactor(list::Array{Any,1})
    if   iseven( length(list) )                                                   error("Wrong number of arguments.")   end
    for  i = 1:2:length(list)
            if  !(typeof(list[i]) == AngularJ64  || typeof(list[i]) == AngularM64)   error("Wrong type of argument $i")    end
            if  i == length(list)                                                    break                                 end
            if  list[i+1] != 1   &&  list[i+1] != -1                error("Wrong type of argument $(i+1); must be +-1")    end
    end
    #
    ja = list[1];        if  ja.den  == 1    jm2 = 2ja.num                  else   jm2 = ja.num                         end 
    for  i = 2:2:length(list)
        ja = list[i+1];   if  ja.den  == 1    jm2 = jm2 + list[i]*2ja.num    else   jm2 = jm2 + list[i]*ja.num           end 
    end
    #
    if rem(jm2,2) != 0    error("Improper combination of angular momenta")  end

    return( (-1.)^(jm2/2) )
end

"""
`AngularMomentum.phaseMultipole(x::ComplexF64, mp::EmMultipole)`  
    ... calculates (x)^p   with   mp = (L,p) and p = 0 (magnetic), p = 1 (electric).
        A  wa ::ComplexF64  is returned.
"""
function  phaseMultipole(x::ComplexF64, mp::EmMultipole)
    if  mp.electric   wa = x    else   wa = ComplexF64(1.0)   end 
    return( wa )
end
    

"""
`AngularMomentum.sigma_reduced_me(suba::Subshell, subb::Subshell)`
    ... calculates the reduced matrix element of the sigma^(1) spherical tensor <suba || sigma^(1) || subb>;
        a value::Float64 is returned.

        NOT IMPLEMENTED, and it now says so instead of returning zero (08-Aug-2026). Until then the body was a
        commented-out formula followed by `redme = 0.`, so the routine returned 0.0 for EVERY pair of subshells
        that passed the parity rule -- a plausible number, silently, from a core module. It has no callers today,
        which is the only reason it never did harm; a documented routine that answers every question with zero is
        a trap left lying in the path of whoever calls it first.

        THE COMMENTED-OUT FORMULA CANNOT SIMPLY BE RE-ENABLED: it references a rank `L` that is never defined in
        this method and appears in no argument, so it never ran in this form. Deciding what `L` was meant to be
        is a physics question, not a repair.

        WORKING SIBLINGS EXIST for the closely related magnetic-quantum-number-resolved elements,
        `sigma_reduced_me_ma` and `sigma_reduced_me_mb`, which ARE implemented and are used by
        `RadialIntegrals`; they are the place to look when this one is completed.
"""
function  sigma_reduced_me(suba::Subshell, subb::Subshell)
    la = Basics.subshell_l(suba);    ja2 = Basics.subshell_2j(suba);
    lb = Basics.subshell_l(subb);    jb2 = Basics.subshell_2j(subb)
    ## the parity rule is a genuine selection rule and its zero is a real answer, not a missing one
    if rem(la + lb, 2) != 0   return 0.     end

    error("\n\nAngularMomentum.sigma_reduced_me():  STOP -- <$suba || sigma^(1) || $subb> is NOT implemented.\n" *
          ">>> This routine returned 0.0 for every allowed pair of subshells until 08-Aug-2026; it now refuses\n"  *
          "    rather than hand back a plausible zero.\n"                                                          *
          ">>> The formula that stood here commented out cannot be re-enabled as it is: it references a rank L\n"  *
          "    that is defined neither in the method nor among its arguments.\n"                                   *
          ">>> See sigma_reduced_me_ma / sigma_reduced_me_mb, which ARE implemented, for the intended structure.\n")
end

function  sigma_reduced_me_ma(mkapa::Int64, kapb::Int64) 
    suba = Subshell(9,-mkapa);   ja2 = Basics.subshell_2j(suba);   ja = Basics.subshell_j(suba);   la = Basics.subshell_l(suba)
    subb = Subshell(9,  kapb);   jb2 = Basics.subshell_2j(subb);   jb = Basics.subshell_j(subb);   lb = Basics.subshell_l(subb)
    
    redme = 0.
    if  ja2-1 == lb
        redme = AngularMomentum.phaseFactor([AngularJ64(1//2), -1, ja]) * sqrt(6*(ja2+1)*(jb2+1)) *
                AngularMomentum.Wigner_6j(AngularJ64(1//2), ja, AngularJ64(ja2-1), jb, AngularJ64(1//2), AngularJ64(1))
    end

    return (redme)
end

function  sigma_reduced_me_mb(kapa::Int64, mkapb::Int64) 
    suba = Subshell(9,  kapa);   ja2 = Basics.subshell_2j(suba);   ja = Basics.subshell_j(suba);   la = Basics.subshell_l(suba)
    subb = Subshell(9,-mkapb);   jb2 = Basics.subshell_2j(subb);   jb = Basics.subshell_j(subb);   lb = Basics.subshell_l(subb)
    
    redme = 0.
    if  jb2-1 == la
        redme = AngularMomentum.phaseFactor([AngularJ64(3//2+la), +1, ja]) * sqrt(6*(ja2+1)*(jb2+1)) *
                AngularMomentum.Wigner_6j(AngularJ64(1//2), ja, AngularJ64(la), jb, AngularJ64(1//2), AngularJ64(1))
    end

    return (redme)
end

"""
`AngularMomentum.sphericalYlm(l::Int64, m::Int64, theta::Float64, phi::Float64)`
    ... calculates the spherical harmonics for low l-values explicitly. A value::Complex{Float64} is returned.
        Ylm = sqrt( (2*l+1) / (two*two*pi) ) * spherical_Clm(l,m,theta,phi).
"""
function sphericalYlm(l::Int64, m::Int64, theta::Float64, phi::Float64)
    # Note (28-Jul-2026): the azimuthal phase must be exp(i*m*phi), not exp(i*phi); the missing
    # "m*" factor here previously made every m != 0 call phi-dependent in the wrong way (e.g. a
    # perfectly octahedral point-charge lattice sum, which must vanish exactly for k=2 at every q,
    # picked up a spurious nonzero q=+-1,+-2 residual -- caught via CrystalField.multipoleLatticeSum,
    # see project memory). Fixed by including "m*" in the exponent.
    iphi  = m*phi*im
    # GSL's sf_legendre_sphPlm only accepts m >= 0; for m < 0 the standard relation
    # Y_l^{-|m|}(theta,phi) = (-1)^|m| * conj(Y_l^{|m|}(theta,phi)) requires this extra phase,
    # which is missing if abs(m) is used without it.
    phase = m < 0 ? (-1.0)^m : 1.0
    ylm   = phase * sf_legendre_sphPlm(l, abs(m), cos(theta)) * exp(iphi)

    return( ylm )
end

"""
`AngularMomentum.triangularDelta(ia2::Int64, ib2::Int64, ic2::Int64)`  
    ... calculates the tringular Delta(ja,jb,jc). The arguments in this integer function are i2a = 2*ja+1, ... 
        The result is 0 if the triangular condition failes and 1 otherwise. 
"""
function triangularDelta(ia2::Int64, ib2::Int64, ic2::Int64)    
    i = ib2 - ic2
    if  ia2 >= abs(i) + 1   &&   ia2 <= ib2 + ic2 - 1    return( 1 )   else    return( 0 )    end 
end

"""
`AngularMomentum.triangularDelta(ja::AngularJ64, jb::AngularJ64, jc::AngularJ64)`  
    ... calculates the tringular Delta(ja,jb,jc). The result is 0 if the triangular condition failes and 1 otherwise. 
"""
function triangularDelta(ja::AngularJ64, jb::AngularJ64, jc::AngularJ64)    
    if  abs(ja.num//ja.den - jb.num//jb.den) <= jc.num//jc.den <= ja.num//ja.den + jb.num//jb.den  return( 1 )   else    return( 0 )    end
end

"""
`AngularMomentum.Wigner_DFunction(j, p, q, alpha::Float64, beta::Float64, gamma::Float64)`  
    ... calculates the value of a Wigner D^j_pq (alpha, beta, gamma) for given quantum numbers and (Euler) angles 
        (alpha, beta, gamma). It makes use of the small Wigner d(beta) matrix as the key part.
        A value::Float64 is returned.
"""
function Wigner_DFunction(j, p, q, alpha::Float64, beta::Float64, gamma::Float64)
    wa = exp(-im*p*alpha -im*q*gamma) * AngularMomentum.Wigner_dmatrix(j, p, q, beta)
end

"""
`AngularMomentum.Wigner_dmatrix(jj, mmp, mm, beta::Float64)`  
    ... calculates the value of the small Wigner d^j_m',m (beta) for given quantum numbers and the angle beta.
        Wigner's formula is applied to evaluate the small Wigner matrix; a value::Float64 is returned.
"""
function Wigner_dmatrix(jj, mmp, mm, beta::Float64)
    j = Float64(jj);      mp = Float64(mmp);    m = Float64(mm)  
    if     mod(j+j + mp+mp + m+m, 2.0) != 0.       error("Inappropriate quantum numbers j=$j, mp=$mp, m=$m")   
    elseif abs(mp) > j    ||   abs(m) > j          return( 0. )
    end
    #
    factor = sqrt( factorial(j+mp) * factorial(j-mp) * factorial(j+m) * factorial(j-m) )
    wa = 0.
    for  s = 0:100
        if  j+m-s < 0  ||  j-mp-s < 0   break       end
        if  mp - m + s < 0              continue    end
        wa = wa + (-1)^(mp-m+s) * cos(beta/2.)^(2j+m-mp-2s) * sin(beta/2.)^(mp-m+2s) / 
                    ( factorial(j+m-s) * factorial(s) * factorial(mp-m+s) * factorial(j-mp-s) )
    end

    return( factor*wa )
end

"""
`AngularMomentum.Wigner_3j(a, b, c, m_a, m_b, m_c)`  
    ... calculates the value of a Wigner 3-j symbol for given quantum numbers as displayed in many texts on the theory of 
        angular momentum (see R. D. Cowan, The Theory of Atomic Structure and Spectra; University of California Press, 1981, p. 142); 
        it calls the corresponding function from the GNU Scientific Library. A value::Float64 is returned.
"""
function Wigner_3j(a, b, c, m_a, m_b, m_c)
    WignerSymbols.wigner3j(Basics.twice(a)/2.0,   Basics.twice(b)/2.0,   Basics.twice(c)/2.0,
                          Basics.twice(m_a)/2.0, Basics.twice(m_b)/2.0, Basics.twice(m_c)/2.0)
end


"""
`AngularMomentum.Wigner_6j(a, b, c, d, e, f)`  
    ... calculates the value of a Wigner 6-j symbol for given quantum numbers as displayed in many texts on the theory of 
        angular momentum (see R. D. Cowan, The Theory of Atomic Structure and Spectra; University of California Press, 1981, p. 142); 
        it calls the corresponding function from the GNU Scientific Library. A value::Float64 is returned.
"""
function Wigner_6j(a, b, c, d, e, f)
    WignerSymbols.wigner6j(Basics.twice(a)/2.0, Basics.twice(b)/2.0, Basics.twice(c)/2.0,
                            Basics.twice(d)/2.0, Basics.twice(e)/2.0, Basics.twice(f)/2.0)
end


    """
    `AngularMomentum.Wigner_9j(a, b, c, d, e, f, g, h, i)`  
        ... calculates the value of a Wigner 3-j symbol for given quantum numbers as displayed in many texts on the theory of 
            angular momentum (see R. D. Cowan, The Theory of Atomic Structure and Spectra; University of California Press, 1981, p. 142); 
            it calls the corresponding function from the GNU Scientific Library. A value::Float64 is returned.
    """
    function Wigner_9j( j11, j12, j13, j21, j22, j23, j31, j32, j33 )
        j11 = Basics.twice(j11) ;   j12 = Basics.twice(j12)   ;   j13 = Basics.twice(j13)
        j21 = Basics.twice(j21) ;   j22 = Basics.twice(j22)   ;   j23 = Basics.twice(j23)
        j31 = Basics.twice(j31) ;   j32 = Basics.twice(j32)   ;   j33 = Basics.twice(j33)
     
        kmin1 = abs(j11 - j33)  ;   kmin2 = abs(j32 - j21)    ;    kmin3 = abs(j23 - j12)
        kmax1 = j11 + j33       ;   kmax2 = j32 + j21         ;    kmax3 = j23 + j12
     
        if (kmin2 > kmin1) kmin1 = kmin2 end
        if (kmin3 > kmin1) kmin1 = kmin3 end
        if (kmax2 < kmax1) kmax1 = kmax2 end
        if (kmax3 < kmax1) kmax1 = kmax3 end
     
        kmin1 = kmin1 + 1
        kmax1 = kmax1 + 1
     
        nineJ = 0.0
     
        if (kmin1 > kmax1) return( nineJ ) end
     
        for k1 = kmin1:2:kmax1
           k = k1 - 1
           s1 = WignerSymbols.wigner6j(j11/2, j21/2, j31/2, j32/2, j33/2, k/2)
           s2 = WignerSymbols.wigner6j(j12/2, j22/2, j32/2, j21/2, k/2, j23/2)
           s3 = WignerSymbols.wigner6j(j13/2, j23/2, j33/2, k/2, j11/2, j12/2)
     
           p = (k+1) * (-1)^k
     
           nineJ += p * s1 * s2 * s3
        end
     
        return( nineJ )
    end


end # module
