#
# This file has mainly been contributed by Gediminas Gaigalas (gediminas.gaigalas@tfai.vu.lt)


"""
`SpinAngular.completlyReducedCfpByIndices(aIndex::Int64, bIndex::Int64)`  
    ... returns the (completely) reduced coefficient of fractional parentage 

        (j alpha Q J||| a^(qj)||| j alpha' Q' J')
        for j = 1/2, 3/2, 5/2, 7/2, and j=9/2 subshell states. A coeff::Float64 is returned.

        aIndex - is the state number of the bra subshell in quasispin representation;
        bIndex - is the state number of the ket subshell in quasispin representation.
"""
function completlyReducedCfpByIndices(aIndex::Int64, bIndex::Int64)
    #
    wa = [0, 1, 1]
    if aIndex <= bIndex   aIndexN = aIndex;    bIndexN = bIndex
    else                  aIndexN = bIndex;    bIndexN = aIndex
    end
    # Distinguish due to different j-values
    if       aIndexN == 1  &&  bIndexN == 2   wa = reducedCfp1half[aIndexN,bIndexN-1]
    elseif   aIndexN == 3  &&  bIndexN  > 3   &&  bIndexN  < 6
        wa = reducedCfp3half[aIndexN-2,bIndexN-3]
    elseif   aIndexN > 5   &&  aIndexN < 9   &&  bIndexN  > 8   &&  bIndexN  < 12
        wa = reducedCfp5half[aIndexN-5,bIndexN-8]
    elseif   aIndexN > 11  &&  aIndexN < 18  &&  bIndexN  > 17  &&  bIndexN  < 26
        wa = reducedCfp7half[aIndexN-11,bIndexN-17]
    elseif   aIndexN > 25  &&  aIndexN < 46  &&  bIndexN  > 45  &&  bIndexN  < 64
        wa = reducedCfp9half[aIndexN-25,bIndexN-45]
    end
    wb = wa[1] * sqrt(wa[2]/wa[3])
    if wa[1] != 0  &&  aIndex >= bIndex
        aT = qspaceTerms(aIndex);  bT = qspaceTerms(bIndex)
        wb = (-1)^Int64((Basics.twice(aT.Q)-Basics.twice(bT.Q)+Basics.twice(aT.J)-Basics.twice(bT.J)+Basics.twice(bT.j)-1)/2) *
                wb
    end
    return( wb )
end


"""
`SpinAngular.completelyReducedWkk(aIndex::Int64, bIndex::Int64, kq::Int64, kj::Int64)`
    ... returns the (completely) reduced coefficient of fractional parentage 

        (j alpha Q J||| W^(kq kj)||| j alpha' Q' J')
        for j = 1/2, 3/2, 5/2, 7/2, and 9/2 subshell states. A coeff::Float64 is returned

        aIndex - is the state number of the bra subshell in quasispin representation;
        bIndex - is the state number of the ket subshell in quasispin representation;
        kq     - is the rank kq in quasispin space q.
        kj     - is the rank kj in angular momentum j space.
"""
function completelyReducedWkk(aIndex::Int64, bIndex::Int64, kq::Int64, kj::Int64)
    #
    # Define a tuple for comparison in the subsequent lists; apply symmetries to reduced the list of coefficients
    #
    wa = [0, 1, 1]
    aT = qspaceTerms(aIndex);  bT = qspaceTerms(bIndex)
    if aT.j ==  bT.j < AngularJ64(9//2)  &&   aT.min_even ==  bT.min_even 
        # Distinguish due to different j-values and different kq, kj
        if aT.j == bT.j == AngularJ64(7//2) && kj > 1
            limits_1 = [0, 5, 9, 12, 14, 15]
            limits_2 = [0, 7, 13, 18, 22, 25, 27, 28]
            if aIndex <= bIndex   aIndexN = aIndex;    bIndexN = bIndex
            else                  aIndexN = bIndex;    bIndexN = aIndex
            end
        end
        if (kq,kj) == (0,0)   &&   aIndex == bIndex
            if aT.j == AngularJ64(1//2)       wa = W_00_one_half[aIndex]
            elseif aT.j == AngularJ64(3//2)   wa = W_00_three_half[aIndex-2]
            elseif aT.j == AngularJ64(5//2)   wa = W_00_five_half[aIndex-5]
            elseif aT.j == AngularJ64(7//2)   wa = W_00_seven_half[aIndex-11]
            end
        elseif (kq,kj) == (1,0)   &&   aIndex == bIndex
            if aT.j == AngularJ64(1//2)       wa = W_10_one_half[aIndex]
            elseif aT.j == AngularJ64(3//2)   wa = W_10_three_half[aIndex-2]
            elseif aT.j == AngularJ64(5//2)   wa = W_10_five_half[aIndex-5]
            elseif aT.j == AngularJ64(7//2)   wa = W_10_seven_half[aIndex-11]
            end
        elseif (kq,kj) == (0,1)   &&   aIndex == bIndex
            if aT.j == AngularJ64(1//2)       wa = W_01_one_half[aIndex]
            elseif aT.j == AngularJ64(3//2)   wa = W_01_three_half[aIndex-2]
            elseif aT.j == AngularJ64(5//2)   wa = W_01_five_half[aIndex-5]
            elseif aT.j == AngularJ64(7//2)   wa = W_01_seven_half[aIndex-11]
            end
        elseif (kq,kj) == (1,2) 
            if aT.j == AngularJ64(3//2) && aT.min_even == 3
                wa = W_12_three_half_odd[aIndex-2,bIndex-2]
            elseif aT.j == AngularJ64(3//2) && aT.min_even == 4
                wa = W_12_three_half_even[aIndex-3,bIndex-3]
            elseif aT.j == AngularJ64(5//2) && aT.min_even == 6
                wa = W_12_five_half_odd[aIndex-5,bIndex-5]
            elseif aT.j == AngularJ64(5//2) && aT.min_even == 9
                wa = W_12_five_half_even[aIndex-8,bIndex-8]
            elseif aT.j == AngularJ64(7//2) && aT.min_even == 12
                wa = W_12_seven_half_odd[limits_1[aIndexN-11]+bIndexN-11]
            elseif aT.j == AngularJ64(7//2) && aT.min_even == 18
                wa = W_12_seven_half_even[limits_2[aIndexN-17]+bIndexN-17]
            end
        elseif (kq,kj) == (0,3) 
            if aT.j == AngularJ64(3//2) && aT.min_even == 3
                wa = W_03_three_half_odd[aIndex-2,bIndex-2]
            elseif aT.j == AngularJ64(3//2) && aT.min_even == 4
                wa = W_03_three_half_even[aIndex-3,bIndex-3]
            elseif aT.j == AngularJ64(5//2) && aT.min_even == 6
                wa = W_03_five_half_odd[aIndex-5,bIndex-5]
            elseif aT.j == AngularJ64(5//2) && aT.min_even == 9
                wa = W_03_five_half_even[aIndex-8,bIndex-8]
            elseif aT.j == AngularJ64(7//2) && aT.min_even  == 12
                wa = W_03_seven_half_odd[limits_1[aIndexN-11]+bIndexN-11]
            elseif aT.j == AngularJ64(7//2) && aT.min_even == 18
                wa = W_03_seven_half_even[limits_2[aIndexN-17]+bIndexN-17]
            end
        elseif (kq,kj) == (1,4) 
            if aT.j == AngularJ64(5//2) && aT.min_even == 6
                wa = W_14_five_half_odd[aIndex-5,bIndex-5]
            elseif aT.j == AngularJ64(5//2) && aT.min_even == 9
                wa = W_14_five_half_even[aIndex-8,bIndex-8]
            elseif aT.j == AngularJ64(7//2) && aT.min_even  == 12
                wa = W_14_seven_half_odd[limits_1[aIndexN-11]+bIndexN-11]
            elseif aT.j == AngularJ64(7//2) && aT.min_even == 18
                wa = W_14_seven_half_even[limits_2[aIndexN-17]+bIndexN-17]
            end
        elseif (kq,kj) == (0,5) 
            if aT.j == AngularJ64(5//2) && aT.min_even == 6
                wa = W_05_five_half_odd[aIndex-5,bIndex-5]
            elseif aT.j == AngularJ64(5//2) && aT.min_even == 9
                wa = W_05_five_half_even[aIndex-8,bIndex-8]
            elseif aT.j == AngularJ64(7//2) && aT.min_even  == 12
                wa = W_05_seven_half_odd[limits_1[aIndexN-11]+bIndexN-11]
            elseif aT.j == AngularJ64(7//2) && aT.min_even == 18
                wa = W_05_seven_half_even[limits_2[aIndexN-17]+bIndexN-17]
            end
        elseif (kq,kj) == (1,6) 
            if aT.j == AngularJ64(7//2) && aT.min_even  == 12
                wa = W_16_seven_half_odd[limits_1[aIndexN-11]+bIndexN-11]
            elseif aT.j == AngularJ64(7//2) && aT.min_even == 18
                wa = W_16_seven_half_even[limits_2[aIndexN-17]+bIndexN-17]
            end
        elseif (kq,kj) == (0,7) 
            if aT.j == AngularJ64(7//2) && aT.min_even  == 12
                wa = W_07_seven_half_odd[limits_1[aIndexN-11]+bIndexN-11]
            elseif aT.j == AngularJ64(7//2) && aT.min_even == 18
                wa = W_07_seven_half_even[limits_2[aIndexN-17]+bIndexN-17]
            end
        end
        wb = wa[1] * sqrt(wa[2]/wa[3])
        if aT.j == AngularJ64(7//2) && wa[1] != 0  &&  aIndex >= bIndex && kj > 1
            wb = (-1)^Int64((Basics.twice(aT.Q)-Basics.twice(bT.Q)+Basics.twice(aT.J)-Basics.twice(bT.J))/2) * wb
        end
    elseif aT.j ==  bT.j == AngularJ64(9//2)  &&   aT.min_even ==  bT.min_even 
        wb  = completelyReducedWkk_calculate(aIndex, bIndex, kq, kj)
    else
        wb = 0.0
    end
    return( wb )
end


"""
`SpinAngular.completelyReducedWkk_calculate(aIndex::Int64, bIndex::Int64, kq::Int64, kj::Int64)`
    ... calculate the (completely) reduced coefficient of fractional parentage 

        (j alpha Q J||| W^(kq kj)||| j alpha' Q' J')
        for j = 9/2 subshell states. A coeff::Float64 is returned.

        aIndex - is the state number of the bra subshell in quasispin representation;
        bIndex - is the state number of the ket subshell in quasispin representation.
        kq     - is the rank kq in quasispin space q.
        kj     - is the rank kj in angular momentum j space.
"""
function completelyReducedWkk_calculate(aIndex::Int64, bIndex::Int64, kq::Int64, kj::Int64)
    #
    # Define a tuple for comparison in the subsequent lists; apply symmetries to reduced the list of coefficients
    #
    wb = 0.0
    aT = qspaceTerms(aIndex);      bT = qspaceTerms(bIndex)
    if aT.min_even ==  bT.min_even
        coeff = 0.0
        for rterm = aT.min_odd:aT.max_odd
            rT = qspaceTerms(rterm)
            coeffJ = AngularMomentum.Wigner_6j(aT.j, aT.j, kj, bT.J, aT.J, rT.J)
            if coeffJ != 0.0
                coeffQ = AngularMomentum.Wigner_6j(AngularJ64(1//2), AngularJ64(1//2), kq, bT.Q, aT.Q, rT.Q)
                if coeffQ != 0.0
                    coeff = coeff + completlyReducedCfpByIndices(aIndex, rterm) * 
                                    completlyReducedCfpByIndices(rterm, bIndex) * coeffJ * coeffQ
                end
            end
        end
        wb = (-1)^Int64((Basics.twice(aT.Q)+Basics.twice(aT.J)+Basics.twice(bT.Q)+Basics.twice(bT.J)+
                            Basics.twice(kq)+Basics.twice(kj))/2) * coeff * sqrt((2*kq+1.0) * (2*kj+1.0))
    end
    return( wb )
end


"""
`SpinAngular.getTermNumber(j::AngularJ64, N::Int64, Q::AngularJ64, J::AngularJ64)`  
    ... returns the internal index for a subshell term which is given by its angular momentum j, the quasispin quantum 
        number Q as well as the total subshell angular momentum J.

        j - is angular momentum j for the subshell;
        N - is number of electrons in the subshell;
        Q - is the subshell total quasispin Q;
        J - is the subshell total angular momentum J.
"""
function  getTermNumber(j::AngularJ64, N::Int64, Q::AngularJ64, J::AngularJ64)
    no_min = [1, 0, 3, 0, 6, 0,12, 0,26]
    no_max = [2, 0, 5, 0,11, 0,25, 0,63]
    I = 0
    if j <= AngularJ64(9//2)
        p_min = no_min[Basics.twice(j)];   p_max = no_max[Basics.twice(j)]
    end
    if j < AngularJ64(9//2)
        for run = p_min:p_max
            if Q == qspaceTerms(run).Q 
                    if J == qspaceTerms(run).J   I = run;  return (I);  end
            end
        end
    elseif j == AngularJ64(9//2)
        if N == 0 || N == 1 || N == 2 
            for run = p_min:p_max
                if Q == qspaceTerms(run).Q 
                    if J == qspaceTerms(run).J   I = run;  return (I);  end
                end
            end
        end
    else
        if N == 0 || N == 1 || N == 2   I = ((Basics.twice(j)*1000) + Basics.twice(Q))*1000 + Basics.twice(J)
        else    error("SpinAngular.getTermNumber: please, contact with G. Gaigalas")
        end
    end
    return( I )
end


"""
`SpinAngular.getTermNumber(j::AngularJ64, Q::AngularJ64, J::AngularJ64, Nr::Int64)`  
    ... returns the internal index for a subshell term which is given by its angular momentum j, the quasispin quantum 
        number Q as well as the total subshell angular momentum J.

        j  - is angular momentum j for the subshell;
        Q  - is the subshell total quasispin Q;
        J  - is the subshell total angular momentum J;
        Nr - is additional parameter for j = 9/2 for one by one indentification. It can be any integere for j = 1/2, 3/2, 5/2, 7/2.
"""
function  getTermNumber(j::AngularJ64, Q::AngularJ64, J::AngularJ64, Nr::Int64)
    no_min = [1, 0, 3, 0, 6, 0,12, 0,26]
    no_max = [2, 0, 5, 0,11, 0,25, 0,63]
    I = 0
    if j <= AngularJ64(9//2)
        p_min = no_min[Basics.twice(j)];   p_max = no_max[Basics.twice(j)]
    end
    if j < AngularJ64(9//2)
        for run = p_min:p_max
            if Q == qspaceTerms(run).Q 
                    if J == qspaceTerms(run).J   I = run;  return (I);  end
            end
        end
    elseif j == AngularJ64(9//2)
        for run = p_min:p_max
            if Q == qspaceTerms(run).Q 
                if J == qspaceTerms(run).J
                    if Nr == qspaceTerms(run).Nr   I = run;  return (I);  end
                end
            end
        end
    else    error("SpinAngular.getTermNumber supporte j= 1/2, 3/2, 5/2, 7/2, 9/2 subshells")
    end
    return( I )
end


"""
` + (jx, Qx, Jx, Nr::Int64)` ... to returns the same but without type binding apart from N.
"""
function  getTermNumber(jx, Qx, Jx, Nr::Int64)
    if  typeof(jx) == AngularJ64      j = jx   else  jz = convert(Int64, 2 * jx);   j  = AngularJ64(jz//2)   end
    if  typeof(Qx) == AngularJ64      Q = Qx   else  Qz = convert(Int64, 2 * Qx);   Q  = AngularJ64(Qz//2)   end
    if  typeof(Jx) == AngularJ64      J = Jx   else  Jz = convert(Int64, 2 * Jx);   J  = AngularJ64(Jz//2)   end
    return( getTermNumber(j::AngularJ64, Q::AngularJ64, J::AngularJ64, Nr::Int64) )
end

    
"""
`SpinAngular.qspaceTerms(iIndex::Int64)`  
    ... returns a list of Q-space term for the given state number; a term::SpinAngular.QspaceTerm,1 is returned.

        iIndex - is the state number of the subshell in quasispin representation.
"""
function  qspaceTerms(iIndex::Int64)
    wa = QTerm[iIndex] 
    return( wa )
end

#######################################################################################################################
#######################################################################################################################

"""
`SpinAngular.QTerm[]`
    ... set-up the list of Q-space term for the given Term number.
"""
    # First define some short-cuts to simplify the subsequent input
    h0  = AngularJ64(0//2);   h1  = AngularJ64(1//2);   h2  = AngularJ64(2//2);   h3  = AngularJ64(3//2);   h4  = AngularJ64(4//2)    
    h5  = AngularJ64(5//2);   h6  = AngularJ64(6//2);   h7  = AngularJ64(7//2);   h8  = AngularJ64(8//2);   h9  = AngularJ64(9//2)    
    h10 = AngularJ64(10//2);  h11 = AngularJ64(11//2);  h12 = AngularJ64(12//2);  h13 = AngularJ64(13//2);  h14 = AngularJ64(14//2)
    h15 = AngularJ64(15//2);  h16 = AngularJ64(16//2);  h17 = AngularJ64(17//2);  h18 = AngularJ64(18//2);  h19 = AngularJ64(19//2)
    h20 = AngularJ64(20//2);  h21 = AngularJ64(21//2);  h22 = AngularJ64(22//2);  h23 = AngularJ64(23//2);  h24 = AngularJ64(24//2)
    h25 = AngularJ64(25//2);  h26 = AngularJ64(26//2);  h27 = AngularJ64(27//2)
    #
    # Distinguish due to different j-values
    # j = 1/2)
    QTerm =  [ QspaceTerm(h1, h0, h1, 0, 2, 2, 1, 1), QspaceTerm(h1, h1, h0, 0, 1, 1, 2, 2),
    # j = 3/2
                QspaceTerm(h3, h1, h3, 0, 4, 5, 3, 3), QspaceTerm(h3, h2, h0, 0, 3, 3, 4, 5), QspaceTerm(h3, h0, h4, 0, 3, 3, 4, 5),
    # j = 5/2
                QspaceTerm(h5, h2, h5, 0, 9,11, 6, 8), QspaceTerm(h5, h0, h3, 0, 9,11, 6, 8), QspaceTerm(h5, h0, h9, 0, 9,11, 6, 8),
                QspaceTerm(h5, h3, h0, 0, 6, 8, 9,11), QspaceTerm(h5, h1, h4, 0, 6, 8, 9,11), QspaceTerm(h5, h1, h8, 0, 6, 8, 9,11),
    # j = 7/2
                QspaceTerm(h7, h1, h3, 0,18,25,12,17), QspaceTerm(h7, h1, h5, 0,18,25,12,17), QspaceTerm(h7, h3, h7, 0,18,25,12,17),
                QspaceTerm(h7, h1, h9, 0,18,25,12,17), QspaceTerm(h7, h1,h11, 0,18,25,12,17), QspaceTerm(h7, h1,h15, 0,18,25,12,17),
                QspaceTerm(h7, h4, h0, 0,12,17,18,25), QspaceTerm(h7, h2, h4, 0,12,17,18,25), QspaceTerm(h7, h0, h4, 0,12,17,18,25),
                QspaceTerm(h7, h2, h8, 0,12,17,18,25), QspaceTerm(h7, h0, h8, 0,12,17,18,25), QspaceTerm(h7, h0,h10, 0,12,17,18,25),
                QspaceTerm(h7, h2,h12, 0,12,17,18,25), QspaceTerm(h7, h0,h16, 0,12,17,18,25),
    # j = 9/2
                QspaceTerm(h9, h0, h1, 0,46,63,26,45), QspaceTerm(h9, h2, h3, 0,46,63,26,45), QspaceTerm(h9, h2, h5, 0,46,63,26,45),
                QspaceTerm(h9, h0, h5, 0,46,63,26,45), QspaceTerm(h9, h2, h7, 0,46,63,26,45), QspaceTerm(h9, h0, h7, 0,46,63,26,45),
                QspaceTerm(h9, h4, h9, 0,46,63,26,45), QspaceTerm(h9, h2, h9, 0,46,63,26,45), QspaceTerm(h9, h0, h9, 0,46,63,26,45),
                QspaceTerm(h9, h2,h11, 0,46,63,26,45), QspaceTerm(h9, h0,h11, 0,46,63,26,45), QspaceTerm(h9, h2,h13, 0,46,63,26,45),
                QspaceTerm(h9, h0,h13, 0,46,63,26,45), QspaceTerm(h9, h2,h15, 0,46,63,26,45), QspaceTerm(h9, h0,h15, 0,46,63,26,45),
                QspaceTerm(h9, h2,h17, 0,46,63,26,45), QspaceTerm(h9, h0,h17, 0,46,63,26,45), QspaceTerm(h9, h0,h19, 0,46,63,26,45),
                QspaceTerm(h9, h2,h21, 0,46,63,26,45), QspaceTerm(h9, h0,h25, 0,46,63,26,45), QspaceTerm(h9, h5, h0, 0,26,45,46,63),
                QspaceTerm(h9, h1, h0, 0,26,45,46,63), QspaceTerm(h9, h3, h4, 0,26,45,46,63), QspaceTerm(h9, h1, h4, 0,26,45,46,63),
                QspaceTerm(h9, h1, h6, 0,26,45,46,63), QspaceTerm(h9, h3, h8, 0,26,45,46,63), QspaceTerm(h9, h1, h8, 1,26,45,46,63),
                QspaceTerm(h9, h1, h8, 2,26,45,46,63), QspaceTerm(h9, h1,h10, 0,26,45,46,63), QspaceTerm(h9, h3,h12, 0,26,45,46,63),
                QspaceTerm(h9, h1,h12, 1,26,45,46,63), QspaceTerm(h9, h1,h12, 2,26,45,46,63), QspaceTerm(h9, h1,h14, 0,26,45,46,63),
                QspaceTerm(h9, h3,h16, 0,26,45,46,63), QspaceTerm(h9, h1,h16, 0,26,45,46,63), QspaceTerm(h9, h1,h18, 0,26,45,46,63),
                QspaceTerm(h9, h1,h20, 0,26,45,46,63), QspaceTerm(h9, h1,h24, 0,26,45,46,63) ]


"""
`SpinAngular.reducedCfp1half[]`
    ... set-up the reduced coefficients of fractional parentage (rcfp) for j=1/2
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
reducedCfp1half = Matrix{Tuple{Int64,Int64,Int64}}(undef,1,1);    for i=1:1  for j=1:1  reducedCfp1half[i,j] = (0,0,0)   end  end  


"""
`SpinAngular.reducedCfp3half[]`
    ... set-up the reduced coefficients of fractional parentage (rcfp) for j=3/2
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
reducedCfp3half = Matrix{Tuple{Int64,Int64,Int64}}(undef,1,2);    for i=1:1  for j=1:2  reducedCfp3half[i,j] = (0,0,0)   end  end


"""
`SpinAngular.reducedCfp5half[]`
    ... set-up the reduced coefficients of fractional parentage (rcfp) for j=5/2
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
reducedCfp5half = Matrix{Tuple{Int64,Int64,Int64}}(undef,6,6);    for i=1:3  for j=1:3  reducedCfp5half[i,j] = (0,0,0)   end  end


"""
`SpinAngular.reducedCfp7half[]`
    ... set-up the reduced coefficients of fractional parentage (rcfp) for j=7/2
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
reducedCfp7half = Matrix{Tuple{Int64,Int64,Int64}}(undef,6,8);  for i=1:6 for j=1:8 reducedCfp7half[i,j] = (0,0,0)   end  end


"""
`SpinAngular.reducedCfp9half[]`
    ... set-up the reduced coefficients of fractional parentage (rcfp) for j=9/2
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
reducedCfp9half = Matrix{Tuple{Int64,Int64,Int64}}(undef,20,18);  for i=1:20 for j=1:18 reducedCfp9half[i,j] = (0,0,0)   end  end

    reducedCfp1half[1,1] = (-1,      4,     1)
    
    reducedCfp3half[1,1] = (-1,     12,     1); reducedCfp3half[1,2] = (-1,     20,     1)

    reducedCfp5half[1,1] = (-1,     24,     1); reducedCfp5half[2,1] = ( 0,      1,     1)
    reducedCfp5half[3,1] = ( 0,      1,     1); reducedCfp5half[1,2] = (-1,     30,     1)
    reducedCfp5half[2,2] = ( 1,    120,     7); reducedCfp5half[3,2] = (-1,     90,     7)
    reducedCfp5half[1,3] = (-1,     54,     1); reducedCfp5half[2,3] = (-1,     48,     7)
    reducedCfp5half[3,3] = ( 1,    330,     7)

    reducedCfp7half[1,1] = ( 0,      1,     1); reducedCfp7half[2,1] = ( 0,      1,     1)
    reducedCfp7half[3,1] = (-1,     40,     1); reducedCfp7half[4,1] = ( 0,      1,     1)
    reducedCfp7half[5,1] = ( 0,      1,     1); reducedCfp7half[6,1] = ( 0,      1,     1)
    reducedCfp7half[1,2] = (-1,     54,     7); reducedCfp7half[2,2] = (-1,     33,     1)
    reducedCfp7half[3,2] = (-1,     40,     1); reducedCfp7half[4,2] = (-1,     65,     7)
    reducedCfp7half[5,2] = ( 1,     30,     1); reducedCfp7half[6,2] = ( 0,      1,     1)
    reducedCfp7half[1,3] = ( 1,     88,     7); reducedCfp7half[2,3] = (-1,      1,     1)
    reducedCfp7half[3,3] = ( 0,      1,     1); reducedCfp7half[4,3] = (-1,   1755,    77)
    reducedCfp7half[5,3] = (-1,     40,    11); reducedCfp7half[6,3] = ( 0,      1,     1)
    reducedCfp7half[1,4] = ( 1,    198,     7); reducedCfp7half[2,4] = (-1,     36,    11)
    reducedCfp7half[3,4] = (-1,     72,     1); reducedCfp7half[4,4] = ( 1,   4500,    77)
    reducedCfp7half[5,4] = (-1,    234,    11); reducedCfp7half[6,4] = (-1,    360,    11)
    reducedCfp7half[1,5] = (-1,     78,    35); reducedCfp7half[2,5] = ( 1,    156,     5)
    reducedCfp7half[3,5] = ( 0,      1,     1); reducedCfp7half[4,5] = ( 1,    108,    91)
    reducedCfp7half[5,5] = (-1,     30,     1); reducedCfp7half[6,5] = ( 1,     96,    13)
    reducedCfp7half[1,6] = ( 1,     66,     5); reducedCfp7half[2,6] = ( 1,     49,     5)
    reducedCfp7half[3,6] = ( 0,      1,     1); reducedCfp7half[4,6] = ( 1,     27,     1)
    reducedCfp7half[5,6] = ( 1,    130,     7); reducedCfp7half[6,6] = ( 1,    136,     7)
    reducedCfp7half[1,7] = ( 0,      1,     1); reducedCfp7half[2,7] = ( 1,    195,    11)
    reducedCfp7half[3,7] = (-1,    104,     1); reducedCfp7half[4,7] = (-1,    245,    11)
    reducedCfp7half[5,7] = (-1,    624,    11); reducedCfp7half[6,7] = ( 1,   1224,    11)
    reducedCfp7half[1,8] = ( 0,      1,     1); reducedCfp7half[2,8] = ( 0,      1,     1)
    reducedCfp7half[3,8] = ( 0,      1,     1); reducedCfp7half[4,8] = ( 1,   2720,   143)
    reducedCfp7half[5,8] = (-1,   2448,    77); reducedCfp7half[6,8] = (-1,   7752,    91)

    # rcfp for term_jj(46)
    reducedCfp9half[ 1,1] = ( 0,      1,     1); reducedCfp9half[ 2,1] = ( 0,      1,     1)
    reducedCfp9half[ 3,1] = ( 0,      1,     1); reducedCfp9half[ 4,1] = ( 0,      1,     1)
    reducedCfp9half[ 5,1] = ( 0,      1,     1); reducedCfp9half[ 6,1] = ( 0,      1,     1)
    reducedCfp9half[ 7,1] = (-1,     60,     1); reducedCfp9half[ 8,1] = ( 0,      1,     1)
    reducedCfp9half[ 9,1] = ( 0,      1,     1); reducedCfp9half[10,1] = ( 0,      1,     1)
    reducedCfp9half[11,1] = ( 0,      1,     1); reducedCfp9half[12,1] = ( 0,      1,     1)
    reducedCfp9half[13,1] = ( 0,      1,     1); reducedCfp9half[14,1] = ( 0,      1,     1)
    reducedCfp9half[15,1] = ( 0,      1,     1); reducedCfp9half[16,1] = ( 0,      1,     1)
    reducedCfp9half[17,1] = ( 0,      1,     1); reducedCfp9half[18,1] = ( 0,      1,     1)
    reducedCfp9half[19,1] = ( 0,      1,     1); reducedCfp9half[20,1] = ( 0,      1,     1)
    # rcfp for term_jj(47)
    reducedCfp9half[ 1,2] = ( 0,      1,     1); reducedCfp9half[ 2,2] = ( 0,      1,     1)
    reducedCfp9half[ 3,2] = ( 0,      1,     1); reducedCfp9half[ 4,2] = ( 0,      1,     1)
    reducedCfp9half[ 5,2] = ( 0,      1,     1); reducedCfp9half[ 6,2] = ( 0,      1,     1)
    reducedCfp9half[ 7,2] = ( 0,      1,     1); reducedCfp9half[ 8,2] = (-1,     12,     1)
    reducedCfp9half[ 9,2] = (-1,      8,     1); reducedCfp9half[10,2] = ( 0,      1,     1)
    reducedCfp9half[11,2] = ( 0,      1,     1); reducedCfp9half[12,2] = ( 0,      1,     1)
    reducedCfp9half[13,2] = ( 0,      1,     1); reducedCfp9half[14,2] = ( 0,      1,     1)
    reducedCfp9half[15,2] = ( 0,      1,     1); reducedCfp9half[16,2] = ( 0,      1,     1)
    reducedCfp9half[17,2] = ( 0,      1,     1); reducedCfp9half[18,2] = ( 0,      1,     1)
    reducedCfp9half[19,2] = ( 0,      1,     1); reducedCfp9half[20,2] = ( 0,      1,     1)
    # rcfp for term_jj(48)
    reducedCfp9half[ 1,3] = ( 0,      1,     1); reducedCfp9half[ 2,3] = ( 0,      1,     1)
    reducedCfp9half[ 3,3] = (-1,     20,     1); reducedCfp9half[ 4,3] = ( 0,      1,     1)
    reducedCfp9half[ 5,3] = (-1,   1664,    33); reducedCfp9half[ 6,3] = ( 0,      1,     1)
    reducedCfp9half[ 7,3] = (-1,     50,     1); reducedCfp9half[ 8,3] = (-1,    130,    33)
    reducedCfp9half[ 9,3] = ( 0,      1,     1); reducedCfp9half[10,3] = (-1,    272,    11)
    reducedCfp9half[11,3] = ( 0,      1,     1); reducedCfp9half[12,3] = (-1,    560,    11)
    reducedCfp9half[13,3] = ( 0,      1,     1); reducedCfp9half[14,3] = ( 0,      1,     1)
    reducedCfp9half[15,3] = ( 0,      1,     1); reducedCfp9half[16,3] = ( 0,      1,     1)
    reducedCfp9half[17,3] = ( 0,      1,     1); reducedCfp9half[18,3] = ( 0,      1,     1)
    reducedCfp9half[19,3] = ( 0,      1,     1); reducedCfp9half[20,3] = ( 0,      1,     1) 
    # rcfp for term_jj(49)
    reducedCfp9half[ 1,4] = ( 0,      1,     1); reducedCfp9half[ 2,4] = ( 0,      1,     1)
    reducedCfp9half[ 3,4] = (-1,    130,     7); reducedCfp9half[ 4,4] = (-1,     48,     7)
    reducedCfp9half[ 5,4] = ( 1,     44,    21); reducedCfp9half[ 6,4] = (-1,    340,    77)
    reducedCfp9half[ 7,4] = ( 0,      1,     1); reducedCfp9half[ 8,4] = ( 1,    140,    33)
    reducedCfp9half[ 9,4] = (-1,     70,    11); reducedCfp9half[10,4] = ( 1,   3808,   143)
    reducedCfp9half[11,4] = (-1,   2280,   143); reducedCfp9half[12,4] = (-1,    110,    13)
    reducedCfp9half[13,4] = (-1,    918,   143); reducedCfp9half[14,4] = ( 0,      1,     1)
    reducedCfp9half[15,4] = ( 0,      1,     1); reducedCfp9half[16,4] = ( 0,      1,     1)
    reducedCfp9half[17,4] = ( 0,      1,     1); reducedCfp9half[18,4] = ( 0,      1,     1)
    reducedCfp9half[19,4] = ( 0,      1,     1); reducedCfp9half[20,4] = ( 0,      1,     1)
    # rcfp for term_jj(50)
    reducedCfp9half[ 1,5] = ( 0,      1,     1); reducedCfp9half[ 2,5] = (-1,     63,     5)
    reducedCfp9half[ 3,5] = (-1,    312,    55); reducedCfp9half[ 4,5] = (-1,      4,    11)
    reducedCfp9half[ 5,5] = ( 1,      9,    11); reducedCfp9half[ 6,5] = ( 1,    255,    11)
    reducedCfp9half[ 7,5] = ( 0,      1,     1); reducedCfp9half[ 8,5] = ( 1,   5040,   143)
    reducedCfp9half[ 9,5] = ( 1,    840,   143); reducedCfp9half[10,5] = (-1,   1428,   143)
    reducedCfp9half[11,5] = (-1,     95,   143); reducedCfp9half[12,5] = ( 1,    504,   715)
    reducedCfp9half[13,5] = (-1,   2856,   143); reducedCfp9half[14,5] = ( 1,  13566,   715)
    reducedCfp9half[15,5] = (-1,    850,   143); reducedCfp9half[16,5] = ( 0,      1,     1)
    reducedCfp9half[17,5] = ( 0,      1,     1); reducedCfp9half[18,5] = ( 0,      1,     1)
    reducedCfp9half[19,5] = ( 0,      1,     1); reducedCfp9half[20,5] = ( 0,      1,     1)
    # rcfp for term_jj(51)
    reducedCfp9half[ 1,6] = ( 0,      1,     1); reducedCfp9half[ 2,6] = (-1,    384,    11)
    reducedCfp9half[ 3,6] = ( 1,    156,    11); reducedCfp9half[ 4,6] = ( 0,      1,     1)
    reducedCfp9half[ 5,6] = (-1,   1920,   143); reducedCfp9half[ 6,6] = ( 0,      1,     1)
    reducedCfp9half[ 7,6] = (-1,     90,     1); reducedCfp9half[ 8,6] = ( 1,   7350,   143)
    reducedCfp9half[ 9,6] = ( 0,      1,     1); reducedCfp9half[10,6] = ( 1,   8160,   143)
    reducedCfp9half[11,6] = ( 0,      1,     1); reducedCfp9half[12,6] = ( 1,   1512,   143)
    reducedCfp9half[13,6] = ( 0,      1,     1); reducedCfp9half[14,6] = (-1,   3648,   143)
    reducedCfp9half[15,6] = ( 0,      1,     1); reducedCfp9half[16,6] = (-1,   9000,   143)
    reducedCfp9half[17,6] = ( 0,      1,     1); reducedCfp9half[18,6] = ( 0,      1,     1)
    reducedCfp9half[19,6] = ( 0,      1,     1); reducedCfp9half[20,6] = ( 0,      1,     1)
    # rcfp for term_jj(52)
    reducedCfp9half[ 1,7] = (-1,   2184,   253); reducedCfp9half[ 2,7] = (-1,     63,    23)
    reducedCfp9half[ 3,7] = (-1,  59904, 19481); reducedCfp9half[ 4,7] = (-1, 302460, 19481)
    reducedCfp9half[ 5,7] = ( 1,   5265,   161); reducedCfp9half[ 6,7] = ( 1, 848691,253253)
    reducedCfp9half[ 7,7] = ( 0,      1,     1); reducedCfp9half[ 8,7] = (-1, 145152, 36179)
    reducedCfp9half[ 9,7] = ( 1, 217728, 36179); reducedCfp9half[10,7] = ( 1,1049580, 36179)
    reducedCfp9half[11,7] = ( 1, 287337, 36179); reducedCfp9half[12,7] = ( 1,   5184,   299)
    reducedCfp9half[13,7] = ( 1, 261120, 36179); reducedCfp9half[14,7] = ( 1, 691866, 36179)
    reducedCfp9half[15,7] = (-1,    750, 36179); reducedCfp9half[16,7] = ( 0,      1,     1)
    reducedCfp9half[17,7] = (-1,  76608,  3289); reducedCfp9half[18,7] = ( 0,      1,     1)
    reducedCfp9half[19,7] = ( 0,      1,     1); reducedCfp9half[20,7] = ( 0,      1,     1)
    # rcfp for term_jj(53)
    reducedCfp9half[ 1,8] = ( 1,   1224,  1265); reducedCfp9half[ 2,8] = ( 1,   2652,   115)
    reducedCfp9half[ 3,8] = ( 1, 188598, 13915); reducedCfp9half[ 4,8] = (-1,  31824,  2783)
    reducedCfp9half[ 5,8] = ( 1,    204,    23); reducedCfp9half[ 6,8] = (-1,  12996, 13915)
    reducedCfp9half[ 7,8] = ( 0,      1,     1); reducedCfp9half[ 8,8] = ( 1,  25500,  2783)
    reducedCfp9half[ 9,8] = (-1,  38250,  2783); reducedCfp9half[10,8] = (-1,    768,  2783)
    reducedCfp9half[11,8] = ( 1,  81396, 13915); reducedCfp9half[12,8] = ( 1,   3213,   115)
    reducedCfp9half[13,8] = (-1,  60543,  2783); reducedCfp9half[14,8] = (-1,3066144,236555)
    reducedCfp9half[15,8] = ( 1, 727776, 47311); reducedCfp9half[16,8] = ( 1,    207,    17)
    reducedCfp9half[17,8] = (-1,  41553, 21505); reducedCfp9half[18,8] = ( 0,      1,     1)
    reducedCfp9half[19,8] = ( 0,      1,     1); reducedCfp9half[20,8] = ( 0,      1,     1)
    # rcfp for term_jj(54)
    reducedCfp9half[ 1,9] = ( 1,     52,     5); reducedCfp9half[ 2,9] = (-1,     36,     5)
    reducedCfp9half[ 3,9] = ( 1,     84,     5); reducedCfp9half[ 4,9] = (-1,     56,    13)
    reducedCfp9half[ 5,9] = ( 1,    126,    13); reducedCfp9half[ 6,9] = (-1,   3570,   325)
    reducedCfp9half[ 7,9] = ( 0,      1,     1); reducedCfp9half[ 8,9] = ( 1,    360,    13)
    reducedCfp9half[ 9,9] = ( 1,     60,    13); reducedCfp9half[10,9] = ( 1,    102,    13)
    reducedCfp9half[11,9] = (-1,   1064,    65); reducedCfp9half[12,9] = (-1,   2142,    65)
    reducedCfp9half[13,9] = ( 1,     42,    13); reducedCfp9half[14,9] = ( 1,    228,    65)
    reducedCfp9half[15,9] = ( 1,    252,    13); reducedCfp9half[16,9] = ( 1,    342,    13)
    reducedCfp9half[17,9] = (-1,     66,    65); reducedCfp9half[18,9] = ( 1,    230,    13)
    reducedCfp9half[19,9] = ( 0,      1,     1); reducedCfp9half[20,9] = ( 0,      1,     1)
    # rcfp for term_jj(55)
    reducedCfp9half[ 1,10] = ( 0,      1,     1); reducedCfp9half[ 2,10] = ( 1,    144,    11)
    reducedCfp9half[ 3,10] = ( 1,    416,    11); reducedCfp9half[ 4,10] = ( 0,      1,     1)
    reducedCfp9half[ 5,10] = (-1,     32,   165); reducedCfp9half[ 6,10] = ( 0,      1,     1)
    reducedCfp9half[ 7,10] = (-1,    130,     1); reducedCfp9half[ 8,10] = (-1,   1922,    33)
    reducedCfp9half[ 9,10] = ( 0,      1,     1); reducedCfp9half[10,10] = ( 1,    896,    55)
    reducedCfp9half[11,10] = ( 0,      1,     1); reducedCfp9half[12,10] = ( 1,    476,    11)
    reducedCfp9half[13,10] = ( 0,      1,     1); reducedCfp9half[14,10] = ( 1,   1344,    11)
    reducedCfp9half[15,10] = ( 0,      1,     1); reducedCfp9half[16,10] = ( 1,   2052,    55)
    reducedCfp9half[17,10] = ( 0,      1,     1); reducedCfp9half[18,10] = ( 0,      1,     1)
    reducedCfp9half[19,10] = (-1,    308,     5); reducedCfp9half[20,10] = ( 0,      1,     1)
    # rcfp for term_jj(56)
    reducedCfp9half[ 1,11] = ( 0,      1,     1); reducedCfp9half[ 2,11] = ( 1,    132,     5)
    reducedCfp9half[ 3,11] = (-1,  52728,  6655); reducedCfp9half[ 4,11] = ( 1,    196,  1331)
    reducedCfp9half[ 5,11] = (-1,     50,    11); reducedCfp9half[ 6,11] = (-1,  24990,  1331)
    reducedCfp9half[ 7,11] = ( 0,      1,     1); reducedCfp9half[ 8,11] = (-1,  21160,  1331)
    reducedCfp9half[ 9,11] = ( 1,  31740,  1331); reducedCfp9half[10,11] = ( 1,  26250,  1331)
    reducedCfp9half[11,11] = ( 1,  12920,  1331); reducedCfp9half[12,11] = (-1,    357,    55)
    reducedCfp9half[13,11] = ( 1,   9583,  1331); reducedCfp9half[14,11] = (-1, 344988,  6655)
    reducedCfp9half[15,11] = ( 1,   5700,  1331); reducedCfp9half[16,11] = ( 1,    171,    11)
    reducedCfp9half[17,11] = (-1,     15,   121); reducedCfp9half[18,11] = (-1,   4830,   121)
    reducedCfp9half[19,11] = (-1,     84,    11); reducedCfp9half[20,11] = ( 0,      1,     1)
    # rcfp for term_jj(57)
    reducedCfp9half[ 1,12] = ( 0,      1,     1); reducedCfp9half[ 2,12] = ( 0,      1,     1)
    reducedCfp9half[ 3,12] = (-1, 209950,  9317); reducedCfp9half[ 4,12] = ( 1,  77520,  9317)
    reducedCfp9half[ 5,12] = ( 1,  12920,   231); reducedCfp9half[ 6,12] = (-1,  25688,  9317)
    reducedCfp9half[ 7,12] = ( 0,      1,     1); reducedCfp9half[ 8,12] = (-1,   4522,  3993)
    reducedCfp9half[ 9,12] = ( 1,   2261,  1331); reducedCfp9half[10,12] = (-1,  48640, 22627)
    reducedCfp9half[11,12] = (-1, 285144,  9317); reducedCfp9half[12,12] = ( 1,    931,    44)
    reducedCfp9half[13,12] = ( 1, 273885, 90508); reducedCfp9half[14,12] = (-1, 112908,429913)
    reducedCfp9half[15,12] = (-1,2138580,158389); reducedCfp9half[16,12] = ( 1, 137781,  3740)
    reducedCfp9half[17,12] = ( 1,6654375,156332); reducedCfp9half[18,12] = (-1,  59616, 39083)
    reducedCfp9half[19,12] = ( 1, 284089, 17765); reducedCfp9half[20,12] = ( 0,      1,     1)
    # rcfp for term_jj(58)
    reducedCfp9half[ 1,13] = ( 0,      1,     1); reducedCfp9half[ 2,13] = ( 0,      1,     1)
    reducedCfp9half[ 3,13] = ( 1,   1530,    77); reducedCfp9half[ 4,13] = ( 1,  13056,  1001)
    reducedCfp9half[ 5,13] = ( 1,  29376,  1001); reducedCfp9half[ 6,13] = ( 1,    720,  1001)
    reducedCfp9half[ 7,13] = ( 0,      1,     1); reducedCfp9half[ 8,13] = (-1,   1890,   143)
    reducedCfp9half[ 9,13] = (-1,    315,   143); reducedCfp9half[10,13] = ( 1, 101124,  2431)
    reducedCfp9half[11,13] = (-1,   4560,  1001); reducedCfp9half[12,13] = (-1,  13965,   572)
    reducedCfp9half[13,13] = (-1,  35131,  9724); reducedCfp9half[14,13] = ( 1,  13500,  2431)
    reducedCfp9half[15,13] = (-1, 685900, 17017); reducedCfp9half[16,13] = (-1,   1197,   884)
    reducedCfp9half[17,13] = (-1,  28875,   884); reducedCfp9half[18,13] = (-1,   5060,   221)
    reducedCfp9half[19,13] = ( 1,    759,    17); reducedCfp9half[20,13] = ( 0,      1,     1)
    # rcfp for term_jj(59)
    reducedCfp9half[ 1,14] = ( 0,      1,     1); reducedCfp9half[ 2,14] = ( 0,      1,     1)
    reducedCfp9half[ 3,14] = ( 0,      1,     1); reducedCfp9half[ 4,14] = ( 0,      1,     1)
    reducedCfp9half[ 5,14] = ( 1,  22848,   715); reducedCfp9half[ 6,14] = ( 0,      1,     1)
    reducedCfp9half[ 7,14] = (-1,    170,     1); reducedCfp9half[ 8,14] = ( 1,    918,   143)
    reducedCfp9half[ 9,14] = ( 0,      1,     1); reducedCfp9half[10,14] = (-1,  32832,   715)
    reducedCfp9half[11,14] = ( 0,      1,     1); reducedCfp9half[12,14] = ( 1,   9044,   143)
    reducedCfp9half[13,14] = ( 0,      1,     1); reducedCfp9half[14,14] = (-1,    576,    13)
    reducedCfp9half[15,14] = ( 0,      1,     1); reducedCfp9half[16,14] = ( 1,   7524,    65)
    reducedCfp9half[17,14] = ( 0,      1,     1); reducedCfp9half[18,14] = ( 0,      1,     1)
    reducedCfp9half[19,14] = ( 1,   1012,     5); reducedCfp9half[20,14] = ( 0,      1,     1)
    # rcfp for term_jj(60)
    reducedCfp9half[ 1,15] = ( 0,      1,     1); reducedCfp9half[ 2,15] = ( 0,      1,     1)
    reducedCfp9half[ 3,15] = ( 0,      1,     1); reducedCfp9half[ 4,15] = ( 0,      1,     1)
    reducedCfp9half[ 5,15] = ( 0,      1,     1); reducedCfp9half[ 6,15] = ( 1,   2128,   143)
    reducedCfp9half[ 7,15] = ( 0,      1,     1); reducedCfp9half[ 8,15] = (-1,   1938,   143)
    reducedCfp9half[ 9,15] = ( 1,   2907,   143); reducedCfp9half[10,15] = (-1,   4860,   143)
    reducedCfp9half[11,15] = (-1,  17136,  2717); reducedCfp9half[12,15] = (-1,   1309,    52)
    reducedCfp9half[13,15] = (-1,   8505,   572); reducedCfp9half[14,15] = ( 1,  15876,   247)
    reducedCfp9half[15,15] = ( 1,    420,    13); reducedCfp9half[16,15] = ( 1,   1287,    20)
    reducedCfp9half[17,15] = (-1,   6075,   988); reducedCfp9half[18,15] = (-1,    132,    19)
    reducedCfp9half[19,15] = (-1,    253,    95); reducedCfp9half[20,15] = (-1,    650,    19)
    # rcfp for term_jj(61)
    reducedCfp9half[ 1,16] = ( 0,      1,     1); reducedCfp9half[ 2,16] = ( 0,      1,     1)
    reducedCfp9half[ 3,16] = ( 0,      1,     1); reducedCfp9half[ 4,16] = ( 0,      1,     1)
    reducedCfp9half[ 5,16] = ( 0,      1,     1); reducedCfp9half[ 6,16] = ( 0,      1,     1)
    reducedCfp9half[ 7,16] = ( 0,      1,     1); reducedCfp9half[ 8,16] = ( 1,    570,    13)
    reducedCfp9half[ 9,16] = ( 1,     95,    13); reducedCfp9half[10,16] = ( 1,   4104,   221)
    reducedCfp9half[11,16] = ( 1,    504,    65); reducedCfp9half[12,16] = (-1,   1463,   260)
    reducedCfp9half[13,16] = (-1,  39501,   884); reducedCfp9half[14,16] = (-1,  60516,  1105)
    reducedCfp9half[15,16] = (-1,   1596,   221); reducedCfp9half[16,16] = ( 1,   3933,    68)
    reducedCfp9half[17,16] = ( 1,    621,  3740); reducedCfp9half[18,16] = (-1,    840,    17)
    reducedCfp9half[19,16] = (-1,    805,    17); reducedCfp9half[20,16] = ( 1,    390,    11)
    # rcfp for term_jj(62)
    reducedCfp9half[ 1,17] = ( 0,      1,     1); reducedCfp9half[ 2,17] = ( 0,      1,     1)
    reducedCfp9half[ 3,17] = ( 0,      1,     1); reducedCfp9half[ 4,17] = ( 0,      1,     1)
    reducedCfp9half[ 5,17] = ( 0,      1,     1); reducedCfp9half[ 6,17] = ( 0,      1,     1)
    reducedCfp9half[ 7,17] = ( 0,      1,     1); reducedCfp9half[ 8,17] = ( 0,      1,     1)
    reducedCfp9half[ 9,17] = ( 0,      1,     1); reducedCfp9half[10,17] = (-1,   5796,   221)
    reducedCfp9half[11,17] = ( 1,  17664,  1235); reducedCfp9half[12,17] = (-1,   5313,    65)
    reducedCfp9half[13,17] = ( 1,   1771,   221); reducedCfp9half[14,17] = (-1,  16632,  1615)
    reducedCfp9half[15,17] = (-1,     88,    17); reducedCfp9half[16,17] = ( 1,    693,    17)
    reducedCfp9half[17,17] = (-1,  94269,  1615); reducedCfp9half[18,17] = ( 1, 192500,  7429)
    reducedCfp9half[19,17] = ( 1,  30030,   323); reducedCfp9half[20,17] = ( 1,  24570,   437)
    # rcfp for term_jj(63)
    reducedCfp9half[ 1,18] = ( 0,      1,     1); reducedCfp9half[ 2,18] = ( 0,      1,     1)
    reducedCfp9half[ 3,18] = ( 0,      1,     1); reducedCfp9half[ 4,18] = ( 0,      1,     1)
    reducedCfp9half[ 5,18] = ( 0,      1,     1); reducedCfp9half[ 6,18] = ( 0,      1,     1)
    reducedCfp9half[ 7,18] = ( 0,      1,     1); reducedCfp9half[ 8,18] = ( 0,      1,     1)
    reducedCfp9half[ 9,18] = ( 0,      1,     1); reducedCfp9half[10,18] = ( 0,      1,     1)
    reducedCfp9half[11,18] = ( 0,      1,     1); reducedCfp9half[12,18] = ( 0,      1,     1)
    reducedCfp9half[13,18] = ( 0,      1,     1); reducedCfp9half[14,18] = (-1,  15000,   323)
    reducedCfp9half[15,18] = ( 1,    280,    17); reducedCfp9half[16,18] = (-1,   1170,    17)
    reducedCfp9half[17,18] = (-1,  48750,  3553); reducedCfp9half[18,18] = (-1,  15600,   437)
    reducedCfp9half[19,18] = ( 1,   3510,    19); reducedCfp9half[20,18] = (-1,  33930,   253) 


"""
`SpinAngular.W_00_one_half[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(00) :: j QJ)  for j = 1/2
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_00_one_half = Array{Tuple{Int64,Int64,Int64}}(undef,2);       for i=1:2  W_00_one_half[i] = (0,0,0)   end


"""
`SpinAngular.W_00_three_half[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(00) :: j QJ)  for j = 3/2
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_00_three_half = Array{Tuple{Int64,Int64,Int64}}(undef,3);       for i=1:3  W_00_three_half[i] = (0,0,0)   end


"""
`SpinAngular.W_00_five_half[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(00) :: j QJ)  for j = 5/2
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_00_five_half = Array{Tuple{Int64,Int64,Int64}}(undef,6);       for i=1:6  W_00_five_half[i] = (0,0,0)   end


"""
`SpinAngular.W_00_seven_half[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(00) :: j QJ)  for j = 7/2
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_00_seven_half = Array{Tuple{Int64,Int64,Int64}}(undef,14);       for i=1:14  W_00_seven_half[i] = (0,0,0)   end

    W_00_one_half[1]    = (-1,     2,     1); W_00_one_half[2]    = (-1,       2,     1)

    W_00_three_half[1]  = (-1,    16,     1); W_00_three_half[2]  = (-1,       6,     1)
    W_00_three_half[3]  = (-1,    10,     1)

    W_00_five_half[1]   = (-1,    54,     1); W_00_five_half[2]   = (-1,      12,     1)
    W_00_five_half[3]   = (-1,    30,     1); W_00_five_half[4]   = (-1,      12,     1)
    W_00_five_half[5]   = (-1,    30,     1); W_00_five_half[6]   = (-1,      54,     1)

    W_00_seven_half[ 1] = (-1,    32,     1); W_00_seven_half[ 2] = (-1,      48,     1)
    W_00_seven_half[ 3] = (-1,   128,     1); W_00_seven_half[ 4] = (-1,      80,     1)
    W_00_seven_half[ 5] = (-1,    96,     1); W_00_seven_half[ 6] = (-1,     128,     1)
    W_00_seven_half[ 7] = (-1,    20,     1); W_00_seven_half[ 8] = (-1,      60,     1)
    W_00_seven_half[ 9] = (-1,    20,     1); W_00_seven_half[10] = (-1,     108,     1)
    W_00_seven_half[11] = (-1,    36,     1); W_00_seven_half[12] = (-1,      44,     1)
    W_00_seven_half[13] = (-1,   156,     1); W_00_seven_half[14] = (-1,      68,     1)


"""
`SpinAngular.W_10_one_half[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(10) :: j QJ)  for j = 1/2
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_10_one_half = Array{Tuple{Int64,Int64,Int64}}(undef,2);       for i=1:2  W_10_one_half[i] = (0,0,0)   end


"""
`SpinAngular.W_10_three_half[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(10) :: j QJ)  for j = 3/2
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_10_three_half = Array{Tuple{Int64,Int64,Int64}}(undef,3);       for i=1:3  W_10_three_half[i] = (0,0,0)   end


"""
`SpinAngular.W_10_five_half[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(10) :: j QJ)  for j = 5/2
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_10_five_half = Array{Tuple{Int64,Int64,Int64}}(undef,6);       for i=1:6  W_10_five_half[i] = (0,0,0)   end


"""
`SpinAngular.W_10_seven_half[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(10) :: j QJ)  for j = 7/2
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_10_seven_half = Array{Tuple{Int64,Int64,Int64}}(undef,14);       for i=1:14  W_10_seven_half[i] = (0,0,0)   end

    W_10_one_half[1] = ( 0,     1,     1);    W_10_one_half[2] = (-1,       6,     1)

    W_10_three_half[1] = (-1,    12,     1);  W_10_three_half[2] = (-1,      12,     1)
    W_10_three_half[3] = ( 0,     1,     1)

    W_10_five_half[1] = (-1,    48,     1);   W_10_five_half[2] = ( 0,       1,     1)
    W_10_five_half[3] = ( 0,     1,     1);   W_10_five_half[4] = (-1,      20,     1)
    W_10_five_half[5] = (-1,    10,     1);   W_10_five_half[6] = (-1,      18,     1)

    W_10_seven_half[ 1] = (-1,     6,     1); W_10_seven_half[ 2] = (-1,       9,     1)
    W_10_seven_half[ 3] = (-1,   120,     1); W_10_seven_half[ 4] = (-1,      15,     1)
    W_10_seven_half[ 5] = (-1,    18,     1); W_10_seven_half[ 6] = (-1,      24,     1)
    W_10_seven_half[ 7] = (-1,    30,     1); W_10_seven_half[ 8] = (-1,      30,     1)
    W_10_seven_half[ 9] = ( 0,     1,     1); W_10_seven_half[10] = (-1,      54,     1)
    W_10_seven_half[11] = ( 0,     1,     1); W_10_seven_half[12] = ( 0,       1,     1)
    W_10_seven_half[13] = (-1,    78,     1); W_10_seven_half[14] = ( 0,       1,     1)


"""
`SpinAngular.W_01_one_half[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(01) :: j QJ)  for j = 1/2
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_01_one_half = Array{Tuple{Int64,Int64,Int64}}(undef,2);       for i=1:2  W_01_one_half[i] = (0,0,0)   end


"""
`SpinAngular.W_01_three_half[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(01) :: j QJ)  for j = 3/2
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_01_three_half = Array{Tuple{Int64,Int64,Int64}}(undef,3);       for i=1:3  W_01_three_half[i] = (0,0,0)   end


"""
`SpinAngular.W_01_five_half[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(01) :: j QJ)  for j = 5/2
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_01_five_half = Array{Tuple{Int64,Int64,Int64}}(undef,6);       for i=1:6  W_01_five_half[i] = (0,0,0)   end


"""
`SpinAngular.W_01_seven_half[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(01) :: j QJ)  for j = 7/2
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_01_seven_half = Array{Tuple{Int64,Int64,Int64}}(undef,14);       for i=1:14  W_01_seven_half[i] = (0,0,0)   end

    W_01_one_half[1]    = (-1,     6,     1); W_01_one_half[2]    = ( 0,       1,     1)

    W_01_three_half[1]  = (-1,    12,     1); W_01_three_half[2]  = ( 0,       1,     1)
    W_01_three_half[3]  = (-1,    12,     1)

    W_01_five_half[1]   = (-1,   126,     7); W_01_five_half[2]   = (-1,      12,     7)
    W_01_five_half[3]   = (-1,   198,     7); W_01_five_half[4]   = ( 0,       1,     1)
    W_01_five_half[5]   = (-1,    48,     7); W_01_five_half[6]   = (-1,     288,     7)

    W_01_seven_half[ 1] = (-1,    10,     7); W_01_seven_half[ 2] = (-1,       5,     1)
    W_01_seven_half[ 3] = (-1,   168,     7); W_01_seven_half[ 4] = (-1,     165,     7)
    W_01_seven_half[ 5] = (-1,   286,     7); W_01_seven_half[ 6] = (-1,     680,     7)
    W_01_seven_half[ 7] = ( 0,     1,     1); W_01_seven_half[ 8] = (-1,      30,     7)
    W_01_seven_half[ 9] = (-1,    10,     7); W_01_seven_half[10] = (-1,     180,     7)
    W_01_seven_half[11] = (-1,    60,     7); W_01_seven_half[12] = (-1,     110,     7)
    W_01_seven_half[13] = (-1,   546,     7); W_01_seven_half[14] = (-1,     408,     7)


"""
`SpinAngular.W_12_three_half_odd[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(12) :: j QJ) j = 3/2  for odd seniority 
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_12_three_half_odd = Matrix{Tuple{Int64,Int64,Int64}}(undef,1,1);  for i=1:1  for j=1:1  W_12_three_half_odd[i,j] = (0,0,0)   end  end


"""
`SpinAngular.W_12_three_half_even[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(12) :: j QJ) j = 3/2  for even seniority 
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_12_three_half_even = Matrix{Tuple{Int64,Int64,Int64}}(undef,2,2); for i=1:2  for j=1:2  W_12_three_half_even[i,j] = (0,0,0)   end  end


"""
`SpinAngular.W_12_five_half_odd[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(12) :: j QJ)  j = 5/2  for odd seniority
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_12_five_half_odd = Matrix{Tuple{Int64,Int64,Int64}}(undef,3,3);   for i=1:3  for j=1:3  W_12_five_half_odd[i,j] = (0,0,0)   end  end


"""
`SpinAngular.W_12_five_half_even[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(12) :: j QJ)  j = 5/2  for even seniority
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_12_five_half_even = Matrix{Tuple{Int64,Int64,Int64}}(undef,3,3);  for i=1:3  for j=1:3  W_12_five_half_even[i,j] = (0,0,0)   end  end


"""
`SpinAngular.W_12_seven_half_odd[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(12) :: j QJ)  j = 7/2  for odd seniority
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_12_seven_half_odd = Array{Tuple{Int64,Int64,Int64}}(undef,21);    for i=1:21  W_12_seven_half_odd[i] = (0,0,0)   end


"""
`SpinAngular.W_12_seven_half_even[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(12) :: j QJ)  j = 7/2  for even seniority
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_12_seven_half_even = Array{Tuple{Int64,Int64,Int64}}(undef,36);       for i=1:36    W_12_seven_half_even[i] = (0,0,0)   end 

    W_12_three_half_odd[1,1]   = ( 1,    60,     1)

    W_12_three_half_even[1, 1] = ( 0,     1,     1); W_12_three_half_even[2, 1] = ( 1,      30,     1)
    W_12_three_half_even[1, 2] = (-1,    30,     1); W_12_three_half_even[2, 2] = ( 0,       1,     1)

    W_12_five_half_odd[1, 1]   = ( 1,   420,     7); W_12_five_half_odd[2, 1]   = ( 1,     360,     7)
    W_12_five_half_odd[3, 1]   = ( 1,   270,     7); W_12_five_half_odd[1, 2]   = ( 1,     360,     7)
    W_12_five_half_odd[2, 2]   = ( 0,     1,     1); W_12_five_half_odd[3, 2]   = ( 0,       1,     1)
    W_12_five_half_odd[1, 3]   = (-1,   270,     7); W_12_five_half_odd[2, 3]   = ( 0,       1,     1)
    W_12_five_half_odd[3, 3]   = ( 0,     1,     1)      

    W_12_five_half_even[1, 1]  = ( 0,     1,     1); W_12_five_half_even[2, 1]  = ( 1,    1960,    49)
    W_12_five_half_even[3, 1]  = ( 0,     1,     1); W_12_five_half_even[1, 2]  = (-1,    1960,    49)
    W_12_five_half_even[2, 2]  = (-1,  1000,    49); W_12_five_half_even[3, 2]  = ( 1,    2430,    49)
    W_12_five_half_even[1, 3]  = ( 0,     1,     1); W_12_five_half_even[2, 3]  = ( 1,    2430,    49)
    W_12_five_half_even[3, 3]  = ( 1,  1980,    49)      

    W_12_seven_half_odd[ 1]    = (-1,   252,    10); W_12_seven_half_odd[ 2]    = ( 1,    1056,    70)
    W_12_seven_half_odd[ 3]    = ( 1,   144,     7); W_12_seven_half_odd[ 4]    = ( 0,       1,     1)
    W_12_seven_half_odd[ 5]    = ( 0,     1,     1); W_12_seven_half_odd[ 6]    = ( 0,       1,     1)
    W_12_seven_half_odd[ 7]    = ( 1,   507,    10); W_12_seven_half_odd[ 8]    = (-1,      88,     1)
    W_12_seven_half_odd[ 9]    = (-1,   325,    14); W_12_seven_half_odd[10]    = ( 0,       1,     1)
    W_12_seven_half_odd[11]    = ( 0,     1,     1); W_12_seven_half_odd[12]    = ( 1,     200,     3)
    W_12_seven_half_odd[13]    = (-1,   520,    21); W_12_seven_half_odd[14]    = ( 1,      80,     1)
    W_12_seven_half_odd[15]    = ( 0,     1,     1); W_12_seven_half_odd[16]    = (-1,    6125,   462)
    W_12_seven_half_odd[17]    = (-1,   560,    11); W_12_seven_half_odd[18]    = ( 0,       1,     1)
    W_12_seven_half_odd[19]    = ( 1,   390,   539); W_12_seven_half_odd[20]    = (-1,    3840,    49)
    W_12_seven_half_odd[21]    = ( 1,  2040,    49)      

    W_12_seven_half_even[ 1]   = ( 0,     1,     1); W_12_seven_half_even[ 2]   = (-1,      50,     1)
    W_12_seven_half_even[ 3]   = ( 0,     1,     1); W_12_seven_half_even[ 4]   = ( 0,       1,     1)
    W_12_seven_half_even[ 5]   = ( 0,     1,     1); W_12_seven_half_even[ 6]   = ( 0,       1,     1)
    W_12_seven_half_even[ 7]   = ( 0,     1,     1); W_12_seven_half_even[ 8]   = ( 0,       1,     1)
    W_12_seven_half_even[ 9]   = (-1,  1280,    49); W_12_seven_half_even[10]   = ( 1,     990,    49)
    W_12_seven_half_even[11]   = ( 1,  2640,    49); W_12_seven_half_even[12]   = ( 1,    1950,    49)
    W_12_seven_half_even[13]   = ( 0,     1,     1); W_12_seven_half_even[14]   = ( 0,       1,     1)
    W_12_seven_half_even[15]   = ( 0,     1,     1); W_12_seven_half_even[16]   = ( 0,       1,     1)
    W_12_seven_half_even[17]   = (-1,   480,    49); W_12_seven_half_even[18]   = ( 0,       1,     1)
    W_12_seven_half_even[19]   = ( 0,     1,     1); W_12_seven_half_even[20]   = ( 0,       1,     1)
    W_12_seven_half_even[21]   = ( 0,     1,     1); W_12_seven_half_even[22]   = (-1,     360,   539)
    W_12_seven_half_even[23]   = ( 1,  1872,    49); W_12_seven_half_even[24]   = ( 1,      42,     1)
    W_12_seven_half_even[25]   = ( 1,   390,    11); W_12_seven_half_even[26]   = ( 0,       1,     1)
    W_12_seven_half_even[27]   = ( 0,     1,     1); W_12_seven_half_even[28]   = ( 0,       1,     1)
    W_12_seven_half_even[29]   = ( 1,    48,     1); W_12_seven_half_even[30]   = ( 0,       1,     1)
    W_12_seven_half_even[31]   = ( 0,     1,     1); W_12_seven_half_even[32]   = ( 1,     234,     7)
    W_12_seven_half_even[33]   = ( 0,     1,     1); W_12_seven_half_even[34]   = ( 1,    1040,    11)
    W_12_seven_half_even[35]   = ( 1,   340,     7); W_12_seven_half_even[36]   = ( 0,       1,     1)
    

"""
`SpinAngular.W_03_three_half_odd[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(03) :: j QJ) j = 3/2  for odd seniority 
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_03_three_half_odd = Matrix{Tuple{Int64,Int64,Int64}}(undef,1,1);  for i=1:1  for j=1:1  W_03_three_half_odd[i,j] = (0,0,0)   end  end


"""
`SpinAngular.W_03_three_half_even[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(03) :: j QJ) j = 3/2  for even seniority 
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_03_three_half_even = Matrix{Tuple{Int64,Int64,Int64}}(undef,2,2); for i=1:2  for j=1:2  W_03_three_half_even[i,j] = (0,0,0)   end  end


"""
`SpinAngular.W_03_five_half_odd[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(03) :: j QJ)  j = 5/2  for odd seniority
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_03_five_half_odd = Matrix{Tuple{Int64,Int64,Int64}}(undef,3,3);   for i=1:3  for j=1:3  W_03_five_half_odd[i,j] = (0,0,0)   end  end


"""
`SpinAngular.W_03_five_half_even[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(03) :: j QJ)  j = 5/2  for even seniority
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_03_five_half_even = Matrix{Tuple{Int64,Int64,Int64}}(undef,3,3);  for i=1:3  for j=1:3  W_03_five_half_even[i,j] = (0,0,0)   end  end


"""
`SpinAngular.W_03_seven_half_odd[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(03) :: j QJ)  j = 7/2  for odd seniority
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_03_seven_half_odd = Array{Tuple{Int64,Int64,Int64}}(undef,21);    for i=1:21  W_03_seven_half_odd[i] = (0,0,0)   end


"""
`SpinAngular.W_03_seven_half_even[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(03) :: j QJ)  j = 7/2  for even seniority
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_03_seven_half_even = Array{Tuple{Int64,Int64,Int64}}(undef,36);   for i=1:36    W_03_seven_half_even[i] = (0,0,0)   end 

    W_03_three_half_odd[1, 1]  = (-1,    28,     1)

    W_03_three_half_even[1, 1] = ( 0,     1,     1); W_03_three_half_even[2, 1] = ( 0,       1,     1)
    W_03_three_half_even[1, 2] = ( 0,     1,     1); W_03_three_half_even[2, 2] = ( 1,      28,     1)

    W_03_five_half_odd[1, 1]   = (-1,   882,    21); W_03_five_half_odd[2, 1]   = ( 0,       1,     1)
    W_03_five_half_odd[3, 1]   = ( 0,     1,     1); W_03_five_half_odd[1, 2]   = ( 0,       1,     1)
    W_03_five_half_odd[2, 2]   = ( 1,   384,    21); W_03_five_half_odd[3, 2]   = (-1,     400,    21)
    W_03_five_half_odd[1, 3]   = ( 0,     1,     1); W_03_five_half_odd[2, 3]   = ( 1,     400,    21)
    W_03_five_half_odd[3, 3]   = ( 1,   286,    21)

    W_03_five_half_even[1, 1]  = ( 0,     1,     1); W_03_five_half_even[2, 1]  = ( 0,       1,     1)
    W_03_five_half_even[3, 1]  = ( 0,     1,     1); W_03_five_half_even[1, 2]  = ( 0,       1,     1)
    W_03_five_half_even[2, 2]  = ( 1,   162,     7); W_03_five_half_even[3, 2]  = (-1,     300,     7)
    W_03_five_half_even[1, 3]  = ( 0,     1,     1); W_03_five_half_even[2, 3]  = (-1,     300,     7)
    W_03_five_half_even[3, 3]  = ( 1,    22,     7)

    W_03_seven_half_odd[ 1]    = (-1,  1188,    70); W_03_seven_half_odd[ 2]    = (-1,     196,    10)
    W_03_seven_half_odd[ 3]    = ( 0,     1,     1); W_03_seven_half_odd[ 4]    = (-1,     234,     7)
    W_03_seven_half_odd[ 5]    = ( 0,     1,     1); W_03_seven_half_odd[ 6]    = ( 0,       1,     1)
    W_03_seven_half_odd[ 7]    = ( 1,   189,   110); W_03_seven_half_odd[ 8]    = ( 0,       1,     1)
    W_03_seven_half_odd[ 9]    = (-1,  1911,   242); W_03_seven_half_odd[10]    = ( 1,    1470,   121)
    W_03_seven_half_odd[11]    = ( 0,     1,     1); W_03_seven_half_odd[12]    = (-1,      56,     1)
    W_03_seven_half_odd[13]    = ( 0,     1,     1); W_03_seven_half_odd[14]    = ( 0,       1,     1)
    W_03_seven_half_odd[15]    = ( 0,     1,     1); W_03_seven_half_odd[16]    = ( 1,  394805, 22022)
    W_03_seven_half_odd[17]    = ( 1,  5250,   121); W_03_seven_half_odd[18]    = ( 1,   53760,  1573)
    W_03_seven_half_odd[19]    = (-1,    78,   847); W_03_seven_half_odd[20]    = ( 1,   17408,   847)
    W_03_seven_half_odd[21]    = ( 1, 12920,  1001)

    W_03_seven_half_even[ 1]   = ( 0,     1,     1); W_03_seven_half_even[ 2]   = ( 0,       1,     1)
    W_03_seven_half_even[ 3]   = ( 0,     1,     1); W_03_seven_half_even[ 4]   = ( 0,       1,     1)
    W_03_seven_half_even[ 5]   = ( 0,     1,     1); W_03_seven_half_even[ 6]   = ( 0,       1,     1)
    W_03_seven_half_even[ 7]   = ( 0,     1,     1); W_03_seven_half_even[ 8]   = ( 0,       1,     1)
    W_03_seven_half_even[ 9]   = ( 1,   110,     7); W_03_seven_half_even[10]   = ( 0,       1,     1)
    W_03_seven_half_even[11]   = (-1,   240,     7); W_03_seven_half_even[12]   = ( 0,       1,     1)
    W_03_seven_half_even[13]   = ( 0,     1,     1); W_03_seven_half_even[14]   = ( 0,       1,     1)
    W_03_seven_half_even[15]   = ( 0,     1,     1); W_03_seven_half_even[16]   = (-1,    1920,    77)
    W_03_seven_half_even[17]   = ( 0,     1,     1); W_03_seven_half_even[18]   = ( 1,      52,     7)
    W_03_seven_half_even[19]   = (-1,   224,    11); W_03_seven_half_even[20]   = ( 0,       1,     1)
    W_03_seven_half_even[21]   = ( 0,     1,     1); W_03_seven_half_even[22]   = ( 1,   32490,   847)
    W_03_seven_half_even[23]   = ( 0,     1,     1); W_03_seven_half_even[24]   = ( 0,       1,     1)
    W_03_seven_half_even[25]   = (-1,  7644,   121); W_03_seven_half_even[26]   = ( 0,       1,     1)
    W_03_seven_half_even[27]   = ( 1,    12,    70); W_03_seven_half_even[28]   = ( 1,     224,    10)
    W_03_seven_half_even[29]   = ( 0,     1,     1); W_03_seven_half_even[30]   = ( 0,       1,     1)
    W_03_seven_half_even[31]   = (-1,    52,    70); W_03_seven_half_even[32]   = ( 0,       1,     1)
    W_03_seven_half_even[33]   = (-1,  2040,    77); W_03_seven_half_even[34]   = (-1,     364,   121)
    W_03_seven_half_even[35]   = ( 0,     1,     1); W_03_seven_half_even[36]   = ( 1,    1292,    77)


"""
`SpinAngular.W_14_five_half_odd[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(14) :: j QJ)  j = 5/2  for odd seniority
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_14_five_half_odd = Matrix{Tuple{Int64,Int64,Int64}}(undef,3,3);   for i=1:3  for j=1:3  W_14_five_half_odd[i,j] = (0,0,0)   end  end


"""
`SpinAngular.W_14_five_half_even[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(14) :: j QJ)  j = 5/2  for even seniority
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_14_five_half_even = Matrix{Tuple{Int64,Int64,Int64}}(undef,3,3);  for i=1:3  for j=1:3  W_14_five_half_even[i,j] = (0,0,0)   end  end


"""
`SpinAngular.W_14_seven_half_odd[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(14) :: j QJ)  j = 7/2  for odd seniority
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_14_seven_half_odd = Array{Tuple{Int64,Int64,Int64}}(undef,21);    for i=1:21  W_14_seven_half_odd[i] = (0,0,0)   end


"""
`SpinAngular.W_14_seven_half_even[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(14) :: j QJ)  j = 7/2  for even seniority
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_14_seven_half_even = Array{Tuple{Int64,Int64,Int64}}(undef,36);       for i=1:36    W_14_seven_half_even[i] = (0,0,0)   end 

    W_14_five_half_odd[1, 1] = ( 1,  3780,    35); W_14_five_half_odd[2, 1] = (-1,     720,    35)
    W_14_five_half_odd[3, 1] = (-1,  4950,    35); W_14_five_half_odd[1, 2] = (-1,     720,    35)
    W_14_five_half_odd[2, 2] = ( 0,     1,     1); W_14_five_half_odd[3, 2] = ( 0,       1,     1)
    W_14_five_half_odd[1, 3] = ( 1,  4950,    35); W_14_five_half_odd[2, 3] = ( 0,       1,     1)
    W_14_five_half_odd[3, 3] = ( 0,     1,     1)      

    W_14_five_half_even[1, 1] = ( 0,     1,     1); W_14_five_half_even[2, 1] = ( 0,       1,     1)
    W_14_five_half_even[3, 1] = ( 1,  3528,    49); W_14_five_half_even[1, 2] = ( 0,       1,     1)
    W_14_five_half_even[2, 2] = ( 1,  2430,    49); W_14_five_half_even[3, 2] = ( 1,    1980,    49)
    W_14_five_half_even[1, 3] = (-1,  3528,    49); W_14_five_half_even[2, 3] = ( 1,    1980,    49)
    W_14_five_half_even[3, 3] = (-1,  7722,    49)

    W_14_seven_half_odd[ 1] = ( 0,     1,     1); W_14_seven_half_odd[ 2] = ( 1,      54,     7)
    W_14_seven_half_odd[ 3] = (-1,   528,     7); W_14_seven_half_odd[ 4] = ( 1,     546,    11)
    W_14_seven_half_odd[ 5] = (-1,   378,    11); W_14_seven_half_odd[ 6] = ( 0,       1,     1)
    W_14_seven_half_odd[ 7] = (-1,  4335,   110); W_14_seven_half_odd[ 8] = (-1,      96,    11)
    W_14_seven_half_odd[ 9] = ( 1, 11271,  1694); W_14_seven_half_odd[10] = ( 1,    3822,   121)
    W_14_seven_half_odd[11] = ( 0,     1,     1); W_14_seven_half_odd[12] = ( 1,     120,     1)
    W_14_seven_half_odd[13] = ( 1, 12000,    77); W_14_seven_half_odd[14] = (-1,     624,    11)
    W_14_seven_half_odd[15] = (-1,   960,    11); W_14_seven_half_odd[16] = (-1,   30345,  3146)
    W_14_seven_half_odd[17] = (-1,   210,   121); W_14_seven_half_odd[18] = (-1,  228480,  1573)
    W_14_seven_half_odd[19] = ( 1,580476,  5929); W_14_seven_half_odd[20] = ( 1,  146880,  5929)
    W_14_seven_half_odd[21] = (-1,627912,  7007)

    W_14_seven_half_even[ 1] = ( 0,     1,     1); W_14_seven_half_even[ 2] = ( 0,       1,     1)
    W_14_seven_half_even[ 3] = ( 0,     1,     1); W_14_seven_half_even[ 4] = (-1,      90,     1)
    W_14_seven_half_even[ 5] = ( 0,     1,     1); W_14_seven_half_even[ 6] = ( 0,       1,     1)
    W_14_seven_half_even[ 7] = ( 0,     1,     1); W_14_seven_half_even[ 8] = ( 0,       1,     1)
    W_14_seven_half_even[ 9] = ( 1,  2640,    49); W_14_seven_half_even[10] = ( 1,     480,    49)
    W_14_seven_half_even[11] = (-1,   360,   539); W_14_seven_half_even[12] = ( 1,   20592,   539)
    W_14_seven_half_even[13] = (-1,    42,     1); W_14_seven_half_even[14] = ( 1,     390,    11)
    W_14_seven_half_even[15] = ( 0,     1,     1); W_14_seven_half_even[16] = ( 0,       1,     1)
    W_14_seven_half_even[17] = ( 1,468180,  5929); W_14_seven_half_even[18] = ( 0,       1,     1)
    W_14_seven_half_even[19] = ( 0,     1,     1); W_14_seven_half_even[20] = ( 1,   21840,   847)
    W_14_seven_half_even[21] = ( 0,     1,     1); W_14_seven_half_even[22] = (-1,  359424,  5929)
    W_14_seven_half_even[23] = (-1,  6750,   539); W_14_seven_half_even[24] = ( 0,       1,     1)
    W_14_seven_half_even[25] = ( 1,917280,  5929); W_14_seven_half_even[26] = (-1,   10710,   121)
    W_14_seven_half_even[27] = ( 0,     1,     1); W_14_seven_half_even[28] = ( 0,       1,     1)
    W_14_seven_half_even[29] = ( 1,    36,    11); W_14_seven_half_even[30] = ( 0,       1,     1)
    W_14_seven_half_even[31] = ( 0,     1,     1); W_14_seven_half_even[32] = (-1,     858,     7)
    W_14_seven_half_even[33] = ( 0,     1,     1); W_14_seven_half_even[34] = (-1,    5304,   121)
    W_14_seven_half_even[35] = (-1, 69768,   847); W_14_seven_half_even[36] = ( 0,       1,     1)


"""
`SpinAngular.W_05_five_half_odd[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(05) :: j QJ)  j = 5/2  for odd seniority
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_05_five_half_odd = Matrix{Tuple{Int64,Int64,Int64}}(undef,3,3);   for i=1:3  for j=1:3  W_05_five_half_odd[i,j] = (0,0,0)   end  end


"""
`SpinAngular.W_05_five_half_even[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(05) :: j QJ)  j = 5/2  for even seniority
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_05_five_half_even = Matrix{Tuple{Int64,Int64,Int64}}(undef,3,3);  for i=1:3  for j=1:3  W_05_five_half_even[i,j] = (0,0,0)   end  end


"""
`SpinAngular.W_05_seven_half_odd[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(05) :: j QJ)  j = 7/2  for odd seniority
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_05_seven_half_odd = Array{Tuple{Int64,Int64,Int64}}(undef,21);    for i=1:21  W_05_seven_half_odd[i] = (0,0,0)   end


"""
`SpinAngular.W_05_seven_half_even[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(05) :: j QJ)  j = 7/2  for even seniority
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_05_seven_half_even = Array{Tuple{Int64,Int64,Int64}}(undef,36);   for i=1:36    W_05_seven_half_even[i] = (0,0,0)   end 

    W_05_five_half_odd[1, 1] = (-1,  1386,    21); W_05_five_half_odd[2, 1] = ( 0,       1,     1)
    W_05_five_half_odd[3, 1] = ( 0,     1,     1); W_05_five_half_odd[1, 2] = ( 0,       1,     1)
    W_05_five_half_odd[2, 2] = ( 0,     1,     1); W_05_five_half_odd[3, 2] = (-1,     440,    21)
    W_05_five_half_odd[1, 3] = ( 0,     1,     1); W_05_five_half_odd[2, 3] = ( 1,     440,    21)
    W_05_five_half_odd[3, 3] = (-1,  1430,    21)

    W_05_five_half_even[1, 1] = ( 0,     1,     1); W_05_five_half_even[2, 1] = ( 0,       1,     1)
    W_05_five_half_even[3, 1] = ( 0,     1,     1); W_05_five_half_even[1, 2] = ( 0,       1,     1)
    W_05_five_half_even[2, 2] = ( 0,     1,     1); W_05_five_half_even[3, 2] = ( 1,     330,     7)
    W_05_five_half_even[1, 3] = ( 0,     1,     1); W_05_five_half_even[2, 3] = ( 1,     330,     7)
    W_05_five_half_even[3, 3] = ( 1,   572,     7)

    W_05_seven_half_odd[ 1] = ( 0,     1,     1); W_05_seven_half_odd[ 2] = ( 0,       1,     1)
    W_05_seven_half_odd[ 3] = ( 0,     1,     1); W_05_seven_half_odd[ 4] = ( 1,     144,     7)
    W_05_seven_half_odd[ 5] = ( 1,    14,     1); W_05_seven_half_odd[ 6] = ( 0,       1,     1)
    W_05_seven_half_odd[ 7] = ( 1,   975,    22); W_05_seven_half_odd[ 8] = ( 0,       1,     1)
    W_05_seven_half_odd[ 9] = (-1, 50421,  2002); W_05_seven_half_odd[10] = ( 1,     336,    11)
    W_05_seven_half_odd[11] = (-1,  7000,   143); W_05_seven_half_odd[12] = (-1,      88,     1)
    W_05_seven_half_odd[13] = ( 0,     1,     1); W_05_seven_half_odd[14] = ( 0,       1,     1)
    W_05_seven_half_odd[15] = ( 0,     1,     1); W_05_seven_half_odd[16] = (-1,    6845, 26026)
    W_05_seven_half_odd[17] = (-1,  3360,   143); W_05_seven_half_odd[18] = ( 1,   14280,  1859)
    W_05_seven_half_odd[19] = (-1,  1836,    77); W_05_seven_half_odd[20] = (-1,  103360,  1001)
    W_05_seven_half_odd[21] = (-1, 28424,143143)

    W_05_seven_half_even[ 1] = ( 0,     1,     1); W_05_seven_half_even[ 2] = ( 0,       1,     1)
    W_05_seven_half_even[ 3] = ( 0,     1,     1); W_05_seven_half_even[ 4] = ( 0,       1,     1)
    W_05_seven_half_even[ 5] = ( 0,     1,     1); W_05_seven_half_even[ 6] = ( 0,       1,     1)
    W_05_seven_half_even[ 7] = ( 0,     1,     1); W_05_seven_half_even[ 8] = ( 0,       1,     1)
    W_05_seven_half_even[ 9] = ( 0,     1,     1); W_05_seven_half_even[10] = ( 0,       1,     1)
    W_05_seven_half_even[11] = ( 1,   390,     7); W_05_seven_half_even[12] = ( 0,       1,     1)
    W_05_seven_half_even[13] = ( 0,     1,     1); W_05_seven_half_even[14] = (-1,      70,     1)
    W_05_seven_half_even[15] = ( 0,     1,     1); W_05_seven_half_even[16] = ( 0,       1,     1)
    W_05_seven_half_even[17] = ( 0,     1,     1); W_05_seven_half_even[18] = ( 1,      32,     7)
    W_05_seven_half_even[19] = ( 1,    14,     1); W_05_seven_half_even[20] = ( 0,       1,     1)
    W_05_seven_half_even[21] = ( 0,     1,     1); W_05_seven_half_even[22] = (-1,     576,    77)
    W_05_seven_half_even[23] = ( 0,     1,     1); W_05_seven_half_even[24] = ( 0,       1,     1)
    W_05_seven_half_even[25] = (-1,   210,    11); W_05_seven_half_even[26] = ( 0,       1,     1)
    W_05_seven_half_even[27] = ( 1, 63888,  1183); W_05_seven_half_even[28] = (-1,     154,    13)
    W_05_seven_half_even[29] = ( 0,     1,     1); W_05_seven_half_even[30] = (-1,    8568,   169)
    W_05_seven_half_even[31] = (-1,   176,     7); W_05_seven_half_even[32] = ( 0,       1,     1)
    W_05_seven_half_even[33] = ( 1,  1938,    91); W_05_seven_half_even[34] = ( 1,    1088,    11)
    W_05_seven_half_even[35] = ( 0,     1,     1); W_05_seven_half_even[36] = (-1,   28424,  1183)

    
"""
`SpinAngular.W_16_seven_half_odd[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(16) :: j QJ)  j = 7/2  for odd seniority
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_16_seven_half_odd = Array{Tuple{Int64,Int64,Int64}}(undef,21);    for i=1:21  W_16_seven_half_odd[i] = (0,0,0)   end


"""
`SpinAngular.W_16_seven_half_even[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(16) :: j QJ)  j = 7/2  for even seniority
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_16_seven_half_even = Array{Tuple{Int64,Int64,Int64}}(undef,36);   for i=1:36    W_16_seven_half_even[i] = (0,0,0)   end 

    W_16_seven_half_odd[ 1] = ( 0,     1,     1); W_16_seven_half_odd[ 2] = ( 0,       1,     1)
    W_16_seven_half_odd[ 3] = ( 0,     1,     1); W_16_seven_half_odd[ 4] = ( 1,     576,    11)
    W_16_seven_half_odd[ 5] = ( 1,   390,    77); W_16_seven_half_odd[ 6] = (-1,     144,     7)
    W_16_seven_half_odd[ 7] = ( 0,     1,     1); W_16_seven_half_odd[ 8] = ( 1,     520,    11)
    W_16_seven_half_odd[ 9] = (-1,    49,   121); W_16_seven_half_odd[10] = (-1,   12480,   121)
    W_16_seven_half_odd[11] = (-1,   408,    11); W_16_seven_half_odd[12] = ( 1,     520,     3)
    W_16_seven_half_odd[13] = (-1,  1960,    33); W_16_seven_half_odd[14] = (-1,    1664,    11)
    W_16_seven_half_odd[15] = ( 1,  3264,    11); W_16_seven_half_odd[16] = ( 1,  552250,  4719)
    W_16_seven_half_odd[17] = (-1, 43520,   847); W_16_seven_half_odd[18] = (-1,   38760, 11011)
    W_16_seven_half_odd[19] = (-1,  2652,   121); W_16_seven_half_odd[20] = ( 1,   15504,   121)
    W_16_seven_half_odd[21] = ( 1, 38760,   143)

    W_16_seven_half_even[ 1] = ( 0,     1,     1); W_16_seven_half_even[ 2] = ( 0,       1,     1)
    W_16_seven_half_even[ 3] = ( 0,     1,     1); W_16_seven_half_even[ 4] = ( 0,       1,     1)
    W_16_seven_half_even[ 5] = ( 0,     1,     1); W_16_seven_half_even[ 6] = ( 0,       1,     1)
    W_16_seven_half_even[ 7] = (-1,   130,     1); W_16_seven_half_even[ 8] = ( 0,       1,     1)
    W_16_seven_half_even[ 9] = ( 0,     1,     1); W_16_seven_half_even[10] = ( 0,       1,     1)
    W_16_seven_half_even[11] = ( 1,   390,    11); W_16_seven_half_even[12] = (-1,      48,     1)
    W_16_seven_half_even[13] = (-1,   234,     7); W_16_seven_half_even[14] = ( 1,    1040,    11)
    W_16_seven_half_even[15] = ( 1,   340,     7); W_16_seven_half_even[16] = ( 0,       1,     1)
    W_16_seven_half_even[17] = ( 1,  3120,   121); W_16_seven_half_even[18] = ( 0,       1,     1)
    W_16_seven_half_even[19] = ( 0,     1,     1); W_16_seven_half_even[20] = (-1,    1170,   121)
    W_16_seven_half_even[21] = ( 0,     1,     1); W_16_seven_half_even[22] = ( 1,   18720,   121)
    W_16_seven_half_even[23] = (-1,    36,    11); W_16_seven_half_even[24] = ( 1,     858,     7)
    W_16_seven_half_even[25] = (-1,  5304,   121); W_16_seven_half_even[26] = (-1,   69768,   847)
    W_16_seven_half_even[27] = ( 0,     1,     1); W_16_seven_half_even[28] = ( 0,       1,     1)
    W_16_seven_half_even[29] = ( 1,  1020,    11); W_16_seven_half_even[30] = ( 0,       1,     1)
    W_16_seven_half_even[31] = ( 0,     1,     1); W_16_seven_half_even[32] = ( 0,       1,     1)
    W_16_seven_half_even[33] = ( 0,     1,     1); W_16_seven_half_even[34] = (-1,   33592,   121)
    W_16_seven_half_even[35] = ( 1, 31654,   121); W_16_seven_half_even[36] = ( 0,       1,     1)
    

"""
`SpinAngular.W_07_seven_half_odd[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(07) :: j QJ)  j = 7/2  for odd seniority
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_07_seven_half_odd = Array{Tuple{Int64,Int64,Int64}}(undef,21);    for i=1:21  W_07_seven_half_odd[i] = (0,0,0)   end


"""
`SpinAngular.W_07_seven_half_even[]`
    ... set-up the ivalue of reduced matrix elements
        ( j QJ :: W(07) :: j QJ)  j = 7/2  for even seniority
        G. Gaigalas, S. Fritzsche, and Z. Rudzikas
        Atomic Data and Nuclear Data Tables 76, 235–269 (2000)
        doi:10.1006/adnd.2000.0844
        Please, contact Gediminas Gaigalas if you have problems/questions.
"""
W_07_seven_half_even = Array{Tuple{Int64,Int64,Int64}}(undef,36);   for i=1:36    W_07_seven_half_even[i] = (0,0,0)   end 

    W_07_seven_half_odd[ 1] = ( 0,     1,     1); W_07_seven_half_odd[ 2] = ( 0,       1,     1)
    W_07_seven_half_odd[ 3] = ( 0,     1,     1); W_07_seven_half_odd[ 4] = ( 0,       1,     1)
    W_07_seven_half_odd[ 5] = ( 1,   162,     7); W_07_seven_half_odd[ 6] = ( 1,     272,     7)
    W_07_seven_half_odd[ 7] = ( 0,     1,     1); W_07_seven_half_odd[ 8] = ( 0,       1,     1)
    W_07_seven_half_odd[ 9] = ( 1, 11025,  1573); W_07_seven_half_odd[10] = ( 1,    4624,   121)
    W_07_seven_half_odd[11] = (-1,  1632,   143); W_07_seven_half_odd[12] = (-1,     120,     1)
    W_07_seven_half_odd[13] = ( 0,     1,     1); W_07_seven_half_odd[14] = ( 0,       1,     1)
    W_07_seven_half_odd[15] = ( 0,     1,     1); W_07_seven_half_odd[16] = ( 1, 1224510, 20449)
    W_07_seven_half_odd[17] = ( 1,306000, 11011); W_07_seven_half_odd[18] = ( 1,12558240,143143)
    W_07_seven_half_odd[19] = (-1,  6460,   121); W_07_seven_half_odd[20] = ( 1,   77520,  1573)
    W_07_seven_half_odd[21] = (-1,297160,  1859)

    W_07_seven_half_even[ 1] = ( 0,     1,     1); W_07_seven_half_even[ 2] = ( 0,       1,     1)
    W_07_seven_half_even[ 3] = ( 0,     1,     1); W_07_seven_half_even[ 4] = ( 0,       1,     1)
    W_07_seven_half_even[ 5] = ( 0,     1,     1); W_07_seven_half_even[ 6] = ( 0,       1,     1)
    W_07_seven_half_even[ 7] = ( 0,     1,     1); W_07_seven_half_even[ 8] = ( 0,       1,     1)
    W_07_seven_half_even[ 9] = ( 0,     1,     1); W_07_seven_half_even[10] = ( 0,       1,     1)
    W_07_seven_half_even[11] = ( 0,     1,     1); W_07_seven_half_even[12] = ( 0,       1,     1)
    W_07_seven_half_even[13] = ( 0,     1,     1); W_07_seven_half_even[14] = ( 1,      60,     1)
    W_07_seven_half_even[15] = ( 0,     1,     1); W_07_seven_half_even[16] = ( 0,       1,     1)
    W_07_seven_half_even[17] = ( 0,     1,     1); W_07_seven_half_even[18] = ( 0,       1,     1)
    W_07_seven_half_even[19] = ( 1,  1600,    77); W_07_seven_half_even[20] = ( 0,       1,     1)
    W_07_seven_half_even[21] = ( 1,  2040,    77); W_07_seven_half_even[22] = ( 1,    4410,   121)
    W_07_seven_half_even[23] = ( 0,     1,     1); W_07_seven_half_even[24] = ( 0,       1,     1)
    W_07_seven_half_even[25] = ( 1, 18360,   121); W_07_seven_half_even[26] = ( 0,       1,     1)
    W_07_seven_half_even[27] = (-1, 18816,   845); W_07_seven_half_even[28] = (-1,   11016,   455)
    W_07_seven_half_even[29] = ( 0,     1,     1); W_07_seven_half_even[30] = ( 1,   11628,  1183)
    W_07_seven_half_even[31] = ( 1,    34,     5); W_07_seven_half_even[32] = ( 0,       1,     1)
    W_07_seven_half_even[33] = ( 1,  7752,   143); W_07_seven_half_even[34] = ( 1,    9690,   121)
    W_07_seven_half_even[35] = ( 0,     1,     1); W_07_seven_half_even[36] = ( 1,  222870,  1859) 
