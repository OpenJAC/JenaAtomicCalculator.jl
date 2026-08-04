#
println("Gc) Symbolic simplification of recoupling coefficients by means of sum rules.")


using SymEngine
j1  = Basic(:j1);    j2  = Basic(:j2);    j3  = Basic(:j3);    j4  = Basic(:j4);    j5  = Basic(:j5)
j6  = Basic(:j6);    j7  = Basic(:j7);    j8  = Basic(:j8);    j9  = Basic(:j9);    j10 = Basic(:j10)
j11 = Basic(:j11);   j12 = Basic(:j12);   j13 = Basic(:j13);   j14 = Basic(:j14);   j15 = Basic(:j15)
j16 = Basic(:j16);   j17 = Basic(:j17);   j18 = Basic(:j18);   j19 = Basic(:j19);   j20 = Basic(:j20)
j21 = Basic(:j21);   j22 = Basic(:j22);   j23 = Basic(:j23);   j24 = Basic(:j24);   j25 = Basic(:j25)

L1  = Basic(:L1);    L2  = Basic(:L2);    L3  = Basic(:L3);    L12 = Basic(:L12);    L123 = Basic(:L123);   
S1  = Basic(:S1);    S2  = Basic(:S2);    S3  = Basic(:S3);    S12 = Basic(:S12);    S123 = Basic(:S123);
T1  = Basic(:T1);    T2  = Basic(:T2);    T3  = Basic(:T3);    T12 = Basic(:T12)
Jm1 = Basic(:Jm1);   Jp1 = Basic(:Jp1);   J1  = Basic(:J1);    Jm2 = Basic(:Jm2);   Jp2 = Basic(:Jp2);   J2  = Basic(:J2)
Jm3 = Basic(:Jm3);   Jp3 = Basic(:Jp3);   J3  = Basic(:J3);    J12 = Basic(:J12)
Xp1 = Basic(:Xp1);   Xm2 = Basic(:Xm2);   Xp2 = Basic(:Xp2);   Xm3 = Basic(:Xm3);   Xp3 = Basic(:Xp3);
J   = Basic(:J)

if  false
    # Last successful:  unknown ...
    # Compute 
    leftCsq  = RacahAlgebra.Csq( j1, j2, j3)
    wa       = RacahAlgebra.evaluate(leftCsq, leftCsq)   
    #
elseif false
    # Last successful:  unknown ...
    # Compute 
    leftCsq  = RacahAlgebra.Csq( RacahAlgebra.Csq( j1, j2, j12), j3, J)
    rightCsq = RacahAlgebra.Csq( j1, RacahAlgebra.Csq( j2, j3, j23), J)
    wa       = RacahAlgebra.evaluate(leftCsq, rightCsq)    
    #
elseif false
    # Last successful:  unknown ...
    # LSJ-recoupling, two open shells
    leftCsq  = RacahAlgebra.Csq( RacahAlgebra.Csq( L1, L2, L12), RacahAlgebra.Csq( S1, S2, S12), J)
    rightCsq = RacahAlgebra.Csq( RacahAlgebra.Csq( L1, S1, T1),  RacahAlgebra.Csq( L2, S2, T2),  J)
    @time  wa       = RacahAlgebra.evaluate(leftCsq, rightCsq)    
    #
elseif false
    # Last successful:  unknown ...
    # jjJ-recoupling, two open shells
    leftCsq  = RacahAlgebra.Csq( RacahAlgebra.Csq( Jm1, Jp1, J1),   RacahAlgebra.Csq( Jm2, Jp2, J2), J)
    rightCsq = RacahAlgebra.Csq( RacahAlgebra.Csq( RacahAlgebra.Csq( Jm1, Jp1, Xp1), Jm2, Xm2), Jp2, J)
    wa       = RacahAlgebra.evaluate(leftCsq, rightCsq) 
    #
elseif false
    # Last successful:  unknown ...
    # LSJ-recoupling, three open shells
    leftCsq  = RacahAlgebra.Csq( RacahAlgebra.Csq( RacahAlgebra.Csq( L1, L2, L12), L3, L123), 
                                 RacahAlgebra.Csq( RacahAlgebra.Csq( S1, S2, S12), S3, S123), J)
    rightCsq = RacahAlgebra.Csq( RacahAlgebra.Csq( RacahAlgebra.Csq( L1, S1, T1),  RacahAlgebra.Csq( L2, S2, T2), T12),
                                                   RacahAlgebra.Csq( L3, S3, T3),   J) 
    wa       = RacahAlgebra.evaluate(leftCsq, rightCsq)    
    #
elseif true
    # Last successful:  unknown ...
    # jjJ-recoupling, three open shells
    leftCsq  = RacahAlgebra.Csq( RacahAlgebra.Csq( RacahAlgebra.Csq( Jm1, Jp1, J1),   RacahAlgebra.Csq( Jm2, Jp2, J2), J12),
                                                   RacahAlgebra.Csq( Jm3, Jp3, J3), J)  
    rightCsq = RacahAlgebra.Csq( RacahAlgebra.Csq( RacahAlgebra.Csq( RacahAlgebra.Csq( 
                                 RacahAlgebra.Csq( Jm1, Jp1, Xp1), Jm2, Xm2), Jp2, Xp2), Jm3, Xm3), Jp3, J)
    wa       = RacahAlgebra.evaluate(leftCsq, rightCsq)   
    #
elseif true
    # Last successful:  unknown ...
    # Compute 
    leftCsq  = RacahAlgebra.Csq( RacahAlgebra.Csq(j1,j2,j5), RacahAlgebra.Csq(j3,j4,j6), j7 )
    rightCsq = RacahAlgebra.Csq( j1, RacahAlgebra.Csq( RacahAlgebra.Csq(j2,j3,j8), j4, j9), j7 )
    wa       = RacahAlgebra.evaluate(leftCsq, rightCsq)
    #
elseif true
    # Last successful:  unknown ...
    # Compute 
    leftCsq  = RacahAlgebra.Csq( RacahAlgebra.Csq( RacahAlgebra.Csq(j1,j2,j9),j3,j10), 
                                 RacahAlgebra.Csq( RacahAlgebra.Csq( RacahAlgebra.Csq(j4,j5,j11),j6,j12), RacahAlgebra.Csq(j7,j8,j13),j14), j15)
    rightCsq = RacahAlgebra.Csq( RacahAlgebra.Csq( RacahAlgebra.Csq(j1,j4,j16), j7, j17), 
                                 RacahAlgebra.Csq( RacahAlgebra.Csq( RacahAlgebra.Csq(j2,j5,j18), RacahAlgebra.Csq(j8,j3,j19), j20), j6,j21), j15)
    wa       = RacahAlgebra.evaluate(leftCsq, rightCsq)
    #
end
