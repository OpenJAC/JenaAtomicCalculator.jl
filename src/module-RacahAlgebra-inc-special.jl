

"""
`RacahAlgebra.specialValue(w3j::RacahAlgebra.W3j)`  
    ... attempts to find a special value for the Wigner 3j symbol w3j. 
        A (istrue, rex)::Tuple{Bool, RacahExpression} is returned, and where istrue determined of whether a 
        special value is returned in rex. For istrue = false, rex has no meaning.
"""
function specialValue(w3j::RacahAlgebra.W3j)
    deltas = Kronecker[];    triangles = Triangle[];   integrals = Integral[]
    w3js = W3j[];   w6js = W6j[];   w9js = W9j[];   ylms = Ylm[];    djpqs = Djpq[]

    rexList = RacahAlgebra.symmetricForms(w3j)
    for  rex in rexList
        ww = rex.w3js[1]
        #
        #  Rule:        ( j1  j2  j3 )        J/2 [(J-2*j1)!(J-2*j2)!(J-2*j3)!] 1/2                   (J/2)!
        #  -----        (            )  = (-1)    [---------------------------]       -----------------------------------
        #               (  0   0   0 )            [          (J+1)!           ]       (J/2 - j1)! (J/2 - j2)! (J/2 - j3)!
        #
        #               with  J = j1+j2+j3  is  even; the 3j symbol is zero if J is odd.
        #
        specialW3j = W3j(ww.ja, ww.jb, ww.jc, 0, 0, 0);     J = ww.ja + ww.jb + ww.jc
        if  ww == specialW3j
            push!( deltas, Kronecker(J, Basic(:even)) )
            wa = RacahExpression( rex.summations, integrals, rex.phase + J/2, 
                    rex.weight * sqrt( factorial(J - 2*ww.ja) * factorial(J - 2*ww.jb) * factorial(J - 2*ww.jb) / factorial(J+1) ) * 
                                    factorial(J/2) / ( factorial(J/2-ww.ja) * factorial(J/2-ww.jb) * factorial(J/2-ww.jc) ), 
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs ) 
            return( (true, wa) )
        end
        #
        #  Rule:        ( j   j   0 )           j-m     -1/2
        #  -----        (           )   =   (-1)      [j]
        #               ( m  -m   0 )
        #
        specialW3j = W3j(ww.ja, ww.ja, 0, ww.ma, -ww.ma, 0)
        if  ww == specialW3j
            wa = RacahExpression( rex.summations, integrals, rex.phase + ww.ja - ww.ma, rex.weight / sqrt(2*ww.ja+1), 
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs ) 
            return( (true, wa) )
        end
        #
        #  Rule:        ( j   j-1/2   1/2 )           j-m-1   (   j - m   ) 1/2
        #  -----        (                 )   =   (-1)        (-----------)
        #               ( m  -m-1/2   1/2 )                   ( 2j (2j+1) )
        #
        specialW3j = W3j(ww.ja, ww.ja-1//2, 1//2, ww.ma, -ww.ma-1//2, 1//2)
        if  ww == specialW3j
            wa = RacahExpression( rex.summations, integrals, rex.phase + ww.ja - ww.ma - 1, 
                                    rex.weight * sqrt( (ww.ja - ww.ma)/(2*ww.ja * (2*ww.ja+1)) ), 
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs ) 
            return( (true, wa) )
        end
        #
        #  Rule:        ( j+1    j    1 )           j-m-1  (    (j-m)(j-m+1)    ) 1/2
        #  -----        (               )   =   (-1)       (--------------------)
        #               (  m   -m-1   1 )                  ( (2j+3)(2j+2)(2j+1) )
        #
        specialW3j = W3j(ww.jb+1, ww.jb, 1, ww.ma, -ww.ma-1, 1)
        if  ww == specialW3j
            wa = RacahExpression( rex.summations, integrals, rex.phase + ww.jb - ww.ma - 1, 
                                    rex.weight * sqrt( (ww.jb - ww.ma)*(ww.jb - ww.ma + 1)/(2*ww.jb+3) * (2*ww.jb+2) * (2*ww.jb+1) ), 
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs ) 
            return( (true, wa) )
        end
        #
        #  Rule:        ( j+1   j   1 )           j-m-1  (  2(j+m+1)(j-m+1)   ) 1/2
        #  -----        (             )   =   (-1)       (--------------------)
        #               (  m   -m   0 )                  ( (2j+3)(2j+2)(2j+1) )
        #
        specialW3j = W3j(ww.jb+1, ww.jb, 1, ww.ma, -ww.ma, 0)
        if  ww == specialW3j
            wa = RacahExpression( rex.summations, integrals, rex.phase + ww.jb - ww.ma - 1, 
                    rex.weight * sqrt( 2*(ww.jb + ww.ma + 1)*(ww.jb - ww.ma + 1)/(2*ww.jb+3) * (2*ww.jb+2) * (2*ww.jb+1) ), 
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs ) 
            return( (true, wa) )
        end
        #
        #  Rule:        ( j    j     1 )           j-m  (  2*(j-m)(j+m+1)  ) 1/2
        #  -----        (              )   =   (-1)     (------------------)
        #               ( m   -m-1   1 )                ( (2j+2)(2j+1)(2j) )
        #
        specialW3j = W3j(ww.ja, ww.ja, 1, ww.ma, -ww.ma-1, 1)
        if  ww == specialW3j
            wa = RacahExpression( rex.summations, integrals, rex.phase + ww.ja - ww.ma, 
                    rex.weight * sqrt( 2*(ww.ja - ww.ma)*(ww.ja + ww.ma + 1)/(2*ww.ja+2) * (2*ww.ja+1) * (2*ww.ja) ), 
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs ) 
            return( (true, wa) )
        end
        #
        #  Rule:        ( j   j   1 )           j-m            2 m
        #  -----        (           )   =   (-1)      -----------------------
        #               ( m  -m   0 )                 ((2j+2)(2j+1)(2j))^(1/2)
        #
        specialW3j = W3j(ww.ja, ww.ja, 1, ww.ma, -ww.ma, 0)
        if  ww == specialW3j
            wa = RacahExpression( rex.summations, integrals, rex.phase + ww.ja - ww.ma, 
                                    rex.weight * 2*ww.ma / sqrt( (2*ww.ja+2) * (2*ww.ja+1) * (2*ww.ja) ), 
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs ) 
            return( (true, wa) )
        end
        #
        #  Rule: 1      ( j+3/2    j    3/2 )           j-m+1/2   ( (j-m-1/2)(j-m+1/2)(j-m+3/2) ) 1/2
        #  -------      (                   )   =   (-1)          (-----------------------------)
        #               (   m   -m-3/2  3/2 )                     (   (2j+4)(2j+3)(2j+2)(2j+1)  )
        #
        specialW3j = W3j(ww.jb+3//2, ww.jb, 3//2, ww.ma, -ww.ma-3//2, 3//2)
        if  ww == specialW3j
            wa = RacahExpression( rex.summations, integrals, rex.phase + ww.jb - ww.ma + 1//2, 
                                    rex.weight * sqrt( ( (ww.jb-ww.ma-1//2) * (ww.jb-ww.ma+1//2) * (ww.jb-ww.ma-3//2) )/
                                                        ( (2*ww.jb+4) * (2*ww.jb+3) * (2*ww.jb+2) * (2*ww.jb+1) ) ), 
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs ) 
            return( (true, wa) )
        end
        #
        #  Rule: 2      ( j+3/2    j    3/2 )           j-m+1/2 ( 3*(j-m+1/2)(j-m+3/2)(j+m+3/2) ) 1/2
        #  -------      (                   )   =   (-1)        (-------------------------------)
        #               (   m   -m-1/2  1/2 )                   (    (2j+4)(2j+3)(2j+2)(2j+1)   )
        #
        specialW3j = W3j(ww.jb+3//2, ww.jb, 3//2, ww.ma, -ww.ma-1//2, 1//2)
        if  ww == specialW3j
            wa = RacahExpression( rex.summations, integrals, rex.phase + ww.jb - ww.ma + 1//2, 
                                    rex.weight * sqrt( ( 3*(ww.jb-ww.ma+1//2) * (ww.jb-ww.ma+3//2) * (ww.jb+ww.ma+3//2) )/
                                                        ( (2*ww.jb+4) * (2*ww.jb+3) * (2*ww.jb+2) * (2*ww.jb+1) ) ), 
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs ) 
            return( (true, wa) )
        end
        #
        #  Rule: 4      ( j+1/2    j    3/2 )            j-m-1/2           (        j-m+1/2         ) 1/2
        #  -------      (                   )    =   (-1)      (j+3*m+3/2) (------------------------)
        #               (   m   -m-1/2  1/2 )                              ( (2j+3)(2j+2)(2j+1)(2j) )
        #
        #
        specialW3j = W3j(ww.jb+1//2, ww.jb, 3//2, ww.ma, -ww.ma-1//2, 1//2)
        if  ww == specialW3j
            wa = RacahExpression( rex.summations, integrals, rex.phase + ww.jb - ww.ma - 1//2, 
                                    rex.weight * (ww.jb + 3*ww.ma + 3//2) *
                                    sqrt( (ww.jb-ww.ma+1//2) / ((2*ww.jb+3) * (2*ww.jb+2) * (2*ww.jb+1) * 2*ww.jb) ), 
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs ) 
            return( (true, wa) )
        end
        #
        #  Rule: 5      ( j+2    j   2 )           j-m  (   (j-m-1)(j-m)(j-m+1)(j-m+2)   ) 1/2
        #  -------      (              )   =   (-1)     (--------------------------------)
        #               ( m    -m-2  2 )                ( (2j+5)(2j+4)(2j+3)(2j+2)(2j+1) )
        #
        specialW3j = W3j(ww.jb+2, ww.jb, 2, ww.ma, -ww.ma-2, 2)
        if  ww == specialW3j
            wa = RacahExpression( rex.summations, integrals, rex.phase + ww.jb - ww.ma, 
                                    rex.weight * sqrt( (ww.jb-ww.ma-1) * (ww.jb-ww.ma) * (ww.jb-ww.ma+1) *(ww.jb-ww.ma-2) / 
                                                        ((2*ww.jb+5) * (2*ww.jb+4) * (2*ww.jb+3) * (2*ww.jb+2)  * (2*ww.jb+1) ) ), 
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs ) 
            return( (true, wa) )
        end
        #
        #  Rule: 7      ( j+2   j   2 )          j-m  ( 6 (j+m+2)(j+m+1)(j-m+2)(j-m+1) ) 1/2
        #  -------      (             )  =   (-1)     (--------------------------------)
        #               ( m    -m   0 )               ( (2j+5)(2j+4)(2j+3)(2j+2)(2j+1) )
        #
        specialW3j = W3j(ww.jb+2, ww.jb, 2, ww.ma, -ww.ma, 0)
        if  ww == specialW3j
            wa = RacahExpression( rex.summations, integrals, rex.phase + ww.jb - ww.ma, 
                                    rex.weight * sqrt( 6* (ww.jb+ww.ma+2) * (ww.jb+ww.ma+1) * (ww.jb-ww.ma+2) * (ww.jb-ww.ma+1) / 
                                                        ((2*ww.jb+5) * (2*ww.jb+4) * (2*ww.jb+3) * (2*ww.jb+2)  * (2*ww.jb+1) ) ), 
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs ) 
            return( (true, wa) )
        end
        #
        #  Rule: 8      ( j+1    j    2 )         j-m+1    (  (j-m-1)(j-m)(j-m+1)(j+m+2)  ) 1/2
        #  -------      (               ) =   (-1)      2* (------------------------------)
        #               ( m    -m-2   2 )                  ( (2j+4)(2j+3)(2j+2)(2j+1)(2j) )
        #
        specialW3j = W3j(ww.jb+1, ww.jb, 2, ww.ma, -ww.ma-2, 2)
        if  ww == specialW3j
            wa = RacahExpression( rex.summations, integrals, rex.phase + ww.jb - ww.ma + 1, 
                                    rex.weight * 2* sqrt( (ww.jb-ww.ma-1) * (ww.jb-ww.ma) * (ww.jb-ww.ma+1) * (ww.jb-ww.ma+2) / 
                                                        ((2*ww.jb+4) * (2*ww.jb+3) * (2*ww.jb+2) * (2*ww.jb+1)  * 2*ww.jb ) ), 
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs ) 
            return( (true, wa) )
        end
        #
        #  Rule: 9      ( j+1    j    2 )         j-m+1            (        (j-m+1)(j-m)          ) 1/2
        #  -------      (               ) =   (-1)     2*(j+2*m+2) (------------------------------)
        #               ( m    -m-1   1 )                          ( (2j+4)(2j+3)(2j+2)(2j+1)(2j) )
        #
        specialW3j = W3j(ww.jb+1, ww.jb, 2, ww.ma, -ww.ma-1, 1)
        if  ww == specialW3j
            wa = RacahExpression( rex.summations, integrals, rex.phase + ww.jb - ww.ma + 1, 
                                    rex.weight * 2* (ww.jb+2*ww.ma+2) * sqrt( (ww.jb-ww.ma-1) * (ww.jb-ww.ma) / 
                                                        ((2*ww.jb+4) * (2*ww.jb+3) * (2*ww.jb+2) * (2*ww.jb+1)  * 2*ww.jb ) ), 
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs ) 
            return( (true, wa) )
        end
        #
        #  Rule: 10     ( j+1   j   2 )          j-m+1     (       6 (j+m+1)(j-m+1)       ) 1/2
        #  --------     (             )  =   (-1)      2*m (------------------------------)
        #               ( m    -m   0 )                    ( (2j+4)(2j+3)(2j+2)(2j+1)(2j) )
        #
        specialW3j = W3j(ww.jb+1, ww.jb, 2, ww.ma, -ww.ma, 0)
        if  ww == specialW3j
            wa = RacahExpression( rex.summations, integrals, rex.phase + ww.jb - ww.ma + 1, 
                                    rex.weight * 2*ww.ma * sqrt( 6* (ww.jb-ww.ma+1) * (ww.jb-ww.ma+1) / 
                                                        ((2*ww.jb+4) * (2*ww.jb+3) * (2*ww.jb+2) * (2*ww.jb+1)  * 2*ww.jb ) ), 
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs ) 
            return( (true, wa) )
        end
        #
        #  Rule: 13     ( j   j   2 )          j-m         2 (3*m*m - j(j+1))
        #  --------     (           )  =   (-1)     ------------------------------------
        #               ( m  -m   0 )               ((2j+3)(2j+2)(2j+1)(2j)(2j-1))^(1/2)
        #
        specialW3j = W3j(ww.jb, ww.jb, 2, ww.ma, -ww.ma, 0)
        if  ww == specialW3j
            wa = RacahExpression( rex.summations, integrals, rex.phase + ww.jb - ww.ma, 
                                    rex.weight * 2* (3*ww.ma*ww.ma - ww.jb*(ww.jb+1)) / 
                                                sqrt( (2*ww.jb+3) * (2*ww.jb+2) * (2*ww.jb+1)  * 2*ww.jb * (2*ww.jb-1) ),
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs ) 
            return( (true, wa) )
        end
    end
    
    return( (false, RacahExpression() ) )
end


"""
`RacahAlgebra.specialValue(w6j::RacahAlgebra.W6j)`  
    ... attempts to find a special value for the Wigner 6-j symbol w6j. 
        A (istrue, rex)::Tuple{Bool, RacahExpression} is returned, and where istrue determined of whether a 
        special value is returned in rex. For istrue = false, rex has no meaning.
"""
function specialValue(w6j::RacahAlgebra.W6j)
    deltas = Kronecker[];    triangles = Triangle[];   integrals = Integral[]
    w3js = W3j[];   w6js = W6j[];   w9js = W9j[];   ylms = Ylm[];    djpqs = Djpq[]

    rexList = RacahAlgebra.symmetricForms(w6j)
    for  rex in rexList
        ww = rex.w6js[1]
        #
        #                                            j1+j2+j3
        #  Rule:         ( j1   j2   j3 )        (-1)
        #  -----        {(              )}   =  --------------  t(j1,j2,j3)  d(j1,l2) d(j2,l1)
        #                ( l1   l2   0  )                   1/2
        #                                         [ j1, j2 ]
        #
        specialW6j = W6j( ww.a, ww.b, ww.c, ww.d, ww.e, 0)
        if  ww == specialW6j
            push!( deltas, Kronecker(ww.a, ww.e) )
            push!( deltas, Kronecker(ww.b, ww.d) )
            push!( triangles, Triangle(ww.a, ww.b, ww.c) )
            wa = RacahExpression( rex.summations, integrals, rex.phase + ww.a + ww.b + ww.c, 
                                    rex.weight / sqrt( (2*ww.a+1)*(2*ww.b+1)), deltas, triangles, w3js, w6js, w9js, ylms, djpqs )
            return( (true, wa) )
        end
        #
        #  Rule:         ( j1     j2      j3   )            j1+j2+j3    (     (j1+j3-j2)(j1+j2-j3+1)     ) 1/2
        #  -----        {(                     )}   =   (-1)            (--------------------------------)
        #                ( 1/2  j3-1/2  j2+1/2 )                        ( (2*j2+1)(2*j2+2)(2*j3)(2*j3+1) )
        #
        specialW6j = W6j( ww.a, ww.b, ww.c, 1//2, ww.c -1//2,  ww.b +1//2)
        if  ww == specialW6j
            wa = RacahExpression( rex.summations, integrals, rex.phase + ww.a + ww.b + ww.c, 
                                    rex.weight * sqrt( (ww.a+ww.c-ww.b)*(ww.a+ww.b-ww.c+1) /
                                                        (2*ww.b+1)*(2*ww.b+2)*2*ww.c*(2*ww.c+1) ), 
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs )
            return( (true, wa) )
        end
        #
        #  Rule:         ( j1     j2      j3   )            j1+j2+j3+1   (     (j1+j2+j3+2)(j2+j3-j1+1)     ) 1/2
        #  -----        {(                     )}   =   (-1)             (----------------------------------)
        #                ( 1/2  j3+1/2  j2+1/2 )                         ( (2*j2+1)(2*j2+2)(2*j3+1)(2*j3+2) )
        #
        specialW6j = W6j( ww.a, ww.b, ww.c, 1//2, ww.c +1//2,  ww.b +1//2)
        if  ww == specialW6j
            wa = RacahExpression( rex.summations, integrals, rex.phase + ww.a + ww.b + ww.c + 1, 
                                    rex.weight * sqrt( (ww.a+ww.b+ww.c+2)*(ww.b+ww.c-ww.a+1) /
                                                        (2*ww.b+1)*(2*ww.b+2)*(2*ww.c+1)*(2*ww.c+2) ), 
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs )
            return( (true, wa) )
        end
        #
        #  Rule:         ( j1   j2    j3  )            S       (           S(S+1)(S-2*j1)(S-2*j1-1)           ) 1/2
        #  -----        {(                )}   =   (-1)        (----------------------------------------------)
        #                (  1  j3-1  j2-1 )                    ( (2*j2-1)(2*j2)(2*j2+1)(2*j3-1)(2*j3)(2*j3+1) )
        #
        #                with  S = j1+j2+j3
        #
        specialW6j = W6j( ww.a, ww.b, ww.c, 1, ww.c -1,  ww.b -1)
        if  ww == specialW6j
            S = ww.a + ww.b + ww.c
            wa = RacahExpression( rex.summations, integrals, rex.phase + S, 
                                    rex.weight * sqrt( S*(S+1)*(S-2*ww.a)*(S-2*ww.a -1) /
                                                        (2*ww.b-1)* 2*ww.b * (2*ww.b+1)* (2*ww.c-1)* 2*ww.c * (2*ww.c+1) ),
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs )
            return( (true, wa) )
        end
        #
        #  Rule:         ( j1   j2   j3  )            S   (      2(S+1)(S-2*j1)(S-2*j2)(S-2*j3+1)        ) 1/2
        #  -----        {(               )}   =   (-1)    (----------------------------------------------)
        #                (  1  j3-1  j2  )                ( (2*j2)(2*j2+1)(2*j2+2)(2*j3-1)(2*j3)(2*j3+1) )
        #
        #                with  S = j1+j2+j3
        #
        specialW6j = W6j( ww.a, ww.b, ww.c, 1, ww.c -1,  ww.b)
        if  ww == specialW6j
            S = ww.a + ww.b + ww.c
            wa = RacahExpression( rex.summations, integrals, rex.phase + S, 
                                    rex.weight * sqrt( 2*(S+1)*(S-2*ww.a)*(S-2*ww.b)*(S-2*ww.c +1) /
                                                        2*ww.b * (2*ww.b+1) * (2*ww.b+2) * (2*ww.c-1) * 2*ww.c * (2*ww.c+1) ),
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs )
            return( (true, wa) )
        end
        #
        #  Rule:         ( j1   j2   j3    )            S   (     (S-2*j2-1)(S-2*j2)(S-2*j3+1)(S-2*j3+2)     ) 1/2
        #  -----        {(                 )}   =   (-1)    (------------------------------------------------)
        #                (  1  j3-1  j2+1  )                ( (2*j2+1)(2*j2+2)(2*j2+3)(2*j3-1)(2*j3)(2*j3+1) )
        #
        #                with  S = j1+j2+j3
        #
        specialW6j = W6j( ww.a, ww.b, ww.c, 1, ww.c -1,  ww.b +1)
        if  ww == specialW6j
            S = ww.a + ww.b + ww.c
            wa = RacahExpression( rex.summations, integrals, rex.phase + S, 
                                    rex.weight * sqrt( (S-2*ww.b-1)*(S-2*ww.b)*(S-2*ww.c +1)*(S-2*ww.c +2) /
                                                        (2*ww.b+1) * (2*ww.b+2) * (2*ww.b+3) * (2*ww.c-1) * 2*ww.c * (2*ww.c+1) ),
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs )
            return( (true, wa) )
        end
        #
        #  Rule:         ( j1  j2  j3  )            S+1         2 [ j2(j2+1) + j3(j3+1) -  j1(j1+1) ]
        #  -----        {(             )}   =   (-1)       ----------------------------------------------  1/2
        #                (  1  j3  j2  )                   ( (2*j2)(2*j2+1)(2*j2+2)(2*j3)(2*j3+1)(2*j3+2) )
        #
        #                with  S = j1+j2+j3
        #
        specialW6j = W6j( ww.a, ww.b, ww.c, 1, ww.c,  ww.b)
        if  ww == specialW6j
            S = ww.a + ww.b + ww.c
            wa = RacahExpression( rex.summations, integrals, rex.phase + S + 1, 
                                    rex.weight * 2 * ( ww.b*(ww.b+1) + ww.c*(ww.c+1) - ww.a*(ww.a+1) ) /
                                                    sqrt( 2*ww.b * (2*ww.b+1) * (2*ww.b+2) * 2*ww.c * (2*ww.c+1) * (2*ww.c+2) ),
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs )
            return( (true, wa) )
        end
        #
        #  Rule: 1       (  a     b      c   )            s   (         (s-1)s(s+1)(s-2*a-2)(s-2*a-1)(s-2*a)         ) 1/2
        #  -----        {(                   )}   =   (-1)    (------------------------------------------------------)
        #                ( 3/2  c-3/2  b-3/2 )                ( (2*b-2)(2*b-1)(2*b)(2*b+1)(2*c-2)(2*c-1)(2*c)(2*c+1) )
        #
        #                with  s = a + b + c
        #
        specialW6j = W6j( ww.a, ww.b, ww.c, 3//2, ww.c-3//2,  ww.b-3//2)
        if  ww == specialW6j
            s = ww.a + ww.b + ww.c
            wa = RacahExpression( rex.summations, integrals, rex.phase + s, 
                                    rex.weight * sqrt( (s-1)* s * (s+1) * (s-2*ww.a-2) * (s-2*ww.a-1) * (s-2*ww.a)/
                                                        (2*ww.b-2)*(2*ww.b-1)*2*ww.b*(2*ww.b+1) * (2*ww.c-2)*(2*ww.c-1)*2*ww.c*(2*ww.c+1) ),
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs )
            return( (true, wa) )
        end
        #
        #  Rule: 2       (  a     b      c   )            s   (       3*s(s+1)(s-2*a-1)(s-2*a)(s-2*b)(s-2*c+1)       ) 1/2
        #  -----        {(                   )}   =   (-1)    (------------------------------------------------------)
        #                ( 3/2  c-3/2  b-1/2 )                ( (2*b-1)(2*b)(2*b+1)(2*b+2)(2*c-2)(2*c-1)(2*c)(2*c+1) )
        #
        #                with  s = a + b + c
        #
        specialW6j = W6j( ww.a, ww.b, ww.c, 3//2, ww.c-3//2,  ww.b-1//2)
        if  ww == specialW6j
            s = ww.a + ww.b + ww.c
            wa = RacahExpression( rex.summations, integrals, rex.phase + s, 
                                    rex.weight * sqrt( 3*s * (s+1) * (s-2*ww.a-1) * (s-2*ww.a) * (s-2*ww.b) * (s-2*ww.c+1)/
                                                        ((2*ww.b-1)*2*ww.b*(2*ww.b+1)*(2*ww.b+2) * (2*ww.c-2)*(2*ww.c-1)*2*ww.c*(2*ww.c+1)) ),
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs )
            return( (true, wa) )
        end
        #
        #  Rule: 3       (  a     b      c   )            s   (   3*(s+1)(s-2*a)(s-2*b-1)(s-2*b)(s-2*c+1)(s-2*c+2)   ) 1/2
        #  -----        {(                   )}   =   (-1)    (------------------------------------------------------)
        #                ( 3/2  c-3/2  b+1/2 )                ( (2*b)(2*b+1)(2*b+2)(2*b+3)(2*c-2)(2*c-1)(2*c)(2*c+1) )
        #
        #                with  s = a + b + c
        #
        specialW6j = W6j( ww.a, ww.b, ww.c, 3//2, ww.c-3//2,  ww.b+1//2)
        if  ww == specialW6j
            s = ww.a + ww.b + ww.c
            wa = RacahExpression( rex.summations, integrals, rex.phase + s, 
                                    rex.weight * sqrt( 3* (s+1) * (s-2*ww.a) * (s-2*ww.b-1) * (s-2*ww.b) * (s-2*ww.c+1) * (s-2*ww.c+2)/
                                                        (2*ww.b *(2*ww.b+1)*(2*ww.b+2)*(2*ww.b+3) * (2*ww.c-2)*(2*ww.c-1)*2*ww.c*(2*ww.c+1)) ),
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs )
            return( (true, wa) )
        end
        #
        #  Rule: 4       (  a     b      c   )            s    (  (s-2*b-2)(s-2*b-1)(s-2*b)(s-2*c+1)(s-2*c+2)(s-2*c+3)  ) 1/2
        #  -----        {(                   )}   =   (-1)     (--------------------------------------------------------)
        #                ( 3/2  c-3/2  b+3/2 )                 ( (2*b+1)(2*b+2)(2*b+3)(2*b+4)(2*c-2)(2*c-1)(2*c)(2*c+1) )
        #
        #                with  s = a + b + c
        #
        specialW6j = W6j( ww.a, ww.b, ww.c, 3//2, ww.c-3//2,  ww.b+3//2)
        if  ww == specialW6j
            s = ww.a + ww.b + ww.c
            wa = RacahExpression( rex.summations, integrals, rex.phase + s, 
                                    rex.weight * sqrt( (s-2*ww.b-2) * (s-2*ww.b-1) * (s-2*ww.b) * (s-2*ww.c+1) * (s-2*ww.c+2) * (s-2*ww.c+3)/
                                                        ((2*ww.b+1)*(2*ww.b+2)*(2*ww.b+3)*(2*ww.b+4) * (2*ww.c-2)*(2*ww.c-1)*2*ww.c*(2*ww.c+1)) ),
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs )
            return( (true, wa) )
        end
        #
        #  Rule: 5       (  a     b      c   )            s                                      (                     (s+1)(s-2*a)                     ) 1/2
        #  -----        {(                   )}   =  (-1)   [2*(s-2*b)(s-2*c) - (s+2)(s-2*a-1)]  (------------------------------------------------------)
        #                ( 3/2  c-1/2  b-1/2 )                                                   ( (2*b-1)(2*b)(2*b+1)(2*b+2)(2*c-1)(2*c)(2*c+1)(2*c+2) )
        #
        #      with  s = a + b + c
        #
        specialW6j = W6j( ww.a, ww.b, ww.c, 3//2, ww.c-1//2,  ww.b-1//2)
        if  ww == specialW6j
            s = ww.a + ww.b + ww.c
            wa = RacahExpression( rex.summations, integrals, rex.phase + s, 
                                    rex.weight * (2*(s-2*ww.b)*(s-2*ww.c) - (s+2)*(s-2*ww.a-1)) *   sqrt( (s+1)*(s-2*ww.a) /
                                                        ((2*ww.b-1)*2*ww.b*(2*ww.b+1)*(2*ww.b+2) * (2*ww.c-1)*2*ww.c*(2*ww.c+1)*(2*ww.c+1)) ),
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs )
            return( (true, wa) )
        end
        #
        #   Rule: 6       (  a     b      c   )            s                                     (                   (s-2*b)(s-2*c+1)                   ) 1/2
        #   -----        {(                   )}  =  (-1)  [ (s-2*b-1)(s-2*c) - 2*(s+2)(s-2*a)]  (------------------------------------------------------)
        #                 ( 3/2  c-1/2  b+1/2 )                                                  ( (2*b)(2*b+1)(2*b+2)(2*b+3)(2*c-1)(2*c)(2*c+1)(2*c+2) )
        #
        #                 with  s = a + b + c
        #
        specialW6j = W6j( ww.a, ww.b, ww.c, 3//2, ww.c-1//2,  ww.b+1//2)
        if  ww == specialW6j
            s = ww.a + ww.b + ww.c
            wa = RacahExpression( rex.summations, integrals, rex.phase + s, 
                                    rex.weight * ((s-2*ww.b-1)*(s-2*ww.c) - 2*(s+2)*(s-2*ww.a)) *   sqrt( (s-2*ww.b) * (s-2*ww.c + 1) /
                                                (2*ww.b * (2*ww.b+1)*(2*ww.b+2)*(2*ww.b+3) * (2*ww.c-1)*2*ww.c*(2*ww.c+1)*(2*ww.c+2)) ),
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs )
            return( (true, wa) )
        end
        #
        #   Rule: 7       ( a   b    c  )            s   (                  (s-2)(s-1)s(s+1)(s-2*a-3)(s-2*a-2)(s-2*a-1)(s-2*a) ) 1/2
        #   -----        {(             )}   =   (-1)    (---------------------------------------------------------------------)
        #                 ( 2  c-2  b-2 )                ( (2*b-3)(2*b-2)(2*b-1)(2*b)(2*b+1)(2*c-3)(2*c-2)(2*c-1)(2*c)(2*c+1)  )
        #
        #                 with  s = a + b + c
        #
        specialW6j = W6j( ww.a, ww.b, ww.c, 2, ww.c-2,  ww.b-2)
        if  ww == specialW6j
            s = ww.a + ww.b + ww.c
            wa = RacahExpression( rex.summations, integrals, rex.phase + s, 
                                    rex.weight * sqrt( (s-2) * (s-1) *s* (s-2*ww.a-3)* (s-2*ww.a-2)* (s-2*ww.a-1)* (s-2*ww.a) /
                        ((2*ww.b-3)*(2*ww.b-2)*(2*ww.b-1)*2*ww.b*(2*ww.b+1) * (2*ww.c-3)*(2*ww.c-2)*(2*ww.c-1)*2*ww.c*(2*ww.c+1)) ),
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs )
            return( (true, wa) )
        end
        #
        #   Rule: 8       ( a   b    c  )            s     (     (s-1)s(s+1)(s-2*a-2)(s-2*a-1)(s-2*a)(s-2*b)(s-2*c+1)            ) 1/2
        #   -----        {(             )}   =   (-1)    2 (---------------------------------------------------------------------)
        #                 ( 2  c-2  b-1 )                  ( (2*b-2)(2*b-1)(2*b)(2*b+1)(2*b+2)(2*c-3)(2*c-2)(2*c-1)(2*c)(2*c+1)  )
        #
        #                 with  s = a + b + c
        #
        specialW6j = W6j( ww.a, ww.b, ww.c, 2, ww.c-2,  ww.b-1)
        if  ww == specialW6j
            s = ww.a + ww.b + ww.c
            wa = RacahExpression( rex.summations, integrals, rex.phase + s, 
                                    rex.weight * 2 + sqrt( (s-1)*s*(s+1) * (s-2*ww.a-2) * (s-2*ww.a-1) * (s-2*ww.a) * (s-2*ww.b) * (s-2*ww.c+1)/
                        ((2*ww.b-2)*(2*ww.b-1)*2*ww.b*(2*ww.b+1)*(2*ww.b+2) * (2*ww.c-3)*(2*ww.c-2)*(2*ww.c-1)*2*ww.c*(2*ww.c+1)) ),
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs )
            return( (true, wa) )
        end
        #
        #   Rule: 9       ( a   b    c  )            s  (       6*s*(s+1)(s-2*a-1)(s-2*b-1)(s-2*a)(s-2*b)(s-2*c+1)(s-2*c+2)  ) 1/2
        #   -----        {(             )}   =   (-1)   (--------------------------------------------------------------------)
        #                 ( 2  c-2   b  )               ( (2*b-1)(2*b)(2*b+1)(2*b+2)(2*b+3)(2*c-3)(2*c-2)(2*c-1)(2*c)(2*c+1) )
        #
        #                 with  s = a + b + c
        #
        specialW6j = W6j( ww.a, ww.b, ww.c, 2, ww.c-2,  ww.b)
        if  ww == specialW6j
            s = ww.a + ww.b + ww.c
            wa = RacahExpression( rex.summations, integrals, rex.phase + s, 
                                    rex.weight * sqrt( 6*s * (s+1)*(s-2*ww.a-1)*(s-2*ww.b-1)*(s-2*ww.a)*(s-2*ww.b)*(s-2*ww.c+1)*(s-2*ww.c+2)/
                        ((2*ww.b-1)*2*ww.b*(2*ww.b+1)*(2*ww.b+2)*(2*ww.b+3) * (2*ww.c-3)*(2*ww.c-2)*(2*ww.c-1)*2*ww.c*(2*ww.c+1)) ),
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs )
            return( (true, wa) )
        end
        #
        #   Rule: 10      ( a   b    c  )            s    (    (s+1)(s-2*a)(s-2*b-2)(s-2*b-1)(s-2*b)(s-2*c+1)(s-2*c+2)(s-2*c+3) ) 1/2
        #   -----        {(             )}   =   (-1)   2 (---------------------------------------------------------------------)
        #                 ( 2  c-2  b+1 )                 ( (2*b)(2*b+1)(2*b+2)(2*b+3)(2*b+4)(2*c-3)(2*c-2)(2*c-1)(2*c)(2*c+1)  )
        #
        #                 with  s = a + b + c
        #
        specialW6j = W6j( ww.a, ww.b, ww.c, 2, ww.c-2,  ww.b+1)
        if  ww == specialW6j
            s = ww.a + ww.b + ww.c
            wa = RacahExpression( rex.summations, integrals, rex.phase + s, 
                                    rex.weight * 2 * sqrt( (s+1)*(s-2*ww.a)*(s-2*ww.b-2)*(s-2*ww.b-1)*(s-2*ww.b)*(s-2*ww.c+1)*(s-2*ww.c+2)*(s-2*ww.c+3) /
                                    ( 2*ww.b*(2*ww.b+1)*(2*ww.b+2)*(2*ww.b+3)*(2*ww.b+4) * (2*ww.c-3)*(2*ww.c-2)*(2*ww.c-1)*2*ww.c*(2*ww.c+1) ) ),
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs )
            return( (true, wa) )
        end
        #
        #   Rule: 11      ( a   b    c  )            s   (  (s-2*b-3)(s-2*b-2)(s-2*b-1)(s-2*b)(s-2*c+1)(s-2*c+2)(s-2*c+3)(s-2*c+4) ) 1/2
        #   -----        {(             )}   =   (-1)    (-------------------------------------------------------------------------)
        #                 ( 2  c-2  b+2 )                ( (2*b+1)(2*b+2)(2*b+3)(2*b+4)(2*b+5)(2*c-3)(2*c-2)(2*c-1)(2*c)(2*c+1)    )
        #
        #                 with  s = a + b + c
        #
        specialW6j = W6j( ww.a, ww.b, ww.c, 2, ww.c-2,  ww.b+2)
        if  ww == specialW6j
            s = ww.a + ww.b + ww.c
            wa = RacahExpression( rex.summations, integrals, rex.phase + s, 
                                    rex.weight * sqrt( (s-2*ww.b-3)*(s-2*ww.b-2)*(s-2*ww.b-1)*(s-2*ww.b) * (s-2*ww.c+1)*(s-2*ww.c+2)*(s-2*ww.c+3)*(s-2*ww.c+4) /
                                    ( (2*ww.b+1)*(2*ww.b+2)*(2*ww.b+3)*(2*ww.b+4)*(2*ww.b+5) * (2*ww.c-3)*(2*ww.c-2)*(2*ww.c-1)*2*ww.c*(2*ww.c+1) ) ),
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs )
            return( (true, wa) )
        end
        #
        #   Rule: 12      ( a   b    c  )        s                                    (                           s*(s+1)(s-2*a-1)(s-2*a)                  ) 1/2
        #   -----        {(             )} = (-1)   4 [ (a+b)(a-b+1) - (c-1)(c-b+1)]  (--------------------------------------------------------------------) 
        #                 ( 2  c-1  b-1 )                                             ( (2*b-2)(2*b-1)(2*b)(2*b+1)(2*b+2)(2*c-2)(2*c-1)(2*c)(2*c+1)(2*c+2) )
        #
        #      with  s = a + b + c
        #
        specialW6j = W6j( ww.a, ww.b, ww.c, 2, ww.c-1,  ww.b-1)
        if  ww == specialW6j
            s = ww.a + ww.b + ww.c
            wa = RacahExpression( rex.summations, integrals, rex.phase + s, 
                                    rex.weight * 4 * ( (ww.a+ww.b)*(ww.a-ww.b+1) - (ww.c-1)*(ww.c-ww.b+1) ) *
                                    sqrt( s*(s+1) * (s-2*ww.a-1)*(s-2*ww.a) /
                                        ( (2*ww.b-2)*(2*ww.b-1)*2*ww.b*(2*ww.b+1)*(2*ww.b+2) * (2*ww.c-2)*(2*ww.c-1)*2*ww.c*(2*ww.c+1)*(2*ww.c+2) ) ),
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs )
            return( (true, wa) )
        end
        #
        #   Rule: 13      ( a   b    c  )        s                                (                    6*(s+1)(s-2*a)(s-2*b)(s-2*c+1)                  ) 1/2
        #   -----        {(             )} = (-1)   2 [ (a+b+1)(a-b) - c*c  + 1 ] (--------------------------------------------------------------------)
        #                 ( 2  c-1   b  )                                         ( (2*b-1)(2*b)(2*b+1)(2*b+2)(2*b+3)(2*c-2)(2*c-1)(2*c)(2*c+1)(2*c+2) )
        #
        #                 with  s = a + b + c
        #
        specialW6j = W6j( ww.a, ww.b, ww.c, 2, ww.c-1, ww.b)
        if  ww == specialW6j
            s = ww.a + ww.b + ww.c
            wa = RacahExpression( rex.summations, integrals, rex.phase + s, 
                                    rex.weight * 2 * ( (ww.a+ww.b+1)*(ww.a-ww.b) - ww.c*ww.c + 1 ) *
                                    sqrt(  6*(s+1) * (s-2*ww.a) * (s-2*ww.b) * (s-2*ww.c+1) /
                                        ( (2*ww.b-1)*2*ww.b*(2*ww.b+1)*(2*ww.b+2)*(2*ww.b+3) * (2*ww.c-2)*(2*ww.c-1)*2*ww.c*(2*ww.c+1)*(2*ww.c+2) ) ),
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs )
            return( (true, wa) )
        end
        #
        #   Rule: 14      ( a   b    c  )        s                                     (                  (s-2*b-1)(s-2*b)(s-2*c+1)(s-2*c+2)                ) 1/2
        #   -----        {(             )} = (-1)  4 [ (a+b+2)(a-b-1) - (c-1)(b+c+2) ] (--------------------------------------------------------------------)
        #                 ( 2  c-1  b+1 )                                              ( (2*b)(2*b+1)(2*b+2)(2*b+3)(2*b+4)(2*c-2)(2*c-1)(2*c)(2*c+1)(2*c+2) )
        #
        #                 with  s = a + b + c
        #
        specialW6j = W6j( ww.a, ww.b, ww.c, 2, ww.c-1,  ww.b+1)
        if  ww == specialW6j
            s = ww.a + ww.b + ww.c
            wa = RacahExpression( rex.summations, integrals, rex.phase + s, 
                                    rex.weight * 4 * ( (ww.a+ww.b+2)*(ww.a-ww.b-1) - (ww.c-1)*(ww.b+ww.c+2) ) *
                                    sqrt( (s-2*ww.b-1)*(s-2*ww.b) * (s-2*ww.c+1)*(s-2*ww.c+2) /
                                        ( 2*ww.b*(2*ww.b+1)*(2*ww.b+2)*(2*ww.b+3)*(2*ww.b+4) * (2*ww.c-2)*(2*ww.c-1)*2*ww.c*(2*ww.c+1)*(2*ww.c+2) ) ),
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs )
            return( (true, wa) )
        end
        #
        #   Rule: 15      ( a   b   c )        s                                     (                                                                    ) 1/2
        #   -----        {(           )} = (-1)  2 [ 3*X(X-1) - 4*b*(b+1)*c*(c+1) ]  (--------------------------------------------------------------------)
        #                 ( 2   c   b )                                              ( (2*b-1)(2*b)(2*b+1)(2*b+2)(2*b+3)(2*c-1)(2*c)(2*c+1)(2*c+2)(2*c+3) )
        #
        #                 with  s = a + b + c   and   X = b*(b+1) + c*(c+1) - a*(a+1)
        #
        specialW6j = W6j( ww.a, ww.b, ww.c, 2, ww.c,  ww.b)
        if  ww == specialW6j
            s = ww.a + ww.b + ww.c;     X = ww.b*(ww.b+1) + ww.c*(ww.c+1) - ww.a*(ww.a+1)
            wa = RacahExpression( rex.summations, integrals, rex.phase + s, 
                                    rex.weight * 2 * ( 3*X*(X-1) - 4*ww.b*(ww.b+1)*ww.c*(ww.c+1) ) *
                                    sqrt( 1 / ( (2*ww.b-1)*2*ww.b*(2*ww.b+1)*(2*ww.b+2)*(2*ww.b+3) * (2*ww.c-1)*2*ww.c*(2*ww.c+1)*(2*ww.c+2)*(2*ww.c+3) ) ),
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs )
            return( (true, wa) )
        end
    end
    
    return( (false, RacahExpression() ) )
end


"""
`RacahAlgebra.specialValue(w9j::RacahAlgebra.W9j)`  
    ... attempts to find a special value for the Wigner 9-j symbol w9j. 
        A (istrue, rex)::Tuple{Bool, RacahExpression} is returned, and where istrue determined of whether a 
        special value is returned in rex. For istrue = false, rex has no meaning.
"""
function specialValue(w9j::RacahAlgebra.W9j)
    deltas = Kronecker[];    triangles = Triangle[];   integrals = Integral[]
    w3js= W3j[];   w6js = W6j[];   w9js = W9j[];   ylms = Ylm[];    djpqs = Djpq[]

    rexList = RacahAlgebra.symmetricForms(w9j)
    for  rex in rexList
        ww = rex.w9js[1]
        #
        #  Rule:         ( ja   ja   0 )
        #  -----         (             )         D(ja,jd,jg)
        #               {( jd   jd   0 )}   =  -------------- 1/2    delta(a,b)  delta(d,e)  delta(g,h)
        #                (             )       [ ja, jd, jg ]
        #                ( jg   jg   0 )
        #
        specialW9j = W9j( ww.a, ww.a, 0, ww.d, ww.d, 0, ww.g, ww.g, 0)
        if  ww == specialW9j
            push!( deltas, Kronecker(ww.a, ww.b) )
            push!( deltas, Kronecker(ww.d, ww.e) )
            push!( deltas, Kronecker(ww.g, ww.h) )
            push!( triangles, Triangle(ww.a, ww.d, ww.g) )
            wa = RacahExpression( rex.summations, integrals, rex.phase, rex.weight / sqrt( (2*ww.a+1)*(2*ww.d+1)*(2*ww.g+1)),
                                    deltas, triangles, w3js, w6js, w9js, ylms, djpqs )
            return( (true, wa) )
        end
        #
        #  Rule:         ( ja   jb   jc )           jb+jc+jd+jg
        #  -----         (              )       (-1)                   ( ja   jb   jc )
        #               {( jd   je   jf )}   =  --------------- 1/2   {(              )}
        #                (              )          [ jc, jg ]          ( je   jd   jg )
        #                ( jg   jh   0  )
        #
        specialW9j = W9j( ww.a, ww.b, ww.c, ww.d, ww.e, ww.f, ww.g, ww.h, 0)
        if  ww == specialW9j
            push!( w6js, W6j(ww.a, ww.b, ww.c, ww.e, ww.d, ww.g) )
            wa = RacahExpression( rex.summations, integrals, rex.phase + ww.b + ww.c + ww.d + ww.g,
                                    rex.weight / sqrt( (2*ww.c+1)*(2*ww.g+1)), deltas, triangles, w3js, w6js, w9js, ylms, djpqs )
            return( (true, wa) )
        end
    end
    
    return( (false, RacahExpression() ) )
end
