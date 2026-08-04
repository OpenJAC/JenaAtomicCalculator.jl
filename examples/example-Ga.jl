#
println("Ga) Symbolic evaluation by means of special values & recursion relations of the Wigner n-j symbols.")
#
using SymEngine

if      false   
    # Last successful:  13May2024
    # Use recursion relations for Wigner 3j symbols
    ja  = Basic(:ja);    jb = Basic(:jb);    jc = Basic(:jc);    ma = Basic(:ma);    mb = Basic(:mb);    mc = Basic(:mc)
    w3j = W3j(ja, jb, jc, ma, mb, mc)
    wa  = RacahAlgebra.recursionW3j(w3j, RacahAlgebra.RecursionW3jMagnetic())
    ## wa  = RacahAlgebra.recursionW3j(w3j, RacahAlgebra.RecursionW3jOneStep())
    ## wa  = RacahAlgebra.recursionW3j(w3j, RacahAlgebra.RecursionW3jHalfStep())
    ## wa  = RacahAlgebra.recursionW3j(w3j, RacahAlgebra.RecursionW3jLouck())
    #
elseif  false 
    # Last successful:  13May2024
    # Look for special value for Wigner 3j symbol 
    j   = Basic(:j);     m = Basic(:m)
    w3j = W3j(j+3//2, j, 3//2, m, -m-3//2, 3//2)
    rex = RacahAlgebra.evaluate(w3j);           @show rex
    rex = RacahAlgebra.equivalentForm(w3j);     @show rex
    RacahAlgebra.evaluate(rex, special=true)
    #
elseif  true    
    # Last successful:  13May2024
    # Look for special value for Wigner 6j symbol 
    a   = Basic(:a);     b  = Basic(:b);    c  = Basic(:c)
    w6j = W6j( a, b, c, 2, c-2, b-2)
    rex = RacahAlgebra.evaluate(w6j);           @show rex
    rex = RacahAlgebra.equivalentForm(w6j);     @show rex
    RacahAlgebra.evaluate(rex, special=true)
    #
elseif  true    
    # Last successful:  unknown ...
    # Look for special value for Wigner 9j symbol 
    a   = Basic(:a);    b = Basic(:b);    c = Basic(:c);    d = Basic(:d);    ee = Basic(:ee);    f = Basic(:f)
    w9j = W9j(a, b, ee, c, d, ee, f, f, Basic(0))
    rex = RacahAlgebra.evaluate(w9j);           @show rex
    rex = RacahAlgebra.equivalentForm(w9j);     @show rex
    RacahAlgebra.evaluate(rex, special=true)
    #
elseif  false    
    # Last successful:  unknown ...
    # Evaluate the recoupling coefficient for 3 angular momenta explicitly from Racah expression; from Racah-02.ps
    j1   = Basic(:j1);    j2 = Basic(:j2);    j3 = Basic(:j3);    J12 = Basic(:J12);    J23 = Basic(:J23);    J = Basic(:J)
    m1   = Basic(:m1);    m2 = Basic(:m2);    m3 = Basic(:m3);    M12 = Basic(:M12);    M23 = Basic(:M23);    M = Basic(:M)
    w3ja = W3j(J12, j3, J, M12, m3, -M);        w3jb = W3j(j1, j2, J12, m1, m2, -M12)       
    w3jc = W3j(j2, j3, J23, m2, m3, -M23);      w3jd = W3j(j1, J23, J, m1, M23, -M)   
    rex  = RacahExpression( [m1, m2, m3, M12, M23], -J12 + 2*j3 - 2*M - 2*j1 - M12 - M23 + J23, 
                            (2*J+1) * sqrt( (2*J12+1)*(2*J23+1) ), Kronecker[], Triangle[], [w3ja, w3jb, w3jc, w3jd], W6j[], W9j[] )
    rex  = RacahAlgebra.evaluate(rex);   @show rex
    rex  = RacahAlgebra.evaluate(rex)
    #
elseif  false   
    # Last successful:  unknown ...
    # Evaluate the Racah expression for the Coulomb interaction strength; from Racah-02.ps
    ja   = Basic(:ja);    jb = Basic(:jb);    jc = Basic(:jc);    jd = Basic(:jd);    Jab = Basic(:Jab);    Jcd = Basic(:Jcd);    L = Basic(:L)
    ma   = Basic(:ma);    mb = Basic(:mb);    mc = Basic(:mc);    md = Basic(:md);    Mab = Basic(:Mab);    Mcd = Basic(:Mcd);    M = Basic(:M)
    rexa = RacahAlgebra.ClebschGordan(ja, ma, jb, mb, Jab, Mab);    rexb = RacahAlgebra.ClebschGordan(jc, mc, jd, md, Jcd, Mcd)
    w3ja = W3j(ja, L, jc, -ma, M, mc);                              w3jb = W3j(jb, L, jd, -mb, -M, md)
    rex  = rexa * rexb * RacahExpression( [ma, mb, mc, md, M], L - M + ja - ma + jb - mb, Basic(1), Kronecker[], Triangle[], [w3ja, w3jb], W6j[], W9j[] )
    rex  = RacahAlgebra.evaluate(rex);   @show rex
    rex  = RacahAlgebra.evaluate(rex)
    #
elseif  true   
    # Last successful:  unknown ...
    # Evaluate the Racah expression for the Feynman-Goldstone diagram; from Racah-02.ps
    ja   = Basic(:ja);    jb = Basic(:jb);    jr = Basic(:jr);    js = Basic(:js);    L1 = Basic(:L1);    L2 = Basic(:L2)
    ma   = Basic(:ma);    mb = Basic(:mb);    mr = Basic(:mr);    ms = Basic(:ms);    M1 = Basic(:M1);    M2 = Basic(:M2)
    w3ja = W3j(jr, L1, ja, -mr, M1, ma);      w3jb = W3j(js, L1, jb, -ms, -M1, mb)
    w3jc = W3j(ja, L2, jr, -ma, M2, mr);      w3jd = W3j(jb, L2, js, -mb, -M2, ms)
    rex  = RacahExpression( [ma, mb, mr, ms, M1, M2], L1 - M1 + L2 - M2 + ja - ma + jb - mb + jr - mr + js - ms, 
                            Basic(1), Kronecker[], Triangle[], [w3ja, w3jb, w3jc, w3jd], W6j[], W9j[] )
    rex  = RacahAlgebra.evaluate(rex);   @show rex
    rex  = RacahAlgebra.evaluate(rex)
    #
elseif  false
    # Last successful:  unknown ...
    # Test for special values of Wigner 3j symbols; selectW3j(1..16)
    w3j = RacahAlgebra.selectW3j(16);                   println("w3j-original        = $w3j")
    rex = RacahAlgebra.equivalentForm(w3j);             println("rex-equivalent      = $rex")
    rex = RacahAlgebra.evaluate(rex, special=true);     println("rex-special value   = $rex")
    #
elseif  false
    # Last successful:  unknown ...
    # Test for special values of the Wigner 6j symbols; selectW6j(1..22)
    w6j = RacahAlgebra.selectW6j(22);                   println("w6j-original        = $w6j")
    rex = RacahAlgebra.equivalentForm(w6j);             println("rex-equivalent      = $rex")
    rex = RacahAlgebra.evaluate(rex, special=true);     println("rex-special value   = $rex")
    ## ws  = [:ja => 2//1, :jb => 2, :jc => 1]
    ## wb  = RacahAlgebra.subs(w6j, ws);                   println("\nw6j-substitution  = $wb")
    ## wc  = RacahAlgebra.subs(wa, ws);                    println("\nrhs-substitution  = $wc")
    ## wd  = RacahAlgebra.evaluateNumerical(wb);           println("\nlhs-numerical     = $wd")
    ## we  = RacahAlgebra.evaluateNumerical(wc);           println("\nrhs-numerical     = $we")
    #
elseif  true
    # Last successful:  unknown ...
    # Test for special values of the Wigner 9-j symbols; selectW6j(1..2)
    w9j = RacahAlgebra.selectW9j(2);                    println("w9j-original        = $w9j")
    rex = RacahAlgebra.equivalentForm(w9j);             println("rex-equivalent      = $rex")
    rex = RacahAlgebra.evaluate(rex, special=true);     println("rex-special value   = $rex")
    #
end
