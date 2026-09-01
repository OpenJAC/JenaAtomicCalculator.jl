#
println("Ge) Symbolic simplification of spherical tensor operators, matrix elements and spherical amplitudes.")

using SymEngine
j1  = Basic(:j1);    j2  = Basic(:j2);    j3  = Basic(:j3);    j4  = Basic(:j4);    j5  = Basic(:j5)
m1  = Basic(:m1);    m2  = Basic(:m2);    m3  = Basic(:m3);    m4  = Basic(:m4);    m5  = Basic(:m5)
ja  = Basic(:ja);    jb  = Basic(:jb);    jc  = Basic(:jc);    jd  = Basic(:jd);    je  = Basic(:je)
ma  = Basic(:ma);    mb  = Basic(:mb);    mc  = Basic(:mc);    md  = Basic(:md);    me  = Basic(:me)
k   = Basic(:k);     q   = Basic(:q );    K   = Basic(:K);     Q   = Basic(:Q);
kk  = Basic(:kk);    qq  = Basic(:qq);    KK  = Basic(:KK);    QQ  = Basic(:QQ);

if  true
    # Last visit:      01-Sep-2026
    # Last successful: 01-Sep-2026 -- the reductions were read and are as described below: the field components
    #                  come out separated from the electronic reduced matrix elements, with the 3-j symbols
    #                  collected into one Racah expression.
    #
    # Branch a: BUILD THE OBJECTS AND REDUCE THEM. Operators (C^(k), the Coulomb and dipole operators, T^(k)), a
    #   scalar constant, a field vector u^(1), tensor products of increasing depth, spherical states, and finally
    #   the matrix elements -- which `expandSphericalMatrixElements` reduces by the Wigner-Eckart theorem, pulling
    #   the constants and the FIELD components out of the electronic reduced matrix elements and collecting the
    #   3-j symbols into one Racah expression for later simplification.
    #
    # WHAT A CORRECT RUN LOOKS LIKE (01-Sep-2026): the last two amplitudes print as a RacahExpression carrying a
    #   sum over magnetic quantum numbers with W3j symbols, a list of constants, the field components separated
    #   out, and the electronic matrix elements left reduced -- e.g. for [u^(1) x D^(1)]^(0) two W3j(1,1,0;...)
    #   factors with the field components u^(1*) and u^(1) pulled clear. That separation IS the point of the
    #   module: the field and the electron parts must come apart, and they do.
    # Define simple electronic and field operators
    Ck      = SphericalTensor.CkOperator(k);                 println("$Ck");         @show Ck
    Coulomb = SphericalTensor.CoulombOperator();             println("$Coulomb");    @show Coulomb
    Dop     = SphericalTensor.DipoleOperator(false);         println("$Dop");        @show Dop
    Tk      = SphericalTensor.TkOperator(k);                 println("$Tk");         @show Tk

    gConst  = SphericalTensor.ScalarConstant(false, Basic(:g));                 println("$gConst");        @show gConst
    uVector = SphericalTensor.UVector(false);                                   println("$uVector");       @show uVector

    # Define spherical tensor product
    prod1   = SphericalTensor.TensorProduct(false, j2, Ck, Tk);                 println("$prod1");         @show prod1
    prod2   = SphericalTensor.TensorProduct(false, Basic(0), uVector, Dop);     println("$prod2");         @show prod2
    prod3   = SphericalTensor.TensorProduct(false, K, prod1, prod2);            println("$prod3");         @show prod3
    prod4   = SphericalTensor.TensorProduct(false, K, uVector, Tk);             println("$prod4");         @show prod4

    # Define spherical and reduced states
    jma     = SphericalTensor.SphericalState(false, ja, ma);                    println("$jma");           @show jma
    jmb     = SphericalTensor.SphericalState(false, jb, mb);                    println("$jmb");           @show jmb
    j2      = SphericalTensor.ReducedState(j2);                                 println("$j2");            @show j2
    j3      = SphericalTensor.ReducedState(j3);                                 println("$j3");            @show j3

    # Define spherical matrix elements
    me1     = SphericalTensor.SphericalMatrixElement(jma, Tk, q, jmb);          println("$me1");           @show me1
    me2     = SphericalTensor.SphericalMatrixElement(jma, prod2, q, jmb);       println("$me2");           @show me2
    me3     = SphericalTensor.SphericalMatrixElement(jma, prod3, q, jmb);       println("$me3");           @show me3
    me4     = SphericalTensor.SphericalMatrixElement(false, jma, gConst, prod3, q, jmb);       
                                                                                println("$me4");           @show me4
    me2s    = SphericalTensor.SphericalMatrixElement(true,  jma, gConst, prod2, q, jmb);       
                                                                                println("$me2s");          @show me2s
    me5     = SphericalTensor.SphericalMatrixElement(false, jma, gConst, prod4, q, jmb);       
                                                                                println("$me5");           @show me5
    me5s    = SphericalTensor.SphericalMatrixElement(true,  jma, gConst, prod4, q, jmb);       
                                                                                println("$me5s");          @show me5s

    # Expand spherical operators
    we1     = SphericalTensor.expandSphericalTensorComponent(Coulomb, q);       println("$Coulomb");       @show we1
    we2     = SphericalTensor.expandSphericalTensorComponent(uVector, q);       println("$uVector");       @show we2
    we3     = SphericalTensor.expandSphericalTensorComponent(prod1, q);         println("$prod1");         @show we3
    we4     = SphericalTensor.expandSphericalTensorComponent(prod2, q);         println("$prod2");         @show we4
    we5     = SphericalTensor.expandSphericalTensorComponent(prod3, q);         println("$prod3");         @show we5

    # Expand spherical operators
    wm1     = SphericalTensor.expandSphericalMatrixElements([me1]);             println("\n $me1");        @show wm1
    wm2     = SphericalTensor.expandSphericalMatrixElements([me2]);             println("\n $me2");        @show wm2
    wm3     = SphericalTensor.expandSphericalMatrixElements([me3]);             println("\n $me3");        @show wm3
    wm4     = SphericalTensor.expandSphericalMatrixElements([me4]);             println("\n $me4");        @show wm4
    wm5     = SphericalTensor.expandSphericalMatrixElements([me2s, me2]);       println("\n $me2s $me2");  @show wm5
    wm6     = SphericalTensor.expandSphericalMatrixElements([me5s, me5]);       println("\n $me5s $me5");  @show wm6

elseif  false
    # Last visit:      01-Sep-2026
    # Last successful: 01-Sep-2026 -- all three identities hold. The Wigner-Eckart factor prints as
    #                  (-1)^(ja+jb-k) sqrt(2ja+1) W3j(jb, k, ja; mb, q, -ja), which IS the Clebsch-Gordan
    #                  <jb mb, k q | ja ma> in 3-j form; conjugating twice returns the identical object; and the
    #                  ranks come out 1 and K.
    #
    # Branch b: THREE IDENTITIES THAT THE SYMBOLIC MACHINERY MUST SATISFY, and none of them needs a reference.
    #   Branch a builds operators and reduces matrix elements; it shows that the machinery RUNS. This branch asks
    #   whether what it produces is RIGHT, using statements that are true by angular-momentum algebra alone.
    #
    # (1) THE WIGNER-ECKART FACTOR IS A CLEBSCH-GORDAN COEFFICIENT, and nothing else. The theorem separates
    #     <a j_a m_a | T^(k)_q | b j_b m_b> into a geometric factor and a reduced matrix element; the geometric
    #     factor is <j_b m_b, k q | j_a m_a>. Asking the module for it and printing it is the cheapest check that
    #     the reduction has not acquired an extra phase or weight.
    wa1 = SphericalTensor.getWignerEckardtFactor( SphericalTensor.SphericalState(false, ja, ma), k, q,
                                                  SphericalTensor.SphericalState(false, jb, mb) )
    println("\n (1) Wigner-Eckart factor <jb mb, k q | ja ma>:");   @show wa1

    # (2) CONJUGATION IS AN INVOLUTION: starring a tensor or a matrix element twice returns the original object.
    #     The `star` flag is carried through products and states, and a sign or a flag dropped on the way would
    #     show up here and nowhere else, since a single conjugation looks perfectly reasonable by itself.
    Tk2   = SphericalTensor.TkOperator(k)
    uV2   = SphericalTensor.UVector(false)
    prodA = SphericalTensor.TensorProduct(false, K, uV2, Tk2)
    twice = SphericalTensor.conjugate( SphericalTensor.conjugate(prodA) )
    println("\n (2) conjugation twice must return the original tensor product:")
    println("     original          : $prodA")
    println("     conjugated twice  : $twice")
    println("     identical         : $(string(prodA) == string(twice))")

    # (3) THE RANK OF A TENSOR PRODUCT IS THE RANK IT WAS COUPLED TO, not something derived from its parts. A
    #     product [A^(k1) x B^(k2)]^(K) has rank K by construction, and getRank must return exactly that -- the
    #     one place where a coupled rank could quietly be replaced by one of the factors' ranks.
    println("\n (3) ranks:")
    println("     getRank(u^(1))                 = $(SphericalTensor.getRank(uV2))       (must be 1)")
    println("     getRank([u^(1) x T^(k)]^(K))   = $(SphericalTensor.getRank(prodA))       (must be K)")

end
