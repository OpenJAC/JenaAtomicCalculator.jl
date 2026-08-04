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
    # Last successful:  unknown ...
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
    # Last successful:  unknown ...

end
