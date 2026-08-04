
println("Aa) Apply & test several radial e-e potentials.")

if  true
    # Last successful:  unknown ... need to be adapted
    # Compute different direct potentials for the charge density of a given level 
    grid      = JAC.Radial.Grid(true)
    multiplet = JAC.ManyElectron.Multiplet("from Ratip2012", "TestAuger/belike-ground-a-csl.inp",
                                                             "TestAuger/belike-ground-a-scf.out", "TestAuger/belike-ground-a-relci.mix") 
    level     = multiplet.levels[1]
    basis     = multiplet.levels[1].basis  
    #
    println("\nCompute the core-Hartree potential.")
    w1 = compute("radial potential: core-Hartree", grid, level)
    println(w1)
    #
    println("\nCompute the Hartree potential.")
    w2 = compute("radial potential: Hartree", grid, level)
    println(w2)
    #
    println("\nCompute the Hartree-Slater potential.")
    w3 = compute("radial potential: Hartree-Slater", grid, level)
    println(w3)
    # 
    println("\nCompute the Kohn-Sham potential.")
    w4 = compute("radial potential: Kohn-Sham", grid, level)
    println(w4)
    # 
    println("\nCompute the Dirac-Fock-Slater potential.")
    w5 = compute("radial potential: Dirac-Fock-Slater", grid, level)
    println(w5)
    error()
    # 
    println("\nCompute the extended-Hartree potential.")
    w6 = compute("radial potential: extended-Hartree", grid, level)
    println(w6)
    # 
    println("\nDisplay all four potentials together graphically.")
    JAC.plot("radial potential", [w1, w2, w3], grid; N = 260)
    #
end



