using PyPlot: plot
## pyplot()

println("Af) Generate, normalize & test for continuum orbitals in a local potential.")

grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 2.0e-2, rbox = 20.0)
nuclearZ = 2.0;   noPoints = grid.NoPoints - 100

if  false
    # Last successful:  unknown 
    # Define a full nuclear + Dirac-Fock-Slater potential
    wa   = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(nuclearZ), 
                              configs=[Configuration("1s^2 2s"), Configuration("1s^2 2p")], 
                              asfSettings=AsfSettings(AsfSettings(), scField=Basics.DFSField())  )

    wb = perform(wa, output=true)

    # Extract the ground level to generate a proper potential
    multiplet = wb["multiplet:"]
    ## grid      = wb["grid:"]
    level     = multiplet.levels[1]
    basis     = multiplet.levels[1].basis  

    #
    println("\n\n>> Compute the nuclear and Hartree potentials.")
    w1  = JAC.Nuclear.pointNucleus(nuclearZ, grid)
    w2  = compute("radial potential: Dirac-Fock-Slater", grid, level)     ## compute("radial potential: Hartree",grid, level)
    ww = JAC.add(w1,w2)                                                   ## ww = Radial.Potential(" ", 0.95*ww.Zr, w2.grid)
    ## ww  = JAC.Nuclear.pointNucleus(nuclearZ, grid)                         ## compute("radial potential: Hartree",grid, level)
    ## plot(grid.r, ww.Zr)
    #
elseif  true
    # Last successful:  unknown 
    # Generate ...
    phases1 = Float64[];   phases2 = Float64[];   phases3 = Float64[];   phases4 = Float64[];   
    norms1  = Float64[];   norms2  = Float64[];   norms3  = Float64[];   norms4  = Float64[]
    # Compute and normalize continuum orbitals for a range of energies
    settings = Continuum.Settings(false, noPoints)
    sh       = Subshell("1000s_1/2")
    energies = [x  for  x = 1.001:2.0:50.0];   
    #
    @time for  energy in energies   
        println("\n*** energy = $energy")
        wavenb   = sqrt( 2energy + energy * Defaults.getDefaults("alpha")^2 );      wavelgth = 2pi / wavenb
        println(">> Generate continuum orbital for kappa = $(sh.kappa) ($sh), energy = $energy a.u. and asympt. wavelength = $wavelgth; hp = $(grid.hp)")
        ## setDefaults("method: continuum, spherical Bessel");    setDefaults("method: normalization, pure sine")
        ## wc1      = Continuum.generateOrbitalLocalPotential(energy, sh, ww, settings)
        ## push!(phases1, wc1[2]);   push!(norms1, abs(wc1[3]))
        #
        setDefaults("method: continuum, Galerkin");            setDefaults("method: normalization, pure sine")
        wc2      = Continuum.generateOrbitalLocalPotential(energy, sh, ww, settings)
        push!(phases2, wc2[2]);   push!(norms2, abs(wc2[3]))
        setDefaults("method: continuum, Galerkin");            setDefaults("method: normalization, pure Coulomb")
        wc3      = Continuum.generateOrbitalLocalPotential(energy, sh, ww, settings)
        push!(phases3, wc3[2]);   push!(norms3, abs(wc3[3]))
        #
        setDefaults("method: continuum, Galerkin");            setDefaults("method: normalization, Ong-Russek")
        wc4      = Continuum.generateOrbitalLocalPotential(energy, sh, ww, settings)
        push!(phases4, wc4[2]);   push!(norms4, abs(wc4[3]))
    end
    #
    ##x plot(energies, norms2);   plot(energies, 1.005*norms3);   plot(energies, norms4) ## ;   plot(energies, phases4)
    ##x plot(energies, phases2);   plot(energies, phases3);   plot(energies, phases4);   plot(energies, norms2)
    plot(energies, phases2);   plot(energies, 1.005*phases3);   plot(energies, 1.01*phases4)
    #
elseif true
    # Last successful:  unknown 
    # Compute and display single continuum orbital
    ##  "method: continuum, pure sine"  "method: continuum, asymptotic Coulomb"  "method: continuum, Galerkin"  "method: continuum, spherical Bessel"
    ##  "method: normalization, pure sine"    "method: normalization, pure Coulomb"  "method: normalization, Ong-Russek"
    @show Defaults.GBL_CONT_SOLUTION, Defaults.GBL_CONT_NORMALIZATION
    settings = Continuum.Settings(false, noPoints)
    sh       = Subshell("1000g_7/2")
    energies = [0.5];   
    for  energy in energies   
        wavenb   = sqrt( 2energy + energy * Defaults.getDefaults("alpha")^2 );      wavelgth = 2pi / wavenb
        println(">> Generate continuum orbital for kappa = $(sh.kappa) ($sh), energy = $energy a.u. and asympt. wavelength = $wavelgth; hp = $(grid.hp)")
        @show " "
        ## setDefaults("method: continuum, pure sine");           setDefaults("method: normalization, pure sine")
        ## @show Defaults.GBL_CONT_SOLUTION, Defaults.GBL_CONT_NORMALIZATION
        ## wc1      = Continuum.generateOrbitalLocalPotential(energy, sh, ww, settings)
        ## plot(grid.r[1:noPoints], wc1[1].P[1:noPoints])
        #
        ## setDefaults("method: continuum, spherical Bessel");    setDefaults("method: normalization, pure sine")
        ## @show " "
        ## @show Defaults.GBL_CONT_SOLUTION, Defaults.GBL_CONT_NORMALIZATION
        ## wc2      = Continuum.generateOrbitalLocalPotential(energy, sh, ww, settings)
        ## plot(grid.r[1:noPoints], wc2[1].P[1:noPoints])
        #
        setDefaults("method: continuum, Galerkin");            setDefaults("method: normalization, pure sine")
        @show " "
        @show Defaults.GBL_CONT_SOLUTION, Defaults.GBL_CONT_NORMALIZATION
        wc3      = Continuum.generateOrbitalLocalPotential(energy, sh, ww, settings)
        plot(grid.r[1:noPoints], wc3[1].P[1:noPoints])
        #
        setDefaults("method: continuum, Galerkin");            setDefaults("method: normalization, pure Coulomb")
        @show " "
        @show Defaults.GBL_CONT_SOLUTION, Defaults.GBL_CONT_NORMALIZATION
        wc4      = Continuum.generateOrbitalLocalPotential(energy, sh, ww, settings)
        plot(grid.r[1:noPoints], wc4[1].P[1:noPoints])
        #
        setDefaults("method: continuum, Galerkin");            setDefaults("method: normalization, Ong-Russek")
        @show " "
        @show Defaults.GBL_CONT_SOLUTION, Defaults.GBL_CONT_NORMALIZATION
        wc5      = Continuum.generateOrbitalLocalPotential(energy, sh, ww, settings)
        plot(grid.r[1:noPoints], wc5[1].P[1:noPoints])
    end
end
