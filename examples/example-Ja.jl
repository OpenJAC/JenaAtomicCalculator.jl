#
println("Ja) Apply & test the average-atom computations.")

if  true
    #
    # Last successful:  08-Oct-2026
    # Branch 1: Plasma.perform(Plasma.Computation(...; scheme=Plasma.AverageAtomScheme(...))) -- a finite-temperature
    #   self-consistent average-atom (DFS) computation, followed by form factors F(q).
    # System: boron (Z=5) at T = 10 eV, rho = 2.463 g/cm^3 (near solid density).
    #
    # THE CELL AND THE BOX ARE TWO DIFFERENT THINGS, and that is the whole shape of this model.  The CELL is the
    #   Wigner-Seitz sphere: it holds Z electrons, the potential is zero outside it, and neutrality is imposed on
    #   what is inside.  The BOX is the B-spline grid, and it only has to be long enough that a bound state can
    #   decay and a continuum state can oscillate -- a few times R^(WS).  Taking the box to BE the cell forces
    #   every state to zero at the wall, which removes about a quarter of the free electrons; taking the cell to be
    #   the box solves a cell of the wrong density.  Branch 4 below measures both against published numbers.
    # Checks:
    #   - Charge neutrality is exact: 2.00000 bound + 3.00000 continuum electrons inside the cell = Z = 5.
    #   - Only the K shell survives at solid density: 1s_1/2 at -5.80740 Ha = -158.03 eV with occupation 1.00000,
    #     while 2s (+1.74 eV) and 2p (+5.25 eV) lie above zero, i.e. are pressure-ionized.  So Z* = 3.0000 is
    #     "boron minus its K shell", and it is exact rather than approximate because f(eps_1s) = 1 to machine
    #     precision and nothing else is bound.
    #   - mu = +0.250241 Ha.
    #   - Form factors F(q) fall monotonically from F(1) = 3.7817 to F(10) = 0.4420 a.u. (q in a_o^-1).  THE
    #     FORM-FACTOR INTEGRAL STOPS AT THE CELL: a form factor belongs to ONE atom, and the orbitals now reach
    #     well past the cell, so integrating to the end of the grid adds the neighbours (it gives F(1) = 2.716).
    #   - The run reports that the subshell-list density accounts for 4.8732 of the 5 electrons inside the cell;
    #     the complete density is results["radial density"] and carries Z over the same range.
    #
    nm          = Nuclear.Model(5.0, 10.82)
    rho         = 2.463      # [g/cm^3]
    temp_au     = Defaults.convertUnits("energy: from eV to atomic", 10.0 * 1)
    temperature = Defaults.convertUnits("temperature: from atomic to Kelvin", temp_au)    # [K]
    radiusWS    = Plasma.determineWignerSeitzRadius(rho, nm) 
    grid        = Radial.generateGrid(Radial.Grid(false), boxSize = 3 * radiusWS)
    settings    = Plasma.Settings(temperature, rho, false)
    qValues     = [ q for q in 1:10 ]
    scheme      = Plasma.AverageAtomScheme(5, 2, Basics.AaDFSField(), false, true, false, Subshell[], Float64[], qValues )
    
    wa          = Plasma.Computation(Plasma.Computation(), scheme=scheme,
                                     nuclearModel=nm, grid=grid, settings=settings)
    @show wa
    wb          = perform(wa, output=true)
    #
elseif  false
    #
    # Last successful:  08-Oct-2026
    # Branch 2: the same boron system at T = 5 eV instead of 10 eV, to see how little the result moves.
    # System: boron (Z=5) at T = 5 eV, rho = 2.463 g/cm^3 (unchanged from Branch 1).
    # Checks:
    #   - Z* = 3.00000, exactly as at 10 eV:  at solid density the K shell is bound and everything else is not,
    #     and a factor of two in temperature does not change which.  mu rises from +0.250241 to +0.375550 Ha.
    #   - 1s_1/2 = -5.84040 Ha = -158.93 eV, i.e. 0.90 eV from its 10 eV value -- the K shell does not feel a
    #     few-eV change in the plasma.  The tabulated NEUTRAL, ISOLATED atom value is 188.0 eV [J. A. Bearden &
    #     A. F. Burr, Rev. Mod. Phys. 39, 125 (1967); LBNL X-Ray Data Booklet Table 1-1; cf.
    #     Empirical.bindingEnergy(5, Shell("1s"), data=PeriodicTable.XrayDataBooklet())].  The ~29 eV gap is the
    #     solid-density environment: the K electron sits in a cell of radius 2.27 a_o holding three free
    #     electrons, which screen it.  Closing that gap needs a LOWER DENSITY, not a lower temperature.
    #   - A PUZZLE RECORDED HERE UNTIL 08-Oct-2026 IS GONE.  This branch used to report Z* = 3.25 at 5 eV against
    #     2.86 at 10 eV -- ionization RISING as the plasma cools -- and left it as "a genuine, only partly
    #     understood finding".  It was neither: both numbers came from a cell 3.0 times too large in volume whose
    #     continuum was summed over about seven of the available states per kappa.  Both temperatures now give
    #     exactly 3.0000.
    #
    nm          = Nuclear.Model(5.0, 10.82)
    rho         = 2.463      # [g/cm^3]
    temp_au     = Defaults.convertUnits("energy: from eV to atomic", 5.0)
    temperature = Defaults.convertUnits("temperature: from atomic to Kelvin", temp_au)    # [K]
    radiusWS    = Plasma.determineWignerSeitzRadius(rho, nm)
    grid        = Radial.generateGrid(Radial.Grid(false), boxSize = 3 * radiusWS)
    settings    = Plasma.Settings(temperature, rho, false)
    scheme      = Plasma.AverageAtomScheme(5, 2, Basics.AaDFSField(), false, false, false, Subshell[], Float64[], Float64[] )

    wa          = Plasma.Computation(Plasma.Computation(), scheme=scheme,
                                     nuclearModel=nm, grid=grid, settings=settings)
    @show wa
    wb          = perform(wa, output=true)

    ## Reference (literature) 1s binding energy for the cold-limit comparison.
    b1sRef      = Defaults.convertUnits("energy: from atomic", Empirical.bindingEnergy(5, Shell("1s"), data=PeriodicTable.XrayDataBooklet()))
    println("\n  Reference (Bearden & Burr 1967, via XrayDataBooklet): neutral B 1s binding energy = $b1sRef eV")
    #
elseif  false
    #
    # Last successful:  08-Oct-2026
    # Branch 3: aluminium at SOLID density, where the Ne-like core (1s^2 2s^2 2p^6) is known to stay bound while
    #   the weakly bound M shell (3s, 3p) ionizes, giving a mean charge Z* close to 3 that changes only slowly
    #   with T from room temperature up to several eV [G. Massacrier, M. Boehme, J. Vorberger, F. Soubiran &
    #   B. Militzer, Phys. Rev. Research 3, 023026 (2021), arXiv:2105.01927 -- the precise Z*(T) values quoted
    #   there were not independently re-verified against the source tables here].
    # System: aluminium (Z=13) at rho = 2.70 g/cm^3 (solid density), T = 1, 5 and 10 eV.
    # Checks:
    #   - Z* = 3.00024, 3.00026 and 3.02594 at T = 1, 5 and 10 eV, with mu = +0.312878, +0.235222 and -0.010768
    #     Ha:  close to 3 and nearly flat in temperature, which is the published behaviour.  1-3 s per point.
    #   - At 10 eV the Ne-like core is bound -- 1s_1/2 -1488.93, 2s_1/2 -92.62, 2p_1/2 -55.15, 2p_3/2 -54.71 eV --
    #     while 3s (+1.08), 3p (+3.24) and 3d (+5.60 eV) all lie above zero.  So the three free electrons are the
    #     M shell, named rather than inferred.
    #   - Neutrality is exact: 9.97406 + 3.02594 = 13 at 10 eV.
    #
    # THIS BRANCH NEVER RAN BEFORE 08-Oct-2026 and was kept only as a record of two failures, both since repaired:
    #   at T = 1 eV the chemical-potential search diverged to ~-2e72 within two iterations and reported a
    #   nonsensical "mean charge = 13.0 = Z" (its Newton derivative had exp(2w) where exp(w) was needed, four
    #   orders too small); and at T = 5 eV it ran for minutes with no result (the O(N^2) DFS potential was rebuilt
    #   once per symmetry block instead of once per iteration).
    #
    for  Tev  in  [1.0, 5.0, 10.0]
        nm          = Nuclear.Model(13.0, 26.98)
        rho         = 2.70       # [g/cm^3]
        temp_au     = Defaults.convertUnits("energy: from eV to atomic", Tev)
        temperature = Defaults.convertUnits("temperature: from atomic to Kelvin", temp_au)    # [K]
        radiusWS    = Plasma.determineWignerSeitzRadius(rho, nm)
        grid        = Radial.generateGrid(Radial.Grid(false), boxSize = 3 * radiusWS)
        settings    = Plasma.Settings(temperature, rho, false)
        scheme      = Plasma.AverageAtomScheme(5, 2, Basics.AaDFSField(), false, false, false, Subshell[], Float64[], Float64[] )

        wa          = Plasma.Computation(Plasma.Computation(), scheme=scheme,
                                         nuclearModel=nm, grid=grid, settings=settings)
        println("\n\n>>> Aluminium at rho = $rho g/cm^3 and T = $Tev eV:")
        wb          = perform(wa, output=true)
    end
    #
elseif  false
    #
    # Last successful:  08-Oct-2026
    # Branch 4: THE PUBLISHED BENCHMARK -- a line-by-line reproduction of Table 1 of W. R. Johnson, C. Guet &
    #   G. F. Bertsch, J. Quant. Spectrosc. Radiat. Transfer 99, 327 (2006): aluminium at 0.27 g/cm^3 (one tenth
    #   of metallic density) and T = 5 eV.  This is the single best test this module has, because every number in
    #   it is published and none of it is fitted.
    # Checks, at lMax = 8 and a box of 3 R^(WS) -- ours against theirs:
    #   R^(WS)      6.4419        6.44
    #   mu         -0.38307      -0.3823      0.20 %
    #   N_bound    11.50892      11.5059      0.03 %
    #   N_free      1.49108       1.4941      0.20 %
    #   2s         -3.9973       -3.980       0.43 %
    #   2p         -2.6196/-2.6033   -2.610    0.37 / 0.26 %
    #   3s         -0.2614       -0.259       0.93 %
    #   3p         -0.0553/-0.0545   -0.054    2.4 / 0.9 %
    #   1s        -55.3180      -55.189       0.23 %   <-- THIS ONE IS RELATIVITY, NOT AN ERROR: their Eq. (1) is
    #                                                      the SCHRODINGER equation, which is also why their table
    #                                                      carries one 2p and one 3p where this one carries two.
    #   All of it is stable between boxes of 2, 3 and 4 R^(WS).
    #
    # TWO THINGS WERE FOUND BY THIS COMPARISON, both of them repaired in src/ on 08-Oct-2026 and both invisible
    #   without a published number to check against:
    #   - the FREE-ELECTRON IDENTITY.  The free box states, counted inside the cell, must reproduce the ideal gas
    #     there;  they now do to 0.07 % (0.83242 against 0.83181), where with the box equal to the cell they were
    #     off by a factor of 3.3.  That identity is what makes their n_0 = 1.3404 meaningful, and our
    #     SelfConsistent.freeElectronDensity at THEIR mu gives 1.3406.
    #   - the ENERGY ZERO.  The electrostatic potential vanishes outside a neutral cell on its own, but the
    #     exchange does not -- it follows the density, which tends to rho_0 and not to zero.  Cutting it at the
    #     wall left a step of 0.0972 Ha there, and mu, 3s and 3p all came out too deep by the SAME 0.09-0.10 Ha
    #     against an exchange at the background density of 0.089185 Ha.  The exchange is now measured from its
    #     background value, so eps = 0 means an electron at rest in the plasma -- the same zero the ideal-gas
    #     density is written in.
    #
    nm          = Nuclear.Model(13.0, 26.98)
    rho         = 0.27       # [g/cm^3], one tenth of metallic density
    temp_au     = Defaults.convertUnits("energy: from eV to atomic", 5.0)
    temperature = Defaults.convertUnits("temperature: from atomic to Kelvin", temp_au)    # [K]
    radiusWS    = Plasma.determineWignerSeitzRadius(rho, nm)
    grid        = Radial.generateGrid(Radial.Grid(false), boxSize = 3 * radiusWS)
    settings    = Plasma.Settings(temperature, rho, false)
    scheme      = Plasma.AverageAtomScheme(11, 8, Basics.AaDFSField(), false, false, false, Subshell[], Float64[], Float64[] )

    wa          = Plasma.Computation(Plasma.Computation(), scheme=scheme,
                                     nuclearModel=nm, grid=grid, settings=settings)
    println("\n\n>>> Johnson, Guet & Bertsch, JQSRT 99, 327 (2006), Table 1:  Al at 0.27 g/cm^3 and T = 5 eV.")
    println(">>> Their values:  R^(WS) 6.44, mu -0.3823, 1s -55.189, 2s -3.980, 2p -2.610, 3s -0.259, 3p -0.054,")
    println(">>>                N_bound 11.5059, N_free 1.4941, free-electron background in the cell 1.3404.")
    wb          = perform(wa, output=true)
    #
end
