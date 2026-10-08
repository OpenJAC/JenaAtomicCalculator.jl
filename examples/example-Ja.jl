#
println("Ja) Apply & test the average-atom computations.")

if  true
    #
    # Last successful:  08-Oct-2026
    # Branch 1: Plasma.perform(Plasma.Computation(...; scheme=Plasma.AverageAtomScheme(...))) -- a finite-temperature
    #   self-consistent average-atom (DFS) computation, followed by form factors F(q).
    # System: boron (Z=5) at T = 10 eV, rho = 2.463 g/cm^3 (near solid density).
    #
    # THE BOX IS THE WIGNER-SEITZ RADIUS, NOT R^(WS) + 1.  Every state of the model, bound or free, is normalised
    #   over the B-spline box, and neutrality is imposed on the electron count in THAT volume, so a box wider than
    #   the cell solves a cell of the wrong density -- here 3.0 times too large in volume.  The code says so now.
    #   Up to 08-Oct-2026 all three branches used radiusWS + 1.0, which is why their recorded numbers changed.
    # Checks:
    #   - Charge neutrality is exact: 2.00000 bound + 3.00000 continuum electrons = Z = 5.
    #   - Only the K shell survives at solid density: 1s_1/2 at -5.5751 Ha = -151.71 eV with occupation 1.00000,
    #     while 2s (+18.19 eV) and 2p (+23.79 eV) are pushed above zero, i.e. pressure-ionized.  Hence the mean
    #     charge is Z* = 3.0000 = "boron minus its K shell", and it is exact rather than approximate because
    #     f(eps_1s) = 1 to machine precision and nothing else is bound.
    #   - mu = +0.61882 Ha, i.e. the cell is degenerate at this density;  an ideal electron gas at the same mu,
    #     temperature and volume would hold 3.676 electrons against the cell's 3.014, the difference being the
    #     density the cell's Dirichlet wall removes.
    #   - Form factors F(q) fall monotonically from F(1) = 4.2104 to F(10) = 0.4401 a.u. (q in a_o^-1).
    #   - The run reports "Mean charge Z* = 3.00000 free and 2.00000 bound electrons", both from the full
    #     per-kappa spectrum.  Until 08-Oct-2026 it printed the BOUND count under the name "mean charge" and
    #     solved the chemical potential a second time from the subshell list;  both are fixed.
    #
    nm          = Nuclear.Model(5.0, 10.82)
    rho         = 2.463      # [g/cm^3]
    temp_au     = Defaults.convertUnits("energy: from eV to atomic", 10.0 * 1)
    temperature = Defaults.convertUnits("temperature: from atomic to Kelvin", temp_au)    # [K]
    radiusWS    = Plasma.determineWignerSeitzRadius(rho, nm) 
    grid        = Radial.generateGrid(Radial.Grid(false), boxSize = radiusWS)
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
    #   - Z* = 3.0000, exactly as at 10 eV:  at solid density the K shell is bound and everything else is not, and
    #     a factor of two in temperature does not change which.  mu rises from +0.6188 to +0.7251 Ha.
    #   - 1s_1/2 = -5.5661 Ha = -151.46 eV, i.e. 0.25 eV from its 10 eV value -- the K shell does not feel a
    #     few-eV change in the plasma.  The tabulated NEUTRAL, ISOLATED atom value is 188.0 eV [J. A. Bearden &
    #     A. F. Burr, Rev. Mod. Phys. 39, 125 (1967); LBNL X-Ray Data Booklet Table 1-1; cf.
    #     Empirical.bindingEnergy(5, Shell("1s"), data=PeriodicTable.XrayDataBooklet())].  The ~36 eV gap is the
    #     solid-density environment: the K electron sits in a cell of radius 2.27 a_o filled with three free
    #     electrons, which screen it.  Closing that gap needs a LOWER DENSITY, not a lower temperature.
    #   - A PUZZLE RECORDED HERE UNTIL 08-Oct-2026 IS GONE.  This branch used to report Z* = 3.25 at 5 eV against
    #     2.86 at 10 eV -- ionization RISING as the plasma cools -- and left it as "a genuine, only partly
    #     understood finding".  It was neither: both numbers came from a cell 3.0 times too large in volume whose
    #     continuum was summed over about seven of the available states per kappa.  On the Wigner-Seitz cell with
    #     the full spectrum both temperatures give exactly 3.0000.
    #
    nm          = Nuclear.Model(5.0, 10.82)
    rho         = 2.463      # [g/cm^3]
    temp_au     = Defaults.convertUnits("energy: from eV to atomic", 5.0)
    temperature = Defaults.convertUnits("temperature: from atomic to Kelvin", temp_au)    # [K]
    radiusWS    = Plasma.determineWignerSeitzRadius(rho, nm)
    grid        = Radial.generateGrid(Radial.Grid(false), boxSize = radiusWS)
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
    # Branch 3: THE LITERATURE COMPARISON -- aluminium at solid density, where the Ne-like core
    #   (1s^2 2s^2 2p^6) is known to stay bound while the weakly bound M shell (3s, 3p) ionizes, giving a mean
    #   charge Z* close to 3 that changes only slowly with T from room temperature up to several eV
    #   [G. Massacrier, M. Boehme, J. Vorberger, F. Soubiran & B. Militzer, Phys. Rev. Research 3, 023026 (2021),
    #   arXiv:2105.01927 -- the precise Z*(T) values quoted there were not independently re-verified against the
    #   source tables here].
    # System: aluminium (Z=13) at rho = 2.70 g/cm^3 (solid density), T = 1, 5 and 10 eV.
    # Checks:
    #   - Z* = 3.0000, 3.0000 and 3.0141 at T = 1, 5 and 10 eV:  close to 3 and nearly flat in temperature, which
    #     is the published behaviour.  Each point takes 1-3 s.
    #   - The Ne-like core is bound with occupation 1 at 10 eV -- 1s_1/2 at -1486.15 eV, 2s_1/2 at -91.25 eV,
    #     2p_1/2 at -53.74 eV, 2p_3/2 at -53.31 eV, f = 1.00000, 0.99995, 0.99776, 0.99766 -- while 3s (+10.57 eV),
    #     3p (+19.76 eV) and 3d (+31.10 eV) all lie above zero.  So the 3 free electrons are the M shell, named.
    #   - Neutrality is exact: 9.98588 + 3.01412 = 13 at 10 eV, 10.00000 + 3.00000 = 13 at 1 and 5 eV.
    #
    # THIS BRANCH NEVER RAN BEFORE 08-Oct-2026 and was kept only as a record of two failures, both since repaired:
    #   at T = 1 eV the chemical-potential search diverged to ~-2e72 within two iterations and reported a
    #   nonsensical "mean charge = 13.0 = Z" (its Newton derivative had exp(2w) where exp(w) was needed, four
    #   orders too small -- repaired 08-Oct-2026); and at T = 5 eV it ran for minutes with no result (the O(N^2)
    #   DFS potential was rebuilt once per symmetry block instead of once per iteration -- it is now built once).
    #
    for  Tev  in  [1.0, 5.0, 10.0]
        nm          = Nuclear.Model(13.0, 26.98)
        rho         = 2.70       # [g/cm^3]
        temp_au     = Defaults.convertUnits("energy: from eV to atomic", Tev)
        temperature = Defaults.convertUnits("temperature: from atomic to Kelvin", temp_au)    # [K]
        radiusWS    = Plasma.determineWignerSeitzRadius(rho, nm)
        grid        = Radial.generateGrid(Radial.Grid(false), boxSize = radiusWS)
        settings    = Plasma.Settings(temperature, rho, false)
        scheme      = Plasma.AverageAtomScheme(5, 2, Basics.AaDFSField(), false, false, false, Subshell[], Float64[], Float64[] )

        wa          = Plasma.Computation(Plasma.Computation(), scheme=scheme,
                                         nuclearModel=nm, grid=grid, settings=settings)
        println("\n\n>>> Aluminium at rho = $rho g/cm^3 and T = $Tev eV:")
        wb          = perform(wa, output=true)
    end
    #
end
