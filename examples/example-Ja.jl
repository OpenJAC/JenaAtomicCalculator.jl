#
println("Ja) Apply & test the average-atom computations.")

if  true
    #
    # Last successful:  18-Jul-2026
    # Branch 1: Plasma.perform(Plasma.Computation(...; scheme=Plasma.AverageAtomScheme(...))) -- a finite-temperature
    #   self-consistent average-atom (DFS) computation, followed by form factors F(q).
    # System: boron (Z=5) at T = 10 eV, rho = 2.463 g/cm^3 (near solid density).
    # Checks:
    #   - Total occupation over all subshells = 5.009, matching Z = 5 to <0.2% (charge neutrality via the
    #     Newton-Raphson chemical-potential search).
    #   - 1s stays tightly bound (-168 eV, close to the isolated-atom value); 2s barely bound (~-1 eV); 2p pushed
    #     just positive (unbound/resonant) -- the expected ordering for boron at this temperature.
    #   - Mean charge state = 2.86, a plausible average ionization for boron at 10 eV.
    #   - Form factors F(q) decrease monotonically and smoothly from F(1)=3.73 to F(10)=0.43 a.u. (q in a_o^-1).
    #
    nm          = Nuclear.Model(5.0, 10.82)
    rho         = 2.463      # [g/cm^3]
    temp_au     = Defaults.convertUnits("energy: from eV to atomic", 10.0 * 1)
    temperature = Defaults.convertUnits("temperature: from atomic to Kelvin", temp_au)    # [K]
    radiusWS    = Plasma.determineWignerSeitzRadius(rho, nm) 
    ## gridx    = Radial.Grid(Radial.Grid(false), rnt = 2.0e-6, h = 5.0e-2, hp = 2.0e-2, rbox = 3.0)
    grid        = Radial.generateGrid(Radial.Grid(false), boxSize = radiusWS +1.0)
    settings    = Plasma.Settings(temperature, rho, false)
    qValues     = [ q for q in 1:10 ]
    scheme      = Plasma.AverageAtomScheme(5, 2, Basics.AaDFSField(), false, true, false, Subshell[], Float64[], qValues )
    
    wa          = Plasma.Computation(Plasma.Computation(), scheme=scheme,
                                     nuclearModel=nm, grid=grid, settings=settings)
    @show wa
    wb          = perform(wa, output=true)

    ## chemMu      = AverageAtom.determineChemicalPotential(orbitals, temp_au, nm::Nuclear.Model, grid)
    #
elseif  false
    #
    # Last successful:  18-Jul-2026
    # Branch 2: a lower-temperature point for the same boron system as Branch 1 -- T = 5 eV instead of 10 eV.
    #   This was originally meant to approach the cold, neutral-atom limit, but exploration this session found
    #   that the plain (undamped) fixed-point SCF iteration in SelfConsistent.solveAverageAtomField, and the
    #   Newton-Raphson chemical-potential search feeding it, become numerically unreliable well before reaching
    #   that limit at *solid* density: at T = 1 eV, chemMu no longer settles down and self-consistency stalls at
    #   ~6e-3 (vs the 1e-6 target); at T = 0.1 eV it oscillates wildly between iterations and never converges at
    #   all. T = 5 eV is the practical lower boundary where this implementation still converges cleanly.
    # System: boron (Z=5) at T = 5 eV, rho = 2.463 g/cm^3 (unchanged from Branch 1).
    # Checks:
    #   - Self-consistency reaches ~7e-6 by the 32nd iteration -- essentially converged, unlike T = 1 eV or below.
    #   - The 1s eigenvalue barely shifts from Branch 1 (-6.16 Hartree here vs -6.18 Hartree at T = 10 eV, i.e.
    #     ~168 eV in both cases) -- expected, since 1s should not feel a few-eV change in the plasma environment.
    #     The tabulated *neutral, isolated* atom value is a much larger 188.0 eV [J. A. Bearden & A. F. Burr,
    #     Rev. Mod. Phys. 39, 125 (1967); LBNL X-Ray Data Booklet Table 1-1, cross-checked against NIST/CRC data
    #     earlier in the Empirical-module work this session -- cf. Empirical.bindingEnergy(5, Shell("1s"),
    #     data=PeriodicTable.XrayDataBooklet())]; the ~20 eV gap is attributed to solid-density confinement
    #     (pressure-ionization/environment) effects, not something this branch expects to close, since reaching
    #     the truly isolated-atom limit would need reducing the *density* as well, not just the temperature (a
    #     separate, considerably more expensive test -- the Wigner-Seitz radius, and hence the required radial
    #     grid, grows quickly as the density drops, and the DFS potential's cost scales as its square).
    #   - Mean charge = 3.25, *higher* than Branch 1's 2.86 at the hotter T = 10 eV -- the opposite of naive
    #     thermal-ionization intuition. Tracing the occupation numbers: this is not driven by 2p (occupation is
    #     essentially unchanged, ~1.71 electrons, between the two temperatures) but by 2s, whose occupation rises
    #     markedly (0.87 -> 1.21) while remaining bound; the net movement in the "positive-energy" (formally
    #     ionized) population comes out larger at the lower temperature at this fixed density. This is reported
    #     as a genuine, only partly understood finding of the self-consistent occupation balance -- not yet
    #     explained away -- rather than a smooth, expected trend.
    #
    nm          = Nuclear.Model(5.0, 10.82)
    rho         = 2.463      # [g/cm^3]
    temp_au     = Defaults.convertUnits("energy: from eV to atomic", 5.0)
    temperature = Defaults.convertUnits("temperature: from atomic to Kelvin", temp_au)    # [K]
    radiusWS    = Plasma.determineWignerSeitzRadius(rho, nm)
    grid        = Radial.generateGrid(Radial.Grid(false), boxSize = radiusWS +1.0)
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
    # Last successful:  NOT WORKING -- kept only to record the diagnosis, see below.
    # Branch 3 (attempted): a qualitative literature comparison for aluminum at solid density, where the Ne-like
    #   core ([Ne] = 1s^2 2s^2 2p^6) is known to stay bound while the weakly-bound M-shell (3s, 3p) ionizes, giving
    #   a mean charge state Z* close to 3 that changes only slowly with T from room temperature up to several eV
    #   [G. Massacrier, M. Boehme, J. Vorberger, F. Soubiran & B. Militzer, "Reconciling ionization energies and
    #   band gaps of warm dense matter derived with ab initio simulations and average atom models", Phys. Rev.
    #   Research 3, 023026 (2021), arXiv:2105.01927 -- the precise Z*(T) values quoted there were not
    #   independently re-verified against the source tables here].
    # System (as originally planned): aluminum (Z=13) at rho = 2.70 g/cm^3 (solid density), T = 1 and 10 eV.
    # Diagnosis (this session): does NOT currently work at either planned temperature.
    #   - T = 1 eV: the initial Newton-Raphson chemical-potential search (run on the starting hydrogenic guess,
    #     before any SCF iteration) diverges catastrophically -- chemMu runs away to ~-2e72 within 2 iterations --
    #     collapsing every orbital occupation to exactly 0 and reporting a nonsensical "mean charge = 13.0 = Z"
    #     (fully ionized, zero density). Far worse than boron's T = 1 eV case (Branch 2's diagnosis), which at
    #     least stayed bounded even though it didn't converge.
    #   - T = 5 eV: still running after several minutes of active CPU work with no result -- abandoned rather than
    #     let it block further; a 19-subshell average-atom SCF should complete in seconds, not minutes, so this
    #     points to a genuine performance problem (most likely candidates: the uncontrolled/undamped Newton-Raphson
    #     iteration, and/or the O(N^2) grid-point double loop in the DFS potential, Basics.computePotential(::
    #     AaDFSField, ...)), separate from -- and probably compounding -- the convergence-robustness gap found for
    #     boron in Branch 2.
    #   - Aluminum (Z=13, more electrons, a shallower/more crowded valence structure than boron) appears to be a
    #     harder case for the current chemical-potential search than boron was, not an easier, more standard one.
    # Conclusion: this branch is left disabled and undeveloped. Revisiting it needs the underlying
    #   chemical-potential search and SCF iteration to be made more robust (and faster) first -- likely worth
    #   testing at a lower level (e.g. Empirical/unit-testing SelfConsistent.determineChemicalPotential and
    #   solveAverageAtomField directly against small, cheap synthetic cases) rather than only through full
    #   Plasma.Computation runs, before attempting more AverageAtomScheme example branches.
    #
end
