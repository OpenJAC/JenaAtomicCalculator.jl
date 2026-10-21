#
println("Sa) Apply & test the StarkBroadening module: the electron-impact Stark width of a line via")
println("    Dimitrijevic & Konjevic's Simplified Modified Semiempirical (SMSE) method.")

setDefaults("print summary: open", "zzz-StarkBroadening.sum")


# THE PHYSICS: free electrons in a plasma repeatedly perturb a radiating ion's upper and lower level,
# interrupting the phase of the emitted/absorbed radiation and broadening the line. SMSE
# (Dimitrijevic & Konjevic, Astron. Astrophys. 172, 345 (1987); formula as restated in Dimitrijevic &
# Popovic, J. Appl. Spectrosc. 68, 893 (2001), Eq. (7)) is the fully closed-form, low-temperature limit
# of their own earlier Modified Semiempirical method -- it needs only the line's own upper and lower
# level, via their effective quantum numbers n* and the orbital angular momentum l of the one active
# (valence) subshell, not an explicit perturber search. See StarkBroadening's own module docstring for
# what is and is not implemented (the shift, and the full perturber-by-perturber MSE, are NOT).


if  true
    # Last successful:  21-Oct-2026   ## C IV 2s-2p (1549 A) at N=1e17 cm^-3: computed W = 3.467e-3,
    #                                    2.452e-3, 1.733e-3, 1.415e-3 A at T=10000/20000/40000/60000 K,
    #                                    against Dimitrijevic & Konjevic (1980) Table 3's own W^SE column
    #                                    0.0038, 0.0027, 0.0019, 0.0016 A -- a CONSISTENT 88-91% across a
    #                                    6x range in T. See this branch's own comment for the residual.
    #
    # Branch a: Li-like C IV (Z=6), the textbook validation case -- the 1980 paper's own introduction
    #   names C III and C IV as "the only reliable [experimental Stark-width] data" available when the
    #   method was built, and its own Table 3 tabulates this EXACT line (2s^2S-2p^2P0, 1549 A) at an
    #   electron density that a garbled OCR exponent made read as 10^11 at first -- running at THAT value
    #   gave a width six orders of magnitude off the table, which is what revealed the true density is
    #   10^17 cm^-3, matching the round number used here.
    #
    # THE RESIDUAL, stated plainly rather than hidden: this branch's own computed width sits a CONSISTENT
    #   9-12% below the paper's own W^SE column at all four tabulated temperatures -- not a random
    #   scatter, and not a different order of magnitude. A genuinely wrong formula would not track a
    #   published table this tightly across a 6x range in T; a small, near-constant fractional offset is
    #   the signature of the SAME physics from a slightly different INPUT, most likely the C IV 2s
    #   ionization energy (64.49 eV, from NIST/standard tables here, not re-derived from the 1980 paper
    #   itself) feeding into n*. This residual is left open rather than tuned away quietly.
    #
    # NIST ASD: C IV (Li-like) 2s ionization energy = 64.49 eV. Dimitrijevic & Konjevic (1980), Table 3:
    #   C IV 2s^2S-2p^2P0 (mult. 1, UV, 1549.1 A), N = 1(10^17) cm^-3.
    ni   = Nuclear.Model(6.)
    grid = Radial.Grid(true)
    asfSettings = AsfSettings()

    wbU = SelfConsistent.performSCF([Configuration("1s^2 2p")], ni, grid, asfSettings; printout=false)
    wbL = SelfConsistent.performSCF([Configuration("1s^2 2s")], ni, grid, asfSettings; printout=false)
    upperMp = Hamiltonian.performCI(wbU.levels[1].basis, ni, grid, asfSettings; printout=false)
    lowerMp = Hamiltonian.performCI(wbL.levels[1].basis, ni, grid, asfSettings; printout=false)

    ev_to_au          = 1.0 / Defaults.convertUnits("energy: from atomic to eV", 1.0)
    ionizationEnergy  = lowerMp.levels[1].energy + 64.49*ev_to_au     ## NIST ASD: C IV 2s -> C V, 64.49 eV
    method            = Basics.SimplifiedModifiedSemiempiricalMethod(ionizationEnergy, 4.0)   ## Z=4 for C IV

    publishedWse = Dict(10000.0 => 0.0038, 20000.0 => 0.0027, 40000.0 => 0.0019, 60000.0 => 0.0016)   ## Angstrom

    cAU      = Defaults.getDefaults("speed of light: c")
    bohr_A   = Defaults.BOHR_RADIUS_SI * 1.0e10

    println("\n>>> C IV 2s-2p (1549 A), N = 1e17 cm^-3 -- against Dimitrijevic & Konjevic (1980) Table 3's own W^SE:")
    for  T  in  (10000.0, 20000.0, 40000.0, 60000.0)
        sbSettings = StarkBroadening.Settings(StarkBroadening.Settings(); multipoles=[E1], gauges=[UseCoulomb],
                          electronDensity=1.0e17, electronTemperature=T, method=method, printBefore=false)
        lines = perform( Atomic.Computation(Atomic.Computation(); name="C IV 2s-2p SMSE", grid=grid, nuclearModel=ni,
                          initialConfigs=[Configuration("1s^2 2p")], finalConfigs=[Configuration("1s^2 2s")],
                          processSettings=sbSettings); output=true )[StarkBroadening.Settings]

        for  ln  in  lines
            lambdaAU   = 2*pi*cAU/ln.omega
            dlambdaAU  = ln.electronWidth * lambdaAU^2/(2*pi*cAU)
            dlambda_A  = dlambdaAU * bohr_A
            @printf("    T=%6.0f K   u=%-4s   W = %.6f A   (published W^SE = %.4f A,  ratio = %.3f)\n",
                    T, string(ln.initialLevel.J), dlambda_A, publishedWse[T], dlambda_A/publishedWse[T])
        end
    end
    println("\n>>> source: Dimitrijevic & Konjevic, J. Quant. Spectrosc. Radiat. Transf. 24, 451 (1980), Table 3;")
    println(">>>         ionization energy from NIST ASD.")
    #
end
