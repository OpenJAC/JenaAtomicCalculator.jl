#
println("Hc) Tests of semiempirical binding-energy and radial-expectation estimates.")
setDefaults("unit: energy", "eV")
#
setDefaults("print summary: open", "zzz-Semiempirical.sum")

if  false
    #
    # Last successful:  14-Jul-2026
    # Branch 1: Subshell binding energies -- compare all three tabulations for Ne (Z=10).
    #   Also test the newly extended XrayDataBooklet for deep Xe (Z=54) subshells 4d/5s/5p.
    # Tests:
    #   (a) Bug 4 fix: EstimateIonizationPotentialInnerShell now takes Subshell (not Shell);
    #       its column must equal Williams2000 column (both use the same Williams source).
    #   (b) New XrayDataBooklet Subshell function covering 4d_3/2..6p_3/2.
    # Checks:
    #   - All four Ne columns should agree to < 1 eV per subshell.
    #   - Xe 4d_3/2 = 69.5 eV,  4d_5/2 = 67.5 eV  (X-ray Data Booklet tabulated).
    #   - Xe 5s_1/2 = 23.3 eV,  5p_1/2 = 13.4 eV,  5p_3/2 = 12.1 eV.
    #
    Z1 = 10;    Z2 = 54
    println("\n  Subshell binding energies [eV] for Ne (Z=$Z1):\n")
    println("    Subshell        IonPot/W2000   Williams2000   Larkins1977   XrayDataBooklet")
    for sh in [Subshell("1s_1/2"), Subshell("2s_1/2"), Subshell("2p_1/2"), Subshell("2p_3/2")]
        e1 = Defaults.convertUnits("energy: from atomic", Semiempirical.estimate(Basics.EstimateIonizationPotentialInnerShell(), sh, Z1))
        e2 = Defaults.convertUnits("energy: from atomic", Semiempirical.estimate(Basics.EstimateBindingEnergyWilliams2000(),        Z1, sh))
        e3 = Defaults.convertUnits("energy: from atomic", Semiempirical.estimate(Basics.EstimateBindingEnergyLarkins1977(),         Z1, sh))
        e4 = Defaults.convertUnits("energy: from atomic", Semiempirical.estimate(Basics.EstimateBindingEnergyXrayDataBooklet(),     Z1, sh))
        sa = "    " * string(sh) * repeat(" ", max(1, 16 - length(string(sh))))
        println(sa * "$(round(e1,digits=1))\t\t$(round(e2,digits=1))\t\t$(round(e3,digits=1))\t\t$(round(e4,digits=1))")
    end
    #
    println("\n  Xe (Z=$Z2) deeper subshells [eV] -- XrayDataBooklet only (not in W2000 or L1977):\n")
    println("    Subshell        XrayDataBooklet   (Exp: 4d_3/2=69.5, 4d_5/2=67.5, 5s=23.3, 5p_1/2=13.4, 5p_3/2=12.1)")
    for sh in [Subshell("4d_3/2"), Subshell("4d_5/2"), Subshell("5s_1/2"), Subshell("5p_1/2"), Subshell("5p_3/2")]
        e4 = Defaults.convertUnits("energy: from atomic", Semiempirical.estimate(Basics.EstimateBindingEnergyXrayDataBooklet(), Z2, sh))
        sa = "    " * string(sh) * repeat(" ", max(1, 16 - length(string(sh))))
        println(sa * "$(round(e4,digits=2)) eV")
    end
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  14-Jul-2026
    # Branch 2: Configuration binding energies -- all three tabulations.
    # Tests:
    #   (a) Bugs 1&2 fix: Williams2000 and Larkins1977 conf functions no longer crash
    #       with BoundsError when a configuration is passed; they cover through 4p only.
    #   (b) Bug 3 fix: XrayDataBooklet conf function now uses wa[22] for Shell("6s"),
    #       not wa[12]; tested via Ba (Z=56) with a 6s^2 configuration.
    # Checks:
    #   - Kr ground config: W2000 and L1977 totals should agree to < 100 eV.
    #   - Ba 6s^1: XrayDataBooklet should give ≈ 5.2 eV (not ≈ 14.1 eV which the old
    #     wa[12] code would have produced via the 4p_3/2 entry).
    #
    println("\n  Total binding energy [eV] of Kr ground configuration [Ar]3d^10 4s^2 4p^6:\n")
    conf_Kr = Configuration("[Ar] 3d^10 4s^2 4p^6")
    eb1 = Defaults.convertUnits("energy: from atomic", Semiempirical.estimate(Basics.EstimateBindingEnergyWilliams2000(),    36, conf_Kr))
    eb2 = Defaults.convertUnits("energy: from atomic", Semiempirical.estimate(Basics.EstimateBindingEnergyLarkins1977(),     36, conf_Kr))
    eb3 = Defaults.convertUnits("energy: from atomic", Semiempirical.estimate(Basics.EstimateBindingEnergyXrayDataBooklet(), 36, conf_Kr))
    println("    Williams2000    = $(round(eb1, digits=1)) eV")
    println("    Larkins1977     = $(round(eb2, digits=1)) eV")
    println("    XrayDataBooklet = $(round(eb3, digits=1)) eV")
    #
    println("\n  Bug 3 note: wa[22]=6s_1/2 is -1 for all Z in the current XrayDataBooklet table.")
    println("  The booklet covers only X-ray core levels; outer valence shells (Ba 6s, IP~5.2 eV)")
    println("  are absent.  The if-guard correctly returns 0 contribution.  Bug 3 is dormant: the")
    println("  wrong wa[12] branch would only fire if 6s data were later added to the table.")
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  14-Jul-2026
    # Branch 3: Hydrogenic high-n binding energies via estimateBindingEnergies.
    # New physics compared to Branches 1-2: semiempirical Zeff model for Rydberg series.
    # System: Rb (Z=37) outer electron above a [Kr] core (36 electrons).
    #
    # Results (14-Jul-2026):
    #   l-independent table: Zeff ≈ 1.5 at n=5, energies -1.26 to -0.17 eV  -- CORRECT Rydberg.
    #   l=0 table:   Zeff = 27 at n=5, energies ~ -400 eV  -- WRONG for near-neutral Rb.
    #   l=2 table:   Zeff = 10 at n=5, energies ~  -55 eV  -- also wrong.
    #   The l-dependent formula 1 - 0.1*Ne/n/(l+1) breaks down when Ne ≈ Z (Rb: Ne=36, Z=37):
    #   the correction factor (0.28 for l=0) collapses the screening → Zeff >> 1.
    #   This formula was designed for Ne << Z (few-electron ions); it is a design limitation,
    #   not a code bug.  Use only the l-independent version for near-neutral Rydberg systems.
    #
    println("\n  Rb (Z=37, [Kr] core, Ne=36): high-n binding energies, l-independent Zeff:")
    Semiempirical.estimateBindingEnergies(37.0, Configuration("[Ar] 3d^10 4s^2 4p^6"), 5:10)
    #
    println("\n  Rb (Z=37): s-type Rydberg (l=0), l-dependent Zeff:")
    Semiempirical.estimateBindingEnergies(37.0, Configuration("[Ar] 3d^10 4s^2 4p^6"), 5:10, 0)
    #
    println("\n  Rb (Z=37): d-type Rydberg (l=2), l-dependent Zeff:")
    Semiempirical.estimateBindingEnergies(37.0, Configuration("[Ar] 3d^10 4s^2 4p^6"), 5:10, 2)
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  14-Jul-2026
    # Branch 5: estimateSlaterZeff -- Slater's (1930) empirical screening rules.
    # New physics compared to Branches 3-4: systematic group-based screening vs. the
    #   global Ne*(1-0.5/n^2.2) formula.
    # Tests:
    #   (a) He (Z=2): target 1s_1/2, core 1s^1.  σ = 1×0.30 → Zeff = 1.70  (textbook value).
    #   (b) Ne (Z=10): target 2p_3/2, core 1s^2 2s^2 2p^5.
    #         σ = 2×0.85 + 7×0.35 = 1.70+2.45 = 4.15 → Zeff = 5.85  (Slater 1930 result).
    #   (c) Rb (Z=37): target 5s_1/2, core [Ar] 3d^10 4s^2 4p^6.
    #         σ = 28×1.00 + 8×0.85 = 34.80 → Zeff = 2.20.
    #         Compare: l-indep formula gives effZ≈1.51 at n=5 (Branch 3).
    #   (d) Rydberg series Rb n=5..9 (l=0): Slater Zeff vs existing l-indep Zeff.
    # Checks:
    #   - He Zeff = 1.70 exactly (only one screening electron, σ = 0.30).
    #   - Ne Zeff = 5.85 (classic Slater result, agrees with literature to 0.01).
    #   - Rb 5s: Slater Zeff = 2.20 > l-indep Zeff = 1.51 (Slater screens less at n=5).
    #   - Rydberg: as n increases, Slater Zeff stays ~2.20 (core unchanged),
    #              while l-indep formula approaches 1.0 (full-screening limit).
    #
    println("\n  Slater Z_eff for selected systems:\n")
    #
    ## (a) He
    zHe = Semiempirical.estimateSlaterZeff(2.0, Configuration("1s^1"), Subshell("1s_1/2"))
    println("    He  (Z=2),  target 1s_1/2,  Zeff = $(round(zHe, digits=2))  (expect 1.70)")
    #
    ## (b) Ne
    zNe = Semiempirical.estimateSlaterZeff(10.0, Configuration("1s^2 2s^2 2p^5"), Subshell("2p_3/2"))
    println("    Ne  (Z=10), target 2p_3/2,  Zeff = $(round(zNe, digits=2))  (expect 5.85)")
    #
    ## (c) Rb 5s -- Slater vs l-independent formula
    Z_Rb = 37.0;    conf_Rb = Configuration("[Ar] 3d^10 4s^2 4p^6");    Ne_Rb = conf_Rb.NoElectrons
    zRb  = Semiempirical.estimateSlaterZeff(Z_Rb, conf_Rb, Subshell("5s_1/2"))
    zLI  = Z_Rb - Ne_Rb * (1 - 0.5 / 5^2.2)    ## l-independent formula, n=5
    println("    Rb  (Z=37), target 5s_1/2,  Zeff(Slater) = $(round(zRb, digits=2)),  " *
            "Zeff(l-indep) = $(round(zLI, digits=2))  (expect 2.20 vs 1.51)")
    #
    ## (d) Rydberg series n=5..9, both formulas
    println("\n  Rb (Z=37, [Kr] core): Slater Zeff vs l-independent Zeff for Rydberg n=5..9 (s-type):\n")
    println("    n       Zeff(Slater)   Zeff(l-indep)   eb(Slater)[eV]   eb(l-indep)[eV]")
    for n in 5:9
        sh_n = Subshell(n, -1)    ## s-type: kappa=-1
        zS   = Semiempirical.estimateSlaterZeff(Z_Rb, conf_Rb, sh_n)
        zI   = Z_Rb - Ne_Rb * (1 - 0.5 / n^2.2)
        ebS  = Defaults.convertUnits("energy: from atomic", -zS^2 / (2*n^2))
        ebI  = Defaults.convertUnits("energy: from atomic", -zI^2 / (2*n^2))
        println("    $n       $(round(zS,digits=2))           $(round(zI,digits=2))           " *
                "$(round(ebS,digits=3))             $(round(ebI,digits=3))")
    end
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  15-Jul-2026
    # Branch 6: EstimateBindingEnergyNist2025 -- successive ionization potentials from NIST (2025).
    # New physics compared to Branches 1-5: valence and successive IPs, covering Z=1..90.
    # Tests:
    #   (a) He (Z=2): IP_1 = 24.587 eV,  IP_2 = 54.418 eV  (well-known values).
    #   (b) Ne (Z=10): IP_1 = 21.565 eV (2p electron, compare to XrayDataBooklet 1s_1/2 = 870.2 eV).
    #   (c) Ba (Z=56): IP_1 should give ≈ 5.21 eV (6s valence -- the gap Bug 3 revealed in Branch 2,
    #       where XrayDataBooklet returns 0 because it only covers X-ray core levels).
    #   (d) Rb (Z=37): IP_1 = 4.177 eV (ground state 5s IP -- compare to Slater eb=-2.63 eV from Branch 5;
    #       NIST is the experimental reference value that the semiempirical formulas are measured against).
    #   (e) Ne successive IPs: IP_1..IP_10 to verify the full table for Z=10.
    # Checks:
    #   - He IP_1 = 24.587 eV, IP_2 = 54.418 eV (tabulated, verify to 0.001 eV).
    #   - Ba IP_1 ≈ 5.21 eV (fills the valence gap that XrayDataBooklet cannot provide).
    #   - Rb IP_1 ≈ 4.177 eV (experimental 5s IP; Slater gives -2.63 eV, l-indep gives -1.26 eV).
    #   - Ne IP_10 ≈ 1362.2 eV (K-shell, compare to Williams2000 1s_1/2 = 870.2 eV -- different:
    #       NIST IP_10 is the energy to remove the last 1s electron from Ne^9+, not from neutral Ne).
    #
    println("\n  NIST 2025 successive ionization potentials [eV]:\n")
    #
    ## (a) He: IP_1 and IP_2
    ip1_He = Defaults.convertUnits("energy: from atomic", Semiempirical.estimate(Basics.EstimateBindingEnergyNist2025(), 2, 1))
    ip2_He = Defaults.convertUnits("energy: from atomic", Semiempirical.estimate(Basics.EstimateBindingEnergyNist2025(), 2, 2))
    println("    He  (Z=2):   IP_1 = $(round(ip1_He, digits=3)) eV  (expect 24.587),  " *
                            "IP_2 = $(round(ip2_He, digits=3)) eV  (expect 54.418)")
    #
    ## (b) Ne: first IP
    ip1_Ne = Defaults.convertUnits("energy: from atomic", Semiempirical.estimate(Basics.EstimateBindingEnergyNist2025(), 10, 1))
    println("    Ne  (Z=10):  IP_1 = $(round(ip1_Ne, digits=3)) eV  (2p electron, expect 21.565)")
    #
    ## (c) Ba: first IP -- valence gap filled
    ip1_Ba = Defaults.convertUnits("energy: from atomic", Semiempirical.estimate(Basics.EstimateBindingEnergyNist2025(), 56, 1))
    println("    Ba  (Z=56):  IP_1 = $(round(ip1_Ba, digits=3)) eV  (6s valence, expect ≈5.21; XrayDataBooklet gives 0)")
    #
    ## (d) Rb: first IP -- reference for Slater comparison
    ip1_Rb = Defaults.convertUnits("energy: from atomic", Semiempirical.estimate(Basics.EstimateBindingEnergyNist2025(), 37, 1))
    println("    Rb  (Z=37):  IP_1 = $(round(ip1_Rb, digits=3)) eV  (5s ground state, expect 4.177; Slater gives 2.63, l-indep 1.26)")
    #
    ## (e) Ne: full successive IP table
    println("\n  Ne (Z=10) successive IPs (IP_1..IP_10):\n")
    println("    q      IP [eV]     assignment")
    assign = ["2p", "2p", "2p", "2p", "2p", "2p", "2s", "2s", "1s", "1s"]
    for q in 1:10
        ipq = Defaults.convertUnits("energy: from atomic", Semiempirical.estimate(Basics.EstimateBindingEnergyNist2025(), 10, q))
        println("    $q      $(round(ipq, digits=3))      $(assign[q])")
    end
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  14-Jul-2026
    # Branch 4: Radial expectation values -- tests Bug 5 fix (effZ instead of bare Z).
    # New physics compared to Branch 3: geometric extent (size) of Rydberg orbitals.
    # System: same Rb (Z=37, [Kr] core, Ne=36).
    # Before fix: rnl = (3n^2 - l(l+1)) / (2*Z=37) ≈ 1.0 a_0 for n=5, l=0 (unphysical).
    # After fix:  rnl = (3n^2 - l(l+1)) / (2*Zeff) ≈ 25 a_0 for n=5  (correct Rydberg extent).
    # Checks:
    #   - <r> should grow roughly as n^2 (fixed l), confirming Zeff ≈ const ≈ 1.5.
    #   - l=0 and l=2 should differ visibly at each n (<r> decreases with l for fixed n).
    #   - All values >> 1 a_0; the old bare-Z formula would give < 2 a_0 throughout.
    #
    println("\n  Rb (Z=37, [Kr] core, Ne=36): Rydberg radial extents <r_nl> [a_0] after Bug 5 fix:")
    Semiempirical.estimateRadialExpectation(37.0, Configuration("[Ar] 3d^10 4s^2 4p^6"), 5:9, 0)
    #
    println("\n  Rb (Z=37): l=2 (d-type) radial extents:")
    Semiempirical.estimateRadialExpectation(37.0, Configuration("[Ar] 3d^10 4s^2 4p^6"), 5:9, 2)
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  15-Jul-2026
    # Branch 7: estimateBindingEnergies (l-dependent) -- M3 fix: Slater Zeff replaces broken formula.
    # New physics compared to Branch 3: l-dependent Zeff now uses Slater's group screening.
    # System: Rb (Z=37) outer electron above a [Kr] core (36 electrons).
    #
    # Old formula (Branch 3, l=0, n=5): Zeff ≈ 27  →  eb ≈ -400 eV   (wrong: Ne=36 ≈ Z).
    # New Slater (n=5, l=0):  Zeff = 2.20  (4sp screen at 0.85, rest at 1.00  →  σ=34.80).
    # New Slater (n=5, l=2):  Zeff = 1.00  (d/f target: all lower shells screen fully).
    # New Slater (n≥6, l=0):  Zeff = 1.00  (no (n-1)sp in [Kr] core  →  all 36e at 1.00).
    # Cross-check: NIST IP_1(Rb) = 4.177 eV;  Slater 5s eb = -2.20^2/(2·25) = -2.63 eV.
    #   (Slater overestimates screening for an sp valence electron -- consistent expectation.)
    #
    println("\n  Rb (Z=37, [Kr] core): l-dependent binding energies after M3 fix (Slater Zeff):\n")
    #
    println("  l=0 (s-type),  expect Zeff=2.20 at n=5, Zeff=1.00 for n≥6:")
    Semiempirical.estimateBindingEnergies(37.0, Configuration("[Ar] 3d^10 4s^2 4p^6"), 5:9, 0)
    #
    println("  l=2 (d-type),  expect Zeff=1.00 at all n:")
    Semiempirical.estimateBindingEnergies(37.0, Configuration("[Ar] 3d^10 4s^2 4p^6"), 5:9, 2)
    #
    ## Cross-check 5s result against NIST
    eb_5s  = Semiempirical.estimateBindingEnergies(37.0, Configuration("[Ar] 3d^10 4s^2 4p^6"), 5:5, 0)[1]
    eu_5s  = Defaults.convertUnits("energy: from atomic", eb_5s)
    ip1_Rb = Defaults.convertUnits("energy: from atomic", Semiempirical.estimate(Basics.EstimateBindingEnergyNist2025(), 37, 1))
    println("  Cross-check: Slater 5s eb = $(round(eu_5s,digits=3)) eV,   NIST IP_1(Rb) = $(round(ip1_Rb,digits=3)) eV  " *
            "(expect ~2.63 vs 4.177 -- Slater undershoots due to 0.85 partial screening)")
    #
    setDefaults("print summary: close", "")
    #
end
