#
println("Na) Tests of empirical binding energies, ionization potentials, and total energies.")
setDefaults("unit: energy", "eV")
#
setDefaults("print summary: open", "zzz-Empirical.sum")

if  false
    #
    # Last successful:  15-Jul-2026
    # Branch 1: Empirical.bindingEnergy(Z, sh::Shell) -- compare Williams2000, Larkins1977, XrayDataBooklet
    #   for Ne (Z=10) and Ar (Z=18).
    # Tests:
    #   - Three datasets for shells 1s..3p of both atoms.
    #   - Checks that all three calls return without error and produce eV values > 0.
    # Checks:
    #   - Ne 1s: W2000 ≈ 870.2 eV; L1977 similar; XDB same source.
    #   - Ar 2p: W2000 ≈ 248.4 eV; all three should agree to < 5%.
    #   - Values decrease along 1s > 2s > 2p > 3s > 3p for both atoms.
    #
    println("\n  Empirical.bindingEnergy(Z, sh::Shell) [eV]:\n")
    println("    Shell         Williams2000   Larkins1977   XrayDataBooklet")
    for  (Z, label, shells)  in  [(10, "Ne", [Shell("1s"), Shell("2s"), Shell("2p")]),
                                  (18, "Ar", [Shell("1s"), Shell("2s"), Shell("2p"), Shell("3s"), Shell("3p")])]
        println("  $label (Z=$Z):")
        for  sh  in  shells
            e1 = Defaults.convertUnits("energy: from atomic", Empirical.bindingEnergy(Z, sh, data=PeriodicTable.Williams2000()))
            e2 = Defaults.convertUnits("energy: from atomic", Empirical.bindingEnergy(Z, sh, data=PeriodicTable.Larkins1977()))
            e3 = Defaults.convertUnits("energy: from atomic", Empirical.bindingEnergy(Z, sh, data=PeriodicTable.XrayDataBooklet()))
            sa = "    " * string(sh) * repeat(" ", max(1, 14 - length(string(sh))))
            println(sa * "$(round(e1,digits=1))           $(round(e2,digits=1))          $(round(e3,digits=1))")
        end
    end
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  15-Jul-2026
    # Branch 2: Empirical.bindingEnergy(Z, subsh::Subshell) -- subshell-resolved comparison.
    # Dataset scope:
    #   Williams2000 / Larkins1977: cover Z up to ~36; error for Z=54 (Xe).
    #   XrayDataBooklet: designed for heavy atoms and X-ray core levels (Z up to ~92).
    # Part (a): Ne (Z=10) -- all three datasets, showing j-splitting of 2p.
    # Part (b): Xe (Z=54) -- XrayDataBooklet only (W2000 and L1977 not available).
    # Checks:
    #   - Ne 2p_1/2 vs 2p_3/2: split of ~0.1 eV (relativistic j-splitting).
    #   - Xe 4d_3/2 ≈ 69.5 eV, 4d_5/2 ≈ 67.5 eV (verified in Hc Branch 1).
    #   - Xe 5s_1/2 ≈ 23.3 eV, 5p_1/2 ≈ 13.4 eV, 5p_3/2 ≈ 12.1 eV.
    #
    println("\n  Empirical.bindingEnergy(Z, subsh::Subshell) [eV]:\n")
    ##  (a) Ne: all three datasets
    println("  Ne (Z=10) -- all three datasets:")
    println("    Subshell       Williams2000   Larkins1977   XrayDataBooklet")
    for  subsh  in  [Subshell("1s_1/2"), Subshell("2s_1/2"), Subshell("2p_1/2"), Subshell("2p_3/2")]
        e1 = Defaults.convertUnits("energy: from atomic", Empirical.bindingEnergy(10, subsh, data=PeriodicTable.Williams2000()))
        e2 = Defaults.convertUnits("energy: from atomic", Empirical.bindingEnergy(10, subsh, data=PeriodicTable.Larkins1977()))
        e3 = Defaults.convertUnits("energy: from atomic", Empirical.bindingEnergy(10, subsh, data=PeriodicTable.XrayDataBooklet()))
        sa = "    " * string(subsh) * repeat(" ", max(1, 15 - length(string(subsh))))
        println(sa * "$(round(e1,digits=1))           $(round(e2,digits=1))          $(round(e3,digits=1))")
    end
    ##  (b) Xe: XrayDataBooklet only
    println("\n  Xe (Z=54) -- XrayDataBooklet only (W2000/L1977 cover Z ≤ ~36):")
    println("    Subshell       XrayDataBooklet [eV]")
    for  subsh  in  [Subshell("4d_3/2"), Subshell("4d_5/2"), Subshell("5s_1/2"), Subshell("5p_1/2"), Subshell("5p_3/2")]
        e3 = Defaults.convertUnits("energy: from atomic", Empirical.bindingEnergy(54, subsh, data=PeriodicTable.XrayDataBooklet()))
        sa = "    " * string(subsh) * repeat(" ", max(1, 15 - length(string(subsh))))
        println(sa * "$(round(e3,digits=1))")
    end
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  15-Jul-2026
    # Branch 3: Empirical.ionizationPotential and Empirical.totalEnergy for neutral atoms.
    # Tests:
    #   - ionizationPotential(Z, conf) for He, Ne, Ar ground configurations.
    #   - totalEnergy(Z, conf) with Williams2000 for the same systems.
    # Checks:
    #   - He IP ≈ 24.587 eV (NIST reference from Hc Branch 6).
    #   - Ne IP ≈ 21.565 eV,  Ar IP ≈ 15.760 eV.
    #   - Total energies (negative) decrease He > Ne > Ar in magnitude.
    #
    println("\n  Empirical.ionizationPotential(Z, conf) and totalEnergy(Z, conf) [eV]:\n")
    println("    System   NoElec   IP [eV]     totalEnergy [eV]  (Williams2000)")
    for  (Z, conf)  in  [(2,  Configuration("1s^2")),
                         (10, Configuration("1s^2 2s^2 2p^6")),
                         (18, Configuration("1s^2 2s^2 2p^6 3s^2 3p^6"))]
        ip = Defaults.convertUnits("energy: from atomic", Empirical.ionizationPotential(Z, conf))
        et = Defaults.convertUnits("energy: from atomic", Empirical.totalEnergy(Z, conf, data=PeriodicTable.Williams2000()))
        println("    Z=$Z       $(conf.NoElectrons)        $(round(ip,digits=3))       $(round(et,digits=1))")
    end
    #
    setDefaults("print summary: close", "")
    #
end
