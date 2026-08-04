#
println("Nk) Tests of empirical excitation-autoionization (E-A), Arnaud & Rothenflug (1985).")
println("    Additive to Empirical.impactIonizationPlasmaAlpha (Lotz, direct ionization) -- does not replace it.")
setDefaults("unit: energy", "eV")
#
setDefaults("print summary: open", "zzz-Empirical.sum")

if  false
    #
    # Last successful:  24-Jul-2026
    # Branch 1: Na-like sequence E-A/D-I ratio vs. Z (cf. Fig. 5 of Arnaud & Rothenflug 1985, though not
    #   reproduced quantitatively here -- see below). Comparison at kT = I_EA (y=1). Both formula sub-ranges
    #   (Z<=16 and 18<=Z<=28) are exercised.
    # System: Na-like ions Mg+, Si3+, S5+ (low-Z range) and Ar7+, Ca9+, Fe15+ (high-Z range); channel
    #   [Ne] 3s^1 -> [Ne] (removing the 3s electron).
    # Checks:
    #   - All E-A and D-I rates positive and finite (this branch originally caught a genuine bug: a hand-
    #     transcribed piece of exponentialIntegralF1 was wrong throughout its range, not just at an isolated
    #     point, giving Inf/negative ratios here at y=1 -- fixed by using SpecialFunctions.expint exactly
    #     instead of the scanned-paper approximation; see that function's docstring).
    #   - Low-Z range (Z<=16): ratio grows monotonically with (Z-10), from ~0.08 to ~0.66 -- matches the
    #     expected qualitative trend.
    #   - High-Z range (18<=Z<=28): ratio instead DECREASES mildly (0.139 -> 0.097) at fixed y=1. Verified
    #     analytically that this is a genuine, self-consistent feature of the fitted a*I_EA Z-scaling in that
    #     sub-range (~(Z-10)^-2.98 at fixed y), not a transcription error -- my initial expectation of continued
    #     growth was an imprecise recollection of Fig. 5's shape, not a hard requirement; not corrected further.
    #
    println("\n  Na-like sequence: E-A/D-I ratio vs. Z, at kT = I_EA (cf. Fig. 5 of Arnaud & Rothenflug 1985):\n")
    ions = [ ("Mg+",  12), ("Si3+", 14), ("S5+",  16), ("Ar7+", 18), ("Ca9+", 20), ("Fe15+", 26) ]
    ratios = Float64[]
    for  (name, Z)  in  ions
        setDefaults("nuclear: charge", Float64(Z))
        iConf = Configuration("1s^2 2s^2 2p^6 3s^1");   fConf = Configuration("1s^2 2s^2 2p^6")
        IEA   = Z <= 16 ? 26.0*(Z-10) : 11.0*(Z-10)^1.5
        local eD  = Distribution.ElectronMaxwell(Defaults.convertUnits("energy: from eV to atomic", IEA))
        aEA   = Empirical.excitationAutoionizationPlasmaAlpha(eD, iConf)
        aDI   = Empirical.impactIonizationPlasmaAlpha(eD, iConf, fConf)
        push!(ratios, aEA/aDI)
        println("    $name (Z=$Z, Z-10=$(Z-10)):   I_EA = $(round(IEA, digits=1)) eV   " *
                "alpha^(EA)/alpha^(DI) = $(round(aEA/aDI, sigdigits=4))")
    end
    allPositive = all(r -> r > 0 && isfinite(r), ratios)
    monotoneLow = ratios[1] < ratios[2] < ratios[3]
    println("\n    All ratios positive & finite: $allPositive" *
            "   Low-Z range (Z<=16) monotonically increasing: $monotoneLow")
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  24-Jul-2026
    # Branch 2: Li-like sequence E-A, exercising the explicit correction factors 0.6/0.8/1.25 for C/N/O -- if
    #   these were accidentally all 1.0 (correction not applied) the three results would not differ from the
    #   uncorrected baseline in the expected proportion.
    # System: C+3, N+4, O+5 (the same 3 ions the paper itself measured and fit the correction factors to),
    #   plus Ne+7 (no correction, factor 1.0) for contrast.
    # Checks:
    #   - All 4 rates positive.
    #   - The ratio (corrected result)/(what the uncorrected formula alone would give) equals the stated factor
    #     (0.6, 0.8, 1.25) for C/N/O, and exactly 1.0 for Ne (no correction applies there).
    #   Verified: applied factors = 0.6/0.8/1.25/1.0 exactly, matching expectation; rates 1.4-1.8e-10 cm^3/s
    #   (physically sensible). This branch originally caught a real bug: the Li-like prefactor was mistranscribed
    #   as 1.60e+7 instead of 1.60e-7 (misread exponent sign), giving absurd rates of order 1e4 cm^3/s; fixed in
    #   Empirical.excitationAutoionizationPlasmaAlpha after re-checking the scanned page directly.
    #
    println("\n  Li-like sequence E-A, correction-factor check (C+3: x0.6, N+4: x0.8, O+5: x1.25, Ne+7: x1.0):\n")
    ions = [ ("C+3", 6, 0.6), ("N+4", 7, 0.8), ("O+5", 8, 1.25), ("Ne+7", 10, 1.0) ]
    for  (name, Z, expectedFactor)  in  ions
        setDefaults("nuclear: charge", Float64(Z))
        local iConf = Configuration("1s^2 2s^1")
        local IEA   = 13.6 * ( (Z-0.835)^2 - 0.25*(Z-1.62)^2 )
        local eD    = Distribution.ElectronMaxwell(Defaults.convertUnits("energy: from eV to atomic", IEA))
        local aEA   = Empirical.excitationAutoionizationPlasmaAlpha(eD, iConf)
        # Recompute the UNCORRECTED baseline directly (same formula, no ion-specific factor) to verify the
        # correction is genuinely applied, not silently a no-op.
        local y     = 1.0;   local Zeff = Z - 0.43;   local b = 1.0/(1.0 + 2.0e-4*Z^3)
        local f1v   = Empirical.exponentialIntegralF1(y)
        local Gy    = 2.22*f1v + 0.67*(1-y*f1v) + 0.49*y*f1v + 1.2*y*(1-y*f1v)
        local aBase = 1.60e-7 * b * 1.2 / (Zeff^2 * sqrt(IEA)) * exp(-y) * Gy
        local fac3  = Defaults.convertUnits("length: from atomic to cm", 1.0)^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)
        appliedFactor = (fac3*aEA) / aBase
        println("    $name (Z=$Z):   alpha^(EA) = $(round(fac3*aEA, sigdigits=6)) cm^3/s   " *
                "applied factor = $(round(appliedFactor, digits=3))   [expected $expectedFactor]")
    end
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  24-Jul-2026
    # Branch 3: total ionization rate alpha^(EII: total) = alpha^(DI) + alpha^(EA) -- the actual physical payoff.
    #   System: Fe15+ (Na-like), same channel and comparison point (kT = I_EA) as branch 1.
    # Checks:
    #   - alpha^(total) > alpha^(DI) (E-A is a genuine addition, not a correction that could go either way).
    #   - E-A/total is consistent with branch 1's independently-computed E-A/DI ratio for this same ion
    #     (0.097 -> 0.097/1.097 = 8.9% of the total), a cross-check between the two branches.
    #   Note: the paper's own statement that their Na-like total rate can exceed Lotz's by up to ~2.8x for iron
    #   is a DIFFERENT comparison (their own choice of comparison temperature/energy across the whole ionization
    #   equilibrium, not the single-point kT=I_EA rate-coefficient ratio used here) -- not reproduced by this
    #   branch, and no attempt is made to match it quantitatively.
    #   Verified: alpha^(DI)=1.978e-9, alpha^(EA)=1.922e-10 cm^3/s, E-A is 8.9% of the total (total/DI=1.1).
    #
    println("\n  Total ionization rate (D-I + E-A) for Na-like Fe15+ at kT = I_EA:\n")
    fac3 = Defaults.convertUnits("length: from atomic to cm", 1.0)^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)
    setDefaults("nuclear: charge", 26.)
    iConf = Configuration("1s^2 2s^2 2p^6 3s^1");   fConf = Configuration("1s^2 2s^2 2p^6")
    IEA   = 11.0*(26-10)^1.5
    eD    = Distribution.ElectronMaxwell(Defaults.convertUnits("energy: from eV to atomic", IEA))
    aDI   = Empirical.impactIonizationPlasmaAlpha(eD, iConf, fConf, printout=true)
    aEA   = Empirical.excitationAutoionizationPlasmaAlpha(eD, iConf, printout=true)
    aTot  = aDI + aEA
    println("    alpha^(DI)    = $(round(fac3*aDI,  sigdigits=4)) cm^3/s")
    println("    alpha^(EA)    = $(round(fac3*aEA,  sigdigits=4)) cm^3/s")
    println("    alpha^(total) = $(round(fac3*aTot, sigdigits=4)) cm^3/s   " *
            "(E-A is $(round(100*aEA/aTot, digits=1))% of the total; total/DI = $(round(aTot/aDI, sigdigits=3)))")
    #
    setDefaults("print summary: close", "")
    #
elseif  true
    #
    # Last successful:  24-Jul-2026
    # Branch 4: scope checks -- two distinct gaps in the paper's own coverage must be rejected cleanly.
    # Systems: (a) Na-like Z=17 (Cl+6, the explicit gap between the two fitted Na-like ranges),
    #          (b) Mg-like Z=11 (Na+, Z<18, where the paper states E-A is negligible and gives no formula).
    # Check: both raise informative errors, not wrong numbers.
    #   Verified: both correctly rejected with clear, specific error messages naming the gap.
    #
    println("\n  Scope checks: Na-like Z=17 gap, and Mg-like Z<18 (both unsupported):\n")
    setDefaults("nuclear: charge", 17.)
    try
        Empirical.excitationAutoionizationPlasmaAlpha(Distribution.ElectronMaxwell(1.0), Configuration("1s^2 2s^2 2p^6 3s^1"))
        println("    ERROR: should have thrown (Na-like Z=17)!")
    catch err
        println("    Correctly rejected (Na-like, Z=17 gap):  ", sprint(showerror, err))
    end
    setDefaults("nuclear: charge", 11.)
    try
        Empirical.excitationAutoionizationPlasmaAlpha(Distribution.ElectronMaxwell(1.0), Configuration("1s^2 2s^2 2p^6 3s^2"))
        println("    ERROR: should have thrown (Mg-like Z=11)!")
    catch err
        println("    Correctly rejected (Mg-like, Z=11 < 18):  ", sprint(showerror, err))
    end
    #
    setDefaults("print summary: close", "")
    #
end
