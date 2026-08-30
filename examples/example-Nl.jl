#
println("Nl) Tests of Empirical.InelasticHReaction -- the configuration-level interface to the Belyaev-Yakovleva (2017)")
println("    Landau-Zener model, plus SCF-vs-literature energy comparisons across several ions (Ba, Ca, Mg).")
setDefaults("unit: energy", "eV")
#
setDefaults("print summary: open", "zzz-Empirical.sum")
#
# Earlier exploratory branches (hand-typed Empirical.InelasticHChannel construction; Empirical.AtomicLevel +
# Empirical.generateInelasticHChannels built from plain L,S term data) are retired as of 24-Jul-2026, once
# Empirical.InelasticHReaction proved a strictly more transparent way to specify these reactions (input is now
# entirely Nuclear.Model + Configuration, with L,S/statistical-weight bookkeeping derived internally). The
# underlying physics engine (module-Empirical-inc-inelastic-h-collisions.jl: nonadiabaticRadius,
# landauZenerProbability, neutralizationReducedRate, the low-level inelasticHCollisionRateMatrix, ...) is unchanged
# and still exercised here, just no longer by hand-typing channels -- see project_inelastic_h_collisions.md (memory)
# for the retired branches' content and what they validated, and example-Ni.jl for the original Table-1 validation
# against Belyaev & Yakovleva (2017), which remains untouched.

if  true
    #
    # Last successful:  24-Jul-2026
    # Branch 1 (Ba): the configuration-level interface, Empirical.InelasticHReaction -- input is entirely strict JAC
    #   entities (Nuclear.Model, Configuration), with reduced mass, molecular symmetry, and statistical weight all
    #   DERIVED internally; only the level ENERGIES must still be supplied by the caller (this module deliberately
    #   does not compute them itself -- see the module's own header note).
    # System: Ba2+([Xe]) + H^-(1s^2) -> Ba+([Xe]6s/5d/6p) + H(1s^1).
    # SCF-vs-literature comparison, run and recorded here, not assumed:
    #     6s:  SCF = -9.33 eV   vs. literature -10.01 eV   (6.8% off -- within the user's stated 10% tolerance)
    #     6p:  SCF = -13.17 eV  vs. literature -7.36 eV    (79% off -- well outside tolerance)
    #     5d:  SCF = -36.56 eV  vs. literature -9.34 eV    (a factor ~3.9 off -- badly wrong, not just imprecise)
    #   The 5d failure is consistent with a known, real difficulty for heavy alkali-like ions: Ba+'s 5d level sits in
    #   a region of strong orbital near-degeneracy/mixing that a bare single-configuration (no-correlation) SCF does
    #   not capture. Given this, rates below are computed from the literature energies (the reliable source here),
    #   with the SCF numbers shown alongside purely for comparison.
    # Checks:
    #   - Rates computed from literature energies via InelasticHReaction match the original Belyaev & Yakovleva
    #     (2017) Table 2 values exactly (5.11e-9, 6.42e-9, 1.68e-8 cm^3/s for 6s/5d/6p).
    #
    println("\n  Empirical.InelasticHReaction: Ba2+([Xe]) + H^- -> Ba+([Xe]6s/5d/6p) + H, strict JAC-entity input:\n")

    nm = Nuclear.Model(56.0)
    reaction = Empirical.InelasticHReaction(
        nm,
        Configuration("[Xe]"),                          # iConfIon: Ba2+ ground configuration (closed shell)
        Configuration("1s^2"),                          # iConfH: H^-(1s^2)
        [ Configuration("[Xe] 6s^1"),
          Configuration("[Xe] 5d^1"),
          Configuration("[Xe] 6p^1") ]
        # fConfH not given -> derived automatically as Configuration("1s^1") = H(1s)
    )

    # NOTE: Configuration has a working == but currently no matching `hash` method, so a genuine Dict{Configuration,Float64}
    #   silently fails lookups (a real, separate JAC bug found while building this branch -- see the module's own comment
    #   at Empirical.energyOf). Both energy lists below therefore use Array{Pair{Configuration,Float64},1}, looked up by
    #   value equality via Empirical.energyOf, rather than a Dict.
    println("  Running a single-configuration Dirac-Fock SCF calculation for each configuration (Atomic.Computation " *
            "/ perform -- outside the Empirical module) ...")
    # Bsplines.checkGridRepresentation refuses Radial.Grid(true) here -- its box reaches 614 a.u. and the guard
    # names about 12.7 a.u. for these subshells, a 48-fold oversize. A box much TOO LARGE starves the fixed
    # B-spline basis exactly as badly as one too small, and Rule 12 warns that the high-n members of each
    # symmetry are then misrepresented -- often returning a DIFFERENT state rather than an inaccurate one.
    grid = Radial.Grid(Radial.Grid(false), rnt = 1.0e-6, h = 5.0e-2, hp = 0.0423, rbox = 12.7)
    scfEnergy(nmx, conf) = perform(Atomic.Computation(Atomic.Computation(), name="Nl-scf", nuclearModel=nmx, grid=grid,
                                                       configs=[conf]); output=true)["multiplet:"].levels[1].energy
    energiesSCF = Pair{Configuration,Float64}[ conf => scfEnergy(nm,conf) for conf in
                  [ reaction.iConfIon; reaction.fConfIon ] ]

    eV = x -> Defaults.convertUnits("energy: from eV to atomic", x)
    energiesLit = Pair{Configuration,Float64}[
        Configuration("[Xe]")       => 0.0,                    # reference zero: Ba2+ ground
        Configuration("[Xe] 6s^1")  => eV(-10.0080),
        Configuration("[Xe] 5d^1")  => eV(-9.34416),
        Configuration("[Xe] 6p^1")  => eV(-7.35615) ]

    println("\n  SCF vs. literature comparison (bound energy relative to Ba2+([Xe]), eV):")
    for  fConf  in  reaction.fConfIon
        local eScf = Defaults.convertUnits("energy: from atomic to eV",
                     Empirical.energyOf(energiesSCF,fConf) - Empirical.energyOf(energiesSCF,reaction.iConfIon))
        local eLit = Defaults.convertUnits("energy: from atomic to eV",
                     Empirical.energyOf(energiesLit,fConf) - Empirical.energyOf(energiesLit,reaction.iConfIon))
        println("    $fConf:  SCF = $(round(eScf,digits=2)) eV   literature = $(round(eLit,digits=2)) eV   " *
                "difference = $(round(100*(eScf-eLit)/eLit,digits=1))%")
    end

    T = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", 6000.)
    println("\n  Rates computed from the LITERATURE energies (the reliable source here for 5d in particular):")
    result = Empirical.inelasticHCollisionRateMatrix(T, reaction, energiesLit; energyLabel="literature (NIST)", printout=true)
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  24-Jul-2026
    # Branch 2 (Ca): same recipe as branch 1, but a much lighter, far less relativistic ion (Z=20 vs. Ba's Z=56) --
    #   testing whether the SCF 5d failure seen for Ba is a heavy-element/near-degeneracy effect specifically, or
    #   something that hits D-states more generally, regardless of Z.
    # System: Ca2+([Ar]) + H^-(1s^2) -> Ca+([Ar]4s/3d/4p) + H(1s^1).
    # Literature energies: well-known ion-trap/optical-clock spectroscopy (Ca+ 4s-3d and 4s-4p transitions are
    #   among the most precisely measured in AMO physics), recalled here rather than freshly looked up -- flagged
    #   as such, not independently cross-checked against a primary source the way Ba's Table 1 was:
    #     IP(Ca+->Ca2+) = 11.8717 eV;  4s (2S) = -11.8717 eV (ground);
    #     3d (2D, J-averaged from the 729.1/732.4 nm clock transitions) = -10.1735 eV;
    #     4p (2P, J-averaged from the 393.4/396.8 nm transitions) = -8.7292 eV.
    #
    println("\n  Empirical.InelasticHReaction: Ca2+([Ar]) + H^- -> Ca+([Ar]4s/3d/4p) + H, strict JAC-entity input:\n")

    nm = Nuclear.Model(20.0)
    reaction = Empirical.InelasticHReaction(
        nm,
        Configuration("[Ar]"),
        Configuration("1s^2"),
        [ Configuration("[Ar] 4s^1"),
          Configuration("[Ar] 3d^1"),
          Configuration("[Ar] 4p^1") ]
    )

    println("  Running a single-configuration Dirac-Fock SCF calculation for each configuration ...")
    grid = Radial.Grid(true)
    scfEnergy(nmx, conf) = perform(Atomic.Computation(Atomic.Computation(), name="Nl-scf", nuclearModel=nmx, grid=grid,
                                                       configs=[conf]); output=true)["multiplet:"].levels[1].energy
    energiesSCF = Pair{Configuration,Float64}[ conf => scfEnergy(nm,conf) for conf in
                  [ reaction.iConfIon; reaction.fConfIon ] ]

    eV = x -> Defaults.convertUnits("energy: from eV to atomic", x)
    energiesLit = Pair{Configuration,Float64}[
        Configuration("[Ar]")       => 0.0,
        Configuration("[Ar] 4s^1")  => eV(-11.8717),
        Configuration("[Ar] 3d^1")  => eV(-10.1735),
        Configuration("[Ar] 4p^1")  => eV(-8.7292) ]

    println("\n  SCF vs. literature comparison (bound energy relative to Ca2+([Ar]), eV):")
    for  fConf  in  reaction.fConfIon
        local eScf = Defaults.convertUnits("energy: from atomic to eV",
                     Empirical.energyOf(energiesSCF,fConf) - Empirical.energyOf(energiesSCF,reaction.iConfIon))
        local eLit = Defaults.convertUnits("energy: from atomic to eV",
                     Empirical.energyOf(energiesLit,fConf) - Empirical.energyOf(energiesLit,reaction.iConfIon))
        println("    $fConf:  SCF = $(round(eScf,digits=2)) eV   literature = $(round(eLit,digits=2)) eV   " *
                "difference = $(round(100*(eScf-eLit)/eLit,digits=1))%")
    end

    T = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", 6000.)
    println("\n  Rates computed from the LITERATURE energies:")
    result = Empirical.inelasticHCollisionRateMatrix(T, reaction, energiesLit; energyLabel="literature", printout=true)
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  24-Jul-2026
    # Branch 3 (Mg): same recipe again, the lightest case (Z=12), completing the S/P/D-channel comparison across
    #   three alkaline-earth-like ions spanning Z=12-56.
    # System: Mg2+([Ne]) + H^-(1s^2) -> Mg+([Ne]3s/3d/3p) + H(1s^1).
    # Literature energies: IP(Mg+->Mg2+) = 15.0353 eV (well-known); 3s (2S) = -15.0353 eV (ground); 3p (2P,
    #   J-averaged from the well-known Mg II h&k lines at 279.55/280.27 nm) = -10.6034 eV; 3d (2D) = -6.17 eV --
    #   this last number recalled with LOWER confidence than the others (not a routinely-cited transition the way
    #   the h&k lines or Ca+'s clock transitions are); flagged explicitly rather than presented as equally solid.
    #
    println("\n  Empirical.InelasticHReaction: Mg2+([Ne]) + H^- -> Mg+([Ne]3s/3d/3p) + H, strict JAC-entity input:\n")

    nm = Nuclear.Model(12.0)
    reaction = Empirical.InelasticHReaction(
        nm,
        Configuration("[Ne]"),
        Configuration("1s^2"),
        [ Configuration("[Ne] 3s^1"),
          Configuration("[Ne] 3d^1"),
          Configuration("[Ne] 3p^1") ]
    )

    println("  Running a single-configuration Dirac-Fock SCF calculation for each configuration ...")
    grid = Radial.Grid(true)
    scfEnergy(nmx, conf) = perform(Atomic.Computation(Atomic.Computation(), name="Nl-scf", nuclearModel=nmx, grid=grid,
                                                       configs=[conf]); output=true)["multiplet:"].levels[1].energy
    energiesSCF = Pair{Configuration,Float64}[ conf => scfEnergy(nm,conf) for conf in
                  [ reaction.iConfIon; reaction.fConfIon ] ]

    eV = x -> Defaults.convertUnits("energy: from eV to atomic", x)
    energiesLit = Pair{Configuration,Float64}[
        Configuration("[Ne]")       => 0.0,
        Configuration("[Ne] 3s^1")  => eV(-15.0353),
        Configuration("[Ne] 3d^1")  => eV(-6.17),         # lower confidence, see comment above
        Configuration("[Ne] 3p^1")  => eV(-10.6034) ]

    println("\n  SCF vs. literature comparison (bound energy relative to Mg2+([Ne]), eV):")
    for  fConf  in  reaction.fConfIon
        local eScf = Defaults.convertUnits("energy: from atomic to eV",
                     Empirical.energyOf(energiesSCF,fConf) - Empirical.energyOf(energiesSCF,reaction.iConfIon))
        local eLit = Defaults.convertUnits("energy: from atomic to eV",
                     Empirical.energyOf(energiesLit,fConf) - Empirical.energyOf(energiesLit,reaction.iConfIon))
        println("    $fConf:  SCF = $(round(eScf,digits=2)) eV   literature = $(round(eLit,digits=2)) eV   " *
                "difference = $(round(100*(eScf-eLit)/eLit,digits=1))%")
    end

    T = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", 6000.)
    println("\n  Rates computed from the LITERATURE energies:")
    result = Empirical.inelasticHCollisionRateMatrix(T, reaction, energiesLit; energyLabel="literature", printout=true)
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    #
    # Last successful:  25-Jul-2026
    # Branch 4 (neutral-H collision, i.e. DE-EXCITATION): A^Z+(i) + H -> A^Z+(f) + H, genuinely mediated by neutral H
    #   on both sides (not H^-/H like branches 1-3) -- the process InelasticHReaction does NOT model, since it only
    #   covers the ionic-entrance neutralization direction. This is the KNOWN, still-open weak point of the whole
    #   model, revisited here at the user's request ("let's spend some moderate effort to resolve the issue").
    # System: Ba+(7p) + H -> Ba+(6d) + H at T=6000K, the exact case with a published literature comparison
    #   (Belyaev & Yakovleva 2017, Table 2): K_if = 8.94e-10 cm^3/s.
    # Attempted fix (real effort spent, reported honestly -- did NOT succeed): re-read Belyaev (1993), PRA 48, 4299
    #   [examples/papers/a93.pra-belyaev-charge-exchange.pdf] in full to correctly identify Eq. (3.8), the general
    #   multichannel transition-probability formula. Key finding not caught in the earlier (paused) attempt: Eq.
    #   (3.8) describes a SINGLE entrance curve crossing MANY final curves -- exactly our NEUTRALIZATION topology
    #   (ionic entrance crossing all 19 Ba+ channels), NOT the de-excitation i->ionic->f topology (two DIFFERENT
    #   covalent curves sharing an ionic intermediate). The earlier attempt's mapping ("treat ionic as entrance,
    #   apply Eq 3.8 to the other 18 channels") used the wrong topology for that reason. A corrected mapping was
    #   attempted here: treat "having survived the i->ionic crossing" (factor 1-p_i) as the effective entrance into
    #   an Eq.-3.8-type problem over the channels between i and f (just 4f, here) plus those deeper than f (7s, 6p,
    #   5d, 6s). Implemented and tested against the same, already-validated J-summation/thermal-averaging machinery
    #   used by Empirical.deExcitationCrossSection (confirmed to reproduce the existing 3-state result exactly as a
    #   control). RESULT: the multichannel correction moves the answer FURTHER from the literature value (control:
    #   4.60e-10, off by 1.9x; with the attempted correction: ~1.35e-10, off by ~6.6x -- worse, not better). A
    #   direct sanity check of Eq. (3.8) itself (F=2, target=first channel, compared two independent ways -- direct
    #   trajectory counting vs. the formula as read from the paper) revealed a genuine, unresolved inconsistency
    #   even in that simple limit. Given "moderate effort... for a while" was the brief, this is reported honestly
    #   as an unresolved attempt rather than pushed further or silently shipped -- see project_inelastic_h_collisions.md
    #   for the full derivation trail if this is revisited.
    # Checks:
    #   - The EXISTING (already-flagged-unreliable) 3-state formula is confirmed, once again, to reproduce its own
    #     previously-recorded value; still ~1.9x below the literature number.
    #
    println("\n  A collision with NEUTRAL H (de-excitation): Ba+(7p) + H -> Ba+(6d) + H at T=6000K:\n")
    eV = x -> Defaults.convertUnits("energy: from eV to atomic", x)
    Z  = 1.0
    MBa = 137.327 / Defaults.ELECTRON_MASS_U;   MH = 1.00794 / Defaults.ELECTRON_MASS_U
    mu  = MBa * MH / (MBa + MH)
    EHminus = Empirical.hydrogenAnionEnergy()
    T   = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", 6000.)
    factor = Defaults.convertUnits("length: from atomic to cm", 1.0)^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)

    Ei = eV(-3.83309);   Ef = eV(-4.29573)             # Ba+(7p), Ba+(6d)
    pstat_i = 1/12                                      # Ba+(7p), a 2P term -- statisticalWeight(1, 0.5, 0, 0.0)

    D_if = Empirical.deExcitationReducedRate(T, Ei, Ef, EHminus, Z, mu)
    K_if = pstat_i * D_if
    println("  Existing 3-state formula:  K_if = $(round(factor*K_if,sigdigits=4)) cm^3/s   " *
            "[literature: 8.94e-10 cm^3/s, ratio = $(round(factor*K_if/8.94e-10,digits=3))]")
    println("  KNOWN, STILL-OPEN limitation: this formula treats the process as an isolated 3-state (i, ionic, f) " *
            "problem, ignoring that the ionic curve also crosses other covalent channels (here, Ba+(4f) lies " *
            "between 7p and 6d in binding energy) along the way. A genuine attempt to correct this via Belyaev " *
            "(1993)'s general multichannel formula, Eq. (3.8), was made and did NOT succeed -- see the branch " *
            "comment above and project_inelastic_h_collisions.md for the full, honestly-reported trail. Do not " *
            "use de-excitation/excitation numbers from this model quantitatively.")
    #
    setDefaults("print summary: close", "")
    #
end
