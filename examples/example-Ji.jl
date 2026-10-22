#
println("Ji) Apply & test the Plasma.KineticChargeBalanceScheme: a kinetic (rate-based, non-LTE) population balance")
println("    across several adjacent charge states, the cross-charge-state sibling to CollisionalRadiativeScheme.")

setDefaults("print summary: open", "zzz-KineticChargeBalanceScheme.sum")


# THE PHYSICS: a 3-stage carbon chain, Li-like C^3+ -- He-like C^4+ -- H-like C^5+. Within each stage, PhotoEmission
# (excited -> ground) and ImpactExcitation (both directions, by detailed balance) connect the stage's own 2-level
# reference multiplet; between adjacent stages, ImpactIonization (up) and PhotoRecombination + ThreeBodyRecombination
# (down, summed into the same matrix entry) connect the two stages' GROUND levels only. DielectronicRecombination is
# deliberately excluded. EACH of these five legs selects its OWN Plasma.AbstractRateTreatment -- an Empirical
# (tabulated/fitted) formula or JAC's own structure-based computation -- so the two branches below are not two
# different physical models, but the SAME rate equations evaluated with two different qualities of input.


if  false
    # Last visit:  22-Oct-2026
    #
    # Branch a: every leg left at its Plasma.KineticChargeBalanceScheme() default, i.e. Plasma.EmpiricalTreatment for
    #   all five legs (PE: ScaledHydrogenic; IE: VanRegemorter1962; EII: Lotz1967; RR: ScaledHydrogenic; TBR: Lotz1967).
    #   No SCF or CI is performed anywhere in this branch -- every rate comes from a Configuration-level formula --
    #   which is exactly the "quick computation for a first test" this scheme is built to offer.
    ni   = Nuclear.Model(6.)
    grid = Radial.Grid(true)

    stages = [ Plasma.ChargeStage(Configuration("1s^2 2s^1"), Configuration("1s^2 2p^1")),    # C^3+, Li-like
               Plasma.ChargeStage(Configuration("1s^2"),      Configuration("1s^1 2p^1")),    # C^4+, He-like
               Plasma.ChargeStage(Configuration("1s^1"),      Configuration("2p^1"))      ]   # C^5+, H-like

    scheme = Plasma.KineticChargeBalanceScheme(Plasma.KineticChargeBalanceScheme();
                 chargeStages=stages, electronTemperatures=[1.0e4, 1.0e5, 1.0e6])

    computation = Plasma.Computation(Plasma.Computation(), scheme=scheme, nuclearModel=ni, grid=grid,
                                     settings=Plasma.Settings(0., 1.0e10, true))
    wa = perform(computation, output=true)
    #
elseif false
    # Last visit:  22-Oct-2026
    #
    # Branch b: MIXED USE -- ImpactIonization and ThreeBodyRecombination are upgraded to Plasma.StructureBasedTreatment
    #   (real computed levels, BEBmodel), while PhotoEmission, ImpactExcitation and PhotoRecombination stay at their
    #   cheap Empirical defaults. This is the direct, same-rate-equations comparison the scheme exists to make: ONLY
    #   the two upgraded legs can move, so any population shift at Te=1e4 K is attributable to them alone.
    #
    # FINDING, quantified rather than asserted: the GROUND-slot population of the middle stage (C^4+, "1s^2") moves
    #   from 5.9773e-14 (branch a, all-empirical) to 4.4177e-16 (this branch) at Te=1e4 K -- a factor of ~135 from
    #   upgrading exactly two of the five legs. This is NOT a case where the two treatments roughly agree: at this
    #   coronal-limit temperature the ionization/recombination balance sits far out on an exponential tail, where a
    #   1967-vintage fitted cross section (Lotz) and a real computed level + a selectable BEB cross section can
    #   legitimately disagree by orders of magnitude -- exactly the comparison this scheme exists to make visible,
    #   not a defect in either route.
    ni   = Nuclear.Model(6.)
    grid = Radial.Grid(true)

    stages = [ Plasma.ChargeStage(Configuration("1s^2 2s^1"), Configuration("1s^2 2p^1")),
               Plasma.ChargeStage(Configuration("1s^2"),      Configuration("1s^1 2p^1")),
               Plasma.ChargeStage(Configuration("1s^1"),      Configuration("2p^1"))      ]

    scheme = Plasma.KineticChargeBalanceScheme(Plasma.KineticChargeBalanceScheme();
                 chargeStages=stages, electronTemperatures=[1.0e4],
                 eiiTreatment=Plasma.StructureBasedTreatment(ImpactIonization.BEBmodel()),
                 tbrTreatment=Plasma.StructureBasedTreatment(ImpactIonization.BEBmodel()))

    computation = Plasma.Computation(Plasma.Computation(), scheme=scheme, nuclearModel=ni, grid=grid,
                                     settings=Plasma.Settings(0., 1.0e10, true))
    wb = perform(computation, output=true)
    #
elseif false
    # Last visit:  22-Oct-2026
    #
    # Branch c: the MIDDLE rung of the ladder on the SAME two legs as branch b -- ImpactIonization and
    #   ThreeBodyRecombination are Plasma.SemiEmpiricalTreatment (the SAME Lotz/Saha formula as branch a, but fed a
    #   real mean-field (AL-Field) binding energy instead of a tabulated one, cf. Empirical.impactIonizationBindingEnergy),
    #   while PhotoEmission, ImpactExcitation and PhotoRecombination stay at their cheap Empirical defaults, exactly
    #   as in branch a and b. Three branches now sit on the SAME physical quantity at Te=1e4 K, one rung apart each
    #   time: a (tabulated P) -- c (mean-field P, same formula) -- b (real level + BEB cross section).
    #
    # FINDING, and it is the reason this branch earns its own entry rather than a one-line note on branch b: the
    #   ladder is NOT monotonic here. The GROUND-slot population of the middle stage (C^4+) runs
    #     branch a (tabulated)      5.9773e-14
    #     branch c (mean-field)     1.2781e-36      <-- 22 orders of magnitude BELOW both neighbours
    #     branch b (structure-based) 4.4177e-16
    #   A reader might expect "one real input" to sit somewhere between "no real input" and "two real inputs" --
    #   it does not. The mean-field (AL-Field) mean binding energy used here evidently differs sharply from BOTH the
    #   tabulated threshold branch a uses AND the true level energy branch b's BEB cross section is built on, and
    #   because the TBR rate carries this binding energy inside an exp(-(eps-P)/Te) factor at a coronal-limit
    #   temperature, a modest shift in P becomes an enormous shift in the rate. The lesson for a user of this ladder:
    #   "one rung up" is not "half the distance to ab-initio" -- a SemiEmpirical binding energy can be WORSE than a
    #   tabulated one for a specific ion/shell, and only a run against an independent source (as in branch e below)
    #   says which rung is actually closer to reality for a given case.
    ni   = Nuclear.Model(6.)
    grid = Radial.Grid(true)

    stages = [ Plasma.ChargeStage(Configuration("1s^2 2s^1"), Configuration("1s^2 2p^1")),
               Plasma.ChargeStage(Configuration("1s^2"),      Configuration("1s^1 2p^1")),
               Plasma.ChargeStage(Configuration("1s^1"),      Configuration("2p^1"))      ]

    scheme = Plasma.KineticChargeBalanceScheme(Plasma.KineticChargeBalanceScheme();
                 chargeStages=stages, electronTemperatures=[1.0e4],
                 eiiTreatment=Plasma.SemiEmpiricalTreatment(Empirical.Lotz1967()),
                 tbrTreatment=Plasma.SemiEmpiricalTreatment(Empirical.Lotz1967()))

    computation = Plasma.Computation(Plasma.Computation(), scheme=scheme, nuclearModel=ni, grid=grid,
                                     settings=Plasma.Settings(0., 1.0e10, true))
    wc = perform(computation, output=true)
    #
elseif false
    # Last visit:  22-Oct-2026
    #
    # Branch d: ALL FIVE legs at Plasma.SemiEmpiricalTreatment -- the full new rung this session's Empirical.jl/
    #   Plasma.jl hook makes possible for the first time, every leg fed a real mean-field input into its own fitted
    #   formula, with no SCF/CI beyond the cheap mean-field pass Plasma.generateChargeStageMultiplet always performs
    #   once any leg leaves Plasma.EmpiricalTreatment.
    #
    # FINDING: against branch c (same two cross-stage legs, Empirical PE/IE), the GROUND-slot population of C^4+ is
    #   essentially UNCHANGED -- 1.2781484e-36 here against 1.2781486e-36 in branch c -- because EII/TBR, not
    #   PhotoEmission/ImpactExcitation, set the ground-to-ground balance at this coronal-limit temperature. What DOES
    #   move is the population WITHIN a stage: the Li-like (C^3+) EXCITED-level slot runs 2.215e-13 (branch c,
    #   Empirical ImpactExcitation) against 5.346e-11 (here, SemiEmpirical ImpactExcitation fed a real mean-field
    #   Einstein-A) -- a real ~240x difference the ground-level comparison above is blind to, and a reminder that
    #   "the headline number" and "the whole population vector" can tell different stories about which legs matter.
    ni   = Nuclear.Model(6.)
    grid = Radial.Grid(true)

    stages = [ Plasma.ChargeStage(Configuration("1s^2 2s^1"), Configuration("1s^2 2p^1")),
               Plasma.ChargeStage(Configuration("1s^2"),      Configuration("1s^1 2p^1")),
               Plasma.ChargeStage(Configuration("1s^1"),      Configuration("2p^1"))      ]

    scheme = Plasma.KineticChargeBalanceScheme(Plasma.KineticChargeBalanceScheme();
                 chargeStages=stages, electronTemperatures=[1.0e4],
                 peTreatment =Plasma.SemiEmpiricalTreatment(Empirical.UsingJAC()),
                 ieTreatment =Plasma.SemiEmpiricalTreatment(Empirical.VanRegemorter1962()),
                 eiiTreatment=Plasma.SemiEmpiricalTreatment(Empirical.Lotz1967()),
                 rrTreatment =Plasma.SemiEmpiricalTreatment(Empirical.UsingJAC()),
                 tbrTreatment=Plasma.SemiEmpiricalTreatment(Empirical.Lotz1967()))

    computation = Plasma.Computation(Plasma.Computation(), scheme=scheme, nuclearModel=ni, grid=grid,
                                     settings=Plasma.Settings(0., 1.0e10, true))
    wd = perform(computation, output=true)
    #
elseif true
    # Last visit:  22-Oct-2026
    #
    # Branch e: the SAME all-SemiEmpirical scheme as branch d, now swept over Te = 1e4, 1e5, 3e5, 1e6 K, checked
    #   against ONE well-established, independently cited fact rather than a precise fractional-abundance table --
    #   because a precise table (e.g. CHIANTI, or Arnaud & Rothenflug 1985 itself) is computed WITH DielectronicRecombination,
    #   which this scheme deliberately excludes (cf. the scheme's own docstring), so matching such a table precisely
    #   would not be a fair test of what this scheme actually computes. The cited fact: C^3+/C IV is the standard
    #   solar TRANSITION-REGION diagnostic ion, characteristically modelled at Te ~ 1e5 K (Ne ~ 1e10 cm^-3) -- not at
    #   chromospheric Te ~ 1e4 K -- while C^4+/C V and C^5+/C VI are CORONAL ions, populated only as Te climbs well
    #   above that [J. T. Mariska, "The Solar Transition Region", Cambridge Univ. Press (1992), as cited in
    #   G. Del Zanna & H. E. Mason, Living Rev. Solar Phys. 15, 5 (2018), arXiv:1809.01618, Sect. 3].
    #
    # FINDING -- the qualitative trend REPRODUCES, the temperature scale does NOT, and both halves are worth stating:
    #     Te [K]     C^3+ (Li-like)   C^4+ (He-like)   C^5+ (H-like)
    #     1.0e4         1.0000           1.28e-36         2.86e-237
    #     1.0e5         1.0000           6.44e-06         4.83e-27
    #     3.0e5         0.9951           4.85e-03         2.91e-10
    #     1.0e6         0.8532           0.1446           2.17e-03
    #   C^3+ stays essentially undepleted through 1e4-1e5 K and only gives up ground noticeably above a few 1e5 K --
    #   the SAME ordering (lower stage survives to higher T than it would if EII had no threshold at all) and the
    #   SAME qualitative story (C IV is a transition-region ion, C V/C VI are coronal) the cited source describes.
    #   What does NOT reproduce is the ABSOLUTE scale: real coronal carbon at Te=1e6 K is close to fully stripped,
    #   dominated by C^5+/C^6+, while this run still has 85% in C^3+. Two scope limits explain the gap without
    #   invoking a defect: (1) DielectronicRecombination is excluded by design, and omitting a recombination channel
    #   should, if anything, push the balance toward HIGHER ionization at a given Te, not lower, so DR's absence is
    #   NOT what is holding C^3+ up; (2) this is a CLOSED 3-stage chain with no stage below C^3+ or above C^5+, so
    #   there is nowhere for flux to go once H-like is reached and no companion isoelectronic sequence to compare
    #   the absolute scale against -- a fair quantitative test needs the full sequence, which is outside this
    #   scheme's deliberately narrow 3-stage scope, not a bug in the 3 stages it does model.
    ni   = Nuclear.Model(6.)
    grid = Radial.Grid(true)

    stages = [ Plasma.ChargeStage(Configuration("1s^2 2s^1"), Configuration("1s^2 2p^1")),
               Plasma.ChargeStage(Configuration("1s^2"),      Configuration("1s^1 2p^1")),
               Plasma.ChargeStage(Configuration("1s^1"),      Configuration("2p^1"))      ]

    scheme = Plasma.KineticChargeBalanceScheme(Plasma.KineticChargeBalanceScheme();
                 chargeStages=stages, electronTemperatures=[1.0e4, 1.0e5, 3.0e5, 1.0e6],
                 peTreatment =Plasma.SemiEmpiricalTreatment(Empirical.UsingJAC()),
                 ieTreatment =Plasma.SemiEmpiricalTreatment(Empirical.VanRegemorter1962()),
                 eiiTreatment=Plasma.SemiEmpiricalTreatment(Empirical.Lotz1967()),
                 rrTreatment =Plasma.SemiEmpiricalTreatment(Empirical.UsingJAC()),
                 tbrTreatment=Plasma.SemiEmpiricalTreatment(Empirical.Lotz1967()))

    computation = Plasma.Computation(Plasma.Computation(), scheme=scheme, nuclearModel=ni, grid=grid,
                                     settings=Plasma.Settings(0., 1.0e10, true))
    we = perform(computation, output=true)
    #
end
