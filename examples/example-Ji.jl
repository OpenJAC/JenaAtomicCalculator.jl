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
elseif true
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
end
