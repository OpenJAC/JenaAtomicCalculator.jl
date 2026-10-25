#
println("De) Apply & test the AutoIonization module with ASF from an internally generated initial- and final-state multiplet.")

setDefaults("print summary: open", "zzz-AutoIonization.sum")
setDefaults("method: continuum, Galerkin")           ## setDefaults("method: continuum, Galerkin")  "method: continuum, asymptotic Coulomb"
                                                     ## setDefaults("method: normalization, Ong-Russek") 
setDefaults("method: normalization, pure sine")      ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")

grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)

## BIORTHOGONAL TRANSFORMATION: AutoIonization.Settings gained a calcBiorthogonal field (2-Aug-2026), analog to
## PhotoIonization -- same reasoning applies essentially verbatim, since Auger decay has exactly the same
## structure as photoionization at the code level (an N-electron autoionizing initial multiplet, an
## (N-1)-electron final-ion multiplet, and a free electron generated fresh per line via
## Continuum.generateOrbitalForLevel with a dummy all-zero placeholder filling the matching slot on the initial
## side). Wired into AutoIonization.computeLines at the same point PhotoIonization does it: transforming the two
## BOUND multiplets before the free electron is generated. See PhotoIonization.Settings' calcBiorthogonal
## docstring for the full discussion; AutoIonization.Settings' own docstring cross-references it rather than
## repeating it. None of the branches below use calcBiorthogonal yet -- a natural next step, not attempted this
## session (today's focus was getting the branches themselves running at all -- see the bug notes below).
##
## FOUR real, pre-existing bugs found and fixed while making calcBiorthogonal usable and running these branches
## (none introduced today, all latent -- consistent with 5 of the original 7 branches never having a real "Last
## successful" date): (1) the keyword copy-constructor's `operator` kwarg was typed Union{Nothing,String},
## while the struct field is ::AbstractEeInteraction -- calling it with operator=CoulombInteraction() (an
## AbstractEeInteraction, not a String) always threw a TypeError, so the copy-constructor could never actually
## set a real operator; (2) three branches passed `process = Auger()` to Atomic.Computation, a keyword that
## no longer exists on its current copy-constructor (process is inferred from processSettings' type) --
## removed; (3) computeAmplitudesProperties's calcTeAuger branch called the pre-1.0-Julia `warn(...)` function,
## long removed from the language -- fixed to the modern `@warn` macro; (4) that calcTeAuger path also needs a
## real resonant Green-function gMultiplet, which none of these branches supplied (it defaults to an empty
## Multiplet()) -- so calcTeAuger was left at its default false for all of them. The LAST branch of this file
## (added 17-Aug-2026) is the one that supplies a real gMultiplet and exercises the path.

if  false
    # Last successful:  28-Aug-2026 -- REPRODUCES. Total K-shell Auger rate 9.461809e-03 a.u. against the recorded
    # 9.464e-3 (0.02 %), lifetime 2.556471e-15 s and total width 2.574689e-01 eV, i.e. 2.556 fs and 0.2575 eV
    # exactly as recorded, still matching the tabulated Ne K-shell width of ~0.24-0.27 eV.
    # Previously: 2-Aug-2026
    # K-LL Auger spectrum of F-like neon (1s-hole, 2s^2 2p^6 spectator shell): Comparison with PhD and related
    # work (an imprecise, unverifiable reference -- kept as-is, not a proper citation).
    #
    #   REPORT: total K-shell Auger rate (level 1, summed over all L-L final channels) = 9.464e-3 a.u., giving a
    #   lifetime of 2.556 fs and a TOTAL Auger width of 0.2575 eV. This matches the well-known, oft-tabulated Ne
    #   K-shell (total, i.e. essentially all-Auger since the fluorescence yield is only ~1.7%) width of
    #   ~0.24-0.27 eV (e.g. Krause & Oliver, J. Phys. Chem. Ref. Data 8, 307 (1979)) very well -- a solid,
    #   independent quantitative check despite the vague reference in the branch's own title.
    augerSettings = AutoIonization.Settings(AutoIonization.Settings(), calcAnisotropy = true, printBefore = true,
                                            lineSelection = LineSelection(true, indexPairs=[(1,0)])  )

    setDefaults("unit: rate", "a.u.")

    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(10.),
                            initialConfigs  =[Configuration("1s 2s^2 2p^6")],
                            finalConfigs    =[Configuration("1s^2 2s^2 2p^4"), Configuration("1s^2 2s 2p^5"), Configuration("1s^2 2p^6")],
                            processSettings = augerSettings )

    wb = perform(wa)
    #
elseif  false
    # Last visit:  18-Oct-2026
    # Last reproduced:  18-Oct-2026 (4 of 4 to 1e-5)
    # Last successful:  28-Aug-2026 -- REPRODUCES ALL FOUR WIDTHS EXACTLY: 7.987094e-02, 7.970066e-02,
    # 7.938559e-02 and 4.612132e-02 eV, against the recorded 0.0799, 0.0797, 0.0794, 0.0461 eV.
    # Previously: 2-Aug-2026
    # K-LL Auger spectrum of Be-like neon (1s-hole, 2s^2 2p spectator shell): Comparison with Bruch (PRA, 1991)
    #
    #   REPORT: 4 initial levels (J^P = 0-,1-,2-,1-), total widths 0.0799, 0.0797, 0.0794, 0.0461 eV
    #   respectively -- narrower than the F-like case above (branch 0), as physically expected: fewer spectator
    #   electrons (2s^2 2p vs. 2s^2 2p^6) means fewer available L-L final-state channels to sum over. All four
    #   values are internally consistent (no NaN/negative/zero rates) and in the right order of magnitude for a
    #   K-hole Auger width at this charge state. No exact Bruch (PRA 1991) numbers available this session to
    #   check absolute magnitudes; dated on physical-consistency + order-of-magnitude grounds (Rule 7).
    ## augerSettings = AutoIonization.Settings(true, true, LineSelection(true, indexPairs=[(1,0)]), 0., 1.0e6, 4, CoulombInteraction())
    augerSettings = AutoIonization.Settings(AutoIonization.Settings(), calcAnisotropy = true, printBefore = true )
    setDefaults("unit: rate", "1/s")   
    
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(10.), 
                            initialConfigs  = [Configuration("1s 2s^2 2p")],
                            finalConfigs    = [Configuration("1s^2 2s"), Configuration("1s^2 2p")], 
                            processSettings = augerSettings )

    wb = perform(wa)
    #
elseif  false
    # Last visit:  20-Oct-2026
    # Last reproduced:  20-Oct-2026 (3 of 3 numbers on the claim line)
    # Last successful:  28-Aug-2026 -- REPRODUCES, INCLUDING THE EXACT ZERO. Level 1 width 8.322900e-02 eV against
    # the recorded 0.0833; levels 2 and 3 at 1.898168e-06 and 5.153823e-06 eV, i.e. the small p_1/2-coupled widths
    # the report describes; and level 4 (J = 5/2-) gives a rate of 0.000000e+00 EXACTLY, which is the selection-rule
    # statement in the report and the most useful single check in this branch -- it is either exactly zero or wrong.
    # Previously: 2-Aug-2026
    # K-LL Auger spectrum of Li-like aluminium: Comparison with Fan et al. (PRA, 2018)
    # (relabeled 2-Aug-2026 from "Be-like" -- initialConfigs below total 3 electrons, i.e. Li-like electron
    #  count with a K-hole, not Be-like; a doubled/broken sibling branch comparing against the same Fan et
    #  al. reference for Al -- with an incompletely-edited config list left over from a Ne K-LL copy-paste,
    #  and a positional AutoIonization.Settings(...) call that no longer matched the current struct even
    #  before today's calcBiorthogonal addition -- was deleted the same day as the most obsolete/doubled of
    #  the two.)
    # NOTE (2-Aug-2026): the old positional AutoIonization.Settings(...) call here had only 7 args -- already
    # inconsistent with the (pre-existing, 10-field) struct before today's calcBiorthogonal addition, i.e. this
    # branch could never actually have run as written; rewritten below using the keyword-based copy-constructor
    # (also needed the process=Auger() and warn() fixes noted at the top of this file to actually run).
    #
    #   REPORT: initial K-hole levels span a wide range of total widths (0.0833 eV for level 1, down to
    #   ~1.9e-6 - 6.7e-7 eV for several p_1/2-coupled levels, with level 4 (J=5/2-) showing EXACTLY zero rate)
    #   -- this spread reflects genuine partial-wave/angular-momentum selection-rule suppression between the
    #   different K-hole fine-structure sublevels (s_1/2 vs. p_1/2/p_3/2/f_5/2 outgoing-electron channels have
    #   very different centrifugal-barrier overlap with the core), not noise or an error. No exact Fan et al.
    #   (PRA 2018) numbers available this session; dated on physical-consistency + selection-rule-plausibility
    #   grounds (Rule 7).
    augerSettings = AutoIonization.Settings(AutoIonization.Settings(), calcAnisotropy=true, printBefore=true,
                                            maxAugerEnergy=1.0e6, maxKappa=4, operator=CoulombInteraction())
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)
    setDefaults("unit: rate", "1/s")

    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(13.),
                            initialConfigs  = [Configuration("1s 2s^2"), Configuration("1s 2s 2p"), Configuration("1s 2p^2")],
                            finalConfigs    = [Configuration("1s^2")], 
                            processSettings = augerSettings )

    wb = perform(wa)
    #
elseif  false
    # Last visit:  20-Oct-2026
    # Last reproduced:  20-Oct-2026 (2 of 2 numbers on the claim line)
    # Last successful:  17-Oct-2026  ## widths 0.102017 and 0.101619 eV, box 33.1 a.u., 4s_1/2 at 0.797
    # Resonant L_23 - M_23 M_23 Auger spectrum of argon (2p-hole with a 4s Rydberg spectator, decaying via a
    # 3p-3p Coster-Kronig-like process while 4s watches): Comparison with Chen (PRA, 1991)
    # NOTE (2-Aug-2026): same stale-positional-call issue as the Al branch above; rewritten with keywords.
    #
    #   REPORT: the 2 selected levels (both J^P=1-, near-degenerate as expected for a spectator-electron
    #   doublet) give total widths 0.102017 and 0.101619 eV -- essentially identical to each other (consistent
    #   with the 4s spectator barely perturbing the L-shell decay) and in the right general range for an Ar
    #   2p-vacancy Auger width (L-shell Auger widths for Ar are commonly cited in the ~0.1-0.2 eV range). No
    #   exact Chen (PRA 1991) numbers available this session; dated on physical-consistency +
    #   order-of-magnitude grounds (Rule 7).
    #
    #   RE-RUN 17-Oct-2026 ON A MATCHED BOX (priority item 77), AND THE RESULT VINDICATES THE SENTENCE ABOVE.
    #   The 2-Aug numbers 0.1022 / 0.1019 eV were taken on a 10.0 a.u. box in which the 4s_1/2 reached 9.97
    #   a.u. -- extent/box = 0.996 -- and its binding energy was -0.049055 against -0.063431 on the matched
    #   box: **29 % out**.  On the matched box the widths come out 0.102017 and 0.101619 eV, i.e. they MOVED
    #   BY 0.18 % AND 0.28 %.  So the spectator is so thoroughly a spectator that a 29 % error in its own
    #   binding barely reaches the L-shell decay, which is exactly what this branch set out to show.  The box
    #   was still worth matching -- the orbital energy was simply wrong before -- but the observable was never
    #   in danger.  Of the eight branches item 77 flagged, this was the ONLY one whose orbital energies moved
    #   at all, and even here the reported numbers did not.  Zero warnings, no grid complaint.
    asfSettings   = AsfSettings(AsfSettings(), scField=Basics.HSField())
    augerSettings = AutoIonization.Settings(AutoIonization.Settings(), calcAnisotropy=true, printBefore=true,
                                            lineSelection=LineSelection(true, indexPairs=[(2,0), (4,0)]), maxAugerEnergy=1.0e6, maxKappa=4,
                                            operator=CoulombInteraction())
    # THE BOX IS MATCHED TO THE 4s SPECTATOR, NOT HAND-TYPED -- repaired 17-Oct-2026, priority item 77.
    # At rbox = 10.0 the 4s_1/2 reached 9.97 of 10.0 a.u. (extent/box = 0.996) and its energy was
    # -0.049055 against -0.063431 on a matched box: 29 % out.  That orbital IS this branch's physics --
    # the 4s Rydberg SPECTATOR whose role is to barely perturb the L-shell decay -- so a 29 % error in
    # its binding is not a detail.  `Basics.recommendedGrid` asks for 33.1 a.u.;  hp is scaled with the
    # box so the point density is unchanged.  Of the eight branches item 77 flagged this was the ONLY
    # one whose numbers actually moved: the others are singly-charged or hydrogenic, where an extent
    # ratio near 0.99 costs under 0.1 %.  What decides it is Z_eff of the flagged orbital, not the ratio.
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 2.65e-2, rbox = 33.1)
    setDefaults("unit: rate", "1/s")   
    
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(18.), 
                            initialConfigs  =[Configuration("1s^2 2s^2 2p^5 3s^2 3p^6 4s")], initialAsfSettings=asfSettings,
                            finalConfigs    =[Configuration("1s^2 2s^2 2p^6 3s^2 3p^4 4s")], finalAsfSettings=asfSettings, 
                            processSettings = augerSettings )

    wb = perform(wa)
    #
elseif  false
    # Last successful:  2-Aug-2026
    # Resonant M_45 - N_23 N_23 Auger spectrum of krypton (4d-hole with a 5p Rydberg spectator): Comparison
    # with Chen (PRA, 1991)
    # NOTE (2-Aug-2026): same stale-positional-call issue as the two branches above; rewritten with keywords.
    #
    #   REPORT: the 3 selected levels (all J^P=1-) give total widths 0.01159, 0.01770, 0.01773 eV -- smaller
    #   than the Ar L-shell case above, consistent with an M-shell (4d) vacancy typically having more
    #   competing (individually narrower) decay channels than an L-shell vacancy. Two of the three levels are
    #   near-degenerate (0.01770/0.01773 eV) while the third is distinctly narrower -- plausible fine-structure
    #   differentiation, not flagged as suspicious. No exact Chen (PRA 1991) numbers available this session;
    #   dated on physical-consistency + order-of-magnitude grounds (Rule 7).
    augerSettings = AutoIonization.Settings(AutoIonization.Settings(), calcAnisotropy=true, printBefore=true,
                                            lineSelection=LineSelection(true, indexPairs=[(4,0), (9,0), (11,0)]), maxAugerEnergy=1.0e6,
                                            maxKappa=4, operator=CoulombInteraction())
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.2e-2, rbox = 15.0)
    setDefaults("unit: rate", "1/s")   
    
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(36.),
                            initialConfigs  =[Configuration("[Ar] 3d^9 4s^2 4p^6 5p")],
                            finalConfigs    =[Configuration("[Ar] 3d^10 4s^2 4p^4 5p")],
                            processSettings = augerSettings )

    wb = perform(wa)
    #
elseif  true
    # Last visit:  23-Oct-2026   ## total 585.39e-7 a.u. over the five 3p^-3 levels (4S 180.7, 2D 266.6, 2P 138.1)
    #   against Zitnik Table I's 228e-7 -- a factor 2.57, and NO LONGER UNIFORM across the terms (3.06/2.67/2.00).
    #   This is the H - E junction of item 98; the V_ee total it replaces, measured earlier the same day, was
    #   9925.871e-7, a factor 43.5 UNIFORM across the terms (42.9/46.1/40.3).  Nothing else moved: electron
    #   energies 455.77, 452.68, 452.66, 450.59, 450.57 eV and Delta E = 9.2882/9.3773/9.4319/18.1406/18.2182 Ha
    #   are unchanged, and the first-order amplitude is still exactly zero.
    # Last successful:  unknown ...
    #
    # FURTHER WORK ON THIS BRANCH IS DELAYED, on the maintainer's decision of 03-Sep-2026, and the reason is a
    # measurement rather than a judgement: THE TEST PROPOSED BELOW DOES NOT FIT. Extending the intermediate series
    # was run on a fast machine that day, as six independent processes, at the commit dcf3a6d:
    #     Green space   levels   maxKappa   total rate
    #     n=3 only        15        2       4664.14e-7 a.u.
    #     n=3 only        15        4       4665.58e-7 a.u.     <-- the reference point, 20.5x the measurement
    #     n=3 only        15        6       4665.58e-7 a.u.
    #     + 3s holes      23        4       4268.08e-7 a.u.     <-- reproduces the -8.5% recorded below, exactly
    #     + 4p           129        4       no result: unfinished after 5 h 20 m
    #     + 4p,5p,4f     419        4       no result: abandoned after 1 h 40 m
    # TWO THINGS FOLLOW. (1) THE PARTIAL-WAVE SUM IS CLOSED, which had never been checked: maxKappa 4 and 6 agree
    # to every digit and maxKappa 2 lies within 0.03%, so the factor 20 is NOT truncation in kappa and that half of
    # the question is settled. (2) The TEA sum scales far worse than linearly in the number of intermediate levels
    # -- 15 levels cost ~210 s, 129 levels did not converge in five hours -- so work/diag-ar-tea.jl cannot be
    # completed as written. It was never "one command away"; it simply does not fit, and saying so is the point of
    # this note. What would revive it is a cheaper intermediate sum (Green-function or pseudo-states instead of an
    # explicit CI over hundreds of levels), or a smaller system in which the series can genuinely be extended.
    # Two-to-one Auger decay of the argon double L vacancy:  Ar^2+ 2p^-2  -->  Ar^3+ 3p^-3 + e-.
    # Comparison target: Zitnik et al., Phys. Rev. A 93, 021401(R) (2016), who resolved this L_23^2 - M^3 spectrum
    # and compared it with second-order perturbation theory. The paper IS in hand (apps/papers/), and its Table I
    # is quoted below; this branch exists to make the TEA path produce a first, checkable amplitude.
    #
    # WHY THIS IS THE RIGHT TEST OF calcTeAuger, and the only branch in this file that uses it. Counting orbital
    # occupations, initial -> final moves THREE bound electrons: two 3p electrons fill the two 2p holes and a third
    # 3p electron is ejected. V^(e-e) is a TWO-body operator, so by the Slater-Condon rules its matrix element
    # between states differing in three spin-orbitals vanishes identically. The channel therefore has NO first-order
    # Auger amplitude at all and exists only at second order -- which is exactly what computeTeaAmplitude sums, and
    # what the @warn in computeAmplitudesProperties checks for. The practical consequence is that the TEA amplitude
    # stands alone here, with no first-order background that could mask a wrong scale.
    #
    # THE GREEN SPACE is chosen by the same counting, not by taste. An intermediate contributes only if it differs
    # from the initial state in at most two spin-orbitals AND from the final state in at most two:
    #     2p^5 3s^2 3p^5   ... one 3p -> 2p from the initial (differs by 1); to the final one 3p -> 2p plus one
    #                          3p -> continuum (differs by 2).                                          CONTRIBUTES
    #     2p^6 3s^2 3p^4   ... two 3p -> 2p from the initial (differs by 2); to the final one 3p -> continuum
    #                          (differs by 1).                                                          CONTRIBUTES
    # The 3s-hole variants 2p^5 3s 3p^6 and 2p^6 3s 3p^5 open no new path of their own: each differs from this final
    # configuration in three or four spin-orbitals, so <f| V^(e-e) |n> vanishes for them. They are NOT harmless all
    # the same -- adding them moves the total by -8.5% (work/diag-ar-tea.jl, 17-Aug-2026). The counting argument is
    # about CSF PAIRS, whereas the intermediate LEVELS are eigenstates of whatever space is given: enlarging it
    # redistributes their mixing coefficients and re-optimises their orbitals, which reweights the paths that do
    # contribute. So the Green space is a convergence parameter here, not a set to be pruned by selection rules.
    #
    # WHY THIS BRANCH DOES NOT YET REPRODUCE THE MEASUREMENT -- an OPEN question, left open deliberately on
    # 17-Aug-2026 rather than guessed at. What is established:
    #   * the DENOMINATORS are right. JAC gets E_i - E_n = 9.16-9.42 Ha for the 2p^5 3p^5 intermediates and
    #     18.19-18.32 Ha for 2p^6 3p^4, against 270 eV and 502 eV from the known binding energies of 2p^-2,
    #     2p^-1 3p^-1 and 3p^-2. Both multiplets are 16-electron total energies, so they share a zero.
    #   * the PREFACTORS are right. ManyElectron.matrixElement_Vee builds the same object as the inner loop of the
    #     validated AutoIonization.amplitude -- same TwoParticleOperator(0, plus, true), same XL_Coulomb, same
    #     coeff.V. There is no convention mismatch between the first- and second-order routes.
    #   * the INTERMEDIATE SPACE is not converged, and that is the live suspect. 2p^5 3s^2 3p^5 is 2p^-1 3p^-2 . 3p,
    #     the n=3 member of the series whose continuum limit is the integral over eps' in Zitnik's Eq. (4); this
    #     branch has that one member and nothing above it. The five contributions printed by one partial wave sum
    #     to -7.2e-4 from terms of individual size up to 3.4e-3 -- a 14x CANCELLATION. Truncating such a sum is not
    #     a small perturbation, and since the rate goes as the square, the 20x in rate is only 4.5x in amplitude.
    # An omitted term can only make a result too LARGE if it would have cancelled what is already there, so the
    # missing continuum cannot by itself explain an excess. Extending the bound series (np, nf -- parity fixes the
    # choice, since symn must equal symi) is therefore the test to run first; work/diag-ar-tea.jl sets it up for
    # n=4,5 and 4f. IT HAS NOT BEEN RUN TO COMPLETION. Until it is, whether the continuum sum with its quadrature
    # weights and resonant denominator is genuinely needed remains UNDECIDED -- note only that the channel
    # 2p^-2 -> 2p^-1 3p^-2 + eps' is energetically OPEN (eps' ~ 235 eV), being the ordinary sequential Auger decay
    # of the first 2p hole, so a real pole does exist somewhere above the 2p^-1 3p^-2 threshold.
    #
    # Enlarging that space was IMPOSSIBLE before 17-Aug-2026: computeTeaAmplitude merged only the initial and
    # continuum subshell lists and then forced the Green multiplet into the result, so any intermediate carrying an
    # orbital the initial and final states lack -- every Rydberg member -- raised "Improper subshell order". The
    # merge now includes the intermediate list. That is why this feature had never produced a checked number: it
    # threw the moment anyone tried to converge it.
    #
    #   MEASURED against Table I of Zitnik et al., which gives the decay rates w_fq directly, in 10^-7 a.u.
    #   RE-MEASURED 23-Oct-2026, AND EVERY FIGURE IN THE PARAGRAPH BELOW THIS ONE IS SUPERSEDED.
    #     Summed over the five 3p^-3 final levels:  JAC 9925.9e-7 a.u. against the paper's 228e-7, a factor 43.5
    #     -- not 20.5. The recorded 4665.58e-7 was computed with a DEFECTIVE spin-angular coefficient and the
    #     change is a REPAIR, not a regression: bisected to 1d05ab3b (13-Sep-2026), "A two-particle coefficient
    #     was labelled with the wrong quadruple, and a rank went missing", where SpinAngular.twoParticleMoveOne
    #     interchanged the direct and crossed labels WHEN TWO CSFs DIFFER BY ONE ELECTRON, so a coefficient
    #     vanished outright. GRASP2018-validated to nine figures and carrying its own suite test. Its trigger --
    #     a ONE-electron difference -- is exactly the <n|V|i> factor of this path and exactly NOT first-order
    #     Auger, which is why branches b and c above reproduce bit-identically while this one doubled.
    #     Verified both ways on 23-Oct-2026: dcf3a6d still returns 4665.585e-7, today returns 9925.871e-7.
    #     PER FINAL TERM the excess is UNIFORM -- 4S 42.9, 2D 46.1, 2P 40.3, a spread of +-7% -- which the old
    #     6x/18x/30x note does NOT reproduce. A truncated, strongly cancelling intermediate sum cannot give a
    #     clean constant across three terms; a defect in the final-state-INDEPENDENT factor can. That is
    #     PRIORITY ITEM 98: V_ee is used where H - E belongs, so the spectator direct+exchange terms return the
    #     whole mean field <b|U_HF|a> instead of cancelling against the one-body -<b|U|a> down to the residual
    #     <b|(U_HF-U)|a>. MEASURED in an isolated clone: dropping every spectator term gives 19.0e-7, i.e. 12x
    #     TOO SMALL, so the two naive operators BRACKET the measurement 43x above and 12x below, a factor 523
    #     apart. Brillouin does not close on this branch and that residual is real -- neither endpoint is a
    #     legitimate approximation. (CORRECTED 23-Oct-2026: the reason is NOT that HSField is unconverged --
    #     it is self-consistent, scfProcedure(HSField()) = :meanFieldIteration -- but that U is not U_HF.)
    #     AND THE RESIDUAL COLLAPSES FOR NO FIELD JAC OFFERS, which was measured rather than assumed and came
    #     out opposite to the expectation: the junction ratio <n|H-E|i>/<n|V_ee|i> is 0.057-0.067 on HSField
    #     and RISES to 0.130-0.141 on Basics.ALField(). The orbital overlaps <initial|Green> are 0.995-1.000
    #     and essentially identical for the two fields, so it is not non-orthogonality between the two SCF
    #     sets either. What the ratio measures is the OFF-DIAGONAL element of whichever mean field made the
    #     orbitals, and Brillouin closes only for the Fock operator OF THE STATE -- which an AVERAGE-LEVEL
    #     field is no more than a local exchange is. So the H - E term must always be carried, not treated as
    #     a correction that a better field would remove. See tools/probeTeaJunctionOperator.jl, whose
    #     NuclearField case returns ratio = 1.0000 exactly and is what verifies the implementation.
    #     THE LEVEL-ASSIGNMENT CLAIM BELOW IS AN ARTIFACT: JAC's 452.68 eV line is 2D and the paper's 452.67 eV
    #     line is 4S, two different transitions that happen to coincide. Aligned by TERM, JAC's splittings are
    #     3.09 and 5.20 eV against the paper's 2.62 and 4.33, which match NIST Ar IV essentially exactly -- so
    #     JAC's term splittings are 18-20% TOO LARGE and the structure side is NOT right.
    #     AND THE CANCELLATION IS NOT "14-FOLD" but channel-dependent, from 0.4x (the terms add) to 390x for
    #     2P_3/2, whose net sits three orders of magnitude below its own largest term.
    #     Five of the fifteen Green levels contribute (J = 2 even): three 2p^-1 3p^-1 at Delta E = 9.29-9.43 Ha
    #     and two 3p^-2 at 18.14-18.22 Ha. The 3p^-2 band is ~2 orders down except for 2P_1/2, where it is
    #     within a factor 3 -- so it is not negligible either.
    #
    #   ITEM 98 WAS THEN BUILT AND RUN, SAME DAY, AND IT CARRIES MOST OF THE DISCREPANCY: the junction
    #   <n| H - E |i> replaces <n| V^(e-e) |i> in AutoIonization.computeTeaAmplitude, and the total falls from
    #   9925.87e-7 to 585.39e-7 a.u. -- A FACTOR 17 -- against the measured 228e-7, i.e. from 43.5x to 2.57x in
    #   rate and from 6.60x to 1.60x in amplitude. The DECAY step keeps V^(e-e) deliberately: there the two states
    #   differ in TWO orbitals, where H - E reduces to V^(e-e) identically, and the final state carries a continuum
    #   orbital whose overlap with the bound set does not vanish, so a one-body term there would need the -E <f|n>
    #   counterterm that is not implemented. See computeTeaAmplitude's own docstring.
    #
    #   THE MECHANISM IS CONFIRMED TERM BY TERM, not inferred from the total. Same partial wave, before and after:
    #       nLevel  Delta E [Ha]  character              V^(e-e)        H - E       ratio
    #          1      18.2182     3p^-2,  two-orbital   -1.637e-08   -2.420e-08     1.48
    #          4      18.1406     3p^-2,  two-orbital   -1.783e-08   -2.590e-08     1.45
    #          8       9.4319     2p^-1 3p^-1, ONE      -1.582e-07   -9.419e-09     0.060
    #         10       9.3773     2p^-1 3p^-1, ONE      +4.524e-07   +2.596e-08     0.057
    #         14       9.2882     2p^-1 3p^-1, ONE      -7.467e-08   -4.965e-09     0.067
    #   The ONE-orbital junctions lose 94% -- that is the mean field that had been counted a second time after it
    #   had already generated both orbitals. The 3p^-2 levels, which the initial state reaches by moving TWO
    #   electrons, move only 1.46x, and only through their CI mixing with the one-orbital CSFs rather than directly.
    #   So the correction lands exactly where the Brillouin argument says it must, and nowhere else.
    #
    #   MEASURED 25-Oct-2026, THREE RUNS THAT NARROW WHAT IS LEFT.  All against the experiment-anchored
    #   251 +- 46 e-7 (Zitnik's Table I column is THEORY; their own 1D_2 branching ratio puts the experiment
    #   at 1.10(20) times it).
    #     * THE ORBITAL SET IS WORTH 1.67x, AND THE BETTER FIELD IS WORSE.  Basics.ALField() in place of
    #       HSField() gives 978.23e-7 against 585.39e-7, i.e. 3.90x the anchor instead of 2.33x, and the
    #       final-state term splittings degrade with it (3.21 / 5.37 eV against HS's 3.09 / 5.18, measured
    #       2.62 / 4.33).  This is the trap the Basics.DFSField docstring records: the local-exchange field
    #       "agrees better by a CANCELLATION OF ERRORS", its more diffuse orbital lowering the Slater
    #       integrals and imitating correlation that a single-configuration calculation does not have.  So
    #       HSField's advantage here is NOT evidence that HSField is right, and tuning the orbital set toward
    #       agreement would be fitting an error -- note that Zitnik's own rule (take the set giving the
    #       LARGEST rates) does not transfer, because their rates sit BELOW experiment and ours above.
    #     * THE BOUND INTERMEDIATE SUM IS CONVERGED TO 3.4 %.  Growing the Green space from 2 configurations
    #       (15 levels) to 4 (23, adding the 3s-hole variants 2p^5 3s 3p^6 and 2p^6 3s 3p^5) moves the total
    #       585.39e-7 -> 565.27e-7, uniformly across the terms (-3.6 / -3.4 / -3.3 %).  Still exactly FIVE
    #       contributing intermediates, but not the same five -- {1,4,8,10,14} becomes {1,4,12,14,18} -- so
    #       the added configurations open NO new path and merely reweight.  THE RECORDED -8.5 % WAS TAKEN
    #       WITH V^(e-e) and is superseded; the sum is far better conditioned now that the operator is right.
    #     * THE SHAKE-DOWN OVERLAP POINTS THE WRONG WAY TO EXPLAIN AN EXCESS.  The cross-shell <2p|3p> between
    #       the SCF sets gives |<2p|3p>|^2 = 4.4e-5 to 1.7e-3 against Zitnik's sudden-approximation X_SA =
    #       2.6e-3 and the measured X = 2.2(4)e-3 -- ours is 0.65x theirs at the largest, and in that picture
    #       the rate is PROPORTIONAL to this quantity.  Within-set orthonormality is exact (1e-16), so
    #       matrixElement_Vee's assumption holds.  See tools/probeTeaJunctionOperator.jl.
    #
    #   TWO THINGS CHANGED CHARACTER, AND BOTH BEAR ON WHAT IS LEFT. The residual is no longer uniform --
    #   3.06/2.67/2.00 against the old 42.9/46.1/40.3 -- so the final-state-INDEPENDENT factor that pointed at an
    #   operator is GONE, which is what one expects of having repaired the right thing. And the sum stopped being
    #   pathological: the largest single term is now 0.67x of the net instead of 2.4x, i.e. the terms add rather
    #   than fight. A 2.57x residual on a bound-only intermediate sum that is still missing the continuum member of
    #   Zitnik's Eq. (5) is the right size and the right sign for what remains absent -- BUT THAT IS THE NEXT TEST,
    #   NOT A CONCLUSION, and 2.57x is close enough that a second, smaller error underneath it is not ruled out.
    #
    #   SUPERSEDED, kept so the earlier numbers can be recognised:
    #     level 1 is Ar^2+ 2p^-2 3P_2, CONFIRMED by its transition energies -- JAC gives 452.68 and 450.57 eV
    #     against the paper's 3P_2 rows at 452.67 and 450.05 eV, i.e. sub-eV. The structure side is right.
    #     Summed over the five 3p^-3 final levels:  JAC 4666e-7 a.u.  against  the paper's 228e-7 a.u.,  a factor
    #     20.5. Per transition the overestimate runs 6x (4S), 18x (2P), 30x (2D). Other anchors in the paper:
    #     X_exp = 2.2(4)e-3 for L23^2(1D_2), X_MBPT = 1.4e-3, and widths Gamma(L_3) = 117(5), Gamma(L23^2) = 323(15) meV.
    #     Convergence is NOT the explanation: maxKappa 4 -> 6 is bit-identical and the Green space moves it 8.5%.
    asfSettings   = AsfSettings(AsfSettings(), scField=Basics.HSField())
    #
    # ---------- STEP 1: the Green-function (intermediate) multiplet ----------
    gComp         = Atomic.Computation(Atomic.Computation(), name="STEP 1: Green-function (intermediate) multiplet",
                                       grid=grid, nuclearModel=Nuclear.Model(18.),
                                       configs = [Configuration("1s^2 2s^2 2p^5 3s^2 3p^5"),
                                                  Configuration("1s^2 2s^2 2p^6 3s^2 3p^4")], asfSettings=asfSettings )
    gMultiplet    = perform(gComp; output=true)[ResultKeys.Multiplet]
    #
    # ---------- STEP 2: the two-to-one Auger decay itself ----------
    augerSettings = AutoIonization.Settings(AutoIonization.Settings(), printBefore=true, maxAugerEnergy=1.0e6,
                                            maxKappa=4, operator=CoulombInteraction(),
                                            lineSelection=LineSelection(true, indexPairs=[(1,0)]),
                                            calcTeAuger=true, gMultiplet=gMultiplet )
    setDefaults("unit: rate", "1/s")

    wa = Atomic.Computation(Atomic.Computation(), name="Ar 2p^-2 --> 3p^-3, two-to-one Auger", grid=grid,
                            nuclearModel    = Nuclear.Model(18.),
                            initialConfigs  =[Configuration("1s^2 2s^2 2p^4 3s^2 3p^6")], initialAsfSettings=asfSettings,
                            finalConfigs    =[Configuration("1s^2 2s^2 2p^6 3s^2 3p^3")], finalAsfSettings=asfSettings,
                            processSettings = augerSettings )

    wb = perform(wa)
    #
elseif  false
    # Last visit:  25-Oct-2026   ## total first-order width 3.606288e-01 eV = 361 meV against the MEASURED
    #   323(15) meV, i.e. 1.12x;  7 final levels, electron energies 214.56-216.61 eV
    # Last successful:  unknown ...
    # THE FIRST-ORDER ANCHOR FOR THE TWO-ELECTRON BRANCH ABOVE, and the reason it exists.
    #
    # Ordinary Auger decay of the SAME Ar^2+ 2p^-2 double L vacancy that branch (f) decays two-to-one: one 2p
    # hole is filled from the M shell and a second M electron is ejected, 2p^-2 -> 2p^-1 3l^-2 + e-.  This is
    # the channel that carries essentially all of the state's width, so its total IS the L23^2 linewidth.
    #
    # WHY IT IS WORTH A BRANCH OF ITS OWN.  Branch (f) sits a factor ~2.3 above the measured two-electron rate
    # and the question is WHERE.  This branch answers half of it by elimination: it has NO intermediate states
    # at all, so it is immune to the difference between our bound Green space and Zitnik's continuum one, and
    # it exercises exactly what the second-order path INHERITS -- these orbitals, this box, this HSField, this
    # continuum normalization, this V^(e-e), this 2*pi*|A|^2.  Everything is identical to branch (f) except
    # the order.  Measured 25-Oct-2026: 361 meV against 323(15) meV, i.e. 1.12x.
    #
    #   SO THE COMMON-MODE MACHINERY IS CLEAN TO 12 %, AND BRANCH (f)'s FACTOR IS NOT IN ANYTHING THE TWO
    #   PATHS SHARE.  For comparison, Zitnik et al.'s own Eq. (3) overshoots this same width by 40 % (~452
    #   meV) with orbitals optimized for 2p^-2, so JAC is closer here than the published calculation.
    #   A second anchor in the same paper, not computed here: Gamma(L_3) = 117(5) meV.
    #
    # It also gives branch (f) a branching ratio that needs NO anchoring assumption: 565.27e-7 / 132528e-7 =
    # 4.3e-3 against the measured X = 2.2(4)e-3.  The initial states differ (3P_2 here, 1D_2 in the measured
    # X), so it is indicative rather than exact -- but it corroborates the 251 +- 46 e-7 anchor independently.
    asfSettings   = AsfSettings(AsfSettings(), scField=Basics.HSField())
    augerSettings = AutoIonization.Settings(AutoIonization.Settings(), printBefore=true, maxAugerEnergy=1.0e6,
                                            maxKappa=4, operator=CoulombInteraction(),
                                            lineSelection=LineSelection(true, indexPairs=[(1,0)]) )
    setDefaults("unit: rate", "1/s")

    wa = Atomic.Computation(Atomic.Computation(), name="Ar 2p^-2 total Auger width, first order", grid=grid,
                            nuclearModel    = Nuclear.Model(18.),
                            initialConfigs  = [Configuration("1s^2 2s^2 2p^4 3s^2 3p^6")], initialAsfSettings=asfSettings,
                            finalConfigs    = [Configuration("1s^2 2s^2 2p^5 3s^2 3p^4"),
                                               Configuration("1s^2 2s^2 2p^5 3s 3p^5"),
                                               Configuration("1s^2 2s^2 2p^5 3p^6")],   finalAsfSettings=asfSettings,
                            processSettings = augerSettings )

    wb = perform(wa)
    #
end
#
setDefaults("print summary: close", "")


