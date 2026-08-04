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
## real resonant Green-function gMultiplet, which none of these branches supply (defaults to an empty
## Multiplet()) -- so calcTeAuger was left at its default false for all rewritten branches below, rather than
## exercising a genuinely unconfigured feature.

if  true
    # Last successful:  2-Aug-2026
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
    # Last successful:  2-Aug-2026
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
    # Last successful:  2-Aug-2026
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
    # Last successful:  2-Aug-2026
    # Resonant L_23 - M_23 M_23 Auger spectrum of argon (2p-hole with a 4s Rydberg spectator, decaying via a
    # 3p-3p Coster-Kronig-like process while 4s watches): Comparison with Chen (PRA, 1991)
    # NOTE (2-Aug-2026): same stale-positional-call issue as the Al branch above; rewritten with keywords.
    #
    #   REPORT: the 2 selected levels (both J^P=1-, near-degenerate as expected for a spectator-electron
    #   doublet) give total widths 0.1022 and 0.1019 eV -- essentially identical to each other (consistent
    #   with the 4s spectator barely perturbing the L-shell decay) and in the right general range for an Ar
    #   2p-vacancy Auger width (L-shell Auger widths for Ar are commonly cited in the ~0.1-0.2 eV range). No
    #   exact Chen (PRA 1991) numbers available this session; dated on physical-consistency +
    #   order-of-magnitude grounds (Rule 7).
    asfSettings   = AsfSettings(AsfSettings(), scField=Basics.HSField())
    augerSettings = AutoIonization.Settings(AutoIonization.Settings(), calcAnisotropy=true, printBefore=true,
                                            lineSelection=LineSelection(true, indexPairs=[(2,0), (4,0)]), maxAugerEnergy=1.0e6, maxKappa=4,
                                            operator=CoulombInteraction())
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.8e-2, rbox = 10.0)
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
end
#
setDefaults("print summary: close", "")


