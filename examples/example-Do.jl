
println("Do) Hyperfine-induced transitions: electronic quenching, nuclear hyperfine mixing, and mixed cases.")

setDefaults("print summary: open", "zzz-HyperfineInduced.sum")
setDefaults("unit: energy", "eV")
setDefaults("unit: rate", "1/s")

## NEW FILE (06-Aug-2026). module-HyperfineInduced.jl had no example at all, and none of the scripts under
## apps/apps-wuwang-hfs-induced/ could run: they build a `Nuclear.Compound` that does not exist anywhere in src/,
## and construct Nuclear.Isomer with seven positional arguments (or scalars where arrays are required) against an
## eight-field struct. The module itself was rewritten the same day -- see its header.
##
## THE THREE KINDS of hyperfine-induced transition, all handled by one amplitude:
##   ELECTRONIC  the nucleus stays put, the electronic level changes, and the decay borrows strength from the
##               hyperfine admixture of a nearby electronic level  [branch a]
##   NUCLEAR     the electronic level stays put and the NUCLEAR level changes; two nuclear states of equal F mix.
##               This is nuclear hyperfine mixing, NHM  [branches b-g]
##   MIXED       both change -- the hyperfine electronic bridge; not covered here, and still marked "??" in the
##               manuscript itself
##
##   o)  OVERVIEW of the Be-like case      -- the cheap first pass: what mixes, and how strongly
##   a)  Be-like 33S12+  3P_0 -> 1S_0        -- classic hyperfine quenching; the regression anchor
##   b)  H-like  229Th89+                    -- the hyperfine matrix elements V_ik vs Shabaev Table I
##   c)  Li-like 229Th87+                    -- the same, at another charge state
##   d)  H-like  229Th89+ isomer lifetime    -- hours -> tens of ms, as a FUNCTION of B(M1)
##   e)  B-like  205Pb77+                    -- Wu et al., PRL: ~15 min -> ~32 ms
##   f)  205Pb q+ for q = 73, 74, 75         -- Wu et al., Atoms: also where electron spins pair
##   g)  235U -- nuclear E3 mixing, E1 photon
##
## REFERENCES, all in examples/papers/ unless noted:
##   b26.pra-hyperfine-induced.tex                    -- the authors' own draft; §II.E and §III
##   2022.pra-shabaev-hyperfine-induced-M1.pdf
##   internal/references/xp-2022.prl-shabaev-thorium-clock.pdf   -- PRL 128, 043001 (2022), Table I
##   2025.prl-wu-nuclear-mixing-boron-like.pdf        -- 205Pb77+
##   2025.atoms-wu-nuclear-mixing-pb.pdf              -- 205Pb, several charge states
##
## NUCLEAR INPUT is supplied through Nuclear.Isomer, which replaces the isotopic information of an ordinary
## structure computation. Use Nuclear.reducedTransitionAmplitude(mp, B, A, spinI) to turn a reduced transition
## probability B(multipole) into the matrix element `elementM`: it applies the Weisskopf unit and the
## sqrt((2 I_i + 1) B) definition in one place, with spinI the spin of the DECAYING nuclear state. Converting by
## hand is how apps/.../job-a-uranium.jl lost a factor sqrt(pi).
##
## NOTE that this is NOT the sqrt(8 pi/(2L+1) B) written in the manuscript. That form is not a reduced matrix
## element -- B already carries the 1/(2 I_i + 1) -- and using it made the nuclear radiation rate too fast by a
## factor of about 9300. See the docstring of Nuclear.reducedTransitionAmplitude for the full argument.
##
## A NOTE ON THE SCF for the H- and Li-like cases: a DFS potential self-interacts badly on a one- or few-electron
## system (it puts H(1s) at -0.194 instead of -0.5 a.u.), which would corrupt exactly the hyperfine matrix
## elements being compared. Basics.NuclearField() is the bare nuclear potential and is the right choice there.

grid = Radial.Grid(Radial.Grid(false), rnt = 1.0e-7, h = 3.0e-2, hp = 1.0e-2, rbox = 6.0)


if  false
    # Last successful:  06-Aug-2026
    # --- Branch o: the OVERVIEW pass, on the same Be-like 33S12+ system as branch a.
    #
    # RUN THIS FIRST, on any new system. It builds the multiplets and the full hyperfine representation, prints
    # every electronic and hyperfine level with the stable index that lineSelection expects, and ranks the
    # admixture channels that will carry the transition -- then stops, without computing one amplitude.
    #
    # THE INDICES IT PRINTS ARE HYPERFINE INDICES, and lineSelection is the ONLY way to select transitions here:
    #     lineSelection = LineSelection(true, indexPairs = [(2,1)])
    # There is deliberately no selection by electronic or by nuclear parent. Neither can express the question:
    # in 205Pb the isomer-based levels sitting on the electronic GROUND level decay by nuclear hyperfine mixing
    # in milliseconds, while isomer-based levels on EXCITED electronic levels decay electronically at 1e11 /s
    # with the nucleus still excited. Both share the same nuclear parent; only the pair identifies the line.
    #
    # WHY IT IS NEEDED. Nothing else tells a user which levels matter. The admixtures this module deals with span
    # 1.2e-2 (nuclear mixing in 229Th89+) to 1.8e-4 (electronic quenching here) to 1e-6 (the E2 channel through
    # 3P_2), and no one can guess that from a configuration list. Before this mode existed, writing the input
    # meant running the job, reading the level table, editing indices, and running again.
    #
    # ELECTRONIC AND NUCLEAR CHANNELS ARE RANKED TOGETHER, which is what makes 229Th89+ stop being a special
    # case: there the electronic list is simply empty -- one 1s level, nothing to mix with -- and a single
    # nuclear channel between the two isomers carries everything.
    nmSo  = Nuclear.Model(Nuclear.Model(16.); spinI = AngularJ64(3//2), mu = 0.6438, Q = -0.0678)
    isoSo = Nuclear.Isomer(Nuclear.Isomer(); spinI = AngularJ64(3//2), parity = Basics.plus,
                           energy = 0.0, mu = 0.6438, Q = -0.0678)
    hiSettings = HyperfineInduced.Settings(HyperfineInduced.Settings();
                        multipoles = [E1], hfMultipoles = [M1, E2], gauges = [UseCoulomb],
                        isomers = Nuclear.Isomer[isoSo],
                        calcOverview  = true,              ## <-- the whole point of this branch
                        lineSelection = LineSelection() )  ## ignored while calcOverview is true
    wo = Atomic.Computation(Atomic.Computation(), name="Do-o: overview of Be-like 33S12+", grid=grid,
                            nuclearModel        = nmSo,
                            initialConfigs      = [Configuration("1s^2 2s^2"), Configuration("1s^2 2s 2p")],
                            finalConfigs        = [Configuration("1s^2 2s^2"), Configuration("1s^2 2s 2p")],
                            processSettings     = hiSettings )
    perform(wo)
    #
elseif  false
    # Last visit:  06-Aug-2026
    # --- Branch a: hyperfine quenching of the 2s2p 3P_0 -> 2s^2 1S_0 line in Be-like 33S12+.
    #
    # PURELY ELECTRONIC, and the regression anchor of this file: no isomer, no nuclear transition, only the
    # nuclear SPIN. For I = 0 (32S) the 3P_0 level cannot decay by a single-photon E1 at all -- J = 0 to J = 0 is
    # strictly forbidden. Replacing one neutron gives 33S with I = 3/2; the hyperfine interaction then admixes the
    # nearby 3P_1 and 1P_1 levels into 3P_0 and opens an effective E1 path. The measured lifetime at a heavy-ion
    # storage ring is tau = 10.4 +- 0.5 s (Schippers/Bernhardt 2012), cited in §II.E.1 of the draft.
    #
    # NOTE HOW LITTLE HAS TO BE SAID. Only the two levels of the LINE are named -- 3P_0 decays, 1S_0 is decayed
    # to. The perturbers 3P_1 and 1P_1 are not listed at all: mixingLevels defaults to the whole multiplet, so
    # they admix automatically. Under the old interface they had to be given by hand as iAddIndices = [3,4,5],
    # which had two traps -- forgetting them gave a rate of exactly ZERO (the correct answer to the wrong
    # question), and naming them ALSO reported their own fully allowed E1 decay at 5e9 /s, swamping the quenched
    # line. The printed "role" table and the composition table below show what actually mixed in.
    nmS  = Nuclear.Model(Nuclear.Model(16.); spinI = AngularJ64(3//2), mu = 0.6438, Q = -0.0678)
    isoS = Nuclear.Isomer(Nuclear.Isomer(); spinI = AngularJ64(3//2), parity = Basics.plus,
                          energy = 0.0, mu = 0.6438, Q = -0.0678)
    hiSettings = HyperfineInduced.Settings(HyperfineInduced.Settings();
                        multipoles = [E1], hfMultipoles = [M1, E2], gauges = [UseCoulomb, UseBabushkin],
                        isomers = Nuclear.Isomer[isoS],
                        calcOverview  = false,
                        ## ONE selector, on HYPERFINE indices, read off a calcOverview run: hyperfine level 2
                        ## is the 3P_0-based F = 3/2 level, 1 is the 1S_0-based one. mixingLevels is omitted, so
                        ## the whole multiplet admixes and 3P_1 lends its strength automatically.
                        lineSelection = LineSelection(true, indexPairs = [(2,1)]),
                        printBefore = true, calcLifetimes = true )
    wa = Atomic.Computation(Atomic.Computation(), name="Do-a: hyperfine quenching in Be-like 33S12+", grid=grid,
                            nuclearModel        = nmS,
                            initialConfigs      = [Configuration("1s^2 2s^2"), Configuration("1s^2 2s 2p"),
                                                   Configuration("1s^2 2p^2")],
                            finalConfigs        = [Configuration("1s^2 2s^2"), Configuration("1s^2 2s 2p"),
                                                   Configuration("1s^2 2p^2")],
                            processSettings     = hiSettings )
    perform(wa)
    #
elseif  false
    # Last successful:  06-Aug-2026
    # --- Branch b: the hyperfine-interaction matrix elements of H-like 229Th89+, against Shabaev Table I.
    #
    # THE PREREQUISITE for everything else in this file, and a direct comparison: Table I of Shabaev et al.,
    # PRL 128, 043001 (2022) tabulates exactly the matrix elements that Hfs.computeHyperfineRepresentation builds.
    #
    #                                Th89+            Th87+ (branch c)
    #     V11/(mu(1)/muN)  [eV]     -1.109 (16)       -0.1833 (27)
    #     V22/(mu(2)/muN)  [eV]      0.783 (14)        0.1293 (25)
    #     V11              [eV]     -0.399 (10)       -0.0660 (16)
    #     V22              [eV]     -0.290 (47)       -0.0478 (77)
    #     V21/d            [eV]     -0.498 (11)       -0.0823 (20)
    #
    # The mu-normalised rows are the primary target: they are purely ELECTRONIC and carry none of the nuclear
    # moment uncertainty. Do not expect better than a few per cent -- Shabaev includes QED, nuclear structure and
    # interelectronic effects that JAC's hyperfine amplitude does not.
    #
    # THE PHYSICS: the ground state I = 5/2+ couples with the 1s_1/2 electron to F = 2, 3; the isomer I = 3/2+ to
    # F = 1, 2. The two F = 2 levels are degenerate but for the 8.356 eV nuclear excitation, and it is precisely
    # these that the hyperfine interaction mixes. Nothing else in the calculation changes.
    #
    # MEASURED against Table I on 06-Aug-2026, through work/diag-nhm-thorium.jl:
    #     V11 = -0.4206 eV  (Shabaev -0.399(10), ratio 1.054)
    #     V22 = -0.3088 eV  (Shabaev -0.290(47), ratio 1.065)
    #     V21 = +0.1263 eV  (non-zero, and V12 = V21 to all printed digits -- the matrix is symmetric)
    # This branch also prints the hyperfine constant A = 0.0088328 a.u. for the ground state.
    #
    # READ THIS TOGETHER WITH BRANCH c, because the pair says much more than either alone. Every ABSOLUTE value
    # is high by 5-7 %, but every RATIO between the two charge states reproduces Shabaev to about 0.05 %:
    #     A(Th89+)/A(Th87+)     = 6.048   vs  V11 ratio from Table I  1.109/0.1833 = 6.050
    #     V21(Th89+)/V21(Th87+) = 6.047   vs  V21/d ratio             0.498/0.0823 = 6.051
    # (V21 fell from 0.1828 to 0.1263 eV when the elementM convention was corrected on 06-Aug-2026 -- by exactly
    #  2/sqrt(8 pi/3) = 0.6911. V11 and V22 come from mu and did NOT move, so the comparisons above are intact,
    #  and the ratio test is normalisation-independent in any case.)
    # and the excess is the SAME 1.054 for V11 at BOTH charge states. Two different charge states and two
    # different nuclear states sharing one common multiplicative factor is the signature of a single missing
    # nuclear-structure effect -- Bohr-Weisskopf, the distribution of nuclear magnetisation, which JAC's point-
    # dipole hyperfine amplitude does not carry -- and NOT of an error in the electronic structure or in the
    # angular algebra, either of which would break the ratios. Do not tune anything toward the absolute numbers.
    BM1   = 0.008                       ## W.u.; theory 0.005-0.008, Tiedau 2024 gives 0.022(2)
    elemM = Nuclear.reducedTransitionAmplitude(M1, BM1, 229, AngularJ64(3//2))
    gsTh  = Nuclear.Isomer(Nuclear.Isomer(); spinI = AngularJ64(5//2), parity = Basics.plus, energy = 0.0,
                           mu =  0.360, multipoleM = [M1], elementM = [elemM])
    isTh  = Nuclear.Isomer(Nuclear.Isomer(); spinI = AngularJ64(3//2), parity = Basics.plus, energy = 8.356,
                           mu = -0.378, multipoleM = [M1], elementM = [elemM])
    asfTh = AsfSettings(AsfSettings(); scField = Basics.NuclearField())
    hiSettings = HyperfineInduced.Settings(HyperfineInduced.Settings();
                        multipoles = [M1], hfMultipoles = [M1], gauges = [UseCoulomb],
                        isomers = Nuclear.Isomer[gsTh, isTh],
                        calcOverview  = false,
                        lineSelection = LineSelection(),      ## all five lines are wanted
                        printBefore = true, calcLifetimes = true )
    wb = Atomic.Computation(Atomic.Computation(), name="Do-b: NHM in H-like 229Th89+", grid=grid,
                            nuclearModel        = Nuclear.Model(90.0, FermiNucleus()),
                            initialConfigs      = [Configuration("1s")],
                            finalConfigs        = [Configuration("1s")],
                            initialAsfSettings  = asfTh, finalAsfSettings = asfTh,
                            processSettings     = hiSettings )
    perform(wb)
    #
elseif  false
    # Last successful:  06-Aug-2026
    # --- Branch c: the same for Li-like 229Th87+, whose unpaired electron is 2s_1/2 rather than 1s_1/2.
    #
    # Table I gives V11/(mu/muN) = -0.1833(27) eV and V11 = -0.0660(16) eV, i.e. about six times smaller than in
    # the H-like ion: the 2s electron produces a correspondingly weaker field at the nucleus. Reproducing that
    # RATIO is a sharper test of the electronic side than either absolute number, because the nuclear input is
    # identical between the two branches and cancels.
    #
    # RESULT, 06-Aug-2026 (work/diag-nhm-th87.jl, and this branch for the rates):
    #     V11 = -0.06954 eV  (Shabaev -0.0660(16), ratio 1.054 -- the SAME excess as Th89+)
    #     V22 = -0.05105 eV  (Shabaev -0.0478(77), ratio 1.068)
    #     V21 = +0.02089 eV,  A = 0.0014604 a.u.   (V21 after the 06-Aug-2026 elementM correction)
    # The ratio to the H-like ion is 6.048 against Table I's 6.050. See branch b for what that pair implies.
    #
    # ONE MORE CHECK falls out here. The ground-state splitting should be 3A if the two nuclear states did not
    # mix: 3A = 0.11922 eV, and the computed 3+ -> 2+ line sits at 0.11933 eV, 0.09 % above it. In the H-like
    # ion the same comparison gives 0.7211 vs 0.7250 eV, 0.5 % -- SIX times larger, in step with the six times
    # stronger hyperfine interaction. The deviation from 3A is the level repulsion of the F = 2 pair (F = 3 has
    # no partner in the isomer and cannot shift), so it must grow with the mixing, and it does.
    #
    # Note this branch uses the default DFS potential, unlike branch b: with three electrons the self-interaction
    # that ruins a one-electron DFS calculation is no longer the dominant error.
    BM1   = 0.008
    elemM = Nuclear.reducedTransitionAmplitude(M1, BM1, 229, AngularJ64(3//2))
    gsTh  = Nuclear.Isomer(Nuclear.Isomer(); spinI = AngularJ64(5//2), parity = Basics.plus, energy = 0.0,
                           mu =  0.360, multipoleM = [M1], elementM = [elemM])
    isTh  = Nuclear.Isomer(Nuclear.Isomer(); spinI = AngularJ64(3//2), parity = Basics.plus, energy = 8.356,
                           mu = -0.378, multipoleM = [M1], elementM = [elemM])
    hiSettings = HyperfineInduced.Settings(HyperfineInduced.Settings();
                        multipoles = [M1], hfMultipoles = [M1], gauges = [UseCoulomb],
                        isomers = Nuclear.Isomer[gsTh, isTh],
                        calcOverview  = false,
                        lineSelection = LineSelection(),
                        printBefore = false, calcLifetimes = true )
    wc = Atomic.Computation(Atomic.Computation(), name="Do-c: NHM in Li-like 229Th87+", grid=grid,
                            nuclearModel        = Nuclear.Model(90.0, FermiNucleus()),
                            initialConfigs      = [Configuration("1s^2 2s")],
                            finalConfigs        = [Configuration("1s^2 2s")],
                            processSettings     = hiSettings )
    perform(wc)
    #
elseif  true
    # Last successful:  06-Aug-2026
    # --- Branch d: the isomeric lifetime of H-like 229Th89+, as a function of B(M1). THE HEADLINE RESULT.
    #
    # Nuclear hyperfine mixing shortens the 229Th isomer lifetime by up to six orders of magnitude, from a few
    # hours for the bare nucleus to a few tens of milliseconds in the H-like ion (Shabaev 2022; §III.C of the
    # draft). That is what makes the transition accessible to laser spectroscopy and is the whole point of the
    # mechanism.
    #
    # QUOTED AS A FUNCTION OF B(M1), never as a single number: B(M1) is genuinely uncertain. Theory gives
    # 0.005-0.008 W.u. (Minkov & Palffy), while a recent estimate from half-life measurements of other 229Th
    # excited states gives 0.022(2) W.u. (Tiedau 2024) -- discrepant by about a factor three. Since the induced
    # rate is proportional to B(M1), the lifetime scales inversely with it, and a single number would hide that.
    #
    # The excitation energy is now known to eight digits, 8.35573(2)(10) eV from direct laser spectroscopy
    # (Zhang et al., Nature 633, 63 (2024)); 8.356 eV is used here.
    #
    # FOUR HYPERFINE LEVELS, because there are two nuclear states: the ground state I = 5/2+ couples with the
    # 1s_1/2 electron to F = 2, 3 and the isomer I = 3/2+ to F = 1, 2. The ladder, from the bottom:
    #
    #     ground F=2   0.000 eV      isomer F=2   8.469 eV
    #     ground F=3   0.725 eV      isomer F=1   9.306 eV
    #
    # Note that F ALONE DOES NOT NAME A LEVEL here -- F = 2 occurs twice -- which is why the output tables carry
    # a "nuclear I^P" column, with * marking the isomer. Note also that the two hyperfine multiplets are
    # INVERTED relative to each other: F=3 lies above F=2 on the ground state (mu = +0.360) while F=1 lies above
    # F=2 on the isomer (mu = -0.378). Each follows the sign of its own mu, which no part of the code is told.
    #
    # RESULT, 06-Aug-2026, for B(M1) = 0.008 W.u. -- five individual transitions, and the level scheme closes:
    #
    #     F_i -> F_f   nuclear I^P    omega [eV]   rate [1/s]   what it is
    #     3+ -> 2+     5/2+ -> 5/2+      0.723        6.236      ground-state hyperfine M1; nucleus UNCHANGED
    #     2+ -> 2+     3/2+* -> 5/2+     8.465        5.576      isomer -> ground, nuclear de-excitation
    #     2+ -> 3+     3/2+* -> 5/2+     7.742        2.390      isomer -> ground
    #     1+ -> 2+     3/2+* -> 5/2+     9.304        4.444      isomer -> ground
    #     1+ -> 2+     3/2+* -> 3/2+*    0.839       14.639      the ISOMER's own hyperfine relaxation
    #
    # Five lines is also the right COUNT: M1 allows dF = 0, +-1, so isomer F=1 -> ground F=3 (dF = 2) is
    # forbidden and correctly absent, and ground F=2 is the lowest level and emits nothing. The energies are
    # internally consistent: 8.465 - 7.742 = 0.723 is the ground-state hyperfine splitting, and
    # 9.304 - 0.839 = 8.465, so the four hyperfine levels form ONE closed scheme.
    #
    # THE FULL SCAN, with the bare-nucleus lifetime for comparison (work/diag-bare-rate2.jl):
    #
    #     B(M1) [W.u.]   tau(F=1)   tau(F=2)   tau(F=3, ground)   tau_bare    enhancement
    #        0.005        57.3 ms   200.9 ms       160.8 ms        3.02 h      5.42e4
    #        0.008        52.4 ms   125.5 ms       160.4 ms        1.89 h      5.42e4
    #        0.022        37.5 ms    45.6 ms       158.2 ms        0.69 h      5.43e4
    #        0.048        24.5 ms    20.9 ms       154.4 ms        0.32 h      5.43e4
    #
    # THREE CONSISTENCY CHECKS come out of that scan, and they are worth more than any single number:
    #   * the ENHANCEMENT is 5.42e4 at every B(M1), to three digits. It must be: the induced rate goes as
    #     b^2 ~ B(M1) through the mixing coefficient, the bare rate goes as B(M1) directly, so the ratio cannot
    #     depend on B(M1) at all. Nothing in the code enforces this -- the two rates are computed by completely
    #     different routes -- so its constancy is a genuine check of the normalisation.
    #   * the F = 3 (ground) row stays at 155-161 ms across the whole scan, i.e. essentially INDEPENDENT of
    #     B(M1), which is what a purely electronic hyperfine transition must be: no nuclear matrix element
    #     enters it. The slight drift is the feedback of the mixing on the level composition.
    #   * the isomeric F = 2 rate is LINEAR in B(M1) to three digits (4.978 -> 47.87 for a factor 9.6).
    #
    # AGAINST SHABAEV. The measured-B(M1) case, 0.022 W.u. (Tiedau 2024), gives tau(F=2) = 45.6 ms -- the "tens
    # of ms" the PRL reports. At the theoretical 0.005-0.008 W.u. it is 125-200 ms. State it as a function of
    # B(M1), never as one number: the enhancement is what the calculation predicts robustly, the absolute
    # lifetime inherits the factor-3 spread in B(M1) directly. A further ~11 % is expected from the 5.4 %
    # Bohr-Weisskopf excess in V21, which enters the rate quadratically; see branch b.
    #
    # THE NORMALISATION WAS WRONG UNTIL 06-Aug-2026 and was found by exactly this bare-nucleus comparison, which
    # is why it is kept in the file. The nuclear radiation term was too fast by ~9300 because
    # Nuclear.reducedTransitionAmplitude followed the manuscript's sqrt(8 pi/(2L+1) B) -- not a reduced matrix
    # element, since B already carries 1/(2 I_i + 1). Correcting it to sqrt((2 I_i + 1) B) collapsed the field
    # factor to sqrt((L+1)/L) (alpha omega)^L/(2L+1)!!, with a further alpha for magnetic multipoles from the
    # 1/c of the magnetic radiation operator. V11 and V22 are untouched by this (they come from mu, not from
    # elementM), so branches b and c did not move; V21 fell by 2/sqrt(8 pi/3) = 0.6911 exactly.
    #
    # READ THE F = 1 ROW WITH CARE, and this is why it is FASTER than F = 2 without being a shorter nuclear
    # lifetime. Its total, 19.08/s, is the sum of a nuclear de-excitation (4.444/s) and the internal hyperfine
    # relaxation F = 1 -> F = 2 WITHIN the isomer (14.64/s), the latter dominating. That relaxation is a real
    # decay of the level but leaves the nucleus excited, so the isomer is NOT destroyed in 52.4 ms: it drops to
    # the isomeric F = 2 in about 68 ms and de-excites from there. Only the F = 2 row is a nuclear lifetime
    # outright; the nuclear channel out of F = 1 alone would be 225 ms.
    BM1s  = [0.005, 0.008, 0.022, 0.048]
    for  BM1 in BM1s
        println("\n", "="^110)
        println("  229Th89+ with B(M1) = $BM1 W.u.")
        println("="^110)
        elemM = Nuclear.reducedTransitionAmplitude(M1, BM1, 229, AngularJ64(3//2))
        gsTh  = Nuclear.Isomer(Nuclear.Isomer(); spinI = AngularJ64(5//2), parity = Basics.plus, energy = 0.0,
                               mu =  0.360, multipoleM = [M1], elementM = [elemM])
        isTh  = Nuclear.Isomer(Nuclear.Isomer(); spinI = AngularJ64(3//2), parity = Basics.plus, energy = 8.356,
                               mu = -0.378, multipoleM = [M1], elementM = [elemM])
        asfTh = AsfSettings(AsfSettings(); scField = Basics.NuclearField())
        hiSettings = HyperfineInduced.Settings(HyperfineInduced.Settings();
                            multipoles = [M1], hfMultipoles = [M1], gauges = [UseCoulomb],
                            isomers = Nuclear.Isomer[gsTh, isTh],
                            calcOverview  = false,
                            lineSelection = LineSelection(),   ## all five: the checks need the ground-state line
                            printBefore = false, calcLifetimes = true )
        wd = Atomic.Computation(Atomic.Computation(), name="Do-d: 229Th89+ isomer lifetime", grid=grid,
                                nuclearModel        = Nuclear.Model(90.0, FermiNucleus()),
                                initialConfigs      = [Configuration("1s")],
                                finalConfigs        = [Configuration("1s")],
                                initialAsfSettings  = asfTh, finalAsfSettings = asfTh,
                                processSettings     = hiSettings )
        perform(wd)
    end
    #
elseif  false
    # Last visit:  06-Aug-2026
    # --- Branch e: B-like 205Pb77+, from Wang/Li/Wang, examples/papers/2025.prl-wu-nuclear-mixing-boron-like.pdf.
    #
    # An INDEPENDENT check of the same mechanism on a different nucleus, a different electronic shell and a
    # different order of magnitude. The 205Pb isomer sits at 2.329 keV -- some 280 times higher than the 229Th
    # one -- and the transition is (I_e = 1/2-, 2p J=1/2) F_i -> (I_g = 5/2-, 2p J=1/2) F_f. The half-life drops
    # from roughly 15 min to ~32 ms, about four orders of magnitude rather than six.
    #
    # The valence electron is 2p_1/2, not an s electron: a useful test that nothing in the machinery quietly
    # assumes an s-wave unpaired electron.
    ## THE NUCLEAR TRANSITION IS E2, NOT M1, and this is forced by angular momentum, not chosen: I_e = 1/2- to
    ## I_g = 5/2- has |dI| = 2 with no parity change, so <I_g||M1||I_e> vanishes identically by the triangle rule
    ## and E2 is the lowest allowed multipole. The same applies to the hyperfine MIXING: only its E2 part can
    ## connect these two nuclear states at all, so hfMultipoles MUST contain E2 or the whole effect is zero.
    ## B(E2) is a PLACEHOLDER here -- take the real value from the paper before quoting any absolute rate.
    elemPb = Nuclear.reducedTransitionAmplitude(E2, 1.0, 205, AngularJ64(1//2))    ## B(E2) = 1 W.u. PLACEHOLDER
    gsPb   = Nuclear.Isomer(Nuclear.Isomer(); spinI = AngularJ64(5//2), parity = Basics.minus, energy = 0.0,
                            mu = 0.7117, multipoleM = [E2], elementM = [elemPb])
    isPb   = Nuclear.Isomer(Nuclear.Isomer(); spinI = AngularJ64(1//2), parity = Basics.minus, energy = 2329.0,
                            mu = -0.7345, multipoleM = [E2], elementM = [elemPb])
    hiSettings = HyperfineInduced.Settings(HyperfineInduced.Settings();
                        multipoles = [M1, E2], hfMultipoles = [M1, E2], gauges = [UseCoulomb],
                        isomers = Nuclear.Isomer[gsPb, isPb],
                        calcOverview  = false,
                        ## Hyperfine levels 3 and 4 are the ISOMER-based ones, 1 and 2 the ground-based ones
                        ## (from calcOverview). Without this selection the ordinary 2p_3/2 -> 2p_1/2 M1 line at
                        ## 1e11 /s is reported alongside and buries the effect by ten orders of magnitude.
                        lineSelection = LineSelection(true, indexPairs = [(3,1), (3,2), (4,1), (4,3)]),
                        printBefore = false, calcLifetimes = true )
    we = Atomic.Computation(Atomic.Computation(), name="Do-e: NHM in B-like 205Pb77+", grid=grid,
                            nuclearModel        = Nuclear.Model(82.0, FermiNucleus()),
                            initialConfigs      = [Configuration("1s^2 2s^2 2p")],
                            finalConfigs        = [Configuration("1s^2 2s^2 2p")],
                            processSettings     = hiSettings )
    perform(we)
    #
elseif  false
    # Last visit:  06-Aug-2026
    # --- Branch f: 205Pb q+ for q = 73, 74, 75, from Wang/Wang, examples/papers/2025.atoms-wu-nuclear-mixing-pb.pdf.
    #
    # The draft notes that these charge states show hyperfine-induced nuclear decay "even if the spin of the
    # electrons can pair", i.e. beyond the single-unpaired-electron picture that the 229Th and B-like 205Pb cases
    # share. This branch therefore tests something the others cannot: that the machinery does not implicitly
    # assume one unpaired electron, and that a closed or paired shell still supports the mixing.
    for  (q, confs) in [ (73, [Configuration("1s^2 2s^2 2p^5")]),
                         (74, [Configuration("1s^2 2s^2 2p^4")]),
                         (75, [Configuration("1s^2 2s^2 2p^3")]) ]
        println("\n", "="^110);   println("  205Pb$(q)+");   println("="^110)
        elemPb = Nuclear.reducedTransitionAmplitude(E2, 1.0, 205, AngularJ64(1//2))   ## PLACEHOLDER, see branch e
        gsPb   = Nuclear.Isomer(Nuclear.Isomer(); spinI = AngularJ64(5//2), parity = Basics.minus, energy = 0.0,
                                mu = 0.7117, multipoleM = [E2], elementM = [elemPb])
        isPb   = Nuclear.Isomer(Nuclear.Isomer(); spinI = AngularJ64(1//2), parity = Basics.minus, energy = 2329.0,
                                mu = -0.7345, multipoleM = [E2], elementM = [elemPb])
        hiSettings = HyperfineInduced.Settings(HyperfineInduced.Settings();
                            multipoles = [M1, E2], hfMultipoles = [M1, E2], gauges = [UseCoulomb],
                            isomers = Nuclear.Isomer[gsPb, isPb],
                            calcOverview  = false,
                            lineSelection = LineSelection(),   ## indices differ per charge state; see the header
                            printBefore = false, calcLifetimes = true )
        wf = Atomic.Computation(Atomic.Computation(), name="Do-f: NHM in 205Pb$(q)+", grid=grid,
                                nuclearModel        = Nuclear.Model(82.0, FermiNucleus()),
                                initialConfigs = confs, finalConfigs = confs, processSettings = hiSettings )
        perform(wf)
    end
    #
elseif  false
    # Last visit:  06-Aug-2026
    # --- Branch g: 235U -- an E3 hyperfine interaction with an E1 photon.
    #
    # THE ONE CASE with a hyperfine multipole other than M1, and therefore the branch that exercises that path
    # end to end. Read the two multipoles carefully, because the draft's own subsection title invites confusion:
    # the BARE nuclear decay I_e = 1/2+ -> I_g = 7/2- is E3 (dI = 3 with a parity change), and it is the E3
    # nuclear moment that mixes the states -- but the photon that is finally emitted comes from the ELECTRONIC E1
    # operator. Hence multipoles = [E1] and hfMultipoles = [E3].
    #
    # dE = 76.737(18) eV, B(E3) = 0.036 W.u. Nuclear.reducedTransitionAmplitude gives
    # <I_g||W(E3)||I_e> = 1.390e-13 a.u., which is the value the draft still leaves as a TODO. Note it is NOT the
    # 7.84e-14 of apps/.../job-a-uranium.jl: that script computes sqrt(8/7 B) instead of sqrt(8 pi/7 B) and is
    # therefore low by sqrt(pi) = 1.772.
    elemU = Nuclear.reducedTransitionAmplitude(E3, 0.036, 235, AngularJ64(1//2))
    gsU   = Nuclear.Isomer(Nuclear.Isomer(); spinI = AngularJ64(7//2), parity = Basics.minus, energy = 0.0,
                           mu = -0.38, multipoleM = [E3], elementM = [elemU])
    isU   = Nuclear.Isomer(Nuclear.Isomer(); spinI = AngularJ64(1//2), parity = Basics.plus, energy = 76.737,
                           mu =  0.0,  multipoleM = [E3], elementM = [elemU])
    hiSettings = HyperfineInduced.Settings(HyperfineInduced.Settings();
                        multipoles = [E1], hfMultipoles = [E3], gauges = [UseCoulomb, UseBabushkin],
                        isomers = Nuclear.Isomer[gsU, isU],
                        ## OVERVIEW FIRST for this branch, and not as a formality: 3d^2 and 3p^5 3d^3 are open
                        ## shells, so the electronic multiplets are large and the full hyperfine basis larger
                        ## still. Computing every allowed line would be an enormous job and almost all of it
                        ## uninteresting -- the E3 mixing is carried by a handful of channels. Run the overview,
                        ## read the ranked channels, then name those few pairs in lineSelection.
                        calcOverview  = true,
                        lineSelection = LineSelection(),
                        printBefore = false, calcLifetimes = true )
    wg = Atomic.Computation(Atomic.Computation(), name="Do-g: hyperfine-induced E1 decay of the 235U isomer", grid=grid,
                            nuclearModel        = Nuclear.Model(92.0, FermiNucleus()),
                            ## THE CONFIGURATION SPACE CANNOT BE MADE SMALL HERE, and trying taught the lesson.
                            ## An E3 hyperfine operator is rank 3, so a non-zero electronic matrix element needs
                            ## |J_r - J_s| <= 3 <= J_r + J_s, i.e. J_r + J_s >= 3, AND -- being parity-odd -- it
                            ## needs the two electronic levels to have OPPOSITE parity. Li-like U89+ was tried as
                            ## a cheap substitute and returned "every hyperfine level is pure": its levels are
                            ## J = 1/2, 1/2, 3/2, whose largest sum is 2, so every E3 element vanishes by the
                            ## triangle rule. The effect genuinely requires a rich open-shell structure, which is
                            ## why this branch is expensive; that is physics, not an inefficiency.
                            initialConfigs      = [Configuration("[Ne] 3s^2 3p^6 3d^2"), Configuration("[Ne] 3s^2 3p^5 3d^3")],
                            finalConfigs        = [Configuration("[Ne] 3s^2 3p^6 3d^2"), Configuration("[Ne] 3s^2 3p^5 3d^3")],
                            processSettings     = hiSettings )
    perform(wg)
    #
end
#
setDefaults("print summary: close", "")
