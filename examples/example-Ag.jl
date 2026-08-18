#
println("Ag) Apply & test the jj-LS transformation of levels from a given multiplet.")

if  false
    # Last successful:  18-Aug-2026
    # Branch a: ONE open shell -- baseline, uses the existing (trusted) LSjj.OneOpenShell code path;
    #   no recoupling coefficient needed at all.
    wa = Atomic.Computation(Atomic.Computation(), name="jj-LS: one open shell", grid=Radial.Grid(true), nuclearModel=Nuclear.Model(10.),
                            configs=[Configuration("[Ne] 3p")],
                            asfSettings=AsfSettings(AsfSettings(); jjLS=LSjjSettings(true, 0.05, 0.001, LevelSelection())) )
    wb = perform(wa)
    #
elseif  false
    # Last successful:  18-Aug-2026
    # Branch b: TWO open shells -- exercises LSjj.GeneralOpenShells(2) (Gaigalas, Atoms 14, 20 (2026),
    #   Eq. 17), routed to by Basics.extractFromConfiguration(Basics.OpenShellNumber(),...). Verified
    #   to reproduce the formerly-trusted, now-retired LSjj.TwoOpenShells code exactly (max|old-new| =
    #   2.2e-16 across all CSFs of this configuration) before that code was removed, so this is also the
    #   regression reference the 3-, 4- and 5-open-shell branches below build confidence from.
    wa = Atomic.Computation(Atomic.Computation(), name="jj-LS: two open shells", grid=Radial.Grid(true), nuclearModel=Nuclear.Model(10.),
                            configs=[Configuration("[Ne] 3s 3p")],
                            asfSettings=AsfSettings(AsfSettings(); jjLS=LSjjSettings(true, 0.05, 0.001, LevelSelection())) )
    wb = perform(wa)
    #
elseif  false
    # Last successful:  18-Aug-2026
    # Branch c: THREE open shells -- exercises LSjj.GeneralOpenShells(3). The formerly-existing
    #   LSjj.ThreeOpenShells code (never independently validated before -- present in the file for years
    #   but never exercised past the 2-open-shell case; now removed) turned out to have a real bug: for CSFs whose
    #   accumulated jj-coupling after the middle shell can take more than one value (e.g. two CSFs with
    #   identical 3s/3d occupation but 3p contributing an accumulated total of either 1 or 2), it silently
    #   summed over both instead of using the one value this particular CSF actually has -- inflating some
    #   physical level weights up to 2.2 instead of the required exactly-1.0. GeneralOpenShells(3) does not
    #   have this bug (fixed by treating the accumulated total as a value read off the CSF itself, not a
    #   free sum, for every shell except the last) and gives exactly weight=1.0 for all 23 physical levels
    #   of [Ne] 3s 3p 3d.
    wa = Atomic.Computation(Atomic.Computation(), name="jj-LS: three open shells", grid=Radial.Grid(true), nuclearModel=Nuclear.Model(26.),
                            configs=[Configuration("[Ne] 3s 3p 3d")],
                            asfSettings=AsfSettings(AsfSettings(); jjLS=LSjjSettings(true, 0.05, 0.001, LevelSelection())) )
    wb = perform(wa)
    #
elseif  false
    # Last successful:  18-Aug-2026
    # Branch d: FOUR open shells -- exercises LSjj.GeneralOpenShells(4); no prior hand-derived code
    #   exists for this shell count, so this is the first-ever 4-open-shell jj-LS transformation in JAC.
    #   Adding the outer 4s shell after the 3s/3p/3d holes surfaced a second real bug: an l=0 shell (only
    #   a single j=1/2 subshell, no minus/plus split at all) folded into the recursion at a non-first
    #   position gave EXACTLY ZERO for every CSF, because the recoupling formula's "diagram B" (6-j) step
    #   -- which exists specifically to reconcile a shell's minus/plus jj-split -- was fed a meaningless
    #   sentinel value for the nonexistent minus subshell. Fixed by skipping diagram B for l=0 shells
    #   (folding in via diagram A, the 9-j, alone, keeping its own normalization prefactor) -- confirmed by
    #   exactly weight=1.0 for all 46 physical levels of [Ne] 3s 3p 3d 4s.
    wa = Atomic.Computation(Atomic.Computation(), name="jj-LS: four open shells", grid=Radial.Grid(true), nuclearModel=Nuclear.Model(26.),
                            configs=[Configuration("[Ne] 3s 3p 3d 4s")],
                            asfSettings=AsfSettings(AsfSettings(); jjLS=LSjjSettings(true, 0.05, 0.001, LevelSelection())) )
    wb = perform(wa)
    #
elseif  false
    # Last successful:  18-Aug-2026
    # Branch e: FIVE open shells -- exercises LSjj.GeneralOpenShells(5), representative of the kind of
    #   configuration that can arise from inner-shell holes in a cascade (multiple simultaneously open
    #   subshells across several n). Confirms the fixes from branches c/d generalize: exactly weight=1.0
    #   for all 249 physical levels of [Ne] 3s 3p 3d 4s 4p.
    wa = Atomic.Computation(Atomic.Computation(), name="jj-LS: five open shells", grid=Radial.Grid(true), nuclearModel=Nuclear.Model(26.),
                            configs=[Configuration("[Ne] 3s 3p 3d 4s 4p")],
                            asfSettings=AsfSettings(AsfSettings(); jjLS=LSjjSettings(true, 0.05, 0.001, LevelSelection())) )
    wb = perform(wa)
    #
elseif  false
    # Last successful:  18-Aug-2026
    # Branch f: GROUND-LEVEL CLASSIFICATION of an open d-shell, Co^2+ (Co III) 3d^7, against NIST ASD.
    #   A d^7 shell is a control case: it is covered by the LS_jj_d_7 table that JAC has always carried, so
    #   this branch tests the CLASSIFICATION itself, not the f-shell machinery added on 10-Aug-2026. It is
    #   the natural companion to branch g below, which does exercise the new work.
    #
    #   NIST ASD (Co III, levels in eV above the ground level) against JAC (DFS, no correlation, no Breit):
    #        term        NIST        JAC       deviation
    #        a ^4F_9/2   0.000000    0.0000    ground level, and JAC assigns it correctly
    #        a ^4F_7/2   0.104315    0.1070      +2.6 %
    #        a ^4F_5/2   0.179942    0.1858      +3.3 %
    #        a ^4F_3/2   0.231540    0.2398      +3.6 %
    #        a ^4P_5/2   1.884877    2.3249     +23   %
    #   JAC gives the ground level as ^4F_9/2 with a weight of 0.99765, the rest being a 0.2 % ^2G_9/2
    #   admixture -- the assignment is unambiguous, which is what makes this configuration a good test.
    #   The pattern of the deviations is the physics: the FINE STRUCTURE inside the ^4F term, which is a
    #   one-electron spin-orbit effect, comes out to 3-4 %, whereas the ^4F-^4P TERM separation, which is an
    #   electron-electron correlation quantity, is 23 % too large. A single-configuration DFS calculation is
    #   expected to behave in exactly this way, and the split verdict is more informative than either number
    #   alone would be.
    grid = Radial.Grid(Radial.Grid(false), rnt=2.0e-6, h=5.0e-2, hp=1.0e-2, rbox=12.0)
    wa   = Atomic.Computation(Atomic.Computation(), name="jj-LS: Co^2+ 3d^7 ground level", grid=grid,
                              nuclearModel=Nuclear.Model(27.0), configs=[Configuration("[Ar] 3d^7")],
                              asfSettings=AsfSettings(AsfSettings(); jjLS=LSjjSettings(true, 0.05, 0.001, LevelSelection())) )
    wb = perform(wa)
    #
elseif  true
    # Last successful:  18-Aug-2026
    # Branch g: GROUND-LEVEL CLASSIFICATION of an open f-shell, Dy^3+ (Dy IV) 4f^9, against NIST ASD.
    #   THIS BRANCH EXERCISES THE WORK OF 10-Aug-2026 and could not have been run before it. Two things had
    #   to be in place. First, the f-shell LS-jj coefficient tables were present in the repository but their
    #   include was commented out, so every open f-shell raised UndefVarError. Second, 4f^9 is MORE than
    #   half filled, and no tables exist beyond f^7 -- in JAC or in GRASP. It is reached through the
    #   electron-hole symmetry of Dyall & Grant, J. Phys. B 15 (1982) L371 (Eq. (A.4) of Gaigalas &
    #   Fritzsche, Comput. Phys. Commun. 149 (2002) 39), which maps f^9 onto the tabulated f^5 by reversing
    #   the j_- occupation and attaching the phase (-1)^(Qm+Qp-Q). So this is a genuine test of that rule
    #   against measured level positions, not merely of the tables.
    #
    #   NIST ASD (Dy IV, levels in eV above the ground level) against JAC (DFS, no correlation, no Breit):
    #        term         NIST      JAC      deviation
    #        ^6H*_15/2    0.000     0.000    ground level, and JAC assigns it correctly, odd parity
    #        ^6H*_13/2    0.429     0.442      +3.0 %
    #        ^6H*_11/2    0.717     0.747      +4.2 %
    #        ^6H*_7/2     1.123     1.164      +3.7 %
    #   JAC gives the ground level as ^6H_15/2 with a weight of 0.93968, the remainder being ^4I_15/2
    #   admixtures of 4.5 % and 1.2 %. That the leading term carries only 94 % rather than the 99.8 % of the
    #   3d^7 case above is NOT a defect: 4f^9 in a heavy lanthanide is a strongly spin-orbit-coupled,
    #   genuinely intermediate-coupling system, and the LSJ label is a leading component rather than a good
    #   quantum number. Reporting that weight alongside the label is the point of the expansion.
    #   The composition now prints as an aligned table of PERCENTAGES (competition item C7, 18-Aug-2026); the
    #   weights above are the same numbers, checked component by component -- 1583 of them across branches a-g,
    #   zero moved.  Each branch passes printWeight = 0.001, i.e. a 0.1 % cut, because the default from
    #   LSjjSettings(makeIt) is 0.1 = 10 % and would drop the very admixtures these branches document.
    #
    #   The two 4.5 % / 1.2 % components above are both ^4I_15/2, and that is NOT a misprint: they are two
    #   different antisymmetric 4f^9 states carrying the same total term, distinguished by their seniority.
    #   LSjj.shortString labels a level by its open shells and their accumulated terms, which separates CSFs
    #   that differ in the COUPLING (510 of 541 levels across these branches were ambiguous before 18-Aug and
    #   200 remain), but it does not print the seniority, which is what these two would need.  CsfNR carries
    #   the quasispin QQ and the classification number w, so the information is available if it is ever wanted.
    grid = Radial.Grid(Radial.Grid(false), rnt=2.0e-6, h=5.0e-2, hp=1.0e-2, rbox=8.0)
    wa   = Atomic.Computation(Atomic.Computation(), name="jj-LS: Dy^3+ 4f^9 ground level", grid=grid,
                              nuclearModel=Nuclear.Model(66.0), configs=[Configuration("[Xe] 4f^9")],
                              asfSettings=AsfSettings(AsfSettings(); jjLS=LSjjSettings(true, 0.05, 0.001, LevelSelection())) )
    wb = perform(wa)
    #
end
