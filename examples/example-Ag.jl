#
println("Ag) Apply & test the jj-LS transformation of levels from a given multiplet.")

if  false
    # Last successful:  30-Jul-2026
    # Branch a: ONE open shell -- baseline, uses the existing (trusted) LSjj.OneOpenShell code path;
    #   no recoupling coefficient needed at all.
    wa = Atomic.Computation(Atomic.Computation(), name="jj-LS: one open shell", grid=Radial.Grid(true), nuclearModel=Nuclear.Model(10.),
                            configs=[Configuration("[Ne] 3p")],
                            asfSettings=AsfSettings(AsfSettings(); jjLS=LSjjSettings(true, 0.05, 0.1, LevelSelection())) )
    wb = perform(wa)
    #
elseif  false
    # Last successful:  30-Jul-2026
    # Branch b: TWO open shells -- exercises LSjj.GeneralOpenShells(2) (Gaigalas, Atoms 14, 20 (2026),
    #   Eq. 17), routed to by Basics.extractFromConfiguration(Basics.OpenShellNumber(),...). Verified
    #   to reproduce the formerly-trusted, now-retired LSjj.TwoOpenShells code exactly (max|old-new| =
    #   2.2e-16 across all CSFs of this configuration) before that code was removed, so this is also the
    #   regression reference the 3-, 4- and 5-open-shell branches below build confidence from.
    wa = Atomic.Computation(Atomic.Computation(), name="jj-LS: two open shells", grid=Radial.Grid(true), nuclearModel=Nuclear.Model(10.),
                            configs=[Configuration("[Ne] 3s 3p")],
                            asfSettings=AsfSettings(AsfSettings(); jjLS=LSjjSettings(true, 0.05, 0.1, LevelSelection())) )
    wb = perform(wa)
    #
elseif  false
    # Last successful:  30-Jul-2026
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
                            asfSettings=AsfSettings(AsfSettings(); jjLS=LSjjSettings(true, 0.05, 0.1, LevelSelection())) )
    wb = perform(wa)
    #
elseif  false
    # Last successful:  30-Jul-2026
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
                            asfSettings=AsfSettings(AsfSettings(); jjLS=LSjjSettings(true, 0.05, 0.1, LevelSelection())) )
    wb = perform(wa)
    #
elseif  true
    # Last successful:  30-Jul-2026
    # Branch e: FIVE open shells -- exercises LSjj.GeneralOpenShells(5), representative of the kind of
    #   configuration that can arise from inner-shell holes in a cascade (multiple simultaneously open
    #   subshells across several n). Confirms the fixes from branches c/d generalize: exactly weight=1.0
    #   for all 249 physical levels of [Ne] 3s 3p 3d 4s 4p.
    wa = Atomic.Computation(Atomic.Computation(), name="jj-LS: five open shells", grid=Radial.Grid(true), nuclearModel=Nuclear.Model(26.),
                            configs=[Configuration("[Ne] 3s 3p 3d 4s 4p")],
                            asfSettings=AsfSettings(AsfSettings(); jjLS=LSjjSettings(true, 0.05, 0.1, LevelSelection())) )
    wb = perform(wa)
    #
end
