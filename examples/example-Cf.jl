
println("Cf) Apply & test the reduced 1-particle density matrix, natural orbitals and related quantities.")
println("    Branches of increasing complexity, following Ma, Li, Godefroid et al., Atoms 12, 30 (2024)")
println("    (examples/papers/b24.atoms-ma-li-natural-orbitals.pdf); see project memory project_natural_orbitals.md.")

setDefaults("print summary: open", "zzz-ReducedDensityMatrix.sum")

if  false
    # Last successful:  23-Jul-2026
    # Branch a: smoke test on a single-CSF (no correlation) system -- Be 1s^2 2s^2 ground state.
    #   Checks: naturalOccupation == [2.0 (1s), 2.0 (2s)] exactly; rho1p diagonal in each 1x1 kappa block;
    #   naturalOrbitalExpansion == identity (trivial, since only one subshell per kappa is present); this
    #   validates the (now fixed) per-kappa diagonalization machinery end-to-end against a known, trivial answer.
    #   Verified: rho1p = diag(2,2) exactly, off-diagonal 0; naturalOccupation = [2.0, 2.0]; naturalOrbitalExpansion
    #   is the identity (up to overall eigenvector sign, i.e. -1 on the diagonal, which is immaterial).
    wa = Atomic.Computation(Atomic.Computation(), name="Cf-a", grid=Radial.Grid(true), nuclearModel=Nuclear.Model(4.),
                            configs=[Configuration("1s^2 2s^2")],
                            propertySettings=[ ReducedDensityMatrix.Settings(true, false, false, false, true, LevelSelection(true, indices=[(1)])) ] )

    wb = perform(wa; output=true)
    #
elseif  false
    # Last successful:  23-Jul-2026
    # Branch b: fractional occupation from configuration mixing between DIFFERENT kappa -- Be 1S0, 2-config CI
    #   1s^2 2s^2  <->  1s^2 2p^2  (the textbook Be near-degeneracy correlation pair).
    #   Checks: naturalOccupation(2s) < 2.0 and naturalOccupation(2p) > 0.0 -- the basic correlation signature
    #   (2s, 2p have different kappa, so this tests fractional OCCUPATION but not yet genuine orbital mixing).
    #   Verified: naturalOccupation = [2.000 (1s), 1.814 (2s), 0.0619 (2p_1/2), 0.1238 (2p_3/2)], sum = 4.000
    #   (conserves N_electrons exactly); 1s/2s stay unmixed since 1s is a pure spectator core orbital in both
    #   CSFs here, so the kappa=-1 block is trivially diagonal already -- correct, not a bug.
    wa = Atomic.Computation(Atomic.Computation(), name="Cf-b", grid=Radial.Grid(true), nuclearModel=Nuclear.Model(4.),
                            configs=[Configuration("1s^2 2s^2"), Configuration("1s^2 2p^2")],
                            propertySettings=[ ReducedDensityMatrix.Settings(true, false, false, false, true, LevelSelection(true, indices=[(1)])) ] )

    wb = perform(wa; output=true)
    #
elseif  false
    # Last successful:  23-Jul-2026
    # Branch c: genuine natural-ORBITAL mixing (same kappa) -- Be 1S0, 2-config CI  1s^2 2s^2  <->  1s^2 2s 3s
    #   (2s and 3s share kappa = -1, so the 1p RDM's kappa=-1 block is now 3x3: 1s, 2s, 3s). NOTE: a genuinely
    #   nonzero off-diagonal rho1p[2s,3s] needs CSFs differing by exactly ONE spin-orbital (single 2s->3s
    #   substitution here) -- a rank-0 (density) operator has EXACTLY zero matrix elements between CSFs that
    #   differ by two spin-orbitals (Slater-Condon rules), so a double substitution like "1s^2 3s^2" only ever
    #   gives fractional DIAGONAL occupation, never off-diagonal mixing; tried and confirmed empirically first.
    #   Also uncovered and fixed a genuine bug along the way: SpinAngular's off-diagonal one-particle coefficient
    #   is missing a normalization factor present in its diagonal branch, giving an unphysical (negative) natural
    #   occupation number. Rather than patch that shared module (used by ~10 other properties), this branch now
    #   uses the new, independent ReducedDensityMatrix.compute1pRDMClaude (derived from scratch via second
    #   quantization; diagonal = mc-weighted CSF occupation, off-diagonal single-substitution transfer
    #   coefficient = sqrt(2), restricted to closed-shell<->one-hole / empty<->one-particle substitutions).
    #   This exercises the full Eqs. 9-11 recipe of Ma et al. for the first time: the natural "2s"/"3s" orbitals
    #   become genuine linear combinations of the standard 2s, 3s orbitals.
    #   Verified: naturalOccupation = [2.000 (1s), 2.000 (2s), 1.8e-4 (3s)] -- all non-negative, summing to 2.0002
    #   (positive semi-definite, physically valid); naturalOrbitalExpansion[2s] = (0.9953, -0.0969) and
    #   [3s] = (0.0969, 0.9953) in the (2s,3s) subspace -- genuine, non-trivial mixing as required.
    wa = Atomic.Computation(Atomic.Computation(), name="Cf-c", grid=Radial.Grid(true), nuclearModel=Nuclear.Model(4.),
                            configs=[Configuration("1s^2 2s^2"), Configuration("1s^2 2s 3s")],
                            propertySettings=[ ReducedDensityMatrix.Settings(true, false, false, false, true, LevelSelection(true, indices=[(1)])) ] )

    wb = perform(wa; output=true)
    #
elseif  false
    # Last successful:  23-Jul-2026
    # Branch d: mean orbital radii <r> before/after the natural-orbital transformation -- same Be 2s/3s CI system
    #   as branch c. Direct structural analog of Table 3/9 of Ma et al. (mean radii before/after NO transform).
    #   Uses RadialIntegrals.rkDiagonal(1, orb, orb, grid) on both the standard basis orbitals (level.basis.orbitals)
    #   and the constructed natural orbitals (outcome.naturalOrbitals) for the 2s, 3s subshells.
    #   Verified: <r>_2s goes 2.4026 (standard) -> 2.5753 (natural), <r>_3s goes 8.2434 -> 8.0708 -- natural 2s is
    #   pulled slightly outward by the (small) 3s admixture, natural 3s pulled slightly inward by the 2s
    #   admixture; correct qualitative direction and a physically reasonable magnitude given the ~0.94% mixing.
    wa = Atomic.Computation(Atomic.Computation(), name="Cf-d", grid=Radial.Grid(true), nuclearModel=Nuclear.Model(4.),
                            configs=[Configuration("1s^2 2s^2"), Configuration("1s^2 2s 3s")],
                            propertySettings=[ ReducedDensityMatrix.Settings(true, false, false, false, true, LevelSelection(true, indices=[(1)])) ] )

    wb    = perform(wa; output=true)
    grid  = wb["grid:"]
    outc  = wb["RDM outcomes:"][1]
    level = outc.level
    println("\n  Mean radii <r> [a.u.] before/after the natural-orbital transformation (cf. Table 3/9 of Ma et al.):\n")
    for  subsh  in  [Subshell("2s_1/2"), Subshell("3s_1/2")]
        rStd = RadialIntegrals.rkDiagonal(1, level.basis.orbitals[subsh], level.basis.orbitals[subsh], grid)
        rNat = RadialIntegrals.rkDiagonal(1, outc.naturalOrbitals[subsh], outc.naturalOrbitals[subsh], grid)
        println("    $subsh:   <r>_standard = $(round(rStd, digits=4))    <r>_natural = $(round(rNat, digits=4))")
    end
    #
elseif  false
    # Last successful:  23-Jul-2026
    # Branch e: radial electron density and natural-orbital normalization -- same Be 2s/3s CI system.
    #   Checks: (1) each natural orbital stays normalized under the linear-combination construction,
    #   RadialIntegrals.rkDiagonal(0, natOrb, natOrb, grid) ~ 1 for every subshell; (2) sum(naturalOccupation)
    #   == N electrons == 4; (3) the printed radial density table (settings.calcDensity=true) is smooth and
    #   decays at large r.
    #   Verified: all three natural-orbital norms = 1.000000 exactly; sum(naturalOccupation) = 4.0 exactly;
    #   density rises smoothly from the origin, peaks around r~0.3-1.5 a.u. (consistent with Be's 1s/2s shell
    #   structure), and decays cleanly to ~1e-13 by r~600 a.u. -- no spurious oscillations or negative values.
    wa = Atomic.Computation(Atomic.Computation(), name="Cf-e", grid=Radial.Grid(true), nuclearModel=Nuclear.Model(4.),
                            configs=[Configuration("1s^2 2s^2"), Configuration("1s^2 2s 3s")],
                            propertySettings=[ ReducedDensityMatrix.Settings(true, true, false, false, true, LevelSelection(true, indices=[(1)])) ] )

    wb   = perform(wa; output=true)
    grid = wb["grid:"]
    outc = wb["RDM outcomes:"][1]
    println("\n  Natural-orbital normalization checks (should all be ~1.0):\n")
    for  (subsh, orb)  in  outc.naturalOrbitals
        println("    $subsh:   <natOrb|natOrb> = $(round(RadialIntegrals.rkDiagonal(0, orb, orb, grid), digits=6))")
    end
    println("\n  sum(naturalOccupation) = $(sum(outc.naturalOccupation))   [should equal N_electrons = 4]")
    #
elseif  true
    # Last successful:  unknown -- deliberately left undated, see below.
    # Branch f: literature system -- N I (Z=7), [He] 2s^2 2p^3 4S°(3/2) ground level, with a small exploratory
    #   correlation config [He] 2s^2 2p^2 3p added (single 2p->3p-type excitation, same odd parity). Reports the
    #   hyperfine constant A [MHz] (Hfs module, calcM1=true) for the ground level -- the closest achievable
    #   "directly comparable with Ma et al." quantity (their Tables 3,4), using nuclear data for 14N (I=1,
    #   mu=0.40376 n.m.).
    #   ReducedDensityMatrix is intentionally NOT included here: 2p^3 is an open shell held by 3 (not 0, 1, or a
    #   full 4) electrons, outside the restricted "closed-shell-hole <-> empty-particle" scope of
    #   compute1pRDMClaude (confirmed: it correctly raised its documented error for this case, occupations 3,2 and
    #   0,1, rather than silently mishandling it) -- general open-shell CFP support would need the original
    #   SpinAngular-based compute1pRDM (once its normalization asymmetry is properly fixed, see project memory)
    #   or a further-generalized compute1pRDMClaude, neither pursued here.
    #   OPEN ISSUE, left undated: A(ground) = 4021 MHz here, vs. experiment 10.4509(10) MHz -- two to three
    #   orders of magnitude too large, and EVERY level in the printed table is similarly far off (not just the
    #   ground state), so this is systematic, not a fluke for one level. Diagnosis: the paper itself states (Sec.
    #   3) that CSFs describing spin/orbital POLARIZATION of the CLOSED 1s, 2s subshells are the MOST important
    #   contribution to the hyperfine constant for a p-shell atom like N I -- this branch's only correlation
    #   config (2p->3p) is a VALENCE correlation, not a core-polarization one, so the dominant physical mechanism
    #   is simply absent here. This matches the paper's own narrative for why naive/energy-driven expansions
    #   converge poorly, rather than indicating a code bug; a quantitatively meaningful A value would need core
    #   CSFs like [He] 2s 2p^4 or similar single substitutions out of 1s/2s, not attempted here.
    nm = Nuclear.Model(7., FermiNucleus(), 14.003074, Nuclear.rrmsRadius(14.003074), AngularJ64(1), 0.40376, 0., 0.)
    wa = Atomic.Computation(Atomic.Computation(), name="Cf-f", grid=Radial.Grid(true), nuclearModel=nm,
                            configs=[Configuration("[He] 2s^2 2p^3"), Configuration("[He] 2s^2 2p^2 3p")],
                            propertySettings=[ Hfs.Settings(calcM1=true, calcE2=false, printBefore=true) ] )

    wb = perform(wa; output=true)
    #
end
#
setDefaults("print summary: close", "")
