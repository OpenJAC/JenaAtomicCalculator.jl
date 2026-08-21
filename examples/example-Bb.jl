
#
println("Bb) Apply & test the WeakInteractionMoment module: the P-odd and P,T-odd one-electron amplitudes")
println("    -- nuclear weak charge, anapole moment and Schiff moment -- between two bound-state levels.")
println("    These are the BARE matrix elements; every observable built from them (a parity-non-conserving")
println("    E1 amplitude, an EDM enhancement factor) needs a sum over intermediate states as well and")
println("    belongs to the property side, not here.")
println("    The module was rebuilt on 21-Aug-2026: until then all three functions returned fabricated")
println("    constants (1+2i through a stub, and a hard-wired 3+3i), which is why every branch below is")
println("    built on a check that a stub could not have satisfied.")

if  true
    # Last successful:  21-Aug-2026
    #   [PROVENANCE: produced on Julia 1.10.9, with a Manifest re-resolved for 1.10 because the checkout's Manifest had been
    #    resolved under 1.12.6 on the maintainer machine.  Re-run there to confirm.]
    # Branch a: THE SELECTION RULES, MADE VISIBLE.  The operators are P-odd, so each must vanish EXACTLY between levels
    #   of the same parity; and the weak charge, being the pseudoscalar rho(r) gamma_5, is of rank ZERO and must in
    #   addition vanish unless J_f = J_i.  The predecessor of this module enforced neither, and its approved test file
    #   asserted non-zero values for three same-parity pairs.
    #
    #   THIS BRANCH EARNED ITS KEEP ON THE DAY IT WAS WRITTEN.  Its anapole column originally came out identically 0.0
    #   for all sixteen pairs -- a pattern that reads like a selection rule and is not one -- which is how a wrong
    #   angular structure in `anapoleAmplitude` was caught.  See branch f, which now covers the anapole in detail.
    #
    #   The test bed is deliberately small and has both parities and two J values: a one-electron ion at Z = 55 in the
    #   BARE nuclear field, giving 1s_1/2+, 2s_1/2+, 2p_1/2- and 2p_3/2-.  Basics.NuclearField() is used rather than the
    #   default DFS so that a one-electron system does not acquire a self-interaction.
    #
    #   RESULT, 21-Aug-2026:  16 ordered pairs, of which 8 have equal parity.  All 8 give exactly 0.0 for all three
    #   operators.  Of the 8 opposite-parity pairs, the weak charge is non-zero for exactly the 4 with J_f = J_i and
    #   exactly zero for the 4 with J_f != J_i, while the rank-1 Schiff and anapole amplitudes survive the J change.
    #   Every zero here is a hard 0.0, not a small number.
    setDefaults("print summary: open", "zzz-WeakInteractionMoment.sum")
    #
    grid   = Radial.Grid(Radial.Grid(false), rnt = 4.0e-7, h = 3.0e-2, hp = 1.0e-2, rbox = 6.0)
    setDefaults("standard grid", grid)
    nModel = Nuclear.Model(55.)
    asfB   = AsfSettings(AsfSettings(); scField = Basics.NuclearField())
    comp   = Atomic.Computation(Atomic.Computation(), name="one-electron test bed, Z=55", grid=grid, nuclearModel=nModel,
                                configs=[Configuration("1s"), Configuration("2s"), Configuration("2p")], asfSettings=asfB)
    mp     = perform(comp; output=true)["multiplet:"]
    #
    println("\n\n  Selection rules of the three P-odd operators (Z = 55, bare nuclear field)\n")
    println("   f <- i     J^P f     J^P i    same par?   dJ?      weak charge        Schiff           anapole")
    nSame = 0;   nViol = 0;   nWeak = 0
    for  a in mp.levels,  b in mp.levels
        wc = WeakInteractionMoment.weakChargeAmplitude(a, b, nModel, grid)
        sm = WeakInteractionMoment.schiffMomentAmplitude(a, b, nModel, grid)
        an = WeakInteractionMoment.anapoleAmplitude(a, b, nModel, grid)
        samePar = (a.parity == b.parity);    sameJ = (a.J == b.J)
        if  samePar
            global nSame += 1
            if  wc != 0. || sm != 0. || an != 0.    global nViol += 1    end
        end
        if  wc != 0.    global nWeak += 1    end
        println("   " * rpad("$(a.index) <- $(b.index)", 10) * rpad(string(LevelSymmetry(a.J,a.parity)), 10) *
                rpad(string(LevelSymmetry(b.J,b.parity)), 9) * rpad(samePar ? "yes" : "no", 12) *
                rpad(sameJ ? "no" : "yes", 9) * rpad(string(round(abs(wc), sigdigits=4)), 18) *
                rpad(string(round(abs(sm), sigdigits=4)), 17) * string(round(abs(an), sigdigits=4)))
    end
    println("\n  same-parity pairs: $nSame,  of which non-vanishing: $nViol   (must be 0)")
    println("  pairs with a non-zero weak charge: $nWeak   (must be the opposite-parity, equal-J ones only)")
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  21-Aug-2026
    #   [PROVENANCE: as branch a.]
    # Branch b: REALITY AND HERMITICITY -- two exact structural facts that no stub can imitate, and the sharpest cheap
    #   check in this file.
    #
    #   (i)  gamma_5 is ANTISYMMETRIC in the large and small Dirac components, so between real radial orbitals the
    #        weak-charge amplitude is PURELY IMAGINARY.  The Schiff operator is a function of position and mixes nothing,
    #        so its amplitude is PURELY REAL.  The old module returned 1+2i for both, i.e. neither.
    #   (ii) H_W is Hermitian, so exchanging the two levels conjugates the reduced matrix element.  For a purely
    #        imaginary rank-0 quantity that means EXACT ANTISYMMETRY, <f||H_W||i> = -<i||H_W||f>.  The rank-1 Schiff
    #        amplitude instead picks up (-1)^(J_f - J_i), so it is symmetric between two J = 1/2 levels and
    #        antisymmetric between J = 1/2 and J = 3/2 -- a sign pattern that is worth reading off, because it is not
    #        something one would put in by hand.
    #
    #   RESULT, 21-Aug-2026:  re/|amp| = 0.0 for every weak-charge amplitude and im/|amp| = 0.0 for every Schiff
    #   amplitude, both exactly.  The antisymmetry sum |<f||H_W||i> + <i||H_W||f>| is 0.000e+00 for every pair.  The
    #   Schiff sign pattern is as predicted: 1<->2 (both J=1/2) symmetric, 1<->4 (J=1/2 against J=3/2) antisymmetric.
    setDefaults("print summary: open", "zzz-WeakInteractionMoment.sum")
    #
    grid   = Radial.Grid(Radial.Grid(false), rnt = 4.0e-7, h = 3.0e-2, hp = 1.0e-2, rbox = 6.0)
    setDefaults("standard grid", grid)
    nModel = Nuclear.Model(55.)
    asfB   = AsfSettings(AsfSettings(); scField = Basics.NuclearField())
    comp   = Atomic.Computation(Atomic.Computation(), name="one-electron test bed, Z=55", grid=grid, nuclearModel=nModel,
                                configs=[Configuration("1s"), Configuration("2s"), Configuration("2p")], asfSettings=asfB)
    mp     = perform(comp; output=true)["multiplet:"]
    #
    println("\n\n  (i) reality:  weak charge must be purely imaginary, Schiff purely real\n")
    println("   f <- i      Re(weak)/|weak|     Im(weak)          Im(Schiff)/|Schiff|   Re(Schiff)")
    for  a in mp.levels,  b in mp.levels
        wc = WeakInteractionMoment.weakChargeAmplitude(a, b, nModel, grid)
        sm = WeakInteractionMoment.schiffMomentAmplitude(a, b, nModel, grid)
        (abs(wc) == 0. && abs(sm) == 0.)  &&  continue
        rw = abs(wc) > 0. ? abs(real(wc))/abs(wc) : NaN
        is = abs(sm) > 0. ? abs(imag(sm))/abs(sm) : NaN
        println("   " * rpad("$(a.index) <- $(b.index)", 11) * rpad(string(rw), 20) *
                rpad(string(round(imag(wc), sigdigits=7)), 18) * rpad(string(is), 22) *
                string(round(real(sm), sigdigits=7)))
    end
    #
    println("\n  (ii) Hermiticity under exchange of the two levels\n")
    println("   pair        weak: |<f|H|i> + <i|H|f>|     Schiff: <f|H|i>/<i|H|f>   expected (-1)^(Jf-Ji)")
    for  a in mp.levels,  b in mp.levels
        a.index >= b.index  &&  continue
        wfi = WeakInteractionMoment.weakChargeAmplitude(a, b, nModel, grid)
        wif = WeakInteractionMoment.weakChargeAmplitude(b, a, nModel, grid)
        sfi = WeakInteractionMoment.schiffMomentAmplitude(a, b, nModel, grid)
        sif = WeakInteractionMoment.schiffMomentAmplitude(b, a, nModel, grid)
        (abs(wfi) == 0. && abs(sfi) == 0.)  &&  continue
        sRatio = abs(sif) > 0. ? real(sfi)/real(sif) : NaN
        expect = (-1.0)^Int(round((Basics.twice(a.J) - Basics.twice(b.J))/2))
        println("   " * rpad("$(a.index)<->$(b.index)", 12) * rpad(string(abs(wfi + wif)), 28) *
                rpad(string(round(sRatio, digits=6)), 24) * string(expect))
    end
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  21-Aug-2026
    #   [PROVENANCE: as branch a.]
    # Branch c: THE RANK-0 NORMALIZATION, measured rather than assumed.  This is the branch to re-run first if anything
    #   about this module is ever doubted, because it is the one that caught a real bug during construction.
    #
    #   `SpinAngular.computeCoefficients` does NOT return its coefficients in a single normalization.  For rank >= 1
    #   (computeCoefficientsNonScalar) the contraction gives the reduced matrix element up to sqrt(2J_f+1); for rank 0
    #   (computeCoefficientsScalar) it gives the ORDINARY matrix element, as it must, since its principal client is the
    #   one-body Hamiltonian -- so the reduced one follows as (2J_f+1) times it.  The weak charge is the rank-0 case.
    #   The first implementation applied sqrt(2J_f+1) to both and came out low by exactly 1/sqrt(2).
    #
    #   For a ONE-electron system the many-electron reduced matrix element IS the one-electron one, so it can be built
    #   independently from the orbitals and the two compared with nothing assumed:
    #        <f||H_W||i>  =  sqrt(2j+1) (G_F/2sqrt2) Q_W  i  INT rho(r)[P_f Q_i - Q_f P_i] dr
    #
    #   RESULT, 21-Aug-2026, Z = 55, 2p_1/2 <- 1s_1/2:  module 6.352592373140062e-9 i against the hand-built
    #   6.3525923731400615e-9 i, ratio 1.0 + 0.0i.  Before the fix the same ratio read 0.7071067811865477 = 1/sqrt(2).
    setDefaults("print summary: open", "zzz-WeakInteractionMoment.sum")
    #
    grid   = Radial.Grid(Radial.Grid(false), rnt = 4.0e-7, h = 3.0e-2, hp = 1.0e-2, rbox = 6.0)
    setDefaults("standard grid", grid)
    nModel = Nuclear.Model(55.)
    asfB   = AsfSettings(AsfSettings(); scField = Basics.NuclearField())
    mp     = SelfConsistent.performSCF([Configuration("1s"), Configuration("2p")], nModel, grid, asfB; printout=false)
    #
    s12 = mp.levels[argmin([l.energy for l in mp.levels])]
    p12 = nothing
    for  l in mp.levels    if  l.parity == Basics.minus && Basics.twice(l.J) == 1    global p12 = l    end    end
    #
    amp    = WeakInteractionMoment.weakChargeAmplitude(p12, s12, nModel, grid)
    rho    = WeakInteractionMoment.nuclearDensity(nModel, grid)
    subS   = first(filter(sh -> sh.kappa == -1, collect(keys(s12.basis.orbitals))))
    subP   = first(filter(sh -> sh.kappa == +1, collect(keys(p12.basis.orbitals))))
    radial = WeakInteractionMoment.radialIntegralPQminus(rho, p12.basis.orbitals[subP], s12.basis.orbitals[subS], grid)
    QW     = WeakInteractionMoment.weakCharge(nModel)
    exact  = sqrt(Basics.twice(p12.J) + 1.0) * WeakInteractionMoment.GF/(2sqrt(2.0)) * QW * im * radial
    #
    println("\n\n  Rank-0 normalization, measured against the exact one-electron reduced matrix element\n")
    println("     Q_W(Z=55, A=125)                    = $QW")
    println("     INT rho [P_p Q_s - Q_p P_s] dr       = $radial")
    println("     module   <2p_1/2||H_W||1s_1/2>       = $amp")
    println("     built by hand from the orbitals      = $exact")
    println("     RATIO (must be exactly 1)            = $(amp/exact)")
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  21-Aug-2026
    #   [PROVENANCE: as branch a.]
    # Branch d: THE NUCLEAR MODEL IS NOT A DETAIL HERE.  These amplitudes are nuclear-overlap quantities: the operator
    #   is proportional to rho(r), which is non-zero only inside the nucleus, so the whole matrix element is fixed by the
    #   electron density over a region of about 1e-5 a.u.  How that region is modelled therefore matters far more than it
    #   does for an ordinary multipole, and this branch measures by how much rather than asserting that it is small.
    #
    #   Three things are varied: the shape at fixed rms radius (uniform sphere against two-parameter Fermi), the rms
    #   radius itself at fixed shape, and the point-nucleus limit -- which is REFUSED rather than approximated, since a
    #   delta density would make every amplitude depend on the orbitals exactly at the origin, where a finite radial grid
    #   says least.  That refusal is itself checked here.
    #
    #   RESULT, 21-Aug-2026, Z = 55, 2p_1/2 <- 1s_1/2:  uniform and Fermi differ by about 2 % at the same rms radius,
    #   and the amplitude falls monotonically as the nucleus is made larger -- both as expected for a quantity weighted by
    #   the electron density near r = 0, and neither of them negligible.  PointNucleus() raises, with an explanation.
    setDefaults("print summary: open", "zzz-WeakInteractionMoment.sum")
    #
    grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-7, h = 3.0e-2, hp = 1.0e-2, rbox = 6.0)
    setDefaults("standard grid", grid)
    asfB = AsfSettings(AsfSettings(); scField = Basics.NuclearField())
    #
    println("\n\n  Dependence on the nuclear model, Z = 55, amplitude <2p_1/2||H_W||1s_1/2>\n")
    println("   model            R_rms [fm]     |amplitude|            relative to the first row")
    ref = 0.0
    for  (mdl, rrms)  in  [ (Nuclear.FermiNucleus(),   4.75), (Nuclear.UniformNucleus(), 4.75),
                            (Nuclear.FermiNucleus(),   4.00), (Nuclear.FermiNucleus(),   5.50) ]
        nm  = Nuclear.Model(55., mdl, 125.125, rrms, AngularJ64(1//2), 0.0, 0.0, 0.0)
        m   = SelfConsistent.performSCF([Configuration("1s"), Configuration("2p")], nm, grid, asfB; printout=false)
        ss  = m.levels[argmin([l.energy for l in m.levels])]
        pp  = nothing
        for  l in m.levels    if  l.parity == Basics.minus && Basics.twice(l.J) == 1    pp = l    end    end
        a   = abs(WeakInteractionMoment.weakChargeAmplitude(pp, ss, nm, grid))
        if  ref == 0.0    global ref = a    end
        println("   " * rpad(string(typeof(mdl).name.name), 17) * rpad(string(rrms), 15) *
                rpad(string(round(a, sigdigits=8)), 24) * string(round(a/ref, digits=6)))
    end
    #
    println("\n  and the point nucleus, which must be refused rather than approximated:")
    try
        nmPt = Nuclear.Model(55., Nuclear.PointNucleus(), 125.125, 4.75, AngularJ64(1//2), 0.0, 0.0, 0.0)
        WeakInteractionMoment.nuclearDensity(nmPt, grid)
        println("     NO ERROR RAISED -- that is a defect")
    catch e
        println("     raised, as intended:  ", first(sprint(showerror, e), 160), " ...")
    end
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  21-Aug-2026
    #   [PROVENANCE: as branch a.]
    # Branch e: THE Z-DEPENDENCE -- CLOSED, and the law is Z^4, not Z^3.
    #
    #   This branch was written to check the amplitudes against the Bouchiat "Z^3 law" and sat undated for a day
    #   because the measured exponent came out at 4.4 rising to 7.3, nothing like 3.  Two things were wrong with the
    #   expectation, and neither was the code:
    #
    #   (i)  Z^3 IS THE LAW FOR A NEUTRAL ATOM'S VALENCE ELECTRON, where the normalization of the valence wave function
    #        near the nucleus combines with Q_W ~ -N ~ -Z in a particular way.  For a hydrogen-like ion the contact
    #        limit gives Q_W * Z^3 * (Z alpha): one factor Z^(3/2) from each of the two wave functions at the origin,
    #        and one more power of Z from the small component, which is what carries gamma_5 between large and small.
    #        Dividing out the trivial Q_W(Z) therefore leaves Z^4.
    #   (ii) The earlier scan let the nuclear radius grow with Z, since Nuclear.Model(Z) sets R from an empirical
    #        A(Z).  The density normalization 1/R^3 then carries its own Z-dependence into the integral.  Here the
    #        radius is held FIXED so that only the electronic Z-dependence is measured.
    #
    #   RESULT, 21-Aug-2026, R_rms = 3.5 fm throughout, 2s_1/2 <-> 2p_1/2:  the local exponent is 4.042 at Z = 6-8,
    #   4.068 at 8-10, 4.117 at 10-14 and 4.224 at 14-20; least squares 4.1198 over the range.  It approaches 4 from
    #   above as Z falls, and the excess is the relativistic enhancement, which must grow with Z alpha.  So the law is
    #   Z^4 and the code obeys it; for contrast, the ORIGINAL scan with a growing radius gave 4.418 to 7.278.
    setDefaults("print summary: open", "zzz-WeakInteractionMoment.sum")
    #
    asfB = AsfSettings(AsfSettings(); scField = Basics.NuclearField())
    Zs   = [6.0, 8.0, 10.0, 14.0, 20.0]
    vals = Float64[]
    println("\n\n  Z-dependence at FIXED R_rms = 3.5 fm, |<2p_1/2||H_W||2s_1/2>| / |Q_W|\n")
    println("      Z      |amp| / |Q_W|             local exponent   (expected 4)")
    for  Zi in Zs
        gr  = Radial.Grid(Radial.Grid(false), rnt = 2.0e-7, h = 3.0e-2, hp = 8.0e-3, rbox = 40.0)
        setDefaults("standard grid", gr)
        nmi = Nuclear.Model(Zi, Nuclear.FermiNucleus(), 2.0*Zi, 3.5, AngularJ64(1//2), 0.0, 0.0, 0.0)
        m   = SelfConsistent.performSCF([Configuration("2s"), Configuration("2p")], nmi, gr, asfB; printout=false)
        sL  = nothing;   pL = nothing
        for  l in m.levels
            if      l.parity == Basics.plus  && Basics.twice(l.J) == 1    sL = l
            elseif  l.parity == Basics.minus && Basics.twice(l.J) == 1    pL = l   end
        end
        v = abs(WeakInteractionMoment.weakChargeAmplitude(pL, sL, nmi, gr)) / abs(WeakInteractionMoment.weakCharge(nmi))
        push!(vals, v)
        sa = length(vals) > 1 ?
             string(round((log(vals[end])-log(vals[end-1]))/(log(Zi)-log(Zs[length(vals)-1])), digits=4)) : "--"
        println("   " * rpad(string(Zi), 8) * rpad(string(round(v, sigdigits=8)), 26) * sa)
    end
    nn = length(Zs);  lx = log.(Zs);  ly = log.(vals)
    println("\n   least-squares exponent = " *
            string(round((nn*sum(lx.*ly) - sum(lx)*sum(ly))/(nn*sum(lx.^2) - sum(lx)^2), digits=4)) * "   (expected 4)")
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  21-Aug-2026
    #   [PROVENANCE: as branch a.]
    # Branch f: THE ANAPOLE, and the mistake that only Hermiticity could catch.  Written last, because the anapole was
    #   the hardest of the three and got its angular structure wrong TWICE before this branch was dated.
    #
    #   The operator is (G_F/sqrt2) kappa_a alpha.I rho(r), rank 1 and P-odd.  alpha = [[0, sigma], [sigma, 0]] connects
    #   the LARGE component of one orbital with the SMALL component of the other, so with Omega_(-kappa) = -(sigma.rhat)
    #   Omega_(kappa) the one-electron reduced matrix element is
    #
    #       <k_a||alpha rho||k_b> = i [ <Om(k_a)||sigma||Om(-k_b)> INT rho P_a Q_b dr
    #                                 - <Om(-k_a)||sigma||Om(k_b)> INT rho Q_a P_b dr ]
    #
    #   FIRST MISTAKE, caught by branch a: giving alpha the MAGNETIC-multipole template
    #   <-k_a||C^1||k_b> (k_a + k_b).  That is parity-EVEN, so it survived for p_1/2 <- p_3/2 and vanished for BOTH
    #   s_1/2 <- p_1/2 and s_1/2 <- p_3/2 -- the anapole column came out identically zero for all sixteen pairs.
    #
    #   SECOND MISTAKE, caught ONLY by Hermiticity: combining the two cross terms above into the symmetric
    #   P_a Q_b + Q_a P_b, on the assumption that their angular factors differ merely by a sign.  That holds for gamma_5,
    #   whose two factors are both delta(k_a, -k_b), and fails here: for s_1/2 <-> p_1/2 they are <s||sigma||s> = sqrt(6)
    #   and <p_1/2||sigma||p_1/2> = -sqrt(6)/3, a ratio of -3.  The wrong version passed the selection rules, passed the
    #   reality check, and gave amplitudes of an entirely plausible size; exchanging the two levels multiplied the
    #   amplitude by -3 instead of -1, and nothing else in this file would have noticed.
    #
    #   RESULT, 21-Aug-2026, Z = 55:
    #     - reachability, from the sigma reduced matrix element: s_1/2 <-> p_1/2 YES (2.44949), p_1/2 <-> d_3/2 YES
    #       (-2.30940), p_3/2 <-> d_3/2 YES (2.58199); the parity-CONSERVING p_1/2 <-> p_3/2 is 0, as it must be;
    #     - the module against the exact one-electron reduced matrix element built from the orbitals: ratio
    #       1.0000000000000002;
    #     - all 8 same-parity pairs exactly 0.0, every amplitude purely imaginary (re/|amp| = 0.0);
    #     - Hermiticity exact: |<f|H|i> + <i|H|f>| = 0.000e+00 for the J_f = J_i pairs and |<f|H|i> - <i|H|f>| = 0.000e+00
    #       for the J = 1/2 <-> 3/2 pairs.
    setDefaults("print summary: open", "zzz-WeakInteractionMoment.sum")
    #
    grid   = Radial.Grid(Radial.Grid(false), rnt = 4.0e-7, h = 3.0e-2, hp = 1.0e-2, rbox = 6.0)
    setDefaults("standard grid", grid)
    nModel = Nuclear.Model(55.)
    asfB   = AsfSettings(AsfSettings(); scField = Basics.NuclearField())
    #
    println("\n\n  (i) which kappa pairs the anapole can reach, from <Om(k_a)||sigma||Om(-k_b)>\n")
    println("   a  <- b            l(k_a)  l(-k_b)   sigma red. m.e.    reachable?")
    for  (ka, kb, nam)  in  [(-1, 1,"s_1/2 <- p_1/2"), (-1,-2,"s_1/2 <- p_3/2"), ( 1,-2,"p_1/2 <- p_3/2"),
                             ( 1, 2,"p_1/2 <- d_3/2"), (-2, 2,"p_3/2 <- d_3/2"), (-2,-3,"p_3/2 <- d_5/2")]
        la = Basics.subshell_l(Subshell(9, ka));    lb = Basics.subshell_l(Subshell(9, -kb))
        sg = WeakInteractionMoment.sigmaReducedMe(ka, -kb)
        println("   " * rpad(nam, 19) * rpad(string(la), 8) * rpad(string(lb), 10) *
                rpad(string(round(sg, digits=6)), 19) * (sg != 0. ? "YES" : "no"))
    end
    println("\n   s_1/2 <-> p_1/2 MUST be reachable: that mixing is how an anapole moment enters a measured")
    println("   PNC transition.  The withdrawn magnetic template gave exactly 0 there.")
    #
    mp   = SelfConsistent.performSCF([Configuration("1s"), Configuration("2p")], nModel, grid, asfB; printout=false)
    s12  = mp.levels[argmin([l.energy for l in mp.levels])]
    p12  = nothing
    for  l in mp.levels    if  l.parity == Basics.minus && Basics.twice(l.J) == 1    global p12 = l    end    end
    rho  = WeakInteractionMoment.nuclearDensity(nModel, grid)
    subS = first(filter(sh -> sh.kappa == -1, collect(keys(s12.basis.orbitals))))
    subP = first(filter(sh -> sh.kappa == +1, collect(keys(p12.basis.orbitals))))
    oP   = p12.basis.orbitals[subP];    oS = s12.basis.orbitals[subS]
    mod  = WeakInteractionMoment.anapoleAmplitude(p12, s12, nModel, grid)
    sgA  = WeakInteractionMoment.sigmaReducedMe( oP.subshell.kappa, -oS.subshell.kappa)
    sgB  = WeakInteractionMoment.sigmaReducedMe(-oP.subshell.kappa,  oS.subshell.kappa)
    r1   = WeakInteractionMoment.radialIntegralPQ(rho, oP, oS, grid)
    r2   = WeakInteractionMoment.radialIntegralPQ(rho, oS, oP, grid)
    hand = WeakInteractionMoment.GF/sqrt(2.0) * im * (sgA*r1 - sgB*r2)
    println("\n  (ii) the module against the exact one-electron reduced matrix element\n")
    println("     subshells used                  : a = $subP,  b = $subS")
    println("     <Om(k_a)||sigma||Om(-k_b)> = $sgA   weights INT rho P_a Q_b = $r1")
    println("     <Om(-k_a)||sigma||Om(k_b)> = $sgB   weights INT rho Q_a P_b = $r2")
    println("     module                          = $mod")
    println("     exact                           = $hand")
    println("     RATIO (must be exactly 1)       = $(mod/hand)")
    #
    println("\n  (iii) Hermiticity -- the check that caught the second mistake\n")
    mp2 = SelfConsistent.performSCF([Configuration("1s"), Configuration("2s"), Configuration("2p")], nModel, grid, asfB; printout=false)
    println("   pair     J_f vs J_i    |<f|H|i> + <i|H|f>|     |<f|H|i> - <i|H|f>|    which must vanish")
    for  a in mp2.levels,  b in mp2.levels
        a.index >= b.index  &&  continue
        x = WeakInteractionMoment.anapoleAmplitude(a, b, nModel, grid)
        y = WeakInteractionMoment.anapoleAmplitude(b, a, nModel, grid)
        abs(x) == 0.  &&  continue
        println("   " * rpad("$(a.index)<->$(b.index)", 9) * rpad(a.J == b.J ? "equal" : "differ", 14) *
                rpad(string(round(abs(x+y), sigdigits=4)), 24) * rpad(string(round(abs(x-y), sigdigits=4)), 23) *
                (a.J == b.J ? "the sum" : "the difference"))
    end
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  21-Aug-2026
    #   [PROVENANCE: as branch a.]
    # Branch g: THE ABSOLUTE SCALE, which until this branch existed was the one thing this module had NOT been tested
    #   on.  Branches a to f are all STRUCTURAL -- selection rules, phases, reality, Hermiticity, and one internal
    #   normalization -- and a module can pass every one of them while being wrong by a constant.  G_F, Q_W, the
    #   1/(2 sqrt 2) of the Hamiltonian and the normalization of rho(r) are each individually defensible and were,
    #   until now, collectively unchecked.
    #
    #   The anchor is the closed-form hydrogenic result for the weak matrix element between ns_1/2 and np_1/2 of the
    #   SAME n, in the contact limit:
    #
    #        M_n  =  i G_F Q_W / (4 pi sqrt2 c) * sqrt(n^2 - 1) / n^4
    #
    #   (cf. the parity-non-conservation-in-hydrogen literature, e.g. arXiv:2605.10321 and Phys. Rev. A 18 (1978) 2421).
    #   It is independent of everything in this module: an analytic result for hydrogenic orbitals, against a numerical
    #   integral over B-splines.  The two share only G_F and Q_W, which are standard constants, so the comparison tests
    #   the grouping 1/(2 sqrt2) against 1/(4 pi sqrt2 c), the normalization of the density, and the radial integral.
    #
    #   TWO CHECKS, and the first needs no constants at all:
    #
    #   (i)  THE n-DEPENDENCE, sqrt(n^2-1)/n^4, in which every prefactor cancels.  Measured at Z = 10 for n = 2,3,4,5,
    #        ratios to n = 2:  0.322642 against 0.322567 (+0.023 %), 0.139770 against 0.139754 (+0.011 %), 0.072404
    #        against 0.072408 (-0.005 %).  That is the orbitals, the radial integral and the density all together, to
    #        better than a part in 4000, with nothing adjustable.
    #
    #   (ii) THE ABSOLUTE VALUE, on hydrogen itself, 2s_1/2 <-> 2p_1/2.  The module's ORDINARY matrix element (the
    #        reduced one divided by sqrt(2J+1) = sqrt 2) comes to 7.426076e-20 a.u. against the closed form's
    #        7.421407e-20 a.u., a ratio of 1.00063.  The residual 0.06 % is the finite nucleus: the closed form is the
    #        point-nucleus contact limit while the module integrates over a uniform sphere of R_rms = 0.8783 fm, and a
    #        correction of that order is exactly what a finite size should cost at Z = 1.
    #
    #   ONE CAUTION ON THE LITERATURE FORM.  It is quoted in some places with 1/c^2 rather than 1/c.  The module
    #   disagrees with that version by a factor of 137.12, which is c to four digits -- so the c^2 version carries one
    #   power of c too many for the convention used here, in which the single 1/c comes from the small component that
    #   gamma_5 reaches.  This branch prints both so that the reader can see which is which rather than take it on
    #   trust.
    setDefaults("print summary: open", "zzz-WeakInteractionMoment.sum")
    #
    asfB = AsfSettings(AsfSettings(); scField = Basics.NuclearField())
    shape(n) = sqrt(n^2 - 1.0) / n^4
    #
    println("\n\n  (i) the n-dependence at Z = 10, expected sqrt(n^2-1)/n^4 -- every prefactor cancels\n")
    println("      n     |<np_1/2||H_W||ns_1/2>|      ratio to n=2     expected       dev [%]")
    gr10 = Radial.Grid(Radial.Grid(false), rnt = 2.0e-7, h = 3.0e-2, hp = 8.0e-3, rbox = 40.0)
    setDefaults("standard grid", gr10)
    nm10 = Nuclear.Model(10.)
    vals = Float64[]
    for  n = 2:5
        m  = SelfConsistent.performSCF([Configuration("$(n)s"), Configuration("$(n)p")], nm10, gr10, asfB; printout=false)
        sL = nothing;   pL = nothing
        for  l in m.levels
            if      l.parity == Basics.plus  && Basics.twice(l.J) == 1    sL = l
            elseif  l.parity == Basics.minus && Basics.twice(l.J) == 1    pL = l   end
        end
        v = abs(WeakInteractionMoment.weakChargeAmplitude(pL, sL, nm10, gr10))
        push!(vals, v)
        r = v/vals[1];    e = shape(n)/shape(2)
        println("   " * rpad(string(n), 8) * rpad(string(round(v, sigdigits=9)), 26) *
                rpad(string(round(r, digits=6)), 17) * rpad(string(round(e, digits=6)), 15) *
                string(round(100*(r/e - 1), digits=3)))
    end
    #
    println("\n  (ii) the absolute value, hydrogen 2s_1/2 <-> 2p_1/2\n")
    grH  = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 3.0e-2, hp = 1.0e-2, rbox = 60.0)
    setDefaults("standard grid", grH)
    nmH  = Nuclear.Model(1., UniformNucleus(), 1.00794, 0.8783, AngularJ64(1//2), 2.7928473, 0.0, 0.0)
    mH   = SelfConsistent.performSCF([Configuration("2s"), Configuration("2p")], nmH, grH, asfB; printout=false)
    sH = nothing;   pH = nothing
    for  l in mH.levels
        if      l.parity == Basics.plus  && Basics.twice(l.J) == 1    global sH = l
        elseif  l.parity == Basics.minus && Basics.twice(l.J) == 1    global pH = l   end
    end
    red  = WeakInteractionMoment.weakChargeAmplitude(pH, sH, nmH, grH)
    ord  = red / sqrt(2.0)
    QW   = WeakInteractionMoment.weakCharge(nmH)
    wc   = Defaults.getDefaults("speed of light: c")
    lit1 = WeakInteractionMoment.GF*QW/(4pi*sqrt(2.0)*wc)    * sqrt(3.0)/16
    lit2 = WeakInteractionMoment.GF*QW/(4pi*sqrt(2.0)*wc^2)  * sqrt(3.0)/16
    println("     Q_W(H), tree level                 = $QW")
    println("     module, reduced                    = $red")
    println("     module, ordinary  = reduced/sqrt2  = $ord")
    println("     closed form with 1/c               = $lit1")
    println("     RATIO ordinary / closed form       = $(abs(ord)/lit1)      <-- the anchor")
    println("     closed form with 1/c^2             = $lit2")
    println("     ratio against that version         = $(abs(ord)/lit2)      <-- i.e. c, so 1/c^2 is one power too many")
    #
    setDefaults("print summary: close", "")
    #
end
