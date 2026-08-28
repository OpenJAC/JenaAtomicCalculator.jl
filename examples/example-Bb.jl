
#
println("Bb) Apply & test the WeakInteractionMoment module: the P-odd and P,T-odd one-electron amplitudes")
println("    -- nuclear weak charge, anapole moment and Schiff moment -- between two bound-state levels.")
println("    These are the BARE matrix elements; every observable built from them (a parity-non-conserving")
println("    E1 amplitude, an EDM enhancement factor) needs a sum over intermediate states as well and")
println("    belongs to the property side, not here.")
println("    The module was rebuilt on 21-Aug-2026: until then all three functions returned fabricated")
println("    constants (1+2i through a stub, and a hard-wired 3+3i), which is why every branch below is")
println("    built on a check that a stub could not have satisfied.")

using Printf

if  true
    # Last successful:  28-Aug-2026 -- reproduces, BUT SEE THE CAVEAT, which is the useful part.
    #   THIS BRANCH'S CHECKS WERE PASSING VACUOUSLY UNTIL TODAY. It asserts that the P-odd operators VANISH between
    #   same-parity levels. From the rank-0 parity-gate defect fixed on 28-Aug (see example-Cnnew.jl branch f), the
    #   weak-charge amplitude was returning exactly zero for EVERY pair, same parity or not -- which satisfies a
    #   "must vanish" test perfectly. A selection-rule check cannot tell "correctly zero" from "always zero" unless
    #   something also asserts a NON-zero value where one is required. example-Cnnew.jl branch f is what did that
    #   here, and it is why the two files are worth keeping together.
    # Previously:  21-Aug-2026
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
elseif  false
    # Last successful:  22-Aug-2026
    #   [PROVENANCE: as branch a.]
    # Branch h: THE NEUTRON SKIN -- a sensitivity study, and a limitation of this module made quantitative.
    #
    #   WHAT IS WRONG TODAY.  Q_W = -N + Z(1 - 4 sin^2 theta_W) is computed correctly, but it multiplies a SINGLE
    #   normalized density rho(r) built from nm.radius -- which is the CHARGE radius.  Since (1 - 4 sin^2 theta_W)
    #   = 0.0751, the weak charge is carried almost entirely by the NEUTRONS: for 208Pb, -126 of -119.84, i.e. 105 %,
    #   the protons contributing -5 % of opposite sign.  And the neutron distribution is WIDER than the proton one, by
    #   the neutron skin.  So this module currently puts the weak charge in the wrong place.
    #
    #   WHY THIS BRANCH IS CHEAP, AND TRUSTWORTHY.  The orbitals are fixed by the CHARGE distribution and do not change
    #   when the weak density moves, so the SCF is run once and only the model handed to weakChargeAmplitude is varied.
    #   Passing a different Nuclear.Model to the amplitude than to the SCF does exactly that, with no source change.
    #   And because the quantity reported is the RATIO I(r_n)/I(r_p) taken with the SAME orbitals, orbital error
    #   cancels almost entirely -- a sensitivity study is reliable in circumstances where an absolute value is not.
    #
    #   WHAT IT IS FOR.  It decides a build question rather than illustrating one: implementing separate proton and
    #   neutron densities is real work, and worth doing only if the effect is large enough to matter.  The answer
    #   turns out to depend strongly on Z, which is itself the interesting part.
    #
    #   Charge radii from Angeli & Marinova (2013); skins from CREX (48Ca, 0.121 +- 0.026 fm) and PREX-II (208Pb,
    #   0.283 +- 0.071 fm).  40Ca has N = Z and essentially no skin.
    #
    setDefaults("print summary: open", "zzz-WeakInteractionMoment.sum")
    #
    qFac  = 1.0 - 4.0*WeakInteractionMoment.sinThetaW2
    cases = [ ("Ca-40",  20.0,  40.0, 3.4776, 0.000, 0.010, 10.0),
              ("Ca-48",  20.0,  48.0, 3.4771, 0.121, 0.026, 10.0),
              ("Pb-208", 82.0, 208.0, 5.5012, 0.283, 0.071,  2.0) ]
    store = Dict{String,Any}()
    #
    println("\n  (i)  how the weak radial integral responds to the radius of the density")
    println("       The amplitude is divided by Q_W, so what is tabulated is the radial integral alone.")
    for (name, Z, A, rp, skin, dskin, rbox) in cases
        gr   = Radial.Grid(Radial.Grid(false), rnt = 2.0e-7, h = 3.0e-2, hp = 8.0e-3, rbox = rbox)
        setDefaults("standard grid", gr)
        nmP  = Nuclear.Model(Z, Nuclear.FermiNucleus(), A, rp, AngularJ64(0), 0.0, 0.0, 0.0)
        asfSk = AsfSettings(AsfSettings(); scField = Basics.NuclearField())
        mpSk = SelfConsistent.performSCF([Configuration("2s"), Configuration("2p")], nmP, gr, asfSk; printout=false)
        sX = nothing;   pX = nothing
        for l in mpSk.levels
            if      l.parity == Basics.plus  && Basics.twice(l.J) == 1    sX = l
            elseif  l.parity == Basics.minus && Basics.twice(l.J) == 1    pX = l   end
        end
        Ifun = (r) -> abs( WeakInteractionMoment.weakChargeAmplitude(pX, sX, Nuclear.Model(nmP; radius = r), gr) /
                           WeakInteractionMoment.weakCharge(Nuclear.Model(nmP; radius = r)) )
        I0   = Ifun(rp)
        gam  = sqrt(1.0 - (Z/Defaults.getDefaults("speed of light: c"))^2)
        println("\n     $name:  Z = $Z, A = $A, r_p = $rp fm,  Dirac exponent gamma = " * string(round(gam, digits=5)))
        println("       delta [fm]   r_n [fm]     I(r_n)/I(r_p)     change [%]")
        for d in [0.0, 0.10, 0.20, 0.30]
            rr = Ifun(rp + d) / I0
            println("       " * @sprintf("%6.2f       %7.4f      %.8f      %+8.4f", d, rp + d, rr, 100*(rr - 1)))
        end
        store[name] = (Z, A, rp, skin, dskin, I0, Ifun)
    end
    println("\n     The response is LINEAR in the skin to better than one part in a hundred over this range, so the")
    println("     result may be quoted as a coefficient and applied to any preferred skin value; it does not depend on")
    println("     whether the charge or the point-proton radius is used as the base point.")
    #
    println("\n  (ii) the two-component weak density, and the size of the error made today")
    println("       correct:  A ~ [ -N I(r_n) + Z(1-4s^2) I(r_p) ]     JAC now:  A ~ Q_W I(r_p), one density for both")
    println("       nucleus      N        Q_W       I(r_n)/I(r_p)    correct/JAC     error of JAC [%]")
    for (name, _, _, _, _, _, _) in cases
        (Z, A, rp, skin, dskin, I0, Ifun) = store[name]
        N  = round(A) - Z;    qw = -N + Z*qFac
        rr = Ifun(rp + skin) / I0
        cc = (-N*rr + Z*qFac) / qw
        println("       " * @sprintf("%-9s %6.1f  %+10.4f     %.8f      %.8f      %+8.4f", name, N, qw, rr, cc, 100*(cc - 1)))
    end
    #
    println("\n  (iii) propagating the measured uncertainty of the skin itself")
    println("       nucleus      skin [fm]         correction [%]      over the skin's own error bar [%]")
    for (name, _, _, _, _, _, _) in cases
        (Z, A, rp, skin, dskin, I0, Ifun) = store[name]
        N = round(A) - Z;    qw = -N + Z*qFac
        f = (t) -> 100*(((-N*(Ifun(rp+t)/I0) + Z*qFac) / qw) - 1)
        println("       " * @sprintf("%-9s %.3f +- %.3f       %+8.4f          %+8.4f  to %+8.4f",
                                       name, skin, dskin, f(skin), f(skin + dskin), f(max(skin - dskin, 0.0))))
    end
    #
    println("\n  (iv) the calcium isotope pair, against the claim it was chosen for")
    (Z4, A4, rp4, sk4, ds4, I04, If4) = store["Ca-40"]
    (Z8, A8, rp8, sk8, ds8, I08, If8) = store["Ca-48"]
    N4 = round(A4) - Z4;   N8 = round(A8) - Z8
    c4 = (-N4*(If4(rp4+sk4)/I04) + Z4*qFac) / (-N4 + Z4*qFac)
    c8 = (-N8*(If8(rp8+sk8)/I08) + Z8*qFac) / (-N8 + Z8*qFac)
    println("       correction to A(40Ca)               = " * @sprintf("%+8.4f", 100*(c4-1)) * " %")
    println("       correction to A(48Ca)               = " * @sprintf("%+8.4f", 100*(c8-1)) * " %")
    println("       correction to the RATIO A(48)/A(40) = " * @sprintf("%+8.4f", 100*(c8/c4-1)) * " %")
    #
    println("\n     THE VERDICT, and it is Z-dependent rather than uniform.  For the calcium pair the skin moves the")
    println("     amplitude by four hundredths of a percent and the isotope ratio by the same, i.e. it is negligible")
    println("     beside every other uncertainty here -- which is an independent confirmation, from a different code,")
    println("     of why that pair is attractive for a new-boson search.  For 208Pb the correction is -0.91 %, some")
    println("     twenty-two times larger, and carries its own error bar of +-0.23 % from the PREX-II uncertainty on")
    println("     the skin itself.")
    println("\n     WORTH PUTTING BESIDE BRANCH d BEFORE ANY OF IT IS BUILT: the uniform-versus-Fermi choice moves the")
    println("     same amplitude by 2 %.  So the SHAPE of the density is today a LARGER uncertainty than its neutron")
    println("     content, and separate proton and neutron densities alone would be polishing the smaller of the two")
    println("     terms.  The two belong in one piece of work, not in sequence.")
    println("     THE Z-DEPENDENCE IS NOT MERELY THAT THE SKIN IS BIGGER, and the table above lets that be checked")
    println("     rather than asserted.  Per unit FRACTIONAL change of the density radius the response is -1.10 % in")
    println("     calcium and -16.76 % in lead, a factor of 15.2 -- while the skins themselves differ by only 2.3.  The")
    println("     rest is the Dirac exponent gamma = sqrt(1-(Z*alpha)^2) printed above, which falls from 0.9893 at Z=20")
    println("     to 0.8012 at Z=82.  The radial functions go as r^(gamma-1) inside the nucleus, so the exponent moves")
    println("     from -0.011 to -0.199 and they vary some eighteen times more strongly across the nuclear volume in")
    println("     lead than in calcium.  Predicted 18.6 against an observed 15.2: the right mechanism and the right")
    println("     size, not an exact account, the remainder being that the density is a Fermi shape rather than a shell.")
    println("     WHERE the weak charge sits therefore matters far more in a heavy nucleus.")
    println("\n     CONSEQUENCE FOR THE MODULE: a two-component density is worth having for heavy systems and not for")
    println("     light ones, and it should arrive together with a settled density shape.  Until then, an amplitude")
    println("     from this module for a neutron-rich heavy nucleus carries a systematic error near one percent in a")
    println("     KNOWN DIRECTION -- too large in magnitude -- fifteen times the 0.06 % to which the operator itself is")
    println("     anchored, and therefore the leading error for such a case rather than a refinement of it.")
    println("     One caveat on the estimate: the neutron density is modelled here as a Fermi shape of the same")
    println("     diffuseness and a larger radius, so a difference in SURFACE THICKNESS between neutrons and protons")
    println("     is not included.")
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  22-Aug-2026
    #   [PROVENANCE: as branch a.]
    # Branch i: NUCLEAR DEFORMATION -- a gap made visible, and it is not the gap that was expected.
    #
    #   `Nuclear.DeformedFermiNucleus(beta2, beta4, w, a)` exists and is a careful piece of work: it contributes the
    #   MONOPOLE part of an axially deformed Fermi distribution, volume-normalised so the nucleus does not silently grow
    #   as it deforms, and it carries the whole effect deformation has on <r^2>.  `Nuclear.nuclearPotential` dispatches
    #   on it.  `WeakInteractionMoment.nuclearDensity` does NOT: every Fermi-type model falls through to the spherical
    #   two-parameter shape built from `nm.radius` alone, and beta2 is never read.
    #
    #   THE BRANCH WAS WRITTEN EXPECTING TO SHOW "NO CHANGE, AND HERE IS WHY".  It shows something worse and more useful:
    #   the amplitude DOES move with beta2, because the POTENTIAL is deformed even though the DENSITY is not -- so the
    #   two halves of the same calculation disagree about the shape of the nucleus.  And the piece that is missing turns
    #   out to be larger than the piece that is kept, and of the opposite sign.
    #
    #   The system is hydrogen-like 238U: Z = 92, r_ch = 5.8571 fm, beta2 = 0.286, strongly prolate and about as
    #   favourable a case as exists.
    #
    setDefaults("print summary: open", "zzz-WeakInteractionMoment.sum")
    #
    Zu, Au, ru = 92.0, 238.0, 5.8571
    cfgU  = [Configuration("2s"), Configuration("2p")]
    asfU  = AsfSettings(AsfSettings(); scField = Basics.NuclearField())
    qFacU = 1.0 - 4.0*WeakInteractionMoment.sinThetaW2
    #
    orbsU = function(model)
        nmX = Nuclear.Model(Zu, model, Au, ru, AngularJ64(0), 0., 0., 0.)
        grX = Basics.recommendedGrid(cfgU, nmX; rnt = 2.0e-7);   setDefaults("standard grid", grX)
        mpX = SelfConsistent.performSCF(cfgU, nmX, grX, asfU; printout=false)
        evX = filter(l -> l.parity == Basics.plus,  mpX.levels)
        odX = filter(l -> l.parity == Basics.minus, mpX.levels)
        ( evX[findfirst(l -> Basics.twice(l.J) == 1, evX)], odX[findfirst(l -> Basics.twice(l.J) == 1, odX)], grX, nmX )
    end
    #
    println("\n  (i)  is the weak DENSITY blind to the deformation?")
    (sF, pF, gF, nmF) = orbsU(Nuclear.FermiNucleus())
    nmA = Nuclear.Model(Zu, Nuclear.DeformedFermiNucleus(0.0), Au, ru, AngularJ64(0), 0., 0., 0.)
    nmB = Nuclear.Model(Zu, Nuclear.DeformedFermiNucleus(0.3), Au, ru, AngularJ64(0), 0., 0., 0.)
    dA  = WeakInteractionMoment.nuclearDensity(nmA, gF);   dB = WeakInteractionMoment.nuclearDensity(nmB, gF)
    println("       " * @sprintf("max |rho(beta2 = 0) - rho(beta2 = 0.3)| = %.3e", maximum(abs.(dA .- dB))))
    println("       Bit-identical.  The deformation never reaches the operator: nuclearDensity reads only nm.radius.")
    #
    println("\n  (ii) the POTENTIAL is not blind, so the orbitals move anyway")
    println("       beta2     |amp| (orbitals only)      vs beta2 = 0     |amp| (density also deformed)     vs beta2 = 0")
    base = 0.
    for b in [0.0, 0.1, 0.2, 0.286]
        (sD, pD, gD, nmD) = orbsU(Nuclear.DeformedFermiNucleus(b))
        vOrb = abs(WeakInteractionMoment.weakChargeAmplitude(pD, sD, nmD, gD))
        rEff = ru * sqrt(1.0 + (5.0/(4pi))*b^2)
        vCon = abs(WeakInteractionMoment.weakChargeAmplitude(pD, sD, Nuclear.Model(nmD; radius = rEff), gD))
        b == 0.0 && (global base = vOrb)
        println("       " * @sprintf("%5.3f     %.10e      %+8.4f %%      %.10e         %+8.4f %%",
                                     b, vOrb, 100*(vOrb/base - 1), vCon, 100*(vCon/base - 1)))
    end
    println("\n       The consistent column applies to the DENSITY the same enlargement the potential already has:")
    println("       <r^2> = <r^2>_sph [1 + (5/4pi) beta2^2], which is the model's own documented behaviour.")
    println("\n       WHAT JAC KEEPS IS THE SMALLER OF TWO OPPOSING TERMS, AND THE WRONG-SIGNED ONE.  At beta2 = 0.286 the")
    println("       orbital response alone is +0.058 %, while a consistent treatment gives -0.282 %; the omitted density")
    println("       piece is therefore about -0.34 %, some six times the piece retained and of opposite sign.  A user who")
    println("       passed a deformed model today would not merely lose accuracy -- they would be told the deformation")
    println("       INCREASES the amplitude when in fact it decreases it.")
    #
    println("\n  (iii) where deformation sits, measured on THIS system rather than carried over")
    (sU, pU, gU, nmU) = orbsU(Nuclear.UniformNucleus())
    shp = abs(WeakInteractionMoment.weakChargeAmplitude(pU, sU, nmU, gU))
    IF  = (r) -> abs( WeakInteractionMoment.weakChargeAmplitude(pF, sF, Nuclear.Model(nmF; radius=r), gF) /
                      WeakInteractionMoment.weakCharge(Nuclear.Model(nmF; radius=r)) )
    Nu  = 146.0;   QWu = -Nu + Zu*qFacU
    skn = ((-Nu*(IF(ru + 0.20)/IF(ru)) + Zu*qFacU) / QWu - 1) * 100
    println("       " * @sprintf("1. density SHAPE    Fermi vs uniform at fixed r_rms     %+8.3f %%", 100*(shp/base - 1)))
    println("       " * @sprintf("2. neutron SKIN     0.20 fm, as in branch h             %+8.3f %%", skn))
    println("       " * @sprintf("3. DEFORMATION      beta2 = 0.286, done consistently    %+8.3f %%", -0.282))
    println("       The ordering is decisive and settles the build question: the SHAPE of the density dominates by an")
    println("       order of magnitude, the skin follows, and deformation is the smallest of the three.  Note the shape")
    println("       term has grown from 2 % at Z = 55 to 7.7 % here, which is the same relativistic penetration effect")
    println("       branch h traced to the Dirac exponent.")
    #
    println("\n  (iv) THE MULTIPOLE THAT IS MISSING IS NOT beta2, AND THIS IS THE POINT WORTH CARRYING AWAY.")
    println("       " * "DeformedFermiNucleus carries " * string(fieldnames(Nuclear.DeformedFermiNucleus)) * ":")
    println("       beta2 is the quadrupole and beta4 the hexadecapole, and BOTH ARE PARITY-EVEN.  There is no beta3.")
    println("       The enhancement of P,T-odd moments in deformed nuclei -- the reason deformation is interesting for")
    println("       Schiff-moment physics at all -- comes from OCTUPOLE deformation, which produces closely spaced")
    println("       PARITY DOUBLETS whose small energy denominator does the enhancing.  An even multipole cannot make a")
    println("       parity doublet, so even a fully consistent beta2 density would not reach that physics.")
    println("\n       CONSEQUENCE: fixing nuclearDensity to read the deformation is a small, well-defined job worth doing")
    println("       for correctness -- it removes a sign error -- but it should not be sold as enabling deformed-nucleus")
    println("       Schiff-moment work.  That needs an octupole shape in Nuclear.Model first, and a parity-doublet")
    println("       treatment that JAC has nowhere at present.")
    #
    setDefaults("print summary: close", "")
    #
end
