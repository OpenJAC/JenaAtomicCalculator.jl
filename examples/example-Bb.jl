
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
    #   THIS BRANCH EARNED ITS KEEP ON THE DAY IT WAS WRITTEN.  It originally printed a third column for the anapole,
    #   which came out identically 0.0 for all sixteen pairs -- a pattern that reads like a selection rule and is not
    #   one.  The cause was a wrong angular structure (the magnetic-multipole template, which is parity-EVEN and so
    #   forbids exactly the s_1/2 <-> p_1/2 mixing through which the anapole reaches a measured PNC transition).
    #   `anapoleAmplitude` now RAISES rather than returning that silent zero, and the column is gone from this table.
    #   See branch d, where the raise is exercised.
    #
    #   The test bed is deliberately small and has both parities and two J values: a one-electron ion at Z = 55 in the
    #   BARE nuclear field, giving 1s_1/2+, 2s_1/2+, 2p_1/2- and 2p_3/2-.  Basics.NuclearField() is used rather than the
    #   default DFS so that a one-electron system does not acquire a self-interaction.
    #
    #   RESULT, 21-Aug-2026:  16 ordered pairs, of which 8 have equal parity.  All 8 give exactly 0.0 for both operators.
    #   Of the 8 opposite-parity pairs, the weak charge is non-zero for exactly the 4 with J_f = J_i and exactly zero for
    #   the 4 with J_f != J_i, while the rank-1 Schiff amplitude survives the J change.  Every zero here is a hard 0.0,
    #   not a small number.
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
    println("\n\n  Selection rules of the P-odd operators (Z = 55, bare nuclear field)\n")
    println("   f <- i     J^P f     J^P i    same par?   dJ?      weak charge        Schiff")
    nSame = 0;   nViol = 0;   nWeak = 0
    for  a in mp.levels,  b in mp.levels
        wc = WeakInteractionMoment.weakChargeAmplitude(a, b, nModel, grid)
        sm = WeakInteractionMoment.schiffMomentAmplitude(a, b, nModel, grid)
        samePar = (a.parity == b.parity);    sameJ = (a.J == b.J)
        if  samePar
            global nSame += 1
            if  wc != 0. || sm != 0.    global nViol += 1    end
        end
        if  wc != 0.    global nWeak += 1    end
        println("   " * rpad("$(a.index) <- $(b.index)", 10) * rpad(string(LevelSymmetry(a.J,a.parity)), 10) *
                rpad(string(LevelSymmetry(b.J,b.parity)), 9) * rpad(samePar ? "yes" : "no", 12) *
                rpad(sameJ ? "no" : "yes", 9) * rpad(string(round(abs(wc), sigdigits=4)), 18) *
                string(round(abs(sm), sigdigits=4)))
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
    println("\n  the anapole amplitude, which is not implemented and must say so rather than return a number:")
    try
        nmA = Nuclear.Model(55.)
        mA  = SelfConsistent.performSCF([Configuration("1s"), Configuration("2p")], nmA, grid, asfB; printout=false)
        WeakInteractionMoment.anapoleAmplitude(mA.levels[1], mA.levels[2], nmA, grid)
        println("     NO ERROR RAISED -- that is a defect")
    catch e
        println("     raised, as intended:  ", first(sprint(showerror, e), 150), " ...")
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
    # Last visit:  21-Aug-2026
    # Last successful:  unknown ...
    # Branch e: THE Z-DEPENDENCE -- DELIBERATELY LEFT OPEN, and the reason is worth reading before re-running it.
    #
    #   This branch was written to check the amplitudes against the Bouchiat "Z^3 law", the standard statement that a
    #   parity-non-conserving amplitude grows as the cube of the nuclear charge.  IT DOES NOT APPLY HERE, and the branch
    #   is kept undated rather than quietly re-targeted at whatever exponent came out.
    #
    #   Why it does not apply: the Z^3 law describes a valence electron of a NEUTRAL heavy atom, where the Z-dependence
    #   of the wave function near the nucleus combines with Q_W ~ -N ~ -Z in a particular way.  What is computed here is
    #   a hydrogen-like ion, whose orbitals scale differently, and whose nuclear radius is itself Z-dependent
    #   (R ~ A^(1/3)), so the density normalization 1/R^3 carries its own Z-dependence into the integral.  There is in
    #   addition a relativistic factor that grows sharply once (Z alpha) is no longer small.
    #
    #   MEASURED, 21-Aug-2026, after dividing out the trivial Q_W(Z): the LOCAL exponent d ln|amp| / d ln Z rises
    #   monotonically from about 4.4 at Z = 20-30 to about 7.3 at Z = 70-80, with a least-squares value of 5.28 over the
    #   whole range.  A rising local exponent is what a relativistic enhancement looks like, so the numbers are not
    #   obviously wrong -- but "not obviously wrong" is not a verification, and NOTHING here is claimed as one.
    #
    #   TO CLOSE THIS BRANCH someone needs the correct hydrogen-like expectation, i.e. the analytic Z-dependence of
    #   INT rho [P_(np) Q_(ns) - Q_(np) P_(ns)] dr for finite-nucleus Dirac orbitals, including the (2 Z R_nuc)^(2gamma-2)
    #   factor with gamma = sqrt(1 - (Z alpha)^2).  With that in hand this becomes a real test; until then it is a scan.
    setDefaults("print summary: open", "zzz-WeakInteractionMoment.sum")
    #
    asfB = AsfSettings(AsfSettings(); scField = Basics.NuclearField())
    Zs   = [20.0, 30.0, 40.0, 55.0, 70.0, 80.0]
    vals = Float64[]
    println("\n\n  Z-dependence of |<2p_1/2||H_W||1s_1/2>| / |Q_W|   -- A SCAN, NOT A TEST\n")
    println("      Z      |amp|/|Q_W|            local exponent d(ln)/d(lnZ)")
    for  Zi in Zs
        gr = Radial.Grid(Radial.Grid(false), rnt = 4.0e-7, h = 3.0e-2, hp = 1.0e-2, rbox = 6.0)
        setDefaults("standard grid", gr)
        nm = Nuclear.Model(Zi)
        m  = SelfConsistent.performSCF([Configuration("1s"), Configuration("2p")], nm, gr, asfB; printout=false)
        ss = m.levels[argmin([l.energy for l in m.levels])]
        pp = nothing
        for  l in m.levels    if  l.parity == Basics.minus && Basics.twice(l.J) == 1    pp = l    end    end
        v  = abs(WeakInteractionMoment.weakChargeAmplitude(pp, ss, nm, gr)) / abs(WeakInteractionMoment.weakCharge(nm))
        push!(vals, v)
        sa = length(vals) > 1 ?
             string(round((log(vals[end])-log(vals[end-1])) / (log(Zi)-log(Zs[length(vals)-1])), digits=4)) : "--"
        println("   " * rpad(string(Zi), 8) * rpad(string(round(v, sigdigits=7)), 22) * sa)
    end
    n = length(Zs);   lx = log.(Zs);   ly = log.(vals)
    slope = (n*sum(lx.*ly) - sum(lx)*sum(ly)) / (n*sum(lx.^2) - sum(lx)^2)
    println("\n   least-squares exponent over the whole range = $(round(slope, digits=4))")
    println("   the Bouchiat Z^3 law is NOT the expectation here; see the comment above this branch.")
    #
    setDefaults("print summary: close", "")
    #
end
