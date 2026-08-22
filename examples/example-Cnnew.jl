#
println("Cn) Apply & test the WeakInteractionEnhancement module: the OBSERVABLES that the P-odd interactions of")
println("    WeakInteractionMoment give rise to -- the parity-non-conserving E1 amplitude between two levels of")
println("    the same parity, and the enhancement factor R = d_atom/d_e of an electron electric-dipole moment.")
println("    Both are sums over intermediate states, supplied by the caller as a gMultiplet; no Green function")
println("    and no pseudo-continuum is used.  The branches below are arranged so that the two observables are")
println("    NOT presented alike: E1_PNC converges usefully and is shown converging, whereas R is a cancelling")
println("    sum that a truncated intermediate set cannot approximate, and branch d demonstrates that rather")
println("    than hiding it.")

using Printf

if  true
    # Last successful:  21-Aug-2026
    #   [PROVENANCE: produced on Julia 1.10.9, with a Manifest re-resolved for 1.10 because the checkout's Manifest had been
    #    resolved under 1.12.6 on the maintainer machine.  Re-run there to confirm.]
    # Branch a: THE NORMALIZATION, MEASURED; AND THE STRUCTURE OF BOTH OPERATORS.  Everything this module computes is a
    #   PRODUCT of a dipole matrix element and a P-odd one, so a normalization error in either factor is an error in the
    #   observable, and the two factors must share one convention.  That is why the dipole here is routed through
    #   WeakInteractionMoment.oneParticleAmplitude rather than through MultipoleMoment: the same contraction, hence the
    #   same convention by construction.  This branch measures what the two routes differ by, rather than assuming it.
    #
    #   Three further facts, none of which a mis-assembled operator would satisfy:
    #     (i)   E1_PNC is PURELY IMAGINARY, inherited from gamma_5 in the weak-charge amplitude;
    #     (ii)  it vanishes exactly between levels of OPPOSITE parity -- it is the same-parity amplitude, by definition;
    #     (iii) both effective operators are HERMITIAN, but with OPPOSITE exchange symmetry, and that difference is the
    #           point: <f||T||i> = (-1)^(J_f-J_i) conj(<i||T||f>) forces the purely imaginary weak-charge amplitude to be
    #           ANTIsymmetric under exchange and the purely real EDM amplitude to be SYMMETRIC.  Getting one of the two
    #           backwards is exactly the error that cost a factor of -3 in the anapole (see example-Bb.jl, branch f).
    #
    #   ONE THING WORTH KNOWING, found here: MultipoleMoment.emmStaticAmplitude cannot span two separately computed
    #   multiplets.  It indexes finalLevel.basis.orbitals with coefficients built from the INITIAL basis and raises a
    #   KeyError when the two bases differ.  oneParticleAmplitude merges the subshell lists first and does not.  Hence the
    #   single mixed-parity multiplet below.
    #
    setDefaults("print summary: open", "zzz-WeakInteractionEnhancement.sum")
    #
    Z    = 20.0
    grid = Radial.Grid(Radial.Grid(false), rnt = 2.0e-7, h = 3.0e-2, hp = 8.0e-3, rbox = 10.0)
    setDefaults("standard grid", grid)
    nm   = Nuclear.Model(Z)
    # a one-electron system, so that a disagreement is the operator's fault and not the atomic model's; the DFS default
    # would add a self-interaction that a single electron must not have (see example-Bb.jl)
    asfB = AsfSettings(AsfSettings(); scField = Basics.NuclearField())
    mp   = SelfConsistent.performSCF([Configuration("2s"), Configuration("3s"), Configuration("2p"), Configuration("3p"),
                                      Configuration("4p")], nm, grid, asfB; printout=false)
    evens = filter(l -> l.parity == Basics.plus,  mp.levels)
    odds  = filter(l -> l.parity == Basics.minus, mp.levels)
    #
    println("\n  (i)  the reduced dipole, this module against MultipoleMoment.emmStaticAmplitude")
    println("        f              i          this module      emmStatic        ratio       sqrt(2J_f+1)")
    for f in odds, i in evens
        d1 = WeakInteractionEnhancement.dipoleReducedMe(f, i, grid)
        d2 = MultipoleMoment.emmStaticAmplitude(1, f, i, grid)
        (abs(d1) == 0. || abs(d2) == 0.) && continue
        println("     " * rpad("$(f.index) " * string(LevelSymmetry(f.J,f.parity)),12) *
                " " * rpad("$(i.index) " * string(LevelSymmetry(i.J,i.parity)),12) *
                "  " * @sprintf("%+.6e", real(d1)) * "  " * @sprintf("%+.6e", real(d2)) *
                "  " * @sprintf("%+.8f", real(d1)/real(d2)) * "   " * @sprintf("%.6f", sqrt(Basics.twice(f.J)+1.0)))
    end
    println("       The ratio must be sqrt(2J_f+1) EXACTLY -- a clean angular factor, not a number near one.")
    #
    println("\n  (ii) E1_PNC: reality, and the opposite-parity zero")
    for f in evens, i in evens
        a = WeakInteractionEnhancement.computePncE1Amplitude(f, i, mp, nm, grid)
        abs(a) == 0. && continue
        println("     " * rpad("$(f.index) <- $(i.index)", 10) * "re = " * @sprintf("%+.4e", real(a)) *
                "   im = " * @sprintf("%+.6e", imag(a)) * "   re/|amp| = " * @sprintf("%.2e", abs(real(a))/abs(a)))
    end
    nViol = 0
    for f in evens, i in odds
        WeakInteractionEnhancement.computePncE1Amplitude(f, i, mp, nm, grid) != 0.  &&  (global nViol += 1)
    end
    println("     opposite-parity pairs with a non-zero E1_PNC: $nViol   (must be 0)")
    #
    println("\n  (iii) exchange symmetry: ANTIsymmetric for the weak charge, SYMMETRIC for the EDM operator")
    for f in evens, i in evens
        f.index >= i.index && continue
        x = WeakInteractionEnhancement.computePncE1Amplitude(f, i, mp, nm, grid)
        y = WeakInteractionEnhancement.computePncE1Amplitude(i, f, mp, nm, grid)
        abs(x) == 0. && continue
        println("     PNC   $(f.index)<->$(i.index)   |<f|i> + <i|f>| = " * @sprintf("%.3e", abs(x+y)) *
                "   |<f|i> - <i|f>| = " * @sprintf("%.3e", abs(x-y)) * "     (the SUM must vanish)")
    end
    for f in odds, i in evens
        x = WeakInteractionEnhancement.edmAmplitude(f, i, nm, grid)
        abs(x) == 0. && continue
        y = WeakInteractionEnhancement.edmAmplitude(i, f, nm, grid)
        println("     EDM   $(f.index)<->$(i.index)   |<f|i> + <i|f>| = " * @sprintf("%.3e", abs(x+y)) *
                "   |<f|i> - <i|f>| = " * @sprintf("%.3e", abs(x-y)) * "     (the DIFFERENCE must vanish)")
    end
    println("\n     For a one-electron system this last pair of facts follows already from the radial integrals, which are")
    println("     antisymmetric and symmetric respectively; it becomes a real test of the angular algebra only for a")
    println("     many-electron case, and is reported here for what it is.")
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  21-Aug-2026
    #   [PROVENANCE: as branch a.]
    # Branch b: DOES THE SUM CONVERGE?  This is the branch that decides whether a number may be quoted at all, and it is
    #   the reason the module takes an explicit gMultiplet instead of hiding one.  E1_PNC is computed for the same pair
    #   with a growing set of intermediate np levels.
    #
    #   The system is LI-LIKE calcium rather than hydrogen-like, deliberately: in a one-electron ion the 2s_1/2 and
    #   2p_1/2 levels are degenerate in Dirac theory and are separated only by QED, which this calculation does not have
    #   (that is branch c).  Screening in a three-electron ion splits them properly, so what is measured here is the
    #   convergence of the sum and not the accuracy of one denominator.
    #
    #   WHAT TO LOOK FOR is not the final value but the FIRST ROW: a single intermediate gives the WRONG SIGN.  The
    #   dominant 2p term is overturned by 3p, after which the tail settles at the per-cent level.  Any calculation that
    #   stops at the leading intermediate -- which is what "the sum is dominated by the lowest np level" is often taken
    #   to license -- gets not merely a poor number but the opposite one.
    #
    setDefaults("print summary: open", "zzz-WeakInteractionEnhancement.sum")
    #
    Z    = 20.0
    grid = Radial.Grid(Radial.Grid(false), rnt = 2.0e-7, h = 3.0e-2, hp = 8.0e-3, rbox = 12.0)
    setDefaults("standard grid", grid)
    nm   = Nuclear.Model(Z)
    base = [Configuration("1s^2 2s"), Configuration("1s^2 3s")]
    sets = [ ("2p",             [Configuration("1s^2 2p")]),
             ("2p 3p",          [Configuration("1s^2 2p"), Configuration("1s^2 3p")]),
             ("2p 3p 4p",       [Configuration("1s^2 2p"), Configuration("1s^2 3p"), Configuration("1s^2 4p")]),
             ("2p 3p 4p 5p",    [Configuration("1s^2 2p"), Configuration("1s^2 3p"), Configuration("1s^2 4p"),
                                 Configuration("1s^2 5p")]) ]
    #
    println("\n  E1_PNC for the Li-like Ca  1s^2 3s <- 1s^2 2s  pair, as the intermediate set grows")
    println("     intermediate set      n(odd levels)     E1_PNC [i e a_0]        change")
    prev = 0.
    for (name, odd) in sets
        mpl = SelfConsistent.performSCF(vcat(base, odd), nm, grid, AsfSettings(); printout=false)
        ev  = filter(l -> l.parity == Basics.plus, mpl.levels)
        od  = filter(l -> l.parity == Basics.minus, mpl.levels)
        s2  = ev[argmin([l.energy for l in ev])];    s3 = ev[argmax([l.energy for l in ev])]
        a   = WeakInteractionEnhancement.computePncE1Amplitude(s3, s2, mpl, nm, grid)
        sa  = prev == 0. ? "     --" : @sprintf("%+8.2f %%", 100*(imag(a)/prev - 1))
        println("     " * rpad(name,20) * "   " * lpad(string(length(od)),5) * "        " *
                @sprintf("%+.8e", imag(a)) * "     " * sa)
        global prev = imag(a)
    end
    println("\n     The first row has the opposite sign to every later one.  The last two rows agree to about 1 %, which is")
    println("     what may be quoted as the truncation uncertainty of this intermediate set -- not a claim of accuracy.")
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  21-Aug-2026
    #   [PROVENANCE: as branch a.]
    # Branch c: A LIMITATION, DEMONSTRATED RATHER THAN DESCRIBED.  For a HYDROGEN-LIKE ion the 2s_1/2 and 2p_1/2 levels
    #   are exactly degenerate in Dirac theory.  What splits them in nature is the Lamb shift, of which this calculation
    #   contains only the finite-nuclear-size part; the self-energy and vacuum-polarization pieces, which dominate at
    #   moderate Z, are simply absent.  The PNC sum has that splitting in a DENOMINATOR, so the error goes straight into
    #   the answer and is not small.
    #
    #   The branch prints, for the same ion, the splitting this calculation produces and the leading-order estimate of
    #   the true one, taken from the measured hydrogen 2s-2p_1/2 Lamb shift of 1057.845 MHz scaled as Z^4.  It then shows
    #   how much of E1_PNC rides on that one denominator.
    #
    #   THE CONCLUSION IS NEGATIVE AND IS THE USEFUL PART: E1_PNC from this module is NOT quantitative for hydrogen-like
    #   ions, and no amount of enlarging the gMultiplet repairs it, because the defect is in an energy the multiplet does
    #   not contain.  A many-electron system, where screening does the splitting, is the regime in which the module can
    #   be believed.  Getting a plausible-looking number out of the H-like case would have been easy and wrong.
    #
    setDefaults("print summary: open", "zzz-WeakInteractionEnhancement.sum")
    #
    Z    = 20.0
    grid = Radial.Grid(Radial.Grid(false), rnt = 2.0e-7, h = 3.0e-2, hp = 8.0e-3, rbox = 10.0)
    setDefaults("standard grid", grid)
    nm   = Nuclear.Model(Z)
    asfB = AsfSettings(AsfSettings(); scField = Basics.NuclearField())
    mp   = SelfConsistent.performSCF([Configuration("2s"), Configuration("3s"), Configuration("2p"), Configuration("3p"),
                                      Configuration("4p")], nm, grid, asfB; printout=false)
    evens = filter(l -> l.parity == Basics.plus,  mp.levels)
    odds  = filter(l -> l.parity == Basics.minus, mp.levels)
    s2    = evens[argmin([l.energy for l in evens])];    s3 = evens[argmax([l.energy for l in evens])]
    p2    = odds[argmin([l.energy for l in odds])]
    #
    dEcalc = s2.energy - p2.energy
    dELamb = 1.6083e-7 * Z^4                    # 1057.845 MHz in a.u., scaled as Z^4; leading order only
    println("\n  (i)  the 2s_1/2 -- 2p_1/2 splitting of this hydrogen-like ion at Z = $Z")
    println("       from this calculation (finite nuclear size only)  = " * @sprintf("%.6e", abs(dEcalc)) * " a.u.")
    println("       leading-order QED estimate, 1.6083e-7 * Z^4       = " * @sprintf("%.6e", dELamb) * " a.u.")
    println("       ratio, true / calculated                          = " * @sprintf("%.1f", dELamb/abs(dEcalc)))
    #
    println("\n  (ii) how much of E1_PNC rides on that denominator")
    tot = WeakInteractionEnhancement.computePncE1Amplitude(s3, s2, mp, nm, grid)
    for n in odds
        a = WeakInteractionEnhancement.computePncE1Amplitude(s3, s2, Multiplet("one", [n]), nm, grid)
        imag(a) == 0. && continue
        println("       n = " * rpad("$(n.index) " * string(LevelSymmetry(n.J,n.parity)),11) *
                " dE(i) = " * @sprintf("%+.3e", s2.energy - n.energy) *
                "   contribution = " * @sprintf("%+.6e", imag(a)) * " i    share = " * @sprintf("%7.2f", 100*imag(a)/imag(tot)) * " %")
    end
    println("       TOTAL = " * @sprintf("%+.6e", imag(tot)) * " i")
    println("\n       The 2p_1/2 share is carried by a denominator too small by the factor printed above, so it is the one")
    println("       term of the sum that this calculation cannot get right.  The levels of j = 3/2 contribute exactly")
    println("       zero, as they must: the weak charge is a rank-0 operator and reaches only J_n = J_i.")
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  21-Aug-2026
    #   [PROVENANCE: as branch a.]
    # Branch d: WHY R IS PRINTED AS UNCONVERGED, AND WHAT THAT DOES AND DOES NOT MEAN.  Schiff's theorem says that the
    #   electric-dipole moment of a POINT, NON-RELATIVISTIC bound system vanishes identically: the sum over states
    #   cancels.  The enhancement of an electron EDM is what survives that cancellation once relativity is restored,
    #   which makes d_atom a DIFFERENCE OF LARGE TERMS rather than a sum of small ones.
    #
    #   THE MEASUREMENT IS MILDER THAN THE THEOREM, and this branch is written to report that rather than the tidier
    #   claim it was set up to make.  Two parts:
    #     (i)  in the hydrogen-like ion the near-degenerate 2s_1/2 and 2p_1/2 pair return values of R equal and opposite
    #          to six digits.  That mirroring IS the cancellation, caught in the act with only two levels.  It also
    #          inflates both of them to ~6e4, for the reason of branch c and not because of Schiff's theorem: the
    #          denominator is the artificial one that a calculation without QED produces.
    #     (ii) in the Li-like ion R is recomputed as the intermediate set grows -- and it CONVERGES, by +59 %, +9.5 %
    #          and +3.2 % over the sequence.  It is merely slower than E1_PNC, which settled to 1.0 % over the very
    #          same sets in branch b.  So a truncated R is not the arbitrary quantity the theorem alone would suggest;
    #          it is a slowly converging one, and the cancellation shows up as a loss of rate, not of meaning.
    #
    #   THE LABEL "UNCONVERGED" THEREFORE RESTS ON THE SECOND REASON, not the first.  Nothing anchors the ABSOLUTE scale
    #   of the EDM operator the way the closed-form hydrogenic matrix element anchors the weak charge (example-Bb.jl,
    #   branch g), so what converges here is an unverified quantity.  Anyone wanting a defensible R should anchor the
    #   operator first; the structural checks of branch a are all it has at present.
    #
    setDefaults("print summary: open", "zzz-WeakInteractionEnhancement.sum")
    #
    Z    = 20.0
    grid = Radial.Grid(Radial.Grid(false), rnt = 2.0e-7, h = 3.0e-2, hp = 8.0e-3, rbox = 10.0)
    setDefaults("standard grid", grid)
    nm   = Nuclear.Model(Z)
    asfB = AsfSettings(AsfSettings(); scField = Basics.NuclearField())
    mp   = SelfConsistent.performSCF([Configuration("2s"), Configuration("3s"), Configuration("2p"), Configuration("3p"),
                                      Configuration("4p")], nm, grid, asfB; printout=false)
    println("\n  (i)  R for the hydrogen-like ion, level by level")
    println("       level      J^P          E [a.u.]            R = d_atom/d_e")
    for l in mp.levels
        r = WeakInteractionEnhancement.computeEdmEnhancement(l, mp, nm, grid)
        println("       " * lpad(string(l.index),5) * "    " * rpad(string(LevelSymmetry(l.J,l.parity)),10) *
                @sprintf("%+.8f", l.energy) * "      " * @sprintf("%+.6e", r))
    end
    println("       The two members of the near-degenerate 2s_1/2 / 2p_1/2 pair mirror each other to six digits; j = 3/2")
    println("       gives exactly zero, the EDM operator being a scalar like the weak charge.  The large magnitude is the")
    println("       artificial denominator of branch c, not Schiff's theorem.")
    #
    grid2 = Radial.Grid(Radial.Grid(false), rnt = 2.0e-7, h = 3.0e-2, hp = 8.0e-3, rbox = 12.0)
    setDefaults("standard grid", grid2)
    base = [Configuration("1s^2 2s"), Configuration("1s^2 3s")]
    sets = [ ("2p",          [Configuration("1s^2 2p")]),
             ("2p 3p",       [Configuration("1s^2 2p"), Configuration("1s^2 3p")]),
             ("2p 3p 4p",    [Configuration("1s^2 2p"), Configuration("1s^2 3p"), Configuration("1s^2 4p")]),
             ("2p 3p 4p 5p", [Configuration("1s^2 2p"), Configuration("1s^2 3p"), Configuration("1s^2 4p"),
                              Configuration("1s^2 5p")]) ]
    println("\n  (ii) R for the Li-like ground level, as the intermediate set grows")
    println("       intermediate set        R (ground)              change")
    prev = 0.
    for (name, odd) in sets
        mpl = SelfConsistent.performSCF(vcat(base, odd), nm, grid2, AsfSettings(); printout=false)
        gl  = mpl.levels[argmin([l.energy for l in mpl.levels])]
        r   = WeakInteractionEnhancement.computeEdmEnhancement(gl, mpl, nm, grid2)
        sa  = prev == 0. ? "     --" : @sprintf("%+8.2f %%", 100*(r/prev - 1))
        println("       " * rpad(name,20) * "  " * @sprintf("%+.8e", r) * "     " * sa)
        global prev = r
    end
    println("\n       R converges, but more slowly than E1_PNC did over the same sets in branch b (3.2 % against 1.0 % at the")
    println("       last step).  That is the cancellation costing rate rather than meaning -- a weaker statement than the")
    println("       theorem alone suggests, and the one the numbers support.")
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  21-Aug-2026
    #   [PROVENANCE: as branch a.]
    # Branch e: THE ORDINARY USER ROUTE, through Atomic.Computation and perform().  The four branches above call the
    #   module's functions directly, which is the right way to test them but not the way anyone will use them; this one
    #   exercises the propertySettings dispatch that was added to Basics.perform, the printBefore path, and both display
    #   tables.  Its numerical content is the same pair as branch b, so the E1_PNC printed here should agree with the
    #   "2p 3p 4p" row there to the accuracy of the slightly different CI basis (branch b diagonalizes even and odd
    #   configurations together, perform() only the even ones).
    #
    setDefaults("print summary: open", "zzz-WeakInteractionEnhancement.sum")
    #
    Z    = 20.0
    grid = Radial.Grid(Radial.Grid(false), rnt = 2.0e-7, h = 3.0e-2, hp = 8.0e-3, rbox = 12.0)
    setDefaults("standard grid", grid)
    nm   = Nuclear.Model(Z)
    gMp  = SelfConsistent.performSCF([Configuration("1s^2 2p"), Configuration("1s^2 3p"), Configuration("1s^2 4p")],
                                     nm, grid, AsfSettings(); printout=false)
    wSet = WeakInteractionEnhancement.Settings(WeakInteractionEnhancement.Settings();
                                               calcPncE1 = true, calcEdmEnhancement = true, gMultiplet = gMp,
                                               printBefore = true)
    comp = Atomic.Computation(Atomic.Computation(); name = "Li-like Ca: PNC E1 and EDM enhancement", grid = grid,
                              nuclearModel = nm, configs = [Configuration("1s^2 2s"), Configuration("1s^2 3s")],
                              propertySettings = [wSet] )
    perform(comp; output=true)
    #
    setDefaults("print summary: close", "")
    #
end
