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
    #   [SUPERSEDED IN PART BY BRANCH h, 22-Aug-2026: the factor 379 below comes from scaling the hydrogen Lamb shift
    #    as Z^4, which OVERSHOOTS the model Hamiltonian by three to eight times; and this branch inspects only ONE of the
    #    sum's two denominators, so it UNDERSTATES the damage.  Branch h puts QedPetersburg in and finds the amplitude
    #    too large by a factor of 110 at Z = 20.  The conclusion below stands; the number does not.]
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
elseif  false
    # Last successful:  28-Aug-2026 -- AND THIS BRANCH CAUGHT A SILENT ZERO, fixed in the same commit.
    #   WHAT WAS WRONG: E1_PNC came out +0.00000000e+00 for EVERY intermediate set and BOTH ions, while the
    #   2s-2p_1/2 splitting printed correctly (36.1 eV in Ca, 235.1 eV in Pb) -- so the levels were fine and only
    #   the amplitude was dead. The branch then died with `KeyError: "Li-like 40Ca"` at the error budget, because
    #   the truncation dictionary is only filled when a change is non-zero.
    #   WHY: `SpinAngular.computeCoefficientsScalar` required leftCsf.parity == rightCsf.parity unconditionally.
    #   That is right for the one-body Hamiltonian, which is parity-EVEN, but the nuclear weak charge is a rank-0
    #   P-ODD operator and connects OPPOSITE parities. Every P-odd rank-0 call therefore got an empty coefficient
    #   list and returned exactly zero. Measured: WeakInteractionMoment.weakChargeAmplitude gave -0.0 - 0.0im.
    #   The gate now tests the OPERATOR's parity, and rank 0 still fixes J_f = J_i -- which is exactly the rule
    #   branch a of example-Bb.jl states in words.
    #   AFTER: E1_PNC = +1.24574758e-13 (Ca) and +6.76245739e-11 (Pb) for the 2p set, and the sign flip on one
    #   intermediate that this branch's own text predicts is back.
    #   [PROVENANCE: as branch a.]
    # Branch f: THE TWO SYSTEMS OF A REAL PAPER -- Li-like 40Ca and 208Pb, and an error budget that turns over with Z.
    #
    #   These are the ions of a talk on atomic parity violation in highly charged Ca and Pb ions, computing the
    #   PV-induced E1 amplitude for 1s^2 3s <- 1s^2 2s.  This branch computes the same quantity independently and, more
    #   usefully than a bare number, says WHICH uncertainty limits it -- which is the question that decides whether a
    #   system is worth pursuing.
    #
    #   WHY LI-LIKE AND NOT H-LIKE, which branch c explains and this branch demonstrates for these very ions.  In a
    #   one-electron ion 2s_1/2 and 2p_1/2 are degenerate in Dirac theory and are split only by a Lamb shift JAC does not
    #   have, which sits in a denominator of the sum.  Screening in a three-electron ion splits them properly, and the
    #   branch prints the splitting so the reader can see it is nowhere near zero.
    #
    #   THE BOX IS LEFT TO Basics.recommendedGrid, which is what Rule 12 asks for.  A hand formula was tried first and
    #   gave 3.3 a.u. where 4p and 5p need about 9.8; Bsplines.checkGridRepresentation REFUSED it rather than returning a
    #   plausible wrong number, which is the guard working exactly as intended.  Only `rnt` is overridden, to 2e-7,
    #   because the weak integral lives INSIDE the nucleus at r ~ 1e-4 a.u. where the default 2e-6 would be coarse.
    #
    setDefaults("print summary: open", "zzz-WeakInteractionEnhancement.sum")
    #
    P4f  = [Configuration("1s^2 2p"), Configuration("1s^2 3p"), Configuration("1s^2 4p")]
    sets = [ ("2p",          [Configuration("1s^2 2p")]),
             ("2p 3p",       [Configuration("1s^2 2p"), Configuration("1s^2 3p")]),
             ("2p 3p 4p",    P4f),
             ("2p 3p 4p 5p", vcat(P4f, [Configuration("1s^2 5p")])) ]
    #
    liPnc = function(Z, A, rr, model, oddCfg)
        nmX = Nuclear.Model(Z, model, A, rr, AngularJ64(0), 0.0, 0.0, 0.0)
        cfg = vcat([Configuration("1s^2 2s"), Configuration("1s^2 3s")], oddCfg)
        grX = Basics.recommendedGrid(cfg, nmX; rnt = 2.0e-7)
        setDefaults("standard grid", grX)
        mpX = SelfConsistent.performSCF(cfg, nmX, grX, AsfSettings(); printout=false)
        evX = filter(l -> l.parity == Basics.plus,  mpX.levels)
        odX = filter(l -> l.parity == Basics.minus, mpX.levels)
        sA  = evX[argmin([l.energy for l in evX])];   sB = evX[argmax([l.energy for l in evX])]
        pA  = odX[argmin([l.energy for l in odX])]
        ( imag(WeakInteractionEnhancement.computePncE1Amplitude(sB, sA, mpX, nmX, grX)), abs(sA.energy - pA.energy) )
    end
    #
    println("\n  (i)  E1_PNC for 1s^2 3s <- 1s^2 2s, as the intermediate set grows")
    trunc = Dict{String,Float64}();   final = Dict{String,Float64}()
    for (nme, Z, A, rr) in [("Li-like 40Ca", 20.0, 40.0, 3.4776), ("Li-like 208Pb", 82.0, 208.0, 5.5012)]
        println("\n     $nme")
        println("       intermediate set      E1_PNC [i e a_0]        change        2s-2p_1/2 [eV]")
        pv = 0.
        for (lbl, oc) in sets
            (a, dE) = liPnc(Z, A, rr, Nuclear.FermiNucleus(), oc)
            sa = pv == 0. ? "     --" : @sprintf("%+8.2f %%", 100*(a/pv - 1))
            println("       " * @sprintf("%-20s %+.8e     %s      %8.1f", lbl, a, sa,
                                          Defaults.convertUnits("energy: from atomic to eV", dE)))
            if  pv != 0.   trunc[nme] = abs(100*(a/pv - 1))   end
            pv = a
        end
        final[nme] = pv
    end
    println("\n     The 2s-2p_1/2 splitting is 36 eV in calcium and 235 eV in lead -- large, screening-dominated, and")
    println("     nothing like the near-degeneracy that makes the hydrogen-like case of branch c untrustworthy.  Note")
    println("     also that ONE intermediate again gives the WRONG SIGN, in both ions, as it did in branch b.")
    #
    println("\n  (ii) the isoelectronic scan: how E1_PNC/Q_W grows along the Li-like sequence")
    println("       Z        A       r [fm]     E1_PNC/Q_W            local exponent")
    Zs = [20.0, 35.0, 50.0, 65.0, 82.0];   vals = Float64[]
    for Z in Zs
        nmD = Nuclear.Model(Z)
        (a, dE) = liPnc(Z, nmD.mass, nmD.radius, Nuclear.FermiNucleus(), P4f)
        v = abs(a) / abs(WeakInteractionMoment.weakCharge(nmD))
        push!(vals, v)
        sa = length(vals) > 1 ?
             @sprintf("%.4f", (log(vals[end])-log(vals[end-1]))/(log(Z)-log(Zs[length(vals)-1]))) : "  --"
        println("       " * @sprintf("%5.1f   %6.1f    %6.3f     %.8e         %s", Z, nmD.mass, nmD.radius, v, sa))
    end
    nz = length(Zs);   lx = log.(Zs);   ly = log.(vals)
    slope = (nz*sum(lx.*ly) - sum(lx)*sum(ly)) / (nz*sum(lx.^2) - sum(lx)^2)
    println("       least-squares exponent = " * @sprintf("%.4f", slope))
    println("\n     About Z^3 on average, but the LOCAL exponent climbs from 2.46 to 4.70 across the sequence, so a single")
    println("     power law describes this badly.  Two Z-dependences are competing: the weak matrix element itself grows")
    println("     as Z^4 (branch e of example-Bb.jl), while the 2s-2p_1/2 denominator grows only as about Z^1.3, being")
    println("     screening-dominated at low Z and relativistic at high Z.  The upward curvature is the relativistic")
    println("     enhancement taking over, and it is why the heaviest systems are the interesting ones.")
    #
    println("\n  (iii) the density SHAPE, at fixed rms radius -- branch d repeated for these two systems")
    println("       system              Fermi                  uniform             spread [%]")
    shape = Dict{String,Float64}()
    for (nme, Z, A, rr) in [("Li-like 40Ca", 20.0, 40.0, 3.4776), ("Li-like 208Pb", 82.0, 208.0, 5.5012)]
        (fF, _) = liPnc(Z, A, rr, Nuclear.FermiNucleus(),   P4f)
        (fU, _) = liPnc(Z, A, rr, Nuclear.UniformNucleus(), P4f)
        shape[nme] = abs(100*(fU/fF - 1))
        println("       " * @sprintf("%-18s  %+.8e   %+.8e   %+9.3f", nme, fF, fU, 100*(fU/fF - 1)))
    end
    #
    println("\n  (iv) THE ERROR BUDGET, and it turns over with Z")
    println("       system              truncation      density shape      neutron skin")
    println("       " * @sprintf("%-18s  %8.2f %%      %8.2f %%       %8.2f %%", "Li-like 40Ca",
                                 trunc["Li-like 40Ca"], shape["Li-like 40Ca"], 0.04))
    println("       " * @sprintf("%-18s  %8.2f %%      %8.2f %%       %8.2f %%", "Li-like 208Pb",
                                 trunc["Li-like 208Pb"], shape["Li-like 208Pb"], 0.91))
    println("       (the skin column is the estimate of branch h in example-Bb.jl, carried over)")
    println("\n     THE DOMINANT UNCERTAINTY SWAPS.  In calcium the ELECTRONIC sum limits the answer and the nucleus is")
    println("     nearly irrelevant; in lead the NUCLEUS limits it and the electronic sum is the best-converged part of")
    println("     the calculation, since the higher np levels are pushed away and the low ones dominate more cleanly.")
    println("     That is a more useful thing to know than either amplitude on its own: it says that improving a light")
    println("     system means adding intermediate states, and improving a heavy one means improving the nuclear model,")
    println("     and that effort spent the other way round is wasted in both cases.")
    println("\n     WHAT MAY AND MAY NOT BE QUOTED.  The two amplitudes above are complete as ELECTRONIC calculations at")
    println("     the stated truncation.  They do NOT contain a neutron skin -- the weak charge is placed on the charge")
    println("     distribution, which branch h shows costs about -0.9 % in lead -- and they carry no QED.  For lead the")
    println("     honest total is therefore a few percent, dominated by the nuclear density, and it would be wrong to")
    println("     present the 8 digits printed above as anything but the arithmetic they are.")
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  22-Aug-2026
    #   [PROVENANCE: as branch a.]
    # Branch g: Ba+, AND THE FIRST REAL EXERCISE OF THE MANY-ELECTRON PATH.  Every verified branch before this one used
    #   a ONE-ELECTRON system on purpose, so that a disagreement would blame the operator rather than the atomic model.
    #   The cost of that discipline is that the spin-angular contraction has never been tested with anything to
    #   contract over.  Ba+ is the mildest possible first case: 55 electrons, but one valence electron over a closed
    #   Xe-like core.
    #
    #   THE SYSTEM IS ALSO THE ONE TWO SEPARATE EXPERIMENTAL GROUPS ARE BUILDING AROUND -- the 6s ^2S_1/2 to 5d ^2D_3/2
    #   transition, whose 80 s metastable lifetime is what makes the Fortson E1_pnc/E2 interference scheme possible.
    #
    #   THE SELECTION RULES SPLIT THE SUM INTO TWO DIFFERENT CHANNELS, which a one-electron test could not exhibit,
    #   because the weak charge is a RANK-0 operator and each term of the sum needs J_n equal to a different one of the
    #   two levels:
    #      term 1, weak acting on the INITIAL level:  <5d||D||n> <n|H_W|6s>   needs J_n = 1/2, kappa = +1  ->  np_1/2
    #      term 2, weak acting on the FINAL level:    <5d|H_W|n> <n||D||6s>   needs J_n = 3/2, kappa = -2  ->  np_3/2
    #   so the two fine-structure partners of 6p enter through DIFFERENT terms and may be compared.
    #
    setDefaults("print summary: open", "zzz-WeakInteractionEnhancement.sum")
    #
    nmBa = Nuclear.Model(56., Nuclear.FermiNucleus(), 138.0, 4.8378, AngularJ64(0), 0.0, 0.0, 0.0)
    cfgB = [Configuration("[Xe] 6s"), Configuration("[Xe] 5d"), Configuration("[Xe] 6p")]
    grB  = Basics.recommendedGrid(cfgB, nmBa; rnt = 2.0e-7)
    rNuc = Defaults.convertUnits("length: from fm to atomic", nmBa.radius)
    nIn  = count(r -> 0 < r <= rNuc, grB.r)
    setDefaults("standard grid", grB)
    #
    println("\n  (i)  one grid, two scales six orders of magnitude apart")
    println("       " * @sprintf("box = %.1f a.u. over %d points;  nuclear radius = %.3e a.u. with %d points inside it",
                                 grB.r[end], grB.NoPoints, rNuc, nIn))
    println("       The valence orbital reaches tens of a.u. while the weak integral lives ENTIRELY inside the nucleus,")
    println("       and both have to be resolved at once.  Basics.recommendedGrid sets the box; only rnt is overridden.")
    #
    mpB = SelfConsistent.performSCF(cfgB, nmBa, grB, AsfSettings(); printout=false)
    gsB = mpB.levels[argmin([l.energy for l in mpB.levels])]
    println("\n       level      J^P        excitation [eV]      experiment [eV]")
    for l in mpB.levels
        ex = Defaults.convertUnits("energy: from atomic to eV", l.energy - gsB.energy)
        sa = Basics.twice(l.J) == 3 && l.parity == Basics.plus  ? "0.604 (5d_3/2)" :
             Basics.twice(l.J) == 1 && l.parity == Basics.minus ? "2.512 (6p_1/2)" : ""
        println("       " * @sprintf("%4d      %-9s   %8.3f            %s", l.index, string(LevelSymmetry(l.J,l.parity)), ex, sa))
    end
    evB = filter(l -> l.parity == Basics.plus,  mpB.levels)
    odB = filter(l -> l.parity == Basics.minus, mpB.levels)
    s6  = evB[findfirst(l -> Basics.twice(l.J) == 1, evB)]
    d5  = evB[findfirst(l -> Basics.twice(l.J) == 3, evB)]
    p1  = odB[findfirst(l -> Basics.twice(l.J) == 1, odB)]
    p3  = odB[findfirst(l -> Basics.twice(l.J) == 3, odB)]
    #
    println("\n  (ii) THE CORE-SPECTATOR TEST -- branch c of example-Bb.jl, now with 54 spectators")
    println("       A closed core cannot contribute to a ONE-PARTICLE P-odd operator, so the 55-electron amplitude must")
    println("       equal the bare valence one EXACTLY.  Same identity as branch c, but the spin-angular contraction now")
    println("       runs over 55 electrons rather than one.  This is the check the earlier branches could not make.")
    rhoB = WeakInteractionMoment.nuclearDensity(nmBa, grB)
    QWB  = WeakInteractionMoment.weakCharge(nmBa)
    wPair = Tuple{Float64,Float64}[]
    for (fL, iL, shF, shI, nme) in [(p1, s6, Subshell(6, 1), Subshell(6,-1), "6p_1/2 <- 6s_1/2   (kappa +1 <- -1)"),
                                    (p3, d5, Subshell(6,-2), Subshell(5, 2), "6p_3/2 <- 5d_3/2   (kappa -2 <- +2)")]
        aM  = WeakInteractionMoment.weakChargeAmplitude(fL, iL, nmBa, grB)
        rad = WeakInteractionMoment.radialIntegralPQminus(rhoB, fL.basis.orbitals[shF], iL.basis.orbitals[shI], grB)
        exa = sqrt(Basics.twice(fL.J) + 1.0) * WeakInteractionMoment.GF/(2sqrt(2.0)) * QWB * im * rad
        push!(wPair, (imag(aM), imag(exa)))
        println("       " * nme)
        println("          " * @sprintf("55-electron amplitude  = %+.10e i", imag(aM)))
        println("          " * @sprintf("bare valence, by hand  = %+.10e i", imag(exa)))
        println("          " * @sprintf("RATIO (must be 1)      = %.12f", real(aM/exa)))
    end
    println("\n       Both to twelve digits.  The many-electron path is not merely untested any more: the closed core is")
    println("       correctly inert, and the angular machinery carries 54 spectator electrons without leaving a trace.")
    #
    println("\n  (iii) the two channels, and why one of them barely exists")
    ampB = WeakInteractionEnhancement.computePncE1Amplitude(d5, s6, mpB, nmBa, grB)
    for n in odB
        a = WeakInteractionEnhancement.computePncE1Amplitude(d5, s6, Multiplet("one", [n]), nmBa, grB)
        println("       " * @sprintf("n = %d [%-9s]  contribution = %+.8e i     share = %9.4f %%",
                                      n.index, string(LevelSymmetry(n.J,n.parity)), imag(a), 100*imag(a)/imag(ampB)))
    end
    println("       " * @sprintf("TOTAL  = %+.8e i  [e a_0]", imag(ampB)))
    println("       " * @sprintf("the two weak matrix elements differ by a factor %.3e", abs(wPair[2][1]/wPair[1][1])))
    println("\n       THE p_3/2 CHANNEL IS A MILLIONTH OF THE p_1/2 ONE, and that is physics rather than an accident.  The")
    println("       weak charge is a CONTACT interaction: it needs both orbitals inside the nucleus, where only s_1/2 and")
    println("       p_1/2 have appreciable amplitude, the others being held out by the centrifugal barrier.  A d_3/2-p_3/2")
    println("       pair barely penetrates at all.  This is why parity-violation experiments are built on transitions in")
    println("       which an s state carries the mixing, and it is visible here as six orders of magnitude.")
    #
    println("\n  (iv) the exchange symmetry, which is the OPPOSITE of branch a and must be")
    revB = WeakInteractionEnhancement.computePncE1Amplitude(s6, d5, mpB, nmBa, grB)
    println("       " * @sprintf("re/|amp| = %.2e   (purely imaginary)", abs(real(ampB))/abs(ampB)))
    println("       " * @sprintf("|<f||i> + <i||f>| / |amp| = %.3e", abs(ampB+revB)/abs(ampB)))
    println("       " * @sprintf("|<f||i> - <i||f>| / |amp| = %.3e     <-- THIS one must vanish here", abs(ampB-revB)/abs(ampB)))
    println("\n       Hermiticity gives <f||T||i> = (-1)^(J_f-J_i) conj(<i||T||f>).  For a purely imaginary amplitude that")
    println("       is ANTIsymmetric when J_f = J_i -- branch a, where the SUM vanished -- and SYMMETRIC when the two J")
    println("       differ by one, as here, where the DIFFERENCE vanishes instead.  The two branches together test the")
    println("       phase factor itself, which neither can do alone.  Worth recording that this branch was first written")
    println("       asserting the sum would vanish, and the calculation said 2.0; the calculation was right.")
    #
    println("\n  (v) the E2 amplitude that the experiments interfere against")
    e2B = MultipoleMoment.emmStaticAmplitude(2, d5, s6, grB) * sqrt(Basics.twice(d5.J) + 1.0)
    println("       " * @sprintf("<5d_3/2||T^(E2)||6s_1/2>  = %+.6e  [e a_0^2]", real(e2B)))
    println("       " * @sprintf("|E1_PNC / E2|             = %.6e", abs(ampB)/abs(real(e2B))))
    #
    println("\n  (vi) WHAT THIS NUMBER IS AND IS NOT.  The structural results above are exact and stand on their own.")
    println("       The AMPLITUDE does not: this is a single configuration per level with no correlation whatever, and the")
    println("       level energies show what that costs -- 5d_3/2 comes out at 0.793 eV against an experimental 0.604,")
    println("       high by 31 %, and 6p_1/2 at 2.321 against 2.512, low by 8 %.  Both sit in denominators of the sum.")
    println("       Adding 7p and 8p moves the amplitude by a further +6.0 % and +2.2 %, so the intermediate set is not")
    println("       converged either, and the largest set already strains the box, which Basics.recommendedGrid pushes")
    println("       past 250 a.u. to accommodate a diffuse 8p while 6s then sits near the edge of what it can carry.")
    println("       So: do NOT read the digits as a prediction.  What has been established is that the machinery runs on")
    println("       a real many-electron system, that the core is correctly inert to twelve digits, and that the two")
    println("       channels behave as the contact nature of the weak interaction requires.  A competitive Ba+ number")
    println("       needs correlation, and should be compared against the published coupled-cluster values rather than")
    println("       against anything in this file.")
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  22-Aug-2026
    #   [PROVENANCE: as branch a.]
    # Branch h: THE MISSING LAMB SHIFT, PUT IN WITH JAC'S OWN QED MODEL -- and branch c turns out to have understated
    #   the problem rather than overstated it.
    #
    #   Branch c estimated the missing 2s-2p_1/2 splitting by scaling the measured hydrogen Lamb shift as Z^4, and
    #   reported a factor of 379 at Z = 20.  Two things are done better here.  FIRST, the splitting is no longer
    #   estimated at all: `AsfSettings(...; qedModel = QedPetersburg())` puts the Shabaev et al. (2013) model
    #   Hamiltonian into the CI, and the corrected level energies flow into the denominators of the sum by themselves.
    #   SECOND, and this is the correction, branch c looked at ONE denominator of a sum that has TWO.
    #
    #   EVERY ns - np_1/2 PAIR OF THE SAME n IS DEGENERATE IN DIRAC THEORY, not only 2s - 2p_1/2.  The PNC sum runs
    #      E1_PNC = SUM_n [ <f||D||n> <n|H_W|i>/(E_i - E_n)  +  <f|H_W|n> <n||D||i>/(E_f - E_n) ]
    #   so an intermediate np_1/2 is near-degenerate with the INITIAL ns through the first denominator and with the
    #   FINAL n's through the second.  For 3s <- 2s the 2p_1/2 term is pathological through E_i - E_n and the 3p_1/2
    #   term is pathological through E_f - E_n -- and those two carry -39 % and +139 % of the total between them.  The
    #   whole amplitude rests on denominators a Dirac calculation gets wrong, not merely one contribution to it.
    #
    #   THE QED OUTPUT IS SUPPRESSED with redirect_stdout(devnull): the model prints every single-electron QED strength
    #   it forms, some five hundred lines per call, which would bury the result.
    #
    setDefaults("print summary: open", "zzz-WeakInteractionEnhancement.sum")
    #
    cfgQ = [Configuration("2s"), Configuration("3s"), Configuration("2p"), Configuration("3p"), Configuration("4p")]
    runQ = function(Z, useQed)
        nmQ = Nuclear.Model(Z)
        stQ = AsfSettings(AsfSettings(); scField = Basics.NuclearField())
        useQed  &&  (stQ = AsfSettings(stQ; qedModel = QedPetersburg()))
        grQ = Basics.recommendedGrid(cfgQ, nmQ; rnt = 2.0e-7);   setDefaults("standard grid", grQ)
        mpQ = redirect_stdout(devnull) do
            SelfConsistent.performSCF(cfgQ, nmQ, grQ, stQ; printout=false)
        end
        evQ = filter(l -> l.parity == Basics.plus,  mpQ.levels)
        odQ = filter(l -> l.parity == Basics.minus, mpQ.levels)
        ( evQ[argmin([l.energy for l in evQ])], evQ[argmax([l.energy for l in evQ])], odQ, mpQ, nmQ, grQ )
    end
    #
    println("\n  (i)  the 2s - 2p_1/2 splitting, three ways")
    println("       Z      Dirac only [a.u.]   QedPetersburg [a.u.]   ratio     1.6083e-7 Z^4     scaling/Qed")
    for Z in [20.0, 30.0, 50.0, 70.0, 92.0]
        (sN, _, odN, _, _, _) = runQ(Z, false);   pN = odN[argmin([l.energy for l in odN])]
        (sQ, _, odQ, _, _, _) = runQ(Z, true);    pQ = odQ[argmin([l.energy for l in odQ])]
        dN  = abs(sN.energy - pN.energy);   dQ = abs(sQ.energy - pQ.energy);   est = 1.6083e-7 * Z^4
        println("       " * @sprintf("%5.1f    %.6e        %.6e      %7.1f    %.6e      %7.2f",
                                     Z, dN, dQ, dQ/dN, est, est/dQ))
    end
    println("\n       The hydrogen-scaled estimate of branch c EXCEEDS the model Hamiltonian by three to eight times, and")
    println("       that is expected rather than alarming: the 2s Lamb shift carries a Bethe logarithm which is large at")
    println("       Z = 1 and falls steadily with Z, so a pure Z^4 extrapolation from hydrogen must overshoot.  The two")
    println("       independent estimates bracket the truth; neither is the truth.")
    #
    println("\n  (ii) BOTH denominators at Z = 20 -- the correction to branch c")
    for useQed in [false, true]
        (s2, s3, od, mp, nm, gr) = runQ(20.0, useQed)
        tot = imag(WeakInteractionEnhancement.computePncE1Amplitude(s3, s2, mp, nm, gr))
        println("       " * (useQed ? "QedPetersburg" : "Dirac only   ") * @sprintf("     TOTAL = %+.6e i", tot))
        println("         n              E_i - E_n       E_f - E_n      contribution         share")
        for n in od
            a = imag(WeakInteractionEnhancement.computePncE1Amplitude(s3, s2, Multiplet("one",[n]), nm, gr))
            a == 0. && continue
            println("         " * @sprintf("%d [%-7s]   %+.4e    %+.4e    %+.6e    %8.2f %%", n.index,
                                           string(LevelSymmetry(n.J,n.parity)), s2.energy-n.energy, s3.energy-n.energy,
                                           a, 100*a/tot))
        end
    end
    println("\n       Read the two denominator columns together.  Without QED the 2p_1/2 row has 6.8e-05 in the FIRST")
    println("       column and the 3p_1/2 row has 2.0e-05 in the SECOND -- two different near-degeneracies, one per term")
    println("       of the sum, carrying -39 % and +139 % of the total.  QED lifts both by about 110.")
    #
    println("\n  (iii) WHAT IT DOES TO THE AMPLITUDE, which is the error bar branch c could not state")
    println("       Z      E1_PNC, Dirac only     E1_PNC with QED       shift")
    for Z in [20.0, 30.0, 50.0, 70.0, 92.0]
        (sN, tN, _, mN, nN, gN) = runQ(Z, false)
        (sQ, tQ, _, mQ, nQ, gQ) = runQ(Z, true)
        aN = imag(WeakInteractionEnhancement.computePncE1Amplitude(tN, sN, mN, nN, gN))
        aQ = imag(WeakInteractionEnhancement.computePncE1Amplitude(tQ, sQ, mQ, nQ, gQ))
        println("       " * @sprintf("%5.1f    %+.8e       %+.8e      %+8.2f %%", Z, aN, aQ, 100*(aQ/aN - 1)))
    end
    println("\n       A PLAIN DIRAC CALCULATION OF H-LIKE PNC IS TOO LARGE BY A FACTOR OF 110 AT Z = 20, not by a")
    println("       percentage.  The error falls with Z as the finite-nuclear-size part of the splitting grows to")
    println("       compete with the QED part, reaching -42 % at Z = 92, where the nucleus does much of the splitting")
    println("       by itself.  Branch c called the hydrogen-like case not quantitative; it is worse than that at low Z")
    println("       and merely bad at high Z.")
    println("\n       AND THE CORRECTED NUMBER IS NOT A RESULT EITHER.  QedPetersburg is a MODEL HAMILTONIAN, not a")
    println("       rigorous bound-state QED evaluation, and the two estimates of the same splitting used above differ")
    println("       by three to eight times.  What may honestly be said is a DIRECTION and an ORDER: the Dirac amplitude")
    println("       is far too large, by about two orders of magnitude at Z = 20 and a factor of two at Z = 92.  A")
    println("       hydrogen-like PNC amplitude worth quoting needs the Lamb shift taken from a dedicated QED")
    println("       calculation for the specific ion, and the Li-like route of branch f -- where screening does the")
    println("       splitting and QED is a correction rather than the whole effect -- remains the better one.")
    #
    setDefaults("print summary: close", "")
    #
end
