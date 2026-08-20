println("Dpnew) Apply & test the GeneralizedOscillatorStrength module: the generalized oscillator")
println("    strength f_n(K) of the Bethe theory of fast charged-particle collisions, cf.")
println("    M. Inokuti, Rev. Mod. Phys. 43 (1971) 297, computed between the levels of an initial-")
println("    and a final-state multiplet.  f_n(K) is the inelastic counterpart of the atomic form")
println("    factor of module-FormFactor.jl (example-Ce.jl): where the form factor is the DIAGONAL")
println("    matrix element of exp(i K.r) for one level, the GOS is the OFF-DIAGONAL one between")
println("    two levels, and it reduces to the ordinary (optical) oscillator strength as K --> 0.")
println("    This version covers bound-bound transitions; the continuum density df(K,E)/dE, and")
println("    with it the full Bethe surface, the Bethe sum rule and the first-Born cross sections,")
println("    are not (yet) included.")

if  true
    # Last successful:  20-Aug-2026
    #   [PROVENANCE: produced on Julia 1.10.9, with a Manifest re-resolved for 1.10 because the checkout's Manifest had been
    #    resolved under 1.12.6 on the maintainer machine.  Every figure quoted below is an ABSOLUTE comparison against a
    #    closed-form result, a published figure or JAC's own independent path -- not a tolerance -- so a genuine change would
    #    be conspicuous rather than marginal.  Re-run on the maintainer machine to confirm; there is nothing to re-derive.]
    # Branch a: ATOMIC HYDROGEN against closed-form theory -- the acceptance test of the module, and
    #   entirely parameter-free.  Inokuti quotes the exact non-relativistic GOS of H in closed form,
    #   and two of those expressions are used here (q = (K a_0)^2 throughout):
    #
    #     Eq. (4.2)   f_2p(K)    = 2^13 * 3^-9 * [1 + (4/9) q]^-6
    #     Eq. (3.9)   f_(n=2)(K) = 24576 * [9 + 4q]^-5
    #     difference  f_2s(K)    = 98304 * q * [9 + 4q]^-6
    #
    #   NOTE on Eq. (4.2): the scanned paper prints the prefactor as 2^18 * 3^-9 = 13.32, which is an
    #   OCR artefact -- it must be 2^13 * 3^-9 = 0.41620, as Fig. 4 of the same paper confirms with its
    #   "opt (0.4162)" label, and as Eq. (4.5), M_2p^2 = 2^15 3^-10 = 0.55493 = f_2p R/E_2p, requires.
    #   The three formulae were checked against each other before use: f_2p + f_2s = f_(n=2) holds to
    #   machine precision at every q, and Eq. (3.9) reproduces the optical limits printed on Figs. 4-7 of
    #   the paper for n = 2,3,4,5 (0.41620, 7.9102e-2, 2.899e-2, 1.3938e-2).
    #
    #   TWO channels are exercised at once, which is the reason for including 2s as well as 2p:
    #     L = 1  carries 1s -> 2p, the optically allowed case;
    #     L = 0  carries 1s -> 2s, the MONOPOLE case.  j_0(Kr) = sin(Kr)/(Kr) is not the identity, so this
    #            matrix element is non-zero for K > 0 and collapses to the (vanishing) overlap only as
    #            K --> 0 -- which is exactly why f_2s(K) starts at zero and rises.  A GOS module that
    #            silently dropped L = 0 would look perfectly healthy on 1s -> 2p and be wrong here.
    #   Relativistically the single non-relativistic 2p is split into 2p_1/2 and 2p_3/2; their GOS must be
    #   SUMMED to recover f_2p(K), and their ratio must be exactly 1 : 2 by statistical weight at every K.
    #
    #   Basics.NuclearField() is essential and not a detail: with the default DFS field a ONE-electron
    #   system acquires a spurious self-interaction, which moved the Lyman-alpha energy from 10.2043 eV to
    #   10.0795 eV, split the two E1 gauges by a factor of 1.93, and -- because initial and final orbitals
    #   then solve DIFFERENT one-body operators -- destroyed the orthogonality that the L = 0 channel lives
    #   on, sending f_2s(K --> 0) to 1.5e3 instead of to zero.  With the bare -Z/r field both multiplets
    #   diagonalise the same operator, so the orbitals are orthogonal across the two SCF runs by
    #   construction, and all of that disappears.
    #
    #   RESULT, 20-Aug-2026, rbox = 40 a.u.:  agreement with the exact formulae to 2-5 significant digits
    #   over five decades of f, for BOTH channels --
    #      (Ka0)^2      f_2p JAC/exact     f_2s JAC/exact
    #        1e-4          0.999977           0.999980
    #        0.09          0.999981           0.999985
    #        1.0           1.000016           1.000032
    #        9.0           1.000156           1.000208
    #   The residual 2e-5 at small K is the relativistic correction, of order (alpha Z)^2 = 5e-5 at Z = 1,
    #   which the exact NON-relativistic formulae cannot contain; the slow drift to 2e-4 at (Ka0)^2 = 9 is
    #   the radial grid beginning to under-resolve the oscillations of j_L(Kr).  Neither is a defect.
    #   Fine-structure split at (Ka0)^2 = 1:  2p_1/2 = 1.5275103e-2, 2p_3/2 = 3.0549698e-2, ratio 2.000000.
    setDefaults("print summary: open", "zzz-GeneralizedOscillatorStrength.sum")
    #
    f2pExact(q) = 2.0^13 * 3.0^-9 * (1 + (4/9)*q)^-6
    f2sExact(q) = 98304 * q * (9 + 4q)^-6.0
    #
    grid   = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 40.0)
    setDefaults("standard grid", grid)
    nModel = Nuclear.Model(1., UniformNucleus(), 1., 0.8783, AngularJ64(1//2), 2.7928473, 0.0, 0.0)
    asfH   = AsfSettings(AsfSettings(); scField = Basics.NuclearField())
    Ks     = [0.01, 0.1, 0.3, 0.5, 1.0, 1.5, 2.0, 3.0]
    #
    gosSettings = GeneralizedOscillatorStrength.Settings(GeneralizedOscillatorStrength.Settings();
                      qValues = Ks, calcOpticalLimit = true, printBefore = true)
    comp = Atomic.Computation(Atomic.Computation(), name = "H 1s -> 2s, 2p generalized oscillator strengths",
              grid = grid, nuclearModel = nModel,
              initialConfigs = [Configuration("1s")],                        initialAsfSettings = asfH,
              finalConfigs   = [Configuration("2s"), Configuration("2p")],   finalAsfSettings   = asfH,
              processSettings = gosSettings)
    wb    = perform(comp; output = true)
    lines = wb["generalized oscillator strengths:"]
    #
    # Sum the two 2p_j components; the 2s line is the even-parity (monopole) one
    f2p = zeros(length(Ks));    f2s = zeros(length(Ks));    fsSplit = Dict{String,Vector{Float64}}()
    for  line in lines
        if  line.finalLevel.parity == Basics.minus
            f2p .+= line.gosValues
            fsSplit[string(LevelSymmetry(line.finalLevel.J, line.finalLevel.parity))] = line.gosValues
        else
            f2s .+= line.gosValues
        end
    end
    println("\n\n  Hydrogen GOS against the exact Bethe formulae of Inokuti (1971):\n")
    println("   (Ka0)^2      f_2p JAC       f_2p exact     ratio        f_2s JAC       f_2s exact     ratio")
    for  (i,K) in enumerate(Ks)
        q = K*K
        println("  " * rpad(round(q, sigdigits=5), 11) * "  " * rpad(round(f2p[i], sigdigits=7), 14) * " " *
                rpad(round(f2pExact(q), sigdigits=7), 14) * " " * rpad(round(f2p[i]/f2pExact(q), digits=6), 12) * " " *
                rpad(round(f2s[i], sigdigits=7), 14) * " " * rpad(round(f2sExact(q), sigdigits=7), 14) * " " *
                string(round(f2s[i]/f2sExact(q), digits=6)))
    end
    println("\n  Fine-structure split of f_2p (must be exactly 1 : 2 at every K):")
    for  (sym, vals) in fsSplit
        println("     " * rpad(sym, 8) * " : " * string(round.(vals, sigdigits=7)))
    end
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  20-Aug-2026
    #   [PROVENANCE: produced on Julia 1.10.9, with a Manifest re-resolved for 1.10 because the checkout's Manifest had been
    #    resolved under 1.12.6 on the maintainer machine.  Every figure quoted below is an ABSOLUTE comparison against a
    #    closed-form result, a published figure or JAC's own independent path -- not a tolerance -- so a genuine change would
    #    be conspicuous rather than marginal.  Re-run on the maintainer machine to confirm; there is nothing to re-derive.]
    # Branch b: HELIUM 1^1S -> 2^1P, the best-studied generalized oscillator strength there is, and the one
    #   to which Inokuti devotes a whole "case history" (Sec. 3.4, p. 317).  Two things are checked, and
    #   they are checks of DIFFERENT things -- which is the point of the branch.
    #
    #   (i) IS THE MODULE RIGHT?  The optical limit of the GOS is, by construction, the LENGTH-form
    #       oscillator strength: the operator exp(i K.r) is a function of position only, so no gauge
    #       freedom enters, and as K --> 0 it must reproduce the Babushkin (length) result of
    #       PhotoExcitation for the very same orbitals.  Measured here: 0.42148351 from the GOS module
    #       against 0.42146859 from PhotoExcitation, a ratio of 1.000035.  Branch a validated the module
    #       on a ONE-electron system; this is the corresponding statement for a many-electron one, where
    #       the spin-angular recoupling of two CSF bases actually does some work.
    #
    #   (ii) IS THE ATOM RIGHT?  Not nearly, and that is Inokuti's point.  The exact optical value is
    #       f = 0.27616 +/- 0.00001 (Schiff & Pekeris 1964); a single-configuration Dirac-Fock helium
    #       gives 0.4215, i.e. 53% too large, and its own two E1 gauges disagree by 11% (0.3780 Coulomb
    #       against 0.4215 Babushkin) -- the standard signature of orbitals that are not good enough. So
    #       the ABSOLUTE GOS of He is limited by the atomic model, not by the Bethe machinery.
    #       An ad-hoc correlation list ([1s^2, 2s^2, 2p^2] -> [1s2p, 2s2p, 2p3d]) was tried and made the
    #       optical limit WORSE, 0.5247; it is not kept here, and is recorded only so that the next reader
    #       does not repeat it.  Reaching 0.27616 is a genuine atomic-structure exercise (an optimised
    #       correlated model for both states), not a GOS one.
    #
    #   What survives the model error is the SHAPE.  f(K)/f(0) is far less sensitive to correlation than
    #   f(0) itself, and it is compared with Lassettre's truncated series, Eq. (3.44),
    #        f(K)/f = (1+x)^-6 [1 + g x/(1+x)],   x = (Ka0)^2/3.391,
    #   here with g = 0.  Measured ratio of the two: 1.000 at (Ka0)^2 = 1e-4, 0.999 at 0.04, 1.003 at 0.25,
    #   1.031 at 0.64, 1.075 at 1.0.  Inokuti quotes the Kim-Inokuti (1968) calculation as good to about
    #   1% for (Ka0)^2 <= 2, and the g = 0 truncation is itself only a fit, so agreement at the few-percent
    #   level up to (Ka0)^2 ~ 1 is the expected outcome and is what is seen.
    setDefaults("print summary: open", "zzz-GeneralizedOscillatorStrength.sum")
    #
    lassettre(q, g) = (1 + q/3.391)^-6 * (1 + g*(q/3.391)/(1 + q/3.391))
    #
    grid   = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 25.0)
    setDefaults("standard grid", grid)
    nModel = Nuclear.Model(2.)
    Ks     = [0.01, 0.2, 0.5, 0.8, 1.0, 1.4, 2.0]
    #
    gosSettings = GeneralizedOscillatorStrength.Settings(GeneralizedOscillatorStrength.Settings();
                      qValues = Ks, calcOpticalLimit = true, printBefore = true)
    comp = Atomic.Computation(Atomic.Computation(), name = "He 1^1S -> 2^1P generalized oscillator strength",
              grid = grid, nuclearModel = nModel,
              initialConfigs = [Configuration("1s^2")], finalConfigs = [Configuration("1s 2p")],
              processSettings = gosSettings)
    wb    = perform(comp; output = true)
    lines = wb["generalized oscillator strengths:"]
    #
    # the same transition through PhotoExcitation, for the gauge comparison of point (i)
    peSettings = PhotoExcitation.Settings(PhotoExcitation.Settings(); multipoles = [E1],
                     gauges = [UseCoulomb, UseBabushkin], printBefore = false)
    comp2 = Atomic.Computation(Atomic.Computation(), name = "He 1^1S -> 2^1P oscillator strengths",
                grid = grid, nuclearModel = nModel,
                initialConfigs = [Configuration("1s^2")], finalConfigs = [Configuration("1s 2p")],
                processSettings = peSettings)
    peLines = perform(comp2; output = true)["photo-excitation lines:"]
    #
    cands = filter(l -> Basics.twice(l.finalLevel.J) == 2  &&  l.finalLevel.parity == Basics.minus, lines)
    line  = cands[ argmax([l.opticalLimit for l in cands]) ]
    println("\n\n  (i) module check -- optical limit of the GOS against JAC's own oscillator strengths:\n")
    for  pl in peLines
        if  pl.finalLevel.index == line.finalLevel.index
            println("     GOS optical limit    = $(line.opticalLimit)")
            println("     PhotoExc f, Coulomb   = $(pl.oscStrength.Coulomb)")
            println("     PhotoExc f, Babushkin = $(pl.oscStrength.Babushkin)   <-- the length form the GOS must match")
            println("     GOS / Babushkin       = $(line.opticalLimit/pl.oscStrength.Babushkin)")
            println("     Coulomb / Babushkin   = $(pl.oscStrength.Coulomb/pl.oscStrength.Babushkin)   (gauge spread = orbital quality)")
        end
    end
    println("\n  (ii) atom check -- optical limit against Schiff-Pekeris 0.27616:  " *
            "ratio = $(line.opticalLimit/0.27616)")
    println("\n  shape of the GOS against the Lassettre series, Eq. (3.44) with g = 0:\n")
    println("   (Ka0)^2     f_n(K) JAC     f(K)/f JAC     Lassettre g=0     ratio")
    for  (i,K) in enumerate(Ks)
        q = K*K;    rel = line.gosValues[i]/line.opticalLimit
        println("  " * rpad(round(q, sigdigits=5), 10) * "  " * rpad(round(line.gosValues[i], sigdigits=7), 14) *
                " " * rpad(round(rel, sigdigits=7), 14) * " " * rpad(round(lassettre(q, 0.0), sigdigits=7), 17) *
                " " * string(round(rel/lassettre(q, 0.0), digits=5)))
    end
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  20-Aug-2026
    #   [PROVENANCE: produced on Julia 1.10.9, with a Manifest re-resolved for 1.10 because the checkout's Manifest had been
    #    resolved under 1.12.6 on the maintainer machine.  Every figure quoted below is an ABSOLUTE comparison against a
    #    closed-form result, a published figure or JAC's own independent path -- not a tolerance -- so a genuine change would
    #    be conspicuous rather than marginal.  Re-run on the maintainer machine to confirm; there is nothing to re-derive.]
    # Branch c: THE GOS MINIMUM -- the sharpest test in this file.  Inokuti's Fig. 12 shows f_n(K) for
    #   He 3^1P <- 2^1S passing through an EXACT ZERO within the first Born approximation, near
    #   ln(Ka0)^2 ~ -2, i.e. (Ka0)^2 ~ 0.14.  This demands more of a GOS code than any smooth curve can:
    #   the amplitude must change SIGN, so the relative sign of competing radial contributions has to be
    #   right, not merely their magnitudes.  Inokuti (p. 319) states the physical reading: "the position
    #   of the minimum is related to the node of the orbitals active in the transition and therefore
    #   provides a stringent test of calculated wave functions".
    #
    #   The initial state is the METASTABLE 2^1S, so the box must hold a 3p orbital of a near-neutral
    #   system: with Z_eff ~ 1 the outer turning point of (n=3, l=1) is r_+ = 16.9 a.u., so rbox = 50 a.u.
    #   is used here rather than the 25 a.u. of branch b (Rule 12).  Note that 1s2s carries BOTH 2^3S_1
    #   and 2^1S_0, so the initial level is selected by J = 0.
    #
    #   RESULT, 20-Aug-2026:  optical limit f = 0.23395; the L = 1 amplitude changes sign between
    #   (Ka0)^2 = 0.1225 (+1.44273e-2) and (Ka0)^2 = 0.16 (-1.27345e-2), placing the zero at
    #   (Ka0)^2 = 0.1424 by linear interpolation, against the ~0.14 of Inokuti Fig. 12.  Across the zero
    #   f_n(K) falls from 0.21787 at (Ka0)^2 = 2.5e-3 to 5.1766e-4 at 0.16, a drop by a factor of 421, and
    #   then rises again to 9.603e-3 at 0.36 -- the characteristic trough-and-recovery of Fig. 12.  (The
    #   overall sign of the amplitude is an arbitrary phase; only the crossing matters.)
    #   That the minimum is not numerically zero is
    #   expected and physical: this is a relativistic calculation, in which (Inokuti, p. 320) the zero of
    #   each j occurs at a slightly different (Ka0)^2, so the summed GOS "can never quite vanish".
    setDefaults("print summary: open", "zzz-GeneralizedOscillatorStrength.sum")
    #
    grid   = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 50.0)
    setDefaults("standard grid", grid)
    nModel = Nuclear.Model(2.)
    Ks     = [0.05, 0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.8, 1.2]
    #
    gosSettings = GeneralizedOscillatorStrength.Settings(GeneralizedOscillatorStrength.Settings();
                      qValues = Ks, calcOpticalLimit = true, printBefore = true)
    comp = Atomic.Computation(Atomic.Computation(), name = "He 2^1S -> 3^1P GOS minimum",
              grid = grid, nuclearModel = nModel,
              initialConfigs = [Configuration("1s 2s")], finalConfigs = [Configuration("1s 3p")],
              processSettings = gosSettings)
    wb    = perform(comp; output = true)
    # the 2^1S_0 initial level, and among its J=1 odd final levels the 3^1P_1 (much the largest strength)
    cands = filter(l -> Basics.twice(l.initialLevel.J) == 0  &&  Basics.twice(l.finalLevel.J) == 2  &&
                        l.finalLevel.parity == Basics.minus, wb["generalized oscillator strengths:"])
    line  = cands[ argmax([l.opticalLimit for l in cands]) ]
    #
    amps = Float64[]
    for  K in Ks
        ch = filter(c -> c.q == K  &&  c.L == 1, line.channels)
        push!(amps, length(ch) > 0 ? real(ch[1].amplitude) : 0.0)
    end
    println("\n\n  He 2^1S -> 3^1P :  optical limit f = $(line.opticalLimit),  dE = $(line.deltaEnergy) a.u.\n")
    println("   (Ka0)^2        f_n(K)           amplitude(L=1)")
    for  (i,K) in enumerate(Ks)
        println("  " * rpad(round(K*K, sigdigits=5), 12) * "  " * rpad(round(line.gosValues[i], sigdigits=7), 16) *
                " " * string(round(amps[i], sigdigits=6)))
    end
    # locate the sign change of the amplitude and interpolate the zero in (Ka0)^2
    for  i = 1:length(Ks)-1
        if  amps[i]*amps[i+1] < 0
            q1 = Ks[i]^2;   q2 = Ks[i+1]^2
            qz = q1 + (q2 - q1) * abs(amps[i]) / (abs(amps[i]) + abs(amps[i+1]))
            println("\n     amplitude changes SIGN between (Ka0)^2 = $q1 and $q2")
            println("     interpolated zero at (Ka0)^2 = $(round(qz, sigdigits=4))" *
                    "     [Inokuti Fig. 12:  about 0.14]")
            println("     f_n at the shallowest scanned point = $(minimum(line.gosValues[1:i+2]))")
        end
    end
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  20-Aug-2026
    #   [PROVENANCE: produced on Julia 1.10.9, with a Manifest re-resolved for 1.10 because the checkout's Manifest had been
    #    resolved under 1.12.6 on the maintainer machine.  Every figure quoted below is an ABSOLUTE comparison against a
    #    closed-form result, a published figure or JAC's own independent path -- not a tolerance -- so a genuine change would
    #    be conspicuous rather than marginal.  Re-run on the maintainer machine to confirm; there is nothing to re-derive.]
    # Branch d: THE GOS THAT ONLY A RELATIVISTIC CODE HAS.  Inokuti notes (p. 303) that the matrix element
    #   eps_n(K) "vanishes for all transitions between states with different spin multiplicities and so
    #   does f_n(K) ... When one takes into account spin-orbit coupling, f_n(K) for such a transition is
    #   finite, as is especially true with heavy atoms."  In LS coupling the intercombination GOS is
    #   therefore identically zero at EVERY K; in a jj-coupled Dirac calculation it is finite, and it must
    #   grow steeply with Z as the singlet-triplet mixing grows.  No non-relativistic GOS code can produce
    #   this branch at all.
    #
    #   He-like ions are the clean testing ground: 1s^2 ^1S_0 -> 1s2p, whose two J=1 final levels are the
    #   resonance line w (1^1P_1) and the intercombination line y (2^3P_1).  The branch computes the GOS
    #   of both at two widely separated Z and reports the RATIO y/w, which is the quantity that carries
    #   the relativistic mixing.  Boxes are matched to the orbitals (Rule 12): the 2p turning point scales
    #   as 1/Z_eff ~ 1/(Z-1), giving r_+ = 0.40 a.u. at Z = 18 and 0.13 a.u. at Z = 54.
    #
    #   RESULT, 20-Aug-2026:
    #        Z = 18 (Ar XVII):  f_w(0) = 0.79124,  f_y(0) = 1.02802e-2,   y/w = 0.01299
    #        Z = 54 (Xe LIII):  f_w(0) = 0.55026,  f_y(0) = 0.21628,      y/w = 0.39306
    #   a THIRTY-fold rise of a quantity that LS coupling puts at exactly zero, which is the effect Inokuti
    #   describes.  The ratio is nearly independent of K at small momentum transfer (0.01299 -> 0.01341
    #   over (Ka0)^2 = 0.25 ... 400 at Z = 18), as it should be: the singlet-triplet mixing coefficient
    #   factorises out of the amplitude there, and only at large K, where the two transitions weight their
    #   radial integrands differently, does the ratio start to drift upwards.
    #   NOTE what is and is not claimed.  The Z-trend and its physical origin are the result; the two sets
    #   of optical f-values have NOT been compared against tabulated He-like w/y oscillator strengths, and
    #   doing so would be the natural next step for this branch.
    setDefaults("print summary: open", "zzz-GeneralizedOscillatorStrength.sum")
    #
    for  (Z, rbox, Ks)  in  [ (18.0, 3.0, [0.5, 2.0, 5.0, 10.0, 20.0]),
                              (54.0, 1.2, [2.0, 5.0, 20.0, 50.0, 80.0]) ]
        println("\n\n>>>>>>  He-like Z = $Z :  1s^2 ^1S_0 -> 1s2p (resonance w and intercombination y)  <<<<<<\n")
        gridZ = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6/Z, h = 5.0e-2, hp = 1.0e-2, rbox = rbox)
        setDefaults("standard grid", gridZ)
        gosSetZ = GeneralizedOscillatorStrength.Settings(GeneralizedOscillatorStrength.Settings();
                          qValues = Ks, calcOpticalLimit = true, printBefore = true)
        compZ = Atomic.Computation(Atomic.Computation(), name = "He-like Z=$Z intercombination GOS",
                  grid = gridZ, nuclearModel = Nuclear.Model(Z),
                  initialConfigs = [Configuration("1s^2")], finalConfigs = [Configuration("1s 2p")],
                  processSettings = gosSetZ)
        wbZ   = perform(compZ; output = true)
        linesZ = filter(l -> Basics.twice(l.finalLevel.J) == 2, wbZ["generalized oscillator strengths:"])
        if  length(linesZ) < 2   println("  fewer than two J=1 final levels found");  continue   end
        # the resonance line w has by far the larger optical limit, the intercombination line y the smaller
        lw = linesZ[ argmax([l.opticalLimit for l in linesZ]) ]
        ly = linesZ[ argmin([l.opticalLimit for l in linesZ]) ]
        println("  resonance w:         optical limit = $(lw.opticalLimit),   dE = $(lw.deltaEnergy) a.u.")
        println("  intercombination y:  optical limit = $(ly.opticalLimit),   dE = $(ly.deltaEnergy) a.u.")
        println("  y/w at K --> 0:  $(ly.opticalLimit/lw.opticalLimit)")
        println("\n   (Ka0)^2        f_w(K)           f_y(K)           y/w")
        for  (i,K) in enumerate(Ks)
            println("  " * rpad(round(K*K, sigdigits=5), 12) * "  " * rpad(round(lw.gosValues[i], sigdigits=7), 16) *
                    " " * rpad(round(ly.gosValues[i], sigdigits=7), 16) * " " *
                    string(round(ly.gosValues[i]/lw.gosValues[i], sigdigits=5)))
        end
    end
    #
    setDefaults("print summary: close", "")
    #
end
