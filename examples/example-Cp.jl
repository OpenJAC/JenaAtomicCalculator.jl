
println("Cp) Apply & test the NuclearShielding module: how the electrons of a closed-shell ion modify, at the")
println("    nucleus, an electric field gradient applied from outside.")

using Printf

# WRITTEN 29-Sep-2026, together with module-NuclearShielding.jl.
#
# WHAT THE MODULE ANSWERS, in one sentence: what does the nucleus actually FEEL, given what is applied from outside?
# It is the mirror image of `Hfs`.  There a NUCLEAR moment couples to the electrons, giving the hyperfine constants
# A, B and C.  Here an EXTERNAL field of the same multipolarity pushes on the electrons, the cloud is distorted, and
# the distorted cloud makes a field of its own AT THE NUCLEUS.  The ratio is the shielding factor, which is why the
# two modules carry the same multipole selectors.
#
# WHY IT MATTERS.  For a heavy ion the electrons do not screen the applied gradient, they AMPLIFY AND INVERT it -- by
# a factor near 185 for Th(4+).  So a nuclear quadrupole splitting measured in a crystal is about 99 % an atomic
# quantity and 1 % crystallography, and the factor is the whole bridge between a lattice geometry and an observable.
# The same bridge, at rank 1 and magnetic, converts an NMR frequency into a nuclear magnetic moment -- which is how a
# wrong shielding constant once produced an apparent 7-sigma failure of strong-field QED in bismuth.
#
# WHAT IS AND IS NOT COMPUTED.  A shielding factor is a LINEAR response of the FREE ION to a field assumed UNIFORM
# over it.  This module returns exactly that, relativistically, with no crystal or molecular content: no charge
# transfer, no bonding, and (so far) no self-consistency -- the induced field is not yet allowed to act back.
# `displayResults` states all of this on every run, beside the number, rather than leaving it to a docstring.


if  true
    # Last visit:      29-Sep-2026
    # Last successful: 29-Sep-2026 -- Th(4+) [Rn] gives gamma_inf = -184.95, hence 1 - gamma_inf = 185.95.
    #   VERIFIED AGAINST AN INDEPENDENT RELATIVISTIC CALCULATION: Feiock & Johnson, Phys. Rev. 187, 39 (1969),
    #   quoted by Peik & Okhapkin, C. R. Physique 16, 516 (2015), give about -178 for Rn-like Th(4+) -- so 4 %,
    #   from a different method, different code and a different decade.  The breakdown printed below also shows
    #   the outermost p shell dominating, which is the expected structure for a closed-shell ion.
    println("\nCp-a)  Th(4+): the quadrupole antishielding factor, and where it comes from.")
    setDefaults("print summary: open", "zzz-NuclearShielding.sum")
    nm   = Nuclear.Model(90.0)
    conf = Configuration("[Rn]")
    grid = Basics.recommendedGrid([conf], nm; printout=false)
    setDefaults("standard grid", grid)
    asf  = AsfSettings(AsfSettings(); scField=Basics.DFSField(1.0))
    wa   = Atomic.Computation(Atomic.Computation(); name="Th4+", grid=grid, nuclearModel=nm,
                              configs=[conf], asfSettings=asf)
    mp   = perform(wa; output=true)["multiplet:"]
    NuclearShielding.computeOutcomes(mp, nm, grid, NuclearShielding.Settings())
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      29-Sep-2026
    # Last successful: 29-Sep-2026 -- the three CLOSED-SHELL ions of Sternheimer, Phys. Rev. 159, 266 (1967),
    #   Table I.  This work against his RAW values (his published ones carry a hand-applied factor 1/1.1 for his
    #   use of neutral-atom orbitals, which we do not need since we use ion orbitals):
    #       Y(3+)   -36.2  against  -34.1    +6 %
    #       In(3+)  -27.4  against  -25.3    +8 %
    #       Bi(3+)  -53.4  against  -42.9   +24 %
    #   THE TREND IS THE RESULT, not the individual numbers: the excess grows with Z, which is the shape of a
    #   relativistic correction -- his calculation is non-relativistic and ours is not.  Note also that his
    #   l -> l+-2 channels were ESTIMATED from a Thomas-Fermi model and are 5-10 % of his total, while we compute
    #   them properly, so agreement much better than ~15 % would not be meaningful in either direction.
    println("\nCp-b)  The validation set: three closed-shell ions with published values.")
    for (name, Z, cstr) in (("Y^3+", 39.0, "[Kr]"), ("In^3+", 49.0, "[Kr] 4d^10"),
                            ("Bi^3+", 83.0, "[Xe] 4f^14 5d^10 6s^2"))
        conf = Configuration(cstr);   nm = Nuclear.Model(Z)
        grid = Basics.recommendedGrid([conf], nm; printout=false);   setDefaults("standard grid", grid)
        asf  = AsfSettings(AsfSettings(); scField=Basics.DFSField(1.0))
        mp   = perform(Atomic.Computation(Atomic.Computation(); name=name, grid=grid, nuclearModel=nm,
                       configs=[conf], asfSettings=asf); output=true)["multiplet:"]
        out  = NuclearShielding.computeOutcomes(mp, nm, grid, NuclearShielding.Settings(); output=true)
        @printf("\n  >>> %-8s gamma_inf = %9.2f\n\n", name, out[1].gammaE2)
    end
    #
elseif  false
    # Last visit:      29-Sep-2026
    # Last successful: 29-Sep-2026 -- A NEGATIVE RESULT, AND THAT IS WHY THE BRANCH EXISTS.  The far-field
    #   coefficient assumes the external charges are infinitely distant.  In CaF2 the nearest neighbours sit
    #   2.4-3.9 A away, comparable with the ion, so the uniform-field assumption was expected to fail and to
    #   explain a factor 1.55-1.69 between the computed 186 and the 110 extracted from measurements.
    #   IT DOES NOT.  Replacing r^2 by the exact q r_<^2 / r_>^3 of a charge at the real distance changes the
    #   answer by 0.0-0.5 %, because Th(4+) is small and tightly bound and almost none of its density lies
    #   beyond 2.4 A.  So the factorisation into lattice x atomic factor is SOUND for this system, and the
    #   110-vs-186 discrepancy must be sought elsewhere -- in the assigned defect charges or in the screening
    #   this module does not yet include.  The sanity check (a charge at 1000 A) reproduces the far field exactly.
    println("\nCp-c)  Does the near neighbour break the uniform-field assumption?  (It does not.)")
    nm   = Nuclear.Model(90.0);   conf = Configuration("[Rn]")
    grid = Basics.recommendedGrid([conf], nm; printout=false);   setDefaults("standard grid", grid)
    asf  = AsfSettings(AsfSettings(); scField=Basics.DFSField(1.0))
    mp   = perform(Atomic.Computation(Atomic.Computation(); name="Th4+", grid=grid, nuclearModel=nm,
                   configs=[conf], asfSettings=asf); output=true)["multiplet:"]
    ANG  = 1.0/0.529177210903
    for model in (NuclearShielding.UniformField(),
                  NuclearShielding.NeighbourField([(1.0, 1000.0*ANG)]),
                  NuclearShielding.NeighbourField([( 2.0, 3.8627*ANG)]),
                  NuclearShielding.NeighbourField([( 1.0, 2.3654*ANG)]))
        out = NuclearShielding.computeOutcomes(mp, nm, grid,
                  NuclearShielding.Settings(NuclearShielding.Settings(); fieldModel=model); output=true)
        @printf("\n  >>> %-46s gamma = %9.2f\n\n", string(model), out[1].gammaE2)
    end
    #
elseif  false
    # Last visit:      29-Sep-2026
    # Last successful: 29-Sep-2026 -- the sum over intermediate states converges FAST and on BOUND virtuals:
    #   keeping only the largest 1, 2, 5, 10, 20 terms per channel reaches the full answer to better than 1 %
    #   by about 20 states out of ~180 available.  That is the practical argument for the sum-over-states
    #   formulation: the answer can be SEEN to come from a handful of low-lying orbitals rather than from the
    #   continuum, so the quality of the pseudo-continuum is not the limiting factor here.
    println("\nCp-d)  How many intermediate states does the response actually need?")
    nm   = Nuclear.Model(90.0);   conf = Configuration("[Rn]")
    grid = Basics.recommendedGrid([conf], nm; printout=false);   setDefaults("standard grid", grid)
    asf  = AsfSettings(AsfSettings(); scField=Basics.DFSField(1.0))
    mp   = perform(Atomic.Computation(Atomic.Computation(); name="Th4+", grid=grid, nuclearModel=nm,
                   configs=[conf], asfSettings=asf); output=true)["multiplet:"]
    for nv in (1, 2, 5, 10, 20, 0)
        out = NuclearShielding.computeOutcomes(mp, nm, grid,
                  NuclearShielding.Settings(NuclearShielding.Settings(); nVirtualMax=nv); output=true)
        @printf("\n  >>> nVirtualMax = %3d   gamma_inf = %9.2f\n\n", nv, out[1].gammaE2)
    end
    #
elseif  false
    # Last visit:      29-Sep-2026
    # Last successful: 29-Sep-2026 -- A GUARD TEST: both calls below MUST raise, and raising is the pass.
    #   (i)  an OPEN-SHELL ion.  A shielding factor is the SCALAR response of a spherical ion; for an open shell
    #        the response is not one number and the factorisation does not apply.  Returning something plausible
    #        would be worse than refusing, so the module refuses and says why.
    #   (ii) calcM1.  The magnetic-dipole shielding sigma -- the NMR analogue, and the quantity behind the
    #        bismuth hyperfine puzzle -- is NOT implemented.  Both rank-1 operators already exist in JAC
    #        (InteractionStrength.zeeman_n1 and hfs_tM1) but the response is not wired up or validated, so the
    #        settings refuse rather than return a number that would be wrong rather than approximate.
    println("\nCp-e)  Two things the module refuses to do, and says why.")
    nm   = Nuclear.Model(90.0)
    for (what, cstr, sets) in (("an OPEN-SHELL ion", "[Rn] 5f^2", NuclearShielding.Settings()),
                               ("calcM1 (not implemented)", "[Rn]",
                                NuclearShielding.Settings(NuclearShielding.Settings(); calcM1=true)))
        conf = Configuration(cstr)
        grid = Basics.recommendedGrid([conf], nm; printout=false);   setDefaults("standard grid", grid)
        asf  = AsfSettings(AsfSettings(); scField=Basics.DFSField(1.0))
        try
            mp = perform(Atomic.Computation(Atomic.Computation(); name="x", grid=grid, nuclearModel=nm,
                         configs=[conf], asfSettings=asf); output=true)["multiplet:"]
            NuclearShielding.computeOutcomes(mp, nm, grid, sets)
            @printf("\n  >>> %-28s DID NOT RAISE -- this is a FAILURE of the guard.\n\n", what)
        catch err
            @printf("\n  >>> %-28s raised as designed:\n       %s\n\n", what,
                    first(split(sprint(showerror, err), "\n"), 4) |> x -> join(x, "\n       "))
        end
    end
    #

elseif  false
    # Last visit:      04-Oct-2026
    # Last successful: 04-Oct-2026 -- THE OVERALL CONSTANT IS DERIVED, AND THIS BRANCH CHECKS EVERY LINK.
    #   `NuclearShielding.GAMMA_PREFACTOR` used to carry a factor of two and a sign that had been fixed by
    #   calibration against Sternheimer's ions, which is a poor state for the one constant that multiplies
    #   every number the module returns.  It is now derived from four factors, and each of them can be checked
    #   here WITHOUT SOLVING ANYTHING: no SCF, no orbitals, no grid.  The branch runs in under a second and it
    #   fails loudly if any link breaks, which is what makes it worth having beside the physics branches.
    println("\nCp-f)  The overall constant of the quadrupole response, link by link.")
    cpOk = Ref(true)
    cpCheck = function (name, got, want, tol)
        good = abs(got - want) <= tol
        @printf("  %-54s %14.9f   expected %12.9f   %s\n", name, got, want, good ? "ok" : "FAILED")
        cpOk[] = cpOk[] && good
    end

    # LINK 1 -- the angular reduction.  Summed over the magnetic substates of a closed subshell,
    #           sum_m |(j_a k j_v; -m_a 0 m_v)|^2 = 1/(2k+1), which is 1/5 at the quadrupole rank k = 2.
    println("\n  link 1:  the 3j orthogonality that reduces a closed subshell, at rank k = 2")
    for  (ja2, jv2)  in  ((1,3), (3,3), (3,5), (5,5), (5,7), (7,7))
        cpSum = 0.
        for  ma2 = -ja2:2:ja2
            cpSum = cpSum + AngularMomentum.Wigner_3j(AngularJ64(ja2//2), AngularJ64(2), AngularJ64(jv2//2),
                                                      AngularM64(-ma2//2), AngularM64(0), AngularM64(ma2//2))^2
        end
        cpCheck("      j_a = $(ja2)/2,  j_v = $(jv2)/2", cpSum, 1/5, 1.0e-10)
    end

    # LINK 2 -- the applied gradient.  `drivingFunction` gives the RADIAL factor of the external rank-2
    #           potential and `appliedGradient` the field gradient that same potential makes at the nucleus.
    #           They are not independent: d^2/dz^2 of r^2 P_2 is 2, so appliedGradient = 2 lim_{r->0} f(r)/r^2.
    #           This is the link that fixes the SIGN, and it needs no outside paper at all.
    println("\n  link 2:  appliedGradient = 2 x lim(r->0) drivingFunction(r)/r^2, for both field models")
    for  (cpWhat, cpModel)  in  (("a uniform field", NuclearShielding.UniformField()),
                                 ("two charges, +1 at 4.0 and -2 at 5.5 a.u.",
                                  NuclearShielding.NeighbourField([(1.0, 4.0), (-2.0, 5.5)])))
        cpR   = 1.0e-4
        cpLim = NuclearShielding.drivingFunction(cpModel, cpR) / (cpR*cpR)
        cpCheck("      $cpWhat", NuclearShielding.appliedGradient(cpModel), 2.0*cpLim, 1.0e-8)
    end

    # LINK 3 -- the angular content, against a validation set that exists independently of us.  4/5 times the
    #           sum of |<kappa_a||C^(2)||kappa'>|^2 over the spin-orbit partners of a closed shell must give
    #           Sternheimer's non-relativistic closed-shell coefficients.  Note that 4/5 is TWICE the constant,
    #           because his coefficients are quoted for a radial integral defined with the other half.
    println("\n  link 3:  4/5 x sum |<kappa||C2||kappa'>|^2 against Sternheimer's closed-shell coefficients")
    for  (cpLab, cpKappas, cpStern)  in  (("p", [1,-2], 48/25), ("d", [2,-3], 16/7), ("f", [3,-4], 224/75))
        cpSq = 0.
        for  ka  in  cpKappas,  kb  in  cpKappas
            cpSq = cpSq + AngularMomentum.CL_reduced_me(Subshell(9, ka), 2, Subshell(9, kb))^2
        end
        @printf("      the %s shell:  sum of squares = %10.7f\n", cpLab, cpSq)
        cpCheck("          4/5 x that", 0.8*cpSq, cpStern, 1.0e-9)
    end

    # LINK 4 -- the product, which is the constant itself, and its sign.
    println("\n  link 4:  the four factors multiplied out")
    println("      2 (linear response)  x  2 (the gradient operator is 2 P_2/r^3)  x  1/5 (rank 2)  /  2 (applied)")
    cpCheck("      the product", 2.0 * 2.0 * (1/5) / 2.0, 0.4, 1.0e-12)
    cpCheck("      |GAMMA_PREFACTOR|", abs(NuclearShielding.GAMMA_PREFACTOR), 0.4, 1.0e-12)
    @printf("  %-54s %14s   expected %12s   %s\n", "      its sign, from gamma = -q_induced/q_applied",
            NuclearShielding.GAMMA_PREFACTOR < 0 ? "negative" : "positive", "negative",
            NuclearShielding.GAMMA_PREFACTOR < 0 ? "ok" : "FAILED")
    cpOk[] = cpOk[] && (NuclearShielding.GAMMA_PREFACTOR < 0)

    println("\n  " * "-"^92)
    println(cpOk[] ? "  ALL LINKS HOLD:  GAMMA_PREFACTOR = -2/5 is derived, not fitted." :
                     "  >>> A LINK FAILED.  The constant is no longer consistent with its own derivation.")
    println("  " * "-"^92)
    #

elseif  false
    # Last visit:      04-Oct-2026
    # Last successful: 04-Oct-2026 -- SELF-CONSISTENT SCREENING, AND A RESULT THAT CONTRADICTS WHAT WAS EXPECTED.
    #   The module's default response is UNCOUPLED: the induced field is not allowed to act back.  With
    #   selfConsistent = true it is, and the loop converges (spread below 0.01 over the last four of 16 steps at
    #   mixing 0.5; UNDAMPED IT OSCILLATES, so do not read six steps and stop).
    #   WHAT WAS EXPECTED was that |gamma| would FALL by 10-30 % everywhere, the induced field opposing the
    #   applied one.  Measured, it falls for the two lighter ions and RISES for thorium:
    #        Y(3+)    -36.18 -> -32.07   -11.4 %
    #        In(3+)   -27.39 -> -18.45   -32.7 %   (outside the expected range)
    #        Th(4+)  -184.95 -> -194.62   +5.2 %   (the WRONG WAY)
    #   Thorium is the one that matters, since -184.95 is to be compared with the 110-120 extracted from the
    #   CaF2 measurements: self-consistency moves it AWAY from them, so it is argued against as the explanation
    #   of that gap rather than for.  The numbers above include the exchange response; direct-only gives
    #   -31.64, -18.09 and -197.99, i.e. exchange is worth about 1.5 percentage points and changes no sign.
    println("\nCp-g)  Uncoupled against self-consistent, for three closed-shell ions.")
    cpSets0 = NuclearShielding.Settings()
    cpSets1 = NuclearShielding.Settings(NuclearShielding.Settings(); selfConsistent=true, scfIterations=16, scfMixing=0.5)
    @printf("\n  %-8s %5s %14s %16s %10s\n", "ion", "Z", "uncoupled", "self-consistent", "change")
    println("  " * "-"^62)
    for  (cpName, cpZ, cpConf)  in  (("Y^3+", 39.0, "[Kr]"), ("In^3+", 49.0, "[Kr] 4d^10"), ("Th^4+", 90.0, "[Rn]"))
        cpC  = Configuration(cpConf);   cpNm = Nuclear.Model(cpZ)
        cpG  = Basics.recommendedGrid([cpC], cpNm; printout=false);   setDefaults("standard grid", cpG)
        cpA  = AsfSettings(AsfSettings(); scField=Basics.DFSField(1.0))
        cpT  = tempname()
        cpMp = open(cpT,"w") do io;  redirect_stdout(io) do
                   perform(Atomic.Computation(Atomic.Computation(); name=cpName, grid=cpG, nuclearModel=cpNm,
                           configs=[cpC], asfSettings=cpA); output=true)["multiplet:"]  end  end
        cpV = Float64[]
        for  cpS  in  (cpSets0, cpSets1)
            cpO = open(cpT,"w") do io;  redirect_stdout(io) do
                      NuclearShielding.computeOutcomes(cpMp, cpNm, cpG, cpS; output=true)  end  end
            push!(cpV, cpO[1].gammaE2)
        end
        rm(cpT, force=true)
        @printf("  %-8s %5.0f %14.2f %16.2f %9.1f %%\n", cpName, cpZ, cpV[1], cpV[2],
                100*(abs(cpV[2])-abs(cpV[1]))/abs(cpV[1]))
        flush(stdout)
    end
    println("  " * "-"^62)
    println("  The rise at Z = 90 is the finding, and it is why priority item 43 was rewritten rather than closed.")
    #

elseif  false
    # Last visit:      04-Oct-2026
    # Last successful: 04-Oct-2026 -- THE GRID, AND THE ONE THING THAT DECIDES WHETHER A SHIELDING FACTOR IS REAL.
    #   gamma_inf is NOT box-dependent, which is the opposite of what was believed for a day.  What decides it is
    #   the SPLINE DENSITY, nsL / r_box.  Hold that fixed and the box may be enlarged tenfold with no change;
    #   let it fall and the number goes wrong smoothly, through the sign and beyond.
    #   The cause is one line in `Basics.recommendedGrid`: hp = r_box/300 keeps the number of outer points
    #   roughly FIXED however large the box, so the B-spline basis does not grow with the box.  Pass hp
    #   explicitly and the problem disappears.
    #   MEASURED HERE, on Y(3+) [Kr]:  r_box 8.7 at 10.7 splines/a.u. gives -36.18, and r_box 92 at 5.6 gives
    #   -36.18 as well -- while r_box 92 on the DEFAULT hp, at 1.07 splines/a.u., gives -41.98.
    #   AND PASSING THE FLOOR PROVES NOTHING, because the density an ion NEEDS belongs to the ion.  Refined by
    #   halving hp four times, the four reference ions converge to -36.14, -102.86, -189.72 and -62.00, and their
    #   own recommended grids are 0.1 %, 1.0 %, 2.5 % and 190 % away from those -- the last being Hg(2+), whose
    #   4f^14 5d^10 shells need 27 splines/a.u. where Y(3+) needs 11.  So the real test is to HALVE hp and
    #   recompute; `tools/probe-nuclearShieldingDensity.jl` does that, and the cancellation ratio A1 is what
    #   flags the Hg(2+) default grid (1.82 against 0.33-0.44 for the three sound ones).
    println("\nCp-h)  What the spline density does to gamma_inf, and what the two guards do and do not catch.")
    cpConf = Configuration("[Kr]");    cpNm = Nuclear.Model(39.0)
    cpAsf  = AsfSettings(AsfSettings(); scField=Basics.DFSField(1.0))
    @printf("\n  %-28s %8s %6s %10s %12s %7s\n", "grid", "r_box", "nsL", "per a.u.", "gamma_inf", "A1")
    println("  " * "-"^80)
    cpRows = Tuple{String,Float64,Union{Nothing,Float64}}[("recommended, default hp", 0.0, nothing),
                   ("tenfold box, hp HELD FIXED", 92.0, 0.0287), ("tenfold box, default hp", 92.0, nothing)]
    for  (cpWhat, cpRbox, cpHp)  in  cpRows
        cpG = cpRbox <= 0.      ?  Basics.recommendedGrid([cpConf], cpNm; printout=false)             :
              isnothing(cpHp)   ?  Basics.recommendedGrid([cpConf], cpNm; rbox=cpRbox, printout=false) :
                                   Basics.recommendedGrid([cpConf], cpNm; rbox=cpRbox, hp=cpHp, printout=false)
        setDefaults("standard grid", cpG)
        cpT = tempname()
        cpO = open(cpT,"w") do io;  redirect_stdout(io) do
                  cpMp = perform(Atomic.Computation(Atomic.Computation(); name="Y^3+", grid=cpG, nuclearModel=cpNm,
                                 configs=[cpConf], asfSettings=cpAsf); output=true)["multiplet:"]
                  NuclearShielding.computeOutcomes(cpMp, cpNm, cpG, NuclearShielding.Settings(); output=true)[1]
              end  end
        rm(cpT, force=true)
        cpA1 = maximum(abs(c.value) for c in cpO.contributions) / abs(cpO.gammaE2)
        @printf("  %-28s %8.1f %6d %10.2f %12.2f %7.2f\n", cpWhat, cpG.tL[end], cpG.nsL,
                cpG.nsL/cpG.tL[end], cpO.gammaE2, cpA1)
        flush(stdout)
    end
    println("  " * "-"^80)
    println("  The first two agree to five figures on boxes a factor of ten apart; the third, whose basis was not")
    println("  allowed to grow with its box, is 16 % wrong -- and its cancellation ratio is as healthy as theirs,")
    println("  which is why the GRID check and the CANCELLATION check are both needed and still not sufficient.")
    #
end
