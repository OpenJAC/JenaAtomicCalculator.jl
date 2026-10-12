#
# probe-cfSecondOrder.jl
#
# A PROBE, not a module:  how does a SECOND-ORDER treatment of a large correlation layer compare with the
# variational CI over the same layer -- in energy, in time and in memory?
#
# Test case:  Cf^17+  [Xe] 4f^14 5d^10 6s^2 5f / 6p,  the 5f -> 6p clock interval, with the n=7 correlation
#             layer opened to 5d CORE excitations.  P is the space a variational calculation can afford
#             (valence singles and doubles plus 5d core SINGLES);  Q is the 5d core DOUBLES.
#
# Usage:  julia --project=. tools/probe-cfSecondOrder.jl [both|variational|pt] [2J]
#
# WHAT IT MEASURED, 12-Oct-2026 (priority item 53), and the three results worth keeping:
#
#   1. THE VARIATIONAL ROUTE IS AFFORDABLE BY THIS ROUTE AND NOT BY THE OBVIOUS ONE.  Driving the same
#      configuration list through Atomic.Computation/perform was killed by the OS at 14.5 GB after 3:42, in
#      CONFIGURATION SETUP, before any CI matrix existed.  Building the basis with Basics.generateBasis and
#      assembling one symmetry block directly costs 3.21 GB and ten minutes for 11 865 CSF.  That gap is a JAC
#      question in its own right.
#
#   2. SECOND ORDER IS GOOD PER LEVEL AND THREE TIMES WORSE ON THE OBSERVABLE.  It recovered 103.9 % and
#      103.5 % of the doubles for the two clock levels, but the levels' corrections very nearly cancel in the
#      interval -- the shift is 4.6 % of the per-level correction -- so the error on the interval is 10.5 %.
#      **A relative error that survives a 95 % cancellation is amplified by the reciprocal of what survives.**
#      For a transition, an isotope shift or a q-factor, validate the DIFFERENCE and not the levels.
#
#   3. NO INTRUDERS, which is why it works at all:  min|D| = 11.5 Ha over 11 442 perturbers, none below
#      0.1 Ha, because a 5d core excitation costs 11-24 Ha.  The same machinery on a VALENCE-valence space
#      is wrong by 3.6-11.5 % (see probe-secondOrderQ.jl), where the coupling is strong.
#
# The orbital construction is NOT optional:  occupied and virtual alike are eigenstates of ONE local potential,
# the reference configuration's own DFS field.  Two earlier constructions in probe-secondOrderQ.jl put ~99 % of
# the second-order ranking into a single excitation that was orbital relaxation masquerading as correlation.
#
# The correction is STATE-SPECIFIC -- one shift per level, added separately, never a diagonalised effective
# Hamiltonian:  that was measured to make splittings WORSE (1021 cm^-1 error against 701), because each P-CSF's
# denominator is level-independent.
#
using JenaAtomicCalculator, LinearAlgebra, Printf
const B  = Basics
const SC = JenaAtomicCalculator.SelfConsistent
const H  = JenaAtomicCalculator.Hamiltonian
# Item 53, step 2:  BOTH LINES on the same block, from the same orbitals.
#   (a) VARIATIONAL -- the full CI over P (+) Q
#   (b) SECOND ORDER -- CI over P alone, plus a state-specific Epstein-Nesbet shift from Q
# P is round 2's SHIPPED space (valence SD + 5d core singles);  Q is the 5d core DOUBLES it could not afford.
# The lean route is used deliberately:  the Atomic.Computation path died at 14.5 GB in SETUP on this space,
# while Basics.generateBasis builds the same basis comfortably, so the matrix is assembled here directly.
mode  = length(ARGS) > 0 ? ARGS[1] : "both"
twoJ  = length(ARGS) > 1 ? parse(Int, ARGS[2]) : 5
sym   = LevelSymmetry(AngularJ64(twoJ, 2), Basics.minus)
Zcf   = 98.0;  core = "[Xe] 4f^14 5d^10"
refs  = [Configuration("$core 6s^2 5f^1"), Configuration("$core 6s^2 6p^1")]
into  = [Shell("7s"), Shell("7p"), Shell("6d"), Shell("6f"), Shell("5g")]
val   = [Shell("5f"), Shell("6s"), Shell("6p")];    cv = vcat([Shell("5d")], val)
nm    = Nuclear.Model(Zcf)
grid  = Basics.recommendedGrid(refs, nm; rnt = 2.0e-7)
set   = AsfSettings(AsfSettings(); eeInteraction=CoulombInteraction(), eeInteractionCI=CoulombInteraction(),
                    gridStopper=false)
@printf(">>> 2J = %d, mode = %s\n", twoJ, mode);  flush(stdout)
# ONE LOCAL POTENTIAL for occupied and virtual alike -- the probe-secondOrderQ construction, and it is not
# optional: two earlier attempts there put ~99 % of the second-order ranking into a single excitation that was
# orbital relaxation masquerading as correlation.  Here the reference's own DFS potential supplies every orbital.
tRef  = @elapsed mpRef = SC.performSCF(refs, nm, grid, AsfSettings(set; scField=Basics.ALField(),
                                        scfRoute=Basics.AverageLevelRoute(60)); printout=false)
refBas = mpRef.levels[1].basis
prim   = Bsplines.generatePrimitives(grid)
pot    = Basics.add( Nuclear.nuclearPotential(nm, grid), Basics.computePotential(Basics.DFSField(), grid, refBas) )
@printf(">>> reference AL + its local potential in %.1f s\n", tRef);  flush(stdout)
# the two spaces
pConfs = unique(vcat(refs, B.generateConfigurations(B.ExciteElectrons(1, val, into), refs),
                           B.generateConfigurations(B.ExciteElectrons(2, val, into), refs),
                           B.generateConfigurations(B.ExciteElectrons(1, cv,  into), refs)))
allConfs = unique(vcat(pConfs, B.generateConfigurations(B.ExciteElectrons(2, cv, into), refs)))
tBas = @elapsed basis0 = redirect_stdout(devnull) do;  B.generateBasis(allConfs, [sym])  end
orbs = redirect_stdout(devnull) do;  Bsplines.generateOrbitals(basis0.subshells, pot, nm, prim; printout=false) end
basis = ManyElectron.Basis(true, basis0.NoElectrons, basis0.subshells, basis0.csfs, basis0.coreSubshells, orbs)
n = length(basis.csfs)
@printf(">>> basis: %d CSF over %d subshells, built in %.1f s\n", n, length(basis.subshells), tBas);  flush(stdout)
# WHICH CSF belong to P?  A CSF is in P when its non-relativistic configuration is one of pConfs.
pKeys = Set{String}()
nlKey(csf, subs) = begin
    d = Dict{Shell,Int64}()
    for k = 1:length(subs)
        if csf.occupation[k] > 0
            sh = Shell(subs[k].n, Basics.subshell_l(subs[k]));   d[sh] = get(d, sh, 0) + csf.occupation[k]
        end
    end
    join(sort([ "$(sh)^$(w)" for (sh,w) in d ]), " ")
end
for  c  in  pConfs    push!(pKeys, join(sort([ "$(sh)^$(w)" for (sh,w) in c.shells if w > 0 ]), " "))    end
keysAll = [ nlKey(csf, basis.subshells) for csf in basis.csfs ]
P = findall(k -> k in pKeys, keysAll);    Q = findall(k -> !(k in pKeys), keysAll)
@printf(">>> P = %d CSF (round 2's space),  Q = %d CSF (the 5d core doubles)\n", length(P), length(Q));  flush(stdout)
cache = InteractionStrength.XLCache()
tMat  = @elapsed Hup = H.setupMatrixKinkAware(sym, basis, nm, grid, set, cache; printout=false)
@printf(">>> matrix (%d x %d) assembled in %.1f s;  %.2f GB\n", n, n, tMat, n^2*8/2^30);  flush(stdout)
Hs = Symmetric(Hup)
if  mode in ["both", "variational"]
    tD = @elapsed eVar = eigvals(Hs)
    @printf("\n>>> (a) VARIATIONAL  E = %.10f Ha   diagonalised in %.1f s\n", eVar[1], tD);  flush(stdout)
end
if  mode in ["both", "pt"]
    tP = @elapsed begin
        FP = eigen(Symmetric(Matrix(Hs[P,P])));    EP = FP.values;   CP = FP.vectors
        Hqq = [ Hup[min(q,q),max(q,q)] for q in Q ]
        Vpq = Matrix(Hs[P,Q])
        c   = CP[:,1];    Viq = vec(c' * Vpq);    D = EP[1] .- Hqq
        dE2 = sum(Viq.^2 ./ D);   c2 = sum((Viq ./ D).^2)
    end
    @printf(">>> (b) P-only      E = %.10f Ha\n", EP[1])
    @printf(">>> (b) SECOND ORDER E = %.10f Ha   dE(2) = %.6e Ha   sum|c|^2 = %.4e   in %.1f s\n",
            EP[1]+dE2, dE2, c2, tP)
    if  mode == "both"
        @printf("\n>>> REMAINDER (exact - 2nd order) = %.6e Ha   i.e. 2nd order recovers %.2f %% of the doubles\n",
                eVar[1] - (EP[1]+dE2), 100*dE2/(eVar[1]-EP[1]))
    end
    @printf(">>> denominators: min|D| = %.4e  median = %.4e  max = %.4e Ha;  %d of %d with |D| < 0.1\n",
            minimum(abs.(D)), sort(abs.(D))[cld(length(D),2)], maximum(abs.(D)), count(abs.(D) .< 0.1), length(D))
end
