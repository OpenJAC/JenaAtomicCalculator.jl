#
# Validation of StarkShift.magicWavelengths, added 20-Oct-2026.
#
# A magic wavelength of a PAIR of levels is a zero of alpha_0(omega, A) - alpha_0(omega, B):  there the trapping
# light shifts both levels equally, so the TRANSITION between them is untouched to first order in the intensity,
# which is what lets an optical trap coexist with a clock or a qubit.  There are generically MANY of them, one or
# more between every pair of consecutive resonances, so the function returns a LIST and never "the" magic
# wavelength.
#
# THE WHOLE DIFFICULTY IS THAT A SIGN CHANGE IS NOT A ZERO.  alpha_0(omega) runs to +infinity as omega approaches
# an excitation energy from below and returns from -infinity above it, so the difference jumps sign across EVERY
# resonance of EITHER level.  A bisection let loose on such a sign change converges contentedly onto the POLE and
# reports a resonance as a magic wavelength -- and that error looks entirely plausible, because what comes back is
# a sharp, reproducible frequency.
#
# THE FIX IS IN TWO PARTS, AND THE SECOND WAS FOUND BY THIS PROBE RATHER THAN FORESEEN.
#
#   PART ONE, the pole partition.  The poles are known EXACTLY -- they are the excitation energies that
#   StarkShift.polarizabilityTerms returns -- so the range is cut at every one of them and each bracket is taken
#   strictly inside a pole-free subinterval.  Each subinterval is SAMPLED rather than merely tested at its ends,
#   because between a pole of A and a pole of B the difference runs from -infinity to -infinity and therefore holds
#   an EVEN number of zeros, which an end-point test misses entirely.
#
#   PART TWO, the residual check, WITHOUT WHICH PART ONE IS NOT ENOUGH.  The first version did the partition alone
#   and returned TWO POLES as magic wavelengths.  THE POLES COME IN NEAR-DEGENERATE PAIRS:  in hydrogen 2p -> 3s
#   and 2p -> 3d sit about 1e-7 apart, and the fine-structure partners of one subshell closer still.  The sliver
#   between such a pair is below numerical resolution, so the sign change inside it is the PAIR OF POLES and not a
#   zero.  Every candidate is therefore verified by |alpha_A - alpha_B| against the scale of the two
#   polarizabilities, and the two populations turn out to be seven orders apart.
#
# CASE A  H 1s against H 2p_1/2 over omega = 0.01 ... 0.25 a.u., with three checks:
#           1. no returned root may sit near a known pole;
#           2. COMPLETENESS against a dense independent scan of 400 000 points, in which every sign change is
#              classified by whether a known pole lies between the two samples -- the count of pole-free sign
#              changes must equal the number of roots returned, with nothing extra;
#           3. the residual |alpha_A - alpha_B| at each root.
#
# RESULTS, 20-Oct-2026 (this file is the validation record):
#
#   alpha_0(1s, 0) = 3.49182825 a.u.,  alpha_0(2p_1/2, 0) = 160.30087308 a.u.;  ten poles in range, in five
#   near-degenerate PAIRS:  0.069446, 0.069446, 0.093752, 0.093752, 0.375005, 0.375006, 0.444450, 0.444451,
#   0.468756, 0.468757.
#
#   BEFORE THE RESIDUAL CHECK -- three "roots" returned, two of them poles:
#        omega = 0.069445794   |alpha_A - alpha_B|/alpha = 1.23e+01   slope -1.20e+12   <- a POLE
#        omega = 0.090592606   |alpha_A - alpha_B|/alpha = 3.45e-09   slope -7.79e+04   <- the zero
#        omega = 0.093751747   |alpha_A - alpha_B|/alpha = 1.06e+01   slope -1.82e+13   <- a POLE
#   The two spurious ones sat 5.4e-08 (relative) from a pole, i.e. essentially ON it.
#
#   AFTER -- ONE root, and it is the right one:
#        omega = 0.090592606 a.u. = 502.948 nm,  alpha_A = alpha_B = 3.697721 a.u. to 3.45e-09,
#        slope -7.79e+04 a.u./a.u., nearest pole 3.37e-02 away (relative), flagged aboveResonance since it lies
#        above 2p's lowest excitation energy of 0.069446 a.u.
#   The dense scan finds exactly 1 pole-free sign change and 2 that straddle a pole;  the function returns 1.
#   COMPLETE, and nothing spurious.
#
# WHAT THIS CASE DOES NOT EXERCISE, said plainly:  every zero here lies ABOVE the lowest resonance of one of the
# two levels, because alpha_0(1s) = 3.49 and alpha_0(2p) = 160.3 do not cross below it -- the difference is
# negative at omega -> 0 and runs to -infinity at 2p's first pole.  A magic wavelength BELOW all resonances, which
# is the case a clock actually wants, needs a pair whose static polarizabilities are closer together;  that belongs
# in an example file with a real alkali or alkaline-earth pair, not here.  The ALGORITHM is what this probe
# validates, and completeness against a dense scan does that independently of which side of a resonance the zeros
# fall on.
#
# Run as:  julia --project=. tools/probe-magicWavelength.jl
#

using JenaAtomicCalculator
const JAC = JenaAtomicCalculator
using Printf

const SS = StarkShift


# the two levels and the one perturber multiplet that serves both:  each level's own selection rules decide which
# members it reaches, so a set of BOTH parities is right here and unreachable members drop out on their own
function setup()
    nm    = Nuclear.Model(1., PointNucleus())
    confs = [Configuration("1s"), Configuration("2p"), Configuration("3p"), Configuration("4p"),
             Configuration("3s"), Configuration("4s"), Configuration("3d"), Configuration("4d")]
    grid  = Basics.recommendedGrid(confs, nm)
    one   = function(name::String, cs::Array{Configuration,1})
        redirect_stdout(devnull) do
            perform( Atomic.Computation(Atomic.Computation(); name = name, grid = grid, nuclearModel = nm,
                                        configs = cs); output = true )[ResultKeys.Multiplet]
        end
    end
    gm   = one("all perturbers", confs)
    levA = one("1s", [Configuration("1s")]).levels[1]
    levB = one("2p", [Configuration("2p")]).levels[1]

    return( (levA, levB, gm, grid) )
end


function caseA(levA::Level, levB::Level, gm::Multiplet, grid::Radial.Grid)
    sA, eA = SS.polarizabilityTerms(levA, gm, grid)
    sB, eB = SS.polarizabilityTerms(levB, gm, grid)
    dAlpha = om -> SS.rebuiltPolarizability(sA, eA, om) - SS.rebuiltPolarizability(sB, eB, om)
    poles  = sort( unique( [dE  for dE in vcat(eA, eB)  if dE > 0.0] ) )
    @printf("\n   alpha_0(A, 0) = %.8f a.u.      alpha_0(B, 0) = %.8f a.u.\n",
            SS.rebuiltPolarizability(sA, eA, 0.0), SS.rebuiltPolarizability(sB, eB, 0.0))
    println("   poles in and beyond the range (a.u.):  ", [round(p, digits=6) for p in poles])

    omegaMin = 0.01;     omegaMax = 0.25
    mws = SS.magicWavelengths(levA, levB, gm, grid; omegaMin = omegaMin, omegaMax = omegaMax)

    println("\n   ", "-"^112)
    @printf("   %-3s %-15s %-13s %-17s %-17s %-12s %-13s %-7s\n",
            "#", "omega [a.u.]", "lambda [nm]", "alpha_A", "alpha_B", "residual", "slope", "above?")
    for  (i, mw) in enumerate(mws)
        aA = SS.rebuiltPolarizability(sA, eA, mw.omega);     aB = SS.rebuiltPolarizability(sB, eB, mw.omega)
        @printf("   %-3d %-15.9f %-13.3f %-17.9f %-17.9f %-12.2e %-13.4e %-7s\n",
                i, mw.omega, mw.wavelength, aA, aB, abs(aA-aB)/max(abs(aA), 1.0e-30), mw.slope, mw.aboveResonance)
    end

    # 1. how close does any returned root come to a pole?
    worst = 1.0e30
    for  mw in mws, p in poles      worst = min(worst, abs(mw.omega - p)/p)      end
    @printf("\n   CHECK 1  the closest any returned root comes to a known pole:  %.3e (relative).  A pole taken\n",
            worst)
    println("            for a root shows up here at ~1e-06, the inset used when bracketing.")

    # 2. completeness against a dense scan, classifying each sign change by whether a pole lies between the samples
    nFine = 400000;      h = (omegaMax - omegaMin)/nFine
    nRoot = 0;      nPole = 0;      fPrev = dAlpha(omegaMin);      xPrev = omegaMin
    for  k = 1:nFine
        x = omegaMin + k*h;      f = dAlpha(x)
        if  fPrev*f < 0.0
            any(p -> xPrev < p < x, poles) ? (nPole += 1) : (nRoot += 1)
        end
        fPrev = f;      xPrev = x
    end
    @printf("\n   CHECK 2  a dense scan of %d points finds %d sign change(s) straddling NO pole (genuine\n",
            nFine, nRoot)
    @printf("            zeros) and %d straddling a pole (which must NOT be reported).  The function returned %d.\n",
            nPole, length(mws))
    println(length(mws) == nRoot ? "            COMPLETE, and nothing spurious." :
                                   "            ** MISMATCH -- read the table above **")

    println("\n   CHECK 3  the printout a user sees:")
    for  mw in mws      print("   ");  show(mw)      end

    return( nothing )
end


function main()
    println("\n", "="^118)
    println("THE MAGIC-WAVELENGTH ROOT FIND -- validation record, 20-Oct-2026")
    println("="^118)
    levA, levB, gm, grid = setup()
    @printf("\n   level A = %s at %.6f a.u.      level B = %s at %.6f a.u.\n",
            string(LevelSymmetry(levA.J, levA.parity)), levA.energy,
            string(LevelSymmetry(levB.J, levB.parity)), levB.energy)
    caseA(levA, levB, gm, grid)
    println("\n", "="^118)

    return( nothing )
end

main()
