#
# Derivation and validation of the ODD-RANK multipole polarizability added to
# module-MultipolePolarizibility.jl on 20-Oct-2026, of which K = 1 -- the VECTOR polarizability alpha_1 of the
# Stark-shift and optical-tweezer literature -- is the member anyone asks for.
#
# WHY IT WAS MISSING, which is the whole point and is not an oversight.  The module's rank coefficient carries
# the Wigner 3-j symbol (L L K; 0 0 0), which vanishes identically whenever L + L + K is odd.  That is not a
# convention to be relaxed: it says that the SYMMETRIC product of two identical multipole operators has no odd
# rank.  A field of helicity q makes the product UNSYMMETRIC, because the two second-order terms carry the
# polarization vectors eps and eps*, hence the operator components q and -q:
#
#     alpha(J,M; q)  =  sum_J'  |<J'||O||J>|^2/(2J'+1) * [ S_(+q)/(dE - omega)  +  S_(-q)/(dE + omega) ]
#
# with S_(q) = C^(J' M+q)_(J M, L q)^2.  Splitting the bracket into the part SYMMETRIC in q and the part
# ANTISYMMETRIC in q separates the two frequency factors exactly:
#
#     symmetric  (even K)  ->  dE / (dE^2 - omega^2)        MultipolePolarizibility.frequencyFactor
#     antisymmetric (odd)  ->  omega / (dE^2 - omega^2)     MultipolePolarizibility.frequencyFactorOdd
#
# SO AN ODD RANK IS PURELY DYNAMIC AND VANISHES AT omega = 0.  A static field has no handedness, and that is
# why a static-first module could be complete and still have no rank 1 in it.
#
# The coefficient is then the one that was already there with the helicity made explicit,
#
#     c(L,K,J,J')  =  2 (-1)^(J+J'+q) C^(JJ)_(JJ,K0) (2K+1)/sqrt(2J+1) * (L L K; q -q 0) * {J J' L; L K J}
#
# and q = 0 reproduces the old expression IDENTICALLY, which is checked below rather than argued.
#
# CASE A  COMPLETENESS.  alpha(J,M) is built DIRECTLY from Clebsch-Gordan coefficients, with no rank
#         decomposition anywhere in it, and then fitted by sum_K A_K(J,M) alpha^(K) over all 2J+1 substates --
#         with J > L throughout, so the fit has spare degrees of freedom and the residual means something.
#         The even rank set must FAIL, and only for a circular field at a non-zero frequency.
# CASE B  THE CLOSED FORM, one perturber at a time, so that alpha^(K) = c * |d|^2 * f_K exactly.
# CASE C  THE FREQUENCY FACTOR.  alpha^(E1,1) divided by each of the two candidate factors, over a sequence
#         of omega.  One of them must be constant at 1.
# CASE D  OUT OF SAMPLE.  Random multi-perturber sets, both helicities and linear, every substate.
# CASE E  THE MODULE ITSELF: bit-identity of the default path, and the substate identity rebuilt from the
#         module's own rankCoefficient / rankShape / frequencyFactor / frequencyFactorOdd.
# CASE F  THE DEGENERACY ANCHOR -- the one PHYSICAL check here, and it needs no literature value.  For an
#         ns_1/2 level the only E1 perturbers are np_1/2 and np_3/2 with |<d>|^2 in the ratio 1 : 2.  Were the
#         two DEGENERATE, the rank-1 sum would have to vanish identically, because the vector light shift of
#         an alkali ground state is carried entirely by the fine-structure splitting.  Then the full property
#         path on the approved H 1s case, where that splitting is tiny and alpha^(E1,1) must be tiny with it.
#
# RESULTS, 20-Oct-2026 (this file is the validation record):
#
#   CASE A   with the EVEN ranks only, the residual is 1e-16 for a linear field at any omega AND for a circular
#            field at omega = 0, and 1.2e-2 ... 1.2e-1 for a circular field at omega = 0.17.  With every rank it
#            is 1.4e-16 ... 6.3e-16 throughout.  Both selection rules therefore come out of the measurement:
#            odd K needs a circular field, and it needs a frequency.
#   CASE B   worst relative deviation of the extracted coefficient from the closed form, over L = 1..3,
#            2J = 1..8, every J' in the triangle, every K, both signs of dE:  3.2e-13 (at L = 3, K = 6, where
#            the least-squares extraction of a rank-6 component is what limits it, not the formula).
#   CASE C   alpha^(E1,1) / [omega/(dE^2-omega^2)] = 1.0000000000 at omega = 0.05, 0.1, 0.2, 0.4, while
#            alpha^(E1,1) / [dE/(dE^2-omega^2)] = 0.0569, 0.1135, 0.2245, 0.4305 -- i.e. nothing.  And
#            alpha^(E1,1)(omega = 0) = -1.5e-17.
#   CASE D   worst relative deviation over 40 random sets, 3 helicities, every substate:  8.9e-15.
#   CASE E   BIT-IDENTITY: 0 of 3888 coefficients differ in any bit at the default helicity = 0, so no number
#            this module has ever produced can move and no approved reference is touched.  The substate identity
#            rebuilt from the module's own functions:  8.8e-15 over 60 random sets and 3 helicities.  And the
#            test suite is 100/100 with no approved reference re-approved.
#            allowedRanks(E1, J=1/2) gives [0] at helicity 0 and [0, 1] at helicity +-1 -- an alkali ground
#            level does have a vector polarizability, which is the case the tweezer literature lives on.
#   CASE F   c(E1,K=1,1/2,1/2) = -1/3 exactly and c(E1,K=1,1/2,3/2) = +1/6 exactly, so the statistical sum
#            1*(-1/3) + 2*(+1/6) = 0.000e+00 -- the cancellation is EXACT, in rational numbers.  The rank-0
#            sum is +1.000000000000000 and does not cancel, as it must not.  On the approved H 1s case:
#            alpha^(E1,0) = 3.491828793108 at BOTH helicities (the approved 3.491829), alpha^(E1,1)(0) = 0.0
#            exactly, and alpha^(E1,1)(0.05) = -7.558e-07 against a 2p fine-structure splitting of 1.664e-06
#            a.u. -- suppressed by just that ratio, which is the physics.
#
# THE POLARIZATION FACTOR that composes a shift from these, MultipolePolarizibility.polarizationFactor:
#            P(E1,0,q) = 1 at every helicity;  P(E1,1,q) = q;  P(E1,2,0) = 1 against P(E1,2,+-1) = -1/2, which
#            is the textbook statement that circular polarization gives MINUS ONE HALF of the linear tensor
#            shift -- the one external check in this file that fell out rather than having to be looked up.
#            For E2:  P(E2,2,+-1) = +1/2 and P(E2,4,+-1) = -2/3.
#
# Run as:  julia --project=. tools/probe-rank1Polarizability.jl           (all cases;  "ph" for the H 1s case only)
#

using JenaAtomicCalculator
const JAC = JenaAtomicCalculator
using Printf, LinearAlgebra, Random

const AM = AngularMomentum
const MP = MultipolePolarizibility

jj(twoJ::Int64) = AngularJ64(twoJ//2)
mm(twoM::Int64) = AngularM64(Rational(twoM, 2))
Ms(Jd::Float64) = [-Jd + i  for i = 0:round(Int64, 2*Jd)]

# A_K(J,M), written out here so that the decomposition does not borrow the module's own rankShape
function shapeA(Jd::Float64, K::Int64, Md::Float64)
    K == 0  &&  return( 1.0 )
    K > round(Int64, 2*Jd)  &&  return( 0.0 )

    return( AM.ClebschGordan(Jd, Md, 1.0*K, 0.0, Jd, Md) / AM.ClebschGordan(Jd, Jd, 1.0*K, 0.0, Jd, Jd) )
end


# the substate polarizability of |J M> in a field of helicity q, built DIRECTLY from Clebsch-Gordan
# coefficients: the eps term carries q and the eps* term carries -q, with the denominators dE -+ omega
function substateAlpha(L::Int64, Jd::Float64, Md::Float64, q::Int64, omega::Float64, perts)
    wa = 0.0
    for (Jnud, dE, d2) in perts
        sp = abs(Md + q) > Jnud + 1.0e-10 ? 0.0 : AM.ClebschGordan(Jd, Md, 1.0*L,  1.0*q, Jnud, Md + q)^2
        sm = abs(Md - q) > Jnud + 1.0e-10 ? 0.0 : AM.ClebschGordan(Jd, Md, 1.0*L, -1.0*q, Jnud, Md - q)^2
        wa = wa + d2 / (2*Jnud + 1.0) * ( sp/(dE - omega) + sm/(dE + omega) )
    end

    return( wa )
end


# the candidate closed form, written out independently of the module
function coeffQ(L::Int64, K::Int64, Jd::Float64, Jnud::Float64, q::Int64)
    if  Jnud < abs(Jd - L) - 1.0e-10  ||  Jnud > Jd + L + 1.0e-10   return( 0.0 )   end
    if  K > 2*L  ||  K > round(Int64, 2*Jd)                         return( 0.0 )   end
    gK  = AM.ClebschGordan(Jd, Jd, 1.0*K, 0.0, Jd, Jd)
    c3j = AM.Wigner_3j(1.0*L, 1.0*L, 1.0*K, 1.0*q, -1.0*q, 0.0)
    c6j = AM.Wigner_6j(Jd, Jnud, 1.0*L, 1.0*L, 1.0*K, Jd)
    ph  = (-1.0)^round(Int64, Jd + Jnud) * (-1.0)^abs(q)

    return( 2.0 * ph * gK * (2*K + 1.0) / sqrt(2*Jd + 1.0) * c3j * c6j )
end


# the OLD coefficient, exactly as it stood before 20-Oct-2026, for the bit-identity check
function oldCoeff(L::Int64, K::Int64, J::AngularJ64, Jnu::AngularJ64)
    Jd = Basics.twice(J) / 2.0;      Jnud = Basics.twice(Jnu) / 2.0
    if  Jnud < abs(Jd - L) - 1.0e-10   ||   Jnud > Jd + L + 1.0e-10    return( 0.0 )    end
    if  isodd(K)   ||   K > 2*L   ||   K > round(Int64, 2*Jd)          return( 0.0 )    end
    gK  = AM.ClebschGordan(Jd, Jd, 1.0*K, 0.0, Jd, Jd)
    c3j = AM.Wigner_3j(1.0*L, 1.0*L, 1.0*K, 0.0, 0.0, 0.0)
    c6j = AM.Wigner_6j(Jd, Jnud, 1.0*L, 1.0*L, 1.0*K, Jd)
    ph  = (-1.0)^round(Int64, Jd + Jnud)

    return( 2.0 * ph * gK * (2*K + 1.0) / sqrt(2*Jd + 1.0) * c3j * c6j )
end


freqFactor(K::Int64, dE::Float64, omega::Float64) =
    iseven(K) ? dE/(dE*dE - omega*omega) : omega/(dE*dE - omega*omega)


# the alpha^(K) that best reproduce the directly-built substate shifts;  (coefficients, relative residual)
function extract(L::Int64, Jd::Float64, q::Int64, omega::Float64, perts, Ks)
    Mlist = Ms(Jd)
    Amat  = [ shapeA(Jd, K, M)  for M in Mlist, K in Ks ]
    targ  = [ substateAlpha(L, Jd, M, q, omega, perts)  for M in Mlist ]
    sol   = Amat \ targ

    return( (sol, norm(Amat*sol - targ) / max(norm(targ), 1.0e-300)) )
end


function randPerts(rng, L::Int64, Jd::Float64, n::Int64)
    out  = Tuple{Float64,Float64,Float64}[]
    cand = [ Jd + d  for d = -L:L  if Jd + d >= 0.0 ]
    for i = 1:n
        push!( out, (cand[rand(rng, 1:length(cand))], 0.5 + 3.0*rand(rng), 0.2 + rand(rng)) )
    end

    return( out )
end


function caseA(rng)
    println("\nCASE A   THE EVEN RANK SET IS NOT COMPLETE FOR A CIRCULAR FIELD AT A NON-ZERO FREQUENCY.")
    @printf("\n   %-4s %-5s %-4s %-8s   %-14s   %-14s\n", "L", "2J", "q", "omega", "even K only", "all K")
    for (L, Jd) in [(1, 2.5), (1, 3.5), (2, 4.0), (2, 5.0), (3, 5.0)]
        perts = randPerts(rng, L, Jd, 6)
        for q in [0, 1], omega in [0.0, 0.17]
            Kall  = [K for K = 0:min(2L, round(Int64, 2*Jd))]
            Keven = [K for K in Kall if iseven(K)]
            rE = extract(L, Jd, q, omega, perts, Keven)[2]
            rA = extract(L, Jd, q, omega, perts, Kall)[2]
            @printf("   %-4d %-5d %-4d %-8.2f   %-14.3e   %-14.3e\n", L, round(Int64, 2*Jd), q, omega, rE, rA)
        end
    end

    return( nothing )
end


function caseB()
    println("\nCASE B   THE EXTRACTED COEFFICIENT MATCHES THE CLOSED FORM AT EVERY RANK, ODD AND EVEN.")
    worst = 0.0;   tag = ""
    for L in [1,2,3], twoJ = 1:8
        Jd = twoJ/2.0
        for dJ = -L:L
            Jnud = Jd + dJ
            Jnud < 0.0 && continue
            for omega in [0.13, 0.41], dE in [1.7, -2.3]
                perts = [(Jnud, dE, 0.73)]
                Kall  = [K for K = 0:min(2L, round(Int64, 2*Jd))]
                sol   = extract(L, Jd, 1, omega, perts, Kall)[1]
                for (i, K) in enumerate(Kall)
                    pred  = coeffQ(L, K, Jd, Jnud, 1) * 0.73 * freqFactor(K, dE, omega)
                    scale = max(abs(pred), abs(sol[i]))
                    scale < 1.0e-14 && continue
                    d = abs(sol[i] - pred)/scale
                    if  d > worst
                        worst = d
                        tag   = @sprintf("L=%d 2J=%d 2J'=%d K=%d omega=%.2f dE=%.1f", L, twoJ,
                                         round(Int64, 2*Jnud), K, omega, dE)
                    end
                end
            end
        end
    end
    @printf("\n   worst relative deviation of extracted vs closed-form coefficient:  %.3e   (%s)\n", worst, tag)

    return( nothing )
end


function caseC()
    println("\nCASE C   THE ODD-K FREQUENCY FACTOR IS  omega/(dE^2-omega^2),  NOT  dE/(dE^2-omega^2).")
    Jd = 2.5;   L = 1;   perts = [(1.5, 1.9, 0.6), (3.5, 2.4, 0.9)]
    @printf("\n   %-10s %-18s %-18s %-18s\n", "omega", "alpha^(E1,1)", "/ [omega/(..)]", "/ [dE/(..)]")
    for omega in [0.0, 0.05, 0.1, 0.2, 0.4]
        a1 = extract(L, Jd, 1, omega, perts, [0,1,2])[1][2]
        rA = sum( coeffQ(L,1,Jd,p[1],1)*p[3]*omega/(p[2]^2 - omega^2)  for p in perts )
        rB = sum( coeffQ(L,1,Jd,p[1],1)*p[3]*p[2]/(p[2]^2 - omega^2)   for p in perts )
        @printf("   %-10.3f %-18.8e %-18s %-18s\n", omega, a1,
                abs(rA) > 1.0e-14 ? @sprintf("%.10f", a1/rA) : "  (0/0)",
                abs(rB) > 1.0e-14 ? @sprintf("%.10f", a1/rB) : "  (0/0)")
    end

    return( nothing )
end


function caseD(rng)
    println("\nCASE D   OUT OF SAMPLE:  random multi-perturber sets, both helicities and linear, every substate.")
    worst = 0.0;   tag = ""
    for trial = 1:40
        L     = rand(rng, 1:3);     Jd = rand(rng, 1:9)/2.0;    n = rand(rng, 1:7)
        perts = randPerts(rng, L, Jd, n)
        omega = 0.05 + 0.3*rand(rng)
        for q in [-1, 0, 1], M in Ms(Jd)
            direct = substateAlpha(L, Jd, M, q, omega, perts)
            pred   = sum( shapeA(Jd, K, M) * sum( coeffQ(L, K, Jd, p[1], q)*p[3]*freqFactor(K, p[2], omega)
                                                  for p in perts )
                          for K = 0:min(2L, round(Int64, 2*Jd)) )
            scale  = max(abs(direct), abs(pred))
            scale < 1.0e-12 && continue
            d = abs(direct - pred)/scale
            if  d > worst
                worst = d
                tag   = @sprintf("L=%d 2J=%d q=%d 2M=%d nPert=%d", L, round(Int64, 2*Jd), q, round(Int64, 2*M), n)
            end
        end
    end
    @printf("\n   worst relative deviation over 40 random sets, 3 helicities, every substate:  %.3e   (%s)\n", worst, tag)

    return( nothing )
end


function caseE(rng)
    println("\nCASE E   THE MODULE ITSELF.")
    bad = 0;   tested = 0
    for L in [1,2,3,4], twoJ = 1:9, twoJnu = 1:12, K = 0:8
        MP.rankCoefficient(L, K, jj(twoJ), jj(twoJnu)) === oldCoeff(L, K, jj(twoJ), jj(twoJnu)) || (bad += 1)
        tested += 1
    end
    @printf("\n   BIT-IDENTITY at the default helicity = 0:  %d of %d coefficients differ in any bit\n", bad, tested)
    @printf("   allowedRanks(E1, J=1/2):  helicity 0 -> %s   helicity +1 -> %s   (an alkali ground level)\n",
            MP.allowedRanks(E1, jj(1)), MP.allowedRanks(E1, jj(1); helicity=1))
    @printf("   allowedRanks(E2, J=5/2):  helicity 0 -> %s   helicity -1 -> %s\n",
            MP.allowedRanks(E2, jj(5)), MP.allowedRanks(E2, jj(5); helicity=-1))

    println("\n   polarizationFactor(L,K,q):")
    @printf("   %-5s %-5s %-11s %-11s %-11s\n", "L", "K", "q=-1", "q=0", "q=+1")
    for L in [1,2], K = 0:2L
        @printf("   %-5d %-5d %-11.6f %-11.6f %-11.6f\n", L, K,
                MP.polarizationFactor(L,K,-1), MP.polarizationFactor(L,K,0), MP.polarizationFactor(L,K,1))
    end

    worst = 0.0;   tag = ""
    for trial = 1:60
        L     = rand(rng, 1:3);   twoJ = rand(rng, 1:9);   Jd = twoJ/2.0
        perts = randPerts(rng, L, Jd, rand(rng, 1:6))
        omega = 0.05 + 0.3*rand(rng)
        for q in [-1, 0, 1], twoM = -twoJ:2:twoJ
            direct = substateAlpha(L, Jd, twoM/2.0, q, omega, perts)
            reb    = 0.0
            for K in MP.allowedRanks(EmMultipole(L, true), jj(twoJ); helicity=q)
                aK = sum( MP.rankCoefficient(L, K, jj(twoJ), jj(round(Int64, 2*p[1])); helicity=q) * p[3] *
                          (isodd(K) ? MP.frequencyFactorOdd(p[2], omega, false)
                                    : MP.frequencyFactor(p[2], omega, false))   for p in perts )
                reb = reb + MP.rankShape(jj(twoJ), K, mm(twoM)) * aK
            end
            scale = max(abs(direct), abs(reb))
            scale < 1.0e-12 && continue
            d = abs(direct - reb)/scale
            if  d > worst
                worst = d;   tag = @sprintf("L=%d 2J=%d q=%d 2M=%d", L, twoJ, q, twoM)
            end
        end
    end
    @printf("\n   the substate identity rebuilt from the MODULE's own functions:  worst %.3e   (%s)\n", worst, tag)

    println("\n   THE GUARDS")
    @printf("   frequencyFactorOdd(1.7, 0.0, false) = %.1f   (an odd rank is zero in the static limit)\n",
            MP.frequencyFactorOdd(1.7, 0.0, false))
    try
        MP.frequencyFactorOdd(1.7, 0.3, true);    println("   imaginary axis: NO ERROR RAISED -- WRONG")
    catch
        println("   imaginary axis: raises as designed")
    end
    try
        MP.allowedRanks(E1, jj(3); helicity=2);    println("   helicity = 2:   NO ERROR RAISED -- WRONG")
    catch
        println("   helicity = 2:   raises as designed")
    end

    return( nothing )
end


function caseF()
    println("\nCASE F   THE DEGENERACY ANCHOR, and then the full property path on the approved H 1s case.")
    c1h = MP.rankCoefficient(1, 1, jj(1), jj(1); helicity=1)
    c3h = MP.rankCoefficient(1, 1, jj(1), jj(3); helicity=1)
    c1s = MP.rankCoefficient(1, 0, jj(1), jj(1); helicity=1)
    c3s = MP.rankCoefficient(1, 0, jj(1), jj(3); helicity=1)
    @printf("\n   c(E1,K=1,1/2,1/2) = %+.15f    c(E1,K=1,1/2,3/2) = %+.15f\n", c1h, c3h)
    @printf("   rank-1 statistical sum  1*c(1/2) + 2*c(3/2) = %+.3e      <-- must be 0\n", c1h + 2*c3h)
    @printf("   rank-0 statistical sum  1*c(1/2) + 2*c(3/2) = %+.15f   <-- must NOT be 0\n", c1s + 2*c3s)

    ni   = Nuclear.Model(1.0, PointNucleus())
    grid = Radial.Grid(Radial.Grid(false), rnt=4.0e-6, h=5.0e-2, hp=4.0e-2, rbox=80.0)
    scf  = Basics.NuclearField();    asf = AsfSettings(AsfSettings(); scField=scf)
    gMp  = redirect_stdout(devnull) do
        generate(Representation("np perturbers", ni, grid,
                   [Configuration("2p"), Configuration("3p"), Configuration("4p")],
                   MeanFieldMultiplet(MeanFieldSettings(scf))), output=true)[ResultKeys.MeanFieldMultiplet]   end

    for  hel in [0, 1]
        set = MP.Settings(MP.Settings(); multipoles=[E1], gMultiplet=gMp, omegas=[0.05], helicity=hel,
                          printBefore=false, levelSelection=LevelSelection())
        wa  = Atomic.Computation(Atomic.Computation(), name="H 1s, helicity $hel", grid=grid, nuclearModel=ni,
                  configs=[Configuration("1s")], asfSettings=asf, propertySettings=[set])
        outcomes = redirect_stdout(devnull) do;   perform(wa; output=true)[MP.Settings]   end
        println("\n   helicity = $hel")
        for  o in outcomes, (mp, K, om, a) in o.alphas
            @printf("      alpha^(%s,%d)  omega = %.3f   =  %+.12e\n", string(mp), K, om, a.Babushkin)
        end
    end
    es = sort( [ (Basics.twice(l.J), l.energy)  for l in gMp.levels ], by = x -> x[2] )
    @printf("\n      2p fine-structure splitting:  %.3e a.u.  --  the scale alpha^(E1,1) must carry\n",
            abs(es[2][2] - es[1][2]))

    return( nothing )
end


function main()
    rng = MersenneTwister(20102026)
    println("\n", "="^118)
    println("ODD-RANK (K = 1, the VECTOR) MULTIPOLE POLARIZABILITY -- derivation and validation.  seed = 20102026")
    println("="^118)
    if  "ph" in ARGS
        caseF()
    else
        caseA(rng);   caseB();   caseC();   caseD(rng);   caseE(rng);   caseF()
    end
    println("\n", "="^118)

    return( nothing )
end

main()
