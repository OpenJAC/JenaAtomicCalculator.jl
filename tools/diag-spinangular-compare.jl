# Compare SpinAngular against SpinAngularNew, coefficient by coefficient, on bases chosen so that the
# comparison CAN FAIL.  Written 27-Aug-2026 after three caller swaps passed a verification set that could not.
#
# THE TWO THINGS A NAIVE CHECK MISSES, both of which cost a day:
#   1. a SAME-KAPPA, DIFFERENT-n substitution with l >= 1.  Two s subshells of different n have no subshell
#      between them, so l = 0 exercises nothing; `3d^2 + 3d 4d` exposes 26 coefficients at once.
#   2. CROSS-SYMMETRY CSF pairs.  SpinAngular ignores its own parity argument and emits rank-0 coefficients
#      between opposite-parity CSFs, which a scalar operator cannot connect; SpinAngularNew omits them and is
#      right.  A same-symmetry-only comparison hides the whole class.
#
# Run:  julia --project=. tools/diag-spinangular-compare.jl
using Printf, JenaAtomicCalculator
const B = JenaAtomicCalculator.Basics
const SA = JenaAtomicCalculator.SpinAngular;  const SAN = JenaAtomicCalculator.SpinAngularNew

const BASES = [("Be   2s^2 + 2p^2",            ["1s^2 2s^2", "1s^2 2p^2"]),
               ("C    2s^2 2p^2 + 2s 2p^3",    ["1s^2 2s^2 2p^2", "1s^2 2s^1 2p^3"]),
               ("Li   2p + 3p        [l>=1]",  ["1s^2 2p^1", "1s^2 3p^1"]),
               ("Ne   2p^6 + 2p^5 3p [l>=1]",  ["1s^2 2s^2 2p^6", "1s^2 2s^2 2p^5 3p^1"]),
               ("3d^2 + 3d 4d        [l>=1]",  ["1s^2 2p^6 3d^2", "1s^2 2p^6 3d^1 4d^1"])]

function csfset(confs)
    rel = JenaAtomicCalculator.ManyElectron.ConfigurationR[]
    for c in confs  append!(rel, B.generateConfigurations(B.RelativisticConfigurations(), Configuration(c)))  end
    sub  = B.generateSubshellList(rel)
    csfs = JenaAtomicCalculator.ManyElectron.CsfR[]
    for rc in rel  append!(csfs, B.generateCsfRs(rc, sub))  end
    return( (csfs, sub) )
end

key1(c) = (c.nu, string(c.a), string(c.b))
key2(c) = (c.nu, string(c.a), string(c.b), string(c.c), string(c.d))

function compare(label, confs)
    csfs, sub = csfset(confs)
    for (what, sameSymOnly) in (("same symmetry", true), ("CROSS symmetry", false))
        miss1 = extra1 = diff1 = miss2 = extra2 = diff2 = raised = pairs = 0
        for r in csfs, s in csfs
            same = (r.J == s.J  &&  r.parity == s.parity)
            same == sameSymOnly || continue
            pairs += 1
            for (kind, oldc, newc, kf, val) in
                (("1p", () -> SA.computeCoefficients(SA.OneParticleOperator(0, B.plus, true), r, s, sub),
                        () -> SAN.computeCoefficients(SAN.OneParticleOperator(0, B.plus), r, s, sub), key1, c -> c.T),
                 ("2p", () -> SA.computeCoefficients(SA.TwoParticleOperator(0, B.plus, true), r, s, sub),
                        () -> SAN.computeCoefficients(SAN.TwoParticleOperator(0, B.plus), r, s, sub), key2, c -> c.V))
                o = oldc()
                n = try newc()  catch;  raised += 1;  continue  end
                ko = Dict(kf(c) => val(c) for c in o if abs(val(c)) > 1.0e-12)
                kn = Dict(kf(c) => val(c) for c in n if abs(val(c)) > 1.0e-12)
                m = length(setdiff(keys(ko), keys(kn)));   e = length(setdiff(keys(kn), keys(ko)))
                d = count(k -> abs(ko[k] - kn[k]) > 1.0e-9 * max(1.0, abs(ko[k])), intersect(keys(ko), keys(kn)))
                kind == "1p" ? (miss1 += m; extra1 += e; diff1 += d) : (miss2 += m; extra2 += e; diff2 += d)
            end
        end
        @printf("  %-30s %-14s pairs %4d | 1p miss %3d extra %3d differ %3d | 2p miss %3d extra %3d differ %3d | raises %d\n",
                label, what, pairs, miss1, extra1, diff1, miss2, extra2, diff2, raised)
    end
end

println("SpinAngular vs SpinAngularNew, rank 0.  'miss' = emitted by SpinAngular only; 'extra' = by SpinAngularNew only.")
println("EXPECTED: all zero on same-symmetry pairs.  On cross-symmetry pairs a non-zero 1p 'miss' is CORRECT --")
println("those are the parity-forbidden terms SpinAngular emits and SpinAngularNew rightly does not.\n")
for (l, c) in BASES    compare(l, c)    end
