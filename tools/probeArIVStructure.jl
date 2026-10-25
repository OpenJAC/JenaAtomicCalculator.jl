#
# tools/probeArIVStructure.jl
#
# Compares correlation models for the Ar IV `3s^2 3p^3` final state of the two-electron Auger case
# (examples/example-De.jl branch f) against the MEASURED term splittings.  Structure only -- no continuum
# orbital, no Green multiplet, no second-order sum -- so it costs minutes and settles which model is worth
# paying for before any TEA run.  Written 25-Oct-2026 for priority item 93.
#
# REFERENCE, double-sourced:  4S*->2D*_3/2 = 2.615 eV and 4S*->2P*_1/2 = 4.321 eV, being the differences of
# the transition energies in Zitnik et al., PRA 93, 021401(R) (2016) Table I, and matching the NIST Ar IV
# levels (21090 and 34855 cm^-1).
#
# WHAT IT ESTABLISHED, 25-Oct-2026 -- A NEGATIVE RESULT, AND THE RATIO IS WHAT PRODUCED IT:
#
#     model                      4S->2D       4S->2P     2P/2D     vs LS 5/3
#     measured                    2.615        4.321     1.652       -0.9 %
#     M0  3s^2 3p^3 only         +18.3 %      +19.8 %    1.674       +0.4 %
#     M1  + 3p^5                 +18.2 %       +0.9 %    1.411      -15.4 %
#     M2  + 3s^2 3p 3d^2         +20.0 %      +10.8 %    1.526       -8.4 %
#     M3  + 3s 3p^3 3d           +21.5 %      +41.2 %    1.920      +15.2 %
#     M4  + all three            +15.7 %       +7.2 %    1.531       -8.1 %
#
# JUDGE BY THE RATIO, NOT BY THE ERROR COLUMN.  For p^3 the Slater algebra fixes E(2D)-E(4S) = 9*F2 and
# E(2P)-E(4S) = 15*F2, so the ratio MUST be 5/3 whatever F2 is.  **Only the UNCORRELATED M0 has it.**  M1
# reads as a triumph at +0.9 % on 2P and is a distortion: `3p^5` is a single-term (2P) configuration and can
# mix with 2P alone.  Every n=3 configuration admitted by parity couples selectively to 2P, and the 2D
# splitting is immovable at +15.7 to +21.5 % in all five models.  On the error columns alone, M1 and M4 would
# have been reported as progress.
#
# THE DIAGNOSIS: `F^2(3p,3p)` is uniformly ~19 % too large with the LS structure otherwise exact, i.e. missing
# screening -- the same effect the Basics.DFSField docstring records when it says a local-exchange field
# "agrees better by a CANCELLATION OF ERRORS".  Repairing it needs core-valence (2p/3s) correlation or n >= 4
# layers, which is a RAS calculation and not a configuration list.
#
# AND IT CANNOT CARRY branch f's factor ~2.3.  The tempting arithmetic is rate ~ F^4 and 1.19^4 = 2.0, but the
# FIRST-ORDER width (example-De.jl branch g, 1.12x measured) pins the 2p-3p integrals that actually drive the
# rate to ~6 % in amplitude, whereas the 19 % sits in `F^2(3p,3p)`, which governs the final-state TERM
# SPACING.  Different integral.
#
# THE BOX IS NOT THE LIMITING FACTOR HERE, measured rather than assumed: re-running at rbox = 24.0 (1582 grid
# points against 1302) reproduces every model to ALL FOUR DIGITS, because JAC's grid extends OUTWARD and
# leaves the inner exponential region -- which is all the 3p sees -- untouched.  The 10.0 a.u. run's own grid
# guard did still correctly flag the 3d-bearing models as under-boxed, so run this at 24.0 if the 3d models
# are to be believed at all.
#
# Run:  julia --project=<JAC-root> tools/probeArIVStructure.jl
#

using JenaAtomicCalculator
const JAC = JenaAtomicCalculator
using Printf

grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.9e-2, rbox = 24.0)
nm   = Nuclear.Model(18.)
asf  = AsfSettings(AsfSettings(), scField=Basics.HSField())

base = Configuration("1s^2 2s^2 2p^6 3s^2 3p^3")
models = [ ("M0  3s^2 3p^3 only (the branch)",        [base]),
           ("M1  + 3p^5",                             [base, Configuration("1s^2 2s^2 2p^6 3p^5")]),
           ("M2  + 3s^2 3p 3d^2",                     [base, Configuration("1s^2 2s^2 2p^6 3s^2 3p 3d^2")]),
           ("M3  + 3s 3p^3 3d",                       [base, Configuration("1s^2 2s^2 2p^6 3s 3p^3 3d")]),
           ("M4  + all three",                        [base, Configuration("1s^2 2s^2 2p^6 3p^5"),
                                                             Configuration("1s^2 2s^2 2p^6 3s^2 3p 3d^2"),
                                                             Configuration("1s^2 2s^2 2p^6 3s 3p^3 3d")]) ]

MEAS_D = 2.615;   MEAS_P = 4.321          # eV, 4S->2D_3/2 and 4S->2P_1/2
println("\n", "="^96)
@printf("  %-34s %6s %9s %9s %9s %9s\n", "model", "nLev", "4S->2D", "err", "4S->2P", "err")
@printf("  %-34s %6s %9.3f %9s %9.3f %9s\n", "MEASURED (Ar IV)", "-", MEAS_D, "-", MEAS_P, "-")
println("="^96)
for  (label, configs) in models
    try
        mp = SelfConsistent.performSCF(configs, nm, grid, asf; printout=false)
        lv = sort(mp.levels, by = x -> x.energy)
        e0 = lv[1].energy
        # the 3p^3 terms in order: 4S_3/2, 2D_3/2, 2D_5/2, 2P_1/2, 2P_3/2 -- take the first level of each J
        eD = nothing;   eP = nothing;   seen = Rational[]
        for  l in lv
            de = Defaults.convertUnits("energy: from atomic to eV", l.energy - e0)
            if  de > 1.0  &&  isnothing(eD)                  eD = de
            elseif  !isnothing(eD)  &&  de > eD + 0.8  &&  isnothing(eP)   eP = de    end
            if  !isnothing(eP)  break   end
        end
        @printf("  %-34s %6d %9.3f %+8.1f%% %9.3f %+8.1f%%\n", label, length(lv), eD, 100*(eD/MEAS_D-1),
                 eP, 100*(eP/MEAS_P-1))
    catch err
        @printf("  %-34s   FAILED: %s\n", label, sprint(showerror, err)[1:min(60,end)])
    end
end
println("\n=== ARIV STRUCTURE PROBE FINISHED ===")
