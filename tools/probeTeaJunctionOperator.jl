#
# tools/probeTeaJunctionOperator.jl
#
# Measures the SECOND-ORDER JUNCTION <n| O |i> of the argon two-electron Auger case with O = V^(e-e) and with
# O = H - E, for several mean fields, and reports the ratio.  Written 23-Oct-2026 for priority items 93 and 98;
# it is a diagnostic, not a test, and it is cheap -- it builds only the two BOUND multiplets and no continuum
# orbital, so it runs in a few minutes where the full branch (examples/example-De.jl, branch f) takes ~5.
#
# WHAT IT ESTABLISHED, 23-Oct-2026, and why each part is here:
#   * NuclearField gives ratio = 1.0000 on all five contributing intermediates.  That is the CALIBRATION and it
#     is the reason to trust the rest: all three configurations carry 16 electrons, so a pure nuclear potential
#     is identical for both multiplets and its orbitals are eigenfunctions of h itself, whence <b|h|a> = 0 and
#     H - E must reduce to V^(e-e) exactly.  A ratio other than 1 here is a defect in the one-body term.
#   * HSField gives 0.057-0.067 and ALField 0.130-0.141 on the one-orbital junctions.  The residual GREW by
#     2.2x on the non-local field, which refuted the expectation that it would collapse there.
#   * The orbital overlaps <initial|Green> are 0.995-1.000 and essentially IDENTICAL for the two fields
#     (3p_3/2: 4.887e-3 against 4.816e-3 from unity), so the growth is NOT non-orthogonality between the two
#     SCF sets.  Both explanations offered on the day were wrong; the overlap block is what killed the second.
#
# THE READING THAT SURVIVES: <b|h|a> = -<b|U|a> is the OFF-DIAGONAL element of whichever mean field made the
# orbitals, and Brillouin closes only for the Fock operator OF THAT STATE.  ALField is an AVERAGE-LEVEL field,
# so it is not that operator for an individual CSF pair any more than a local Slater exchange is.  The H - E
# term is therefore field-dependent and vanishes for NO field JAC offers, which argues that it must always be
# carried rather than that it is a small correction -- the opposite of the conclusion the run was expected to
# support.
#
# Run:  julia --project=<JAC-root> tools/probeTeaJunctionOperator.jl
#

using JenaAtomicCalculator
const JAC = JenaAtomicCalculator
using Printf

# [AR-JUNCTION-FIELDS]  ONE question: does the H - E residual at the second-order junction collapse when the
# mean field's exchange becomes NON-LOCAL?  The residual is <b|(U_HF - U)|a>, which is identically zero only
# for U = U_HF.  ONE variable per run: the mean field.  Everything else is the branch-f setup verbatim.
#
# PREDICTION, WRITTEN BEFORE THE RUN AND KEPT BECAUSE IT WAS REFUTED:  ratio = <n|H-E|i> / <n|V_ee|i> is
# 0.057-0.067 on HSField, whose Slater exchange is LOCAL;  on ALField, taken for "the non-local Dirac-Fock
# field", it was to fall substantially toward zero.  **IT ROSE, to 0.130-0.141.**  The error in the prediction
# was the premise, not the arithmetic: ALField is an AVERAGE-LEVEL field, and Brillouin's theorem closes only
# for the Fock operator OF THE STATE IN QUESTION, which an average over a configuration's levels is no more
# than a local exchange is.  Read the ratio here as a measurement of the off-diagonal element of whichever
# mean field made the orbitals, NOT as a distance from some ideal at which it would vanish.

grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)
nm   = Nuclear.Model(18.)
pot  = Nuclear.nuclearPotential(nm, grid)
opV  = CoulombInteraction()
opH  = Basics.HminusE(CoulombInteraction())

gConfigs = [Configuration("1s^2 2s^2 2p^5 3s^2 3p^5"), Configuration("1s^2 2s^2 2p^6 3s^2 3p^4")]
iConfigs = [Configuration("1s^2 2s^2 2p^4 3s^2 3p^6")]

function junction(scField, label)
    asf  = AsfSettings(AsfSettings(), scField=scField)
    gMp  = SelfConsistent.performSCF(gConfigs, nm, grid, asf; printout=false)
    iMp  = SelfConsistent.performSCF(iConfigs, nm, grid, asf; printout=false)
    subshells = Basics.merge(iMp.levels[1].basis.subshells, gMp.levels[1].basis.subshells)
    iLevel    = Level(iMp.levels[1], subshells)
    nMult     = Multiplet(gMp, subshells)
    symi = LevelSymmetry(iLevel.J, iLevel.parity);   eni = iLevel.energy
    ni   = length(iLevel.basis.csfs);                nn  = length(nMult.levels[1].basis.csfs)
    println("\n=== $label :  initial level J^P = $symi,  energy = $eni Ha ===")
    @printf("  %-4s %-9s %14s %14s %10s\n", "nLev", "Delta E", "V^(e-e)", "H - E", "ratio")
    for  nLevel in nMult.levels
        symn = LevelSymmetry(nLevel.J, nLevel.parity)
        if  symn != symi                               continue    end
        if  abs(nLevel.energy - eni) <= 1.0e-8         continue    end
        jv = 0.;   jh = 0.
        for  tp = 1:nn
            if  nLevel.mc[tp] == 0.   continue    end
            for  s = 1:ni
                if  LevelSymmetry(iLevel.basis.csfs[s].J, iLevel.basis.csfs[s].parity) != symi   continue   end
                if  iLevel.mc[s] == 0.   continue    end
                wv = ManyElectron.matrixElement_Vee(opV, nLevel.basis, tp, iLevel.basis, s, grid)
                wh = ManyElectron.matrixElement_HminusE(opH, nLevel.basis, tp, iLevel.basis, s, grid, pot)
                jv = jv + nLevel.mc[tp] * wv * iLevel.mc[s]
                jh = jh + nLevel.mc[tp] * wh * iLevel.mc[s]
            end
        end
        if  abs(jv) > 1.0e-14
            @printf("  %-4d %9.4f %14.5e %14.5e %10.4f\n", nLevel.index, eni - nLevel.energy, jv, jh, jh/jv)
        end
    end
end



println("\n=== PROBE FINISHED ===")

# [AR-JUNCTION-CALIBRATION]  ONE question: is matrixElement_HminusE itself correct, or is the growth of the
# residual from HSField to ALField a bug?  NuclearField decides it.  All three configurations here carry 16
# electrons, so a pure nuclear potential is IDENTICAL for the initial and the Green multiplet: one common,
# orthonormal, hydrogenic orbital set.  Those orbitals are eigenfunctions of h itself, so <b|h|a> = 0 for
# b != a and the one-body term must VANISH -- i.e. H - E must reproduce V^(e-e) EXACTLY, ratio = 1.0000.
# This is the U = 0 limit of <b|(U_HF - U)|a>, which is then the whole mean field, which is V^(e-e)'s
# spectator sum.  A ratio that is NOT 1 here is a defect in the one-body term, and nothing else.
junction(Basics.HSField(), "HSField  -- self-consistent, LOCAL Slater exchange")
junction(Basics.ALField(), "ALField  -- average-level, non-local exchange")
junction(Basics.NuclearField(), "NuclearField -- common hydrogenic set, MUST give ratio = 1.0000")

# And the size of the non-orthogonality itself, which is the quantity the dropped -E<n|i> term would carry.
for  (scField, label) in [(Basics.HSField(), "HSField"), (Basics.ALField(), "ALField")]
    asf = AsfSettings(AsfSettings(), scField=scField)
    gMp = SelfConsistent.performSCF(gConfigs, nm, grid, asf; printout=false)
    iMp = SelfConsistent.performSCF(iConfigs, nm, grid, asf; printout=false)
    go  = gMp.levels[1].basis.orbitals;   io = iMp.levels[1].basis.orbitals
    println("\n=== orbital overlap <initial | Green>, $label  (1.0000 = one common set) ===")
    for  sh in sort(collect(keys(io)), by = x -> (x.n, x.kappa))
        if  haskey(go, sh)
            ov = RadialIntegrals.overlap(io[sh], go[sh], grid)
            @printf("   %-9s  <i|n> = %12.8f   deviation = %10.3e\n", string(sh), ov, abs(ov) - 1.0)
        end
    end
end
println("\n=== CALIBRATION FINISHED ===")
