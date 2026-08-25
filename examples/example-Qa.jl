
println("Qa) Muonic atoms: muon binding energies, muonic X-rays, and the two screenings.")

using Printf

# WRITTEN 25-Aug-2026, first implementation of module-MuonicAtom.jl.  The module that stood here before computed
# nothing -- nine functions, every one returning zero -- and was never included into JAC.  It is replaced entirely.
#
# THE MODEL, agreed with the maintainer and the reason the module is small.  A muon is an electron with 206.77
# times the mass: same charge, same spin, the same Dirac equation.  It is treated as a SINGLE PARTICLE in the
# field of the nucleus plus, if wanted, the mean field of the electrons.  Muon and electrons never appear together
# in one amplitude and there is no exchange between them; each sees only the other's average charge.  That is what
# makes the whole thing a one-particle problem, and it is a good approximation precisely because the muon orbits
# some two hundred times closer to the nucleus than any electron does.
#
# WHAT HAD TO CHANGE IN JAC, and it was less than expected.  Written out, the Dirac matrix in a B-spline basis is
#
#         [   V          c D_kappa^-  ]                      [   V          c D_kappa^-  ]
#         [ c D_kappa^+   V - 2 c^2   ]     becomes          [ c D_kappa^+   V - 2 m c^2 ]
#
# -- the kinetic blocks carry the speed of light but NOT the mass, so the particle mass enters in exactly one
# place.  One companion had to move with it: Bsplines.findPositiveBranchStart separates real bound states from the
# unphysical negative-energy branch using a threshold of -1.999 c^2, and that threshold scales with the mass too.
# It matters more than it looks.  A muon 1s level in lead lies near 10 MeV = 3.7e5 Hartree, far BELOW the electron
# threshold of 3.75e4, so with the electron value every muon bound state of a heavy atom is discarded as spurious
# -- silently, because that function returns an index and not an error.


if  true
    # Last visit:      25-Aug-2026
    # Last successful: 25-Aug-2026
    #
    # Branch a: THE KNOWN-ANSWER TEST, which needs no literature at all.  For a POINT nucleus the Dirac energy is
    #   m c^2 f(Z alpha), and the fine-structure constant does not depend on which particle is bound -- so a muon
    #   level must be EXACTLY 206.7683 times the corresponding electron level at the same Z.  That is an identity,
    #   not an approximation, and it tests the mass threading end to end: the analytic formula, the B-spline
    #   matrix, and the branch selection that decides which eigenvalue is a bound state at all.
    #
    # REPORT (25-Aug-2026): at Z = 10 the analytic ratio is 206.768283 for every subshell, to all printed digits.
    #   The SOLVER then reproduces the analytic value with a relative deviation of
    #        electron   1s 5.8e-09   2s 1.3e-07   2p_1/2 3.1e-08   2p_3/2 4.5e-08
    #        muon       1s 6.3e-09   2s 2.0e-07   2p_1/2 5.2e-08   2p_3/2 6.7e-08
    #   i.e. the same accuracy for both particles, which is what must happen for one equation with a scaled mass.
    #   THE FIRST ATTEMPT FAILED THIS TEST and the failure is worth recording: the muon energies came out POSITIVE
    #   and meaningless.  The mass keyword had been added to every function that needed it, but one CALL SITE
    #   inside Bsplines.generateOrbitals did not pass it on, so the branch selection ran with the electron mass and
    #   picked an index in the positive continuum.  Neither a parse check nor loading the package showed anything;
    #   only running it did.  Adding a parameter and passing a parameter are different claims.
    #
    #   NOTE THE GRID.  A muon orbit is ~207 times smaller than the electron one, so Rule 12 -- the box must match
    #   the orbitals -- applies with the same force to a different scale.  The electron here uses rbox = 6 a.u.,
    #   the muon rbox = 0.05 a.u.  Using the electron's box for the muon starves the basis exactly as a box that
    #   is too large always does.
    setDefaults("print summary: open", "zzz-MuonicAtom-Qa-scaling.sum")

    mmu = Defaults.getDefaults("mass: muon")
    Z   = 10.0;    shs = [Subshell("1s_1/2"), Subshell("2s_1/2"), Subshell("2p_1/2"), Subshell("2p_3/2")]
    println("\n  (i) the analytic point-nucleus energies, and their ratio [Hartree]\n")
    @printf("      %-10s %16s %16s %12s\n", "subshell", "electron", "muon", "ratio")
    for sh in shs
        ee = Basics.computeDiracEnergy(sh, Z);   em = Basics.computeDiracEnergy(sh, Z; mass=mmu)
        @printf("      %-10s %16.8f %16.4f %12.6f\n", string(sh), ee, em, em/ee)
    end
    #
    println("\n  (ii) the B-spline SOLVER against that formula, point nucleus\n")
    for (label, mass, rbox) in [("electron", 1.0, 6.0), ("muon", mmu, 0.05)]
        nm   = Nuclear.Model(Z, Nuclear.PointNucleus(), 0., 0., AngularJ64(0), 0., 0., 0.)
        grid = Radial.Grid(Radial.Grid(false); rnt=1.0e-9, h=5.0e-2, hp=0., rbox=rbox)
        prim = Bsplines.generatePrimitives(grid)
        orbs = Bsplines.generateOrbitalsHydrogenic(shs, nm, prim; printout=false, mass=mass)
        @printf("      %-8s  box = %-8.3g %14s %18s %12s\n", label, rbox, "solver", "analytic", "rel.dev.")
        for sh in shs
            es = orbs[sh].energy;   ea = Basics.computeDiracEnergy(sh, Z; mass=mass)
            @printf("      %-10s %26.8f %18.8f %12.2e\n", string(sh), es, ea, abs(es-ea)/abs(ea))
        end
        println("")
    end
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      25-Aug-2026
    # Last successful: 25-Aug-2026
    #
    # Branch b: MUONIC X-RAYS, AND WHY THEY MEASURE NUCLEAR RADII.  The muon 1s orbit has a radius of about
    #   256/Z fm, while a nucleus of mass number A has a radius of about 1.2 A^(1/3) fm.  Above Z ~ 40 the muon
    #   therefore orbits INSIDE the nuclear charge, where it no longer feels the full -Z/r, and its binding is
    #   reduced enormously.  This branch computes the 1s and 2p levels twice, once for a point nucleus and once
    #   for a Fermi charge distribution, and forms the 2p -> 1s muonic X-ray from the second.
    #
    # REPORT (25-Aug-2026): the shift grows from a rounding error to a factor of two across the periodic table.
    #        Z   el     1s point     1s Fermi     2p->1s [keV]    finite-size shift     measured 2p->1s
    #        6   C        -101.3       -100.9          75.6            0.4 %              ~75.3 keV
    #       12   Mg       -405.9       -397.6         296.3            2.0 %              ~296  keV
    #       26   Fe      -1919.2      -1723.1        1246.8           10.2 %              ~1257 keV
    #       82   Pb     -21003.8     -10507.7        5908.3           50.0 %              ~5963 keV
    #   (energies in keV; the measured column is the standard muonic K-alpha, quoted from memory of the muonic
    #   X-ray literature and NOT from a table checked here -- see the caveat below.)
    #
    #   ALL FOUR AGREE WITH MEASUREMENT TO ABOUT 1%, which for a computation with no free parameter is a real
    #   result: the only inputs are Z, the mass number, and the muon mass.  IN LEAD THE FINITE NUCLEUS HALVES THE
    #   BINDING -- -21004 keV becomes -10508 keV -- so a muonic X-ray energy is far more sensitive to the nuclear
    #   charge distribution than to anything atomic, which is exactly why the technique exists.
    #
    #   THE RESIDUAL IS THE RIGHT SIZE AND THE RIGHT SIGN.  The computed energies lie ~1% BELOW the measured ones
    #   for the heavy cases, and the largest omission here is VACUUM POLARISATION, which is the dominant QED
    #   correction in muonic atoms and increases the binding.  That is consistent rather than confirmed: no
    #   vacuum-polarisation calculation was done, so this is an argument about direction and magnitude, not a
    #   demonstration.  NOT dated on the 1% agreement alone -- the measured column needs a proper reference before
    #   any of these numbers is quoted, and the nuclear radius used is a default from a mass formula rather than a
    #   measured charge radius.
    setDefaults("print summary: open", "zzz-MuonicAtom-Qb-xrays.sum")

    shs = [Subshell("1s_1/2"), Subshell("2p_1/2"), Subshell("2p_3/2")]
    println("\n      Z   el      1s point      1s Fermi      2p->1s        shift   [all in keV]\n")
    for (Z, el, A) in [(6.,"C",12.), (12.,"Mg",24.), (26.,"Fe",56.), (82.,"Pb",208.)]
        grid = MuonicAtom.recommendedGrid(Z, 2)
        nmP  = Nuclear.Model(Z, Nuclear.PointNucleus(), A, 0., AngularJ64(0), 0., 0., 0.)
        nmF  = Nuclear.Model(Z, A)
        oP   = MuonicAtom.computeOrbitals(shs, nmP, grid)
        oF   = MuonicAtom.computeOrbitals(shs, nmF, grid)
        keV(e) = Defaults.convertUnits("energy: from atomic to eV", e) * 1.0e-3
        e1P = keV(oP[shs[1]].energy);   e1F = keV(oF[shs[1]].energy)
        @printf("     %2.0f   %-3s %12.1f %13.1f %12.1f %11.1f %%\n",
                Z, el, e1P, e1F, keV(oF[shs[3]].energy) - e1F, 100*(e1F-e1P)/abs(e1P))
    end
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      25-Aug-2026
    # Last successful: 25-Aug-2026
    #
    # Branch c: THE OTHER SCREENING -- what the MUON does to the ELECTRONS.  A muonic atom is chemically almost
    #   the element one place LOWER in the periodic table, and the reason is geometric rather than subtle: the
    #   muon sits some two hundred times closer to the nucleus than any electron, so the electrons cannot resolve
    #   it from the nucleus and simply see one unit of charge glued on.  MuonicAtom.screeningPotential forms that
    #   contribution properly, as the Hartree potential of the muon's own charge density, and the Z-1 statement is
    #   then its LARGE-r LIMIT rather than an assumption put in by hand.
    #
    # REPORT (25-Aug-2026): for the 1s muon of lead the screening charge Zr(r) climbs from essentially nothing at
    #   the origin to exactly 1 at large r:
    #        r = 2.6e-10 a.u.   Zr = 0.000002          r = 2.6e-07 a.u.   Zr = 0.002227
    #        r = 1.5e-08 a.u.   Zr = 0.000127          r = 1.7e-05 a.u. (0.91 fm)  Zr = 0.148070
    #                                                  r = 1.1e-03 a.u. (61 fm)    Zr = 1.000000
    #   The 1s ELECTRON of lead has its maximum near 1.2e-2 a.u., an order of magnitude beyond the last row, so at
    #   every distance an electron actually occupies the muon screening is a full unit of charge.  Hence Z-1.
    #
    #   THE LARGE-r LIMIT IS THE CHECK, and it needs no reference: a bound particle of unit charge must screen
    #   exactly one unit at infinity, and any other value means the density is not normalised or the integral is
    #   wrong.  IT FAILED THE FIRST TIME, at 1.000374, and the 0.04% was not physics -- it was a hand-rolled
    #   trapezoidal rule on an exponential mesh.  Replacing it by the grid's own Gauss-Legendre weights, grid.wr,
    #   gives 1.000000.  A quadrature error and a physics error look identical in the answer; only an exact limit
    #   separates them, which is the whole reason to have one.
    setDefaults("print summary: open", "zzz-MuonicAtom-Qc-screening.sum")

    Z = 82.;  A = 208.
    grid = MuonicAtom.recommendedGrid(Z, 2)
    nm   = Nuclear.Model(Z, A)
    orbs = MuonicAtom.computeOrbitals([Subshell("1s_1/2")], nm, grid)
    scr  = MuonicAtom.screeningPotential(orbs[Subshell("1s_1/2")], grid)
    println("\n  the muon's screening charge Zr(r), lead, muon in 1s:\n")
    println("        r [a.u.]        r [fm]       Zr(r)")
    for  frac  in [0.02, 0.05, 0.1, 0.2, 0.4, 0.7, 1.0]
        i = max(2, round(Int64, frac*length(grid.r)))
        @printf("     %12.4e %14.4e %11.6f\n", grid.r[i], grid.r[i]*0.529177e5, scr.Zr[i])
    end
    @printf("\n     large-r limit  Zr = %.6f   (must be 1: one bound particle screens one unit of charge)\n",
            scr.Zr[end])
    setDefaults("print summary: close", "")
    #
end
