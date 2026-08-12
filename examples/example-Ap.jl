#
println("Ap) Apply & test the three-parameter Fermi nucleus: the shape of the nuclear charge distribution.")

using Printf

# WRITTEN 12-Aug-2026.  JAC's two-parameter Fermi nucleus, rho(r) = rho_0 / (1 + exp((r-c)/a)), is flat in
# the nuclear interior by construction.  The three-parameter form adds one dimensionless parameter,
#
#     rho(r) = rho_0 (1 + w r^2/c^2) / (1 + exp((r-c)/a)),
#
# with w > 0 a central bump and w < 0 a CENTRAL DEPRESSION.  In medium-mass nuclei the depression is a
# shell effect (the 34Si bubble); in the superheavy region it is driven by the Coulomb repulsion pushing
# protons outward, which makes it robust rather than model-dependent -- an appreciable semibubble in the
# proton density is predicted for 294Og.  This is also the form tabulated by the elastic electron-scattering
# compilations (de Vries, de Jager & de Vries, At.Data Nucl.Data Tables 36 (1987) 495) and used in the
# analysis of muonic atoms.
#
# THE POINT OF THE COMPARISON BELOW.  Each branch holds the rms radius FIXED and varies only w.  That is
# deliberate: the leading finite-nuclear-size shift depends on <r^2>, so holding <r^2> fixed cancels it and
# what remains isolates the SHAPE of the distribution -- essentially the fourth moment.  It is therefore a
# small effect for an electron, and the branches report it both absolutely and as a fraction of the total
# finite-size shift, so that the size of the residual is visible rather than implied.


if  true
    # Last visit:  12-Aug-2026
    # Last successful:  12-Aug-2026
    #
    ## (a)  The closed-form rms radius, checked against the tabulated three-parameter nuclei.
    ##
    ## The moments of the Fermi factor are analytic, M_n = -n! a^(n+1) Li_(n+1)(-exp(c/a)), so
    ##      <r^2> = [M_4 + (w/c^2) M_6] / [M_2 + (w/c^2) M_4]
    ## in closed form; Nuclear.threeParameterFermiRrms evaluates it through Math.polylogExp.  Reproducing
    ## the tabulated R_rms from the tabulated (c, a, w) is the one test of this model that uses measured
    ## numbers rather than JAC's own.
    ##
    ## CAVEAT ON THESE THREE ROWS (12-Aug-2026).  The (c, a, w) triples and the reference R_rms below were
    ## written down from recollection of the de Vries compilation, NOT read off the table itself.  The
    ## formula is separately verified to 1e-15 against direct numerical quadrature of the same density, so
    ## the 0.1-0.3 % seen here is a statement about the input numbers, not about the closed form.  Anyone
    ## using this as a literature validation should re-enter the triples from the published table first.
    println("\n  Tabulated three-parameter Fermi nuclei (de Vries et al. 1987):\n")
    println("     nucleus     c [fm]    a [fm]      w        R_rms (JAC)   R_rms (table)   difference")
    for (nam, c, a, w, Rtab) in [("40Ca", 3.7660, 0.5860, -0.1610, 3.478),
                                 ("48Ca", 3.7369, 0.5245, -0.0300, 3.474),
                                 ("58Ni", 4.3092, 0.5169, -0.1308, 3.776)]
        R = Nuclear.threeParameterFermiRrms(c, a, w)
        @printf("     %-8s  %7.4f   %7.4f   %8.4f    %8.5f      %8.3f      %+7.2f %%\n",
                nam, c, a, w, R, Rtab, 100*(R-Rtab)/Rtab)
    end

    ## w = 0 must recover the two-parameter model identically -- the regression test this model brings with it
    println("\n  w = 0 against the two-parameter closed form Nuclear.fermiRrms:\n")
    println("     R_rms [fm]     c [fm]        2-par R_rms          3-par R_rms (w=0)     rel. difference")
    for R in [3.0, 4.5, 5.5, 7.0]
        b  = Nuclear.computeFermiBParameter(R)
        R2 = Nuclear.fermiRrms(b)
        R3 = Nuclear.threeParameterFermiRrms(b, Nuclear.fermiA, 0.0)
        @printf("     %8.3f    %10.6f    %18.14f   %18.14f      %9.2e\n", R, b, R2, R3, abs(R3-R2)/R2)
    end


elseif  false
    # Last visit:  12-Aug-2026
    # Last successful:  12-Aug-2026
    #
    ## (b)  The shape effect on the 1s level of a hydrogen-like heavy and superheavy ion.
    ##
    ## Rule 12: the box is matched to the orbital.  A 1s electron at Z ~ 100 turns over at r ~ 1/Z a.u., so
    ## the 614 a.u. default grid would starve the B-spline basis completely; rbox = 25/Z is used instead.
    println("\n  1s_1/2 of a hydrogen-like ion, at FIXED rms radius, for several shape parameters w:\n")

    function level1s(Z::Float64, model, Rrms::Float64, prim)
        nm  = Nuclear.Model(Z, model, 2.5*Z, Rrms, AngularJ64(0), 0., 0., 0.)
        orb = Bsplines.generateOrbitalsHydrogenic([Subshell("1s_1/2")], nm, prim; printout=false)
        return( orb[Subshell("1s_1/2")].energy )
    end
    toCm(x) = Defaults.convertUnits("energy: from atomic to Kayser", x)

    for Z in [92.0, 120.0]
        grid = Radial.Grid(Radial.Grid(false), rnt=1.0e-8, h=3.0e-2, hp=5.0e-3, rbox=25.0/Z)
        prim = Bsplines.generatePrimitives(grid)
        Rrms = Nuclear.rrmsRadius(2.5*Z)
        @printf("\n     Z = %5.1f,  R_rms = %6.3f fm,  box = %8.4f a.u.,  %d grid points\n",
                Z, Rrms, 25.0/Z, grid.NoPoints)
        ePoint = level1s(Z, PointNucleus(), Rrms, prim)
        eFermi = level1s(Z, FermiNucleus(), Rrms, prim)
        @printf("       point nucleus            1s = %16.8f a.u.\n", ePoint)
        @printf("       2-par Fermi  (w =  0.00)  1s = %16.8f a.u.   finite-size shift = %13.2f cm^-1\n",
                eFermi, toCm(eFermi-ePoint))
        for w in [-0.30, -0.20, -0.10, 0.20]
            e3 = level1s(Z, ThreeParameterFermiNucleus(w), Rrms, prim)
            @printf("       3-par Fermi  (w = %5.2f)  1s = %16.8f a.u.   shape effect      = %13.2f cm^-1  (%5.2f %% of it)\n",
                    w, e3, toCm(e3-eFermi), 100*(e3-eFermi)/(eFermi-ePoint))
        end
    end

    ## WHAT THE NUMBERS SAY (12-Aug-2026).  For U^91+ the finite-size shift is 1.537e6 cm^-1 = 190.5 eV, in
    ## the right place for the ~198 eV of the literature once the slightly different R_rms is allowed for.
    ## The shape effect at fixed R_rms is -288 cm^-1 for w = -0.2, i.e. 0.02 % of that, and it is LINEAR in w
    ## (-144, -288, -422 cm^-1 for w = -0.1, -0.2, -0.3) exactly as a term first order in <r^4> should be.
    ## At Z = 120 the same shape effect is -8921 cm^-1 = -1.11 eV, 0.025 % of a much larger finite-size
    ## shift: the shape term grows FASTER with Z than the leading one, which is again what a higher radial
    ## moment does.  A central depression (w < 0) makes the level MORE bound, since at fixed <r^2> it moves
    ## charge into the surface and so raises <r^4>.
    ##
    ## For an electron this is a fraction of a wavenumber to an eV -- real but small.  The same quantity is
    ## what muonic-atom spectroscopy measures directly, because the muon orbits largely INSIDE the nucleus
    ## and samples the shape rather than only the second moment.


elseif  false
    # Last visit:  12-Aug-2026
    # Last successful:  12-Aug-2026
    #
    ## (c)  Where the three-parameter form stops being meaningful, and how it says so.
    ##
    ## Two limits are worth knowing about, and both raise rather than return a quietly wrong number.
    println("\n  Limits of the three-parameter Fermi form:\n")

    ## (i) The weight 1 + w r^2/c^2 must not compete with the Fermi factor.  It does once |w| a^2/c^2 grows
    ##     to O(1), i.e. once c falls towards a, and there the second moment ceases to be positive.
    ##
    ##     READ THE c = 0.50 fm LINE CAREFULLY.  It does NOT raise -- it returns 3.30 fm, a positive number
    ##     for a nucleus whose half-density radius is half a femtometre.  The guard in
    ##     Nuclear.threeParameterFermiRrms catches only the case where the ratio of moments turns negative,
    ##     which is a sufficient condition for nonsense and not a necessary one.  This whole region is
    ##     outside the model's domain, which is why Nuclear.computeThreeParameterFermiC never enters it:
    ##     it brackets above c = 4a and reports an unreachable radius instead of walking into this.
    for (c, w) in [(5.5, -0.20), (2.5, -0.20), (1.0, -0.20), (0.5, -0.20)]
        wa = try  @sprintf("%10.5f fm", Nuclear.threeParameterFermiRrms(c, Nuclear.fermiA, w))
             catch e   "raises (the weight outweighs the Fermi factor)"   end
        @printf("     c = %5.2f fm, a = %6.4f fm, w = %5.2f   ->   R_rms = %s\n", c, Nuclear.fermiA, w, wa)
    end

    ## (ii) A requested rms radius that no (c, a, w) can produce is reported, not silently approximated.
    println()
    for R in [5.5, 3.0, 1.5]
        wa = try  @sprintf("c = %10.6f fm", Nuclear.computeThreeParameterFermiC(R, Nuclear.fermiA, -0.1))
             catch e   "raises (radius below the representable floor)"   end
        @printf("     R_rms = %5.2f fm, w = -0.10   ->   %s\n", R, wa)
    end
end
