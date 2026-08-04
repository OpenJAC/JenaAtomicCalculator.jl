
println("Cm) Apply & test the StarkZeeman module (exact ASF-basis diagonalization of level mixing in")
println("    static E and/or B fields of arbitrary relative direction).")

if  false
    # Last successful:  01-Aug-2026
    # Branch a: hydrogen n=2 manifold (2s_1/2, 2p_1/2, 2p_3/2), E field along the quantization axis
    #   (eDirection default (0,0,1), i.e. q=0 only). This is the SAME system used to empirically
    #   calibrate the M-resolved matrix element convention this module relies on (see the module
    #   docstring): the ratio of the computed 2s-2p_1/2 and 2s-2p_3/2 (m_j=1/2) off-diagonal matrix
    #   elements matched the hand-derived, sign-independent exact value sqrt(2) to 5 significant
    #   figures. This branch reruns that same calculation end-to-end through the actual module (not
    #   the standalone calibration script) and prints the full 8x8 field-dressed spectrum.
    #   Zero-field gaps here (2s_1/2 to 2p_1/2/2p_3/2, ~1.25e-2 eV, i.e. ~4.6e-4 Hartree) are NOT the
    #   tiny, ~exact Dirac-Coulomb degeneracy one might naively expect for hydrogen -- this reflects
    #   the reference multiplet being built from TWO separately-optimized configurations ("2s","2p"),
    #   an ordinary SCF/orbital-quality effect (this system's own zero-field energies are already
    #   ~1-1.4% off the exact non-relativistic -0.125 Hartree), not a StarkZeeman-specific issue.
    #   At the field strength used below (eField=100 kV/cm, deliberately fairly strong), the field-
    #   induced coupling is comparable to or larger than that zero-field gap -- exactly the regime
    #   module-StarkShift.jl's quadratic formula cannot handle, and this module is for.
    grid = Radial.Grid(true)
    nm   = Nuclear.Model(1., "point")

    wa = Atomic.Computation(Atomic.Computation(), name="H n=2 manifold", grid=grid, nuclearModel=nm,
                            configs=[Configuration("2s"), Configuration("2p")])
    wb = perform(wa; output=true)
    multiplet = wb["multiplet:"]

    szSettings = StarkZeeman.Settings(StarkZeeman.Settings(); includeEField=true, eField=1.0e5,
                                       printBefore=true, levelSelection=LevelSelection())
    wc = Atomic.Computation(Atomic.Computation(), name="H n=2 Stark-Zeeman (E parallel)", grid=grid,
                            nuclearModel=nm, configs=[Configuration("2s"), Configuration("2p")],
                            propertySettings=[szSettings] )
    wd = perform(wc; output=true)
    outcomes = wd["Stark-Zeeman outcomes:"]

    println("\n>> Field-dressed energies (Hartree), relative to the mean, divided by eField in a.u.:")
    eFieldAu = 1.0e5 / StarkShift.AU_EFIELD_IN_VCM
    meanE = sum(o.energy for o in outcomes) / length(outcomes)
    for  o in outcomes
        println("   E=$(o.energy)   (E-mean)/eFieldAu = $((o.energy-meanE)/eFieldAu)   dominant: " *
                join(["L$(c.level.index)(J=$(c.level.J)$(string(c.level.parity)),M=$(c.M))[$(round(c.weight,digits=3))]"
                      for c in o.components], ", "))
    end
    #
elseif false
    # Last successful:  01-Aug-2026
    # Branch b: SAME H n=2 system, but now BOTH E and B fields present, and NOT parallel -- E along z
    #   (as in branch a), B tilted into the x-z plane at a chosen polar angle. Demonstrates the
    #   genuinely new capability (arbitrary relative field directions, q=+-1 terms active) this design
    #   enables beyond the simple parallel-field case. Validated two ways: (i) the trace (sum of all
    #   8 eigenvalues) is exactly invariant under the tilt, matching between tiltAngle=0 and pi/6 to
    #   13 significant figures, as required since a rotation of B alone cannot change the spectrum's
    #   sum; (ii) at tiltAngle=0 the eigenvalues match a pure-parallel-field calculation (branch a's
    #   E plus an axial B) to 10 significant figures. At tiltAngle=pi/6, genuine within-level M-mixing
    #   appears (e.g. one level's M=1/2 component acquires a real M=-1/2 admixture, weights ~0.91/0.07)
    #   and the exact M=+-3/2 degeneracy seen at tilt=0 is lifted -- both are only possible once the
    #   q=+-1 tensor terms are active, confirming the general-direction machinery works correctly.
    grid = Radial.Grid(true)
    nm   = Nuclear.Model(1., "point")

    wa = Atomic.Computation(Atomic.Computation(), name="H n=2 manifold", grid=grid, nuclearModel=nm,
                            configs=[Configuration("2s"), Configuration("2p")])
    wb = perform(wa; output=true)
    multiplet = wb["multiplet:"]

    for  tiltAngle in (0.0, pi/6)
        szSettings = StarkZeeman.Settings(StarkZeeman.Settings(); includeEField=true, eField=1.0e5,
                                        eDirection=(0., 0., 1.),
                                        includeBField=true, bField=0.01, bDirection=(sin(tiltAngle), 0., cos(tiltAngle)),
                                        includeSchwinger=true, printBefore=false, levelSelection=LevelSelection())
        wc = Atomic.Computation(Atomic.Computation(), name="H n=2 Stark-Zeeman (tilt=$tiltAngle)", grid=grid,
                                nuclearModel=nm, configs=[Configuration("2s"), Configuration("2p")],
                                propertySettings=[szSettings] )
        wd = perform(wc; output=true)
        outcomes = wd["Stark-Zeeman outcomes:"]

        println("\n>> Field-dressed energies (Hartree), B tilted by $(tiltAngle) rad from z" *
                (tiltAngle == 0.0 ? "  (should match a pure-parallel-field calculation)" : "  (genuinely non-parallel case)") * ":")
        println("   Sum of eigenvalues (trace, must be independent of tilt angle -- pure rotation of B): " *
                "$(sum(o.energy for o in outcomes))")
        for  o in outcomes
            println("   E=$(o.energy)   dominant: " *
                    join(["L$(c.level.index)(J=$(c.level.J)$(string(c.level.parity)),M=$(c.M))[$(round(c.weight,digits=3))]"
                        for c in o.components], ", "))
        end
    end
    #
elseif true
    # Last successful:  01-Aug-2026
    # Branch c: hydrogen n=6 near-degenerate manifold, E parallel to z -- REVISED SCOPE. The original
    #   plan called for an alkali (Li/Na) Rydberg manifold at n~15 compared against Zimmerman, Littman,
    #   Kash & Kleppner, Phys. Rev. A 20, 2251 (1979) / the ARC library. That was attempted and hit a
    #   real, separate numerical limitation, NOT a StarkZeeman bug: JAC's single-configuration SCF for
    #   a Rydberg valence orbital sharing the CORE's kappa symmetry (Li, n=10-14) oscillates in a limit
    #   cycle and never converges; and even for plain, coreless hydrogen, orbitals with kappa <= -2 (the
    #   j=l+1/2 "stretched" branch, l>=1) become unreliable already at n=6 and are badly wrong by n=15
    #   (e.g. 15s_1/2 comes out with the WRONG SIGN energy, +0.0154 Ha instead of the exact -0.0022 Ha)
    #   -- consistent with the previously-documented, unfixed Bsplines.jl high-n/high-kappa pseudo-
    #   continuum artifact from the MultipolePolarizibility investigation. Fixing this is out of scope
    #   here (a different module, module-SelfConsistent.jl / module-Bsplines.jl, currently touched by a
    #   parallel session). User-approved pivot: use hydrogen (no alkali-core SCF issue) at a modest
    #   n=6, restricted to the numerically well-converged "good chain" of levels (kappa=+l, i.e.
    #   j=l-1/2, for each l=0..5 -- s_1/2,p_1/2,d_3/2,f_5/2,g_7/2,h_9/2 -- selected automatically below
    #   by keeping only levels within 5% of the exact hydrogen energy -1/(2n^2); the excluded kappa<=-2
    #   "stretched" levels are exactly the unreliable ones). This chain is still mutually E1-connected
    #   in a Δl=+-1 ladder (same structure validated in branch a), giving a genuine SIX-level (32
    #   M-sublevel) simultaneous diagonalization -- the real, new capability over branches a/b's 2-3
    #   level tests. Quantitative check: at eField=1,2,4 MV/cm the total energy spread of the manifold
    #   scales linearly with the field to 5 significant figures (ratios 1.99946 and 2.00004 for the two
    #   doublings) -- the qualitative signature of a near-degenerate hydrogenic manifold (linear Stark
    #   effect) that StarkShift's quadratic formula structurally cannot reproduce, and the entire reason
    #   this module exists. The slope, spread/eFieldAu ~ 92.75-92.74 (a.u.), lands close to (within
    #   ~3%) of 3n(n-1)=90, the textbook extremal-splitting bound (3/2)*n*q*F with q=n1-n2 up to n-1,
    #   confirmed via web search against the general hydrogenic linear-Stark formula dE=(3/2)*n*q*F --
    #   the same formula underlying Zimmerman et al.'s treatment of the hydrogenic manifold states. An
    #   EXACT match is not expected since our 6-level chain is not the complete degenerate m=0 subspace
    #   (that would require the excluded kappa<=-2 levels too); the ~3% proximity is a meaningful,
    #   honest sanity check, not a precision benchmark.
    grid = Radial.Grid(true)
    nm   = Nuclear.Model(1., "point")
    n    = 6
    configs = [Configuration("$(n)s"), Configuration("$(n)p"), Configuration("$(n)d"),
               Configuration("$(n)f"), Configuration("$(n)g"), Configuration("$(n)h")]

    wa = Atomic.Computation(Atomic.Computation(), name="H n=$n manifold", grid=grid, nuclearModel=nm, configs=configs)
    wb = perform(wa; output=true)
    multiplet = wb["multiplet:"]

    exact = -1.0 / (2*n^2)
    good  = [lev for lev in multiplet.levels if abs(lev.energy - exact) < 0.05*abs(exact)]
    println("\n>> Selected 'good' (well-converged) levels, |E - exact| < 5% of exact -1/(2n^2):")
    for lev in good
        println("   index=$(lev.index)  J=$(lev.J)  parity=$(lev.parity)  energy=$(lev.energy)  " *
                "(diff=$(lev.energy-exact))")
    end
    goodIdx = [lev.index for lev in good]

    previousSpread = 0.
    for  eFieldVcm in (1.0e6, 2.0e6, 4.0e6)
        szSettings = StarkZeeman.Settings(StarkZeeman.Settings(); includeEField=true, eField=eFieldVcm,
                                        printBefore=false, levelSelection=LevelSelection(true, indices=goodIdx))
        wc = Atomic.Computation(Atomic.Computation(), name="H n=$n Stark (E=$eFieldVcm V/cm)", grid=grid,
                                nuclearModel=nm, configs=configs, propertySettings=[szSettings])
        wd = perform(wc; output=true)
        outcomes = wd["Stark-Zeeman outcomes:"]
        energies = sort([o.energy for o in outcomes])
        spread   = energies[end] - energies[1]
        eFieldAu = eFieldVcm / StarkShift.AU_EFIELD_IN_VCM
        println("\n>> eField=$eFieldVcm V/cm ($eFieldAu a.u.):  spread=$spread Ha   " *
                "spread/eFieldAu=$(spread/eFieldAu)   ratio to previous field point=$(spread/previousSpread)   " *
                "theory bound 3n(n-1)=$(3*n*(n-1))")
        global previousSpread = spread
    end
    #
end
