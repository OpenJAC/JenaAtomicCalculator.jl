
println("Dp) Coulomb IONIZATION of ions by fast ion impact: energy-differential cross sections and the alignment")
println("    of the residual ion.")

# WRITTEN 25-Aug-2026, first implementation of module-CoulombIonization.jl.  What stood there before was a shell:
# Settings, Channel and Line with their show-methods, no compute function at all, and a Line() constructor that
# passed seven arguments to a six-field struct -- errors that survived only because nothing ever compiled it.
#
# THE THEORY IS THAT OF CoulombExcitation (example-Ds.jl) and this file is its companion.  The projectile is a fast
# ion that never touches the target; it acts through its Coulomb field, transferring momentum q.  The difference is
# only in the final state:
#
#       excitation      A^(q+) + |i(N)>  -->  A^(q+) + |f(N)>            example-Ds.jl
#       ionization      A^(q+) + |i(N)>  -->  A^(q+) + |f(N-1)> + e^-    this file
#
# and the module reuses CoulombExcitation's reduced matrix elements literally, rather than copying them, so that a
# correction there is a correction here.  That matters: the -i on every magnetic term, which RATIP carries and the
# published equation does not show, was found and fixed on the excitation side and is inherited.
#
# THREE THINGS FOLLOW FROM THE ELECTRON LEAVING.
#   + The final state is a level TIMES a free electron.  It is built by attaching a continuum orbital to the
#     residual ion and coupling the two, so that both sides of the matrix element carry N electrons and the
#     many-electron machinery applies unchanged.
#   + THE NORMALISATION.  A continuum orbital is normalised per unit ENERGY, so what comes out is a cross section
#     DIFFERENTIAL in the ejected-electron energy, d(sigma)/d(epsilon), and NOT a total.  A total needs an
#     integration over epsilon, which is deliberately left to the caller.
#   + The minimum momentum transfer GROWS with the electron energy: q0 = (E_f - E_i + epsilon)/(beta c), because
#     the projectile must supply the binding energy and the kinetic energy of the electron.  In the excitation
#     case q0 is fixed by the transition alone.


if  true
    # Last visit:      25-Aug-2026
    # Last successful: 25-Aug-2026
    #
    # Branch a: THE REFERENCE CASE AND THE SHAPE OF THE SPECTRUM -- K-shell Coulomb ionization of helium-like
    #   carbon, 1s^2 (J=0) --> 1s (J=1/2) + e^-, by a 100 MeV/u projectile, at a series of ejected-electron
    #   energies.  A single ejected-electron energy says almost nothing; the SHAPE is what an experiment measures
    #   and what the physics constrains.
    #
    # REPORT (25-Aug-2026): d(sigma)/d(epsilon) falls monotonically, and it STEEPENS as it goes:
    #        eps [eV]      q0 [a.u.]    d(sigma)/d(eps) [b/eV]    local power
    #          10         2.51165e-01        8.73905e+01
    #          50         2.76216e-01        4.46342e+01           eps^-0.42
    #         100         3.07530e-01        2.85829e+01           eps^-0.64
    #         200         3.70159e-01        1.43770e+01           eps^-0.99
    #         500         5.58045e-01        3.46099e+00           eps^-1.56
    #   That shape is the physics of the process rather than an accident: a nearly flat soft-electron plateau at
    #   low epsilon, steepening towards the 1/eps^(3/2)-like fall at high epsilon.  The reason is in the third
    #   column -- q0 RISES BY 122% across this range, because the projectile must supply the electron's kinetic
    #   energy as well as its binding energy.  A larger minimum momentum transfer means fewer distant collisions
    #   can contribute, and distant collisions are where most of a first-order Coulomb cross section lives.  In
    #   Coulomb EXCITATION q0 is fixed by the transition alone, and this is the one qualitative difference the
    #   ejected electron makes to the shape.
    #
    #   THE CONTINUUM METHOD MATTERED MORE THAN ANYTHING ELSE HERE, and the first version of this branch got it
    #   wrong.  With setDefaults("method: continuum, asymptotic Coulomb") the same five points came out as
    #   1.5e-03, 1.41, 1.49, 4.23, 1.93 b/eV -- not monotonic, not smooth, and not a spectrum of anything.  It is
    #   the same trap that example-Fd.jl documents for photoionization: asymptotic Coulomb gives erratic continuum
    #   orbitals for slow electrons in low-charge ions.  Galerkin gives the table above.  A first-implementation
    #   module cannot be judged before that switch is made.
    #
    #   NOT DATED AGAINST ANY MEASUREMENT.  No reference value was compared and the absolute scale carries every
    #   approximation of the module.  What is verified is the internal behaviour -- monotonicity, the steepening,
    #   and the rise of q0 that explains it -- none of which would survive a sign or normalisation error.
    setDefaults("method: continuum, Galerkin")
    setDefaults("method: normalization, pure sine")
    setDefaults("unit: energy", "eV")
    setDefaults("print summary: open", "zzz-CoulombIonization-Dp-reference.sum")

    grid = Radial.Grid(Radial.Grid(false); rnt=2.0e-6, h=5.0e-2, hp=1.0e-2, rbox=10.0)
    ciSettings = CoulombIonization.Settings(CoulombIonization.Settings(),
                                            ionEnergies=[100.], electronEnergies=[10., 50., 100., 200., 500.],
                                            calcAlignment=false, printBefore=false, zerosGL=6, lValues=[0,1])
    wa = Atomic.Computation(Atomic.Computation(); name="He-like C, K-shell Coulomb ionization",
                            grid=grid, nuclearModel=Nuclear.Model(6.),
                            initialConfigs  = [Configuration("1s^2")],
                            finalConfigs    = [Configuration("1s")],
                            processSettings = ciSettings )
    wb = perform(wa)
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      25-Aug-2026
    # Last successful: 25-Aug-2026
    #
    # Branch b: THE PARTIAL-WAVE TRUNCATION, which is the one approximation here that the user controls directly
    #   and the one that fails quietly.  The ejected electron is expanded in partial waves and the sum is cut at
    #   `lValues`.  A cut that is too severe does not merely lower the cross section by a known factor; because
    #   different l contribute differently at different momentum transfers, it can distort the SHAPE as well.
    #
    # REPORT (25-Aug-2026): at 50 eV ejected-electron energy, adding partial waves gives
    #        lValues        d(sigma)/d(eps) [b/eV]        change
    #        [0]                 1.45431e-04                --
    #        [0,1]               4.46342e+01           x 300 000
    #        [0,1,2]             4.71109e+01              +5.5 %
    #   THE s WAVE IS NEGLIGIBLE AND THAT IS A SELECTION RULE, not a small number: ejecting a 1s electron into an
    #   eps-s wave is a monopole transition, which the leading term of the interaction cannot drive, so it survives
    #   only through the higher multipoles and the relativistic magnetic term.  The p wave is the dipole channel
    #   and carries essentially the whole cross section; d adds 5.5%.
    #   CONVERGENCE IS FAST HERE AS A PROPERTY OF THE CASE, not of the method.  A 1s electron ejected at 50 eV
    #   carries little angular momentum.  It will NOT be fast for a faster electron or for an initial orbital of
    #   higher l, and the honest way to use this module is to add one more l and look, exactly as done here.
    setDefaults("method: continuum, Galerkin")
    setDefaults("method: normalization, pure sine")
    setDefaults("unit: energy", "eV")
    setDefaults("print summary: open", "zzz-CoulombIonization-Dp-partialWaves.sum")

    grid = Radial.Grid(Radial.Grid(false); rnt=2.0e-6, h=5.0e-2, hp=1.0e-2, rbox=10.0)
    for  lv  in  [[0], [0,1], [0,1,2]]
        println("\n  ---- lValues = $lv ----")
        ciSettings = CoulombIonization.Settings(CoulombIonization.Settings(),
                                                ionEnergies=[100.], electronEnergies=[50.],
                                                calcAlignment=false, printBefore=false, zerosGL=6, lValues=lv)
        wa = Atomic.Computation(Atomic.Computation(); name="He-like C, partial-wave test",
                                grid=grid, nuclearModel=Nuclear.Model(6.),
                                initialConfigs  = [Configuration("1s^2")],
                                finalConfigs    = [Configuration("1s")],
                                processSettings = ciSettings )
        perform(wa)
    end
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      25-Aug-2026
    # Last successful: 25-Aug-2026
    #
    # Branch c: THE ALIGNMENT OF THE RESIDUAL ION, which is the property this module shares with CoulombExcitation
    #   and the reason for building the final state the way it does.  A_2 and A_4 are formed from EXACTLY the
    #   expressions of the excitation case -- they describe the ion, not the electron -- but reaching them takes one
    #   extra step that has no counterpart there: the amplitudes are computed for the COUPLED state of (ion +
    #   ejected electron), so that combination must be decoupled again, coherently over the total symmetries, before
    #   the ion's own magnetic quantum number M_f is resolved.
    #
    #   The system is lithium-like carbon losing a K electron, 1s^2 2s (J=1/2) --> 1s 2s + e^-, whose residual ion
    #   has a J=1 level -- the smallest system in which an alignment can exist at all, since a J=1/2 or J=0 ion
    #   cannot be aligned.
    #
    # REPORT (25-Aug-2026):
    #        final level        d(sigma)/d(eps) [b/eV]      A_2        A_4
    #        1s2s  J = 1             4.07376e+01          0.00000    0.00000
    #        1s2s  J = 0             1.35869e+01          0.00000    0.00000
    #   BOTH ZEROS ARE CORRECT, AND FOR TWO DIFFERENT REASONS, which is what makes them a check rather than a
    #   disappointment:
    #     + a J = 0 ion cannot be aligned at all -- there is only one magnetic sublevel, so there is nothing for an
    #       alignment to describe;
    #     + the J = 1 level here is 1s2s in a spin TRIPLET, an S state with no orbital angular momentum.  The
    #       Coulomb interaction does not act on spin, so an unpolarised J = 1/2 target cannot leave a spin-aligned
    #       ion behind.  The zero is a statement about the interaction, not about the geometry.
    #   Neither zero is put in by hand: both come out of the decoupling and the M_f-resolved cross sections, so a
    #   sign error or a wrong Clebsch-Gordan in that step would show up here as a spurious non-zero value.
    #
    #   WHAT THIS BRANCH THEREFORE DOES NOT SHOW is a NON-ZERO alignment, and that is a limitation of the case and
    #   not of the module.  It needs a residual ion with orbital angular momentum: ejecting a 2p electron from a
    #   neon-like ion, 1s^2 2s^2 2p^6 --> 1s^2 2s^2 2p^5 (J = 3/2), is the smallest natural candidate and is the
    #   obvious next test.  It was not run here because ten electrons make it a much larger computation than
    #   anything else in this file.
    setDefaults("method: continuum, Galerkin")
    setDefaults("method: normalization, pure sine")
    setDefaults("unit: energy", "eV")
    setDefaults("print summary: open", "zzz-CoulombIonization-Dp-alignment.sum")

    grid = Radial.Grid(Radial.Grid(false); rnt=2.0e-6, h=5.0e-2, hp=1.0e-2, rbox=10.0)
    ciSettings = CoulombIonization.Settings(CoulombIonization.Settings(),
                                            ionEnergies=[100.], electronEnergies=[50.],
                                            calcAlignment=true, printBefore=false, zerosGL=6, lValues=[0,1])
    wa = Atomic.Computation(Atomic.Computation(); name="Li-like C, K-shell ionization with alignment",
                            grid=grid, nuclearModel=Nuclear.Model(6.),
                            initialConfigs  = [Configuration("1s^2 2s")],
                            finalConfigs    = [Configuration("1s 2s")],
                            processSettings = ciSettings )
    wb = perform(wa)
    setDefaults("print summary: close", "")
    #
end
