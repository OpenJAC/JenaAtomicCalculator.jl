
println("Dp) Photo-ionization with subsequent fluorescence: strength, angular distribution and polarization.")

setDefaults("print summary: open", "zzz-PhotoIonizationFluores.sum")
setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, pure sine")
setDefaults("unit: energy", "eV")

# WRITTEN 19-Aug-2026, the first implementation of Cascade-free photo-ionization-fluorescence.
#
# THE PROCESS.  gamma(omega, P1 P2 P3) + A(J_0) --> A^+(J_m) + e_p --> A^+(J_f) + e_p + gamma'.  Light of
# well-defined Stokes parameters ionizes the atom; the residual ion is thereby left ALIGNED (and ORIENTED, if the
# light is circularly polarized), so its subsequent radiative decay is neither isotropic nor unpolarized.
#
# HOW THE MODULE IS BUILT, because it explains what these branches can and cannot check.  The two steps are
# computed INDEPENDENTLY -- PhotoIonization.computeLines for i --> m + e_p and PhotoEmission.computeLines for
# m --> f + gamma' -- and only then combined, so the partial-wave sum is carried once per (i,m) pair however many
# final levels follow.  The strength of a pathway is
#
#     S(i,m,f)  =  sigma_ion(i --> m) * A_r(m --> f) / Gamma_r(m) ,
#
# the exact analogue of the resonance strength of DielectronicRecombination.  The statistical tensors rho_kq of the
# intermediate ion come from PhotoIonization.computeStatisticalTensorUnpolarized, which already handles arbitrary
# (P1, P2, P3); the angular distribution and the Stokes parameters of the fluorescence are then read off its
# 2x2 helicity density matrix.
#
# WHAT MUST BE CHECKED BEFORE ANY NUMBER IS BELIEVED, and why branch a comes first: a level with J_m = 0 or 1/2
# CANNOT BE ALIGNED.  Its fluorescence must therefore come out isotropic and unpolarized whatever the incident
# light.  That is a property of angular momentum alone, independent of every amplitude in the calculation, so it
# tests the machinery without needing an external reference.

grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 20.0)

if  true
    # Last visit:       19-Aug-2026
    # Last successful:  unknown ... see the verdict written in after the first run.
    #
    # Branch a: THE NULL TEST.  Neon, 2s photo-ionization followed by 2p --> 2s fluorescence:
    #
    #     gamma + Ne(1s^2 2s^2 2p^6, J=0)  -->  Ne^+(1s^2 2s 2p^6, ^2S_1/2) + e_p  -->  Ne^+(1s^2 2s^2 2p^5, ^2P) + gamma'
    #
    #   The intermediate level is ^2S_1/2, so J_m = 1/2 and NO alignment is possible: rho_2q must vanish, A_20 must
    #   be zero, and the fluorescence must be isotropic and unpolarized for EVERY choice of incident Stokes
    #   parameters.  The branch therefore runs the same case under unpolarized, linearly and circularly polarized
    #   light and requires the same (null) answer from all three.  A non-zero anisotropy here would mean the
    #   machinery is manufacturing alignment out of nothing.
    #   The system is otherwise the one of example-Dc.jl branch a, which is a working photo-ionization case.
    angles = [ SolidAngle(0.0, 0.0), SolidAngle(pi/4, 0.0), SolidAngle(pi/2, 0.0) ]
    for  (label, stokes)  in  [ ("unpolarized", ExpStokes(0., 0., 0.)),
                                ("linear 0 deg", ExpStokes(1., 0., 0.)),
                                ("circular",     ExpStokes(0., 0., 1.)) ]
        println("\n\n*** Branch a -- null test, incident light: $label ***\n")
        pifSettings = PhotoIonizationFluores.Settings(PhotoIonizationFluores.Settings();
                                                      multipoles     = [E1],
                                                      gauges         = [UseCoulomb, UseBabushkin],
                                                      photonEnergies = [80.0],
                                                      lValues        = [0, 1, 2],
                                                      incidentStokes = stokes,
                                                      calcAngular    = true,
                                                      calcStokes     = true,
                                                      solidAngles    = angles )
        wa = Atomic.Computation(Atomic.Computation(), name="Dp-a: Ne 2s photoionization + fluorescence, $label",
                                grid=grid, nuclearModel=Nuclear.Model(10.),
                                initialConfigs      = [Configuration("1s^2 2s^2 2p^6")],
                                intermediateConfigs = [Configuration("1s^2 2s^1 2p^6")],
                                finalConfigs        = [Configuration("1s^2 2s^2 2p^5")],
                                processSettings     = pifSettings )
        perform(wa)
    end
    #
elseif  false
    # Last visit:       19-Aug-2026
    # Last successful:  unknown ...
    #
    # Branch b: THE STRENGTH, and the one internal identity it must satisfy.  Magnesium, 2p photo-ionization
    #   followed by 3s --> 2p fluorescence:
    #
    #     gamma + Mg(1s^2 2s^2 2p^6 3s^2, J=0) --> Mg^+(1s^2 2s^2 2p^5 3s^2, ^2P_1/2,3/2) + e_p
    #                                          --> Mg^+(1s^2 2s^2 2p^6 3s, ^2S_1/2) + e_p + gamma'
    #
    #   Here the intermediate ^2P_3/2 HAS J_m = 3/2 and can be aligned, which is what branches c and d need.
    #   The check available without any external reference is the branching identity: summed over ALL final levels
    #   of a given intermediate level m, the strengths must reproduce the photo-ionization cross section into m,
    #   since sum_f A_r(m-->f)/Gamma_r(m) = 1 by construction.  With a single final level, as here, that means the
    #   strength must EQUAL sigma_ion(i-->m) to machine precision -- a check on the assembly, not on the physics.
    pifSettings = PhotoIonizationFluores.Settings(PhotoIonizationFluores.Settings();
                                                  multipoles     = [E1],
                                                  gauges         = [UseCoulomb, UseBabushkin],
                                                  photonEnergies = [80.0],
                                                  lValues        = [0, 1, 2],
                                                  incidentStokes = ExpStokes(0., 0., 0.),
                                                  calcAngular    = false,
                                                  calcStokes     = false )
    wa = Atomic.Computation(Atomic.Computation(), name="Dp-b: Mg 2p photoionization + 3s-2p fluorescence",
                            grid=grid, nuclearModel=Nuclear.Model(12.),
                            initialConfigs      = [Configuration("1s^2 2s^2 2p^6 3s^2")],
                            intermediateConfigs = [Configuration("1s^2 2s^2 2p^5 3s^2")],
                            finalConfigs        = [Configuration("1s^2 2s^2 2p^6 3s")],
                            processSettings     = pifSettings )
    perform(wa)
    #
elseif  false
    # Last visit:       19-Aug-2026
    # Last successful:  unknown ...
    #
    # Branch c: THE ANGULAR DISTRIBUTION, and its dependence on the incident polarization.  The Mg case of branch b,
    #   run under unpolarized, linearly and circularly polarized light.  Two things must hold, and both follow from
    #   the structure of the photon statistical tensors rather than from any amplitude:
    #     + A_20 IS INDEPENDENT OF P3.  Circular polarization enters only through rho_1q, so it produces ORIENTATION
    #       and cannot change the ALIGNMENT.  The unpolarized and circular runs must give the same A_20.
    #     + LINEAR POLARIZATION CHANGES IT.  rho_2,+-2 is carried by P1 -+ i P2, so a linearly polarized beam adds a
    #       component that an unpolarized one does not have, and W(theta) must differ.
    #   An unpolarized beam nevertheless still defines an axis, so its A_20 is NOT zero -- alignment without
    #   handedness.  That is the point most easily got wrong.
    angles = [ SolidAngle(0.0, 0.0), SolidAngle(pi/4, 0.0), SolidAngle(pi/2, 0.0), SolidAngle(3pi/4, 0.0) ]
    for  (label, stokes)  in  [ ("unpolarized", ExpStokes(0., 0., 0.)),
                                ("linear 0 deg", ExpStokes(1., 0., 0.)),
                                ("circular",     ExpStokes(0., 0., 1.)) ]
        println("\n\n*** Branch c -- angular distribution, incident light: $label ***\n")
        pifSettings = PhotoIonizationFluores.Settings(PhotoIonizationFluores.Settings();
                                                      multipoles     = [E1],
                                                      gauges         = [UseCoulomb, UseBabushkin],
                                                      photonEnergies = [80.0],
                                                      lValues        = [0, 1, 2],
                                                      incidentStokes = stokes,
                                                      calcAngular    = true,
                                                      calcStokes     = false,
                                                      solidAngles    = angles )
        wa = Atomic.Computation(Atomic.Computation(), name="Dp-c: Mg angular distribution, $label",
                                grid=grid, nuclearModel=Nuclear.Model(12.),
                                initialConfigs      = [Configuration("1s^2 2s^2 2p^6 3s^2")],
                                intermediateConfigs = [Configuration("1s^2 2s^2 2p^5 3s^2")],
                                finalConfigs        = [Configuration("1s^2 2s^2 2p^6 3s")],
                                processSettings     = pifSettings )
        perform(wa)
    end
    #
elseif  false
    # Last visit:       19-Aug-2026
    # Last successful:  unknown ...
    #
    # Branch d: THE POLARIZATION of the fluorescence, which separates alignment from orientation more sharply than
    #   the angular distribution does.  Same Mg case, same three incident polarizations, now reading the Stokes
    #   parameters of the emitted line:
    #     + P1 of the fluorescence is fed by the ALIGNMENT, so it is non-zero already for UNPOLARIZED incident light
    #       and changes with linear polarization.
    #     + P3 of the fluorescence is fed by the ORIENTATION, i.e. by rho_1q, which exists ONLY for circularly
    #       polarized incidence.  It must therefore be zero for the unpolarized and linear runs and non-zero for the
    #       circular one.  That rule is exact and is the sharpest check in this file.
    angles = [ SolidAngle(pi/4, 0.0), SolidAngle(pi/2, 0.0) ]
    for  (label, stokes)  in  [ ("unpolarized", ExpStokes(0., 0., 0.)),
                                ("linear 0 deg", ExpStokes(1., 0., 0.)),
                                ("circular",     ExpStokes(0., 0., 1.)) ]
        println("\n\n*** Branch d -- fluorescence polarization, incident light: $label ***\n")
        pifSettings = PhotoIonizationFluores.Settings(PhotoIonizationFluores.Settings();
                                                      multipoles     = [E1],
                                                      gauges         = [UseCoulomb, UseBabushkin],
                                                      photonEnergies = [80.0],
                                                      lValues        = [0, 1, 2],
                                                      incidentStokes = stokes,
                                                      calcAngular    = true,
                                                      calcStokes     = true,
                                                      solidAngles    = angles )
        wa = Atomic.Computation(Atomic.Computation(), name="Dp-d: Mg fluorescence polarization, $label",
                                grid=grid, nuclearModel=Nuclear.Model(12.),
                                initialConfigs      = [Configuration("1s^2 2s^2 2p^6 3s^2")],
                                intermediateConfigs = [Configuration("1s^2 2s^2 2p^5 3s^2")],
                                finalConfigs        = [Configuration("1s^2 2s^2 2p^6 3s")],
                                processSettings     = pifSettings )
        perform(wa)
    end
    #
end

setDefaults("print summary: close", "")
