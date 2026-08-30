
println("Ci) Apply & test the AlphaVariation module: q- and Q-sensitivity coefficients for the ground 2s^2 2p^2")
println("    fine-structure multiplet (3P0, 3P1, 3P2, 1D2, 1S0) of Carbon-like ions, by the finite-difference method.")
println("    Reference Q-values (Kozlov, Tupitsyn & Reimers, arXiv:0812.3210, Table I -- relative to the 3P0 ground level):")
println("       C I  (Z=6):   3P1  Q=1.0086    3P2  Q=1.0004    1D2  Q=0.0052    1S0  Q=0.0031")
println("       O III(Z=8):   3P1  Q=1.0197    3P2  Q=1.0040    1D2  Q=0.0099    1S0  Q=0.0058")
println("       Si IX(Z=14):  3P1  Q=1.1403    3P2  Q=1.0254    1D2  Q=0.0990    1S0  Q=0.0573")
println("    Nd13+/Sm15+: single valence electron above a Pd-like ([Kr] 4d^10) closed core, with a near 5s-4f level")
println("    crossing that gives a huge, opposite-sign q for the two ions (Berengut, Dzuba & Flambaum, arXiv:1208.4157):")
println("       Nd13+(Z=60):  ground 5s_1/2,   4f_5/2 at  58897 cm^-1,  q(5s-4f5/2)  = +106000 cm^-1")
println("       Sm15+(Z=62):  ground 4f_5/2,   5s_1/2 at  55675 cm^-1,  q(4f5/2-5s)  = -136000 cm^-1")

if true
    # Last successful:  21-Jul-2026 -- 3P1 Q=1.0035, 3P2 Q=1.0003, 1D2 Q=0.0033, 1S0 Q=0.0012 (vs literature above)
    # C I -- the lightest, simplest member of the sequence; small but nonzero relativistic sensitivity.
    setDefaults("print summary: open", "zzz-AlphaVariation.sum")
    wa = Atomic.Computation(Atomic.Computation(), name="CI-alpha", grid=Radial.Grid(true), nuclearModel=Nuclear.Model(6.),
                            configs=[Configuration("1s^2 2s^2 2p^2")],
                            propertySettings=[AlphaVariation.Settings(true, 0.125, true, LevelSelection())] )
    wb = perform(wa)
    setDefaults("print summary: close", "")
    #
elseif true
    # Last successful:  21-Jul-2026 -- 3P1 Q=1.0174, 3P2 Q=1.0028, 1D2 Q=0.0115, 1S0 Q=0.0046 (vs literature above)
    # O III -- doubly-ionized, moderate Z; noticeably larger Q than C I for the same 2s^2 2p^2 configuration.
    setDefaults("print summary: open", "zzz-AlphaVariation.sum")
    wa = Atomic.Computation(Atomic.Computation(), name="OIII-alpha", grid=Radial.Grid(true), nuclearModel=Nuclear.Model(8.),
                            configs=[Configuration("1s^2 2s^2 2p^2")],
                            propertySettings=[AlphaVariation.Settings(true, 0.125, true, LevelSelection())] )
    wb = perform(wa)
    setDefaults("print summary: close", "")
    #
elseif true
    # Last successful:  21-Jul-2026 -- 3P1 Q=1.1240, 3P2 Q=1.0166, 1D2 Q=0.1039, 1S0 Q=0.0454 (vs literature above)
    # Si IX -- highest Z of the three, clearly enhanced relativistic sensitivity; still a single configuration.
    setDefaults("print summary: open", "zzz-AlphaVariation.sum")
    wa = Atomic.Computation(Atomic.Computation(), name="SiIX-alpha", grid=Radial.Grid(true), nuclearModel=Nuclear.Model(14.),
                            configs=[Configuration("1s^2 2s^2 2p^2")],
                            propertySettings=[AlphaVariation.Settings(true, 0.125, true, LevelSelection())] )
    wb = perform(wa)
    setDefaults("print summary: close", "")
    #
elseif true
    # Last successful:  21-Jul-2026 -- q(5s-4f5/2) = 105889 cm^-1 (0.10% from literature 106000 cm^-1); excitation
    # 6.256 eV vs literature 7.302 eV (~14% low, expected for a bare single-configuration calculation without
    # core-valence correlation -- q itself is far more robust, being dominated by near-nucleus one-electron
    # relativistic contraction of the s vs f orbital).
    # Nd13+ -- single valence electron (5s or 4f) above a Pd-like [Kr]4d^10 core; an open-f-shell/near-crossing case.
    setDefaults("print summary: open", "zzz-AlphaVariation.sum")
    wa = Atomic.Computation(Atomic.Computation(), name="Nd13-alpha", grid=Radial.Grid(true), nuclearModel=Nuclear.Model(60.),
                            configs=[Configuration("[Kr] 4d^10 5s"), Configuration("[Kr] 4d^10 4f")],
                            propertySettings=[AlphaVariation.Settings(true, 0.125, true, LevelSelection())] )
    wb = perform(wa)
    setDefaults("print summary: close", "")
    #
elseif true
    # Last successful:  21-Jul-2026 -- q(4f5/2-5s) = -136121 cm^-1 (0.09% from literature -136000 cm^-1); excitation
    # 8.235 eV vs literature 6.903 eV (~19% high, same expected cause as Nd13+ above). Ground state comes out as
    # 4f_5/2 (not 5s_1/2), correctly reproducing the level crossing between the Nd13+ and Sm15+ members of this
    # isoelectronic sequence -- and the SIGN of q flips accordingly, exactly as in the literature.
    # Sm15+ -- same recipe as Nd13+, two charge states further along the sequence, past the 5s-4f crossing.
    setDefaults("print summary: open", "zzz-AlphaVariation.sum")
    # Bsplines.checkOrbitalConsistency refuses Radial.Grid(true) here: its two spin-orbit partners come out
    # describing DIFFERENT states, which Rule 12 names as the signature of a box not matched to the orbitals.
    # Radial.Grid(true) reaches 614 a.u.; Basics.recommendedGrid sizes Sm15+ at 7.4 a.u., an 83-fold oversize.
    grid = Basics.recommendedGrid([Configuration("[Kr] 4d^10 5s"), Configuration("[Kr] 4d^10 4f")],
                                  Nuclear.Model(62.), printout=false)
    wa = Atomic.Computation(Atomic.Computation(), name="Sm15-alpha", grid=grid, nuclearModel=Nuclear.Model(62.),
                            configs=[Configuration("[Kr] 4d^10 5s"), Configuration("[Kr] 4d^10 4f")],
                            propertySettings=[AlphaVariation.Settings(true, 0.125, true, LevelSelection())] )
    wb = perform(wa)
    setDefaults("print summary: close", "")
    #
end
