
println("Pc) RESONANT INELASTIC SCATTERING (RIXS): gamma + |i> --> |f> + gamma through a resonant intermediate level.")

setDefaults("print summary: open", "zzz-PhotonScattering.sum")
setDefaults("unit: energy", "eV")

if  true
    # Last visit:      21-Aug-2026
    # Last successful:  unknown ... the Lorentzian IS reproduced to 0.6 %, but the absolute cross sections rest on
    #                              an underived prefactor, Gamma was chosen rather than computed, and nothing here
    #                              constrains a PHASE.  See the report below.
    #
    # Branch a (the resonance profile): beryllium-like neon, scanning the incoming photon energy ACROSS the
    #   1s^2 2s^2 ^1S_0 --> 1s^2 2s 2p (1^-) resonance at 25.003 eV, on the resonant-ELASTIC channel.
    #
    #   WHY THIS IS THE RIGHT FIRST TEST, and why it needs no prefactor.  On resonance the amplitude is governed by
    #   the regularised denominator E_i + omega_in - E_nu + i*Gamma/2, so the cross section traces a LORENTZIAN in
    #   the detuning delta:
    #
    #       sigma(delta) / sigma(0)  =  (Gamma^2/4) / (delta^2 + Gamma^2/4)
    #
    #   That is a RATIO, so the underived constant in rayleighCrossSection cancels completely -- the same trick as
    #   the omega^4 law in example-Pb.jl and the Z^5 scan in example-Of.jl.  It yields one sharp, falsifiable
    #   number: the HALF MAXIMUM must fall at a detuning of exactly Gamma/2, i.e. the profile's FWHM must equal the
    #   width that was put in.  With Gamma = 0.1 eV the points at delta = +-0.05 eV must sit at half the peak.
    #
    #   THE WIDTH IS CHOSEN, NOT COMPUTED, and that is worth stating plainly.  Gamma = 0.1 eV is picked so the
    #   resonance is resolvable on a practical scan.  The PHYSICAL width of a 2s2p level in Be-like neon is set by
    #   its radiative lifetime and is many orders of magnitude smaller, which would demand a scan in micro-eV
    #   steps.  A physical Gamma would come from a PhotoEmission (or AutoIonization) lifetime calculation for the
    #   intermediate level and be fed in here; nothing in this branch should be read as a prediction of a real
    #   line shape.  What IS tested is that the code produces the Lorentzian it claims to.
    #
    #   WHAT MAKES THIS RIXS RATHER THAN RAYLEIGH.  PhotonScattering.rayleighAmplitude SKIPS any intermediate level
    #   whose denominator falls below selfTolerance, because the perturbative expression is undefined there.  That
    #   is exactly the level RIXS is tuned onto, so resonantAmplitude keeps it and lets the width regularise it.
    #   The two files therefore make OPPOSITE decisions about the same situation, deliberately.
    #
    #   The intermediate levels come from gMultiplet, as everywhere else in this module -- 1s^2 2s 2p supplies the
    #   two 1^- levels at 10.816 and 25.003 eV.  This implementation is LINE-shaped; the older ResonantInelastic
    #   module is pathway-shaped and takes its intermediates from Basics.perform's intermediateMultiplet instead.
    #
    #   TWO DESIGN ERRORS OF MINE, recorded because neither was a code fault and both are easy to repeat.
    #     (1) THE FIRST ATTEMPT SCANNED THE WRONG RESONANCE, the 1^- level at 10.816 eV, and showed NO resonance
    #         at all -- the cross section merely rose 37 % across the scan, which was almost entirely the
    #         omega_in*omega_out^3 prefactor.  That level is the INTERCOMBINATION line, amplitude 7.56e-07 against
    #         1.06e-03 for the 25.003 eV line.  The resonant numerator goes as the SQUARE of the vertex, so a weak
    #         intermediate is doubly penalised: the resonant term came to ~3e-10 while the strong level 14 eV OFF
    #         resonance contributed ~2e-06, beating it by four orders of magnitude.  The same smallness that gave
    #         this level the worst gauge agreement in example-Pb.jl is why it cannot carry a resonance -- one cause,
    #         two symptoms.
    #     (2) THE FIRST ATTEMPT ALSO MIXED ENERGY SCALES.  finalConfigs listed 1s^2 2s^2 AND 1s^2 2p^2 while
    #         initialConfigs listed only the first, so the two multiplets came from SEPARATE SCF runs and the same
    #         physical level had different absolute energies in each.  E_f - E_i was then wrong by 3.09 eV and the
    #         nominally ELASTIC line printed omega_out = 13.90 eV against omega_in = 10.82 eV.  Absolute energies
    #         from separate SCF runs are not on a common scale, so any DIFFERENCE taken across them is meaningless.
    #         Both configurations lists are now identical and omega_out = omega_in exactly, as the table shows.
    #
    #   REPORT (21-Aug-2026, after retargeting to 25.003 eV):
    #
    #     detuning [eV]   sigma(Coulomb)   obs/peak   Lorentzian   Babushkin obs/peak
    #        -0.400        9.254062e-17     0.0147      0.0154          0.0147
    #        -0.200        3.625535e-16     0.0574      0.0588          0.0574
    #        -0.050        3.138101e-15     0.4970      0.5000          0.4970
    #        +0.000        6.313628e-15     1.0000      1.0000          1.0000
    #        +0.050        3.175594e-15     0.5030      0.5000          0.5030
    #        +0.200        3.803564e-16     0.0602      0.0588          0.0602
    #        +0.400        1.018614e-16     0.0161      0.0154          0.0161
    #
    #   THE HALF MAXIMUM FALLS AT delta = +-Gamma/2 TO WITHIN 0.6 % -- 0.4970 and 0.5030 against 0.5000 -- so the
    #   profile's FWHM equals the width that was put in, which is what the branch set out to test.  Peak
    #   enhancement over the far wing is 68.2x against 65.0x predicted.
    #
    #   TWO THINGS TO READ CORRECTLY.  First, BOTH GAUGES GIVE IDENTICAL PROFILE RATIOS, and that is not a bonus:
    #   the profile is set by the DENOMINATOR, which carries no gauge, so this test validates the resonance
    #   machinery and says NOTHING about the amplitude's gauge behaviour.  For that, see example-Pb.jl.  Second,
    #   the small asymmetry -- the + side high by ~9 % -- is the omega_in*omega_out^3 prefactor, which rises by
    #   13.7 % across the scan; it is not a Fano asymmetry and should not be read as one.
    #
    #   NOT DATED.  The SHAPE is verified to sub-percent, but the ABSOLUTE cross sections rest on the same
    #   underived prefactor as everything else in this module, Gamma was chosen rather than computed, and nothing
    #   here constrains a PHASE.
    #
    grid       = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, rbox = 10.0)
    nm         = Nuclear.Model(10.)
    gMultiplet = SelfConsistent.performSCF([Configuration("1s^2 2s 2p")], nm, grid, AsfSettings())

    settings   = PhotonScattering.Settings(PhotonScattering.Settings();
                        process        = PhotonScattering.ResonantScattering(),
                        approximation  = PhotonScattering.SecondOrderGreen(),
                        # symmetric about the 25.003 eV resonance; +-0.05 eV are the half-maximum points for Gamma = 0.1
                        photonEnergies = [24.603, 24.803, 24.953, 25.003, 25.053, 25.203, 25.403],
                        multipoles     = [E1],
                        gMultiplet     = gMultiplet,
                        width          = Defaults.convertUnits("energy: to atomic", 0.1),
                        selfTolerance  = 1.0e-6,
                        printBefore    = true )

    wa = Atomic.Computation(Atomic.Computation(), name="RIXS on Be-like neon", grid=grid, nuclearModel=nm,
                            initialConfigs  = [Configuration("1s^2 2s^2")],
                            finalConfigs    = [Configuration("1s^2 2s^2")],
                            processSettings = settings )
    wb = perform(wa; output=true)
    #
elseif  false
    # Last visit:      21-Aug-2026
    # Last successful:  unknown ... not yet run; it is a guard test, not a physics branch.
    #
    # Branch b (the zero-width refusal): the same computation with width = 0, to check that the module REFUSES
    #   rather than quietly dividing by whatever small number the energy grid happens to supply.
    #
    #   WHY THIS DESERVES A BRANCH.  With Gamma = 0 the resonant denominator is E_i + omega_in - E_nu, and its
    #   value on or near resonance is then decided by how close a scan point lands to an intermediate level --
    #   an artefact of the chosen energy grid rather than a property of the atom.  The result would look like a
    #   number and be one, and would change completely if a scan point moved by a micro-eV.  A process that
    #   cannot be evaluated should say so, following the corePolarization.doApply pattern in
    #   module-PhotoEmission.jl: an error that explains, rather than a wrong result or a crash at a random
    #   undefined name.
    #
    #   The expected behaviour is an error naming the width, saying why it is needed, and pointing at
    #   RayleighScattering() as the thing to use if the NON-resonant amplitude was what was wanted.
    #
    grid       = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, rbox = 10.0)
    nm         = Nuclear.Model(10.)
    gMultiplet = SelfConsistent.performSCF([Configuration("1s^2 2s 2p")], nm, grid, AsfSettings())

    settings   = PhotonScattering.Settings(PhotonScattering.Settings();
                        process        = PhotonScattering.ResonantScattering(),
                        photonEnergies = [10.816], multipoles = [E1],
                        gMultiplet     = gMultiplet, width = 0., printBefore = true )

    wa = Atomic.Computation(Atomic.Computation(), name="RIXS with a zero width", grid=grid, nuclearModel=nm,
                            initialConfigs  = [Configuration("1s^2 2s^2")],
                            finalConfigs    = [Configuration("1s^2 2s^2")],
                            processSettings = settings )
    wb = perform(wa; output=true)
    #
end
#
setDefaults("print summary: close", "")
