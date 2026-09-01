#
println("Ei) Two-photon IONIZATION of a one-electron system: generalized cross sections in GM, the convergence of")
println("    the intermediate sum, and the resonance guard -- all driven through Atomic.Computation.")
#
# MOVED HERE FROM example-Dv.jl ON 01-Sep-2026, on the maintainer's instruction, and the move is the point rather
# than a tidying.  Two-photon ionization needs an INTERMEDIATE spectrum, so it belongs in the E line -- the one for
# processes with initial, intermediate and final states -- and not in the D line, which is for single-step
# processes.  And it is now driven through `Atomic.Computation` like every other process in JAC, where before it
# called `MultiPhotonIonization.computeLines2pOneElectron` directly.
#
# WHAT THE MOVE CHANGES FOR THE READER, and it is an improvement rather than a formality: the intermediate spectrum
# is now passed as `intermediateStates`, the same convention `MultiPhotonTransition`, `PhotoExcitationFluores` and
# `PhotonScattering` use.  A second-order amplitude is a sum over a spectrum, and WHICH spectrum changes the answer;
# making it an explicit argument puts that choice in front of the user instead of leaving it to whatever orbitals
# happened to be lying around.  Branch b below is exactly a study of that choice, so it reads better this way.
#
# THE MODULE REMAINS A ONE-ELECTRON ONE.  `MultiPhotonIonization.computeLines` is an ADAPTER: it checks that the
# initial state really has a single electron and refuses anything else with a message that says why.  A
# many-electron two-photon ionization is real physics that nobody has written, not a keyword away.
#
using SymEngine

setDefaults("method: continuum, Galerkin");   setDefaults("method: normalization, pure sine")

if  true
    # Last visit:      01-Sep-2026
    # Last successful: unknown ... ROUTE CHANGED 01-Sep-2026; the numbers quoted below were obtained on the
    #                  direct route in example-Dv.jl and have NOT yet been re-verified through Atomic.Computation.
    #
    # Branch a: THE REFERENCE COMPUTATION.  Hydrogen 1s, two linearly polarized photons of equal energy, well below
    #   the 1s -> 2p resonance at 0.375 a.u.  sigma^(2) has units of cm^4 s rather than of an area, because a
    #   two-photon rate goes as the flux SQUARED; it is quoted in Goeppert-Mayer, 1 GM = 1e-50 cm^4 s, and
    #   1 a.u. = 1.8968 GM.
    setDefaults("print summary: open", "zzz-MultiPhotonIonization-Ei-reference.sum")
    # THE BOX HOLDS THE INTERMEDIATE SET, and that is a change the move forced into the open. The retired
    # example-Dv.jl summed np up to n = 8 in a 90 a.u. box by building the orbitals directly, which BYPASSED
    # Bsplines.checkGridRepresentation; through Atomic.Computation the guard runs and refuses that grid -- 8p at
    # Z = 1 turns over near 127 a.u. and does not fit. The reference branch therefore stops at 4p, which the same
    # 90 a.u. box holds properly (4p turns over near 32 a.u.). Going higher is branch b's business, and there the
    # high np are PSEUDO-states rather than physical ones, so the box must grow with them -- and a box GROWN for
    # the bound states then runs into the CONTINUUM normalisation, which refused a 150 a.u. box here. The two
    # demands pull in opposite directions, which is the whole difficulty of a second-order continuum process.
    grid = Radial.Grid(Radial.Grid(false); rnt=2.0e-6, h=5.0e-2, hp=1.0e-2, rbox=90.0)
    nm   = Nuclear.Model(1.0, Nuclear.PointNucleus(), 0., 0., AngularJ64(0), 0., 0., 0.)
    scfN = Basics.NuclearField()
    asfN = AsfSettings(AsfSettings(); scField=scfN, gridStopper=false)
    #
    # the intermediate spectrum, supplied explicitly: the np series is what an E1 step from 1s can reach
    interConfs = [Configuration("$(n)p")  for n = 2:4]
    interRep   = Representation("intermediate np states", nm, grid, interConfs,
                                MeanFieldMultiplet(MeanFieldSettings(scfN)))
    interMp    = generate(interRep, output=true)["mean-field multiplet"]
    #
    mpiSettings = MultiPhotonIonization.Settings(MultiPhotonIonization.Settings();
                      scheme = MultiPhotonIonization.TwoPhotonOneElectronScheme([0.30, 0.32, 0.35], [E1], 1.0e-3),
                      intermediateStates = interMp, printBefore = true)
    wa = Atomic.Computation(Atomic.Computation(), name="Ei-a: two-photon ionization of H", grid=grid, nuclearModel=nm,
                            initialConfigs = [Configuration("1s")],    initialAsfSettings = asfN,
                            finalConfigs   = [Configuration("1s^0")],  finalAsfSettings   = asfN,
                            processSettings = mpiSettings)
    wb = perform(wa; output=true)
    setDefaults("print summary: close", "")
    #
elseif false
    # Last visit:      01-Sep-2026
    # Last successful: unknown ... ROUTE CHANGED 01-Sep-2026; see branch a.
    #
    # Branch b: DOES THE INTERMEDIATE SUM CONVERGE?  A second-order amplitude sums over the WHOLE spectrum, bound
    #   AND continuum, and this branch adds bound states one shell at a time and watches.  With the new route the
    #   study is simply a loop over `intermediateStates`, which is what makes the question visible.
    #
    # REPORT (25-Aug-2026, on the DIRECT route of the retired example-Dv.jl), hydrogen 1s at omega = 0.30 a.u.,
    #   90 a.u. box, sigma^(2) linear in ATOMIC units:
    #   THE SUM CONVERGES, and slowly: the gauge ratio climbs monotonically towards 1 but is still 9 % away with
    #   44 bound intermediate states, which is the expected behaviour of a bound-only sum -- length and velocity
    #   forms agree exactly only for a COMPLETE set, continuum included.
    #   THIS BRANCH FIRST GAVE THE OPPOSITE ANSWER, and the mistake is worth keeping.  Run in a 60 a.u. box the
    #   ratios marched monotonically AWAY from 1 -- 0.775, 0.553, 0.396, 0.198 -- which reads as "the bound sum
    #   does not converge" and was half written up as exactly that.  In a 90 a.u. box the same study converges:
    #   0.793, 0.803, 0.873, 0.908.  The corrupted continuum orbital grew worse as the sum grew, so the artifact
    #   and the physics carried the SAME signature, and it was the DIRECTION of the march that made it convincing.
    #   See Rule 12: a box that is merely too small fails quietly and progressively.
    setDefaults("print summary: open", "zzz-MultiPhotonIonization-Ei-convergence.sum")
    grid = Radial.Grid(Radial.Grid(false); rnt=2.0e-6, h=5.0e-2, hp=1.0e-2, rbox=90.0)
    nm   = Nuclear.Model(1.0, Nuclear.PointNucleus(), 0., 0., AngularJ64(0), 0., 0., 0.)
    scfN = Basics.NuclearField();    asfN = AsfSettings(AsfSettings(); scField=scfN, gridStopper=false)
    for  nmax  in  [3, 4, 6, 8]
        println("\n  ---- intermediate np states up to n = $nmax ----")
        interConfs = [Configuration("$(n)p")  for n = 2:nmax]
        interMp    = generate(Representation("inter", nm, grid, interConfs,
                              MeanFieldMultiplet(MeanFieldSettings(scfN))), output=true)["mean-field multiplet"]
        mpiSettings = MultiPhotonIonization.Settings(MultiPhotonIonization.Settings();
                          scheme = MultiPhotonIonization.TwoPhotonOneElectronScheme([0.30], [E1], 1.0e-3),
                          intermediateStates = interMp)
        wa = Atomic.Computation(Atomic.Computation(), name="Ei-b: n <= $nmax", grid=grid, nuclearModel=nm,
                                initialConfigs = [Configuration("1s")],    initialAsfSettings = asfN,
                                finalConfigs   = [Configuration("1s^0")],  finalAsfSettings   = asfN,
                                processSettings = mpiSettings)
        perform(wa)
    end
    setDefaults("print summary: close", "")
    #
elseif false
    # Last visit:      01-Sep-2026
    # Last successful: unknown ... ROUTE CHANGED 01-Sep-2026; see branch a.
    #
    # Branch c: THE RESONANCE GUARD, which REFUSES rather than skips.  If a photon energy brings a real
    #   intermediate state onto the energy shell, the perturbative denominator vanishes and the expression is not
    #   defined there; dropping the term instead would remove most of the amplitude and hand back a small number
    #   in place of no number.
    #
    # REPORT (25-Aug-2026): asking for hydrogen at omega = 0.375 a.u., the 1s -> 2p resonance, the guard trips and
    #   names 2p_1/2.  Below it, at omega = 0.35, the same computation runs normally.
    #   THE TOLERANCE HAD TO BE MADE RELATIVE, and the first version was useless.  With an ABSOLUTE tolerance of
    #   1e-6 a.u. the guard did NOT trip at omega = 0.375: the computed 1s -> 2p transition sits at 0.3750046, so
    #   the denominator was 4.6e-6 -- "safe" by that test -- while inflating that one term by 2e5.  It is now
    #   relative to the photon energy, refusing at |denominator| < tolerance * omega with a default of 1e-3.
    setDefaults("print summary: open", "zzz-MultiPhotonIonization-Ei-resonance.sum")
    grid = Radial.Grid(Radial.Grid(false); rnt=2.0e-6, h=5.0e-2, hp=1.0e-2, rbox=90.0)
    nm   = Nuclear.Model(1.0, Nuclear.PointNucleus(), 0., 0., AngularJ64(0), 0., 0., 0.)
    scfN = Basics.NuclearField();    asfN = AsfSettings(AsfSettings(); scField=scfN, gridStopper=false)
    interConfs = [Configuration("$(n)p")  for n = 2:8]
    interMp    = generate(Representation("inter", nm, grid, interConfs,
                          MeanFieldMultiplet(MeanFieldSettings(scfN))), output=true)["mean-field multiplet"]
    for  omega  in  [0.35, 0.375]
        println("\n  ---- omega = $omega a.u. " * (omega == 0.375 ? "(the 1s -> 2p resonance; must be REFUSED)" :
                                                                   "(below the resonance; must run)") * " ----")
        mpiSettings = MultiPhotonIonization.Settings(MultiPhotonIonization.Settings();
                          scheme = MultiPhotonIonization.TwoPhotonOneElectronScheme([omega], [E1], 1.0e-3),
                          intermediateStates = interMp)
        wa = Atomic.Computation(Atomic.Computation(), name="Ei-c: omega = $omega", grid=grid, nuclearModel=nm,
                                initialConfigs = [Configuration("1s")],    initialAsfSettings = asfN,
                                finalConfigs   = [Configuration("1s^0")],  finalAsfSettings   = asfN,
                                processSettings = mpiSettings)
        try     perform(wa)
        catch e print("\n  >> REFUSED, as it should be: ");  println(split(sprint(showerror, e), "\n")[3])
        end
    end
    setDefaults("print summary: close", "")
    #
end
