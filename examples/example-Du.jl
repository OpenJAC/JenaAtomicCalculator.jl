
println("Du) Apply & test the MultiPhotonTransition module: THREE-photon processes between bound levels.")

setDefaults("print summary: open", "zzz-MultiPhotonTransition-3p.sum")
setDefaults("unit: energy", "eV")
setDefaults("unit: rate", "1/s")

## THIS FILE IS THE THREE-PHOTON HALF of the MultiPhotonTransition module; example-Dh.jl carries the two-photon
## half. The split is by PHOTON NUMBER, not by direction (emission vs absorption), because bichromatic absorption
## IS two-photon absorption and only the photon-number axis changes what has to be built: a third-order amplitude,
## a two-dimensional sharing simplex, and the coupling of three multipoles rather than two.
##
## THE TWO HALVES OF THREE-PHOTON ARE IN DIFFERENT STATES, and the asymmetry is deliberate rather than an
## oversight:
##   a)  three-photon EMISSION      -- only the ENERGY SHARINGS exist, and they are verified exactly. The
##                                     third-order amplitude is NOT implemented and a full run refuses.
##   b)  three-photon ABSORPTION    -- implemented, in an ELEMENTARY formulation (see below), and it runs.
##
## WHY ABSORPTION COULD BE DONE FIRST AND EMISSION COULD NOT. Absorption fixes the three photon energies, so it
## needs no sharings at all. Emission fixes only their SUM, so the rate is differential over a two-dimensional
## simplex, and that quadrature -- not the amplitude -- was the self-contained piece worth finishing first.
##
## WHAT "ELEMENTARY" MEANS FOR THE ABSORPTION SCHEME, stated up front so that nothing here is over-read:
##   * THREE BEAMS, ALL LINEARLY POLARIZED ALONG THE SAME AXIS. One geometry; no polarization observables, no
##     Stokes parameters, no rank-K decomposition. With every polarization along z only the q = 0 component
##     contributes and the magnetic quantum number is conserved along the whole chain, so the four-fold m-sum
##     collapses to a single M -- which is exactly what makes this the case to do first.
##   * A STRENGTH, NOT A CROSS SECTION:  S^(3) = 1/(2J_i+1) sum_M |A(M)|^2 in atomic units. A generalized
##     three-photon cross section needs a normalisation of order F^3, and the module's TWO-photon absorption
##     normalisation has never been derived either -- so inventing a three-photon one would add a second
##     undetermined constant wearing the units of a measured quantity.
##
## WHY THE AMPLITUDE IS BUILT BY SUMMING OVER MAGNETIC QUANTUM NUMBERS rather than by coupling the three photon
## multipoles to a total rank K. The coupled route needs the six time orderings re-expressed in one common
## coupling order, which brings in recoupling coefficients (6-j, and 9-j for unequal multipoles) whose phases are
## easy to get confidently wrong and hard to test -- and getting exactly that kind of phase wrong is what blocker
## A1 was, in the two-photon file. The m-sum route needs no recoupling at all: each ordering is a plain product
## of three Wigner-Eckart 3-j symbols and three reduced matrix elements. It costs an m-loop and buys checkability.
##
## AND NOTE WHAT IS THEREFORE *NOT* A TEST HERE. For two photons, invariance under exchanging the colours was a
## real check, because the two orderings were combined with a relative phase that could be (and was) wrong. Here
## all six orderings are summed explicitly with no relative phase to get wrong, so permuting the three colours is
## invariant BY CONSTRUCTION and proves nothing. The checks with teeth are in branch b.
##
## REFERENCES for a later verification stage (no three-photon benchmark is used or claimed here):
##   Drake, PRA 34, 2871 (1986)                     -- two-photon rates; the anchor for the two-photon half
##   Fritzsche, Molecules 26, 2660 (2021)           -- the intermediate-state summation itself

grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 20.0)


##
## =====================================================================================================
##  THREE-PHOTON EMISSION -- the energy sharings only; there is no amplitude
## =====================================================================================================
if  false
    # Last successful:  08-Aug-2026        ## for the SIMPLEX SHARINGS only; there is still no amplitude
    #
    # --- Branch a: the ENERGY SHARINGS on the two-dimensional simplex, and their sum rule.
    #
    # MOVED HERE 08-Aug-2026 from example-Dh.jl, where it had been placed as branch j against the agreed
    # organisation -- three-photon and beyond belong in this file.
    #
    # WHY THIS PIECE FIRST. Three-photon emission needs a third-order amplitude (a double sum over two
    # intermediate sets, all 3! = 6 orderings, three coupled multipoles) which is not implemented. But ONE piece
    # of it is self-contained and has an EXACT answer to check against: the sharings on the simplex
    # omega1 + omega2 + omega3 = E_i - E_f. Get that quadrature wrong and the eventual three-photon rate is wrong
    # by a factor that no amount of checking the amplitude would ever reveal -- the silent-normalisation class of
    # error that cost this module a factor 4 in the two-photon totals, found only on 08-Aug-2026.
    #
    # The moments of a simplex are known in closed form,
    #     Int w1^a w2^b w3^c dw1 dw2  =  E^(a+b+c+2) a! b! c! / (a+b+c+2)!
    # so the weights are verified against exact numbers before any physics is attached to them.
    #
    # MEASURED (work/diag-3p-sharings.jl), relative error against the exact moments:
    #     n = 2 (4 points):   1, w1, w1w2 exact to 1e-16;  w1w2w3 off by 1.1e-1
    #     n = 4 (16 points):  ALL FOUR exact to 1e-15
    #     n = 6, 8:           all four exact to 1e-15
    # The n = 2 failure is not a defect but the confirmation: with the Jacobian the (1,1,1) integrand is degree 4
    # in the collapsed coordinate, and a 2-point Gauss-Legendre rule is exact only to degree 3. The quadrature
    # fails exactly where degree counting says it must, and nowhere else.
    #
    # ONE LIMITATION, documented rather than discovered later: the collapsed-coordinate grid is EXACT but not
    # PERMUTATION SYMMETRIC -- omega1 is distinguished by the transform. For the integrated rate that is
    # irrelevant, but a spectrum tabulated on these nodes will not LOOK symmetric under exchanging the photons
    # even though the function is. The two-photon symmetry check must therefore be applied to the function at
    # permuted points, not to the node table.
    #
    # HOW TO SEE IT: calcOverview = true prints the sharings and the sum-rule table. A full run still REFUSES,
    # and its message lists what is missing AND states that the sharings are already done.
    ni          = Nuclear.Model(1.0, PointNucleus())   ## Fermi cannot represent Z = 1
    mpSettings  = MultiPhotonTransition.Settings(MultiPhotonTransition.Settings();
                        scheme = MultiPhotonTransition.ThreePhotonEmissionScheme(
                                     MultiPhotonTransition.AbstractMultiPhotonProperty[
                                         MultiPhotonTransition.EnergyDiffCs()], 4 ),  ## points PER DIRECTION
                        multipoles = [E1], gauges = [UseCoulomb],
                        intermediateStates = Multiplet(), calcOverview = true,
                        lineSelection = LineSelection() )
    wa = Atomic.Computation(Atomic.Computation(), name="Du-a: three-photon energy sharings on the simplex", grid=grid,
                            nuclearModel   = ni,
                            initialConfigs = [Configuration("3s")],
                            finalConfigs   = [Configuration("1s")],
                            processSettings= mpSettings )
    perform(wa)
    #
##
## =====================================================================================================
##  THREE-PHOTON ABSORPTION -- elementary, and it runs
## =====================================================================================================
elseif  true
    # Last successful:  08-Aug-2026        ## for the SELECTION RULE, the ordering count and the Z-scaling
    #
    # DATED for what it verifies -- three exact, parameter-free structural checks -- and NOT for the magnitude,
    # which carries no derived normalisation at all (see the header).
    #
    # --- Branch b: THREE-PHOTON ABSORPTION in hydrogen, 1s -> 2p, three parallel linearly-polarized beams.
    #
    # 1s -> 2p is deliberately a transition that ONE photon can also drive. Three-photon absorption is not being
    # proposed here as the way to reach 2p; the point is that the third-order machinery can be exercised on the
    # simplest possible system, where the intermediate spectrum is understood and the checks below are sharp.
    #
    # THREE CHECKS, ALL PASSED (work/diag-3p-absorption.jl), none of which needs a benchmark or a normalisation:
    #
    #   (1) THE PARITY SELECTION RULE, which is the sharp one. Three E1 photons change the parity, so an
    #       even -> even transition cannot be driven at all:
    #           1s -> 2s :  S^(3) = 0.000000e+00 in BOTH gauges, exactly
    #           1s -> 2p :  S^(3) = 3.747137e-18 (Cou) / 5.930890e-18 (Bab) for J_f = 1/2
    #                       S^(3) = 1.594001e-17 (Cou) / 1.637483e-17 (Bab) for J_f = 3/2
    #       An exact zero from a routine that is demonstrably non-zero elsewhere is worth more than a plausible
    #       number: it says the parity rule is enforced by the machinery and not by a filter written around it.
    #
    #   (2) THE SIX TIME ORDERINGS. For equal colours all six coincide, so the amplitude must be exactly six
    #       times the single-ordering one and the strength exactly 36 times. MEASURED: 36.000000 for both J_f.
    #       This is the check that the ordering bookkeeping is complete -- five missing orderings would still
    #       have produced a perfectly plausible number.
    #
    #   (3) HYDROGENIC Z-SCALING. Each length-form E1 amplitude in JAC carries a factor omega (measured
    #       explicitly, work/diag-omega-dependence.jl), so it scales as omega * r ~ Z^2 * Z^-1 = Z; the two
    #       energy denominators scale as Z^2 each. Hence A ~ Z^3/Z^4 = Z^-1 and S^(3) ~ Z^-2.
    #       MEASURED over Z = 1, 2, 4:  power -2.0020 and -2.0011, the residual being the relativistic
    #       correction. A clean integer power confirms that the three amplitudes and the two denominators are
    #       combined consistently.
    #
    # THE GAUGE RATIO IS NOT UNIFORM -- Cou/Bab is 0.63 for J_f = 1/2 and 0.97 for J_f = 3/2 -- and that is
    # expected rather than alarming here: the velocity form carries the LEVEL gap instead of the photon energy
    # in an off-shell intermediate step (found 08-Aug-2026, see the note at the bichromatic ranking display), so
    # a third-order sum with TWO off-shell steps is exactly where the two gauges should be expected to part
    # company. Babushkin is the form to read.
    #
    # ON THE POTENTIAL: NuclearField throughout, so that the initial, final and intermediate states are
    # eigenstates of the SAME one-body Hamiltonian -- the consistency requirement that was worth a factor 6 in
    # the two-photon work.
    ni          = Nuclear.Model(1.0, PointNucleus())
    scfN        = Basics.NuclearField()
    asfN        = AsfSettings(AsfSettings(); scField = scfN)
    interConfs  = [Configuration("2s"), Configuration("2p"), Configuration("3s"), Configuration("3p"),
                   Configuration("3d"), Configuration("4p")]
    interRep    = Representation("intermediate levels", ni, grid, interConfs,
                                 MeanFieldMultiplet(MeanFieldSettings(scfN)))
    interMp     = generate(interRep, output=true)["mean-field multiplet"]
    #
    ## omega1 = omega2 = 0 selects the MONOCHROMATIC case, all three photons carrying (E_f - E_i)/3. Give two
    ## energies in the user-selected units instead to drive it with three different colours; the third then
    ## follows from energy conservation.
    mpSettings  = MultiPhotonTransition.Settings(MultiPhotonTransition.Settings();
                        scheme = MultiPhotonTransition.ThreePhotonAbsorptionScheme(
                                     MultiPhotonTransition.AbstractMultiPhotonProperty[
                                         MultiPhotonTransition.TotalCsUnpolarized()], 0., 0. ),
                        multipoles = [E1], gauges = [UseCoulomb, UseBabushkin],
                        intermediateStates = interMp, calcOverview = false,
                        lineSelection = LineSelection() )
    wb = Atomic.Computation(Atomic.Computation(), name="Du-b: three-photon absorption of H, 1s -> 2p", grid=grid,
                            nuclearModel   = ni,
                            initialConfigs = [Configuration("1s")],
                            finalConfigs   = [Configuration("2p")],
                            initialAsfSettings = asfN, finalAsfSettings = asfN,
                            processSettings= mpSettings )
    perform(wb)
    #
elseif  false
    # Last visit:  08-Aug-2026
    #
    # --- Branch c: THE SAME TRANSITION WITH THREE DIFFERENT COLOURS, to show that the scheme is not restricted
    #     to the monochromatic case.
    #
    # NOT A TEST, and labelled so: permuting the three colours is invariant by construction here (all six
    # orderings are summed explicitly, with no relative phase to get wrong), so this branch demonstrates the
    # capability rather than checking it. What it does show is the physically interesting part -- how strongly
    # the strength depends on how close one of the partial sums E_i + omega1 or E_i + omega1 + omega2 comes to a
    # real intermediate level. That near-resonant enhancement is the whole reason multi-colour schemes are used
    # in experiments, and it is where the non-resonant expression computed here eventually stops applying.
    ni          = Nuclear.Model(1.0, PointNucleus())
    scfN        = Basics.NuclearField()
    asfN        = AsfSettings(AsfSettings(); scField = scfN)
    interConfs  = [Configuration("2s"), Configuration("2p"), Configuration("3s"), Configuration("3p"),
                   Configuration("3d"), Configuration("4p")]
    interRep    = Representation("intermediate levels", ni, grid, interConfs,
                                 MeanFieldMultiplet(MeanFieldSettings(scfN)))
    interMp     = generate(interRep, output=true)["mean-field multiplet"]
    #
    ## the H 1s -> 2p transition energy is 10.2 eV; 2.0 + 3.0 + 5.2 eV is a strongly asymmetric sharing
    mpSettings  = MultiPhotonTransition.Settings(MultiPhotonTransition.Settings();
                        scheme = MultiPhotonTransition.ThreePhotonAbsorptionScheme(
                                     MultiPhotonTransition.AbstractMultiPhotonProperty[
                                         MultiPhotonTransition.TotalCsUnpolarized()], 2.0, 3.0 ),  ## [eV]
                        multipoles = [E1], gauges = [UseCoulomb, UseBabushkin],
                        intermediateStates = interMp, calcOverview = false,
                        lineSelection = LineSelection() )
    wc = Atomic.Computation(Atomic.Computation(), name="Du-c: three-colour three-photon absorption of H", grid=grid,
                            nuclearModel   = ni,
                            initialConfigs = [Configuration("1s")],
                            finalConfigs   = [Configuration("2p")],
                            initialAsfSettings = asfN, finalAsfSettings = asfN,
                            processSettings= mpSettings )
    perform(wc)
    #
end
#
