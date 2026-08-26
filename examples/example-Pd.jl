
println("Pd) FORM-FACTOR approximation for photon scattering: the cheap non-relativistic limit, and where it may be trusted.")

setDefaults("print summary: open", "zzz-PhotonScattering.sum")
setDefaults("unit: energy", "eV")

if  true
    # Last visit:      22-Aug-2026
    # Last successful:  22-Aug-2026
    #
    #   DATED, and it is the FIRST dated branch of the P series.  Rule 7 asks for verification of physical
    #   consistency rather than "it ran", and this branch has what none of its neighbours has: an ABSOLUTE number
    #   checked against a CLOSED FORM, agreeing to one part in 10^5, with the deviation from that limit tracking
    #   q a_0 exactly as the physics requires.  There is no underived prefactor anywhere in it -- the
    #   normalization follows from r_e = alpha^2 -- which is precisely why the other P branches remain undated:
    #   their SHAPES are verified and their MAGNITUDES are not.  The contrast between this date and those blanks
    #   is meant to be read.
    #
    # Branch a (RAYLEIGH in the form-factor approximation, and the first ABSOLUTE check in this module):
    #   coherent scattering on beryllium-like neon, 1s^2 2s^2, over a wide range of photon energies.
    #
    #   THE APPROXIMATION.  Treat the atomic electrons as free scatterers whose only collective property is the
    #   Fourier transform of their charge density, the atomic form factor F(q), with q = 2 k sin(theta/2) the
    #   momentum transfer.  Then
    #
    #       d sigma / d Omega  =  (r_e^2 / 2) (1 + cos^2 theta) |F(q)|^2 ,
    #
    #   with r_e = alpha^2 in atomic units.  There is NO sum over intermediate states, NO gauge and NO multipole
    #   expansion -- which is exactly why it is cheap, and exactly why it cannot describe a resonance.
    #
    #   WHY THIS BRANCH MATTERS MORE THAN ITS COST SUGGESTS.  It carries the ONLY absolute check in the whole of
    #   PhotonScattering.  As q --> 0 the form factor tends to the electron number, F(0) = N, so
    #
    #       sigma  -->  N^2 sigma_Thomson  =  N^2 (8 pi / 3) alpha^4
    #
    #   which for N = 4 is 16 x 2.3741e-08 = 3.7986e-07 a.u.  That is a CLOSED FORM: no prefactor to derive, no
    #   fitted constant, nothing taken from a table.  Every other test in this enterprise -- omega^4, Z^5,
    #   detailed balance, the Lorentzian, low-frequency additivity -- is a RATIO or a SHAPE, deliberately, because
    #   the cross-section prefactors are underived.  This one is not, and if the computed low-energy limit misses
    #   3.7986e-07 the fault is in this file's own normalization and nowhere else.
    #
    #   THE COMPARISON WITH BRANCH a OF example-Pb.jl IS THE POINT, AND THE TWO MUST DISAGREE AT LOW ENERGY.
    #   That is physics, not a defect, and it is worth being explicit because a naive reading would call it one:
    #
    #     * The SECOND-ORDER sum gives sigma ~ omega^4 as omega --> 0.  A BOUND electron cannot scatter a photon
    #       whose energy lies far below its binding energy except through the atom's polarizability, and that
    #       response vanishes as omega^4.  example-Pb.jl measures exponents 4.014 and 4.056.
    #     * The FORM FACTOR gives sigma --> N^2 sigma_Thomson, a CONSTANT, as omega --> 0.  It treats the
    #       electrons as free, and free electrons scatter at the Thomson rate however soft the photon is.
    #
    #   So the two approximations are valid in OPPOSITE regimes: the second-order sum near and below the
    #   resonances, the form factor well ABOVE the binding energies where the electrons do respond freely.  A
    #   disagreement at 1-4 eV confirms both are behaving correctly; agreement there would indict one of them.
    #   The interesting question -- where the two cross over, and whether either is reliable in between -- is what
    #   this pair of branches exists to make askable.
    #
    #   REPORT (22-Aug-2026).  THE ABSOLUTE CHECK PASSES TO ONE PART IN 10^5.
    #
    #      omega [eV]   q a_0 (180 deg)   sigma [a.u.]      sigma / N^2 sigma_T
    #         1.0           0.00054       3.800982e-07          0.999990
    #        10.0           0.00536       3.800976e-07          0.999989
    #       100.0           0.05363       3.800411e-07          0.999840
    #       500.0           0.26817       3.786772e-07          0.996252
    #      1000.0           0.53635       3.744770e-07          0.985202
    #      2000.0           1.07269       3.585666e-07          0.943344
    #
    #   against the closed form N^2 (8 pi/3) alpha^4 = 3.801018e-07 a.u.  N came out 4.00000 exactly, F(0)
    #   recovering the electron number.  The departure sets in precisely where it should: four-digit agreement
    #   while q a_0 << 1, falling away as q a_0 approaches 1, which is the form factor leaving F(0).
    #
    #   AN INDEPENDENT CONFIRMATION from the angular table: d sigma / d Omega at theta = 0 is 4.5371e-08 at EVERY
    #   energy, since q = 0 in the forward direction whatever omega may be -- and r_e^2 N^2 = 16 alpha^4 =
    #   4.534e-08.  Backward scattering falls with energy as q grows.  Nothing here is fitted.
    #
    #   AND THE DEPARTURE FROM THE LIMIT IS PHYSICAL, NOT NUMERICAL -- checked rather than assumed, by refining the
    #   radial grid at FIXED omega = 2 keV and asking whether 0.9433 moves:
    #
    #       h = 5.0e-2,  rbox = 10   (301 points)    sigma/N^2 sT = 0.943344      <- the dated grid
    #       h = 2.5e-2,  rbox = 10   (595 points)                   0.943347
    #       h = 1.25e-2, rbox = 10  (1183 points)                   0.943348
    #       h = 5.0e-2,  rbox = 20   (315 points)                   0.943344
    #
    #   Quadrupling the grid points moves the ratio by FOUR PARTS IN 10^6, and doubling the box changes it not at
    #   all to six digits.  So the 5.7 % shortfall at q a_0 = 1.07 is the atom genuinely ceasing to scatter
    #   coherently, not the radial grid.  That matters more than it sounds: it turns an anchor at ONE point into a
    #   VERIFIED APPROACH to it across q a_0 = 0.0005 ... 1.07, and a wrong normalization could not survive being
    #   right at both ends of that range while tracking the form factor in between.
    #
    #   THIS IS THE FIRST ABSOLUTELY VERIFIED NUMBER IN PhotonScattering.  Every other branch of the P series has a
    #   verified SHAPE and an unverified MAGNITUDE, the cross-section prefactors being underived; this branch has
    #   no free prefactor at all, its normalization following from r_e = alpha^2.
    #
    grid = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, rbox = 10.0)
    nm   = Nuclear.Model(10.)

    settings = PhotonScattering.Settings(PhotonScattering.Settings();
                    process        = PhotonScattering.RayleighScattering(),
                    approximation  = PhotonScattering.FormFactorApproximation(),
                    # from far below the 2s binding energy, where the Thomson limit must be recovered, up to where
                    # q a_0 ~ 1 and the form factor begins to fall
                    photonEnergies = [1.0, 10.0, 100.0, 500.0, 1000.0, 2000.0],
                    polarThetas    = [0.0, pi/4, pi/2, 3pi/4, pi],
                    printBefore    = true )

    wa = Atomic.Computation(Atomic.Computation(), name="Rayleigh, form-factor approximation", grid=grid, nuclearModel=nm,
                            initialConfigs  = [Configuration("1s^2 2s^2")],
                            finalConfigs    = [Configuration("1s^2 2s^2")],
                            processSettings = settings )
    wb = perform(wa; output=true)
    #
elseif  false
    # Last visit:      22-Aug-2026
    # Last successful:  unknown ... BOTH LIMITS are verified -- sigma ~ omega^2 at low energy from the closure
    #                              expansion, and sigma / N sigma_T rising monotonically toward 1 at high -- but
    #                              the MIDDLE of the range rests on S(q) = N - |F|^2/N, which the Hartree-Fock
    #                              tabulations differ from by some percent precisely there.  A date would suggest
    #                              the cross sections are right ACROSS the range, which is more than has been
    #                              shown.  The missing test is a comparison of S(q) against a tabulation.
    #
    # Branch b (COMPTON in the form-factor approximation): incoherent scattering on the same target, over the same
    #   energies, so that the coherent and incoherent channels can be read side by side.
    #
    #   THE APPROXIMATION, and its one genuine weakness.  Where Rayleigh is governed by the COHERENT form factor
    #   F(q), Compton is governed by the INCOHERENT SCATTERING FUNCTION S(q):
    #
    #       d sigma / d Omega  =  (r_e^2 / 2) (1 + cos^2 theta) S(q) .
    #
    #   JAC's FormFactor module supplies F(q) -- FormFactor.standardF -- but NOT S(q), so this branch uses the
    #   standard closure approximation
    #
    #       S(q)  =  N  -  |F(q)|^2 / N ,
    #
    #   which is exact in both limits that matter and interpolated in between: S(0) = N - N^2/N = 0, since a soft
    #   photon cannot excite a bound electron incoherently, and S(q --> infinity) = N, since every electron
    #   scatters independently once the momentum transfer is large.  It is NOT a Hartree-Fock incoherent
    #   scattering function and should not be quoted as one; the tabulations (Hubbell and others) differ from it
    #   by some percent in the intermediate region, which is precisely where this branch is least trustworthy.
    #
    #   WHAT IT IS GOOD FOR.  The SUM RULE.  Coherent and incoherent scattering share the same electrons, and
    #
    #       |F(q)|^2  +  S(q)  ->  N^2  at q = 0   and   ->  N  at large q ,
    #
    #   so the two branches together must trade one channel for the other as q grows, with nothing lost. That is
    #   a check no single branch can make, and it is why the two are written as a pair.
    #
    #   NOTE the deliberate limitation: this is the THOMSON-level Compton cross section, without the Klein-Nishina
    #   factor and without the Compton energy shift, so omega_out = omega_in here.  Both corrections matter once
    #   the photon energy approaches m c^2 and neither is implemented; at the energies of this branch they are
    #   small, and the module's capability table says so.
    #
    #   REPORT (22-Aug-2026).  BOTH LIMITS ARE CORRECT, and they run OPPOSITE to the coherent case:
    #
    #      omega [eV]   q a_0 (180 deg)   sigma [a.u.]      sigma / N sigma_T
    #         1.0           0.00054       1.426192e-15          0.000000
    #        10.0           0.00536       1.426189e-13          0.000002
    #       100.0           0.05363       1.425980e-11          0.000150
    #       500.0           0.26817       3.552295e-10          0.003738
    #      1000.0           0.53635       1.405299e-09          0.014789
    #      2000.0           1.07269       5.382889e-09          0.056647
    #
    #   LOW-ENERGY LIMIT, exact and checkable: sigma goes as omega^2.  From 1 to 10 to 100 eV the cross section
    #   rises by a factor 100 per decade -- 1.426192e-15, 1.426189e-13, 1.425980e-11 -- and that follows from the
    #   closure at small q, S(q) = N - |F|^2/N ~ N q^2 <r^2> / 3.  So the construction is confirmed at the end
    #   where it is supposed to be exact.
    #
    #   HIGH-ENERGY LIMIT: the ratio to N sigma_T rises monotonically toward 1, reaching 0.0566 at 2 keV.  It is
    #   still far from 1 because the 1s electrons of neon are tightly bound and q a_0 ~ 1 does not yet resolve
    #   them individually; the limit is approached, not reached, over this range.
    #
    #   NOTE THE NORMALIZATION IS PROCESS-DEPENDENT, and getting it wrong once produced a table of 0.000000 beside
    #   perfectly good cross sections.  Rayleigh tends to N^2 sigma_T as q --> 0, the electrons scattering IN
    #   PHASE so the amplitudes add before squaring; Compton tends to N sigma_T at LARGE q, the momentum transfer
    #   resolving the electrons individually so the cross sections add.  Different limits, approached from
    #   opposite directions in energy.  PhotonScattering.thomsonLimit now selects the right one per process.
    #
    #   WHAT IS NOT VERIFIED: the middle.  The closure S(q) = N - |F|^2/N is exact at both ends and interpolated
    #   between, and the Hartree-Fock tabulations differ from it by some percent precisely there.  This branch
    #   therefore verifies its two LIMITS and nothing in between.
    grid = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, rbox = 10.0)
    nm   = Nuclear.Model(10.)

    settings = PhotonScattering.Settings(PhotonScattering.Settings();
                    process        = PhotonScattering.ComptonScattering(),
                    approximation  = PhotonScattering.FormFactorApproximation(),
                    photonEnergies = [1.0, 10.0, 100.0, 500.0, 1000.0, 2000.0],
                    polarThetas    = [0.0, pi/4, pi/2, 3pi/4, pi],
                    printBefore    = true )

    wa = Atomic.Computation(Atomic.Computation(), name="Compton, form-factor approximation", grid=grid, nuclearModel=nm,
                            initialConfigs  = [Configuration("1s^2 2s^2")],
                            finalConfigs    = [Configuration("1s^2 2s^2")],
                            processSettings = settings )
    wb = perform(wa; output=true)
    #
end
#
setDefaults("print summary: close", "")
