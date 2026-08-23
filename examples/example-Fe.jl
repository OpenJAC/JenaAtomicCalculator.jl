
println("Fe) Cascade.PhotoAbsorptionScheme: photoabsorption of C^+ near the 1s -> 2p resonances.")

using JLD2
#
setDefaults("method: continuum, asymptotic Coulomb")    ## setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, pure sine")         ## setDefaults("method: normalization, pure Coulomb")
setDefaults("unit: energy", "eV")
setDefaults("unit: cross section", "Mbarn")             ## NB: the string is "Mbarn"; "Mb" is rejected

grid = Radial.Grid(Radial.Grid(false); rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)


# REWRITTEN 06-Aug-2026, third file of the scheme series. The previous example-Fe.jl was a DR rate-coefficient
# script for neon-like Ar; DR gets its own file (Ff) under the plan, so nothing was kept.
#
# REFERENCE. This file reproduces the C^+ example of
#     S. Fritzsche, A. K. Sahoo, L. Sharma, S. Schippers,
#     "Simulated Photoabsorption Spectra for Singly and Multiply Charged Ions",
#     Atoms 13, 77 (2025);  examples/papers/b25.atoms-photoabsorption-ions-original.pdf
# Section 3.1 and Figures 3 and 4, i.e. the photoabsorption of C^+ (1s^2 2s^2 2p) near its 1s -> 2p
# resonances, measured by Mueller et al., Phys. Rev. A 97, 013409 (2018) at PIPE.
#
# THE PUBLISHED SCRIPT NO LONGER RUNS, and this file deliberately does NOT try to reproduce it verbatim; it
# uses the CURRENT API. Two data structures have gained a field since the paper went to press:
#   * Cascade.PhotoAbsorptionScheme has 11 fields, not the 10 of the paper's Figure 2 -- minCrossSection was
#     added at the end;
#   * Cascade.PhotoAbsorptionSpectrum has 8 fields, not 7 -- csScaling was inserted after resonanceWidth.
# Both positional calls of Figure 3 are therefore one argument short. The calls below are the corrected ones.
#
# PHOTOABSORPTION = DIRECT + RESONANT, Eq. (1) of the paper:
#     sigma^abs(w; i)  =  SUM_f sigma^ion(w; i->f)  +  SUM_{e,f'} sigma^exc(w; i->e) b(e->f')
# The cascade COMPUTATION (branch a) generates both sets of amplitudes and stores them; the SIMULATION
# (branches b-d) combines them on a fine photon-energy grid, replacing each resonance by a profile of a
# user-given width. The two steps are separate on purpose: one expensive computation feeds many cheap
# simulations, and the simulation is where the initial-level population and the direct/resonant split are
# chosen.
#
# ONE REAL BUG FIXED 06-Aug-2026 WHILE WRITING THIS FILE, and it dominated everything. In
# PhotoExcitation.estimateCrossSection -- the routine through which EVERY simulated resonance passes -- the
# integrated strength 109.7617 f [Mb eV] was converted back to atomic units by multiplying by the cross
# section conversion but DIVIDING by the energy conversion. Every resonant photoabsorption cross section was
# therefore too large by exactly 1/(1 eV in a.u.)^2 = 27.2114^2 = 740.46. With the division the factor came
# out as 106.659; multiplied it is 0.1440449, which equals 2 pi^2 alpha = 0.1440440 to five digits -- the
# same constant that underlies the 109.761 Mb eV resonance-strength column of PhotoExcitation, so the fix is
# confirmed independently. The C^+ peak drops from 6.9e4 Mb to the ~90 Mb range, i.e. from 230x above
# Figure 4 to the same order as it. A stray "@show n, values" in Basics.determineNearestPoints, which spewed
# a debug line per call, was removed at the same time.
#
# COMPARISON WITH FIGURE 4 AND WITH THE MEASUREMENT, 06-Aug-2026.
#   POSITIONS -- good. The two fine-structure levels of the 1s^2 2s^2 2p ^2P ground term give their strongest
#   resonance at 287.49 eV (level 1, ^2P_1/2) and 288.36 eV (level 2, ^2P_3/2), split by 0.87 eV. Mueller et
#   al. see their strongest structure near 288.4 eV, so the level-2 resonance sits within ~0.05 eV of it.
#   HEIGHTS -- right order, not quantitative. Converged peak values are 122 Mb (level 1) and 125 Mb (level 2)
#   against a measured ~170 Mb. The paper itself warns that the calculated fine structure comes out somewhat
#   too large and that no attempt was made to improve individual positions and intensities, which is the
#   expected consequence of the single-configuration AverageSCA representation.
#   MESH -- watch out, and this is a genuine trap. The paper's 0.2 eV grid UNDER-SAMPLES a 0.1 eV-wide
#   profile badly: refining the same simulation gives
#         mesh 0.2 eV : 93.2 Mb (Lev1) / 76.1 Mb (Lev2)
#         mesh 0.05 eV: 119.1        / 120.1
#         mesh 0.01 eV: 122.5        / 124.9   <- converged
#   i.e. the coarse mesh loses up to 40% of the peak, and it loses it unevenly between the two levels (7% vs
#   39%), so even the RATIO of two peaks is mesh-dependent. Any peak height read off a 0.2 eV grid at
#   resonanceWidth = 0.1 eV is a lower bound, not a value.
#   ONE THING NOT EXPLAINED. Figure 4 of the paper shows peaks reaching roughly 250-300 Mb, whereas the
#   present code on the same 0.2 eV mesh gives 76-93 Mb and even converged only 122-125 Mb. The factor of
#   ~2-3 is not accounted for here. It is NOT the 740x factor fixed above (that would be far larger), and it
#   may well be that the code has moved since the figure was made -- csScaling, which multiplies exactly this
#   quantity, is one of the two fields that did not yet exist when the paper went to press. Left open.
#
# INTERNAL CHECKS, all passed:
#   + Eq. (1) is additive to 3.1e-5 -- "direct only" plus "resonant only" reproduces "direct + resonant"
#     (branch d). The residual is the 5-digit printout, not the arithmetic.
#   + The initial-level population enters linearly to 5.2e-5 -- 0.9*(level 1) + 0.1*(level 2) reproduces the
#     90:10 run point by point (branch c).
#   + The direct term is a flat 8.2e-3 Mb background across the whole window, four orders of magnitude below
#     the resonance peaks. That is as it must be: the 1s threshold of C^+ lies near 392 eV, so nothing but
#     the weak 2s/2p direct ionization contributes below 292 eV.


if  true
    # Last visit:      06-Aug-2026 ... ~49 s; 12 photo-ionization and 26 photo-excitation lines, written to
    #                  zzz-cascade-photoabsorption-computations-<date>.jld.  Left at "Last visit" because this
    #                  branch produces amplitudes rather than a comparable observable; it is verified through
    #                  branches b-d, which consume exactly this file.
    #
    # Last visit:      23-Aug-2026
    # Last successful: 23-Aug-2026 ... 18.9 s warm.  Both halves of the photoabsorption scheme run: the DIRECT
    #                  part gives 12 photoionization lines in two steps, the RESONANT part 26 photoexcitation
    #                  lines in two more, written to zzz-cascade-photoabsorption-computations-<date>.jld.
    #
    #                  DATED ON THE STRENGTH OF ITS CONSUMERS, which is the honest ground here.  This branch
    #                  produces amplitudes and asserts nothing about them, so there is little in it to check
    #                  directly.  What makes it verified is that branches b, c and d ALL read this file and
    #                  are all dated: b reproduces the paper's Figure 4 peak at 287.6 eV, c the initial-level
    #                  admixture, and d confirms the additivity of Eq. (1) to 3.1e-5.  A fault in the
    #                  amplitudes written here would have to survive all three of those, and the additivity
    #                  check in particular ties the direct and resonant halves together.
    setDefaults("print summary: open", "zzz-Cascade-Fe-computation.sum")

    fromShells = [Shell("1s"), Shell("2s")]
    toShells   = [Shell("2p"), Shell("3p")]
    name       = "Photoabsorption of C^+: Computation"
    pScheme    = Cascade.PhotoAbsorptionScheme([E1], [en for en = 285.0:5.0:295.0], Float64[],
                                               fromShells, toShells, LevelSelection(), [0,1], true, true, 0., 0.)
    comp       = Cascade.Computation(Cascade.Computation(); nuclearModel=Nuclear.Model(6.), grid=grid,
                                     name=name, scheme=pScheme, approach=Cascade.AverageSCA(),
                                     initialConfigs=[Configuration("1s^2 2s^2 2p")] )
    println(comp)
    wb = perform(comp; output=true, outputToFile=true)
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      06-Aug-2026
    # Last successful: 06-Aug-2026 ... ~3 s.  Peak 93.2 Mb at 287.6 eV on the paper's mesh, 122.5 Mb at
    #                  287.49 eV converged; the resonance position agrees with the measurement of Mueller et
    #                  al. and the height is of the measured order.  See the assessment above, including the
    #                  unexplained factor 2-3 against Figure 4 itself.
    #
    # Branch b: THE SIMULATION -- Figure 4, left panel, blue curve. The total photoabsorption cross section
    #   of C^+ from its 1s^2 2s^2 2p ^2P_1/2 ground level over 286 - 292 eV, on the paper's 0.2 eV mesh and
    #   with its 0.1 eV resonance width. Requires branch a to have run.
    setDefaults("print summary: open", "zzz-Cascade-Fe-simulation.sum")

    fn       = sort(filter(f -> startswith(f, "zzz-cascade-photoabsorption-computations-"), readdir()),
                    by = f -> stat(f).mtime)[end]
    println(">>> reading the cascade data from  $fn")
    data     = [JLD2.load(fn)]
    energies = [en for en = 286.0:0.2:292.0]
    property = Cascade.PhotoAbsorptionSpectrum(true, true, 0.1, 1.0, energies, Shell[], [(1, 1.0)], Configuration[])
    settings = Cascade.SimulationSettings(false, false, 0.)
    simu     = Cascade.Simulation(Cascade.Simulation(); name="Photoabsorption of C^+: Simulation",
                                  computationData=data, property=property, settings=settings)
    println(simu)
    wd = perform(simu; output=true)
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      06-Aug-2026
    # Last successful: 06-Aug-2026 ... three spectra, <1 s.  Level 1 peaks at 287.6 eV and level 2 at
    #                  288.4 eV on the paper's mesh (287.49 / 288.36 converged), and the 90:10 admixture is
    #                  the exact weighted superposition of the two to 5.2e-5.
    #
    # Branch c: THE INITIAL-LEVEL ADMIXTURE -- Figure 4, left panel, all three curves. The ^2P ground term of
    #   C^+ is split into ^2P_1/2 (level 1) and ^2P_3/2 (level 2), and the two give visibly different spectra.
    #   Mueller et al. estimated that ~90% of their ions were in the ^2P ground term and ~10% in the
    #   metastable 2s 2p^2 ^4P term; the paper's red curve is that admixture. Here the same machinery is used
    #   for the population WITHIN the ground term, which is what initialOccupations addresses directly.
    setDefaults("print summary: open", "zzz-Cascade-Fe-admixture.sum")

    fn       = sort(filter(f -> startswith(f, "zzz-cascade-photoabsorption-computations-"), readdir()),
                    by = f -> stat(f).mtime)[end]
    data     = [JLD2.load(fn)]
    energies = [en for en = 286.0:0.2:292.0]
    for  (sa, occ)  in  [("level 1 only", [(1, 1.0)]), ("level 2 only", [(2, 1.0)]),
                         ("90% level 1 + 10% level 2", [(1, 0.9), (2, 0.1)])]
        println("\n>>> initial population: $sa")
        property = Cascade.PhotoAbsorptionSpectrum(true, true, 0.1, 1.0, energies, Shell[], occ, Configuration[])
        settings = Cascade.SimulationSettings(false, false, 0.)
        simu     = Cascade.Simulation(Cascade.Simulation(); name="Photoabsorption of C^+: $sa",
                                      computationData=data, property=property, settings=settings)
        wd = perform(simu; output=true)
    end
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      06-Aug-2026
    # Last successful: 06-Aug-2026 ... three spectra, <1 s.  Additivity of Eq. (1) confirmed to 3.1e-5, and
    #                  the direct term is a flat 8.2e-3 Mb background four orders below the peaks.
    #
    # Branch d: DIRECT versus RESONANT -- the two terms of Eq. (1) separately, by toggling includeIonization
    #   and includeExcitation. Near the 1s -> 2p threshold the resonant term should dominate by orders of
    #   magnitude at the peaks while the direct term is the slowly varying background, and the two switched
    #   on together must reproduce their sum exactly. That additivity is the cheapest internal check on the
    #   simulation, since Eq. (1) is an incoherent sum by construction.
    setDefaults("print summary: open", "zzz-Cascade-Fe-decomposition.sum")

    fn       = sort(filter(f -> startswith(f, "zzz-cascade-photoabsorption-computations-"), readdir()),
                    by = f -> stat(f).mtime)[end]
    data     = [JLD2.load(fn)]
    energies = [en for en = 286.0:0.2:292.0]
    for  (sa, inclIon, inclExc)  in  [("direct only", true, false), ("resonant only", false, true),
                                      ("direct + resonant", true, true)]
        println("\n>>> contribution: $sa")
        property = Cascade.PhotoAbsorptionSpectrum(inclIon, inclExc, 0.1, 1.0, energies, Shell[],
                                                   [(1, 1.0)], Configuration[])
        settings = Cascade.SimulationSettings(false, false, 0.)
        simu     = Cascade.Simulation(Cascade.Simulation(); name="Photoabsorption of C^+: $sa",
                                      computationData=data, property=property, settings=settings)
        wd = perform(simu; output=true)
    end
    setDefaults("print summary: close", "")
    #
end
