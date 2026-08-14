println("Fk) Cascade.ExpansionOpacityScheme: bound-bound expansion opacities of Sr^+ in kilonova ejecta.")

using JLD2, Printf
#
setDefaults("unit: energy", "eV")

grid = Radial.Grid(Radial.Grid(false); rnt = 1.0e-5, h = 5.0e-2, hp = 2.0e-2, rbox = 25.0)


# WRITTEN 08-Aug-2026, last stage of the cascade-scheme series.  The previous example-Fk.jl held the P^2+
# expansion-opacity material moved out of the old example-Fg.jl; none of it ran, and nothing was kept.
#
# WHAT AN EXPANSION OPACITY IS.  In ejecta that expand homologously (v proportional to r), a photon travelling
# outwards is continuously red-shifted in the local comoving frame.  It therefore sweeps across one bound-bound
# line after another, and a line that is optically thick in a narrow frequency interval nevertheless removes
# photons from a WIDE band, because every photon in that band is eventually shifted into resonance with it.  The
# effective opacity is obtained by binning the lines (Eastman & Pinto 1993),
#
#     kappa_exp(lambda) = 1/(c t_exp rho) * SUM_l (lambda_l / Delta-lambda) * (1 - exp(-tau_l))
#
# with the Sobolev optical depth of a line, in atomic units,
#
#     tau_l = pi alpha * n_l * lambda_l * t_exp * f_l ,        n_l = population of the LOWER level.
#
# The two densities play different roles and are NOT interchangeable: rho [g/cm^3] is a MASS density and only
# sets the prefactor -- it is what makes kappa come out in cm^2/g -- while n_l [1/cm^3] is a NUMBER density and
# decides how optically thick each line is.  Cascade.ExpansionOpacities now asks for both explicitly.
#
# WHY Sr^+.  The 1.0 micron absorption feature of AT2017gfo, the kilonova of GW170817, is attributed to the
# Sr II 4d - 5p triplet; Sr^+ is the single most-discussed ion in kilonova spectroscopy.  It is also cheap here:
# a Kr core plus one valence electron, so the whole cascade is seven blocks of one or two levels each.  Its
# strong lines are the 5s - 5p resonance doublet near 4100 A and the 4d - 5p triplet near 1.0 micron.
#
# SCOPE.  A single ion stage with LTE (Boltzmann) level populations, which is how the kilonova literature
# reports expansion opacities.  No Saha ionization balance across charge states -- that is a different
# calculation, and JAC has Plasma.SahaBoltzmannScheme for it.
#
# NINE DEFECTS were fixed to get this far; the scheme had never run end to end.  Three had broken the
# computation -> simulation chain outright, all three being leftovers of the unfinished migration from the
# per-scheme Cascade.ExcitationData to the general Cascade.Data{T}:
#   1. perform stored the lines under "photoexcitation line data:" while Cascade.extractPhotoExcitationData
#      reads "photoexcitation lines:", so the simulation always aborted with "No photoexcitationData provided.";
#   2. computeSteps returned Cascade.ExcitationData (field linesE) where the rest of the module expects
#      Cascade.Data{T} (field lines);
#   3. simulateExpansionOpacities read excData.linesE and raised a FieldError.
# The other three concern the physics inputs:
#   4. the number density entering tau was HARD-CODED as  ne = 1.0  with a literal "## number density ??"
#      comment in the source, so the input density never reached tau at all and the absolute scale of every
#      opacity ever produced by this scheme was arbitrary;
#   5. property.levelPopulation was read zero times.  The population factor was  ge/g0 * exp(-omega/kT):  no
#      partition function, the TRANSITION energy in place of the excitation energy of the absorbing level, and
#      the statistical weights the wrong way round for a lower-level population.  It is now a genuine Boltzmann
#      population normalised over the levels that occur in the line list;
#   6. the binning factor for the frequency- and temperature-dependent cases returned Delta-omega/omega instead
#      of omega/Delta-omega.  Since |Delta-lambda|/lambda = |Delta-omega|/omega exactly, all three dependences
#      must use the same factor; the two wrong ones were off by (omega/Delta-omega)^2.  The temperature bin was
#      additionally not centred on its own bin centre.
# And three stopped the computation itself, each surfacing only once the previous one was out of the way:
#   7. two API drifts in generateBlocks -- Bsplines.generateOrbitalsHydrogenic called with its arguments in
#      the wrong order (MethodError on the very first block), and SelfConsistent.performSCF used as if it
#      returned a Basis when it returns a Multiplet;  PhotoExcitation.computeLinesCascade had meanwhile
#      gained a LevelSelection argument.  Here NO selection may be applied: under a Boltzmann distribution
#      every level of the lower block absorbs, and the metastable ones carry the interesting lines;
#   8. the INITIAL configurations never became cascade blocks.  Basics.generateConfigurations returns the
#      excited configurations only, so the ground configuration -- which holds ~95% of the population and
#      hence the resonance lines, the strongest contributors to any opacity -- was absent from every
#      expansion-opacity cascade ever run;
#   9. displayExpansionOpacities labelled the binning [Hartree] whatever the dependence (it is nm for a
#      wavelength dependence) and printed the bin centres in eV, so a wavelength-dependent opacity could not
#      be plotted against wavelength without converting by hand.
# Two dead fields went with them: ExpansionOpacityScheme.meanEnergyShift (read nowhere, and duplicating
# ExpansionOpacities.transitionEnergyShift, which IS used), and printTransitions, which reached only a
# println("... not yet !!").  The latter is now implemented and produces the line list of branch a.
#
# ACCURACY, HONESTLY.  Each cascade block is a single configuration, so there is no correlation between them.
# The resonance doublet comes out well -- 4451/4585 A against NIST 4078/4216 A (9% long), f = 0.83/0.40
# against 0.71/0.34 (~15% high), and the ratio f(3/2)/f(1/2) = 2.09 exactly as measured.  The 4d - 5p
# separation does NOT: 4d lands 9% too high and 5p 8% too low, and since the transition is their DIFFERENCE
# the two errors compound to 34%.  The famous 1.0 micron triplet therefore appears here near 1.5 micron.
# Its oscillator strengths disagree between the gauges by a factor of 3-4, which is the honest signal that
# these near-cancelling weak lines are not described.  That limitation is a matter of the atomic structure,
# not of the opacity machinery, and it does not affect the checks in branches c and d, which are exact.


## The same cascade computation feeds all four branches, so it is defined once.  Sr^+ is cheap enough
## (a few minutes) that re-running it per branch is preferable to a run-dated .jld filename.
function srCascade(; printTransitions=false)
    scheme = Cascade.ExpansionOpacityScheme([E1], 0.02, 0.5, 1, [Shell("5s"), Shell("4d"), Shell("5p")],
                                                                [Shell("5p"), Shell("4d"), Shell("6s"), Shell("5d"), Shell("6p")],
                                            printTransitions)
    ## Cascade.SCA(), not AverageSCA():  under AverageSCA every block after the first re-uses the orbitals of
    ## block 1 and fills the rest with HYDROGENIC orbitals of the bare Z = 38 nucleus.  For valence excitations
    ## of a 37-electron ion that is hopeless -- it placed 4d some 350 eV above 5p, where the true splitting is
    ## about 2 eV, and every line then fell outside the photon-energy window, leaving an empty line list.
    wa = Cascade.Computation(Cascade.Computation(); name="Expansion opacity cascade for Sr^+", grid=grid,
                             nuclearModel=Nuclear.Model(38.), approach=Cascade.SCA(), scheme=scheme,
                             initialConfigs=[Configuration("[Kr] 5s")] )
    println(wa)
    return( perform(wa; output=true) )
end

## Bin centres are handed to the simulation as photon ENERGIES in Hartree, while the wavelength binning itself
## is in nm -- a mixture that is easy to get wrong, so the conversion is written out here.
omegaOfLambda(lambda_A) = convertUnits("energy: from eV to atomic", 12398.42 / lambda_A)


if  true
    # Last successful:  09-Aug-2026
    #
    # a) The line list.  Computation only: which bound-bound transitions of Sr^+ exist in the window, at what
    #    wavelengths, and with what oscillator strengths.  This is the whole input to every opacity below, so
    #    it is worth looking at before trusting any kappa.  Expect the 5s - 5p resonance doublet near 4100 A
    #    with f of order 1, and the much weaker 4d - 5p triplet near 1.0 micron; the latter carries the
    #    1 micron feature of AT2017gfo only because 4d is metastable and gets populated.
    #
    #    RESULT (7 blocks, 36 steps, ~4 min):  17 bound-bound lines from 1919 to 16823 A.
    #      + resonance doublet 5s - 5p at 4451/4585 A with f = 0.827/0.396 (Coulomb), gauge agreement 7%;
    #        NIST gives 4078/4216 A and f = 0.71/0.34, and the statistical ratio f(3/2)/f(1/2) = 2.09 is
    #        reproduced exactly.  This is the pair that dominates the opacity.
    #      + 4d - 5p triplet at 15141/15502/16823 A -- the 1.0 micron feature, displaced by the 34% error in
    #        the 4d - 5p separation discussed in the header.  Gauges differ by 3-4x on these weak lines.
    #      + the two 5s - 6p lines at 1919/1929 A have f ~ 1e-4 with a factor-12 gauge spread: near-cancelling
    #        and not to be trusted individually, though they are saturated at kilonova density anyway.
    #    Dated for what it verifies -- that the scheme yields a real, physically ordered line list -- with the
    #    4d - 5p displacement stated rather than hidden.
    setDefaults("print summary: open", "zzz-Fk-a-expansion-opacity-lines.sum")
    wb    = srCascade(printTransitions=true)
    lines = wb["photoexcitation lines:"]
    println("\n*** $(length(lines)) bound-bound lines were computed.")
    lambdas = [ convertUnits("energy: from atomic to Angstrom", line.omega)  for line in lines ]
    @printf("*** Wavelength range:  %.1f  ...  %.1f  Angstrom. \n", minimum(lambdas), maximum(lambdas))
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  09-Aug-2026
    #
    # b) kappa_exp(lambda) at kilonova conditions:  t_exp = 1 day, T = 5000 K, rho = 1.5e-13 g/cm^3.
    #    The two densities are chosen CONSISTENTLY for pure Sr ejecta:  rho = n * A * m_u  with A = 87.62,
    #    m_u = 1.6605e-24 g, so n = 1.0e9 /cm^3 gives rho = 1.455e-13 g/cm^3.  The literature quotes of order
    #    1 cm^2/g for lanthanide-poor and 10 cm^2/g for lanthanide-rich ejecta at about one day.
    #
    #    RESULT:  kappa_exp peaks at 5.31e-2 cm^2/g in the 4000 A bin (the resonance doublet), with 2.4e-2 at
    #    5000 A, 2.2e-2 at 3000 A, 1.0e-2 at 2000 A and a 3-4e-2 group from 12000 to 17000 A.  Bins from 6000
    #    to 11000 A are empty: 17 lines cannot fill a spectrum, and an empty bin here means "no line", not
    #    "no opacity".  The total is one to two orders below the ~1 cm^2/g of real ejecta, which is expected --
    #    published kilonova opacities are built from line lists of 1e5 to 1e7 transitions, several ion stages
    #    and many elements, against 17 lines of one ion here.
    #
    #    CHECKED BY HAND.  Five lines fall in the 4000 A bin (3637, 3713, 3726, 4451, 4492 A); each is deeply
    #    saturated, so each contributes its full lambda/Delta-lambda.  Sum = 20.02, and the prefactor is
    #    1/(rho c t_exp) = 1/(1.455e-13 * 2.9979e10 * 86400) = 2.654e-3 cm^2/g.  Product = 5.314e-2 against
    #    5.311656e-2 computed -- four digits, with the small residue being the one line that is not quite
    #    saturated in the Coulomb gauge.
    #
    #    Coulomb and Babushkin agree to all digits in almost every bin.  That is NOT a statement about the
    #    quality of the oscillator strengths: once tau >> 1, (1 - exp(-tau)) = 1 whatever f is, and saturation
    #    hides the gauge difference entirely.  Branch c shows it reappearing in the thin limit.
    setDefaults("print summary: open", "zzz-Fk-b-expansion-opacity.sum")
    wb        = srCascade()
    #
    ## CONTIGUOUS bins: centres every 1000 A with a 100 nm = 1000 A binning, so the bins tile the range of
    ## the line list without gaps or overlap.  1000 A is a coarse bin -- 25% of lambda at 4000 A, where the
    ## literature uses a few per cent -- but with only 17 lines a finer grid is mostly empty bins, and the
    ## Eastman & Pinto sum is meaningful only where a bin actually contains lines.
    lambdas   = collect(2000.:1000.:17000.)                                                   ## [A], bin centres
    depValues = omegaOfLambda.(lambdas)
    binning   = 100.                                                                          ## [nm] = 1000 A
    prop      = Cascade.ExpansionOpacities(Basics.BoltzmannLevelPopulation(),
                                           Cascade.WavelengthOpacityDependence(binning),
                                           1.0e9, 1.455e-13, 5000., 86400., 0., depValues)
    wc        = Cascade.Simulation(Cascade.Simulation(); name="kappa_exp(lambda) for Sr^+ at 1 day", property=prop,
                                   settings=Cascade.SimulationSettings(false, false, 0.),
                                   computationData=Dict{String,Any}[ Dict{String,Any}("results" => wb) ] )
    println(wc)
    wd = perform(wc; output=true)
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  09-Aug-2026
    #
    # c) THE MAIN INTERNAL CHECK: the two limits of (1 - exp(-tau)).  No external number is needed, because
    #    both limits are exact.
    #
    #      optically THIN  (tau << 1):   1 - exp(-tau) -> tau  is proportional to n, and n is proportional to
    #                                    rho, so the rho in the prefactor cancels:  kappa is INDEPENDENT of rho.
    #      optically THICK (tau >> 1):   1 - exp(-tau) -> 1, so kappa is proportional to 1/rho.
    #
    #    The pair (n, rho) is scaled together, as it must be for one and the same material.  Sr^+ at kilonova
    #    densities is deeply thick (tau of order 1e7 for the resonance lines), so the thin limit is only
    #    reached far below any astrophysical density -- which is the point: the check is on the formula, not
    #    on the conditions.  A tenfold drop in density should leave kappa unchanged at the thin end and raise
    #    it tenfold at the thick end.
    #
    #    RESULT -- both limits come out exactly.  kappa [cm^2/g], Coulomb gauge:
    #
    #        scale     4000 A        15000 A       regime
    #        1e-11     9.474e+05     1.137e+03     thin
    #        1e-10     9.440e+05     1.137e+03     thin
    #        1e-9      9.111e+05     1.137e+03
    #        1e-8      6.538e+05     1.137e+03
    #        1e-6      1.545e+04     1.121e+03
    #        1e-4      5.168e+02     3.780e+02
    #        1e-2      5.312e+00     4.017e+00     thick
    #        1e0       5.312e-02     4.017e-02     thick  <- kilonova conditions of branch b
    #        1e2       5.312e-04     4.017e-04     thick
    #
    #      + THICK end:  each factor of 100 in density divides kappa by exactly 100, to all printed digits,
    #        in both bins.  kappa proportional to 1/rho, confirmed.
    #      + THIN end:   kappa stops changing.  The 15000 A bin is already flat to 5 digits by scale 1e-9;
    #        the 4000 A bin, whose lines are ~300x stronger, needs 1e-10 and then changes by only 0.36% for a
    #        further factor of ten.  kappa independent of rho, confirmed.
    #      + the crossover sits near scale 1e-4 ... 1e-6, where the two bins part company because their lines
    #        saturate at different densities -- exactly the behaviour the formula must show.
    #      + THE GAUGES SEPARATE IN THE THIN LIMIT:  at 15000 A, 1.137e3 (Coulomb) against 3.791e3
    #        (Babushkin), a factor of 3.3, whereas in the thick regime they agree to all digits.  This is the
    #        cleanest illustration of the point made in branch b: saturation conceals the quality of the
    #        oscillator strengths, and any thin-limit opacity inherits their full uncertainty.
    setDefaults("print summary: open", "zzz-Fk-c-opacity-limits.sum")
    wb        = srCascade()
    #
    lambdas   = [4000., 15000.]
    depValues = omegaOfLambda.(lambdas)
    binning   = 100.
    scale     = [1.0e-11, 1.0e-10, 1.0e-9, 1.0e-8, 1.0e-6, 1.0e-4, 1.0e-2, 1.0, 1.0e2]   ## multiplies BOTH densities
    for  sc in scale
        prop = Cascade.ExpansionOpacities(Basics.BoltzmannLevelPopulation(),
                                          Cascade.WavelengthOpacityDependence(binning),
                                          1.0e9*sc, 1.455e-13*sc, 5000., 86400., 0., depValues)
        wc   = Cascade.Simulation(Cascade.Simulation(); name="density scaling factor $sc", property=prop,
                                  settings=Cascade.SimulationSettings(false, false, 0.),
                                  computationData=Dict{String,Any}[ Dict{String,Any}("results" => wb) ] )
        perform(wc; output=true)
    end
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  09-Aug-2026
    #
    # d) The two parameters a light-curve model varies: the expansion time and the temperature.
    #
    #    t_exp enters twice -- through the 1/t prefactor and through tau, which is proportional to t.  In the
    #    thick limit the prefactor wins and kappa is proportional to 1/t; in the thin limit the two cancel and
    #    kappa does not depend on t at all.  Sr^+ here is thick, so expect the 1/t behaviour.
    #
    #    T does NOT change the line list, only how the population is distributed over the lower levels.  But
    #    that is invisible at kilonova density, where every line is saturated and (1 - exp(-tau)) = 1 no matter
    #    how small the population is: a saturated line is population-blind.  The temperature scan is therefore
    #    ALSO run in the thin limit (both densities scaled by 1e-10), where kappa is proportional to n_l and
    #    the redistribution shows.  4d lies about 1.97 eV above 5s here, so at 3000 K (kT = 0.26 eV) it is
    #    almost empty and the 1.5 micron group is negligible; by 12000 K it should have gained by orders of
    #    magnitude relative to the resonance bin.
    #
    #    RESULT (1)  t_exp, thick.  kappa halves when t doubles and falls by exactly 5 when t is multiplied
    #    by 5:  1.0623e-1, 5.3117e-2, 2.6558e-2, 1.0623e-2 cm^2/g at 0.5, 1, 2 and 5 days in the 4000 A bin.
    #    kappa proportional to 1/t to all digits, as the thick limit requires -- the t in tau is irrelevant
    #    once the line is saturated, so only the 1/t prefactor survives.
    #
    #    RESULT (2)  T, thin.  kappa [cm^2/g], Coulomb:
    #
    #        T [K]      4000 A        15000 A       15000/4000
    #         3000      9.899e+05     5.707e+01     5.8e-5
    #         5000      9.440e+05     1.137e+03     1.2e-3
    #         8000      7.752e+05     4.961e+03     6.4e-3
    #        12000      5.915e+05     8.730e+03     1.5e-2
    #        20000      4.329e+05     1.068e+04     2.5e-2
    #
    #    The total does not grow without bound: opacity MOVES from the resonance bin, which loses a factor
    #    2.3 as the ground level is depopulated, into the 4d - 5p bin, which gains a factor 187.  That is the
    #    redistribution the Boltzmann factor must produce, and it is the reason the 1 micron feature exists
    #    at all -- 4d is metastable, and only a warm enough plasma populates it.
    #
    #    CHECKED INDEPENDENTLY.  In the thin limit kappa is proportional to the population of the lower level,
    #    so the 15000 A bin (all three of its lines start from 4d) must scale as the 4d fraction
    #      [4 exp(-1.9669/kT) + 6 exp(-1.98597/kT)] / Z .
    #    Evaluated by hand with the seven lowest levels, that fraction rises from 2.36e-3 at 3000 K to 0.476
    #    at 20000 K, a factor of 201, against the computed 187 -- the residue being the levels omitted from
    #    the hand-computed Z, which are exactly the ones that matter at 20000 K.  The Boltzmann population,
    #    which this scheme previously did not have at all, is therefore right to within the accuracy of the
    #    check.
    setDefaults("print summary: open", "zzz-Fk-d-opacity-time-temperature.sum")
    wb        = srCascade()
    #
    lambdas   = [4000., 15000.]
    depValues = omegaOfLambda.(lambdas)
    binning   = 100.
    #
    println("\n***  t_exp scan at kilonova density (optically thick: expect kappa proportional to 1/t):")
    for  texp  in  [43200., 86400., 172800., 432000.]
        prop = Cascade.ExpansionOpacities(Basics.BoltzmannLevelPopulation(),
                                          Cascade.WavelengthOpacityDependence(binning),
                                          1.0e9, 1.455e-13, 5000., texp, 0., depValues)
        wc   = Cascade.Simulation(Cascade.Simulation(); name="t_exp = $texp s, T = 5000 K, thick", property=prop,
                                  settings=Cascade.SimulationSettings(false, false, 0.),
                                  computationData=Dict{String,Any}[ Dict{String,Any}("results" => wb) ] )
        perform(wc; output=true)
    end
    #
    println("\n***  Temperature scan in the THIN limit (densities scaled by 1e-10), where populations matter:")
    for  temp  in  [3000., 5000., 8000., 12000., 20000.]
        prop = Cascade.ExpansionOpacities(Basics.BoltzmannLevelPopulation(),
                                          Cascade.WavelengthOpacityDependence(binning),
                                          1.0e-1, 1.455e-23, temp, 86400., 0., depValues)
        wc   = Cascade.Simulation(Cascade.Simulation(); name="T = $temp K, thin", property=prop,
                                  settings=Cascade.SimulationSettings(false, false, 0.),
                                  computationData=Dict{String,Any}[ Dict{String,Any}("results" => wb) ] )
        perform(wc; output=true)
    end
    #
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last successful:  14-Aug-2026:
    #
    #    contributions                    Planck [cm^2/g]     Rosseland [cm^2/g]
    #    bound-bound only                 1.42252e-02         0            (2 of 8 bins empty)
    #    + scattering                     1.87975e-02         8.17795e-03
    #    + scattering + bound-free        2.80495e+01         8.18263e-03
    #
    #    READ THE LAST ROW TWICE.  Adding bound-free multiplies the PLANCK mean by a factor 1500 and moves
    #    the ROSSELAND mean by 0.06 %.  That is not a bug; it is the definition doing its work.  Bound-free
    #    opens in only two of the eight bins, so an arithmetic mean -- which is dominated by whatever is most
    #    opaque -- is transformed by it, while a harmonic mean -- dominated by whatever is most transparent --
    #    barely notices, because the other six bins are unchanged and they are what a Rosseland mean is
    #    about.  Anyone who took the pre-14-Aug-2026 code at its word, where the arithmetic sum was LABELLED
    #    Rosseland, would have concluded that bound-free raises the Rosseland opacity of Sr^+ by three orders
    #    of magnitude.  It raises it by 0.06 %.
    #
    #    Checks that hold in the table: the Planck mean is exactly additive (1.42252e-02 + 4.57236e-03 =
    #    1.87975e-02, the scattering floor being kappa_es = n_e sigma_T/rho = 4.57236e-03), the Rosseland
    #    mean is not; and kappa_R < kappa_P in every row, as harmonic <= arithmetic requires.
    #
    #    CAVEAT ON THE BOUND-FREE NUMBER, which is why 2.80495e+01 is NOT a physical Sr^+ opacity.
    #    Empirical.photoionizationCrossSection scales Kramers by a binding energy drawn from X-ray tables
    #    compiled for NEUTRAL atoms: for Sr it returns 5.573 eV, essentially the neutral Sr I -> Sr II
    #    potential of 5.695 eV, whereas the true edge of Sr^+ is the Sr II -> Sr III potential of 11.030 eV.
    #    The edge is therefore placed at about half the right energy and opens inside the 0.07 - 9.85 eV band
    #    that T = 5000 K spans, where a real Sr^+ ion cannot be ionised at all.  With the correct threshold
    #    the bound-free contribution here would be identically zero and the last row would equal the second.
    #    The cross section itself is sound -- hydrogen 1s gives 6.29992 Mb at threshold against the
    #    Gaunt-corrected Kramers value of 6.30 Mb, which is what TestFrames.testMethod_Opacities asserts.
    #    It is the THRESHOLD, not the formula, that does not know about ion stages.
    #
    # e) THE TWO MEAN OPACITIES, and what each contribution does to them.  Branches a-d all report the
    #    SPECTRAL opacity kappa(lambda); this one reduces it to a single number, which is what a
    #    radiative-transfer or light-curve model actually consumes.  Two means are available and they answer
    #    different questions:
    #      Rosseland, 1/kappa_R = int (1/kappa_nu)(dB/dT) / int (dB/dT)  -- HARMONIC, so dominated by the most
    #        TRANSPARENT frequencies, because radiation leaks through the windows between lines.  This is the
    #        mean for the optically thick, diffusion regime.
    #      Planck,    kappa_P   = int kappa_nu B / int B                 -- ARITHMETIC, so dominated by the
    #        most OPAQUE frequencies.  This is the mean for optically thin emission.
    #    They are not interchangeable and on a line spectrum they differ by orders of magnitude.
    #
    #    THE POINT OF RUNNING IT WITH THREE DIFFERENT CONTRIBUTION LISTS: a Rosseland mean of a bound-bound
    #    line list ALONE is exactly zero.  An empty bin means "no line", kappa_nu = 0 there, and a single
    #    perfectly transparent window short-circuits a harmonic mean.  That is physics, not a numerical
    #    accident, and it is why a continuum contribution is not optional for a Rosseland mean.
    setDefaults("print summary: open", "zzz-Fk-e-mean-opacities.sum")
    wb    = srCascade()
    rho   = 1.455e-13                                    ## [g/cm^3], as branch b
    nIon  = rho / (87.62 * 1.6605e-24)                   ## [1/cm^3]; A = 87.62 for Sr
    nElec = nIon                                         ## singly ionised: one free electron per ion
    #
    boundBound = Cascade.BoundBoundOpacity(Basics.BoltzmannLevelPopulation())
    scattering = Cascade.ScatteringOpacity(nElec)
    boundFree  = Cascade.BoundFreeOpacity(38.0, [ (Configuration("[Kr] 5s"), Configuration("[Kr]")) ])
    #
    for  (label, contributions)  in
            [ ("bound-bound only",             Cascade.AbstractOpacityContribution[boundBound]),
              ("+ scattering",                 Cascade.AbstractOpacityContribution[boundBound, scattering]),
              ("+ scattering + bound-free",    Cascade.AbstractOpacityContribution[boundBound, scattering, boundFree]) ]
        for  mean  in  [Cascade.PlanckMean(), Cascade.RosselandMean()]
            prop = Cascade.MeanOpacities(mean, contributions, Cascade.TemperatureOpacityDependence(0.05),
                                         [nIon], [rho], [5000.], 86400., 0.)
            wc   = Cascade.Simulation(Cascade.Simulation(); name="$label", property=prop,
                                      settings=Cascade.SimulationSettings(false, false, 0.),
                                      computationData=Dict{String,Any}[ Dict{String,Any}("results" => wb) ] )
            println("\n***  $label:")
            perform(wc; output=true)
        end
    end
    setDefaults("print summary: close", "")
    #
end
