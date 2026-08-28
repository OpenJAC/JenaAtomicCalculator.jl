# Two-photon absorption by BI-CHROMATIC photons, i.e. two beams of different frequency.
# NEW FILE, 06-Aug-2026; IMPLEMENTED 08-Aug-2026 (Phase C). `TwoPhotonAbsorptionBichromatic` had existed as a
# process type since the module was written, was offered in the docstring of the abstract type and had its own
# default constructor -- but no `-inc-` file and therefore no `computeLines` method at all, so selecting it died
# in a MethodError with nothing to say why. Until 08-Aug-2026 this file held only a scaffold that named what was
# missing; the physics below now replaces it.
# THE PHYSICS DIFFERENCE from the monochromatic case, and it is the whole point of this file:
#     omega1 = scheme.omegaLess [user units]        omega2 = (E_f - E_i) - omega1
# are FIXED and unequal, so the two photons are DISTINGUISHABLE. Three consequences follow, and each of them is
# a test rather than merely a feature:
#   (1) THE TWO TIME ORDERINGS ARE PHYSICALLY DISTINCT. Absorbing omega1 first reaches the intermediate level at
#       E_i + omega1; absorbing omega2 first reaches E_i + omega2. Both orderings contribute to the same final
#       state and must be added coherently, with their OWN denominators. In the monochromatic case the two
#       orderings coincide term by term, which is why the single-beam routine gets away with one of them.
#   (2) ODD K IS ALLOWED, where a single beam forbids it. Two photons from one beam are identical bosons, so
#       their polarization state must be symmetric and odd K cancels (blocker A2, 07-Aug-2026). Two beams of
#       different colour are not identical particles, and nothing forces that cancellation. This tests the A2
#       fix FROM THE OPPOSITE SIDE: if odd K came out suppressed here as well, the mechanism claimed in A2 --
#       exchange symmetry of identical bosons, not some accidental property of the angular algebra -- would be
#       wrong. The cancellation must switch on CONTINUOUSLY as omega1 -> omega2, and it does; see below.
#       ONE REFINEMENT, MEASURED HERE AND WORTH STATING PRECISELY: what odd K needs is an ANTISYMMETRIC
#       two-photon POLARIZATION state, not merely distinguishable photons. Two beams in the SAME pure
#       polarization -- which is all Settings can describe today, since it carries one Stokes vector -- still
#       give exactly zero for odd K, because the coherent helicity sum then contains the symmetric combination
#       only. The angular weight sum_q |sum_(lambda1 lambda2) (-lambda1)(-lambda2) 3j(1,1,K;l1,l2,q)|^2 is
#       1.3333 for K = 0, EXACTLY 0 for K = 1 and 0.5333 for K = 2, independently of any amplitude. Odd K
#       therefore shows up in the UNPOLARIZED cross section, which is an incoherent mixture over helicities and
#       does contain the antisymmetric part -- and that is where it is measured below.
#   (3) THE CROSS SECTION IS EXACTLY INVARIANT under omega1 <-> omega2 when both beams carry the same
#       polarization. Exchanging the two beams exchanges the two orderings and multiplies the rank-K amplitude
#       by the phase (-1)^(L1+L2-K), whose modulus is one; since the modulus is taken per K, every observable is
#       unchanged. This is the same free symmetry check that found the K-interference bug in the emission
#       spectrum, and it costs nothing: run the same line twice with omegaLess and (E_f - E_i) - omegaLess.
# THE EXCHANGE PHASE (-1)^(L1+L2-K) is the SAME one that blocker A1 corrected in the emission file: rewriting
# [O(mp1) (x) O(mp2)]^K as [O(mp2) (x) O(mp1)]^K carries it. It appears here for a different reason -- not to
# symmetrize identical photons, but to express the second time ordering in the same coupling order as the first.
# WHAT THIS FIXES AND WHAT IT DOES NOT. The monochromatic limit of the amplitude is unambiguous: as
# omega1 -> omega2 the two orderings become term-by-term equal, so the rank-K amplitude tends to
# (1 + (-1)^K) U, i.e. 2 U for even K and 0 for odd K, where U is exactly the single-ordering quantity the
# monochromatic file computes and squares. The bichromatic cross section therefore approaches FOUR times the
# monochromatic one, and that factor is a statement about the AMPLITUDE, which is derived, not fitted. What it
# does NOT settle by itself is the remaining flux convention -- whether a single beam should be counted with
# F^2 (ordered pairs) or F^2/2 (unordered) -- which is a photon-counting question on top of the amplitude and
# has to be decided when the absorption normalisation is finally derived. Both are recorded rather than merged
# into one number; see the note at computeTotalCsUnpolarized_2pAbsorptionBichromatic.


"""
`struct  MultiPhotonTransition.Channel_2pAbsorptionBichromatic`
    ... defines a type for a two-photon absorption channel for two beams of DIFFERENT frequency and with
        well-defined multipolarities.

    + K              ::AngularJ64             ... Rank K of the channel.
    + omega1         ::Float64                ... Energy of the photon from beam 1.
    + omega2         ::Float64                ... Energy of the photon from beam 2.
    + multipole1     ::EmMultipole            ... Multipole of the photon from beam 1.
    + multipole2     ::EmMultipole            ... Multipole of the photon from beam 2.
    + gauge          ::EmGauge                ... Gauge for dealing with the (coupled) radiation field.
    + Jsym           ::LevelSymmetry          ... Symmetry of the intermediate levels used in the summation.
    + amplitude1     ::Complex{Float64}       ... reduced amplitude for the ordering in which the photon of
                                                  beam 1 is absorbed FIRST, i.e. with denominator E_i + omega1 - E_nu.
    + amplitude2     ::Complex{Float64}       ... reduced amplitude for the ordering in which the photon of
                                                  beam 2 is absorbed FIRST, i.e. with denominator E_i + omega2 - E_nu.

    TWO AMPLITUDES, NOT ONE. This is the single structural difference from Channel_2pAbsorptionMonochromatic,
    and it is deliberate: for two distinguishable beams the two time orderings are different numbers, so storing
    them separately is what lets the exchange phase be applied where it belongs and lets the monochromatic limit
    be checked term by term rather than only in the total.
"""
struct  Channel_2pAbsorptionBichromatic
    K                ::AngularJ64
    omega1           ::Float64
    omega2           ::Float64
    multipole1       ::EmMultipole
    multipole2       ::EmMultipole
    gauge            ::EmGauge
    Jsym             ::LevelSymmetry
    amplitude1       ::Complex{Float64}
    amplitude2       ::Complex{Float64}
end


"""
`struct  MultiPhotonTransition.Line_2pAbsorptionBichromatic`
    ... defines a type for a two-photon absorption line driven by two beams of different frequency.

    + initialLevel     ::Level          ... initial-(state) level
    + finalLevel       ::Level          ... final-(state) level
    + omega1           ::Float64        ... Energy of the photon from beam 1.
    + omega2           ::Float64        ... Energy of the photon from beam 2.
    + csLinear         ::EmProperty     ... Total cross section for linearly-polarized incident light.
    + csRightCircular  ::EmProperty     ... Total cross section for right-circularly polarized incident light.
    + csUnpolarized    ::EmProperty     ... Total cross section for unpolarized incident light.
    + channels         ::Array{MultiPhotonTransition.Channel_2pAbsorptionBichromatic,1}
                                        ... List of MultiPhotonTransition.Channel_2pAbsorptionBichromatic's of this line.
"""
struct  Line_2pAbsorptionBichromatic
    initialLevel       ::Level
    finalLevel         ::Level
    omega1             ::Float64
    omega2             ::Float64
    csLinear           ::EmProperty
    csRightCircular    ::EmProperty
    csUnpolarized      ::EmProperty
    channels           ::Array{MultiPhotonTransition.Channel_2pAbsorptionBichromatic,1}
end


# `Base.show(io::IO, line::MultiPhotonTransition.Line_2pAbsorptionBichromatic)`
#   ... prepares a proper printout of the variable line::MultiPhotonTransition.Line_2pAbsorptionBichromatic.
function Base.show(io::IO, line::MultiPhotonTransition.Line_2pAbsorptionBichromatic)
    println(io, "initialLevel:      $(line.initialLevel)  ")
    println(io, "finalLevel:        $(line.finalLevel)  ")
    println(io, "omega1:            $(line.omega1)  ")
    println(io, "omega2:            $(line.omega2)  ")
    println(io, "csLinear:          $(line.csLinear)  ")
    println(io, "csRightCircular:   $(line.csRightCircular)  ")
    println(io, "csUnpolarized:     $(line.csUnpolarized)  ")
    println(io, "channels:          $(line.channels)  ")
end


"""
`MultiPhotonTransition.computeChannelAmplitudes_2pAbsorptionBichromatic(
                            line::MultiPhotonTransition.Line_2pAbsorptionBichromatic, grid::Radial.Grid,
                            settings::MultiPhotonTransition.Settings)`
    ... to compute both ordering amplitudes of every channel of the given line; a
        line::MultiPhotonTransition.Line_2pAbsorptionBichromatic is returned for which the amplitudes are now
        evaluated.
"""
function  computeChannelAmplitudes_2pAbsorptionBichromatic(line::MultiPhotonTransition.Line_2pAbsorptionBichromatic,
                                                           grid::Radial.Grid, settings::MultiPhotonTransition.Settings)
    newChannels = MultiPhotonTransition.Channel_2pAbsorptionBichromatic[]
    for channel in line.channels
        # Ordering (a): the photon of beam 1 is absorbed FIRST, so the intermediate level is reached at
        # E_i + omega1 and the second (left-standing) operator is multipole2 at omega2.
        amplitude1 = MultiPhotonTransition.computeReducedAmplitudeBichromatic(channel.K, line.finalLevel,
                            channel.omega2, channel.multipole2, channel.Jsym, channel.omega1, channel.multipole1,
                            line.initialLevel, channel.gauge, grid, settings.intermediateStates, settings.selfTolerance)
        # Ordering (b): the photon of beam 2 is absorbed first -- the two roles simply exchange.
        amplitude2 = MultiPhotonTransition.computeReducedAmplitudeBichromatic(channel.K, line.finalLevel,
                            channel.omega1, channel.multipole1, channel.Jsym, channel.omega2, channel.multipole2,
                            line.initialLevel, channel.gauge, grid, settings.intermediateStates, settings.selfTolerance)
        push!( newChannels, MultiPhotonTransition.Channel_2pAbsorptionBichromatic(channel.K, channel.omega1,
                                    channel.omega2, channel.multipole1, channel.multipole2, channel.gauge,
                                    channel.Jsym, amplitude1, amplitude2) )
    end
    line = MultiPhotonTransition.Line_2pAbsorptionBichromatic( line.initialLevel, line.finalLevel, line.omega1,
                                    line.omega2, EmProperty(0.), EmProperty(0.), EmProperty(0.), newChannels)

    return( line )
end


"""
`MultiPhotonTransition.computeLines(scheme::TwoPhotonAbsorptionBichromaticScheme, finalMultiplet::Multiplet,
                            initialMultiplet::Multiplet, grid::Radial.Grid,
                            settings::MultiPhotonTransition.Settings; output=true)`
    ... to compute the two-photon absorption lines for two beams of different frequency; a list of
        lines::Array{MultiPhotonTransition.Line_2pAbsorptionBichromatic,1} is returned.
"""
function  computeLines(scheme::TwoPhotonAbsorptionBichromaticScheme, finalMultiplet::Multiplet,
                       initialMultiplet::Multiplet, grid::Radial.Grid,
                       settings::MultiPhotonTransition.Settings; output=true)
    println("")
    printstyled("MultiPhotonTransition.computeLines(::TwoPhotonAbsorptionBichromaticScheme): The computation of amplitudes starts now ... \n",
                color=:light_green)
    printstyled("---------------------------------------------------------------------------------------------------------------------- \n",
                color=:light_green)
    println("")
    lines = MultiPhotonTransition.determineLines_2pAbsorptionBichromatic(finalMultiplet, initialMultiplet, settings)
    # Display all selected lines before the computations start
    if  settings.printBefore    MultiPhotonTransition.displayLines_2pAbsorptionBichromatic(lines)    end
    # OVERVIEW MODE: rank the intermediate levels for both beams and stop before any amplitude is formed
    if  settings.calcOverview
        MultiPhotonTransition.displayIntermediateRanking_2pAbsorptionBichromatic(stdout, lines, grid, settings)
        println("\n>>> Overview only (calcOverview = true); no amplitude was computed.")
        if  output    return( MultiPhotonTransition.Line_2pAbsorptionBichromatic[] )
        else          return( nothing )
        end
    end
    # Calculate all amplitudes and requested properties
    newLines = MultiPhotonTransition.Line_2pAbsorptionBichromatic[]
    for  line in lines
        newLine = MultiPhotonTransition.computeChannelAmplitudes_2pAbsorptionBichromatic(line, grid, settings)
        newLine = MultiPhotonTransition.computeProperties_2pAbsorptionBichromatic(newLine, grid, settings)
        push!( newLines, newLine)
    end
    # Print all results to screen
    MultiPhotonTransition.displayResults_2pAbsorptionBichromatic(stdout, settings.scheme.properties, newLines)
    MultiPhotonTransition.displayKContributions_2pAbsorptionBichromatic(stdout, newLines, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    MultiPhotonTransition.displayResults_2pAbsorptionBichromatic(iostream, settings.scheme.properties, newLines)
                        MultiPhotonTransition.displayKContributions_2pAbsorptionBichromatic(iostream, newLines, settings)     end
    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`MultiPhotonTransition.computeProperties_2pAbsorptionBichromatic(
                            line::MultiPhotonTransition.Line_2pAbsorptionBichromatic, grid::Radial.Grid,
                            settings::MultiPhotonTransition.Settings)`
    ... to compute all requested properties of the given line; a
        line::MultiPhotonTransition.Line_2pAbsorptionBichromatic is returned for which the properties are now
        evaluated.
"""
function  computeProperties_2pAbsorptionBichromatic(line::MultiPhotonTransition.Line_2pAbsorptionBichromatic,
                                                    grid::Radial.Grid, settings::MultiPhotonTransition.Settings)
    csLinear        = EmProperty(0., 0.)
    csRightCircular = EmProperty(0., 0.)
    csUnpolarized   = EmProperty(0., 0.)
    for property in settings.scheme.properties
        if      typeof(property) == TotalCsLinear
            if  Basics.UseCoulomb  in  settings.gauges
                    totalCs_Cou = MultiPhotonTransition.computeTotalCsLinear_2pAbsorptionBichromatic(line, EmGauge("Coulomb"), settings)
            else    totalCs_Cou = 0.
            end
            if  Basics.UseBabushkin  in  settings.gauges
                    totalCs_Bab = MultiPhotonTransition.computeTotalCsLinear_2pAbsorptionBichromatic(line, EmGauge("Babushkin"), settings)
            else    totalCs_Bab = 0.
            end
            csLinear    = EmProperty( totalCs_Cou,  totalCs_Bab)
        elseif  typeof(property) == TotalCsRightCircular
            if  Basics.UseCoulomb  in  settings.gauges
                    totalCs_Cou = MultiPhotonTransition.computeTotalCsRightCircular_2pAbsorptionBichromatic(line, EmGauge("Coulomb"), settings)
            else    totalCs_Cou = 0.
            end
            if  Basics.UseBabushkin  in  settings.gauges
                    totalCs_Bab = MultiPhotonTransition.computeTotalCsRightCircular_2pAbsorptionBichromatic(line, EmGauge("Babushkin"), settings)
            else    totalCs_Bab = 0.
            end
            csRightCircular = EmProperty( totalCs_Cou,  totalCs_Bab)
        elseif  typeof(property) == TotalCsUnpolarized
            if  Basics.UseCoulomb  in  settings.gauges
                    totalCs_Cou = MultiPhotonTransition.computeTotalCsUnpolarized_2pAbsorptionBichromatic(line, EmGauge("Coulomb"), settings)
            else    totalCs_Cou = 0.
            end
            if  Basics.UseBabushkin  in  settings.gauges
                    totalCs_Bab = MultiPhotonTransition.computeTotalCsUnpolarized_2pAbsorptionBichromatic(line, EmGauge("Babushkin"), settings)
            else    totalCs_Bab = 0.
            end
            csUnpolarized = EmProperty( totalCs_Cou,  totalCs_Bab)
        else
            # REFUSED RATHER THAN INVENTED. TotalAlpha0 divides by omega^3, and TotalCsDensityMatrix was derived
            # for two photons drawn from ONE beam (rho (x) rho with a single Stokes vector). Neither carries over
            # to two beams of different colour without a definition that does not yet exist -- alpha_0 would need
            # a convention for which omega to use, and the density matrix would need a Stokes vector per beam.
            @warn("MultiPhotonTransition: $(typeof(property)) is not defined for the bichromatic scheme " *
                  "(alpha_0 would need one omega, the density matrix one Stokes vector per beam); it is skipped.")
        end
    end
    line = MultiPhotonTransition.Line_2pAbsorptionBichromatic( line.initialLevel, line.finalLevel, line.omega1,
                                    line.omega2, csLinear, csRightCircular, csUnpolarized, line.channels)
    return( line )
end


"""
`MultiPhotonTransition.computeReducedAmplitudeBichromatic(K::AngularJ64, finalLevel::Level, omegaB::Float64,
                            multipoleB::EmMultipole, Jsym::LevelSymmetry, omegaA::Float64, multipoleA::EmMultipole,
                            initialLevel::Level, gauge::EmGauge, grid::Radial.Grid,
                            intermediateStates::Union{Multiplet,Array{AtomicState.GreenChannel,1}},
                            selfTolerance::Float64=1.0e-8)`
    ... to compute ONE time ordering of the reduced two-photon absorption amplitude for two beams of different
        frequency: the photon (omegaA, multipoleA) is absorbed FIRST and (omegaB, multipoleB) SECOND. An
        amplitude::Complex{Float64} is returned.

        THE ONLY DIFFERENCE from the monochromatic `computeReducedAmplitudeAbsorption` is that the two photons
        carry their own energies: the denominator belongs to the photon absorbed FIRST, E_i + omegaA - E_nu,
        while each matrix element is evaluated at the energy of ITS OWN photon. Calling this routine twice with
        the two roles exchanged gives the two orderings, which for a single beam coincide and here do not.
"""
function computeReducedAmplitudeBichromatic(K::AngularJ64, finalLevel::Level, omegaB::Float64, multipoleB::EmMultipole,
                                            Jsym::LevelSymmetry, omegaA::Float64, multipoleA::EmMultipole,
                                            initialLevel::Level, gauge::EmGauge, grid::Radial.Grid,
                                            intermediateStates::Union{Multiplet,Array{AtomicState.GreenChannel,1}},
                                            selfTolerance::Float64=1.0e-8)
    U = Complex(0.);    nuLevels = MultiPhotonTransition.intermediateLevels(intermediateStates, Jsym)
    found = length(nuLevels) > 0
    for  nuLevel in nuLevels
        # A vanishing denominator is a RESONANT intermediate level, where the non-resonant perturbative
        # expression does not apply; it is skipped rather than silently producing an enormous cross section.
        # NOTE that with two beams the near-resonant case is not an accident but the experimentally interesting
        # regime -- one tunes omega1 onto an intermediate level -- so calcOverview prints BOTH denominators.
        denom = initialLevel.energy + omegaA - nuLevel.energy
        if  abs(denom) < selfTolerance    continue    end
        U = U + PhotoEmission.amplitude(Absorption(), multipoleB, gauge, omegaB, finalLevel, nuLevel, grid,
                                        display=false, printout=false) *
                PhotoEmission.amplitude(Absorption(), multipoleA, gauge, omegaA, nuLevel, initialLevel, grid,
                                        display=false, printout=false) / denom
    end

    if    found
            U = U * AngularMomentum.Wigner_6j(initialLevel.J, finalLevel.J, K, AngularJ64(multipoleB.L),
                                              AngularJ64(multipoleA.L), Jsym.J)
    else  @warn("No intermediate level of symmetry $Jsym for U^{K, 2gamma bichromatic}; " *
                "the intermediate basis does not span this symmetry.")
    end

    return( U )
end


"""
`MultiPhotonTransition.getReducedAmplitudeBichromatic(K::AngularJ64, multipole1::EmMultipole,
                            multipole2::EmMultipole, Jsym::LevelSymmetry, gauge::EmGauge,
                            channels::Array{MultiPhotonTransition.Channel_2pAbsorptionBichromatic,1})`
    ... to get/return the FULL reduced amplitude of the given channel, i.e. the coherent sum of the two time
        orderings; an amplitude::Complex{Float64} is returned.

        THE EXCHANGE PHASE (-1)^(L1+L2-K) is applied to the second ordering. It is not a symmetrization -- the
        two photons here are distinguishable -- but the price of writing both orderings in the SAME coupling
        order [O(mp2) (x) O(mp1)]^K; interchanging two coupled tensors carries exactly this phase. It is the same
        factor that blocker A1 corrected in the emission file, where it appears for the other reason.

        THE MONOCHROMATIC LIMIT IS BUILT INTO THIS LINE. As omega1 -> omega2 the two orderings become equal term
        by term, so the sum tends to (1 + (-1)^(L1+L2-K)) times one ordering: 2 U for even K (E1E1) and exactly
        zero for odd K. The vanishing of odd K for a single beam is therefore not imposed here by a rule, as it
        is in the monochromatic routine, but emerges continuously as the two colours merge.
"""
function getReducedAmplitudeBichromatic(K::AngularJ64, multipole1::EmMultipole, multipole2::EmMultipole,
                                        Jsym::LevelSymmetry, gauge::EmGauge,
                                        channels::Array{MultiPhotonTransition.Channel_2pAbsorptionBichromatic,1})
    U = Complex(0.);    found = false
    for channel in channels
        if  K == channel.K            &&  Jsym == channel.Jsym  &&  multipole1 == channel.multipole1  &&
            multipole2 == channel.multipole2                    &&  gauge == channel.gauge
            U = channel.amplitude1 +
                (-1.0)^( multipole1.L + multipole2.L - Basics.twice(K)/2 ) * channel.amplitude2
            found = true;    break
        end
    end

    if  !found
        @warn("No stored U^{K, 2gamma bichromatic} for ($gauge, K = $K, $multipole1, $multipole2, " *
              "Jsym = $Jsym); returning zero.")
    end

    return( U )
end


"""
`MultiPhotonTransition.computeTotalCsLinear_2pAbsorptionBichromatic(
                            line::MultiPhotonTransition.Line_2pAbsorptionBichromatic, gauge::EmGauge,
                            settings::MultiPhotonTransition.Settings)`
    ... to compute the total cross section for two linearly-polarized beams. A tcs::Float64 is returned.
"""
function computeTotalCsLinear_2pAbsorptionBichromatic(line::MultiPhotonTransition.Line_2pAbsorptionBichromatic,
                                                      gauge::EmGauge, settings::MultiPhotonTransition.Settings)
    tcs  = 0.;   Klist = oplus(line.finalLevel.J, line.initialLevel.J)
    symi = LevelSymmetry(line.initialLevel.J, line.initialLevel.parity)
    symf = LevelSymmetry(line.finalLevel.J,   line.finalLevel.parity)

    for  K in Klist
        # NO ODD-K RESTRICTION IS IMPOSED, and that is the physics of this file: the monochromatic routine must
        # skip odd K because two photons from one beam are identical bosons, whereas two beams of different
        # colour are not identical particles and nothing forbids it a priori.
        #
        # FOR LINEARLY POLARIZED LIGHT IT NEVERTHELESS VANISHES, and by itself rather than by a rule. Both beams
        # carry the SAME polarization here (Settings holds one Stokes vector), so the helicity sum below is
        # coherent and symmetric, while the 3-j is antisymmetric under lambda1 <-> lambda2 for odd K with
        # L1 = L2: the angular weight is exactly zero for K = 1, whatever the amplitude. Odd K is therefore
        # visible in the UNPOLARIZED cross section, an incoherent mixture that does contain the antisymmetric
        # component -- see computeTotalCsUnpolarized_2pAbsorptionBichromatic, where it is measured.
        qList = AngularMomentum.m_values(K)
        for  q in qList
            amp = ComplexF64(0.)
            for  lambda1  in [-1, 1]
                for  lambda2  in [-1, 1]
                    for  mp1 in settings.multipoles
                        for  mp2 in settings.multipoles
                            if   mp1.electric   p1 = 1    else    p1 = 0    end
                            if   mp2.electric   p2 = 1    else    p2 = 0    end
                            symmetries = MultiPhotonTransition.allowedSymmetries_2pAbsorptionBichromatic(symf, mp2, mp1, symi)
                            for Jsym in symmetries
                                amp = amp + (1.0im)^(mp1.L - p1 + mp2.L - p2) * (-lambda1)^p1 * (-lambda2)^p2 *
                                            sqrt( (2*mp1.L + 1)*(2*mp2.L + 1) ) * (Basics.twice(K) + 1)       *
                                            AngularMomentum.Wigner_3j(mp1.L, mp2.L, K, lambda1, lambda2, q)   *
                                            MultiPhotonTransition.getReducedAmplitudeBichromatic(K, mp1, mp2, Jsym,
                                                                                                 gauge, line.channels)
                            end
                        end
                    end
                end
            end
            tcs = tcs + abs( amp )^2
        end
    end

    tcs = tcs * 8*pi^5 * Defaults.getDefaults("alpha")^2 / (Basics.twice(line.initialLevel.J) + 1) /
                (line.omega1 * line.omega2)

    return( tcs )
end


"""
`MultiPhotonTransition.computeTotalCsRightCircular_2pAbsorptionBichromatic(
                            line::MultiPhotonTransition.Line_2pAbsorptionBichromatic, gauge::EmGauge,
                            settings::MultiPhotonTransition.Settings)`
    ... to compute the total cross section for two right-circularly polarized beams. A tcs::Float64 is returned.

        THE FACTOR 2 IN THE AMPLITUDE IS A NORMALISATION, not a physical factor, and it is derived rather than
        guessed. A photon of linear polarization along x is the helicity combination c_lambda = -lambda/sqrt(2),
        while the routines above write the coherent helicity sum WITHOUT the 1/sqrt(2) -- so both the linear and
        the unpolarized cross sections carry a common factor 2 per photon, i.e. 4 in the cross section. Right
        circular light is the single term c_(+1) = 1, so to be expressed in the SAME units it must be multiplied
        by sqrt(2) per photon, i.e. by 2 in the amplitude. Only then are the three polarization cross sections
        mutually comparable, which is the entire value of computing them together.

        VERIFIED, and deliberately not on a case where zero is the answer. For H 1s -> 2s this routine returns
        exactly 0, which is correct -- J = 1/2 -> 1/2 admits K <= 1 while two photons of equal helicity require
        |q| = 2 and hence K >= 2 -- but it is also what a routine computing nothing would return. The test is
        therefore H 1s -> 3d, where K = 2 is open: right-circular comes out NON-ZERO, and the ratio to linear
        light is 1.500000 for both J_f = 3/2 and J_f = 5/2, against the value 0.8/0.5333 = 1.5 that the 3-j
        weights alone predict for a pure K = 2 channel. That fixes the factor 2 above independently of the
        amplitude and of the absolute normalisation.

        NOTE that the corresponding MONOCHROMATIC routine, computeTotalCsRightCircular, has never been filled in:
        its body is a comment reading "Need to be filled" and it returns zero for every input. Any earlier
        observation that right-circular light gives exactly 0.0000 there is therefore vacuous, whatever the
        physics of the case; it is left untouched here because this file may change only its own scheme.
"""
function computeTotalCsRightCircular_2pAbsorptionBichromatic(line::MultiPhotonTransition.Line_2pAbsorptionBichromatic,
                                                             gauge::EmGauge, settings::MultiPhotonTransition.Settings)
    tcs  = 0.;   Klist = oplus(line.finalLevel.J, line.initialLevel.J);    lambda1 = 1;    lambda2 = 1
    symi = LevelSymmetry(line.initialLevel.J, line.initialLevel.parity)
    symf = LevelSymmetry(line.finalLevel.J,   line.finalLevel.parity)

    for  K in Klist
        qList = AngularMomentum.m_values(K)
        for  q in qList
            amp = ComplexF64(0.)
            for  mp1 in settings.multipoles
                for  mp2 in settings.multipoles
                    if   mp1.electric   p1 = 1    else    p1 = 0    end
                    if   mp2.electric   p2 = 1    else    p2 = 0    end
                    symmetries = MultiPhotonTransition.allowedSymmetries_2pAbsorptionBichromatic(symf, mp2, mp1, symi)
                    for Jsym in symmetries
                        amp = amp + 2.0 * (1.0im)^(mp1.L - p1 + mp2.L - p2) * (-lambda1)^p1 * (-lambda2)^p2 *
                                    sqrt( (2*mp1.L + 1)*(2*mp2.L + 1) ) * (Basics.twice(K) + 1)             *
                                    AngularMomentum.Wigner_3j(mp1.L, mp2.L, K, lambda1, lambda2, q)         *
                                    MultiPhotonTransition.getReducedAmplitudeBichromatic(K, mp1, mp2, Jsym,
                                                                                         gauge, line.channels)
                    end
                end
            end
            tcs = tcs + abs( amp )^2
        end
    end

    tcs = tcs * 8*pi^5 * Defaults.getDefaults("alpha")^2 / (Basics.twice(line.initialLevel.J) + 1) /
                (line.omega1 * line.omega2)

    return( tcs )
end


"""
`MultiPhotonTransition.computeTotalCsUnpolarized_2pAbsorptionBichromatic(
                            line::MultiPhotonTransition.Line_2pAbsorptionBichromatic, gauge::EmGauge,
                            settings::MultiPhotonTransition.Settings)`
    ... to compute the total cross section for two unpolarized beams. A tcs::Float64 is returned.

        THE PREFACTOR is taken over from the monochromatic routine with the single generalisation
        1/omega^2 -> 1/(omega1 * omega2), which is the only replacement consistent with the monochromatic limit.
        IT IS NOT DERIVED, and neither is the monochromatic one: the absorption normalisation of this module has
        never been established, unlike the emission prefactor whose residual discrepancy is at least localised to
        the constant 6.679. The POLARIZATION RATIOS, the odd-K behaviour and the omega1 <-> omega2 invariance are
        independent of it and are what this scheme can legitimately claim today.

        ON THE SINGLE-BEAM CONVENTION, which this scheme was meant to fix. The limit omega1 -> omega2 gives
        exactly FOUR times the monochromatic cross section, and that factor decomposes into two parts that must
        not be confused: (i) a factor 4 = 2^2 because the monochromatic routine keeps only ONE of the two time
        orderings in its amplitude, which is a genuine amplitude-level deficiency and not a convention; and
        (ii) the flux counting, F_1 F_2 versus F^2 or F^2/2 for one beam, which IS a convention and is untouched
        by anything computed here. Part (i) is measured and reproducible; part (ii) still needs the normalisation
        to be derived. The monochromatic routine is deliberately NOT patched from this observation, because
        changing it would move every absorption number in the module on the strength of an argument that has not
        been carried through to the absolute scale.

        MEASURED, H 1s -> 2s, Babushkin (work/diag-bichromatic.jl):

            x = omega1/(E_f-E_i)    unpolarized     linear/unpol    odd-K fraction    bichromatic/monochromatic
                 0.1                1.607734e-09      0.211197        8.944013e-01           145.060
                 0.25               3.006394e-10      0.436600        7.817002e-01            27.126
                 0.4                7.474408e-11      1.254942        3.725292e-01             6.744
                 0.49               4.462238e-11      1.988122        5.938897e-03             4.026
                 0.499              4.433576e-11      1.999881        5.974393e-05             4.000261
                 0.4999             4.433290e-11      1.999999        5.974750e-07             4.000003

        THREE THINGS ARE VISIBLE AT ONCE, and none of them needs a normalisation. The ratio to the monochromatic
        result converges to 4.000003, the derived value. The odd-K fraction falls by exactly 100 for every factor
        10 in (0.5 - x), i.e. QUADRATICALLY in omega1 - omega2 -- as it must, since the odd-K amplitude is the
        DIFFERENCE of the two orderings and that difference is linear in the detuning. And linear/unpolarized
        approaches 2.0000, the value the monochromatic scheme gives, only once odd K has switched off; away from
        the merge point it differs precisely because odd K contributes to one of the two and not to the other.
"""
function computeTotalCsUnpolarized_2pAbsorptionBichromatic(line::MultiPhotonTransition.Line_2pAbsorptionBichromatic,
                                                           gauge::EmGauge, settings::MultiPhotonTransition.Settings)
    tcs  = 0.;   Klist = oplus(line.finalLevel.J, line.initialLevel.J)
    symi = LevelSymmetry(line.initialLevel.J, line.initialLevel.parity)
    symf = LevelSymmetry(line.finalLevel.J,   line.finalLevel.parity)

    for  K in Klist
        qList = AngularMomentum.m_values(K)
        for  q in qList
            for  lambda1  in [-1, 1]
                for  lambda2  in [-1, 1]
                    amp = ComplexF64(0.)
                    for  mp1 in settings.multipoles
                        for  mp2 in settings.multipoles
                            if   mp1.electric   p1 = 1    else    p1 = 0    end
                            if   mp2.electric   p2 = 1    else    p2 = 0    end
                            symmetries = MultiPhotonTransition.allowedSymmetries_2pAbsorptionBichromatic(symf, mp2, mp1, symi)
                            for Jsym in symmetries
                                amp = amp + (1.0im)^(mp1.L - p1 + mp2.L - p2) * (-lambda1)^p1 * (-lambda2)^p2 *
                                            sqrt( (2*mp1.L + 1)*(2*mp2.L + 1) ) * (Basics.twice(K) + 1)       *
                                            AngularMomentum.Wigner_3j(mp1.L, mp2.L, K, lambda1, lambda2, q)   *
                                            MultiPhotonTransition.getReducedAmplitudeBichromatic(K, mp1, mp2, Jsym,
                                                                                                 gauge, line.channels)
                            end
                        end
                    end
                    tcs = tcs + abs( amp )^2
                end
            end
        end
    end

    tcs = tcs * 8*pi^5 * Defaults.getDefaults("alpha")^2 / (Basics.twice(line.initialLevel.J) + 1) /
                (line.omega1 * line.omega2)

    return( tcs )
end


"""
`MultiPhotonTransition.computeCsPerK_2pAbsorptionBichromatic(
                            line::MultiPhotonTransition.Line_2pAbsorptionBichromatic, gauge::EmGauge,
                            settings::MultiPhotonTransition.Settings)`
    ... to compute the unpolarized cross section RESOLVED BY the rank K; a
        contributions::Array{Tuple{AngularJ64,Float64},1} is returned.

        THIS IS THE ACCEPTANCE TEST OF THE WHOLE SCHEME, not a display convenience. Odd K must be present for two
        colours and must vanish as omega1 -> omega2; K-resolving the cross section is the only way to see it
        happen, and it is the same technique that exposed the emission exchange-phase bug (blocker A1), where
        K = 1 was carrying 48 % of a rate in which it must vanish.
"""
function computeCsPerK_2pAbsorptionBichromatic(line::MultiPhotonTransition.Line_2pAbsorptionBichromatic,
                                               gauge::EmGauge, settings::MultiPhotonTransition.Settings)
    contributions = Tuple{AngularJ64,Float64}[];    Klist = oplus(line.finalLevel.J, line.initialLevel.J)
    symi = LevelSymmetry(line.initialLevel.J, line.initialLevel.parity)
    symf = LevelSymmetry(line.finalLevel.J,   line.finalLevel.parity)

    for  K in Klist
        tcs   = 0.;    qList = AngularMomentum.m_values(K)
        for  q in qList
            for  lambda1  in [-1, 1]
                for  lambda2  in [-1, 1]
                    amp = ComplexF64(0.)
                    for  mp1 in settings.multipoles
                        for  mp2 in settings.multipoles
                            if   mp1.electric   p1 = 1    else    p1 = 0    end
                            if   mp2.electric   p2 = 1    else    p2 = 0    end
                            symmetries = MultiPhotonTransition.allowedSymmetries_2pAbsorptionBichromatic(symf, mp2, mp1, symi)
                            for Jsym in symmetries
                                amp = amp + (1.0im)^(mp1.L - p1 + mp2.L - p2) * (-lambda1)^p1 * (-lambda2)^p2 *
                                            sqrt( (2*mp1.L + 1)*(2*mp2.L + 1) ) * (Basics.twice(K) + 1)       *
                                            AngularMomentum.Wigner_3j(mp1.L, mp2.L, K, lambda1, lambda2, q)   *
                                            MultiPhotonTransition.getReducedAmplitudeBichromatic(K, mp1, mp2, Jsym,
                                                                                                 gauge, line.channels)
                            end
                        end
                    end
                    tcs = tcs + abs( amp )^2
                end
            end
        end
        tcs = tcs * 8*pi^5 * Defaults.getDefaults("alpha")^2 / (Basics.twice(line.initialLevel.J) + 1) /
                    (line.omega1 * line.omega2)
        push!( contributions, (K, tcs) )
    end

    return( contributions )
end


"""
`MultiPhotonTransition.allowedSymmetries_2pAbsorptionBichromatic(symf::LevelSymmetry, multipole2::EmMultipole,
                            multipole1::EmMultipole, symi::LevelSymmetry)`
    ... to return those intermediate symmetries that are allowed for EITHER time ordering of the two photons; an
        Array{LevelSymmetry,1} is returned.

        BOTH ORDERINGS, hence the union. `AngularMomentum.allowedTotalSymmetries(symf, mp2, mp1, symi)` describes
        the path i -> (mp1) -> nu -> (mp2) -> f only; with two distinguishable photons the reversed path
        i -> (mp2) -> nu -> (mp1) -> f is a physically different contribution to the SAME final state and may
        admit intermediate symmetries the first does not. For equal multipoles -- the E1E1 case that is all this
        module is exercised with today -- the two sets coincide and the union changes nothing; it matters as soon
        as mixed multipoles such as E1M1 are requested.
"""
function allowedSymmetries_2pAbsorptionBichromatic(symf::LevelSymmetry, multipole2::EmMultipole,
                                                   multipole1::EmMultipole, symi::LevelSymmetry)
    symmetries = AngularMomentum.allowedTotalSymmetries(symf, multipole2, multipole1, symi)
    for  sym in AngularMomentum.allowedTotalSymmetries(symf, multipole1, multipole2, symi)
        if  !(sym in symmetries)    push!(symmetries, sym)    end
    end
    return( symmetries )
end


"""
`MultiPhotonTransition.determineChannels_2pAbsorptionBichromatic(omega1::Float64, omega2::Float64,
                            finalLevel::Level, initialLevel::Level, settings::MultiPhotonTransition.Settings)`
    ... to determine a list of MultiPhotonTransition.Channel_2pAbsorptionBichromatic for a transition from the
        initial to the final level and by taking the particular settings of this computation into account; an
        Array{MultiPhotonTransition.Channel_2pAbsorptionBichromatic,1} is returned.
"""
function determineChannels_2pAbsorptionBichromatic(omega1::Float64, omega2::Float64, finalLevel::Level,
                                                   initialLevel::Level, settings::MultiPhotonTransition.Settings)
    channels = MultiPhotonTransition.Channel_2pAbsorptionBichromatic[]
    symi     = LevelSymmetry(initialLevel.J, initialLevel.parity)
    symf     = LevelSymmetry(finalLevel.J,   finalLevel.parity)
    Klist    = oplus(symf.J, symi.J)
    for  mp1 in settings.multipoles
        for  mp2 in settings.multipoles
            symmetries = MultiPhotonTransition.allowedSymmetries_2pAbsorptionBichromatic(symf, mp2, mp1, symi)
            for  symn in symmetries
                for  gauge in settings.gauges
                    # Include further restrictions if appropriate
                    if     string(mp1)[1] == 'E' && string(mp2)[1] == 'E'  &&   gauge == Basics.UseCoulomb
                        for K in Klist  push!(channels, MultiPhotonTransition.Channel_2pAbsorptionBichromatic(K, omega1, omega2,
                                                            mp1, mp2, Basics.Coulomb, symn, 0., 0.) )     end
                    elseif string(mp1)[1] == 'E' && string(mp2)[1] == 'E'  &&   gauge == Basics.UseBabushkin
                        for K in Klist  push!(channels, MultiPhotonTransition.Channel_2pAbsorptionBichromatic(K, omega1, omega2,
                                                            mp1, mp2, Basics.Babushkin, symn, 0., 0.) )   end
                    elseif string(mp1)[1] == 'M' && string(mp2)[1] == 'M'
                        for K in Klist  push!(channels, MultiPhotonTransition.Channel_2pAbsorptionBichromatic(K, omega1, omega2,
                                                            mp1, mp2, Basics.Magnetic, symn, 0., 0.) )    end
                    end
                end
            end
        end
    end

    return( channels )
end


"""
`MultiPhotonTransition.determineLines_2pAbsorptionBichromatic(finalMultiplet::Multiplet,
                            initialMultiplet::Multiplet, settings::MultiPhotonTransition.Settings)`
    ... to determine a list of MultiPhotonTransition.Line_2pAbsorptionBichromatic's for transitions between the
        levels of the given initial- and final-state multiplets and by taking the particular selections and
        settings of this computation into account; an
        Array{MultiPhotonTransition.Line_2pAbsorptionBichromatic,1} is returned. Apart from the level
        specification, all physical properties are set to zero during this initialization process.
"""
function  determineLines_2pAbsorptionBichromatic(finalMultiplet::Multiplet, initialMultiplet::Multiplet,
                                                 settings::MultiPhotonTransition.Settings)
    lines = MultiPhotonTransition.Line_2pAbsorptionBichromatic[]
    for  iLevel  in  initialMultiplet.levels
        for  fLevel  in  finalMultiplet.levels
            if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                energy = fLevel.energy - iLevel.energy + settings.photonEnergyShift
                # omegaLess is given in the USER-SELECTED energy units, following the convention of the photon
                # energies of the other process modules, and is converted here -- once, and in the only place
                # it enters.
                omega1 = Defaults.convertUnits("energy: to atomic", settings.scheme.omegaLess)
                omega2 = energy - omega1
                # GUARDED, because an omegaLess that does not fit the transition is a user error that would
                # otherwise show up as a negative photon energy deep inside an amplitude.
                if  omega1 <= 0.  ||  omega2 <= 0.
                    error("\n\nMultiPhotonTransition: omegaLess = $(settings.scheme.omegaLess) "                *
                          TableStrings.inUnits("energy") * " does not fit the transition "                      *
                          "$(iLevel.index) -> $(fLevel.index), whose total energy is "                          *
                          "$(Defaults.convertUnits("energy: from atomic", energy)) "                            *
                          TableStrings.inUnits("energy") * ".\n"                                                *
                          ">>> Both photon energies must be positive, i.e. 0 < omegaLess < E_f - E_i.\n")
                end
                channels = MultiPhotonTransition.determineChannels_2pAbsorptionBichromatic(omega1, omega2, fLevel,
                                                                                           iLevel, settings)
                push!( lines, MultiPhotonTransition.Line_2pAbsorptionBichromatic(iLevel, fLevel, omega1, omega2,
                                            EmProperty(0., 0.), EmProperty(0., 0.), EmProperty(0., 0.), channels) )
            end
        end
    end
    return( lines )
end


"""
`MultiPhotonTransition.displayLines_2pAbsorptionBichromatic(
                            lines::Array{MultiPhotonTransition.Line_2pAbsorptionBichromatic,1})`
    ... to display a list of lines and channels that have been selected due to the prior settings. A neat table
        of all selected transitions and energies is printed but nothing is returned otherwise.
"""
function  displayLines_2pAbsorptionBichromatic(lines::Array{MultiPhotonTransition.Line_2pAbsorptionBichromatic,1})
    nx = 185
    println(" ")
    println("  Selected two-photon absorption lines for two beams of different frequency:")
    println(" ")
    println("  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=3);                         sb = sb * TableStrings.hBlank(21)
    sa = sa * TableStrings.center(12, "Energy"; na=2)
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=2)
    sa = sa * TableStrings.center(12, "omega1"; na=2)
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=2)
    sa = sa * TableStrings.center(12, "omega2"; na=4)
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.flushleft(90, "List of multipoles & intermediate level symmetries"; na=4)
    sb = sb * TableStrings.flushleft(90, "(K-rank, multipole_1, Jsym, multipole_2, gauge), ..."; na=4)
    println(sa);    println(sb);    println("  ", TableStrings.hLine(nx))
    for  line in lines
        sa   = "";   isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                     fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        energy = line.finalLevel.energy - line.initialLevel.energy
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", energy))       * "    "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.omega1))  * "    "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.omega2))  * "     "
        mpGaugeList = Tuple{AngularJ64, Basics.EmMultipole, LevelSymmetry, Basics.EmMultipole, Basics.EmGauge}[]
        for  channel in  line.channels
            push!( mpGaugeList, (channel.K, channel.multipole1, channel.Jsym, channel.multipole2, channel.gauge) )
        end
        wa = TableStrings.twoPhotonGaugeTupels(105, mpGaugeList)
        if  length(wa) > 0    sb = sa * wa[1];    println( sb )    end
        for  i = 2:length(wa)
            sb = TableStrings.hBlank( length(sa) );    sb = sb * wa[i];    println( sb )
        end
    end
    println("  ", TableStrings.hLine(nx))
    return( nothing )
end


"""
`MultiPhotonTransition.displayResults_2pAbsorptionBichromatic(stream::IO,
                            properties::Array{MultiPhotonTransition.AbstractMultiPhotonProperty,1},
                            lines::Array{MultiPhotonTransition.Line_2pAbsorptionBichromatic,1})`
    ... to display all results, energies, cross sections, etc. of the selected lines. A neat table is printed but
        nothing is returned otherwise.
"""
function  displayResults_2pAbsorptionBichromatic(stream::IO,
                            properties::Array{MultiPhotonTransition.AbstractMultiPhotonProperty,1},
                            lines::Array{MultiPhotonTransition.Line_2pAbsorptionBichromatic,1})
    nx = 88
    wx = Defaults.convertUnits("length: from atomic to cm", 1.0)^4 / Defaults.convertUnits("energy: from atomic to Ws", 1.0)
    println(stream, " ")
    println(stream, "  Two-photon absorption by two beams of different frequency (bichromatic):")
    println(stream, " ")
    println(stream, "  Cross sections [cm^4/Ws] are given for (absolute scale NOT yet verified):")
    noCs = 0  # Number of cross sections to be printed
    for property in properties
        if      typeof(property) == TotalCsLinear
            noCs = noCs + 1;   println(stream, "    + total cross sections for two linearly-polarized beams ($noCs)")
        elseif  typeof(property) == TotalCsRightCircular
            noCs = noCs + 1;   println(stream, "    + total cross sections for two right-circularly polarized " *
                                               "beams ($noCs); MUST vanish for J = 0 -> J = 0")
        elseif  typeof(property) == TotalCsUnpolarized
            noCs = noCs + 1;   println(stream, "    + total cross sections for two unpolarized beams ($noCs)")
        end
    end
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx + 34noCs))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(10, "omega1";  na=4)
    sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(10, "omega2";  na=4)
    sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=7)
    for no = 1:noCs
        sa = sa * TableStrings.center(28, "Cou -- cross section ($no) -- Bab"; na=4)
        sb = sb * TableStrings.center(28, "[cm^4/Ws]" * "         " * "[cm^4/Ws]"; na=8)
    end

    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx + 34noCs))
    for  line in lines
        sa   = "";   isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                     fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=4)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.omega1))    * "    "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.omega2))    * "        "
        for property in properties
            if      typeof(property) == TotalCsLinear
                sa = sa * @sprintf("%.4e", line.csLinear.Coulomb          * wx)   * "      "
                sa = sa * @sprintf("%.4e", line.csLinear.Babushkin        * wx)   * "          "
            elseif  typeof(property) == TotalCsRightCircular
                sa = sa * @sprintf("%.4e", line.csRightCircular.Coulomb   * wx)   * "      "
                sa = sa * @sprintf("%.4e", line.csRightCircular.Babushkin * wx)   * "          "
            elseif  typeof(property) == TotalCsUnpolarized
                sa = sa * @sprintf("%.4e", line.csUnpolarized.Coulomb     * wx)   * "      "
                sa = sa * @sprintf("%.4e", line.csUnpolarized.Babushkin   * wx)   * "          "
            end
        end
        println(stream, sa )
    end
    println(stream, "  ", TableStrings.hLine(nx + 34noCs))
    return( nothing )
end


"""
`MultiPhotonTransition.displayKContributions_2pAbsorptionBichromatic(stream::IO,
                            lines::Array{MultiPhotonTransition.Line_2pAbsorptionBichromatic,1},
                            settings::MultiPhotonTransition.Settings)`
    ... to display the unpolarized cross section resolved by the rank K, together with the fraction carried by
        ODD K. Nothing is returned.

        THE ODD-K FRACTION IS THE POINT OF THIS TABLE. For a single beam it must be exactly zero (two identical
        bosons, blocker A2); here it must be non-zero, and it must fall to zero continuously as omega1 -> omega2.
        Printing it with every run makes that the default observation rather than a special diagnostic, since it
        is the one statement this scheme can make that no absolute normalisation is needed for.
"""
function  displayKContributions_2pAbsorptionBichromatic(stream::IO,
                            lines::Array{MultiPhotonTransition.Line_2pAbsorptionBichromatic,1},
                            settings::MultiPhotonTransition.Settings)
    nx = 96
    println(stream, " ")
    println(stream, "  Unpolarized cross section resolved by the rank K (relative contributions):")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, "     i-level-f     gauge         K       fraction of the total      odd-K fraction")
    println(stream, "  ", TableStrings.hLine(nx))
    for  line in lines
        for  gauge in settings.gauges
            if      gauge == Basics.UseCoulomb      emGauge = EmGauge("Coulomb")
            elseif  gauge == Basics.UseBabushkin    emGauge = EmGauge("Babushkin")
            else                                    continue
            end
            contributions = MultiPhotonTransition.computeCsPerK_2pAbsorptionBichromatic(line, emGauge, settings)
            total = 0.;    odd = 0.
            for  (K, cs) in contributions
                total = total + cs
                if  isodd( Int(Basics.twice(K)/2) )    odd = odd + cs    end
            end
            if  total == 0.    continue    end
            for  (K, cs) in contributions
                sa = "     " * TableStrings.center(12, TableStrings.levels_if(line.initialLevel.index,
                                                                             line.finalLevel.index); na=2)
                sa = sa * TableStrings.center(12, string(emGauge); na=2) * TableStrings.center(6, string(K); na=4)
                sa = sa * @sprintf("%.6e", cs/total) * "              "
                if  Basics.twice(K) == Basics.twice(contributions[1][1])   sa = sa * @sprintf("%.6e", odd/total)   end
                println(stream, sa)
            end
        end
    end
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, ">>> Odd K is FORBIDDEN for a single beam (identical bosons) and ALLOWED here; it must " *
                    "fall to zero\n>>> continuously as omega1 -> omega2.")
    return( nothing )
end


"""
`MultiPhotonTransition.displayIntermediateRanking_2pAbsorptionBichromatic(stream::IO,
                            lines::Array{MultiPhotonTransition.Line_2pAbsorptionBichromatic,1},
                            grid::Radial.Grid, settings::MultiPhotonTransition.Settings; nShow::Int64=15)`
    ... to rank, for every selected line, the intermediate levels |nu> by their contribution to the second-order
        sum and to print them together with BOTH energy denominators; nothing is returned.

        TWO DENOMINATORS, ONE PER BEAM, and that is what makes this table more than the emission analogue. With
        two colours the near-resonant case is not an accident to be avoided but the regime an experiment actually
        works in: one tunes omega1 close to a real intermediate level to gain orders of magnitude. Seeing
        E_i + omega1 - E_nu and E_i + omega2 - E_nu side by side is what tells the user how close the chosen
        colours already are to that boundary -- and beyond it the non-resonant expression computed here no longer
        applies at all.

        READ THIS TABLE IN THE LENGTH (BABUSHKIN) GAUGE, and the gauge is printed for that reason. Measured with
        work/diag-omega-dependence.jl, for H at fixed levels and omega varied over a factor 40:

            <2s| d(omega) |nu>       omega = 0.5      2.0      6.0      20.0 eV
            nu = 2p, Babushkin        1.884e-4  7.534e-4  2.260e-3  7.533e-3      <- exactly proportional to omega
            nu = 2p, Coulomb          5.971e-6  5.970e-6  5.960e-6  5.851e-6      <- independent of omega
            nu = 3p, Babushkin        1.337e-4  5.350e-4  1.605e-3  5.347e-3
            nu = 3p, Coulomb          4.045e-4  4.045e-4  4.045e-4  4.043e-4

        The LENGTH form carries the photon energy it is given. The VELOCITY form does not: by the commutator
        identity <f|p|i> = i m (E_f - E_i) <f|r|i> it carries the LEVEL difference instead, which is the same
        thing only ON SHELL. In a second-order sum the intermediate step is off shell by construction, and the
        consequence here is severe rather than cosmetic: for H 1s -> 2s the NEARLY DEGENERATE 2p intermediate --
        physically the dominant one -- is suppressed in the Coulomb ranking by a factor 250, which is precisely
        the ratio of the level gaps (2s-3p over 2s-2p = 276). It then ranks BELOW 3p and 4p, which is backwards.
        In Babushkin, 2p correctly comes out largest.

        THIS IS NOT SPECIFIC TO THE BICHROMATIC SCHEME. Every second-order sum in this module calls the same
        PhotoEmission.amplitude, so the observation applies to the emission and monochromatic-absorption paths
        as well; it is recorded here because this is where it was found, and it is an independent reason -- with
        a mechanism -- for the long-standing observation that Babushkin is the trustworthy gauge in this module.
"""
function  displayIntermediateRanking_2pAbsorptionBichromatic(stream::IO,
                            lines::Array{MultiPhotonTransition.Line_2pAbsorptionBichromatic,1},
                            grid::Radial.Grid, settings::MultiPhotonTransition.Settings; nShow::Int64=15)
    for  line in lines
        nx = 134
        println(stream, " ")
        println(stream, "  Intermediate levels ranked by their contribution, for the line  " *
                "$(line.initialLevel.index) - $(line.finalLevel.index)  with omega1 = " *
                "$(round(Defaults.convertUnits("energy: from atomic", line.omega1), sigdigits=6)) and omega2 = " *
                "$(round(Defaults.convertUnits("energy: from atomic", line.omega2), sigdigits=6)) " *
                TableStrings.inUnits("energy") * ":")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        println(stream, "     contribution     nu   J^P (nu)     E_i + omega1 - E_nu     E_i + omega2 - E_nu " *
                TableStrings.inUnits("energy") * "     multipoles      Jsym      gauge")
        println(stream, "  ", TableStrings.hLine(nx))
        entries = NamedTuple{(:c,:idx,:jsym,:d1,:d2,:mps,:sym,:gauge),
                             Tuple{Float64,Int64,String,Float64,Float64,String,String,String}}[]
        # Deduplicated over K, as in the emission ranking: the quantity ranked here does not depend on K, so
        # without this every intermediate level would be listed once per K value.
        seen = Set{Tuple{Int64,String,String,EmGauge}}()
        for  ch in line.channels
            nuLevels = MultiPhotonTransition.intermediateLevels(settings.intermediateStates, ch.Jsym)
            for  nuLevel in nuLevels
                key = (nuLevel.index, "$(ch.multipole1),$(ch.multipole2)", string(ch.Jsym), ch.gauge)
                if  key in seen    continue    end
                push!(seen, key)
                denom1 = line.initialLevel.energy + line.omega1 - nuLevel.energy
                denom2 = line.initialLevel.energy + line.omega2 - nuLevel.energy
                if  abs(denom1) < settings.selfTolerance    continue    end
                wa = PhotoEmission.amplitude(Absorption(), ch.multipole2, ch.gauge, line.omega2, line.finalLevel,
                                             nuLevel, grid, display=false, printout=false) *
                     PhotoEmission.amplitude(Absorption(), ch.multipole1, ch.gauge, line.omega1, nuLevel,
                                             line.initialLevel, grid, display=false, printout=false) / denom1
                push!(entries, (c = abs(wa), idx = nuLevel.index,
                                jsym = string(LevelSymmetry(nuLevel.J, nuLevel.parity)),
                                d1 = Defaults.convertUnits("energy: from atomic", denom1),
                                d2 = Defaults.convertUnits("energy: from atomic", denom2),
                                mps = "$(ch.multipole1),$(ch.multipole2)", sym = string(ch.Jsym),
                                gauge = string(ch.gauge)))
            end
        end
        sort!(entries, by = e -> -e.c)
        if  length(entries) == 0
            println(stream, "     -- none: the intermediate basis spans none of the required symmetries, so " *
                    "every amplitude would be zero.")
        end
        for  (n, e) in enumerate(entries)
            if  n > nShow
                println(stream, "     ... and $(length(entries)-nShow) further contributions, all smaller.")
                break
            end
            println(stream, "     " * @sprintf("%.6e", e.c) * "   " * TableStrings.center(6, string(e.idx); na=1) *
                    TableStrings.center(10, e.jsym; na=2) * @sprintf("%+14.6e", e.d1) * "        " *
                    @sprintf("%+14.6e", e.d2) * "     " *
                    TableStrings.center(10, e.mps; na=2) * TableStrings.center(10, e.sym; na=1) *
                    TableStrings.center(10, e.gauge; na=1))
        end
        println(stream, "  ", TableStrings.hLine(nx))
        println(stream, ">>> The ranking uses the ordering in which beam 1 is absorbed first; the second " *
                "ordering carries the\n>>> other denominator and is what makes odd K possible for two colours. " *
                "READ THE BABUSHKIN ROWS: the velocity\n>>> form carries the level gap rather than the photon " *
                "energy and mis-ranks near-degenerate intermediates.")
    end
    return( nothing )
end
