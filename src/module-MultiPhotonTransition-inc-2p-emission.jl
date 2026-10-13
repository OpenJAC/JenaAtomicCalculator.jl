# Two-photon emission of initially-unpolarized atoms.
"""
`struct  MultiPhotonTransition.ReducedChannel_2pEmission`  
    ... defines a type for a two-photon emission channel for the emission of photons with well-defined 
        multipolarities.

    + K              ::AngularJ64             ... Rank K of the channel.
    + omega1         ::Float64                ... omega1.
    + omega2         ::Float64                ... omega2.
    + multipole1     ::EmMultipole            ... Multipole M1.
    + multipole2     ::EmMultipole            ... Multipole M2.
    + gauge          ::EmGauge                ... Gauge for dealing with the (coupled) radiation field.
    + Jsym           ::LevelSymmetry          ... Symmetry of the Green function channel used in the summation.
    + amplitude      ::Complex{Float64}       ... reduced two-photon emission amplitude U^(K, 2gamma, emission) (..)
                                                    associated with the given channel.
"""
struct  ReducedChannel_2pEmission
    K                ::AngularJ64 
    omega1           ::Float64
    omega2           ::Float64 
    multipole1       ::EmMultipole
    multipole2       ::EmMultipole
    gauge            ::EmGauge
    Jsym             ::LevelSymmetry
    amplitude        ::Complex{Float64}
end


"""
`struct  MultiPhotonTransition.Sharing_2pEmission`  
    ... defines a type for a 2pEmission sharing to help characterize energy sharing between the two emitted photons.

    + omega1           ::Float64         ... Energy of the emitted photon 1.
    + omega2           ::Float64         ... Energy of the emitted photon 2.
    + weight           ::Float64         ... Gauss-Lengendre weight of this sharing for energy-integrated quantities.
    + differentialRate ::EmProperty      ... differential rate of this energy sharing.
    + channels         ::Array{MultiPhotonTransition.ReducedChannel_2pEmission,1}  
                                            ... List of 2pEmission channels of this sharing.
"""
struct  Sharing_2pEmission
    omega1             ::Float64
    omega2             ::Float64
    weight             ::Float64
    differentialRate   ::EmProperty
    channels           ::Array{MultiPhotonTransition.ReducedChannel_2pEmission,1}
end


"""
`struct  MultiPhotonTransition.Line_2pEmission`  
    ... defines a type for a two-photon absorption line by monochromatic light that may include the definition of channels.

    + initialLevel     ::Level                ... initial-(state) level
    + finalLevel       ::Level                ... final-(state) level
    + totalRate        ::EmProperty           ... Total rate for the two-photon transition. 
    + sharings         ::Array{MultiPhotonTransition.Sharing_2pEmission,1}  
                                                ... List of MultiPhotonTransition.Sharing_2pEmission's of this line.
"""
struct  Line_2pEmission
    initialLevel       ::Level
    finalLevel         ::Level
    totalRate          ::EmProperty
    sharings           ::Array{MultiPhotonTransition.Sharing_2pEmission,1} 
end


# `Base.show(io::IO, line::MultiPhotonTransition.Line_2pEmission)`  
#   ... prepares a proper printout of the variable line::MultiPhotonTransition.Line_2pEmission.
function Base.show(io::IO, line::MultiPhotonTransition.Line_2pEmission) 
    println(io, "initialLevel:      $(line.initialLevel)  ")
    println(io, "finalLevel:        $(line.finalLevel)  ")
    println(io, "totalRate:         $(line.totalRate)  ")
    println(io, "sharings:          $(line.sharings)  ")
end


"""
`MultiPhotonTransition.computeAmplitudes_2pEmission(line::MultiPhotonTransition.Line_2pEmission, 
                                                        grid::Radial.Grid, settings::MultiPhotonTransition.Settings)` 
    ... to compute all amplitudes of the given line; a line::MultiPhotonTransition.Line_2pEmission 
        is returned for which the amplitudes are now evaluated.
"""
function  computeAmplitudes_2pEmission(line::MultiPhotonTransition.Line_2pEmission, 
                                        grid::Radial.Grid, settings::MultiPhotonTransition.Settings)
    newSharings = MultiPhotonTransition.Sharing_2pEmission[]
    for sharing in line.sharings
        newChannels = MultiPhotonTransition.ReducedChannel_2pEmission[]
        for ch in sharing.channels
            amplitude = MultiPhotonTransition.computeReducedAmplitudeEmission(ch.K, line.finalLevel, ch.omega2, ch.multipole2, 
                                                            ch.Jsym, ch.omega1, ch.multipole1, line.initialLevel, ch.gauge, grid,
                                                            settings.intermediateStates, settings.selfTolerance)
            push!( newChannels, MultiPhotonTransition.ReducedChannel_2pEmission(ch.K, ch.omega1, ch.omega2, ch.multipole1, ch.multipole2, 
                                                                                    ch.gauge, ch.Jsym, amplitude) )
        end
        push!( newSharings, MultiPhotonTransition.Sharing_2pEmission( sharing.omega1, sharing.omega2, sharing.weight, EmProperty(0.), newChannels) )
    end
    newLine = MultiPhotonTransition.Line_2pEmission( line.initialLevel, line.finalLevel, EmProperty(0.), newSharings)
    return( newLine )
end


"""
`MultiPhotonTransition.computeEnergyDiffCs(sharing::MultiPhotonTransition.Sharing_2pEmission, line::MultiPhotonTransition.Line_2pEmission, 
                                                gauge::EmGauge, settings::MultiPhotonTransition.Settings)`  
    ... to compute the energy differential cross section for the given line and the energy sharing omega1. A dcs::Float64 is returned.
"""
function computeEnergyDiffCs(sharing::MultiPhotonTransition.Sharing_2pEmission, line::MultiPhotonTransition.Line_2pEmission, 
                                gauge::EmGauge, settings::MultiPhotonTransition.Settings)
    omega1 = sharing.omega1;    omega2 = sharing.omega2;    dcs = 0.
    for  mp1 in settings.multipoles
        for  mp2 in settings.multipoles
            symi        = LevelSymmetry(line.initialLevel.J, line.initialLevel.parity);    symf = LevelSymmetry(line.finalLevel.J, line.finalLevel.parity) 
            symmetries  = AngularMomentum.allowedTotalSymmetries(symf, mp2, mp1, symi)
            Klist       = oplus(line.finalLevel.J, line.initialLevel.J)
            # THE SUM OVER THE INTERMEDIATE SYMMETRY Jsym IS COHERENT, and sits INSIDE the modulus (corrected
            # 03-Sep-2026). It used to be the outer loop, with abs()^2 taken per Jsym and the squares added --
            # the exact mistake the K-sum note below describes, made on the other index.
            #
            # THE TWO INDICES ARE NOT ALIKE. K is the total angular-momentum transfer between the initial and
            # final LEVEL: different K are distinguishable channels for an unpolarized atom and an
            # angle-integrated rate, so their squares add, which is what the note below establishes. Jsym is the
            # symmetry of the VIRTUAL intermediate state -- nothing observes it, it is summed over inside one
            # second-order sum, and np_1/2 and np_3/2 therefore INTERFERE.
            #
            # THE TEST THAT SETTLED IT, and it was not fitted. For H 2s -> 1s the 2E1 operator is spin-independent
            # and orbitally rank 0 (l = 0 -> 0), hence a total scalar, so K = 1 must vanish up to (Z alpha)^2 ~
            # 1e-5. That cancellation runs BETWEEN the J_nu = 1/2 and J_nu = 3/2 intermediates, so squaring them
            # separately makes it impossible. Measured at n_max = 12, central sharing:
            #
            #     K = 1 share of the rate      incoherent Jsym   21.05 %      coherent Jsym   1.2e-11
            #     total, coherent/incoherent   1.4211 in Coulomb AND 1.4211 in Babushkin, identically
            #
            # Eleven digits of cancellation in both gauges independently is not something a wrong grouping
            # produces. The residual 7 % that the 31-Aug-2026 note recorded for K = 1, and left unexplained, was
            # this. NOTE the rate RISES by 1.4211: it does not close the remaining factor ~15 against the exact
            # 8.2206 1/s, and the gauge ratio is untouched because both gauges scale by the same number.
            for  K in Klist
                wk = Complex(0.)
                for Jsym in symmetries
                # THE MODULUS IS TAKEN PER K (corrected 06-Aug-2026). The K-sum used to sit INSIDE abs(), i.e.
                # dcs += |sum_K M_K|^2 instead of sum_K |M_K|^2, and that is wrong twice over.
                #
                # PHYSICALLY: different K are different total angular-momentum transfers between the initial and
                # final level. For an initially unpolarized atom and an angle-integrated rate they are
                # distinguishable channels and cannot interfere; their squares add.
                #
                # AND IT BROKE THE omega1 <-> omega2 SYMMETRY, which is how it was found. Writing the two photon
                # orderings as A_K and B_K with phase p_K, the old expression was sum_K (A_K + p_K B_K); under
                # omega1 <-> omega2 the orderings exchange, giving sum_K (B_K + p_K A_K). Those two agree only if
                # p_K is the SAME for every K. For H 2s -> 1s, J_i = J_f = 1/2 gives Klist = {0,1} and
                # p_K = (-1)^(K+1) = -1, +1 -- so they differ, and the computed spectrum came out asymmetric by
                # 60-90 % even though the two emitted photons are indistinguishable. Per K, p_K = +-1 gives
                # |wk'| = |wk| identically, so taking the modulus first restores the symmetry exactly.
                    wk = wk +
                          MultiPhotonTransition.getReducedAmplitudeEmission(K, line.finalLevel, omega2, mp2, Jsym, omega1, mp1,
                                                                                        line.initialLevel, gauge, line.sharings) +
                            # EXCHANGE PHASE (-1)^(L1+L2-K), not (-1)^(K+J_f+J_i)  -- A1, 07-Aug-2026.
                            # The two emitted photons are identical bosons, so the amplitude must be symmetric
                            # under exchanging them. Exchanging two coupled multipoles carries
                            # (-1)^(L1+L2-K); for E1E1 that is (-1)^K. The previous factor was (-1)^(K+J_f+J_i),
                            # which for J_i = J_f = 1/2 equals -(-1)^K -- an overall sign that SWAPS which K
                            # adds its two orderings and which cancels.
                            # MEASURED CONSEQUENCE: for H 2s -> 1s the 2E1 operator is spin-independent and
                            # orbitally rank 0 (l = 0 -> 0), hence a total scalar, so K = 1 must vanish up to
                            # (Z alpha)^2 ~ 1e-5. With the old phase K = 1 carried 48 % of the rate; with this
                            # one it carries 7 %.
                            (-1.0)^( mp1.L + mp2.L - Basics.twice(K)/2 ) *
                                MultiPhotonTransition.getReducedAmplitudeEmission(K, line.finalLevel, omega1, mp1, Jsym, omega2, mp2,
                                                                                        line.initialLevel, gauge, line.sharings)
                end
                dcs = dcs + abs( wk )^2
            end
        end
    end
    # THE PREFACTOR WAS WRONG BY A CONSTANT FACTOR OF 16, FOUND AND FIXED 13-Oct-2026 (challenge S5).
    #
    # Months of measurement (Z-scans, n_max scans, box/mesh convergence, a gauge-ratio trend) had excluded a
    # grid artefact, an off-shell-omega effect and basis incompleteness one at a time, but never checked the
    # prefactor itself against an INDEPENDENT reference formula -- which is what finally settled it. Two
    # hand-algebra-free comparisons, neither resting on a hand-evaluated 6j phase (an earlier same-day attempt
    # that did rest on one was retracted after it turned out to be a sign error in that derivation, not a code
    # defect -- kept in the S5 record as a cautionary example, not repeated here):
    #   * Drake (Phys. Rev. A 34, 2871 (1986), Eqs 17/18/22) gives the H 2s->1s two-photon rate as an explicit
    #     sum over exact hydrogenic np states, with NO 6j symbols (nonrelativistic, single np per n). Summing
    #     Drake's own formula numerically (exact radial integrals via GSL's generalised Laguerre polynomials)
    #     over the SAME bound n=2..20 this module's own intermediate states use gives 10.0466 /s; this module's
    #     own already-measured, dated n_max=20 Babushkin total is 0.62789 /s -- ratio 16.00.
    #   * Independently, `DalgarnoLewis.twoPhotonAmplitude` (the complete, non-truncated chi solve, see
    #     `example-As.jl`) gives 0.514 /s against the exact 8.2206 /s -- ratio also 16.00. `example-As.jl`'s own
    #     header already named this shared prefactor as the leading suspect; this is that check, carried out.
    # This also matches a hypothesis from 31-Aug-2026 (below, kept for the record) that was tried and abandoned
    # as a coincidence because the THEN-current code (pre several fixes since) measured 12.98, not 16; with
    # those fixes in place the match is now exact.
    #
    # THE 31-Aug-2026 HISTORY, KEPT BECAUSE IT IS WHAT LED HERE: a long Z-independence/n_max-dependence
    # argument ruled OUT several wrong explanations in turn (grid, off-shell omega, box size, two measurement
    # artefacts at n_max=12) without ever finding the right one. The one candidate it floated and set aside --
    # matching JAC's one-photon normalisation `A_1 = 8 pi alpha omega/(2J_i+1)|T|^2` to the golden rule and
    # naively squaring it for two photons gives `32 pi`, i.e. a factor 16 above the `2 pi` used below -- was
    # rejected then because the measured gap (12.98) didn't match and because the reduced amplitude's own 6-j
    # coupling was assumed to invalidate the naive squaring. The squaring argument's K=0 reduction turns out
    # not to change the overall scale after all; full reasoning in `[[project_challenges_list]]`, item S5.
    dcs = dcs * 32pi * Defaults.getDefaults("alpha")^2 / (Basics.twice(line.initialLevel.J) + 1) * omega1 * omega2

    return( dcs )
end



"""
`MultiPhotonTransition.computeLines(scheme::TwoPhotonEmissionScheme, finalMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                        grid::Radial.Grid, settings::MultiPhotonTransition.Settings)` 
    ... to compute the multiphoton transition amplitudes and all properties as requested by the given settings. A list of 
        lines::Array{MultiPhotonTransition.Line_2pEmission,1} is returned.
"""
function  computeLines(scheme::TwoPhotonEmissionScheme, finalMultiplet::Multiplet, initialMultiplet::Multiplet, 
                        grid::Radial.Grid, settings::MultiPhotonTransition.Settings; output=true)
    println("")
    printstyled("MultiPhotonTransition.computeLines(::TwoPhotonEmissionScheme): The computation of amplitudes starts now ... \n", color=:light_green)
    printstyled("------------------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    lines = MultiPhotonTransition.determineLines_2pEmission(finalMultiplet, initialMultiplet, settings)
    # Display all selected lines before the computations start
    if  settings.printBefore    MultiPhotonTransition.displayLines_2pEmission(lines)    end
    # OVERVIEW MODE: rank the intermediate levels and stop before any amplitude is formed
    if  settings.calcOverview
        MultiPhotonTransition.displayIntermediateRanking(stdout, lines, grid, settings)
        println("\n>>> Overview only (calcOverview = true); no amplitude was computed.")
        if  output    return( MultiPhotonTransition.Line_2pEmission[] )
        else          return( nothing )
        end
    end
    # Calculate all amplitudes
    newLines = MultiPhotonTransition.Line_2pEmission[]
    for  line in lines
        newLine = MultiPhotonTransition.computeAmplitudes_2pEmission(line, grid, settings) 
        newLine = MultiPhotonTransition.computeProperties_2pEmission(newLine, grid, settings)
        push!( newLines, newLine)
    end
    # Print all results to screen
    MultiPhotonTransition.displayTotalRates_2pEmission(stdout, newLines, settings)
    MultiPhotonTransition.displayDifferentialRates_2pEmission(stdout, newLines, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   MultiPhotonTransition.displayTotalRates_2pEmission(iostream, newLines, settings)
                        MultiPhotonTransition.displayDifferentialRates_2pEmission(iostream, newLines, settings)     end
    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`MultiPhotonTransition.computeProperties_2pEmission(line::MultiPhotonTransition.Line_2pEmission, 
                                                        grid::Radial.Grid, settings::MultiPhotonTransition.Settings)` 
    ... to compute all properties of the given line; a line::MultiPhotonTransition.Line_2pEmission 
        is returned for which the properties are now evaluated.
"""
function  computeProperties_2pEmission(line::MultiPhotonTransition.Line_2pEmission, 
                                        grid::Radial.Grid, settings::MultiPhotonTransition.Settings)
    newSharings = MultiPhotonTransition.Sharing_2pEmission[];     totalRate = EmProperty(0.)
    for sharing in line.sharings
        diffRate = EmProperty(0.)
        if  Basics.UseCoulomb  in  settings.gauges
                dr_Cou = MultiPhotonTransition.computeEnergyDiffCs(sharing, line, EmGauge("Coulomb"), settings)
        else    dr_Cou = 0.
        end
        if  Basics.UseBabushkin  in  settings.gauges
                dr_Bab = MultiPhotonTransition.computeEnergyDiffCs(sharing, line, EmGauge("Babushkin"), settings)
        else    dr_Bab = 0.
        end
        diffRate = diffRate + EmProperty(dr_Cou, dr_Bab)
        push!( newSharings, MultiPhotonTransition.Sharing_2pEmission( sharing.omega1, sharing.omega2, sharing.weight, diffRate, sharing.channels) )
        totalRate = totalRate + sharing.weight * diffRate
    end
    # THE FACTOR 1/2 IS THE INDISTINGUISHABLE-PHOTON FACTOR, applied here explicitly (08-Aug-2026) and written
    # down rather than left inside a quadrature weight.
    #
    # `computeEnergyDiffCs` already contains BOTH time orderings in its amplitude, so the differential rate at a
    # sharing (omega1, omega2) is the rate density for "one photon at omega1 and its partner at omega2" -- and
    # the sharings run over the FULL interval [0, E], where (omega1, omega2) and (omega2, omega1) are the SAME
    # physical event. Integrating over the full range therefore counts every event twice:
    #
    #     A  =  1/2 Int_0^E domega1 (dA/domega1)  =  Int_0^(E/2) domega1 (dA/domega1)
    #
    # the two being equal because the spectrum is symmetric -- which this module verifies to all printed digits
    # (branch d of example-Dh.jl) and which is exactly the symmetry the exchange-phase bug A1 was found through.
    #
    # TOGETHER WITH THE CORRECTED QUADRATURE WEIGHTS this changes every two-photon emission TOTAL by a factor 4
    # and no differential rate at all. `Basics.determineEnergySharings` used to return weights that were larger
    # by exactly 2 than the quadrature weights of Int_0^E (the 1/2 of the interval map was missing), so the old
    # totalRate was 2 * Int, where the physical answer is (1/2) * Int. The two factors are independent -- one is
    # arithmetic, the other is physics -- and are now separated: the arithmetic lives in Basics, the statistics
    # lives here.
    #
    # CONSEQUENCE FOR THE OPEN CONSTANT, stated so that nobody compares against a stale number: the deficit
    # measured against the exact H 2s -> 1s rate was 6.679 with the old totals, hence 4 * 6.679 = 26.7 with
    # these. THE DISCREPANCY IS THEREFORE LARGER, NOT SMALLER, and no factor has been introduced anywhere to
    # make it look better. The constant remains underived and must be derived, not fitted; see the long note at
    # the prefactor in computeEnergyDiffCs.
    totalRate = 0.5 * totalRate
    # Calculate the totalRate
    newLine = MultiPhotonTransition.Line_2pEmission( line.initialLevel, line.finalLevel, totalRate, newSharings)
    return( newLine )
end


"""
`MultiPhotonTransition.computeReducedAmplitudeEmission(K::AngularJ64, finalLevel::Level, 
                                                            omega2::Float64, multipole2::EmMultipole, Jsym::LevelSymmetry, 
                                                            omega1::Float64, multipole1::EmMultipole, initialLevel::Level,
                                            gauge::EmGauge, grid::Radial.Grid, greenChannels::Array{AtomicState.GreenChannel,1})`  
    ... to compute the reduced amplitude U^{K, 2gamma emission} (K, Jf, omega2, multipole2, Jsym, omega1, multipole1, Ji) by means of the
        given Green function channels. An amplitude::Complex{Float64} is returned.
"""
function computeReducedAmplitudeEmission(K::AngularJ64, finalLevel::Level, omega2::Float64, multipole2::EmMultipole, Jsym::LevelSymmetry, 
                                                                            omega1::Float64, multipole1::EmMultipole, initialLevel::Level,
                                        gauge::EmGauge, grid::Radial.Grid,
                                        intermediateStates::Union{Multiplet,Array{AtomicState.GreenChannel,1}},
                                        selfTolerance::Float64=1.0e-8)
    U = Complex(0.);    nuLevels = MultiPhotonTransition.intermediateLevels(intermediateStates, Jsym)
    found = length(nuLevels) > 0
    for  nuLevel in nuLevels
        # Skip a resonant intermediate level: the denominator vanishes and the perturbative expression is not
        # defined there. This is the boundary between non-resonant and resonant two-photon decay, and it is
        # guarded rather than hidden -- calcOverview reports the merely-small denominators.
        # MINUS omega1, not plus (corrected 07-Aug-2026). For EMISSION the intermediate state is reached AFTER
        # the first photon has left, so it lies at E_i - omega1; the "+" here was the ABSORPTION denominator,
        # which is correct in the absorption file and wrong in this one.
        #
        # The sign produced SPURIOUS POLES. For H 2s -> 1s, E_2s + omega1 - E_3p = omega1 - 1.889 eV vanishes at
        # omega1 = 1.889 eV, and the np series gives poles accumulating at the limit omega1 = 3.40 eV. The
        # differential spectrum showed exactly that: a spike at 3.225 eV and a MINIMUM at the centre, where the
        # true spectrum is a smooth dome peaked in the middle. With the correct sign the denominator is strictly
        # negative for every bound np and no pole exists.
        #
        # It also explains three things that had been puzzling for days: why ADDING intermediate states made the
        # result worse (more spurious poles, accumulating toward the series limit), why the gauge ratio was
        # Z-INDEPENDENT in the Z-scan (a structural error, not a relativistic one), and why the total sat ~6 %
        # above 8.2206 with a stubborn factor-1.7 gauge split.
        denom = initialLevel.energy - omega1 - nuLevel.energy
        if  abs(denom) < selfTolerance    continue    end
        U = U + PhotoEmission.amplitude(Emission(), multipole2, gauge, omega2, finalLevel, nuLevel, grid,
                                        display=false, printout=false) *
                PhotoEmission.amplitude(Emission(), multipole1, gauge, omega1, nuLevel, initialLevel, grid,
                                        display=false, printout=false) / denom
    end 
    
    if    found
        wx = AngularMomentum.Wigner_6j(initialLevel.J, finalLevel.J, K, AngularJ64(multipole2.L), AngularJ64(multipole1.L), Jsym.J)
        U  = U * wx
    else  @warn("No intermediate level of symmetry $Jsym for U^{K, 2gamma emission}; " *
                "the intermediate basis does not span this symmetry.")
    end 
    
    return( U )
end


"""
`MultiPhotonTransition.getReducedAmplitudeEmission(K::AngularJ64, finalLevel::Level, 
                                                        omega2::Float64, multipole2::EmMultipole, Jsym::LevelSymmetry, 
                                                        omega1::Float64, multipole1::EmMultipole, initialLevel::Level, 
                                        gauge::EmGauge, sharings::Array{MultiPhotonTransition.Sharing_2pEmission,1})`  
    ... to get/return the reduced amplitude U^{K, 2gamma emission} (K, Jf, omega2, multipole2, Jsym, omega1, multipole1, Ji) 
        from the calculated list of channels. An amplitude::Complex{Float64} is returned.
"""
function getReducedAmplitudeEmission(K::AngularJ64, finalLevel::Level, omega2::Float64, multipole2::EmMultipole, Jsym::LevelSymmetry, 
                                                                        omega1::Float64, multipole1::EmMultipole, initialLevel::Level, 
                                                        gauge::EmGauge, sharings::Array{MultiPhotonTransition.Sharing_2pEmission,1})
    U = Complex(0.);    found = false
    for sharing in sharings
        for channel in sharing.channels
            # omega1 == channel.omega1  &&  omega2 == channel.omega2  &&  
            if  K == channel.K         &&  abs(omega1 - channel.omega1) < 1.e-10  &&   abs(omega2 - channel.omega2) < 1.e-10 && 
                Jsym == channel.Jsym   &&  multipole1 == channel.multipole1  &&  multipole2 == channel.multipole2            && 
                gauge == channel.gauge
                # (gauge == channel.gauge  ||  EmGauge("Magnetic")  == channel.gauge)
                U = channel.amplitude;    found = true;     break
            end
        end
    end 
    
    if  !found
        @warn("No stored U^{K, 2gamma emission} for ($gauge, K = $K, $omega2, $multipole2, Jsym = $Jsym, " *
              "$omega1, $multipole1); returning zero.")
    end 
    
    return( U )
end


"""
`MultiPhotonTransition.determineLines_2pEmission(finalMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                                    settings::MultiPhotonTransition.Settings)`
    ... to determine a list of MultiPhotonTransition.Line_2pEmission's for transitions between the levels from the given 
        initial- and final-state multiplets and by taking into account the particular selections and settings for this computation; 
        an Array{MultiPhotonTransition.Line_2pEmission,1} is returned. Apart from the level specification, all physical 
        properties are set to zero during this initialization process.  
"""
function  determineLines_2pEmission(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::MultiPhotonTransition.Settings)
    lines = MultiPhotonTransition.Line_2pEmission[]
    for  iLevel  in  initialMultiplet.levels
        for  fLevel  in  finalMultiplet.levels
            if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                energy   = iLevel.energy - fLevel.energy
                sharings = MultiPhotonTransition.determineSharingsAndChannels(fLevel, iLevel, energy, settings) 
                push!( lines, MultiPhotonTransition.Line_2pEmission(iLevel, fLevel, EmProperty(0., 0.,), sharings) )
            end
        end
    end
    return( lines )
end


"""
`MultiPhotonTransition.determineSharingsAndChannels(finalLevel::Level, initialLevel::Level, energy::Float64, settings::MultiPhotonTransition.Settings)`  
    ... to determine a list of MultiPhotonTransition 2pEmission Sharing's and Channel's for a transitions from the initial to 
        final level and by taking into account the particular settings of for this computation; 
        an Array{MultiPhotonTransition.Sharing_2pEmission,1} is returned.
"""
function determineSharingsAndChannels(finalLevel::Level, initialLevel::Level, energy::Float64, settings::MultiPhotonTransition.Settings)
    sharings = MultiPhotonTransition.Sharing_2pEmission[];    eSharings = Basics.determineEnergySharings(energy, settings.scheme.noSharings) 
    for  es in eSharings
        omega1    = es[1];    omega2 = es[2];    weight = es[3] 
        channels  = MultiPhotonTransition.ReducedChannel_2pEmission[];   
        symi      = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
        for  mp1 in settings.multipoles
            for  mp2 in settings.multipoles
                symmetries  = AngularMomentum.allowedTotalSymmetries(symf, mp2, mp1, symi)
                Klist       = oplus(symf.J, symi.J)
                for  symx in symmetries
                    for  gauge in settings.gauges
                        # Include further restrictions if appropriate
                        if     (string(mp1)[1] == 'E' || string(mp2)[1] == 'E')  &&   gauge == Basics.UseCoulomb
                            for K in Klist  push!(channels, ReducedChannel_2pEmission(K, omega1, omega2, mp1, mp2, Basics.Coulomb, symx, 0.) )    end 
                        elseif (string(mp1)[1] == 'E' || string(mp2)[1] == 'E')  &&   gauge == Basics.UseBabushkin    
                            for K in Klist  push!(channels, ReducedChannel_2pEmission(K, omega1, omega2, mp1, mp2, Basics.Babushkin, symx, 0.) )  end
                        elseif (string(mp1)[1] == 'M' && string(mp2)[1] == 'M')
                            for K in Klist  push!(channels, ReducedChannel_2pEmission(K, omega1, omega2, mp1, mp2, Basics.Magnetic, symx, 0.) )   end
                        end
                    end
                end
            end
        end
        push!(sharings, MultiPhotonTransition.Sharing_2pEmission(omega1, omega2, weight, EmProperty(0., 0.), channels) )
    end
    return( sharings )  
end


"""
`MultiPhotonTransition.displayDifferentialRates_2pEmission(stream::IO, lines::Array{MultiPhotonTransition.Line_2pEmission,1}, 
                                                                settings::MultiPhotonTransition.Settings)`  
    ... to display all differential rates, etc. of the selected lines. A neat table is printed but nothing is returned otherwise.
"""
function  displayDifferentialRates_2pEmission(stream::IO, lines::Array{MultiPhotonTransition.Line_2pEmission,1}, 
                                                settings::MultiPhotonTransition.Settings)
    # First, print lines and sharings
    nx = 130
    println(stream, " ")
    println(stream, "  Energy-differential rates of selected two-photon emission lines:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.flushleft(38, "Energies (all in " * TableStrings.inUnits("energy") * ")"; na=4);              
    sb = sb * TableStrings.flushleft(38, "  i -- f          omega1        omega2"; na=4)
    sa = sa * TableStrings.center(14, "Weight"; na=0);                            sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(34, "Cou -- diff. rate -- Bab"; na=3)      
    sb = sb * TableStrings.center(34, TableStrings.inUnits("rate") * "        " * 
                                        TableStrings.inUnits("rate"); na=3)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    for  line in lines
        sa  = "";      isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4) 
        energy = line.finalLevel.energy - line.initialLevel.energy
        sa = sa * @sprintf("%.5e", Defaults.convertUnits("energy: from atomic", energy)) * "    "
        for  (is, sharing)  in  enumerate(line.sharings)
            if  is == 1     sb = sa     else    sb = TableStrings.hBlank( length(sa) )    end
            sb = sb * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.omega1))                   * "    "
            sb = sb * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.omega2))                   * "    "
            sb = sb * @sprintf("%.4e",                                              sharing.weight)                    * "       "
            sb = sb * @sprintf("%.5e", Defaults.convertUnits("rate: from atomic", sharing.differentialRate.Coulomb))   * "   "
            sb = sb * @sprintf("%.5e", Defaults.convertUnits("rate: from atomic", sharing.differentialRate.Babushkin)) * "   "
            println(stream,  sb )
        end
    end
    println(stream, "  ", TableStrings.hLine(nx))
    return( nothing )
end


"""
`MultiPhotonTransition.displayLines_2pEmission(lines::Array{MultiPhotonTransition.Line_2pEmission,1})`  
    ... to display a list of lines, sharings & (reduced) channels that have been selected due to the prior settings. 
        A neat table of all selected transitions and energies is printed but nothing is returned otherwise.
"""
function  displayLines_2pEmission(lines::Array{MultiPhotonTransition.Line_2pEmission,1})
    nx = 157
    println(" ")
    println("  Selected two-photon emission lines with given photon splitting:")
    println(" ")
    println("  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.flushleft(34, "Energies (all in " * TableStrings.inUnits("energy") *")"; na=5);              
    sb = sb * TableStrings.flushleft(34, "  i -- f      omega1      omega2  "; na=5)
    sa = sa * TableStrings.flushleft(77, "List of multipoles, gauges & intermediate level symmetries"; na=4)  
    sb = sb * TableStrings.flushleft(77, "(K-rank, multipole_1, Jsym, multipole_2, gauge)           "; na=4)
    println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
    for  line in lines
        sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        energy = line.initialLevel.energy - line.finalLevel.energy
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", energy)) * "  "
        for  sharing  in  line.sharings
            sb =      @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.omega1))   * "  "
            sb = sb * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", sharing.omega2))   * "     "
            mpGaugeList = Tuple{AngularJ64, Basics.EmMultipole, LevelSymmetry, Basics.EmMultipole, Basics.EmGauge}[]
            for  channel in  sharing.channels
                push!( mpGaugeList, (channel.K, channel.multipole1, channel.Jsym, channel.multipole2, channel.gauge) )
            end
            wa = TableStrings.twoPhotonGaugeTupels(75, mpGaugeList)
            sc = sa * sb * wa[1];    println( sc )  
            for  i = 2:length(wa)
                sc = TableStrings.hBlank( length(sa*sb) ) * wa[i];    println( sc )
            end
        end
    end
    println("  ", TableStrings.hLine(nx))
    return( nothing )
end


"""
`MultiPhotonTransition.displayTotalRates_2pEmission(stream::IO, lines::Array{MultiPhotonTransition.Line_2pEmission,1}, 
                                                        settings::MultiPhotonTransition.Settings)`  
    ... to display all total rates, etc. of the selected lines. A neat table is printed but nothing is returned otherwise.
"""
function  displayTotalRates_2pEmission(stream::IO, lines::Array{MultiPhotonTransition.Line_2pEmission,1}, settings::MultiPhotonTransition.Settings)
    # First, print lines and sharings
    nx = 88
    println(stream, " ")
    println(stream, "  Total rates of selected two-photon emission lines:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=2);                         sb = sb * TableStrings.hBlank(21)
    sa = sa * TableStrings.center(12, "i--Energy--f"; na=3)               
    sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=3)
    sa = sa * TableStrings.center(30, "Cou --   total rate   -- Bab"; na=4)      
    sb = sb * TableStrings.center(30, TableStrings.inUnits("rate") * "          " * 
                                        TableStrings.inUnits("rate"); na=3)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    for  line in lines
        sa  = "";      isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4) 
        energy = line.initialLevel.energy - line.finalLevel.energy
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", energy))                 * "      "
        sb = sa * @sprintf("%.5e", Defaults.convertUnits("rate: from atomic", line.totalRate.Coulomb))   * "     "
        sb = sb * @sprintf("%.5e", Defaults.convertUnits("rate: from atomic", line.totalRate.Babushkin)) * "   "
        println(stream,  sb )
    end
    println(stream, "  ", TableStrings.hLine(nx))
    return( nothing )
end


"""
`MultiPhotonTransition.displayIntermediateRanking(stream::IO, lines::Array{MultiPhotonTransition.Line_2pEmission,1},
                            grid::Radial.Grid, settings::MultiPhotonTransition.Settings; nShow::Int64=15)`
    ... to rank, for every selected line, the intermediate levels |nu> by their contribution to the second-order
        sum, and to print them together with their energy denominators; nothing is returned.

        THIS IS WHAT calcOverview EXISTS FOR. The intermediate sum is both the expensive part of a two-photon
        calculation and the delicate one, and three questions about it cannot be answered from a configuration
        list -- only from the numbers:

          * WHICH intermediate levels actually carry the transition. Usually a handful dominate and the rest are
            noise, so `intermediateStates` can be chosen from evidence instead of guessed.
          * WHETHER any denominator is near-resonant. A small (E_i + omega1 - E_nu) makes one term dominate and
            puts the whole non-resonant perturbative treatment in doubt; that is the physical boundary between
            non-resonant and resonant two-photon decay. `selfTolerance` removes the exactly singular term, but a
            merely SMALL denominator is far more dangerous, because it produces a large finite number that looks
            like a result.
          * WHETHER the basis is big enough. If the largest contributions sit on the highest levels included,
            the sum is truncated too early and the answer is not converged, however smooth it looks.

        Only the first sharing of each line is scanned: the ranking is a diagnostic of the intermediate basis,
        which barely changes across sharings, and scanning all of them would defeat the purpose of a cheap pass.
"""
function  displayIntermediateRanking(stream::IO, lines::Array{MultiPhotonTransition.Line_2pEmission,1},
                                     grid::Radial.Grid, settings::MultiPhotonTransition.Settings; nShow::Int64=15)
    for  line in lines
        if  length(line.sharings) == 0    continue    end
        sharing = line.sharings[1];    omega1 = sharing.omega1;    omega2 = sharing.omega2
        nx = 124
        println(stream, " ")
        println(stream, "  Intermediate levels ranked by their contribution, for the line  " *
                "$(line.initialLevel.index) - $(line.finalLevel.index)  at the first sharing " *
                "(omega1 = $(round(Defaults.convertUnits("energy: from atomic", omega1), sigdigits=6)) " *
                TableStrings.inUnits("energy") * "):")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        println(stream, "     contribution     nu   J^P (nu)     E_i - omega1 - E_nu " *
                TableStrings.inUnits("energy") * "      multipoles      Jsym")
        println(stream, "  ", TableStrings.hLine(nx))
        entries = NamedTuple{(:c,:idx,:jsym,:dE,:mps,:sym),Tuple{Float64,Int64,String,Float64,String,String}}[]
        # DEDUPLICATE OVER K (07-Aug-2026): sharing.channels carries one entry per K value, but the quantity
        # ranked here -- <f|O|nu><nu|O|i>/denominator -- does not depend on K at all, so every intermediate
        # level was being listed once per K. For J_i = J_f = 1/2 that meant each line printed twice, which
        # reads as though two distinct channels contributed where there is only one.
        seen = Set{Tuple{Int64,String,String,EmGauge}}()
        for  ch in sharing.channels
            nuLevels = MultiPhotonTransition.intermediateLevels(settings.intermediateStates, ch.Jsym)
            for  nuLevel in nuLevels
                key = (nuLevel.index, "$(ch.multipole1),$(ch.multipole2)", string(ch.Jsym), ch.gauge)
                if  key in seen    continue    end
                push!(seen, key)
                denom = line.initialLevel.energy - omega1 - nuLevel.energy   ## emission: see computeReducedAmplitudeEmission
                if  abs(denom) < settings.selfTolerance    continue    end
                wa = PhotoEmission.amplitude(Emission(), ch.multipole2, ch.gauge, omega2, line.finalLevel, nuLevel,
                                             grid, display=false, printout=false) *
                     PhotoEmission.amplitude(Emission(), ch.multipole1, ch.gauge, omega1, nuLevel,
                                             line.initialLevel, grid, display=false, printout=false) / denom
                push!(entries, (c = abs(wa), idx = nuLevel.index,
                                jsym = string(LevelSymmetry(nuLevel.J, nuLevel.parity)),
                                dE = Defaults.convertUnits("energy: from atomic", denom),
                                mps = "$(ch.multipole1),$(ch.multipole2)", sym = string(ch.Jsym)))
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
                    TableStrings.center(10, e.jsym; na=2) * @sprintf("%+14.6e", e.dE) * "        " *
                    TableStrings.center(10, e.mps; na=2) * TableStrings.center(10, e.sym; na=1))
        end
        println(stream, "  ", TableStrings.hLine(nx))
        println(stream, ">>> The rate scales as the SQUARE of the summed amplitude; a near-resonant denominator " *
                "makes one term dominate\n>>> and puts the non-resonant treatment itself in doubt.")
    end
    return( nothing )
end
