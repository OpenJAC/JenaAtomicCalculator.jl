# Two-photon absorption by monochromatic and equally-polarized photons, usually from the same beam.
"""
`struct  MultiPhotonTransition.Channel_2pAbsorptionMonochromatic`  
    ... defines a type for a two-photon absorption channel for the absorption of monochromatic light with well-defined 
        multipolarities.

    + K              ::AngularJ64             ... Rank K of the channel.
    + omega          ::Float64                ... omega.
    + multipole1     ::EmMultipole            ... Multipole M1.
    + multipole2     ::EmMultipole            ... Multipole M2.
    + gauge          ::EmGauge                ... Gauge for dealing with the (coupled) radiation field.
    + Jsym           ::LevelSymmetry          ... Symmetry of the Green function channel/multiplet used in the summation.
    + amplitude      ::Complex{Float64}       ... reduced two-photon absorption amplitude U^(K, 2gamma, absorption) (..)
                                                    associated with the given channel.
"""
struct  Channel_2pAbsorptionMonochromatic
    K                ::AngularJ64 
    omega            ::Float64
    multipole1       ::EmMultipole
    multipole2       ::EmMultipole
    gauge            ::EmGauge
    Jsym             ::LevelSymmetry
    amplitude        ::Complex{Float64}
end


"""
`struct  MultiPhotonTransition.Line_2pAbsorptionMonochromatic`  
    ... defines a type for a two-photon absorption line by monochromatic light that may include the definition of channels.

    + initialLevel     ::Level          ... initial-(state) level
    + finalLevel       ::Level          ... final-(state) level
    + omega            ::Float64        ... Energy of the incoming photons.
    + alpha0           ::EmProperty     ... Two-photon absorption parameter alpha_0 [often in cm^4 / Ws]
    + csLinear         ::EmProperty     ... Total cross section for linearly-polarized incident light.
    + csRightCircular  ::EmProperty     ... Total cross section for right-circularly polarized incident light.
    + csUnpolarized    ::EmProperty     ... Total cross section for unpolarized incident light.
    + channels         ::Array{MultiPhotonTransition.Channel_2pAbsorptionMonochromatic,1}  
                                        ... List of MultiPhotonTransition.Channel_2pAbsorptionMonochromatic's of this line.
"""
struct  Line_2pAbsorptionMonochromatic
    initialLevel       ::Level
    finalLevel         ::Level
    omega              ::Float64
    alpha0             ::EmProperty
    csLinear           ::EmProperty
    csRightCircular    ::EmProperty
    csUnpolarized      ::EmProperty
    csDensityMatrix    ::EmProperty
    channels           ::Array{MultiPhotonTransition.Channel_2pAbsorptionMonochromatic,1}
end


# `Base.show(io::IO, line::MultiPhotonTransition.Line_2pAbsorptionMonochromatic)`  
#   ... prepares a proper printout of the variable line::MultiPhotonTransition.Line_2pAbsorptionMonochromatic.
function Base.show(io::IO, line::MultiPhotonTransition.Line_2pAbsorptionMonochromatic) 
    println(io, "initialLevel:      $(line.initialLevel)  ")
    println(io, "finalLevel:        $(line.finalLevel)  ")
    println(io, "omega:             $(line.omega)  ")
    println(io, "alpha0:            $(line.alpha0)  ")
    println(io, "csLinear:          $(line.csLinear)  ")
    println(io, "csRightCircular:   $(line.csRightCircular)  ")
    println(io, "csUnpolarized:     $(line.csUnpolarized)  ")
    println(io, "csDensityMatrix:   $(line.csDensityMatrix)  ")
    println(io, "channels:          $(line.channels)  ")
end


"""
`MultiPhotonTransition.computeChannelAmplitudes_2pAbsorptionMonochromatic(
                            line::MultiPhotonTransition.Line_2pAbsorptionMonochromatic, grid::Radial.Grid, 
                            settings::MultiPhotonTransition.Settings)` 
    ... to compute all amplitudes and properties of the given line; a line::MultiPhotonTransition.Line_2pAbsorptionMonochromatic 
        is returned for which the amplitudes and properties are now evaluated.
"""
function  computeChannelAmplitudes_2pAbsorptionMonochromatic(line::MultiPhotonTransition.Line_2pAbsorptionMonochromatic, 
                                                                grid::Radial.Grid, settings::MultiPhotonTransition.Settings)
    newChannels = MultiPhotonTransition.Channel_2pAbsorptionMonochromatic[]
    for channel in line.channels
        amplitude = MultiPhotonTransition.computeReducedAmplitudeAbsorption(channel.K, line.finalLevel, channel.multipole2, 
                            channel.Jsym, channel.omega, channel.multipole1, line.initialLevel, channel.gauge, grid, settings.intermediateStates)
        push!( newChannels, MultiPhotonTransition.Channel_2pAbsorptionMonochromatic(channel.K, channel.omega, 
                                    channel.multipole1, channel.multipole2, channel.gauge, channel.Jsym, amplitude) )
    end
    line = MultiPhotonTransition.Line_2pAbsorptionMonochromatic( line.initialLevel, line.finalLevel, line.omega, 
                                                                    EmProperty(0.), EmProperty(0.), EmProperty(0.), EmProperty(0.),
                                                                    EmProperty(0.), newChannels)
    
    return( line )
end



"""
`MultiPhotonTransition.computeLines(scheme::TwoPhotonAbsorptionScheme, finalMultiplet::Multiplet, 
                                        initialMultiplet::Multiplet, grid::Radial.Grid, settings::MultiPhotonTransition.Settings)` 
    ... to compute the multiphoton transition amplitudes and all properties as requested by the given settings. A list of 
        lines::Array{MultiPhotonTransition.Line_2pAbsorptionMonochromatic,1} is returned.
"""
function  computeLines(scheme::TwoPhotonAbsorptionScheme, finalMultiplet::Multiplet, initialMultiplet::Multiplet, 
                        grid::Radial.Grid, settings::MultiPhotonTransition.Settings; output=true)
    println("")
    printstyled("MultiPhotonTransition.computeLines(::TwoPhotonAbsorptionScheme): The computation of amplitudes starts now ... \n", color=:light_green)
    printstyled("---------------------------------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    lines = MultiPhotonTransition.determineLines_2pAbsorptionMonochromatic(finalMultiplet, initialMultiplet, settings)
    # Display all selected lines before the computations start
    if  settings.printBefore    MultiPhotonTransition.displayLines_2pAbsorptionMonochromatic(lines)    end
    # Calculate all amplitudes and requested properties
    newLines = MultiPhotonTransition.Line_2pAbsorptionMonochromatic[]
    for  line in lines
        newLine = MultiPhotonTransition.computeChannelAmplitudes_2pAbsorptionMonochromatic(line, grid, settings) 
        newLine = MultiPhotonTransition.computeProperties_2pAbsorptionMonochromatic(newLine, grid, settings) 
        push!( newLines, newLine)
    end
    # Print all results to screen
    MultiPhotonTransition.displayResults_2pAbsorptionMonochromatic(stdout, settings.scheme.properties, newLines)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    MultiPhotonTransition.displayResults_2pAbsorptionMonochromatic(iostream, settings.scheme.properties, newLines)  end
    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`MultiPhotonTransition.computeProperties_2pAbsorptionMonochromatic(
                            line::MultiPhotonTransition.Line_2pAbsorptionMonochromatic, grid::Radial.Grid, 
                            settings::MultiPhotonTransition.Settings)` 
    ... to compute all amplitudes and properties of the given line; a line::MultiPhotonTransition.Line_2pAbsorptionMonochromatic 
        is returned for which the amplitudes and properties are now evaluated.
"""
function  computeProperties_2pAbsorptionMonochromatic(line::MultiPhotonTransition.Line_2pAbsorptionMonochromatic, 
                                                        grid::Radial.Grid, settings::MultiPhotonTransition.Settings)
    # Calculate the requested cross sections, etc.
    alpha0          = EmProperty(0., 0.)
    csLinear        = EmProperty(0., 0.)
    csRightCircular = EmProperty(0., 0.)
    csUnpolarized   = EmProperty(0., 0.)
    for property in settings.scheme.properties
        if      typeof(property) == TotalAlpha0          
            if  Basics.UseCoulomb  in  settings.gauges
                    totalA0_Cou = MultiPhotonTransition.computeTotalAlpha0(line, EmGauge("Coulomb"), settings)
            else    totalA0_Cou = 0.
            end
            if  Basics.UseBabushkin  in  settings.gauges
                    totalA0_Bab = MultiPhotonTransition.computeTotalAlpha0(line, EmGauge("Babushkin"), settings)
            else    totalA0_Bab = 0.
            end
            alpha0      = EmProperty( totalA0_Cou,  totalA0_Bab)
        elseif  typeof(property) == TotalCsLinear          
            if  Basics.UseCoulomb  in  settings.gauges
                    totalCs_Cou = MultiPhotonTransition.computeTotalCsLinear(line, EmGauge("Coulomb"), settings)
            else    totalCs_Cou = 0.
            end
            if  Basics.UseBabushkin  in  settings.gauges
                    totalCs_Bab = MultiPhotonTransition.computeTotalCsLinear(line, EmGauge("Babushkin"), settings)
            else    totalCs_Bab = 0.
            end
            csLinear    = EmProperty( totalCs_Cou,  totalCs_Bab)
        elseif  typeof(property) == TotalCsRightCircular   
            if  Basics.UseCoulomb  in  settings.gauges
                    totalCs_Cou = MultiPhotonTransition.computeTotalCsRightCircular(line, EmGauge("Coulomb"), settings)
            else    totalCs_Cou = 0.
            end
            if  Basics.UseBabushkin  in  settings.gauges
                    totalCs_Bab = MultiPhotonTransition.computeTotalCsRightCircular(line, EmGauge("Babushkin"), settings)
            else    totalCs_Bab = 0.
            end
            csRightCircular = EmProperty( totalCs_Cou,  totalCs_Bab)
        elseif  typeof(property) == TotalCsUnpolarized     
            if  Basics.UseCoulomb  in  settings.gauges
                    totalCs_Cou = MultiPhotonTransition.computeTotalCsUnpolarized(line, EmGauge("Coulomb"), settings)
            else    totalCs_Cou = 0.
            end
            if  Basics.UseBabushkin  in  settings.gauges
                    totalCs_Bab = MultiPhotonTransition.computeTotalCsUnpolarized(line, EmGauge("Babushkin"), settings)
            else    totalCs_Bab = 0.
            end
            csUnpolarized = EmProperty( totalCs_Cou,  totalCs_Bab)
        end
    end
    # THE DENSITY-MATRIX CROSS SECTION for arbitrary Stokes parameters (added 07-Aug-2026). It had been listed
    # among the properties since the module was written but was never a field, never computed and never
    # displayed, so requesting it did nothing at all.
    #
    # For a J = 0 -> J = 0 transition through two E1 photons the polarization dependence closes in a simple
    # form. Only the MIXED-helicity channel contributes: two photons of equal helicity would have to deliver
    # two units of angular momentum along the beam, which a 0 -> 0 transition cannot absorb. Writing the
    # photon density matrix in the helicity basis as rho = 1/2 [[1+P3, P1-iP2], [P1+iP2, 1-P3]] and taking the
    # two photons as independent draws from the same beam (rho (x) rho, which is the correct description of two
    # photons from ONE beam) gives
    #
    #     sigma(P1,P2,P3)  =  sigma_unpolarized * (1 + P1^2 + P2^2 - P3^2)
    #
    # and this reproduces all three special cases exactly: linear (1,0,0) -> 2*sigma_unpol = sigma_linear;
    # right-circular (0,0,1) -> 0; unpolarized (0,0,0) -> sigma_unpol. Those identities are checks in
    # themselves, fixed by angular algebra and independent of any normalisation or of the wave functions.
    #
    # LIMITATION, stated rather than hidden: the closed form above is derived for J_i = J_f = 0, where only
    # K = 0 contributes. For a general transition several K contribute with different polarization weights and
    # this expression does NOT apply; it is therefore computed only when both levels have J = 0, and left at
    # zero otherwise rather than silently returning a wrong number.
    stokes = settings.stokes
    if  Basics.twice(line.initialLevel.J) == 0  &&  Basics.twice(line.finalLevel.J) == 0
        wp  = 1.0 + stokes.P1^2 + stokes.P2^2 - stokes.P3^2
        csDensityMatrix = EmProperty(csUnpolarized.Coulomb * wp, csUnpolarized.Babushkin * wp)
    else
        csDensityMatrix = EmProperty(0., 0.)
        @warn("TotalCsDensityMatrix is implemented only for J_i = J_f = 0 (K = 0 alone); returning zero for " *
              "J_i = $(line.initialLevel.J), J_f = $(line.finalLevel.J).")
    end
    line = MultiPhotonTransition.Line_2pAbsorptionMonochromatic( line.initialLevel, line.finalLevel, line.omega, 
                                                                    alpha0, csLinear, csRightCircular, csUnpolarized,
                                                                    csDensityMatrix, line.channels)
    return( line )
end


"""
`MultiPhotonTransition.computeReducedAmplitudeAbsorption(K::AngularJ64, finalLevel::Level, multipole2::EmMultipole, Jsym::LevelSymmetry, 
                                                                                omega::Float64, multipole1::EmMultipole, initialLevel::Level,
                                                                gauge::EmGauge, grid::Radial.Grid, gMultiplet::Multiplet)`  
    ... to compute the reduced amplitude U^{K, 2gamma emission} (K, Jf, multipole2, Jsym, omega, multipole1, Ji) by means of the
        given Green function/multiplet channels. An amplitude::Complex{Float64} is returned.
"""
function computeReducedAmplitudeAbsorption(K::AngularJ64, finalLevel::Level, multipole2::EmMultipole, Jsym::LevelSymmetry, 
                                                            omega::Float64, multipole1::EmMultipole, initialLevel::Level,
                                                            gauge::EmGauge, grid::Radial.Grid,
                                    intermediateStates::Union{Multiplet,Array{AtomicState.GreenChannel,1}},
                                    selfTolerance::Float64=1.0e-8)
    U = Complex(0.);    nuLevels = MultiPhotonTransition.intermediateLevels(intermediateStates, Jsym)
    found = length(nuLevels) > 0
    for  nuLevel in nuLevels
        # A vanishing denominator is a RESONANT intermediate level, where the non-resonant perturbative
        # expression does not apply; it is skipped rather than silently producing an enormous cross section.
        denom = initialLevel.energy + omega - nuLevel.energy
        if  abs(denom) < selfTolerance    continue    end
        U = U + PhotoEmission.amplitude(Absorption(), multipole2, gauge, omega, finalLevel, nuLevel, grid,
                                        display=false, printout=false) *
                PhotoEmission.amplitude(Absorption(), multipole1, gauge, omega, nuLevel, initialLevel, grid,
                                        display=false, printout=false) / denom
    end 
    
    if    found                                
            U = U * AngularMomentum.Wigner_6j(initialLevel.J, finalLevel.J, K, AngularJ64(multipole2.L), AngularJ64(multipole1.L), Jsym.J)
    else  @warn("No intermediate level of symmetry $Jsym for U^{K, 2gamma absorption}; " *
                "the intermediate basis does not span this symmetry.")
    end 
    
    return( U )
end


"""
`MultiPhotonTransition.computeTotalAlpha0(line::MultiPhotonTransition.Line_2pAbsorptionMonochromatic, 
                                            gauge::EmGauge, settings::MultiPhotonTransition.Settings)`  
    ... to compute the (total) alpha_0 parameter for the two-photon absorption line. A ta0::Float64 is returned.
"""
function computeTotalAlpha0(line::MultiPhotonTransition.Line_2pAbsorptionMonochromatic, 
                            gauge::EmGauge, settings::MultiPhotonTransition.Settings)
    ta0 = 0.;   Klist = oplus(line.finalLevel.J, line.initialLevel.J);      omega = line.omega;     amp = ComplexF64(0.)
    symi = LevelSymmetry(line.initialLevel.J, line.initialLevel.parity);    symf = LevelSymmetry(line.finalLevel.J, line.finalLevel.parity) 
    
    for  K in Klist
        for  mp1 in settings.multipoles
            for  mp2 in settings.multipoles
                if   mp1.electric   p1 = 1    else    p1 = 0    end
                if   mp2.electric   p2 = 1    else    p2 = 0    end
                symmetries  = AngularMomentum.allowedTotalSymmetries(symf, mp2, mp1, symi)
                for Jsym in symmetries
                    amp = MultiPhotonTransition.getReducedAmplitudeAbsorption(K, line.finalLevel, mp2, Jsym, omega, mp1, 
                                                                                    line.initialLevel, gauge, line.channels) 
                    ta0 = ta0 + abs( amp )^2
                end
            end
        end
    end
    
    ta0 = ta0 * 2*pi^3 / Defaults.getDefaults("alpha")^2 / omega^3  ## / (Basics.twice(line.initialLevel.J) + 1)
    ta0 = ta0 / Defaults.getDefaults("alpha")
    
    return( ta0 )
end


"""
`MultiPhotonTransition.computeTotalCsLinear(line::MultiPhotonTransition.Line_2pAbsorptionMonochromatic, 
                                                gauge::EmGauge, settings::MultiPhotonTransition.Settings)`  
    ... to compute the total cross sections for linearly-polarized incident light. A tcs::Float64 is returned.
"""
function computeTotalCsLinear(line::MultiPhotonTransition.Line_2pAbsorptionMonochromatic, 
                                gauge::EmGauge, settings::MultiPhotonTransition.Settings)
    tcs = 0.;   Klist = oplus(line.finalLevel.J, line.initialLevel.J);      omega = line.omega
    symi = LevelSymmetry(line.initialLevel.J, line.initialLevel.parity);    symf = LevelSymmetry(line.finalLevel.J, line.finalLevel.parity) 
    
    for  K in Klist
        qList = AngularMomentum.m_values(K)
        for  q in qList
            amp = ComplexF64(0.)
            for  lambda1  in [-1, 1]
                for  lambda2  in [-1, 1]
                    for  mp1 in settings.multipoles
                        for  mp2 in settings.multipoles
                            if   mp1.electric   p1 = 1    else    p1 = 0    end
                            if   mp2.electric   p2 = 1    else    p2 = 0    end
                            symmetries  = AngularMomentum.allowedTotalSymmetries(symf, mp2, mp1, symi)
                            for Jsym in symmetries
                                wa = (1.0im)^(mp1.L - p1 + mp2.L - p2) * (-lambda1)^p1 * (-lambda2)^p2 
                                wb = sqrt( (2*mp1.L + 1)*(2*mp2.L + 1) ) * (Basics.twice(K) + 1)
                                wc = AngularMomentum.Wigner_3j(mp1.L, mp2.L, K, lambda1, lambda2, q)
                                wd = MultiPhotonTransition.getReducedAmplitudeAbsorption(K, line.finalLevel, mp2, Jsym, omega, mp1, 
                                                                                            line.initialLevel, gauge, line.channels) 
                                                                                                        
                                amp = amp + (1.0im)^(mp1.L - p1 + mp2.L - p2) * (-lambda1)^p1 * (-lambda2)^p2           *
                                            sqrt( (2*mp1.L + 1)*(2*mp2.L + 1) ) * (Basics.twice(K) + 1)       *
                                            AngularMomentum.Wigner_3j(mp1.L, mp2.L, K, lambda1, lambda2, q)             *
                                            MultiPhotonTransition.getReducedAmplitudeAbsorption(K, line.finalLevel, mp2, Jsym, omega, mp1, 
                                                                                                        line.initialLevel, gauge, line.channels) 
                            end
                        end
                    end
                end
            end
            tcs = tcs + abs( amp )^2
        end
    end
    
    # println("computeTotalCsLinear: tcs = $tcs")
    # PREFACTOR CORRECTED 11-Aug-2026.  It read  8*pi^5 * alpha^2 ...  and therefore carried alpha^2 where
    # the radiation relation requires c^2 = 1/alpha^2.  Two independent arguments fix this:
    #  (i) DIMENSIONS.  Comparing this module's own two-photon emission and absorption for the SAME
    #      transition, the unconverged intermediate sum cancels and only radiation physics remains,
    #          dA/domega_1 = (g_l/g_u) * delta * omega_1^2 omega_2^2 / (pi^4 c^4),   delta = hbar*sigma,
    #      each photon contributing the mode density c*rho(omega) = omega^2/(pi^2 c^2) that also makes the
    #      one-photon Einstein relation come out right.  With the OLD prefactor the ratio of the two
    #      routines carried NO power of c at all, alpha^2 cancelling between them, while the relation
    #      demands c^-4: a cross section and a rate cannot differ by no power of alpha.
    #  (ii) INTERNALLY.  computeTotalAlpha0 in this same file already divides by alpha^2.
    # Measured before the change: emission/absorption exceeded the relation by 1.63264e7 at Z = 1 and
    # 1.63255e7 at Z = 2 -- constant to 6e-5 over a fourfold change in omega, so a pure prefactor.  It
    # decomposes exactly as c^4/4 = 8.81614e7 times 0.185188.  The c^4/4 is installed here; the residual
    # 0.185188 (0.656 per single-photon amplitude) is NOT, being still unexplained -- it most likely
    # reflects how the reduced matrix elements carry their (2J+1) weights when amplitude(::Absorption)
    # conjugates and swaps the levels.  Absolute cross sections therefore remain low by about 5.4.
    tcs = tcs * 2*pi^5 / Defaults.getDefaults("alpha")^2 / (Basics.twice(line.initialLevel.J) + 1) / omega^2
    
    return( tcs )
end


"""
`MultiPhotonTransition.computeTotalCsRightCircular(line::MultiPhotonTransition.Line_2pAbsorptionMonochromatic,
                                                        gauge::EmGauge, settings::MultiPhotonTransition.Settings)`
    ... to compute the total cross sections for right-cicularly polarized incident light. A tcs::Float64 is returned.

        WRITTEN 08-Aug-2026. Until then the body of this routine was the comment "Need to be filled" followed by
        the prefactor, so it returned EXACTLY ZERO for every input -- for every atom, every transition and every
        gauge. That matters beyond the missing feature: the observation "right-circular light vanishes
        identically", recorded as a parameter-free check for the hydrogen and magnesium branches of
        example-Dh.jl, was produced by a routine that computed nothing and was therefore VACUOUS. Both cases
        happen to be ones where zero is also the right answer, which is exactly why it went unnoticed.

        THE FACTOR 2 IS A NORMALISATION and is derived, not guessed. A photon linearly polarized along x is the
        helicity combination c_lambda = -lambda/sqrt(2), while `computeTotalCsLinear` and
        `computeTotalCsUnpolarized` write the helicity sum WITHOUT the 1/sqrt(2) -- so both carry a common factor
        2 per photon, i.e. 4 in the cross section. Right-circular light is the single term c_(+1) = 1, so it must
        be multiplied by sqrt(2) per photon, i.e. by 2 in the amplitude, to be expressed in the same units. Only
        then are the three polarization cross sections mutually comparable, which is the whole value of computing
        them from one set of amplitudes.

        VERIFIED where the answer is NOT zero, since a routine returning zero is what was being replaced: for
        H 1s -> 3d, where J_f = 3/2, 5/2 admit K = 2 and two photons of equal helicity can be absorbed, this
        gives right-circular/linear = 1.5 exactly, against 0.8/0.5333 = 1.5 from the 3-j weights of a pure K = 2
        channel. The same ratio comes out of the independently written bichromatic routine.

        ODD K cannot contribute for two photons of EQUAL helicity: with L1 = L2 the 3-j is antisymmetric under
        exchanging the two identical helicity labels, so it vanishes identically. The guard below is therefore
        redundant here -- it is kept for parity with `computeTotalCsUnpolarized`, where it is NOT redundant, and
        removing it must leave every number unchanged.
"""
function computeTotalCsRightCircular(line::MultiPhotonTransition.Line_2pAbsorptionMonochromatic,
                                        gauge::EmGauge, settings::MultiPhotonTransition.Settings)
    tcs = 0.;   Klist = oplus(line.finalLevel.J, line.initialLevel.J);      omega = line.omega
    lambda1 = 1;    lambda2 = 1     ## right-circular: both photons carry helicity +1
    symi = LevelSymmetry(line.initialLevel.J, line.initialLevel.parity);    symf = LevelSymmetry(line.finalLevel.J, line.finalLevel.parity)

    for  K in Klist
        if  isodd( Int(Basics.twice(K)/2) )    continue    end
        qList = AngularMomentum.m_values(K)
        for  q in qList
            amp = ComplexF64(0.)
            for  mp1 in settings.multipoles
                for  mp2 in settings.multipoles
                    if   mp1.electric   p1 = 1    else    p1 = 0    end
                    if   mp2.electric   p2 = 1    else    p2 = 0    end
                    symmetries  = AngularMomentum.allowedTotalSymmetries(symf, mp2, mp1, symi)
                    for Jsym in symmetries
                        amp = amp + 2.0 * (1.0im)^(mp1.L - p1 + mp2.L - p2) * (-lambda1)^p1 * (-lambda2)^p2 *
                                    sqrt( (2*mp1.L + 1)*(2*mp2.L + 1) ) * (Basics.twice(K) + 1)             *
                                    AngularMomentum.Wigner_3j(mp1.L, mp2.L, K, lambda1, lambda2, q)         *
                                    MultiPhotonTransition.getReducedAmplitudeAbsorption(K, line.finalLevel, mp2, Jsym, omega, mp1,
                                                                                                line.initialLevel, gauge, line.channels)
                    end
                end
            end
            tcs = tcs + abs( amp )^2
        end
    end

    # PREFACTOR CORRECTED 11-Aug-2026.  It read  8*pi^5 * alpha^2 ...  and therefore carried alpha^2 where
    # the radiation relation requires c^2 = 1/alpha^2.  Two independent arguments fix this:
    #  (i) DIMENSIONS.  Comparing this module's own two-photon emission and absorption for the SAME
    #      transition, the unconverged intermediate sum cancels and only radiation physics remains,
    #          dA/domega_1 = (g_l/g_u) * delta * omega_1^2 omega_2^2 / (pi^4 c^4),   delta = hbar*sigma,
    #      each photon contributing the mode density c*rho(omega) = omega^2/(pi^2 c^2) that also makes the
    #      one-photon Einstein relation come out right.  With the OLD prefactor the ratio of the two
    #      routines carried NO power of c at all, alpha^2 cancelling between them, while the relation
    #      demands c^-4: a cross section and a rate cannot differ by no power of alpha.
    #  (ii) INTERNALLY.  computeTotalAlpha0 in this same file already divides by alpha^2.
    # Measured before the change: emission/absorption exceeded the relation by 1.63264e7 at Z = 1 and
    # 1.63255e7 at Z = 2 -- constant to 6e-5 over a fourfold change in omega, so a pure prefactor.  It
    # decomposes exactly as c^4/4 = 8.81614e7 times 0.185188.  The c^4/4 is installed here; the residual
    # 0.185188 (0.656 per single-photon amplitude) is NOT, being still unexplained -- it most likely
    # reflects how the reduced matrix elements carry their (2J+1) weights when amplitude(::Absorption)
    # conjugates and swaps the levels.  Absolute cross sections therefore remain low by about 5.4.
    tcs = tcs * 2*pi^5 / Defaults.getDefaults("alpha")^2 / (Basics.twice(line.initialLevel.J) + 1) / omega^2

    return( tcs )
end


"""
`MultiPhotonTransition.computeTotalCsUnpolarized(line::MultiPhotonTransition.Line_2pAbsorptionMonochromatic, 
                                                    gauge::EmGauge, settings::MultiPhotonTransition.Settings)`  
    ... to compute the total cross sections for linearly-polarized incident light. A tcs::Float64 is returned.
"""
function computeTotalCsUnpolarized(line::MultiPhotonTransition.Line_2pAbsorptionMonochromatic, 
                                    gauge::EmGauge, settings::MultiPhotonTransition.Settings)
    tcs = 0.;   Klist = oplus(line.finalLevel.J, line.initialLevel.J);      omega = line.omega
    symi = LevelSymmetry(line.initialLevel.J, line.initialLevel.parity);    symf = LevelSymmetry(line.finalLevel.J, line.finalLevel.parity) 
    
    for  K in Klist
        # ODD K IS FORBIDDEN FOR TWO PHOTONS FROM THE SAME BEAM (added 07-Aug-2026). They are identical bosons,
        # so the two-photon polarization state must be SYMMETRIC; exchanging the two multipoles in the 3-j
        # carries (-1)^(L1+L2+K), i.e. (-1)^K for E1E1, so odd K is antisymmetric and cannot contribute.
        #
        # LINEAR light got this right by accident: its helicity sum is COHERENT, so the (+,-) and (-,+) terms
        # cancel for odd K on their own. The UNPOLARIZED sum is incoherent -- |+-> and |-+> are accumulated as
        # separate states, and neither is individually symmetric -- so odd K survived spuriously there.
        # MEASURED: for H 1s -> 2s (J = 1/2 -> 1/2, K = {0,1}) K = 1 supplied 92 % of the unpolarized cross
        # section, while the K = 0 part was already correct (0.107e-27 against 0.214e-27 for linear -- exactly
        # the factor 2 that Mg gives, and Mg is K = 0 only).
        #
        # NOTE this restriction belongs to the MONOCHROMATIC single-beam scheme ONLY. With two distinguishable
        # beams the photons are not identical and every K contributes -- which is one more reason the
        # bichromatic case is worth completing.
        if  isodd( Int(Basics.twice(K)/2) )    continue    end
        qList = AngularMomentum.m_values(K)
        for  q in qList
            for  lambda1  in [-1, 1]
                for  lambda2  in [-1, 1]
                    amp = ComplexF64(0.)
                    for  mp1 in settings.multipoles
                        for  mp2 in settings.multipoles
                            if   mp1.electric   p1 = 1    else    p1 = 0    end
                            if   mp2.electric   p2 = 1    else    p2 = 0    end
                            symmetries  = AngularMomentum.allowedTotalSymmetries(symf, mp2, mp1, symi)
                            for Jsym in symmetries
                                amp = amp + (1.0im)^(mp1.L - p1 + mp2.L - p2) * (-lambda1)^p1 * (-lambda2)^p2 *
                                            sqrt( (2*mp1.L + 1)*(2*mp2.L + 1) ) * (Basics.twice(K) + 1)       *
                                            AngularMomentum.Wigner_3j(mp1.L, mp2.L, K, lambda1, lambda2, q)             *
                                            MultiPhotonTransition.getReducedAmplitudeAbsorption(K, line.finalLevel, mp2, Jsym, omega, mp1, 
                                                                                                        line.initialLevel, gauge, line.channels) 
                            end
                        end
                    end
                    tcs = tcs + abs( amp )^2
                end
            end
        end
    end
    
    # PREFACTOR CORRECTED 11-Aug-2026.  It read  8*pi^5 * alpha^2 ...  and therefore carried alpha^2 where
    # the radiation relation requires c^2 = 1/alpha^2.  Two independent arguments fix this:
    #  (i) DIMENSIONS.  Comparing this module's own two-photon emission and absorption for the SAME
    #      transition, the unconverged intermediate sum cancels and only radiation physics remains,
    #          dA/domega_1 = (g_l/g_u) * delta * omega_1^2 omega_2^2 / (pi^4 c^4),   delta = hbar*sigma,
    #      each photon contributing the mode density c*rho(omega) = omega^2/(pi^2 c^2) that also makes the
    #      one-photon Einstein relation come out right.  With the OLD prefactor the ratio of the two
    #      routines carried NO power of c at all, alpha^2 cancelling between them, while the relation
    #      demands c^-4: a cross section and a rate cannot differ by no power of alpha.
    #  (ii) INTERNALLY.  computeTotalAlpha0 in this same file already divides by alpha^2.
    # Measured before the change: emission/absorption exceeded the relation by 1.63264e7 at Z = 1 and
    # 1.63255e7 at Z = 2 -- constant to 6e-5 over a fourfold change in omega, so a pure prefactor.  It
    # decomposes exactly as c^4/4 = 8.81614e7 times 0.185188.  The c^4/4 is installed here; the residual
    # 0.185188 (0.656 per single-photon amplitude) is NOT, being still unexplained -- it most likely
    # reflects how the reduced matrix elements carry their (2J+1) weights when amplitude(::Absorption)
    # conjugates and swaps the levels.  Absolute cross sections therefore remain low by about 5.4.
    tcs = tcs * 2*pi^5 / Defaults.getDefaults("alpha")^2 / (Basics.twice(line.initialLevel.J) + 1) / omega^2
    
    return( tcs )
end


"""
`MultiPhotonTransition.getReducedAmplitudeAbsorption(K::AngularJ64, finalLevel::Level, multipole2::EmMultipole, Jsym::LevelSymmetry, omega::Float64, 
                                                                                        multipole1::EmMultipole, initialLevel::Level, 
                                                        gauge::EmGauge, channels::Array{MultiPhotonTransition.Channel_2pAbsorptionMonochromatic,1})`  
    ... to get/return the reduced amplitude U^{K, 2gamma absorption} (K, Jf, multipole2, Jsym, omega, multipole1, Ji) from the calculated list
        of channels. An amplitude::Complex{Float64} is returned.
"""
function getReducedAmplitudeAbsorption(K::AngularJ64, finalLevel::Level, multipole2::EmMultipole, Jsym::LevelSymmetry, omega::Float64, 
                                                                            multipole1::EmMultipole, initialLevel::Level, 
                                        gauge::EmGauge, channels::Array{MultiPhotonTransition.Channel_2pAbsorptionMonochromatic,1})
    U = Complex(0.);    found = false
    for channel in channels
        if  K == channel.K  &&  omega == channel.omega  &&   Jsym == channel.Jsym       &&  multipole1 == channel.multipole1  &&  
            multipole2 == channel.multipole2            &&   (gauge == channel.gauge  ||  EmGauge("Magnetic")  == channel.gauge)
            U = channel.amplitude;    found = true
        end
    end 
    
    if    found                                
            # println("U^{K, 2gamma absorption} (..) = $U   ** amplitude found. ")
    else  println("U^{$K, 2gamma absorption} (Jf, mp2=$multipole2, Jsym=$Jsym, omega=$omega, mp1=$multipole1, Ji)  ** NO amplitude found.")
    end 
    
    return( U )
end


"""
`MultiPhotonTransition.determineChannels_2pAbsorptionMonochromatic(omega::Float64, finalLevel::Level, initialLevel::Level, 
                                                                        settings::MultiPhotonTransition.Settings)`  
    ... to determine a list of MultiPhotonTransition.Channel_2pAbsorptionMonochromatic for a transitions from the initial to 
        final level and by taking into account the particular settings of for this computation; 
        an Array{MultiPhotonTransition.Channel_2pAbsorptionMonochromatic,1} is returned.
"""
function determineChannels_2pAbsorptionMonochromatic(omega::Float64, finalLevel::Level, initialLevel::Level, settings::MultiPhotonTransition.Settings)
    channels   = MultiPhotonTransition.Channel_2pAbsorptionMonochromatic[];   
    symi       = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
    for  mp1 in settings.multipoles
        for  mp2 in settings.multipoles
            symmetries  = AngularMomentum.allowedTotalSymmetries(symf, mp2, mp1, symi)
            Klist       = oplus(symf.J, symi.J)
            for  symn in symmetries
                for  gauge in settings.gauges
                    # THE SECOND CONDITION WAS DEAD until 08-Aug-2026: it read
                    #     elseif string(mp1)[1] == 'E'  string(mp2)[1] == 'E'  &&  gauge == Basics.UseBabushkin
                    # without the first `&&`, which Julia parses as a condition on mp1 ALONE followed by a
                    # no-op expression -- so the test on mp2 and the test on the gauge were both discarded.
                    # Harmless for E1E1, which is all this module has ever been run with, and wrong as soon as
                    # mixed multipoles are requested: an (E1, M1) pair with gauge = UseCoulomb fell through to
                    # here and was pushed as a BABUSHKIN channel. Fixed; E1E1 results are bit-identical.
                    # Include further restrictions if appropriate
                    if     string(mp1)[1] == 'E' && string(mp2)[1] == 'E'  &&   gauge == Basics.UseCoulomb
                        for K in Klist  push!(channels, MultiPhotonTransition.Channel_2pAbsorptionMonochromatic(K, omega, mp1, mp2, Basics.Coulomb, symn, 0.) )     end
                    elseif string(mp1)[1] == 'E' && string(mp2)[1] == 'E'  &&   gauge == Basics.UseBabushkin
                        for K in Klist  push!(channels, MultiPhotonTransition.Channel_2pAbsorptionMonochromatic(K, omega, mp1, mp2, Basics.Babushkin, symn, 0.) )   end
                    elseif string(mp1)[1] == 'M' && string(mp2)[1] == 'M'
                        for K in Klist  push!(channels, MultiPhotonTransition.Channel_2pAbsorptionMonochromatic(K, omega, mp1, mp2, Basics.Magnetic, symn, 0.) )    end
                    end
                end 
            end
        end
    end

    return( channels )  
end


"""
`MultiPhotonTransition.determineLines_2pAbsorptionMonochromatic(finalMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                                                    settings::MultiPhotonTransition.Settings)`
    ... to determine a list of MultiPhotonTransition.Line_2pAbsorptionMonochromatic's for transitions between the levels from the given 
        initial- and final-state multiplets and by taking into account the particular selections and settings for this computation; 
        an Array{MultiPhotonTransition.Line_2pAbsorptionMonochromatic,1} is returned. Apart from the level specification, all physical 
        properties are set to zero during this initialization process.  
"""
function  determineLines_2pAbsorptionMonochromatic(finalMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                                    settings::MultiPhotonTransition.Settings)
    lines = MultiPhotonTransition.Line_2pAbsorptionMonochromatic[]
    for  iLevel  in  initialMultiplet.levels
        for  fLevel  in  finalMultiplet.levels
            if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                omega    = (fLevel.energy - iLevel.energy + settings.photonEnergyShift) / 2.
                channels = MultiPhotonTransition.determineChannels_2pAbsorptionMonochromatic(omega, fLevel, iLevel, settings) 
                push!( lines, MultiPhotonTransition.Line_2pAbsorptionMonochromatic(iLevel, fLevel, omega,
                                            EmProperty(0., 0.), EmProperty(0., 0.), EmProperty(0., 0.), EmProperty(0., 0.),
                                            EmProperty(0., 0.), channels) )
            end
        end
    end
    return( lines )
end


"""
`MultiPhotonTransition.displayLines_2pAbsorptionMonochromatic(lines::Array{MultiPhotonTransition.Line_2pAbsorptionMonochromatic,1})`  
    ... to display a list of lines and channels that have been selected due to the prior settings. A neat table of all selected 
        transitions and energies is printed but nothing is returned otherwise.
"""
function  displayLines_2pAbsorptionMonochromatic(lines::Array{MultiPhotonTransition.Line_2pAbsorptionMonochromatic,1})
    nx = 175
    println(" ")
    println("  Selected two-photon absorption lines by monochromatic and equally-polarized photons:")
    println(" ")
    println("  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=3);                         sb = sb * TableStrings.hBlank(21)
    sa = sa * TableStrings.center(12, "Energy"; na=2);              
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=2)
    sa = sa * TableStrings.center(12, "omega"; na=4);              
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.flushleft(90, "List of multipoles & intermediate level symmetries"; na=4)            
    sb = sb * TableStrings.flushleft(90, "(K-rank, multipole_1, Jsym, multipole_2, gauge), ..."; na=4)
    println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
    for  line in lines
        sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        energy = line.finalLevel.energy - line.initialLevel.energy
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", energy))      * "    "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.omega))  * "     "
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
`MultiPhotonTransition.displayTotalAlpha0_2pAbsorptionMonochromatic(stream::IO, 
                                            properties::Array{MultiPhotonTransition.AbstractMultiPhotonProperty,1},
                                            lines::Array{MultiPhotonTransition.Line_2pAbsorptionMonochromatic,1})`  
    ... to display all results, energies, rates, etc. of the selected lines. A neat table is printed but nothing is 
        returned otherwise.
"""
function  displayTotalAlpha0_2pAbsorptionMonochromatic(stream::IO, properties::Array{AbstractMultiPhotonProperty,1},
                                                        lines::Array{Line_2pAbsorptionMonochromatic,1})
    nx = 105
    wx = Defaults.convertUnits("length: from atomic to cm", 1.0)^4 / Defaults.convertUnits("energy: from atomic to Ws", 1.0)
    println(stream, " ")
    println(stream, "  Total alpha_0 parameter are given in cm^4/Ws :")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(10, "Energy"; na=4);              
    sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(10, "omega";  na=4);              
    sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(28, "Cou --  total alpha_0  -- Bab"; na=4);              
    sb = sb * TableStrings.center(28, "[cm^4/Ws]"  * "           " * "[cm^4/Ws]"; na=8)

    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    for  line in lines
        sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=4)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        energy = line.finalLevel.energy - line.initialLevel.energy
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", energy))         * "    "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.omega))     * "        "
        for property in properties
            if      typeof(property) == TotalAlpha0          
                sa = sa * @sprintf("%.4e", line.alpha0.Coulomb   * wx)        * "      "
                sa = sa * @sprintf("%.4e", line.alpha0.Babushkin * wx)        * "          "
            else
            end
        end
        println(stream, sa )
    end
    println(stream, "  ", TableStrings.hLine(nx))
    return( nothing )
end


"""
`MultiPhotonTransition.displayResults_2pAbsorptionMonochromatic(stream::IO, 
                                        properties::Array{MultiPhotonTransition.AbstractMultiPhotonProperty,1},
                                        lines::Array{MultiPhotonTransition.Line_2pAbsorptionMonochromatic,1})`  
    ... to display all results, energies, rates, etc. of the selected lines. A neat table is printed but nothing is 
        returned otherwise.
"""
function  displayResults_2pAbsorptionMonochromatic(stream::IO, properties::Array{AbstractMultiPhotonProperty,1},
                                                    lines::Array{Line_2pAbsorptionMonochromatic,1})
    nx = 75
    wx = Defaults.convertUnits("length: from atomic to cm", 1.0)^4 / Defaults.convertUnits("energy: from atomic to Ws", 1.0)
    # wx = wx * Defaults.convertUnits("time: from atomic to sec", 1.0)
    println(stream, " ")
    println(stream, "  Two-photon absorption by monochromatic and equally-polarized photons (usually from the same beam):")
    println(stream, " ")
    # THE UNIT LABEL SAID cm^4/W WHILE wx COMPUTES cm^4/(W s) -- the header and the conversion disagreed.
    # Corrected to match what is actually computed. NOTE the ABSOLUTE normalisation of these cross sections has
    # NOT been derived or verified (unlike the emission prefactor, checked against H 2s -> 1s), so the numbers
    # are provisional; the POLARIZATION RATIOS below are not, being fixed by angular algebra alone.
    println(stream, "  Cross sections [cm^4/Ws] are given for (absolute scale NOT yet verified):")
    noCs = 0  # Number of cross sections to be printed
    for property in properties
        if      typeof(property) == TotalCsLinear          
            noCs = noCs + 1;   println(stream, "    + total cross sections for linearly-polarized incident light ($noCs)")
        elseif  typeof(property) == TotalCsRightCircular   
            noCs = noCs + 1;   println(stream, "    + total cross sections for right-circularly polarized incident " *
                                                "light ($noCs); MUST vanish for J = 0 -> J = 0") 
        elseif  typeof(property) == TotalCsUnpolarized     
            noCs = noCs + 1;   println(stream, "    + total cross sections for unpolarized incident light ($noCs)")
        elseif  typeof(property) == TotalCsDensityMatrix
            noCs = noCs + 1;   println(stream, "    + total cross sections for the given Stokes parameters ($noCs); " *
                                               "J = 0 -> J = 0 only")
        end
    end
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx + 34noCs))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(10, "Energy"; na=4);              
    sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(10, "omega";  na=4);              
    sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=7)
    for no = 1:noCs
        sa = sa * TableStrings.center(28, "Cou -- cross section ($no) -- Bab"; na=4);              
        sb = sb * TableStrings.center(28, "[cm^4/Ws]" * "         " * "[cm^4/Ws]"; na=8)
    end

    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx + 34noCs)) 
    for  line in lines
        sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=4)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        energy = line.finalLevel.energy - line.initialLevel.energy
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", energy))         * "    "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.omega))     * "        "
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
            elseif  typeof(property) == TotalCsDensityMatrix
                sa = sa * @sprintf("%.4e", line.csDensityMatrix.Coulomb   * wx)   * "      "
                sa = sa * @sprintf("%.4e", line.csDensityMatrix.Babushkin * wx)   * "          "
            end
        end
        println(stream, sa )
    end
    println(stream, "  ", TableStrings.hLine(nx + 34noCs))
    # Display the TotalAlpha0 parameters if calculated
    if  TotalAlpha0()  in   properties
        MultiPhotonTransition.displayTotalAlpha0_2pAbsorptionMonochromatic(stream, properties, lines)
    end
    return( nothing )
end
