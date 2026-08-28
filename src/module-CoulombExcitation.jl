


#######################################################################################################################
#######################################################################################################################

"""
`module  JAC.CoulombExcitation`  
... a submodel of JAC that contains all methods for computing Coulomb excitation cross sections and alignment parameters 
    for the excitation of target or projectile electrons by fast ion impact.
"""
module CoulombExcitation


using  Printf, ..AngularMomentum, ..Basics, ..Defaults, ..ManyElectron, ..Nuclear, ..Radial, ..RadialIntegrals, ..SpinAngular, ..TableStrings

"""
`struct  CoulombExcitation.Settings  <:  AbstractProcessSettings`  
    ... defines a type for the details and parameters of computing Coulomb-excitation amplitudes as well as 
        partial and total cross of selected lines. All cross sections are always calculated in momentum space and by
        assuming the excitation by a proton with given ion energy. For other targets, these cross sections need to be 
        multipled with (Z_target)^2.

    + ionEnergies             ::Array{Float64,1}    ... List of ion energies [MeV/u].
    + calcAlignment           ::Bool                ... True, if alignment parameters to be calculated and false otherwise.
    + printBefore             ::Bool                ... True, if all energies and lines are printed before their evaluation.
    + lineSelection           ::LineSelection       ... Specifies the selected levels, if any.
    + zerosGL                 ::Int64               ... Number of Gauss-Legendre zeros in the integration over the momentum transfer.
"""
struct Settings  <:  AbstractProcessSettings 
    ionEnergies               ::Array{Float64,1} 
    calcAlignment             ::Bool 
    printBefore               ::Bool 
    lineSelection             ::LineSelection
    zerosGL                   ::Int64
end 


"""
`CoulombExcitation.Settings()`  ... constructor for the default values of Coulomb-excitation line computations.
"""
function Settings()
    Settings(Float64[], false, false, LineSelection(), 0)
end


"""
`CoulombExcitation.Settings(set::CoulombExcitation.Settings;`

        ionEnergies=..,         calcAlignment=..,           printBefore=..,         lineSelection=..,
        zerosGL=..)
                    
    ... constructor for modifying the given CoulombExcitation.Settings by 'overwriting' the previously selected parameters.
"""
function Settings(set::CoulombExcitation.Settings;    
    ionEnergies::Union{Nothing,Array{Float64,1}}=nothing,           calcAlignment::Union{Nothing,Bool}=nothing,    
    printBefore::Union{Nothing,Bool}=nothing,                       lineSelection::Union{Nothing,LineSelection}=nothing,
    zerosGL::Union{Nothing,Int64}=nothing)  
    
    ## `ionEnergiesx`, not `ionEnergies` (fixed 08-Aug-2026): the if-branch assigned the PARAMETER instead of the
    ## local that the constructor call below uses, so every copy-construction that did not pass ionEnergies
    ## explicitly died in an UndefVarError -- i.e. the normal case. Found by the new structural test in TestFrames.
    if  isnothing(ionEnergies)           ionEnergiesx         = set.ionEnergies             else  ionEnergiesx         = ionEnergies       end
    if  isnothing(calcAlignment)         calcAlignmentx       = set.calcAlignment           else  calcAlignmentx       = calcAlignment     end 
    if  isnothing(printBefore)           printBeforex         = set.printBefore             else  printBeforex         = printBefore       end 
    if  isnothing(lineSelection)         lineSelectionx       = set.lineSelection           else  lineSelectionx       = lineSelection     end 
    if  isnothing(zerosGL)               zerosGLx             = set.zerosGL                 else  zerosGLx             = zerosGL           end 
    
    Settings( ionEnergiesx, calcAlignmentx, printBeforex, lineSelectionx, zerosGLx)
end


# `Base.show(io::IO, settings::CoulombExcitation.Settings)`  ... prepares a proper printout of the variable settings::CoulombExcitation.Settings.
function Base.show(io::IO, settings::CoulombExcitation.Settings) 
    println(io, "ionEnergies:              $(settings.ionEnergies)  ")
    println(io, "calcAlignment:            $(settings.calcAlignment)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "lineSelection:            $(settings.lineSelection)  ")
    println(io, "zerosGL:                  $(settings.zerosGL)  ")
end


"""
`struct  CoulombExcitation.Channel`  
    ... defines a type for a Coulomb-excitation channel to characterize a single q- and w-dependent (complex) amplitude.
        Each of these amplitudes are built on a number of (L, M)-dependent matrix elements which are all computed "on fly".
        
    + q              ::Float64              ... q-value of the magnetic subline.
    + w              ::Float64              ... Gauss-Legendre weight of this amplitude in the integration.
    + amplitude      ::Complex{Float64}     
        ... Coulomb-excitation amplitude K^(Coulex) (vec{q}; alpha_i J_i M_i --> alpha_f J_f M_f) associated with the given channel.
"""
struct  Channel
    q                ::Float64  
    w                ::Float64 
    amplitude        ::Complex{Float64} 
end


# `Base.show(io::IO, channel::CoulombExcitation.Channel)`  ... prepares a proper printout of the variable channel::CoulombExcitation.Channel.
function Base.show(io::IO, channel::CoulombExcitation.Channel) 
    println(io, "q:             $(channel.q)  ")
    println(io, "w:             $(channel.w)  ")
    println(io, "amplitude:     $(channel.amplitude)  ")
end


"""
`struct  CoulombExcitation.MagneticLine`  
    ... defines a type for a Coulomb-excitation magnetic lines to characterize a pair of Mi and Mf values, along with an
        (sub-) cross section and a list of channels
        Each of these amplitudes are built on a number of (L, M)-dependent matrix elements which are all computed "on fly".
        
    + Mi             ::AngularM64           ... magnetic quantum number of initial level.
    + Mf             ::AngularM64           ... magnetic quantum number of final level.
    + partialCs      ::Float64              ... partial cross section (alpha_i J_i M_i --> alpha_f J_f M_f) of this magnetic line
    + channels       ::Array{CoulombExcitation.Channel,1}  ... channels of the magnetic line.
"""
struct  MagneticLine
    Mi               ::AngularM64
    Mf               ::AngularM64
    partialCs        ::Float64
    channels         ::Array{CoulombExcitation.Channel,1}
end


# `Base.show(io::IO, channel::CoulombExcitation.MagneticLine)`  ... prepares a proper printout of the variable channel::CoulombExcitation.MagneticLine.
function Base.show(io::IO, mLine::CoulombExcitation.MagneticLine) 
    println(io, "Mi:            $(mLine.Mi)  ")
    println(io, "Mf:            $(mLine.Mf)  ")
    println(io, "partialCs:     $(mLine.partialCs)  ")
    println(io, "channels:      $(mLine.channels)  ")
end


"""
`struct  CoulombExcitation.Line`  ... defines a type for a Coulomb-excitation line that may include the definition of channels.

    + initialLevel   ::Level                  ... initial-(state) level
    + finalLevel     ::Level                  ... final-(state) level
    + ionEnergy      ::Float64                ... ion energy [MeV/u].
    + q0             ::Float64                ... minimum momentum transfer q0 that, for a given transition, is equivalent to ionEnergy.
    + totalCs        ::Float64                ... total cross section (alpha_i J_i --> alpha_f J_f) for this line.
    + alignmentA2    ::Float64                ... Alignment A_2 parameter.
    + alignmentA4    ::Float64                ... Alignment A_4 parameter.
    + mLines         ::Array{MagneticLine,1}  ... List of CoulombExcitation.MagneticLine's of this line.
"""
struct  Line
    initialLevel     ::Level
    finalLevel       ::Level
    ionEnergy        ::Float64
    q0               ::Float64
    totalCs          ::Float64
    alignmentA2      ::Float64
    alignmentA4      ::Float64
    mLines           ::Array{CoulombExcitation.MagneticLine,1}
end


# `Base.show(io::IO, line::CoulombExcitation.Line)`  ... prepares a proper printout of the variable line::CoulombExcitation.Line.
function Base.show(io::IO, line::CoulombExcitation.Line) 
    println(io, "initialLevel:      $(line.initialLevel)  ")
    println(io, "finalLevel:        $(line.finalLevel)  ")
    println(io, "ionEnergy:         $(line.ionEnergy)  ")
    println(io, "q0:                $(line.q0)  ")
    println(io, "totalCs:           $(line.totalCs)  ")
    println(io, "alignmentA2:       $(line.alignmentA2)  ")
    println(io, "alignmentA4:       $(line.alignmentA4)  ")
    println(io, "mLines:            $(line.mLines)  ")
end

#######################################################################################################################
#######################################################################################################################


"""
`CoulombExcitation.betaProjectile(ionEnergy)`  
    ... returns the (projectile) beta = v/c for ions of given ionEnergy [MeV/u]; a beta::Float64 is returned.
"""
function  betaProjectile(ionEnergy)
    gamma  = 1.0 + ionEnergy / 938.272;    beta = sqrt(1.0 - 1.0/gamma^2)
    return( beta )
end


"""
`CoulombExcitation.computeKjYme(finalLevel::Level, L::AngularJ64, initialLevel::Level, q::Float64, grid::Radial.Grid)`
    ... computes the (many-electron) reduced matrix element <alpha_f J_f || Sum_i j_L(q r_i) Y_L(n_i) || alpha_i J_i>,
        following Eq. (8) of Surzhykov, Jentschura, Stohlker, Gumberidze, Fritzsche, Phys. Rev. A 77, 042722 (2008);
        a me::ComplexF64 is returned.
"""
function  computeKjYme(finalLevel::Level, L::AngularJ64, initialLevel::Level, q::Float64, grid::Radial.Grid)
    if  initialLevel.basis.subshells == finalLevel.basis.subshells
        iLevel = initialLevel;   fLevel = finalLevel
    else
        subshells = Basics.merge(initialLevel.basis.subshells, finalLevel.basis.subshells)
        iLevel    = Level(initialLevel, subshells)
        fLevel    = Level(finalLevel, subshells)
    end

    Lint = Int64( AngularMomentum.oneJ(L) )
    nf   = length(fLevel.basis.csfs);    ni = length(iLevel.basis.csfs)
    matrix = zeros(ComplexF64, nf, ni)
    #
    for  r = 1:nf
        if  fLevel.basis.csfs[r].J != fLevel.J  ||  fLevel.basis.csfs[r].parity != fLevel.parity    continue    end
        for  s = 1:ni
            if  iLevel.basis.csfs[s].J != iLevel.J  ||  iLevel.basis.csfs[s].parity != iLevel.parity    continue    end
            subshellList = fLevel.basis.subshells
            opa = SpinAngular.OneParticleOperator(Lint, Basics.multipoleParity(EmMultipole(Lint, true)))
            wa  = SpinAngular.computeCoefficients(opa, fLevel.basis.csfs[r], iLevel.basis.csfs[s], subshellList)
            me  = 0.
            for  coeff in wa
                orbf = fLevel.basis.orbitals[coeff.a];   orbi = iLevel.basis.orbitals[coeff.b]
                ja   = Basics.subshell_2j(orbf.subshell)
                kernel = sqrt( (2Lint+1) / (4pi) ) * AngularMomentum.CL_reduced_me(orbf.subshell, Lint, orbi.subshell) *
                         RadialIntegrals.GrantJL(Lint, q, orbf, orbi, grid) / sqrt(ja + 1)
                me = me + coeff.T * kernel * sqrt( Basics.twice(fLevel.J) + 1 )
            end
            matrix[r,s] = me
        end
    end
    amplitude = transpose(fLevel.mc) * matrix * iLevel.mc
    return( amplitude )
end


"""
`CoulombExcitation.computeKjTme(finalLevel::Level, t::AngularJ64, L::AngularJ64, initialLevel::Level, q::Float64, grid::Radial.Grid)`
    ... computes the (many-electron) reduced matrix element <alpha_f J_f || Sum_i j_L(q r_i) alpha_i . T_tL(n_i) || alpha_i J_i>,
        following Eq. (8) of Surzhykov, Jentschura, Stohlker, Gumberidze, Fritzsche, Phys. Rev. A 77, 042722 (2008);
        a me::ComplexF64 is returned.
"""
function  computeKjTme(finalLevel::Level, t::AngularJ64, L::AngularJ64, initialLevel::Level, q::Float64, grid::Radial.Grid)
    if  initialLevel.basis.subshells == finalLevel.basis.subshells
        iLevel = initialLevel;   fLevel = finalLevel
    else
        subshells = Basics.merge(initialLevel.basis.subshells, finalLevel.basis.subshells)
        iLevel    = Level(initialLevel, subshells)
        fLevel    = Level(finalLevel, subshells)
    end

    Lint = Int64( AngularMomentum.oneJ(L) );    tint = Int64( AngularMomentum.oneJ(t) )
    nf   = length(fLevel.basis.csfs);    ni = length(iLevel.basis.csfs)
    matrix = zeros(ComplexF64, nf, ni)
    #
    for  r = 1:nf
        if  fLevel.basis.csfs[r].J != fLevel.J  ||  fLevel.basis.csfs[r].parity != fLevel.parity    continue    end
        for  s = 1:ni
            if  iLevel.basis.csfs[s].J != iLevel.J  ||  iLevel.basis.csfs[s].parity != iLevel.parity    continue    end
            subshellList = fLevel.basis.subshells
            opa = SpinAngular.OneParticleOperator(tint, Basics.multipoleParity(EmMultipole(Lint, false)))
            wa  = SpinAngular.computeCoefficients(opa, fLevel.basis.csfs[r], iLevel.basis.csfs[s], subshellList)
            me  = 0.
            for  coeff in wa
                orbf = fLevel.basis.orbitals[coeff.a];   orbi = iLevel.basis.orbitals[coeff.b]
                kapa = orbf.subshell.kappa;   kapb = orbi.subshell.kappa
                ja   = Basics.subshell_2j(orbf.subshell)
                kernel = ( AngularMomentum.sigma_TtL_reduced_me(-kapa, Lint, tint,  kapb) * RadialIntegrals.GrantIL0(Lint, q, orbi, orbf, grid)  -
                           AngularMomentum.sigma_TtL_reduced_me( kapa, Lint, tint, -kapb) * RadialIntegrals.GrantIL0(Lint, q, orbf, orbi, grid) ) / sqrt(ja + 1)
                me = me + coeff.T * kernel * sqrt( Basics.twice(fLevel.J) + 1 )
            end
            matrix[r,s] = me
        end
    end
    amplitude = transpose(fLevel.mc) * matrix * iLevel.mc
    return( amplitude )
end


"""
`CoulombExcitation.computeAmplitude(channel::CoulombExcitation.Channel, Mi::AngularM64, Mf::AngularM64, 
                                    line::CoulombExcitation.Line, grid::Radial.Grid, 
                                    settings::CoulombExcitation.Settings; printout::Bool=true)`  
    ... to compute the amplitudes K^(Coulex) (vec{q}; alpha_i J_i M_i --> alpha_f J_f M_f);
        an amplitude::ComplexF64 is returned.
"""
function  computeAmplitude(channel::CoulombExcitation.Channel, Mi::AngularM64, Mf::AngularM64, 
                                    line::CoulombExcitation.Line, grid::Radial.Grid, 
                                    settings::CoulombExcitation.Settings; printout::Bool=true)
    Ji  = line.initialLevel.J;        Jf  = line.finalLevel.J;    amplitude = ComplexF64(0.)
    Jix = AngularMomentum.oneJ(Ji);   Jfx = AngularMomentum.oneJ(Jf)
    Mix = AngularMomentum.oneM(Mi);   Mfx = AngularMomentum.oneM(Mf)
    
    beta = CoulombExcitation.betaProjectile(line.ionEnergy)
    
    # Simply sum over all summations; determine the range of t and L, when it becomes relevant
    ts = AngularMomentum.j_values(Ji, Jf)
    for  t  in  ts
        tx = AngularMomentum.oneJ(t)
        wa = AngularMomentum.ClebschGordan(Jix, Mix, tx, Mfx-Mix, Jfx, Mfx) / sqrt(2*Jfx + 1.0)
        Ls = AngularMomentum.j_values(t, AngularJ64(1))
        Mval = Mfx - Mix
        wb = ComplexF64(0.)
        for  L  in Ls
            Lx = AngularMomentum.oneJ(L);   Lint = Int64(Lx);   Mint = Int64(Mval)
            if  abs(Mint) > Lint   continue   end
            wc = ComplexF64(0.)
            if  t == L   wc = wc + CoulombExcitation.computeKjYme(line.finalLevel, L, line.initialLevel, channel.q, grid)   end
            ## RATIP's coulex_pure_matrix() multiplies every magnetic (F23-based) contribution by an explicit
            ## cmplx(zero,-one) = -i prefactor (rabs_coulex.f90, all three of terms ii/iii/iv) that Eq. (8) of
            ## Surzhykov et al. does not show explicitly -- this was missing here and was the confirmed root cause
            ## of the sigma(Mi,Mf) != sigma(-Mi,-Mf) asymmetry found and diagnosed earlier: with this factor, the
            ## L=t bracket becomes A + i*beta*CG(Mf)*B, and since CG(-Mf)=-CG(Mf), bracket(-Mf)=conj(bracket(Mf)),
            ## so |bracket(-Mf)|=|bracket(Mf)| exactly. Confirmed both algebraically and empirically (exact
            ## Mf<->-Mf symmetry now holds without any post-hoc symmetrization).
            wc = wc + im * beta * AngularMomentum.ClebschGordan(Lx, Mval, 1., 0., tx, Mval) *
                      CoulombExcitation.computeKjTme(line.finalLevel, t, L, line.initialLevel, channel.q, grid)
            wc = im^Lx * conj( AngularMomentum.sphericalYlm(Lint, Mint, acos(line.q0/channel.q), 0.) ) * wc
            wb = wb + wc
        end
        amplitude = amplitude + wa * wb
    end
    
    return( amplitude )
end


"""
`CoulombExcitation.computeAmplitudesProperties(mLine::CoulombExcitation.MagneticLine, line::CoulombExcitation.Line,
                                               grid::Radial.Grid, settings::CoulombExcitation.Settings; printout::Bool=true)`  
    ... to compute all amplitudes and properties of the given magnetic line by using various parameters from line;
        a mline::CoulombExcitation.MagneticLine is returned for which the amplitudes and properties have now been evaluated.
"""
function  computeAmplitudesProperties(mLine::CoulombExcitation.MagneticLine, line::CoulombExcitation.Line,
                                      grid::Radial.Grid, settings::CoulombExcitation.Settings; printout::Bool=true)
    # Compute the amplitudes K^(Coulex) for the given magnetic line
    newChannels = CoulombExcitation.Channel[]
    for  channel in mLine.channels
        amplitude = CoulombExcitation.computeAmplitude(channel, mLine.Mi, mLine.Mf, line, grid, settings; printout=printout)
        push!(newChannels, CoulombExcitation.Channel(channel.q, channel.w, amplitude))
    end

    # Compute the partial cross section; collect parameters for line
    Ji2 = Basics.twice(line.initialLevel.J);    beta = CoulombExcitation.betaProjectile(line.ionEnergy);    wSum = 0.;
    for ch in newChannels
        wa = (ch.q^2 - (line.q0^2 * beta^2) )
        wSum = wSum + ch.q * ch.w / (wa^2) * real( ch.amplitude * conj(ch.amplitude) )
    end
    partialCs = wSum * 2pi * (8pi * Defaults.getDefaults("alpha") / beta)^2 / (Ji2 + 1)

    newmLine = CoulombExcitation.MagneticLine(mLine.Mi, mLine.Mf, partialCs, newChannels)
    return( newmLine )
end


"""
`CoulombExcitation.computeAmplitudesProperties(line::CoulombExcitation.Line, grid::Radial.Grid, 
                                               settings::CoulombExcitation.Settings; printout::Bool=true)`  
    ... to compute the amplitudes, cross sections and alignment parameters of the given line; a line::CoulombExcitation.Line is 
        returned for which the amplitudes and properties have now been evaluated.
"""
function  computeAmplitudesProperties(line::CoulombExcitation.Line, grid::Radial.Grid, 
                                      settings::CoulombExcitation.Settings; printout::Bool=true)
    newmLines = CoulombExcitation.MagneticLine[]
    for  mLine in line.mLines
        newmLine = CoulombExcitation.computeAmplitudesProperties(mLine, line, grid, settings; printout=printout)
        push!(newmLines, newmLine)
    end

    # Compute the total cross section
    totalCs     = 0.;    for  mLine in newmLines    totalCs = totalCs + mLine.partialCs   end
    alignmentA2 = 0.
    alignmentA4 = 0.

    # Compute the alignment parameters A2, A4 following Eqs. (9)-(10) of Surzhykov et al., Phys. Rev. A 77, 042722 (2008);
    # this requires the cross section sigma(Jf,Mf), summed over all initial sublevels Mi for a fixed final Mf.
    if  settings.calcAlignment  &&  totalCs > 0.
        Jf  = line.finalLevel.J;    Jfx = AngularMomentum.oneJ(Jf)
        MfList  = AngularMomentum.m_values(Jf)
        sigmaMf = Dict{AngularM64, Float64}( Mf => 0. for Mf in MfList )
        for  mLine in newmLines   sigmaMf[mLine.Mf] = sigmaMf[mLine.Mf] + mLine.partialCs   end

        wa2 = 0.;   wa4 = 0.
        for  Mf in MfList
            Mfx    = AngularMomentum.oneM(Mf)
            phase  = (-1)^Int64( round(Jfx - Mfx) )
            cg2    = AngularMomentum.ClebschGordan(Jfx, Mfx, Jfx, -Mfx, 2., 0.)
            cg4    = AngularMomentum.ClebschGordan(Jfx, Mfx, Jfx, -Mfx, 4., 0.)
            wa2    = wa2 + phase * cg2 * sigmaMf[Mf]
            wa4    = wa4 + phase * cg4 * sigmaMf[Mf]
        end
        alignmentA2 = sqrt(2*Jfx + 1) * wa2 / totalCs
        alignmentA4 = sqrt(2*Jfx + 1) * wa4 / totalCs
    end

    newLine = CoulombExcitation.Line(line.initialLevel, line.finalLevel, line.ionEnergy, line.q0, totalCs, alignmentA2, alignmentA4, newmLines)
    return( newLine )
end


"""
`CoulombExcitation.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                                settings::CoulombExcitation.Settings; output=true)`  
    ... to compute the Coulomb excitation amplitudes and all properties as requested by the given settings. A list of 
        lines::Array{CoulombExcitation.Lines} is returned.
"""
function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, settings::CoulombExcitation.Settings; output=true)
    println("")
    printstyled("CoulombExcitation.computeLines(): The computation of Coulomb excitation cross sections starts now ... \n", color=:light_green)
    printstyled("----------------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
        
    sa =    "\n* Coulomb excitation cross sections for many-electron atoms and ions can be computed within different representations " *
            "and methods: " *
            "\n  they are all rather tricky and only approximate. The following assumptions and approximations are presently made: \n" *
            "\n    + All (projectile) ion energies are given in [MeV/u] which are converted into relative velocities beta = v/c. " *
            "\n    + All cross sections are computed for a (target) proton Z_t=1; multiply with Z_t^2 to find the correct cross sections. \n"
    println(sa)
    
    lines = CoulombExcitation.determineLines(finalMultiplet, initialMultiplet, settings)
    # Display all selected lines before the computations start
    if  settings.printBefore    CoulombExcitation.displayLines(stdout, lines)    end
    # Calculate all amplitudes and requested properties
    newLines = CoulombExcitation.Line[]
    for  line in lines
        newLine = CoulombExcitation.computeAmplitudesProperties(line, grid, settings) 
        push!( newLines, newLine)
    end
    # Print all results to screen
    CoulombExcitation.displayCrossSections(stdout, newLines, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    CoulombExcitation.displayCrossSections(iostream, newLines, settings)   end
    #
    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`CoulombExcitation.determineChannels(finalLevel::Level, initialLevel::Level, q0::Float64, settings::CoulombExcitation.Settings)`  
    ... to determine a list of CoulombExcitation.Channel for a magnetic line (mLine) from the initial to final level 
        and by taking into account the particular settings of for this computation; an Array{CoulombExcitation.Channel,1} is returned.
"""
function determineChannels(finalLevel::Level, initialLevel::Level, q0::Float64, settings::CoulombExcitation.Settings)
    channels = CoulombExcitation.Channel[];  
    # Compute the q's and associated weights in the interval [q0, 10*q0]
    gaussLegendre = Radial.GridGL(Radial.GridGaussLegendreFinite(), q0, 10*q0, settings.zerosGL);     qs = gaussLegendre.t;     ws = gaussLegendre.wt
    for  (iq, q)  in  enumerate(qs)
        push!(channels, CoulombExcitation.Channel(q, ws[iq], ComplexF64(0.)) )
    end

    return( channels )  
end


"""
`CoulombExcitation.determineMagneticLines(finalLevel::Level, initialLevel::Level, q0::Float64, settings::CoulombExcitation.Settings)`  
    ... to determine a list of CoulombExcitation.MagneticLine's for line from the initial to final level 
        and by taking into account the particular settings of for this computation; an Array{CoulombExcitation.Channel,1} 
        is returned.
"""
function determineMagneticLines(finalLevel::Level, initialLevel::Level, q0::Float64, settings::CoulombExcitation.Settings)
    mLines = CoulombExcitation.MagneticLine[]
    
    for  Mi in AngularMomentum.m_values(initialLevel.J)
        for  Mf in AngularMomentum.m_values(finalLevel.J)
            channels = CoulombExcitation.determineChannels(finalLevel, initialLevel, q0, settings)
            if   length(channels) == 0   continue   end
            push!( mLines, CoulombExcitation.MagneticLine(Mi, Mf, 0., channels) )
        end
    end

    return( mLines )  
end


"""
`CoulombExcitation.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::CoulombExcitation.Settings)`  
    ... to determine a list of Coulomb-excitation Line's for transitions between the levels from the given initial- and 
        final-state multiplets and by taking into account the particular selections and settings for this computation; 
        an Array{CoulombExcitation.Line,1} is returned. Apart from the level specification, all physical properties are set to 
        zero during the initialization process.  
"""
function  determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::CoulombExcitation.Settings)
    lines = CoulombExcitation.Line[]
    for  iLevel  in  initialMultiplet.levels
        for  fLevel  in  finalMultiplet.levels
            if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                if   fLevel.energy - iLevel.energy == 0.   continue   end  
                for  ionEnergy  in  settings.ionEnergies
                    # Determine q0 for the given transition and ion energy     
                    beta = CoulombExcitation.betaProjectile(ionEnergy)
                    q0   = (fLevel.energy - iLevel.energy) / (beta * Defaults.getDefaults("speed of light: c"))
                    mLines = CoulombExcitation.determineMagneticLines(fLevel, iLevel, q0, settings)
                    if   length(mLines) == 0   continue   end
                    push!( lines, CoulombExcitation.Line(iLevel, fLevel, ionEnergy, q0, 0., 0., 0., mLines) )
                end 
            end
        end
    end
    
    return( lines )
end


"""
`CoulombExcitation.displayCrossSections(stream::IO, lines::Array{CoulombExcitation.Line,1}, settings::CoulombExcitation.Settings)`  
    ... to display a list of lines, magnetic lines and channels that have been selected due to the prior settings. 
        A neat table of all selected transitions and energies is printed but nothing is returned otherwise.
"""
function  displayCrossSections(stream::IO, lines::Array{CoulombExcitation.Line,1}, settings::CoulombExcitation.Settings)
    nx = 120
    println(stream, " ")
    println(stream, "  Partial and total Coulomb-excitation lines, magnetic lines and channel parameters:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(14, "Ion energy"; na=4);              
    sb = sb * TableStrings.center(14, " [MeV/u]"; na=4)
    sa = sa * TableStrings.center(16, "Total CS"; na=2);       
    sb = sb * TableStrings.center(16, TableStrings.inUnits("cross section"); na=4)
    sa = sa * TableStrings.center(16, "Mi    Mf"; na=2);                          sb = sb * TableStrings.hBlank(16)   
    sa = sa * TableStrings.center(16, "Partial CS"; na=2);       
    sb = sb * TableStrings.center(16, TableStrings.inUnits("cross section"); na=4)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                     fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        sa = sa * @sprintf("%.8e", line.ionEnergy)  * "    "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", line.totalCs)) * "    "
        println(stream, sa )        
        for  mL in line.mLines 
            sc = "      " * string(mL.Mi) * "    " * string(mL.Mf) * "           "
            sc = sc[1:22]
            sc = sc * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", mL.partialCs))
            println(stream, TableStrings.hBlank(length(sa)) * sc )   
        end
     end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`CoulombExcitation.displayLines(stream::IO, lines::Array{CoulombExcitation.Line,1})`  
    ... to display a list of lines, magnetic lines and channels that have been selected due to the prior settings. 
        A neat table of all selected transitions and energies is printed but nothing is returned otherwise.
"""
function  displayLines(stream::IO, lines::Array{CoulombExcitation.Line,1})
    nx = 160
    println(stream, " ")
    println(stream, "  Selected Coulomb-excitation lines, magnetic lines and channel parameters:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=0);                         sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(14, "Ion energy"; na=4);              
    sb = sb * TableStrings.center(14, " [MeV/u]"; na=4)
    sa = sa * TableStrings.flushleft(50, "List of (Mi, Mf, q1, w1) for mLines + channel_1"; na=4);       
    sb = sb * TableStrings.hBlank(34)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                     fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        sa = sa * @sprintf("%.8e", line.ionEnergy)  * "    "
        MMtuples = Tuple{AngularM64, AngularM64, Float64, Float64}[]
        for  mLine  in  line.mLines
            push!( MMtuples, (mLine.Mi, mLine.Mf, mLine.channels[1].q, mLine.channels[1].w) )
        end
        wa = TableStrings.MMffTupels(100, MMtuples, "mL")
        println(stream, sa * wa[1] )        
        for  ia = 2:length(wa)   println(stream, TableStrings.hBlank(length(sa)) * wa[ia] )   end
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end

end # module
