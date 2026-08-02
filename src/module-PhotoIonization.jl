
"""
`module  JAC.PhotoIonization`  
... a submodel of JAC that contains all methods for computing photoionization properties between some initial 
    and final-state multiplets.
"""
module PhotoIonization


using Printf, ..AngularMomentum, ..Basics, ..BiOrthogonal, ..Continuum, ..Defaults, ..Radial, ..Nuclear, ..ManyElectron, ..PhotoEmission,
              ..TableStrings

"""
`struct  PhotoIonization.Settings  <:  AbstractProcessSettings`  ... defines a type for the details and parameters of computing photoionization lines.

    + multipoles                    ::Array{EmMultipole}  ... Specifies the multipoles of the radiation field that are to be included.
    + gauges                        ::Array{UseGauge}     ... Specifies the gauges to be included into the computations.
    + photonEnergies                ::Array{Float64,1}    ... List of photon energies [in user-selected units].  
    + electronEnergies              ::Array{Float64,1}    ... List of electron energies; usually only one of these lists are utilized. 
    + thetas                        ::Array{Float64,1}    ... List of theta-values if angle-differential CS are calculated explicitly. 
    + calcAnisotropy                ::Bool                ... True, if the beta anisotropy parameters are to be calculated and false otherwise (o/w). 
    + calcPartialCs                 ::Bool                ... True, if partial cross sections are to be calculated and false otherwise.  
    + calcTimeDelay                 ::Bool                ... True, if time-delays are to be calculated and false otherwise.  
    + calcNonE1AngleDifferentialCS  ::Bool                ... True, if non-E1 angle-differential CS are be calculated and false otherwise.  
    + calcTensors                   ::Bool                ... True, if statistical tensors of the excited atom are to be calculated and false o/w.
    + calcBiorthogonal              ::Bool                ... True, if the (bound) initial- and final-state multiplets are first brought into
                                                            a bi-orthogonal representation (`BiOrthogonal.computeTransformation`) before the
                                                            free-electron (continuum) orbital is generated and the transition amplitudes are
                                                            evaluated, and false (the default) if the two multiplets are used as they are. This
                                                            transforms ONLY the bound spectator orbitals shared between the N-electron neutral
                                                            (initial) and (N-1)-electron ionic (final) multiplets -- i.e. exactly the kind of
                                                            core relaxation upon ionization the bi-orthogonal method is designed to correct.
                                                            The continuum orbital itself is generated afterward, per line, directly in the field
                                                            of the (possibly now bi-orthogonally-transformed) ionic core via
                                                            `Continuum.generateOrbitalForLevel`, and is NEVER itself passed through the
                                                            bi-orthogonal machinery: applying it any later (after the continuum orbital and its
                                                            placeholder counterpart on the initial side already exist) would fail, since the
                                                            initial-side placeholder for the free-electron subshell
                                                            (`Basics.generateLevelWithExtraSubshell`) is a dummy, all-zero orbital with no
                                                            physical overlap to speak of -- the per-kappa overlap matrix used by
                                                            `BiOrthogonal.computeTransformationMatrices` would be singular. Although
                                                            `BiOrthogonal.computeTransformation`'s own docstring says the two multiplets "must
                                                            have the same number of electrons", nothing in its actual implementation enforces
                                                            or requires this: `computeTransformationMatrices` only needs matching per-kappa
                                                            orbital COUNTS (not matching NoElectrons), and `generateCounterRotatingCiMatrices`
                                                            builds its counter-rotation matrix entirely from one side's OWN CSF list, with no
                                                            cross-reference to the other side's electron count at all -- the N-vs-(N-1) case
                                                            used here appears mathematically sound on inspection, but (unlike the PhotoEmission/
                                                            PhotoExcitation cases) has not yet been checked against an independent known-answer
                                                            test; treat quantitative results with appropriate caution until such a test exists.
    + printBefore                   ::Bool                ... True, if all energies and lines are printed before their evaluation.
    + lineSelection                 ::LineSelection       ... Specifies the selected levels, if any.
    + stokes                        ::ExpStokes           ... Stokes parameters of the incident radiation.
    + freeElectronShift             ::Float64             ... An overall energy shift of all free-electron energies [user-specified units].
    + lValues                       ::Array{Int64,1}      ... Orbital angular momentum of free-electrons, for which partial waves are considered.
"""
struct Settings  <:  AbstractProcessSettings 
    multipoles                      ::Array{EmMultipole}
    gauges                          ::Array{UseGauge}
    photonEnergies                  ::Array{Float64,1} 
    electronEnergies                ::Array{Float64,1} 
    thetas                          ::Array{Float64,1}
    calcAnisotropy                  ::Bool 
    calcPartialCs                   ::Bool 
    calcTimeDelay                   ::Bool 
    calcNonE1AngleDifferentialCS    ::Bool  
    calcTensors                     ::Bool
    calcBiorthogonal                ::Bool
    printBefore                     ::Bool
    lineSelection                   ::LineSelection
    stokes                          ::ExpStokes
    freeElectronShift               ::Float64 
    lValues                         ::Array{Int64,1}
end 


"""
`PhotoIonization.Settings()`  ... constructor for the default values of photoionization line computations
"""
function Settings()
    Settings(Basics.EmMultipole[E1], Basics.UseGauge[Basics.UseCoulomb, Basics.UseBabushkin], Float64[], Float64[], Float64[],
                false, false, false, false, false, false, false, LineSelection(), Basics.ExpStokes(), 0., [0,1,2,3,4,5])
end


"""
`PhotoIonization.Settings(set::PhotoIonization.Settings;`

        multipoles=..,                      gauges=..,                  photonEnergies=..,          electronEnergies=..,
        thetas=..,                          calcAnisotropy=..,          calcPartialCs..,            calcTimeDelay=..,
        calcNonE1AngleDifferentialCS=..,    calcTensors=..,             calcBiorthogonal=..,        printBefore=..,
        lineSelection=..,                   stokes=..,                  freeElectronShift=..,       lValues=.. )
                    
    ... constructor for modifying the given PhotoIonization.Settings by 'overwriting' the previously selected parameters.
"""
function Settings(set::PhotoIonization.Settings;    
    multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,                gauges::Union{Nothing,Array{UseGauge,1}}=nothing,  
    photonEnergies::Union{Nothing,Array{Float64,1}}=nothing,                electronEnergies::Union{Nothing,Array{Float64,1}}=nothing, 
    thetas::Union{Nothing,Array{Float64,1}}=nothing,                        calcAnisotropy::Union{Nothing,Bool}=nothing,
    calcPartialCs::Union{Nothing,Bool}=nothing,                             calcTimeDelay::Union{Nothing,Bool}=nothing,  
    calcNonE1AngleDifferentialCS::Union{Nothing,Bool}=nothing,              calcTensors::Union{Nothing,Bool}=nothing,
    calcBiorthogonal::Union{Nothing,Bool}=nothing,                          printBefore::Union{Nothing,Bool}=nothing,
    lineSelection::Union{Nothing,LineSelection}=nothing,
    stokes::Union{Nothing,ExpStokes}=nothing,                               freeElectronShift::Union{Nothing,Float64}=nothing,
    lValues::Union{Nothing,Array{Int64,1}}=nothing)  
    
    if  isnothing(multipoles)          multipolesx        = set.multipoles        else  multipolesx        = multipoles         end 
    if  isnothing(gauges)              gaugesx            = set.gauges            else  gaugesx            = gauges             end 
    if  isnothing(photonEnergies)      photonEnergiesx    = set.photonEnergies    else  photonEnergiesx    = photonEnergies     end 
    if  isnothing(electronEnergies)    electronEnergiesx  = set.electronEnergies  else  electronEnergiesx  = electronEnergies   end 
    if  isnothing(thetas)              thetasx            = set.thetas            else  thetasx            = thetas             end 
    if  isnothing(calcAnisotropy)      calcAnisotropyx    = set.calcAnisotropy    else  calcAnisotropyx    = calcAnisotropy     end 
    if  isnothing(calcPartialCs)       calcPartialCsx     = set.calcPartialCs     else  calcPartialCsx     = calcPartialCs      end 
    if  isnothing(calcTimeDelay)       calcTimeDelayx     = set.calcTimeDelay     else  calcTimeDelayx     = calcTimeDelay      end 
    if  isnothing(calcNonE1AngleDifferentialCS)     calcNonE1AngleDifferentialCSx = set.calcNonE1AngleDifferentialCS        else  
        calcNonE1AngleDifferentialCSx  = calcNonE1AngleDifferentialCS                                                           end 
    if  isnothing(calcTensors)         calcTensorsx       = set.calcTensors       else  calcTensorsx       = calcTensors        end
    if  isnothing(calcBiorthogonal)    calcBiorthogonalx  = set.calcBiorthogonal  else  calcBiorthogonalx  = calcBiorthogonal   end
    if  isnothing(printBefore)         printBeforex       = set.printBefore       else  printBeforex       = printBefore        end
    if  isnothing(lineSelection)       lineSelectionx     = set.lineSelection     else  lineSelectionx     = lineSelection      end 
    if  isnothing(stokes)              stokesx            = set.stokes            else  stokesx            = stokes             end 
    if  isnothing(freeElectronShift)   freeElectronShiftx = set.freeElectronShift else  freeElectronShiftx = freeElectronShift  end 
    if  isnothing(lValues)             lValuesx           = set.lValues           else  lValuesx           = lValues            end 

    Settings( multipolesx, gaugesx, photonEnergiesx, electronEnergiesx, thetasx, calcAnisotropyx, calcPartialCsx, calcTimeDelayx,
                calcNonE1AngleDifferentialCSx, calcTensorsx, calcBiorthogonalx, printBeforex, lineSelectionx, stokesx, freeElectronShiftx, lValuesx)
end


# `Base.show(io::IO, settings::PhotoIonization.Settings)`  ... prepares a proper printout of the variable settings::PhotoIonization.Settings.
function Base.show(io::IO, settings::PhotoIonization.Settings) 
    println(io, "multipoles:                    $(settings.multipoles)  ")
    println(io, "gauges:                        $(settings.gauges)  ")
    println(io, "photonEnergies:                $(settings.photonEnergies)  ")
    println(io, "electronEnergies:              $(settings.electronEnergies)  ")
    println(io, "thetas:                        $(settings.thetas)  ")
    println(io, "calcAnisotropy:                $(settings.calcAnisotropy)  ")
    println(io, "calcPartialCs:                 $(settings.calcPartialCs)  ")
    println(io, "calcTimeDelay:                 $(settings.calcTimeDelay)  ")
    println(io, "calcNonE1AngleDifferentialCS:  $(settings.calcNonE1AngleDifferentialCS)  ")
    println(io, "calcTensors:                   $(settings.calcTensors)  ")
    println(io, "calcBiorthogonal:              $(settings.calcBiorthogonal)  ")
    println(io, "printBefore:                   $(settings.printBefore)  ")
    println(io, "lineSelection:                 $(settings.lineSelection)  ")
    println(io, "stokes:                        $(settings.stokes)  ")
    println(io, "freeElectronShift:             $(settings.freeElectronShift)  ")
    println(io, "lValues:                       $(settings.lValues)  ")
end


"""
`struct  PhotoIonization.PlasmaSettings  <:  Basics.AbstractLineShiftSettings`  
    ... defines a type for the details and parameters of computing photoionization rates with plasma interactions.

    + multipoles             ::Array{Basics.EmMultipole}     ... Specifies the multipoles of the radiation field that are to be included.
    + gauges                 ::Array{Basics.UseGauge}        ... Specifies the gauges to be included into the computations.
    + photonEnergies         ::Array{Float64,1}              ... List of photon energies.  
    + printBefore            ::Bool                          ... True, if all energies and lines are printed before their evaluation.
    + lineSelection          ::LineSelection                 ... Specifies the selected levels, if any.
"""
struct PlasmaSettings  <:  Basics.AbstractLineShiftSettings 
    multipoles               ::Array{Basics.EmMultipole}
    gauges                   ::Array{Basics.UseGauge}  
    photonEnergies           ::Array{Float64,1} 
    printBefore              ::Bool 
    lineSelection            ::LineSelection
end 


"""
`PhotoIonization.PlasmaSettings()`  ... constructor for a standard instance of PhotoIonization.PlasmaSettings.
"""
function PlasmaSettings()
    PlasmaSettings([E1], [Basics.UseCoulomb], Float64[], true, LineSelection() )
end


# `Base.show(io::IO, settings::PhotoIonization.PlasmaSettings)`  ... prepares a proper printout of the settings::PhotoIonization.PlasmaSettings.
function Base.show(io::IO, settings::PhotoIonization.PlasmaSettings)
    println(io, "multipoles:              $(settings.multipoles)  ")
    println(io, "gauges:                  $(settings.gauges)  ")
    println(io, "photonEnergies:          $(settings.photonEnergies)  ")
    println(io, "printBefore:             $(settings.printBefore)  ")
    println(io, "lineSelection:           $(settings.lineSelection)  ")
end


"""
`struct  PhotoIonization.Channel`  
    ... defines a type for a photoionization channel to help characterize a single multipole and scattering (continuum) state 
        of many electron-states with a single free electron.

    + multipole      ::EmMultipole          ... Multipole of the photon absorption.
    + gauge          ::EmGauge              ... Gauge for dealing with the (coupled) radiation field.
    + kappa          ::Int64                ... partial-wave of the free electron
    + symmetry       ::LevelSymmetry        ... total angular momentum and parity of the scattering state
    + phase          ::Float64              ... phase of the partial wave
    + amplitude      ::Complex{Float64}     ... Photoionization amplitude associated with the given channel.
"""
struct  Channel
    multipole        ::EmMultipole
    gauge            ::EmGauge
    kappa            ::Int64
    symmetry         ::LevelSymmetry
    phase            ::Float64
    amplitude        ::Complex{Float64}
end


"""
`struct  Line`  ... defines a type for a photoionization line that may include the definition of channels.

    + initialLevel   ::Level                  ... initial-(state) level
    + finalLevel     ::Level                  ... final-(state) level
    + electronEnergy ::Float64                ... Energy of the (outgoing free) electron.
    + photonEnergy   ::Float64                ... Energy of the absorbed photon.
    + crossSection   ::EmProperty             ... Cross section for this photoionization.
    + angularBeta    ::EmProperty             ... beta -parameter for unpolarized targets with J=0, 1/2, 1
    + coherentDelay  ::EmProperty             ... coherent time-delay due to the selected averaging of phases.
    + incoherentDelay::EmProperty             ... incoherent time-delay due to the selected averaging of phases.
    + channels       ::Array{PhotoIonization.Channel,1}  ... List of PhotoIonization.Channels of this line.
"""
struct  Line
    initialLevel     ::Level
    finalLevel       ::Level
    electronEnergy   ::Float64
    photonEnergy     ::Float64
    crossSection     ::EmProperty
    angularBeta      ::EmProperty 
    coherentDelay    ::EmProperty  
    incoherentDelay  ::EmProperty  
    channels         ::Array{PhotoIonization.Channel,1}
end


"""
`PhotoIonization.Line(initialLevel::Level, finalLevel::Level, crossSection::Float64)`  
    ... constructor for an photoionization line between a specified initial and final level.
"""
function Line(initialLevel::Level, finalLevel::Level, crossSection::EmProperty)
    Line(initialLevel, finalLevel, totalRate, 0., 0., crossSection, EmProperty(0.), EmProperty(0.), EmProperty(0.), 
            EmProperty(0.), PhotoChannel[] )
end


# `Base.show(io::IO, line::PhotoIonization.Line)`  ... prepares a proper printout of the variable line::PhotoIonization.Line.
function Base.show(io::IO, line::PhotoIonization.Line) 
    println(io, "initialLevel:      $(line.initialLevel)  ")
    println(io, "finalLevel:        $(line.finalLevel)  ")
    println(io, "electronEnergy:    $(line.electronEnergy)  ")
    println(io, "photonEnergy:      $(line.photonEnergy)  ")
    println(io, "crossSection:      $(line.crossSection)  ")
    println(io, "angularBeta:       $(line.angularBeta)  ")
    println(io, "coherentDelay:     $(line.coherentDelay)  ")
    println(io, "incoherentDelay:   $(line.incoherentDelay)  ")
end


"""
`PhotoIonization.amplitude(kind::String, channel::PhotoIonization.Channel, omega::Float64, continuumLevel::Level, 
                                initialLevel::Level, grid::Radial.Grid)`  
    ... to compute the kind = (photoionization) amplitude  <(alpha_f J_f, epsilon kappa) J_t || O^(photoionization) || alpha_i J_i>  
        due to the electron-photon interaction for the given final and initial level, the partial wave of the outgoing 
        electron as well as the given multipole and gauge. A value::ComplexF64 is returned.
"""
function amplitude(kind::String, channel::PhotoIonization.Channel, omega::Float64, continuumLevel::Level, initialLevel::Level, grid::Radial.Grid)
    if      kind in [ "photoionization"]
    #-----------------------------------
        amp = PhotoEmission.amplitude(Absorption(), channel.multipole, channel.gauge, omega, continuumLevel, initialLevel, grid, 
                                        display=false, printout=false)
        l         = Basics.subshell_l(Subshell(101, channel.kappa))
        amplitude = (1.0im)^(-l) * exp( -im*channel.phase ) * amp
        
    else    error("stop b")
    end
    
    return( amplitude )
end


"""
`PhotoIonization.angularFunctionK(L1::Int64, L2::Int64, X::Int64, Ji::AngularJ64, Jf::AngularJ64, 
                                  kappa1::Int64, J1::AngularJ64, kappa2::Int64, J2::AngularJ64)`  
    ... to compute angular function K(...) as defined for the non-E1 angle-differential cross sections by
        Nishita Hosea (2025). No tests are made that the triangular conditions of the quantum numbers
        are fulfilled. A wa::Float64 is returned.
"""
function angularFunctionK(L1::Int64, L2::Int64, X::Int64, Ji::AngularJ64, Jf::AngularJ64, 
                          kappa1::Int64, J1::AngularJ64, kappa2::Int64, J2::AngularJ64)
    s1 = Subshell(20,kappa1);   j1 = Basics.subshell_j(s1)
    s2 = Subshell(20,kappa2);   j2 = Basics.subshell_j(s2)
    wb = (2L1+1) * (Basics.twice(j1)+1) * (Basics.twice(J1)+1) * (2L2+1) * (Basics.twice(j2)+1) * (Basics.twice(J2)+1)
    wa = (2X+1) / (Basics.twice(Ji)+1)  * sqrt(wb) * AngularMomentum.phaseFactor([Ji, -1, Jf, 1, AngularJ64(1//2)])
    wa = wa * AngularMomentum.Wigner_6j(J2, J1, X, j1, j2, Jf) * AngularMomentum.Wigner_6j(J2, J1, X, L2, L1, Ji)
    
    return( wa )
end


"""
`PhotoIonization.angularFunctionW(theta::Float64, L1::Int64, L2::Int64, X::Int64, lambda1::Int64, lambda2::Int64,
                                  kappa1::Int64, mu1::Rational{Int64}, kappa2::Int64, mu2::Rational{Int64})`  
    ... to compute angular function W(theta; ...) as defined for the non-E1 angle-differential cross sections by
        Nishita Hosea (2025). No tests are made that the triangular conditions of the quantum numbers
        are fulfilled. A  wa::Float64 is returned.
"""
function angularFunctionW(theta::Float64, L1::Int64, L2::Int64, X::Int64, lambda1::Int64, lambda2::Int64,
                          kappa1::Int64, mu1::Rational{Int64}, kappa2::Int64, mu2::Rational{Int64})
    if  abs(lambda1) > L1  ||   abs(lambda2) > L2  ||   abs(lambda2-lambda1) > X   return( 0.)    end
    s1 = Subshell(20,kappa1);   j1 = Basics.subshell_j(s1)
    s2 = Subshell(20,kappa2);   j2 = Basics.subshell_j(s2)
    wa = # AngularMomentum.Wigner_dmatrix(X, lambda2-lambda1, Int64(mu2-mu1), theta) *
         # AngularMomentum.Wigner_3j(j2, j1, X, -mu2, mu1, mu2-mu1) * 
         AngularMomentum.Wigner_3j(L2, L1, X, -lambda2, lambda1, lambda2-lambda1)
    
    return( wa )
end


"""
`PhotoIonization.computeDisplayNonE1AngleDifferentialCS(stream::IO, lines::Array{PhotoIonization.Line,1}, 
                                                        settings::PhotoIonization.Settings)`  
    ... to compute & display the non-E1 angle-differential photoionization cross sections for all PhotoIonization.Line's 
        and at all angles theta as defined in the settings. The general formula by Nishita Hosea (2025) is applied here.
        A neat table is printed for each line but nothing is returned otherwise. 
"""
function computeDisplayNonE1AngleDifferentialCS(stream::IO, lines::Array{PhotoIonization.Line,1}, settings::PhotoIonization.Settings)
    function spinDensityMatrix(lambda1::Int64, lambda2::Int64, stokes::ExpStokes)
        # Convert the Stokes parameters of the incoming light into a spin-density matrix on the indices lambda = +-1
        if      lambda1 == lambda2  == 1               return( (1.0 + stokes.P3)/2. )
        elseif  lambda1 ==  1   &&   lambda2  == -1    return( (stokes.P1 - stokes.P1*im)/2. )
        elseif  lambda1 == -1   &&   lambda2  ==  1    return( (stokes.P1 + stokes.P1*im)/2. )
        elseif  lambda1 == lambda2  == -1              return( (1.0 - stokes.P3)/2. )
        else    error("stop a")
        end
    end
    #
    angCS = Tuple{Float64, ComplexF64, ComplexF64}[]  # theta, angCs.Coulomb, angCs.Babushkin
    nx = 50
    # Define the 2x2 spin
    # Loop about all lines; a table is printed independently for each line
    for  line in lines
        # Loop over all angles theta
        for  theta in settings.thetas
            #
            # Loop twice about all channels but distinguish the two gauges
            csCoulomb = 0.;   csBabushkin = 0.
            for  cha in line.channels, chb in line.channels
                if  cha.gauge == Basics.Coulomb  &&   chb.gauge == Basics.Babushkin   continue    end  
                if  chb.gauge == Basics.Coulomb  &&   cha.gauge == Basics.Babushkin   continue    end
                s1 = Subshell(20,cha.kappa);   j1 = Basics.subshell_j(s1)
                s2 = Subshell(20,chb.kappa);   j2 = Basics.subshell_j(s2)                
                for  X = 0:20  # Test for triangular conditions for X and continue otherwise
                    if  AngularMomentum.isTriangle(cha.multipole.L, chb.multipole.L, X)             &&
                        AngularMomentum.isTriangle(cha.symmetry.J,  chb.symmetry.J, AngularJ64(X) ) &&   
                        AngularMomentum.isTriangle(j1,  j2, AngularJ64(X) )
                        K = PhotoIonization.angularFunctionK(cha.multipole.L, chb.multipole.L, X, line.initialLevel.J, line.finalLevel.J,
                                                             cha.kappa, cha.symmetry.J, chb.kappa, chb.symmetry.J)
                        # Compute the summation over lambda's and mu's
                        W = 0.
                        for  lambda1 = -1:2:1,   lambda2 = -1:2:1,   mu = -1//2:1: 1//2
                            W = W + spinDensityMatrix(lambda1, lambda2, settings.stokes) *
                                    PhotoIonization.angularFunctionW(theta, cha.multipole.L, chb.multipole.L, X, lambda1, lambda2,
                                                                     cha.kappa, mu, chb.kappa, mu) *
                                    (1.0im)^(chb.multipole.L - cha.multipole.L) * (lambda1*lambda2) / 2. *
                                     AngularMomentum.phaseMultipole(1.0im*lambda1, cha.multipole) *
                                     AngularMomentum.phaseMultipole(1.0im*lambda2, chb.multipole)    
                        end
                        #
                        if  cha.gauge != Basics.Coulomb     &&   chb.gauge != Basics.Coulomb 
                            csBabushkin =  csBabushkin + K * W * cha.amplitude * conj(chb.amplitude)      end
                        if  cha.gauge != Basics.Babushkin   &&   chb.gauge != Basics.Babushkin 
                            csCoulomb   =  csCoulomb   + K * W * cha.amplitude * conj(chb.amplitude)      end
                    end
                end 
            end
            push!(angCS, (theta, csCoulomb, csBabushkin) ) 
        end 
        #
        # Prepare and printout a table for the angle-differential cross sections
        nx = 69;   symi = LevelSymmetry(line.initialLevel.J, line.initialLevel.parity)
                   symf = LevelSymmetry(line.finalLevel.J, line.finalLevel.parity)
        println(stream, " ")
        println(stream, "  Non-E1 angle-differential cross sections for line:" *
                        "  $(line.initialLevel.index) [$symi] -- $(line.finalLevel.index) [$symf] "    )
        println(stream, " ")
        println(stream, "  + Photon energy:   $(line.photonEnergy)    [Hartree]")
        println(stream, "  + Multipoles:      $(settings.multipoles)")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(14, "theta" ; na=6);               
        sb = sb * TableStrings.center(14, "[rad]" ; na=6)
        sa = sa * TableStrings.center(44, "Coulomb -- cross sections -- Babushkin"; na=3);                        
        sb = sb * TableStrings.center(44, "    [Mb]     "; na=3)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
        for  cs in angCS
            sa = "     " * @sprintf("%.2e", cs[1])      * "      "
            sa = sa      * @sprintf("%.4e", cs[2].re)   * "  " * @sprintf("%.4e", cs[2].im)   * "    "
            sa = sa      * @sprintf("%.4e", cs[3].re)   * "  " * @sprintf("%.4e", cs[3].im)
            println(stream, sa)
        end
        println(stream, "  ", TableStrings.hLine(nx)) 
        #   
    end 
    
    return( nothing )
end


"""
`PhotoIonization.computeAmplitudesProperties(line::PhotoIonization.Line, nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64, 
                                                    settings::PhotoIonization.Settings; printout::Bool=false)`  
    ... to compute all amplitudes and properties of the given line; a line::PhotoIonization.Line is returned for which the amplitudes and 
        properties are now evaluated.
"""
function  computeAmplitudesProperties(line::PhotoIonization.Line, nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64, 
                                        settings::PhotoIonization.Settings; printout::Bool=false)
    nChannels    = PhotoIonization.Channel[];   nxChannels    = PhotoIonization.Channel[];   
    contSettings = Continuum.Settings(false, nrContinuum);      csC = csB = dtC = dtB = 0.
    for channel in line.channels
        newiLevel = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, line.initialLevel.basis.subshells)
        newfLevel = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel, newiLevel.basis.subshells)
        newiLevel = Basics.generateLevelWithExtraSubshell(Subshell(101, channel.kappa), newiLevel)
        cOrbital, phase  = Continuum.generateOrbitalForLevel(line.electronEnergy, Subshell(101, channel.kappa), newfLevel, nm, grid, contSettings)
        #
        newcLevel  = Basics.generateLevelWithExtraElectron(cOrbital, channel.symmetry, newfLevel)
        nChannel   = PhotoIonization.Channel(channel.multipole, channel.gauge, channel.kappa, channel.symmetry, phase, 0.)
        amplitude  = PhotoIonization.amplitude("photoionization", nChannel, line.photonEnergy, newcLevel, newiLevel, grid)
        push!( nChannels, PhotoIonization.Channel(nChannel.multipole, nChannel.gauge, nChannel.kappa, nChannel.symmetry, 
                                                    nChannel.phase, amplitude) )
        if       channel.gauge == Basics.Coulomb     csC = csC + abs(amplitude)^2
        elseif   channel.gauge == Basics.Babushkin   csB = csB + abs(amplitude)^2
        elseif   channel.gauge == Basics.Magnetic    csB = csB + abs(amplitude)^2;   csC = csC + abs(amplitude)^2
        end
        #
        if  settings.calcTimeDelay
            cOrbitalx, phasex  = Continuum.generateOrbitalForLevel(line.electronEnergy+0.01, Subshell(101, channel.kappa), 
                                                                    newfLevel, nm, grid, contSettings)
            newcLevelx = Basics.generateLevelWithExtraElectron(cOrbitalx, channel.symmetry, newfLevel)
            nxChannel  = PhotoIonization.Channel(channel.multipole, channel.gauge, channel.kappa, channel.symmetry, phasex, 0.)
            amplitude  = PhotoIonization.amplitude("photoionization", nxChannel, line.photonEnergy+0.01, newcLevelx, newiLevel, grid)
            push!( nxChannels, PhotoIonization.Channel(nxChannel.multipole, nxChannel.gauge, nxChannel.kappa, nxChannel.symmetry, 
                                                       nxChannel.phase, amplitude) )
        end
    end
    Ji2 = Basics.twice(line.initialLevel.J)
    csFactor     = 8 * pi^3 / Defaults.getDefaults("alpha") / line.photonEnergy
    ##  csFactor     = csFactor / 2.   # Not fully clear, arises likely from the Rydberg normalization
    ##  Correct for energy normalization 
    ##  if  line.electronEnergy < 2.0   csFactor = csFactor * (line.electronEnergy/2.0)^1.5     end
    crossSection = EmProperty(csFactor * csC, csFactor * csB)
    if    settings.calcAnisotropy
            angularBeta  = PhotoIonization.computeAngularBeta(line.initialLevel, line.finalLevel, nChannels)
    else  angularBeta  = EmProperty(0.)
    end
    if    settings.calcTimeDelay
          coherentDelay, incoherentDelay = PhotoIonization.computeTimeDelays(nChannels, nxChannels, 0.01, line.finalLevel.J)
    else  coherentDelay = EmProperty(0.);     incoherentDelay = EmProperty(0.)
    end
    #
    nLine = PhotoIonization.Line( line.initialLevel, line.finalLevel, line.electronEnergy, line.photonEnergy, 
                                    crossSection, angularBeta, coherentDelay, incoherentDelay, nChannels)

    return( nLine )
end


"""
`PhotoIonization.computeAmplitudesPropertiesPlasma(line::PhotoIonization.Line, nm::Nuclear.Model, grid::Radial.Grid,
                                                   settings::PhotoIonization.PlasmaSettings, plasmaModel::Basics.AbstractPlasmaModel)`
    ... to compute all amplitudes and properties of the given line but for the given plasma model. The photon-electron
        multipole operator itself carries no e-e Coulomb interaction to screen (unlike the Auger case), so the plasma
        dependence enters here through the continuum (photoelectron) orbital only, cf.
        Continuum.generateOrbitalForLevel(...,plasmaModel). A line::PhotoIonization.Line is returned for which the
        amplitudes and properties are now evaluated.
"""
function  computeAmplitudesPropertiesPlasma(line::PhotoIonization.Line, nm::Nuclear.Model, grid::Radial.Grid,
                                            settings::PhotoIonization.PlasmaSettings, plasmaModel::Basics.AbstractPlasmaModel)
    newChannels = PhotoIonization.Channel[];   contSettings = Continuum.Settings(false, grid.NoPoints-50);    csC = csB = 0.
    subshellList = Basics.generate(OrderedSubshellList(), line.finalLevel.basis, line.initialLevel.basis)
    Defaults.setDefaults("relativistic subshell list", subshellList; printout=false)
    #
    for channel in line.channels
        newiLevel = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, subshellList)
        newfLevel = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel, subshellList)
        newiLevel = Basics.generateLevelWithExtraSubshell(Subshell(101, channel.kappa), newiLevel)
        cOrbital, phase  = Continuum.generateOrbitalForLevel(line.electronEnergy, Subshell(101, channel.kappa), newfLevel,
                                                              nm, grid, contSettings, plasmaModel)
        newcLevel  = Basics.generateLevelWithExtraElectron(cOrbital, channel.symmetry, newfLevel)
        newChannel = PhotoIonization.Channel(channel.multipole, channel.gauge, channel.kappa, channel.symmetry, phase, 0.)
        amplitude  = PhotoIonization.amplitude("photoionization", newChannel, line.photonEnergy, newcLevel, newiLevel, grid)
        push!( newChannels, PhotoIonization.Channel(newChannel.multipole, newChannel.gauge, newChannel.kappa, newChannel.symmetry,
                                                    newChannel.phase, amplitude) )
        if       channel.gauge == Basics.Coulomb     csC = csC + abs(amplitude)^2
        elseif   channel.gauge == Basics.Babushkin   csB = csB + abs(amplitude)^2
        elseif   channel.gauge == Basics.Magnetic    csB = csB + abs(amplitude)^2;   csC = csC + abs(amplitude)^2
        end
    end
    Ji2 = Basics.twice(line.initialLevel.J)
    csFactor     = 4 * pi^2 * Defaults.getDefaults("alpha") * line.photonEnergy / (2*(Ji2 + 1))
    crossSection = EmProperty(csFactor * csC, csFactor * csB)
    newline = PhotoIonization.Line( line.initialLevel, line.finalLevel, line.electronEnergy, line.photonEnergy,
                                    crossSection, EmProperty(0.), EmProperty(0.), EmProperty(0.), newChannels)
    return( newline )
end


"""
`PhotoIonization.computeAngularBeta(iLevel::Level, fLevel::Level, channels::Array{PhotoIonization.Channel,1})`  
    ... to compute the beta anisotropy parameter for the photoionization transition i -> f with the given channels;
        here, the formula from Balashov (1994, Eq. 2.135) has been utilized. A beta::EmProperty parameter is returned.
        These (gauge-dependent) beta parameters are set to -9., if no amplitudes are calculated for the given gauge.
"""
function  computeAngularBeta(iLevel::Level, fLevel::Level, channels::Array{PhotoIonization.Channel,1})
    wnC = wnB = waC = waB = 0.
    for  ch in channels    
        if  ch.multipole != E1   continue    end      # These beta parameters are valid only in E1 approximation
        if       ch.gauge == Basics.Coulomb     wnC = wnC + conj(ch.amplitude) * ch.amplitude   
        elseif   ch.gauge == Basics.Babushkin   wnB = wnB + conj(ch.amplitude) * ch.amplitude
        else
        end
    end
    Ji = iLevel.J;    Jf = fLevel.J;
    for  ch  in channels
        if  ch.multipole != E1   continue    end  # These beta parameters are valid only in E1 approximation
        j = AngularMomentum.kappa_j(ch.kappa);    l = AngularMomentum.kappa_l(ch.kappa);   Jt = ch.symmetry.J
        for  chp  in channels 
            if  ch.gauge !=  chp.gauge    continue     end
            jp = AngularMomentum.kappa_j(chp.kappa);    lp = AngularMomentum.kappa_l(chp.kappa);   Jtp = chp.symmetry.J
            wa = AngularMomentum.phaseFactor([Jf, -1, Ji, -1, AngularJ64(1//2)]) *
                    sqrt( AngularMomentum.bracket([Jt, Jtp, j, jp, l, lp]) ) *  
                        AngularMomentum.ClebschGordan(l, AngularM64(0), lp, AngularM64(0), AngularJ64(2), AngularM64(0)) *
                        AngularMomentum.Wigner_6j(j, l, AngularJ64(1//2), lp, jp, AngularJ64(2)) * 
                        AngularMomentum.Wigner_6j(j, Jt, Jf, Jtp, jp, AngularJ64(2)) * 
                        AngularMomentum.Wigner_6j(AngularJ64(1), Jt, Ji, Jtp, AngularJ64(1), AngularJ64(2)) * 
                        ch.amplitude * conj(chp.amplitude)
            if      chp.gauge == Basics.Coulomb      waC = waC + wa
            elseif  chp.gauge == Basics.Babushkin    waB = waB + wa
            else    error("stop a")
            end
        end    
    end
    
    if  wnC == 0.   waC = ComplexF64(-9.0)    else    waC = sqrt(6.0) * waC / wnC      end
    if  wnB == 0.   waB = ComplexF64(-9.0)    else    waB = sqrt(6.0) * waB / wnB      end
    @show waC, waB
        
    return( EmProperty(waC.re, waB.re) )
end


"""
`PhotoIonization.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                    settings::PhotoIonization.Settings; output::Bool=true)`  
    ... to compute the photoIonization transition amplitudes and all properties as requested by the given settings. 
        A list of lines::Array{PhotoIonization.Lines} is returned.
"""
function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                        settings::PhotoIonization.Settings; output::Bool=true)
    if  settings.calcBiorthogonal
        initialMultiplet, finalMultiplet = BiOrthogonal.computeTransformation(initialMultiplet, finalMultiplet, grid)
    end
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")

    println("")
    printstyled("PhotoIonization.computeLines(): The computation of photo-ionization and properties starts now ... \n", color=:light_green)
    printstyled("------------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    lines = PhotoIonization.determineLines(finalMultiplet, initialMultiplet, settings)
    # Display all selected lines before the computations start
    if  settings.printBefore    PhotoIonization.displayLines(stdout, lines)    end
    # Determine maximum energy and check for consistency of the grid
    maxEnergy = 0.;   for  line in lines   maxEnergy = max(maxEnergy, line.electronEnergy)   end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
    # Calculate all amplitudes and requested properties
    newLines = PhotoIonization.Line[]
    for  line in lines
        println("\n>> Calculate photoionization amplitudes and properties for line: $(line.initialLevel.index) - $(line.finalLevel.index) " *
                "for the photon energy $(Defaults.convertUnits("energy: from atomic", line.photonEnergy)) " * Defaults.GBL_ENERGY_UNIT)
        newLine = PhotoIonization.computeAmplitudesProperties(line, nm, grid, nrContinuum, settings) 
        push!( newLines, newLine)
    end
    # Print all results to screen
    PhotoIonization.displayPhases(newLines)
    PhotoIonization.displayResults(stdout, newLines, settings)
    ## PhotoIonization.displayTimeDelay(stdout, newLines, settings)
    if  printSummary   PhotoIonization.displayResults(iostream, newLines, settings)     
                        ## PhotoIonization.displayTimeDelay(iostream, newLines, settings)    
    end
    #
    # Add printout about non-E1 angle-differential cross sections, if required
    if  settings.calcNonE1AngleDifferentialCS
        PhotoIonization.computeDisplayNonE1AngleDifferentialCS(stdout, newLines, settings)
    end
    #
    #
    if    output    return( newLines )
    else            return( nothing )
    end
end



"""
`PhotoIonization.computeLinesCascade(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                     settings::PhotoIonization.Settings, initialLevelSelection::LevelSelection; 
                                     output=true, printout::Bool=true)`  
    ... to compute the photoionization transition amplitudes and all properties as requested by the given settings. The computations
        and printout is adapted for large cascade computations by including only lines with at least one channel and by sending
        all printout to a summary file only. A list of lines::Array{PhotoIonization.Lines} is returned.
"""
function  computeLinesCascade(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                              settings::PhotoIonization.Settings, initialLevelSelection::LevelSelection; output=true, printout::Bool=true)
    
    lines = PhotoIonization.determineLines(finalMultiplet, initialMultiplet, settings)
    # Display all selected lines before the computations start
    # if  settings.printBefore    PhotoIonization.displayLines(stdout, lines)    end  
    # Determine maximum energy and check for consistency of the grid
    maxEnergy = 0.;   for  en in settings.photonEnergies     maxEnergy = max(maxEnergy, Defaults.convertUnits("energy: to atomic", en))   end
                      for  en in settings.electronEnergies   maxEnergy = max(maxEnergy, Defaults.convertUnits("energy: to atomic", en))   end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
    # Calculate all amplitudes and requested properties
    newLines = PhotoIonization.Line[]
    for  (i,line)  in  enumerate(lines)
        if  rem(i,10) == 0    println("> Photo line $i:  ... calculated ")    end
        # Do not compute line if initial level is not in initialLevelSelection()
        ## @show Basics.selectLevel(line.initialLevel, initialLevelSelection), line.initialLevel.index
        if  !Basics.selectLevel(line.initialLevel, initialLevelSelection)   continue   ## @show "jump photoioization line";    continue  
        end
        #
        newLine = PhotoIonization.computeAmplitudesProperties(line, nm, grid, nrContinuum, settings, printout=printout) 
        push!( newLines, newLine)
    end
    # Print all results to a summary file, if requested
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   PhotoIonization.displayResults(iostream, newLines, settings)   end
    #
    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`PhotoIonization.computeLinesPlasma(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                                    settings::PhotoIonization.PlasmaSettings, plasmaModel::Basics.AbstractPlasmaModel; output::Bool=true)`
    ... to compute the photoIonization transition amplitudes and all properties as requested by the given settings and
        plasma model. A list of lines::Array{PhotoIonization.Lines} is returned.
"""
function  computeLinesPlasma(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                                settings::PhotoIonization.PlasmaSettings, plasmaModel::Basics.AbstractPlasmaModel; output::Bool=true)
    println("")
    printstyled("PhotoIonization.computeLinesPlasma(): The computation of photo-ionization cross sections starts now ... \n", color=:light_green)
    printstyled("------------------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    photoSettings = PhotoIonization.Settings(settings.multipoles, settings.gauges, settings.photonEnergies, Float64[], Float64[],
                        false, false, false, false, false, settings.printBefore, settings.lineSelection, Basics.ExpStokes(), 0., [0,1,2,3,4,5])

    lines = PhotoIonization.determineLines(finalMultiplet, initialMultiplet, photoSettings)
    # Display all selected lines before the computations start
    if  settings.printBefore    PhotoIonization.displayLines(stdout, lines)    end
    # Calculate all amplitudes and requested properties
    newLines = PhotoIonization.Line[]
    for  line in lines
        newLine = PhotoIonization.computeAmplitudesPropertiesPlasma(line, nm, grid, settings, plasmaModel)
        push!( newLines, newLine)
    end
    # Print all results to screen
    PhotoIonization.displayResults(stdout, newLines, photoSettings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   PhotoIonization.displayResults(iostream, newLines, photoSettings)     end
    #
    if    output    return( newLines )
    else            return( nothing )
    end
end



"""
`PhotoIonization.computePartialCrossSectionUnpolarized(gauge::EmGauge, Mf::AngularM64, line::PhotoIonization.Line)`  
    ... to compute the partial photoionization cross section for initially unpolarized atoms by unpolarized plane-wave photons.
        A value::Float64 is returned.
"""
function  computePartialCrossSectionUnpolarized(gauge::EmGauge, Mf::AngularM64, line::PhotoIonization.Line)
    # Define an internal Racah expression for the summation over t, lambda
    function Racahexpr(kappa::Int64, Ji::AngularJ64, Jf::AngularJ64, Mf::AngularM64, J::AngularJ64, Jp::AngularJ64, 
                        L::Int64, Lp::Int64, p::Int64, pp::Int64)
        # Determine the allowed values of t
        t1 = Basics.oplus( AngularJ64(Lp), Jf);    t2 = Basics.oplus( AngularJ64(L), Jf);    tList = intersect(t1, t2)
        wb = 0.
        for  t  in tList
            for  lambda = -1:2:1
                j = AngularMomentum.kappa_j(kappa);   Mf_lambda = Basics.add(AngularM64(lambda), Mf)
                wb = wb + (1.0im * lambda)^p * (-1.0im * lambda)^pp * 
                        AngularMomentum.ClebschGordan( AngularJ64(Lp), AngularM64(lambda), Jf, Mf, t, Mf_lambda) *
                        AngularMomentum.ClebschGordan( AngularJ64(L),  AngularM64(lambda), Jf, Mf, t, Mf_lambda) *
                        AngularMomentum.Wigner_9j(j, Jp, Jf, J, Ji, AngularJ64(L), Jf, AngularJ64(Lp), t)
            end
        end
        
        return( wb )
    end
    
        wa = 0.0im;    Ji = line.initialLevel.J;    Jf = line.finalLevel.J;
    kappaList = PhotoIonization.getLineKappas(line)
    for  kappa in kappaList
        for  cha  in line.channels
            if  kappa != cha.kappa  ||  (gauge != cha.gauge  &&  gauge != Basics.Magnetic)   continue    end
            J = cha.symmetry.J;    L = cha.multipole.L;    if  cha.multipole.electric   p = 1   else    p = 0   end
            for  chp  in line.channels  
                if  kappa != chp.kappa  ||  (gauge != cha.gauge  &&  gauge != Basics.Magnetic)    continue    end
                Jp = chp.symmetry.J;    Lp = chp.multipole.L;    if  chp.multipole.electric   pp = 1   else    pp = 0   end
                wa = wa + 1.0im^(L - Lp) * (-1)^(L + Lp) * AngularMomentum.bracket([AngularJ64(L), AngularJ64(Lp), J, Jp]) *  
                        Racahexpr(kappa, Ji, Jf, Mf, J, Jp, L, Lp, p, pp) * cha.amplitude * conj(chp.amplitude)
            end
        end
    end
    csFactor = 8 * pi^3 * Defaults.getDefaults("alpha") / (2*line.photonEnergy * (Basics.twice(Ji) + 1))
    wa       = csFactor * wa

    return( wa )
end



"""
`PhotoIonization.computeStatisticalTensorUnpolarized(k::Int64, q::Int64, gauge::EmGauge, line::PhotoIonization.Line, 
                                                            settings::PhotoIonization.Settings)`  
    ... to compute the statistical tensor of the photoion in its final level after the photoionization of initially unpolarized atoms 
        by plane-wave photons with given Stokes parameters (density matrix). A value::ComplexF64 is returned.
"""
function  computeStatisticalTensorUnpolarized(k::Int64, q::Int64, gauge::EmGauge, line::PhotoIonization.Line, 
                                                settings::PhotoIonization.Settings)
    wa = 0.0im;    Ji = line.initialLevel.J;    Jf = line.finalLevel.J   
    kappaList = PhotoIonization.getLineKappas(line);    P1 = settings.stokes.P1;   P2 = settings.stokes.P2;   P3 = settings.stokes.P3
    for  kappa in kappaList
        j = AngularMomentum.kappa_j(kappa)
        for  cha  in line.channels
            if  kappa != cha.kappa  ||  (gauge != cha.gauge  &&  gauge != Basics.Magnetic)   continue    end
            J = cha.symmetry.J;    L = cha.multipole.L;    if  cha.multipole.electric   p = 1   else    p = 0   end
            for  chp  in line.channels  
                if  kappa != chp.kappa  ||  (gauge != cha.gauge  &&  gauge != Basics.Magnetic)    continue    end
                Jp = chp.symmetry.J;    Lp = chp.multipole.L;    if  chp.multipole.electric   pp = 1   else    pp = 0   end
                #
                for  lambda = -1:2:1
                    for  lambdap = -1:2:1
                        if  lambda == lambdap   wb = (1.0 + 0.0im +lambda*P3)    else    wb = P1 - lambda * P2 * im    end
                        wa = wa + wb * 1.0im^(L - Lp + p - pp) * lambda^p * lambdap^pp *
                                sqrt( AngularMomentum.bracket([AngularJ64(L), AngularJ64(Lp), J, Jp]) ) *  
                                AngularMomentum.phaseFactor([J, +1, Jp, +1, Jf, +1, Ji, +1, j, +1, AngularJ64(1)]) *
                                AngularMomentum.ClebschGordan( AngularJ64(L),  AngularM64(lambda), AngularJ64(Lp),  AngularM64(-lambda), 
                                                                AngularJ64(k),  AngularM64(q)) *
                                AngularMomentum.Wigner_6j(Jf, j, Jp, J, AngularJ64(k), Jf) * 
                                AngularMomentum.Wigner_6j(Jp, Ji, AngularJ64(Lp), AngularJ64(L), AngularJ64(k), J) * 
                                cha.amplitude * conj(chp.amplitude)
                    end
                end
            end
        end
    end
    
    wa = pi / (Basics.twice(Ji) + 1) * wa
    return( wa )
end


"""
`PhotoIonization.computeTimeDelays(channels::Array{PhotoIonization.Channel,1}, xchannels::Array{PhotoIonization.Channel,1}, 
                                    deltaE::Float64, Jf::AngularJ64)`  
    ... to compute the -- coherent and incoherent -- time delay from the channels as calculated for two neighboured photon 
        energies (deltaE = xE - E). Two tuple of two time delays (coherentDelay::EmProperty, incoherentDelay::EmProperty)
        is returned.
"""
function  computeTimeDelays(channels::Array{PhotoIonization.Channel,1}, xchannels::Array{PhotoIonization.Channel,1}, 
                            deltaE::Float64, Jf::AngularJ64)
    #== Calculate the coherent and incorent time delays separately; start with the coherent delays.
    # First attempt
    MeffC = MeffCx = MeffB = MeffBx = ComplexF64(0.)
    for  ch in channels
        if       ch.gauge == Basics.Coulomb     MeffC  = MeffC  + ch.amplitude
        elseif   ch.gauge == Basics.Babushkin   MeffB  = MeffB  + ch.amplitude
        end
    end
    for  ch in xchannels
        if       ch.gauge == Basics.Coulomb     MeffCx = MeffCx + ch.amplitude
        elseif   ch.gauge == Basics.Babushkin   MeffBx = MeffBx + ch.amplitude
        end
    end
    DeffC = angle(MeffC);    DeffCx = angle(MeffCx);    DeffB = angle(MeffB);    DeffBx = angle(MeffBx)
    coherentDelay = EmProperty( (DeffCx - DeffC) / deltaE,  - (DeffBx - DeffB) / deltaE )
    #
    # Calculate the coherent time delays
    # Second attempt
    nomTauC = denTauC = nomTauB = denTauB= ComplexF64(0.)
    for  (ic, ch) in  enumerate(channels)
        @show  "***", ch.gauge, ch.amplitude, xchannels[ic].amplitude, (xchannels[ic].amplitude - ch.amplitude) / deltaE
        if       ch.gauge == Basics.Coulomb                 
            nomTauC = nomTauC + conj(ch.amplitude) * (xchannels[ic].amplitude - ch.amplitude) / deltaE
            denTauC = denTauC + conj(ch.amplitude) * ch.amplitude
        elseif   ch.gauge == Basics.Babushkin
            nomTauB = nomTauB + conj(ch.amplitude) * (xchannels[ic].amplitude - ch.amplitude) / deltaE
            denTauB = denTauB + conj(ch.amplitude) * ch.amplitude
        end
    end
    @warn "Multiply coherentTauC by 20. for mean energy calibration !!!"
    coherentTauC = -im * nomTauC / denTauC / deltaE * 20.;   coherentTauB = im * nomTauB / denTauB
    @show coherentTauC, coherentTauB
    coherentDelay = EmProperty( coherentTauC.im, coherentTauB.im )    ==#
    #
    # Calculate the coherent time delays
    # Third attempt, explicit derivation by Nikolay, April 2024
    ## @warn "l0 = 1 ... for p_1/2, 3/2 splitting"
    @warn "l0 = 2 ... for d_3/2, 5/2 splitting"
    ampC = ampCx = ampB = ampBx = ComplexF64(0.)
    for  (ic, ch) in  enumerate(channels)
        j  = AngularMomentum.kappa_j(ch.kappa);   l  = AngularMomentum.kappa_l(ch.kappa);   l0 = AngularJ64(2);   Jc = Jf
        @show  "***", j, l, l0, Jc
        @show  "***", ch.gauge, ch.amplitude, xchannels[ic].amplitude, (xchannels[ic].amplitude - ch.amplitude) / deltaE
        ## factor = (1.0im)^( Basics.twice(l)/2 ) * sqrt(3/(4pi))  * AngularMomentum.phaseFactor([j, +1, l, +1, AngularJ64(1//2)]) *
        factor = sqrt(3/(4pi))  * AngularMomentum.phaseFactor([j, +1, l, +1, AngularJ64(1//2)]) *
                 AngularMomentum.ClebschGordan( l0, AngularM64(0), AngularJ64(1), AngularM64(0), l,  AngularM64(0)) *
                 sqrt(Basics.twice(j)+1) * AngularMomentum.Wigner_6j(Jc, AngularJ64(1//2), l0, l, AngularJ64(1), j)
        if       ch.gauge == Basics.Coulomb
            ampC  = ampC  + factor * ch.amplitude
            ampCx = ampCx + factor * xchannels[ic].amplitude
        elseif   ch.gauge == Basics.Babushkin
            ampB  = ampB  + factor * ch.amplitude
            ampBx = ampBx + factor * xchannels[ic].amplitude
        end
    end
    coherentTauC = (ampCx - ampC) / deltaE / ampC  
    coherentTauB = (ampBx - ampB) / deltaE / ampB
    ## phiEffB      = log(ampB);       phiEffBx = log(ampBx);    coherentTaulnB = (phiEffBx.re - phiEffB.re) / deltaE
    @show coherentTauC, coherentTauB
    coherentDelay = EmProperty( coherentTauC.im, coherentTauB.im )
    #
    println("\n\nChannel amplitudes M_lj for photon energy: \n")
    for channel in channels
        println("   $(Subshell(11, channel.kappa))   $(channel.gauge)     phase=$(channel.phase)    M_lj=$(channel.amplitude) ")
    end
    # Incoherent time delays
    nomDeffC = nomDeffCx = nomDeffB = nomDeffBx = 0.
    denDeffC = denDeffCx = denDeffB = denDeffBx = 0.
    for  ch in channels
        if       ch.gauge == Basics.Coulomb     nomDeffC  = nomDeffC  + abs(ch.amplitude)^2 * ch.phase
                                                denDeffC  = denDeffC  + abs(ch.amplitude)^2
        elseif   ch.gauge == Basics.Babushkin   nomDeffB  = nomDeffB  + abs(ch.amplitude)^2 * ch.phase
                                                denDeffB  = denDeffB  + abs(ch.amplitude)^2
        end
    end
    for  ch in xchannels
        if       ch.gauge == Basics.Coulomb     nomDeffCx = nomDeffCx + abs(ch.amplitude)^2 * ch.phase
                                                denDeffCx = denDeffCx + abs(ch.amplitude)^2
        elseif   ch.gauge == Basics.Babushkin   nomDeffBx = nomDeffBx + abs(ch.amplitude)^2 * ch.phase
                                                denDeffBx = denDeffBx + abs(ch.amplitude)^2
        end
    end
    DeffC  = nomDeffC  / denDeffC;     DeffB  = nomDeffB  / denDeffB
    DeffCx = nomDeffCx / denDeffCx;    DeffBx = nomDeffBx / denDeffBx
    incoherentDelay = EmProperty( (DeffCx - DeffC) / deltaE,  (DeffBx - DeffB) / deltaE )
    
    return( coherentDelay, incoherentDelay)
end


"""
`PhotoIonization.determineChannels(finalLevel::Level, initialLevel::Level, settings::PhotoIonization.Settings)`  
    ... to determine a list of photoionization Channel for a transitions from the initial to final level and by taking into account 
        the particular settings of for this computation; an Array{PhotoIonization.Channel,1} is returned.
"""
function determineChannels(finalLevel::Level, initialLevel::Level, settings::PhotoIonization.Settings)
    channels = PhotoIonization.Channel[];   
    symi = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
    if  Basics.UseCoulomb  in  settings.gauges   gaugeM = Basics.UseCoulomb    else   gaugeM = Basics.UseBabushkin    end
    for  mp in settings.multipoles
        symList = AngularMomentum.allowedMultipoleSymmetries(symi, mp)
        for  symt in symList
            kappaList = AngularMomentum.allowedKappaSymmetries(symt, symf)
            for  kappa in kappaList
                if  !(Basics.subshell_l(Subshell(10,kappa)) in  settings.lValues)    continue     end
                for  gauge in settings.gauges
                    # Include further restrictions if appropriate
                    if     string(mp)[1] == 'E'  &&   gauge == Basics.UseCoulomb      
                        push!(channels, PhotoIonization.Channel(mp, Basics.Coulomb,   kappa, symt, 0., Complex(0.)) )
                    elseif string(mp)[1] == 'E'  &&   gauge == Basics.UseBabushkin    
                        push!(channels, PhotoIonization.Channel(mp, Basics.Babushkin, kappa, symt, 0., Complex(0.)) )  
                    elseif string(mp)[1] == 'M'  &&   gauge == gaugeM                               
                        push!(channels, PhotoIonization.Channel(mp, Basics.Magnetic,  kappa, symt, 0., Complex(0.)) ) 
                    end 
                end
            end
        end
    end
    return( channels )  
end


"""
`PhotoIonization.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::PhotoIonization.Settings)`  
    ... to determine a list of PhotoIonization.Line's for transitions between levels from the initial- and final-state multiplets, 
        and  by taking into account the particular selections and settings for this computation; an Array{PhotoIonization.Line,1} 
        is returned. Apart from the level specification, all physical properties are set to zero during the initialization process.
"""
function  determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::PhotoIonization.Settings)
    lines    = PhotoIonization.Line[]
    shift_au = Defaults.convertUnits("energy: to atomic", settings.freeElectronShift)
    for  iLevel  in  initialMultiplet.levels
        for  fLevel  in  finalMultiplet.levels
            if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                # Add lines for all photon energies
                for  omega in settings.photonEnergies
                    # Photon energies are still in 'pre-defined' units; convert to Hartree
                    omega_au = Defaults.convertUnits("energy: to atomic", omega)
                    energy   = omega_au - (fLevel.energy - iLevel.energy) + shift_au
                    if  energy < 0.    continue   end  
                    channels = PhotoIonization.determineChannels(fLevel, iLevel, settings) 
                    push!( lines, PhotoIonization.Line(iLevel, fLevel, energy, omega_au, EmProperty(0.), EmProperty(0.), 
                                                        EmProperty(0.), EmProperty(0.), channels) )
                end
                # Add lines for all electron energies
                for  en in settings.electronEnergies
                    # Electron energies are still in 'pre-defined' units; convert to Hartree
                    energy_au = Defaults.convertUnits("energy: to atomic", en) + shift_au
                    omega     = energy_au + (fLevel.energy - iLevel.energy)
                    if  energy_au < 0.    continue   end  
                    channels = PhotoIonization.determineChannels(fLevel, iLevel, settings) 
                    push!( lines, PhotoIonization.Line(iLevel, fLevel, energy_au, omega, EmProperty(0.), EmProperty(0.), 
                                                        EmProperty(0.), EmProperty(0.), channels) )
                end
            end
        end
    end
    return( lines )
end


"""
`PhotoIonization.displayLineData(stream::IO, lines::Array{PhotoIonization.Line,1})`  
    ... to display the calculated data, ordered by the initial levels and the photon energies involved.
        Neat tables of all initial levels and photon energies as well as all associated cross sections are printed
        but nothing is returned otherwise.
"""
function  displayLineData(stream::IO, lines::Array{PhotoIonization.Line,1})
    # Extract and display all initial levels by their total energy and symmetry
    energies = Float64[]
    for  line  in  lines    push!(energies, line.initialLevel.energy)  end     
    energies = unique(energies);    energies = sort(energies)
    println(stream, "\n  Initial levels, available in the given photoionization line data:")
    println(stream, "\n  ", TableStrings.hLine(42))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                              sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                              sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(12, "Level energy"   ; na=3);               
    sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=3)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(42))
    for  energy in energies
        isNew = true
        for  line in lines
            if  line.initialLevel.energy == energy  &&  isNew
                sa  = "  ";    sym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity )
                sa = sa * TableStrings.center(10, TableStrings.level(line.initialLevel.index); na=2)
                sa = sa * TableStrings.center(10, string(sym); na=4)
                sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.initialLevel.energy)) * "    "
                println(stream, sa);   isNew = false
            end
        end
    end
    #
    # Extract and display photon energies for which line data are available
    omegas = Float64[]
    for  line  in  lines    push!(omegas, line.photonEnergy)   end;    
    omegas = unique(omegas);    omegas = sort(omegas)
    sa  = "    ";   
    for   omega in omegas  sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", omega)) * "   "      end
    #
    println(stream, "\n  Photon energies $(TableStrings.inUnits("energy")) for which given photoionization line data are given:")
    println(stream, "\n" * sa)
    #
    # Extract and display total cross sections, ordered by the initial levels and photon energies, for which line data are available
    println(stream, "\n  Initial level, photon energies and total cross sections for given photoionization line data: \n")
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                              sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                              sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(12, "Omega"   ; na=3);               
    sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=3)
    sa = sa * TableStrings.center(26, "Cou--Total CS--Bab"   ; na=5);               
    sb = sb * TableStrings.center(26,TableStrings.inUnits("cross section"); na=5)
    sa = sa * TableStrings.center(18, "Electron energies"   ; na=3);               
    sb = sb * TableStrings.center(18,TableStrings.inUnits("energy"); na=3)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(107))
    #
    if  length(lines) == 0    println("  >>> No photoionization data in lines here ... return ");    return(nothing)    end 
    #
    wLine = lines[1]
    for  energy in energies
        for  omega in omegas
            cs = Basics.EmProperty(0.);   electronEnergies = Float64[]
            for  line in lines
                if  line.initialLevel.energy == energy  &&   line.photonEnergy == omega     cs = cs + line.crossSection
                    push!(electronEnergies, line.electronEnergy); wLine = line   end
            end
            sa  = "  ";    sym = LevelSymmetry( wLine.initialLevel.J, wLine.initialLevel.parity )
            sa = sa * TableStrings.center(10, TableStrings.level(wLine.initialLevel.index); na=2)
            sa = sa * TableStrings.center(10, string(sym); na=4)
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", omega)) * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", cs.Coulomb))   * "  " *
                        @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", cs.Babushkin)) * "     "
            for en in electronEnergies  sa = sa * @sprintf("%.3e", Defaults.convertUnits("energy: from atomic", en)) * "  " end
            println(stream, sa)
        end
    end
    
    
    return( nothing )
end


"""
`PhotoIonization.displayLines(stream::IO, lines::Array{PhotoIonization.Line,1})`  
    ... to display a list of lines and channels that have been selected due to the prior settings. A neat table of all selected 
        transitions and energies is printed but nothing is returned otherwise.
"""
function  displayLines(stream::IO, lines::Array{PhotoIonization.Line,1})
    nx = 175
    println(stream, " ")
    println(stream, "  Selected photoionization lines:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"   ; na=0);                       sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"   ; na=2);                       sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(10, "Energy_fi"; na=3);              
    sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=3)
    sa = sa * TableStrings.center(10, "omega"; na=3);              
    sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=2)
    sa = sa * TableStrings.center(12, "Energy e_p"; na=3);              
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.flushleft(57, "List of multipoles, gauges, kappas and total symmetries"; na=4)  
    sb = sb * TableStrings.flushleft(57, "partial (multipole, gauge, total J^P)                  "; na=4)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    nchannels = 0
    for  line in lines
        sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=3)
        energy = line.finalLevel.energy - line.initialLevel.energy
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", energy))              * "   "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.photonEnergy))   * "   "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.electronEnergy)) * "    "
        kappaMultipoleSymmetryList = Tuple{Int64,EmMultipole,EmGauge,LevelSymmetry}[]
        for  i in 1:length(line.channels)
            push!( kappaMultipoleSymmetryList, (line.channels[i].kappa, line.channels[i].multipole, line.channels[i].gauge, 
                                                line.channels[i].symmetry) )
            nchannels = nchannels + 1
        end
        wa = TableStrings.kappaMultipoleSymmetryTupels(85, kappaMultipoleSymmetryList)
        sb = sa * wa[1];    println(stream,  sb )  
        for  i = 2:length(wa)
            sb = TableStrings.hBlank( length(sa) ) * wa[i];    println(stream,  sb )
        end
    end
    println(stream, "  ", TableStrings.hLine(nx), "\n")
    println(stream, "  A total of $nchannels channels need to be calculated. \n")
    #
    return( nothing )
end


"""
`PhotoIonization.displayPhases(lines::Array{PhotoIonization.Line,1})`  
    ... to display a list of lines, channels and phases of the continuum wave that have been selected due to the prior settings.
        A neat table of all selected transitions and energies is printed but nothing is returned otherwise.
"""
function  displayPhases(lines::Array{PhotoIonization.Line,1})
    nx = 185
    println(" ")
    println("  Selected photoionization lines and phases:")
    println(" ")
    println("  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"   ; na=0);                       sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"   ; na=2);                       sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(10, "Energy_fi"; na=3);              
    sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=3)
    sa = sa * TableStrings.center(10, "omega"; na=3);              
    sb = sb * TableStrings.center(10, TableStrings.inUnits("energy"); na=2)
    sa = sa * TableStrings.center(12, "Energy e_p"; na=3);              
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.flushleft(57, "List of multipoles, gauges, kappas and total symmetries"; na=4)  
    sb = sb * TableStrings.flushleft(57, "partial (multipole, gauge, total J^P, phase)           "; na=4)
    println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
    #   
    nchannels = 0
    for  line in lines
        sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=3)
        energy = line.finalLevel.energy - line.initialLevel.energy
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", energy))              * "   "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.photonEnergy))   * "   "
        sa = sa * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.electronEnergy)) * "    "
        kappaMultipoleSymmetryPhaseList = Tuple{Int64,EmMultipole,EmGauge,LevelSymmetry,Float64}[]
        for  i in 1:length(line.channels)
            push!( kappaMultipoleSymmetryPhaseList, (line.channels[i].kappa, line.channels[i].multipole, line.channels[i].gauge, 
                                                        line.channels[i].symmetry, line.channels[i].phase) )
            nchannels = nchannels + 1
        end
        wa = TableStrings.kappaMultipoleSymmetryPhaseTupels(85, kappaMultipoleSymmetryPhaseList)
        sb = sa * wa[1];    println( sb )  
        for  i = 2:length(wa)
            sb = TableStrings.hBlank( length(sa) ) * wa[i];    println( sb )
        end
    end
    println("  ", TableStrings.hLine(nx), "\n")
    println("  A total of $nchannels channels has been calculated. \n")
    #
    return( nothing )
end



"""
`PhotoIonization.displayResults(stream::IO, lines::Array{PhotoIonization.Line,1}, settings::PhotoIonization.Settings)`  
    ... to list all results, energies, cross sections, etc. of the selected lines. A neat table is printed but nothing 
        is returned otherwise.
"""
function  displayResults(stream::IO, lines::Array{PhotoIonization.Line,1}, settings::PhotoIonization.Settings)
    nx = 130
    println(stream, " ")
    println(stream, "  Total photoionization cross sections for initially unpolarized atoms by unpolarized plane-wave photons:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"   ; na=0);                       sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(18, "i--J^P--f"   ; na=2);                       sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(12, "f--Energy--i"; na=4)               
    sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(12, "omega"     ; na=4)             
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(12, "Energy e_p"; na=3)             
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=3)
    sa = sa * TableStrings.center(10, "Multipoles"; na=1);                              sb = sb * TableStrings.hBlank(13)
    sa = sa * TableStrings.center(30, "Cou -- Cross section -- Bab"; na=3)      
    sb = sb * TableStrings.center(30, TableStrings.inUnits("cross section") * "          " * 
                                            TableStrings.inUnits("cross section"); na=3)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    # 
    wx = 1.0 ## wx = 1.0  (Schippers, August'23; wx = 2.0; wx = pi/2)
    for  line in lines
        sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=3)
        en = line.finalLevel.energy - line.initialLevel.energy
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", en))                  * "    "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.photonEnergy))   * "    "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.electronEnergy)) * "    "
        multipoles = EmMultipole[]
        for  ch in line.channels
            multipoles = push!( multipoles, ch.multipole)
        end
        multipoles = unique(multipoles);   mpString = TableStrings.multipoleList(multipoles) * "          "
        sa = sa * TableStrings.flushleft(11, mpString[1:10];  na=2)
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", wx * line.crossSection.Coulomb))     * "    "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", wx * line.crossSection.Babushkin))   * "                 "
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    #==
    # Print, if useful, the sum of all cross sections (Schippers, August'23)
    tcs = Basics.EmProperty(0.);    for  line in lines   tcs = tcs + line.crossSection   end 
    sa = repeat(" ", 102)
    sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", wx * tcs.Coulomb))     * "    "
    sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", wx * tcs.Babushkin))   * "                 "
    println(stream, sa)  ==#
    #
    #
    if  settings.calcAnisotropy
        nx = 120
        println(stream, " ")
        println(stream, "  Angular beta-parameters in E1 approximation for unpolarized target atoms with Ji = 0, 1/2, 1:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"   ; na=0);                       sb = sb * TableStrings.hBlank(18)
        sa = sa * TableStrings.center(18, "i--J^P--f"   ; na=2);                       sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.center(12, "f--Energy--i"; na=4)               
        sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(12, "omega"     ; na=4)             
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(12, "Energy e_p"; na=3)             
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=3)
        sa = sa * TableStrings.center(30, "Cou -- angular beta_2 -- Bab"; na=3);       sb = sb * TableStrings.hBlank(33)  
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
        #   
        for  line in lines
            sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                            fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=3)
            en = line.finalLevel.energy - line.initialLevel.energy
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", en))                  * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.photonEnergy))   * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.electronEnergy)) * "   "
            sa = sa * @sprintf("% .6e", line.angularBeta.Coulomb)     * "   "
            sa = sa * @sprintf("% .6e", line.angularBeta.Babushkin)   * "   "
            println(stream, sa)
        end
        println(stream, "  ", TableStrings.hLine(nx))
    end
    #
    #
    if  settings.calcTimeDelay
        nx = 150
        println(stream, " ")
        println(stream, "  (Averaged) Time-delays of individual lines:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"   ; na=0);                       sb = sb * TableStrings.hBlank(18)
        sa = sa * TableStrings.center(18, "i--J^P--f"   ; na=2);                       sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.center(12, "f--Energy--i"; na=4)               
        sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(12, "omega"     ; na=4)             
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(12, "Energy e_p"; na=3)             
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=3)
        sa = sa * TableStrings.center(30, "Cou -- coherent delay -- Bab"; na=3)      
        sb = sb * TableStrings.center(30, TableStrings.inUnits("time"); na=4)  
        sa = sa * TableStrings.center(30, "Cou -- incoherent delay -- Bab"; na=3);    
        sb = sb * TableStrings.center(30, TableStrings.inUnits("time"); na=4)  
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
        #   
        for  line in lines
            sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                            fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=3)
            en = line.finalLevel.energy - line.initialLevel.energy
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", en))                  * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.photonEnergy))   * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.electronEnergy)) * "   "
            sa = sa * @sprintf("% .6e", Defaults.convertUnits("time: from atomic", line.coherentDelay.Coulomb))     * "   "
            sa = sa * @sprintf("% .6e", Defaults.convertUnits("time: from atomic", line.coherentDelay.Babushkin))   * "     "
            sa = sa * @sprintf("% .6e", Defaults.convertUnits("time: from atomic", line.incoherentDelay.Coulomb))   * "   "
            sa = sa * @sprintf("% .6e", Defaults.convertUnits("time: from atomic", line.incoherentDelay.Babushkin)) * "     "
            println(stream, sa)
        end
        println(stream, "  ", TableStrings.hLine(nx))
    end
    #
    #
    if  settings.calcPartialCs  
        nx = 144 
        println(stream, " ")
        println(stream, "  Partial cross sections for initially unpolarized atoms by unpolarized plane-wave photons:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
            sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"   ; na=0);                       sb = sb * TableStrings.hBlank(18)
        sa = sa * TableStrings.center(18, "i--J^P--f"   ; na=2);                       sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.center(12, "f--Energy--i"; na=4)               
        sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(12, "omega"     ; na=4)             
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(12, "Energy e_p"; na=3)             
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=3)
        sa = sa * TableStrings.center(10, "Multipoles"; na=1);                              sb = sb * TableStrings.hBlank(13)
        sa = sa * TableStrings.center( 7, "M_f"; na=1);                                     sb = sb * TableStrings.hBlank(11)
        sa = sa * TableStrings.center(30, "Cou -- Partial cross section -- Bab"; na=3)      
        sb = sb * TableStrings.center(30, TableStrings.inUnits("cross section") * "          " * 
                                            TableStrings.inUnits("cross section"); na=3)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
        #   
        for  line in lines
            sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                            fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=3)
            en = line.finalLevel.energy - line.initialLevel.energy
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", en))                  * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.photonEnergy))   * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.electronEnergy)) * "    "
            multipoles = EmMultipole[]
            for  ch in line.channels
                multipoles = push!( multipoles, ch.multipole)
            end
            multipoles = unique(multipoles);   mpString = TableStrings.multipoleList(multipoles) * "          "
            sa = sa * TableStrings.flushleft(11, mpString[1:10];  na=2)
            println(stream, sa)
            MfList = Basics.projections(line.finalLevel.J)
            for  Mf in MfList
                sb  = TableStrings.hBlank(97)
                wac = PhotoIonization.computePartialCrossSectionUnpolarized(Basics.Coulomb, Mf, line)
                wab = PhotoIonization.computePartialCrossSectionUnpolarized(Basics.Babushkin, Mf, line)
                sb  = sb * TableStrings.flushright( 8, string(Mf))                             * "       "
                sb  = sb * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", wac.re))     * "    "
                sb  = sb * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", wab.re))     * "    "
                println(stream, sb)
            end
        end
        println(stream, "  ", TableStrings.hLine(nx))
    end
    #
    #
    if  settings.calcTensors  
        nx = 144 
        println(stream, " ")
        println(stream, "  Reduced statistical tensors of the photoion in its final level after the photoionization ")
        println(stream, "  of initially unpolarized atoms by plane-wave photons with given Stokes parameters (density matrix):")
        println(stream, "\n     + tensors are printed for k = 0, 1, 2 and if non-zero only.")
        println(stream,   "     + Stokes parameters are:  P1 = $(settings.stokes.P1),  P2 = $(settings.stokes.P2),  P3 = $(settings.stokes.P3) ")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
            sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"   ; na=0);                       sb = sb * TableStrings.hBlank(18)
        sa = sa * TableStrings.center(18, "i--J^P--f"   ; na=2);                       sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.center(12, "f--Energy--i"; na=4)               
        sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(12, "omega"     ; na=4)             
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(12, "Energy e_p"; na=3)             
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=3)
        sa = sa * TableStrings.center(10, "Multipoles"; na=1);                              sb = sb * TableStrings.hBlank(13)
        sa = sa * TableStrings.center(10, "k    q"; na=4);                                  sb = sb * TableStrings.hBlank(11)
        sa = sa * TableStrings.center(30, "Cou --  rho_kq (J_f)  -- Bab"; na=3)      
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
        #   
        for  line in lines
            sa  = "";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                            fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=3)
            en = line.finalLevel.energy - line.initialLevel.energy
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", en))                  * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.photonEnergy))   * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.electronEnergy)) * "    "
            multipoles = EmMultipole[]
            for  ch in line.channels
                multipoles = push!( multipoles, ch.multipole)
            end
            multipoles = unique(multipoles);   mpString = TableStrings.multipoleList(multipoles) * "          "
            sa = sa * TableStrings.flushleft(11, mpString[1:10];  na=2)
            println(stream, sa)
            for  k = 0:2
                for q = -k:k
                    sb   = TableStrings.hBlank(102)
                    rhoc = PhotoIonization.computeStatisticalTensorUnpolarized(k, q, Basics.Coulomb,   line, settings)
                    rhob = PhotoIonization.computeStatisticalTensorUnpolarized(k, q, Basics.Babushkin, line, settings)
                    sb   = sb * string(k) * " " * TableStrings.flushright( 4, string(q))             * "       "
                    sb   = sb * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", rhoc.re))     * "    "
                    sb   = sb * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", rhob.re))     * "    "
                    println(stream, sb)
                end
            end
        end
        println(stream, "  ", TableStrings.hLine(nx))
    end
    #
    return( nothing )
end

#==
"""
`PhotoIonization.displayResultsDetailed(stream::IO, line::PhotoIonization.Line, settings::PhotoIonization.Settings)`  
    ... to list the detailed results, energies, etc. for the given line. A neat table is printed but nothing 
        is returned otherwise.
"""
function  displayResultsDetailed(stream::IO, line::PhotoIonization.Line, settings::PhotoIonization.Settings)
    symi = LevelSymmetry(line.initialLevel.J, line.initialLevel.parity)
    symf = LevelSymmetry(line.finalLevel.J, line.finalLevel.parity)
    aeffC = ComplexF64(0.);     aeffB = ComplexF64(0.);    
    beffC = Float64(0.);        beffB = Float64(0.);       deffC = Float64(0.);     deffB = Float64(0.)
    for  channel in line.channels   
        if      channel.gauge == Basics.Coulomb     aeffC = aeffC + channel.amplitude
                                                    beffC = beffC + abs(channel.amplitude)^2
                                                    deffC = deffC + abs(channel.amplitude)^2 * atan(channel.amplitude.im, channel.amplitude.re)
                                                    ## deffC = deffC + abs(channel.amplitude)^2 * angle(channel.amplitude)
        elseif  channel.gauge == Basics.Babushkin   aeffB = aeffB + channel.amplitude  
                                                    beffB = beffB + abs(channel.amplitude)^2
                                                    deffB = deffB + abs(channel.amplitude)^2 * atan(channel.amplitude.im, channel.amplitude.re)
                                                    ## deffB = deffB + abs(channel.amplitude)^2 * angle(channel.amplitude)
        else
        end
    end
    deffC = deffC / beffC;    deffB = deffB / beffB  
    #
    sa   = "\n  Results for PI line from the transition  $(line.initialLevel.index) -  $(line.finalLevel.index):  " *
            "  $symi  - $symf "
    
    println(stream, sa, "\n  ", TableStrings.hLine(length(sa)-3))
    println(stream, "\n  Photon energy               = " * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.photonEnergy)) *
                        "  " * TableStrings.inUnits("energy") )
    println(stream,   "  Electron energy             = " * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.electronEnergy)) *
                        "  " * TableStrings.inUnits("energy") )
    println(stream,   "  Total cross section         = " * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", line.crossSection.Coulomb)) *
                        "  " * TableStrings.inUnits("cross section") * " (Coulomb gauge)" )
    println(stream,   "                                " * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", line.crossSection.Babushkin)) *
                        "  " * TableStrings.inUnits("cross section") * " (Babushkin gauge)" )
    println(stream,   "  A^eff (Re, Im, abs, phi)    = " * @sprintf("%.4e", aeffC.re)   * "   " * @sprintf("%.4e", aeffC.im) * "   " * 
                                                            @sprintf("%.4e", abs(aeffC)) * "   " * @sprintf("%.4e", angle(aeffC)) * "  (Coulomb gauge)" )
    println(stream,   "                              = " * @sprintf("%.4e", aeffB.re)   * "   " * @sprintf("%.4e", aeffB.im) * "   " * 
                                                            @sprintf("%.4e", abs(aeffB)) * "   " * @sprintf("%.4e", angle(aeffB)) * "  (Babushkin gauge)" )
    println(stream,   "  A^r_eff, del_eff (non-coh)  = " * @sprintf("%.4e", sqrt(beffC))* "   " * @sprintf("%.4e", deffC)        * "  (Coulomb gauge)" )
    println(stream,   "                              = " * @sprintf("%.4e", sqrt(beffB))* "   " * @sprintf("%.4e", deffB)        * "  (Babushkin gauge)" )
    println(stream,   "  A^r_eff, del_eff (coherent) = " * @sprintf("%.4e", abs(aeffC)) * "   " * @sprintf("%.4e", atan(aeffC.im, aeffC.re)) * "  (Coulomb gauge)" )
    println(stream,   "                              = " * @sprintf("%.4e", abs(aeffB)) * "   " * @sprintf("%.4e", atan(aeffB.im, aeffB.re)) * "  (Babushkin gauge)" )
    println(stream, "\n  Kappa    Total J^P  Mp     Gauge               Amplitude          Real-Amplitude  Cross section (b)   Phase    atan()" *
                    "\n  " * TableStrings.hLine(112) )
                
    for  channel in line.channels 
        Ji2 = Basics.twice(line.initialLevel.J)
        csFactor = 4 * pi^3 / Defaults.getDefaults("alpha") / line.photonEnergy
        # csFactor     = 4 * pi^2 * Defaults.getDefaults("alpha") * line.photonEnergy / (2*(Ji2 + 1))
        cs = csFactor * abs(channel.amplitude)^2
        sg = string(channel.gauge) * "      "
        sk = "    " * string(channel.kappa)
        sb = "    " * @sprintf("%.4e", channel.amplitude.re)
        sc = "    " * @sprintf("%.4e", channel.amplitude.im)
        sp = "    " * @sprintf("%.4e", channel.phase)
        sx = "    " * @sprintf("%.4e", atan(channel.amplitude.im, channel.amplitude.re) )
        sa = " " * sk[end-3:end] * "          " * string(channel.symmetry) * "    " * string(channel.multipole) * "     " * sg[1:11] * 
                sb[end-12:end] * sc[end-12:end] * "    " * @sprintf("%.4e", abs(channel.amplitude)) *
                "       " * @sprintf("%.4e", Defaults.convertUnits("cross section: from atomic", cs)) * "   " * sp[end-12:end] * 
                "       " * sx[end-12:end]
        println(stream, sa) 
    end
    
    return( nothing )
end   ==#



"""
`PhotoIonization.extractPhotonEnergies(lines::Array{PhotoIonization.Line,1})`  
    ... to extract all photon energies for which photoionization data and cross sections are provided by lines;
        an list of energies::Array{Float64,1} is returned.
"""
function  extractPhotonEnergies(lines::Array{PhotoIonization.Line,1})
    pEnergies = Float64[]
    for  line  in  lines    push!(pEnergies, line.photonEnergy)  end     
    pEnergies = unique(pEnergies)
    pEnergies = sort(pEnergies)
    
    return( pEnergies )
end


"""
`PhotoIonization.extractLines(lines::Array{PhotoIonization.Line,1}, omega::Float64)`  
    ... to extract from lines all those that refer to the given omega;
        a reduced list rLines::Array{PhotoIonization.Line,1} is returned.
"""
function  extractLines(lines::Array{PhotoIonization.Line,1}, omega::Float64)
    rLines = PhotoIonization.Line[]
    for  line  in  lines    
        if  line.photonEnergy == omega   @show  "aa", omega;    push!(rLines, line)  end  
    end
    
    return( rLines )
end


"""
`PhotoIonization.extractCrossSection(lines::Array{PhotoIonization.Line,1}, omega::Float64, initialLevel)`  
    ... to extract from lines the total PI cross section that refer to the given omega and initial level;
        a cross section cs::EmProperty is returned.
"""
function  extractCrossSection(lines::Array{PhotoIonization.Line,1}, omega::Float64, initialLevel)
    cs = Basics.EmProperty(0.)
    for  line  in  lines    
        if  line.initialLevel.index == initialLevel.index  &&  line.initialLevel.energy == initialLevel.energy  &&  
            line.photonEnergy       == omega 
            cs = cs + line.crossSection
        end  
    end
    
    return( cs )
end


"""
`PhotoIonization.extractCrossSection(lines::Array{PhotoIonization.Line,1}, omega::Float64, shell::Shell, initialLevel)`  
    ... to extract from lines the total PI cross section that refer to the given omega and initial level and to 
        to the ionization of an electron from shell; a cross section cs::EmProperty is returned.
"""
function  extractCrossSection(lines::Array{PhotoIonization.Line,1}, omega::Float64, shell::Shell, initialLevel)
    cs = Basics.EmProperty(0.)
    for  line  in  lines    
        if  line.initialLevel.index == initialLevel.index  &&  line.initialLevel.energy == initialLevel.energy  &&  
            line.photonEnergy       == omega 
            # Now determined of whether the photoionization refers to the given shell
            confi     = Basics.extractConfiguration(Basics.LeadingConfiguration(), line.initialLevel)
            conff     = Basics.extractConfiguration(Basics.LeadingConfiguration(), line.finalLevel)
            shellOccs = Basics.extractFromConfigurations(Basics.OccupationDifference(), confi, conff)
            if  length(shellOccs) > 1  ||  shellOccs[1][2] < 0      error("stop a")   end
            if  shellOccs[1][1] == shell   cs = cs + line.crossSection      end  
        end  
    end
    
    return( cs )
end


"""
`PhotoIonization.interpolateCrossSection(lines::Array{PhotoIonization.Line,1}, omega::Float64, initialLevel)`  
    ... to interpolate (or extrapolate) from lines the total PI cross section for any given omega and 
        initial level. The procedure applies a linear interpolation/extrapolation by just using the 
        cross sections from the two nearest (given) omega points; a cross section cs::EmProperty is returned.
"""
function  interpolateCrossSection(lines::Array{PhotoIonization.Line,1}, omega::Float64, initialLevel)
    # First determine for which omegas cross sections are available and which associated cross sections
    # are to be applied in the interpolation
    omegas     = PhotoIonization.extractPhotonEnergies(lines)
    oms, diffs = Basics.determineNearestPoints(omega, 2, omegas)
    if  oms[1] < oms[2]     om1 = oms[1];   om2 = oms[2];   diff = - diffs[1]  
    else                    om1 = oms[2];   om2 = oms[1];   diff = - diffs[2] 
    end
    cs1        = PhotoIonization.extractCrossSection(lines, om1, initialLevel)
    cs2        = PhotoIonization.extractCrossSection(lines, om2, initialLevel)
    if  omega < om1  &&  omega < om2   ||  omega > om1  &&  omega > om2
        @warn("No extrapolation of cross sections; cs = 0.");       
        return( Basics.EmProperty(0.) )
    end
    # Interpolate/extrapolate linearly with the cross section data for the two omegas
    wm         = 1 / (om2 - om1) * (cs2 - cs1)
    cs         = cs1  +  diff * wm
    ## @show om1, om2, diff, cs1.Coulomb, cs2.Coulomb, cs.Coulomb
    
    return( cs )
end


"""
`PhotoIonization.getLineKappas(line::PhotoIonization.Line)`  
    ... returns a list of kappa-values (partial waves) which contribute to the given line, to which one or several channels are 
        assigned. An kappaList::Array{Int64,1} is returned.
"""
function getLineKappas(line::PhotoIonization.Line)
    kappaList = Int64[]
    for  ch  in line.channels
        kappaList = union(kappaList, ch.kappa)
    end
    
    return( kappaList )
end

end # module
