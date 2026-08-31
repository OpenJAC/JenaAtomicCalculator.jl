
"""
`module  JAC.PhotoIonization`  
... a submodel of JAC that contains all methods for computing photoionization properties between some initial and final-state multiplets.
"""
module PhotoIonization


using Printf, ..AngularMomentum, ..Basics, ..BiOrthogonal, ..Bsplines, ..Continuum, ..Defaults, ..Radial, ..Nuclear, ..ManyElectron,
              ..PhotoEmission, ..TableStrings


"""
`struct  PhotoIonization.Settings  <:  AbstractProcessSettings`
    ... defines a type for the details and parameters of computing photoionization lines.

    + multipoles                    ::Array{EmMultipole}  ... Specifies the multipoles of the radiation field that are to be included.
    + gauges                        ::Array{UseGauge}     ... Specifies the gauges to be included into the computations.
    + photonEnergies                ::Array{Float64,1}    ... List of photon energies [in user-selected units].
    + electronEnergies              ::Array{Float64,1}    ... List of electron energies; usually only one of these lists are utilized.
    + thetas                        ::Array{Float64,1}    ... List of theta-values if angle-differential CS are calculated explicitly.
    + calcAnisotropy                ::Bool                ... True, if the beta anisotropy parameters are to be calculated and false
                                                              otherwise (o/w).
    + calcPartialCs                 ::Bool                ... True, if partial cross sections are to be calculated and false otherwise.
    + calcTimeDelay                 ::Bool                ... True, if time-delays are to be calculated and false otherwise.
    + calcNonE1AngleDifferentialCS  ::Bool                ... True, if non-E1 angle-differential CS are be calculated and false otherwise.
    + calcTensors                   ::Bool                ... True, if statistical tensors of the excited atom are to be calculated and
                                                              false o/w.
    + calcBiorthogonal              ::Bool                ... True, if the (bound) initial- and final-state multiplets are first brought
                                                              into a bi-orthogonal representation (`BiOrthogonal.computeTransformation`)
                                                              before the free-electron (continuum) orbital is generated and the transition
                                                              amplitudes are evaluated, and false (the default) if the two multiplets are
                                                              used as they are. This transforms ONLY the bound spectator orbitals shared
                                                              between the N-electron neutral (initial) and (N-1)-electron ionic (final)
                                                              multiplets -- i.e. exactly the kind of core relaxation upon ionization the
                                                              bi-orthogonal method is designed to correct. The continuum orbital itself is
                                                              generated afterward, per line, directly in the field of the (possibly now
                                                              bi-orthogonally-transformed) ionic core via
                                                              `Continuum.generateOrbitalForLevel`, and is NEVER itself passed through the
                                                              bi-orthogonal machinery: applying it any later (after the continuum orbital
                                                              and its placeholder counterpart on the initial side already exist) would fail,
                                                              since the initial-side placeholder for the free-electron subshell
                                                              (`Basics.generateLevelWithExtraSubshell`) is a dummy, all-zero orbital with no
                                                              physical overlap to speak of -- the per-kappa overlap matrix used by
                                                              `BiOrthogonal.computeTransformationMatrices` would be singular. Although
                                                              `BiOrthogonal.computeTransformation`'s own docstring says the two multiplets
                                                              "must have the same number of electrons", nothing in its actual implementation
                                                              enforces or requires this: `computeTransformationMatrices` only needs matching
                                                              per-kappa orbital COUNTS (not matching NoElectrons), and
                                                              `generateCounterRotatingCiMatrices` builds its counter-rotation matrix
                                                              entirely from one side's OWN CSF list, with no cross-reference to the other
                                                              side's electron count at all -- the N-vs-(N-1) case used here appears
                                                              mathematically sound on inspection, but (unlike the PhotoEmission/
                                                              PhotoExcitation cases) has not yet been checked against an independent
                                                              known-answer test; treat quantitative results with appropriate caution until
                                                              such a test exists.
    + printBefore                   ::Bool                ... True, if all energies and lines are printed before their evaluation.
    + lineSelection                 ::LineSelection       ... Specifies the selected levels, if any.
    + stokes                        ::ExpStokes           ... Stokes parameters of the incident radiation.
    + freeElectronShift             ::Float64             ... An overall energy shift of all free-electron energies [user-specified units].
    + lValues                       ::Array{Int64,1}      ... Orbital angular momentum of free-electrons, for which partial waves are
                                                              considered.
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


# `Base.show(io::IO, settings::PhotoIonization.Settings)`
# ... prepares a proper printout of the variable settings::PhotoIonization.Settings.
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


# `Base.show(io::IO, settings::PhotoIonization.PlasmaSettings)`
# ... prepares a proper printout of the settings::PhotoIonization.PlasmaSettings.
function Base.show(io::IO, settings::PhotoIonization.PlasmaSettings)
    println(io, "multipoles:              $(settings.multipoles)  ")
    println(io, "gauges:                  $(settings.gauges)  ")
    println(io, "photonEnergies:          $(settings.photonEnergies)  ")
    println(io, "printBefore:             $(settings.printBefore)  ")
    println(io, "lineSelection:           $(settings.lineSelection)  ")
end


"""
`struct  PhotoIonization.Channel`
    ... ONE asymptotic scattering state: the final ion plus a free electron, coupled to a total symmetry.

    + symmetry       ::LevelSymmetry                            ... total J^parity of the scattering state.
    + amplitudes     ::Array{MultipoleAmplitude,1}        ... one entry per contributing multipole.
"""
struct  Channel
    symmetry         ::LevelSymmetry
    amplitudes       ::Array{MultipoleAmplitude,1}
end


"""
`struct  PhotoIonization.PartialWave`
    ... ONE partial wave of the free electron, and the channels it serves.

    + kappa          ::Int64                            ... partial wave of the free electron.
    + energy         ::Float64                          ... energy of the free electron.
    + phase          ::Float64                          ... scattering phase; a property of (energy, kappa).
    + channels       ::Array{Channel,1}           ... the total symmetries this partial wave serves.

        The radial orbital and the phase belong HERE, not to a channel: they do not depend on the total symmetry, which is precisely why one
        kappa can serve several of them.
"""
struct  PartialWave
    kappa            ::Int64
    energy           ::Float64
    phase            ::Float64
    channels         ::Array{Channel,1}
end


"""
`struct  PhotoIonization.Line`
    ... defines a type for a photoionization line between an initial and a final level, for one photon energy, and carries the outgoing
        electron as a list of partial waves; the cross section and the angular parameters of the line are stored alongside them.

    + initialLevel   ::Level                            ... initial-(state) level
    + finalLevel     ::Level                            ... final-(state) level
    + electronEnergy ::Float64                          ... Energy of the (outgoing free) electron.
    + photonEnergy   ::Float64                          ... Energy of the absorbed photon.
    + crossSection   ::EmProperty                       ... Cross section for this photoionization.
    + angularBeta    ::EmProperty                       ... beta_2 anisotropy parameter.
    + coherentDelay  ::EmProperty                       ... coherent time delay.
    + incoherentDelay::EmProperty                       ... incoherent time delay.
    + partialWaves   ::Array{PartialWave,1}       ... partial waves, each with the channels it serves.
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
    partialWaves     ::Array{PartialWave,1}
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
    println(io, "partialWaves:      $(line.partialWaves)  ")
end


#################################################################################################################################
#################################################################################################################################

"""
`PhotoIonization.amplitude(kind::String, mp::EmMultipole, gauge::EmGauge, kappa::Int64, phase::Float64,
                                omega::Float64, continuumLevel::Level, initialLevel::Level, grid::Radial.Grid)`
    ... to compute the kind = (photoionization) amplitude  <(alpha_f J_f, epsilon kappa) J_t || O^(photoionization) || alpha_i J_i> due to
        the electron-photon interaction for the given final and initial level, the partial wave of the outgoing electron as well as the
        given multipole and gauge. A value::ComplexF64 is returned.
"""
function amplitude(kind::String, mp::EmMultipole, gauge::EmGauge, kappa::Int64, phase::Float64, omega::Float64,
                    continuumLevel::Level, initialLevel::Level, grid::Radial.Grid)
    if      kind in [ "photoionization"]
    #-----------------------------------
        amp = PhotoEmission.amplitude(Absorption(), mp, gauge, omega, continuumLevel, initialLevel, grid,
                                        display=false, printout=false)
        l         = Basics.subshell_l(Subshell(101, kappa))
        amplitude = (1.0im)^(-l) * exp( -im*phase ) * amp
        
    else    error("stop b")
    end
    
    return( amplitude )
end


"""
`PhotoIonization.angularFunctionK(L1::Int64, L2::Int64, X::Int64, Ji::AngularJ64, Jf::AngularJ64, 
                                  kappa1::Int64, J1::AngularJ64, kappa2::Int64, J2::AngularJ64)`  
    ... to compute angular function K(...) as defined for the non-E1 angle-differential cross sections by Nishita Hosea (2025). No tests are
        made that the triangular conditions of the quantum numbers are fulfilled. A wa::Float64 is returned.
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
    ... to compute angular function W(theta; ...) as defined for the non-E1 angle-differential cross sections by Nishita Hosea (2025). No
        tests are made that the triangular conditions of the quantum numbers are fulfilled. A  wa::Float64 is returned.
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
`PhotoIonization.computeAmplitudesProperties(line::PhotoIonization.Line, nm::Nuclear.Model,
        grid::Radial.Grid, nrContinuum::Int64, settings::PhotoIonization.Settings; printout::Bool=false,
        nuclearPot::Union{Nothing,Radial.Potential}=nothing,
        primitives::Union{Nothing,Bsplines.Primitives}=nothing)`
    ... computes all amplitudes and properties of the given line; a Line is returned in which the amplitudes, the phases and the cross
        section are filled.

        The physics is unchanged: the amplitudes come from the very same PhotoIonization.amplitude as the flat path. Only the ORDER of the
        work differs, and that is the point -- the continuum orbital and the initial level carrying the extra continuum subshell are formed
        ONCE PER PARTIAL WAVE, because they depend on (energy, kappa) and on nothing else. The flat path had to re-derive them per channel
        until a cache was added; here there is nothing to cache.

        An electric multipole is evaluated twice, once per gauge, into one EmPropertyC; a magnetic multipole once, into an EmPropertyC with
        equal components.
"""
function computeAmplitudesProperties(line::PhotoIonization.Line, nm::Nuclear.Model, grid::Radial.Grid,
                                           nrContinuum::Int64, settings::PhotoIonization.Settings; printout::Bool=false,
                                           nuclearPot::Union{Nothing,Radial.Potential}=nothing,
                                           primitives::Union{Nothing,Bsplines.Primitives}=nothing)
    contSettings = Continuum.Settings(false, nrContinuum)
    redILevel = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, line.initialLevel.basis.subshells)
    newfLevel = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel, redILevel.basis.subshells)
    newPartialWaves = PhotoIonization.PartialWave[]

    for  pw in line.partialWaves
        # ONE orbital and ONE initial level per partial wave -- structurally, not by a cache.
        cSubshell = Subshell(101, pw.kappa)
        newiLevel = Basics.generateLevelWithExtraSubshell(cSubshell, redILevel)
        cOrbital, phase = Continuum.generateOrbitalForLevel(line.electronEnergy, cSubshell, newfLevel, nm, grid,
                                                            contSettings; nuclearPot=nuclearPot, primitives=primitives)
        newChannels = PhotoIonization.Channel[]
        for  ch in pw.channels
            newcLevel = Basics.generateLevelWithExtraElectron(cOrbital, ch.symmetry, newfLevel)
            newAmps   = MultipoleAmplitude[]
            for  ma in ch.amplitudes
                mp = ma.multipole
                if  string(mp)[1] == 'E'
                    ampC = PhotoIonization.amplitude("photoionization", mp, Basics.Coulomb, pw.kappa, phase, line.photonEnergy, newcLevel, newiLevel, grid)
                    ampB = PhotoIonization.amplitude("photoionization", mp, Basics.Babushkin, pw.kappa, phase, line.photonEnergy, newcLevel, newiLevel, grid)
                    push!(newAmps, MultipoleAmplitude(mp, EmPropertyC(ampC, ampB)))
                else
                    # A magnetic multipole does not depend on the gauge; one evaluation, equal components.
                    ampM = PhotoIonization.amplitude("photoionization", mp, Basics.Magnetic, pw.kappa, phase, line.photonEnergy, newcLevel, newiLevel, grid)
                    push!(newAmps, MultipoleAmplitude(mp, EmPropertyC(ampM)))
                end
            end
            push!(newChannels, PhotoIonization.Channel(ch.symmetry, newAmps))
        end
        push!(newPartialWaves, PhotoIonization.PartialWave(pw.kappa, line.electronEnergy, phase, newChannels))
    end

    crossSection = PhotoIonization.computeCrossSection(newPartialWaves, line.photonEnergy)
    if    settings.calcAnisotropy
          angularBeta = PhotoIonization.computeAngularBeta(line.initialLevel, line.finalLevel, newPartialWaves)
    else  angularBeta = EmProperty(0.)
    end
    newLine = PhotoIonization.Line(line.initialLevel, line.finalLevel, line.electronEnergy, line.photonEnergy,
                                         crossSection, angularBeta, EmProperty(0.), EmProperty(0.), newPartialWaves)

    return( newLine )
end


"""
`PhotoIonization.computeAmplitudesPropertiesPlasma(line::PhotoIonization.Line, nm::Nuclear.Model, grid::Radial.Grid,
                                                   settings::PhotoIonization.PlasmaSettings, plasmaModel::Basics.AbstractPlasmaModel;
                                                   nuclearPot::Union{Nothing,Radial.Potential}=nothing,
                                                   primitives::Union{Nothing,Bsplines.Primitives}=nothing)`
    ... to compute all amplitudes and properties of the given line but for the given plasma model. The photon-electron multipole operator
        itself carries no e-e Coulomb interaction to screen (unlike the Auger case), so the plasma dependence enters here through the
        continuum (photoelectron) orbital only, cf. Continuum.generateOrbitalForLevel(...,plasmaModel). A line::PhotoIonization.Line is
        returned for which the amplitudes and properties are now evaluated.
"""
function  computeAmplitudesPropertiesPlasma(line::PhotoIonization.Line, nm::Nuclear.Model, grid::Radial.Grid,
                                            settings::PhotoIonization.PlasmaSettings, plasmaModel::Basics.AbstractPlasmaModel;
                                            nuclearPot::Union{Nothing,Radial.Potential}=nothing,
                                            primitives::Union{Nothing,Bsplines.Primitives}=nothing)
    contSettings = Continuum.Settings(false, grid.NoPoints-50);    cs = EmProperty(0., 0.)
    subshellList = Basics.generate(OrderedSubshellList(), line.finalLevel.basis, line.initialLevel.basis)
    # Display-only; set once by the driver, not here.  See Defaults.setStandardSubshellList. As above: the two symmetry-reduced levels do
    # not depend on the partial wave.
    redILevel = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, subshellList)
    newfLevel = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel,   subshellList)
    newPartialWaves = PhotoIonization.PartialWave[]

    # ONE continuum orbital per PARTIAL WAVE, as in computeAmplitudesProperties; the plasma model is the only difference between the two,
    # and it enters Continuum.generateOrbitalForLevel.
    for  pw in line.partialWaves
        cSubshell = Subshell(101, pw.kappa)
        newiLevel = Basics.generateLevelWithExtraSubshell(cSubshell, redILevel)
        cOrbital, phase = Continuum.generateOrbitalForLevel(line.electronEnergy, cSubshell, newfLevel,
                                                            nm, grid, contSettings, plasmaModel;
                                                            nuclearPot=nuclearPot, primitives=primitives)
        newChannels = PhotoIonization.Channel[]
        for  ch in pw.channels
            newcLevel = Basics.generateLevelWithExtraElectron(cOrbital, ch.symmetry, newfLevel)
            newAmps   = MultipoleAmplitude[]
            for  ma in ch.amplitudes
                mp = ma.multipole
                if  string(mp)[1] == 'E'
                    ampC = PhotoIonization.amplitude("photoionization", mp, Basics.Coulomb,   pw.kappa, phase,
                                                      line.photonEnergy, newcLevel, newiLevel, grid)
                    ampB = PhotoIonization.amplitude("photoionization", mp, Basics.Babushkin, pw.kappa, phase,
                                                      line.photonEnergy, newcLevel, newiLevel, grid)
                    amp  = EmPropertyC(ampC, ampB)
                else
                    ampM = PhotoIonization.amplitude("photoionization", mp, Basics.Magnetic,  pw.kappa, phase,
                                                      line.photonEnergy, newcLevel, newiLevel, grid)
                    amp  = EmPropertyC(ampM)
                end
                cs = cs + abs2(amp)
                push!(newAmps, MultipoleAmplitude(mp, amp))
            end
            push!(newChannels, PhotoIonization.Channel(ch.symmetry, newAmps))
        end
        push!(newPartialWaves, PhotoIonization.PartialWave(pw.kappa, line.electronEnergy, phase, newChannels))
    end
    # CORRECTED 31-Aug-2026. This line read  4 pi^2 alpha omega / (2(2J_i+1))  and was wrong by 4 pi/(alpha omega)^2
    # -- 174.7 at omega = 1 keV, and growing as omega falls. It now uses the SAME factor as
    # PhotoIonization.computeCrossSection, which is what the field-free path uses and which is validated against
    # Stobbe's closed-form 1s cross section.
    #
    # WHY THE OLD FACTOR LOOKED RIGHT: 4 pi^2 alpha omega |d|^2 IS the textbook photoionization cross section --
    # for a DIPOLE MATRIX ELEMENT. JAC's multipole amplitude is not one: it is built on j_L(alpha omega r), whose
    # leading term is (alpha omega r)^L/(2L+1)!!, so an E1 amplitude ALREADY carries one power of alpha omega and
    # |M|^2 carries two. The prefactor must therefore go as 1/(alpha omega), not as alpha omega. The old line was a
    # correct formula applied to the wrong object.
    #
    # MEASURED, both routes in one session on the same case (Ne 1s^2 2s^2 2p^6 -> 1s 2s^2 2p^6, omega = 1000 eV,
    # E1, Coulomb, NoPlasmaModel so that only the route differs):
    #        route                                  channels   sum |amp|^2      sigma/sum
    #        PhotoIonization.computeLines              2        8.84431072e-06   9.249642e+02  = 8 pi^3/(alpha w)
    #        PhotoIonization.computeLinesPlasma        2        8.84430289e-06   5.293518e+00  = 4 pi^2 alpha w/2
    # The summed amplitudes agree to six figures and the channel counts are equal, so the ENTIRE difference was the
    # prefactor. In barn: 0.2291 Mb against 0.001311 Mb, and 0.229 Mb is the right order for Ne K-shell just above
    # its 870 eV edge, which the old value was not.
    #
    # WHY IT SURVIVED SO LONG: examples/example-Jb.jl branch b is the only caller, and it runs ALL FOUR of its
    # plasma models -- including the "field-free" NoPlasmaModel() -- through THIS function. Its checks are
    # monotonicity in screening strength, absence of discontinuities, and convergence to the field-free limit;
    # every one of them is scale-invariant, so a factor common to all four is invisible to all three.
    csFactor     = 8 * pi^3 / Defaults.getDefaults("alpha") / line.photonEnergy
    crossSection = csFactor * cs
    newline = PhotoIonization.Line( line.initialLevel, line.finalLevel, line.electronEnergy, line.photonEnergy,
                                    crossSection, EmProperty(0.), EmProperty(0.), EmProperty(0.), newPartialWaves)

    return( newline )
end


"""
`PhotoIonization.computeAngularBeta(iLevel::Level, fLevel::Level, partialWaves::Array{PhotoIonization.PartialWave,1})`
    ... computes the beta_2 anisotropy parameter of the photoelectron angular distribution from the given partial waves, valid in the E1
        approximation; a beta::EmProperty is returned.

        This is the only place in the module where PAIRS of amplitudes are needed rather than one at a time, and therefore the only one
        where the two gauges could in principle be mixed. No guard against that is needed: the product `amp * conj(ampp)` of two EmPropertyC
        is componentwise, so Coulomb meets Coulomb and Babushkin meets Babushkin by construction.

        The partial wave supplies kappa -- one kappa may serve several total symmetries -- and the channel supplies the total symmetry,
        which is exactly the pairing this formula needs.
"""
function computeAngularBeta(iLevel::Level, fLevel::Level, partialWaves::Array{PhotoIonization.PartialWave,1})
    # Flatten to (kappa, symmetry, E1 amplitude) once; every pair below is then a plain double loop.
    entries = Tuple{Int64,LevelSymmetry,EmPropertyC}[]
    for  pw in partialWaves,  ch in pw.channels,  ma in ch.amplitudes
        ma.multipole == E1   ||   continue          # these beta parameters are valid only in E1 approximation
        push!(entries, (pw.kappa, ch.symmetry, ma.amplitude))
    end
    wn = EmProperty(0., 0.)
    for  (kappa, symt, amp) in entries    wn = wn + abs2(amp)    end

    Ji = iLevel.J;    Jf = fLevel.J
    waC = ComplexF64(0.);    waB = ComplexF64(0.)
    for  (kappa, symt, amp) in entries
        j = AngularMomentum.kappa_j(kappa);      l = AngularMomentum.kappa_l(kappa);      Jt  = symt.J
        for  (kappap, symtp, ampp) in entries
            jp = AngularMomentum.kappa_j(kappap);    lp = AngularMomentum.kappa_l(kappap);    Jtp = symtp.J
            wb = AngularMomentum.phaseFactor([Jf, -1, Ji, -1, AngularJ64(1//2)]) *
                    sqrt( AngularMomentum.bracket([Jt, Jtp, j, jp, l, lp]) ) *
                        AngularMomentum.ClebschGordan(l, AngularM64(0), lp, AngularM64(0), AngularJ64(2), AngularM64(0)) *
                        AngularMomentum.Wigner_6j(j, l, AngularJ64(1//2), lp, jp, AngularJ64(2)) *
                        AngularMomentum.Wigner_6j(j, Jt, Jf, Jtp, jp, AngularJ64(2)) *
                        AngularMomentum.Wigner_6j(AngularJ64(1), Jt, Ji, Jtp, AngularJ64(1), AngularJ64(2))
            # Componentwise: no gauge can leak into the other, and no guard is needed to say so.
            wc  = (amp * conj(ampp)) * wb
            waC = waC + wc.Coulomb;    waB = waB + wc.Babushkin
        end
    end

    if  wn.Coulomb   == 0.   waC = ComplexF64(-9.0)    else    waC = sqrt(6.0) * waC / wn.Coulomb      end
    if  wn.Babushkin == 0.   waB = ComplexF64(-9.0)    else    waB = sqrt(6.0) * waB / wn.Babushkin    end

    return( EmProperty(waC.re, waB.re) )
end


"""
`PhotoIonization.computeCrossSection(partialWaves::Array{PhotoIonization.PartialWave,1}, photonEnergy::Float64)`
    ... computes the total photoionization cross section from the given partial waves; an EmProperty is returned.

        THE SUMMATION IS THE SAME AS THE FLAT PATH'S, deliberately: |amplitude|^2 is summed INCOHERENTLY over multipoles as well as over
        channels. The new structure would make a coherent multipole sum easy to write, and that would be a change of physics -- a question
        to be settled separately against a reference, not as a side effect of restructuring. What the new form buys is that the rule now
        lives in ONE named place instead of being spelled out wherever channels are summed.

        No gauge appears: abs2(::EmPropertyC) returns an EmProperty, so the two gauges are carried along untouched, and a magnetic multipole
        -- whose two components are equal -- reproduces the flat path's rule of adding it to both sums without any special case.
"""
function computeCrossSection(partialWaves::Array{PhotoIonization.PartialWave,1}, photonEnergy::Float64)
    cs = EmProperty(0., 0.)
    for  pw in partialWaves,  ch in pw.channels,  ma in ch.amplitudes
        cs = cs + abs2(ma.amplitude)
    end
    # THIS IS THE CORRECT NORMALISATION, settled 28-Aug-2026 against Stobbe's closed-form 1s cross section rather
    # than by preference. Hydrogen 1s, point nucleus, E1 only, at x = omega/omega_th = 1.5, 2 and 3: this factor
    # gives 2.51913e+06, 1.02584e+06 and 3.22474e+05 barn in the Babushkin (length) gauge against Stobbe's
    # 2.09127e+06, 9.31429e+05 and 2.88397e+05 -- ratios 1.205, 1.101, 1.118, i.e. FLAT to ~10 %, which is the
    # agreement examples/example-Dd.jl already recorded. See the warning in computeAmplitudesPropertiesPlasma.
    csFactor = 8 * pi^3 / Defaults.getDefaults("alpha") / photonEnergy
    return( EmProperty(csFactor * cs.Coulomb, csFactor * cs.Babushkin) )
end


"""
`PhotoIonization.computeDisplayNonE1AngleDifferentialCS(stream::IO, lines::Array{PhotoIonization.Line,1},
                                                              settings::PhotoIonization.Settings)`
    ... computes and displays the non-E1 angle-differential photoionization cross sections for all PhotoIonization.Line's and at all angles
        theta as defined in the settings; the general formula by Nishita Hosea (2025) is applied, exactly as in
        PhotoIonization.computeDisplayNonE1AngleDifferentialCS. A table is printed for each line but nothing is returned otherwise.

        THIS IS THE ONE CONSUMER THAT PAIRS AMPLITUDES ACROSS PARTIAL WAVES. computeAngularBeta and the partial cross section pair only
        within one kappa; here kappa and kappa' are independent, and the interference between two DIFFERENT partial waves is the whole
        content of the non-dipole terms. The loop therefore runs over the flattened (kappa, symmetry, multipole, amplitude) entries of the
        line.

        Of the three flat functions that must keep the gauges apart, this is the one that does it correctly: it guards the pair from BOTH
        sides (`cha == Coulomb && chb == Babushkin` and the reverse) and then books a pair into the Coulomb sum unless either partner is
        Babushkin, and vice versa -- so a magnetic multipole, whose flat gauge label is Basics.Magnetic, enters both sums. Every one of
        those three rules is reproduced here by the single product `ma.amplitude * conj(mp.amplitude)`: it is componentwise, so the mixed
        pairs cannot form, and a magnetic amplitude has equal components, so it contributes to both. Four lines of bookkeeping become none.

        MIRRORED, NOT REPAIRED: `angCS` is filled once per line but never emptied between lines, so the table printed for the second line
        repeats the rows of the first. That is a defect of the flat function; it is reproduced here unchanged, so that the two paths can be
        compared line for line, and it is to be fixed in the flat function -- where it belongs -- rather than silently here.

        NO CALLER, AND NEVER EXECUTED -- treat that as a warning, not a reassurance. It is a user entry point for non-E1
        angle-differential cross sections. On 28-Aug-2026 the two OTHER never-executed functions in this module,
        `computePartialCrossSectionUnpolarized` and `computeStatisticalTensorUnpolarized`, were run for the first time
        and held THREE defects between them: a retired call signature, a Clebsch-Gordan coefficient evaluated where it
        had to be skipped, and a normalisation that still does not sum to the total. Expect the same class here until
        somebody runs it, and do not quote a number from it before then.
"""
function computeDisplayNonE1AngleDifferentialCS(stream::IO, lines::Array{PhotoIonization.Line,1},
                                                      settings::PhotoIonization.Settings)
    function spinDensityMatrix(lambda1::Int64, lambda2::Int64, stokes::ExpStokes)
        # Convert the Stokes parameters of the incoming light into a spin-density matrix on the indices lambda = +-1
        if      lambda1 == lambda2  == 1               return( (1.0 + stokes.P3)/2. )
        elseif  lambda1 ==  1   &&   lambda2  == -1    return( (stokes.P1 - stokes.P1*im)/2. )
        elseif  lambda1 == -1   &&   lambda2  ==  1    return( (stokes.P1 + stokes.P1*im)/2. )
        elseif  lambda1 == lambda2  == -1              return( (1.0 - stokes.P3)/2. )
        else    error("stop a")
        end
    end

    angCS = Tuple{Float64, ComplexF64, ComplexF64}[]  # theta, angCs.Coulomb, angCs.Babushkin
    # Loop about all lines; a table is printed independently for each line
    for  line in lines
        # Flatten once: one entry per (partial wave, total symmetry, multipole), each carrying both gauges.
        entries = Tuple{Int64,LevelSymmetry,EmMultipole,EmPropertyC}[]
        for  pw in line.partialWaves,  ch in pw.channels,  ma in ch.amplitudes
            push!(entries, (pw.kappa, ch.symmetry, ma.multipole, ma.amplitude))
        end
        # Loop over all angles theta
        for  theta in settings.thetas
            cs = EmPropertyC(0.0im)
            for  (kapa, syma, mpa, ampa)  in entries,   (kapb, symb, mpb, ampb)  in entries
                s1 = Subshell(20, kapa);   j1 = Basics.subshell_j(s1)
                s2 = Subshell(20, kapb);   j2 = Basics.subshell_j(s2)
                for  X = 0:20  # Test for triangular conditions for X and continue otherwise
                    if  AngularMomentum.isTriangle(mpa.L, mpb.L, X)                    &&
                        AngularMomentum.isTriangle(syma.J,  symb.J, AngularJ64(X) )    &&
                        AngularMomentum.isTriangle(j1,  j2, AngularJ64(X) )
                        K = PhotoIonization.angularFunctionK(mpa.L, mpb.L, X, line.initialLevel.J, line.finalLevel.J,
                                                             kapa, syma.J, kapb, symb.J)
                        # Compute the summation over lambda's and mu's
                        W = 0.
                        for  lambda1 = -1:2:1,   lambda2 = -1:2:1,   mu = -1//2:1: 1//2
                            W = W + spinDensityMatrix(lambda1, lambda2, settings.stokes) *
                                    PhotoIonization.angularFunctionW(theta, mpa.L, mpb.L, X, lambda1, lambda2,
                                                                     kapa, mu, kapb, mu) *
                                    (1.0im)^(mpb.L - mpa.L) * (lambda1*lambda2) / 2. *
                                     AngularMomentum.phaseMultipole(1.0im*lambda1, mpa) *
                                     AngularMomentum.phaseMultipole(1.0im*lambda2, mpb)
                        end
                        # No gauge bookkeeping: the product is componentwise, and a magnetic amplitude has equal components and so enters
                        # both sums by itself.
                        cs = cs + K * W * (ampa * conj(ampb))
                    end
                end
            end
            push!(angCS, (theta, cs.Coulomb, cs.Babushkin) )
        end

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

    end

    return( nothing )
end


"""
`PhotoIonization.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model,
        grid::Radial.Grid, settings::PhotoIonization.Settings; output::Bool=true)`
    ... computes the photoionization amplitudes and all requested properties for the lines between the levels of the two given multiplets
        and for every photon energy of the settings; this is the standard entry point of this module. The nuclear potential and the B-spline
        basis are built once and shared by all lines, which are then computed in parallel. An Array{PhotoIonization.Line,1} is returned if
        output=true, and nothing otherwise.
"""
function computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model,
                            grid::Radial.Grid, settings::PhotoIonization.Settings; output::Bool=true)
    println("")
    printstyled("PhotoIonization.computeLines(): The computation of photoionization amplitudes starts now ... \n", color=:light_green)
    printstyled("-------------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    lines = PhotoIonization.determineLines(finalMultiplet, initialMultiplet, settings)
    if  settings.printBefore    PhotoIonization.displayLines(stdout, lines)    end
    maxEnergy = 0.;   for  line in lines   maxEnergy = max(maxEnergy, line.electronEnergy)   end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
    nuclearPot = Nuclear.nuclearPotential(nm, grid)
    primitives = Bsplines.generatePrimitives(grid)
    newLines = Vector{PhotoIonization.Line}(undef, length(lines))
    Threads.@threads for  i  in  eachindex(lines)
        newLines[i] = PhotoIonization.computeAmplitudesProperties(lines[i], nm, grid, nrContinuum, settings;
                                                                        nuclearPot=nuclearPot, primitives=primitives)
    end
    PhotoIonization.displayPhases(newLines)
    PhotoIonization.displayResults(stdout, newLines, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   PhotoIonization.displayResults(iostream, newLines, settings)   end
    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`PhotoIonization.computeLinesCascade(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                     settings::PhotoIonization.Settings, initialLevelSelection::LevelSelection; 
                                     output=true, printout::Bool=true)`  
    ... to compute the photoionization transition amplitudes and all properties as requested by the given settings. The computations and
        printout is adapted for large cascade computations by including only lines with at least one channel and by sending all printout to
        a summary file only. A list of lines::Array{PhotoIonization.Lines} is returned.
"""
function  computeLinesCascade(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                              settings::PhotoIonization.Settings, initialLevelSelection::LevelSelection; output=true, printout::Bool=true)
    
    lines = PhotoIonization.determineLines(finalMultiplet, initialMultiplet, settings)
    # Display-only and the same for every line; set ONCE per computation, never from inside the line loop.
    Defaults.setStandardSubshellList(Basics.generate(OrderedSubshellList(), finalMultiplet.levels[1].basis,
                                                     initialMultiplet.levels[1].basis); printout=false)
    # Display all selected lines before the computations start if  settings.printBefore    PhotoIonization.displayLines(stdout, lines)
    # end Determine maximum energy and check for consistency of the grid
    maxEnergy = 0.;   for  en in settings.photonEnergies     maxEnergy = max(maxEnergy, Defaults.convertUnits("energy: to atomic", en))   end
                      for  en in settings.electronEnergies   maxEnergy = max(maxEnergy, Defaults.convertUnits("energy: to atomic", en))   end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
    # Calculate all amplitudes and requested properties Constants of the whole computation; see the note in computeLines above.
    nuclearPot = Nuclear.nuclearPotential(nm, grid)
    primitives = Bsplines.generatePrimitives(grid)
    # Threaded as computeLines is, but this loop SKIPS lines whose initial level is not selected, so the results are collected into a
    # Union{Nothing,Line} vector by index and the gaps are dropped afterwards -- the same shape ImpactExcitation.computeLinesCascade uses,
    # and for the same reason: push! is neither thread-safe nor order-preserving.
    tmpLines = Vector{Union{Nothing, PhotoIonization.Line}}(nothing, length(lines))
    doPrint  = printout  &&  Threads.nthreads() == 1
    Threads.@threads for  i  in  eachindex(lines)
        if  doPrint  &&  rem(i,10) == 0    println("> Photo line $i:  ... calculated ")    end
        # Do not compute line if initial level is not in initialLevelSelection()
        if  !Basics.selectLevel(lines[i].initialLevel, initialLevelSelection)   continue
        end

        tmpLines[i] = PhotoIonization.computeAmplitudesProperties(lines[i], nm, grid, nrContinuum, settings;
                                                                  printout=printout, nuclearPot=nuclearPot,
                                                                  primitives=primitives)
    end
    newLines = PhotoIonization.Line[ l  for l in tmpLines  if l !== nothing ]
    # Print all results to a summary file, if requested
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   PhotoIonization.displayResults(iostream, newLines, settings)   end

    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`PhotoIonization.computeLinesPlasma(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                                    settings::PhotoIonization.PlasmaSettings, plasmaModel::Basics.AbstractPlasmaModel; output::Bool=true)`
    ... to compute the photoIonization transition amplitudes and all properties as requested by the given settings and plasma model. A list
        of lines::Array{PhotoIonization.Lines} is returned.
"""
function  computeLinesPlasma(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                                settings::PhotoIonization.PlasmaSettings, plasmaModel::Basics.AbstractPlasmaModel; output::Bool=true)
    println("")
    printstyled("PhotoIonization.computeLinesPlasma(): The computation of photo-ionization cross sections starts now ... \n", color=:light_green)
    printstyled("------------------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    photoSettings = PhotoIonization.Settings(settings.multipoles, settings.gauges, settings.photonEnergies, Float64[], Float64[],
                        false, false, false, false, false, false, settings.printBefore, settings.lineSelection,
                        Basics.ExpStokes(), 0., [0,1,2,3,4,5])

    lines = PhotoIonization.determineLines(finalMultiplet, initialMultiplet, photoSettings)
    # Display-only and the same for every line; set ONCE per computation, never from inside the line loop.
    Defaults.setStandardSubshellList(Basics.generate(OrderedSubshellList(), finalMultiplet.levels[1].basis,
                                                     initialMultiplet.levels[1].basis); printout=false)
    # Display all selected lines before the computations start
    if  settings.printBefore    PhotoIonization.displayLines(stdout, lines)    end
    # Calculate all amplitudes and requested properties Constants of the whole computation; see the note in computeLines above.
    nuclearPot = Nuclear.nuclearPotential(nm, grid)
    primitives = Bsplines.generatePrimitives(grid)
    newLines = Vector{PhotoIonization.Line}(undef, length(lines))
    Threads.@threads for  i  in  eachindex(lines)
        newLines[i] = PhotoIonization.computeAmplitudesPropertiesPlasma(lines[i], nm, grid, settings, plasmaModel;
                                                                        nuclearPot=nuclearPot, primitives=primitives)
    end
    # Print all results to screen
    PhotoIonization.displayResults(stdout, newLines, photoSettings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   PhotoIonization.displayResults(iostream, newLines, photoSettings)     end

    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`PhotoIonization.computePartialCrossSectionUnpolarized(Mf::AngularM64, line::PhotoIonization.Line)`
    ... computes the partial photoionization cross section for a given magnetic quantum number Mf of the final level; an EmProperty is
        returned, holding BOTH gauges at once.

        RUN FOR THE FIRST TIME ON 28-Aug-2026, and it took two fixes to get there -- see the caller at
        `displayResults`, which was still passing the retired flat signature, and the |M| > J guard in `Racahexpr`
        below. Both are corrected.

        IT DOES NOT YET SUM TO THE TOTAL, and this is the open part. Measured on Ne-like Kr at 3000 eV: with E1
        and M1 only, `sum over M_f` divided into `line.crossSection` gives 2086.5406 (J_f = 3/2) and 2086.5404
        (J_f = 1/2), which is 1/(9 alpha^2) to eight figures. That looked like a clean missing scalar and IS NOT
        one: adding E2 and M2 moves the two ratios to 2063.9997 and 2065.5892, so they are no longer equal to each
        other and no longer 1/(9 alpha^2). The 9 is an artefact of dipole-only channels. What IS certain is the
        alpha: this factor carries `* alpha` while `computeCrossSection` carries `/ alpha`. Do not "fix" this by
        multiplying in 9 alpha^2 -- the L-dependence has to be derived.

        NO GAUGE CAN BE SELECTED WRONGLY HERE. The product `ma.amplitude * conj(mp.amplitude)` of two EmPropertyC
        is componentwise, so Coulomb meets Coulomb and Babushkin meets Babushkin whatever anyone writes, and one
        call returns both.
"""
function computePartialCrossSectionUnpolarized(Mf::AngularM64, line::PhotoIonization.Line)
    function Racahexpr(kappa::Int64, Ji::AngularJ64, Jf::AngularJ64, Mf::AngularM64, J::AngularJ64, Jp::AngularJ64,
                        L::Int64, Lp::Int64, p::Int64, pp::Int64)
        t1 = Basics.oplus( AngularJ64(Lp), Jf);    t2 = Basics.oplus( AngularJ64(L), Jf);    tList = intersect(t1, t2)
        wb = 0.
        for  t  in tList
            for  lambda = -1:2:1
                j = AngularMomentum.kappa_j(kappa);   Mf_lambda = Basics.add(AngularM64(lambda), Mf)
                # A Clebsch-Gordan coefficient with |M| > J is ZERO BY DEFINITION and must be skipped, not
                # evaluated: WignerSymbols raises DomainError rather than returning 0. The case is reached
                # whenever the coupled t is smaller than |M_f + lambda| -- e.g. J_f = 3/2, M_f = 3/2, lambda = +1
                # asks for t = 1/2 with M = 5/2. Never seen before 28-Aug-2026 because this function had never run.
                if  abs(Basics.twice(Mf_lambda)) > Basics.twice(t)    continue    end
                wb = wb + (1.0im * lambda)^p * (-1.0im * lambda)^pp *
                        AngularMomentum.ClebschGordan( AngularJ64(Lp), AngularM64(lambda), Jf, Mf, t, Mf_lambda) *
                        AngularMomentum.ClebschGordan( AngularJ64(L),  AngularM64(lambda), Jf, Mf, t, Mf_lambda) *
                        AngularMomentum.Wigner_9j(j, Jp, Jf, J, Ji, AngularJ64(L), Jf, AngularJ64(Lp), t)
            end
        end
        return( wb )
    end

    Ji = line.initialLevel.J;    Jf = line.finalLevel.J
    waC = 0.0im;    waB = 0.0im
    # Pairs are formed WITHIN one partial wave; two amplitudes of different kappa never multiply each other.
    for  pw in line.partialWaves
        for  cha in pw.channels,  ma in cha.amplitudes
            J = cha.symmetry.J;    L = ma.multipole.L;    p  = ma.multipole.electric ? 1 : 0
            for  chp in pw.channels,  mp in chp.amplitudes
                Jp = chp.symmetry.J;   Lp = mp.multipole.L;   pp = mp.multipole.electric ? 1 : 0
                wc = (1.0im)^(L - Lp) * (-1)^(L + Lp) * AngularMomentum.bracket([AngularJ64(L), AngularJ64(Lp), J, Jp]) *
                        Racahexpr(pw.kappa, Ji, Jf, Mf, J, Jp, L, Lp, p, pp) * (ma.amplitude * conj(mp.amplitude))
                waC = waC + wc.Coulomb;    waB = waB + wc.Babushkin
            end
        end
    end
    # THIS FUNCTION DOES NOT REPRODUCE THE TOTAL, AND THE DISCREPANCY IS NOT A CONSTANT. Measured 30-Aug-2026,
    # E1, omega = 80 eV, summing this function over all M_f and dividing by PhotoIonization.computeCrossSection --
    # which IS validated, against Stobbe's closed form (see the note at computeCrossSection):
    #
    #     system          J_i    J_f                 sum/total, in units of alpha^2
    #     H 1s            1/2    0                   5.0002
    #     Ne 2p^6         0      3/2, 1/2            9.0000, 9.0000
    #     Ne+ 2p^5        3/2    2, 1, 0, 2, 0       3.8382, 4.0199, 4.2463, 3.0982, 4.2469
    #     Ne+ 2p^5        1/2    2, 1, 0, 2, 0       4.5545, 5.7424, 5.8912, 4.5544, 5.8917
    #
    # TWO SEPARATE FAULTS, and only the first is understood. (i) The alpha below sits in the NUMERATOR where
    # computeCrossSection has it in the DENOMINATOR. That one is certainly wrong: JAC's multipole amplitude already
    # carries the photon-momentum factor from j_L(alpha omega r), and the total's 1/alpha is there to cancel it --
    # the partial uses the SAME amplitudes, so it needs the same cancellation. alpha^2 is accordingly common to
    # every row above. (ii) The residue after removing alpha^2 is NOT constant: it varies from 3.10 to 9.00, and
    # TWO LINES WITH THE SAME J_i AND THE SAME J_f differ (3.8382 against 3.0982). No prefactor in J_i can absorb
    # that, so the ANGULAR expression above is wrong as well, and by an amount that is not established.
    #
    # DO NOT "FIX" THIS BY MULTIPLYING BY 9 alpha^2. That number is the Ne value only, and Ne is the special case:
    # J_i = 0 with a single initial level. Patching the call site with a compensating factor is exactly the drift
    # Rule 18 exists to stop. The correct prefactor and the correct angular expression have to be derived together
    # and checked against a J_i != 0 case; that is carried as challenge 63.
    #
    # NOTHING IN THE PACKAGE CONSUMES THIS. It has never been used in an application nor checked against a known
    # value, and computeCrossSection is unaffected -- so the numbers this returns are the only ones at risk.
    csFactor = 8 * pi^3 * Defaults.getDefaults("alpha") / (2*line.photonEnergy * (Basics.twice(Ji) + 1))
    Defaults.warn(AddWarning(), "PhotoIonization.computePartialCrossSectionUnpolarized(): these partial cross " *
                  "sections do NOT sum to the total (they are low by a level-dependent 3.1-9.0 times alpha^2); " *
                  "see the note at this function and challenge 63. Do not use them quantitatively.")

    return( EmProperty(real(csFactor * waC), real(csFactor * waB)) )
end


"""
`PhotoIonization.computeStatisticalTensorUnpolarized(k::Int64, q::Int64, line::PhotoIonization.Line,
                                                           settings::PhotoIonization.Settings)`
    ... computes the statistical tensor rho_kq of the photoion; an EmPropertyC is returned, holding both gauges. As above, no gauge argument
        is taken because the gauge is carried by the amplitude.

        CARRIES THE SAME CORRECTION as computePartialCrossSectionUnpolarized: the flat version's inner loop tests `cha.gauge` where it means
        `chp.gauge`, so it multiplies Coulomb amplitudes with Babushkin ones. That cannot happen here. The two versions therefore differ by
        more than rounding, and the difference is the point rather than a defect of this one.

        A GREEN SUITE SAYS NOTHING ABOUT THIS FUNCTION UNDER LINEARLY POLARIZED LIGHT, and that is worth knowing before trusting it. The
        Clebsch-Gordan coupling the two helicities carried the projection -lambda where it needs -lambdap, so both projections summed to zero
        and every lambda /= lambda' term -- exactly the terms carrying LINEAR polarization -- was deposited at q = 0 rather than at
        q = +-2. Under P1 that cancelled the alignment exactly, A_20 coming out 0.000000 instead of 0.128571; under P2 it produced an
        imaginary orientation that Hermiticity forbids; and q = +-2 was always zero. For UNPOLARIZED and CIRCULAR light those terms vanish
        outright, so the wrong expression was RIGHT there -- which is why every approved case passed throughout. The defect lived precisely
        in the configurations nobody had approved, so the next reader should not take a green suite as coverage of this one.
"""
function computeStatisticalTensorUnpolarized(k::Int64, q::Int64, line::PhotoIonization.Line,
                                                   settings::PhotoIonization.Settings)
    Ji = line.initialLevel.J;    Jf = line.finalLevel.J
    P1 = settings.stokes.P1;     P2 = settings.stokes.P2;     P3 = settings.stokes.P3
    waC = 0.0im;    waB = 0.0im
    for  pw in line.partialWaves
        j = AngularMomentum.kappa_j(pw.kappa)
        for  cha in pw.channels,  ma in cha.amplitudes
            J = cha.symmetry.J;    L = ma.multipole.L;    p  = ma.multipole.electric ? 1 : 0
            for  chp in pw.channels,  mp in chp.amplitudes
                Jp = chp.symmetry.J;   Lp = mp.multipole.L;   pp = mp.multipole.electric ? 1 : 0
                for  lambda = -1:2:1
                    for  lambdap = -1:2:1
                        if  lambda == lambdap   wb = (1.0 + 0.0im + lambda*P3)    else    wb = P1 - lambda * P2 * im    end
                        wc = wb * (1.0im)^(L - Lp + p - pp) * lambda^p * lambdap^pp *
                                sqrt( AngularMomentum.bracket([AngularJ64(L), AngularJ64(Lp), J, Jp]) ) *
                                AngularMomentum.phaseFactor([J, +1, Jp, +1, Jf, +1, Ji, +1, j, +1, AngularJ64(1)]) *
                                AngularMomentum.ClebschGordan( AngularJ64(L),  AngularM64(lambda), AngularJ64(Lp),
                                                                AngularM64(-lambdap), AngularJ64(k),  AngularM64(q)) *
                                AngularMomentum.Wigner_6j(Jf, j, Jp, J, AngularJ64(k), Jf) *
                                AngularMomentum.Wigner_6j(Jp, Ji, AngularJ64(Lp), AngularJ64(L), AngularJ64(k), J) *
                                (ma.amplitude * conj(mp.amplitude))
                        waC = waC + wc.Coulomb;    waB = waB + wc.Babushkin
                    end
                end
            end
        end
    end
    fc = pi / (Basics.twice(Ji) + 1)

    return( EmPropertyC(fc * waC, fc * waB) )
end


"""
`PhotoIonization.computeTimeDelays(partialWaves::Array{PhotoIonization.PartialWave,1},
        xPartialWaves::Array{PhotoIonization.PartialWave,1}, deltaE::Float64, Jf::AngularJ64)`
    ... computes the coherent and incoherent time delays from two sets of partial waves that have been computed at energies differing by
        deltaE; a tuple (coherentDelay, incoherentDelay)::Tuple{EmProperty,EmProperty} is returned. The phase of each partial wave is taken
        from the wave itself, which is what it is a property of.

        ONLY THE INCOHERENT DELAY IS EVALUATED: it is the amplitude-weighted mean of the partial-wave phases, differenced over deltaE.

        The COHERENT delay is NOT evaluated and EmProperty(0.) is returned in its place, so a caller must not read it as a computed zero. No
        general expression for it is available here; a formula restricted to a single initial subshell would look general once written into
        this return value, and that is why none is given. Completing it is an open task of this function.
"""
function computeTimeDelays(partialWaves::Array{PhotoIonization.PartialWave,1},
                                 xPartialWaves::Array{PhotoIonization.PartialWave,1},
                                 deltaE::Float64, Jf::AngularJ64)
    function weightedPhase(pws)
        nom = EmProperty(0., 0.);    den = EmProperty(0., 0.)
        for  pw in pws,  ch in pw.channels,  ma in ch.amplitudes
            w   = abs2(ma.amplitude)                       # an EmProperty: |amp|^2 in both gauges
            nom = nom + w * pw.phase
            den = den + w
        end
        return( EmProperty(nom.Coulomb / den.Coulomb, nom.Babushkin / den.Babushkin) )
    end
    Deff  = weightedPhase(partialWaves);    Deffx = weightedPhase(xPartialWaves)
    incoherentDelay = EmProperty( (Deffx.Coulomb - Deff.Coulomb) / deltaE, (Deffx.Babushkin - Deff.Babushkin) / deltaE )
    coherentDelay   = EmProperty(0.)       # not evaluated; see the docstring

    return( coherentDelay, incoherentDelay )
end


"""
`PhotoIonization.determineChannels(finalLevel::Level, initialLevel::Level, settings::PhotoIonization.Settings)`
    ... determines the partial waves of the outgoing electron and, within each of them, the scattering channels of the given transition, by
        applying AngularMomentum.allowedMultipoleSymmetries and AngularMomentum.allowedKappaSymmetries. The distinct pairs (kappa, symt) are
        collected first, and every multipole that reaches such a pair is then attached to it; the gauge is not iterated over at all, since
        an electric multipole yields one EmPropertyC holding both gauges and a magnetic one an EmPropertyC with equal components. An
        Array{PhotoIonization.PartialWave,1} is returned, with every amplitude still zero.
"""
function determineChannels(finalLevel::Level, initialLevel::Level, settings::PhotoIonization.Settings)
    symi = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity)
    # (1) Collect, for every physical pair (kappa, symt), the multipoles that can reach it.
    mpsFor = Dict{Tuple{Int64,LevelSymmetry}, Array{EmMultipole,1}}()
    order  = Tuple{Int64,LevelSymmetry}[]                    # to keep a reproducible sequence
    for  mp in settings.multipoles
        for  symt in AngularMomentum.allowedMultipoleSymmetries(symi, mp)
            for  kappa in AngularMomentum.allowedKappaSymmetries(symt, symf)
                if  !(Basics.subshell_l(Subshell(10,kappa)) in settings.lValues)    continue     end
                key = (kappa, symt)
                if  haskey(mpsFor, key)   push!(mpsFor[key], mp)
                else                      mpsFor[key] = EmMultipole[mp];   push!(order, key)
                end
            end
        end
    end
    # (2) Group the pairs by kappa: one partial wave per kappa, carrying the symmetries it serves.
    kappas = Int64[];   for (kappa, symt) in order    if !(kappa in kappas)   push!(kappas, kappa)   end    end
    partialWaves = PartialWave[]
    for  kappa in kappas
        channels = Channel[]
        for  (ka, symt) in order
            ka == kappa   ||   continue
            amps = MultipoleAmplitude[]
            for  mp in mpsFor[(ka, symt)]
                push!(amps, MultipoleAmplitude(mp, EmPropertyC(Complex(0.), Complex(0.))))
            end
            push!(channels, Channel(symt, amps))
        end
        push!(partialWaves, PartialWave(kappa, 0., 0., channels))
    end

    return( partialWaves )
end


"""
`PhotoIonization.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::PhotoIonization.Settings)`
    ... selects the photoionization lines to be computed by forming all pairs of an initial and a final level that pass the line
        selection, once for every photon energy requested in the settings, and assigns the free-electron energy that follows from the
        transition energy and settings.freeElectronShift. An Array{PhotoIonization.Line,1} is returned, with all amplitudes still zero.
"""
function determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::PhotoIonization.Settings)
    lines    = PhotoIonization.Line[]
    shift_au = Defaults.convertUnits("energy: to atomic", settings.freeElectronShift)
    for  iLevel  in  initialMultiplet.levels
        for  fLevel  in  finalMultiplet.levels
            if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                for  omega in settings.photonEnergies
                    omega_au = Defaults.convertUnits("energy: to atomic", omega)
                    energy   = omega_au - (fLevel.energy - iLevel.energy) + shift_au
                    if  energy < 0.    continue   end
                    pws = PhotoIonization.determineChannels(fLevel, iLevel, settings)
                    push!( lines, PhotoIonization.Line(iLevel, fLevel, energy, omega_au, EmProperty(0.),
                                                             EmProperty(0.), EmProperty(0.), EmProperty(0.), pws) )
                end
                for  en in settings.electronEnergies
                    energy_au = Defaults.convertUnits("energy: to atomic", en) + shift_au
                    omega     = energy_au + (fLevel.energy - iLevel.energy)
                    if  energy_au < 0.    continue   end
                    pws = PhotoIonization.determineChannels(fLevel, iLevel, settings)
                    push!( lines, PhotoIonization.Line(iLevel, fLevel, energy_au, omega, EmProperty(0.),
                                                             EmProperty(0.), EmProperty(0.), EmProperty(0.), pws) )
                end
            end
        end
    end

    return( lines )
end


#################################################################################################################################
#################################################################################################################################

"""
`PhotoIonization.displayLineData(stream::IO, lines::Array{PhotoIonization.Line,1})`
    ... lists the initial levels of the given lines and, for each of them, the photon energies and the total cross sections that have been
        obtained, so that the data of a photoionization computation can be inspected in compact form. A neat table is printed to stream;
        nothing::Nothing is returned.
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

    # Extract and display photon energies for which line data are available
    omegas = Float64[]
    for  line  in  lines    push!(omegas, line.photonEnergy)   end;    
    omegas = unique(omegas);    omegas = sort(omegas)
    sa  = "    ";   
    for   omega in omegas  sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", omega)) * "   "      end

    println(stream, "\n  Photon energies $(TableStrings.inUnits("energy")) for which given photoionization line data are given:")
    println(stream, "\n" * sa)

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

    if  length(lines) == 0    println("  >>> No photoionization data in lines here ... return ");    return(nothing)    end 

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
    ... lists the photoionization lines that have been selected, together with the free-electron energy and the multipoles and gauges that
        will be computed for each of them. One row is printed per (multipole, gauge), an expansion made here because a single amplitude
        holds both gauges: the gauge is a property of this table, not of the partial wave. A neat table is printed to stream;
        nothing::Nothing is returned.
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
        # One row per (multipole, GAUGE), which is a property of the PRESENTATION: an electric multipole is one amplitude holding two
        # gauges, a magnetic one has no gauge freedom.  Expanding it here keeps the table exactly as it was when the channel itself carried
        # a gauge label.
        for  pw in line.partialWaves,  ch in pw.channels,  ma in ch.amplitudes
            if  string(ma.multipole)[1] == 'E'
                push!( kappaMultipoleSymmetryList, (pw.kappa, ma.multipole, Basics.Coulomb,   ch.symmetry) )
                push!( kappaMultipoleSymmetryList, (pw.kappa, ma.multipole, Basics.Babushkin, ch.symmetry) )
                nchannels = nchannels + 2
            else
                push!( kappaMultipoleSymmetryList, (pw.kappa, ma.multipole, Basics.Magnetic,  ch.symmetry) )
                nchannels = nchannels + 1
            end
        end
        wa = TableStrings.kappaMultipoleSymmetryTupels(85, kappaMultipoleSymmetryList)
        sb = sa * wa[1];    println(stream,  sb )  
        for  i = 2:length(wa)
            sb = TableStrings.hBlank( length(sa) ) * wa[i];    println(stream,  sb )
        end
    end
    println(stream, "  ", TableStrings.hLine(nx), "\n")
    println(stream, "  A total of $nchannels channels need to be calculated. \n")

    return( nothing )
end


"""
`PhotoIonization.displayPhases(lines::Array{PhotoIonization.Line,1})`
    ... lists the selected photoionization lines together with the kappa and the phase of every outgoing partial wave. The phase is a
        property of the partial wave and is taken from it directly, so that each kappa contributes a single value. A neat table is printed
        to stdout; nothing::Nothing is returned.
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
        # As in displayLines, one row per (multipole, gauge); the phase is read from the partial wave.
        for  pw in line.partialWaves,  ch in pw.channels,  ma in ch.amplitudes
            if  string(ma.multipole)[1] == 'E'
                push!( kappaMultipoleSymmetryPhaseList, (pw.kappa, ma.multipole, Basics.Coulomb,   ch.symmetry, pw.phase) )
                push!( kappaMultipoleSymmetryPhaseList, (pw.kappa, ma.multipole, Basics.Babushkin, ch.symmetry, pw.phase) )
                nchannels = nchannels + 2
            else
                push!( kappaMultipoleSymmetryPhaseList, (pw.kappa, ma.multipole, Basics.Magnetic,  ch.symmetry, pw.phase) )
                nchannels = nchannels + 1
            end
        end
        wa = TableStrings.kappaMultipoleSymmetryPhaseTupels(85, kappaMultipoleSymmetryPhaseList)
        sb = sa * wa[1];    println( sb )  
        for  i = 2:length(wa)
            sb = TableStrings.hBlank( length(sa) ) * wa[i];    println( sb )
        end
    end
    println("  ", TableStrings.hLine(nx), "\n")
    println("  A total of $nchannels channels has been calculated. \n")

    return( nothing )
end


"""
`PhotoIonization.displayResults(stream::IO, lines::Array{PhotoIonization.Line,1}, settings::PhotoIonization.Settings)`
    ... lists the total photoionization cross sections for initially unpolarized atoms and unpolarized plane-wave photons, in both gauges,
        and appends whatever further tables the settings request: the angular beta parameters, the time delays, the partial cross sections
        and the reduced statistical tensors of the photoion. This is the main result table of the module. A neat table is printed to stream;
        nothing::Nothing is returned.
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

    wx = 1.0 # wx = 1.0  (Schippers, August'23; wx = 2.0; wx = pi/2)
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
        for  pw in line.partialWaves,  ch in pw.channels,  ma in ch.amplitudes
            multipoles = push!( multipoles, ma.multipole)
        end
        multipoles = unique(multipoles);   mpString = TableStrings.multipoleList(multipoles) * "          "
        sa = sa * TableStrings.flushleft(11, mpString[1:10];  na=2)
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", wx * line.crossSection.Coulomb))     * "    "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", wx * line.crossSection.Babushkin))   * "                 "
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))


    # Total (summed) cross sections.  The table above is resolved into the individual final levels f, which is rarely what one wants to
    # quote; the sum over f at a fixed initial level and photon energy is.
    PhotoIonization.displayTotalCrossSections(stream, lines, settings)


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
            for  pw in line.partialWaves,  ch in pw.channels,  ma in ch.amplitudes
                multipoles = push!( multipoles, ma.multipole)
            end
            multipoles = unique(multipoles);   mpString = TableStrings.multipoleList(multipoles) * "          "
            sa = sa * TableStrings.flushleft(11, mpString[1:10];  na=2)
            println(stream, sa)
            MfList = Basics.projections(line.finalLevel.J)
            for  Mf in MfList
                sb  = TableStrings.hBlank(97)
                # ONE call returns BOTH gauges, as it does for the statistical tensor below. This site still passed a
                # `gauge` and read `.re` -- the FLAT signature, retired on 15-Aug-2026 -- so it raised a MethodError
                # the moment `calcPartialCs` was switched on, which nobody had ever done.
                wa  = PhotoIonization.computePartialCrossSectionUnpolarized(Mf, line)
                sb  = sb * TableStrings.flushright( 8, string(Mf))                             * "       "
                sb  = sb * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", wa.Coulomb))   * "    "
                sb  = sb * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", wa.Babushkin)) * "    "
                println(stream, sb)
            end
        end
        println(stream, "  ", TableStrings.hLine(nx))
    end


    if  settings.calcTensors  
        nx = 178 
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
        sa = sa * TableStrings.center(52, "Cou --  rho_kq (J_f)  (re, im)  -- Bab"; na=3)      
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 

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
            for  pw in line.partialWaves,  ch in pw.channels,  ma in ch.amplitudes
                multipoles = push!( multipoles, ma.multipole)
            end
            multipoles = unique(multipoles);   mpString = TableStrings.multipoleList(multipoles) * "          "
            sa = sa * TableStrings.flushleft(11, mpString[1:10];  na=2)
            println(stream, sa)
            for  k = 0:2
                for q = -k:k
                    sb   = TableStrings.hBlank(102)
                    # One call returns both gauges, and the component is COMPLEX: with linearly polarized light the
                    # photoion carries coherence between sublevels two units apart, so q = +-2 has an imaginary part
                    # whenever P2 does not vanish.  Both parts are printed.
                    rho  = PhotoIonization.computeStatisticalTensorUnpolarized(k, q, line, settings)
                    if  abs(rho.Coulomb) == abs(rho.Babushkin) == 0.    continue    end
                    cnv(x) = Defaults.convertUnits("cross section: from atomic", x)
                    sb   = sb * string(k) * " " * TableStrings.flushright( 4, string(q))             * "       "
                    sb   = sb * @sprintf("%.5e %+.5e", cnv(real(rho.Coulomb)),   cnv(imag(rho.Coulomb)))    * "    "
                    sb   = sb * @sprintf("%.5e %+.5e", cnv(real(rho.Babushkin)), cnv(imag(rho.Babushkin)))  * "    "
                    println(stream, sb)
                end
            end
        end
        println(stream, "  ", TableStrings.hLine(nx))
    end

    return( nothing )
end


"""
`PhotoIonization.displayTotalCrossSections(stream::IO, lines::Array{PhotoIonization.Line,1}, settings::PhotoIonization.Settings)`
    ... lists the total photoionization cross sections summed over all final levels, i.e. one entry per initial level and photon energy,
        which is the quantity to be compared with a measured total cross section. Nothing is printed if no line was selected. A neat table
        is printed to stream; nothing::Nothing is returned.
"""
function  displayTotalCrossSections(stream::IO, lines::Array{PhotoIonization.Line,1}, settings::PhotoIonization.Settings)
    if  length(lines) == 0      return( nothing )    end
    nx = 96
    println(stream, " ")
    println(stream, "  Total photoionization cross sections, summed over all final levels:")
    println(stream, " ")
    println(stream, "    Each row sums the line-resolved cross sections above over ALL final levels f that belong to")
    println(stream, "    the given initial level i and photon energy, and over the partial waves of the photo-electron.")
    println(stream, "    It is the direct (non-resonant) cross section only: no resonant excitation-autoionization, no")
    println(stream, "    shake-off and no multiple ionization are contained.  If the lines cover more than one ionized")
    println(stream, "    subshell, this sum runs over those subshells as well and is then the grand total.")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(12, "i-level"    ; na=2);                        sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(10, "i--J^P"     ; na=3);                        sb = sb * TableStrings.hBlank(13)
    sa = sa * TableStrings.center(12, "omega"      ; na=4)
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(10, "No lines"   ; na=3);                        sb = sb * TableStrings.hBlank(13)
    sa = sa * TableStrings.center(30, "Cou -- Total cross section -- Bab"; na=3)
    sb = sb * TableStrings.center(30, TableStrings.inUnits("cross section") * "          " *
                                            TableStrings.inUnits("cross section"); na=3)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))

    wx = 1.0
    # Collect the distinct (initial level, photon energy) pairs in the order in which they first occur
    keys = Tuple{Int64,Float64}[]
    for  line in lines
        wa = (line.initialLevel.index, line.photonEnergy)
        if  !(wa in keys)   push!(keys, wa)    end
    end
    for  (idx, omega)  in  keys
        tcs = Basics.EmProperty(0.);   nl = 0;   isym = LevelSymmetry(AngularJ64(0), Basics.plus)
        for  line in lines
            if  line.initialLevel.index == idx  &&  line.photonEnergy == omega
                tcs = tcs + line.crossSection;    nl = nl + 1
                isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
            end
        end
        sa = "  " * TableStrings.center(12, string(idx); na=2)
        sa = sa * TableStrings.center(10, string(isym); na=3)
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", omega))  * "    "
        sa = sa * TableStrings.center(10, string(nl); na=3)
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", wx * tcs.Coulomb))    * "    "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", wx * tcs.Babushkin))
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))

    return( nothing )
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
    ... to extract from lines the total PI cross section that refer to the given omega and initial level and to to the ionization of an
        electron from shell; a cross section cs::EmProperty is returned.
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
`PhotoIonization.extractLines(lines::Array{PhotoIonization.Line,1}, omega::Float64)`  
    ... to extract from lines all those that refer to the given omega;
        a reduced list rLines::Array{PhotoIonization.Line,1} is returned.
"""
function  extractLines(lines::Array{PhotoIonization.Line,1}, omega::Float64)
    rLines = PhotoIonization.Line[]
    for  line  in  lines    
        if  line.photonEnergy == omega   push!(rLines, line)  end
    end
    
    return( rLines )
end


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
`PhotoIonization.getLineKappas(line::PhotoIonization.Line)`  
    ... returns a list of kappa-values (partial waves) which contribute to the given line, to which one or several channels are assigned. An
        kappaList::Array{Int64,1} is returned.
"""
function getLineKappas(line::PhotoIonization.Line)
    # There is one partial wave per kappa, so the list is distinct by construction.
    return( Int64[pw.kappa  for pw in line.partialWaves] )
end


"""
`PhotoIonization.interpolateCrossSection(lines::Array{PhotoIonization.Line,1}, omega::Float64, initialLevel)`  
    ... to interpolate (or extrapolate) from lines the total PI cross section for any given omega and initial level. The procedure applies a
        linear interpolation/extrapolation by just using the cross sections from the two nearest (given) omega points; a cross section
        cs::EmProperty is returned.
"""
function  interpolateCrossSection(lines::Array{PhotoIonization.Line,1}, omega::Float64, initialLevel)
    # First determine for which omegas cross sections are available and which associated cross sections are to be applied in the
    # interpolation
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
    
    return( cs )
end

end # module
