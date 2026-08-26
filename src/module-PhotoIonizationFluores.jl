
"""
`module  JAC.PhotoIonizationFluores`
... a submodule of JAC for photo-ionization followed by fluorescence,

        gamma(omega, P1 P2 P3)  +  A(J_0)   -->   A^+(J_m)  +  e_p(eps, kappa)   -->   A^+(J_f)  +  e_p  +  gamma'(omega')

    An atom is ionized by light of well-defined Stokes parameters; the residual ion is thereby left ALIGNED and, for
    circularly polarized light, ORIENTED, and its subsequent radiative decay is therefore neither isotropic nor
    unpolarized. The three observables are the STRENGTH of each fluorescence line, its ANGULAR DISTRIBUTION, and its
    POLARIZATION.

    The module works from well-defined initial, intermediate and final multiplets, in the manner of
    `DielectronicRecombination`, and is NOT a Cascade computation.

    ## The physics, in three stages

    **Stage A -- the incident light.**  A photon travelling along z with Stokes parameters (P1, P2, P3) has, in the
    helicity basis lambda = +-1, the density matrix  rho(lambda,lambda') = 1/2 [1+P3, -P1+iP2; -P1-iP2, 1-P3].
    Its statistical tensors follow by contracting that with <1 lambda, 1 -lambda'| k q>, and only k = 0, 1, 2 survive.
    The structure is what makes such an experiment interpretable: k = 0 carries the intensity, **k = 1 is non-zero only
    for CIRCULAR polarization and is the sole source of orientation**, and **k = 2 is non-zero even for UNPOLARIZED
    light** -- an unpolarized beam still defines an axis, hence alignment, but no handedness.

    **Stage B -- the residual ion.**  Coupling those photon tensors to the photo-ionization amplitudes and summing over
    the UNOBSERVED photoelectron leaves the ion with tensors rho_kq(J_m).  THIS ALREADY EXISTS IN JAC:
    `PhotoIonization.computeStatisticalTensorUnpolarized(k, q, line, settings)` returns rho_kq for both gauges, and
    despite its name it reads `settings.stokes` and builds the helicity matrix of stage A inline, so arbitrary
    (P1, P2, P3) are handled.  Its "Unpolarized" refers to the unobserved photoelectron, not to the incident light.
    A level with J_m = 0 or 1/2 cannot be aligned at all, so its fluorescence must come out isotropic and unpolarized
    whatever the incident light -- the null test that `examples/example-Dp.jl` runs first.

    **Stage C -- the fluorescence.**  The photon density matrix of the emitted line at (theta, phi) is the
    rho_kq-weighted sum over pairs of fluorescence multipoles, with Wigner D-functions, a Clebsch-Gordan coefficient
    and a 6j symbol; `computePhotonDm` below evaluates it.  Its intensity and Stokes parameters are then read off the
    2x2 helicity matrix directly, which avoids having to fix a separate convention for an intrinsic anisotropy
    coefficient alpha_2.

    ## How it is put together

    The two steps are computed INDEPENDENTLY of any pathway and only then combined:
    `PhotoIonization.computeLines` for i --> m + e_p, and `PhotoEmission.computeLines` for m --> f + gamma'.  Each
    (initial, intermediate) pair therefore carries its partial-wave sum once, however many final levels follow it.
    The strength of a pathway is the ionization cross section times the radiative branching ratio,

        S(i,m,f)  =  sigma_ion(i --> m)  *  A_r(m --> f) / Gamma_r(m) ,

    the exact analogue of the resonance strength of `DielectronicRecombination`, and both gauges are carried
    throughout as `EmProperty`.

    STATUS, 19-Aug-2026: first implementation.  Until 18-Aug-2026 this module computed nothing at all: it printed a
    "computation starts now" banner and returned the STRING "Not yet implemented !" as its list of pathways.
"""
module PhotoIonizationFluores

using Printf, ..AngularMomentum, ..Basics, ..Defaults, ..ManyElectron, ..Nuclear, ..PhotoEmission, ..PhotoIonization,
      ..Radial, ..TableStrings


"""
`struct  PhotoIonizationFluores.Settings`
    ... defines a type for the details and parameters of computing photo-ionization-fluorescence pathways
        |i(N)>  -->  |m(N-1)> + e_p  -->  |f(N-1)> + e_p + gamma'.

    + multipoles              ::Array{EmMultipole,1}   ... Multipoles of the IONIZING radiation field.
    + gauges                  ::Array{UseGauge,1}      ... Gauges to be included into the computations.
    + photonEnergies          ::Array{Float64,1}       ... List of ionizing photon energies [user-specified units].
    + lValues                 ::Array{Int64,1}         ... Orbital angular momenta of the photoelectron partial waves.
    + incidentStokes          ::ExpStokes              ... Stokes parameters (P1, P2, P3) of the incident radiation.
    + calcAngular             ::Bool                   ... Calculate the angular distribution of the fluorescence photon.
    + calcStokes              ::Bool                   ... Calculate the Stokes parameters of the fluorescence photon.
    + solidAngles             ::Array{SolidAngle,1}    ... Solid angles [(theta_1, phi_1), ...] at which they are wanted.
    + printBefore             ::Bool                   ... True, if all pathways are printed before their evaluation.
    + pathwaySelection        ::PathwaySelection       ... Specifies the selected levels/pathways, if any.
"""
struct Settings  <:  AbstractProcessSettings
    multipoles                ::Array{EmMultipole,1}
    gauges                    ::Array{UseGauge,1}
    photonEnergies            ::Array{Float64,1}
    lValues                   ::Array{Int64,1}
    incidentStokes            ::ExpStokes
    calcAngular               ::Bool
    calcStokes                ::Bool
    solidAngles               ::Array{SolidAngle,1}
    printBefore               ::Bool
    pathwaySelection          ::PathwaySelection
end


"""
`PhotoIonizationFluores.Settings()`  ... constructor for the default values of photo-ionization-fluorescence settings.
"""
function Settings()
    Settings( Basics.EmMultipole[E1], Basics.UseGauge[Basics.UseCoulomb, Basics.UseBabushkin], Float64[], Int64[0, 1, 2],
              Basics.ExpStokes(), true, true, SolidAngle[], false, PathwaySelection() )
end


"""
`PhotoIonizationFluores.Settings(set::PhotoIonizationFluores.Settings;`

        multipoles=..,        gauges=..,          photonEnergies=..,    lValues=..,        incidentStokes=..,
        calcAngular=..,       calcStokes=..,      solidAngles=..,       printBefore=..,    pathwaySelection=.. )

    ... constructor for modifying the given PhotoIonizationFluores.Settings by 'overwriting' the previously selected
        parameters. A settings::PhotoIonizationFluores.Settings is returned.
"""
function Settings(set::PhotoIonizationFluores.Settings;
    multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,      gauges::Union{Nothing,Array{UseGauge}}=nothing,
    photonEnergies::Union{Nothing,Array{Float64,1}}=nothing,      lValues::Union{Nothing,Array{Int64,1}}=nothing,
    incidentStokes::Union{Nothing,ExpStokes}=nothing,             calcAngular::Union{Nothing,Bool}=nothing,
    calcStokes::Union{Nothing,Bool}=nothing,                      solidAngles::Union{Nothing,Array{SolidAngle,1}}=nothing,
    printBefore::Union{Nothing,Bool}=nothing,                     pathwaySelection::Union{Nothing,PathwaySelection}=nothing)

    if  isnothing(multipoles)        multipolesx       = set.multipoles       else   multipolesx       = multipoles       end
    if  isnothing(gauges)            gaugesx           = set.gauges           else   gaugesx           = gauges           end
    if  isnothing(photonEnergies)    photonEnergiesx   = set.photonEnergies   else   photonEnergiesx   = photonEnergies   end
    if  isnothing(lValues)           lValuesx          = set.lValues          else   lValuesx          = lValues          end
    if  isnothing(incidentStokes)    incidentStokesx   = set.incidentStokes   else   incidentStokesx   = incidentStokes   end
    if  isnothing(calcAngular)       calcAngularx      = set.calcAngular      else   calcAngularx      = calcAngular      end
    if  isnothing(calcStokes)        calcStokesx       = set.calcStokes       else   calcStokesx       = calcStokes       end
    if  isnothing(solidAngles)       solidAnglesx      = set.solidAngles      else   solidAnglesx      = solidAngles      end
    if  isnothing(printBefore)       printBeforex      = set.printBefore      else   printBeforex      = printBefore      end
    if  isnothing(pathwaySelection)  pathwaySelectionx = set.pathwaySelection else   pathwaySelectionx = pathwaySelection end

    Settings( multipolesx, gaugesx, photonEnergiesx, lValuesx, incidentStokesx, calcAngularx, calcStokesx,
              solidAnglesx, printBeforex, pathwaySelectionx )
end


# `Base.show(io::IO, settings::PhotoIonizationFluores.Settings)`
#        ... prepares a proper printout of the variable settings::PhotoIonizationFluores.Settings.
function Base.show(io::IO, settings::PhotoIonizationFluores.Settings)
    println(io, "multipoles:              $(settings.multipoles)  ")
    println(io, "gauges:                  $(settings.gauges)  ")
    println(io, "photonEnergies:          $(settings.photonEnergies)  ")
    println(io, "lValues:                 $(settings.lValues)  ")
    println(io, "incidentStokes:          $(settings.incidentStokes)  ")
    println(io, "calcAngular:             $(settings.calcAngular)  ")
    println(io, "calcStokes:              $(settings.calcStokes)  ")
    println(io, "solidAngles:             $(settings.solidAngles)  ")
    println(io, "printBefore:             $(settings.printBefore)  ")
    println(io, "pathwaySelection:        $(settings.pathwaySelection)  ")
end


"""
`struct  PhotoIonizationFluores.Pathway`
    ... defines a type for one photo-ionization-fluorescence pathway i --> m --> f, assembled from one
        `PhotoIonization.Line` and one `PhotoEmission.Line`. No `Channel` type of its own is needed: the photoelectron
        partial waves live in the ionization line and the fluorescence multipole amplitudes in the emission line.

    + initialLevel        ::Level      ... initial-(state) level i of the neutral atom.
    + intermediateLevel   ::Level      ... intermediate-(state) level m of the residual ion.
    + finalLevel          ::Level      ... final-(state) level f of the residual ion.
    + photonEnergy        ::Float64    ... energy of the ionizing photon.
    + electronEnergy      ::Float64    ... energy of the emitted photoelectron.
    + fluorEnergy         ::Float64    ... energy of the fluorescence photon.
    + crossSection        ::EmProperty ... photo-ionization cross section sigma(i --> m).
    + photonRate          ::EmProperty ... radiative rate A_r(m --> f) of this line alone.
    + totalPhotonRate     ::EmProperty ... Gamma_r(m), summed over ALL final levels of the given multiplet.
    + strength            ::EmProperty ... sigma(i --> m) * A_r(m --> f) / Gamma_r(m).
    + alignmentA2         ::EmProperty ... alignment parameter A_20 = rho_20 / rho_00 of the intermediate level.
    + statisticalTensors  ::Dict{Tuple{Int64,Int64},EmPropertyC}  ... rho_kq of the intermediate level, keyed by (k,q).
"""
struct  Pathway
    initialLevel          ::Level
    intermediateLevel     ::Level
    finalLevel            ::Level
    photonEnergy          ::Float64
    electronEnergy        ::Float64
    fluorEnergy           ::Float64
    crossSection          ::EmProperty
    photonRate            ::EmProperty
    totalPhotonRate       ::EmProperty
    strength              ::EmProperty
    alignmentA2           ::EmProperty
    statisticalTensors    ::Dict{Tuple{Int64,Int64},EmPropertyC}
end


"""
`PhotoIonizationFluores.Pathway()`  ... constructor for an `empty` photo-ionization-fluorescence pathway.
"""
function Pathway()
    Pathway(Level(), Level(), Level(), 0., 0., 0., EmProperty(0.), EmProperty(0.), EmProperty(0.), EmProperty(0.),
            EmProperty(0.), Dict{Tuple{Int64,Int64},EmPropertyC}() )
end


# `Base.show(io::IO, pathway::PhotoIonizationFluores.Pathway)`
#        ... prepares a proper printout of the variable pathway::PhotoIonizationFluores.Pathway.
function Base.show(io::IO, pathway::PhotoIonizationFluores.Pathway)
    println(io, "initialLevel:            $(pathway.initialLevel)  ")
    println(io, "intermediateLevel:       $(pathway.intermediateLevel)  ")
    println(io, "finalLevel:              $(pathway.finalLevel)  ")
    println(io, "photonEnergy:            $(pathway.photonEnergy)  ")
    println(io, "electronEnergy:          $(pathway.electronEnergy)  ")
    println(io, "fluorEnergy:             $(pathway.fluorEnergy)  ")
    println(io, "crossSection:            $(pathway.crossSection)  ")
    println(io, "photonRate:              $(pathway.photonRate)  ")
    println(io, "totalPhotonRate:         $(pathway.totalPhotonRate)  ")
    println(io, "strength:                $(pathway.strength)  ")
    println(io, "alignmentA2:             $(pathway.alignmentA2)  ")
    println(io, "statisticalTensors:      $(pathway.statisticalTensors)  ")
end


"""
`PhotoIonizationFluores.computePhotonDm(pathway::PhotoIonizationFluores.Pathway, fluorLine::PhotoEmission.Line,
                                        solidAngle::SolidAngle, gauge::EmGauge)`
    ... computes the 2x2 helicity density matrix of the fluorescence photon emitted into the given solid angle, from
        the statistical tensors of the intermediate level and the multipole amplitudes of the fluorescence line. This
        is stage C of the module docstring; the expression is the one of `PhotoExcitationFluores.computePhotonDm`,
        here evaluated with the REAL rho_kq rather than with the literal 2.0 carried there.
        A `Dict{Tuple{Int64,Int64},ComplexF64}` keyed by (lambda, lambda') is returned.
"""
function computePhotonDm(pathway::PhotoIonizationFluores.Pathway, fluorLine::PhotoEmission.Line,
                         solidAngle::SolidAngle, gauge::EmGauge)
    Jm = pathway.intermediateLevel.J;    Jf = pathway.finalLevel.J
    dm = Dict{Tuple{Int64,Int64},ComplexF64}()
    for  lambda in [1, -1],  lambdap in [1, -1]
        wa = ComplexF64(0.)
        for  ((k, q), rho)  in  pathway.statisticalTensors
            rhoValue = (gauge == Basics.Babushkin) ? rho.Babushkin : rho.Coulomb
            qp = lambdap - lambda
            ## The photon helicities give qp = lambda' - lambda in {0, +-2}, which need not be carried by rank k:
            ## the Clebsch-Gordan <L lambda, L' -lambda' | k -qp> vanishes unless |qp| <= k, and asking for it
            ## outside that range is a DomainError rather than a zero.
            if  abs(qp) > k    ||    abs(q) > k    continue    end
            for  ma in fluorLine.amplitudes
                L  = ma.multipole.L;    p  = ma.multipole.electric ? 1 : 0
                ampa = (gauge == Basics.Babushkin) ? ma.amplitude.Babushkin : ma.amplitude.Coulomb
                for  mp in fluorLine.amplitudes
                    Lp = mp.multipole.L;   pp = mp.multipole.electric ? 1 : 0
                    if  !AngularMomentum.isTriangle(AngularJ64(L), AngularJ64(Lp), AngularJ64(k))         continue    end
                    if  !AngularMomentum.isTriangle(Jm, Jm, AngularJ64(k))                                continue    end
                    ampp = (gauge == Basics.Babushkin) ? mp.amplitude.Babushkin : mp.amplitude.Coulomb
                    wa = wa + AngularMomentum.Wigner_DFunction(k, -q, qp, solidAngle.phi, solidAngle.theta, 0.0) * rhoValue *
                              im^(Lp + pp - L - p) * lambda^p * lambdap^pp * sqrt((2L+1)*(2Lp+1))                          *
                              (-1)^( Int(Basics.twice(Jf)/2 + Basics.twice(Jm)/2) + k + q + 1 )                            *
                              AngularMomentum.ClebschGordan(AngularJ64(L), AngularM64(lambda), AngularJ64(Lp),
                                                            AngularM64(-lambdap), AngularJ64(k), AngularM64(-qp))         *
                              AngularMomentum.Wigner_6j(AngularJ64(L), AngularJ64(Lp), AngularJ64(k), Jm, Jm, Jf)         *
                              ampa * conj(ampp)
                end
            end
        end
        dm[(lambda, lambdap)] = wa
    end

    return( dm )
end


"""
`PhotoIonizationFluores.computeObservables(pathway::PhotoIonizationFluores.Pathway, fluorLine::PhotoEmission.Line,
                                           solidAngle::SolidAngle, gauge::EmGauge)`
    ... reads the intensity and the Stokes parameters of the fluorescence line off its helicity density matrix at the
        given solid angle. A tuple `(intensity::Float64, P1::Float64, P2::Float64, P3::Float64)` is returned; the
        Stokes parameters are normalized to the intensity and are zero whenever the ion carries no alignment.
"""
function computeObservables(pathway::PhotoIonizationFluores.Pathway, fluorLine::PhotoEmission.Line,
                            solidAngle::SolidAngle, gauge::EmGauge)
    dm  = PhotoIonizationFluores.computePhotonDm(pathway, fluorLine, solidAngle, gauge)
    tr  = real( dm[(1,1)] + dm[(-1,-1)] )
    if  abs(tr) < 1.0e-30    return( (0., 0., 0., 0.) )    end
    P1  = -2.0 * real( dm[(1,-1)] ) / tr
    P2  = -2.0 * imag( dm[(1,-1)] ) / tr
    P3  = real( dm[(1,1)] - dm[(-1,-1)] ) / tr

    return( (tr, P1, P2, P3) )
end


"""
`PhotoIonizationFluores.computePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet,
                                        nm::Nuclear.Model, grid::Radial.Grid, settings::PhotoIonizationFluores.Settings; output=true)`
    ... computes the photo-ionization-fluorescence pathways for the given multiplets and settings. The two steps are
        evaluated INDEPENDENTLY -- `PhotoIonization.computeLines` for i --> m + e_p and `PhotoEmission.computeLines`
        for m --> f + gamma' -- and only then combined into pathways, so that the expensive partial-wave sum is carried
        once per (i, m) pair however many final levels follow it. `settings.multipoles` is applied to BOTH steps.
        A list of pathways::Array{PhotoIonizationFluores.Pathway,1} is returned.
"""
function  computePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet,
                          nm::Nuclear.Model, grid::Radial.Grid, settings::PhotoIonizationFluores.Settings; output=true)
    println("")
    printstyled("PhotoIonizationFluores.computePathways(): The computation of photo-ionization-fluorescence pathways starts now ... \n",
                color=:light_green)
    printstyled("------------------------------------------------------------------------------------------------------------------ \n",
                color=:light_green)
    println("")
    if  length(settings.photonEnergies) == 0
        error("PhotoIonizationFluores: settings.photonEnergies is empty; give at least one ionizing photon energy.")
    end
    #
    ## STEP 1, independent of any pathway:  i --> m + e_p, with the statistical tensors switched on.
    ## calcTensors is left FALSE on purpose. It controls only PhotoIonization's own DISPLAY of the tensors, and that
    ## path is broken: `displayLines` calls computeStatisticalTensorUnpolarized(k, q, gauge, line, settings) with five
    ## arguments (module-PhotoIonization.jl:1404) while only the four-argument method exists, so `calcTensors = true`
    ## raises a MethodError and has evidently never been exercised. The tensors themselves are fine and are obtained
    ## here by calling that four-argument method directly, below.
    piSettings = PhotoIonization.Settings(PhotoIonization.Settings(); multipoles=settings.multipoles, gauges=settings.gauges,
                                          photonEnergies=settings.photonEnergies, lValues=settings.lValues,
                                          calcTensors=false, stokes=settings.incidentStokes, printBefore=settings.printBefore)
    piLines    = PhotoIonization.computeLines(intermediateMultiplet, initialMultiplet, nm, grid, piSettings; output=true)
    #
    ## STEP 2, likewise independent:  m --> f + gamma'.
    peSettings = PhotoEmission.Settings(PhotoEmission.Settings(); multipoles=settings.multipoles, gauges=settings.gauges,
                                        printBefore=settings.printBefore)
    peLines    = PhotoEmission.computeLines(finalMultiplet, intermediateMultiplet, grid, peSettings; output=true)
    #
    ## Gamma_r(m), the branching denominator, summed over ALL final levels of the given multiplet.
    totalRates = Dict{Int64,EmProperty}()
    for  peLine in peLines
        m  = peLine.initialLevel.index
        totalRates[m] = get(totalRates, m, EmProperty(0.)) + peLine.photonRate
    end
    #
    pathways = PhotoIonizationFluores.Pathway[]
    for  piLine in piLines
        mIndex = piLine.finalLevel.index
        ## The statistical tensors of the intermediate level; k = 0, 1, 2 is all a dipole-ionized ion can carry.
        tensors = Dict{Tuple{Int64,Int64},EmPropertyC}()
        for  k = 0:2,  q = -k:k
            rho = PhotoIonization.computeStatisticalTensorUnpolarized(k, q, piLine, piSettings)
            if  abs(rho.Coulomb) > 1.0e-15  ||  abs(rho.Babushkin) > 1.0e-15    tensors[(k,q)] = rho    end
        end
        rho00 = get(tensors, (0,0), EmPropertyC(0., 0.));    rho20 = get(tensors, (2,0), EmPropertyC(0., 0.))
        aC    = abs(rho00.Coulomb)   > 1.0e-30 ? real(rho20.Coulomb   / rho00.Coulomb)   : 0.
        aB    = abs(rho00.Babushkin) > 1.0e-30 ? real(rho20.Babushkin / rho00.Babushkin) : 0.
        alignment = EmProperty(aC, aB)
        #
        for  peLine in peLines
            if  peLine.initialLevel.index != mIndex    continue    end
            if  !Basics.selectLevelTriple(piLine.initialLevel, piLine.finalLevel, peLine.finalLevel, settings.pathwaySelection)
                continue
            end
            gamma_r  = totalRates[mIndex]
            strength = EmProperty( abs(gamma_r.Coulomb)   > 1.0e-30 ?
                                       piLine.crossSection.Coulomb   * peLine.photonRate.Coulomb   / gamma_r.Coulomb   : 0.,
                                   abs(gamma_r.Babushkin) > 1.0e-30 ?
                                       piLine.crossSection.Babushkin * peLine.photonRate.Babushkin / gamma_r.Babushkin : 0. )
            push!( pathways, PhotoIonizationFluores.Pathway( piLine.initialLevel, piLine.finalLevel, peLine.finalLevel,
                                                             piLine.photonEnergy, piLine.electronEnergy, peLine.omega,
                                                             piLine.crossSection, peLine.photonRate, gamma_r, strength,
                                                             alignment, tensors ) )
        end
    end
    #
    PhotoIonizationFluores.displayResults(stdout, pathways, peLines, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    PhotoIonizationFluores.displayResults(iostream, pathways, peLines, settings)    end

    if    output    return( pathways )
    else            return( nothing )
    end
end


"""
`PhotoIonizationFluores.displayResults(stream::IO, pathways::Array{PhotoIonizationFluores.Pathway,1},
                                       fluorLines::Array{PhotoEmission.Line,1}, settings::PhotoIonizationFluores.Settings)`
    ... displays the strengths and, if requested, the angular distribution and the Stokes parameters of the fluorescence
        lines. Nothing is returned.
"""
function  displayResults(stream::IO, pathways::Array{PhotoIonizationFluores.Pathway,1},
                         fluorLines::Array{PhotoEmission.Line,1}, settings::PhotoIonizationFluores.Settings)
    nx = 137
    println(stream, " ")
    println(stream, "  Photo-ionization with subsequent fluorescence: strengths and alignment")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, "     i -- m -- f        J^P(i) -- J^P(m) -- J^P(f)      omega        e_p          omega'      " *
                    "  sigma(i-m)      A_r(m-f)        Strength        A_20  ")
    println(stream, "                                                      " * TableStrings.inUnits("energy") *
                    "                                     [barn]          [1/s]           [barn]              ")
    println(stream, "  ", TableStrings.hLine(nx))
    for  pw in pathways
        si = LevelSymmetry(pw.initialLevel.J, pw.initialLevel.parity)
        sm = LevelSymmetry(pw.intermediateLevel.J, pw.intermediateLevel.parity)
        sf = LevelSymmetry(pw.finalLevel.J, pw.finalLevel.parity)
        sa = "  " * @sprintf("%3d -- %3d -- %3d", pw.initialLevel.index, pw.intermediateLevel.index, pw.finalLevel.index)
        sa = sa * "     " * @sprintf("%-8s %-8s %-8s", string(si), string(sm), string(sf))
        sa = sa * @sprintf("%.4e  %.4e  %.4e  ", Defaults.convertUnits("energy: from atomic", pw.photonEnergy),
                           Defaults.convertUnits("energy: from atomic", pw.electronEnergy),
                           Defaults.convertUnits("energy: from atomic", pw.fluorEnergy))
        sa = sa * @sprintf("%.4e  %.4e  %.4e  %8.5f", pw.crossSection.Coulomb, pw.photonRate.Coulomb,
                           pw.strength.Coulomb, pw.alignmentA2.Coulomb)
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    if  (settings.calcAngular || settings.calcStokes)  &&  length(settings.solidAngles) > 0
        println(stream, " ")
        println(stream, "  Angular distribution and Stokes parameters of the fluorescence photon (Coulomb gauge)")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        println(stream, "     i -- m -- f       theta       phi          W(theta,phi)        P1            P2            P3   ")
        println(stream, "  ", TableStrings.hLine(nx))
        for  pw in pathways
            fluorLine = nothing
            for  fl in fluorLines
                if  fl.initialLevel.index == pw.intermediateLevel.index  &&  fl.finalLevel.index == pw.finalLevel.index
                    fluorLine = fl;   break
                end
            end
            isnothing(fluorLine)  &&  continue
            for  sAngle in settings.solidAngles
                (w, p1, p2, p3) = PhotoIonizationFluores.computeObservables(pw, fluorLine, sAngle, Basics.Coulomb)
                sa = "  " * @sprintf("%3d -- %3d -- %3d", pw.initialLevel.index, pw.intermediateLevel.index, pw.finalLevel.index)
                sa = sa * @sprintf("   %8.4f   %8.4f     %.6e   %11.6f   %11.6f   %11.6f", sAngle.theta, sAngle.phi, w, p1, p2, p3)
                println(stream, sa)
            end
        end
        println(stream, "  ", TableStrings.hLine(nx))
    end

    return( nothing )
end

end # module
