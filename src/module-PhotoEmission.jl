
"""
`module  JAC.PhotoEmission`  
... a submodel of JAC that contains all methods for computing Einstein coefficients, oscillator strength, etc. between 
    some initial and final-state multiplets.    ***done***
"""
module PhotoEmission


using Printf, ..AngularMomentum, ..Basics, ..BiOrthogonal, ..Defaults, ..InteractionStrength, ..ManyElectron, ..Radial,
                ..SpinAngular, ..TableStrings


"""
`struct  PhotoEmission.Settings  <:  AbstractProcessSettings`  
    ... defines a type for the details and parameters of computing radiative lines.

    + multipoles              ::Array{EmMultipoles}     ... Specifies the (radiat. field) multipoles to be included.
    + gauges                  ::Array{UseGauge}         ... Gauges to be included into the computations.
    + calcAnisotropy          ::Bool                    ... True, if the anisotropy (structure) functions are to be 
                                                            calculated and false otherwise 
    + printBefore             ::Bool                    ... True, if all energies and lines are printed before comput.
    + corePolarization        ::CorePolarization        ... Parametrization of the core-polarization potential/contribution.
    + lineSelection           ::LineSelection           ... Specifies the selected levels, if any.
    + photonEnergyShift       ::Float64                 ... An overall energy shift for all photon energies.
    + mimimumPhotonEnergy     ::Float64                 ... minimum transition energy for which (photon) transitions
                                                            are included into the computation.
    + maximumPhotonEnergy     ::Float64                 ... maximum transition energy for which (photon) transitions
                                                            are included.
    + calcBiorthogonal        ::Bool                    ... True, if the initial- and final-state multiplets are first
                                                            brought into a bi-orthogonal representation
                                                            (`BiOrthogonal.computeTransformation`) before the transition
                                                            amplitudes are evaluated, and false (the default) if the
                                                            two multiplets are used as they are.
"""
struct Settings  <:  AbstractProcessSettings
    multipoles                ::Array{EmMultipole,1}
    gauges                    ::Array{UseGauge}
    calcAnisotropy            ::Bool
    printBefore               ::Bool
    corePolarization          ::CorePolarization
    lineSelection             ::LineSelection
    photonEnergyShift         ::Float64
    mimimumPhotonEnergy       ::Float64
    maximumPhotonEnergy       ::Float64
    calcBiorthogonal          ::Bool
end


"""
`PhotoEmission.Settings()`  ... constructor for the default values of radiative line computations
"""
function Settings()
    Settings(EmMultipole[E1], UseGauge[Basics.UseCoulomb], false, false, CorePolarization(), LineSelection(), 0., 0., 10000., false)
end


"""
`PhotoEmission.Settings(set::PhotoEmission.Settings;`

        multipoles::=..,        gauges=..,                calcAnisotropy=..,          printBefore=..,
        corePolarization=..,    lineSelection=..,         photonEnergyShift=..,
        mimimumPhotonEnergy=.., maximumPhotonEnergy=..,   calcBiorthogonal=..)

    ... constructor for modifying the given PhotoEmission.Settings by 'overwriting' the previously selected parameters.
"""
function Settings(set::PhotoEmission.Settings;
    multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,    gauges::Union{Nothing,Array{UseGauge}}=nothing,
    calcAnisotropy::Union{Nothing,Bool}=nothing,                printBefore::Union{Nothing,Bool}=nothing,
    corePolarization::Union{Nothing,CorePolarization}=nothing,  lineSelection::Union{Nothing,LineSelection}=nothing,
    photonEnergyShift::Union{Nothing,Float64}=nothing,          mimimumPhotonEnergy::Union{Nothing,Float64}=nothing,
    maximumPhotonEnergy::Union{Nothing,Float64}=nothing,        calcBiorthogonal::Union{Nothing,Bool}=nothing)

    if  isnothing(multipoles)            multipolesx          = set.multipoles              else  multipolesx          = multipoles            end
    if  isnothing(gauges)                gaugesx              = set.gauges                  else  gaugesx              = gauges                end
    if  isnothing(calcAnisotropy)        calcAnisotropyx      = set.calcAnisotropy          else  calcAnisotropyx      = calcAnisotropy        end
    if  isnothing(printBefore)           printBeforex         = set.printBefore             else  printBeforex         = printBefore           end
    if  isnothing(corePolarization)      corePolarizationx    = set.corePolarization        else  corePolarizationx    = corePolarization      end
    if  isnothing(lineSelection)         lineSelectionx       = set.lineSelection           else  lineSelectionx       = lineSelection         end
    if  isnothing(photonEnergyShift)     photonEnergyShiftx   = set.photonEnergyShift       else  photonEnergyShiftx   = photonEnergyShift     end
    if  isnothing(mimimumPhotonEnergy)   mimimumPhotonEnergyx = set.mimimumPhotonEnergy     else  mimimumPhotonEnergyx = mimimumPhotonEnergy   end
    if  isnothing(maximumPhotonEnergy)   maximumPhotonEnergyx = set.maximumPhotonEnergy     else  maximumPhotonEnergyx = maximumPhotonEnergy   end
    if  isnothing(calcBiorthogonal)      calcBiorthogonalx    = set.calcBiorthogonal        else  calcBiorthogonalx    = calcBiorthogonal      end

    Settings( multipolesx, gaugesx, calcAnisotropyx, printBeforex, corePolarizationx, lineSelectionx,
              photonEnergyShiftx, mimimumPhotonEnergyx, maximumPhotonEnergyx, calcBiorthogonalx)
end


# `Base.show(io::IO, settings::PhotoEmission.Settings)`  ... prepares a proper printout of the variable settings::PhotoEmissionSettings.
function Base.show(io::IO, settings::PhotoEmission.Settings)
    println(io, "multipoles:             $(settings.multipoles)  ")
    println(io, "gauges:                 $(settings.gauges)  ")
    println(io, "calcAnisotropy:         $(settings.calcAnisotropy)  ")
    println(io, "printBefore:            $(settings.printBefore)  ")
    println(io, "corePolarization:       $(settings.corePolarization)  ")
    println(io, "lineSelection:          $(settings.lineSelection)  ")
    println(io, "photonEnergyShift:      $(settings.photonEnergyShift)  ")
    println(io, "mimimumPhotonEnergy:    $(settings.mimimumPhotonEnergy)  ")
    println(io, "maximumPhotonEnergy:    $(settings.maximumPhotonEnergy)  ")
    println(io, "calcBiorthogonal:       $(settings.calcBiorthogonal)  ")
end










#####################################################################################################################
## THE PHYSICAL FORM, built BESIDE the flat one above -- the third module, after PhotoIonization and
## PhotoRecombination, and the one where the CHANNEL CONCEPT DISAPPEARS ALTOGETHER.
##
## Photo emission is BOUND-BOUND: there is no free electron, hence no kappa, no scattering phase and no total
## symmetry of a scattering state.  A PhotoEmission.Channel is (multipole, gauge, amplitude) -- that is,
## ENTIRELY a label of the interaction OPERATOR, with no state content whatsoever.  In the other two modules the
## channel had to be reorganized, separating the state (kappa, symmetry) from the operator (multipole, gauge);
## here there is nothing to separate, so the whole intermediate layer goes away and a line carries a flat list of
## multipole amplitudes.  No PartialWave, no Channel: a struct holding exactly one field would describe nothing.
##
## For a line with [E1, M1, E2] and both gauges requested:
##     flat  ->  5 Channel objects:   E1-Coulomb, E1-Babushkin, M1-Magnetic, E2-Coulomb, E2-Babushkin
##     new   ->  3 amplitude entries: E1, M1, E2 -- each carrying BOTH gauges
#####################################################################################################################




"""
`struct  PhotoEmission.Line`
    ... as PhotoEmission.Line, but carrying one entry per MULTIPOLE instead of one per (multipole, gauge).

    + initialLevel   ::Level                            ... initial-(state) level
    + finalLevel     ::Level                            ... final-(state) level
    + omega          ::Float64                          ... Transition frequency of this line.
    + photonRate     ::EmProperty                       ... Total rate of this line.
    + angularBeta    ::EmProperty                       ... Angular beta_2 coefficient.
    + amplitudes     ::Array{MultipoleAmplitude,1} ... one entry per contributing multipole.

        `hasSublines` is deliberately NOT carried over. In the flat Line it is declared, stored and printed and
        never read -- not here, not in PhotoExcitation, which is the only other survivor of the field, and not
        in any of the sixteen modules that use PhotoEmission.Line.
"""
struct  Line
    initialLevel     ::Level
    finalLevel       ::Level
    omega            ::Float64
    photonRate       ::EmProperty
    angularBeta      ::EmProperty
    amplitudes       ::Array{MultipoleAmplitude,1}
end




"""
`PhotoEmission.amplitude(::Emission, Mp::EmMultipole, gauge::EmGauge, omega::Float64, finalLevel::Level, initialLevel::Level,
                            grid::Radial.Grid; display::Bool=false, printout::Bool=false)`
    ... to compute the photon emission amplitude  <alpha_f J_f || O^(Mp) || alpha_i J_i> for the
        interaction with a photon of multipolarity Mp and for the given transition energy and gauge. A value::ComplexF64 is
        returned. The amplitude value is printed to screen if display=true.
"""
function amplitude(::Emission, Mp::EmMultipole, gauge::EmGauge, omega::Float64, finalLevel::Level, initialLevel::Level,
                    grid::Radial.Grid; display::Bool=false, printout::Bool=false)
    if  initialLevel.basis.subshells == finalLevel.basis.subshells
        iLevel = initialLevel;   fLevel = finalLevel
    else
        subshells = Basics.merge(initialLevel.basis.subshells, finalLevel.basis.subshells)
        iLevel    = Level(initialLevel, subshells)
        fLevel    = Level(finalLevel, subshells)
    end

    nf = length(fLevel.basis.csfs);    ni = length(iLevel.basis.csfs)
    if  printout   printstyled("Compute radiative $(Mp) matrix of dimension $nf x $ni in the initial- and final-state bases " *
                                "for the transition [$(iLevel.index)-$(fLevel.index)] ... ", color=:light_green)    end
    matrix = zeros(ComplexF64, nf, ni)
    #
    for  r = 1:nf
        if  fLevel.basis.csfs[r].J != fLevel.J      ||  fLevel.basis.csfs[r].parity  != fLevel.parity    continue    end
        for  s = 1:ni
            if  iLevel.basis.csfs[s].J != iLevel.J  ||  iLevel.basis.csfs[s].parity  != iLevel.parity    continue    end
            subshellList = fLevel.basis.subshells
            opa = SpinAngular.OneParticleOperator(Mp.L, plus, true)
            wa  = SpinAngular.computeCoefficients(opa, fLevel.basis.csfs[r], iLevel.basis.csfs[s], subshellList)
            me = 0.
            for  coeff in wa
                MabEm = InteractionStrength.MabEmission(Mp, gauge, omega, fLevel.basis.orbitals[coeff.a],
                                                                                        iLevel.basis.orbitals[coeff.b], grid)
                ja = Basics.subshell_2j(fLevel.basis.orbitals[coeff.a].subshell)
                ## jb = Basics.subshell_2j(iLevel.basis.orbitals[coeff.b].subshell)
                me = me + coeff.T * MabEm / sqrt( ja + 1) * sqrt( (Basics.twice(fLevel.J) + 1))      ## * sqrt( jb + 1)
            end
            matrix[r,s] = me
        end
    end
    if  printout   printstyled("done. \n", color=:light_green)    end
    amplitude = transpose(fLevel.mc) * matrix * iLevel.mc
    # Multiply with the multipolarity factors to keep different multipoles on the same footings; this factor need to be better understood
    # amplitude = amplitude * sqrt( (2Mp.L+1)*(Mp.L+1)/Mp.L )

    if  display
        println("    < level=$(finalLevel.index) [J=$(finalLevel.J)$(string(finalLevel.parity))] ||" *
                " O^($Mp, emission) ($omega a.u., $gauge) ||" *
                " $(initialLevel.index) [$(initialLevel.J)$(string(initialLevel.parity))] >  = $amplitude  ")
    end

    return( amplitude )
end


"""
` + amplitude(::Absorption, Mp::EmMultipole, gauge::EmGauge, omega::Float64, finalLevel::Level, initialLevel::Level,
                grid::Radial.Grid; display::Bool=false, printout::Bool=false)`
    ... to compute the photon absorption amplitude as the conjugate of the emission amplitude with swapped levels,
        i.e. conj( <initialLevel || O^(Mp) || finalLevel> ).  A value::ComplexF64 is returned.
"""
function amplitude(::Absorption, Mp::EmMultipole, gauge::EmGauge, omega::Float64, finalLevel::Level, initialLevel::Level,
                    grid::Radial.Grid; display::Bool=false, printout::Bool=false)
    amplitude = conj( PhotoEmission.amplitude(Emission(), Mp, gauge, omega, initialLevel, finalLevel, grid; printout=printout) )
    if  display
        println("    < level=$(finalLevel.index) [J=$(finalLevel.J)$(string(finalLevel.parity))] ||" *
                " O^($Mp, absorption) ($omega a.u., $gauge) ||" *
                " $(initialLevel.index) [$(initialLevel.J)$(string(initialLevel.parity))] >  = $amplitude  ")
    end
    return( amplitude )
end






















#####################################################################################################################
## STAGE A of the physical form: the core.  Everything below is ADDITIVE -- nothing above is altered.
#####################################################################################################################


"""
`PhotoEmission.determineChannels(finalLevel::Level, initialLevel::Level, settings::PhotoEmission.Settings)`
    ... as PhotoEmission.determineChannels, but returning one entry per MULTIPOLE; an
        Array{MultipoleAmplitude,1} is returned with all amplitudes still zero.

        The selection rule is the same AngularMomentum.isAllowedMultipole. What disappears is everything that
        was only there to service the gauge: the `for gauge in settings.gauges` loop, the three-way push, and
        the `hasMagnetic` flag that stopped a magnetic multipole being emitted once per requested gauge. The
        flag was necessary in the flat form and is unrepresentable here -- there is no gauge loop to guard.
"""
function determineChannels(finalLevel::Level, initialLevel::Level, settings::PhotoEmission.Settings)
    amplitudes = MultipoleAmplitude[]
    symi = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity)
    for  mp in settings.multipoles
        if   AngularMomentum.isAllowedMultipole(symi, mp, symf)
            push!(amplitudes, MultipoleAmplitude(mp, EmPropertyC(Complex(0.), Complex(0.))))
        end
    end
    return( amplitudes )
end


"""
`PhotoEmission.computeAmplitudesProperties(line::PhotoEmission.Line, grid::Radial.Grid,
                                                 settings::PhotoEmission.Settings; printout::Bool=true)`
    ... as PhotoEmission.computeAmplitudesProperties, but on the physical form; a Line with all
        amplitudes and the photon rate evaluated is returned.

        An electric multipole is evaluated twice, once per gauge, into one EmPropertyC; a magnetic multipole
        once, into an EmPropertyC with equal components. The three-way gauge `if` that accumulated rateC and
        rateB by hand becomes a single `rate = rate + abs2(ma.amplitude)`, since abs2 of an EmPropertyC is an
        EmProperty and a magnetic amplitude enters both components by itself.

        The rate is summed INCOHERENTLY over multipoles, exactly as in the flat version.
"""
function computeAmplitudesProperties(line::PhotoEmission.Line, grid::Radial.Grid,
                                           settings::PhotoEmission.Settings; printout::Bool=true)
    if  settings.corePolarization.doApply
        ## RETIRED 09-Aug-2026, as in computeAmplitudesProperties above; refuse rather than return a number
        ## nobody can defend.
        error("\n\nPhotoEmission: corePolarization.doApply = true, but the core-polarization amplitude was "  *
              "RETIRED on 09-Aug-2026.\n"                                                                     *
              ">>> It was based on MbaEmissionMigdalek, which never worked in a useful form.\n"                *
              ">>> Set doApply = false, or implement a new core-polarization correction from scratch.\n")
    end
    newAmplitudes = MultipoleAmplitude[];    rate = EmProperty(0., 0.)
    for  ma in line.amplitudes
        mp = ma.multipole
        if  string(mp)[1] == 'E'
            ampC = PhotoEmission.amplitude(Emission(), mp, Basics.Coulomb,   line.omega, line.finalLevel,
                                            line.initialLevel, grid; printout=printout)
            ampB = PhotoEmission.amplitude(Emission(), mp, Basics.Babushkin, line.omega, line.finalLevel,
                                            line.initialLevel, grid; printout=printout)
            amp  = EmPropertyC(ampC, ampB)
        else
            ampM = PhotoEmission.amplitude(Emission(), mp, Basics.Magnetic,  line.omega, line.finalLevel,
                                            line.initialLevel, grid; printout=printout)
            amp  = EmPropertyC(ampM)
        end
        rate = rate + abs2(amp)
        push!(newAmplitudes, MultipoleAmplitude(mp, amp))
    end
    wa          = 8pi * Defaults.getDefaults("alpha") * line.omega / (Basics.twice(line.initialLevel.J) + 1)
    photonrate  = wa * rate
    angularBeta = EmProperty(-9., -9.)      ## the sentinel of the flat version; never computed there either
    return( PhotoEmission.Line(line.initialLevel, line.finalLevel, line.omega, photonrate, angularBeta,
                                     newAmplitudes) )
end


"""
`PhotoEmission.computeAnisotropyFunctions(line::PhotoEmission.Line)`
    ... computes the anisotropy (structure) functions f_2 and f_4 of the given line; a tuple
        (f2, f4)::Tuple{EmPropertyC,EmPropertyC} is returned, each holding BOTH gauges.

        THIS IS THE FUNCTION THAT JUSTIFIES THE NEW FORM IN THIS MODULE. In displayAnisotropies the Coulomb
        and the Babushkin sums are written out as two blocks of twenty-two lines that differ in nothing but
        the gauge name -- the only place in the three converted modules where the gauge branching is literally
        duplicated code. Written once against EmPropertyC they are one block: the product
        conj(ampa) * ampb is componentwise, so Coulomb meets Coulomb and Babushkin meets Babushkin, and a
        magnetic amplitude has equal components and so enters both.

        BOTH defects fixed in 8bf17cb lived in those blocks and neither is expressible here: there is no gauge
        to test on either partner, and the pair is destructured into (mpa, ampa) and (mpb, ampb), so writing
        the outer multipole where the inner is meant is a visible mistake rather than an invisible one.

        The angular algebra, the phase factor and the parity factor are reproduced verbatim.
"""
function computeAnisotropyFunctions(line::PhotoEmission.Line)
    f2 = EmPropertyC(0.0im);   f4 = EmPropertyC(0.0im);   norm = EmProperty(0., 0.)
    angJi = line.initialLevel.J;    angJf = line.finalLevel.J
    for  ma in line.amplitudes    norm = norm + abs2(ma.amplitude)    end
    #
    for  ma in line.amplitudes
        mpa = ma.multipole;    angL  = AngularJ64(mpa.L);   p  = Basics.multipole_p(mpa)
        for  mb in line.amplitudes
            mpb = mb.multipole;    angLp = AngularJ64(mpb.L);   pp = Basics.multipole_p(mpb)
            #
            wa = (1.0im)^(mpb.L + pp - mpa.L - p)                                                    *
                 sqrt( (2mpa.L+1) * (2mpb.L+1) )                                                     *
                 (1. + (-1.)^(mpa.L + p + mpb.L + pp - 2) )                                          *
                 (conj(ma.amplitude) * mb.amplitude)
            f2 = f2 + wa * AngularMomentum.phaseFactor([angJf, +1, angJi, +1, AngularJ64(3)])        *
                 AngularMomentum.ClebschGordan(angL, AngularM64(1), angLp, AngularM64(-1), AngularJ64(2), AngularM64(0)) *
                 AngularMomentum.Wigner_6j(angL, angLp, AngularJ64(2), angJi, angJi, angJf)
            f4 = f4 + wa * AngularMomentum.phaseFactor([angJf, +1, angJi, +1, AngularJ64(5)])        *
                 AngularMomentum.ClebschGordan(angL, AngularM64(1), angLp, AngularM64(-1), AngularJ64(4), AngularM64(0)) *
                 AngularMomentum.Wigner_6j(angL, angLp, AngularJ64(4), angJi, angJi, angJf)
        end
    end
    wn = sqrt(Basics.twice(line.initialLevel.J) + 1) / 2.
    f2 = EmPropertyC(f2.Coulomb / norm.Coulomb * wn, f2.Babushkin / norm.Babushkin * wn)
    f4 = EmPropertyC(f4.Coulomb / norm.Coulomb * wn, f4.Babushkin / norm.Babushkin * wn)
    return( f2, f4 )
end


#####################################################################################################################
## STAGE B: the bridge back to the flat form, and the displays.
#####################################################################################################################






"""
`PhotoEmission.displayLines(stream::IO, lines::Array{PhotoEmission.Line,1})`
`PhotoEmission.displayRates(stream::IO, lines::Array{PhotoEmission.Line,1}, settings)`
`PhotoEmission.displayLifetimes(stream::IO, lines::Array{PhotoEmission.Line,1}, settings)`
`PhotoEmission.displayAnisotropies(stream::IO, lines::Array{PhotoEmission.Line,1}, settings)`
    ... the display layer.  These are the tables printed before the flat channel form was retired, reading the
        line's multipole amplitudes instead of its channels.

        displayLines prints one row per (multipole, GAUGE), which is no longer how an amplitude is stored, so
        that expansion is written where the table is built.  displayAnisotropies no longer computes anything:
        it calls computeAnisotropyFunctions, which is the same physics written once instead of twice.
"""
function  displayLines(stream::IO, lines::Array{PhotoEmission.Line,1})
    nx = 95
    println(stream, " ")
    println(stream, "  Selected radiative lines:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(14, "Energy"; na=4);              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.flushleft(30, "List of multipoles"; na=4);             sb = sb * TableStrings.hBlank(34)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", line.omega)) * "    "
        mpGaugeList = Tuple{EmMultipole,EmGauge}[]
        ## One row per (multipole, GAUGE).  An electric multipole is one amplitude holding two gauges, a
        ## magnetic one has no gauge freedom; the gauge is a property of this TABLE, not of the amplitude.
        for  ma in line.amplitudes
            if  string(ma.multipole)[1] == 'E'
                push!( mpGaugeList, (ma.multipole, Basics.Coulomb) )
                push!( mpGaugeList, (ma.multipole, Basics.Babushkin) )
            else
                push!( mpGaugeList, (ma.multipole, Basics.Magnetic) )
            end
        end
        sa = sa * TableStrings.multipoleGaugeTupels(50, mpGaugeList)
        println(stream,  sa )
    end
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, " ")
    #
    return( nothing )
end


function  displayRates(stream::IO, lines::Array{PhotoEmission.Line,1}, settings::PhotoEmission.Settings)
    nx = 161
    println(stream, " ")
    println(stream, "  Einstein coefficients, transition rates and oscillator strengths:")
    println(stream, " ")
    if  settings.corePolarization.doApply
        println(stream, "  + Calculate E1 amplitude with core-polarization corrections.")
        println(stream, "  + alpha_c = $(settings.corePolarization.coreAlpha) a.u.")
        println(stream, "  + r_c     = $(settings.corePolarization.coreRadius) a.u.")
        println(stream, " ")
    end
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(12, "Energy"   ; na=4);               
    sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center( 9, "Multipole"; na=0);                         sb = sb * TableStrings.hBlank(10)
    sa = sa * TableStrings.center(11, "Gauge"    ; na=4);                         sb = sb * TableStrings.hBlank(17)
    sa = sa * TableStrings.center(26, "A--Einstein--B"; na=3);       
    sb = sb * TableStrings.center(26, TableStrings.inUnits("rate")*"          "*TableStrings.inUnits("rate"); na=2)
    sa = sa * TableStrings.center(11, "Osc. strength"    ; na=3);                 sb = sb * TableStrings.hBlank(17)
    sa = sa * TableStrings.center(12, "Decay widths"; na=3);       
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(13, "Line strength"; na=4);       
    sb = sb * TableStrings.center(12, "[a.u.]"       ; na=4)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        ## This table has a GAUGE column, so one amplitude prints as two rows for an electric multipole and
        ## one for a magnetic one.  The gauge belongs to the presentation here; the amplitude carries both.
        for  ma in line.amplitudes
            if  string(ma.multipole)[1] == 'E'
                gaugeValues = [(Basics.Coulomb, ma.amplitude.Coulomb), (Basics.Babushkin, ma.amplitude.Babushkin)]
            else
                gaugeValues = [(Basics.Magnetic, ma.amplitude.Coulomb)]
            end
            for  (gauge, amplitude) in gaugeValues
                sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                               fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
                sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
                sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
                sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.omega)) * "    "
                sa = sa * TableStrings.center(9,  string(ma.multipole); na=4)
                sa = sa * TableStrings.flushleft(11, string(gauge);  na=2)
                chRate =  8pi * Defaults.getDefaults("alpha") * line.omega / (Basics.twice(line.initialLevel.J) + 1) * (abs(amplitude)^2)
                sa = sa * @sprintf("%.6e", Basics.recast(RecastRateToEinsteinA(),    line, chRate)) * "  "
                sa = sa * @sprintf("%.6e", Basics.recast(RecastRateToEinsteinB(),    line, chRate)) * "    "
                sa = sa * @sprintf("%.6e", Basics.recast(RecastRateToOscillatorGf(),           line, chRate)) * "    "
                sa = sa * @sprintf("%.6e", Basics.recast(RecastRateToDecayWidth(),   line, chRate)) * "    "
                if  ma.multipole == E1
                        sa = sa * @sprintf("%.6e", Basics.recast(RecastRateToLineStrengthS(),     line, chRate)) * "    "
                else    sa = sa * "  --  "
                end
                println(stream, sa)
            end
        end
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


function  displayLifetimes(stream::IO, lines::Array{PhotoEmission.Line,1}, settings::PhotoEmission.Settings)
    # Determine all initial levels (and their level information) for printing the lifetimes
    ilevels = Int64[];   istr = String[]
    for  i = 1:length(lines)
        ii = lines[i].initialLevel.index;   
        if  !(ii in ilevels)   
            sa  = "  ";    sym = LevelSymmetry( lines[i].initialLevel.J, lines[i].initialLevel.parity )
            sa = sa * TableStrings.center(10, TableStrings.level(lines[i].initialLevel.index); na=2)
            sa = sa * TableStrings.center(10, string(sym); na=4)
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", lines[i].initialLevel.energy)) * "    "
            push!( ilevels, ii);    push!( istr, sa )   
        end
    end
    # Determine the lifetime (in a.u.) of the selected initial levels
    irates = Basics.EmProperty[]
    for  ii in  ilevels
        waCoulomb = waBabushkin = 0.
        for  i = 1:length(lines)
            if   lines[i].initialLevel.index == ii    
                waCoulomb   = waCoulomb   + lines[i].photonRate.Coulomb
                waBabushkin = waBabushkin + lines[i].photonRate.Babushkin
            end
        end
        push!(irates, EmProperty( waCoulomb, waBabushkin) )
    end
    
    nx = 105
    println(stream, " ")
    println(stream, "  PhotoEmission lifetimes (as derived from these computations):")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                              sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                              sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(12, "Level energy"   ; na=3);               
    sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=3)
    sa = sa * TableStrings.center(12, "Used Gauge"    ; na=6);                     sb = sb * TableStrings.hBlank(18)
    sa = sa * TableStrings.center(26, "Lifetime"; na=3);       
    sb = sb * TableStrings.center(26, "[a.u.]"*"          "*TableStrings.inUnits("time"); na=5)
    sa = sa * TableStrings.center(12, "Decay widths"; na=4);       
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #    
    for  ii = 1:length(ilevels)
        sa = istr[ii]
        sa = sa * "Coulomb          " * @sprintf("%.6e",              1.0/irates[ii].Coulomb)     * "  "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("time: from atomic",   1.0/irates[ii].Coulomb) )   * "    "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic",     irates[ii].Coulomb) )
        println(stream, sa)
        sa = repeat(" ", length(istr[ii]) )
        sa = sa * "Babushkin        " * @sprintf("%.6e",              1.0/irates[ii].Babushkin)   * "  "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("time: from atomic",   1.0/irates[ii].Babushkin) ) * "    "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic",     irates[ii].Babushkin) )
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


function  displayAnisotropies(stream::IO, lines::Array{PhotoEmission.Line,1}, settings::PhotoEmission.Settings)
    nx = 153
    println(stream, " ")
    println(stream, "  Anisotropy (structure) functions:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(14, "Energy"   ; na=4);               
    sb = sb * TableStrings.center(14,TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.flushleft(20, "Multipoles";   na=0);              
    sa = sa * TableStrings.center(14, "f_2 (Coulomb)";   na=3);       
    sa = sa * TableStrings.center(14, "f_2 (Babushkin)"; na=3);       
    sa = sa * TableStrings.center(14, "f_4 (Coulomb)";   na=3);       
    sa = sa * TableStrings.center(14, "f_4 (Babushkin)"; na=3);       
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  line in lines
        ## THE COMPUTATION IS NO LONGER DONE HERE.  Until the flat channel form was retired this function
        ## carried it inline, as two blocks of twenty-two lines differing in nothing but the gauge name -- the
        ## only place in the five converted modules where the gauge branching was literally duplicated code,
        ## and where both defects fixed in 8bf17cb lived.  computeAnisotropyFunctions is that physics written
        ## once against EmPropertyC, and this table now simply asks for it.
        f2, f4 = PhotoEmission.computeAnisotropyFunctions(line)
        f2Coulomb = f2.Coulomb;   f2Babushkin = f2.Babushkin
        f4Coulomb = f4.Coulomb;   f4Babushkin = f4.Babushkin
        mpList    = EmMultipole[]
        for  ma in line.amplitudes    if  !(ma.multipole in mpList)    push!(mpList, ma.multipole)    end    end
        sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                       fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", line.omega)) * "    "
        #
        sa = sa * TableStrings.flushleft( 17, TableStrings.multipoleList(mpList);  na=3)
        sa = sa * TableStrings.flushright(15, @sprintf("%.8e", f2Coulomb.re);    na=3) 
        sa = sa * TableStrings.flushright(15, @sprintf("%.8e", f2Babushkin.re);  na=3)
        sa = sa * TableStrings.flushright(15, @sprintf("%.8e", f4Coulomb.re);    na=3)
        sa = sa * TableStrings.flushright(15, @sprintf("%.8e", f4Babushkin.re);  na=3)
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


#####################################################################################################################
## STAGE C: the drivers, so that a whole computation can be run either way and compared end to end.
#####################################################################################################################


"""
`PhotoEmission.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet,
                                    settings::PhotoEmission.Settings)`
    ... as PhotoEmission.determineLines, but producing Line; an Array{PhotoEmission.Line,1}
        is returned.
"""
function determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::PhotoEmission.Settings)
    lines = PhotoEmission.Line[]
    for  iLevel  in  initialMultiplet.levels
        for  fLevel  in  finalMultiplet.levels
            if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                omega = iLevel.energy - fLevel.energy   + settings.photonEnergyShift
                ## `mimimumPhotonEnergy` is spelled so in the Settings struct; kept verbatim.
                if  omega <= settings.mimimumPhotonEnergy  ||  omega > settings.maximumPhotonEnergy    continue   end
                amplitudes = PhotoEmission.determineChannels(fLevel, iLevel, settings)
                if   length(amplitudes) == 0   continue   end
                push!( lines, PhotoEmission.Line(iLevel, fLevel, omega, EmProperty(0., 0.),
                                                       EmProperty(0., 0.), amplitudes) )
            end
        end
    end
    return( lines )
end


"""
`PhotoEmission.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid,
                                  settings::PhotoEmission.Settings; output=true)`
    ... as PhotoEmission.computeLines, but on the physical form; a list of
        lines::Array{PhotoEmission.Line,1} is returned.
"""
function computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid,
                            settings::PhotoEmission.Settings; output=true)
    if  settings.calcBiorthogonal
        initialMultiplet, finalMultiplet = BiOrthogonal.computeTransformation(initialMultiplet, finalMultiplet, grid)
    end
    subshellList = Basics.generate(OrderedSubshellList(), finalMultiplet.levels[1].basis, initialMultiplet.levels[1].basis)
    Defaults.setDefaults("relativistic subshell list", subshellList; printout=true)
    println("")
    printstyled("PhotoEmission.computeLines(): The computation of the transition amplitudes and properties starts now ... \n", color=:light_green)
    printstyled("-------------------------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    lines = PhotoEmission.determineLines(finalMultiplet, initialMultiplet, settings)
    if  settings.printBefore    PhotoEmission.displayLines(stdout, lines)    end
    newLines = PhotoEmission.Line[]
    for  line in lines
        push!( newLines, PhotoEmission.computeAmplitudesProperties(line, grid, settings) )
    end
    PhotoEmission.displayRates(stdout, newLines, settings)
    if  settings.calcAnisotropy    PhotoEmission.displayAnisotropies(stdout, newLines, settings)    end
    PhotoEmission.displayLifetimes(stdout, newLines, settings)
    #
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   PhotoEmission.displayRates(iostream, newLines, settings)
        if  settings.calcAnisotropy    PhotoEmission.displayAnisotropies(iostream, newLines, settings)    end
                       PhotoEmission.displayLifetimes(iostream, newLines, settings)
    end
    #
    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`PhotoEmission.computeLinesCascade(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid,
                                         settings::PhotoEmission.Settings; output=true, printout::Bool=true)`
    ... as PhotoEmission.computeLinesCascade, but on the physical form; a list of
        lines::Array{PhotoEmission.Line,1} is returned.

        The "drop a line that does not contribute" test loses its explicit loop: abs2 of an EmPropertyC is an
        EmProperty, so summing it gives both gauges at once and either component answers the question.
"""
function computeLinesCascade(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid,
                                   settings::PhotoEmission.Settings; output=true, printout::Bool=true)
    subshellList = Basics.generate(OrderedSubshellList(), finalMultiplet.levels[1].basis, initialMultiplet.levels[1].basis)
    Defaults.setDefaults("relativistic subshell list", subshellList; printout=false)
    lines = PhotoEmission.determineLines(finalMultiplet, initialMultiplet, settings)
    newLines = PhotoEmission.Line[]
    for  (i,line)  in  enumerate(lines)
        if  rem(i,500) == 0    println("> Radiative line $i:")   end
        newLine = PhotoEmission.computeAmplitudesProperties(line, grid, settings, printout=printout)
        #
        wa = EmProperty(0., 0.);    for  ma in newLine.amplitudes    wa = wa + abs2(ma.amplitude)    end
        if   wa.Coulomb == 0.  &&  wa.Babushkin == 0.    continue    end
        push!( newLines, newLine)
    end
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   PhotoEmission.displayRates(iostream, newLines, settings)    end
    #
    if    output    return( newLines )
    else            return( nothing )
    end
end

end # module

