
"""
`module  JAC.AutoIonization`  
	... a submodel of JAC that contains all methods for computing Auger properties between two given initial- and final-state multiplets.
"""
module AutoIonization


using  Printf, ..AngularMomentum, ..Basics, ..BiOrthogonal, ..Bsplines, ..Continuum, ..Defaults, ..InteractionStrength,
               ..ManyElectron, ..Nuclear, ..Radial, ..SpinAngular, ..TableStrings


"""
`struct  Settings  <:  AbstractProcessSettings`  ... defines a type for the details and parameters of computing Auger lines.

    + calcAnisotropy      ::Bool               ... True, if the intrinsic alpha_2,4 angular parameters are to be calculated, and false
                                                   otherwise.
    + calcTeAuger         ::Bool               
        ... True, if contributions of the two-electron Auger transitions are to be calculated, and false otherwise;
            this flag requires a proper (resonant) Green function that supports the TEA transitions.
    + calcBiorthogonal    ::Bool
        ... True, if the (bound) initial- and final-state multiplets are first brought into a bi-orthogonal
            representation (`BiOrthogonal.computeTransformation`) before the free (Auger) electron orbital is
            generated and the transition amplitudes are evaluated, and false (the default) otherwise. As for
            PhotoIonization, this transforms ONLY the bound spectator orbitals shared between the N-electron
            autoionizing initial multiplet and the (N-1)-electron final-ion multiplet -- the core relaxation
            upon Auger decay -- and must be applied BEFORE the free electron is generated: the free-electron
            partial wave is built fresh per line via `Continuum.generateOrbitalForLevel`, with the matching
            "empty slot" on the initial side filled by `Basics.generateLevelWithExtraSubshell`'s dummy,
            all-zero placeholder orbital; a transformation applied any later would LU-decompose a singular
            per-kappa overlap matrix built partly from that placeholder. See PhotoIonization.Settings'
            calcBiorthogonal docstring for the same reasoning in full, including the still-open question of
            whether BiOrthogonal.computeTransformation is rigorously valid for the resulting N-vs-(N-1)
            electron-count mismatch (empirically it runs and gives plausible results, but has not yet been
            checked against an independent known-answer test).
    + printBefore         ::Bool               ... True, if all energies and lines are printed before their evaluation.
    + lineSelection       ::LineSelection      ... Specifies the selected levels, if any.
    + augerEnergyShift    ::Float64            ... An overall energy shift for all Auger (free-electron) energies.
    + minAugerEnergy      ::Float64            ... Minimum energy of free (Auger) electrons to be included.
    + maxAugerEnergy      ::Float64            ... Maximum energy of free (Auger) electrons to be included.
    + maxKappa            ::Int64              ... Maximum kappa value of partial waves to be included.
    + operator            ::AbstractEeInteraction   
        ... Auger operator that is to be used for evaluating the Auger amplitudes; allowed values are: 
            CoulombInteraction(), BreitInteraction(), ...
    + gMultiplet          ::Multiplet      
        ... Mean-field multiplet of intermediate levels in the computations, sometimes referred to as
            (resonant) Green function.
"""
struct Settings  <:  AbstractProcessSettings
    calcAnisotropy        ::Bool         
    calcTeAuger           ::Bool
    calcBiorthogonal      ::Bool
    printBefore           ::Bool
    lineSelection         ::LineSelection
    augerEnergyShift      ::Float64
    minAugerEnergy        ::Float64
    maxAugerEnergy        ::Float64
    maxKappa              ::Int64
    operator              ::AbstractEeInteraction 
    gMultiplet            ::Multiplet      
end 


"""
`AutoIonization.Settings()`  ... constructor for the default values of AutoIonization line computations.
"""
function Settings()
    Settings(false, false, false, false, LineSelection(), 0., 0., 10e5, 100, CoulombInteraction(), Multiplet())
end


"""
`AutoIonization.Settings(set::AutoIonization.Settings;`

        calcAnisotropy=..,      calcTeAuger..,              calcBiorthogonal=..,    printBefore=..,
        augerEnergyShift=..,    minAugerEnergy=..,          maxAugerEnergy=..,      maxKappa=..,
        operator=..,            gMultiplet=.. )
                    
    ... constructor for modifying the given AutoIonization.Settings by 'overwriting' the previously selected parameters.
"""
function Settings(set::AutoIonization.Settings;
    calcAnisotropy::Union{Nothing,Bool}=nothing,            calcTeAuger::Union{Nothing,Bool}=nothing,
    calcBiorthogonal::Union{Nothing,Bool}=nothing,          printBefore::Union{Nothing,Bool}=nothing,
    lineSelection::Union{Nothing,LineSelection}=nothing,
    augerEnergyShift::Union{Nothing,Float64}=nothing,
    minAugerEnergy::Union{Nothing,Float64}=nothing,         maxAugerEnergy::Union{Nothing,Float64}=nothing,
    maxKappa::Union{Nothing,Int64}=nothing,                 operator::Union{Nothing,AbstractEeInteraction}=nothing,
    gMultiplet::Union{Nothing,Multiplet}=nothing)

    if  isnothing(calcAnisotropy)     calcAnisotropyx    = set.calcAnisotropy    else  calcAnisotropyx    = calcAnisotropy    end
    if  isnothing(calcTeAuger)        calcTeAugerx       = set.calcTeAuger       else  calcTeAugerx       = calcTeAuger       end
    if  isnothing(calcBiorthogonal)   calcBiorthogonalx  = set.calcBiorthogonal  else  calcBiorthogonalx  = calcBiorthogonal  end
    if  isnothing(printBefore)        printBeforex       = set.printBefore       else  printBeforex       = printBefore       end
    if  isnothing(lineSelection)      lineSelectionx     = set.lineSelection     else  lineSelectionx     = lineSelection     end
    if  isnothing(augerEnergyShift)   augerEnergyShiftx  = set.augerEnergyShift  else  augerEnergyShiftx  = augerEnergyShift  end
    if  isnothing(minAugerEnergy)     minAugerEnergyx    = set.minAugerEnergy    else  minAugerEnergyx    = minAugerEnergy    end
    if  isnothing(maxAugerEnergy)     maxAugerEnergyx    = set.maxAugerEnergy    else  maxAugerEnergyx    = maxAugerEnergy    end
    if  isnothing(maxKappa)           maxKappax          = set.maxKappa          else  maxKappax          = maxKappa          end
    if  isnothing(operator)           operatorx          = set.operator          else  operatorx          = operator          end
    if  isnothing(gMultiplet)         gMultipletx        = set.gMultiplet        else  gMultipletx        = gMultiplet        end

    Settings( calcAnisotropyx, calcTeAugerx, calcBiorthogonalx, printBeforex, lineSelectionx, augerEnergyShiftx,
              minAugerEnergyx, maxAugerEnergyx, maxKappax, operatorx, gMultipletx)
end


# `Base.show(io::IO, settings::AutoIonization.Settings)`  ... prepares a proper printout of the variable settings::AutoIonization.Settings.
function Base.show(io::IO, settings::AutoIonization.Settings) 
    println(io, "calcAnisotropy:                $(settings.calcAnisotropy)  ")
    println(io, "calcTeAuger:                   $(settings.calcTeAuger)  ")
    println(io, "calcBiorthogonal:              $(settings.calcBiorthogonal)  ")
    println(io, "printBefore:                   $(settings.printBefore)  ")
    println(io, "lineSelection:                 $(settings.lineSelection)  ")
    println(io, "augerEnergyShift:              $(settings.augerEnergyShift)  ")
    println(io, "minAugerEnergy:                $(settings.minAugerEnergy)  ")
    println(io, "maxAugerEnergy:                $(settings.maxAugerEnergy)  ")
    println(io, "maxKappa:                      $(settings.maxKappa)  ")
    println(io, "operator:                      $(settings.operator)  ")
    println(io, "gMultiplet.name:               $(settings.gMultiplet.name)  ")
end


"""
`struct  AutoIonization.PlasmaSettings  <:  Basics.AbstractLineShiftSettings`  
    ... defines a type for the details and parameters of computing Auger rates with plasma interactions.

    + printBefore         ::Bool             ... True, if all energies and lines are printed before their evaluation.
    + lineSelection       ::LineSelection    ... Specifies the selected levels, if any.
"""
struct PlasmaSettings  <:  Basics.AbstractLineShiftSettings
    printBefore           ::Bool 
    lineSelection         ::LineSelection
end 


"""
`AutoIonization.PlasmaSettings()`  ... constructor for a standard instance of AutoIonization.PlasmaSettings.
"""
function PlasmaSettings()
    PlasmaSettings(true, LineSelection() )
end


# `Base.show(io::IO, settings::AutoIonization.PlasmaSettings)`  ... prepares a proper printout of the settings::AutoIonization.PlasmaSettings.
function Base.show(io::IO, settings::AutoIonization.PlasmaSettings)
    println(io, "printBefore:             $(settings.printBefore)  ")
    println(io, "lineSelection:           $(settings.lineSelection)  ")
end


"""
`struct  AutoIonization.PartialWave`
    ... ONE partial wave of the emitted Auger electron, and its amplitude.

    + kappa          ::Int64                ... partial wave of the free electron.
    + energy         ::Float64              ... energy of the free electron.
    + phase          ::Float64              ... scattering phase; a property of (energy, kappa).
    + amplitude      ::Complex{Float64}     ... Auger amplitude of this partial wave.

        `amplitude` is a plain Complex{Float64} and deliberately NOT an EmPropertyC. The Coulomb-Breit operator has no gauge;
        `settings.operator` selects CoulombInteraction() or BreitInteraction(), which is a choice of OPERATOR, not a pair of representations
        computed together. An EmPropertyC here would assert that two gauges exist. (Measured: this module contains no occurrence of
        `EmProperty` or `gauge` at all.)
"""
struct  PartialWave
    kappa            ::Int64
    energy           ::Float64
    phase            ::Float64
    amplitude        ::Complex{Float64}
end


"""
`struct  AutoIonization.Line`
    ... defines a type for an Auger transition between an initial (autoionizing) and a final level, and carries the emitted electron as a
        list of partial waves; the total rate and the intrinsic angular parameter of the line are stored alongside them.

    + initialLevel   ::Level                            ... initial-(state) level
    + finalLevel     ::Level                            ... final-(state) level
    + electronEnergy ::Float64                          ... Energy of the emitted Auger electron.
    + totalRate      ::Float64                          ... Total rate of this line.
    + angularAlpha   ::Float64                          ... Intrinsic angular parameter alpha_2.
    + partialWaves   ::Array{PartialWave,1}       ... one entry per partial wave of the free electron.

        The total symmetry of the scattering state is NOT stored: it is LevelSymmetry(initialLevel.J, initialLevel.parity) for every partial
        wave, and it is recomputed where needed as LevelSymmetry(line.initialLevel.J, line.initialLevel.parity).
"""
struct  Line
    initialLevel     ::Level
    finalLevel       ::Level
    electronEnergy   ::Float64
    totalRate        ::Float64
    angularAlpha     ::Float64
    partialWaves     ::Array{PartialWave,1}
end


"""
`AutoIonization.amplitude(kind::AbstractEeInteraction, kappa::Int64, phase::Float64, continuumLevel::Level, initialLevel::Level, 
                            grid::Radial.Grid; printout::Bool=true)`  
    ... to compute the kind in  CoulombInteraction(), BreitInteraction(), CoulombBreit(), CoulombGaunt()   Auger amplitude <(alpha_f J_f,
        kappa) J_i || O^(Auger, kind) || alpha_i J_i>  due to the interelectronic interaction for the given final and initial level. A
        value::ComplexF64 is returned.
"""
function amplitude(kind::AbstractEeInteraction, kappa::Int64, phase::Float64, continuumLevel::Level, initialLevel::Level, grid::Radial.Grid; 
                    printout::Bool=true)
    nt = length(continuumLevel.basis.csfs);    ni = length(initialLevel.basis.csfs);    partial = Subshell(9,kappa)
    if  printout  printstyled("Compute ($kind) Auger matrix of dimension $nt x $ni in the continuum- and initial-state bases " *
                                "for the transition [$(initialLevel.index)- ...] and for partial wave $(string(partial)[2:end]) ... ", 
                                color=:light_green)    end
    matrix = zeros(ComplexF64, nt, ni)

    if  initialLevel.basis.subshells == continuumLevel.basis.subshells
        iLevel = initialLevel;   cLevel = continuumLevel
    else
        subshells = Basics.merge(initialLevel.basis.subshells, continuumLevel.basis.subshells)
        Defaults.setDefaults("relativistic subshell list", subshells; printout=false)
        iLevel    = Level(initialLevel, subshells)
        cLevel    = Level(continuumLevel, subshells)
    end

    if      typeof(kind) in [ CoulombInteraction, BreitInteraction, CoulombBreit, CoulombGaunt]        # pure V^Coulomb interaction
    #------------------------------------------------------------------------------------------
        for  r = 1:nt
            for  s = 1:ni
                if  iLevel.basis.csfs[s].J != iLevel.J  ||  iLevel.basis.csfs[s].parity != iLevel.parity      continue    end 
                    # Calculate the spin-angular coefficients
                if  Defaults.saRatip()
                    waR = compute(AngularCoeffsEeRatip2013(), continuumLevel.basis.csfs[r], initialLevel.basis.csfs[s])
                    wa  = waR       
                end
                if  Defaults.saGG()
                        subshellList = cLevel.basis.subshells
                    opa  = SpinAngular.TwoParticleOperator(0, plus, true)
                    waG2 = SpinAngular.computeCoefficients(opa, cLevel.basis.csfs[r], iLevel.basis.csfs[s], subshellList)
                    wa   = [1.0, waG2]
                end
                if  Defaults.saRatip() && Defaults.saGG() && true
                    if  length(waR[2]) != 0     println("\n>> Angular coeffients from Ratip2013   = $(waR[2]) ")    end
                    if  length(waG2)   != 0     println(  ">> Angular coeffients from SpinAngular = $waG2 ")        end
                end

                me = 0.
                for  coeff in wa[2]
                    if   typeof(kind) in [ CoulombInteraction, CoulombBreit, CoulombGaunt]    
                        me = me + coeff.V * InteractionStrength.XL_Coulomb(coeff.nu, 
                                                cLevel.basis.orbitals[coeff.a], cLevel.basis.orbitals[coeff.b],
                                                iLevel.basis.orbitals[coeff.c], iLevel.basis.orbitals[coeff.d], grid)   end
                    if   typeof(kind) in [ BreitInteraction, CoulombBreit, CoulombGaunt]    
                        me = me + coeff.V * InteractionStrength.XL_Breit(coeff.nu, 
                                                cLevel.basis.orbitals[coeff.a], cLevel.basis.orbitals[coeff.b],
                                                iLevel.basis.orbitals[coeff.c], iLevel.basis.orbitals[coeff.d], grid, kind)   end
                    end
                matrix[r,s] = me
            end
        end 
        if  printout  printstyled("done. \n", color=:light_green)    end
        amplitude = transpose(cLevel.mc) * matrix * iLevel.mc 
        amplitude = im^Basics.subshell_l(Subshell(101, kappa)) * exp( -im*phase ) * amplitude


        elseif  kind == "H-E"
    #--------------------
        iLevel = finalLevel;   fLevel = initialLevel
        amplitude = 0.;    error("stop a")
    else    error("stop b")
    end

    return( amplitude )
end


"""
`AutoIonization.amplitude(kind::AbstractEeInteraction, kappa::Int64, phase::Float64, continuumLevel::Level,
                          initialLevel::Level, grid::Radial.Grid, plasmaModel::Basics.AbstractPlasmaModel; printout::Bool=true)`
    ... to compute the Auger amplitude as AutoIonization.amplitude(...) above, but with the electron-electron Coulomb interaction
        Debye-Hueckel-screened according to the given plasma model (InteractionStrength.XL_Coulomb_DH in place of XL_Coulomb) -- the plasma
        parameters thereby enter the transition operator itself, not only the level energies. Only kind = CoulombInteraction() is supported,
        consistent with Basics.compute(...,plasmaModel) which also excludes the Breit interaction from plasma computations. For plasmaModel
        = Basics.NoPlasmaModel(), this delegates directly to the field-free method above, so that the computation of the standard
        (field-free) Auger process is never touched by this plasma-specific code path. A value::ComplexF64 is returned.
"""
function amplitude(kind::AbstractEeInteraction, kappa::Int64, phase::Float64, continuumLevel::Level, initialLevel::Level,
                    grid::Radial.Grid, plasmaModel::Basics.AbstractPlasmaModel; printout::Bool=true)
    if  typeof(plasmaModel) == Basics.NoPlasmaModel
        return( AutoIonization.amplitude(kind, kappa, phase, continuumLevel, initialLevel, grid; printout=printout) )
    end
    typeof(plasmaModel) == Basics.DebyeHueckelModel  ||
        error("Unsupported plasma model = $(plasmaModel)  (only Basics.DebyeHueckelModel is currently supported " *
              "for the plasma-screened Auger amplitude).")
    typeof(kind) == CoulombInteraction  ||
        error("No Breit interaction supported for plasma computations; use kind = CoulombInteraction() for the " *
              "plasma-screened Auger amplitude.")

    nt = length(continuumLevel.basis.csfs);    ni = length(initialLevel.basis.csfs);    partial = Subshell(9,kappa)
    if  printout  printstyled("Compute ($kind, Debye-Hueckel-screened) Auger matrix of dimension $nt x $ni in the " *
                                "continuum- and initial-state bases for the transition [$(initialLevel.index)- ...] " *
                                "and for partial wave $(string(partial)[2:end]) ... ", color=:light_green)    end
    matrix = zeros(ComplexF64, nt, ni)

    if  initialLevel.basis.subshells == continuumLevel.basis.subshells
        iLevel = initialLevel;   cLevel = continuumLevel
    else
        subshells = Basics.merge(initialLevel.basis.subshells, continuumLevel.basis.subshells)
        Defaults.setDefaults("relativistic subshell list", subshells; printout=false)
        iLevel    = Level(initialLevel, subshells)
        cLevel    = Level(continuumLevel, subshells)
    end

    for  r = 1:nt
        for  s = 1:ni
            if  iLevel.basis.csfs[s].J != iLevel.J  ||  iLevel.basis.csfs[s].parity != iLevel.parity      continue    end
            if  Defaults.saRatip()
                waR = compute(AngularCoeffsEeRatip2013(), continuumLevel.basis.csfs[r], initialLevel.basis.csfs[s])
                wa  = waR
            end
            if  Defaults.saGG()
                subshellList = cLevel.basis.subshells
                opa  = SpinAngular.TwoParticleOperator(0, plus, true)
                waG2 = SpinAngular.computeCoefficients(opa, cLevel.basis.csfs[r], iLevel.basis.csfs[s], subshellList)
                wa   = [1.0, waG2]
            end
            me = 0.
            for  coeff in wa[2]
                me = me + coeff.V * InteractionStrength.XL_Coulomb_DH(coeff.nu,
                                        cLevel.basis.orbitals[coeff.a], cLevel.basis.orbitals[coeff.b],
                                        iLevel.basis.orbitals[coeff.c], iLevel.basis.orbitals[coeff.d], grid,
                                        1/plasmaModel.debyeLength)
            end
            matrix[r,s] = me
        end
    end
    if  printout  printstyled("done. \n", color=:light_green)    end
    amplitude = transpose(cLevel.mc) * matrix * iLevel.mc
    amplitude = im^Basics.subshell_l(Subshell(101, kappa)) * exp( -im*phase ) * amplitude

    return( amplitude )
end


"""
`AutoIonization.computeAmplitudesProperties(line::AutoIonization.Line, nm::Nuclear.Model,
        grid::Radial.Grid, nrContinuum::Int64, settings::AutoIonization.Settings; printout::Bool=true,
        nuclearPot::Union{Nothing,Radial.Potential}=nothing, primitives::Union{Nothing,Bsplines.Primitives}=nothing)`
    ... computes the Auger amplitude of every partial wave of the given line and, from them, the total rate and -- if the settings ask for
        it -- the intrinsic angular parameter alpha_2. The rate is summed INCOHERENTLY over the partial waves. The two keywords carry
        the nuclear potential and the B-spline primitives, which are constant for a whole computation and are therefore passed in rather
        than rebuilt per line. A newLine::AutoIonization.Line with all amplitudes and properties evaluated is returned.
"""
function computeAmplitudesProperties(line::AutoIonization.Line, nm::Nuclear.Model, grid::Radial.Grid,
                                           nrContinuum::Int64, settings::AutoIonization.Settings; printout::Bool=true,
                                           nuclearPot::Union{Nothing,Radial.Potential}=nothing,
                                           primitives::Union{Nothing,Bsplines.Primitives}=nothing)
    newPartialWaves = AutoIonization.PartialWave[];   contSettings = Continuum.Settings(false, nrContinuum)
    rate = 0.
    subshellList = Basics.generate(OrderedSubshellList(), line.finalLevel.basis, line.initialLevel.basis)
    nucPot    = isnothing(nuclearPot) ? Nuclear.nuclearPotential(nm, grid) : nuclearPot
    redILevel = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, subshellList)
    newfLevel = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel,   subshellList)
    # The one total symmetry of the scattering state; see the note at PartialWave.
    symi      = LevelSymmetry(line.initialLevel.J, line.initialLevel.parity)

    for  pw in line.partialWaves
        newiLevel = Basics.generateLevelWithExtraSubshell(Subshell(101, pw.kappa), redILevel)
        cOrbital, phase = Continuum.generateOrbitalForLevel(line.electronEnergy, Subshell(101, pw.kappa), newfLevel,
                                                            nm, grid, contSettings; nuclearPot=nucPot,
                                                            primitives=primitives)
        newcLevel  = Basics.generateLevelWithExtraElectron(cOrbital, symi, newfLevel)
        amplitude  = AutoIonization.amplitude(settings.operator, pw.kappa, phase, newcLevel, newiLevel, grid, printout=printout)
        # Calculate two-electron Auger (TEA) contributions if requested; write an extra note if the amplitude is non-zero
        if  settings.calcTeAuger
            if  amplitude != ComplexF64(0.)     @warn ">>> TEA contributions start from non-zero amplitude = $amplitude"    end
            println(">>> Normal Auger for ($symi --> $symi) transition with amplitude = $amplitude")
            amp = AutoIonization.computeTeaAmplitude(settings.operator, newcLevel, settings.gMultiplet,
                                                     newiLevel, grid, printout=printout)
            amplitude = amplitude + amp
        end

        rate = rate + conj(amplitude) * amplitude
        push!(newPartialWaves, AutoIonization.PartialWave(pw.kappa, line.electronEnergy, phase, amplitude))
    end
    totalRate = 2pi * rate;    angularAlpha = 0.
    newLine   = AutoIonization.Line(line.initialLevel, line.finalLevel, line.electronEnergy, totalRate,
                                          angularAlpha, newPartialWaves)

    if  settings.calcAnisotropy    angularAlpha = AutoIonization.computeIntrinsicAlpha(2, newLine)
        newLine = AutoIonization.Line(line.initialLevel, line.finalLevel, line.electronEnergy, totalRate,
                                            real(angularAlpha), newPartialWaves)
    end

    return( newLine )
end


"""
`AutoIonization.computeAmplitudesPropertiesPlasma(line::AutoIonization.Line, nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64,
                                                  settings::AutoIonization.PlasmaSettings, plasmaModel::Basics.AbstractPlasmaModel;
                                                  printout::Bool=true)`
    ... to compute all amplitudes and properties of the given line but for the given plasma model; both the continuum orbital and the Auger
        (e-e Coulomb) transition operator are screened according to plasmaModel, cf. Continuum.generateOrbitalForLevel(...,plasmaModel) and
        AutoIonization.amplitude(...,plasmaModel). A line::AutoIonization.Line is returned for which the amplitudes and properties are now
        evaluated.
"""
function computeAmplitudesPropertiesPlasma(line::AutoIonization.Line, nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64,
                                           settings::AutoIonization.PlasmaSettings, plasmaModel::Basics.AbstractPlasmaModel;
                                           printout::Bool=true)
    newPartialWaves = AutoIonization.PartialWave[];   contSettings = Continuum.Settings(false, nrContinuum);   rate = 0.
    # Define a common subshell list for both multiplets, as in the field-free computeAmplitudesProperties
    subshellList = Basics.generate(OrderedSubshellList(), line.finalLevel.basis, line.initialLevel.basis)
    # Display-only; set once by the driver, not here.  See Defaults.setStandardSubshellList. The one total symmetry of the scattering state:
    # that of the level which autoionizes.
    symi = LevelSymmetry(line.initialLevel.J, line.initialLevel.parity)

    for  pw in line.partialWaves
        newiLevel = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, subshellList)
        newfLevel = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel, subshellList)
        newiLevel = Basics.generateLevelWithExtraSubshell(Subshell(101, pw.kappa), newiLevel)
        cOrbital, phase  = Continuum.generateOrbitalForLevel(line.electronEnergy, Subshell(101, pw.kappa), newfLevel,
                                                              nm, grid, contSettings, plasmaModel)
        newcLevel  = Basics.generateLevelWithExtraElectron(cOrbital, symi, newfLevel)
        amplitude  = AutoIonization.amplitude(CoulombInteraction(), pw.kappa, phase, newcLevel, newiLevel, grid, plasmaModel; printout=printout)
        rate       = rate + conj(amplitude) * amplitude
        push!( newPartialWaves, AutoIonization.PartialWave(pw.kappa, line.electronEnergy, phase, amplitude) )
    end
    totalRate = 2pi* rate;   angularAlpha = 0.
    newLine   = AutoIonization.Line(line.initialLevel, line.finalLevel, line.electronEnergy, totalRate, angularAlpha, newPartialWaves)

    return( newLine )
end


"""
`AutoIonization.computeIntrinsicAlpha(k::Int64, line::AutoIonization.Line)`
    ... computes the intrinsic angular parameter alpha_k of the given line, which determines the anisotropy of the emitted Auger electrons;
        a value::ComplexF64 is returned.

        Pairs are formed ACROSS partial waves -- kappa and kappa' are independent, and that interference is the physical content of the
        parameter. Within a pair, l and j belong to the outer partial wave and lp, jp to the inner one.

        The `if false` block of the flat version -- an unfinished comparison against M. H. Chen, PRA 47 (1993) 3733, which "does not agree
        exactly so far" -- is NOT carried over. It is dead code with a println, and the place to finish it is the flat function.
"""
function computeIntrinsicAlpha(k::Int64, line::AutoIonization.Line)
    wn = 0.;    for  pw in line.partialWaves    wn = wn + conj(pw.amplitude) * pw.amplitude    end
    wa = 0.;    Ji = line.initialLevel.J;    Jf = line.finalLevel.J
    for  pwa in line.partialWaves
        j = AngularMomentum.kappa_j(pwa.kappa);    l = AngularMomentum.kappa_l(pwa.kappa)
        for  pwb in line.partialWaves
            jp = AngularMomentum.kappa_j(pwb.kappa);    lp = AngularMomentum.kappa_l(pwb.kappa)
            wa = wa + sqrt( AngularMomentum.bracket([l, lp, j, jp]) ) *
                      AngularMomentum.ClebschGordan(l, AngularM64(0), lp, AngularM64(0), AngularJ64(k), AngularM64(0)) *
                      AngularMomentum.Wigner_6j(Ji, j, Jf, jp, Ji, AngularJ64(k)) *
                      AngularMomentum.Wigner_6j(l,  j, AngularJ64(1//2), jp, lp, AngularJ64(k)) *
                        pwa.amplitude * conj(pwb.amplitude)
        end
    end

    return( AngularMomentum.phaseFactor([Ji, +1, Jf, +1, AngularJ64(k), -1, AngularJ64(1//2)]) *
            sqrt(Basics.twice(Ji) + 1) * wa / wn )
end


"""
`AutoIonization.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model,
                                   grid::Radial.Grid, settings::AutoIonization.Settings; output=true, printout::Bool=true)`
    ... computes the Auger amplitudes and all requested properties for the lines between the levels of the two given multiplets; this is the
        standard entry point of this module. The continuum orbitals are generated on a shared nuclear potential and B-spline basis. A list
        of lines::Array{AutoIonization.Line,1} is returned if output=true, and nothing otherwise.
"""
function computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                            settings::AutoIonization.Settings; output=true, printout::Bool=true)
    if  settings.calcBiorthogonal
        initialMultiplet, finalMultiplet = BiOrthogonal.computeTransformation(initialMultiplet, finalMultiplet, grid)
    end
    println("")
    printstyled("AutoIonization.computeLines(): The computation of Auger rates and properties starts now ... \n", color=:light_green)
    printstyled("-------------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    lines = AutoIonization.determineLines(finalMultiplet, initialMultiplet, settings)
    # Display-only and the same for every line; set ONCE per computation, never from inside the line loop.
    Defaults.setStandardSubshellList(Basics.generate(OrderedSubshellList(), finalMultiplet.levels[1].basis,
                                                     initialMultiplet.levels[1].basis); printout=false)
    if  settings.printBefore    AutoIonization.displayLines(stdout, lines)    end
    maxEnergy = 0.;   for  line in lines   maxEnergy = max(maxEnergy, line.electronEnergy)   end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
    # Both are constants of the whole computation and are READ-ONLY below, so they are shared safely.
    nuclearPot = Nuclear.nuclearPotential(nm, grid)
    primitives = Bsplines.generatePrimitives(grid)
    newLines = Vector{AutoIonization.Line}(undef, length(lines))
    doPrint  = Threads.nthreads() == 1
    Threads.@threads for  nl  in  eachindex(lines)
        if  doPrint  &&  rem(nl,10) == 0    println("\n>> Auger computations for line No = $nl   ...")     end
        newLines[nl] = AutoIonization.computeAmplitudesProperties(lines[nl], nm, grid, nrContinuum, settings;
                                                                        nuclearPot=nuclearPot, primitives=primitives)
    end
    AutoIonization.displayRates(stdout, newLines, settings)
    AutoIonization.displayLifetimes(stdout, newLines)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   AutoIonization.displayRates(iostream, newLines, settings)
                       AutoIonization.displayLifetimes(iostream, newLines)     end

    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`AutoIonization.computeLinesCascade(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                    settings::AutoIonization.Settings; output::Bool=true, printout::Bool=true)`  
    ... to compute the Auger transition amplitudes and all properties as requested by the given settings. The computations and printout is
        adapted for large cascade computations by including only lines with at least one channel and by sending all printout to a summary
        file only. A list of lines::Array{AutoIonization.Lines} is returned.
"""
function  computeLinesCascade(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                settings::AutoIonization.Settings; output::Bool=true, printout::Bool=true)
    
    lines = AutoIonization.determineLines(finalMultiplet, initialMultiplet, settings)
    # Display-only and the same for every line; set ONCE per computation, never from inside the line loop.
    Defaults.setStandardSubshellList(Basics.generate(OrderedSubshellList(), finalMultiplet.levels[1].basis,
                                                     initialMultiplet.levels[1].basis); printout=false)
    # Display all selected lines before the computations start if  settings.printBefore    AutoIonization.displayLines(stdout, lines)    end
    # Determine maximum energy and check for consistency of the grid
    maxEnergy = 0.;   for  line in lines   maxEnergy = max(maxEnergy, line.electronEnergy)   end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
    # Calculate all amplitudes and requested properties Threaded as computeLines is; see the note there.
    newLines = Vector{AutoIonization.Line}(undef, length(lines))
    doPrint  = printout  &&  Threads.nthreads() == 1
    Threads.@threads for  i  in  eachindex(lines)
        if  doPrint  &&  rem(i,10) == 0    println("> Auger line $i:  ... calculated ")    end
        newLines[i] = AutoIonization.computeAmplitudesProperties(lines[i], nm, grid, nrContinuum, settings,
                                                                 printout=printout)
    end
    # Print all results to a summary file, if requested
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   AutoIonization.displayRates(iostream, newLines, settings)     end

    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`AutoIonization.computeLinesFromOrbitals(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                            settings::AutoIonization.Settings, contOrbitals::Dict{Subshell, Orbital},
                                            contPhases::Dict{Subshell, Float64}; output::Bool=true, printout::Bool=true)`  
    ... to compute the Auger transition amplitudes and all properties as requested by the given settings but by using the given set of
        continuum orbitals and their scattering phases, both keyed on the continuum subshell. The computations and printout is adapted for
        large cascade computations by including only lines with at least one channel and by sending all printout to a summary file only.
        The phase cancels in the total rate, which is an incoherent sum over the partial waves, but NOT in the intrinsic angular
        parameter, where amplitudes of different kappa interfere. A list of lines::Array{AutoIonization.Line,1} is returned.
"""
function  computeLinesFromOrbitals(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                    settings::AutoIonization.Settings, contOrbitals::Dict{Subshell, Orbital},
                                    contPhases::Dict{Subshell, Float64}; output::Bool=true, printout::Bool=true)

    lines = AutoIonization.determineLines(finalMultiplet, initialMultiplet, settings)
    # Display-only and the same for every line; set ONCE per computation, never from inside the line loop.
    Defaults.setStandardSubshellList(Basics.generate(OrderedSubshellList(), finalMultiplet.levels[1].basis,
                                                     initialMultiplet.levels[1].basis); printout=false)
    # Calculate all amplitudes and requested properties
    newLines = AutoIonization.Line[]
    
    for  (i,line)  in  enumerate(lines)
        # Define a common subshell list for both multiplets
        subshellList = Basics.generate(OrderedSubshellList(), line.finalLevel.basis, line.initialLevel.basis)
        # Display-only; set once by the driver, not here.  See Defaults.setStandardSubshellList.
    
        if  rem(i,500) == 0    println("> Auger line $i:  ... calculated ")    end
        # Calculate the individual channels with the given orbitals
        newPartialWaves = AutoIonization.PartialWave[];   rate = 0.
        # The one total symmetry of the scattering state: that of the level which autoionizes.
        symi = LevelSymmetry(line.initialLevel.J, line.initialLevel.parity)

        for  pw in line.partialWaves
            newiLevel  = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, subshellList)
            newfLevel  = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel, subshellList)
            sh         = Subshell(101, pw.kappa)
            if haskey(contOrbitals, sh)   cOrbital = contOrbitals[sh]      else    println(">>> skip Auger channel for $sh");   continue    end
            phase      = contPhases[sh]
            newiLevel  = Basics.generateLevelWithExtraSubshell(sh, newiLevel)
            newcLevel  = Basics.generateLevelWithExtraElectron(cOrbital, symi, newfLevel)
            amplitude  = AutoIonization.amplitude(settings.operator, pw.kappa, phase, newcLevel, newiLevel, grid, printout=printout)
            rate       = rate + conj(amplitude) * amplitude
            push!( newPartialWaves, AutoIonization.PartialWave(pw.kappa, line.electronEnergy, phase, amplitude) )
        end
        totalRate = 2pi* rate;   angularAlpha = 0.
        newLine   = AutoIonization.Line(line.initialLevel, line.finalLevel, line.electronEnergy, totalRate, angularAlpha, newPartialWaves)
        ##
        # if  settings.calcAnisotropy    angularAlpha = AutoIonization.computeIntrinsicAlpha(2, newLine)
        #  newLine   = AutoIonization.Line(line.initialLevel, line.finalLevel, line.electronEnergy, totalRate, real(angularAlpha), newChannels)
        # end
        push!( newLines, newLine)
    end
    # Print all results to a summary file, if requested
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   AutoIonization.displayRates(iostream, newLines, settings)     end

    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`AutoIonization.computeLinesPlasma(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                                   settings::AutoIonization.PlasmaSettings, plasmaModel::Basics.AbstractPlasmaModel; output=true)`
    ... to compute the Auger transition amplitudes and all properties as requested by the given settings and plasma model. A list of
        lines::Array{AutoIonization.Lines} is returned.
"""
function  computeLinesPlasma(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                             settings::AutoIonization.PlasmaSettings, plasmaModel::Basics.AbstractPlasmaModel; output=true)
    println("")
    printstyled("AutoIonization.computeLinesPlasma(): The computation of Auger rates starts now ... \n", color=:light_green)
    printstyled("---------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    augerSettings = AutoIonization.Settings(AutoIonization.Settings(), printBefore=settings.printBefore, lineSelection=settings.lineSelection)
    lines         = AutoIonization.determineLines(finalMultiplet, initialMultiplet, augerSettings)
    # Display-only and the same for every line; set ONCE per computation, never from inside the line loop.
    Defaults.setStandardSubshellList(Basics.generate(OrderedSubshellList(), finalMultiplet.levels[1].basis,
                                                     initialMultiplet.levels[1].basis); printout=false)
    # Display all selected lines before the computations start
    if  settings.printBefore    AutoIonization.displayLines(stdout, lines)    end
    # Determine maximum energy and check for consistency of the grid
    maxEnergy = 0.;   for  line in lines   maxEnergy = max(maxEnergy, line.electronEnergy)   end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
    # Calculate all amplitudes and requested properties
    newLines = AutoIonization.Line[]
    for  line in lines
        newLine = AutoIonization.computeAmplitudesPropertiesPlasma(line, nm, grid, nrContinuum, settings, plasmaModel)
        push!( newLines, newLine)
    end
    # Print all results to screen
    AutoIonization.displayRates(stdout, newLines, augerSettings)
    AutoIonization.displayLifetimes(stdout, newLines)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   AutoIonization.displayRates(iostream, newLines, augerSettings);   AutoIonization.displayLifetimes(iostream, newLines)     end

    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`AutoIonization.computeTeaAmplitude(kind::AbstractEeInteraction, continuumLevel::Level, 
                                    gMultiplet::Multiplet, initialLevel::Level, grid::Radial.Grid; printout::Bool=true)`  
    ... to compute the kind in  CoulombInteraction(), BreitInteraction(), CoulombBreit() total Auger amplitude for the two-electron Auger
        transitions via the gMultiplet as the resonant Green function:
        
        <(alpha_f J_f, kappa) J_i || O^(TEA, kind) || alpha_i J_i> 
        
                    <(alpha_f J_f, kappa) J_i || O^(Auger, kind) || alpha_in J_n> <J_n || O^(e-e, kind) || alpha_i J_i>
                = ---------------------------------------------------------------------------------------------------
                                                                    E_i  -  E_n
                                                                    
        due to the interelectronic interaction as well as the given initial, intermediate (gMultiplet), and final (continuum) levels and the
        given kind of interaction. A value::ComplexF64 is returned.
"""
function computeTeaAmplitude(kind::AbstractEeInteraction, continuumLevel::Level,
                                gMultiplet::Multiplet, initialLevel::Level, grid::Radial.Grid; printout::Bool=true)

    # Always ensure the same subshell list for all initial, intermediate and final (continuum) levels
    subshells  = Basics.merge(initialLevel.basis.subshells, continuumLevel.basis.subshells)
    iLevel     = Level(initialLevel, subshells)
    fLevel     = Level(continuumLevel, subshells)
    nMultiplet = Multiplet(gMultiplet, subshells)
    
    nf = length(fLevel.basis.csfs);    symf = LevelSymmetry(fLevel.J, fLevel.parity)
    ni = length(iLevel.basis.csfs);    symi = LevelSymmetry(iLevel.J, iLevel.parity);    eni = iLevel.energy
    nn = length(nMultiplet.levels[1].basis.csfs)
    
    if  printout   printstyled("Compute two-electron Auger amplitude for the transition [$(iLevel.index)-$(fLevel.index)] ... \n", 
                                color=:light_green)    end
    amplitude = ComplexF64(0.)

    for  nLevel in nMultiplet.levels
        amp  = ComplexF64(0.)
        symn = LevelSymmetry(nLevel.J, nLevel.parity);    enn = nLevel.energy
        if  symf != symn  ||  symn != symi     continue    end

        for  r = 1:nf
            symr = LevelSymmetry(fLevel.basis.csfs[r].J, fLevel.basis.csfs[r].parity);      if  symr != symf    continue    end
            for  s = 1:ni
                syms = LevelSymmetry(iLevel.basis.csfs[s].J, iLevel.basis.csfs[s].parity);  if  syms != symi    continue    end

                #   Compute <alpha_f J_i || V^(e-e) || alpha_n J_i> <alpha_n J_i || V^(e-e) || alpha_i J_i>
                for  t = 1:nn
                    if  nLevel.mc[t] == 0.  continue    end
                    Vee = ManyElectron.matrixElement_Vee(kind, nLevel.basis, t, iLevel.basis, s, grid)
                    Vae = ManyElectron.matrixElement_Vee(kind, fLevel.basis, r, nLevel.basis, t, grid)
                    amp = amp + fLevel.mc[r] * Vae * nLevel.mc[t]^2 * Vee * iLevel.mc[s] / (eni - enn) / 4.
                end
            end
        end
        # Display the amplitude for nLevel, if desired
        if  true   
            println(">>>> TEA amplitude < fLevel=$(continuumLevel.index) [J=$(continuumLevel.J)$(string(continuumLevel.parity))] ||" *
                    " { nLevel=$(nLevel.index) [J=$(nLevel.J)$(string(nLevel.parity))] } ||" *
                    " iLevel=$(initialLevel.index) [$(initialLevel.J)$(string(initialLevel.parity))] >  = $(amp.re)  " *
                    " with  Delta E=$(eni - enn)")    
        end
        amplitude = amplitude + amp
    end
    # if  printout   printstyled("done. \n", color=:light_green)    end
    
    if  printout  
        println(">>>  Total TEA amplitude < fLevel=$(continuumLevel.index) [J=$(continuumLevel.J)$(string(continuumLevel.parity))] ||" *
                " TAE^($kind) ||" *
                " iLevel=$(initialLevel.index) [$(initialLevel.J)$(string(initialLevel.parity))] >  = $amplitude  ")
    end
    
    return( amplitude )
end


"""
`AutoIonization.determineChannels(finalLevel::Level, initialLevel::Level, settings::AutoIonization.Settings)`
    ... determines the partial waves of the emitted Auger electron for the given transition, by applying
        AngularMomentum.allowedKappaSymmetries(symi, symf) and discarding every kappa with abs(kappa) > settings.maxKappa. No total symmetry
        is stored with a partial wave: it is symi for all of them and hence a property of the line. An Array{AutoIonization.PartialWave,1}
        is returned, with all amplitudes still zero.
"""
function determineChannels(finalLevel::Level, initialLevel::Level, settings::AutoIonization.Settings)
    partialWaves = AutoIonization.PartialWave[]
    symi         = LevelSymmetry(initialLevel.J, initialLevel.parity)
    symf         = LevelSymmetry(finalLevel.J,   finalLevel.parity)
    for  kappa in AngularMomentum.allowedKappaSymmetries(symi, symf)
        if  abs(kappa) > settings.maxKappa      continue    end
        push!(partialWaves, AutoIonization.PartialWave(kappa, 0., 0., Complex(0.)))
    end
    return( partialWaves )
end


"""
`AutoIonization.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet,
                                     settings::AutoIonization.Settings)`
    ... selects the Auger lines to be computed by forming all pairs of an initial and a final level that pass the line selection and for
        which the electron energy -- the transition energy shifted by settings.augerEnergyShift -- is positive. An
        Array{AutoIonization.Line,1} is returned, with all amplitudes still zero.
"""
function determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::AutoIonization.Settings)
    lines = AutoIonization.Line[]
    augerEnergyShift = Defaults.convertUnits("energy: to atomic", settings.augerEnergyShift)
    for  iLevel  in  initialMultiplet.levels
        for  fLevel  in  finalMultiplet.levels
            if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                energy = iLevel.energy - fLevel.energy + augerEnergyShift
                if   energy < 0.01                                                             continue   end
                if   energy < settings.minAugerEnergy  ||  energy > settings.maxAugerEnergy    continue   end
                partialWaves = AutoIonization.determineChannels(fLevel, iLevel, settings)
                push!( lines, AutoIonization.Line(iLevel, fLevel, energy, 0., 0., partialWaves) )
            end
        end
    end

    return( lines )
end


"""
`AutoIonization.displayLifetimes(stream::IO, lines::Array{AutoIonization.Line,1})`
    ... derives the lifetime of every autoionizing level from the total Auger rate of all lines that start from it, and lists it together
        with the corresponding level width. A neat table is printed to stream; nothing::Nothing is returned.
"""
function  displayLifetimes(stream::IO, lines::Array{AutoIonization.Line,1})
    nx = 104
    println(stream, " ")
    println(stream, "  Auger lifetimes, total rates and widths:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level";    na=2);                           sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center( 8, "J^P";      na=4);                           sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(12, "Lifetime"; na=4);               
    sb = sb * TableStrings.center(12,TableStrings.inUnits("time"); na=4)
    sa = sa * TableStrings.center(14, "Total rate"; na=6);               
    sb = sb * TableStrings.center(14,TableStrings.inUnits("rate"); na=4)
    sa = sa * TableStrings.center(42, "Widths"; na=2);       
    sb = sb * TableStrings.center(42, "  Hartrees         Kaysers           eV"; na=2)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 

    notYetDone = trues(1000)
    for  line in lines
        if  notYetDone[line.initialLevel.index]
            notYetDone[line.initialLevel.index] = false
            totalRate = 0.
            for  ln in lines
                if  ln.initialLevel.index == line.initialLevel.index   totalRate = totalRate + ln.totalRate    end
            end
            sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
            sa = sa * TableStrings.center(10, TableStrings.level(line.initialLevel.index); na=2)
            sa = sa * TableStrings.center( 8, string(isym); na=4)
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("time: from atomic",  1/totalRate))            * "     "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("rate: from atomic",    totalRate))            * "      "
            sa = sa * @sprintf("%.6e", totalRate)                                                           * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic to Kayser",  totalRate))  * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic to eV",      totalRate))  * "    "
            println(stream, sa)
        end
    end
    println(stream, "  ", TableStrings.hLine(nx))

    return( nothing )
end


"""
`AutoIonization.displayLines(stream::IO, lines::Array{AutoIonization.Line,1})`
    ... lists the Auger lines that have been selected, together with the electron energy and the kappa values of the partial waves that will
        be computed for each of them, and counts the lines and partial waves involved. The total symmetry is taken from the autoionizing
        level, which is where it belongs. A neat table is printed to stream; nothing::Nothing is returned.
"""
function  displayLines(stream::IO, lines::Array{AutoIonization.Line,1})
    nx = 150;   noLines = 0;    noChannels = 0
    println(stream, " ")
    println(stream, "  Selected Auger lines:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=2);                                sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                                sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(14, "Energy"; na=4);              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(14, "Energy e_A"; na=4);              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.flushleft(37, "List of kappas and total symmetries"; na=4)  
    sb = sb * TableStrings.flushleft(37, "partial (total J^P)                "; na=4)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 

    for  line in lines
        noLines = noLines + 1;    noChannels = noChannels + length(line.partialWaves)
        sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", line.initialLevel.energy))  * "    "
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", line.electronEnergy))       * "   "
        kappaSymmetryList = Tuple{Int64,LevelSymmetry}[]
        # One entry per partial wave; the total symmetry is that of the autoionizing level, the same for all.
        symi = LevelSymmetry(line.initialLevel.J, line.initialLevel.parity)
        for  pw in line.partialWaves
            push!( kappaSymmetryList, (pw.kappa, symi) )
        end
        sa = sa * TableStrings.kappaSymmetryTupels(80, kappaSymmetryList)
        println(stream,  sa )
    end
    println(stream, "  ", TableStrings.hLine(nx), "\n\n  A total of $noLines lines with $noChannels Auger channels will be compute. \n")

    return( nothing )
end


"""
`AutoIonization.displayRates(stream::IO, lines::Array{AutoIonization.Line,1}, settings::AutoIonization.Settings)`
    ... lists, for every Auger line, the transition energy, the energy of the emitted electron and the total Auger rate; the intrinsic
        angular parameter alpha_2 is shown in addition whenever settings.calcAnisotropy is true. A neat table is printed to stream;
        nothing::Nothing is returned.
"""
function  displayRates(stream::IO, lines::Array{AutoIonization.Line,1}, settings::AutoIonization.Settings)
    nx = 106
    println(stream, " ")
    if  settings.calcAnisotropy    println(stream, "  Auger rates and intrinsic angular parameters: \n")
    else                           println(stream, "  Auger rates (without angular parameters): \n")        end
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=2);                         sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                         sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(12, "Energy"   ; na=2);               
    sb = sb * TableStrings.center(12,TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(14, "Electron energy"   ; na=2);               
    sb = sb * TableStrings.center(14,TableStrings.inUnits("energy"); na=2)
    sa = sa * TableStrings.center(14, "Auger rate"; na=2);       
    sb = sb * TableStrings.center(14, TableStrings.inUnits("rate"); na=2)
    sa = sa * TableStrings.center(15, "alpha_2"; na=2);                           sb = sb * TableStrings.hBlank(18)     
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 

    for  line in lines
        sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.initialLevel.energy))  * "    "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.electronEnergy))       * "    "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("rate: from atomic", line.totalRate))              * "    "
        sa = sa * TableStrings.flushright(13, @sprintf("%.4e", line.angularAlpha))            * "    "
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))

    return( nothing )
end

end # module
