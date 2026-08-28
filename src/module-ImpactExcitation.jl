
"""
`module  JAC.ImpactExcitation`  
... a submodel of JAC that contains all methods for computing electron impact excitation cross sections and collision strengths.
"""
module ImpactExcitation 

using   Distributed, FastGaussQuadrature, GSL, SpecialFunctions, Printf,
        ..AngularMomentum, ..Basics, ..Bsplines, ..Continuum, ..Defaults, ..InteractionStrength, ..ManyElectron,
        ..Nuclear, ..Radial, ..SpinAngular, ..TableStrings, ..RadialIntegrals

## Relative change of sum |amplitude|^2 contributed by the last partial wave, below which the kappa summation
## counts as converged. Used BOTH to terminate that summation in computeAmplitudesProperties() and to flag the
## lines that never reached it in displayResults() -- the two must not be allowed to drift apart.
const convergenceCriterion = 1.0e-5



"""
`struct  ImpactExcitation.RateCoefficients`     ... Defines a type for the output results from excitation rate or
                                                    effective collision strengths calculations
                                                    
    + initialLevel        ::Level               ... initial- (bound-state) level
    + finalLevel          ::Level               ... final- (bound-state) level
    + temperatures        ::Array{Float64,1}    ... Temperatures in [K] to calculate excitation rates and effective collision strengths
    + alphas              ::Array{Float64,1}    ... Excitation rate coefficients in [cm^3/s] for the input temperatures
    + effOmegas           ::Array{Float64,1}    ... Effective collision strengths for the input temperatures
"""
struct  RateCoefficients   
    initialLevel        ::Level 
    finalLevel          ::Level 
    temperatures        ::Array{Float64,1}
    alphas              ::Array{Float64,1}
    effOmegas           ::Array{Float64,1}
end


"""
`ImpactExcitation.RateCoefficients()`  ... constructor for the default values of electron-impact excitation rate coefficient or
                                            effective collision strength computations.
"""
function RateCoefficients()
    RateCoefficients(Level(), Level(), Float64[], Float64[], Float64[])
end


"""
`struct  ImpactExcitation.Settings  <:  AbstractProcessSettings`  ... defines a type for the details and parameters of computing 
                                                                                                electron-impact excitation lines.

    + lineSelection            ::LineSelection      ... Specifies the selected levels, if any.
    + electronEnergies         ::Array{Float64,1}   ... List of impact-energies of the incoming elecgtrons (in user-defined units).
    + energyShift              ::Float64            ... An overall energy shift for all transitions |i> --> |f>.
    + maxKappa                 ::Int64              ... Maximum kappa value of partial waves to be included.
    + calcRateCoefficient      ::Bool               ... True, if the plasma rate coefficients to be calculated, false otherwise.
    + maxEnergyMultiplier      ::Float64            ... Maximum initial electron energy for eff. collision strength integration.
                                     (maxEnergyMultiplier * Excitation threshold energy), after this assymptotic limit is applied.
    + numElectronEnergies      ::Int64              ... No. of different electron energy points at which collision strengths to compute
                                                        in electron energy range [0, (maxEnergyMultiplier * Excitation threshold energy)]
    + temperatures             ::Array{Float64, 1}  ... Electron temperatures [K] for eff. collision strengths and rate coefficients.
    + printBefore              ::Bool               ... True, if all energies and lines are printed before their evaluation.
    + operator                 ::AbstractEeInteraction   
        ... Interaction operator that is to be used for evaluating the e-e interaction amplitudes; allowed values are: 
            CoulombInteraction(), BreitInteraction(), ...
"""
struct Settings  <:  AbstractProcessSettings
    lineSelection             ::LineSelection  
    electronEnergies          ::Array{Float64,1}
    energyShift               ::Float64   
    maxKappa                  ::Int64
    calcRateCoefficient       ::Bool 
    maxEnergyMultiplier       ::Float64
    numElectronEnergies       ::Int64
    temperatures              ::Array{Float64, 1}
    printBefore               ::Bool  
    operator                  ::AbstractEeInteraction 
end 


"""
`ImpactExcitation.Settings()`  ... constructor for the default values of electron-impact excitation line computations.
"""
function Settings()
    Settings( LineSelection(), Float64[], 0., 300, false, 30.0, 6, [5e4, 1e5, 2e5, 3e5, 5e5, 7e5, 1e6], false, CoulombInteraction())
end


"""
`ImpactExcitation.Settings(set::ImpactExcitation.Settings;`

    lineSelection...,  electronEnergies..., energyShift..., maxKappa...,
    calcRateCoefficient..., maxEnergyMultiplier..., numElectronEnergies..., temperatures...,
    printBefore..., operator...)

    ... constructor for modifying the given ImpactExcitation.Settings by 'overwriting' the previously selected parameters.
"""
function Settings(set::ImpactExcitation.Settings;
    lineSelection::Union{Nothing, LineSelection} = nothing, 
    electronEnergies::Union{ Nothing, Array{Float64,1} } = nothing,
    energyShift::Union{Nothing, Float64} = nothing,
    maxKappa::Union{Nothing, Int64} = nothing,
    calcRateCoefficient::Union{Nothing, Bool} = nothing,
    maxEnergyMultiplier::Union{Nothing, Float64} = nothing,
    numElectronEnergies::Union{Nothing, Int64} = nothing,
    temperatures::Union{ Nothing, Array{Float64,1} } = nothing,
    printBefore::Union{Nothing, Bool} = nothing,
    operator::Union{Nothing, AbstractEeInteraction} = nothing )

    if isnothing(lineSelection)     lineSelectionx = set.lineSelection       else lineSelectionx = lineSelection        end
    if isnothing(electronEnergies)  electronEnergiesx = set.electronEnergies else electronEnergiesx = electronEnergies  end
    if isnothing(energyShift)       energyShiftx = set.energyShift           else energyShiftx = energyShift            end
    if isnothing(maxKappa)          maxKappax = set.maxKappa                 else maxKappax = maxKappa                  end
    if isnothing(calcRateCoefficient)   calcRateCoefficientx = set.calcRateCoefficient else 
                                                                                calcRateCoefficientx = calcRateCoefficient end
    if isnothing(maxEnergyMultiplier) maxEnergyMultiplierx = set.maxEnergyMultiplier else maxEnergyMultiplierx = maxEnergyMultiplier end
    if isnothing(numElectronEnergies) numElectronEnergiesx = set.numElectronEnergies else 
                                                                             numElectronEnergiesx = numElectronEnergies end
    if isnothing(temperatures)      temperaturesx = set.temperatures         else temperaturesx = temperatures          end
    if isnothing(printBefore)       printBeforex = set.printBefore           else printBeforex = printBefore            end
    if isnothing(operator)          operatorx = set.operator                 else operatorx = operator                  end

    Settings( lineSelectionx, electronEnergiesx, energyShiftx, maxKappax, calcRateCoefficientx, maxEnergyMultiplierx, 
                                                            numElectronEnergiesx, temperaturesx, printBeforex, operatorx )

end


# `Base.show(io::IO, settings::ImpactExcitation.Settings)`  ... prepares a proper printout of the variable settings::ImpactExcitation.Settings.
function Base.show(io::IO, settings::ImpactExcitation.Settings) 
    println(io, "lineSelection:              $(settings.lineSelection)  ")
    println(io, "electronEnergies:           $(settings.electronEnergies)  ")
    println(io, "energyShift:                $(settings.energyShift)  ")
    println(io, "maxKappa:                   $(settings.maxKappa)  ")
    println(io, "calcRateCoefficient:        $(settings.calcRateCoefficient)  ")
    println(io, "maxEnergyMultiplier:        $(settings.maxEnergyMultiplier)  ")
    println(io, "numElectronEnergies:        $(settings.numElectronEnergies)  ")
    println(io, "temperatures:               $(settings.temperatures)  ")
    println(io, "printBefore:                $(settings.printBefore)  ")
    println(io, "operator:                   $(settings.operator)  ")
end


"""
`struct  ImpactExcitation.Channel`  
    ... defines a type for a electron-impact excitaiton channel to help characterize the incoming and outgoing (continuum) states of 
        many electron-states with a single free electron

    + initialKappa     ::Int64              ... partial-wave of the incoming free electron
    + finalKappa       ::Int64              ... partial-wave of the outgoing free electron
    + symmetry         ::LevelSymmetry      ... total angular momentum and parity of the scattering state
    + initialPhase     ::Float64            ... phase of the incoming partial wave
    + finalPhase       ::Float64            ... phase of the outgoing partial wave
    + amplitude        ::Complex{Float64}   ... Collision amplitude associated with the given channel.
"""
struct  Channel
    initialKappa       ::Int64 
    finalKappa         ::Int64 
    symmetry           ::LevelSymmetry
    initialPhase       ::Float64
    finalPhase         ::Float64
    amplitude          ::Complex{Float64}
end


# `Base.show(io::IO, channel::ImpactExcitation.Channel)`  ... prepares a proper printout of the variable channel::ImpactExcitation.Channel.
function Base.show(io::IO, channel::ImpactExcitation.Channel) 
    println(io, "initialKappa:       $(channel.initialKappa)  ")
    println(io, "finalKappa:         $(channel.finalKappa)  ")
    println(io, "symmetry:           $(channel.symmetry)  ")
    println(io, "initialPhase:       $(channel.initialPhase)  ")
    println(io, "finalPhase:         $(channel.finalPhase)  ")
    println(io, "amplitude:          $(channel.amplitude)  ")
end


"""
`struct  ImpactExcitation.Line`  
    ... defines a type for a electron-impact excitation line that may include the definition of channels and their corresponding
                                                                                                                        amplitudes.

    + initialLevel           ::Level         ... initial- (bound-state) level
    + finalLevel             ::Level         ... final- (bound-state) level
    + initialElectronEnergy  ::Float64       ... energy of the incoming (initial-state) free-electron
    + finalElectronEnergy    ::Float64       ... energy of the outgoing (final-state) free-electron
    + crossSection           ::Float64       ... total cross section of this line
    + collisionStrength      ::Float64       ... total collision strength of this line
    + channels               ::Array{ImpactExcitation.Channel,1}  ... List of ImpactExcitation channels of this line.
    + convergence            ::Float64       ... convergence of calculation
"""
struct  Line
    initialLevel             ::Level
    finalLevel               ::Level
    initialElectronEnergy    ::Float64
    finalElectronEnergy      ::Float64
    crossSection             ::Float64 
    collisionStrength        ::Float64 
    channels                 ::Array{ImpactExcitation.Channel,1}
    convergence              ::Float64
end 


"""
`ImpactExcitation.Line()`  ... 'empty' constructor for an electron-impact excitation line between a specified initial and final level.
"""
function Line()
    Line(Level(), Level(), 0., 0., 0., 0., ImpactExcitation.Channel[], 1.0 )
end


"""
`ImpactExcitation.Line(initialLevel::Level, finalLevel::Level, crossSection::Float64)`  
    ... constructor for an electron-impact excitation line between a specified initial and final level.
"""
function Line(initialLevel::Level, finalLevel::Level, crossSection::Float64)
    Line(initialLevel, finalLevel, 0., 0., crossSection, 0., ImpactExcitation.Channel[], 1.0 )
end


# `Base.show(io::IO, line::ImpactExcitation.Line)`  ... prepares a proper printout of the variable line::ImpactExcitation.Line.
function Base.show(io::IO, line::ImpactExcitation.Line) 
    println(io, "initialLevel:            $(line.initialLevel)  ")
    println(io, "finalLevel:              $(line.finalLevel)  ")
    println(io, "initialElectronEnergy:   $(line.initialElectronEnergy)  ")
    println(io, "finalElectronEnergy:     $(line.finalElectronEnergy)  ")
    println(io, "crossSection:            $(line.crossSection)  ")
    println(io, "collisionStrength:       $(line.collisionStrength)  ")
    println(io, "channels:                $(line.channels)  ")
end


"""
`ImpactExcitation.amplitude(kind::AbstractEeInteraction, channel::ImpactExcitation.Channel, cFinalLevel::Level, cInitialLevel::Level, 
                            grid::Radial.Grid; printout::Bool=true)`  
    ... to compute the kind in  CoulombInteraction(), BreitInteraction(), CoulombBreit() electron-impact interaction amplitude 
        <(alpha_f J_f, kappa_f) J_t || O^(e-e, kind) || (alpha_i J_i, kappa_i) J_t>  due to the interelectronic interaction for 
        the given final and initial (continuum) level. A value::ComplexF64 is returned.
"""
function amplitude(kind::AbstractEeInteraction, channel::ImpactExcitation.Channel, cFinalLevel::Level, cInitialLevel::Level, 
                    grid::Radial.Grid; printout::Bool=true)
    nf = length(cFinalLevel.basis.csfs);      fPartial = Subshell(9,channel.finalKappa)         
    ni = length(cInitialLevel.basis.csfs);    iPartial = Subshell(9,channel.initialKappa)    
    
    if  printout  printstyled("Compute ($kind) e-e matrix of dimension $nf x $ni in the final- and initial-state (continuum) bases " *
                                "for the transition [$(cInitialLevel.index)- ...] " * 
                                "and for partial waves $(string(fPartial)[2:end]),  $(string(iPartial)[2:end])... ", color=:light_green) end
    matrix = zeros(Float64, nf, ni)
    #
    if  cInitialLevel.basis.subshells == cFinalLevel.basis.subshells
        iLevel = cInitialLevel;   fLevel = cFinalLevel
    else
        subshells = Basics.merge(cInitialLevel.basis.subshells, cFinalLevel.basis.subshells)
        iLevel    = Level(cInitialLevel, subshells)
        fLevel    = Level(cFinalLevel,   subshells)
    end
    #
    if      typeof(kind) in [ CoulombInteraction, BreitInteraction, CoulombBreit]        ## pure V^Coulomb interaction
    #--------------------------------------------------------------------------
        for  r = 1:nf
            for  s = 1:ni
                #if ( fLevel.mc[r] == 0.0 && iLevel.mc[s] == 0.0 ) continue end
                if  iLevel.basis.csfs[s].J != iLevel.J  ||  iLevel.basis.csfs[s].parity != iLevel.parity      continue    end 
                subshellList = fLevel.basis.subshells
                op2  = SpinAngular.TwoParticleOperator(0, plus)
                wa   = SpinAngular.computeCoefficients(op2, fLevel.basis.csfs[r], iLevel.basis.csfs[s], subshellList)
                #
                me = 0.

                for  coeff in wa
                    # The exchnage T-matrices are considered up to 40 partial waves
                    if ( (coeff.c.n == 101) && (coeff.nu > 40) )  continue end

                    if   typeof(kind) in [ CoulombInteraction, CoulombBreit]    
                        me = me + coeff.V * InteractionStrength.XL_Coulomb(coeff.nu, 
                                                fLevel.basis.orbitals[coeff.a], fLevel.basis.orbitals[coeff.b],
                                                iLevel.basis.orbitals[coeff.c], iLevel.basis.orbitals[coeff.d], grid)   end
                    if   typeof(kind) in [ BreitInteraction, CoulombBreit]    
                        me = me + coeff.V * InteractionStrength.XL_Breit(coeff.nu, 
                                                fLevel.basis.orbitals[coeff.a], fLevel.basis.orbitals[coeff.b],
                                                iLevel.basis.orbitals[coeff.c], iLevel.basis.orbitals[coeff.d], grid,
                                                                                                   CoulombBreit(1.0))   end
                end
                matrix[r,s] = me
            end
        end 
        if  printout  printstyled("done. \n", color=:light_green)    end
        amplitude = transpose(fLevel.mc) * matrix * iLevel.mc 
        amplitude = im^( -1.0 * Basics.subshell_l(Subshell(102, channel.finalKappa)) )   * exp(  im*channel.finalPhase )   * 
                    im^Basics.subshell_l(Subshell(101, channel.initialKappa)) * exp(  im*channel.initialPhase ) * amplitude
        #
        #
        elseif  kind == "H-E"
    #--------------------
        amplitude = 0.;    error("stop a")
    else    error("stop b")
    end
    
    return( amplitude )
end


"""
`ImpactExcitation.computeAmplitudesProperties(line::ImpactExcitation.Line, nm::Nuclear.Model, grid::Radial.Grid,
                                                nrContinuum::Int64, settings::ImpactExcitation.Settings; printout::Bool=true,
                                                nuclearPot::Union{Nothing,Radial.Potential}=nothing,
                                                primitives::Union{Nothing,Bsplines.Primitives}=nothing)`
    ... to compute all amplitudes and properties of the given line; a line::ImpactExcitation.Line is returned for which
    the amplitudes and properties are now evaluated. The two keywords carry quantities that are CONSTANT for a whole
    computation and are otherwise rebuilt for every line and every partial wave: the nuclear potential depends only on
    the nuclear model and the grid, the B-spline primitives only on the grid. Nothing is cached; omitting them
    reproduces the previous behaviour exactly.
"""
function  computeAmplitudesProperties(line::ImpactExcitation.Line, nm::Nuclear.Model, grid::Radial.Grid,
                                        nrContinuum::Int64, settings::ImpactExcitation.Settings; printout::Bool=true,
                                        nuclearPot::Union{Nothing,Radial.Potential}=nothing,
                                        primitives::Union{Nothing,Bsplines.Primitives}=nothing)
    newChannels = ImpactExcitation.Channel[];
    contSettings = Continuum.Settings(false, nrContinuum);   cross = 0.;   coll = 0.; convergence = 1.
    conv = 0.; conv0 = 0. ; n = 0

    # Define a common subshell list for both multiplets
    subshellList = Basics.generate(OrderedSubshellList(), line.finalLevel.basis, line.initialLevel.basis)
    ## The standard subshell list is display-only and identical for every line (verified over 50 level pairs),
    ## so it is set once in the drivers rather than here -- this function runs inside a @threads loop, and a
    ## write to a global from there is a race.  See Defaults.setStandardSubshellList.
    ## Both symmetry-reduced levels are properties of the line alone; see the note in the innermost loop.
    redILevel = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, subshellList)
    redFLevel = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel,   subshellList)

    # First determine a common set of continuum orbitals for the incoming and outgoing electron
    ciOrbitals = Dict{Subshell, Orbital}();     ciPhases = Dict{Subshell, Float64}()
    cfOrbitals = Dict{Subshell, Orbital}();     cfPhases = Dict{Subshell, Float64}()
    #
    symi = LevelSymmetry(line.initialLevel.J, line.initialLevel.parity);    
    symf = LevelSymmetry(line.finalLevel.J, line.finalLevel.parity) 
    #
    kappa = 1
    while kappa <= settings.maxKappa
        for  inKappa = -kappa:kappa
            if  inKappa == 0    continue    end
            insymtList = AngularMomentum.allowedTotalSymmetries(symi, inKappa)
            for  outKappa = -(kappa):kappa
                if  outKappa == 0  || ( abs(inKappa) < kappa && abs(outKappa) < kappa )  continue    end
                outsymtList = AngularMomentum.allowedTotalSymmetries(symf, outKappa)
                for  symt in insymtList
                    if symt in outsymtList
                        channel = ImpactExcitation.Channel( inKappa, outKappa, symt, 0., 0.,Complex(0.))
                        # Generate the continuum orbitals if they were not yet generated before
                        iSubshell  = Subshell(101, channel.initialKappa)
                        if  !haskey(ciOrbitals, iSubshell)
                            ciOrbital, iPhase     = Continuum.generateOrbitalForLevel(line.initialElectronEnergy, iSubshell,
                                                                                            redILevel, nm, grid, contSettings;
                                                                                            nuclearPot=nuclearPot,
                                                                                            primitives=primitives)
                            ciOrbitals[iSubshell] = ciOrbital;      ciPhases[iSubshell] = iPhase
                            println("\n Initial channel")
                        end

                        fSubshell  = Subshell(102, channel.finalKappa)
                        if  !haskey(cfOrbitals, fSubshell)
                            cfOrbital, fPhase     = Continuum.generateOrbitalForLevel(line.finalElectronEnergy,   fSubshell,
                                                                                            redFLevel, nm, grid, contSettings;
                                                                                            nuclearPot=nuclearPot,
                                                                                            primitives=primitives)
                            cfOrbitals[fSubshell] = cfOrbital;      cfPhases[fSubshell] = fPhase
                            println("\n Final channel")
                        end
                    #end
                    
                    #for channel in line.channels
                        ## The two symmetry-reduced levels depend on the line but on neither partial wave nor
                        ## total symmetry, so they are formed ONCE above rather than in this innermost loop,
                        ## which runs O(maxKappa^3) times.  The four calls below rebind newiLevel/newfLevel to
                        ## fresh Levels (generateLevelWithExtraElectron/-Subshell copy their argument), so
                        ## redILevel and redFLevel are never modified.
                        newiLevel  = redILevel
                        newfLevel  = redFLevel
                        iSubshell  = Subshell(101, channel.initialKappa)
                        fSubshell  = Subshell(102, channel.finalKappa)
                        ciOrbital  = ciOrbitals[iSubshell];     iPhase = ciPhases[iSubshell]
                        cfOrbital  = cfOrbitals[fSubshell];     fPhase = cfPhases[fSubshell]
                        newiLevel = Basics.generateLevelWithExtraElectron(ciOrbital, channel.symmetry, newiLevel)
                        newiLevel = Basics.generateLevelWithExtraSubshell(fSubshell,   newiLevel)
                        newfLevel = Basics.generateLevelWithExtraSubshell(iSubshell, newfLevel)
                        newfLevel = Basics.generateLevelWithExtraElectron(cfOrbital, channel.symmetry, newfLevel)
                        newChannel = ImpactExcitation.Channel(channel.initialKappa, channel.finalKappa, channel.symmetry, 
                                                                                                            iPhase, fPhase, 0.)
                        #
                        amplitude  = ImpactExcitation.amplitude(settings.operator, newChannel, newfLevel, newiLevel, grid, 
                                                                                                                printout=printout)

                        coll       = coll + AngularMomentum.bracket([channel.symmetry.J]) * ( conj(amplitude) * amplitude ).re

                        conv += (conj(amplitude) * amplitude).re

                        push!( newChannels, ImpactExcitation.Channel(newChannel.initialKappa, newChannel.finalKappa, 
                                                                            newChannel.symmetry, iPhase, fPhase, amplitude) )
                    end
                end
            end
        end
        
        # Checking convergence
        convergence = abs(conv - conv0)/conv
        if convergence < convergenceCriterion
            printstyled("\nConvergence Achieved "; bold=true, underline=true, color=:light_red); 
            if n == 5 
                printstyled("\nConvergence Achieved "; bold=true, underline=true, color=:light_red); break
            end
            n += 1
        end
        conv0 = conv
        kappa += 1
    end

    wc = Defaults.getDefaults("speed of light: c")

    qa = sqrt( (line.initialElectronEnergy/wc)^2 + 2*line.initialElectronEnergy )
    qb = sqrt( (line.finalElectronEnergy/wc)^2 + 2*line.finalElectronEnergy )

    initialEnergy = sqrt( (qa*wc)^2 + wc^4 );    finalEnergy = sqrt( (qb*wc)^2 + wc^4 )

    N = 4. * sqrt( (initialEnergy + wc^2) * ( finalEnergy + wc^2 ) / ( initialEnergy * finalEnergy ) ) / qa / qb
    N = N / ( sqrt( qa / pi / line.initialElectronEnergy) * sqrt( qb / pi / line.finalElectronEnergy) )^2

    coll = coll * N 
    cross = coll * pi / ( (qa^2) * (Basics.twice(line.initialLevel.J) + 1) )

    newLine = ImpactExcitation.Line( line.initialLevel, line.finalLevel, line.initialElectronEnergy, line.finalElectronEnergy,
                                        cross, coll, newChannels, convergence )
    ImpactExcitation.displayResults(newLine)
    return( newLine )
end


"""
`ImpactExcitation.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                settings::ImpactExcitation.Settings; output=true)`  
    ... to compute the electron-impact excitation transition amplitudes and all properties as requested by the given settings. 
        A list of lines::Array{ImpactExcitation.Lines} is returned.
"""
function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::ImpactExcitation.Settings; 
                        output=true)
    println("")
    printstyled("ImpactExcitation.computePathways(): The computation of electron-impact excitation cross sections starts now ... \n", color=:light_green)
    printstyled("--------------------------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    #
    lines = ImpactExcitation.determineLines(finalMultiplet, initialMultiplet, settings)
    # Display all selected lines before the computations start
    if  settings.printBefore    ImpactExcitation.displayLines(lines)    end
    ## Determine the maximum (continuum) electron energy and check the grid for consistency with it. This module used to
    ## skip this check -- alone among the continuum-based process modules -- and to normalize at grid.NoPoints-5 instead
    ## of the nrContinuum returned here. Both were wrong: an under-resolved continuum orbital does not fail loudly, it
    ## returns cross sections that drift by tens of percent and then diverge by 30+ orders of magnitude once the grid can
    ## no longer represent the oscillation. Continuum.gridConsistency() refuses such a grid up front and says why.
    maxEnergy = 0.;   for  line in lines   maxEnergy = max(maxEnergy, line.initialElectronEnergy)   end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
    # Calculate all amplitudes and requested properties
    if  Distributed.nprocs() > 1
        ## NOT given the two keywords below: pmap would have to SERIALIZE both objects to a worker process for
        ## every line, and rebuilding the basis there costs only ~0.03 s. Omitting them is the documented
        ## fallback, so this branch simply behaves as before. Untested against a measurement -- with one
        ## process the threaded branch is what runs.
        newLines = Distributed.pmap(line -> ImpactExcitation.computeAmplitudesProperties(line, nm, grid, nrContinuum, settings), lines)
    else
        newLines = Vector{ImpactExcitation.Line}(undef, length(lines))
        ## Constants of the whole computation, built once and SHARED across the threads below; both are
        ## read-only (nothing writes into primitives.bsplinesL/bsplinesS), so sharing introduces no race.
        ## Display-only and the same for every line; set ONCE here rather than from inside the threaded loop.
        Defaults.setStandardSubshellList(Basics.generate(OrderedSubshellList(), finalMultiplet.levels[1].basis,
                                                         initialMultiplet.levels[1].basis); printout=false)
        nuclearPot = Nuclear.nuclearPotential(nm, grid)
        primitives = Bsplines.generatePrimitives(grid)
        Threads.@threads for l in eachindex(lines)
            newLines[l] = ImpactExcitation.computeAmplitudesProperties(lines[l], nm, grid, nrContinuum, settings;
                                                                       nuclearPot=nuclearPot, primitives=primitives)
        end
    end
    # Print all results to screen
    ImpactExcitation.displayResults(newLines)
    
    if settings.calcRateCoefficient 
        rates = ImpactExcitation.computeRateCoefficients(newLines, settings)
        ImpactExcitation.displayResults(rates)
    end
    #
    if    ( output  &&  settings.calcRateCoefficient) return( newLines, rates );
    elseif output   return( newLines );
    else            return( nothing );
    end
end


"""
`ImpactExcitation.determineChannels(finalLevel::Level, initialLevel::Level, settings::ImpactExcitation.Settings)`  
    ... to determine a list of electron-impact excitation Channels for a transitions from the initial to the final level and by 
        taking into account the particular settings of for this computation; an Array{ImpactExcitation.Channel,1} is returned.
"""
function determineChannels(finalLevel::Level, initialLevel::Level, settings::ImpactExcitation.Settings)
    channels = ImpactExcitation.Channel[];   
    symi = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
    #
    for  inKappa = -settings.maxKappa:settings.maxKappa
        if  inKappa == 0    continue    end
        symtList = AngularMomentum.allowedTotalSymmetries(symi, inKappa)
        for  symt in symtList
            outKappaList = AngularMomentum.allowedKappaSymmetries(symt, symf)
            for  outKappa in outKappaList
                push!(channels, ImpactExcitation.Channel( inKappa, outKappa, symt, 0., 0.,Complex(0.)) )
            end
        end
    end
    return( channels )  
end


"""
`ImpactExcitation.computeLinesCascade(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model,
                                      grid::Radial.Grid, settings::ImpactExcitation.Settings;
                                      output::Bool=true, printout::Bool=true)`
    ... to compute the electron-impact excitation amplitudes and cross sections as requested by the given settings.
        The computations and printout are adapted for larger cascade computations by including only lines with at
        least one channel and by reducing the printout. A list of lines::Array{ImpactExcitation.Line,1} is returned.

        This is the cascade counterpart of `ImpactExcitation.computeLines` and follows the same pattern as
        `PhotoExcitation.computeLinesCascade` and `AutoIonization.computeLinesCascade`. It was CALLED by
        Cascade.computeSteps(::ImpactExcitationScheme, ...) but had never been defined, so every electron-impact
        excitation cascade raised an UndefVarError at its first step.
"""
function  computeLinesCascade(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model,
                              grid::Radial.Grid, settings::ImpactExcitation.Settings;
                              output::Bool=true, printout::Bool=true)
    lines = ImpactExcitation.determineLines(finalMultiplet, initialMultiplet, settings)
    if  length(lines) == 0    return( ImpactExcitation.Line[] )    end
    ## The grid must be able to carry the fastest continuum orbital that occurs; see the note in computeLines().
    maxEnergy   = 0.;   for  line in lines   maxEnergy = max(maxEnergy, line.initialElectronEnergy)   end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
    #
    ## The per-line work is independent -- each line builds its own continuum orbital and its own amplitudes --
    ## so it is spread over the available Julia threads.  This mirrors ImpactExcitation.computeLines() above,
    ## which already threads the very same computeAmplitudesProperties(); the cascade variant was serial and
    ## thereby left the expensive part of every cascade on a single thread.  Start Julia with `julia -t N` to
    ## use it; with the default of one thread the loop below behaves exactly as the serial one did.
    ##
    ## NB the results are written BY INDEX into a pre-allocated vector rather than push!-ed: push! is not
    ## thread-safe, and indexing also preserves the order of the lines, which the rate-coefficient machinery
    ## relies on (ImpactExcitation.groupLines reshapes the list into an energies x transitions block).
    ## Lines without a single channel are dropped afterwards rather than skipped inside the loop.
    tmpLines = Vector{Union{Nothing, ImpactExcitation.Line}}(nothing, length(lines))
    doPrint  = printout  &&  Threads.nthreads() == 1     ## interleaved per-line printout would be unreadable
    ## Constants of the whole computation, built once and SHARED across the threads below; see computeLines.
    ## Display-only and the same for every line; set ONCE here rather than from inside the threaded loop.
    Defaults.setStandardSubshellList(Basics.generate(OrderedSubshellList(), finalMultiplet.levels[1].basis,
                                                     initialMultiplet.levels[1].basis); printout=false)
    nuclearPot = Nuclear.nuclearPotential(nm, grid)
    primitives = Bsplines.generatePrimitives(grid)
    Threads.@threads for  i  in  eachindex(lines)
        newLine = ImpactExcitation.computeAmplitudesProperties(lines[i], nm, grid, nrContinuum, settings; printout=doPrint,
                                                               nuclearPot=nuclearPot, primitives=primitives)
        if  length(newLine.channels) > 0    tmpLines[i] = newLine    end
    end
    newLines = ImpactExcitation.Line[ l  for l in tmpLines  if l !== nothing ]
    # Print the results; ImpactExcitation.displayResults() writes to stdout only
    if  printout   ImpactExcitation.displayResults(newLines)    end
    ## No rate coefficients here: a cascade computation returns LINES.  Effective collision strengths and plasma
    ## rate coefficients are formed from those lines by a cascade simulation (Cascade.EieRateCoefficients).
    #
    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`ImpactExcitation.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::ImpactExcitation.Settings)`
    ... to determine a list of ImpactExcitation.Line's for transitions between levels from the initial- and final-state multiplets,
        and by taking into account the particular selections and settings for this computation; an Array{ImpactExcitation.Line,1} is 
        returned. Apart from the level specification, all physical properties are set to zero during the initialization process.
"""
function  determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::ImpactExcitation.Settings)
    energyShift = Defaults.convertUnits("energy: to atomic", settings.energyShift)
    lines = ImpactExcitation.Line[]
    
    # The rate-coefficient path builds its OWN logarithmic energy grid below, from the threshold of each pair
    # out to maxEnergyMultiplier * dE, because the rate integral needs that range and that spacing.  Until
    # 28-Aug-2026 it did so while silently ignoring settings.electronEnergies: a user who asked for [500., 900.]
    # eV got 945.2, 937.2, 937.2 and no word about it, and the run then failed further downstream.  Rather than
    # quietly overriding a setting the caller made deliberately, the combination is refused.
    if  settings.calcRateCoefficient  &&  length(settings.electronEnergies) > 0
        error("\n\nImpactExcitation: calcRateCoefficient = true and electronEnergies = " *
              "$(settings.electronEnergies) cannot both be given.\n"                                            *
              "The rate coefficient is integrated over an energy grid that this module builds itself, from the\n" *
              "threshold of each level pair out to maxEnergyMultiplier times it, with numElectronEnergies\n"      *
              "points; your energies would be discarded.\n\n"                                                    *
              "    For rate coefficients:  leave electronEnergies empty and set maxEnergyMultiplier and\n"        *
              "                            numElectronEnergies instead.\n"                                       *
              "    For cross sections at YOUR energies:  set calcRateCoefficient = false.\n")
    end
    if settings.calcRateCoefficient
        for  iLevel  in  initialMultiplet.levels
            for  fLevel  in  finalMultiplet.levels
                if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                    ΔE = fLevel.energy - iLevel.energy
                    #enGrid = Radial.GridGL("Finite", ΔE, ΔE * settings.maxEnergyMultiplier, settings.numElectronEnergies, printout=false)
                    electronEnergies = exp.( LinRange(log(ΔE+0.0003), log(max(5., settings.maxEnergyMultiplier*ΔE)), 
                                                                                                settings.numElectronEnergies) ) 
                    for  initialElectronEnergy in electronEnergies         # initial electron energies
                        finalElectronEnergy    = initialElectronEnergy - (fLevel.energy - iLevel.energy) + energyShift
                        if  finalElectronEnergy < 0    continue   end  
                        #channels = ImpactExcitation.determineChannels(fLevel, iLevel, settings) 
                        channels = ImpactExcitation.Channel[]
                        # order channels as per Kappa to check convergence
                        push!( lines, ImpactExcitation.Line(iLevel, fLevel, initialElectronEnergy, finalElectronEnergy, 0., 0.,
                                                                                                                    channels, 1.0) )
                    end
                end
            end
        end
    else
        for  iLevel  in  initialMultiplet.levels
            for  fLevel  in  finalMultiplet.levels
                if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                    for  en in settings.electronEnergies
                        initialElectronEnergy  = Defaults.convertUnits("energy: to atomic", en)
                        finalElectronEnergy    = initialElectronEnergy - (fLevel.energy - iLevel.energy) + energyShift
                        if  finalElectronEnergy < 0    continue   end  
                        #channels = ImpactExcitation.determineChannels(fLevel, iLevel, settings) 
                        channels = ImpactExcitation.Channel[]
                        # order channels as per Kappa to check convergence
                        push!( lines, ImpactExcitation.Line(iLevel, fLevel, initialElectronEnergy, finalElectronEnergy, 0., 0.,
                                                                                                                    channels, 1.0) )
                    end
                end
            end
        end
    end

    return( lines )
end


"""
`ImpactExcitation.groupLines(lines::Array{ImpactExcitation.Line,1})`
    ... groups lines having the same initial and final level but different energies
        returns an Array{Array{ImpactExcitation.Line,1},1}
"""
function groupLines(lines::Array{ImpactExcitation.Line,1})
    indexPairs = Set{Tuple{Level, Level}}()

    for line in lines
        push!(indexPairs,(line.initialLevel, line.finalLevel))
    end

    groupLines = Array{ImpactExcitation.Line,1}[]

    for indexPair in indexPairs
        gl = ImpactExcitation.Line[]
        for line in lines
            if indexPair[1] == line.initialLevel && indexPair[2] == line.finalLevel
                push!(gl, line)
            end
        end
        push!(groupLines, gl)
    end

    return( groupLines )
end


"""
`ImpactExcitation.groupLines(lines::Array{ImpactExcitation.Line,1}, settings::ImpactExcitation.Settings)`
    ... groups lines having the same initial and final level but different energies
        returns an Array{Array{ImpactExcitation.Line,1},1}
"""
function groupLines(lines::Array{ImpactExcitation.Line,1}, settings::ImpactExcitation.Settings)

    n = settings.numElectronEnergies
    q, r = divrem(length(lines), n)
    if  r != 0
        error("ImpactExcitation.groupLines(): the $(length(lines)) lines do not divide evenly into groups of " *
              "settings.numElectronEnergies = $n.  This list has to be a rectangular (energies x transitions) " *
              "block; set numElectronEnergies to the number of DISTINCT incident electron energies actually " *
              "computed.")
    end

    groupLines = Array{ImpactExcitation.Line,1}[]

    lines = reshape(lines, n, q)

    for i = 1:q
        push!(groupLines, lines[:,i])
    end

    return( groupLines )
end


"""
`ImpactExcitation.interpolateCS(x::Float64, xa::Vector{Float64}, ya::Vector{Float64})`
    ... interpolates the cross section or collision stregths for a given initial electron Energy
        returns a Float64.
"""
function interpolateCS(x::Float64, xa::Vector{Float64}, ya::Vector{Float64})
    n = length(xa)

    # Alloc and setup
    obj = GSL.interp_alloc(gsl_interp_linear, n)
    #obj = GSL.interp_alloc(gsl_interp_cspline, n)
    acc = GSL.interp_accel_alloc()
    GSL.interp_init(obj, xa, ya, n)

    # Interpolate
    y = GSL.interp_eval(obj, xa, ya, x, acc)
    # Free
    GSL.interp_accel_free(acc)
    GSL.interp_free(obj)

    return( y )
end    


"""
`ImpactExcitation.interpolateCS(x::Float64, xa::Vector{Float64}, ya::Vector{Float64}, isE1Allowed::Bool)`
    ... interpolates the cross section or collision stregths for a given initial electron Energy
        if the energy is beyond the upper bound then uses asymptotic approximation for extrapolation,
        depending on the transition is E1 allowed or not.
        Returns a Float64.
"""
function interpolateCS(x::Float64, xa::Vector{Float64}, ya::Vector{Float64}, isE1Allowed::Bool)
    ## Guard the LOWER bound as well.  Only the upper one was handled, but the Gauss-Laguerre nodes of the
    ## thermal average are xt = t*kT with t starting near 0.137, so they reach far below the sampled energies
    ## whenever the energy grid does not begin at threshold -- and GSL then aborts with a bare
    ## "interpolation error at interp.c:150".  That never showed up as long as calcRateCoefficient = true
    ## auto-generated its own grid starting AT threshold (see ImpactExcitation.determineLines), but a cascade
    ## supplies its own energies and need not start there.  Below the sampled range the collision strength is
    ## simply not known, so a constant extrapolation is used -- the same conservative choice already made above
    ## the range for non-E1 transitions.  NB this makes rate coefficients at kT well below the lowest sampled
    ## energy unreliable; sample from near threshold if the low-temperature end matters.
    if      x < minimum(xa)     return( ya[1] )
    elseif  x <= maximum(xa)
        n = length(xa)
        # Alloc and setup
        obj = GSL.interp_alloc(gsl_interp_linear, n)
        #obj = GSL.interp_alloc(gsl_interp_cspline, n)
        acc = GSL.interp_accel_alloc()
        GSL.interp_init(obj, xa, ya, n)
        # Interpolate
        y = GSL.interp_eval(obj, xa, ya, x, acc)
        # Free
        GSL.interp_accel_free(acc)
        GSL.interp_free(obj)
    else  
        if isE1Allowed
            a = ya[end] / log(xa[end])
            y = a * log(x)
        else
            y = ya[end]
        end
    end
    return y
end 


"""
`ImpactExcitation.computeEffStrengths(lines::Array{ImpactExcitation.Line, 1}, settings::ImpactExcitation.Settings)`
    ... computes Effective collision strengths from the calculated line collision strengths at temperature(s) [K].
        Returns an Array{ImpactExcitation.RateCoefficients,1}.
"""
function computeEffStrengths(lines::Array{ImpactExcitation.Line, 1}, settings::ImpactExcitation.Settings)
    gLines = ImpactExcitation.groupLines(lines, settings)

    allEffOmegas = ImpactExcitation.RateCoefficients[]

    for lines in gLines

        symi = LevelSymmetry(lines[1].initialLevel.J, lines[1].initialLevel.parity);    
        symf = LevelSymmetry(lines[1].finalLevel.J, lines[1].finalLevel.parity)
        isE1Allowed = AngularMomentum.isAllowedMultipole(symi, EmMultipole("E1"), symf)

        # Interplating collision strengths for GL Quadrature
        energies = [ line.finalElectronEnergy for line in lines]
        collisionStrengths = [ line.collisionStrength for line in lines]
        
        enGrid = Radial.GridGL(Radial.GridGaussLegendreFinite(), minimum(energies), maximum(energies) , 25, printout=false)
        cs     = [ ImpactExcitation.interpolateCS(en, energies, collisionStrengths) for en in enGrid.t ]

        # For low temperature Gauss-Laguerre integration
        t, wt = FastGaussQuadrature.gausslaguerre(10)

        effOmegas = Float64[]
        for temp in settings.temperatures
            temp_au  = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", temp)
            wa = 0.0
            xt = t .* temp_au
            for i in eachindex(t)
                wa += ImpactExcitation.interpolateCS(xt[i], energies, collisionStrengths, isE1Allowed) * wt[i]
            end

            push!(effOmegas, wa)
        end
        push!(allEffOmegas, ImpactExcitation.RateCoefficients(lines[1].initialLevel, lines[1].finalLevel, settings.temperatures, 
                                                                    zeros(length(settings.temperatures)), effOmegas ) )
    end
    
    return( allEffOmegas )
end


"""
`ImpactExcitation.computeRateCoefficients(lines::Array{ImpactExcitation.Line,1}, settings::ImpactExcitation.Settings)`
    ... computes Exitation rate coefficients from the calculated collsion strengths at a temperature(s) [K].
        The rate coefficients are returned in [cm^3/s]. Returns an Array{ImpactExcitation.RateCoefficients,1}.
"""
function computeRateCoefficients(lines::Array{ImpactExcitation.Line, 1}, settings::ImpactExcitation.Settings)

    effStrengths = ImpactExcitation.computeEffStrengths(lines, settings)
    allAlphas = ImpactExcitation.RateCoefficients[]
    
    for effOmegas in effStrengths

        initialLevel = effOmegas.initialLevel ; finalLevel = effOmegas.finalLevel
        ΔE = finalLevel.energy - initialLevel.energy
        gi = 2.0 * ( initialLevel.J.num / initialLevel.J.den )  + 1

        alphas = Float64[]
        for (i, temp) in enumerate(effOmegas.temperatures)
            temp_au  = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", temp)
            wa = 8.6291269e-6 * effOmegas.effOmegas[i] * exp(-ΔE / temp_au) / gi / sqrt(temp)
            push!(alphas, wa)
        end
        push!(allAlphas, RateCoefficients( initialLevel, finalLevel, effOmegas.temperatures, alphas, effOmegas.effOmegas) )
    end

    return( allAlphas )
end


"""
`ImpactExcitation.displayLines(lines::Array{ImpactExcitation.Line,1})`  
    ... to display a list of lines and channels that have been selected due to the prior settings. A neat table of all 
        selected transitions and energies is printed but nothing is returned otherwise.
"""
function  displayLines(lines::Array{ImpactExcitation.Line,1})
    nx = 180
    println(" ")
    println("  Selected electron-impact ionization lines:")
    println(" ")
    println("  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=2);                                sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=3);                                sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(12, "f--Energy--i"; na=4)               
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=5)
    sa = sa * TableStrings.center(12, "Energy e_in"; na=3);              
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=3)
    sa = sa * TableStrings.center(12, "Energy e_out"; na=4);              
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.flushleft(57, "List of partial waves and total symmetries"; na=4)  
    sb = sb * TableStrings.flushleft(57, "partial-in [total J^P] partial-out        "; na=4)
    println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
    #
    NoLines = 0   
    for  line in lines
        NoLines = NoLines + 1
        sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        en = line.finalLevel.energy - line.initialLevel.energy
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", en))                          * "    "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.initialElectronEnergy))  * "    "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.finalElectronEnergy))    * "    "
        kappaInOutSymmetryList = Tuple{Int64,Int64,LevelSymmetry}[]
        for  i in 1:length(line.channels)
            push!( kappaInOutSymmetryList, (line.channels[i].initialKappa, line.channels[i].finalKappa, line.channels[i].symmetry) ) 
        end
        wa = TableStrings.kappaKappaSymmetryTupels(95, kappaInOutSymmetryList)
        # printBefore = true calls this BEFORE the channels exist, so the tuple list is empty and so is wa.
        # Until 28-Aug-2026 the `wa[1]` below then raised a BoundsError and killed the run on a documented
        # setting.  An empty list is the normal state here, not an error: the point of printing before the
        # computation is to show WHICH lines were selected, and the partial waves are not yet known.
        if  length(wa) == 0    println( sa );    continue    end
        sb = sa * wa[1];    println( sb )  
        for  i = 2:length(wa)
            NoLines = NoLines + 1
            sb = TableStrings.hBlank( length(sa) ) * wa[i];    println( sb )
        end
        # Avoid long table printouts
        if NoLines > 100  println("\n  ImpactExcitation.displayLines():  A maximum of 100 lines are printed in this table. \n")
            break   
        end
    end
    println("  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`ImpactExcitation.displayResults(lines::Array{ImpactExcitation.Line,1})`  
    ... to list all results, energies, cross sections, etc. of the selected lines. A neat table is printed but nothing is 
        returned otherwise.
"""
function  displayResults(lines::Array{ImpactExcitation.Line,1})
    nx = 142
    println(" ")
    println("  Electron-impact excitation cross sections:")
    println(" ")
    println("  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=2);                                sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=3);                                sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(12, "f--Energy--i"; na=4)               
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=5)
    sa = sa * TableStrings.center(12, "Energy e_in"; na=3);              
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=3)
    sa = sa * TableStrings.center(12, "Energy e_out"; na=3);              
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=3)
    sa = sa * TableStrings.center(15, "Cross section"; na=3)      
    #sb = sb * TableStrings.center(15, TableStrings.inUnits("cross section")*"*1E-8"; na=3)
    sb = sb * TableStrings.center(15, "10^(-20) m^2"; na=3)
    sa = sa * TableStrings.center(15, "Collision strength"; na=3)      
    sb = sb * TableStrings.center(15, " ";                  na=3)
    sa = sa * TableStrings.center(15, "Convergence"; na=3)      
    sb = sb * TableStrings.center(15, " ";                  na=3)
    println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
    #
    for  line in lines
        sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        en = line.finalLevel.energy - line.initialLevel.energy
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", en))                          * "    "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.initialElectronEnergy))  * "    "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.finalElectronEnergy))    * "     "
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", line.crossSection)*1e-8)    * "        "
        sa = sa * @sprintf("%.6e", line.collisionStrength)                                                    * "          "
        sa = sa * @sprintf("%.2e", line.convergence)                                                          * "    "
        if  line.convergence > convergenceCriterion   sa = sa * "<== NOT CONVERGED"   end
        println(sa)
    end
    println("  ", TableStrings.hLine(nx))
    #
    ## Say so explicitly if the partial-wave sum was still growing when it ran into settings.maxKappa. Unlike a
    ## too-coarse grid -- which Continuum.gridConsistency() now refuses outright -- this failure is silent: the
    ## collision strength simply comes out too SMALL, while remaining smooth and entirely plausible-looking.
    nNotConverged = count(line -> line.convergence > convergenceCriterion, lines)
    if  nNotConverged > 0
        println("")
        printstyled(">>> WARNING: $nNotConverged of $(length(lines)) lines did NOT reach the partial-wave " *
                    "convergence criterion ($convergenceCriterion).\n", color=:light_red)
        println("    The `convergence` column is the relative change of sum |amplitude|^2 contributed by the last ")
        println("    partial wave. A value above the criterion means the sum over kappa was still growing when it ")
        println("    hit settings.maxKappa, so the collision strengths of the flagged lines are too SMALL -- by an ")
        println("    amount this computation cannot bound. Nothing else signals this: the numbers stay smooth and ")
        println("    look entirely plausible. For a dipole-allowed transition the tell-tale symptom is that Omega ")
        println("    stops growing like ln(E) and flattens or even turns over at the upper end of the energy range.")
        println("    The requirement grows with impact energy, because the Bethe logarithm is built from large ")
        println("    impact parameters b ~ kappa/k, i.e. precisely from the highest partial waves.")
        printstyled("    >>> Increase settings.maxKappa and re-run.\n", color=:light_red)
        println("")
    end
    #
    return( nothing )
end


"""
`ImpactExcitation.displayResults(line::ImpactExcitation.Line)`
    ... to display energies, cross sections, collision stregths, convergence etc. for a single line immediatly after
        the calculation for that line is completed.
"""
function displayResults(line::ImpactExcitation.Line)
    nx = 142
    println(" ")
    println("  Electron-impact excitation cross sections:")
    println(" ")
    println("  ", TableStrings.hLine(nx))
    sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
    fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
    sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
    sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
    en = line.finalLevel.energy - line.initialLevel.energy
    sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", en))                          * "    "
    sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.initialElectronEnergy))  * "    "
    sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.finalElectronEnergy))    * "     "
    sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", line.crossSection)*1e-8)    * "        "
    sa = sa * @sprintf("%.6e", line.collisionStrength)                                                    * "          "
    sa = sa * @sprintf("%.2e", line.convergence)                                                          * "    "
    println(sa)
    println("  ", TableStrings.hLine(nx))
    return( nothing )
end


"""
`ImpactExcitation.displayResults(allRates::Array{RateCoefficients,1})`  
    ... to list the excitation rate coefficients and effective collision strengths for the selected lines at the selected temperatures.
        A neat table is printed but nothing is returned otherwise.
"""
function displayResults(allRates::Array{RateCoefficients,1})
    ImpactExcitation.displayResults(stdout, allRates)
end


"""
`ImpactExcitation.displayResults(stream::IO, allRates::Array{ImpactExcitation.RateCoefficients,1})`
    ... the same table, but written to the given stream so that a cascade can also place it into the summary
        file.  A rate coefficient is the quantity such a computation is run for and should not go to the
        screen alone.
"""
function displayResults(stream::IO, allRates::Array{RateCoefficients,1})
    nx = 115
    println(stream, " ")
    println(stream, "  Electron-impact excitation rate coefficients / effective collision strengths:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "      ";   sb = "  "
    sa = sa * TableStrings.center(14, "i-level-f"; na=2);                                sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=3);                                sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(12, "f--Energy--i"; na=4)               
    sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(14, "Temperature"; na=2) 
    sb = sb * TableStrings.center(12, "[K]"; na=5)
    sa = sa * TableStrings.center(12, "Rate Coeffs"; na=2) 
    sb = sb * TableStrings.center(12, "[cm^3/s]"; na=5)
    sa = sa * TableStrings.center(15, "Eff Coll Strength"; na=4) 
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))

    for rates in allRates
        for (i, temp) in enumerate(rates.temperatures)
            sa  = "  ";    isym = LevelSymmetry( rates.initialLevel.J, rates.initialLevel.parity)
                            fsym = LevelSymmetry( rates.finalLevel.J,   rates.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(rates.initialLevel.index, rates.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
            en = rates.finalLevel.energy - rates.initialLevel.energy
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", en))                          * "      "
            sa = sa * @sprintf("%.4e", temp)         * "      "
            sa = sa * @sprintf("%.4e", rates.alphas[i])         * "       "
            sa = sa * @sprintf("%.4e", rates.effOmegas[i])         * "      "
            println(stream, sa)
        end
    end

    println(stream, "  ", TableStrings.hLine(nx))
end


end # module end
