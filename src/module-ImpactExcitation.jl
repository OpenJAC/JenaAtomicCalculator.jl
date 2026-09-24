
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
    + continuumExchange        ::Union{Nothing,Bool,Symbol}
        ... True/false to force exchange between the impacting/scattered electron and the bound electrons on or off;
            `nothing` -- the default -- lets `Basics.perform` decide automatically from the final-ion `AsfSettings.scField`
            (`Continuum.resolveExchangeDefault`, item 99).  A Symbol instead selects the exchange model of the continuum
            potential itself, and is then used as given: :Slater, :FurnessMcCarthy or :none (Continuum.continuumPotential).
    + calcDCS                  ::Bool               ... True, if differential cross sections are to be calculated.
    + angles                   ::Array{Float64,1}   ... Angles in degrees for calculating differential cross sections.
    + cacheAngularCoefficients ::Bool               ... Cache SpinAngular coefficient tables across energies in this run.
    + cacheRadialIntegrals     ::Bool               ... Cache radial integrals within each individual line/energy.
    + orthogonalizeContinuum   ::Bool               ... Orthogonalize each continuum orbital to the bound orbitals of the same kappa
                                                        (Continuum.orthogonalizeToBasis), as the spin-angular coefficients assume.
    + asymptoticTail           ::Bool               ... Add the part of the direct Coulomb integrals beyond the radial box
                                                        (ImpactExcitation.tailCorrection).
    + topUp                    ::Bool               ... Add the partial waves beyond maxKappa by geometric extrapolation of the
                                                        last two kappa shells (ImpactExcitation.topUpCollisionStrength).
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
    continuumExchange         ::Union{Nothing,Bool,Symbol}
    calcDCS                   ::Bool
    angles                    ::Array{Float64,1}
    cacheAngularCoefficients  ::Bool
    cacheRadialIntegrals      ::Bool
    orthogonalizeContinuum    ::Bool
    asymptoticTail            ::Bool
    topUp                     ::Bool
end


"""
`ImpactExcitation.Settings()`  ... constructor for the default values of electron-impact excitation line computations.
"""
function Settings()
    Settings( LineSelection(), Float64[], 0., 300, false, 30.0, 6,
              [5e4, 1e5, 2e5, 3e5, 5e5, 7e5, 1e6], false, CoulombInteraction(), nothing,
              false, Float64[], true, true, false, false, false)
end


"""
`ImpactExcitation.continuumExchangeModel(settings::ImpactExcitation.Settings)`
    ... returns the exchange model of the continuum potential selected by settings.continuumExchange, which is :Slater unless that
        field names another one; a model::Symbol is returned.
"""
function continuumExchangeModel(settings::ImpactExcitation.Settings)
    return( settings.continuumExchange isa Symbol ? settings.continuumExchange : :Slater )
end


"""
`ImpactExcitation.Settings(set::ImpactExcitation.Settings;`

    lineSelection...,  electronEnergies..., energyShift..., maxKappa...,
    calcRateCoefficient..., maxEnergyMultiplier..., numElectronEnergies..., temperatures...,
    printBefore..., operator..., continuumExchange..., calcDCS..., angles..., cacheAngularCoefficients..., cacheRadialIntegrals...,
    orthogonalizeContinuum..., asymptoticTail..., topUp...)

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
    operator::Union{Nothing, AbstractEeInteraction} = nothing,
    continuumExchange::Union{Nothing, Bool, Symbol} = nothing,
    calcDCS::Union{Nothing, Bool} = nothing,
    angles::Union{Nothing, Array{Float64,1}} = nothing,
    cacheAngularCoefficients::Union{Nothing, Bool} = nothing,
    cacheRadialIntegrals::Union{Nothing, Bool} = nothing,
    orthogonalizeContinuum::Union{Nothing, Bool} = nothing,
    asymptoticTail::Union{Nothing, Bool} = nothing,
    topUp::Union{Nothing, Bool} = nothing )

    if  isnothing(lineSelection)            lineSelectionx            = set.lineSelection            else  lineSelectionx            = lineSelection            end
    if  isnothing(electronEnergies)         electronEnergiesx         = set.electronEnergies         else  electronEnergiesx         = electronEnergies         end
    if  isnothing(energyShift)              energyShiftx              = set.energyShift              else  energyShiftx              = energyShift              end
    if  isnothing(maxKappa)                 maxKappax                 = set.maxKappa                 else  maxKappax                 = maxKappa                 end
    if  isnothing(calcRateCoefficient)      calcRateCoefficientx      = set.calcRateCoefficient      else  calcRateCoefficientx      = calcRateCoefficient      end
    if  isnothing(maxEnergyMultiplier)      maxEnergyMultiplierx      = set.maxEnergyMultiplier      else  maxEnergyMultiplierx      = maxEnergyMultiplier      end
    if  isnothing(numElectronEnergies)      numElectronEnergiesx      = set.numElectronEnergies      else  numElectronEnergiesx      = numElectronEnergies      end
    if  isnothing(temperatures)             temperaturesx             = set.temperatures             else  temperaturesx             = temperatures             end
    if  isnothing(printBefore)              printBeforex              = set.printBefore              else  printBeforex              = printBefore              end
    if  isnothing(operator)                 operatorx                 = set.operator                 else  operatorx                 = operator                 end
    if  isnothing(continuumExchange)        continuumExchangex        = set.continuumExchange        else  continuumExchangex        = continuumExchange        end
    if  isnothing(calcDCS)                  calcDCSx                  = set.calcDCS                  else  calcDCSx                  = calcDCS                  end
    if  isnothing(angles)                   anglesx                   = set.angles                   else  anglesx                   = angles                   end
    if  isnothing(cacheAngularCoefficients) cacheAngularCoefficientsx = set.cacheAngularCoefficients else  cacheAngularCoefficientsx = cacheAngularCoefficients end
    if  isnothing(cacheRadialIntegrals)     cacheRadialIntegralsx     = set.cacheRadialIntegrals     else  cacheRadialIntegralsx     = cacheRadialIntegrals     end
    if  isnothing(orthogonalizeContinuum)   orthogonalizeContinuumx   = set.orthogonalizeContinuum   else  orthogonalizeContinuumx   = orthogonalizeContinuum   end
    if  isnothing(asymptoticTail)           asymptoticTailx           = set.asymptoticTail           else  asymptoticTailx           = asymptoticTail           end
    if  isnothing(topUp)                    topUpx                    = set.topUp                    else  topUpx                    = topUp                    end

    Settings( lineSelectionx, electronEnergiesx, energyShiftx, maxKappax, calcRateCoefficientx, maxEnergyMultiplierx,
              numElectronEnergiesx, temperaturesx, printBeforex, operatorx, continuumExchangex, calcDCSx, anglesx,
              cacheAngularCoefficientsx, cacheRadialIntegralsx, orthogonalizeContinuumx, asymptoticTailx,
              topUpx )

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
    println(io, "continuumExchange:          $(settings.continuumExchange)  ")
    println(io, "calcDCS:                    $(settings.calcDCS)  ")
    println(io, "angles:                     $(settings.angles)  ")
    println(io, "cacheAngularCoefficients:   $(settings.cacheAngularCoefficients)  ")
    println(io, "cacheRadialIntegrals:       $(settings.cacheRadialIntegrals)  ")
    println(io, "orthogonalizeContinuum:     $(settings.orthogonalizeContinuum)  ")
    println(io, "asymptoticTail:             $(settings.asymptoticTail)  ")
    println(io, "topUp:                      $(settings.topUp)  ")
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
    + angles                 ::Array{Float64,1}   ... angles for differential cross section calculation.
    + dcs                    ::Array{Float64,1}   ... differential cross sections for these angles.
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
    angles                   ::Array{Float64,1}
    dcs                      ::Array{Float64,1}
    channels                 ::Array{ImpactExcitation.Channel,1}
    convergence              ::Float64
end 


"""
`ImpactExcitation.Line()`  ... 'empty' constructor for an electron-impact excitation line between a specified initial and final level.
"""
function Line()
    Line(Level(), Level(), 0., 0., 0., 0., Float64[], Float64[], ImpactExcitation.Channel[], 1.0 )
end


"""
`ImpactExcitation.Line(initialLevel::Level, finalLevel::Level, crossSection::Float64)`  
    ... constructor for an electron-impact excitation line between a specified initial and final level.
"""
function Line(initialLevel::Level, finalLevel::Level, crossSection::Float64)
    Line(initialLevel, finalLevel, 0., 0., crossSection, 0., Float64[], Float64[], ImpactExcitation.Channel[], 1.0 )
end


# `Base.show(io::IO, line::ImpactExcitation.Line)`  ... prepares a proper printout of the variable line::ImpactExcitation.Line.
function Base.show(io::IO, line::ImpactExcitation.Line) 
    println(io, "initialLevel:            $(line.initialLevel)  ")
    println(io, "finalLevel:              $(line.finalLevel)  ")
    println(io, "initialElectronEnergy:   $(line.initialElectronEnergy)  ")
    println(io, "finalElectronEnergy:     $(line.finalElectronEnergy)  ")
    println(io, "crossSection:            $(line.crossSection)  ")
    println(io, "collisionStrength:       $(line.collisionStrength)  ")
    println(io, "angles:                  $(line.angles)  ")
    println(io, "dcs:                     $(line.dcs)  ")
    println(io, "channels:                $(line.channels)  ")
end


"""
`struct  ImpactExcitation.AngularCoefficientCache`
    ... memoises, for ONE run of computeLines() or computeLinesCascade(), the table of spin-angular coefficients of every
        (final CSF, initial CSF) pair for a given line and channel. These coefficients depend on the CSF topology -- the two
        target levels, the two continuum kappas and the total symmetry -- but not on the electron energy, so a run over
        several impact energies reuses each table once per energy. The cache is created by the driver and goes out of
        scope with it, like InteractionStrength.XLCache, so its keys never outlive the multiplets they refer to. The lock
        makes it safe to share between the threads of one run.

    + tables    ::Dict{Tuple, Matrix{Any}}   ... the memoised tables; an entry `nothing` marks a CSF pair of wrong symmetry.
    + lock      ::ReentrantLock              ... guards reads and writes of tables.
"""
struct  AngularCoefficientCache
    tables      ::Dict{Tuple, Matrix{Any}}
    lock        ::ReentrantLock
end


"""
`ImpactExcitation.AngularCoefficientCache()`
    ... constructor for an empty cache; a cache::ImpactExcitation.AngularCoefficientCache is returned.
"""
AngularCoefficientCache() = AngularCoefficientCache( Dict{Tuple, Matrix{Any}}(), ReentrantLock() )


# SpinAngular keeps process-wide memo caches of its own that are not thread-safe (Rule 19), so every call into
# SpinAngular.computeCoefficients() from the threads of this module is serialized through this one lock.
const SPIN_ANGULAR_LOCK = ReentrantLock()


"""
`ImpactExcitation.angularCoefficientKey(channel::ImpactExcitation.Channel, fLevel::Level, iLevel::Level)`
    ... identifies an energy-independent table of spin-angular coefficients within one run: the two target levels, the two
        continuum kappas, the total symmetry and the common subshell list. A key::Tuple is returned.
"""
function angularCoefficientKey(channel::ImpactExcitation.Channel, fLevel::Level, iLevel::Level)
    Jt = channel.symmetry.J
    return( (iLevel.index, fLevel.index, channel.initialKappa, channel.finalKappa, Jt.num, Jt.den, string(channel.symmetry.parity),
             length(iLevel.basis.csfs), length(fLevel.basis.csfs), Tuple(fLevel.basis.subshells)) )
end


"""
`ImpactExcitation.computeAngularCoefficientTable(fLevel::Level, iLevel::Level)`
    ... computes the spin-angular coefficients of the scalar two-particle operator for every (final CSF, initial CSF) pair of
        the two (continuum) levels; the entry is `nothing` for an initial CSF whose symmetry differs from the level's. The call
        is serialized through SPIN_ANGULAR_LOCK. A table::Matrix{Any} is returned.
"""
function computeAngularCoefficientTable(fLevel::Level, iLevel::Level)
    nf = length(fLevel.basis.csfs);   ni = length(iLevel.basis.csfs)
    table = Matrix{Any}(undef, nf, ni)
    op2   = SpinAngular.TwoParticleOperator(0, plus)
    lock(SPIN_ANGULAR_LOCK) do
        for  r = 1:nf,  s = 1:ni
            if  iLevel.basis.csfs[s].J != iLevel.J  ||  iLevel.basis.csfs[s].parity != iLevel.parity
                table[r,s] = nothing
            else
                table[r,s] = SpinAngular.computeCoefficients(op2, fLevel.basis.csfs[r], iLevel.basis.csfs[s], fLevel.basis.subshells)
            end
        end
    end

    return( table )
end


"""
`ImpactExcitation.angularCoefficientTable(channel::ImpactExcitation.Channel, fLevel::Level, iLevel::Level,
                                           cache::Union{Nothing,ImpactExcitation.AngularCoefficientCache})`
    ... returns the table of spin-angular coefficients for the given channel, from the cache when it holds one and computing
        and storing it otherwise; with cache = nothing the table is always computed. A table::Matrix{Any} is returned.
"""
function angularCoefficientTable(channel::ImpactExcitation.Channel, fLevel::Level, iLevel::Level,
                                 cache::Union{Nothing,ImpactExcitation.AngularCoefficientCache})
    isnothing(cache)   &&   return( ImpactExcitation.computeAngularCoefficientTable(fLevel, iLevel) )
    key    = ImpactExcitation.angularCoefficientKey(channel, fLevel, iLevel)
    cached = lock(cache.lock) do;   get(cache.tables, key, nothing)   end
    isnothing(cached)  ||  return( cached )
    # Two threads may both miss and both compute the same table; they produce identical tables, so the second write is harmless.
    table  = ImpactExcitation.computeAngularCoefficientTable(fLevel, iLevel)
    lock(cache.lock) do;   cache.tables[key] = table   end

    return( table )
end


"""
`struct  ImpactExcitation.ContinuumTail`
    ... the analytic continuation of one continuum orbital beyond the last point r0 of its grid representation. Outside the
        bound orbitals the potential is the pure Coulomb tail -Z_ion/r, so the large component is  P(r) = A F_l(eta,qr) +
        B G_l(eta,qr)  with eta = -Z_ion (1 + E/c^2)/q; A and B are fitted to the grid orbital at r0 and a quarter wavelength
        inside it. The small component is not continued: in the tail it contributes at relative order (alpha q)^2.

    + r0          ::Float64   ... last grid point of the orbital, where the grid representation ends.
    + q           ::Float64   ... (relativistic) wave number.
    + eta         ::Float64   ... Sommerfeld parameter, negative for an attractive ion.
    + l           ::Int64     ... orbital angular momentum of the large component.
    + A           ::Float64   ... coefficient of the regular Coulomb function F_l.
    + B           ::Float64   ... coefficient of the irregular Coulomb function G_l.
    + fitError    ::Float64   ... |P_fit - P_grid| / sqrt(A^2+B^2) at an intermediate point; the quality of the continuation.
"""
struct  ContinuumTail
    r0            ::Float64
    q             ::Float64
    eta           ::Float64
    l             ::Int64
    A             ::Float64
    B             ::Float64
    fitError      ::Float64
end


"""
`struct  ImpactExcitation.TailCache`
    ... memoises, for ONE line (i.e. one energy), the continuations of its continuum orbitals, the bound multipole moments
        and the tail integrals, all of which recur across channels and total symmetries.

    + Zion        ::Float64                                           ... charge of the target ion seen by the free electron.
    + tails       ::Dict{Subshell, ImpactExcitation.ContinuumTail}    ... continuation of each continuum orbital.
    + moments     ::Dict{Tuple{Int64,Subshell,Subshell}, Float64}     ... bound multipole moments M_L(a,c).
    + integrals   ::Dict{Tuple{Int64,Subshell,Subshell}, Float64}     ... tail integrals int_r0^infinity P_b P_d r^-(L+1) dr.
"""
struct  TailCache
    Zion          ::Float64
    tails         ::Dict{Subshell, ImpactExcitation.ContinuumTail}
    moments       ::Dict{Tuple{Int64,Subshell,Subshell}, Float64}
    integrals     ::Dict{Tuple{Int64,Subshell,Subshell}, Float64}
end


"""
`ImpactExcitation.TailCache(Zion::Float64)`
    ... constructor for an empty tail cache of a line with target-ion charge Zion; a cache::ImpactExcitation.TailCache is returned.
"""
TailCache(Zion::Float64) = TailCache(Zion, Dict{Subshell, ContinuumTail}(), Dict{Tuple{Int64,Subshell,Subshell}, Float64}(),
                                     Dict{Tuple{Int64,Subshell,Subshell}, Float64}())


"""
`ImpactExcitation.coulombFG(t::ImpactExcitation.ContinuumTail, r::Float64)`
    ... evaluates the continued large component and its r-derivative at r >= t.r0. A tuple (P::Float64, dP::Float64) is returned.
"""
function coulombFG(t::ImpactExcitation.ContinuumTail, r::Float64)
    f, fp, g, gp = GSL.sf_coulomb_wave_FG_e(t.eta, t.q * r, Float64(t.l), 0, 0.0, 0.0)

    return( (t.A * f.val + t.B * g.val,  t.q * (t.A * fp.val + t.B * gp.val)) )
end


"""
`ImpactExcitation.continuumTail(orbital::Orbital, grid::Radial.Grid, Zion::Float64)`
    ... fits the Coulomb continuation of a continuum orbital beyond its last grid point; see ImpactExcitation.ContinuumTail.
        A tail::ImpactExcitation.ContinuumTail is returned.
"""
function continuumTail(orbital::Orbital, grid::Radial.Grid, Zion::Float64)
    wc  = Defaults.getDefaults("speed of light: c");    E = orbital.energy
    q   = sqrt(E * (E + 2wc^2)) / wc;                  eta = -Zion * (1.0 + E/wc^2) / q
    l   = Basics.subshell_l(orbital.subshell)
    m0  = length(orbital.P);   r0 = grid.r[m0]
    m1  = m0;   while  m1 > 2  &&  r0 - grid.r[m1] < 0.5pi/q    m1 = m1 - 1   end
    mh  = div(m0 + m1, 2)
    fg(r) = (f = GSL.sf_coulomb_wave_FG_e(eta, q*r, Float64(l), 0, 0.0, 0.0);   (f[1].val, f[3].val))
    F0, G0 = fg(r0);   F1, G1 = fg(grid.r[m1])
    det = F0*G1 - F1*G0
    A   = (orbital.P[m0]*G1 - orbital.P[m1]*G0) / det
    B   = (F0*orbital.P[m1] - F1*orbital.P[m0]) / det
    Fh, Gh = fg(grid.r[mh])
    err = abs(A*Fh + B*Gh - orbital.P[mh]) / sqrt(A^2 + B^2)

    return( ContinuumTail(r0, q, eta, l, A, B, err) )
end


"""
`ImpactExcitation.tailIntegral(tb::ImpactExcitation.ContinuumTail, td::ImpactExcitation.ContinuumTail, n::Int64)`
    ... computes  int_r0^infinity P_b(r) P_d(r) r^(-n) dr  for two continued continuum orbitals with a common r0: by Simpson
        on a uniform mesh of 40 points per shortest wavelength out to X, where both waves are past their centrifugal
        turning points and at least four times r0; then analytically from X to infinity, with each wave written as
        amp sin(theta(X) + theta'(X)(r-X)) so that the remainder is a pair of generalized exponential integrals E_n. This
        is the treatment of the asymptotic region of the Fortran RDW code (Purcell et al. 1987), with the Coulomb phase
        carried by theta' = q - eta/X. A value::Float64 is returned.
"""
function tailIntegral(tb::ImpactExcitation.ContinuumTail, td::ImpactExcitation.ContinuumTail, n::Int64)
    r0 = tb.r0
    X  = max( 4r0, 3.0 * (max(tb.l, td.l) + 1 + abs(tb.eta) + abs(td.eta)) / min(tb.q, td.q) )
    h  = 2pi / max(tb.q, td.q) / 40;    N = 2 * ceil(Int64, (X - r0) / h / 2);    h = (X - r0) / N
    wa = 0.
    for  i = 0:N
        r  = r0 + i*h;    w = (i == 0 || i == N) ? 1.0 : (isodd(i) ? 4.0 : 2.0)
        wa = wa + w * coulombFG(tb, r)[1] * coulombFG(td, r)[1] / r^n
    end
    wa = wa * h / 3

    function ampPhase(t::ImpactExcitation.ContinuumTail)
        u, du = coulombFG(t, X);    om = t.q - t.eta / X
        return( (sqrt(u^2 + (du/om)^2), atan(u, du/om), om) )
    end
    Ab, thb, omb = ampPhase(tb);    Ad, thd, omd = ampPhase(td)
    I(om) = abs(om) < 1.0e-12 ? X^(1-n) / (n-1) : exp(-im*om*X) * X^(1-n) * SpecialFunctions.expint(n, -im*om*X)
    rem = Ab * Ad / 2 * real( exp(im*(thb - thd)) * I(omb - omd) - exp(im*(thb + thd)) * I(omb + omd) )

    return( wa + rem )
end


"""
`ImpactExcitation.tailCorrection(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid,
                                  cache::ImpactExcitation.TailCache)`
    ... the part of the Coulomb strength X^L(abcd) that lies beyond the radial box, for a DIRECT term in which one of the
        pairs (a,c), (b,d) is bound-bound and the other continuum-continuum: outside the bound orbitals the potential of
        the bound transition density is exactly M_L / r^(L+1), so the missing piece is  xc * M_L * int_r0^inf P P' r^-(L+1),
        with the same angular prefactor xc as in InteractionStrength.XL_Coulomb. Zero for every other term -- exchange
        and bound-continuum densities vanish beyond the bound orbitals -- and for L = 0, where M_0 = <a|c> = 0 by
        orthogonality. A value::Float64 is returned.
"""
function tailCorrection(L::Int64, a::Orbital, b::Orbital, c::Orbital, d::Orbital, grid::Radial.Grid,
                        cache::ImpactExcitation.TailCache)
    L == 0   &&   return( 0. )
    if      a.isBound  &&  c.isBound  &&  !b.isBound  &&  !d.isBound     bnd1, bnd2, con1, con2 = a, c, b, d
    elseif  b.isBound  &&  d.isBound  &&  !a.isBound  &&  !c.isBound     bnd1, bnd2, con1, con2 = b, d, a, c
    else    return( 0. )
    end
    ja2 = Basics.subshell_2j(a.subshell);   jb2 = Basics.subshell_2j(b.subshell)
    jc2 = Basics.subshell_2j(c.subshell);   jd2 = Basics.subshell_2j(d.subshell)
    if  AngularMomentum.triangularDelta(ja2+1,jc2+1,L+L+1) * AngularMomentum.triangularDelta(jb2+1,jd2+1,L+L+1) == 0   ||
        rem(Basics.subshell_l(a.subshell) + Basics.subshell_l(c.subshell) + L, 2) == 1   ||
        rem(Basics.subshell_l(b.subshell) + Basics.subshell_l(d.subshell) + L, 2) == 1
        return( 0. )
    end
    xc = AngularMomentum.CL_reduced_me(a.subshell, L, c.subshell) * AngularMomentum.CL_reduced_me(b.subshell, L, d.subshell)
    if  rem(L,2) == 1    xc = - xc    end

    M  = get!(cache.moments, (L, bnd1.subshell, bnd2.subshell)) do
             n = min(length(bnd1.P), length(bnd2.P), grid.NoPoints)
             sum( (bnd1.P[i]*bnd2.P[i] + bnd1.Q[i]*bnd2.Q[i]) * grid.r[i]^L * grid.wr[i] for i = 1:n )
         end
    tb = get!(cache.tails, con1.subshell) do;   ImpactExcitation.continuumTail(con1, grid, cache.Zion)   end
    td = get!(cache.tails, con2.subshell) do;   ImpactExcitation.continuumTail(con2, grid, cache.Zion)   end
    abs(tb.r0 - td.r0) < 1.0e-10   ||   error("ImpactExcitation.tailCorrection: the two continua end at different radii.")
    T  = get!(cache.integrals, (L, con1.subshell, con2.subshell)) do
             ImpactExcitation.tailIntegral(tb, td, L+1)
         end

    return( xc * M * T )
end


"""
`ImpactExcitation.topUpCollisionStrength(shellColl::Array{Float64,1})`
    ... estimates the contribution of all kappa shells beyond the last one computed, from the (unnormalized) partial collision
        strengths of the kappa shells, by the geometric series  Omega_K r/(1-r)  with  r = Omega_K/Omega_(K-1). This is the
        usual top-up for NON-dipole transitions and a lower bound for dipole ones, whose partial strengths fall off more
        slowly than geometrically (the Bethe logarithm); a warning is issued when r > 0.9, where the estimate is least
        reliable. Zero is returned if fewer than two shells exist or the last two do not decrease. A value::Float64 is
        returned, in the units of shellColl.
"""
function topUpCollisionStrength(shellColl::Array{Float64,1})
    length(shellColl) < 2   &&   return( 0. )
    r = shellColl[end] / shellColl[end-1]
    if  !(0. < r < 1.)    return( 0. )    end
    if  r > 0.9
        @warn("ImpactExcitation.topUpCollisionStrength: the partial collision strengths fall off slowly (ratio $(round(r, digits=4)))" *
              " -- the geometric top-up is a lower bound here; raise maxKappa.")
    end

    return( shellColl[end] * r / (1. - r) )
end


"""
`ImpactExcitation.amplitude(kind::AbstractEeInteraction, channel::ImpactExcitation.Channel, cFinalLevel::Level, cInitialLevel::Level, 
                            grid::Radial.Grid; printout::Bool=true)`  
    ... to compute the kind in  CoulombInteraction(), BreitInteraction(), CoulombBreit() electron-impact interaction amplitude 
        <(alpha_f J_f, kappa_f) J_t || O^(e-e, kind) || (alpha_i J_i, kappa_i) J_t>  due to the interelectronic interaction for 
        the given final and initial (continuum) level. A value::ComplexF64 is returned.
"""
function amplitude(kind::AbstractEeInteraction, channel::ImpactExcitation.Channel, cFinalLevel::Level, cInitialLevel::Level,
                    grid::Radial.Grid; printout::Bool=true,
                    radialCache::Union{Nothing,InteractionStrength.XLCache}=nothing,
                    angularCache::Union{Nothing,ImpactExcitation.AngularCoefficientCache}=nothing,
                    tailCache::Union{Nothing,ImpactExcitation.TailCache}=nothing)
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
        angularTable = ImpactExcitation.angularCoefficientTable(channel, fLevel, iLevel, angularCache)
        for  r = 1:nf
            for  s = 1:ni
                wa = angularTable[r,s]
                wa === nothing && continue
                me = 0.

                for  coeff in wa
                    # The exchange T-matrices are considered up to 40 partial waves
                    if ( (coeff.c.n == 101) && (coeff.nu > 40) )  continue end

                    oa = fLevel.basis.orbitals[coeff.a];   ob = fLevel.basis.orbitals[coeff.b]
                    oc = iLevel.basis.orbitals[coeff.c];   od = iLevel.basis.orbitals[coeff.d]
                    if  typeof(kind) in [CoulombInteraction, CoulombBreit]
                        me = me + coeff.V * ( isnothing(radialCache) ? InteractionStrength.XL_Coulomb(coeff.nu, oa, ob, oc, od, grid) :
                                              InteractionStrength.XL_Coulomb(coeff.nu, oa, ob, oc, od, grid, radialCache) )
                        if  !isnothing(tailCache)
                            me = me + coeff.V * ImpactExcitation.tailCorrection(coeff.nu, oa, ob, oc, od, grid, tailCache)
                        end
                    end
                    if  typeof(kind) in [BreitInteraction, CoulombBreit]
                        me = me + coeff.V * ( isnothing(radialCache) ?
                                              InteractionStrength.XL_Breit(coeff.nu, oa, ob, oc, od, grid, CoulombBreit(1.0)) :
                                              InteractionStrength.XL_Breit(coeff.nu, oa, ob, oc, od, grid, CoulombBreit(1.0), radialCache) )
                    end
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
`ImpactExcitation.checkLowEnergyGrid(minEnergy::Float64, grid::Radial.Grid)`
    ... warns when the SLOWEST continuum electron of a run has a de Broglie wavelength too long to complete an oscillation
        inside the radial box, so that its normalization and phase cannot be read off the asymptotic region. Nothing is
        returned; the run continues.
"""
function checkLowEnergyGrid(minEnergy::Float64, grid::Radial.Grid)
    if  minEnergy <= 0.    return( nothing )    end
    wavenb   = sqrt( 2minEnergy + minEnergy * Defaults.getDefaults("alpha")^2 )
    wavelgth = 2pi / wavenb
    rbox = grid.r[end]
    if  wavelgth > rbox
        @warn("ImpactExcitation: the slowest continuum electron of this run does not fit inside the box.\n" *
              @sprintf("    lowest continuum energy  = %.6e Hartree  (%.4f eV)\n", minEnergy,
                       Defaults.convertUnits("energy: from atomic to eV", minEnergy))                        *
              @sprintf("    de Broglie wavelength    = %.4f a.u.\n", wavelgth)                              *
              @sprintf("    box radius grid.r[end]   = %.4f a.u.\n", rbox)                                  *
              @sprintf("    wavelength / box         = %.2f      (should be well below 1)\n", wavelgth/rbox) *
              "Normalization and phase are taken near the outer boundary, and an orbital that has not completed one\n"  *
              "oscillation by then has no asymptotic region to be read from. The cross sections at the lowest energies\n" *
              "are therefore unreliable by an amount this run cannot bound; higher energies are unaffected.\n"           *
              ">>> Enlarge the box, or raise the lowest energy (with calcRateCoefficient = true the grid starts AT\n"       *
              "    threshold and only a larger box will help).")
    end

    return( nothing )
end




"""
`ImpactExcitation.dcsYlmVector!(cache::Dict{Tuple{Int64,Int64},Vector{ComplexF64}}, l::Int64, m::Int64, angles::Array{Float64,1})`
    ... returns Y_lm(theta, phi=0) on the given angle grid (in degrees) from AngularMomentum.sphericalYlm(), memoised in the
        cache, which belongs to ONE angle grid. A ylm::Vector{ComplexF64} is returned.
"""
function dcsYlmVector!(cache::Dict{Tuple{Int64,Int64},Vector{ComplexF64}}, l::Int64, m::Int64, angles::Array{Float64,1})
    return get!(cache, (l,m)) do
        ComplexF64[ AngularMomentum.sphericalYlm(l, m, deg2rad(theta), 0.0) for theta in angles ]
    end
end


"""
`ImpactExcitation.accumulateChannelDCS!(fAmp::Array{ComplexF64,5}, line::ImpactExcitation.Line, channel::ImpactExcitation.Channel,
                                         angles::Array{Float64,1},
                                         ylmCache::Union{Nothing,Dict{Tuple{Int64,Int64},Vector{ComplexF64}}}=nothing)`
    ... adds the (already computed) amplitude of one channel to the coherent scattering amplitude
        fAmp[theta, M_i, M_f, mu_a, mu_b] for an incoming electron along z and an outgoing one at angle theta in the x-z plane.
        The projection m_b of the outgoing partial wave is not summed over, since the Clebsch-Gordan coefficients fix it:
        M_t = M_i + mu_a, m_b = M_t - M_f, m_lb = m_b - mu_b. The channel amplitude is the ordinary (M-independent) matrix
        element of the scalar e-e interaction, as in the collision strength. Nothing is returned; fAmp is modified in place.
"""
function accumulateChannelDCS!(fAmp::Array{ComplexF64,5}, line::ImpactExcitation.Line, channel::ImpactExcitation.Channel,
                               angles::Array{Float64,1},
                               ylmCache::Union{Nothing,Dict{Tuple{Int64,Int64},Vector{ComplexF64}}}=nothing)
    Ji = line.initialLevel.J.num // line.initialLevel.J.den
    Jf = line.finalLevel.J.num // line.finalLevel.J.den

    iSubshell = Subshell(101, channel.initialKappa)
    fSubshell = Subshell(102, channel.finalKappa)
    la = Basics.subshell_l(iSubshell)
    lb = Basics.subshell_l(fSubshell)
    ja = Basics.subshell_j(iSubshell).num // Basics.subshell_j(iSubshell).den
    jb = Basics.subshell_j(fSubshell).num // Basics.subshell_j(fSubshell).den
    Jt = channel.symmetry.J.num // channel.symmetry.J.den

    for (mubIdx, mub) in enumerate(-1//2:1//2)
        for (muaIdx, mua) in enumerate(-1//2:1//2)
            cgA = AngularMomentum.ClebschGordan(la, 0, 1//2, mua, ja, mua)
            cgA == 0.0 && continue

            for (mfIdx, Mf) in enumerate(-Jf:Jf)
                for (miIdx, Mi) in enumerate(-Ji:Ji)
                    Mt = Mi + mua
                    abs(Mt) > Jt && continue

                    # CG(Jf Mf, jb mb | Jt Mt) forces mb = Mt - Mf.
                    mb = Mt - Mf
                    abs(mb) > jb && continue

                    mlb = mb - mub
                    denominator(mlb) != 1 && continue
                    mlbInt = Int64(numerator(mlb))
                    abs(mlbInt) > lb && continue

                    cgI = AngularMomentum.ClebschGordan(Ji, Mi, ja, mua, Jt, Mt)
                    cgI == 0.0 && continue

                    cgF = AngularMomentum.ClebschGordan(Jf, Mf, jb, mb, Jt, Mt)
                    cgF == 0.0 && continue

                    cgB = AngularMomentum.ClebschGordan(lb, mlbInt, 1//2, mub, jb, mb)
                    cgB == 0.0 && continue

                    ylm = isnothing(ylmCache) ? ComplexF64[ AngularMomentum.sphericalYlm(lb, mlbInt, deg2rad(theta), 0.0) for theta in angles ] :
                                                ImpactExcitation.dcsYlmVector!(ylmCache, lb, mlbInt, angles)

                    factor = cgA * cgB * sqrt((2la + 1) / 4π) * cgI * cgF * channel.amplitude
                    @views fAmp[:, miIdx, mfIdx, muaIdx, mubIdx] .+= factor .* ylm
                end
            end
        end
    end

    return( nothing )
end


"""
`ImpactExcitation.computeDCS(line::ImpactExcitation.Line, angles::Array{Float64,1}=line.angles)`
    ... recomputes the differential cross section of a line from its stored channel amplitudes, with the same accumulator and
        normalization as the on-the-fly computation in computeAmplitudesProperties(); an empty angle list means 0:1:180 degrees.
        A dcs::Array{Float64,1} in a_0^2/sr is returned.
"""
function computeDCS(line::ImpactExcitation.Line, angles::Array{Float64,1}=line.angles)
    anglesx  = isempty(angles) ? collect(0.0:1.0:180.0) : angles
    numMi    = Basics.twice(line.initialLevel.J) + 1;    numMf = Basics.twice(line.finalLevel.J) + 1
    fAmp     = zeros(ComplexF64, length(anglesx), numMi, numMf, 2, 2)
    ylmCache = Dict{Tuple{Int64,Int64},Vector{ComplexF64}}()
    for  channel in line.channels
        ImpactExcitation.accumulateChannelDCS!(fAmp, line, channel, anglesx, ylmCache)
    end

    return( ImpactExcitation.dcsFromAmplitudeTensor(line, fAmp) )
end


"""
`ImpactExcitation.dcsFromAmplitudeTensor(line::ImpactExcitation.Line, fAmp::Array{ComplexF64,5})`
    ... turns the coherent amplitude tensor into the differential cross section, averaged over M_i and mu_a and summed over
        M_f and mu_b. The normalization is tied to that of the integrated cross section of the same line: integrating the result
        over the sphere reproduces line.crossSection to the accuracy of the angle quadrature (3e-5 with a 0.5 degree trapezoid,
        C3+ 2s --> 2p at 50 eV, measured 24-Sep-2026). A factor pi^2/4 had been missing here before that check was made.
        A dcs::Array{Float64,1} in a_0^2/sr is returned.
"""
function dcsFromAmplitudeTensor(line::ImpactExcitation.Line, fAmp::Array{ComplexF64,5})
    matx = zeros(Float64, size(fAmp, 1))
    for  mubIdx = 1:2,  muaIdx = 1:2,  mfIdx = axes(fAmp,3),  miIdx = axes(fAmp,2)
        matx .+= abs2.(fAmp[:, miIdx, mfIdx, muaIdx, mubIdx])
    end

    wc = Defaults.getDefaults("speed of light: c")
    qa = sqrt((line.initialElectronEnergy / wc)^2 + 2 * line.initialElectronEnergy)
    qb = sqrt((line.finalElectronEnergy   / wc)^2 + 2 * line.finalElectronEnergy)
    Ea = sqrt((qa * wc)^2 + wc^4);    Eb = sqrt((qb * wc)^2 + wc^4)
    Ji = line.initialLevel.J.num // line.initialLevel.J.den

    dcsCoeff = 0.5 * sqrt((Ea + wc^2) * (Eb + wc^2) / (Ea * Eb)) / qa / qb
    normDCS  = (2π)^4 * (qb / qa) / (2 * (2Ji + 1))
    dcs      = matx .* normDCS .* dcsCoeff^2 .* qa .* qb

    return( dcs )
end


"""
`ImpactExcitation.computeAmplitudesProperties(line::ImpactExcitation.Line, nm::Nuclear.Model, grid::Radial.Grid,
                                                nrContinuum::Int64, settings::ImpactExcitation.Settings; printout::Bool=true,
                                                nuclearPot::Union{Nothing,Radial.Potential}=nothing,
                                                primitives::Union{Nothing,Bsplines.Primitives}=nothing,
                                                angularCache::Union{Nothing,ImpactExcitation.AngularCoefficientCache}=nothing)`
    ... to compute all amplitudes and properties of the given line; a line::ImpactExcitation.Line is returned for which
    the amplitudes and properties are now evaluated. The two keywords carry quantities that are CONSTANT for a whole
    computation and are otherwise rebuilt for every line and every partial wave: the nuclear potential depends only on
    the nuclear model and the grid, the B-spline primitives only on the grid. Nothing is cached; omitting them
    reproduces the previous behaviour exactly.
"""
function  computeAmplitudesProperties(line::ImpactExcitation.Line, nm::Nuclear.Model, grid::Radial.Grid,
                                        nrContinuum::Int64, settings::ImpactExcitation.Settings; printout::Bool=true,
                                        nuclearPot::Union{Nothing,Radial.Potential}=nothing,
                                        primitives::Union{Nothing,Bsplines.Primitives}=nothing,
                                        angularCache::Union{Nothing,ImpactExcitation.AngularCoefficientCache}=nothing)
    newChannels = ImpactExcitation.Channel[];
    contSettings = Continuum.Settings(settings.continuumExchange === true, nrContinuum)
    cross = 0.;   coll = 0.; convergence = 1.
    # The TTp B-spline matrices depend on the grid and on kappa but never on the potential or the energy, so all the continuum
    # orbitals of ONE line -- both directions and every partial wave, ~35 of them at the maxKappa of a typical case -- can share a
    # single storage dictionary instead of each rebuilding its own. It is created HERE, per line, and deliberately NOT hoisted
    # beside nuclearPot and primitives above the `Threads.@threads` loop that calls this function: a Dict is mutable and is written
    # on every first use of a block, so one shared across threads would be a data race. Per line it is thread-local by construction.
    storage = Dict{String,Array{Float64,2}}()
    conv = 0.; conv0 = 0. ; n = 0

    # Define a common subshell list for both multiplets
    subshellList = Basics.generate(OrderedSubshellList(), line.finalLevel.basis, line.initialLevel.basis)
    ## The standard subshell list is display-only and identical for every line (verified over 50 level pairs),
    ## so it is set once in the drivers rather than here -- this function runs inside a @threads loop, and a
    ## write to a global from there is a race.  See Defaults.setStandardSubshellList.
    ## Both symmetry-reduced levels are properties of the line alone; see the note in the innermost loop.
    redILevel = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, subshellList)
    redFLevel = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel,   subshellList)

    # The radial integrals depend on the continuum orbitals and hence on the energy: one cache per line, owned by one thread.
    radialCache = settings.cacheRadialIntegrals ? InteractionStrength.XLCache() : nothing
    # The free electron sees the target ion's charge beyond the bound orbitals; the tail is continued in that Coulomb field.
    Zion        = nm.Z - sum(line.initialLevel.basis.csfs[1].occupation)
    tailCache   = settings.asymptoticTail ? ImpactExcitation.TailCache(Float64(Zion)) : nothing

    angles = settings.calcDCS ? (isempty(settings.angles) ? collect(0.0:1.0:180.0) : settings.angles) : Float64[]
    numMi = Basics.twice(line.initialLevel.J) + 1
    numMf = Basics.twice(line.finalLevel.J) + 1
    fAmpDCS = settings.calcDCS ? zeros(ComplexF64, length(angles), numMi, numMf, 2, 2) : nothing
    dcsYlmCache = settings.calcDCS ? Dict{Tuple{Int64,Int64},Vector{ComplexF64}}() : nothing

    # First determine a common set of continuum orbitals for the incoming and outgoing electron
    ciOrbitals = Dict{Subshell, Orbital}();     ciPhases = Dict{Subshell, Float64}()
    cfOrbitals = Dict{Subshell, Orbital}();     cfPhases = Dict{Subshell, Float64}()
    #
    symi = LevelSymmetry(line.initialLevel.J, line.initialLevel.parity);    
    symf = LevelSymmetry(line.finalLevel.J, line.finalLevel.parity) 
    #
    kappa = 1;   shellColl = Float64[];   collBefore = 0.;   converged = false
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
                                                                                            primitives=primitives, storage=storage,
                                                                                            orthogonalize=settings.orthogonalizeContinuum,
                                                                                            exchange=ImpactExcitation.continuumExchangeModel(settings))
                            ciOrbitals[iSubshell] = ciOrbital;      ciPhases[iSubshell] = iPhase
                            println("\n Initial channel")
                        end

                        fSubshell  = Subshell(102, channel.finalKappa)
                        if  !haskey(cfOrbitals, fSubshell)
                            cfOrbital, fPhase     = Continuum.generateOrbitalForLevel(line.finalElectronEnergy,   fSubshell,
                                                                                            redFLevel, nm, grid, contSettings;
                                                                                            nuclearPot=nuclearPot,
                                                                                            primitives=primitives, storage=storage,
                                                                                            orthogonalize=settings.orthogonalizeContinuum,
                                                                                            exchange=ImpactExcitation.continuumExchangeModel(settings))
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
                        amplitude  = ImpactExcitation.amplitude(settings.operator, newChannel, newfLevel, newiLevel, grid;
                                                               printout=printout, radialCache=radialCache, angularCache=angularCache,
                                                               tailCache=tailCache)

                        coll       = coll + AngularMomentum.bracket([channel.symmetry.J]) * ( conj(amplitude) * amplitude ).re

                        conv += (conj(amplitude) * amplitude).re

                        chan = ImpactExcitation.Channel(newChannel.initialKappa, newChannel.finalKappa,
                                                        newChannel.symmetry, iPhase, fPhase, amplitude)
                        push!(newChannels, chan)
                        if settings.calcDCS
                            ImpactExcitation.accumulateChannelDCS!(fAmpDCS, line, chan, angles, dcsYlmCache)
                        end
                    end
                end
            end
        end
        
        push!(shellColl, coll - collBefore);   collBefore = coll
        # Checking convergence
        convergence = abs(conv - conv0)/conv
        if convergence < convergenceCriterion
            printstyled("\nConvergence Achieved "; bold=true, underline=true, color=:light_red); 
            if n == 5 
                printstyled("\nConvergence Achieved "; bold=true, underline=true, color=:light_red);   converged = true;   break
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

    if  settings.topUp  &&  !converged    coll = coll + ImpactExcitation.topUpCollisionStrength(shellColl)    end
    coll = coll * N 
    cross = coll * pi / ( (qa^2) * (Basics.twice(line.initialLevel.J) + 1) )

    dcs = settings.calcDCS ? ImpactExcitation.dcsFromAmplitudeTensor(line, fAmpDCS) : Float64[]

    newLine = ImpactExcitation.Line( line.initialLevel, line.finalLevel, line.initialElectronEnergy, line.finalElectronEnergy,
                                     cross, coll, angles, dcs, newChannels, convergence )
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
    angularCache = settings.cacheAngularCoefficients ? ImpactExcitation.AngularCoefficientCache() : nothing
    # Display all selected lines before the computations start
    if  settings.printBefore    ImpactExcitation.displayLines(lines)    end
    ## Determine the maximum (continuum) electron energy and check the grid for consistency with it. This module used to
    ## skip this check -- alone among the continuum-based process modules -- and to normalize at grid.NoPoints-5 instead
    ## of the nrContinuum returned here. Both were wrong: an under-resolved continuum orbital does not fail loudly, it
    ## returns cross sections that drift by tens of percent and then diverge by 30+ orders of magnitude once the grid can
    ## no longer represent the oscillation. Continuum.gridConsistency() refuses such a grid up front and says why.
    ## gridConsistency() tests the FAST end only -- whether grid.hp can still sample the shortest oscillation.
    ## The slow end fails the opposite way and was unguarded until 28-Aug-2026: near threshold the wavelength grows
    ## without bound, and an orbital that cannot complete one oscillation inside rbox has no asymptotic region for
    ## its normalization to be read from. example-Dl.jl branch c) meets exactly this -- 0.008 eV out, lambda ~ 256
    ## a.u. against rbox = 30. It WARNS rather than refusing, because how wrong the near-threshold cross section
    ## then is has not been measured, and refusing on an unmeasured criterion would stop runs that may be fine.
    maxEnergy = 0.;   minEnergy = Inf
    for  line in lines
        maxEnergy = max(maxEnergy, line.initialElectronEnergy);   minEnergy = min(minEnergy, line.finalElectronEnergy)
    end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
    ImpactExcitation.checkLowEnergyGrid(minEnergy, grid)
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
                                                                       nuclearPot=nuclearPot, primitives=primitives,
                                                                       angularCache=angularCache)
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
    angularCache = settings.cacheAngularCoefficients ? ImpactExcitation.AngularCoefficientCache() : nothing
    if  length(lines) == 0    return( ImpactExcitation.Line[] )    end
    ## The grid must be able to carry the fastest continuum orbital that occurs; see the note in computeLines().
    maxEnergy   = 0.;   minEnergy = Inf
    for  line in lines
        maxEnergy = max(maxEnergy, line.initialElectronEnergy);   minEnergy = min(minEnergy, line.finalElectronEnergy)
    end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
    ImpactExcitation.checkLowEnergyGrid(minEnergy, grid)
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
                                                               nuclearPot=nuclearPot, primitives=primitives, angularCache=angularCache)
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
                                                                                  Float64[], Float64[], channels, 1.0) )
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
                                                                                  Float64[], Float64[], channels, 1.0) )
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

    if !isempty(line.dcs)
        println(" ")
        println("  Differential Cross Sections (DCS):")
        println("  ", TableStrings.hLine(50))
        println("    Theta [deg]         DCS [10^-20 m^2/sr]")
        println("  ", TableStrings.hLine(50))
        for (idx, θ) in enumerate(line.angles)
            dcsVal = line.dcs[idx] * 0.5291772109 * 0.5291772109
            @printf("     %8.2f           %.6e\n", θ, dcsVal)
        end
        println("  ", TableStrings.hLine(50))
    end
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
