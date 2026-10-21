

"""
`module  JAC.ThreeBodyRecombination`
    ... a submodel of JAC that computes the rate coefficient for three-body recombination, e + e + A^(z+1) ->
        A^z(m) + e, via the Saha/Milne detailed-balance relation to its inverse process, electron-impact
        ionization. THIS IS NOT A FIRST-PRINCIPLES CALCULATION -- no genuine two-free-electron continuum
        scattering amplitude is computed anywhere here, and the module's name should not be read as promising
        one; a real ab-initio three-body rate needs a three-body continuum that essentially no production
        atomic-structure code computes.

        THE PHYSICS. Under LTE (a Maxwellian electron distribution), detailed balance between three-body
        recombination into a specific bound level m and electron-impact ionization of that same level fixes
        the one rate coefficient in terms of the other:

            alpha^TBR(T_e; i->m) = k^EII(T_e; m->i) * (g_m/(2 g_i)) * (2 pi/T_e)^(3/2) * exp(P/T_e)

        (atomic units; g = 2J+1 is a level's statistical weight, the factor 2 is the free electron's own spin
        degeneracy, P = E_i - E_m is the capture energy, and k^EII is the Maxwellian-averaged rate coefficient
        of the inverse ionization). This is the ordinary Saha equation applied to one specific bound level
        rather than to an ensemble average -- see e.g. Mihalas, "Stellar Atmospheres" (2nd ed.), the chapter on
        LTE populations, or any standard plasma-physics text on detailed balance; nothing about the RELATION
        itself is new or particular to this module.

        WHAT THIS MODULE ADDS is the two INPUTS to that relation. P comes from two actually computed
        `Level`s, reached the same way every other JAC process reaches them -- `SelfConsistent`/`Hamiltonian`
        on explicit configurations -- rather than from a scaled or tabulated binding energy. And k^EII is
        evaluated with a SELECTABLE `ImpactIonization.AbstractModel` (BEB, BED, a relativistic variant, ...),
        rather than one fixed formula. See `ThreeBodyRecombination.eiiRateCoefficient` for exactly how the one
        subshell-dependent binding energy that model uses internally (from the orbitals of the bound level's
        OWN basis) relates to P -- the two are independent numbers that should agree closely for a
        well-converged structure, and nothing here forces them to.

        WHAT IS NOT HERE. A real two-free-electron continuum amplitude; shake-off or autoionizing-cascade
        capture; anything beyond the single-subshell capture a `ThreeBodyRecombination.Line` models (see
        `ThreeBodyRecombination.capturedSubshell` for what "single capture" means precisely and what happens
        otherwise).

        STATUS, 21-Oct-2026: a first design, agreed with the maintainer and being built stepwise -- Settings,
        Line and the Saha/Milne core land in this commit; `examples/example-Ta.jl` and a TestFrames test
        follow next.
"""
module ThreeBodyRecombination


using  Printf, FastGaussQuadrature, ..Basics, ..Defaults, ..ImpactIonization, ..ManyElectron, ..Nuclear, ..Radial


"""
`struct  ThreeBodyRecombination.Settings  <:  AbstractProcessSettings`
    ... defines a type for the details and parameters of a three-body-recombination rate-coefficient
        computation.

    + eiiModel             ::ImpactIonization.AbstractModel   ... which empirical electron-impact-ionization
                                                                    cross-section model supplies the inverse
                                                                    process; see `ImpactIonization.AbstractModel`
                                                                    for the available choices.
    + electronTemperature  ::Float64                          ... electron temperature of the surrounding
                                                                    plasma, in Kelvin; three-body recombination
                                                                    is only meaningful under the LTE
                                                                    (Maxwellian-electron) assumption this module
                                                                    makes throughout.
    + lineSelection         ::LineSelection                    ... Specifies the selected (initialLevel,
                                                                    finalLevel) pairs, if any.
    + printBefore            ::Bool                           ... True, if the lines are to be printed before
                                                                    the rate-coefficient computation returns.
"""
struct  Settings  <:  AbstractProcessSettings
    eiiModel                ::ImpactIonization.AbstractModel
    electronTemperature      ::Float64
    lineSelection             ::LineSelection
    printBefore                ::Bool
end


"""
`ThreeBodyRecombination.Settings()`  ... constructor for the default values of a three-body-recombination
    computation.
"""
function  Settings()
    Settings(ImpactIonization.BEBmodel(), 0., LineSelection(), false)
end


"""
`ThreeBodyRecombination.Settings(set::ThreeBodyRecombination.Settings;`

        eiiModel=..,   electronTemperature=..,   lineSelection=..,   printBefore=..)

    ... constructor for modifying the given ThreeBodyRecombination.Settings by 'overwriting' the previously
        selected parameters.
"""
function  Settings(set::ThreeBodyRecombination.Settings;
    eiiModel::Union{Nothing,ImpactIonization.AbstractModel}=nothing,   electronTemperature::Union{Nothing,Float64}=nothing,
    lineSelection::Union{Nothing,LineSelection}=nothing,                printBefore::Union{Nothing,Bool}=nothing)

    if  isnothing(eiiModel)              eiiModelx              = set.eiiModel              else  eiiModelx              = eiiModel              end
    if  isnothing(electronTemperature)   electronTemperaturex   = set.electronTemperature   else  electronTemperaturex   = electronTemperature   end
    if  isnothing(lineSelection)         lineSelectionx         = set.lineSelection         else  lineSelectionx         = lineSelection         end
    if  isnothing(printBefore)           printBeforex           = set.printBefore           else  printBeforex           = printBefore           end

    Settings( eiiModelx, electronTemperaturex, lineSelectionx, printBeforex )
end


function  Base.show(io::IO, settings::ThreeBodyRecombination.Settings)
    println(io, "eiiModel:                 $(settings.eiiModel)  ")
    println(io, "electronTemperature [K]:  $(settings.electronTemperature)  ")
    println(io, "lineSelection:            $(settings.lineSelection)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
end


"""
`struct  ThreeBodyRecombination.Line`
    ... one (initialLevel, finalLevel) capture pair and its three-body-recombination rate coefficient.

    + initialLevel        ::Level      ... the ion BEFORE capture, with one electron fewer than `finalLevel`.
    + finalLevel          ::Level      ... the bound level the two free electrons recombine into.
    + capturedSubshell    ::Subshell   ... the one subshell gained between `initialLevel` and `finalLevel`;
                                            see `ThreeBodyRecombination.capturedSubshell`.
    + bindingEnergy       ::Float64    ... P = initialLevel.energy - finalLevel.energy, in Hartree -- the REAL,
                                            computed capture energy, not a tabulated or scaled value.
    + eiiRateCoefficient  ::Float64    ... the Maxwellian-averaged rate coefficient of the INVERSE process,
                                            k^EII(T_e), in cm^3/s, from `settings.eiiModel`.
    + rateCoefficient     ::Float64    ... the three-body-recombination rate coefficient alpha^TBR(T_e;
                                            initialLevel -> finalLevel), in cm^6/s, via the Saha/Milne
                                            detailed-balance relation -- see this module's own docstring for
                                            the formula.
"""
struct  Line
    initialLevel          ::Level
    finalLevel            ::Level
    capturedSubshell       ::Subshell
    bindingEnergy           ::Float64
    eiiRateCoefficient       ::Float64
    rateCoefficient           ::Float64
end


function  Base.show(io::IO, line::ThreeBodyRecombination.Line)
    println(io, "ThreeBodyRecombination.Line  i = $(line.initialLevel.J)$(line.initialLevel.parity)  ->  " *
                "m = $(line.finalLevel.J)$(line.finalLevel.parity),  capturedSubshell = $(line.capturedSubshell)")
    println(io, "  bindingEnergy P: $(line.bindingEnergy) Hartree")
    println(io, "  k^EII(T_e):      $(line.eiiRateCoefficient) cm^3/s")
    println(io, "  alpha^TBR(T_e):  $(line.rateCoefficient) cm^6/s")
end


"""
`ThreeBodyRecombination.capturedSubshell(initialLevel::Level, finalLevel::Level)`
    ... to find the ONE subshell whose mean occupation in `finalLevel` exceeds its mean occupation in
        `initialLevel` by (approximately) one electron -- the subshell the two free electrons recombine into.
        Both occupations are the CI-weighted means of the two specific levels (`Basics.
        computeMeanSubshellOccupation(sh, [level])`), not an unweighted average over a whole basis, so a
        level dominated by a different CSF than its basis as a whole is still read correctly. Raises an
        informative error if there is not EXACTLY one such subshell, since `finalLevel` is then not related
        to `initialLevel` by a single capture and the pair does not belong in one `ThreeBodyRecombination.
        Line` at all. A sh::Subshell is returned.
"""
function  capturedSubshell(initialLevel::Level, finalLevel::Level)
    gained = Subshell[]

    for  sh  in  finalLevel.basis.subshells
        qFinal = Basics.computeMeanSubshellOccupation(sh, [finalLevel])
        qInit  = sh in initialLevel.basis.subshells ? Basics.computeMeanSubshellOccupation(sh, [initialLevel]) : 0.
        abs( (qFinal - qInit) - 1.0 ) < 0.3   &&   push!(gained, sh)
    end

    length(gained) != 1   &&
        error("ThreeBodyRecombination.capturedSubshell(): found $(length(gained)) subshells gained between " *
              "initialLevel and finalLevel, not exactly one. A ThreeBodyRecombination.Line models a single " *
              "capture and does not apply to this pair as it stands.")

    return( gained[1] )
end


"""
`ThreeBodyRecombination.eiiRateCoefficient(sh::Subshell, finalLevel::Level, nm::Nuclear.Model, grid::Radial.Grid,
                                            settings::ThreeBodyRecombination.Settings; nPoints::Int64=32)`
    ... the Maxwellian-averaged electron-impact-ionization rate coefficient of subshell `sh`, out of the basis
        belonging to `finalLevel`, at `settings.electronTemperature`, in `settings.eiiModel`:

            k^EII(T_e) = sqrt(8/pi) * T_e^(-1/2) * exp(-E_th/T_e) *
                         integral_0^infinity  sigma(E_th + s*T_e) * (E_th + s*T_e) * exp(-s)  ds

        (atomic units throughout). E_th is the subshell binding energy `settings.eiiModel` uses INTERNALLY,
        read from the orbital energies of `finalLevel.basis` -- a Koopmans-like approximation to the true
        capture energy P this module gets instead from two actually computed level energies (see
        `ThreeBodyRecombination.Line.bindingEnergy`); the two should be close for a converged structure, and a
        caller comparing them is comparing the orbital picture against the correlated one, which is
        informative rather than a defect to silence. The substitution E = E_th + s*T_e removes the 0 x
        infinity cancellation that a direct integral from E_th to infinity hits at low T_e, the same trap
        every Maxwellian rate integral in this codebase must guard against. The integral is evaluated by a
        `nPoints`-point `FastGaussQuadrature.gausslegendre` rule mapped onto s in [0, 40] (exp(-40) is
        negligible at double precision); `nPoints` is exposed so a caller can check the quadrature has
        converged rather than having to trust a single fixed order. A value::Float64 (atomic units,
        bohr^3/time_au) is returned.
"""
function  eiiRateCoefficient(sh::Subshell, finalLevel::Level, nm::Nuclear.Model, grid::Radial.Grid,
                               settings::ThreeBodyRecombination.Settings; nPoints::Int64=32)
    basis = finalLevel.basis
    Eth   = - basis.orbitals[sh].energy
    Te_au = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", settings.electronTemperature)

    iiSettings = ImpactIonization.Settings(ImpactIonization.Settings(); model=settings.eiiModel)

    sMax             = 40.0
    (xNodes, wNodes) = FastGaussQuadrature.gausslegendre(nPoints)
    sNodes           = 0.5*sMax .* (xNodes .+ 1.0)
    sWeights         = 0.5*sMax .* wNodes

    integral = 0.0
    for  (s, w)  in  zip(sNodes, sWeights)
        E     = Eth + s*Te_au
        iE_eV = Defaults.convertUnits("energy: from atomic to eV", E)
        cs    = ImpactIonization.CrossSection(sh, iE_eV, 0., Float64[], Float64[])
        newCs = ImpactIonization.computeCrossSections(settings.eiiModel, cs, basis, nm, grid, iiSettings)
        integral += w * newCs.partialCS * E * exp(-s)
    end

    rate = sqrt(8/pi) / sqrt(Te_au) * exp(-Eth/Te_au) * integral

    return( rate )
end


"""
`ThreeBodyRecombination.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model,
                                      grid::Radial.Grid, settings::ThreeBodyRecombination.Settings;
                                      output::Bool=true)`
    ... to determine, for every (initialLevel, finalLevel) pair selected by `settings.lineSelection` from
        `initialMultiplet` (the ion BEFORE capture, one electron fewer) and `finalMultiplet` (the bound level
        the two free electrons recombine into), the three-body-recombination rate coefficient at
        `settings.electronTemperature` via the Saha/Milne detailed-balance relation against the
        electron-impact-ionization rate coefficient of the inverse process, in `settings.eiiModel` -- see this
        module's own docstring for the formula. A PAIR NOT RELATED BY A SINGLE CAPTURE (zero or more than one
        subshell gained, see `ThreeBodyRecombination.capturedSubshell`) or a pair with `finalLevel.energy >=
        initialLevel.energy` (not actually bound relative to this initial level, so P would not be positive)
        is skipped rather than raising, since a general (initialMultiplet, finalMultiplet) pair of multiplets
        is expected to contain such non-matching combinations alongside the real capture channels. An
        Array{ThreeBodyRecombination.Line,1} is returned if `output=true`, otherwise `nothing`.
"""
function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                        settings::ThreeBodyRecombination.Settings; output::Bool=true)
    lines     = ThreeBodyRecombination.Line[]
    toCm3     = Defaults.convertUnits("length: from atomic to cm", 1.0)^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)
    toCm6     = Defaults.convertUnits("length: from atomic to cm", 1.0)^6 / Defaults.convertUnits("time: from atomic to sec", 1.0)
    Te_au     = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", settings.electronTemperature)

    for  i  in  initialMultiplet.levels
        for  m  in  finalMultiplet.levels
            Basics.selectLevelPair(i, m, settings.lineSelection)   ||   continue
            P = i.energy - m.energy
            P <= 0.   &&   continue

            sh = try   ThreeBodyRecombination.capturedSubshell(i, m)   catch;   nothing   end
            isnothing(sh)   &&   continue

            kEII_au  = ThreeBodyRecombination.eiiRateCoefficient(sh, m, nm, grid, settings)
            gRatio   = (Basics.twice(m.J) + 1) / (2*(Basics.twice(i.J) + 1))
            alpha_au = kEII_au * gRatio * (2*pi/Te_au)^1.5 * exp(P/Te_au)

            push!( lines, ThreeBodyRecombination.Line(i, m, sh, P, kEII_au*toCm3, alpha_au*toCm6) )
        end
    end

    if  settings.printBefore
        for  line in lines   println(line)   end
    end

    output   &&   return( lines )

    return( nothing )
end


end # module
