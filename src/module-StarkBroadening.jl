
"""
`module  JAC.StarkBroadening`
    ... a submodel of JAC that computes the ELECTRON-IMPACT contribution to the Stark broadening of a
        spectral line -- the width and shift a line picks up because free electrons in a plasma repeatedly
        perturb the radiating ion's upper and lower level, interrupting the phase of the emitted/absorbed
        radiation (Lindholm-Foley / Anderson-Baranger impact theory, in its semiclassical-perturbation form,
        as used throughout the non-hydrogenic Stark-broadening literature -- Griem, "Spectral Line
        Broadening by Plasmas" (1974); Seaton, Mon. Not. R. Astron. Soc. 118, 504 (1958); Dimitrijevic &
        Konjevic's modified-semiempirical papers of the 1980s build on exactly this chain).

        WHY THIS NEEDS NO NEW SCATTERING ENGINE. The semiclassical route does not touch an R-matrix or a
        close-coupling calculation at all -- it needs only the SAME structure quantities JAC already
        computes to publication quality: level energies and reduced E1 (optionally E2) matrix elements to
        every OTHER level the upper and lower level of the line are dipole-connected to. Each such
        "perturbing transition" contributes to the width in proportion to its own oscillator strength and
        inversely to its own energy gap, summed over a Maxwellian electron-velocity distribution. This is
        the entire reason the competition-list item names "computed from JAC's own structure" rather than a
        new capability: `StarkBroadening` calls `PhotoExcitation.computeLines` directly for every candidate
        perturbing pair, reusing its already cross-checked-against-Einstein-A oscillator strength, rather
        than re-deriving a transition-moment formula a second time.

        TWO PERTURBER POOLS, READ OFF THE SAME TWO MULTIPLETS `PhotoEmission`/`PhotoExcitation` ALREADY
        TAKE. For a line between an upper level u (`initialMultiplet`) and a lower level l
        (`finalMultiplet`), the E1 selection rule means every dipole-connected perturber of u has the
        SAME parity as l, and every dipole-connected perturber of l has the SAME parity as u -- i.e. the
        perturbers of EACH endpoint of the line already sit inside the OTHER multiplet (plus, in general,
        further levels of that same parity inside either multiplet if more than one configuration was
        given). So no third, separate "perturber" configuration list is introduced: the pool searched for
        perturbers of either level is simply the union of `initialMultiplet.levels` and
        `finalMultiplet.levels`, excluding the level itself. A caller wanting a richer perturber set gets
        it the same way `AdjustedLevels`' own examples reach a richer CSF space -- by handing in a richer
        configuration list upstream, not by a separate argument here.

        STATUS, 21-Oct-2026:
        1. The SCAFFOLD (Settings, the perturber-gathering, the oscillator-strength-based
           `PerturberContribution`s) is COMPLETE and exercised end to end through `perform()`.
        2. ONE WIDTH METHOD SHIPS: `Basics.SimplifiedModifiedSemiempiricalMethod` (SMSE), Dimitrijevic &
           Konjevic, Astron. Astrophys. 172, 345 (1987), formula as restated in Dimitrijevic & Popovic, J.
           Appl. Spectrosc. 68, 893 (2001), Eq. (7). This formula was deliberately NOT written down from
           memory -- the primary sources were fetched and cross-checked FIRST (the Z-dependent Gaunt-factor
           threshold "0.9-1.1/Z" in that equation was independently reproduced by combining a SEPARATE,
           single-source 1980 table value with a SEPARATE, independently-sourced Z-correction formula,
           which is as strong a confirmation as was practical to obtain without the original 1968 Phys.
           Rev. reference). SMSE needs no perturber search at all -- only the line's own two levels, via
           `StarkBroadening.effectiveQuantumNumber` and `StarkBroadening.activeSubshell`.
        3. SMSE's SHIFT IS NOT COMPUTED -- the companion Eq. (8)/(9) in the same source carries a sign
           factor this module's reading could not pin down with confidence; see
           `StarkBroadening.electronWidth`'s own docstring. `Line.electronShift` stays `nothing`
           unconditionally until that is resolved.
        4. NOT YET VALIDATED AGAINST A PUBLISHED NUMBER. The formula is now implemented and runs, but no
           `examples/*.jl` branch carries a `Last successful` date yet -- that needs checking a real case
           against a STARK-B/literature value first, the same discipline every other dated branch in this
           project follows.
        5. The FULL Modified Semiempirical method (perturber-by-perturber, via the Gaunt-factor table
           rather than SMSE's closed-form limit) and the ion/quasi-static (microfield) half of Stark
           broadening are BOTH out of scope here and were never started.
        A session picking this up should read points 1-5 above, and should treat SMSE's own numbers as
        "runs correctly, not yet checked against a published width" rather than as settled.
"""
module StarkBroadening


using Printf, ..Basics, ..Defaults, ..ManyElectron, ..PhotoExcitation, ..Radial


"""
`struct  StarkBroadening.Settings  <:  AbstractProcessSettings`
    ... defines a type for the details and parameters of an electron-impact Stark-broadening computation.

    + multipoles              ::Array{EmMultipole,1}   ... which multipoles couple a level to its perturbers;
                                                             E1 alone is the standard choice, since E2
                                                             perturbers are a small, same-parity correction.
    + gauges                  ::Array{UseGauge,1}      ... Gauges to be included for every perturber's own
                                                             oscillator strength.
    + electronDensity         ::Float64                ... electron number density, in cm^-3.
    + electronTemperature     ::Float64                ... electron temperature, in Kelvin.
    + perturberCutoffEnergy   ::Float64                ... perturbers farther than this from the reference
                                                             level, in Hartree, are excluded; 0. (the
                                                             default) means NO cutoff -- every dipole-
                                                             connected level found in the two given
                                                             multiplets is kept.
    + lineSelection           ::LineSelection           ... Specifies the selected lines, if any.
    + method                  ::Union{Nothing,Basics.AbstractStarkBroadeningMethod}
                                                         ... which closed-form relation turns the gathered
                                                             perturbers (or the line's own levels) into an
                                                             actual width/shift; `nothing` (the default)
                                                             means the scaffold runs -- perturbers and
                                                             oscillator strengths are gathered, but
                                                             `Line.electronWidth` stays `nothing` too. See
                                                             `Basics.AbstractStarkBroadeningMethod`.
    + printBefore             ::Bool                   ... True, if the lines (with their perturbers) are
                                                             to be printed before the width computation.
"""
struct  Settings  <:  AbstractProcessSettings
    multipoles                ::Array{EmMultipole,1}
    gauges                    ::Array{UseGauge,1}
    electronDensity            ::Float64
    electronTemperature        ::Float64
    perturberCutoffEnergy      ::Float64
    lineSelection               ::LineSelection
    method                      ::Union{Nothing,Basics.AbstractStarkBroadeningMethod}
    printBefore                  ::Bool
end


"""
`StarkBroadening.Settings()`  ... constructor for the default values of an electron-impact Stark-broadening
    computation.
"""
function  Settings()
    Settings(EmMultipole[E1], UseGauge[Basics.UseCoulomb], 0., 0., 0., LineSelection(), nothing, false)
end


"""
`StarkBroadening.Settings(set::StarkBroadening.Settings;`

        multipoles=..,            gauges=..,                   electronDensity=..,
        electronTemperature=..,   perturberCutoffEnergy=..,    lineSelection=..,
        method=..,                printBefore=..)

    ... constructor for modifying the given StarkBroadening.Settings by 'overwriting' the previously
        selected parameters.
"""
function  Settings(set::StarkBroadening.Settings;
    multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,       gauges::Union{Nothing,Array{UseGauge,1}}=nothing,
    electronDensity::Union{Nothing,Float64}=nothing,                electronTemperature::Union{Nothing,Float64}=nothing,
    perturberCutoffEnergy::Union{Nothing,Float64}=nothing,          lineSelection::Union{Nothing,LineSelection}=nothing,
    method::Union{Nothing,Basics.AbstractStarkBroadeningMethod,Missing}=missing,
    printBefore::Union{Nothing,Bool}=nothing)

    if  isnothing(multipoles)              multipolesx              = set.multipoles              else  multipolesx              = multipoles              end
    if  isnothing(gauges)                  gaugesx                  = set.gauges                  else  gaugesx                  = gauges                  end
    if  isnothing(electronDensity)         electronDensityx         = set.electronDensity         else  electronDensityx         = electronDensity         end
    if  isnothing(electronTemperature)     electronTemperaturex     = set.electronTemperature     else  electronTemperaturex     = electronTemperature     end
    if  isnothing(perturberCutoffEnergy)   perturberCutoffEnergyx   = set.perturberCutoffEnergy   else  perturberCutoffEnergyx   = perturberCutoffEnergy   end
    if  isnothing(lineSelection)           lineSelectionx           = set.lineSelection           else  lineSelectionx           = lineSelection           end
    ## `missing`, not `nothing`, is the sentinel here: `nothing` is itself a meaningful VALUE for this field
    ## (scaffold-only), so it must stay distinguishable from "the caller didn't mention `method` at all".
    if  ismissing(method)                  methodx                  = set.method                  else  methodx                  = method                  end
    if  isnothing(printBefore)             printBeforex             = set.printBefore             else  printBeforex             = printBefore             end

    Settings( multipolesx, gaugesx, electronDensityx, electronTemperaturex, perturberCutoffEnergyx, lineSelectionx, methodx, printBeforex )
end


function  Base.show(io::IO, settings::StarkBroadening.Settings)
    println(io, "multipoles:               $(settings.multipoles)  ")
    println(io, "gauges:                   $(settings.gauges)  ")
    println(io, "electronDensity [cm^-3]:  $(settings.electronDensity)  ")
    println(io, "electronTemperature [K]:  $(settings.electronTemperature)  ")
    println(io, "perturberCutoffEnergy:    $(settings.perturberCutoffEnergy)  ")
    println(io, "lineSelection:            $(settings.lineSelection)  ")
    println(io, "method:                   $(settings.method)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
end


"""
`struct  StarkBroadening.PerturberContribution`
    ... one OTHER level that a reference level (the upper or lower level of a line) is dipole-connected to,
        together with the oscillator strength a width/shift formula turns an energy gap into a contribution
        from -- see `Basics.AbstractStarkBroadeningMethod` for which formula, if any, actually uses it.

    + level         ::Level        ... the perturbing level.
    + omega         ::Float64      ... |E_perturber - E_reference|, in Hartree.
    + oscStrength   ::EmProperty   ... the absorption oscillator strength between the reference level and
                                        this perturber, from `PhotoExcitation.computeLines` directly (both
                                        gauges) -- not re-derived here.
"""
struct  PerturberContribution
    level          ::Level
    omega          ::Float64
    oscStrength    ::EmProperty
end


"""
`struct  StarkBroadening.LineProfile`
    ... a sampled line-shape: a detuning grid and the intensity at each point, normalized to a peak of 1.
        When `Line.electronWidth` comes from `Basics.SimplifiedModifiedSemiempiricalMethod`, this is a plain
        LORENTZIAN of that width (no shift, no ion/quasi-static contribution -- see that method's own
        docstring for exactly what it does and does not include); a method combining several broadening
        mechanisms would return something else here, which is the whole reason this is its own sampled
        array rather than a single number.

    + omega        ::Array{Float64,1}   ... detuning from line center, in Hartree.
    + intensity    ::Array{Float64,1}   ... intensity at each point, peak-normalized to 1.
"""
struct  LineProfile
    omega        ::Array{Float64,1}
    intensity    ::Array{Float64,1}
end


"""
`struct  StarkBroadening.Line`
    ... one spectral line and its electron-impact broadening information.

    + initialLevel     ::Level                             ... the line's upper level u.
    + finalLevel       ::Level                              ... the line's lower level l.
    + omega            ::Float64                            ... the line's own transition energy, Hartree.
    + oscStrength      ::EmProperty                          ... the line's OWN absorption oscillator
                                                                   strength (u<-l), from `PhotoExcitation`
                                                                   directly -- also what makes this a REAL
                                                                   line in the first place: `computeLines`
                                                                   never returns a (u,l) pair this is zero
                                                                   for, since a dipole-forbidden pair (e.g.
                                                                   J=0 to J=0, or J=0 to J=2) is not a
                                                                   spectral line to broaden at all.
    + perturbersUpper  ::Array{PerturberContribution,1}      ... every level u is dipole-connected to
                                                                   (other than l), found in the combined
                                                                   initial+final multiplet pool.
    + perturbersLower  ::Array{PerturberContribution,1}      ... the same for l.
    + electronWidth    ::Union{Nothing,Float64}              ... the electron-impact FWHM width, in Hartree
                                                                    (angular frequency), from
                                                                    `settings.method` -- `nothing` if
                                                                    `settings.method === nothing` (the
                                                                    scaffold-only case). NOT an `EmProperty`:
                                                                    this is a semi-empirical/semiclassical
                                                                    number, not a JAC transition amplitude,
                                                                    so it carries no Coulomb/Babushkin gauge
                                                                    at all.
    + electronShift    ::Union{Nothing,Float64}              ... NOT YET COMPUTED by any method shipped so
                                                                    far -- see `Basics.
                                                                    SimplifiedModifiedSemiempiricalMethod`'s
                                                                    own docstring for why the shift was left
                                                                    for a later pass.
    + profile          ::Union{Nothing,LineProfile}          ... a Lorentzian built from `electronWidth`
                                                                    (zero shift) when a method computed one,
                                                                    `nothing` otherwise.
"""
struct  Line
    initialLevel      ::Level
    finalLevel        ::Level
    omega             ::Float64
    oscStrength        ::EmProperty
    perturbersUpper    ::Array{PerturberContribution,1}
    perturbersLower    ::Array{PerturberContribution,1}
    electronWidth       ::Union{Nothing,Float64}
    electronShift       ::Union{Nothing,Float64}
    profile             ::Union{Nothing,LineProfile}
end


function  Base.show(io::IO, line::StarkBroadening.Line)
    println(io, "StarkBroadening.Line  u = $(line.initialLevel.J)$(line.initialLevel.parity)  ->  " *
                "l = $(line.finalLevel.J)$(line.finalLevel.parity),  omega = $(line.omega) Hartree")
    println(io, "  oscStrength (u<-l): $(line.oscStrength)")
    println(io, "  perturbers of u: $(length(line.perturbersUpper)),  perturbers of l: $(length(line.perturbersLower))")
    if  isnothing(line.electronWidth)
        println(io, "  electronWidth: nothing (settings.method === nothing -- scaffold only)")
    else
        println(io, "  electronWidth: $(line.electronWidth) Hartree   (shift: NOT YET COMPUTED by any method)")
    end
end


"""
`StarkBroadening.gatherPerturbers(reference::Level, pool::Array{Level,1}, grid::Radial.Grid, settings::StarkBroadening.Settings)`
    ... to find every level in `pool` (other than `reference` itself) that carries a nonzero oscillator
        strength to `reference` under the requested `settings.multipoles`/`settings.gauges`, via
        `PhotoExcitation.computeLines` -- the SAME, already cross-checked-against-Einstein-A oscillator
        strength every other JAC module uses, not a new formula. `settings.perturberCutoffEnergy` (if
        nonzero) excludes a perturber farther away in energy than that bound. An
        Array{PerturberContribution,1} is returned, sorted by ascending `omega`.
"""
function  gatherPerturbers(reference::Level, pool::Array{Level,1}, grid::Radial.Grid, settings::StarkBroadening.Settings)
    contributions = PerturberContribution[]
    pxSettings    = PhotoExcitation.Settings(PhotoExcitation.Settings(); multipoles=settings.multipoles, gauges=settings.gauges)

    for  lev in pool
        lev === reference   &&   continue
        omega = abs(lev.energy - reference.energy)
        omega == 0.   &&   continue
        settings.perturberCutoffEnergy > 0.   &&   omega > settings.perturberCutoffEnergy   &&   continue

        upperLev, lowerLev = lev.energy >= reference.energy ? (lev, reference) : (reference, lev)
        pxLines = PhotoExcitation.computeLines( Multiplet("upper", [upperLev]), Multiplet("lower", [lowerLev]), grid, pxSettings;
                                                  output=true )
        isempty(pxLines)   &&   continue
        osc = pxLines[1].oscStrength
        ( osc.Coulomb == 0.   &&   osc.Babushkin == 0. )   &&   continue

        push!( contributions, PerturberContribution(lev, omega, osc) )
    end

    return( sort(contributions, by = pc -> pc.omega) )
end


"""
`StarkBroadening.activeSubshell(level::Level)`
    ... to find the ONE subshell whose occupation is neither zero nor full in `level`'s own CSFs -- the
        single "active" (valence) subshell both `Basics.SimplifiedModifiedSemiempiricalMethod` and the full
        MSE method are built around. Raises an informative error if there is not EXACTLY one such subshell,
        since that is exactly the physical scope these semi-empirical formulas assume (a single dominant
        valence electron) -- not a limitation of this function to be worked around, but the method's own
        boundary made visible rather than silently misapplied. A sh::Subshell is returned.
"""
function  activeSubshell(level::Level)
    idx = findfirst( c -> c.J == level.J  &&  c.parity == level.parity, level.basis.csfs )
    isnothing(idx)   &&   error("StarkBroadening.activeSubshell(): no CSF of this level's own symmetry found in its basis.")
    csf  = level.basis.csfs[idx]
    open = [ sh   for  (i, sh)  in  enumerate(level.basis.subshells)
             if   0 < csf.occupation[i] < Basics.subshell_2j(sh) + 1 ]

    length(open) != 1   &&
        error("StarkBroadening.activeSubshell(): found $(length(open)) partially-occupied subshells for this " *
              "level, not exactly one. The modified-semiempirical methods are built for a single dominant " *
              "valence electron and do not apply here as they stand.")

    return( open[1] )
end


"""
`StarkBroadening.effectiveQuantumNumber(level::Level, method::Basics.SimplifiedModifiedSemiempiricalMethod)`
    ... to compute n* = Z*sqrt(E_H/(ionizationEnergy - level.energy)), E_H = 0.5 Hartree, from
        `method.ionizationEnergy` and `method.residualCharge` -- the standard Rydberg-series effective
        principal quantum number, exactly as used in Dimitrijevic & Konjevic (1987)/Dimitrijevic & Popovic
        (2001) Eq. (7). Raises if `level.energy` is not below `method.ionizationEnergy`, since n* is not
        real otherwise and that signals a zero-point mismatch between the two, not a physical result. A
        value::Float64 is returned.
"""
function  effectiveQuantumNumber(level::Level, method::Basics.SimplifiedModifiedSemiempiricalMethod)
    bindingEnergy = method.ionizationEnergy - level.energy
    bindingEnergy <= 0.   &&
        error("StarkBroadening.effectiveQuantumNumber(): level.energy is not below method.ionizationEnergy -- " *
              "check that both are on the SAME absolute energy scale.")

    return( method.residualCharge * sqrt(0.5/bindingEnergy) )
end


"""
`StarkBroadening.electronWidth(line::StarkBroadening.Line, method::Basics.SimplifiedModifiedSemiempiricalMethod,
                                 settings::StarkBroadening.Settings)`
    ... the Simplified Modified Semiempirical (SMSE) width, Dimitrijevic & Konjevic, Astron. Astrophys. 172,
        345 (1987), Eq. (7) as restated in Dimitrijevic & Popovic, J. Appl. Spectrosc. 68, 893 (2001):
        needing only the line's own upper and lower level (via their effective quantum numbers and active-
        subshell orbital angular momenta), NOT the gathered perturbers -- this is the method's own defining
        simplification, valid when the nearest real dipole-allowed perturber is far enough that its own
        Gaunt factor has reached the x<=2 threshold constant. `settings.electronDensity` (cm^-3) and
        `.electronTemperature` (K) enter exactly as in the cited formula. THE SHIFT IS NOT COMPUTED: the
        paper's own Eq. (8)/(9) for it carries a sign factor (their epsilon_j) this module's reading of the
        source could not pin down with confidence, and guessing a sign is worse than leaving it undone. A
        value::Float64 (Hartree, angular frequency) is returned.
"""
function  electronWidth(line::StarkBroadening.Line, method::Basics.SimplifiedModifiedSemiempiricalMethod,
                          settings::StarkBroadening.Settings)
    Z     = method.residualCharge
    nStarU = effectiveQuantumNumber(line.initialLevel, method)
    nStarL = effectiveQuantumNumber(line.finalLevel,   method)
    lU     = Basics.subshell_l( activeSubshell(line.initialLevel) )
    lL     = Basics.subshell_l( activeSubshell(line.finalLevel)  )

    sumTerm = (3*nStarU^2/(2*Z^2)) * (nStarU^2 - lU*(lU+1) - 1) +
              (3*nStarL^2/(2*Z^2)) * (nStarL^2 - lL*(lL+1) - 1)

    cAU        = Defaults.getDefaults("speed of light: c")
    lambdaAU   = 2*pi*cAU / line.omega                                 ## standard omega = 2 pi c / lambda, atomic units
    bohrInCm   = Defaults.BOHR_RADIUS_SI * 100.0
    lambda_cm  = lambdaAU * bohrInCm

    widthAngstrom = 2.2151e-8 * lambda_cm^2 * settings.electronDensity / sqrt(settings.electronTemperature) *
                    (0.9 - 1.1/Z) * sumTerm

    ## Convert the WIDTH, a wavelength interval (Angstrom), into an angular-frequency width (Hartree):
    ## |domega/dlambda| = 2 pi c / lambda^2, evaluated at the line's own wavelength, same convention
    ## throughout (atomic units in, atomic units out) -- an ordinary unit conversion, not broadening theory.
    deltaLambda_au = (widthAngstrom * 1.0e-8) / bohrInCm
    widthAU        = (2*pi*cAU / lambdaAU^2) * deltaLambda_au

    return( widthAU )
end


"""
`StarkBroadening.lorentzianProfile(width::Float64; nPoints::Int64=201, span::Float64=20.0)`
    ... a plain, peak-normalized Lorentzian of FWHM `width` (Hartree), sampled over `+-span/2 * width` around
        line center -- the profile `Basics.SimplifiedModifiedSemiempiricalMethod` produces, since it supplies
        no shift and no ion/quasi-static contribution (see that method's own docstring). A
        lp::StarkBroadening.LineProfile is returned.
"""
function  lorentzianProfile(width::Float64; nPoints::Int64=201, span::Float64=20.0)
    halfWidth = width/2
    omegaGrid = collect( range(-span*halfWidth, span*halfWidth, length=nPoints) )
    intensity = [ halfWidth^2 / (halfWidth^2 + om^2)   for  om  in  omegaGrid ]

    return( LineProfile(omegaGrid, intensity) )
end


"""
`StarkBroadening.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid,
                                settings::StarkBroadening.Settings; output::Bool=true)`
    ... to determine, for every (u,l) level pair selected by `settings.lineSelection` from `initialMultiplet`
        (upper levels u) and `finalMultiplet` (lower levels l), the perturbers of u and of l within the
        combined pool of both multiplets' levels -- see this module's own docstring for why no separate
        perturber configuration is needed. A (u,l) PAIR THAT IS ITSELF DIPOLE-FORBIDDEN (zero oscillator
        strength under `settings.multipoles` -- e.g. J=0 to J=0, or J=0 to J=2) is NOT a spectral line and
        is skipped outright, never returned as a `Line` with an empty or zero rate. If `settings.method` is
        given, each `Line.electronWidth` and `.profile` are filled in by that method (see
        `Basics.AbstractStarkBroadeningMethod`); otherwise they stay `nothing` (the scaffold-only case). An
        Array{StarkBroadening.Line,1} is returned if `output=true`, otherwise `nothing`.
"""
function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid,
                        settings::StarkBroadening.Settings; output::Bool=true)
    pool       = vcat(initialMultiplet.levels, finalMultiplet.levels)
    lines      = StarkBroadening.Line[]
    pxSettings = PhotoExcitation.Settings(PhotoExcitation.Settings(); multipoles=settings.multipoles, gauges=settings.gauges)

    for  u in initialMultiplet.levels
        for  l in finalMultiplet.levels
            Basics.selectLevelPair(u, l, settings.lineSelection)   ||   continue
            omega = u.energy - l.energy
            omega <= 0.   &&   continue

            pxLines = PhotoExcitation.computeLines( Multiplet("u", [u]), Multiplet("l", [l]), grid, pxSettings; output=true )
            isempty(pxLines)   &&   continue
            osc = pxLines[1].oscStrength
            ( osc.Coulomb == 0.   &&   osc.Babushkin == 0. )   &&   continue   ## not a real (dipole-forbidden) line

            perturbersUpper = gatherPerturbers(u, pool, grid, settings)
            perturbersLower = gatherPerturbers(l, pool, grid, settings)
            line = StarkBroadening.Line(u, l, omega, osc, perturbersUpper, perturbersLower, nothing, nothing, nothing)

            if  !isnothing(settings.method)
                width = electronWidth(line, settings.method, settings)
                line  = StarkBroadening.Line(u, l, omega, osc, perturbersUpper, perturbersLower,
                                               width, nothing, lorentzianProfile(width))
            end
            push!( lines, line )
        end
    end

    if  settings.printBefore
        for  line in lines   println(line)   end
    end

    output   &&   return( lines )

    return( nothing )
end


end # module
