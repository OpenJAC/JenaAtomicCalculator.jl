
"""
`module  JAC.ResonantImpactIonization`
... a submodel of JAC that provides the RESONANT channels of electron-impact ionization: an incident electron is CAPTURED into a
    doubly-excited state of the next lower charge state, and that resonance then sheds TWO electrons, so that the ion ends up one
    charge state HIGHER than it began.  Two routes are distinguished, and they differ only in how the resonance decays:

    + resonant-electron-capture-with-sequential-double-autoionization
        ... the two electrons leave ONE AFTER THE OTHER, through an intermediate level of the original charge state that lies above
            the ionization threshold and therefore autoionizes in its turn.  Selected by `SequentialAuger()`.
            Known in the older literature as REDA.
    + resonant-electron-capture-with-simultaneous-double-autoionization
        ... the two electrons leave TOGETHER, in a single double-Auger transition.  Selected by `SimultaneousAuger()`.
            Known in the older literature as READI.

    In electron counting, with N the electrons of the initial ion,

        capture     i(N)  +  e-   -->   d(N+1)                        the doubly-excited resonance
        sequential  d(N+1)  -->  n(N)  +  e-   -->  f(N-1)  +  2e-     two Auger steps
        simultaneous  d(N+1)  -->  f(N-1)  +  2e-                      one double-Auger step

    so both are RESONANT SINGLE IONIZATION, adding sharp structure on top of the smooth direct channel and on top of the
    impact-excitation-with-subsequent-autoionization channel of `Cascade.ElectronIonizationScheme`.

    WHY THIS MATTERS BEYOND ADDING A PROCESS.  A doubly-excited resonance either RADIATES -- which is dielectronic recombination --
    or AUTOIONIZES, which is what this module computes.  They are competing decays of the SAME resonance, so a computation that
    generates the resonance obtains the branching between recombination and ionization at the same time, and that branching is what
    a plasma ionization balance needs.

    WHAT THIS MODULE IS, AND IS NOT.  It contains the two channel selectors, an estimate of the double-Auger probability, and the
    strength formulas -- nothing else.  It computes NO amplitudes and NO angular coefficients of its own: every rate it uses is
    handed to it, having been computed by `AutoIonization` and `PhotoEmission` through the cascade.  The strength functions
    therefore take RATES AS ARGUMENTS rather than levels, which is what makes them checkable on numbers written by hand.  It is
    also deliberately independent of `DoubleAutoIonization`: the simultaneous route is estimated here, not computed from
    amplitudes.

    APPROXIMATIONS, stated so that a user knows what a computation actually delivers:

    + isolated resonances
        ... each resonance is treated on its own and the strengths of different resonances are ADDED -- no interference between
            them, and no overlap of their widths.  Sound where the resonances are narrow against their spacing, which holds for
            the low-n captures; it degrades for high n, where the Rydberg resonances crowd together.
    + capture by detailed balance
        ... the capture rate is the Auger rate of the time-reversed transition, which is exact rather than approximate.  But the
            capture therefore inherits whatever partial-wave truncation the Auger calculation used: a `maxKappa` too small lowers
            every strength here in the same proportion, silently.
    + branchings from the rates that were generated
        ... every branching ratio is formed from TOTAL rates summed over the decay channels the cascade actually built.  The
            branchings therefore sum to one BY CONSTRUCTION, and that sum can never reveal a decay route that was left out of the
            configuration lists.  It checks the arithmetic, not the completeness.
    + the two sequential Auger steps are independent
        ... the intermediate is taken to be formed and then to decay, with no memory of its formation and no interference between
            the two time orderings.  This is the isolated-resonance picture applied twice.
    + the simultaneous route is a SHAKE-OFF ESTIMATE, not a computed rate
        ... see `ResonantImpactIonization.shakeProbability`.  Its error has BOTH SIGNS and it must not be quoted as a bound.

    THE SIMULTANEOUS ROUTE DESERVES A WARNING OF ITS OWN, because it is the weakest thing here.  A double-Auger width cannot be had
    without amplitudes, so it is estimated from shake-off in the sudden approximation -- real physics, computed from orbital
    overlaps, but incomplete in two opposite directions at once: a finite bound-orbital set lets some shake-UP leak in and inflates
    it, while KNOCKOUT, which is comparable to or larger than shake-off for double Auger, is missing altogether and deflates it.
    An estimate with errors of both signs is not a bound, and the printout says so rather than offering a number to quote.
"""
module ResonantImpactIonization


using  Printf, ..Basics, ..Defaults, ..Radial, ..RadialIntegrals


"""
`struct  ResonantImpactIonization.SequentialAuger  <:  Basics.AbstractProcess`
    ... selects resonant-electron-capture-with-sequential-double-autoionization, in which the doubly-excited resonance sheds its two
        electrons ONE AFTER THE OTHER, through an intermediate level that autoionizes in its turn.  Known as REDA in the older
        literature.  It subtypes `Basics.AbstractProcess` so that it can stand directly in the `processes` list of a cascade scheme,
        rather than needing a second set of marker types beside it.
"""
struct   SequentialAuger    <:  Basics.AbstractProcess   end


"""
`struct  ResonantImpactIonization.SimultaneousAuger  <:  Basics.AbstractProcess`
    ... selects resonant-electron-capture-with-simultaneous-double-autoionization, in which the doubly-excited resonance sheds both
        electrons TOGETHER, in a single double-Auger transition.  Known as READI in the older literature.  Its rate is ESTIMATED
        from shake-off, see the module docstring and `ResonantImpactIonization.shakeProbability`; it is not computed from
        amplitudes.
"""
struct   SimultaneousAuger  <:  Basics.AbstractProcess   end


# `Base.string(kind::ResonantImpactIonization.SequentialAuger)`  ... provides a String notation for the variable kind.
function Base.string(kind::ResonantImpactIonization.SequentialAuger)
    return( "resonant-electron-capture-with-sequential-double-autoionization" )
end


# `Base.string(kind::ResonantImpactIonization.SimultaneousAuger)`  ... provides a String notation for the variable kind.
function Base.string(kind::ResonantImpactIonization.SimultaneousAuger)
    return( "resonant-electron-capture-with-simultaneous-double-autoionization" )
end


"""
`ResonantImpactIonization.doubleAugerProbability(passiveOrbitals::Array{Radial.Orbital,1}, occupations::Array{Int64,1},`
                                                 `finalOrbitals::Array{Radial.Orbital,1}, grid::Radial.Grid)`
    ... to estimate the probability that a SECOND electron is emitted together with the Auger electron, i.e. the branching fraction
        of the double-Auger route, by combining the one-electron shake probabilities of all passive subshells; a value::Float64 in
        [0,1] is returned.

        With P(a) the shake probability of subshell a and N(a) its occupation, the probability that NO electron is shaken is the
        product of (1 - P(a))^N(a) over the passive subshells, so that

            P(double Auger)  =  1  -  PROD_a  (1 - P(a))^N(a)

        This treats the passive electrons as shaking INDEPENDENTLY of one another, which is the usual sudden-approximation counting
        and is the third approximation stacked on top of the two carried by `shakeProbability` itself.  Read that function's
        docstring before using the number: the estimate is wrong in both directions and is not a bound.

        PASSIVE MEANS PASSIVE, and it is the easiest thing to get wrong here.  The two electrons that PARTICIPATE in the Auger
        transition -- the one that fills the hole and the one that is ejected -- must NOT appear in `passiveOrbitals`, or they are
        counted twice: once as the Auger itself and once as a shake.  For a KLL transition the passive set is everything except the
        two L electrons involved.

        Measured behaviour, for orientation rather than as a claim: with 2s and 2p passive, neon between its neutral and
        singly-ionized orbitals gives 6.2 %, and the one-electron probabilities are 1.0 % for 2p in Ne, 0.77 % for 3p in Ar and
        0.76 % for 4p in Kr, while the 1s is unaffected to five decimals -- a deeply bound orbital hardly notices the loss of an
        outer electron.  Those are the right order for light elements.  Note that a change of NUCLEAR charge is NOT the right
        comparison and gives numbers several times larger; an Auger raises the IONIC charge at fixed Z.
"""
function doubleAugerProbability(passiveOrbitals::Array{Radial.Orbital,1}, occupations::Array{Int64,1},
                                finalOrbitals::Array{Radial.Orbital,1}, grid::Radial.Grid)
    if  length(passiveOrbitals) != length(occupations)
        error("ResonantImpactIonization.doubleAugerProbability(): $(length(passiveOrbitals)) orbitals were given but " *
              "$(length(occupations)) occupations; there must be one occupation per passive subshell.")
    end
    wa = 1.0
    for  (i, orbital)  in  enumerate(passiveOrbitals)
        pa = ResonantImpactIonization.shakeProbability(orbital, finalOrbitals, grid)
        wa = wa * (1.0 - pa)^occupations[i]
    end

    return( 1.0 - wa )
end


"""
`ResonantImpactIonization.resonanceStrength(energy::Float64, captureRate::Float64, twiceJd::Int64, twiceJi::Int64)`
    ... to compute the energy-integrated strength of ONE resonance for its FORMATION alone, i.e. before any decay branching is
        applied; a value::Float64 in atomic units is returned.

            S_0  =  pi^2 / k^2  *  A(capture)  *  (2J_d + 1) / (2J_i + 1)

        with k the wave number of the incident electron at the resonance energy.  The capture rate follows from the Auger rate of
        the time-reversed transition by detailed balance, so no separate capture calculation enters.  The wave number is taken from
        `Defaults.convertUnits("kinetic energy to wave number: atomic units", ...)`, i.e. relativistically and by exactly the same
        route as `Cascade.simulateDrRateCoefficients` uses for the dielectronic-recombination strength -- which is what allows the
        recombination and ionization strengths of the SAME resonance to be compared without a conversion between them.
"""
function resonanceStrength(energy::Float64, captureRate::Float64, twiceJd::Int64, twiceJi::Int64)
    if  energy <= 0.    return( 0. )    end
    k = Defaults.convertUnits("kinetic energy to wave number: atomic units", energy)

    return( pi*pi / (k*k) * captureRate * (twiceJd + 1.0) / (twiceJi + 1.0) )
end


"""
`ResonantImpactIonization.sequentialStrength(energy::Float64, captureRate::Float64, twiceJd::Int64, twiceJi::Int64,`
                                             `augerRateToIntermediate::Float64, totalWidthResonance::Float64,`
                                             `augerRateIntermediate::Float64, totalWidthIntermediate::Float64)`
    ... to compute the energy-integrated strength of ONE resonance for the sequential double-autoionization route, i.e. the
        formation strength multiplied by the two decay branchings it must pass through; a value::Float64 in atomic units is
        returned.

            S  =  S_0  *  [ A(d -> n) / Gamma(d) ]  *  [ A(n, total Auger) / Gamma(n) ]

        The first bracket is the branching of the resonance d into the autoionizing intermediate n, the second is the branching of
        n to autoionize rather than to radiate; Gamma is the TOTAL width, Auger plus radiative, of the level in question.  Zero is
        returned if either total width vanishes, which happens for a level whose decays were not generated.

        The two steps are treated as INDEPENDENT, with no interference between the time orderings; see the module docstring.
"""
function sequentialStrength(energy::Float64, captureRate::Float64, twiceJd::Int64, twiceJi::Int64,
                            augerRateToIntermediate::Float64, totalWidthResonance::Float64,
                            augerRateIntermediate::Float64, totalWidthIntermediate::Float64)
    if  totalWidthResonance <= 0.   ||   totalWidthIntermediate <= 0.    return( 0. )    end
    wa = ResonantImpactIonization.resonanceStrength(energy, captureRate, twiceJd, twiceJi)

    return( wa * (augerRateToIntermediate/totalWidthResonance) * (augerRateIntermediate/totalWidthIntermediate) )
end


"""
`ResonantImpactIonization.shakeProbability(initialOrbital::Radial.Orbital, finalOrbitals::Array{Radial.Orbital,1},`
                                           `grid::Radial.Grid)`
    ... to estimate, in the SUDDEN APPROXIMATION, the probability that a passive electron in `initialOrbital` fails to remain in
        any bound orbital of the ion left behind once the Auger electron has gone; a value::Float64 in [0,1] is returned.

            P  =  1  -  SUM_(n' bound, same kappa)  |< phi(n' kappa) | phi(initial) >|^2

        The sudden change of the potential is spherically symmetric, so it cannot change kappa: only final orbitals of the SAME
        kappa are summed, and `RadialIntegrals.overlap` supplies each overlap.  `finalOrbitals` must be the orbitals of the ion
        AFTER the Auger transition, i.e. one charge state higher than those of `initialOrbital`.

        THE ESTIMATE IS WRONG IN BOTH DIRECTIONS AND IS NOT A BOUND, which is the single most important thing about it:

        + a finite bound set INFLATES it
            ... the sum runs over the bound orbitals that were actually computed.  What the projection loses to bound orbitals
                outside that set is counted here as if the electron had been ejected, so shake-UP leaks in as shake-off.
        + the missing knockout DEFLATES it
            ... shake-off is only one of the two mechanisms of double Auger.  The other, KNOCKOUT -- the departing electron
                striking a second on its way out -- is a collision and cannot be reached from overlaps at all.  It is comparable
                to, and often larger than, shake-off.

        The two errors do not cancel in any controlled way.  A negative value from numerical noise, which an over-complete bound
        set can produce, is clipped to zero.
"""
function shakeProbability(initialOrbital::Radial.Orbital, finalOrbitals::Array{Radial.Orbital,1}, grid::Radial.Grid)
    kappa = initialOrbital.subshell.kappa
    wa    = 0.
    for  orbital  in  finalOrbitals
        if  orbital.subshell.kappa != kappa    continue    end
        wa = wa + RadialIntegrals.overlap(orbital, initialOrbital, grid)^2
    end

    return( max(0., min(1., 1.0 - wa)) )
end


"""
`ResonantImpactIonization.simultaneousStrength(energy::Float64, captureRate::Float64, twiceJd::Int64, twiceJi::Int64,`
                                               `doubleAugerRate::Float64, totalWidthResonance::Float64)`
    ... to compute the energy-integrated strength of ONE resonance for the simultaneous double-autoionization route, i.e. the
        formation strength multiplied by the branching into the double-Auger channel; a value::Float64 in atomic units is returned.

            S  =  S_0  *  A(double Auger) / Gamma(d)

        with Gamma(d) the total width of the resonance.  Zero is returned if that width vanishes.  The double-Auger rate is
        normally `ResonantImpactIonization.doubleAugerProbability(...)` times the resonance's total Auger rate, and is therefore an
        ESTIMATE whose error has both signs; a strength from this function is an order of magnitude, not a value to quote.
"""
function simultaneousStrength(energy::Float64, captureRate::Float64, twiceJd::Int64, twiceJi::Int64,
                              doubleAugerRate::Float64, totalWidthResonance::Float64)
    if  totalWidthResonance <= 0.    return( 0. )    end
    wa = ResonantImpactIonization.resonanceStrength(energy, captureRate, twiceJd, twiceJi)

    return( wa * doubleAugerRate/totalWidthResonance )
end

end # module
