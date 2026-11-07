
"""
`module  JAC.PhotoIonizationPseudostates`
... a submodel of JAC that treats resonant photoionization by DIAGONALIZING ONCE a CI matrix that mixes the
    bound autoionizing ("doubly-excited") configurations together with the ejected electron's BOX-DISCRETIZED
    continuum -- a ladder of pseudostate orbitals of one fixed partial wave kappa, all mutually orthogonal
    members of ONE generalized one-electron eigenproblem (`Bsplines.diagonalizeLocalMatrix`) -- rather than
    computing the direct and the resonant amplitude SEPARATELY and adding them through Fano's closed-form
    q-formula, as `PhotoIonizationInterference` does.

    THIS IS THE DIRECT ALTERNATIVE TO `PhotoIonizationInterference`, NOT A REPLACEMENT OF IT, built for the
    SAME physical observable after that module's own q-extraction got stuck on a phase-convention mismatch
    among three separately-computed amplitudes (challenges-list S29). Diagonalizing once removes the need for
    the q-formula, and with it the phase convention it assumes, entirely: the interference between the direct
    and the resonant route lives IN THE EIGENVECTORS of the combined CI, not in a sum of two amplitudes that
    must be put on a common phase footing by hand.

    WHERE THE BOUND, AUTOIONIZING CHARACTER COMES FROM -- deliberately NOT a `Settings` field here. A caller
    hands `computeSpectra` an `initialLevel::Level` whose own `.basis` may ALREADY be a multi-configuration CI
    spanning both the ordinary bound configuration and any doubly-excited ("autoionizing") ones -- exactly how
    `PhotoIonizationInterference`'s own `intermediateMultiplet` is built, just folded into ONE multiplet
    instead of two. Everything already sitting in `initialLevel.basis.csfs` is carried, unchanged, into the
    combined basis; only the PSEUDOSTATE part is new. This is simpler than an earlier sketch that gave this
    module its own `autoionizingConfigs` field: that field would have duplicated what `initialLevel.basis`
    already knows, and merging it at the `Configuration`/`Shell` level would have reintroduced exactly the
    kappa-splitting problem the pseudostate ladder itself had to avoid (see below).

    WHY THE PSEUDOSTATE LADDER CANNOT USE `AtomicState.RasStep`/`Basics.Configuration`. A pseudostate ladder is,
    by construction, ONE fixed kappa -- the ejected electron's own partial wave. `Basics.shellSplitOccupation`
    always splits a non-relativistic `l > 0` shell into BOTH kappa branches (`j = l-1/2` and `j = l+1/2`), so
    routing the ladder through the ordinary `Configuration`/`Shell` machinery would silently double every
    pseudostate into two unwanted kappa families. The ladder is therefore built directly at the `ConfigurationR`
    level -- one relativistic configuration per retained pseudostate orbital, each the final ion's own core
    occupation plus one electron in `Subshell(1000+k, kappa)` -- merged with `initialLevel.basis`'s own CSFs
    through the SAME `Basics.generateSubshellList`/`Basics.generateCsfRs`/`Basis(...)` sequence
    `Basics.generateBasis` itself uses, so no new CSF machinery is written, only composed differently.

    ONE CONSISTENT ONE-ELECTRON BASIS, SO NO `BiOrthogonal` IS NEEDED. The pseudostate ladder is generated in
    the FROZEN local potential of `finalIonLevel` (nuclear + `Basics.DFSField(1.0)` mean field, the same
    combination `Continuum.generateOrbitalForLevel` already uses), via the SAME generalized eigenproblem
    (`Bsplines.diagonalizeLocalMatrix`) that an ordinary hydrogenic SCF start already solves -- the pseudostate
    orbitals are mutually orthogonal to each other by construction. They are explicitly Gram-Schmidt
    orthogonalized against any OCCUPIED bound orbital of the SAME kappa already in `finalIonLevel.basis`
    (`removeOccupiedComponent`, mirroring `RandomPhase.projectOntoUnoccupied`'s own, already-validated
    algorithm) -- a core electron cannot be excited into a filled orbital, and without this step a pseudostate
    whose kappa coincides with an occupied core shell would not be orthogonal to it. `checkFrozenOrbitalConsistency`
    then RAISES, rather than silently mixing inconsistent orbitals, if `initialLevel.basis` and `finalIonLevel.basis`
    disagree (overlap below `1.0 - 1.0e-6`) on any subshell they share -- the caller's responsibility is to have
    generated both multiplets against ONE consistent core, exactly as every frozen-orbital CI in JAC already
    assumes elsewhere.

    **THE ONE HONEST LIMITATION THIS BUYS, STATED RATHER THAN HIDDEN**: the pseudostate ladder is FROZEN, built
    once in a fixed local potential, never relaxed or EOL-optimized against the full combined CI. This is the
    same kind of approximation a DFS-based CI already makes elsewhere in JAC, not a new one invented here --
    but it IS an approximation, and a future round that compares against a relaxed treatment should expect some
    residual from it.

    THE NORMALIZATION: `crossSectionAt` turns the discrete "stick spectrum" of dipole strengths into a
    continuous cross section via the density-of-states recipe `sigma(E_k) = dipoleStrength_k / localSpacing_k`,
    MULTIPLIED by the SAME physical prefactor `PhotoIonization.computeCrossSection` applies to an energy-
    normalized continuum amplitude (`8*pi^3/(alpha*omega*gi)`) -- the prefactor is NOT optional and was missing
    from the first working version: Phase 0's toy validation only ever checked a RATIO (the Fano lineshape
    shape, and a bare box against its own defining formula), in which any constant prefactor cancels and so was
    never caught. `examples/example-Dz.jl` branch a's own PI-limit check is what found this, giving ratios of
    order 1e-4 before the fix and order 1 after it.

    **A SECOND LIMITATION, FOUND THE SAME WAY: the pseudostate ladder used to share the SAME grid the bound
    SCF used, with no mechanism to give it its OWN box -- PARTLY FIXED 05-Nov-2026, AND THE FIX EXPOSES A
    DEEPER, QUANTIFIED COMPUTATIONAL CEILING RATHER THAN SOLVING THE PROBLEM OUTRIGHT.** `settings.pseudostateRbox`
    now lets the ladder's own grid EXTEND the caller's (`ladderGridFor`), verified to need no basis
    reconciliation at all: `Radial.Grid(gr; rbox=...)` reproduces the smaller grid's own mesh points EXACTLY
    (checked numerically, max|diff| = 0.0), so a bound orbital's existing array is already a valid, exact
    prefix of the larger grid and ordinary zero-padding (already how every radial-integral routine here treats
    differing array lengths) is all that is needed.

    **MEASURED, NOT ASSUMED, WHAT THIS ACTUALLY BUYS.** Fe K-LL (branch b): going from rbox=4 (Fe's own bound
    box, ZERO pseudostates in a 5 Ha window at the resonance's ~171 Ha electron energy) to rbox=300 (nsL~2200,
    a dense eigensolve of several thousand that still completes) gives REAL pseudostates there -- but the cross
    section comes out ~22x smaller than `example-Dy.jl`'s own validated peak, because the LOCAL spacing even at
    rbox=300 (~0.19 Ha, estimated from the free-particle law pi*sqrt(2E)/rbox) is still far wider than the
    resonance's own physical width. He 2s2p (branch c): rbox=200 (nsL~2900) still returns a smooth, structureless
    RAMP, not a lineshape -- reaching Gamma = 38 meV-level local spacing there would need rbox of order 17000
    a.u. by the same estimate, which would need a dense eigensolve of a basis roughly 10x larger again, a cubic-
    scaling cost increase that is not a reasonable next box size to simply try.

    **06-Nov-2026 -- INSTRUMENTING THE REAL STICK STATES DIRECTLY (NOT THE PRINTED SCAN CURVE) SHOWS THE ABOVE
    TWO PARAGRAPHS WERE PRECISE ABOUT THE SYMPTOM BUT IMPRECISE ABOUT WHAT A "SCAN" ACTUALLY MEASURES, AND THE
    CORRECTED PICTURE FOUND A REAL, ENCOURAGING SIGNAL.** A printed scan point is one `crossSectionAt` call,
    which LINEARLY INTERPOLATES between whichever two real stick states happen to bracket that one energy -- it
    is not independent physics, however many points are printed. Checked directly against `Spectrum.states`:
    Fe's "~22x too small at resonance" is a SINGLE interpolated point between two sticks 2.86 eV apart with
    ZERO real sticks within +-0.3 eV of the resonance -- a background value, not a diluted resonance. He DOES
    have two real sticks inside its scanned 0.6 eV window, and they already show the right physics: one
    (omega=56.741 eV) is background-like (dipole strength 3.09e-4, 1.6% weight on the 2s2p reference); the other
    (omega=56.889 eV, almost exactly the literature 56.8856 eV) carries dipole strength 6.34e-6 -- 50x SMALLER --
    with 98% weight on the 2s2p reference. A sharp dipole-strength suppression exactly at the state carrying the
    resonance's own character is the textbook Fano destructive-interference signature. **The "smooth ramp" was
    the interpolated curve through 2-3 real points; the real points underneath it already show the effect, just
    far too sparsely (one point per ~0.35-0.4 eV against a 38 meV width) to trace a lineshape or extract q.**

    **06-Nov-2026 -- A SECOND CHECK, PROMPTED DIRECTLY: IS THE FROZEN ION CORE ITSELF A GOOD SINGLE CSF?**
    `buildPseudostateLadder` uses only `finalIonLevel`'s DOMINANT CSF's occupation to build every pseudostate
    configuration (see `computeSpectrum`'s own `coreCsf`/`coreOccupation`), discarding any other CSF the ion
    might carry -- a real approximation IF the ion has its own configuration mixing, separate from the box-
    resolution problem above. Tested directly for the Fe Li-like ion ("1s^2 2s", used in branch b): adding a
    realistic correlating partner (`1s 2p^2`, same frozen 1s/2s/2p orbitals, same J=1/2 even symmetry) and
    re-diagonalizing gives a ground-level weight of 0.999995 on the pure "1s^2 2s" CSF (0.0005% mixing, 0.03 eV
    energy shift). For He's "1s" ion core this question doesn't even arise -- a one-electron ion cannot mix.
    So, against the one plausible partner tried, the single-CSF ion-core approximation is MEASURED accurate
    here, not merely assumed for convenience; this is not exhaustive (only one correlating configuration was
    tried), and a different final-ion system should not inherit this result without its own check.

    **06-Nov-2026 -- A THIRD, STRUCTURAL LIMITATION, FOUND BY A DIRECT QUESTION ABOUT WHY `computeSpectra` RAN
    ONE DIAGONALIZATION PER KAPPA RATHER THAN ONE DIAGONALIZATION TOTAL -- FIXED, SAME DAY.** The atomic
    Hamiltonian conserves total J and parity exactly, so CSFs of DIFFERENT (J,parity) are exactly block-diagonal
    and splitting those costs nothing. But checked directly via `PhotoIonization.determineChannels` for both
    test systems: kappa=-2 (p3/2) and kappa=+1 (p1/2) BOTH reach the SAME symmetry, `LevelSymmetry(1,-)`, for
    this E1 process from a J=0 even initial level -- the GENERIC case for single-photon ionization, since the
    dipole selection rule fixes the final (J,parity) independent of which partial wave carries it. The earlier
    version of `computeSpectra` built a SEPARATE `Spectrum` (a separate `Hamiltonian.performCI` call) per
    `(finalIonLevel, kappa)` pair and simply ADDED the resulting cross sections afterward -- dropping the real
    Hamiltonian matrix elements connecting a p1/2-built CSF to a p3/2-built CSF of that SAME (J=1,-) block, an
    approximation rather than an exact symmetry decomposition.

    **FIXED, 06-Nov-2026: `computeSpectra` now inverts `determineChannels`' own (kappa -> symmetries) map into
    (symmetry -> kappas), and `computeSpectrum` takes the FULL SET of kappas reaching one symmetry, pools every
    one of their pseudostate ladders into ONE `combinedBasis`, and diagonalizes ONCE** (`Spectrum.kappas`,
    `Spectrum.symmetry` replace the old scalar `kappa` field; `StickState.weightByKappa` is the new, genuinely
    useful by-product -- each eigenstate's own weight broken down by WHICH partial wave's pseudostate CSFs it
    sits on, computed for free from the CSF index ranges already tracked while building `combinedCsfs`). The
    single-kappa case (e.g. bare hydrogen, where kappa alone fixes (J,parity) one-to-one) is simply
    `length(kappas) == 1`, not a separately-maintained code path -- confirmed by re-running `example-Dz.jl`
    branch a: ratio unchanged at 1.0068, byte-identical to before this fix.

    **MEASURED, NOT ASSUMED, WHAT COUPLING ACTUALLY CHANGES.** Fe K-LL (branch b), at the resonance energy:
    coupled sigma(Coulomb) = 2.930310e-02, against the OLD split-and-add sum of 6.003433e-02 -- roughly HALVED
    by destructive interference between the two partial waves, at a point with zero real sticks nearby (pure
    background; the resolution gap to `example-Dy.jl`'s 1.336798 peak remains, but is now attributable entirely
    to pseudostate DENSITY, not to this approximation). He 2s2p (branch c) is the striking case: **the coupled
    scan shows a SHARP, asymmetric feature that was completely ABSENT from the old split-and-add scan at the
    SAME box and density** -- rising to sigma=0.133 at omega=56.7656 eV, then falling to a near-total
    destructive ZERO, sigma=2.5e-7, at omega=56.9056 eV, close to the resonance's own 56.8856 eV. Coupling the
    partial waves was not merely a quantitative correction here; it is what makes the resonance feature VISIBLE
    AT ALL at this density. The scan is still far too coarse (one point per 0.02 eV against Gamma ~ 38 meV) to
    read off q or Gamma -- that is the separate, still-open density problem below.

    **SO THE RECOMMENDED NEXT STEPS, IN ORDER (held at challenges-list S29):**
    1. ~~Couple same-symmetry kappas into one diagonalization~~ -- DONE, 06-Nov-2026, above.
    2. **Raise the local pseudostate DENSITY near one target energy** via a sparse/shift-invert eigensolver
       targeting a narrow window (cost scales with the window's own state count, not the whole basis) rather
       than a bigger dense box -- reuses the existing B-spline primitive basis with a different, targeted solve,
       not new theory. Brute-force box enlargement under DENSE diagonalization has a genuine, quantified
       computational ceiling (rbox~300 for Fe, ~200 for He already near-practical; reaching either resonance's
       own width would need an order of magnitude more, at cubic cost, now compounded by step 1's own larger
       pooled basis) -- NOT a reasonable next box size to simply try.
    3. **Re-validate the ion-core single-CSF approximation (point above) for any NEW final-ion system** before
       trusting its pseudostate ladder -- cheap to check (one extra `performCIwithFrozenOrbitals` call with a
       plausible correlating configuration), not yet shown to generalize beyond the two cases tested.
    Alternative, not attempted: genuine R-matrix outer-region propagation, which handles the asymptotic
    continuum analytically rather than by sampling it densely -- a larger, separate undertaking than 2-3 above.

    VALIDATED, before any `src/` code existed: two independent toy cases (see `work/photoionization-pseudostates-
    toy/`, Phase 0, 05-Nov-2026) -- a closed-form Fano two-level-plus-flat-continuum model (exact in the trivial
    uncoupled limit, converging toward the exact Fano lineshape as the box widens at fixed pseudostate density)
    and a direct check that `Bsplines.diagonalizeLocalMatrix`'s own positive-energy branch reproduces the exact
    free-particle box density of states to 1e-3..1e-4 relative accuracy, improving with box size. Both toys
    validated the SHAPE of the normalization, which is why neither caught the missing prefactor above.

    **A THIRD THING FOUND AND FIXED THE SAME WAY: `crossSectionAt` LINEARLY INTERPOLATES between the two
    sticks bracketing the requested energy, rather than reading off the single nearest one.** Measured against
    hydrogen's own EXACT Stobbe cross section (one-electron, no SCF or frozen-core approximation anywhere, so
    nothing else could be responsible) -- the first, nearest-neighbour version gave ratio = 1.0636 (6.4% high).
    Printing the individual stick-state sigma_k showed them SMOOTH and monotonically falling with energy, not
    noisy, so the error was never statistical scatter to be averaged away: it was evaluating at the nearest
    AVAILABLE stick's own energy rather than the one actually requested, on a cross section that falls off
    steeply here. Interpolating between the two bracketing sticks closed the gap to ratio = 1.0068 (0.7%).
    This is a DIFFERENT fix from the missing-prefactor one above, found the same way (a real example, not the
    toys) and independent of it. It does NOT touch the second limitation (the shared-grid one): a window with
    zero pseudostates in it, or one far coarser than Gamma, has nothing to interpolate between either.
"""
module PhotoIonizationPseudostates


using  Printf, SparseArrays, ..Basics, ..Bsplines, ..Continuum, ..Defaults, ..Hamiltonian, ..ManyElectron,
       ..Nuclear, ..PhotoEmission, ..PhotoIonization, ..Radial, ..RadialIntegrals, ..TableStrings


"""
`struct  PhotoIonizationPseudostates.Settings  <:  AbstractProcessSettings`
    ... defines a type for the details and parameters of the diagonalize-once, pseudostate-continuum treatment
        of resonant photoionization.

    + multipoles               ::Array{EmMultipole,1}   ... Multipoles of the absorbed photon that are taken into account.
    + gauges                   ::Array{UseGauge,1}      ... Gauges to be included into the computations.
    + photonEnergies           ::Array{Float64,1}       ... Photon energies [a.u.] at which the SMOOTHED, density-of-states
                                                            normalized cross section is evaluated -- ONE diagonalization
                                                            per (initial level, final-ion level, kappa) serves the WHOLE
                                                            scan, unlike `PhotoIonizationInterference`'s per-energy solve.
    + lValues                  ::Array{Int64,1}         ... Orbital angular momenta of the outgoing electron that are
                                                            considered, exactly as in `PhotoIonization.Settings`.
    + pseudostateWindow        ::Tuple{Float64,Float64}  ... Energy window [a.u.], in the ejected electron's OWN one-
                                                            electron energy (i.e. directly comparable to photonEnergy minus
                                                            the ionization energy), of box-discretized pseudostates to
                                                            retain per kappa. Should comfortably bracket the electron-energy
                                                            range `photonEnergies` implies, with margin on both sides --
                                                            `crossSectionAt` is only trustworthy away from the window's own
                                                            edges, exactly as the Phase-0 toy validations were restricted.
    + pseudostateRbox          ::Union{Nothing,Float64} ... Box radius [a.u.] for the pseudostate ladder's OWN grid,
                                                            built by EXTENDING the caller's own grid (same rnt/h/hp,
                                                            larger rbox only) -- see `buildPseudostateLadder`'s own
                                                            docstring for why this is exact, not an approximation, over
                                                            the smaller grid's own range. `nothing` -- the default --
                                                            reuses the caller's grid unchanged (its own rbox), which is
                                                            adequate only near threshold; a resonance far above
                                                            threshold or narrower than that grid's local level spacing
                                                            needs this set explicitly, larger than the bound SCF's own
                                                            box ever needed to be.
    + continuumExchange        ::Union{Nothing,Bool}    ... True/false to add the frozen nonlocal Dirac-Fock exchange
                                                            operator (`Bsplines.setupExchangeMatrix`, item 99's mechanism)
                                                            to the pseudostate ladder's own generating potential;
                                                            `nothing` -- the default -- resolves to `false`. Unlike other
                                                            process modules this does NOT resolve automatically from a
                                                            `finalAsfSettings.scField`, since this entry point is handed a
                                                            `Level` rather than an `AsfSettings`.
    + printBefore              ::Bool                   ... True, if the combined-basis CI is printed before diagonalization.
"""
struct Settings  <:  AbstractProcessSettings
    multipoles                 ::Array{EmMultipole,1}
    gauges                      ::Array{UseGauge,1}
    photonEnergies               ::Array{Float64,1}
    lValues                       ::Array{Int64,1}
    pseudostateWindow              ::Tuple{Float64,Float64}
    pseudostateRbox                  ::Union{Nothing,Float64}
    continuumExchange                  ::Union{Nothing,Bool}
    printBefore                          ::Bool
end


"""
`PhotoIonizationPseudostates.Settings()`
    ... constructor for an `empty` instance of PhotoIonizationPseudostates.Settings; a settings::Settings is returned.
"""
function Settings()
    Settings( EmMultipole[E1], UseGauge[UseCoulomb], Float64[], [0,1,2,3,4,5], (0., 1.), nothing, nothing, false )
end


"""
`PhotoIonizationPseudostates.Settings(set::PhotoIonizationPseudostates.Settings;
                                      multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,
                                      gauges::Union{Nothing,Array{UseGauge,1}}=nothing,
                                      photonEnergies::Union{Nothing,Array{Float64,1}}=nothing,
                                      lValues::Union{Nothing,Array{Int64,1}}=nothing,
                                      pseudostateWindow::Union{Nothing,Tuple{Float64,Float64}}=nothing,
                                      pseudostateRbox::Union{Nothing,Float64}=nothing,
                                      continuumExchange::Union{Nothing,Bool}=nothing,
                                      printBefore::Union{Nothing,Bool}=nothing)`

        multipoles=..,          gauges=..,              photonEnergies=..,       lValues=..,
        pseudostateWindow=..,   pseudostateRbox=..,     continuumExchange=..,   printBefore=..)

    ... keyword copy-constructor for re-defining selected values of a settings::PhotoIonizationPseudostates.Settings; a
        settings::Settings is returned.
"""
function Settings(set::PhotoIonizationPseudostates.Settings;
    multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,        gauges::Union{Nothing,Array{UseGauge,1}}=nothing,
    photonEnergies::Union{Nothing,Array{Float64,1}}=nothing,        lValues::Union{Nothing,Array{Int64,1}}=nothing,
    pseudostateWindow::Union{Nothing,Tuple{Float64,Float64}}=nothing,
    pseudostateRbox::Union{Nothing,Float64}=nothing,                 continuumExchange::Union{Nothing,Bool}=nothing,
    printBefore::Union{Nothing,Bool}=nothing)

    if  isnothing(multipoles)          multipolesx         = set.multipoles         else   multipolesx         = multipoles         end
    if  isnothing(gauges)              gaugesx             = set.gauges             else   gaugesx             = gauges             end
    if  isnothing(photonEnergies)      photonEnergiesx     = set.photonEnergies     else   photonEnergiesx     = photonEnergies     end
    if  isnothing(lValues)             lValuesx            = set.lValues            else   lValuesx            = lValues            end
    if  isnothing(pseudostateWindow)   pseudostateWindowx  = set.pseudostateWindow  else   pseudostateWindowx  = pseudostateWindow  end
    if  isnothing(pseudostateRbox)     pseudostateRboxx    = set.pseudostateRbox    else   pseudostateRboxx    = pseudostateRbox    end
    if  isnothing(continuumExchange)   continuumExchangex  = set.continuumExchange  else   continuumExchangex  = continuumExchange  end
    if  isnothing(printBefore)         printBeforex        = set.printBefore        else   printBeforex        = printBefore        end

    Settings( multipolesx, gaugesx, photonEnergiesx, lValuesx, pseudostateWindowx, pseudostateRboxx, continuumExchangex,
              printBeforex )
end


# `Base.show(io::IO, settings::PhotoIonizationPseudostates.Settings)`  ... prepares a proper printout of settings.
function Base.show(io::IO, settings::PhotoIonizationPseudostates.Settings)
    println(io, "multipoles:               $(settings.multipoles)  ")
    println(io, "gauges:                   $(settings.gauges)  ")
    println(io, "photonEnergies:           $(settings.photonEnergies)  ")
    println(io, "lValues:                  $(settings.lValues)  ")
    println(io, "pseudostateWindow:        $(settings.pseudostateWindow)  ")
    println(io, "pseudostateRbox:          $(settings.pseudostateRbox)  ")
    println(io, "continuumExchange:        $(settings.continuumExchange)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
end


"""
`struct  PhotoIonizationPseudostates.StickState`
    ... defines a type for ONE eigenstate of the combined (bound + pseudostate-ladder) diagonalization.

    + energy                   ::Float64                ... Total energy of this eigenstate [a.u.].
    + dipoleStrength            ::EmProperty             ... |<this state|D|initialLevel>|^2, both gauges, summed
                                                            incoherently over the requested multipoles.
    + localSpacing               ::Float64                ... (E_{k+1}-E_{k-1})/2 -- the density-of-states normalization
                                                            divisor `crossSectionAt` uses.
    + weightOnReference            ::Float64                ... |c|^2 summed over the CSFs this state shares with
                                                            `initialLevel.basis` (i.e. everything that is NOT a newly
                                                            added pseudostate-ladder CSF) -- a resonant-vs-background
                                                            diagnostic: close to 1 marks a state dominated by the bound,
                                                            autoionizing character, close to 0 a smooth background state.
    + weightByKappa                ::Dict{Int64,Float64}    ... |c|^2 summed over the pseudostate-ladder CSFs of EACH
                                                            contributing partial wave separately (06-Nov-2026, added
                                                            when same-symmetry kappas were pooled into one combined
                                                            diagonalization rather than diagonalized apart and added --
                                                            see `Spectrum`'s own docstring) -- a genuine partial-channel
                                                            diagnostic: how much of this eigenstate's own character sits
                                                            in each outgoing partial wave, not merely whether it is
                                                            bound-like or continuum-like.
"""
struct StickState
    energy                      ::Float64
    dipoleStrength                ::EmProperty
    localSpacing                    ::Float64
    weightOnReference                 ::Float64
    weightByKappa                       ::Dict{Int64,Float64}
end


"""
`PhotoIonizationPseudostates.StickState()`
    ... constructor for an `empty` instance of PhotoIonizationPseudostates.StickState; a state::StickState is returned.
"""
function StickState()
    StickState( 0., EmProperty(0., 0.), 0., 0., Dict{Int64,Float64}() )
end


# `Base.show(io::IO, state::PhotoIonizationPseudostates.StickState)`  ... prepares a proper printout of state.
function Base.show(io::IO, state::PhotoIonizationPseudostates.StickState)
    println(io, "energy:                   $(state.energy)  ")
    println(io, "dipoleStrength:           $(state.dipoleStrength)  ")
    println(io, "localSpacing:             $(state.localSpacing)  ")
    println(io, "weightOnReference:        $(state.weightOnReference)  ")
    println(io, "weightByKappa:            $(state.weightByKappa)  ")
end


"""
`struct  PhotoIonizationPseudostates.Spectrum`
    ... defines a type for the full stick spectrum of ONE combined diagonalization, for ONE (initial level,
        final-ion level, symmetry) combination -- spanning EVERY partial wave that reaches that symmetry,
        pooled into the SAME `Hamiltonian.performCI` call.

        WHY POOLED BY SYMMETRY, NOT BY KAPPA (06-Nov-2026) -- the atomic Hamiltonian conserves total J and
        parity exactly, so CSFs of DIFFERENT (J,parity) are exactly block-diagonal and splitting those costs
        nothing; but two DIFFERENT partial waves of the SAME (J,parity) are NOT exactly block-diagonal from
        each other -- their mutual Hamiltonian matrix elements are, in general, nonzero. Checked directly for
        both test systems in `example-Dz.jl` (He 2s2p, Fe K-LL): kappa=-2 (p3/2) and kappa=+1 (p1/2) both reach
        the SAME symmetry there, the generic case for single-photon ionization since the dipole selection rule
        fixes the final (J,parity) independent of which partial wave carries it. An earlier version of this
        module diagonalized each kappa SEPARATELY and added the resulting cross sections afterward -- dropping
        those real matrix elements, an approximation rather than an exact decomposition. `computeSpectrum`
        now takes the FULL SET of kappas reaching one symmetry and diagonalizes them together; a system where
        only one kappa ever reaches a given symmetry (e.g. bare hydrogen, where kappa alone fixes (J,parity)
        one-to-one) is simply the special case `length(kappas) == 1`, not a separate code path.

    + symmetry                   ::LevelSymmetry          ... Total (J,parity) this spectrum's combined CI was
                                                            restricted to.
    + kappas                       ::Array{Int64,1}         ... Every partial wave of the ejected electron pooled into
                                                            this ONE diagonalization (all reaching `symmetry`).
    + initialLevel                   ::Level                  ... Initial (bound, possibly correlated) level.
    + finalIonLevel                     ::Level                  ... Final-ion level, one electron fewer than initialLevel.
    + states                               ::Array{StickState,1}   ... Every eigenstate of the combined diagonalization,
                                                            sorted by energy.
"""
struct Spectrum
    symmetry                      ::LevelSymmetry
    kappas                           ::Array{Int64,1}
    initialLevel                        ::Level
    finalIonLevel                          ::Level
    states                                    ::Array{StickState,1}
end


"""
`PhotoIonizationPseudostates.Spectrum()`
    ... constructor for an `empty` instance of PhotoIonizationPseudostates.Spectrum; a spectrum::Spectrum is returned.
"""
function Spectrum()
    Spectrum( LevelSymmetry(AngularJ64(0), Basics.plus), Int64[], Level(), Level(), StickState[] )
end


# `Base.show(io::IO, spectrum::PhotoIonizationPseudostates.Spectrum)`  ... prepares a proper printout of spectrum.
function Base.show(io::IO, spectrum::PhotoIonizationPseudostates.Spectrum)
    println(io, "symmetry:                 $(spectrum.symmetry)  ")
    println(io, "kappas:                   $(spectrum.kappas)  ")
    println(io, "initialLevel:             $(spectrum.initialLevel)  ")
    println(io, "finalIonLevel:            $(spectrum.finalIonLevel)  ")
    println(io, "states:                   $(length(spectrum.states)) stick state(s)  ")
end


"""
`PhotoIonizationPseudostates.crossSectionAt(spectrum::PhotoIonizationPseudostates.Spectrum, omega::Float64)`
    ... to evaluate the density-of-states-normalized cross section at photon energy omega [a.u.]. Each stick
        state k already carries `sigma_k = dipoleStrength_k / localSpacing_k`, the normalization this module's
        own docstring describes and Phase 0 validated; this function then LINEARLY INTERPOLATES sigma_k between
        the two sticks that bracket `spectrum.initialLevel.energy + omega`, rather than reading off the single
        nearest one. A value::EmProperty is returned.

        **WHY INTERPOLATION, NOT NEAREST-NEIGHBOUR -- FOUND 05-Nov-2026, read before reverting this.** The
        first working version used nearest-neighbour and gave a stable ~6-8% discrepancy against hydrogen's
        EXACT Stobbe cross section (`examples/example-Dz.jl` branch a). Printing the individual sigma_k around
        the evaluation point showed them to be SMOOTH and monotonic, not noisy -- so the discrepancy was never
        statistical scatter to be averaged away, it was simply evaluating at the nearest AVAILABLE stick's own
        energy instead of the requested one, on a function (hydrogen's photoionization cross section) that
        falls off steeply with energy. A plain linear interpolation between the two bracketing sticks closed
        the hydrogen residual from 6.4% to under 1%. At the array's own edges, where no upper or lower
        bracketing stick exists, this falls back to the nearest endpoint.
"""
function crossSectionAt(spectrum::PhotoIonizationPseudostates.Spectrum, omega::Float64)
    length(spectrum.states) == 0   &&   error("PhotoIonizationPseudostates.crossSectionAt(): spectrum carries no stick states.")
    totalEnergy = spectrum.initialLevel.energy + omega
    states = spectrum.states                                    # already sorted by energy (computeSpectrum)
    sigmaOf(st) = EmProperty(st.dipoleStrength.Coulomb/st.localSpacing, st.dipoleStrength.Babushkin/st.localSpacing)

    if      totalEnergy <= states[1].energy       return( sigmaOf(states[1]) )
    elseif  totalEnergy >= states[end].energy     return( sigmaOf(states[end]) )
    end
    hi = findfirst(st -> st.energy >= totalEnergy, states);   lo = hi - 1
    frac = (totalEnergy - states[lo].energy) / (states[hi].energy - states[lo].energy)
    sigmaLo = sigmaOf(states[lo]);   sigmaHi = sigmaOf(states[hi])

    return( EmProperty( sigmaLo.Coulomb   + frac*(sigmaHi.Coulomb   - sigmaLo.Coulomb),
                        sigmaLo.Babushkin + frac*(sigmaHi.Babushkin - sigmaLo.Babushkin) ) )
end


"""
`PhotoIonizationPseudostates.removeOccupiedComponent(orb::Radial.Orbital, occupiedOrbitals::Dict{Subshell,Orbital},
                                                      grid::Radial.Grid)`
    ... removes from orb every component along an OCCUPIED, BOUND orbital of orb's own kappa (mirroring
        `RandomPhase.projectOntoUnoccupied`'s own, already-validated algorithm, written locally rather than
        reused since that function's input type is `DalgarnoLewis.PerturbedOrbital`, not a plain
        `Radial.Orbital`), and renormalizes the result to unit norm in this grid's quadrature. A core electron
        cannot be excited into a filled orbital, so without this step a pseudostate whose kappa coincides with
        an occupied ion-core shell would not be orthogonal to it. An orbital::Orbital is returned.
"""
function removeOccupiedComponent(orb::Radial.Orbital, occupiedOrbitals::Dict{Subshell,Orbital}, grid::Radial.Grid)
    P = copy(orb.P);   Q = copy(orb.Q)
    for  (sh, occ)  in  occupiedOrbitals
        (sh.kappa == orb.subshell.kappa  &&  occ.energy < 0.0)   ||   continue
        ov = RadialIntegrals.overlap(occ.P, P, grid) + RadialIntegrals.overlap(occ.Q, Q, grid)
        nn = RadialIntegrals.overlap(occ.P, occ.P, grid) + RadialIntegrals.overlap(occ.Q, occ.Q, grid)
        nn <= 0.0   &&   continue
        for  i = 1:min(length(P), length(occ.P))   P[i] = P[i] - ov/nn * occ.P[i]   end
        for  i = 1:min(length(Q), length(occ.Q))   Q[i] = Q[i] - ov/nn * occ.Q[i]   end
    end
    norm = RadialIntegrals.overlap(P, P, grid) + RadialIntegrals.overlap(Q, Q, grid)
    norm <= 0.0   &&   error("PhotoIonizationPseudostates.removeOccupiedComponent(): a pseudostate was entirely " *
                              "consumed by the occupied-orbital projection -- the box or the window is too small.")
    P = P / sqrt(norm);   Q = Q / sqrt(norm)

    return( Orbital(orb.subshell, orb.isBound, orb.useStandardGrid, orb.energy, P, Q, orb.Pprime, orb.Qprime, orb.grid) )
end


"""
`PhotoIonizationPseudostates.ladderGridFor(grid::Radial.Grid, settings::PhotoIonizationPseudostates.Settings)`
    ... builds the grid the pseudostate ladder itself is generated on -- the caller's own `grid` UNCHANGED if
        `settings.pseudostateRbox` is `nothing` (the default), or that same grid EXTENDED to the larger box via
        `Radial.Grid(grid; rbox=settings.pseudostateRbox)` otherwise.

        **WHY EXTENDING rbox, WITH rnt/h/hp UNCHANGED, NEEDS NO BASIS RECONCILIATION -- VERIFIED NUMERICALLY,
        05-Nov-2026, NOT ASSUMED.** `Radial.Grid(gr; rbox=...)` recomputes only the NUMBER of points needed to
        reach the new box at the SAME point-spacing law; the mesh points inside the OLD box are reproduced
        EXACTLY (checked: two grids built with `hp=0.02` and `rbox=4.0` vs `rbox=40.0` agree on their first 497
        points to `max|diff| = 0.0`). An orbital tabulated on the smaller grid is therefore already a valid,
        exact prefix of the larger grid's own mesh -- extending it is zero-padding, not interpolation, and
        every radial-integral routine already used here (`RadialIntegrals.overlap`, `removeOccupiedComponent`)
        already loops to `min(length(...))`, so nothing further is needed to mix the two. A `grid::Radial.Grid`
        is returned.
"""
function ladderGridFor(grid::Radial.Grid, settings::PhotoIonizationPseudostates.Settings)
    isnothing(settings.pseudostateRbox)   &&   return( grid )

    return( Radial.Grid(grid; rbox=settings.pseudostateRbox) )
end


"""
`PhotoIonizationPseudostates.solveWindowBySegments(kappa::Int64, matrixA::AbstractMatrix{Float64},
        matrixB::AbstractMatrix{Float64}, primitives::Bsplines.Primitives, window::Tuple{Float64,Float64})`
    ... covers the WHOLE energy `window` with real eigenpairs via several small shift-invert calls
        (`Bsplines.diagonalizeLocalMatrixNearEnergy`) rather than one call asking for however many eigenpairs
        the whole window needs -- Arnoldi/Lanczos converges best for eigenvalues NEAR the shift, so a wide
        window is split into `nSegments` equal-width sub-ranges, each solved by its OWN shift-invert call
        centred on its own midpoint. `nSegments` is estimated from the free-particle local-spacing formula
        `pi*sqrt(2E)/rbox` (the same estimate used throughout `Bsplines`' own development today) at the
        window's own midpoint, targeting `targetStatesPerSegment` states per segment -- a heuristic, not a
        derived bound -- and NATURALLY degenerates to 1 for a window already narrow relative to the local
        spacing, so there is no separate code path for that common, near-threshold case. Each segment's own
        `nev` is estimated the SAME way from its OWN midpoint, since spacing itself varies with E, and
        doubled, up to 5 times, if the returned eigenvalues do not cover the segment -- fewer than 3 kept, or
        the kept range falls short of either edge by more than one local spacing. A
        `(values::Array{Float64,1}, vectors::Array{Vector{Float64},1})` pair, sorted by value and covering the
        whole window, is returned.
"""
function solveWindowBySegments(kappa::Int64, matrixA::AbstractMatrix{Float64}, matrixB::AbstractMatrix{Float64},
                                primitives::Bsplines.Primitives, window::Tuple{Float64,Float64})
    windowLo, windowHi     = window
    rbox                   = primitives.grid.r[end]
    targetStatesPerSegment = 30

    midE         = max( (windowLo+windowHi)/2, 1.0e-6 )
    spacingAtMid = pi * sqrt(2*midE) / rbox
    nSegments    = max( 1, ceil(Int64, (windowHi-windowLo) / (targetStatesPerSegment*spacingAtMid)) )
    edges        = range(windowLo, windowHi, length=nSegments+1)

    allValues  = Float64[];    allVectors = Vector{Float64}[]
    for  seg = 1:nSegments
        eLo = edges[seg];   eHi = edges[seg+1]
        sigma        = (eLo+eHi)/2
        localSpacing = pi * sqrt(2*max(sigma,1.0e-6)) / rbox
        nev          = max( 6, ceil(Int64, 1.5*(eHi-eLo)/localSpacing) )

        attempt = 0
        while true
            eig        = Bsplines.diagonalizeLocalMatrixNearEnergy(kappa, matrixA, matrixB, primitives, sigma, nev)
            keepRightEdge = seg == nSegments
            keptIdx    = findall( E -> (eLo <= E < eHi)  ||  (keepRightEdge && E == eHi),  eig.values )
            covered    = length(keptIdx) >= 3   &&   minimum(eig.values[keptIdx]) <= eLo + localSpacing   &&
                                                       maximum(eig.values[keptIdx]) >= eHi - localSpacing
            if  covered
                append!(allValues, eig.values[keptIdx]);   append!(allVectors, eig.vectors[keptIdx])
                break
            end
            attempt = attempt + 1
            attempt > 5   &&   error("PhotoIonizationPseudostates.solveWindowBySegments(): segment [$eLo, " *
                                      "$eHi) for kappa = $kappa did not reach full coverage after $attempt " *
                                      "attempts (last nev = $nev) -- check pseudostateRbox and " *
                                      "pseudostateWindow units, or that the box is matched to this energy range.")
            nev = nev * 2
        end
    end

    order = sortperm(allValues)

    return( allValues[order], allVectors[order] )
end


"""
`PhotoIonizationPseudostates.buildPseudostateLadder(kappa::Int64, finalIonLevel::Level, nm::Nuclear.Model,
                                                      grid::Radial.Grid, boundGrid::Radial.Grid,
                                                      settings::PhotoIonizationPseudostates.Settings)`
    ... builds the box-discretized continuum pseudostate ladder of the given kappa, in the frozen local
        potential of finalIonLevel (nuclear + `Basics.DFSField(1.0)` mean field, the same combination
        `Continuum.generateOrbitalForLevel` already uses), keeping every eigenstate inside
        settings.pseudostateWindow. Each retained state is labelled `Subshell(1000+k, kappa)` (k = 1, 2, ...,
        purely a bookkeeping label, chosen far above any physical principal quantum number so it can never
        collide with an ordinary bound subshell) and Gram-Schmidt orthogonalized against finalIonLevel's own
        occupied orbitals of the same kappa via `removeOccupiedComponent`. `grid` here is already
        `ladderGridFor`'s own result, not necessarily the caller's original grid -- see that function's
        docstring; `boundGrid` is that ORIGINAL, small grid, needed for one thing only -- see next.

        **WHY THE MEAN FIELD IS BUILT ON `boundGrid`, THEN EXTENDED ONTO `grid`, RATHER THAN BUILT ON `grid`
        DIRECTLY.** `Basics.computePotential(DFSField(...), ...)` integrates via a brute-force O(npoints^2)
        double loop -- harmless on the small `boundGrid` (npoints ~ hundreds) but, measured 06-Nov-2026,
        ~7e11 operations and effectively unaffordable on the large ladder `grid` a resolved resonance needs
        (npoints ~ 850000 at rbox~17000) -- see priority-list item 7, where fixing this PROPERLY (inside
        `Basics` itself) is held as its own, separate, deliberately deferred task. This function instead
        computes the mean field on the cheap, small `boundGrid` (exactly as before today's work), then hands
        `Basics.add` a `Radial.Potential` whose `Zr` array is that short result but whose `grid` field points
        at the LARGE ladder `grid` -- `Basics.add`'s own existing logic (`Zr[nmin+1:nmax] .= Zr[nmin]`,
        `module-BasicsAZ-inc-AG.jl`) already extends a shorter potential as a CONSTANT past its own end. This
        is not a new approximation invented here: it is the EXACT physical behaviour for a mean field from a
        charge density that is itself exactly zero beyond `boundGrid`'s own extent (the core's charge seen
        from far away is just its total enclosed charge, a constant, no further integration needed) -- it is
        the existing convention `Basics.add` was already written to exploit, not a shortcut specific to this
        module. `boundGrid` and `grid` share the same `rnt`/`h`/`hp` by construction of `ladderGridFor`, so
        `Basics.add`'s own grid-consistency check passes exactly.

        **TWO PATHS, SPLIT ON `settings.continuumExchange`, NOT ONE.** Exchange (`Bsplines.setupExchangeMatrix`)
        is genuinely NONLOCAL, so it cannot be represented by the sparse/banded assembly the default path
        uses -- `Bsplines.diagonalizeLocalMatrixNearEnergy`'s own guard would refuse it anyway, for the right
        reason. When `continuumExchange == true` this function falls back to the ORIGINAL dense
        `generateTTpMatrix!`/`setupLocalMatrix`/`diagonalizeLocalMatrix` route, filtering the full positive-
        energy branch to the window exactly as before -- that is the only path that can carry exchange, not
        a kept-around fallback for its own sake. Otherwise (the common, default case) it builds the matrices
        directly as sparse (`Bsplines.generateTTpMatrixSparse!`/`setupLocalMatrixSparse`) and covers the
        window via `PhotoIonizationPseudostates.solveWindowBySegments`'s segmented shift-invert solve, rather
        than one dense diagonalization of the whole positive-energy branch -- this is what makes a box
        genuinely matched to a resonance's own width (rbox in the thousands, not capped near rbox~300 by a
        dense diagonalization's memory cost) practical.
"""
function buildPseudostateLadder(kappa::Int64, finalIonLevel::Level, nm::Nuclear.Model, grid::Radial.Grid,
                                boundGrid::Radial.Grid, settings::PhotoIonizationPseudostates.Settings)
    nuclearPot    = Nuclear.nuclearPotential(nm, grid)
    meanPotSmall  = Basics.computePotential(Basics.DFSField(1.0), boundGrid, finalIonLevel)
    meanPot       = Radial.Potential(meanPotSmall.name, meanPotSmall.Zr, grid)
    pot           = Basics.add(nuclearPot, meanPot)
    windowLo, windowHi = settings.pseudostateWindow

    if  settings.continuumExchange == true
        # The ONLY path that can represent exchange -- see this function's own docstring. Unchanged from
        # before today's sparse/shift-invert work.
        primitives = Bsplines.generatePrimitives(grid)
        storage    = Dict{String,Array{Float64,2}}()
        nsL        = primitives.grid.nsL;   nsS = primitives.grid.nsS
        overlap    = zeros(nsL+nsS, nsL+nsS)
        overlap[1:nsL,1:nsL]                 = Bsplines.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
        overlap[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = Bsplines.generateTTpMatrix!("SS-overlap", 0, primitives, storage)

        localMatrix     = Bsplines.setupLocalMatrix(kappa, primitives, pot, storage)
        coreOccupations = Basics.extractMeanOccupation(finalIonLevel.basis)
        localMatrix     = localMatrix + Bsplines.setupExchangeMatrix(kappa, primitives, finalIonLevel.basis.orbitals,
                                                                      coreOccupations)
        eig = Bsplines.diagonalizeLocalMatrix(kappa, localMatrix, overlap, primitives)

        idxWindow = findall(E -> windowLo <= E <= windowHi, eig.values)
        length(idxWindow) < 3   &&   error("PhotoIonizationPseudostates.buildPseudostateLadder(): only " *
                                            "$(length(idxWindow)) pseudostate(s) fall inside pseudostateWindow = " *
                                            "$(settings.pseudostateWindow) for kappa = $kappa -- widen the window, " *
                                            "set settings.pseudostateRbox larger, or check that the window is given " *
                                            "in the ejected electron's OWN energy units.")
        values  = eig.values[idxWindow];    vectors = eig.vectors[idxWindow]
    else
        # The default path, sparse/shift-invert throughout -- see this function's own docstring.
        primitives = Bsplines.generatePrimitives(grid)
        storage    = Dict{String,SparseMatrixCSC{Float64,Int64}}()
        nsL        = primitives.grid.nsL;   nsS = primitives.grid.nsS
        oLL        = Bsplines.generateTTpMatrixSparse!("LL-overlap", 0, primitives, storage)
        oSS        = Bsplines.generateTTpMatrixSparse!("SS-overlap", 0, primitives, storage)
        overlap    = [ oLL  spzeros(nsL,nsS) ;  spzeros(nsS,nsL)  oSS ]

        localMatrix = Bsplines.setupLocalMatrixSparse(kappa, primitives, pot, storage)
        values, vectors = PhotoIonizationPseudostates.solveWindowBySegments(kappa, localMatrix, overlap, primitives,
                                                                             settings.pseudostateWindow)
    end

    ladder = Dict{Subshell,Orbital}()
    mtp    = primitives.grid.NoPoints
    for  (k, (value, vector))  in  enumerate(zip(values, vectors))
        sh  = Subshell(1000+k, kappa)
        raw = Bsplines.generateOrbitalFromPrimitives(sh, value, mtp, vector, primitives)
        ladder[sh] = PhotoIonizationPseudostates.removeOccupiedComponent(raw, finalIonLevel.basis.orbitals, grid)
    end

    return( ladder )
end


"""
`PhotoIonizationPseudostates.checkFrozenOrbitalConsistency(initialLevel::Level, finalIonLevel::Level, grid::Radial.Grid)`
    ... RAISES, rather than silently mixing inconsistent orbitals, if initialLevel.basis and finalIonLevel.basis
        disagree -- overlap below 1.0 - 1.0e-6 -- on any subshell they share. The combined CI this module builds
        assumes ordinary (non-BiOrthogonal) Slater-Condon rules, which require every orbital shared by the two
        sides to be LITERALLY the same radial function; the caller's responsibility is to have generated both
        multiplets against one consistent frozen core. Nothing is returned on success.
"""
function checkFrozenOrbitalConsistency(initialLevel::Level, finalIonLevel::Level, grid::Radial.Grid)
    shared = intersect(keys(initialLevel.basis.orbitals), keys(finalIonLevel.basis.orbitals))
    for  sh  in  shared
        oi = initialLevel.basis.orbitals[sh];    of = finalIonLevel.basis.orbitals[sh]
        ov = RadialIntegrals.overlap(oi.P, of.P, grid) + RadialIntegrals.overlap(oi.Q, of.Q, grid)
        if  ov < 1.0 - 1.0e-6
            error("PhotoIonizationPseudostates.checkFrozenOrbitalConsistency(): initialLevel and finalIonLevel " *
                  "disagree on subshell $sh -- overlap = $ov, expected ~1.0.  The combined CI needs ONE " *
                  "consistent frozen core; re-generate both multiplets against the same orbitals.")
        end
    end

    return( nothing )
end


"""
`PhotoIonizationPseudostates.computeSpectrum(initialLevel::Level, finalIonLevel::Level, kappas::Array{Int64,1},
        symmetry::LevelSymmetry, nm::Nuclear.Model, grid::Radial.Grid,
        settings::PhotoIonizationPseudostates.Settings)`
    ... builds the pseudostate ladder of EVERY kappa in `kappas` (one `buildPseudostateLadder` call per kappa --
        this is the ONLY place that calls it), merges ALL of them together with initialLevel's own CSFs into
        ONE combined basis restricted to the SINGLE `symmetry` (re-expanding initialLevel.basis onto the
        extended subshell list via `ManyElectron.Basis(basis, subshells)`'s existing zero-padding constructor,
        then appending the pseudostate CSFs of every pooled kappa), diagonalizes it ONCE via
        `Hamiltonian.performCI`, and computes the dipole amplitude from initialLevel onto every resulting
        eigenstate via `PhotoEmission.amplitude(Absorption(), ...)` -- which already tolerates the two levels
        carrying different subshell lists/orbital dicts, so no further basis surgery is needed here. A
        spectrum::Spectrum is returned.

        POOLING EVERY KAPPA REACHING `symmetry` INTO ONE DIAGONALIZATION IS THE POINT, NOT AN OPTION -- see
        `Spectrum`'s own docstring for why diagonalizing them apart (an earlier version of this function) drops
        real Hamiltonian matrix elements between same-symmetry partial waves. `length(kappas) == 1` (e.g. bare
        hydrogen, where kappa alone fixes (J,parity)) is simply the special case where that pooling has nothing
        to do, not a separate code path.
"""
function computeSpectrum(initialLevel::Level, finalIonLevel::Level, kappas::Array{Int64,1}, symmetry::LevelSymmetry,
                         nm::Nuclear.Model, grid::Radial.Grid, settings::PhotoIonizationPseudostates.Settings)
    # ladderGrid is grid itself unless settings.pseudostateRbox asks for a larger box -- see ladderGridFor's own
    # docstring for why extending rbox this way needs no basis reconciliation. EVERYTHING below uses ladderGrid,
    # not grid, so the pseudostate ladder's full extent is correctly integrated over wherever it actually reaches.
    ladderGrid = PhotoIonizationPseudostates.ladderGridFor(grid, settings)
    PhotoIonizationPseudostates.checkFrozenOrbitalConsistency(initialLevel, finalIonLevel, ladderGrid)
    # ONE ladder PER kappa, kept separate here (not yet merged) so the CSF index range each kappa contributes to
    # combinedCsfs, below, can be recorded for weightByKappa -- Subshell(1000+k, kappa) never collides ACROSS
    # kappas, since kappa itself is part of a Subshell's identity, so merging the ladders themselves is safe.
    laddersByKappa = Dict{Int64,Dict{Subshell,Orbital}}( kappa => PhotoIonizationPseudostates.buildPseudostateLadder(
                                                          kappa, finalIonLevel, nm, ladderGrid, grid, settings) for kappa in kappas )
    ladderAll = merge( values(laddersByKappa)... )

    # The final ion's own DOMINANT CSF supplies the core occupation every pseudostate configuration shares.
    dominantIdx = argmax( abs.(finalIonLevel.mc) )
    dominantW2  = finalIonLevel.mc[dominantIdx]^2
    dominantW2 < 0.9   &&
        @warn("PhotoIonizationPseudostates.computeSpectrum(): finalIonLevel's dominant CSF carries only " *
              "$dominantW2 of the weight -- the single-CSF ion-core picture the pseudostate ladder assumes " *
              "is a poorer approximation here than usual.")
    coreCsf          = finalIonLevel.basis.csfs[dominantIdx]
    coreOccupation   = Dict{Subshell,Int64}( finalIonLevel.basis.subshells[i] => coreCsf.occupation[i]
                                              for i = 1:length(finalIonLevel.basis.subshells)
                                              if  coreCsf.occupation[i] > 0 )

    # MUST be in Base.isless(::Subshell,...)'s own canonical (n, then |kappa|, then kappa) order -- not a bare
    # vcat/Dict-key order -- since PhotoEmission.amplitude's own internal Basics.merge(...) produces exactly
    # that sorted order and then requires THIS list to appear, in order, as a subsequence of it.
    extendedSubshells = Base.sort( vcat(initialLevel.basis.subshells, collect(keys(ladderAll))), lt=Base.isless )
    paddedBasis       = ManyElectron.Basis(initialLevel.basis, extendedSubshells)
    boundCsfs         = CsfR[];   for  csf in paddedBasis.csfs
                                      LevelSymmetry(csf.J, csf.parity) == symmetry   &&   push!(boundCsfs, csf)   end

    # APPEND each kappa's own pseudostate CSFs in turn, recording the index RANGE it occupies in combinedCsfs --
    # this range is what lets weightByKappa, below, read off a per-eigenstate partial-channel weight directly
    # from the already-computed mixing vector, with no further matrix element machinery.
    combinedCsfs = copy(boundCsfs);   nBound = length(boundCsfs)
    kappaRanges  = Dict{Int64,UnitRange{Int64}}()
    for  kappa  in  kappas
        pseudoRelConfs = ConfigurationR[]
        for  sh  in  keys(laddersByKappa[kappa])
            subshells = merge( coreOccupation, Dict(sh => 1) )
            push!( pseudoRelConfs, ConfigurationR(subshells, finalIonLevel.basis.NoElectrons + 1) )
        end
        pseudoCsfs = CsfR[]
        for  relconf in pseudoRelConfs
            newCsfs = Basics.generateCsfRs(relconf, extendedSubshells)
            for  csf in newCsfs    LevelSymmetry(csf.J, csf.parity) == symmetry   &&   push!(pseudoCsfs, csf)    end
        end
        length(pseudoCsfs) == 0   &&   error("PhotoIonizationPseudostates.computeSpectrum(): no pseudostate CSF " *
                                              "of kappa = $kappa reaches symmetry $symmetry from this ion core, " *
                                              "although PhotoIonization.determineChannels said it should.")
        startIdx = length(combinedCsfs) + 1
        append!( combinedCsfs, pseudoCsfs )
        kappaRanges[kappa] = startIdx:length(combinedCsfs)
    end

    combinedOrbitals  = merge(initialLevel.basis.orbitals, ladderAll)
    coreSubshellList  = Subshell[]
    for  (k, sh)  in  enumerate(extendedSubshells)
        mocc = Basics.subshell_2j(sh) + 1;    isFilled = true
        for  csf in combinedCsfs    if  csf.occupation[k] != mocc    isFilled = false;   break   end   end
        isFilled   &&   push!(coreSubshellList, sh)
    end
    combinedBasis = Basis(true, initialLevel.basis.NoElectrons, extendedSubshells, combinedCsfs, coreSubshellList,
                          combinedOrbitals)

    mp = Hamiltonian.performCI(combinedBasis, nm, ladderGrid, AsfSettings(); printout=settings.printBefore, writeSummary=false)

    sortedLevels = sort(mp.levels, by = lev -> lev.energy)
    gi = Basics.twice(initialLevel.J) + 1.0
    states = StickState[]
    for  (n, lev)  in  enumerate(sortedLevels)
        weightOnReference = sum( lev.mc[1:nBound] .^ 2 )
        weightByKappa      = Dict{Int64,Float64}( kappa => sum( lev.mc[kappaRanges[kappa]] .^ 2 ) for kappa in kappas )
        omega = lev.energy - initialLevel.energy
        dipoleStrength = EmProperty(0., 0.)
        for  mp_ in settings.multipoles
            if  string(mp_)[1] == 'E'
                aC = PhotoEmission.amplitude(Absorption(), mp_, Basics.Coulomb,   omega, lev, initialLevel, ladderGrid; printout=false)
                aB = PhotoEmission.amplitude(Absorption(), mp_, Basics.Babushkin, omega, lev, initialLevel, ladderGrid; printout=false)
                dipoleStrength = dipoleStrength + EmProperty(abs2(aC), abs2(aB))
            else
                aM = PhotoEmission.amplitude(Absorption(), mp_, Basics.Magnetic, omega, lev, initialLevel, ladderGrid; printout=false)
                dipoleStrength = dipoleStrength + EmProperty(abs2(aM), abs2(aM))
            end
        end
        # THE SAME PHYSICAL PREFACTOR PhotoIonization.computeCrossSection APPLIES TO AN ENERGY-NORMALIZED
        # continuum amplitude -- 8*pi^3/(alpha*omega*gi) -- converts a squared amplitude into a cross section;
        # it is MISSING here on purpose until this point, since nothing before it needs it, and omega varies
        # per stick state so it cannot be folded into a single constant. Found missing (05-Nov-2026) when
        # branch (a)'s own PI-limit check came out ~10^4 too small -- Phase 0's toy validated only the SHAPE
        # (a ratio, in which any missing constant factor cancels), never the absolute scale.
        if  omega > 0.
            csFactor = 8 * pi^3 / Defaults.getDefaults("alpha") / omega / gi
            dipoleStrength = EmProperty(csFactor * dipoleStrength.Coulomb, csFactor * dipoleStrength.Babushkin)
        end
        if      n == 1               localSpacing = sortedLevels[2].energy - sortedLevels[1].energy
        elseif  n == length(sortedLevels)   localSpacing = sortedLevels[end].energy - sortedLevels[end-1].energy
        else                           localSpacing = (sortedLevels[n+1].energy - sortedLevels[n-1].energy) / 2
        end
        push!( states, StickState(lev.energy, dipoleStrength, localSpacing, weightOnReference, weightByKappa) )
    end

    return( Spectrum(symmetry, kappas, initialLevel, finalIonLevel, states) )
end


"""
`PhotoIonizationPseudostates.computeSpectra(finalMultiplet::Multiplet, initialLevel::Level, nm::Nuclear.Model,
        grid::Radial.Grid, settings::PhotoIonizationPseudostates.Settings; output::Bool=true)`
    ... the standard entry point of this module: for every final-ion level of finalMultiplet, groups every
        partial wave `PhotoIonization.determineChannels` finds reachable from initialLevel (under
        settings.multipoles/lValues) by the TOTAL SYMMETRY it reaches -- not by kappa -- and builds and
        diagonalizes ONE combined CI per (finalIonLevel, symmetry) pooling every kappa that shares it
        (`computeSpectrum`; see its own and `Spectrum`'s docstring for why pooling by symmetry, not kappa, is
        the point). A list of spectra::Array{PhotoIonizationPseudostates.Spectrum,1} is returned if output =
        true, and nothing otherwise.

        TAKES ONE initialLevel, NOT A WHOLE initialMultiplet -- deliberately, unlike `PhotoIonizationInterference
        .computePathways`. Looping over every level of a multiplet that may ALSO contain the autoionizing
        levels (see this module's own docstring) would have no way to tell a genuine below-threshold initial
        state apart from one of the higher-lying autoionizing ones; the caller already knows which level of
        their own multiplet is the one to ionize FROM, and passes it directly.
"""
function computeSpectra(finalMultiplet::Multiplet, initialLevel::Level, nm::Nuclear.Model, grid::Radial.Grid,
                        settings::PhotoIonizationPseudostates.Settings; output::Bool=true)
    println("")
    printstyled("PhotoIonizationPseudostates.computeSpectra(): the diagonalize-once pseudostate-continuum " *
                "computation starts now ... \n", color=:light_green)
    printstyled("---------------------------------------------------------------------------------------- \n",
                color=:light_green)
    println("")
    piSettings = PhotoIonization.Settings(PhotoIonization.Settings(); multipoles=settings.multipoles,
                                          gauges=settings.gauges, lValues=settings.lValues)
    spectra = PhotoIonizationPseudostates.Spectrum[]
    for  finalIonLevel  in  finalMultiplet.levels
        partialWaves = PhotoIonization.determineChannels(finalIonLevel, initialLevel, piSettings)
        # INVERT the (kappa -> symmetries) map determineChannels returns into (symmetry -> kappas): every
        # kappa reaching the SAME symmetry must be pooled into ONE computeSpectrum call, not one call each --
        # see Spectrum's own docstring for why.
        kappasBySymmetry = Dict{LevelSymmetry,Array{Int64,1}}()
        for  pw  in  partialWaves,  ch  in  pw.channels
            kappasBySymmetry[ch.symmetry] = push!( get(kappasBySymmetry, ch.symmetry, Int64[]), pw.kappa )
        end
        for  (symmetry, kappas)  in  kappasBySymmetry
            push!( spectra, PhotoIonizationPseudostates.computeSpectrum(initialLevel, finalIonLevel,
                                                                        unique(kappas), symmetry, nm, grid, settings) )
        end
    end
    PhotoIonizationPseudostates.displayResults(stdout, spectra, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    PhotoIonizationPseudostates.displayResults(iostream, spectra, settings)    end

    if    output    return( spectra )
    else            return( nothing )
    end
end


"""
`PhotoIonizationPseudostates.displayResults(stream::IO, spectra::Array{PhotoIonizationPseudostates.Spectrum,1},
        settings::PhotoIonizationPseudostates.Settings)`
    ... to display, for every spectrum and every requested photon energy, the density-of-states normalized
        cross section from `crossSectionAt`. A neat table is printed but nothing is returned otherwise.
"""
function displayResults(stream::IO, spectra::Array{PhotoIonizationPseudostates.Spectrum,1},
                        settings::PhotoIonizationPseudostates.Settings)
    nx = 100
    println(stream, " ")
    println(stream, "  Diagonalize-once pseudostate-continuum photoionization:  cross sections [a.u.], Coulomb gauge:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  "
    sa = sa * TableStrings.center(10, "i-level"; na=2)
    sa = sa * TableStrings.center(10, "f-level"; na=2)
    sa = sa * TableStrings.center(10, "symmetry"; na=2)
    sa = sa * TableStrings.center(10, "kappas"; na=2)
    sa = sa * TableStrings.center(16, "E_photon"; na=2)
    sa = sa * TableStrings.center(16, "sigma (Coulomb)"; na=2)
    println(stream, sa);    println(stream, "  ", TableStrings.hLine(nx))

    for  spectrum in spectra,  omega in settings.photonEnergies
        omegaAu = Defaults.convertUnits("energy: to atomic", omega)
        sigma   = PhotoIonizationPseudostates.crossSectionAt(spectrum, omegaAu)
        sa = "  " * TableStrings.center(10, string(spectrum.initialLevel.index); na=2)
        sa = sa * TableStrings.center(10, string(spectrum.finalIonLevel.index); na=2)
        sa = sa * TableStrings.center(10, string(spectrum.symmetry); na=2)
        sa = sa * TableStrings.center(10, string(spectrum.kappas); na=2)
        sa = sa * @sprintf("%14.6e", omega) * "  "
        sa = sa * @sprintf("%14.6e", sigma.Coulomb)
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx), "\n")

    return( nothing )
end

end # module
