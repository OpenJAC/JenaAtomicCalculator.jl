
"""
`module  JAC.WeakInteractionEnhancement`
... a submodel of JAC that computes the OBSERVABLES which the P-odd and P,T-odd interactions of `WeakInteractionMoment` give rise to, and
    which -- unlike those bare amplitudes -- all require a SUM OVER INTERMEDIATE STATES.  Two are provided:

        E1_PNC  ... the parity-non-conserving electric-dipole amplitude between two levels of the SAME parity, which is what a PNC
                    experiment on a forbidden transition such as the caesium 6s-7s line actually measures.  It comes in two kinds, and
                    they are DIFFERENT OBSERVABLES rather than a leading term and a correction:
                    -- the NUCLEAR-SPIN-INDEPENDENT part, driven by the weak charge, which is a scalar and therefore lives between
                       fine-structure levels;  `computePncE1Amplitude`, re-coupled to hyperfine levels by `computePncE1AmplitudeNsi`;
                    -- the NUCLEAR-SPIN-DEPENDENT part, driven by the nuclear anapole moment, whose operator is a scalar only in the
                       COUPLED electron-plus-nucleus space and which therefore exists only between hyperfine levels and only for a
                       nucleus of non-zero spin;  `computePncE1AmplitudeNsd`.
                    The distinction is not academic.  The weak charge cannot change J, so a transition such as barium 6s_1/2 - 5d_5/2
                    has NO spin-independent amplitude at all and is a clean anapole probe; caesium 6s-7s is the opposite case, where the
                    weak charge dominates and the anapole is the few-per-cent difference between two hyperfine components.
        R       ... the enhancement factor d_atom / d_e of an electron electric-dipole moment in the atomic state.

    Both are second-order quantities of the same shape, differing only in which P-odd operator is admixed and in how the dipole is then
    projected:

        E1_PNC(f<-i) = SUM_n [ <f||D||n> <n|H_W|i> / (E_i - E_n)  +  <f|H_W|n> <n||D||i> / (E_f - E_n) ]
        d_atom       = 2 SUM_n <0|D_z|n> <n|H_EDM|0> / (E_0 - E_n)

    THE INTERMEDIATE STATES ARE SUPPLIED BY THE CALLER as an ordinary `Multiplet` in `Settings.gMultiplet`, obtained from a separate
    `Atomic.Computation` over whatever opposite-parity configurations are thought to matter.  No Green function and no pseudo-continuum is
    attempted.  That is a deliberate choice and its price differs sharply between the two observables, which is the single most important
    thing to know before reading a number out of this module:

    -- E1_PNC IS DOMINATED BY A FEW LOW-LYING np LEVELS, so a named, inspectable set of intermediates is close to how the calculation is
       actually done.  It is still convergence-limited, and the example accordingly runs more than one `gMultiplet` and shows the spread
       rather than quoting one number as though it were converged.

    -- R IS A CANCELLING SUM AND MUST BE READ DIFFERENTLY.  Schiff's theorem states that the electric-dipole moment of a POINT,
       NON-RELATIVISTIC bound system vanishes identically; the whole enhancement is what survives that cancellation once relativity is
       restored, so a truncated sum retains part of what was meant to cancel.  What that costs in practice was MEASURED rather than
       assumed, in branch d of `example-Cnnew.jl`, and the answer is more nuanced than the theorem alone suggests: R does converge for a
       Li-like ion, but more slowly than E1_PNC over the very same sequence of intermediate sets -- 3.2 % against 1.0 % at the last step --
       which is the signature of a cancellation eating into the accuracy without destroying it.  In a hydrogen-like ion, by contrast, the
       two members of a near-degenerate pair return values equal and opposite to six digits: the cancellation caught in the act with two
       levels, and a warning that R there is dominated by the same wrong denominator as E1_PNC (see the next paragraph).

       R is therefore printed as UNCONVERGED rather than as a result.  That is the honest label: the observed convergence is encouraging,
       but nothing anchors the ABSOLUTE scale of the EDM operator the way the closed-form hydrogenic matrix element now anchors the weak
       charge, so a converged R here is a converged value of an unverified quantity.

    A HYDROGEN-LIKE ION IS THE WRONG TEST SYSTEM FOR BOTH OBSERVABLES, which is worth knowing before a plausible-looking number is taken
    from one.  There the 2s_1/2 and 2p_1/2 levels are exactly degenerate in Dirac theory and are split only by the Lamb shift, of which a
    JAC calculation contains the finite-nuclear-size part alone; at Z = 20 the missing QED is some 380 times larger than the splitting that
    remains.  That splitting sits in a DENOMINATOR of both sums, so the error passes straight into the answer and no enlargement of the
    gMultiplet repairs it.  A many-electron system, in which screening does the splitting, is the regime in which this module can be
    believed.

    THE EDM OPERATOR LIVES HERE RATHER THAN IN `WeakInteractionMoment`, where it belongs by kind.  `H_EDM = -d_e beta Sigma.E` is a bare
    one-electron amplitude and would sit naturally beside the weak charge, the anapole and the Schiff moment; it was placed here only to
    keep a single source module per task.  It should be MOVED to `WeakInteractionMoment` when that module is next opened, together with
    `radialIntegralPPminus`, which is of the same family as the three integrals already there.
"""
module WeakInteractionEnhancement


using  Printf, ..AngularMomentum, ..Basics, ..Defaults, ..InteractionStrength, ..ManyElectron, ..Nuclear, ..Radial, ..TableStrings,
       ..WeakInteractionMoment


"""
`struct  WeakInteractionEnhancement.Settings  <:  AbstractPropertySettings`
    ... defines a type for the details and parameters of computing the parity-non-conserving E1 amplitude and the electron-EDM enhancement
        factor.

    + calcPncE1                ::Bool             ... True if the PNC E1 amplitudes between the selected levels are to be computed.
    + calcEdmEnhancement       ::Bool             ... True if the electron-EDM enhancement factor R = d_atom/d_e is to be computed; see the
                                                        module docstring for why R from a truncated sum is reported as unconverged.
    + gMultiplet               ::Multiplet        ... The intermediate levels of the sum over states, supplied by the caller; only those of
                                                        parity opposite to the level in question contribute, and the others are ignored.
    + printBefore              ::Bool             ... True if a list of selected levels is printed before the actual computations start.
    + levelSelection           ::LevelSelection   ... Specifies the selected levels, if any.
"""
struct Settings  <:  AbstractPropertySettings
    calcPncE1                  ::Bool
    calcEdmEnhancement         ::Bool
    gMultiplet                 ::Multiplet
    printBefore                ::Bool
    levelSelection             ::LevelSelection
end


"""
`WeakInteractionEnhancement.Settings()`
    ... constructor for an `empty` instance of WeakInteractionEnhancement.Settings; a settings::WeakInteractionEnhancement.Settings is
        returned.  The empty `gMultiplet` must be replaced before anything can be computed, since the sum over states has no intermediates
        without it.
"""
function Settings()
    Settings(false, false, Multiplet(), false, LevelSelection() )
end


"""
`WeakInteractionEnhancement.Settings(set::WeakInteractionEnhancement.Settings;`

        calcPncE1=.., calcEdmEnhancement=.., gMultiplet=.., printBefore=.., levelSelection=..)

    ... keyword copy-constructor for re-defining selected values of a settings::WeakInteractionEnhancement.Settings; a
        settings::WeakInteractionEnhancement.Settings is returned.
"""
function Settings(set::WeakInteractionEnhancement.Settings;
        calcPncE1::Union{Nothing,Bool}=nothing,            calcEdmEnhancement::Union{Nothing,Bool}=nothing,
        gMultiplet::Union{Nothing,Multiplet}=nothing,      printBefore::Union{Nothing,Bool}=nothing,
        levelSelection::Union{Nothing,LevelSelection}=nothing)
    if  isnothing(calcPncE1)           calcPncE1x          = set.calcPncE1          else   calcPncE1x          = calcPncE1          end
    if  isnothing(calcEdmEnhancement)  calcEdmEnhancementx = set.calcEdmEnhancement else   calcEdmEnhancementx = calcEdmEnhancement end
    if  isnothing(gMultiplet)          gMultipletx         = set.gMultiplet         else   gMultipletx         = gMultiplet         end
    if  isnothing(printBefore)         printBeforex        = set.printBefore        else   printBeforex        = printBefore        end
    if  isnothing(levelSelection)      levelSelectionx     = set.levelSelection     else   levelSelectionx     = levelSelection     end

    Settings( calcPncE1x, calcEdmEnhancementx, gMultipletx, printBeforex, levelSelectionx )
end


# `Base.show(io::IO, settings::WeakInteractionEnhancement.Settings)`  ... prepares a proper printout of settings.
function Base.show(io::IO, settings::WeakInteractionEnhancement.Settings)
    println(io, "calcPncE1:                $(settings.calcPncE1)  ")
    println(io, "calcEdmEnhancement:       $(settings.calcEdmEnhancement)  ")
    println(io, "gMultiplet:               $(length(settings.gMultiplet.levels)) intermediate levels  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "levelSelection:           $(settings.levelSelection)  ")
end


"""
`struct  WeakInteractionEnhancement.Outcome`
    ... defines a type to keep the outcome of a weak-interaction-enhancement computation for ONE level.

    + level                     ::Level                            ... Atomic level to which the outcome refers to.
    + edmEnhancement            ::Float64                          ... R = d_atom/d_e for this level, from the truncated sum over
                                                                        gMultiplet; UNCONVERGED rather than an estimate, and of an
                                                                        operator whose absolute scale is unverified -- see the module
                                                                        docstring before quoting it.
    + pncAmplitudes             ::Array{Tuple{Level,ComplexF64},1} ... List of (partner level, E1_PNC) pairs, in which the amplitude
                                                                        <level || D_PNC || partner> is given in atomic units (e a_0) and is
                                                                        purely imaginary for real orbitals.
"""
struct Outcome
    level                       ::Level
    edmEnhancement              ::Float64
    pncAmplitudes               ::Array{Tuple{Level,ComplexF64},1}
end


"""
`WeakInteractionEnhancement.Outcome()`
    ... constructor for an `empty` instance of WeakInteractionEnhancement.Outcome; an outcome::WeakInteractionEnhancement.Outcome is
        returned.
"""
function Outcome()
    Outcome(Level(), 0., Tuple{Level,ComplexF64}[] )
end


# `Base.show(io::IO, outcome::WeakInteractionEnhancement.Outcome)`  ... prepares a proper printout of outcome.
function Base.show(io::IO, outcome::WeakInteractionEnhancement.Outcome)
    println(io, "level:                   $(outcome.level)  ")
    println(io, "edmEnhancement:          $(outcome.edmEnhancement)  ")
    println(io, "pncAmplitudes:           $(outcome.pncAmplitudes)  ")
end


"""
`WeakInteractionEnhancement.computeEdmEnhancement(level::Level, gMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid)`
    ... to compute the enhancement factor R = d_atom/d_e of an electron electric-dipole moment in the given level, from the sum over the
        opposite-parity intermediate levels of gMultiplet.  A value::Float64 is returned.

        The atomic moment is the expectation of D_z in the perturbed state of the stretched magnetic substate M = J,

            d_atom = 2 SUM_n <0|D_z|n> <n|H_EDM|0> / (E_0 - E_n),

        in which H_EDM is a scalar, so that only intermediates of the same J and opposite parity contribute.  Using the projection
        <j j|T^1_0|j j> = sqrt( j / ((j+1)(2j+1)) ) <j||T^1||j> the whole sum reduces to reduced matrix elements.  R is returned for a unit
        electron EDM and is therefore dimensionless.  Read the module docstring before using the number: Schiff's theorem makes this a
        CANCELLING sum, which converges more slowly than E1_PNC does over the same intermediate sets, and in a hydrogen-like ion it is
        dominated by the 2s-2p_1/2 denominator that a calculation without QED gets wrong by some two orders of magnitude.
"""
function computeEdmEnhancement(level::Level, gMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid)
    jj = Basics.twice(level.J) / 2.0
    if  jj == 0.    return( 0. )    end
    wa = ComplexF64(0.)
    for  nLevel  in  gMultiplet.levels
        if  nLevel.parity == level.parity   ||   nLevel.J != level.J    continue    end
        dE = level.energy - nLevel.energy
        if  abs(dE) < 1.0e-12                                           continue    end
        hOrd = WeakInteractionEnhancement.edmAmplitude(nLevel, level, nm, grid) / sqrt( Basics.twice(level.J) + 1.0 )
        wa   = wa + WeakInteractionEnhancement.dipoleReducedMe(level, nLevel, grid) * hOrd / dE
    end

    return( 2.0 * sqrt( jj / ((jj + 1.0) * (2.0*jj + 1.0)) ) * real(wa) )
end


"""
`WeakInteractionEnhancement.computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,`
                                            `settings::WeakInteractionEnhancement.Settings; output=true)`
    ... to compute (as selected) the parity-non-conserving E1 amplitudes and the electron-EDM enhancement factors for the levels of the
        given multiplet.  The intermediate states of both sums are taken from settings.gMultiplet.  The results are printed in neat tables
        to screen and, if output, an Array{WeakInteractionEnhancement.Outcome,1} is returned; nothing otherwise.
"""
function computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                          settings::WeakInteractionEnhancement.Settings; output=true)
    println("")
    printstyled("WeakInteractionEnhancement.computeOutcomes(): The computation of PNC amplitudes and EDM enhancements starts now ... \n",
                color=:light_green)
    printstyled("--------------------------------------------------------------------------------------------------------------------- \n",
                color=:light_green)
    #
    if  length(settings.gMultiplet.levels) == 0
        error("WeakInteractionEnhancement.computeOutcomes(): settings.gMultiplet is empty. Both observables of this module are sums over " *
              "intermediate states, so the caller must supply the opposite-parity levels to be summed over; there is no default.")
    end
    #
    outcomes = WeakInteractionEnhancement.determineOutcomes(multiplet, settings)
    if  settings.printBefore    WeakInteractionEnhancement.displayOutcomes(stdout, outcomes, settings)    end
    #
    newOutcomes = WeakInteractionEnhancement.Outcome[]
    for  outcome in outcomes
        edm = 0.;   pncs = Tuple{Level,ComplexF64}[]
        if  settings.calcEdmEnhancement
            edm = WeakInteractionEnhancement.computeEdmEnhancement(outcome.level, settings.gMultiplet, nm, grid)
        end
        if  settings.calcPncE1
            # each same-parity pair is reported once, from the higher level of the two, so that a pair does not appear twice
            for  iLevel  in  multiplet.levels
                if  iLevel.index >= outcome.level.index   ||   iLevel.parity != outcome.level.parity              continue    end
                if  !AngularMomentum.isTriangle(outcome.level.J, AngularJ64(1), iLevel.J)                         continue    end
                amp = WeakInteractionEnhancement.computePncE1Amplitude(outcome.level, iLevel, settings.gMultiplet, nm, grid)
                push!( pncs, (iLevel, amp) )
            end
        end
        push!( newOutcomes, WeakInteractionEnhancement.Outcome(outcome.level, edm, pncs) )
    end
    #
    WeakInteractionEnhancement.displayResults(stdout, newOutcomes, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    WeakInteractionEnhancement.displayResults(iostream, newOutcomes, settings)   end
    #
    if    output    return( newOutcomes )
    else            return( nothing )
    end
end


"""
`WeakInteractionEnhancement.computePncE1Amplitude(finalLevel::Level, initialLevel::Level, gMultiplet::Multiplet, nm::Nuclear.Model,`
                                                  `grid::Radial.Grid)`
    ... to compute the parity-non-conserving electric-dipole amplitude <finalLevel || D_PNC || initialLevel> between two levels of the SAME
        parity, from the sum over the opposite-parity intermediate levels of gMultiplet.  A value::ComplexF64 is returned, in atomic units
        (e a_0).

        The weak-charge Hamiltonian is a scalar and so admixes opposite-parity states without altering the angular structure; the amplitude
        therefore remains an ordinary rank-1 reduced matrix element,

            <f||D_PNC||i> = SUM_n [ <f||D||n> <n|H_W|i>/(E_i - E_n)  +  <f|H_W|n> <n||D||i>/(E_f - E_n) ],

        in which <n|H_W|i> is the ORDINARY matrix element, i.e. the reduced one divided by sqrt(2J+1), and is nonzero only for J_n = J_i.
        The result is purely imaginary for real radial orbitals, since the weak-charge amplitude is.
"""
function computePncE1Amplitude(finalLevel::Level, initialLevel::Level, gMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid)
    # a PNC amplitude exists only between levels of EQUAL parity; the dipole itself supplies the usual triangle rule
    if      finalLevel.parity != initialLevel.parity                                    return( ComplexF64(0.) )    end
    if      !AngularMomentum.isTriangle(finalLevel.J, AngularJ64(1), initialLevel.J)     return( ComplexF64(0.) )    end
    wa = ComplexF64(0.)
    for  nLevel  in  gMultiplet.levels
        if  nLevel.parity == initialLevel.parity    continue    end
        #
        # first term: the weak interaction acts on the INITIAL level, the dipole on the admixture
        if  nLevel.J == initialLevel.J
            dE = initialLevel.energy - nLevel.energy
            if  abs(dE) > 1.0e-12
                wOrd = WeakInteractionMoment.weakChargeAmplitude(nLevel, initialLevel, nm, grid) /
                       sqrt( Basics.twice(initialLevel.J) + 1.0 )
                wa   = wa + WeakInteractionEnhancement.dipoleReducedMe(finalLevel, nLevel, grid) * wOrd / dE
            end
        end
        #
        # second term: the weak interaction acts on the FINAL level
        if  nLevel.J == finalLevel.J
            dE = finalLevel.energy - nLevel.energy
            if  abs(dE) > 1.0e-12
                wOrd = WeakInteractionMoment.weakChargeAmplitude(finalLevel, nLevel, nm, grid) / sqrt( Basics.twice(nLevel.J) + 1.0 )
                wa   = wa + wOrd * WeakInteractionEnhancement.dipoleReducedMe(nLevel, initialLevel, grid) / dE
            end
        end
    end

    return( wa )
end


"""
`WeakInteractionEnhancement.computePncE1AmplitudeNsd(finalLevel::Level, Ff::AngularJ64, initialLevel::Level, Fi::AngularJ64,`
                                                     `gMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid;`
                                                     `kappaAnapole::Float64=1.0)`
    ... to compute the NUCLEAR-SPIN-DEPENDENT parity-non-conserving electric-dipole amplitude <f J_f I F_f || D_PNC || i J_i I F_i> between
        two HYPERFINE levels of the same electronic parity, from the sum over the opposite-parity intermediate levels of gMultiplet.  A
        value::ComplexF64 is returned, in atomic units (e a_0).

        THIS IS A DIFFERENT OBSERVABLE FROM `computePncE1Amplitude` AND NOT A REFINEMENT OF IT.  The weak charge is a scalar, so the
        amplitude it drives is an ordinary rank-1 reduced matrix element between FINE-STRUCTURE levels and the hyperfine structure only
        re-couples it.  The anapole moment is not: its operator is the scalar product `alpha.I rho(r)` of an electronic rank-1 tensor with
        the nuclear spin, so it is a scalar only in the COUPLED electron-plus-nucleus space, and the amplitude it drives exists only between
        hyperfine levels.  A nucleus of spin zero has no anapole amplitude at all, and this function returns zero for one.

        The consequence that matters for an experiment is a selection rule.  Because the weak charge cannot change J, the intermediate level
        must carry J_n = J_i on one side and J_n = J_f on the other, and a transition such as barium 6s_1/2 - 5d_5/2 then has NO
        nuclear-spin-independent amplitude whatever: no single E1 step reaches from J = 1/2 to J = 5/2.  The anapole operator obeys only the
        triangle |J_n - J_i| <= 1, so p_3/2 intermediates open and the amplitude is non-zero.  Such a transition is therefore a CLEAN
        anapole probe rather than a weak-charge measurement with a small spin-dependent correction, which is the reverse of the caesium
        6s-7s situation.

        The two geometrical factors are `hyperfineDipoleFactor` and `hyperfineScalarFactor`; the electronic reduced matrix element of the
        anapole operator is `WeakInteractionMoment.anapoleAmplitude`, which already carries the `G_F/sqrt(2) kappa` prefactor.  Since the
        anapole Hamiltonian is diagonal in F, the admixture on the initial side carries F_i and that on the final side F_f.

        + kappaAnapole  ::Float64   ... the dimensionless nuclear anapole constant, left at 1.0 so that the returned amplitude is the
                                        coefficient OF kappa and may be scaled by whatever value the user adopts.  The convention is the one
                                        of Flambaum and Khriplovich, `H = (G_F/sqrt2) kappa (alpha.I/I) rho(r)`, i.e. with the 1/I that
                                        makes kappa comparable with the caesium value 0.112(16) of Wood et al., Science 275, 1759 (1997).
"""
function computePncE1AmplitudeNsd(finalLevel::Level, Ff::AngularJ64, initialLevel::Level, Fi::AngularJ64,
                                  gMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid; kappaAnapole::Float64=1.0)
    # a PNC amplitude exists only between levels of EQUAL electronic parity, and the dipole supplies the triangle rule in the F's
    if      finalLevel.parity != initialLevel.parity                        return( ComplexF64(0.) )    end
    if      !AngularMomentum.isTriangle(Ff, AngularJ64(1), Fi)              return( ComplexF64(0.) )    end
    if      Basics.twice(nm.spinI) == 0                                     return( ComplexF64(0.) )    end
    wa = ComplexF64(0.)
    for  nLevel  in  gMultiplet.levels
        if  nLevel.parity == initialLevel.parity    continue    end
        #
        # first term: the anapole acts on the INITIAL level, which fixes the admixture to F_i, and the dipole on the admixture
        dE = initialLevel.energy - nLevel.energy
        if  abs(dE) > 1.0e-12
            wg = WeakInteractionEnhancement.hyperfineScalarFactor(nLevel.J, initialLevel.J, Fi, nm.spinI) *
                 WeakInteractionEnhancement.hyperfineDipoleFactor(finalLevel.J, Ff, nLevel.J, Fi, nm.spinI)
            if  wg != 0.
                wa = wa + wg * WeakInteractionEnhancement.dipoleReducedMe(finalLevel, nLevel, grid) *
                          WeakInteractionMoment.anapoleAmplitude(nLevel, initialLevel, nm, grid; kappaAnapole=kappaAnapole) / dE
            end
        end
        #
        # second term: the anapole acts on the FINAL level, which fixes the admixture to F_f
        dE = finalLevel.energy - nLevel.energy
        if  abs(dE) > 1.0e-12
            wg = WeakInteractionEnhancement.hyperfineScalarFactor(finalLevel.J, nLevel.J, Ff, nm.spinI) *
                 WeakInteractionEnhancement.hyperfineDipoleFactor(nLevel.J, Ff, initialLevel.J, Fi, nm.spinI)
            if  wg != 0.
                wa = wa + wg * WeakInteractionMoment.anapoleAmplitude(finalLevel, nLevel, nm, grid; kappaAnapole=kappaAnapole) *
                          WeakInteractionEnhancement.dipoleReducedMe(nLevel, initialLevel, grid) / dE
            end
        end
    end

    return( wa )
end


"""
`WeakInteractionEnhancement.computePncE1AmplitudeNsi(finalLevel::Level, Ff::AngularJ64, initialLevel::Level, Fi::AngularJ64,`
                                                     `gMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid)`
    ... to re-couple the nuclear-spin-independent amplitude of `computePncE1Amplitude` into the hyperfine basis, giving
        <f J_f I F_f || D_PNC || i J_i I F_i>.  A value::ComplexF64 is returned, in atomic units (e a_0).

        It carries no physics of its own: the weak charge is a purely electronic scalar, so the F-dependence is the single geometrical
        factor `hyperfineDipoleFactor` and nothing else.  The function exists so that the two contributions to one measured hyperfine line
        can be added, which is what an experiment sees, and so that their RATIO can be formed on the same footing -- the caesium anapole
        moment was extracted from exactly that ratio, as the small difference between the F = 3 -> 4 and F = 4 -> 3 components of one line.
"""
function computePncE1AmplitudeNsi(finalLevel::Level, Ff::AngularJ64, initialLevel::Level, Fi::AngularJ64,
                                  gMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid)
    if      !AngularMomentum.isTriangle(Ff, AngularJ64(1), Fi)              return( ComplexF64(0.) )    end
    wg = WeakInteractionEnhancement.hyperfineDipoleFactor(finalLevel.J, Ff, initialLevel.J, Fi, nm.spinI)
    if      wg == 0.                                                        return( ComplexF64(0.) )    end

    return( wg * WeakInteractionEnhancement.computePncE1Amplitude(finalLevel, initialLevel, gMultiplet, nm, grid) )
end


"""
`WeakInteractionEnhancement.determineOutcomes(multiplet::Multiplet, settings::WeakInteractionEnhancement.Settings)`
    ... to determine a list of Outcome's for the computation of the enhancement properties for the given multiplet, taking the level
        selection of the settings into account.  An Array{WeakInteractionEnhancement.Outcome,1} is returned, in which all physical
        properties are still set to zero.
"""
function determineOutcomes(multiplet::Multiplet, settings::WeakInteractionEnhancement.Settings)
    outcomes = WeakInteractionEnhancement.Outcome[]
    for  level  in  multiplet.levels
        if  Basics.selectLevel(level, settings.levelSelection)
            push!( outcomes, WeakInteractionEnhancement.Outcome(level, 0., Tuple{Level,ComplexF64}[]) )
        end
    end

    return( outcomes )
end


"""
`WeakInteractionEnhancement.dipoleReducedMe(finalLevel::Level, initialLevel::Level, grid::Radial.Grid)`
    ... to compute the reduced electric-dipole matrix element <alpha_f J_f || D || alpha_i J_i> in the LENGTH form and in atomic units
        (e a_0).  A value::ComplexF64 is returned.

        It is deliberately built from `WeakInteractionMoment.oneParticleAmplitude`, the same contraction that produces the weak-charge
        amplitude, rather than from `MultipoleMoment.dipoleAmplitude`.  The two observables of this module are PRODUCTS of a dipole and a
        weak matrix element, so what matters is that both factors carry the same normalization; routing them through one contraction makes
        that automatic and leaves only an overall constant that the example branch measures.
"""
function dipoleReducedMe(finalLevel::Level, initialLevel::Level, grid::Radial.Grid)
    if      finalLevel.parity == initialLevel.parity                                    return( ComplexF64(0.) )    end
    if      !AngularMomentum.isTriangle(finalLevel.J, AngularJ64(1), initialLevel.J)     return( ComplexF64(0.) )    end
    kernel = (orba, orbb) -> ComplexF64( InteractionStrength.eMultipole(1, orba, orbb, grid) )

    return( WeakInteractionMoment.oneParticleAmplitude(1, kernel, finalLevel, initialLevel) )
end


"""
`WeakInteractionEnhancement.displayOutcomes(stream::IO, outcomes::Array{WeakInteractionEnhancement.Outcome,1},`
                                            `settings::WeakInteractionEnhancement.Settings)`
    ... to display the levels that have been selected, together with the size of the intermediate multiplet, before the computation starts.
        Nothing is returned.
"""
function displayOutcomes(stream::IO, outcomes::Array{WeakInteractionEnhancement.Outcome,1}, settings::WeakInteractionEnhancement.Settings)
    nx = 43
    println(stream, " ")
    println(stream, "  Selected weak-interaction-enhancement levels:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "Energy"; na=4)
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))
    #
    for  outcome in outcomes
        sa  = "  ";    sym = LevelSymmetry( outcome.level.J, outcome.level.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(outcome.level.index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", outcome.level.energy)) * "    "
        println(stream, sa )
    end
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, "  The sum over intermediate states uses $(length(settings.gMultiplet.levels)) levels of the given gMultiplet; only")
    println(stream, "  those of opposite parity, and of the total angular momentum required by the operator, actually contribute.")
    println(stream, " ")

    return( nothing )
end


"""
`WeakInteractionEnhancement.displayResults(stream::IO, outcomes::Array{WeakInteractionEnhancement.Outcome,1},`
                                           `settings::WeakInteractionEnhancement.Settings)`
    ... to display the PNC E1 amplitudes and the EDM enhancement factors that have been computed, each in its own table.  Nothing is
        returned.
"""
function displayResults(stream::IO, outcomes::Array{WeakInteractionEnhancement.Outcome,1}, settings::WeakInteractionEnhancement.Settings)
    if  settings.calcPncE1
        nx = 84
        println(stream, " ")
        println(stream, "  Parity-non-conserving E1 amplitudes  <f || D_PNC || i>,  in atomic units (e a_0):")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  "
        sa = sa * TableStrings.center(10, "f-level"; na=2)
        sa = sa * TableStrings.center(10, "i-level"; na=2)
        sa = sa * TableStrings.center(18, "J^P (f) -- J^P (i)"; na=4)
        sa = sa * TableStrings.center(30, "Re  and  Im  of  E1_PNC"; na=2)
        println(stream, sa);    println(stream, "  ", TableStrings.hLine(nx))
        #
        for  outcome in outcomes
            for  (iLevel, amp)  in  outcome.pncAmplitudes
                symf = LevelSymmetry( outcome.level.J, outcome.level.parity);    symi = LevelSymmetry( iLevel.J, iLevel.parity)
                sa   = "  "
                sa   = sa * TableStrings.center(10, TableStrings.level(outcome.level.index); na=2)
                sa   = sa * TableStrings.center(10, TableStrings.level(iLevel.index); na=2)
                sa   = sa * TableStrings.center(18, string(symf) * " -- " * string(symi); na=4)
                sa   = sa * @sprintf("%.6e", amp.re) * "   " * @sprintf("%.6e", amp.im)
                println(stream, sa )
            end
        end
        println(stream, "  ", TableStrings.hLine(nx))
        println(stream, "  The amplitude is purely imaginary for real radial orbitals; a nonzero real part is a numerical residue.")
        println(stream, "  It is convergence-limited by the choice of gMultiplet -- run more than one and compare before quoting a value.")
        println(stream, " ")
    end
    #
    if  settings.calcEdmEnhancement
        nx = 60
        println(stream, " ")
        println(stream, "  Electron-EDM enhancement factors  R = d_atom / d_e,  dimensionless:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  "
        sa = sa * TableStrings.center(10, "Level";  na=2)
        sa = sa * TableStrings.center(10, "J^P";    na=4)
        sa = sa * TableStrings.center(20, "R (unconverged)"; na=2)
        println(stream, sa);    println(stream, "  ", TableStrings.hLine(nx))
        #
        for  outcome in outcomes
            sym = LevelSymmetry( outcome.level.J, outcome.level.parity)
            sa  = "  "
            sa  = sa * TableStrings.center(10, TableStrings.level(outcome.level.index); na=2)
            sa  = sa * TableStrings.center(10, string(sym); na=4)
            sa  = sa * @sprintf("%.6e", outcome.edmEnhancement)
            println(stream, sa )
        end
        println(stream, "  ", TableStrings.hLine(nx))
        println(stream, "  UNCONVERGED, AND OF AN UNANCHORED OPERATOR.  Schiff's theorem makes d_atom a CANCELLING sum, which converges")
        println(stream, "  more slowly with the intermediate set than E1_PNC does; and nothing here fixes the absolute scale of the EDM")
        println(stream, "  operator the way a closed-form matrix element fixes the weak charge.  For a hydrogen-like ion the value is in")
        println(stream, "  addition dominated by the 2s-2p_1/2 denominator, which is set by QED and is absent from this calculation.")
        println(stream, " ")
    end

    return( nothing )
end


"""
`WeakInteractionEnhancement.edmAmplitude(finalLevel::Level, initialLevel::Level, nm::Nuclear.Model, grid::Radial.Grid;`
                                         `display::Bool=false)`
    ... to compute the electron-EDM amplitude <alpha_f J_f || H^(EDM) || alpha_i J_i> for a UNIT electron electric-dipole moment, in the
        nuclear field of the given model.  A value::ComplexF64 is returned.

        The operator is `H_EDM = -d_e beta Sigma.E` with `E = -grad phi` the internal electric field.  Sigma is an axial vector and E a
        polar one, so their contraction is a PSEUDOSCALAR: rank 0, P-odd and, since Sigma changes sign under time reversal while E does not,
        also T-odd.  It therefore obeys the same selection rules as the weak charge -- equal J, opposite parity, kappa_b = -kappa_a -- but a
        different radial structure, because beta Sigma.r does not mix the large and small components and instead weights them with opposite
        signs:

            <a|H_EDM|b>  =  d_e  INT E_r(r) [P_a P_b - Q_a Q_b] dr ,    kappa_b = -kappa_a

        using  <Om(kappa_a)|sigma.r|Om(kappa_b)> = -delta(kappa_a, -kappa_b)  for both the upper and the lower pair.  The amplitude is REAL
        for real radial orbitals, in contrast to the weak charge, whose gamma_5 supplies an explicit factor of i.
"""
function edmAmplitude(finalLevel::Level, initialLevel::Level, nm::Nuclear.Model, grid::Radial.Grid; display::Bool=false)
    if      finalLevel.parity == initialLevel.parity   amplitude = ComplexF64(0.)
    elseif  finalLevel.J      != initialLevel.J        amplitude = ComplexF64(0.)
    else
        # the electrostatic potential of the nucleus is phi(r) = Z(r)/r, and the radial field is E_r = -d phi/dr
        pot = Nuclear.nuclearPotential(nm, grid)
        phi = zeros( grid.NoPoints )
        mtp = min( grid.NoPoints, length(pot.Zr) )
        for  i = 2:mtp    phi[i] = pot.Zr[i] / grid.r[i]    end
        eField = -WeakInteractionMoment.radialDerivative(phi, grid)
        kernel = (orba, orbb) -> begin
            if  orba.subshell.kappa != -orbb.subshell.kappa    return( ComplexF64(0.) )    end
            ComplexF64( WeakInteractionEnhancement.radialIntegralPPminus(eField, orba, orbb, grid) )
        end
        amplitude = WeakInteractionMoment.oneParticleAmplitude(0, kernel, finalLevel, initialLevel)
    end
    #
    if  display
        sa = @sprintf("%.5e", amplitude.re) * "  " * @sprintf("%.5e", amplitude.im)
        println("    < level=$(finalLevel.index) [J=$(finalLevel.J)$(string(finalLevel.parity))] || H^(EDM) ||" *
                " $(initialLevel.index) [$(initialLevel.J)$(string(initialLevel.parity))] >  = " * sa)
    end

    return( amplitude )
end


"""
`WeakInteractionEnhancement.hyperfineDipoleFactor(Ja::AngularJ64, Fa::AngularJ64, Jb::AngularJ64, Fb::AngularJ64, spinI::AngularJ64)`
    ... to compute the geometrical factor that turns a purely ELECTRONIC rank-1 reduced matrix element <a J_a || T^(1) || b J_b> into the
        reduced matrix element <a J_a I F_a || T^(1) || b J_b I F_b> of the same operator in the coupled hyperfine basis.  A value::Float64
        is returned, and it is zero unless the F's satisfy the triangle rule.

            factor = (-1)^(J_a + I + F_b + 1) sqrt( (2F_a+1)(2F_b+1) ) { J_a  F_a  I ; F_b  J_b  1 }

        This is the standard reduction of an operator acting on one part of a coupled system (Edmonds, Eq. 7.1.7).  It applies to the
        ordinary electric dipole and, unchanged, to the nuclear-spin-independent PNC amplitude, both of which are electronic rank-1 tensors.
"""
function hyperfineDipoleFactor(Ja::AngularJ64, Fa::AngularJ64, Jb::AngularJ64, Fb::AngularJ64, spinI::AngularJ64)
    # a level whose J cannot couple with I to the given F is not part of this hyperfine manifold and contributes nothing; the test is made
    # on the doubled momenta because AngularMomentum.isTriangle RAISES rather than returns false when they cannot be coupled at all
    if  rem(Basics.twice(Ja) + Basics.twice(spinI) + Basics.twice(Fa), 2) != 0                       return( 0. )    end
    if  rem(Basics.twice(Jb) + Basics.twice(spinI) + Basics.twice(Fb), 2) != 0                       return( 0. )    end
    if  rem(Basics.twice(Fa) + Basics.twice(Fb), 2) != 0                                             return( 0. )    end
    if  !AngularMomentum.isTriangle(Fa, AngularJ64(1), Fb)                                           return( 0. )    end
    if  !AngularMomentum.isTriangle(Ja, spinI, Fa)  ||  !AngularMomentum.isTriangle(Jb, spinI, Fb)   return( 0. )    end
    n2 = Basics.twice(Ja) + Basics.twice(spinI) + Basics.twice(Fb) + 2
    wa = (-1.)^div(n2, 2) * sqrt( (Basics.twice(Fa) + 1.0) * (Basics.twice(Fb) + 1.0) ) *
         AngularMomentum.Wigner_6j(Ja, Fa, spinI, Fb, Jb, AngularJ64(1))

    return( wa )
end


"""
`WeakInteractionEnhancement.hyperfineScalarFactor(Ja::AngularJ64, Jb::AngularJ64, F::AngularJ64, spinI::AngularJ64)`
    ... to compute the geometrical factor that turns the ELECTRONIC rank-1 reduced matrix element <a J_a || alpha rho || b J_b> of the
        anapole operator into the matrix element of the full scalar Hamiltonian `H = (G_F/sqrt2) kappa (alpha.I/I) rho(r)` between the
        hyperfine states |a J_a I F> and |b J_b I F>.  A value::Float64 is returned; it is zero for a nucleus of spin zero.

            factor = (-1)^(I + J_b + F) { F  J_a  I ; 1  I  J_b }  sqrt( I(I+1)(2I+1) ) / I

        The first three factors are the reduction of a scalar product of two commuting tensors (Edmonds, Eq. 7.1.6), in the SAME form and
        with the same phase convention as `Hfs.computeHyperfineRepresentation` uses for the hyperfine matrix itself, which was verified
        there against the Casimir formula.  The last factor is the nuclear reduced matrix element <I||I||I> = sqrt(I(I+1)(2I+1)), divided by
        I because the anapole constant kappa is conventionally defined with the unit vector I/I; carrying that 1/I here is what makes a
        kappa extracted from a computed amplitude comparable with the published caesium value.

        Note that the operator is diagonal in F but NOT in J: the triangle |J_a - J_b| <= 1 is what lets the anapole reach intermediate
        levels that the weak charge cannot, and is the whole reason a spin-dependent amplitude survives where the spin-independent one
        vanishes.
"""
function hyperfineScalarFactor(Ja::AngularJ64, Jb::AngularJ64, F::AngularJ64, spinI::AngularJ64)
    floatI = Basics.twice(spinI) / 2.0
    if  floatI == 0.                                                                                 return( 0. )    end
    if  rem(Basics.twice(Ja) + Basics.twice(Jb), 2) != 0                                             return( 0. )    end
    if  rem(Basics.twice(Ja) + Basics.twice(spinI) + Basics.twice(F), 2) != 0                        return( 0. )    end
    if  rem(Basics.twice(Jb) + Basics.twice(spinI) + Basics.twice(F), 2) != 0                        return( 0. )    end
    if  !AngularMomentum.isTriangle(Ja, AngularJ64(1), Jb)                                           return( 0. )    end
    if  !AngularMomentum.isTriangle(Ja, spinI, F)  ||  !AngularMomentum.isTriangle(Jb, spinI, F)     return( 0. )    end
    n2 = Basics.twice(spinI) + Basics.twice(Jb) + Basics.twice(F)
    wa = (-1.)^div(n2, 2) * AngularMomentum.Wigner_6j(F, Ja, spinI, AngularJ64(1), spinI, Jb) *
         sqrt( floatI * (floatI + 1.0) * (2floatI + 1.0) ) / floatI

    return( wa )
end


"""
`WeakInteractionEnhancement.radialIntegralPPminus(weight::Array{Float64,1}, a::Orbital, b::Orbital, grid::Radial.Grid)`
    ... to compute INT weight(r) [P_a P_b - Q_a Q_b] dr, the combination that the Dirac matrix beta produces by weighting the large and the
        small components with opposite signs; a value::Float64 is returned.

        It completes the family of `WeakInteractionMoment.radialIntegralPPplus`, `radialIntegralPQminus` and `radialIntegralPQ`, and belongs
        beside them; it sits here only because the EDM operator does, see the module docstring.
"""
function radialIntegralPPminus(weight::Array{Float64,1}, a::Orbital, b::Orbital, grid::Radial.Grid)
    mtp = min( size(a.P, 1), size(b.P, 1), length(weight) )
    wa  = 0.
    for  i = 2:mtp   wa = wa + weight[i] * (a.P[i]*b.P[i] - a.Q[i]*b.Q[i]) * grid.wr[i]    end

    return( wa )
end

end # module
