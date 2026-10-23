
"""
`module  JAC.InternalConversion`
    ... a submodel of JAC that contains all methods for computing internal conversion coefficients
        (ICC) and rates: the ejection of a bound atomic electron, driven by a nuclear gamma-ray
        transition of given energy and multipolarity, into the continuum.

        This is the TRANSITION/PROCESS side of a bound-electron-to-continuum-electron process (two
        configurations: `initialConfigs`, the atom/ion before conversion, and `finalConfigs`, the ion
        after the converted electron is ejected) -- architecturally identical in topology to
        `AutoIonization` (bound N-electron state -> one continuum channel of given kappa, no external
        photon field, single fixed transition energy, sum over allowed continuum kappa), NOT to
        `PhotoEmission`/`PhotoIonization` (no photon-energy sweep, no gauge -- there is no "gauge"
        concept at all for this near-field nuclear interaction). The actual one-electron matrix
        element, however, is a genuine rank-L ONE-PARTICLE multipole operator (like
        `PhotoEmission.amplitude`'s `SpinAngular.OneParticleOperator`), not `AutoIonization`'s
        two-particle Coulomb operator -- so this module combines `AutoIonization`'s pipeline shape
        with `PhotoEmission`'s one-particle spin-angular machinery.

        CORE FORMULAS FOLLOW Rysavy & Dragoun, Report NPI Rez-TH-03/2000 (arXiv:nucl-ex/0006009), whose
        Eqs. (1)-(8) were read verbatim from the source and are the equation numbers used throughout this
        module. **The citation was WRONG until 23-Oct-2026 and is corrected here**: every "Raman Eq. ..." in
        this file pointed at Raman, Nestor, Ichihara & Trzhaskovskaya, Phys. Rev. C 66, 044312 (2002), which
        is the MEASUREMENT this module is compared against (its Table VII gives the 207Pb alpha_K) and not
        the source of the formalism. Nothing in the code changed with the correction; only the labels.

          Eq. (1)      alpha_i(tau L) = alpha*pi*omega * (2j_i+1)/(L(L+1)) * sum_kappa_f |C * (R + T)|^2,
                       with C the Clebsch-Gordan C^(j_f,-1/2)_(j_i,-1/2;L,0);
          Eqs. (2)/(3) the magnetic and electric radial factors R^(ML) and R^(EL);
          Eqs. (4)/(5)/(6) the bilinear orbital combinations F^(0), F and F^(-1) those integrate;
          Eqs. (7)/(8) the intranuclear PENETRATION term T, which this module does NOT implement --
                       the "no-penetration" (Rose) model, X_Lambda(kr) = h_Lambda(kr) everywhere. Section 3.7
                       of that report puts T's effect at 1e-2 % to 10 % for non-hindered transitions, so it is
                       a real but bounded omission and cannot account for a large discrepancy.

        In this module's own naming, the summand of Eq. (1) is split as B_if^(tau L) * R_if^(tau L): B_if a
        Clebsch-Gordan x Racah-W geometric factor, R_if the radial factor of Eq. (2) or (3). Because the kernel
        is a spherical HANKEL function (h_L = j_L + i*y_L), reflecting the ejected electron's outgoing-wave
        boundary condition, R_if -- and hence InternalConversion.Channel.amplitude -- is genuinely COMPLEX;
        not because of any photon helicity/gauge ambiguity, but because of that outgoing-wave near-field
        character.

        Note on validation (28/29-Jul-2026): the B_if^(tau L) angular factor was validated by implementing the
        published angular factor -- the Clebsch-Gordan of Eq. (1) together with the standard Racah-W recoupling
        -- literally with JAC's own AngularMomentum.ClebschGordan/Wigner_6j and comparing numerically, across
        92 (kappa_i,kappa_f,L) combinations, against AngularMomentum.CL_reduced_me: the magnitude ratio collapsed
        EXACTLY to 1/sqrt(L(L+1))
        (independent of kappa_i,kappa_f) in every case -- strong evidence the literal transcription
        (Clebsch-Gordan order, and the standard Racah W(j1j2j3j4;j5j6) = (-1)^(j1+j2+j3+j4) *
        {j1 j2 j5; j4 j3 j6} relation to the Wigner 6-j symbol) is correct. The residual per-(kappa_i,
        kappa_f) SIGN could not be matched to any of 14 simple candidate phase formulas against it -- and at
        the time that function's source carried an acknowledged open phase ambiguity ("this factor is not
        clear"), so this module uses its own literal, magnitude-validated B_if implementation directly
        rather than forcing a match to that convention.

        REDONE 10-Aug-2026, after the three CL_reduced_me variants were consolidated into one, the empirical
        phase comment removed and the missing parity selection rule added -- the premise of the comparison
        above had changed, so it was repeated rather than merely re-flagged. Over kappa_i, kappa_f in
        {-1,+1,-2,+2,-3,+3,-4} and L = 1..3: of the 55 non-trivial combinations the magnitude ratio is
        EXACTLY 1/sqrt(L(L+1)) in all 55, so that conclusion stands unchanged. The comparison is now in fact
        STRONGER: a further 92 combinations are parity-forbidden and both functions return zero for every one
        of them, with no case where one vanishes and the other does not. Before the parity fix CL_reduced_me
        returned spurious nonzero values on all 92, while reducedBif -- which has always applied its own
        parity check -- returned zero, so the two disagreed there. The SIGN still matches in only 28 of the
        55, i.e. the residual phase question is untouched by the consolidation, and this module's decision to
        use its own literal implementation stands.

        STAGE 1 (implemented here): single-active-electron, closed-shell-plus-one case (Eqs. (1)-(6) above)
        -- sufficient for the two literature anchors that examples/example-Dq.jl actually compares against,
        which are Raman et al. (207Pb, branch b) and Nica et al. (119Sn, branch c). STAGE 2 (explicitly NOT
        implemented): open-shell
        multi-electron spectator coupling (an extra 6-j recoupling factor, needed only for IC from
        excited/open-shell configurations as in Bilous et al., Phys. Rev. A 95, 032503 (2017)).

        STATUS, 23-Oct-2026: THE ABSOLUTE COEFFICIENTS NOW REPRODUCE TWO INDEPENDENT LITERATURE ANCHORS and the
        45-order discrepancy of the note below is explained. Whether the module counts as un-parked is the maintainer's
        call: the hard refusal that the 15-Aug note installed has been replaced by a printed note on scope, and he may
        put the refusal back.

        What was wrong -- the radial integral, exactly as the superseded note predicted it would be:
          * the sum over grid points STARTED AT THE FIRST POINT. The irregular kernel h_L(kr) ~ r^(-L-1) amplifies the
            innermost points instead of suppressing them, so it was multiplying the continuum orbital's own round-off
            floor; for the 207Pb M4 case h_4 = 2.9e+23 there against an orbital that was sign-oscillating noise. That
            single mechanism carried all 45 orders. The sum now starts at InternalConversion.icInnermostReliablePoint,
            which reads the floor off the orbitals themselves. A FIXED radius does not do it: the nuclear radius was
            measured to lie INSIDE the floor for 119Sn and OUTSIDE it for 207Pb -- the same multipole, either side of
            the only boundary that matters -- which is why an R_nuc cut looked stable on one case and not the other.
          * the Hankel argument was the ELECTRON momentum where Eq. (2) has the photon wavenumber omega/c, 374.2
            against 285.2 a.u. for 207Pb. h_L(omega r) is the field outside the nuclear current; the electron's own
            momentum has no place in it.
          * the prefactor alpha*pi*omega was missing altogether, and with it the (2L+1) of the reduced amplitude.

        Validated, computed value against the literature one:
          * 207Pb, 1063.6 keV M4: alpha_K = 0.09602 against the measured 0.0945(22) of Raman et al., Phys. Rev. C
            66, 044312 (2002), Table VII -- a ratio of 1.016, inside the measurement's own uncertainty band.
          * 119Sn, 65.660 keV M4: alpha_K = 1530.7 against Nica et al., Phys. Rev. C 89, 014303 (2014), Table IV,
            which quotes THREE calculations differing only in the treatment of the K vacancy: no vacancy 1544(1),
            SCF of the ion 1603(1), frozen orbitals 1618(1), against a measurement of 1621(25). We land 0.9 % from
            the NO-VACANCY value and 5.6 % below the measurement. THE VACANCY IS THE OPEN QUESTION: the continuum
            orbital here IS generated in the potential of the final, 1s-hole level, so the vacancy is present, yet
            its effect comes out under a per cent where that table puts it at +3.8 %. Isolating it needs the ion and
            the neutral potential compared at the SAME electron energy, which the perform route cannot express today
            because the choice lives inside Continuum.generateOrbitalForLevel. See examples/example-Dq.jl branch c.
          * both carry a +-2 % systematic from the inner match point, measured by moving it across its plateau
            rather than assumed, and both figures are what the SHIPPED pipeline returns rather than a hand assembly.
            A SECOND uncertainty of the same size comes from the grid, measured 23-Oct-2026 on the 119Sn case: the
            step refinement is settling (1530.7, 1529.7, 1519.3 for hp = 3.0e-3, 1.5e-3, 1.0e-3) but the box is NOT
            (1558.4, 1530.7, 1524.0 for rbox = 15, 20, 30, still drifting at 30). So neither anchor is converged
            better than one or two per cent, which is the same size as the vacancy effect above -- and that is why
            the vacancy question cannot be settled on these grids however the potential is chosen.

        WHAT THE TWO ANCHORS DO NOT SPEAK FOR, since both are K-shell MAGNETIC multipoles from a J = 0 closed-shell
        atom: an open-shell initial level, an L-shell coefficient, and any electric multipole. The (2L+1) was measured
        on that same closed-shell shape, so an open-shell case needs its own check before its numbers are used. The
        internal-consistency checks of example-Dq.jl branch a remain what they were and still say nothing about
        absolute values.
"""
module InternalConversion


using Printf, GSL, ..AngularMomentum, ..Basics, ..Continuum, ..Defaults, ..ManyElectron, ..Nuclear, ..Radial, ..SpinAngular, ..TableStrings


"""
`struct  InternalConversion.Settings  <:  Basics.AbstractProcessSettings`
    ... defines the settings for computing internal conversion coefficients between an initial
        (bound) and a final (ion + ejected electron) multiplet.

    + multipoles     ::Array{EmMultipole,1}   ... nuclear transition multipolarities. **ONE PER RUN: this
                                                   is NOT a mixed-multipolarity facility, and the docstring
                                                   claimed it was until 23-Oct-2026.** Every multipole given
                                                   has its channels summed INCOHERENTLY into one ICC, whereas
                                                   a genuinely mixed transition has
                                                   alpha = (alpha(M1) + delta^2 alpha(E2))/(1 + delta^2) with
                                                   the NUCLEAR mixing ratio delta, which this module neither
                                                   receives nor could know. So an array of more than one entry
                                                   returns a sum with no physical meaning; compute each
                                                   multipole separately and combine them where delta is known
                                                   (examples/example-Dq.jl branch d). No "gauge" concept
                                                   applies to this near-field nuclear interaction, unlike a
                                                   real photon field.
    + gammaEnergy    ::Float64                ... the nuclear transition energy (Hartree) -- an
                                                   explicit, externally-known input, NOT derived
                                                   from atomic level energies.
    + maxKappa       ::Int64                  ... maximum |kappa| of continuum partial waves to
                                                   sum over.
    + printBefore    ::Bool                   ... True if all lines are printed before evaluation.
    + lineSelection  ::LineSelection          ... Specifies the selected levels, if any.
    + acceptUnvalidated ::Bool                ... True silences the printed note on what the absolute
                                                   coefficients have been validated against; see the STATUS
                                                   note in this module's docstring for that scope.
"""
struct Settings  <:  Basics.AbstractProcessSettings
    multipoles      ::Array{EmMultipole,1}
    gammaEnergy     ::Float64
    maxKappa        ::Int64
    printBefore     ::Bool
    lineSelection   ::LineSelection
    acceptUnvalidated ::Bool
end


"""
`InternalConversion.Settings()`  ... constructor for the default values of internal-conversion computations.
"""
function Settings()
    Settings(EmMultipole[], 0., 100, false, LineSelection(), false)
end


"""
`InternalConversion.Settings(set::InternalConversion.Settings;`

        multipoles=.., gammaEnergy=.., maxKappa=.., printBefore=.., lineSelection=..,
        acceptUnvalidated=..)

    ... keyword copy-constructor for re-defining selected values of a settings::InternalConversion.Settings.
"""
function Settings(set::InternalConversion.Settings;
        multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,   gammaEnergy::Union{Nothing,Float64}=nothing,
        maxKappa::Union{Nothing,Int64}=nothing,                    printBefore::Union{Nothing,Bool}=nothing,
        lineSelection::Union{Nothing,LineSelection}=nothing,       acceptUnvalidated::Union{Nothing,Bool}=nothing)
    if  isnothing(multipoles)      multipolesx     = set.multipoles     else   multipolesx     = multipoles     end
    if  isnothing(gammaEnergy)     gammaEnergyx    = set.gammaEnergy    else   gammaEnergyx    = gammaEnergy    end
    if  isnothing(maxKappa)        maxKappax       = set.maxKappa       else   maxKappax       = maxKappa       end
    if  isnothing(printBefore)     printBeforex    = set.printBefore    else   printBeforex    = printBefore    end
    if  isnothing(lineSelection)   lineSelectionx  = set.lineSelection  else   lineSelectionx  = lineSelection  end
    if  isnothing(acceptUnvalidated)  acceptUnvalidatedx = set.acceptUnvalidated  else  acceptUnvalidatedx = acceptUnvalidated  end

    Settings( multipolesx, gammaEnergyx, maxKappax, printBeforex, lineSelectionx, acceptUnvalidatedx )
end


# `Base.show(io::IO, settings::InternalConversion.Settings)`  ... prepares a proper printout of settings::InternalConversion.Settings.
function Base.show(io::IO, settings::InternalConversion.Settings)
    println(io, "multipoles:               $(settings.multipoles)  ")
    println(io, "gammaEnergy:              $(settings.gammaEnergy)  ")
    println(io, "maxKappa:                 $(settings.maxKappa)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "lineSelection:            $(settings.lineSelection)  ")
    println(io, "acceptUnvalidated:        $(settings.acceptUnvalidated)  ")
end


"""
`struct  InternalConversion.Channel`
    ... one nuclear-multipolarity contribution to one internal-conversion line, mirroring
        PhotoIonization.PartialWave (a rank-L operator couples initialLevel's symmetry to an
        INTERMEDIATE symmetry `symmetry` -- the total J/parity of the final-ion-plus-continuum-
        electron coupled state -- before that intermediate symmetry connects on to finalLevel; this
        differs from AutoIonization.PartialWave, whose rank-0 operator makes the coupled-state symmetry
        trivially equal to the initial level's own symmetry).

    + multipole  ::EmMultipole      ... nuclear transition multipolarity (tau, L) of this contribution.
    + kappa      ::Int64            ... relativistic kappa of the ejected (continuum) electron.
    + symmetry   ::LevelSymmetry    ... total J/parity of the (final ion + continuum electron) coupled state.
    + amplitude  ::ComplexF64       ... the MANY-ELECTRON REDUCED amplitude of this channel, i.e. the CSF-pair
                                        sum of the one-particle B_if^(tau L) * R_if^(tau L) products of Eq. (1)
                                        weighted by the spin-angular coefficients and the mixing coefficients --
                                        NOT the bare one-particle product, as this entry said until 23-Oct-2026.
                                        Being REDUCED, it still needs the (2L+1) that
                                        InternalConversion.computeAmplitudesProperties applies before squaring.
                                        It is complex because the kernel h_L(omega r) is a spherical Hankel
                                        function (outgoing-wave near field), not because of any
                                        photon-helicity/gauge ambiguity.
"""
struct  Channel
    multipole   ::EmMultipole
    kappa       ::Int64
    symmetry    ::LevelSymmetry
    amplitude   ::ComplexF64
end


"""
`struct  InternalConversion.Line`
    ... one internal-conversion line between an initial (bound, N-electron) and final (ion + ejected
        electron, N-1-electron) level, mirroring AutoIonization.Line.

    + initialLevel    ::Level                                 ... initial bound-electron level.
    + finalLevel      ::Level                                 ... final ion level (electron vacancy).
    + electronEnergy  ::Float64                                ... E_k = gammaEnergy - bindingEnergy.
    + ICC             ::Float64                                ... total (multipole-summed) internal
                                                                     conversion coefficient alpha for this line.
    + channels        ::Array{InternalConversion.Channel,1}    ... per-multipole/kappa contributions.
"""
struct  Line
    initialLevel     ::Level
    finalLevel       ::Level
    electronEnergy   ::Float64
    ICC              ::Float64
    channels         ::Array{Channel,1}
end


"""
`InternalConversion.Line(initialLevel::Level, finalLevel::Level, electronEnergy::Float64)`
    ... constructor for an internal-conversion line with as-yet unevaluated channels/ICC.
"""
function Line(initialLevel::Level, finalLevel::Level, electronEnergy::Float64)
    Line(initialLevel, finalLevel, electronEnergy, 0., Channel[])
end


# `Base.show(io::IO, line::InternalConversion.Line)`  ... prepares a proper printout of line::InternalConversion.Line.
function Base.show(io::IO, line::InternalConversion.Line)
    println(io, "initialLevel:           $(line.initialLevel)  ")
    println(io, "finalLevel:             $(line.finalLevel)  ")
    println(io, "electronEnergy:         $(line.electronEnergy)  ")
    println(io, "ICC:                    $(line.ICC)  ")
    println(io, "channels:               $(line.channels)  ")
end


#################################################################################################################################
#################################################################################################################################
## Angular (B_if) and radial (R_if) machinery, Rysavy & Dragoun, arXiv:nucl-ex/0006009, Eqs. (1)-(6).


"""
`InternalConversion.kappaToLJ(kappa::Int64)`
    ... returns (l::Int64, j::AngularJ64) for a given relativistic quantum number kappa.
"""
function kappaToLJ(kappa::Int64)
    if  kappa > 0    l = kappa;        j = AngularJ64(2*kappa-1, 2)
    else             l = -kappa - 1;   j = AngularJ64(-2*kappa-1, 2)
    end
    return (l, j)
end


"""
`InternalConversion.racahW(j1::AngularJ64, j2::AngularJ64, j3::AngularJ64, j4::AngularJ64, j5::AngularJ64, j6::AngularJ64)`
    ... the Racah W-coefficient, related to the Wigner 6-j symbol by the standard
        (Rotenberg-Bivins-Metropolis-Wooten) relation
        W(j1 j2 j3 j4; j5 j6) = (-1)^(j1+j2+j3+j4) {j1 j2 j5; j4 j3 j6}.
"""
function racahW(j1::AngularJ64, j2::AngularJ64, j3::AngularJ64, j4::AngularJ64, j5::AngularJ64, j6::AngularJ64)
    phase = (-1.0)^Int64( Basics.twice(j1)/2 + Basics.twice(j2)/2 + Basics.twice(j3)/2 + Basics.twice(j4)/2 )
    return  phase * AngularMomentum.Wigner_6j(j1, j2, j5, j4, j3, j6)
end


"""
`InternalConversion.reducedBif(mp::EmMultipole, kappaInitial::Int64, kappaFinal::Int64)`
    ... computes the angular factor B_if^(tau L) for a nuclear multipole `mp` between a bound orbital of
        kappaInitial and a continuum orbital of kappaFinal: the Clebsch-Gordan of Eq. (1) times a Racah-W
        recoupling, magnitude-validated against AngularMomentum.CL_reduced_me as the module docstring records.
        For magnetic multipoles l_f is replaced by l_bar_f = 2*j_f - l_f, which is the standard device for
        writing the magnetic factor in the same Clebsch-Gordan form as the electric one; it is this module's
        own transcription and no equation number is claimed for it. The normalisation carried here INCLUDES
        the published (2j_i+1)/(L(L+1)) of Eq. (1), which is why
        InternalConversion.computeAmplitudesProperties must not apply it again. A value::Float64 is returned
        (real; only the radial factor R_if is complex).
"""
function reducedBif(mp::EmMultipole, kappaInitial::Int64, kappaFinal::Int64)
    li, ji = kappaToLJ(kappaInitial)
    lf, jf = kappaToLJ(kappaFinal)
    lUse   = mp.electric ? lf : (Basics.twice(jf) - lf)
    if  rem(li + lUse + mp.L, 2) != 0   return 0.0   end
    #
    CG    = AngularMomentum.ClebschGordan(AngularJ64(li), AngularM64(0), AngularJ64(lUse), AngularM64(0), AngularJ64(mp.L), AngularM64(0))
    W     = racahW(AngularJ64(li), ji, AngularJ64(lUse), jf, AngularJ64(1//2), AngularJ64(mp.L))
    phase = (-1.0)^Int64( Basics.twice(jf)/2 + 1//2 + mp.L )
    norm  = sqrt( (Basics.twice(ji)+1) * (2li+1) * (Basics.twice(jf)+1) * (2lUse+1) / (mp.L*(mp.L+1)*(2mp.L+1)) )
    return phase * CG * W * norm
end


"""
`InternalConversion.sphericalHankel1(L::Int64, x::Float64)`
    ... the spherical Hankel function of the first kind, h_L(x) = j_L(x) + i*y_L(x), built from
        GSL's regular (sf_bessel_jl) and irregular (sf_bessel_yl) spherical Bessel functions -- both
        already used elsewhere in JAC (module-RadialIntegrals.jl, module-InteractionStrength.jl's
        Breit density kernel). A value::ComplexF64 is returned.
"""
function sphericalHankel1(L::Int64, x::Float64)
    return  GSL.sf_bessel_jl(L, x) + im * GSL.sf_bessel_yl(L, x)
end


"""
`InternalConversion.icInnermostReliablePoint(a::Radial.Orbital, b::Radial.Orbital)`
    ... determines the innermost grid point from which a radial integral of this module may be accumulated, i.e. the first
        point at which both orbitals still carry signal rather than their own round-off floor. The distinction matters here
        and nowhere else in JAC: the irregular kernel h_L(kr) ~ r^(-L-1) AMPLIFIES whatever sits at the innermost points
        instead of suppressing it, so a continuum orbital that has fallen below its floor is magnified rather than
        neglected. Nothing physical is lost by starting further out, since P_a ~ r^(l_a+1) and P_b ~ r^(l_b+1) make the
        integrand rise from zero as a positive power of r; the omitted piece therefore scales as a power of the match
        radius and is negligible against an integral that is built at r ~ 1/Z. The floor is recognised from the orbitals
        alone -- it neither keeps its sign nor grows monotonically, whereas an analytic power law does both -- and NOT from
        any fixed radius such as the nuclear one, which was measured to fall on opposite sides of the floor for two cases
        of the same multipole. An index::Int64 is returned.
"""
function icInnermostReliablePoint(a::Radial.Orbital, b::Radial.Orbital)
    mtp    = min(size(a.P,1), size(b.P,1))
    window = 8
    for  m = 2:mtp-window
        if  icIsAnalyticFrom(a, m, window)  &&  icIsAnalyticFrom(b, m, window)   return( m )   end
    end

    return( 2 )
end


"""
`InternalConversion.icIsAnalyticFrom(orbital::Radial.Orbital, m::Int64, window::Int64)`
    ... tests whether `orbital` still follows its analytic small-r behaviour at grid point `m`, by requiring both the large
        and the small component to keep their sign and to grow in modulus over the next `window` points. A power law
        r^(l+1) satisfies both; a round-off floor, which is flat and sign-indifferent, satisfies neither.
        A value::Bool is returned.
"""
function icIsAnalyticFrom(orbital::Radial.Orbital, m::Int64, window::Int64)
    for  comp in (orbital.P, orbital.Q)
        if  m + window > size(comp,1)   return( false )   end
        sgn = sign(comp[m])
        if  sgn == 0.0                  return( false )   end
        for  mm = m:m+window-1
            if  sign(comp[mm]) != sgn  ||  abs(comp[mm+1]) <= abs(comp[mm])   return( false )   end
        end
    end

    return( true )
end


"""
`InternalConversion.icRadialIntegral(kind::Int64, L::Int64, k::Float64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid)`
    ... computes one of the three radial integrals out of which the R_if of Eqs. (2) and (3) are assembled;
        they are the bilinear orbital combinations F^(0), F and F^(-1) of Eqs. (4), (5) and (6) integrated
        against the Hankel kernel. This is the "no-penetration" (Rose) model, X_Lambda(kr) = h_Lambda(kr)
        everywhere, i.e. the T term of Eqs. (7)/(8) is omitted. JAC's orbital.P/.Q are used directly as the
        report's G(r) = r*g(r) and F(r) = r*f(r), its v and u respectively:

        kind=1: R_1 = int G_a(r) F_b(r) h_L(kr) dr
        kind=2: R_2 = int F_a(r) G_b(r) h_L(kr) dr
        kind=3: R_3 = int [G_a(r) G_b(r) + F_a(r) F_b(r)] h_L(kr) dr

        A value::ComplexF64 is returned. Mirrors RadialIntegrals.GrantILplus's structure, but with
        the irregular Hankel h_L(kr) (electron-continuum near field) in place of the regular photon
        Bessel j_L(qr).

        The sum starts at InternalConversion.icInnermostReliablePoint rather than at the first grid point: h_L(kr) grows
        as r^(-L-1) towards the origin and so amplifies the round-off floor of the continuum orbital, which is where the
        45 orders of magnitude of this module's old 207Pb M4 result came from. See there for why no fixed radius does.
"""
function icRadialIntegral(kind::Int64, L::Int64, k::Float64, a::Radial.Orbital, b::Radial.Orbital, grid::Radial.Grid)
    mtp = min(size(a.P,1), size(b.P,1))
    m0  = icInnermostReliablePoint(a, b)
    wa  = ComplexF64(0)
    for  m = m0:mtp
        if      kind == 1   fr = a.P[m]*b.Q[m]
        elseif  kind == 2   fr = a.Q[m]*b.P[m]
        else                fr = a.P[m]*b.P[m] + a.Q[m]*b.Q[m]
        end
        wa = wa + fr * sphericalHankel1(L, k*grid.r[m]) * grid.wr[m]
    end

    return( wa )
end


"""
`InternalConversion.electronMomentum(energy::Float64)`
    ... the relativistic electron momentum (atomic units) for kinetic energy `energy` (Hartree),
        q = sqrt(energy*(energy + 2c^2))/c, the same expression Continuum uses when it generates a continuum
        orbital. A value::Float64 is returned.

        **NOTHING IN THIS MODULE CALLS IT any more, as of 23-Oct-2026**, and that is deliberate rather than an
        oversight: the Hankel kernel takes the PHOTON wavenumber omega/c, and using the electron momentum here
        was one of the three defects repaired that day. It is kept only so that the distinction stays visible,
        and whether to retire it is the maintainer's call -- but it must NOT be reintroduced into
        InternalConversion.amplitude.
"""
function electronMomentum(energy::Float64)
    wc = Defaults.getDefaults("speed of light: c")
    return  sqrt( energy*(energy + 2*wc^2) ) / wc
end


"""
`InternalConversion.reducedRif(mp::EmMultipole, kappaInitial::Int64, kappaFinal::Int64, k::Float64, orbitalInitial::Radial.Orbital, orbitalFinal::Radial.Orbital, grid::Radial.Grid)`
    ... computes the radial factor R_if^(tau L) from the InternalConversion.icRadialIntegral building blocks:
        Eq. (2) for a magnetic multipole, which is the single integral of F^(0), and Eq. (3) for an electric
        one, which is the two-term combination carrying h_(L-1) beside h_L. A value::ComplexF64 is returned.
"""
function reducedRif(mp::EmMultipole, kappaInitial::Int64, kappaFinal::Int64, k::Float64,
                     orbitalInitial::Radial.Orbital, orbitalFinal::Radial.Orbital, grid::Radial.Grid)
    L = mp.L
    if  mp.electric
        # Lambda = L-1 (which is 0 when L==1 -- that h_0 branch still exists, not omitted)
        R1m = icRadialIntegral(1, L-1, k, orbitalInitial, orbitalFinal, grid)
        R2m = icRadialIntegral(2, L-1, k, orbitalInitial, orbitalFinal, grid)
        R3  = icRadialIntegral(3, L, k, orbitalInitial, orbitalFinal, grid)
        return  (kappaInitial-kappaFinal)*(R1m+R2m) + L*(R2m-R1m+R3)
    else
        R1 = icRadialIntegral(1, L, k, orbitalInitial, orbitalFinal, grid)
        R2 = icRadialIntegral(2, L, k, orbitalInitial, orbitalFinal, grid)
        return  (kappaInitial+kappaFinal)*(R1+R2)
    end
end


#################################################################################################################################
#################################################################################################################################
## Many-electron amplitude and standard process pipeline (mirrors AutoIonization.jl's shape,
## PhotoEmission.jl's one-particle spin-angular machinery).


"""
`InternalConversion.amplitude(mp::EmMultipole, gammaEnergy::Float64, finalLevel::Level, initialLevel::Level, grid::Radial.Grid)`
    ... computes the many-electron reduced amplitude for one requested nuclear multipole `mp`
        between `finalLevel` (the ion, with the continuum orbital of a specific kappa already
        attached as an extra "subshell" -- mirroring AutoIonization.amplitude/computeAmplitudesProperties)
        and `initialLevel` (the bound N-electron level), by summing the one-particle contributions
        InternalConversion.reducedBif(...) * InternalConversion.reducedRif(...) over CSF pairs via
        SpinAngular.OneParticleOperator(mp.L,...), exactly as PhotoEmission.amplitude does for a real
        photon field (including the same sqrt(2 j_a+1) "undo" of SpinAngular's internal GRASP-like
        convention, cf. the note in Hfs.amplitude). A value::ComplexF64 is returned.

        The argument of the Hankel kernel is the PHOTON wavenumber omega/c and not the electron momentum: h_L(omega r)
        is the multipole field outside the nuclear transition current, so the electron's own momentum never enters it.
        The two differ by 30 % for a 1 MeV transition, and the electron energy is still needed -- but only to decide
        whether a subshell is open at all.
"""
function amplitude(mp::EmMultipole, gammaEnergy::Float64, finalLevel::Level, initialLevel::Level, grid::Radial.Grid)
    nf = length(finalLevel.basis.csfs);   ni = length(initialLevel.basis.csfs)
    matrix = zeros(ComplexF64, nf, ni)
    opa    = SpinAngular.OneParticleOperator(mp.L, Basics.multipoleParity(mp))
    for  r = 1:nf
        for  s = 1:ni
            if  initialLevel.basis.csfs[s].J != initialLevel.J  ||  initialLevel.basis.csfs[s].parity != initialLevel.parity   continue    end
            subshellList = finalLevel.basis.subshells
            wa = SpinAngular.computeCoefficients(opa, finalLevel.basis.csfs[r], initialLevel.basis.csfs[s], subshellList)
            me = ComplexF64(0)
            for  coeff in wa
                orbFinal    = finalLevel.basis.orbitals[coeff.a]     # bra/final side (may be the continuum orbital)
                orbInitial  = initialLevel.basis.orbitals[coeff.b]   # ket/initial side (bound orbital)
                ja2         = Basics.subshell_2j(coeff.a)
                bindingEnergy = -orbInitial.energy
                electronEnergy = gammaEnergy - bindingEnergy
                if  electronEnergy <= 0.   continue    end
                k    = gammaEnergy / Defaults.getDefaults("speed of light: c")   ## the PHOTON wavenumber, not the electron's
                Bif  = reducedBif(mp, orbInitial.subshell.kappa, orbFinal.subshell.kappa)
                if  Bif == 0.   continue    end
                Rif  = reducedRif(mp, orbInitial.subshell.kappa, orbFinal.subshell.kappa, k, orbInitial, orbFinal, grid)
                me   = me + coeff.T / sqrt(ja2 + 1) * Bif * Rif
            end
            matrix[r,s] = me
        end
    end
    return  transpose(finalLevel.mc) * matrix * initialLevel.mc
end


"""
`InternalConversion.determineChannels(finalLevel::Level, initialLevel::Level, settings::InternalConversion.Settings)`
    ... determines the allowed (multipole, kappa) channels for one level pair, mirroring
        AutoIonization.determineChannels but generalized to rank-L multipole coupling (mirroring
        PhotoIonization.determineChannels's use of AngularMomentum.allowedMultipoleSymmetries +
        AngularMomentum.allowedKappaSymmetries, since a rank-L operator -- unlike AutoIonization's
        rank-0 Coulomb interaction -- couples initialLevel's symmetry to an INTERMEDIATE symmetry
        first). An Array{InternalConversion.Channel,1} is returned (amplitudes not yet evaluated).
"""
function determineChannels(finalLevel::Level, initialLevel::Level, settings::InternalConversion.Settings)
    channels = InternalConversion.Channel[]
    symi = LevelSymmetry(initialLevel.J, initialLevel.parity);   symf = LevelSymmetry(finalLevel.J, finalLevel.parity)
    for  mp in settings.multipoles
        symtList = AngularMomentum.allowedMultipoleSymmetries(symi, mp)
        for  symt in symtList
            kappaList = AngularMomentum.allowedKappaSymmetries(symt, symf)
            for  kappa in kappaList
                if  abs(kappa) > settings.maxKappa   continue    end
                push!(channels, InternalConversion.Channel(mp, kappa, symt, ComplexF64(0)))
            end
        end
    end
    return  channels
end


"""
`InternalConversion.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::InternalConversion.Settings)`
    ... forms the Cartesian product of initialMultiplet.levels x finalMultiplet.levels, filtered by
        settings.lineSelection, mirroring AutoIonization.determineLines/PhotoIonization.determineLines.
        The electron energy for a line is settings.gammaEnergy minus the energy difference between
        initial and final ATOMIC levels (the vacancy-creation energy), NOT simply an atomic level
        energy difference (unlike ordinary EM transitions) -- gammaEnergy is the externally-given
        nuclear transition energy. Lines with non-positive electron energy are skipped.
        An Array{InternalConversion.Line,1} is returned.
"""
function determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::InternalConversion.Settings)
    lines = InternalConversion.Line[]
    for  iLevel in initialMultiplet.levels
        for  fLevel in finalMultiplet.levels
            if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                electronEnergy = settings.gammaEnergy - (fLevel.energy - iLevel.energy)
                if  electronEnergy <= 0.   continue    end
                channels = InternalConversion.determineChannels(fLevel, iLevel, settings)
                if  length(channels) == 0   continue    end
                push!(lines, InternalConversion.Line(iLevel, fLevel, electronEnergy, 0., channels))
            end
        end
    end
    return  lines
end


"""
`InternalConversion.computeAmplitudesProperties(line::InternalConversion.Line, nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64, settings::InternalConversion.Settings)`
    ... fills in a Line's channels and total ICC (Eq. (1), summed over multipoles and kappa),
        mirroring AutoIonization.computeAmplitudesProperties: generates the continuum orbital for
        each channel's kappa via Continuum.generateOrbitalForLevel(...), attaches it to the final
        level, and calls InternalConversion.amplitude(...). A line::InternalConversion.Line is returned.

        Two factors multiply the squared amplitude here and both are easy to apply twice by mistake. alpha*pi*omega is
        the prefactor of Eq. (1) and is the WHOLE of it: the published (2j_i+1)/(L(L+1)) is already carried inside
        InternalConversion.reducedBif, which was measured term by term against the published angular factor. (2L+1)
        undoes the Wigner-Eckart normalisation of the rank-L REDUCED amplitude that InternalConversion.amplitude
        returns; measured as exactly 1/(2L+1) over ten channels spanning M1, E2, M2, E3 and M3, each to six digits.
"""
function computeAmplitudesProperties(line::InternalConversion.Line, nm::Nuclear.Model, grid::Radial.Grid, nrContinuum::Int64,
                                      settings::InternalConversion.Settings)
    newChannels = InternalConversion.Channel[];   contSettings = Continuum.Settings(false, nrContinuum);   icc = 0.
    alfa        = Defaults.getDefaults("alpha")
    subshellList = Basics.generate(OrderedSubshellList(), line.finalLevel.basis, line.initialLevel.basis)
    Defaults.setDefaults("relativistic subshell list", subshellList; printout=false)
    #
    for  channel in line.channels
        newiLevel = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, subshellList)
        newfLevel = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel,   subshellList)
        newiLevel = Basics.generateLevelWithExtraSubshell(Subshell(101, channel.kappa), newiLevel)
        cOrbital, _ = Continuum.generateOrbitalForLevel(line.electronEnergy, Subshell(101, channel.kappa), newfLevel, nm, grid, contSettings)
        newcLevel = Basics.generateLevelWithExtraElectron(cOrbital, channel.symmetry, newfLevel)
        amp = InternalConversion.amplitude(channel.multipole, settings.gammaEnergy, newcLevel, newiLevel, grid)
        icc = icc + alfa * pi * settings.gammaEnergy * (2*channel.multipole.L + 1) * real(conj(amp)*amp)
        push!(newChannels, InternalConversion.Channel(channel.multipole, channel.kappa, channel.symmetry, amp))
    end
    return  InternalConversion.Line(line.initialLevel, line.finalLevel, line.electronEnergy, icc, newChannels)
end


"""
`InternalConversion.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::InternalConversion.Settings)`
    ... the top-level driver, called automatically by Basics.perform(::Atomic.Computation) once
        computation.processSettings isa InternalConversion.Settings. `nm` is the nuclear model of the
        computation (Z, charge distribution) -- it fixes the potential in which the ejected electron's
        continuum orbital is generated, exactly as for AutoIonization/PhotoIonization; only the nuclear
        TRANSITION energy (settings.gammaEnergy) is external, not the nuclear charge model itself.
        An Array{InternalConversion.Line,1} is returned.
"""
function computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::InternalConversion.Settings)
    if  !settings.acceptUnvalidated
        # A note on scope rather than a refusal; see the STATUS block in the module docstring for the measurements.
        printstyled("\n>>> InternalConversion: the absolute coefficients are validated for K-shell MAGNETIC "       *
                    "multipoles of a J = 0 closed-shell atom only --\n"                                              *
                    "    207Pb 1063.6 keV M4: alpha_K = 0.09602 against the measured 0.0945(22), inside the "         *
                    "measurement's own band;\n"                                                                      *
                    "    119Sn 65.660 keV M4: alpha_K = 1530.7, 0.9 % from the NO-VACANCY 1544(1) but 5.6 % below "   *
                    "the measured 1621(25),\n"                                                                       *
                    "    the gap being the K vacancy, whose effect comes out here far smaller than published. An "    *
                    "open-shell initial level, an\n"                                                                 *
                    "    L-shell coefficient and any electric multipole are UNTESTED against the literature. See "    *
                    "examples/example-Dq.jl, branches b-e.\n"                                                        *
                    "    Set  acceptUnvalidated = true  to silence this note.\n", color=:cyan)
    end
    println("")
    printstyled("InternalConversion.computeLines(): The computation of internal conversion coefficients starts now ... \n", color=:light_green)
    printstyled("----------------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    lines = InternalConversion.determineLines(finalMultiplet, initialMultiplet, settings)
    if  settings.printBefore    InternalConversion.displayLines(stdout, lines)    end
    maxEnergy = 0.;   for line in lines   maxEnergy = max(maxEnergy, line.electronEnergy)   end
    nrContinuum = Continuum.gridConsistency(maxEnergy, grid)
    newLines = InternalConversion.Line[]
    for  line in lines
        push!(newLines, InternalConversion.computeAmplitudesProperties(line, nm, grid, nrContinuum, settings))
    end
    InternalConversion.displayLines(stdout, newLines)
    return  newLines
end


"""
`InternalConversion.displayLines(stream::IO, lines::Array{InternalConversion.Line,1})`
    ... prints the computed ICC (and per-multipole/kappa breakdown) in a table, mirroring
        AutoIonization.displayRates. Nothing is returned.
"""
function displayLines(stream::IO, lines::Array{InternalConversion.Line,1})
    nx = 100
    println(stream, " ")
    println(stream, "  Internal conversion coefficients:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(18, "i-level-f"; na=2);                          sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(18, "i--J^P--f"; na=4);                          sb = sb * TableStrings.hBlank(22)
    sa = sa * TableStrings.center(16, "Electron energy"; na=2);
    sb = sb * TableStrings.center(16, TableStrings.inUnits("energy"); na=2)
    sa = sa * TableStrings.center(16, "ICC (alpha)"; na=4);                        sb = sb * TableStrings.hBlank(20)
    println(stream, sa);   println(stream, sb);   println(stream, "  ", TableStrings.hLine(nx))
    for  line in lines
        sa   = "  ";    isym = LevelSymmetry(line.initialLevel.J, line.initialLevel.parity)
                        fsym = LevelSymmetry(line.finalLevel.J,   line.finalLevel.parity)
        sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
        sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
        sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.electronEnergy)) * "    "
        sa = sa * @sprintf("%.6e", line.ICC) * "    "
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    return  nothing
end


end # module InternalConversion
