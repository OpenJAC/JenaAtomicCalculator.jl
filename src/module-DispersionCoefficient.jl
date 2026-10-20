
"""
`module  JAC.DispersionCoefficient`
... a submodel of JAC that turns multipole polarizabilities at IMAGINARY frequency into the long-range
    dispersion coefficients C_n of the atom-atom interaction V(R) = - C_6/R^6 - C_8/R^8 - C_10/R^10 - ... and of
    the atom-wall interaction U(z) = - C_3/z^3. It computes no structure of its own: every number here is a
    re-summation of the terms that `MultipolePolarizibility` already stores per perturber, and the module's whole
    content is the ANGULAR/COMBINATORIAL weight that turns those terms into a coefficient.

    THE CASIMIR-POLDER RELATION IS ONE BINOMIAL, not a table of special cases. With alpha^(E l,0)(i omega) the
    scalar 2^l-pole polarizability at imaginary frequency and n = l1 + l2 + 1,

        C_2n = (1/2 pi) sum over ORDERED pairs (l1,l2) with l1 + l2 = n - 1  of
                        binomial(2 l1 + 2 l2, 2 l1) * integral_0^inf alpha^(E l1,0)_A(i w) alpha^(E l2,0)_B(i w) dw

    which reproduces the four prefactors the literature writes out by hand: C_6 = (3/pi) int a1 a1,
    C_8 = (15/2 pi) [int a1 a2 + int a2 a1], and C_10 = (28/pi) int a1 a3 + (35/pi) int a2 a2 for identical
    partners. `DispersionCoefficient.casimirPolderPrefactor` is that weight and nothing else.

    AND NO NUMERICAL omega QUADRATURE IS NEEDED, which is what makes this module small. On a DISCRETE perturber
    spectrum each polarizability is a sum of single poles, alpha(i w) = sum_a A_a dE_a/(dE_a^2 + w^2), and

        integral_0^inf [a/(a^2+w^2)] [b/(b^2+w^2)] dw = pi / (2 (a+b))

    exactly. Every C_n is therefore a DOUBLE SUM over perturber pairs weighted by 1/(dE_a + dE_b) -- Dalgarno's
    oscillator-strength form -- and the only inputs are the per-perturber strengths A_a and energies dE_a that
    `MultipolePolarizibility.Contribution` keeps. `DispersionCoefficient.overlapIntegral` is that double sum.

    THE ATOM-WALL COEFFICIENT IS THE SAME INTEGRAL WITH ONE POLARIZABILITY INSTEAD OF TWO. For a perfectly
    conducting wall the non-retarded (Lennard-Jones) energy is U(z) = -C_3/z^3 with

        C_3 = (1/4 pi) integral_0^inf alpha^(E1,0)(i omega) d omega = (1/8) sum_a A_a

    the second form again being exact on a discrete spectrum. A wall of CONSTANT dielectric constant eps simply
    scales this by the reflection factor (eps-1)/(eps+1). The closed form is worth stating because it is
    independently checkable:
    with A_a = (2/3) |<a|d|0>|^2 the sum rule sum_a A_a = (2/3) <d^2> gives C_3 = <d^2>/12, which is exactly what
    the electrostatic image dipole yields, and for H(1s), <r^2> = 3, so C_3 = 1/4 EXACTLY.

    C_3 IS ALSO THE QUANTITY THAT A POOR PERTURBER SET DAMAGES FIRST, because that sum carries NO energy
    denominator at all, while alpha(0) = sum_a A_a/dE_a and every C_n weight 1/(dE_a + dE_b) suppress the high
    perturbers. Measured on ONE np n = 2..16 set of H(1s) at two radial boxes: on a 64 a.u. box alpha^(E1,0) is
    1.00000 of the reference and C_3 is 0.99995 of the exact 1/4, while on a 775 a.u. box they are 0.86054 and
    0.77026 -- the error is 1.7 times larger on C_3, and NOTHING in alpha(0) says so.
    `DispersionCoefficient.checkSaturation` is the guard for it.

    A FREQUENCY-DEPENDENT WALL COSTS NO QUADRATURE EITHER, AND IT SHOWS THAT THE CONSTANT-eps OPTION IS A LIMIT
    RATHER THAN A MODEL. One undamped Lorentz oscillator, eps(i omega) = 1 + omega_p^2/(omega_0^2 + omega^2), has a
    reflection factor that is ITSELF a single pole,

        r(i omega) = (eps-1)/(eps+1) = R/(Omega^2 + omega^2),   R = omega_p^2/2,   Omega^2 = omega_0^2 + omega_p^2/2

    which is the same shape as one polarizability term -- so the wall integral stays the closed form the rest of
    this module uses, and `DispersionCoefficient.wallReflectionPole` is that one line. Writing the oscillator by its
    STATIC dielectric constant eps(0) = 1 + omega_p^2/omega_0^2 and its resonance omega_0 gives

        C_3 = (R/8) sum_a A_a / [Omega (dE_a + Omega)],   Omega = omega_0 sqrt((eps(0)+1)/2),  R = (eps(0)-1) omega_0^2/2

    THE CONSTANT-eps ANSWER IS THE Omega >> dE LIMIT OF THAT, i.e. a wall whose electrons respond instantaneously on
    the atom's own scale, and the error of using it is one factor per pole, exactly 1 + dE_a/Omega. That is not a
    small correction for a real wall: a silica-like oscillator (eps(0) = 3.9, omega_0 = 10 eV) puts Omega = 0.575
    a.u. against dE = 0.375 a.u. for H(1s) -> 2p, so the constant-eps shortcut OVERSTATES C_3 by 65 % on that one
    pole and, MEASURED over the whole np n = 2..16 set, by 81 % -- 0.147952 against 0.081850 a.u., the ratio 0.5532
    (`tools/probe-dielectricWall.jl`, 20-Oct-2026). A perfect conductor needs no resonance at all, reflecting fully
    at every frequency, and is the default.

    AND THE RADIAL BOX IS WHAT DECIDES THAT, IN THE OPPOSITE DIRECTION TO RULE 12'S USUAL ADVICE. A sum over
    states needs the upper members of its perturber series to stand in for the CONTINUUM, and a box matched to
    the highest orbital of the series leaves them as genuine bound states instead -- so the set becomes
    bound-only and the continuum is missing entirely, not approximated badly. Measured on the same three series,
    `Basics.recommendedGrid` of the union gives 775-1342 a.u. and alpha^(E1,0), alpha^(E2,0), alpha^(E3,0) come
    out at 0.811, 0.349 and 0.085 of their references, which is exactly what a few BOUND perturbers give; a
    64-100 a.u. box gives 1.00000, 1.00000 and 1.00004. So choose the box small enough that the series runs into
    the continuum, and beware the other end too: a box that is far too large multiplies a high multipole instead
    of starving it, since the matrix elements weight r^L.

    THE ONE GUARD THAT IS NOT OPTIONAL: BOTH PARTNERS MUST BE IN THEIR LOWEST LEVEL. The derivation above
    integrates along the imaginary axis because all the denominators then have one sign. If a partner is EXCITED,
    at least one dE_a is negative, the contour argument fails, and the true interaction acquires an additional
    REAL-frequency (resonant) term -- see Lepers & Dulieu, arXiv:1703.02833, Eq. (75). Excited levels are exactly
    what JAC is for, so here that is the common case and not an edge one, and
    `DispersionCoefficient.checkGroundLevel` therefore RAISES rather than returning a plausible number.

    VALIDATED against Tang, Zhang, Zhang, Jiang & Mitroy, arXiv:1205.6883, which is the same paper that anchored
    the polarizabilities themselves: H(1s)-H(1s) C_6 = 6.498 392 245 754, C_8 = 124.384 030 658,
    C_10 = 3285.370 788 289 (Dirac). See `tools/probe-dispersionCoefficient.jl` for the measured numbers.
"""
module DispersionCoefficient


using Printf, ..Basics, ..Defaults, ..ManyElectron, ..MultipolePolarizibility, ..Nuclear, ..Radial, ..TableStrings


"""
`struct  DispersionCoefficient.MultipoleTerms`
    ... defines a type to keep ONE atom's complete 2^L-pole response as the list of single poles it is made of, which
        is the only input any dispersion coefficient needs. The scalar polarizability at imaginary frequency is
        alpha^(E L,0)(i w) = sum_a strengths[a] * energies[a] / (energies[a]^2 + w^2), and the static value is
        alpha^(E L,0)(0) = sum_a strengths[a] / energies[a].

    + multipole   ::EmMultipole        ... The multipole E L whose response these terms carry.
    + strengths   ::Array{Float64,1}   ... The weights A_a = c(L,0,J,J_a) |<J_a || O^(E L) || J>|^2 of each perturber,
                                           with c the rank-0 angular coefficient `MultipolePolarizibility.rankCoefficient`.
    + energies    ::Array{Float64,1}   ... The excitation energies dE_a = W_a - W_level (atomic units), all POSITIVE
                                           for a lowest level, which is what the Casimir-Polder form requires.
"""
struct MultipoleTerms
    multipole     ::EmMultipole
    strengths     ::Array{Float64,1}
    energies      ::Array{Float64,1}
end


# `Base.show(io::IO, terms::DispersionCoefficient.MultipoleTerms)`
#       ... prepares a proper printout of the variable terms::DispersionCoefficient.MultipoleTerms.
function Base.show(io::IO, terms::DispersionCoefficient.MultipoleTerms)
    println(io, "multipole:                 $(terms.multipole)  ")
    println(io, "number of poles:           $(length(terms.strengths))  ")
    println(io, "alpha(0):                  $(DispersionCoefficient.staticPolarizability(terms))  ")
end


"""
`struct  DispersionCoefficient.Partner`
    ... defines a type to keep one of the two atoms (or the single atom, for an atom-wall coefficient) of a dispersion
        calculation, together with its multipole response. A Partner is the unit that is checked by
        `DispersionCoefficient.checkGroundLevel`, since the Casimir-Polder form is valid only if its level is the
        lowest one its own perturber set can see.

    + name        ::String                           ... A label for the printout, e.g. "H(1s)"; it carries no physics.
    + level       ::Level                            ... The level whose dispersion coefficients are wanted.
    + terms       ::Array{DispersionCoefficient.MultipoleTerms,1}
                                                     ... One entry per multipole available for this atom; a C_n that
                                                         needs a multipole absent from this list is reported as
                                                         INCOMPLETE rather than silently computed too small.
"""
struct Partner
    name          ::String
    level         ::Level
    terms         ::Array{DispersionCoefficient.MultipoleTerms,1}
end


"""
`DispersionCoefficient.Partner(name::String, level::Level, gMultiplet::Multiplet, grid::Radial.Grid,
                                multipoles::Array{EmMultipole,1})`
    ... constructor that builds a partner::DispersionCoefficient.Partner directly from a level and a perturber
        multiplet, by asking `MultipolePolarizibility` for the terms of each requested multipole in turn; a
        `partner::DispersionCoefficient.Partner` is returned. This is the convenient path and uses ONE perturber set
        for every multipole. Where the multipoles need different sets -- and they do, since the discrete series
        saturates the more slowly the higher the multipole, and the radial box that suits one may ruin another -- build
        each `MultipoleTerms` separately with `DispersionCoefficient.multipoleTerms` and use the plain constructor
        above.
"""
function Partner(name::String, level::Level, gMultiplet::Multiplet, grid::Radial.Grid,
                 multipoles::Array{EmMultipole,1})
    terms = DispersionCoefficient.MultipoleTerms[]
    for  multipole in multipoles
        push!( terms, DispersionCoefficient.multipoleTerms(multipole, level, gMultiplet, grid) )
    end

    return( DispersionCoefficient.Partner(name, level, terms) )
end


"""
`DispersionCoefficient.Partner()`
    ... constructor for an `empty` instance of DispersionCoefficient.Partner. It carries no multipole response at all,
        and is the value that tells `DispersionCoefficient.computeOutcomes` to treat the pair as HOMONUCLEAR; a
        `partner::DispersionCoefficient.Partner` is returned.
"""
function Partner()
    Partner( "", Level(), DispersionCoefficient.MultipoleTerms[] )
end


# `Base.show(io::IO, partner::DispersionCoefficient.Partner)`
#       ... prepares a proper printout of the variable partner::DispersionCoefficient.Partner.
function Base.show(io::IO, partner::DispersionCoefficient.Partner)
    println(io, "name:                      $(partner.name)  ")
    println(io, "level:                     $(partner.level)  ")
    println(io, "multipoles available:      $([t.multipole  for t in partner.terms])  ")
end


"""
`struct  DispersionCoefficient.Coefficient`
    ... defines a type to keep ONE atom-atom dispersion coefficient C_n together with the (l1,l2) pieces it was built
        from, because the pieces are what a reader needs: they say which multipole response the coefficient actually
        rests on, and hence which one has to be converged.

    + n           ::Int64                                     ... The order, i.e. the n of C_n in -C_n/R^n; even and
                                                                  at least 6.
    + value       ::Float64                                   ... The coefficient itself, in atomic units.
    + pieces      ::Array{Tuple{Int64,Int64,Float64},1}       ... (l1, l2, contribution) for every ordered multipole
                                                                  pair with l1 + l2 = n/2 - 1, in the order computed.
    + complete    ::Bool                                      ... False if some (l1,l2) pair was skipped because a
                                                                  partner carried no terms for that multipole, in
                                                                  which case `value` is a LOWER bound and nothing
                                                                  more.
"""
struct Coefficient
    n             ::Int64
    value         ::Float64
    pieces        ::Array{Tuple{Int64,Int64,Float64},1}
    complete      ::Bool
end


# `Base.show(io::IO, coefficient::DispersionCoefficient.Coefficient)`
#       ... prepares a proper printout of the variable coefficient::DispersionCoefficient.Coefficient.
function Base.show(io::IO, coefficient::DispersionCoefficient.Coefficient)
    println(io, "C_$(coefficient.n):                     $(coefficient.value)  ")
    println(io, "pieces (l1,l2,value):      $(coefficient.pieces)  ")
    println(io, "complete:                  $(coefficient.complete)  ")
end


"""
`struct  DispersionCoefficient.Outcome`
    ... defines a type to keep the outcome of one dispersion-coefficient computation, i.e. of ONE pair of partners.

    + partnerA        ::DispersionCoefficient.Partner               ... The first atom.
    + partnerB        ::DispersionCoefficient.Partner               ... The second atom; identical to partnerA for a
                                                                        homonuclear pair.
    + coefficients    ::Array{DispersionCoefficient.Coefficient,1}  ... The atom-atom coefficients C_n requested, in
                                                                        ascending order of n.
    + atomWallC3      ::Float64                                     ... The atom-wall coefficient C_3 of partnerA in
                                                                        atomic units, or 0. if it was not requested.
    + wallDielectric  ::Float64                                     ... The STATIC dielectric constant the wall was
                                                                        given; Inf for a perfect conductor.
    + wallResonance   ::Float64                                     ... The resonance omega_0 of that wall's Lorentz
                                                                        oscillator, in atomic units; Inf for the
                                                                        instantaneous-response limit, in which the
                                                                        constant-eps reflection factor is exact.
"""
struct Outcome
    partnerA          ::DispersionCoefficient.Partner
    partnerB          ::DispersionCoefficient.Partner
    coefficients      ::Array{DispersionCoefficient.Coefficient,1}
    atomWallC3        ::Float64
    wallDielectric    ::Float64
    wallResonance     ::Float64
end


"""
`DispersionCoefficient.Outcome()`
    ... constructor for an `empty` instance of DispersionCoefficient.Outcome.
"""
function Outcome()
    Outcome( Partner("", Level(), DispersionCoefficient.MultipoleTerms[]),
             Partner("", Level(), DispersionCoefficient.MultipoleTerms[]),
             DispersionCoefficient.Coefficient[], 0., Inf, Inf )
end


# `Base.show(io::IO, outcome::DispersionCoefficient.Outcome)`
#       ... prepares a proper printout of the variable outcome::DispersionCoefficient.Outcome.
function Base.show(io::IO, outcome::DispersionCoefficient.Outcome)
    println(io, "partnerA:                  $(outcome.partnerA.name)  ")
    println(io, "partnerB:                  $(outcome.partnerB.name)  ")
    println(io, "coefficients:              $([ (c.n, c.value)  for c in outcome.coefficients ])  ")
    println(io, "atomWallC3:                $(outcome.atomWallC3)  ")
end


"""
`struct  DispersionCoefficient.Settings`
    ... defines a type for the details and parameters of computing long-range dispersion coefficients. It IS an
        `AbstractPropertySettings`, so the module is reached through `propertySettings` of an `Atomic.Computation`
        exactly as every other property module is, and the module's own functions remain available for the cases the
        property route cannot express. A dispersion coefficient is a property of a PAIR of atoms, and
        `perform(::Atomic.Computation)` carries only ONE multiplet -- so the computed atom is partner A, taken level
        by level from that multiplet, while partner B is carried IN THESE SETTINGS, already built. Leaving
        `partnerB` empty asks for a HOMONUCLEAR pair, which is what a published C_n of a single species means.

    + orders            ::Array{Int64,1}        ... The orders n of the atom-atom coefficients C_n wanted; even and at
                                                      least 6. The default [6, 8, 10] is what the literature quotes.
    + multipoles        ::Array{EmMultipole,1}  ... The multipoles whose response is to be built when a Partner is
                                                      constructed from a level and a perturber multiplet. C_6 needs
                                                      E1, C_8 needs E1 and E2, and C_10 needs E1, E2 and E3.
    + calcAtomWall      ::Bool                  ... True if the atom-wall coefficient C_3 of partnerA is wanted
                                                      beside the atom-atom ones.
    + wallDielectric    ::Float64               ... The STATIC dielectric constant eps(0) of that wall; Inf, the
                                                      default, is a perfect conductor at every frequency.
    + wallResonance     ::Float64               ... The resonance omega_0 of the wall's one Lorentz oscillator, in
                                                      atomic units, which fixes eps(i omega) = 1 + omega_p^2 /
                                                      (omega_0^2 + omega^2) together with wallDielectric. Inf, the
                                                      default, is the instantaneous-response limit, where the
                                                      constant-eps reflection factor (eps-1)/(eps+1) is exact; a
                                                      FINITE omega_0 reduces C_3 by 1 + dE_a/Omega per perturber,
                                                      which is not a small factor for a real wall -- see the module
                                                      docstring.
    + gMultiplet        ::Multiplet             ... The perturber multiplet of the atom whose levels are computed,
                                                      i.e. the intermediate ASF over which the response is summed.
                                                      As in `MultipolePolarizibility`, the caller generates it ahead
                                                      of time; see the module docstring.
    + partnerB          ::DispersionCoefficient.Partner
                                                ... The SECOND atom, fully built by the caller, including its own
                                                      perturber response. An empty partner -- the default -- means a
                                                      HOMONUCLEAR pair, and partner B is then taken to be partner A
                                                      itself, which is the common case and the one every published
                                                      C_n of a single species refers to.
    + printBefore       ::Bool                  ... True if the partners and their response are printed before the
                                                      coefficients are computed.
    + levelSelection    ::LevelSelection        ... Specifies the levels to be computed. The Casimir-Polder form
                                                      holds only for the LOWEST level a perturber set can see, so a
                                                      selection that admits an excited level makes
                                                      `DispersionCoefficient.checkGroundLevel` raise rather than
                                                      return a plausible number.
"""
struct Settings <: Basics.AbstractPropertySettings
    orders              ::Array{Int64,1}
    multipoles          ::Array{EmMultipole,1}
    calcAtomWall        ::Bool
    wallDielectric      ::Float64
    wallResonance       ::Float64
    gMultiplet          ::Multiplet
    partnerB            ::DispersionCoefficient.Partner
    printBefore         ::Bool
    levelSelection      ::LevelSelection
end


"""
`DispersionCoefficient.Settings()`
    ... constructor for an `empty` instance of DispersionCoefficient.Settings.
"""
function Settings()
    Settings( [6, 8, 10], EmMultipole[E1, E2, E3], false, Inf, Inf, Multiplet(), DispersionCoefficient.Partner(),
              false, LevelSelection() )
end


"""
`DispersionCoefficient.Settings(set::DispersionCoefficient.Settings;`

        orders=.., multipoles=.., calcAtomWall=.., wallDielectric=.., wallResonance=.., gMultiplet=..,
        partnerB=.., printBefore=.., levelSelection=..)

    ... keyword copy-constructor for re-defining selected values of a settings::DispersionCoefficient.Settings.
"""
function Settings(set::DispersionCoefficient.Settings;
        orders::Union{Nothing,Array{Int64,1}}=nothing,       multipoles::Union{Nothing,Array{EmMultipole,1}}=nothing,
        calcAtomWall::Union{Nothing,Bool}=nothing,           wallDielectric::Union{Nothing,Float64}=nothing,
        wallResonance::Union{Nothing,Float64}=nothing,       gMultiplet::Union{Nothing,Multiplet}=nothing,
        partnerB::Union{Nothing,DispersionCoefficient.Partner}=nothing,
        printBefore::Union{Nothing,Bool}=nothing,            levelSelection::Union{Nothing,LevelSelection}=nothing)
    if  isnothing(orders)           ordersx         = set.orders         else   ordersx         = orders         end
    if  isnothing(multipoles)       multipolesx     = set.multipoles     else   multipolesx     = multipoles     end
    if  isnothing(calcAtomWall)     calcAtomWallx   = set.calcAtomWall   else   calcAtomWallx   = calcAtomWall   end
    if  isnothing(wallDielectric)   wallDielectricx = set.wallDielectric else   wallDielectricx = wallDielectric end
    if  isnothing(wallResonance)    wallResonancex  = set.wallResonance  else   wallResonancex  = wallResonance  end
    if  isnothing(gMultiplet)       gMultipletx     = set.gMultiplet     else   gMultipletx     = gMultiplet     end
    if  isnothing(partnerB)         partnerBx       = set.partnerB       else   partnerBx       = partnerB       end
    if  isnothing(printBefore)      printBeforex    = set.printBefore    else   printBeforex    = printBefore    end
    if  isnothing(levelSelection)   levelSelectionx = set.levelSelection else   levelSelectionx = levelSelection end

    Settings( ordersx, multipolesx, calcAtomWallx, wallDielectricx, wallResonancex, gMultipletx, partnerBx,
              printBeforex, levelSelectionx )
end


# `Base.show(io::IO, settings::DispersionCoefficient.Settings)`
#       ... prepares a proper printout of the variable settings::DispersionCoefficient.Settings.
function Base.show(io::IO, settings::DispersionCoefficient.Settings)
    println(io, "orders:                   $(settings.orders)  ")
    println(io, "multipoles:               $(settings.multipoles)  ")
    println(io, "calcAtomWall:             $(settings.calcAtomWall)  ")
    println(io, "wallDielectric:           $(settings.wallDielectric)  ")
    println(io, "wallResonance:            $(settings.wallResonance)  ")
    println(io, "gMultiplet:               $(settings.gMultiplet)  ")
    println(io, "partnerB:                 $(settings.partnerB.name)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "levelSelection:           $(settings.levelSelection)  ")
end


#####################################################################################################################
#####################################################################################################################


"""
`DispersionCoefficient.casimirPolderPrefactor(l1::Int64, l2::Int64)`
    ... returns the weight with which the overlap integral of a 2^l1-pole and a 2^l2-pole polarizability enters the
        dispersion coefficient of order n = 2(l1 + l2 + 1), i.e. binomial(2 l1 + 2 l2, 2 l1)/(2 pi); a
        `prefactor::Float64` is returned. This single expression is the whole of the Casimir-Polder combinatorics and
        reproduces every prefactor the literature writes out: summed over ORDERED pairs it gives 3/pi for C_6,
        15/2 pi for each of the two C_8 terms, and 28/pi (l1,l2 = 1,3 and 3,1 together) and 35/pi (2,2) for C_10.
"""
function casimirPolderPrefactor(l1::Int64, l2::Int64)
    return( binomial(2*l1 + 2*l2, 2*l1) / (2*pi) )
end


"""
`DispersionCoefficient.checkGroundLevel(partner::DispersionCoefficient.Partner)`
    ... verifies that `partner`'s level is the LOWEST one its own perturber set can reach, which is the condition the
        whole Casimir-Polder construction rests on, and RAISES an informative error if it is not; nothing is returned.

        THIS IS THE ONE GUARD THAT IS NOT OPTIONAL, and it refuses rather than warns because the number it would
        otherwise return is finite, smooth and wrong. The imaginary-frequency representation is reached by a contour
        argument that needs every excitation energy dE_a = W_a - W_level to have the SAME sign. For an excited level at
        least one perturber lies BELOW it, that dE_a is negative, and two things happen at once: the integrand acquires
        a pole on the contour, and the true interaction gains an additional REAL-frequency (resonant) term that no
        imaginary-axis integral contains -- Lepers & Dulieu, arXiv:1703.02833, Eq. (75). Since excited levels are
        exactly what JAC is for, this is the common case here and not an edge one.
"""
function checkGroundLevel(partner::DispersionCoefficient.Partner)
    for  terms in partner.terms
        for  (a, dE) in enumerate(terms.energies)
            if  dE <= 0.0
                error("\nDispersionCoefficient: partner '$(partner.name)' is NOT in its lowest level -- its " *
                      "$(terms.multipole) perturber no. $a lies $(abs(dE)) a.u. BELOW it.\n" *
                      "A Casimir-Polder coefficient is NOT DEFINED by the imaginary-frequency integral in that " *
                      "case: the contour argument needs every excitation energy to carry the same sign, and an " *
                      "excited partner adds a REAL-frequency (resonant) term that this integral does not contain " *
                      "(Lepers & Dulieu, arXiv:1703.02833, Eq. 75).\n" *
                      "This is refused rather than warned about because the number would be finite, smooth and " *
                      "wrong. Use the lowest level of each partner, or extend the module by that term.")
            end
        end
    end

    return( nothing )
end


"""
`DispersionCoefficient.checkSaturation(partner::DispersionCoefficient.Partner)`
    ... inspects whether each of `partner`'s multipole responses rests on a perturber set that reaches into the
        CONTINUUM, and WARNS LOUDLY, to stdout, where it does not; nothing is returned.

        THE CONDITION BEING TESTED IS THE ONE THAT MAKES A SUM A LOWER BOUND. A bound Rydberg series accumulates at
        the ionization threshold, so its excitation energies crowd into a narrow band and its highest member lies
        barely above its lowest; a set that continues into the continuum, as a box-confined pseudo-state series does,
        spreads over orders of magnitude. The span max(dE_a)/min(dE_a) therefore separates the two without needing the
        threshold, the box or anything else from outside. Calibrated on one H(1s) np n = 2..16 set at four boxes:

            span 54.6 (64 a.u.)    alpha^(E1,0) 1.00000 of reference,   C_3 0.99995 of exact
            span  6.7 (200 a.u.)                0.99992                     0.99908
            span  2.2 (400 a.u.)                0.98888                     0.96403
            span  1.4 (775 a.u.)                0.86054                     0.77026

        and on nd n = 3..16, where a span of 1.15 goes with alpha^(E2,0) at 0.443 of reference. A span below 3 is
        therefore warned about, and the warning names C_3 because that is the quantity which suffers most.

        WHAT THIS GUARD DOES NOT CATCH, said here so that it is not relied upon for it: a box that is far too LARGE
        inflates a high multipole rather than starving it, and the span is then huge while the number is wrong -- E3
        on a 775 a.u. box gives a span of 10.4 and 1.275 of its reference. That failure belongs to
        `MultipolePolarizibility.checkPerturberSet` and `Bsplines.checkOrbitalConsistency`, which are run when the
        polarizability itself is computed.

        ONE CANDIDATE DETECTOR WAS MEASURED AND REJECTED, and its numbers are kept so that it is not rebuilt: the
        SHARE of sum_a A_a carried by the highest-lying pole runs 0.0000 on the perfect E1 set and 0.0306 on the
        23 %-short one -- discriminating in the wrong direction and by far too little -- and it is 0.0000 for E3 in
        every box, good and bad alike.
"""
function checkSaturation(partner::DispersionCoefficient.Partner)
    for  terms in partner.terms
        isempty(terms.energies)   &&   continue
        dEmin = minimum(terms.energies);    dEmax = maximum(terms.energies)
        dEmin <= 0.0              &&   continue      ## an excited partner;  checkGroundLevel refuses it anyway
        span  = dEmax / dEmin
        if  span < 3.0
            println(">>> DispersionCoefficient: '$(partner.name)', $(terms.multipole): the $(length(terms.energies)) " *
                    "perturbers span only " * @sprintf("%.2f", span) * " in excitation energy " *
                    "(" * @sprintf("%.4f", dEmin) * " to " * @sprintf("%.4f", dEmax) * " a.u.), so the series has " *
                    "ACCUMULATED AT A THRESHOLD and represents no continuum.")
            println(">>>   Every quantity here is then a LOWER bound, and C_3 the weakest of them: measured on " *
                    "H(1s) np n = 2..16 at a span of 1.4, alpha^(E1,0) came out 0.861 of its reference and C_3 " *
                    "0.770 of the exact value.  REDUCE the radial box until the upper members run into the " *
                    "continuum;  a box matched to the highest orbital is what causes this.")
        end
    end

    return( nothing )
end


"""
`DispersionCoefficient.checkWall(settings::DispersionCoefficient.Settings)`
    ... refuses a wall that is not a dielectric BEFORE any perturber set is built; `nothing` is returned. Neither
        eps(0) < 1, which has no reflection factor on the imaginary axis, nor omega_0 <= 0, which has no pole, is an
        extreme case to be clamped; both are typing errors. WHAT THIS BUYS IS THE TIMING, since the wall enters only
        at the END of `DispersionCoefficient.computeOutcome`: a mistyped eps(0) was discovered after the multipole
        sums of the first selected level had been assembled, and the remaining levels never ran at all. The algebra
        is `DispersionCoefficient.wallReflectionPole` and costs nothing, so it is done once here as well. A perfect
        conductor reflects fully at every frequency and so carries no resonance to check.
"""
function checkWall(settings::DispersionCoefficient.Settings)
    settings.calcAtomWall              ||   return( nothing )
    isinf(settings.wallDielectric)     &&   return( nothing )
    if  settings.wallDielectric < 1.0
        error("\nDispersionCoefficient: a wall needs a static dielectric constant eps(0) >= 1, but eps(0) = " *
              "$(settings.wallDielectric) was given;  Inf asks for a perfect conductor.  Below 1 the reflection " *
              "factor (eps-1)/(eps+1) comes out NEGATIVE, which no dielectric does, so this is refused here " *
              "rather than silently carried into C_3.")
    end
    isinf(settings.wallResonance)      ||   DispersionCoefficient.wallReflectionPole(settings.wallDielectric,
                                                                                    settings.wallResonance)

    return( nothing )
end


"""
`DispersionCoefficient.computeAtomWallCoefficient(partner::DispersionCoefficient.Partner, wallDielectric::Float64)`
    ... computes the non-retarded atom-wall coefficient C_3 of `partner` against a wall that responds
        INSTANTANEOUSLY on the atom's own frequency scale, i.e. in the limit omega_0 >> dE_a; a `C3::Float64` in
        atomic units is returned. It is `DispersionCoefficient.computeAtomWallCoefficient(partner, wallDielectric,
        Inf)` and is kept as its own method because that limit is the only wall most callers want.
"""
function computeAtomWallCoefficient(partner::DispersionCoefficient.Partner, wallDielectric::Float64)
    return( DispersionCoefficient.computeAtomWallCoefficient(partner, wallDielectric, Inf) )
end


"""
`DispersionCoefficient.computeAtomWallCoefficient(partner::DispersionCoefficient.Partner, wallDielectric::Float64,
                                                   wallResonance::Float64)`
    ... computes the non-retarded atom-wall coefficient C_3 of `partner` in U(z) = - C_3/z^3, for a wall of STATIC
        dielectric constant `wallDielectric` (Inf for a perfect conductor) whose one Lorentz oscillator sits at
        `wallResonance` (Inf for the instantaneous-response limit); a `C3::Float64` in atomic units is returned, and
        0. if the partner carries no E1 response.

        The Lennard-Jones form C_3 = (1/4 pi) int_0^inf alpha^(E1,0)(i w) dw collapses on a discrete spectrum to
        C_3 = (1/8) sum_a A_a, since int_0^inf dE/(dE^2 + w^2) dw = pi/2 for every pole alike. The wall enters as the
        reflection factor r(i w), which for one oscillator is itself a single pole R/(Omega^2 + w^2) -- see
        `DispersionCoefficient.wallReflectionPole` -- so the integral stays closed and

        C_3 = (R/8) sum_a A_a / [Omega (dE_a + Omega)],   collapsing to (eps-1)/(eps+1) times (1/8) sum_a A_a as
        Omega >> dE_a.

        SO THE CONSTANT-eps ANSWER IS A LIMIT AND NOT A SEPARATE MODEL, and its error is one factor per pole,
        exactly 1 + dE_a/Omega. For a silica-like wall (eps(0) = 3.9, omega_0 = 10 eV, Omega = 0.575 a.u.) against
        H(1s) that is 65 % on the 2p pole alone and 81 % measured over the np n = 2..16 set, so a real dielectric is
        NOT well described by its static constant alone. A perfect conductor reflects fully at every frequency and
        ignores `wallResonance` altogether.

        The perfectly conducting closed form is exactly checkable: with A_a = (2/3) |<a|d|0>|^2 the dipole sum rule
        gives C_3 = <d^2>/12, which is what an electrostatic image dipole yields, so C_3(H 1s) = 1/4 exactly.

        BEWARE THAT A POOR PERTURBER SET DAMAGES THIS BEFORE ANYTHING ELSE, because sum_a A_a carries no energy
        denominator while alpha(0) and every C_n suppress the high poles: measured on one H(1s) np n = 2..16 set,
        moving the radial box from 64 to 775 a.u. costs alpha^(E1,0) 14 % and C_3 23 %. See
        `DispersionCoefficient.checkSaturation`.
"""
function computeAtomWallCoefficient(partner::DispersionCoefficient.Partner, wallDielectric::Float64,
                                    wallResonance::Float64)
    terms = DispersionCoefficient.selectTerms(partner, 1)
    if  isnothing(terms)
        println(">>> DispersionCoefficient: '$(partner.name)' carries no E1 response, so no atom-wall C_3.")
        return( 0.0 )
    end
    ## A perfect conductor has r = 1 at EVERY frequency, so no resonance can enter; the same closed form serves
    ## the instantaneous-response limit, where r is the constant (eps-1)/(eps+1).
    if  isinf(wallDielectric)   ||   isinf(wallResonance)
        reflection = isinf(wallDielectric) ? 1.0 : (wallDielectric - 1.0) / (wallDielectric + 1.0)
        return( 0.125 * sum(terms.strengths) * reflection )
    end
    strength, Omega = DispersionCoefficient.wallReflectionPole(wallDielectric, wallResonance)
    wa = 0.0
    for  a = 1:length(terms.strengths)
        wa = wa + terms.strengths[a] / ( Omega * (terms.energies[a] + Omega) )
    end

    return( 0.125 * strength * wa )
end


"""
`DispersionCoefficient.computeCoefficient(n::Int64, partnerA::DispersionCoefficient.Partner,
                                           partnerB::DispersionCoefficient.Partner)`
    ... computes the atom-atom dispersion coefficient C_n of the pair (partnerA, partnerB) in
        V(R) = - ... - C_n/R^n - ..., by summing the Casimir-Polder weight times the polarizability overlap integral
        over every ORDERED multipole pair (l1,l2) with l1 + l2 = n/2 - 1; a
        `coefficient::DispersionCoefficient.Coefficient` in atomic units is returned, carrying its own (l1,l2) pieces.
        Both partners are first checked by `DispersionCoefficient.checkGroundLevel`, which raises for an excited level.
        A pair whose multipole either partner does not carry is SKIPPED with a warning and marks the result
        incomplete, so that a coefficient missing a term is never mistaken for a converged one.
"""
function computeCoefficient(n::Int64, partnerA::DispersionCoefficient.Partner,
                            partnerB::DispersionCoefficient.Partner)
    if  isodd(n)  ||  n < 6
        error("DispersionCoefficient: C_$n is not defined here -- n must be even and at least 6.")
    end
    DispersionCoefficient.checkGroundLevel(partnerA);    DispersionCoefficient.checkGroundLevel(partnerB)

    lSum     = div(n, 2) - 1
    value    = 0.0;     pieces = Tuple{Int64,Int64,Float64}[];     complete = true
    for  l1 = 1:lSum-1
        l2      = lSum - l1
        termsA  = DispersionCoefficient.selectTerms(partnerA, l1)
        termsB  = DispersionCoefficient.selectTerms(partnerB, l2)
        if  isnothing(termsA)  ||  isnothing(termsB)
            complete = false
            println(">>> DispersionCoefficient: C_$n needs the (E$l1, E$l2) pair, and " *
                    "'$(isnothing(termsA) ? partnerA.name : partnerB.name)' carries no " *
                    "E$(isnothing(termsA) ? l1 : l2) response;  this term is MISSING and C_$n is a lower bound.")
            continue
        end
        wa      = DispersionCoefficient.casimirPolderPrefactor(l1, l2) *
                  DispersionCoefficient.overlapIntegral(termsA, termsB)
        value   = value + wa
        push!( pieces, (l1, l2, wa) )
    end

    return( DispersionCoefficient.Coefficient(n, value, pieces, complete) )
end


"""
`DispersionCoefficient.computeOutcome(partnerA::DispersionCoefficient.Partner, partnerB::DispersionCoefficient.Partner,
                                       settings::DispersionCoefficient.Settings)`
    ... computes every dispersion coefficient asked for in `settings` for the pair (partnerA, partnerB), i.e. the
        atom-atom C_n of `settings.orders` and, if requested, the atom-wall C_3 of partnerA; an
        `outcome::DispersionCoefficient.Outcome` is returned. The saturation of both partners is reported first, since
        a coefficient is only as good as the perturber set under it.
"""
function computeOutcome(partnerA::DispersionCoefficient.Partner, partnerB::DispersionCoefficient.Partner,
                        settings::DispersionCoefficient.Settings)
    DispersionCoefficient.checkSaturation(partnerA)
    partnerB === partnerA   ||   DispersionCoefficient.checkSaturation(partnerB)

    coefficients = DispersionCoefficient.Coefficient[]
    for  n in sort(settings.orders)
        push!( coefficients, DispersionCoefficient.computeCoefficient(n, partnerA, partnerB) )
    end
    c3 = settings.calcAtomWall ?
         DispersionCoefficient.computeAtomWallCoefficient(partnerA, settings.wallDielectric,
                                                          settings.wallResonance) : 0.0

    return( DispersionCoefficient.Outcome(partnerA, partnerB, coefficients, c3, settings.wallDielectric,
                                          settings.wallResonance) )
end


"""
`DispersionCoefficient.computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                                        settings::DispersionCoefficient.Settings; output=true)`
    ... computes the dispersion coefficients for the selected levels of the given multiplet, which is the entry point
        `perform(::Atomic.Computation)` uses, so that the module is reached exactly as every other property module is;
        an `outcomes::Array{DispersionCoefficient.Outcome,1}` is returned. The computed atom is partner A and its
        response is built from `settings.gMultiplet`; partner B is `settings.partnerB` where the caller supplied one
        and partner A itself otherwise, which is the homonuclear case every published C_n of a single species means.
        Each selected level passes `DispersionCoefficient.checkGroundLevel` FIRST, so a selection that admits an
        excited level raises here rather than returning a number the Casimir-Polder form cannot support.
"""
function computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
            settings::DispersionCoefficient.Settings; output=true)
    println("")
    printstyled("DispersionCoefficient.computeOutcomes(): The computation of dispersion coefficients starts now ... \n",
                color=:light_green)
    printstyled("------------------------------------------------------------------------------------------------- \n",
                color=:light_green)
    #
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    DispersionCoefficient.checkWall(settings)
    outcomes = DispersionCoefficient.Outcome[]
    for  level  in  multiplet.levels
        Basics.selectLevel(level, settings.levelSelection)   ||   continue
        partnerA = DispersionCoefficient.Partner( "level $(level.index), $(LevelSymmetry(level.J, level.parity))",
                                                  level, settings.gMultiplet, grid, settings.multipoles )
        DispersionCoefficient.checkGroundLevel(partnerA)
        partnerB = isempty(settings.partnerB.terms) ? partnerA : settings.partnerB
        if  settings.printBefore
            println("");   println(partnerA);   partnerB === partnerA   ||   println(partnerB)
        end
        outcome = DispersionCoefficient.computeOutcome(partnerA, partnerB, settings)
        push!( outcomes, outcome )
        DispersionCoefficient.displayResults(stdout, outcome)
        if  printSummary    DispersionCoefficient.displayResults(iostream, outcome)   end
    end

    if    output    return( outcomes )
    else            return( nothing )
    end
end


"""
`DispersionCoefficient.displayResults(stream::IO, outcome::DispersionCoefficient.Outcome)`
    ... prints the dispersion coefficients of `outcome` in a neat table, each with the (l1,l2) multipole pairs it is
        built from and the static polarizabilities that carry it, so that a reader can see at once which response
        would have to be improved; nothing is returned.
"""
function displayResults(stream::IO, outcome::DispersionCoefficient.Outcome)
    nx = 103
    println(stream, " ")
    println(stream, "  Long-range dispersion coefficients of  $(outcome.partnerA.name) -- $(outcome.partnerB.name)  " *
                    "(atomic units):")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  "
    sa = sa * TableStrings.center(10, "Coefficient"; na=4)
    sa = sa * TableStrings.center(20, "Value [a.u.]";  na=4)
    sa = sa * TableStrings.center( 8, "complete";  na=4)
    sa = sa * TableStrings.center(44, "(l1,l2) pieces [a.u.]";  na=2)
    println(stream, sa);    println(stream, "  ", TableStrings.hLine(nx))
    #
    for  coefficient in outcome.coefficients
        sa = "  " * TableStrings.center(10, "C_$(coefficient.n)"; na=4)
        sa = sa * TableStrings.flushright(20, @sprintf("%.8e", coefficient.value) )                       * "    "
        sa = sa * TableStrings.center( 8, coefficient.complete ? "yes" : "NO"; na=4)
        sb = ""
        for  (l1, l2, wa) in coefficient.pieces
            sb = sb * "(E$l1,E$l2) " * @sprintf("%.6e", wa) * "   "
        end
        println(stream, sa * sb)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    if  outcome.atomWallC3 != 0.0
        wall = isinf(outcome.wallDielectric) ? "a perfectly conducting wall" :
               isinf(outcome.wallResonance)    ? "a wall of constant eps = " * @sprintf("%.4f", outcome.wallDielectric) :
               "a wall with eps(0) = " * @sprintf("%.4f", outcome.wallDielectric) * " at omega_0 = " *
               @sprintf("%.4f", outcome.wallResonance) * " a.u."
        println(stream, " ")
        println(stream, "  Atom-wall coefficient of $(outcome.partnerA.name) against " * wall *
                        ",  U(z) = - C_3/z^3:")
        println(stream, "      C_3 = " * @sprintf("%.8e", outcome.atomWallC3) * " a.u.")
        println(stream, "  ")
        if  !isinf(outcome.wallDielectric)  &&  !isinf(outcome.wallResonance)
            om = DispersionCoefficient.wallReflectionPole(outcome.wallDielectric, outcome.wallResonance)[2]
            println(stream, "  ")
            println(stream, "  >>> This wall is frequency dependent, Omega = " * @sprintf("%.4f", om) * " a.u., so " *
                            "C_3 is SMALLER than its constant-eps value")
            println(stream, "  >>> by the factor 1 + dE_a/Omega on every perturber.  The constant-eps number is the " *
                            "Omega >> dE limit, not a model of its own.")
        end
        println(stream, "  >>> C_3 = (1/8) sum_a A_a carries NO energy denominator, so a poor perturber set damages")
        println(stream, "  >>> it before anything else:  measured on one H(1s) np n = 2..16 set, moving the radial")
        println(stream, "  >>> box from 64 to 775 a.u. costs alpha^(E1,0) 14 % and C_3 23 %, and alpha(0) does not")
        println(stream, "  >>> show it.  Read the saturation warnings above before quoting C_3.")
    end
    #
    println(stream, " ")
    println(stream, "  The static polarizabilities these coefficients rest on:")
    for  partner in (outcome.partnerA, outcome.partnerB)
        for  terms in partner.terms
            println(stream, "      $(partner.name)   alpha^($(terms.multipole),0)(0) = " *
                            @sprintf("%.8e", DispersionCoefficient.staticPolarizability(terms)) *
                            " a.u.   from $(length(terms.strengths)) poles, lowest dE = " *
                            @sprintf("%.6f", isempty(terms.energies) ? 0.0 : minimum(terms.energies)) * " a.u.")
        end
        partner === outcome.partnerB   &&   break
    end
    println(stream, " ")

    return( nothing )
end


"""
`DispersionCoefficient.multipoleTerms(multipole::EmMultipole, level::Level, gMultiplet::Multiplet, grid::Radial.Grid)`
    ... builds the single-pole representation of `level`'s 2^L-pole response from the perturbers of `gMultiplet`, by
        asking `MultipolePolarizibility.computeMultipolePolarizabilities` for its per-perturber contributions and
        keeping, for each, the rank-0 strength A_a = c(L,0,J,J_a) |<J_a || O^(E L) || J>|^2 and the excitation energy
        dE_a = W_a - W_level; a `terms::DispersionCoefficient.MultipoleTerms` is returned. Poles of vanishing strength
        -- a perturber the multipole cannot reach, or one outside the rank-0 triangle -- are dropped, so the returned
        lists hold only terms that contribute.

        NO NEW STRUCTURE CALCULATION HAPPENS HERE, and that is the point of taking this route rather than evaluating
        polarizabilities at many frequencies: the stored reduced matrix element and energy denominator are all a
        dispersion coefficient ever needs, so one pass over the perturbers serves every C_n and the atom-wall C_3 at
        once.
"""
function multipoleTerms(multipole::EmMultipole, level::Level, gMultiplet::Multiplet, grid::Radial.Grid)
    alphas, contributions = MultipolePolarizibility.computeMultipolePolarizabilities(multipole, level, gMultiplet, grid)
    strengths = Float64[];      energies = Float64[]
    for  contribution in contributions
        coeff = MultipolePolarizibility.rankCoefficient(multipole.L, 0, level.J, contribution.nuLevel.J)
        wa    = coeff * abs2(contribution.reducedMe)
        wa == 0.0   &&   continue
        push!( strengths, wa );     push!( energies, -contribution.deltaE )
    end

    return( DispersionCoefficient.MultipoleTerms(multipole, strengths, energies) )
end


"""
`DispersionCoefficient.overlapIntegral(termsA::DispersionCoefficient.MultipoleTerms,
                                        termsB::DispersionCoefficient.MultipoleTerms)`
    ... computes the imaginary-frequency overlap integral int_0^inf alpha_A(i w) alpha_B(i w) dw of two multipole
        responses in CLOSED FORM; an `integral::Float64` in atomic units is returned.

        No quadrature is involved and none is wanted: on a discrete spectrum each polarizability is a sum of single
        poles, and int_0^inf [a/(a^2+w^2)] [b/(b^2+w^2)] dw = pi/(2(a+b)) exactly, so the integral is the double sum
        (pi/2) sum_ab A_a A_b/(dE_a + dE_b) -- Dalgarno's oscillator-strength form. Reading the weight 1/(dE_a + dE_b)
        is also how one sees which perturbers a dispersion coefficient actually rests on: the low-lying ones, though
        less exclusively so than for alpha(0), whose weight is 1/dE_a.
"""
function overlapIntegral(termsA::DispersionCoefficient.MultipoleTerms, termsB::DispersionCoefficient.MultipoleTerms)
    integral = 0.0
    for  a = 1:length(termsA.strengths),  b = 1:length(termsB.strengths)
        integral = integral + termsA.strengths[a] * termsB.strengths[b] /
                              (termsA.energies[a] + termsB.energies[b])
    end

    return( 0.5 * pi * integral )
end


"""
`DispersionCoefficient.polarizability(terms::DispersionCoefficient.MultipoleTerms, omega::Float64)`
    ... returns the scalar multipole polarizability alpha^(E L,0)(i omega) of `terms` at the IMAGINARY frequency
        i*omega, i.e. sum_a A_a dE_a/(dE_a^2 + omega^2); an `alpha::Float64` in atomic units is returned. It is
        provided because it is what a user plots and what makes this module's closed-form integrals checkable by
        quadrature, not because any coefficient here needs it: the imaginary axis has no poles, so this function is
        finite, positive and monotone decreasing for every omega.
"""
function polarizability(terms::DispersionCoefficient.MultipoleTerms, omega::Float64)
    alpha = 0.0
    for  a = 1:length(terms.strengths)
        alpha = alpha + terms.strengths[a] * MultipolePolarizibility.frequencyFactor(terms.energies[a], omega, true)
    end

    return( alpha )
end


"""
`DispersionCoefficient.selectTerms(partner::DispersionCoefficient.Partner, L::Int64)`
    ... picks out `partner`'s electric 2^L-pole response; a `terms::DispersionCoefficient.MultipoleTerms` is returned,
        or `nothing` if the partner carries no response for that multipole, which the callers report as an incomplete
        coefficient rather than treat as zero.
"""
function selectTerms(partner::DispersionCoefficient.Partner, L::Int64)
    for  terms in partner.terms
        if  terms.multipole.L == L  &&  terms.multipole.electric    return( terms )    end
    end

    return( nothing )
end


"""
`DispersionCoefficient.staticPolarizability(terms::DispersionCoefficient.MultipoleTerms)`
    ... returns the static scalar multipole polarizability alpha^(E L,0)(0) = sum_a A_a/dE_a of `terms`; an
        `alpha::Float64` in atomic units is returned. It is the omega = 0 member of
        `DispersionCoefficient.polarizability` and is given its own method because it is the number a reader compares
        against the literature to judge whether a dispersion coefficient rests on anything.
"""
function staticPolarizability(terms::DispersionCoefficient.MultipoleTerms)
    alpha = 0.0
    for  a = 1:length(terms.strengths)
        alpha = alpha + terms.strengths[a] / terms.energies[a]
    end

    return( alpha )
end


"""
`DispersionCoefficient.wallReflectionPole(wallDielectric::Float64, wallResonance::Float64)`
    ... writes the imaginary-frequency reflection factor of a wall with one undamped Lorentz oscillator as a SINGLE
        POLE, which is what keeps the atom-wall integral closed; a `(strength, energy)::Tuple{Float64,Float64}` in
        atomic units is returned, meaning r(i w) = strength/(energy^2 + w^2).

        With eps(i w) = 1 + omega_p^2/(omega_0^2 + w^2) the factor (eps-1)/(eps+1) is omega_p^2 / (omega_p^2 +
        2 omega_0^2 + 2 w^2), hence strength = omega_p^2/2 and energy^2 = omega_0^2 + omega_p^2/2. Expressed through
        the STATIC constant eps(0) = 1 + omega_p^2/omega_0^2, which is what a wall is tabulated by,

            strength = (eps(0) - 1) omega_0^2 / 2,     energy = omega_0 sqrt( (eps(0) + 1) / 2 )

        and the three limits are worth knowing because each is a physical wall: energy >> dE returns the constant
        (eps(0)-1)/(eps(0)+1), eps(0) -> Inf returns 1 at every frequency (a perfect conductor), and omega_0 -> 0
        returns 0 -- a wall whose charges are too slow to follow the atom does not reflect at the frequencies that
        carry C_3. The function RAISES for eps(0) < 1 or omega_0 <= 0, neither of which is a dielectric.
"""
function wallReflectionPole(wallDielectric::Float64, wallResonance::Float64)
    if  wallDielectric < 1.0   ||   wallResonance <= 0.0
        error("A wall needs a static dielectric constant eps(0) >= 1 and a resonance omega_0 > 0, but " *
              "eps(0) = $wallDielectric and omega_0 = $wallResonance were given;  Inf for either means " *
              "a perfect conductor or the instantaneous-response limit and is handled without this function.")
    end
    strength = (wallDielectric - 1.0) * wallResonance^2 / 2.0
    energy   = wallResonance * sqrt( (wallDielectric + 1.0) / 2.0 )

    return( (strength, energy) )
end

end # module
