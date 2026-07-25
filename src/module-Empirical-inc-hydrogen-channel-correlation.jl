
#################################################################################################################################
### Molecular-symmetry channel correlation for atom/ion -- hydrogen collisions ##################################################
##
##  Automates Steps 1-2 of A. K. Belyaev & S. A. Yakovleva, A&A 608, A33 (2017) [examples/papers/b17.aa-belyav-inelastic.pdf]:
##  given the ground term of the ionic entrance state A^(Z+1)+ + H^- and a list of candidate atomic/ionic levels A^Z+(j), decide
##  which levels j correlate (by molecular-term symmetry) to the ionic entrance state, and what statistical weight p_stat_j each
##  gets -- i.e. it turns a plain list of atomic levels into the Array{Empirical.InelasticHChannel,1} already consumed, unchanged,
##  by Empirical.neutralizationReducedRate / Empirical.inelasticHCollisionRateMatrix (module-Empirical-inc-inelastic-h-collisions.jl).
##
##  Physics in one paragraph: only atomic/ionic states A^Z+(j) whose combination with a hydrogen partner produces a molecular
##  term of the SAME symmetry (total spin S, orbital projection Lambda) as the ionic entrance state A^(Z+1)+ + H^- can couple to
##  it via the Landau-Zener ionic-covalent crossing (Belyaev & Yakovleva 2017, Sec. 2.1). Since one partner is always hydrogen --
##  H^-(1s^2, ^1S, L=0, S=0) on the ionic side, neutral H(1s, ^2S, L=0, S=1/2) on every covalent side -- the standard atom-to-
##  molecule (Wigner-Witmer) term correlation collapses to simple closed-form combinatorics: combining an atomic (L,S) term with
##  a partner of L=0 spreads the atomic term's (2L+1)(2S+1) magnetic sublevels over molecular Lambda = 0,1,...,L (multiplicity 1
##  for Lambda=0, 2 for each Lambda>0) and, independently, over molecular S_mol = S +/- S_partner (multiplicity 2*S_mol+1 each).
##
##  statisticalWeight below is this module's own combinatorial reconstruction of p_stat_j -- the paper states only that "the
##  statistical probabilities p_stat_j can be readily calculated" (Sec. 2.1) without giving an explicit formula. It was derived
##  from first principles (the paragraph above) and then checked, level by level, against every one of the 19 tabulated p_stat_j
##  values in the paper's own Table 1 (Ba+ + H, Z=1): it reproduces all of them exactly (1/4 for every ^2S channel, 1/12 for every
##  ^2P, 1/20 for ^2D, 1/28 for ^2F, 1/36 for ^2G) -- see example-Nl.jl for the side-by-side check.
##
##  Known scope limitations (not silently glossed over):
##  - Reflection parity (Sigma+/Sigma-) for Lambda=0 terms is NOT resolved separately; the formula below treats a shared
##    (Lambda=0, S_mol) pair as always matching. This is exact for the Ba+H case (both partners' relevant terms are always
##    Sigma+ there) but has not been checked against a case where it could matter (e.g. a candidate level whose Lambda=0
##    component might carry Sigma- character); flagged here for future verification, not assumed away.
##  - This module does NOT decide *which* levels of A^Z+ belong in the candidate list in the first place -- selecting states
##    that correspond to genuine one-electron transitions from the same ionic core/parent as A^(Z+1)+ (paper's Step 2, second
##    clause) is a configuration/parentage judgment this module has no information to make; the caller must supply a candidate
##    list already restricted to single-active-electron states sharing that parent (true automatically for simple one-valence-
##    electron ions like the alkalis, alkaline earths once ionized, and Ba+ -- not necessarily for open-shell/transition-metal
##    ions with several possible parent cores).


"""
`struct  Empirical.AtomicLevel`
    ... a single atomic/ionic term (or, optionally, one fine-structure level within a term) of the covalent species A^Z+(j),
        used as the input candidate for Empirical.generateInelasticHChannels. L and S are the standard LS-coupling quantum
        numbers of the term (e.g. from NIST, or from the dominant LS parentage of a JAC structure level); J, if given, selects
        one fine-structure level within that term for a finer-grained (state-resolved) treatment instead of a term-averaged one.

    + name  ::String                 ... a short label, e.g. "Ba+(6d 2D)", for display purposes only.
    + L     ::Int64                  ... orbital angular momentum quantum number of the term.
    + S     ::Float64                ... spin quantum number of the term (0, 1/2, 1, 3/2, ...).
    + J     ::Union{Nothing,Float64} ... total angular momentum of one fine-structure level within the term, or `nothing` for
                                          a term-averaged level (see Empirical.statisticalWeight for how J is used).
    + E     ::Float64                ... electronic bound energy [a.u.], negative, referenced to the ionization limit of A^Z+
                                          (same convention as Empirical.InelasticHChannel.E).
"""
struct  AtomicLevel
    name    ::String
    L       ::Int64
    S       ::Float64
    J       ::Union{Nothing,Float64}
    E       ::Float64
end


# `Base.show(io::IO, level::AtomicLevel)`  ... provides a String notation for the variable level::AtomicLevel.
function Base.show(io::IO, level::AtomicLevel)
    Jstring = isnothing(level.J) ? "term-averaged" : "J = $(level.J)"
    print(io, "$(level.name):  L = $(level.L),  S = $(level.S),  $Jstring,  E = $(level.E) [a.u.]")
end


"""
`Empirical.molecularOrbitalFraction(Lj::Int64, Lion::Int64)`
    ... to compute the fraction of channel j's own (2*Lj+1) orbital magnetic sublevels that fall on a molecular Lambda value
        also reachable by the ionic entrance term (Lion combined with a partner of L=0, i.e. Lambda_ion = 0,1,...,Lion). Since
        both Lj and Lion are combined with an L=0 partner (H or H^-), their own Lambda ranges are 0,1,...,Lj and 0,1,...,Lion
        respectively, so the shared range is 0,1,...,min(Lj,Lion) -- always including Lambda=0. A value::Float64 in (0,1] is
        returned.
"""
function molecularOrbitalFraction(Lj::Int64, Lion::Int64)
    sharedMultiplicity = 1.0 + 2.0 * min(Lj, Lion)
    return( sharedMultiplicity / (2*Lj + 1) )
end


"""
`Empirical.molecularSpinFraction(Sj::Float64, Sion::Float64)`
    ... to compute the fraction of channel j's own 2*(2*Sj+1) spin magnetic sublevels (from combining the atomic spin Sj with
        neutral hydrogen's S=1/2) that fall on the molecular spin S_mol = Sion required by the ionic entrance term (combining
        the ionic term's own spin Sion with H^-'s S=0, so the ionic side has the single, unsplit value S_mol = Sion). Returns
        0.0 if Sion is not among {Sj+1/2, Sj-1/2} (no one-electron-transfer path conserves spin between the two terms). A
        value::Float64 in [0,1] is returned.
"""
function molecularSpinFraction(Sj::Float64, Sion::Float64)
    total = 2.0 * (2*Sj + 1)
    if      Sion == Sj + 0.5                  matched = 2*Sion + 1
    elseif  Sj > 0.0   &&   Sion == Sj - 0.5   matched = 2*Sion + 1
    else                                       matched = 0.0
    end
    return( matched / total )
end


"""
`Empirical.statisticalWeight(Lj::Int64, Sj::Float64, Lion::Int64, Sion::Float64)`
    ... to compute the statistical probability p_stat_j (Belyaev & Yakovleva 2017, Eqs. 1-2) that channel j = A^Z+(j) + H
        correlates, by molecular-term symmetry, to the ionic entrance state A^(Z+1)+ + H^- of term (Lion,Sion), as the product
        of Empirical.molecularOrbitalFraction and Empirical.molecularSpinFraction. Returns 0.0 (no correlation, channel j is
        not accessible from this ionic entrance) if the spin fraction is zero; see the module note above for the derivation
        and its verification against Table 1 of the paper. A value::Float64 in [0,1) is returned.
"""
function statisticalWeight(Lj::Int64, Sj::Float64, Lion::Int64, Sion::Float64)
    return( Empirical.molecularOrbitalFraction(Lj, Lion) * Empirical.molecularSpinFraction(Sj, Sion) )
end


"""
`Empirical.statisticalWeight(level::Empirical.AtomicLevel, Lion::Int64, Sion::Float64)`
    ... convenience method that takes the term-level statistical weight from Empirical.statisticalWeight(level.L, level.S,
        Lion, Sion) and, if level.J is given (a specific fine-structure level rather than a term average), redistributes it
        over the term's J-multiplet in proportion to (2*level.J+1)/[(2*level.L+1)(2*level.S+1)] -- exact in the degenerate
        (negligible fine-structure splitting) limit, since summing (2J+1) over a full J-multiplet reproduces (2L+1)(2S+1)
        exactly, and a strict refinement otherwise since each J-level then keeps its own true bound energy level.E for the
        crossing-radius calculation instead of a single term-averaged energy. A value::Float64 in [0,1) is returned.
"""
function statisticalWeight(level::Empirical.AtomicLevel, Lion::Int64, Sion::Float64)
    pTerm = Empirical.statisticalWeight(level.L, level.S, Lion, Sion)
    if  isnothing(level.J)   return( pTerm )   end
    return( pTerm * (2*level.J + 1) / ((2*level.L + 1) * (2*level.S + 1)) )
end


"""
`Empirical.generateInelasticHChannels(levels::Array{Empirical.AtomicLevel,1}, Lion::Int64, Sion::Float64;
                                      printout::Bool=false)`
    ... to build the Array{Empirical.InelasticHChannel,1} needed by Empirical.inelasticHCollisionRateMatrix (and, level by
        level, by Empirical.neutralizationReducedRate/deExcitationReducedRate) directly from a candidate list of atomic/ionic
        term (or fine-structure) levels of A^Z+ plus the (Lion,Sion) term of the ionic entrance state A^(Z+1)+. Every level
        is kept only if Empirical.statisticalWeight(level, Lion, Sion) > 0.0 (a genuine molecular-symmetry correlation);
        levels that do not correlate are dropped silently unless printout=true, in which case they are listed explicitly
        together with the kept channels and their statistical weights. A channels::Array{Empirical.InelasticHChannel,1} is
        returned.
        Quantity: this function performs no rate-coefficient calculation itself -- it only turns level data (L,S,J,E) into
            the (E,p_stat) pairs that Empirical.InelasticHChannel carries; hand the result to
            Empirical.inelasticHCollisionRateMatrix as usual.
"""
function generateInelasticHChannels(levels::Array{Empirical.AtomicLevel,1}, Lion::Int64, Sion::Float64;
                                    printout::Bool=false)
    channels  = Empirical.InelasticHChannel[]
    discarded = Empirical.AtomicLevel[]
    for  level  in  levels
        p = Empirical.statisticalWeight(level, Lion, Sion)
        if  p > 0.0   push!(channels,  Empirical.InelasticHChannel(level.name, level.E, p))
        else          push!(discarded, level)
        end
    end

    if  printout
        println("\n* Generate inelastic-H-collision channels by molecular-symmetry correlation, for the ionic entrance " *
                "term A^(Z+1)+ + H^- with (Lion,Sion) = ($Lion, $Sion): " *
                "\n    + A level A^Z+(j) is KEPT as a channel only if combining it with neutral H produces a molecular " *
                "state of the SAME symmetry as the ionic entrance -- i.e. only if p_stat_j > 0. Without a shared " *
                "symmetry, the two potential curves may still cross, but nothing can couple them (Belyaev & Yakovleva " *
                "2017, Steps 1-2). " *
                "\n    + p_stat_j = the FRACTION of level j's own magnetic sublevels that land on that shared symmetry. " *
                "It is NOT a factor in the forward (neutralization) rate coefficient -- that only needs p_stat_ionic " *
                "(passed to Empirical.inelasticHCollisionRateMatrix, normally 1 since the ionic entrance term is " *
                "unique). p_stat_j instead governs the REVERSE (ion-pair-formation/excitation) rate, via " *
                "Empirical.detailedBalanceRate. " *
                "\n    + The crossing radius Rx for a kept channel (where the ionic and covalent potential curves " *
                "meet) is NOT printed here -- it is a separate diagnostic, available via " *
                "Empirical.nonadiabaticRadius(level.E, EHminus, Z) for anyone wanting to sanity-check that deeper-" *
                "bound levels get a smaller Rx. ")
        println("    Correlated channels (kept), with their statistical weight p_stat:")
        for  ch  in  channels   println("      $ch")   end
        if  !isempty(discarded)
            println("    Discarded (no shared molecular symmetry with the ionic entrance term -- p_stat = 0, no " *
                    "Landau-Zener coupling possible for these levels from this ionic entrance):")
            for  level  in  discarded   println("      $level")   end
        end
        println()
    end

    return( channels )
end
