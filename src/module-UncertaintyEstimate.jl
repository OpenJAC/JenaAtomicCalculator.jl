using Printf

"""
`module  JAC.UncertaintyEstimate`
    ... a submodel of JAC that collects a handful of DELIBERATELY, PHYSICALLY named alternative computations
        of the SAME quantity -- e.g. the same transition rate computed once with the DFS field and once with
        the AL field, or once from a reference-only correlation space and once with one more correlation
        layer added -- and reports the spread across them.

        WHAT THIS IS NOT, STATED PLAINLY BEFORE ANYTHING ELSE. This is NOT a calibrated statistical
        uncertainty in the metrology sense (no claim of a confidence level, no sampling of a probability
        distribution), and it is NOT the Kramida-style Monte-Carlo-over-many-random-perturbations idea in
        full: that would need hundreds of full SCF+CI re-runs per quantity, which is not affordable for a
        code built the way JAC is (unlike the semi-empirical, cheap-to-rerun codes the method was designed
        for). What IS reported is the spread across a SMALL set of variants the caller chose because each one
        is a real, physically understood, already-available JAC alternative (a different mean field, a
        different correlation extent, ...) -- and because only a few of the many possible sources of error
        are ever sampled this way, the honest reading of the result is a LOWER BOUND on the true uncertainty,
        not a full characterisation of it: the real uncertainty is at least this large, very possibly larger.
        Say so whenever this number is quoted; the module's own `displayResults` does.

        WHAT THIS MODULE DOES NOT DO: it runs NO physics of its own and depends on no process module, nor is
        it depended on by one (the same one-way relationship `NistCalibration` and the `Export*` family
        keep). The caller runs JAC exactly as they normally would, once per variant, under whichever settings
        (field choice, correlation space, grid, ...) define that variant, and hands the resulting values to
        this module purely for collection and summary. Because each variant costs one ordinary JAC
        calculation, the natural use is occasional and deliberate -- a handful of hand-picked, important
        results -- never a routine, automatic pass over every quantity a calculation produces.

        EXTENDING TO A NEW KIND OF VARIANT NEEDS NO NEW CODE HERE. A `Variant` is just a label and a value;
        whether the label names a field choice, a correlation layer, a box size, or something not yet tried
        is entirely the caller's business. A new jitter type is a new calling pattern in an example file, not
        a new code path in this module.
"""
module UncertaintyEstimate

using Printf


"""
`struct  UncertaintyEstimate.Variant`
    ... one alternative computation of a quantity, under one explicit, physically named choice.

    + label   ::String    ... a short, physically meaningful name for this variant, e.g. "DFS", "AL",
                                "ref-only", "ref+1layer" -- never a generic index or a random-draw number.
    + value   ::Float64   ... the quantity's computed value under this variant, in whatever unit the caller
                                is working in (the same unit for every Variant of one Estimate).
"""
struct  Variant
    label     ::String
    value     ::Float64
end


"""
`struct  UncertaintyEstimate.Estimate`
    ... the spread across a set of Variants of one quantity. See this module's own docstring for what this
        is and -- more importantly -- is NOT.

    + quantity    ::String             ... a short description of what was varied, e.g.
                                             "A(1s2p 1P1 -> 1s^2 1S0)".
    + variants    ::Array{Variant,1}   ... the variants this estimate was built from, in the order supplied.
    + spread      ::Float64            ... max(value) - min(value) over the variants.
    + relSpread   ::Float64            ... spread / mean(value) over the variants; dimensionless.
"""
struct  Estimate
    quantity      ::String
    variants      ::Array{Variant,1}
    spread        ::Float64
    relSpread     ::Float64
end


"""
`UncertaintyEstimate.estimate(quantity::String, variants::Array{UncertaintyEstimate.Variant,1})`
    ... to build an Estimate from at least two Variants of the same quantity. An est::Estimate is returned.
"""
function  estimate(quantity::String, variants::Array{Variant,1})
    length(variants) < 2   &&   error("UncertaintyEstimate.estimate(): need at least two variants to have a spread at all.")
    values    = [v.value for v in variants]
    spread    = maximum(values) - minimum(values)
    meanValue = sum(values) / length(values)
    meanValue == 0.0   &&   error("UncertaintyEstimate.estimate(): mean value is exactly zero; relSpread is undefined.")

    return( Estimate(quantity, variants, spread, spread/abs(meanValue)) )
end


"""
`UncertaintyEstimate.displayResults(stream::IO, est::UncertaintyEstimate.Estimate)`
    ... to print `est`'s variants and spread onto `stream`, with the lower-bound caveat this module's own
        docstring explains repeated alongside the number, not left to a reader's memory. Nothing is returned.
"""
function  displayResults(stream::IO, est::Estimate)
    nx = 70
    println(stream, " ");   println(stream, "  Quantity: ", est.quantity);   println(stream, " ")
    println(stream, "  ", "-"^nx)
    println(stream, "  ", rpad("Variant",20), "Value")
    println(stream, "  ", "-"^nx)
    for  v in est.variants
        println(stream, "  ", rpad(v.label,20), @sprintf("%.6e", v.value))
    end
    println(stream, "  ", "-"^nx)
    @printf(stream, "  spread = %.6e   (relative: %.4f %%)\n", est.spread, 100*est.relSpread)
    println(stream, "  THIS IS A LOWER BOUND on the true uncertainty, not a full statistical characterisation --")
    println(stream, "  only the $(length(est.variants)) variants named above were sampled. See this module's own docstring.")

    return( nothing )
end

end # module
