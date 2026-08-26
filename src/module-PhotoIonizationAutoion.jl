
"""
`module  JAC.PhotoIonizationAutoion`
... a submodule of JAC that was intended for photo-ionization followed by autoionization,
    |i(N)>  -->  |m(N-1)> + e_p  -->  |f(N-2)> + e_p + e_a.

    STATUS, 18-Aug-2026:  NOT IMPLEMENTED, AND NEVER HAS BEEN.  This is a declaration of intent and nothing more.
    The decision to leave it so is the maintainer's; do not begin an implementation here without asking.

    What was removed on that date, because it was dead rather than unfinished: a `Channel` type that nothing ever
    constructed -- carrying an `excitationChannel ::MultipoleAmplitude` field that was never written and never read --
    together with a `Pathway` type, its constructor and its `Base.show`, none of which any code path reached.  The
    `Settings` docstring also advertised an `incidentStokes ::ExpStokes` field that the struct did not have.

    What was NOT removed, and why: `Settings` is constructed by `TestFrames`, by `examples/example-Em.jl` and by the
    type list in `perform(::Atomic.Computation)`, and `computePathways` is called from that same dispatch chain.  Both
    therefore stay, but `computePathways` now RAISES instead of pretending.  Until 18-Aug-2026 it printed a green
    "The computation ... starts now" banner, then "Not yet implemented", and returned the STRING
    "Not yet implemented !" as its list of pathways -- which `perform` merged into its results under the ordinary key
    "photo-ionization-autoionization pathways:".  A caller got success-shaped output for a computation that had not
    happened, which is worse than a failure; cf. the `corePolarization.doApply` pattern of `module-PhotoEmission.jl`.

    Where the physics can be done today: the CASCADE route.  `Cascade.PhotoIonizationScheme` populates the ionic
    states and `Cascade.StepwiseDecayScheme` with `Auger()` carries the autoionization; both are implemented and
    exercised by examples.  See also `PhotoExcitationAutoion` for the excitation analogue.
"""
module PhotoIonizationAutoion

using ..Basics, ..ManyElectron, ..Radial


"""
`struct  PhotoIonizationAutoion.Settings`
    ... defines a type for the details and parameters of computing photo-ionization-autoionization pathways
        |i(N)>  -->  |m(N-1)>  -->  |f(N-2)>.  KEPT ONLY so that the dispatch chain, the struct tests and
        `examples/example-Em.jl` continue to resolve; no field of it is read by any computation.

    + multipoles              ::Array{EmMultipole,1}               ... Multipoles of the radiation field to be included.
    + gauges                  ::Array{UseGauge,1}                  ... Gauges to be included into the computations.
    + printBefore             ::Bool                               ... True, if all energies and lines are printed before evaluation.
    + selectPathways          ::Bool                               ... True if particular pathways are selected for the computations.
    + selectedPathways        ::Array{Tuple{Int64,Int64,Int64},1}  ... Selected pathways, as tuples (initial, intermediate, final).
    + maxKappa                ::Int64                              ... Maximum kappa value of partial waves to be included.
"""
struct Settings
    multipoles                ::Array{EmMultipole,1}
    gauges                    ::Array{UseGauge,1}
    printBefore               ::Bool
    selectPathways            ::Bool
    selectedPathways          ::Array{Tuple{Int64,Int64,Int64},1}
    maxKappa                  ::Int64
end


"""
`PhotoIonizationAutoion.Settings()`  ... constructor for the default values of photo-ionization-autoionization settings.
"""
function Settings()
    Settings( Basics.EmMultipole[], Basics.UseGauge[], false, false, Tuple{Int64,Int64,Int64}[], 0)
end


# `Base.show(io::IO, settings::PhotoIonizationAutoion.Settings)`
#        ... prepares a proper printout of the variable settings::PhotoIonizationAutoion.Settings.
function Base.show(io::IO, settings::PhotoIonizationAutoion.Settings)
    println(io, "multipoles:              $(settings.multipoles)  ")
    println(io, "gauges:                  $(settings.gauges)  ")
    println(io, "printBefore:             $(settings.printBefore)  ")
    println(io, "selectPathways:          $(settings.selectPathways)  ")
    println(io, "selectedPathways:        $(settings.selectedPathways)  ")
    println(io, "maxKappa:                $(settings.maxKappa)  ")
end


"""
`PhotoIonizationAutoion.computePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet,
                                        grid::Radial.Grid, settings::PhotoIonizationAutoion.Settings; output=true)`
    ... this module is not implemented; the call raises an error that says so. Nothing is returned.
"""
function  computePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid,
                            settings::PhotoIonizationAutoion.Settings; output=true)
    error("\n\nPhotoIonizationAutoion.computePathways() is NOT implemented and never has been.\n"                     *
          ">>> The module holds a Settings type and nothing else; no photo-ionization or autoionization\n"            *
          ">>> amplitude is computed for the pathway.  Until 18-Aug-2026 this function returned the STRING\n"          *
          ">>> \"Not yet implemented !\" as its pathway list, which reached the results dictionary looking\n"          *
          ">>> like an ordinary result.\n"                                                                            *
          ">>> For photo-ionization followed by autoionization, use the CASCADE route instead:\n"                     *
          ">>> Cascade.PhotoIonizationScheme to populate the ionic states, then Cascade.StepwiseDecayScheme\n"         *
          ">>> with Auger() for the autoionization -- both are implemented and exercised by examples.\n")
end

end # module
