
# Functions and methods for the two RESONANT channels of scheme::Cascade.ElectronIonizationScheme, in which the incident
# electron is CAPTURED into a doubly-excited resonance which then sheds two electrons, so that the ion ends up one charge
# state higher than it began:
#
#   resonant-electron-capture-with-sequential-double-autoionization     ResonantImpactIonization.SequentialAuger()
#       ...  i(N) + e-  ->  d(N+1)  ->  n(N) + e-  ->  f(N-1) + 2e-     two Auger steps through an autoionizing n
#   resonant-electron-capture-with-simultaneous-double-autoionization   ResonantImpactIonization.SimultaneousAuger()
#       ...  i(N) + e-  ->  d(N+1)  ->  f(N-1) + 2e-                    one double-Auger step
#
# The impact-excitation channel of the same scheme lives in module-Cascade-inc-electron-ionization.jl; these two are kept
# apart from it because they share nothing but the configuration lists: a capture is the time reverse of an Auger and
# needs no partial-wave sum over impact energies, which is why they cost seconds where that one costs hours.
#
# THE SIMULTANEOUS CHANNEL GENERATES NO STEPS OF ITS OWN.  It uses the resonance's total Auger width -- which the
# branching denominator needs anyway -- together with the shake-off estimate of ResonantImpactIonization, so requesting
# it costs nothing beyond the sequential channel's resonances.  It does, however, require those resonances to exist, so
# SimultaneousAuger() alone still builds the capture and decay steps below.


"""
`Cascade.determineResonantSteps(scheme::Cascade.ElectronIonizationScheme, comp::Cascade.Computation,
                                initialList::Array{Cascade.Block,1}, capturedList::Array{Cascade.Block,1},
                                excitedList::Array{Cascade.Block,1}, ionizedList::Array{Cascade.Block,1},
                                decayList::Array{Cascade.Block,1}; withSecondAuger::Bool=true)`
    ... determines all step::Cascade.Step's of the resonant channels of an electron-ionization cascade; a
        stepList::Array{Cascade.Step,1} is returned.

        FOUR FAMILIES OF STEP are set up, and each is needed for a different reason:

        + Auger, captured(N+1) -> initial(N)
            ... the resonance decaying back whence it came.  Its rate is what gives the CAPTURE rate, by detailed
                balance, and it is also the largest single term in the resonance's total width.
        + Auger, captured(N+1) -> excited(N)
            ... the first of the two sequential steps, into an intermediate that lies above the ionization threshold
                and can therefore autoionize again.
        + Auger, excited(N) -> ionized(N-1)
            ... the second sequential step, which completes the ionization.  The impact-excitation channel sets up the
                VERY SAME family for its own autoionization, so `withSecondAuger=false` suppresses it here whenever both
                channels are requested together.  (`Cascade.modifySteps` does NOT deduplicate -- it is a placeholder that
                copies every step through -- so this has to be handled at the point where the steps are built.)
        + Radiative, captured(N+1) -> decay(N+1)
            ... radiative stabilization of the resonance.  This is NOT an ionization route -- it is dielectronic
                recombination -- but it belongs in the resonance's TOTAL WIDTH, and leaving it out would make every
                ionization branching too large.  It is the reason a resonant-ionization computation also yields the
                recombination branching of the same resonances.

        Each step is set up only where it is energetically possible at all, i.e. where the mean energy of the upper
        block exceeds that of the lower one.

    + withSecondAuger  ::Bool
        ... whether family (3) is emitted.  Pass false when the impact-excitation channel runs alongside, since it emits
            the identical family and the amplitudes would otherwise be computed twice.
"""
function determineResonantSteps(scheme::Cascade.ElectronIonizationScheme, comp::Cascade.Computation,
                                initialList::Array{Cascade.Block,1}, capturedList::Array{Cascade.Block,1},
                                excitedList::Array{Cascade.Block,1}, ionizedList::Array{Cascade.Block,1},
                                decayList::Array{Cascade.Block,1}; withSecondAuger::Bool=true)
    stepList = Cascade.Step[]
    if  comp.approach  in  [Cascade.AverageSCA(), Cascade.SCA()]
        maxKappa = length(scheme.lValues) > 0 ? maximum(scheme.lValues) + 1 : 7
        aSettings = AutoIonization.Settings(AutoIonization.Settings(), maxKappa=maxKappa)
        rSettings = PhotoEmission.Settings(PhotoEmission.Settings(), multipoles=[E1], gauges=[UseCoulomb, UseBabushkin])
        #
        # (1) the resonance decaying back to the initial block; this supplies the capture rate by detailed balance
        for  capturedBlock in capturedList
            for  initialBlock in initialList
                if  capturedBlock.NoElectrons == initialBlock.NoElectrons + 1
                    push!( stepList, Cascade.Step(Basics.Auger(), aSettings, capturedBlock.confs,     initialBlock.confs,
                                                                             capturedBlock.multiplet, initialBlock.multiplet) )
                end
            end
        end
        # (2) the first sequential Auger step, into an autoionizing intermediate of the original charge state
        for  capturedBlock in capturedList
            for  excitedBlock in excitedList
                if  capturedBlock.NoElectrons == excitedBlock.NoElectrons + 1   &&
                    Basics.determineMeanEnergy(capturedBlock.multiplet) - Basics.determineMeanEnergy(excitedBlock.multiplet) > 0.
                    push!( stepList, Cascade.Step(Basics.Auger(), aSettings, capturedBlock.confs,     excitedBlock.confs,
                                                                             capturedBlock.multiplet, excitedBlock.multiplet) )
                end
            end
        end
        # (3) the second sequential Auger step, which completes the ionization; suppressed when the impact-excitation
        #     channel is running too, since that one emits the identical family
        for  excitedBlock in (withSecondAuger ? excitedList : Cascade.Block[])
            for  ionizedBlock in ionizedList
                if  excitedBlock.NoElectrons == ionizedBlock.NoElectrons + 1   &&
                    Basics.determineMeanEnergy(excitedBlock.multiplet) - Basics.determineMeanEnergy(ionizedBlock.multiplet) > 0.
                    push!( stepList, Cascade.Step(Basics.Auger(), aSettings, excitedBlock.confs,     ionizedBlock.confs,
                                                                             excitedBlock.multiplet, ionizedBlock.multiplet) )
                end
            end
        end
        # (4) radiative stabilization of the resonance; not an ionization route, but part of its total width
        for  capturedBlock in capturedList
            for  decayBlock in decayList
                if  capturedBlock.NoElectrons == decayBlock.NoElectrons   &&
                    Basics.determineMeanEnergy(capturedBlock.multiplet) - Basics.determineMeanEnergy(decayBlock.multiplet) > 0.
                    push!( stepList, Cascade.Step(Basics.Radiative(), rSettings, capturedBlock.confs,     decayBlock.confs,
                                                                                 capturedBlock.multiplet, decayBlock.multiplet) )
                end
            end
        end
        #
    else  error("Unsupported cascade approach.")
    end

    return( stepList )
end


"""
`Cascade.generateConfigurationsForResonantIonization(multiplets::Array{Multiplet,1}, scheme::Cascade.ElectronIonizationScheme,
                                                     nm::Nuclear.Model)`
    ... generates the DOUBLY-EXCITED configurations of the resonant channels, i.e. those reached by capturing the incident
        electron into one of `scheme.intoShells` while exciting a further electron from `scheme.excitationFromShells` into
        `scheme.excitationToShells`, each carrying one electron MORE than the initial configurations.

        A Tuple(capturedConfList, decayConfList) is returned: the doubly-excited resonances themselves, and the
        (N+1)-electron configurations they can radiate INTO, which are needed for the resonance's radiative width.

        The generation is delegated to `Basics.ForDielectronicRecombination`, the same configuration theme that
        `Cascade.DielectronicRecombinationScheme` uses -- these resonances ARE dielectronic-capture resonances, and only
        their subsequent decay distinguishes the two schemes, so building them by a second route would be an invitation
        for the two to drift apart.

        The intermediate and final configurations of these channels are NOT generated here: they are the very same
        inner-shell excited and ionized lists that the impact-excitation channel already builds, which is what lets the
        two channels share one scheme at all.

        An explanatory error is raised if no intoShells were given -- a resonant channel without a capture shell has
        nothing to capture into, and would otherwise fail later and less clearly.
"""
function generateConfigurationsForResonantIonization(multiplets::Array{Multiplet,1}, scheme::Cascade.ElectronIonizationScheme,
                                                     nm::Nuclear.Model)
    if  length(scheme.intoShells) == 0
        error("Cascade.generateConfigurationsForResonantIonization(): a resonant channel was requested but scheme.intoShells " *
              "is empty.  The incident electron has to be captured into some shell; give at least one, e.g. " *
              "intoShells = [Shell(\"3s\"), Shell(\"3p\"), Shell(\"3d\")].")
    end
    initialConfList = Basics.extractConfigurations(Basics.FromMultiplet(), multiplets)
    # The decay shells are taken to be the capture shells: the excited core electron drops back while the captured
    # electron stays where it is, which is the dominant radiative stabilization of such a resonance.
    theme = Basics.ForDielectronicRecombination(scheme.excitationFromShells, scheme.excitationToShells,
                                                scheme.intoShells, scheme.intoShells)
    capturedConfList, decayConfList = Basics.generateConfigurations(theme, initialConfList)

    return( unique(capturedConfList), unique(decayConfList) )
end


"""
`Cascade.hasResonantChannel(scheme::Cascade.ElectronIonizationScheme)`
    ... returns true::Bool if either of the two resonant channels appears in `scheme.processes`, and false otherwise.
"""
function hasResonantChannel(scheme::Cascade.ElectronIonizationScheme)
    return( ResonantImpactIonization.SequentialAuger()   in  scheme.processes   ||
            ResonantImpactIonization.SimultaneousAuger() in  scheme.processes )
end
