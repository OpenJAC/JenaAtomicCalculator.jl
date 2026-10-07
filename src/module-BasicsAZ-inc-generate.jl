
using   QuadGK, ..Hamiltonian
export  generate


"""
`Basics.generate(representation::AtomicState.Representation)`  
    ... to generate an atomic representation as specified by the representation.repType::AbstractRepresentationType.
        All relevant intermediate and final results are printed to screen (stdout). Nothing is returned.

`Basics.generate(representation::AtomicState.Representation; output::Bool=false)`  
    ... to generate the same but to return the complete output in a dictionary; the particular output depends on the type and 
        specifications of the representation but can easily accessed by the keys of this dictionary.
"""
function Basics.generate(representation::AtomicState.Representation; output::Bool=false)
    results = Basics.generate(representation.repType, representation; output=output)
    
    Defaults.warn(PrintWarnings())
    Defaults.warn(ResetWarnings())
    println(" ")
    return( results isa Dict{String,Any} ? Basics.PerformResults(results) : results )
end



"""
`Basics.generate(repType::AtomicState.MeanFieldBasis, representation::AtomicState.Representation)`  
    ... to generate a mean-field basis (representation) for a set of reference configurations; all relevant intermediate 
        and final results are printed to screen (stdout). Nothing is returned.

`Basics.generate(repType::AtomicState.MeanFieldBasis, representation::AtomicState.Representation; output::Bool=false)`  
    ... to generate the same but to return the complete output in a dictionary; the particular output depends on the type and 
        specifications of the representation but can easily accessed by the keys of this dictionary.
"""
function Basics.generate(repType::AtomicState.MeanFieldBasis, rep::AtomicState.Representation; output::Bool=false)
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    nModel    = rep.nuclearModel

    # The asfSettings only define the SCF part and are partly derived from the MeanFieldSettings
    asfSettings = AsfSettings(AsfSettings(); scField=repType.settings.scField) 
    
    multiplet      = SelfConsistent.performSCF(rep.refConfigs, nModel, rep.grid, asfSettings; printout=true)
    basis          = multiplet.levels[1].basis
    if output    results = Base.merge( results, Dict("mean-field basis" => basis) )          end
    
    return( results isa Dict{String,Any} ? Basics.PerformResults(results) : results )
end


"""
`Basics.generate(repType::AtomicState.MeanFieldMultiplet, representation::AtomicState.Representation)`  
    ... to generate a mean-field basis (representation) for a set of reference configurations; all relevant intermediate 
        and final results are printed to screen (stdout). Nothing is returned.

`Basics.generate(repType::AtomicState.MeanFieldMultiplet, representation::AtomicState.Representation; output::Bool=false)`  
    ... to generate the same but to return the complete output in a dictionary; the particular output depends on the type and 
        specifications of the representation but can easily accessed by the keys of this dictionary.
"""
function Basics.generate(repType::AtomicState.MeanFieldMultiplet, rep::AtomicState.Representation; output::Bool=false)
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    nModel    = rep.nuclearModel

    # The asfSettings only define the SCF part and are partly derived from the MeanFieldSettings
    asfSettings = AsfSettings(AsfSettings(); scField=repType.settings.scField) 
    
    multiplet  = SelfConsistent.performSCF(rep.refConfigs, nModel, rep.grid, asfSettings; printout=true)
    if output    
        ## `levels`, not `level` (fixed 06-Aug-2026): this typo made generate(MeanFieldMultiplet, ...) raise a
        ## FieldError for every call with output=true, i.e. the only useful call. It survived because the one
        ## example that exercised this path -- the old branch 1 of example-Dh.jl -- had never been run.
        results = Base.merge( results, Dict("mean-field basis" => multiplet.levels[1].basis) )
        results = Base.merge( results, Dict("mean-field multiplet" => multiplet) )          
    end
    
    return( results isa Dict{String,Any} ? Basics.PerformResults(results) : results )
end




"""
`Basics.generate(repType::AtomicState.OneElectronSpectrum, rep::AtomicState.Representation)`  
    ... to generate a one-electron spectrum for the atomic potential from the (given) levels, based on a set of reference 
        configurations as well as for given settings. Relevant intermediate and final results are printed to screen (stdout). 
        Nothing is returned in this case.

    + `(repType::AtomicState.OneElectronSpectrum, representation::AtomicState.Representation; output=true)`  
    ... to generate the same but to return the complete output in a orbitals::Dict{Subshell, Orbital}.
"""
function Basics.generate(repType::AtomicState.OneElectronSpectrum, rep::AtomicState.Representation; output::Bool=true)
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    nModel    = rep.nuclearModel
    settings  = repType.settings
    # First perform a SCF+CI computations for the reference configurations below to generate a spectrum of start orbitals
    asfSettings   = AsfSettings()   ## Use default settings to define a first multiplet from the reference configurations;
                                    ## all further details are specified for each step
    refMultiplet  = SelfConsistent.performSCF(rep.refConfigs, nModel, rep.grid, asfSettings; printout=true)
    refBasis      = refMultiplet.levels[1].basis
    Basics.display(stdout, refBasis.orbitals, rep.grid; longTable=false)
    nuclearPot    = Nuclear.nuclearPotential(nModel, rep.grid)
    @warn("The potential is generated for the mean-field basis but not (yetc hosen for the selected levels.")
    ## A NuclearField carries no electronic part at all; adding a DFS screening term on top of it would silently
    ## reintroduce the very mismatch this field exists to remove.
    if  typeof(settings.scField) == Basics.NuclearField
        electronicPot = Radial.Potential("zero", zeros(length(rep.grid.r)), rep.grid)
    else
        electronicPot = Basics.computePotential(Basics.DFSField(1.0), rep.grid, refBasis)
    end
    meanPot       = Basics.add(nuclearPot, electronicPot)
    
    println("")
    printstyled("Compute an one-electron spectrum for the selected partial waves ... \n", color=:light_green)
    printstyled("------------------------------------------------------------------- \n", color=:light_green)
    
    # Generate all non-relativistic and relativistic subshells and the associated single-electron spectrum for this potential
    shellList = Basics.generateShellList(1, settings.nMax, settings.lValues)
    subshellList = Subshell[]
    for  shell in shellList     append!(subshellList, Basics.shellSplitIntoSubshells(shell))    end
    primitives = Bsplines.generatePrimitives(rep.grid)
    orbitals   = Bsplines.generateOrbitals(subshellList, meanPot, nModel, primitives; printout=true)
    
    # Print all results to screen
    Basics.display(stdout, orbitals, rep.grid; longTable=true)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    Basics.display(iostream, orbitals, rep.grid; longTable=false)         end
    
    if output    
        results = Base.merge( results, Dict("mean potential" => meanPot) )              
        results = Base.merge( results, Dict("orbitals" => orbitals) )              
    end
    
    return( results isa Dict{String,Any} ? Basics.PerformResults(results) : results )
end



"""
`Basics.generate(repType::AtomicState.CiExpansion, representation::AtomicState.Representation)`  
    ... to generate a configuration-interaction expansion for a single level symmetry and based on a set of reference configurations
        and a number of pre-specified steps. All relevant intermediate and final results are printed to screen (stdout). 
        Nothing is returned.

`Basics.generate(repType::AtomicState.CiExpansion, representation::AtomicState.Representation; output::Bool=false)`  
    ... to generate the same but to return the complete output in a dictionary; the particular output depends on the type and 
        specifications of the computations but can easily accessed by the keys of this dictionary.
"""
function Basics.generate(repType::AtomicState.CiExpansion, rep::AtomicState.Representation; output::Bool=false)
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    nModel    = rep.nuclearModel
    orbitals  = repType.applyOrbitals
    
    # Generate a list of relativistic configurations and  CSF's for the given subshell list
    relconfList = ConfigurationR[]
    for  conf in rep.refConfigs
        wa = Basics.generateConfigurations(Basics.RelativisticConfigurations(), conf)
        append!( relconfList, wa)
    end
    subshellList = Basics.generateSubshellList(relconfList)
    csfList = CsfR[]
    for  relconf in relconfList
        newCsfs = Basics.generateCsfRs(relconf, subshellList)
        append!( csfList, newCsfs)
    end
    
    # Generate all level symmetries that need to be taken into account
    if     length(repType.settings.levelSelectionCI.symmetries) != 0    symmetries = repType.settings.levelSelectionCI.symmetries
    else   symmetries = LevelSymmetry[]
        for csf in csfList
            if      LevelSymmetry(csf.J, csf.parity)  in  symmetries
            else    push!( symmetries, LevelSymmetry(csf.J, csf.parity) )
            end
        end
    end
    println("*** Level symmetries = $symmetries ")

    # The asfSettings only define the CI part and are partly derived from the CiSettings
    asfSettings = AsfSettings(true, CoulombInteraction(), Basics.DFSField(), StartFromHydrogenic(), Basics.AutomaticRoute(), 0., Subshell[], Subshell[], 
                                1.0e-3, true, repType.settings.eeInteractionCI, NoneQed(), LSjjSettings(false), 
                                repType.settings.levelSelectionCI) 
    
    basis      = Basics.generateBasis(rep.refConfigs, symmetries, repType.excitations)
    basis      = Basis( true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, orbitals )  
    # Test that all orbitals are available 
    for  subshell in  basis.subshells   
        if   !(haskey(orbitals, subshell))  error("No orbital available for subshell $subshell ")     end
    end
    multiplet  = Hamiltonian.performCI(basis,  nModel, rep.grid, asfSettings; printout=true) 
    if output    results = Base.merge( results, Dict("CI multiplet" => Multiplet("CI multiplet:", multiplet.levels)) )              end
    
    return( results isa Dict{String,Any} ? Basics.PerformResults(results) : results )
end



"""
`Basics.printRasStepDiagnostic(istep::Int64, multiplet::Multiplet, basis::Basis, refConfigs::Array{Configuration,1},`
                              `frozenShells::Array{Shell,1}, grid::Radial.Grid)`
    ... prints, for one step of a RAS expansion, the two numbers that say what the layer actually did: the weight
        of each of the lowest levels on the REFERENCE CSFs, and the mean radius of every subshell the step newly
        introduced. Nothing is returned. See priority item 25 for why a step that prints only its multiplet is
        not enough.
"""
function Basics.printRasStepDiagnostic(istep::Int64, multiplet::Multiplet, basis::Basis,
                                       refConfigs::Array{Configuration,1}, frozenShells::Array{Shell,1},
                                       grid::Radial.Grid)
    levels = sort(multiplet.levels, by = l -> l.energy)
    isempty(levels)  &&  return( nothing )
    refIdx = SelfConsistent.referenceCsfIndices(levels[1].basis, refConfigs)
    println("\n>> [RAS step $istep] what this layer did:")
    if  isempty(refIdx)
        println(">>   NO CSF of this basis belongs to the reference configurations -- the reference space does " *
                "not describe this step, and every weight below would be zero.")
    else
        for  (i, level)  in  enumerate(levels[1:min(5, end)])
            w = sum( level.mc[r]^2  for r in refIdx )
            println(">>   level $i   J^P = " * string(LevelSymmetry(level.J, level.parity)) *
                    "   E = " * @sprintf("%.8f", level.energy) * " Ha   weight on the reference CSFs = " *
                    @sprintf("%.5f", w) * (w < 0.5 ? "   <-- NOT a reference level" : ""))
        end
    end
    # THE SHELLS THIS STEP ACTUALLY VARIED are the ones it did NOT freeze, which is the same rule the driver
    # uses to build AsfSettings.frozenSubshells.  Taking it from the frozen list rather than from a "new
    # shells" field keeps this correct for any RasStep, however its excitation lists are written.
    orbitals = levels[1].basis.orbitals
    for  subsh  in  sort(collect(keys(orbitals)), by = t -> (t.n, abs(t.kappa)))
        isFrozen = any(sh -> sh.n == subsh.n  &&  Basics.subshell_l(subsh) == sh.l, frozenShells)
        if  !isFrozen  &&  length(orbitals[subsh].P) > 1
            r0 = RadialIntegrals.rkDiagonal(0, orbitals[subsh], orbitals[subsh], grid)
            r1 = RadialIntegrals.rkDiagonal(1, orbitals[subsh], orbitals[subsh], grid)
            println(">>   varied subshell " * string(subsh) * "   <r> = " * @sprintf("%.4f", r1/r0) *
                    " a.u.   norm = " * @sprintf("%.6f", r0))
        end
    end
    return( nothing )
end


"""
`Basics.generate(repType::AtomicState.RasExpansion, representation::AtomicState.Representation)`  
    ... to generate a restricted active-space expansion for a single level symmetry and based on a set of reference configurations
        and a number of pre-specified steps. All relevant intermediate and final results are printed to screen (stdout). 
        Nothing is returned.

`Basics.generate(repType::AtomicState.RasExpansion, representation::AtomicState.Representation; output::Bool=false)`  
    ... to generate the same but to return the complete output in a dictionary; the particular output depends on the type and 
        specifications of the computations but can easily accessed by the keys of this dictionary.
"""
function Basics.generate(repType::AtomicState.RasExpansion, rep::AtomicState.Representation; output::Bool=false)
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    nModel   = rep.nuclearModel
    # First perform a SCF+CI computations for the reference configurations below to generate a spectrum of start orbitals
    asfSettings    = AsfSettings()  ## Use default settings to define a first multiplet from the reference configurations;
                                    ## all further details are specified for each step
    priorMultiplet = SelfConsistent.performSCF(rep.refConfigs, nModel, rep.grid, asfSettings; printout=true)
    nuclearPot     = Nuclear.nuclearPotential(nModel, rep.grid)
    # ITEM 23, FIXED 01-Sep-2026.  The two lines that build a SCREENED potential sat here commented out, and the
    # start orbitals were generated in the BARE nuclear potential -- so every correlation orbital a layer
    # introduces began life as a hydrogenic function of the FULL nuclear charge.  MEASURED at Z = 76 with 60
    # electrons: the generated spectrum reproduced the point-Dirac energies to 1.9e-05 relative for 2p_1/2 and
    # 1e-10 for 2p_3/2, i.e. the screening of sixty electrons was exactly absent, and a 5f correlation orbital
    # started far too compact for the valence region it is meant to correlate.
    #   THE SAME DEFECT WAS FIXED IN THE GreenExpansion DRIVER ON 07-Aug-2026, a hundred lines below in this
    # file, and the reasoning recorded there -- that a spectrum must match the Hamiltonian of the computation
    # that consumes it -- applies here word for word.  This is the wiring job that reasoning implies, using the
    # reference multiplet already computed just above as the density for the screening.
    electronicPot  = Basics.computePotential(Basics.DFSField(1.0), rep.grid, priorMultiplet.levels[1].basis)
    meanPot        = Basics.add(nuclearPot, electronicPot)
    subshellList   = Basics.extractRelativisticSubshellList(rep)             ## extract all subshells that occur in the RAS computation
    primitives     = Bsplines.generatePrimitives(rep.grid)
    startOrbitals  = Bsplines.generateOrbitals(subshellList, meanPot, nModel, primitives, printout=true)  ## generate a spectrum of sufficient size
    if output    results = Base.merge( results, Dict("reference multiplet" => Multiplet("Reference multiplet:", priorMultiplet.levels) ) )  end

    # Now, cycle over all steps of the RasExpansion. Each step re-optimizes ONLY the shells newly introduced
    # at this layer (EOL, so the SAME target level(s) specified once in repType.settings.levelSelectionCI stay
    # the optimization target as the CSF space grows -- an AL/configuration-average functional would let
    # higher-lying states entering the growing space dilute the correlation for the levels actually wanted),
    # freezing every shell already present/optimized in an earlier step -- GRASP's own default behavior
    # (LFIX(:NW)=.TRUE. unless explicitly varied).
    #
    # WHAT THE AL ALTERNATIVE ACTUALLY COSTS, measured 11-Sep-2026, because "dilute" understates it. An average-
    # level field over a correlation expansion optimizes the core for the configurations that carry a core HOLE:
    # on a Cl III basis whose reference is 2s^2 it returns a mean occupation of 1.807 for 2s, so the reference is
    # not a stationary point of its own orbitals. A single excitation between orbitals of EQUAL kappa is then no
    # longer suppressed -- the one-particle operator is rank 0, so it connects only equal kappa, and 2s -> 4s is
    # an orbital ROTATION rather than a correlation excitation. Ranked by a second-order estimate it took 99.7 %
    # of the whole core contribution, against 71.6 % for the genuine double 2p^2 -> 3d^2 once it was removed.
    # Under EOL the target stays the reference-dominated level and that channel is absorbed into the orbitals.
    # Same-kappa correlation orbitals cannot be avoided in any case -- layers n = 4,5,6,7 give 4s,5s,6s,7s -- so
    # what protects the expansion is the CHOICE OF FUNCTIONAL here, not the choice of shells.
    # priorIsLayer says whether priorMultiplet is a previous LAYER or still the reference SCF of line 254:
    # the boundary check below is meaningful only between two layers, since the reference SCF runs a
    # different field and a first layer may legitimately sit above it.
    priorIsLayer = false
    for (istep, step)  in  enumerate(repType.steps)
        println("")
        printstyled("++ Compute the orbitals, orbitals and multiplet for step $istep ... \n", color=:light_green)
        printstyled("--------------------------------------------------------------      \n", color=:light_green)
        # THE FREEZING POLICY THIS DRIVER ALREADY CLAIMS, NOW ACTUALLY APPLIED -- 06-Oct-2026, closing the cause of
        # challenge S18.  The comment above says a step re-optimizes only its new shells, "freezing every shell
        # already present/optimized in an earlier step -- GRASP's own default behavior".  That was NOT what
        # happened: AtomicState.RasExpansion accumulates its frozen set from coreShells, fromShells and each
        # layer's newShells ONLY, so a reference shell named in NEITHER coreShells NOR fromShells was never frozen
        # and was re-optimized at every layer -- by performSCF's AVERAGE-LEVEL pass, i.e. for the configuration
        # average rather than for the level energy the layer minimizes.
        #   MEASURED: 43Ca+ [Ar] 4s with 3s -> 5s,6s singles handed layer 2 a point 36.1 mHa ABOVE where layer 1
        # had finished, and layer 2 then ended 2.8 mHa WORSE than its own reference -- which a variational layer
        # cannot be.  Freezing those five shells puts the entry back AT the previous layer's energy and the ladder
        # descends again.  See the boundary check further down, which reports the symptom when it survives.
        #   WHY IT WENT UNSEEN: in examples/example-Ai.jl the Be reference 1s^2 2s^2 is exactly coreShells=[1s]
        # plus fromShells=[2s], so the claimed policy and the implemented one coincide there.
        #   THE FIRST STEP IS LEFT FREE, deliberately: it is the step that optimizes the reference, and there is no
        # earlier step whose work could be lost.
        #   IT IS DERIVED HERE AND NOT IN RasExpansion because that constructor never receives refConfigs, and the
        # reference's occupied shells are what the policy is about.  A hand-built RasStep still carries whatever
        # frozen set it is given, so the other policies Rule 21 names ("core frozen", "all free") remain reachable
        # that way; what changes is only the DEFAULT of the layer API, which now matches its own documentation.
        stepFrozenShells = step.frozenShells
        if  istep > 1
            refShells = Shell[]
            for  conf in rep.refConfigs,  (sh, occ) in conf.shells
                if  occ > 0  &&  !(sh in refShells)    push!(refShells, sh)    end
            end
            extraFrozen = Shell[ sh  for sh in refShells
                                 if !any(f -> f.n == sh.n  &&  f.l == sh.l, step.frozenShells) ]
            if  !isempty(extraFrozen)
                stepFrozenShells = vcat(step.frozenShells, extraFrozen)
                println(">> [RAS] step $istep also freezes the reference shells " *
                        join(string.(extraFrozen), ", ") * ", which an earlier step already optimized;  without " *
                        "this they would be re-optimized here for the configuration average and the layer could " *
                        "end ABOVE the layer it sits on.")
            end
        end
        basis      = Basics.generateBasis(rep.refConfigs, repType.symmetries, step)
        orbitals   = Basics.generateOrbitalsForBasis(basis, stepFrozenShells, priorMultiplet.levels[1].basis, startOrbitals)
        basis      = Basis( true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, orbitals )

        # step.frozenShells is a list of non-relativistic Shell(n,l); translate to the concrete, relativistic
        # Subshell(n,kappa) instances actually present in THIS step's basis, as AsfSettings.frozenSubshells needs.
        # Subshell[...] and not [...]:  a step with NO frozen shells makes this comprehension empty, Julia then
        # infers Vector{Any}, and AsfSettings rejects it.  Unreachable through RasLayer, which always passes a
        # non-empty core, and reached at once by building the steps by hand -- which is what a step that adds new
        # EXCITATIONS rather than new SHELLS has to do.
        frozenSubshellsThisStep = Subshell[ sh  for  shell in stepFrozenShells  for sh in basis.subshells
                                                if  sh.n == shell.n  &&  Basics.subshell_l(sh) == shell.l ]

        # ITEM 22, FIXED 01-Sep-2026.  RasSettings.levelsScf was declared, documented and PRINTED, and never read:
        # the EOL target came from levelSelectionCI alone, so a user who wrote RasSettings([1], ...) to optimize
        # the layer on the ground level got, silently, a functional over every level the CI had selected.  It is
        # now read -- but only when levelSelectionCI is INACTIVE, because levelSelectionCI is the more expressive
        # and the safer selector: it can name levels by their weight on the reference configurations, which is
        # the only one of the three that survives an intruder sinking below the reference, whereas levelsScf can
        # only name indices into the energy-sorted multiplet.  Silently replacing a reference-weight selection by
        # an index one would be a downgrade, so when both are given levelSelectionCI wins and this says so.
        stepLevelSelection = repType.settings.levelSelectionCI
        if  !stepLevelSelection.active  &&  !isempty(repType.settings.levelsScf)
            stepLevelSelection = LevelSelection(true, indices=repType.settings.levelsScf)
            if  istep == 1
                println(">> [RAS] the EOL target of every step is levelsScf = $(repType.settings.levelsScf), by index.")
            end
        elseif  stepLevelSelection.active  &&  !isempty(repType.settings.levelsScf)  &&  istep == 1
            println(">> [RAS] BOTH levelsScf and levelSelectionCI were given; levelSelectionCI is in force as the " *
                    "EOL target and levelsScf = $(repType.settings.levelsScf) is ignored.")
        end
        # THE ROUTE IS PASSED THROUGH, NOT REBUILT.  Until 04-Oct-2026 this line read
        #     scfRoute = Basics.RotationRoute(repType.settings.maxIterationsScf)
        # which converted an Int field into a route here, one line deep in the driver -- so a RAS ladder could name only the
        # ROTATION route and only its iteration COUNT, while RotationRoute also carries nVirtual and stepping and FockRoute
        # carries unscaledOffDiagonal.  Four of the six things a route can say were unreachable from a RAS computation, and the
        # route hierarchy exists precisely so a run may be repeated on a dearer route when a cheaper one does not converge.
        # A STEP MAY OVERRIDE THE LADDER, because a ladder's steps are not alike: a reference layer is converged in about twelve
        # iterations where a correlation layer never converges at all, so one budget cannot serve both.
        stepRoute    = isnothing(step.scfRoute) ? repType.settings.scfRoute : step.scfRoute
        if  !isnothing(step.scfRoute)
            println(">> [RAS] step $istep takes its OWN scf route, $(step.scfRoute), instead of the ladder's " *
                    "$(repType.settings.scfRoute).")
        end
        stepSettings = AsfSettings( AsfSettings();  scField=Basics.EOLField(),  frozenSubshells=frozenSubshellsThisStep,
                                     eeInteractionCI=repType.settings.eeInteractionCI,  levelSelectionCI=stepLevelSelection,
                                     scfRoute=stepRoute,
                                     accuracyScf=repType.settings.accuracyScf )

        # A VARIATIONAL step optimizes its new shells and puts every configuration it generated into the CI, which
        # is what a step has always done.  A SECOND-ORDER step does neither: its orbitals are the ones this layer
        # was given -- earlier layers frozen, the new shells taken from the spectrum generated once in the mean
        # potential above -- and its configurations are ranked, then promoted, folded in or dropped.  Skipping the
        # SCF is the point rather than an economy: the layer exists to RECOVER what P leaves out, so optimizing
        # orbitals on the whole of it would pay for Q in the one place the treatment is meant to avoid, and the
        # virtuals of a perturbation belong in the field of the reference in any case.
        if  typeof(step.treatment) == Basics.Variational
            multiplet  = SelfConsistent.performSCF(basis, nModel, rep.grid, stepSettings; printout=true)
            # THE RELEASE PHASE -- RasSettings.freezing, 06-Oct-2026.  A released policy runs the layer TWICE and
            # that is what makes it trustworthy rather than merely richer.
            #   Phase 1 above WARMED this layer's new shells with the earlier ones frozen, so the new correlation
            # orbitals have a sensible shape.  Phase 2 now RELEASES, and starts the optimized-level solver from those
            # warm orbitals through StartFromPrevious -- which (since 06-Oct-2026) skips the average-level pre-pass.
            # Skipping it is the whole point: that pass optimizes the configuration AVERAGE, a different functional
            # from the level energy, and would re-derive the inherited orbitals before the solver ever saw them.
            #   MEASURED without the two-step entry, 43Ca+ [Ar] 4s with 3s -> 5s,6s singles: a released layer BEGAN
            # 36.1 mHa above the layer it sits on -- more than the whole correlation effect -- and its increment then
            # had a budget-dependent SIGN, +11.2, +5.9, +1.0 mHa at budgets 20, 60, 120.  With the two-step entry the
            # same layer gives -0.013, -0.053, -0.063 mHa: monotone, correctly signed, and saturating.
            #   THE POLICY GOVERNS WHAT A LAYER DOES WITH WHAT EARLIER LAYERS OPTIMIZED, so it applies from step 2 on.
            # Step 1 optimizes the reference and has nothing to inherit, and its coreShells stay the caller's
            # explicit declaration under every policy.
            if  repType.settings.freezing != :previousLayersFrozen  &&  istep > 1
                eWarm        = minimum( lv.energy  for lv in multiplet.levels )
                warmOrbitals = multiplet.levels[1].basis.orbitals
                # :coreFrozen keeps the expansion's own coreShells, which are exactly the frozen set of step 1;
                # :allFree keeps nothing.
                releaseShells = repType.settings.freezing == :coreFrozen ? repType.steps[1].frozenShells : Shell[]
                releaseSubshells = Subshell[ sh  for  shell in releaseShells  for sh in basis.subshells
                                             if  sh.n == shell.n  &&  Basics.subshell_l(sh) == shell.l ]
                printstyled(">> [RAS] step $istep RELEASE phase (freezing = $(repType.settings.freezing)): the " *
                            "warm-up reached " * @sprintf("%.9f", eWarm) * " Ha and the solver now restarts from " *
                            "those orbitals with " *
                            (isempty(releaseSubshells) ? "nothing frozen" :
                             "only " * join(string.(releaseSubshells), ", ") * " frozen") * ". \n",
                            color=:light_cyan)
                warmBasis   = Basis( true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells,
                                     warmOrbitals )
                relSettings = AsfSettings( stepSettings;  frozenSubshells = releaseSubshells,
                                           startScfFrom = ManyElectron.StartFromPrevious(warmOrbitals) )
                multiplet   = SelfConsistent.performSCF(warmBasis, nModel, rep.grid, relSettings; printout=true)
                eRel        = minimum( lv.energy  for lv in multiplet.levels )
                printstyled(">> [RAS] step $istep RELEASE gained " * @sprintf("%.4f mHa", (eRel - eWarm)*1000) *
                            " over its own warm-up" *
                            (eRel > eWarm ? " -- which is POSITIVE and should not be: the release began above the " *
                                            "warm-up it started from, so suspect the grid or the level selection." :
                                            ".") * " \n", color=:light_cyan)
            end
        else
            printstyled(">> step $istep is treated to SECOND ORDER in its Q space; the orbitals of this layer are " *
                        "NOT re-optimized. \n", color=:light_yellow)
            (multiplet, partitions) = Hamiltonian.performCIwithSecondOrderQ(basis, priorMultiplet.levels[1].basis,
                                                rep.refConfigs, nModel, rep.grid, stepSettings, step.treatment;
                                                printout=true)
            if  !all(p -> p.isJustified, partitions)
                printstyled(">> step $istep: the second-order treatment was NOT justified in every symmetry; see " *
                            "the configurations named above. \n", color=:light_red)
            end
        end
        # ITEM 25, ADDED 01-Sep-2026.  Until now a step printed one Multiplet and nothing else, so nothing in
        # the output said whether the layer just added is a CORRELATION layer or has collapsed onto a
        # spectroscopic orbital -- a distinction that is not academic: a 5f treated as occupied once came out
        # at 4f's radius, mixed in at 22 %, and reordered a multiplet by 3.2 eV, caught only because the
        # application computed the reference weight by hand.  Two numbers per step say it.  The weight is on
        # the REFERENCE CSFs, so a level that drifts away from the reference space shows up as w_ref falling;
        # the radii say whether a new shell contracted into the valence region (a correlation orbital) or
        # expanded out of it (a Rydberg orbital, which is not what a layer is for).
        # THE DOUBLING TEST, REPORTED PER LAYER AND COSTING NOTHING -- 07-Oct-2026.  A correlation layer never truly
        # converges, so its iteration count is a cost dial and the honest question is whether DOUBLING the effort
        # would still move the layer's increment.  That used to need a second computation at half the budget;  the
        # solver now carries the descent at half its own iterations, so the answer comes from the run in hand -- and
        # more correctly, since both numbers lie on ONE trajectory rather than on two runs whose earlier layers may
        # have stopped in different places.
        #   IT IS ADVISORY AND NEVER GATING, as the verdict itself is.  A layer whose increment is still moving is
        # not wrong, it is unfinished, and the reader is the one who decides whether that matters.
        if  repType.settings.scfRoute isa Basics.RotationRoute
            vd = SelfConsistent.lastScfVerdict()
            if  !isnothing(vd)  &&  vd.descent > 0.  &&  vd.descentAtHalf > 0.
                drift = abs(vd.descent - vd.descentAtHalf) / vd.descent
                # TWO DIFFERENT NUMBERS, AND THEY MUST NOT BE CONFLATED -- which this print did on the day it was
                # written.  A layer's INCREMENT is everything it gains over the previous layer, and most of it
                # arrives at the ENTRY, from mixing the new CSFs in at the inherited orbitals: that part is one
                # diagonalization, exact, and cannot be unconverged.  The solver's DESCENT is only the part the
                # orbital optimization then adds, and it is the only part an iteration count can affect.  Measured
                # on C II 1s^2 2s^2 2p + 3s3p3d(SD): the increment is -80.42 mHa and the orbital optimization
                # contributes -3.26 mHa of it, i.e. 4 %.  So the doubling test applies to the DESCENT -- correctly,
                # since that is the only thing still moving -- and a tolerance of 5 % on it is worth about 1 % on
                # the increment.  Printing both keeps that visible instead of inviting the reader to read the
                # tolerance as an error bar on the layer.
                layerInc = istep > 1 ? (minimum(lv.energy for lv in multiplet.levels) -
                                        minimum(lv.energy for lv in priorMultiplet.levels)) : NaN
                println(">> [RAS] step $istep: increment " *
                        (isnan(layerInc) ? "n/a (first step)" : @sprintf("%+.4f mHa", 1000*layerInc)) *
                        ", of which the orbital optimization gave " * @sprintf("%+.4f mHa", -1000*vd.descent) *
                        " in $(vd.iterations) iterations;  the rest arrived at entry, from the CI, and needs no " *
                        "iteration.")
                # A RELATIVE TEST ON A NEGLIGIBLE QUANTITY SAYS NOTHING, and saying it anyway is worse than silence.
                # Seen on the first run that printed this: a REFERENCE layer whose orbitals were already converged
                # gained 0.0002 mHa by optimizing them, and a 20 % change in THAT was reported as "not converged" --
                # alarming, and about nothing.  So the test is only applied where the optimized part is large enough
                # to matter: above a MICRO-HARTREE in absolute terms -- a thousandth of a milli-Hartree, below which
                # no quoted atomic energy is affected -- and above a thousandth of the layer's own increment in
                # relative ones.  The absolute floor was tried at 1e-9 Ha first and was useless: the two reference
                # layers that provoked this gained 1e-7 and 2e-7 Ha, so they sailed over it and were still judged.
                negligible = vd.descent < 1.0e-6  ||
                             (!isnan(layerInc)  &&  vd.descent < 1.0e-3 * abs(layerInc))
                if  negligible
                    println(">> [RAS] step $istep: the orbital optimization moved this layer by only " *
                            @sprintf("%+.4f mHa", -1000*vd.descent) * ", so there was nothing left for it to do " *
                            "and the doubling test does not apply.")
                else
                println(">> [RAS] step $istep: the DOUBLING TEST on that optimized part -- " *
                        @sprintf("%+.4f mHa", -1000*vd.descent) * " against " *
                        @sprintf("%+.4f mHa", -1000*vd.descentAtHalf) * " at half the iterations, a change of " *
                        @sprintf("%.2f %%", 100*drift) * " -- " *
                        (drift <= 0.10 ? "so this layer is CONVERGED in its budget." :
                                         "so this layer is NOT converged;  raise the budget, or quote the " *
                                         "increment WITH it."))
                end
            end
        end
        Basics.printRasStepDiagnostic(istep, multiplet, basis, rep.refConfigs, stepFrozenShells, rep.grid)
        # A LAYER THAT ENDS ABOVE THE LAYER IT SITS ON IS NOT PHYSICS, AND UNTIL 06-Oct-2026 NOTHING SAID SO.
        # A variational layer's CSF space CONTAINS the previous layer's, so its lowest level cannot lie above the
        # previous one: taking the previous orbitals and the mixing vector (1,0,...) is a point of this layer's own
        # parameter space.  When it does lie above, something has gone wrong in the HANDOVER rather than in the
        # minimization, and this is the one place that can see it.
        #
        # MEASURED, and it is what this check was written from.  43Ca+ [Ar] 4s with 3s -> 5s,6s singles and nothing
        # frozen: layer 1 converged to -679.5212307 and layer 2 ENTERED at -679.4850961, i.e. 36.1 mHa above it,
        # ending at -679.5184105 -- a correlation layer 2.8 mHa WORSE than its own reference.  The cause is named
        # below and is NOT the optimizer: the second-order route drives |grad| from 1.04 to 0.0016, a factor 640
        # better than the first-order one, and lands HIGHER still (-679.5156563), because a tighter solve converges
        # to the stationary point of whatever basin it is handed.
        #
        # WHY THE HANDOVER LOSES IT.  performSCF runs an AVERAGE-LEVEL pass before the optimized-level solver, and
        # that pass re-optimizes every subshell not named in AsfSettings.frozenSubshells -- for the configuration
        # AVERAGE, which is a different functional from the level energy the layer is supposed to minimize.  So a
        # reference shell that no earlier step froze is re-optimized at every layer, discarding what the previous
        # layer achieved.  AtomicState.RasExpansion accumulates its frozen set from coreShells, fromShells and each
        # layer's newShells only, so a reference shell in NEITHER coreShells NOR fromShells is never frozen.  In
        # examples/example-Ai.jl the Be reference 1s^2 2s^2 is exactly coreShells=[1s] plus fromShells=[2s], the
        # two coincide, and this is invisible -- which is why it went unseen.
        #
        # THE CHECK MOVES NO NUMBER.  It compares two energies that were computed anyway and prints;  the remedy is
        # the caller's, and naming the unfrozen reference shells is what makes it actionable, since those are
        # exactly the shells to pass as coreShells.
        if  typeof(step.treatment) == Basics.Variational  &&  priorIsLayer
            eNow  = minimum( lv.energy  for lv in multiplet.levels )
            ePrev = minimum( lv.energy  for lv in priorMultiplet.levels )
            if  eNow > ePrev
                refShells = Shell[]
                for  conf in rep.refConfigs,  (sh, occ) in conf.shells
                    if  occ > 0  &&  !(sh in refShells)    push!(refShells, sh)    end
                end
                unfrozen = Shell[ sh  for sh in refShells
                                  if !any(f -> f.n == sh.n  &&  f.l == sh.l, stepFrozenShells) ]
                sa = ">> [RAS] STEP $istep ENDED ABOVE STEP $(istep-1) BY " *
                     @sprintf("%.3f mHa", (eNow - ePrev)*1000) * " -- a variational layer CANNOT do that, since " *
                     "its CSF space contains the previous layer's.  The layer was handed the previous orbitals and " *
                     "then lost ground before it began to descend."
                printstyled(sa * "\n", color=:light_red)
                if  !isempty(unfrozen)
                    sb = ">> [RAS]   These reference shells are frozen in NO earlier step, so they are " *
                         "re-optimized here -- and performSCF's average-level pass does that for the configuration " *
                         "AVERAGE, not for the level energy this layer minimizes:  " *
                         join(string.(unfrozen), ", ") * ".\n" *
                         ">> [RAS]   Pass them as coreShells to keep what the earlier layers achieved, or read " *
                         "this layer's increment as a lower bound rather than as the correlation energy."
                    printstyled(sb * "\n", color=:light_red)
                else
                    sb = ">> [RAS]   Every reference shell IS frozen here, so the loss is not the freezing policy; " *
                         "suspect the grid (a correlation shell that does not fit cannot lower the energy) or the " *
                         "level selection (a different level may have become the target)."
                    printstyled(sb * "\n", color=:light_red)
                end
                Defaults.warn(AddWarning(), "Basics.generate(RasExpansion): step $istep ended " *
                              @sprintf("%.3f mHa", (eNow - ePrev)*1000) * " ABOVE step $(istep-1), which a " *
                              "variational layer cannot do.  Unfrozen reference shells: " *
                              (isempty(unfrozen) ? "none" : join(string.(unfrozen), ", ")) * ".")
            end
        end
        if output    results = Base.merge( results, Dict("step"*string(istep) => Multiplet("Multiplet:", multiplet.levels)) )              end
        # A PERTURBATIVE STEP IS A LEAF AND MUST NOT BECOME THE TRUNK.  It reports what the excitations left out of
        # the previous step are still worth, over that step's own orbitals; a later variational step therefore
        # inherits from the last VARIATIONAL step, not from this one.  Were the leaf carried forward, the next step
        # would take the whole folded space as its P -- including exactly the configurations that were deliberately
        # kept out of the CI.  So a chain may read step2, step2PT, step3, step3PT, ... and each PT branches off its
        # own layer.
        if  typeof(step.treatment) == Basics.Variational    priorMultiplet = multiplet;   priorIsLayer = true    end
    end
    
    return( results isa Dict{String,Any} ? Basics.PerformResults(results) : results )
end



"""
`Basics.generate(repType::AtomicState.GreenExpansion, representation::AtomicState.Representation)`  
    ... to generate a Green (function) expansion for a given approach and excitation scheme of the electron,
        based on a set of reference configurations, a list of level symmetries as well as for given settings.
        All relevant intermediate and final results are printed to screen (stdout). Nothing is returned.

`Basics.generate(repType::AtomicState.GreenExpansion, representation::AtomicState.Representation; output::Bool=false)`  
    ... to generate the same but to return the complete output in a dictionary; the particular output depends on the type and 
        specifications of the representation but can easily accessed by the keys of this dictionary.
"""
function Basics.generate(repType::AtomicState.GreenExpansion, rep::AtomicState.Representation; output::Bool=false)
    if  output    results = Dict{String, Any}()    else    results = nothing    end
    nModel    = rep.nuclearModel
    settings  = repType.settings
    # First perform a SCF+CI computations for the reference configurations below to generate a spectrum of start orbitals
    ## THE POTENTIAL IS NOW TAKEN FROM settings.scField (07-Aug-2026), not hardcoded. Both places below used to
    ## be fixed to DFS -- the reference AsfSettings and the mean potential -- so a Green expansion could never be
    ## matched to the potential of the computation that consumes it. That matters because a Green expansion is
    ## normally used as the INTERMEDIATE spectrum of a second-order calculation, where gauge invariance requires
    ## the same one-body Hamiltonian throughout.
    asfSettings   = AsfSettings(AsfSettings(); scField=settings.scField)
    refMultiplet  = SelfConsistent.performSCF(rep.refConfigs, nModel, rep.grid, asfSettings; printout=true)
    refBasis      = refMultiplet.levels[1].basis
    nuclearPot    = Nuclear.nuclearPotential(nModel, rep.grid)
    ## A NuclearField carries no electronic part at all; adding a DFS screening term on top of it would silently
    ## reintroduce the very mismatch this field exists to remove.
    if  typeof(settings.scField) == Basics.NuclearField
        electronicPot = Radial.Potential("zero", zeros(length(rep.grid.r)), rep.grid)
    else
        electronicPot = Basics.computePotential(Basics.DFSField(1.0), rep.grid, refBasis)
    end
    meanPot       = Basics.add(nuclearPot, electronicPot)
    
    println("")
    printstyled("Compute an approximate Green function expansion ... \n", color=:light_green)
    printstyled("--------------------------------------------------- \n", color=:light_green)
    
    channels = AtomicState.GreenChannel[]
    
    # Generate all (non-relativistic) configurations from the bound configurations due to the given excitation scheme 
    confList = Basics.generateConfigurationsForExcitationScheme(rep.refConfigs, repType.excitationScheme, settings.nMax, settings.lValues)
    # Print (if required) information about the generated configuration list
    if  settings.printBefore    Basics.displayConfigurations(stdout, confList)    end
    
    # Generate shell list abd a full single-electron spectrum for this potential
    subshellList = Basics.extractRelativisticSubshellList(confList)                      ## extract all subshells that occur in confList
    primitives   = Bsplines.generatePrimitives(rep.grid)
    orbitals     = Bsplines.generateOrbitals(subshellList, meanPot, nModel, primitives, printout=true) ## generate a spectrum of sufficient size
    Defaults.setDefaults("relativistic subshell list", subshellList; printout=true)
    Basics.display(stdout, orbitals, rep.grid)

    # The asfSettings only define the CI part of the Green channels and are partly derived from the GreenSettings
    asfSettings = AsfSettings(true, CoulombInteraction(), Basics.DFSField(), StartFromHydrogenic(), Basics.AutomaticRoute(), 0., Subshell[], Subshell[], 
                                1.0e-3, true, CoulombInteraction(), NoneQed(), LSjjSettings(false), settings.levelSelection ) 
    
    # Cycle over all selected level symmetries to generate the requested channels
    for  levelSymmetry  in  repType.levelSymmetries
        basis      = Basics.generateBasis(confList, [levelSymmetry])
        basis      = Basis( true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, orbitals )  
        multiplet  = Basics.computeMultipletForGreenApproach(repType.approach, basis, nModel, rep.grid, asfSettings, settings; printout=true) 
        push!( channels, GreenChannel( levelSymmetry, multiplet) )
    end        
    
    # Print all results to screen
    Basics.display(stdout, channels)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    Basics.display(iostream, channels)         end
    
    if output    results = Base.merge( results, Dict("Green channels" => channels) )   end
    return( results isa Dict{String,Any} ? Basics.PerformResults(results) : results )
end



"""
`Basics.generate(::CondensedMultiplet, multiplet::Multiplet)`
    ... to condense/reduce the number of CSF in the basis of the given multiplet due to a single 'weight';
        a multiplet::Multiplet is returned.  **Not yet implemented !**
"""
function Basics.generate(::CondensedMultiplet, multiplet::Multiplet)
    error("Not yet implemented !")
end



"""
`Basics.generate(::ConfigurationListNRFromBasis, basis::Basis)`
    ... to (re-) generate the list of NR configurations from the given basis; a confList::Array{Configuration,1} is returned.
"""
function Basics.generate(::ConfigurationListNRFromBasis, basis::Basis)
    confList    = Configuration[]
    NoElectrons = sum( basis.csfs[1].occupation )
    for  csf in basis.csfs
        shellList = Dict{Shell,Int64}[]
        for  i = 1:length(csf.occupation)
            n = csf.subshell[i].n;    l = Basics.subshell_l( csf.subshell[i] );    sh = Shell(n, l)
            if     haskey(shellList, sh)    shellList(sh) = shellList(sh) + csf.occupation[i]
            else   shellList = merge( shellList, Dict(sh => csf.occupation[i]) )
            end
        end
        confNew = Configuration( shellList, NoElectrons )
        push!( confList, confNew )
    end
    confListNew = ConfigurationsExcludeDoubles(confList)
    
    return( confListNew )
end


"""
`Basics.generate(::ConfigurationListNRFromConfiguration, refConf::Configuration, NoExcitations::Int64, fromShells::Array{Shell,1},`
                    toShells::Array{Shell,1})
    ... to generate a non-relativistic configuration list, including the given reference configuration (refConf) and with
        all configurations that differ by NoExcitations from the fromShells into the toShells; an Array{Configuration,1}
        is returned.
"""
function Basics.generate(::ConfigurationListNRFromConfiguration, refConf::Configuration, NoExcitations::Int64, fromShells::Array{Shell,1}, toShells::Array{Shell,1})
    confList = [refConf]
    # First prepare a proper reference configuration that also includes all fromShells and toShells with zero occupation
    shellDict = deepcopy(refConf.shells)
    allShells = deepcopy(fromShells)
    for  sh in toShells
        if     sh in allShells
        else   append!(allShells, [sh] )
        end
    end
    
    for  sh in allShells
    if  !haskey(shellDict, sh )    shellDict = merge( shellDict, Dict( sh => 0 ) )    end
    end

    rconf = Configuration(shellDict, refConf.NoElectrons)
    println("aa", rconf)

    # Add single excitations
    if      NoExcitations == 0
    elseif  NoExcitations == 1
        for  fromsh in fromShells
            for  tosh in toShells
                confnew = deepcopy(rconf)
                confnew.shells[fromsh] = confnew.shells[fromsh] - 1
                confnew.shells[tosh]   = confnew.shells[tosh]   + 1
                println("bb", confnew)

                # Check that this is valid configuration and append, if appropriate
                add = true
                for (k,v) in  confnew.shells
                    if  v <= -1  ||  Parity(rconf) != Parity(confnew)    add = false    end
                end
                if  add 
                    for (k,v) in confnew.shells   
                    if  v == 0   delete!(confnew.shells, k)   end
                    end
                    push!(confList, confnew);   println("fromsh = $fromsh  tosh = $tosh")   
                end
            end
        end

    # Add double excitations
    elseif  NoExcitations == 2
        for  fromsha in fromShells
        for  fromshb in fromShells
            for  tosha in toShells
            for  toshb in toShells
                confnew = deepcopy(rconf)
                confnew.shells[fromsha] = confnew.shells[fromsha] - 1
                confnew.shells[fromshb] = confnew.shells[fromshb] - 1
                confnew.shells[tosha]   = confnew.shells[tosha]   + 1
                confnew.shells[toshb]   = confnew.shells[toshb]   + 1
                println("bb", confnew)

                # Check that this is valid configuration and append, if appropriate
                add = true
                for (k,v) in  confnew.shells
                    if  v <= -1  ||  Parity(rconf) != Parity(confnew)    add = false    end
                end
                if  add    
                    for (k,v) in confnew.shells   
                    if  v == 0   delete!(confnew.shells, k)   end
                    end
                    push!(confList, confnew)    
                end
            end
            end
        end
        end

    else
        error("stop a")
    end

    # Exclude all configurations that appear twice (and more) in the list
    confList = Basics.ConfigurationsExcludeDoubles(confList)

    return( confList )
end


"""
`Basics.generateCsfRs(conf::ConfigurationR, subshellList::Array{Subshell,1})` 
    ... to construct from a given (relativistic) configuration all possible CSF with regard to the subshell order as specified 
        by subshellList; a list::Array{CsfR,1} is returned.
"""
function Basics.generateCsfRs(conf::ConfigurationR, subshellList::Array{Subshell,1})
    parity  = Basics.extractFromConfiguration(Basics.GetParity(), conf)
    csfList = CsfR[];   useStandardSubshells = true;    first = true;    previousCsfs = CsfR[]
    # subhshellList = Subshell[];   
    for  subsh in subshellList
        if   subsh in keys(conf.subshells)    occ = conf.subshells[subsh]    else    occ = 0    end
        if   first
            stateList   = ManyElectron.provideSubshellStates(subsh, occ)
            currentCsfs = CsfR[]
            for  state in stateList
                push!( currentCsfs, CsfR( true, AngularJ64(state.Jsub2//2), parity, [state.occ], [state.nu],
                                            [AngularJ64(state.Jsub2//2)], [AngularJ64(state.Jsub2//2)], Subshell[]) )
            end
            previousCsfs = copy(currentCsfs)
            first        = false
        else
            # Now support also all couplings of the subshell states with the CSFs that were built-up so far
            stateList   = ManyElectron.provideSubshellStates(subsh, occ)
            currentCsfs = CsfR[]
            for  csf in  previousCsfs
                for  state in stateList
                    occupation = deepcopy(csf.occupation);    seniorityNr = deepcopy(csf.seniorityNr);    
                    subshellJ  = deepcopy(csf.subshellJ);     subshells = deepcopy(csf.subshells)
                    push!(occupation, state.occ);   push!(seniorityNr, state.nu);   push!(subshellJ, AngularJ64(state.Jsub2//2) ) 
                    push!(subshells, subsh)
                    newXList = oplus( csf.subshellX[end], AngularJ64(state.Jsub2//2) )
                    for  newX in newXList
                        subshellX = deepcopy(csf.subshellX);   push!(subshellX, newX) 
                        push!( currentCsfs, CsfR( true, subshellX[end], parity, occupation, seniorityNr, subshellJ, subshellX, Subshell[]) ) 
                    end
                end
            end
            previousCsfs = copy(currentCsfs)
        end
    end
    
    return( previousCsfs )
end
    


"""
`Basics.generate(::OrderedShellList, confs::Array{Configuration,1})`
    ... to generate for confs, i.e. all the given (non-relativistic) configurations, a common and ordered shell list;
        a list::Array{Shell,1} is returned.
"""
function Basics.generate(::OrderedShellList, confs::Array{Configuration,1})
    shells = Shell[]

    wa = Defaults.getDefaults("ordered shell list: non-relativistic", 11)
    ## wa = Defaults.getDefaults("ordered shell list: non-relativistic", 19)
    for  a in wa
        for  cf in 1:length(confs)
            ks = keys(confs[cf].shells)
            if  a in ks   push!(shells, a);    break    end
        end 
    end

    return( shells )
end


"""
`Basics.generateSubshellList(confs::Array{ConfigurationR,1})`  
    ... to generate for confs, i.e. all the given (relativistic) configurations, common and ordered subshell list; 
        a list::Array{Subshell,1} is returned.
"""
function Basics.generateSubshellList(confs::Array{ConfigurationR,1})
    subshells = Subshell[]   

    for  conf in confs
        for  subsh in keys(conf.subshells)      push!(subshells, subsh)     end
    end
    subshells = Base.unique(subshells)
    subshells = Base.sort( subshells, lt=Base.isless)

    return( subshells )
end



"""
`Basics.generate(::OrderedSubshellList, basisA::Basis, basisB::Basis)`
    ... to generate common and ordered subshell list for the two basis A and B; a list::Array{Subshell,1} is returned.
"""
function Basics.generate(::OrderedSubshellList, basisA::Basis, basisB::Basis)
    function areEqual(nn::Int64, sha::Array{Subshell,1}, shb::Array{Subshell,1})
        # Determines whether the first nn subshells are equal in sha and shb (true) or not (false)
        for  i = 1:nx
            if    sha[i] != shb[i]    return( false )   end
        end
        return( true )
    end

    subshells = Subshell[]

    nx = min(length(basisA.subshells), length(basisB.subshells))
    if  areEqual(nx, basisA.subshells, basisB.subshells)
        # If subshell order is equal, any subshell order is accepted for those subshells that occur in both basis
        for  i = 1:nx   
            if    basisA.subshells[i] != basisB.subshells[i]   error("Inconsistent subshells of two bases.")
            else  push!( subshells, basisA.subshells[i])
            end
        end
        #
        if       length(basisA.subshells) > nx   
            for  i = nx+1:length(basisA.subshells)    push!(subshells, basisA.subshells[i])    end
        elseif   length(basisB.subshells) > nx   
            for  i = nx+1:length(basisB.subshells)    push!(subshells, basisB.subshells[i])    end
        end
    else
        ## println("basisA.subshells = $(basisA.subshells)")
        ## println("basisB.subshells = $(basisB.subshells)")
        # If subshell order is NOT equal, all subshell in basisA and basisB must follow standard order
        standardList = Defaults.getDefaults("ordered subshell list: relativistic", 7)
        na = 0;  nb = 0
        for  sh in standardList   
            if  sh in basisA.subshells  ||   sh in basisB.subshells   push!( subshells, sh)    end    
        end
        nn = 0;   
        for sh in basisA.subshells   
            if    !(sh in subshells) push!( subshells, sh)    
            else  wb = findall(x->x==sh, subshells);    if wb[1] <= nn   error("stop a")   else   nn = wb[1]  end    
            end
        end
        nn = 0;   
        for sh in basisB.subshells   
            if    !(sh in subshells) push!( subshells, sh)    
            else  wb = findall(x->x==sh, subshells);    if wb[1] <= nn   error("stop a")   else   nn = wb[1]  end    
            end
        end
        println(">>> Extended subshells from two basis = $subshells")
    end

    return( subshells )
end


"""
`Basics.generate(::SlaterTypeSpectrum, N::Int64, potential::Radial.Potential, grid::Radial.Grid; N_0::Int64=30, alpha_0::Float64=1.0,`
                    beta_0::Float64=1.1)
    ... to generate a complete one-electron spectrum with N positive and N negative states, and by using even-tempered Slater-type
        orbitals (STO) with parameters ``\\alpha_i = \\alpha_0 \\beta_0^i``; a spectrum::SingleElecSpectrum is returned where just
        N0 positive and N_0 negative are kept for later use.  **Not yet implemented !**
"""
function Basics.generate(::SlaterTypeSpectrum, N::Int64, potential::Radial.Potential, grid::Radial.Grid; N_0::Int64=30, alpha_0::Float64=1.0, beta_0::Float64=1.1)
    error("Not yet implemented !")
    return( nothing )
end


"""
`Basics.generate(::SlaterTypeSpectrumPositive, N::Int64, potential::Radial.Potential, grid::Radial.Grid; N_0::Int64=30, alpha_0::Float64=1.0,`
                    beta_0::Float64=1.1)
    ... to generate the same but to return only the N_0 positive states.  **Not yet implemented !**
"""
function Basics.generate(::SlaterTypeSpectrumPositive, N::Int64, potential::Radial.Potential, grid::Radial.Grid; N_0::Int64=30, alpha_0::Float64=1.0, beta_0::Float64=1.1)
    error("Not yet implemented !")
    return( nothing )
end



"""
`Basics.generateBasis(confList::Array{Configuration,1}, symmetries::Array{LevelSymmetry,1})`  
    ... generates the CSF basis for the given configuration list and level symmetries; a basis::Basis is returned but 
        without a valid representation of the radial orbitals.
"""
function Basics.generateBasis(confList::Array{Configuration,1}, symmetries::Array{LevelSymmetry,1})
    #
    relconfList = ConfigurationR[]
    for  conf in confList
        wa = Basics.generateConfigurations(Basics.RelativisticConfigurations(), conf)
        append!( relconfList, wa)
    end
    subshellList = Basics.generateSubshellList(relconfList)
    Defaults.setDefaults("relativistic subshell list", subshellList; printout=true)

    # Generate the relativistic CSF's for the given subshell list
    csfList = CsfR[]
    for  relconf in relconfList
        newCsfs = Basics.generateCsfRs(relconf, subshellList)
        for  csf in newCsfs     if  LevelSymmetry(csf.J, csf.parity) in symmetries   push!(csfList, csf)   end   end
    end
    #
    if length(csfList) == 0     error("There are no CSF with $symmetries in the given configuration list.")      end

    # Determine the number of electrons and the list of coreSubshells
    NoElectrons      = sum( csfList[1].occupation )
    coreSubshellList = Subshell[]
    for  k in 1:length(subshellList)
        mocc = Basics.subshell_2j(subshellList[k]) + 1;    is_filled = true
        for  csf in csfList
            if  csf.occupation[k] != mocc    is_filled = false;    break   end
        end
        if   is_filled    push!( coreSubshellList, subshellList[k])    end
    end
    
    basis = Basis(true, NoElectrons, subshellList, csfList, coreSubshellList, Dict{Subshell, Orbital}())
    
    println("Construct a basis with $(length(basis.csfs)) CSF for J^P = $symmetries with $(length(basis.subshells)) subshells: " *
            "$(basis.subshells[1])  $(basis.subshells[2]) ...  $(basis.subshells[end-1])  $(basis.subshells[end])")
    return( basis )
end


"""
`Basics.generateBasis(refConfigs::Array{Configuration,1}, symmetries::Array{LevelSymmetry,1}, step::AtomicState.RasStep)`  
    ... generates the CSF basis for the given reference configurations and single, double, ... excitations of
        electrons from the corresponding fromShells --> toShells. A basis::Basis is returned but without a valid
        representation of the radial orbitals
"""
function Basics.generateBasis(refConfigs::Array{Configuration,1}, symmetries::Array{LevelSymmetry,1}, step::AtomicState.RasStep)
    # Single excitations
    if      step.seFrom == Shell[]  ||  step.seTo == Shell[]    confSingles = Configuration[]
    else    confSingles = Basics.generateConfigurations(refConfigs, step.seFrom, step.seTo)
    end
    # Double excitations
    if      step.deFrom == Shell[]  ||  step.deTo == Shell[]    confDoubles = Configuration[]
    else    confDoubles = Basics.generateConfigurations(refConfigs, step.deFrom, step.deTo)
            confDoubles = Basics.generateConfigurations(confDoubles,     step.deFrom, step.deTo)
    end
    # Triple excitations
    if      step.teFrom == Shell[]  ||  step.teTo == Shell[]    confTriples = Configuration[]
    else    confTriples = Basics.generateConfigurations(refConfigs,      step.teFrom, step.teTo)
            confTriples = Basics.generateConfigurations(confTriples,     step.teFrom, step.teTo)
            confTriples = Basics.generateConfigurations(confTriples,     step.teFrom, step.teTo)
    end
    # Quadruple excitations.  The guard tests qeFrom/qeTo and NOT teFrom/teTo -- it tested the triples' fields
    # until 11-Sep-2026, so a step asking for quadruples ALONE silently produced none, while a step asking for
    # triples ran this branch with empty shell lists.  Invisible because nothing in the repository sets qeFrom.
    if      step.qeFrom == Shell[]  ||  step.qeTo == Shell[]    confQuadruples = Configuration[]
    else    confQuadruples = Basics.generateConfigurations(refConfigs, step.qeFrom, step.qeTo)
            confQuadruples = Basics.generateConfigurations(confQuadruples,  step.qeFrom, step.qeTo)
            confQuadruples = Basics.generateConfigurations(confQuadruples,  step.qeFrom, step.qeTo)
            confQuadruples = Basics.generateConfigurations(confQuadruples,  step.qeFrom, step.qeTo)
    end
    
    # Now get a unique set of configurations and generate the relativistic configurations
    # configurations = Base.unique((==), configurations) does not work here
    confList       = Configuration[]  
    configurations = deepcopy(refConfigs);      append!(configurations, confSingles);    append!(configurations, confDoubles)
                                                append!(configurations, confTriples);    append!(configurations, confQuadruples)
    for  confa in configurations
        addTo = true
        for confb in confList    if   confa == confb    addTo = false;    break     end     end
        if  addTo    push!(confList, confa)     end
    end
    for  i = 1:length(confList)    println(">> include ", confList[i])    end
    #
    relconfList = ConfigurationR[]
    for  conf in confList
        wa = Basics.generateConfigurations(Basics.RelativisticConfigurations(), conf)
        append!( relconfList, wa)
    end
    subshellList = Basics.generateSubshellList(relconfList)
    Defaults.setDefaults("relativistic subshell list", subshellList; printout=true)

    # Generate the relativistic CSF's for the given subshell list
    csfList = CsfR[]
    for  relconf in relconfList
        newCsfs = Basics.generateCsfRs(relconf, subshellList)
        for  csf in newCsfs     if  LevelSymmetry(csf.J, csf.parity) in symmetries   push!(csfList, csf)   end   end
    end

    # Determine the number of electrons and the list of coreSubshells
    NoElectrons      = sum( csfList[1].occupation )
    coreSubshellList = Subshell[]
    for  k in 1:length(subshellList)
        mocc = Basics.subshell_2j(subshellList[k]) + 1;    is_filled = true
        for  csf in csfList
            if  csf.occupation[k] != mocc    is_filled = false;    break   end
        end
        if   is_filled    push!( coreSubshellList, subshellList[k])    end
    end
    
    basis = Basis(true, NoElectrons, subshellList, csfList, coreSubshellList, Dict{Subshell, Orbital}())
    
    println("Construct a basis with $(length(basis.csfs)) CSF for J^P = $symmetries with $(length(basis.subshells)) subshells: " *
            "$(basis.subshells[1])  $(basis.subshells[2]) ...  $(basis.subshells[end-1])  $(basis.subshells[end])")
    return( basis )
end


"""
`Basics.generateConfigurations(refConfigs::Array{Configuration,1}, fromShells::Array{Shell,1}, toShells::Array{Shell,1})`  
    ... generates all nonrelativistic configurations due to the (single) excitation of an electron fromShells --> toShells
        for all given reference configurations. A confList::Array{Configuration,1} is returned.
"""
function Basics.generateConfigurations(refConfigs::Array{Configuration,1}, fromShells::Array{Shell,1}, toShells::Array{Shell,1})
    confList = Configuration[]
    for  config in refConfigs
        for  fromShell in fromShells
            for  toShell in toShells
                # Cycle through all shells of config to built-up a new configuration with a single excitation if fromShell occurs
                if  haskey(config.shells, fromShell )
                    newShells = Dict{Shell,Int64}();   addShell = true;   addConfiguration = true
                    for (k,v) in config.shells
                        if        k == fromShell == toShell       newShells = Base.merge( newShells, Dict( k => v))
                                    addShell = false   
                        elseif    k == fromShell    if  v-1 > 0   newShells = Base.merge( newShells, Dict( fromShell => v-1))  end
                        elseif    k == toShell      
                            if    v+1 <= 2*(2*toShell.l + 1)      newShells = Base.merge( newShells, Dict( toShell => v+1))
                                    addShell = false   
                            else  addConfiguration = false
                            end
                        else      newShells = Base.merge( newShells, Dict( k => v))
                        end
                        #
                        if  addShell                              newShells = Base.merge( newShells, Dict( toShell => 1))      end
                    end
                    if  addConfiguration   push!(confList, Configuration(newShells, config.NoElectrons ))          end
                else    # No configuration to add if fromShell does not occur in the reference configuration
                end
            end
        end
    end
    
    confList = Base.unique(confList)
    
    return( confList )
end



"""
`Basics.generateConfigurations(refConfigs::Array{Configuration,1}, fromShells::Array{Shell,1}, toShells::Array{Shell,1}, 
                                noex::Int64; restrictions::Array{AbstractConfigurationRestriction,1}=AbstractConfigurationRestriction[])`  
    ... generates all nonrelativistic configurations with excitation of up to noex electrons fromShells --> toShells
        for all given reference configurations. Moreover, the list of restrictions (if any is given) is finally applied 
        to restrict the configurations to a given set of limitations. It remains the reponsibility of the user to make sure that
        the given restrictions do not contradict each other and are consistent with what is to be achieved. The given set of 
        restrictions can be easily extended if this need arises by the users. A confList::Array{Configuration,1} is returned.
"""
function Basics.generateConfigurations(refConfigs::Array{Configuration,1}, fromShells::Array{Shell,1}, toShells::Array{Shell,1}, 
                                        noex::Int64; restrictions::Array{AbstractConfigurationRestriction,1}=AbstractConfigurationRestriction[])
    if     noex == 0
        newConfigs = refConfigs
    elseif noex > 0
        newConfigs = Basics.generateConfigurations(refConfigs, fromShells, toShells)
        newConfigs = Base.unique(newConfigs)
        for ne = 2:noex    
            newConfigs = Basics.generateConfigurations(newConfigs, fromShells, toShells)    
            newConfigs = Base.unique(newConfigs)
        end
    end
    
    # Now apply in turn all given restrictions, if any, and append if no restriction is violated
    confList = Configuration[]
    for  conf  in  newConfigs
        addConf = true
        for res in restrictions
            if  Basics.isViolated(conf,res)     addConf = false;    break   end
        end
        if  addConf     push!(confList, conf)   end
    end
    
    return( confList )
end



"""
`Basics.generateConfigurationsForExcitationScheme(confs::Array{Configuration,1}, exScheme::Basics.NoExcitationScheme, 
                                                    nMax::Int64, lValues::Array{Int64,1})`  
    ... generates a list of non-relativistic configurations for the given (reference) confs and the excitation scheme. 
        All orbitals in standard order are considered up to maxShell = (n_max, l_max).
"""
function Basics.generateConfigurationsForExcitationScheme(confs::Array{Configuration,1}, exScheme::Basics.NoExcitationScheme, 
                                                        nMax::Int64, lValues::Array{Int64,1})
    @warn(  "No excitations are included if exScheme::NoExcitationScheme. ")
    println("No excitations are included if exScheme::NoExcitationScheme. ")
    #
    newConfList = confs
    return( newConfList )
end



"""
`Basics.generateConfigurationsForExcitationScheme(confs::Array{Configuration,1}, exScheme::Basics.DeExciteSingleElectron, 
                                                    nMax::Int64, lValues::Array{Int64,1})`  
    ... generates a list of non-relativistic configurations for the given (reference) confs and the excitation scheme. 
        All orbitals in standard order are considered up to maxShell = (n_max, l_max).
"""
function Basics.generateConfigurationsForExcitationScheme(confs::Array{Configuration,1}, exScheme::Basics.DeExciteSingleElectron, 
                                                            nMax::Int64, lValues::Array{Int64,1})
    confList  = Configuration[];    NoElectrons = confs[1].NoElectrons
    shellList = Basics.generateShellList(confs, nMax, lValues)
    
    # Create all configurations with single excitations from the given list of configurations
    for  conf in confs
        # Take one electron fromShell and add toShell
        for  (fromShell,occ)  in  conf.shells
            for  toShell  in  shellList
                newShells = deepcopy( conf.shells )
                if      fromShell == toShell    
                elseif  haskey(conf.shells, toShell )   fromOcc = occ;   toOcc = conf.shells[toShell]
                        if  fromOcc - 1 < 0                        println("..");    continue    end
                        if  toOcc   + 1 > 2*(2*toShell.l + 1)      println("..");    continue    end
                        newShells[fromShell] = newShells[fromShell] - 1
                        newShells[toShell]   = newShells[toShell] + 1
                else    fromOcc = occ
                        if  fromOcc - 1 < 0                        println("..");    continue    end
                        newShells[fromShell] = newShells[fromShell] - 1
                        newShells = Base.merge( newShells, Dict( toShell => 1))
                end
                # Add a new configuration
                push!( confList, Configuration( newShells, NoElectrons))
                if  true  println(">> Generate $(Configuration( newShells, NoElectrons))")   end
            end
        end
    end

    nbefore     = length(confList)
    confList = unique(confList)
    ## newConfList = Configuration[]
    ## for  conf  in  confList
    ##     if  conf in newConfList    continue;    else    push!( newConfList,  conf)    end
    ## end
    nafter      = length(confList)
    println(">> Number of generated configurations for $exScheme is: $nbefore (before) and $nafter (after).")

    return( confList )
end



"""
`Basics.generateConfigurationsForExcitationScheme(confs::Array{Configuration,1}, exScheme::Basics.DeExciteTwoElectrons, 
                                                    nMax::Int64, lValues::Array{Int64,1})`  
    ... generates a list of non-relativistic configurations for the given (reference) confs and the excitation scheme. 
        All orbitals in standard order are considered up to maxShell = (n_max, l_max).
"""
function Basics.generateConfigurationsForExcitationScheme(confs::Array{Configuration,1}, exScheme::Basics.DeExciteTwoElectrons, 
                                                            nMax::Int64, lValues::Array{Int64,1})
    newConfList = Basics.generateConfigurationsForExcitationScheme(confs, Basics.DeExciteSingleElectron(), nMax, lValues)
    newConfList = Basics.generateConfigurationsForExcitationScheme(newConfList, Basics.DeExciteSingleElectron(), nMax, lValues)
    newConfList = unique(newConfList)
    return( newConfList )
end



"""
`Basics.generateConfigurationsForExcitationScheme(confs::Array{Configuration,1}, exScheme::Basics.ExciteByCapture, 
        fromShells::Array{Shell,1}, toShells::Array{Shell,1}, intoShells::Array{Shell,1}, noex::Int64)`  
    ... generates a list of non-relativistic configurations for the given (reference) confs with the excitation
        of upto noex electrons fromShell --> toShells and the additional capture of one electron --> intoShells.
        All configuration will therefore contain N+1 electrons.
"""
function Basics.generateConfigurationsForExcitationScheme(confs::Array{Configuration,1}, exScheme::Basics.ExciteByCapture, 
                fromShells::Array{Shell,1}, toShells::Array{Shell,1}, intoShells::Array{Shell,1}, noex::Int64)
    NoElectrons = confs[1].NoElectrons
    confList    = Basics.generateConfigurations(confs, fromShells, toShells, noex)
    shellList   = Basics.extractNonrelativisticShellList(confList)
    #
    # Now add one electron to the intoShells
    newConfList = Configuration[]
    for  conf in confList
        for intoShell in intoShells
            newShells = deepcopy( conf.shells )
            if   haskey(conf.shells, intoShell )   occ = conf.shells[intoShell]
                if  occ   + 1 > 2*(2*intoShell.l + 1)      println("..");    continue    end
                newShells[intoShell] = newShells[intoShell] + 1
            else  
                newShells = Base.merge( newShells, Dict( intoShell => 1))
            end
            # Add a new configuration
            push!( newConfList, Configuration( newShells, NoElectrons+1))
            if  true  println(">> Generate $(Configuration( newShells, NoElectrons+1))")   end
        end
    end
    newConfList = unique(newConfList)
    
    nbefore     = length(confList)
    nafter      = length(newConfList)
    println(">> Number of generated configurations for $exScheme is: $nbefore (before) and $nafter (after).")

    return( newConfList )
end


"""
`Basics.generateFieldCoordinates(mesh::Basics.AbstractMesh)`  
    ... generates a list of field values of proper type due to the specification by the given mesh; 
        this specification determines both, the number and kind of coordinates as well as the type of the array.
        A list::Array{...FieldValue{Type},1} is returned.
"""
function Basics.generateFieldCoordinates(mesh::Basics.AbstractMesh)
    #
    if      typeof(mesh) == Basics.Cartesian2DMesh
        fValues = Basics.Cartesian2DFieldValue{Float64}[]
        xs      = Basics.generateMeshCoordinates(mesh.xMesh)
        ys      = Basics.generateMeshCoordinates(mesh.yMesh)
        for  x in xs, y in ys             push!(fValues, Basics.Cartesian2DFieldValue{Float64}(x,y, 0.))     end
    elseif  typeof(mesh) == Basics.PolarMesh 
        fValues = Basics.PolarFieldValue{Float64}[]
        rhos    = Basics.generateMeshCoordinates(mesh.rhoMesh)
        phis    = Basics.generateMeshCoordinates(mesh.phiMesh)
        for  rho in rhos, phi in phis     push!(fValues, Basics.PolarFieldValue{Float64}(rho,phi, 0.))       end
    elseif  typeof(mesh) == Basics.SphericalMesh
        fValues = Basics.SphericalFieldValue{Float64}[]
        rs      = Basics.generateMeshCoordinates(mesh.rMesh)
        thetas  = Basics.generateMeshCoordinates(mesh.thetaMesh)
        phis    = Basics.generateMeshCoordinates(mesh.phiMesh)
        for  r in rs, theta in thetas, phi in phis     push!(fValues, Basics.SphericalFieldValue{Float64}(r,theta,phi, 0.))    end
    else    error("Unknown mesh type: typeof(mesh) = $(typeof(mesh))")
    end 
    
    return( fValues )
end
        


"""
`Basics.generateLevelWithExtraElectron(newOrbital, symt::LevelSymmetry, level::Level)`  
    ... generates a (new) level with one extra electron in subshell newOrbital.subshell and with overall symmetry symt. The function 
        assumes that all CSF in the basis of level have the same symmetry as the level itself, that the new subshell is not yet part 
        of the basis and that CSF with the total symmetry can be constructed; it terminates with an error if one of this assumptions
        is violated. A newLevel::Level with the same number of CSF and the same representation (eigenvector) as the given level is 
        returned. From a physics viewpoint, newLevel refers to a antisymmetrized product state of [level x orbital(sh, energy)] J^P 
        (symt) and with an total energy = level.energy + newOrbital.energy.
"""
function Basics.generateLevelWithExtraElectron(newOrbital, symt::LevelSymmetry, level::Level)
    basis = level.basis;    newSubshells = copy(basis.subshells);    newCsfs = CsfR[];   J = level.J;   parity = level.parity
    
    push!(newSubshells, newOrbital.subshell)
    Defaults.setDefaults("relativistic subshell list", newSubshells; printout=false)
    newOrbitals = Base.merge( copy(basis.orbitals), Dict( newOrbital.subshell => newOrbital ))

    for  i = 1:length(basis.csfs)
        stateList   = ManyElectron.provideSubshellStates(newOrbital.subshell, 1)
        if  basis.csfs[i].J != J   ||   basis.csfs[i].parity != parity    error("Improper symmetry of CSF.")             end
        if  length(stateList) != 1                                        error("Improper number of subshell states.")   end
        if  AngularMomentum.triangularDelta(J, Basics.subshell_j(newOrbital.subshell), symt.J) != 1    
                                                                        error("Improper coupling of (new) subshell & total J.")   end
        substate   = stateList[1]
        occupation = copy(basis.csfs[i].occupation);    push!(occupation, substate.occ)
        seniorityNr  = copy(basis.csfs[i].seniorityNr);     push!(seniorityNr,  substate.nu)
        subshellJ  = copy(basis.csfs[i].subshellJ);     push!(subshellJ,  AngularJ64(substate.Jsub2//2) )
        subshellX  = copy(basis.csfs[i].subshellX);     push!(subshellX,  symt.J)
        push!(newCsfs, CsfR(true, symt.J, symt.parity, occupation, seniorityNr, subshellJ, subshellX, Subshell[] ) )
    end

    newBasis = Basis(true, basis.NoElectrons+1, newSubshells, newCsfs, copy(basis.coreSubshells), newOrbitals)
    newLevel = Level(symt.J, AngularM64(symt.J), symt.parity, level.index, level.energy+newOrbital.energy, level.relativeOcc, 
                    true, newBasis, level.mc)
    return( newLevel )
end


"""
`Basics.generateLevelWithExtraTwoElectrons(orb1::Orbital, symx::LevelSymmetry, orb2::Orbital, symt::LevelSymmetry, level::Level)`  
    ... generates a (new) level with two extra electrons in subshells orb1.subshell and orb2.subshell and with intermediate
        symmetry symx and overall symmetry symt. The function assumes that all CSF in the basis of level have the 
        same symmetry as the level itself, that the new subshells are not yet part of the basis and that CSF with the 
        intermediate and total symmetries can be constructed; it terminates with an error if one of this assumptions
        is violated. A newLevel::Level with the same number of CSF and the same representation (eigenvector) as the given level is 
        returned. From a physics viewpoint, newLevel refers to a antisymmetrized product state of 
        [ [level x orbital1(sh, energy)] Jx^Px  x orbital1(sh, energy)] J^P  (symt) and with an 
        total energy = level.energy + orb1.energy + orb2.energy.
"""
function Basics.generateLevelWithExtraTwoElectrons(orb1::Orbital, symx::LevelSymmetry, orb2::Orbital, symt::LevelSymmetry, level::Level)
    basis = level.basis;    newSubshells = copy(basis.subshells);    newCsfs = CsfR[];   J = level.J;   parity = level.parity
    
    push!(newSubshells, orb1.subshell);     push!(newSubshells, orb2.subshell)
    Defaults.setDefaults("relativistic subshell list", newSubshells; printout=false)
    newOrbitals = Base.merge( copy(basis.orbitals), Dict( orb1.subshell => orb1, orb2.subshell => orb2 ))

    for  i = 1:length(basis.csfs)
        stateList1   = ManyElectron.provideSubshellStates(orb1.subshell, 1)
        stateList2   = ManyElectron.provideSubshellStates(orb2.subshell, 1)
        if  basis.csfs[i].J != J     ||   basis.csfs[i].parity != parity    error("Improper symmetry of CSF.")             end
        if  length(stateList1) != length(stateList2) != 1                   error("Improper number of subshell states.")   end
        if  AngularMomentum.triangularDelta(J,      Basics.subshell_j(orb1.subshell), symx.J) != 1    
                                                            error("Improper coupling of (new) subshell & intermediate J.")  end
        if  AngularMomentum.triangularDelta(symx.J, Basics.subshell_j(orb2.subshell), symt.J) != 1    
                                                            error("Improper coupling of (new) subshell & total J.")         end
        substate1   = stateList1[1];    substate2  = stateList2[1]
        occupation  = copy(basis.csfs[i].occupation);    push!(occupation, substate1.occ);    push!(occupation, substate2.occ)
        seniorityNr = copy(basis.csfs[i].seniorityNr);   push!(seniorityNr,substate1.nu);     push!(seniorityNr,substate2.nu)
        subshellJ   = copy(basis.csfs[i].subshellJ);     push!(subshellJ,  AngularJ64(substate1.Jsub2//2) )
                                                            push!(subshellJ,  AngularJ64(substate2.Jsub2//2) )
        subshellX   = copy(basis.csfs[i].subshellX);     push!(subshellX,  symx.J);           push!(subshellX,  symt.J)
        push!(newCsfs, CsfR(true, symt.J, symt.parity, occupation, seniorityNr, subshellJ, subshellX, Subshell[] ) )
    end

    newBasis = Basis(true, basis.NoElectrons+2, newSubshells, newCsfs, copy(basis.coreSubshells), newOrbitals)
    newLevel = Level(symt.J, AngularM64(symt.J), symt.parity, level.index, level.energy+orb1.energy+orb2.energy, level.relativeOcc, 
                    true, newBasis, level.mc)
    return( newLevel )
end


"""
`Basics.generateLevelWithExtraSubshell(sh::Subshell, level::Level)`  
    ... generates a (new) level with one extra subshell sh but with the same overall symmetry as before. The function assumes that the 
        new subshell is not yet part of the basis; it terminates with an error if one of this assumptions is violated. 
        A newLevel::Level with the same symmetry, number of CSF, energy and the same representation (eigenvector) as the given 
        level is returned. 
"""
function Basics.generateLevelWithExtraSubshell(sh::Subshell, level::Level)
    basis = level.basis;    newSubshells = copy(basis.subshells);    newCsfs = CsfR[];   J = level.J;   parity = level.parity
    newOrbitals = copy(basis.orbitals)

    push!(newSubshells, sh)
    newOrbitals = Base.merge( newOrbitals, Dict( sh => Orbital(sh, -10000.) ))
    Defaults.setDefaults("relativistic subshell list", newSubshells; printout=false)

    for  i = 1:length(basis.csfs)
        stateList   = ManyElectron.provideSubshellStates(sh, 0)
        if  length(stateList) != 1                                        error("Improper number of subshell states.")   end
        substate   = stateList[1]
        occupation = copy(basis.csfs[i].occupation);    push!(occupation, substate.occ)
        seniorityNr  = copy(basis.csfs[i].seniorityNr);     push!(seniorityNr,  substate.nu)
        subshellJ  = copy(basis.csfs[i].subshellJ);     push!(subshellJ,  AngularJ64(substate.Jsub2//2) )
        subshellX  = copy(basis.csfs[i].subshellX);     push!(subshellX,  subshellX[end])
        push!(newCsfs, CsfR(true, basis.csfs[i].J, basis.csfs[i].parity, occupation, seniorityNr, subshellJ, subshellX, Subshell[] ) )
    end

    newBasis = Basis(true, basis.NoElectrons, newSubshells, newCsfs, copy(basis.coreSubshells), newOrbitals)
    newLevel = Level(level.J, level.M, level.parity, level.index, level.energy, level.relativeOcc, true, newBasis, level.mc)
    return( newLevel )
end


"""
`Basics.generateLevelWithExtraSubshells(subshells::Array{Subshell,1}, level::Level)`  
    ... generates a (new) level with one or several extra subshells but with the same overall symmetry as before. 
        The function assumes that the new subshell is not yet part of the basis. A newLevel::Level with the same symmetry, 
        number of CSF, energy and the same representation (eigenvector) as the given level is returned. 
"""
function Basics.generateLevelWithExtraSubshells(subshells::Array{Subshell,1}, level::Level)
    newLevel = deepcopy(level)
    for  sh in subshells    newLevel = Basics.generateLevelWithExtraSubshell(sh::Subshell, newLevel::Level)    end
    return( newLevel )
end


"""
`Basics.generateLevelWithExtraSubshell(sh::Subshell, level::Level, ns::Int64)`  
    ... generates a (new) level with one extra subshell sh but with the same overall symmetry as before. The function assumes that the 
        new subshell is not yet part of the basis and terminates if this assumptions is not fulfilled. 
        It builts the subshell sh at position ns and moves all other subshells further.
        A newLevel::Level with the same symmetry, number of CSF, energy and the same representation (eigenvector) as the given 
        level is returned. 
"""
function Basics.generateLevelWithExtraSubshell(sh::Subshell, level::Level, ns::Int64)
    basis = level.basis;    newCsfs = CsfR[];   J = level.J;   parity = level.parity
    newOrbitals = copy(basis.orbitals)
    
    # Build the new subshell list and add the associated orbital
    newSubshells = basis.subshells[1:ns-1];   push!(newSubshells, sh);    append!(newSubshells, basis.subshells[ns:end])
    newOrbitals = Base.merge( newOrbitals, Dict( sh => Orbital(sh, -10000.) ))
    Defaults.setDefaults("relativistic subshell list", newSubshells; printout=false)

    for  i = 1:length(basis.csfs)
        stateList   = ManyElectron.provideSubshellStates(sh, 0)
        if    length(stateList) != 1        error("Improper number of subshell states.")   
        else                                substate = stateList[1]
        end
        
        occupation  = copy(basis.csfs[i].occupation[1:ns-1]);  push!(occupation, substate.occ);    append!(occupation, basis.csfs[i].occupation[ns:end]) 
        seniorityNr = copy(basis.csfs[i].seniorityNr[1:ns-1]); push!(seniorityNr,  substate.nu);   append!(seniorityNr, basis.csfs[i].seniorityNr[ns:end]) 
        subshellJ   = copy(basis.csfs[i].subshellJ[1:ns-1]);   push!(subshellJ,  AngularJ64(substate.Jsub2//2) );  append!(subshellJ, basis.csfs[i].subshellJ[ns:end]) 
        if  ns  > 1    subshellX   = copy(basis.csfs[i].subshellX[1:ns-1]);    push!(subshellX,  subshellX[ns-1])
        else           subshellX   = AngularJ64[AngularJ64(0)]                           
        end
        append!(subshellX, basis.csfs[i].subshellX[ns:end])
        push!(newCsfs, CsfR(true, basis.csfs[i].J, basis.csfs[i].parity, occupation, seniorityNr, subshellJ, subshellX, Subshell[] ) )
    end

    newBasis = Basis(true, basis.NoElectrons, newSubshells, newCsfs, copy(basis.coreSubshells), newOrbitals)
    newLevel = Level(level.J, level.M, level.parity, level.index, level.energy, level.relativeOcc, true, newBasis, level.mc)
    return( newLevel )
end


"""
`Basics.generateLevelWithSymmetryReducedBasis(level::Level, subshells::Array{Subshell,1})`  
    ... generates a level with a new basis and representation that only includes the CSF with the same symmetry as the level 
        itself. It also ensures that the new basis is based on the given subshells so that it can be later used with levels
        which have one or several subshells more (at it often appears in ionization processes). Both, the underlying CSF basis 
        and the eigenvectors are adopted accordingly. The procedures only supports the 'addition' of subshells but it stops 
        if (1) the individual sequence differs or if the basis of level contains more subshells than the given subshell-list.
        A (new) level::levelR is returned.
"""
function Basics.generateLevelWithSymmetryReducedBasis(level::Level, subshells::Array{Subshell,1})
    testBasis = level.basis;    newCsfs = CsfR[];    newMc = Float64[]

    # Decide of whether the basis need to be extended by some subshells
    nt = length(testBasis.subshells)
    if      testBasis.subshells == subshells      basis = testBasis;    nLevel = level
    elseif  nt >= length(subshells)               error("stop a")
    elseif  nt + 1 == length(subshells)  &&  testBasis.subshells == subshells[1:nt]
            nLevel = Basics.generateLevelWithExtraSubshell(subshells[end], level)
            basis  = nLevel.basis
    elseif  testBasis.subshells == subshells[1:nt]
        nLevel  = level
        for  ns = nt+1:length(subshells)
            nLevel = Basics.generateLevelWithExtraSubshell(subshells[ns], nLevel)
        end
        basis  = nLevel.basis
    else
        nLevel  = level
        for  (ns, sh) in enumerate(subshells)
            if  sh in testBasis.subshells    continue   
            else      nLevel = Basics.generateLevelWithExtraSubshell(sh, nLevel, ns)
            end
        end
        basis  = nLevel.basis
    end
        
    
    for  i = 1:length(basis.csfs)
        if  basis.csfs[i].J == nLevel.J  &&   basis.csfs[i].parity == nLevel.parity        &&    abs(nLevel.mc[i]) > 1.0e-7
            push!(newCsfs, deepcopy(basis.csfs[i]));    push!(newMc, deepcopy(nLevel.mc[i]))
        end
    end
    newBasis = Basis(true, basis.NoElectrons, deepcopy(basis.subshells), newCsfs, deepcopy(basis.coreSubshells), deepcopy(basis.orbitals))
    newLevel = Level(nLevel.J, nLevel.M, nLevel.parity, nLevel.index, nLevel.energy, nLevel.relativeOcc, true, newBasis, newMc)
    return( newLevel )
end


"""
`Basics.generateMeshCoordinates(mesh::Basics.AbstractMesh)`  
    ... generates a list of coordinates for the given  (1D)mesh; a list::Array{Float64,1} is returned.
        Note: No 'weights' for integration are returned in the present form ... though this could be done readily.
"""
## function Basics.generateMeshCoordinates(mesh::Basics.GLegenreMesh)
function Basics.generateMeshCoordinates(mesh::Basics.AbstractMesh)
    #
    ## using  QuadGK
    if      typeof(mesh) == Basics.GLegenreMesh
        # The weights of the GL coordinates are determined here ... but not returned in the present version.
        wax  = QuadGK.gauss(mesh.NoZeros);    t = wax[1];     wt = wax[2]        
        fac1 = 0.5 * ( mesh.b - mesh.a );   fac2 = 0.5 * ( mesh.b + mesh.a )
        for  j = 1:mesh.NoZeros
            t[j] = fac1 * t[j] + fac2;    wt[j] = fac1 * wt[j]
        end
        coords = t
    elseif  typeof(mesh) == Basics.LinearMesh 
        coords = Float64[];    x0 = (mesh.b - mesh.a) / (mesh.NoPoints-1);   x = mesh.a
        for  i = 1:mesh.NoPoints    x = x + x0;    push!(coords, x)     end
    else    error("Unknown mesh type for generating coordinates: typeof(mesh) = $(typeof(mesh))")
    end 
    
    return( coords )
end


"""
`Basics.generateOrbitalsForBasis(basis::Basis, frozenShells::Array{Shell,1}, priorBasis::Basis, startOrbitals::Dict{Subshell, Orbital})`  
    ... generates a dict of (relativistic) orbitals as specificed by subshells of the basis. These orbitals are taken from priorBasis
        if contained in frozen-shells, and from spectrum (start orbitals) otherwises. An error message is issued if some
        requested orbital is not found. A orbList::Dict{Subshell, Orbital} is returned.
"""
function Basics.generateOrbitalsForBasis(basis::Basis, frozenShells::Array{Shell,1}, priorBasis::Basis, 
                                            startOrbitals::Dict{Subshell, Orbital})
    orbitals = Dict{Subshell, Orbital}()
    # ITEM 24, FIXED 01-Sep-2026.  `subsh in frozenShells` compared a Subshell(n,kappa) against an
    # Array{Shell,1}: there is no == method for that pair, so the fallback === returned false ALWAYS and both
    # frozen branches were dead code.  The consequence looked benign, because the next branch takes the same
    # orbital from priorBasis anyway -- but the error below could never fire, so a shell declared frozen and
    # MISSING from the prior basis was given a hydrogenic start orbital in silence, which at high Z (and with
    # item 23 still open, in the BARE nuclear potential) is not a small error.  The test now compares on (n,l),
    # which is how the driver itself builds AsfSettings.frozenSubshells a few lines below.
    isFrozen(subsh) = any(sh -> sh.n == subsh.n  &&  Basics.subshell_l(subsh) == sh.l, frozenShells)
    for  subsh in basis.subshells
        if  isFrozen(subsh)    &&   haskey(priorBasis.orbitals, subsh)
            orbitals = Base.merge( orbitals, Dict( subsh => priorBasis.orbitals[subsh]) )
        elseif  isFrozen(subsh)           error("Frozen orbital $subsh not found in prior basis.")
        elseif  haskey(priorBasis.orbitals, subsh)
            println(">> Start orbital $subsh is taken from prior basis")
            orbitals = Base.merge( orbitals, Dict( subsh => priorBasis.orbitals[subsh]) )
        else
            println(">> Start orbital $subsh is taken from hydrogenic orbitals")
            orbitals = Base.merge( orbitals, Dict( subsh => startOrbitals[subsh]) )
        end
    end
    return( orbitals )
end



"""
`Basics.generateOrbitalSuperposition(a::Orbital, b::Orbital, cx::Float64, grid::Radial.Grid)`  
    ... generates a superposition  a + cx * b of two given orbitals; the function just takes the linear combination
        of the large and small components as well as the energy. The function assumes that both orbitals are defined on the same
        grid. A re-normalized newOrbital::Orbital is returned for which the derivatives are not defined. 
        This function is used to accelerate the convergence (hopefully).
"""
function Basics.generateOrbitalSuperposition(a::Orbital, b::Orbital, cx::Float64, grid::Radial.Grid)
    if  a.subshell != b.subshell  ||  !(a.isBound)  ||  !(b.isBound)   error("stop a")     end
    mtp  = max(size(a.P, 1), size(b.P, 1));     newP = zeros(mtp);  newQ = zeros(mtp);   newEnergy = (a.energy + cx * b.energy) / (1 + cx)
    mtpa = size(a.P, 1);    newP[1:mtpa] = a.P[1:mtpa];    newQ[1:mtpa] = a.Q[1:mtpa]
    mtpb = size(b.P, 1);    newP[1:mtpb] = newP[1:mtpb] + cx * b.P[1:mtpb]    
                            newQ[1:mtpb] = newQ[1:mtpb] + cx * b.Q[1:mtpb]
    newOrbital = Orbital(a.subshell, a.isBound, a.useStandardGrid, newEnergy, newP, newQ, Float64[], Float64[], a.grid)
    norm       = RadialIntegrals.overlap(newOrbital, newOrbital, grid)
    newP       = newP / sqrt(norm);      newQ = newQ / sqrt(norm)
    newOrbital = Orbital(a.subshell, a.isBound, a.useStandardGrid, newEnergy, newP, newQ, Float64[], Float64[], a.grid)

    return( newOrbital )
end


"""
`Basics.generateShellList(confs::Array{Configuration,1}, nMax::Int64, lValues::Array{Int64,1})`  
    ... generates a list of (non-relativistic) shells in standard order that contains all shells from the given list of
        configurations as well as the shells up to the principal quantum number nMax and the list of orbital quantum
        numbers lValues. A shellList::Array{Shell,1} is returned.
"""
function Basics.generateShellList(confs::Array{Configuration,1}, nMax::Int64, lValues::Array{Int64,1})
    shellList = Shell[]
    # Add all shells from the given configurations
    for  conf in confs
        for (k,v) in conf.shells    push!( shellList, k)    end
    end
    # Add all shells for the given quantum numbers
    for  n = 1:nMax
        for  ll in  lValues  
            if  n >= ll + 1     push!( shellList, Shell(n,ll))      end    
        end
    end
    
    # Now bring the shells in standard order
    newShellList = Shell[]
    for  n = 1:nMax
        for  l = 0:30  
            if  n >= l + 1  && Shell(n,l) in shellList    push!( newShellList, Shell(n,l))     end    
        end
    end
    println(">> From configurations generated shell list $newShellList ")
    
    return( newShellList )
end


"""
`Basics.generateShellList(nMin::Int64, nMax::Int64, lValues::Array{Int64,1})`  
    ... generates a list of (non-relativistic) shells in standard order that contains all shells with principal quantum
        numbers from nMin .. nMaxthe and orbital angular momenta from lValues. A shellList::Array{Shell,1} is returned.
"""
function Basics.generateShellList(nMin::Int64, nMax::Int64, lValues::Array{Int64,1})
    shellList = Shell[]
    # Add all shells for the given quantum numbers
    for  n = nMin:nMax
        for  ll in  lValues  
            if  n >= ll + 1     push!( shellList, Shell(n,ll))      end    
        end
    end
    
    # Now bring the shells in standard order
    newShellList = Shell[]
    for  n = 1:nMax
        for  l = 0:30  
            if  n >= l + 1  && Shell(n,l) in shellList    push!( newShellList, Shell(n,l))     end    
        end
    end
    println(">> Generated shell list $newShellList ")
    
    return( newShellList )
end


"""
`Basics.generateShellList(nMin::Int64, nMax::Int64, lMax::Int64)`  
    ... generates a list of (non-relativistic) shells in standard order that contains all shells with principal quantum
        numbers from nMin .. nMaxthe and orbital angular momenta l <= lMax. A shellList::Array{Shell,1} is returned.
"""
function Basics.generateShellList(nMin::Int64, nMax::Int64, lMax::Int64)
    shellList = Shell[]
    # Add all shells for the given quantum numbers
    for  n = nMin:nMax
        for  ll = 0:lMax 
            if  n >= ll + 1     push!( shellList, Shell(n,ll))      end    
        end
    end
    
    # Now bring the shells in standard order
    newShellList = Shell[]
    for  n = 1:nMax
        for  l = 0:30  
            if  n >= l + 1  && Shell(n,l) in shellList    push!( newShellList, Shell(n,l))     end    
        end
    end
    ## println(">> Generated shell list $newShellList ")
    
    return( newShellList )
end


"""
`Basics.generateShellList(nMin::Int64, nMax::Int64, symMax::String)`  
    ... generates a list of (non-relativistic) shells in standard order that contains all shells with principal quantum
        numbers from nMin .. nMaxthe and orbital angular momenta l <= l(symMax), and where symMax is the symmetry character
        string ("s", "p", ...,"z") of the corresponding orbital angular momentum. A shellList::Array{Shell,1} is returned.
"""
function Basics.generateShellList(nMin::Int64, nMax::Int64, symMax::String)
    lMax = Basics.shellNotation(symMax)
    shellList = Shell[]
    # Add all shells for the given quantum numbers
    for  n = nMin:nMax
        for  ll = 0:lMax 
            if  n >= ll + 1     push!( shellList, Shell(n,ll))      end    
        end
    end
    
    # Now bring the shells in standard order
    newShellList = Shell[]
    for  n = 1:nMax
        for  l = 0:30  
            if  n >= l + 1  && Shell(n,l) in shellList    push!( newShellList, Shell(n,l))     end    
        end
    end
    println(">> Generated shell list $newShellList ")
    
    return( newShellList )
end


"""
`Basics.generateSubshellList(shells::Array{Shell,1})`  
    ... generates a list of relativistic subshells from the given shell list, and by keeping the same order.
        A subshellList::Array{Subshell,1} is returned.
"""
function Basics.generateSubshellList(shells::Array{Shell,1})
    subshellList = Subshell[]
    for  sh in shells
        if  sh.l == 0   push!(subshellList, Subshell(sh.n, -1))
        else            push!(subshellList, Subshell(sh.n, sh.l))
                        push!(subshellList, Subshell(sh.n, -sh.l - 1))
        end
    end
        
    return( subshellList )
end


"""
`Basics.generateSpectrumLorentzian(xIntensities::Array{Tuple{Float64,Float64},1}, widths::Float64; 
                                    energyShift::Float64=0., resolution::Int64=500)`  
    ... to generate the Lorentzian spectrum for the given intensities (x-position, y-intensity) and the (constant) width.
        A tuple of an xList::Array{Float64,1} (positions) and yList::Array{Float64,1} (valiues) is returned that can immediately
        be plotted. The energy shift is provided already in 'current' units and is simply added to all energies.
"""
function Basics.generateSpectrumLorentzian(xIntensities::Array{Tuple{Float64,Float64},1}, widths::Float64; 
                                            energyShift::Float64=0., resolution::Int64=500)
    function normalLorentzian(en::Float64, ga::Float64)
        gaover2 = ga /2.0;  wa = gaover2 / (en^2 + gaover2^2);      wa = wa / pi
    end
    #
    # Convert energies and determine their range 
    wc = Defaults.convertUnits("energy: from atomic", 1.0);  xwidths = widths * wc
    energies = Float64[];     intensities = Float64[]
    for xInt in xIntensities    push!(energies, xInt[1]*wc + energyShift);    push!(intensities, xInt[2])   end
    xMin = minimum(energies) - 3.0xwidths
    xMax = maximum(energies) + 3.0xwidths
    res = 1 / resolution;   nn = resolution * Int64( floor(xMax-xMin) + 1.0 )
    x = zeros(nn);    y = zeros(nn)
    for  n = 1:nn     x[n] = xMin +n*res    end
    for  n = 1:nn 
        for ie = 1:length(energies)
            en = energies[ie];   int = intensities[ie]
            y[n] = y[n] + int * normalLorentzian(x[n] - en, xwidths)
        end
    end
    
    return(x,y)
end


