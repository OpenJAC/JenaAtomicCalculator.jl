
println("Fj) Cascade.HollowIonScheme: formation and decay of a hollow carbon ion.")

using JLD2, Printf
#
setDefaults("method: continuum, Galerkin")              ## see the note in example-Fg.jl on asymptotic Coulomb
setDefaults("method: normalization, Alok")
setDefaults("unit: energy", "eV")

grid = Radial.Grid(Radial.Grid(true); rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 20.0)


# WRITTEN 08-Aug-2026, Stage 5 of the cascade-scheme series.  The previous example-Fj.jl held the Si^-
# photo-ionization scenario rescued from the old example-Fb.jl; example-Fd.jl covers PhotoIonizationScheme
# properly since, so nothing was kept.
#
# WHAT A HOLLOW ION IS.  An ion whose INNER shells are empty while outer ones are occupied -- the state left
# behind when a highly charged ion captures several electrons into high nl shells, as happens at a surface or
# in a dense plasma.  Such a configuration is far above the ground state of its own charge state and decays by
# a cascade of two competing processes: radiative transitions that fill the inner vacancy, and autoionization
# (Auger) that ejects an electron.  Cascade.HollowIonScheme follows both.
#
# The reference case below is the simplest possible: a BARE carbon nucleus that captures two electrons into
# n = 2, giving 2s^2 / 2s2p / 2p^2 with an entirely empty K shell.  Those states lie above the 1s threshold of
# the one-electron ion, so both channels are open -- 2p -> 1s radiative decay and Auger emission to 1s.
#
# SIX DEFECTS were fixed in module-Cascade-inc-hollow-ion.jl to get this far; the scheme had never run.  The
# first five are of the classes already familiar from the other schemes:
#   1. five stray "@show" statements, three of them appended to live code lines;
#   2. SelfConsistent.performSCF returns a Multiplet, but generateBlocks used it as a Basis (the same stale
#      assumption as in module-Cascade-inc-impact-excitation.jl);
#   3. Cascade.DecayData(linesR, linesA) -- that type no longer exists, Cascade.Data{T} replaced it;
#   4. the results were stored ONLY under "hollow-ion line data:", while every simulation reads
#      "cascade data:" via Cascade.reviewData.  No Cascade.Simulation could consume a hollow-ion cascade at
#      all; branch d now does.
#   5. the configuration generator emitted a 0-electron configuration -- the fully stripped ion, produced by
#      the repeated RemoveElectrons -- which is not a cascade block and which aborted Hamiltonian.performCI
#      inside Basics.merge with a bare error("stop a").
#   6. Basics.computePotentialDFS, called to build the potential for the Auger continuum, DOES NOT EXIST.  The
#      working pattern is that of the dielectronic-recombination scheme,
#      Basics.computePotential(comp.asfSettings.scField, grid, level), which also honours a non-DFS choice.
#      The same dead call still sits in module-Cascade-inc-stepwise-decay.jl:174 and in
#      module-PhotoDoubleIonization.jl:356 and is left there for those files' own tasks.
#
# AND ONE CHANGE OF SUBSTANCE: decayShells now means what its docstring says.  It is documented as "shells into
# which the electrons can decay", but generateConfigurationsForHollowIons only ever REMOVED electrons from
# those shells (two RemoveElectrons loops) and never filled them.  Starting from an empty K shell, no
# configuration with an occupied 1s was generated at all -- so the dominant 2p -> 1s radiative channel and the
# Auger decay into the ground state of the next ion were both unreachable, and the scheme computed only the
# redistribution among the captured shells.  The generator now also produces the configurations reached by
# moving an electron from an into-shell down into a decay shell.  For the case below that takes the list from
# three configurations to nine:
#       2 electrons :  2s^2   2s2p   2p^2        (hollow)
#                      1s2s   1s2p   1s^2        (reached by radiative decay into the K shell)
#       1 electron  :  2s     2p     1s          (reached by autoionization)
# A TRAP inside that fix, worth remembering: intoShells and decayShells normally overlap, and
# Basics.generateConfigurations then also returns configurations with an electron ADDED rather than moved.
# Filtering those out on Configuration.NoElectrons DOES NOT WORK -- that field is not kept in step with the
# shells by these generators, which is why the routine's own final loop recomputes it.  Counting the shell
# occupations directly is what reduced 26 generated configurations to the 9 above.


if  true
    # Last visit:      08-Aug-2026
    # Last successful: 08-Aug-2026 ... 41 s.  14 steps; 38 PhotoEmission and 10 AutoIonization lines, written
    #                  to zzz-cascade-hollow-ion-computations-<date>.jld.  Both channels are open, which is the
    #                  whole point of the scheme and was not the case before the decayShells fix above.
    #                  Physics checks: the radiative steps come out at 354 - 365 eV against the C VI Lyman-alpha
    #                  of 367 eV, and the Auger electron at 269 eV.  Dated on that agreement together with the
    #                  channel-union check of branch c.
    #
    # Branch a: THE REFERENCE CASCADE -- a bare carbon nucleus captures two electrons into n = 2 and decays by
    #   both channels.  Radiative and Auger steps are requested together, which is the point of the scheme.
    #   The lines are written to a .jld that branch d reads.
    setDefaults("print summary: open", "zzz-Cascade-Fj-reference.sum")

    name   = "Hollow carbon: K-shell empty, two electrons captured into n = 2"
    scheme = Cascade.HollowIonScheme([Radiative(), Auger()], [E1], 2, [Shell("2s"), Shell("2p")],
                                     [Shell("1s"), Shell("2s"), Shell("2p")])
    wa     = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(6.), grid=grid,
                                 approach=Cascade.AverageSCA(), scheme=scheme,
                                 initialConfigs=[Configuration("1s^0")] )
    println(wa)
    wb = perform(wa; output=true, outputToFile=true)
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      08-Aug-2026
    # Last successful: 08-Aug-2026 ... 26 s for both runs.
    #                      decayShells            radiative lines    Auger lines
    #                      2s, 2p (no 1s)               11                0
    #                      1s, 2s, 2p                   38               10
    #                  Without the 1s the K shell can never be filled, so there is no 2p -> 1s photon and no
    #                  autoionization whatever: the cascade collapses to the redistribution among the captured
    #                  n = 2 shells.  This is the measurement behind the generator change described above.
    #
    # Branch b: WHAT decayShells CONTROLS -- the same capture, run once WITHOUT the 1s among the decay shells
    #   and once with it.  Without 1s no configuration with an occupied K shell can be generated, so the
    #   cascade is reduced to the redistribution among the captured n = 2 shells: no 2p -> 1s photon and no
    #   Auger decay at all.  This is the branch that shows why the generator had to be extended, and it is
    #   also the practical warning: decayShells must list the shells the electrons are to decay INTO, not only
    #   the ones they start in.
    setDefaults("print summary: open", "zzz-Cascade-Fj-decayshells.sum")

    for  (sa, dec)  in  [("without 1s", [Shell("2s"), Shell("2p")]),
                         ("with 1s",    [Shell("1s"), Shell("2s"), Shell("2p")])]
        println("\n>>> decayShells $sa")
        scheme = Cascade.HollowIonScheme([Radiative(), Auger()], [E1], 2, [Shell("2s"), Shell("2p")], dec)
        wa = Cascade.Computation(Cascade.Computation(); name="Hollow C, decayShells $sa",
                                 nuclearModel=Nuclear.Model(6.), grid=grid,
                                 approach=Cascade.AverageSCA(), scheme=scheme,
                                 initialConfigs=[Configuration("1s^0")] )
        wb = perform(wa; output=true, outputToFile=false)
        for  x  in  wb["cascade data:"]
            println(">>>   ", eltype(x.lines), " : ", length(x.lines), " lines")
        end
    end
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      08-Aug-2026
    # Last successful: 08-Aug-2026 ... 48 s for the three runs.
    #                      processes           radiative lines    Auger lines
    #                      [Radiative()]             38                0
    #                      [Auger()]                  0               10
    #                      both                      38               10
    #                  Requesting both gives EXACTLY the union of the two, 38 and 10 in each case -- so the
    #                  step generation neither double-counts a block pair nor drops one when both processes
    #                  are active.  That is a free internal check and it passes exactly.
    #
    # Branch c: THE TWO CHANNELS SEPARATELY -- the same cascade with processes = [Radiative()], [Auger()] and
    #   both.  scheme.processes is one of the few process lists in the Cascade module that is genuinely read
    #   (determineSteps loops over it), so this is a real switch and not decoration.  Requesting both must give
    #   exactly the union of the two, which is the cheapest available check that the step generation is not
    #   double-counting or dropping pairs.
    setDefaults("print summary: open", "zzz-Cascade-Fj-channels.sum")

    for  (sa, procs)  in  [("radiative only", [Radiative()]), ("Auger only", [Auger()]),
                           ("both",           [Radiative(), Auger()])]
        println("\n>>> processes: $sa")
        scheme = Cascade.HollowIonScheme(procs, [E1], 2, [Shell("2s"), Shell("2p")],
                                         [Shell("1s"), Shell("2s"), Shell("2p")])
        wa = Cascade.Computation(Cascade.Computation(); name="Hollow C, $sa", nuclearModel=Nuclear.Model(6.),
                                 grid=grid, approach=Cascade.AverageSCA(), scheme=scheme,
                                 initialConfigs=[Configuration("1s^0")] )
        wb = perform(wa; output=true, outputToFile=false)
        for  x  in  wb["cascade data:"]
            println(">>>   ", eltype(x.lines), " : ", length(x.lines), " lines")
        end
    end
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      08-Aug-2026 ... runs, and the machinery is verified: Cascade.reviewData now finds the
    #                  lines (it could not before defect 4 was fixed), sorts 21 levels across all charge states,
    #                  propagates, and returns a normalised distribution summing to 1.00000.
    #                  NOT "Last successful", and the reason is a TRAP worth knowing.  With
    #                  initialOccupations = [(1, 1.0)] the run reported
    #                      No. electrons   Lev-No   J^P      Energy [eV]      Rel. occ.
    #                             1           1    3/2 -    -1.224164e+02     1.00000
    #                  i.e. level 1 is a ONE-electron level (2p of C^5+), not a hollow two-electron state at
    #                  all.  The cascade therefore started at a dead end and the resulting "100% at one
    #                  electron" says nothing: we began there.  Level numbers in initialOccupations refer to the
    #                  SORTED, CASCADE-WIDE list that spans every charge state, so the obvious choice of 1 can
    #                  land on a FINAL state rather than an initial one.  Read the "Initial level occupation"
    #                  table that the simulation prints and pick a level number whose electron count matches the
    #                  hollow configuration before quoting any distribution from this branch.
    #
    # Last visit:      23-Aug-2026
    # Last successful: 23-Aug-2026 ... 4.6 s for the three simulations, after branch a (60 s warm) has written
    #                  the data file.  Ion distribution, per starting level:
    #                      level  4  (0+, -213.6 eV)   2 electrons 1.72e-02   1 electron 9.83e-01
    #                      level  9  (1+, -218.9 eV)   2 electrons 1.00e+00   1 electron 5.62e-05
    #                      level 13  (0-, -223.8 eV)   2 electrons 6.14e-02   1 electron 9.39e-01
    #                  Each column sums to 1 to five digits, which is the internal check.
    #
    # Branch d: THE SIMULATION -- the ion distribution that the hollow ion decays into.  This is the physical
    #   question a hollow-ion cascade is usually asked: starting from the doubly-captured state, how is the
    #   final charge distributed between the one- and two-electron ions?  Radiative decay keeps the electron
    #   number, Auger decay reduces it, so the answer is set by the competition that branch c separates.
    #   Requires branch a to have run.
    #
    #   A DEFECT FOUND ON 23-Aug-2026, WHICH IS WHY THIS BRANCH WAS NEVER DATED.  It used to start the
    #   propagation from Cascade.IonDistribution([(1, 1.0)], ...), i.e. from level 1 -- and the level numbers
    #   of a simulation run over the SORTED list of all levels of all charge states together, where level 1 is
    #   a ONE-ELECTRON level at -122 eV.  The cascade therefore began where it should have ended, propagated
    #   nothing, and reported "1 electron: 1.0" -- a result that is perfectly self-consistent and answers no
    #   question at all.  The same trap as any level index that means something different from what the caller
    #   assumes.
    #
    #   WHICH LEVELS ARE THE HOLLOW ONES.  Sorted, the 21 generated levels fall out as: 1-3 and 14 with ONE
    #   electron; 4-13 with two, and these are exactly the 2l2l' manifold (2s^2 gives one level, 2s2p four,
    #   2p^2 five); 15-20 the 1s2l states; 21 the 1s^2 ground state.  The doubly-captured hollow ion is the
    #   2l2l' manifold, i.e. levels 4-13.
    #
    #   AND THE ANSWER DEPENDS STRONGLY ON WHICH OF THEM IS POPULATED, which is why three are run rather than
    #   one being chosen.  Level 9 keeps 99.99 % of its population at two electrons -- its Auger channel is
    #   all but shut -- while levels 4 and 13 lose 98 % and 94 % to Auger.  A single representative level
    #   would have made the hollow ion look either almost entirely Auger-decaying or almost entirely
    #   radiative, depending only on which was picked.
    setDefaults("print summary: open", "zzz-Cascade-Fj-simulation.sum")

    # The sibling cascade files -- Fb, Fe and Fi -- test for an empty list here and say what to do.
    # This line indexed it with [end] instead, so a missing cascade file surfaced as
    # `BoundsError: attempt to access 0-element Vector{String} at index [0]`, which names neither the
    # cause nor the remedy. Index [0] is never valid in Julia, so the message was doubly misleading.
    fnList     = sort(filter(f -> startswith(f, "zzz-cascade-hollow-ion-"), readdir()), by = f -> stat(f).mtime)
    if  isempty(fnList)   error("Run branch a of this file first; no zzz-cascade-hollow-ion-* file found in the " *
                                "working directory.")   end
    fn    = fnList[end]
    println(">>> reading the cascade data from  $fn")
    data = [JLD2.load(fn)]
    for  levelNo  in  [4, 9, 13]
        println("\n>>> starting the propagation from level $levelNo of the 2l2l' hollow manifold")
        prop = Cascade.IonDistribution([(levelNo, 1.0)], Configuration[])
        simu = Cascade.Simulation(Cascade.Simulation(); name="Hollow-carbon decay from level $levelNo",
                                  computationData=data, property=prop,
                                  settings=Cascade.SimulationSettings(false, false, 0.) )
        wd = perform(simu; output=true)
    end
    setDefaults("print summary: close", "")
    #
end
