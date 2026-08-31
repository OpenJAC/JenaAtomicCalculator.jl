
"""
`module  JAC.SelfConsistent`  
	... a submodel of JAC that contains all structs and methods to generate self-consistent fields of different 
	    kind and complexity.
"""
module SelfConsistent

using  Printf, LinearAlgebra, ..AngularMomentum, ..Basics, ..Bsplines, ..Defaults, ..Hamiltonian, ..InteractionStrength, ..ManyElectron, ..Nuclear, ..Radial,
       ..RadialIntegrals, ..SpinAngular

# The coefficient carriers, aliased so that the signatures below stay readable. The KIND parameter is the point:
# a rank-0 coefficient multiplies an ORDINARY one-electron matrix element and a rank-k one a REDUCED one, and
# pairing the wrong two raises a MethodError naming both rather than returning a wrong number (Rule 18).
const Coefficient1p = SpinAngular.Coefficient1p{SpinAngular.OrdinaryKind}
const Coefficient2p = SpinAngular.Coefficient2p{SpinAngular.EffectiveStrengthKind}


# Module-level defaults. Each is set through its own SelfConsistent.set... function below.

# EXPERIMENTAL SWITCH (09-Aug-2026), default false = the behaviour that has always been in place.
# When true, solveOptimizedLevelField builds the Fock matrix as
#      h  +  (1/occ) * V(diagonal CSF pairs)  +  V(off-diagonal CSF pairs)
# instead of  h + (1/occ) * V(all pairs).  This tests whether the off-diagonal (CI-coupling) part should
# be scaled by the generalized occupation at all: it carries weight ~ c_r c_s while occ carries ~ c_r^2,
# so dividing it by occ introduces a c_1/c_2 factor that grows without bound as a correlating CSF's own
# coefficient shrinks -- the suspected mechanism behind the winner-take-all collapse documented in
# solveOptimizedLevelField's KNOWN LIMITATION.  Set with SelfConsistent.setEolUnscaledOffDiagonal(true).
GBL_EOL_UNSCALED_OFFDIAGONAL = false

# Re-orthonormalize the same-kappa orbitals after the damping step. DEFAULT TRUE since 10-Aug-2026;
# the switch is kept only so the two behaviours can still be compared.
#
# Hamiltonian.projectHamiltonian makes each raw eigenvector orthogonal to the already-processed same-kappa
# orbitals, but the damping that follows, mixed = 0.5*old + 0.5*raw, mixes it back with the PREVIOUS
# iteration's vector, which is NOT orthogonal to them -- and nothing restored it. The CSF expansion assumes
# an orthonormal orbital set, so the resulting energies were not legitimate variational numbers.
# Measured on Li 1s^2 2s + 1s^2 3s + 1s^2 3d (three s-orbitals in kappa = -1), converged:
#
#                     <2s|3s>      <1s|2s>      E
#     as it was      -1.128e-03   -5.558e-05   -7.4335291982
#     re-orthonorm.  -1.373e-11   -5.358e-13   -7.4335284248
#
# The energy RISES by 7.7e-07 Ha, which is the honest direction: the non-orthogonal set was giving a
# slightly-too-low number. The whole approved test suite is blind to this (44/44 either way), which is why
# it went unnoticed -- TestFrames.testMethod_OrbitalOrthonormality now asserts it directly.
# The fix needs no new code: orthonormalizeSameKappa was written for exactly this and had never
# been called from anywhere.
GBL_SCF_REORTHONORMALIZE = true

# Anderson depth for the AVERAGE-LEVEL field, separate from the mean-field one above because the iterate is
# different: there it is the screening potential, here the orbitals themselves.  0 = the plain damped
# iteration exactly as before.
#
# ON since 17-Aug-2026, at the same depth 2 the mean-field driver uses.  It reaches the SAME solution --
# Ar 3s^2 3p^6 agrees to 4.3e-9 once accuracyScf is tight enough to converge at all -- and is 1.4x to 1.9x
# faster, the gain GROWING with the accuracy demanded (1.47x at 1e-6, 1.79x at 1e-9, 1.89x at 1e-12).
# On Be 4-config it is also the more STABLE of the two: across accuracyScf = 1e-6, 1e-9, 1e-12 it drifts by
# 1e-5 where the plain iteration swings 6.5e-5 NON-MONOTONICALLY, and it gets there in 211 s against 1093 s.
#
# It could not be switched on until 18aaf5f.  Anderson perturbs the orbital tails just enough to move the
# old mtp cut, which made TestFrames.testMethod_OrbitalOrthonormality report 6.5e-08 where the plain
# iteration gave 6.8e-10 -- an artefact of the truncated integral, not of the orbitals, which were
# orthonormal to 1e-17 throughout.  With the tails kept, both give ~1e-16.
#
# STILL OPEN, and NOT fixed by any of this: plain AL does not converge for a multi-configuration basis with
# near-degenerate CSFs (Be 4-config), and the default accuracyScf = 1e-6 hides it by stopping early.
# Tightening the tolerance is therefore not a general cure -- right for Ar-like cases, worse for Be.
GBL_AL_ANDERSON_DEPTH = 2

# Anderson depth for the mean-field (DFS/HS) SCF.  0 = the plain iteration exactly as before; a positive
# value routes performSCF to SelfConsistent.solveMeanFieldBasisAnderson, which reaches the SAME self-consistent
# solution in fewer iterations.  DEFAULT 0 so that nothing changes unless it is asked for.
#
# The plain iteration converges linearly, the residual shrinking by a constant factor r per step; measured
# 12-Aug-2026, r = 0.44 (Ar 1s^2..3p^6), 0.57 (Ne 1s^2 2s^2 2p^6), 0.69 (Fe [Ar] 3d^6 4s^2).  Anderson
# mixing builds the next screening potential from a least-squares combination of the last few iterates and
# their residuals, cancelling the slowest-decaying error rather than waiting for it to decay:
#
#                        plain        depth 2      agreement of the converged orbital energies
#     Ne  2s^2 2p^6      28 it        13 it        7.9e-07
#     Ar  3s^2 3p^6      23 it        14 it        6.9e-09
#     Fe  3d^6 4s^2      45 it        18 it        6.5e-07
#     Ne+ 1s hole        17 it        11 it        1.2e-07
#     Ar+ 1s hole        14 it        11 it        1.2e-07   (added 30-Aug-2026)
#
# THE 30-Aug ROW WAS TAKEN WITH maxIterationsScf = 200, AND THAT IS NOT A DETAIL. AsfSettings() defaults the
# ceiling to 24, and the plain Ne iteration needs 28: a first attempt at this measurement reported "24 it" for
# Ne, which is the CEILING and not a convergence count, beside a NOT-CONVERGED flag easy to read as a property
# of the case rather than of the settings. Ne was therefore re-run as a CONTROL for the new row and reproduces
# 28 exactly. Any iteration count taken here without raising the ceiling first measures the default.
# Ne's depth-2 count came out 12 against the 13 recorded above; the grid differs (Basics.recommendedGrid was
# used in 2026-08-30, and the 12-Aug measurement does not record what it used), so this is most likely the
# convergence path rather than a change. It is left standing rather than smoothed over.
# Ar+ is the SMALLEST gain of the five (1.27x, against 2.15x Ne, 1.64x Ar, 2.50x Fe, 1.55x Ne+), consistent
# with Anderson helping most where the plain iteration converges slowest -- a 1s-hole system already converges
# quickly.
#
# Depth 2 is the measured optimum; 3 is nearly equal, and LARGER IS WORSE (Ne: 24 it at depth 5, 36 at 12),
# the usual ill-conditioning of a long Anderson history.  Note that depth 0 in the Anderson driver itself is a
# JACOBI sweep and does NOT converge in 60 iterations -- the Gauss-Seidel ordering of the original driver is
# what makes the plain iteration viable at all.
#
# STANDARD SINCE 12-Aug-2026 (was 0 = the plain iteration when this was first added).  Setting it to 0
# restores the old path exactly, which is how the two were compared.  Making it the default is a deliberate
# editorial act: both iterations reach the same self-consistent solution, but only to within accuracyScf,
# so results can move by ~1e-6.  What that costs was measured by regenerating every approved reference with
# it on: 27 of 29 came out BITWISE IDENTICAL, one moved by 3.0e-09, and the only large apparent change --
# test-Cascade-StepwiseDecay -- was not numerical at all, but two DEGENERATE levels swapping index labels
# (all 606 transition rows the same set).
GBL_SCF_ANDERSON_DEPTH = 2


"""
`SelfConsistent.setScfAndersonDepth(depth::Int64)`
    ... sets the Anderson-mixing depth of the mean-field SCF; 0 restores the plain iteration. Nothing is returned.
"""
function setScfAndersonDepth(depth::Int64)
    global GBL_SCF_ANDERSON_DEPTH = depth
    return( nothing )
end


"""
`SelfConsistent.setEolUnscaledOffDiagonal(flag::Bool)`
    ... sets the experimental switch that keeps the off-diagonal CSF-pair contributions OUT of the
        (1/occ) scaling in the EOL Fock matrix; nothing is returned.
"""
function setEolUnscaledOffDiagonal(flag::Bool)
    global GBL_EOL_UNSCALED_OFFDIAGONAL = flag
    return( nothing )
end


"""
`SelfConsistent.setScfReorthonormalize(flag::Bool)`  
    ... sets the experimental switch that re-orthonormalizes same-kappa orbitals after the SCF damping
        step; nothing is returned.
"""
function setScfReorthonormalize(flag::Bool)
    global GBL_SCF_REORTHONORMALIZE = flag
    return( nothing )
end


"""
`SelfConsistent.initializeBasis(configs::Array{Configuration,1}, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives,
                                settings::AsfSettings; levelSymmetries::Array{LevelSymmetry,1}=LevelSymmetry[], printout::Bool=false)` 
    ... Initialized a many-electron basis from the given list of configurations, the nuclear model as well as ASF settings.
        It assumes that a proper set of primitives::Primitives has been initialized before. The initial set of orbitals in this
        basis is determined by the settings::AsfSettings.  A basis::Basis is returned.
"""
function initializeBasis(configs::Array{Configuration,1}, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives, 
                         settings::AsfSettings; levelSymmetries::Array{LevelSymmetry,1}=LevelSymmetry[], printout::Bool=true)
    NoElectrons = configs[1].NoElectrons;   subshells = Subshell[];   coreSubshells = Subshell[];     csfs = CsfR[] 
    orbitals    = Dict{Subshell, Orbital}()
    
    # Perform some simple tests: Number of electrons must be equal in all configurations
    for  conf in configs   if  conf.NoElectrons != NoElectrons    error("stop a")   end     end
    
    # Generate a full set of relativistic CSF from the given configurations and collect the associated level symmetries
    relconfList = ConfigurationR[]
    for  conf in configs
        wa = Basics.generateConfigurations(Basics.RelativisticConfigurations(), conf)
        append!( relconfList, wa)
    end
    if  printout    for  i = 1:length(relconfList)    println(">>> include ", relconfList[i])    end   end
    subshells = Basics.generateSubshellList(relconfList)
    Defaults.setDefaults("relativistic subshell list", subshells; printout=printout)

    # Generate the relativistic CSF's for the given subshell list
    csfList = CsfR[]
    for  relconf in relconfList
        newCsfs = Basics.generateCsfRs(relconf, subshells)
        append!( csfList, newCsfs)
    end
    
    # Select CSF with requested symmetry if needed
    if  length(levelSymmetries) == 0
        csfs = csfList          # Take all relativistic CSF into account
    else
        for  csf in csfList
            if  LevelSymmetry(csf.J, csf.parity)  in  levelSymmetries   push!(csfs, csf)    end
        end
    end

    # Determine the number of electrons and the list of coreSubshells
    for  (k,sh)  in  enumerate(subshells)
        mocc = Basics.subshell_2j(sh) + 1;    is_filled = true
        for  csf in csfList
            if  csf.occupation[k] != mocc     is_filled = false;           break   end
        end
        if   is_filled    push!( coreSubshells, subshells[k])      else    break   end
    end
        
    # Check that the radial grid is able to represent these subshells at all, BEFORE any orbital is generated.
    # This sits here rather than inside Bsplines.generateOrbitalsHydrogenic so that it applies whichever way the
    # orbitals are seeded -- StartFromPrevious inherits a grid just as much as StartFromHydrogenic does. The
    # occupations are handed over so that each subshell is tested at the charge it actually feels; at the bare
    # charge the check rejects the valence orbital of any heavy near-neutral system, whose box must be matched
    # to a screened orbital some thirty times more extended than the bare-Z one.
    occupations = Dict{Shell,Int64}()
    for  conf in configs
        for  (sh, occ)  in conf.shells
            if  occ > 0     occupations[sh] = max( get(occupations, sh, 0), occ )    end
        end
    end
    Bsplines.checkGridRepresentation(subshells, nuclearModel.Z, primitives; occupations = occupations,
                                     accuracy = settings.gridAccuracy, stopper = settings.gridStopper)

    # Initialize the orbitals
    if  typeof(settings.startScfFrom) == StartFromHydrogenic
        if  printout   println("> Start SCF process with hydrogenic orbitals.")   end
        # Generate start orbitals for the SCF field by using B-splines
        orbitals  = Bsplines.generateOrbitalsHydrogenic(subshells, nuclearModel, primitives; printout=printout)
    elseif  typeof(settings.startScfFrom) == StartFromThomasFermi
        if  printout   println("> Start SCF process with orbitals in a Thomas-Fermi potential.")   end
        # The nucleus screened by a statistical model of the electron cloud.  Unlike every self-consistent
        # field this needs no density, so it is available before any orbital exists -- which is the point of
        # a start potential.  Bsplines.generateOrbitals then does what it does for any other potential.
        tfPot     = Basics.add( Nuclear.nuclearPotential(nuclearModel, primitives.grid),
                                Basics.computePotential(Basics.ThomasFermiField(), primitives.grid,
                                                        nuclearModel.Z, NoElectrons) )
        orbitals  = Bsplines.generateOrbitals(subshells, tfPot, nuclearModel, primitives; printout=printout)
    elseif  typeof(settings.startScfFrom) == StartFromPrevious
        if  printout   println("> Start SCF process from given list of orbitals.energy")    end
        # Take what the given set provides and fall back to a hydrogenic orbital for anything it does not.
        # THE FALLBACK NEVER RAN BEFORE (fixed 13-Aug-2026): it called HydrogenicIon.radialOrbital(subsh, ...,
        # grid) where neither `subsh` nor `grid` exists in this method -- the loop variable is `sh` and only
        # `primitives` is passed -- and radialOrbital takes a Shell rather than a Subshell in any case, so the
        # branch could only ever have thrown.  It went unnoticed because StartFromPrevious had exactly one
        # caller, which always supplied a complete set.  Warm-starting one cascade block from another does
        # not: consecutive blocks differ in which subshells are occupied.  The missing ones are now generated
        # by the same B-spline routine the StartFromHydrogenic branch above uses, in ONE call.
        orbitals = Dict{Subshell, Orbital}()
        missingSubshells = Subshell[]
        for  sh in subshells
            if  haskey(settings.startScfFrom.orbitals, sh)   orbitals[sh] = settings.startScfFrom.orbitals[sh]
            else                                             push!(missingSubshells, sh)
            end
        end
        if  length(missingSubshells) > 0
            if  printout   println("> Start orbitals do not cover $missingSubshells; these are taken hydrogenic.")   end
            hydrogenic = Bsplines.generateOrbitalsHydrogenic(missingSubshells, nuclearModel, primitives; printout=false)
            for  sh in missingSubshells   orbitals[sh] = hydrogenic[sh]   end
        end
    else  error("stop b")
    end
    
    basis = Basis(true, NoElectrons, subshells, csfs, coreSubshells, orbitals)
    return( basis )
end


"""
`SelfConsistent.checkOneElectronSelfInteraction(configs::Array{Configuration,1}, scField::Basics.AbstractScField)`
    ... answers whether a ONE-ELECTRON system is about to be solved in a mean-field potential built from its own
        density, printing the explanation when it is; a `value::Bool` is returned, `true` meaning that the field must
        be replaced. `performSCF` acts on that answer by switching to `Basics.NuclearField()`, so the condition is
        CORRECTED rather than merely reported -- see the note there for why the correction is exact and not a guess.

        WHY THIS IS WORTH SAYING OUT LOUD RATHER THAN A NOTE SOMEWHERE. A mean field such as DFS is built from the total
        electron density, so with a SINGLE electron it contains that electron's own charge: the electron is repelled
        by itself. The bound orbital and the continuum orbital then solve DIFFERENT one-body operators, which is not
        a small error and does not look like one. Measured 30-Aug-2026 on H(1s) -> continuum: the Coulomb and
        Babushkin photoionization cross sections, which are EQUAL for exact hydrogenic wavefunctions, come out in
        the ratio 1.238, 1.299 and 1.294 at omega/I = 1.5, 2 and 3 under `DFSField()`, and 1.000000 at every one of
        them under `NuclearField()` -- unchanged across boxes of 30 to 150 a.u., meshes of 0.05 to 0.01, and the
        pure-sine, pure-Coulomb and Ong-Russek normalisations.

        That 24-30 % discrepancy was carried on the priority list for weeks as a suspected defect of the CONTINUUM
        machinery, and the three suspects it named were all innocent. `AsfSettings()` defaults to `DFSField()`, so
        the default is the trap; the remark that would have prevented it existed only inside a test for another
        module, where no user would meet it.
"""
function checkOneElectronSelfInteraction(configs::Array{Configuration,1}, scField::Basics.AbstractScField)
    NoElectrons = length(configs) == 0  ?  0  :  sum( values(configs[1].shells) )
    if  NoElectrons != 1   ||   typeof(scField) == Basics.NuclearField   return( false )   end

    sa = "SelfConsistent.performSCF(): a ONE-ELECTRON system was requested in $(nameof(typeof(scField))), " *
         "which is built from the electron's OWN density, so the electron would be repelled by itself; " *
         "solving in Basics.NuclearField() instead, where the one-electron problem is exact."
    printstyled("\n>> " * sa * "\n", color=:light_red)
    printstyled(">> This would not be a small error: on H(1s) it puts the two photoionization gauges, which are " *
                "EQUAL\n>> for exact hydrogenic wavefunctions, in the ratio 1.24-1.30 rather than 1.000000.\n",
                color=:light_red)
    Defaults.warn(AddWarning(), sa)

    return( true )
end


"""
`SelfConsistent.checkScFieldIsSupported(scField::Basics.AbstractScField)`
    ... verifies that the given field is one that performSCF can actually iterate, and raises an explanatory
        error if it is not. Several members of Basics.AbstractScField are POTENTIALS rather than fields: they
        answer Basics.providesPotential but not Basics.providesScfDriver, and one of them is not
        self-consistent even in principle. Checked HERE, before a grid, a set of primitives and a many-electron
        basis have been built, so that an unsupported choice costs nothing and says what to do instead.
        Nothing is returned.
"""
function checkScFieldIsSupported(scField::Basics.AbstractScField)
    Basics.providesScfDriver(scField)   &&   return( nothing )
    sa = "\n\nSelfConsistent.performSCF(): $(nameof(typeof(scField))) is not a self-consistent field that " *
         "this driver can iterate.\n"
    if      typeof(scField) == Basics.ThomasFermiField
        sa = sa * ">>> It is not self-consistent even in principle: it needs the nuclear charge and the "     *
                  "electron number only, and no\n    density at all. That is exactly what makes it useful "   *
                  "as a STARTING potential -- use ManyElectron.StartFromThomasFermi,\n    or ask for the "    *
                  "potential itself with Basics.computePotential.\n"
    elseif  typeof(scField) in [Basics.AaDFSField, Basics.AaHSField]
        sa = sa * ">>> It is an AVERAGE-ATOM potential for finite temperature: its Basics.computePotential "  *
                  "takes a chemical potential\n    and a temperature, so it belongs to the plasma line "      *
                  "rather than to this bound-state SCF. Use Plasma.AverageAtomScheme\n    (see "              *
                  "examples/example-Ja.jl).\n"
    else
        sa = sa * ">>> It is a screened POTENTIAL, not a field: it owns a Basics.computePotential method but " *
                  "no SCF driver of its own.\n    Use Basics.computePotential(...) to obtain it, or choose "  *
                  "one of the fields listed below.\n"
    end
    sa = sa * ">>> The fields performSCF iterates are:  " * join(Basics.scfDriverFields(), ", ") *
              "  (NuclearField with a hydrogenic start only).\n"
    error(sa)
end


"""
`SelfConsistent.performSCF(configs::Array{Configuration,1}, nm::Nuclear.Model,
                           settings::AsfSettings; levelSymmetries::Array{LevelSymmetry,1}=LevelSymmetry[], printout::Bool=true)`
    ... performs a SCF computation for which NO grid is given, so that the radial box is derived from the
        configurations themselves by Basics.recommendedGrid; a multiplet::Multiplet is returned.

        This exists because choosing the box is the step most likely to be got wrong, and getting it wrong
        does not look like a grid problem: the record attributes four separate "bugs" -- an E3 rate 1000x too
        small, a Zeeman kappa <= -3 failure, a MultipolePolarizibility defect and a Breit sign flip -- to a
        box that did not match the orbitals, and each was first blamed on the angular machinery.  The derived
        box beats JAC's hand-chosen default grid for every system it has been measured on (see
        Basics.recommendedGrid), by 2.9e-3 Ha for argon and 2.3e-2 Ha for Ti+.

        ONE CASE NEEDS THE GRID GIVEN BY HAND.  The estimate cannot tell a spectroscopic Rydberg shell from a
        CORRELATION shell of the same n and l, and reads both as diffuse: a beryllium basis carrying 3s and 3d
        for correlation is given 67 a.u., which is right for a real 1s^2 2s 3s state and far too generous for
        a correlation orbital that contracts onto the valence region.  Where the high-n shells are there to
        correlate rather than to be occupied, pass a grid, or pass `rbox` to Basics.recommendedGrid.
"""
function performSCF(configs::Array{Configuration,1}, nm::Nuclear.Model,
                    settings::AsfSettings; levelSymmetries::Array{LevelSymmetry,1}=LevelSymmetry[], printout::Bool=true)
    grid = Basics.recommendedGrid(configs, nm, printout=printout)

    return( SelfConsistent.performSCF(configs, nm, grid, settings,
                                      levelSymmetries=levelSymmetries, printout=printout) )
end


"""
`SelfConsistent.performSCF(configs::Array{Configuration,1}, nm::Nuclear.Model, grid::Radial.Grid,
                           settings::AsfSettings; levelSymmetries::Array{LevelSymmetry,1}=LevelSymmetry[], printout::Bool=false)`
    ... Performs a SCF computation for the given list of configurations, the nuclear model as well as ASF settings.
        If explicit levelSymmetries are given, only these symmetries are considered. Internally, a proper set of primitives::Primitives 
        is initialized and used in the computations. The generated SCF field is controlled by the settings::AsfSettings.  
        A multiplet::Multiplet is returned.
"""
function performSCF(configs::Array{Configuration,1}, nm::Nuclear.Model, grid::Radial.Grid, 
                    settings::AsfSettings; levelSymmetries::Array{LevelSymmetry,1}=LevelSymmetry[], printout::Bool=true)
    
    SelfConsistent.checkScFieldIsSupported(settings.scField)
    # A system with ONE electron has no electron-electron interaction to average, so every self-consistent field
    # degenerates to the bare nuclear one -- except that a mean field built from the density would include this
    # electron's own charge and repel it from itself. The substitution is therefore not a fallback or a guess: it
    # is the exact answer, and the hydrogenic orbitals it starts from ARE the converged ones, which is why
    # startScfFrom is set with it (a caller's StartFromPrevious would otherwise reach the :hydrogenicStartOnly
    # branch below and raise). Found 31-Aug-2026 in four example branches, three of them dated -- and in three of
    # the four the one-electron multiplet is built inside cascade or MultiPhotonTransition machinery, where no
    # keyword the caller could set reaches this point. That is what makes correcting it here right and warning
    # about it wrong: the person who would have to act on the warning has no way to.
    if  SelfConsistent.checkOneElectronSelfInteraction(configs, settings.scField)
        settings = AsfSettings(settings; scField = Basics.NuclearField(), startScfFrom = StartFromHydrogenic())
    end

    # Generate primitives and initialize the many-electron basis
    Defaults.setDefaults("standard grid", grid)
    primitives = Bsplines.generatePrimitives(grid)    
    basis      = SelfConsistent.initializeBasis(configs, nm, primitives, settings; levelSymmetries, printout)
    
    # Solve a self-consistent field for this basis
    scfProc = Basics.scfProcedure(settings.scField)
    if   scfProc == :meanFieldIteration
        # GBL_SCF_ANDERSON_DEPTH = 0 keeps the plain iteration; a positive depth reaches the SAME
        # self-consistent solution in fewer iterations (see the note at the switch).
        if  GBL_SCF_ANDERSON_DEPTH > 0
            basis = SelfConsistent.solveMeanFieldBasisAnderson(basis, nm, primitives, settings; printout=printout,
                                                            andersonDepth=GBL_SCF_ANDERSON_DEPTH)
        else
            basis = SelfConsistent.solveMeanFieldBasis(basis, nm, primitives, settings; printout=printout)
        end 
    elseif   scfProc == :hydrogenicStartOnly  &&  settings.startScfFrom == StartFromHydrogenic()
        # Return the basis as already generated.
    elseif   scfProc == :averageLevel
        basis     = SelfConsistent.solveAverageLevelField(basis, nm, primitives, settings; printout=printout)
    elseif   scfProc == :optimizedLevel
        # EOL is done by ORBITAL ROTATION.  The older solveOptimizedLevelField, which folds the off-diagonal
        # CSF-pair terms into the same (1/occ)-scaled Fock matrix, converges to a DEGENERATE stationary point
        # whenever two near-degenerate CSFs compete for one correlation channel: the correlating weight runs
        # to zero, the correlation orbital then no longer enters the energy, and its gradient vanishes for a
        # trivial reason.  Measured on Be 1s^2 2s^2 + 1s^2 2p^2 it lands 19.4 mHa ABOVE the average-level
        # field on the very level it is asked to optimize, with the 2p weight collapsed from 0.25 to 0.0001.
        # Rotating the orbitals instead escapes that point and reaches 5.3 mHa BELOW AL.  See example-Ao.jl.
        # The rotation is a LOCAL optimizer, so it starts from an average-level basis rather than from the
        # initial guess -- that is how it was validated, and a hydrogenic start has no reason to lie in its
        # basin.  Both solvers return a complete, correctly (kink-aware) diagonalized multiplet, so return it
        # directly; falling through would re-diagonalize with the bare, non-kink-aware Hamiltonian.performCI.
        alSettings = AsfSettings(settings; scField = Basics.ALField())
        basis      = SelfConsistent.solveAverageLevelField(basis, nm, primitives, alSettings; printout=printout)
        return( SelfConsistent.solveOptimizedLevelFieldByRotation(basis, nm, primitives, settings; printout=printout) )
    else  error("stop a")
    end

    # Now that the orbitals are final, check that no symmetry has converged onto the wrong state. This
    # catches what Bsplines.checkGridRepresentation cannot: that check tests HYDROGENIC orbitals -- since
    # 17-Aug-2026 at the SCREENED charge rather than the bare one, which is what lets a heavy near-neutral
    # system through at all -- and a hydrogenic test can never see a symmetry that has converged onto a
    # different state, Ge II 4f on a 614 a.u. box being the case that motivated it. Note the EOLField branch
    # above returns early and is therefore not covered here.
    Bsplines.checkOrbitalConsistency(basis.orbitals, grid; stopper = settings.gridStopper)
    # ... and that the box was in fact large enough for them.  This is the only test made on the CONVERGED
    # orbitals rather than on a hydrogenic stand-in, so it is the one that says whether the estimate the
    # grid was built from turned out to be right.
    (boxOk, boxReport) = Bsplines.checkOrbitalBox(basis.orbitals, grid; stopper = false)
    if  printout   println(">> Radial box: " * boxReport * (boxOk ? "  -- adequate." : "  -- NOT adequate."))   end

    # Setup and diagonalize the Hamiltonian matrix; assign mixing coefficients
    if   scfProc == :averageLevel
        mp = Hamiltonian.performCIKinkAware(basis, nm, grid, settings, printout=printout)
    else
        mp = Hamiltonian.performCI(basis, nm, grid, settings, printout=printout)
    end

    return( mp )
end


"""
`SelfConsistent.performSCF(basis::Basis, nm::Nuclear.Model, grid::Radial.Grid,
                           settings::AsfSettings; levelSymmetries::Array{LevelSymmetry,1}=LevelSymmetry[], printout::Bool=false)`
    ... Performs a SCF computation for the given list of configurations, the nuclear model as well as ASF settings.
        If explicit levelSymmetries are given, only these symmetries are considered. Internally, a proper set of primitives::Primitives 
        is initialized and used in the computations. The generated SCF field is controlled by the settings::AsfSettings.  
        A multiplet::Multiplet is returned.
"""
function performSCF(basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, 
                    settings::AsfSettings; levelSymmetries::Array{LevelSymmetry,1}=LevelSymmetry[], printout::Bool=false)
    
    SelfConsistent.checkScFieldIsSupported(settings.scField)

    # Generate primitives
    primitives = Bsplines.generatePrimitives(grid)    
    
    # Solve a self-consistent field for this basis
    scfProc = Basics.scfProcedure(settings.scField)
    if   scfProc == :meanFieldIteration
        # GBL_SCF_ANDERSON_DEPTH = 0 keeps the plain iteration; a positive depth reaches the SAME
        # self-consistent solution in fewer iterations (see the note at the switch).
        if  GBL_SCF_ANDERSON_DEPTH > 0
            basis = SelfConsistent.solveMeanFieldBasisAnderson(basis, nm, primitives, settings; printout=printout,
                                                            andersonDepth=GBL_SCF_ANDERSON_DEPTH)
        else
            basis = SelfConsistent.solveMeanFieldBasis(basis, nm, primitives, settings; printout=printout)
        end 
    elseif   scfProc == :hydrogenicStartOnly  &&  settings.startScfFrom == StartFromHydrogenic()
        # Return the basis as already generated.
    elseif   scfProc == :averageLevel
        basis     = SelfConsistent.solveAverageLevelField(basis, nm, primitives, settings; printout=printout)
    elseif   scfProc == :optimizedLevel
        # See the note in the other performSCF overload just above: EOL is done by orbital rotation, started
        # from an average-level basis, and returns a complete, correctly (kink-aware) diagonalized multiplet.
        alSettings = AsfSettings(settings; scField = Basics.ALField())
        basis      = SelfConsistent.solveAverageLevelField(basis, nm, primitives, alSettings; printout=printout)
        return( SelfConsistent.solveOptimizedLevelFieldByRotation(basis, nm, primitives, settings; printout=printout) )
    else  error("stop a")
    end

    # Now that the orbitals are final, check that no symmetry has converged onto the wrong state. This
    # catches what Bsplines.checkGridRepresentation cannot: that check tests HYDROGENIC orbitals -- since
    # 17-Aug-2026 at the SCREENED charge rather than the bare one, which is what lets a heavy near-neutral
    # system through at all -- and a hydrogenic test can never see a symmetry that has converged onto a
    # different state, Ge II 4f on a 614 a.u. box being the case that motivated it. Note the EOLField branch
    # above returns early and is therefore not covered here.
    Bsplines.checkOrbitalConsistency(basis.orbitals, grid; stopper = settings.gridStopper)
    # ... and that the box was in fact large enough for them.  This is the only test made on the CONVERGED
    # orbitals rather than on a hydrogenic stand-in, so it is the one that says whether the estimate the
    # grid was built from turned out to be right.
    (boxOk, boxReport) = Bsplines.checkOrbitalBox(basis.orbitals, grid; stopper = false)
    if  printout   println(">> Radial box: " * boxReport * (boxOk ? "  -- adequate." : "  -- NOT adequate."))   end

    # Setup and diagonalize the Hamiltonian matrix; assign mixing coefficients
    if   scfProc == :averageLevel
        mp = Hamiltonian.performCIKinkAware(basis, nm, grid, settings, printout=printout)
    else
        mp = Hamiltonian.performCI(basis, nm, grid, settings, printout=printout)
    end

    return( mp )
end


# The four self-consistent-field paths and the machinery they share. Split out of this file on 27-Aug-2026;
# the dispatcher (performSCF), the basis setup and the module-level defaults stay here.

include("module-SelfConsistent-inc-orbitals.jl")
include("module-SelfConsistent-inc-averagelevel.jl")
include("module-SelfConsistent-inc-optimizedlevel.jl")
include("module-SelfConsistent-inc-meanfield.jl")
include("module-SelfConsistent-inc-averageatom.jl")

end # module
