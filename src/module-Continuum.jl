
"""
`module  JAC.Continuum`  
... a submodel of JAC that contains all functions and methods for computing, either in a given potential or with
    regard to a well-defined bound state.

    THE GUARD ORDER, agreed 03-Nov-2026, and a NEW SOLUTION METHOD OR NORMALIZATION CONVENTION IS EXPECTED TO FIT IT.
    This module supplies the free electron to every process in the package and is called at a low level, often from
    inside a cascade or a second-order sum, so a caller never names it and a wrong orbital here is silent and reaches
    every observable downstream.  The guards are therefore LOUD and sit at four fixed stages, each raising with the
    measurement that justifies it rather than with a marker:

    1. BEFORE SOLVING -- can this grid represent an orbital at this energy at all?  `gridConsistency` tests the step
       against the de Broglie wavelength (15 points per oscillation) and `checkNormalizationRadius` tests that the
       point where the normalization is fixed lies outside the bound charge cloud, in the dimensionless Zeff*r.
       Both are A-PRIORI: they examine the grid and the energy, never the orbital that comes out.
    2. AT EVERY DISPATCH on one of the two session globals -- `GBL_CONT_SOLUTION` and `GBL_CONT_NORMALIZATION` --
       the `else` branch calls `unsupportedScheme`, which names the value in force, the path that has no branch for
       it, what that path does implement and how to change the setting.  The three dispatchers do NOT implement the
       same set, so this is per site and not once for the module.
    3. AFTER NORMALIZING -- does the orbital that was actually produced deserve the scaling just applied to it?
       `reportContinuumAmplitude` asks two questions in this order, because they fail independently:
         (a) the FORM -- the rms misfit of the large component against the asymptotic form the convention assumed,
             as a FRACTION of the amplitude, from `asymptoticAmplitude`.  Only this measure sees an orbital that has
             not reached its oscillation inside the box.
         (b) the MAGNITUDE -- `checkContinuumAmplitude`'s interior/asymptotic ratio, which sees a mesh that cannot
             follow the oscillation it did reach.
       ALL FOUR CONVENTIONS PASS THROUGH IT, BUT A MEASURE MAY ONLY GUARD THE ASSUMPTION ITS OWN ROUTINE MAKES.
       Three of them normalize against the ASYMPTOTIC FORM of the wave, so 3(a) is precisely their validity
       criterion.  `normalizeOrbitalCoulombFunctions` instead matches the EXACT Coulomb or Bessel functions, which hold at any
       radius, and is legitimately usable where the asymptotic expansion is not -- so it supplies no residual and
       3(a) is skipped for it.  Uniformity of the MEASURE was tried on 03-Nov-2026 and withdrawn the same day:  on
       Ne-like Kr at 16.7 eV in a 10 a.u. box, eta ~ 24 makes the asymptotic form meaningless and the residual 4.77,
       which refused four suite tests whose numbers were sound.  The consistency that matters is of the PRINCIPLE.
       AND 3(a) ITSELF MUST FOLLOW THE SOLUTION METHOD, not only the potential:  `asymptoticCharge` drops the
       logarithmic Coulomb phase for the two free-wave solvers, which otherwise read a hundred times too large.
    4. EVERY CHECK HAS TWO TIERS -- `@warn` where the result is suspect and `error` where it is unusable -- with the
       thresholds DERIVED from a stated measurement and quoted in the message, never chosen to make a case pass.

    So a new scheme needs: a branch at all three dispatchers (stage 2 names the sites if one is forgotten), and either
    a residual for stage 3(a) or an explicit `NaN` with the reason stated, as `normalizeOrbitalCoulombFunctions` once did.
"""
module Continuum


using  GSL, Printf, SpecialFunctions, DelimitedFiles
using  ..Basics, ..Bsplines, ..Defaults, ..ManyElectron, ..Radial, ..Nuclear


"""
`struct  Continuum.Settings`  ... defines a type for the parameters for computing continuum orbitals.

    + includeExchange         ::Bool    ... True, if the exchange with the bound electrons of the level is to be included,
                                            false otherwise. SEE THE NOTE BELOW: `true` is honoured in exactly ONE place.
    + mtp                     ::Int64   ... No of grid points for which the continuum orbital(s) are to be computed.

    DELIBERATELY POSITIONAL, AND DELIBERATELY OUTSIDE THE `Settings` CONVENTION -- decided 28-Aug-2026. Every other
    `Settings` in the package carries a zero-argument constructor and a keyword copy-constructor, because every other one is
    something a USER assembles and then varies. This one is not: it is an internal argument-holder that a process module
    fills in on its way to `Continuum.generateOrbitalForLevel`, and it is built positionally at all 27 sites that construct it in
    `src/` and `examples/` -- 20 of them the identical `Continuum.Settings(false, nrContinuum)`. Adding the convention here
    would give a user a constructor for a struct they have no decision to make about. Do not "fix" this; it is the decision.

    **`includeExchange = true` IS NOW REAL, AS OF 29-Oct-2026 (item 99), BUT ONLY UNDER `Continuum.generateOrbitalForLevel`
    WITH THE GALERKIN SOLUTION METHOD** (`setDefaults("method: continuum, Galerkin")`). There it adds
    `Bsplines.setupExchangeMatrix` -- the frozen nonlocal Dirac-Fock exchange operator built for item 75/S30, validated to
    1.1e-16 against the whole-atom exchange energy -- against every occupied subshell of `level.basis`, with occupations from
    `Basics.extractMeanOccupation(level.basis)`. **No iteration is needed**: with the core frozen, the exchange operator is
    LINEAR in the unknown continuum orbital, so it simply adds to the Galerkin matrix that is already being solved once; this
    is the same one-shot assembly a frozen-core bound Fock matrix uses, not a new convergence loop. Under any OTHER solution
    method (`ContBessel`, `ContSine`, `AsymptoticCoulomb`, `NonrelativisticCoulomb`) `generateOrbitalForLevel` RAISES rather
    than silently ignoring the flag -- those integrate directly against a local potential function and have no place to put a
    nonlocal term without a genuine fixed-point iteration, which does not exist here (see the module docstring). And
    `generateOrbitalLocalPotential` still raises on `true` unconditionally, for a different reason: it is handed a bare
    `pot::Radial.Potential` with no `level` at all, so it has no bound orbitals to exchange with, Galerkin or not. All 27
    EXISTING construction sites pass `false` and are completely unaffected; the field remains the way a caller OPTS IN.
"""
struct Settings
    includeExchange           ::Bool
    mtp                       ::Int64
end


"""
`Continuum.resolveExchangeDefault(explicit::Union{Nothing,Bool}, scField::Basics.AbstractScField)`
    ... decides whether a continuum orbital built for THIS computation should carry exchange with the bound electrons:
        `explicit`, whenever it is given, always wins; otherwise the decision is made FROM THE SCF FIELD alone, so that the
        question never has to be put to a user who never asked it. `true` for `ALField` and `EOLField` -- the two fields
        whose bound orbitals come from a genuine nonlocal Fock SCF, so the continuum orbital built in `Continuum`'s own
        LOCAL potential is otherwise missing exactly the interaction its own bound orbitals were built with -- and `false`
        for every other `scField`: a local potential (`DFSField`, `HSField`, `KSField`, `CHField`, `ThomasFermiField`,
        `AaDFSField`, `AaHSField`) has no exchange of its own to match, and `NuclearField` has no bound electrons at all.
        A `value::Bool` is returned. Item 99, route 1; this is the one function every process module's `continuumExchange`
        field is resolved through, via `Basics.perform`.
"""
function resolveExchangeDefault(explicit::Union{Nothing,Bool}, scField::Basics.AbstractScField)
    !isnothing(explicit)   &&   return( explicit )

    return( scField isa Basics.ALField  ||  scField isa Basics.EOLField )
end


"""
`Continuum.unsupportedScheme(fn::String, setting::String, scheme, implemented::String, howToSet::String)`
    ... raises because a dispatch on one of the two session-wide continuum choices -- the SOLUTION method
        `Defaults.GBL_CONT_SOLUTION` or the NORMALIZATION convention `Defaults.GBL_CONT_NORMALIZATION` -- reached its `else`
        branch, i.e. the value in force is one that this code path has no branch for.  The message names the value, the path,
        what that path does implement and how to change the setting;  nothing is returned, since it always raises.

        WHY THIS IS A FUNCTION AND NOT SIX `error("stop a")` CALLS, which is what stood here until 03-Nov-2026.  These two
        settings are SESSION GLOBALS that a user sets once and that are then read far below, inside a cascade or a
        second-order sum, so an unhandled value surfaces as a bare `stop a` from a place the caller never named.  Worse, the
        three dispatch sites do NOT implement the same set:  `generateOrbitalLocalPotential` handles `ContSine()` and the two
        `generateOrbitalForLevel` methods do not, so `setDefaults("method: continuum, pure sine")` followed by the level
        route is a reachable trap that said only `stop a`.  AND THE POINT IS THE NEXT SCHEME, not this one:  a scheme added
        to `Basics.AbstractContinuumSolutions` or `Basics.ContinuumNormalization` compiles and runs until it reaches a
        dispatcher, and it must then fail loudly, naming itself and the site that needs the new branch.
"""
function unsupportedScheme(fn::String, setting::String, scheme, implemented::String, howToSet::String)

    error("\n\nContinuum." * fn * "():  STOP -- " * setting * " = " * string(scheme) *
          " has no branch on this code path.\n" *
          "    this path implements   : " * implemented * "\n" *
          "    how to change it       : " * howToSet * "\n" *
          "The two continuum choices are session globals, set once and read far below -- inside a cascade or a second-order\n" *
          "sum -- so this is raised where the dispatch happens rather than passed on as a wrong number.  If you have just\n" *
          "ADDED a scheme to Basics.AbstractContinuumSolutions or Basics.ContinuumNormalization, this is the site that needs\n" *
          "the new branch;  the other dispatchers in this module need one too, and they do not all implement the same set.\n")
end


#  The wording of the six dispatch guards, written once so that they cannot drift apart.  Note that the two
#  generateOrbitalForLevel methods implement a SMALLER set than generateOrbitalLocalPotential:  there is no
#  level-route branch for ContSine, and a caller who selects it gets a named refusal rather than a marker.
const  CONT_IMPLEMENTED_SOLUTION_LEVEL = "ContBessel, AsymptoticCoulomb, NonrelativisticCoulomb, BsplineGalerkin  " *
                                         "(NOT ContSine -- the level route has no branch for it)"
const  CONT_IMPLEMENTED_SOLUTION_POT   = "ContBessel, ContSine, AsymptoticCoulomb, NonrelativisticCoulomb, BsplineGalerkin"
const  CONT_IMPLEMENTED_NORMALIZATION  = "PureSineNorm, CoulombSineNorm, OngRussekNorm, CoulombFunctionNorm"
const  CONT_HOW_TO_SET_SOLUTION        = "setDefaults(\"method: continuum, ...\")  with  spherical Bessel | pure sine | " *
                                         "asymptotic Coulomb | nonrelativistic Coulomb | Galerkin"
const  CONT_HOW_TO_SET_NORMALIZATION   = "setDefaults(\"method: normalization, ...\")  with  pure sine | pure Coulomb | " *
                                         "Ong-Russek | Coulomb functions"


#  THRESHOLDS FOR THE MATCHING-RADIUS CHECK of the Coulomb-function normalization, in RADIANS of inner phase shift.
#  MEASURED 05-Nov-2026 by tools/probeContinuumMatchingRadius.jl over Z = 2, 18, 36, 92 at 0.05, 0.5, 5 and 50 Ha
#  for kappa = -1, 1, -2, each on a box of 10 wavelengths with 40 mesh points per oscillation.
#    33 SOUND orbitals spread between 6.0e-05 and 5.4e-04 rad -- so the worst sound case sits a factor 3.7 below
#  the warn level, and all 33 reproduce the pre-guard normalization constant and phase EXACTLY.
#    THE SPREAD IS CONVERTIBLE INTO WHAT A USER CARES ABOUT, which is why these numbers can be defended:  across
#  the whole set the normalization constant's own variation over the probed radii is 2.4x the phase spread for a
#  sound orbital and 3.5x for a bad one.  So `warn` is the point at which the normalization becomes uncertain at
#  about the half-percent level, and `fatal` the point at which it exceeds ~30 % and is no longer a property of
#  the solution at all.
const  CONT_MATCH_WARN                 = 2.0e-3
const  CONT_MATCH_FATAL                = 1.0e-1
#  AND THE PHASE ALONE IS HALF-BLIND, which the mutation pass of 05-Nov-2026 found and no amount of looking at
#  sound orbitals could have:  the inner phase shift is extracted as the ANGLE of a complex number built from `P`
#  and `Pprime`, so multiplying BOTH by the same factor leaves it EXACTLY unchanged.  A 30 % multiplicative
#  distortion of the outer half of an orbital, and 2 % point-to-point noise, each left the phase spread identical
#  to the tenth digit and passed as sound -- while the NORMALIZATION CONSTANT extracted at those same radii moved
#  by 4.8 % and 4.0 % against 1.5e-04 for the sound orbital, a separation of more than a factor 300.
#    So the normalization constant's own spread across the probed radii is the second, independent measure, and it
#  is the one that answers the question a user actually asks:  is the normalization well defined?  Sound orbitals
#  run to 1.2e-03;  the ruined ones to between 2.0e-02 and 1.8e-01.
const  CONT_SCALE_WARN                 = 5.0e-3
const  CONT_SCALE_FATAL                = 1.0e-1

#  THE MESH DENSITY BELOW WHICH MOVING THE MATCHING RADIUS IS NOT A REPAIR, in mesh points COUNTED over one
#  oscillation of the de Broglie wavelength at the matching point.
#    WHY THIS IS A BRACKET AND NOT A CROSSING, and the measurement is worth reading before anyone tightens it.
#  MEASURED 06-Nov-2026 over four systems -- Ar and Kr at 5 Ha, U at 50 Ha, Ar at 0.5 Ha, each kappa = -1 -- with
#  the box held fixed and grid.hp swept over thirteen values, comparing the normalization constant against the
#  finest mesh of its own system.  THE ERROR IS NOT MONOTONE IN THE MESH DENSITY:

#      Ar  5 Ha    counted  16     19     21     21     22     23     24     27     30     36     48
#                  error   +11.5  +9.4   +1.5   -1.9   +1.7   -0.9   +0.1   +0.1   -0.1   -0.5   -0.3  %
#      U  50 Ha    counted  15     20     20     21     22     23     25     29     35     46     62
#                  error    +0.0  +3.4   +1.8   -1.2   +2.3   -0.3   +0.0   +0.1   +0.0   +0.0   -0.1  %

#  TWO MESHES AT THE SAME COUNTED DENSITY, 21, DIFFER BY 3.3 % (+1.5 and -1.9), AND U AT 15 POINTS IS EXACT WHILE U
#  AT 20 IS 3.4 % OUT.  So the counted density does NOT determine the error, and no single value of it can be a
#  calibrated threshold -- which is the same non-monotone signature `Continuum.checkNormalizationRadius` already
#  records for the normalization RADIUS (0.085 then 3.78 on a FINER grid), and has the same cause:  the matching
#  point landing at an arbitrary phase of the oscillation.
#    WHAT THE SAME SWEEP DOES ESTABLISH, strongly, is that THE VERDICT IS THE RELIABLE SIGNAL AND THE DENSITY IS
#  NOT.  Every mesh judged `:sound` reproduced its system's converged constant to better than 0.5 %, and most to
#  better than 0.2 %, in all four systems;  the meshes judged `:poor` scatter by up to 12 %.  That is what the
#  phase-settling measure was built to do and it does it.
#    HENCE WHAT THIS CONSTANT IS FOR, narrowly:  it decides only whether a MOVE of the matching radius may be
#  taken, never whether an orbital is good.  Two measured points bracket it.  A move at 30 counted points returned
#  -0.09 %, i.e. the repair worked;  the move that the dielectronic-recombination test case took at about 18
#  counted points turned a +6.0 % error into a -5.6 % one, i.e. it relocated the error rather than removing it and
#  reported success while doing so.  28 lies inside that bracket.  The crossing is not resolved more finely, and
#  the scatter above shows why it cannot be from this quantity alone.
const  CONT_MESH_SOUND                 = 28.0


#  HOW MANY OSCILLATIONS OF THE OUTGOING WAVE MUST FIT IN THE BOX, for a convention that FITS AN ASYMPTOTIC SINE.
#  Both numbers were measured 06-Nov-2026, and the experiment is `tools/probeSineBoxWaves.jl`:  a bare POINT proton,
#  so the exact solution is the Coulomb function and `CoulombFunctionNorm` is correct at any radius and serves as
#  the CONTROL.  The statistic is the norm of the NORMALIZED large component over a radial window both boxes
#  contain, which is scale-free -- the Galerkin eigenvector's arbitrary scale cancels, and after normalization the
#  same physical orbital must come out of any adequate box.  `hp` is equal in the two boxes of a pair, so the
#  grids NEST and the window carries identical mesh points and identical weights.

#      series A, box fixed at 12 a.u. and the ENERGY varied (eta = 5.0 ... 0.33):
#          waves  0.54   0.76   1.05   1.43   1.91   2.67   3.82   5.73
#          error  14.6%  61.2%  23.1%  11.3%   5.6%   9.1%   3.4%   0.8%      control 0.5 ... 3.6 %
#      series B, ENERGY fixed (eta = 1.818) and only the box varied -- the confound removed:
#          waves  0.54   0.71   0.97   1.23   1.76   2.46   3.51   5.27
#          error  74.6%  57.7%  26.2%  10.2%  19.3%   2.6%  13.2%   6.0%      control 7e-8 above 8 a.u.

#  THE CONTROL IS WHAT MAKES THIS A MEASUREMENT.  In series B the default convention reproduces its own value to
#  7e-8 across a factor seven in box size, so the box, the mesh and the solve are all innocent and the error
#  belongs to the sine fit alone.  A first version of this experiment compared the normalization CONSTANT between
#  boxes and was discarded:  N absorbs the eigenvector's arbitrary scale, so that comparison has no meaning, and
#  its "control" drifted by factors of two.

#  SO THE FATAL LEVEL IS ONE FULL OSCILLATION, and it is where the two series agree:  below one wave the error is
#  15-75 %, which squares into a cross section wrong by a factor of 1.3 to 3.  The WARN level stays at four waves
#  but claims LESS than it looks:  above one wave the error is an erratic 2.6-19.3 % at fixed eta and does NOT
#  settle by 5.3 waves, and the 240 a.u. reference itself carries a fit residual of 0.337 -- so four waves is the
#  level above which this experiment cannot certify anything, NOT a level above which the answer is known good.

#  AND THE FIT RESIDUAL CANNOT DO THIS JOB, which is the finding that decided there should be a box tier at all.
#  It is ANTI-correlated at the dangerous end:  74.6 % error at residual 0.038 and 61.2 % at 0.140, both BELOW its
#  own 0.20 warning, while 2.6 % error reads 0.387 and warns.  The mechanism is plain once seen -- with less than
#  one wave in the box the fit has almost nothing to satisfy, so it matches the stub nearly perfectly.  A small
#  residual there means the fit is unconstrained, not that it is right.
const  CONT_BOX_WAVES_FATAL            = 1.0
const  CONT_BOX_WAVES_WARN             = 4.0


"""
`Continuum.occupiedCoreOf(level::Level, includeExchange::Bool)`  
    ... selects the orbitals that the frozen nonlocal exchange operator is taken with, together with their mean occupations.
        THE CORE IS THE OCCUPIED PART OF THE BASIS AND NOT ALL OF ITS ORBITALS. `level.basis.orbitals` carries every orbital
        the basis was built with, which for a configuration interaction on a shared or extended orbital set includes VIRTUALS
        that no CSF of this level occupies; `Basics.extractMeanOccupation` runs over `basis.subshells` and so returns no entry
        for them, and `Bsplines.setupExchangeMatrix` -- whose contract is one occupation per core orbital -- then raised "no
        occupation given for core subshell ...". That raise fired on a perfectly legitimate state: a shared orbital set is what
        an explicit second-order amplitude REQUIRES, and such a set carries subshells the final ion does not occupy. A virtual
        contributes nothing to the exchange sum, so the core is restricted here to the subshells that HAVE an occupation. This
        removes the raise and changes no number that was computable before it, since every orbital dropped carried no occupation
        and hence no weight. A tuple (coreOrbitals::Dict{Subshell,Orbital}, coreOccupations::Dict{Subshell,Float64}) is
        returned, or (nothing, nothing) where exchange is not wanted.
"""
function occupiedCoreOf(level::Level, includeExchange::Bool)
    if  !includeExchange    return( (nothing, nothing) )    end
    occupations  = Basics.extractMeanOccupation(level.basis)
    coreOrbitals = Dict{Subshell,Orbital}( sh => orb  for (sh, orb) in level.basis.orbitals  if haskey(occupations, sh) )

    return( (coreOrbitals, occupations) )
end


"""
`Continuum.generateOrbitalForLevel(energy::Float64, sh::Subshell, level::Level, nm::Nuclear.Model, grid::Radial.Grid,
                                    settings::Continuum.Settings; nuclearPot::Union{Nothing,Radial.Potential}=nothing, 
                                    primitives::Union{Nothing,Bsplines.Primitives}=nothing, 
                                    storage::Union{Nothing,Dict{String,Array{Float64,2}}}=nothing)`
    ... to generate a continuum orbital for the (continuum) subshell sh, the energy and the effective charge within the given
        potential. The continuum orbital is generated orthogonal with regard to all subshells of the same symmetry in the basis.
        All further specifications about this generations are made by proper settings. All three keywords carry quantities that are
        CONSTANT for a whole computation -- the nuclear potential depends only on the nuclear model and the grid, the B-spline
        primitives only on the grid, and the storage dictionary holds the TTp matrices, which depend on the grid and on kappa but
        never on the potential or the energy -- so that a caller which generates many continuum orbitals need not rebuild them for
        every line and every partial wave; omitting them reproduces the previous behaviour exactly. A tupel of a (continuum)
        (orbital::Orbital, phase::Float64) is returned.

        NOTE ON THE POTENTIAL: every continuum orbital is generated in a LOCAL (DFS) potential, whatever field the bound
        orbitals of `level` actually came from -- that remains true of every run and is stated here rather than emitted as
        a per-orbital warning as it was until 12-Aug-2026 (which produced one identical entry in jac-warn.report for every
        partial wave, and said nothing about the run). **AS OF 29-Oct-2026 (item 99) exchange with the bound electrons CAN
        be added on top of that local potential**, under `settings.includeExchange = true` and the Galerkin solution method
        only -- see `Continuum.Settings`'s own docstring for what that does and why every other combination raises.
"""
function generateOrbitalForLevel(energy::Float64, sh::Subshell, level::Level, nm::Nuclear.Model, grid::Radial.Grid,
                                 settings::Continuum.Settings; nuclearPot::Union{Nothing,Radial.Potential}=nothing, 
                                 primitives::Union{Nothing,Bsplines.Primitives}=nothing, 
                                 storage::Union{Nothing,Dict{String,Array{Float64,2}}}=nothing)

    # Generate a (local) potential for the given level.  The NUCLEAR part depends only on the nuclear model
    # and the grid, so it is constant for a whole computation -- but this function is called once per line
    # AND per partial wave, so rebuilding it here is pure repeated work.  Profiling an AutoIonization run
    # (10-Aug-2026) put ~19% of the wall clock in the nuclear charge density and a further ~16% in the
    # adaptive quadrature integrating it, of which this call site was responsible for 57.6%.
    # Callers may therefore pass a potential built once; nothing is cached here, and omitting the keyword
    # reproduces the previous behaviour exactly.
    nuclearPotential  = isnothing(nuclearPot) ? Nuclear.nuclearPotential(nm, grid) : nuclearPot
    # wp1 = compute("radial potential: core-Hartree", grid, wLevel)
    # wp2 = compute("radial potential: Hartree-Slater", grid, wLevel)
    # wp3 = compute("radial potential: Kohn-Sham", grid, wLevel)
    # wp4 = compute("radial potential: Dirac-Fock-Slater", grid, wLevel)
    # wp  = Basics.computePotential(Basics.DFSField(0.42), grid, level)   
    wp  = Basics.computePotential(Basics.DFSField(1.0), grid, level);
    pot = Basics.add(nuclearPotential, wp)

    # `includeExchange` is honoured in exactly one place -- the Galerkin solution method, where it costs no new
    # iteration (see Continuum.Settings' docstring, item 99). Elsewhere it has no mechanism to act through, and
    # raising here is cheaper to debug than a flag that was silently ignored.
    settings.includeExchange  &&  Defaults.GBL_CONT_SOLUTION != BsplineGalerkin()  &&
        error("Continuum.generateOrbitalForLevel(): settings.includeExchange = true needs the Galerkin solution " *
              "method (setDefaults(\"method: continuum, Galerkin\")); the current method has no nonlocal term to add it to.")
    coreOrbitals, coreOccupations = Continuum.occupiedCoreOf(level, settings.includeExchange)

    # Generate a continuum orbital due to the given solution method
    if      Defaults.GBL_CONT_SOLUTION  ==  ContBessel()
        cOrbital = Continuum.generateOrbitalBessel(energy, sh, grid::Radial.Grid, settings)
    elseif  Defaults.GBL_CONT_SOLUTION  ==  AsymptoticCoulomb()
        cOrbital = Continuum.generateOrbitalAsymptoticCoulomb(energy, sh, pot, settings)
    elseif  Defaults.GBL_CONT_SOLUTION  ==  NonrelativisticCoulomb()
        cOrbital = Continuum.generateOrbitalNonrelativisticCoulomb(energy, sh, pot, settings)
    elseif  Defaults.GBL_CONT_SOLUTION  ==  BsplineGalerkin()
        cOrbital = Continuum.generateOrbitalGalerkin(energy, sh, pot, settings; primitives=primitives, storage=storage,
                                                       coreOrbitals=coreOrbitals, occupations=coreOccupations)
    else    Continuum.unsupportedScheme("generateOrbitalForLevel", "Defaults.GBL_CONT_SOLUTION", Defaults.GBL_CONT_SOLUTION,
                                        CONT_IMPLEMENTED_SOLUTION_LEVEL, CONT_HOW_TO_SET_SOLUTION)
    end


    # Normalize the continuum orbital and determine its phase
    if      Defaults.GBL_CONT_NORMALIZATION  ==  PureSineNorm()
        cOrbital, phase = Continuum.normalizeOrbitalPureSine(cOrbital, pot.grid, settings; Zbar = Radial.determineZbar(pot))
    elseif  Defaults.GBL_CONT_NORMALIZATION  ==  CoulombSineNorm()
        cOrbital, phase = Continuum.normalizeOrbitalCoulombSine(cOrbital, pot, settings)
    elseif  Defaults.GBL_CONT_NORMALIZATION  ==  OngRussekNorm()
        cOrbital, phase = Continuum.normalizeOrbitalOngRussek(cOrbital, pot, settings)
    elseif  Defaults.GBL_CONT_NORMALIZATION  ==  CoulombFunctionNorm()
        cOrbital, phase = Continuum.normalizeOrbitalCoulombFunctions(cOrbital, pot, settings)
    else    Continuum.unsupportedScheme("generateOrbitalForLevel", "Defaults.GBL_CONT_NORMALIZATION", Defaults.GBL_CONT_NORMALIZATION,
                                        CONT_IMPLEMENTED_NORMALIZATION, CONT_HOW_TO_SET_NORMALIZATION)
    end

    #=open("continuum.txt","w") do file
        write(file,"energy $energy , n $(sh.n), kappa $(sh.kappa), phase $phase ")
        write(file,"potential")
        writedlm(file, pot.Zr)
        write(file,"radial grid")
        writedlm(file, pot.grid.r)
    end=#
    
    return( cOrbital, phase )
end


"""
`Continuum.generateOrbitalForLevel(energy::Float64, sh::Subshell, level::Level, nm::Nuclear.Model, grid::Radial.Grid,
                                    settings::Continuum.Settings, plasmaModel::Basics.AbstractPlasmaModel;
                                    nuclearPot::Union{Nothing,Radial.Potential}=nothing, 
                                    primitives::Union{Nothing,Bsplines.Primitives}=nothing, 
                                    storage::Union{Nothing,Dict{String,Array{Float64,2}}}=nothing)`
    ... to generate a continuum orbital as Continuum.generateOrbitalForLevel(...) above, but with the combined local
        potential (nuclear + bound-electron mean field) Debye-Hueckel-screened according to the given plasma model.
        The WHOLE combined potential is screened, not just its nuclear part: screening only the (attractive) nuclear
        term while leaving the (repulsive) mean-field term unscreened leaves an uncancelled repulsive tail at large r
        once the Debye length becomes comparable to the box radius, giving an unphysical (negative) asymptotic
        effective charge for the continuum electron. For plasmaModel = Basics.NoPlasmaModel(), this delegates
        directly to the field-free method above, so that the computation of the standard (field-free) processes is
        never touched by this plasma-specific code path. A tupel of a (continuum) (orbital::Orbital, phase::Float64)
        is returned.
"""
function generateOrbitalForLevel(energy::Float64, sh::Subshell, level::Level, nm::Nuclear.Model, grid::Radial.Grid,
                                  settings::Continuum.Settings, plasmaModel::Basics.AbstractPlasmaModel;
                                  nuclearPot::Union{Nothing,Radial.Potential}=nothing, 
                                  primitives::Union{Nothing,Bsplines.Primitives}=nothing, 
                                  storage::Union{Nothing,Dict{String,Array{Float64,2}}}=nothing)
    if  typeof(plasmaModel) == Basics.NoPlasmaModel
        return( Continuum.generateOrbitalForLevel(energy, sh, level, nm, grid, settings; nuclearPot=nuclearPot,
                                                  primitives=primitives, storage=storage) )
    end
    typeof(plasmaModel) == Basics.DebyeHueckelModel  ||
        error("Unsupported plasma model = $(plasmaModel)  (only Basics.DebyeHueckelModel is currently supported " *
              "for the plasma-screened continuum orbital).")

    # Generate the (local) potential for the given level exactly as in the field-free case, then Debye-Hueckel
    # screen the combined (nuclear + mean-field) potential as a whole -- see the docstring note above.
    nuclearPotential = isnothing(nuclearPot) ? Nuclear.nuclearPotential(nm, grid) : nuclearPot
    wp  = Basics.computePotential(Basics.DFSField(1.0), grid, level)
    pot = Basics.add(nuclearPotential, wp)
    lambda = 1/plasmaModel.debyeLength
    screenedZr = deepcopy(pot.Zr)
    for  i = 1:length(screenedZr)   screenedZr[i] = screenedZr[i] * exp(-lambda * pot.grid.r[i])    end
    pot = Radial.Potential(pot.name * "+ Debye-Hueckel screening", screenedZr, pot.grid)

    # `includeExchange` is honoured in exactly one place -- see the field-free method above. The exchange operator
    # is built from the ACTUAL (unscreened) bound orbitals and is deliberately NOT itself Debye-screened: it is a
    # short-range, finite-overlap operator (Bsplines.setupExchangeMatrix), not a long-range Coulomb tail, so it is
    # not the term the screening in this function exists to fix.
    settings.includeExchange  &&  Defaults.GBL_CONT_SOLUTION != BsplineGalerkin()  &&
        error("Continuum.generateOrbitalForLevel(): settings.includeExchange = true needs the Galerkin solution " *
              "method (setDefaults(\"method: continuum, Galerkin\")); the current method has no nonlocal term to add it to.")
    coreOrbitals, coreOccupations = Continuum.occupiedCoreOf(level, settings.includeExchange)

    # Generate a continuum orbital due to the given solution method
    if      Defaults.GBL_CONT_SOLUTION  ==  ContBessel()
        cOrbital = Continuum.generateOrbitalBessel(energy, sh, grid::Radial.Grid, settings)
    elseif  Defaults.GBL_CONT_SOLUTION  ==  AsymptoticCoulomb()
        cOrbital = Continuum.generateOrbitalAsymptoticCoulomb(energy, sh, pot, settings)
    elseif  Defaults.GBL_CONT_SOLUTION  ==  NonrelativisticCoulomb()
        cOrbital = Continuum.generateOrbitalNonrelativisticCoulomb(energy, sh, pot, settings)
    elseif  Defaults.GBL_CONT_SOLUTION  ==  BsplineGalerkin()
        cOrbital = Continuum.generateOrbitalGalerkin(energy, sh, pot, settings; primitives=primitives, storage=storage,
                                                       coreOrbitals=coreOrbitals, occupations=coreOccupations)
    else    Continuum.unsupportedScheme("generateOrbitalForLevel", "Defaults.GBL_CONT_SOLUTION", Defaults.GBL_CONT_SOLUTION,
                                        CONT_IMPLEMENTED_SOLUTION_LEVEL, CONT_HOW_TO_SET_SOLUTION)
    end

    # Normalize the continuum orbital and determine its phase
    if      Defaults.GBL_CONT_NORMALIZATION  ==  PureSineNorm()
        cOrbital, phase = Continuum.normalizeOrbitalPureSine(cOrbital, pot.grid, settings; Zbar = Radial.determineZbar(pot))
    elseif  Defaults.GBL_CONT_NORMALIZATION  ==  CoulombSineNorm()
        cOrbital, phase = Continuum.normalizeOrbitalCoulombSine(cOrbital, pot, settings)
    elseif  Defaults.GBL_CONT_NORMALIZATION  ==  OngRussekNorm()
        cOrbital, phase = Continuum.normalizeOrbitalOngRussek(cOrbital, pot, settings)
    elseif  Defaults.GBL_CONT_NORMALIZATION  ==  CoulombFunctionNorm()
        cOrbital, phase = Continuum.normalizeOrbitalCoulombFunctions(cOrbital, pot, settings)
    else    Continuum.unsupportedScheme("generateOrbitalForLevel", "Defaults.GBL_CONT_NORMALIZATION", Defaults.GBL_CONT_NORMALIZATION,
                                        CONT_IMPLEMENTED_NORMALIZATION, CONT_HOW_TO_SET_NORMALIZATION)
    end

    return( cOrbital, phase )
end


"""
`Continuum.generateOrbitalLocalPotential(energy::Float64, sh::Subshell, pot::Radial.Potential, settings::Continuum.Settings)`
    ... to generate a continuum orbital for the (continuum) subshell sh, the energy and the given potential. The effective
        charge for the normalization of the continuum orbital is derived from the potential. All further specifications about
        this generations are made by proper settings; however, the function terminates if settings.includeExchange = true,
        UNCONDITIONALLY and for a different reason than in `Continuum.generateOrbitalForLevel` (item 99): this function is
        handed a bare `pot::Radial.Potential` with no `level` at all, so it has no bound orbitals to build an exchange
        operator from, whatever the solution method.
        A tupel of a (continuum) (orbital::Orbital, phase::Float64, normFactor::Float64) is returned.
"""
function generateOrbitalLocalPotential(energy::Float64, sh::Subshell, pot::Radial.Potential, settings::Continuum.Settings)
    Continuum.checkNormalizationRadius(pot)
    settings.includeExchange   &&   error("Continuum orbital for local potential does not allow 'exchange'.")
    
    # Generate a continuum orbital due to the given solution method
    if      Defaults.GBL_CONT_SOLUTION  ==  ContBessel() 
        cOrbital = Continuum.generateOrbitalBessel(energy, sh, pot.grid, settings)
    elseif  Defaults.GBL_CONT_SOLUTION  ==  ContSine()
        cOrbital = Continuum.generateOrbitalPureSine(energy, sh, pot.grid, settings)
    elseif  Defaults.GBL_CONT_SOLUTION  ==  AsymptoticCoulomb()
        cOrbital = Continuum.generateOrbitalAsymptoticCoulomb(energy, sh, pot, settings)
    elseif  Defaults.GBL_CONT_SOLUTION  ==  NonrelativisticCoulomb()
        cOrbital = Continuum.generateOrbitalNonrelativisticCoulomb(energy, sh, pot, settings)
    elseif  Defaults.GBL_CONT_SOLUTION  ==  BsplineGalerkin()
        cOrbital = Continuum.generateOrbitalGalerkin(energy, sh, pot, settings)
    else    Continuum.unsupportedScheme("generateOrbitalLocalPotential", "Defaults.GBL_CONT_SOLUTION", Defaults.GBL_CONT_SOLUTION,
                                        CONT_IMPLEMENTED_SOLUTION_POT, CONT_HOW_TO_SET_SOLUTION)
    end

    # Normalize the continuum orbital and determine its phase
    if      Defaults.GBL_CONT_NORMALIZATION  ==  PureSineNorm()
        cOrbital, phase, normFactor = Continuum.normalizeOrbitalPureSine(cOrbital, pot.grid, settings; Zbar = Radial.determineZbar(pot))
    elseif  Defaults.GBL_CONT_NORMALIZATION  ==  CoulombSineNorm()
        cOrbital, phase, normFactor = Continuum.normalizeOrbitalCoulombSine(cOrbital, pot, settings)
    elseif  Defaults.GBL_CONT_NORMALIZATION  ==  OngRussekNorm()
        cOrbital, phase, normFactor = Continuum.normalizeOrbitalOngRussek(cOrbital, pot, settings)
    elseif  Defaults.GBL_CONT_NORMALIZATION  ==  CoulombFunctionNorm()
        cOrbital, phase             = Continuum.normalizeOrbitalCoulombFunctions(cOrbital, pot, settings)
        normFactor                  = 1.0
    else    Continuum.unsupportedScheme("generateOrbitalLocalPotential", "Defaults.GBL_CONT_NORMALIZATION", Defaults.GBL_CONT_NORMALIZATION,
                                        CONT_IMPLEMENTED_NORMALIZATION, CONT_HOW_TO_SET_NORMALIZATION)
    end

    return( cOrbital, phase, normFactor )
end


"""
`Continuum.generateOrbitalAsymptoticCoulomb(energy::Float64, sh::Subshell, pot::Radial.Potential, settings::Continuum.Settings)`   
    ... to generate a simple relativistic and non-normalized continuum orbital but by assuming an asymptotic sin/cos behaviour 
        for both components at all mesh points. For this behaviour, both components and their derivatives can be calculated 
        analytically. The effective charge Zbar is determined for the given potential and the non-Coulombic phase shifts is set 
        to zero. An non-normalized cOrbital::Orbital is returned.
        
        Warning: At present, only the real part eta.re is taken into account; check for consistency !!
"""
function generateOrbitalAsymptoticCoulomb(energy::Float64, sh::Subshell, pot::Radial.Potential, settings::Continuum.Settings)  
    Continuum.checkNormalizationRadius(pot)
    P = zeros(settings.mtp);    Q = zeros(settings.mtp);   Pprime = zeros(settings.mtp);    Qprime = zeros(settings.mtp)
    
    phi0 = 0.3
    wc   = Defaults.getDefaults("speed of light: c");    q = sqrt( 2*energy + energy/(wc*wc));    kappa = sh.kappa
    Zbar = Radial.determineZbar(pot)
    y    = Zbar * (energy + wc^2) / (wc^2 * q);   gammaBar = sqrt(kappa^2 - Zbar^2 / wc^2)
    eta  = - (kappa - im*y) / (energy + wc^2) / (gammaBar + im*y) / (2im)
    NP   = sqrt( (energy + 2*wc^2)/(pi*wc^2*q) );   NQ = -sqrt( energy / (pi*wc^2*q) )
    
    # println("Continuum-asymptotic Coulomb orbital for q = $q,  kappa = $kappa,  Zbar=$Zbar  y = $y,   gammaBar = $gammaBar,  eta = $eta")
    if  Defaults.GBL_PRINT_DEBUG  println("Continuum-asymptotic Coulomb orbital for q = $q,  kappa = $kappa,  Zbar=$Zbar,  eta = $eta")    end

    
    for  i = 2:settings.mtp  
        theta = q * pot.grid.r[i]  +  y * log(2q * pot.grid.r[i])  - angle( SpecialFunctions.gamma(gammaBar + im*y) )  -  pi*gammaBar/2  +  eta.re
        P[i] = NP * cos(theta + phi0);      Pprime[i] = - NP * (q + y/pot.grid.r[i]) * sin(theta + phi0)       
        Q[i] = - NQ * sin(theta + phi0);    Qprime[i] = - NQ * (q + y/pot.grid.r[i]) * cos(theta + phi0)       
    end

    cOrbital = Orbital( sh, false, true, energy, P, Q, Pprime, Qprime, Radial.Grid())
    
    return( cOrbital )
end


"""
`Continuum.generateOrbitalBessel(energy::Float64, sh::Subshell, grid::Radial.Grid, settings::Continuum.Settings)`   
    ... to generate a simple, non-relativistic and non-normalized Bessel wave for the large component of the continuum
        orbital and a small component due to the kinetic-balance condition. A (non-normalized) orbital::Orbital is returned.
        While Pprime is obtained from the analytical expression of the large component, Qprime is set to zero here.
"""
function generateOrbitalBessel(energy::Float64, sh::Subshell, grid::Radial.Grid, settings::Continuum.Settings)  
    P = zeros(settings.mtp);    Q = zeros(settings.mtp);   Pprime = zeros(settings.mtp);    Qprime = zeros(settings.mtp)
    
    q = sqrt( 2*energy );    l = Basics.subshell_l(sh);   wc = Defaults.getDefaults("speed of light: c")
    for  i = 2:settings.mtp   P[i] = GSL.sf_bessel_jl(l, q * grid.r[i]) * grid.r[i]    end
    for  i = 2:settings.mtp   qz = q * grid.r[i];   Pprime[i] = qz * (l/qz * GSL.sf_bessel_jl(l, qz) - GSL.sf_bessel_jl(l+1, qz)) + P[i]/grid.r[i]   end
    for  i = 2:settings.mtp   Q[i] = -1/(2 * wc) * Pprime[i]  +  sh.kappa/grid.r[i] * P[i]         end
    
    println("Continuum-Bessel orbital for q = $q,  l = $l")

    cOrbital = Orbital( sh, false, false, energy, P, Q, Pprime, Qprime, grid)
    
    return( cOrbital )
end


"""
`Continuum.generateOrbitalGalerkin(energy::Float64, sh::Subshell, pot::Radial.Potential, settings::Continuum.Settings;
                                    primitives::Union{Nothing,Bsplines.Primitives}=nothing,
                                    storage::Union{Nothing,Dict{String,Array{Float64,2}}}=nothing,
                                    coreOrbitals::Union{Nothing,Dict{Subshell,Orbital}}=nothing,
                                    occupations::Union{Nothing,Dict{Subshell,Float64}}=nothing)`
    ... to generate a non-normalized continuum orbital within the given local potential by using the Galerkin method and a
        given B-spline basis. Both keywords carry quantities that are constant for a whole computation on one grid: the
        B-spline primitives depend on the radial grid alone, and the storage dictionary holds the TTp matrices, which depend
        on the grid and on kappa but never on the potential or the energy. A caller that generates many continuum orbitals
        should build both once and pass them in. Omitting either reproduces the previous behaviour exactly. A
        (non-normalized) orbital::Orbital is returned.

        **`coreOrbitals`/`occupations`, ADDED 29-Oct-2026 (item 99), ADD EXCHANGE WITH A FROZEN CORE, AND NEED NO ITERATION
        HERE.** When both are given, `Bsplines.setupExchangeMatrix` -- the same frozen nonlocal Dirac-Fock exchange operator
        validated for item 75/S30 -- is added to the Galerkin matrix before it is solved. The core is frozen and the operator
        is LINEAR in the unknown continuum orbital, so this is exactly the one-shot assembly `A = setupLocalMatrix +
        setupExchangeMatrix` already uses for a frozen-core bound Fock matrix, not a new convergence loop: the matrix is built
        once, and the same single generalized-eigenvalue solve below finds the continuum orbital, now of the Dirac-FOCK
        problem instead of the Dirac-Hartree one. Omitting either keyword (the default) reproduces the previous behaviour
        exactly -- this is the switch `Continuum.Settings.includeExchange` turns on, and its own docstring says why every
        other solution method must raise instead of reaching this keyword.
"""
function generateOrbitalGalerkin(energy::Float64, sh::Subshell, pot::Radial.Potential, settings::Continuum.Settings;
                                  primitives::Union{Nothing,Bsplines.Primitives}=nothing,
                                  storage::Union{Nothing,Dict{String,Array{Float64,2}}}=nothing,
                                  coreOrbitals::Union{Nothing,Dict{Subshell,Orbital}}=nothing,
                                  occupations::Union{Nothing,Dict{Subshell,Float64}}=nothing)
    Continuum.checkNormalizationRadius(pot)
    P = zeros(settings.mtp);   Q = zeros(settings.mtp);   Pprime = zeros(settings.mtp);    Qprime = zeros(settings.mtp)
    nsL = pot.grid.nsL - 1;    nsS = pot.grid.nsS - 1
    # The B-spline basis depends only on the grid, but this function is called once per line AND per partial wave;
    # profiling the Auger rates alone (11-Aug-2026) put 34.3% of them into generatePrimitives.  Nothing is cached
    # here: a caller that knows its grid is fixed builds the basis once and hands it in.
    if  isnothing(primitives)   wa = Bsplines.generatePrimitives(pot.grid)
    else                        wa = primitives
        wa.grid.nsL == pot.grid.nsL  &&  wa.grid.nsS == pot.grid.nsS  &&  wa.grid.NoPoints == pot.grid.NoPoints  ||
            error("The given primitives belong to a different grid than the potential.")
    end
    wb = Bsplines.generateGalerkinMatrix(sh, energy, pot, wa; storage=storage)
    if  !isnothing(coreOrbitals)
        wb = wb + Bsplines.setupExchangeMatrix(sh.kappa, wa, coreOrbitals, occupations)
    end
    wc = adjoint(wb) * wb
    
    # Test for 'real-symmetric matrix' ... and symmetrize otherwise
    nx = 0
    for  i = 1:nsL+nsS    
        for  j = i+1:nsL+nsS    
            if  abs(  (wc[i,j] - wc[j,i])/(wc[i,j] + wc[j,i]) ) > 1.0e-7   nx = nx + 1    
            end
        end
    end
    # Only the LOWEST eigenpair is used below -- the null-space vector of wc = wb'*wb, which is the Galerkin
    # solution.  Asking LAPACK for that one pair (syevr) instead of decomposing the whole matrix is ~2.7x
    # faster at this size (671x671, measured 12-Aug-2026); the eigenvector agrees with the full solve to
    # ~1e-12, and the eigenvalue -- which is ~0 by construction and enters only the printout below -- to the
    # same absolute accuracy.  Profiling the Auger rates put 18.7% of them in this call.
    wd = Basics.diagonalize(MatrixWithLinearAlgebra(), wc; range=1:1)
    # println(">>> Galerkin-eigenvalues = $(wd.values[1]), $(wd.values[2]) for  $sh  with  energy = $energy")
    
    cOrbital = Bsplines.generateOrbitalFromPrimitives(sh, energy, settings.mtp, wd.vectors[1], wa)  
    mtp      = size(cOrbital.P,1)
    println(">> Continuum B-spline-Galerkin orbital for energy=" * @sprintf("%.4e",energy) * ",  kappa=$(sh.kappa) " *
            "[mpt=$mtp, r[mtp]=" * @sprintf("%.4e",pot.grid.r[mtp]) * ", smallest eigenvalue=" * @sprintf("%.4e",wd.values[1]) * "].")
    
    return( cOrbital )
end


"""
`Continuum.generateOrbitalNonrelativisticCoulomb(energy::Float64, sh::Subshell, Zeff::Float64, grid::Radial.Grid, 
                                                    settings::Continuum.Settings)`   
    ... RAISES.  It was intended to generate a simple, non-relativistic and non-normalized Coulomb wave for the (continuum)
        subshell sh, with the small component from the kinetic-balance condition, but it has never produced an orbital and no
        `orbital::Orbital` is returned:  the Coulomb wave needs a confluent hypergeometric function 1F1 with COMPLEX arguments,
        `1F1(l+1+i*eta, 2l+2, 2iqr)`, and `GSL.hypergeom` accepts real arguments only.

        WHY IT RAISES RATHER THAN WARNS, since 03-Nov-2026.  That limitation was already recorded here and at the head of the
        function, but as a `println`, which goes to stdout and is suppressed by every `redirect_stdout(devnull)` path in the
        package;  execution then continued into the loops, which assign the complex result into a `Float64` array and die on
        an `InexactError` naming neither this function nor the setting that selected it.  Until the same date all three
        dispatchers compounded it by passing a `Radial.Grid` into the `Radial.Potential` slot of the four-argument method, so
        the branch raised a `MethodError` on a signature before ever reaching the physics.  Both are now a single explained
        refusal that names the alternatives:  `AsymptoticCoulomb` for an analytic Coulomb wave, `BsplineGalerkin` to solve in
        the actual potential.  The formulas are left in place for whoever supplies a complex 1F1.
"""
function generateOrbitalNonrelativisticCoulomb(energy::Float64, sh::Subshell, Zeff::Float64, grid::Radial.Grid, settings::Continuum.Settings) 
    #  A LOUD STOP, WHERE A SWALLOWED println USED TO BE.  This was already the author's intent -- the line that
    #  stood here said "presently, NO Julia implementation of the hypergeometric function with complex arguments
    #  available" -- but as a `println` it goes to stdout, is suppressed by every `redirect_stdout(devnull)` path in
    #  the package, and execution then continued into the loops below, which pass COMPLEX arguments to GSL.hypergeom
    #  (a real-only 1F1) and assign a complex result into a Float64 array.  Measured 03-Nov-2026, that surfaces as
    #  `InexactError: Float64(1.0 + 0.9027im)` from inside this function, with nothing to connect it to the setting
    #  the user chose.  The formulas below are left in place as the intended implementation, for whoever supplies a
    #  complex 1F1;  they have never produced an orbital.
    error("\n\nContinuum.generateOrbitalNonrelativisticCoulomb():  STOP -- this solution method is NOT implemented.\n" *
          "    requested by           : setDefaults(\"method: continuum, nonrelativistic Coulomb\")\n" *
          "    what is missing        : a confluent hypergeometric function 1F1 with COMPLEX arguments;  GSL.hypergeom\n" *
          "                             accepts real arguments only, and the Coulomb wave needs 1F1(l+1+i*eta, 2l+2, 2iqr)\n" *
          "    use instead            : AsymptoticCoulomb -- an analytic Coulomb wave, fitted to its own asymptotic form\n" *
          "                             to 7e-16 -- or BsplineGalerkin, which solves in the actual potential\n" *
          "This method has produced no orbital since the first commit of the public history:  until 03-Nov-2026 all three\n" *
          "dispatchers also called it with a Radial.Grid in the Radial.Potential slot, so it raised a MethodError on a\n" *
          "signature before reaching this point.  That call is now correct and the refusal is deliberate and explained.\n")
    
    # Define the derivative of the confluent hypergeometric function
    function hgPrime(a,b,z)    wa = a/z * (GSL.hypergeom(a+1,b,z) - GSL.hypergeom(a,b,z));   return(wa)   end
    
    P = zeros(settings.mtp);    Q = zeros(settings.mtp);   Pprime = zeros(settings.mtp);    Qprime = zeros(settings.mtp)
    
    q  = sqrt( 2*energy );    l = Basics.subshell_l(sh);   wc = Defaults.getDefaults("speed of light: c");   np = Zeff/q;   ws = 1
    wa = 2 * sqrt(Zeff) / sqrt(1 - Base.MathConstants.e^(-2pi*np)) / factorial(2l + 1)
    for  s = 1:l    ws = ws * sqrt(s^2 + np^2)    end
    
    for  i = 2:settings.mtp   P[i] = wa * ws * (2q*grid.r[i])^l * Base.MathConstants.e^(-im*q*grid.r[i]) * 
                                        GSL.hypergeom(l+1+im*np, 2l+2., 2im*q*grid.r[i])                     end
    for  i = 2:settings.mtp   
        Pprime[i] = 2q*l * (2q*grid.r[i])^(l-1) * Base.MathConstants.e^(-im*q*grid.r[i]) * GSL.hypergeom(l+1+im*np, 2l+2., 2im*q*grid.r[i]) -
                    im*q * (2q*grid.r[i])^l * Base.MathConstants.e^(-im*q*grid.r[i]) * GSL.hypergeom(l+1+im*np, 2l+2., 2im*q*grid.r[i])     +
                    (2q*grid.r[i])^l * Base.MathConstants.e^(-im*q*grid.r[i]) * hgPrime(l+1+im*np, 2l+2., 2im*q*grid.r[i])
        Pprime[i] = wa * ws * Pprime[i]          
    end
    for  i = 2:settings.mtp   Q[i] = -1/(2 * wc) * Pprime[i]  +  sh.kappa/grid.r[i] * P[i]                end

    cOrbital = Orbital( sh, false, true, energy, P, Q, Pprime, Qprime, Radial.Grid())
    
    return( cOrbital )
end


"""
`Continuum.generateOrbitalNonrelativisticCoulomb(energy::Float64, sh::Subshell, pot::Radial.Potential, settings::Continuum.Settings)`   
    ... to generate a simple, non-relativistic and non-normalized Coulomb wave for the (continuum) subshell sh, 
        the small component is simply obtained by the kinetic-balance condition. A orbital::Orbital is returned.
"""
function generateOrbitalNonrelativisticCoulomb(energy::Float64, sh::Subshell, pot::Radial.Potential, settings::Continuum.Settings) 
    Zbar     = Radial.determineZbar(pot)
    cOrbital = Continuum.generateOrbitalNonrelativisticCoulomb(energy, sh, Zbar, pot.grid, settings)
    
    return( cOrbital )
end


"""
`Continuum.generateOrbitalPureSine(energy::Float64, sh::Subshell, grid::Radial.Grid, settings::Continuum.Settings)`   
    ... to generate a simple, non-relativistic and non-normalized Bessel wave for the large component of the continuum
        orbital and a small component due to the kinetic-balance condition. A (non-normalized) orbital::Orbital is returned.
        While Pprime is obtained from the analytical expression of the large component, Qprime is set to zero here.
"""
function generateOrbitalPureSine(energy::Float64, sh::Subshell, grid::Radial.Grid, settings::Continuum.Settings)  
    P = zeros(settings.mtp);    Q = zeros(settings.mtp);   Pprime = zeros(settings.mtp);    Qprime = zeros(settings.mtp)
    
    phi0 = 0.0
    q = sqrt( 2*energy );    l = Basics.subshell_l(sh);   wc = Defaults.getDefaults("speed of light: c");    N = sqrt(2/(pi*q))
    for  i = 1:settings.mtp   P[i]      = N * sin(q * grid.r[i] - l*pi/2 + phi0)                  end
    for  i = 1:settings.mtp   Pprime[i] = N * q * cos(q * grid.r[i] - l*pi/2 + phi0)              end
    for  i = 2:settings.mtp   Q[i] = -1/(2 * wc) * Pprime[i]  +  sh.kappa/grid.r[i] * P[i]        end
    
    println(">> Continuum-pure sine orbital for q = $q,  l = $l")

    cOrbital = Orbital( sh, false, false, energy, P, Q, Pprime, Qprime, grid)
    
    return( cOrbital )
end


"""
`Continuum.checkContinuumAmplitude(P::Array{Float64,1}, energy::Float64, grid::Radial.Grid)`
    ... checks that a NORMALIZED continuum orbital's large component does not blow up in the interior relative to its own
        asymptotic amplitude, which is how a radial mesh too coarse for the requested energy makes itself visible in the orbital
        it produced rather than in the grid it was produced on.  A NamedTuple
        `(ok, checked, interior, asymptotic, ratio, required, rIndex, note)` is returned; nothing is raised here, so that the
        CALLER judges -- the same contract as `checkNormalizationRadius` below.

        WHY A POSTERIOR CHECK IS NEEDED AT ALL, when three guards already look at the grid.  `gridConsistency` compares
        `grid.hp` with the shortest de Broglie wavelength and requires 15 points across it; `checkNormalizationRadius` asks
        whether the normalization point lies outside the atom; `Bsplines.checkGridRepresentation` and `checkOrbitalBox` test
        BOUND orbitals.  **ALL THREE ARE A PRIORI** -- they judge the mesh against the energy, before any orbital exists -- and
        a mesh can satisfy every one of them and still fail.  Measured 22-Oct-2026 on the 3.583e+04 Hartree electron of the
        207Pb K-shell internal-conversion case, Z = 82, kappa = -4, where BOTH meshes below pass `gridConsistency`:

            points per oscillation   15.65        46.94
            max|P| anywhere          2.03e+08     6.59e-02
            max|P| asymptotically    6.12e-02     5.77e-02
            RATIO                    3.32e+09     1.14

        The 15.65 mesh is the one `example-Dq.jl` used, and the orbital it returns is wrong by ~4e+09 in amplitude and ~1e+19 in
        amplitude SQUARED while every a-priori guard is satisfied.  The spike sits at r = 6e-05 a.u., i.e. in the innermost
        region, where a continuum orbital goes as r^gamma and must therefore be SMALL -- which is what makes the ratio a sharp
        instrument rather than a heuristic.

        THE THRESHOLD IS MEASURED, AND THE RATIO IS CONVENTION-INDEPENDENT.  Over 144 cases spanning Z = 2, 18, 54, 92,
        energies 0.01 to 50 Hartree and kappa from -1 to 8 -- including near-threshold and high centrifugal barriers, the two
        regimes where a legitimately large interior amplitude might be expected -- the WORST ratio was 1.269, while the failing
        cases give 3.3e+09 and 2.3e+10.  So `required = 1.0e3` sits three orders above anything legitimate yet measured and six
        orders below the defect.  And because the ratio divides one amplitude by another, the normalization CONSTANT cancels:
        the same good case gives 1.141 under all four of `PureSineNorm`, `CoulombSineNorm`, `OngRussekNorm` and `CoulombFunctionNorm`, so
        one threshold serves every convention and the check cannot be evaded by changing the global.

        NOT CHECKED IS NOT PASSED, and `checked` is the field that keeps the two apart.  A grid with fewer than 400 points has
        no asymptotic region to measure against, and an identically vanishing asymptotic amplitude leaves the ratio undefined;
        both return `ok = true` because nothing failed, and a caller reading only `ok` would take either for a clean bill of
        health.
"""
function checkContinuumAmplitude(P::Array{Float64,1}, energy::Float64, grid::Radial.Grid)
    required = 2.0;     fatal = 1.0e2
    if  length(P) < 400
        return( (ok = true, checked = false, interior = NaN, asymptotic = NaN, ratio = NaN, required = required, fatal = fatal, rIndex = 0,
                 nWindow = 0, note = "not checked: fewer than 400 grid points, so there is no asymptotic region to compare against") )
    end
    # THE WINDOW IS SIZED BY THE WAVELENGTH, NOT BY A FIXED POINT COUNT, and that is the whole correctness of the
    # measure rather than a refinement of it.  A window shorter than half an oscillation samples the outgoing wave at
    # one arbitrary phase, so landing near a node makes the "asymptotic amplitude" arbitrarily small and the ratio
    # arbitrarily large.  Measured 22-Oct-2026: a fixed 300-point window at hp = 5e-3 and E = 0.05 Hartree spans
    # 1.5 a.u. of a 19.87 a.u. wavelength -- 0.075 of one oscillation -- and inflated 8 of 120 SOUND screened orbitals
    # above 2, the worst to 3.96.  Half a wavelength always contains one antinode, so it cannot miss the amplitude.
    wavenb   = sqrt( 2energy + energy * Defaults.getDefaults("alpha")^2 )
    nHalf    = wavenb > 0.0  &&  grid.hp > 0.0  ?  Int64(ceil( (pi/wavenb) / grid.hp ))  :  300
    nWindow  = min( max(nHalf, 300), div(length(P), 3) )
    if  nWindow < 300
        return( (ok = true, checked = false, interior = maximum( abs.(P) ), asymptotic = NaN, ratio = NaN,
                 required = required, fatal = fatal, rIndex = argmax( abs.(P) ), nWindow = nWindow,
                 note = "not checked: the box holds less than three asymptotic windows, so interior and asymptotic " *
                        "regions cannot be separated") )
    end
    asymptotic = maximum( abs.(P[end-nWindow+1:end]) )
    interior   = maximum( abs.(P) );     rIndex = argmax( abs.(P) )
    if  asymptotic <= 0.0
        return( (ok = true, checked = false, interior = interior, asymptotic = asymptotic, ratio = NaN, required = required,
                 fatal = fatal, rIndex = rIndex, nWindow = nWindow,
                 note = "not checked: the asymptotic amplitude vanishes identically, so the ratio is undefined") )
    end
    ratio = interior / asymptotic

    return( (ok = ratio <= required, checked = true, interior = interior, asymptotic = asymptotic, ratio = ratio,
             required = required, fatal = fatal, rIndex = rIndex, nWindow = nWindow, note = "") )
end


"""
`Continuum.reportContinuumAmplitude(P::Array{Float64,1}, cOrbital::Orbital, grid::Radial.Grid, sa::String; residual::Float64=NaN)`
    ... validates the normalized large component `P` of `cOrbital` in the two stages of the module's guard order, naming the
        convention `sa` that produced it:  first its FORM, against the `residual` of the asymptotic fit below, then its
        MAGNITUDE, through `Continuum.checkContinuumAmplitude`'s interior/asymptotic ratio.  They fail independently, which is
        why both are asked and in that order.  Each stage warns where the result is suspect and raises where it is unusable;
        `nothing` is returned when neither fires.

        + residual ::Float64   ... the rms misfit of the orbital against its fitted asymptotic form, as a FRACTION of the
                                   amplitude, from `Continuum.asymptoticAmplitude`.  It is checked in two tiers before the
                                   amplitude ratio below, because it is the only one of the two that sees a non-asymptotic
                                   orbital at all.  `NaN`, the default, means the caller has no fit to offer -- which is the
                                   case for `normalizeOrbitalCoulombFunctions`, matching analytic functions instead -- and the tier is
                                   then skipped rather than guessed at.

        THE WARNING GOES TO `stderr` AND IS THEREFORE NOT SUPPRESSED by a `printout = false` path or by a
        `redirect_stdout(devnull)`, which is deliberate: the orbital it concerns is unusable for every observable, and the
        failure is otherwise entirely silent.  It WARNS rather than raises, following `checkNormalizationRadius` below, whose
        own note records why a normalization check is not promoted to an `error` -- there the penalty for violating the
        criterion is observable-dependent by four orders.  That reasoning is respected here and is also why the threshold is
        set six orders below the measured defect rather than at the edge of the legitimate range.

        THE GRID IS AN ARGUMENT AND MUST BE, although an `Orbital` carries a `grid` field of its own: on this path that field
        comes back EMPTY -- `length(cOrbital.grid.r) == 0` for an orbital produced by `generateOrbitalLocalPotential` -- so
        reading the spike radius from it yields `NaN` and the most useful line of the warning says nothing.  Measured
        22-Oct-2026 while calibrating this check.
"""
function reportContinuumAmplitude(P::Array{Float64,1}, cOrbital::Orbital, grid::Radial.Grid, sa::String;
                                  residual::Float64=NaN)
    #  THE FIT RESIDUAL IS THE ONLY ONE OF THESE MEASURES THAT SEES A NON-ASYMPTOTIC ORBITAL, measured 03-Nov-2026
    #  and it is why `abs(N) > 30` could be deleted rather than merely corrected.  Over sixteen SOUND He+ orbitals,
    #  0.5 to 200 eV in 24 and 40 a.u. boxes, the interior/asymptotic ratio read 1.000 or 1.001 -- and it read
    #  1.000 for all seven DELIBERATELY ruined ones too, slow electrons in boxes far too small for them, so the
    #  ratio cannot separate the two cases at all.  The residual separates them by orders of magnitude:

    #      sound     2.2e-04 (200 eV, 40 a.u.)  ...  1.8e-01 (0.5 eV, 24 a.u.)
    #      ruined    1.9e-01 (1.0 eV,  8 a.u.)  ...  1.8e+01 (0.05 eV, 12 a.u.)

    #  Hence `required` 0.2 and `fatal` 1.0 below:  at a residual of 1.0 the misfit is as large as the amplitude it
    #  is a fraction of, so the asymptotic form the normalization assumed does not describe this orbital.  A NaN
    #  means the caller has no fit to offer -- normalizeOrbitalCoulombFunctions, which matches analytic functions instead --
    #  and the tier is then skipped rather than guessed at.
    resRequired = 0.2;    resFatal = 1.0
    if       !isnan(residual)  &&  residual > resFatal
        error("\n\nContinuum.normalizeOrbital" * sa * "():  STOP -- the orbital does not have the asymptotic form " *
              "that this normalization assumes.\n" * Continuum.asymptoticFitReport(residual, resRequired, resFatal, cOrbital) )
    elseif   !isnan(residual)  &&  residual > resRequired
        @warn("\n\nContinuum.normalizeOrbital" * sa * "():  the orbital only approximately has the asymptotic form " *
              "that this normalization assumes.\n" * Continuum.asymptoticFitReport(residual, resRequired, resFatal, cOrbital) )
    end

    wa = Continuum.checkContinuumAmplitude(P, cOrbital.energy, grid)
    if  wa.checked  &&  wa.ratio > wa.fatal
        # THE TWO TIERS DO DIFFERENT JOBS AND THE SPLIT IS MEASURED, not a hedge.  The normalization is fixed by the
        # ASYMPTOTIC amplitude, so a ratio of a few means the observable is off by percent at worst -- refusing such a
        # run is the mistake `checkNormalizationRadius` below records as its own reason for warning rather than
        # raising.  A ratio of 1e2 or more cannot be read that way: the orbital carries no usable amplitude at all,
        # and every number built from it is wrong by that factor or its square.  This note USED to add that the
        # three sibling conventions refuse the same grid on their own `abs(N) > 30` test;  that test was deleted on
        # 03-Nov-2026 as uncorrelated with usability, and the residual tier above is what now carries that duty.
        error("\n\nContinuum.normalizeOrbital" * sa * "():  STOP -- the continuum orbital is not represented on this grid.\n" *
              Continuum.continuumAmplitudeReport(wa, cOrbital, grid) )
    elseif  wa.checked  &&  !wa.ok
        @warn("\n\nContinuum.normalizeOrbital" * sa * "():  the continuum orbital may not be represented on this grid.\n" *
              Continuum.continuumAmplitudeReport(wa, cOrbital, grid) )
    end

    return( nothing )
end


"""
`Continuum.continuumAmplitudeReport(wa::NamedTuple, cOrbital::Orbital, grid::Radial.Grid)`
    ... formats the measurement from `Continuum.checkContinuumAmplitude` into the explanatory block that both the warning and
        the refusal in `Continuum.reportContinuumAmplitude` print, so that the two cannot drift apart; a `report::String` is
        returned.
"""
function continuumAmplitudeReport(wa::NamedTuple, cOrbital::Orbital, grid::Radial.Grid)
    energy   = cOrbital.energy
    wavenb   = sqrt( 2energy + energy * Defaults.getDefaults("alpha")^2 );      wavelgth = 2pi / wavenb
    rSpike   = wa.rIndex <= length(grid.r) ? grid.r[wa.rIndex] : NaN

    return( @sprintf("    subshell, energy         = %s,  %.4e Hartree\n", string(cOrbital.subshell), energy) *
            @sprintf("    max|P| anywhere          = %.4e   at r = %.3e a.u.\n", wa.interior, rSpike) *
            @sprintf("    max|P| asymptotically    = %.4e   (over the outer %d points = %.2f oscillations)\n",
                     wa.asymptotic, wa.nWindow, wa.nWindow * grid.hp / wavelgth) *
            @sprintf("    interior / asymptotic    = %.4e   (above %.1f is suspect, above %.1e unusable)\n",
                     wa.ratio, wa.required, wa.fatal) *
            "A continuum orbital goes as r^gamma near the origin and oscillates with a nearly constant amplitude outside " *
            "the atom, so\nthis ratio is 1.27 or below for every case it has been measured on -- 264 of them, Z = 2 to 92, " *
            "0.01 to 3.6e+04\nHartree, l up to 6, bare and screened potentials, all four normalization conventions.  A " *
            "larger value means the radial\nmesh cannot follow the orbital's own oscillation, so the amplitude -- and any " *
            "rate, cross section or coefficient built\nfrom it -- is wrong by about that factor or its square.\n" *
            "NOTHING ELSE CATCHES THIS: the mesh satisfies Continuum.gridConsistency's 15-points-per-oscillation rule, " *
            "which is an\na-priori test of the step against the energy and cannot see the orbital that was produced.\n" *
            ">>> Reduce grid.hp to at most " * @sprintf("%.3e", wavelgth/40.0) * " a.u. (about 40 points per oscillation " *
            "at this\n    energy), or lower the highest continuum energy.\n" )
end


"""
`Continuum.asymptoticFitReport(residual::Float64, required::Float64, fatal::Float64, cOrbital::Orbital)`
    ... formats the measurement behind the residual tier of `Continuum.reportContinuumAmplitude` into the explanatory block that
        both its warning and its refusal print, so that the two cannot drift apart;  a `report::String` is returned.
"""
function asymptoticFitReport(residual::Float64, required::Float64, fatal::Float64, cOrbital::Orbital)

    return( @sprintf("    subshell, energy         = %s,  %.4e Hartree\n", string(cOrbital.subshell), cOrbital.energy) *
            @sprintf("    rms misfit / amplitude   = %.4e   (above %.2f is suspect, above %.2f unusable)\n",
                     residual, required, fatal) *
            "The normalization is fixed by requiring the orbital's asymptotic oscillation to carry a prescribed amplitude, so " *
            "it is\nmeaningless unless the orbital HAS reached that oscillation inside the box.  This number is the rms " *
            "departure of the\nlarge component from the fitted asymptotic form, as a fraction of the amplitude, and it is free " *
            "of the arbitrary scale\nthe Galerkin eigenvector carries.  Measured 03-Nov-2026 on He+: 2.2e-04 at 200 eV to " *
            "1.8e-01 at 0.5 eV for sound\norbitals, against 1.9e-01 to 1.8e+01 for orbitals in boxes far too small for the " *
            "electron's wavelength.\n" *
            "NOTE THAT A LARGER BOX IS THE CURE AND NEVER THE CAUSE HERE, since this measure is a FRACTION of the " *
            "amplitude:\nenlarging the box lets the orbital reach its oscillation and the number falls -- measured 1.8e-01 " *
            "to 1.4e-01 at\n0.5 eV on going from 24 to 40 a.u.  (The superseded `abs(N) > 30` test behaved the opposite way " *
            "and is described in\nthis function's docstring.)\n" *
            ">>> Enlarge rbox until at least one full oscillation of this electron fits outside the atom, or raise the " *
            "energy.\n" )
end


"""
`Continuum.checkNormalizationRadius(pot::Radial.Potential)`
    ... checks that the point at which a continuum orbital's normalization and phase are fixed lies far enough
        outside the atom, and raises with the measurement behind it otherwise; `nothing` is returned.

        WHY THIS EXISTS SEPARATELY FROM `gridConsistency`, which asks a related question and cannot answer this
        one. `gridConsistency` receives only `(maxEnergy, grid)` and therefore knows nothing about the atom, so
        the best it can do is an ABSOLUTE test -- `grid.r[nrContinuum] < 2.0` -- and an absolute length cannot be
        right for every charge: 2 a.u. is deep inside hydrogen and far outside Ne^9+. What matters is the
        DIMENSIONLESS product Zeff * r.

        AND Zeff NEEDS NO ARGUMENT, which is the point of putting the check here. `Radial.Potential` carries
        `Zr`, with `V(r) = -Zr(r)/r`, so `pot.Zr[n]` IS the effective charge the escaping electron still feels at
        radius r[n] -- the residual ion's charge, seen from where the orbital is normalized. Every generator in
        this module that solves in a potential already has it, so all twenty call sites of `gridConsistency`
        across the package are protected without any of them being changed.

        THE CRITERION, measured 01-Sep-2026 on H-like photoionization at Z = 1 and Z = 10, where the two gauges
        are EQUAL for exact wavefunctions so every departure from 1.00000 is error and nothing else:

            Zeff * r    2.12    4.09    4.84    6.50    9.77    10.74   14.70   21.69   48.32
            Cou/Bab     0.004   0.085   4.15    3.78    1.072   0.9996  0.9997  1.00001 1.00000

        The two charges fall on the SAME curve in this variable, which is what identifies it as the right one:
        Zeff*r = 21.69 gives 1.00001 at Z = 1 and 1.00000 at Z = 10. Below 10 the errors are not small -- a
        factor of 200 at the worst point -- and the failure is NON-MONOTONE (0.085 and then 3.78 on a FINER
        grid), which is the signature of the normalization point landing at an arbitrary phase of the oscillation
        rather than of insufficient resolution.

        THE MANY-ELECTRON CASE, measured 02-Sep-2026, which is what the criterion had been waiting for. For a
        many-electron ion the two gauges are NOT equal even for exact wavefunctions, so the H-like error measure
        does not exist; CONVERGENCE does. Neon with a 1s hole -- the system of `testModule_Cascade_StepwiseDecay`,
        i.e. the one case a refusal here would block -- decaying by Auger emission to Ne^2+, total rate against box:

            Zeff * r    6.81       10.46      13.21      18.25      22.13      30.13      44.18      62.04
            rate        9.432734   9.427362   9.427437   9.427417   9.427461   9.427415   9.427445   9.427469   (e-03)

        THE THRESHOLD IS IN THE RIGHT PLACE: the rate settles to 9.42742e-03 within ~5e-6 from Zeff*r = 10.46 on,
        and the only point that misses it is the one BELOW the criterion. **And the local-vs-asymptotic ambiguity
        that held this up turns out not to be material**: the local Zr equals the asymptotic charge (2.000) to 0.1 %
        at every radius where the check is evaluated, and differs (2.099) only at r = 3.2 a.u., which the check
        refuses anyway.

        BUT THE PENALTY FOR VIOLATING IT IS OBSERVABLE-DEPENDENT BY ORDERS OF MAGNITUDE, and that is the reason
        this is still a warning. At Zeff*r = 6.8 the Auger rate above is wrong by 5.7e-4 RELATIVE; at Zeff*r = 4.8
        the H-like photoionization gauge ratio is wrong by a FACTOR OF FOUR. A ratio of two amplitudes exposes the
        continuum phase directly, while a rate built from one matrix element is far more forgiving of it. So the
        number 10 is right for the quantity it was measured on, and very conservative for others -- which is why
        promoting this to an unconditional `error` would refuse calculations whose error is 1e-4.

        A FREE PARTICLE IS EXEMPT, and must be. When the residual system is neutral the electron sees no Coulomb
        tail, Zeff -> 0, and `Zeff * r >= 10` could never be satisfied by any box; there is also nothing for it to
        guard, since a plane wave has no Coulomb phase to get wrong. The check therefore applies only where a
        Coulomb tail exists.

        A NamedTuple is returned, carrying the MEASUREMENT rather than a verdict -- priority item 40, whose rule is
        that a diagnostic hands back a measured quantity with its context so the CALLER judges:
        `(ok, checked, Zeff, r, ZeffTimesR, required, rWanted, note)`.  **`checked` is the field that matters most
        and the reason this is not a bare Bool**: two paths leave the criterion UNTESTED -- a grid too short for
        the probe point to exist, and a free particle with no Coulomb tail -- and both return `ok = true` because
        nothing failed.  NOT CHECKED IS NOT PASSED, and a caller reading only `ok` would take either for a clean
        bill of health.  `ZeffTimesR` against `required` is the quantity to judge on, and `rWanted` is the radius
        that would satisfy it.  Nothing raises: this warns, deliberately, for the reason given above.
"""
function checkNormalizationRadius(pot::Radial.Potential)
    grid = pot.grid;    n = grid.NoPoints - 200
    # NOT CHECKED IS NOT THE SAME AS PASSED, and the returned `checked` field is what keeps the two apart.  Both
    # early exits below leave the criterion UNTESTED -- the first because the grid is too short for the probe
    # point to exist at all, the second because a free particle has no Coulomb tail to guard -- and a caller that
    # read only `ok` would take either for a clean bill of health.  Priority item 40: the diagnostic returns what
    # it measured, with its context, and the caller judges.
    if  n < 1  ||  n > length(pot.Zr)
        return( (ok = true, checked = false, Zeff = NaN, r = NaN, ZeffTimesR = NaN, required = 10.0,
                 rWanted = NaN, note = "not checked: the grid has fewer than 200 points beyond the probe index") )
    end
    Zeff = pot.Zr[n]
    if  Zeff < 0.5
        return( (ok = true, checked = false, Zeff = Zeff, r = grid.r[n], ZeffTimesR = Zeff * grid.r[n],
                 required = 10.0, rWanted = NaN,
                 note = "not checked: Zeff < 0.5 at the probe point, i.e. a free particle with no Coulomb tail") )
    end
    if  Zeff * grid.r[n] < 10.0
        # WARNS RATHER THAN REFUSES, and the distinction is deliberate. Promoting this to an `error` is a
        # ONE-WORD change and is the intended end state; what it waits on is a MANY-ELECTRON case that validates
        # the criterion. The number 10 was measured on hydrogen-like systems, where Zeff is the bare nuclear
        # charge and is CONSTANT in r. In a many-electron ion Zr(r) falls towards the residual charge, so
        # "Zeff" at the normalization point and the asymptotic charge are different numbers, and the data behind
        # the 10 cannot say which of them the criterion should use. Refusing on that basis would be exactly the
        # "a finding from one scan is a hypothesis" mistake this list has twice paid for.

        # WHAT IT COSTS TODAY, measured 01-Sep-2026: as an ERROR this refuses exactly ONE case in the whole test
        # suite -- testModule_Cascade_StepwiseDecay, at Zeff*r = 9.36 against 10, i.e. 6 % short. Repairing that
        # needs rbox 9.5 -> 14.0, a 47 % larger box, and it MOVES THE NUMBERS, so the approved reference would
        # have to be re-approved. Note also that enlarging the box from 9.5 to 11.5 makes the product WORSE
        # (9.36 -> 8.63), because Zr falls faster than r grows over that range; the requirement is really
        # r >= 10/Z_asymptotic, i.e. about 10 a.u. for a singly charged residual ion.
        Defaults.warn(AddWarning(), "Continuum.checkNormalizationRadius(): Zeff * r = " *
                      "$(round(Zeff * grid.r[n], digits=2)) at the normalization point, below the 10 that " *
                      "hydrogen-like measurements require; the continuum phase and normalization may be wrong.")
        sa = @sprintf("r = %.4f a.u., Zeff = %.3f, Zeff*r = %.2f against a required 10.\n",
                      grid.r[n], Zeff, Zeff * grid.r[n])
        sb = @sprintf(">> Enlarge rbox until grid.r[NoPoints-200] exceeds %.3f a.u.\n", 10.0/Zeff)
        printstyled("\n>> WARNING: the continuum normalization point lies close to the atom: " * sa * sb *
                    ">> See the docstring of Continuum.checkNormalizationRadius for the measurement.\n",
                    color=:light_red)
        return( (ok = false, checked = true, Zeff = Zeff, r = grid.r[n], ZeffTimesR = Zeff * grid.r[n],
                 required = 10.0, rWanted = 10.0/Zeff, note = "the normalization point lies too close to the atom") )
    end

    return( (ok = true, checked = true, Zeff = Zeff, r = grid.r[n], ZeffTimesR = Zeff * grid.r[n],
             required = 10.0, rWanted = 10.0/Zeff, note = "") )
end


"""
`Continuum.gridConsistency(maxEnergy::Float64, grid::Radial.Grid; minEnergy::Union{Nothing,Float64}=nothing,
                            Zeff::Union{Nothing,Float64}=nothing)`   
    ... to check the consistency of the given grid with the energies of the required continuum electrons; 
        an error message is issued if the grid.hp = 0.   or  15 * grid.hp < wavelength(maxEnergy)   or  if the grid has 
        less than 600 grid points. The function also returns the recommended grid point where the normalization and phase 
        is to be determined. This number is currently set to nrContinuum = grid.NoPoints - 200  ... to correct for the wrong
        'phase behaviour' at large r-values.

        TWO QUANTITIES ARE TESTED, NOT ONE, AND THEY ARE INDEPENDENT.  The first is the familiar one, the asymptotic step
        against the de Broglie wavelength -- how finely one oscillation is SAMPLED.  The second is how many oscillations
        the grid holds at all, and it is the criterion this function lacked until 06-Nov-2026.  Measured 01-Sep-2026 on H
        1s photoionization at eps = 0.5 Ha, where the two gauges are equal for exact wavefunctions so every departure
        from 1.00000 is error, with grid.hp held FIXED at 0.060 and only the box changed:

            rbox      wavelengths in box    points per oscillation    Cou/Bab
            18.09     2.88                  104.7                     1.067
            30.00     4.77                  104.7                     1.00001

        THE SAMPLING IS IDENTICAL ON BOTH ROWS and one is 6.7 % wrong, so no criterion on grid.hp alone can separate
        them.  The phase and the normalization are read near the outer boundary, and reading them needs several complete
        oscillations to read FROM;  2.88 is not enough and 4.77 is.  Below 4 this function warns.

        BUT THAT SECOND CRITERION APPLIES ONLY TO THE ASYMPTOTIC-SINE CONVENTIONS, and the warning is withheld for
        `CoulombFunctionNorm`, the shipped default.  Measured 06-Nov-2026 on Mg^10+ (1s^2) with a self-consistent
        screened potential and a 10 eV electron:  enlarging rbox from 6.5 to 30.0 a.u., i.e. from 0.89 to 4.09
        oscillations, left the inner phase shift, the normalization constant and the orbital identical to six decimals
        for kappa = -1, 1 and -2 alike.  Exact Coulomb functions hold wherever the potential has reached its -Zbar/r
        tail rather than only asymptotically, and there Zeff*r = 46 at the matching point, so the tail condition is met
        well inside one oscillation.  A measure may only guard the assumption its own routine makes.

        THE TESTS WANT OPPOSITE ENDS OF THE ENERGY RANGE, which is what `minEnergy` is for.  The mesh must sample
        the SHORTEST oscillation, so it is tested at the fastest electron, `maxEnergy`;  the box must CONTAIN the
        longest, so that belongs to the slowest.  A caller that knows its range passes `minEnergy` and is told
        about both ends;  one that passes only a maximum is told about the fast end alone, which is how every
        call site behaved before 06-Nov-2026 and how most still do.

        AND THE TWO ENDS CARRY DIFFERENT SEVERITIES, which is a decision and not an oversight -- the maintainer's,
        06-Nov-2026.  The fast end may be REFUSED, because at a fixed energy "enlarge the box to four wavelengths"
        is an instruction that can be carried out.  The slow end only WARNS, because the wavelength grows as
        1/sqrt(E) without bound:  no box of any size holds one wave at threshold, and several of JAC's own
        calculations place their lowest energy AT threshold deliberately, a rate coefficient being the obvious one.
        A refusal there would forbid the calculation rather than ask for a larger box.  Nor can the offending
        energies be dismissed as negligible:  the colder the plasma, the more the slowest electrons contribute to
        the very average they would be dropped from.  So the slow end reports WHICH PART of the output is affected
        and leaves the run to proceed.  A `minEnergy` of zero or below, and an empty line list arriving as
        `maxEnergy = 0`, are NOT ASKED rather than failed, since there is then no oscillation to count.

        AND A FIXED POINTS-PER-OSCILLATION NUMBER CANNOT BE CALIBRATED, which is why the sampling test refuses at 15 and
        only WARNS at 40 rather than refusing there.  Three convergence studies in this repository demand three widely
        different densities:  34 points converge a 4.66 keV Auger rate to 0.02 % (06-Nov-2026), 104 are not enough for
        neutral H in a small box (01-Sep-2026), and 463 are needed for radiative recombination into H-like Bi at 1 eV
        (17-Sep-2026).  The ordering follows the Sommerfeld parameter eta = Zbar/q -- about 1 for the first two and 306
        for the last -- because the asymptotic wavelength is the LONGEST one in the problem and therefore flatters a
        heavy ion at low energy.  So 40 is the value below which the published measurements show percent-level error and
        is worth saying out loud;  it is not a value above which anything is guaranteed, and only halving grid.hp until
        the answer stops moving establishes that.

        OPEN POINT, stated here rather than emitted as a per-call warning as it was until 12-Aug-2026: the grid
        point at which the normalization and phase are fixed deserves a better choice than the present fixed
        offset. It is a property of the method, not of any particular run.

        MEASURED 01-Sep-2026, and it turns that open point into the sharpest of these checks. The offset is FIXED
        AT 200 GRID POINTS from the end, so where it lands physically depends on hp and NoPoints and not on the
        atom at all. On H-like photoionization, where the two gauges are equal for exact wavefunctions:

            Zeff * r[NoPoints-200]   2.12    4.09   4.84   6.50   9.77    10.74   14.70   21.69   48.32
            Cou/Bab                  0.004   0.085  4.15   3.78   1.072   0.9996  0.9997  1.00001 1.00000

        Z = 1 and Z = 10 fall on the SAME curve in the dimensionless variable Zeff*r, which is what identifies it
        as the right one, and the errors below 10 are not small -- a factor of 200 at the worst point. THE
        FAILURE IS NON-MONOTONE (0.085 and then 3.78 on a FINER grid), which is the signature of the
        normalization point landing at an arbitrary phase of the oscillation and NOT of insufficient resolution.

        Hence the optional `Zeff`: when it is given, this function refuses a grid with Zeff*r < 10. It is optional
        only because the existing callers pass `(maxEnergy, grid)` and none of them can be changed from here; each
        should be migrated to supply it, one module at a time. Until a caller does, the absolute `2.0` test is all
        that guards it, and the table above shows that 2.0 establishes nothing: at Z = 1 it admits grids that are
        wrong by two orders of magnitude.
"""
function gridConsistency(maxEnergy::Float64, grid::Radial.Grid; minEnergy::Union{Nothing,Float64}=nothing,
                            Zeff::Union{Nothing,Float64}=nothing)
    wavenb      = sqrt( 2maxEnergy + maxEnergy * Defaults.getDefaults("alpha")^2 )
    wavelgth    = 2pi / wavenb
    nrContinuum = grid.NoPoints - 200
    # The tests below want OPPOSITE extrema of the energy range, so the slow end carries its own wavelength:
    # the mesh must sample the SHORTEST oscillation and the box must contain the LONGEST.  A caller that passes
    # no minimum reproduces the behaviour of every version before 06-Nov-2026 exactly.
    slowEnergy   = isnothing(minEnergy) ? maxEnergy : minEnergy
    slowWavenb   = slowEnergy  > 0.  ?  sqrt( 2slowEnergy + slowEnergy * Defaults.getDefaults("alpha")^2 )  :  0.
    slowWavelgth = slowWavenb  > 0.  ?  2pi / slowWavenb  :  0.
    
    ptsPerOsc = wavelgth / grid.hp

    if      grid.hp == 0.
        error("\n\nContinuum.gridConsistency():  STOP -- this grid has grid.hp = 0. and therefore no asymptotic region " *
                "with a constant step size.\n" *
                "A continuum orbital oscillates out to the box boundary and cannot be represented on such a grid.\n" *
                ">>> Build the grid with an explicit hp, e.g.  Radial.Grid(Radial.Grid(true), rnt=4.0e-6, h=5.0e-2, " *
                "hp=2.5e-2, rbox=30.0),\n    or let Radial.generateGrid(grid; maximumFreeElectronEnergy = ...) choose hp " *
                "from the shortest de Broglie wavelength.\n")
    elseif  15 * grid.hp > wavelgth
        error("\n\nContinuum.gridConsistency():  STOP -- the radial grid is too coarse for the requested continuum energy.\n" *
                @sprintf("    continuum energy         = %.6e  Hartree\n", maxEnergy) *
                @sprintf("    de Broglie wavelength    = %.6e  a.u.\n",    wavelgth)  *
                @sprintf("    asymptotic step grid.hp  = %.6e  a.u.\n",    grid.hp)   *
                @sprintf("    points per oscillation   = %.2f      (a minimum of 15 is required)\n", ptsPerOsc) *
                "The continuum orbital would be sampled too sparsely to resolve its own oscillation. This does NOT fail " *
                "loudly by itself:\ncross sections first drift by tens of percent while still looking entirely plausible, " *
                "and then diverge by many orders of\nmagnitude once the grid can no longer follow the phase at all. " *
                "The computation is therefore stopped here rather than\nreturning numbers that cannot be trusted.\n" *
                @sprintf(">>> Reduce grid.hp to at most %.3e a.u., or lower the highest continuum energy;\n", wavelgth/15) *
                "    Radial.generateGrid(grid; maximumFreeElectronEnergy = ...) sets hp from this same wavelength.\n" )
    elseif  grid.NoPoints < 600
        error("\n\nContinuum.gridConsistency():  STOP -- too few grid points for continuum processes; " *
                "grid.NoPoints = $(grid.NoPoints) < 600.\n" *
                "The normalization and phase of the continuum orbital are determined near the outer boundary and need a " *
                "well-resolved\nasymptotic region.\n>>> Enlarge the grid (smaller h/hp, or a larger rbox).\n")
    elseif  grid.r[nrContinuum] < 2.0
        error("\n\nContinuum.gridConsistency():  STOP -- the grid does not extend far enough for continuum processes; " *
                @sprintf("grid.r[nrContinuum] = %.4e < 2.\n", grid.r[nrContinuum]) *
                "The normalization point (grid.NoPoints-200) still lies inside the atomic core region, where the orbital " *
                "has not yet\nreached its asymptotic form.\n>>> Enlarge rbox.\n")
    elseif  !isnothing(Zeff)  &&  Zeff * grid.r[nrContinuum] < 10.0
        # THE 2.0 ABOVE IS AN ABSOLUTE LENGTH AND SO CANNOT BE RIGHT FOR EVERY CHARGE: the atom's own size scales
        # as 1/Zeff, so a normalization point 2 a.u. out is deep inside hydrogen and far outside Ne^9+.  What
        # matters is the DIMENSIONLESS Zeff * r, and this clause is the one with the measurement behind it.
        error("\n\nContinuum.gridConsistency():  STOP -- the normalization point lies too close to the atom.\n" *
                @sprintf("    normalization radius     = %.4f a.u.   (grid point NoPoints-200)\n", grid.r[nrContinuum]) *
                @sprintf("    Zeff * that radius       = %.2f       (a minimum of 10 is required)\n",
                         Zeff * grid.r[nrContinuum]) *
                "The normalization and phase of the continuum orbital are fixed at that ONE point, so it must lie well " *
                "outside\nthe bound-state charge cloud, whose extent scales as 1/Zeff.  Measured 01-Sep-2026 on H-like " *
                "photoionization,\nwhere the two gauges are EQUAL for exact wavefunctions so every departure from " *
                "1.00000 is error:\n\n" *
                "    Zeff*r   2.12    4.09    4.84    6.50    9.77    10.74   14.70   21.69   48.32\n" *
                "    Cou/Bab  0.004   0.085   4.15    3.78    1.072   0.9996  0.9997  1.00001 1.00000\n\n" *
                "The two charges Z = 1 and Z = 10 fall on the SAME curve in this variable, which is why it is the " *
                "right one:\nZeff*r = 21.69 gives 1.00001 at Z = 1 and 1.00000 at Z = 10.  Note also that the failure " *
                "is NON-MONOTONE\n(0.085 then 3.78 on a FINER grid), which is the signature of a normalization point " *
                "landing at an arbitrary\nplace in the oscillation rather than of insufficient resolution.\n" *
                @sprintf(">>> Enlarge rbox until grid.r[NoPoints-200] exceeds %.3f a.u.\n", 10.0/Zeff))
    end

    # TWO WARNING TIERS SIT ABOVE THE REFUSALS, and they are warnings rather than refusals on purpose.  The refusals
    # above mark grids on which nothing can be computed;  these two mark grids on which the published measurements
    # show percent-level error, which is a different statement and one the user may legitimately accept.  Raising
    # either refusal to the warning level would reject grids in current use -- and, for the sampling test, would
    # still not be sufficient, since 104 points per oscillation is measurably too few for neutral H in a small box.
    #   `maxlog=1` is deliberate:  the cascade modules call this once per decay step with the SAME grid, so an
    # unguarded warning would print dozens of identical lines and be read past.  The AddWarning() record is what
    # carries it into the summary file, where it cannot be scrolled away.
    # THE DENSITY IS COUNTED HERE RATHER THAN DIVIDED, for the reason given at Continuum.checkMatchingRadius:  the
    # mesh is a B-spline knot sequence sampled at Gauss-Legendre points, so no single step represents it.  The mean
    # step over a wavelength is 0.87 to 1.00 times grid.hp on the systems measured 06-Nov-2026, so `wavelength/hp`
    # and the counted density agree to a few percent and the refusal above needs no change on this account.
    #   AND THIS TIER WARNS RATHER THAN REFUSES BECAUSE THE QUANTITY IS NOT A PREDICTOR.  Over four systems and
    # thirteen meshes each, the error in the normalization constant is NON-MONOTONE in the density:  two meshes at
    # 21 counted points differ by 3.3 %, and uranium at 15 points is exact while at 20 points it is 3.4 % out.  So a
    # low density is a reason to doubt a number, not a proof that it is wrong -- and a high one is no guarantee.  The
    # reliable signal is the posterior verdict of Continuum.checkMatchingRadius, which this cannot anticipate.
    # THE BOX IS MEASURED IN WAVELENGTHS OF ITS OWN FULL EXTENT, not of the normalization radius.  The measurement
    # behind this test varied grid.rbox, so rbox is the quantity it calibrates;  r[NoPoints-200] sits short of the
    # box end and testing THAT against the same number is a harsher criterion than anything measured.
    #   AND IT IS A WARNING, NOT A REFUSAL, DELIBERATELY.  Written first as a refusal at 2 oscillations on
    # 06-Nov-2026, it stopped four tests that carry approved references -- PhotoRecombination at 0.63 oscillations
    # of the normalization radius, Cascade_StepwiseDecay at 0.51, PhotoIonizationInterference at 1.47 -- which says
    # that a substantial part of JAC's low-energy continuum work sits below that line.  Whether JAC should refuse it
    # is a decision about what the code accepts, and it rests on ONE measurement;  that makes it a question for the
    # maintainer rather than a stop to be introduced here.
    # A non-positive energy has no oscillation to count, so the box tier is NOT ASKED rather than failed:
    # a caller whose energy range starts AT threshold, and one whose line list came out empty, both reach
    # here with a zero energy and must not be refused for it.
    boxWaves  = wavenb       > 0.  ?  grid.r[grid.NoPoints] / wavelgth      :  Inf
    slowWaves = slowWavelgth > 0.  ?  grid.r[grid.NoPoints] / slowWavelgth  :  Inf
    iOne     = something( findfirst(k -> grid.r[k] >= grid.r[nrContinuum] - wavelgth, 1:nrContinuum), 1 )
    ptsCount = Float64(nrContinuum - iOne + 1)
    if  ptsCount < 40.0
        Defaults.warn(AddWarning(), "Continuum.gridConsistency(): " *
                      @sprintf("%.1f", ptsCount) * " counted points per oscillation at " *
                      @sprintf("%.4e", maxEnergy) * " Ha, below the 40 that measurements show to be sound; " *
                      "continuum amplitudes may carry percent-level error.")
        @warn("Continuum.gridConsistency(): the mesh holds " * @sprintf("%.0f", ptsCount) *
              " points across one oscillation at the highest continuum energy " * @sprintf("%.4e", maxEnergy) *
              " Ha, below the 40 that measurements show to be sound.  Measured 06-Nov-2026 on three systems with " *
              "the box held fixed, the normalization constant is 4 to 6 % out at 13 counted points and about 1 % " *
              "out at 15, while a 4.66 keV Auger rate converges to 0.02 % by 34 points of wavelength/grid.hp.  " *
              "Reduce grid.hp by a factor of about " * @sprintf("%.1f", 40.0/max(ptsCount,1.0)) * ", from " *
              @sprintf("%.3e", grid.hp) * " to about " * @sprintf("%.3e", grid.hp*max(ptsCount,1.0)/40.0) *
              " a.u., or accept the error knowingly.  Note that 40 is NOT a guarantee -- 104 points per " *
              "oscillation is measurably too few for neutral hydrogen in a small box -- and that halving grid.hp " *
              "until the answer stops moving is the only real check.  For reference, wavelength/grid.hp here is " *
              @sprintf("%.1f", ptsPerOsc) * ", which is the quantity the refusal above tests.", maxlog=3)
    end
    # AND THE BOX CRITERION BELONGS TO THE ASYMPTOTIC-SINE CONVENTIONS ONLY.  This is a CORRECTION to the first form
    # of this warning rather than a refinement of it.  MEASURED 06-Nov-2026 on Mg^10+ (1s^2) with a self-consistent
    # SCREENED potential, a 10 eV electron and kappa = -1, 1, -2:  enlarging rbox from 6.5 to 30.0 a.u., i.e. from
    # 0.89 to 4.09 oscillations, left the inner phase shift, the normalization constant AND the orbital itself
    # identical to six decimals in every case.
    #   THE REASON IS THE CONVENTION.  normalizeOrbitalCoulombFunctions matches EXACT Coulomb functions, which are
    # valid wherever the potential has reached its -Zbar/r tail and not merely asymptotically;  there Zeff*r = 46 at
    # the matching point, so the tail condition is amply met INSIDE ONE OSCILLATION.  "A few oscillations in the box"
    # is the requirement of a routine that FITS AN ASYMPTOTIC SINE, which is what normalizeOrbitalPureSine does.  A
    # measure may only guard the assumption its own routine makes, so the warning is raised for those conventions
    # and withheld for this one.
    #   NOTE ALSO, because it bears on whether this test adds anything at all:  the 01-Sep-2026 box measurement it
    # was built on is better explained by the variable checkNormalizationRadius ALREADY uses.  Varying rbox at fixed
    # grid.hp moves r[NoPoints-200] and hence Zeff*r, and that docstring records Z = 1 and Z = 10 falling on the SAME
    # curve in Zeff*r -- which is the evidence that Zeff*r is the right parameterisation and oscillations-in-box a
    # disguised version of it.  This tier is therefore kept narrow deliberately.
    # AND THIS TIER IS ONLY AS GOOD AS THE ENERGY IT IS HANDED, which is a real limitation and not a quibble.
    # `maxEnergy` is the right extremum for the MESH tier above -- the FASTEST electron has the shortest wavelength
    # and is the hardest to resolve -- and it is the WRONG one here, because the SLOWEST electron has the longest
    # wavelength and is the one that may not fit in the box.  The two tiers want opposite extrema and this function
    # receives only one.
    #   SO THE BOX TIER BITES ONLY WHERE THE CALLER HAPPENS TO PASS A SLOW ENERGY, and whether it does is a property
    # of the caller rather than of the grid.  `Cascade-inc-stepwise-decay` calls this per STEP with that step's own
    # mean energy, which is how a 0.1 Ha Auger channel of neon was caught;  a photoionization run at a fixed photon
    # energy passes one fast energy and a slow partner hidden in an energy SHARING -- photo-double-ionization is the
    # case -- is never seen at all.  Measured 06-Nov-2026: `example-Fb.jl` branch (a), the same scheme and a 10 a.u.
    # box, only WARNS because Mg's slowest step sits at 0.997 Ha, while the neon cascade of
    # `testModule_Cascade_StepwiseDecay` refuses at 0.1 Ha in a 9.5 a.u. box.  The difference is the level structure,
    # not the grid.
    #   HENCE THE OPTIONAL `minEnergy`, added 06-Nov-2026, which is the fix for exactly that:  when it is given, the
    # box tier is evaluated at the SLOWEST electron of the run while the mesh tier stays at the fastest.  It is optional
    # for the same reason `Zeff` below is -- 27 call sites in 21 modules pass only a maximum and cannot be changed from
    # here, so they migrate one module at a time.  UNTIL A CALLER PASSES IT, READ A SILENT BOX TIER AS "NOT ASKED"
    # RATHER THAN AS "CHECKED".
    #   `module-ImpactExcitation.jl` IS WHERE TO LOOK FIRST:  it already computes the minimum and already tests it, in
    # a module-local `checkLowEnergyGrid` applying the same `lambda > rbox` criterion -- but only WARNS, and only there.
    # Its warning predates the measurement behind the refusal below.
    #   THE MIGRATION IS THEREFORE SAFE AND MECHANICAL, now that the slow end only warns:  a module passes the
    # minimum it already has beside the maximum, and no run that works today can be stopped by it.  That was NOT
    # true of the first design, in which `minEnergy` fed the refusal:  `Cascade-inc-photorecombination` builds a
    # Gauss-Legendre energy grid from 0.0 upward and `ImpactExcitation` says the same of `calcRateCoefficient =
    # true`, so the honest minimum would have refused nearly every rate-coefficient run under a sine convention.
    #   AND HERE IS WHY THE SLOW END IS WORTH REPORTING AT ALL, measured 06-Nov-2026 rather than argued.
    # `examples/example-Dd.jl` branch (b) -- bare Z = 18, K-shell radiative recombination, `rbox` = 10 a.u.,
    # electron energies 1, 10, 100, 1000 and 10000 eV -- hands this function its MAXIMUM, 367 Ha, where the
    # box holds 43 waves.  Every tier passes.  The run then dies inside `normalizeOrbitalPureSine` on the
    # asymptotic-form residual, with nothing pointing at the energy that caused it.  The box holds 0.43 waves
    # at 1 eV and 1.37 at 10 eV;  dropping those two points makes the branch complete, keeping either makes it
    # raise, and from 100 eV up (4.34 waves) it is fine.  So the slow end was the whole story and the box
    # tier, asked about the fast end, had nothing to say.  The branch carries a `Last successful` date of
    # 28-Aug-2026 QUOTING ITS 1 eV CROSS SECTION, which the residual guard now says was never trustworthy.
    #   NOTE ALSO, because two instruments agreeing is worth more than either alone:  `CONT_BOX_WAVES_WARN`
    # was calibrated at 4 waves on a BARE POINT PROTON, knowing nothing of this example, and the independent
    # asymptotic-form residual starts refusing between 1.37 and 4.34 waves on it.
    if  boxWaves < CONT_BOX_WAVES_FATAL  &&  Defaults.GBL_CONT_NORMALIZATION != CoulombFunctionNorm()
        error("Continuum.gridConsistency(): STOP -- this box holds only " * @sprintf("%.2f", boxWaves) *
              " oscillations of the continuum orbital at " * @sprintf("%.4e", maxEnergy) * " Ha, i.e. NOT ONE " *
              "COMPLETE WAVE, and the selected normalization fixes the orbital's scale by fitting an asymptotic " *
              "sine.  There is no complete oscillation inside the box to fit it to.  Measured 06-Nov-2026 on a " *
              "bare point proton, where the exact solution is known: below one wave the normalized orbital is " *
              "15 to 75 % wrong, which squares into a cross section wrong by a factor of 1.3 to 3, while the " *
              "Coulomb-function convention reproduces itself to 7e-8 over the same boxes.  THE FIT RESIDUAL DOES " *
              "NOT WARN YOU HERE: it reads 0.04 at the 75 % error, below its own 0.20 limit, because with less " *
              "than one wave the fit has almost nothing to satisfy.  REMEDY: enlarge grid.rbox to at least " *
              @sprintf("%.2f", CONT_BOX_WAVES_WARN*wavelgth) * " a.u. (" * @sprintf("%.0f", CONT_BOX_WAVES_WARN) *
              " waves at this energy; " * @sprintf("%.2f", wavelgth) * " a.u. would be the bare minimum of one), " *
              "or lower the continuum energy, or select " *
              "setDefaults(\"method: normalization, Coulomb functions\"), which matches exact Coulomb functions " *
              "instead of an asymptotic sine and was measured insensitive to the box.")
    end
    if  boxWaves < CONT_BOX_WAVES_WARN  &&  Defaults.GBL_CONT_NORMALIZATION != CoulombFunctionNorm()
        Defaults.warn(AddWarning(), "Continuum.gridConsistency(): the box holds only " *
                      @sprintf("%.2f", boxWaves) * " oscillations at " * @sprintf("%.4e", maxEnergy) *
                      " Ha, below the 4 above which the box can be ruled out as a cause; the continuum phase and " *
                      "normalization may be wrong by several percent to 20 %.")
        @warn("Continuum.gridConsistency(): the box reaches " * @sprintf("%.2f", boxWaves) *
              " oscillations of the continuum orbital at " * @sprintf("%.4e", maxEnergy) * " Ha.  THIS IS " *
              "INDEPENDENT OF THE STEP SIZE and is not implied by it: the phase and the normalization are read " *
              "near the outer boundary, so there must be complete oscillations out there to read them FROM, " *
              "however finely each one is sampled.  Measured 06-Nov-2026 on a bare point proton, at a FIXED " *
              "energy so that only the box varied: the normalized orbital came out 26 % wrong at 0.97 waves and " *
              "10 to 19 % wrong between 1.2 and 1.8, against 7e-8 for the Coulomb-function convention over the " *
              "same boxes.  Above 4 waves this experiment can no longer attribute the residual few percent to " *
              "the box, which is why 4 is the level here -- it is where the box stops being demonstrably the " *
              "cause, NOT a level above which the answer is known to be converged.  Enlarge grid.rbox to about " *
              @sprintf("%.2f", CONT_BOX_WAVES_WARN*wavelgth) * " a.u., or lower the continuum energy.  This " *
              "applies to the asymptotic-sine conventions;  " *
              "setDefaults(\"method: normalization, Coulomb functions\") matches " *
              "exact Coulomb functions instead and was measured insensitive to the box.", maxlog=3)
    end
    # THE SLOW END WARNS AND NEVER REFUSES, and the asymmetry with the two tiers above is the whole point --
    # the maintainer's decision, 06-Nov-2026.  A refusal has to name a remedy that can be carried out.  At the
    # PRINCIPAL energy it can:  the wavelength is fixed, so `enlarge grid.rbox to 4 lambda` is a finite
    # instruction.  At the SLOW end it cannot:  the wavelength grows as 1/sqrt(E) without bound, so NO FINITE
    # BOX holds one wave at threshold, and several of JAC's own calculations put their lowest energy AT zero on
    # purpose -- a rate coefficient needs the slow electrons.  Refusing there would not ask for a larger box, it
    # would forbid the calculation.
    #   AND THE OBVIOUS ESCAPE DOES NOT EXIST:  one cannot drop the offending energies as negligible, because
    # the colder the plasma the more the slowest electrons carry -- in a cold case they dominate the very
    # average they would be dropped from.  So what is owed to the caller is the TRUTH ABOUT WHICH PART OF THE
    # OUTPUT IS AFFECTED, not a stop.  A caller holding a list of energies should additionally say HOW MANY of
    # them fall below the line;  this function sees only the extremes and reports the extreme.
    if  !isnothing(minEnergy)  &&  slowWaves < CONT_BOX_WAVES_WARN  &&  slowWaves < boxWaves  &&
        Defaults.GBL_CONT_NORMALIZATION != CoulombFunctionNorm()
        if  slowWaves < CONT_BOX_WAVES_FATAL
            sa = "NOT ONE COMPLETE WAVE, so its normalization is read from a wave the box never shows;  " *
                 "measured 15 to 75 % wrong on a bare point proton, a factor 1.3 to 3 in a cross section"
        else
            sa = "fewer than the 4 waves above which the box can be ruled out as a cause;  percent to " *
                 "20 % in the phase and the normalization"
        end
        Defaults.warn(AddWarning(), "Continuum.gridConsistency(): THE SLOW END OF THIS RUN IS NOT "             *
                      "REPRESENTED. The slowest continuum electron, at " * @sprintf("%.4e", slowEnergy)        *
                      " Ha, reaches " * @sprintf("%.2f", slowWaves) * " oscillations in this box -- " * sa     *
                      ". The faster energies of the same run are unaffected.")
        @warn("Continuum.gridConsistency(): THE SLOW END OF THIS RUN IS NOT REPRESENTED.  The slowest "         *
              "continuum electron, at " * @sprintf("%.4e", slowEnergy) * " Ha, reaches only "                  *
              @sprintf("%.2f", slowWaves) * " oscillations inside a box of "                                   *
              @sprintf("%.2f", grid.r[grid.NoPoints]) * " a.u., where the selected normalization fits an "      *
              "asymptotic sine -- " * sa * ".  The faster energies of the same run are unaffected, so this "   *
              "is a statement about PART of the output.  THIS IS A WARNING AND NOT A REFUSAL BY DESIGN: the "   *
              "wavelength grows without bound as the energy falls, so no box of any size holds one wave at "    *
              "threshold, and a calculation that deliberately starts at threshold -- a rate coefficient -- "    *
              "cannot be given a box that would satisfy it.  WHAT DOES HELP: "                                  *
              "setDefaults(\"method: normalization, Coulomb functions\"), which matches exact Coulomb "        *
              "functions rather than an asymptotic sine and was measured insensitive to the box, or raising "   *
              "the lowest energy of the run if the slow end is not needed.", maxlog=3)
    end

    return( nrContinuum )
end


"""
`Continuum.asymptoticCharge(Zbar::Float64)`
    ... the residual charge that the ASYMPTOTIC FORM of the computed orbital actually carries, which is not always the charge in
        the potential:  a FREE-wave solution method produces an orbital with no logarithmic Coulomb phase at all, whatever the
        potential says.  `Zbar` is returned unchanged for the solution methods that solve in the potential, and 0.0 for those
        that do not, so that `Continuum.asymptoticAmplitude` is handed the form the orbital really has.  A `Float64` is returned.

        MEASURED 03-Nov-2026 on one He+ orbital at 16.7 eV, as the fit residual with and without the logarithmic term:

            solution method           with Zbar     with 0
            ContBessel                7.83e-02    2.49e-05
            ContSine                  7.83e-02    2.49e-05
            AsymptoticCoulomb         6.96e-16    7.85e-02

        -- so each solver is described by exactly one of the two forms, and `AsymptoticCoulomb` is fitted to MACHINE PRECISION
        by the right one.  Taking the phase model from the potential alone is therefore wrong for two of the five methods, and
        the error does not announce itself:  it appears as a residual a hundred times too large on a perfectly sound orbital.
"""
function asymptoticCharge(Zbar::Float64)

    return( Defaults.GBL_CONT_SOLUTION == ContBessel()  ||  Defaults.GBL_CONT_SOLUTION == ContSine()  ?  0.0  :  Zbar )
end


"""
`Continuum.asymptoticAmplitude(P::Array{Float64,1}, grid::Radial.Grid, mtp::Int64, energy::Float64, Zbar::Float64)`
    ... to determine the asymptotic amplitude and phase of a continuum orbital's large component by FITTING it to its known
        asymptotic form over the outer third of the mesh, instead of reading a maximum off it.  Asymptotically
        `P(r) --> A sin(theta(r) + delta)` with `theta(r) = q*r + y*ln(2*q*r)`, where q is the relativistic wavenumber and y the
        Coulomb parameter formed from the residual charge Zbar that the outgoing electron still feels;  every CONSTANT phase --
        centrifugal, Coulomb gamma-function, relativistic -- is absorbed into delta and so need not be known here.  Writing the
        form as `P = a*sin(theta) + b*cos(theta)` makes the extraction a two-parameter LINEAR least squares, solved below through
        its 2x2 normal equations.  A named tuple `(A, delta, residual, nFitted)` is returned.

        Why a fit and not a window.  Reading a maximum requires a crest inside the window, loses up to (q*hp)^2/8 to discrete
        sampling, and must fail once half a wavelength exceeds the box, however the window is sized;  the fit requires no crest,
        has no node to divide by, averages the discretization away and carries the logarithmic Coulomb phase exactly.  Measured
        03-Nov-2026 against `Continuum.generateOrbitalPureSine`, whose amplitude is `sqrt(2/(pi*q))` in closed form, the fit is
        exact to five or six digits from 68 eV down to 0.054 eV, whereas a fixed 300-point window returned 0.4462 of the true
        amplitude at 0.136 eV -- a normalization 2.24 times too large -- while being exact at 0.272 eV.  That the error is
        ERRATIC rather than monotone in the energy is the signature of the defect, and is why it survived so long.

        `residual` is the rms misfit as a FRACTION of the amplitude, so it is free of the arbitrary scale that the Galerkin
        eigenvector carries.  It measures whether the orbital has actually reached the form that was fitted:  measured on He+ it
        falls from 7.0e-2 at 1.449 eV to 1.5e-3 at 120 eV, and fitting the same Coulomb wave with a FREE phase model (Zbar = 0)
        raises it by a factor of 5 to 20 at every energy while biasing the amplitude by up to 14 %.  It is therefore the
        scale-free quantity on which a refusal should rest.
"""
function asymptoticAmplitude(P::Array{Float64,1}, grid::Radial.Grid, mtp::Int64, energy::Float64, Zbar::Float64)
    wc   = Defaults.getDefaults("speed of light: c")
    q    = sqrt( 2*energy + energy/(wc*wc) )
    y    = Zbar * (energy + wc^2) / (wc^2 * q)
    #  the outer third of the mesh, which keeps the inner, non-asymptotic region out of the fit
    iLow = max(2, mtp - div(mtp, 3));     nFitted = mtp - iLow + 1
    if  energy <= 0.0   ||   nFitted < 10   ||   grid.r[iLow] <= 0.0
        return( (A = 0.0, delta = 0.0, residual = Inf, nFitted = nFitted) )
    end

    Saa = 0.0;   Sab = 0.0;   Sbb = 0.0;   va = 0.0;   vb = 0.0
    for  i = iLow:mtp
        theta = q*grid.r[i] + y*log(2q*grid.r[i]);    sn = sin(theta);    cs = cos(theta)
        Saa = Saa + sn*sn;      Sab = Sab + sn*cs;    Sbb = Sbb + cs*cs
        va  = va  + P[i]*sn;    vb  = vb  + P[i]*cs
    end
    det = Saa*Sbb - Sab*Sab
    if  abs(det) < 1.0e-30   return( (A = 0.0, delta = 0.0, residual = Inf, nFitted = nFitted) )   end
    a = (Sbb*va - Sab*vb) / det;    b = (Saa*vb - Sab*va) / det;    A = sqrt(a*a + b*b)
    if  A <= 0.0             return( (A = 0.0, delta = 0.0, residual = Inf, nFitted = nFitted) )   end

    wa = 0.0
    for  i = iLow:mtp
        theta = q*grid.r[i] + y*log(2q*grid.r[i])
        wa    = wa + (P[i] - a*sin(theta) - b*cos(theta))^2
    end

    return( (A = A, delta = atan(b, a), residual = sqrt(wa/nFitted)/A, nFitted = nFitted) )
end


"""
`Continuum.normalizeOrbitalPureSine(cOrbital::Orbital, grid::Radial.Grid, settings::Continuum.Settings; Zbar::Float64=0.0)`
    ... to normalize the given continuum orbital with regard to a (asymptotic) pure sine-function.
        An (on-energy-scale-normalized) orbital::Orbital and its relative phase phi w.r.t.  sin(kr + phi) is returned.

        + Zbar  ::Float64   ... the residual charge the outgoing electron still feels, which enters only the PHASE that
                                `Continuum.asymptoticAmplitude` fits;  the target amplitude keeps this convention's own
                                `sqrt(2/(pi*q))`.  It defaults to 0.0, the genuinely free wave, so that a caller holding
                                only a grid -- as the test suite does -- gets the free-wave model;  the three callers in
                                this module all have a potential and pass `Radial.determineZbar(pot)`.  Measured
                                03-Nov-2026, fitting a Coulomb wave with the free model biases the amplitude by up to
                                14 % at 1.449 eV and raises the fit residual by a factor of 5 to 20.
"""
function normalizeOrbitalPureSine(cOrbital::Orbital, grid::Radial.Grid, settings::Continuum.Settings; Zbar::Float64=0.0) 
    mtp = size( cOrbital.P, 1);    q = sqrt( 2*cOrbital.energy );        l = Basics.subshell_l(cOrbital.subshell)

    PPprime = cOrbital.P[mtp] / cOrbital.Pprime[mtp];    kr  = q * grid.r[mtp];   at = atan( PPprime * q )
    phi = at - kr + l*pi/2; 
    phi = rem(phi+1000pi, pi)    # to bring phi in the interval 0. <= phi < pi
    # For normalization, FIT the asymptotic form rather than read a maximum off it;  see Continuum.asymptoticAmplitude
    # for why no window can do this job at every energy.  Zbar is the residual charge the outgoing electron still feels:
    # it enters only the fitted PHASE, while the target amplitude keeps this convention's own non-relativistic q.
    wa = Continuum.asymptoticAmplitude(cOrbital.P, grid, mtp, cOrbital.energy, Continuum.asymptoticCharge(Zbar))
    A  = wa.A;    N   = sqrt(2/(pi*q)) / A
    
    println(">> Pure-sine normalized continuum orbital with normalization constant N=" * @sprintf("%.4e",N) *
            " and phase phi=" * @sprintf("%.4e",phi) * " at r=" * @sprintf("%.4e",grid.r[mtp])   * " a.u." )

    
    P = N .* cOrbital.P;   Q = N .* cOrbital.Q;   Pprime = N .* cOrbital.Pprime;   Qprime = N .* cOrbital.Qprime;   
    Continuum.reportContinuumAmplitude(P, cOrbital, grid, "PureSine"; residual = wa.residual)
    newOrbital = Orbital( cOrbital.subshell, cOrbital.isBound, cOrbital.useStandardGrid, cOrbital.energy, 
                            P, Q, Pprime, Qprime, cOrbital.grid)

    normF = N * wa.A   # the fitted asymptotic amplitude of the P component after scaling, as a measure of the amplitude

    return( newOrbital, phi, normF )
end


"""
`Continuum.normalizeOrbitalCoulombSine(cOrbital::Orbital, pot::Radial.Potential, settings::Continuum.Settings)`   
    ... to normalize the given continuum orbital with regard to a (asymptotic) pure sine-function. 
        An (on-energy-scale-normalized) orbital::Orbital is returned.
"""
function normalizeOrbitalCoulombSine(cOrbital::Orbital, pot::Radial.Potential, settings::Continuum.Settings) 
    mtp = size( cOrbital.P, 1);   energy = cOrbital.energy    
    
    if  false  # use non-relativistic Coulomb-sine wave representation
        wc   = Defaults.getDefaults("speed of light: c");    q = sqrt( 2*cOrbital.energy );        l = Basics.subshell_l(cOrbital.subshell)
        Zbar = Radial.determineZbar(pot);    
        PPprime = cOrbital.P[mtp] / cOrbital.Pprime[mtp];    kr  = q * pot.grid.r[mtp];   at = atan( PPprime * q )
        sigmac  = SpecialFunctions.gamma(l + 1 - im*Zbar/q)
        phi = at - kr + l*pi/2 - Zbar/q * log(2*q*pot.grid.r[mtp]) - atan(sigmac.im / sigmac.re); 
        phi = rem(phi+1000pi, pi)    # to bring phi in the interval 0. <= phi < pi
        # For normalization, determine the maximum of the large component within the last 300 point
        A = maximum( abs.(cOrbital.P[end-300:end]) );    N   = sqrt(2/(pi*q)) / A
        
        println(">> (Non-relativistic) Coulomb-sine normalized continuum orbital with normalization constant N=" * @sprintf("%.4e",N) *
                " and phase phi=" * @sprintf("%.4e",phi) * " at r=" * @sprintf("%.4e",pot.grid.r[mtp])   * " a.u." )

    else      # use relativistic Coulomb-sine wave representation
        wc   = Defaults.getDefaults("speed of light: c");    q = sqrt( 2*energy + energy/(wc*wc));      kappa = cOrbital.subshell.kappa
        Zbar = Radial.determineZbar(pot);    y = Zbar * (energy + wc^2) / (wc^2 * q);   gammaBar = sqrt(kappa^2 - Zbar^2 / wc^2)
        eta  = log( - (kappa - im*y) / (energy + wc^2) / (gammaBar + im*y)) / (2*im)
        NP   = sqrt( (energy + 2*wc^2)/(pi*wc^2*q) );   NQ = -sqrt( energy / (pi*wc^2*q) )
        
        # Select the 'ratio' on which the thetaprimephase is determined
        if  false   PprimeP = cOrbital.Pprime[mtp] / cOrbital.P[mtp];        at = atan( - PprimeP / thetaprime );   son = "P'/P"    
        else        PprimeP = cOrbital.Q[mtp] * NP / cOrbital.P[mtp] / NQ;   at = atan( PprimeP );                  son = "Q/P"     end

        argG  = SpecialFunctions.gamma(gammaBar + im*y)
        theta = q * pot.grid.r[mtp]  +  y * log(2q * pot.grid.r[mtp])  - atan(argG.im / argG.re)  -  pi*gammaBar/2 #  +  eta.re
        phi = at - theta
        phi = rem(phi+1000pi, pi)    # to bring phi in the interval 0. <= phi < pi
        #  A windowless form was tried here and REFUSED by measurement, 03-Nov-2026:  with phi = at - theta it
        #  reduces to A = sqrt(P^2 + (Q*NP/NQ)^2), the envelope built from the two components of the spinor, which
        #  are in quadrature only asymptotically.  Calibrated against this window ABOVE 14 eV, where the window is
        #  known to read the true amplitude, it missed by 1.4 to 4.2 % and non-monotonically, because at the radii
        #  used here -- about 22 a.u. -- the Coulomb tail has not yet left that quadrature exact.  The FIT below
        #  needs no such assumption:  it absorbs the phase as a free parameter and reads the amplitude from hundreds
        #  of points, which is why it succeeds where that one-point form failed.
        wa = Continuum.asymptoticAmplitude(cOrbital.P, pot.grid, mtp, energy, Continuum.asymptoticCharge(Zbar))
        A  = wa.A;    N   = NP / A
        
        println(">> (Relativistic) Coulomb-sine normalized continuum orbital with normalization constant N=" * @sprintf("%.4e",N) *
                " and phase phi=" * @sprintf("%.4e",phi) * " at r=" * @sprintf("%.4e",pot.grid.r[mtp])   * " a.u." )
    end

    
    P = N .* cOrbital.P;   Q = N .* cOrbital.Q;   Pprime = N .* cOrbital.Pprime;   Qprime = N .* cOrbital.Qprime;   
    Continuum.reportContinuumAmplitude(P, cOrbital, pot.grid, "CoulombSine"; residual = wa.residual)
    newOrbital = Orbital( cOrbital.subshell, cOrbital.isBound, cOrbital.useStandardGrid, cOrbital.energy, 
                            P, Q, Pprime, Qprime, cOrbital.grid)

    normF = N * wa.A   # the fitted asymptotic amplitude of the P component after scaling, as a measure of the amplitude

    return( newOrbital, phi, normF )
end


"""
`Continuum.normalizeOrbitalOngRussek(cOrbital::Orbital, pot::Radial.Potential, settings::Continuum.Settings)`   
    ... to normalize the given continuum orbital with regard to a (asymptotic) pure sine-function. 
        An (on-energy-scale-normalized) orbital::Orbital is returned.
"""
function normalizeOrbitalOngRussek(cOrbital::Orbital, pot::Radial.Potential, settings::Continuum.Settings) 
    mtp = size( cOrbital.P, 1);                         energy = cOrbital.energy;      kappa = cOrbital.subshell.kappa
    wc   = Defaults.getDefaults("speed of light: c");   E = wc^2 + cOrbital.energy;    q = sqrt( 2*energy + energy/(wc*wc))     
    qq = sqrt(E^2/wc^2 - wc^2)
    
    dVdr = pot.Zr[mtp] / (pot.grid.r[mtp]^2);    V  = - pot.Zr[mtp] / pot.grid.r[mtp]
    Pe = cOrbital.P[mtp];   Pep = cOrbital.Pprime[mtp];       r0 = pot.grid.r[mtp]
    # Compute U-function
    U = wc /2. * ((E-V)*dVdr - wc^2 *kappa*(kappa+1) /r0^3) / (( (E-V)^2 -wc^4 - wc^2 *kappa*(kappa+1) / r0^2)^(3/2)) * Pe
    # U = U - wc * (Pep + dVdr * Pe / (2*(E-V+ 2* wc^2)) ) / (( (E-V)^2 -wc^4 - wc^2 *kappa*(kappa+1)/r0^2)^(1/2))
    U = U - wc * (Pep + dVdr * Pe / (2*(E-V+ 2* wc^2)) ) / (( (E-V)^2 -wc^4 - wc^2 *kappa*(kappa+1)/r0^2)^(1/2))
    # Compute A-function
    A = (( (E-V)^2 -wc^4 - wc^2 *kappa*(kappa+1)/r0^2)^(1/2)) / (E-V+wc^2) * (Pe^2 + U^2)
    A = sqrt(A);              phir0 = atan(U, Pe)
    phi = phir0 - q*r0;       N  = 1 / (A * sqrt(pi*wc))
    phi = rem(phi+1000pi, pi)    # to bring phi in the interval 0. <= phi < pi
    
    println(">> WKB (Ong-Russek) normalized continuum orbital with normalization constant N=" * @sprintf("%.4e",N) *
            " and phase phi=" * @sprintf("%.4e",phi) *
            " at r=" * @sprintf("%.4e",pot.grid.r[mtp])   * " a.u.  U/Pe = $(U/Pe)" )

    
    P = N .* cOrbital.P;   Q = N .* cOrbital.Q;   Pprime = N .* cOrbital.Pprime;   Qprime = N .* cOrbital.Qprime;   
    #  this routine's WKB envelope is phase-complete and stays;  the fit supplies the REPORTED amplitude and, through
    #  its residual, the check that the orbital has reached the asymptotic form at all
    wa = Continuum.asymptoticAmplitude(P, pot.grid, mtp, cOrbital.energy,
                                       Continuum.asymptoticCharge(Radial.determineZbar(pot)))
    Continuum.reportContinuumAmplitude(P, cOrbital, pot.grid, "OngRussek"; residual = wa.residual)
    newOrbital = Orbital( cOrbital.subshell, cOrbital.isBound, cOrbital.useStandardGrid, cOrbital.energy, 
                            P, Q, Pprime, Qprime, cOrbital.grid)

    normF = wa.A

    return( newOrbital, phi, normF )
end


"""
`Continuum.coulombMatchingPair(i::Int64, pot::Radial.Potential, energy::Float64, kappa::Int64, l::Int64, Zbar::Float64)`
    ... evaluates, AT THE SINGLE MESH POINT `i`, the pair of reference solutions onto which the Coulomb-function
        normalization matches the numerical orbital:  the regular and irregular relativistic Coulomb functions where the
        residual charge is appreciable, and the spherical Bessel pair where it is not.  Both large components and their
        radial derivatives are returned, together with `x = q r_i`, the classical turning point `xTP` in the same
        variable, and the Coulomb phase -- which depends on `eta`, `lambda`, `l` and `kappa` but NOT on the radius, and
        is therefore the same whichever point is matched at.  Nothing is printed, so the function may be called at many
        radii in turn.  A `(fu, gu, fup, gup, cPhase, kind, x, xTP)::NamedTuple` is returned.
"""
function coulombMatchingPair(i::Int64, pot::Radial.Potential, energy::Float64, kappa::Int64, l::Int64, Zbar::Float64)
    α = Defaults.getDefaults("alpha")  ;      wc = 1/α
    q  = sqrt( energy * (energy + 2 * wc^2) ) / wc  ;    x  = q * pot.grid.r[i]

    if abs(Zbar) > 0.1
        sa = "Normalization with Coulomb functions"
        λ  = sqrt(kappa^2 - Zbar^2 / wc^2)  ; λm1 = λ - 1.0
        η  = Zbar * α * (energy + wc^2) / sqrt( energy * (energy + 2 * wc^2) ) 

        xTP = η + sqrt(η^2 + λ*(λ+1.0))

        Δ = angle(SpecialFunctions.gamma(λm1 + 1 + im * η))
        if Δ >= 0.0 Δ = mod(Δ, 2pi) else Δ = -mod(-Δ, 2pi) end
    
        θ  = x - λm1*pi/2 - η*log(2x) + Δ                               #Eq 7.3
        if θ > 1e4   θ = mod(θ, 2pi) end
    
        hgfλ = twoFzero(im*η - λm1, im*η + λm1 + 1, im*2*x)
        hgfλm1 = twoFzero(im*η - λm1 + 1, im*η + λm1 + 2, im*2*x)
        hgfλm1 = im * hgfλm1 * (im*η - λm1) * (im*η + λm1 + 1) / ( 2.0 * x^2 )
    
        GiFλm1 = hgfλ * exp(im*θ)  ;   GPiFPλm1 = ( hgfλm1 + im *(1.0 - η/x) * hgfλ ) * exp(im*θ)
        gm_1 = GiFλm1.re ; fm_1 = GiFλm1.im   ;   gpm_1 = GPiFPλm1.re ; fpm_1 = GPiFPλm1.im

        if abs(gm_1*fpm_1 - fm_1*gpm_1 - 1.0) > 1e-15 
            ef= 0.0 ; eg = 0.0
            f, fp, g, gp = GSL.sf_coulomb_wave_FG_e(η, x, λm1, 0, ef, eg)
            gm_1 = g.val ; fm_1 = f.val ; gpm_1 = gp.val ; fpm_1 = fp.val
        end
    
        f = λ * ((λ/x + η/λ)*fm_1 - fpm_1) / sqrt(λ^2 + η^2)
        g = λ * ((λ/x + η/λ)*gm_1 - gpm_1) / sqrt(λ^2 + η^2)
    
        N  =  ( α^2 * Zbar^2 * (energy + 2 * wc^2)^2 + (kappa + λ)^2 * wc^2 * q^2 )^(-0.5) / λ

        # Coulomb phase shift for λ
        rnu = angle(Zbar * α * (energy + 2 * wc^2) - im * (kappa+λ) * sqrt(energy*(energy+2*wc^2)))
        Δ =  rnu - (λ-l-1.0)*pi/2 + Δ
        if ( Zbar < 0.0 && kappa < 0)  Δ = Δ - pi ; N = -N end
        if Δ >= 0.0 Δ = mod(Δ, 2pi) else Δ = -mod(-Δ, 2pi) end
    
        fu = N * ( (kappa + λ) * sqrt(λ^2 + η^2) * wc * q * f + α * Zbar * (λ * wc^2 - kappa * (energy + wc^2)) * fm_1 )
        gu = N * ( (kappa + λ) * sqrt(λ^2 + η^2) * wc * q * g + α * Zbar * (λ * wc^2 - kappa * (energy + wc^2)) * gm_1 )
    
        fl = -N * ( α * Zbar * sqrt(λ^2 + η^2) * wc * q * f + (kappa + λ) * (λ * wc^2 - kappa * (energy + wc^2)) * fm_1 )
        gl = -N * ( α * Zbar * sqrt(λ^2 + η^2) * wc * q * g + (kappa + λ) * (λ * wc^2 - kappa * (energy + wc^2)) * gm_1 )
    
        fup = -kappa * fu / pot.grid.r[i] + fl * (energy + pot.Zr[i]/pot.grid.r[i] + 2.0 * wc^2 ) / wc
        gup = -kappa * gu / pot.grid.r[i] + gl * (energy + pot.Zr[i]/pot.grid.r[i] + 2.0 * wc^2 ) / wc

    else
        
        sa = "Normalization with Bessel functions"

        if kappa < 0 ksign = 1 else ksign = -1 end

        fu =  x * SpecialFunctions.sphericalbesselj(l, x)
        gu = -x * SpecialFunctions.sphericalbessely(l, x)
        fl = -sqrt(energy / ( energy + 2*wc^2 )) * ksign * x * SpecialFunctions.sphericalbesselj(l + ksign, x)
        gl =  sqrt(energy / ( energy + 2*wc^2 )) * ksign * x * SpecialFunctions.sphericalbessely(l + ksign, x)

        fup = -kappa * fu / pot.grid.r[i] + fl * (energy + pot.Zr[i]/pot.grid.r[i] + 2.0 * wc^2 ) / wc
        gup = -kappa * gu / pot.grid.r[i] + gl * (energy + pot.Zr[i]/pot.grid.r[i] + 2.0 * wc^2 ) / wc

        Δ   = 0.0
        # A free wave has no classical turning point;  0.0 makes the `x < xTP` test below vacuous, which is correct.
        xTP = 0.0
    end

    return( (fu=fu, gu=gu, fup=fup, gup=gup, cPhase=Δ, kind=sa, x=x, xTP=xTP) )
end


"""
`Continuum.matchAtRadius(i::Int64, cOrbital::Orbital, pot::Radial.Potential, l::Int64, Zbar::Float64)`
    ... extracts, from the single mesh point `i`, the two numbers the Coulomb-function normalization is built on:  the
        inner (non-Coulomb) phase shift `iPhase` and the normalization constant `N`, obtained by matching the orbital's
        large component and its derivative to the reference pair at that point.  THE POINT OF HAVING THIS AS ITS OWN
        FUNCTION is that `iPhase` must be INDEPENDENT of `i` wherever the match is legitimate -- the numerical solution
        is `cos(iPhase) F + sin(iPhase) G` exactly, once the potential has reached its `-Zbar/r` tail -- so calling it
        at several radii turns an unverifiable one-point match into a measurement.  A
        `(N, iPhase, cPhase, kind, x, xTP)::NamedTuple` is returned.
"""
function matchAtRadius(i::Int64, cOrbital::Orbital, pot::Radial.Potential, l::Int64, Zbar::Float64)
    energy = cOrbital.energy  ;   kappa = cOrbital.subshell.kappa
    α = Defaults.getDefaults("alpha")  ;   wc = 1/α
    q = sqrt( energy * (energy + 2 * wc^2) ) / wc
    pair = Continuum.coulombMatchingPair(i, pot, energy, kappa, l, Zbar)

    p = cOrbital.P[i]    ;   pPrime = cOrbital.Pprime[i]

    δ = angle((p*pair.gup - pPrime*pair.gu) + im*(pPrime*pair.fu - p*pair.fup))
    if  ( abs(δ) > pi/2 )   δ = δ * ( 1.0 - pi/abs(δ))     end

    if  ( abs(p) > 1e-15 )  N = ( cos(δ) * pair.fu + sin(δ) * pair.gu ) / p
    else                    N = ( cos(δ) * pair.fup + sin(δ) * pair.gup ) / pPrime    end

    # Normalize wrt energy
    N = N * sqrt(q / pi / energy)

    return( (N=N, iPhase=δ, cPhase=pair.cPhase, kind=pair.kind, x=pair.x, xTP=pair.xTP) )
end


"""
`Continuum.checkMatchingRadius(cOrbital::Orbital, pot::Radial.Potential, l::Int64, Zbar::Float64, mtp::Int64;
                             warn::Float64=CONT_MATCH_WARN,   fatal::Float64=CONT_MATCH_FATAL,
                             warnN::Float64=CONT_SCALE_WARN,  fatalN::Float64=CONT_SCALE_FATAL,
                             meshSound::Float64=CONT_MESH_SOUND)`

    ... asks whether the ONE radius at which `normalizeOrbitalCoulombFunctions` matches the orbital was a fair place to
        do it, and moves it if not.  TWO quantities are extracted by `Continuum.matchAtRadius` over the outer two
        wavelengths and BOTH must be constant there, because each is blind to what the other catches:

        * the INNER PHASE SHIFT, a constant of the solution wherever the potential has reached its `-Zbar/r` tail.  Its
          spread catches an orbital whose `P` and `Pprime` are mutually inconsistent, i.e. one that does not satisfy the
          equation it is supposed to;
        * the NORMALIZATION CONSTANT itself, whose RELATIVE spread catches a distorted SHAPE.  This second measure is
          necessary rather than decorative:  the phase is the ANGLE of a complex number built from `P` and `Pprime`, so
          scaling both by the same factor leaves it EXACTLY unchanged.  Measured 05-Nov-2026 -- a 30 % multiplicative
          distortion of the outer half of an orbital, and 2 % point-to-point noise, each moved the phase spread by
          NOTHING and the normalization spread by 4.8 % and 4.0 % against 1.5e-04 for the sound orbital.

        A sound orbital has both near zero.  Either one large says the matching point carries something that is not the
        Coulomb solution, and the constant returned is then untrustworthy however confidently it was produced.

        WHAT THIS DOES NOT CATCH, stated because a guard claiming to catch everything is the suspicious kind:  a single
        enormous INTERIOR spike leaves both measures untouched, since neither looks inside the matching region.  That is
        what `Continuum.checkContinuumAmplitude` is for, and the two tiers are complementary, not alternatives.
        
        THE MEASURE IS SCALE-FREE AND NEEDS NO REFERENCE DATA.  `iPhase` is an angle, so its spread is in radians
        regardless of the orbital's arbitrary eigenvector scale, and the comparison is between radii of ONE orbital
        rather than against a tabulated value.
        
        IT ALSO SAYS WHICH GRID IS AT FAULT, which a single number cannot.  The ratio of the total variation of
        `iPhase` to its net change over the window separates the two failure modes:  near 1 the phase drifts
        SMOOTHLY, meaning the potential has not reached its tail inside the box and the BOX IS TOO SMALL;  far above 1
        it jitters from point to point, meaning the mesh does not resolve the oscillation and the MESH IS TOO COARSE.
        
        AND THE REPAIR MOVES INWARD, NOT OUTWARD, which is worth stating because it is the opposite of what one
        expects:  the matching point is `size(cOrbital.P, 1)`, the OUTERMOST mesh point of the orbital, so there is no
        larger radius to move to -- and that point is precisely where a B-spline Galerkin solution is most constrained
        by the boundary condition at the box wall.  So the repair is to step AWAY FROM THE WALL to the outermost radius
        at which the phase has settled.  Where no such radius exists on this grid, nothing in the code can repair it
        and the function RAISES, naming the grid change needed, rather than returning a number.
        
        AND THE MOVE IS GATED ON MESH DENSITY, which is a CORRECTION to the first form of this check rather than a
        refinement of it.  A move is taken only where the mesh carries at least `meshSound` points, counted, across one
        oscillation.  Measured 06-Nov-2026:  a move at 30 counted points reproduced the converged normalization to
        0.09 %, while the move the dielectronic-recombination test case took at about 18 turned a +6.0 % error into a
        -5.6 % one and reported success.  Below the gate the move is therefore DECLINED -- the ORIGINAL radius is kept,
        `declined` is set, and the verdict is `:poor`, so the caller is told loudly with the uncertainty quoted.  It
        does not raise:  the defect measured was a move that relocated an error while reporting success, and declining
        the move removes exactly that and nothing more.
        An orbital already sound at the box wall is untouched by this gate;  it bites only where a radius would
        otherwise have been moved.  Note what the gate does NOT claim:  the same sweep shows the error is non-monotone
        in the density, so a coarse mesh is a reason to distrust a MOVE rather than proof that the orbital is bad.

        A `(index, moved, declined, spread, jitter, scaleDrift, ppo, nProbes, inside, rTP, verdict)::NamedTuple` is
        returned, with `verdict` one of `:sound`, `:moved`, `:poor`, `:unfair` or `:tooFewPoints`.
            A `checkMatchingRadius::Continuum.checkMatchingRadius` is returned.
"""
function checkMatchingRadius(cOrbital::Orbital, pot::Radial.Potential, l::Int64, Zbar::Float64, mtp::Int64;
                             warn::Float64=CONT_MATCH_WARN,   fatal::Float64=CONT_MATCH_FATAL,
                             warnN::Float64=CONT_SCALE_WARN,  fatalN::Float64=CONT_SCALE_FATAL,
                             meshSound::Float64=CONT_MESH_SOUND)
    energy = cOrbital.energy
    α  = Defaults.getDefaults("alpha")  ;   wc = 1/α
    q  = sqrt( energy * (energy + 2 * wc^2) ) / wc
    λr = 2pi / q                                             # the de Broglie wavelength, in a.u.
    rOuter = pot.grid.r[mtp]

    # POINTS PER OSCILLATION IS COUNTED OVER THE OUTER WAVELENGTH, NOT OBTAINED BY DIVIDING BY A LOCAL STEP, and the
    # difference is a factor of two rather than a nicety.  JAC's radial mesh is a B-spline knot sequence sampled at
    # Gauss-Legendre points, so consecutive differences are NOT equal:  measured 06-Nov-2026 on a 25 a.u. box with
    # grid.hp = 0.0497, individual steps near the boundary run 0.70, 1.13 and 1.37 times grid.hp depending on which
    # sub-point of the pattern they span.  Dividing the wavelength by ONE such step therefore reports anything from
    # 0.7x to 1.4x the true sampling according to where the matching point happens to fall, which is noise rather
    # than a measurement.  Counting is what the name says and is immune to it.
    #   THE MEAN step over a whole wavelength is a steadier 0.87 to 1.00 times grid.hp across the systems measured,
    # so `wavelength/grid.hp` -- the quantity Continuum.gridConsistency tests and Radial.generateGrid constructs from
    # -- is a fair proxy for the counted density, running a few percent BELOW it.  It is the per-point quotient, not
    # that proxy, which is unusable.
    iOne = something( findfirst(k -> pot.grid.r[k] >= rOuter - λr, 1:mtp), 1 )
    ppo  = Float64(mtp - iOne + 1)                           # counted points per oscillation, carried to the caller

    # THE WINDOW IS THE OUTER TWO WAVELENGTHS, not a fixed number of mesh points.  A fixed point count is the error
    # that cost weeks on the sibling guard (see the amplitude note in Continuum.asymptoticAmplitude):  at a fine mesh
    # it spans a fraction of one oscillation and samples the wave at an arbitrary phase.
    rInner = max( pot.grid.r[2], rOuter - 2.0 * λr )
    iInner = something( findfirst(k -> pot.grid.r[k] >= rInner, 1:mtp), 1 )
    nAvail = mtp - iInner + 1
    if  nAvail < 16
        # Too few points in the outer two wavelengths to say anything;  that is itself a coarse-mesh statement,
        # and it is reported rather than silently passed.
        return( (index=mtp, moved=false, declined=false, spread=NaN, jitter=NaN, scaleDrift=NaN, ppo=ppo,
                 nProbes=nAvail, inside=false, rTP=0.0, verdict=:tooFewPoints) )
    end

    nProbes = min(24, nAvail)
    idx     = unique( round.(Int, range(iInner, mtp, length=nProbes)) )
    ms      = [ Continuum.matchAtRadius(i, cOrbital, pot, l, Zbar)  for i in idx ]
    δs      = [ m.iPhase  for m in ms ]
    Ns      = [ abs(m.N)  for m in ms ]

    # `iPhase` is folded into (-pi/2, pi/2] by matchAtRadius, so it lives on a circle of period pi and its spread
    # must be measured there.  A plain standard deviation reports a spurious pi/2 whenever the true phase sits near
    # the fold, which is a value the wave itself never has.
    spreadOf(v) = begin
        z = sum(exp.(2im .* v)) / length(v)
        μ = angle(z) / 2
        maximum( abs(rem(d - μ, pi, RoundNearest))  for d in v )
    end
    totalVariation(v) = sum( abs(rem(v[k+1] - v[k], pi, RoundNearest))  for k = 1:length(v)-1 )
    # The normalization constant's spread is taken RELATIVE to its own mean, so it is scale-free in the same sense
    # the phase is:  the orbital's arbitrary eigenvector scale divides out and what is left is the fractional
    # uncertainty the user would carry into a cross section.
    relSpread(v)      = (maximum(v) - minimum(v)) / max(sum(v)/length(v), 1.0e-300)

    worst  = spreadOf(δs)
    # `jitter` is the total variation measured AGAINST THE SPREAD, not against the net change.  Dividing by the net
    # change was the first form and it is useless:  for a sound orbital the net change is ~0 by construction, so the
    # ratio ran to 400 on orbitals that were perfectly fine.  Against the spread a monotone drift gives ~2 and
    # point-to-point noise gives much more, which is the distinction wanted.
    jitter = totalVariation(δs) / max(worst, 1.0e-12)
    sDrift = (maximum(Ns) - minimum(Ns)) / max(sum(Ns)/length(Ns), 1.0e-300)
    # WHETHER THE OUTERMOST MATCHING POINT LIES INSIDE THE CLASSICAL TURNING POINT is the physically unambiguous
    # half of the diagnosis, and it is what every badly-behaved case in the calibration set had in common:  inside
    # the turning point the wave does not oscillate at all, so no radius there can carry an asymptotic match.  It
    # is reported rather than inferred from the phase behaviour.
    inside = ms[end].x < ms[end].xTP
    rTP    = ms[end].xTP / q

    # THE REPAIR, IN THREE TIERS.  Walk inward from the box wall and take the OUTERMOST window of 8 consecutive probes
    # over which the phase has settled;  8 of 24 spans a quarter of the window, enough that one bad point cannot carry
    # it.  Below `warn` the match is fair and the orbital is used silently.  Between `warn` and `fatal` it is used but
    # said out loud, because the error it implies is percent-level rather than ruinous.  Above `fatal` there is no fair
    # radius on this grid at all and the caller must be stopped.
    # Below `meshSound` points per oscillation a moved radius is not a repair but a differently wrong answer; see
    # the measurement at CONT_MESH_SOUND.  The gate bites only on a MOVE, so an orbital that is sound at the wall is
    # unaffected however coarse the mesh.
    coarse = ppo < meshSound
    win  = 8
    # BOTH measures are judged per window, and `bad` expresses each in units of its OWN tolerance so that the worse
    # of the two decides.  That keeps the two thresholds independent -- they guard different failures and were
    # measured separately -- while leaving one number to compare windows by.
    best = (b=0, sp=Inf, sc=Inf, bad=Inf)
    for  b  in  length(idx):-1:win
        sp  = spreadOf(δs[b-win+1:b])
        sc  = relSpread(Ns[b-win+1:b])
        bad = max(sp/warn, sc/warnN)
        if  bad < best.bad    best = (b=b, sp=sp, sc=sc, bad=bad)    end
        if  bad < 1.0
            chosen = idx[b]
            # THE GATE DECLINES THE MOVE AND KEEPS THE ORIGINAL RADIUS;  it does not raise.  Raising was the first
            # form and it stopped six tests carrying approved references, which makes it a decision about what JAC
            # accepts rather than a correction.  Declining is the whole of what the measurement supports:  the
            # defect found was a move that RELOCATED an error and reported success, and not moving removes exactly
            # that.  The caller is told loudly, with the uncertainty quoted, through the `:poor` verdict.
            if  coarse  &&  chosen != mtp
                return( (index=mtp, moved=false, declined=true, spread=sp, jitter=jitter, scaleDrift=sc, ppo=ppo,
                         nProbes=length(idx), inside=inside, rTP=rTP, verdict=:poor) )
            end
            return( (index=chosen, moved=(chosen != mtp), declined=false, spread=sp, jitter=jitter, scaleDrift=sc,
                     ppo=ppo, nProbes=length(idx), inside=inside, rTP=rTP,
                     verdict=(chosen == mtp ? :sound : :moved)) )
        end
    end

    if  best.b > 0  &&  best.sp < fatal  &&  best.sc < fatalN
        chosen = coarse ? mtp : idx[best.b]
        return( (index=chosen, moved=(chosen != mtp), declined=(coarse && idx[best.b] != mtp), spread=best.sp,
                 jitter=jitter, scaleDrift=best.sc, ppo=ppo, nProbes=length(idx), inside=inside, rTP=rTP,
                 verdict=:poor) )
    end

    return( (index=mtp, moved=false, declined=false, spread=min(worst, best.sp), jitter=jitter,
             scaleDrift=(best.b > 0 ? best.sc : sDrift), ppo=ppo, nProbes=length(idx), inside=inside, rTP=rTP,
             verdict=:unfair) )
end


"""
`Continuum.matchingRadiusReport(wa::NamedTuple, cOrbital::Orbital, pot::Radial.Potential, mtp::Int64)`
    ... says out loud what `Continuum.checkMatchingRadius` found, and RAISES where the matching radius cannot be made
        fair on the given grid;  `nothing` is returned when the orbital may be used.  A sound orbital prints nothing,
        so the check is silent in the normal case and unmissable otherwise.
        
        WHY IT RAISES RATHER THAN WARNS IN THE LAST CASE.  A warning on a continuum orbital is read past:  the orbital
        is consumed inside an amplitude loop that may run thousands of times, and the resulting cross section carries
        no trace of it.  Where the phase never settles, the normalization constant is not a property of the solution
        at all, and every number built on it is wrong by an unknown factor -- which is the failure mode that reads as
        box dependence and has twice been mistaken for physics in this code.  The grid is NOT repaired here:  it
        belongs to the whole `Atomic.Computation`, and a silent re-solve would contract a continuum orbital from one
        grid against bound orbitals from another.
"""
function matchingRadiusReport(wa::NamedTuple, cOrbital::Orbital, pot::Radial.Potential, mtp::Int64)
    wa.verdict == :sound   &&   return( nothing )

    energy = cOrbital.energy
    α  = Defaults.getDefaults("alpha")  ;   wc = 1/α
    q  = sqrt( energy * (energy + 2 * wc^2) ) / wc  ;   λr = 2pi / q
    box = pot.grid.r[mtp] / λr                               # box, in wavelengths
    sh  = string(cOrbital.subshell)

    # WHICH GRID IS AT FAULT, and the order of these tests is the measured one rather than a guess.  EVERY badly
    # behaved case in the 05-Nov-2026 calibration set had its outermost matching point INSIDE the classical turning
    # point, which is physically decisive:  there the wave does not oscillate, so no radius in the box can carry the
    # match and only a larger box can help.  That is tested first and quotes the radius actually needed.  Only when
    # the turning point is comfortably inside the box is a coarse mesh the remaining explanation, and then the points
    # per oscillation says so directly -- ~40 is the measured sound value and 15.7 fails catastrophically.
    sPpo = @sprintf("%.1f", wa.ppo)   ;   sBox = @sprintf("%.1f", box)   ;   sJit = @sprintf("%.1f", wa.jitter)
    if      wa.inside
        cause  = "THE BOX IS TOO SMALL -- the outermost matching point lies INSIDE the classical turning point "  *
                 "(q*r = " * @sprintf("%.1f", pot.grid.r[mtp]*q) * " against a turning point at "                 *
                 @sprintf("%.1f", wa.rTP*q) * "), so the electron is still in the tunnelling region there and no " *
                 "radius in this box carries an oscillating wave to match onto"
        remedy = "enlarge grid.rbox to at least " * @sprintf("%.0f", 1.5*wa.rTP) * " a.u. (the turning point is "  *
                 "at " * @sprintf("%.0f", wa.rTP) * " a.u. and this box reaches "                                  *
                 @sprintf("%.0f", pot.grid.r[mtp]) * " a.u.)"
    elseif  wa.ppo < CONT_MESH_SOUND
        cause  = "THE MESH IS TOO COARSE -- this mesh carries only " * sPpo * " points per oscillation at this "   *
                 "energy, below the " * @sprintf("%.0f", CONT_MESH_SOUND) * " below which a MOVED matching "      *
                 "radius has been measured to relocate the error rather than remove it, and the extracted phase "  *
                 "varies " * sJit * "x its own spread from point to point"
        remedy = "reduce grid.hp by a factor of about " * @sprintf("%.1f", 40.0/max(wa.ppo, 1.0)) * ", from "     *
                 @sprintf("%.3e", pot.grid.hp) * " to about " * @sprintf("%.3e", pot.grid.hp*max(wa.ppo,1.0)/40.0) *
                 " a.u., which brings the counted density to about 40.  The factor is scaled from the COUNTED "   *
                 "density rather than assumed from grid.hp, since the mesh is sampled at Gauss-Legendre points and " *
                 "no single step represents it"
    else
        cause  = "NEITHER STANDARD CAUSE APPLIES -- the turning point is inside the box and the mesh carries "     *
                 sPpo * " points per oscillation, yet the phase still does not settle (varying " * sJit         *
                 "x its own spread)"
        remedy = "treat this orbital as suspect and report it:  the box holds " * sBox * " wavelengths, so the "   *
                 "usual two remedies do not obviously apply and the cause is worth finding"
    end

    if      wa.verdict == :tooFewPoints
        @warn("Continuum.checkMatchingRadius(): the outer two wavelengths of $sh at " *
              @sprintf("%.4e Ha hold only %d mesh points", energy, wa.nProbes) * ", too few to test whether the " *
              "normalization's matching radius is fair.  The orbital is used as it stands.  " * remedy * ".")
    elseif  wa.verdict == :moved
        printstyled("\n>> Continuum: THE MATCHING RADIUS WAS MOVED for $sh at " * @sprintf("%.4e", energy) *
                    " Ha.\n", color=:light_yellow)
        printstyled(">>   The outermost mesh point (index " * string(mtp) * ", r = " *
                    @sprintf("%.4f", pot.grid.r[mtp]) * " a.u.) did not give a settled phase;  matching moved "  *
                    "inward to index " * string(wa.index) * ", r = " * @sprintf("%.4f", pot.grid.r[wa.index])   *
                    " a.u., where the spread over 8 probes is " * @sprintf("%.2e", wa.spread) * " rad.\n",
                    color=:light_yellow)
        printstyled(">>   There the inner phase is constant to " * @sprintf("%.2e", wa.spread) *
                    " rad and the normalization constant to " * @sprintf("%.2f", 100*wa.scaleDrift) * " %.\n",
                    color=:light_yellow)
        printstyled(">>   This is normal close to a box wall, where a B-spline solution is constrained by its " *
                    "boundary condition.  The orbital IS usable;  the phase and normalization below are the moved " *
                    "ones.\n", color=:light_yellow)
    elseif  wa.verdict == :poor
        printstyled("\n>> Continuum: THE NORMALIZATION OF $sh AT " * @sprintf("%.4e", energy) *
                    " Ha IS ONLY APPROXIMATELY DEFINED.\n", color=:light_red)
        printstyled(">>   On the best window this grid offers, the inner phase spreads by " *
                    @sprintf("%.2e", wa.spread) * " rad (limit " * @sprintf("%.2e", CONT_MATCH_WARN) *
                    ") and the normalization constant by " * @sprintf("%.2f", 100*wa.scaleDrift) * " % (limit " *
                    @sprintf("%.2f", 100*CONT_SCALE_WARN) * " %).\n", color=:light_red)
        printstyled(">>   SO TREAT EVERY CROSS SECTION BUILT ON THIS ORBITAL AS CARRYING ABOUT " *
                    @sprintf("%.1f", 100*wa.scaleDrift) * " % UNCERTAINTY FROM THE NORMALIZATION ALONE.\n",
                    color=:light_red)
        if  wa.declined
            printstyled(">>   A RADIUS AT WHICH THE PHASE SETTLES DOES EXIST HERE AND WAS DELIBERATELY NOT USED, " *
                        "because this mesh carries only " * sPpo * " counted points per oscillation, below the " *
                        @sprintf("%.0f", CONT_MESH_SOUND) * " at which moving the radius has been measured to be " *
                        "a repair rather than a relocation of the error.\n", color=:light_red)
            printstyled(">>   The ORIGINAL matching radius is kept, so this number is the one the mesh gives;  it " *
                        "is not improved and not disguised.\n", color=:light_red)
        end
        printstyled(">>   DIAGNOSIS: " * cause * ".\n", color=:light_red)
        printstyled(">>   REMEDY: " * remedy * ".\n", color=:light_red)
        printstyled(">>   The orbital IS used.  Treat any cross section built on it as carrying this uncertainty.\n",
                    color=:light_red)
    elseif  wa.verdict == :unfair
        # WHY THERE ARE TWO WORDINGS HERE.  The gate at CONT_MESH_SOUND refuses a grid on which a settled radius DOES
        # exist, and saying "no radius works" in that case would be false and would send the reader looking for the
        # wrong thing.  The two endings are therefore kept apart:  one reports that nothing fair was found, the other
        # that something fair was found and deliberately not used.
        if  wa.ppo < CONT_MESH_SOUND
            sWhy = "The mesh carries only " * sPpo * " counted points per oscillation AND no radius in the outer "  *
                   "two wavelengths gives a settled phase, so the two failures coincide:  there is nothing fair to " *
                   "move to, and the mesh would not make a move trustworthy if there were.  The inner phase "       *
                   "spreads by " * @sprintf("%.2e", wa.spread) * " rad and the normalization constant by "          *
                   @sprintf("%.1f", 100*wa.scaleDrift) * " % over " * string(wa.nProbes) * " probes."
        else
            sWhy = "The inner phase shift, which is a CONSTANT of the solution wherever the match is legitimate, "  *
                   "spreads by " * @sprintf("%.2e", wa.spread) * " rad over the outer two wavelengths ("           *
                   string(wa.nProbes) * " probes) against a limit of " * @sprintf("%.2e", CONT_MATCH_WARN)         *
                   ", and the normalization constant itself varies by " * @sprintf("%.1f", 100*wa.scaleDrift)      *
                   " % across them against a limit of " * @sprintf("%.2f", 100*CONT_SCALE_FATAL) * " %.  So there " *
                   "is NO radius on this grid at which the orbital is the Coulomb solution it is being matched "    *
                   "to, and the constant returned would be an arbitrary number rather than a property of the "      *
                   "solution."
        end
        error("Continuum.normalizeOrbitalCoulombFunctions(): the normalization of $sh at " *
              @sprintf("%.4e", energy) * " Ha CANNOT BE TRUSTED on this grid.  " * sWhy * "  " *
              "DIAGNOSIS: " * cause * ".  REMEDY: " * remedy * ".  The grid is not changed here because it belongs " *
              "to the whole computation;  change it in the Radial.Grid you pass in.  To proceed anyway and accept " *
              "an unverified normalization, select another convention with " *
              "setDefaults(\"method: normalization, pure sine\"), which carries its own asymptotic-form guard.")
    end

    return( nothing )
end


"""
`Continuum.normalizeOrbitalCoulombFunctions(cOrbital::Orbital, pot::Radial.Potential, settings::Continuum.Settings)`   
    ... to normalize the given continuum orbital with regard to a (asymptotic) wave function as per Salvat Code. 
        The orbitasl are normalized to unit amplitude.
        An ( orbital::Orbital, (δ + Δ)::Float64, N::Float64 ) is returned.
        
        THE MATCHING RADIUS IS NOW TESTED BEFORE IT IS TRUSTED (05-Nov-2026).  This convention pins the orbital to the
        exact Coulomb functions at ONE radius, and for a long time nothing asked whether that radius was a fair place
        to do it:  where the potential had not reached its `-Zbar/r` tail there, or where the orbital was corrupted,
        the match still succeeded and returned a confident number.  `Continuum.checkMatchingRadius` now extracts the
        inner phase shift over the outer two wavelengths and requires it to be CONSTANT, which it is wherever the
        match is legitimate;  the matching point is moved away from the box wall when that helps, and the routine
        RAISES when no radius on the grid will do.  See that function for the measure and `Continuum.matchAtRadius`
        for the one-radius extraction this was factored out of.
"""
function normalizeOrbitalCoulombFunctions(cOrbital::Orbital, pot::Radial.Potential, settings::Continuum.Settings)
    mtp   = size( cOrbital.P, 1)     ;    energy = cOrbital.energy
    kappa = cOrbital.subshell.kappa  ;    l      = Basics.subshell_l(cOrbital.subshell) 
    Zbar  = -Radial.determineZbar(pot)

    # THE GUARD RUNS FIRST and may move the matching point or raise;  `matchAtRadius` then does at the chosen radius
    # exactly what this routine used to do at `mtp`, so a sound orbital follows the old path number for number.
    wa = Continuum.checkMatchingRadius(cOrbital, pot, l, Zbar, mtp)
    Continuum.matchingRadiusReport(wa, cOrbital, pot, mtp)
    m  = Continuum.matchAtRadius(wa.index, cOrbital, pot, l, Zbar)
    if  m.x < m.xTP    println("The kr is less than the Coulomb turning point")    end

    N = m.N  ;    δ = m.iPhase  ;    Δ = m.cPhase  ;    sa = m.kind

    P = N .* cOrbital.P  ;   Q = N .* cOrbital.Q  ;   Pprime = N .* cOrbital.Pprime  ;   Qprime = N .* cOrbital.Qprime  ;   
    #  NO ASYMPTOTIC-FIT TIER HERE, AND THE REASON IS THE CONVENTION ITSELF, not an oversight -- it was added on
    #  03-Nov-2026 for uniformity and WITHDRAWN the same day, which is worth recording.  The other three conventions
    #  normalize against the ASYMPTOTIC FORM of the wave, so a fit to that form is exactly their validity criterion.
    #  This one matches the orbital to the EXACT Coulomb (or Bessel) functions, which hold at any radius, so it is
    #  legitimately usable where the asymptotic expansion is not.  Measured: Ne-like Kr, Zbar ~ 27, 16.7 eV in a
    #  10 a.u. box gives eta ~ 24, where the logarithmic term eta*ln(2qr) ~ 76 dwarfs q*r ~ 11 and the asymptotic
    #  form is meaningless -- the fit residual is 4.77 and refused four suite tests, while this normalization is
    #  sound there and always has been.  A measure may only guard the assumption its own routine makes.
    #    AND THAT IS WHY THE GUARD ABOVE TAKES THE FORM IT DOES (05-Nov-2026).  The assumption THIS routine makes is
    #  not that the wave is asymptotic but that ONE radius suffices, so the measure asks exactly that: the inner
    #  phase shift is radius-independent wherever the match is legitimate, and its spread over the outer two
    #  wavelengths is therefore the convention's own validity criterion -- needing no asymptotic form, no reference
    #  data, and no absolute scale.
    Continuum.reportContinuumAmplitude(P, cOrbital, pot.grid, "CoulombFunctions")
    newOrbital = Orbital( cOrbital.subshell, cOrbital.isBound, cOrbital.useStandardGrid, cOrbital.energy, 
                            P, Q, Pprime, Qprime, cOrbital.grid)
    
    # The CONVENTION and its constant belong in this line, not only the phases.  The absolute scale of every
    # continuum-derived quantity depends on which of the four `Defaults.GBL_CONT_NORMALIZATION` conventions is in
    # force, and that is a session global a user never sees; a number printed without it cannot be read back.
    println(">> Coulomb-function-normalized continuum orbital (" * sa * ") with normalization constant N=" * @sprintf("%.4e",N) *
            ":   r = $(pot.grid.r[wa.index]),   iPhase = $δ,   cPhase = $Δ ")

    return( newOrbital, δ + Δ, N )
end


"""
`Continuum.twoFzero(CA::ComplexF64, CB::ComplexF64, CZ::ComplexF64)`
    ... Calculates the Hypergeometric function 2F0(CA,CB;1/CZ) hypergeometric asymptotic series.
        Taken from Radial package by Salvat et al.
        A ComplexF64 value is returned.  
"""
function twoFzero(CA::ComplexF64, CB::ComplexF64, CZ::ComplexF64)

    EPS=1.0E-16; ACCUR=0.5E-15; NTERM=75

    RRP=1.0
    RRN=0.0
    RIP=0.0
    RIN=0.0
    CDF=1.0 + 0.0im
    ERR2=0.0
    ERR3=1.0
    AR=0.0
    AF=0.0
    CF=0.0 + 0.0im

    for I = 1 : NTERM
        J=I-1
        CDF=CDF*(CA+J)*(CB+J)/(I*CZ)
        ERR1=ERR2
        ERR2=ERR3
        ERR3=abs(CDF)
        if (ERR1 > ERR2 && ERR2 < ERR3) break end
        AR=CDF.re
        if(AR > 0.0) 
        RRP=RRP+AR
        else
        RRN=RRN+AR
        end
        AI=(-im*CDF).re
        if(AI > 0.0) 
        RIP=RIP+AI
        else
        RIN=RIN+AI
        end
        CF=complex(RRP+RRN,RIP+RIN)
        AF=abs(CF)
        if(AF > 1.0e25) 
        CF=0.0 + 0im
        ERR=1.0
        break
        end
        if(ERR3 < 1.0e-25*AF || ERR3 < EPS) 
        ERR=EPS
        break
        end    
    end

    # ****  Round off error.

    TR=abs(RRP+RRN)
    if(TR > 1.0e-25) 
    ERRR=(RRP-RRN)*ACCUR/TR
    else
    ERRR=1.0e0
    end
    TI=abs(RIP+RIN)
    if(TI > 1.0e-25) 
    ERRI=(RIP-RIN)*ACCUR/TI
    else
    ERRI=1.0
    end

    #  ****  ... and truncation error.
    if(AF > 1.0e-25) 
    ERR=max(ERRR,ERRI)+ERR2/AF
    else
    ERR=max(ERRR,ERRI)
    end

    return CF

end


end # module
