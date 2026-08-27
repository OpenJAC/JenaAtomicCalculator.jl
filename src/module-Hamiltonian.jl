
"""
`module  JAC.Hamiltonian`  
	... a submodel of JAC that contains all structs and methods to efficiently set-up and work
	    with (large) Hamiltonian matrices under different conditions.
"""
module Hamiltonian

using  Printf, ..Basics, ..Bsplines, ..Defaults, ..InteractionStrength, ..InteractionStrengthQED, ..ManyElectron,
       ..Nuclear, ..Radial, ..RadialIntegrals, ..SpinAngularNew


"""
`Hamiltonian.CI_PARTIAL_DIAG_MIN_DIM`
    ... minimum symmetry-block CI matrix dimension above which Hamiltonian.diagonalizeCiMatrix considers a
        partial (LAPACK syevr-based) diagonalization instead of a full one, when settings.levelSelectionCI
        allows it (see diagonalizeCiMatrix). Below this size, the partial solver's own overhead can make it
        slightly SLOWER than a full diagonalization (confirmed 28-Jul-2026: dense random matrices showed a
        clean, monotonic crossover between not-worth-it and worth-it around n=10-20, with ~1.5x+ speedup
        comfortably established by n=20); chosen as a round, conservative value with margin above that
        crossover, not a hard-measured optimum.
"""
const  CI_PARTIAL_DIAG_MIN_DIM = 20


"""
`Hamiltonian.diagonalizeCiMatrix(matrix::Array{Float64,2}, levelSelectionCI::LevelSelection)`
    ... diagonalizes a single symmetry-block CI matrix, automatically choosing between a full and a partial
        (LAPACK syevr-based, via Basics.diagonalize's range keyword) diagonalization. Partial diagonalization
        is used only when levelSelectionCI selects a bounded set of global level INDICES (never symmetries --
        those already keep only whole, fully-wanted blocks via Basics.selectSymmetry, so there is no partial-
        eigenpair opportunity within a kept block) and the matrix is at least CI_PARTIAL_DIAG_MIN_DIM large.
        When used, only the lowest min(n, maximum(levelSelectionCI.indices)) eigenpairs are requested -- this
        is always sufficient and safe: the GLOBAL (energy-sorted, merged-across-blocks) rank of any level can
        never be smaller than its own LOCAL (within-block) rank, so if a level's global index is
        <= maximum(indices), its local index must be too; the true global top-`maximum(indices)` levels are
        therefore guaranteed to survive in the collected per-block partial spectra before the final
        merge-and-sort in performCI/performCIKinkAware. An  eigen::Basics.Eigen  is returned.
"""
function diagonalizeCiMatrix(matrix::Array{Float64,2}, levelSelectionCI::LevelSelection)
    n = size(matrix, 1)
    if  levelSelectionCI.active  &&  !isempty(levelSelectionCI.indices)  &&  n >= CI_PARTIAL_DIAG_MIN_DIM
        k = min(n, maximum(levelSelectionCI.indices))
        return( Basics.fixEigenvectorPhase!(Basics.diagonalize(MatrixWithLinearAlgebra(), matrix; range=1:k)) )
    else
        return( Basics.fixEigenvectorPhase!(Basics.diagonalize(MatrixWithLinearAlgebra(), matrix)) )
    end
end


"""
`Hamiltonian.projectHamiltonian(subshell::Subshell, matrix::Array{Float64,2}, matrixB::Array{Float64,2}, 
                                bVectors::Dict{Subshell, Vector{Float64}})` 
    ... projects the (single-electron DHF) Hamiltonian matrix of subshell with regard to the B-vectors of the same
        symmetry. This projection is necessary to ensure orthogonality among the various orbitals.
        A (nsL+nsS) x (nsL+nsS) matrix::Array{Float64,2} is returned.
"""
function projectHamiltonian(subshell::Subshell, matrix::Array{Float64,2}, matrixB::Array{Float64,2}, 
                            bVectors::Dict{Subshell, Vector{Float64}})
    nn = size(matrix, 1);   matrixP = deepcopy(matrix);           
    
    matrixI = zeros(nn, nn);   for i = 1:nn   matrixI[i,i] = 1.0   end
    
    for  (subsh, bVector)  in  bVectors
        if  subsh != subshell  &&  subshell.kappa == subsh.kappa
            ## println(">>> Project $subshell  upon $subsh  vectors")
            bbtMatrix = bVector * transpose(bVector)
            matrixP   = (matrixI - matrixB * bbtMatrix) * matrixP * (matrixI - bbtMatrix * matrixB)
        end
    end

    return( matrixP )
end


"""
`Hamiltonian.performCI(basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings; printout::Bool=false,
                       writeSummary::Bool=true)
    ... Computes and diagonalizes the Hamiltonian matrix for all CSF in the given basis. It also assigns the
        mixing coefficients to the individual levels. The individual contributions from the Breit or diagonal interaction as well as
        from QED to this matrix are controlled by the settings. A  multiplet::Multiplet  is returned.
        With writeSummary=false, the level table is not written to the summary file even when a summary file is
        open; this is needed by callers that invoke performCI inside a per-configuration loop (Cascade, Plasma)
        and would otherwise append one table per configuration to the .sum file.
"""
function performCI(basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings; printout::Bool=false,
                   writeSummary::Bool=true)
    
    # First determine the number of CSF in each J^P symmetry block
    symmetries = Dict{LevelSymmetry,Int64}()
    for  csf in basis.csfs
        sym = LevelSymmetry(csf.J, csf.parity)
        if     haskey(symmetries, sym)    symmetries[sym] = symmetries[sym] + 1
        else                              symmetries[sym] = 1
        end
    end

    # Test the total number of CSF
    NoCsf = 0;   for (sym,v) in symmetries   NoCsf = NoCsf + v   end
    if  NoCsf != length(basis.csfs)   error("stop b; NoCsf = $NoCsf ")   end

    # The radial-integral cache is created ONCE per performCI call and handed to every symmetry block, not
    # once per block as setupMatrix used to do (28-Jul-2026): its key is rank + subshell labels, which does
    # not depend on which block is being built, so an integral shared across blocks -- very common, since
    # many blocks reference the same subshell set -- is legitimately reusable across the whole computation.
    # OWNING IT HERE is what makes the labels-only key safe (15-Aug-2026): every block of this call shares
    # one BASIS, so a label identifies one orbital throughout, and the cache goes out of scope with the
    # call, so it cannot leak into the next performCI (a successive RAS layer, say) or into another task.
    # Until 15-Aug this was a module global wiped by an XL_*_reset_storage call here, and its correctness
    # rested on that wipe rather than on anything a reader of setupMatrix could see.
    xlCache = InteractionStrength.XLCache()

    # Calculate for each symmetry block the corresponding CI matrix, diagonalize it and append a Multiplet for this block
    multiplets = Multiplet[]
    for  (sym,v) in  symmetries
        # Skip the symmetry block if it not selected
        if  !Basics.selectSymmetry(sym, settings.levelSelectionCI)     continue    end
        matrix = Hamiltonian.setupMatrix(sym, basis, nm, grid, settings, xlCache; printout=printout)
        eigen  = Hamiltonian.diagonalizeCiMatrix(matrix, settings.levelSelectionCI)

        # Reassign state vectors to levels
        levels = Level[]
        for  ev = 1:length(eigen.values)
            # Construct the eigenvector with regard to the given basis (not w.r.t the symmetry block)
            evSym = eigen.vectors[ev];    vector = zeros( length(basis.csfs) );   ns = 0
            for  r = 1:length(basis.csfs) 
                if LevelSymmetry(basis.csfs[r].J, basis.csfs[r].parity) == sym    ns = ns + 1;   vector[r] = evSym[ns]   end
            end
            newlevel = Level( sym.J, AngularM64(sym.J.num//sym.J.den), sym.parity, 0, eigen.values[ev], 0., true, basis, vector ) 
            push!( levels, newlevel)
        end
        wa = Multiplet(string(sym) * "+", levels)
        push!( multiplets, wa)
    end
    
    # Merge all multiplets into a single one
    mp = Basics.merge(multiplets)
    mp = Basics.sortByEnergy(mp)
    
    # Determine the level list to be printed out
    levelNos = Int64[]
    for (ilev, level) in  enumerate(mp.levels)
        if  Basics.selectLevel(level, settings.levelSelectionCI)    push!(levelNos, ilev)    end
    end
    
    # Display all level energies and energy splittings
    if  printout
        Basics.tabulate(stdout, mp, levelNos)
    end
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary  &&  writeSummary
        Basics.tabulate(iostream, mp, levelNos)
    end

    return( mp )
end


"""
`Hamiltonian.performCIwithFrozenOrbitals(configs::Array{Configuration,1}, frozenOrbitals::Dict{Subshell, Orbital},
                                         nuclearModel::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings; printout::Bool=true)`
    ... to generate, from the given frozen orbitals, the CI multiplet for the given configurations: the CSF basis
        is set up for these configurations and the frozen orbitals, and is then diagonalized by Hamiltonian.performCI
        without any further orbital optimization. The e-e interaction that enters the CI matrix is whatever
        settings.eeInteractionCI specifies -- CoulombInteraction() (the AsfSettings() default), CoulombBreit(),
        ..., or DiagonalCoulomb() for the cheap, single-CSF-level approximation this function was once hard-wired
        to. A multiplet::Multiplet, sorted by energy, is returned.
"""
function performCIwithFrozenOrbitals(configs::Array{Configuration,1}, frozenOrbitals::Dict{Subshell, Orbital},
                                     nuclearModel::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings; printout::Bool=true)
    if  printout    println("\n... in Hamiltonian.perform...: a multiplet from frozen orbitals, CI with $(settings.eeInteractionCI) ...")    end

    # Generate a list of relativistic configurations and determine an ordered list of subshells for these configurations
    relconfList = ConfigurationR[]
    for  conf in configs
        wa = Basics.generateConfigurations(Basics.RelativisticConfigurations(), conf)
        append!( relconfList, wa)
    end
    if  printout    for  i = 1:length(relconfList)    println(">> include ", relconfList[i])    end   end
    subshellList = Basics.generateSubshellList(relconfList)
    Defaults.setDefaults("relativistic subshell list", subshellList; printout=printout)

    # Generate the relativistic CSF's for the given subshell list
    csfList = CsfR[]
    for  relconf in relconfList
        newCsfs = Basics.generateCsfRs(relconf, subshellList)
        append!( csfList, newCsfs)
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
    
    # Set-up a basis for calculating the Hamiltonian matrix
    basis = Basis(true, NoElectrons, subshellList, csfList, coreSubshellList, frozenOrbitals)

    # Hand the frozen-orbital basis to the standard CI machinery. NOTE (04-Aug-2026): this function used to
    # build only the DIAGONAL of the Hamiltonian matrix itself and therefore always returned single-CSF levels,
    # silently ignoring settings.eeInteractionCI (while erroring out on the one value, DiagonalCoulomb(), that
    # actually described its behaviour). That cost the genuine CSF mixing: for He-like carbon the two J=1 odd
    # levels came back as the bare jj CSFs (1s2p_1/2)_1 and (1s2p_3/2)_1 split by 1.330 eV, instead of the
    # near-pure LS states 1P1/3P1 split by 3.902 eV -- so the resonance line's strength was shared 1:2 between
    # them and the intercombination line was not reproduced at all (see examples/example-Je.jl, branch a).
    # Hamiltonian.performCI does the same job properly: it blocks the basis by J^P symmetry, assigns each
    # level's J and parity from its symmetry block rather than from a single dominant CSF, and honours
    # settings.eeInteractionCI throughout. The cheap original behaviour is NOT lost -- it is now reached
    # explicitly via settings.eeInteractionCI = DiagonalCoulomb(), which Hamiltonian.setupMatrix implements by
    # skipping every off-diagonal element BEFORE any spin-angular or radial-integral work is done, so its cost
    # is unchanged. Callers that need the cheap path (e.g. the Saha-Boltzmann level generation) simply say so.
    mp = Hamiltonian.performCI(basis, nuclearModel, grid, settings; printout=printout, writeSummary=printout)

    return( mp )
end



"""
`Hamiltonian.setupMatrix(sym::LevelSymmetry, basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings,
                          cache::InteractionStrength.XLCache; printout::Bool=false)
    ... Set-up (computes) the Hamiltonian matrix for all CSF with symmetry sym in the given basis. The individual contributions
        from the Breit or diagonal interaction as well as from QED to this matrix are controlled by the settings.
        A  matrix::Arrays{Float64,2}  is returned.
"""
function setupMatrix(sym::LevelSymmetry, basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings,
                      cache::InteractionStrength.XLCache; printout::Bool=false)

    # Determine the dimension of the CI matrix and the indices of the CSF with J^P symmetry in the basis
    idx_csf = Int64[]
    for  idx = 1:length(basis.csfs)
        if  basis.csfs[idx].J ==sym.J   &&   basis.csfs[idx].parity == sym.parity    push!(idx_csf, idx)    end
    end
    n = length(idx_csf)
    if printout    print("> Compute CI matrix of dimension $n x $n for the symmetry $(string(sym.J))^$(string(sym.parity)) ...")    end

    # Generate an effective nuclear charge Z(r) on the given grid to add QED contributions, if requested
    potential = Nuclear.nuclearPotential(nm, grid)
    if  settings.qedModel in [QedPetersburg(), QedSydney()]    
        meanPot = potential
        ## meanPot = compute("radial potential: Dirac-Fock-Slater", grid, basis)
        ## meanPot = Basics.add(potential, meanPot)   
    end   

    # Compute the Coulomb-(Breit-) interaction matrix. NOTE (28-Jul-2026): the radial-integral cache
    # reset/lifetime is owned by the CALLER (Hamiltonian.performCI), not here -- see that function's own
    # note. Resetting once per BLOCK here would discard cross-block-shared integrals (the same subshell
    # combination very often recurs across different symmetry blocks) for no benefit, since the cache key
    # never depends on which block is being built.
    #
    # Hermitian-symmetry shortcut (28-Jul-2026): only the UPPER triangle (r<=s) is ever computed below. This
    # is not merely safe but exact -- Basics.diagonalize(MatrixWithLinearAlgebra(),...) wraps this matrix in
    # LinearAlgebra.Symmetric(matrix), whose DEFAULT uplo=:U already reads ONLY the upper triangle and
    # completely ignores the lower one (confirmed: Symmetric([1 2 3; 999 5 6; 999 999 9]) == the honest
    # symmetric matrix, garbage lower-left included) -- so the lower triangle was always discarded even
    # before this change; no mirroring step is needed, only skipping the redundant r>s work.
    matrix = zeros(Float64, n, n)
    for  r = 1:n
        for  s = r:n
            if  settings.eeInteractionCI == DiagonalCoulomb()  &&  r != s    continue    end
            # Calculate the spin-angular coefficients
            subshellList = basis.subshells
            opa  = SpinAngularNew.OneParticleOperator(0, plus)
            waG1 = SpinAngularNew.computeCoefficients(opa, basis.csfs[idx_csf[r]], basis.csfs[idx_csf[s]], subshellList)
            opa  = SpinAngularNew.TwoParticleOperator(0, plus)
            waG2 = SpinAngularNew.computeCoefficients(opa, basis.csfs[idx_csf[r]], basis.csfs[idx_csf[s]], subshellList)
            wa   = [waG1, waG2]
            #
            me = 0.
            for  coeff in waG1
                me = me + coeff.T * RadialIntegrals.GrantIab(basis.orbitals[coeff.a], basis.orbitals[coeff.b], grid, potential)
                if  settings.qedModel != NoneQed()
                    me = me + InteractionStrengthQED.qedLocal(basis.orbitals[coeff.a], basis.orbitals[coeff.b], nm, settings.qedModel, meanPot, grid)
                end
            end

            for  coeff in waG2
                if  typeof(settings.eeInteractionCI) in [DiagonalCoulomb, CoulombInteraction, CoulombBreit, CoulombGaunt]
                    me = me + coeff.V * InteractionStrength.XL_Coulomb(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                                 basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid, cache)
                end
                                                                                        
                if      typeof(settings.eeInteractionCI) in [BreitInteraction, CoulombBreit, CoulombGaunt]
                    me = me + coeff.V * InteractionStrength.XL_Breit(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                               basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid,
                                                                               settings.eeInteractionCI, cache)     
                end
            end
            matrix[r,s] = me
        end
    end 
    if printout    println("   ... done.")    end

    return( matrix )
end


"""
`Hamiltonian.performCIKinkAware(basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings; printout::Bool=false)
    ... computes and diagonalizes the same CI Hamiltonian matrix as performCI, but via setupMatrixKinkAware (kink-aware
        two-electron Slater integral) instead of setupMatrix. Isolated from performCI; reached from
        SelfConsistent.performSCF when settings.scField = Basics.ALField(), and directly from
        SelfConsistent.solveOptimizedLevelField (EOLField), so that the FINAL reported level energies reflect
        the same kink-aware integral used during SCF orbital optimization, not just the SCF's own internal
        energy tracking. A  multiplet::Multiplet  is returned.
"""
function performCIKinkAware(basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings; printout::Bool=false)

    # First determine the number of CSF in each J^P symmetry block
    symmetries = Dict{LevelSymmetry,Int64}()
    for  csf in basis.csfs
        sym = LevelSymmetry(csf.J, csf.parity)
        if     haskey(symmetries, sym)    symmetries[sym] = symmetries[sym] + 1
        else                              symmetries[sym] = 1
        end
    end

    # Test the total number of CSF
    NoCsf = 0;   for (sym,v) in symmetries   NoCsf = NoCsf + v   end
    if  NoCsf != length(basis.csfs)   error("stop b; NoCsf = $NoCsf ")   end

    # One cache for the whole call, shared by every symmetry block; see performCI for why owning it here
    # is what makes a key of rank + subshell labels legitimate.  All three routes share it, the route being
    # part of the key.
    xlCache = InteractionStrength.XLCache()

    # Calculate for each symmetry block the corresponding CI matrix, diagonalize it and append a Multiplet for this block
    multiplets = Multiplet[]
    for  (sym,v) in  symmetries
        # Skip the symmetry block if it not selected
        if  !Basics.selectSymmetry(sym, settings.levelSelectionCI)     continue    end
        matrix = Hamiltonian.setupMatrixKinkAware(sym, basis, nm, grid, settings, xlCache; printout=printout)
        eigen  = Hamiltonian.diagonalizeCiMatrix(matrix, settings.levelSelectionCI)

        # Reassign state vectors to levels
        levels = Level[]
        for  ev = 1:length(eigen.values)
            # Construct the eigenvector with regard to the given basis (not w.r.t the symmetry block)
            evSym = eigen.vectors[ev];    vector = zeros( length(basis.csfs) );   ns = 0
            for  r = 1:length(basis.csfs)
                if LevelSymmetry(basis.csfs[r].J, basis.csfs[r].parity) == sym    ns = ns + 1;   vector[r] = evSym[ns]   end
            end
            newlevel = Level( sym.J, AngularM64(sym.J.num//sym.J.den), sym.parity, 0, eigen.values[ev], 0., true, basis, vector )
            push!( levels, newlevel)
        end
        wa = Multiplet(string(sym) * "+", levels)
        push!( multiplets, wa)
    end

    # Merge all multiplets into a single one
    mp = Basics.merge(multiplets)
    mp = Basics.sortByEnergy(mp)

    # Determine the level list to be printed out
    levelNos = Int64[]
    for (ilev, level) in  enumerate(mp.levels)
        if  Basics.selectLevel(level, settings.levelSelectionCI)    push!(levelNos, ilev)    end
    end

    # Display all level energies and energy splittings
    if  printout
        Basics.tabulate(stdout, mp, levelNos)
    end
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary
        Basics.tabulate(iostream, mp, levelNos)
    end

    return( mp )
end


#
# ===========================================================================================================
# MEASURED 10-Aug-2026 -- why there is NO iterative (Krylov/Davidson) eigensolver here, and why the caching
# below is already as good as it gets.  Recorded so that neither question has to be re-opened from scratch.
#
# Ne-like 1s^2 2s^2 2p^6, single + double excitations into n <= 5: 234 configurations, 19478 CSFs, largest
# symmetry block J = 2+ of dimension 2037 (2075703 upper-triangle pairs):
#
#     BUILD           158.1 s
#     DIAGONALIZE       1.1 s      dense LAPACK, ALL 2037 eigenvalues
#     ratio           138 : 1
#     memory           31.7 MB     sparsity 378903 of 4149369 nonzero (9.1%)
#
# 1) AN ITERATIVE EIGENSOLVER WOULD NOT HELP AT THESE SIZES, and would usually hurt.  Diagonalization is
#    0.7% of the cost.  A matrix-free Davidson needs H*v per iteration, and the matrix ELEMENTS are exactly
#    what is expensive: ~30 iterations touching the 9.1% nonzero elements costs roughly 2.7 builds, against
#    one build plus a 1.1 s LAPACK call -- unless every element is cached, which is building the matrix.
#    Build scales as n^2 (pairs), dense diagonalization as n^3; they cross only near n ~ 300000, where the
#    matrix would need ~700 GB.  So diagonalization time never becomes the binding constraint in practice.
#    Memory is not binding either at these sizes (31.7 MB here, ~3.2 GB at n = 20000) -- but the BUILD at
#    n = 20000 would take ~4.2 hours, so the build wall arrives long before the memory wall.
#    An iterative solver becomes a CAPABILITY enabler for n >~ 20000, where dense storage fails outright;
#    it is not a speed optimization for anything JAC reaches today.  KrylovKit.jl was considered and
#    deliberately NOT added as a dependency (user decision, 10-Aug-2026): BLAS is enough for now.
#    Basics.AbstractDiagonalizeTheme is the seam if this is ever revisited -- a new member of that family
#    keeps dense LAPACK available for the many JAC paths that want the WHOLE spectrum.
#
# 2) THE BUILD-LOCAL CACHING OF SLATER INTEGRALS IS ALREADY ~99% EFFECTIVE, so there is nothing to gain
#    there.  Measured over one symmetry block:
#
#        n <= 3, dim 118:   17632 coefficient terms ->   257 distinct integrals, 98.54% hits (68.6x reuse)
#        n <= 4, dim 684:  286452 coefficient terms ->  3259 distinct integrals, 98.86% hits (87.9x reuse)
#
#    The cost is therefore (number of DISTINCT integrals) x (cost of one integral), ~5.4 ms each, of which
#    the kink-aware screened potential is the bulk.  Making the CI build faster means making ONE Slater
#    integral cheaper -- not caching more of them, and not changing how the matrix is diagonalized.
# ===========================================================================================================
#


"""
`Hamiltonian.setupMatrixKinkAware(sym::LevelSymmetry, basis::Basis, nm::Nuclear.Model, grid::Radial.Grid,
                                   settings::AsfSettings, cache::InteractionStrength.XLCache; printout::Bool=false)
    ... sets up the same Hamiltonian matrix as setupMatrix, but using InteractionStrength.XL_CoulombKinkAware (kink-aware)
        instead of InteractionStrength.XL_Coulomb for the Coulomb two-electron contributions. The Breit contribution
        (XL_Breit) and the QED contribution are left untouched -- this bug is specific to the Coulomb Slater
        integral's quadrature, not to those other terms. Isolated from setupMatrix; only used by performCIKinkAware.
        A  matrix::Arrays{Float64,2}  is returned.
"""
function setupMatrixKinkAware(sym::LevelSymmetry, basis::Basis, nm::Nuclear.Model, grid::Radial.Grid,
                               settings::AsfSettings, cache::InteractionStrength.XLCache; printout::Bool=false)

    # Determine the dimension of the CI matrix and the indices of the CSF with J^P symmetry in the basis
    idx_csf = Int64[]
    for  idx = 1:length(basis.csfs)
        if  basis.csfs[idx].J ==sym.J   &&   basis.csfs[idx].parity == sym.parity    push!(idx_csf, idx)    end
    end
    n = length(idx_csf)
    if printout    print("> Compute CI matrix of dimension $n x $n for the symmetry $(string(sym.J))^$(string(sym.parity)) ...")    end

    # Generate an effective nuclear charge Z(r) on the given grid to add QED contributions, if requested
    potential = Nuclear.nuclearPotential(nm, grid)
    if  settings.qedModel in [QedPetersburg(), QedSydney()]
        meanPot = potential
    end

    # Compute the Coulomb-(Breit-) interaction matrix. NOTE (28-Jul-2026): the radial-integral cache
    # reset/lifetime is owned by the CALLER (Hamiltonian.performCIKinkAware), not here -- identical reasoning to
    # setupMatrix's own note.
    #
    # Hermitian-symmetry shortcut (28-Jul-2026): only the UPPER triangle (r<=s) is ever computed below --
    # see setupMatrix's identical note for why this is exact, not just safe (Symmetric(matrix)'s default
    # uplo=:U already discards the lower triangle in Basics.diagonalize).
    matrix = zeros(Float64, n, n)
    for  r = 1:n
        for  s = r:n
            if  settings.eeInteractionCI == DiagonalCoulomb()  &&  r != s    continue    end
            # Calculate the spin-angular coefficients
            subshellList = basis.subshells
            opa  = SpinAngularNew.OneParticleOperator(0, plus)
            waG1 = SpinAngularNew.computeCoefficients(opa, basis.csfs[idx_csf[r]], basis.csfs[idx_csf[s]], subshellList)
            opa  = SpinAngularNew.TwoParticleOperator(0, plus)
            waG2 = SpinAngularNew.computeCoefficients(opa, basis.csfs[idx_csf[r]], basis.csfs[idx_csf[s]], subshellList)
            wa   = [waG1, waG2]
            #
            me = 0.
            for  coeff in waG1
                me = me + coeff.T * RadialIntegrals.GrantIab(basis.orbitals[coeff.a], basis.orbitals[coeff.b], grid, potential)
                if  settings.qedModel != NoneQed()
                    me = me + InteractionStrengthQED.qedLocal(basis.orbitals[coeff.a], basis.orbitals[coeff.b], nm, settings.qedModel, meanPot, grid)
                end
            end

            for  coeff in waG2
                if  typeof(settings.eeInteractionCI) in [DiagonalCoulomb, CoulombInteraction, CoulombBreit, CoulombGaunt]
                    me = me + coeff.V * InteractionStrength.XL_CoulombKinkAware(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                                        basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid, cache)
                end

                if      typeof(settings.eeInteractionCI) in [BreitInteraction, CoulombBreit, CoulombGaunt]
                    me = me + coeff.V * InteractionStrength.XL_Breit(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                               basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid,
                                                                               settings.eeInteractionCI, cache)
                end
            end
            matrix[r,s] = me
        end
    end
    if printout    println("   ... done.")    end

    return( matrix )
end

end # module
