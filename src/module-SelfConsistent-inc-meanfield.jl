

# The mean-field path, DFSField and HSField: the plain iteration and its Anderson-accelerated form.

"""
`SelfConsistent.solveMeanFieldBasis(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives, 
                                    settings::AsfSettings; printout::Bool=true)` 
    ... solves the self-consistent field for the given orbitals (from basis), the nuclear model as well as
        the (local) mean-field potential as specified by the settings::AsfSettings. A (new) basis::Basis is returned.
"""
function solveMeanFieldBasis(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives, 
                             settings::AsfSettings; printout::Bool=true)
    # Defaults.setDefaults("standard grid", primitives.grid; printout=printout)
    # Define the storage for the calculations of matrices
    if  printout    println(">> (Re-) Define a storage array for dealing with single-electron TTp B-spline matrices:")    end
    storage  = Dict{String,Array{Float64,2}}()
    
    # Set-up the overlap matrix; compute or fetch the diagonal 'overlap' blocks
    nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS
    wb  = zeros( nsL+nsS, nsL+nsS )
    wb[1:nsL,1:nsL]                 = Bsplines.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
    wb[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = Bsplines.generateTTpMatrix!("SS-overlap", 0, primitives, storage)
    
    # Determine te nuclear potential once at the beginning
    nuclearPotential  = Nuclear.nuclearPotential(nuclearModel, primitives.grid)
    
    # Determine the symmetry block of this basis and define storage for the kappa blocks and orbitals from the last iteration
    kappas   = Int64[];   for sh in basis.subshells  push!(kappas, sh.kappa)   end;   kappas = unique(kappas)

    bsplineBlock = Dict{Int64,Basics.Eigen}();   previousOrbitals = deepcopy(basis.orbitals)
    for  kappa  in  kappas           bsplineBlock[kappa]  = Basics.Eigen( zeros(2), [zeros(2), zeros(2)])   end
    
    # Start the SCF procedure for all symmetries
    isNotSCF = true;   NoIteration = 0;   accuracyScf = 0.
    while  isNotSCF
        NoIteration = NoIteration + 1;   go_on = false 
        if  NoIteration >  settings.maxIterationsScf
                println(">> Maximum number of SCF iterations = $(settings.maxIterationsScf) is reached at accuracy " * 
                        @sprintf("%.4e", accuracyScf) * " ... computations proceed.")
                # Collected as well; see the note at solveAverageAtomField above.
                Defaults.warn(AddWarning(), "SelfConsistent.solveMeanFieldBasis(): the SCF did NOT converge for " *
                              string(basis.subshells) * " -- stopped at accuracy " * @sprintf("%.1e", accuracyScf) *
                              " after $(settings.maxIterationsScf) iterations.")
            break
        end
        if  printout    println("\nIteration $NoIteration for symmetries ... ")    end
        for kappa in kappas
            # (1) First re-define an (arbitrary) 'level' that represents the mean occupation for the local potential
            wBasis = Basis(true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, previousOrbitals)
            NoCsf  = length(wBasis.csfs)
            wmc    = zeros( NoCsf );   wN = 0.
            for i = 1:NoCsf   wmc[i] = Basics.twice(wBasis.csfs[i].J) + 1.0;   wN = wN + abs(wmc[i])^2    end
            for i = 1:NoCsf   wmc[i] = wmc[i] / sqrt(wN)   end
            wLevel = Level( AngularJ64(0), AngularM64(0), Basics.plus, 0, -1., 0., true, wBasis, wmc)
            # (2) Re-compute the local potential
            wp  = Basics.computePotential(settings.scField, primitives.grid, wLevel)
            pot = Basics.add(nuclearPotential, wp)
            # (3) Set-up the diagonal part of the Hamiltonian matrix
            wa = Bsplines.setupLocalMatrix(kappa, primitives, pot, storage)
            # (4) Solve the generalized eigenvalue problem
            wc = Bsplines.diagonalizeLocalMatrix(kappa, wa, wb, primitives)
            # (5) Analyse and print information about the convergence of the symmetry blocks and the occupied orbitals
            wcBlock = Basics.analyzeConvergence(bsplineBlock[kappa], wc)
            if  wcBlock > 1.000 * settings.accuracyScf   go_on = true   end
            for  sh in basis.subshells
                if      sh in settings.frozenSubshells   ## do nothing
                elseif  sh.kappa == kappa
                    newOrbital = Bsplines.generateOrbitalFromPrimitives(sh, wc, primitives)
                    wcOrbital  = Basics.analyzeConvergence(previousOrbitals[sh], newOrbital)
                    if  wcOrbital > settings.accuracyScf   accuracyScf = wcOrbital;   go_on = true   end
                        sa = "  $sh::  en [a.u.] = " * @sprintf("%.7e", newOrbital.energy) * ";   self-cons'cy = "  
                        sa = sa * @sprintf("%.4e", wcOrbital)   * "  ["
                        sa = sa * @sprintf("%.4e", wcBlock)             * " for sym-block kappa = $kappa]"
                        if  printout    println(sa)    end
                    previousOrbitals[sh] = newOrbital
                end
            end
            # (6) Re-define the bsplineBlock
            bsplineBlock[kappa] = wc
        end
        if  go_on   nothing   else   break   end
    end
    
    analyzedOrbitals = Basics.analyze(previousOrbitals, printout=true)

    newBasis = Basis(true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, analyzedOrbitals)
    return( newBasis )
end


"""
`SelfConsistent.solveMeanFieldBasisAnderson(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives,
                                          settings::AsfSettings; printout::Bool=true, andersonDepth::Int64=2)`
    ... solves the same mean-field (DFS/HS) self-consistent problem as SelfConsistent.solveMeanFieldBasis, but
        accelerates the iteration with Anderson mixing. The physics is identical -- the same functional, the same
        fixed point; only the route to it differs, and the returned basis::Basis is the SAME self-consistent
        solution, reached in fewer iterations.

        WHY. The plain iteration is an undamped fixed-point (Picard) map: build the screening potential from the
        current orbitals, solve the one-electron Dirac equation in it, take the resulting orbitals as the new
        ones. Such a map converges LINEARLY, the residual shrinking by a roughly constant factor r each step, so
        the cost to reach a given accuracy scales as log(tol)/log(r). Measured on this code (12-Aug-2026), r is
        0.44 for Ar 1s^2..3p^6, 0.57 for Ne 1s^2 2s^2 2p^6 and 0.69 for Fe [Ar] 3d^6 4s^2 -- 17, 25 and 37
        iterations respectively to reach 1e-6, and the open-3d case additionally sits at residual ~1 for seven
        iterations before descending at all. Anderson mixing forms the next screening potential from a
        least-squares combination of the last few iterates and their residuals, which cancels the slowest-decaying
        error components instead of waiting for them to decay geometrically.

        WHAT IS MIXED. The fixed-point variable is the ELECTRONIC screening potential on the radial grid, not the
        orbitals: it is one vector, the nuclear part is constant and drops out, and the map
        potential -> orbitals -> potential is exactly the self-consistency condition. andersonDepth = 0 recovers
        the plain iteration (and is then Jacobi rather than Gauss-Seidel, see below); the default 2 is the measured optimum here (3 is nearly equal, larger is WORSE: 24 iterations at 5 and 36 at 12 for Ne), and it
        choice and rarely needs changing.

        ONE DELIBERATE DIFFERENCE FROM solveMeanFieldBasis. There the potential is rebuilt inside the kappa loop,
        so each symmetry is solved in a field the earlier symmetries have already improved (a Gauss-Seidel
        sweep); here it is built once per iteration and all symmetries see the same field (a Jacobi sweep), which
        is what makes a single well-defined residual per iteration -- and hence Anderson -- possible. Both have
        the SAME fixed point, so a converged result is unaffected; only an iteration stopped short of
        convergence could differ.

        A basis::Basis with the self-consistent orbitals is returned.
"""
function solveMeanFieldBasisAnderson(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives,
                                   settings::AsfSettings; printout::Bool=true, andersonDepth::Int64=2)
    storage  = Dict{String,Array{Float64,2}}()
    nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS
    wb  = zeros( nsL+nsS, nsL+nsS )
    wb[1:nsL,1:nsL]                 = Bsplines.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
    wb[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = Bsplines.generateTTpMatrix!("SS-overlap", 0, primitives, storage)
    nuclearPotential  = Nuclear.nuclearPotential(nuclearModel, primitives.grid)
    kappas   = Int64[];   for sh in basis.subshells  push!(kappas, sh.kappa)   end;   kappas = unique(kappas)
    bsplineBlock = Dict{Int64,Basics.Eigen}();   previousOrbitals = deepcopy(basis.orbitals)
    for  kappa  in  kappas           bsplineBlock[kappa]  = Basics.Eigen( zeros(2), [zeros(2), zeros(2)])   end

    # The screening potential built from a given set of orbitals -- the right-hand side of the self-consistency
    # condition, and the vector Anderson works on.
    function screeningZr(orbitals)
        wBasis = Basis(true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, orbitals)
        NoCsf  = length(wBasis.csfs)
        wmc    = zeros( NoCsf );   wN = 0.
        for i = 1:NoCsf   wmc[i] = Basics.twice(wBasis.csfs[i].J) + 1.0;   wN = wN + abs(wmc[i])^2    end
        for i = 1:NoCsf   wmc[i] = wmc[i] / sqrt(wN)   end
        wLevel = Level( AngularJ64(0), AngularM64(0), Basics.plus, 0, -1., 0., true, wBasis, wmc)
        return( Basics.computePotential(settings.scField, primitives.grid, wLevel).Zr )
    end

    xUsed    = screeningZr(previousOrbitals)          ## the screening potential the next solve will use
    xHistory = Vector{Float64}[];   fHistory = Vector{Float64}[]
    isNotSCF = true;   NoIteration = 0;   accuracyScf = 0.
    while  isNotSCF
        NoIteration = NoIteration + 1;   go_on = false
        if  NoIteration >  settings.maxIterationsScf
                println(">> Maximum number of SCF iterations = $(settings.maxIterationsScf) is reached at accuracy " *
                        @sprintf("%.4e", accuracyScf) * " ... computations proceed.")
                # Collected as well; see the note at solveAverageAtomField above.
                Defaults.warn(AddWarning(), "SelfConsistent.solveMeanFieldBasisAnderson(): the SCF did NOT converge for " *
                              string(basis.subshells) * " -- stopped at accuracy " * @sprintf("%.1e", accuracyScf) *
                              " after $(settings.maxIterationsScf) iterations.")
            break
        end
        if  printout    println("\nIteration $NoIteration for symmetries ... ")    end
        pot = Basics.add(nuclearPotential, Radial.Potential("mean field", xUsed, primitives.grid))
        for kappa in kappas
            wa = Bsplines.setupLocalMatrix(kappa, primitives, pot, storage)
            wc = Bsplines.diagonalizeLocalMatrix(kappa, wa, wb, primitives)
            wcBlock = Basics.analyzeConvergence(bsplineBlock[kappa], wc)
            if  wcBlock > 1.000 * settings.accuracyScf   go_on = true   end
            for  sh in basis.subshells
                if      sh in settings.frozenSubshells   ## do nothing
                elseif  sh.kappa == kappa
                    newOrbital = Bsplines.generateOrbitalFromPrimitives(sh, wc, primitives)
                    wcOrbital  = Basics.analyzeConvergence(previousOrbitals[sh], newOrbital)
                    if  wcOrbital > settings.accuracyScf   accuracyScf = wcOrbital;   go_on = true   end
                        sa = "  $sh::  en [a.u.] = " * @sprintf("%.7e", newOrbital.energy) * ";   self-cons'cy = "
                        sa = sa * @sprintf("%.4e", wcOrbital)   * "  ["
                        sa = sa * @sprintf("%.4e", wcBlock)             * " for sym-block kappa = $kappa]"
                        if  printout    println(sa)    end
                    previousOrbitals[sh] = newOrbital
                end
            end
            bsplineBlock[kappa] = wc
        end
        if  go_on   nothing   else   break   end
        # The self-consistency residual: what the orbitals just obtained say the screening potential should be,
        # minus what was actually used to obtain them. It vanishes exactly at the self-consistent solution.
        gNew = screeningZr(previousOrbitals);    fNew = gNew - xUsed
        if  andersonDepth <= 0
            xUsed = gNew                                                        ## plain (Jacobi) iteration
        else
            push!(xHistory, copy(xUsed));   push!(fHistory, fNew)
            if  length(xHistory) > andersonDepth + 1
                popfirst!(xHistory);   popfirst!(fHistory)
            end
            m = length(xHistory) - 1
            if  m < 1     xUsed = gNew                                          ## no history yet
            else
                # Least squares: choose the combination of the last m residual DIFFERENCES that best cancels
                # the current residual, then apply it to the iterates as well (standard Anderson, beta = 1).
                dF = zeros( length(fNew), m );   dX = zeros( length(fNew), m )
                for j = 1:m   dF[:,j] = fHistory[j+1] - fHistory[j];   dX[:,j] = xHistory[j+1] - xHistory[j]   end
                local gamma
                try
                    gamma = dF \ fNew
                catch
                    gamma = zeros(m)      ## a rank-deficient window: fall back to the plain step
                end
                if  any(!isfinite, gamma)   gamma = zeros(m)   end
                xUsed = xUsed + fNew - (dX + dF) * gamma
            end
        end
    end

    analyzedOrbitals = Basics.analyze(previousOrbitals, printout=true)
    newBasis = Basis(true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, analyzedOrbitals)
    return( newBasis )
end
