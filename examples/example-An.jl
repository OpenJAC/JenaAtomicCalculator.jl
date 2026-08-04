#
println("An) Apply & test the bi-orthogonal transformation of two independently-generated multiplets.")

using LinearAlgebra

if  true
    # Last successful:  30-Jul-2026
    # Branch a: Stage 1 -- orbital biorthogonalization. Generates an initial- and a final-state multiplet
    #   via two INDEPENDENT SCF processes (exactly as example-Da.jl does today), then confirms two things:
    #   (1) the two multiplets' orbitals, sharing subshell labels like 3p_3/2, are indeed slightly
    #   non-orthogonal to each other today (this is the real, present gap BiOrthogonal.jl addresses, not a
    #   hypothetical one); (2) BiOrthogonal.generateBiorthogonalShellMatrices makes them exactly
    #   biorthonormal (to floating-point precision), and correctly reduces to the identity transformation
    #   when both sides already share the same orbitals.
    grid = Radial.Grid(true)
    nm   = Nuclear.Model(26.)
    photoSettings = PhotoEmission.Settings(PhotoEmission.Settings(), multipoles=[E1], gauges=[UseCoulomb,UseBabushkin], printBefore=false)
    comp = Atomic.Computation(Atomic.Computation(), name="An: bi-orthogonal Stage 1", grid=grid, nuclearModel=nm,
                              initialConfigs = [Configuration("[Ne] 3s 3p^6")],
                              finalConfigs   = [Configuration("[Ne] 3s^2 3p^5")],
                              processSettings = photoSettings )
    wb = perform(comp; output=true)
    initialMultiplet = wb["initialMultiplet"];   finalMultiplet = wb["finalMultiplet"]
    iBasis = initialMultiplet.levels[1].basis;   fBasis = finalMultiplet.levels[1].basis

    println("\n>> Orbital overlaps BEFORE the bi-orthogonal transformation (should deviate slightly from 1):")
    for sh in intersect(iBasis.subshells, fBasis.subshells)
        ov = RadialIntegrals.overlap(iBasis.orbitals[sh], fBasis.orbitals[sh], grid)
        println("     <$sh (initial)|$sh (final)> = $ov")
    end

    newLeft, newRight = BiOrthogonal.generateBiorthogonalShellMatrices(iBasis, fBasis, grid)
    println("\n>> Orbital overlaps AFTER the bi-orthogonal transformation (should be exactly 1):")
    newOverlaps = [ RadialIntegrals.overlap(newLeft[sh], newRight[sh], grid)  for sh in iBasis.subshells ]
    for  (sh, ov)  in  zip(iBasis.subshells, newOverlaps)
        println("     <$sh (new-left)|$sh (new-right)> = $ov")
    end
    println(">> worst |overlap - 1| = $(maximum(abs.(newOverlaps .- 1.0)))")

    newLeft2, newRight2 = BiOrthogonal.generateBiorthogonalShellMatrices(iBasis, iBasis, grid)
    worst2 = maximum( abs.(newLeft2[sh].P[1:min(length(newLeft2[sh].P),length(iBasis.orbitals[sh].P))] -
                            iBasis.orbitals[sh].P[1:min(length(newLeft2[sh].P),length(iBasis.orbitals[sh].P))])
                        for sh in iBasis.subshells ) |> maximum
    println(">> trivial-identity check (same basis on both sides): worst |P_new - P_old| = $worst2")
    #
elseif  true
    # Last successful:  30-Jul-2026
    # Branch b: Stage 3 -- CI-vector counter-rotation, validated in isolation from Stage 1's specific
    #   biorthogonalizing choice. A minimal complete active space (1 electron among {3s,4s}, so the full CI
    #   here is exactly invariant under any orthogonal rotation mixing 3s and 4s) is rotated by a known
    #   angle theta; BiOrthogonal.generateCounterRotatingCiMatrices is fed this rotation DIRECTLY (bypassing
    #   Stage 1's LU-decomposition-derived transformation matrix) and used to predict the mixing
    #   coefficients in the rotated basis. This prediction is compared against an INDEPENDENT
    #   re-diagonalization of the same Hamiltonian using the (genuinely self-orthonormal, since the
    #   rotation is orthogonal) rotated orbitals directly -- a real, external ground truth, not something
    #   derived from BiOrthogonal.jl's own machinery.
    grid = Radial.Grid(true)
    nm   = Nuclear.Model(26.)
    comp = Atomic.Computation(Atomic.Computation(), name="An: bi-orthogonal Stage 3", grid=grid, nuclearModel=nm,
                              configs = [Configuration("[Ne] 3s"), Configuration("[Ne] 4s")],
                              asfSettings = AsfSettings() )
    wb = perform(comp; output=true)
    multiplet = wb["multiplet:"];   basis = multiplet.levels[1].basis
    nm2 = Nuclear.Model(26.);       asfSettings = AsfSettings()

    targetLevel = multiplet.levels[1]
    JP          = LevelSymmetry(targetLevel.J, targetLevel.parity)
    idxCsf      = [i for (i,c) in enumerate(basis.csfs) if c.J == JP.J && c.parity == JP.parity]

    theta = 0.29
    Q11, Q12 =  cos(theta), sin(theta)
    Q21, Q22 = -sin(theta), cos(theta)
    sh3s = Subshell("3s_1/2");   sh4s = Subshell("4s_1/2")
    o3   = basis.orbitals[sh3s]; o4 = basis.orbitals[sh4s]
    len  = max(length(o3.P), length(o4.P))
    padTo(v,n) = length(v) >= n ? v[1:n] : vcat(v, zeros(n-length(v)))
    P3,P4   = padTo(o3.P,len),      padTo(o4.P,len);        Q3,Q4   = padTo(o3.Q,len),      padTo(o4.Q,len)
    Pp3,Pp4 = padTo(o3.Pprime,len), padTo(o4.Pprime,len);   Qp3,Qp4 = padTo(o3.Qprime,len), padTo(o4.Qprime,len)
    new3s = Radial.Orbital(o3.subshell, o3.isBound, o3.useStandardGrid, o3.energy, Q11*P3+Q21*P4, Q11*Q3+Q21*Q4, Q11*Pp3+Q21*Pp4, Q11*Qp3+Q21*Qp4, grid)
    new4s = Radial.Orbital(o4.subshell, o4.isBound, o4.useStandardGrid, o4.energy, Q12*P3+Q22*P4, Q12*Q3+Q22*Q4, Q12*Pp3+Q22*Pp4, Q12*Qp3+Q22*Qp4, grid)

    rotatedOrbitals = copy(basis.orbitals)
    rotatedOrbitals[sh3s] = new3s;   rotatedOrbitals[sh4s] = new4s
    rotatedBasis = Basis(basis.isDefined, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, rotatedOrbitals)

    println("\n>> INDEPENDENT reference: re-diagonalize H directly using the (self-orthonormal) rotated orbitals")
    Hmat_rot        = Basics.compute(Basics.CImatrixWithSymmetryJP(), JP, rotatedBasis, nm2, grid, asfSettings; printout=false)
    vals_rot, vecs_rot = eigen(Symmetric(Hmat_rot))
    iMin            = argmin(abs.(vals_rot .- targetLevel.energy))
    refMcSub        = vecs_rot[:, iMin]
    println("     matched eigenvalue = $(vals_rot[iMin])   (original level energy = $(targetLevel.energy))")

    println("\n>> Stage 3 prediction, with C = the rotation matrix directly (bypassing Stage 1)")
    shellList  = sort([sh for sh in basis.subshells if sh.kappa == sh3s.kappa], by = sh -> sh.n)
    n          = length(shellList)
    Qfull      = Matrix{Float64}(I, n, n)
    i3 = findfirst(==(sh3s), shellList);   i4 = findfirst(==(sh4s), shellList)
    Qfull[i3,i3] = Q11;   Qfull[i3,i4] = Q12;   Qfull[i4,i3] = Q21;   Qfull[i4,i4] = Q22
    transformationDirect = Dict(sh3s.kappa => (shellList, shellList, Qfull, Qfull))
    Mdirect         = BiOrthogonal.generateCounterRotatingCiMatrices(basis, transformationDirect, :left)
    predictedMcFull = Mdirect * targetLevel.mc
    predictedMcSub  = predictedMcFull[idxCsf];   predictedMcSub = predictedMcSub / norm(predictedMcSub)

    diffPlus  = maximum(abs.(predictedMcSub - refMcSub))
    diffMinus = maximum(abs.(predictedMcSub + refMcSub))
    println("     predicted mc = $predictedMcSub")
    println("     reference mc = $refMcSub")
    println("     best match (allowing for the usual overall-sign ambiguity of an eigenvector) = $(min(diffPlus,diffMinus))")
    #
end
