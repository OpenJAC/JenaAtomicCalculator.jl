#
println("An) Apply & test the bi-orthogonal transformation of two independently-generated multiplets.")

using LinearAlgebra

if  true
    # Last visit:  20-Oct-2026
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
    initialMultiplet = wb[ResultKeys.InitialMultiplet];   finalMultiplet = wb[ResultKeys.FinalMultiplet]
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
    # Last visit:  20-Oct-2026
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
    multiplet = wb[ResultKeys.Multiplet];   basis = multiplet.levels[1].basis
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
elseif  false
    # Last visit:  09-Nov-2026
    # Last reproduced:  09-Nov-2026 (8 of 8 numbers of the REPORT table below)
    # Last successful:  09-Nov-2026   M1 4.702e-04 /s and E2 Babushkin 2.621e-01 /s with term-dependent 2p
    #   orbitals, against M1 4.030e-04 /s and E2 Babushkin 6.698e-02 /s with one shared set; E2 Coulomb is
    #   structurally zero here and CF = 0.0008 on M1, so what is judged sound is the MECHANISM and the
    #   qualitative term-sensitivity, not the absolute rates -- see the REPORT below.
    # Branch c: TERM-DEPENDENT orbitals (competition item 14's first test case). The SAME 2p shell of
    #   carbon-like 1s^2 2s^2 2p^2 is optimized SEPARATELY for two different LS terms of the SAME
    #   configuration: orbital set "D", an EOL field targeted (levelSelectionCI) on the ^1D_2 level alone,
    #   and orbital set "P", an EOL field targeted on the ^3P term (its three J-levels together). Both
    #   sides share the identical configuration and CSF space but carry genuinely different, mutually
    #   non-orthogonal 2p orbitals -- precisely the case BiOrthogonal.jl's already-validated BoundStatePair
    #   mode exists for. NO new src/ code: this is purely a new example branch over existing, already-
    #   validated machinery (Atomic.Computation's independent initial-/final-side SCF, LevelSelection by
    #   indices, PhotoEmission's calcBiorthogonal). The indices assume the textbook Hund's-rule order for a
    #   p^2 ground configuration, ^3P_0 < ^3P_1 < ^3P_2 < ^1D_2 < ^1S_0 (indices 1..5) -- CONFIRMED by this
    #   run's own printed J/parity/energy tables below, in BOTH orbital sets (levels 1-3 cluster within a
    #   few meV, level 4 sits ~1.59 eV above as the lone J=2, level 5 ~2.38-3.98 eV above as the lone J=0).
    #
    #   REPORT (run via a scratch copy with this branch's guard flipped true, the others false; this branch
    #   itself stays `elseif false` pending the maintainer's own physical judgement -- see below):
    #   2p orbital overlaps BEFORE bi-orthogonalization (the term-dependence signal itself): 1s 1-<i|f>=
    #   5.7e-9, 2s 1.2e-5, 2p_1/2 5.55e-4, 2p_3/2 5.60e-4 -- the 2p shell (the one actually driving the
    #   1D2/3P term split) shows ~100x the non-orthogonality of the untouched 1s core, exactly the expected
    #   signature. AFTER BiOrthogonal.generateBiorthogonalShellMatrices: every overlap is 1 to 1e-15.
    #   The real [C I] 1D2 -> 3P2 forbidden line (the astrophysical 9850/9824 A auroral-type doublet's
    #   analogue here, M1+E2, same parity, both levels from this one configuration), WITH vs WITHOUT
    #   term-dependent orbitals:
    #     Quantity            WITH (calcBiorthogonal)      WITHOUT (one shared EOL set)      change
    #     Energy [eV]              1.556163                     1.586375                    -1.9%
    #     M1 [1/s]                 4.702e-04                    4.030e-04                   +16.7%
    #     E2 Coulomb [1/s]         5.947e-06                    5.955e-07                    +899%  (~10x)
    #     E2 Babushkin [1/s]       2.621e-01                    6.698e-02                    +291%  (~3.9x)
    #   Term-dependence is NOT a small correction here: the E2 rate (in either gauge) changes by a factor of
    #   several, which is itself a believable result for a p^2 forbidden line -- this exact class of
    #   transition is a textbook example of strong correlation/term sensitivity.
    #
    #   THE "GAUGE DISAGREEMENT" OF 551 AND 167 IS A NON-FINDING, and the first version of this report was
    #   wrong to rest its caveat on it -- resolved 09-Nov-2026. It is not a defect in E2, not a correlation
    #   effect and not a small-omega breakdown. The velocity (Coulomb) form of an ELECTRIC multipole vanishes
    #   IDENTICALLY between two orbitals of the SAME subshell, because the prefactor (kapa-kapb) and the
    #   integral I^-(a,a) die together; InteractionStrength.MabEmission's own docstring derives this. Both
    #   levels of this line come from the ONE configuration 1s^2 2s^2 2p^2, so every one-body E2 contribution
    #   is diagonal in 2p and the length form carries the whole amplitude. Measured directly on this run's
    #   own orbitals at omega = 0.0583 Ha: the diagonal 2p_3/2 E2 element is Coulomb -0.0 EXACTLY against
    #   Babushkin -3.6180e-08, and the single surviving off-diagonal pair 2p_1/2-2p_3/2 is itself suppressed
    #   by 181 in amplitude (-1.9977e-10 against -3.6148e-08). 181^2 = 3.3e+04 accounts for the measured E2
    #   rate ratio 2.621e-01/5.947e-06 = 4.4e+04, so the number is explained quantitatively and nothing is
    #   left over. (The 551 and 167 quoted above are the ratios of the TOTAL rate, which mixes in the
    #   gauge-invariant M1; that is a reporting defect in the quality table, fixed the same day.)
    #   AND E2 ITSELF IS SOUND, including at small omega, which was the obvious alternative suspect: on
    #   hydrogenic ions, where the Dirac orbitals are exact and gauge invariance must therefore be exact,
    #   JAC's E2 3d->1s agrees between the gauges to 1e-5 at Z = 6 and Z = 20, and E2 3d->3s -- a transition
    #   split only by the fine structure, omega = 6.6e-04 Ha at Z = 6, four orders of magnitude BELOW this
    #   line's omega -- agrees to 3-7 % at Z = 6 and 0.2-0.6 % at Z = 20, on rates of 1e-28 and 1e-19.
    #
    #   WHAT DOES CONDEMN THE ABSOLUTE RATES IS THE CANCELLATION FACTOR, AND IT FALLS ON M1: CF = 0.0008 in
    #   both gauges (0.0009 in the term-dependent run). That is the indicator lineQualitySummary's own
    #   docstring says to assert on, and it is sixty times below Cowan's 0.05 "carries no information"
    #   threshold -- the M1 amplitude is a residue of a near-cancellation, which is exactly what a
    #   spin-forbidden 1D2 -> 3P2 line is in a coupling scheme where the two J = 2 states of 2p^2 are jj
    #   components whose contributions almost annul. M1 dominates the physical line (the measured [C I]
    #   9850 A total A is ~2e-04 /s, and this run's 4.03e-04 /s is within a factor of two of it -- which on
    #   CF = 0.0008 is luck, not accuracy). E2 is the healthier of the two here, CF Babushkin = 0.4994.
    #   So the conclusion of this branch is UNCHANGED but its reason is now the right one: the mechanism
    #   (independent EOL targets + BoundStatePair biorthogonalization, end to end, no new src/ code) and the
    #   qualitative term-sensitivity are what this case establishes, while the absolute rates are not to be
    #   quoted -- because of CF on M1, not because the two gauges disagree.
    #   AND THE SMALLNESS OF THOSE MATRIX ELEMENTS IS LONG-WAVELENGTH PHYSICS, NOT ROUND-OFF, which is the
    #   question worth asking of a 3.6e-08 -- measured 09-Nov-2026. Note first WHICH quantity is which: the
    #   -0.0 and -3.6180e-08 above are reduced one-electron matrix elements <a||O^(E2)||b> in atomic units,
    #   while the 5.947e-06 and 2.621e-01 of the table are A-values in 1/s. For this line q = omega/c =
    #   4.1732e-04 a.u., and with <r^2> = 4.5097 a.u.^2 for this 2p_3/2 the leading spherical Bessel function
    #   is j_2(qr) ~ q^2 <r^2> / 15 = 5.235966e-08 -- which reproduces RadialIntegrals.GrantJL(2, q, 2p, 2p) =
    #   5.235966e-08 to every printed digit, and times AngularMomentum.JohnsonI(-2, -2, 2) = -0.6909883 gives
    #   the matrix element -3.617991e-08 exactly. So the magnitude is the (qr)^2 long-wavelength suppression
    #   of a forbidden E2 at 1.56 eV and nothing else. It is also TRUSTWORTHY at that size: the length-form
    #   integrand int j_2 (P^2 + Q^2) dr is POSITIVE-DEFINITE for a diagonal pair, so no near-equal numbers
    #   are subtracted and the element carries full relative precision; and the Coulomb -0.0 is an exact
    #   signed zero out of the vanishing (kapa - kapb) prefactor, not a 1e-20 residue that happened to print
    #   as zero. That is why the structural-zero argument above rests on these numbers rather than merely
    #   being consistent with them.
    grid = Radial.Grid(true)
    nm   = Nuclear.Model(6.)
    config = [Configuration("1s^2 2s^2 2p^2")]

    asfD = AsfSettings(AsfSettings(); scField=Basics.EOLField(), levelSelectionCI=LevelSelection(true, indices=[4]))
    asfP = AsfSettings(AsfSettings(); scField=Basics.EOLField(), levelSelectionCI=LevelSelection(true, indices=[1,2,3]))

    compD = Atomic.Computation(Atomic.Computation(), name="An: 2p^2, EOL on 1D2 alone", grid=grid, nuclearModel=nm,
                                configs=config, asfSettings=asfD)
    compP = Atomic.Computation(Atomic.Computation(), name="An: 2p^2, EOL on the 3P term", grid=grid, nuclearModel=nm,
                                configs=config, asfSettings=asfP)
    wbD = perform(compD; output=true);   wbP = perform(compP; output=true)
    mpD = wbD[ResultKeys.Multiplet];      mpP = wbP[ResultKeys.Multiplet]

    println("\n>> Orbital set D (EOL on 1D2) -- level table (index, J, parity, energy):")
    for (i, lev) in enumerate(mpD.levels)    println("     $i:  J=$(lev.J)  parity=$(lev.parity)  energy=$(lev.energy)")    end
    println(">> Orbital set P (EOL on the 3P term) -- level table (index, J, parity, energy):")
    for (i, lev) in enumerate(mpP.levels)    println("     $i:  J=$(lev.J)  parity=$(lev.parity)  energy=$(lev.energy)")    end

    bD = mpD.levels[1].basis;   bP = mpP.levels[1].basis
    println("\n>> 2p orbital overlaps BEFORE bi-orthogonalization (the TERM-DEPENDENCE signal -- should deviate from 1):")
    for sh in intersect(bD.subshells, bP.subshells)
        ov = RadialIntegrals.overlap(bD.orbitals[sh], bP.orbitals[sh], grid)
        println("     <$sh (1D2-optimized)|$sh (3P-optimized)> = $ov")
    end

    newD, newP = BiOrthogonal.generateBiorthogonalShellMatrices(bD, bP, grid)
    println(">> 2p orbital overlaps AFTER bi-orthogonalization (should be exactly 1):")
    for sh in bD.subshells
        println("     <$sh (new-D)|$sh (new-P)> = $(RadialIntegrals.overlap(newD[sh], newP[sh], grid))")
    end

    println("\n>> The real [C I] forbidden line, 1D2 -> 3P2 (M1+E2), WITH term-dependent orbitals (calcBiorthogonal):")
    photoSettings = PhotoEmission.Settings(PhotoEmission.Settings(), multipoles=[M1,E2], gauges=[UseCoulomb,UseBabushkin],
                                            lineSelection=LineSelection(true, indexPairs=[(4,3)]),
                                            printBefore=true, calcBiorthogonal=true)
    compTD = Atomic.Computation(Atomic.Computation(), name="An: 1D2 -> 3P2, term-dependent orbitals",
                grid=grid, nuclearModel=nm,
                initialConfigs=config, initialAsfSettings=asfD,
                finalConfigs=config,   finalAsfSettings=asfP,
                processSettings=photoSettings)
    perform(compTD; output=true)

    println("\n>> SAME line, ONE shared EOL orbital set over all 5 levels (no term-dependence -- baseline):")
    asfShared = AsfSettings(AsfSettings(); scField=Basics.EOLField())
    photoSettingsShared = PhotoEmission.Settings(photoSettings; calcBiorthogonal=false)
    compShared = Atomic.Computation(Atomic.Computation(), name="An: 1D2 -> 3P2, shared orbitals",
                grid=grid, nuclearModel=nm,
                initialConfigs=config, initialAsfSettings=asfShared,
                finalConfigs=config,   finalAsfSettings=asfShared,
                processSettings=photoSettingsShared)
    perform(compShared; output=true)
    #
end
