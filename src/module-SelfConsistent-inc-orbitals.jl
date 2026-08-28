

# Shared machinery: orbital orthonormalization, the AL/EOL energy functional and the Fock matrix,
# reached by more than one of the SCF paths below.

"""
`SelfConsistent.orthonormalizeSameKappa(newOrbitals::Dict{Subshell, Orbital}, newbVectors::Dict{Subshell, Vector{Float64}},
                                              subshells::Array{Subshell,1}, primitives::Bsplines.Primitives,
                                              matrixB::Array{Float64,2})`
    ... applies a Loewdin SYMMETRIC orthogonalization to every group of subshells sharing the same kappa, treating
        them as a SIMULTANEOUSLY-valid set rather than the sequence in which they happen to have been refined --
        this is the explicit "orthonormalize the orbital set" step shown as its OWN, separate box (after, not
        instead of, the per-orbital "apply projection" step) in the block diagram of Zatsarinny \\& Froese
        Fischer, Comput. Phys. Commun. 202, 287-303 (2016), Fig. 1. The per-orbital
        Hamiltonian.projectHamiltonian step alone does not guarantee a simultaneously-consistent same-kappa set
        when subshells are refined one at a time with lagged (Jacobi-style) data from the previous outer
        iteration -- this step corrects that afterwards, once per outer iteration, independent of processing
        order.

        For a kappa-group of n>1 subshells with B-spline expansion coefficient vectors b_1...b_n (columns of a
        (nsL+nsS) x n matrix B) and S-metric overlap matrix M (M_ij = b_i' * matrixB * b_j), the orthonormalized
        set is B_new = B * M^(-1/2), computed via an eigendecomposition of the small (n x n) matrix M. Unlike
        Gram-Schmidt, this treats every member of the group symmetrically -- no member is privileged by
        processing order -- matching the paper's "simultaneous" variation of same-symmetry orbitals directly.
        Groups of size 1 (no same-kappa partner) are left untouched, already trivially orthonormal.
        A (orbitals::Dict{Subshell,Orbital}, bVectors::Dict{Subshell,Vector{Float64}}) tuple is returned.
"""
function orthonormalizeSameKappa(newOrbitals::Dict{Subshell, Orbital}, newbVectors::Dict{Subshell, Vector{Float64}},
                                       subshells::Array{Subshell,1}, primitives::Bsplines.Primitives,
                                       matrixB::Array{Float64,2})
    orbitalsOut = deepcopy(newOrbitals)
    bVectorsOut = deepcopy(newbVectors)

    for  kappa  in  unique( [sh.kappa for sh in subshells] )
        group = [sh for sh in subshells if sh.kappa == kappa]
        if  length(group) < 2    continue    end

        B = hcat( [newbVectors[sh] for sh in group]... )        # (nsL+nsS) x n
        M = transpose(B) * matrixB * B                          # n x n, symmetric overlap matrix of the group

        eig  = Basics.diagonalize(MatrixWithLinearAlgebra(), M)
        n    = length(group)
        Mm12 = zeros(n, n)
        for  k = 1:n
            vk    = eig.vectors[k]
            Mm12 .+= (1.0/sqrt(eig.values[k])) .* (vk * transpose(vk))
        end

        Bnew = B * Mm12
        for  (idx, sh)  in  enumerate(group)
            vec             = Bnew[:,idx]
            en              = newOrbitals[sh].energy
            bVectorsOut[sh] = vec
            orbitalsOut[sh] = Bsplines.generateOrbitalFromVector(sh, en, vec, primitives)
        end
    end

    return( orbitalsOut, bVectorsOut )
end


"""
`SelfConsistent.computeFunctional(coeffs1p::Array{Coefficient1p,1}, coeffs2p::Array{Coefficient2p,1},
                                        orbitals::Dict{Subshell, Orbital}, grid::Radial.Grid, potential::Radial.Potential)`
    ... computes the MCDHF energy functional using InteractionStrength.XL_CoulombKinkAware (kink-aware
        two-electron Slater integral) for the two-electron term. Used by the ALField/EOLField code line, cf.
        solveAverageLevelField. An energy::Float64 is returned.
"""
function computeFunctional(coeffs1p::Array{Coefficient1p,1}, coeffs2p::Array{Coefficient2p,1},
                                 orbitals::Dict{Subshell, Orbital}, grid::Radial.Grid, potential::Radial.Potential)
    energy = 0.

    # Collect one-electron contributions -- unchanged, no kink in this integral
    for  cf  in  coeffs1p
        energy = energy + cf.T * RadialIntegrals.GrantIab(orbitals[cf.a], orbitals[cf.b], grid, potential)
    end

    # Collect two-electron contributions via the kink-aware integral.  The integral depends only on the rank and
    # the four SUBSHELLS, while distinct angular coefficients -- one per contributing CSF pair -- repeatedly ask
    # for the same one with a different weight, so it is memoised for the duration of this call.  The orbitals
    # are fixed here by construction, which is what makes that exact rather than approximate.  The gain is modest
    # and confined to MULTI-configuration bases, which is where cross-CSF coefficients repeat a quadruple: 4.3%
    # on a four-configuration beryllium EOL (429.9 s against 449.3 s) and nothing measurable on single-
    # configuration Ti+ or W+, whose quadruples are nearly all distinct.  It is kept because the RAS line is
    # multi-configuration by construction, and because the memo is demonstrably exact -- both runs returned
    # E = -14.6142687294, identical to every digit.
    rkCache = Dict{NTuple{5,Any}, Float64}()
    for  cf  in  coeffs2p
        rk = get!(rkCache, (cf.nu, cf.a, cf.b, cf.c, cf.d)) do
                 InteractionStrength.XL_CoulombKinkAware(cf.nu, orbitals[cf.a], orbitals[cf.b],
                                                               orbitals[cf.c], orbitals[cf.d], grid)
             end
        energy = energy + cf.V * rk
    end

    return( energy )
end


"""
`SelfConsistent.rotateOrbitals(subshellList::Array{Subshell,1}, orbitals::Dict{Subshell, Orbital}, grid::Radial.Grid,
                               settings::AsfSettings)`
    ... rotates pairwise the orbitals to enhance the covergence of the SCF procedures; it follows that 
        rotation schemes of C.F. Fischer ... . The pairwise rotation of the orbitals of the same symmetry 
        is done from the inner-to-outer subshells. A rotation is not necessary if two subshells are filled
        of if both are fixed. For each rotation, a rotation parameter beta is determined.
        A pair (newOrbitals::Dict{Subshell, Orbital}, betaMax::Float64) us returned which contains the rotated
        orbitals as well as the maximum beta parameter that occurred during the pairwise rotation of the orbital set.
"""
function rotateOrbitals(subshellList::Array{Subshell,1}, orbitals::Dict{Subshell, Orbital}, grid::Radial.Grid, 
                        settings::AsfSettings)
    betaMax = 0.;    newOrbitals = orbitals
    
    @warn("Not yet properly implemented.")
    
    return( newOrbitals, betaMax)
end


"""
`SelfConsistent.computeTwoElectronV(subshell::Subshell, coeffs2p::Array{Coefficient2p,1},
                                           bVectors::Dict{Subshell, Vector{Float64}}, primitives::Bsplines.Primitives,
                                           tensorCaches::Dict{Int64, NTuple{3,RadialIntegrals.ScreenedPotentialCache}})`
    ... computes the direct and exchange two-electron potential matrix for the given subshell, taking
        bVectors::Dict{Subshell,Vector{Float64}} -- B-spline expansion coefficients -- as the SOLE, canonical
        per-subshell state, rather than a persistent Dict{Subshell,Orbital}. Every coeffs2p entry relevant to
        subshell involves exactly one "partner" subshell (the diagonal self-term is its own partner); a
        partner's tabulated form is built only as a disposable, read-only byproduct
        (Bsplines.generateOrbitalFromVector), cached within this call so repeated coefficients sharing
        the same partner do not re-evaluate it, and discarded once this function returns -- never re-fit back
        into bVectors, since bVectors are never derived FROM a tabulated form here in the first place; they
        come directly from diagonalization (solveAverageLevelField/solveOptimizedLevelField). Shared
        by the ALField and EOLField code lines. A (nsL+nsS) x (nsL+nsS) matrixV::Array{Float64,2} is
        returned.
"""
function computeTwoElectronV(subshell::Subshell, coeffs2p::Array{Coefficient2p,1},
                                    bVectors::Dict{Subshell, Vector{Float64}}, primitives::Bsplines.Primitives,
                                    tensorCaches::Dict{Int64, NTuple{3,RadialIntegrals.ScreenedPotentialCache}};
                                    directKernels::Dict{Tuple{Int64,Subshell,Subshell},Array{Float64,2}} =
                                                   Dict{Tuple{Int64,Subshell,Subshell},Array{Float64,2}}(),
                                    exchangeKernels::Dict{Tuple{Int64,Subshell},Array{Float64,2}} =
                                                     Dict{Tuple{Int64,Subshell},Array{Float64,2}}())
    nsL     = primitives.grid.nsL;        nsS = primitives.grid.nsS
    matrixV = zeros( nsL+nsS, nsL+nsS )
    partnerOrbitals = Dict{Subshell, Orbital}()
    partnerCVectors = Dict{Subshell, Vector{Float64}}()

    for  cf  in  coeffs2p
        if      subshell == cf.a  &&  subshell == cf.b  &&  subshell == cf.c  &&  subshell == cf.d
            partner = subshell;   pattern = :diagonal
        elseif  subshell == cf.a  &&  subshell == cf.c  &&  cf.b == cf.d
            partner = cf.b;       pattern = :direct
        elseif  subshell == cf.a  &&  subshell == cf.d  &&  cf.b == cf.c
            partner = cf.b;       pattern = :exchange
        elseif  subshell == cf.b  &&  subshell == cf.d  &&  cf.a == cf.c
            partner = cf.a;       pattern = :direct
        elseif  subshell == cf.b  &&  subshell == cf.c  &&  cf.a == cf.d
            partner = cf.a;       pattern = :exchange
        elseif  subshell == cf.a  &&  subshell == cf.b  &&  cf.c == cf.d  &&  cf.c != subshell
            # "Pair-excitation" pattern (a=b=subshell, c=d=partner, subshell!=partner): arises only from
            # genuine cross-CSF (off-diagonal) two-particle matrix elements, e.g. an EOL coefficient
            # combination mixing two different configurations (never produced by AL's single-CSF-average
            # scheme, hence never needed before). Algebraically IDENTICAL to the :exchange pattern above --
            # R^k(a,b,c,d) with (a,b,c,d)=(X,X,Y,Y) equals R^k(X,Y,Y,X) exactly, since the second
            # integration variable's density is just a product of two scalar functions (Y(r)*X(r) times
            # X(r)*Y(r) are the same product, whichever argument slot each factor nominally occupies) --
            # so it is handled by the identical XL_CoulombTensor call, unweighted, as :exchange.
            partner = cf.c;       pattern = :exchange
        elseif  subshell == cf.c  &&  subshell == cf.d  &&  cf.a == cf.b  &&  cf.a != subshell
            partner = cf.a;       pattern = :exchange
        else
            continue
        end
        # SpinAngular does not pre-filter its coefficients on the C^L parity rule, so combinations whose
        # reduced matrix element <a||C^nu||c> vanishes arrive here carrying a NONZERO V -- seven of them even
        # for Be 1s^2 2s^2 + 1s^2 2p^2. The radial kernels annihilate them, so skipping is numerically inert;
        # it only avoids building a matrix already known to be zero, and stops the kernels' own guards from
        # firing thousands of times per SCF. The rule depends on the pattern, which is exactly what the
        # 26-Jul note below got wrong: the direct call passes a = c and b = d, so it collapses to `nu even`,
        # whereas the exchange call pairs (subshell, partner) and needs l_sh + l_partner + nu even.
        if      pattern in [:diagonal, :direct]   &&   isodd(cf.nu)                                     continue
        elseif  pattern == :exchange              &&
                isodd( Basics.subshell_l(subshell) + Basics.subshell_l(partner) + cf.nu )               continue
        end
        # NOTE (26-Jul-2026): a previous `iseven(subshell_l + partner_l + cf.nu)` filter here, intended to
        # drop Breit-only coefficients, was REMOVED -- it used the parity of two unrelated orbitals' l-values
        # plus nu, which is not a real Breit/Coulomb selection rule and was silently discarding legitimate
        # direct (nu=0) Coulomb terms between different-l shells (e.g. 1s-2p, 2s-2p). This is why Ne/Ar were
        # wrong while He/Be (pure s-shells, where the bogus filter never triggered) validated fine. coeffs2p
        # here comes from Basics.ALField()'s own (Coulomb-only) angular coefficients, so no Breit terms are
        # present to filter in the first place. Confirmed: summing ALL entries reproduces Ne's literature
        # total energy to ~3e-3 Ha (residual consistent with import precision), vs -157.6 Ha when filtered.

        partnerOrb = get!(partnerOrbitals, partner) do
            Bsplines.generateOrbitalFromVector(partner, 0.0, bVectors[partner], primitives)
        end

        if      pattern == :diagonal
            matrixV = matrixV + 2 * cf.V *
                      InteractionStrength.XL_CoulombKinkAware(cf.nu, subshell, partnerOrb, subshell, partnerOrb,
                                                              primitives, directKernels)
        elseif  pattern == :direct
            matrixV = matrixV + cf.V *
                      InteractionStrength.XL_CoulombKinkAware(cf.nu, subshell, partnerOrb, subshell, partnerOrb,
                                                              primitives, directKernels)
        else    # :exchange
            # partnerOrb was built by Bsplines.generateOrbitalFromVector, which silently canonicalizes its
            # sign so that P is positive at small r (see that function's `wSign` step) -- a convention applied
            # to the TABULATED reconstruction only, never fed back into bVectors[partner] itself. Whenever that
            # canonicalization actually flips the sign, partnerOrb and the raw bVectors[partner] end up
            # oppositely signed, even though XL_CoulombTensor's "b"/"c" arguments (from partnerOrb) and its
            # "cVector" argument (from bVectors[partner]) are meant to represent the identical orbital -- an odd
            # total power of partner's sign then leaks into the exchange matrix element (confirmed empirically
            # this session: exchange terms flip sign in exact lockstep with this mismatch, while diagonal/direct
            # terms -- which only ever use partnerOrb, an even number of times -- stay sign-invariant as they
            # should). Fix: mirror generateOrbitalFromVector's own small-r wSign test on the raw bVector,
            # and pass a sign-matched cVector so it always agrees with partnerOrb's canonical sign.
            cVector = get!(partnerCVectors, partner) do
                bV      = bVectors[partner]
                wSign   = 0.0
                for  i = 1:nsL
                    Bi    = primitives.bsplinesL[i];   add = 1 - Bi.lower
                    for  j = Bi.lower:min(Bi.upper,30)   wSign = wSign + bV[i] * Bi.bs[j+add]   end
                end
                wSign < 0.  ?  -bV  :  bV
            end
            (cacheLL, cacheLS, cacheSS) = tensorCaches[cf.nu]
            matrixV = matrixV + cf.V *
                      InteractionStrength.XL_CoulombTensor(cf.nu, subshell, partnerOrb, partnerOrb,
                                                                  cVector, subshell,
                                                                  cacheLL, cacheLS, cacheSS, primitives,
                                                                  exchangeKernels)
        end
    end

    return( matrixV )
end


"""
`SelfConsistent.computeFockMatrix(subshell::Subshell, coeffs2p::Array{Coefficient2p,1},
                                         bVectors::Dict{Subshell, Vector{Float64}}, primitives::Bsplines.Primitives,
                                         nucPot::Radial.Potential, storage::Dict{String,Array{Float64,2}},
                                         occ::Float64, tensorCaches::Dict{Int64, NTuple{3,RadialIntegrals.ScreenedPotentialCache}})`
    ... computes the full per-subshell Hamiltonian matrix (one-electron + direct + exchange), matching the
        role of DBSR_HF's hf_matrix.f90: Bsplines.setupLocalMatrix provides the one-electron (kinetic +
        nuclear + rest-mass) block natively in B-spline basis (reused unchanged), and
        computeTwoElectronV provides the two-electron block, divided by the subshell's own mean
        occupation -- algebraically identical to DBSR_HF's own hfm=qsum(i)*dhl (one-electron scaled up)
        followed by hfm=hfm/qsum(i) (whole matrix scaled down) once expanded, given the two-electron
        coefficients already carry their full (not per-electron) occupation weighting (verified this
        session against the DBSR_HF av_energy_coef reference for all of Ne's coefficients). Shared by the
        ALField and EOLField code lines. A (nsL+nsS) x (nsL+nsS) matrix::Array{Float64,2} is returned.
"""
function computeFockMatrix(subshell::Subshell, coeffs2p::Array{Coefficient2p,1},
                                  bVectors::Dict{Subshell, Vector{Float64}}, primitives::Bsplines.Primitives,
                                  nucPot::Radial.Potential, storage::Dict{String,Array{Float64,2}},
                                  occ::Float64, tensorCaches::Dict{Int64, NTuple{3,RadialIntegrals.ScreenedPotentialCache}};
                                  coeffs2pUnscaled::Array{Coefficient2p,1}=Coefficient2p[],
                                  directKernels::Dict{Tuple{Int64,Subshell,Subshell},Array{Float64,2}} =
                                                 Dict{Tuple{Int64,Subshell,Subshell},Array{Float64,2}}(),
                                  exchangeKernels::Dict{Tuple{Int64,Subshell},Array{Float64,2}} =
                                                   Dict{Tuple{Int64,Subshell},Array{Float64,2}}())
    # occ == 0 would give 1/occ = Inf and, against the zero two-electron matrix such a subshell has,
    # Inf * 0 = NaN in every element -- a whole Fock matrix of NaN that only surfaces much later, and
    # then as a completely unrelated-looking complaint about the negative-energy continuum.  A subshell
    # with no occupation has no mean field to be refined in; the caller must decide what to do with it
    # (solveOptimizedLevelField carries it forward unchanged), so refuse here rather than return NaN.
    if  abs(occ) < 1.0e-12
        error("SelfConsistent.computeFockMatrix(): subshell $subshell has occupation $occ. There is no " *
              "mean field to define for an unoccupied subshell; the caller must skip it instead.")
    end
    matrix  = Bsplines.setupLocalMatrix(subshell.kappa, primitives, nucPot, storage)
    matrixV = computeTwoElectronV(subshell, coeffs2p, bVectors, primitives, tensorCaches;
                                  directKernels=directKernels, exchangeKernels=exchangeKernels)
    # coeffs2pUnscaled, when given, contributes WITHOUT the 1/occ scaling. Empty by default, so the
    # returned matrix is bit-for-bit what it always was unless a caller opts in.
    if  length(coeffs2pUnscaled) > 0
        matrixU = computeTwoElectronV(subshell, coeffs2pUnscaled, bVectors, primitives, tensorCaches;
                                      directKernels=directKernels, exchangeKernels=exchangeKernels)
        return( matrix + (1.0/occ) * matrixV + matrixU )
    end
    return( matrix + (1.0/occ) * matrixV )
end
