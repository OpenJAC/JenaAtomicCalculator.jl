#
# probe-secondOrderQ.jl
#
# A PROBE, not a module:  does a second-order Epstein-Nesbet treatment of the Q space reproduce the
# energy lowering that an exact CI in P (+) Q gives?  The difference between the two IS the higher-order
# remainder, so this measures directly whether "the coupling to Q is weak enough that second order suffices".
#
# Test case:  Cl III  1s^2 2s^2 2p^6 3s^2 3p^3,  J = 3/2 odd  --  the system of
#             G. Gaigalas, P. Rynkun, L. Kitoviene, Lith. J. Phys. 64(1), 20-39 (2024), Section 4.1.
#             Their orbitals cannot be reproduced here, so their numbers serve only as a scale check.
#
# Path (Rule 21):  AL . ref(3s^2 3p^3) + 3d,4s . all free . Coulomb
#
using JenaAtomicCalculator, LinearAlgebra, Printf
const JAC = JenaAtomicCalculator

println("="^110)
println("PROBE:  second-order Q-space correction against the exact CI,  Cl III  3s^2 3p^3  J = 3/2 odd")
println("="^110)

nm       = Nuclear.Model(17.)
refConf  = Configuration("1s^2 2s^2 2p^6 3s^2 3p^3")
sym      = LevelSymmetry(AngularJ64(3//2), Basics.minus)
# Which shells may be excited FROM decides which of Gaigalas' three papers we are reproducing:
#   "VV"  3s,3p        -> valence-valence           (his Lith.J.Phys. 65(1) 32 (2025))
#   "C"   2s,2p        -> core and core-core        (his 64(2) 73 and 64(3) 139 (2024))
#   "ALL" 2s,2p,3s,3p  -> all four classes at once, separated afterwards by where the holes sit
mode     = length(ARGS) > 0 ? ARGS[1] : "VV"
coreSh   = [Shell("2s"), Shell("2p")]
valSh    = [Shell("3s"), Shell("3p")]
fromSh   = mode == "VV" ? valSh : (mode == "C" ? coreSh : vcat(coreSh, valSh))
toSh     = length(ARGS) > 1 && ARGS[2] == "3d" ? [Shell("3d")] : [Shell("3d"), Shell("4s")]
println(">>> mode = $mode ;  exciting from $fromSh into $toSh")
step     = RasStep(RasStep(), seFrom=fromSh, seTo=deepcopy(toSh), deFrom=fromSh, deTo=deepcopy(toSh))

settings = AsfSettings(AsfSettings(); scField = Basics.ALField(), eeInteraction = CoulombInteraction(),
                                      eeInteractionCI = CoulombInteraction(), gridStopper = false)

println("\n>>> generating the P (+) Q configuration list ...")
# built exactly as Basics.generateBasis does it, de-duplicated by == rather than by hash
confS    = Basics.generateConfigurations([refConf], fromSh, toSh)
confD    = Basics.generateConfigurations(Basics.generateConfigurations([refConf], fromSh, toSh), fromSh, toSh)
confList = Configuration[]
for  confa in vcat([refConf], confS, confD)
    addTo = true
    for  confb in confList    if  confa == confb    addTo = false;  break   end   end
    if  addTo    push!(confList, confa)   end
end
println(">>> $(length(confList)) configurations")

# THE BOX (Rule 12).  Taken from the WHOLE list the estimate is set by 4s at Zeff = 1 and gives 96 a.u.,
# because it cannot tell a spectroscopic Rydberg shell from a CORRELATION shell of the same n and l (see
# its own docstring).  Here 3d and 4s are correlation orbitals and belong in the valence region, so the box
# is set by hand.  The 9.7 and 96 a.u. runs already agreed on the remainders to 0.1 %, so the conclusion
# does not depend on this choice; 20 a.u. simply holds both the valence and the correlation orbitals.
grid     = Basics.recommendedGrid([refConf], nm; rbox = 20.0, printout=true)

# ---------------------------------------------------------------------------------------------------------
# ORBITALS, IN TWO LAYERS.  An average-level field over the WHOLE correlation space optimizes the core for
# configurations that already carry a core hole (mean occ of 2s came out as 1.807), so the reference is not
# a stationary point, Brillouin's theorem does not protect the single excitations, and the ranking is then
# dominated by singles that are compensating for the orbitals rather than describing correlation.
#   Layer 1:  orbitals optimized for the REFERENCE configuration alone.
#   Layer 2:  those are FROZEN and only the correlation orbitals relax in their field.
# ---------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------
# ORBITALS -- THE MBPT CONSTRUCTION.  Two earlier attempts both put ~99 % of the ranking into a 2s -> 4s
# SINGLE excitation, which is orbital relaxation masquerading as correlation:
#   (a) an average-level field over the WHOLE space optimizes the core for configurations that already
#       carry a core hole (mean occ of 2s = 1.807);
#   (b) freezing the reference core and relaxing only 3d,4s made it WORSE (97.6 % -> 98.8 %), because the
#       AL field still shapes 4s to serve the core-hole configurations that dominate the average.
# The standard remedy is to let the occupied AND the virtual orbitals be eigenstates of ONE one-particle
# operator -- here the reference's own potential.  Then a single-excitation matrix element is only the
# difference between that local potential and the true Fock operator, rather than a full orbital rotation.
# ---------------------------------------------------------------------------------------------------------
println("\n>>> layer 1: orbitals optimized for the reference configuration alone ...")
primitives = Bsplines.generatePrimitives(grid)
mpRef      = SelfConsistent.performSCF([refConf], nm, grid, settings, printout=false)
refBas     = mpRef.levels[1].basis

println(">>> building the reference's own one-particle potential, and its FULL spectrum ...")
pot        = Basics.add( Nuclear.nuclearPotential(nm, grid),
                         Basics.computePotential(Basics.DFSField(), grid, refBas) )
basis0     = Basics.generateBasis(confList, [sym])
allOrbs    = Bsplines.generateOrbitals(basis0.subshells, pot, nm, primitives; printout=false)
basis      = Basis(true, basis0.NoElectrons, basis0.subshells, basis0.csfs, basis0.coreSubshells, allOrbs)
println(">>> orbital energies:  " * join([ @sprintf("%s %.4f", sh, allOrbs[sh].energy) for sh in basis0.subshells ], "   "))

# ---------------------------------------------------------------------------------------------------------
# The Hamiltonian over the whole P (+) Q space, in the CSF basis.
# ---------------------------------------------------------------------------------------------------------
cache    = InteractionStrength.XLCache()
Hup      = Hamiltonian.setupMatrix(sym, basis, nm, grid, settings, cache, printout=true)
n        = size(Hup, 1)
H        = Symmetric(Hup) * Matrix(1.0I, n, n)          # materialise the full symmetric matrix

# indices of the CSF carrying this symmetry, in the order setupMatrix used
idx      = [k for k = 1:length(basis.csfs) if basis.csfs[k].J == sym.J && basis.csfs[k].parity == sym.parity]
@assert length(idx) == n

# ---------------------------------------------------------------------------------------------------------
# Label every CSF by its NON-RELATIVISTIC configuration:  the jj partners of one nl shell are kept together,
# so a selection threshold can never split a spin-orbit multiplet.
# ---------------------------------------------------------------------------------------------------------
function nlKey(csf, subshells)
    d = Dict{Shell,Int64}()
    for  k = 1:length(subshells)
        occ = csf.occupation[k]
        if  occ > 0
            sh = Shell(subshells[k].n, Basics.subshell_l(subshells[k]))
            d[sh] = get(d, sh, 0) + occ
        end
    end
    return( join(sort([ "$(sh)^$(w)" for (sh,w) in d ]), " ") )
end

keys_ = [ nlKey(basis.csfs[i], basis.subshells) for i in idx ]
refKey = nlKey(basis.csfs[idx[1]], basis.subshells)   # provisional; fixed below
# the reference key is the one built from the reference configuration itself
refKey = join(sort([ "$(sh)^$(w)" for (sh,w) in refConf.shells if w > 0 ]), " ")

isP    = [ k == refKey for k in keys_ ]
P      = findall(isP);    Q = findall(.!isP)
println("\n>>> P = $(length(P)) CSF of  $refKey ;   Q = $(length(Q)) CSF in $(length(unique(keys_[Q]))) further nl-configurations")

if length(P) == 0    error("The reference configuration produced no CSF of this symmetry -- check refKey: $refKey vs $(unique(keys_))")   end

# ---------------------------------------------------------------------------------------------------------
# (a) exact CI in P (+) Q     (b) CI in P alone     (c) P + second-order Epstein-Nesbet correction from Q
# ---------------------------------------------------------------------------------------------------------
Eexact = eigvals(Symmetric(H))
FP     = eigen(Symmetric(H[P,P]));   EP = FP.values;   CP = FP.vectors
Hqq    = [ H[q,q] for q in Q ]
Vpq    = H[P,Q]                                    # the coupling block, |P| x |Q|

nLev   = min(3, length(EP))
println("\n" * "-"^110)
@printf("%-4s %18s %18s %18s %14s %14s %12s\n",
        "lev", "E(P only) [a.u.]", "E(P)+dE(2) [a.u.]", "E exact [a.u.]", "dE(2) [a.u.]", "remainder", "sum|c|^2")
println("-"^110)

results = []
for i = 1:nLev
    c      = CP[:,i]
    Viq    = vec(c' * Vpq)                          # <Psi_i | V | q >  for every q
    D      = EP[i] .- Hqq                           # Epstein-Nesbet denominators, per CSF
    dE     = sum(Viq.^2 ./ D)
    c2     = sum((Viq ./ D).^2)
    rem    = Eexact[i] - (EP[i] + dE)
    @printf("%-4d %18.9f %18.9f %18.9f %14.6e %14.6e %12.4e\n", i, EP[i], EP[i]+dE, Eexact[i], dE, rem, c2)
    push!(results, (EP[i], dE, Eexact[i], rem, c2, D, Viq))
end
println("-"^110)

# ---------------------------------------------------------------------------------------------------------
# The denominator distribution -- the second watch-item:  a near-zero D is an intruder, and says
# that its configuration belongs in P rather than in the perturbation.
# ---------------------------------------------------------------------------------------------------------
D1 = results[1][6]
println(@sprintf("\n>>> denominators for level 1:  min|D| = %.4e   median|D| = %.4e   max|D| = %.4e  [a.u.]",
                 minimum(abs.(D1)), sort(abs.(D1))[cld(length(D1),2)], maximum(abs.(D1))))
nSmall = count(abs.(D1) .< 0.1)
println(">>> CSF with |D| < 0.1 a.u. (intruder candidates):  $nSmall of $(length(D1))")

# ---------------------------------------------------------------------------------------------------------
# Ranking by nl-configuration, for level 1
# ---------------------------------------------------------------------------------------------------------
Viq1 = results[1][7]
contrib = Dict{String,Float64}()
for (j,q) in enumerate(Q)
    k = keys_[q];   contrib[k] = get(contrib, k, 0.0) + Viq1[j]^2 / D1[j]
end
total = sum(values(contrib))
println(@sprintf("\n>>> level 1:  total second-order lowering from Q = %.6e a.u.", total))
println("    NOT comparable with Gaigalas' -2.4830e-04 a.u.: that is one K' of a CORE-valence space (2s,2p -> 3d),")
println("    whereas this Q is VALENCE-valence (3s,3p -> 3d,4s).  Different physics, two orders of magnitude apart.")
# how many holes has this configuration made in the core, and how many in the valence?
refOcc = Dict{Shell,Int64}(sh => w for (sh,w) in refConf.shells)
function holes(k::String)
    nc = 0;  nv = 0
    for tok in split(k)
        m = match(r"^(\d+)([spdfg])\^(\d+)$", tok);   m === nothing && continue
        sh = Shell(parse(Int64, m[1]), findfirst(isequal(m[2][1]), ['s','p','d','f','g']) - 1)
        d  = get(refOcc, sh, 0) - parse(Int64, m[3])
        if  d > 0
            if      sh in coreSh   nc += d
            elseif  sh in valSh    nv += d
            end
        end
    end
    return( (nc, nv) )
end
label(nc,nv) = nc == 0 ? "VV" : (nc == 1 && nv == 0 ? "C " : (nc == 1 ? "CV" : "CC"))

byClass = Dict{String,Float64}()
for (k,v) in contrib
    (nc,nv) = holes(k);   cl = label(nc,nv);   byClass[cl] = get(byClass, cl, 0.0) + v
end
println("\n>>> by correlation class  (C = one core hole, CV = one core + valence, CC = two core holes, VV = valence only):")
for (cl,v) in sort(collect(byClass), by = x -> -abs(x[2]))
    println(@sprintf("    %-4s %16.6e a.u.  %9.2f%%", cl, v, 100*v/total))
end

println("\n>>> ranking by nl-configuration (descending |contribution|):")
println(@sprintf("    %-40s %16s %10s", "nl-configuration", "dE [a.u.]", "% of total"))
for (k,v) in sort(collect(contrib), by = x -> -abs(x[2]))
    (nc,nv) = holes(k)
    println(@sprintf("    %-38s %-4s %16.6e %9.2f%%", k, label(nc,nv), v, 100*v/total))
end

# ---------------------------------------------------------------------------------------------------------
# THE EFFECTIVE HAMILTONIAN IN P  (Gaigalas' mode 1 proper, his Eq. (37) with the symmetrised denominator
# of his Eq. (30) -- but with EPSTEIN-NESBET denominators: each Q-CSF's own diagonal element in place of a
# configuration average, and each P-CSF's own diagonal element in place of Ebar(K)).
#
#   Heff(a,b) = H(a,b)  -  SUM_q  V(a,q) V(q,b) * 1/2 [ 1/(Hqq - Haa) + 1/(Hqq - Hbb) ]
#
# The state-specific treatment above corrects each level in isolation, so the levels do not stay mutually
# orthogonal and every one of them claims the full lowering from the same Q-CSF.  Diagonalising Heff
# restores that orthogonality within P, and the question is how much of the remainder it recovers.
# ---------------------------------------------------------------------------------------------------------
Hpp   = H[P,P]
Hdiag = [ H[p,p] for p in P ]
Heff  = copy(Hpp)
for (j,q) in enumerate(Q)
    for a = 1:length(P),  b = 1:length(P)
        Heff[a,b] -= Vpq[a,j]*Vpq[b,j] * 0.5*( 1/(Hqq[j] - Hdiag[a]) + 1/(Hqq[j] - Hdiag[b]) )
    end
end
Eeff = eigvals(Symmetric(Heff))

println("\n" * "-"^110)
println("EFFECTIVE HAMILTONIAN IN P  vs  the state-specific second order  vs  the exact CI")
println("-"^110)
@printf("%-4s %17s %17s %17s | %13s %13s\n", "lev", "state-specific", "Heff in P", "exact CI",
        "rem (state)", "rem (Heff)")
for i = 1:nLev
    (EPi, dE, Ex, rm, c2, D, Viq) = results[i]
    @printf("%-4d %17.9f %17.9f %17.9f | %13.6e %13.6e\n", i, EPi+dE, Eeff[i], Ex, Ex-(EPi+dE), Ex-Eeff[i])
end
println("-"^110)
for i = 1:nLev
    (EPi, dE, Ex, rm, c2, D, Viq) = results[i]
    low = Ex - EPi
    @printf("   lev %d:  state-specific %+7.2f %%   Heff %+7.2f %%   of the exact lowering\n",
            i, 100*((EPi+dE)-Ex)/abs(low), 100*(Eeff[i]-Ex)/abs(low))
end

# ---------------------------------------------------------------------------------------------------------
# What does the AVERAGING cost?  Gaigalas must replace the per-CSF denominator by one configuration
# average per K', because his implicit route pulls it out in front of the orbital sums.  We build Q
# explicitly and need not.  Here both are computed on the SAME matrix, so the difference is the
# averaging alone -- everything else is identical.
#
# Ebar(K) is the (2J+1)-weighted mean of the diagonal elements over ALL CSF of the configuration, which
# is the trace identity;  restricted here to this symmetry block, so it is their option 1 as far as this
# block can see it.
# ---------------------------------------------------------------------------------------------------------
Ebar = Dict{String,Float64}()
wsum = Dict{String,Float64}()
for  (j,i) in enumerate(idx)
    k = keys_[j];   w = Basics.twice(basis.csfs[i].J) + 1.0
    Ebar[k] = get(Ebar, k, 0.0) + w * H[j,j];    wsum[k] = get(wsum, k, 0.0) + w
end
for  k in keys(Ebar)    Ebar[k] = Ebar[k] / wsum[k]    end

println("\n" * "-"^110)
println("THE COST OF THE CONFIGURATION AVERAGE  (same matrix, same couplings; only the denominator differs)")
println("-"^110)
@printf("%-4s %20s %20s %16s %12s\n", "lev", "dE(2) per-CSF [a.u.]", "dE(2) averaged [a.u.]", "difference", "rel.")
for i = 1:nLev
    (EPi, dE, Ex, rm, c2, D, Viq) = results[i]
    dEav = 0.0
    for (j,q) in enumerate(Q)
        dEav += Viq[j]^2 / (Ebar[refKey] - Ebar[keys_[q]])
    end
    @printf("%-4d %20.9e %20.9e %16.3e %11.2f%%\n", i, dE, dEav, dEav-dE, 100*(dEav-dE)/abs(dE))
end
println("-"^110)

println("\n" * "="^110)
println("DONE.  The 'remainder' column is what second order does NOT capture:  the higher-order")
println("       contribution of Q, which is exactly the quantity that decides whether P is large enough.")
println("="^110)
