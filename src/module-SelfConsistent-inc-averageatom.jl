

# The finite-temperature AVERAGE-ATOM path. Not reached by SelfConsistent.performSCF -- it is driven
# from Plasma.AverageAtomScheme (module-Plasma-inc-average-atom.jl), and kept here beside the other
# SCF drivers because it iterates orbitals in the same way.

"""
`SelfConsistent.freeElectronDensity(mu::Float64, temp::Float64)`
    ... returns the number density of a FREE (ideal, uniform) electron gas at chemical potential mu and temperature
        temp, both in atomic units;  a value::Float64 [electrons per a_o^3] is returned.

        `n = (2 temp)^(3/2) / (2 pi^2) F_(1/2)(mu/temp)` with the Fermi-Dirac integral
        `F_(1/2)(eta) = int_0^inf sqrt(u) du / (1 + exp(u-eta))`.  This is the "exactly known free-particle
        density" that the Blenski-Ishikawa prescription adds back after the free reference has been subtracted
        from every partial wave.

        **THE INTEGRAL IS TAKEN IN u = x^2 AND NOT IN u**, which costs nothing and buys five digits:  the
        integrand carries a sqrt(u) whose DERIVATIVE is singular at the origin, so a quadrature in u converges
        slowly and misleadingly -- measured 08-Oct-2026, Simpson in u reached only 1.3e-04 and QuadGK in u with
        rtol = 1e-10 was itself wrong by 1.4e-05, i.e. it reported success.  In x the integrand is smooth and
        twelve digits are stable.
"""
function freeElectronDensity(mu::Float64, temp::Float64)
    eta  = mu / temp
    xMax = sqrt( 60. + max(eta, 0.) )
    wa,_ = QuadGK.quadgk(x -> 2x*x / (1. + exp(x*x - eta)), 0., xMax; rtol=1.0e-12)

    return( (2temp)^1.5 / (2pi^2) * wa )
end


"""
`SelfConsistent.averageAtomFreeWaveLimit(mu::Float64, temp::Float64, radiusBox::Float64)`
    ... returns the partial wave `lMaxFree::Int64` beyond which the free box states of the cell can no longer hold
        `tolerance` electrons, so that the sum over partial waves may stop there.

        The lowest state of partial wave l in a box of radius R lies at `eps_1 = j_(l,1)^2 / (2 R^2)` with
        `j_(l,1)` the first zero of the spherical Bessel function j_l, taken here from McMahon's expansion --
        accurate to 0.3 % for l >= 2 and 4 % at l = 0, far more than a stopping criterion needs.  The whole wave
        can hold at most `(2l+2) f(eps_1, mu, temp)` electrons, and the sum stops where that falls below the
        tolerance.

        **THE LIMIT IS PHYSICAL AND CANNOT BE GUESSED FROM lMax**, which is why it is computed rather than
        offered as a parameter with a default.  It scales as `k R` with the thermal momentum, so it is short for
        a cold dense cell and long for a hot dilute one:  measured 08-Oct-2026 for silicon, solid density at
        10 eV needs no free wave beyond l = 5, while 0.1 g/cm^3 at 100 eV was still gaining +0.09 electrons
        between l = 25 and l = 30.  A fixed default would have been silently too small for exactly the
        warm-dense conditions the model exists for.
"""
function averageAtomFreeWaveLimit(mu::Float64, temp::Float64, radiusBox::Float64; tolerance::Float64=1.0e-6,
                                  lCap::Int64=160)
    for  l = 0:lCap
        nu = l + 0.5
        jl = nu + 1.8557571 * nu^(1/3) + 1.0331500 * nu^(-1/3)
        if  (2l + 2) * Basics.FermiDirac( jl^2 / (2 * radiusBox^2), mu, temp ) < tolerance    return( l )    end
    end

    return( lCap )
end


"""
`SelfConsistent.averageAtomSpectra(pot::Radial.Potential, lMax::Int64, primitives::Bsplines.Primitives)`
    ... diagonalizes the one-electron problem of every kappa with l <= lMax in the given potential, and the
        FREE-PARTICLE problem of every kappa with l <= lMaxFree in the same box, and returns the NamedTuple
        `(kappas, kappasFree, lMax, lMaxFree, trueEigen, trueStart, freeEigen, freeStart)`:  the `Basics.Eigen`
        of each problem, keyed by kappa, and the index at which its physical (electron-branch) states begin.

        **THE WHOLE SPECTRUM OF EACH KAPPA IS KEPT, which is the point of this routine.**  The average-atom
        density used to be built from `Basics.generateShellList(1, nMax, lMax)`, i.e. `nMax - l` states per l --
        about SEVEN of the ~78 the B-spline basis offers for each kappa.  Johnson's continuum term is an integral
        over eps which in a box becomes a plain sum over the box states, the energy normalisation cancelling the
        level spacing, but ONLY if the sum runs over the WHOLE spectrum.

        **THE FREE SPECTRA ARE WHAT MAKES THE PARTIAL-WAVE SUM CONVERGE**, and they are computed to a HIGHER
        lMaxFree than the potential is solved to, deliberately.  A partial wave with l above lMax barely reaches
        into the potential, so its states are those of a free particle in the same box and need no
        diagonalization in the field;  the sum therefore runs over the true spectra up to lMax and over the free
        ones from lMax+1 to lMaxFree, which is exact as both grow and converges in lMax as fast as the difference
        `f(eps) - f(eps0)` falls away.  Measured 08-Oct-2026 on silicon, this is the whole difference between a
        usable and an unusable l-sum:  the plain sum of true states alone still moved by +0.40 electrons between
        lMax = 6 and 7 for a hot dilute plasma.

        **THE IDEAL-GAS CLOSED FORM IS NOT USED, and that was measured rather than assumed.**  The textbook
        Blenski-Ishikawa prescription subtracts the free partial waves WITHIN lMax and adds back `n_free(mu,T) V`
        in closed form.  That is exact only where the box gas reproduces the uniform gas, and a Wigner-Seitz cell
        is too small for it:  measured on silicon at four conditions, the cell's continuum holds 0.97, 0.82, 0.91
        and 0.94 of an ideal gas at the same mu, temperature and volume -- the Dirichlet depletion at the cell
        wall, 3 to 18 %, and NOT a small correction.  The closed form therefore converges in lMax to the wrong
        number rather than slowly to the right one, which is the dangerous of the two failures:  for solid-density
        silicon at 10 eV it gives Z* = 3.22 with 10.78 bound electrons, against 4.00 and 9.999 here, and 4.00 is
        the answer the shell structure demands -- the L shell closed and the M shell gone.  Summing the free box
        states instead costs one diagonalization per extra kappa and is right by construction.

        **THE FREE SPECTRA ARE THE SAME AT EVERY SCF ITERATION**, the potential being zero, so a converged run
        diagonalizes them once:  pass a previous call's result back as `freeSpectra` and only the true half is
        redone.
"""
function averageAtomSpectra(pot::Radial.Potential, lMax::Int64, primitives::Bsplines.Primitives;
                            lMaxFree::Int64=lMax, freeSpectra::Union{Nothing,NamedTuple}=nothing)
    grid    = primitives.grid;      npoints = grid.NoPoints
    nsL     = grid.nsL;             nsS     = grid.nsS
    if  lMaxFree < lMax    error("SelfConsistent.averageAtomSpectra(): lMaxFree = $lMaxFree is below lMax = $lMax.")   end
    storage = Dict{String,Array{Float64,2}}()
    matrixB = zeros(nsL+nsS, nsL+nsS)
    matrixB[1:nsL,1:nsL]                 = Bsplines.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
    matrixB[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = Bsplines.generateTTpMatrix!("SS-overlap", 0, primitives, storage)
    zeroPot = Radial.Potential("free particle", zeros(npoints), grid)
    kappas    = Int64[];                      kappasFree = Int64[]
    trueEigen = Dict{Int64,Basics.Eigen}();   trueStart  = Dict{Int64,Int64}()
    freeEigen = Dict{Int64,Basics.Eigen}();   freeStart  = Dict{Int64,Int64}()
    for  l = 0:lMaxFree,  kappa  in  (l == 0 ? (-1,) : (-l-1, l))
        push!(kappasFree, kappa)
        if  l <= lMax
            push!(kappas, kappa)
            wc = Bsplines.diagonalizeLocalMatrix(kappa, Bsplines.setupLocalMatrix(kappa, primitives, pot, storage),
                                                 matrixB, primitives)
            trueEigen[kappa] = wc;    trueStart[kappa] = Bsplines.findPositiveBranchStart(wc.values)
        end
        if  isnothing(freeSpectra)  ||  !haskey(freeSpectra.freeEigen, kappa)
            w0 = Bsplines.diagonalizeLocalMatrix(kappa, Bsplines.setupLocalMatrix(kappa, primitives, zeroPot, storage),
                                                 matrixB, primitives)
            freeEigen[kappa] = w0
            freeStart[kappa] = something( findfirst(e -> e > 0., w0.values), length(w0.values) + 1 )
        else
            freeEigen[kappa] = freeSpectra.freeEigen[kappa];    freeStart[kappa] = freeSpectra.freeStart[kappa]
        end
    end
    # THE SPURIOUS kappa > 0 STATE.  The Dirac problem in a finite B-spline basis carries one unphysical state for
    # each kappa > 0 -- the classical kinetic-balance artifact -- and for l >= 8 that state RISES ABOVE the
    # negative-continuum threshold, so `Bsplines.findPositiveBranchStart` hands back ITS index as the first
    # physical one.  Measured 08-Oct-2026 on a free-particle box of R = 3.218 a.u.:  kappa = -9 begins 8.0320,
    # 13.5931, 19.9907 while its partner kappa = 8 reads 0.0000, 8.0317, 13.5926 -- the same spectrum with one
    # extra state at zero energy prepended.  Counted as a state it is occupied with f(0) ~ 1/2, and the l-sum then
    # JUMPS by six electrons between lMax = 7 and 8 instead of converging.
    # It is found by the SPIN-ORBIT PARTNER rather than by a threshold:  kappa > 0 and kappa' = -(l+1) share the
    # same l and therefore very nearly the same radial spectrum, and the partner carries no artifact, so whichever
    # of the two leading states of kappa agrees with the partner's leading state is the physical one.  That test
    # needs no tolerance to be chosen and works in a real potential as well as in a box.
    for  (eigen, start)  in  ( (trueEigen, trueStart), (freeEigen, freeStart) )
        for  kappa  in  sort( collect(keys(eigen)) )
            if  kappa > 0  &&  haskey(start, -kappa-1)
                i0 = start[kappa];    ep = eigen[-kappa-1].values[ start[-kappa-1] ]
                if  i0 + 1 <= length(eigen[kappa].values)  &&
                    abs(eigen[kappa].values[i0+1] - ep) < abs(eigen[kappa].values[i0] - ep)
                    start[kappa] = i0 + 1
                end
            end
        end
    end

    return( (kappas=kappas, kappasFree=kappasFree, lMax=lMax, lMaxFree=lMaxFree, trueEigen=trueEigen,
             trueStart=trueStart, freeEigen=freeEigen, freeStart=freeStart) )
end


"""
`SelfConsistent.averageAtomElectronCount(mu::Float64, spectra::NamedTuple, temp::Float64)`
    ... counts the electrons of the Wigner-Seitz cell at chemical potential mu;  a triple
        `(nBound, nCont, nTotal)`::Tuple{Float64,Float64,Float64} is returned.

        `n = sum_(l<=lMax) (2j+1) sum_n f(eps_n)  +  sum_(lMax<l<=lMaxFree) (2j+1) sum_n f(eps0_n)`, i.e. the
        true box states where the potential was solved and the FREE box states of the same box above that -- the
        statement that a partial wave of high l does not reach into the atom.  A state is split into `nBound` or
        `nCont` by the sign of its energy, so that `nCont` is Johnson's mean charge Z* and, at neutrality,
        `nBound = Z - Z*`.

        Every box state is normalised over the cell, so the count needs no radial integral:  this is what the
        `int d eps` of the continuum term reduces to in a B-spline box.
"""
function averageAtomElectronCount(mu::Float64, spectra::NamedTuple, temp::Float64)
    nBound = 0.;    nCont = 0.
    for  kappa  in  spectra.kappas
        occ = 2 * abs(kappa);    wc = spectra.trueEigen[kappa]
        for  i = spectra.trueStart[kappa]:length(wc.values)
            wf = occ * Basics.FermiDirac(wc.values[i], mu, temp)
            if  wc.values[i] < 0.   nBound = nBound + wf   else   nCont = nCont + wf   end
        end
    end
    for  kappa  in  spectra.kappasFree
        if  abs(kappa) <= spectra.lMax + 1  &&  haskey(spectra.trueEigen, kappa)    continue    end
        occ = 2 * abs(kappa);    w0 = spectra.freeEigen[kappa]
        for  i = spectra.freeStart[kappa]:length(w0.values)
            nCont = nCont + occ * Basics.FermiDirac(w0.values[i], mu, temp)
        end
    end

    return( (nBound, nCont, nBound + nCont) )
end


"""
`SelfConsistent.averageAtomChemicalPotential(spectra::NamedTuple, temp::Float64, nm::Nuclear.Model)`
    ... solves the neutrality condition `nTotal(mu) = Z` on the full per-kappa spectrum;  a `chemMu::Float64` is
        returned.

        The spectra do not depend on mu, so this root solve is free compared with the diagonalizations that
        produced them:  one Fermi factor per box state per trial mu.  The count is strictly increasing in mu,
        which makes the root unique and bracketable, and the exit is on the RESIDUAL in electrons -- a quantity
        with a meaning -- rather than on the step.  The derivative is taken analytically from `df/dmu = f(1-f)/T`,
        which is the one place a Fermi-factor derivative has ever gone wrong here (see the note in
        `determineChemicalPotential`), so it is written once and used by both sums.
"""
function averageAtomChemicalPotential(spectra::NamedTuple, temp::Float64, nm::Nuclear.Model)
    g(mu::Float64)  = SelfConsistent.averageAtomElectronCount(mu, spectra, temp)[3] - nm.Z
    function gprime(mu::Float64)
        wa = 0.
        for  kappa  in  spectra.kappas
            occ = 2 * abs(kappa);    wc = spectra.trueEigen[kappa]
            for  i = spectra.trueStart[kappa]:length(wc.values)
                wf = Basics.FermiDirac(wc.values[i], mu, temp);    wa = wa + occ * wf * (1. - wf) / temp
            end
        end
        for  kappa  in  spectra.kappasFree
            if  abs(kappa) <= spectra.lMax + 1  &&  haskey(spectra.trueEigen, kappa)    continue    end
            occ = 2 * abs(kappa);    w0 = spectra.freeEigen[kappa]
            for  i = spectra.freeStart[kappa]:length(w0.values)
                wf = Basics.FermiDirac(w0.values[i], mu, temp);    wa = wa + occ * wf * (1. - wf) / temp
            end
        end
        return( wa )
    end
    epsLo = minimum( spectra.trueEigen[k].values[spectra.trueStart[k]]  for k in spectra.kappas )
    muLo  = epsLo - 60temp - 1.0;    muHi = max(1.0, 60temp);    nExp = 0
    while  g(muLo) > 0.  &&  nExp < 200    muLo = muLo - max(1.0, abs(muLo));    nExp = nExp + 1    end
    while  g(muHi) < 0.  &&  nExp < 400    muHi = muHi + max(1.0, abs(muHi));    nExp = nExp + 1    end
    if  g(muLo) > 0.  ||  g(muHi) < 0.
        error("SelfConsistent.averageAtomChemicalPotential(): no bracket for the neutrality condition between " *
              @sprintf("%.3e", muLo) * " and " * @sprintf("%.3e", muHi) * " Ha;  the electron count runs from " *
              @sprintf("%.3e", g(muLo) + nm.Z) * " to " * @sprintf("%.3e", g(muHi) + nm.Z) * " against Z = $(nm.Z).")
    end
    chemMu = 0.5 * (muLo + muHi);    gNow = g(chemMu);    nx = 0
    for  it = 1:200
        nx = it
        if  abs(gNow) < 1.0e-10 * max(1.0, nm.Z)    break    end
        if  gNow > 0.    muHi = chemMu    else    muLo = chemMu    end
        gp    = gprime(chemMu)
        newMu = gp > 0. ? chemMu - gNow / gp : 0.5 * (muLo + muHi)
        if  !(muLo < newMu < muHi)    newMu = 0.5 * (muLo + muHi)    end
        if  abs(newMu - chemMu) < 1.0e-14 * max(1.0, abs(chemMu))    chemMu = newMu;    break    end
        chemMu = newMu;    gNow = g(chemMu)
    end
    if  abs(gNow) > 1.0e-6 * max(1.0, nm.Z)
        Defaults.warn(AddWarning(), "SelfConsistent.averageAtomChemicalPotential(): the neutrality condition is " *
                      "satisfied only to " * @sprintf("%.2e", gNow) * " electrons after $nx iterations.")
    end

    return( chemMu )
end


"""
`SelfConsistent.averageAtomDensity(mu::Float64, spectra::NamedTuple, temp::Float64, primitives::Bsplines.Primitives)`
    ... accumulates the radial electron density of the cell from the full per-kappa spectrum;  a
        `rhot::Array{Float64,1}` is returned, normalised as `rhot(r) = 4 pi r^2 rho(r)`, so that
        `int_0^R rhot dr` is the electron number -- the convention `Basics.computePotential` uses for its
        average-atom fields.

        It is the r-space form of exactly the sum `averageAtomElectronCount` evaluates for the number:  the true
        box states up to lMax and the free box states of the same box above it.  The free states carry the
        Dirichlet depletion at the cell boundary, which is why this tail is right where a uniform gas added as a
        constant is not.

        **A STATE IS BUILT ONLY WHERE ITS OCCUPATION MATTERS.**  With ~50 kappa's and ~78 box states each,
        building every orbital would mean some 4000 B-spline reconstructions per iteration, almost all of them
        weighted by a Fermi factor below 1e-12;  the cutoff leaves the few hundred that carry the density.
"""
function averageAtomDensity(mu::Float64, spectra::NamedTuple, temp::Float64, primitives::Bsplines.Primitives)
    grid = primitives.grid;    npoints = grid.NoPoints;    nsL = grid.nsL;    nsS = grid.nsS
    rhot = zeros(npoints);     P = zeros(npoints);         Q = zeros(npoints)
    fMin = 1.0e-12
    branches = Tuple{Basics.Eigen,Int64,Int64}[]
    for  kappa  in  spectra.kappas        push!(branches, (spectra.trueEigen[kappa], spectra.trueStart[kappa], kappa))   end
    for  kappa  in  spectra.kappasFree
        if  abs(kappa) <= spectra.lMax + 1  &&  haskey(spectra.trueEigen, kappa)    continue    end
        push!(branches, (spectra.freeEigen[kappa], spectra.freeStart[kappa], kappa))
    end
    for  (eigen, i0, kappa)  in  branches
        occ = 2 * abs(kappa)
        for  i = i0:length(eigen.values)
            wf = Basics.FermiDirac(eigen.values[i], mu, temp)
            if  wf < fMin    continue    end
            ev = eigen.vectors[i];    fill!(P, 0.);    fill!(Q, 0.)
            for  k = 1:nsL
                bs = primitives.bsplinesL[k];   add = 1 - bs.lower
                for  j = bs.lower:bs.upper    P[j] = P[j] + ev[k] * bs.bs[j+add]    end
            end
            for  k = 1:nsS
                bs = primitives.bsplinesS[k];   add = 1 - bs.lower
                for  j = bs.lower:bs.upper    Q[j] = Q[j] + ev[nsL+k] * bs.bs[j+add]    end
            end
            wo = occ * wf
            for  j = 1:npoints    rhot[j] = rhot[j] + wo * (P[j]^2 + Q[j]^2)    end
        end
    end

    return( rhot )
end


"""
`SelfConsistent.averageAtomPotential(scField::Basics.AbstractScField, rhot::Array{Float64,1}, grid::Radial.Grid)`
    ... builds the ELECTRONIC part of the average-atom potential from a given radial density rhot, and returns a
        `pot::Radial.Potential` holding `Z(r) = -r V(r)`, as everywhere else in JAC.

        It is the density-space counterpart of `Basics.computePotential(::AaDFSField, grid, orbitals, mu, temp)`
        and uses the same two terms -- the direct `int dr' rhot(r')/r_>` and Slater exchange
        `-(3 rhot / (4 pi^2 r^2))^(1/3)` for `AaDFSField`, the Hartree-Slater pair for `AaHSField`.  It exists
        because the spectrum formulation has a density but no subshell list, and because the density no longer
        needs rebuilding once per kappa:  the quadratic direct term is now evaluated ONCE per SCF iteration
        rather than once per symmetry block, which for lMax = 7 is sixteen times less work.
"""
function averageAtomPotential(scField::Basics.AbstractScField, rhot::Array{Float64,1}, grid::Radial.Grid)
    npoints = grid.NoPoints;    wb = zeros(npoints);    wx = zeros(npoints)
    for  i = 1:npoints
        for  j = 1:npoints    wx[j] = rhot[j] / max( grid.r[i], grid.r[j] )    end
        wb[i] = RadialIntegrals.V0(wx, npoints, grid)
    end
    if      scField isa Basics.AaDFSField
        for  i = 1:npoints   wb[i] = wb[i] - (3 / (4pi^2 * grid.r[i]^2) * max(rhot[i], 0.))^(1/3)             end
    elseif  scField isa Basics.AaHSField
        for  i = 1:npoints   wb[i] = wb[i] - (3/2) * (3 / (4pi^2 * grid.r[i]^2) * max(rhot[i], 0.))^(1/3)     end
    else    error("SelfConsistent.averageAtomPotential(): no density-space form for $(typeof(scField)); the " *
                  "average-atom schemes are Basics.AaDFSField() and Basics.AaHSField().")
    end
    for  i = 2:npoints    wb[i] = - wb[i] * grid.r[i]    end;    wb[1] = 0.

    return( Radial.Potential("average-atom (density space)", wb, grid) )
end


"""
`SelfConsistent.determineChemicalPotential(orbitals::Dict{Subshell, Orbital}, temp::Float64, radiusWS::Float64,
                                           nm::Nuclear.Model, grid::Radial.Grid)`
    ... determines the chemical potential mu from the neutrality condition  Sum_i (2j_i+1) f(eps_i, mu, temp) = Z,
        where f is the Fermi-Dirac factor; a safeguarded Newton iteration inside a bracket is used, and a
        chemMu::Float64 is returned.

        Note: this general finite-temperature Fermi-Dirac root-finding utility was moved here from module Plasma
              (where it originated as `determineChemicalPotential`), since Plasma.perform(::AverageAtomScheme,
              ...) needs SelfConsistent.solveAverageAtomField below, and solveAverageAtomField itself needs this
              function internally at every SCF iteration; keeping it in Plasma would have made the two modules
              depend on each other circularly. Nothing here is Plasma-scheme-specific.
"""
function determineChemicalPotential(orbitals::Dict{Subshell, Orbital}, temp::Float64, radiusWS::Float64, nm::Nuclear.Model,
                                    grid::Radial.Grid)
    # THE FERMI FACTOR IS EVALUATED IN THE STABLE BRANCH, which is also what makes the derivative provable.
    # Writing 1/(exp(w)+1) directly overflows for large positive w, and the previous version avoided that by CLAMPING
    # w at 300 -- which silently changes the function whose root is being sought.  Choosing the branch instead is
    # exact everywhere, and then  df/dmu = f (1-f) / temp  follows in one line and cannot be mistyped.
    fermi(w::Float64)   = w > 0. ? exp(-w) / (1. + exp(-w))  :  1. / (1. + exp(w))
    function g(mu::Float64)
        wa = -nm.Z
        for  (k,v)  in orbitals
            wa = wa + (Basics.twice(Basics.subshell_j(k)) + 1) * fermi( (v.energy - mu) / temp )
        end
        return( wa )
    end
    function gprime(mu::Float64)
        wa = 0.
        for  (k,v)  in orbitals
            wf = fermi( (v.energy - mu) / temp )
            wa = wa + (Basics.twice(Basics.subshell_j(k)) + 1) * wf * (1. - wf) / temp
        end
        return( wa )
    end
    # THE ROOT IS UNIQUE AND BRACKETABLE: every Fermi factor increases with mu, so g is strictly increasing, g -> -Z
    # as mu -> -inf and g -> (sum of all 2j+1) - Z > 0 as mu -> +inf whenever the subshell set can hold Z electrons.
    # A bracket therefore always exists, and finding it first is what stops the iteration running away.
    #
    # WHY THAT MATTERS HERE, measured 08-Oct-2026 on the version this replaces.  Its Newton derivative read
    # `occ * wc^2 / temp / (wc+1)^2` with wc = exp(w) -- exp(2w) where the derivative needs exp(w).  Against a
    # central difference of g itself, at its own starting point mu = -0.1 for Si at T = 10 eV: finite difference
    # 6.267236e-04, correct analytic 6.267246e-04, and that expression 2.406385e-08 -- FOUR ORDERS too small.  The
    # first step went to mu = -2.49e+08 instead of -9.57e+03 and never recovered, ending near -6.4e+24, whereupon
    # every occupation is zero, the electron sum is 0 instead of Z, and the mean charge comes out as exactly Z.
    # The error was worst where it mattered: for a deeply bound orbital exp(w) is tiny and squaring it annihilates
    # the derivative.  A fudge of -0.0011 had been added to the result "for stability", which is what a misbehaving
    # root-finder looks like from the outside.
    epsLo = minimum( v.energy  for (k,v) in orbitals );   epsHi = maximum( v.energy  for (k,v) in orbitals )
    muLo  = epsLo - 60temp - 1.0;    muHi = epsHi + 60temp + 1.0
    nExp  = 0
    while  g(muLo) > 0.  &&  nExp < 200    muLo = muLo - max(1.0, abs(muLo));    nExp = nExp + 1    end
    while  g(muHi) < 0.  &&  nExp < 400    muHi = muHi + max(1.0, abs(muHi));    nExp = nExp + 1    end
    if  g(muLo) > 0.  ||  g(muHi) < 0.
        error("SelfConsistent.determineChemicalPotential(): no bracket for the neutrality condition between " *
              @sprintf("%.3e", muLo) * " and " * @sprintf("%.3e", muHi) * " Ha;  g = " * @sprintf("%.3e", g(muLo)) *
              " and " * @sprintf("%.3e", g(muHi)) * ".  The subshell set holds " *
              string(sum(Basics.twice(Basics.subshell_j(k)) + 1  for (k,v) in orbitals)) * " electrons against " *
              "Z = $(nm.Z);  raise nMax or lMax if that is less than Z.")
    end
    # A SAFEGUARDED NEWTON: the step is taken when it stays inside the bracket and bisection is taken when it does
    # not, so the iteration inherits Newton's speed without its ability to leave the interval.  The exit is on the
    # RESIDUAL of the neutrality condition itself -- electrons, a quantity with a meaning -- and not on the step.
    chemMu = 0.5 * (muLo + muHi);     nx = 0;    gNow = g(chemMu)
    for  it = 1:200
        nx = it
        if  abs(gNow) < 1.0e-10 * max(1.0, nm.Z)    break    end
        if  gNow > 0.    muHi = chemMu    else    muLo = chemMu    end
        gp    = gprime(chemMu)
        newMu = gp > 0. ? chemMu - gNow / gp : 0.5 * (muLo + muHi)
        if  !(muLo < newMu < muHi)    newMu = 0.5 * (muLo + muHi)    end
        if  abs(newMu - chemMu) < 1.0e-14 * max(1.0, abs(chemMu))    chemMu = newMu;    break    end
        chemMu = newMu;    gNow = g(chemMu)
    end
    if  abs(gNow) > 1.0e-6 * max(1.0, nm.Z)
        Defaults.warn(AddWarning(), "SelfConsistent.determineChemicalPotential(): the neutrality condition is " *
                      "satisfied only to " * @sprintf("%.2e", gNow) * " electrons after $nx iterations.")
    end
    println(">>> chemical potential: mu = " * @sprintf("%.8f", chemMu) * " Ha after $nx safeguarded Newton steps;  " *
            "the neutrality residual is " * @sprintf("%.2e", gNow) * " electrons.")

    return( chemMu )
end


"""
`SelfConsistent.solveAverageAtomSpectra(orbitals::Dict{Subshell, Orbital}, nuclearModel::Nuclear.Model,
                                         scField::Basics.AbstractScField, temp::Float64, radiusWS::Float64,
                                         primitives::Bsplines.Primitives; printout::Bool=true)`
    ... solves the average-atom field self-consistently IN DENSITY SPACE, i.e. on the full per-kappa B-spline
        spectrum rather than on a list of subshells, and returns the NamedTuple
        `(orbitals, chemMu, nBound, nCont, rhot, pot, spectra, lMax, radiusBox, accuracy, iterations)`.

        The cycle is  spectra -> chemical potential -> density -> potential -> spectra, each step the one named in
        Johnson's average-atom equations:  `averageAtomSpectra` for the box states of every kappa,
        `averageAtomChemicalPotential` for neutrality, `averageAtomDensity` for `4 pi r^2 rho(r)`, and
        `averageAtomPotential` for the field it generates.  The density is mixed linearly between iterations,
        which an average atom needs:  the continuum responds to the potential far more strongly than a bound
        orbital does, and an unmixed cycle oscillates between a too-ionized and a too-neutral solution instead of
        settling.

        **THE CELL RADIUS IS THE B-SPLINE BOX, NOT A SEPARATE PARAMETER**, and that is a physical statement rather
        than an implementation detail:  every state of the model, bound or free, is a box state normalised over
        the cell, and the electron count that neutrality is imposed on is the count in THAT volume.  A box larger
        than the Wigner-Seitz cell therefore solves a different cell and shifts the neutrality condition bodily;
        this routine computes with the box it is given and WARNS when the two disagree by more than 5 %.

        `orbitals` is used only for its subshell list -- which fixes lMax and which orbitals are handed back for
        the properties that follow -- and for the starting density of the first iteration.  `lMaxFree` is the
        partial wave to which the FREE box states stand in for the true ones;  left unset it is derived at every
        iteration from `averageAtomFreeWaveLimit`, i.e. from the temperature and the cell, and the printed table
        carries the value actually used.
"""
function solveAverageAtomSpectra(orbitals::Dict{Subshell, Orbital}, nuclearModel::Nuclear.Model,
                                 scField::Basics.AbstractScField, temp::Float64, radiusWS::Float64,
                                 primitives::Bsplines.Primitives; printout::Bool=true, lMaxFree::Int64=-1,
                                 freeWaveTolerance::Float64=1.0e-6)
    grid    = primitives.grid;      npoints = grid.NoPoints
    lMax    = maximum( Basics.subshell_l(k)  for (k,v) in orbitals )
    adaptFree = lMaxFree < lMax;    lMaxFree = max(lMaxFree, lMax)
    rBox    = grid.tL[end]
    if  abs(rBox - radiusWS) > 0.05 * radiusWS
        # Printed as well as collected: a collected warning is shown only where a caller flushes them, and this
        # one changes the cell the model solves rather than merely degrading it.
        println(">> The B-spline box R = " * @sprintf("%.4f", rBox) * " a.u. is NOT the Wigner-Seitz radius " *
                @sprintf("%.4f", radiusWS) * " a.u.:  the average atom is solved in a cell of volume " *
                @sprintf("%.1f", 4pi/3*rBox^3) * " a_o^3 instead of " * @sprintf("%.1f", 4pi/3*radiusWS^3) *
                " a_o^3, so its density is not the one asked for.  Use boxSize = R^(WS) when generating the grid.")
        Defaults.warn(AddWarning(), "SelfConsistent.solveAverageAtomSpectra(): the B-spline box R = " *
                      @sprintf("%.3f", rBox) * " a.u. is not the Wigner-Seitz radius " *
                      @sprintf("%.3f", radiusWS) * " a.u.;  every state is normalised over the BOX, so the " *
                      "cell this solves is the box.  Generate the grid with boxSize = R^(WS).")
    end
    nuclearPotential = Nuclear.nuclearPotential(nuclearModel, grid)
    # The starting density is the subshell-list one, which is all the hydrogenic orbitals can give; every later
    # iteration uses the spectrum.
    chemMu = determineChemicalPotential(orbitals, temp, radiusWS, nuclearModel, grid)
    rhot   = zeros(npoints)
    for  (k,v)  in  orbitals
        occ = (Basics.twice(Basics.subshell_j(k)) + 1) * Basics.FermiDirac(v.energy, chemMu, temp)
        for  i = 1:length(v.P)    rhot[i] = rhot[i] + occ * (v.P[i]^2 + v.Q[i]^2)    end
    end
    pot = Basics.add( nuclearPotential, averageAtomPotential(scField, rhot, grid) )
    #
    if  printout
        println("\n>> Average-atom SCF in density space:  lMax = $lMax, box R = " * @sprintf("%.4f", rBox) *
                " a.u., cell volume " * @sprintf("%.2f", 4pi/3*rBox^3) * " a_o^3" *
                (adaptFree ? ";  free partial waves to the tolerance " * @sprintf("%.0e", freeWaveTolerance) *
                             " electrons." : ";  free partial waves to l = $lMaxFree, as asked for."))
        println("   iter  lFree      mu [Ha]     n(bound)    n(cont)     n(total)    accuracy")
    end
    spectra = nothing;   nBound = 0.;   nCont = 0.;   accuracy = 1.0;   nx = 0;   eta = 0.4
    for  it = 1:60
        nx      = it
        spectra = averageAtomSpectra(pot, lMax, primitives; lMaxFree=lMaxFree, freeSpectra=spectra)
        chemMu  = averageAtomChemicalPotential(spectra, temp, nuclearModel)
        # The free tail is as long as the temperature and the cell make it, and mu is what says so, so the limit
        # is re-derived at every iteration.  It never shrinks: the spectra already computed are reused, and a sum
        # that lost a partial wave between two iterations would stop the SCF converging at all.
        if  adaptFree
            lNeeded = averageAtomFreeWaveLimit(chemMu, temp, rBox; tolerance=freeWaveTolerance)
            if  lNeeded > lMaxFree
                lMaxFree = lNeeded
                spectra  = averageAtomSpectra(pot, lMax, primitives; lMaxFree=lMaxFree, freeSpectra=spectra)
                chemMu   = averageAtomChemicalPotential(spectra, temp, nuclearModel)
            end
        end
        (nBound, nCont, nTotal) = averageAtomElectronCount(chemMu, spectra, temp)
        rhotNew  = averageAtomDensity(chemMu, spectra, temp, primitives)
        # The self-consistency is measured on the DENSITY, in electrons: it is the quantity the potential is built
        # from, and an error of 1e-6 electrons out of Z means the same thing at every temperature, which neither a
        # change in mu nor one in an orbital energy does.
        wx       = [ abs(rhotNew[i] - rhot[i])  for i = 1:npoints ]
        accuracy = RadialIntegrals.V0(wx, npoints, grid) / nuclearModel.Z
        for  i = 1:npoints    rhot[i] = (1. - eta) * rhot[i] + eta * rhotNew[i]    end
        pot      = Basics.add( nuclearPotential, averageAtomPotential(scField, rhot, grid) )
        if  printout
            println("   " * @sprintf("%4d", it) * "  " * @sprintf("%5d", lMaxFree) * " " *
                    @sprintf("%12.6f", chemMu) * " " *
                    @sprintf("%11.5f", nBound) * " " * @sprintf("%10.5f", nCont) * " " *
                    @sprintf("%11.6f", nTotal) * "  " * @sprintf("%.3e", accuracy))
        end
        if  accuracy < 1.0e-7    break    end
    end
    if  accuracy > 1.0e-5
        println(">> The average-atom field did NOT converge: it stopped at accuracy " *
                @sprintf("%.3e", accuracy) * " after $nx iterations ... computations proceed.")
        Defaults.warn(AddWarning(), "SelfConsistent.solveAverageAtomSpectra(): the density is self-consistent " *
                      "only to " * @sprintf("%.1e", accuracy) * " electrons per Z after $nx iterations.")
    end
    # Hand back the subshells that were asked for, reconstructed from the converged spectrum of their own kappa
    newOrbitals = Dict{Subshell, Orbital}()
    for  (k,v)  in  orbitals
        newOrbitals[k] = Bsplines.generateOrbitalFromPrimitives(k, spectra.trueEigen[k.kappa], primitives)
    end
    if  printout
        wIdeal = SelfConsistent.freeElectronDensity(chemMu, temp) * 4pi/3 * rBox^3
        println(">> mu = " * @sprintf("%.8f", chemMu) * " Ha;  " * @sprintf("%.5f", nBound) * " bound and " *
                @sprintf("%.5f", nCont) * " continuum electrons, i.e. a mean charge Z* = " *
                @sprintf("%.4f", nCont) * " after $nx iterations.")
        # How far the cell is from an ideal gas, which is what decides whether a closed-form free-electron term
        # would have done: a ratio near one says the continuum is essentially uniform, and the small ratios a
        # Wigner-Seitz cell actually gives are the Dirichlet surface term, not a correction.
        println("   An ideal electron gas at this mu and temperature would hold " * @sprintf("%.5f", wIdeal) *
                " electrons in the same volume, i.e. " * @sprintf("%.2f", nCont/max(wIdeal,1.0e-30)) *
                " times fewer than the cell does.")
    end

    return( (orbitals=newOrbitals, chemMu=chemMu, nBound=nBound, nCont=nCont, rhot=rhot, pot=pot,
             spectra=spectra, lMax=lMax, lMaxFree=lMaxFree, radiusBox=rBox, accuracy=accuracy,
             iterations=nx) )
end


"""
`SelfConsistent.solveAverageAtomField(orbitals::Dict{Subshell, Orbital}, nuclearModel::Nuclear.Model, scField::Basics.AbstractScField,
                                      temp::Float64, radiusWS::Float64, primitives::Bsplines.Primitives; printout::Bool=true)`
    ... solves the self-consistent field for a given local average-atom potential as specified by scField;
        a (new) set of `orbitals::Dict{Subshell, Orbital}` is returned.

        A thin wrapper over `solveAverageAtomSpectra`, which does the work and also returns the chemical
        potential, the electron counts and the converged density;  this form is kept because the subshell
        orbitals are all that the properties computed afterwards -- form factors, photoionization cross
        sections -- need from the field.
"""
function solveAverageAtomField(orbitals::Dict{Subshell, Orbital}, nuclearModel::Nuclear.Model, scField::Basics.AbstractScField,
                               temp::Float64, radiusWS::Float64, primitives::Bsplines.Primitives; printout::Bool=true)
    wa = solveAverageAtomSpectra(orbitals, nuclearModel, scField, temp, radiusWS, primitives; printout=printout)

    return( wa.orbitals )
end
