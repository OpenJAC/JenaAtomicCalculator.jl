
"""
`struct  Plasma.ConductivityLine`
    ... one dipole-allowed transition of the average atom, together with the conductivity it carries.

    + omega        ::Float64   ... transition energy eps_b - eps_a > 0 [Hartree].
    + sigma        ::Float64   ... the INTEGRATED conductivity of this transition [a.u.], i.e. its contribution
                                   to `int sigma(omega) d omega` and not a value of `sigma` at a point.
    + electrons    ::Float64   ... the number of electrons it carries in the conductivity sum rule,
                                   `(2 Omega/pi) sigma`.
    + kind         ::Symbol    ... `:boundBound`, `:boundFree` or `:freeFree`, by the signs of the two energies.
    + kappaA       ::Int64     ... kappa of the lower state.
    + kappaB       ::Int64     ... kappa of the upper state.
    + energyA      ::Float64   ... eps_a [Hartree].
    + energyB      ::Float64   ... eps_b [Hartree].
"""
struct  ConductivityLine
    omega          ::Float64
    sigma          ::Float64
    electrons      ::Float64
    kind           ::Symbol
    kappaA         ::Int64
    kappaB         ::Int64
    energyA        ::Float64
    energyB        ::Float64
end


# `Base.show(io::IO, line::Plasma.ConductivityLine)`  ... prepares a proper printout of line::Plasma.ConductivityLine.
function Base.show(io::IO, line::Plasma.ConductivityLine)
    println(io, "ConductivityLine($(line.kind)):  omega = " * @sprintf("%.5e", line.omega) * " Ha, sigma = " *
                @sprintf("%.5e", line.sigma) * " a.u., " * @sprintf("%.5e", line.electrons) * " electrons;  " *
                "kappa " * string(line.kappaA) * " -> " * string(line.kappaB) * ", eps " *
                @sprintf("%.4e", line.energyA) * " -> " * @sprintf("%.4e", line.energyB) * " Ha.")
end


"""
`Plasma.computeConductivityLines(scfResult::NamedTuple, temp::Float64, radiusWS::Float64, grid::Radial.Grid,
                                  primitives::Bsplines.Primitives)`
    ... computes every dipole-allowed transition between the box states of a converged average atom and the
        conductivity it carries;  an `Array{Plasma.ConductivityLine,1}` is returned, sorted by omega.

        This is the average-atom Kubo-Greenwood formula of Johnson, Guet & Bertsch, JQSRT 99, 327 (2006), Eq.
        (13), evaluated in the B-spline box.  Their Eqs. (14), (16) and (17) are the free-free, bound-bound and
        bound-free parts of the same sum;  here they are not separated, because in a box they are the same object
        -- a list of transitions between box states -- and only the SIGNS of the two energies say which name a
        transition goes by.  **THAT IS THE SIMPLIFICATION THE BOX BUYS.**  Johnson's free-free and bound-free
        terms carry an integral over the continuum energy with energy-normalised continuum functions;  in a box
        the energy normalisation cancels the level spacing and the integral becomes a plain sum over box states,
        exactly as it does for the density.

        Each transition contributes `sigma_ab = (pi/3 Omega) (f_a - f_b) omega_ab |<kappa_a||C^(1)||kappa_b>|^2
        R_ab^2`, with `R_ab = int (P_a P_b + Q_a Q_b) r dr` the length-form radial integral and `f` the Fermi
        factor.  The value stored is the INTEGRATED strength of the line, i.e. its contribution to
        `int sigma d omega`;  `Plasma.computeConductivitySpectrum` folds the lines into a spectrum.

        **THE BOX WALL SCATTERS, SO THE FREE-FREE PART NEEDS THE SAME SUBTRACTION THE DENSITY NEEDS.**  A free
        particle has no free-free absorption at all:  between two exact free waves normalised on the energy scale
        the dipole matrix element is a delta function in the energy, so it vanishes for any omega > 0, and
        Johnson's free-free l-sum converges because a high-l continuum pair barely feels the potential.  A BOX
        state is not a free wave -- the wall supplies the momentum -- so free-free transitions between box states
        are non-zero however high l goes, and the sum rule then grows without limit:  measured 11-Oct-2026 on
        aluminium at 0.27 g/cm^3, N_ff came out 3.13 at lMax = 3 and 6.82 at lMax = 5.  Called with
        `useFreeSpectra = true` this routine returns the same list for the FREE (V = 0) box, and the difference of
        the two is the physical free-free conductivity;  `Plasma.computeConductivitySumRules` and
        `Plasma.computeConductivitySpectrum` both take the two lists and subtract.

        **THE RELATIVISTIC ANGULAR FACTOR REPLACES JOHNSON'S `l_>`, AND REDUCES TO IT.**  Johnson's model is
        non-relativistic -- his Eq. (1) is the Schroedinger equation -- so his line strength is `l_> R^2`.  Summed
        over the two spin-orbit partners of each l, JAC's `AngularMomentum.CL_reduced_me` gives exactly `2 l_>`
        (checked for l = 0...4), and the factor of two is the spin that his prefactor carries separately.  So the
        form above is his wherever fine structure does not matter, and better where it does.

        Validated 11-Oct-2026 against their Table 2, aluminium at 0.27 g/cm^3 and T = 5 eV:  the four
        bound-bound lines come out at -0.24 %, -0.14 %, -0.42 % and +3.4 % (3s -> 3p, the smallest omega and the
        one most sensitive to the 3s and 3p energies), and the sum rule gives N_bb = 0.504 against their 0.493.
"""
function computeConductivityLines(scfResult::NamedTuple, temp::Float64, radiusWS::Float64, grid::Radial.Grid,
                                  primitives::Bsplines.Primitives; fMin::Float64=1.0e-8, sigmaMin::Float64=1.0e-14,
                                  useFreeSpectra::Bool=false)
    spectra = scfResult.spectra;    chemMu = scfResult.chemMu;    volume = 4pi/3 * radiusWS^3
    eigenOf = useFreeSpectra ? spectra.freeEigen : spectra.trueEigen
    startOf = useFreeSpectra ? spectra.freeStart : spectra.trueStart
    lines   = ConductivityLine[]
    # The radial functions are needed for every state that can take part, so they are built once and kept: a
    # transition needs both partners and every state appears in many pairs.
    orbs = Dict{Tuple{Int64,Int64}, Radial.Orbital}()
    occs = Dict{Tuple{Int64,Int64}, Float64}()
    for  kappa  in  spectra.kappas
        wc = eigenOf[kappa]
        for  i = startOf[kappa]:length(wc.values)
            wf = Basics.FermiDirac(wc.values[i], chemMu, temp)
            if  wf < fMin    break    end          ## the energies increase, so nothing beyond can contribute
            orbs[(kappa,i)] = Bsplines.generateOrbitalFromPrimitives(Subshell(100, kappa), wc.values[i],
                                                                     grid.NoPoints, wc.vectors[i], primitives)
            occs[(kappa,i)] = wf
        end
    end
    # ... and the upper partner of a transition may be empty, so it is built on demand and cached as well.
    function orbitalOf(kappa::Int64, i::Int64)
        if  !haskey(orbs, (kappa,i))
            wc = eigenOf[kappa]
            orbs[(kappa,i)] = Bsplines.generateOrbitalFromPrimitives(Subshell(100, kappa), wc.values[i],
                                                                     grid.NoPoints, wc.vectors[i], primitives)
            occs[(kappa,i)] = Basics.FermiDirac(wc.values[i], chemMu, temp)
        end
        return( orbs[(kappa,i)] )
    end
    #
    for  (keyA, orbA)  in  collect(orbs)
        kappaA, ia = keyA;    epsA = eigenOf[kappaA].values[ia];    fA = occs[keyA]
        lA = Basics.subshell_l( Subshell(100, kappaA) )
        for  kappaB  in  spectra.kappas
            lB = Basics.subshell_l( Subshell(100, kappaB) )
            if  abs(lA - lB) != 1    continue    end        ## E1 parity and triangle rule
            cl = AngularMomentum.CL_reduced_me(Subshell(100, kappaA), 1, Subshell(100, kappaB))
            if  cl == 0.    continue    end
            wc = eigenOf[kappaB]
            for  ib = startOf[kappaB]:length(wc.values)
                epsB = wc.values[ib];    omega = epsB - epsA
                if  omega <= 0.    continue    end
                fB = Basics.FermiDirac(epsB, chemMu, temp)
                if  fA - fB < fMin    continue    end
                orbB  = orbitalOf(kappaB, ib)
                rInt  = RadialIntegrals.rkDiagonal(1, orbA, orbB, grid)
                sigma = pi / (3 * volume) * (fA - fB) * omega * cl^2 * rInt^2
                if  sigma < sigmaMin    continue    end
                kind  = epsA < 0. ? (epsB < 0. ? :boundBound : :boundFree) : :freeFree
                push!(lines, ConductivityLine(omega, sigma, 2 * volume / pi * sigma, kind, kappaA, kappaB,
                                              epsA, epsB))
            end
        end
    end
    sort!(lines, by = x -> x.omega)

    return( lines )
end


"""
`Plasma.computeConductivitySumRules(lines::Array{Plasma.ConductivityLine,1}, nFree::Float64, radiusWS::Float64)`
    ... evaluates the conductivity sum rule of Johnson, Guet & Bertsch (2006), Eq. (15), separately for the three
        kinds of transition;  a NamedTuple `(nBoundBound, nBoundFree, nFreeFree, nTotal, nFree, gamma)` is
        returned.

        `int sigma(omega) d omega = pi N_free / (2 Omega)` in atomic units, so each transition carries
        `(2 Omega/pi) sigma_ab` electrons and the three partial sums are the effective electron numbers that
        Johnson quotes -- `N_bb = 0.493` for aluminium at 0.27 g/cm^3 and 5 eV.

        **THE SUM RULE IS WHAT FIXES THE ONE FREE PARAMETER OF THE MODEL.**  The free-free conductivity diverges
        as `1/omega^2` towards zero frequency, because the radial matrix element between two nearly degenerate
        continuum states grows without bound;  in reality a finite lifetime cuts it off.  Johnson mocks that up
        with a Drude factor `omega^2/(gamma^2 + omega^2)` and chooses `gamma` so that the FREE-FREE part of the
        sum rule returns the number of free electrons.  `gamma` is then the mean electron-ion collision rate, and
        the returned value is that choice:  `gamma` such that `sum_ff sigma_ab omega_ab^2/(gamma^2+omega_ab^2)`
        equals `pi n_free/(2 Omega)`.  For their aluminium case `gamma = 0.0466 a.u. = 1.27 eV`.
"""
function computeConductivitySumRules(lines::Array{Plasma.ConductivityLine,1}, nFree::Float64, radiusWS::Float64;
                                     freeLines::Array{Plasma.ConductivityLine,1}=ConductivityLine[])
    nbb = 0.;    nbf = 0.;    nff = 0.
    for  line  in  lines
        if      line.kind == :boundBound    nbb = nbb + line.electrons
        elseif  line.kind == :boundFree     nbf = nbf + line.electrons
        else                                nff = nff + line.electrons
        end
    end
    # the box-wall artifact, removed by the same free-particle reference the density uses
    nffFree = sum( [ line.electrons  for line in freeLines  if line.kind == :freeFree ], init=0. )
    nff     = nff - nffFree
    # gamma is the root of  sum_ff electrons_ab omega^2/(gamma^2+omega^2) = nFree, which decreases monotonically
    # in gamma from nff to zero;  so a bracket exists exactly when the undamped free-free sum exceeds nFree,
    # which is the divergence the damping is there to remove.
    ffLines = [ line  for line in lines      if line.kind == :freeFree ]
    ffFree  = [ line  for line in freeLines  if line.kind == :freeFree ]
    damped(g::Float64) = sum( [ line.electrons * line.omega^2 / (g^2 + line.omega^2)  for line in ffLines ], init=0. ) -
                         sum( [ line.electrons * line.omega^2 / (g^2 + line.omega^2)  for line in ffFree ], init=0. )
    gamma = 0.
    if  !isempty(ffLines)  &&  damped(0.) > nFree
        gLo = 0.;    gHi = 1.0
        while  damped(gHi) > nFree  &&  gHi < 1.0e4    gHi = 2gHi    end
        for  it = 1:200
            gamma = 0.5 * (gLo + gHi)
            if  damped(gamma) > nFree    gLo = gamma    else    gHi = gamma    end
            if  gHi - gLo < 1.0e-12 * max(1.0, gHi)    break    end
        end
    end

    return( (nBoundBound=nbb, nBoundFree=nbf, nFreeFree=nff, nTotal=nbb+nbf+nff, nFree=nFree, gamma=gamma,
             nFreeFreeBox=nff+nffFree, nFreeFreeArtifact=nffFree) )
end


"""
`Plasma.displayConductivitySumRules(sumRules::NamedTuple, lines::Array{Plasma.ConductivityLine,1})`
    ... displays the three partial sum rules, the collision rate taken from them and the strongest transitions of
        each kind;  nothing is returned.
"""
function displayConductivitySumRules(sumRules::NamedTuple, lines::Array{Plasma.ConductivityLine,1})
    println("\n  Conductivity sum rule, int sigma dw = pi N/(2 Omega), by kind of transition:\n")
    println("  " * "-"^74)
    println("   kind            transitions      N (electrons)")
    println("  " * "-"^74)
    for  (kind, nx)  in  ( (:boundBound, sumRules.nBoundBound), (:boundFree, sumRules.nBoundFree),
                           (:freeFree,   sumRules.nFreeFree) )
        nLines = count(x -> x.kind == kind, lines)
        println("   " * rpad(string(kind), 16) * lpad(string(nLines), 8) * "      " * @sprintf("%12.5f", nx))
    end
    println("  " * "-"^74)
    println("   total          " * lpad(string(length(lines)), 9) * "      " * @sprintf("%12.5f", sumRules.nTotal))
    println("   free electrons of the cell, for comparison:    " * @sprintf("%12.5f", sumRules.nFree))
    if  sumRules.gamma > 0.
        println("   Drude collision rate from the free-free sum rule:  gamma = " * @sprintf("%.5f", sumRules.gamma) *
                " a.u. = " * @sprintf("%.3f", Defaults.convertUnits("energy: from atomic to eV", sumRules.gamma)) *
                " eV")
    else
        println("   No Drude damping was needed:  the free-free sum rule does not exceed the free-electron count.")
    end
    println("  " * "-"^74)
    #
    for  kind  in  (:boundBound, :boundFree, :freeFree)
        wa = sort( [ line  for line in lines  if line.kind == kind ], by = x -> -x.sigma )
        if  isempty(wa)    continue    end
        println("\n   the strongest " * string(kind) * " transitions:")
        println("     omega [Ha]    omega [eV]     sigma [a.u.]    electrons    kappa_a -> kappa_b")
        for  line  in  wa[1:min(6,length(wa))]
            println("    " * @sprintf("%11.5f", line.omega) * "  " *
                    @sprintf("%11.3f", Defaults.convertUnits("energy: from atomic to eV", line.omega)) * "   " *
                    @sprintf("%14.7e", line.sigma) * "  " * @sprintf("%11.6f", line.electrons) * "      " *
                    lpad(string(line.kappaA), 4) * " -> " * lpad(string(line.kappaB), 4))
        end
    end
    println(" ")

    return( nothing )
end


"""
`Plasma.computeConductivitySpectrum(lines::Array{Plasma.ConductivityLine,1}, omegas::Array{Float64,1},
                                    broadening::Float64)`
    ... folds the conductivity lines into a spectrum on the given photon energies;  a
        `sigmas::Array{Float64,1}` of the real part of the conductivity [a.u.] is returned.

        Every transition of the box is a delta line, and a spectrum needs a width.  Each line is folded with a
        normalised Lorentzian of half-width `broadening`, `L(w) = broadening/(pi ((w-w_ab)^2 + broadening^2))`,
        which preserves `int sigma dw` exactly and so leaves the sum rules of
        `Plasma.computeConductivitySumRules` untouched.  Johnson, Guet & Bertsch (2006) do the same for their
        bound-bound lines, with the width taken to be the electron-ion collision rate.

        **THE WIDTH IS A MODEL AND NOT A RESULT**, which is the one thing to keep in view when reading a
        spectrum from this routine.  A broadening much smaller than the box level spacing shows the comb of box
        states rather than a continuum;  one much larger washes out real structure, such as the opening of an
        ionization edge.  The box spacing of the continuum near an energy eps is about `pi^2 n/R^2` with
        `n = R sqrt(2 eps)/pi`, so a sensible choice is a few times that.

        `freeLines`, if given, are folded in the same way and SUBTRACTED, which removes the box-wall artifact
        from the free-free part;  see `Plasma.computeConductivityLines`.
"""
function computeConductivitySpectrum(lines::Array{Plasma.ConductivityLine,1}, omegas::Array{Float64,1},
                                     broadening::Float64;
                                     freeLines::Array{Plasma.ConductivityLine,1}=ConductivityLine[],
                                     kinds::Array{Symbol,1}=[:boundBound, :boundFree, :freeFree])
    if  broadening <= 0.    error("Plasma.computeConductivitySpectrum(): the broadening must be positive.")   end
    sigmas = zeros( length(omegas) )
    for  (iw, omega)  in  enumerate(omegas)
        wa = 0.
        for  line  in  lines
            if  line.kind in kinds
                wa = wa + line.sigma * broadening / (pi * ((omega - line.omega)^2 + broadening^2))
            end
        end
        for  line  in  freeLines
            if  line.kind in kinds  &&  line.kind == :freeFree
                wa = wa - line.sigma * broadening / (pi * ((omega - line.omega)^2 + broadening^2))
            end
        end
        sigmas[iw] = wa
    end
    # A NEGATIVE CONDUCTIVITY IS NOT A ROUNDING ERROR AND IS NOT SILENCED HERE.  The free-free subtraction
    # removes the box wall by differencing two combs of lines, and the true potential shifts the continuum
    # levels, so the two combs have different spacings and their difference changes sign where neither is
    # smooth compared with the folding width.  That happens at small omega -- measured 11-Oct-2026 on aluminium
    # at 0.27 g/cm^3, Re sigma came out -1.2e-03 a.u. below 0.05 Ha -- and it is exactly the region the static
    # conductivity comes from.  So the free-free part of this spectrum is NOT usable towards zero frequency;
    # use `Plasma.drudeConductivity` there, which is a model with the sum rule fixing its one amplitude.
    if  minimum(sigmas) < 0.
        iBad = findall(x -> x < 0., sigmas)
        println(">> The conductivity is NEGATIVE at " * string(length(iBad)) * " of " * string(length(omegas)) *
                " photon energies, from " * @sprintf("%.4e", omegas[iBad[1]]) * " to " *
                @sprintf("%.4e", omegas[iBad[end]]) * " Ha (most negative " * @sprintf("%.3e", minimum(sigmas)) *
                " a.u.):  the free-free subtraction differences two combs of box lines and does not survive " *
                "towards zero frequency.  Read the bound-bound and bound-free parts there, and take the " *
                "free-electron term from Plasma.drudeConductivity.")
        Defaults.warn(AddWarning(), "Plasma.computeConductivitySpectrum(): the conductivity is negative below " *
                      @sprintf("%.3e", omegas[iBad[end]]) * " Ha;  the free-free subtraction does not survive " *
                      "towards zero frequency.")
    end

    return( sigmas )
end


"""
`Plasma.transportCrossSection(scfResult::NamedTuple, energy::Float64, grid::Radial.Grid; lMax::Int64=8)`
    ... returns the momentum-transfer (transport) cross section an electron of the given kinetic energy sees in
        the converged average-atom potential;  a `value::Float64` [a.u.^2] is returned, or 0 where the continuum
        orbital cannot be formed.

        `sigma_tr(k) = (4 pi/k^2) sum_l (l+1) sin^2(delta_{l+1} - delta_l)`, the standard phase-shift form, with
        the phases taken from `Continuum.generateOrbitalLocalPotential` in the self-consistent potential.

        **ONLY THE DIFFERENCE OF NEIGHBOURING PHASES ENTERS, AND THAT IS WHAT MAKES THIS WELL POSED.**  A phase
        shift is defined modulo pi and the routines that compute one do wrap:  on aluminium at 0.27 g/cm^3 and
        5 eV, delta_0 runs -1.447 at 0.2 Ha and +1.310 at 0.35 Ha, a jump of very nearly pi.  Because `sin^2` is
        pi-periodic the cross section does not notice, so no unwrapping is needed and a constant offset in the
        phase convention cannot matter.  Contrast the free-free dipole amplitude, where the box boundary enters
        the integrand itself;  that is the difference between this route to the DC conductivity and the one that
        failed.
"""
function transportCrossSection(scfResult::NamedTuple, energy::Float64, grid::Radial.Grid; lMax::Int64=8)
    energy <= 0.  &&  return( 0. )
    cSet  = Continuum.Settings(false, grid.NoPoints - 10)
    phase = fill(NaN, lMax + 2)
    for  l = 0:lMax+1
        try
            phase[l+1] = redirect_stdout(devnull) do
                wa, wp = Continuum.generateOrbitalLocalPotential(energy, Subshell(200, -(l+1)), scfResult.pot, cSet)
                wp
            end
        catch;    phase[l+1] = NaN
        end
    end
    k = sqrt(2*energy);    wa = 0.
    for  l = 0:lMax
        if  isnan(phase[l+1])  ||  isnan(phase[l+2])    continue    end
        wa = wa + (l + 1) * sin( phase[l+2] - phase[l+1] )^2
    end

    return( 4pi / k^2 * wa )
end


"""
`Plasma.dcConductivity(scfResult::NamedTuple, temp::Float64, radiusWS::Float64, grid::Radial.Grid;
                        lMax::Int64=8, energyMin::Float64=0.001, energyMax::Float64=2.5, noEnergies::Int64=48)`
    ... computes the DC electrical conductivity of the average atom from its own scattering phase shifts, with no
        adjustable parameter;  a NamedTuple `(sigma0, gamma, crossSections, energies)` is returned, `sigma0` in
        a.u. and `gamma` the collision rate the same sum rule assigns to it.

        **WHY NOT FROM THE FREE-FREE CONDUCTIVITY, which is where one would look first.**  A single-site model has
        no relaxation time:  one scattering centre cannot give a finite DC conductivity, so the small-frequency
        divergence of its free-free Kubo-Greenwood sum -- and the l-sum that grows with it -- is the physics and
        not a defect to be removed.  A finite `sigma(0)` needs the DENSITY of scatterers, which enters here as
        `n_ion = 1/Omega` through the relaxation time and nowhere else.  Johnson, Guet & Bertsch (2006) reach
        their own `sigma(0)` the same way, cross-checking it against the Ziman formula.

        In the relaxation-time approximation, with one ion per cell and `v = k` in atomic units,

            tau(k)   = 1 / ( n_ion k sigma_tr(k) )
            sigma(0) = (1/3 pi^2) int dk k^4 tau(k) (-df/d eps)|_{eps = k^2/2}

        and the prefactor is fixed rather than fitted:  for a constant `tau` the integral reduces to
        `sigma(0) = n_e tau` exactly, which is how it was checked.

        Validated against Johnson, Guet & Bertsch (2006) for aluminium at 0.27 g/cm^3 and T = 5 eV:  this gives
        `sigma(0) = 0.02763 a.u.` against their 0.0262, i.e. 5.5 % high, and the collision rate it implies is
        0.0486 a.u. against their 0.0466.  Converged to 0.1 % in `lMax` (6 -> 8), to 0.01 % in the mesh (24 -> 48
        points) and to 0.02 % in the lower energy cut (0.005 -> 0.001 Ha).  **The agreement is the same few per
        cent the Drude form reached while being GIVEN their collision rate;  what has changed is that nothing is
        given.**
"""
function dcConductivity(scfResult::NamedTuple, temp::Float64, radiusWS::Float64, grid::Radial.Grid;
                        lMax::Int64=8, energyMin::Float64=0.001, energyMax::Float64=2.5, noEnergies::Int64=48)
    volume = 4pi/3 * radiusWS^3;    nIon = 1 / volume;    chemMu = scfResult.chemMu
    # the mesh is uniform in k, since the integrand is written in k and carries k^4
    kMin = sqrt(2*energyMin);       kMax = sqrt(2*energyMax)
    ks   = [ kMin + (kMax - kMin) * (i-1) / (noEnergies - 1)   for i = 1:noEnergies ]
    es   = [ k^2/2  for k in ks ];      sigTr = zeros(noEnergies);      wa = zeros(noEnergies)
    for  i = 1:noEnergies
        sigTr[i] = Plasma.transportCrossSection(scfResult, es[i], grid; lMax=lMax)
        wx       = (es[i] - chemMu) / temp
        dfde     = wx > 60. ? 0. : exp(wx) / (temp * (1 + exp(wx))^2)
        wa[i]    = sigTr[i] <= 0. ? 0. : ks[i]^4 / (nIon * ks[i] * sigTr[i]) * dfde
    end
    wb = 0.;    for i = 1:noEnergies-1   wb = wb + 0.5*(wa[i] + wa[i+1])*(ks[i+1] - ks[i])   end
    sigma0 = wb / (3 * pi^2)
    gamma  = sigma0 > 0. ? scfResult.nCont / (volume * sigma0) : 0.

    return( (sigma0=sigma0, gamma=gamma, crossSections=sigTr, energies=es) )
end


"""
`Plasma.drudeConductivity(nFree::Float64, radiusWS::Float64, gamma::Float64, omegas::Array{Float64,1})`
    ... returns the free-electron (Drude) conductivity `sigma(w) = sigma(0) gamma^2/(gamma^2 + w^2)` of the cell
        on the given photon energies;  a `sigmas::Array{Float64,1}` [a.u.] is returned.

        **THIS IS A MODEL WITH ONE PARAMETER, AND IT IS HERE BECAUSE THE BOX CANNOT DO THE JOB.**  The free-free
        conductivity of a plasma carries the whole free-electron oscillator strength, and it carries it at photon
        energies of order the electron-ion collision rate -- below the level spacing of any affordable box, and
        in a region where the free-free subtraction of `Plasma.computeConductivitySpectrum` changes sign.  So the
        free-electron response is supplied here in closed form instead of computed.

        Its AMPLITUDE is not free:  the conductivity sum rule `int sigma dw = pi N_free/(2 Omega)` gives
        `sigma(0) = N_free/(Omega gamma)` with no adjustable factor left.  Only `gamma` is a choice -- and since
        11-Oct-2026 it need not be a choice at all:  `Plasma.dcConductivity` computes it from the average atom's
        own scattering phase shifts, so the whole optical chain can be run with nothing supplied from outside.  Johnson, Guet & Bertsch (2006) obtain
        `gamma = 0.0466 a.u. = 1.27 eV` for aluminium at 0.27 g/cm^3 and 5 eV by imposing that same sum rule on
        their (divergent) free-free integral, and quote `sigma(0) = 0.0262 a.u. = 0.121 (mOhm m)^-1` against
        0.125 from the Ziman formula.  The closed form above gives 0.0286 at their gamma -- 9 % high, the
        difference being that their value comes from the damped Kubo-Greenwood integrand and not from a pure
        Drude shape.
"""
function drudeConductivity(nFree::Float64, radiusWS::Float64, gamma::Float64, omegas::Array{Float64,1})
    if  gamma <= 0.    error("Plasma.drudeConductivity(): the collision rate gamma must be positive.")   end
    volume = 4pi/3 * radiusWS^3;    sigma0 = nFree / (volume * gamma)

    return( [ sigma0 * gamma^2 / (gamma^2 + omega^2)  for omega in omegas ] )
end


"""
`Plasma.computeOpticalProperties(omegas::Array{Float64,1}, sigmas::Array{Float64,1})`
    ... turns the real part of the conductivity into the optical constants of the plasma;  a NamedTuple
        `(sigmaImag, epsReal, epsImag, refIndex, extinction, absorption)` of `Array{Float64,1}` is returned, each
        on the given photon energies.

        The imaginary part of the conductivity follows from the real part by the Kramers-Kronig relation,
        `Im sigma(w) = -(2w/pi) P int_0^inf Re sigma(w') dw'/(w'^2 - w^2)`, and then, in atomic units and
        following Johnson, Guet & Bertsch (2006), Section 4.1,

            Re eps = 1 - 4 pi Im sigma / w,     Im eps = 4 pi Re sigma / w,     n + i k = sqrt(eps),

        with the refractive index and the extinction coefficient
        `n = sqrt((|eps| + Re eps)/2)`, `k = sqrt((|eps| - Re eps)/2)`, and the absorption coefficient
        `alpha = 2 w k / c` [1/a_o].

        **THE PRINCIPAL VALUE IS TAKEN BY SUBTRACTION AND IS THEREFORE EXACT RATHER THAN REGULARISED.**
        `P int_0^inf dw'/(w'^2-w^2)` vanishes identically, so subtracting `Re sigma(w)` from the numerator
        leaves an integrand with no pole and changes nothing;  at `w' = w` the integrand is the derivative
        `sigma'(w)/(2w)`, which is taken from the neighbouring points.

        **AND THE ANSWER IS ONLY AS GOOD AS THE RANGE OF `omegas`.**  The dispersion integral runs to infinity,
        so a grid that stops below the K edge of the element silently omits the bound-free strength above it and
        the refractive index comes out wrong at every frequency, not only near the edge.  Carry `omegas` well
        past the deepest ionization threshold, and log-spaced, since the integrand spans decades.
"""
function computeOpticalProperties(omegas::Array{Float64,1}, sigmas::Array{Float64,1})
    nw = length(omegas)
    if  nw != length(sigmas)    error("Plasma.computeOpticalProperties(): omegas and sigmas differ in length.")  end
    alpha = Defaults.getDefaults("alpha");    cLight = 1.0/alpha
    sigmaImag = zeros(nw);   epsReal = zeros(nw);   epsImag = zeros(nw)
    refIndex  = zeros(nw);   extinct = zeros(nw);   absorb  = zeros(nw)
    # the subtracted principal-value integral, on the trapezoidal rule of the given (possibly log-spaced) grid
    for  iw = 1:nw
        w  = omegas[iw];    wa = 0.
        for  jw = 1:nw-1
            w1 = omegas[jw];      w2 = omegas[jw+1]
            g1 = jw   == iw ? 0. : (sigmas[jw]   - sigmas[iw]) / (w1^2 - w^2)
            g2 = jw+1 == iw ? 0. : (sigmas[jw+1] - sigmas[iw]) / (w2^2 - w^2)
            # at the pole the subtracted integrand is sigma'(w)/(2w);  a one-sided difference is enough, the
            # neighbouring interval carrying the same value with the opposite offset
            if      jw   == iw  &&  iw < nw
                g1 = (sigmas[iw+1] - sigmas[iw]) / (omegas[iw+1] - omegas[iw]) / (2w)
            elseif  jw+1 == iw  &&  iw > 1
                g2 = (sigmas[iw] - sigmas[iw-1]) / (omegas[iw] - omegas[iw-1]) / (2w)
            end
            wa = wa + 0.5 * (g1 + g2) * (w2 - w1)
        end
        sigmaImag[iw] = -(2w/pi) * wa
        epsReal[iw]   = 1.0 - 4pi * sigmaImag[iw] / w
        epsImag[iw]   =       4pi * sigmas[iw]    / w
        wabs          = sqrt( epsReal[iw]^2 + epsImag[iw]^2 )
        refIndex[iw]  = sqrt( max(0., (wabs + epsReal[iw]) / 2) )
        extinct[iw]   = sqrt( max(0., (wabs - epsReal[iw]) / 2) )
        absorb[iw]    = 2 * w * extinct[iw] / cLight
    end

    return( (sigmaImag=sigmaImag, epsReal=epsReal, epsImag=epsImag, refIndex=refIndex, extinction=extinct,
             absorption=absorb) )
end


"""
`Plasma.displayOpticalProperties(omegas::Array{Float64,1}, sigmas::Array{Float64,1}, optical::NamedTuple)`
    ... displays the conductivity and the optical constants of the plasma against the photon energy;  nothing is
        returned.
"""
function displayOpticalProperties(omegas::Array{Float64,1}, sigmas::Array{Float64,1}, optical::NamedTuple)
    println("\n  Conductivity and optical constants of the plasma:\n")
    println("  " * "-"^112)
    println("   omega [eV]    omega [Ha]    Re sigma      Im sigma      Re eps        Im eps        n" *
            "             k")
    println("  " * "-"^112)
    for  (iw, omega)  in  enumerate(omegas)
        println("  " * @sprintf("%11.4e", Defaults.convertUnits("energy: from atomic to eV", omega)) * "  " *
                @sprintf("%11.4e", omega) * "  " * @sprintf("%12.5e", sigmas[iw]) * "  " *
                @sprintf("%12.5e", optical.sigmaImag[iw]) * "  " * @sprintf("%12.5e", optical.epsReal[iw]) * "  " *
                @sprintf("%12.5e", optical.epsImag[iw]) * "  " * @sprintf("%12.6f", optical.refIndex[iw]) * "  " *
                @sprintf("%12.6f", optical.extinction[iw]))
    end
    println("  " * "-"^112)
    println("   Re sigma, Im sigma in a.u.;  n and k are the real and imaginary parts of the refractive index.")
    println("   The absorption coefficient alpha = 2 w k/c is available as optical.absorption [1/a_o].\n")

    return( nothing )
end
