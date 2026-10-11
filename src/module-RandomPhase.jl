
"""
`module  JAC.RandomPhase`
    ... a submodel of JAC that orchestrates one `DalgarnoLewis.solve` per ACTIVE electron, which is where
        the many-electron physics this module is named for belongs -- `DalgarnoLewis` itself knows nothing
        about more than one electron and nothing about self-consistency. The split is deliberate: a direct
        inhomogeneous-equation solve (`DalgarnoLewis`) is a numerical primitive, at the same layer as
        `Bsplines` or `RadialIntegrals`; coupling several such solves into a shared response is a different
        kind of code, and the two are kept in separate modules so neither grows unreadable mixing both.

        THE COUPLING TERM -- IMPLEMENTED 11-Oct-2026, DIRECT (RING) TERM ONLY, NOT YET INDEPENDENTLY
        VALIDATED THE WAY DalgarnoLewis's OWN COULOMB-GAUGE SOLVE WAS. `settings.couplingEnabled = true` now
        iterates: for every active electron a, the induced DIRECT (Hartree, "ring") potential at the external
        operator's own multipole rank L, felt from every OTHER active electron b (b excluded from its own
        induced field -- this module builds no exchange term, so there is no self-interaction to subtract
        beyond that), is built from the induced density delta_rho_b(r) = (2j_b+1) * 2 * [Pb(r)*Pchi_b(r) +
        Qb(r)*Qchi_b(r)] -- the (2j_b+1) a STATED CHOICE, treating each active subshell statistically as a
        filled shell, exactly as `Basics.computePotential(::Basics.AaHSField, ...)` already does for its own
        (L=0) Hartree term; the factor 2 is the psi*chi + chi*psi symmetrization for a real response -- via
        `RadialIntegrals.Yk_ab(L, r, Sum_(b!=a) delta_rho_b, mtp, grid)`, projected onto electron a's source
        vector through `DalgarnoLewis.solve`'s new `extraSource` keyword, and the whole process iterated
        (`Settings.maxIterations`, `Settings.accuracy`, measured on the induced potential's own relative
        change between iterations) until it stops moving.

        WHAT HAS BEEN CHECKED, AND WHAT HAS NOT -- READ THIS BEFORE TRUSTING A COUPLED NUMBER. Checked, and
        EXACT: with only ONE active electron, couplingEnabled = true reproduces couplingEnabled = false to
        the last bit (there is no "other" electron to induce anything, so the induced source is identically
        zero) -- the one rigorous, decisive check this first version admits. Checked and NOT reassuring: on a
        toy two-electron-channel case (H's own 1s and 2s treated as two independently active subshells, E1,
        omega = 0.3 a.u., `Basics.NuclearField` orbitals -- see example-At.jl branch c), the induced SOURCE
        starts small (0.3 % of the external one on the first pass) but the DAMPED fixed-point iteration
        (Settings.mixing) still only creeps toward self-consistency, of order 6 % relative change per
        iteration still after 30 -- not divergent, but not practically convergent either within a modest
        budget. The most likely reading, stated rather than hidden: `Basics.NuclearField` orbitals carry NO
        mean-field interaction between 1s and 2s AT ALL, so this toy case asks the coupling term to supply the
        ENTIRE electron-electron interaction by itself, with the statistical (2j_b+1) = 2 weight on each
        subshell on top -- an unusually stiff fixed point, quite possibly not representative of a genuine
        multi-electron SCF system (where the orbitals already carry a mean field and the induced term is a
        CORRECTION to it, not the whole interaction). NOT YET checked against any independent reference
        either way: unlike DalgarnoLewis's Coulomb gauge (validated to 6 figures against a full-eigenbasis
        explicit sum before being trusted anywhere), the (2j_b+1) weight, the factor-of-2 symmetrization and
        the choice of rank L have no literature RPA/RPAE cross-check yet. TREAT A COUPLED NUMBER AS A
        DEMONSTRATION THAT THE MACHINERY RUNS AND REDUCES TO THE RIGHT LIMIT, NOT YET AS A VALIDATED PHYSICAL
        RESULT. The natural next step is a genuine multi-electron SCF system (a real mean-field potential, not
        `NuclearField`) with a known polarizability or RPA correction to check against -- exactly the
        direction of `MultipolePolarizibility`, should that module be unparked.

        WHAT IS STILL OUT OF SCOPE EVEN NOW. Only the DIRECT term is built: the EXCHANGE contribution to the
        induced field needs the same two-electron angular machinery (`SpinAngular`) the CI path already
        carries, a substantially larger undertaking -- "Random Phase Approximation" (no exchange) is the
        right target for this module, not "RPA with Exchange".
"""
module RandomPhase


using  LinearAlgebra, Printf, ..Basics, ..Bsplines, ..DalgarnoLewis, ..ManyElectron, ..Radial, ..RadialIntegrals,
       JenaAtomicCalculator


"""
`struct  RandomPhase.Settings`
    ... defines a type for the details and parameters of a many-electron Dalgarno-Lewis orchestration,
        independent-particle or (direct-term-only) self-consistently coupled.

    + activeSubshells   ::Array{Subshell,1}        ... the electrons whose response is solved for; an empty
                                                         array (the default) means every subshell of the
                                                         `activeOrbitals` dictionary `RandomPhase.solve` is
                                                         called with.
    + targetSubshells   ::Dict{Subshell,Subshell}   ... for each active subshell, the kappa-channel chi is
                                                         solved in -- the same `targetSubshell`
                                                         `DalgarnoLewis.solve` itself requires explicitly;
                                                         this module does not derive it from a selection rule
                                                         any more than `DalgarnoLewis.solve` does.
    + couplingEnabled    ::Bool                      ... True to iterate the self-consistent (direct-term-
                                                         only) induced-potential coupling; False (the default)
                                                         for the independent-particle limit. See the module's
                                                         own docstring for what the coupled path has, and has
                                                         not yet, been checked against.
    + maxIterations      ::Int64                     ... the coupled iteration's budget; unused (every
                                                         electron solved once) while couplingEnabled is false.
    + accuracy           ::Float64                   ... the coupled iteration stops once the induced
                                                         potential's relative change between iterations falls
                                                         below this; unused while couplingEnabled is false.
    + mixing             ::Float64                    ... the LINEAR DAMPING of the induced potential between
                                                         iterations (newInduced = mixing*justComputed +
                                                         (1-mixing)*previous); the plain, undamped iteration
                                                         (mixing = 1.0) was measured NOT to converge on even a
                                                         small two-electron case (11-Oct-2026 -- the induced
                                                         SOURCE itself is small, under 0.3% of the external
                                                         one, so the instability is in the ITERATION, not an
                                                         oversized physical term); unused while couplingEnabled
                                                         is false.
    + dalgarnoLewis       ::DalgarnoLewis.Settings    ... the settings passed through to every individual
                                                         `DalgarnoLewis.solve` call.
    + printBefore         ::Bool                      ... True if a short overview is to be printed before
                                                         the orchestration starts.
"""
struct Settings
    activeSubshells        ::Array{Subshell,1}
    targetSubshells        ::Dict{Subshell,Subshell}
    couplingEnabled         ::Bool
    maxIterations           ::Int64
    accuracy                ::Float64
    mixing                  ::Float64
    dalgarnoLewis            ::DalgarnoLewis.Settings
    printBefore              ::Bool
end


"""
`RandomPhase.Settings()`  ... constructor for an `empty` instance of RandomPhase.Settings.
"""
function Settings()
    Settings( Subshell[], Dict{Subshell,Subshell}(), false, 15, 1.0e-6, 0.3, DalgarnoLewis.Settings(), false )
end


"""
`RandomPhase.Settings(settings::RandomPhase.Settings; activeSubshells=.., targetSubshells=.., couplingEnabled=..,
                       maxIterations=.., accuracy=.., dalgarnoLewis=.., printBefore=..)`
    ... the standard JAC keyword copy-constructor.
"""
function Settings(settings::RandomPhase.Settings;
    activeSubshells::Union{Nothing,Array{Subshell,1}}=nothing,   targetSubshells::Union{Nothing,Dict{Subshell,Subshell}}=nothing,
    couplingEnabled::Union{Nothing,Bool}=nothing,                 maxIterations::Union{Nothing,Int64}=nothing,
    accuracy::Union{Nothing,Float64}=nothing,                     mixing::Union{Nothing,Float64}=nothing,
    dalgarnoLewis::Union{Nothing,DalgarnoLewis.Settings}=nothing, printBefore::Union{Nothing,Bool}=nothing)

    if  isnothing(activeSubshells)   activeSubshellsx = settings.activeSubshells   else   activeSubshellsx = activeSubshells   end
    if  isnothing(targetSubshells)   targetSubshellsx = settings.targetSubshells   else   targetSubshellsx = targetSubshells   end
    if  isnothing(couplingEnabled)   couplingEnabledx = settings.couplingEnabled   else   couplingEnabledx = couplingEnabled   end
    if  isnothing(maxIterations)     maxIterationsx   = settings.maxIterations     else   maxIterationsx   = maxIterations     end
    if  isnothing(accuracy)          accuracyx        = settings.accuracy          else   accuracyx        = accuracy          end
    if  isnothing(mixing)            mixingx          = settings.mixing            else   mixingx          = mixing            end
    if  isnothing(dalgarnoLewis)     dalgarnoLewisx   = settings.dalgarnoLewis     else   dalgarnoLewisx   = dalgarnoLewis      end
    if  isnothing(printBefore)       printBeforex     = settings.printBefore       else   printBeforex     = printBefore       end

    Settings( activeSubshellsx, targetSubshellsx, couplingEnabledx, maxIterationsx, accuracyx, mixingx, dalgarnoLewisx, printBeforex )
end


# `Base.show(io::IO, settings::RandomPhase.Settings)`  ... prepares a proper printout of settings::RandomPhase.Settings.
function Base.show(io::IO, settings::RandomPhase.Settings)
    println(io, "activeSubshells:          $(settings.activeSubshells)  ")
    println(io, "targetSubshells:          $(settings.targetSubshells)  ")
    println(io, "couplingEnabled:          $(settings.couplingEnabled)  ")
    println(io, "maxIterations:            $(settings.maxIterations)  ")
    println(io, "accuracy:                 $(settings.accuracy)  ")
    println(io, "mixing:                   $(settings.mixing)  ")
    println(io, "dalgarnoLewis:            $(settings.dalgarnoLewis)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
end


"""
`struct  RandomPhase.CoupledResponse`
    ... defines a type for the result of a many-electron Dalgarno-Lewis orchestration: one
        `DalgarnoLewis.PerturbedOrbital` per active electron, driven by the same multipole operator at the
        same frequency, independent-particle or coupled depending on the `Settings` it was solved with.

    + multipole          ::EmMultipole                                       ... the multipole rank Mp of
                                                                                  the external operator.
    + gauge               ::EmGauge                                           ... the gauge it was evaluated
                                                                                  in; see `DalgarnoLewis`'s
                                                                                  own docstring for which
                                                                                  gauges are currently usable.
    + omega                ::Float64                                          ... the perturbing frequency,
                                                                                  shared by every electron.
    + perturbedOrbitals    ::Dict{Subshell,DalgarnoLewis.PerturbedOrbital}    ... chi, keyed by the ACTIVE
                                                                                  (source) subshell it was
                                                                                  solved for.
    + nIterations          ::Int64                                            ... always 1 in the
                                                                                  independent-particle limit
                                                                                  (couplingEnabled = false),
                                                                                  which needs no iteration;
                                                                                  the iteration count reached
                                                                                  otherwise.
    + converged             ::Bool                                            ... always true in the
                                                                                  independent-particle limit;
                                                                                  otherwise whether the
                                                                                  induced potential's relative
                                                                                  change fell below
                                                                                  settings.accuracy before
                                                                                  settings.maxIterations.
"""
struct CoupledResponse
    multipole               ::EmMultipole
    gauge                    ::EmGauge
    omega                     ::Float64
    perturbedOrbitals         ::Dict{Subshell,DalgarnoLewis.PerturbedOrbital}
    nIterations               ::Int64
    converged                  ::Bool
end


# `Base.show(io::IO, response::RandomPhase.CoupledResponse)`  ... prepares a proper printout of response::RandomPhase.CoupledResponse.
function Base.show(io::IO, response::RandomPhase.CoupledResponse)
    println(io, "multipole, gauge:         $(response.multipole), $(response.gauge)  ")
    println(io, "omega:                    $(response.omega)  ")
    println(io, "active subshells:         $(collect(keys(response.perturbedOrbitals)))  ")
    println(io, "nIterations, converged:   $(response.nIterations), $(response.converged)  ")
end


"""
`RandomPhase.inducedPotentialZr(multipole::EmMultipole, activeOrbitals::Dict{Subshell,Orbital},
                                 chis::Dict{Subshell,DalgarnoLewis.PerturbedOrbital}, excludeSubshell::Subshell,
                                 grid::Radial.Grid)`
    ... builds the DIRECT (Hartree, "ring") induced potential electron `excludeSubshell` feels from every OTHER
        active electron's CURRENT chi, at the multipole's own rank L, in the SAME "Zr" (= -r*V(r)) convention
        `Radial.Potential` and `Bsplines.computeVlocal` use throughout JAC -- so that the minus sign below is
        the ordinary nuclear-potential-style sign, not a new convention. Each contributing electron b is
        weighted by (2j_b+1), treating it statistically as a filled shell (the same choice
        `Basics.computePotential(::Basics.AaHSField, ...)` makes for its own L=0 Hartree term), and the
        induced density itself is 2*[Pb*Pchi_b + Qb*Qchi_b] -- the psi*chi + chi*psi symmetrization for a real
        response. NOT YET validated against a literature RPA result -- see the module's own docstring.
        A VindZr::Array{Float64,1} of length grid.NoPoints is returned.
"""
function inducedPotentialZr(multipole::EmMultipole, activeOrbitals::Dict{Subshell,Orbital},
                             chis::Dict{Subshell,DalgarnoLewis.PerturbedOrbital}, excludeSubshell::Subshell,
                             grid::Radial.Grid)
    rho = zeros(grid.NoPoints);   mtpRho = 1
    for  (b, chiB)  in  chis
        b == excludeSubshell   &&   continue
        orbB = activeOrbitals[b];   occB = Basics.twice(Basics.subshell_j(b)) + 1
        mtpB = min(length(orbB.P), length(chiB.P), length(orbB.Q), length(chiB.Q))
        for  i = 1:mtpB   rho[i] = rho[i] + occB * 2.0 * (orbB.P[i]*chiB.P[i] + orbB.Q[i]*chiB.Q[i])   end
        mtpRho = max(mtpRho, mtpB)
    end

    return( [ -RadialIntegrals.Yk_ab(multipole.L, grid.r[i], rho, mtpRho, grid)   for i = 1:grid.NoPoints ] )
end


"""
`RandomPhase.inducedSourceVector(VindZr::Array{Float64,1}, sourceOrbital::Orbital, primitives::Bsplines.Primitives)`
    ... projects the induced potential VindZr (in the SAME "Zr" convention as `inducedPotentialZr` returns)
        onto every B-spline primitive against sourceOrbital, exactly reproducing `Bsplines.computeVlocal`'s own
        formula and sign -- the one difference being that computeVlocal's SECOND argument is itself a single
        primitive, where here it is a general tabulated orbital, so the shared formula is written out again
        rather than called. This is `DalgarnoLewis.solve`'s `extraSource` keyword's own input.
        An extraSource::Array{Float64,1} of length primitives.grid.nsL + primitives.grid.nsS is returned.
"""
function inducedSourceVector(VindZr::Array{Float64,1}, sourceOrbital::Orbital, primitives::Bsplines.Primitives)
    grid = primitives.grid;   nsL = grid.nsL;   nsS = grid.nsS
    extraSource = zeros(nsL+nsS)
    mtpP = length(sourceOrbital.P);   mtpQ = length(sourceOrbital.Q)
    for  i = 1:nsL
        bsp = primitives.bsplinesL[i];   lower = bsp.lower;   upper = min(bsp.upper, mtpP)
        upper < lower   &&   continue
        wa = 0.0
        for  r = lower:upper   wa = wa - bsp.bs[r-bsp.lower+1] * VindZr[r] / grid.r[r] * sourceOrbital.P[r] * grid.wr[r]   end
        extraSource[i] = wa
    end
    for  i = 1:nsS
        bsp = primitives.bsplinesS[i];   lower = bsp.lower;   upper = min(bsp.upper, mtpQ)
        upper < lower   &&   continue
        wa = 0.0
        for  r = lower:upper   wa = wa - bsp.bs[r-bsp.lower+1] * VindZr[r] / grid.r[r] * sourceOrbital.Q[r] * grid.wr[r]   end
        extraSource[nsL+i] = wa
    end

    return( extraSource )
end


"""
`RandomPhase.solve(multipole::EmMultipole, gauge::EmGauge, omega::Float64, activeOrbitals::Dict{Subshell,Orbital},
                    potential::Radial.Potential, primitives::Bsplines.Primitives;
                    settings::RandomPhase.Settings=RandomPhase.Settings())`
    ... solves one `DalgarnoLewis` equation per active electron named in `settings.activeSubshells` (or, if
        that is empty, per subshell of `activeOrbitals`), all driven by the SAME external multipole operator
        at the SAME frequency omega. With settings.couplingEnabled == false (the default), this is the
        independent-particle limit: each chi is exactly what a bare `DalgarnoLewis.solve` call for that
        electron alone would return. With settings.couplingEnabled == true, the electrons are iterated
        together, each one's source gaining the DIRECT induced-potential contribution of every OTHER active
        electron's CURRENT chi (`RandomPhase.inducedPotentialZr` + `RandomPhase.inducedSourceVector`), until
        the induced potential's relative change falls below settings.accuracy or settings.maxIterations is
        reached -- see the module's own docstring for what this coupled path has, and has not yet, been
        checked against.
        A response::RandomPhase.CoupledResponse is returned.
"""
function solve(multipole::EmMultipole, gauge::EmGauge, omega::Float64, activeOrbitals::Dict{Subshell,Orbital},
                potential::Radial.Potential, primitives::Bsplines.Primitives;
                settings::RandomPhase.Settings=RandomPhase.Settings())
    subshells = isempty(settings.activeSubshells) ? collect(keys(activeOrbitals)) : settings.activeSubshells
    targets   = Dict{Subshell,Subshell}()
    for  a in subshells
        haskey(activeOrbitals, a)  ||
            error("RandomPhase.solve(): settings.activeSubshells names $a, which activeOrbitals does not carry.")
        target = get(settings.targetSubshells, a, nothing)
        isnothing(target)  &&
            error("RandomPhase.solve(): no target subshell given for active subshell $a -- " *
                  "settings.targetSubshells must name one, exactly as DalgarnoLewis.solve itself requires.")
        targets[a] = target
    end

    if  settings.printBefore
        @printf(">> RandomPhase.solve():  %d active electron(s), multipole = %s, gauge = %s, omega = %.8f, coupled = %s\n",
                length(subshells), string(multipole), string(gauge), omega, string(settings.couplingEnabled))
    end

    # THE INDEPENDENT-PARTICLE LIMIT -- one DalgarnoLewis.solve per electron, no coupling term, no iteration.
    if  !settings.couplingEnabled
        perturbedOrbitals = Dict{Subshell,DalgarnoLewis.PerturbedOrbital}()
        for  a in subshells
            perturbedOrbitals[a] = DalgarnoLewis.solve(multipole, gauge, omega, targets[a], activeOrbitals[a],
                                                         potential, primitives; settings=settings.dalgarnoLewis)
        end

        return( CoupledResponse(multipole, gauge, omega, perturbedOrbitals, 1, true) )
    end

    # THE COUPLED ITERATION. Start from the independent-particle chi's (iteration 0, no induced potential yet),
    # then repeatedly rebuild each electron's induced source from the OTHERS' current chi and re-solve.
    grid = primitives.grid
    chis = Dict{Subshell,DalgarnoLewis.PerturbedOrbital}()
    for  a in subshells
        chis[a] = DalgarnoLewis.solve(multipole, gauge, omega, targets[a], activeOrbitals[a], potential, primitives;
                                       settings=settings.dalgarnoLewis)
    end

    # THE DAMPED (LINEARLY MIXED) ITERATION -- see Settings.mixing's own docstring for why the plain,
    # undamped iteration (mixing = 1.0) does not converge.
    converged = false;   nIterations = 0;   previousExtraSource = Dict{Subshell,Array{Float64,1}}()
    for  iteration = 1:settings.maxIterations
        nIterations = iteration
        newChis = Dict{Subshell,DalgarnoLewis.PerturbedOrbital}();   maxRelChange = 0.0
        for  a in subshells
            VindZr      = RandomPhase.inducedPotentialZr(multipole, activeOrbitals, chis, a, grid)
            justComputed = RandomPhase.inducedSourceVector(VindZr, activeOrbitals[a], primitives)
            priorSource  = get(previousExtraSource, a, zeros(length(justComputed)))
            dampedSource = settings.mixing * justComputed + (1.0 - settings.mixing) * priorSource

            newChis[a] = DalgarnoLewis.solve(multipole, gauge, omega, targets[a], activeOrbitals[a], potential,
                                              primitives; settings=settings.dalgarnoLewis, extraSource=dampedSource)

            denom = LinearAlgebra.norm(priorSource)
            # denom == 0 on the FIRST iteration just means "nothing to compare against yet" (Inf, keep
            # iterating) -- UNLESS dampedSource is ALSO zero (e.g. a single active electron, which has no
            # coupling partner at all and so no induced source ever), which already IS the fixed point.
            relChange = denom > 0.0 ? LinearAlgebra.norm(dampedSource - priorSource) / denom :
                        (LinearAlgebra.norm(dampedSource) > 0.0 ? Inf : 0.0)
            maxRelChange = max(maxRelChange, relChange)
            previousExtraSource[a] = dampedSource
        end
        chis = newChis

        if  settings.printBefore
            @printf(">> RandomPhase.solve():  iteration %d, max relative change of the (damped) induced source = %.3e\n",
                    iteration, maxRelChange)
        end
        if  maxRelChange < settings.accuracy
            converged = true
            break
        end
    end

    return( CoupledResponse(multipole, gauge, omega, chis, nIterations, converged) )
end


end # module
