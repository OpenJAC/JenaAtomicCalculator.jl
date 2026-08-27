
"""
`module  JAC.ReducedDensityMatrix`
    ... a submodel of JAC that contains all methods for computing reduced density matrices, natural orbitals, radial density distributions
        and related information for some level(s).
"""
module ReducedDensityMatrix


using Printf, ..Basics, ..Defaults, ..InteractionStrength, ..ManyElectron, ..Nuclear, ..Radial, ..RadialIntegrals,
                ..SpinAngular, ..TableStrings


"""
`struct  ReducedDensityMatrix.Settings  <:  AbstractPropertySettings`
    ... defines a type for the parameters that control the computation of reduced density matrices and the quantities derived from them.

    + calcNatural              ::Bool             ... True if natural orbitals need to be calculated, and false otherwise.
    + calcDensity              ::Bool             ... True if the radial density need to be calculated, and false otherwise.
    + calcIpq                  ::Bool             ... True if orbital interaction I_pq need to be calculated, and false otherwise.
    + calc2pRDM                ::Bool             ... True if the 2p RDM need to be calculated, and false otherwise;
                                                        the 1pRDM is calculated by default.
    + printBefore              ::Bool             ... True if a list of selected levels is printed before the actual computations start. 
    + levelSelection           ::LevelSelection   ... Specifies the selected levels, if any.
"""
struct Settings  <:  AbstractPropertySettings 
    calcNatural                ::Bool 
    calcDensity                ::Bool 
    calcIpq                    ::Bool 
    calc2pRDM                  ::Bool 
    printBefore                ::Bool 
    levelSelection             ::LevelSelection
end 


"""
`ReducedDensityMatrix.Settings(; calcNatural::Bool=true, calcDensity::Bool=true, calcIpq::Bool=false,
                                calc2pRDM::Bool=false, printBefore::Bool=true, levelSelection::LevelSelection=LevelSelection())`
    ... keyword constructor to overwrite selected values of the reduced-density-matrix computations.
        A settings::ReducedDensityMatrix.Settings is returned.
"""
function Settings(; calcNatural::Bool=true, calcDensity::Bool=true, calcIpq::Bool=false,
                    calc2pRDM::Bool=false, printBefore::Bool=true, levelSelection::LevelSelection=LevelSelection())
    Settings(calcNatural, calcDensity, calcIpq, calc2pRDM, printBefore, levelSelection)
end


"""
`Base.show(io::IO, settings::ReducedDensityMatrix.Settings)`
    ... prepares a proper printout of the variable settings::ReducedDensityMatrix.Settings; nothing is returned.
"""
function Base.show(io::IO, settings::ReducedDensityMatrix.Settings) 
    println(io, "calcNatural:              $(settings.calcNatural)  ")
    println(io, "calcDensity:              $(settings.calcDensity)  ")
    println(io, "calcIpq:                  $(settings.calcIpq)  ")
    println(io, "calc2pRDM:                $(settings.calc2pRDM)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "levelSelection:           $(settings.levelSelection)  ")
end


"""
`struct  ReducedDensityMatrix.Outcome`  
    ... defines a type to keep the outcome of a reduced density-matrix and natural-orbital computation.

    + level                     ::Level                 ... Atomic level to which the outcome refers to.
    + naturalSubshells          ::Array{Subshell,1}     ... List of natural orbitals (subshells).
    + naturalOccupation         ::Array{Float64,1}      ... Occupation numbers of natural orbitals.
    + naturalOrbitalExpansion   ::Dict{Subshell, Array{Float64,1}} 
        ... Dictionary of the expansion of (one-electron) natural orbitals in terms of the standard orbitals
            from the given basis.
    + naturalOrbitals           ::Dict{Subshell, Orbital} ... Dictionary of (one-electron) natural orbitals.
    + orbitalInteraction        ::Array{Float64,2}      ... Orbital interaction I_pq = I[ip,iq].
    + rho1p                     ::Array{Float64,2}      ... One-particle RDM rho^(1p) [ip,iq].
    + rho2p                     ::Array{Float64,4}      ... Two-particle RDM rho^(2p) [ip,iq,ir,is].
    + electronDensity           ::Radial.Density        ... Radial density distribution.
"""
struct Outcome 
    level                       ::Level 
    naturalSubshells            ::Array{Subshell,1}
    naturalOccupation           ::Array{Float64,1}
    naturalOrbitalExpansion     ::Dict{Subshell, Array{Float64,1}}
    naturalOrbitals             ::Dict{Subshell, Orbital}
    orbitalInteraction          ::Array{Float64,2}
    rho1p                       ::Array{Float64,2}
    rho2p                       ::Array{Float64,4}
    electronDensity             ::Radial.Density
end 


"""
`ReducedDensityMatrix.Outcome()`
    ... constructor for an `empty` instance of ReducedDensityMatrix.Outcome, i.e. with the level and every derived quantity unset.
        An outcome::ReducedDensityMatrix.Outcome is returned.
"""
function Outcome()
    Outcome(Level(), Subshell[], Float64[], Dict{Subshell, Array{Float64,1}}, Dict{Subshell, Orbital}(), 
            zeros(1,1), zeros(1,1), zeros(1,1,1,1), Radial.Density() )
end


"""
`Base.show(io::IO, outcome::ReducedDensityMatrix.Outcome)`
    ... prepares a proper printout of the variable outcome::ReducedDensityMatrix.Outcome; nothing is returned.
"""
function Base.show(io::IO, outcome::ReducedDensityMatrix.Outcome) 
    println(io, "level:                   $(outcome.level)  ")
    println(io, "naturalSubshells:        $(outcome.naturalSubshells)  ")
    println(io, "naturalOccupation:       $(outcome.naturalOccupation)  ")
    println(io, "naturalOrbitalExpansion: $(outcome.naturalOrbitalExpansion)  ")
    println(io, "naturalOrbitals:         $(outcome.naturalOrbitals)  ")
    println(io, "orbitalInteraction:      $(outcome.orbitalInteraction)  ")
    println(io, "rho1p:                   $(outcome.rho1p)  ")
    println(io, "rho2p:                   $(outcome.rho2p)  ")
    println(io, "electronDensity:         $(outcome.electronDensity)  ")
end


"""
`ReducedDensityMatrix.compute1pRDM(level::Level)`  
    ... to compute 1p RDM; a rdm::Array{Float64,2} is returned.
"""
function  compute1pRDM(level::Level)
    subshellList = level.basis.subshells
    lenNO = length(subshellList);    rho_pq     = zeros(lenNO,lenNO)
    # Cycle over all matrix elements of the CSF basis
    opa          = SpinAngular.OneParticleOperator(0, plus, true)
    for  (ir, rcsf) in enumerate(level.basis.csfs)
        for  (is, scsf) in enumerate(level.basis.csfs)
            # Calculate angular coefficient for rank-0 operator
            wa = SpinAngular.computeCoefficients(opa, rcsf, scsf, subshellList) 
            # Cycle over the pair of natural subshells in rho_pq
            for (ip,p)  in  enumerate(subshellList)
                for (iq,q)  in  enumerate(subshellList)
                    for  coeff in wa
                        if  (p == coeff.a   &&  q == coeff.b)  ||  (p == coeff.b   &&  q == coeff.a)
                            jj = Basics.subshell_2j(level.basis.orbitals[coeff.a].subshell)
                            rho_pq[ip,iq] = rho_pq[iq,ip] = rho_pq[ip,iq] + level.mc[ir] * coeff.T * sqrt( jj + 1) * level.mc[is]
                        end
                    end
                end
            end
        end
    end

    return( rho_pq )
end


"""
`ReducedDensityMatrix.compute1pRDMDirect(level::Level)`
    ... to compute the 1p RDM directly from the CI mixing coefficients and CSF occupation numbers, following Eq. (7) of Ma et al., Atoms 12,
        30 (2024), rho^kappa_nn' = sum_ij c_i v^ij_nn' c_j -- but derived here from elementary second quantization, independently of
        SpinAngular.computeCoefficients.

        Diagonal elements need no angular-momentum recoupling at all: the expectation value of a number operator in a normalized state is
        just the mc-weighted sum of that CSF's own occupation number.

        Off-diagonal elements are supported for the simplest, and here relevant, case only: two CSFs that differ by moving ONE electron
        between two subshells of the SAME kappa, each holding 0 or 1 electron (e.g. a single 2s->3s substitution), with all other subshells
        in identical coupling. For this case the transfer coefficient is the universal value sqrt(2), independent of kappa/j: writing the
        singly-occupied pair coupled to a scalar (J=0) as |ab;0> = (1/sqrt(2j+1)) sum_m (-1)^(j-m) a+_{am} a+_{b,-m} |0> (the standard
        pair-coupling identity) and the closed-shell pair as |aa;0> = (1/sqrt(2*(2j+1))) sum_m (-1)^(j-m) a+_{am} a+_{a,-m} |0> (note the
        extra factor of 2 under the root here, from the (m,-m) double counting that only affects the same-shell case), direct evaluation of
        N_ab = sum_m a+_{am} a_{bm} gives N_ab |ab;0> = sqrt(2) |aa;0> for ANY j. CSF pairs with a more general occupation pattern (more
        than one electron differing per subshell) are NOT handled and raise an error -- those need the full
        coefficient-of-fractional-parentage machinery in SpinAngular.
        A rho_pq::Array{Float64,2} is returned.
"""
function  compute1pRDMDirect(level::Level)
    subshellList = level.basis.subshells;    lenNO = length(subshellList)
    rho_pq       = zeros(lenNO, lenNO)
    csfs         = level.basis.csfs

    # Diagonal expectation values: mc-weighted CSF occupation numbers, no coupling needed.
    for  (i, csf)  in  enumerate(csfs)
        for  ip = 1:lenNO    rho_pq[ip,ip] = rho_pq[ip,ip] + level.mc[i]^2 * csf.occupation[ip]    end
    end

    # Off-diagonal (cross-CSF) contributions from single-electron substitutions between subshells of the same kappa.
    for  (ir, rcsf)  in  enumerate(csfs)
        for  (is, scsf)  in  enumerate(csfs)
            if  ir == is    continue    end
            diffSubshells = Int64[];    validSingleSubstitution = true
            for  i = 1:lenNO
                diff = rcsf.occupation[i] - scsf.occupation[i]
                if       diff == 0          continue
                elseif   abs(diff) == 1     push!(diffSubshells, i)
                else     validSingleSubstitution = false;    break
                end
            end
            if  !validSingleSubstitution  ||  length(diffSubshells) != 2    continue    end

            ia, ib = diffSubshells
            if  subshellList[ia].kappa != subshellList[ib].kappa    continue    end
            # Determine, order-independently, whether each of the two differing subshells oscillates between
            # {0,1} (empty-type) or {full-1,full} (full-type) across the CSF pair -- direction (which CSF is
            # bra/ket) does not matter since the sqrt(2) transfer coefficient is accumulated symmetrically below.
            fullA = Basics.subshell_2j(subshellList[ia]) + 1;    fullB = Basics.subshell_2j(subshellList[ib]) + 1
            occAset = Set([rcsf.occupation[ia], scsf.occupation[ia]])
            occBset = Set([rcsf.occupation[ib], scsf.occupation[ib]])
            if      occAset == Set([0,1])        &&  occBset == Set([fullB-1,fullB])
                emptySite, fullSite = ia, ib
            elseif  occAset == Set([fullA-1,fullA])  &&  occBset == Set([0,1])
                emptySite, fullSite = ib, ia
            else
                error("ReducedDensityMatrix.compute1pRDMDirect: only single substitutions between a fully-closed " *
                        "shell (one hole created) and an empty shell (one electron added) are supported; got "     *
                        "occupations $(rcsf.occupation[ia]),$(scsf.occupation[ia]) and $(rcsf.occupation[ib]),"    *
                        "$(scsf.occupation[ib]).")
            end
            # All other subshells must be in identical coupling for the plain sqrt(2) transfer coefficient to apply.
            samecoupling = true
            for  i = 1:lenNO
                if  i == ia  ||  i == ib    continue    end
                if  rcsf.occupation[i]  != scsf.occupation[i]   ||  rcsf.subshellJ[i] != scsf.subshellJ[i]  ||
                    rcsf.subshellX[i]   != scsf.subshellX[i]     ||  rcsf.seniorityNr[i] != scsf.seniorityNr[i]
                    samecoupling = false;    break
                end
            end
            if  !samecoupling    continue    end

            # The operator a+_fullSite a_emptySite has a NONZERO matrix element in only ONE (ir,is) direction:
            # bra (rcsf) = the closed configuration (full at fullSite, empty at emptySite), ket (scsf) = the open
            # one. The opposite ordering is an exact zero (it would need to annihilate an electron at emptySite
            # that the closed-shell ket does not have) and must be skipped here -- accumulating it too would
            # double the physical value. Hermiticity/reality already mirrors this single contribution into
            # rho_pq[emptySite,fullSite] below.
            fullCapacity = Basics.subshell_2j(subshellList[fullSite]) + 1
            if  rcsf.occupation[fullSite] == fullCapacity  &&  rcsf.occupation[emptySite] == 0
                rho_pq[fullSite,emptySite] = rho_pq[emptySite,fullSite] =
                    rho_pq[fullSite,emptySite] + level.mc[ir] * sqrt(2.) * level.mc[is]
            end
        end
    end

    return( rho_pq )
end


"""
`ReducedDensityMatrix.compute2pRDM(level::Level)`
    ... to compute 2p RDM for all pairs (p,q;r,s) of subshells; a rdm::Array{Float64,4} is returned.
"""
function  compute2pRDM(level::Level)
    subshellList = level.basis.subshells
    lenNO = length(subshellList);    rdm     = zeros(lenNO,lenNO,lenNO,lenNO)
    # Cycle over all matrix elements of the CSF basis
    opa = SpinAngular.TwoParticleOperator(0, plus, true)
    for  (i, icsf) in enumerate(level.basis.csfs)
        for  (j, jcsf) in enumerate(level.basis.csfs)
            # Calculate angular coefficient for rank-0 operator
            wa = SpinAngular.computeCoefficients(opa, icsf, jcsf, subshellList) 
            # Cycle over the pair of natural subshells in rho_pqrs
            for (ip,p)  in  enumerate(subshellList)
                for (iq,q)  in  enumerate(subshellList)
                    for (ir,r)  in  enumerate(subshellList)
                        for (is,s)  in  enumerate(subshellList)
                            for  coeff in wa
                                if  (p == coeff.a   &&  q == coeff.b   &&  r == coeff.c   &&  s == coeff.d)  ||  
                                    (p == coeff.c   &&  q == coeff.d   &&  r == coeff.a   &&  s == coeff.b)
                                    rdm[ip,iq,ir,is] = rdm[ir,is,ip,iq] = rdm[ip,iq,ir,is] + level.mc[i] * coeff.V * level.mc[j]
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    return( rdm )
end


"""
`ReducedDensityMatrix.computeNaturalOrbitalExpansion(rho1p::Array{Float64,2}, level::Level)`
    ... to perform the expansion of the natural orbitals for the given level in terms of the standard orbitals as defined by the basis. A
        tuple (naturalOcc::Array{Float64,1}, naturalExp::Dict{Subshell, Array{Float64,1}}) is returned which provides the natural occpuation
        numbers and expansion coefficients with regard to the given list of natural orbitals. This method makes use of the 1-particle RDM to
        extract and diagonalize the expansion matrix. Since the (rank-0) 1p RDM has no matrix elements between subshells of different
        relativistic symmetry kappa, it is block-diagonal in kappa; each kappa-block is diagonalized separately (Eqs. 7-8 of Ma et al.,
        Atoms 12, 30 (2024)) and the resulting natural orbitals (sorted by decreasing occupation within the block) are associated, in turn,
        with the standard subshells of that same kappa.
"""
function  computeNaturalOrbitalExpansion(rho1p::Array{Float64,2}, level::Level)
    subshellList = level.basis.subshells;    lenNO = length(subshellList)
    naturalOcc   = zeros(Float64, lenNO);     naturalExp = Dict{Subshell, Array{Float64,1}}()

    kappaGroups = Dict{Int64, Array{Int64,1}}()
    for  (i, sh)  in  enumerate(subshellList)     push!( get!(kappaGroups, sh.kappa, Int64[]), i )     end

    for  (_, idxList)  in  kappaGroups
        subMatrix = rho1p[idxList, idxList]
        eigen     = Basics.fixEigenvectorPhase!( Basics.diagonalize(MatrixWithLinearAlgebra(), subMatrix) )
        order     = sortperm(eigen.values, rev=true)                     # decreasing natural occupation
        for  (k, ord)  in  enumerate(order)
            i             = idxList[k]
            labelSubshell = subshellList[i]
            naturalOcc[i] = eigen.values[ord]
            fullExp       = zeros(Float64, lenNO)
            fullExp[idxList] = eigen.vectors[ord]
            naturalExp[labelSubshell] = fullExp
        end
    end

    return( naturalOcc, naturalExp )
end


"""
`ReducedDensityMatrix.computeNaturalOrbitals(naturalSubshells::Array{Subshell,1}, naturalExp::Dict{Subshell, Array{Float64,1}},
                                            level::Level)`
    ... to compute the natural orbitals as superposition of the standard orbitals; for each subshell in naturalExp, an orbital is computed,
        and naturalOrbitals::Dict{Subshell, Orbital} returned. Follows Eqs. 10-11 of Ma et al., Atoms 12, 30 (2024): P-tilde_n'kappa(r) =
        sum_n u^kappa_{n,n'} P_nkappa(r), and likewise for Q.
"""
function  computeNaturalOrbitals(naturalSubshells::Array{Subshell,1}, naturalExp::Dict{Subshell, Array{Float64,1}}, level::Level)
    subshellList    = level.basis.subshells
    naturalOrbitals = Dict{Subshell, Orbital}()

    for  subsh  in  naturalSubshells
        exp = naturalExp[subsh]
        mtp = 0
        for  (i, s)  in  enumerate(subshellList)
            if  exp[i] != 0.    mtp = max(mtp, length(level.basis.orbitals[s].P))    end
        end
        P = zeros(mtp);   Q = zeros(mtp);   Pprime = zeros(mtp);   Qprime = zeros(mtp);   energy = 0.

        for  (i, s)  in  enumerate(subshellList)
            c = exp[i];    if  c == 0.    continue    end
            orb = level.basis.orbitals[s];   n = length(orb.P)
            P[1:n]      = P[1:n]      + c * orb.P
            Q[1:n]      = Q[1:n]      + c * orb.Q
            if  !isempty(orb.Pprime)
                Pprime[1:n] = Pprime[1:n] + c * orb.Pprime;   Qprime[1:n] = Qprime[1:n] + c * orb.Qprime
            end
            energy = energy + c^2 * orb.energy
        end

        parentOrb = level.basis.orbitals[subsh]
        naturalOrbitals[subsh] = Orbital(subsh, parentOrb.isBound, parentOrb.useStandardGrid, energy, P, Q, Pprime, Qprime, parentOrb.grid)
    end

    return( naturalOrbitals )
end


"""
`ReducedDensityMatrix.computeOrbitalInteractions(naturalSubshells::Array{Subshell,1}, naturalOrbitals::Dict{Subshell, Orbital},
                                                level::Level)`
    ... to compute the orbital interaction I_pq between pairs of natural orbitals; a matrix orbitalInteraction::Array{Float64,2} is
        returned. Not yet implemented -- this quantity is not part of the Ma et al. (2024) recipe used elsewhere in this module and still
        needs its own definition/derivation.
"""
function  computeOrbitalInteractions(naturalSubshells::Array{Subshell,1}, naturalOrbitals::Dict{Subshell, Orbital}, level::Level)
    ns = length(naturalSubshells)
    orbitalInteraction = zeros(ns,ns)    
    println(">> computeOrbitalInteractions() ... not yet implemented !!")

    return( orbitalInteraction )
end


"""
`ReducedDensityMatrix.computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                        settings::ReducedDensityMatrix.Settings; output=true)` 
    ... to compute (as selected) the 1-particle and 2-particle reduced density matrices, natural orbitals or other requested information for
        the levels of interest. This is the driver of the module: it determines the outcomes, computes the properties of each and displays
        the results to screen and to the summary file.
        A newOutcomes::Array{ReducedDensityMatrix.Outcome,1} is returned for output = true, and nothing otherwise.
"""
function computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                            settings::ReducedDensityMatrix.Settings; output=true)
    println("")
    printstyled("ReducedDensityMatrix.computeOutcomes(): The computation of the reduced density matrices starts now ... \n", color=:light_green)
    printstyled("------------------------------------------------------------------------------------------------------ \n", color=:light_green)

    outcomes = ReducedDensityMatrix.determineOutcomes(multiplet, settings)
    # Display all selected levels before the computations start
    if  settings.printBefore    ReducedDensityMatrix.displayOutcomes(outcomes, settings)    end
    # Calculate all requested properties for each outcome
    newOutcomes = ReducedDensityMatrix.Outcome[]
    for  outcome in outcomes
        newOutcome = ReducedDensityMatrix.computeProperties(outcome, nm, grid, settings) 
        push!( newOutcomes, newOutcome)
    end
    # Print all results to screen
    ReducedDensityMatrix.displayResults(stdout, newOutcomes, grid, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    ReducedDensityMatrix.displayResults(iostream, newOutcomes, grid, settings)   end

    if    output    return( newOutcomes )
    else            return( nothing )
    end
end


"""
`ReducedDensityMatrix.computeProperties(outcome::ReducedDensityMatrix.Outcome, nm::Nuclear.Model, grid::Radial.Grid,
                                        settings::ReducedDensityMatrix.Settings)`
    ... to compute all properties for a given level; an outcome::ReducedDensityMatrix.Outcome is returned for which the properties are now
        evaluated explicitly.
"""
function  computeProperties(outcome::ReducedDensityMatrix.Outcome, nm::Nuclear.Model, grid::Radial.Grid,
                            settings::ReducedDensityMatrix.Settings)
    naturalOccupation = outcome.naturalOccupation;        naturalOrbitalExp  = outcome.naturalOrbitalExpansion;  
    naturalOrbitals   = outcome.naturalOrbitals;          orbitalInteraction = outcome.orbitalInteraction;       
    electronDensity   = outcome.electronDensity;          rho1p              = outcome.rho1p;             rho2p = outcome.rho2p
    
    naturalSubshells   = copy(outcome.level.basis.subshells)
    rho1p              = compute1pRDMDirect(outcome.level)
    naturalOccupation, naturalOrbitalExp = computeNaturalOrbitalExpansion(rho1p, outcome.level)

    if  settings.calcNatural  ||  settings.calcDensity  ||  settings.calcIpq
        naturalOrbitals = computeNaturalOrbitals(naturalSubshells, naturalOrbitalExp, outcome.level)
    end

    if  settings.calcDensity
        electronDensity = computeRadialDistribution(naturalOccupation, naturalOrbitals, naturalSubshells, grid)
    end

    if  settings.calcIpq
        orbitalInteraction = computeOrbitalInteractions(naturalSubshells, naturalOrbitals, outcome.level)
    end

    if  settings.calc2pRDM
        rho2p = compute2pRDM(outcome.level)
    end

    newOutcome = ReducedDensityMatrix.Outcome( outcome.level, naturalSubshells,  naturalOccupation, naturalOrbitalExp,
                                                naturalOrbitals, orbitalInteraction, rho1p, rho2p, electronDensity )

    return( newOutcome )
end


"""
`ReducedDensityMatrix.computeRadialDistribution(naturalOcc::Array{Float64,1}, naturalOrbitals::Dict{Subshell, Orbital},
                                                    naturalSubshells::Array{Subshell,1}, grid::Radial.Grid)`
    ... to compute the total radial electron density from the (occupation-weighted) natural orbitals, D(r) = sum_(n'kappa) occ_(n'kappa) [
        P-tilde_(n'kappa)(r)^2 + Q-tilde_(n'kappa)(r)^2 ]. A Radial.Density is returned w.r.t. the given grid.
"""
function  computeRadialDistribution(naturalOcc::Array{Float64,1}, naturalOrbitals::Dict{Subshell, Orbital},
                                        naturalSubshells::Array{Subshell,1}, grid::Radial.Grid)
    Dr = zeros(grid.NoPoints)
    for  (i, subsh)  in  enumerate(naturalSubshells)
        occ = naturalOcc[i];    if  occ == 0.    continue    end
        orb = naturalOrbitals[subsh];    mtp = length(orb.P)
        for  k = 1:mtp    Dr[k] = Dr[k] + occ * (orb.P[k]^2 + orb.Q[k]^2)    end
    end

    return( Radial.Density("Radial density from natural orbitals", Dr, grid) )
end


"""
`ReducedDensityMatrix.determineOutcomes(multiplet::Multiplet, settings::ReducedDensityMatrix.Settings)`  
    ... to determine a list of Outcomes's for the computation of reduced density matrices and other information for the given multiplet. It
        takes into account the particular selections and settings.
        An Array{ReducedDensityMatrix.Outcome,1} is returned. Apart from the level specification, all physical 
        properties are set to zero during the initialization process.
"""
function  determineOutcomes(multiplet::Multiplet, settings::ReducedDensityMatrix.Settings) 
    outcomes = ReducedDensityMatrix.Outcome[]
    for  level  in  multiplet.levels
        if  Basics.selectLevel(level, settings.levelSelection)
            push!( outcomes, ReducedDensityMatrix.Outcome(level, Subshell[], Float64[], Dict{Subshell, Array{Float64,1}}(), 
                                    Dict{Subshell, Orbital}(), zeros(1,1), zeros(1,1), zeros(1,1,1,1), Radial.Density()) )
        end
    end
    return( outcomes )
end


"""
`ReducedDensityMatrix.displayOutcomes(outcomes::Array{ReducedDensityMatrix.Outcome,1}, settings::ReducedDensityMatrix.Settings)`  
    ... to display a list of levels that have been selected for the computations. A small neat table of all selected levels and their
        energies is printed but nothing is returned otherwise. Moreover, the selected properties and information is printed as well, though
        not yet calculated.
"""
function  displayOutcomes(outcomes::Array{ReducedDensityMatrix.Outcome,1}, settings::ReducedDensityMatrix.Settings)
    nx = 43
    println(" ")
    println("  Results for the natural occupation number, natural orbitals and kp RDM are printed in turn for the following levels:")
    println(" ")
    println("  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "Energy"; na=4);              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 

    for  outcome in outcomes
        sa  = "  ";    sym = Basics.LevelSymmetry( outcome.level.J, outcome.level.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(outcome.level.index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", outcome.level.energy)) * "    "
        println( sa )
    end
    println("  ", TableStrings.hLine(nx))

    println("\n\n  The following properties are calculated: \n")
    println("  + Natural occupations numbers")
    println("  + Natural orbital expansion in terms of the standard orbitals.")
    if  settings.calcNatural   println("  + Representation of the natural orbitals.")  end
    if  settings.calcDensity   println("  + Radial electron density distribution.")    end
    if  settings.calcIpq       println("  + Orbital interaction I_pq.")                end
    
    return( nothing )
end


"""
`ReducedDensityMatrix.displayResults(stream::IO, outcomes::Array{ReducedDensityMatrix.Outcome,1}, grid::Radial.Grid,
                                     settings::ReducedDensityMatrix.Settings)`
    ... to display the computed results for all outcomes on the given stream. The selected levels and their energies are tabulated first,
        and then, for each level in turn, the 1-particle RDM, the natural occupation numbers and the natural orbital expansion; the natural
        orbitals, the radial electron density, the orbital interactions I_pq and the 2-particle RDM follow for whichever of them the
        settings requested. Nothing is returned.
"""
function  displayResults(stream::IO, outcomes::Array{ReducedDensityMatrix.Outcome,1}, grid::Radial.Grid, 
                            settings::ReducedDensityMatrix.Settings)
    nx = 43
    println(stream, " ")
    println(stream, "  Results for the natural occupation number, natural orbitals and kp RDM are printed in turn for the following levels:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "Energy"; na=4);              
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 

    for  outcome in outcomes
        sa  = "  ";    sym = Basics.LevelSymmetry( outcome.level.J, outcome.level.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(outcome.level.index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", outcome.level.energy)) * "    "
        println(stream,  sa )
    end
    println(stream, "  ", TableStrings.hLine(nx))

    # Now print all selected result for each outcome
    for  outcome in outcomes
        nx = 120;       sym = Basics.LevelSymmetry( outcome.level.J, outcome.level.parity) 
        println(stream, " ")
        println(stream, " ")
        println(stream, "  =============================")
        println(stream, "  Level: $(outcome.level.index) with symmetry $sym  ")
        println(stream, "  =============================")

        println(stream, " ")
        println(stream, "  1-particle RDM:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        wa = TableStrings.subshellList(10, outcome.level.basis.subshells)
        for  (i, sa)  in  enumerate(wa)
            if  i == 1      sb = "     Orb | Orb  " * sa
            else            sb = "                " * sa      end
            println(stream,  sb)
        end
        println(stream, "  ", TableStrings.hLine(nx))

        for  (i, subsh)  in  enumerate(outcome.level.basis.subshells)
            sb = string(subsh)
            wa = TableStrings.floatList(10, outcome.rho1p[i,:])
            for  (j, sa)  in  enumerate(wa)
                if       j == 1  sb = "     " * sb * "     " * sa
                else             sb = "                    " * sa      end
                println(stream,  sb)
            end
        end
        println(stream, "  ", TableStrings.hLine(nx))

        println(stream, " ")
        println(stream, "  Natural occupation numbers:  ")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        wa = TableStrings.subshellList(10, outcome.naturalSubshells)
        for  (i, sa)  in  enumerate(wa)
            if       i == 1  sb = "     NO   " * sa
            else             sb = "          " * sa      end
            println(stream,  sb)
        end
        println(stream, "  ", TableStrings.hLine(nx))

        wa = TableStrings.floatList(10, outcome.naturalOccupation)
        for  (i, sa)  in  enumerate(wa)
            println(stream,  "          " * sa)
        end
        println(stream, "  ", TableStrings.hLine(nx))

        nx = 120
        println(stream, " ")
        println(stream, "  Natural orbital expansion:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        wa = TableStrings.subshellList(10, outcome.naturalSubshells)
        for  (i, sa)  in  enumerate(wa)
            if       i == 1  sb = "     Orbitals: " * sa
            else             sb = "               " * sa      end
            println(stream,  sb)
        end
        println(stream, "  ", TableStrings.hLine(nx))

        for  subsh  in  outcome.naturalSubshells
            sb = string(subsh)
            wa = TableStrings.floatList(10, outcome.naturalOrbitalExpansion[subsh])
            for  (j, sa)  in  enumerate(wa)
                if       j == 1  sb = "     " * sb * "   " * sa
                else             sb = "                  " * sa      end
                println(stream,  sb)
            end
        end
        println(stream, "  ", TableStrings.hLine(nx))

        if  settings.calcNatural
            println(stream, "\n  Natural orbitals are calculated for the following subshells and are available by outcome.naturalOrbitals:\n")
            wa = TableStrings.subshellList(10, outcome.naturalSubshells)
            nx = 120
            println(stream, "  ", TableStrings.hLine(nx))
            for  (i, sa)  in  enumerate(wa)
                if       i == 1  sb = "     NO:  " * sa
                else             sb = "          " * sa      end
                println(stream,  sb)
            end
            println(stream, "  ", TableStrings.hLine(nx))
        end

        if  settings.calcDensity
            println(stream, "\n  Radial electron distribution:\n")
            nx = 42;    sa = "         i)         r[i]       density[i]"
            println(stream, "  ", TableStrings.hLine(nx))
            println(stream,  sa)
            println(stream, "  ", TableStrings.hLine(nx))
            for  i = 1:30:outcome.electronDensity.grid.NoPoints
                println(stream, "   " * ("      "*string(i))[end-6:end] * ")      " *
                                @sprintf("%.4e", grid.r[i]) * "    " * @sprintf("%.4e", outcome.electronDensity.Dr[i])  )
            end
            println(stream, "  ", TableStrings.hLine(nx))
        end

        if  settings.calcIpq
            nx = 120
            println(stream, "\n  Orbital interactions I_pq:\n")
            println(stream, "  ", TableStrings.hLine(nx))
            wa = TableStrings.subshellList(10, outcome.naturalSubshells)
            for  (i, sa)  in  enumerate(wa)
                if  i == 1      sb = "     Orb | Orb  " * sa
                else            sb = "                " * sa      end
                println(stream,  sb)
            end
            println(stream, "  ", TableStrings.hLine(nx))

            for  (i, subsh)  in  enumerate(outcome.naturalSubshells)
                sb = string(subsh)
                wa = TableStrings.floatList(10, outcome.orbitalInteraction[i,:])
                for  (j, sa)  in  enumerate(wa)
                    if       j == 1  sb = "     " * sb * "     " * sa
                    else             sb = "                    " * sa      end
                    println(stream,  sb)
                end
            end
            println(stream, "  ", TableStrings.hLine(nx))
        end

        if  settings.calc2pRDM
            println(stream, "\n  Reduced 2-particle RDM were calculated for the following subshells and are available by outcome.rho2p:\n")
            wa = TableStrings.subshellList(10, outcome.level.basis.subshells)
            nx = 120
            println(stream, "  ", TableStrings.hLine(nx))
            for  (i, sa)  in  enumerate(wa)
                if       i == 1  sb = "     Orb:  " * sa
                else             sb = "           " * sa      end
                println(stream,  sb)
            end
            println(stream, "  ", TableStrings.hLine(nx))

            # Compute the diagonal part of rho2p
            lenNO = length(outcome.level.basis.subshells);   rho_pprr = zeros(lenNO,lenNO)
            for  ip = 1:lenNO   for  ir = 1:lenNO     rho_pprr[ip,ir] = outcome.rho2p[ip,ip, ir,ir]     end   end

            println(stream, "\n  The diagonal part of the 2-particle RDM refers the mean (product) of occuaption numbers " *
                            "rho_pp,rr = (n_p * n_r)_av:\n")
            wa = TableStrings.subshellList(10, outcome.level.basis.subshells)
            nx = 120
            println(stream, "  ", TableStrings.hLine(nx))
            for  (i, sa)  in  enumerate(wa)
                if       i == 1  sb = "     Orb:  " * sa
                else             sb = "           " * sa      end
                println(stream,  sb)
            end
            println(stream, "  ", TableStrings.hLine(nx))

            for  (i, subsh)  in  enumerate(outcome.level.basis.subshells)
                sb = string(subsh)
                wa = TableStrings.floatList(10, outcome.orbitalInteraction[i,:])
                for  (j, sa)  in  enumerate(wa)
                    if       j == 1  sb = "     " * sb * "     " * sa
                    else             sb = "                    " * sa      end
                    println(stream,  sb)
                end
            end
            println(stream, "  ", TableStrings.hLine(nx))
        end

    end # outcomes

    return( nothing )
end


end # module
