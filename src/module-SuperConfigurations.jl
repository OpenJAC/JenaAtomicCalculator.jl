
"""
`module  JAC.SuperConfigurations`
    ... a submodel of JAC that gives the Gaussian (two-moment) statistical description of ONE configuration's
        UNRESOLVED level spectrum -- a mean energy, a variance, and a total statistical weight -- reported
        WITHOUT diagonalizing the configuration's CI matrix. This is the building block a superconfiguration
        opacity (and, later, ionization-balance or NLTE-kinetics) scheme needs for an ion whose open-shell
        CSF count is too large for ordinary, level-by-level JAC machinery to resolve one level at a time.

        ONE-WAY DEPENDENCY, the same discipline every statistical/empirical module this project has built
        keeps: no process or `Plasma` module is called FROM here, and this module is never called FROM any
        process module either. A caller -- a future `Plasma.SuperConfigurationOpacityScheme` -- builds a
        `Group` and reads its `meanEnergy`/`variance`/`degeneracy` fields directly.

        THE MEAN IS REUSED, NOT REIMPLEMENTED. `BasicsAZ.meanEnergyInClosedForm` already gives a
        configuration's (2J+1)-weighted average energy in closed form, touching no CSF at all, independently
        validated to 1.6e-9 relative against the exact trace on nine configurations through `[Xe] 4f^7`
        (`module-BasicsAZ-inc-AG.jl`). This module calls it directly for `Group.meanEnergy` -- which is also
        why this module is included AFTER `module-BasicsAZ.jl` in `JenaAtomicCalculator.jl`, not alongside
        the other statistical modules built earlier.

        THE VARIANCE IS DELIBERATELY NOT YET THE SAME KIND OF CLOSED FORM. A true shell-pair-level closed
        form for the SECOND moment needs genuine Racah-algebra sum-rule reduction -- 6-j-symbol-level
        identities for sums of PRODUCTS of two-body matrix elements -- which this module's own first version
        chose not to derive from memory without a reference to check it against: the risk of a subtly wrong
        combinatorial factor, one that passes every test that does not happen to probe it, was judged higher
        than the value of rushing it. WHAT IS IMPLEMENTED INSTEAD IS EXACT AND STILL USEFUL:
        `varianceFromCsfMatrix` builds the full CSF-pair Hamiltonian matrix (the SAME basis
        `Basics.determineMeanEnergy` builds) and reads the variance from Tr(H^2), avoiding the O(n^3)
        DIAGONALIZATION step but not the O(n^2) matrix CONSTRUCTION -- exactly the role
        `Basics.determineMeanEnergy` played before `meanEnergyInClosedForm` replaced it, and the natural
        internal validation target for whatever shell-pair closed form is derived later. For a configuration
        too large even for this, the module is not yet usable; that is Stage 2, named here rather than hidden.

        SCOPE, DELIBERATELY NARROW, the same choice `AdjustedLevels` made for the same reason: ONE ordinary
        (non-relativistic) `Configuration` at a time, pure Coulomb (no Breit/QED).
"""
module SuperConfigurations

using ..Basics, ..BasicsAZ, ..ManyElectron, ..Nuclear, ..Radial, ..InteractionStrength, ..Hamiltonian


"""
`struct  SuperConfigurations.Group`
    ... the Gaussian (two-moment) statistical description of one configuration's UNRESOLVED level spectrum:
        its mean energy, variance, and total (2J+1)-weighted degeneracy (the statistical weight the group
        represents), computed without diagonalizing its CI matrix. See this module's own docstring for
        exactly how `meanEnergy` and `variance` are obtained, and what each still costs.

    + conf          ::Configuration   ... the configuration this group describes.
    + meanEnergy    ::Float64          ... (2J+1)-weighted mean energy, in Hartree.
    + variance      ::Float64          ... (2J+1)-weighted variance of the level energies, in Hartree^2.
    + degeneracy    ::Float64          ... total (2J+1)-weighted statistical weight, SUM_csf (2*J_csf+1).
"""
struct  Group
    conf            ::Configuration
    meanEnergy      ::Float64
    variance        ::Float64
    degeneracy      ::Float64
end


"""
`Base.show(io::IO, group::SuperConfigurations.Group)`
    ... prints the configuration, its mean energy, its width (sqrt(variance), the natural quantity to quote
        for a spread) and its degeneracy. Nothing is returned.
"""
function  Base.show(io::IO, group::Group)
    println(io, "SuperConfigurations.Group: ", string(group.conf))
    println(io, "  mean energy           = $(group.meanEnergy) Hartree")
    println(io, "  width (sqrt variance) = $(sqrt(max(group.variance, 0.0))) Hartree")
    println(io, "  degeneracy (SUM 2J+1) = $(group.degeneracy)")
end


"""
`SuperConfigurations.buildBasisAndCsfs(conf::Configuration, orbitals::Dict{Subshell,Orbital})`
    ... a private helper that builds the SAME relativistic-CSF basis `Basics.determineMeanEnergy` builds, so
        that `varianceFromCsfMatrix` and `BasicsAZ.meanEnergyInClosedForm` describe EXACTLY the same set of
        states. A (basis::Basis, csfList::Array{CsfR,1}) tuple is returned.
"""
function  buildBasisAndCsfs(conf::Configuration, orbitals::Dict{Subshell,Orbital})
    relconfList  = Basics.generateConfigurations(Basics.RelativisticConfigurations(), conf)
    subshellList = Basics.generateSubshellList(relconfList)
    csfList      = CsfR[]
    for  relconf in relconfList    append!(csfList, Basics.generateCsfRs(relconf, subshellList))    end
    NoElectrons  = sum( csfList[1].occupation )

    coreSubshellList = Subshell[]
    for  k = 1:length(subshellList)
        mocc = Basics.subshell_2j(subshellList[k]) + 1;    isFilled = true
        for  csf in csfList
            if  csf.occupation[k] != mocc    isFilled = false;    break    end
        end
        isFilled   &&   push!(coreSubshellList, subshellList[k])
    end
    basis = Basis(true, NoElectrons, subshellList, csfList, coreSubshellList, orbitals)

    return( (basis, csfList) )
end


"""
`SuperConfigurations.varianceFromCsfMatrix(conf::Configuration, orbitals::Dict{Subshell,Orbital},
                                             nm::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings,
                                             meanEnergy::Float64)`
    ... to determine the (2J+1)-weighted variance of a configuration's level energies from the full CSF-pair
        Hamiltonian matrix, via Tr(H^2), WITHOUT diagonalizing it. THIS IS NOT YET THE CLOSED FORM this
        module's own docstring names as future work: it still touches every CSF PAIR (O(n^2) matrix
        elements), only the O(n^3) diagonalization is avoided. Pairs of CSFs in different J^P blocks are
        skipped outright rather than relied on to return zero, since the Hamiltonian's own J/parity
        selection rule makes them exactly zero by construction. A value::Float64 (in Hartree^2) is returned.
"""
function  varianceFromCsfMatrix(conf::Configuration, orbitals::Dict{Subshell,Orbital}, nm::Nuclear.Model,
                                  grid::Radial.Grid, settings::AsfSettings, meanEnergy::Float64)
    basis, csfList = buildBasisAndCsfs(conf, orbitals)
    potential = Nuclear.nuclearPotential(nm, grid)
    cache     = InteractionStrength.XLCache()

    num2 = 0.;   den = 0.
    for  r = 1:length(csfList)
        wr  = Basics.twice(csfList[r].J) + 1
        den = den + wr
        for  s = 1:length(csfList)
            ( csfList[s].J != csfList[r].J   ||   csfList[s].parity != csfList[r].parity )   &&   continue
            Hrs  = Hamiltonian.matrixElement(basis, r, s, nm, grid, settings, potential, cache)
            num2 = num2 + wr * Hrs^2
        end
    end

    return( num2/den - meanEnergy^2 )
end


"""
`SuperConfigurations.buildGroup(conf::Configuration, orbitals::Dict{Subshell,Orbital}, nm::Nuclear.Model,
                                  grid::Radial.Grid, settings::AsfSettings)`
    ... to build the Gaussian statistical description of `conf`: the mean energy from
        `BasicsAZ.meanEnergyInClosedForm` (closed form, no CSF touched) and the variance from
        `SuperConfigurations.varianceFromCsfMatrix` (exact, CSF-pair-level, not yet closed form -- see this
        module's own docstring for why). A group::SuperConfigurations.Group is returned.
"""
function  buildGroup(conf::Configuration, orbitals::Dict{Subshell,Orbital}, nm::Nuclear.Model,
                       grid::Radial.Grid, settings::AsfSettings)
    meanE          = BasicsAZ.meanEnergyInClosedForm(conf, orbitals, nm, grid)
    varE           = varianceFromCsfMatrix(conf, orbitals, nm, grid, settings, meanE)
    basis, csfList = buildBasisAndCsfs(conf, orbitals)
    degeneracy     = sum( Basics.twice(csf.J) + 1   for csf in csfList )

    return( Group(conf, meanE, varE, degeneracy) )
end

end # module
