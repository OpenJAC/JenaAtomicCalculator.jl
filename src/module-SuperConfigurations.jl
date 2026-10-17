
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

        STAGE 2 STARTS HERE. `SuperConfigurations.TransitionArray` gives the SAME two-moment description for
        the UNRESOLVED electric-dipole (E1) transition array BETWEEN two configurations -- the building
        block `Plasma.SuperConfigurationOpacityScheme` needs to turn a population weighting into an actual
        opacity. It is built with the SAME discipline as `Group.variance`: EXACT, at CSF-pair level (no
        diagonalization of either configuration, and no closed-form Bauche-Arnoult moment formula reached
        for without a reference to check it against), via `SuperConfigurations.buildTransitionArray`. See
        that function's own docstring for the trace identities it uses and why no extra statistical-weight
        factor is needed beyond the dipole reduced matrix element's own Wigner-Eckart normalization.

        SCOPE, DELIBERATELY NARROW, the same choice `AdjustedLevels` made for the same reason: ONE ordinary
        (non-relativistic) `Configuration` at a time, pure Coulomb (no Breit/QED).
"""
module SuperConfigurations

using ..Basics, ..BasicsAZ, ..ManyElectron, ..Nuclear, ..Radial, ..InteractionStrength, ..Hamiltonian, ..SpinAngular,
      LinearAlgebra


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
`struct  SuperConfigurations.TransitionArray`
    ... the Gaussian (two-moment) statistical description of the UNRESOLVED electric-dipole (E1) transition
        array between two configurations' level spectra: the line-strength-weighted mean transition energy,
        its variance, and the array's total line strength, computed without diagonalizing either
        configuration's CI matrix. See `SuperConfigurations.buildTransitionArray` for exactly how these are
        obtained, and this module's own docstring for why they are EXACT rather than a closed form.

    + confFinal     ::Configuration   ... the configuration the transition ends in.
    + confInitial   ::Configuration   ... the configuration the transition starts from.
    + meanEnergy    ::Float64          ... line-strength-weighted mean of (E_final - E_initial), in Hartree;
                                            POSITIVE for an allowed absorption transition to a higher complex.
    + variance      ::Float64          ... line-strength-weighted variance of the transition energy, in
                                            Hartree^2.
    + totalStrength ::Float64          ... the array's total E1 line strength, SUM over every CSF pair of
                                            |<confFinal CSF || D^(1) || confInitial CSF>|^2, in atomic units.
"""
struct  TransitionArray
    confFinal       ::Configuration
    confInitial     ::Configuration
    meanEnergy      ::Float64
    variance        ::Float64
    totalStrength   ::Float64
end


"""
`Base.show(io::IO, ta::SuperConfigurations.TransitionArray)`
    ... prints the two configurations, the mean transition energy, its width (sqrt(variance)) and the
        array's total line strength. Nothing is returned.
"""
function  Base.show(io::IO, ta::TransitionArray)
    println(io, "SuperConfigurations.TransitionArray:  ", string(ta.confInitial), "  ->  ", string(ta.confFinal))
    println(io, "  mean transition energy = $(ta.meanEnergy) Hartree")
    println(io, "  width (sqrt variance)  = $(sqrt(max(ta.variance, 0.0))) Hartree")
    println(io, "  total line strength    = $(ta.totalStrength)")
end


"""
`SuperConfigurations.buildBasisAndCsfs(conf::Configuration, orbitals::Dict{Subshell,Orbital};
                                         subshellList::Union{Nothing,Array{Subshell,1}}=nothing)`
    ... a private helper that builds the SAME relativistic-CSF basis `Basics.determineMeanEnergy` builds, so
        that `varianceFromCsfMatrix` and `BasicsAZ.meanEnergyInClosedForm` describe EXACTLY the same set of
        states. If `subshellList` is given, `conf`'s CSFs are built on THAT list instead of generating one
        from `conf` alone -- the one case this is needed is `buildTransitionArray`, where two DIFFERENT
        configurations' CSFs must be indexed against their COMBINED subshell list before a transition
        operator can connect them. A (basis::Basis, csfList::Array{CsfR,1}) tuple is returned.
"""
function  buildBasisAndCsfs(conf::Configuration, orbitals::Dict{Subshell,Orbital};
                              subshellList::Union{Nothing,Array{Subshell,1}}=nothing)
    relconfList = Basics.generateConfigurations(Basics.RelativisticConfigurations(), conf)
    subshells   = isnothing(subshellList) ? Basics.generateSubshellList(relconfList) : subshellList
    csfList     = CsfR[]
    for  relconf in relconfList    append!(csfList, Basics.generateCsfRs(relconf, subshells))    end
    NoElectrons = sum( csfList[1].occupation )

    coreSubshellList = Subshell[]
    for  k = 1:length(subshells)
        mocc = Basics.subshell_2j(subshells[k]) + 1;    isFilled = true
        for  csf in csfList
            if  csf.occupation[k] != mocc    isFilled = false;    break    end
        end
        isFilled   &&   push!(coreSubshellList, subshells[k])
    end
    basis = Basis(true, NoElectrons, subshells, csfList, coreSubshellList, orbitals)

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


"""
`SuperConfigurations.buildHamiltonianMatrix(csfList::Array{CsfR,1}, basis::Basis, nm::Nuclear.Model,
                                              grid::Radial.Grid, settings::AsfSettings)`
    ... a private helper that fills the FULL CSF-pair Hamiltonian matrix of `basis` -- zero between different
        J^P blocks, exactly as `varianceFromCsfMatrix`'s own loop skips them -- needed here because
        `buildTransitionArray` contracts this matrix with the transition operator rather than only reading
        its trace. A matrix::Array{Float64,2} is returned.
"""
function  buildHamiltonianMatrix(csfList::Array{CsfR,1}, basis::Basis, nm::Nuclear.Model, grid::Radial.Grid,
                                   settings::AsfSettings)
    potential = Nuclear.nuclearPotential(nm, grid)
    cache     = InteractionStrength.XLCache()
    n         = length(csfList)
    H         = zeros(Float64, n, n)
    for  r = 1:n
        for  s = 1:n
            ( csfList[s].J != csfList[r].J   ||   csfList[s].parity != csfList[r].parity )   &&   continue
            H[r,s] = Hamiltonian.matrixElement(basis, r, s, nm, grid, settings, potential, cache)
        end
    end

    return( H )
end


"""
`SuperConfigurations.buildDipoleMatrix(csfListFinal::Array{CsfR,1}, basisFinal::Basis,
                                         csfListInitial::Array{CsfR,1}, basisInitial::Basis,
                                         subshellList::Array{Subshell,1}, grid::Radial.Grid)`
    ... a private helper that fills the CSF-pair electric-dipole (E1) REDUCED matrix element
        <csfListFinal[s] || D^(1) || csfListInitial[r]> for every pair, the two-configuration analogue of
        `buildHamiltonianMatrix` for a rank-1 operator instead of the rank-0 Hamiltonian. Reuses the SAME
        `InteractionStrength.eMultipole(1, ...)` radial integral `MultipoleMoment.emmStaticAmplitude` pairs
        with `coeff.T` (GRASP-dressed, Rule 18) -- NOT `InteractionStrength.MabEmission`, which is the
        retarded photon-emission operator and needs a transition energy this module does not yet have. The
        residual sqrt(2J_f+1) factor Rule 18 places AFTER the CSF loop is applied HERE PER ROW, using that
        row's OWN csf.J, since `csfListFinal` spans several J blocks at once rather than one level's fixed J
        -- the direct generalization of `varianceFromCsfMatrix`'s per-CSF `wr` weight to a rank-1 operator.
        A matrix::Array{Float64,2} of size (length(csfListFinal), length(csfListInitial)) is returned.
"""
function  buildDipoleMatrix(csfListFinal::Array{CsfR,1}, basisFinal::Basis, csfListInitial::Array{CsfR,1},
                              basisInitial::Basis, subshellList::Array{Subshell,1}, grid::Radial.Grid)
    nF  = length(csfListFinal);   nI = length(csfListInitial)
    D   = zeros(Float64, nF, nI)
    opa = SpinAngular.OneParticleOperator(1, Basics.multipoleParity(EmMultipole(1, true)))
    for  s = 1:nF
        wf = sqrt( Basics.twice(csfListFinal[s].J) + 1.0 )
        for  r = 1:nI
            wa = SpinAngular.computeCoefficients(opa, csfListFinal[s], csfListInitial[r], subshellList)
            me = 0.
            for  coeff in wa
                me = me + coeff.T * InteractionStrength.eMultipole(1, basisFinal.orbitals[coeff.a],
                                                                       basisInitial.orbitals[coeff.b], grid)
            end
            D[s,r] = me * wf
        end
    end

    return( D )
end


"""
`SuperConfigurations.buildTransitionArray(confFinal::Configuration, confInitial::Configuration,
                                            orbitals::Dict{Subshell,Orbital}, nm::Nuclear.Model,
                                            grid::Radial.Grid, settings::AsfSettings)`
    ... to build the Gaussian statistical description of the UNRESOLVED E1 transition array from
        `confInitial` to `confFinal`: the line-strength-weighted mean energy, variance and total line
        strength, EXACT at CSF-pair level, with no diagonalization of either configuration's Hamiltonian and
        no closed-form Bauche-Arnoult moment formula. The two configurations' CSFs are built on their
        COMBINED subshell list (`Basics.merge`), exactly as every two-configuration transition amplitude
        elsewhere in JAC already does, so the dipole operator can connect them.

        THE ALGEBRA, so it can be checked rather than trusted. Let H1, H2 be the ordinary CSF-pair
        Hamiltonian matrices of confInitial, confFinal (each block-diagonal in J, parity, with EXACT
        eigenvalues equal to that configuration's true level energies), and D the CSF-pair E1 REDUCED matrix
        element (Rule 18: already Wigner-Eckart-normalized, so |D_sr|^2 IS the standard line strength between
        the levels D's CSF indices expand into -- no extra statistical weight is needed anywhere below).
        Inserting the completeness relation of each configuration's OWN exact eigenvectors on either side of
        D, with NO approximation:

            S_total = Tr(D^T D)                                           -- the array's total line strength
            mean    = [Tr(H2 D D^T) - Tr(H1 D^T D)] / S_total              -- strength-weighted <E_final-E_initial>
            var     = [Tr(H2^2 D D^T) - 2 Tr(H1 D^T H2 D) + Tr(H1^2 D^T D)] / S_total - mean^2

        See `TestFrames.testModule_SuperConfigurations` for the brute-force check -- diagonalizing H1, H2
        directly and summing over the resulting level pairs -- that these traces reproduce the double sum
        rather than being a tidy-looking guess. A transitionArray::SuperConfigurations.TransitionArray is
        returned.
"""
function  buildTransitionArray(confFinal::Configuration, confInitial::Configuration,
                                 orbitals::Dict{Subshell,Orbital}, nm::Nuclear.Model, grid::Radial.Grid,
                                 settings::AsfSettings)
    Basics.extractFromConfiguration(Basics.GetParity(), confFinal) ==
        Basics.extractFromConfiguration(Basics.GetParity(), confInitial)   &&
        error("SuperConfigurations.buildTransitionArray(): confFinal and confInitial have the SAME parity; " *
              "an E1 transition array between them is identically zero.")

    relconfFinal   = Basics.generateConfigurations(Basics.RelativisticConfigurations(), confFinal)
    relconfInitial = Basics.generateConfigurations(Basics.RelativisticConfigurations(), confInitial)
    subshellList   = Basics.merge( Basics.generateSubshellList(relconfInitial),
                                     Basics.generateSubshellList(relconfFinal) )

    basisFinal,   csfListFinal   = buildBasisAndCsfs(confFinal,   orbitals; subshellList=subshellList)
    basisInitial, csfListInitial = buildBasisAndCsfs(confInitial, orbitals; subshellList=subshellList)

    H2 = buildHamiltonianMatrix(csfListFinal,   basisFinal,   nm, grid, settings)
    H1 = buildHamiltonianMatrix(csfListInitial, basisInitial, nm, grid, settings)
    D  = buildDipoleMatrix(csfListFinal, basisFinal, csfListInitial, basisInitial, subshellList, grid)

    DDt    = D * transpose(D)                 # (nF x nF)
    DtD    = transpose(D) * D                 # (nI x nI)
    Stotal = LinearAlgebra.tr(DtD)

    Stotal < 1.0e-14   &&
        error("SuperConfigurations.buildTransitionArray(): the total E1 line strength between $confInitial " *
              "and $confFinal is zero to machine precision -- check that they are connected by a genuine " *
              "single-electron E1 transition.")

    meanE = ( LinearAlgebra.tr(H2 * DDt) - LinearAlgebra.tr(H1 * DtD) ) / Stotal

    secondMoment = ( LinearAlgebra.tr(H2^2 * DDt)   -   2*LinearAlgebra.tr(H1 * transpose(D) * H2 * D)   +
                      LinearAlgebra.tr(H1^2 * DtD) ) / Stotal
    varE = secondMoment - meanE^2

    return( TransitionArray(confFinal, confInitial, meanE, varE, Stotal) )
end

end # module
