
"""
`module  JAC.CrystalField`
    ... a submodel of JAC that contains all methods for computing the crystal-field (Stark)
        splitting of atomic levels in the field of an external point-charge lattice, following the
        point-charge model of Gaigalas & Kato, Comput. Phys. Commun. 261 (2021) 107772.

        The crystal-field potential is expanded in rank-k spherical tensors (Eqs. (1)-(2) of that
        paper) and its matrix elements are evaluated between the CSF/ASF representation already
        available from a JAC computation (Eq. (3)), reusing
        AngularMomentum.CL_reduced_me_sms(...) (one-electron reduced matrix elements of C^(k),
        with the l_a+l_b+k parity selection rule of Eq. (5) built in),
        SpinAngular.computeCoefficients(...) (many-electron spin-angular coefficients of a
        one-particle tensor operator) and RadialIntegrals.rkDiagonal(...) (the radial moment
        <r^k>_ab, valid for point charges outside the electron shell, i.e. the r<R branch of
        Eq. (2)), following the same architecture as Hfs.jl.

        The crystal-field source (currently: a point-charge lattice) is represented by
        CrystalField.AbstractCrystalFieldModel, so that a future model -- e.g. a superposition
        model, or a point-charge model with an empirical covalency/hybridization rescaling as in
        Uldry, Vernay & Delley, Phys. Rev. B 85, 125133 (2012) -- can be added later as a new
        subtype without changing the signatures of the functions below.

        A crystal field breaks the (2J+1)-fold degeneracy of a free-ion level and, if J-mixing is
        requested, also couples levels of different total J. The resulting Stark sublevels are no
        longer eigenstates of J^2 and therefore cannot be represented by a plain ManyElectron.Level
        (which assumes one definite J and M); instead they are collected in a
        CrystalField.CfMultiplet of CrystalField.CfLevel's, each of which is a mixture over the
        CrystalField.CfBasisVector's (one per M_J of every parent level under consideration) -- see
        CrystalField.computeRepresentation(...).

        Note (28-Jul-2026): in the standard |JM> basis, the crystal-field Hamiltonian matrix is
        Hermitian but, for anything less symmetric than a centrosymmetric/axially-aligned lattice,
        genuinely COMPLEX off-diagonal (its diagonal is always real) -- Gaigalas & Kato's own code
        diagonalizes it with LAPACK's general-complex solver ZGEEV for exactly this reason.
        CrystalField.computeInteractionMatrix therefore returns a Hermitian ComplexF64 matrix, and
        CrystalField.computeRepresentation diagonalizes it with
        LinearAlgebra.eigen(Hermitian(matrix)) rather than Basics.diagonalize(...) (which only
        accepts real Float64 matrices, appropriate for e.g. Hfs.jl but not here); real eigenvalues
        are still guaranteed by Hermiticity, but the eigenvectors (CfLevel.mc) are complex.
"""
module CrystalField


using LinearAlgebra, Printf, ..AngularMomentum, ..Basics, ..Defaults, ..ManyElectron, ..Radial, ..RadialIntegrals, ..SpinAngular


#################################################################################################################################
#################################################################################################################################

"""
`abstract type  CrystalField.AbstractCrystalFieldModel`
    ... the physical model used to represent the external crystal-field source acting on the
        atom/ion. Concrete subtypes carry only what differs between models; every computational
        function below dispatches on `model::AbstractCrystalFieldModel`, so a new model (e.g. a
        superposition model, or an Uldry-style covalency/hybridization correction) can be added
        later as a new subtype without touching these signatures.
"""
abstract type  AbstractCrystalFieldModel                                              end


"""
`struct  CrystalField.PointChargeModel  <:  CrystalField.AbstractCrystalFieldModel`
    ... the classical point-charge crystal-field model of Gaigalas & Kato, CPC 261 (2021) 107772.

    + scaleField   ::Float64   ... empirical overall scaling of the crystal-field strength;
                                   1.0 = unscaled ab-initio point-charge value. Provided so an
                                   Uldry-style (Phys. Rev. B 85, 125133 (2012)) empirical
                                   rescaling can be applied later without introducing a new model
                                   type.
"""
struct   PointChargeModel  <:  AbstractCrystalFieldModel
    scaleField     ::Float64
end


"""
`CrystalField.PointChargeModel()`  ... constructor for the unscaled (scaleField = 1.0) point-charge model.
"""
function PointChargeModel()
    PointChargeModel(1.0)
end


"""
`struct  CrystalField.PointCharge`
    ... one neighboring ion of the crystal lattice, in a spherical frame centered on the atom/ion
        under study.

    + charge   ::Float64   ... charge of the ion, in units of the elementary charge.
    + rho      ::Float64   ... distance of the ion from the central atom/ion, in Bohr radii.
    + theta    ::Float64   ... polar angle of the ion, in radian.
    + phi      ::Float64   ... azimuthal angle of the ion, in radian.
"""
struct   PointCharge
    charge     ::Float64
    rho        ::Float64
    theta      ::Float64
    phi        ::Float64
end


"""
`struct  CrystalField.Lattice`
    ... the set of point charges that together define the external crystal field.

    + ions            ::Array{CrystalField.PointCharge,1}   ... the point charges of the lattice.
    + symmetryLabel   ::String                              ... optional point-group label (e.g. "Oh", "C2"), for
                                                                 documentation/display only; no symmetry is
                                                                 exploited internally in the first implementation.
"""
struct   Lattice
    ions              ::Array{PointCharge,1}
    symmetryLabel     ::String
end


"""
`CrystalField.Lattice()`  ... constructor for an empty lattice (no external field).
"""
function Lattice()
    Lattice(PointCharge[], "")
end


# `Base.show(io::IO, lattice::CrystalField.Lattice)`  ... prepares a proper printout of lattice::CrystalField.Lattice.
function Base.show(io::IO, lattice::Lattice)
    println(io, "Crystal-field lattice \"$(lattice.symmetryLabel)\" of $(length(lattice.ions)) point charges.")
end


"""
`struct  CrystalField.CfBasisVector`
    ... one M_J sublevel of one parent (free-ion) level; the diagonalization basis of a crystal-field
        computation is the list of all CfBasisVector's belonging to the selected parent level(s) --
        just one level if no J-mixing is requested, several (possibly of different J) otherwise.

    + M             ::AngularM64   ... M_J projection of this sublevel.
    + parentLevel   ::Level        ... the (unperturbed) free-ion level this M_J value belongs to.
"""
struct  CfBasisVector
    M             ::AngularM64
    parentLevel   ::Level
end


"""
`struct  CrystalField.CfLevel`
    ... one crystal-field (Stark) sublevel, i.e. one eigenstate of the diagonalized crystal-field
        Hamiltonian. Since the crystal field breaks the M_J degeneracy of (and, under J-mixing,
        possibly also couples) the parent free-ion level(s), a CfLevel is in general a mixture over
        several CfBasisVector's and therefore has no single well-defined J or M -- unlike a plain
        ManyElectron.Level.

    + energy    ::Float64                             ... Stark sublevel energy.
    + cfBasis   ::Array{CrystalField.CfBasisVector,1}   ... the common basis of M_J sublevels for the whole CfMultiplet.
    + mc        ::Vector{ComplexF64}                   ... mixing coefficients of this sublevel with regard to cfBasis
                                                            (complex in general, see the module docstring note on
                                                            why the crystal-field matrix is not real for a generic
                                                            lattice).
"""
struct  CfLevel
    energy      ::Float64
    cfBasis     ::Array{CfBasisVector,1}
    mc          ::Vector{ComplexF64}
end


"""
`struct  CrystalField.CfMultiplet`
    ... the ordered list of crystal-field (Stark) sublevels that arise from splitting one or
        several parent free-ion levels in an external point-charge lattice.

    + name      ::String                       ... a name associated with this crystal-field multiplet.
    + cfLevels  ::Array{CrystalField.CfLevel,1}  ... the list of Stark sublevels.
"""
struct  CfMultiplet
    name        ::String
    cfLevels    ::Array{CfLevel,1}
end


# `Base.show(io::IO, cfLevel::CrystalField.CfLevel)`  ... prepares a proper printout of cfLevel::CrystalField.CfLevel.
function Base.show(io::IO, cfLevel::CfLevel)
    println(io, "Crystal-field sublevel:  energy = $(cfLevel.energy)")
end


# `Base.show(io::IO, multiplet::CrystalField.CfMultiplet)`  ... prepares a proper printout of multiplet::CrystalField.CfMultiplet.
function Base.show(io::IO, multiplet::CfMultiplet)
    println(io, "Crystal-field multiplet \"$(multiplet.name)\" of $(length(multiplet.cfLevels)) Stark sublevels.")
end


"""
`struct  Settings  <:  Basics.AbstractPropertySettings`
    ... defines the settings for computing the crystal-field (Stark) splitting of one or several
        atomic levels.

    + lattice          ::CrystalField.Lattice               ... the external point-charge lattice.
    + model            ::CrystalField.AbstractCrystalFieldModel  ... the crystal-field model to be used.
    + maxRank          ::Int64                               ... maximum tensor rank k retained in the multipole
                                                                   expansion of the crystal-field potential
                                                                   (k=0 is omitted throughout: it is a common
                                                                   energy shift of all levels and does not
                                                                   contribute to the splitting pattern).
    + includeJmixing   ::Bool                                 ... if true, allow mixing between the selected levels
                                                                   even if they belong to different total J.
    + printBefore      ::Bool                                 ... True if a list of selected levels is printed
                                                                   before the actual computations start.
    + levelSelection   ::LevelSelection                       ... specifies the levels for which the crystal-field
                                                                   splitting is to be computed.
"""
struct  Settings  <:  Basics.AbstractPropertySettings
    lattice          ::Lattice
    model            ::AbstractCrystalFieldModel
    maxRank          ::Int64
    includeJmixing   ::Bool
    printBefore      ::Bool
    levelSelection   ::LevelSelection
end


"""
`CrystalField.Settings()`  ... constructor for a Settings with an empty lattice and otherwise default values.
"""
function Settings()
    Settings( Lattice(), PointChargeModel(), 6, false, false, LevelSelection() )
end


"""
`CrystalField.Settings(set::CrystalField.Settings;`

        lattice=.., model=.., maxRank=.., includeJmixing=.., printBefore=.., levelSelection=..)

    ... keyword copy-constructor for re-defining selected values of a settings::CrystalField.Settings.
"""
function Settings(set::CrystalField.Settings;
        lattice::Union{Nothing,Lattice}=nothing,                     model::Union{Nothing,AbstractCrystalFieldModel}=nothing,
        maxRank::Union{Nothing,Int64}=nothing,                       includeJmixing::Union{Nothing,Bool}=nothing,
        printBefore::Union{Nothing,Bool}=nothing,                    levelSelection::Union{Nothing,LevelSelection}=nothing)
    if  isnothing(lattice)           latticex        = set.lattice        else   latticex        = lattice        end
    if  isnothing(model)             modelx          = set.model          else   modelx          = model          end
    if  isnothing(maxRank)           maxRankx        = set.maxRank        else   maxRankx        = maxRank        end
    if  isnothing(includeJmixing)    includeJmixingx = set.includeJmixing else   includeJmixingx = includeJmixing end
    if  isnothing(printBefore)       printBeforex    = set.printBefore    else   printBeforex    = printBefore    end
    if  isnothing(levelSelection)    levelSelectionx = set.levelSelection else   levelSelectionx = levelSelection end

    Settings( latticex, modelx, maxRankx, includeJmixingx, printBeforex, levelSelectionx )
end


# `Base.show(io::IO, settings::CrystalField.Settings)`  ... prepares a proper printout of settings::CrystalField.Settings.
function Base.show(io::IO, settings::CrystalField.Settings)
    println(io, "lattice:                  $(settings.lattice)  ")
    println(io, "model:                    $(settings.model)  ")
    println(io, "maxRank:                  $(settings.maxRank)  ")
    println(io, "includeJmixing:           $(settings.includeJmixing)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "levelSelection:           $(settings.levelSelection)  ")
end


"""
`struct  CrystalField.Outcome`
    ... defines a type to keep the outcome of a crystal-field computation for one parent level (or,
        under J-mixing, one group of jointly-diagonalized parent levels).

    + Jlevel        ::Level               ... (one of) the parent free-ion level(s) this outcome refers to.
    + cfMultiplet   ::CrystalField.CfMultiplet  ... the resulting Stark sublevels.
"""
struct  Outcome
    Jlevel          ::Level
    cfMultiplet     ::CfMultiplet
end


"""
`CrystalField.Outcome()`  ... constructor for an `empty` instance of CrystalField.Outcome.
"""
function Outcome()
    Outcome( Level(), CfMultiplet("", CfLevel[]) )
end


# `Base.show(io::IO, outcome::CrystalField.Outcome)`  ... prepares a proper printout of outcome::CrystalField.Outcome.
function Base.show(io::IO, outcome::CrystalField.Outcome)
    println(io, "Jlevel:                   J = $(outcome.Jlevel.J), parity = $(outcome.Jlevel.parity)  ")
    println(io, "cfMultiplet:              $(outcome.cfMultiplet)  ")
end


#################################################################################################################################
#################################################################################################################################


"""
`CrystalField.multipoleLatticeSum(lattice::CrystalField.Lattice, k::Int64, q::Int64)`
    ... computes the purely geometric lattice-sum coefficient of rank k, component q,

            A_kq = sum_i  charge_i / rho_i^(k+1) * (-1)^(q+1) sqrt(4pi/(2k+1)) Y^k_{-q}(theta_i, phi_i)

        summed over the ions of `lattice` (Eq. (1)-(2) of Gaigalas & Kato, CPC 261 (2021) 107772,
        with the r<R branch of Eq. (2) since every lattice ion lies outside the electron shell of
        the central atom/ion); independent of the atomic system. Uses AngularMomentum.sphericalYlm(...).
        A value::ComplexF64 is returned.

    + lattice   ::CrystalField.Lattice   ... the external point-charge lattice.
    + k         ::Int64                 ... tensor rank of the multipole term.
    + q         ::Int64                 ... tensor component, -k <= q <= k.
"""
function multipoleLatticeSum(lattice::Lattice, k::Int64, q::Int64)
    wa = ComplexF64(0)
    for  ion in lattice.ions
        wa = wa + ion.charge / ion.rho^(k+1) * AngularMomentum.sphericalYlm(k, -q, ion.theta, ion.phi)
    end
    return  (-1.0)^(q+1) * sqrt(4pi/(2k+1)) * wa
end


"""
`CrystalField.electrostaticIntegral(k::Int64, orbitalA::Radial.Orbital, orbitalB::Radial.Orbital, grid::Radial.Grid)`
    ... computes the one-electron reduced matrix element  [kappa_a||C^(k)||kappa_b] * <r^k>_ab
        (Eq. (4) of Gaigalas & Kato) of the rank-k crystal-field tensor operator between two
        relativistic orbitals, by combining AngularMomentum.CL_reduced_me_sms(...) (which already
        enforces the l_a+l_b+k parity selection rule of Eq. (5)) with
        RadialIntegrals.rkDiagonal(k, orbitalA, orbitalB, grid). A value::Float64 is returned.

    + k          ::Int64            ... tensor rank.
    + orbitalA   ::Radial.Orbital   ... one relativistic orbital.
    + orbitalB   ::Radial.Orbital   ... the other relativistic orbital.
    + grid       ::Radial.Grid      ... the radial grid.
"""
function electrostaticIntegral(k::Int64, orbitalA::Radial.Orbital, orbitalB::Radial.Orbital, grid::Radial.Grid)
    return  AngularMomentum.CL_reduced_me_sms(orbitalA.subshell, k, orbitalB.subshell) *
            RadialIntegrals.rkDiagonal(k, orbitalA, orbitalB, grid)
end


"""
`CrystalField.reducedMatrixElement(k::Int64, rLevel::Level, sLevel::Level, grid::Radial.Grid)`
    ... computes the rank-k reduced matrix element <rLevel||C^(k)||sLevel> between two ASF levels,
        by looping over all CSF pairs of the two levels' own basis, obtaining the many-electron
        spin-angular coefficients from SpinAngular.computeCoefficients(...) for a
        SpinAngular.OneParticleOperator(k, ...), contracting each with
        CrystalField.electrostaticIntegral(...) and finally sandwiching the resulting CSF x CSF
        matrix between the two levels' mixing vectors -- following exactly the pattern of
        Hfs.amplitude (module-Hfs.jl). As in Hfs.amplitude, each coeff.T is divided by
        sqrt(2*j_a+1) to undo the internal "GRASP-like" normalization that
        SpinAngular.computeCoefficients applies for rank>0 one-particle operators (see the note on
        Hfs.amplitude). No same-parity restriction is imposed here: for odd k the l_a+l_b+k
        selection rule already built into CL_reduced_me_sms allows rLevel and sLevel to belong to
        different overall parities, which is exactly the "ASF mixing with different parities"
        capability described in the Program Summary of Gaigalas & Kato. A value::Float64 is returned.

    + k        ::Int64   ... tensor rank.
    + rLevel   ::Level   ... bra level.
    + sLevel   ::Level   ... ket level.
    + grid     ::Radial.Grid   ... the radial grid.
"""
function reducedMatrixElement(k::Int64, rLevel::Level, sLevel::Level, grid::Radial.Grid)
    nr = length(rLevel.basis.csfs);   ns = length(sLevel.basis.csfs)
    matrix = zeros(Float64, nr, ns)
    opa    = SpinAngular.OneParticleOperator(k, Basics.plus, true)
    for  r = 1:nr
        for  s = 1:ns
            subshellList = sLevel.basis.subshells
            wa = SpinAngular.computeCoefficients(opa, rLevel.basis.csfs[r], sLevel.basis.csfs[s], subshellList)
            me = 0.
            for  coeff in wa
                ja2  = Basics.subshell_2j(coeff.a)
                tamp = electrostaticIntegral(k, rLevel.basis.orbitals[coeff.a], sLevel.basis.orbitals[coeff.b], grid)
                me   = me + coeff.T / sqrt(ja2 + 1) * tamp
            end
            matrix[r,s] = me
        end
    end
    return  transpose(rLevel.mc) * matrix * sLevel.mc
end


"""
`CrystalField.computeInteractionMatrix(levels::Array{Level,1}, lattice::CrystalField.Lattice, model::CrystalField.AbstractCrystalFieldModel, grid::Radial.Grid, maxRank::Int64)`
    ... builds the full crystal-field Hamiltonian matrix over the M_J-resolved basis of all
        `levels` (one CrystalField.CfBasisVector per M_J of every level in `levels`), following
        the Wigner-Eckart decomposition of Eq. (3): for every tensor rank k = 1, .., maxRank the
        M_J-independent reduced matrix element <levelP||C^(k)||levelQ> is computed once (via
        CrystalField.reducedMatrixElement) and reused for every pair of M_J values, each weighted
        by CrystalField.multipoleLatticeSum(lattice,k,q) and the appropriate 3-j symbol and phase.
        The off-diagonal (crystal-field) part is scaled by model.scaleField if `model` is a
        CrystalField.PointChargeModel. The matrix is generally COMPLEX (see the module docstring
        note) since CrystalField.multipoleLatticeSum(lattice,k,q) is complex for any lattice
        without the right symmetry to make it real; the returned matrix is explicitly
        Hermitized ((matrix+matrix')/2, adjoint not transpose) as a safeguard against residual
        floating-point asymmetry, since the upper and lower triangles are summed independently
        rather than mirrored.
        A tuple (matrix::Array{ComplexF64,2}, cfBasis::Array{CrystalField.CfBasisVector,1}) is returned.

    + levels   ::Array{Level,1}                             ... the parent (unperturbed) levels to be diagonalized jointly.
    + lattice  ::CrystalField.Lattice                        ... the external point-charge lattice.
    + model    ::CrystalField.AbstractCrystalFieldModel      ... the crystal-field model to be applied.
    + grid     ::Radial.Grid                                ... the radial grid.
    + maxRank  ::Int64                                      ... maximum tensor rank k to be included.
"""
function computeInteractionMatrix(levels::Array{Level,1}, lattice::Lattice, model::AbstractCrystalFieldModel, grid::Radial.Grid, maxRank::Int64)
    # Build the (level, M) basis
    cfBasis = CfBasisVector[]
    for  lev in levels
        J2 = Basics.twice(lev.J)
        for  m2 = -J2:2:J2
            push!(cfBasis, CfBasisVector(AngularM64(m2//2), lev))
        end
    end
    n = length(cfBasis)
    #
    # Pre-compute the M_J-independent reduced matrix elements <levelP||C^(k)||levelQ>, once per (rank, levelP, levelQ)
    redme = Dict{Tuple{Int64,Int64,Int64},Float64}()
    for  k = 1:maxRank
        for  levP in levels
            for  levQ in levels
                redme[(k, levP.index, levQ.index)] = reducedMatrixElement(k, levP, levQ, grid)
            end
        end
    end
    #
    scaleField = model isa PointChargeModel  ?  model.scaleField  :  1.0
    #
    matrix = zeros(ComplexF64, n, n)
    for  bp = 1:n
        for  bq = 1:n
            vp = cfBasis[bp];   vq = cfBasis[bq]
            if  vp.parentLevel.index == vq.parentLevel.index  &&  vp.M == vq.M
                matrix[bp,bq] = matrix[bp,bq] + vp.parentLevel.energy
            end
            Jp = vp.parentLevel.J;   Jq = vq.parentLevel.J;   wa = ComplexF64(0)
            for  k = 1:maxRank
                rme = redme[(k, vp.parentLevel.index, vq.parentLevel.index)]
                if  rme == 0.   continue    end
                for  q = -k:k
                    latt   = multipoleLatticeSum(lattice, k, q)
                    threej = AngularMomentum.Wigner_3j(Jp, AngularJ64(k), Jq, AngularM64(-vp.M.num//vp.M.den), AngularM64(q//1), vq.M)
                    if  threej == 0.   continue    end
                    phase = (-1.0)^Int64( (Basics.twice(Jp) - Basics.twice(vp.M)) ÷ 2 )
                    wa    = wa + latt * phase * sqrt(Basics.twice(Jp)+1) * threej * rme
                end
            end
            matrix[bp,bq] = matrix[bp,bq] + scaleField * wa
        end
    end
    matrix = (matrix + matrix') / 2
    #
    return (matrix, cfBasis)
end


"""
`CrystalField.computeRepresentation(levels::Array{Level,1}, lattice::CrystalField.Lattice, model::CrystalField.AbstractCrystalFieldModel, grid::Radial.Grid, maxRank::Int64)`
    ... builds the crystal-field representation for the given (possibly J-mixed) set of parent
        levels: calls CrystalField.computeInteractionMatrix(...) to obtain the full (complex
        Hermitian) Hamiltonian matrix over the M_J-resolved basis, diagonalizes it with
        LinearAlgebra.eigen(Hermitian(matrix)) -- NOT Basics.diagonalize(MatrixWithLinearAlgebra(), ...),
        which only accepts real matrices and is therefore unsuitable here, see the module docstring
        note -- and repacks the eigenpairs into a CrystalField.CfMultiplet of CrystalField.CfLevel's.
        Hermiticity guarantees real eigenvalues even though the matrix and eigenvectors are complex.
        A cfMultiplet::CrystalField.CfMultiplet is returned.

    + levels   ::Array{Level,1}                             ... the parent (unperturbed) levels to be diagonalized jointly.
    + lattice  ::CrystalField.Lattice                        ... the external point-charge lattice.
    + model    ::CrystalField.AbstractCrystalFieldModel      ... the crystal-field model to be applied.
    + grid     ::Radial.Grid                                ... the radial grid.
    + maxRank  ::Int64                                      ... maximum tensor rank k to be included.
"""
function computeRepresentation(levels::Array{Level,1}, lattice::Lattice, model::AbstractCrystalFieldModel, grid::Radial.Grid, maxRank::Int64)
    matrix, cfBasis = computeInteractionMatrix(levels, lattice, model, grid, maxRank)
    eigen = LinearAlgebra.eigen(LinearAlgebra.Hermitian(matrix))
    cfLevels = CfLevel[]
    for  ev = 1:length(eigen.values)
        push!(cfLevels, CfLevel(real(eigen.values[ev]), cfBasis, eigen.vectors[:,ev]))
    end
    return  CfMultiplet("crystal-field", cfLevels)
end


"""
`CrystalField.computeOutcomes(multiplet::Multiplet, lattice::CrystalField.Lattice, grid::Radial.Grid, settings::CrystalField.Settings)`
    ... the standard top-level driver of the module. Selects the levels of `multiplet` specified by
        settings.levelSelection (or just the lowest level if no selection is active) and computes
        their crystal-field splitting: if settings.includeJmixing == false, each selected level is
        diagonalized separately (its own M_J manifold only, one CrystalField.Outcome per level); if
        settings.includeJmixing == true, all selected levels are diagonalized jointly in a single,
        common M_J-resolved basis via CrystalField.computeRepresentation(...), and a single
        CrystalField.Outcome is returned. An Array{CrystalField.Outcome,1} is returned.

    + multiplet   ::Multiplet              ... the (unperturbed) atomic levels.
    + lattice     ::CrystalField.Lattice    ... the external point-charge lattice.
    + grid        ::Radial.Grid            ... the radial grid.
    + settings    ::CrystalField.Settings   ... the crystal-field computation settings.
"""
function computeOutcomes(multiplet::Multiplet, lattice::Lattice, grid::Radial.Grid, settings::Settings)
    if  settings.levelSelection.active
        levels = [ lev  for lev in multiplet.levels  if  lev.index in settings.levelSelection.indices ]
    else
        levels = [ multiplet.levels[1] ]
    end
    if  settings.printBefore
        println("CrystalField.computeOutcomes(): crystal-field splitting requested for $(length(levels)) level(s), " *
                 "J-mixing = $(settings.includeJmixing).")
    end
    #
    outcomes = Outcome[]
    if  settings.includeJmixing
        cfMultiplet = computeRepresentation(levels, lattice, settings.model, grid, settings.maxRank)
        push!(outcomes, Outcome(levels[1], cfMultiplet))
    else
        for  lev in levels
            cfMultiplet = computeRepresentation([lev], lattice, settings.model, grid, settings.maxRank)
            push!(outcomes, Outcome(lev, cfMultiplet))
        end
    end
    #
    return  outcomes
end


"""
`CrystalField.characteristicSplitting(cfMultiplet::CrystalField.CfMultiplet)`
    ... computes the "characteristic crystal-field splitting" (CXS) of Uldry, Vernay & Delley,
        Phys. Rev. B 85, 125133 (2012), Sec. II F: a single scalar descriptor of the overall
        strength of a crystal-field multiplet, defined there as the energy difference between the
        barycenter of the group of Stark sublevels lying above the largest gap in the spectrum and
        that of the group lying below it (their own examples use case-specific variants, e.g.
        "two lowest vs. three highest" for Cu2+/CuO, or "highest negative vs. average of positive"
        for Ti3+/LaTiO3; splitting the spectrum at its single largest gap generalizes both of these
        into one algorithm, since that gap is exactly where their examples split it by hand).
        A value::Float64 (in Hartree) is returned; 0.0 if cfMultiplet has fewer than 2 sublevels.

    + cfMultiplet   ::CrystalField.CfMultiplet   ... the Stark sublevels of one crystal-field outcome.
"""
function characteristicSplitting(cfMultiplet::CfMultiplet)
    energies = sort([ lev.energy  for lev in cfMultiplet.cfLevels ])
    n = length(energies)
    if  n < 2   return 0.0   end
    #
    gaps = [ energies[i+1] - energies[i]  for i = 1:n-1 ]
    imax = argmax(gaps)
    lowerGroup = energies[1:imax];   upperGroup = energies[imax+1:end]
    #
    return  sum(upperGroup)/length(upperGroup) - sum(lowerGroup)/length(lowerGroup)
end


"""
`CrystalField.characteristicSplitting(outcome::CrystalField.Outcome)`
    ... convenience method, equivalent to CrystalField.characteristicSplitting(outcome.cfMultiplet).
"""
function characteristicSplitting(outcome::Outcome)
    return  characteristicSplitting(outcome.cfMultiplet)
end


"""
`CrystalField.fitScaleField(levels::Array{Level,1}, lattice::CrystalField.Lattice, grid::Radial.Grid, targetSplitting::Float64, maxRank::Int64;`

        tolerance::Float64=1.0e-6, maxIterations::Int64=20)

    ... finds the CrystalField.PointChargeModel.scaleField value -- the S_xtal empirical rescaling
        parameter of Uldry, Vernay & Delley, Phys. Rev. B 85, 125133 (2012), Sec. IV -- for which
        CrystalField.characteristicSplitting(...) of the crystal-field representation of `levels`
        matches `targetSplitting` (in Hartree). Whenever `levels` contains a single parent level
        (no J-mixing), the interaction matrix's diagonal is a multiple of the identity (the shared
        level energy), so only the scaleField-proportional off-diagonal part affects the eigenvalue
        differences: CXS is then EXACTLY proportional to scaleField, and a single evaluation at
        scaleField=1.0 already gives the exact answer. Under J-mixing this proportionality is only
        approximate, so a short secant iteration (seeded from that exact linear estimate) refines it
        further. A scaleField::Float64 is returned.

    + levels            ::Array{Level,1}      ... the parent (unperturbed) levels to be diagonalized jointly.
    + lattice           ::CrystalField.Lattice ... the external point-charge lattice.
    + grid              ::Radial.Grid         ... the radial grid.
    + targetSplitting   ::Float64             ... the target characteristic splitting, in Hartree.
    + maxRank           ::Int64               ... maximum tensor rank k to be included.
    + tolerance         ::Float64             ... relative convergence tolerance on targetSplitting.
    + maxIterations     ::Int64               ... maximum number of secant iterations.
"""
function fitScaleField(levels::Array{Level,1}, lattice::Lattice, grid::Radial.Grid, targetSplitting::Float64, maxRank::Int64;
                        tolerance::Float64=1.0e-6, maxIterations::Int64=20)
    evaluate(s) = characteristicSplitting( computeRepresentation(levels, lattice, PointChargeModel(s), grid, maxRank) )
    #
    cxsAtUnitScale = evaluate(1.0)
    if  cxsAtUnitScale == 0.
        error("CrystalField.fitScaleField: characteristic splitting vanishes at scaleField=1.0 " *
              "(no allowed multipole term for this lattice/level combination); cannot fit.")
    end
    #
    s0 = targetSplitting / cxsAtUnitScale;   f0 = evaluate(s0) - targetSplitting
    s1 = s0 * 1.01;                          f1 = evaluate(s1) - targetSplitting
    #
    for  _ = 1:maxIterations
        if  abs(f1) < tolerance * abs(targetSplitting)   break    end
        slope = (f1 - f0) / (s1 - s0)
        sNew  = s1 - f1/slope
        s0, f0 = s1, f1
        s1 = sNew;   f1 = evaluate(s1) - targetSplitting
    end
    #
    return  s1
end


"""
`CrystalField.displayResults(stream::IO, outcomes::Array{CrystalField.Outcome,1})`
    ... prints the crystal-field (Stark) splitting pattern for every outcome to `stream`, sorted by
        increasing sublevel energy and given relative to the lowest sublevel of each outcome, both
        in Hartree and in cm^-1 (directly comparable to Table 1 of Gaigalas & Kato), followed by the
        characteristic crystal-field splitting (CXS) of Uldry, Vernay & Delley (Sec. II F). Nothing
        is returned.

    + stream     ::IO                              ... the output stream.
    + outcomes   ::Array{CrystalField.Outcome,1}    ... the computed crystal-field outcomes.
"""
function displayResults(stream::IO, outcomes::Array{Outcome,1})
    println(stream, "\n  Crystal-field (Stark) splitting:\n")
    for  outcome in outcomes
        println(stream, "  Parent level:  J = $(outcome.Jlevel.J), parity = $(outcome.Jlevel.parity), " *
                          "energy = $(outcome.Jlevel.energy) Hartree")
        sortedLevels = sort(outcome.cfMultiplet.cfLevels, by = lev -> lev.energy)
        e0 = sortedLevels[1].energy
        for  (i, lev)  in enumerate(sortedLevels)
            enCm = Defaults.convertUnits("energy: from atomic to Kayser", lev.energy - e0)
            @printf(stream, "    sublevel %3i:   E - E0 = %14.6f cm^-1   (%.10f Hartree)\n", i, enCm, lev.energy)
        end
        cxs      = characteristicSplitting(outcome)
        cxsCm    = Defaults.convertUnits("energy: from atomic to Kayser", cxs)
        @printf(stream, "    characteristic crystal-field splitting (CXS) = %.6f cm^-1  (%.10f Hartree)\n", cxsCm, cxs)
        println(stream, "")
    end
    return  nothing
end


end # module CrystalField
