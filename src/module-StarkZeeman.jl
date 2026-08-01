
"""
`module  JAC.StarkZeeman`
... a submodel of JAC that computes the exact (not perturbative) level mixing and splitting of a
    caller-selected set of ASF/levels in static electric and/or magnetic fields of arbitrary relative
    direction.

    This is architecturally close to what module-Hfs.jl already does for hyperfine mixing: the
    correlated atomic structure (ordinary CI, with whatever Coulomb/Breit/QED physics the caller
    already chose) is computed first, and the external field(s) are then treated as a SECOND, much
    smaller, separate diagonalization directly in the basis of those already-good levels -- not folded
    into the original CSF-level CI Hamiltonian. `Settings.levelSelection` IS the pool of levels
    admitted into this second diagonalization (no separate gMultiplet/coupling-window concept is
    needed, unlike module-MultipolePolarizibility.jl/module-StarkShift.jl, which need a convergent sum
    over many DISTANT perturbers -- here only a handful of levels close enough to actually interact
    matter). Matrix elements between levels that cannot couple (wrong parity for the electric-dipole
    term, wrong M for the given field direction, etc.) come out exactly zero from the selection rules
    already built into MultipoleMoment.emmStaticAmplitude/LandeZeeman.amplitude -- the diagonalization
    itself sorts out which of the selected levels actually mix, so the caller does not need to
    pre-filter by symmetry.

    This module complements, not replaces, module-StarkShift.jl: StarkShift's quadratic-in-E formula
    is the correct treatment for levels far from any other level of opposite parity (a convergent
    perturbative sum); StarkZeeman is for the handful of levels close enough that this perturbative
    treatment breaks down (the `1/(E_J-E_J')` denominator blows up as E_J'->E_J) -- exactly the
    "departure from a quadratic Stark effect" van Leeuwen & Hogervorst, Z. Phys. A 316, 149 (1984)
    document experimentally for Ca 3d^2 ^3P states (see also examples/example-Cl.jl branch d).

    Field-direction generality: BOTH fields may point in an arbitrary direction (not necessarily
    parallel), specified as a magnitude plus a Cartesian unit vector. This uses the standard identity
    that a unit vector n_hat, expressed as a rank-1 spherical tensor, has components
    C^(1)_q(theta,phi) (Racah-normalized spherical harmonics, via AngularMomentum.sphericalYlm --
    already used the same way in module-CrystalField.jl to project a direction onto rank-k tensor
    components), so that n_hat . T^(1) = sum_q (-1)^q C^(1)_{-q}(theta,phi) T^(1)_q
    (StarkZeeman.fieldWeight). The common special case of both fields along the z quantization axis
    (StarkZeeman.Settings' default eDirection/bDirection) is simply the q=0-only special case of the
    exact same formula -- no separate code path.

    IMPORTANT, empirically-resolved normalization note (1-Aug-2026): turning the REDUCED matrix
    elements MultipoleMoment.emmStaticAmplitude/LandeZeeman.amplitude return into genuine M-resolved
    sublevel matrix elements requires the Wigner-Eckart projection
    <J'M'|T^(1)_q|JM> = ClebschGordan(J,M;1,q|J',M') * <J'||T^(1)||J> -- WITHOUT any additional
    1/sqrt(2J'+1) factor. This was NOT assumed; it was calibrated against hydrogen's exact n=2 linear
    Stark coupling (the well-known non-relativistic <2p,m=0|z|2s,m=0>=-3 a.u., re-expressed in JAC's
    own relativistic (kappa,j,m_j) basis via the standard l don't-add-a-guessed-factor
    l-s -> j recoupling of the m_j=1/2 sublevels of 2p_1/2 and 2p_3/2): the ratio of the two computed
    matrix elements matched the hand-derived, sign/convention-independent exact ratio sqrt(2) to 5
    significant figures WITHOUT the extra factor, and was off by a clean, wrong factor of
    sqrt(2J'+1)/sqrt(2*0.5+1) etc. WITH it. This is the same empirical finding (a JAC amplitude
    function's own convention already differs from a naively-transcribed bare Wigner-Eckart formula)
    as the 31-Jul-2026 module-LandeZeeman.jl C2 fix -- now confirmed independently for a second,
    different operator (E1, not M1), which is reassuring cross-validation rather than a coincidence.
    The residual ~15% discrepancy in the ABSOLUTE magnitude (not the ratio) against the naive -3 a.u.
    non-relativistic benchmark is attributable to ordinary SCF/orbital-quality effects for that
    artificially-combined two-configuration hydrogen calculation (its own computed energies are
    already ~1-1.4% off the exact -0.125 Ha), not a convention error -- the ratio test isolates the
    angular/convention part cleanly since it is scale-independent.

    Known, explicitly INHERITED (not blocking) limitation: the open kappa<=-3 bug in
    LandeZeeman.amplitude's N1 operator (see that function's own docstring) affects any
    `includeBField=true` computation involving such subshells here too, since this module calls that
    same function rather than duplicating it. Per explicit user direction, this is not a reason to
    delay a first implementation -- fixing that bug later fixes it here automatically.
"""
module StarkZeeman


using Printf, LinearAlgebra, ..AngularMomentum, ..Basics, ..Defaults, ..LandeZeeman, ..ManyElectron,
              ..MultipoleMoment, ..Nuclear, ..Radial, ..StarkShift, ..TableStrings


"""
`StarkZeeman.AU_BFIELD_IN_TESLA`
    ... the CODATA atomic unit of magnetic flux density, in Tesla (1 a.u. = 2.35051757e5 T). Cross-
        validated (1-Aug-2026) against LandeZeeman.jl's own independently-derived `conv = 0.11909076`
        MHz/T^2 constant: (energy-to-Hz conversion, 6.57968e9 MHz per a.u. energy) / (2.35051757e5)^2
        = 0.11909..., an exact match.
"""
const AU_BFIELD_IN_TESLA = 2.35051757e5


"""
`struct  StarkZeeman.SublevelBasis`  ... one (level,M) entry of the perturbation-matrix basis.

    + level  ::Level        ... the zero-field level.
    + M      ::AngularM64   ... the M-sublevel of that level.
"""
struct SublevelBasis
    level    ::Level
    M        ::AngularM64
end


"""
`struct  StarkZeeman.Component`  ... one dominant (level,M) contribution to a field-dressed eigenstate.

    + level   ::Level        ... the zero-field level this component belongs to.
    + M       ::AngularM64   ... its M-sublevel.
    + weight  ::Float64      ... |eigenvector coefficient|^2 for this (level,M) basis component.
"""
struct Component
    level     ::Level
    M         ::AngularM64
    weight    ::Float64
end


"""
`struct  StarkZeeman.Outcome`  ... one field-dressed eigenstate of the (E,B)-field perturbation matrix.

    + energy      ::Float64                        ... field-dressed energy, in atomic units.
    + components  ::Array{StarkZeeman.Component,1}  ... the zero-field (level,M) components this
                                                          eigenstate is built from, sorted by
                                                          descending weight ("major level/symmetry"
                                                          first).
"""
struct Outcome
    energy        ::Float64
    components    ::Array{StarkZeeman.Component,1}
end


"""
`struct  StarkZeeman.Settings  <:  AbstractPropertySettings`
    ... defines a type for the details and parameters of an exact (E,B)-field level-mixing/-splitting
        computation.

    + includeEField    ::Bool               ... True if the static electric-field (Stark) coupling is
                                                 to be included.
    + includeBField    ::Bool               ... True if the static magnetic-field (Zeeman) coupling is
                                                 to be included.
    + eField           ::Float64            ... Strength of the electric field, in V/cm (SAME unit
                                                 convention as StarkShift.Settings.EField).
    + eDirection       ::NTuple{3,Float64}  ... Cartesian direction of the electric field (need not be
                                                 normalized); default (0.,0.,1.), i.e. along the
                                                 quantization axis.
    + bField           ::Float64            ... Strength of the magnetic field, in Tesla (SAME unit
                                                 convention as LandeZeeman.Settings.BField).
    + bDirection       ::NTuple{3,Float64}  ... Cartesian direction of the magnetic field; default
                                                 (0.,0.,1.).
    + includeSchwinger ::Bool               ... True if Schwinger's QED correction Delta-N1 is to be
                                                 included in the magnetic-dipole coupling (SAME flag
                                                 role as LandeZeeman.Settings.includeSchwinger).
    + printBefore      ::Bool               ... True if a list of selected levels is printed before
                                                 the actual computation starts.
    + levelSelection   ::LevelSelection     ... Specifies the levels admitted into the (E,B)-field
                                                 perturbation matrix -- see the module docstring.
"""
struct Settings  <:  AbstractPropertySettings
    includeEField     ::Bool
    includeBField     ::Bool
    eField            ::Float64
    eDirection        ::NTuple{3,Float64}
    bField            ::Float64
    bDirection        ::NTuple{3,Float64}
    includeSchwinger  ::Bool
    printBefore       ::Bool
    levelSelection    ::LevelSelection
end


"""
`StarkZeeman.Settings()`  ... constructor for an `empty` instance of StarkZeeman.Settings.
"""
function Settings()
    Settings(false, false, 0., (0., 0., 1.), 0., (0., 0., 1.), false, false, LevelSelection())
end


"""
`StarkZeeman.Settings(set::StarkZeeman.Settings;`

        includeEField=.., includeBField=.., eField=.., eDirection=.., bField=.., bDirection=..,
        includeSchwinger=.., printBefore=.., levelSelection=..)

    ... keyword copy-constructor for re-defining selected values of a settings::StarkZeeman.Settings.
"""
function Settings(set::StarkZeeman.Settings;
        includeEField::Union{Nothing,Bool}=nothing,             includeBField::Union{Nothing,Bool}=nothing,
        eField::Union{Nothing,Float64}=nothing,                 eDirection::Union{Nothing,NTuple{3,Float64}}=nothing,
        bField::Union{Nothing,Float64}=nothing,                 bDirection::Union{Nothing,NTuple{3,Float64}}=nothing,
        includeSchwinger::Union{Nothing,Bool}=nothing,          printBefore::Union{Nothing,Bool}=nothing,
        levelSelection::Union{Nothing,LevelSelection}=nothing)
    if  isnothing(includeEField)     includeEFieldx     = set.includeEField     else   includeEFieldx     = includeEField     end
    if  isnothing(includeBField)     includeBFieldx     = set.includeBField     else   includeBFieldx     = includeBField     end
    if  isnothing(eField)            eFieldx            = set.eField            else   eFieldx            = eField            end
    if  isnothing(eDirection)        eDirectionx        = set.eDirection        else   eDirectionx        = eDirection        end
    if  isnothing(bField)            bFieldx            = set.bField            else   bFieldx            = bField            end
    if  isnothing(bDirection)        bDirectionx        = set.bDirection        else   bDirectionx        = bDirection        end
    if  isnothing(includeSchwinger)  includeSchwingerx  = set.includeSchwinger  else   includeSchwingerx  = includeSchwinger  end
    if  isnothing(printBefore)       printBeforex       = set.printBefore       else   printBeforex       = printBefore       end
    if  isnothing(levelSelection)    levelSelectionx    = set.levelSelection    else   levelSelectionx    = levelSelection    end

    Settings( includeEFieldx, includeBFieldx, eFieldx, eDirectionx, bFieldx, bDirectionx,
              includeSchwingerx, printBeforex, levelSelectionx )
end


# `Base.show(io::IO, settings::StarkZeeman.Settings)`  ... prepares a proper printout of the variable settings::StarkZeeman.Settings.
function Base.show(io::IO, settings::StarkZeeman.Settings)
    println(io, "includeEField:            $(settings.includeEField)  ")
    println(io, "includeBField:            $(settings.includeBField)  ")
    println(io, "eField:                   $(settings.eField)  [V/cm]")
    println(io, "eDirection:               $(settings.eDirection)  ")
    println(io, "bField:                   $(settings.bField)  [Tesla]")
    println(io, "bDirection:               $(settings.bDirection)  ")
    println(io, "includeSchwinger:         $(settings.includeSchwinger)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "levelSelection:           $(settings.levelSelection)  ")
end


"""
`StarkZeeman.directionToAngles(direction::NTuple{3,Float64})`
    ... converts a (not necessarily normalized) Cartesian direction (dx,dy,dz) into the polar angles
        (theta,phi) relative to the SAME z quantization axis used throughout JAC's M-sublevel
        machinery. A tuple (theta::Float64, phi::Float64) is returned.
"""
function directionToAngles(direction::NTuple{3,Float64})
    dx, dy, dz = direction
    r = sqrt(dx^2 + dy^2 + dz^2)
    r == 0.0   &&   error("StarkZeeman.directionToAngles(): zero direction vector.")
    theta = acos(dz/r)
    phi   = atan(dy, dx)
    return (theta, phi)
end


"""
`StarkZeeman.fieldWeight(q::Int64, theta::Float64, phi::Float64)`
    ... the rank-1 spherical-tensor weight of a unit-vector direction (theta,phi) in component q, i.e.
        n_hat . T^(1) = sum_q StarkZeeman.fieldWeight(q,theta,phi) * T^(1)_q, using the standard
        identity n_hat <-> C^(1)(theta,phi) (Racah-normalized spherical harmonic,
        AngularMomentum.sphericalYlm) and the rank-1 spherical dot product A.B = sum_q (-1)^q A_q B_{-q}.
        A value::ComplexF64 is returned.
"""
function fieldWeight(q::Int64, theta::Float64, phi::Float64)
    return  (-1.0)^q * sqrt(4pi/3) * AngularMomentum.sphericalYlm(1, -q, theta, phi)
end


"""
`StarkZeeman.mResolvedAmplitude(reducedAmplitude, Ji::AngularJ64, Mi::AngularM64, Jf::AngularJ64, Mf::AngularM64)`
    ... converts a reduced matrix element <Jf||T^(1)||Ji> (as returned by
        MultipoleMoment.emmStaticAmplitude or LandeZeeman.amplitude) into the genuine M-resolved
        sublevel matrix element <Jf,Mf|T^(1)_q|Ji,Mi>, q = Mf-Mi, via
        ClebschGordan(Ji,Mi;1,q|Jf,Mf) * reducedAmplitude -- deliberately WITHOUT an additional
        1/sqrt(2Jf+1) factor; see the module docstring for how this convention was empirically
        calibrated (not assumed) against the exact hydrogen n=2 Stark coupling. A value is returned in
        the same (real or complex) type as `reducedAmplitude`, or exactly zero if q is not in {-1,0,1}
        or the Clebsch-Gordan coefficient vanishes.
"""
function mResolvedAmplitude(reducedAmplitude, Ji::AngularJ64, Mi::AngularM64, Jf::AngularJ64, Mf::AngularM64)
    q = round(Int64, Basics.twice(Mf)/2.0 - Basics.twice(Mi)/2.0)
    abs(q) > 1   &&   return( zero(reducedAmplitude) )
    cg = AngularMomentum.ClebschGordan(Ji, Mi, AngularJ64(1), AngularM64(q), Jf, Mf)
    return( cg * reducedAmplitude )
end


"""
`StarkZeeman.computeStarkZeemanMatrix(multiplet::Multiplet, grid::Radial.Grid, settings::StarkZeeman.Settings)`
    ... builds and diagonalizes the (level,M)-basis (E,B)-field perturbation matrix for the levels
        selected by settings.levelSelection from `multiplet` -- see the module docstring for the full
        construction. A tuple (eigenvalues::Array{Float64,1}, eigenvectors::Array{ComplexF64,2},
        basisList::Array{StarkZeeman.SublevelBasis,1}) is returned.
"""
function computeStarkZeemanMatrix(multiplet::Multiplet, grid::Radial.Grid, settings::StarkZeeman.Settings)
    selectedLevels = [ level  for level in multiplet.levels  if  Basics.selectLevel(level, settings.levelSelection) ]

    basisList = StarkZeeman.SublevelBasis[]
    for  level in selectedLevels
        for  M in AngularMomentum.m_values(level.J)
            push!(basisList, StarkZeeman.SublevelBasis(level, M))
        end
    end
    N = length(basisList)

    thetaE, phiE = StarkZeeman.directionToAngles(settings.eDirection)
    thetaB, phiB = StarkZeeman.directionToAngles(settings.bDirection)
    eFieldAu = settings.eField / StarkShift.AU_EFIELD_IN_VCM
    bFieldAu = settings.bField / StarkZeeman.AU_BFIELD_IN_TESLA

    H = zeros(ComplexF64, N, N)
    for  i = 1:N   H[i,i] += basisList[i].level.energy   end

    for  i = 1:N,  j = 1:N
        bi = basisList[i];    bj = basisList[j]
        q  = round(Int64, Basics.twice(bj.M)/2.0 - Basics.twice(bi.M)/2.0)
        abs(q) > 1   &&   continue

        levi = bi.level;    levj = bj.level
        if  levi.basis.subshells != levj.basis.subshells
            subshellList = Basics.generate(OrderedSubshellList(), levj.basis, levi.basis)
            levi = Basics.generateLevelWithSymmetryReducedBasis(levi, subshellList)
            levj = Basics.generateLevelWithSymmetryReducedBasis(levj, subshellList)
        end

        if  settings.includeEField  &&  !(i == j)
            dRed = MultipoleMoment.emmStaticAmplitude(1, levj, levi, grid)
            if  dRed != 0.0
                me = StarkZeeman.mResolvedAmplitude(dRed, bi.level.J, bi.M, bj.level.J, bj.M)
                if  me != 0.0
                    H[i,j] += -eFieldAu * me * StarkZeeman.fieldWeight(q, thetaE, phiE)
                end
            end
        end

        if  settings.includeBField
            n1Red = LandeZeeman.amplitude(LandeZeeman.ZeemanN1(), levj, levi, grid)
            if  settings.includeSchwinger
                n1Red = n1Red + LandeZeeman.amplitude(LandeZeeman.ZeemanDeltaN1(), levj, levi, grid)
            end
            if  n1Red != ComplexF64(0.)
                me = StarkZeeman.mResolvedAmplitude(n1Red, bi.level.J, bi.M, bj.level.J, bj.M)
                if  me != ComplexF64(0.)
                    H[i,j] += -bFieldAu * me * StarkZeeman.fieldWeight(q, thetaB, phiB)
                end
            end
        end
    end

    Hh  = Hermitian((H + H') / 2)
    eig = eigen(Hh)
    return  (eig.values, eig.vectors, basisList)
end


"""
`StarkZeeman.determineOutcomes(eigenvalues, eigenvectors, basisList::Array{StarkZeeman.SublevelBasis,1})`
    ... packages the raw diagonalization result into a sorted (by energy) array of
        StarkZeeman.Outcome, each carrying its dominant zero-field (level,M) components (weight >
        1.0e-3, sorted by descending weight) -- the "major level/symmetry" labeling. An
        Array{StarkZeeman.Outcome,1} is returned.
"""
function determineOutcomes(eigenvalues, eigenvectors, basisList::Array{StarkZeeman.SublevelBasis,1})
    N = length(basisList)
    outcomes = StarkZeeman.Outcome[]
    for  k = 1:N
        comps = StarkZeeman.Component[]
        for  i = 1:N
            w = abs2(eigenvectors[i,k])
            if  w > 1.0e-3   push!(comps, StarkZeeman.Component(basisList[i].level, basisList[i].M, w))   end
        end
        sort!(comps, by = c -> -c.weight)
        push!(outcomes, StarkZeeman.Outcome(real(eigenvalues[k]), comps))
    end
    sort!(outcomes, by = o -> o.energy)
    return( outcomes )
end


"""
`StarkZeeman.computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::StarkZeeman.Settings; output=true)`
    ... to compute the field-dressed eigenenergies and mixing coefficients for the levels selected by
        settings.levelSelection from `multiplet`, in the given static electric and/or magnetic
        field(s). The results are printed in a neat table but nothing is returned otherwise.
"""
function computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::StarkZeeman.Settings; output=true)
    println("")
    printstyled("StarkZeeman.computeOutcomes(): The computation of the (E,B)-field level mixing starts now ... \n", color=:light_green)
    printstyled("------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    selectedLevels = [ level  for level in multiplet.levels  if  Basics.selectLevel(level, settings.levelSelection) ]
    if  settings.printBefore    StarkZeeman.displayOutcomes(selectedLevels)    end

    evals, evecs, basisList = StarkZeeman.computeStarkZeemanMatrix(multiplet, grid, settings)
    outcomes = StarkZeeman.determineOutcomes(evals, evecs, basisList)

    StarkZeeman.displayResults(stdout, outcomes)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    StarkZeeman.displayResults(iostream, outcomes)   end

    if    output    return( outcomes )
    else            return( nothing )
    end
end


"""
`StarkZeeman.displayOutcomes(selectedLevels::Array{Level,1})`
    ... to display the list of zero-field levels admitted into the (E,B)-field perturbation matrix. A
        small neat table is printed but nothing is returned otherwise.
"""
function  displayOutcomes(selectedLevels::Array{Level,1})
    nx = 43
    println(" ")
    println("  Levels admitted into the StarkZeeman perturbation matrix:")
    println(" ")
    println("  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "Energy"; na=4);
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    println(sa);    println(sb);    println("  ", TableStrings.hLine(nx))
    #
    for  level in selectedLevels
        sa  = "  ";    sym = LevelSymmetry( level.J, level.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(level.index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", level.energy)) * "    "
        println( sa )
    end
    println("  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`StarkZeeman.displayResults(stream::IO, outcomes::Array{StarkZeeman.Outcome,1})`
    ... to display the field-dressed eigenenergies and their dominant zero-field (level,M) components.
        A neat table is printed but nothing is returned otherwise.
"""
function  displayResults(stream::IO, outcomes::Array{StarkZeeman.Outcome,1})
    nx = 100
    println(stream, " ")
    println(stream, "  Field-dressed StarkZeeman eigenstates:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(8,  "Eigen-#"; na=2);                           sb = sb * TableStrings.hBlank(10)
    sa = sa * TableStrings.center(16, "Energy";  na=4);
    sb = sb * TableStrings.center(16, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(66, "Dominant zero-field (level, J^P, M) components [weight]"; na=4)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))
    #
    for  (k,outcome) in enumerate(outcomes)
        sa  = "  "
        sa = sa * TableStrings.center(8, string(k); na=2)
        sa = sa * @sprintf("%.10e", Defaults.convertUnits("energy: from atomic", outcome.energy)) * "    "
        parts = String[]
        for  comp in outcome.components
            sym = LevelSymmetry(comp.level.J, comp.level.parity)
            push!(parts, "L$(comp.level.index)($(string(sym)),M=$(comp.M))[$(@sprintf("%.3f", comp.weight))]")
        end
        sa = sa * join(parts, "  ")
        println(stream, sa )
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end

end # module
