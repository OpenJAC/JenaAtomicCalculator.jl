
"""
`module  JAC.HfsPaschenBack`
... a submodel of JAC that computes the exact (not perturbative) hyperfine Paschen-Back splitting of a
    caller-selected set of hyperfine (F,M) levels in a static magnetic and/or electric field of arbitrary
    direction -- the Breit-Rabi problem, generalized beyond its textbook two-level closed form to an
    arbitrary manifold, possibly spanning several F and several fine-structure J at once.

    ARCHITECTURE. This module is to Hfs.jl what StarkZeeman.jl already is to the plain electronic (ASF)
    levels: a SECOND, much smaller diagonalization performed directly in the basis of already-good levels,
    not folded into the original CSF-level CI. Where StarkZeeman diagonalizes a field Hamiltonian in the
    basis of ASF (Level,M) sublevels, this module diagonalizes the SAME KIND of Hamiltonian in the basis of
    Hfs.HfLevel (F,M) sublevels -- each of which is ITSELF already a J-mixture within its own F, the ordinary
    output of Hfs.computeHyperfineRepresentation. Hfs.jl's own zero-field Hamiltonian is built block-diagonal
    in F on purpose (level mixing within one F only); an external field is exactly what breaks that block
    structure, which is the whole content of "hyperfine Paschen-Back".

    INPUT, DELIBERATELY NOT THREADED THROUGH `perform()`. Like StarkShift.Settings.gMultiplet /
    MultipolePolarizibility.Settings.gMultiplet, Settings.hfMultiplets is a list of Hfs.HfMultiplet built
    INDEPENDENTLY, ahead of time, by the caller (via Hfs.computeHyperfineMultiplet or Hfs.computeOutcomes),
    one call per fine-structure level whose hyperfine manifold is to participate. This is necessary rather
    than a stylistic choice: a hyperfine multiplet is itself the diagonalized output of a SEPARATE
    computation (Hfs's own zero-field Hamiltonian), and `perform()` only threads the bare ELECTRONIC
    multiplet of one Atomic.Computation into a propertySettings' computeOutcomes -- it has no mechanism to
    chain a second module's diagonalized output into a third. `computeOutcomes(multiplet, nm, grid,
    settings)` therefore accepts `multiplet`/`nm` only for the uniform dispatch signature every
    AbstractPropertySettings method shares; neither is read.

    THE TWO FIELD TERMS, AND WHERE EACH REDUCED MATRIX ELEMENT COMES FROM. The field Hamiltonian is
    H_field = g_J-type electronic Zeeman/Stark term + g_I-type nuclear Zeeman term, and the two are built by
    the SAME recoupling idea run in OPPOSITE directions:

    + ELECTRONIC term (magnetic N1/DeltaN1 from LandeZeeman.amplitude, or electric E1 from
      MultipoleMoment.emmStaticAmplitude): the operator acts on the electrons alone and leaves the nucleus
      untouched, so Hfs.recouplingElectronicOperator(spinI,Ja,Fa,Jb,Fb,L) -- already written, and already
      independently VERIFIED there via a 454-case sum rule -- turns the ordinary (J-basis) reduced matrix
      element into the needed (F-basis) one, with the nucleus as the spectator in its own 6-j.
    + NUCLEAR term (the small g_I mu_N B I_z Breit-Rabi correction): the mirror image,
      Hfs.recouplingNuclearOperator(J,Ia,Fa,Ib,Fb,L), with the roles of nucleus and electrons swapped -- the
      electronic LEVEL, not merely its J, must be IDENTICAL on both sides (the operator leaves the electrons
      untouched, which for a multi-level HfMultiplet means more than equal J). Its reduced matrix element is
      already available as Hfs.computeInteractionAmplitudeM(Basics.M1, isomer, isomer) -- the same nuclear
      magnetic-dipole reduced matrix element Hfs.jl already uses for the hyperfine A-constant, since mu^(1)
      is one physical operator whether it couples to the electrons' own field (hyperfine) or to an external
      uniform one (nuclear Zeeman).

    Both HfLevel objects entering a matrix element may be J-MIXTURES within their own F (not bare (I,J)F
    product states), so the actual matrix element is a double sum over each HfLevel's own `mc`-weighted
    `hfBasisVectors`, not a single recoupling call -- see HfsPaschenBack.recoupledElectronicMe/
    recoupledNuclearZeemanMe.

    ARBITRARY FIELD GEOMETRY IS NOT RE-IMPLEMENTED HERE. StarkZeeman.fieldWeight/directionToAngles already
    turn a Cartesian field direction into the q=-1,0,+1 spherical-tensor weights this module needs, and
    StarkZeeman.mResolvedAmplitude already projects a reduced matrix element onto M via plain Wigner-Eckart
    -- both work for F exactly as they do for J (F is an ordinary AngularJ64), so both are called directly,
    unmodified, from `..StarkZeeman`.

    VALIDATION TARGET. The textbook two-level hyperfine Paschen-Back case -- one J, F = I+-1/2 -- has the
    exact closed-form Breit-Rabi formula, independent of everything else in JAC; see examples/example-Cr.jl
    branch (a), the calibration this module is built against, the same role hydrogen's exact n=2 Stark ratio
    plays for StarkZeeman.
"""
module HfsPaschenBack


using Printf, LinearAlgebra, ..AngularMomentum, ..Basics, ..Defaults, ..Hfs, ..LandeZeeman, ..ManyElectron,
              ..MultipoleMoment, ..Nuclear, ..Radial, ..StarkShift, ..StarkZeeman, ..TableStrings


"""
`struct  HfsPaschenBack.Component`  ... one dominant zero-field Hfs.HfLevel contribution to a field-dressed eigenstate.

    + hfLevel  ::Hfs.HfLevel  ... the zero-field hyperfine level this component belongs to.
    + weight   ::Float64      ... |eigenvector coefficient|^2 for this hfLevel basis component.
"""
struct Component
    hfLevel    ::Hfs.HfLevel
    weight     ::Float64
end


"""
`struct  HfsPaschenBack.Outcome`  ... one field-dressed eigenstate of the hyperfine (B,E)-field perturbation matrix.

    + energy      ::Float64                           ... field-dressed energy, in atomic units.
    + components  ::Array{HfsPaschenBack.Component,1}  ... the zero-field Hfs.HfLevel components this eigenstate is
                                                             built from, sorted by descending weight.
"""
struct Outcome
    energy        ::Float64
    components    ::Array{HfsPaschenBack.Component,1}
end


"""
`struct  HfsPaschenBack.Settings  <:  AbstractPropertySettings`
    ... defines a type for the details and parameters of an exact hyperfine Paschen-Back (Breit-Rabi) computation.

    + includeBField        ::Bool                     ... True if the static magnetic-field coupling is to be included.
    + bField               ::Float64                  ... Strength of the magnetic field, in Tesla (same convention as
                                                             LandeZeeman.Settings.BField/StarkZeeman.Settings.bField).
    + bDirection            ::NTuple{3,Float64}        ... Cartesian direction of the magnetic field; default (0.,0.,1.).
    + includeNuclearZeeman  ::Bool                     ... True if the small nuclear-Zeeman term (g_I mu_N B I_z) is to be
                                                             added to the electronic one -- the term that distinguishes a
                                                             genuine Breit-Rabi diagram from its leading-order approximation.
    + includeEField         ::Bool                     ... True if the static electric-field (Stark) coupling is to be included.
    + eField                ::Float64                  ... Strength of the electric field, in V/cm (same convention as
                                                             StarkShift.Settings.EField/StarkZeeman.Settings.eField).
    + eDirection             ::NTuple{3,Float64}        ... Cartesian direction of the electric field; default (0.,0.,1.).
    + includeSchwinger       ::Bool                     ... True if Schwinger's QED correction Delta-N1 is included in the
                                                             electronic magnetic-dipole coupling (same role as
                                                             LandeZeeman.Settings.includeSchwinger/StarkZeeman.Settings.includeSchwinger).
    + hfMultiplets           ::Array{Hfs.HfMultiplet,1} ... The zero-field hyperfine multiplet(s) admitted into the
                                                             perturbation matrix -- see the module docstring; built
                                                             independently, ahead of time, by the caller.
    + printBefore            ::Bool                     ... True if a list of selected hyperfine levels is printed before
                                                             the actual computation starts.
"""
struct Settings  <:  AbstractPropertySettings
    includeBField        ::Bool
    bField               ::Float64
    bDirection           ::NTuple{3,Float64}
    includeNuclearZeeman ::Bool
    includeEField        ::Bool
    eField               ::Float64
    eDirection           ::NTuple{3,Float64}
    includeSchwinger     ::Bool
    hfMultiplets         ::Array{Hfs.HfMultiplet,1}
    printBefore          ::Bool
end


"""
`HfsPaschenBack.Settings()`  ... constructor for an `empty` instance of HfsPaschenBack.Settings.
"""
function Settings()
    Settings(false, 0., (0., 0., 1.), false, false, 0., (0., 0., 1.), false, Hfs.HfMultiplet[], false)
end


"""
`HfsPaschenBack.Settings(set::HfsPaschenBack.Settings;`

        includeBField=.., bField=.., bDirection=.., includeNuclearZeeman=.., includeEField=.., eField=..,
        eDirection=.., includeSchwinger=.., hfMultiplets=.., printBefore=..)

    ... keyword copy-constructor for re-defining selected values of a settings::HfsPaschenBack.Settings.
"""
function Settings(set::HfsPaschenBack.Settings;
        includeBField::Union{Nothing,Bool}=nothing,             bField::Union{Nothing,Float64}=nothing,
        bDirection::Union{Nothing,NTuple{3,Float64}}=nothing,    includeNuclearZeeman::Union{Nothing,Bool}=nothing,
        includeEField::Union{Nothing,Bool}=nothing,              eField::Union{Nothing,Float64}=nothing,
        eDirection::Union{Nothing,NTuple{3,Float64}}=nothing,    includeSchwinger::Union{Nothing,Bool}=nothing,
        hfMultiplets::Union{Nothing,Array{Hfs.HfMultiplet,1}}=nothing,   printBefore::Union{Nothing,Bool}=nothing)
    if  isnothing(includeBField)         includeBFieldx        = set.includeBField        else   includeBFieldx        = includeBField        end
    if  isnothing(bField)                bFieldx               = set.bField               else   bFieldx               = bField               end
    if  isnothing(bDirection)            bDirectionx           = set.bDirection            else   bDirectionx           = bDirection           end
    if  isnothing(includeNuclearZeeman)  includeNuclearZeemanx = set.includeNuclearZeeman  else   includeNuclearZeemanx = includeNuclearZeeman  end
    if  isnothing(includeEField)         includeEFieldx        = set.includeEField         else   includeEFieldx        = includeEField        end
    if  isnothing(eField)                eFieldx               = set.eField               else   eFieldx               = eField               end
    if  isnothing(eDirection)            eDirectionx           = set.eDirection            else   eDirectionx           = eDirection           end
    if  isnothing(includeSchwinger)      includeSchwingerx     = set.includeSchwinger      else   includeSchwingerx     = includeSchwinger      end
    if  isnothing(hfMultiplets)          hfMultipletsx         = set.hfMultiplets          else   hfMultipletsx         = hfMultiplets          end
    if  isnothing(printBefore)           printBeforex          = set.printBefore           else   printBeforex          = printBefore           end

    Settings( includeBFieldx, bFieldx, bDirectionx, includeNuclearZeemanx, includeEFieldx, eFieldx, eDirectionx,
              includeSchwingerx, hfMultipletsx, printBeforex )
end


# `Base.show(io::IO, settings::HfsPaschenBack.Settings)`  ... prepares a proper printout of the variable settings::HfsPaschenBack.Settings.
function Base.show(io::IO, settings::HfsPaschenBack.Settings)
    println(io, "includeBField:            $(settings.includeBField)  ")
    println(io, "bField:                   $(settings.bField)  [Tesla]")
    println(io, "bDirection:               $(settings.bDirection)  ")
    println(io, "includeNuclearZeeman:     $(settings.includeNuclearZeeman)  ")
    println(io, "includeEField:            $(settings.includeEField)  ")
    println(io, "eField:                   $(settings.eField)  [V/cm]")
    println(io, "eDirection:               $(settings.eDirection)  ")
    println(io, "includeSchwinger:         $(settings.includeSchwinger)  ")
    println(io, "hfMultiplets:             $(length(settings.hfMultiplets)) given  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
end


"""
`HfsPaschenBack.expandMSublevels(hfMultiplets::Array{Hfs.HfMultiplet,1})`
    ... expands each Hfs.HfLevel of every given Hfs.HfMultiplet into its full set of 2F+1 (F,M) sublevels.

        Hfs.computeHyperfineRepresentation deliberately stores only ONE HfLevel per F, carrying the placeholder
        M = F rather than a real projection -- its own comment explains why: at zero field the 2F+1 sublevels are
        exactly degenerate and M cancels by the Wigner-Eckart theorem in every quantity Hfs.jl itself computes, so
        storing all of them would be pure duplication. An external field is exactly the thing that "becomes
        meaningful" once M stops being a mere label, which is this module's reason to exist, so the expansion
        happens HERE rather than in Hfs.jl: the energy, hfBasisVectors and mc of a given F are UNCHANGED across its
        2F+1 copies (the zero-field Hamiltonian is a scalar in M), and the only thing that varies is the M label
        itself -- the reduced matrix elements HfsPaschenBack.recoupledElectronicMe/recoupledNuclearZeemanMe compute
        from hfBasisVectors/mc are themselves M-independent; only the final Wigner-Eckart projection
        (StarkZeeman.mResolvedAmplitude) needs the true M of each copy. An Array{Hfs.HfLevel,1} is returned, with
        fresh, sequential `index` fields (the original Hfs-assigned index is no longer unique once a level has been
        expanded into several M-copies).
"""
function expandMSublevels(hfMultiplets::Array{Hfs.HfMultiplet,1})
    basisList = Hfs.HfLevel[];    idx = 0
    for  hfMultiplet in hfMultiplets,  hfLevel in hfMultiplet.hfLevels,  M in AngularMomentum.m_values(hfLevel.F)
        idx = idx + 1
        push!(basisList, Hfs.HfLevel(hfLevel.F, M, hfLevel.parity, idx, hfLevel.energy,
                                      hfLevel.hfBasisVectors, hfLevel.mc))
    end

    return( basisList )
end


"""
`HfsPaschenBack.recoupledElectronicMe(reducedMeOf::Function, levi::Hfs.HfLevel, levj::Hfs.HfLevel, grid::Radial.Grid)`
    ... the F-coupled reduced matrix element <levj||T^(1)||levi> of a rank-1 operator that acts on the ELECTRONS alone
        (`reducedMeOf(bra::Level, ket::Level)` is the caller-supplied ordinary J-basis reduced matrix element, e.g.
        LandeZeeman.amplitude(ZeemanN1(),bra,ket,grid) or MultipoleMoment.emmStaticAmplitude(1,bra,ket,grid)), obtained
        by recoupling that J-basis element into the (I,J)->F basis via Hfs.recouplingElectronicOperator and summing
        over each HfLevel's own mc-weighted mixture of hfBasisVectors -- see the module docstring. A value::ComplexF64
        is returned, exactly 0. if no basis-vector pair of `levi`/`levj` shares a common nuclear state (the operator
        conserves the isomer, since it does not touch the nucleus).
"""
function recoupledElectronicMe(reducedMeOf::Function, levi::Hfs.HfLevel, levj::Hfs.HfLevel, grid::Radial.Grid)
    me = ComplexF64(0.)
    for  (k, vk)  in  enumerate(levi.hfBasisVectors),  (l, vl)  in  enumerate(levj.hfBasisVectors)
        mcw = levi.mc[k] * levj.mc[l]
        if  abs(mcw) < 1.0e-10  ||  vk.isomer != vl.isomer     continue    end

        lhsLevel = vk.levelJ;    rhsLevel = vl.levelJ
        if  lhsLevel.basis.subshells != rhsLevel.basis.subshells
            subshellList = Basics.generate(OrderedSubshellList(), rhsLevel.basis, lhsLevel.basis)
            lhsLevel = Basics.generateLevelWithSymmetryReducedBasis(lhsLevel, subshellList)
            rhsLevel = Basics.generateLevelWithSymmetryReducedBasis(rhsLevel, subshellList)
        end
        redMe = reducedMeOf(rhsLevel, lhsLevel)                        ## <rhsLevel(bra) || op || lhsLevel(ket)>
        redMe == zero(redMe)   &&   continue

        recpl = Hfs.recouplingElectronicOperator(vk.isomer.spinI, vk.levelJ.J, vk.F, vl.levelJ.J, vl.F, 1)
        me = me + mcw * recpl * redMe
    end

    return( me )
end


"""
`HfsPaschenBack.recoupledNuclearZeemanMe(levi::Hfs.HfLevel, levj::Hfs.HfLevel)`
    ... the F-coupled reduced matrix element <levj||mu^(1)||levi> of the nuclear-Zeeman operator, which acts on the
        nucleus alone and therefore requires the ELECTRONIC level -- not merely its J -- to be literally identical on
        both sides (Hfs.recouplingNuclearOperator's spectator), summed over each HfLevel's own mc-weighted mixture of
        hfBasisVectors exactly as HfsPaschenBack.recoupledElectronicMe does for the electronic term. The nuclear
        reduced matrix element itself is Hfs.computeInteractionAmplitudeM(Basics.M1,isomer,isomer), the same one
        Hfs.jl already uses for the hyperfine A-constant, converted here from nuclear magnetons to atomic units. A
        value::Float64 is returned, exactly 0. if no basis-vector pair shares both the same electronic level and the
        same nuclear state.
"""
function recoupledNuclearZeemanMe(levi::Hfs.HfLevel, levj::Hfs.HfLevel)
    me = 0.0
    for  (k, vk)  in  enumerate(levi.hfBasisVectors),  (l, vl)  in  enumerate(levj.hfBasisVectors)
        mcw = levi.mc[k] * levj.mc[l]
        if  abs(mcw) < 1.0e-10  ||  vk.levelJ != vl.levelJ  ||  vk.isomer != vl.isomer     continue    end

        muRed = Hfs.computeInteractionAmplitudeM(Basics.M1, vk.isomer, vl.isomer) *
                Defaults.convertUnits("moment: from nuclear magneton to atomic", 1.0)
        recpl = Hfs.recouplingNuclearOperator(vk.levelJ.J, vk.isomer.spinI, vk.F, vl.isomer.spinI, vl.F, 1)
        me = me + mcw * recpl * muRed
    end

    return( me )
end


"""
`HfsPaschenBack.computeHfsZeemanMatrix(settings::HfsPaschenBack.Settings, grid::Radial.Grid)`
    ... builds and diagonalizes the Hfs.HfLevel-basis hyperfine Paschen-Back perturbation matrix for every HfLevel of
        every Hfs.HfMultiplet in settings.hfMultiplets -- see the module docstring for the full construction. A tuple
        (eigenvalues::Array{Float64,1}, eigenvectors::Array{ComplexF64,2}, basisList::Array{Hfs.HfLevel,1}) is returned.
"""
function computeHfsZeemanMatrix(settings::HfsPaschenBack.Settings, grid::Radial.Grid)
    basisList = HfsPaschenBack.expandMSublevels(settings.hfMultiplets)
    N = length(basisList)

    thetaB, phiB = StarkZeeman.directionToAngles(settings.bDirection)
    thetaE, phiE = StarkZeeman.directionToAngles(settings.eDirection)
    bFieldAu = settings.bField / StarkZeeman.AU_BFIELD_IN_TESLA
    eFieldAu = settings.eField / StarkShift.AU_EFIELD_IN_VCM

    H = zeros(ComplexF64, N, N)
    for  i = 1:N    H[i,i] += basisList[i].energy    end

    for  i = 1:N,  j = 1:N
        levi = basisList[i];    levj = basisList[j]
        q  = round(Int64, Basics.twice(levj.M)/2.0 - Basics.twice(levi.M)/2.0)
        abs(q) > 1   &&   continue

        if  settings.includeBField
            n1Me = HfsPaschenBack.recoupledElectronicMe(
                       (bra,ket) -> settings.includeSchwinger ?
                           LandeZeeman.amplitude(LandeZeeman.ZeemanN1(), bra, ket, grid) +
                           LandeZeeman.amplitude(LandeZeeman.ZeemanDeltaN1(), bra, ket, grid) :
                           LandeZeeman.amplitude(LandeZeeman.ZeemanN1(), bra, ket, grid),
                       levi, levj, grid)
            if  n1Me != ComplexF64(0.)
                me = StarkZeeman.mResolvedAmplitude(n1Me, levi.F, levi.M, levj.F, levj.M)
                if  me != ComplexF64(0.)    H[i,j] += -bFieldAu * me * StarkZeeman.fieldWeight(q, thetaB, phiB)    end
            end

            if  settings.includeNuclearZeeman
                nucMe = HfsPaschenBack.recoupledNuclearZeemanMe(levi, levj)
                if  nucMe != 0.0
                    me = StarkZeeman.mResolvedAmplitude(nucMe, levi.F, levi.M, levj.F, levj.M)
                    if  me != 0.0    H[i,j] += -bFieldAu * me * StarkZeeman.fieldWeight(q, thetaB, phiB)    end
                end
            end
        end

        if  settings.includeEField  &&  !(i == j)
            dMe = HfsPaschenBack.recoupledElectronicMe(
                      (bra,ket) -> MultipoleMoment.emmStaticAmplitude(1, bra, ket, grid), levi, levj, grid)
            if  dMe != ComplexF64(0.)
                me = StarkZeeman.mResolvedAmplitude(dMe, levi.F, levi.M, levj.F, levj.M)
                if  me != ComplexF64(0.)    H[i,j] += -eFieldAu * me * StarkZeeman.fieldWeight(q, thetaE, phiE)    end
            end
        end
    end

    Hh  = Hermitian((H + H') / 2)
    eig = eigen(Hh)
    return  (eig.values, eig.vectors, basisList)
end


"""
`HfsPaschenBack.determineOutcomes(eigenvalues, eigenvectors, basisList::Array{Hfs.HfLevel,1})`
    ... packages the raw diagonalization result into a sorted (by energy) array of HfsPaschenBack.Outcome, each
        carrying its dominant zero-field Hfs.HfLevel components (weight > 1.0e-3, sorted by descending weight). An
        Array{HfsPaschenBack.Outcome,1} is returned.
"""
function determineOutcomes(eigenvalues, eigenvectors, basisList::Array{Hfs.HfLevel,1})
    N = length(basisList)
    outcomes = HfsPaschenBack.Outcome[]
    for  k = 1:N
        comps = HfsPaschenBack.Component[]
        for  i = 1:N
            w = abs2(eigenvectors[i,k])
            if  w > 1.0e-3   push!(comps, HfsPaschenBack.Component(basisList[i], w))   end
        end
        sort!(comps, by = c -> -c.weight)
        push!(outcomes, HfsPaschenBack.Outcome(real(eigenvalues[k]), comps))
    end
    sort!(outcomes, by = o -> o.energy)
    return( outcomes )
end


"""
`HfsPaschenBack.computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::HfsPaschenBack.Settings; output=true)`
    ... to compute the field-dressed hyperfine Paschen-Back eigenenergies and mixing coefficients for every HfLevel of
        every Hfs.HfMultiplet in settings.hfMultiplets, in the given static magnetic and/or electric field(s).
        `multiplet`/`nm` are accepted only for the uniform AbstractPropertySettings dispatch signature and are not
        read -- see the module docstring. The results are printed in a neat table; an Array{HfsPaschenBack.Outcome,1}
        is returned if `output`, nothing otherwise.
"""
function computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::HfsPaschenBack.Settings; output=true)
    println("")
    printstyled("HfsPaschenBack.computeOutcomes(): The computation of the hyperfine Paschen-Back mixing starts now ... \n", color=:light_green)
    printstyled("---------------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    if  settings.printBefore    HfsPaschenBack.displayOutcomes(settings.hfMultiplets)    end

    evals, evecs, basisList = HfsPaschenBack.computeHfsZeemanMatrix(settings, grid)
    outcomes = HfsPaschenBack.determineOutcomes(evals, evecs, basisList)

    HfsPaschenBack.displayResults(stdout, outcomes)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    HfsPaschenBack.displayResults(iostream, outcomes)   end

    if    output    return( outcomes )
    else            return( nothing )
    end
end


"""
`HfsPaschenBack.displayOutcomes(hfMultiplets::Array{Hfs.HfMultiplet,1})`
    ... to display the list of zero-field hyperfine levels admitted into the Paschen-Back perturbation matrix. A
        small neat table is printed but nothing is returned otherwise.
"""
function  displayOutcomes(hfMultiplets::Array{Hfs.HfMultiplet,1})
    nx = 43
    println(" ")
    println("  Hyperfine levels admitted into the HfsPaschenBack perturbation matrix:")
    println(" ")
    println("  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "HfLevel"; na=2);                          sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "F^P, M";  na=4);                          sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "Energy"; na=4);
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    println(sa);    println(sb);    println("  ", TableStrings.hLine(nx))
    #
    for  hfLevel in HfsPaschenBack.expandMSublevels(hfMultiplets)
        sa  = "  "
        sa = sa * TableStrings.center(10, string(hfLevel.index); na=2)
        sa = sa * TableStrings.center(10, "$(hfLevel.F)$(string(hfLevel.parity)), M=$(hfLevel.M)"; na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", hfLevel.energy)) * "    "
        println( sa )
    end
    println("  ", TableStrings.hLine(nx))
    #
    return( nothing )
end


"""
`HfsPaschenBack.displayResults(stream::IO, outcomes::Array{HfsPaschenBack.Outcome,1})`
    ... to display the field-dressed hyperfine Paschen-Back eigenenergies and their dominant zero-field HfLevel
        components. A neat table is printed but nothing is returned otherwise.
"""
function  displayResults(stream::IO, outcomes::Array{HfsPaschenBack.Outcome,1})
    nx = 100
    println(stream, " ")
    println(stream, "  Field-dressed HfsPaschenBack eigenstates:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(8,  "Eigen-#"; na=2);                           sb = sb * TableStrings.hBlank(10)
    sa = sa * TableStrings.center(16, "Energy";  na=4);
    sb = sb * TableStrings.center(16, TableStrings.inUnits("energy"); na=4)
    sa = sa * TableStrings.center(66, "Dominant zero-field (HfLevel, F^P, M) components [weight]"; na=4)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))
    #
    for  (k,outcome) in enumerate(outcomes)
        sa  = "  "
        sa = sa * TableStrings.center(8, string(k); na=2)
        sa = sa * @sprintf("%.10e", Defaults.convertUnits("energy: from atomic", outcome.energy)) * "    "
        parts = String[]
        for  comp in outcome.components
            push!(parts, "HF$(comp.hfLevel.index)(F=$(comp.hfLevel.F)$(string(comp.hfLevel.parity)),M=$(comp.hfLevel.M))" *
                         "[$(@sprintf("%.3f", comp.weight))]")
        end
        sa = sa * join(parts, "  ")
        println(stream, sa )
    end
    println(stream, "  ", TableStrings.hLine(nx))
    #
    return( nothing )
end

end # module
