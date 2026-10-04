"""
`module  JAC.NuclearShielding`  
... a submodel of JAC that computes how the electrons of a closed-shell ion MODIFY, at the nucleus, an electric or magnetic
    field that is applied from outside.  The question it answers is: *what does the nucleus actually feel?*

    This is the mirror image of `Hfs`.  There, a NUCLEAR moment couples to the electrons and one obtains the hyperfine
    constants A, B, C.  Here, an EXTERNAL field of the same multipolarity acts on the electrons, the electron cloud is
    distorted, and the distorted cloud produces a field of its own AT THE NUCLEUS.  The ratio of what is produced to what
    was applied is the shielding factor, and the two modules therefore carry the same multipole selectors.

    + the ELECTRIC QUADRUPOLE case (calcE2) returns gamma_inf, the quadrupole ANTISHIELDING factor.  For heavy ions it is
      large and NEGATIVE, i.e. the electrons amplify and invert the applied gradient instead of screening it, so that the
      gradient actually felt at the nucleus is (1 - gamma_inf) times the applied one.  This is the quantity that converts
      a lattice field gradient into an observed nuclear quadrupole splitting.
    + the MAGNETIC DIPOLE case (calcM1) would return sigma, the analogue of the NMR chemical shielding constant, which
      converts an applied magnetic field into the field felt at the nucleus and hence a measured NMR frequency into a
      nuclear magnetic moment.  IT IS NOT IMPLEMENTED YET and the settings refuse it rather than return a wrong number.

    WHAT THE NUMBERS MEAN, AND THE TWO CONDITIONS THEY REST ON.  A shielding factor is a LINEAR RESPONSE of the FREE ION:
    it is defined as the first-order reaction of the electron cloud to a field that is (i) weak and (ii) UNIFORM over the
    ion.  Condition (i) is always satisfied by crystal fields.  Condition (ii) is the subscript `inf`: it assumes the
    charges producing the field are infinitely far away.  In a real crystal the nearest neighbours sit only a little
    beyond the ion's own outer shells, and then the uniform-field value is NOT the right number -- which is why
    `NeighbourField` exists and why `displayResults` says so on every run.

    WHAT IS *NOT* COMPUTED HERE: any crystal or molecular effect.  No neighbour overlap beyond the electrostatics of the
    given charges, no charge transfer, no chemical bonding.  Those are the business of the application that calls this.
"""
module NuclearShielding


using  Printf, LinearAlgebra, ..AngularMomentum, ..Basics, ..Bsplines, ..Defaults, ..ManyElectron, ..Nuclear,
       ..Radial, ..RadialIntegrals, ..TableStrings


"""
`abstract type NuclearShielding.AbstractFieldModel`  
    ... defines an abstract type to specify WHAT is applied to the ion from outside, and hence what the returned factor
        means.  The choice is not a numerical parameter: it changes the physical content of the answer.

    + UniformField()          ... a strictly uniform external field gradient, i.e. the sources are infinitely far away.
                                  This returns the transferable, tabulated coefficient gamma_inf.
    + NeighbourField(charges) ... the exact rank-2 potential of point charges at finite distances, which falls off as
                                  R^2/r^3 outside each neighbour instead of growing as r^2.  This returns the EFFECTIVE
                                  factor for that particular arrangement, which is smaller in magnitude.
"""
abstract type  AbstractFieldModel                           end


"""
`struct  NuclearShielding.UniformField  <:  AbstractFieldModel`  
    ... the far-field limit; the perturbation grows as r^2 everywhere and the returned factor is the tabulated gamma_inf.
"""
struct   UniformField  <:  AbstractFieldModel                end


"""
`struct  NuclearShielding.NeighbourField  <:  AbstractFieldModel`  
    ... the potential of point charges at FINITE distance, whose rank-2 part is q * r_<^2 / r_>^3 with r_< and r_> the
        smaller and larger of (r, R).  Only the distances matter for the radial response, since the angular placement
        of the neighbours affects the LATTICE gradient and not the ion's response to it.

    + charges  ::Array{Tuple{Float64,Float64},1}  ... list of (charge [e], distance [a.u.]) of the neighbours.
"""
struct   NeighbourField  <:  AbstractFieldModel
    charges                    ::Array{Tuple{Float64,Float64},1}
end


"""
`struct  NuclearShielding.Contribution`  
    ... keeps what one occupied subshell contributes through one excitation channel, so that a total can be read rather
        than merely trusted.  This breakdown is the reason a sum-over-states formulation was chosen.

    + subshell     ::Subshell   ... the occupied subshell that is distorted.
    + kappaP       ::Int64      ... the kappa of the states it is distorted into.
    + value        ::Float64    ... its contribution to gamma_inf.
    + nStates      ::Int64      ... how many intermediate states entered the sum.
    + nOccupied    ::Int64      ... how many were excluded because they are occupied (Pauli).
"""
struct   Contribution
    subshell                   ::Subshell
    kappaP                     ::Int64
    value                      ::Float64
    nStates                    ::Int64
    nOccupied                  ::Int64
end


"""
`struct  NuclearShielding.Settings  <:  AbstractPropertySettings`  
    ... defines a type for the details and parameters of computing nuclear shielding factors.

    + calcE2        ::Bool               ... True if the electric-quadrupole antishielding factor gamma_inf is wanted.
    + calcM1        ::Bool               ... True if the magnetic-dipole shielding sigma is wanted; NOT YET IMPLEMENTED.
    + fieldModel    ::AbstractFieldModel ... what is applied from outside; see AbstractFieldModel.  This decides what
                                             the returned number MEANS and is echoed in the output for that reason.
    + nVirtualMax   ::Int64              ... maximum number of intermediate states per channel, ordered by |contribution|;
                                             0 keeps all of them.  Only for convergence studies -- the default keeps all.
    + selfConsistent::Bool               ... True if the induced field is allowed to act back on the electrons, i.e. if the
                                             response is iterated to convergence.  DEFAULT FALSE, so that the uncoupled
                                             number every earlier result was computed with stays what the module returns
                                             unless it is asked for otherwise.  See `computeInducedPotential` for what the
                                             feedback contains and, just as important, for what it does NOT.
    + scfIterations ::Int64              ... iterations of that loop.  The undamped iteration OSCILLATES, so this is used
                                             with `scfMixing` below; 16 is ample for the ions measured so far.
    + scfMixing     ::Float64            ... linear mixing of the new driving term into the old, in (0,1].  1.0 is no
                                             damping and oscillates; 0.5 converges to a spread of 0.002 in 16 steps.
    + printBefore   ::Bool               ... True if a list of selected levels is printed before the computations start.
    + levelSelection::LevelSelection     ... Specifies the selected levels, if any.
"""
struct Settings  <:  AbstractPropertySettings
    calcE2                     ::Bool
    calcM1                     ::Bool
    fieldModel                 ::AbstractFieldModel
    nVirtualMax                ::Int64
    selfConsistent             ::Bool
    scfIterations              ::Int64
    scfMixing                  ::Float64
    printBefore                ::Bool
    levelSelection             ::LevelSelection
end


"""
`NuclearShielding.Settings()`  ... constructor for an `empty` instance of NuclearShielding.Settings.
"""
function Settings()
    Settings(true, false, UniformField(), 0, false, 16, 0.5, false, LevelSelection() )
end


"""
`NuclearShielding.Settings(set::NuclearShielding.Settings;`

        calcE2=.., calcM1=.., fieldModel=.., nVirtualMax=.., selfConsistent=.., scfIterations=..,
        scfMixing=.., printBefore=.., levelSelection=..)

    ... keyword copy-constructor for re-defining selected values of a settings::NuclearShielding.Settings.
"""
function Settings(set::NuclearShielding.Settings;
        calcE2::Union{Nothing,Bool}=nothing,                calcM1::Union{Nothing,Bool}=nothing,
        fieldModel::Union{Nothing,AbstractFieldModel}=nothing,  nVirtualMax::Union{Nothing,Int64}=nothing,
        selfConsistent::Union{Nothing,Bool}=nothing,        scfIterations::Union{Nothing,Int64}=nothing,
        scfMixing::Union{Nothing,Float64}=nothing,
        printBefore::Union{Nothing,Bool}=nothing,           levelSelection::Union{Nothing,LevelSelection}=nothing)
    if  isnothing(calcE2)           calcE2x         = set.calcE2         else   calcE2x         = calcE2         end
    if  isnothing(calcM1)           calcM1x         = set.calcM1         else   calcM1x         = calcM1         end
    if  isnothing(fieldModel)       fieldModelx     = set.fieldModel     else   fieldModelx     = fieldModel     end
    if  isnothing(nVirtualMax)      nVirtualMaxx    = set.nVirtualMax    else   nVirtualMaxx    = nVirtualMax    end
    if  isnothing(selfConsistent)   selfConsistentx = set.selfConsistent else   selfConsistentx = selfConsistent end
    if  isnothing(scfIterations)    scfIterationsx  = set.scfIterations  else   scfIterationsx  = scfIterations  end
    if  isnothing(scfMixing)        scfMixingx      = set.scfMixing      else   scfMixingx      = scfMixing      end
    if  isnothing(printBefore)      printBeforex    = set.printBefore    else   printBeforex    = printBefore    end
    if  isnothing(levelSelection)   levelSelectionx = set.levelSelection else   levelSelectionx = levelSelection end

    Settings( calcE2x, calcM1x, fieldModelx, nVirtualMaxx, selfConsistentx, scfIterationsx, scfMixingx,
              printBeforex, levelSelectionx )
end


# `Base.show(io::IO, settings::NuclearShielding.Settings)`  ... prepares a proper printout of settings::NuclearShielding.Settings.
function Base.show(io::IO, settings::NuclearShielding.Settings)
    println(io, "calcE2:                   $(settings.calcE2)  ")
    println(io, "calcM1:                   $(settings.calcM1)  ")
    println(io, "fieldModel:               $(settings.fieldModel)  ")
    println(io, "nVirtualMax:              $(settings.nVirtualMax)  ")
    println(io, "selfConsistent:           $(settings.selfConsistent)  ")
    println(io, "scfIterations:            $(settings.scfIterations)  ")
    println(io, "scfMixing:                $(settings.scfMixing)  ")
    println(io, "printBefore:              $(settings.printBefore)  ")
    println(io, "levelSelection:           $(settings.levelSelection)  ")
end


"""
`struct  NuclearShielding.Outcome`  
    ... defines a type to keep the outcome of a nuclear-shielding computation.

    + level          ::Level                       ... Atomic level to which the outcome refers to.
    + gammaE2        ::Float64                     ... the quadrupole antishielding factor gamma_inf.
    + sigmaM1        ::Float64                     ... the magnetic-dipole shielding sigma; 0. while unimplemented.
    + fieldModel     ::AbstractFieldModel          ... what was applied, kept WITH the number so that a later script
                                                       cannot separate a value from its meaning.
    + contributions  ::Array{Contribution,1}       ... the per-subshell, per-channel breakdown.
"""
struct Outcome
    level                      ::Level
    gammaE2                    ::Float64
    sigmaM1                    ::Float64
    fieldModel                 ::AbstractFieldModel
    contributions              ::Array{Contribution,1}
end


"""
`NuclearShielding.Outcome()`  ... constructor for an `empty` instance of NuclearShielding.Outcome.
"""
function Outcome()
    Outcome( Level(), 0., 0., UniformField(), Contribution[] )
end


# `Base.show(io::IO, outcome::NuclearShielding.Outcome)`  ... prepares a proper printout of outcome::NuclearShielding.Outcome.
function Base.show(io::IO, outcome::NuclearShielding.Outcome)
    println(io, "level:                    $(outcome.level)  ")
    println(io, "gammaE2:                  $(outcome.gammaE2)  ")
    println(io, "sigmaM1:                  $(outcome.sigmaM1)  ")
    println(io, "fieldModel:               $(outcome.fieldModel)  ")
    println(io, "contributions:            $(length(outcome.contributions)) entries  ")
end


"""
`NuclearShielding.drivingFunction(model::NuclearShielding.UniformField, r::Float64)`  
    ... returns the radial factor f(r) = r^2 of a strictly uniform external field gradient, whose rank-2 potential grows
        as r^2 everywhere; a value::Float64 is returned.
"""
function drivingFunction(model::NuclearShielding.UniformField, r::Float64)
    return( r*r )
end


"""
`NuclearShielding.drivingFunction(model::NuclearShielding.NeighbourField, r::Float64)`  
    ... returns the radial factor f(r) = sum_k q_k r_<^2 / r_>^3 of the EXACT rank-2 potential of point charges at finite
        distance, the standard r_<^L / r_>^(L+1) kernel with one coordinate frozen at the neighbour's distance.  It grows
        as r^2 only INSIDE the nearest neighbour and falls as 1/r^3 outside it, which is the whole difference between an
        effective and a far-field shielding factor; a value::Float64 is returned.
"""
function drivingFunction(model::NuclearShielding.NeighbourField, r::Float64)
    wa = 0.
    for  (q, R)  in  model.charges
        rs = min(r, R);   rg = max(r, R);   wa = wa + q * rs*rs / (rg*rg*rg)
    end
    return( wa )
end


"""
`NuclearShielding.appliedGradient(model::NuclearShielding.UniformField)`  
    ... returns the field gradient V_zz that the driving function of this model produces at the nucleus, so that the
        computed response can be divided by it; a value::Float64 is returned.  For f(r) = r^2 it is 2.
"""
function appliedGradient(model::NuclearShielding.UniformField)
    return( 2.0 )
end


"""
`NuclearShielding.appliedGradient(model::NuclearShielding.NeighbourField)`  
    ... returns the field gradient V_zz produced at the nucleus by the given point charges, sum_k 2 q_k / R_k^3, which is
        the quantity a lattice sum reports; a value::Float64 is returned.
"""
function appliedGradient(model::NuclearShielding.NeighbourField)
    wa = 0.
    for  (q, R)  in  model.charges      wa = wa + 2.0 * q / (R*R*R)      end
    return( wa )
end


"""
`NuclearShielding.allowedKappas(kappa::Int64)`  
    ... returns every kappa' that a RANK-2 operator can reach from the given kappa, i.e. those for which the reduced
        matrix element <kappa || C^(2) || kappa'> does not vanish by the triangle and parity rules; this is the list of
        channels through which a shell can be distorted.  An array::Array{Int64,1} is returned.
"""
function allowedKappas(kappa::Int64)
    kappas = Int64[];    sa = Subshell(9, kappa)
    for  kp = -9:9
        kp == 0  &&  continue
        abs(AngularMomentum.CL_reduced_me(sa, 2, Subshell(9, kp))) > 1.0e-12   &&   push!(kappas, kp)
    end
    return( kappas )
end


"""
`NuclearShielding.projectOperator(f::Function, orb::Radial.Orbital, primitives::Bsplines.Primitives, grid::Radial.Grid)`  
    ... returns the vector <B_i | f(r) | orb> of the given radial operator acting on an orbital and projected onto the
        B-spline primitives, large and small component alike, since f(r) is multiplicative and acts on both the same way.
        An array::Array{Float64,1} of length nsL+nsS is returned.
"""
function projectOperator(f::Function, orb::Radial.Orbital, primitives::Bsplines.Primitives, grid::Radial.Grid)
    rhs = zeros(grid.nsL + grid.nsS)
    for  (offset, splines, comp)  in  ((0, primitives.bsplinesL, orb.P), (grid.nsL, primitives.bsplinesS, orb.Q))
        for  i = 1:length(splines)
            Bi = splines[i];    add = 1 - Bi.lower;    wa = 0.
            for  r = max(2, Bi.lower):min(Bi.upper, length(comp))
                wa = wa + Bi.bs[r+add] * f(grid.r[r]) * comp[r] * grid.wr[r]
            end
            rhs[offset+i] = wa
        end
    end
    return( rhs )
end


"""
`NuclearShielding.computeChannelResponse(orb::Radial.Orbital, kappaP::Int64, bDrive::Array{Float64,1}, bObserve::Array{Float64,1},`
                                        ` pot::Radial.Potential, primitives::Bsplines.Primitives, occupiedEnergies::Array{Float64,1},`
                                        ` storage::Dict{String,Array{Float64,2}}, overlap::Array{Float64,2}, nVirtualMax::Int64)`  
    ... performs the sum over intermediate states of one excitation channel: the occupied orbital is pushed by the external
        field into the states of symmetry kappaP, and each of them produces a field gradient back at the nucleus.  States
        that are themselves OCCUPIED are excluded, since a rotation among occupied orbitals is not a physical excitation.
        A tuple tpl(value::Float64, nStates::Int64, nOccupied::Int64, dvector::Array{Float64,1}) is returned, the last
        being the B-spline expansion of the PERTURBED orbital of this channel, sum_v |v> <v|f_drive|a>/(eps_a - eps_v).
        It is what the self-consistent loop needs and it costs nothing extra to accumulate, since every term of it is
        already formed for the scalar.
"""
function computeChannelResponse(orb::Radial.Orbital, kappaP::Int64, bDrive::Array{Float64,1}, bObserve::Array{Float64,1},
                                pot::Radial.Potential, primitives::Bsplines.Primitives, occupiedEnergies::Array{Float64,1},
                                storage::Dict{String,Array{Float64,2}}, overlap::Array{Float64,2}, nVirtualMax::Int64)
    matrixA = Bsplines.setupLocalMatrix(kappaP, primitives, pot, storage)
    eigen   = Bsplines.diagonalizeLocalMatrix(kappaP, matrixA, overlap, primitives)
    epsA    = orb.energy;    terms = Float64[];    nOccupied = 0
    dvector = zeros(length(bDrive))

    for  (i, vector)  in  enumerate(eigen.vectors)
        epsV = eigen.values[i]
        # exclude the occupied states of this symmetry, and with them the vanishing denominator
        isOcc = false
        for  eOcc  in  occupiedEnergies
            abs(epsV - eOcc) < 1.0e-6 * max(1.0, abs(eOcc))   &&   (isOcc = true;  break)
        end
        if  isOcc    nOccupied = nOccupied + 1;   continue    end
        abs(epsA - epsV) < 1.0e-10   &&   continue
        coeff = dot(vector, bDrive) / (epsA - epsV)
        push!(terms, coeff * dot(vector, bObserve))
        dvector = dvector + coeff * vector
    end

    if  nVirtualMax > 0  &&  length(terms) > nVirtualMax
        terms = sort(terms, by = t -> -abs(t))[1:nVirtualMax]
    end

    return( (sum(terms), length(terms), nOccupied, dvector) )
end


"""
`NuclearShielding.GAMMA_PREFACTOR`
    ... the single overall constant of the quadrupole response, -2/5.

        IT IS DERIVED, SIGN AND ALL -- four factors, each a line of electrostatics or angular-momentum algebra,
        and nothing in it is fitted (derived 04-Oct-2026; until then the magnitude's factor of two and the sign
        had been calibrated against the validation set, and this docstring said so):

            2            LINEAR RESPONSE.  The induced gradient is <0|Q|v><v|H'|0> + c.c., and the two terms are
                         equal for real orbitals and Hermitian operators.
            x 2          THE OBSERVABLE.  The field gradient at the nucleus from one electron is 2 P_2(cos t)/r^3,
                         the 2 being d^2/dz^2 of r^2 P_2.  `fObserve` in `computeAmplitudesProperties` carries
                         only the 1/r^3, SO THAT 2 LIVES HERE.  That is what made this constant look arbitrary:
                         the factor was real and had simply been moved.
            x 1/5        THE ANGULAR REDUCTION.  Summed over the magnetic substates of a closed subshell,
                         sum_m |(j_a 2 j_v; -m_a 0 m_v)|^2 = 1/(2k+1) = 1/5 at rank k = 2, which is what turns
                         the reduced matrix elements into the response of the whole shell.
            / 2          THE APPLIED GRADIENT.  `drivingFunction` uses f(r) = r^2, whose own gradient at the
                         nucleus is `appliedGradient` = 2, and the ratio is divided by it to make a shielding
                         FACTOR rather than a field.
            = 2/5,  and the MINUS from  gamma_inf = -q_induced/q_applied,  the two electron charges (-1 in the
                    perturbation and -1 in the observable) having already cancelled against each other.

        THE SIGN NEVER NEEDED STERNHEIMER'S PAPER.  It is fixed internally, by the relation between
        `drivingFunction` and `appliedGradient`, both of which are in this module; the earlier note that it
        "depends on the sign convention of his driving potential, which reaches us as a poor scan" was looking
        for the answer in the wrong place.

        AND THE ANGULAR CONTENT IS A CHECK RATHER THAN THE DEFINITION.  4/5 times the sum of
        |<kappa_a||C^(2)||kappa'>|^2 over the spin-orbit partners reproduces Sternheimer's non-relativistic
        closed-shell coefficients 48/25 (p), 16/7 (d) and 224/75 (f) with ratio 1.0000000 in all three cases --
        the sums themselves being 12/5, 20/7 and 56/15.  Note that 4/5 is TWICE the constant above, because
        Sternheimer's coefficients are quoted for a radial integral defined with the other half.

        `example-Cp.jl` branch (f) checks every one of these links, and needs no orbitals to do it.
        Do not re-tune this constant to improve an individual ion.
"""
const GAMMA_PREFACTOR = -0.4


"""
`NuclearShielding.multipolePotential(density::Array{Float64,1}, grid::Radial.Grid)`
    ... returns the rank-2 potential of a radial density, Y_2[f](r) = (1/r^3) int_0^r r'^2 f dr' + r^2 int_r^inf f/r'^3 dr',
        tabulated on the grid; an array::Array{Float64,1} is returned.

        THE SECOND TERM IS WHY THIS FUNCTION IS WORTH ITS OWN NAME.  Its coefficient, int_0^inf f/r'^3 dr', is the
        gradient-carrying part: as r -> 0 the first term goes to a constant and only the second behaves as r^2.  That is
        what lets the induced potential be checked against the gradient the module computes directly -- see
        `computeInducedPotential`.
"""
function multipolePotential(density::Array{Float64,1}, grid::Radial.Grid)
    n = min(length(density), length(grid.r));    inner = zeros(n);    outer = zeros(n);    acc = 0.
    for  i = 2:n     acc = acc + grid.r[i]^2 * density[i] * grid.wr[i];    inner[i] = acc    end
    acc = 0.
    for  i = n:-1:2  acc = acc + density[i] / grid.r[i]^3 * grid.wr[i];    outer[i] = acc    end
    potential = zeros(length(grid.r))
    for  i = 2:n     potential[i] = inner[i]/grid.r[i]^3 + grid.r[i]^2 * outer[i]    end

    return( potential )
end


"""
`NuclearShielding.computeInducedPotential(density::Array{Float64,1}, rhoTotal::Array{Float64,1}, grid::Radial.Grid)`
    ... returns the radial factor of the rank-2 potential that the induced electron density makes, i.e. the term that
        turns a first-order response into a self-consistent one; an array::Array{Float64,1} is returned.

        TWO PIECES, AND THE SECOND IS OFTEN LEFT OUT.
        (i)  THE DIRECT (HARTREE) TERM, `INDUCED_SCALE` times the rank-2 potential of the induced density.  The scale is
             not fitted: the gradient that Y_2[density] makes at the nucleus was measured against the gradient the module
             computes directly, and the ratio came out as 5/2 EXACTLY on three ions spanning Z = 39 to 90 and gamma from
             -27 to -185.  A constant rational across that range is bookkeeping, not physics -- it is the same 1/5 and 2
             that `GAMMA_PREFACTOR` carries, met once more on the way back from density to potential.
        (ii) THE EXCHANGE RESPONSE, the local derivative of the SAME Slater term the DFS mean field is built with:
             V_x = -(3 rho_t/(4 pi^2 r^2))^(1/3) goes as rho_t^(1/3), so dV_x = (1/3)(V_x/rho_t) drho_2.  No new
             functional and no new parameter.  Measured, it is a correction of about 1.5 percentage points on the
             self-consistent shift and it does NOT change its sign anywhere.

        WHAT IS STILL MISSING, stated because it bounds what the self-consistent number means: the exchange response is
        LOCAL, as the mean field it derives from is; a non-local (Dirac-Fock) exchange response is a different and larger
        calculation.  And the response remains that of a free, spherical, closed-shell ion.
"""
function computeInducedPotential(density::Array{Float64,1}, rhoTotal::Array{Float64,1}, grid::Radial.Grid)
    n       = length(grid.r)
    induced = INDUCED_SCALE .* multipolePotential(density, grid)
    drho2   = INDUCED_SCALE .* density
    for  i = 2:min(n, length(rhoTotal), length(drho2))
        rhoTotal[i] <= 1.0e-30   &&   continue
        vX = -(3*rhoTotal[i]/(4*pi^2*grid.r[i]^2))^(1/3)
        induced[i] = induced[i] + (1/3) * vX * drho2[i] / rhoTotal[i]
    end

    return( induced )
end


"""
`NuclearShielding.INDUCED_SCALE`
    ... the one constant of the feedback, 2/5, and it is measured rather than chosen -- see `computeInducedPotential`.
        It is the same 2/5 as `GAMMA_PREFACTOR` because it is the same chain read backwards.
"""
const INDUCED_SCALE = 0.4


"""
`NuclearShielding.computeAmplitudesProperties(outcome::NuclearShielding.Outcome, nm::Nuclear.Model, grid::Radial.Grid,`
                                             ` settings::NuclearShielding.Settings)`  
    ... computes the shielding factor of the given level by summing, over every occupied subshell and every rank-2
        excitation channel, the response of the electron cloud to the external field.  The angular weight of a channel is
        |<kappa_a || C^(2) || kappa'>|^2 and the overall constant is `GAMMA_PREFACTOR`, which is DERIVED -- see its own
        docstring for the four factors -- and checked against Sternheimer's closed-shell coefficients 48/25 (p->p),
        16/7 (d->d) and 224/75 (f->f).  NOTE that `fObserve` below is 1/r^3 and NOT the full 2/r^3 of the gradient
        operator; that factor of two sits in `GAMMA_PREFACTOR`.
        A new outcome::NuclearShielding.Outcome is returned.
"""
function computeAmplitudesProperties(outcome::NuclearShielding.Outcome, nm::Nuclear.Model, grid::Radial.Grid,
                                     settings::NuclearShielding.Settings)
    settings.calcM1   &&   error("\n\nNuclearShielding: calcM1 = true, but the MAGNETIC-DIPOLE shielding sigma is NOT yet "  *
                                 "implemented.\n"                                                                            *
                                 ">>> The rank-1 machinery exists in JAC -- InteractionStrength.zeeman_n1 for the applied "   *
                                 "field and hfs_tM1 for the field at the nucleus -- but the response has not been wired up "  *
                                 "or validated, and a number returned now would be wrong rather than approximate.\n"          *
                                 ">>> Set calcM1 = false.\n")
    settings.calcE2   ||   return( outcome )

    level = outcome.level;    basis = level.basis
    # the mean field the electrons move in, and the B-spline primitives for the intermediate states
    pot     = Basics.add( Nuclear.nuclearPotential(nm, grid), Basics.computePotential(Basics.DFSField(1.0), grid, basis) )
    prims   = Bsplines.generatePrimitives(grid)
    storage = Dict{String,Array{Float64,2}}()
    nsL     = grid.nsL;    nsS = grid.nsS
    overlap = zeros(nsL+nsS, nsL+nsS)
    overlap[1:nsL, 1:nsL]                 = Bsplines.generateTTpMatrix!("LL-overlap", 0, prims, storage)
    overlap[nsL+1:nsL+nsS, nsL+1:nsL+nsS] = Bsplines.generateTTpMatrix!("SS-overlap", 0, prims, storage)

    # the occupied subshells, and their energies gathered per kappa so that Pauli-blocked states can be skipped
    occupied = Subshell[];    occEnergies = Dict{Int64,Array{Float64,1}}()
    for  sh  in  basis.subshells
        haskey(basis.orbitals, sh)  ||  continue
        Basics.computeMeanSubshellOccupation(sh, [level]) < 1.0e-6   &&   continue
        push!(occupied, sh)
        push!( get!(occEnergies, sh.kappa, Float64[]), basis.orbitals[sh].energy )
    end

    model = settings.fieldModel;    applied = appliedGradient(model)
    abs(applied) < 1.0e-30   &&   error("\n\nNuclearShielding: the chosen fieldModel applies NO gradient at the nucleus, "  *
                                        "so a shielding RATIO is undefined.  Check the charges.\n")
    nr          = length(grid.r)
    fObserve(r) = 1.0 / (r*r*r)
    fExt        = [ drivingFunction(model, grid.r[i]) for i = 1:nr ]
    fEff        = copy(fExt)

    # the total density of the ion, needed only by the exchange response of the self-consistent loop
    rhoTotal = zeros(nr)
    if  settings.selfConsistent
        for  sh  in  occupied
            orb = basis.orbitals[sh];    occ = Basics.computeMeanSubshellOccupation(sh, [level])
            for  i = 1:min(length(orb.P), nr)    rhoTotal[i] += occ * (orb.P[i]^2 + orb.Q[i]^2)    end
        end
    end

    # ONE PASS IF UNCOUPLED, AND THE SAME PASS REPEATED IF NOT.  The driving term is the only thing that changes
    # between iterations, which is why it sits behind `drivingFunction` in the first place; the B-spline spectrum of
    # each symmetry does not depend on it and is rebuilt per pass only because `computeChannelResponse` owns it.
    nPass         = settings.selfConsistent ? max(1, settings.scfIterations) : 1
    contributions = Contribution[];    gamma = 0.
    for  pass = 1:nPass
        contributions = Contribution[];    gamma = 0.;    density = zeros(nr)
        for  sh  in  occupied
            orb      = basis.orbitals[sh]
            bDrive   = projectOperator(r -> fEff[max(2, searchsortedfirst(grid.r, r))], orb, prims, grid)
            bObserve = projectOperator(fObserve, orb, prims, grid)
            for  kappaP  in  allowedKappas(sh.kappa)
                weight = AngularMomentum.CL_reduced_me(sh, 2, Subshell(9, kappaP))^2
                occE   = get(occEnergies, kappaP, Float64[])
                value, nStates, nOcc, dvec = computeChannelResponse(orb, kappaP, bDrive, bObserve, pot, prims, occE,
                                                                    storage, overlap, settings.nVirtualMax)
                wa     = GAMMA_PREFACTOR * (2.0/applied) * weight * value
                gamma  = gamma + wa
                push!(contributions, Contribution(sh, kappaP, wa, nStates, nOcc))
                if  settings.selfConsistent
                    dOrb = Bsplines.generateOrbitalFromPrimitives(Subshell(9, kappaP), 0.0, nr, dvec, prims)
                    for  i = 2:min(length(orb.P), length(dOrb.P), nr)
                        density[i] += weight * (orb.P[i]*dOrb.P[i] + orb.Q[i]*dOrb.Q[i])
                    end
                end
            end
        end
        settings.selfConsistent   ||   break
        induced = computeInducedPotential(density, rhoTotal, grid)
        fEff    = (1 - settings.scfMixing) .* fEff .+ settings.scfMixing .* (fExt .+ induced)
    end

    return( Outcome(level, gamma, 0., model, contributions) )
end


"""
`NuclearShielding.determineOutcomes(multiplet::Multiplet, settings::NuclearShielding.Settings)`  
    ... determines which levels of the given multiplet should be considered; an array of outcomes::Array{Outcome,1} is
        returned, one for each selected level, with the shielding factors still unset.
"""
function  determineOutcomes(multiplet::Multiplet, settings::NuclearShielding.Settings)
    outcomes = NuclearShielding.Outcome[]
    for  level  in  multiplet.levels
        if  Basics.selectLevel(level, settings.levelSelection)
            push!( outcomes, NuclearShielding.Outcome(level, 0., 0., settings.fieldModel, Contribution[]) )
        end
    end
    return( outcomes )
end


"""
`NuclearShielding.displayOutcomes(outcomes::Array{NuclearShielding.Outcome,1})`  
    ... lists the selected levels before the computation starts, so that a long run can be interrupted early if the wrong
        levels were picked; nothing is returned.
"""
function  displayOutcomes(outcomes::Array{NuclearShielding.Outcome,1})
    nx = 43
    println(" ");    println("  Selected levels for nuclear-shielding computations:");    println(" ")
    println("  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2);                             sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "J^P";   na=4);                             sb = sb * TableStrings.hBlank(14)
    sa = sa * TableStrings.center(14, "Energy"; na=4)
    sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=4)
    println(sa);    println(sb);    println("  ", TableStrings.hLine(nx))
    for  outcome  in  outcomes
        sa  = "  ";    sym = LevelSymmetry( outcome.level.J, outcome.level.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(outcome.level.index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=4)
        sa = sa * @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", outcome.level.energy)) * "    "
        println( sa )
    end
    println("  ", TableStrings.hLine(nx))

    return( nothing )
end


"""
`NuclearShielding.displayResults(stream::IO, outcomes::Array{NuclearShielding.Outcome,1}, settings::NuclearShielding.Settings)`  
    ... prints the shielding factors together with the STATEMENT OF WHAT THEY MEAN, which is not decoration: a shielding
        factor is a conditional quantity and a number quoted without its field model and its free-ion caveat invites
        exactly the misuse this module exists to prevent.  Both gamma_inf and (1 - gamma_inf) are printed because the
        published sign conventions for this quantity disagree with one another.  Nothing is returned.
"""
function  displayResults(stream::IO, outcomes::Array{NuclearShielding.Outcome,1}, settings::NuclearShielding.Settings)
    nx = 76
    println(stream, " ");    println(stream, "  Nuclear shielding factors:");    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  "
    sa = sa * TableStrings.center(10, "Level"; na=2)
    sa = sa * TableStrings.center(10, "J^P";   na=3)
    sa = sa * TableStrings.center(18, "gamma_inf";     na=2)
    sa = sa * TableStrings.center(18, "1 - gamma_inf"; na=2)
    println(stream, sa);    println(stream, "  ", TableStrings.hLine(nx))
    for  outcome  in  outcomes
        sa  = "  ";    sym = LevelSymmetry( outcome.level.J, outcome.level.parity)
        sa = sa * TableStrings.center(10, TableStrings.level(outcome.level.index); na=2)
        sa = sa * TableStrings.center(10, string(sym); na=3)
        sa = sa * @sprintf("%16.4f", outcome.gammaE2) * "    "
        sa = sa * @sprintf("%16.4f", 1.0 - outcome.gammaE2)
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))

    # --- what the number means; printed with it, never only in a docstring
    model = settings.fieldModel
    println(stream, "")
    if      typeof(model) == NuclearShielding.UniformField
        println(stream, "  field model   UniformField   ... the external charges are assumed FAR compared with the ion.")
        println(stream, "                This is the transferable, free-ion coefficient: it multiplies a LATTICE field")
        println(stream, "                gradient to give the gradient felt at the nucleus, in ANY host.")
        println(stream, "                For neighbours closer than about twice the ion radius it OVERESTIMATES the")
        println(stream, "                response; use fieldModel = NeighbourField(...) and compare.")
    else
        println(stream, "  field model   NeighbourField ... the EXACT rank-2 potential of the given charges at finite")
        println(stream, "                distance.  The number above is the EFFECTIVE factor for that arrangement and is")
        println(stream, "                NOT transferable to another host; divide it by the UniformField value to obtain")
        println(stream, "                the penetration correction.")
    end
    println(stream, "  response      uncoupled      ... the induced field is NOT yet allowed to act back on the")
    println(stream, "                electrons (no RPA/self-consistency).  That correction is expected to be of order")
    println(stream, "                10-30 % and to REDUCE the magnitude; see the module docstring.")
    println(stream, "  THIS IS A FREE-ION NUMBER: no neighbour overlap beyond the given point charges, no charge")
    println(stream, "  transfer, no chemical bonding.")

    # --- the breakdown, which is why a sum-over-states form was chosen
    for  outcome  in  outcomes
        length(outcome.contributions) == 0   &&   continue
        println(stream, "\n  Breakdown for level $(outcome.level.index), largest channels first:")
        println(stream, "  ", TableStrings.hLine(62))
        @printf(stream, "  %-12s %8s %16s %10s %10s\n", "subshell", "kappa'", "contribution", "states", "blocked")
        println(stream, "  ", TableStrings.hLine(62))
        for  c  in  first(sort(outcome.contributions, by = x -> -abs(x.value)), 12)
            @printf(stream, "  %-12s %8d %16.5f %10d %10d\n", string(c.subshell), c.kappaP, c.value, c.nStates, c.nOccupied)
        end
        println(stream, "  ", TableStrings.hLine(62))
    end

    return( nothing )
end


"""
`NuclearShielding.computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::NuclearShielding.Settings; output=true)`  
    ... computes the nuclear shielding factors for the selected levels of the given multiplet.  The ion must be CLOSED
        SHELL: the response of an open-shell ion is not a scalar and the factorisation this module rests on does not
        apply, so such a case is refused rather than approximated.  An array of outcomes::Array{Outcome,1} is returned if
        output = true, and nothing otherwise.
"""
function computeOutcomes(multiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::NuclearShielding.Settings;
                         output=true)
    println("")
    printstyled("NuclearShielding.computeOutcomes(): The computation of the shielding factors starts now ... \n", color=:light_green)
    printstyled("--------------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    outcomes = NuclearShielding.determineOutcomes(multiplet, settings)
    if  settings.printBefore    NuclearShielding.displayOutcomes(outcomes)    end

    newOutcomes = NuclearShielding.Outcome[]
    for  outcome  in  outcomes
        if  Basics.twice(outcome.level.J) != 0   ||   length(outcome.level.basis.csfs) != 1
            error("\n\nNuclearShielding: the level $(outcome.level.index) is not a CLOSED-SHELL, J = 0 state "  *
                  "(2J = $(Basics.twice(outcome.level.J)), $(length(outcome.level.basis.csfs)) CSF).\n"          *
                  ">>> A shielding factor is the SCALAR response of a spherical ion.  For an open-shell ion the "*
                  "response is not a single number and this module does not apply.\n")
        end
        push!( newOutcomes, NuclearShielding.computeAmplitudesProperties(outcome, nm, grid, settings) )
    end

    NuclearShielding.displayResults(stdout, newOutcomes, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    NuclearShielding.displayResults(iostream, newOutcomes, settings)    end

    if    output    return( newOutcomes )
    else            return( nothing )
    end
end

end # module
