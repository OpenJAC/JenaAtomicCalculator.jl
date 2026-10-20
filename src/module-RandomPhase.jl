
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
        RESULT.

        **THAT NEXT STEP HAS NOW BEEN TAKEN, IT FOUND A DEFECT, AND THE DEFECT IS REPAIRED -- 14-Oct-2026.**
        The genuine multi-electron SCF system the paragraph above asks for is Ne 1s^2 2s^2 2p^6 in its own
        DFS potential, E1 at omega = 0.1, three active subshells. The coupled iteration diverged there under
        any damping, and measuring the iteration map's spectral radius directly -- by scaling the induced term
        by a factor lambda and tracking successive residual norms -- gave rho = 5.03, where a closed shell's
        tens-of-percent screening calls for 0.2-0.4. The cause was that THIS MODULE APPLIED NO ANGULAR
        COEFFICIENT AT ALL, while a rank-L induced potential needs a reduced matrix element of C^L at each
        end, and weighted each electron by (2j_b+1) where the reduced convention wants 2*sqrt(2j_b+1). Both
        are now fixed, the radial route is verified against `InteractionStrength.XL_Coulomb` to ten digits
        including its sign, and the response is confirmed to SCREEN rather than anti-screen.

        **AND WHY THE COUPLED PATH STILL CANNOT BE ITERATED -- THE ANSWER, NOT A SUSPICION.** The angular
        repair improved the spectral radius 2.17x, from 5.03 to 2.32, and left three candidates for the rest.
        All three were then measured and two are eliminated: widening the channel set to every symmetry the
        multipole reaches (`Settings.allChannels`, 3 channels to 7 on this case) does NOT help, 2.33 against
        2.49; and the `2.0` density symmetrization is not the cause either, since sweeping omega down to
        0.001 -- where that factor becomes EXACT -- leaves rho RISING to 2.72 rather than falling. What
        remains is the approximation itself: **only the DIRECT (ring) term is built, a direct-only induced
        field OVER-screens because exchange opposes it, and rho ~ 2.7 is what that looks like.**

        **SO THE POSITION IS NOW CLEAR, AND IT IS A QUESTION OF SCOPE.** rho ~ 2.7 is a property of the
        operator rather than of the stepping, so no damping, acceleration or channel choice can reach the
        fixed point -- a DIRECT SOLVE could, since (I - J) stays invertible, and that is the right method for
        it. But solving this term exactly would make the module CONVERGE without making it RIGHT: a
        direct-only screening correction for a closed shell is not trustworthy, and exchange is not a small
        correction in this regime. Whether to build the direct solve anyway, to extend to exchange (which the
        scope note below puts out of reach), or to leave `couplingEnabled = true` as a documented
        demonstration is the maintainer's decision. The single-electron limit remains exact throughout and is
        unaffected by any of this. The full account, with every number, is at `RandomPhase.solve`.

        WHAT IS STILL OUT OF SCOPE EVEN NOW. Only the DIRECT term is built: the EXCHANGE contribution to the
        induced field needs the same two-electron angular machinery (`SpinAngular`) the CI path already
        carries, a substantially larger undertaking -- "Random Phase Approximation" (no exchange) is the
        right target for this module, not "RPA with Exchange".
"""
module RandomPhase


using  LinearAlgebra, Printf, ..AngularMomentum, ..Basics, ..Bsplines, ..DalgarnoLewis, ..ManyElectron, ..Radial,
       ..RadialIntegrals, JenaAtomicCalculator


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
    + andersonDepth       ::Int64                     ... depth of the Anderson history used to accelerate the
                                                         coupled iteration.  **THE DEFAULT IS 0, i.e. plain
                                                         linear mixing, and that is a MEASURED choice rather
                                                         than a cautious one:** Anderson was implemented for
                                                         priority item 62 and does NOT repair this fixed point
                                                         -- on the Ne case in `solve`'s own note it is WORSE
                                                         than heavy damping at every depth tried.  It is kept
                                                         because it is sound, costs nothing when unused, and
                                                         may help a better-posed case;  do not switch it on
                                                         without checking that it converges.
    + directSolve         ::Bool                      ... True if the coupled equations are SOLVED as one linear
                                                         system, (I - M) s = b on the concatenated channel
                                                         sources, instead of being iterated.  **This is the
                                                         only route that can reach the fixed point at all**:
                                                         the iteration's spectral radius was measured at
                                                         about 2.7 on a closed shell (14-Oct-2026), so no
                                                         damping or acceleration converges, while (I - M)
                                                         stays invertible.  Implemented with
                                                         `RandomPhase.gmres` over `RandomPhase.sourceMap`,
                                                         so no matrix is ever formed.  `maxIterations` then
                                                         bounds the Krylov dimension and `accuracy` is the
                                                         relative residual demanded.  **THE DEFAULT IS TRUE
                                                         since 14-Oct-2026**, because an iterating default
                                                         cannot converge at all;  and on this route
                                                         `nIterations` of the returned `CoupledResponse`
                                                         counts KRYLOV STEPS rather than sweeps.
    + allChannels         ::Bool                      ... True if the induced field is built from EVERY target
                                                         symmetry the multipole can reach from each active
                                                         subshell, rather than from the single one named in
                                                         `targetSubshells`.  **This is what an RPA actually
                                                         asks for**, and the one-target truncation was measured
                                                         on 14-Oct-2026 to change the iteration's spectral
                                                         radius between 1.94 and 3.66 on nothing but the choice
                                                         of target -- because the two spin-orbit partners of a
                                                         subshell carry OPPOSITE angular signs and cancel, so a
                                                         truncated set breaks a cancellation rather than merely
                                                         dropping a term.  `targetSubshells` is still required
                                                         and still names each active subshell's PRIMARY channel,
                                                         which is the one whose chi is returned;  the extra
                                                         channels are solved too and enter the induced density.
                                                         **THE DEFAULT IS TRUE since 14-Oct-2026.** It was
                                                         false until then, which made the default coupled path
                                                         the one measured to give a sign-reversed total; see
                                                         `RandomPhase.solve`'s note for the numbers.
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
    andersonDepth            ::Int64
    allChannels              ::Bool
    directSolve              ::Bool
end


"""
`RandomPhase.Settings()`  ... constructor for an `empty` instance of RandomPhase.Settings.
"""
function Settings()
    ## allChannels DEFAULTS TO TRUE -- changed 14-Oct-2026, on the maintainer's instruction and on a
    ## measurement: with one target per subshell the total coherent response comes out SIGN-REVERSED
    ## (ratio -0.793 against the uncoupled total), while the complete channel set gives a physical +0.835.
    ## A default that produces the wrong sign is not a conservative default.
    ## directSolve DEFAULTS TO TRUE for the same reason -- 14-Oct-2026, the maintainer's instruction.  The
    ## iteration's spectral radius is about 2.7 and belongs to the OPERATOR, so an iterating default could not
    ## converge for any setting of anything else.  NOTE what this changes for a caller: on the direct route
    ## `maxIterations` bounds the KRYLOV DIMENSION and `nIterations` reports Krylov steps, not sweeps.
    Settings( Subshell[], Dict{Subshell,Subshell}(), false, 15, 1.0e-6, 0.3, DalgarnoLewis.Settings(), false, 0,
               true, true )
end


"""
`RandomPhase.Settings(settings::RandomPhase.Settings; activeSubshells=.., targetSubshells=.., couplingEnabled=..,
                       maxIterations=.., accuracy=.., mixing=.., dalgarnoLewis=.., printBefore=.., andersonDepth=..,
                       allChannels=.., directSolve=..)`
    ... the standard JAC keyword copy-constructor.
"""
function Settings(settings::RandomPhase.Settings;
    activeSubshells::Union{Nothing,Array{Subshell,1}}=nothing,   targetSubshells::Union{Nothing,Dict{Subshell,Subshell}}=nothing,
    couplingEnabled::Union{Nothing,Bool}=nothing,                 maxIterations::Union{Nothing,Int64}=nothing,
    accuracy::Union{Nothing,Float64}=nothing,                     mixing::Union{Nothing,Float64}=nothing,
    dalgarnoLewis::Union{Nothing,DalgarnoLewis.Settings}=nothing, printBefore::Union{Nothing,Bool}=nothing,
    andersonDepth::Union{Nothing,Int64}=nothing,                  allChannels::Union{Nothing,Bool}=nothing,
    directSolve::Union{Nothing,Bool}=nothing)

    if  isnothing(activeSubshells)   activeSubshellsx = settings.activeSubshells   else   activeSubshellsx = activeSubshells   end
    if  isnothing(targetSubshells)   targetSubshellsx = settings.targetSubshells   else   targetSubshellsx = targetSubshells   end
    if  isnothing(couplingEnabled)   couplingEnabledx = settings.couplingEnabled   else   couplingEnabledx = couplingEnabled   end
    if  isnothing(maxIterations)     maxIterationsx   = settings.maxIterations     else   maxIterationsx   = maxIterations     end
    if  isnothing(accuracy)          accuracyx        = settings.accuracy          else   accuracyx        = accuracy          end
    if  isnothing(mixing)            mixingx          = settings.mixing            else   mixingx          = mixing            end
    if  isnothing(dalgarnoLewis)     dalgarnoLewisx   = settings.dalgarnoLewis     else   dalgarnoLewisx   = dalgarnoLewis      end
    if  isnothing(printBefore)       printBeforex     = settings.printBefore       else   printBeforex     = printBefore       end
    if  isnothing(andersonDepth)     andersonDepthx   = settings.andersonDepth     else   andersonDepthx   = andersonDepth     end
    if  isnothing(allChannels)       allChannelsx     = settings.allChannels       else   allChannelsx     = allChannels       end
    if  isnothing(directSolve)       directSolvex     = settings.directSolve       else   directSolvex     = directSolve       end

    Settings( activeSubshellsx, targetSubshellsx, couplingEnabledx, maxIterationsx, accuracyx, mixingx, dalgarnoLewisx, printBeforex,
              andersonDepthx, allChannelsx, directSolvex )
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
    println(io, "andersonDepth:            $(settings.andersonDepth)  ")
    println(io, "allChannels:              $(settings.allChannels)  ")
    println(io, "directSolve:              $(settings.directSolve)  ")
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
    + channelOrbitals     ::Dict{Tuple{Subshell,Subshell},PerturbedOrbital}  ... EVERY channel's perturbed
                                                                                  orbital, keyed by the
                                                                                  (active, target) pair.
                                                                                  **Added 14-Oct-2026 because
                                                                                  `perturbedOrbitals` alone
                                                                                  makes the physics
                                                                                  unreachable:** with
                                                                                  `allChannels` one electron
                                                                                  carries several channels, the
                                                                                  observable is the COHERENT
                                                                                  SUM over all of them, and a
                                                                                  caller given only the primary
                                                                                  ones cannot form it -- it
                                                                                  gets a 3-channel sum of a
                                                                                  7-channel solve, which is a
                                                                                  different quantity.  Empty on
                                                                                  the uncoupled path.
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
    channelOrbitals      ::Dict{Tuple{Subshell,Subshell},DalgarnoLewis.PerturbedOrbital}
end


# `Base.show(io::IO, response::RandomPhase.CoupledResponse)`  ... prepares a proper printout of response::RandomPhase.CoupledResponse.
function Base.show(io::IO, response::RandomPhase.CoupledResponse)
    println(io, "multipole, gauge:         $(response.multipole), $(response.gauge)  ")
    println(io, "omega:                    $(response.omega)  ")
    println(io, "active subshells:         $(collect(keys(response.perturbedOrbitals)))  ")
    println(io, "nIterations, converged:   $(response.nIterations), $(response.converged)  ")
end


"""
`RandomPhase.angularCoefficient(L::Int64, tSub::Subshell, aSub::Subshell)`
    ... the one-electron angular factor the rank-L induced interaction carries on the channel aSub -> tSub,
        namely `AngularMomentum.CL_reduced_me(tSub, L, aSub)`, or 0.0 where the triangular and parity
        selection rules forbid the channel. It exists as a named function, rather than inline, because the
        convention is the thing most easily got wrong here and Rule 18 asks for it to be visible.

        **THE ARGUMENT ORDER IS THE CALLER'S BUSINESS, AND THE TWO ENDS OF THE INDUCED TERM USE DIFFERENT ONES.**
        This function is the bare reduced matrix element with its selection rules attached, and nothing else;  it does
        not know which end of the induced interaction it is serving.  `RandomPhase.inducedPotentialZr` calls it
        SOURCE-FIRST, as `angularCoefficient(L, b, t_b)`, because the weight it needs is the sum over the magnetic
        substates of a filled subshell and that sum is CL_reduced_me(b, L, t_b)^2 / (2L+1) -- positive definite, so the
        order is not free.  The receiving end does NOT call this function at all:  a one-body potential's matrix element
        carries `InteractionStrength.eMultipole`'s own reduction, which is `RandomPhase.receivingCoefficient`.  The two
        orders differ by (-1)^(j_1 - j_2), i.e. by a SIGN on every channel with |Delta j| = 1, so getting one of them
        wrong corrupts the BALANCE between channels rather than the scale.

        **AND THERE IS NO `(-1)^L` ANYWHERE ON THE DIRECT PATH -- amended 18-Oct-2026, priority item 75.**  It was
        applied once in `RandomPhase.solve`, on the reading that the induced direct term is the two-electron matrix
        element <t_a t_b| 1/r12 |a b> that `InteractionStrength.XL_Coulomb` assembles.  The direct term is not assembled
        that way here:  a charge density is built, Poisson is solved for a potential, and that potential's ONE-BODY
        matrix element is taken at the other end, so the two-electron product's phase never enters.  EXCHANGE is a
        different object and does go through XL_Coulomb's convention, which is why
        `RandomPhase.exchangeAngularCoefficient` keeps target-first at both of its ends.
        A coefficient::Float64 is returned, which is 0.0 for a channel the selection rules forbid.
"""
function angularCoefficient(L::Int64, tSub::Subshell, aSub::Subshell)
    lt = Basics.subshell_l(tSub);   jt2 = Basics.subshell_2j(tSub)
    la = Basics.subshell_l(aSub);   ja2 = Basics.subshell_2j(aSub)
    if  AngularMomentum.triangularDelta(jt2+1, ja2+1, L+L+1) == 0   ||   rem(lt+la+L, 2) == 1
        return( 0.0 )
    end

    return( AngularMomentum.CL_reduced_me(tSub, L, aSub) )
end


"""
`RandomPhase.receivingCoefficient(L::Int64, aSub::Subshell, tSub::Subshell)`
    ... the one-electron angular factor the rank-L induced potential carries where it ACTS, i.e. at the active end of
        the channel aSub -> tSub, namely `AngularMomentum.CL_reduced_me(aSub, L, tSub) / sqrt(2*j_a+1)`, or 0.0 where
        the triangular and parity selection rules forbid the channel. It is the counterpart of
        `RandomPhase.angularCoefficient`, which is the factor the induced DENSITY carries at the other end, and they are
        deliberately two functions because they are two different reductions of the same algebra.

        **WHY IT IS NOT `angularCoefficient`, AND HOW THAT WAS DECIDED -- 18-Oct-2026, priority item 75.** The induced
        potential and the EXTERNAL drive are both ordinary one-body potentials, acting on the same electron, in the same
        Dalgarno-Lewis equation; and a one-body potential has ONE matrix element, so the same radial function v(r) must
        produce the same source vector down either path. `DalgarnoLewis.solve` reduces its drive with
        `InteractionStrength.eMultipole`, which is CL_reduced_me(a, L, t)/sqrt(2j_a+1) times the bare radial integral --
        source first, and DIVIDED. Pushing v(r) = r down both paths on the nine E1 channels of Ne, the induced path as
        it stood returned -sqrt(2j_a+1)*CL(t,L,a)/CL(a,L,t) times the drive instead of 1.0, i.e. sqrt(2j_a+1) in
        magnitude with a channel-dependent sign, reproducible to 2.0e-15.

        **AND THE DRIVE IS THE ANCHORED SIDE, which is what makes that an arbitration rather than one more comparison:**
        `eMultipole`'s normalisation is fixed by the EXACT alpha(H 1s) = 4.5 a.u. through the polarizability assembly
        alpha = (2/3)(2j_a+1) eMultipole(1, phi_a, chi), which reproduces 4.49972686 on this grid. So the drive's
        reduction was never in question and the induced term's was.
        A coefficient::Float64 is returned, which is 0.0 for a channel the selection rules forbid.
"""
function receivingCoefficient(L::Int64, aSub::Subshell, tSub::Subshell)
    la = Basics.subshell_l(aSub);   ja2 = Basics.subshell_2j(aSub)
    lt = Basics.subshell_l(tSub);   jt2 = Basics.subshell_2j(tSub)
    if  AngularMomentum.triangularDelta(ja2+1, jt2+1, L+L+1) == 0   ||   rem(la+lt+L, 2) == 1
        return( 0.0 )
    end

    return( AngularMomentum.CL_reduced_me(aSub, L, tSub) / sqrt(ja2 + 1.0) )
end


"""
`RandomPhase.reachableChannels(L::Int64, aSub::Subshell)`
    ... every target symmetry a rank-L multipole can reach from the subshell `aSub`, i.e. every kappa allowed
        by the triangular condition on j and by the parity rule l_t + l_a + L even. This is what an RPA's
        induced field has to be built from, and it needs no input from the caller because the selection rules
        determine it completely.

        **ONLY THE kappa OF THE RETURNED SUBSHELLS MATTERS.** `DalgarnoLewis.solve` takes a `targetSubshell`
        purely to read `targetSubshell.kappa` off it -- the channel is a symmetry, not a particular orbital --
        so each entry here is built with the LOWEST principal quantum number that symmetry admits, n = l+1,
        and that n carries no meaning. A caller must not read anything into it.

        Worked, because it is easy to under-count: from 2p_3/2 with L = 1 the reachable symmetries are
        s_1/2, d_3/2 and d_5/2 -- THREE channels, where `Settings.targetSubshells` can name only one. The p
        symmetries are excluded by parity and j = 7/2 by the triangular condition. An
        `Array{Subshell,1}` is returned, ordered by l and then by j, and it is empty only for a rank and
        symmetry that admit nothing at all. The order is by l and then by kappa, so the j = l+1/2 partner comes
        before the j = l-1/2 one.
"""
function reachableChannels(L::Int64, aSub::Subshell)
    la = Basics.subshell_l(aSub);   ja2 = Basics.subshell_2j(aSub)
    channels = Subshell[]
    for  lt = max(0, la-L):(la+L)
        rem(lt + la + L, 2) == 1   &&   continue
        for  kappa  in  (lt == 0 ? (-1,) : (-lt-1, lt))
            jt2 = abs(2*kappa) - 1
            AngularMomentum.triangularDelta(jt2+1, ja2+1, L+L+1) == 0   &&   continue
            push!(channels, Subshell(lt+1, kappa))
        end
    end

    return( channels )
end


"""
`RandomPhase.projectOntoUnoccupied(chi::DalgarnoLewis.PerturbedOrbital, targetSubshell::Subshell,
                                    occupiedOrbitals::Dict{Subshell,Orbital}, grid::Radial.Grid)`
    ... removes from chi every component along an OCCUPIED orbital of chi's own kappa, and returns the result
        as an `orbital::Orbital`. A channel label such as `5p_3/2 -> 3d_5/2` names a kappa SYMMETRY and not one
        target orbital, so the solved chi spans that whole symmetry -- filled members included. A core electron
        cannot be excited into a filled orbital, so those components are forbidden and must not appear in any
        observable built from chi.

        **WHY THIS EXISTS AS ITS OWN FUNCTION -- 16-Oct-2026, item 73(c).** The projection had been applied to
        the POLARIZABILITY assembly since 14-Oct and never to the induced DENSITY, although both are built from
        the same chi. The consequence was measured and is large: the density's rank-1 moment, which must equal
        the polarizability computed from the same chi, instead ran from -886 to +2081 times it across 42 core
        channels of Ba+. With the projection applied the ratio is **exactly 3.00000 on every one of them**, and
        that residue is the (2L+1) of `RandomPhase.inducedPotentialZr`. Two separately-correct pieces again
        failed to compose because nothing required them to.
"""
function projectOntoUnoccupied(chi::DalgarnoLewis.PerturbedOrbital, targetSubshell::Subshell,
                                occupiedOrbitals::Dict{Subshell,Orbital}, grid::Radial.Grid)
    P = copy(chi.P);   Q = copy(chi.Q)
    for  (sh, orb)  in  occupiedOrbitals
        (sh.kappa == targetSubshell.kappa  &&  orb.energy < 0.0)   ||   continue
        ov = RadialIntegrals.overlap(orb.P, P, grid) + RadialIntegrals.overlap(orb.Q, Q, grid)
        nn = RadialIntegrals.overlap(orb.P, orb.P, grid) + RadialIntegrals.overlap(orb.Q, orb.Q, grid)
        nn <= 0.0   &&   continue
        for  i = 1:min(length(P), length(orb.P))   P[i] = P[i] - ov/nn * orb.P[i]   end
        for  i = 1:min(length(Q), length(orb.Q))   Q[i] = Q[i] - ov/nn * orb.Q[i]   end
    end

    return( Orbital(chi.subshell, false, true, chi.targetEnergy, P, Q,
                    zeros(length(P)), zeros(length(Q)), Radial.Grid()) )
end


"""
`RandomPhase.inducedPotentialZr(multipole::EmMultipole, activeOrbitals::Dict{Subshell,Orbital},
                                 chis::Dict{Tuple{Subshell,Subshell},DalgarnoLewis.PerturbedOrbital},
                                 excludeSubshell::Subshell, grid::Radial.Grid;
                                 occupiedOrbitals::Union{Nothing,Dict{Subshell,Orbital}}=nothing,
                                 occupations::Union{Nothing,Dict{Subshell,Float64}}=nothing)`
    ... builds the DIRECT (Hartree, "ring") induced potential ONE electron of `excludeSubshell` feels from every
        OTHER active electron's CURRENT chi -- its `(N-1)` PARTNERS IN ITS OWN SUBSHELL INCLUDED, whenever
        `occupations` is given -- at the multipole's own rank L, in the SAME "Zr" (= -r*V(r)) convention
        `Radial.Potential` and `Bsplines.computeVlocal` use throughout JAC -- so that the minus sign below is
        the ordinary nuclear-potential-style sign, not a new convention. The induced density itself is
        2*[Pb*Pchi_b + Qb*Qchi_b] -- the psi*chi + chi*psi symmetrization for a real response -- which is
        exact only at omega = 0, since at finite omega the two terms are the two time orderings and are
        different functions.

        **THE WEIGHT IS sqrt(2j_b+1) TIMES THE CHANNEL'S REDUCED MATRIX ELEMENT, AND THE sqrt IS NOT A TYPO.**
        Until 14-Oct-2026 each electron b was weighted by the full (2j_b+1), treating it statistically as a
        filled shell, with no angular factor at all -- which made the coupled iteration diverge (see
        `RandomPhase.solve`). The weight is fixed, without any literature, by the L = 0 HARTREE LIMIT:

          * at L = 0 with t_b = b the induced potential must be the ordinary first-order Hartree potential,
            2*(2j_b+1)*Y_0[rho]/r, which is what the old code had and what is not in dispute;
          * the source this potential feeds is `DalgarnoLewis.solve`'s, which is REDUCED -- it comes from
            `DalgarnoLewis.amplitudeKernel`, a reduced one-electron matrix element -- so the induced term must
            be reduced to match:  `RandomPhase.angularCoefficient` at this end, and
            `RandomPhase.receivingCoefficient`, which divides by sqrt(2j_a+1), where the potential acts;
          * and `AngularMomentum.CL_reduced_me(b, 0, b) = sqrt(2j_b+1)` EXACTLY, checked for s, p, d and f.

        Equating the two therefore gives weight * sqrt(2j_b+1) = 2*(2j_b+1), i.e. a weight of
        2*sqrt(2j_b+1): the occupation splits in half, one sqrt(2j_b+1) staying here and the other being
        carried by the reduced matrix element's own normalisation. This is precisely the sqrt(2j+1) Rule 18
        is about, and it is written out because the two quantities it sits between -- an ordinary potential
        and a reduced source -- are otherwise denoted alike.

        **AND THE WEIGHT IS POSITIVE DEFINITE, WHICH AN EARLIER VERSION OF THIS DOCSTRING DENIED -- 18-Oct-2026,
        priority item 75.**  Written out as the explicit sum over the magnetic substates of a filled subshell, the
        weight of the channel b -> t_b is sum_m <b m|C^L_0|t_b m> <t_b m|C^L_0|b m> = CL_reduced_me(b, L, t_b)^2/(2L+1),
        where the (2L+1) is 3j orthogonality at fixed middle projection and not a calibration.  It cannot change sign.
        The text here previously argued the opposite -- that the two spin-orbit partners of a subshell enter the induced
        field with OPPOSITE signs, that the cancellation is "physical", and that it is "a large part of why the
        uncorrected term was far too strong".  THAT CLAIM IS WITHDRAWN:  the alternation was an artefact of calling
        `RandomPhase.angularCoefficient` target-first at this end, which differs from source-first by
        (-1)^(j_b - j_t_b) and so flips exactly the |Delta j| = 1 channels -- 25 of Ba+'s 42 core channels.

        **THE RADIAL PART IS VERIFIED AGAINST JAC'S OWN, TO TEN DIGITS.** `InteractionStrength.XL_Coulomb`
        builds the corresponding two-electron strength from `RadialIntegrals.SlaterRkComponent`, so the two
        can be compared entry by entry if the source orbital is chosen to BE a single B-spline, which makes
        its coefficient vector a unit vector. Done on the Ne case for three spline pairs, this route and
        JAC's agree to every digit printed -- ratio 1.000000, values -7.3537752978e-03 against
        -7.3537752978e-03 -- SIGN INCLUDED. So the radial machinery and the "Zr" sign convention here are the
        ones JAC uses elsewhere, and the angular factor was the only thing this module lacked.
        A VindZr::Array{Float64,1} of length grid.NoPoints is returned.
"""
function inducedPotentialZr(multipole::EmMultipole, activeOrbitals::Dict{Subshell,Orbital},
                             chis::Dict{Tuple{Subshell,Subshell},DalgarnoLewis.PerturbedOrbital},
                             excludeSubshell::Subshell, grid::Radial.Grid;
                             occupiedOrbitals::Union{Nothing,Dict{Subshell,Orbital}}=nothing,
                             occupations::Union{Nothing,Dict{Subshell,Float64}}=nothing,
                             withExchange::Bool=false)
    rho = zeros(grid.NoPoints);   mtpRho = 1
    # EVERY CHANNEL OF EVERY OTHER ELECTRON CONTRIBUTES, EACH WITH ITS OWN ANGULAR FACTOR.  `chis` is keyed
    # by the (active, target) PAIR, so the target is read off the key rather than looked up -- which is what
    # lets one electron carry several channels, as an RPA requires.
    #   **AN ELECTRON IS SCREENED BY ITS OWN SUBSHELL PARTNERS, AND UNTIL 17-Oct-2026 IT WAS NOT.**  The line
    # here was `b == excludeSubshell && continue`, which skipped the WHOLE source subshell and so denied every
    # electron its (N-1) partners -- 5 of 6 for a p^6 shell, and for a TWO-ELECTRON atom all of them, so He
    # came out with EXACTLY ZERO screening although two electrons must certainly screen each other.
    #   The self-term needs the ACTUAL occupation and cannot be guessed from 2j+1:  a one-electron subshell has
    # no partners at all, so assuming a closed shell would invent self-screening for hydrogen and break both
    # the exact alpha(H 1s) = 4.5 anchor and `example-At.jl` branch (a), whose whole assertion is that a single
    # active electron induces nothing.  Hence `occupations`, and hence the CONSERVATIVE default: with no
    # occupations given the partners are omitted exactly as before, which is incomplete but never wrong.
    for  ((b, tb), chiB)  in  chis
        isnothing(occupations)  &&  b == excludeSubshell   &&   continue
        orbB   = activeOrbitals[b];   jb2 = Basics.subshell_2j(b)
        # The weight is sqrt(2j_b+1), NOT (2j_b+1):  see this method's docstring for the Hartree calibration
        # that fixes it.  The remaining sqrt(2j_b+1) of the occupation is carried by CL_reduced_me's own
        # normalisation at this end, since CL_reduced_me(b, 0, b) = sqrt(2j_b+1) exactly.
        # THE (2L+1), DERIVED AND NOT FITTED -- 16-Oct-2026, item 73(c).  The density's rank-L moment
        # int r^L rho dr is the induced multipole moment, which for L = 1 is the quantity the POLARIZABILITY
        # assembly computes from this very chi -- and that assembly reproduces the exact 4.5 a.u. of H 1s.
        # Measured on 42 core channels of Ba+ spanning j = 1/2...7/2 and l = 0...3, the two agreed to
        # EXACTLY 3.00000 = 2L+1, with zero spread, once the projection below was applied.  The L = 0
        # Hartree calibration that fixes the sqrt(2j_b+1) is silent about this factor, since 2L+1 = 1 there.
        # **THE ARGUMENT ORDER IS SOURCE-FIRST, AND THE MINUS IS PART OF IT -- derived and measured 18-Oct-2026,
        # priority item 75.**  This end turns a channel's chi into the induced CHARGE DENSITY, and the weight that does
        # so is the explicit sum over the magnetic substates of a filled subshell,
        #     w(b -> t_b)  =  sum_m <b m|C^L_0|t_b m> <t_b m|C^L_0|b m>  =  CL_reduced_me(b, L, t_b)^2 / (2L+1) ,
        # the (2L+1) coming from 3j orthogonality at fixed middle projection rather than from a calibration.  The sum is
        # POSITIVE DEFINITE, and that is what decides the argument order:  the sign of a channel's contribution to the
        # density is the sign of its chi and nothing else.
        #   **SO ITEM 73(c)'s TARGET-FIRST ORDER IS WITHDRAWN.**  It was taken from `InteractionStrength.XL_Coulomb`'s
        # two-electron convention, X^L(abcd) = (-1)^L CL(a,L,c) CL(b,L,d) R^L(abcd), on the reading that the induced
        # direct term is <t_a t_b| 1/r12 |a b>.  That reading is right about the two-electron matrix element and wrong
        # here, because this module does not assemble one:  it builds a density, solves Poisson for a potential, and
        # takes that potential's ONE-BODY matrix element at the other end (`RandomPhase.receivingCoefficient`).
        # Target-first differs from source-first by (-1)^(j_b - j_t_b), so it FLIPPED the 25 of Ba+'s 42 core channels
        # that have |Delta j| = 1 and made the weight alternate in sign between spin-orbit partners.
        #   The MINUS is JAC's chi convention made visible:  `DalgarnoLewis` returns the NEGATIVE of the first-order
        # correction, so a density built from phi*chi carries the opposite sign to the physical one.  With the minus
        # here the induced field OPPOSES the drive, which is what makes the response screen;  without it the sign of
        # every screening factor this module produces is wrong.
        weightB = -sqrt(jb2 + 1.0) * RandomPhase.angularCoefficient(multipole.L, b, tb) / (2*multipole.L + 1)
        # ONE electron is removed from its own subshell, never the whole subshell.  The weight above carries
        # the occupation N_b (half here as sqrt(2j_b+1), half through CL_reduced_me's normalisation), so an
        # electron of `excludeSubshell` feels (N-1)/N of that subshell's induced density: its partners screen
        # it, only IT does not.
        #   **THE REPAIR MAKES THE AGREEMENT WORSE ON ITS OWN, AND THAT IS EXPECTED AND NOT A REGRESSION:** it
        # strengthens the direct term and so LOWERS alpha, while measured alpha is already too low for He and
        # Ne.  The partners it restores were an uncontrolled stand-in for the self-interaction that EXCHANGE
        # cancels properly, so the two belong together and a number measured before this repair is not a
        # baseline for the exchange work.
        #   **AND WITH EXCHANGE SWITCHED ON THE SELF TERM IS KEPT IN FULL**, which is not an inconsistency but
        # the whole point: in a proper RPAE the direct term runs over ALL N_b electrons, self included, and the
        # spurious self-interaction is cancelled by EXCHANGE rather than by hand.  **DO NOT CALIBRATE THAT ON A
        # ONE-ELECTRON SYSTEM, 20-Oct-2026, item 75:** hydrogen is the one case where this method's own
        # N_b/(2j_b+1) is not 1, and its nuclear-only left-hand side drops a term the identity needs -- see
        # `RandomPhase.exchangeSourceVector` for the two measurements and for what IS settled.  Removing the
        # self term here as well would remove it twice.
        if  !isnothing(occupations)
            nb = get(occupations, b, jb2 + 1.0)
            # A PARTLY FILLED SUBSHELL CONTRIBUTES ONLY ITS SHARE.  The weight above is built for a FULL
            # subshell (2j_b+1 electrons), so a general occupation enters as the ratio N_b/(2j_b+1).  For a
            # closed subshell the ratio is exactly 1, which is why this cannot disturb a closed-shell result.
            # Without it a one-electron subshell would screen as though it held two, and hydrogen is the case
            # that matters: it is where the exchange normalisation is fixed.
            weightB = weightB * nb / (jb2 + 1.0)
            # and ONE electron of its own subshell is removed, unless exchange is there to do it properly
            if  b == excludeSubshell  &&  !withExchange   weightB = weightB * (nb - 1.0) / nb   end
        end
        weightB == 0.0   &&   continue
        # THE PAULI PROJECTION.  A channel label names a kappa SYMMETRY, so chi spans the filled members of
        # it too, and a core electron cannot be excited into a filled orbital.  Without this the forbidden
        # components swamp the density -- the moment ratio above ran from -886 to +2081 instead of 3.
        chiP = RandomPhase.projectOntoUnoccupied(chiB, tb, isnothing(occupiedOrbitals) ? activeOrbitals :
                                                 occupiedOrbitals, grid)
        mtpB = min(length(orbB.P), length(chiP.P), length(orbB.Q), length(chiP.Q))
        for  i = 1:mtpB   rho[i] = rho[i] + weightB * 2.0 * (orbB.P[i]*chiP.P[i] + orbB.Q[i]*chiP.Q[i])   end
        mtpRho = max(mtpRho, mtpB)
    end

    return( [ -RadialIntegrals.Yk_ab(multipole.L, grid.r[i], rho, mtpRho, grid)   for i = 1:grid.NoPoints ] )
end


"""
`RandomPhase.exchangeAngularCoefficient(L::Int64, Lp::Int64, tSubA::Subshell, aSub::Subshell,
                                        tSubB::Subshell, bSub::Subshell, term::Symbol)`
    ... the angular coefficient the EXCHANGE part of the induced interaction carries on the channel
        aSub -> tSubA, from the partner channel bSub -> tSubB, at Coulomb rank `Lp` and perturbation rank `L`.
        `term` is `:E1` or `:E2`, the two exchange terms (below). A coefficient::Float64 is returned, 0.0 where
        the selection rules forbid it.

        **WHY THERE ARE TWO TERMS AND WHY THEY DIFFER.** Varying the Hartree-Fock potential
        `V_HF phi_a = sum_b [ <b|b> phi_a - <b|phi_a> phi_b ]` to first order gives two exchange pieces with
        DIFFERENT pairings of the four subshells:

            E1 = -sum_b <chi_b|1/r12|phi_a> phi_b      pairs (t_a, b)   and (t_b, a)
            E2 = -sum_b <phi_b|1/r12|phi_a> chi_b      pairs (t_a, t_b) and (b,   a)

        **WHY A 6j APPEARS AT ALL, where the direct term needs none.** The direct term carries only the
        perturbation's rank L at both ends, so one `CL_reduced_me` per end suffices. Exchange carries a COULOMB
        rank Lp that differs from L, so the free index must be recoupled from Lp back to L. That recoupling is
        the 6j, and the two terms use the same one with `j_b` and `j_tb` swapped in the bottom row:

            E1:  (-1)^(j_tb + j_b + L + Lp) * 6j{ j_ta  j_a  L ;  j_tb  j_b   Lp }
            E2:  (-1)^(j_tb + j_b + L + Lp) * 6j{ j_ta  j_a  L ;  j_b   j_tb  Lp }

        **NONE OF THIS WAS RECALLED FROM A TEXTBOOK, AND THAT IS DELIBERATE -- Rule 18.** The m-resolved form
        uses only Wigner-Eckart, which is unambiguous once the 3j convention is fixed, and the coefficient was
        then READ OFF numerically and matched against enumerated 6j argument orderings.

        **COVARIANCE IS A COMPLETE TEST OF THE STRUCTURE AND BLIND TO THE SCALE -- amended 20-Oct-2026, item
        75.** It fixes the 6j, the pairing and the phase, and it caught both of the author's errors below. What
        it cannot see is a factor depending only on (j_b, L, Lp), an overall normalisation included, since such a
        factor leaves every magnetic quantum number untouched. The scale is fixed elsewhere and by other means:
        for the E2 branch by the frozen-Hartree comparison recorded in `RandomPhase.exchangeSourceVector`, which
        holds for a closed subshell; for the E1 branch not at all, because the frozen NONLOCAL exchange it would
        have to be compared against does not exist in JAC. Do not read the 1.11e-16 below as a check on the
        magnitude.

        **COVARIANCE, THEN, IN ITS OWN RIGHT.** The m-resolved term must be
        PROPORTIONAL to `ME(t_a, m_ta, L, q, a, m_a)`, i.e. the ratio must be independent of every magnetic
        quantum number, because the problem is rotationally covariant. Only a correct form satisfies that, and
        it caught both of the author's errors: a first E1 built with the FORWARD chi amplitude showed
        m-dependence of 3.3e-1, as large as the coefficient itself -- `<chi_b|` is a BRA and needs the
        bra/ket-SWAPPED amplitude -- and E2 had been omitted altogether, which only became visible because it
        is covariant on its own to 1e-17.

        **THE PHASE WAS VALIDATED OUT OF SAMPLE.** It was read off 11 cases, so it was then checked on an
        exhaustive sweep of j = 1/2...7/2 and L, Lp = 0...3: **970 informative cases per term, worst deviation
        1.11e-16**. A formula tested on the data that produced it is a fit and not a check. A further 3126
        combinations vanish on both sides and carry no information; they were counted and excluded rather than
        allowed to inflate the agreement.
"""
function exchangeAngularCoefficient(L::Int64, Lp::Int64, tSubA::Subshell, aSub::Subshell,
                                    tSubB::Subshell, bSub::Subshell, term::Symbol)
    jta2 = Basics.subshell_2j(tSubA);   ja2 = Basics.subshell_2j(aSub)
    jtb2 = Basics.subshell_2j(tSubB);   jb2 = Basics.subshell_2j(bSub)
    # The PAIRING of the two reduced C^Lp matrix elements is what distinguishes the two terms, and it is the
    # same pairing the radial assembly in `exchangeSourceVector` uses -- the two must be read together.
    if       term == :E1
        cc = AngularMomentum.CL_reduced_me(tSubA, Lp, bSub)  * AngularMomentum.CL_reduced_me(tSubB, Lp, aSub)
        w6 = AngularMomentum.Wigner_6j(AngularJ64(jta2//2), AngularJ64(ja2//2), AngularJ64(L),
                                       AngularJ64(jtb2//2), AngularJ64(jb2//2), AngularJ64(Lp))
    elseif   term == :E2
        cc = AngularMomentum.CL_reduced_me(tSubA, Lp, tSubB) * AngularMomentum.CL_reduced_me(bSub, Lp, aSub)
        w6 = AngularMomentum.Wigner_6j(AngularJ64(jta2//2), AngularJ64(ja2//2), AngularJ64(L),
                                       AngularJ64(jb2//2), AngularJ64(jtb2//2), AngularJ64(Lp))
    else
        error("\n\nRandomPhase.exchangeAngularCoefficient: `term` must be :E1 or :E2, not $(term).\n")
    end
    (cc == 0.0  ||  w6 == 0.0)   &&   return( 0.0 )

    return( (-1.0)^((jtb2 + jb2)//2 + L + Lp) * w6 * cc )
end


"""
`RandomPhase.exchangeSourceVector(multipole::EmMultipole, activeOrbitals::Dict{Subshell,Orbital},`
                                 `chis::Dict{Tuple{Subshell,Subshell},DalgarnoLewis.PerturbedOrbital},`
                                 `aSub::Subshell, tSubA::Subshell, grid::Radial.Grid,`
                                 `primitives::Bsplines.Primitives;`
                                 `occupations::Union{Nothing,Dict{Subshell,Float64}}=nothing, maxRank::Int64=8)`
    ... builds the EXCHANGE induced source for the channel aSub -> tSubA, summed over every partner channel in
        `chis` and over the Coulomb ranks Lp, in the same B-spline projection and sign convention
        `RandomPhase.inducedSourceVector` uses. An extraSource::Array{Float64,1} of length
        primitives.grid.nsL + primitives.grid.nsS is returned.

        **IT IS NOT A POTENTIAL, which is why it needs its own function rather than a weight.** The direct term
        returns a radial potential multiplying `phi_a`; exchange returns `phi_b` and `chi_b` instead, each
        multiplied by a two-orbital Yukawa function -- so there is nothing to hand to `inducedSourceVector`.
        The two terms are, with the angular coefficients from
        `RandomPhase.exchangeAngularCoefficient`:

            E1:  Y^Lp[chi_b, phi_a](r) * phi_b(r)
            E2:  Y^Lp[phi_b, phi_a](r) * chi_b(r)

        **AND IT SUMS OVER RANKS, which changes the SHAPE of the computation and not just a weight.** The
        direct term carries only the perturbation's own rank L; this one needs every Lp the selection rules
        allow, so a case with 42 channels and about five ranks costs some thousands of Yukawa sweeps per
        iteration rather than tens.

        **THE OCCUPATION IS ALREADY IN THE COEFFICIENT, so do not apply it twice.** The derivation summed over
        every m_b of the partner subshell, which for a filled subshell IS its occupation 2j_b+1; `occupations`
        therefore enters only as the RATIO N_b/(2j_b+1) for a subshell that is not full. For the same reason
        the self term b = a is KEPT HERE IN FULL: it is exactly what cancels the direct term's spurious
        self-interaction, and `inducedPotentialZr` must then be called with `withExchange=true` so that it does
        not remove that self-interaction a second time.

        **HYDROGEN IS NOT THE TEST OF THIS NORMALISATION, AND THE CLAIM THAT IT IS HAS BEEN WITHDRAWN --
        measured 20-Oct-2026, priority item 75.** The appeal is obvious and it is a trap: for one electron
        V_HF = J - K vanishes identically, so a reader concludes that the induced direct and induced exchange
        sources must be the same vector and that comparing them needs no reference value. Two measurements say
        otherwise.

        FIRST, THE IDENTITY HAS THREE TERMS AND THE THIRD ONE CANNOT BE FORMED HERE. Varying
        (J[phi] - K[phi]) phi = 0, which holds for ANY single orbital, gives delta(J-K)[chi] phi_a + (J-K) chi = 0:
        the FROZEN potential acting on chi is part of it. The problem carries exactly two radial shapes,
        v1 = phi_a Y_1[phi_a, chi] and v2 = chi Y_0[phi_a, phi_a], and the identity splits into
        `v1:  A - B1 - C_K = 0` and `v2:  -B2 + C_J = 0`, where A is the induced direct coefficient, B1 and B2
        the two induced exchange coefficients, and C_J, C_K the frozen Hartree and frozen exchange ones. JAC has
        no NONLOCAL exchange operator, so C_K does not exist and the v1 equation is not checkable. Comparing the
        direct source against this function's output tests neither equation: measured on H 1s in a 120 a.u. box
        the ratios are 1.402708 (2p_3/2) and 1.811708 (2p_1/2) with cos 0.826089 and 0.899776, and both are
        REPRODUCED EXACTLY by the split below -- they are the angle between a pure v1 and a v1/v2 mixture, not a
        defect.

        SECOND, HYDROGEN IS THE ONE CASE WHERE THE FACTOR UNDER TEST IS NOT 1. `scale` above is N_b/(2j_b+1),
        which is exactly 1 for every closed subshell and 1/2 for a single electron -- so the open-shell
        occupation rule, the only thing a hydrogen comparison could calibrate, is also the only thing hydrogen
        cannot be calibrated against. Its left-hand side is nuclear-only as well, so C_J = 0 there while a
        many-electron mean field supplies it.

        **WHAT IS SETTLED, AND IT IS THE PART THAT MATTERS FOR A CLOSED SHELL.** The v2 equation IS checkable
        and it holds: C_J is the frozen Hartree coefficient of the partner subshell, N_b, and B2 comes out as
        `1.0` EXACTLY for every subshell and every multipole rank -- analytically, since at Lp = 0 the 6j with a
        zero argument is 1/sqrt((2j_a+1)(2j_t+1)) and cancels both reduced matrix elements. With `scale` = 1 for
        a closed subshell, B2 = 1 = C_J/N_b per electron and the E2 branch is correctly normalised. The radial
        assembly is sound too: the E1 term is collinear with the direct term to cos = -1.000000000000, the E2
        term lies along v2 to a residual of 2e-16, and E1 + E2 reproduces this function to 3e-16.
        `tools/mResolvedPolarizability.jl` CASE L carries the measurement.
"""
function exchangeSourceVector(multipole::EmMultipole, activeOrbitals::Dict{Subshell,Orbital},
                              chis::Dict{Tuple{Subshell,Subshell},DalgarnoLewis.PerturbedOrbital},
                              aSub::Subshell, tSubA::Subshell, grid::Radial.Grid,
                              primitives::Bsplines.Primitives;
                              occupations::Union{Nothing,Dict{Subshell,Float64}}=nothing, maxRank::Int64=8)
    np   = grid.NoPoints;   nsL = grid.nsL;   nsS = grid.nsS
    fP   = zeros(np);   fQ = zeros(np);   rho = zeros(np);   mtpF = 1
    orbA = activeOrbitals[aSub]
    for  ((bSub, tSubB), chiB)  in  chis
        !haskey(activeOrbitals, bSub)   &&   continue
        orbB  = activeOrbitals[bSub];   jb2 = Basics.subshell_2j(bSub)
        scale = isnothing(occupations) ? 1.0 : get(occupations, bSub, jb2 + 1.0) / (jb2 + 1.0)
        scale == 0.0   &&   continue
        for  Lp = 0:maxRank
            c1 = RandomPhase.exchangeAngularCoefficient(multipole.L, Lp, tSubA, aSub, tSubB, bSub, :E1)
            c2 = RandomPhase.exchangeAngularCoefficient(multipole.L, Lp, tSubA, aSub, tSubB, bSub, :E2)
            (c1 == 0.0  &&  c2 == 0.0)   &&   continue
            if  c1 != 0.0
                mtp = min(length(chiB.P), length(orbA.P), length(chiB.Q), length(orbA.Q))
                fill!(rho, 0.0)
                for  i = 1:mtp   rho[i] = chiB.P[i]*orbA.P[i] + chiB.Q[i]*orbA.Q[i]   end
                mtpB = min(length(orbB.P), length(orbB.Q))
                for  i = 2:mtpB
                    wy = RadialIntegrals.Yk_ab(Lp, grid.r[i], rho, mtp, grid) / grid.r[i]
                    fP[i] = fP[i] - c1 * scale * wy * orbB.P[i]
                    fQ[i] = fQ[i] - c1 * scale * wy * orbB.Q[i]
                end
                mtpF = max(mtpF, mtpB)
            end
            if  c2 != 0.0
                mtp = min(length(orbB.P), length(orbA.P), length(orbB.Q), length(orbA.Q))
                fill!(rho, 0.0)
                for  i = 1:mtp   rho[i] = orbB.P[i]*orbA.P[i] + orbB.Q[i]*orbA.Q[i]   end
                mtpC = min(length(chiB.P), length(chiB.Q))
                for  i = 2:mtpC
                    wy = RadialIntegrals.Yk_ab(Lp, grid.r[i], rho, mtp, grid) / grid.r[i]
                    fP[i] = fP[i] - c2 * scale * wy * chiB.P[i]
                    fQ[i] = fQ[i] - c2 * scale * wy * chiB.Q[i]
                end
                mtpF = max(mtpF, mtpC)
            end
        end
    end
    # The projection is `inducedSourceVector`'s own, with the general function in place of V*phi_a.
    extraSource = zeros(nsL+nsS)
    for  i = 1:nsL
        bsp = primitives.bsplinesL[i];   lower = bsp.lower;   upper = min(bsp.upper, mtpF)
        upper < lower   &&   continue
        wa = 0.0
        for  r = lower:upper   wa = wa + bsp.bs[r-bsp.lower+1] * fP[r] * grid.wr[r]   end
        extraSource[i] = wa
    end
    for  i = 1:nsS
        bsp = primitives.bsplinesS[i];   lower = bsp.lower;   upper = min(bsp.upper, mtpF)
        upper < lower   &&   continue
        wa = 0.0
        for  r = lower:upper   wa = wa + bsp.bs[r-bsp.lower+1] * fQ[r] * grid.wr[r]   end
        extraSource[nsL+i] = wa
    end

    return( extraSource )
end


"""
`RandomPhase.inducedSourceVector(VindZr::Array{Float64,1}, sourceOrbital::Orbital, primitives::Bsplines.Primitives)`
    ... projects the induced potential VindZr (in the SAME "Zr" convention as `inducedPotentialZr` returns)
        onto every B-spline primitive against sourceOrbital, exactly reproducing `Bsplines.computeVlocal`'s own
        formula and sign -- the one difference being that computeVlocal's SECOND argument is itself a single
        primitive, where here it is a general tabulated orbital, so the shared formula is written out again
        rather than called. This is `DalgarnoLewis.solve`'s `extraSource` keyword's own input.

        **IT CARRIES NO ANGULAR COEFFICIENT, DELIBERATELY, AND THE CALLER MUST SUPPLY ONE.** What this
        function returns is the bare radial projection;  the reduction belonging to the active channel is applied by
        `RandomPhase.solve` and `RandomPhase.sourceMap` through `RandomPhase.receivingCoefficient`. Until 14-Oct-2026
        nothing applied any angular factor at all, which is what made the coupled iteration diverge. A caller that uses
        this function directly, as a diagnostic does, must multiply by `receivingCoefficient(L, active, target)` --
        ACTIVE first, divided by sqrt(2j_active+1), and with no (-1)^L.
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
`RandomPhase.sourceMap(multipole::EmMultipole, gauge::EmGauge, omega::Float64, activeOrbitals::Dict{Subshell,Orbital},
                        channels::Array{Tuple{Subshell,Subshell},1}, sources::Array{Float64,1},
                        potential::Radial.Potential, primitives::Bsplines.Primitives,
                        dlSettings::DalgarnoLewis.Settings;
                        occupiedOrbitals::Union{Nothing,Dict{Subshell,Orbital}}=nothing,
                        occupations::Union{Nothing,Dict{Subshell,Float64}}=nothing,
                        staticDrive::Bool=false)`
    ... one application of the coupled map: given the induced source of EVERY channel, concatenated into a
        single vector, solve each channel's Dalgarno-Lewis equation and rebuild every channel's source from
        the resulting chi's. This is exactly one sweep of the fixed-point iteration, isolated so that it can
        serve as the matrix-vector product of a linear solver as well.

        **THE MAP IS AFFINE, AND THAT IS WHAT MAKES THE DIRECT SOLVE POSSIBLE.** Writing it as
        F(s) = b + M s, the drive b is simply F(0) -- the sources rebuilt from the UNCOUPLED chi's -- and the
        action of the linear part follows as M s = F(s) - b, with no matrix ever formed. The fixed point is
        the solution of (I - M) s = b, which exists whenever (I - M) is invertible, whether or not the
        spectral radius of M is below one.
        **`staticDrive` IS PASSED STRAIGHT THROUGH TO `DalgarnoLewis.solve`, and at omega = 0 it is not
        optional** -- the dynamic kernel is identically zero there, so every channel would be driven by
        nothing. `RandomPhase.solve` raises rather than let that happen; this function is the lower level and
        trusts its caller.
        A newSources::Array{Float64,1} of the same length is returned, together with the chis it passed
        through, as a tuple (newSources, chis).
"""
function sourceMap(multipole::EmMultipole, gauge::EmGauge, omega::Float64, activeOrbitals::Dict{Subshell,Orbital},
                    channels::Array{Tuple{Subshell,Subshell},1}, sources::Array{Float64,1},
                    potential::Radial.Potential, primitives::Bsplines.Primitives,
                    dlSettings::DalgarnoLewis.Settings;
                    occupiedOrbitals::Union{Nothing,Dict{Subshell,Orbital}}=nothing,
                    occupations::Union{Nothing,Dict{Subshell,Float64}}=nothing,
                    staticDrive::Bool=false)
    grid = primitives.grid;   m = grid.nsL + grid.nsS
    chis = Dict{Tuple{Subshell,Subshell},DalgarnoLewis.PerturbedOrbital}()
    for  (k, ch)  in  enumerate(channels)
        a = ch[1];   t = ch[2];   slice = sources[(k-1)*m+1 : k*m]
        chis[ch] = DalgarnoLewis.solve(multipole, gauge, omega, t, activeOrbitals[a], potential, primitives;
                                        settings=dlSettings, extraSource=slice, staticDrive=staticDrive)
    end
    newSources = zeros(length(sources))
    for  (k, ch)  in  enumerate(channels)
        a = ch[1];   t = ch[2]
        VindZr = RandomPhase.inducedPotentialZr(multipole, activeOrbitals, chis, a, grid;
                                                 occupiedOrbitals=occupiedOrbitals, occupations=occupations)
        # the receiving end's reduction, `InteractionStrength.eMultipole`'s own -- see
        # `RandomPhase.receivingCoefficient` and the same lines in `RandomPhase.solve`
        aCoeff = RandomPhase.receivingCoefficient(multipole.L, a, t)
        newSources[(k-1)*m+1 : k*m] = aCoeff * RandomPhase.inducedSourceVector(VindZr, activeOrbitals[a], primitives)
    end

    return( newSources, chis )
end


"""
`RandomPhase.gmres(applyA::Function, b::Array{Float64,1}, maxIterations::Int64, tolerance::Float64)`
    ... solves A x = b for a linear operator available only through its action `applyA`, by the standard
        Arnoldi/GMRES construction: build an orthonormal Krylov basis of b, A b, A^2 b, ..., and pick the
        combination minimising the residual over it.

        **WHY A KRYLOV SOLVER HERE AND NOT A MATRIX FACTORISATION.** The coupled system's operator is
        (I - M) on the concatenated channel sources, and M is never formed: assembling it in the B-spline
        basis would need one `RadialIntegrals.Yk_ab` sweep per column, which is O(N^2) per column and so
        hours of work for a case the map itself sweeps in seconds. GMRES needs only the action, and one
        action is one `RandomPhase.sourceMap` call. It is also indifferent to the spectral radius -- which
        is the whole point, since that radius is about 2.7 here and no fixed-point iteration can converge.

        This is unrelated to the Krylov/Lanczos CI eigensolver that was measured and DECLINED in Aug-2026:
        that was an EIGENVALUE problem for CI matrices, rejected because diagonalisation was a small share of
        wall time, and it proposed a new package dependency. This is a LINEAR SYSTEM whose iteration provably
        cannot converge, and it uses nothing beyond `LinearAlgebra`.

        A tuple (x, iterations, relativeResidual) is returned. The residual is the true GMRES residual
        ||b - A x|| / ||b||, so a caller can tell a converged solve from an exhausted one.
"""
function gmres(applyA::Function, b::Array{Float64,1}, maxIterations::Int64, tolerance::Float64)
    n = length(b);   beta = LinearAlgebra.norm(b)
    beta == 0.0   &&   return( zeros(n), 0, 0.0 )
    kMax = min(maxIterations, n)
    Q = zeros(n, kMax+1);   H = zeros(kMax+1, kMax);   Q[:,1] = b / beta
    for  k = 1:kMax
        w = applyA(Q[:,k])
        for  i = 1:k
            H[i,k] = LinearAlgebra.dot(Q[:,i], w);   w = w - H[i,k] * Q[:,i]
        end
        H[k+1,k] = LinearAlgebra.norm(w)
        e1 = zeros(k+1);   e1[1] = beta
        y  = H[1:k+1, 1:k] \ e1
        relative = LinearAlgebra.norm(H[1:k+1,1:k] * y - e1) / beta
        # a vanishing subdiagonal means the Krylov space is exhausted and y is then EXACT, so stop either way
        if  relative < tolerance   ||   H[k+1,k] < 1.0e-13 * beta
            return( Q[:,1:k] * y, k, relative )
        end
        Q[:,k+1] = w / H[k+1,k]
    end
    e1 = zeros(kMax+1);   e1[1] = beta
    y  = H[1:kMax+1, 1:kMax] \ e1

    return( Q[:,1:kMax] * y, kMax, LinearAlgebra.norm(H[1:kMax+1,1:kMax] * y - e1) / beta )
end


"""
`RandomPhase.solve(multipole::EmMultipole, gauge::EmGauge, omega::Float64, activeOrbitals::Dict{Subshell,Orbital},
                    potential::Radial.Potential, primitives::Bsplines.Primitives;
                    settings::RandomPhase.Settings=RandomPhase.Settings(),
                    occupiedOrbitals::Union{Nothing,Dict{Subshell,Orbital}}=nothing,
                    occupations::Union{Nothing,Dict{Subshell,Float64}}=nothing,
                    staticDrive::Bool=false)`
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
        **WHAT THIS FUNCTION HAS BEEN VALIDATED TO DO, AND WHAT IT STILL CANNOT DO -- 18-Oct-2026, priority
        item 62 closed.** At omega = 0 with `staticDrive = true` the coupled response has been checked against an
        INDEPENDENT m-resolved implementation -- explicit magnetic substates, no angular reduction, no 6j, no
        CL convention -- function by function rather than in a total: on the nine E1 channels of Ne the chi
        returned here is CL(a,1,t_a)/sqrt(2j_a+1) times that calculation's, the whole difference being the
        normalisation `InteractionStrength.eMultipole` puts on the drive, to 9.9e-11 against a GMRES tolerance of
        1e-8. **Checking the TOTAL would not have shown this:** the two angular conventions tried before
        `cbdedacb` agree to 0.7 % on Ne's screened total while differing by 0.644 to 2.194 PER CHANNEL, which is
        how a convention error survived two repairs here.
        **WHAT IS STILL ABSENT IS EXCHANGE** (item 75, HELD rather than open since 20-Oct-2026: its one blocker is
        a frozen NONLOCAL exchange operator, which JAC does not have anywhere; the verified lines are in this file
        and deliberately unwired). A direct-only induced field over-screens, so a number from this function is a
        direct-only result in a strongly-coupled regime and NOT a validated screening correction. Do not read the
        paragraph above as more than it says: it establishes that two routes compute the SAME thing, not that the
        thing is the measured one.
        A response::RandomPhase.CoupledResponse is returned.
"""
function solve(multipole::EmMultipole, gauge::EmGauge, omega::Float64, activeOrbitals::Dict{Subshell,Orbital},
                potential::Radial.Potential, primitives::Bsplines.Primitives;
                settings::RandomPhase.Settings=RandomPhase.Settings(),
                occupiedOrbitals::Union{Nothing,Dict{Subshell,Orbital}}=nothing,
                occupations::Union{Nothing,Dict{Subshell,Float64}}=nothing,
                staticDrive::Bool=false)
    # THE STATIC CASE MUST SAY SO, AND UNTIL 18-Oct-2026 IT COULD NOT -- item 69.  `DalgarnoLewis.solve` gained
    # `staticDrive` on 14-Oct, which swaps the photon-normalised multipole field for `eMultipole`'s
    # frequency-independent r^L C_L;  this module had no pass-through, so a call at omega = 0 drove every
    # channel with an `amplitudeKernel` that is identically zero there, returned chi = 0 on every channel and
    # reported `converged = true`.  The module as shipped had therefore never produced a static screened
    # response, and nothing said so.  Raising is the right answer rather than warning: there is no useful
    # omega = 0 calculation with the dynamic drive, so a caller that asks for one has made a mistake.
    if  omega == 0.0  &&  !staticDrive
        error("RandomPhase.solve(): omega = 0 with staticDrive = false drives every channel with a kernel " *
              "that is identically zero at omega = 0, so the response would be zero for a reason that has " *
              "nothing to do with the physics. Pass staticDrive = true for the static case.")
    end
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
                                                         potential, primitives; settings=settings.dalgarnoLewis,
                                                         staticDrive=staticDrive)
        end

        return( CoupledResponse(multipole, gauge, omega, perturbedOrbitals, 1, true,
                                 Dict{Tuple{Subshell,Subshell},DalgarnoLewis.PerturbedOrbital}()) )
    end

    # THE COUPLED ITERATION. Start from the independent-particle chi's (iteration 0, no induced potential yet),
    # then repeatedly rebuild each electron's induced source from the OTHERS' current chi and re-solve.
    grid = primitives.grid
    # THE CHANNEL LIST.  With settings.allChannels the induced field is built from every target symmetry the
    # multipole can reach from each active subshell (`RandomPhase.reachableChannels`), which is what an RPA
    # asks for;  otherwise it is the single channel named in `targetSubshells`, which is the old behaviour and
    # remains the default.  The PRIMARY channel -- the one `targetSubshells` names -- is listed first for each
    # active subshell and is the one whose chi is returned, so `CoupledResponse` keeps its meaning either way.
    channels = Tuple{Subshell,Subshell}[]
    for  a in subshells
        push!(channels, (a, targets[a]))
        if  settings.allChannels
            for  t  in  RandomPhase.reachableChannels(multipole.L, a)
                # COMPARE ON kappa, NOT ON THE SUBSHELL.  `reachableChannels` labels each symmetry with the
                # lowest principal quantum number it admits, so the kappa = -1 channel comes back as 1s_1/2
                # while `targetSubshells` may name it 2s_1/2.  Those are unequal as Subshells and identical as
                # CHANNELS -- `DalgarnoLewis.solve` reads only the kappa -- so testing the Subshell would add
                # the primary channel a second time and double-count it in the induced field.
                t.kappa == targets[a].kappa   &&   continue
                push!(channels, (a, t))
            end
        end
    end
    if  settings.printBefore   &&   settings.allChannels
        @printf(">> RandomPhase.solve():  %d channel(s) over %d active subshell(s):  %s\n",
                length(channels), length(subshells),
                join([ "$(c[1])->$(c[2])" for c in channels ], ", "))
    end
    # TWO WARNINGS A CALLER MUST NOT BE ABLE TO MISS, and neither is suppressible by printBefore: both say the
    # returned number is untrustworthy, which is not a matter of verbosity.  `maxlog=1` because each states a
    # standing fact about the module rather than something about one call -- repeating it per call would bury
    # it in its own noise, which is how a warning stops being read.
    if  !settings.allChannels
        @warn("RandomPhase.solve(): allChannels = false keeps ONE target channel per active subshell, and that " *
              "truncation was measured on 14-Oct-2026 to give a SIGN-REVERSED total coherent response " *
              "(-0.793 of the uncoupled total, against +0.835 with every reachable channel). The spin-orbit " *
              "partners carry OPPOSITE angular signs, so a truncated set breaks a cancellation rather than " *
              "merely dropping a term. Use it for diagnostics, not for a physical number.", maxlog=1)
    end
    if  !settings.directSolve
        @warn("RandomPhase.solve(): directSolve = false iterates a map whose spectral radius was measured at " *
              "about 2.7, which belongs to the operator and not to the stepping -- so NO damping, Anderson " *
              "depth or channel choice converges, and whatever this call returns is an unconverged iterate. " *
              "Set directSolve = true.", maxlog=1)
    end

    chis = Dict{Tuple{Subshell,Subshell},DalgarnoLewis.PerturbedOrbital}()
    for  (a, t)  in  channels
        chis[(a, t)] = DalgarnoLewis.solve(multipole, gauge, omega, t, activeOrbitals[a], potential, primitives;
                                            settings=settings.dalgarnoLewis, staticDrive=staticDrive)
    end

    # THE DIRECT SOLVE.  The coupled map is AFFINE in the concatenated channel sources, F(s) = b + M s, so the
    # drive b is F(0) and the action of M is F(s) - b.  The fixed point therefore solves (I - M) s = b, which
    # a Krylov solver reaches whether or not the spectral radius of M is below one -- and it is about 2.7 on a
    # closed shell, which is why NO iteration of this map can work (see the long note below).  Nothing is
    # assembled: one operator application is one `RandomPhase.sourceMap` call.
    if  settings.directSolve
        m = (primitives.grid.nsL + primitives.grid.nsS) * length(channels)
        drive, _ = RandomPhase.sourceMap(multipole, gauge, omega, activeOrbitals, channels, zeros(m),
                                          potential, primitives, settings.dalgarnoLewis;
                                          occupiedOrbitals=occupiedOrbitals, occupations=occupations,
                                          staticDrive=staticDrive)
        applyIminusM = function(v::Array{Float64,1})
            mapped, _ = RandomPhase.sourceMap(multipole, gauge, omega, activeOrbitals, channels, v,
                                               potential, primitives, settings.dalgarnoLewis;
                                               occupiedOrbitals=occupiedOrbitals, occupations=occupations,
                                               staticDrive=staticDrive)
            return( v - (mapped - drive) )
        end
        sources, nKrylov, relResidual =
            RandomPhase.gmres(applyIminusM, drive, settings.maxIterations, settings.accuracy)
        # one last application turns the solved sources into the chis themselves
        checkSources, chis = RandomPhase.sourceMap(multipole, gauge, omega, activeOrbitals, channels, sources,
                                                    potential, primitives, settings.dalgarnoLewis;
                                                    occupiedOrbitals=occupiedOrbitals, occupations=occupations,
                                                    staticDrive=staticDrive)
        # the SELF-CONSISTENCY of the answer, measured on the sources rather than taken on trust from the solver
        selfConsistency = LinearAlgebra.norm(checkSources - sources) /
                          max(LinearAlgebra.norm(sources), 1.0e-300)
        if  settings.printBefore
            @printf(">> RandomPhase.solve():  DIRECT solve over %d channel(s), %d unknowns, %d Krylov step(s)\n",
                    length(channels), m, nKrylov)
            @printf(">>                       relative residual %.3e, induced-source self-consistency %.3e\n",
                    relResidual, selfConsistency)
        end
        perturbedOrbitals = Dict{Subshell,DalgarnoLewis.PerturbedOrbital}()
        for  a in subshells   perturbedOrbitals[a] = chis[(a, targets[a])]   end

        return( CoupledResponse(multipole, gauge, omega, perturbedOrbitals, nKrylov,
                                 relResidual < settings.accuracy  &&  selfConsistency < sqrt(settings.accuracy),
                                 chis) )
    end

    # THE ITERATION.  The fixed-point variable is the INDUCED SOURCE, one vector per active subshell, and the
    # map is "given the sources, solve for the chis;  from the chis, rebuild the sources".  Two schemes are
    # available and the choice is `settings.andersonDepth`.
    #
    # **PLAIN LINEAR MIXING HAS NO USABLE WINDOW, AND THAT WAS MEASURED RATHER THAN SUSPECTED (13-Oct-2026,
    # priority item 62).**  On Ne 1s^2 2s^2 2p^6 in its own DFS potential, E1 at omega = 0.1, three active
    # subshells, the ratio |chi|_coupled/|chi|_uncoupled after 40 iterations came out
    #       mixing 0.300 -> 9.4e+08     0.100 -> 2.3e+03     0.030 -> 7.8     0.010 -> 1.9     0.003 -> 1.0
    # and `converged` was false at every one:  DIVERGENT above about 0.03 and merely CREEPING below it.  The
    # ratio FALLING with the mixing is what rules out a sign error -- positive feedback would be
    # mixing-independent -- so the induced term is sound and it is the SCHEME that fails.
    #
    # **ANDERSON ACCELERATION WAS TRIED AND DOES NOT REPAIR IT -- MEASURED 13-Oct-2026, SAME CASE.**  It is the
    # construction `SelfConsistent.solveMeanFieldBasisAnderson` uses for the mean-field SCF:  keep a short
    # history of the iterate x and the residual f = x_new - x, solve gamma = dF \ f, step to
    #       x <- x + f - (dX + dF) gamma ,
    # with a rank-deficient window falling back to the plain step.  On the Ne case above, at 40 iterations:
    #       depth 2 -> 18.1 / 3.22 / 1.38      depth 3 -> 77.8 / 101 / 140      depth 5 -> 124 / 48.2 / 52.9
    # none converged, and ALL ARE WORSE THAN PLAIN MIXING AT 0.003.  So the default is andersonDepth = 0 and
    # the code is kept available rather than imposed.
    #
    # WHAT THAT TOGETHER IMPLIES, and it is the useful part:  a fixed point that neither damping nor Anderson
    # can reach is not a step-size problem at all.  Anderson accelerates a map whose iteration matrix has
    # spectral radius near or below one;  it cannot rescue one well above it.  So the question is WHY the
    # spectral radius is large, and that was measured next rather than argued about.
    #
    # **THE ANGULAR ALGEBRA WENT IN ON 14-Oct-2026, AND THIS NOTE RECORDS WHAT IT DID AND DID NOT FIX.**
    # Until then the module applied no angular coefficient anywhere, and weighted each contributing electron
    # by its full (2j_b+1).  Both are corrected -- a weight of 2*sqrt(2j_b+1) times the channel's reduced matrix
    # element at the density end, and `InteractionStrength.eMultipole`'s own reduction where the potential acts -- with
    # the derivation in `RandomPhase.inducedPotentialZr` and the two conventions in `RandomPhase.angularCoefficient`
    # and `RandomPhase.receivingCoefficient`.  **THE (-1)^L THIS NOTE CLAIMED IS GONE, 18-Oct-2026, item 75:** it
    # belonged to a two-electron product this formulation never assembles.
    #
    # THREE THINGS ARE NOW ESTABLISHED, each by a check that could have failed.
    #
    #   (1) THE RADIAL ROUTE IS JAC'S OWN, TO TEN DIGITS AND INCLUDING THE SIGN.  Compared entry by entry
    #       against `RadialIntegrals.SlaterRkComponent`, the object `InteractionStrength.XL_Coulomb`
    #       contracts: ratio 1.000000 on three spline pairs.  So nothing radial was ever wrong here.
    #
    #   (2) THE SIGN OF THE WHOLE ASSEMBLY IS RIGHT -- IT SCREENS.  rho depends on |eigenvalue| and is blind
    #       to this, so it was checked where the iteration converges (lambda = 1/8, residual 2.95e-05 ->
    #       4.46e-21 over 40 sweeps).  All three active subshells come out SMALLER than their uncoupled
    #       values -- 0.949, 0.982, 0.630 -- i.e. the induced field OPPOSES the external one, which is what a
    #       screening response must do and what a reversed sign would have shown as a ratio above 1.
    #
    #   (3) THE SPECTRAL RADIUS IMPROVED BY 2.17x BUT IS STILL ABOVE ONE:  rho_1 = 5.03 before, 2.32 after,
    #       both as the geometric mean of the per-iteration residual ratios (the ratios cycle, so a single
    #       one is not the number).  The undamped Jacobi iteration therefore still fails at full strength,
    #       and `couplingEnabled = true` still cannot be trusted for a converged number.
    #
    # **AND ONE CLAIM THIS NOTE PREVIOUSLY MADE IS WITHDRAWN:  THE CHANNEL TRUNCATION IS IMPLICATED AFTER
    # ALL.**  It had been argued that a narrower channel set could only make the induced term WEAKER, so the
    # one-target-per-subshell truncation could not explain an EXCESS.  That holds only while every
    # contribution carries the same sign, and the angular factors destroy exactly that premise: the two
    # spin-orbit partners of a subshell contribute with OPPOSITE signs -- CL(2p_1/2, 1, 2s_1/2) = -0.816
    # against CL(2p_3/2, 1, 2s_1/2) = +1.155 -- so dropping a channel removes a CANCELLATION and can make the
    # net term larger.  Measured, varying only which target each active subshell is given:
    #       2p1/2->2s, 2p3/2->2s, 2s->2p1/2   rho = 2.43
    #       2p1/2->2s, 2p3/2->2s, 2s->2p3/2   rho = 3.66
    #       only the two 2p partners active   rho = 1.94
    # Non-monotonic in the number of channels, which is the signature of cancellation rather than of a sum.
    #
    # **AND BOTH REMAINING CANDIDATES ARE NOW ELIMINATED -- 14-Oct-2026, SAME DAY.**  The angular repair left
    # rho at 2.32 and named three things that might account for the rest.  All three have been measured.
    #
    #   (a) THE CHANNEL TRUNCATION DOES NOT EXPLAIN IT.  `settings.allChannels` widens the Ne case from 3
    #       channels to 7 -- 2p_1/2 -> s_1/2, d_3/2;  2p_3/2 -> s_1/2, d_5/2, d_3/2;  2s_1/2 -> p_1/2, p_3/2 --
    #       and rho does NOT fall:  2.33 on three channels against 2.49 on seven.  So the cancellation between
    #       the spin-orbit partners is real (rho moves 1.94-3.66 on WHICH single target is chosen) but
    #       restoring all of it does not bring the iteration into range.  The hypothesis is dead; the setting
    #       stays, because the 7-channel field is the physically correct one regardless.
    #
    #   (b) THE FACTOR 2 DOES NOT EXPLAIN IT EITHER, and the test is clean because that factor becomes EXACT
    #       at omega = 0.  Sweeping omega down with all seven channels:
    #           omega   0.100   0.030   0.010   0.003   0.001
    #           rho     2.433   2.623   2.690   2.714   2.721
    #       rho RISES slightly and saturates rather than falling towards one, so at the point where the
    #       symmetrization is no longer an approximation the iteration is at its WORST.  The factor 2 is
    #       therefore not the defect, and if anything it slightly understates the term.
    #
    #   (c) WHAT IS LEFT IS THE APPROXIMATION ITSELF:  ONLY THE DIRECT (ring) TERM IS BUILT.  A direct-only
    #       induced field is known to OVER-screen, because exchange opposes it, and rho ~ 2.7 is what that
    #       looks like.  This is a property of "RPA without exchange" for a closed shell, not a defect in the
    #       code -- and it is omega-independent in character, which is exactly what (b) measured.
    #
    # **SO TWO CONCLUSIONS, AND THE SECOND MATTERS MORE THAN THE FIRST.**
    #
    # FIRST, A JACOBI ITERATION CANNOT WORK HERE AND A DIRECT SOLVE IS THE RIGHT METHOD.  rho = 2.7 is a
    # property of the operator, not of the stepping, so no damping, acceleration or channel choice will reach
    # the fixed point -- while (I - J) remains invertible and a direct solve would find it in one step.  This
    # RESTORES the conclusion of 13-Oct that an earlier version of this note withdrew.  The withdrawal was
    # wrong for a reason worth stating: it assumed a correctly-built induced term for a closed shell would
    # have rho ~ 0.3, extrapolating from the tens-of-percent screening the FULL RPA is known to give.  That
    # extrapolation does not survive the direct-only approximation, and nothing had measured it.
    #
    # SECOND, AND IT IS A QUESTION OF SCOPE RATHER THAN OF CODE:  rho ~ 2.7 means this term OVER-SCREENS
    # badly, so even solving it exactly would not yield a trustworthy screening correction for a closed shell.
    # Exchange is not a small correction in this regime, and the module's own scope note puts it out of
    # reach ("RPA", not "RPA with exchange").  **So a direct solve makes the module converge without making
    # it right**, and that remains true now that the direct solve exists.
    #
    # **THE DIRECT SOLVE IS BUILT AND VERIFIED -- 14-Oct-2026, `settings.directSolve`.**  It is GMRES over
    # `RandomPhase.sourceMap`, so no matrix is assembled, and it converges in about ten Krylov steps where no
    # iteration converges at all:  11 steps on 1395 unknowns (3 channels) and 10 on 3255 (all 7), with
    # relative residuals of 1.7e-12 and 8.7e-11 and the self-consistency of the induced source independently
    # re-measured at 4.8e-12 and 1.7e-10 by applying the map once more.
    #
    # **IT IS VERIFIED AGAINST THE ITERATION IN THE ONE REGIME WHERE BOTH WORK, AND THAT IS THE WHOLE PROOF.**
    # At 1/8 coupling the plain iteration converges (rho = 0.456), so an independent answer exists.  The two
    # methods agree to THIRTEEN DIGITS there -- worst |direct/iteration - 1| = 1.1e-13 over the three channels,
    # with both giving 0.94915 / 0.98168 / 0.63017 of the uncoupled response.  At full strength the iteration
    # has overflowed to 1e18 while the direct solve returns finite, self-consistent values in 12 steps.
    #
    # AND WHAT THE CONVERGED ANSWER SHOWS, which is the point of having it:  at full coupling the response is
    # NOT a screened version of the uncoupled one.  The 2p_1/2 channel comes out ENHANCED, by 1.22x over three
    # channels and 1.82x over seven, where at 1/8 coupling every channel was screened.  That is what rho ~ 2.7
    # means from the inside -- the induced term is comparable to the response rather than a correction to it.
    # **A per-channel |chi| ratio is therefore NOT a screening diagnostic at full strength**:  screening is a
    # statement about the total coherent response, and the channels mix strongly here.  Read a coupled number
    # from this module as a direct-only result in a strongly-coupled regime, not as a screening correction.
    converged = false;   nIterations = 0
    previousExtraSource = Dict{Tuple{Subshell,Subshell},Array{Float64,1}}()
    xHistory = Dict{Tuple{Subshell,Subshell},Vector{Vector{Float64}}}()
    fHistory = Dict{Tuple{Subshell,Subshell},Vector{Vector{Float64}}}()
    for  ch in channels   xHistory[ch] = Vector{Float64}[];   fHistory[ch] = Vector{Float64}[]   end
    for  iteration = 1:settings.maxIterations
        nIterations = iteration
        newChis = Dict{Tuple{Subshell,Subshell},DalgarnoLewis.PerturbedOrbital}();   maxRelChange = 0.0
        for  ch in channels
            a = ch[1];   t = ch[2]
            VindZr      = RandomPhase.inducedPotentialZr(multipole, activeOrbitals, chis, a, grid;
                                                      occupiedOrbitals=occupiedOrbitals, occupations=occupations)
            # THE a-END ANGULAR FACTOR.  `inducedSourceVector` returns the bare radial projection;  what belongs
            # here is the reduction `InteractionStrength.eMultipole` applies to the EXTERNAL drive of this very
            # equation, because the induced potential is a one-body potential exactly as the drive is -- see
            # `RandomPhase.receivingCoefficient`, which on 18-Oct-2026 (item 75) replaced a target-first
            # `angularCoefficient` and a (-1)^L here.  It is THIS channel's own t, not targets[a], since one active
            # subshell may carry several channels.
            aCoeff      = RandomPhase.receivingCoefficient(multipole.L, a, t)
            justComputed = aCoeff * RandomPhase.inducedSourceVector(VindZr, activeOrbitals[a], primitives)
            priorSource  = get(previousExtraSource, ch, zeros(length(justComputed)))
            if  settings.andersonDepth <= 0  ||  iteration == 1
                dampedSource = settings.mixing * justComputed + (1.0 - settings.mixing) * priorSource
            else
                # x is the source actually used last time, f its residual against what the chis then implied
                fNew = justComputed - priorSource
                push!(xHistory[ch], copy(priorSource));   push!(fHistory[ch], copy(fNew))
                while  length(xHistory[ch]) > settings.andersonDepth + 1
                    popfirst!(xHistory[ch]);   popfirst!(fHistory[ch])
                end
                m = length(xHistory[ch]) - 1
                if  m < 1
                    dampedSource = priorSource + fNew
                else
                    dF = zeros(length(fNew), m);   dX = zeros(length(fNew), m)
                    for  j = 1:m
                        dF[:,j] = fHistory[ch][j+1] - fHistory[ch][j];   dX[:,j] = xHistory[ch][j+1] - xHistory[ch][j]
                    end
                    local gamma
                    try
                        gamma = dF \ fNew
                    catch
                        gamma = zeros(m)            ## a rank-deficient window: fall back to the plain step
                    end
                    if  any(!isfinite, gamma)   gamma = zeros(m)   end
                    dampedSource = priorSource + fNew - (dX + dF) * gamma
                    if  any(!isfinite, dampedSource)   dampedSource = priorSource + settings.mixing * fNew   end
                end
            end

            newChis[ch] = DalgarnoLewis.solve(multipole, gauge, omega, t, activeOrbitals[a], potential,
                                               primitives; settings=settings.dalgarnoLewis,
                                               extraSource=dampedSource, staticDrive=staticDrive)

            denom = LinearAlgebra.norm(priorSource)
            # denom == 0 on the FIRST iteration just means "nothing to compare against yet" (Inf, keep
            # iterating) -- UNLESS dampedSource is ALSO zero (e.g. a single active electron, which has no
            # coupling partner at all and so no induced source ever), which already IS the fixed point.
            relChange = denom > 0.0 ? LinearAlgebra.norm(dampedSource - priorSource) / denom :
                        (LinearAlgebra.norm(dampedSource) > 0.0 ? Inf : 0.0)
            maxRelChange = max(maxRelChange, relChange)
            previousExtraSource[ch] = dampedSource
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

    # BACK TO ONE chi PER ACTIVE SUBSHELL.  `CoupledResponse` carries the PRIMARY channel's perturbed orbital,
    # the one `targetSubshells` names, so its type and meaning are the same whether or not allChannels was set.
    # The extra channels are solved and do enter the induced field;  they are not returned, because an active
    # subshell's full first-order change is a SUM over its channels and that sum is not a PerturbedOrbital.
    perturbedOrbitals = Dict{Subshell,DalgarnoLewis.PerturbedOrbital}()
    for  a in subshells   perturbedOrbitals[a] = chis[(a, targets[a])]   end

    return( CoupledResponse(multipole, gauge, omega, perturbedOrbitals, nIterations, converged, chis) )
end


end # module
