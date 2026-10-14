
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
                                                         The default is false, so an existing caller is
                                                         unaffected.
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
end


"""
`RandomPhase.Settings()`  ... constructor for an `empty` instance of RandomPhase.Settings.
"""
function Settings()
    Settings( Subshell[], Dict{Subshell,Subshell}(), false, 15, 1.0e-6, 0.3, DalgarnoLewis.Settings(), false, 0, false )
end


"""
`RandomPhase.Settings(settings::RandomPhase.Settings; activeSubshells=.., targetSubshells=.., couplingEnabled=..,
                       maxIterations=.., accuracy=.., mixing=.., dalgarnoLewis=.., printBefore=.., andersonDepth=..,
                       allChannels=..)`
    ... the standard JAC keyword copy-constructor.
"""
function Settings(settings::RandomPhase.Settings;
    activeSubshells::Union{Nothing,Array{Subshell,1}}=nothing,   targetSubshells::Union{Nothing,Dict{Subshell,Subshell}}=nothing,
    couplingEnabled::Union{Nothing,Bool}=nothing,                 maxIterations::Union{Nothing,Int64}=nothing,
    accuracy::Union{Nothing,Float64}=nothing,                     mixing::Union{Nothing,Float64}=nothing,
    dalgarnoLewis::Union{Nothing,DalgarnoLewis.Settings}=nothing, printBefore::Union{Nothing,Bool}=nothing,
    andersonDepth::Union{Nothing,Int64}=nothing,                  allChannels::Union{Nothing,Bool}=nothing)

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

    Settings( activeSubshellsx, targetSubshellsx, couplingEnabledx, maxIterationsx, accuracyx, mixingx, dalgarnoLewisx, printBeforex,
              andersonDepthx, allChannelsx )
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
`RandomPhase.angularCoefficient(L::Int64, tSub::Subshell, aSub::Subshell)`
    ... the one-electron angular factor the rank-L induced interaction carries on the channel aSub -> tSub,
        namely `AngularMomentum.CL_reduced_me(tSub, L, aSub)`, or 0.0 where the triangular and parity
        selection rules forbid the channel. It exists as a named function, rather than inline, because the
        convention is the thing most easily got wrong here and Rule 18 asks for it to be visible.

        **WHICH CONVENTION, AND HOW IT WAS FIXED -- 14-Oct-2026.** `DalgarnoLewis.solve` builds its EXTERNAL
        source from `DalgarnoLewis.amplitudeKernel`, i.e. from a REDUCED one-electron matrix element on the
        channel aSub -> tSub. An induced source added to that same equation must therefore be reduced the
        same way, which is exactly one factor of `CL_reduced_me` per end and no division by sqrt(2j+1)
        anywhere. The argument ORDER matters and is not free: `CL_reduced_me(tSub, L, aSub)` and
        `CL_reduced_me(aSub, L, tSub)` differ by (-1)^(j_t - j_a), so they have OPPOSITE signs for
        2s_1/2 <-> 2p_3/2 and the same sign for 2s_1/2 <-> 2p_1/2. The order used here is the one
        `InteractionStrength.XL_Coulomb(L, a::Subshell, b, c::Subshell, d, primitives)` uses for the
        corresponding two-electron strength -- target (row) first, source (column) second -- so that this
        module's radial route and JAC's own already-validated angular factor agree term by term.

        **THE `(-1)^L` IS NOT APPLIED HERE.** `XL_Coulomb` carries it once for the whole two-electron
        product, and this function is called at BOTH ends -- once where the induced density is formed and
        once where it acts -- so applying it here would square it away. `RandomPhase.solve` applies it once,
        and says so at the line that does.
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
`RandomPhase.inducedPotentialZr(multipole::EmMultipole, activeOrbitals::Dict{Subshell,Orbital},
                                 chis::Dict{Tuple{Subshell,Subshell},DalgarnoLewis.PerturbedOrbital},
                                 excludeSubshell::Subshell, grid::Radial.Grid)`
    ... builds the DIRECT (Hartree, "ring") induced potential electron `excludeSubshell` feels from every OTHER
        active electron's CURRENT chi, at the multipole's own rank L, in the SAME "Zr" (= -r*V(r)) convention
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
            be reduced to match, one factor of `CL_reduced_me` per end (see `RandomPhase.angularCoefficient`);
          * and `AngularMomentum.CL_reduced_me(b, 0, b) = sqrt(2j_b+1)` EXACTLY, checked for s, p, d and f.

        Equating the two therefore gives weight * sqrt(2j_b+1) = 2*(2j_b+1), i.e. a weight of
        2*sqrt(2j_b+1): the occupation splits in half, one sqrt(2j_b+1) staying here and the other being
        carried by the reduced matrix element's own normalisation. This is precisely the sqrt(2j+1) Rule 18
        is about, and it is written out because the two quantities it sits between -- an ordinary potential
        and a reduced source -- are otherwise denoted alike.

        **AND THE ANGULAR FACTOR CHANGES THE PHYSICS, NOT ONLY THE SCALE.** `CL_reduced_me(b, L, t_b)` carries
        a SIGN, and for the two spin-orbit partners of a subshell those signs are OPPOSITE: for L = 1,
        CL(2p_1/2, 1, 2s_1/2) = -0.816 against CL(2p_3/2, 1, 2s_1/2) = +1.155. So the partners partially
        cancel in the induced field, where the old weighting added them with one sign and with weights 2 and
        4. That cancellation is physical and is a large part of why the uncorrected term was far too strong.

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
                             excludeSubshell::Subshell, grid::Radial.Grid)
    rho = zeros(grid.NoPoints);   mtpRho = 1
    # EVERY CHANNEL OF EVERY OTHER ELECTRON CONTRIBUTES, EACH WITH ITS OWN ANGULAR FACTOR.  `chis` is keyed
    # by the (active, target) PAIR, so the target is read off the key rather than looked up -- which is what
    # lets one electron carry several channels, as an RPA requires.  Channels of the excluded electron are
    # skipped whatever their target: b is excluded from its own induced field, not merely one channel of it.
    for  ((b, tb), chiB)  in  chis
        b == excludeSubshell   &&   continue
        orbB   = activeOrbitals[b];   jb2 = Basics.subshell_2j(b)
        # The weight is sqrt(2j_b+1), NOT (2j_b+1):  see this method's docstring for the Hartree calibration
        # that fixes it.  The remaining sqrt(2j_b+1) of the occupation is carried by CL_reduced_me's own
        # normalisation at this end, since CL_reduced_me(b, 0, b) = sqrt(2j_b+1) exactly.
        weightB = sqrt(jb2 + 1.0) * RandomPhase.angularCoefficient(multipole.L, b, tb)
        weightB == 0.0   &&   continue
        mtpB = min(length(orbB.P), length(chiB.P), length(orbB.Q), length(chiB.Q))
        for  i = 1:mtpB   rho[i] = rho[i] + weightB * 2.0 * (orbB.P[i]*chiB.P[i] + orbB.Q[i]*chiB.Q[i])   end
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

        **IT CARRIES NO ANGULAR COEFFICIENT, DELIBERATELY, AND THE CALLER MUST SUPPLY ONE.** What this
        function returns is the bare radial projection; the active channel's own reduced matrix element of
        C^L, and the `(-1)^L` of the two-electron product, are applied by `RandomPhase.solve` through
        `RandomPhase.angularCoefficient`. Until 14-Oct-2026 nothing applied them at all, which is what made
        the coupled iteration diverge. A caller that uses this function directly, as a diagnostic does, must
        multiply by `angularCoefficient(L, target, active)` and the sign itself -- see `RandomPhase.solve`.
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

    chis = Dict{Tuple{Subshell,Subshell},DalgarnoLewis.PerturbedOrbital}()
    for  (a, t)  in  channels
        chis[(a, t)] = DalgarnoLewis.solve(multipole, gauge, omega, t, activeOrbitals[a], potential, primitives;
                                            settings=settings.dalgarnoLewis)
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
    # by its full (2j_b+1).  Both are now corrected -- one reduced matrix element of C^L per end, a weight of
    # 2*sqrt(2j_b+1), and the (-1)^L of the two-electron product applied once -- with the derivation in
    # `RandomPhase.inducedPotentialZr` and the convention in `RandomPhase.angularCoefficient`.
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
    # reach ("RPA", not "RPA with exchange").  **So a direct solve would make the module converge without
    # making it right.**  Whether to build it anyway, to extend to exchange, or to leave the coupled path
    # as a documented demonstration is the maintainer's decision and is not taken here.
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
            VindZr      = RandomPhase.inducedPotentialZr(multipole, activeOrbitals, chis, a, grid)
            # THE a-END ANGULAR FACTOR AND THE SIGN.  `inducedSourceVector` returns the bare radial
            # projection;  THIS channel's own reduced matrix element belongs here -- t, not targets[a], since
            # one active subshell may carry several channels -- and the (-1)^L that
            # `InteractionStrength.XL_Coulomb` carries for the two-electron product is applied ONCE, here,
            # rather than at both ends (see `RandomPhase.angularCoefficient`).
            aCoeff      = RandomPhase.angularCoefficient(multipole.L, t, a)
            if  rem(multipole.L, 2) == 1    aCoeff = -aCoeff    end
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
                                               primitives; settings=settings.dalgarnoLewis, extraSource=dampedSource)

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

    return( CoupledResponse(multipole, gauge, omega, perturbedOrbitals, nIterations, converged) )
end


end # module
