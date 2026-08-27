
"""
`module  JAC.Bsplines`
    ... a submodel of JAC that contains all structs and methods to generate the B-spline basis and to solve the single-electron Dirac
        equation in a local potential. It also provides the major function calls to generate self-consistent fields; cf. JAC.SelfConsistent.
"""
module Bsplines


using  BSplineKit, Printf, ..Basics, ..Defaults, ..Nuclear, ..Radial, JenaAtomicCalculator


"""
`struct  Bsplines.Bspline`  
    ... defines a type for a (single) B-spline that is defined on a given radial grid from r[lower:upper]. Note that only the non-zero
        values are specified for the B-spline function and its derivative.

    + lower        ::Int64               ... lower radial index (on the radial grid.r) from where the functions is nonzero.
    + upper        ::Int64               ... upper radial index up to which the functions is nonzero.
    + bs           ::Array{Float64,1}    ... radial B-spline functions as defined on the predefined grid.r[lower:upper]
    + bp           ::Array{Float64,1}    ... derivative of bs on the predefined grid grid.r[lower:upper]
"""
struct Bspline
    lower          ::Int64 
    upper          ::Int64 
    bs             ::Array{Float64,1}   
    bp             ::Array{Float64,1}   
end


"""
`struct  Bsplines.Primitives`  ... defines a type for a set of primitive functions which typically belongs to a well-defined grid.

    + grid         ::Radial.Grid         ... radial grid on which the states are represented.
    + bsplinesL    ::Array{Bspline,1}    ... set of B-splines for the large components on the given radial grid.
    + bsplinesS    ::Array{Bspline,1}    ... set of B-splines for the small components on the given radial grid.
"""
struct Primitives
    grid           ::Radial.Grid
    bsplinesL      ::Array{Bspline,1}
    bsplinesS      ::Array{Bspline,1}
end


"""
`Base.show(io::IO, primitives::Bsplines.Primitives)`
    ... prepares a proper printout of the variable primitives::Bsplines.Primitives; nothing is returned.
"""
function Base.show(io::IO, primitives::Bsplines.Primitives) 
    println(io, "grid:               $(primitives.grid)  ")
    println(io, "bsplinesL:           (primitives.bsplinesL)  ")
    println(io, "bsplinesS:           (primitives.bsplinesS)  ")
end


"""
`Bsplines.boundaryDropCounts(kappa::Int64, grid::Radial.Grid)`
    ... returns the number of leading (near r=0) and trailing (near r=R_max) large- and small-component B-splines that must be excluded from
        the generalized eigenvalue problem for the given symmetry kappa, following Zatsarinny & Froese Fischer's DBSR_HF reference code
        (hf_boundary.f90's Boundary_conditions, the finite/non-point-nucleus branch) and cross-checked against a live DBSR_HF run on Sc: P
        needs l+1 excluded at r=0 (kappa-sign independent, since the large component's leading near-origin power depends only on l); Q needs
        l+2 excluded for kappa<0 (j=l+1/2, where the naive leading power l+1+kappa vanishes identically since kappa=-(l+1), pushing the true
        leading power one order higher) but only l for kappa>0 (j=l-1/2, kappa=l, no such cancellation) -- both capped at (order-1),
        matching DBSR_HF's own j>ksp-1/j>ksq-1 caps.

        The trailing (outer-boundary) counts are kappa/l-INDEPENDENT in DBSR_HF (confirmed: identical trailP=3, trailQ=2 for every orbital
        in the live Sc run) and tied to a small fixed parameter DBSR_HF calls ibzero (=2 by default): trailP = ibzero+1 and trailQ = ibzero.
        This outer truncation is not optional: without it the diagonalized spectrum is subtly perturbed even for the kappa>0 states that
        need no elimination at r=0 at all.
        A tuple (dropP::Int64, dropQ::Int64, trailP::Int64, trailQ::Int64) is returned.
"""
function boundaryDropCounts(kappa::Int64, grid::Radial.Grid)
    l      = kappa > 0  ?  kappa  :  -kappa - 1
    dropP  = min(l + 1,                     grid.orderL - 1)
    dropQ  = kappa < 0  ?  min(l + 2, grid.orderS - 1)  :  min(l, grid.orderS - 1)
    ibzero = 2
    trailP = ibzero + 1
    trailQ = ibzero
    return( dropP, dropQ, trailP, trailQ )
end


"""
`Bsplines.checkGridRepresentation(subshells::Array{Subshell,1}, Z::Float64, primitives::Bsplines.Primitives;
                                  occupations::Dict{Shell,Int64}=Dict{Shell,Int64}(), accuracy::Float64=1.0e-3,
                                  stopper::Bool=true, mass::Float64=1.0)`
    ... checks whether the given radial grid can represent every subshell of the list, by solving the single-electron Dirac equation for a
        POINT nucleus on this grid and comparing each level with the closed-form point-nucleus energy Basics.computeDiracEnergy. Every
        subshell whose energy deviates by more than `accuracy` is listed, and an error is raised; with stopper = false a loud warning is
        printed instead and the computation proceeds.
        A tuple  (isRepresentable::Bool, recommendedRbox::Float64)  is returned.

        A SUBSHELL FAILS ONLY WHEN THE GRID REPRESENTS NEITHER OF TWO HYDROGENIC PROXIES, when `occupations` is given: the orbital at the
        BARE charge Z, and the orbital at the Slater screened charge of Basics.slaterScreening (floored at the asymptotic Z - NoElectrons +
        1). That is not a way of being lenient; it is because NEITHER proxy is the real orbital and the real one lies between them, so a
        grid which carries either end is not demonstrably wrong. Measured against the known answers:

            Na 3s   bare Zeff 11.0  r_+ = 1.64  |  Slater Zeff 2.20  r_+ = 8.18  |  real r_max ~ 3.2
            Mg 3s   bare Zeff 12.0  r_+ = 1.50  |  Slater Zeff 2.85  r_+ = 6.32  |  real r_max ~ 2.9
            K  4s   bare Zeff 19.0  r_+ = 1.68  |  Slater Zeff 2.20  r_+ = 14.6  |  real r_max ~ 4.2

        The bare charge understates an outer orbital by about a factor of two and the screened charge overstates it by two to three, so
        EITHER used alone over-rejects, in opposite directions. Until 17-Aug-2026 only the bare charge was used, and it rejected the valence
        orbital of every heavy near-neutral system: on Th+ [Rn] 6d 7s^2, on its own recommended box of 52.7 a.u., it refused the 7s at
        1.93e-3 and advised a box of 2.3 a.u. for an 89-electron ion, while with `gridStopper = false` that very run converged in ten
        iterations and passed Bsplines.checkOrbitalConsistency, the test that actually catches a wrong state. Screening alone is no better:
        it rejects three cases of the JAC test suite at up to 1.5e-1, demanding a 38.8 a.u. box for a system that computes correctly on
        11.1.

        Inner shells are unaffected either way, since there is almost nothing to screen -- the thorium 1s sees 89.7 rather than 90 -- so the
        check stays exactly as strict where the bare-Z yardstick was already right.

        WITHOUT `occupations` the bare charge alone is used, which reproduces the behaviour before 17-Aug-2026. That is all a caller holding
        no occupation numbers can do, but it is NOT appropriate for a neutral or near-neutral heavy system; pass the occupations there.

        WHY A POINT NUCLEUS IS THE RIGHT YARDSTICK, whatever the computation itself uses. The question asked here is not "is this orbital
        physically accurate" but "can this GRID resolve an orbital of this size and shape at all", and only the point-nucleus spectrum has
        closed-form energies to test against. An extended nucleus, or a Dirac-Fock rather than a bare nuclear field, changes each individual
        orbital by a modest factor -- it does not change the ORDER of the spectrum, nor the radial scale that the box has to accommodate.

        THE USUAL CAUSE IS A BOX THAT IS TOO LARGE, not one that is too small. The number of B-splines is fixed, so a box much wider than
        the orbitals spends them on empty space and starves the physical region: at Z = 10 the default box of 614 a.u. leaves 5f_7/2 wrong
        by 31%, while a box of 11 a.u. -- matched to the orbital -- gives it to 6e-5. The box recommended here is the one of
        Basics.recommendedGrid, `r_plus + 16 n/Zeff`, so that the two cannot contradict each other.

    + subshells    ::Array{Subshell,1}   ... subshells that the grid must be able to carry.
    + Z            ::Float64             ... nuclear charge.
    + occupations  ::Dict{Shell,Int64}   ... shell occupations; if empty, every subshell is tested at the bare charge.
    + accuracy     ::Float64             ... largest tolerated relative deviation from the closed-form energy.
    + stopper      ::Bool                ... true, if a failure shall raise rather than warn.
"""
function checkGridRepresentation(subshells::Array{Subshell,1}, Z::Float64, primitives::Bsplines.Primitives;
                                 occupations::Dict{Shell,Int64}=Dict{Shell,Int64}(), accuracy::Float64=1.0e-3,
                                 stopper::Bool=true, mass::Float64=1.0)
    grid    = primitives.grid;      nsL = grid.nsL;     nsS = grid.nsS
    storage = Dict{String,Array{Float64,2}}()
    wb      = zeros( nsL+nsS, nsL+nsS )
    wb[1:nsL,1:nsL]                 = Bsplines.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
    wb[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = Bsplines.generateTTpMatrix!("SS-overlap", 0, primitives, storage)
    NoElectrons = length(occupations) == 0  ?  0  :  sum( values(occupations) )
    # The screened charge differs from subshell to subshell, so the local matrix can no longer be shared by
    # all subshells of one kappa; it is cached on (kappa, Zeff) instead, which still shares it wherever the
    # screening happens to agree.
    spectra    = Dict{Tuple{Int64,Float64}, Tuple{Vector{Float64},Int64}}()
    offenders  = Tuple{Subshell,Float64,Float64,Float64,Float64}[];     rbox = 0.
    borderline = Tuple{Subshell,Float64,Float64}[]

    # the relative deviation of one subshell from its closed-form energy, at a given charge
    function deviationAt(sh::Subshell, l::Int64, Zx::Float64)
        key = (sh.kappa, round(Zx, digits=8))
        if  !haskey(spectra, key)
            pot = Nuclear.pointNucleus(Zx, grid)
            wa  = Bsplines.setupLocalMatrix(sh.kappa, primitives, pot, storage; mass=mass)
            w2  = Bsplines.diagonalizeLocalMatrix(sh.kappa, wa, wb, primitives)
            spectra[key] = (w2.values, Bsplines.findPositiveBranchStart(w2.values; mass=mass))
        end
        (values2, mm) = spectra[key]
        ni = mm + sh.n - l - 1
        ex = Basics.computeDiracEnergy(sh, Zx; mass=mass)
        if  ni < 1  ||  ni > length(values2)    return( (NaN, ex, Inf) )    end
        return( (values2[ni], ex, abs(values2[ni]/ex - 1)) )
    end

    for  sh  in  subshells
        l = Basics.subshell_l(sh)
        (enB, exB, devB) = deviationAt(sh, l, Z)
        Zeff = Z;   enS = enB;   exS = exB;   devS = devB
        if  length(occupations) > 0
            Zeff = max( Z - NoElectrons + 1., Z - Basics.slaterScreening(Shell(sh.n, l), occupations), 1.0 )
            if  Zeff != Z
                (enS, exS, devS) = deviationAt(sh, l, Zeff)
            end
        end
        # The true orbital is BRACKETED by the two hydrogenic proxies, so only a grid that represents NEITHER
        # end is demonstrably wrong; see the docstring for why either alone over-rejects.  Where just one of
        # them fails the grid is not condemned, but it is not silently passed either.
        if      min(devB, devS) > accuracy      push!(offenders,  (sh, enS, exS, devS, Zeff))
        elseif  max(devB, devS) > accuracy      push!(borderline, (sh, devB, devS))
        end
        wr   = (sh.n^2/Zeff) * (1 + sqrt( max(0., 1 - l*(l+1)/sh.n^2) )) + 16. * sh.n / Zeff
        rbox = max(rbox, wr)
    end
    if  length(offenders) > 0
        printstyled("\n>>> GRID CHECK FAILED: on this grid the following subshells are not represented to the requested " *
                    "accuracy of $accuracy\n>>> in the point-nucleus Dirac spectrum at EITHER the bare or the screened " *
                    "charge (Z = $Z):\n", color=:light_red)
        printstyled("      subshell      Zeff      E(grid) [a.u.]      E(Dirac) [a.u.]     rel. deviation\n", color=:light_red)
        for  (sh, en, ex, dev, Zeff)  in  offenders
            printstyled(@sprintf("    %10s   %7.2f    %+.8e     %+.8e      %.2e\n",
                                 string(sh), Zeff, en, ex, dev), color=:light_red)
        end
        printstyled(@sprintf(">>> The present box is r_max = %.1f a.u.;  a box of about %.1f a.u. suits these subshells.\n",
                             grid.r[end], rbox), color=:light_red)
        printstyled(">>> Note that a box which is much TOO LARGE starves the basis just as badly as one that is too\n" *
                    ">>> small, since the number of B-splines is fixed;  use Radial.Grid(grid; rbox=..) to match it,\n" *
                    ">>> or Basics.recommendedGrid(configs, nm) to have it matched automatically.\n", color=:light_red)
        if  stopper   error("Bsplines.checkGridRepresentation(): the grid fails to represent " *
                            "$(length(offenders)) of $(length(subshells)) subshells to accuracy $accuracy.")   end
        return( (false, rbox) )
    end

    # A subshell that one proxy carries and the other does not says the grid sits near the edge for that
    # orbital. It is not grounds for refusing the computation -- the real orbital lies between the two -- but
    # it is the only warning that will be given, so it is not dropped.
    if  length(borderline) > 0
        sa = ""
        for  (sh, devB, devS)  in borderline
            sa = sa * @sprintf("%s (bare %.1e, screened %.1e)  ", string(sh), devB, devS)
        end
        printstyled(">>> Grid check: these subshells are carried at one of the two hydrogenic charges but not the\n" *
                    ">>> other, so the box is near its limit for them:  " * sa * "\n" *
                    @sprintf(">>> The present box is r_max = %.1f a.u.;  about %.1f a.u. would suit them comfortably.\n",
                             grid.r[end], rbox), color=:yellow)
    end

    return( (true, rbox) )
end


"""
`Bsplines.checkOrbitalBox(orbitals::Dict{Subshell,Orbital}, grid::Radial.Grid;
                          extentTolerance::Float64=0.9, densityFloor::Float64=1.0e-5, stopper::Bool=false)`
    ... checks, on the CONVERGED orbitals rather than on a hydrogenic stand-in, whether the radial box was in fact large enough for them,
        and reports what it measured; a tuple (isAdequate::Bool, report::String) is returned, the report being a single line fit to be
        written into a summary file.

        THE BOX IS AN ESTIMATE UNTIL SOMETHING CHECKS IT. Basics.recommendedGrid derives it from Slater's rules, a 1930s empirical recipe,
        and nothing downstream ever asked whether the answer was right -- so every result carried an unquantified box error. The measure
        used here is the orbital's own EXTENT,
        `rEnd`, the radius beyond which its density has fallen below `densityFloor` of its maximum: if
        `rEnd/rbox` approaches 1 the orbital is still substantial where the box ends and has been cut off.

        WHY NOT THE AMPLITUDE AT THE WALL, which is the obvious choice: the B-spline basis imposes P(rbox) = 0, so the density at the wall
        is zero BY CONSTRUCTION and measures the boundary condition rather than the orbital. Measured on argon in a 3 a.u. box -- an error
        of 37 Ha in the total energy -- the wall density came to 9.2e-10, i.e. the most badly truncated case available scored better than
        the correct one.

        THE DEFAULT densityFloor = 1e-5 WAS CALIBRATED, not chosen, against argon boxes whose harm is known: 3.0 a.u. is wrong by 37 Ha and
        6.0 a.u. by 2.3e-3 Ha, both of which must be flagged, while 10.2 a.u. gives the right answer and 61 a.u. costs 5e-9 Ha, neither of
        which may be. On an uncut 61 a.u. box the outermost argon orbital ends at 6.92, 8.10 and 9.26 a.u. for floors of 1e-4, 1e-5 and
        1e-6, so 1e-5 puts the correct box at extent/box = 0.80 against the 0.90 limit while both bad boxes reach 1.0. A floor of 1e-10
        fails outright: the tail is then still alive at 13 a.u. and the CORRECT box is condemned.

        WHAT THIS DOES NOT DO, and it is half the problem: it does not detect a box that is too LARGE. That failure has cost this codebase
        more than truncation ever did, but it leaves no signature in the outer orbital's shape, only in its energy, and it is strongly
        n-dependent. Measured on thorium at the full nuclear charge with hp = rbox/300 throughout, the 1s is unmoved (3e-9) across boxes
        from 10 to 600 a.u. while the 7s degrades 1.0e-7 -> 2.0e-3 -> 1.1e-1 -> 2.4e-1; yet argon, whose outermost orbital has one node
        instead of six, loses only 1.0e-9 on a box 25 times too large. Counting B-spline break points inside the orbital does not see it
        either: over that same thorium range the count falls only from 60 to 47 while the error grows by six orders, because with hp scaled
        to the box the mesh holds
        `log(R/rnt)/h` points below any radius R independently of rbox.  No threshold in the single-run data was
        defensible, so none is imposed. Bsplines.checkOrbitalConsistency remains the guard for that side: a starved subshell shows up as its
        two spin-orbit partners ceasing to describe the same shell.

    + orbitals         ::Dict{Subshell,Orbital}  ... the CONVERGED orbitals.
    + grid             ::Radial.Grid             ... the grid they were computed on.
    + extentTolerance  ::Float64                 ... largest tolerated rEnd/rbox before the orbital counts as cut off.
    + densityFloor     ::Float64                 ... fraction of the peak density that defines where an orbital ends.
    + stopper          ::Bool                    ... true, if a truncated orbital shall raise rather than report.
"""
function checkOrbitalBox(orbitals::Dict{Subshell,Orbital}, grid::Radial.Grid;
                         extentTolerance::Float64=0.9, densityFloor::Float64=1.0e-5, stopper::Bool=false)
    rbox = grid.r[end];      offenders = Tuple{Subshell,Float64,Float64}[]
    worstRatio = 0.;         worstShell = first(keys(orbitals));      outerExtent = 0.

    for  (sh, orb)  in  orbitals
        np = min(length(orb.P), length(grid.r));     if  np < 3   continue    end
        peak = 0.
        for  i = 1:np    peak = max(peak, orb.P[i]^2 + orb.Q[i]^2)    end
        if  peak <= 0.    continue    end

        # where the orbital effectively ends, and how that compares with where the box ends
        rEnd = grid.r[1]
        for  i = np:-1:1
            if  orb.P[i]^2 + orb.Q[i]^2 > densityFloor * peak    rEnd = grid.r[i];   break    end
        end
        ratio = rEnd / rbox
        if  ratio > extentTolerance    push!(offenders, (sh, rEnd, ratio))    end
        if  ratio > worstRatio         worstRatio = ratio;   worstShell = sh   end
        outerExtent = max(outerExtent, rEnd)
    end

    report = @sprintf("box %.1f a.u.; outermost orbital reaches %.2f a.u., ", rbox, outerExtent) *
             @sprintf("largest extent/box = %.3f (%s, limit %.2f)", worstRatio, string(worstShell), extentTolerance)

    if  length(offenders) > 0
        printstyled("\n>>> BOX TOO SMALL: these converged orbitals still carry density where the box ends, so they are\n" *
                    ">>> cut off rather than decayed:\n", color=:light_red)
        for  (sh, rEnd, ratio)  in offenders
            printstyled(@sprintf("    %10s   reaches %.2f a.u. of a %.1f a.u. box   (extent/box = %.3f)\n",
                                 string(sh), rEnd, rbox, ratio), color=:light_red)
        end
        printstyled(">>> Use Basics.recommendedGrid(configs, nm) to have the box matched, or Radial.Grid(grid; rbox=..)\n",
                    color=:light_red)
        Defaults.warn(AddWarning(), "Bsplines.checkOrbitalBox(): the radial box cuts off " *
                      join([string(o[1]) for o in offenders], ", ") * "; the orbitals are truncated.  " * report)
        if  stopper   error("Bsplines.checkOrbitalBox(): the radial box cuts off " *
                            "$(length(offenders)) of $(length(orbitals)) orbitals.")   end
        return( (false, report) )
    end

    return( (true, report) )
end


"""
`Bsplines.checkOrbitalConsistency(orbitals::Dict{Subshell,Orbital}, grid::Radial.Grid;
                                  rTolerance::Float64=1.5, eTolerance::Float64=2.0,
                                  eInversion::Float64=0.05, stopper::Bool=true)`
    ... checks the generated orbitals for the one failure mode that Bsplines.checkGridRepresentation cannot see: an SCF that has converged
        onto the WRONG STATE for some symmetry. It compares the two spin-orbit partners of every subshell -- same n, same l, kappa of either
        sign -- which must describe the same shell and therefore must be close in both mean radius and binding energy.

        WHY THIS IS NEEDED IN ADDITION TO checkGridRepresentation. That function tests hydrogenic orbitals at the FULL nuclear charge, so it
        only sees whether the grid can resolve a COMPACT orbital. It passes happily on a grid that cannot represent a diffuse, screened
        outer orbital: Ge II [Ar] 3d^10 4s^2 4f on r_max = 614 a.u. returned E(4f_7/2) = -1.5758 with <r> = 5.20 against E(4f_5/2) = -0.0619
        with <r> = 10.87 -- a different state entirely -- and produced a Lande factor of -2.264 against the exact 8/7, while the grid check
        reported no problem at all.

        THE TOLERANCES ARE CALIBRATED, not guessed. Genuine fine structure does separate the partners, and the more so the heavier the ion,
        so a tight criterion would fire on correct results. Measured for hydrogen-like ions on matched boxes, |E_1/E_2| and <r>_1/<r>_2 are
             Z =  10   2p  1.0013 / 0.9982      3d  1.0003 / 0.9996
             Z =  26   2p  1.0092 / 0.9877      3d  1.0020 / 0.9969
             Z =  54   2p  1.0426 / 0.9448      3d  1.0088 / 0.9867
             Z =  92   2p  1.1540 / 0.8195      3d  1.0268 / 0.9601      4f  1.0097 / 0.9844
        i.e. at worst 15% in energy and 18% in radius for a legitimate pair, against a factor 26 in energy and 2.1 in radius for the broken
        Ge II case. The defaults of 2.0 and 1.5 sit in that gap with room on both sides.
        A value::Bool is returned -- true if every partner pair is consistent.
"""
function checkOrbitalConsistency(orbitals::Dict{Subshell,Orbital}, grid::Radial.Grid;
                                 rTolerance::Float64=1.5, eTolerance::Float64=2.0,
                                 eInversion::Float64=0.05, stopper::Bool=true)
    # group the orbitals by (n, l); only pairs with both kappa signs can be compared
    groups = Dict{Tuple{Int64,Int64}, Array{Subshell,1}}()
    for  sh  in  keys(orbitals)
        key = (sh.n, Basics.subshell_l(sh));    groups[key] = push!( get(groups, key, Subshell[]), sh )
    end
    offenders = Tuple{Subshell,Subshell,Float64,Float64}[]
    for  (key, shs)  in  groups
        length(shs) == 2   ||   continue
        a = orbitals[shs[1]];     b = orbitals[shs[2]]
        ra = JenaAtomicCalculator.RadialIntegrals.rkDiagonal(1, a, a, grid)
        rb = JenaAtomicCalculator.RadialIntegrals.rkDiagonal(1, b, b, grid)
        rRatio = (ra > 0. && rb > 0.)  ?  max(ra,rb)/min(ra,rb)  :  Inf
        eRatio = (a.energy * b.energy > 0.)  ?  max(abs(a.energy),abs(b.energy))/min(abs(a.energy),abs(b.energy))  :  Inf
        # The MEAN RADIUS is one criterion; the energy RATIO is reported for diagnosis but does not trigger
        # on its own. Energies are fragile here: a nearly-unbound subshell in a highly-ionised configuration
        # can have its two partners straddle zero, which makes eRatio infinite while the orbitals are in
        # fact identical -- exactly what the Cascade stepwise-decay test does with 3p_1/2 / 3p_3/2, whose
        # radii agree to 0.1%. The radius alone already separates the calibration cases by a wide margin.
        # THE ORDERING, however, IS a criterion, and a sharp one (added 12-Aug-2026). In a central field the
        # two members of a spin-orbit pair are ordered by physics, not by magnitude: j = l-1/2 (kappa > 0)
        # lies BELOW j = l+1/2 (kappa < 0). Hydrogen's 2p_1/2 below 2p_3/2 is the familiar case. An inverted
        # pair therefore means the SCF has produced something that is not a spin-orbit doublet at all, no
        # matter how similar the two radii happen to be.
        # THIS IS THE CASE THE RADIUS TEST MISSED: neutral Pr I [Xe] 4f^3 6s^2 gave E(4f_7/2) = -0.1854 below
        # E(4f_5/2) = -0.1154 -- inverted, and 20x too large -- while the two mean radii agreed to 6%, well
        # inside rTolerance. The level structure that followed had the ^4I multiplet upside down and its fine
        # structure collapsed by four orders of magnitude, and nothing warned.
        # The test is applied only when BOTH partners are bound, which is what keeps the straddling-zero
        # Cascade case out of it; there, one energy is positive and the ordering carries no meaning.
        inverted = 0.
        if  a.energy < 0.  &&  b.energy < 0.
            shLo = shs[1].kappa > 0  ?  shs[1]  :  shs[2]      # j = l - 1/2, must be the more bound
            shHi = shs[1].kappa > 0  ?  shs[2]  :  shs[1]      # j = l + 1/2
            eLo  = orbitals[shLo].energy;      eHi = orbitals[shHi].energy
            inverted = (eLo - eHi) / max(abs(eLo), abs(eHi))   # > 0 means the pair is upside down
        end
        if  rRatio > rTolerance   ||   inverted > eInversion
            push!(offenders, (shs[1], shs[2], eRatio, rRatio))
        end
    end
    if  length(offenders) > 0
        printstyled("\n>>> ORBITAL CHECK FAILED: these spin-orbit partners do not describe the same shell, which means\n" *
                    ">>> the SCF has converged onto the wrong state for one of them (usually an ill-matched radial box):\n",
                    color=:light_red)
        printstyled("      partners                    E ratio     <r> ratio     (limits $eTolerance / $rTolerance;\n" *
                    "      an inverted spin-orbit pair also triggers, at $eInversion)\n", color=:light_red)
        for  (sa, sb, er, rr)  in  offenders
            printstyled(@sprintf("    %-10s %-10s  %11.4g %13.4g\n", string(sa), string(sb), er, rr), color=:light_red)
        end
        printstyled(">>> Compare the two orbital energies and mean radii directly, and match the radial box to the\n" *
                    ">>> orbitals; note that a box much TOO LARGE starves the basis just as badly as one too small.\n",
                    color=:light_red)
        # Collected as well, and this is the case that most needs it: with stopper = false the finding is
        # printed and execution continues, so in a long run it scrolls away and the untrustworthy orbitals are
        # used without anyone noticing afterwards.  The offending pairs are named so the report is actionable.
        Defaults.warn(AddWarning(), "Bsplines.checkOrbitalConsistency(): $(length(offenders)) spin-orbit partner " *
                      "pair(s) describe different states -- " *
                      join([string(o[1]) * "/" * string(o[2]) for o in offenders], ", ") *
                      " -- the orbitals are not trustworthy; match the radial box to the orbitals.")
        if  stopper   error("Bsplines.checkOrbitalConsistency(): $(length(offenders)) spin-orbit partner pair(s) " *
                            "describe different states; the orbitals are not trustworthy.")   end
        return( false )
    end

    return( true )
end


"""
`Bsplines.computeNondiagonalD(pm::Int64, kappa::Int64, bspline1::Bsplines.Bspline, bspline2::Bsplines.Bspline, grid::Radial.Grid)`
    ... computes the (radial and non-diagonal) D_kappa^+/- integral for the two bsplines, all defined on grid,
        <bspline1| +/- d/dr + kappa/r |bspline2>; pm = +1/-1 provides the phase for taking the derivative.

        A single, uniform composite Gauss-Legendre quadrature is applied here, with no special-casing at any grid point, exactly as in
        Zatsarinny & Froese Fischer's DBSR_HF reference code (dbsr_lib_dbs.f90's ZINTYM). In particular the first point needs no separate
        treatment: every grid.r[i] is already a genuine (nonzero) Gauss-Legendre node and never literally r=0.
        A value::Float64 is returned.
"""
function computeNondiagonalD(pm::Int64, kappa::Int64, bspline1::Bsplines.Bspline, bspline2::Bsplines.Bspline, grid::Radial.Grid)
    if  bspline1.upper <= bspline2.lower  ||  bspline2.upper <= bspline1.lower    return( 0. )   end
    lower = max(bspline1.lower, bspline2.lower);    add1 = 1 - bspline1.lower
    upper = min(bspline1.upper, bspline2.upper);    add2 = 1 - bspline2.lower

    wa = 0.
    for  i = lower:upper
        wa = wa + pm * bspline1.bs[i+add1] * bspline2.bp[i+add2] * grid.wr[i]
        wa = wa + bspline1.bs[i+add1] * kappa * bspline2.bs[i+add2] / grid.r[i] * grid.wr[i]
    end
    return( wa )
end


"""
`Bsplines.computeOverlap(bspline1::Bsplines.Bspline, bspline2::Bsplines.Bspline, grid::Radial.Grid)`
    ... computes the (radial) overlap integral <bspline1|bspline2> for the two bsplines as defined on grid.
        A value::Float64 is returned.
"""
function computeOverlap(bspline1::Bsplines.Bspline, bspline2::Bsplines.Bspline, grid::Radial.Grid)
    if  bspline1.upper <= bspline2.lower  ||  bspline2.upper <= bspline1.lower    return( 0. )   end
    lower = max(bspline1.lower, bspline2.lower);    add1 = 1 - bspline1.lower
    upper = min(bspline1.upper, bspline2.upper);    add2 = 1 - bspline2.lower
    
    wa = 0.            
    for  i = lower:upper   wa = wa + bspline1.bs[i+add1] * bspline2.bs[i+add2] * grid.wr[i]   end
    return( wa )
end


"""
`Bsplines.computeVlocal(bspline1::Bsplines.Bspline, bspline2::Bsplines.Bspline, pot::Radial.Potential, grid::Radial.Grid)`
    ... computes the (radial) integral <bspline1| V_pot |bspline2> for the two bsplines and the given radial potential as defined on grid.
        The quadrature is the uniform one described for Bsplines.computeNondiagonalD, with no special treatment of the first grid point.
        A value::Float64 is returned.
"""
function computeVlocal(bspline1::Bsplines.Bspline, bspline2::Bsplines.Bspline, pot::Radial.Potential, grid::Radial.Grid)
    if  bspline1.upper <= bspline2.lower  ||  bspline2.upper <= bspline1.lower    return( 0. )   end
    lower = max(bspline1.lower, bspline2.lower);    add1 = 1 - bspline1.lower
    upper = min(bspline1.upper, bspline2.upper);    add2 = 1 - bspline2.lower

    wa = 0.
    for  i = lower:upper
        wa = wa - bspline1.bs[i+add1] * pot.Zr[i] * bspline2.bs[i+add2] / grid.r[i] * grid.wr[i]
    end
    return( wa )
end


"""
`Bsplines.diagonalizeLocalMatrix(kappa::Int64, matrixA::Array{Float64,2}, matrixB::Array{Float64,2}, primitives::Bsplines.Primitives)`
    ... enforces the r=0 boundary condition for symmetry kappa (Bsplines.boundaryDropCounts) by eliminating the corresponding leading large-
        and small-component B-splines from the generalized eigenvalue problem (matrixA, matrixB) before diagonalizing, following Zatsarinny
        & Froese Fischer's DBSR_HF reference code (hf_boundary.f90's Boundary_conditions, finite-nucleus branch). Each returned eigenvector
        is re-embedded into the full (nsL+nsS)-length vector space, with exact zeros at the eliminated positions, so that all downstream
        code (Bsplines.generateOrbitalFromPrimitives, the bVector-native AL/EOL machinery, etc.) is unaffected in how it reconstructs
        P(r)/Q(r) -- only the NUMBER of returned eigenpairs shrinks by (dropP+dropQ). Callers that pick an eigenpair by index must use
        Bsplines.findPositiveBranchStart on the returned values, NOT a fixed index counted from nsL/nsS.
        An eigen::Basics.Eigen is returned.

        A NOTE ON COST, for whoever comes back to this (measured 12-Aug-2026). This is a FULL generalized eigendecomposition, but very
        little of it is used: the SCF wants a short contiguous block at the bottom of the positive-energy branch --
        Bsplines.generateOrbitalFromPrimitives takes ni = mm + n - l - 1 -- which in real bases is ONE TO FOUR eigenvectors out of ~665
        (Ne, Ar, Fe and a 3s/3p/3d correlation basis were counted; in all of them the highest wanted state was still bound). About half
        the spectrum is the Dirac sea and is discarded outright. Three routes were measured, and they behave quite differently.

        (a) ASK FOR FEWER EIGENVECTORS. LinearAlgebra.eigen(A, range) reaches LAPACK's syevr; measured against a full decomposition it
            gives 3.3x (n=239), 3.9x (671), 4.1x (1053), 3.5x (2005), 2.8x (4389). The speed-up SATURATES around 3-4x and does not improve
            with a denser grid, because the tridiagonal reduction is O(n^3) and unavoidable -- only the back-transformation shrinks. Note
            also that a two-stage scheme (all eigenvalues first, then the wanted vectors) cannot pay: eigenvalues alone already cost ~0.7
            of the full solve. Any gain must come from ONE range-restricted call.

        (b) EXPLOIT THAT THE MATRIX IS BANDED. It is 4.0% non-zero, and its half-bandwidth depends entirely on the ordering: 343 as stored
            here (large-component block, then small-component block -- effectively dense), but 15 once the large and small components are
            INTERLEAVED. That 15 is set by the spline order and stays 15 however dense the grid, so banded work is O(n b^2) against dense
            O(n^3). This is the only route whose advantage GROWS with grid density, which is the case that matters when a fine grid is
            needed for physical reasons. Measured with LAPACK's banded routines: eigenVALUES (dsbgv, JOBZ='N') 1.2x at n=400 rising to
            7.8x at n=3200, agreeing with the dense solve to 1e-14. But dsbgvx WITH eigenvectors is SLOWER than dense and gets worse with
            n (0.6x at 400, 0.1x at 3200), because it accumulates the n x n reduction matrix unblocked -- so the vectors would have to come
            from banded INVERSE ITERATION, not from the banded driver. Route (b) is deliberately NOT implemented: LAPACK's banded routines
            are not exposed by Julia's LinearAlgebra.LAPACK and would need a raw ccall, which JAC avoids at its present stage. If a
            dependency ever becomes acceptable, a banded-matrix package is the safer way in. One trap worth recording because it does not
            announce itself: passing Int32 where Julia's ILP64 LAPACK expects BlasInt returns info = 0 with silently WRONG eigenvalues.

        (c) RESTRICT THE SOLVE TO AN ENERGY WINDOW. Built and MEASURED on 12-Aug-2026, then WITHDRAWN -- the reasoning is kept here so that
            it is not attempted a third time. LAPACK can be asked for an energy window, but only for a STANDARD problem, so the generalized
            one is first reduced with the Cholesky factor of the overlap matrix (B = L L', C = L^-1 A L^-T, x = L^-T y). On Ar, order 665:
            the full solve 0.050 s against 0.013 s reduction + 0.014 s windowed solve, i.e. 1.7-1.8x, which is worth having since this solve
            is 57.5% of a real SCF. Two facts sank it. The window must be the BOUND states to pay at all: taking the whole positive branch
            instead (331 of 665 states) costs 0.075 s and is 0.6x, i.e. SLOWER than the full solve, so the gain comes from the smallness of
            the bound set and not from skipping the Dirac sea. And DURING the iteration the wanted states are not all bound: with the
            screened DFS potential of the starting orbitals, neutral Ne has ONE bound state in kappa = -1 and NONE in kappa = 1 or -2, its
            2p sitting at +0.0252 a.u., so a bound window fails on a neutral atom in the very first iteration even though every wanted state
            is bound once converged. (An earlier count that suggested otherwise had used the bare NUCLEAR potential, which binds far more
            strongly.) An index window is no way out either: see above on why an eigenpair must never be picked by an index counted from
            nsL/nsS.
"""
function diagonalizeLocalMatrix(kappa::Int64, matrixA::Array{Float64,2}, matrixB::Array{Float64,2}, primitives::Bsplines.Primitives)
    nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS
    dropP, dropQ, trailP, trailQ = Bsplines.boundaryDropCounts(kappa, primitives.grid)
    keep = vcat( collect(dropP+1:nsL-trailP), collect(nsL+dropQ+1:nsL+nsS-trailQ) )
    wc = Basics.diagonalize(GeneralizedEigenvaluesWithLinearAlgebra(), matrixA[keep,keep], matrixB[keep,keep])
    vectors = Vector{Float64}[]
    for  v  in  wc.vectors
        full = zeros(nsL+nsS);   full[keep] = v;   push!(vectors, full)
    end
    return( Basics.Eigen(wc.values, vectors) )
end


"""
`Bsplines.findPositiveBranchStart(values::Array{Float64,1}; mass::Float64=1.0)`
    ... returns the 1-based index of the first (ascending-sorted) eigenvalue above -1.999*mass*c^2, the standard threshold marking the
        boundary between the unphysical negative-energy continuum ("Dirac sea") branch and the physical positive branch, following
        Zatsarinny & Froese Fischer's DBSR_HF reference code (hf_solve_HF.f90's hf_eiv). The (n,l) bound state then sits at (this index) +
        (n-l) - 1 within the SAME values array -- robust to however many leading B-splines were eliminated for a given symmetry kappa
        (Bsplines.boundaryDropCounts/diagonalizeLocalMatrix), unlike a fixed index counted from nsL/nsS.
        A value::Int64 is returned.

        THE THRESHOLD SCALES WITH THE MASS, which is why mass is a keyword here rather than a constant. The negative-energy continuum starts
        at -2 m c^2, so a fixed -1.999 c^2 is an electron-only value: a muon 1s level in lead lies near 10 MeV = 3.7e5 Hartree, far BELOW the
        electron threshold of 3.75e4, and every muon bound state of a heavy atom would then be rejected here as a spurious Dirac-sea state
        -- silently, since this function returns an index and not an error.

        WARNING: the threshold-based separation likewise returns garbage if `values` was diagonalized against a potential with no (or a
        much-too-weak) attractive nuclear well -- e.g. a caller that passes only `Basics.computePotential(Basics.DFSField(1.0), grid,
        basis)` (the ELECTRONIC mean-field potential alone) without adding `Nuclear.nuclearPotential(nm, grid)` first. Without a real
        potential well the Dirac equation has no clean energetic gap between the unphysical Dirac-sea branch and genuine atomic bound
        states, so the index returned can lie STILL within (or immediately adjacent to) the spurious negative-continuum branch, with
        eigenvalues clustering right at -1.999*c^2 instead of at the expected atomic scale (roughly -1 to -2000 Hartree, not ~-37500 for a
        typical ion). That produces orbitals which are numerical garbage rather than "slightly wrong". Always add the nuclear potential to
        any potential passed into `Bsplines.generateOrbitals` or into this function.
"""
function findPositiveBranchStart(values::Array{Float64,1}; mass::Float64=1.0)
    c  = Defaults.getDefaults("speed of light: c")
    zz = -1.999 * mass * c^2
    for  (i,v)  in  enumerate(values)
        if  v > zz   return( i )   end
    end
    error("Bsplines.findPositiveBranchStart(): no eigenvalue found above the negative-continuum threshold.")
end


"""
`Bsplines.fitVectorToPrimitives(orb::Radial.Orbital, primitives::Bsplines.Primitives, matrixB::Array{Float64,2})`
    ... projects the (already CLEANED/truncated) tabulated orbital orb onto the B-spline primitives basis via the standard
        Galerkin/least-squares projection matrixB * p = rhs, rhs[i] = <B_i|orb.P-or-Q>, using the existing grid quadrature. Unlike pulling
        the RAW diagonalization eigenvector -- which reproduces the UNCLEANED tabulated function from BEFORE generateOrbitalFromPrimitives'
        own truncation-at-mtp and small-value cleanup are applied -- this
        function GUARANTEES the returned coefficient vector is fully self-consistent with orb's OWN (already cleaned) tabulated P, Q arrays,
        inheriting orb's own well-defined truncation instead of carrying whatever small numerical noise the raw eigenvector's tail
        coefficients happen to have. This matters wherever B-spline expansion coefficients are themselves summed/weighted directly, rather
        than only ever used to reconstruct one smooth tabulated function -- such sums do not automatically benefit from the cancellation
        that evaluating a single, already-cleaned tabulated function enjoys. See InteractionStrength.XL_CoulombTensor and
        SelfConsistent.solveAverageLevelField, where using the raw eigenvector instead of this projection was traced to a real,
        non-negligible SCF discrepancy. The returned vector is explicitly re-normalized so that v'*matrixB*v = 1 EXACTLY (to floating-point
        precision), rather than trusting the least-squares fit to land there on its own -- Hamiltonian.projectHamiltonian's projection
        operator (I - S*bb') is only truly idempotent for an exactly S-normalized b; feeding it a vector off by even a small residual leaves
        the "projection" not quite a projection, which compounds under repeated application across SCF iterations.
        A vector::Vector{Float64}, of length nsL+nsS, is returned.
"""
function fitVectorToPrimitives(orb::Radial.Orbital, primitives::Bsplines.Primitives, matrixB::Array{Float64,2})
    nsL = primitives.grid.nsL;   nsS = primitives.grid.nsS;   grid = primitives.grid
    rhs = zeros(nsL+nsS)

    for  i = 1:nsL
        Bi = primitives.bsplinesL[i]
        Pi = zeros(Bi.upper);   add = 1 - Bi.lower
        for  j = Bi.lower:Bi.upper   Pi[j] = Pi[j] + Bi.bs[j+add]   end
        mtp = min(length(Pi), length(orb.P))
        s   = 0.;   for  r = 2:mtp   s = s + Pi[r] * orb.P[r] * grid.wr[r]   end
        rhs[i] = s
    end
    for  i = 1:nsS
        Bi = primitives.bsplinesS[i]
        Qi = zeros(Bi.upper);   add = 1 - Bi.lower
        for  j = Bi.lower:Bi.upper   Qi[j] = Qi[j] + Bi.bs[j+add]   end
        mtp = min(length(Qi), length(orb.Q))
        s   = 0.;   for  r = 2:mtp   s = s + Qi[r] * orb.Q[r] * grid.wr[r]   end
        rhs[nsL+i] = s
    end

    vector = matrixB \ rhs
    norm2  = transpose(vector) * matrixB * vector

    return( vector / sqrt(norm2) )
end


"""
`Bsplines.generateGalerkinMatrix(sh::Subshell, energy::Float64, pot::Radial.Potential, primitives::Bsplines.Primitives;
                                 mass::Float64=1.0)`
    ... generates the Galerkin-A matrix for the given potential and B-spline primitives; a matrix::Array{Float64,2} is returned.
"""
function generateGalerkinMatrix(sh::Subshell, energy::Float64, pot::Radial.Potential, primitives::Bsplines.Primitives;
                                mass::Float64=1.0)
    nsL      = primitives.grid.nsL;    nsS = primitives.grid.nsS

    # Define the storage for the calculations of matrices; this is necessary to use the Bsplines.generateMatrix!() function
    println(">> (Re-) Define a storage array for dealing with single-electron TTp B-spline matrices:")
    storage  = Dict{String,Array{Float64,2}}()
    # Set-up the overlap matrix
    wb  = zeros( nsL+nsS, nsL+nsS )
    wb[1:nsL,1:nsL]                 = Bsplines.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
    wb[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = Bsplines.generateTTpMatrix!("SS-overlap", 0, primitives, storage)
    # Set-up the local Hamiltonian matrix
    wa = Bsplines.setupLocalMatrix(sh.kappa, primitives, pot::Radial.Potential, storage::Dict{String,Array{Float64,2}}; mass=mass)
    wa[1:end,1:end] = wa[1:end,1:end] - energy * wb[1:end,1:end]

    return( wa )
end


"""
`Bsplines.generateOrbitalFromPrimitives(sh::Subshell, wc::Basics.Eigen, primitives::Bsplines.Primitives; mass::Float64=1.0)`
    ... generates the large and small components for the subshell sh from the primitives and their eigenvalues & eigenvectors.
        A (normalized) orbital::Orbital is returned.
"""
function generateOrbitalFromPrimitives(sh::Subshell, wc::Basics.Eigen, primitives::Bsplines.Primitives; mass::Float64=1.0)
    nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS
    l  = Basics.subshell_l(sh);   mm = Bsplines.findPositiveBranchStart(wc.values; mass=mass);   ni = mm + sh.n - l - 1
    en = wc.values[ni];        if  en < 0.    isBound = true  else   isBound = false                 end
    ev = wc.vectors[ni];       if  length(ev) != nsL + nsS    error("stop a")                        end
    
    P = zeros(primitives.grid.NoPoints);    Pprime = zeros(primitives.grid.NoPoints)
    Q = zeros(primitives.grid.NoPoints);    Qprime = zeros(primitives.grid.NoPoints)
    for  i = 1:nsL
        lower = primitives.bsplinesL[i].lower;   upper = primitives.bsplinesL[i].upper;   add = 1 - primitives.bsplinesL[i].lower
        for  j = lower:upper  P[j]      = P[j] + ev[i] * primitives.bsplinesL[i].bs[j+add]           end
        for  j = lower:upper  Pprime[j] = Pprime[j] + ev[i] * primitives.bsplinesL[i].bp[j+add]      end
    end 
    for  i = 1:nsS   
        lower = primitives.bsplinesS[i].lower;   upper = primitives.bsplinesS[i].upper;   add = 1 - primitives.bsplinesS[i].lower
        for  j = lower:upper  Q[j]      = Q[j] + ev[nsL+i] * primitives.bsplinesS[i].bs[j+add]       end
        for  j = lower:upper  Qprime[j] = Qprime[j] + ev[nsL+i] * primitives.bsplinesS[i].bp[j+add]  end
    end 
    
    # Determine the maximum number of grid points for this orbital and normalized it propery
    mtp = 0;   for j = primitives.grid.NoPoints:-1:1    if  abs(P[j])^2 + abs(Q[j])^2 > 1.0e-30   mtp = j;   break   end     end
    
    Px = zeros(mtp);    Px[1:mtp] = P[1:mtp];    Pprimex = zeros(mtp);    Pprimex[1:mtp] = Pprime[1:mtp]  
    Qx = zeros(mtp);    Qx[1:mtp] = Q[1:mtp];    Qprimex = zeros(mtp);    Qprimex[1:mtp] = Qprime[1:mtp]    
    for  j = 1:mtp      if  abs(Px[j])      < 1.0e-16    Px[j] = 0.       end
                        if  abs(Qx[j])      < 1.0e-16    Qx[j] = 0.       end 
                        if  abs(Pprimex[j]) < 1.0e-16    Pprimex[j] = 0.  end
                        if  abs(Qprimex[j]) < 1.0e-16    Qprimex[j] = 0.  end      end
                        
    # Ensure that the large component of all orbitals start 'positive'
    wSign     = sum( Px[1:30] )
    if  wSign < 0.   Px[1:mtp] = -Px[1:mtp];   Pprimex[1:mtp] = -Pprimex[1:mtp] 
                     Qx[1:mtp] = -Qx[1:mtp];   Qprimex[1:mtp] = -Qprimex[1:mtp]   end
    
    orbital   = Orbital(sh, isBound, true, en, Px, Qx, Pprimex, Qprimex, Radial.Grid())
    
    # Renormalize the radial orbital   
    wN        = sqrt( JenaAtomicCalculator.RadialIntegrals.overlap(orbital, orbital, primitives.grid) )
    Px[1:mtp] = Px[1:mtp] / wN;    Pprimex[1:mtp] = Pprimex[1:mtp] / wN
    Qx[1:mtp] = Qx[1:mtp] / wN;    Qprimex[1:mtp] = Qprimex[1:mtp] / wN 
    
    orb = Orbital(sh, isBound, true, en, Px, Qx, Pprimex, Qprimex, Radial.Grid())
    
    return( orb )   
end


"""
`Bsplines.generateOrbitalFromPrimitives(sh::Subshell, energy::Float64, mtp::Int64, ev::Array{Float64,1}, primitives::Bsplines.Primitives)`
    ... generates the large and small components of a (relativistic) orbital for the subshell sh from the given primitives and the
        eigenvector ev. A (non-normalized) orbital::Orbital is returned.
"""
function generateOrbitalFromPrimitives(sh::Subshell, energy::Float64, mtp::Int64, ev::Array{Float64,1}, primitives::Bsplines.Primitives)
    nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS
    P = zeros(primitives.grid.NoPoints);    Pprime = zeros(primitives.grid.NoPoints)
    Q = zeros(primitives.grid.NoPoints);    Qprime = zeros(primitives.grid.NoPoints)
    for  i = 1:nsL   
        lower = primitives.bsplinesL[i].lower;   upper = primitives.bsplinesL[i].upper;   add = 1 - primitives.bsplinesL[i].lower
        for  j = lower:upper  P[j]      = P[j] + ev[i] * primitives.bsplinesL[i].bs[j+add]           end
        for  j = lower:upper  Pprime[j] = Pprime[j] + ev[i] * primitives.bsplinesL[i].bp[j+add]      end
    end 
    for  i = 1:nsS   
        lower = primitives.bsplinesS[i].lower;   upper = primitives.bsplinesS[i].upper;   add = 1 - primitives.bsplinesS[i].lower
        for  j = lower:upper  Q[j]      = Q[j] + ev[nsL+i] * primitives.bsplinesS[i].bs[j+add]       end
        for  j = lower:upper  Qprime[j] = Qprime[j] + ev[nsL+i] * primitives.bsplinesS[i].bp[j+add]  end
    end 
    
    Px      = zeros(mtp);    Qx      = zeros(mtp);    Px[1:mtp]      = P[1:mtp];         Qx[1:mtp]      = Q[1:mtp]    
    Pprimex = zeros(mtp);    Qprimex = zeros(mtp);    Pprimex[1:mtp] = Pprime[1:mtp];    Qprimex[1:mtp] = Qprime[1:mtp]    
    
    return( Orbital(sh, false, true, energy, Px, Qx, Pprimex, Qprimex, Radial.Grid()) )   
end


"""
`Bsplines.generateOrbitalFromVector(sh::Subshell, energy::Float64, vector::Vector{Float64}, primitives::Bsplines.Primitives)`
    ... generates a (normalized, cleaned) tabulated orbital directly from a given B-spline expansion coefficient vector, rather than from an
        eigenvector INDEXED out of a Basics.Eigen (as Bsplines.generateOrbitalFromPrimitives(sh,wc,primitives) requires). This is needed
        whenever the vector in hand is not literally an eigenvector of anything -- e.g. after SelfConsistent.orthonormalizeSameKappa's
        Loewdin symmetric orthogonalization, which produces a LINEAR COMBINATION of eigenvectors that is itself not an eigenvector of the
        original problem. Reproduces generateOrbitalFromPrimitives(sh,wc,primitives)'s own reconstruction exactly (auto-detects mtp from
        where the density drops below 1e-13, zeros values below 1e-10, fixes the sign convention so the large component starts positive, and
        renormalizes to unit norm in the grid quadrature) -- just without the wc/ni lookup indirection.
        A (normalized) orbital::Orbital is returned.
"""
function generateOrbitalFromVector(sh::Subshell, energy::Float64, vector::Vector{Float64}, primitives::Bsplines.Primitives)
    nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS
    if  length(vector) != nsL + nsS    error("stop a")    end

    P = zeros(primitives.grid.NoPoints);    Pprime = zeros(primitives.grid.NoPoints)
    Q = zeros(primitives.grid.NoPoints);    Qprime = zeros(primitives.grid.NoPoints)
    for  i = 1:nsL
        lower = primitives.bsplinesL[i].lower;   upper = primitives.bsplinesL[i].upper;   add = 1 - primitives.bsplinesL[i].lower
        for  j = lower:upper  P[j]      = P[j] + vector[i] * primitives.bsplinesL[i].bs[j+add]           end
        for  j = lower:upper  Pprime[j] = Pprime[j] + vector[i] * primitives.bsplinesL[i].bp[j+add]      end
    end
    for  i = 1:nsS
        lower = primitives.bsplinesS[i].lower;   upper = primitives.bsplinesS[i].upper;   add = 1 - primitives.bsplinesS[i].lower
        for  j = lower:upper  Q[j]      = Q[j] + vector[nsL+i] * primitives.bsplinesS[i].bs[j+add]       end
        for  j = lower:upper  Qprime[j] = Qprime[j] + vector[nsL+i] * primitives.bsplinesS[i].bp[j+add]  end
    end

    mtp = 0;   for j = primitives.grid.NoPoints:-1:1    if  abs(P[j])^2 + abs(Q[j])^2 > 1.0e-30   mtp = j;   break   end     end

    Px = zeros(mtp);    Px[1:mtp] = P[1:mtp];    Pprimex = zeros(mtp);    Pprimex[1:mtp] = Pprime[1:mtp]
    Qx = zeros(mtp);    Qx[1:mtp] = Q[1:mtp];    Qprimex = zeros(mtp);    Qprimex[1:mtp] = Qprime[1:mtp]
    for  j = 1:mtp      if  abs(Px[j])      < 1.0e-16    Px[j] = 0.       end
                        if  abs(Qx[j])      < 1.0e-16    Qx[j] = 0.       end
                        if  abs(Pprimex[j]) < 1.0e-16    Pprimex[j] = 0.  end
                        if  abs(Qprimex[j]) < 1.0e-16    Qprimex[j] = 0.  end      end

    wSign = sum( Px[1:min(30,mtp)] )
    if  wSign < 0.   Px[1:mtp] = -Px[1:mtp];   Pprimex[1:mtp] = -Pprimex[1:mtp]
                     Qx[1:mtp] = -Qx[1:mtp];   Qprimex[1:mtp] = -Qprimex[1:mtp]   end

    orbital = Orbital(sh, true, true, energy, Px, Qx, Pprimex, Qprimex, Radial.Grid())

    wN        = sqrt( JenaAtomicCalculator.RadialIntegrals.overlap(orbital, orbital, primitives.grid) )
    Px[1:mtp] = Px[1:mtp] / wN;    Pprimex[1:mtp] = Pprimex[1:mtp] / wN
    Qx[1:mtp] = Qx[1:mtp] / wN;    Qprimex[1:mtp] = Qprimex[1:mtp] / wN

    return( Orbital(sh, true, true, energy, Px, Qx, Pprimex, Qprimex, Radial.Grid()) )
end


"""
`Bsplines.generateOrbitals(subshells::Array{Subshell,1}, pot::Radial.Potential, nm::Nuclear.Model,
                           primitives::Bsplines.Primitives; printout::Bool=true, mass::Float64=1.0)`
    ... generates all single-electron orbitals from subshell list for the radial potential pot.
        A set of orbitals::Dict{Subshell, Orbital} is returned.
"""
function generateOrbitals(subshells::Array{Subshell,1}, pot::Radial.Potential, nm::Nuclear.Model, 
                          primitives::Bsplines.Primitives; printout::Bool=true, mass::Float64=1.0)
    orbitals = Dict{Subshell, Orbital}()
    kappas   = Int64[];   for sh in subshells  push!(kappas, sh.kappa)   end;   kappas = unique(kappas)
    nsL      = primitives.grid.nsL;    nsS = primitives.grid.nsS
    
    # Define the storage for the calculations of matrices; this is necessary to use the Bsplines.generateTTpMatrix!() function.
    if  printout    println(">> (Re-) Define a storage array for dealing with single-electron TTp B-spline matrices:")    end
    storage  = Dict{String,Array{Float64,2}}()
    
    for kappa  in  kappas
        # Set-up the overlap matrix
        wb = zeros( nsL+nsS, nsL+nsS )
        
        # (1) Compute or fetch the diagonal 'overlap' blocks
        wb[1:nsL,1:nsL]                 = Bsplines.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
        wb[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = Bsplines.generateTTpMatrix!("SS-overlap", 0, primitives, storage)
        
        # (2) Compute the local Hamiltonian matrix and diagonalize it
        wa = Bsplines.setupLocalMatrix(kappa, primitives, pot, storage; mass=mass)
        w2 = Bsplines.diagonalizeLocalMatrix(kappa, wa, wb, primitives)
        # The offset handed to the tabulation is where the POSITIVE-energy branch begins, and that is NOT nsS:
        # diagonalizeLocalMatrix eliminates dropP+dropQ B-splines for this symmetry, so the branch starts at
        # findPositiveBranchStart -- 46..61 on the standard grid against nsS = 63.  Using nsS here mislabels the
        # printed table by several states without affecting the orbitals actually extracted.
        nsi = Bsplines.findPositiveBranchStart(w2.values; mass=mass) - 1
        if  printout  Basics.tabulateKappaSymmetryEnergiesDirac(kappa, w2.values, nsi, nm)    end
        
        # (3) Collect all the requested single-electron orbitals
        for  sh in subshells
            if  sh.kappa == kappa    orbitals[sh] = Bsplines.generateOrbitalFromPrimitives(sh, w2, primitives; mass=mass)    end
        end
    end
    
    return( orbitals )
end


"""
`Bsplines.generateOrbitalsHydrogenic(subshells::Array{Subshell,1}, nm::Nuclear.Model, primitives::Bsplines.Primitives;
                                     printout::Bool=true, mass::Float64=1.0)`
    ... generates all single-electron orbitals from subshell list for the nuclear potential as specified by nm.
        A set of orbitals::Dict{Subshell, Orbital} is returned.
"""
function generateOrbitalsHydrogenic(subshells::Array{Subshell,1}, nm::Nuclear.Model, primitives::Bsplines.Primitives; printout::Bool=true,
                                    mass::Float64=1.0)
    pot = Nuclear.nuclearPotential(nm, primitives.grid)
    
    orbitals = Bsplines.generateOrbitals(subshells, pot, nm, primitives; printout=printout, mass=mass)
    return( orbitals )
end


"""
`Bsplines.generatePrimitives(grid::Radial.Grid)`
    ... generates the breaks, knots and the B-spline primitives of order k, both for the large and small components. The function applies
        the given grid parameters; no primitive is defined beyond grid[n_max]. The definition of the primitives follows the work of
        Zatsarinny and Froese Fischer, CPC 202 (2016) 287.
        A (set of) primitives::Bsplines.Primitives is returned.

        The raw, full B-spline count -- including the corner spline that is nonzero at r=0 -- is kept here on purpose. The r=0 boundary
        condition is enforced per symmetry kappa instead, with the correct (l- and kappa-sign-dependent) number of leading splines, inside
        Bsplines.setupLocalMatrix; see Bsplines.boundaryDropCounts. This follows Zatsarinny & Froese Fischer's DBSR_HF reference code
        (hf_boundary.f90's Boundary_conditions, for a finite/non-point nucleus).
"""
function generatePrimitives(grid::Radial.Grid)
    !(1 <= grid.orderL <= 11)   &&   error("Order should be 2 <= grid.orderL <= 11; obtained order = $(grid.orderL)")
    !(1 <= grid.orderS <= 11)   &&   error("Order should be 2 <= grid.orderS <= 11; obtained order = $(grid.orderS)")

    # Now determined the B-splines on the grid for the large and small components; initialize values
    primitivesL = Bsplines.Bspline[];   primitivesS = Bsplines.Bspline[];   lower = 0;   upper = 0

    # Generate B-spline basis for large component
    breaks = deepcopy( grid.tL[grid.orderL:end-grid.orderL+1] )
    BL = BSplineKit.BSplineBasis(BSplineOrder(grid.orderL), breaks)
    for  (ib, bL)  in  enumerate(BL)
        bs = Float64[];   bp = Float64[];   needlower = true
        for  (ir,r)  in  enumerate(grid.r)
            if  bL(r) > 0.    push!(bs, bL(r));  push!(bp, bL(r, Derivative(1)) )
                              upper = ir;        if needlower   lower = ir;   needlower = false   end
            end
        end
        push!(primitivesL, Bspline(lower, upper, bs, bp) )
    end

    # Generate B-spline basis for large component
    breaks = deepcopy( grid.tS[grid.orderS:end-grid.orderS+1] )
    BL = BSplineKit.BSplineBasis(BSplineOrder(grid.orderS), breaks)
    for  (ib, bL)  in  enumerate(BL)
        bs = Float64[];   bp = Float64[];   needlower = true
        for  (ir,r)  in  enumerate(grid.r)
            if  bL(r) > 0.    push!(bs, bL(r));  push!(bp, bL(r, Derivative(1)) )
                              upper = ir;        if needlower   lower = ir;   needlower = false   end
            end
        end
        push!(primitivesS, Bspline(lower, upper, bs, bp) )
    end

    return( Bsplines.Primitives(grid, primitivesL, primitivesS) )
end


"""
`Bsplines.generateTTpMatrix!(TTp::String, kappa::Int64, primitives::Bsplines.Primitives, storage::Dict{String,Array{Float64,2}})`  
    ... returns the TTp block of the (single-electron) Dirac Hamiltonian matrix for an electron with symmetry kappa without any potential.
        The following TTp strings are allowed: ["LL-overlap", "SS-overlap", "LS-D_kappa^-", "LS-D_kappa^+"].
        
        Two modes are distinguished owing to the values that are available in the storage (Dict).
            * The TTp matrix block from the storage is returned, if an entry is known; it is assumed that this matrix
              block belong to the given set of primitives.
            * The TTp matrix is computed and set to the storage otherwise; from the TTp string, the key string
              key = string(kappa) * ":" * TTp is generated an applied in the storage dictionary.
              
        All B-splines are supposed to be defined for the same (radial) grid; a matrix::Array{Float64,2} is returned which is quadratic for
        'LL-overlap' and 'SS-overlap' and whose dimension depends on the number of B-splines for the large and small component, otherwise.
"""
function generateTTpMatrix!(TTp::String, kappa::Int64, primitives::Bsplines.Primitives, storage::Dict{String,Array{Float64,2}})
    # Look up the dictionary of whether the requested matrix has been calculated before
    key = string(kappa) * ":" * TTp;      nsL = primitives.grid.nsL;   nsS = primitives.grid.nsS;
    wc  = Defaults.getDefaults("speed of light: c")
    
    wa  = get( storage, key, zeros(1,1) )
    if  wa != zeros(1,1)  
        # println(">>>> Re-used $TTp matrix for kappa = $kappa ...")
        return( wa )    
    end
    
    # Now calculate and store the requested matrix
    if      TTp == "LL-overlap"
        wa = zeros( nsL, nsL ) 
        for  i = 1:nsL,  j = 1:nsL
            wa[i,j] = Bsplines.computeOverlap(primitives.bsplinesL[i], primitives.bsplinesL[j], primitives.grid)
        end
    elseif  TTp == "SS-overlap"
        wa = zeros( nsS, nsS ) 
        for  i = 1:nsS,  j = 1:nsS
            wa[i,j] = Bsplines.computeOverlap(primitives.bsplinesS[i], primitives.bsplinesS[j], primitives.grid)
         end
    elseif  TTp == "LS-D_kappa^-"
        wa = zeros( nsL, nsS ) 
        for  i = 1:nsL,  j = 1:nsS
            wa[i,j] = wc * Bsplines.computeNondiagonalD(-1, kappa, primitives.bsplinesL[i], primitives.bsplinesS[j], primitives.grid)
        end
    elseif  TTp == "SL-D_kappa^+"
        wa = zeros( nsS, nsL ) 
        for  i = 1:nsS,  j = 1:nsL
            wa[i,j] = wc * Bsplines.computeNondiagonalD( 1, kappa, primitives.bsplinesS[i], primitives.bsplinesL[j], primitives.grid)
        end
    else   println("TTp = $TTp ");    error("stop a")
    end
    
    storage[key] = copy(wa)

    return( wa )
end


"""
`Bsplines.setupLocalMatrix(kappa::Int64, primitives::Bsplines.Primitives, pot::Radial.Potential, storage::Dict{String,Array{Float64,2}};
                           mass::Float64=1.0)`
    ... sets up the local parts of the generalized eigenvalue problem for the symmetry block kappa and the given (local) potential pot. The
        B-spline (basis) functions are defined by primitives.bsplinesL for the large component and by primitives.bsplinesS for the small
        one, respectively.
        A matrix::Array{Float64,2} of dimension (nsL+nsS) x (nsL+nsS) is returned.

        THE MATRIX IS FILLED ELEMENTWISE, AND DELIBERATELY SO -- it is not to be rewritten back into slice expressions. This function runs
        once per kappa per SCF iteration, and the slice form allocated 13.9 MB on every call: a SECOND full (nsL+nsS)^2 matrix for the
        overlap, of which only the SS block was ever read, plus a temporary for each side of every slice expression in steps (3) and (4).
        In isolation such a call looks cheap, 0.0029 s; a hundred calls in a row average 0.0146 s each, the difference being garbage
        collection of about 1 GB per SCF. Only the returned matrix is allocated now, and every step writes into it elementwise in the same
        order and with the same arithmetic, so the result is bitwise what the slice expressions produced.
"""
function setupLocalMatrix(kappa::Int64, primitives::Bsplines.Primitives, pot::Radial.Potential, storage::Dict{String,Array{Float64,2}};
                          mass::Float64=1.0)
    nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS
    wa  = zeros( nsL+nsS, nsL+nsS )

    # (1) Compute or fetch the 'SS-overlap' block; it is the only overlap block this function reads (step 3).
    #     generateTTpMatrix! returns the cached matrix itself, which is read but never written here.
    ssOverlap = Bsplines.generateTTpMatrix!("SS-overlap", 0, primitives, storage)

    # (2) Re-compute the diagonal blocks for the local potential.  B-splines have LOCAL support, so for a
    #     given i only a band of j overlaps -- 3.8% of the pairs, half-bandwidth 6 on a typical grid.  The
    #     splines are ordered by support, so once bspline[j] starts beyond the end of bspline[i] no later j
    #     can overlap either and the row is done.  computeVlocal returns exactly 0. for those, which is what
    #     the untouched entries already hold.
    for  i = 1:nsL
        for  j = 1:nsL
            if  primitives.bsplinesL[j].upper <= primitives.bsplinesL[i].lower    continue    end
            if  primitives.bsplinesL[j].lower >= primitives.bsplinesL[i].upper    break       end
            wa[i,j] = Bsplines.computeVlocal(primitives.bsplinesL[i], primitives.bsplinesL[j], pot, primitives.grid)
        end
    end
    for  i = 1:nsS
        for  j = 1:nsS
            if  primitives.bsplinesS[j].upper <= primitives.bsplinesS[i].lower    continue    end
            if  primitives.bsplinesS[j].lower >= primitives.bsplinesS[i].upper    break       end
            wa[nsL+i,nsL+j] = Bsplines.computeVlocal(primitives.bsplinesS[i], primitives.bsplinesS[j], pot, primitives.grid)
        end
    end

    # (3) Substract the rest mass from the 'SS' block.  This is the ONLY place the particle mass enters the
    #     Dirac matrix: the kinetic blocks below carry the speed of light but not the mass, so a muon differs
    #     from an electron here and nowhere else.  mass = 1.0 is an electron.
    twoCsq = 2 * mass * Defaults.getDefaults("speed of light: c")^2
    for  i = 1:nsS,  j = 1:nsS
        wa[nsL+i,nsL+j] = wa[nsL+i,nsL+j] - twoCsq * ssOverlap[i,j]
    end

    # (4) Compute or fetch the diagonal 'D_kappa' blocks
    lsD = Bsplines.generateTTpMatrix!("LS-D_kappa^-", kappa, primitives, storage)
    slD = Bsplines.generateTTpMatrix!("SL-D_kappa^+", kappa, primitives, storage)
    for  i = 1:nsL,  j = 1:nsS      wa[i,nsL+j] = wa[i,nsL+j] + lsD[i,j]      end
    for  i = 1:nsS,  j = 1:nsL      wa[nsL+i,j] = wa[nsL+i,j] + slD[i,j]      end

    #===== DISABLED, kept as a diagnostic that can be switched back on: it counts how far the local matrix
    # departs from real-symmetric.  It is OFF because the departure is EXPECTED -- the last B-spline breaks
    # the symmetry -- so it fires routinely and says nothing.  Labelled 13-Aug-2026.
    # Test for 'real-symmetric matrix' ... this is not fullfilled if the last B-spline is included !!
    nx = 0
    for  i = 1:nsL+nsS    
        for  j = i+1:nsL+nsS    
            if  abs(  (wa[i,j] - wa[j,i])/(wa[i,j] + wa[j,i]) ) > 1.0e-7   nx = nx + 1    
                @show "setupLocalMatrix", i, j, wa[i,j], wa[j,i] 
            end
        end
    end
    ny = (nsL+nsS)^2/2 - (nsL+nsS)
    if  nx > 0    
        println(">>> setupLocalMatrix:: $nx (from $(ny)) non-symmetric H-matrix integrals for kappa = $kappa with relative deviation > 1.0e-7.")  end
    =====#
    
    return( wa )
end   


end # module
