# SUPERSEDED 04-Oct-2026, THE SAME DAY IT WAS WRITTEN.  Its premise -- the first line below -- is FALSE:
# gamma_inf is NOT box-dependent.  The box scan that motivated this probe used `Basics.recommendedGrid`'s
# default hp = rbox/300, which keeps the spline basis at nsL ~ 95 however large the box, so what was varying was
# the SPLINE DENSITY and not the box at all.  At fixed density a tenfold box gives the same gamma to five
# figures.  See `tools/probe-nuclearShieldingDensity.jl`, which carries the prescription that replaced this one.
# KEPT because the two cancellation ratios it measures became `NuclearShielding.checkCancellation`, and because
# its own rejected candidate (the locality ratio) is worth not retrying.
#
# probe-nuclearShieldingBoxCriterion.jl -- WHY gamma_inf RUNS AWAY WITH THE BOX, AND WHAT WOULD WARN US EARLY.
#
# Measured 04-Oct-2026: gamma_inf is not converged in the radial box for ANY ion, and for Th(4+) it changes sign
# between 14 and 61 a.u.  Before a grid can be prescribed, the mechanism has to be named.  The hypothesis this
# probe tests is the one the pattern suggests -- the damage grows with box size and is worst for the ions whose
# OUTERMOST occupied orbital lies closest to threshold:
#
#   gamma is a sum of  <a|r^2|v><v|r^-3|a> / (eps_a - eps_v)  over the B-spline pseudo-spectrum.  Enlarging the
#   box does two things at once.  It adds pseudo-states just ABOVE eps_a for a weakly bound outer orbital, so
#   1/(eps_a - eps_v) grows without bound; and the driving operator is r^2, which weights exactly the large-r
#   region those states live in.  A deeply bound orbital is immune -- nothing crowds its energy -- which is why
#   the light, highly charged ions misbehave least.
#
# So the probe prints, per channel and per box: the smallest |eps_a - eps_v| actually used, the single largest
# term, and how much of the channel that one term carries.  IF ONE TERM CARRIES THE CHANNEL, the sum is not a
# sum any more and the number is an accident of where the nearest pseudo-state happened to fall.
using JenaAtomicCalculator, Printf, LinearAlgebra
const JAC = JenaAtomicCalculator
const NS  = JenaAtomicCalculator.NuclearShielding

function analyse(name, Z, confs, rbox)
    conf = Configuration(confs);   nm = Nuclear.Model(Z)
    grid = rbox > 0 ? Basics.recommendedGrid([conf], nm; rbox=rbox, printout=false) :
                      Basics.recommendedGrid([conf], nm; printout=false)
    setDefaults("standard grid", grid)
    asf  = AsfSettings(AsfSettings(); scField=Basics.DFSField(1.0))
    tmp  = tempname()
    mp   = open(tmp,"w") do io;  redirect_stdout(io) do
               perform(Atomic.Computation(Atomic.Computation(); name=name, grid=grid, nuclearModel=nm,
                       configs=[conf], asfSettings=asf); output=true)["multiplet:"]  end  end
    rm(tmp, force=true)
    level = mp.levels[1];   basis = level.basis
    pot   = Basics.add( Nuclear.nuclearPotential(nm, grid), Basics.computePotential(Basics.DFSField(1.0), grid, basis) )
    prims = Bsplines.generatePrimitives(grid);   storage = Dict{String,Array{Float64,2}}()
    nsL, nsS = grid.nsL, grid.nsS
    ov = zeros(nsL+nsS, nsL+nsS)
    ov[1:nsL,1:nsL]                 = Bsplines.generateTTpMatrix!("LL-overlap", 0, prims, storage)
    ov[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = Bsplines.generateTTpMatrix!("SS-overlap", 0, prims, storage)

    occ = Subshell[];   occE = Dict{Int64,Array{Float64,1}}()
    for  sh in basis.subshells
        haskey(basis.orbitals, sh) || continue
        Basics.computeMeanSubshellOccupation(sh, [level]) < 1.0e-6 && continue
        push!(occ, sh);   push!(get!(occE, sh.kappa, Float64[]), basis.orbitals[sh].energy)
    end
    outermost = occ[argmax([basis.orbitals[s].energy for s in occ])]
    eOuter    = basis.orbitals[outermost].energy

    model = NS.UniformField();   applied = NS.appliedGradient(model)
    gamma = 0.;   worst = ("", 0.0, 0.0, 0.0, 0.0)      # channel, |term|, denom, share, contribution
    for  sh in occ
        orb = basis.orbitals[sh]
        bD  = NS.projectOperator(r -> r*r,        orb, prims, grid)
        bO  = NS.projectOperator(r -> 1/(r*r*r),  orb, prims, grid)
        for kp in NS.allowedKappas(sh.kappa)
            w  = AngularMomentum.CL_reduced_me(sh, 2, Subshell(9,kp))^2
            mA = Bsplines.setupLocalMatrix(kp, prims, pot, storage)
            eg = Bsplines.diagonalizeLocalMatrix(kp, mA, ov, prims)
            terms = Float64[];  dens = Float64[]
            for (i,vec) in enumerate(eg.vectors)
                ev = eg.values[i];   isO = false
                for eO in get(occE,kp,Float64[])   abs(ev-eO) < 1.0e-6*max(1.0,abs(eO)) && (isO=true;break)  end
                isO && continue
                abs(orb.energy-ev) < 1.0e-10 && continue
                push!(terms, dot(vec,bD)*dot(vec,bO)/(orb.energy-ev));   push!(dens, abs(orb.energy-ev))
            end
            isempty(terms) && continue
            contrib = NS.GAMMA_PREFACTOR*(2.0/applied)*w*sum(terms)
            gamma  += contrib
            k = argmax(abs.(terms))
            share = abs(terms[k])/max(1e-30, sum(abs.(terms)))
            if  abs(w*terms[k]) > worst[2]
                worst = ("$(sh) -> $(kp)", abs(w*terms[k]), dens[k], share, contrib)
            end
        end
    end
    return( grid.tL[end], gamma, string(outermost), eOuter, worst )
end

println("\n  WHAT RUNS AWAY, AND WHAT WOULD HAVE WARNED US")
println("  " * "="^118)
@printf("  %-7s %8s %11s %-9s %11s   %-16s %11s %11s %8s\n",
        "ion", "r_box", "gamma", "outermost", "eps_outer", "worst channel", "|term|", "denom", "share")
println("  " * "-"^118)
for  (name, Z, confs)  in  (("Y^3+", 39.0, "[Kr]"), ("Ba^2+", 56.0, "[Xe]"), ("Th^4+", 90.0, "[Rn]"))
    for rb in (0.0, 30.0, 60.0, 90.0)
        rbox, g, out, eo, w = analyse(name, Z, confs, rb)
        @printf("  %-7s %8.1f %11.2f %-9s %11.3f   %-16s %11.3f %11.2e %7.0f %%\n",
                rb == 0.0 ? name : "", rbox, g, out, eo, w[1], w[2], w[3], 100*w[4])
        flush(stdout)
    end
    println()
end
println("  " * "-"^118)
println("  'share' is how much of its own channel the single largest term carries.  A sum in which one term")
println("  carries most of the channel is not converged, whatever its value happens to be.")
