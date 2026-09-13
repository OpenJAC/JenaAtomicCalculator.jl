#
# probe-lsBlockDiagonality.jl   --   proves, WITHOUT an external oracle, which two-particle coefficients are
#                                    wrong.  Evidence for priority item 27.
#
# The electron-electron interaction commutes with L and with S, so in LS coupling the 3s3p -- 3p3d block must be
# STRICTLY DIAGONAL in (L,S).  Build the jj matrix from a code's own coefficients, transform with the jj->LS 9-j,
# and read off the elements that must vanish.  Measured 12-Sep-2026:
#
#                  J = 1                     J = 2
#     GRASP2018    3.3e-16  obeys L,S        1.4e-16  obeys L,S
#     JAC          9.5e-01  VIOLATES         2.4e-01  VIOLATES
#
# so JAC's coefficients for CSF pairs differing by one electron across two open shells are the ones in error.
# This does not assume GRASP is right; it only assumes V_ee commutes with L and S.  The 9-j is self-checked by
# the orthogonality of the transformation, 5e-16.
#
# TWO CONSTRAINTS ARE ESSENTIAL, AND WITHOUT EITHER THE TEST CONDEMNS BOTH CODES -- which is how a broken version
# of this file announces itself, and both mistakes were made here before the test worked:
#   (1) THE LS LIMIT.  R^k must depend on (n,l) and NOT on j, because L and S are good quantum numbers only when
#       the spin-orbit partners share a radial function.  Keeping 3p- and 3p radially distinct is the
#       relativistic case, where L and S are not good and the criterion has nothing to say.
#   (2) THE SLATER SYMMETRY  R^k(abcd) = R^k(badc) = R^k(cdab) = R^k(dcba).
# Random draws of the radial integrals make the test structural: a cancellation holding for one draw will not
# hold for ten.
#
# INPUT is the two coefficient dumps produced by `tools/probe-twoParticleVsGrasp.jl` (JAC) and
# `tools/grasp-drivers/angdump2.f90` (GRASP); edit the two paths below to point at them.
#




using JenaAtomicCalculator, Printf, Random, LinearAlgebra
const AM = JenaAtomicCalculator.AngularMomentum

# 9-j from the standard single sum over 6-j's
function ninej(a,b,c, d,e,f, g,h,j)
    s = 0.0
    xlo = max(abs(a-j), abs(b-f), abs(d-h));   xhi = min(a+j, b+f, d+h)
    x = xlo
    while x <= xhi + 1e-9
        s += (-1)^round(Int, 2x) * (2x+1) *
             AM.Wigner_6j(a,b,c, f,j,x) * AM.Wigner_6j(d,e,f, b,x,h) * AM.Wigner_6j(g,h,j, x,a,d)
        x += 1.0
    end
    return( s )
end

# transformation |(l1 l2)L (ss)S; J>  ->  |(l1 s)j1 (l2 s)j2; J>
tmat(l1, l2, jjs, lss, J) =
    [ sqrt((2L+1)*(2S+1)*(2j1+1)*(2j2+1)) * ninej(l1, 0.5, j1, l2, 0.5, j2, L, S, J)
      for (j1,j2) in jjs, (L,S) in lss ]

# --- the bases, in GRASP's CSF order
bas = Dict(
 0 => (jjA=[(0.5,0.5)],                          lsA=[(1.0,1.0)],
       jjB=[(1.5,1.5)],                          lsB=[(1.0,1.0)],
       idA=[1], idB=[2]),
 1 => (jjA=[(0.5,1.5),(0.5,0.5)],                lsA=[(1.0,0.0),(1.0,1.0)],
       jjB=[(1.5,2.5),(1.5,1.5),(0.5,1.5)],      lsB=[(1.0,0.0),(1.0,1.0),(2.0,1.0)],
       idA=[3,4], idB=[5,6,7]),
 2 => (jjA=[(0.5,1.5)],                          lsA=[(1.0,1.0)],
       jjB=[(1.5,2.5),(1.5,1.5),(0.5,2.5),(0.5,1.5)],
                                                 lsB=[(2.0,0.0),(1.0,1.0),(2.0,1.0),(3.0,1.0)],
       idA=[8], idB=[9,10,11,12]) )

function readcoeff(fn, remap)
    d = Dict{Tuple{Int,Int}, Vector{Tuple{Int,Int,Int,Int,Int,Float64}}}()
    for line in eachline(fn)
        p = split(line);   length(p) == 8 || continue
        all(x -> occursin(r"^-?\d+$", x), p[1:7]) || continue
        ic, ir, a, b, c, dd, k = parse.(Int, p[1:7]);   v = parse(Float64, p[8])
        ic = remap(ic);  ir = remap(ir)
        push!(get!(d, (ic,ir), []), (k,a,b,c,dd,v))
    end
    return( d )
end

S = "/tmp/claude-1000/-home-fritzsch-fri-JAC-jl-examples/d2e5b6ed-950b-4f37-9076-3b4e7a4d92a8/scratchpad/"
j2g = Dict(1=>3,2=>8,3=>1,4=>4,5=>5,6=>9,7=>2,8=>6,9=>10,10=>11,11=>7,12=>12)
G = readcoeff(S*"grasp/run/grasp2p.log", identity)
J = readcoeff(S*"jac2p.log", i -> j2g[i])

Random.seed!(20260912)
println("\nThe 3s3p -- 3p3d block in LS coupling.  V_ee commutes with L and S, so every element with")
println("(L,S) != (L',S') MUST vanish.  Ten random draws of the radial integrals R^k(abcd).\n")
for  Jv in [1, 2]
    b = bas[Jv]
    TA = tmat(0.0, 1.0, b.jjA, b.lsA, Float64(Jv))
    TB = tmat(1.0, 2.0, b.jjB, b.lsB, Float64(Jv))
    @printf("J = %d :  orthogonality |TA'TA - I| = %.2e, |TB'TB - I| = %.2e\n", Jv,
            maximum(abs.(TA'*TA - I(size(TA,2)))), maximum(abs.(TB'*TB - I(size(TB,2)))))
    for (name, D) in [("GRASP2018", G), ("JAC", J)]
        worst = 0.0
        for trial = 1:10
            # THE SLATER INTEGRAL IS SYMMETRIC: R^k(abcd) = R^k(badc) = R^k(cdab) = R^k(dcba).  Giving each
            # label tuple its own random value breaks that and destroys the LS structure by construction --
            # it made BOTH codes appear to violate L,S on the first run of this test.
            # TWO CONSTRAINTS, AND THE TEST IS MEANINGLESS WITHOUT EITHER.
            # (1) THE LS LIMIT.  L and S are good quantum numbers only when the spin-orbit partners share a
            #     radial function, so R^k must depend on (n,l) and NOT on j:  subshell 1 = 3s -> l 0,
            #     2,3 = 3p-,3p -> l 1,  4,5 = 3d-,3d -> l 2.  Keeping them distinct is the relativistic case,
            #     where L and S are NOT good and the test has nothing to say -- that is what made both codes
            #     look like violators on the first two runs of this file.
            # (2) THE SLATER SYMMETRY:  R^k(abcd) = R^k(badc) = R^k(cdab) = R^k(dcba).
            ell = Dict(1=>0, 2=>1, 3=>1, 4=>2, 5=>2)
            Rv  = Dict{NTuple{5,Int},Float64}()
            function getR(k,a,bb,c,d)
                (la,lb,lc,ld) = (ell[a], ell[bb], ell[c], ell[d])
                key = (k, min((la,lb,lc,ld), (lb,la,ld,lc), (lc,ld,la,lb), (ld,lc,lb,la))...)
                return( get!(Rv, key, randn()) )
            end
            M = zeros(length(b.idA), length(b.idB))
            for (ii, ia) in enumerate(b.idA), (jj, jb) in enumerate(b.idB)
                for (k,a,bb,c,d,v) in get(D, (ia,jb), [])
                    M[ii,jj] += v * getR(k,a,bb,c,d)
                end
            end
            MLS = TA' * M * TB
            for  (p, (L,Sp)) in enumerate(b.lsA), (q, (Lq,Sq)) in enumerate(b.lsB)
                (L == Lq && Sp == Sq) && continue
                worst = max(worst, abs(MLS[p,q]))
            end
        end
        @printf("          %-10s  worst FORBIDDEN element over 10 draws = %.4e   %s\n",
                name, worst, worst < 1e-10 ? "<-- obeys L,S" : "<-- VIOLATES L,S")
    end
    println()
end
