# Writes JAC orbitals in the text format the relativistic distorted-wave (RDW) Fortran code reads (routine INPUT, iold = 0):
#
#     G92RWF
#     n  kappa  -eps  npoint  p0            format (2i4,d24.16,i4,d24.16)
#     r  P  Q                                format (3d24.16), npoint lines, the first at r = 0
#
# on the mesh that code builds itself, R(i) = 2e-6 (exp(0.05 i) - 1): the r column is read and then DISCARDED, so the values
# must sit exactly on that mesh. JAC's P, Q are interpolated there by cubic Hermite polynomials using JAC's own P', Q'.
#
# Used for the Na 3s -> 3p comparison of the non-local exchange (work/rdw-na): both codes then run on the SAME orbitals, so the
# atomic structure drops out of the comparison.
#
#   julia --project=. tools/export-rdw-orbitals.jl <output file>

using JenaAtomicCalculator, Printf
const JAC = JenaAtomicCalculator

function hermite(r::Vector{Float64}, y::Vector{Float64}, yp::Vector{Float64}, x::Float64)
    x <= r[1] && return( y[1] )
    k = searchsortedlast(r, x);   k >= length(r) && return( 0. )
    h = r[k+1] - r[k];   t = (x - r[k]) / h
    return( (2t^3 - 3t^2 + 1) * y[k] + (t^3 - 2t^2 + t) * h * yp[k] + (-2t^3 + 3t^2) * y[k+1] + (t^3 - t^2) * h * yp[k+1] )
end

"""
`writeRdwOrbitals(file, orbitals)`  ... orbitals is a list of JAC Orbitals in the order the RDW input expects.
"""
function writeRdwOrbitals(file::String, orbitals, grid::Radial.Grid; nmax::Int=390, cut::Float64=1.0e-14)
    R(i) = 2.0e-6 * (exp(0.05 * i) - 1.0)
    open(file, "w") do io
        println(io, " G92RWF")
        for  o in orbitals
            n  = min(length(o.P), grid.NoPoints);   r = grid.r[1:n]
            P  = [hermite(r, o.P[1:n], o.Pprime[1:n], R(i)) for i = 1:nmax]
            Q  = [hermite(r, o.Q[1:n], o.Qprime[1:n], R(i)) for i = 1:nmax]
            np = something(findlast(i -> abs(P[i]) > cut, 1:nmax), nmax) + 1
            np = min(np, nmax)
            @printf(io, "%4d%4d%24.16E%4d%24.16E\n", o.subshell.n, o.subshell.kappa, -o.energy, np + 1, 0.0)
            @printf(io, "%24.16E%24.16E%24.16E\n", 0.0, 0.0, 0.0)
            for  i = 1:np    @printf(io, "%24.16E%24.16E%24.16E\n", R(i), P[i], Q[i])    end
        end
    end
    return( nothing )
end

if  abspath(PROGRAM_FILE) == @__FILE__
    out  = length(ARGS) > 0 ? ARGS[1] : "naorb.in"
    grid = Radial.Grid(Radial.Grid(true), rnt = 2.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 60.0)
    nm   = Nuclear.Model(11.)
    lev(c) = redirect_stdout(devnull) do
        perform(Atomic.Computation(Atomic.Computation(), name="Na", grid=grid, nuclearModel=nm, configs=[Configuration(c)]),
                output=true)["multiplet:"].levels[1]
    end
    ls = lev("[Ne] 3s");   lp = lev("[Ne] 3p")
    o(l, s) = l.basis.orbitals[Subshell(s)]
    orbs = [o(ls,"1s_1/2"), o(ls,"2s_1/2"), o(ls,"2p_1/2"), o(ls,"2p_3/2"), o(ls,"3s_1/2"), o(lp,"3p_1/2"), o(lp,"3p_3/2")]
    b = orbs[1];   i = findfirst(x -> x > 1.0e-4, grid.r)
    @printf("sign check 1s at r = %.2e:  P = % .3e  Q = % .3e   (GRASP: P > 0, Q < 0)\n", grid.r[i], b.P[i], b.Q[i])
    writeRdwOrbitals(out, orbs, grid)
    @printf("excitation energies [eV] (JAC, these orbitals' levels): 3p_1/2 %.6f   3p_3/2 from multiplet\n",
            Defaults.convertUnits("energy: from atomic to eV", lp.energy - ls.energy))
    println("wrote ", out)
end
