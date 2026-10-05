#
# tools/probe-hfDrivenCaptureScaling.jl
#
# THE n-SCALING OF A HYPERFINE-DRIVEN DIELECTRONIC CAPTURE, AND THE EXTRAPOLATION TO n = 133.
#
# PRIORITY ITEM 34.  A dielectronic resonance whose only excitation is the core's own hyperfine flip exists in
# H-like Bi only where the Rydberg electron is bound by LESS than the hyperfine splitting:  JAC gives
# E_HFS(1s, Bi82+) = 5.1958 eV against the measured 5.084, the captured electron sees q = 82, and so
#
#       n_min = sqrt( q^2 / 2 E_HFS ) = 133 .
#
# A 133s orbital reaches ~244 a.u., and a grid that must ALSO resolve 1s at Z = 83 cannot carry it at the spline
# density Rule 46 asks for.  So the rate at n = 133 is obtained by measuring reachable n and extrapolating, which
# the maintainer approved on 05-Oct-2026.
#
# HOW REACHABLE n ARE USED AT ALL.  At those n the channel is CLOSED (E_res < 0) and `determineHfCaptureLines`
# rightly discards it.  But the capture rate A(n; eps) is a well-defined function of n and eps -- the resonance
# condition only decides WHICH eps is realised -- so the energy is supplied directly to
# `DielectronicRecombination.computeCaptureAmplitudes` instead.  What is computed here is the ELECTRONIC capture
# rate;  the hyperfine recoupling coefficient multiplying it is n-INDEPENDENT, so it carries the n-scaling
# unchanged.
#
# AND THE EXPONENT DEPENDS ON eps, WHICH NEARLY PRODUCED AN ANSWER A FACTOR OF FOUR TOO LOW.  At eps = 50 eV the
# local exponent CLIMBS with n -- 1.73, 1.82, 1.97, 2.17, 2.36, 2.49 over n = 12..55, fitting 3 - p = 6.1 n^-0.59
# -- and an extrapolation built on that drift gives 2.7e-04 a.u.  At the physical energy it is FLAT instead.  The
# scaling must be measured AT the energy wanted;  a drift fitted at a convenient energy does not transfer.  The
# two regimes differ in whether the continuum wavelength (3.3 a.u. at 50 eV, 52 at 0.2 eV, 150 at 0.024 eV) is
# small or comparable to the Rydberg extent (88 a.u. at n = 60, 244 at n = 100).
#
# WHAT THIS MEASURED, 05-Oct-2026, DFS mean field, Coulomb capture operator, eps = 0.2 eV where A(eps) has
# saturated (1.1248, 1.1379, 1.1406 x its 50 eV value at 5, 1, 0.2 eV -- no threshold suppression, as an
# attractive Coulomb field requires):
#
#       n =  20   2.661489e-02        n =  60   4.320246e-03   p = 1.670
#       n =  24   1.972982e-02        n =  80   2.677647e-03   p = 1.663
#       n =  28   1.528247e-02        n = 100   1.838243e-03   p = 1.686
#       n =  34   1.110026e-02
#       n =  40   8.479111e-03        p = 1.64 to 1.69 across n = 20..100, flat to +/- 0.02
#       n =  50   5.857671e-03
#
#       ==>  A_capture(1s 133s, eps -> 0) = 1.14e-03 a.u. = 4.7e+13 1/s
#
# THE GRID CONTROL IS EXACT, and without it none of the above would be evidence: A(60) = 4.320246e-03 on a
# 300 a.u. box with 1287 splines against 4.320206e-03 on a 120 a.u. box with 631 -- ratio 1.0000 to five figures.
#
# ONE QUESTION IS STILL OPEN, and it is physics rather than numerics: p = 1.67 is NOT the 1/n^3 that a genuine
# two-electron Auger rate shows.  That is consistent with the mechanism -- the core is electronically UNCHANGED,
# so the amplitude is a monopole-screened one-body transition eps s -> ns whose integrand lives over the whole
# Rydberg orbital rather than near the core -- but it has not been confirmed against an independent value.
#
# Usage:   julia --project=. tools/probe-hfDrivenCaptureScaling.jl            # n = 20..60, 120 a.u. box
#          julia --project=. tools/probe-hfDrivenCaptureScaling.jl --long     # adds n = 80, 100, 300 a.u. box
#
using JenaAtomicCalculator, Printf
const DR = JenaAtomicCalculator.DielectronicRecombination

const REFERENCE = Dict(20 => 2.661489e-02, 24 => 1.972982e-02, 28 => 1.528247e-02, 34 => 1.110026e-02,
                       40 => 8.479111e-03, 50 => 5.857671e-03, 60 => 4.320246e-03, 80 => 2.677647e-03,
                       100 => 1.838243e-03)


"""
`captureRateAt(n::Int64, epsEV::Float64, grid::Radial.Grid, nm::Nuclear.Model)`  
    ... computes the electronic capture rate for  1s + e(eps) --> 1s ns  in H-like bismuth, summed over the
        intermediate levels and the capture partial waves, with the electron energy supplied DIRECTLY rather than
        taken from the level difference; a rate::Float64 in atomic units is returned.

        Supplying the energy is what makes a reachable n usable: at n well below 133 the hyperfine-driven channel
        is closed and `determineHfCaptureLines` discards it, while A(n; eps) remains a perfectly well-defined
        function whose n-dependence is what this probe is after.
"""
function captureRateAt(n::Int64, epsEV::Float64, grid::Radial.Grid, nm::Nuclear.Model)
    asfSettings = AsfSettings(AsfSettings(); scField = Basics.DFSField())
    drSettings  = DR.Settings(DR.Settings(); multipoles=[Basics.E1], gauges=[Basics.UseCoulomb],
                              calcOnlyPassages=true)
    mult(confs, name) = redirect_stdout(devnull) do
        perform(Atomic.Computation(Atomic.Computation(); name=name, grid=grid, nuclearModel=nm,
                configs=confs, asfSettings=asfSettings); output=true)["multiplet:"]
    end
    iMultiplet = mult([Configuration("1s")],       "Bi82+")
    mMultiplet = mult([Configuration("1s $(n)s")], "Bi81+ 1s$(n)s")
    eps        = Defaults.convertUnits("energy: to atomic", epsEV)
    nrContinuum = Continuum.gridConsistency(eps, grid)
    nuclearPot  = Nuclear.nuclearPotential(nm, grid);   primitives = Bsplines.generatePrimitives(grid)
    rate = 0.
    for  iLevel in iMultiplet.levels,  mLevel in mMultiplet.levels
        pws = DR.determineCaptureChannels(mLevel, iLevel, drSettings)
        if  isempty(pws)    continue    end
        cLine = DR.CaptureLine(iLevel, mLevel, eps, 0., 0., EmProperty(0., 0.), EmProperty(0., 0.), pws)
        rate  = rate + DR.computeCaptureAmplitudes(cLine, nm, grid, nrContinuum, drSettings;
                                                   nuclearPot=nuclearPot, primitives=primitives).captureRate
    end

    return( rate )
end


"""
`scanOverN(ns::Array{Int64,1}, epsEV::Float64, rbox::Float64, hp::Float64, nm::Nuclear.Model)`  
    ... measures the capture rate over the given principal quantum numbers on one grid and prints the local
        exponent p of each consecutive pair, beside the value this probe recorded on 05-Oct-2026; a tuple
        (ns::Array{Int64,1}, rates::Array{Float64,1}) is returned.
"""
function scanOverN(ns::Array{Int64,1}, epsEV::Float64, rbox::Float64, hp::Float64, nm::Nuclear.Model)
    grid = Radial.Grid(Radial.Grid(false), rnt=1.0e-6, h=5.0e-2, hp=hp, rbox=rbox)
    setDefaults("standard grid", grid)
    @printf("\n  grid: rbox = %.1f a.u., %d points, nsL = %d, %.2f splines per a.u.;   eps = %.3f eV\n",
            grid.tL[end], grid.NoPoints, grid.nsL, grid.nsL/grid.tL[end], epsEV)
    @printf("  %5s %12s %18s %10s %16s\n", "n", "r_max(ns)", "A_capture [a.u.]", "local p", "vs 05-Oct")
    println("  " * "-"^68)
    rates = Float64[]
    for  (j, n) in enumerate(ns)
        rate = captureRateAt(n, epsEV, grid, nm);    push!(rates, rate)
        p    = j == 1 ? NaN : -log(rates[j]/rates[j-1]) / log(ns[j]/ns[j-1])
        ref  = get(REFERENCE, n, NaN)
        @printf("  %5d %12.1f %18.6e %10.3f %16s\n", n, 2.0*n^2/82.0, rate, p,
                isnan(ref) ? "--" : @sprintf("%.5f", rate/ref))
        flush(stdout)
    end
    println("  " * "-"^68)

    return( (ns, rates) )
end


"""
`extrapolate(ns::Array{Int64,1}, rates::Array{Float64,1}, nTarget::Int64)`  
    ... extrapolates the measured power law to nTarget and prints the value together with the bracket that the
        exponent's own spread implies; a value::Float64 in atomic units is returned.
"""
function extrapolate(ns::Array{Int64,1}, rates::Array{Float64,1}, nTarget::Int64)
    ps = [ -log(rates[j]/rates[j-1]) / log(ns[j]/ns[j-1])  for j = 2:length(ns) ]
    pLo = minimum(ps);    pHi = maximum(ps);    pMid = (pLo + pHi)/2
    lnr = log(nTarget/ns[end]);    auToPerSec = 4.1341373e16
    aMid = rates[end]*exp(-pMid*lnr);   aLo = rates[end]*exp(-pHi*lnr);   aHi = rates[end]*exp(-pLo*lnr)
    @printf("\n  the exponent runs %.3f to %.3f over n = %d..%d, so extrapolating to n = %d (a factor %.2f):\n",
            pLo, pHi, ns[1], ns[end], nTarget, nTarget/ns[end])
    @printf("      A(%d) = %.4e a.u. = %.3e 1/s,   bracket %.3e .. %.3e 1/s\n",
            nTarget, aMid, aMid*auToPerSec, aLo*auToPerSec, aHi*auToPerSec)

    return( aMid )
end


nm   = Nuclear.Model(Nuclear.Model(83.); spinI=AngularJ64(9//2), mu=4.1106, Q=-0.516)
long = "--long" in ARGS
println("\n  HYPERFINE-DRIVEN DIELECTRONIC CAPTURE IN H-LIKE Bi:  the n-scaling and n = 133")
println("  " * "="^88)
ns, rates = scanOverN([20, 24, 28, 34, 40, 50, 60], 0.2, 120.0, 3.0e-2, nm)
extrapolate(ns, rates, 133)
if  long
    println("\n  THE SAME ON A 300 a.u. BOX, which reaches n = 100 and repeats n = 60 as the grid CONTROL:")
    nsL, ratesL = scanOverN([60, 80, 100], 0.2, 300.0, 3.5e-2, nm)
    extrapolate(nsL, ratesL, 133)
    @printf("\n  GRID CONTROL: A(60) = %.6e here against %.6e on the 120 a.u. box  (ratio %.5f)\n",
            ratesL[1], rates[end], ratesL[1]/rates[end])
end
println("\n  " * "="^88)
