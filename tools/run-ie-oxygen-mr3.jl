# Electron-impact excitation of neutral O from the 2p^4 3P_J ground levels to 2p^3 3s 3S and 2p^3 3d 3D (and the unresolvable 3d 5D),
# fine-structure resolved, with the structure from the 79-CSF GRASP2018 expansion MR_3 (1s core; 2s, 2p, 3s, 3p, 3d; Breit + QED).
# Every J_i -> J_f line is computed; the term cross sections that the (unresolved) experiments measure are then formed by averaging over
# the initial 3P_J with weights 2J+1 and summing over the final J of each term.
#
# Environment (all optional):
#   GDIR      directory with MR_3.c, MR_3.w, MR_3.cm     (default <JAC-root>/O_I/MR_3)
#   ENERGIES  comma-separated electron energies in eV
#   EXCH      comma-separated exchange models: FurnessMcCarthy, nonlocal, none
#   DCSE      energies (eV) at which the DCS is computed  (default 30,50,100, where DCS measurements exist; "all" for every energy)
#   CBKAPPA   upper |kappa| of the Coulomb-Bethe partial waves (default 100)
#   LINES     "i:f,i:f" level-index pairs instead of all (convergence tests)
#   MAXKAPPA  upper |kappa| of the distorted partial waves (default 20)
#   TOPUP     true: geometric top-up of the integrated cross section beyond the last partial wave (the DCS is not topped up)
#   SHIFT     "term=eV,term=eV": energy added to the outgoing electron of that term's lines, e.g. 3s3S=1.101,3d3D=1.004
#             moves the MR_3 thresholds (10.62, 13.09 eV) onto the measured ones (9.521, 12.087 eV)
#   OUT       output prefix                               (default oxygen-mr3)
#
# Run from <JAC-root>:  OPENBLAS_NUM_THREADS=2 julia --project=. -t 8 tools/run-ie-oxygen-mr3.jl

using JenaAtomicCalculator, Printf, Dierckx
const JAC = JenaAtomicCalculator
using JenaAtomicCalculator: ImpactExcitation, ManyElectron

setDefaults("unit: energy", "eV")
setDefaults("method: continuum, Galerkin")
setDefaults("method: normalization, Alok")

const GDIR     = get(ENV, "GDIR", joinpath(@__DIR__, "..", "O_I", "MR_3"))
const ENERGIES = parse.(Float64, split(get(ENV, "ENERGIES", "15,20,25,30,40,50,100,200"), ","))
const MODELS   = Symbol.(split(get(ENV, "EXCH", "FurnessMcCarthy,nonlocal"), ","))
const OUT      = get(ENV, "OUT", "oxygen-mr3")
const DCSE     = get(ENV, "DCSE", "30,50,100") == "all" ? ENERGIES : parse.(Float64, split(get(ENV, "DCSE", "30,50,100"), ","))
const CBKAPPA  = parse(Int, get(ENV, "CBKAPPA", "100"))
const MAXKAPPA = parse(Int, get(ENV, "MAXKAPPA", "20"))
const TOPUP    = get(ENV, "TOPUP", "false") == "true"
const SHIFT    = Dict(String(split(x, "=")[1]) => parse(Float64, split(x, "=")[2]) for x in split(get(ENV, "SHIFT", ""), ",") if !isempty(x))
const GRID     = Radial.Grid(Radial.Grid(true), rnt = 2.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 45.0)
const ANGLES   = collect(0.0:1.0:180.0)

"""
`readMR3(dir::String)`  ... the GRASP multiplet of `dir`; the CSF file is read through a copy whose one-subshell core line is padded to
    five columns, since `Basics.readCslFileGrasp92` slices that line in fixed 5-character fields. A multiplet::Multiplet is returned.
"""
function readMR3(dir::String)
    tmp = mktempdir()
    lines = readlines(joinpath(dir, "MR_3.c"))
    lines[2] = rpad(lines[2], 5 * cld(length(lines[2]), 5))
    write(joinpath(tmp, "MR_3.c"), join(lines, "\n") * "\n")
    mp = cd(tmp) do
        redirect_stdout(devnull) do
            Basics.readFilesGrasp18(GRID, joinpath(tmp, "MR_3.c"), joinpath(dir, "MR_3.w"), joinpath(dir, "MR_3.cm"))
        end
    end
    return( mp )
end


"""
`correctOrbitals(mp::Multiplet)`  ... a copy of the multiplet whose bound orbitals carry isBound = true, a negative energy and derivatives
    of P and Q taken from a spline through the tabulated functions; a multiplet::Multiplet is returned.
"""
function correctOrbitals(mp::Multiplet)
    b    = mp.levels[1].basis
    r    = GRID.r
    orbs = Dict{Subshell,Radial.Orbital}()
    for  (sh, o) in b.orbitals
        Pp = Dierckx.derivative(Dierckx.Spline1D(r, o.P; k=3), r)
        Qp = Dierckx.derivative(Dierckx.Spline1D(r, o.Q; k=3), r)
        orbs[sh] = Radial.Orbital(sh, true, o.useStandardGrid, -abs(o.energy), o.P, o.Q, Pp, Qp, o.grid)
    end
    basis  = ManyElectron.Basis(true, b.NoElectrons, b.subshells, b.csfs, b.coreSubshells, orbs)
    levels = [ManyElectron.Level(l.J, l.M, l.parity, l.index, l.energy, l.relativeOcc, l.hasStateRep, basis, l.mc) for l in mp.levels]
    return( ManyElectron.Multiplet(mp.name, levels) )
end


mp  = correctOrbitals(readMR3(GDIR))
E0  = mp.levels[1].energy
eV(i) = Defaults.convertUnits("energy: from atomic to eV", mp.levels[i].energy - E0)

# Terms by level index in the energy-ordered MR_3 multiplet (J and energy checked on reading, see below).
initial = [1, 2, 3]
# 3d 5D (levels 14-17) lies 20 meV below 3d 3D and is unresolved in the experiments, but at 50 eV (Furness-McCarthy) it carries only
# 1.5 % of the 3d 3D cross section; set WITH5D=true to include it.
terms   = [("3s3S", [7]), ("3d3D", [18, 19, 20])]
get(ENV, "WITH5D", "false") == "true"  &&  push!(terms, ("3d5D", [14, 15, 16, 17]))
for (name, fi) in vcat([("3P", initial)], terms)
    println(rpad(name, 6), join([@sprintf("%d:J=%s%s %.4f eV", i, string(mp.levels[i].J), string(mp.levels[i].parity), eV(i)) for i in fi], "   "))
end
pairs = haskey(ENV, "LINES") ? [Tuple(parse.(Int, split(x, ":"))) for x in split(ENV["LINES"], ",")] :
                               [(i, f) for i in initial for (_, fi) in terms for f in fi]
w(l)  = 2 * Float64(l.initialLevel.J) + 1
sw    = sum(2 * Float64(mp.levels[i].J) + 1 for i in initial)
integral(l) = (th = deg2rad.(l.angles);   f = l.dcs .* sin.(th);   2pi * sum((th[i+1]-th[i]) * (f[i]+f[i+1]) / 2 for i in 1:length(th)-1))
const A0SQ = 0.28002852   # a0^2 in 1e-16 cm^2

ics = open(OUT * "-ics.csv", "w");     println(ics, "model,energy,", join([t[1] for t in terms], ","), "   # 1e-18 cm^2")
fsr = open(OUT * "-fs.csv", "w");      println(fsr, "model,energy,Ji_level,Jf_level,Ji,Jf,sigma_1e-18cm2,intdcs_over_sigma")
# All energies go into ONE computeLines call per exchange model: the spin-angular tables and the non-local exchange matrices depend
# only on the bound levels, and their caches live for one call, so they are then built once instead of once per energy; the
# (energy, line) pairs also give the threads enough work to stay busy. The radial integrals are cached per line and energy anyway,
# since the continuum orbitals change with the energy. The price is that nothing is written until the whole model has finished.
# The DCS costs the angular sum at every partial wave, so it is computed only at the energies in DCSE: one call with and one without.
for model in MODELS
    results = ImpactExcitation.Line[]
    # The energy shift is one number per call, so the lines of each term go into their own call when SHIFT is given.
    groups = isempty(SHIFT) ? [(pairs, 0.0)] :
             [([p for p in pairs if p[2] in fi], get(SHIFT, name, 0.0)) for (name, fi) in terms if any(p -> p[2] in fi, pairs)]
    for (energies, dcs) in ((filter(E -> E in DCSE, ENERGIES), true), (filter(E -> !(E in DCSE), ENERGIES), false)), (gpairs, shift) in groups
        isempty(energies)  &&  continue
        set = ImpactExcitation.Settings(ImpactExcitation.Settings(); electronEnergies = energies, calcRateCoefficient = false, maxKappa = MAXKAPPA,
                                        coulombBetheKappa = CBKAPPA, asymptoticTail = true, printBefore = false, continuumExchange = model,
                                        orthogonalizeContinuum = model == :nonlocal, calcDCS = dcs, angles = dcs ? ANGLES : Float64[],
                                        topUp = TOPUP, energyShift = shift, lineSelection = LineSelection(true, indexPairs = gpairs))
        st    = redirect_stdout(devnull) do;   @timed ImpactExcitation.computeLines(mp, mp, Nuclear.Model(8.), GRID, set, output=true);   end
        append!(results, st.value)
        @printf("%-16s %d energies (DCS %s, shift %.3f eV) in %.1f s\n", model, length(energies), dcs, shift, st.time);   flush(stdout)
    end
    Ea    = [Defaults.convertUnits("energy: to atomic", E) for E in ENERGIES]
    for (k, E) in enumerate(ENERGIES)
        lines = [l for l in results if argmin(abs.(Ea .- l.initialElectronEnergy)) == k]
        for l in sort(lines, by = l -> (l.initialLevel.index, l.finalLevel.index))
            println(fsr, join([model, E, l.initialLevel.index, l.finalLevel.index, l.initialLevel.J, l.finalLevel.J,
                               100 * A0SQ * l.crossSection, isempty(l.dcs) ? NaN : integral(l) / l.crossSection], ","))
        end
        sig = [sum(w(l) * l.crossSection for l in lines if l.finalLevel.index in fi; init = 0.0) / sw for (_, fi) in terms]
        println(ics, model, ",", E, ",", join(100 * A0SQ .* sig, ","));   flush(ics);   flush(fsr)
        E in DCSE  &&  open(@sprintf("%s-dcs-%s-%geV.csv", OUT, model, E), "w") do io
            println(io, "theta,", join([t[1] for t in terms], ","), "   # 1e-18 cm^2/sr")
            cols = [sum((w(l) / sw) .* l.dcs for l in lines if l.finalLevel.index in fi) for (_, fi) in terms]
            for (j, t) in enumerate(ANGLES)   println(io, t, ",", join([100 * A0SQ * c[j] for c in cols], ","))   end
        end
        @printf("%-16s %6.1f eV  %2d lines   %s  [1e-18 cm^2]\n", model, E, length(lines),
                join([@sprintf("%s %8.3f", t[1], 100 * A0SQ * s) for (t, s) in zip(terms, sig)], "   "))
    end
    flush(stdout)
end
close(ics);   close(fsr)
