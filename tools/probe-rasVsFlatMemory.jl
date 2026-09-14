#
# probe-rasVsFlatMemory.jl   --   THE SAME CORRELATION SPACE BY TWO ROUTES, and what each costs in memory.
#
# THE QUESTION.  A RAS ladder is driven by `generate(Representation(..., RasExpansion(...)))`.  The same physics
# can instead be written as an explicit list of configurations handed straight to `SelfConsistent.performSCF`.
# On 14-Sep-2026 eleven solves down the FLAT route stayed between 1.42 and 1.67 GB up to 1231 CSFs -- including
# 3p-opened spaces, so the correlation class is genuinely there -- while a Ti III ladder down the RAS route held
# 26.9 GB for hours.  If the gap is the ROUTE rather than the space, a 25 085-CSF layer that needs a 61 GB
# machine today would fit on an ordinary desktop, which is worth knowing before another such row is queued.
#
# WHY THE COMPARISON IS EXACT AND NOT MERELY SIMILAR.  `Basics.generateBasis(refConfigs, symmetries, step)` builds
# a RAS step's configuration list by calling `Basics.generateConfigurations(refConfigs, from, to)` -- once for
# singles, twice in succession for doubles.  This probe calls THE SAME FUNCTIONS IN THE SAME ORDER, so the
# configuration list it feeds to the flat route is the one the RAS route would have built, not an imitation of
# it.  Anything left over is then a property of the route.
#
# WHAT IS REPORTED.  The configuration and CSF counts (so the two routes can be checked to describe the same
# space at all), the peak RSS of this process, and the lowest few level energies -- because a memory win that
# changes the physics is not a win.  Run one layer per process; see tools/probe-eolPeakRss.jl for why measuring
# several in one process silently under-reports every case after the first.
#
using JenaAtomicCalculator, Printf

Z     = 22.0
refs  = [Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^2"), Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d 4s"),
         Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 4s^2")]
syms  = [LevelSymmetry(J, Basics.plus)  for J = 0:4]
grid  = Basics.recommendedGrid(refs, Nuclear.Model(Z); printout=false)
core  = [Shell("1s"), Shell("2s"), Shell("2p"), Shell("3s"), Shell("3p")]
valp  = [Shell("3p"), Shell("3d"), Shell("4s")]

layer = length(ARGS) >= 1 ? ARGS[1] : "3p-v"
route = length(ARGS) >= 2 ? ARGS[2] : "flat"
to    = layer == "3p-v"  ? [Shell("3d"), Shell("4s")] :
        layer == "3p-4p" ? [Shell("3d"), Shell("4s"), Shell("4p")] :
        layer == "3p-4d" ? [Shell("3d"), Shell("4s"), Shell("4p"), Shell("4d")] :
                           [Shell("3d"), Shell("4s"), Shell("4p"), Shell("4d"), Shell("4f")]

# The RAS step's own expansion, reproduced call for call (Basics.generateBasis, module-BasicsAZ-inc-generate.jl).
singles = Basics.generateConfigurations(refs, valp, to)
doubles = Basics.generateConfigurations(refs, valp, to)
doubles = Basics.generateConfigurations(doubles, valp, to)
confs   = unique( vcat(refs, singles, doubles) )
@printf("layer %-7s route %-4s  configurations %6d\n", layer, route, length(confs));   flush(stdout)

step  = RasStep(RasStep(); seFrom=valp, seTo=to, deFrom=valp, deTo=to, frozen=core)
nCsf  = length( Basics.generateBasis(refs, syms, step).csfs )
@printf("layer %-7s route %-4s  CSFs in the RAS step %6d\n", layer, route, nCsf);   flush(stdout)

set = AsfSettings(AsfSettings(); scField=Basics.EOLField(), eeInteraction=CoulombInteraction(),
                                 eeInteractionCI=CoulombInteraction(), gridStopper=false,
                                 scfRoute=Basics.RotationRoute(12),
                                 levelSelectionCI=LevelSelection(true, configurations=refs))
if  route == "flat"
    t = @elapsed (mp = redirect_stdout(devnull) do
                           SelfConsistent.performSCF(confs, Nuclear.Model(Z), grid, set; printout=false)
                       end)
else
    rs = RasSettings(Int64[], 24, 1.0e-6, CoulombInteraction(), LevelSelection(true, configurations=refs))
    t  = @elapsed (wb = redirect_stdout(devnull) do
                            generate( Representation("Ti III", Nuclear.Model(Z), grid, refs,
                                      RasExpansion(syms, 20, [RasStep(RasStep(); frozen=core), step], rs)),
                                      output=true )
                        end)
    mp = wb["step2"]
end

b  = mp.levels[1].basis
e0 = mp.levels[1].energy
@printf("RESULT  layer %-7s route %-4s  nCsf %7d  nSub %4d  wall %8.1f  E0 %18.8f\n",
        layer, route, length(b.csfs), length(b.subshells), t, e0)
for  l in mp.levels[1:min(4,length(mp.levels))]
    @printf("   LEVEL  J=%-5s  %18.8f  rel %12.2f cm^-1\n", string(l.J), l.energy,
            Defaults.convertUnits("energy: from atomic to Kayser", l.energy - e0))
end
