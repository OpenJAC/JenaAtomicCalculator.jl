
#
println("Ff) Computations for RR plasma rate coefficients with module Cascade (atoms, 2022).")

using JLD2
#
setDefaults("unit: energy", "Hartree")
setDefaults("method: continuum, Galerkin")            ## setDefaults("method: continuum, Galerkin") setDefaults("method: continuum, asymptotic Coulomb") 
setDefaults("method: normalization, pure sine")       ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")
grid          = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 0.6e-2, rbox = 10.0)
temps         = [1.0e+3,   1.0e+4,   1.0e+5,   1.0e+6,   1.0e+7,  1.0e+8]
# temps_au    = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", 1.0) .* temps;     @show temps_au


if  true
    # Last successful:  unknown ...
    # RR plasma rate coefficients for initially Li-like ion: Radiative recombination cascade computations
    setDefaults("print summary: open", "zzz-RadiativeRecombinationCascade.sum")
    name                 = "Radiative recombination cascade for Li-like ions"
    intoShells           = Basics.generateShellList(2,  8, "f")
    excitationScheme     = Cascade.RadiativeRecombinationScheme([E1], [0, 1], 10, 1.0, -0.0, 0., intoShells)
    wa   = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(26.), 
                               grid=grid, approach=Cascade.AverageSCA(), scheme=excitationScheme,
                               initialConfigs=[Configuration("1s^2")] )
    perform(wa; output=true)
    setDefaults("print summary: close", "")
    #
elseif true
    # Last successful:  unknown ...
    # RR plasma rate coefficients for initially Li-like ion: Radiative recombination cascade simulations
    setDefaults("print summary: open", "zzz-RadiativeRecombinationSimulationsum")
    simulationSettings = Cascade.SimulationSettings(true, false, 0.)

    data = [JLD2.load("zzz-cascade-rr-rate-computations-2022-01-08T20.jld")]
    name = "Radiative recombination simulations for Li-like ions"
    prop = Cascade.RrRateCoefficients(1, temps)

    wb   = Cascade.Simulation(Cascade.Simulation(), name=name, property=prop, settings=simulationSettings, computationData=data )
    perform(wb; output=true)
    setDefaults("print summary: close", "")
    #
elseif  true
    # Last successful:  unknown ...
    # Compute 
    function  alphaRR(T, p)
        s0 = sqrt(T/p[3]);   s1 = sqrt(T/p[4])
        wa = p[1] / (s0 * (1+s0)^(1-p[2]) * (1+s1)^(1+p[2]))
    end
    setDefaults("unit: energy", "eV")
    setDefaults("unit: rate", "1/s");   setDefaults("unit: strength", "cm^2 eV")
    using Plots; pyplot()
    # Figure 5: Plasma rate coefficients as function of temperature
    ## xList        =    [1.0e+3,   3.0e+3,  1.0e+4,   3.0e+4,   1.0e+5,   3.0e+5,   1.0e+6,   3.0e+6,  1.0e+7,  3.0e+7,  1.0e+8]
    xList        =    [10.0^T  for T =  3.5:0.25:4.0 ]
    # Verner etal (1996) values, O VII
    abT0T1xx     =    (4.897e-10,  0.7048,  1.906e+2,  4.093e+7)
    alphaRRxx    =    [ alphaRR(T, abT0T1xx)  for  T in xList ]
    # Verner etal (1996) values, C IV
    abT0T1yy     =    (8.540e-11,  0.5247,  5.014e+2,  1.479e+7)
    alphaRRyy    =    [ alphaRR(T, abT0T1yy)  for  T in xList ]
    # Verner etal (1996) values, C IV
    abT0T1Fe24p  =    (1.198e-9,   0.6443,  3.789e+3,  1.437e+8)
    alphaRRFe24p =    [ alphaRR(T, abT0T1Fe24p)  for  T in xList ]
    
    plot()
    ## plot!(xList, alphaRRxx, smooth=false,      lw=3.0, seriestype =:path, label="O VII")
    ## plot!(xList, alphaRRyy, smooth=false,      lw=3.0, seriestype =:path, label="C IV")
    plot!(xList, alphaRRFe24p, smooth=false,   lw=3.0, seriestype =:path, label="Fe^24+")
    plot!(xlabel="Plasma temperature [K]", ylabel="Rate coefficients [cm^3/s]", gridlw=0.1, 
          tickfontsize=12, labelfontsize=15, xscale=:log10, yscale=:log10)
    #
end
