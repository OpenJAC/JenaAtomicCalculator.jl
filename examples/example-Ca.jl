
println("Ca) Apply & test the Einstein module with ASF from an internally generated multiplet.")

setDefaults("print summary: open", "zzz-Einstein.sum")
setDefaults("unit: energy", "Kayser")
grid=Radial.Grid(true)
## grid = Radial.Grid(Radial.Grid(true), rnt = 2.0e-5,h = 1.0e-2, hp = 1.0e-2, NoPoints = 2000)


if  true
    # Last successful:  unknown ...
    # Compute 
    asfSettings   = AsfSettings(AsfSettings(), scField=Basics.DFSField())
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(26.), asfSettings=asfSettings,
                            configs=[Configuration("[Ne] 3s 3p^6"), Configuration("[Ne] 3s^2 3p^4 4s"), Configuration("[Ne] 3s^2 3p^5")],
                            propertySettings=[ Einstein.Settings([E1], true, LineSelection(false), 0., 0., 10000. ) ])

    wb = perform(wa; output=true)
    #
elseif  true
    # Last successful:  unknown ...
    # Compute 
    ## grid = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 1.0e-2, hp = 0., rbox = 2.0)
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(6.), 
                            configs=[Configuration("1s^2 2s"), Configuration("1s 2s^2"), Configuration("1s 2s 2p"), Configuration("1s 2p^2")],
                            propertySettings=[ Einstein.Settings([E1], true, LineSelection(true, indexPairs=[(13,0), (15,0)]), 0., 0., 10000. )] )
                            ## Einstein.Settings([M2], true, LineSelection(true, indexPairs=[(13,2), (15,5), (15,4)]), 0., 0., 10000. ) )

    wb = perform(wa; output=true)
    #
elseif  false
    # Last successful:  unknown ...
    # Compute 
    ## grid = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 1.0e-2, hp = 0., rbox = 2.0)
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(36.), 
                            configs=[Configuration("1s^2 2s^2")],
                            ## configs=[Configuration("1s^2 2s"), Configuration("1s^2 2p"), Configuration("1s^2 3d"), Configuration("1s^2 4f")],
                            ##x configs=[Configuration("1s 2s^2"), Configuration("1s 2s 2p"), Configuration("1s^2 2s")],
                            propertySettings=[ Einstein.Settings([E1, M1, E2], true, LineSelection(), 0., 0., 10000. ) ])
                            ##                 Einstein.Settings([M2], true, LineSelection(true, indexPairs=[(13,2), (15,5), (15,4)]), 0., 0., 10000. ) )

    wb         = perform(wa; output=true)
    frozenOrbs = wb["multiplet:"].levels[1].basis.orbitals
    #
elseif  true
    # Last successful:  unknown ...
    # Compute 
    ## grid = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 1.0e-2, hp = 0., rbox = 2.0)
    ## asfSettings   = AsfSettings(AsfSettings(), frozenSubshells=[Subshell("1s_1/2"), Subshell("2s_1/2"), Subshell("2p_1/2"), Subshell("2p_3/2"), Subshell("3d_3/2"), Subshell("3d_5/2")],
    ## asfSettings   = AsfSettings(AsfSettings(), frozenSubshells=[Subshell("1s_1/2"), Subshell("2s_1/2")],
    ##                                            startScfFrom=StartFromPrevious(frozenOrbs), generateScf=true)
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(6.), ## asfSettings=asfSettings,
                            properties=[JAC.EinsteinX()], 
                            ##x configs=[Configuration("1s^2 2s"), Configuration("1s^2 2p"), Configuration("1s^2 3d"), Configuration("1s^2 4f")],
                            configs=[Configuration("1s^2 2s^2"), Configuration("1s^2 2s 3p")],
                            einsteinSettings=Einstein.Settings([E1], true, LineSelection(), 0., 0., 10000. ) )
                            ## einsteinSettings=Einstein.Settings([E1,M1,E2], true, LineSelection(true, indexPairs=[(3,1), (5,1)]), 0., 0., 10000. ) )

    wb = perform(wa; output=true)
    #
end
#
setDefaults("print summary: close", "")
