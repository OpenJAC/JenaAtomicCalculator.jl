
println("Ca) Apply & test the Einstein module with ASF from an internally generated multiplet.")

setDefaults("print summary: open", "zzz-Einstein.sum")
setDefaults("unit: energy", "Kayser")
grid=Radial.Grid(true)
## grid = Radial.Grid(Radial.Grid(true), rnt = 2.0e-5,h = 1.0e-2, hp = 1.0e-2, NoPoints = 2000)


if  true
    # Last successful:  30-Aug-2026 -- Fe (Z=26), three configurations, 16 E1 transitions, and the
    #   two gauges agree reasonably: Coulomb/Babushkin has median 1.039 over all 16, range 0.544 to 1.083.
    #   The single 0.544 sits on a weak line, where the two gauges are differences of large cancelling terms
    #   and disagreement is expected. NOT DATED only because Rule 7 leaves that to a human; this branch is a
    #   Dated on the maintainer's decision of 30-Aug-2026.
    # Compute 
    asfSettings   = AsfSettings(AsfSettings(), scField=Basics.DFSField())
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(26.), asfSettings=asfSettings,
                            configs=[Configuration("[Ne] 3s 3p^6"), Configuration("[Ne] 3s^2 3p^4 4s"), Configuration("[Ne] 3s^2 3p^5")],
                            propertySettings=[ Einstein.Settings([E1], true, LineSelection(false), 0., 0., 10000. ) ])

    wb = perform(wa; output=true)
    #
elseif  true
    # Last visit:  30-Aug-2026 -- RUNS (exit 0), and WHAT STOPS A DATE IS THE GAUGES. C (Z=6), four
    #   configurations, 13 E1 transitions. The strong lines to levels 6 and 7 (A ~ 4e8 1/s) agree to 25 %,
    #   which is ordinary; but the equally strong lines 13-->11 and 13-->12 (A ~ 1e8 and 7e7 1/s) come out
    #   with Coulomb/Babushkin = 6.28 and 6.31, and 15-->11, 15-->12 with 4.01 and 4.12. A factor of four to
    #   six between the gauges on transitions that are NOT weak is not something to date. (The 2.9e4 ratio on
    #   13-->3 is a different matter and is fine: A = 3e2 against 1e-2, i.e. a near-cancelling weak line.)
    #   Whether this is correlation, the basis, or something in the module is NOT established here.
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
    # Last visit:  30-Aug-2026 -- RUNS (exit 0) AND COMPUTES NOTHING, which is why it cannot be dated.
    #   The configuration list is the single closed shell 1s^2 2s^2 at Z=36, so the multiplet holds ONE level
    #   (J=0, even) and there is no pair of levels to connect: the Einstein table comes out empty. The branch
    #   is therefore correct and useless as it stands -- it would need a second configuration to test
    #   anything. A does-it-run sweep would report this branch as healthy.
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
    # Last successful:  30-Aug-2026 -- C (Z=6), 1s^2 2s^2 + 1s^2 2s 3p, two E1 transitions, and the
    #   gauge ratio is 1.1898 and 1.1919 -- consistent between the two lines, so the 19 % is systematic rather
    #   than erratic, as a length-against-velocity difference in a small CI usually is. NOT DATED only because
    #   Rule 7 leaves that to a human; a date here would record a 19 % gauge spread, which is worth saying
    #   out loud rather than leaving implicit. Dated on the maintainer's decision of 30-Aug-2026.
    # Compute 
    ## grid = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 1.0e-2, hp = 0., rbox = 2.0)
    ## asfSettings   = AsfSettings(AsfSettings(), frozenSubshells=[Subshell("1s_1/2"), Subshell("2s_1/2"), Subshell("2p_1/2"), Subshell("2p_3/2"), Subshell("3d_3/2"), Subshell("3d_5/2")],
    ## asfSettings   = AsfSettings(AsfSettings(), frozenSubshells=[Subshell("1s_1/2"), Subshell("2s_1/2")],
    ##                                            startScfFrom=StartFromPrevious(frozenOrbs), generateScf=true)
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(6.), ## asfSettings=asfSettings,
                            # `properties=[JAC.EinsteinX()]` with a separate `einsteinSettings=` stood here until
                            # 29-Aug-2026. Atomic.Computation carries ONE list now -- propertySettings -- and the
                            # settings object itself says which property is wanted, so the type is no longer named.
                            configs=[Configuration("1s^2 2s^2"), Configuration("1s^2 2s 3p")],
                            propertySettings=[Einstein.Settings([E1], true, LineSelection(), 0., 0., 10000. )] )
                            ## einsteinSettings=Einstein.Settings([E1,M1,E2], true, LineSelection(true, indexPairs=[(3,1), (5,1)]), 0., 0., 10000. ) )

    wb = perform(wa; output=true)
    #
end
#
setDefaults("print summary: close", "")
