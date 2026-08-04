
println("Ea)  Test of the BeamPhotoExcitation module with ASF from an internally generated initial- and final-state multiplet.")

setDefaults("print summary: open", "zzz-BeamPhotoExcitation.sum")


if  true
    # Last successful:  unknown ...
    # Test of beam-photo-excitation computatations of dominant-multipole regions and intensity pattern
    ## grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 20.0)
    grid = Radial.Grid(true)
    ## Basics.Cartesian2DMesh, Basics.GLegenreMesh, Basics.LinearMesh, Basics.PolarMesh, Basics.SphercialMesh
    obsMesh    = Basics.Cartesian2DMesh(Basics.GLegenreMesh(0., 15., 20), Basics.GLegenreMesh())            
    beamType   = Beam.LaguerreGauss(2, 0, 0.4)      ## Beam.PlaneWave, Beam.BesselBeam, Beam.LaguerreGauss, Beam.Component
    observable = Beam.DominantMultipoles(obsMesh)   ## Beam.DominantMultipoles,  Beam.IntensityPattern, Beam.AnisotropyParameter
    bSettings  = BeamPhotoExcitation.Settings(beamType, observable, [E1, M1], [UseCoulomb, UseBabushkin], true, 
                                              LineSelection(true, [(1,1), (1,2)], Tuple{LevelSymmetry,LevelSymmetry}[]), 0. )
    wa = Atomic.Computation(Atomic.Computation(), name="BeamPhotoExcitation computations", grid=grid, nuclearModel=Nuclear.Model(6.), 
                            initialConfigs = [Configuration("1s^2 2s")],
                            finalConfigs   = [Configuration("1s^2 2p")], 
                            processSettings= bSettings )

    wb = perform(wa)
    #
end
#
setDefaults("print summary: close", "")

#==
Further examples:
-----------------
++ Dominant multipole regions (M1, E2) for 
   Be  2s^2  ^1S_0  -->  2s 2p  ^3P_1,2
   Mg  3s^2  ^1S_0  -->  3s 3p  ^3P_1,2
   
++ Dominant multipole regions (E1, M2) for 
   Na  3s  -->  3p
   compare polarized, partly polarized and unpolarized sodium atoms
   LG (lambda=1) and LG  (lambda=-1) for the same set of polarizations
   LG-based radial vector beams == Beam.SuperposedBeam

++ Intensity pattern (only E2)
   Na  3s  -->  3d_3/2, 3d_5/2
   compare three different beam waists
   compare full intensity pattern vs. dominant multipole (M1, E2) multipole regions.
   
++ Anisotropy parameter 
   Na  3s  -->  3p, analog to your HG study but by considering different superposition
   of beams.
   
==#
