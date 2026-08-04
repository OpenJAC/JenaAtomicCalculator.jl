#
println("Id) Tests of the StrongField2 module to compute in ATI photoelectron momentum distributions.")

#Laser pulse
wavelength      = 800.                             # nm
intensity       = 2e14                             # W/cm^2
CEP             = pi/4                             # carrier-envelope phase 0 <= CEP < 2*pi
envelope        = Pulse.SinSquaredEnvelope(2)      # InfiniteEnvelope(), SinSquaredEnvelope(np), GaussianEnvelope(np) with number of optical cycles np
polarization    = Basics.RightCircular()	       # RightCircular(), LeftCircular(); for SinSquared also: RightElliptical(epsilon), LeftElliptical(epsilon)

#Electronic states                                                                  
hydrogenic      = false  #true: Use hydrogenic wave functions for the initial state (quantum numbers n,l,m and ionization potential are taken from initialLevel)
hydrogenic1s    = false  #true: If both are true, a hydrogen-like 1s initial state is used with modified ionization potential
volkov          = StrongField2.DistortedVolkov() #FreeVolkov(), CoulombVolkov(Z) with charge Z of ion, or DistortedVolkov()
mAverage        = true   #true: Average initialState projections mj or ml and sum final state spin projection msp

#Derived quantities
omega           = convertUnits("energy: from wavelength [nm] to atomic", wavelength)
intensity       = convertUnits("intensity: from W/cm^2 to atomic", intensity)
A0              = Pulse.computeFieldAmplitude(intensity, omega)
beam            = Pulse.PlaneWaveBeam(A0, omega, CEP)  

#Strong-field observable
observable      = StrongField2.SfaMomentumDistribution(pi/2, 40, 20, 8*omega)  #Arguments: theta, No azimuthal angles, No energies, maximum energy

if true
    # Last successful:  unknown ...
    # Preprare the atomic target
	asfSettings      = AsfSettings(AsfSettings(), generateScf=true)
	rGrid            = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 2.0e-2, rbox = 20.0)

	initialName      = "Li ground-state"
	nuclearModel     = Nuclear.Model(3.39336) #Nuclear.Model(Z): Z effective nuclear charge to match the ionization potential
	refConfigInitial = [Configuration("[He] 2s")]

	finalName        = "He-like ground-state"
	refConfigFinal   = [Configuration("[He]")]

	AtomicComp       = Atomic.Computation(Atomic.Computation(), name=initialName, grid=rGrid, nuclearModel=nuclearModel, 
	                                      configs=refConfigInitial,  asfSettings=asfSettings )
	AtomicData       = perform(AtomicComp, output=true)
	initialLevel     = AtomicData["multiplet:"].levels[1]

	AtomicComp       = Atomic.Computation(Atomic.Computation(), name=finalName, grid=rGrid, nuclearModel=nuclearModel, 
	                                      configs=refConfigFinal,  asfSettings=asfSettings )
	AtomicData       = perform(AtomicComp, output=true)
	finalLevel       = AtomicData["multiplet:"].levels[1]

elseif true
    # Last successful:  unknown ...
    # Perform the StrongField2 computation
	sfaSettings = StrongField2.Settings([E1], "VelocityGauge", true, true, hydrogenic, hydrogenic1s, mAverage)
	SFIComp	    = StrongField2.Computation(observable, nuclearModel, rGrid, initialLevel, finalLevel, beam, envelope, polarization, volkov, sfaSettings)
	SFIData     = StrongField2.perform(SFIComp, output=true)
    StrongField2.exportData([SFIComp], [SFIData], "example-Id")
    #
end


