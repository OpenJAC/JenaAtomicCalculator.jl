#
println("Ic) Tests of the StrongField module to demonstrate the Coulomb asymmetry in ATI azimuthal angular distributions.")

#Laser pulse
wavelength      = 800.                             # nm
intensity       = 2e14                             # W/cm^2
CEP             = 0.                               # carrier-envelope phase 0 <= CEP < 2*pi
envelope        = Pulse.SinSquaredEnvelope(4)      # InfiniteEnvelope(), SinSquaredEnvelope(np), GaussianEnvelope(np) with number of optical cycles np
polarization    = Basics.RightElliptical(0.8)	   # RightCircular(), LeftCircular(); for SinSquared also: RightElliptical(epsilon), LeftElliptical(epsilon)

#Derived quantities
omega           = convertUnits("energy: from wavelength [nm] to atomic", wavelength)
intensity       = convertUnits("intensity: from W/cm^2 to atomic", intensity)
A0              = Pulse.computeFieldAmplitude(intensity, omega)
beam            = Pulse.PlaneWaveBeam(A0, omega, CEP)  

#Strong-field observable
observable      = StrongField.SfaAzimuthalAngularDistribution(pi/2, 100, 2.21*omega)  # Arguments: theta, number of azimuthal angles, energy

if true
    # Last successful:  unknown ...
    # Prepare the atomic target
	asfSettings     = AsfSettings(AsfSettings(), generateScf=true)
	rGrid           = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 2.0e-2, rbox = 20.0)

	initialName     = "Li ground-state"
	nuclearModel    = Nuclear.Model(3.39336) # Nuclear.Model(Z): Z effective nuclear charge to match the ionization potential
	refConfigInitial= [Configuration("[He] 2s")]

	finalName       = "He-like ground-state"
	refConfigFinal  = [Configuration("[He]")]

	AtomicComp      = Atomic.Computation(Atomic.Computation(), name=initialName, grid=rGrid, nuclearModel=nuclearModel, 
	                                     configs=refConfigInitial,  asfSettings=asfSettings )
	AtomicData      = perform(AtomicComp, output=true)
	initialLevel    = AtomicData["multiplet:"].levels[1]

	AtomicComp      = Atomic.Computation(Atomic.Computation(), name=finalName, grid=rGrid, nuclearModel=nuclearModel, 
	                                     configs=refConfigFinal,  asfSettings=asfSettings )
	AtomicData      = perform(AtomicComp, output=true)
	finalLevel      = AtomicData["multiplet:"].levels[1]

elseif true
    # Last successful:  unknown ...
    # Perform the StrongField computation; use a hydrogen-like 1s initial state and FreeVolkov continuum
	hydrogenic 	    = true
	hydrogenic1s    = true
	volkov          = StrongField.FreeVolkov()

	sfaSettings     = StrongField.Settings([E1], "VelocityGauge", true, true, hydrogenic, hydrogenic1s, true)
	SFIComp	        = StrongField.Computation(observable, nuclearModel, rGrid, initialLevel, finalLevel, beam, envelope, polarization, volkov, sfaSettings)
	SFIData         = StrongField.perform(SFIComp, output=true)
	
	#Atomic initial state and CoulombVolkov continuum
	hydrogenic 	    = false
	hydrogenic1s    = false
	volkov          = StrongField.CoulombVolkov(1.0) 

	sfaSettings     = StrongField.Settings([E1], "VelocityGauge", true, true, hydrogenic, hydrogenic1s, true)
	SFIComp2	    = StrongField.Computation(observable, nuclearModel, rGrid, initialLevel, finalLevel, beam, envelope, polarization, volkov, sfaSettings)
	SFIData2        = StrongField.perform(SFIComp2, output=true)

elseif true
    # Last successful:  unknown ...
    # Export data and plots
    StrongField.exportData([SFIComp,SFIComp2], [SFIData,SFIData2], "example-Ic")
    StrongField.plot([SFIComp,SFIComp2], [SFIData,SFIData2], "omega", "linear", "example-Ic")
end
