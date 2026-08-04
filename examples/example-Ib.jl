#
println("Ib) Tests of the StrongField module to calculate ATI energy distributions in the SFA with hydrogen-like initial states.")

#Laser pulse
wavelength      = 800.                             # nm
intensity       = 2e14                             # W/cm^2
CEP             = 0.                               # carrier-envelope phase 0 <= CEP < 2*pi
envelope        = Pulse.SinSquaredEnvelope(4)      # InfiniteEnvelope(), SinSquaredEnvelope(np), GaussianEnvelope(np) 
                                                   # with number of optical cycles np
polarization    = Basics.RightCircular()	       # RightCircular(), LeftCircular(); for SinSquared also: 
                                                   # RightElliptical(epsilon), LeftElliptical(epsilon)
              
#Electronic states                                                                  
hydrogenic      = true  #true: Use hydrogenic wave functions for the initial state (quantum numbers n,l,m and ionization potential are taken from initialLevel)
hydrogenic1s    = true  #true: If both are true, a hydrogen-like 1s initial state is used with modified ionization potential
volkov          = StrongField.FreeVolkov()  #F reeVolkov(), CoulombVolkov(Z) with charge Z of ion, or DistortedVolkov()
mAverage        = false #true: Average initialState projections mj or ml and sum final state spin projection msp

#Derived quantities
omega           = convertUnits("energy: from wavelength [nm] to atomic", wavelength)
intensity       = convertUnits("intensity: from W/cm^2 to atomic", intensity)
A0              = Pulse.computeFieldAmplitude(intensity, omega)
beam            = Pulse.PlaneWaveBeam(A0, omega, CEP)  

#Strong-field observable
observable      = StrongField.SfaEnergyDistribution(pi/2, pi/2, 100, 10*omega)   #Arguments: theta, phi, number of energies, maximum energy


if true
    # Last successful:  unknown ...
    # Prepare the atomic target
	asfSettings      = AsfSettings(AsfSettings(), generateScf=true)
	rGrid            = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 2.0e-2, rbox = 20.0)

	initialName      = "Li ground-state"
	nuclearModel     = Nuclear.Model(3.39336)  # Nuclear.Model(Z): Z effective nuclear charge to match the ionization potential
	refConfigInitial = [Configuration("[He] 2s")]

	finalName        = "He-like ground-state"
	refConfigFinal   = [Configuration("[He]")]

	AtomicComp       = Atomic.Computation(Atomic.Computation(), name=initialName, grid=rGrid, 
	                                     nuclearModel=nuclearModel, configs=refConfigInitial,  asfSettings=asfSettings )
	AtomicData       = perform(AtomicComp, output=true)
	initialLevel     = AtomicData["multiplet:"].levels[1]

	AtomicComp       = Atomic.Computation(Atomic.Computation(), name=finalName, grid=rGrid, 
	                                     nuclearModel=nuclearModel, configs=refConfigFinal,  asfSettings=asfSettings )
	AtomicData       = perform(AtomicComp, output=true)
	finalLevel       = AtomicData["multiplet:"].levels[1]

elseif true
    # Last successful:  unknown ...
    # Perform the StrongField computation
	sfaSettings      = StrongField.Settings([E1], "VelocityGauge", true, true, hydrogenic, hydrogenic1s, mAverage)
	SFIComp	         = StrongField.Computation(observable, nuclearModel, rGrid, initialLevel, finalLevel, 
	                                           beam, envelope, polarization, volkov, sfaSettings)
	SFIData          = StrongField.perform(SFIComp, output=true)

elseif true
    # Last successful:  unknown ...
    # Export the data and plots
    StrongField.exportData([SFIComp], [SFIData], "example-Ib")
    StrongField.plot([SFIComp], [SFIData], "omega", "linear", "example-Ib")
end
