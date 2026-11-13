
##  Autoionization
```@autodocs
Modules = [AutoIonization]
Order   = [:type, :function]
```

##  Crystal-field-resolved emission
```@autodocs
Modules = [CrystalFieldEmission]
Order   = [:type, :function]
```

##  Decay yields
```@autodocs
Modules = [DecayYield]
Order   = [:type, :function]
```

##  Dielectronic Recombination
```@autodocs
Modules = [DielectronicRecombination]
Order   = [:type, :function]
```

##  Hyperfine-induced transitions
```@autodocs
Modules = [HyperfineInduced]
Order   = [:type, :function]
```

##  Electron Impact Excitation
```@autodocs
Modules = [ImpactExcitation]
Order   = [:type, :function]
```

### Non-Maxwellian excitation rates

Energy-resolved `ImpactExcitation.Line` results can be averaged over isotropic
κ or Druyvesteyn distributions without recalculating their scattering amplitudes:

```julia
using JenaAtomicCalculator: ImpactExcitation

# lines contains computed collision strengths for each transition at several energies.
settings = ImpactExcitation.Settings(ImpactExcitation.Settings(); temperatures=[1.e4, 5.e4, 1.e5])
kappaRates = ImpactExcitation.computeRateCoefficients(lines, settings; distribution=:Kappa, kappa=3.0)
druyRates = ImpactExcitation.computeRateCoefficients(lines, settings; distribution=:Druyvesteyn)
ImpactExcitation.displayResults(kappaRates)
```

Both models use temperatures in Kelvin with mean electron energy `3kBT/2`;
κ must be finite and greater than `3/2`. The public
`electronEnergyDistribution(energy, temperature; distribution=:Kappa, kappa=3.0)`
returns a normalized probability density in Hartree⁻¹ for incoming energy in Hartree.
For κ, the energy scale is `(κ-3/2)kBT`; for Druyvesteyn it is
`(3/2)kBT Γ(3/4)/Γ(5/4)` in `sqrt(E) exp(-(E/b)^2)`.

The non-Maxwellian result is `DistributionRateCoefficients`, containing rates in
cm³/s, the model/κ and the lower/upper extrapolated fractions of each rate.
It does not report Maxwellian effective collision strengths. Omitting the
distribution keyword preserves the existing Maxwellian result type and algorithm;
the existing `Settings` layout and automatic `computeLines` rate calculation are
unchanged. The new models are selected explicitly in the postprocessing call above.

Integration uses the actual shifted threshold and the same declared interpolation
and extrapolation of collision strengths. A warning identifies rates with over 1%
contribution from outside the sampled energy range. This fraction is a diagnostic,
not an error estimate: especially for small κ, extend and converge the high-energy
table rather than relying on the assumed tail. Rates use the conventional
nonrelativistic velocity/collision-strength relation, not a relativistic thermal model.

Distribution references: [Dzifčáková et al., KAPPA](https://arxiv.org/abs/1502.00853)
and [Druyvesteyn/mean-energy conventions](https://www.comsol.com/blogs/electron-energy-distribution-function).

##  Photoemission
```@autodocs
Modules = [PhotoEmission]
Order   = [:type, :function]
```

##  Photoexcitation
```@autodocs
Modules = [PhotoExcitation]
Order   = [:type, :function]
```

##  Photoionization
```@autodocs
Modules = [PhotoIonization]
Order   = [:type, :function]
```

##  Photorecombination
```@autodocs
Modules = [PhotoRecombination]
Order   = [:type, :function]
```

##  Two-electron one-photon transitions
```@autodocs
Modules = [TwoElectronOnePhoton]
Order   = [:type, :function]
```

##  Multi-photon transitions
```@autodocs
Modules = [MultiPhotonTransition]
Order   = [:type, :function]
```


##  Particle scattering
```@autodocs
Modules = [ParticleScattering]
Order   = [:type, :function]
```


##  Rayleigh-Compton scattering
```@autodocs
Modules = [RayleighCompton]
Order   = [:type, :function]
```

##  Electron-impact ionization
```@autodocs
Modules = [ImpactIonization]
Order   = [:type, :function]
```


##  Coulomb ionization by ion impact
```@autodocs
Modules = [CoulombIonization]
Order   = [:type, :function]
```

##  Multi-photon ionization
```@autodocs
Modules = [MultiPhotonIonization]
Order   = [:type, :function]
```

##  Photorecombination with RR-DR interference
```@autodocs
Modules = [PhotoRecombinationInterference]
Order   = [:type, :function]
```

##  Resonant electron-impact ionization
```@autodocs
Modules = [ResonantImpactIonization]
Order   = [:type, :function]
```

##  Generalized oscillator strengths
```@autodocs ; canonical=false
Modules = [GeneralizedOscillatorStrength]
Order   = [:type, :function]
```
