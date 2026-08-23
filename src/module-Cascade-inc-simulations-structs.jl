

"""
`struct  Cascade.Data{T}`  ... defines a type for communicating different line lists to the cascade simulations

    + lines    ::Array{T,1}       ... List of lines to type T.
"""  
struct  Data{T}
    lines      ::Array{T,1}
end 


#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################



#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################


"""
`abstract type  Cascade.AbstractSimulationProperty`  
    ... defines an abstract and various singleton types for the different properties that can be obtained from the simulation of 
        cascade data.

    + struct DecayPathes                ... determine the major 'decay pathes' of the cascade.
    + struct DrRateCoefficients         ... simulate the DR (plasma) rate coefficients for given plasma temperatures. 
    + struct ElectronCoincidence        ... simulate electron-coincidence spectra (not yet).
    + struct FinalLevelDistribution     ... simulate the 'final-level distribution' as it is found after all cascade 
                                            processes are completed.
    + struct IonDistribution            ... simulate the 'ion distribution' as it is found after all cascade processes are completed.
    + struct MeanLineWidths             ... simulate the mean line widths of a line near to a given energy (not yet). 
    + struct RelaxationCurve            ... simulate the relaxation CURVE, i.e. the times by which 70%, 80%, 90%
                                            of the initially occupied levels have decayed to the ground
                                            configuration.        
    + struct ElectronIntensities        ... simulate the electron-line intensities as function of electron energy.
    + struct PhotoAbsorptionCS          ... simulate the (total) photoabsorption cross sections for a given set of photo-excitation 
                                            and ionization processes.
    + struct PhotoResonances            ... simulate the position and strength of photo-resonances for a given set of photo-excitation 
                                            amplitudes & cross sections (not yet).
    + struct PhotonIntensities          ... simulate the photon-line intensities as function of electron energy. 
    + struct PiRateCoefficients         ... simulate the PI (plasma) rate coefficients for given plasma temperatures (not yet). 
    + struct TimeBinnedPhotonIntensity  ... simulate the photon-line intensities as function of electron energy 
                                            in a given time interval (not yet).
"""
abstract type  AbstractSimulationProperty                              end


"""
`struct  Cascade.DrRateCoefficients   <:  Cascade.AbstractSimulationProperty`  
    ... defines a type for simulating the DR plasma rate coefficients as function of the (free) electron energy and
        plasma temperature.

    + initialLevelNo      ::Int64       ... Level No of initial level for which rate coefficients are to be computed.
    + electronEnergyShift ::Float64     
        ... (total) energy shifts that apply to all resonances when alpha^(DR) is computed.
    + temperatures        ::Array{Float64,1}
        ... temperatures [K] for which the DR plasma rate coefficieints to be calculated.
    + nDetailed           ::Int64       
        ... principal quantum of the `last' shell for which Auger and radiatiative amplitudes have been calculated
            by the Cascade.Computation; all contributions of this shell are scaled for nDetailed < n <= nMax
            also for higher shells by a simple scaling rule.
    + nMax                ::Int64 
        ... Maximum n (principal quantum number), for which contributions are scaled; NO scaling is taken
            into account for nMax <= nDetailed.
    + resonanceSelection  ::DielectronicRecombination.ResonanceSelection
"""  
struct  DrRateCoefficients   <:  Cascade.AbstractSimulationProperty
    initialLevelNo        ::Int64 
    electronEnergyShift   ::Float64 
    temperatures          ::Array{Float64,1}
    nDetailed             ::Int64 
    nMax                  ::Int64 
    resonanceSelection    ::DielectronicRecombination.ResonanceSelection
end 


"""
`Cascade.DrRateCoefficients()`  ... (simple) constructor for cascade DrRateCoefficients.
"""
function DrRateCoefficients()
    DrRateCoefficients(1, 0.,  Float64[], 0, 0, DielectronicRecombination.ResonanceSelection() )
end


# `Base.show(io::IO, dist::Cascade.DrRateCoefficients)`  ... prepares a proper printout of the variable data::Cascade.DrRateCoefficients.
function Base.show(io::IO, dist::Cascade.DrRateCoefficients) 
    println(io, "initialLevelNo:           $(dist.initialLevelNo)  ")
    println(io, "electronEnergyShift:      $(dist.electronEnergyShift)  ")
    println(io, "temperatures:             $(dist.temperatures)  ")
    println(io, "nDetailed:                $(dist.nDetailed)  ")
    println(io, "nMax:                     $(dist.nMax)  ")
    println(io, "resonanceSelection:       $(dist.resonanceSelection)  ")
end


"""
`struct  Cascade.RrRateCoefficients   <:  Cascade.AbstractSimulationProperty`  
    ... defines a type for simulating the RR plasma rate coefficients as function of the (free) electron energy and
        plasma temperature.

    + initialLevelNo      ::Int64       ... Level No of initial level for which rate coefficients are to be computed.
    + temperatures        ::Array{Float64,1}
        ... temperatures [K] for which the DR plasma rate coefficieints to be calculated.
    + multipoles          ::Array{EmMultipole}           
        ... Multipoles of the radiation field that are to be included into the excitation processes.
    + finalLevelSelection ::LevelSelection    
        ... Specifies the selected final levels of some given final-state configurations; these final level numbers/
            symmetries always refer to single configurations.
    + finalConfigurations ::Array{Configuration,1}
    

"""  
struct  RrRateCoefficients   <:  Cascade.AbstractSimulationProperty
    initialLevelNo        ::Int64 
    temperatures          ::Array{Float64,1}
    multipoles            ::Array{EmMultipole} 
    finalLevelSelection   ::LevelSelection
    finalConfigurations   ::Array{Configuration,1}
end 


"""
`Cascade.RrRateCoefficients()`  ... (simple) constructor for cascade RrRateCoefficients.
"""
function RrRateCoefficients()
    RrRateCoefficients(1, Float64[], EmMultipole[], LevelSelection(), Configuration[])
end


# `Base.show(io::IO, dist::Cascade.RrRateCoefficients)`  ... prepares a proper printout of the variable data::Cascade.RrRateCoefficients.
function Base.show(io::IO, dist::Cascade.RrRateCoefficients) 
    println(io, "initialLevelNo:           $(dist.initialLevelNo)  ")
    println(io, "temperatures:             $(dist.temperatures)  ")
    println(io, "multipoles:               $(dist.multipoles)  ")
    println(io, "finalLevelSelection:      $(dist.finalLevelSelection)  ")
    println(io, "finalConfigurations:      $(dist.finalConfigurations)  ")
end


"""
`struct  Cascade.EieRateCoefficients   <:  Cascade.AbstractSimulationProperty`
    ... defines a type for simulating electron-impact excitation (EIE) plasma rate coefficients and effective
        collision strengths from the ImpactExcitation.Line's of a previous cascade computation.

        This is the simulation counterpart of Cascade.ImpactExcitationScheme.  The computation produces LINES --
        one per (initial level, final level, incident electron energy), each carrying a cross section and a
        collision strength; forming the effective collision strength Upsilon(T) and the rate coefficient
        alpha(T) from them requires an interpolation over energy and an integration against a Maxwellian, and
        that is simulation work.  Compare Cascade.DrRateCoefficients and Cascade.RrRateCoefficients.

    + initialLevelNo      ::Int64             ... Level No of the initial level; 0 selects all initial levels.
    + temperatures        ::Array{Float64,1}  ... Temperatures [K] at which the rate coefficients are formed.
    + finalLevelSelection ::LevelSelection    ... Specifies the selected final levels, if any.
"""
struct  EieRateCoefficients   <:  Cascade.AbstractSimulationProperty
    initialLevelNo        ::Int64
    temperatures          ::Array{Float64,1}
    finalLevelSelection   ::LevelSelection
end


"""
`Cascade.EieRateCoefficients()`  ... (simple) constructor for cascade EieRateCoefficients.
"""
function EieRateCoefficients()
    EieRateCoefficients(1, Float64[], LevelSelection())
end


# `Base.show(io::IO, dist::Cascade.EieRateCoefficients)`  ... prepares a proper printout of the variable data::Cascade.EieRateCoefficients.
function Base.show(io::IO, dist::Cascade.EieRateCoefficients)
    println(io, "initialLevelNo:           $(dist.initialLevelNo)  ")
    println(io, "temperatures:             $(dist.temperatures)  ")
    println(io, "finalLevelSelection:      $(dist.finalLevelSelection)  ")
end


"""
`struct  Cascade.EaCrossSections   <:  Cascade.AbstractSimulationProperty`
    ... defines a type for simulating the excitation-autoionization (EA) contribution to the electron-impact
        ionization cross section from the lines of a previous Cascade.ElectronIonizationScheme computation.

        The EA cross section at a given impact energy is

            sigma^EA(E; i) = SUM_e  sigma^exc(E; i -> e) * B_a(e)

        i.e. the electron-impact excitation cross section into each autoionizing level e, weighted by the
        branching ratio B_a(e) = A_a / (A_a + A_r) with which that level decays by autoionization rather than
        radiatively.  The ElectronIonizationScheme computes no radiative rates, so A_r = 0 and B_a = 1 for every
        level that carries at least one Auger line; the sum then runs over exactly those levels, which is what
        the autoionization lines of the computation identify.  That is an upper bound on the EA cross section,
        and a good one for inner-shell excited levels of light ions, where autoionization dominates radiative
        decay by orders of magnitude -- but it IS an upper bound and should be read as such.

    + initialLevelNo      ::Int64             ... Level No of the initial level; 0 selects all initial levels.
    + electronEnergies    ::Array{Float64,1}  ... Impact energies at which sigma^EA is reported; these must be
                                                  among the energies of the underlying computation.
"""
struct  EaCrossSections   <:  Cascade.AbstractSimulationProperty
    initialLevelNo        ::Int64
    electronEnergies      ::Array{Float64,1}
end


"""
`Cascade.EaCrossSections()`  ... (simple) constructor for cascade EaCrossSections.
"""
function EaCrossSections()
    EaCrossSections(1, Float64[])
end


# `Base.show(io::IO, dist::Cascade.EaCrossSections)`  ... prepares a proper printout of the variable data::Cascade.EaCrossSections.
function Base.show(io::IO, dist::Cascade.EaCrossSections)
    println(io, "initialLevelNo:           $(dist.initialLevelNo)  ")
    println(io, "electronEnergies:         $(dist.electronEnergies)  ")
end


"""
`struct  Cascade.ResonantIonizationStrengths   <:  Cascade.AbstractSimulationProperty`
    ... defines a type for simulating the RESONANT contribution to electron-impact ionization from the lines of a
        previous Cascade.ElectronIonizationScheme computation that requested one or both of the resonant capture
        channels.

        An incident electron is CAPTURED into a doubly-excited resonance d, which then sheds two electrons, so that
        the ion ends one charge state HIGHER than it began.  What is reported is the energy-integrated RESONANCE
        STRENGTH of each resonance, the natural quantity for an isolated resonance, together with its sum:

            S_0(d)   = pi^2/k^2 * A(capture) * (2J_d+1)/(2J_i+1)              formation alone
            S_seq(d) = S_0 * SUM_n [A(d->n)/Gamma(d)] * [A_a(n)/Gamma(n)]     sequential, through each autoionizing n
            S_sim(d) = S_0 * [P_dbl * A_a(d)/Gamma(d)]                        simultaneous, one double-Auger step

        with Gamma the total width, Auger plus radiative.  The prefactor and the relativistic wave number are those
        of Cascade.simulateDrRateCoefficients, deliberately and identically, so that the RECOMBINATION and IONIZATION
        strengths of the same resonance may be placed side by side without a conversion between them -- which is the
        main thing this property is good for, since the two are competing decays of one resonance.

    APPROXIMATIONS, so that it is clear what these numbers are:

    + isolated resonances
        ... strengths are ADDED over resonances, with no interference and no overlap of widths.  Sound while the
            resonances stay narrow against their spacing; it degrades for high capture shells.
    + branchings from the rates that were generated
        ... every Gamma is a sum over the decay steps the cascade actually built.  A decay route absent from the
            configuration lists is absent from Gamma too, so the branchings still sum to one.  That checks the
            arithmetic, never the completeness.
    + an autoionizing intermediate is one that HAS an Auger line
        ... the sequential route needs the intermediate n to autoionize in its turn, and n is admitted exactly when
            the computation produced at least one Auger line out of it.  A level that would autoionize into a final
            configuration nobody generated therefore does not count, and its route is missing rather than wrong.
    + the simultaneous route is NOT computed here
        ... its double-Auger probability is supplied by the user through `dblAugerProbability`, because the choice of
            which subshells count as PASSIVE is a physics judgement the cascade cannot make for itself: the two
            electrons that participate in the Auger transition must be excluded, and which those are depends on the
            transition.  Use ResonantImpactIonization.doubleAugerProbability to obtain a value, read its docstring on
            the two-sided error first, and leave this at 0. to omit the channel entirely rather than guess.

    + initialLevelNo        ::Int64             ... Level No of the initial level from which the capture starts.
    + electronEnergyShift   ::Float64           ... Shift applied to every resonance energy, in the user-selected
                                                     units; the initial level is shifted by the negative amount.
    + dblAugerProbability   ::Float64           ... Probability that a resonance sheds BOTH electrons at once rather
                                                     than one at a time; 0. omits the simultaneous channel.
"""
struct  ResonantIonizationStrengths   <:  Cascade.AbstractSimulationProperty
    initialLevelNo        ::Int64
    electronEnergyShift   ::Float64
    dblAugerProbability   ::Float64
end


"""
`Cascade.ResonantIonizationStrengths()`  ... (simple) constructor for cascade ResonantIonizationStrengths.
"""
function ResonantIonizationStrengths()
    ResonantIonizationStrengths(1, 0., 0.)
end


# `Base.show(io::IO, dist::Cascade.ResonantIonizationStrengths)`  ... prepares a proper printout of the variable dist.
function Base.show(io::IO, dist::Cascade.ResonantIonizationStrengths)
    println(io, "initialLevelNo:           $(dist.initialLevelNo)  ")
    println(io, "electronEnergyShift:      $(dist.electronEnergyShift)  ")
    println(io, "dblAugerProbability:      $(dist.dblAugerProbability)  ")
end


"""
`struct  Cascade.PiRateCoefficients   <:  Cascade.AbstractSimulationProperty`  
    ... defines a type for simulating the PI plasma rate coefficients as function of the (free) photon energy distribution and
        the plasma temperature. These coefficients included a convolution of the direct photoionization cross sections over the 
        Maxwell distribution of electron and (if requested) add the pre-evaluated resonant parts afterwards.
        
        For the implementation: (1) Perform GL integration over free-electron contributions analog to the photoabsorption;
        (2) for the resonant part, simply add convoluted summation terms for all final states (to be worked out in detail);
        (3) use a electronEnergies grid ... instead of the photonEnergies grid.
        (4) in the direct part, sum over the final states before convolution

    + includeResonantPart ::Bool        ... True, if the resonant contributions are to be included.
    + initialLevelNo      ::Int64       ... Level No of initial level for which rate coefficients are to be computed.
    + temperatures        ::Array{Float64,1}
        ... temperatures [K] for which the DR plasma rate coefficieints to be calculated.
    + multipoles          ::Array{EmMultipole}           
        ... Multipoles of the radiation field that are to be included into the photoionization processes.
    + finalConfigurations ::Array{Configuration,1}
    

"""  
struct  PiRateCoefficients   <:  Cascade.AbstractSimulationProperty
    includeResonantPart   ::Bool
    initialLevelNo        ::Int64 
    temperatures          ::Array{Float64,1}
    multipoles            ::Array{EmMultipole} 
    finalConfigurations   ::Array{Configuration,1}
end 


"""
`Cascade.PiRateCoefficients()`  ... (simple) constructor for cascade PiRateCoefficients.
"""
function PiRateCoefficients()
    PiRateCoefficients(false, 1, Float64[], EmMultipole[], Configuration[])
end


# `Base.show(io::IO, dist::Cascade.PiRateCoefficients)`  ... prepares a proper printout of the variable data::Cascade.PiRateCoefficients.
function Base.show(io::IO, dist::Cascade.PiRateCoefficients) 
    println(io, "includeResonantPart:      $(dist.includeResonantPart)  ")
    println(io, "initialLevelNo:           $(dist.initialLevelNo)  ")
    println(io, "temperatures:             $(dist.temperatures)  ")
    println(io, "multipoles:               $(dist.multipoles)  ")
    println(io, "finalConfigurations:      $(dist.finalConfigurations)  ")
end


"""
abstract type  Cascade.AbstractOpacityDependence`  
    ... defines an abstract type to distinguish different dependencies for the opacity; see also:

    + struct FrequencyOpacityDependence     ... to deal with omega-dependence opacities [omega].
    + struct WavelengthOpacityDependence    ... to deal with wavelength-dependence opacities [lambda].
    + struct TemperatureOpacityDependence   ... to deal with temperature-normalized dependent opacities [u].
"""
abstract type  AbstractOpacityDependence      end


"""
`struct  Cascade.FrequencyOpacityDependence   <:  Cascade.AbstractOpacityDependence`  
    ... to deal with omega-dependence opacities [omega] and a binning that need to be given in Hartree.
"""  
struct  FrequencyOpacityDependence   <:  Cascade.AbstractOpacityDependence
    binning    ::Float64
end


"""
`struct  Cascade.WavelengthOpacityDependence   <:  Cascade.AbstractOpacityDependence`  
    ... to deal with wavelength-dependence opacities [lambda] and a binning that need always to be given in nm.
"""  
struct  WavelengthOpacityDependence   <:  Cascade.AbstractOpacityDependence
    binning    ::Float64
end


"""
`struct  Cascade.TemperatureOpacityDependence   <:  Cascade.AbstractOpacityDependence`  
    ... to deal with temperature-normalized dependent opacities [u] and a binning that need to be given Hartree.
"""  
struct  TemperatureOpacityDependence   <:  Cascade.AbstractOpacityDependence
    binning    ::Float64
end


"""
`struct  Cascade.ExpansionOpacities   <:  Cascade.AbstractSimulationProperty`  
    ... defines a type for simulating the expansion opacity as function of the wavelength as well as (parametrically) the density
        and expansion time.

        The expansion opacity is formed in the line-binned Sobolev picture of Eastman & Pinto (1993),

            kappa_exp(lambda) = 1/(c * t_exp * rho) * SUM_l (lambda_l/Delta-lambda) * (1 - exp(-tau_l))

        with the Sobolev optical depth of each bound-bound line l

            tau_l = pi * alpha * n_l * lambda_l * t_exp * f_l         (atomic units; n_l = LOWER level)

        NOTE THAT TWO DIFFERENT DENSITIES ENTER, and they are given as two separate fields on purpose:
        a MASS density rho fixes the prefactor and hence the unit of kappa [cm^2/g], while a NUMBER density
        of the absorbing level enters the optical depth.  A single "ionDensity" field previously served as
        the first while the second was hard-wired to 1.0 in the code, which made the absolute scale of every
        opacity meaningless.  Stating both explicitly keeps an inconsistent pair visible in the printout
        rather than silent.

    + levelPopulation        ::Basics.AbstractLevelPopulation  
        ... the level population used to obtain n_l from ionNumberDensity; BoltzmannLevelPopulation() gives
            the LTE population of the lower level, normalised by the partition function over the levels that
            occur in the given line list.
    + opacityDependence      ::Cascade.AbstractOpacityDependence    
        ... to specify the dependence of the opacities [omega, lambda, (temperature-normalized) u]
    + ionNumberDensity       ::Float64       ... number density n of the ions [in 1/cm^3]; enters tau_l
    + massDensity            ::Float64       ... mass density rho [in g/cm^3]; enters the 1/(rho c t) prefactor
    + temperature            ::Float64       ... temperature [in K]
    + expansionTime          ::Float64       ... (expansion/observation) time [in sec]
    + transitionEnergyShift  ::Float64     
        ... (total) energy shifts that apply to all transition energies; the amplitudes are re-scaled accordingly.
    + dependencyValues       ::Array{Float64,1}
        ... the bin CENTRES [in a.u.] at which the expansion opacity is reported; the bin width is taken from
            opacityDependence.binning.
"""  
struct  ExpansionOpacities      <:  Cascade.AbstractSimulationProperty
    levelPopulation          ::Basics.AbstractLevelPopulation
    opacityDependence        ::Cascade.AbstractOpacityDependence
    ionNumberDensity         ::Float64
    massDensity              ::Float64
    temperature              ::Float64
    expansionTime            ::Float64
    transitionEnergyShift    ::Float64
    dependencyValues         ::Array{Float64,1}
end 


"""
`Cascade.ExpansionOpacities()`  ... (simple) constructor for expansion opacity simulations.
"""
function ExpansionOpacities()
    ExpansionOpacities(Basics.BoltzmannLevelPopulation(), WavelengthOpacityDependence(0.01), 1.0e10, 1.0e-13, 5000., 86400.,  0., Float64[])
end


# `Base.show(io::IO, opacities::Cascade.ExpansionOpacities)`  ... prepares a proper printout of the opacities::Cascade.ExpansionOpacities.
function Base.show(io::IO, opacities::Cascade.ExpansionOpacities) 
    println(io, "levelPopulation:            $(opacities.levelPopulation)  ")
    println(io, "opacityDependence:          $(opacities.opacityDependence)  ")
    println(io, "ionNumberDensity:           $(opacities.ionNumberDensity)  [1/cm^3]")
    println(io, "massDensity:                $(opacities.massDensity)  [g/cm^3]")
    println(io, "temperature:                $(opacities.temperature)  [K]")
    println(io, "expansionTime:              $(opacities.expansionTime)  [sec]")
    println(io, "transitionEnergyShift:      $(opacities.transitionEnergyShift)  ")
    println(io, "dependencyValues:           $(opacities.dependencyValues)  ")
end


"""
`abstract type  Cascade.AbstractOpacityMean`
    ... defines an abstract and a number of singleton types for the MEAN that is taken over a spectral
        opacity kappa_nu. A mean is a FUNCTIONAL of kappa_nu and is applied after the spectral opacity has
        been assembled; it is therefore kept separate from both the contributions that build kappa_nu and
        the dependence (binning) axis on which it is reported.

        Both means below carry the same units as kappa_nu itself, [cm^2/g].

    + struct RosselandMean   ... the HARMONIC mean weighted by dB_nu/dT; use in the optically thick,
                                 diffusion regime.
    + struct PlanckMean      ... the ARITHMETIC mean weighted by B_nu; use in the optically thin,
                                 emission regime.
"""
abstract type  AbstractOpacityMean      end


"""
`struct  Cascade.RosselandMean  <:  Cascade.AbstractOpacityMean`
    ... the Rosseland mean opacity, defined by the HARMONIC average

            1 / kappa_R  =  int (1/kappa_nu) (dB_nu/dT) dnu  /  int (dB_nu/dT) dnu

        In the temperature-normalised variable u = omega/kT the normalised weight is
        w(u) = 15/(4 pi^4) u^4 e^-u / (1 - e^-u)^2, whose prefactor is exact because
        int_0^inf u^4 e^u/(e^u-1)^2 du = 4 Gamma(4) zeta(4) = 4 pi^4/15.

        BECAUSE IT IS HARMONIC, THE ROSSELAND MEAN IS DOMINATED BY THE MOST TRANSPARENT FREQUENCIES --
        radiation leaks through the windows between the lines. One consequence must be understood before
        the number is used: if any bin of the spectral opacity is EMPTY, kappa_nu = 0 there and the mean is
        exactly zero. That is physically correct, not a numerical artefact, and it means a Rosseland mean of
        a bound-bound line list ALONE is meaningless -- a continuum floor (electron scattering, bound-free)
        is required. Cascade.simulateMeanOpacities therefore reports the empty bins rather than returning a
        number that merely looks usable.
"""
struct   RosselandMean          <:  Cascade.AbstractOpacityMean   end


"""
`struct  Cascade.PlanckMean  <:  Cascade.AbstractOpacityMean`
    ... the Planck (emission) mean opacity, defined by the ARITHMETIC average

            kappa_P  =  int kappa_nu B_nu dnu  /  int B_nu dnu

        Being arithmetic it is dominated by the most OPAQUE frequencies, i.e. by the lines, and it stays
        finite on an incomplete line list. It is the appropriate mean for optically thin emission, and it is
        what was computed -- under the name "Rosseland" and with the Rosseland weight -- by every JAC version
        before 14-Aug-2026.
"""
struct   PlanckMean             <:  Cascade.AbstractOpacityMean   end


"""
`abstract type  Cascade.AbstractOpacityContribution`
    ... defines an abstract and a number of singleton types for the CONTRIBUTIONS that build up a spectral
        opacity. Contributions are ADDITIVE, kappa_nu = sum over contributions, and each carries its own
        physical assumptions in its docstring, so that a user reading the type sees what is being computed.
        All contributions are mass opacities in [cm^2/g].

    + struct BoundBoundOpacity   ... lines, in the Sobolev/expansion treatment.
    + struct BoundFreeOpacity    ... photoionization of the bound shells; a true continuum.
    + struct ScatteringOpacity   ... Thomson scattering off free electrons; frequency-independent.

        NOT provided, deliberately rather than as a stub: free-free (bremsstrahlung). See
        Cascade.RosselandMean for why a continuum contribution is needed at all.
"""
abstract type  AbstractOpacityContribution      end


"""
`struct  Cascade.BoundBoundOpacity  <:  Cascade.AbstractOpacityContribution`
    ... the bound-bound (line) contribution in the Sobolev/expansion treatment of Eastman & Pinto,

            kappa_nu^bb = 1/(rho c t_exp) * sum_lines (lambda_l/Delta-lambda) (1 - exp(-tau_l))

        with the Sobolev optical depth tau_l = pi alpha n_l lambda_l t_exp f_l.

        ASSUMPTIONS: a homologously expanding medium in which each line is a sharp resonance swept through
        by the velocity gradient (Sobolev); LTE level populations; and a line list that is complete enough
        for the chosen binning. The last is not a formality -- an empty bin means "no line in this bin", NOT
        "no opacity", and it is what makes a Rosseland mean of lines alone vanish.

    + levelPopulation  ::Basics.AbstractLevelPopulation
        ... the level population used to obtain n_l; BoltzmannLevelPopulation() gives the LTE population of
            the lower level, normalised by the partition function over the levels of the given line list.
"""
struct   BoundBoundOpacity      <:  Cascade.AbstractOpacityContribution
    levelPopulation          ::Basics.AbstractLevelPopulation
end


"""
`struct  Cascade.BoundFreeOpacity  <:  Cascade.AbstractOpacityContribution`
    ... the bound-free (photoionization) contribution,

            kappa_nu^bf = n_ion * sum_shells sigma^PI_shell(omega) / rho ,

        with the cross sections taken from Empirical.photoionizationCrossSection in the ScaledHydrogenic
        approximation: Kramers' (1923) formula scaled by the tabulated binding energy of the ionized shell
        and corrected by a bound-free Gaunt factor, which for H 1s turns the Kramers 7.91 Mb at threshold
        into the exact 6.30 Mb and steepens the tail from omega^-3 to the exact omega^-7/2.

        UNLIKE the bound-bound contribution this is a genuine CONTINUUM: it is evaluated at each node
        directly and carries no lambda/Delta-lambda bin factor, since there is no line to be smeared over a
        bin. Each shell contributes only above its own threshold, so the bound-free opacity has edges.

        ASSUMPTIONS: an independent-particle, hydrogenic-scaled description of a single active electron per
        shell; tabulated binding energies; and every ion in the configuration given, i.e. no ionization
        balance is solved here -- the caller states which ion stage is present through the configurations.

        KNOW THIS BEFORE USING IT ON AN ION. Empirical.scaledBindingEnergy draws on X-ray tables compiled for
        NEUTRAL atoms, so the threshold is the neutral binding energy whatever ion stage the configurations
        describe. Measured 14-Aug-2026: for Sr it returns 5.573 eV, essentially the neutral Sr I -> Sr II
        potential of 5.695 eV, where the true edge of Sr^+ is the Sr II -> Sr III potential of 11.030 eV --
        low by a factor of about two. The edge is therefore placed at too LOW a photon energy for any ion,
        which lets bound-free opacity appear in a band where the ion cannot in fact be ionized. For neutral
        species the tables are exactly what is wanted.

        AND THE OBVIOUS ESCAPE DOES NOT WORK. Empirical.photoionizationCrossSection has a second
        approximation, UsingJAC, which computes the cross section ab initio from mean-field orbitals and E1
        amplitudes -- but its own docstring records that ITS thresholds are the DFS orbital energies and
        likewise lie below the true ionization potentials (Ne 2p at 12.0 eV against the true 21.6 eV). Both
        routes therefore open their edge too early, for different reasons: a table that ignores the ion
        stage, and a Koopmans-like eigenvalue that is not an ionization energy. A threshold that is right
        for an ion needs the DIFFERENCE of two total energies, E(ion^q+1) - E(ion^q), which JAC can compute
        but neither cross-section route uses. Until then, treat the position of a bound-free edge here as
        approximate and the opacity within one edge-width of it as unreliable.

    + nuclearCharge          ::Float64
        ... the nuclear charge Z. THIS IS CARRIED EXPLICITLY ON PURPOSE. Empirical.photoionizationCrossSection
            reads Z from the global Defaults, whose GBL_NUCLEAR_CHARGE is 1.0 unless something sets it, and
            no JAC computation does. A contribution relying on that global would silently return HYDROGEN
            cross sections for every element; here Z is taken from this field, set around the call and
            restored afterwards.
    + ionizationPairs        ::Array{Tuple{Configuration,Configuration},1}
        ... the (initial, final) configuration pairs whose photoionization is to be summed; the two must
            differ by exactly one electron in one shell, which is the shell that is ionized.
"""
struct   BoundFreeOpacity       <:  Cascade.AbstractOpacityContribution
    nuclearCharge            ::Float64
    ionizationPairs          ::Array{Tuple{Configuration,Configuration},1}
end


"""
`struct  Cascade.ScatteringOpacity  <:  Cascade.AbstractOpacityContribution`
    ... the Thomson scattering contribution of free electrons,

            kappa^sc = n_e sigma_T / rho ,      sigma_T = 6.6524587321e-25 cm^2

        which is INDEPENDENT of frequency and therefore contributes the same value to every bin. It is small
        beside a strong line, but it is the CONTINUUM FLOOR that makes a Rosseland mean well defined at all:
        without it a single empty bin sends the harmonic mean to zero.

        ASSUMPTIONS: non-relativistic, unpolarised Thomson scattering (kT << m_e c^2, i.e. T << 6e9 K), and
        scattering treated as pure extinction. For a fully ionised hydrogen plasma n_e = rho/m_H and this
        reduces to the textbook kappa_es = 0.20 (1 + X) = 0.40 cm^2/g, which is the cleanest known-answer
        check available for any mean-opacity implementation.

    + electronNumberDensity  ::Float64      ... free-electron number density n_e [in 1/cm^3]
"""
struct   ScatteringOpacity      <:  Cascade.AbstractOpacityContribution
    electronNumberDensity    ::Float64
end


"""
`struct  Cascade.MeanOpacities   <:  Cascade.AbstractSimulationProperty`
    ... defines a type for simulating a MEAN opacity as function of the temperature or density as well as
        (parametrically) the expansion time. Usually, different values are given for either the temperature
        or density. A value::EmProperty in [cm^2/g] is reported for each (rho, T) pair.

        This type replaces the former `RosselandOpacities` (14-Aug-2026), which had never been run: it took
        only ONE density list, documented as a mass density, and derived a number density from it as
        rho/1.6605e-24 -- i.e. one ion per nucleon mass, A = 1 for every element, wrong by roughly A. Both
        densities are now carried explicitly, exactly as Cascade.ExpansionOpacities already does and for the
        same reason. It also fixed the mean to a single formula that was in fact the arithmetic one; the
        mean is now named, and is a field.

    + opacityMean            ::Cascade.AbstractOpacityMean
        ... the mean to be taken over the spectral opacity: RosselandMean() or PlanckMean().
    + contributions          ::Array{Cascade.AbstractOpacityContribution,1}
        ... what the spectral opacity is made of; the contributions are summed bin by bin. A list holding
            only BoundBoundOpacity reproduces the pre-14-Aug-2026 behaviour, and cannot support a Rosseland
            mean; adding ScatteringOpacity supplies the continuum floor that can.
    + opacityDependence      ::Cascade.TemperatureOpacityDependence
        ... to specify the dependence of the opacities in terms of (the temperature-normalized) u
    + ionNumberDensities     ::Array{Float64,1}      ... list of ion NUMBER densities [in 1/cm^3]; enters tau_l
    + massDensities          ::Array{Float64,1}      ... list of MASS densities [in g/cm^3]; enters 1/(rho c t)
    + temperatures           ::Array{Float64,1}      ... list of temperatures [in K]
    + expansionTime          ::Float64               ... (expansion/observation) time [in sec]
    + transitionEnergyShift  ::Float64
        ... (total) energy shifts that apply to all transition energies; the amplitudes are re-scaled accordingly.
"""
struct  MeanOpacities           <:  Cascade.AbstractSimulationProperty
    opacityMean              ::Cascade.AbstractOpacityMean
    contributions            ::Array{Cascade.AbstractOpacityContribution,1}
    opacityDependence        ::Cascade.TemperatureOpacityDependence
    ionNumberDensities       ::Array{Float64,1}
    massDensities            ::Array{Float64,1}
    temperatures             ::Array{Float64,1}
    expansionTime            ::Float64
    transitionEnergyShift    ::Float64
end


"""
`Cascade.MeanOpacities()`  ... (simple) constructor for mean-opacity simulations.
"""
function MeanOpacities()
    MeanOpacities(RosselandMean(), Cascade.AbstractOpacityContribution[ BoundBoundOpacity(Basics.BoltzmannLevelPopulation()) ],
                  TemperatureOpacityDependence(0.01), [1.0e10], [1.0e-13], [5000.], 86400.,  0.)
end


# `Base.show(io::IO, opacities::Cascade.MeanOpacities)`
#   ... prepares a proper printout of opacities::Cascade.MeanOpacities.
function Base.show(io::IO, opacities::Cascade.MeanOpacities)
    println(io, "opacityMean:                $(opacities.opacityMean)  ")
    println(io, "contributions:              $(opacities.contributions)  ")
    println(io, "opacityDependence:          $(opacities.opacityDependence)  ")
    println(io, "ionNumberDensities:         $(opacities.ionNumberDensities)  [1/cm^3]")
    println(io, "massDensities:              $(opacities.massDensities)  [g/cm^3]")
    println(io, "temperatures:               $(opacities.temperatures)  [K]")
    println(io, "expansionTime:              $(opacities.expansionTime)  [sec]")
    println(io, "transitionEnergyShift:      $(opacities.transitionEnergyShift)  ")
end


"""
`struct  Cascade.ElectronIntensities   <:  Cascade.AbstractSimulationProperty`  
    ... defines a type for simulating the electron-line intensities as function of electron energy.

    + minElectronEnergy   ::Float64     ... Minimum electron energy for the simulation of electron spectra.
    + maxElectronEnergy   ::Float64     ... Maximum electron energy for the simulation of electron spectra.
    + initialOccupations  ::Array{Tuple{Int64,Float64},1}   
        ... List of one or several (tupels of) levels in the overall cascade tree together with their relative population.
"""  
struct  ElectronIntensities   <:  Cascade.AbstractSimulationProperty
    minElectronEnergy     ::Float64
    maxElectronEnergy     ::Float64
    initialOccupations    ::Array{Tuple{Int64,Float64},1} 
end 


"""
`Cascade.ElectronIntensities()`  ... (simple) constructor for cascade ElectronIntensities.
"""
function ElectronIntensities()
    ElectronIntensities(0., 1.0e6,  [(1, 1.0)])
end


# `Base.show(io::IO, data::Cascade.ElectronIntensities)`  ... prepares a proper printout of the variable data::Cascade.ElectronIntensities.
function Base.show(io::IO, dist::Cascade.ElectronIntensities) 
    println(io, "minElectronEnergy:        $(dist.minElectronEnergy)  ")
    println(io, "maxElectronEnergy:        $(dist.maxElectronEnergy)  ")
    println(io, "initialOccupations:       $(dist.initialOccupations)  ")
end


"""
`struct  Cascade.FinalLevelDistribution   <:  Cascade.AbstractSimulationProperty`  
    ... defines a type for simulating the 'final-level distribution' as it is found after all cascade processes are completed.

    + initialOccupations  ::Array{Tuple{Int64,Float64},1}   
        ... List of one or several (tupels of) levels in the overall cascade tree together with their relative population.
    + leadingConfigs      ::Array{Configuration,1}   
        ... List of leading configurations whose levels are equally populated, either initially or ....
    + finalConfigs        ::Array{Configuration,1}   
        ... List of final configurations whose level population are to be "listed" finally; these configuration
            only determine the printout but not the propagation of the probabilities. If this list is empty, all 
            the levels are shown.
"""  
struct  FinalLevelDistribution   <:  Cascade.AbstractSimulationProperty
    initialOccupations    ::Array{Tuple{Int64,Float64},1} 
    leadingConfigs        ::Array{Configuration,1}
    finalConfigs          ::Array{Configuration,1}
end 


"""
`Cascade.FinalLevelDistribution()`  ... (simple) constructor for cascade FinalLevelDistribution.
"""
function FinalLevelDistribution()
    FinalLevelDistribution([(1, 1.0)], Configuration[], Configuration[])
end


# `Base.show(io::IO, data::Cascade.FinalLevelDistribution)`  
#       ... prepares a proper printout of the variable data::Cascade.FinalLevelDistribution.
function Base.show(io::IO, dist::Cascade.FinalLevelDistribution) 
    println(io, "initialOccupations:       $(dist.initialOccupations)  ")
    println(io, "leadingConfigs:           $(dist.leadingConfigs)  ")
    println(io, "finalConfigs:             $(dist.finalConfigs)  ")
end


"""
`struct  Cascade.IonDistribution   <:  Cascade.AbstractSimulationProperty`  
    ... defines a type for simulating the 'ion distribution' as it is found after all cascade processes are completed.

    + initialOccupations  ::Array{Tuple{Int64,Float64},1}   
        ... List of one or several (tupels of) levels in the overall cascade tree together with their relative population.
    + leadingConfigs      ::Array{Configuration,1}   
        ... List of leading configurations whose levels are equally populated, either initially or ....
"""  
struct  IonDistribution   <:  Cascade.AbstractSimulationProperty
    initialOccupations    ::Array{Tuple{Int64,Float64},1} 
    leadingConfigs        ::Array{Configuration,1}
end 


"""
`Cascade.IonDistribution()`  ... (simple) constructor for cascade IonDistribution data.
"""
function IonDistribution()
    IonDistribution([(1, 1.0)], Configuration[])
end


# `Base.show(io::IO, data::Cascade.IonDistribution)`  ... prepares a proper printout of the variable data::Cascade.IonDistribution.
function Base.show(io::IO, dist::Cascade.IonDistribution) 
    println(io, "initialOccupations:       $(dist.initialOccupations)  ")
    println(io, "leadingConfigs:           $(dist.leadingConfigs)  ")
end


"""
`struct  Cascade.RelaxationCurve   <:  Cascade.AbstractSimulationProperty`  
    ... defines a type for simulating the mean relaxation times in which 70%, 80%, of the occupied levels decay to the 
        ground configuration.

    + timeStep            ::Float64     ... Time-step for following the decay of levels [a.u.]
    + initialOccupations  ::Array{Tuple{Int64,Float64},1}   
        ... List of one or several (tupels of) levels in the overall cascade tree together with their relative population.
    + leadingConfigs      ::Array{Configuration,1}   
        ... List of leading configurations whose levels are equally populated, either initially or ....
    + groundConfigs       ::Array{Configuration,1}   
        ... List of ground configurations into which the decay is considered.
"""  
struct  RelaxationCurve   <:  Cascade.AbstractSimulationProperty
    timeStep              ::Float64
    initialOccupations    ::Array{Tuple{Int64,Float64},1} 
    leadingConfigs        ::Array{Configuration,1}
    groundConfigs         ::Array{Configuration,1}
end 


"""
`Cascade.RelaxationCurve()`  ... (simple) constructor for cascade RelaxationCurve.
"""
function RelaxationCurve()
    RelaxationCurve(10., [(1, 1.0)], Configuration[], Configuration[])
end


# `Base.show(io::IO, dist::Cascade.RelaxationCurve)`  ... prepares a proper printout of the variable data::Cascade.RelaxationCurve.
function Base.show(io::IO, dist::Cascade.RelaxationCurve) 
    println(io, "timeStep:                 $(dist.timeStep)  ")
    println(io, "initialOccupations:       $(dist.initialOccupations)  ")
    println(io, "leadingConfigs:           $(dist.leadingConfigs)  ")
    println(io, "groundConfigs:            $(dist.groundConfigs)  ")
end


"""
`struct  Cascade.PhotonIntensities   <:  Cascade.AbstractSimulationProperty`  
    ... defines a type for simulating the photon-line intensities as function of photon energy.

    + minPhotonEnergy     ::Float64     ... Minimum photon energy for the simulation of photon spectra.
    + maxPhotonEnergy     ::Float64     ... Maximum photon energy for the simulation of photon spectra.
    + initialOccupations  ::Array{Tuple{Int64,Float64},1}   
        ... List of one or several (tupels of) levels in the overall cascade tree together with their relative population.
    + leadingConfigs      ::Array{Configuration,1}   
        ... List of leading configurations whose levels are equally populated, either initially or ....
"""  
struct  PhotonIntensities   <:  Cascade.AbstractSimulationProperty
    minPhotonEnergy       ::Float64
    maxPhotonEnergy       ::Float64
    initialOccupations    ::Array{Tuple{Int64,Float64},1} 
    leadingConfigs        ::Array{Configuration,1}
end 


"""
`Cascade.PhotonIntensities()`  ... (simple) constructor for cascade PhotonIntensities.
"""
function PhotonIntensities()
    PhotonIntensities(0., 1.0e6,  [(1, 1.0)], Configuration[])
end


# `Base.show(io::IO, dist::Cascade.PhotonIntensities)`  ... prepares a proper printout of the variable data::Cascade.PhotonIntensities.
function Base.show(io::IO, dist::Cascade.PhotonIntensities) 
    println(io, "minPhotonEnergy:          $(dist.minPhotonEnergy)  ")
    println(io, "maxPhotonEnergy:          $(dist.maxPhotonEnergy)  ")
    println(io, "initialOccupations:       $(dist.initialOccupations)  ")
    println(io, "leadingConfigs:           $(dist.leadingConfigs)  ")
end


"""
`struct  Cascade.PhotoAbsorptionSpectrum   <:  Cascade.AbstractSimulationProperty`  
    ... defines a type for simulating the total photo-absorption cross sections in a given interval of photon energies
        as well as for a given set of photo-ionization and photo-excitation cross sections

    + includeIonization   ::Bool             ... True, if photo-ionization cross sections are to be considered.
    + includeExcitation   ::Bool             ... True, if photo-excitation lines are to be considered.
    + resonanceWidth      ::Float64          ... Widths of the resonances (user-defined units)
    + csScaling           ::Float64          ... Scaling factor do enhance the strengths of resonances (default=1.0)
    + photonEnergies      ::Array{Float64,1} ... Photon energies (in user-selected units) for the simulation of photon spectra
                                                    to describe the interval and resolution of the absorption cross sections.
    + shells              ::Array{Shell,1}   
        ... Shells that should be included for a partial absorption cross sections; a cross section contribution is considered,
            if the occupation of one of these shells is lowered in the leading configurations.
    + initialOccupations  ::Array{Tuple{Int64,Float64},1}   
        ... List of one or several (tupels of) levels in the overall cascade tree together with their relative population;
            at least one of these tuples must be given.
    + leadingConfigs      ::Array{Configuration,1}   
        ... List of leading configurations whose levels are equally populated, either initially or ....
"""  
struct  PhotoAbsorptionSpectrum   <:  Cascade.AbstractSimulationProperty
    includeIonization     ::Bool
    includeExcitation     ::Bool
    resonanceWidth        ::Float64
    csScaling             ::Float64
    photonEnergies        ::Array{Float64,1}
    shells                ::Array{Shell,1}
    initialOccupations    ::Array{Tuple{Int64,Float64},1} 
    leadingConfigs        ::Array{Configuration,1}
end 


"""
`Cascade.PhotoAbsorptionSpectrum()`  ... (simple) constructor for cascade PhotoAbsorptionSpectrum.
"""
function PhotoAbsorptionSpectrum()
    PhotoAbsorptionSpectrum(true, false, 0.1, 1., [1.0],  Shell[], [(1, 1.0)], Configuration[])
end


# `Base.show(io::IO, dist::Cascade.PhotoAbsorptionSpectrum)`  ... prepares a proper printout of the variable data::Cascade.PhotoAbsorptionSpectrum.
function Base.show(io::IO, dist::Cascade.Cascade.PhotoAbsorptionSpectrum) 
    println(io, "includeIonization:        $(dist.includeIonization)  ")
    println(io, "includeExcitation:        $(dist.includeExcitation)  ")
    println(io, "resonanceWidth:           $(dist.resonanceWidth)  ")
    println(io, "csScaling:                $(dist.csScaling)  ")
    println(io, "photonEnergies:           $(dist.photonEnergies)  ")
    println(io, "shells:                   $(dist.shells)  ")
    println(io, "initialOccupations:       $(dist.initialOccupations)  ")
    println(io, "leadingConfigs:           $(dist.leadingConfigs)  ")
end


struct   DecayPathes                  <:  AbstractSimulationProperty   end
struct   ElectronCoincidence          <:  AbstractSimulationProperty   end
struct   TimeBinnedPhotonIntensity    <:  AbstractSimulationProperty   end
struct   MeanLineWidth                <:  AbstractSimulationProperty   end
struct   PhotoResonances              <:  AbstractSimulationProperty   end



"""
abstract type  Cascade.AbstractSimulationMethod`  
    ... defines a abstract and a list of singleton data types for the properties that can be 'simulated' from a given
        list of lines.

    + struct ProbPropagation     ... to propagate the (occupation) probabilites of the levels until no further changes occur.
    + struct MonteCarlo          ... to simulate the cascade decay by a Monte-Carlo approach of possible pathes (not yet considered).
    + struct RateEquations       ... to solve the cascade by a set of rate equations (not yet considered).
"""
abstract type  AbstractSimulationMethod                  end
struct   ProbPropagation  <:  AbstractSimulationMethod   end
struct   MonteCarlo       <:  AbstractSimulationMethod   end
struct   RateEquations    <:  AbstractSimulationMethod   end


"""
`struct  Cascade.SimulationSettings`  ... defines settings for performing the simulation of some cascade (data).

    + printTree           ::Bool        ... Print the cascade tree in a short form
    + printLongTree       ::Bool        ... Print the cascade tree in a long form
    + initialPhotonEnergy ::Float64     ... Photon energy for which photoionization data are considered. 
"""
struct  SimulationSettings
    printTree             ::Bool
    printLongTree         ::Bool 
    initialPhotonEnergy   ::Float64
end 


"""
`Cascade.SimulationSettings()`  ... constructor for an 'empty' instance of a Cascade.Block.
"""
function SimulationSettings()
    SimulationSettings(false, false, 0.)
end


# `Base.show(io::IO, settings::SimulationSettings)`  ... prepares a proper printout of the variable settings::SimulationSettings.
function Base.show(io::IO, settings::SimulationSettings) 
    println(io, "printTree:                $(settings.printTree)  ")
    println(io, "printLongTree:            $(settings.printLongTree)  ")
    println(io, "initialPhotonEnergy:      $(settings.initialPhotonEnergy)  ")
end


"""
`struct  Cascade.Simulation`  ... defines a simulation on some given cascade (data).

    + name            ::String                              ... Name of the simulation
    + property        ::Cascade.AbstractSimulationProperty 
        ... Property that is to be considered in this simulation of the cascade (data).
    + method          ::Cascade.AbstractSimulationMethod    
        ... Method that is used in the cascade simulation; cf. Cascade.SimulationMethod.
    + settings        ::Cascade.SimulationSettings          ... Settings for performing these simulations.
    + computationData ::Array{Dict{String,Any},1}           ... Date on which the simulations are performed
"""
struct  Simulation
    name              ::String
    property          ::Cascade.AbstractSimulationProperty
    method            ::Cascade.AbstractSimulationMethod
    settings          ::Cascade.SimulationSettings 
    computationData   ::Array{Dict{String,Any},1}
end 


"""
`Cascade.Simulation()`  ... constructor for an 'default' instance of a Cascade.Simulation.
"""
function Simulation()
    Simulation("Default cascade simulation", Cascade.PhotoAbsorptionSpectrum(), Cascade.ProbPropagation(), 
                Cascade.SimulationSettings(), Array{Dict{String,Any},1}[] )
end


"""
`Cascade.Simulation(sim::Cascade.Simulation;`
    
            name=..,               property=..,             method=..,              settings=..,     computationData=.. )
            
    ... constructor for re-defining the computation::Cascade.Simulation.
"""
function Simulation(sim::Cascade.Simulation;                              
    name::Union{Nothing,String}=nothing,                                  property::Union{Nothing,Cascade.AbstractSimulationProperty}=nothing,
    method::Union{Nothing,Cascade.AbstractSimulationMethod}=nothing,      settings::Union{Nothing,Cascade.SimulationSettings}=nothing,    
    computationData::Union{Nothing,Array{Dict{String,Any},1}}=nothing )

    if  isnothing(name)              namex            = sim.name              else  namex            = name                end 
    if  isnothing(property)          propertyx        = sim.property          else  propertyx        = property            end 
    if  isnothing(method)            methodx          = sim.method            else  methodx          = method              end 
    if  isnothing(settings)          settingsx        = sim.settings          else  settingsx        = settings            end 
    if  isnothing(computationData)   computationDatax = sim.computationData   else  computationDatax = computationData     end 
    
    Simulation(namex, propertyx, methodx, settingsx, computationDatax)
end


# `Base.show(io::IO, simulation::Cascade.Simulation)`  ... prepares a proper printout of the variable simulation::Cascade.Simulation.
function Base.show(io::IO, simulation::Cascade.Simulation) 
    println(io, "name:              $(simulation.name)  ")
    println(io, "property:          $(simulation.property)  ")
    println(io, "method:            $(simulation.method)  ")
    println(io, "computationData:   ... based on $(length(simulation.computationData)) cascade data sets ")
    println(io, "> settings:      \n$(simulation.settings)  ")
end


"""
`struct  Cascade.LineReference{T}`  ... refers to one line within one of the line lists of a Cascade.Data set: it keeps
        the list itself, the atomic process it belongs to, and the position of the line within that list. Renamed
        from LineIndex (05-Aug-2026), which suggested a bare integer index rather than the reference it is.

    + lines        ::Array{T,1}              ... refers to the line list for which this index is defined.
    + process      ::Basics.AbstractProcess  ... refers to the particular lineList of cascade (data).
    + index        ::Int64                   ... index of the corresponding line.
"""
struct  LineReference{T}
    lines          ::Array{T,1}
    process        ::Basics.AbstractProcess
    index          ::Int64
end


"""
`abstract type Cascade.AbstractCoincidenceGate`
    ... specifies one condition that an emitted particle (or the resulting ion) has to satisfy for a cascade
        event to be counted in a coincidence measurement. The gates are named after the PARTICLE that has to
        be detected, since that is what distinguishes them; ElectronGate and PhotonGate then differ only in
        which of the two an experiment's analyser accepted.

    + struct ElectronGate     ... an emitted electron within [minEnergy, maxEnergy].
    + struct PhotonGate       ... an emitted photon within [minEnergy, maxEnergy].
    + struct ChargeStateGate  ... the ion left with exactly NoElectrons electrons at the end of the cascade;
                                  the odd one out, since it gates on the ion rather than on an emitted
                                  particle, and therefore carries a different field.
"""
abstract type  AbstractCoincidenceGate                                                            end
struct  ElectronGate     <:  Cascade.AbstractCoincidenceGate;  minEnergy::Float64;  maxEnergy::Float64   end
struct  PhotonGate       <:  Cascade.AbstractCoincidenceGate;  minEnergy::Float64;  maxEnergy::Float64   end
struct  ChargeStateGate  <:  Cascade.AbstractCoincidenceGate;  NoElectrons::Int64                        end


"""
`struct  Cascade.ParticleCoincidences   <:  Cascade.AbstractSimulationProperty`
    ... defines a type for simulating a coincidence measurement on a cascade: the spectrum of one emitted
        particle, recorded only for those cascade events in which a set of further conditions was met.

    + gates               ::Array{Cascade.AbstractCoincidenceGate,1}
        ... the coincidence conditions. They are matched IN EMISSION ORDER along each decay pathway: gates[1]
            against the first emission of the pathway, gates[2] against the second, and so on. Ordered
            matching is chosen because a cascade pathway is intrinsically ordered and the rule is then
            unambiguous; for the usual cases it makes no difference, since the gated groups are separated far
            beyond any analyser window (a Mg KLL electron near 1100 eV against an LMM electron of a few tens
            of eV).
    + observed            ::Cascade.AbstractCoincidenceGate
        ... the particle whose spectrum is reported, together with the window it is observed in. It is matched
            against the emission FOLLOWING the last gate. Giving a ChargeStateGate here is meaningless -- that
            is an ion distribution, not a spectrum -- and is rejected.
    + initialOccupations  ::Array{Tuple{Int64,Float64},1}
        ... levels of the cascade tree together with their initial relative population.

        One type covers every two- and multi-fold combination, because the gates and the observed particle are
        described in the same way:

            Auger-Auger      gates = [ElectronGate(1050., 1200.)]     observed = ElectronGate(0., 100.)
            photon-photon    gates = [PhotonGate(...)]                observed = PhotonGate(...)
            electron-photon  gates = [ElectronGate(...)]              observed = PhotonGate(...)
            electron-ion     gates = [ChargeStateGate(9)]             observed = ElectronGate(0., 2000.)
            triple           gates = [ElectronGate(...), ElectronGate(...)]  observed = ElectronGate(...)
            singles          gates = []                               observed = ElectronGate(0., 1.0e6)

        The last line reproduces Cascade.ElectronIntensities and is worth running as a check on the machinery.
        Cascade.PhotonIntensities and Cascade.ElectronIntensities are deliberately kept as the simple ungated
        cases: they are the common request, and having both lets each verify the other.

        LIMITATION, and it is a physical one rather than an implementation gap. The intensities computed here
        are an INCOHERENT sum over decay pathways, so they are energy correlations only. They are not angular
        correlations between the emitted particles, and they do not describe interference between pathways
        that reach the same final state by different routes. Fritzsche et al., Symmetry 13, 520 (2021),
        Sect. 2.3 makes the point explicitly: once the fine-structure splitting of intermediate levels becomes
        comparable to their natural widths, "the angular emission of the emitted Auger electrons can no longer
        be described by an incoherent summation over individual pathways". Coincidence experiments are often
        precisely where that matters.
"""
struct  ParticleCoincidences   <:  Cascade.AbstractSimulationProperty
    gates                 ::Array{Cascade.AbstractCoincidenceGate,1}
    observed              ::Cascade.AbstractCoincidenceGate
    initialOccupations    ::Array{Tuple{Int64,Float64},1}
end


"""
`Cascade.ParticleCoincidences()`  ... (simple) constructor for cascade ParticleCoincidences.
"""
function ParticleCoincidences()
    ParticleCoincidences( Cascade.AbstractCoincidenceGate[], Cascade.ElectronGate(0., 1.0e6), [(1, 1.0)] )
end


# `Base.show(io::IO, prop::Cascade.ParticleCoincidences)`  ... prepares a proper printout.
function Base.show(io::IO, prop::Cascade.ParticleCoincidences)
    println(io, "gates:                    $(prop.gates)  ")
    println(io, "observed:                 $(prop.observed)  ")
    println(io, "initialOccupations:       $(prop.initialOccupations)  ")
end


"""
`struct  Cascade.LineFlux`
    ... records how much probability has flowed through ONE line during a cascade propagation, together with
        everything needed to interpret it afterwards. This is what separates the propagation from the
        extraction of observables: Cascade.propagateProbability! records a LineFlux for every line it passes
        probability through, and each simulation property then derives what it needs from that record.

    + process       ::Basics.AbstractProcess  ... the process of the line (Radiative(), Auger(), Photo()).
    + energy        ::Float64                 ... energy carried away by the emitted particle, i.e. the
                                                  difference of the initial and final level energies.
    + flux          ::Float64                 ... probability that has flowed through this line.
    + initialIndex  ::Int64                   ... index of the initial level within the propagated level list.
    + finalIndex    ::Int64                   ... index of the final level within the propagated level list.

        Before 05-Aug-2026 the propagation collected observables itself, through the keywords
        collectPhotonIntensities and collectElectronIntensities, and refused to collect both at once
        (error("stop a")). That made a photon-electron coincidence inexpressible and forced a second full
        propagation whenever a different observable was wanted. Recording the flux once, and keeping the level
        indices with it, removes both restrictions -- and the indices are what a later coincidence simulation
        needs in order to chain individual emissions into pathways.
"""
struct  LineFlux
    process        ::Basics.AbstractProcess
    energy         ::Float64
    flux           ::Float64
    initialIndex   ::Int64
    finalIndex     ::Int64
end


# `Base.show(io::IO, lf::Cascade.LineFlux)`  ... prepares a proper printout of the variable lf::Cascade.LineFlux.
function Base.show(io::IO, lf::Cascade.LineFlux)
    println(io, "process:       $(lf.process)  ")
    println(io, "energy:        $(lf.energy)  ")
    println(io, "flux:          $(lf.flux)  ")
    println(io, "initialIndex:  $(lf.initialIndex)  ")
    println(io, "finalIndex:    $(lf.finalIndex)  ")
end 


# `Base.show(io::IO, index::Cascade.LineReference)`  ... prepares a proper printout of the variable index::Cascade.LineReference.
function Base.show(io::IO, index::Cascade.LineReference) 
    println(io, "lines (typeof):        $(typeof(lines))  ")
    println(io, "process:               $(index.process)  ")
    println(io, "index:                 $(index.index)  ")
end


"""
`mutable struct  Cascade.Level`  ... defines a level specification for dealing with cascade transitions.

    + energy       ::Float64                     ... energy of the level.
    + J            ::AngularJ64                  ... total angular momentum of the level
    + parity       ::Basics.Parity               ... total parity of the level
    + NoElectrons  ::Int64                       ... total number of electrons of the ion to which this level belongs.
    + majorConfig  ::Configuration               ... major (dominant) configuration of this level.
    + relativeOcc  ::Float64                     ... relative occupation  
    + parents      ::Array{Cascade.LineReference,1}  ... list of parent lines that (may) populate the level.     
    + daughters    ::Array{Cascade.LineReference,1}  ... list of daughter lines that (may) de-populate the level.     
"""
mutable struct  Level
    energy         ::Float64 
    J              ::AngularJ64 
    parity         ::Basics.Parity 
    NoElectrons    ::Int64 
    majorConfig    ::Configuration 
    relativeOcc    ::Float64 
    parents        ::Array{Cascade.LineReference,1} 
    daughters      ::Array{Cascade.LineReference,1} 
end 


# `Base.show(io::IO, level::Cascade.Level)`  ... prepares a proper printout of the variable level::Cascade.Level.
function Base.show(io::IO, level::Cascade.Level) 
    println(io, "energy:        $(level.energy)  ")
    println(io, "J:             $(level.J)  ")
    println(io, "parity:        $(level.parity)  ")
    println(io, "NoElectrons:   $(level.NoElectrons)  ")
    println(io, "majorConfig:   $(level.majorConfig)  ")
    println(io, "relativeOcc:   $(level.relativeOcc)  ")
    println(io, "parents:       $(level.parents)  ")
    println(io, "daughters:     $(level.daughters)  ")
end


"""
`Base.:(==)(leva::Cascade.Level, levb::Cascade.Level)`  ... returns true if both levels are equal and false otherwise.
"""
function  Base.:(==)(leva::Cascade.Level, levb::Cascade.Level)
    if  leva.energy == levb.energy   &&   leva.J == levb.J  && leva.parity == levb.parity  && leva.NoElectrons == levb.NoElectrons 
            return( true )
    else    return( false )
    end
end

"""
`struct  Cascade.AbsorptionCrossSection`  
    ... defines the absorption cross section for a particular (incident) photon energy in terms of its discrete and 
        (direct photoionization) contributions. Of course, this absorption cross section depends on the relative population
        of the initial levels.

    + photonEnergy ::Float64              ... incident photon energy/photon-energy dependence of the absorption spectrum.
    + excitationCS ::Basics.EmProperty    ... contribution due to discrete excitation processes.
    + ionizationCS ::Basics.EmProperty    ... contribution due to contineous ionization processes.
"""
struct  AbsorptionCrossSection
    photonEnergy   ::Float64
    excitationCS   ::Basics.EmProperty
    ionizationCS   ::Basics.EmProperty
end 


# `Base.show(io::IO, cs::Cascade.AbsorptionCrossSection)`  ... prepares a proper printout of the variable cs::Cascade.AbsorptionCrossSection.
function Base.show(io::IO, cs::Cascade.AbsorptionCrossSection) 
    println(io, "photonEnergy:      $(cs.photonEnergy)  ")
    println(io, "excitationCS:      $(cs.excitationCS)  ")
    println(io, "ionizationCS:      $(cs.ionizationCS)  ")
end
