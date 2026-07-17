
#################################################################################################################################
### Common helpers ##############################################################################################################


"""
`Empirical.boundFreeGauntFactor(n::Int64, x::Float64)`
    ... to provide the bound-free Gaunt factor g^(bf) for the photoionization of shell n at the photon energy
        x = omega/omega_threshold, i.e. the quantum correction to Kramers' semiclassical cross section,
        sigma = sigma^(Kramers) * g^(bf). For n = 1 the factor is exact (nonrelativistic): it is the ratio of
        Stobbe's (1930) closed-form 1s cross section,
            sigma^(Stobbe)(x) = 2^9 pi^2/(3 e^4) alpha a_o^2 x^(-4) exp(4 - 4 arctan(kappa)/kappa) / (1 - e^(-2 pi/kappa))
        with kappa = sqrt(x-1), to the Kramers form: at threshold, sigma^(Stobbe) = 6.304 Mb (Z = 1) analytically,
        so that g(1) = 6.304/7.907 = 0.7973, while at large x the exact cross section falls off as x^(-7/2) against
        Kramers' x^(-3), so that g ~ x^(-1/2). For n >= 2 the deviation from unity is scaled with the leading-order
        n-dependence of the Menzel & Pekeris (1935) expansion, g = 1 + [g_1s(x) - 1] n^(-2/3), which preserves the
        exact n = 1 limit and approaches the correct semiclassical limit g -> 1 for large n, where Kramers' formula
        becomes exact. A gbf::Float64 is returned; g = 1 is returned for x <= 1, where no cross section remains to
        be corrected.
"""
function boundFreeGauntFactor(n::Int64, x::Float64)
    if  x <= 1.0    return( 1.0 )   end
    kappa  = sqrt( max(x - 1.0, 1.0e-14) )
    ## Stobbe's exact 1s cross section over the Kramers form; the common factor alpha a_o^2 cancels in the ratio.
    sigmaS = 2^9 * pi^2 / 3 * x^(-4) * exp(4 - 4*atan(kappa)/kappa) / (1 - exp(-2pi/kappa)) / exp(4.0)
    sigmaK = 64 * pi / (3 * sqrt(3.0)) * x^(-3)
    g1s    = sigmaS / sigmaK

    return( 1.0 + (g1s - 1.0) * n^(-2/3) )
end


"""
`Empirical.scaledBindingEnergy(Z::Float64, sh::Shell, conf::Configuration, data::PeriodicTable.AbstractEnergyData)`
    ... to provide the binding energy of an electron in shell sh of the given configuration conf; this energy scales all
        hydrogenic (Kramers) estimates below and, hence, fixes the photoionization threshold as well as the transition
        energies. The semi-empirical binding energies of the *neutral* atom are taken from the tabulation data whenever
        this data set covers Z and sh *and* conf is not itself a genuine few-electron ion (see below). If no tabulated
        value applies, a Slater-screened hydrogenic estimate bE = Z^(eff)^2 / (2 n^2) is applied instead and a warning
        is issued; Z^(eff) then follows from Slater's (1930) rules for the remaining electrons of conf. An
        energy::Float64 > 0. [a.u.] is returned.

        Note: A pure hydrogenic estimate is quite unreliable for valence shells; for neutral Ne, for instance, it places
              the 2p threshold near 85 eV, while the tabulated (true) binding energy is 21.6 eV. The tabulated values
              should therefore be preferred whenever they are available.

        Note: The tabulation is indexed by Z alone and, hence, always returns the *neutral*-atom value regardless of
              conf.NoElectrons. This is correct for the common "spectator-omitted shorthand" usage, where conf lists
              only the shells taking part in a transition of an otherwise near-neutral atom (e.g. a K-alpha analog
              Configuration("1s^1 2p^6") for Ne, with the 2s^2 spectator shell left out); it silently fails, however,
              for a conf that literally describes a stripped few-electron ion, such as He+'s Configuration("1s^1").
              A neutral atom with only 1 or 2 electrons is simply H or He, for which the tabulation already gives the
              correct (neutral) value; conf.NoElectrons < Z with conf.NoElectrons <= 2 therefore unambiguously signals
              a genuine few-electron ion rather than a spectator-omitted shorthand, and is routed directly to the
              Slater-screened hydrogenic estimate. This estimate is exact for H-like ions (no other electron to
              screen); for He-like ions, Slater's constant same-group screening (sigma = 0.30) systematically
              overestimates the binding energy, from about +31% at Z=3 (Li+) down to a few percent by Z ~ 20, as
              the ion becomes increasingly hydrogenic and the fixed screening constant matters proportionally
              less (checked against NIST/CRC data for the He isoelectronic sequence, Z = 3...18).
"""
function scaledBindingEnergy(Z::Float64, sh::Shell, conf::Configuration, data::PeriodicTable.AbstractEnergyData)
    fewElectronIon = conf.NoElectrons <= 2  &&  conf.NoElectrons < round(Int64, Z)

    bEnergy = 0.
    if  !fewElectronIon
        ## The tabulated data sets cover only a limited range of Z and shells; they raise an error or return 0. otherwise.
        try
            bEnergy = Empirical.bindingEnergy(round(Int64, Z), sh, data=data)
        catch
            bEnergy = 0.
        end
        if  bEnergy > 0.    return( bEnergy )    end
    end

    ## No tabulated value applies: screen the nuclear charge by all *other* electrons of conf with Slater's rules.
    shells = deepcopy(conf.shells)
    if      haskey(shells, sh)  &&  shells[sh] > 1    shells[sh] = shells[sh] - 1
    elseif  haskey(shells, sh)                        delete!(shells, sh)
    end
    rConf = Configuration(shells, conf.NoElectrons - 1)
    Zeff  = Semiempirical.estimateSlaterZeff(Z, rConf, Subshell(sh.n, -sh.l - 1))
    ## An electron of a neutral atom or positive ion sees at least a unit charge asymptotically.
    if  Zeff < 1.0      Zeff = 1.0      end
    if  fewElectronIon
        ## For a genuine few-electron ion this branch is exact for H-like ions; for He-like ions, Slater's 0.30
        ## same-group screening constant overestimates the binding energy by ~31% at Z=3, shrinking to a few
        ## percent by Z ~ 20 (checked against NIST/CRC data), not merely a fallback, so the warning says so.
        sa = "conf = $conf has only $(conf.NoElectrons) electron(s) at Z = $Z, i.e. it is a genuine few-electron ion, " *
             "for which the tabulated *neutral*-atom binding energy in $data would be wrong. A Slater-screened " *
             "hydrogenic estimate is used instead (exact for H-like ions; for He-like ions it overestimates the " *
             "binding energy by ~5-30%, decreasing with Z; this warning is shown at most 5 times)."
    else
        ## For (Rydberg) shells outside the tabulation this estimate is the designed behavior -- for a Rydberg electron
        ## it is even exact -- so the warning is limited; a summation over many capture channels would otherwise flood
        ## the output with one warning per shell.
        sa = "No tabulated binding energy for Z = $Z and $sh in $data; a Slater-screened hydrogenic estimate is used " *
             "(exact for Rydberg shells; this warning is shown at most 5 times)."
    end
    @warn sa maxlog=5

    return( Zeff^2 / (2 * sh.n^2) )
end


#################################################################################################################################
### Photoemission (PE) ##########################################################################################################


"""
`Empirical.photoemissionEinsteinA(iConf::Configuration, fConf::Configuration, approx::Empirical.ScaledHydrogenic;
                                  printout::Bool=false, data::PeriodicTable.AbstractEnergyData=PeriodicTable.XrayDataBooklet())`
    ... to estimate empirically the Einstein-A value for a transition from iConf -> fConf. The transition energy is
        obtained as the difference of the binding energies of the two shells involved, cf.
        Empirical.scaledBindingEnergy(). A named triple (multipole::EmMultipole=, energy::Float64=, rate::Float64=)
        is returned.
        Quantity: a spontaneous rate [1/s] -- an intrinsic property of the ion, independent of the plasma.

        Rates are estimated for E1 and M1 transitions only:
        + E1 (Delta-l = 1): from a hydrogenic scaling of the r-expectation values; good to a factor ~1.5
          (Ne K-alpha: 6.0e13 vs. the literature 8.8e13 1/s).
        + M1 (Delta-l = 0, always Delta-n != 0 here): such transitions are non-relativistically forbidden and
          proceed only via the relativistically induced M1 amplitude; the estimate scales the hydrogenic
          2s -> 1s anchor W(M1) = 2.496e-6 Z^10 1/s (Breit & Teller 1940; Johnson 1972) with the transition
          energy, (Delta-E/0.375 Hartree)^5, and is good to an order of magnitude only.
        + E2 and higher multipoles raise an error: the crude <r>-proxies misjudge their line strengths by
          orders of magnitude, and no silent zero is returned; use GivenEinsteinA(..) or the PhotoEmission
          module instead.
"""
function photoemissionEinsteinA(iConf::Configuration, fConf::Configuration, approx::Empirical.ScaledHydrogenic;
                                printout::Bool=false, data::PeriodicTable.AbstractEnergyData=PeriodicTable.XrayDataBooklet())
    multipole = Basics.E1;   tEnergy = rate = 0.;   Z = Defaults.getDefaults("nuclear: charge")
    iShell = fShell = Shell(0,0);    diff = 0
    
    # Determine initial and final shells involved in the transition; no multipole transition is assigned if the 
    # occupation of configurations differ by more than 1.
    wa = Basics.extractFromConfigurations(Basics.OccupationDifference(), iConf, fConf)
    if length(wa) > 2   multipole = missing   end 
    for  (k,v) in wa
        diff = diff + v
        if      v == -1    fShell = k
        elseif  v ==  1    iShell = k
        else    multipole = missing;   break
        end 
    end
    if  diff != 0   error("Incompatible initial and final configurations for a radiative transition.")   end

    # Determine the transition energy from the binding energies of the two shells involved; the electron moves from
    # iShell into the more tightly bound fShell, so that tEnergy = bE(fShell) - bE(iShell) > 0. For the Ne K-alpha
    # analogue this yields 870.2 - 21.6 = 848.6 eV, in agreement with the tabulated K-alpha energy.
    tEnergy = Empirical.scaledBindingEnergy(Z, fShell, fConf, data) - Empirical.scaledBindingEnergy(Z, iShell, iConf, data)

    # Determine the effective charges Zi and Zf felt by an electron in the initial and final shells, and apply the
    # standard hydrogenic formula <r> = (3 n^2 - l (l+1)) / (2 Z^(eff)) for the radial expectation values;
    # Just substract all inner-shell electrons (assume a charge 0.5 due to the other electrons in iShell or fShell).
    # Note that these radial Z^(eff) are kept distinct from the charge implied by the binding energy above; the two
    # differ considerably for valence shells, and the r-expectation values need the less-screened (radial) charge.
    ce      = 0;
    for  (shell, v) in iConf.shells
        if  shell.n < iShell.n  ||  (shell.n == iShell.n  &&  shell.l < iShell.l)   ce = ce + v   end
    end
    Zi      = Z - ce - 0.5;    rxi = (3*iShell.n^2 - iShell.l * (iShell.l+1)) / (2*Zi)
    ce      = 0;
    for  (shell, v) in fConf.shells
        if  shell.n < fShell.n  ||  (shell.n == fShell.n  &&  shell.l < fShell.l)   ce = ce + v   end
    end
    Zf      = Z - ce - 0.5;    rxf = (3*fShell.n^2 - fShell.l * (fShell.l+1)) / (2*Zf)

    # Determine the multipolarity of the transition
    if  multipole === missing    ||   diff != 0
    else
        if       abs(iShell.l - fShell.l) == 0        multipole = M1
        elseif   abs(iShell.l - fShell.l) == 1        multipole = E1
        elseif   abs(iShell.l - fShell.l) == 2        multipole = E2
        elseif   abs(iShell.l - fShell.l) == 3        multipole = E3
        else     error("Shell structure not supported.")
        end 
    end

    # Assign the energies and rates for the different multipoles
    if      multipole === missing        triple = (missing, 0., 0.)
    elseif  multipole == E1
        #  Determine an effective Z and apply the standard formula for the r-expectation values
        c      = Defaults.getDefaults("speed of light: c")
        rate   = 4/3 * tEnergy^3 / c^3 * abs( rxi * rxf )
        triple = (multipole, tEnergy, rate)
    elseif  multipole == M1
        #  Any M1 assigned here has Delta-n != 0 and Delta-l = 0 (equal shells never appear in a configuration
        #  difference); such transitions are forbidden non-relativistically -- the radial orbitals are orthogonal --
        #  and proceed only through the relativistically induced M1 amplitude. The hydrogenic anchor is the famous
        #  2s -> 1s decay, W(M1) = 1/972 m alpha (alpha Z)^10 ~ 2.496e-6 Z^10 1/s (Breit & Teller 1940; Johnson 1972);
        #  with Delta-E = 3/8 Z^2 Hartree, the Z^10 scaling becomes (Delta-E / 0.375 Hartree)^5. This energy-scaled
        #  form is applied as an order-of-magnitude estimate for all (induced) M1 transitions of this kind.
        rate   = 2.496e-6 * Defaults.convertUnits("time: from atomic to sec", 1.0) * (tEnergy / 0.375)^5
        triple = (multipole, tEnergy, rate)
    else
        error("No empirical rate estimate is supported for a $multipole transition ($iShell -> $fShell): the " *
              "hydrogenic <r>-expectation proxies of this ScaledHydrogenic approximation misjudge E2 (and higher) " *
              "line strengths by several orders of magnitude. Please provide the rate via GivenEinsteinA(..) or " *
              "compute it with the PhotoEmission module.")
    end
    
    # Report about this estimate
    if  printout
        unRate = Defaults.getDefaults("unit: rate");   unEnergy = Defaults.getDefaults("unit: energy")
        energyx = Defaults.convertUnits("energy: from atomic to " * unEnergy, tEnergy)
        ratex   = Defaults.convertUnits("rate: from atomic to "   * unRate,   rate)
        sa = "\n* Estimate empirically the Einstein-A value for a given transition i -> f with the " *
             "following assumptions/simplifications: " *
             "\n    + Use the tabulated binding energies ($data) for the energy and a hydrogenic scaling for the A-value. " *
             "\n    + iConf = $iConf  -->  fConf = $fConf " *
             "\n    + Extract transition energy from binding energies of shells $iShell -> $fShell " *
             "\n    + $multipole transition with energy [$unEnergy] = $energyx  " *
             "\n    + Rate [$unRate]                     = $ratex " *
             "\n    + Quantity: a spontaneous rate [$unRate] -- an intrinsic property of the ion, independent of the plasma. " * "\n"
        println(sa)
    end
    
    return( (multipole = multipole, energy = tEnergy, rate = rate) )
end

        
"""
`Empirical.photoemissionEinsteinA(iConf::Configuration, fConf::Configuration, approx::Empirical.UsingJAC;
                                  printout::Bool=false)`  
    ... to estimate empirically the Einstein-A value for a transition from iConf -> fConf by using simple 
        JAC computations of transition and Einstein rates. A named triple
        (multipole::EmMultipole=, energy::Float64=, rate::Float64=) is returned.
        Quantity: a spontaneous rate [1/s] -- an intrinsic property of the ion, independent of the plasma.
"""
function photoemissionEinsteinA(iConf::Configuration, fConf::Configuration, approx::Empirical.UsingJAC;
                                printout::Bool=false) 
    multipole = Basics.E1;   tEnergy = rate = 0.;   Z = Defaults.getDefaults("nuclear: charge")
    iShell = fShell = Shell(0,0);    diff = 0
    
    # Determine initial and final shells involved in the transition; no multipole transition is assigned if the 
    # occupation of configurations differ by more than 1.
    wa = Basics.extractFromConfigurations(Basics.OccupationDifference(), iConf, fConf)
    if length(wa) > 2   multipole = missing   end 
    for  (k,v) in wa
        diff = diff + v
        if      v == -1    fShell = k
        elseif  v ==  1    iShell = k
        else    multipole = missing;   break
        end 
    end
    if  diff != 0   error("Incompatible initial and final configurations for a radiative transition.")   end 
    
    # Determine the multipolarity of the transition
    if  multipole === missing    ||   diff != 0
    else
        if       abs(iShell.l - fShell.l) == 0        multipole = M1
        elseif   abs(iShell.l - fShell.l) == 1        multipole = E1
        elseif   abs(iShell.l - fShell.l) == 2        multipole = E2
        elseif   abs(iShell.l - fShell.l) == 3        multipole = E3
        else     error("Shell structure not supported.")
        end
    end

    # Generate mean-field orbitals in order to extract the transition energies and amplitudes
    grid        = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, rbox = 10.0) 
    mfSettings  = AtomicState.MeanFieldSettings(Basics.DFSField(1.0))
    meanField   = Representation("Internal", Nuclear.Model(Z), grid, [iConf], MeanFieldBasis(mfSettings) )
    mfrep       = generate(meanField; output=true)
    iOrbitals   = mfrep["mean-field basis"].orbitals

    meanField   = Representation("Internal", Nuclear.Model(Z), grid, [fConf], MeanFieldBasis(mfSettings) )
    mfrep       = generate(meanField; output=true)
    fOrbitals   = mfrep["mean-field basis"].orbitals
        

    # Assign the energies and rates for the different multipoles
    c       = Defaults.getDefaults("speed of light: c")
    iSubsh  = Subshell(iShell.n, -iShell.l -1)
    fSubsh  = Subshell(fShell.n, -fShell.l -1)
    tEnergy = iOrbitals[iSubsh].energy - fOrbitals[fSubsh].energy
    
    if      multipole == E1
        amp  = InteractionStrength.MabEmissionJohnsony(Basics.E1, Basics.Babushkin, tEnergy, fOrbitals[fSubsh],  
                                                       iOrbitals[iSubsh], grid)
        wa   = 8.0pi * Defaults.getDefaults("alpha") * tEnergy / 
               Basics.extractFromConfiguration(Basics.Multiplicity(), iConf) *
               Basics.extractFromConfiguration(Basics.Multiplicity(), fConf) 
        rate = wa * abs(amp)^2
        triple = (multipole, tEnergy, rate)
    elseif  multipole == M1
        amp  = InteractionStrength.MabEmissionJohnsony(Basics.M1, Basics.Magnetic, tEnergy, fOrbitals[fSubsh],  
                                                       iOrbitals[iSubsh], grid)
        wa   = 8.0pi * Defaults.getDefaults("alpha") * tEnergy / 
               Basics.extractFromConfiguration(Basics.Multiplicity(), iConf) *
               Basics.extractFromConfiguration(Basics.Multiplicity(), fConf) 
        rate = wa * abs(amp)^2
        triple = (multipole, tEnergy, rate)
    else                                 
        triple = (multipole, 0., 0.) 
    end
    
    # Report about this estimate
    if  printout
        unRate = Defaults.getDefaults("unit: rate");   unEnergy = Defaults.getDefaults("unit: energy")
        energyx = Defaults.convertUnits("energy: from atomic to " * unEnergy, tEnergy)
        ratex   = Defaults.convertUnits("rate: from atomic to "   * unRate,   rate)
        sa = "\n* Estimate empirically the Einstein-A value for a given transition i -> f with the " *
             "following assumptions/simplifications: " *
             "\n    + Use a simple JAC computations for the energy and A-value. " * 
             "\n    + iConf = $iConf  -->  fConf = $fConf " * 
             "\n    + Extract transition energy from binding energies of subshells $iSubsh -> $fSubsh " * 
             "\n    + $multipole transition with energy [$unEnergy] = $energyx  " * 
             "\n    + Rate [$unRate]                     = $ratex " *
             "\n    + Quantity: a spontaneous rate [$unRate] -- an intrinsic property of the ion, independent of the plasma. " * "\n"
        println(sa)
    end
    
    return( (multipole = multipole, energy = tEnergy, rate = rate) )
end


"""
`Empirical.photoemissionEinsteinA(approx::Empirical.GivenEinsteinA; printout::Bool=false)`  
    ... to simply return the Einstein-A value for the given transition. A named
        (multipole::EmMultipole=, energy::Float64=, rate::Float64=) is returned.
        Quantity: a spontaneous rate [1/s] -- an intrinsic property of the ion, independent of the plasma.
"""
function photoemissionEinsteinA(approx::Empirical.GivenEinsteinA; printout::Bool=false) 
    
    # Report about this estimate
    if  printout
        unRate = Defaults.getDefaults("unit: rate");   unEnergy = Defaults.getDefaults("unit: energy")
        energyx = Defaults.convertUnits("energy: from atomic to " * unEnergy, approx.energy)
        ratex   = Defaults.convertUnits("rate: from atomic to "   * unRate,   approx.rate)
        sa = "\n* User-given Einstein-A value for a given transition i -> f: "   *
             "\n    + $(approx.multipole) transition with energy = $energyx  " * unEnergy *
             "\n    + rate = $ratex " * unRate *
             "\n    + Quantity: a spontaneous rate [" * unRate * "] -- an intrinsic property of the ion, independent of the plasma. " * "\n"
        println(sa)
    end
    
    return( (multipole = approx.multipole, energy = approx.energy, rate = approx.rate) )
end


#################################################################################################################################
### Photoexcitation (PX) ########################################################################################################

"""
`Empirical.photoexcitationPlasmaRatePerIon(dist::Distribution.AbstractPhotonDistribution,
                                           iConf::Configuration, fConf::Configuration; 
                                           approx::Empirical.AbstractEmpiricalApproximation=UsingJAC(), printout::Bool=false)` 
                                     
    ... to estimate the photoexcitation plasma rate per ion R^(PX: per ion) (T; i -> f) for a transition from 
        iConf -> fConf by applying a given photon distribution dist. A rate::Float64 is returned.
        Quantity: a rate per ion [1/s] in the given field -- multiply by the ion number density n_ion [1/cm^3] to obtain
            the volumetric rate [1/(cm^3 s)].
"""
function photoexcitationPlasmaRatePerIon(dist::Distribution.AbstractPhotonDistribution,
                                         iConf::Configuration, fConf::Configuration; 
                                         approx::Empirical.AbstractEmpiricalApproximation=UsingJAC(), printout::Bool=false)
    EinsteinA = Empirical.photoemissionEinsteinA(iConf, fConf, approx, printout=false)
    c         = Defaults.getDefaults("speed of light: c")
    pnDensity = Distribution.photonNumberDensity(dist, EinsteinA.energy)
    occNumber = pnDensity * pi^2 * c^3 / EinsteinA.energy^2    ## n̄ = 1/(exp(ω/T)-1); photonNumberDensity returns n(ω) = ω²/(π²c³) × n̄
    rate      = Basics.extractFromConfiguration(Basics.Multiplicity(), iConf) /
                Basics.extractFromConfiguration(Basics.Multiplicity(), fConf) * EinsteinA.rate * occNumber

    # Report about this estimate
    if  printout
        Tx      = Defaults.convertUnits("temperature: from atomic to Kelvin", dist.T)
        unRate  = Defaults.getDefaults("unit: rate");   unEnergy = Defaults.getDefaults("unit: energy")
        energyx = Defaults.convertUnits("energy: from atomic to " * unEnergy, EinsteinA.energy)
        ratex   = Defaults.convertUnits("rate: from atomic to "   * unRate,   rate)
        sa = "\n* Estimate empirically the photoexcitation plasma rate R^(PX: per ion) (T; i -> f) for a given transition " *
             "i -> f with the following assumptions/simplifications: " *
             "\n    + Photon field: $dist,  i.e. T [K] = $(Tx). " * 
             "\n    + Einstein-A values are determined in the $(SubString(string(approx), 22)) approximation " * 
             "\n    + iConf = $iConf  -->  fConf = $fConf " * 
             "\n    + Transition energy [$unEnergy]                                = $energyx " *
             "\n    + Plasma rate per ion R^(PX: per ion) (T; i -> f) [$unRate] = $ratex " *
             "\n    + Quantity: a rate per ion [$unRate] in this field -- multiply by the ion number density n_ion [1/cm^3] for the volumetric rate [1/(cm^3 s)]. " * "\n"
        println(sa)
    end
    
    return( rate )
end



#################################################################################################################################
### Photodeexcitation (PD) ######################################################################################################

"""
`Empirical.photodeexcitationPlasmaRatePerIon(dist::Distribution.AbstractPhotonDistribution,
                                             iConf::Configuration, fConf::Configuration; 
                                             approx::Empirical.AbstractEmpiricalApproximation=UsingJAC(), printout::Bool=false)` 
                                     
    ... to estimate the total photodeexcitation plasma rate R^(PD: total, per ion) (T; i -> f) for a transition from 
        iConf -> fConf by applying a given photon distribution dist. A rate::Float64 is returned.
        Quantity: a rate per ion [1/s] in the given field -- multiply by the ion number density n_ion [1/cm^3] to obtain
            the volumetric rate [1/(cm^3 s)].
"""
function photodeexcitationPlasmaRatePerIon(dist::Distribution.AbstractPhotonDistribution,
                                           iConf::Configuration, fConf::Configuration; 
                                           approx::Empirical.AbstractEmpiricalApproximation=UsingJAC(), printout::Bool=false)
    EinsteinA = Empirical.photoemissionEinsteinA(iConf, fConf, approx, printout=false)
    c         = Defaults.getDefaults("speed of light: c")
    pnDensity = Distribution.photonNumberDensity(dist, EinsteinA.energy)
    occNumber = pnDensity * pi^2 * c^3 / EinsteinA.energy^2    ## n̄ = 1/(exp(ω/T)-1)
    rate      = EinsteinA.rate * (1.0 + occNumber)

    # Report about this estimate
    if  printout
        Tx      = Defaults.convertUnits("temperature: from atomic to Kelvin", dist.T)
        unRate  = Defaults.getDefaults("unit: rate");   unEnergy = Defaults.getDefaults("unit: energy")
        energyx = Defaults.convertUnits("energy: from atomic to " * unEnergy, EinsteinA.energy)
        ratex   = Defaults.convertUnits("rate: from atomic to "   * unRate,   rate)
        sa = "\n* Estimate empirically the (total) photodeexcitation plasma rate R^(PD: total, per ion) (T; i -> f) " *
             "for a given transition i -> f" *
             "\n  with the following assumptions/simplifications: " *
             "\n    + Photon field: $dist,  i.e. T [K] = $(Tx). " * 
             "\n    + Einstein-A values are determined in the $(SubString(string(approx), 22)) approximation " * 
             "\n    + iConf = $iConf  -->  fConf = $fConf " * 
             "\n    + Transition energy [$unEnergy]                               = $energyx " *
             "\n    + Plasma rate R^(PD: total, per ion) (T; i -> f) [$unRate] = $ratex " *
             "\n    + Quantity: a rate per ion [$unRate] in this field -- multiply by the ion number density n_ion [1/cm^3] for the volumetric rate [1/(cm^3 s)]. " * "\n"
        println(sa)
    end
    
    return( rate )
end



#################################################################################################################################
### Photoionization (PI) ########################################################################################################


"""
`Empirical.photoionizationCrossSection(omega::Float64, iConf::Configuration, fConf::Configuration; 
                                       approx::Empirical.AbstractEmpiricalApproximation=ScaledHydrogenic(), printout::Bool=false)` 
                                     
    ... to estimate empirically the PI cross section for a transition from iConf -> fConf by applying some simple approximation 
        as determined by approx. For printout=true, basic information are printed about the input parameters, approximation 
        as well as the results in user-defined units. A cs::Float64 [in a.u.] is returned.
        Quantity: a cross section [a.u.] -- a property of the ion alone, independent of the plasma; fold it with the
            photon/electron field, or pass the configurations to a *PlasmaRatePerIon / *PlasmaAlpha function, to obtain a rate.
"""
function photoionizationCrossSection(omega::Float64, iConf::Configuration, fConf::Configuration; 
                                     approx::Empirical.AbstractEmpiricalApproximation=ScaledHydrogenic(), printout::Bool=false)
    
    cs = Empirical.photoionizationCrossSection([omega], iConf, fConf, approx, printout=printout) 
    
    return( cs[1] )
end


"""
`Empirical.photoionizationCrossSection(omegas::Array{Float64,1}, iConf::Configuration, fConf::Configuration,
                                       approx::Empirical.ScaledHydrogenic; printout::Bool=false,
                                       data::PeriodicTable.AbstractEnergyData=PeriodicTable.XrayDataBooklet())`
    ... to estimate empirically the photoionization cross sections for a transition from iConf -> fConf by using the binding
        energy (ionization potential) of the ionized shell, cf. Empirical.scaledBindingEnergy(), and Kramer's (1923)
        empirical formula. A css::Array{Float64,1} [a.u.] is returned.
        Quantity: a cross section [a.u.] -- a property of the ion alone, independent of the plasma; fold it with the
            photon/electron field, or pass the configurations to a *PlasmaRatePerIon / *PlasmaAlpha function, to obtain a rate.
"""
function photoionizationCrossSection(omegas::Array{Float64,1}, iConf::Configuration, fConf::Configuration,
                                     approx::Empirical.ScaledHydrogenic; printout::Bool=false,
                                     data::PeriodicTable.AbstractEnergyData=PeriodicTable.XrayDataBooklet())
    Z = Defaults.getDefaults("nuclear: charge");    iShell = Shell(0,0);    diff = 0

    # Determine the initial shell and its binding (threshold) energy.
    wa = Basics.extractFromConfigurations(Basics.OccupationDifference(), iConf, fConf)
    if length(wa) > 1   error("Incompatible initial and final configurations for a photoionization cross section.")   end

    for  (k,v) in wa
        diff = diff + v
        if      v == 1     iShell = k    end
    end
    if  diff != 1   error("Incompatible initial and final configurations for a photoionization cross section.")   end

    # Scale the hydrogenic formula by the binding energy (photoionization threshold) of the ionized shell.
    bEnergy = Empirical.scaledBindingEnergy(Z, iShell, iConf, data)

    # Kramers' (1923) bound-free cross section, sigma = 32 pi alpha N_e / (3 sqrt(3) n) * bE^2 / omega^3 * g^(bf),
    # where N_e is the number of equivalent electrons in iShell: the photon may eject any one of them and, in an
    # independent-particle picture, their cross sections simply add. The bound-free Gaunt factor g^(bf), cf.
    # Empirical.boundFreeGauntFactor(), corrects the semiclassical form: for H 1s it turns the Kramers 7.91 Mb at
    # threshold into the exact 6.30 Mb and steepens the tail from omega^(-3) to the exact omega^(-7/2).
    nEquiv  = iConf.shells[iShell]
    factor  = 32 * pi * Defaults.getDefaults("alpha") * nEquiv / (3 * sqrt(3) * iShell.n)

    # Compute all cross sections in the ScaledHydrogenic approximation
    css = Float64[]
    for  omega in omegas
        if  omega > bEnergy
            push!(css, factor * bEnergy^2 / omega^3 * Empirical.boundFreeGauntFactor(iShell.n, omega/bEnergy) )
        else
            push!(css, 0.)
        end
    end

    # Report about these estimates
    if  printout
        unCs    = Defaults.getDefaults("unit: cross section");   unEnergy = Defaults.getDefaults("unit: energy")
        energyx = Defaults.convertUnits("energy: from atomic to " * unEnergy, bEnergy)
        omegasx = Float64[];   cssx = Float64[];   
        for omega in omegas  push!(omegasx, Defaults.convertUnits("energy: from atomic to " * unEnergy, omega))   end
        for cs    in css     push!(cssx,    Defaults.convertUnits("cross section: from atomic to " * unCs, cs))   end
        sa = "\n* Estimate empirically the photoionization cross section for a given transition i -> f with the " *
             "following assumptions/simplifications: " *
             "\n    + Use Kramers' (1923) formula, scaled by the binding energy of the ionized shell. " *
             "\n    + Binding energy (photoionization threshold) taken from $data " *
             "\n    + iConf = $iConf  -->  fConf = $fConf " *
             "\n    + Binding energy of $iShell   = $energyx  " * unEnergy *
             "\n    + Omegas [$unEnergy]            = $omegasx " *
             "\n    + Cross sections [$unCs]  = $cssx     " *
             "\n    + Quantity: cross section [$unCs] -- a property of the ion alone, independent of the plasma; fold it with the photon/electron field, or pass iConf/fConf to a *PlasmaRatePerIon / *PlasmaAlpha function, to obtain a rate. " * "\n"
        println(sa)
    end
    
    return( css )
end


"""
`Empirical.photoionizationCrossSection(omegas::Array{Float64,1}, iConf::Configuration, fConf::Configuration,
                                       approx::Empirical.UsingJAC; printout::Bool=false)`
    ... to estimate the photoionization cross section for a transition from iConf -> fConf from mean-field (DFS)
        orbitals and single-electron E1 amplitudes,
            sigma^(PI) (omega) = N_shell/(4l+2) * 8 pi^3/(alpha omega) *
                                 sum_(subshells j = l+-1/2) sum_(E1 channels kappa_c) |<eps kappa_c || O^(E1) || n kappa>|^2,
        with energy-normalized continuum orbitals in the local (nuclear + DFS) potential; the prefactor
        8 pi^3/(alpha omega) belongs to the Johnson (2007) convention of the reduced emission matrix elements and is
        consistent with the (well-tested) bound-bound Einstein rate A = 8 pi alpha omega/(2J_i+1) |<f||O||i>|^2 by
        smearing a line into the continuum. A css::Array{Float64,1} [a.u.] is returned.
        Quantity: a cross section [a.u.] -- a property of the ion alone, independent of the plasma; fold it with the
            photon/electron field, or pass the configurations to a *PlasmaRatePerIon / *PlasmaAlpha function, to obtain a rate.

        Note: Validated against experiment: He 1s^2 gives 6.3/4.9/1.7/0.32 Mb at 26/30/50/100 eV against the measured
              ~7.4/5.5/2/0.3 Mb (within ~15%); Ne 2p^6 reproduces shape and magnitude to ~30-50%. The *thresholds* are
              the mean-field (DFS) orbital energies and lie below the true ionization potentials (Ne 2p: 12.0 eV
              instead of 21.6 eV) -- a property of this approximation that has to be kept in mind near threshold.
              The continuum requires a radial grid with a linear tail; the exponential grid formerly used here
              produced spurious box resonances (92 Mb at 100 eV for Ne 2p) that invalidated all earlier results.
              The amplitudes further depend on the globally selected continuum method; the validation above refers
              to JAC's defaults, i.e. setDefaults("method: continuum, Galerkin") and
              setDefaults("method: normalization, Alok").
"""
function photoionizationCrossSection(omegas::Array{Float64,1}, iConf::Configuration, fConf::Configuration,
                                     approx::Empirical.UsingJAC; printout::Bool=false)
    Z = Defaults.getDefaults("nuclear: charge");    iShell = Shell(0,0);    diff = 0

    # Determine the ionized shell.
    wa = Basics.extractFromConfigurations(Basics.OccupationDifference(), iConf, fConf)
    if length(wa) > 1   error("Incompatible initial and final configurations for a photoionization cross section.")   end
    for  (k,v) in wa
        diff = diff + v
        if      v == 1     iShell = k    end
    end
    if  diff != 1   error("Incompatible initial and final configurations for a photoionization cross section.")   end

    # The E1 continuum channels kappa_c of an ionized subshell (l, 2j): parity change l_c = l +- 1 and |j_c - j| <= 1.
    function e1Channels(li::Int64, ji2::Int64)
        chans = Int64[]
        for  lc in (li-1, li+1)
            lc < 0    &&  continue
            for  jc2 in (2lc-1, 2lc+1)
                jc2 < 1              &&  continue
                abs(jc2 - ji2) > 2   &&  continue
                push!(chans, jc2 == 2lc+1 ? -lc-1 : lc)
            end
        end
        return( chans )
    end

    # Generate mean-field orbitals in order to extract the threshold energies and amplitudes; the radial grid carries
    # a linear tail (hp), which the oscillating continuum orbitals require.
    grid        = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.5e-2, rbox = 20.0)
    mfSettings  = AtomicState.MeanFieldSettings(Basics.DFSField(1.0))
    meanField   = Representation("Internal", Nuclear.Model(Z), grid, [iConf], MeanFieldBasis(mfSettings) )
    mfrep       = generate(meanField; output=true)
    iBasis      = mfrep["mean-field basis"]

    contSettings = Continuum.Settings(false, size(grid.r,1) - 11);    nm = Nuclear.Model(Z)
    nucPot   = Nuclear.nuclearPotential(nm, grid)
    dfsPot   = Basics.computePotential(Basics.DFSField(1.0), grid, iBasis)
    localPot = Basics.add(nucPot, dfsPot)

    li        = iShell.l;    nShell = iConf.shells[iShell];    alfa = Defaults.getDefaults("alpha")
    subshells = li == 0  ?  [Subshell(iShell.n, -1)]  :  [Subshell(iShell.n, li), Subshell(iShell.n, -li-1)]
    bEnergy   = minimum( -iBasis.orbitals[s].energy  for s in subshells if haskey(iBasis.orbitals, s) )

    css = Float64[]
    for  omega in omegas
        sigma = 0.
        for  subsh in subshells
            haskey(iBasis.orbitals, subsh)   ||  continue
            eps = omega + iBasis.orbitals[subsh].energy         # orbital energy is negative
            eps <= 0.   &&  continue
            ji2 = Basics.twice(Basics.subshell_j(subsh))
            for  kapc in e1Channels(li, ji2)
                cOrbital, phase, norm = Continuum.generateOrbitalLocalPotential(eps, Subshell(101, kapc), localPot, contSettings)
                amp   = InteractionStrength.MabEmissionJohnsony(Basics.E1, Basics.Babushkin, omega, cOrbital,
                                                                iBasis.orbitals[subsh], grid)
                sigma = sigma + nShell / (2*(2li+1)) * 8 * pi^3 / (alfa * omega) * abs(amp)^2
            end
        end
        push!(css, sigma)
    end

    # Report about these estimates
    if  printout
        unCs    = Defaults.getDefaults("unit: cross section");   unEnergy = Defaults.getDefaults("unit: energy")
        energyx = Defaults.convertUnits("energy: from atomic to " * unEnergy, bEnergy)
        omegasx = Float64[];   cssx = Float64[];   
        for omega in omegas  push!(omegasx, Defaults.convertUnits("energy: from atomic to " * unEnergy, omega))   end
        for cs    in css     push!(cssx,    Defaults.convertUnits("cross section: from atomic to " * unCs, cs))   end
        sa = "\n* Estimate empirically the photoionization cross section for a given transition i -> f with the " *
             "following assumptions/simplifications: " *
             "\n    + Use a JAC mean-field (DFS) approach: one-electron E1 amplitudes with energy-normalized continuum " *
             "orbitals, summed over the subshells j = l +- 1/2 and their E1 channels. " *
             "\n    + The threshold energies are the mean-field orbital energies; they lie below the true ionization " *
             "potentials. " *
             "\n    + iConf = $iConf  -->  fConf = $fConf " *
             "\n    + Mean-field threshold of $iShell   = $energyx  " * unEnergy *
             "\n    + Omegas [$unEnergy]            = $omegasx " *
             "\n    + Cross sections [$unCs]  = $cssx     " *
             "\n    + Quantity: cross section [$unCs] -- a property of the ion alone, independent of the plasma; fold it with the photon/electron field, or pass iConf/fConf to a *PlasmaRatePerIon / *PlasmaAlpha function, to obtain a rate. " * "\n"
        println(sa)
    end
    
    return( css )
end


"""
`Empirical.photoionizationPlasmaRatePerIon(dist::Distribution.AbstractPhotonDistribution,
                                      iConf::Configuration, fConf::Configuration;
                                      approx::Empirical.AbstractEmpiricalApproximation=ScaledHydrogenic(), printout::Bool=false,
                                      data::PeriodicTable.AbstractEnergyData=PeriodicTable.XrayDataBooklet())`
    ... to estimate empirically the photoionization plasma rate per ion R^(PI: per ion) (T; i -> f) for a transition from
        iConf -> fConf and a photon field dist at the radiation temperature T. A rate::Float64 [a.u.] is returned.
        Quantity: a rate per ion [1/s] in the given field -- multiply by the ion number density n_ion [1/cm^3] to obtain
            the volumetric rate [1/(cm^3 s)].

        Note: In contrast to the photorecombination coefficient alpha^(PR) [cm^3/s], this quantity is a *rate* [1/s] and
              not a rate coefficient. The convolution R^(PI: per ion) = int d(omega) n(omega; T) c sigma^(PI)(omega)
              already contains the photon number density of the radiation field and, hence, needs no further
              multiplication by a density; it is the number of photoionization events per ion and per unit time.
              The free-electron density, in contrast, does not enter the PI process at all.
"""
function photoionizationPlasmaRatePerIon(dist::Distribution.AbstractPhotonDistribution,
                                    iConf::Configuration, fConf::Configuration;
                                    approx::Empirical.AbstractEmpiricalApproximation=ScaledHydrogenic(), printout::Bool=false,
                                    data::PeriodicTable.AbstractEnergyData=PeriodicTable.XrayDataBooklet())
    if  typeof(approx) != Empirical.ScaledHydrogenic
        error("The photoionization plasma rate is presently supported only for approx = ScaledHydrogenic(); the " *
              "threshold energy, which is needed to place the integration mesh, is not available for $approx.")
    end
    Z = Defaults.getDefaults("nuclear: charge");   alpha = 0.;   iShell = Shell(0,0);   diff = 0

    # Determine the ionized shell and its binding energy; the PI cross section vanishes below this threshold.
    wa = Basics.extractFromConfigurations(Basics.OccupationDifference(), iConf, fConf)
    if length(wa) > 1   error("Incompatible initial and final configurations for a photoionization plasma rate.")   end
    for  (k,v) in wa
        diff = diff + v
        if      v == 1     iShell = k    end
    end
    if  diff != 1   error("Incompatible initial and final configurations for a photoionization plasma rate.")   end
    bEnergy = Empirical.scaledBindingEnergy(Z, iShell, iConf, data)
    c       = Defaults.getDefaults("speed of light: c")

    # Perform a Gauss-Legendre integration of  R^(PI: per ion) = int d(omega) n(omega; T) c sigma^(PI)(omega).
    # The mesh starts at the threshold: sigma^(PI) jumps from zero to its maximum at omega = bEnergy, and a mesh that
    # straddles this step converges badly (> 50% error for a mesh 0 ... 30 T). It ends at bEnergy + 30 T, since all
    # supported photon fields fall off exponentially above their temperature. This quadrature is converged to < 0.2%.
    gridGL    = Radial.GridGL(Radial.GridGaussLegendreFinite(), bEnergy, bEnergy + 30*dist.T, 96; printout=true)
    css       = Empirical.photoionizationCrossSection(gridGL.t, iConf, fConf, approx; printout=false, data=data)
    for  n = 1:gridGL.nt
        alpha = alpha + css[n] * c * Distribution.photonNumberDensity(dist, gridGL.t[n]) * gridGL.wt[n]
    end

    # Report about this estimate
    if  printout
        unRate = Defaults.getDefaults("unit: rate");   unEnergy = Defaults.getDefaults("unit: energy")
        Tx     = Defaults.convertUnits("temperature: from atomic to Kelvin", dist.T)
        alphax = Defaults.convertUnits("rate: from atomic to " * unRate, alpha)
        bEx    = Defaults.convertUnits("energy: from atomic to " * unEnergy, bEnergy)
        sa = "\n* Estimate empirically the photoionization plasma rate per ion for a given transition i -> f with the " *
             "following assumptions/simplifications: " *
             "\n    + Photon field: $dist,  i.e. T [K] = $(Tx). " *
             "\n    + PI cross sections are generated in the $(SubString(string(approx), 22)) approximation. " *
             "\n    + iConf = $iConf  -->  fConf = $fConf " *
             "\n    + Threshold energy of $iShell [$unEnergy] = $bEx " *
             "\n    + Plasma rate per ion R^(PI: per ion) (T; i -> f) [$unRate] = $alphax " *
             "\n    + Quantity: a rate per ion [$unRate] in this field -- multiply by the ion number density n_ion [1/cm^3] for the volumetric rate [1/(cm^3 s)]. " * "\n"
        println(sa)
    end

    return( alpha )
end
    

#################################################################################################################################
### Photorecombination (PR) #####################################################################################################


"""
`Empirical.photorecombinationCrossSection(energy::Float64, iConf::Configuration, fConf::Configuration; 
                                          approx::Empirical.AbstractEmpiricalApproximation=ScaledHydrogenic(), printout::Bool=false)` 
                                     
    ... to estimate empirically the PR cross section for a transition from iConf -> fConf by applying some simple approximation 
        as determined by approx. For printout=true, basic information are printed about the input parameters, approximation 
        as well as the results in user-defined units. A cs::Float64 [in a.u.] is returned.
        Quantity: a cross section [a.u.] -- a property of the ion alone, independent of the plasma; fold it with the
            photon/electron field, or pass the configurations to a *PlasmaRatePerIon / *PlasmaAlpha function, to obtain a rate.
"""
function photorecombinationCrossSection(energy::Float64, iConf::Configuration, fConf::Configuration; 
                                        approx::Empirical.AbstractEmpiricalApproximation=ScaledHydrogenic(), printout::Bool=false)
    
    cs = Empirical.photorecombinationCrossSection([energy], iConf, fConf, approx, printout=printout) 
    
    return( cs )
end


"""
`Empirical.photorecombinationCrossSection(energies::Array{Float64,1}, iConf::Configuration, fConf::Configuration,
                                          approx::Empirical.ScaledHydrogenic; printout::Bool=false,
                                          data::PeriodicTable.AbstractEnergyData=PeriodicTable.XrayDataBooklet())`
    ... to estimate empirically the (spontaneous) photorecombination cross section for a transition from iConf -> fConf
        by using the Einstein-Milne relation and the binding energy (ionization potential) of the captured electron,
        cf. Empirical.scaledBindingEnergy(), together with Kramer's (1923) empirical formula.
        A css::Array{Float64,1} [a.u.] is returned.
        Quantity: a cross section [a.u.] -- a property of the ion alone, independent of the plasma; fold it with the
            photon/electron field, or pass the configurations to a *PlasmaRatePerIon / *PlasmaAlpha function, to obtain a rate.
"""
function photorecombinationCrossSection(energies::Array{Float64,1}, iConf::Configuration, fConf::Configuration,
                                        approx::Empirical.ScaledHydrogenic; printout::Bool=false,
                                        data::PeriodicTable.AbstractEnergyData=PeriodicTable.XrayDataBooklet())
    Z = Defaults.getDefaults("nuclear: charge");    fShell = Shell(0,0);    diff = 0

    # Determine the final shell into which the electron is captured.
    wa = Basics.extractFromConfigurations(Basics.OccupationDifference(), iConf, fConf)
    if length(wa) > 1   error("Incompatible initial and final configurations for a photorecombination cross section.")   end

    for  (k,v) in wa
        diff = diff + v
        if     v == -1     fShell = k    end
    end
    if  diff != -1   error("Incompatible initial and final configurations for a photorecombination cross section.")   end

    # Binding energy of the captured electron in fShell of the recombined ion fConf; this fixes the photon energies
    # omega = eps + bE and must be the same threshold as seen by the inverse (photoionization) process.
    bEnergy = Empirical.scaledBindingEnergy(Z, fShell, fConf, data)
    omegas  = energies .+ bEnergy

    # Compute the PI cross sections of the inverse process directly from the same bEnergy; cf. Kramers' formula and
    # the bound-free Gaunt factor in Empirical.photoionizationCrossSection(..., ScaledHydrogenic). The inverse process
    # ionizes fShell of fConf, so that N_e is here the occupation of fShell in the *recombined* ion fConf.
    nEquiv  = fConf.shells[fShell]
    factor  = 32 * pi * Defaults.getDefaults("alpha") * nEquiv / (3 * sqrt(3) * fShell.n)
    piCss   = [factor * bEnergy^2 / omega^3 * Empirical.boundFreeGauntFactor(fShell.n, omega/bEnergy)   for omega in omegas]

    # Apply the Einstein-Milne relation:  sigma^(PR) = omega^2 / (2 eps c^2) * (g_f/g_i) * sigma^(PI)
    c       = Defaults.getDefaults("speed of light: c")
    gf_gi   = Basics.extractFromConfiguration(Basics.Multiplicity(), fConf) /
              Basics.extractFromConfiguration(Basics.Multiplicity(), iConf)
    prCss   = [gf_gi * omegas[im]^2 / (2 * energies[im] * c^2) * piCss[im] for im in eachindex(omegas)]
    
    # Report about this estimate
    if  printout
        unCs    = Defaults.getDefaults("unit: cross section");   unEnergy = Defaults.getDefaults("unit: energy")
        energyx = Defaults.convertUnits("energy: from atomic to " * unEnergy, bEnergy)
        omegasx = Float64[];   cssx = Float64[];   
        for omega in omegas  push!(omegasx, Defaults.convertUnits("energy: from atomic to " * unEnergy, omega))   end
        for cs    in prCss   push!(cssx,    Defaults.convertUnits("cross section: from atomic to " * unCs, cs))   end
        sa = "\n* Estimate empirically the photorecombination cross section for a given transition i -> f with the " *
             "following assumptions/simplifications: " *
             "\n    + Use the Einstein-Milne relation to obtain PR cross sections from PI cross sections from the " * 
             "ground level of $fConf " *
             "\n    + PI cross sections are determined in the $(SubString(string(approx), 22)) approximation " * 
             "\n    + iConf = $iConf  -->  fConf = $fConf " * 
             "\n    + Binding energy of $fShell  [$unEnergy]  = $energyx  " *
             "\n    + Omegas [$unEnergy]                 = $omegasx " *
             "\n    + PR Cross sections [$unCs]    = $cssx     " *
             "\n    + Quantity: cross section [$unCs] -- a property of the ion alone, independent of the plasma; fold it with the photon/electron field, or pass iConf/fConf to a *PlasmaRatePerIon / *PlasmaAlpha function, to obtain a rate. " * "\n"
        println(sa)
    end
    
    return( prCss )
end


"""
`Empirical.photorecombinationCrossSection(energies::Array{Float64,1}, iConf::Configuration, fConf::Configuration,
                                          approx::Empirical.UsingJAC; printout::Bool=false)`
    ... to estimate the (spontaneous) photorecombination cross section for a transition from iConf -> fConf by applying
        the Einstein-Milne relation to the mean-field (UsingJAC) PI cross sections of the inverse process; the binding
        energy of the captured electron is taken from the mean-field orbital of fConf and is, therefore, consistent
        with the thresholds of those PI cross sections. A css::Array{Float64,1} [a.u.] is returned.
        Quantity: a cross section [a.u.] -- a property of the ion alone, independent of the plasma; fold it with the
            photon/electron field, or pass the configurations to a *PlasmaRatePerIon / *PlasmaAlpha function, to obtain a rate.

        Note: cf. the validation notes of Empirical.photoionizationCrossSection(..., UsingJAC); in particular, the
              mean-field thresholds lie below the true ionization potentials.
"""
function photorecombinationCrossSection(energies::Array{Float64,1}, iConf::Configuration, fConf::Configuration,
                                        approx::Empirical.UsingJAC; printout::Bool=false)
    Z = Defaults.getDefaults("nuclear: charge");    fShell = Shell(0,0);    diff = 0

    # Determine the final shell into which the electron is captured.
    wa = Basics.extractFromConfigurations(Basics.OccupationDifference(), iConf, fConf)
    if length(wa) > 1   error("Incompatible initial and final configurations for a photorecombination cross section.")   end

    for  (k,v) in wa
        diff = diff + v
        if     v == -1     fShell = k    end
    end
    if  diff != -1   error("Incompatible initial and final configurations for a photorecombination cross section.")   end

    # Generate mean-field orbitals in order to extract the transition energies and amplitudes
    grid        = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, rbox = 10.0)
    mfSettings  = AtomicState.MeanFieldSettings(Basics.DFSField(1.0))
    meanField   = Representation("Internal", Nuclear.Model(Z), grid, [fConf], MeanFieldBasis(mfSettings) )
    mfrep       = generate(meanField; output=true)
    fOrbitals   = mfrep["mean-field basis"].orbitals
    fSubsh      = Subshell(fShell.n, -fShell.l -1)
    bEnergy     = - fOrbitals[fSubsh].energy
    omegas  = energies .+ bEnergy

    # Determine the photoionization cross section
    piCss   = Empirical.photoionizationCrossSection(omegas, fConf, iConf, approx, printout=false)

    # Apply the Einstein-Milne relation:  sigma^(PR) = omega^2 / (2 eps c^2) * (g_f/g_i) * sigma^(PI); cf. the same
    # relation in Empirical.photorecombinationCrossSection(..., ScaledHydrogenic). The 1/eps factor is essential:
    # it produces the divergence sigma^(PR) ~ 1/eps as eps -> 0 that characterizes radiative recombination.
    c       = Defaults.getDefaults("speed of light: c")
    gf_gi   = Basics.extractFromConfiguration(Basics.Multiplicity(), fConf) /
              Basics.extractFromConfiguration(Basics.Multiplicity(), iConf)
    prCss   = [gf_gi * omegas[im]^2 / (2 * energies[im] * c^2) * piCss[im] for im in eachindex(omegas)]

    # Report about this estimate
    if  printout
        unCs    = Defaults.getDefaults("unit: cross section");   unEnergy = Defaults.getDefaults("unit: energy")
        energyx = Defaults.convertUnits("energy: from atomic to " * unEnergy, bEnergy)
        omegasx = Float64[];   cssx = Float64[];   
        for omega in omegas  push!(omegasx, Defaults.convertUnits("energy: from atomic to " * unEnergy, omega))   end
        for cs    in prCss   push!(cssx,    Defaults.convertUnits("cross section: from atomic to " * unCs, cs))   end
        sa = "\n* Estimate empirically the photorecombination cross section for a given transition i -> f with the " *
             "following assumptions/simplifications: " *
             "\n    + Use the Einstein-Milne relation to obtain PR cross sections from PI cross sections from the " * 
             "ground level of $fConf " *
             "\n    + PI cross sections are determined in the $(SubString(string(approx), 22)) approximation " * 
             "\n    + iConf = $iConf  -->  fConf = $fConf " * 
             "\n    + Binding energy of $fShell  [$unEnergy]  = $energyx  " *
             "\n    + Omegas [$unEnergy]                 = $omegasx " *
             "\n    + PR Cross sections [$unCs]    = $cssx     " *
             "\n    + Quantity: cross section [$unCs] -- a property of the ion alone, independent of the plasma; fold it with the photon/electron field, or pass iConf/fConf to a *PlasmaRatePerIon / *PlasmaAlpha function, to obtain a rate. " * "\n"
        println(sa)
    end
    
    return( prCss )
end


"""
`Empirical.photorecombinationPlasmaAlpha(eDist::Distribution.AbstractElectronDistribution,
                                         pDist::Distribution.AbstractPhotonDistribution,
                                         iConf::Configuration, fConf::Configuration;
                                         approx::Empirical.AbstractEmpiricalApproximation=ScaledHydrogenic(), printout::Bool=false,
                                         data::PeriodicTable.AbstractEnergyData=PeriodicTable.XrayDataBooklet())`
    ... to estimate empirically the total PR plasma rate coefficient alpha^(PR: total), i.e. the sum of the spontaneous and
        the stimulated contribution, for a transition from iConf -> fConf. The free electrons follow the distribution eDist
        at the electron temperature T_e, while the ambient radiation field is described by pDist at the (independent)
        radiation temperature T. An alpha::Float64 [a.u.] is returned.
        Quantity: a rate coefficient [cm^3/s] -- multiply by the electron number density n_e [1/cm^3] to obtain the rate
            per ion [1/s].

        Note: The stimulated enhancement [1 + nbar(omega)] is evaluated *inside* the integral over the electron energies.
              Photorecombination is a bound-free process, for which the emitted photon energy omega = eps + bEnergy is not
              a single number but varies across the whole electron distribution; there is, therefore, no unique omega at
              which a common factor [1 + nbar] could be taken out of the integral. Such a factorization is exact only for
              bound-bound processes (PD, PX), where omega = E_f - E_i is fixed. It would introduce an error that grows with
              the departure of the radiation field from the electron temperature, i.e. just in the non-LTE regime.
              For pDist = PhotonVacuumField, nbar = 0 and the spontaneous coefficient is recovered.
"""
function photorecombinationPlasmaAlpha(eDist::Distribution.AbstractElectronDistribution,
                                       pDist::Distribution.AbstractPhotonDistribution,
                                       iConf::Configuration, fConf::Configuration;
                                       approx::Empirical.AbstractEmpiricalApproximation=ScaledHydrogenic(), printout::Bool=false,
                                       data::PeriodicTable.AbstractEnergyData=PeriodicTable.XrayDataBooklet())
    # The threshold energy is needed only to place the emitted photon, i.e. only for the stimulated term. A vacuum
    # field carries no photons and, hence, works for every approximation.
    isVacuum = typeof(pDist)  == Distribution.PhotonVacuumField
    isScaled = typeof(approx) == Empirical.ScaledHydrogenic
    if  !isVacuum  &&  !isScaled
        error("A stimulated PR contribution is presently supported only for approx = ScaledHydrogenic(); the threshold " *
              "energy of the capture channel, which fixes the emitted photon energy omega = eps + bEnergy, is not " *
              "available for $approx.")
    end
    Z = Defaults.getDefaults("nuclear: charge");   alpha = 0.;   fShell = Shell(0,0);   diff = 0

    # Determine the shell into which the electron is captured; its binding energy fixes the photon energy
    # omega = eps + bEnergy of each individual recombination channel.
    wa = Basics.extractFromConfigurations(Basics.OccupationDifference(), iConf, fConf)
    if length(wa) > 1   error("Incompatible initial and final configurations for a photorecombination rate coefficient.")   end
    for  (k,v) in wa
        diff = diff + v
        if     v == -1     fShell = k    end
    end
    if  diff != -1   error("Incompatible initial and final configurations for a photorecombination rate coefficient.")   end
    if  isVacuum   bEnergy = 0.
    else           bEnergy = Empirical.scaledBindingEnergy(Z, fShell, fConf, data)
    end
    c       = Defaults.getDefaults("speed of light: c")

    # Perform a Gauss-Legendre integration of
    #     alpha^(PR: total) = int d(eps) f_e(eps; T_e) v(eps) sigma^(PR: spontaneous)(eps) [1 + nbar(eps + bEnergy; T)]
    # with the electron velocity v(eps) = sqrt(2 eps) in atomic units. The integration range is scaled to the electron
    # temperature: a fixed upper limit misses the thermal peak altogether for a cool plasma (at T_e = 0.1 eV, a range
    # of 0 ... 100 a.u. returns a vanishing coefficient) and wastes all mesh points where f_e is negligible.
    # With eenMax = 30 T_e and 96 mesh points, this quadrature is converged to < 0.01% for T_e = 0.1 ... 1000 eV.
    if      hasproperty(eDist, :T)       Te = eDist.T
    elseif  hasproperty(eDist, :Tpar)    Te = max(eDist.Tpar, eDist.Tperp)
    else    error("No electron temperature can be extracted from $eDist.")
    end
    gridGL    = Radial.GridGL(Radial.GridGaussLegendreFinite(), 0., 30*Te, 96; printout=true)
    eEns      = gridGL.t
    if  isScaled   css = Empirical.photorecombinationCrossSection(gridGL.t, iConf, fConf, approx; printout=false, data=data)
    else           css = Empirical.photorecombinationCrossSection(gridGL.t, iConf, fConf, approx; printout=false)
    end
    for  n = 1:gridGL.nt
        ## nbar(omega) follows from the spectral photon number density by dividing out the density of photon modes
        ## omega^2/(pi^2 c^3); this conversion applies to any (isotropic) photon distribution.
        if  isVacuum   occNumber = 0.
        else
            omega     = eEns[n] + bEnergy
            occNumber = Distribution.photonNumberDensity(pDist, omega) * pi^2 * c^3 / omega^2
        end
        alpha     = alpha + css[n] * sqrt(2*eEns[n]) * Distribution.electronEnergyDistribution(eDist, eEns[n]) *
                            (1.0 + occNumber) * gridGL.wt[n]
    end

    # Report about this estimate
    if  printout
        Tex    = Defaults.convertUnits("temperature: from atomic to Kelvin", eDist.T)
        Tpx    = Defaults.convertUnits("temperature: from atomic to Kelvin", pDist.T)
        factor = Defaults.convertUnits("length: from atomic to cm", 1.0)
        factor = factor^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)
        alphax = factor * alpha
        sa = "\n* Estimate empirically the (total) photorecombination plasma rate coefficient alpha for a given transition " *
             "i -> f with the following assumptions/simplifications: " *
             "\n    + Electron field: $eDist,  i.e. T_e [K] = $(Tex). " *
             "\n    + Photon field: $pDist,  i.e. T [K] = $(Tpx). " *
             "\n    + Spontaneous PR cross sections are generated in the $approx approximation. " *
             "\n    + The stimulated enhancement [1 + nbar(eps + bEnergy)] is applied inside the electron-energy integral. " *
             "\n    + iConf = $iConf  -->  fConf = $fConf " *
             "\n    + Plasma rate coefficient alpha^(PR: total) [cm^3/s] = $alphax   " *
             "\n    + Quantity: a rate coefficient [cm^3/s] -- multiply by the electron number density n_e [1/cm^3] for the rate per ion [1/s]. " * "\n"
        println(sa)
    end

    return( alpha )
end


"""
`Empirical.photorecombinationPlasmaAlpha(dist::Distribution.AbstractElectronDistribution,
                                         iConf::Configuration, fConf::Configuration;
                                         approx::Empirical.AbstractEmpiricalApproximation=ScaledHydrogenic(), printout::Bool=false)`
    ... to estimate empirically the spontaneous PR plasma rate coefficient alpha^(PR: spontaneous) for a transition
        from iConf -> fConf by applying some simple approximation as determined by approx. For printout=true, basic information
        are printed about the approximation as well as the results. An alpha::Float64 is returned. This is a short-cut to the
        method above for a vanishing (vacuum) photon field.
        Quantity: a rate coefficient [cm^3/s] -- multiply by the electron number density n_e [1/cm^3] to obtain the rate
            per ion [1/s].
"""
function photorecombinationPlasmaAlpha(dist::Distribution.AbstractElectronDistribution,
                                       iConf::Configuration, fConf::Configuration;
                                       approx::Empirical.AbstractEmpiricalApproximation=ScaledHydrogenic(), printout::Bool=false)
    # The vacuum field has n(omega) = 0 and, hence, no stimulated contribution; this just returns the spontaneous coefficient.
    alpha = Empirical.photorecombinationPlasmaAlpha(dist, Distribution.PhotonVacuumField(0.), iConf, fConf,
                                                    approx=approx, printout=false)

    # Report about this estimate
    if  printout
        Tx     = Defaults.convertUnits("temperature: from atomic to Kelvin", dist.T)
        factor = Defaults.convertUnits("length: from atomic to cm", 1.0)
        factor = factor^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)
        alphax = factor * alpha
        sa = "\n* Estimate empirically the (spontaneous) photorecombination plasma rate coefficient alpha for a given transition " *
             "i -> f with the following assumptions/simplifications: " *
             "\n    + Electron field: $dist,  i.e. T_e [K] = $(Tx). " * 
             "\n    + Spontaneous PR cross sections are generated in the $approx approximation. " * 
             "\n    + iConf = $iConf  -->  fConf = $fConf " * 
             "\n    + Plasma rate coefficient alpha^(PR: spontaneous) [cm^3/s] = $alphax   " *
             "\n    + Quantity: a rate coefficient [cm^3/s] -- multiply by the electron number density n_e [1/cm^3] for the rate per ion [1/s]. " * "\n"
        println(sa)
    end
    
    return( alpha )
end

    
#################################################################################################################################
### Three-Body Recombination (TBR) ##############################################################################################


"""
`Empirical.recombinationConfigurations(conf::Configuration; nLayers::Int64=10)`
    ... to enumerate the configurations that arise from conf by the capture of one additional electron. These are the
        capture channels of (three-body or radiative) recombination: every shell of conf with a vacancy as well as
        all empty (Rydberg) shells up to n_max = n_valence + nLayers, where n_valence is the largest principal
        quantum number occupied in conf. The enumeration truncates there -- the contribution of a capture channel
        falls off with n, but the infinite Rydberg sum would diverge logarithmically, and very high Rydberg states
        are pressure-ionized in any real plasma. Only shells with l <= 3 (s, p, d, f) are included, since Slater's
        screening rules, which fix the binding energies of untabulated shells, are defined up to f electrons.
        A confs::Array{Configuration,1} is returned.
"""
function recombinationConfigurations(conf::Configuration; nLayers::Int64=10)
    confs    = Configuration[]
    nValence = conf.NoElectrons == 0 ? 0 : maximum( sh.n for (sh, occ) in conf.shells if occ > 0 )
    for  n = 1:(nValence + nLayers)
        for  l = 0:min(n-1, 3)
            sh  = Shell(n, l)
            occ = haskey(conf.shells, sh) ? conf.shells[sh] : 0
            if  occ < 2*(2l + 1)
                shells     = deepcopy(conf.shells);    shells[sh] = occ + 1
                push!(confs, Configuration(shells, conf.NoElectrons + 1))
            end
        end
    end

    return( confs )
end


"""
`Empirical.threeBodyRecombinationPlasmaAlpha(eDist::Distribution.ElectronMaxwell,
                                             iConf::Configuration, fConf::Configuration;
                                             approx::Empirical.AbstractEmpiricalApproximation=Lotz1967(), printout::Bool=false,
                                             data::PeriodicTable.AbstractEnergyData=PeriodicTable.XrayDataBooklet())`
    ... to estimate empirically the three-body recombination (TBR) plasma rate coefficient alpha^(TBR) for the capture
        channel iConf + e + e -> fConf + e, by detailed balance with the (Lotz) electron-impact ionization of fConf.
        In LTE, the Saha equation fixes the ratio of the forward and backward rates, so that
            alpha^(TBR) (T_e; i -> f) = alpha^(EII) (T_e; f -> i) * g_f/(2 g_i) * (2 pi/T_e)^(3/2) * exp(P/T_e),
        with P the binding energy of the captured electron and the factor 2 the spin degeneracy of the free electron.
        An alpha::Float64 [a.u.] is returned.
        Quantity: a three-body rate coefficient [cm^6/s] -- multiply by the electron number density *squared*
            n_e^2 [1/cm^6] to obtain the rate per ion [1/s]; TBR therefore dominates over radiative recombination
            at high electron densities.

        Note: The two exponentials exp(+P/T_e) (Saha) and exp(-P/T_e) (threshold suppression of the EII coefficient)
              cancel analytically but *not* numerically: for T_e << P the EII coefficient underflows and the product
              becomes 0 * Inf. The Saha exponential is therefore folded into the integrand, where it combines to the
              well-behaved exp(-(eps - P)/T_e) on the threshold-started mesh. This cancellation requires a Maxwellian
              electron distribution -- detailed balance in this form holds only in (local) thermodynamic equilibrium --
              and the method is restricted accordingly.
"""
function threeBodyRecombinationPlasmaAlpha(eDist::Distribution.ElectronMaxwell,
                                           iConf::Configuration, fConf::Configuration;
                                           approx::Empirical.AbstractEmpiricalApproximation=Lotz1967(), printout::Bool=false,
                                           data::PeriodicTable.AbstractEnergyData=PeriodicTable.XrayDataBooklet())
    if  typeof(approx) != Empirical.Lotz1967
        error("The TBR plasma rate coefficient is presently supported only for approx = Lotz1967().")
    end
    Z = Defaults.getDefaults("nuclear: charge");   alpha = 0.;   cShell = Shell(0,0);   diff = 0

    # Determine the shell into which the electron is captured; its binding energy P is both the threshold of the
    # inverse (EII) process and the energy released in the capture.
    wa = Basics.extractFromConfigurations(Basics.OccupationDifference(), iConf, fConf)
    if length(wa) > 1   error("Incompatible initial and final configurations for a TBR rate coefficient.")   end
    for  (k,v) in wa
        diff = diff + v
        if     v == -1     cShell = k    end
    end
    if  diff != -1   error("Incompatible initial and final configurations for a TBR rate coefficient.")   end
    P      = Empirical.scaledBindingEnergy(Z, cShell, fConf, data)
    xi     = fConf.shells[cShell]
    Te     = eDist.T

    # Lotz's constant A = 4.5e-14 cm^2 eV^2 in atomic units; cf. Empirical.impactIonizationCrossSection.
    Aau    = 4.5e-14 / Defaults.convertUnits("length: from atomic to cm", 1.0)^2 *
             Defaults.convertUnits("energy: from eV to atomic", 1.0)^2

    # Saha (statistical) factor without its exponential; the latter is kept inside the integrand for stability.
    gf_gi  = Basics.extractFromConfiguration(Basics.Multiplicity(), fConf) /
             Basics.extractFromConfiguration(Basics.Multiplicity(), iConf)
    saha   = gf_gi / 2 * (2pi/Te)^1.5

    # Gauss-Legendre integration of the exponentially shifted EII coefficient,
    #     alpha^(TBR) = saha * int d(eps) [f_e(eps) e^(P/T_e)] v(eps) sigma^(EII)(eps)
    # with  f_e(eps) e^(P/T_e) = 2/sqrt(pi) sqrt(eps) T_e^(-3/2) exp(-(eps - P)/T_e),  well-behaved for all T_e.
    gridGL = Radial.GridGL(Radial.GridGaussLegendreFinite(), P, P + 30*Te, 96; printout=false)
    for  n = 1:gridGL.nt
        eps    = gridGL.t[n]
        sigma  = Aau * xi * log(eps/P) / (eps * P)
        fShift = 2/sqrt(pi) * sqrt(eps) / Te^1.5 * exp(-(eps - P)/Te)
        alpha  = alpha + fShift * sqrt(2*eps) * sigma * gridGL.wt[n]
    end
    alpha  = saha * alpha

    # Report about this estimate
    if  printout
        unEnergy = Defaults.getDefaults("unit: energy")
        Tex    = Defaults.convertUnits("temperature: from atomic to Kelvin", Te)
        Px     = Defaults.convertUnits("energy: from atomic to " * unEnergy, P)
        factor = Defaults.convertUnits("length: from atomic to cm", 1.0)
        factor = factor^6 / Defaults.convertUnits("time: from atomic to sec", 1.0)
        alphax = factor * alpha
        sa = "\n* Estimate empirically the three-body recombination plasma rate coefficient alpha for a given capture " *
             "channel i -> f with the following assumptions/simplifications: " *
             "\n    + Electron field: $eDist,  i.e. T_e [K] = $(Tex). " *
             "\n    + Detailed balance (Saha) with the simplified Lotz (1967) EII cross section of the inverse process. " *
             "\n    + iConf = $iConf  -->  fConf = $fConf " *
             "\n    + Binding energy of the captured electron in $cShell [$unEnergy] = $Px " *
             "\n    + Plasma rate coefficient alpha^(TBR) [cm^6/s] = $alphax " *
             "\n    + Quantity: a three-body rate coefficient [cm^6/s] -- multiply by the electron number density squared n_e^2 [1/cm^6] for the rate per ion [1/s]. " * "\n"
        println(sa)
    end

    return( alpha )
end


"""
`Empirical.threeBodyRecombinationPlasmaAlpha(eDist::Distribution.ElectronMaxwell, iConf::Configuration;
                                             nLayers::Int64=10, approx::Empirical.AbstractEmpiricalApproximation=Lotz1967(),
                                             printout::Bool=false,
                                             data::PeriodicTable.AbstractEnergyData=PeriodicTable.XrayDataBooklet())`
    ... to estimate empirically the *total* three-body recombination plasma rate coefficient alpha^(TBR) (T_e; i) of
        the ion iConf, i.e. the sum over all capture channels iConf + e + e -> fConf + e as enumerated by
        Empirical.recombinationConfigurations(iConf; nLayers): every vacancy of iConf and all empty (Rydberg) shells
        with l <= 3 up to n_max = n_valence + nLayers, where the sum is truncated. An alpha::Float64 [a.u.] is returned.
        Quantity: a three-body rate coefficient [cm^6/s] -- multiply by the electron number density *squared*
            n_e^2 [1/cm^6] to obtain the rate per ion [1/s].

        Note: The truncation is not merely numerical but physical: capture into shells with binding energies well
              below T_e is undone by the next collision, so that including ever higher Rydberg shells would inflate
              the coefficient without describing effective recombination; n_max should be chosen in accordance with
              the plasma density (pressure ionization). As a consequence, this *truncated direct-capture* sum is not
              the classical collisional-radiative coefficient of Mansbach & Keck (1969): their scaling
              alpha^(TBR) ~ T_e^(-9/2) rests on the cascade through levels with binding energy ~ T_e, which lie above
              the truncation at low T_e. For H^+ with nLayers = 10, the present sum scales only ~ T_e^(-1.2 ... -1.4)
              between 0.05 and 0.4 eV and lies about an order of magnitude below the classical estimate at 0.1 eV.
"""
function threeBodyRecombinationPlasmaAlpha(eDist::Distribution.ElectronMaxwell, iConf::Configuration;
                                           nLayers::Int64=10, approx::Empirical.AbstractEmpiricalApproximation=Lotz1967(),
                                           printout::Bool=false,
                                           data::PeriodicTable.AbstractEnergyData=PeriodicTable.XrayDataBooklet())
    alpha  = 0.;    confs = Empirical.recombinationConfigurations(iConf, nLayers=nLayers)
    for  fConf in confs
        alpha = alpha + Empirical.threeBodyRecombinationPlasmaAlpha(eDist, iConf, fConf, approx=approx,
                                                                    printout=false, data=data)
    end

    # Report about this estimate
    if  printout
        nValence = iConf.NoElectrons == 0 ? 0 : maximum( sh.n for (sh, occ) in iConf.shells if occ > 0 )
        Tex    = Defaults.convertUnits("temperature: from atomic to Kelvin", eDist.T)
        factor = Defaults.convertUnits("length: from atomic to cm", 1.0)
        factor = factor^6 / Defaults.convertUnits("time: from atomic to sec", 1.0)
        alphax = factor * alpha
        sa = "\n* Estimate empirically the total three-body recombination plasma rate coefficient alpha for a given " *
             "ion i with the following assumptions/simplifications: " *
             "\n    + Electron field: $eDist,  i.e. T_e [K] = $(Tex). " *
             "\n    + Detailed balance (Saha) with the simplified Lotz (1967) EII cross section of each inverse process. " *
             "\n    + iConf = $iConf " *
             "\n    + Sum over $(length(confs)) capture channels: all vacancies of iConf and the empty (Rydberg) shells " *
             "with l <= 3 up to n_max = $(nValence + nLayers); the sum is truncated there. " *
             "\n    + Plasma rate coefficient alpha^(TBR: total) [cm^6/s] = $alphax " *
             "\n    + Quantity: a three-body rate coefficient [cm^6/s] -- multiply by the electron number density squared n_e^2 [1/cm^6] for the rate per ion [1/s]. " * "\n"
        println(sa)
    end

    return( alpha )
end


#################################################################################################################################
### Autoionization (AI) #########################################################################################################


"""
`Empirical.fluorescenceYield(Z::Int64; data::PeriodicTable.AbstractYieldData=PeriodicTable.KrauseAdopted2016())`
    ... to provide the K-shell fluorescence yield omega_K of the neutral element with nuclear charge Z, i.e. the
        probability that a K-shell vacancy is filled by a radiative (K X-ray) rather than a radiationless (Auger)
        transition. A yield::Float64 in [0,1] is returned; cf. Empirical.augerYield for its complement.
"""
function fluorescenceYield(Z::Int64; data::PeriodicTable.AbstractYieldData=PeriodicTable.KrauseAdopted2016())
    if  typeof(data) == PeriodicTable.KrauseAdopted2016
        return( PeriodicTable.fluorescenceYields_KrauseAdopted2016(Z) )
    else
        error("Unsupported fluorescence-yield data set = $data")
    end
end


"""
`Empirical.augerYield(Z::Int64; data::PeriodicTable.AbstractYieldData=PeriodicTable.KrauseAdopted2016())`
    ... to provide the K-shell Auger yield a_K = 1 - omega_K of the neutral element with nuclear charge Z, i.e. the
        probability that a K-shell vacancy decays by a (radiationless) Auger transition. A yield::Float64 in [0,1]
        is returned.
"""
function augerYield(Z::Int64; data::PeriodicTable.AbstractYieldData=PeriodicTable.KrauseAdopted2016())
    return( 1.0 - Empirical.fluorescenceYield(Z, data=data) )
end


"""
`Empirical.augerRate(iConf::Configuration, fConf::Configuration, approx::Empirical.ScaledHydrogenic;
                     printout::Bool=false, data::PeriodicTable.AbstractEnergyData=PeriodicTable.XrayDataBooklet(),
                     yieldData::PeriodicTable.AbstractYieldData=PeriodicTable.KrauseAdopted2016())`
    ... to estimate empirically the (total) K-shell Auger rate A_auger for filling the K-shell vacancy of the transition
        iConf -> fConf, from the fluorescence yield omega_K and the radiative K-shell rate. The two decay channels of a
        K-shell hole share the same total rate Gamma = A_rad + A_auger, and their branching is fixed by the yield,
        A_rad = omega_K Gamma, so that
            A_auger = A_rad (1 - omega_K) / omega_K,
        where A_rad is the (empirical) radiative rate of iConf -> fConf, cf. Empirical.photoemissionEinsteinA. A
        rate::Float64 [a.u.] is returned.
        Quantity: a spontaneous rate [1/s] -- an intrinsic property of the ion, independent of the plasma.

        Note: This estimate exploits that K-shell Auger rates are nearly independent of Z along an isoelectronic
              sequence, whereas radiative rates scale approximately as Z^4; the whole Z dependence is therefore
              carried by omega_K. It presumes that iConf -> fConf refills a *K-shell* (1s) vacancy; omega_K is taken
              for the neutral element and is only a rough guide for an ion.

        Note: The branching ratio A_auger/A_rad = (1 - omega_K)/omega_K is accurate (it is the definition of the
              yield), but the *absolute* Auger rate is only as good as the radiative rate A_rad that feeds it. In the
              ScaledHydrogenic approximation A_rad is itself good to about a factor 1.5, and this uncertainty is
              amplified by the large factor (1 - omega_K)/omega_K for light elements: for neutral Ne the implied
              K-shell width A_rad + A_auger comes out near 2.6 eV, against a measured natural width of ~0.27 eV.
              The result should therefore be read as an order-of-magnitude estimate, most reliable near omega_K ~ 0.5.
"""
function augerRate(iConf::Configuration, fConf::Configuration, approx::Empirical.ScaledHydrogenic;
                   printout::Bool=false, data::PeriodicTable.AbstractEnergyData=PeriodicTable.XrayDataBooklet(),
                   yieldData::PeriodicTable.AbstractYieldData=PeriodicTable.KrauseAdopted2016())
    Z      = round(Int64, Defaults.getDefaults("nuclear: charge"))
    omegaK = Empirical.fluorescenceYield(Z, data=yieldData)
    Arad   = Empirical.photoemissionEinsteinA(iConf, fConf, approx, printout=false, data=data).rate
    Aauger = Arad * (1.0 - omegaK) / omegaK

    # Report about this estimate
    if  printout
        unRate = Defaults.getDefaults("unit: rate")
        Aradx  = Defaults.convertUnits("rate: from atomic to " * unRate, Arad)
        Aaugx  = Defaults.convertUnits("rate: from atomic to " * unRate, Aauger)
        sa = "\n* Estimate empirically the K-shell Auger rate for a given transition i -> f with the " *
             "following assumptions/simplifications: " *
             "\n    + A_auger = A_rad (1 - omega_K)/omega_K, with omega_K the K-shell fluorescence yield ($yieldData). " *
             "\n    + Radiative rate A_rad is estimated in the $approx approximation. " *
             "\n    + iConf = $iConf  -->  fConf = $fConf " *
             "\n    + K-shell fluorescence yield omega_K (Z = $Z) = $omegaK " *
             "\n    + Radiative rate A_rad [$unRate] = $Aradx " *
             "\n    + Auger rate A_auger [$unRate] = $Aaugx " *
             "\n    + Quantity: a spontaneous rate [$unRate] -- an intrinsic property of the ion, independent of the plasma. " * "\n"
        println(sa)
    end

    return( Aauger )
end


#################################################################################################################################
### Electron-impact excitation (EIE) ############################################################################################


"""
`Empirical.effectiveGauntFactor(x2::Float64, isIon::Bool)`
    ... to provide Van Regemorter's (1962) effective Gaunt factor gbar for electron-impact excitation of optically
        allowed transitions; x2 = (eps - deltaE)/deltaE is the *final* electron energy in units of the transition
        energy. For positive ions, the paper prescribes gbar = 0.2 up to x2 = 2 and the hydrogen 1s -> 2p curve
        beyond; its high-energy form [Eq. (10)] gbar = sqrt(3)/(2 pi) ln(x2) reaches 0.2 just at x2 = 2.06, so that
        both branches join seamlessly in gbar = max(0.2, sqrt(3)/(2 pi) ln(x2)). For neutral atoms, the
        near-threshold behavior [Eq. (9)] gbar = 0.074 x (1 + x2) is applied and capped by the ion branch; the
        cross section then vanishes at threshold, as it must for a neutral target. A gbar::Float64 is returned.
"""
function effectiveGauntFactor(x2::Float64, isIon::Bool)
    if  x2 <= 0.    return( isIon ? 0.2 : 0. )   end
    gLog = sqrt(3.0) / (2pi) * log(x2)
    if  isIon    return( max(0.2, gLog) )
    else         x = sqrt(x2);    return( min( 0.074 * x * (1.0 + x2), max(0.2, gLog) ) )
    end
end


"""
`Empirical.impactExcitationCrossSection(energies::Array{Float64,1}, iConf::Configuration, fConf::Configuration,
                                        approx::Empirical.VanRegemorter1962; printout::Bool=false,
                                        aSource::Empirical.AbstractEmpiricalApproximation=Empirical.ScaledHydrogenic(),
                                        data::PeriodicTable.AbstractEnergyData=PeriodicTable.XrayDataBooklet())`
    ... to estimate empirically the electron-impact excitation (EIE) cross section for the optically allowed
        transition iConf -> fConf at the given (free-electron) energies, by using Van Regemorter's (1962) formula
            sigma^(EIE) (eps) = 8 pi / sqrt(3) * pi a_o^2 * (E_Ryd^2 / (eps deltaE)) * f * gbar(x),
        where f is the absorption oscillator strength of the transition, obtained from the (empirical) Einstein-A
        value via f = (g_f/g_i) A / (2 alpha^3 omega^2), and gbar the effective Gaunt factor, cf.
        Empirical.effectiveGauntFactor(). The approx = VanRegemorter1962() argument selects this formula (as opposed
        to, e.g., the forbidden-transition treatment of a different approx type); the separate aSource keyword fixes
        how the Einstein-A value itself is estimated, with both ScaledHydrogenic() and GivenEinsteinA(..) supported.
        A css::Array{Float64,1} [a.u.] is returned.
        Quantity: a cross section [a.u.] -- a property of the ion alone, independent of the plasma; fold it with the
            electron field, or pass the configurations to Empirical.impactExcitationPlasmaAlpha, to obtain a rate.

        Note: Van Regemorter's formula applies to optically allowed (E1) transitions only and is accurate to about
              a factor 2; it fails for forbidden transitions, for which no estimate is returned (an error is raised).
"""
function impactExcitationCrossSection(energies::Array{Float64,1}, iConf::Configuration, fConf::Configuration,
                                      approx::Empirical.VanRegemorter1962; printout::Bool=false,
                                      aSource::Empirical.AbstractEmpiricalApproximation=Empirical.ScaledHydrogenic(),
                                      data::PeriodicTable.AbstractEnergyData=PeriodicTable.XrayDataBooklet())
    Z = Defaults.getDefaults("nuclear: charge")

    # Obtain the (downward) Einstein-A value of the transition fConf -> iConf; iConf is here the *lower* level.
    if      typeof(aSource) == Empirical.GivenEinsteinA
        ea = Empirical.photoemissionEinsteinA(aSource, printout=false)
    elseif  typeof(aSource) == Empirical.ScaledHydrogenic
        ea = Empirical.photoemissionEinsteinA(fConf, iConf, aSource, printout=false, data=data)
    else
        error("Unsupported aSource $aSource for the EIE cross section; use ScaledHydrogenic() or GivenEinsteinA(..).")
    end
    if  ea.multipole != Basics.E1
        error("Van Regemorter's formula applies to optically allowed (E1) transitions only; " *
              "the transition iConf -> fConf was assigned the multipole $(ea.multipole).")
    end
    deltaE = ea.energy
    if  deltaE <= 0.    error("Non-positive transition energy; fConf must lie above iConf.")   end

    # Absorption oscillator strength from the Einstein-A value:  f = (g_f/g_i) A / (2 alpha^3 omega^2)
    alfa   = Defaults.getDefaults("alpha")
    gf_gi  = Basics.extractFromConfiguration(Basics.Multiplicity(), fConf) /
             Basics.extractFromConfiguration(Basics.Multiplicity(), iConf)
    fosc   = gf_gi * ea.rate / (2 * alfa^3 * deltaE^2)

    # Van Regemorter's cross section; in atomic units 8 pi/sqrt(3) E_Ryd^2 = 2 pi/sqrt(3), since E_Ryd = 1/2 Hartree.
    # A positive ion is recognized from the net charge of the target configuration.
    isIon  = Z - iConf.NoElectrons >= 1
    css    = Float64[]
    for  eps in energies
        if  eps > deltaE    x2 = (eps - deltaE) / deltaE
                            push!(css, 2pi / sqrt(3.0) * fosc * Empirical.effectiveGauntFactor(x2, isIon) / (eps * deltaE))
        else                push!(css, 0.)
        end
    end

    # Report about these estimates
    if  printout
        unCs    = Defaults.getDefaults("unit: cross section");   unEnergy = Defaults.getDefaults("unit: energy")
        energyx = Defaults.convertUnits("energy: from atomic to " * unEnergy, deltaE)
        epsx    = Float64[];   cssx = Float64[]
        for eps in energies  push!(epsx, Defaults.convertUnits("energy: from atomic to " * unEnergy, eps))   end
        for cs  in css       push!(cssx, Defaults.convertUnits("cross section: from atomic to " * unCs, cs))   end
        sa = "\n* Estimate empirically the electron-impact excitation cross section for a given transition i -> f with the " *
             "following assumptions/simplifications: " *
             "\n    + Use Van Regemorter's (1962) formula with the effective Gaunt factor for a " *
             (isIon ? "positive ion. " : "neutral atom. ") *
             "\n    + Einstein-A values are determined in the $aSource approximation. " *
             "\n    + iConf = $iConf  -->  fConf = $fConf " *
             "\n    + Transition energy [$unEnergy] = $energyx " *
             "\n    + Oscillator strength f = $fosc " *
             "\n    + Electron energies [$unEnergy] = $epsx " *
             "\n    + EIE cross sections [$unCs] = $cssx " *
             "\n    + Quantity: cross section [$unCs] -- a property of the ion alone, independent of the plasma; fold it with the photon/electron field, or pass iConf/fConf to a *PlasmaRatePerIon / *PlasmaAlpha function, to obtain a rate. " * "\n"
        println(sa)
    end

    return( css )
end


"""
`Empirical.impactExcitationPlasmaAlpha(eDist::Distribution.AbstractElectronDistribution,
                                       iConf::Configuration, fConf::Configuration;
                                       approx::Empirical.VanRegemorter1962=Empirical.VanRegemorter1962(), printout::Bool=false,
                                       aSource::Empirical.AbstractEmpiricalApproximation=Empirical.ScaledHydrogenic(),
                                       data::PeriodicTable.AbstractEnergyData=PeriodicTable.XrayDataBooklet())`
    ... to estimate empirically the electron-impact excitation plasma rate coefficient alpha^(EIE) for the optically
        allowed transition iConf -> fConf and an electron distribution eDist at the electron temperature T_e, by
        folding Van Regemorter's cross section with the electron distribution,
            alpha^(EIE) (T_e; i -> f) = int d(eps) f_e(eps; T_e) v(eps) sigma^(EIE)(eps).
        The approx = VanRegemorter1962() argument selects this formula; the separate aSource keyword fixes how the
        Einstein-A value itself is estimated, with both ScaledHydrogenic() and GivenEinsteinA(..) supported.
        An alpha::Float64 [a.u.] is returned.
        Quantity: a rate coefficient [cm^3/s] -- multiply by the electron number density n_e [1/cm^3] to obtain the rate
            per ion [1/s].
"""
function impactExcitationPlasmaAlpha(eDist::Distribution.AbstractElectronDistribution,
                                     iConf::Configuration, fConf::Configuration;
                                     approx::Empirical.VanRegemorter1962=Empirical.VanRegemorter1962(), printout::Bool=false,
                                     aSource::Empirical.AbstractEmpiricalApproximation=Empirical.ScaledHydrogenic(),
                                     data::PeriodicTable.AbstractEnergyData=PeriodicTable.XrayDataBooklet())
    alpha = 0.
    # The EIE cross section of an ion is finite at threshold; the mesh therefore starts at the threshold, so that
    # the Gauss-Legendre quadrature does not straddle this step, and follows the electron temperature.
    if      typeof(aSource) == Empirical.GivenEinsteinA   deltaE = aSource.energy
    else    ea = Empirical.photoemissionEinsteinA(fConf, iConf, aSource, printout=false, data=data);   deltaE = ea.energy
    end
    if      hasproperty(eDist, :T)       Te = eDist.T
    elseif  hasproperty(eDist, :Tpar)    Te = max(eDist.Tpar, eDist.Tperp)
    else    error("No electron temperature can be extracted from $eDist.")
    end
    gridGL    = Radial.GridGL(Radial.GridGaussLegendreFinite(), deltaE, deltaE + 30*Te, 96; printout=true)
    css       = Empirical.impactExcitationCrossSection(gridGL.t, iConf, fConf, approx; printout=false, aSource=aSource, data=data)
    for  n = 1:gridGL.nt
        alpha = alpha + css[n] * sqrt(2*gridGL.t[n]) * Distribution.electronEnergyDistribution(eDist, gridGL.t[n]) * gridGL.wt[n]
    end

    # Report about this estimate
    if  printout
        Tex    = Defaults.convertUnits("temperature: from atomic to Kelvin", Te)
        factor = Defaults.convertUnits("length: from atomic to cm", 1.0)
        factor = factor^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)
        alphax = factor * alpha
        sa = "\n* Estimate empirically the electron-impact excitation plasma rate coefficient alpha for a given transition " *
             "i -> f with the following assumptions/simplifications: " *
             "\n    + Electron field: $eDist,  i.e. T_e [K] = $(Tex). " *
             "\n    + EIE cross sections follow Van Regemorter (1962); Einstein-A values from the $aSource approximation. " *
             "\n    + iConf = $iConf  -->  fConf = $fConf " *
             "\n    + Plasma rate coefficient alpha^(EIE) [cm^3/s] = $alphax " *
             "\n    + Quantity: a rate coefficient [cm^3/s] -- multiply by the electron number density n_e [1/cm^3] for the rate per ion [1/s]. " * "\n"
        println(sa)
    end

    return( alpha )
end


"""
`Empirical.forbiddenExcitationCrossSection(energies::Array{Float64,1}, iConf::Configuration, fConf::Configuration,
                                           approx::Empirical.ConstantCollisionStrength; printout::Bool=false,
                                           data::PeriodicTable.AbstractEnergyData=PeriodicTable.XrayDataBooklet())`
    ... to estimate empirically the electron-impact excitation (EIE) cross section for the optically forbidden
        (M1 or E2) transition iConf -> fConf at the given (incident free-electron) energies, by assuming a
        constant, energy-independent collision strength Omega = approx.Omega (of order unity by convention when
        no detailed calculation is available), cf. Empirical.ConstantCollisionStrength. The standard
        collision-strength cross section [e.g. D. Osterbrock & G. Ferland, Astrophysics of Gaseous Nebulae and
        Active Galactic Nuclei, 2nd ed., University Science Books (2006)]
            sigma^(EIE) (eps) = pi / (2 eps) * Omega / g_i,
        with eps the incident electron energy and g_i the statistical weight of iConf, cf.
        Basics.extractFromConfiguration(Basics.Multiplicity(), iConf). A css::Array{Float64,1} [a.u.] is returned.
        Quantity: a cross section [a.u.] -- a property of the ion alone, independent of the plasma; fold it with the
            electron field, or pass the configurations to Empirical.forbiddenExcitationPlasmaAlpha, to obtain a rate.

        Note: the transition multipole is determined from Delta-l between the shells involved, exactly as in
              Empirical.photoemissionEinsteinA: Delta-l = 0 -> M1, Delta-l = 2 -> E2. An optically allowed (E1,
              Delta-l = 1) transition is rejected -- use Empirical.impactExcitationCrossSection(energies, iConf,
              fConf, VanRegemorter1962()) instead -- and Delta-l >= 3 (E3 and higher) raises an error, since no
              useful collision-strength estimate is attempted for those higher multipoles.
"""
function forbiddenExcitationCrossSection(energies::Array{Float64,1}, iConf::Configuration, fConf::Configuration,
                                         approx::Empirical.ConstantCollisionStrength; printout::Bool=false,
                                         data::PeriodicTable.AbstractEnergyData=PeriodicTable.XrayDataBooklet())
    Z = Defaults.getDefaults("nuclear: charge");    iShell = fShell = Shell(0,0);    diff = 0

    # Determine the initial and final shells; iConf is here the *lower* (bound) level.
    wa = Basics.extractFromConfigurations(Basics.OccupationDifference(), iConf, fConf)
    if length(wa) > 2   error("Incompatible initial and final configurations for a forbidden EIE cross section.")   end
    for  (k,v) in wa
        diff = diff + v
        if      v == -1    fShell = k
        elseif  v ==  1    iShell = k
        else    error("Incompatible initial and final configurations for a forbidden EIE cross section.")
        end
    end
    if  diff != 0   error("Incompatible initial and final configurations for a forbidden EIE cross section.")   end

    # Determine the multipole from Delta-l, exactly as in Empirical.photoemissionEinsteinA, and restrict to M1/E2.
    dl = abs(iShell.l - fShell.l)
    if      dl == 1    error("The transition iConf -> fConf is optically allowed (E1); use " *
                              "Empirical.impactExcitationCrossSection(energies, iConf, fConf, VanRegemorter1962()) instead.")
    elseif  dl == 0    multipole = M1
    elseif  dl == 2    multipole = E2
    else                error("Forbidden EIE is supported for M1 and E2 transitions only (Delta-l = 0 or 2); " *
                              "the transition iConf -> fConf has Delta-l = $dl.")
    end

    # Unlike Empirical.photoemissionEinsteinA (where iConf is the *upper*, excited level), iConf here is the *lower*
    # (bound) level and fConf the *upper* (excited) one -- the same EIE convention as impactExcitationCrossSection --
    # so iShell (found in iConf) is the deeper, more tightly bound shell and the subtraction order is reversed.
    deltaE = Empirical.scaledBindingEnergy(Z, iShell, iConf, data) - Empirical.scaledBindingEnergy(Z, fShell, fConf, data)
    if  deltaE <= 0.    error("Non-positive transition energy; fConf must lie above iConf.")   end

    gi  = Basics.extractFromConfiguration(Basics.Multiplicity(), iConf)
    css = Float64[]
    for  eps in energies
        if  eps > deltaE    push!(css, pi / (2*eps) * approx.Omega / gi)
        else                push!(css, 0.)
        end
    end

    # Report about these estimates
    if  printout
        unCs    = Defaults.getDefaults("unit: cross section");   unEnergy = Defaults.getDefaults("unit: energy")
        energyx = Defaults.convertUnits("energy: from atomic to " * unEnergy, deltaE)
        epsx    = Float64[];   cssx = Float64[]
        for eps in energies  push!(epsx, Defaults.convertUnits("energy: from atomic to " * unEnergy, eps))   end
        for cs  in css       push!(cssx, Defaults.convertUnits("cross section: from atomic to " * unCs, cs))   end
        sa = "\n* Estimate empirically the electron-impact excitation cross section for a forbidden ($multipole) " *
             "transition i -> f with the following assumptions/simplifications: " *
             "\n    + Use a constant collision strength Omega = $(approx.Omega) (order unity by convention). " *
             "\n    + iConf = $iConf  -->  fConf = $fConf " *
             "\n    + Transition energy [$unEnergy] = $energyx " *
             "\n    + Statistical weight g_i = $gi " *
             "\n    + Electron energies [$unEnergy] = $epsx " *
             "\n    + EIE cross sections [$unCs] = $cssx " *
             "\n    + Quantity: cross section [$unCs] -- a property of the ion alone, independent of the plasma; fold it with the photon/electron field, or pass iConf/fConf to Empirical.forbiddenExcitationPlasmaAlpha, to obtain a rate. " * "\n"
        println(sa)
    end

    return( css )
end


"""
`Empirical.forbiddenExcitationPlasmaAlpha(eDist::Distribution.AbstractElectronDistribution,
                                          iConf::Configuration, fConf::Configuration;
                                          approx::Empirical.ConstantCollisionStrength=Empirical.ConstantCollisionStrength(),
                                          printout::Bool=false,
                                          data::PeriodicTable.AbstractEnergyData=PeriodicTable.XrayDataBooklet())`
    ... to estimate empirically the electron-impact excitation plasma rate coefficient alpha^(EIE) for the
        optically forbidden (M1 or E2) transition iConf -> fConf and an electron distribution eDist at the electron
        temperature T_e, by folding Empirical.forbiddenExcitationCrossSection with the electron distribution,
            alpha^(EIE) (T_e; i -> f) = int d(eps) f_e(eps; T_e) v(eps) sigma^(EIE)(eps).
        For a constant collision strength Omega and a Maxwellian eDist, this reduces analytically to the textbook
        formula alpha^(EIE) (T_e; i -> f) = 8.629e-6 / (g_i sqrt(T_e[K])) * Omega * exp(-deltaE/T_e)  [cm^3/s]
        [e.g. D. Osterbrock & G. Ferland, Astrophysics of Gaseous Nebulae and Active Galactic Nuclei, 2nd ed.,
        University Science Books (2006)], which served as an independent numerical cross check of the
        Gauss-Legendre folding used here. An alpha::Float64 [a.u.] is returned.
        Quantity: a rate coefficient [cm^3/s] -- multiply by the electron number density n_e [1/cm^3] to obtain the
            rate per ion [1/s].
"""
function forbiddenExcitationPlasmaAlpha(eDist::Distribution.AbstractElectronDistribution,
                                        iConf::Configuration, fConf::Configuration;
                                        approx::Empirical.ConstantCollisionStrength=Empirical.ConstantCollisionStrength(),
                                        printout::Bool=false,
                                        data::PeriodicTable.AbstractEnergyData=PeriodicTable.XrayDataBooklet())
    alpha = 0.
    Z = Defaults.getDefaults("nuclear: charge");    iShell = fShell = Shell(0,0);    diff = 0
    wa = Basics.extractFromConfigurations(Basics.OccupationDifference(), iConf, fConf)
    if length(wa) > 2   error("Incompatible initial and final configurations for a forbidden EIE plasma rate coefficient.")   end
    for  (k,v) in wa
        diff = diff + v
        if      v == -1    fShell = k
        elseif  v ==  1    iShell = k
        else    error("Incompatible initial and final configurations for a forbidden EIE plasma rate coefficient.")
        end
    end
    if  diff != 0   error("Incompatible initial and final configurations for a forbidden EIE plasma rate coefficient.")   end
    # iConf is the *lower* (bound) level and fConf the *upper* (excited) one; cf. the note in
    # Empirical.forbiddenExcitationCrossSection.
    deltaE = Empirical.scaledBindingEnergy(Z, iShell, iConf, data) - Empirical.scaledBindingEnergy(Z, fShell, fConf, data)

    if      hasproperty(eDist, :T)       Te = eDist.T
    elseif  hasproperty(eDist, :Tpar)    Te = max(eDist.Tpar, eDist.Tperp)
    else    error("No electron temperature can be extracted from $eDist.")
    end
    gridGL    = Radial.GridGL(Radial.GridGaussLegendreFinite(), deltaE, deltaE + 30*Te, 96; printout=true)
    css       = Empirical.forbiddenExcitationCrossSection(gridGL.t, iConf, fConf, approx; printout=false, data=data)
    for  n = 1:gridGL.nt
        alpha = alpha + css[n] * sqrt(2*gridGL.t[n]) * Distribution.electronEnergyDistribution(eDist, gridGL.t[n]) * gridGL.wt[n]
    end

    # Report about this estimate
    if  printout
        Tex    = Defaults.convertUnits("temperature: from atomic to Kelvin", Te)
        factor = Defaults.convertUnits("length: from atomic to cm", 1.0)
        factor = factor^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)
        alphax = factor * alpha
        sa = "\n* Estimate empirically the electron-impact excitation plasma rate coefficient alpha for a forbidden " *
             "transition i -> f with the following assumptions/simplifications: " *
             "\n    + Electron field: $eDist,  i.e. T_e [K] = $(Tex). " *
             "\n    + EIE cross sections follow a constant collision strength Omega = $(approx.Omega). " *
             "\n    + iConf = $iConf  -->  fConf = $fConf " *
             "\n    + Plasma rate coefficient alpha^(EIE) [cm^3/s] = $alphax " *
             "\n    + Quantity: a rate coefficient [cm^3/s] -- multiply by the electron number density n_e [1/cm^3] for the rate per ion [1/s]. " * "\n"
        println(sa)
    end

    return( alpha )
end


#################################################################################################################################
### Electron-impact ionization (EII) ############################################################################################


"""
`Empirical.impactIonizationCrossSection(energies::Array{Float64,1}, iConf::Configuration, fConf::Configuration,
                                        approx::Empirical.Lotz1967; printout::Bool=false,
                                        data::PeriodicTable.AbstractEnergyData=PeriodicTable.XrayDataBooklet())`
    ... to estimate empirically the electron-impact ionization (EII) cross section for a transition iConf -> fConf
        at the given (free-electron) energies, by using the simplified formula of Lotz (1967, 1968)
            sigma^(EII) (eps) = A xi ln(eps/P) / (eps P),        A = 4.5e-14 cm^2 eV^2,
        where P is the binding energy of the ionized shell, cf. Empirical.scaledBindingEnergy(), and xi the number
        of its equivalent electrons. Lotz's full formula carries subshell constants a_i = 2.9 ... 4.5e-14 cm^2 eV^2
        and two further parameters (b_i, c_i) that lower the near-threshold cross section of neutral atoms; the
        simplified form with A = 4.5e-14 applies best to ions and overestimates light neutral atoms by ~40%
        (accuracy +40/-30%). A css::Array{Float64,1} [a.u.] is returned.
        Quantity: a cross section [a.u.] -- a property of the ion alone, independent of the plasma; fold it with the
            electron field, or pass the configurations to Empirical.impactIonizationPlasmaAlpha, to obtain a rate.
"""
function impactIonizationCrossSection(energies::Array{Float64,1}, iConf::Configuration, fConf::Configuration,
                                      approx::Empirical.Lotz1967; printout::Bool=false,
                                      data::PeriodicTable.AbstractEnergyData=PeriodicTable.XrayDataBooklet())
    Z = Defaults.getDefaults("nuclear: charge");   iShell = Shell(0,0);   diff = 0

    # Determine the ionized shell, its binding energy and its number of equivalent electrons.
    wa = Basics.extractFromConfigurations(Basics.OccupationDifference(), iConf, fConf)
    if length(wa) > 1   error("Incompatible initial and final configurations for an EII cross section.")   end
    for  (k,v) in wa
        diff = diff + v
        if      v == 1     iShell = k    end
    end
    if  diff != 1   error("Incompatible initial and final configurations for an EII cross section.")   end
    P      = Empirical.scaledBindingEnergy(Z, iShell, iConf, data)
    xi     = iConf.shells[iShell]

    # Lotz's constant A = 4.5e-14 cm^2 eV^2, expressed in atomic units [a_o^2 Hartree^2].
    Aau    = 4.5e-14 / Defaults.convertUnits("length: from atomic to cm", 1.0)^2 *
             Defaults.convertUnits("energy: from eV to atomic", 1.0)^2

    css    = Float64[]
    for  eps in energies
        if  eps > P    push!(css, Aau * xi * log(eps/P) / (eps * P))   else   push!(css, 0.)   end
    end

    # Report about these estimates
    if  printout
        unCs    = Defaults.getDefaults("unit: cross section");   unEnergy = Defaults.getDefaults("unit: energy")
        energyx = Defaults.convertUnits("energy: from atomic to " * unEnergy, P)
        epsx    = Float64[];   cssx = Float64[]
        for eps in energies  push!(epsx, Defaults.convertUnits("energy: from atomic to " * unEnergy, eps))   end
        for cs  in css       push!(cssx, Defaults.convertUnits("cross section: from atomic to " * unCs, cs))   end
        sa = "\n* Estimate empirically the electron-impact ionization cross section for a given transition i -> f with the " *
             "following assumptions/simplifications: " *
             "\n    + Use the simplified Lotz (1967) formula with A = 4.5e-14 cm^2 eV^2 (best for ions; ~40% high " *
             "for light neutral atoms). " *
             "\n    + Binding energy (ionization threshold) taken from $data " *
             "\n    + iConf = $iConf  -->  fConf = $fConf " *
             "\n    + Binding energy of $iShell [$unEnergy] = $energyx;  equivalent electrons xi = $xi " *
             "\n    + Electron energies [$unEnergy] = $epsx " *
             "\n    + EII cross sections [$unCs] = $cssx " *
             "\n    + Quantity: cross section [$unCs] -- a property of the ion alone, independent of the plasma; fold it with the photon/electron field, or pass iConf/fConf to a *PlasmaRatePerIon / *PlasmaAlpha function, to obtain a rate. " * "\n"
        println(sa)
    end

    return( css )
end


"""
`Empirical.impactIonizationPlasmaAlpha(eDist::Distribution.AbstractElectronDistribution,
                                       iConf::Configuration, fConf::Configuration;
                                       approx::Empirical.AbstractEmpiricalApproximation=Lotz1967(), printout::Bool=false,
                                       data::PeriodicTable.AbstractEnergyData=PeriodicTable.XrayDataBooklet())`
    ... to estimate empirically the electron-impact ionization plasma rate coefficient alpha^(EII) for a transition
        iConf -> fConf and an electron distribution eDist at the electron temperature T_e, by folding the Lotz cross
        section with the electron distribution,
            alpha^(EII) (T_e; i -> f) = int d(eps) f_e(eps; T_e) v(eps) sigma^(EII)(eps).
        An alpha::Float64 [a.u.] is returned.
        Quantity: a rate coefficient [cm^3/s] -- multiply by the electron number density n_e [1/cm^3] to obtain the rate
            per ion [1/s].
"""
function impactIonizationPlasmaAlpha(eDist::Distribution.AbstractElectronDistribution,
                                     iConf::Configuration, fConf::Configuration;
                                     approx::Empirical.AbstractEmpiricalApproximation=Lotz1967(), printout::Bool=false,
                                     data::PeriodicTable.AbstractEnergyData=PeriodicTable.XrayDataBooklet())
    if  typeof(approx) != Empirical.Lotz1967
        error("The EII plasma rate coefficient is presently supported only for approx = Lotz1967().")
    end
    Z = Defaults.getDefaults("nuclear: charge");   alpha = 0.;   iShell = Shell(0,0);   diff = 0

    # Determine the threshold in order to start the mesh there; cf. photoionizationPlasmaRatePerIon.
    wa = Basics.extractFromConfigurations(Basics.OccupationDifference(), iConf, fConf)
    if length(wa) > 1   error("Incompatible initial and final configurations for an EII rate coefficient.")   end
    for  (k,v) in wa
        diff = diff + v
        if      v == 1     iShell = k    end
    end
    if  diff != 1   error("Incompatible initial and final configurations for an EII rate coefficient.")   end
    P = Empirical.scaledBindingEnergy(Z, iShell, iConf, data)

    if      hasproperty(eDist, :T)       Te = eDist.T
    elseif  hasproperty(eDist, :Tpar)    Te = max(eDist.Tpar, eDist.Tperp)
    else    error("No electron temperature can be extracted from $eDist.")
    end
    gridGL    = Radial.GridGL(Radial.GridGaussLegendreFinite(), P, P + 30*Te, 96; printout=true)
    css       = Empirical.impactIonizationCrossSection(gridGL.t, iConf, fConf, approx; printout=false, data=data)
    for  n = 1:gridGL.nt
        alpha = alpha + css[n] * sqrt(2*gridGL.t[n]) * Distribution.electronEnergyDistribution(eDist, gridGL.t[n]) * gridGL.wt[n]
    end

    # Report about this estimate
    if  printout
        Tex    = Defaults.convertUnits("temperature: from atomic to Kelvin", Te)
        factor = Defaults.convertUnits("length: from atomic to cm", 1.0)
        factor = factor^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)
        alphax = factor * alpha
        sa = "\n* Estimate empirically the electron-impact ionization plasma rate coefficient alpha for a given transition " *
             "i -> f with the following assumptions/simplifications: " *
             "\n    + Electron field: $eDist,  i.e. T_e [K] = $(Tex). " *
             "\n    + EII cross sections follow the simplified Lotz (1967) formula. " *
             "\n    + iConf = $iConf  -->  fConf = $fConf " *
             "\n    + Plasma rate coefficient alpha^(EII) [cm^3/s] = $alphax " *
             "\n    + Quantity: a rate coefficient [cm^3/s] -- multiply by the electron number density n_e [1/cm^3] for the rate per ion [1/s]. " * "\n"
        println(sa)
    end

    return( alpha )
end


