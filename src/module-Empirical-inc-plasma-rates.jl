
#################################################################################################################################
### Photoemission (PE) ##########################################################################################################


"""
`Empirical.photoemissionEinsteinA(iConf::Configuration, fConf::Configuration, approx::Empirical.ScaledHydrogenic;
                                  printout::Bool=false)`  
    ... to estimate empirically the Einstein-A value for a transition from iConf -> fConf by using a simple
        hydrogenic scaling with an effective nuclear charge. A named triple
        (multipole::EmMultipole=, energy::Float64=, rate::Float64=) is returned. 
"""
function photoemissionEinsteinA(iConf::Configuration, fConf::Configuration, approx::Empirical.ScaledHydrogenic;
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

    # Determine the effective charges Zi and Zf felt by an electron in the initial and final shells;
    # Just substract all inner-shell electrons and use a pure hydrogenic scaling of the binding energy 
    # (assume a charge 0.5 due to the other electrons in iShell or fShell)
    ce      = 0;   
    for  (shell, v) in iConf.shells   
        if  shell.n < iShell.n  ||  (shell.n == iShell.n  &&  shell.l < iShell.l)   ce = ce + v   end
    end
    Zi      = Z - ce - 0.5;    rxi = (3*iShell.n^2 - iShell.l + (iShell.l+1)) / (2*Zi)
    ce      = 0;   
    for  (shell, v) in fConf.shells   
        if  shell.n < fShell.n  ||  (shell.n == fShell.n  &&  shell.l < fShell.l)   ce = ce + v   end
    end
    Zf      = Z - ce - 0.5;    rxf = (3*fShell.n^2 - fShell.l + (fShell.l+1)) / (2*Zf)
    tEnergy = Zi^2 / (iShell.n - 0.07)^2 - Zf^2 / fShell.n^2
    
    # Determine the multipolarity of the transition
    if  multipole === missing    ||   diff != 0  
    else
        if       abs(iShell.l - fShell.l) == 0        multipole = M1
        elseif   abs(iShell.l - fShell.l) == 1        multipole = E1
        elseif   abs(iShell.l - fShell.l) == 2        multipole = E2
        elseif   abs(iShell.l - fShell.l) == 1        multipole = E3
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
    else                                triple = (multipole, 0., 0.) 
    end
    
    # Report about this estimate
    if  printout
        unRate = Defaults.getDefaults("unit: rate");   unEnergy = Defaults.getDefaults("unit: energy")
        energyx = Defaults.convertUnits("energy: from atomic to " * unEnergy, tEnergy)
        ratex   = Defaults.convertUnits("rate: from atomic to "   * unRate,   rate)
        sa = "\n* Estimate empirically the Einstein-A value for a given transition i -> f with the " *
             "following assumptions/simplifications: " *
             "\n    + Use a simple hydrogenic scaling for the energy and A-value. " * 
             "\n    + iConf = $iConf  -->  fConf = $fConf " * 
             "\n    + Extract transition energy from binding energies of shells $iShell -> $fShell " * 
             "\n    + $multipole transition with energy [$unEnergy] = $energyx  " * 
             "\n    + Rate [$unRate]                     = $ratex " * "\n"
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
        elseif   abs(iShell.l - fShell.l) == 1        multipole = E3
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
             "\n    + Rate [$unRate]                     = $ratex " * "\n"
        println(sa)
    end
    
    return( (multipole = multipole, energy = tEnergy, rate = rate) )
end


"""
`Empirical.photoemissionEinsteinA(approx::Empirical.GivenEinsteinA; printout::Bool=false)`  
    ... to simply return the Einstein-A value for the given transition. A named
        (multipole::EmMultipole=, energy::Float64=, rate::Float64=) is returned. 
"""
function photoemissionEinsteinA(approx::Empirical.GivenEinsteinA; printout::Bool=false) 
    
    # Report about this estimate
    if  printout
        unRate = Defaults.getDefaults("unit: rate");   unEnergy = Defaults.getDefaults("unit: energy")
        energyx = Defaults.convertUnits("energy: from atomic to " * unEnergy, approx.energy)
        ratex   = Defaults.convertUnits("rate: from atomic to "   * unRate,   approx.rate)
        sa = "\n* User-given Einstein-A value for a given transition i -> f: "   *
             "\n    + $(approx.multipole) transition with energy = $energyx  " * unEnergy *
             "\n    + rate = $ratex " * unRate * "\n"
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
"""
function photoexcitationPlasmaRatePerIon(dist::Distribution.AbstractPhotonDistribution,
                                         iConf::Configuration, fConf::Configuration; 
                                         approx::Empirical.AbstractEmpiricalApproximation=UsingJAC(), printout::Bool=false)
    EinsteinA = Empirical.photoemissionEinsteinA(iConf, fConf, approx, printout=false)
    pnDensity = Distribution.photonNumberDensity(dist, EinsteinA.energy)
    rate      = Basics.extractFromConfiguration(Basics.Multiplicity(), fConf) / 
                Basics.extractFromConfiguration(Basics.Multiplicity(), iConf) * EinsteinA.rate * pnDensity

    # Report about this estimate
    if  printout
        Tx      = Defaults.convertUnits("temperature: from atomic to Kelvin", dist.T)
        unRate  = Defaults.getDefaults("unit: rate");   unEnergy = Defaults.getDefaults("unit: energy")
        energyx = Defaults.convertUnits("energy: from atomic to " * unEnergy, EinsteinA.energy)
        ratex   = Defaults.convertUnits("rate: from atomic to "   * unRate,   rate)
        sa = "\n* Estimate empirically the photoexcitation plasma rate R^(PX: per ion) (T; i -> f) for a given transition " *
             "i -> f with the following assumptions/simplifications: " *
             "\n    + Photon field follows a $(SubString(string(dist), 22)) at temperature T [K] = $(Tx). " * 
             "\n    + Einstein-A values are determined in the $(SubString(string(approx), 22)) approximation " * 
             "\n    + iConf = $iConf  -->  fConf = $fConf " * 
             "\n    + Transition energy [$unEnergy]                                = $energyx " *
             "\n    + Plasma rate per ion R^(PX: per ion) (T; i -> f) [$unRate] = $ratex \n"
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
"""
function photodeexcitationPlasmaRatePerIon(dist::Distribution.AbstractPhotonDistribution,
                                           iConf::Configuration, fConf::Configuration; 
                                           approx::Empirical.AbstractEmpiricalApproximation=UsingJAC(), printout::Bool=false)
    EinsteinA = Empirical.photoemissionEinsteinA(iConf, fConf, approx, printout=false)
    pnDensity = Distribution.photonNumberDensity(dist, EinsteinA.energy)
    rate      = EinsteinA.rate * (1.0 + pnDensity)

    # Report about this estimate
    if  printout
        Tx      = Defaults.convertUnits("temperature: from atomic to Kelvin", dist.T)
        unRate  = Defaults.getDefaults("unit: rate");   unEnergy = Defaults.getDefaults("unit: energy")
        energyx = Defaults.convertUnits("energy: from atomic to " * unEnergy, EinsteinA.energy)
        ratex   = Defaults.convertUnits("rate: from atomic to "   * unRate,   rate)
        sa = "\n* Estimate empirically the (total) photodeexcitation plasma rate R^(PD: total, per ion) (T; i -> f) " *
             "for a given transition i -> f" *
             "\n  with the following assumptions/simplifications: " *
             "\n    + Photon field follows a $(SubString(string(dist), 22)) at temperature T [K] = $(Tx). " * 
             "\n    + Einstein-A values are determined in the $(SubString(string(approx), 22)) approximation " * 
             "\n    + iConf = $iConf  -->  fConf = $fConf " * 
             "\n    + Transition energy [$unEnergy]                               = $energyx " *
             "\n    + Plasma rate R^(PD: total, per ion) (T; i -> f) [$unRate] = $ratex \n"
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
"""
function photoionizationCrossSection(omega::Float64, iConf::Configuration, fConf::Configuration; 
                                     approx::Empirical.AbstractEmpiricalApproximation=ScaledHydrogenic(), printout::Bool=false)
    
    cs = Empirical.photoionizationCrossSection([omega], iConf, fConf, approx, printout=printout) 
    
    return( cs[1] )
end


"""
`Empirical.photoionizationCrossSection(omegas::Array{Float64,1}, iConf::Configuration, fConf::Configuration, 
                                       approx::Empirical.ScaledHydrogenic; printout::Bool=false)`  
    ... to estimate empirically the photoionization cross sections for a transition from iConf -> fConf by using the binding
        energy (ionization potential) of the ionized shell from JAC and Kramer's (1923) empirical formula.
        A css::Array{Float64,1} [a.u.] is returned. 
"""
function photoionizationCrossSection(omegas::Array{Float64,1}, iConf::Configuration, fConf::Configuration, 
                                     approx::Empirical.ScaledHydrogenic; printout::Bool=false) 
    Z = Defaults.getDefaults("nuclear: charge");    iShell = Shell(0,0);    diff = 0;   zeroCss = false
    
    # Determine the initial shell and its binding (threshold) energy; set all css = 0, if the occupation of 
    # configurations differ by more than 1.
    wa = Basics.extractFromConfigurations(Basics.OccupationDifference(), iConf, fConf)
    if length(wa) > 1   zeroCss = true   end 
 
    for  (k,v) in wa
        diff = diff + v
        if      v == 1     iShell = k    end 
    end
    if  diff != 1   error("Incompatible initial and final configurations for a photoionization cross section.")   end
    
    # Just substract all inner-shell electrons and use a pure hydrogenic scaling of the binding energy 
    # (assume a charge 0.5 due to the other electrons in iShell)
    cShells = Basics.extractFromConfiguration(Basics.ClosedShells(), fConf)
    ce      = 0;   for  cShell in cShells   ce = ce + 2 * (2*cShell.l + 1)   end
    Zi      = Z - ce - 1.0
    bEnergy = Zi^2 / iShell.n^2
    factor  = 64 * pi^2 * Defaults.getDefaults("alpha") / (3 * sqrt(3) * iShell.n)
    
    # Compute all cross sections in the ScaledHydrogenic approximation
    css = Float64[]
    for  omega in omegas   
        if  omega > bEnergy    push!(css, factor * (bEnergy / omega)^3 )   else    push!(css, 0.)  end 
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
             "\n    + Use a simple hydrogenic scaling of the PI cross section following Bethe-Salpeter (1957). " * 
             "\n    + Hydrogenic estimate of the binding energy (photoionization threshold) with Z^(eff) = $Zi " * 
             "\n    + iConf = $iConf  -->  fConf = $fConf " * 
             "\n    + Binding energy of $iShell   = $energyx  " * unEnergy *
             "\n    + Omegas [$unEnergy]            = $omegasx " *
             "\n    + Cross sections [$unCs]  = $cssx     \n"
        println(sa)
    end
    
    return( css )
end


"""
`Empirical.photoionizationCrossSection(omegas::Array{Float64,1}, iConf::Configuration, fConf::Configuration, 
                                       approx::Empirical.UsingJAC; printout::Bool=false)`  
    ... to estimate empirically the photoionization cross section for a transition from iConf -> fConf by using the binding
        energy (ionization potential) of the ionized shell from JAC and Kramer's (1923) empirical formula.
        A cs::Float64 is returned. 
"""
function photoionizationCrossSection(omegas::Array{Float64,1}, iConf::Configuration, fConf::Configuration, 
                                     approx::Empirical.UsingJAC; printout::Bool=false) 
    Z = Defaults.getDefaults("nuclear: charge");    iShell = Shell(0,0);    diff = 0;   zeroCss = false
    
    # Determine the initial shell and its binding (threshold) energy; set all css = 0, if the occupation of 
    # configurations differ by more than 1.
    wa = Basics.extractFromConfigurations(Basics.OccupationDifference(), iConf, fConf)
    if length(wa) > 1   zeroCss = true   end 
 
    for  (k,v) in wa
        diff = diff + v
        if      v == 1     iShell = k    end 
    end
    if  diff != 1   error("Incompatible initial and final configurations for a photoionization cross section.")   end

    # Generate mean-field orbitals in order to extract the transition energies and amplitudes
    grid        = Radial.Grid(Radial.Grid(true), rnt = 4.0e-6, h = 5.0e-2, rbox = 10.0) 
    mfSettings  = AtomicState.MeanFieldSettings(Basics.DFSField(1.0))
    meanField   = Representation("Internal", Nuclear.Model(Z), grid, [iConf], MeanFieldBasis(mfSettings) )
    mfrep       = generate(meanField; output=true)
    iBasis      = mfrep["mean-field basis"]
    iSubsh      = Subshell(iShell.n, -iShell.l -1)
    cSubsh      = Subshell(101, -iShell.l -1 -1)
    bEnergy     = -iBasis.orbitals[iSubsh].energy
    
    # Compute single-electron photoionization cross sections
    css      = Float64[];    contSettings = Continuum.Settings(false, 290);    nm = Nuclear.Model(Z)  
    nucPot   = Nuclear.nuclearPotential(nm, grid)
    dfsPot   = Basics.computePotential(Basics.DFSField(1.0), grid, iBasis)
    localPot = Basics.add(nucPot, dfsPot)
    
    for  omega in omegas
        cOrbital, phase, norm = Continuum.generateOrbitalLocalPotential(omega-bEnergy, cSubsh, localPot, contSettings)
        amp      = InteractionStrength.MabEmissionJohnsony(Basics.E1, Basics.Babushkin, omega, cOrbital,  
                                                           iBasis.orbitals[iSubsh], grid)
        wa       = 8 * pi^3 / Defaults.getDefaults("alpha") / omega / 
                   Basics.extractFromConfiguration(Basics.Multiplicity(), iConf) *
                   Basics.extractFromConfiguration(Basics.Multiplicity(), fConf) / (2*(2*iShell.l+1))
        push!(css, wa * 1.0 * abs(amp)^2)  # A contribution 0.0 is considered due to the l-1 partial waves
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
             "\n    + Using a JAC mean-field approach for binding energy (threshold energy) and the PI cross sections. " * 
             "\n    + Mean-field computations for inital configuration & one-electron PI matrix elements." * 
             "\n    + iConf = $iConf  -->  fConf = $fConf " * 
             "\n    + Binding energy of $iShell   = $energyx  " * unEnergy *
             "\n    + Omegas [$unEnergy]            = $omegasx " *
             "\n    + Cross sections [$unCs]  = $cssx     \n"
        println(sa)
    end
    
    return( css )
end


"""
`Empirical.photoionizationPlasmaAlpha(dist::Distribution.AbstractPhotonDistribution, 
                                      iConf::Configuration, fConf::Configuration; 
                                      approx::Empirical.AbstractEmpiricalApproximation=ScaledHydrogenic(), printout::Bool=false)` 
                                     
    ... to estimate empirically the PI plasma rate coefficient alpha^(PI) for a transition from iConf -> fConf by applying some simple
        approximation as determined by approx. For printout=true, basic information are printed about the approximation as well
        as the results. An alpha::Float64 is returned.
"""
function photoionizationPlasmaAlpha(dist::Distribution.AbstractPhotonDistribution,
                                    iConf::Configuration, fConf::Configuration; 
                                    approx::Empirical.AbstractEmpiricalApproximation=ScaledHydrogenic(), printout::Bool=false)
    alpha = 0.;    omegaMax = 100.
    # Perform a Gauss-Legendre integration by selecting the photon number distribution and cross sectins at the required energies
    gridGL    = Radial.GridGL(Radial.GridGaussLegendreFinite(), 0., omegaMax, 24; printout=true)
    css       = Empirical.photoionizationCrossSection(gridGL.t, iConf, fConf, approx; printout=false)
    for  n = 1:gridGL.nt
        alpha = alpha + css[n] * Distribution.photonNumberDensity(dist, gridGL.t[n]) * gridGL.wt[n]    
    end
    
    # Report about this estimate
    if  printout
        Tx     = Defaults.convertUnits("temperature: from atomic to Kelvin", dist.T)
        factor = Defaults.convertUnits("length: from atomic to cm", 1.0)
        factor = factor^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)
        alphax = factor * alpha
        sa = "\n* Estimate empirically the photoionization plasma rate coefficient alpha for a given transition i -> f with the " *
             "following assumptions/simplifications: " *
             "\n    + Photon field follows a $(SubString(string(dist), 22)) at temperature T [K] = $(Tx). " * 
             "\n    + PI cross sections are generated in the $(SubString(string(approx), 22)) approximation. " * 
             "\n    + iConf = $iConf  -->  fConf = $fConf " * 
             "\n    + Plasma rate coefficient alpha^(PI) [cm^3/s] = $alphax \n"
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
"""
function photorecombinationCrossSection(energy::Float64, iConf::Configuration, fConf::Configuration; 
                                        approx::Empirical.AbstractEmpiricalApproximation=ScaledHydrogenic(), printout::Bool=false)
    
    cs = Empirical.photorecombinationCrossSection([energy], iConf, fConf, approx, printout=printout) 
    
    return( cs )
end


"""
`Empirical.photorecombinationCrossSection(energies::Array{Float64,1}, iConf::Configuration, fConf::Configuration, 
                                          approx::Empirical.ScaledHydrogenic; printout::Bool=false)`  
    ... to estimate empirically the (spontaneous) photorecombination cross section for a transition from iConf -> fConf 
        by using the Einstein-Milne relation and the binding energy (ionization potential) of the ionized shell 
        from JAC and Kramer's (1923) empirical formula.A css::Array{Float64,1} [a.u.] is returned. 
"""
function photorecombinationCrossSection(energies::Array{Float64,1}, iConf::Configuration, fConf::Configuration, 
                                        approx::Empirical.ScaledHydrogenic; printout::Bool=false) 
    Z = Defaults.getDefaults("nuclear: charge");    fShell = Shell(0,0);    diff = 0;   zeroCss = false
    
    # Determine the initial shell and its binding (threshold) energy; set all css = 0, if the occupation of 
    # configurations differ by more than 1.
    wa = Basics.extractFromConfigurations(Basics.OccupationDifference(), iConf, fConf)
    if length(wa) > 1   zeroCss = true   end 
 
    for  (k,v) in wa
        diff = diff + v
        if     v == -1     fShell = k    end 
    end
    if  diff != -1   error("Incompatible initial and final configurations for a photorecombination cross section.")   end
    
    # Just substract all inner-shell electrons and use a pure hydrogenic scaling of the binding energy 
    # (assume a charge 0.5 due to the other electrons in iShell)
    ce      = iConf.NoElectrons   
    Zf      = Z - iConf.NoElectrons + 0.1   # effective Z felt by the captured electron in fShell
    bEnergy = Zf^2 / fShell.n^2
    omegas  = energies .+ bEnergy
    
    # Determine the photoionization cross section
    piCss   = Empirical.photoionizationCrossSection(omegas, fConf, iConf, approx, printout=false)
    
    # Apply the Einstein-Milne relation
    prCss   = Float64[]
    factor  = pi^2 * Defaults.getDefaults("speed of light: c")^2 *
              Basics.extractFromConfiguration(Basics.Multiplicity(), fConf) / 
              Basics.extractFromConfiguration(Basics.Multiplicity(), iConf)
    for (im, omega) in enumerate(omegas)    push!(prCss, factor * piCss[im] / omega^2)   end 
    
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
             "\n    + PR Cross sections [$unCs]    = $cssx     \n"
        println(sa)
    end
    
    return( prCss )
end


"""
`Empirical.photorecombinationCrossSection(energies::Array{Float64,1}, iConf::Configuration, fConf::Configuration, 
                                          approx::Empirical.UsingJAC; printout::Bool=false)`  
    ... to estimate empirically the (spontaneous) photorecombination cross section for a transition from iConf -> fConf 
        by using the Einstein-Milne relation and the binding energy (ionization potential) of the ionized shell 
        from JAC and Kramer's (1923) empirical formula.A css::Array{Float64,1} [a.u.] is returned. 
"""
function photorecombinationCrossSection(energies::Array{Float64,1}, iConf::Configuration, fConf::Configuration, 
                                        approx::Empirical.UsingJAC; printout::Bool=false) 
    Z = Defaults.getDefaults("nuclear: charge");    fShell = Shell(0,0);    diff = 0
    
    # Determine the initial shell and its binding (threshold) energy; set all css = 0, if the occupation of 
    # configurations differ by more than 1.
    wa = Basics.extractFromConfigurations(Basics.OccupationDifference(), iConf, fConf)
    if length(wa) > 1   zeroCss = true   end 
 
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
    
    # Apply the Einstein-Milne relation
    prCss   = Float64[]
    factor  = pi^2 * Defaults.getDefaults("speed of light: c")^2 *
              Basics.extractFromConfiguration(Basics.Multiplicity(), fConf) / 
              Basics.extractFromConfiguration(Basics.Multiplicity(), iConf)
    for (im, omega) in enumerate(omegas)    push!(prCss, factor * piCss[im] / omega^2)   end 
    
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
             "\n    + PR Cross sections [$unCs]    = $cssx     \n"
        println(sa)
    end
    
    return( prCss )
end


"""
`Empirical.photorecombinationPlasmaAlpha(dist::Distribution.AbstractElectronDistribution, 
                                         iConf::Configuration, fConf::Configuration; 
                                         approx::Empirical.AbstractEmpiricalApproximation=ScaledHydrogenic(), printout::Bool=false)` 
    ... to estimate empirically the spontaneous PR plasma rate coefficient alpha^(PR: spontaneous) for a transition 
        from iConf -> fConf by applying some simple approximation as determined by approx. For printout=true, basic information 
        are printed about the approximation as well as the results. An alpha::Float64 is returned.
"""
function photorecombinationPlasmaAlpha(dist::Distribution.Distribution.AbstractElectronDistribution, 
                                       iConf::Configuration, fConf::Configuration; 
                                       approx::Empirical.AbstractEmpiricalApproximation=ScaledHydrogenic(), printout::Bool=false)
    alpha = 0.;    eenMax = 100.
    # Perform a Gauss-Legendre integration by selecting the electron energy distribution and cross sections at the required energies
    gridGL    = Radial.GridGL(Radial.GridGaussLegendreFinite(), 0., eenMax, 24; printout=true)
    eEns      = gridGL.t
    css       = Empirical.photorecombinationCrossSection(gridGL.t, iConf, fConf, approx; printout=false)
    for  n = 1:gridGL.nt
        alpha = alpha + css[n] * eEns[n] * Distribution.electronEnergyDistribution(dist, eEns[n]) * gridGL.wt[n]    
    end
    
    # Report about this estimate
    if  printout
        Tx     = Defaults.convertUnits("temperature: from atomic to Kelvin", dist.T)
        factor = Defaults.convertUnits("length: from atomic to cm", 1.0)
        factor = factor^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)
        alphax = factor * alpha
        sa = "\n* Estimate empirically the (spontaneous) photorecombination plasma rate coefficient alpha for a given transition " *
             "i -> f with the following assumptions/simplifications: " *
             "\n    + Electron field follows a $(SubString(string(dist), 22)) at temperature T [K] = $(Tx). " * 
             "\n    + Spontaneous PR cross sections are generated in the $approx approximation. " * 
             "\n    + iConf = $iConf  -->  fConf = $fConf " * 
             "\n    + Plasma rate coefficient alpha^(PR: spontaneous) [cm^3/s] = $alphax   \n"
        println(sa)
    end
    
    return( alpha )
end

    
#################################################################################################################################
### Three-Body Recombination (TBR) ##############################################################################################


#################################################################################################################################
### Autoionization (AI) #########################################################################################################


#################################################################################################################################
### Electron-impact excitation (EIE) ############################################################################################


#################################################################################################################################
### Electron-impact ionization (EII) ############################################################################################


