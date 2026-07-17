

"""
`Empirical.bindingEnergy(Z::Int64, sh::Shell; data::PeriodicTable.AbstractEnergyData=XrayDataBooklet() )`  
    ... to estimate the binding) energy of the given shell for a neutral atom; these binding energies 
        are taken from semi-empirical tabulations due to data by different data sets.
        An energy::Float64 > 0. is returned. 
        ...... Perhaps, introduce an optional argument scaleLast::Bool to scale values, if no explicit binding
               energies are available.
"""
function  bindingEnergy(Z::Int64, sh::Shell; data::PeriodicTable.AbstractEnergyData=XrayDataBooklet() )
    if         typeof(data)== PeriodicTable.Williams2000
        wa = PeriodicTable.bindingEnergies_Williams2000(Z)    
    elseif     typeof(data)== PeriodicTable.Larkins1977
        wa = PeriodicTable.bindingEnergies_Larkins1977(Z) 
    elseif     typeof(data)== PeriodicTable.XrayDataBooklet
        wa = PeriodicTable.bindingEnergies_XrayDataBooklet(Z)    
    else   error("Unsupported energy data set = $data")
    end
    
    # Now extract the energies from the given subshell
    shellIndex = Dict{Shell, Int64}()
    shellIndex[Shell("1s")] =  1;    shellIndex[Shell("2s")] =  2;    shellIndex[Shell("2p")] =  4
    shellIndex[Shell("3s")] =  5;    shellIndex[Shell("3p")] =  7;    shellIndex[Shell("3d")] =  9
    shellIndex[Shell("4s")] = 10;    shellIndex[Shell("4p")] = 12;    shellIndex[Shell("4d")] = 14;    shellIndex[Shell("4f")] = 16
    shellIndex[Shell("5s")] = 17;    shellIndex[Shell("5p")] = 19;    shellIndex[Shell("5d")] = 21
    shellIndex[Shell("6s")] = 22;    shellIndex[Shell("6p")] = 24
    
    # If no useful value is found, raise an error; callers such as Empirical.scaledBindingEnergy catch it and
    # apply a Slater-screened hydrogenic estimate instead. The former behavior -- a silent sentinel of 2.0 eV --
    # poisoned every Rydberg shell (e.g. H 2s: 2.0 eV instead of the exact 3.4 eV).
    if  sh in [Shell("5g"), Shell("6d"), Shell("6f"), Shell("6g")]   error("Not defined shell $sh in data set $data")   end

    wb = wa[ shellIndex[sh] ]
    if  wb == -1.   error("No useful binding energy found for shell $sh in data set $data")   end
    wc = Defaults.convertUnits("energy: from eV to atomic", wb)

    return( wc )
end


"""
`Empirical.bindingEnergy(Z::Int64, subsh::Subshell; data::PeriodicTable.AbstractEnergyData=XrayDataBooklet() )`  
    ... to estimate the binding) energy of the given subshell for a neutral atom; these binding energies 
        are taken from semi-empirical tabulations due to data by different data sets.
        An energy::Float64 > 0. is returned.
"""
function  bindingEnergy(Z::Int64, subsh::Subshell; data::PeriodicTable.AbstractEnergyData=XrayDataBooklet() )
    if         typeof(data)== PeriodicTable.Williams2000
        wa = PeriodicTable.bindingEnergies_Williams2000(Z)    
    elseif     typeof(data)== PeriodicTable.Larkins1977
        wa = PeriodicTable.bindingEnergies_Larkins1977(Z) 
    elseif     typeof(data)== PeriodicTable.XrayDataBooklet
        wa = PeriodicTable.bindingEnergies_XrayDataBooklet(Z)    
    else   error("Unsupported energy data set = $data")
    end
    
    # Now extract the energies from the given subshell
    subshellIndex = Dict{Subshell, Int64}()
    subshellIndex[Subshell("1s_1/2")] =  1;    subshellIndex[Subshell("2s_1/2")] =  2    
    subshellIndex[Subshell("2p_1/2")] =  3;    subshellIndex[Subshell("2p_3/2")] =  4
    subshellIndex[Subshell("3s_1/2")] =  5;    subshellIndex[Subshell("3p_1/2")] =  6    
    subshellIndex[Subshell("3p_3/2")] =  7;    subshellIndex[Subshell("3d_3/2")] =  8    
    subshellIndex[Subshell("3d_5/2")] =  9;    subshellIndex[Subshell("4s_1/2")] = 10
    subshellIndex[Subshell("4p_1/2")] = 11;    subshellIndex[Subshell("4p_3/2")] = 12
    subshellIndex[Subshell("4d_3/2")] = 13;    subshellIndex[Subshell("4d_5/2")] = 14
    subshellIndex[Subshell("4f_5/2")] = 15;    subshellIndex[Subshell("4f_7/2")] = 16
    subshellIndex[Subshell("5s_1/2")] = 17;    subshellIndex[Subshell("5p_1/2")] = 18
    subshellIndex[Subshell("5p_3/2")] = 19;    subshellIndex[Subshell("5d_3/2")] = 20
    subshellIndex[Subshell("5d_5/2")] = 21;    subshellIndex[Subshell("6s_1/2")] = 22
    subshellIndex[Subshell("6p_1/2")] = 23;    subshellIndex[Subshell("6p_3/2")] = 24
    
    # If no useful value is found, issue a warning and return 2.0 eV; this can be readily improved
    wb = wa[ subshellIndex[subsh] ];       
    if  wb == -1.   wb = 2.0;  @warn("No useful binding energy found for subshell $subsh in data set $data")   end
    wc = Defaults.convertUnits("energy: from eV to atomic", wb)
    
    return( wc )
end
    
  

"""
`Empirical.bindingEnergy(Z::Int64, subsh::Subshell, conf::Configuration;
                         data::PeriodicTable.AbstractEnergyData=XrayDataBooklet())`
    ... to evaluate the binding energy of a subshell electron if part of the given configuration.
        The tabulated binding energies refer to the neutral atom; this procedure assumes that each missing
        electron (compared with the neutral system) adds 0.3 / n Hartree to the binding energy, where n is
        the principal quantum number of the subshell. It also requires that the shell associated with subsh
        is occupied in the given configuration. An energy::Float64 is returned.
"""
function  bindingEnergy(Z::Int64, subsh::Subshell, conf::Configuration;
                        data::PeriodicTable.AbstractEnergyData=XrayDataBooklet())
    sh = Shell(subsh.n, Basics.subshell_l(subsh))
    if !(sh in keys(conf.shells))   error("Shell $sh not part of $conf")   end
    wc = Empirical.bindingEnergy(Z, subsh, data=data)
    nm = Z - conf.NoElectrons
    if  nm  >  0   wc = wc + nm * 0.3 / subsh.n   end

    return( wc )
end


"""
`Empirical.ionizationPotential(Z::Int64, conf::Configuration)`  
    ... to evaluate the ionization potential I_p that is needed to ionize a valence electron from the 
        given configuration. An energy::Float64 is returned.
"""
function ionizationPotential(Z::Int64, conf::Configuration)
    if  conf != Configuration( conf.NoElectrons )
        error("Ionization potential for non-standard filling required; conf = $conf")
    end

    wa = PeriodicTable.ionizationPotentials_Nist2025(Z)
    wi = Z - conf.NoElectrons
    wb = wa[ wi + 1 ]
    wc = Defaults.convertUnits("energy: from eV to atomic", wb)
    
    return( wc )
end


"""
`Empirical.meanCharge(Z::Int64, subshell::Subshell, conf::Configuration)`
    ... to estimate the mean charge  Z as seen by an electron in subshell for an ion with the given configuation.
        This charge is obtained by inverting the hydrogenic relation bE = Z^(eff)^2 / (2 n^2) for the binding
        energy bE > 0. of the given subshell, i.e. Z^(eff) = n sqrt(2 bE). A charge::Float64 is returned.
"""
function meanCharge(Z::Int64, subshell::Subshell, conf::Configuration)
    we   = Empirical.bindingEnergy(Z, subshell, conf)
    Zeff = sqrt( 2 * subshell.n^2 * we )

    return( Zeff )
end


"""
`Empirical.totalEnergy(Z::Int64, conf::Configuration; data::PeriodicTable.AbstractEnergyData=Nist2025() )`  
    ... to estimate the total energy of the given configuration for an atom or ions with nuclear charge Z; 
        this total energy is determined either from the (negative) sum of all binding energies
        or the (negative) sum of ionization potentials that are needed to remove the remaining electrons
        from the configuration. A warning is issued if electrons occur in excited shells with regard to
        the ground configuration of the neutral atom. An energy::Float64 is returned.
        
        Note: This is an estimate of the total energy which need not to be so accurate.
"""
function  totalEnergy(Z::Int64, conf::Configuration; data::PeriodicTable.AbstractEnergyData = PeriodicTable.XrayDataBooklet() )
    if         typeof(data) == PeriodicTable.Nist2025
        wa = PeriodicTable.ionizationPotentials_Nist2025(Z) 
        # Use the ionization potentials to make a realistic estimate of the total energy ... and return the value
        # Warn if the givrn configuration does not follow the standard filling
        if  conf != Configuration( conf.NoElectrons )
            @warn("Ionization potentials for non-standard filling of shells applied; conf = $conf")
        end

        wi = Z - conf.NoElectrons
        wb = 0.;   for  ni = wi+1:length(wa)   wb = wb + wa[ni]   end
        wc = Defaults.convertUnits("energy: from eV to atomic", -wb)
    
        return( wc )
        
    elseif     typeof(data) == PeriodicTable.Williams2000
        wa = PeriodicTable.bindingEnergies_Williams2000(Z)    
    elseif     typeof(data) == PeriodicTable.Larkins1977
        wa = PeriodicTable.bindingEnergies_Larkins1977(Z) 
    elseif     typeof(data) == PeriodicTable.XrayDataBooklet
        wa = PeriodicTable.bindingEnergies_XrayDataBooklet(Z)    
    else   error("Unsupported energy data set = $data")
    end

    shellIndex = Dict{Shell, Int64}()
    shellIndex[Shell("1s")] =  1;    shellIndex[Shell("2s")] =  2;    shellIndex[Shell("2p")] =  4
    shellIndex[Shell("3s")] =  5;    shellIndex[Shell("3p")] =  7;    shellIndex[Shell("3d")] =  9
    shellIndex[Shell("4s")] = 10;    shellIndex[Shell("4p")] = 12;    shellIndex[Shell("4d")] = 14;    shellIndex[Shell("4f")] = 16
    shellIndex[Shell("5s")] = 17;    shellIndex[Shell("5p")] = 19;    shellIndex[Shell("5d")] = 21
    shellIndex[Shell("6s")] = 22;    shellIndex[Shell("6p")] = 24
    
    # Use the binding energies of the neutral atom to obtain a useful total binding energy
    totalE = 0.;   wb = 0.
    for  (sh, occ) in conf.shells
        if  sh in [Shell("5f"), Shell("5g"), Shell("6d"), Shell("6f"), Shell("6g"), Shell("6h")]  || sh.n >= 7
            wb = 1.0;  @warn("Not defined shell $sh in data set $data"); continue   
        elseif   shellIndex[sh] == -1
            wb = 2.0;  @warn("No useful binding energy found for shell $sh in data set $data")
        else 
            wb = wa[ shellIndex[sh] ]
        end
        totalE = totalE - occ * wb
    end
        
    wc = Defaults.convertUnits("energy: from eV to atomic", totalE)
    
    return( wc )
end
