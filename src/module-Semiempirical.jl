
"""
`module  JAC.Semiempirical`  
... a submodel of JAC that contains all methods to set-up and process (simple)  semiempirical estimations of 
    atomic properties that cannot be calculated so easily in a rigorious manner.
"""
module Semiempirical


using Printf,
      ..Basics, ..Defaults, ..ManyElectron, ..PeriodicTable, ..TableStrings


"""
`Semiempirical.estimate(::EstimateIonizationPotentialInnerShell, sh::Subshell, Z::Int64)`
    ... estimates the ionization potential = binding energy of an electron in the given subshell;
        values are taken from Williams et al. (2000); covers 1s_1/2 through 4p_3/2 (12 subshells);
        an energy value::Float64 in Hartree is returned.
"""
function estimate(::EstimateIonizationPotentialInnerShell, sh::Subshell, Z::Int64)
    shellIndex = Dict( Subshell("1s_1/2") => 1, Subshell("2s_1/2") => 2, Subshell("2p_1/2") => 3, Subshell("2p_3/2") => 4,
                       Subshell("3s_1/2") => 5, Subshell("3p_1/2") => 6, Subshell("3p_3/2") => 7,
                       Subshell("3d_3/2") => 8, Subshell("3d_5/2") => 9,
                       Subshell("4s_1/2") => 10, Subshell("4p_1/2") => 11, Subshell("4p_3/2") => 12 )
    idx = shellIndex[sh]

    return( Defaults.convertUnits("energy: from eV to atomic", PeriodicTable.bindingEnergies_Williams2000(Z)[idx]) )
end


"""
`Semiempirical.estimate(::EstimateBindingEnergyWilliams2000, Z::Int64, sh::Subshell)`
    ... provides the binding energy of a subshell electron from the Williams et al. (2000) tabulation,
        https://userweb.jlab.org/~gwyn/ebindene.html; an energy::Float64 in Hartree is returned.
"""
function estimate(::EstimateBindingEnergyWilliams2000, Z::Int64, sh::Subshell)
    wa = PeriodicTable.bindingEnergies_Williams2000(Z)
    if      sh == Subshell("1s_1/2")    wb = wa[1]
    elseif  sh == Subshell("2s_1/2")    wb = wa[2]
    elseif  sh == Subshell("2p_1/2")    wb = wa[3]
    elseif  sh == Subshell("2p_3/2")    wb = wa[4]
    elseif  sh == Subshell("3s_1/2")    wb = wa[5]
    elseif  sh == Subshell("3p_1/2")    wb = wa[6]
    elseif  sh == Subshell("3p_3/2")    wb = wa[7]
    elseif  sh == Subshell("3d_3/2")    wb = wa[8]
    elseif  sh == Subshell("3d_5/2")    wb = wa[9]
    elseif  sh == Subshell("4s_1/2")    wb = wa[10]
    elseif  sh == Subshell("4p_1/2")    wb = wa[11]
    elseif  sh == Subshell("4p_3/2")    wb = wa[12]
    else    error("No binding energy available for Z = $Z and subshell $sh ")
    end
    if  wb == -1.   error("No binding energy available for Z = $Z and subshell $sh ")   end

    return( Defaults.convertUnits("energy: from eV to atomic", wb) )
end


"""
`Semiempirical.estimate(::EstimateBindingEnergyLarkins1977, Z::Int64, sh::Subshell)`
    ... provides the binding energy of a subshell electron from the Larkins (1977) tabulation;
        an energy::Float64 in Hartree is returned.
"""
function estimate(::EstimateBindingEnergyLarkins1977, Z::Int64, sh::Subshell)
    wa = PeriodicTable.bindingEnergies_Larkins1977(Z)
    if      sh == Subshell("1s_1/2")    wb = wa[1]
    elseif  sh == Subshell("2s_1/2")    wb = wa[2]
    elseif  sh == Subshell("2p_1/2")    wb = wa[3]
    elseif  sh == Subshell("2p_3/2")    wb = wa[4]
    elseif  sh == Subshell("3s_1/2")    wb = wa[5]
    elseif  sh == Subshell("3p_1/2")    wb = wa[6]
    elseif  sh == Subshell("3p_3/2")    wb = wa[7]
    elseif  sh == Subshell("3d_3/2")    wb = wa[8]
    elseif  sh == Subshell("3d_5/2")    wb = wa[9]
    elseif  sh == Subshell("4s_1/2")    wb = wa[10]
    elseif  sh == Subshell("4p_1/2")    wb = wa[11]
    elseif  sh == Subshell("4p_3/2")    wb = wa[12]
    else    error("No binding energy available for Z = $Z and subshell $sh ")
    end
    if  wb == -1.   error("No binding energy available for Z = $Z and subshell $sh ")   end

    return( Defaults.convertUnits("energy: from eV to atomic", wb) )
end


"""
`Semiempirical.estimate(::EstimateBindingEnergyXrayDataBooklet, Z::Int64, sh::Subshell)`
    ... provides the binding energy of a subshell electron from the X-ray Data Booklet tabulation;
        covers 1s_1/2 through 6p_3/2 (24 subshells); an energy::Float64 in Hartree is returned.
"""
function estimate(::EstimateBindingEnergyXrayDataBooklet, Z::Int64, sh::Subshell)
    wa = PeriodicTable.bindingEnergies_XrayDataBooklet(Z)
    if      sh == Subshell("1s_1/2")    wb = wa[1]
    elseif  sh == Subshell("2s_1/2")    wb = wa[2]
    elseif  sh == Subshell("2p_1/2")    wb = wa[3]
    elseif  sh == Subshell("2p_3/2")    wb = wa[4]
    elseif  sh == Subshell("3s_1/2")    wb = wa[5]
    elseif  sh == Subshell("3p_1/2")    wb = wa[6]
    elseif  sh == Subshell("3p_3/2")    wb = wa[7]
    elseif  sh == Subshell("3d_3/2")    wb = wa[8]
    elseif  sh == Subshell("3d_5/2")    wb = wa[9]
    elseif  sh == Subshell("4s_1/2")    wb = wa[10]
    elseif  sh == Subshell("4p_1/2")    wb = wa[11]
    elseif  sh == Subshell("4p_3/2")    wb = wa[12]
    elseif  sh == Subshell("4d_3/2")    wb = wa[13]
    elseif  sh == Subshell("4d_5/2")    wb = wa[14]
    elseif  sh == Subshell("4f_5/2")    wb = wa[15]
    elseif  sh == Subshell("4f_7/2")    wb = wa[16]
    elseif  sh == Subshell("5s_1/2")    wb = wa[17]
    elseif  sh == Subshell("5p_1/2")    wb = wa[18]
    elseif  sh == Subshell("5p_3/2")    wb = wa[19]
    elseif  sh == Subshell("5d_3/2")    wb = wa[20]
    elseif  sh == Subshell("5d_5/2")    wb = wa[21]
    elseif  sh == Subshell("6s_1/2")    wb = wa[22]
    elseif  sh == Subshell("6p_1/2")    wb = wa[23]
    elseif  sh == Subshell("6p_3/2")    wb = wa[24]
    else    error("No binding energy available for Z = $Z and subshell $sh ")
    end
    if  wb == -1.   error("No binding energy available for Z = $Z and subshell $sh ")   end

    return( Defaults.convertUnits("energy: from eV to atomic", wb) )
end


"""
`Semiempirical.estimate(::EstimateBindingEnergyWilliams2000, Z::Int64, conf::Configuration)`
    ... provides an approximate total binding energy of a configuration from Williams et al. (2000);
        no relaxation effects between hole states are included; an energy::Float64 in Hartree is returned.
"""
function estimate(::EstimateBindingEnergyWilliams2000, Z::Int64, conf::Configuration)
    wa = PeriodicTable.bindingEnergies_Williams2000(Z)
    wb = 0.
    for (sh,v) in  conf.shells
        if      sh == Shell("1s")    if wa[1]  == -1.   error("stop aa")   else   wb = wb + v * wa[1]    end
        elseif  sh == Shell("2s")    if wa[2]  == -1.   error("stop ab")   else   wb = wb + v * wa[2]    end
        elseif  sh == Shell("2p")    if wa[4]  == -1.   error("stop ac")   else   wb = wb + v * wa[4]    end
        elseif  sh == Shell("3s")    if wa[5]  == -1.   error("stop ad")   else   wb = wb + v * wa[5]    end
        elseif  sh == Shell("3p")    if wa[7]  == -1.   error("stop ae")   else   wb = wb + v * wa[7]    end
        elseif  sh == Shell("3d")    if wa[9]  == -1.   error("stop af")   else   wb = wb + v * wa[9]    end
        elseif  sh == Shell("4s")    if wa[10] == -1.   error("stop ag")   else   wb = wb + v * wa[10]   end
        elseif  sh == Shell("4p")    if wa[12] == -1.   error("stop ah")   else   wb = wb + v * wa[12]   end
        end
    end

    return( Defaults.convertUnits("energy: from eV to atomic", wb) )
end


"""
`Semiempirical.estimate(::EstimateBindingEnergyLarkins1977, Z::Int64, conf::Configuration)`
    ... provides an approximate total binding energy of a configuration from Larkins (1977);
        no relaxation effects between hole states are included; an energy::Float64 in Hartree is returned.
"""
function estimate(::EstimateBindingEnergyLarkins1977, Z::Int64, conf::Configuration)
    wa = PeriodicTable.bindingEnergies_Larkins1977(Z)
    wb = 0.
    for (sh,v) in  conf.shells
        if      sh == Shell("1s")    if wa[1]  == -1.   error("stop aa")   else   wb = wb + v * wa[1]    end
        elseif  sh == Shell("2s")    if wa[2]  == -1.   error("stop ab")   else   wb = wb + v * wa[2]    end
        elseif  sh == Shell("2p")    if wa[4]  == -1.   error("stop ac")   else   wb = wb + v * wa[4]    end
        elseif  sh == Shell("3s")    if wa[5]  == -1.   error("stop ad")   else   wb = wb + v * wa[5]    end
        elseif  sh == Shell("3p")    if wa[7]  == -1.   error("stop ae")   else   wb = wb + v * wa[7]    end
        elseif  sh == Shell("3d")    if wa[9]  == -1.   error("stop af")   else   wb = wb + v * wa[9]    end
        elseif  sh == Shell("4s")    if wa[10] == -1.   error("stop ag")   else   wb = wb + v * wa[10]   end
        elseif  sh == Shell("4p")    if wa[12] == -1.   error("stop ah")   else   wb = wb + v * wa[12]   end
        end
    end

    return( Defaults.convertUnits("energy: from eV to atomic", wb) )
end


"""
`Semiempirical.estimate(::EstimateBindingEnergyXrayDataBooklet, Z::Int64, conf::Configuration)`
    ... provides an approximate total binding energy of a configuration from the X-ray Data Booklet tabulation;
        no relaxation effects between hole states are included; an energy::Float64 in Hartree is returned.
"""
function estimate(::EstimateBindingEnergyXrayDataBooklet, Z::Int64, conf::Configuration)
    wa = PeriodicTable.bindingEnergies_XrayDataBooklet(Z)
    wb = 0.
    for (sh,v) in  conf.shells
        if      sh == Shell("1s")    if wa[1]  == -1.   error("stop aa")   else   wb = wb + v * wa[1]    end
        elseif  sh == Shell("2s")    if wa[2]  == -1.   error("stop ab")   else   wb = wb + v * wa[2]    end
        elseif  sh == Shell("2p")    if wa[4]  == -1.   error("stop ac")   else   wb = wb + v * wa[4]    end
        elseif  sh == Shell("3s")    if wa[5]  == -1.   error("stop ad")   else   wb = wb + v * wa[5]    end
        elseif  sh == Shell("3p")    if wa[7]  == -1.   error("stop ae")   else   wb = wb + v * wa[7]    end
        elseif  sh == Shell("3d")    if wa[9]  == -1.   error("stop af")   else   wb = wb + v * wa[9]    end
        elseif  sh == Shell("4s")    if wa[10] == -1.   error("stop ag")   else   wb = wb + v * wa[10]   end
        elseif  sh == Shell("4p")    if wa[12] == -1.   error("stop ah")   else   wb = wb + v * wa[12]   end
        elseif  sh == Shell("4d")    if wa[14] == -1.   error("stop ai")   else   wb = wb + v * wa[14]   end
        elseif  sh == Shell("4f")    if wa[16] == -1.   error("stop aj")   else   wb = wb + v * wa[16]   end
        elseif  sh == Shell("5s")    if wa[17] == -1.   error("stop ag")   else   wb = wb + v * wa[17]   end
        elseif  sh == Shell("5p")    if wa[19] == -1.   error("stop ah")   else   wb = wb + v * wa[19]   end
        elseif  sh == Shell("5d")    if wa[21] == -1.                      else   wb = wb + v * wa[21]   end
        elseif  sh == Shell("6s")    if wa[22] == -1.                      else   wb = wb + v * wa[22]   end
        elseif  sh == Shell("6p")    if wa[24] == -1.                      else   wb = wb + v * wa[24]   end
        end
    end

    return( Defaults.convertUnits("energy: from eV to atomic", wb) )
end


"""
`Semiempirical.estimateBindingEnergies(Z::Float64, coreConf::Configuration, nRange::UnitRange{Int64})`  
    ... to estimate the binding energies of the high-n shells with n = nRange for a (multiply-charged) 
        ion with nuclear charge Z and core configuration coreConfiguration.
        A list of energies::Array{Float64} [in Hartree] is returned.
"""
function estimateBindingEnergies(Z::Float64, coreConf::Configuration, nRange::UnitRange{Int64})
    energies = Float64[];    Ne = coreConf.NoElectrons
    # Terminate if  nRange.start  <=  nMaxCore,  i.e. if any core-shell has a higher n-value
    nMaxCore = 0;     for (k,v) in coreConf.shells   if  k.n > nMaxCore   nMaxCore = k.n    end    end
    if  nRange.start  <=  nMaxCore   error("Unsupported n-values:  nMin = $(nRange.start)  <=  nMaxCore = $nMaxCore")    end
    
    # Compute an effective Zeff for each n separately (later by considering the (mean) overlap with the core-shell electrons) 
    sn = "";    sZ = "";    se = "";    su = "";    
    for  n = nRange
        effZ = Z - Ne * (1 - 0.5/n^2.2)
        eb   = -effZ^2/(2*n^2)
        push!(energies, eb );     eu = Defaults.convertUnits("energy: from atomic", eb) 
        sx = "           " * string(n);                   sn = sn * sx[end-10:end]
        sx = "           " * @sprintf("% 2.2f", effZ);    sZ = sZ * sx[end-10:end]
        sx = "           " * @sprintf("% 0.2e", eb);      se = se * sx[end-10:end]
        sx = "           " * @sprintf("% 0.2e", eu);      su = su * sx[end-10:end]
    end
    
    sa = TableStrings.inUnits("energy")
    println("  Estimated high-n binding energies for Z = $Z and the given core configuration $coreConf: \n\n" *
            "    n                   = " * sn *"\n" *
            "    Zeff                = " * sZ *"\n" * 
            "    epsilon_b [Hartree] = " * se *"\n" *
            "    epsilon_b $sa      = "  * su *"\n"    )
    
    return( energies )
end



"""
`Semiempirical.estimateBindingEnergies(Z::Float64, coreConf::Configuration, nRange::UnitRange{Int64}, l::Int64)`  
    ... to estimate the binding energies of the high-n shells with n = nRange and orbital angular momentum l
        for a (multiply-charged) ion with nuclear charge Z and core configuration coreConfiguration.
        A list of energies::Array{Float64} [in Hartree] is returned.
"""
function estimateBindingEnergies(Z::Float64, coreConf::Configuration, nRange::UnitRange{Int64}, l::Int64)
    energies = Float64[];    Ne = coreConf.NoElectrons
    # Terminate if  nRange.start  <=  nMaxCore,  i.e. if any core-shell has a higher n-value
    nMaxCore = 0;     for (k,v) in coreConf.shells   if  k.n > nMaxCore   nMaxCore = k.n    end    end
    if  nRange.start  <=  nMaxCore   error("Unsupported n-values:  nMin = $(nRange.start)  <=  nMaxCore = $nMaxCore")    end
    if  l > nRange.start - 1         error("Unsupported l-value:   l = $l  >  nMin = $(nRange.start) - 1")               end
    
    # Compute an effective Zeff for each n separately (later by considering the (mean) overlap with the core-shell electrons) 
    sn = "";    sZ = "";    se = "";    su = "";    
    for  n = nRange
        effZ = Z - Ne * (1 - 0.5/n^2.2) * (1 - 0.1 * Ne / n /  (l+1) )
        eb   = -effZ^2/(2*n^2);     eu = Defaults.convertUnits("energy: from atomic", eb)
        push!(energies, eb ) 
        sx = "           " * string(n);                   sn = sn * sx[end-10:end]
        sx = "           " * @sprintf("% 2.2f", effZ);    sZ = sZ * sx[end-10:end]
        sx = "           " * @sprintf("% 0.2e", eb);      se = se * sx[end-10:end]
        sx = "           " * @sprintf("% 0.2e", eu);      su = su * sx[end-10:end]
    end
    
    sa = TableStrings.inUnits("energy")
    println("  Estimated (n, l=$l) binding energies for Z = $Z and the given core configuration $coreConf: \n\n" *
            "    n                   = " * sn *"\n" *
            "    Zeff                = " * sZ *"\n" * 
            "    epsilon_b [Hartree] = " * se *"\n" *
            "    epsilon_b $sa      = "  * su *"\n"    )
    
    return( energies )
end



"""
`Semiempirical.estimateRadialExpectation(Z::Float64, coreConf::Configuration, nRange::UnitRange{Int64}, l::Int64)`  
    ... to estimate the radial extent (r-expectation value) of the high-n shells with n = nRange and orbital
        angular momentum l for a (multiply-charged) ion with nuclear charge Z and core configuration 
        coreConfiguration. A list of radii::Array{Float64} [in a_o] is returned.
"""
function estimateRadialExpectation(Z::Float64, coreConf::Configuration, nRange::UnitRange{Int64}, l::Int64)
    rValues = Float64[];    Ne = coreConf.NoElectrons
    # Terminate if  nRange.start  <=  nMaxCore,  i.e. if any core-shell has a higher n-value
    nMaxCore = 0;     for (k,v) in coreConf.shells   if  k.n > nMaxCore   nMaxCore = k.n    end    end
    if  nRange.start  <=  nMaxCore   error("Unsupported n-values:  nMin = $(nRange.start)  <=  nMaxCore = $nMaxCore")    end
    
    # Compute an effective Zeff for each n separately (later by considering the (mean) overlap with the core-shell electrons) 
    sn = "";    sZ = "";    sr = ""    
    for  n = nRange
        effZ = Z - Ne * (1 - 0.5/n^2.2)
        rnl  = (3*n^2 - l*(l+1)) / (2*effZ)
        push!(rValues, rnl)
        sx = "           " * string(n);                   sn = sn * sx[end-10:end]
        sx = "           " * @sprintf("% 2.2f", effZ);    sZ = sZ * sx[end-10:end]
        sx = "           " * @sprintf("% 0.2e", rnl);     sr = sr * sx[end-10:end]
    end
    
    println("  Estimated high-n <r_nl> expectation values for Z = $Z and the given core configuration $coreConf: \n\n" *
            "    n            = " * sn *"\n" *
            "    Zeff         = " * sZ *"\n" * 
            "    <r_nl> [a_o] = " * sr *"\n" )
    
    return( rValues )
end


"""
`Semiempirical.estimateSlaterZeff(Z::Float64, conf::Configuration, sh::Subshell)`
    ... estimates the effective nuclear charge Z_eff for a single electron in subshell sh,
        given the remaining electrons in conf and the nuclear charge Z.  Slater's (1930)
        empirical screening rules are applied: same-group electrons contribute σ = 0.35
        (0.30 for [1s]), electrons with n' = n−1 contribute σ = 0.85 (sp targets only),
        and all inner electrons contribute σ = 1.00.  For d and f electrons all inner
        electrons contribute σ = 1.00.  The configuration conf must represent the other
        N−1 electrons only (not include the target electron itself); a Zeff::Float64 is returned.
"""
function estimateSlaterZeff(Z::Float64, conf::Configuration, sh::Subshell)
    ##  Slater groups in ascending order (used by nested helpers below):
    ##  [1s]  [2sp] [3sp] [3d]  [4sp] [4d]  [4f]  [5sp] [5d]  [5f]  [6sp] [6d]  [7sp]
    ##    1     2     3     4     5     6     7     8     9    10    11    12    13
    function groupIndex(n::Int64, l::Int64)
        if      l <= 1    return( n <= 7 ? [1, 2, 3, 5, 8, 11, 13][n] : 13 + 2*(n - 7) )
        elseif  l == 2    return( [4, 6, 9, 12][n - 2] )
        elseif  l == 3    return( [7, 10][n - 3] )
        else    error("No Slater group defined for n = $n, l = $l")
        end
    end
    function screen(n_t::Int64, l_t::Int64, n_o::Int64, l_o::Int64)
        gi_t = groupIndex(n_t, l_t)
        gi_o = groupIndex(n_o, l_o)
        if      gi_o >  gi_t                     return( 0.0 )              ## outer group: no shielding
        elseif  gi_o == gi_t                     return( n_t == 1 ? 0.30 : 0.35 )   ## same group
        elseif  l_t <= 1  &&  n_o == n_t - 1    return( 0.85 )             ## sp target, one shell below
        else                                     return( 1.00 )             ## sp target inner / d or f target
        end
    end
    n_t   = sh.n
    l_t   = sh.kappa > 0 ? sh.kappa : -sh.kappa - 1
    sigma = 0.0
    for (shell, occ) in conf.shells
        sigma = sigma + occ * screen(n_t, l_t, shell.n, shell.l)
    end
    return( Z - sigma )
end


end # module
