
using  QuadGK, ..Basics,  ..Defaults


"""
`struct  Plasma.IonicLevel`  
    ... a struct to comprise the (necessary) information about a single ionic level in Saha-Boltzmann computations.

    + energy          ::Float64        ... total energy E (Z, A, q) of the ionic level.
    + g               ::Float64        ... degeneracy g (Z, A, q) of the ionic level.
    + nDensity        ::Float64        ... number density of ions I(Z, A, q) in this level.
    + uppermostShell  ::Shell          ... uppermost shell in the underlying configuration.
"""
struct   IonicLevel
    energy            ::Float64
    g                 ::Float64
    nDensity          ::Float64
    uppermostShell    ::Shell  
end  


# `Base.string(lev::Plasma.IonicLevel)`  ... provides a String notation for the variable lev::Plasma.IonicLevel.
function Base.string(lev::Plasma.IonicLevel)
    sa = "Ionic level (energy=$(lev.energy), g=$(lev.g), nDensity=$(lev.nDensity), uppermostShell=$(lev.uppermostShell)):"
    return( sa )
end


# `Base.show(io::IO, lev::Plasma.IonicLevel)`  ... prepares a proper printout of the lev::Plasma.IonicLevel.
function Base.show(io::IO, lev::Plasma.IonicLevel)
    sa = Base.string(lev);                   print(io, sa)
end


"""
`struct  Plasma.IonicClass`  
    ... a struct to comprise the (necessary) information about a single ionic charge state in 
        Saha-Boltzmann computations.

    + q               ::Int64                       ... charge state q of the ions.
    + groundEnergy    ::Float64                     ... ground level energy of the ion.
    + ionLevels       ::Array{Plasma.IonicLevel,1}  ... ionic levels of this charge state
"""
struct   IonicClass
    q                 ::Int64 
    groundEnergy      ::Float64  
    ionLevels         ::Array{Plasma.IonicLevel,1}
end  


# `Base.string(ion::Plasma.IonicClass)`  ... provides a String notation for the variable ion::Plasma.IonicClass.
function Base.string(ion::Plasma.IonicClass)
    sa = "Ionic class (q=$(ion.q), groundEnergy=$(ion.groundEnergy), No of ionic levels=$(length(ion.ionLevels))):"
    return( sa )
end


# `Base.show(io::IO, ion::Plasma.IonicClass)`  ... prepares a proper printout of the ion::Plasma.IonicClass.
function Base.show(io::IO, ion::Plasma.IonicClass)
    sa = Base.string(ion);                   print(io, sa)
end


"""
`struct  Plasma.IsotopeClass`  
    ... a struct to comprise the (necessary) information about a single isotopic fraction in Saha-Boltzmann computations.

    + isotopicDensity    ::Float64                    ... (number) density n (Z,A) of all ions of this isotope class.
    + Lambda             ::Float64                    ... thermal length Lambda of the given isotope.
    + dominantEnergy     ::Float64                    
        ... energy of the "dominant" ion level of this isotope; this depends on energy and help improve convergence
    + isotopicFraction   ::Basics.IsotopicFraction    ... isotopic fraction in the Saha-Boltzmann mixture.
    + ionClasses         ::Array{Plasma.IonicClass,1} ... list of all ions of the isotope I(Z,A).
"""
struct   IsotopeClass
    isotopicDensity      ::Float64
    Lambda               ::Float64 
    dominantEnergy       ::Float64                    
    isotopicFraction     ::Basics.IsotopicFraction 
    ionClasses           ::Array{Plasma.IonicClass,1}
end  


# `Base.string(isotope::Plasma.IsotopeClass)`  ... provides a String notation for the variable isotope::Plasma.IsotopeClass.
function Base.string(isotope::Plasma.IsotopeClass)
    sa = "Isotope class (Z=$(isotope.isotopicFraction.Z), A=$(isotope.isotopicFraction.A)):"
    return( sa )
end


# `Base.show(io::IO, isotope::Plasma.IsotopeClass)`  ... prepares a proper printout of the isotope::Plasma.IsotopeClass.
function Base.show(io::IO, isotope::Plasma.IsotopeClass)
    sa = Base.string(isotope);             print(io, sa)
    println(io, "isotopicDensity:          $(isotope.isotopicDensity)  ")
    println(io, "Lambda:                   $(isotope.Lambda)  ")     
    println(io, "dominantEnergy:           $(isotope.dominantEnergy)  ")     
    println(io, "isotopicFraction:         $(isotope.isotopicFraction)  ")
    println(io, "ionClasses:               $(isotope.ionClasses)  ")
end

#################################################################################################################################
#################################################################################################################################


"""
`Plasma.computeDominantIsotopeEnergy(isoClass::IsotopeClass)`  
    ... computes the (dominant) energy of all ions of the given isotope class, i.e. an energy at which
        ionic levels contributes most to the overall mixture; this (dominant) energy is utilized to make 
        the computation more stable; an energy::Float64 is returned.   ... Not used in the present code.
"""
function computeDominantIsotopeEnergy(isoClass::IsotopeClass)
    energy = 0.;   nDensity = 0.
    for  ionClass  in  isoClass.ionClasses
        for  level  in  ionClass.ionLevels
            if  level.nDensity > nDensity    nDensity = level.nDensity;    energy = level.energy    end
        end
    end
    
    return ( energy )
end


"""
`Plasma.computeElectronChemicalPotential(temp::Float64, ne::Float64)`  
    ... computes the electron chemical potential for the given electron (number) density; a chemMu::Float64 is returned.
"""
function computeElectronChemicalPotential(temp::Float64, ne::Float64)
    eLambda = sqrt(2pi) / sqrt(temp)
    chemMu  = - temp * log( 2.0 / (ne * eLambda^3) )
    
    return ( chemMu )
end


"""
`Plasma.computeElectronNumberDensity(temp::Float64, isoClasses::Array{IsotopeClass,1})`  
    ... computes the electron (number) density n_e (T) from the density of all ionic levels, charge states
        and isotopic mixture at temperature T; an nDensity::Float64 is returned.
"""
function computeElectronNumberDensity(temp::Float64, isoClasses::Array{IsotopeClass,1})
    nDensity = 0.;
    for  isoClass  in  isoClasses
        for  ionClass  in  isoClass.ionClasses
            for  level  in  ionClass.ionLevels
                nDensity = nDensity + level.nDensity * ionClass.q
            end
        end
    end
        
    return ( nDensity )
end


"""
`Plasma.computeIonicLevelChemicalPotential(temp::Float64, ionLevel::IonicLevel, Lambda::Float64)`  
    ... computes the ionic-level chemical potential for the given level (number) density; 
        a chemMu::Float64 is returned.   ... Not used in the present code.
"""
function computeIonicLevelChemicalPotential(temp::Float64, ionLevel::IonicLevel, Lambda::Float64)
    chemMu = - temp * log( ionLevel.g / (ionLevel.nDensity * Lambda^3) ) + ionLevel.energy
    
    return ( chemMu )
end


"""
`Plasma.computeIsotopicPartitionFunction(temp::Float64, isoClass::IsotopeClass, chemMuE::Float64, dominantEnergy::Float64)`  
    ... computes the (unnormalized) partition function of the given isotopic class in terms of a summation 
        over all charge states and ionic levels of these class. A pf::Float64 is returned.
"""
function computeIsotopicPartitionFunction(temp::Float64, isoClass::IsotopeClass, chemMuE::Float64, dominantEnergy::Float64)
    beta = 1 / temp;    pf = 0.
    for  ionClass in isoClass.ionClasses
        for level in ionClass.ionLevels
            pf = pf + level.g * exp( -beta * (chemMuE * ionClass.q + level.energy - dominantEnergy) ) 
        end
    end
        
    return ( pf )
end

#==  20.7.24
"""
`Plasma.computeIsotopicNumberDensity(temp::Float64, isoClass::IsotopeClass)`  
    ... computes the (number) density n (Z,A) of a given isotope class at temperature T;
        an nDensity::Float64 is returned.   ... Not used in the present code.
"""
function computeIsotopicNumberDensity(temp::Float64, isoClass::IsotopeClass)
    nDensity = 0.
    for  ionClass  in  isoClass.ionClasses
        for  level  in  ionClass.ionLevels      nDensity = nDensity + level.nDensity    end
    end
        
    return ( isoClass.isotopicFraction.x * nDensity )
end  ==#


"""
`Plasma.computeIonLevelNumberDensity(temp::Float64, q::Int64, ionLevel::IonicLevel, pfIsoClass::Float64,
                                     isoDensity::Float64, chemMuE::Float64, dominantEnergy::Float64)`  
    ... computes the (number) density of an individual ionLevel of ions of charge state q at temperature T, 
        the electron chemical potential chemMuE and for a given partion function of the associated
        IsotopeClass pfIonClass; an nDensity::Float64 is returned.
"""
function computeIonLevelNumberDensity(temp::Float64, q::Int64, ionLevel::IonicLevel, pfIsoClass::Float64,
                                      isoDensity::Float64, chemMuE::Float64, dominantEnergy::Float64)
    beta     = 1 / temp
    nDensity = isoDensity * ionLevel.g * exp( -beta * (chemMuE * q + ionLevel.energy - dominantEnergy) ) / pfIsoClass
    return ( nDensity )
end


"""
`Plasma.computeIonLevelNumberDensityTotal(isoClass::IsotopeClass)`  
    ... computes the total (number) density of all charge states and individual ionLevel's just by summation.
        an nDensity::Float64 is returned.
"""
function computeIonLevelNumberDensityTotal(isoClass::IsotopeClass)
    nDensity = 0.
    for  ionClass in isoClass.ionClasses
        for level in ionClass.ionLevels    nDensity = nDensity + level.nDensity    end
    end
    return ( nDensity )
end


"""
`Plasma.computeMeanIsotopeChargeState(isoClass::IsotopeClass)`  
    ... computes the (mean) charge states <q> of all ions of the given isotope; a qbar::Float64 is returned.
"""
function computeMeanIsotopeChargeState(isoClass::IsotopeClass)
    qnDensity = 0.
    for  ionClass  in  isoClass.ionClasses
        for  level  in  ionClass.ionLevels    qnDensity = qnDensity + level.nDensity * ionClass.q    end
    end
     
    return ( qnDensity / isoClass.isotopicDensity )
end


"""
`Plasma.determineIsotopeClasses(scheme::Plasma.SahaBoltzmannScheme, temp::Float64)`  
    ... determines which isotope classes which need to be considered for the given Saha-Boltzmann mixture;
        it first normalizes all isotopic fractions to 1 and then selects the relevant charge states in the 
        mixture due to the given temperature and the number of charge states involved the self-consistent 
        treatment of the mixture. However, this procedure does neither determine the number nor the detailed 
        properties of the ionic levels of these ion classes. These ionic levels and all further properties are 
        first left out (empty) for later specification. A list of isoClasses::Array{Plasma.IsotopicClass,1}
"""
function determineIsotopeClasses(scheme::Plasma.SahaBoltzmannScheme, temp::Float64)
    isoClasses = Plasma.IsotopeClass[];    wa = 0.
    newIsoFractions = Basics.IsotopicFraction[]
    # Re-normalize the ionic mixture to make sure the total ionic and isotopic densities agree
    for  isoFraction in scheme.isotopicMixture    wa = wa + isoFraction.x     end   # Determine the normalization
    for  isoFraction in scheme.isotopicMixture    
        push!( newIsoFractions,  Basics.IsotopicFraction(isoFraction.Z, isoFraction.A, isoFraction.x/wa) )
    end
    
    # Determine the relevant ionic classes for each isotope fraction of the mixture
    for  isoFraction in newIsoFractions 
        ionClasses = Plasma.determineIonicClasses(scheme, temp, isoFraction)
        Lambda     = Plasma.thermalLength(temp, isoFraction)
        push!(isoClasses, Plasma.IsotopeClass(0., Lambda, 0., isoFraction, ionClasses) )
    end
    
    return ( isoClasses )
end


"""
`Plasma.determineIonicClasses(scheme::Plasma.SahaBoltzmannScheme, temp::Float64, isoFraction::Basics.IsotopicFraction)`  
    ... determines all relevant ionic classes of the given isotope with fraction x at temperature temp 
        which (sufficiently) contribute to the ionic mixture. The procedure selects the relevant charge states 
        but does neither determine the number nor the properties of the ionic levels of these ion classes. 
        These levels need to be later initialized by either explicit computations or reading them from a file.
        A list of ionClasses::Array{Plasma.IonicClass,1} is returned.
"""
function determineIonicClasses(scheme::Plasma.SahaBoltzmannScheme, temp::Float64, isoFraction::Basics.IsotopicFraction)
    ionClasses = Plasma.IonicClass[];    nZ = round(Int, isoFraction.Z);    nA = round(Int, isoFraction.A);    
    closestNe = 0;   deltaIp = 1.0e6     
    # Determine the charge state with an ionization potential closest to kT; 
    for  ne = 1:nZ
        groundConf = Plasma.determineReferenceConfiguration(ne)
        subshells  = Basics.extractRelativisticSubshellList([groundConf])
        nMax = 0;    for  subsh in subshells    if  subsh.n  > nMax  nMax = subsh.n                                          end   end
        kMax = 0;    for  subsh in subshells    if  subsh.n == nMax  &&  abs(subsh.kappa) > abs(kMax)   kMax = subsh.kappa   end   end
        bindingEn  = Semiempirical.estimate("binding energy", nZ, Subshell(nMax,kMax), useLarkins=true )
        if  abs(bindingEn - temp)  <  deltaIp   deltaIp = abs(bindingEn - temp);   closestNe = ne   end
    end
    
    neMin = max(0, nZ-scheme.qRange.stop);   neMax = nZ-scheme.qRange.start
    for  ne = neMin:neMax
        if  ne > nZ   continue    end
        push!(ionClasses, Plasma.IonicClass(nZ-ne, -1000., Plasma.IonicLevel[]) )
    end
    
    println(">> Selected ionic charge states for isotope  I(Z=$nZ,A=$nA)^q+  with  q = $(scheme.qRange)  and  " * 
            "ne = $(neMin:neMax)   ... number of electrons")
    
    return ( ionClasses )
end


"""
`Plasma.determineInitialIonDensitiesPropterties(scheme::Plasma.SahaBoltzmannScheme, temp::Float64, ionDensity::Float64,
                                                isoClass::IsotopeClass)`  
    ... determines the initial level densities and all temperature-dependent properties which are not yet defined
        for solving the Saha-Boltzmann equations. The procedure assumes however that the ionic level information has been
        provided before. A list of isoClasses::Array{Plasma.IsotopeClass,1} is returned.
"""
function determineInitialIonDensitiesPropterties(scheme::Plasma.SahaBoltzmannScheme, temp::Float64, ionDensity::Float64, 
                                                 isoClass::IsotopeClass)
    # Determine the mean charge state in the mixture
    nq = 0;     qbar = 0.
    for  ionClass in isoClass.ionClasses    nq = nq + 1;   qbar = qbar + ionClass.q     end
    qbar            = qbar / nq
    newIonClasses   = Plasma.IonicClass[]
    ne              = qbar * ionDensity                                  # Electron number density; assume qbar-fold ionized ions
    chemMuE         = Plasma.computeElectronChemicalPotential(temp, ne)  # Electron chemical potential

    # Initial the isotope class data
    isotopicDensity = ionDensity * isoClass.isotopicFraction.x
    Lambda          = Plasma.thermalLength(temp, isoClass.isotopicFraction)
    dominantEnergy  = -4*temp   ## Plasma.computeDominantIsotopeEnergy(isoClass)  -32.3703 
    
    # Compute factor to re-normalize the partition function
    nDensityTotal = Plasma.computeIonLevelNumberDensityTotal(isoClass)
    pfIsoClass    = Plasma.computeIsotopicPartitionFunction(temp, isoClass, chemMuE, dominantEnergy) 
    pfIsoClass    = pfIsoClass * isotopicDensity / nDensityTotal
    
    for  ionClass  in  isoClass.ionClasses
        # Initialize the ion level data; first compute the (re-normalized) partition function of the ionic class
        groundEnergy = 0.
        newIonLevels = Plasma.IonicLevel[]
        for  ionLevel  in  ionClass.ionLevels
            if   ionLevel.energy < groundEnergy     groundEnergy = ionLevel.energy     end
            nDensity = Plasma.computeIonLevelNumberDensity(temp, ionClass.q, ionLevel, pfIsoClass, 
                                                           isotopicDensity, chemMuE, dominantEnergy)
            push!(newIonLevels, Plasma.IonicLevel(ionLevel.energy, ionLevel.g, nDensity, ionLevel.uppermostShell) )
        end
        # Initialize the ion class data
        push!(newIonClasses, Plasma.IonicClass(ionClass.q, groundEnergy, newIonLevels) )
    end
    newIsoClass     = Plasma.IsotopeClass(isotopicDensity, Lambda, dominantEnergy, isoClass.isotopicFraction, newIonClasses)
        
    return( newIsoClass )
end


"""
`Plasma.determineIonNumberDensity(rho::Float64, mixture::Array{IsotopicFraction,1})`  
    ... determines the ion number density for a given matter density [g/cm^3] and mixture of isotopic fractions.
        It re-normalizes the isotopic fractions and performs the conversion by taking the mass of ^12C as the basis.
        A ion number density ni [ions/a_o^3] is returned.
"""
function determineIonNumberDensity(rho::Float64, mixture::Array{IsotopicFraction,1})
    xTotal = 0.;   Axtotal = 0.
    for  isoFraction in mixture    xTotal = xTotal + isoFraction.x;    Axtotal = Axtotal + isoFraction.A * isoFraction.x   end 
    factor = (0.529_177_210_67e-8)^3 * 12 / 1.66054e-24
    ni     = rho * factor * xTotal / Axtotal
    return( ni )
end


"""
`Plasma.determineIpShifts(pm::Basics.AbstractPlasmaModel, temp::Float64, ni::Float64, isotopeClass::IsotopeClass)`  
    ... determines the ground-state energie as well as the shifts Delta I_p (q --> q+1) for all charge states q = 0 ... qmax 
        of the ions defined by isotopeClass. The procedure assumes that the plasma model pm brings in the current plasma
        parameters. In some models, the ion number density ni is required as well.
        Two dictionaries groundE[q]::Dict{Int64,Float64}, deltaIp[q]::Dict{Int64,Float64} 
        are returned with all the energies (shifts) in [Hartree].
"""
function determineIpShifts(pm::Basics.AbstractPlasmaModel, temp::Float64, ni::Float64, isoClass::IsotopeClass)
    groundE = Dict{Int64,Float64}();     deltaIp = Dict{Int64,Float64}()
    
    # Determine the ground-state energies of the different charge states
    for ionClass  in  isoClass.ionClasses   groundE[ionClass.q] = ionClass.groundEnergy   end
    
    # Determine the ionization potential shifts that depend on the given plasma model
    if      typeof(pm) == Basics.NoPlasmaModel
        for ionClass  in  isoClass.ionClasses   deltaIp[ionClass.q] = -1000.   end
    elseif  typeof(pm) == Basics.WithoutAutoionizationModel
        for ionClass  in  isoClass.ionClasses   deltaIp[ionClass.q] = 0.       end
    elseif  typeof(pm) == Basics.DebyeHueckelModel
        for ionClass  in  isoClass.ionClasses   
            deltaIp[ionClass.q] = - (ionClass.q + 1)^2 / pm.debyeLength    
        end
    elseif  typeof(pm) == Basics.IonSphereModel 
        for ionClass  in  isoClass.ionClasses
            R = ( 3 * ionClass.q / (4pi * pm.electronDensity) )^(1/3)
            deltaIp[ionClass.q] = - 3/2 * ionClass.q / R
        end
    elseif  typeof(pm) == Basics.StewartPyattModel
        for ionClass  in  isoClass.ionClasses
            R  = ( 3 * ionClass.q / (4pi * pm.electronDensity) )^(1/3)
            wa = ( (R / pm.lambda)^3 + 1. )^(2/3)
            wb = - 1 / (2temp) / (pm.electronDensity / ni + 1.) 
            deltaIp[ionClass.q] = wb * (wa - 1.) 
        end
    else
        error("Undefined plasma model.")
    end
        
    return( groundE, deltaIp )
end



"""
`Plasma.determineReferenceConfiguration(ne::Int64)`  
    ... determines the reference configuration for an ion with ne electrons. A refConfig::Configuration is returned.
"""
function determineReferenceConfiguration(ne::Int64)
    if ne > 20  error("Re-check this procedure ...")    end
    
    # These reference configurations are independent of Z
    if      ne == 1         refConfig = Configuration("1s")
    elseif  ne == 2         refConfig = Configuration("1s^2")
    elseif  ne == 3         refConfig = Configuration("1s^2 2s")
    elseif  ne == 4         refConfig = Configuration("1s^2 2s^2")
    elseif  3 <= ne < 11    refConfig = Configuration("1s^2 2s^2 2p^$(ne-4)")
    elseif  ne == 11        refConfig = Configuration("1s^2 2s^2 2p^6 3s")
    elseif  ne == 12        refConfig = Configuration("1s^2 2s^2 2p^6 3s^2")
    elseif  13 <= ne < 19   refConfig = Configuration("1s^2 2s^2 2p^6 3s^2 3p^$(ne-12)")
    elseif  ne == 19        refConfig = Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 4s")
    elseif  ne == 20        refConfig = Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 4s^2")
    else    error("stop a")
    end
    
    return ( refConfig )
end


"""
`Plasma.determineValenceShells(ne::Int64)`  
    ... determines the valance shells for an ion with ne electrons. A shellList::Array{Shell,1} is returned.
"""
function determineValenceShells(ne::Int64)
    if ne > 20  error("Re-check this procedure ...")    end
    
    # These reference configurations are independent of Z
    if      1 <= ne <=  2    shellList = Shell[]
    elseif  3 <= ne <=  4    shellList = [Shell("2s")]
    elseif  5 <= ne <=  6    shellList = [Shell("2s"), Shell("2p")]
    elseif  7 <= ne <= 10    shellList = [Shell("2p")]
    else    error("stop a")
    end
    
    return ( shellList )
end


"""
`Plasma.displayIsotopeClasses(stream::IO, isoClasses::Array{IsotopeClass,1})`  
    ... list all basic information about the ionic mixture, i.e. the isotopic fractions, the thermal length, 
        the isotopic density as well as the number of charge states and ionic levels involved into the computations.
        A neat table is printed but nothing is returned otherwise.
"""
function  displayIsotopeClasses(stream::IO, isoClasses::Array{IsotopeClass,1})
    nx = 105
    println(stream, " ")
    println(stream, "  Isotope classes of the ionic mixture:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(22, "Isotope (Z, A, x)"   ; na=2);                       sb = sb * TableStrings.hBlank(24)
    sa = sa * TableStrings.center(18, "Lambda (T,A)"        ; na=2);                       sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(10, "n (Z,A)"             ; na=2);                       sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "Dominant E"          ; na=2);                       
    sb = sb * TableStrings.center(10, " [Hartree]"          ; na=2)
    sa = sa * TableStrings.center(10, "No (ions)"           ; na=2);                       sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(16, "q    No (levels)"    ; na=2);                       sb = sb * TableStrings.hBlank(18)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  isoClass in isoClasses
        sa  = "";    sc  = "                 ";    
        sc = sc * string(round(isoClass.isotopicFraction.Z, digits=2)) * ", "
        sc = sc * string(round(isoClass.isotopicFraction.A, digits=2)) * ", "
        sc = sc * string(round(isoClass.isotopicFraction.x, digits=2))
        sa = sa * sc[end-20:end]                                * "          "
        sa = sa * @sprintf("%.2e", isoClass.Lambda)             * "        "
        sa = sa * @sprintf("%.2e", isoClass.isotopicDensity)    * "    "
        sa = sa * @sprintf("%.4e", isoClass.dominantEnergy)     * "   "
        sc = "      " * string(length(isoClass.ionClasses))
        sa = sa * sc[end-4:end] * "   "    
        sb = repeat(" ", length(sa))
        println(stream, sa)
        for  ionClass  in  isoClass.ionClasses   
            sc = "      " * string(ionClass.q);                     sb = sb * sc[end-4:end] * "  "; 
            sc = "      " * string(length(ionClass.ionLevels));     sb = sb * sc[end-6:end]
            println(stream, sb);    
            sb = repeat(" ", length(sa))   
        end
    end
    println(stream, "  ", TableStrings.hLine(nx)) 
     
    return( nothing )
end


"""
`Plasma.displaySahaBoltzmannEquilibrium(stream::IO, isoClasses::Array{IsotopeClass,1}, 
                                        pm::Basics.AbstractPlasmaModel; printLevelDensities::Bool=false)`  
    ... list some basic information about the ionic mixture, along with the charge states and the ionic densities;
        it also prints the ionic level number densities if this is required.
        A neat table is printed but nothing is returned otherwise.
"""
function  displaySahaBoltzmannEquilibrium(stream::IO, isoClasses::Array{IsotopeClass,1}, 
                                          pm::Basics.AbstractPlasmaModel; printLevelDensities::Bool=true)
    nx = 115
    println(stream, " ")
    println(stream, "  Saha-Boltzmann equilibrium of the ionic mixture:")
    println(stream, " ")
    println(stream, "    + Plasma model used for IPD: $(typeof(pm)) ")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  ";   sb = "  "
    sa = sa * TableStrings.center(22, "Isotope (Z, A, x)"   ; na=2);                       sb = sb * TableStrings.hBlank(24)
    sa = sa * TableStrings.center(18, "Lambda (T,A)"        ; na=2);                       sb = sb * TableStrings.hBlank(20)
    sa = sa * TableStrings.center(10, "n (Z,A)"             ; na=2);                       sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(10, "Dominant E"          ; na=2);                       
    sb = sb * TableStrings.center(10, " [Hartree]"          ; na=2)
    sa = sa * TableStrings.center(10, "No (ions)"           ; na=2);                       sb = sb * TableStrings.hBlank(12)
    sa = sa * TableStrings.center(30, "q    Energy     g    nDensity"    ; na=2);             sb = sb * TableStrings.hBlank(18)
    println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx)) 
    #   
    for  isoClass in isoClasses
        sa  = "";    sc  = "                 ";    
        sc = sc * string(round(isoClass.isotopicFraction.Z, digits=2)) * ", "
        sc = sc * string(round(isoClass.isotopicFraction.A, digits=2)) * ", "
        sc = sc * string(round(isoClass.isotopicFraction.x, digits=2))
        sa = sa * sc[end-20:end]                                * "         "
        sa = sa * @sprintf("% .2e", isoClass.Lambda)             * "       "
        sa = sa * @sprintf("% .2e", isoClass.isotopicDensity)    * "   "
        sa = sa * @sprintf("% .4e", isoClass.dominantEnergy)     * "   "
        sc = "      " * string(length(isoClass.ionClasses))
        sa = sa * sc[end-4:end] * "   "    
        sb = repeat(" ", length(sa))
        println(stream, sa)
        for  ionClass  in  isoClass.ionClasses   
            sc = "         " * string(ionClass.q);         sb = sb * sc[end-4:end] * "  " 
            println(stream, sb);    
            sb = repeat(" ", length(sa))
            sd = repeat(" ", length(sb)) * "       "   
            if  printLevelDensities
                for  ionLevel in ionClass.ionLevels
                    se = sd *  @sprintf("% .3e", ionLevel.energy) * " " * @sprintf("% 3.0f", ionLevel.g) * 
                                                                    " " * @sprintf("% .3e",  ionLevel.nDensity)
                    println(stream, se)
                end
            end
        end
    end
    println(stream, "  ", TableStrings.hLine(nx)) 

    # Display charge-state distributions for all given isotopes in the mixture; loop through all isotope-classes
    nx = 95
    println(stream, " ")
    println(stream, "  Charge-state distributions of the ionic mixture:")
    println(stream, " ")
    for  isoClass  in  isoClasses
        sa = "  ";   sb = "  ";    sc  = "   ";    
        sc = sc * string(round(isoClass.isotopicFraction.Z, digits=2)) * ", "
        sc = sc * string(round(isoClass.isotopicFraction.A, digits=2)) * ", "
        sc = sc * string(round(isoClass.isotopicFraction.x, digits=2))
        sa = sa * TableStrings.center(22, "Isotope (Z, A, x)"   ; na=2);                  
        sb = sb * sc
        sa = sa * "  q =   ";                                                            
        sb = sb * "        |    "
        for  ionClass  in  isoClass.ionClasses
            ntot = 0.
            for level in ionClass.ionLevels   ntot = ntot + level.nDensity    end;     
            ntot = ntot / isoClass.isotopicDensity
            sa   = sa * TableStrings.center(10, string(ionClass.q);       na=2);                  
            sb   = sb * TableStrings.center(10, @sprintf("% .3e", ntot) ; na=2);  
        end
        println(stream, sa);    println(stream, "  ", TableStrings.hLine(nx));  println(stream, sb) 
        println(stream, " ")
    end
     
    return( nothing )
end


"""
`Plasma.displayThermodynamicProperties(stream::IO, temp::Float64, totalIonDensity::Float64, isoClasses::Array{IsotopeClass,1})`  
    ... evaluates, converts units and displays the thermodynamic properties of the ionic mixture, based on the given 
        isotopic-class data. A neat table is printed but nothing is returned otherwise.
"""
function  displayThermodynamicProperties(stream::IO, temp::Float64, totalIonDensity::Float64, isoClasses::Array{IsotopeClass,1})
    #
    ne         = Plasma.computeElectronNumberDensity(temp, isoClasses)
    betamu     = Plasma.computeElectronChemicalPotential(temp, ne) / temp
    meanCharge = Plasma.computeElectronNumberDensity(temp, isoClasses) / totalIonDensity
    pressure   = Plasma.pressure(temp, totalIonDensity, isoClasses)
    internalEn = Plasma.internalEnergyDensity(temp, totalIonDensity, isoClasses)
    freeEn     = Plasma.freeEnergyDensity(temp, totalIonDensity, isoClasses)
    #
    println(stream, "\n  Thermodynamic properties of the ionic mixture: \n")
    println(stream, "  + betamu:                  " * @sprintf("% .4e", betamu)                  )
    println(stream, "  + mean charge state <q>:   " * @sprintf("% .4e", meanCharge)              )
    println(stream, "  + pressure:                " * @sprintf("% .4e", pressure)     * " [a.u.]")
    println(stream, "  + internal energy:         " * @sprintf("% .4e", internalEn)   * " [a.u.]")
    println(stream, "  + free-energy:             " * @sprintf("% .4e", freeEn)       * " [a.u.]")
     
    return( nothing )
end



"""
`Plasma.freeEnergyDensity(temp::Float64, totalIonDensity::Float64, isoClasses::Array{IsotopeClass,1})`  
    ... computes the internal energy density u(T) of an isotopic mixture at temperature T; a 
        u::Float64 is returned.
"""
function freeEnergyDensity(temp::Float64, totalIonDensity::Float64, isoClasses::Array{IsotopeClass,1})
    wf = 0.;   beta = 1.0/temp
    
    # Make use of the mean charge
    for  isoClass  in  isoClasses
        Lambda = isoClass.Lambda 
        for  ionClass in isoClass.ionClasses
            for  level in ionClass.ionLevels 
                br = level.g / (level.nDensity * Lambda^3)
                wf = wf - level.nDensity * temp * (log(br) - beta * level.energy + 1.0)   
            end
        end
    end
    # Add the contributions from the free-electron density
    ne = Plasma.computeElectronNumberDensity(temp, isoClasses) 
    wf = wf - ne * temp * (log(2.0/ne) + 1.0)        

    println(">>>> Free energy density f(T) = $wf")
        
    return ( wf )
end
    

"""
`Plasma.internalEnergyDensity(temp::Float64, totalIonDensity::Float64, isoClasses::Array{IsotopeClass,1})`  
    ... computes the internal energy density u(T) of an isotopic mixture at temperature T; a 
        u::Float64 is returned.
"""
function internalEnergyDensity(temp::Float64, totalIonDensity::Float64, isoClasses::Array{IsotopeClass,1})
    wu = 0.;   qbar = 0.;   wa = 0.
    
    # Make use of the mean charge
    for  isoClass  in  isoClasses
        qbar = isoClass.isotopicFraction.x * Plasma.computeMeanIsotopeChargeState(isoClass)
        for  ionClass in isoClass.ionClasses
            for  level in ionClass.ionLevels    wa = wa + level.nDensity * level.energy   end
        end
    end
    wu = 3/2 * totalIonDensity * temp * (1. + qbar) + wa
    
    println(">>>> Internal energy density u(T) = $wu")
        
    return ( wu )
end


"""
`Plasma.perform(scheme::Plasma.SahaBoltzmannScheme, computation::Plasma.Computation; output::Bool=true)`  
    ... to perform a Saha-Boltzmann equilibrium computation for a given ion mixture. For output=true, a dictionary 
        is returned from which the relevant results can be can easily accessed by proper keys.
"""
function  perform(scheme::Plasma.SahaBoltzmannScheme, computation::Plasma.Computation; output::Bool=true)
    if  output    results = Dict{String, Any}()    else    results = nothing    end
        
    wa         = 0.;    for  fraction  in  scheme.isotopicMixture     wa = wa + fraction.x    end
    isoMixture = IsotopicFraction[]; 
    for  fraction  in  scheme.isotopicMixture   push!(isoMixture, IsotopicFraction(fraction.Z, fraction.A, fraction.x/wa) )   end
    #
    if  computation.settings.useNumberDensity    totalIonDensity = computation.settings.density
    else                                         error("stop a")
    end
    
    # Return results if required
    if  output   
        results["temperature"]    = computation.settings.temperature               
        results["number density"] = totalIonDensity               
        results["(normalized) isotope mixture"] = isoMixture                 
    end
    println(" ")

    # Determine the isotope and ion classes for which a Saha-Boltzmann LTE need to be considered
    isoClasses = Plasma.determineIsotopeClasses(scheme, computation.settings.temperature)
    # Read or generate ion-level data for each isotope class and initialize the number densities
    newIsoClasses = IsotopeClass[]
    for  isoClass  in  isoClasses
        newIsoClass = Plasma.readEvaluateIonLevelData(scheme, isoClass, computation.grid)
        newIsoClass = Plasma.restrictIonLevelData(scheme, newIsoClass) 
        newIsoClass = Plasma.determineInitialIonDensitiesPropterties(scheme, computation.settings.temperature, totalIonDensity, newIsoClass)
        push!(newIsoClasses, newIsoClass)  
    end

    Plasma.displayIsotopeClasses(stdout, newIsoClasses)
    # Solve for the equilibrium densities
    newIsoClasses = Plasma.solveSahaBoltzmannEquilibrium(computation.settings.temperature, newIsoClasses, scheme.plasmaModel)
    # Display equilibrium densities and thermodynamic properties of the given mixture
    if  scheme.printIonLevels     Plasma.displaySahaBoltzmannEquilibrium(stdout, newIsoClasses, scheme.plasmaModel)   end 
    Plasma.displayThermodynamicProperties(stdout, computation.settings.temperature, totalIonDensity, newIsoClasses)   
    
    println("\n> Saha-Boltzmann LTE computation complete ...")
    
    Defaults.warn(PrintWarnings())
    Defaults.warn(ResetWarnings())
    return( results )
end


"""
`Plasma.pressure(temp::Float64, totalIonDensity::Float64, isoClasses::Array{IsotopeClass,1})`  
    ... computes the pressure of an isotopic mixture at temperature T; a pressure::Float64 is returned.
"""
function pressure(temp::Float64, totalIonDensity::Float64, isoClasses::Array{IsotopeClass,1})
    pressure = 0.;      wp = 0.;   qbar = 0.
    # Compute (if desired) the pressure from the level population explicitly ... for comparison
    if true
        for  isoClass  in  isoClasses
            for  ionClass in isoClass.ionClasses
                for  level in ionClass.ionLevels    wp = wp + level.nDensity * temp   end
            end
        end
        wp = wp + Plasma.computeElectronNumberDensity(temp, isoClasses) * temp
    end
    
    # Make use of the mean charge
    for  isoClass  in  isoClasses
        qbar = qbar + isoClass.isotopicFraction.x * Plasma.computeMeanIsotopeChargeState(isoClass)
    end
    pressure = totalIonDensity * temp * (1. + qbar)
    
    println(">>>> Pressure = $pressure  =!= wp = $wp")
        
    return ( pressure )
end


"""
`Plasma.solveSahaBoltzmannEquilibrium(temp::Float64, isoClasses::Array{IsotopeClass,1}, pm::Basics.AbstractPlasmaModel)`  
    ... solves the Saha-Boltzmann equations for a given ionic mixture and temperature; it iterates
        the ion-level number densities until no relevant changes occur. 
        A newIsoClasses::Array{IsotopeClass,1} with the proper number densities is returned.
"""
function solveSahaBoltzmannEquilibrium(temp::Float64, isoClasses::Array{IsotopeClass,1}, pm::Basics.AbstractPlasmaModel)
    noIteration   = 0;   isNotConverged = true;     newIsoClasses = Plasma.IsotopeClass[]
    oldIsoClasses = deepcopy(isoClasses);           betamuOld     = 1.0;                   opm = pm
    ionDensity    = 0. 
    for isoClass  in  isoClasses    ionDensity = ionDensity + isoClass.isotopicDensity * isoClass.isotopicFraction.x    end 
    
    while  isNotConverged
        noIteration = noIteration + 1
        println("\n> Saha-Boltzmann iteration $noIteration ... \n  ------------------------------")
        #
        if  noIteration == 1   ne = 2.0e-5 * 6     ## ion number density * Z
        else                   ne = Plasma.computeElectronNumberDensity(temp, oldIsoClasses)
        end
        chemMuE       = Plasma.computeElectronChemicalPotential(temp, ne) 
        
        # Update plasma model for ionization-potential depression (IPD)
        npm = Plasma.updateIpdPlasmaModel(opm, temp, ne, oldIsoClasses)
        
        # Update all isotopic classes in the ionic mixture
        newIsoClasses = Plasma.IsotopeClass[]
        for  isoClass  in  oldIsoClasses
            # Extract the ground energies and Delta I_p for the given plasma model
            groundE, deltaIp = Plasma.determineIpShifts(npm, temp, ionDensity, isoClass)
            qMax = 0;    for ionClass  in  isoClass.ionClasses   qMax = max(qMax, ionClass.q)   end
            
            dominantEnergy = -4*temp  ## Plasma.computeDominantIsotopeEnergy(isoClass)  -32.3703
                
            # Compute factor to re-normalize the partition function
            pfIsoClass    = Plasma.computeIsotopicPartitionFunction(temp, isoClass, chemMuE, dominantEnergy) 

            newIonClasses  = Plasma.IonicClass[]
            for  ionClass  in  isoClass.ionClasses
                groundEnergy = 0.
                newIonLevels = Plasma.IonicLevel[]
                for  ionLevel  in  ionClass.ionLevels
                    if   ionLevel.energy < groundEnergy     groundEnergy = ionLevel.energy     end
                    # Assign a non-zero level density only if the level is bound under the given IPD plasma model
                    if  ionLevel.energy == 0.                                          ||
                        ionClass.q      == qMax                                        ||
                        ionLevel.energy + deltaIp[ionClass.q] < groundE[ionClass.q+1]
                           nDensity = Plasma.computeIonLevelNumberDensity(temp, ionClass.q, ionLevel, pfIsoClass, 
                                                                          isoClass.isotopicDensity, chemMuE, dominantEnergy)
                    else   nDensity = 0.
                    end
                    push!(newIonLevels, Plasma.IonicLevel(ionLevel.energy, ionLevel.g, nDensity, ionLevel.uppermostShell) )
                end
                push!(newIonClasses, Plasma.IonicClass(ionClass.q, groundEnergy, newIonLevels) )
            end
            newIsoClass = Plasma.IsotopeClass(isoClass.isotopicDensity, isoClass.Lambda, dominantEnergy, 
                                              isoClass.isotopicFraction, newIonClasses)
            push!(newIsoClasses, newIsoClass)
                
            
            
        end
        # Determine and report about convergence
        ne     = Plasma.computeElectronNumberDensity(temp, newIsoClasses)
        betamu = Plasma.computeElectronChemicalPotential(temp, ne) / temp
        if  abs(betamu - betamuOld) < 1.0e-10    isNotConverged = false   else  betamuOld = betamu   end 
        println(">> betamu = $betamu ")
        #
        if     isNotConverged  &&  noIteration < 20   oldIsoClasses = deepcopy(newIsoClasses);   opm = npm
        else   break
        end
    end
    
    return( newIsoClasses )
    
end


"""
`Plasma.thermalLength(temp::Float64, isoFraction::Basics.IsotopicFraction)`  
    ... returns the thermal length Lambda(T, M) for a particle with mass M at temperature T.
        A length Lambda::Float64 in {a_o] is returned.
"""
function thermalLength(temp::Float64, isoFraction::Basics.IsotopicFraction)
    M = isoFraction.A / Defaults.ELECTRON_MASS_U
    Lambda = sqrt(2pi) / sqrt(M * temp)
    
    return ( Lambda )
end


#################################################################################################################################
#################################################################################################################################


"""
`Plasma.generateIonLevelData(scheme::Plasma.SahaBoltzmannScheme, isoClass::IsotopeClass, q::Int64, grid::Radial.Grid)`  
    ... generates the ionic-level data for the given isotope class and charge state q of the ionic mixture.
        These ionic-level data are generated for the requested number in scheme.NoIonLevels and 
        excitations in scheme.NoExcitations of levels. The procedure presently assumes that the electron
        shells in the reference configuration are filled in standard order with the Z-q electrons. 
        A special treatment of the reference configuration need to be introduced if this is not the case. 
        The energy levels are generated by including scheme.NoExcitations excitations w.r.t. the reference 
        configuration, and all levels with an excitation  energy < maxEn are taken into account, up
        to scheme.NoIonLevels.         
        An ionLevels::Array{IonicLevels,1} is returned where all ionic-level data are properly places
        but all (number) densities are still set to 0. 
"""
function generateIonLevelData(scheme::Plasma.SahaBoltzmannScheme, isoClass::IsotopeClass, q::Int64, grid::Radial.Grid)
    ionLevels = IonicLevel[];  Zint = round(Int, isoClass.isotopicFraction.Z)
    
    # Return a single ionic level for bare ions
    if  Zint - q == 0   return( [Plasma.IonicLevel(0., 1., 0., Shell(0,0))] )   end
        
    # Determine the reference and the associated configurations for an atomic computation
    refConfig  = Plasma.determineReferenceConfiguration( Zint - q )
    valShells  = Plasma.determineValenceShells( Zint - q )
    
    # Check the grid for a proper presentation of all orbitals
    ## Radial.checkGrid(grid, Shell(scheme.upperShellNo, 0), isoClass.isotopicFraction.Z)
    @warn("Check the grid in Plasma: re-install !!")
    
    # The generation of the ionic level data requires special care, especially if double or multiple 
    # excitations and excitations into high-n shells will be considered. In practice, this requires a 
    # specific procedure (abstract tpye) to distinguish between core and valence electrons and the 
    # procedure how the wave functions for individual ionic levels will be optimized.
    
    # At present, the following  selections are made 
    if true
        println(">> Rules for generating the ionic level data:                                       \n" *
                "   (1) Single excitations from all shells of the reference configuration $refConfig \n" * 
                "   (2) Double & higher excitations only from the valence shells          $valShells \n" *
                "   (3) Generate a mean-field basis (orbitals) for all configurations together;      \n" *
                "   (4) Generate the ionic levels separately for each configuration from the list below.")
    end
    
    # Generate from refConfig a set of configurations in line with the given scheme
    fromShells = Basics.extractShellList([refConfig])
    toShells   = Basics.generateShellList([refConfig], scheme.upperShellNo, [0,1,2,3])
    #
    sConfigs   = Configuration[];    dtConfigs   = Configuration[]
    if  scheme.NoExcitations >= 1    sConfigs   = Basics.generateConfigurations([refConfig], fromShells, toShells, 1)  end
    if  scheme.NoExcitations >  1    
        dtConfigs  = Basics.generateConfigurations([refConfig], valShells,  toShells, scheme.NoExcitations) 
    end
    allConfigs = append!(sConfigs, dtConfigs)
    allConfigs = unique(allConfigs)
    #
    println(">> Generation of ion-level data for Z = $Zint and charge state q = $q  with configurations \n   $allConfigs")
    
    if  length(valShells) == 0      repConfigs = allConfigs
    else                            repConfigs = dtConfigs
    end
    
    
    # Generate a mean-field multiplet for these configurations and compute all atomic levels in Dirac-Coulomb approximation
    sZ           = trunc(isoClass.isotopicFraction.Z, digits=2);
    nm           = Nuclear.Model(isoClass.isotopicFraction.Z, isoClass.isotopicFraction.A)
    name         = "Mean multiplet for Z=$sZ and q=$q"
    mfSettings   = AtomicState.MeanFieldSettings()
    repBasis     = Representation(name, nm, grid, repConfigs, MeanFieldBasis(mfSettings) )
    repOutput    = generate(repBasis, output=true)
    
    orbitals = repOutput["mean-field basis"].orbitals
    
    
    # Extract levels from output
    repMultiplet = Multiplet()
    for  conf in allConfigs
        mp           = Basics.performCI([conf], orbitals, nm, grid, AsfSettings(), printout=false)
        repMultiplet = Multiplet("append", append!(repMultiplet.levels, mp.levels) )
    end
    Basics.displayLevels(stdout, [repMultiplet], N=100)
    
    for  level in repMultiplet.levels
        # Determine the uppermost shell of the given level
        confs          = Basics.extractConfigurations(Basics.FromBasis(), level.basis)
        shellList      = Basics.extractNonrelativisticShellList(confs)
        uppermostShell = Shell(0,0)
        for  shell in shellList     if  shell > uppermostShell    uppermostShell = shell    end     end
            
        push!(ionLevels, Plasma.IonicLevel(level.energy, Basics.twice(level.J)+1, 0., uppermostShell) )
    end
    #
    println("\n >> Generation of $(length(ionLevels)) ion-levels for Z = $Zint and charge state q = $q.")

    # Print information to summury file
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if printSummary
        println(iostream, "\n\n>> Generation of ion-level data for Z = $Zint and charge state q = $q" *
                          " with $(length(allConfigs)) configurations.\n")
        for  conf in allConfigs   println(iostream, "   $conf ")     end
        println(iostream, "\n>> SCF field based on shells: \n" *
                          "   valShells  = $valShells       \n" *
                          "   fromShells = $fromShells       \n" *
                          "   toShells   = $toShells       \n")
        Basics.displayLevels(iostream, [repMultiplet], N=100)
        println(iostream, "\n>> Generation of $(length(ionLevels)) ion-levels for Z = $Zint and charge state q = $q.")
    end 
    
    return( ionLevels )
end


"""
`Plasma.readEvaluateIonLevelData(filename::String, isoClass::IsotopeClass, NoExcitations::Int64, upperShellNo::Int64)`  
    ... reads in, if available, the ionic-level data for the given isotope class from filename.
        The ionic-level data are accepted for return, if (1) they belong to isotope (Z,A) in isoClasses
        and if they fullfill (2) upperShellNo <= filename:upperShellNo, (3) NoExcitations <= filename:NoExcitations.
        If proper data are found, an updated newIsoClass::IsotopeClass returned and missing otherwise.
"""
function readEvaluateIonLevelData(filename::String, isoClass::IsotopeClass, NoExcitations::Int64, upperShellNo::Int64)
    found = false;  chargeStates = Int64[];   newIonClasses = Plasma.IonicClass[]
    # Open the file and try to read in the directory
    data    = JLD2.load(filename)
    isoData = data["isoData"]

    @show trunc(isoData["Z"], digits=3), trunc(isoClass.isotopicFraction.Z, digits=3) 
    @show trunc(isoData["A"], digits=3), trunc(isoClass.isotopicFraction.A, digits=3)
    @show trunc(isoData["Z"], digits=3)  ==  trunc(isoClass.isotopicFraction.Z, digits=3)   &&
        trunc(isoData["A"], digits=3)  ==  trunc(isoClass.isotopicFraction.A, digits=3)   &&
        isoData["NoExcitations"]       >=  NoExcitations                                  &&                      
        isoData["upperShellNo"]        >=  upperShellNo
    
    if  trunc(isoData["Z"], digits=3)  ==  trunc(isoClass.isotopicFraction.Z, digits=3)   &&
        trunc(isoData["A"], digits=3)  ==  trunc(isoClass.isotopicFraction.A, digits=3)   &&
        isoData["NoExcitations"]       >=  NoExcitations                                  &&                      
        isoData["upperShellNo"]        >=  upperShellNo  
        # The given data look suitable to be used; now check for ion-level data for all requested charge states
        for  ionClass in isoClass.ionClasses
            qkey = "q" * string(ionClass.q) * "Levels"
            if   haskey(isoData, qkey)  push!(newIonClasses, isoData[qkey] )
                                        push!(chargeStates,  ionClass.q)
            else                        push!(newIonClasses, ionClass)
            end
        end
        #
        newIsoClass = Plasma.IsotopeClass(isoClass.isotopicDensity, isoClass.Lambda, isoClass.dominantEnergy, 
                                          isoClass.isotopicFraction, newIonClasses)
        found       = true
        println(">>> Ion-level data are found in $filename  for $(chargeStates)+." )
        #
    elseif  false  # Use this branch to read-in external data that were generated for test purposes
        for  ionClass in isoClass.ionClasses
            qkey = "q" * string(ionClass.q) * "Levels"
            if   haskey(isoData, qkey)  push!(newIonClasses, isoData[qkey] )
                                        push!(chargeStates,  ionClass.q)
            else                        push!(newIonClasses, ionClass)
            end
        end
        #
        newIsoClass = Plasma.IsotopeClass(isoClass.isotopicDensity, isoClass.Lambda, isoClass.dominantEnergy, 
                                          isoClass.isotopicFraction, newIonClasses)
        found       = true
        println(">>> Ion-level data are found in $filename  for $(chargeStates)+. \n" )
        #
    else
        println(">>> No ion-level data found in $filename." )
    end
     
    if  found   return( newIsoClass )   else    return( missing )     end
end


"""
`Plasma.readEvaluateIonLevelData(scheme::Plasma.SahaBoltzmannScheme, isoClass::IsotopeClass, grid::Radial.Grid)`  
    ... reads in or evaluates the ionic-level data for the given isotope class of the ionic mixture.
        This steps basically generates/provides all atomic level data upon which Saha-Boltzmann ionic mixture 
        is based on. To generate these data it analyses of whether ion-level data are provided by one of the 
        files in scheme.isotopeFilenames::{String,1}; if this is the case, it further analyzes of whether 
        the provided data fullfill the requested number scheme.NoIonLevels and excitations scheme.NoExcitations
        of levels for each charge state of the isotope. If no suitable ion-level data are found for an
        isotope (Z,A) in isoClasses, the routine generates a (new) set of ion-level data for the current 
        computation. In addition, a (new) file with name NewIonLevelsZxxxAyyyy.... is written to disk 
        and can be utilized in any subsequent computation. It is typically assumed that these file are renamed 
        with the proper physical settings for their generation in mind.
        
        A newIsoClass::IsotopeClass is returned where all ionic-level data are properly places
        but all (number) densities are still set to 0. All other fields remain unchanged.
"""
function readEvaluateIonLevelData(scheme::Plasma.SahaBoltzmannScheme, isoClass::IsotopeClass, grid::Radial.Grid)
    
    writeData = false;   incomplete = false;   newIsoClass = missing;  
    for  filename in scheme.isotopeFilenames
        newIsoClass = readEvaluateIonLevelData(filename, isoClass, scheme.NoExcitations, scheme.upperShellNo)
        if  !(typeof(newIsoClass) == Missing)   break   end
    end
    
    # Analyze whether newIsoClass has a complete set of ionic level data and generate new data, if incomplete
    if  typeof(newIsoClass) == Missing         
        incomplete = true;    newIsoClass = deepcopy(isoClass)                             
    else
        for ionClass in newIsoClass.ionClasses
            if  length(ionClass.ionLevels) == 0    incomplete = true    end
        end
    end
    #
    if incomplete
        newIonClasses = Plasma.IonicClass[] 
        for  ionClass in newIsoClass.ionClasses
            if  length(ionClass.ionLevels) > 0
                push!(newIonClasses, ionClass)
            else
                println("\n\n>>> Generate new ion-level data for Z = $(trunc(isoClass.isotopicFraction.Z, digits=3)), " *
                                                                "A = $(trunc(isoClass.isotopicFraction.A, digits=3)) and " *
                        "charge state $(ionClass.q) ... ")
                writeData = true
                ionLevels = Plasma.generateIonLevelData(scheme, newIsoClass, ionClass.q, grid)
                push!(newIonClasses, Plasma.IonicClass(ionClass.q, ionClass.groundEnergy, ionLevels) )
            end
        end
        newIsoClass = Plasma.IsotopeClass( newIsoClass.isotopicDensity, newIsoClass.Lambda, newIsoClass.dominantEnergy,
                                           newIsoClass.isotopicFraction, newIonClasses)
    end

    # Write out newIsoClass if needed
    if  writeData
        sa = string(round(Int, isoClass.isotopicFraction.Z));    sb = string(round(Int, isoClass.isotopicFraction.A))
        filename = "newIonicLevelDataZ" * sa * "A" * sb * "-" * string(Dates.now())[1:13] * ".jld"
        Plasma.writeIonLevelData(filename, newIsoClass, scheme.NoExcitations, scheme.upperShellNo)
    end
    
    return ( newIsoClass )
end


"""
`Plasma.readUpdateIonLevelDataObsolete(filename::String, isoClass::IsotopeClass, NoExcitations::Int64, upperShellNo::Int64,
                               newIonClasses::Array{Plasma.IonicClass,1})`  
    ... reads in, if available, the ionic-level data for the given isotope class from filename.
        The ionic-level data are accepted for return, if (1) they belong to isotope (Z,A) in isoClasses
        and if they fullfill (2) upperShellNo <= filename:upperShellNo, (3) NoExcitations <= filename:NoExcitations.
        The newIonClasses are updated/appended in the given ionClasses data and a new ionic-level data file is written
        out. The procedure terminates with an error if no proper ionic level data are found in filename.
        Nothing is returned. ... This procedures is not used a present.
"""
function readUpdateIonLevelDataObsolete(filename::String, isoClass::IsotopeClass, noExcitations::Int64, upperShellNo::Int64,
                                newIonClasses::Array{Plasma.IonicClass,1})
    found = false;  chargeStates = Int64[];    updatedIonClasses = Plasma.IonicClass[]
        
    # Open the file and try to read in the directory
    data = JLD2.load(filename)
    @show data
    
    if  trunc(isoData["Z"], digits=3)  ==  trunc(isoClass.isotopicFraction.Z, digits=3)   &&
        trunc(isoData["A"], digits=3)  ==  trunc(isoClass.isotopicFraction.A, digits=3)   &&
        isoData["NoExcitations"]       >=  NoExcitations                                  &&                      
        isoData["upperShellNo"]        >=  upperShellNo  
        
        # Update or append all new data
        for  (qkey,v)  in  newIonLevelData    isoData[qkey] = v   end
            
        # The given data look suitable to be used; now construct a newIsoClass to be printed out
        for  ionClass in isoClass.ionClasses
            qkey = "q" * string(ionClass.q) * "Levels"
            if   haskey(isoData, qkey)     push!(updatedIonClasses, IonicClass(ionClass.q, ionClass.groundEnergy, isoData[qkey]) )
            else                           push!(updatedIonClasses, ionClass)
                println(">> No ionic-level data are found for q = $(ionClass.q) ... ")
            end
        end
        #
        newIsoClass = Plasma.IsotopeClass( newIsoClass.isotopicDensity, newIsoClass.Lambda, newIsoClass.dominantEnergy,
                                           newIsoClass.isotopicFraction, updatedIonClasses)
    end

    # Write out newIsoClass if needed
    if  writeData
        sa = "   " * round(Int, isoClass.isotopicFraction.Z);    sb = "   " * round(Int, isoClass.isotopicFraction.A)
        filename = "newIonicLevelDataZ" * sa[end-1:end] * sa[end-2:end]
        Plasma.writeIonLevelData(filename, newIsoClass)
        println(">>> New ionic-level data are printed to file $filename.")
    end
    
    return ( nothing )
end


"""
`Plasma.restrictIonLevelData(scheme::Plasma.SahaBoltzmannScheme, isoClass::IsotopeClass)`  
    ... restricts the compiled ion-level-data tp those levels that need to be included into the Saha-Boltzmann
        equilibrium iteration. a newIsoClass::IsotopeClass is returned.
"""
function restrictIonLevelData(scheme::Plasma.SahaBoltzmannScheme, isoClass::IsotopeClass)
    println(">>> Restrict ion-level data to uppermost n <= $(scheme.upperShellNo) ... ")
    newIonClasses = Plasma.IonicClass[]
    for  ionClass  in  isoClass.ionClasses
        newIonLevels = Plasma.IonicLevel[]
        for  level  in  ionClass.ionLevels
            if  level.uppermostShell.n <= scheme.upperShellNo    push!(newIonLevels, level)    end
        end
        push!(newIonClasses, Plasma.IonicClass(ionClass.q, ionClass.groundEnergy, newIonLevels))
    end
    newIsoClass = Plasma.IsotopeClass(isoClass.isotopicDensity, isoClass.Lambda, isoClass.dominantEnergy, 
                                      isoClass.isotopicFraction, newIonClasses)
    
    return ( newIsoClass )
end


"""
`Plasma.updateIpdPlasmaModel(pm::Basics.AbstractPlasmaModel, temp::Float64, ne::Float64, isotopeClasses::Array{IsotopeClass,1})`  
    ... updates the (parameters) of the plasma model pm due to the given temperature temp, electron density ne and the 
        level population of the mixture. A new plasma model npm::type(pm) is returned.
"""
function updateIpdPlasmaModel(pm::Basics.AbstractPlasmaModel, temp::Float64, ne::Float64, isotopeClasses::Array{IsotopeClass,1})
    
    # Determine the ionization potential shifts that depend on the given plasma model
    if      typeof(pm) == Basics.NoPlasmaModel                    npm = pm
    elseif  typeof(pm) == Basics.WithoutAutoionizationModel       npm = pm
    elseif  typeof(pm) == Basics.DebyeHueckelModel
        nt = ne
        for  isotopeClass in isotopeClasses
            for ionClass  in  isotopeClass.ionClasses
                for  level in ionClass.ionLevels
                    nt = nt + level.nDensity * ionClass.q^2
                end 
            end 
        end 
        debyeLength = sqrt( nt/ temp )
        npm = Basics.DebyeHueckelModel(debyeLength, pm.debyeRadius)
    elseif  typeof(pm) == Basics.IonSphereModel 
        npm = Basics.IonSphereModel(pm.radius, ne)
    elseif  typeof(pm) == Basics.StewartPyattModel
        ni = 0.;    for isoClass  in  isotopeClasses    ni = ni + isoClass.isotopicDensity * isoClass.isotopicFraction.x   end 
        lambda = 4pi / temp * ne * (ne / ni + 1.)
        lambda = 1 / sqrt(lambda)
        npm = Basics.StewartPyattModel(pm.radius, ne, lambda)
    else
        error("Undefined plasma model.")
    end
        
    return( npm )
end


"""
`Plasma.writeIonLevelData(filename::String, isoClass::IsotopeClass, NoExcitations::Int64, upperShellNo::Int64)`  
    ... writes-out the ionic-level data for the given isotope class to filename. Nothing is returned.
"""
function writeIonLevelData(filename::String, isoClass::IsotopeClass, NoExcitations::Int64, upperShellNo::Int64)
    # Generate a proper dictionary to be printed out
    isoData = Dict{String, Any}()
    isoData["Z"]              =  trunc(isoClass.isotopicFraction.Z, digits=3)
    isoData["A"]              =  trunc(isoClass.isotopicFraction.A, digits=3)
    isoData["NoExcitations"]  =  NoExcitations            
    isoData["upperShellNo"]   =  upperShellNo  
    # The given data look suitable to be used; now construct a newIsoClass to be printed out
    for  ionClass in isoClass.ionClasses
        qkey = "q" * string(ionClass.q) * "Levels"
        isoData[qkey]         = ionClass
    end
    
    # Open the file and dump the directory
    println(">>> Save ion level data to file  $filename")
    JLD2.@save filename isoData
    
    return ( nothing )
end


"""
`Plasma.writeIonLevelDataRobin(filenameRobin::String, filenameJac::String)`  
    ... reads-in ion-level data from filenameRobin (in tabular form) and writes a proper isoData file
        for the given isotope class to filenameJac. Nothing is returned.
"""
function writeIonLevelDataRobin(filenameRobin::String, filenameJac::String)
    # Generate a proper dictionary to be printed out; Z and A must be commuicated and set explicitly
    wZ = 6.0;   wA = 12.2
    isoData = Dict{String, Any}();     ionClasses = Plasma.IonicClass[]
    isoData["Z"]              =  trunc(wZ, digits=3)
    isoData["A"]              =  trunc(wA, digits=3)
    isoData["NoExcitations"]  =  0            
    isoData["upperShellNo"]   =  6  
    # The given data look suitable to be used; now construct a newIsoClass to be printed out
    f = open(filenameRobin)
    readline(f);    readline(f)  # Jump over the header
    
    sLines = String[]
    while  !eof(f)
        push!(sLines, readline(f))
    end
    @show  length(sLines), sLines[1]
    
    # Test that all data are read-in properly
    #== for  sLine  in sLines
        nele = Base.parse(Int64,   strip(sLine[17:20]) );   q = round(Int, wZ) - nele
        en   = Base.parse(Float64, strip(sLine[33:52]) )
        g    = Base.parse(Int64,   strip(sLine[54:57]) )
        @show nele, q, en, g
    end
    error("Check the test !!")  ==#

    # Determine all q in the data set
    qs = Int64[]
    for  sLine  in sLines
        nele = Base.parse(Int64,   strip(sLine[17:20]) );   q = round(Int, wZ) - nele
        push!(qs, q)
    end
    qs = unique(qs)
    println(">> Charge state in the input file are: $qs ")
    
    # Create an IonClass for each q
    for  q  in qs
        ionLevels = Plasma.IonicLevel[];  groundEnergy = 1.0e6
        for  sLine  in sLines
            nele = Base.parse(Int64,   strip(sLine[17:20]) );   qx = round(Int, wZ) - nele
            if  qx == q
                en   = Base.parse(Float64, strip(sLine[33:52]) )
                g    = Base.parse(Int64,   strip(sLine[54:57]) )
                push!(ionLevels, Plasma.IonicLevel(en, g, 0.) )
                if  en < groundEnergy   groundEnergy = en   end
            end
        end
        push!(ionClasses, Plasma.IonicClass(q, groundEnergy, ionLevels) )
    end
    println(">> Number of ionic classes in the input file is: $(length(ionClasses)) ")
    
    for  ionClass in ionClasses
        qkey = "q" * string(ionClass.q) * "Levels"
        isoData[qkey]         = ionClass
    end
    
    # Open the file and dump the directory
    println(">> Save ion level data to file  $filenameJac")
    JLD2.@save filenameJac isoData
    
    return ( nothing )
end

