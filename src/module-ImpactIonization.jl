
#== August 2025, the following replacements need to be made and tested properly:
++ module-ImpactIonization.jl:   newConfs  = Basics.generateConfigurationsWithElectronLoss(confs::Array{Configuration,1}, [valenceShell])
++ Basics.generateConfigurationsWithElectronLoss(confs::Array{Configuration,1}, fromShells::Array{Shell,1}) 
   ... generates a list of non-relativistic configurations for the given (reference) confs and with one 
       removed (ionized) electron fromShells.
++ See Basics.generateConfigurations(RemoveElectrons(), confs)
==#


"""
`module  JAC.ImpactIonization`  
... a submodel of JAC that contains all methods for computing electron impact ionization cross sections.
    This module was contributed by Yuancheng Wang (2023) and especially also by Bowen Li (2025).
"""
module ImpactIonization


using  Printf, ..Basics, ..Cubature, 
       ..Defaults, ..ManyElectron, ..Nuclear, ..PeriodicTable, ..Radial, ..RadialIntegrals, ..TableStrings


"""
`abstract type ImpactIonization.AbstractModel` 
    ... defines an abstract and a number of singleton types for the computation of empirical EII cross sections in different 
        models.
    
    + struct BEBmodel        
        ... to apply a generalized BEB model for non-relativistic impact energies due to
            Kim et al. PRA (2001) & Wang  et al. JPB (2024).
    + struct BEDmodel        
        ... to apply a improved BED model due to Huo PRA (2001) for non-relativistic impact energies.  
    + struct RelativisticBEBmodel        
        ... to apply a generalized BEB model for non-relativistic impact energies due to
            Kim et al. PRA (2001) & Wang  et al. JPB (2024).
    + struct RelativisticBEDmodel
        ... to apply a modified relativistic BEB model due to Uddin et al. PRA (2005) for relativistic impact energies.
    + struct FittedBEDmodel 
        ... to apply a modified relativistic BEB model due to Haque et al. AQC (2017) for non-relativistic & 
            relativistic impact energies, and by using fitting coefficients for K, L, and M subshells
            (also known as MUIBED). 
    + struct DirectMultipleModel        
        ... to apply a generalized DirectMultiple model for electron-impact multiple ionization
            Belenger et al. JPB (1997).
    + struct InDirectDoubleModel        
        ... to apply a generalized InDirectDouble model for electron-impact double ionization
            Shevelko et al. JPB (2005).
    + struct DoubleExperimentModel        
        ... to apply a improved InDirectDouble model based on experimental data fitting for electron-impact double ionization
            Shevelko et al. JPB (2006).
    + struct LotzMultipleModel
        ... a generalized LotzMultiple model for electron-impact multiple ionization
            Hahn et al. ApJ (2015).
"""
abstract type  AbstractModel  
end

struct         BEBmodel                  <:  AbstractModel   end
struct         BEDmodel                  <:  AbstractModel   end
struct         RelativisticBEBmodel      <:  AbstractModel   end
struct         RelativisticBEDmodel      <:  AbstractModel   end
struct         FittedBEDmodel            <:  AbstractModel   end
struct         DirectMultipleModel       <:  AbstractModel   end   # modified for EIMI
struct         InDirectDoubleModel       <:  AbstractModel   end   # modified for EIMI
struct         DoubleExperimentModel     <:  AbstractModel   end   # modified for EIMI
struct         LotzMultipleModel         <:  AbstractModel   end   # modified for EIMI

export BEBmodel, BEDmodel, RelativisticBEBmodel, RelativisticBEDmodel, FittedBEDmodel,
       DirectMultipleModel, InDirectDoubleModel, DoubleExperimentModel, LotzMultipleModel


"""
`struct  ImpactIonization.Settings  <:  AbstractEmpiricalSettings`  
    ... defines a type for the settings of empirical electron-impact ionization cross sections.
    
    + model                      ::ImpactIonization.AbstractModel 
        ... model for computing empirical EII cross sections.
    + multipleN                  ::Int64  # modified for EIMI
        ... N of multiply charge   # modified for EIMI
    + impactEnergies             ::Array{Float64,1}                       
        ... List of impact-energies of the incoming elecgtrons.
    ## + electronEnergies           ::Array{Float64,1}   ... List of continuum-electron energies.
    ## + calcDifferentialCs         ::Bool                                   
    ##     ... True if energy-differential cross sections are to be calculated, and false otherwise.
    + calcPartialCs              ::Bool                                   
        ... True if partial (shell-dependent) cross sections are to be calculated, and false otherwise.
    + calcTotalCs                ::Bool                                   
        ... True if total cross sections are to be calculated, and false otherwise.
    + shellSelection             ::ShellSelection      ... Describe the selected shells, if any.                               ... 
"""    
struct Settings   <:  AbstractEmpiricalSettings
    model                        ::ImpactIonization.AbstractModel
    multipleN                    ::Int64   # modified for EIMI
    impactEnergies               ::Array{Float64,1}
    calcPartialCs                ::Bool 
    calcTotalCs                  ::Bool 
    shellSelection               ::ShellSelection
end      


"""
`JAC.ImpactIonization.Settings()`  ... constructor for the default values of empirical electron-impact ionization cross sections.
"""
function Settings()
    Settings(ImpactIonization.AbstractModel(), Int64[], Float64[], false, false, ShellSelection())  # Int64 for multipleN added, modified for EIMI
end


## `Base.show(io::IO, settings::ImpactIonization.Settings)`  
##  ... prepares a proper printout of the variable settings::ImpactIonization.Settings.
function Base.show(io::IO, settings::ImpactIonization.Settings) 
    println(io, "model:                       $(settings.model)  ")
    println(io, "multipleN:                   $(settings.multipleN)  ")  # modified for EIMI
    println(io, "impactEnergies:              $(settings.impactEnergies)  ")
    println(io, "calcPartialCs:               $(settings.calcPartialCs)  ")
    println(io, "calcTotalCs:                 $(settings.calcTotalCs)  ")
    println(io, "shellSelection:              $(settings.shellSelection)  ")     
end


"""
`struct  ImpactIonization.CrossSection`  
    ... defines a type for dealing with (shell-dependent) partial EII cross sections that are based on empirical computations.

    + subshell          ::Subshell          ... subshell for which the (differential and partial) cross sections are calculated.
    + impactEnergy      ::Float64           ... energy of impact electron [Hartree].
    + partialCS         ::Float64           ... partial cross sections for this shell.
    + electronEnergies  ::Array{Float64,1}  ... energies of the ionized electron.
    + differentialCS    ::Array{Float64,1}  ... differential cross sections for this shell.
"""
struct  CrossSection
    subshell            ::Subshell
    impactEnergy        ::Float64
    partialCS           ::Float64
    electronEnergies    ::Array{Float64,1}
    differentialCS      ::Array{Float64,1}
end 

struct  MultipleCrossSection
    impactEnergy        ::Float64
    totalCS             ::Float64
end 


# `Base.show(io::IO, cs::ImpactIonization.CrossSection)`  ... prepares a proper printout of the variable line::ImpactIonization.CrossSection.
function Base.show(io::IO, cs::ImpactIonization.CrossSection) 
    println(io, "subshell:           $(cs.subshell)  ")
    println(io, "impactEnergy:       $(cs.impactEnergy)  ")
    println(io, "partialCS:          $(cs.partialCS)  ") 
    println(io, "electronEnergies:   $(cs.electronEnergies)  ")
    println(io, "differentialCS:     $(cs.differentialCS)  ")
end

    
"""
`ImpactIonization.computeCrossSections(basis::Basis, grid::Radial.Grid, nm::Nuclear.Model, settings::ImpactIonization.Settings; output::Bool=true)`  
    ... to compute the EII cross sections for an atom or ion in chosen model. While the basis describes the
        configuration of the atom, the settings provide access to the selected subshells and energies; 
        a list of crossSections::Array{ImpactIonization.CrossSection,1} is returned.
"""    
function  computeCrossSections(basis::Basis, grid::Radial.Grid, nm::Nuclear.Model, settings::ImpactIonization.Settings; output::Bool=true)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    #
    println("")
    printstyled("ImpactIonization.computeCrossSections(): Use model $(settings.model) ... \n", color=:light_green)
    printstyled("---------------------------------------------------------------------------------------- \n", color=:light_green)
    println("")
    #
    newCrossSections = ImpactIonization.CrossSection[]
    crossSections    = ImpactIonization.determineCrossSections(basis, settings)
    # 
    for  cs in crossSections
        newCs = ImpactIonization.computeCrossSections(settings.model, cs, basis, nm, grid, settings) 
        push!( newCrossSections, newCs)
    end
    # Print all results to screen
    ImpactIonization.displayCrossSections(stdout, newCrossSections, settings)
    if  printSummary   ImpactIonization.displayCrossSections(iostream, newCrossSections, settings)     end
    #
    if    output    return( newCrossSections )
    else            return( nothing )
    end
end

    
"""
`ImpactIonization.computeCrossSections(model::BEBmodel, cs::ImpactIonization.CrossSection,
                                        basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, settings::ImpactIonization.Settings)`  
    ... to compute the particular EII cross sections (cs) for an atom or ion in BEBmodel. While the basis describes the
        configuration of the atom, the settings provide access to the selected subshells and energies; 
        a cs::ImpactIonization.CrossSection is returned.
"""    
function  computeCrossSections(model::BEBmodel, cs::ImpactIonization.CrossSection,
                                basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, settings::ImpactIonization.Settings)
    bindingEnergy = kineticEnergy = q = t = u = 0.0
    partialCS     = 0.0 
    Ryd           = 13.6
    # Determine the binding energy, kinetic energy and occupation number for the subshell
    bindingEnergy    = - basis.orbitals[cs.subshell].energy
    orb              = basis.orbitals[cs.subshell]
    kineticEnergy    = RadialIntegrals.isotope_nms(orb, orb, 0.0, grid)
    occupationNumber = Basics.computeMeanSubshellOccupation(cs.subshell, basis)
    eC               = ImpactIonization.effectiveCharge(cs.subshell, nm, basis)
    # Transform the imput energies into atomic unit and get the default values
    iE = Defaults.convertUnits("energy: from eV to atomic", cs.impactEnergy)
    ryd = Defaults.convertUnits("energy: from eV to atomic", Ryd)  # epislon/R
    # Define parameters using symbols in BEB & BED formulas
    t = iE / bindingEnergy     
    u = kineticEnergy / bindingEnergy  
    q = occupationNumber 
    @show cs.subshell, eC
    # Calculate the partial and total cross sections
    csFactor = 4 * pi * q * ryd^2 * (1/bindingEnergy)^2/(t + (u + 1)/(eC+1))
    wc       = 1 - 1/t - log(t)/(1+t)
    wd       = log(t) * 1/2 * (1 - 1/t^2 )            
    if t < 1   else    partialCS = csFactor * (wc + wd)     end 
    newCs = ImpactIonization.CrossSection(cs.subshell, cs.impactEnergy, partialCS, Float64[], Float64[])
    return( newCs )
end

    
"""
`ImpactIonization.computeCrossSections(model::BEDmodel, cs::ImpactIonization.CrossSection,
                                        basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, settings::ImpactIonization.Settings)`  
    ... to compute the EII cross sections (cs) for an atom or ion in BEDmodel. While the basis describes the
        configuration of the atom, the settings provide access to the selected subshells and energies; 
        a cs::ImpactIonization.CrossSection is returned.
"""    
function  computeCrossSections(model::BEDmodel, cs::ImpactIonization.CrossSection,
                                basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, settings::ImpactIonization.Settings)
    bindingEnergy = kineticEnergy = q = t = u = 0.0
    partialCS     = 0.0 
    # Determine the binding energy, kinetic energy and occupation number for the subshell
    bindingEnergy    = - basis.orbitals[cs.subshell].energy
    orb              = basis.orbitals[cs.subshell]
    kineticEnergy    = RadialIntegrals.isotope_nms(orb, orb, 0.0, grid)
    occupationNumber = Basics.computeMeanSubshellOccupation(cs.subshell, basis)
    eC               = ImpactIonization.effectiveCharge(cs.subshell, nm, basis)
    # Transform the imput energies into atomic unit and get the default values
    iE = Defaults.convertUnits("energy: from eV to atomic", cs.impactEnergy); @show iE
    # Define parameters using symbols in BEB & BED formulas
    t = iE / bindingEnergy; @show t     
    u = kineticEnergy / bindingEnergy  
    q = occupationNumber 
    # Calculate the dipole term
    k = 2*bindingEnergy
    f = x -> (sqrt(2*x[2])*(2*x[2]+k))/(x[1]*((x[1]+sqrt(2*x[2]))^2+k)^3*(((x[1])-sqrt(2*x[2]))^2+k)^3) 
    xmin,xmax = (bindingEnergy/sqrt(2*iE), 0), (sqrt(8*iE)*(1+t), iE-bindingEnergy) 
    (val,err) = hcubature(f, xmin, xmax; reltol=1e-8, abstol=0, maxevals=0) 
    # Calculate the partial and total cross sections
    Gfactor = val
    Hfactor = 0.5 * (1/bindingEnergy - 1/iE - log(t)/(iE + bindingEnergy))
    Sfactor = 4 * pi * q /(iE +(bindingEnergy + kineticEnergy)/(eC + 1))
    Ffactor = 32 * q * (sqrt(2*bindingEnergy))^3 /iE            
    if t < 1   else    partialCS = Sfactor *Hfactor + Ffactor *Gfactor     end   
    newCs = ImpactIonization.CrossSection(cs.subshell, cs.impactEnergy, partialCS, Float64[], Float64[])
    return( newCs )
end

    
"""
`ImpactIonization.computeCrossSections(model::FittedBEDmodel, cs::ImpactIonization.CrossSection,
                                        basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, settings::ImpactIonization.Settings)`  
    ... to compute the EII cross sections (cs) for an atom or ion in FittedBEDmodel. While the basis describes the
        configuration of the atom, the settings provide access to the selected subshells and energies; 
        a cs::ImpactIonization.CrossSection is returned.
"""    
function  computeCrossSections(model::FittedBEDmodel, cs::ImpactIonization.CrossSection,
                                basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, settings::ImpactIonization.Settings)
    subshell    = cs.subshell
    # Determine a scaling factor in the FittedBEDmodel
    function scalingFactor(subshell::Subshell, nm::Nuclear.Model)
        parameters = [ (Subshell("1s_1/2"),  0.,    0.), 
                       (Subshell("2s_1/2"), -0.02,  0.001),  (Subshell("2p_1/2"),  0.13, 0.11),  (Subshell("2p_3/2"), 0.04, 0.32), 
                       (Subshell("3s_1/2"), -0.01,  0.01),   (Subshell("3p_1/2"), -0.07, 0.01),  (Subshell("3p_3/2"), 0.05, 0.04), 
                       (Subshell("3d_3/2"),  0.45,  0.05),   (Subshell("3d_5/2"),  0.40, 0.01)  ]
        hasfound = false
        m = 0.0 
        lambda = 0.0
        for  par in parameters   
            if par[1] == subshell  m = par[2];  lambda = par[3];   hasfound = true;    break    end
        end
        if  !hasfound   error("Can only calculate K- L- and M-shells!")     end
        return  1 + m * nm.Z^lambda
    end
        scale_factor_value = scalingFactor(subshell, nm)
       
    bindingEnergy = kineticEnergy = q = t = u = 0.0
    partialCS     = 0.0 
    # Determine the binding energy, kinetic energy and occupation number for the subshell
    bindingEnergy    = - basis.orbitals[cs.subshell].energy
    orb              = basis.orbitals[cs.subshell]
    kineticEnergy    = RadialIntegrals.isotope_nms(orb, orb, 0.0, grid)
    occupationNumber = Basics.computeMeanSubshellOccupation(cs.subshell, basis)
    eC               = ImpactIonization.effectiveCharge(cs.subshell, nm, basis)
    # Transform the imput energies into atomic unit and get the default values
    iE    = Defaults.convertUnits("energy: from eV to atomic", cs.impactEnergy)
    mC2   = Defaults.getDefaults("mc^2") 
    # Define parameters using symbols in BEB & BED formulas
    t = iE / bindingEnergy     
    u = kineticEnergy / bindingEnergy  
    q = occupationNumber 
    # Define parameters using symbols in relativisticBEB & BED formulas
    tP      = iE / mC2;          uP      = kineticEnergy / mC2; bP      = bindingEnergy / mC2
    relTpra = 1 - 1/(1 + tP)^2;  relUpra = 1 - 1/(1 + uP)^2;    relBpra = 1 - 1/(1 + bP)^2
    # Calculate the dipole term
    k = 2*bindingEnergy
    f = x -> (sqrt(2*x[2])*(2*x[2]+k))/(x[1]*((x[1]+sqrt(2*x[2]))^2+k)^3*(((x[1])-sqrt(2*x[2]))^2+k)^3) 
    xmin,xmax = (bindingEnergy/sqrt(2*iE), 0), (sqrt(8*iE)*(1+tP), iE-bindingEnergy) 
    (val,err) = hcubature(f, xmin, xmax; reltol=1e-8, abstol=0, maxevals=0) 
    Gfactor = val
    # Calculate the partial and total cross sections 
    Hfactor   = 0.5 * (1/bindingEnergy - 1/iE - log(t)/(iE + bindingEnergy))
    Sfactor   = 4 * pi * q /(iE +(bindingEnergy + kineticEnergy)/(eC + 1))
    Ffactor   = 32 * q * (sqrt(2*bindingEnergy))^3 /iE 
    RfactorB  = (((1 + t)*(t + 2/bP)*(1 + 1/bP)^2)/((1 + 2/bP)/(bP^2) + t*(t + 2/bP)*(1 + 1/bP)^2))^1.5
    RfactorLM = (1 + 2/bP)/(0.6 * t +2/bP) * ((t + 1/bP)/(1 + 1/bP))^2.11
    if t < 1   else    partialCS = scale_factor_value * (RfactorLM *Sfactor *Hfactor + RfactorB *Ffactor *Gfactor)    end  
    newCs = ImpactIonization.CrossSection(cs.subshell, cs.impactEnergy, partialCS, Float64[], Float64[])
    return( newCs )
end


"""
`ImpactIonization.computeCrossSections(model::RelativisticBEBmodel, cs::ImpactIonization.CrossSection,
                                        basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, settings::ImpactIonization.Settings)`  
    ... to compute the particular EII cross sections (cs) for an atom or ion in BEBmodel. While the basis describes the
        configuration of the atom, the settings provide access to the selected subshells and energies; 
        a cs::ImpactIonization.CrossSection is returned.
"""    
function  computeCrossSections(model::RelativisticBEBmodel, cs::ImpactIonization.CrossSection,
                                basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, settings::ImpactIonization.Settings)
    bindingEnergy = kineticEnergy = q = t = u = 0.0
    partialCS     = 0.0 
    # Determine the binding energy, kinetic energy and occupation number for the subshell
    bindingEnergy    = - basis.orbitals[cs.subshell].energy
    orb              = basis.orbitals[cs.subshell]
    kineticEnergy    = RadialIntegrals.isotope_nms(orb, orb, 0.0, grid)
    occupationNumber = Basics.computeMeanSubshellOccupation(cs.subshell, basis)
    eC               = ImpactIonization.effectiveCharge(cs.subshell, nm, basis)
    # Transform the imput energies into atomic unit and get the default values
    iE    = Defaults.convertUnits("energy: from eV to atomic", cs.impactEnergy)
    mC2   = Defaults.getDefaults("mc^2") 
    alpha = Defaults.getDefaults("alpha")
    # Define parameters using symbols in BEB & BED formulas
    t = iE / bindingEnergy     
    u = kineticEnergy / bindingEnergy  
    q = occupationNumber 
    # Define parameters using symbols in relativisticBEB & BED formulas
    tP      = iE / mC2;          uP      = kineticEnergy / mC2; bP      = bindingEnergy / mC2
    relTpra = 1 - 1/(1 + tP)^2;  relUpra = 1 - 1/(1 + uP)^2;    relBpra = 1 - 1/(1 + bP)^2
    # Calculate the partial and total cross sections
    csFactor       = 4 * pi * alpha^4 * q /( 2 * bP) * 1 /(relTpra + (relBpra + relUpra)/(eC +1))
    wa             = log( relTpra / (1 - relTpra) ) - relTpra - log(2 * bP) 
    wb             = 1 - 1/t - log(t)/(1 + t) * (1 + 2 * tP)/(1 + tP / 2)^2 + bP^2/(1 + tP / 2)^2 *(t - 1)/2        
    analyticDipole = 0.5 * (1 - 1 / t^2)            
    if t < 1   else    partialCS = csFactor * (wa * analyticDipole +  wb)     end 
    newCs = ImpactIonization.CrossSection(cs.subshell, cs.impactEnergy, partialCS, Float64[], Float64[])
    return( newCs )
end    


"""
`ImpactIonization.computeCrossSections(model::RelativisticBEDmodel, cs::ImpactIonization.CrossSection,
                                        basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, settings::ImpactIonization.Settings)`  
    ... to compute the particular EII cross sections (cs) for an atom or ion in BEBmodel. While the basis describes the
        configuration of the atom, the settings provide access to the selected subshells and energies; 
        a cs::ImpactIonization.CrossSection is returned.
"""    
function  computeCrossSections(model::RelativisticBEDmodel, cs::ImpactIonization.CrossSection,
                                basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, settings::ImpactIonization.Settings)
    bindingEnergy = kineticEnergy = q = t = u = 0.0
    partialCS     = 0.0 
    # Determine the binding energy, kinetic energy and occupation number for the subshell
    bindingEnergy    = - basis.orbitals[cs.subshell].energy
    orb              = basis.orbitals[cs.subshell]
    kineticEnergy    = RadialIntegrals.isotope_nms(orb, orb, 0.0, grid)
    occupationNumber = Basics.computeMeanSubshellOccupation(cs.subshell, basis)
    eC               = ImpactIonization.effectiveCharge(cs.subshell, nm, basis)
    # Transform the imput energies into atomic unit and get the default values
    iE    = Defaults.convertUnits("energy: from eV to atomic", cs.impactEnergy)
    mC2   = Defaults.getDefaults("mc^2") 
    alpha = Defaults.getDefaults("alpha")
    # Define parameters using symbols in BEB & BED formulas
    t = iE / bindingEnergy     
    u = kineticEnergy / bindingEnergy  
    q = occupationNumber 
    # Define parameters using symbols in relativisticBEB & BED formulas
    tP      = iE / mC2;          uP      = kineticEnergy / mC2; bP      = bindingEnergy / mC2
    relTpra = 1 - 1/(1 + tP)^2;  relUpra = 1 - 1/(1 + uP)^2;    relBpra = 1 - 1/(1 + bP)^2
    # Calculate the dipole term
    k = 2*bindingEnergy
    f = x -> (sqrt(2*x[2])*(2*x[2]+k))/(x[1]*((x[1]+sqrt(2*x[2]))^2+k)^3*(((x[1])-sqrt(2*x[2]))^2+k)^3) 
    xmin,xmax = (bindingEnergy/sqrt(2*iE), 0), (sqrt(8*iE)*(1+tP), iE-bindingEnergy) 
    (val,err) = hcubature(f, xmin, xmax; reltol=1e-8, abstol=0, maxevals=0) 
    # Calculate the partial and total cross sections
    Gfactor = val
    Hfactor = 0.5 * (1/bindingEnergy - 1/iE - log(t)/(iE + bindingEnergy))
    Sfactor = 4 * pi * alpha^2 * q / (relTpra + (relBpra + relUpra)/(eC +1))
    Ffactor = 64 * q * (sqrt(relUpra))^3 /relTpra/alpha            
    if t < 1   else    partialCS = Sfactor *Hfactor + Ffactor *Gfactor     end 
    newCs = ImpactIonization.CrossSection(cs.subshell, cs.impactEnergy, partialCS, Float64[], Float64[])
    return( newCs )
end  


"""
`ImpactIonization.computeMultipleCrossSections(model::DirectMultipleModel, cs::ImpactIonization.CrossSection,
                                               basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, settings::ImpactIonization.Settings)`  
    ... to compute the particular EIMI cross sections (cs) for an atom or ion in DirectMultipleModel. While the basis describes the
        configuration of the atom, the settings provide access to the selected subshells and energies; 
        a cs::ImpactIonization.MultipleCrossSection is returned.
"""    
function  computeCrossSections(model::DirectMultipleModel, cs::ImpactIonization.CrossSection,
                               basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, settings::ImpactIonization.Settings)
    # The parameters fitted by Belenger et al, a, b for multiply charged N, c=1 for neutral while 0.75 for other ions.
    param_dir_mul_a = [14.0, 6.30, 0.50, 0.14, 0.049, 0.021, 0.0096, 0.0049, 0.0027]
    param_dir_mul_b = [1.08, 1.20, 1.73, 1.85, 1.96, 2.00, 2.00, 2.00, 2.00]
    param_dir_mul_c = [1, 0.75]                            
    totalEnergy     = u = 0.0
    partialCS       = 0.0
    Ryd             = 13.6

    atomicNumber = Int.(nm.Z); 
    multipleN    = settings.multipleN; 

    if  multipleN  >  10
        a1 = 1350 * multipleN^(-5.7)
        b1 = 2
    else
        a1 = param_dir_mul_a[multipleN - 1]
        b1 = param_dir_mul_b[multipleN - 1]
    end
    
    # Here we calculate total electron number and charge state.
    conf          = Basics.extractConfigurations(Basics.FromBasis(), basis)[1];   
    totalElectron = conf.NoElectrons; @show totalElectron
    chargeState   = atomicNumber - totalElectron; @show chargeState
   
    if chargeState == 0
        c1 = param_dir_mul_c[1]
    else
        c1 = param_dir_mul_c[2]
    end

    if totalElectron < multipleN
        printstyled("inputError: Target electron number must be more than ionization order!\n", color=:light_green)
        error("inputError: InDirectDoubleModel only works for double ionization!")
    end


    # This part is left for calculating the ionization potential from JAC, to be added
    if false
        valenceShell = Basics.extractValenceShell(basis);                    @show valenceShell
        
        confs        = Basics.extractConfigurations(Basics.FromBasis(), basis);   @show confs
        #
        # Now generate N-1 (N = multipleN from your settings) new bases to extract improved one-particle energies
        newBases = ManyElectron.Basis[];  @show newBases
        totalEnergy = Float64[];
        for  N = 1:multipleN
            newConfs      = Basics.generateConfigurations(Basics.RemoveElectrons([valenceShell]), confs)
            @show confs
            multipletX    = SelfConsistent.performSCF(computation.configs, nm, computation.grid, asfSettings)
            newBasis      = multipletX.levels[1].basis
            push!(newBases, newBasis)
            #
            confs        = Basics.extractConfigurations(Basics.FromBasis(), newBasis);   @show confs
            valenceShell = Basics.extractValenceShell(newBasis);                    @show valenceShell
        end    
        for  i = 1:multipleN
            BIP = basis[i].orbitals[subshell]
            push!(totalEnergy,BIP)
        end
        totalEnergyNew = sum(totalEnergy[chargeState+1:chargeState+multipleN]);     @show totalEnergyNew
    end

    EIP = PeriodicTable.ionizationPotentials_Nist2025(atomicNumber);    @show EIP

    totalEnergy = sum(EIP[chargeState+1 : chargeState+multipleN]);      @show totalEnergy
    ## error("xx")

    # Impact energy and parameter u
    iE  = cs.impactEnergy
    u   = iE / totalEnergy;                                             @show u     
    
    # Calculate the total cross section   
    term11 = (a1 * totalElectron^b1) / (totalEnergy / Ryd)^2
    term12 = ((u + 1) / u)^c1
    term13 = log(u) / u 
      
    Barn2Au = 0.280_028_520_292_481_56e8  
    if  u < 1   else   partialCS = term11 * term12 * term13 / Barn2Au     end 
    newCs = ImpactIonization.CrossSection(cs.subshell, cs.impactEnergy, partialCS, Float64[], Float64[])
    
    return( newCs )
end

    
"""
`ImpactIonization.computeMultipleCrossSections(model::InDirectDoubleModel, cs::ImpactIonization.CrossSection,
                                               basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, settings::ImpactIonization.Settings)`  
    ... to compute the particular EIMI cross sections (cs) for an atom or ion in InDirectDoubleModel. While the basis describes the
        configuration of the atom, the settings provide access to the selected subshells and energies; 
        a cs::ImpactIonization.MultipleCrossSection is returned.
"""    
function  computeCrossSections(model::InDirectDoubleModel, cs::ImpactIonization.CrossSection,
                                basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, settings::ImpactIonization.Settings)
    # The parameters fitted by shevelko et al.: 
    # [atomic number; charge state; order; A_dir; [C]; [PHI]; [electron number of inner shell]; [threshold of inner shell]; 
    # [branching     ratio of autoionization]]       
    param_shevelko_fitting = [ [2, 0, 2, 5.0, [0], [0], [2], [24.587], [1.0]],
                               [3, 0, 2, 12.0, [0], [0], [2], [58], [1.0]],
                               [3, 1, 2, 5.0, [0], [0], [2], [75.638], [1.0]],
                               [5, 0, 2, 10.0, [3.6], [5.0], [2], [200.33], [1.0]],
                               [5, 1, 2, 1.8, [3.6], [4.5], [2], [217.66], [1 - 7129 / 1e4]],
                               [5, 3, 2, 5.0, [0], [0], [2], [259.368], [1.0]]                 ]

    totalEnergy = u = EIP = 0.0 
    Ryd         = 13.6
    partialCS   = 0.0
    
    atomicNumber = Int.(nm.Z); 
    multipleN = settings.multipleN; 

    # InDirectDoubleModel only works for double ionzation
    if multipleN != 2
        printstyled("inputError: InDirectDoubleModel only works for double ionization!\n", color=:light_green)
        error("inputError: InDirectDoubleModel only works for double ionization!")
    end
    
    all_atom_shev   = []
    all_charge_shev = []

    # Get the Z and charge_state for Shevelko
    for js2 in 1:length(param_shevelko_fitting)
        push!(all_atom_shev, param_shevelko_fitting[js2][1])
        push!(all_charge_shev, param_shevelko_fitting[js2][2])
    end
    
    # Here we calculate the total electron number and charge state.
    conf          = Basics.extractConfigurations(Basics.FromBasis(), basis)[1];   
    totalElectron = conf.NoElectrons; @show totalElectron
    chargeState   = atomicNumber - totalElectron; @show chargeState
   
    if totalElectron < multipleN
        printstyled("inputError: Target electron number must be more than ionization order!\n", color=:light_green)
        error("inputError: Target electron number must be more than ionization order!")
    end

    # Calculate total ionization potential
    EIP         = PeriodicTable.ionizationPotentials_Nist2025(atomicNumber); @show EIP
    totalEnergy = sum(EIP[chargeState+1 : chargeState+multipleN]); @show totalEnergy
     
    A_dir = C_GAMMA = PHI_GAMMA =0.0
    inner_electron_number = 1
    inner_threshold = f_BR = 0.0
     
    if !(atomicNumber in all_atom_shev) || !(chargeState in all_charge_shev) 
        printstyled("InputError: Input atomic number or charge state out of range!\n Parameters are available for limited ions!" *
                    "\n To be added in near future. \n", color=:light_green)
        error("InputError: Input atomic number or charge state out of range!")
    else
        for sf in 1:length(param_shevelko_fitting)
            if param_shevelko_fitting[sf][1] == atomicNumber && param_shevelko_fitting[sf][2] == chargeState && 
               param_shevelko_fitting[sf][3] == multipleN
                A_dir                 = param_shevelko_fitting[sf][4]
                C_GAMMA               = param_shevelko_fitting[sf][5]
                PHI_GAMMA             = param_shevelko_fitting[sf][6]
                inner_electron_number = param_shevelko_fitting[sf][7]
                inner_threshold       = param_shevelko_fitting[sf][8]
                f_BR                  = param_shevelko_fitting[sf][9]                
            end
        end
    end

    iE          = cs.impactEnergy
    u           = iE / totalEnergy

    # Initialize the direct and indirect cross section
    sigma_dir   = 0.0
    sigma_indir = 0.0

    # Calculate the direct ionization part
    term21    = 1 - exp(-3 * (u - 1))
    term22    = A_dir * (u - 1) / totalEnergy^3 / (u + 0.5)^2
    sigma_dir = term21 * term22
    # Calculate the indirect ionization part        
    for epsilon in iE
        for i2 in 1:length(inner_threshold)
            if epsilon < inner_threshold[i2]
                IA = 0
            else
                alpha_gamma = f_BR[i2]
                c_gamma     = C_GAMMA[i2]
                I_gamma     = inner_threshold[i2]
                phi_gamma   = PHI_GAMMA[i2]
                x2          = epsilon / I_gamma
                IA          = alpha_gamma * c_gamma * (x2 - 1) / I_gamma^2 / x2 / (x2 + phi_gamma)
            end
            sigma_indir += IA
        end
        end
     
    # Calculate the total cross section
    Barn2Au = 0.280_028_520_292_481_56e8 
    if  u < 1   else    partialCS = (sigma_dir + sigma_indir) * 1e5 / Barn2Au     end 
    newCs = ImpactIonization.CrossSection(cs.subshell, cs.impactEnergy, partialCS, Float64[], Float64[])
    
    return( newCs )
end


"""
`ImpactIonization.computeMultipleCrossSections(model::DoubleExperimentModel, cs::ImpactIonization.CrossSection,
                                               basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, settings::ImpactIonization.Settings)`  
    ... to compute the particular EIMI cross sections (cs) for an atom or ion in DoubleExperimentModel. While the basis describes the
        configuration of the atom, the settings provide access to the selected subshells and energies; 
        a cs::ImpactIonization.MultipleCrossSection is returned.
"""    
function  computeCrossSections(model::DoubleExperimentModel, cs::ImpactIonization.CrossSection,
                               basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, settings::ImpactIonization.Settings)
    # The parameters fitted by Shevelko et al.: 
    # [atomic number; charge state; order; A_dir; [C]; [PHI]; [electron number of inner shell]; [threshold of inner shell]; 
    # [branching ratio of autoionization]]
    param_shevelko_fitting = [  [2, 0, 2, 5.0, [0], [0], [2], [24.587], [1.0]],
                                [3, 0, 2, 12.0, [0], [0], [2], [58], [1.0]],
                                [3, 1, 2, 5.0, [0], [0], [2], [75.638], [1.0]],
                                [5, 0, 2, 10.0, [3.6], [5.0], [2], [200.33], [1.0]],
                                [5, 1, 2, 1.8, [3.6], [4.5], [2], [217.66], [1 - 7129 / 1e4]],
                                [5, 3, 2, 5.0, [0], [0], [2], [259.368], [1.0]]              ]                           
    totalEnergy = u = 0.0
    partialCS   = 0.0
    Ryd         = 13.6
    
    atomicNumber = Int.(nm.Z); 
    multipleN    = settings.multipleN; 
    
    # InDirectDoubleModel only works for double ionzation
    if multipleN != 2
        printstyled("inputError: DoubleExperimentModel only works for double ionization!\n", color=:light_green)
        error("inputError: DoubleExperimentModel only works for double ionization!")
    end
    
    
    all_atom_shev = []
    all_charge_shev = []
    # get the Z and charge_state for shevelko
    for js2 in 1:length(param_shevelko_fitting)
        push!(all_atom_shev, param_shevelko_fitting[js2][1])
        push!(all_charge_shev, param_shevelko_fitting[js2][2])
    end
    
    # Here we calculate the total electron number and charge state.
    conf          = Basics.extractConfigurations(Basics.FromBasis(), basis)[1];   
    totalElectron = conf.NoElectrons;               @show totalElectron
    chargeState   = atomicNumber - totalElectron;   @show chargeState
   
    if totalElectron < multipleN
        printstyled("inputError: Target electron number must be more than ionization order!\n", color=:light_green)
        error("inputError: Target electron number must be more than ionization order!")
    end

   
    # Calculate the total ionization potential
    EIP         = PeriodicTable.ionizationPotentials_Nist2025(atomicNumber);   @show EIP
    totalEnergy = sum(EIP[chargeState+1 : chargeState+multipleN]);             @show totalEnergy

    # parameters used in the formula
    A_dir = C_GAMMA = PHI_GAMMA = 0.0
    inner_electron_number  = 1
    inner_threshold = f_BR = 0.0
     
    if !(atomicNumber in all_atom_shev)  ||  !(chargeState in all_charge_shev) 
        printstyled("InputError: Input atomic number or charge state out of range!\n Parameters are available for limited ions! " *
                    "\n To be added in near future. \n", color=:light_green)
        error("InputError: Input atomic number or charge state out of range!")
    else
        for sf in 1:length(param_shevelko_fitting)
            if param_shevelko_fitting[sf][1] == atomicNumber && param_shevelko_fitting[sf][2] == chargeState && param_shevelko_fitting[sf][3] == multipleN
                A_dir                 = param_shevelko_fitting[sf][4]
                C_GAMMA               = param_shevelko_fitting[sf][5]
                PHI_GAMMA             = param_shevelko_fitting[sf][6]
                inner_electron_number = param_shevelko_fitting[sf][7]
                inner_threshold       = param_shevelko_fitting[sf][8]
                f_BR                  = param_shevelko_fitting[sf][9]                
            end
        end
    end
        
    # Impact energy
    iE = cs.impactEnergy
    u  = iE / totalEnergy

    # Initialize the direct and indirect cross section
    sigma_dir   = 0.0
    sigma_indir = 0.0

    # Calculate direct ionization part 
    term31    = 1 - exp(-3 * (u - 1))
    term32    = A_dir * (u - 1) / totalEnergy^3 / (u + 0.5)^2
    term33    = 1 + 0.1 * log(4 * u + 1)
    sigma_dir = term31 * term32 * term33
            
    # Calculate indirect ionization part
    for epsilon in iE
        for i3 in 1:length(inner_threshold)
            if epsilon < inner_threshold[i3]
                IA          = 0
            else
                t           = epsilon / totalEnergy
                alpha_gamma = f_BR[i3]
                c_gamma     = C_GAMMA[i3]
                I_gamma     = inner_threshold[i3]
                phi_gamma   = PHI_GAMMA[i3]
                x3          = epsilon / I_gamma 
                IA1         = alpha_gamma * c_gamma * (x3 - 1) / I_gamma^2 / x3 / (x3 + 5)
                IA2         = 1 + 0.3 / phi_gamma / log(4 * t + 1)
                IA          = IA1 * IA2
            end
            sigma_indir += IA
        end
    end
     
    # Calculate the total cross section     
    Barn2Au = 0.280_028_520_292_481_56e8 
    
    if  u < 1   else    partialCS = (sigma_dir + sigma_indir) * 1e5 / Barn2Au    end 
    newCs = ImpactIonization.CrossSection(cs.subshell, cs.impactEnergy, partialCS, Float64[], Float64[])
    
    return( newCs )
end


"""
`ImpactIonization.computeMultipleCrossSections(model::LotzMultipleModel, cs::ImpactIonization.CrossSection,
                                               basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, settings::ImpactIonization.Settings)`  
    ... to compute the particular EIMI cross sections (cs) for an atom or ion in LotzMultipleModel. While the basis describes the
        configuration of the atom, the settings provide access to the selected subshells and energies; 
        a cs::ImpactIonization.MultipleCrossSection is returned.
"""    
function  computeCrossSections(model::LotzMultipleModel, cs::ImpactIonization.CrossSection,
                                basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, settings::ImpactIonization.Settings)
                        
    param_shevelko_fitting = [  [2, 0, 2, 5.0, [0], [0], [2], [24.587], [1.0]],
                                [3, 0, 2, 12.0, [0], [0], [2], [58], [1.0]],
                                [3, 1, 2, 5.0, [0], [0], [2], [75.638], [1.0]],
                                [5, 0, 2, 10.0, [3.6], [5.0], [2], [200.33], [1.0]],
                                [5, 1, 2, 1.8, [3.6], [4.5], [2], [217.66], [1 - 7129 / 1e4]],
                                [5, 3, 2, 5.0, [0], [0], [2], [259.368], [1.0]]     ]
                         
    totalEnergy  = u = 0.0
    partialCS    = 0.0
    
    atomicNumber = Int.(nm.Z); 
    multipleN    = settings.multipleN; 
 
    if multipleN != 2
        printstyled("inputError: LotzMultipleModel currently only works for double ionization!" *
                    "\n Branching ratio parameters are available for limited ions! \n To be added in near future. \n", color=:light_green)
        error("inputError: LotzMultipleModel currentl only works for double ionization!")
    end
    
    all_atom_shev   = []
    all_charge_shev = []
    # Get the Z and charge_state for shevelko
    for js2 in 1:length(param_shevelko_fitting)
        push!(all_atom_shev, param_shevelko_fitting[js2][1])
        push!(all_charge_shev, param_shevelko_fitting[js2][2])
    end
    
    # Calculate the total electron numbers and charge state
    conf          = Basics.extractConfigurations(Basics.FromBasis(), basis)[1];   
    totalElectron = conf.NoElectrons; @show totalElectron
    chargeState   = atomicNumber - totalElectron; @show chargeState
   
    if totalElectron < multipleN
       println("inputError: Target electron number must be more than ionization order!")
    end
   
    # Calculate the total ionization potential
    EIP         = PeriodicTable.ionizationPotentials_Nist2025(atomicNumber);     @show EIP
    totalEnergy = sum(EIP[chargeState+1 : chargeState+multipleN]);               @show totalEnergy

    # Initialize parameter
    A_dir = C_GAMMA = PHI_GAMMA =0.0
    inner_electron_number  = 1
    inner_threshold = f_BR = 0.0
     
    if !(atomicNumber in all_atom_shev)  ||  !(chargeState in all_charge_shev) 
        println("InputError: Input atomic number or charge state out of range!")
        exit(1)
    else
        for sf in 1:length(param_shevelko_fitting)
            if param_shevelko_fitting[sf][1] == atomicNumber && param_shevelko_fitting[sf][2] == chargeState && 
               param_shevelko_fitting[sf][3] == multipleN
                A_dir                 = param_shevelko_fitting[sf][4]
                C_GAMMA               = param_shevelko_fitting[sf][5]
                PHI_GAMMA             = param_shevelko_fitting[sf][6]
                inner_electron_number = param_shevelko_fitting[sf][7]
                inner_threshold       = param_shevelko_fitting[sf][8]
                f_BR                  = param_shevelko_fitting[sf][9]                
            end
        end
    end
     
    # Impact energy
    iE = cs.impactEnergy
    u  = iE / totalEnergy
    # Define parameters using symbols in BEB & BED formulas
    sigma_indir = 0.0
    for epsilon in iE
        t = epsilon / totalEnergy
        for i4 in 1:length(inner_threshold)
            if epsilon < inner_threshold[i4]
                IA = 0
            else
                alpha_gamma = f_BR[i4]
                Ns          = inner_electron_number[i4]
                I_gamma     = inner_threshold[i4]
                IA = 4.5 * alpha_gamma * Ns * log(t) / I_gamma^2 / t
            end
            sigma_indir += IA
        end
    end
    
    Barn2Au = 0.280_028_520_292_481_56e8    
    if u < 1   else    partialCS = sigma_indir * 1e4 / Barn2Au    end 
    newCs = ImpactIonization.CrossSection(cs.subshell, cs.impactEnergy, partialCS, Float64[], Float64[])
    
    return( newCs )
end
    

"""
`ImpactIonization.determineCrossSections(basis::Basis, settings::ImpactIonization.Settings)`  
    ... to determine an array of EII cross sections for all requested shells of the given multiplett that need to be calculated in these
        computations; an Array{ImpactIonization.CrossSection,1} is returned.
"""
function  determineCrossSections(basis::Basis, settings::ImpactIonization.Settings)
    crossSections = ImpactIonization.CrossSection[]
    if  settings.shellSelection.active
        subshells = Basics.generateSubshellList(settings.shellSelection.shells) 
    else
        confs     = Basics.extractConfigurations(Basics.FromBasis(), basis)
        shells    = Basics.extractShellList(confs)
        subshells = Basics.generateSubshellList(shells)
    end
    # Generate the list of cross sections
    for  subsh in subshells
        for  iE in  settings.impactEnergies
            push!(crossSections, ImpactIonization.CrossSection(subsh, iE, 0., Float64[], Float64[]))
        end
    end
    
    return( crossSections )
end


"""
`ImpactIonization.displayCrossSections(stream::IO, crossSections::Array{ImpactIonization.CrossSection,1}, settings::ImpactIonization.Settings)`  
    ... displays the EII cross sections in a neat tabular form; nothing is returned.
"""
function  displayCrossSections(stream::IO, crossSections::Array{ImpactIonization.CrossSection,1}, settings::ImpactIonization.Settings)
    # Display the partial (shell-dependent) cross sections
    if  settings.calcPartialCs 
        nx = 120 
        println(stream, " ")
        println(stream, "  Partial ionization cross sections for atoms by electron impact:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "subshell"         ; na=0);                       
        sb = sb * TableStrings.hBlank(18)
        sa = sa * TableStrings.center(18, "binding energy"   ; na=2);                       
        sb = sb * TableStrings.center(18, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.center(18, "kinetic energy"   ; na=4)               
        sb = sb * TableStrings.center(14, TableStrings.inUnits("energy"); na=0)
        sa = sa * TableStrings.center(12, "occupation number"  ; na=4)             
        sb = sb * TableStrings.hBlank(30) 
        sa = sa * TableStrings.center(12, "impact energies"  ; na=4)             
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)            
        sa = sa * TableStrings.center(28, "Cross section"; na=3)      
        sb = sb * TableStrings.center(30, TableStrings.inUnits("cross section"); na=3)
        println(stream, sa);    println(stream, sb);    println(stream, "  ", TableStrings.hLine(nx))
        #
        sa  = ""; 
        sa = sa * TableStrings.center(18, string(crossSections[1].subshell); na=6)
        println(stream, sa)            
        for  cs  in crossSections
            sb = TableStrings.hBlank(86)
            sb = sb * @sprintf("%.6e", cs.impactEnergy) * "             "    
            wa = cs.partialCS
            sb = sb * @sprintf("%.10e", Defaults.convertUnits("cross section: from atomic to barn", wa))
            println(sb)
        end
        println(stream, "  ", TableStrings.hLine(nx))    
    end
    #
    # Display the total EII cross sections
    if  settings.calcTotalCs 
        nx = 120 
        println(stream, " ")
        println(stream, "  Total ionization cross sections for atoms by electron impact:")
        println(stream, " ")
        println(stream, "  ", TableStrings.hLine(nx))
        sa = "  "
        for  iE in  settings.impactEnergies
        totalCs = 0.0
            for  cs  in crossSections
                if  cs.impactEnergy == iE   totalCs = totalCs + cs.partialCS      end
            end
            sa = TableStrings.hBlank(86)
            sa = sa * @sprintf("%.6e", iE) * "             "    
            wa = totalCs
            sa = sa * @sprintf("%.10e", Defaults.convertUnits("cross section: from atomic to barn", wa))
            println(stream, sa)            
        end            
        println(stream, "  ", TableStrings.hLine(nx))    
    end
    #
    return( nothing )
end


"""
`ImpactIonization.effectiveCharge(subshell::Subshell, nm::Nuclear.Model, basis::Basis)`  
    ... determined the effective charge that is felt by an electron in the given subshell;
        a zEff::Float64 is returned.
"""
function  effectiveCharge(subshell::Subshell, nm::Nuclear.Model, basis::Basis)
    wz        = nm.Z
    confs     = Basics.extractConfigurations(Basics.FromBasis(), basis)
    shells    = Basics.extractShellList(confs)
    subshells = Basics.generateSubshellList(shells)
    enSub     = basis.orbitals[subshell].energy
    for  subsh  in subshells   
        if  basis.orbitals[subsh].energy <= enSub
            wz = wz - Basics.computeMeanSubshellOccupation(subsh, basis) 
        end 
    end
    
    return(wz)
end

end # module
