
export compute


"""
`Basics.compute("angular coefficients: e-e, Ratip2013", csfa::CsfR, csfb::CsfR)`  
    ... to compute the angular coefficients in the decomposition of a (reduced) many-electron matrix element with a general 
        rank-0 electron-electron interaction operator ``⟨csf_a ||V(e-e)|| csfb⟩ = \\sum_t  T(a_t, b_t, c_t, d_t) * R^k (a_t, b_t, c_t, d_t)``
        by a call to the Fortran procedure `anco_calculate_csf_pair` of the RATIP program; a Tuple{Array{AngularTcoeff,1},Array{AngularVcoeff,1}}` 
        is returned.
"""
function Basics.compute(sa::String, csfa::CsfR, csfb::CsfR)
    if sa == "angular coefficients: e-e, Ratip2013"
        if csfa == csfb
            AngularCoefficientsRatip2013.load_csl(csfa)
            t_coeffs, v_coeffs = AngularCoefficientsRatip2013.angular_coefficients_pair(1,1)
        else
            AngularCoefficientsRatip2013.load_csl(csfa, csfb)
            t_coeffs, v_coeffs = AngularCoefficientsRatip2013.angular_coefficients_pair(1,2)
        end
    else 
        error("Unsupported keystring = ", sa)
    end
    return t_coeffs, v_coeffs
end


"""
`Basics.compute("angular coefficients: 1-p, Ratip2013", rank::Int64, csfa::CsfR, csfb::CsfR)`  
    ... to compute the the angular coefficients in the decomposition of a (reduced) many-electron matrix element with a general 
        single-particle operator of the given rank ``⟨csf_a ||O^rank|| csfb⟩ = \\sum_t  T(a_t, b_t) * R (a_t, b_t)``  by a call 
        to the Fortran procedure `anco_calculate_csf_pair_1p` of the RATIP program; an `Array{AngularTcoeff,1}` is returned.
"""
function Basics.compute(sa::String, rank, csfa::CsfR, csfb::CsfR)
    if      sa == "angular coefficients: 1-p, Ratip2013"
        if csfa == csfb
            AngularCoefficientsRatip2013.load_csl(csfa)
            t_coeffs = AngularCoefficientsRatip2013.angular_coefficients_pair_1p(rank,1,1)
        else
            AngularCoefficientsRatip2013.load_csl(csfa, csfb)
            t_coeffs = AngularCoefficientsRatip2013.angular_coefficients_pair_1p(rank,1,2)
        end
    else    error("Unsupported keystring = ", sa)
    end
    return t_coeffs
end


"""
`Basics.compute("angular coefficients: 1-p, Grasp92", parity, rank::Integer, csfa::CsfR, csfb::CsfR)`  
    ... to compute the the angular coefficients in the decomposition of a (reduced) many-electron matrix element with a general 
        single-particle operator of the given parity and rank by a call to the Fortran procedure `mct_generate_coefficients` 
        of the RATIP program; an `Array{AngularTcoeff,1}` is returned.
"""
function Basics.compute(sa::String, parity, rank::Integer, csfa::CsfR, csfb::CsfR)
    if      sa == "angular coefficients: 1-p, Grasp92"
        if csfa == csfb
            subshells  = AngularCoefficientsRatip2013.load_csl(csfa)
            mct_coeffs = AngularCoefficientsRatip2013.mct_generate_coefficients(1, 1, 1, 1, Int32(parity), rank)
        else
            # Add 'zeros' to the fields if the length of occupation does not agree
            if       ( nz = length(csfa.occupation) - length(csfb.occupation) ) <  0   csfa = Basics.addZerosToCsfR( -nz, csfa)
            elseif   ( nz = length(csfa.occupation) - length(csfb.occupation) ) >  0   csfb = Basics.addZerosToCsfR(  nz, csfb)    end
            subshells  = AngularCoefficientsRatip2013.load_csl(csfa, csfb)
            mct_coeffs = AngularCoefficientsRatip2013.mct_generate_coefficients(1, 1, 2, 2, Int32(parity), rank)
        end
    else    error("Unsupported keystring = ", sa)
    end
    return map(t -> AngularCoefficientsRatip2013.AngularTcoeff(t, subshells), mct_coeffs) # Convert Fmctcoefficient to AngularTcoeff
end



"""
`Basics.compute("matrix: CI, J^P symmetry", JP::LevelSymmetry, basis::Basis, nuclearModel::Nuclear.Model, grid::Radial.Grid,
                                            settings::AsfSettings; printout::Bool=true)`  
    ... to compute the CI matrix for a given J^P symmetry block of basis and by making use of the nuclear model and the grid; 
        a matrix::Array{Float64,2} is returned.
"""
function Basics.compute(sa::String, JP::LevelSymmetry, basis::Basis, nuclearModel::Nuclear.Model, grid::Radial.Grid,
                        settings::AsfSettings; printout::Bool=true)    
    !(sa == "matrix: CI, J^P symmetry")   &&   error("Not supported keystring")

    # Determine the dimension of the CI matrix and the indices of the CSF with J^P symmetry in the basis
    idx_csf = Int64[]
    for  idx = 1:length(basis.csfs)
        if  basis.csfs[idx].J == JP.J   &&   basis.csfs[idx].parity == JP.parity    push!(idx_csf, idx)    end
    end
    n = length(idx_csf)
    if printout    print("Compute CI matrix of dimension $n x $n for the symmetry $(string(JP.J))^$(string(JP.parity)) ...")    end

    # Generate an effective nuclear charge Z(r) on the given grid
    potential = Nuclear.nuclearPotential(nuclearModel, grid)
    if  settings.qedModel in [QedPetersburg(), QedSydney()]    
        meanPot = potential
        ## meanPot = compute("radial potential: Dirac-Fock-Slater", grid, basis)
        ## meanPot = Basics.add(potential, meanPot)   
    end   

    matrix = zeros(Float64, n, n)
    keep = true
    InteractionStrength.XL_Coulomb_reset_storage(keep, printout=false)
    InteractionStrength.XL_Breit_reset_storage(keep, printout=false)
    for  r = 1:n
        for  s = 1:n
            #
            if  settings.eeInteractionCI == DiagonalCoulomb()  &&  r != s    continue    end
            # Calculate the spin-angular coefficients
            if  Defaults.saRatip()
                waR = compute("angular coefficients: e-e, Ratip2013", basis.csfs[idx_csf[r]], basis.csfs[idx_csf[s]])
                wa  = waR       
            end
            if  Defaults.saGG()
                subshellList = basis.subshells
                opa  = SpinAngular.OneParticleOperator(0, plus, true)
                waG1 = SpinAngular.computeCoefficients(opa, basis.csfs[idx_csf[r]], basis.csfs[idx_csf[s]], subshellList) 
                opa  = SpinAngular.TwoParticleOperator(0, plus, true)
                waG2 = SpinAngular.computeCoefficients(opa, basis.csfs[idx_csf[r]], basis.csfs[idx_csf[s]], subshellList)
                wa   = [waG1, waG2]
            end
            if  Defaults.saRatip() && Defaults.saGG() && true
                if  length(waR[1]) != 0     println(  ">> Angular coeffients from Ratip2013   = $(waR[1]) ")    end
                if  length(waG1)   != 0     println("\n>> Angular coeffients from SpinAngular = $waG1 ")        end
                if  length(waR[2]) != 0     println(  ">> Angular coeffients from Ratip2013   = $(waR[2]) ")    end
                if  length(waG2)   != 0     println("\n>> Angular coeffients from SpinAngular = $waG2 ")        end
            end
            #
            me = 0.
            for  coeff in wa[1]
                jj = Basics.subshell_2j(basis.orbitals[coeff.a].subshell)
                me = me + coeff.T * sqrt( jj + 1) * RadialIntegrals.GrantIab(basis.orbitals[coeff.a], basis.orbitals[coeff.b], grid, potential)
                if  settings.qedModel != NoneQed()  
                    me = me + InteractionStrengthQED.qedLocal(basis.orbitals[coeff.a], basis.orbitals[coeff.b], nuclearModel, 
                                                                settings.qedModel, meanPot, grid)  
                end
            end

            for  coeff in wa[2]
                if  typeof(settings.eeInteractionCI) in [DiagonalCoulomb, CoulombInteraction, CoulombBreit, CoulombGaunt]
                    me = me + coeff.V * InteractionStrength.XL_Coulomb(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                                    basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid, keep=keep)
                elseif  false
                    xl1 = InteractionStrength.XL_Coulomb(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                                            basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid, keep=false)
                    xl2 = InteractionStrength.XL_Coulomb_WO(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                                                basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid)
                    if abs(xl1 - xl2) > 1.0e-12  println("XL_Coulomb differ: $xl1   $xl2")      end
                elseif  false 
                    me = me + coeff.V * InteractionStrength.XL_Coulomb_WO(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                                    basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid)   
                end
                                                                                        
                if      typeof(settings.eeInteractionCI) in [BreitInteraction, CoulombBreit, CoulombGaunt]
                    me = me + coeff.V * InteractionStrength.XL_Breit(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                                basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid,
                                                                                settings.eeInteractionCI, keep=keep)     
                elseif  false
                    xl1 = InteractionStrength.XL_Breit(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                                            basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid, keep=false)
                    xl2 = InteractionStrength.XL_Breit_WO(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                                                basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid)
                    if abs(xl1 - xl2) > 1.0e-12  println("XL_Breit differ: $xl1   $xl2")      end
                elseif  false 
                    me = me + coeff.V * InteractionStrength.XL_Breit_WO(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                                    basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid)
                end
            end
            matrix[r,s] = me
        end
    end 
    if printout    println("   ... done.")    end

    return( matrix )
end


"""
`Basics.compute(JP::LevelSymmetry, basis::Basis, nuclearModel::Nuclear.Model, grid::Radial.Grid,
                settings::AsfSettings, plasmaModel::Basics.AbstractPlasmaModel; printout::Bool=true)`  
    ... to compute the CI matrix for a given J^P symmetry block of basis and by making use of the nuclear model and the grid; 
        a matrix::Array{Float64,2} is returned.
"""
function Basics.compute(JP::LevelSymmetry, basis::Basis, nuclearModel::Nuclear.Model, grid::Radial.Grid,
                        settings::AsfSettings, plasmaModel::Basics.AbstractPlasmaModel; printout::Bool=true)    

    # Determine the dimension of the CI matrix and the indices of the CSF with J^P symmetry in the basis
    idx_csf = Int64[]
    for  idx = 1:length(basis.csfs)
        if  basis.csfs[idx].J == JP.J   &&   basis.csfs[idx].parity == JP.parity    push!(idx_csf, idx)    end
    end
    n = length(idx_csf)
    
    if  typeof(settings.eeInteractionCI) in [BreitInteraction, CoulombBreit]              
                                         error("No Breit interaction supported for plasma computations; use breitCI=false  in the asfSettings.")    end   
    if  settings.qedModel != NoneQed()   error("No QED estimates supported for plasma computations; use qedModel=NoneQed()  in the asfSettings.")   end   
    
    # Now distinguis the CI matrix for different plasma models
    if  typeof(plasmaModel) == Basics.DebyeHueckel
        if printout    print("Compute DebyeHueckel-CI matrix of dimension $n x $n for the symmetry $(string(JP.J))^$(string(JP.parity)) ...")    end

        # Generate an effective nuclear charge Z(r) for a screened Debye-Hueckel potential on the given grid
        potential = Nuclear.nuclearPotentialDH(nuclearModel, grid, 1/plasmaModel.debyeLength)

        matrix = zeros(Float64, n, n)
        for  r = 1:n
            for  s = 1:n
                # Calculate the spin-angular coefficients
                if  Defaults.saRatip()
                    waR = compute("angular coefficients: e-e, Ratip2013", basis.csfs[idx_csf[r]], basis.csfs[idx_csf[s]])
                    wa  = waR       
                end
                if  Defaults.saGG()
                    subshellList = basis.subshells
                    opa  = SpinAngular.OneParticleOperator(0, plus, true)
                    waG1 = SpinAngular.computeCoefficients(opa, basis.csfs[idx_csf[r]], basis.csfs[idx_csf[s]], subshellList) 
                    opa  = SpinAngular.TwoParticleOperator(0, plus, true)
                    waG2 = SpinAngular.computeCoefficients(opa, basis.csfs[idx_csf[r]], basis.csfs[idx_csf[s]], subshellList)
                    wa   = [waG1, waG2]
                end
                if  Defaults.saRatip() && Defaults.saGG() && true
                    if  length(waR[1]) != 0     println(  ">> Angular coeffients from Ratip2013   = $(waR[1]) ")    end
                    if  length(waG1)   != 0     println("\n>> Angular coeffients from SpinAngular = $waG1 ")        end
                    if  length(waR[2]) != 0     println(  ">> Angular coeffients from Ratip2013   = $(waR[2]) ")    end
                    if  length(waG2)   != 0     println("\n>> Angular coeffients from SpinAngular = $waG2 ")        end
                end
                #
                me = 0.
                for  coeff in wa[1]
                    jj = Basics.subshell_2j(basis.orbitals[coeff.a].subshell)
                    me = me + coeff.T * sqrt( jj + 1) * RadialIntegrals.GrantIab(basis.orbitals[coeff.a], basis.orbitals[coeff.b], grid, potential)
                end

                for  coeff in wa[2]
                    if  typeof(settings.eeInteractionCI)  in  [CoulombInteraction, CoulombBreit]
                        me = me + coeff.V * InteractionStrength.XL_Coulomb_DH(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                                        basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid, 
                                                                                        1/plasmaModel.debyeLength)   end
                end
                matrix[r,s] = me
            end
        end 
    else
        error("Unsupported plasma model = $(plasmaSettings.plasmaModel)")
    end
    
    if printout    println("   ... done.")    end

    return( matrix )
end



"""
`Basics.compute("radial orbital: NR, Bunge (1993)", subshell::Subshell, Z::Int64)`  
    ... to compute a radial orbital::Orbital for the given subshell and nuclear charge by using the Roothan-Hartree-Fock data by 
        Bunge et al., Atomic Data and Nuclear Data Tables 53 (1993) 113, as obtained for a non-relativistic RHF computation of the 
        neutral atom. These functions are used a large component, and the small component is obtained from the kinetic balance 
        condition. Radial orbitals can be obtained for the ground-state configuration for all elements with 2 <= Z <= 54. 
`Basics.compute("radial orbital: NR, McLean (1981)", subshell::Subshell, Z::Int64)`  
    ... to compute a radial orbital::Orbital for the given subshell and nuclear charge by using the Roothan-Hartree-Fock data by 
        McLean and McLean, Atomic Data and Nuclear Data Tables 26 (1981) 197., as obtained for a non-relativistic RHF computation of 
        the neutral atom. These functions are used a large component, and the small component is obtained from the kinetic balance 
        condition. Radial orbitals can be obtained for the ground-state configuration for all elements with 55 <= Z <= 92. 
"""
function Basics.compute(sa::String, subshell::Subshell, Z::Int64)
    if      sa == "radial orbital: NR, Bunge (1993)"          wa = Radial.OrbitalBunge1993(subshell,Z)
    elseif  sa == "radial orbital: NR, McLean (1981)"         wa = Radial.OrbitalMcLean1981(subshell,Z)
    elseif  sa == "radial orbital: hydrogenic"                wa = Radial.OrbitalHydrogenic(subshell,Z)
    elseif  sa == "radial orbital: Thomas-Fermi"              wa = Radial.OrbitalThomasFermi(subshell,Z)
    else    error("Unsupported keystring = $sa ")
    end

    return( wa )
end


#==
"""
`Basics.compute(scField::Basics.AbstractScField, grid::Radial.Grid, level::Level)`  
    ... to compute a (radial) SCF potential of type scField::AbstractScField from the given list of orbitals. 
        A potential::RadialPotential is returned. 
"""
function Basics.compute(scField::Basics.AbstractScField, grid::Radial.Grid, level::Level)
    
    if      sa == "radial potential: core-Hartree"               wa = Basics.computePotentialCoreHartree(grid, level)
    elseif  sa == "radial potential: Hartree"                    wa = Basics.computePotentialHartree(grid, level)
    elseif  sa == "radial potential: Hartree-Slater"             wa = Basics.computePotentialHartreeSlater(grid, level)
    elseif  sa == "radial potential: Kohn-Sham"                  wa = Basics.computePotentialKohnSham(grid, level)
    elseif  sa == "radial potential: Dirac-Fock-Slater"          wa = Basics.computePotentialDFS(grid, level) 
    elseif  sa == "radial potential: extended-Hartree"           wa = Basics.computePotentialExtendedHartree(grid, level)
    else    error("Unsupported keystring = $sa ")
    end

    return( wa )
end

function Basics.compute(sa::String, grid::Radial.Grid, basis::Basis)
    if      sa == "radial potential: Dirac-Fock-Slater"          wa = Basics.computePotentialDFS(grid, basis) 
    else    error("Unsupported keystring = $sa ")
    end

    return( wa )
end  ==#


"""
`Basics.computeDensity(level::Level, grid::Radial.Grid)`  
    ... computes the electronic density of level; a rho::Array{Float64,1} is returned.
"""
function Basics.computeDensity(level::Level, grid::Radial.Grid)
    basis = level.basis;    npoints = grid.NoPoints;    wx = zeros( npoints );    rho = zeros( npoints )
    # Sum_a ...
    for  a in basis.subshells
        occa = Basics.computeMeanSubshellOccupation(a, [level])
        orba = basis.orbitals[a];   nrho = length(orba.P);      rhoaa  = zeros(nrho)
        for  i = 1:nrho    rhoaa[i] = orba.P[i]^2 + orba.Q[i]^2    end
        for  i = 1:nrho    rho[i]   = rho[i] + occa * rhoaa[i]     end
    end
    
    return( rho )
end



"""
`Basics.computeDiracEnergy(sh::Subshell, Z::Float64)`  
    ... computes the Dirac energy for the hydrogenic subshell sh and for a point-like nucleus with nuclear charge Z; 
        a energy::Float64 in atomic units and without the rest energy of the electron is returned. That is the binding 
        energy of a 1s_1/2 electron for Z=1 is -0.50000665.
"""
function Basics.computeDiracEnergy(sh::Subshell, Z::Float64)
    if  Z <= 0.1    error("Requires nuclear charge Z >= 0.1")    end
    # Compute the energy from the Dirac formula
    jPlusHalf = (Basics.subshell_2j(sh) + 1) / 2;   nr = sh.n - jPlusHalf;    alpha = Defaults.getDefaults("alpha") 
    wa = sqrt(jPlusHalf^2 - Z^2 * alpha^2 )  
    wa = sqrt(1.0 + Z^2 * alpha^2 / (nr + wa)^2)
    wa = Defaults.getDefaults("speed of light: c")^2 * (1/wa - 1.0)
    return( wa )
end



"""
`Basics.computeMeanSubshellOccupation(sh::Subshell, levels::Array{Level,1})`  
    ... computes the mean subshell occupation for the subshell sh and for the given levels; a q::Float64 is returned.
"""
function Basics.computeMeanSubshellOccupation(sh::Subshell, levels::Array{Level,1})
    q = 0.
    if  length(levels) < 1   error("stop a")    end
    for  level in levels
        subshells = level.basis.subshells;    nsh = 0
        for  i = 1:length(level.basis.subshells)    if  sh == subshells[i]   nsh = i;   break   end    end
        if   nsh == 0   error("Subshell not found in CSF basis.")   end
        #
        for  i = 1:length(level.basis.csfs)   q = q + abs(level.mc[i])^2 * level.basis.csfs[i].occupation[nsh]    end
    end
    return( q/length(levels) )
end



"""
`Basics.computeMeanSubshellOccupation(sh::Subshell, basis::Basis)`  
    ... computes the mean subshell occupation for the subshell sh and for the given CSF in the basis; a q::Float64 is returned.
"""
function Basics.computeMeanSubshellOccupation(sh::Subshell, basis::Basis)
    q = 0.
    for  csf in basis.csfs
        subshells = basis.subshells;    nsh = 0
        for  i = 1:length(basis.subshells)    if  sh == subshells[i]   nsh = i;   break   end    end
        if   nsh == 0   error("Subshell not found in CSF basis.")   end
        #
        for  i = 1:length(basis.csfs)   q = q + basis.csfs[i].occupation[nsh]    end
    end
    return( q/length(basis.csfs) )
end



"""
`Basics.computeMultipletForGreenApproach(approach::AtomicState.SingleCSFwithoutCI, basis::Basis, nModel::Nuclear.Model, grid::Radial.Grid, 
                                            asfSettings::AsfSettings, greenSettings::GreenSettings; printout::Bool=false)`  
    ... computes the (Green channel) multiplet from the given basis with the SingleCSFwithoutCI approach.
"""
function Basics.computeMultipletForGreenApproach(approach::AtomicState.SingleCSFwithoutCI, basis::Basis, nModel::Nuclear.Model, 
                                                    grid::Radial.Grid, asfSettings::AsfSettings, greenSettings::GreenSettings; printout::Bool=false)
    print("Compute a Green function multiplet in $approach approach ... ")
    # In the SingleCSFwithoutCI, only the diagonal ME are included into the Hamiltonian matrix
    # (1) Check that all CSF have same (level) symmetry; issue an error if not
    sym = LevelSymmetry(0, Basics.plus)
    for  (r, csf)  in enumerate(basis.csfs)
        if     r == 1   sym  = LevelSymmetry( csf.J, csf.parity )
        elseif          sym != LevelSymmetry( csf.J, csf.parity )  error("stop a")
        end
    end
    
    # (2) Compute Hamiltonian matrix with only diagonal ME with given settings; issue a message which flags are considered.
    if  asfSettings.qedModel != NoneQed()   println("   ++ No QED terms included for this Green function approach.")               end
    if  asfSettings.jjLS.makeIt             println("   ++ No jj-LS transformation included for this Green function approach.")    end
    potential = Nuclear.nuclearPotential(nModel, grid)
    ncsf      = length(basis.csfs);    matrix = zeros(Float64, ncsf, ncsf)
    for  (r, csf)  in enumerate(basis.csfs)
        # Calculate the spin-angular coefficients
        #== 
        wa = compute("angular coefficients: e-e, Ratip2013", basis.csfs[r], basis.csfs[r])
        if  Defaults.saRatip()
            waR = compute("angular coefficients: e-e, Ratip2013", basis.csfs[r], basis.csfs[r])
            wa  = waR       
        end  ==#
        if  Defaults.saGG()
            subshellList = basis.subshells
            opa  = SpinAngular.OneParticleOperator(0, plus, true)
            waG1 = SpinAngular.computeCoefficients(opa, basis.csfs[r], basis.csfs[r], subshellList) 
            opa  = SpinAngular.TwoParticleOperator(0, plus, true)
            waG2 = SpinAngular.computeCoefficients(opa, basis.csfs[r], basis.csfs[r], subshellList)
            wa   = [waG1, waG2]
        end
        #==
        if  Defaults.saRatip() && Defaults.saGG() && true
            if  length(waR[1]) != 0     println(  ">> Angular coeffients from Ratip2013   = $(waR[1]) ")    end
            if  length(waG1)   != 0     println("\n>> Angular coeffients from SpinAngular = $waG1 ")        end
            if  length(waR[2]) != 0     println(  ">> Angular coeffients from Ratip2013   = $(waR[2]) ")    end
            if  length(waG2)   != 0     println("\n>> Angular coeffients from SpinAngular = $waG2 ")        end
        end  ==#
        #
        me = 0.
        for  coeff in wa[1]
            jj = Basics.subshell_2j(basis.orbitals[coeff.a].subshell)
            me = me + coeff.T * sqrt( jj + 1) * RadialIntegrals.GrantIab(basis.orbitals[coeff.a], basis.orbitals[coeff.b], grid, potential)
        end

        for  coeff in wa[2]
            if  asfSettings.coulombCI    
                me = me + coeff.V * InteractionStrength.XL_Coulomb(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                                basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid)   end
            if  asfSettings.breitCI
                me = me + coeff.V * InteractionStrength.XL_Breit(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                            basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid,
                                                                            Basics.CoulombBreit(0.))     end
        end
        matrix[r,r] = me
    end
    
    # (3) Diagonalize matrix with Julia;   assign a multiplet 
    eigen  = Basics.diagonalize("matrix: LinearAlgebra", matrix)
    levels = Level[]
    for  ev = 1:length(eigen.values)
        level = Level( sym.J, AngularM64(sym.J.num//sym.J.den), sym.parity, 0, eigen.values[ev], 0., true, basis, eigen.vectors[ev] ) 
        push!( levels, level)
    end
    println("done with $(length(levels)) levels")
    
    multiplet = Multiplet("SingleCSFwithoutCI multiplet for $sym", levels)
    return( multiplet )
end



"""
`Basics.computeMultipletForGreenApproach(approach::AtomicState.CoreSpaceCI, basis::Basis, nModel::Nuclear.Model, grid::Radial.Grid, 
                                            asfSettings::AsfSettings, greenSettings::GreenSettings; printout::Bool=false)`  
    ... computes the (Green channel) multiplet from the given basis with the CoreSpaceCI approach in which the electron-electron 
        interaction is taken into account only between the bound-state orbitals.
"""
function Basics.computeMultipletForGreenApproach(approach::AtomicState.CoreSpaceCI, basis::Basis, nModel::Nuclear.Model, 
                                                    grid::Radial.Grid, asfSettings::AsfSettings, greenSettings::GreenSettings; printout::Bool=false)
    print("Compute a Green function multiplet in $approach approach ... ")
    # In the SingleCSFwithoutCI, only the diagonal ME are included into the Hamiltonian matrix
    # (1) Check that all CSF have same (level) symmetry; issue an error if not
    sym = LevelSymmetry(0, Basics.plus)
    for  (r, csf)  in enumerate(basis.csfs)
        if     r == 1   sym  = LevelSymmetry( csf.J, csf.parity )
        elseif          sym != LevelSymmetry( csf.J, csf.parity )  error("stop a")
        end
    end
    
    # (2) Compute Hamiltonian matrix with only diagonal ME with given settings; issue a message which flags are considered.
    if  asfSettings.qedModel != NoneQed()   println("   ++ No QED terms included for this Green function approach.")               end
    if  asfSettings.jjLS.makeIt             println("   ++ No jj-LS transformation included for this Green function approach.")    end
    potential = Nuclear.nuclearPotential(nModel, grid)
    ncsf      = length(basis.csfs);    matrix = zeros(Float64, ncsf, ncsf)
    for  (r, rcsf)  in enumerate(basis.csfs)
        for  (s, scsf)  in enumerate(basis.csfs)
            # Calculate the spin-angular coefficients
            if  Defaults.saGG()
                subshellList = basis.subshells
                opa  = SpinAngular.OneParticleOperator(0, plus, true)
                waG1 = SpinAngular.computeCoefficients(opa, basis.csfs[r], basis.csfs[s], subshellList) 
                opa  = SpinAngular.TwoParticleOperator(0, plus, true)
                waG2 = SpinAngular.computeCoefficients(opa, basis.csfs[r], basis.csfs[s], subshellList)
                wa   = [waG1, waG2]
            end
            #==
            if  Defaults.saRatip() && Defaults.saGG() && true
                if  length(waR[1]) != 0     println(  ">> Angular coeffients from Ratip2013   = $(waR[1]) ")    end
                if  length(waG1)   != 0     println("\n>> Angular coeffients from SpinAngular = $waG1 ")        end
                if  length(waR[2]) != 0     println(  ">> Angular coeffients from Ratip2013   = $(waR[2]) ")    end
                if  length(waG2)   != 0     println("\n>> Angular coeffients from SpinAngular = $waG2 ")        end
            end  ==#
            #
            me = 0.
            for  coeff in wa[1]
                jj = Basics.subshell_2j(basis.orbitals[coeff.a].subshell)
                if  basis.orbitals[coeff.a].isBound   &&   basis.orbitals[coeff.b].isBound
                    me = me + coeff.T * sqrt( jj + 1) * RadialIntegrals.GrantIab(basis.orbitals[coeff.a], basis.orbitals[coeff.b], grid, potential)
                end
            end

            for  coeff in wa[2]
                if  basis.orbitals[coeff.a].isBound   &&   basis.orbitals[coeff.b].isBound   && 
                    basis.orbitals[coeff.c].isBound   &&   basis.orbitals[coeff.d].isBound
                    if  asfSettings.coulombCI    
                        me = me + coeff.V * InteractionStrength.XL_Coulomb(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                                        basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid)   end
                    if  asfSettings.breitCI
                        me = me + coeff.V * InteractionStrength.XL_Breit(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                                    basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid,
                                                                                    Basics.CoulombBreit(0.))                                    end
                end
            end
            matrix[r,s] = me
        end
    end
    
    # (3) Diagonalize matrix with Julia;   assign a multiplet 
    eigen  = Basics.diagonalize("matrix: LinearAlgebra", matrix)
    levels = Level[]
    for  ev = 1:length(eigen.values)
        level = Level( sym.J, AngularM64(sym.J.num//sym.J.den), sym.parity, 0, eigen.values[ev], 0., true, basis, eigen.vectors[ev] ) 
        push!( levels, level)
    end
    println("done with $(length(levels)) levels")
    
    multiplet = Multiplet("CoreSpaceCI multiplet for $sym", levels)
    return( multiplet )
end



"""
`Basics.computeMultipletForGreenApproach(approach::AtomicState.DampedSpaceCI, basis::Basis, nModel::Nuclear.Model, grid::Radial.Grid,
                                            asfSettings::AsfSettings, greenSettings::GreenSettings; printout::Bool=false)`  
    ... computes the (Green channel) multiplet from the given basis with the DampedSpaceCI approach in which the electron-electron 
        interaction strength are all damped by the factor exp(- dampingTau * r).
"""
function Basics.computeMultipletForGreenApproach(approach::AtomicState.DampedSpaceCI, basis::Basis, nModel::Nuclear.Model, 
                                                    grid::Radial.Grid, asfSettings::AsfSettings, greenSettings::GreenSettings; 
                                                    printout::Bool=false)
    print("Compute a Green function multiplet in $approach approach ... ")
    # In the SingleCSFwithoutCI, only the diagonal ME are included into the Hamiltonian matrix
    # (1) Check that all CSF have same (level) symmetry; issue an error if not
    sym = LevelSymmetry(0, Basics.plus)
    for  (r, csf)  in enumerate(basis.csfs)
        if     r == 1   sym  = LevelSymmetry( csf.J, csf.parity )
        elseif          sym != LevelSymmetry( csf.J, csf.parity )  error("stop a")
        end
    end
    
    # (2) Compute Hamiltonian matrix with only diagonal ME with given settings; issue a message which flags are considered.
    if  asfSettings.qedModel != NoneQed()   println("   ++ No QED terms included for this Green function approach.")               end
    if  asfSettings.jjLS.makeIt             println("   ++ No jj-LS transformation included for this Green function approach.")    end
    potential = Nuclear.nuclearPotential(nModel, grid)
    ncsf      = length(basis.csfs);    matrix = zeros(Float64, ncsf, ncsf);    tau = greenSettings.dampingTau
    for  (r, rcsf)  in enumerate(basis.csfs)
        for  (s, scsf)  in enumerate(basis.csfs)
            # Calculate the spin-angular coefficients
            if  Defaults.saGG()
                subshellList = basis.subshells
                opa  = SpinAngular.OneParticleOperator(0, plus, true)
                waG1 = SpinAngular.computeCoefficients(opa, basis.csfs[r], basis.csfs[s], subshellList) 
                opa  = SpinAngular.TwoParticleOperator(0, plus, true)
                waG2 = SpinAngular.computeCoefficients(opa, basis.csfs[r], basis.csfs[s], subshellList)
                wa   = [waG1, waG2]
            end
            #
            me = 0.
            for  coeff in wa[1]
                jj = Basics.subshell_2j(basis.orbitals[coeff.a].subshell)
                me = me + coeff.T * sqrt( jj + 1) * RadialIntegrals.GrantIabDamped(tau, basis.orbitals[coeff.a], basis.orbitals[coeff.b], grid, potential)
            end

            for  coeff in wa[2]
                if  typeof(asfSettings.eeInteraction) in [CoulombInteraction, CoulombBreit]   
                    me = me + coeff.V * InteractionStrength.XL_CoulombDamped(tau, coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                                            basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid)   end
                if  typeof(asfSettings.eeInteraction) in [BreitInteraction, CoulombBreit]
                    me = me + coeff.V * InteractionStrength.XL_BreitDamped(tau, coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                                            basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid,
                                                                                            asfSettings.eeInteraction)                                  end
            end
            matrix[r,s] = me
        end
    end
    
    # (3) Diagonalize matrix with Julia;   assign a multiplet 
    eigen  = Basics.diagonalize("matrix: LinearAlgebra", matrix)
    levels = Level[]
    for  ev = 1:length(eigen.values)
        level = Level( sym.J, AngularM64(sym.J.num//sym.J.den), sym.parity, 0, eigen.values[ev], 0., true, basis, eigen.vectors[ev] ) 
        push!( levels, level)
    end
    println("done with $(length(levels)) levels")
    
    multiplet = Multiplet("DampedSpaceCI multiplet for $sym", levels)
    return( multiplet )
end


"""
`Basics.computePotential(scField::Basics.AaHSField, grid::Radial.Grid, orbitals::Dict{Subshell, Orbital}, mu::Float64, temp::Float64)`  
    ... to compute a (radial) Hartree-Slater-type potential for the density as modeled by the given set of orbital
        and their occupation numbers in the atomic-average model; a potential::RadialPotential is returned. 
        Here, the atomic-average HS potential is defined by
                                                    
                                        (2ja+1)                       3  [     3          ]^1/3      r
            V_HS(r) = = - Sum_a  --------------------  Y0_aa(r)  +  -  [ ---------  rho ]       *  -
                                    [exp(eps-mu/kT) + 1]               2  [ 4pi^2 r^2      ]          2

        with rho(r) = Sum_a  (2j_a+1) * (Pa^2(r) + Qa^2(r)) / [1 + exp(epsilon-mu)/kT]   ... charge density of all electrons.  
        An Radial.Potential with V_aa_HS(r)  is returned that is consistent with an effective charge Z(r).
"""
function Basics.computePotential(scField::Basics.AaHSField, grid::Radial.Grid, orbitals::Dict{Subshell, Orbital}, mu::Float64, temp::Float64)
    npoints = grid.NoPoints;    wx = zeros( npoints );    rho = zeros( npoints )
    # Sum_a ...
    for (k,v)  in  orbitals 
        occa = Basics.twice( Basics.subshell_j(k)) + 1
        occa = occa * Basics.FermiDirac(v.energy, mu, temp)
        nrho = length(v.P);      rhoaa  = zeros(nrho)
        for  i = 1:nrho    rhoaa[i] = v.P[i]^2 + v.Q[i]^2    end
        for  i = 1:nrho    wx[i]    = wx[i]  - occa * RadialIntegrals.Yk_ab(0, grid.r[i], rhoaa, nrho, grid)   end
        for  i = 1:nrho    rho[i]   = rho[i] + occa * rhoaa[i]     end
        Yk = RadialIntegrals.Yk_ab(0, grid.r[nrho], rhoaa, nrho, grid)
        for  i = nrho+1:npoints       wx[i] = wx[i] - occa * Yk    end
    end
    ## for  i = 2:npoints     rho[i]   = rho[i] / (4pi * grid.r[i])   end;       rho[1] = 0.
    for  i = 2:npoints     wx[i]    = wx[i] + (3/2) * (3/(4pi^2 * grid.r[i]^2) * rho[i])^(1/3) * grid.r[i] / 2   end
    #
    wc = Radial.Potential("average-atom HS", wx, grid)
    return( wc )
end



"""
`Basics.computePotential(scField::Basics.AaDFSField, grid::Radial.Grid, orbitals::Dict{Subshell, Orbital}, mu::Float64, temp::Float64)`  
    ... to compute a (radial) average-atom DFS-type potential for  the density as modeled by the given set of orbital
        and their occupation numbers in the atomic-average model; a potential::RadialPotential is returned. 
        The atomic-average DFS potential is defined by

                                                            [   3                   ]^(1/3)
            V_DFS(r) = int_0^infty dr'  rho_t(r') / r_>  -- [-----------   rho_t(r) ]
                                                            [ 4 pi^2 r^2            ]

        but with the occupation numbers (2j_a+1) / [exp(eps-mu/kT) + 1]
        with r_> = max(r,r')  and  rho_t(r) = sum_a (Pa^2(r) + Qa^2(r))   ... charge density of all electrons.  
        An Radial.Potential with -r * V_DFS(r)  is returned to be consistent with an effective charge Z(r).
"""
function Basics.computePotential(scField::Basics.AaDFSField, grid::Radial.Grid, orbitals::Dict{Subshell, Orbital}, mu::Float64, temp::Float64)
    npoints = grid.NoPoints;    rhot = zeros( npoints );     wb = zeros( npoints );    wx = zeros( npoints )
    # Compute the charge density with the average-atom occupation numbers
    for (k,v)  in  orbitals 
        occ  = Basics.twice( Basics.subshell_j(k)) + 1
        occ  = occ * Basics.FermiDirac(v.energy, mu, temp)
        nrho = length(v.P)
        for    i = 1:nrho   rhot[i] = rhot[i] + occ * (v.P[i]^2 + v.Q[i]^2)    end
    end
    # Define the integrant and take the integrals
    for i = 1:npoints
        for j = 1:npoints    rg = max( grid.r[i],  grid.r[j]);    wx[j] = rhot[j]/rg   end
        wb[i] = RadialIntegrals.V0(wx, npoints, grid::Radial.Grid)
        wb[i] = wb[i] - (3 /(4*pi^2 * grid.r[i]^2) * rhot[i])^(1/3)
    end
    # Define the potential with regard to Z(r)
    for i = 2:npoints    wb[i] = - wb[i] * grid.r[i]   end;    wb[1] = 0.
    #
    wc = Radial.Potential("average-atom DFS", wb, grid)
    return( wc )
end



"""
`Basics.computePotential(scField::Basics.CHField, grid::Radial.Grid, level::Level)`  
    ... to compute a (radial) core-Hartree potential for the given level; a potential::RadialPotential is returned. 
        The core-Hartree potential is defined by

            V_CH(r) = int_0^infty dr'  rho_c(r') / r_>                 
            
        with r_> = max(r,r')  and  rho_c(r) =  sum_a (Pa^2(r) + Qa^2(r)) the charge density of the core-electrons.  
        An Radial.Potential with -r * V_CH(r)  is returned to be consistent with an effective charge Z(r).
"""
function Basics.computePotential(scField::Basics.CHField, grid::Radial.Grid, level::Level)
    basis = level.basis;    npoints = grid.NoPoints
    println("coreSubshells = $(basis.coreSubshells);    subshells = $(basis.subshells)")
    rhoc = zeros( npoints );    wx = zeros( npoints );    wb = zeros( npoints )
    # Compute the charge density of the core orbitals for the given level
    for  sh in basis.coreSubshells
    ## for  sh in basis.subshells
        orb  = basis.orbitals[sh]
        occ  = Basics.computeMeanSubshellOccupation(sh, [level])    ## (Basics.subshell_2j(sh) + 1.0)  for closed-shells
        nrho = length(orb.P)
        for    i = 1:nrho   rhoc[i] = rhoc[i] + occ * (orb.P[i]^2 + orb.Q[i]^2)    end
    end
    # Define the integrant and take the integrals
    for i = 1:npoints
        for j = 1:npoints    rg = max( grid.r[i],  grid.r[j]);    wx[j] = rhoc[j]/rg   end
        wb[i] = RadialIntegrals.V0(wx, npoints, grid::Radial.Grid)
    end
    # Define the potential with regard to Z(r)
    for i = 1:npoints    wb[i] = - wb[i] * grid.r[i]   end
    #
    wc = Radial.Potential("CoreHartree", wb, grid)
    return( wc )
end



"""
`Basics.computePotential(scField::Basics.HartreeField, grid::Radial.Grid, level::Level)`  
    ... to compute a (radial) Hartree potential for the given level; a potential::RadialPotential is returned. 
        Here, the Hartree potential is simply defined by

            V_H(r) = - Sum_a  q_a^bar  Y0_aa(r)

        An Radial.Potential with V_H(r)  is returned that is consistent with an effective charge Z(r).
"""
function Basics.computePotential(scField::Basics.HartreeField, grid::Radial.Grid, level::Level)
    basis = level.basis;    npoints = grid.NoPoints;    wx = zeros( npoints )
    # Sum_a ...
    for  a in basis.subshells
        occa = Basics.computeMeanSubshellOccupation(a, [level])
        orba = basis.orbitals[a];   nrho = length(orba.P);      rhoaa  = zeros(nrho)
        for  i = 1:nrho    rhoaa[i] = orba.P[i]^2 + orba.Q[i]^2    end
        for  i = 1:nrho    wx[i]    = wx[i] - occa * RadialIntegrals.Yk_ab(0, grid.r[i], rhoaa, nrho, grid)   end
        Yk = RadialIntegrals.Yk_ab(0, grid.r[nrho], rhoaa, nrho, grid)
        for  i = nrho+1:npoints       wx[i] = wx[i] - occa * Yk    end
    end
    
    wc = Radial.Potential("Hartree", wx, grid)
    return( wc )
end


"""
`Basics.computePotential(scField::Basics.HSField, grid::Radial.Grid, level::Level)`  
    ... to compute a (radial) Hartree-Slater potential for the given level; a potential::RadialPotential is returned. 
        Here, the Hartree-Slater potential is defined by
                                                    
                                                    3  [     3          ]^1/3      r
            V_HS(r) = = - Sum_a  q_a^bar  Y0_aa(r)  +  -  [ ---------  rho ]       *  -
                                                    2  [ 4pi^2 r^2      ]          2

        with rho(r) = Sum_a  q_a^bar * (Pa^2(r) + Qa^2(r))   ... charge density of all electrons.  
        An Radial.Potential with V_HS(r)  is returned that is consistent with an effective charge Z(r).
"""
function Basics.computePotential(scField::Basics.HSField, grid::Radial.Grid, level::Level)
    basis = level.basis;    npoints = grid.NoPoints;    wx = zeros( npoints );    rho = zeros( npoints )
    # Sum_a ...
    for  a in basis.subshells
        occa = Basics.computeMeanSubshellOccupation(a, [level])
        orba = basis.orbitals[a];   nrho = length(orba.P);      rhoaa  = zeros(nrho)
        for  i = 1:nrho    rhoaa[i] = orba.P[i]^2 + orba.Q[i]^2    end
        for  i = 1:nrho    wx[i]    = wx[i]  - occa * RadialIntegrals.Yk_ab(0, grid.r[i], rhoaa, nrho, grid)   end
        for  i = 1:nrho    rho[i]   = rho[i] + occa * rhoaa[i]     end
        Yk = RadialIntegrals.Yk_ab(0, grid.r[nrho], rhoaa, nrho, grid)
        for  i = nrho+1:npoints       wx[i] = wx[i] - occa * Yk    end
    end
    ## for  i = 2:npoints     rho[i]   = rho[i] / (4pi * grid.r[i])   end;       rho[1] = 0.
    for  i = 2:npoints     wx[i]    = wx[i] + (3/2) * (3/(4pi^2 * grid.r[i]^2) * rho[i])^(1/3) * grid.r[i] / 2   end
    #
    wc = Radial.Potential("Hartree-Slater", wx, grid)
    return( wc )
end



"""
`Basics.computePotential(scField::Basics.KSField, grid::Radial.Grid, level::Level)`  
    ... to compute a (radial) Kohn-Sham potential for the given level; a potential::RadialPotential is returned. 
        The Kohn-Sham potential is defined by

                                                            2     [   81                   ]^(1/3)
            V_KS(r) =  int_0^infty dr'  rho_t(r') / r_>  --  ---  [--------   r * rho_t(r) ]
                                                            3 r   [ 32 pi^2                ]

        with r_> = max(r,r')  and  rho_t(r) = sum_a (Pa^2(r) + Qa^2(r))   ... charge density of all electrons.  
        An Radial.Potential with -r * V_KS(r)  is returned to be consistent with an effective charge Z(r).
"""
function Basics.computePotential(scField::Basics.KSField, grid::Radial.Grid, level::Level)
    basis = level.basis;       npoints = grid.NoPoints
    rhot = zeros( npoints );   wb = zeros( npoints );   wx = zeros( npoints )
    # Compute the charge density of the core orbitals for the given level
    for  sh in basis.subshells
        orb  = basis.orbitals[sh]
        occ  = Basics.computeMeanSubshellOccupation(sh, [level])
        nrho = length(orb.P)
        for    i = 1:nrho   rhot[i] = rhot[i] + occ * (orb.P[i]^2 + orb.Q[i]^2)    end
    end
    # Define the integrant and take the integrals (without alpha)
    for i = 1:npoints
        for j = 1:npoints    rg = max( grid.r[i],  grid.r[j]);    wx[j] = rhot[j]/rg   end
        wb[i] = RadialIntegrals.V0(wx, npoints, grid::Radial.Grid)
        wb[i] = wb[i] - 2 / (3grid.r[i]) * (81 /(32 * pi^2) * grid.r[i] * rhot[i])^(1/3)
    end
    # Define the potential with regard to Z(r)
    for i = 2:npoints    wb[i] = - wb[i] * grid.r[i]   end;    wb[1] = 0.
    #
    wc = Radial.Potential("Kohn-Sham", wb, grid)
    return( wc )
end



"""
`Basics.computePotential(scField::Basics.DFSField, grid::Radial.Grid, level::Level)`  
    ... to compute a (radial) Dirac-Fock-Slater potential for the given level; a potential::RadialPotential is returned. 
        The Dirac-Fock-Slater potential is defined by

                                                            [   3                   ]^(1/3)
            V_DFS(r) = int_0^infty dr'  rho_t(r') / r_>  -- [-----------   rho_t(r) ]
                                                            [ 4 pi^2 r^2            ]

        with r_> = max(r,r')  and  rho_t(r) = sum_a (Pa^2(r) + Qa^2(r))   ... charge density of all electrons.  
        An Radial.Potential with -r * V_DFS(r)  is returned to be consistent with an effective charge Z(r).
"""
function Basics.computePotential(scField::Basics.DFSField, grid::Radial.Grid, level::Level)
    basis = level.basis;    npoints = grid.NoPoints
    rhot = zeros( npoints );    wb = zeros( npoints );    wx = zeros( npoints )
    # Compute the charge density of the core orbitals for the given level
    for  sh in basis.subshells
        orb  = basis.orbitals[sh]
        occ  = Basics.computeMeanSubshellOccupation(sh, [level])
        nrho = length(orb.P)
        for    i = 1:nrho   rhot[i] = rhot[i] + occ * (orb.P[i]^2 + orb.Q[i]^2)    end
    end
    # Define the integrant and take the integrals
    ## println(">>>>> wy = $(scField.strength) * DFS potential");     
    wy = scField.strength
    for i = 1:npoints
        for j = 1:npoints    rg = max( grid.r[i],  grid.r[j]);    wx[j] = rhot[j]/rg   end
        wb[i] = RadialIntegrals.V0(wx, npoints, grid::Radial.Grid)
        wb[i] = wb[i] - (3 /(4*pi^2 * grid.r[i]^2) * rhot[i])^(1/3) * wy
    end
    ## wb = 1.02 * wb
    # Define the potential with regard to Z(r)
    for i = 2:npoints    wb[i] = - wb[i] * grid.r[i]   end;    wb[1] = 0.
    #
    wc = Radial.Potential("DFS", wb, grid)
    return( wc )
end



"""
`Basics.computePotential(scField::Basics.DFSField, grid::Radial.Grid, basis::Basis)`  
    ... to compute the same but for the mean occupation of the orbitals in the given basis.
"""
function Basics.computePotential(scField::Basics.DFSField, grid::Radial.Grid, basis::Basis)
    npoints = grid.NoPoints
    rhot = zeros( npoints );    wb = zeros( npoints );    wx = zeros( npoints )
    # Compute the charge density of the core orbitals for the given level
    for  sh in basis.subshells
        orb  = basis.orbitals[sh]
        occ  = Basics.computeMeanSubshellOccupation(sh, basis)
        nrho = length(orb.P)
        for    i = 1:nrho   rhot[i] = rhot[i] + occ * (orb.P[i]^2 + orb.Q[i]^2)    end
    end
    # Define the integrant and take the integrals
    for i = 1:npoints
        for j = 1:npoints    rg = max( grid.r[i],  grid.r[j]);    wx[j] = rhot[j]/rg   end
        wb[i] = RadialIntegrals.V0(wx, npoints, grid::Radial.Grid)
        wb[i] = wb[i] - (3 /(4*pi^2 * grid.r[i]^2) * rhot[i])^(1/3)
    end
    # Define the potential with regard to Z(r)
    for i = 2:npoints    wb[i] = - wb[i] * grid.r[i]   end;    wb[1] = 0.
    #
    wc = Radial.Potential("DFS for CSF basis", wb, grid)
    return( wc )
end


"""
`Basics.computePotential(scField::Basics.DFSwCPField, kappa::Int64, cp::CorePolarization, grid::Radial.Grid, level::Level)`  
    ... to compute a (radial) Dirac-Fock-Slater potential for the given level; 
        a potential::RadialPotential is returned. The Dirac-Fock-Slater potential with core-polarization is defined by
                                            
            V_DFS_wCP(r) = V_DFS + V_p (r) delta_a,valence        with
            
                                alpha_d  r^2
            V_p (r)      = - ---------------
                            2 (r^2 + r_c)^3
            
        and where alpha_d is the core-polarizibility and r_c the core radius.
        An Radial.Potential with -r * V_DFS(r)  is returned to be consistent with an effective charge Z(r).
"""
function Basics.computePotential(scField::Basics.DFSwCPField, kappa::Int64, cp::CorePolarization, grid::Radial.Grid, level::Level)
    basis = level.basis;    npoints = grid.NoPoints
    rhot = zeros( npoints );    wb = zeros( npoints );    wx = zeros( npoints )
    # Compute the charge density of the core orbitals for the given level
    for  sh in basis.subshells
        orb  = basis.orbitals[sh]
        occ  = Basics.computeMeanSubshellOccupation(sh, [level])
        nrho = length(orb.P)
        for    i = 1:nrho   rhot[i] = rhot[i] + occ * (orb.P[i]^2 + orb.Q[i]^2)    end
    end
    # Define the integrant and take the integrals
    for i = 1:npoints
        for j = 1:npoints    rg = max( grid.r[i],  grid.r[j]);    wx[j] = rhot[j]/rg   end
        wb[i] = RadialIntegrals.V0(wx, npoints, grid::Radial.Grid)
        wb[i] = wb[i] - (3 /(4*pi^2 * grid.r[i]^2) * rhot[i])^(1/3)
    end
    # Add the core-polarization contributions from the valence shells
    if  cp.doApply
        println(">>> Add core-polarization potential for kappa = $kappa.")
        for i = 1:npoints
            denom = 2.0 * (grid.r[i]^2 + cp.coreRadius^2)^3
            wb[i] = wb[i] - cp.coreAlpha * grid.r[i]^2 / denom
        end
    else
        println(">>> No core-polarization potential added for kappa = $kappa; this is an ERROR in the present implementation.")
    end
    #
    # Define the potential with regard to Z(r)
    for i = 2:npoints    wb[i] = - wb[i] * grid.r[i]   end;    wb[1] = 0.
    #
    wc = Radial.Potential("DFS", wb, grid)
    return( wc )
end



"""
`Basics.computePotential(scField::Basics.EHField, grid::Radial.Grid, level::Level)`  
    ... to compute a (radial) extended-Hartree potential for the given level; a potential::RadialPotential is returned. 
        The extended-Hartreepotential is defined by Gu (2008, Eq.(8)) and is based on the Y^k_ab (r) function and various 
        direct and exchange coefficients as the arise from average energy of the configuration.
"""
function Basics.computePotential(scField::Basics.EHField, grid::Radial.Grid, level::Level)
    basis = level.basis;    npoints = grid.NoPoints
    rhoab = zeros( npoints );    wx = zeros( npoints );    wg = zeros( npoints )
    # First term Sum_ab ...
    for  a in basis.subshells
        occa = Basics.computeMeanSubshellOccupation(a, [level])
        orba = basis.orbitals[a]
        nrho = length(orba.P);      rhoa  = zeros(nrho)
        for  i = 1:nrho    rhoa[i] = orba.P[i]^2 + orba.Q[i]^2    end
        for  b in basis.subshells
            occb = Basics.computeMeanSubshellOccupation(b, [level])
            orbb = basis.orbitals[b]
            nrho = length(orbb.P);   rhobb = zeros(nrho)
            if  a != b  wc = occa*occb   else   wc = occa * (occb -1.0)   end
            if  wc < 0.  error("Improper potential for subshells with mean occupation < 1; stop a")    end
            for  i = 1:nrho    rhobb[i] = orbb.P[i]^2 + orbb.Q[i]^2       end
            for  i = 1:length(orba.P)    wx[i] = wx[i] + wc * RadialIntegrals.Yk_ab(0, grid.r[i], rhobb, nrho, grid) * rhoa[i]   end
        end
    end
    # Second term Sum_a ...
    for  a in basis.subshells
        occa = Basics.computeMeanSubshellOccupation(a, [level])
        orba = basis.orbitals[a]
        nrho = length(orba.P);      rhoa  = zeros(nrho)
        for  i = 1:nrho    rhoa[i] = orba.P[i]^2 + orba.Q[i]^2   end
        wc = occa * (occa - 1.0)
        if  wc < 0.  error("Improper potential for subshells with mean occupation < 1; stop b")    end
        for  k = 1:10
            fk = AngularMomentum.fk_ab(a,a);    if  fk == 0   cycle    end
            for  i = 1:length(orba.P)    wx[i] = wx[i] + wc * fk * RadialIntegrals.Yk_ab(k, grid.r[i], rhoa, nrho, grid) * rhoa[i]   end
        end
    end
    # Third term Sum_a != b ...
    for  a in basis.subshells,   b in basis.subshells
        if  a == b   continue   end 
        occa = Basics.computeMeanSubshellOccupation(a, [level])
        occb = Basics.computeMeanSubshellOccupation(b, [level])
        orba = basis.orbitals[a]
        orbb = basis.orbitals[b]
        nrho = min(length(orba.P), length(orbb.P));      rhoab  = zeros(nrho)
        for  i = 1:nrho    rhoab[i] = orba.P[i]*orbb.P[i] + orba.Q[i]*orbb.Q[i]   end
        wc = occa * occb
        for  k = 0:10
            gk = AngularMomentum.gk_ab(a,b);    if  gk == 0   continue    end
            for  i = 1:nrho    wx[i] = wx[i] + wc * gk * RadialIntegrals.Yk_ab(k, grid.r[i], rhoab, nrho, grid) * rhoab[i]   end
        end
    end
    # Weight factor
    for  a in basis.subshells
        occa = Basics.computeMeanSubshellOccupation(a, [level])
        orba = basis.orbitals[a]
        nrho = length(orba.P);      rhoa  = zeros(nrho)
        for  i = 1:nrho    rhoa[i] = orba.P[i]^2 + orba.Q[i]^2   end
        wc = occa
        for  i = 1:length(orba.P)    wg[i] = wg[i] + wc * rhoa[i] * grid.r[i]   end
    end
    #
    for  i = 1:npoints   if  wg[i] > 0.    wx[i] = wx[i]/wg[i]   else   wx[i] = 0.   end    end
    return( wx )
end



"""
`Basics.computeScfCoefficients(field::Basics.ALField, basis::Basis, subsh::Subshell)`  
    ... to compute and sort out all angular coefficients that contribute to the self-consistent field of the subshell sh.
        For an average-level field, the energy functional is given by the sum of all diagonal matrix elements of the given
        basis; the function also adds all coefficients that give rise to the same XL_Coulomb interaction strengths.
        A tuple of two lists Tuple{Array{AngularTcoeff,1},Array{AngularVcoeff,1}} is returned.
"""
function Basics.computeScfCoefficients(field::Basics.ALField, basis::Basis, subsh::Subshell)
    #  Select only angular coefficients for the Coulomb interaction
    function isCoulomb(L::Int64, sha::Subshell, shb::Subshell, shc::Subshell, shd::Subshell)
        la = Basics.subshell_l(sha);    ja2 = Basics.subshell_2j(sha)
        lb = Basics.subshell_l(shb);    jb2 = Basics.subshell_2j(shb)
        lc = Basics.subshell_l(shc);    jc2 = Basics.subshell_2j(shc)
        ld = Basics.subshell_l(shd);    jd2 = Basics.subshell_2j(shd)

        if  AngularMomentum.triangularDelta(ja2+1,jc2+1,L+L+1) * AngularMomentum.triangularDelta(jb2+1,jd2+1,L+L+1) == 0   ||   
            rem(la+lc+L,2) == 1   ||   rem(lb+ld+L,2) == 1
            return( false )
        end
        return( true )
    end
    #
    Tcoeffs = AngularCoefficientsRatip2013.AngularTcoeff[];   allTcoeffs = AngularCoefficientsRatip2013.AngularTcoeff[]
    Vcoeffs = AngularCoefficientsRatip2013.AngularVcoeff[];   allVcoeffs = AngularCoefficientsRatip2013.AngularVcoeff[]
    # Compute all angular coefficients for the diagonal matrix elements
    for  csf in basis.csfs
        # Calculate the spin-angular coefficients
        if  Defaults.saGG()
            subshellList = basis.subshells
            opa  = SpinAngular.OneParticleOperator(0, plus, true)
            waG1 = SpinAngular.computeCoefficients(opa, csf, csf, subshellList) 
            opa  = SpinAngular.TwoParticleOperator(0, plus, true)
            waG2 = SpinAngular.computeCoefficients(opa, csf, csf, subshellList)
            wa   = [waG1, waG2]
        end
        #
        append!(allTcoeffs, wa[1])
        append!(allVcoeffs, wa[2])
    end
    # Determine a list reduced angular that just depend on subshell subsh
    tcoeffs = Tuple{Int64,Subshell,Subshell}[]
    vcoeffs = Tuple{Int64,Subshell,Subshell,Subshell,Subshell}[]
    for  coeff in allTcoeffs
        if   subsh == coeff.a
            wa = (coeff.nu, coeff.a, coeff.b)
            if      wa in tcoeffs
            else    push!(tcoeffs, wa)
            end
        end
    end
    for  coeff in allVcoeffs
        if   subsh == coeff.a   &&   isCoulomb(coeff.nu, coeff.a, coeff.b, coeff.c, coeff.d)
            wa = (coeff.nu, coeff.a, coeff.b, coeff.c, coeff.d)
            if      wa in vcoeffs
            else    push!(vcoeffs, wa)
            end
        end
        if   subsh == coeff.b   &&   isCoulomb(coeff.nu, coeff.b, coeff.a, coeff.d, coeff.c)
            wa = (coeff.nu, coeff.b, coeff.a, coeff.d, coeff.c)
            if      wa in vcoeffs
            else    push!(vcoeffs, wa)
            end
        end
    end
    # Now collect all contributions that belong to the selected reduced coefficients
    for  redCoeff in tcoeffs
        redT = 0.
        for  coeff  in  allTcoeffs
            wa = (coeff.nu, coeff.a, coeff.b)
            if  wa == redCoeff      redT = redT + coeff.T   end
        end
        push!( Tcoeffs, AngularCoefficientsRatip2013.AngularTcoeff( redCoeff[1], redCoeff[2], redCoeff[3], redT) )
    end
    for  redCoeff in vcoeffs
        redV = 0.
        for coeff  in  allVcoeffs
            ## V = coeff.V
            if      coeff.a == coeff.d   &&  coeff.a != coeff.b     V = coeff.V / 2.    
            elseif  coeff.b == coeff.c   &&  coeff.a != coeff.b     V = coeff.V / 2.    
            else                                                                V = coeff.V
            end
            wa = (coeff.nu, coeff.a, coeff.b, coeff.c, coeff.d);    
            if  wa == redCoeff      redV = redV + V    end
            wa = (coeff.nu, coeff.b, coeff.a, coeff.d, coeff.c)
            if  wa == redCoeff      redV = redV + V    end
        end
        push!( Vcoeffs, AngularCoefficientsRatip2013.AngularVcoeff( redCoeff[1], redCoeff[2], redCoeff[3], redCoeff[4], redCoeff[5],  redV) )
    end
    
    return( Tcoeffs, Vcoeffs )
end
