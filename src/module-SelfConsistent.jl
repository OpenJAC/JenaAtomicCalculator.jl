
"""
`module  JAC.SelfConsistent`  
	... a submodel of JAC that contains all structs and methods to generate self-consistent fields of different 
	    kind and complexity.
"""
module SelfConsistent

using  Printf, ..Basics, ..Bsplines, ..Defaults, ..Hamiltonian, ..ManyElectron, ..Nuclear, ..Radial, ..SpinAngular


"""
`SelfConsistent.computeAngularCoefficients(scField::Basics.ALField, basis::Basis)` 
    ... computes all spin-angular coefficients for the average-level (AL) functional as obtained for the given
        basis. In the AL functional, only the diagonal matrix elements (from the trace of the Hamiltonian matrix)
        are considered, and the energy is averaged with regard to the number of CSF in the basis. This makes the total 
        energy of the system comparable and independent of the size of the basis. A Tuple of two lists with one- and 
        two-particle coefficients  tpl::Tuple{coeffs1p::Array{Coefficient1p,1}, coeffs2p::Array{Coefficient2p,1}} 
        is returned.
"""
function computeAngularCoefficients(scField::Basics.ALField, basis::Basis)
    ncsf = length(basis.csfs);    coeffs1p = SpinAngular.Coefficient1p[];     coeffs2p = SpinAngular.Coefficient2p[] 
    
    # Compute angular coefficients in turn for all diagonal ME
    for  csf  in  basis.csfs
        coeffs = SpinAngular.computeCoefficientsScalar(SpinAngular.OneParticleOperator(0, Basics.plus, true), 
                                                       csf, csf, basis.subshells)
        # Add to the existing list
        for  cf in coeffs   push!(coeffs1p, SpinAngular.Coefficient1p(cf.nu, cf.a, cf.b, cf.T / ncsf) )   end
        
        coeffs = SpinAngular.computeCoefficientsScalar(SpinAngular.TwoParticleOperator(0, Basics.plus, true), 
                                                       csf, csf, basis.subshells)
        # Add to the existing lists
        for  cf in coeffs   push!(coeffs2p, SpinAngular.Coefficient2p(cf.nu, cf.a, cf.b, cf.c, cf.d, cf.V / ncsf) )   end
    end 

    # Condense angular coefficients if they refer to the same set of orbital; 
    # include symmetry <ab||cd> == <ba||dc> for symmetric interactions
    coeffs1px = SpinAngular.Coefficient1p[];     coeffs2px = SpinAngular.Coefficient2p[]
    
    hasConsidered = falses( length(coeffs1p) );   T = 0.
    for  (ic, cf) in enumerate(coeffs1p)
        if    hasConsidered[ic]   
        else  nu = cf.nu;   a = cf.a;   b = cf.b
            T = T + cf.T;    hasConsidered[ic] = true
            for   (icx, cfx) in enumerate(coeffs1p)
                if    hasConsidered[icx]   
                elseif  nu == cfx.nu  &&  a == cfx.a  &&  b == cf.b     T = T + cfx.T;    hasConsidered[icx] = true
                end 
            end
            push!(coeffs1px, SpinAngular.Coefficient1p(nu, a, b, T) );  T = 0.
        end 
    end
    
    hasConsidered = falses( length(coeffs2p) );   V = 0.
    for  (ic, cf) in enumerate(coeffs2p)
        if    hasConsidered[ic]   
        else  
            nu = cf.nu;       a = cf.a;   b = cf.b;   c = cf.c;   d = cf.d
            V  = V + cf.V;    hasConsidered[ic] = true
            for   (icx, cfx) in enumerate(coeffs2p)
                if    hasConsidered[icx]   
                elseif  nu == cfx.nu &&    a == cfx.a  &&  b == cfx.b  &&  c == cfx.c &&  d == cfx.d    
                        V  = V + cfx.V;    hasConsidered[icx] = true
                ## elseif  nu == cfx.nu &&  a == cfx.b && b == cf.a &&  c == cfx.d &&  d == cf.c    
                ##         V = V + cfx.V;    hasConsidered[icx] = true
                end 
            end
            push!(coeffs2px, SpinAngular.Coefficient2p(nu, a, b, c, d, V) );   V = 0.
        end 
    end
        
    return( (coeffs1px, coeffs2px) )
end


"""
`SelfConsistent.computeConvergence(aOrbitals::Dict{Subshell, Orbital}, bOrbitals::Dict{Subshell, Orbital}, grid::Radial.Grid)` 
    ... computes a measure for the orbital convergence, i.e. the convergence of all generated orbitals.
        Such a measure is given by the modulus of all products of overlap integrals 
        vonv = | Prod_a   <a|a>  =!= 1, if all orbitals are the same. The procedure assumes that all orbitals
        in aOrbitals and bOrbitals are normalized. The procedure terminates with an error message if any
        orbital in aOrbitals is not found also in aOrbitals. A value conv::Float64 is returned.
"""
function computeConvergence(aOrbitals::Dict{Subshell, Orbital}, bOrbitals::Dict{Subshell, Orbital}, grid::Radial.Grid)
    conv = 1.
    
    for (asubsh, aorb)  in  aOrbitals
        borb = bOrbitals[asubsh]   # Deal with situations that no orbital is found
        conv = conv * RadialIntegrals.overlap(aorb, borb, grid)
    end
    
    return( abs(conv) )
end


"""
`SelfConsistent.computeConvergence(aVectors::Dict{Subshell, Vector{Float64}}, bVectors::Dict{Subshell, Vector{Float64}},
                                   matrixB::Array{Float64,2})` 
    ... computes a measure for the orbital convergence due to the B-spline expansion coefficients, i.e. the convergence 
        of all generated orbitals in terms of their expansion coefficients. Such a measure is given by the modulus 
        of all products of overlap integrals 
        
            conv = | Prod_a   aVector[a] * B * bVector[a] | =!= 1, if all orbitals are the same, and where 
        
        matrixB contains the overlap matrix of the primitives. Formally, this should be always equivalent to
        convergence criterium that is directly based on the orbital representation and, therefore, enables an
        independent test of the generation of orbitals. The procedure assumes that all coefficient vectors 
        are properly normalized. The procedure terminates with an error message if any vector for an orbital 
        in aVectors is not found also in bVectors. A value conv::Float64 is returned.
"""
function computeConvergence(aVectors::Dict{Subshell, Vector{Float64}}, bVectors::Dict{Subshell, Vector{Float64}},
                            matrixB::Array{Float64,2})
    conv = 1.
    
    for (asubsh, aVector)  in  aVectors
        bVector = bVectors[asubsh]   # Deal with situations that no orbital is found
        conv    = conv * (transpose(aVector) * matrixB * bVector)
    end
    
    return( abs(conv) )
end


"""
`SelfConsistent.computeDirectExchangeV(subshell::Subshell, coeffs2p::Array{SpinAngular.Coefficient2p,1}, 
                                       primitives::Bsplines.Primitives, orbitals::Dict{Subshell, Orbital})`
    ... computes the direct and exchange contributions to the one-electron Hamiltonian matrix of the given subshelll.
        These contributions and their position in the Hamiltonian matrix can be derived from the position of subshell
        in the individual coefficient. The coefficient need to be "doubled" if the interaction refer to the same
        subshell due to the variation of the functional. Likely, moreover, the coefficient need to be multiplied
        by sqrt(2j + 1) where j refers to subshell (??)
        A (nsL+nsS) x (nsL+nsS) matrixV::Array{Float64,2} is returned.
"""
function computeDirectExchangeV(subshell::Subshell, coeffs2p::Array{SpinAngular.Coefficient2p,1}, 
                                primitives::Bsplines.Primitives, orbitals::Dict{Subshell, Orbital})
    nsL     = primitives.grid.nsL;        nsS = primitives.grid.nsS;    grid = primitives.grid
    matrixV = zeros( nsL+nsS, nsL+nsS )
    
    for  cf  in  coeffs2p
        ## @show cf
        if      subshell == cf.a  &&  subshell == cf.b  &&  subshell == cf.c  &&  subshell == cf.d  
            if     iseven( Basics.subshell_l(subshell) + Basics.subshell_l(cf.b) + cf.nu )
                matrixV = matrixV + 2 * cf.V *
                          InteractionStrength.XL_Coulomb(cf.nu, subshell, orbitals[cf.b], subshell, orbitals[cf.d], primitives)
            else   continue  # Breit coefficient
            end
        elseif  subshell == cf.a  &&  subshell == cf.c  &&  cf.b == cf.d                            
            if     iseven( Basics.subshell_l(subshell) + Basics.subshell_l(cf.b) + cf.nu )
                matrixV = matrixV + cf.V *
                          InteractionStrength.XL_Coulomb(cf.nu, subshell, orbitals[cf.b], subshell, orbitals[cf.d], primitives)
             else   continue  # Breit coefficient
            end
        elseif  subshell == cf.a  &&  subshell == cf.d  &&  cf.b == cf.c                            
            if     iseven( Basics.subshell_l(subshell) + Basics.subshell_l(cf.b) + cf.nu ) 
                matrixV = matrixV + cf.V *
                          InteractionStrength.XL_Coulomb(cf.nu, subshell, orbitals[cf.b], orbitals[cf.c], subshell, primitives)
            else   continue  # Breit coefficient
            end
        elseif  subshell == cf.b  &&  subshell == cf.d  &&  cf.a == cf.c                            
            if     iseven( Basics.subshell_l(subshell) + Basics.subshell_l(cf.a) + cf.nu )
                matrixV = matrixV + cf.V *
                          InteractionStrength.XL_Coulomb(cf.nu, subshell, orbitals[cf.a], subshell, orbitals[cf.c], primitives)
            else   continue  # Breit coefficient
            end
        elseif  subshell == cf.b  &&  subshell == cf.c  &&  cf.a == cf.d                            
            if     iseven( Basics.subshell_l(subshell) + Basics.subshell_l(cf.a) + cf.nu ) 
                matrixV = matrixV + cf.V *
                          InteractionStrength.XL_Coulomb(cf.nu, subshell, orbitals[cf.a], orbitals[cf.d], subshell, primitives)
            else   continue  # Breit coefficient
            end
        elseif  subshell != cf.a  &&  subshell != cf.b  &&  subshell != cf.c  &&  subshell != cf.d  continue
        else    error("No proper interpretation of coeff = $cf")
        end
    end
    
    return( matrixV )
end


"""
`SelfConsistent.computeFunctional(coeffs1p::Array{SpinAngular.Coefficient1p,1}, coeffs2p::Array{SpinAngular.Coefficient2p,1}, 
                                  orbitals::Dict{Subshell, Orbital}, grid::Radial.Grid, potential::Radial.Potential)` 
    ... computes the value of the energy functional upon which this MCDHF optimization is based on.
        The procedure assumes that this functional is simply given by the form 
        
             EF = Sum_t  coeff1p(a_t,b_t) * I(a_t,b_t)  +   
                  Sum_t  coeff2p(k_t; a_t,b_t,c_t,d_t) * <a_t,b_t | R^(k_t) | c_t,d_t>
                  
        where  <a_t,b_t | R^(k_t) | c_t,d_t>  denotes the rank_k two-electron interaction strengths of the orbitals
        a, ,b , c, and d. An energy::Float64 is returned. 
"""
function computeFunctional(coeffs1p::Array{SpinAngular.Coefficient1p,1}, coeffs2p::Array{SpinAngular.Coefficient2p,1}, 
                           orbitals::Dict{Subshell, Orbital}, grid::Radial.Grid, potential::Radial.Potential)
    energy = 0.

    # Collect one-electron contributions
    for  cf  in  coeffs1p
        jj     = Basics.subshell_2j(cf.a)
        energy = energy + cf.T * sqrt( jj + 1) * RadialIntegrals.GrantIab(orbitals[cf.a], orbitals[cf.b], grid, potential)
    end
    
    # Collect two-electron contributions
    for  cf  in  coeffs2p
        energy = energy + cf.V * InteractionStrength.XL_Coulomb(cf.nu, orbitals[cf.a], orbitals[cf.b], 
                                                                       orbitals[cf.c], orbitals[cf.d], grid)
    end 
    
    return( energy )
end


"""
`SelfConsistent.initializeBasis(configs::Array{Configuration,1}, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives, 
                                settings::AsfSettings; levelSymmetries::Array{LevelSymmetry,1}=LevelSymmetry[], printout::Bool=false)` 
    ... Initialized a many-electron basis from the given list of configurations, the nuclear model as well as ASF settings.
        It assumes that a proper set of primitives::Primitives has been initialized before. The initial set of orbitals in this
        basis is determined by the settings::AsfSettings.  A basis::Basis is returned.
"""
function initializeBasis(configs::Array{Configuration,1}, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives, 
                         settings::AsfSettings; levelSymmetries::Array{LevelSymmetry,1}=LevelSymmetry[], printout::Bool=true)
    NoElectrons = configs[1].NoElectrons;   subshells = Subshell[];   coreSubshells = Subshell[];     csfs = CsfR[] 
    orbitals    = Dict{Subshell, Orbital}()
    
    # Perform some simple tests: Number of electrons must be equal in all configurations
    for  conf in configs   if  conf.NoElectrons != NoElectrons    error("stop a")   end     end
    
    # Generate a full set of relativistic CSF from the given configurations and collect the associated level symmetries
    relconfList = ConfigurationR[]
    for  conf in configs
        wa = Basics.generateConfigurations(Basics.RelativisticConfigurations(), conf)
        append!( relconfList, wa)
    end
    if  printout    for  i = 1:length(relconfList)    println(">>> include ", relconfList[i])    end   end
    subshells = Basics.generateSubshellList(relconfList)
    Defaults.setDefaults("relativistic subshell list", subshells; printout=printout)

    # Generate the relativistic CSF's for the given subshell list
    csfList = CsfR[]
    for  relconf in relconfList
        newCsfs = Basics.generateCsfRs(relconf, subshells)
        append!( csfList, newCsfs)
    end
    
    # Select CSF with requested symmetry if needed
    if  length(levelSymmetries) == 0
        csfs = csfList          # Take all relativistic CSF into account
    else
        for  csf in csfList
            if  LevelSymmetry(csf.J, csf.parity)  in  levelSymmetries   push!(csfs, csf)    end
        end
    end

    # Determine the number of electrons and the list of coreSubshells
    for  (k,sh)  in  enumerate(subshells)
        mocc = Basics.subshell_2j(sh) + 1;    is_filled = true
        for  csf in csfList
            if  csf.occupation[k] != mocc     is_filled = false;           break   end
        end
        if   is_filled    push!( coreSubshells, subshells[k])      else    break   end
    end
        
    # Initialize the orbitals
    if  typeof(settings.startScfFrom) == StartFromHydrogenic
        if  printout   println("> Start SCF process with hydrogenic orbitals.")   end
        # Generate start orbitals for the SCF field by using B-splines
        orbitals  = Bsplines.generateOrbitalsHydrogenic(subshells, nuclearModel, primitives; printout=printout)
    elseif  typeof(settings.startScfFrom) == StartFromPrevious
        if  printout   println("> Start SCF process from given list of orbitals.energy")    end
        # Taking starting orbitals for the given dictionary; non-relativistic orbitals with a proper nuclear charge
        # are adapted if no orbital is found
        orbitals = Dict{Subshell, Orbital}()
        for  sh in subshells
            if  haskey(settings.startScfFrom.orbitals, sh)  
                orbitals[sh] = settings.startScfFrom.orbitals[sh]
            else
                println("Start orbitals do not contain an Orbital for subshell $sh ")
                orb          = HydrogenicIon.radialOrbital(subsh, nuclearModel.Z, grid)
                orb          = Orbital(orb.subshell, orb.isBound, true, orb.energy, orb.P, orb.Q, orb.Pprime, orb.Qprime, Radial.Grid())
                orbitals[sh] = orb 
            end
        end
    else  error("stop b")
    end
    
    basis = Basis(true, NoElectrons, subshells, csfs, coreSubshells, orbitals)
    return( basis )
end


"""
`SelfConsistent.performSCF(configs::Array{Configuration,1}, nm::Nuclear.Model, grid::Radial.Grid, 
                           settings::AsfSettings; levelSymmetries::Array{LevelSymmetry,1}=LevelSymmetry[], printout::Bool=false)` 
    ... Performs a SCF computation for the given list of configurations, the nuclear model as well as ASF settings.
        If explicit levelSymmetries are given, only these symmetries are considered. Internally, a proper set of primitives::Primitives 
        is initialized and used in the computations. The generated SCF field is controlled by the settings::AsfSettings.  
        A multiplet::Multiplet is returned.
"""
function performSCF(configs::Array{Configuration,1}, nm::Nuclear.Model, grid::Radial.Grid, 
                    settings::AsfSettings; levelSymmetries::Array{LevelSymmetry,1}=LevelSymmetry[], printout::Bool=true)
    
    # Generate primitives and initialize the many-electron basis
    Defaults.setDefaults("standard grid", grid)
    primitives = Bsplines.generatePrimitives(grid)    
    basis      = SelfConsistent.initializeBasis(configs, nm, primitives, settings; levelSymmetries, printout)
    
    # Solve a self-consistent field for this basis
    if   typeof(settings.scField)  in  [Basics.DFSField, Basics.DFSwCPField, Basics.HSField]
        basis = SelfConsistent.solveMeanFieldBasis(basis, nm, primitives, settings; printout=printout) 
    elseif   settings.scField in [Basics.NuclearField()]  && settings.startScfFrom == StartFromHydrogenic() 
        # Return the basis as already generated.
    elseif   settings.scField in [Basics.ALField()] 
        basis     = SelfConsistent.solveAverageLevelField(basis, nm, primitives, settings; printout=printout)
        multiplet = Hamiltonian.performCI(basis, nm, primitives.grid, settings, printout=true)
    elseif   typeof(settings.scField) == Basics.EOLField 
        multiplet = SelfConsistent.solveOptimizedLevelField(basis, nm, primitives, settings; printout=printout)
        basis     = multiplet.levels[1].basis
    else  error("stop a")
    end    
    
    # Setup and diagonalize the Hamiltonian matrix; assign mixing coefficients
    mp = Hamiltonian.performCI(basis, nm, grid, settings, printout=printout)

    return( mp )
end


"""
`SelfConsistent.performSCF(basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, 
                           settings::AsfSettings; levelSymmetries::Array{LevelSymmetry,1}=LevelSymmetry[], printout::Bool=false)` 
    ... Performs a SCF computation for the given list of configurations, the nuclear model as well as ASF settings.
        If explicit levelSymmetries are given, only these symmetries are considered. Internally, a proper set of primitives::Primitives 
        is initialized and used in the computations. The generated SCF field is controlled by the settings::AsfSettings.  
        A multiplet::Multiplet is returned.
"""
function performSCF(basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, 
                    settings::AsfSettings; levelSymmetries::Array{LevelSymmetry,1}=LevelSymmetry[], printout::Bool=false)
    
    # Generate primitives
    primitives = Bsplines.generatePrimitives(grid)    
    
    # Solve a self-consistent field for this basis
    if   typeof(settings.scField)  in  [Basics.DFSField, Basics.DFSwCPField, Basics.HSField]
        basis = SelfConsistent.solveMeanFieldBasis(basis, nm, primitives, settings; printout=printout) 
    elseif   settings.scField in [Basics.NuclearField()]  && settings.startScfFrom == StartFromHydrogenic() 
        # Return the basis as already generated.
    elseif   settings.scField in [Basics.ALField()] 
        basis     = SelfConsistent.solveAverageLevelField(basis, nm, primitives, settings; printout=printout)
        multiplet = Hamiltonian.performCI(basis, nm, primitives.grid, settings, printout=true)
    elseif   typeof(settings.scField) == Basics.EOLField 
        multiplet = SelfConsistent.solveOptimizedLevelField(basis, nm, primitives, settings; printout=printout)
        basis     = multiplet.levels[1].basis
    else  error("stop a")
    end    
    
    # Setup and diagonalize the Hamiltonian matrix; assign mixing coefficients
    mp = Hamiltonian.performCI(basis, nm, settings, grid, printout=printout)

    return( mp )
end


"""
`SelfConsistent.rotateOrbitals(subshellList::Array{Subshell,1}, orbitals::Dict{Subshell, Orbital}, grid::Radial.Grid, 
                               settings::AsfSettings)` 
    ... rotates pairwise the orbitals to enhance the covergence of the SCF procedures; it follows that 
        rotation schemes of C.F. Fischer ... . The pairwise rotation of the orbitals of the same symmetry 
        is done from the inner-to-outer subshells. A rotation is not necessary if two subshells are filled
        of if both are fixed. For each rotation, a rotation parameter beta is determined.
        A pair (newOrbitals::Dict{Subshell, Orbital}, betaMax::Float64) us returned which contains the rotated
        orbitals as well as the maximum beta parameter that occurred during the pairwise rotation of the orbital set.
"""
function rotateOrbitals(subshellList::Array{Subshell,1}, orbitals::Dict{Subshell, Orbital}, grid::Radial.Grid, 
                        settings::AsfSettings)
    betaMax = 0.;    newOrbitals = orbitals
    
    @warn("Not yet properly implemented.")
    
    return( newOrbitals, betaMax)
end


"""
`SelfConsistent.determineChemicalPotential(orbitals::Dict{Subshell, Orbital}, temp::Float64, radiusWS::Float64,
                                           nm::Nuclear.Model, grid::Radial.Grid)`
    ... determines the chemical potential so that Sum_i f(epsilon_i, mu, temp) = Z.
        The Newton-Raphson methods is used to iterate to the chemical potential; a chemMu::Float64 is returned.

        Note: this general finite-temperature Fermi-Dirac root-finding utility was moved here from module Plasma
              (where it originated as Plasma.determineChemicalPotential), since Plasma.perform(::AverageAtomScheme,
              ...) needs SelfConsistent.solveAverageAtomField below, and solveAverageAtomField itself needs this
              function internally at every SCF iteration; keeping it in Plasma would have made the two modules
              depend on each other circularly. Nothing here is Plasma-scheme-specific.
"""
function determineChemicalPotential(orbitals::Dict{Subshell, Orbital}, temp::Float64, radiusWS::Float64, nm::Nuclear.Model,
                                    grid::Radial.Grid)
    function g(mu::Float64, orbitals::Dict{Subshell, Orbital}, temp::Float64, nm::Nuclear.Model)
        wa = - nm.Z
        for  (k,v)  in orbitals
            occ = Basics.twice( Basics.subshell_j(k)) + 1
            wb  = (v.energy - mu) /temp
            if  wb > 300.   wb = 300.   end
            wa  = wa + occ / (exp(wb) + 1)
        end
        return( wa )
    end
    #
    function gprime(mu::Float64, orbitals::Dict{Subshell, Orbital}, temp::Float64, nm::Nuclear.Model)
        wa = 0.
        for  (k,v)  in orbitals
            occ = Basics.twice( Basics.subshell_j(k)) + 1
            wb  = (v.energy - mu) /temp
            if  wb > 300.   wb = 300.   end
            wc  = exp( wb )
            wa  = wa + occ * wc^2 / temp / (wc+1)^2
        end
        return( wa )
    end
    # Iterate for the chemical potential
    chemMu = -0.1;     newMu = 0.;     nx = 0
    while true
        nx = nx + 1
        newMu = chemMu - g(chemMu, orbitals, temp, nm) / gprime(chemMu, orbitals, temp, nm)
        if  abs(newMu - chemMu) < 1.0e-4  break
        else    chemMu = newMu
        end
    end

    chemMu = chemMu - 0.0011  ## Seems to bring better stability in the SCF computations

    println(">>> Newton-Raphson: $nx)  chemMu = $chemMu  g = $(g(chemMu, orbitals, temp, nm)) ")
    return ( chemMu )
end


"""
`SelfConsistent.solveAverageAtomField(orbitals::Dict{Subshell, Orbital}, nuclearModel::Nuclear.Model, scField::Basics.AbstractScField,
                                      temp::Float64, radiusWS::Float64, primitives::Bsplines.Primitives; printout::Bool=true)`
    ... solves the self-consistent field for a given local average-atom potential as specified by scField
        A (new) set of orbitals::Dict{Subshell, Orbital} is returned.
"""
function solveAverageAtomField(orbitals::Dict{Subshell, Orbital}, nuclearModel::Nuclear.Model, scField::Basics.AbstractScField,
                               temp::Float64, radiusWS::Float64, primitives::Bsplines.Primitives; printout::Bool=true)
    # Determine the chemical potential
    chemMu    = determineChemicalPotential(orbitals, temp, radiusWS, nuclearModel, primitives.grid);        @show chemMu
    
    # Extract the kappa's from orbitals
    kappas = Int64[];     for (k,v)  in  orbitals     push!(kappas, k.kappa)    end;    kappas = unique(kappas);   @show kappas

    ## Defaults.setDefaults("standard grid", primitives.grid; printout=printout)
    # Define the storage for the calculations of matrices
    if  printout    println(">> (Re-) Define a storage array for dealing with single-electron TTp B-spline matrices:")    end
    storage  = Dict{String,Array{Float64,2}}()
    
    # Set-up the overlap matrix; compute or fetch the diagonal 'overlap' blocks
    nsL = primitives.grid.nsL;        nsS = primitives.grid.nsS;    grid = primitives.grid
    wb  = zeros( nsL+nsS, nsL+nsS )
    wb[1:nsL,1:nsL]                 = Bsplines.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
    wb[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = Bsplines.generateTTpMatrix!("SS-overlap", 0, primitives, storage)
    
    # Determine the symmetry block of this basis and define storage for the kappa blocks and orbitals from the last iteration
    bsplineBlock = Dict{Int64,Basics.Eigen}();   previousOrbitals = deepcopy(orbitals)
    for  kappa  in  kappas           bsplineBlock[kappa]  = Basics.Eigen( zeros(2), [zeros(2), zeros(2)])   end
    # Determine te nuclear potential once at the beginning
    nuclearPotential  = Nuclear.nuclearPotential(nuclearModel, grid)
            
    # Start the SCF procedure for all symmetries
    isNotSCF = true;   NoIteration = 0;   accuracyScf = 0.
    while  isNotSCF
        NoIteration = NoIteration + 1;   go_on = false 
        if  NoIteration >  32
                println(">> Maximum number of SCF iterations = 32 is reached at accuracy " * 
                        @sprintf("%.4e", accuracyScf) * " ... computations proceed.")
            break
        end
        if  printout    println("\nIteration $NoIteration for symmetries ... ")    end
        #
        for kappa in kappas
            # (1) Re-compute the local potential
            wp  = Basics.computePotential(scField, grid, previousOrbitals, chemMu, temp)
            pot = Basics.add(nuclearPotential, wp)
            
            # (2) Set-up the diagonal part of the Hamiltonian matrix
            wa = Bsplines.setupLocalMatrix(kappa, primitives, pot, storage)
            # (3) Solve the generalized eigenvalue problem
            wc = Basics.diagonalize(GeneralizedEigenvaluesWithLinearAlgebra(), wa, wb)
            
            # (4) Analyse and print information about the convergence of the symmetry blocks and the occupied orbitals
            wcBlock = Basics.analyzeConvergence(bsplineBlock[kappa], wc)
            if  wcBlock > 1.0e-6   go_on = true   end     ## accuracyScf
            for  (k,v)  in  orbitals
                if      k.kappa == kappa
                    newOrbital = Bsplines.generateOrbitalFromPrimitives(k, wc, primitives)
                    wcOrbital  = Basics.analyzeConvergence(previousOrbitals[k], newOrbital)
                    if  wcOrbital > 1.0e-6   accuracyScf = wcOrbital;   go_on = true   end     ## accuracyScf
                        sa = "  $k::  en [a.u.] = " * @sprintf("%.7e", newOrbital.energy) * ";   self-cons'cy = "  
                        sa = sa * @sprintf("%.4e", wcOrbital)   * "  ["
                        sa = sa * @sprintf("%.4e", wcBlock)             * " for sym-block kappa = $kappa]"
                        if  printout    println(sa)    end
                    ## println("  $sh  en [a.u.] = $(newOrbital.energy)   self-consistency = $(wcOrbital), $(wcBlock) [kappa=$kappa] ") 
                    previousOrbitals[k] = newOrbital
                end
            end
            # (5) Re-define the bsplineBlock
            bsplineBlock[kappa] = wc
        end
        chemMu              = determineChemicalPotential(previousOrbitals, temp, radiusWS, nuclearModel, primitives.grid)
        if  go_on   nothing   else   break   end
    end
    
    analyzedOrbitals = Basics.analyze(previousOrbitals, printout=true)    
    return( analyzedOrbitals )
end


"""
`SelfConsistent.solveAverageLevelField(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives, 
                                       settings::AsfSettings; printout::Bool=true)` 
    ... solves the self-consistent field for the given orbitals (from basis), the nuclear model as well as
        for the average-level (AL) functional. In addition, the settings::AsfSettings are taken into account.
        A (new) basis::Basis is returned.
"""
function solveAverageLevelField(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives, 
                                settings::AsfSettings; printout::Bool=true)
    nsL    = primitives.grid.nsL;    nsS = primitives.grid.nsS;    grid = primitives.grid
    rotate = false
    
    # (1) Initialize storage and important arrays; determine nuclear potential and mean occupation once at the beginning   
    if  printout    println(">> (Re-) Define a storage array for dealing with single-electron TTp B-spline matrices:")    end
    storage = Dict{String,Array{Float64,2}}()
    matrixB = zeros( nsL+nsS, nsL+nsS )
    matrixB[1:nsL,1:nsL]                 = Bsplines.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
    matrixB[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = Bsplines.generateTTpMatrix!("SS-overlap", 0, primitives, storage)

    nucPot   = Nuclear.nuclearPotential(nuclearModel, primitives.grid)
    meanOcc  = Basics.extractMeanOccupation(basis)
    orbitals = deepcopy( basis.orbitals )
    bVectors = Dict{Subshell, Vector{Float64}}()
    v = zeros(nsL+nsS);     for  (k, orb)  in  orbitals   bVectors[k] = v   end
    
    # (2) Generate angular coefficients for all diagonal matrix elements including the division by nCSF 
    # to obtain mean energy functional
    (coeffs1p, coeffs2p) = SelfConsistent.computeAngularCoefficients(settings.scField, basis::Basis)
    for  cf in  coeffs1p   wx = Basics.twice(Basics.subshell_j(cf.a)) + 1;   @show cf, cf.T * sqrt(wx)   end
    
    # (3) Iterate the orbital set of the given basis until convergence is achieved; initialize bVectors in order to keep 
    # the B-spline coefficients of each orbital for the orbital projections
 
    for iter = 1:settings.maxIterationsScf
        println("\n> SCF interation $(iter): ")
        newOrbitals = Dict{Subshell, Orbital}();   newbVectors = Dict{Subshell, Vector{Float64}}()
        
        # (4) Rotate orbitals pairwise if required; such rotations are not needed if both subshells are filled
        #     determine a maximum rotation parameter
        if rotate
            orbitals = SelfConsistent.rotateOrbitals(basis.subshells, orbitals, grid, settings)
        end 
        
        # (5) Cycle through all subshells from the inner-to-outer subshells; 
        for subshell  in  basis.subshells
            occ      = meanOcc[subshell];  occm1 = 1.0 / occ   
            orb      = orbitals[subshell] 
            print(">> Refine $subshell orbital with mean occ = $occ ... ")
            
            # (6) Set-up the Hamiltonian matrix for this shell, including the Dirac Hamiltonian + direct + exchange potentials
            matrix  = Bsplines.setupLocalMatrix(subshell.kappa, primitives, nucPot, storage)
            matrixV = SelfConsistent.computeDirectExchangeV(subshell, coeffs2p, primitives, orbitals)
            matrix  = matrix + occm1 * matrixV
            
            # (7) Apply projections of matrix upon the other orbitals of the same symmetry (accounts for Langrange multipliers)
            ## matrix = Hamiltonian.projectHamiltonian(subshell, matrix, matrixB, bVectors)
            
            # (8) Diagonalize the Hamiltonian matrix and refine the orbital and bVectors
            wc = Basics.diagonalize(GeneralizedEigenvaluesWithLinearAlgebra(), matrix, matrixB)
            newOrb                = Bsplines.generateOrbitalFromPrimitives(subshell, wc, primitives)
            newOrbitals[subshell] = newOrb
            newbVectors[subshell] = Bsplines.extractVectorFromPrimitives(subshell,   wc, primitives)
            
            # (9) Report about the new orbital
            ovlap = abs( RadialIntegrals.overlap(orb, newOrb, grid) )
            println("     overlap = $ovlap   acc = $(1.0 - ovlap)  ... ")
        end
        
        # (10) Compute energy functional and analyze the convergence of orbitals
        eFunctional = SelfConsistent.computeFunctional(coeffs1p, coeffs2p, newOrbitals, grid, nucPot)
        orbitalConv = SelfConsistent.computeConvergence(orbitals, newOrbitals, grid)
        bVectorConv = SelfConsistent.computeConvergence(bVectors, newbVectors, matrixB)
        
        println(">> Total energy = $(eFunctional*1)   orbital-conv = $orbitalConv   orbital-acc = $(1.0 - orbitalConv)   " *
                "B-conv = $bVectorConv")
        
        if  abs(1.0 - orbitalConv) < settings.accuracyScf    break   end
        orbitals = deepcopy(newOrbitals);    bVectors = newbVectors
    end
    
    newBasis = Basis(true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, orbitals)
    return( newBasis )
end

    
function normX(sa::String, matrix::Array{Float64,2})
    absD = absND = 0.
    for  i = 1:size(matrix,1)                    absD  = absD  + abs(matrix[i,i])   
        for  j = 1:size(matrix,2)   if i != j    absND = absND + abs(matrix[i,j])  end   end
    end 
    @show  sa, size(matrix), absD, absND
    return(nothing)
end


"""
`SelfConsistent.solveMeanFieldBasis(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives, 
                                    settings::AsfSettings; printout::Bool=true)` 
    ... solves the self-consistent field for the given orbitals (from basis), the nuclear model as well as
        the (local) mean-field potential as specified by the settings::AsfSettings. A (new) basis::Basis is returned.
"""
function solveMeanFieldBasis(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives, 
                             settings::AsfSettings; printout::Bool=true)
    ## Defaults.setDefaults("standard grid", primitives.grid; printout=printout)
    # Define the storage for the calculations of matrices
    if  printout    println(">> (Re-) Define a storage array for dealing with single-electron TTp B-spline matrices:")    end
    storage  = Dict{String,Array{Float64,2}}()
    
    # Set-up the overlap matrix; compute or fetch the diagonal 'overlap' blocks
    nsL = primitives.grid.nsL;    nsS = primitives.grid.nsS
    wb  = zeros( nsL+nsS, nsL+nsS )
    wb[1:nsL,1:nsL]                 = Bsplines.generateTTpMatrix!("LL-overlap", 0, primitives, storage)
    wb[nsL+1:nsL+nsS,nsL+1:nsL+nsS] = Bsplines.generateTTpMatrix!("SS-overlap", 0, primitives, storage)
    
    # Determine te nuclear potential once at the beginning
    nuclearPotential  = Nuclear.nuclearPotential(nuclearModel, primitives.grid)
    
    # Determine the symmetry block of this basis and define storage for the kappa blocks and orbitals from the last iteration
    kappas   = Int64[];   for sh in basis.subshells  push!(kappas, sh.kappa)   end;   kappas = unique(kappas)

    bsplineBlock = Dict{Int64,Basics.Eigen}();   previousOrbitals = deepcopy(basis.orbitals)
    for  kappa  in  kappas           bsplineBlock[kappa]  = Basics.Eigen( zeros(2), [zeros(2), zeros(2)])   end
    
    # Start the SCF procedure for all symmetries
    isNotSCF = true;   NoIteration = 0;   accuracyScf = 0.
    while  isNotSCF
        NoIteration = NoIteration + 1;   go_on = false 
        if  NoIteration >  settings.maxIterationsScf
                println(">> Maximum number of SCF iterations = $(settings.maxIterationsScf) is reached at accuracy " * 
                        @sprintf("%.4e", accuracyScf) * " ... computations proceed.")
            break
        end
        if  printout    println("\nIteration $NoIteration for symmetries ... ")    end
        #
        for kappa in kappas
            # (1) First re-define an (arbitrary) 'level' that represents the mean occupation for the local potential
            wBasis = Basis(true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, previousOrbitals)
            NoCsf  = length(wBasis.csfs)
            wmc    = zeros( NoCsf );   wN = 0.
            for i = 1:NoCsf   wmc[i] = Basics.twice(wBasis.csfs[i].J) + 1.0;   wN = wN + abs(wmc[i])^2    end
            for i = 1:NoCsf   wmc[i] = wmc[i] / sqrt(wN)   end
            wLevel = Level( AngularJ64(0), AngularM64(0), Basics.plus, 0, -1., 0., true, wBasis, wmc)
            # (2) Re-compute the local potential
            wp  = Basics.computePotential(settings.scField, primitives.grid, wLevel)
            pot = Basics.add(nuclearPotential, wp)
            # (3) Set-up the diagonal part of the Hamiltonian matrix
            wa = Bsplines.setupLocalMatrix(kappa, primitives, pot, storage)
            # (4) Solve the generalized eigenvalue problem
            wc = Basics.diagonalize(GeneralizedEigenvaluesWithLinearAlgebra(), wa, wb)
            # (5) Analyse and print information about the convergence of the symmetry blocks and the occupied orbitals
            wcBlock = Basics.analyzeConvergence(bsplineBlock[kappa], wc)
            if  wcBlock > 1.000 * settings.accuracyScf   go_on = true   end
            for  sh in basis.subshells
                if      sh in settings.frozenSubshells   ## do nothing
                elseif  sh.kappa == kappa
                    newOrbital = Bsplines.generateOrbitalFromPrimitives(sh, wc, primitives)
                    wcOrbital  = Basics.analyzeConvergence(previousOrbitals[sh], newOrbital)
                    if  wcOrbital > settings.accuracyScf   accuracyScf = wcOrbital;   go_on = true   end
                        sa = "  $sh::  en [a.u.] = " * @sprintf("%.7e", newOrbital.energy) * ";   self-cons'cy = "  
                        sa = sa * @sprintf("%.4e", wcOrbital)   * "  ["
                        sa = sa * @sprintf("%.4e", wcBlock)             * " for sym-block kappa = $kappa]"
                        if  printout    println(sa)    end
                    previousOrbitals[sh] = newOrbital
                end
            end
            # (6) Re-define the bsplineBlock
            bsplineBlock[kappa] = wc
        end
        if  go_on   nothing   else   break   end
    end
    
    analyzedOrbitals = Basics.analyze(previousOrbitals, printout=true)
    
    newBasis = Basis(true, basis.NoElectrons, basis.subshells, basis.csfs, basis.coreSubshells, analyzedOrbitals)
    return( newBasis )
end
    

"""
`SelfConsistent.solveOptimizedLevelField(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives, 
                                         settings::AsfSettings; printout::Bool=true)` 
    ... solves the self-consistent field for the given orbitals (from basis), the nuclear model as well as
        for the average-level (AL) functional. In addition, the settings::AsfSettings are taken into account.
        A (new) multiplet::Multiplet is returned.
"""
function solveOptimizedLevelField(basis::Basis, nuclearModel::Nuclear.Model, primitives::Bsplines.Primitives, 
                                  settings::AsfSettings; printout::Bool=true)
    error("Not yet implemented; this was never done so far.")
    
    multiplet = Multiplet()
    return( multiplet )
end

end # module
