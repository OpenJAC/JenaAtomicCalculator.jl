

# The finite-temperature AVERAGE-ATOM path. Not reached by SelfConsistent.performSCF -- it is driven
# from Plasma.AverageAtomScheme (module-Plasma-inc-average-atom.jl), and kept here beside the other
# SCF drivers because it iterates orbitals in the same way.

"""
`SelfConsistent.determineChemicalPotential(orbitals::Dict{Subshell, Orbital}, temp::Float64, radiusWS::Float64,
                                           nm::Nuclear.Model, grid::Radial.Grid)`
    ... determines the chemical potential so that Sum_i f(epsilon_i, mu, temp) = Z.
        The Newton-Raphson methods is used to iterate to the chemical potential; a chemMu::Float64 is returned.

        Note: this general finite-temperature Fermi-Dirac root-finding utility was moved here from module Plasma
              (where it originated as `determineChemicalPotential`), since Plasma.perform(::AverageAtomScheme,
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
    chemMu    = determineChemicalPotential(orbitals, temp, radiusWS, nuclearModel, primitives.grid);
    
    # Extract the kappa's from orbitals
    kappas = Int64[];     for (k,v)  in  orbitals     push!(kappas, k.kappa)    end;    kappas = unique(kappas);

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
                ## Collected as well: in a long run this line scrolls away, and nobody learns afterwards that a
                ## field never converged.  The accuracy is rounded so that repeated identical failures collapse
                ## into one counted entry; see Defaults.warn.
                Defaults.warn(AddWarning(), "SelfConsistent.solveAverageAtomField(): the SCF did NOT converge -- " *
                              "stopped at accuracy " * @sprintf("%.1e", accuracyScf) * " after 32 iterations.")
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
            wc = Bsplines.diagonalizeLocalMatrix(kappa, wa, wb, primitives)
            
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
