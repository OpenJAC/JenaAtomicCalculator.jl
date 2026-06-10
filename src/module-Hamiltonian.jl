
"""
`module  JAC.Hamiltonian`  
	... a submodel of JAC that contains all structs and methods to efficiently set-up and work
	    with (large) Hamiltonian matrices under different conditions.
"""
module Hamiltonian

using  Printf, ..Basics, ..BsplinesN, ..Defaults, ..InteractionStrength, ..ManyElectron, ..Nuclear, ..Radial, 
       ..RadialIntegrals, ..SpinAngular


"""
`Hamiltonian.projectHamiltonian(subshell::Subshell, matrix::Array{Float64,2}, matrixB::Array{Float64,2}, 
                                bVectors::Dict{Subshell, Vector{Float64}})` 
    ... projects the (single-electron DHF) Hamiltonian matrix of subshell with regard to the B-vectors of the same
        symmetry. This projection is necessary to ensure orthogonality among the various orbitals.
        A (nsL+nsS) x (nsL+nsS) matrix::Array{Float64,2} is returned.
"""
function projectHamiltonian(subshell::Subshell, matrix::Array{Float64,2}, matrixB::Array{Float64,2}, 
                            bVectors::Dict{Subshell, Vector{Float64}})
    nn = size(matrix, 1);   matrixP = deepcopy(matrix);           
    
    matrixI = zeros(nn, nn);   for i = 1:nn   matrixI[i,i] = 1.0   end
    
    for  (subsh, bVector)  in  bVectors
        if  subsh != subshell  &&  subshell.kappa == subsh.kappa
            ## println(">>> Project $subshell  upon $subsh  vectors")
            bbtMatrix = bVector * transpose(bVector)
            matrixP   = (matrixI - matrixB * bbtMatrix) * matrixP * (matrixI - bbtMatrix * matrixB)
        end 
    end 
    
    return( matrix )
end


"""
`Hamiltonian.performCI(basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings; printout::Bool=false)
    ... Computes and diagonalizes the Hamiltonian matrix for all CSF in the given basis. It also assigns the
        mixing coefficients to the individual levels. The individual contributions from the Breit or diagonal interaction as well as 
        from QED to this matrix are controlled by the settings. A  multiplet::Multiplet  is returned.
"""
function performCI(basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings; printout::Bool=false)
    
    # First determine the number of CSF in each J^P symmetry block
    symmetries = Dict{LevelSymmetry,Int64}()
    for  csf in basis.csfs
        sym = LevelSymmetry(csf.J, csf.parity)
        if     haskey(symmetries, sym)    symmetries[sym] = symmetries[sym] + 1
        else                              symmetries[sym] = 1
        end
    end

    # Test the total number of CSF
    NoCsf = 0;   for (sym,v) in symmetries   NoCsf = NoCsf + v   end
    if  NoCsf != length(basis.csfs)   error("stop b; NoCsf = $NoCsf ")   end

    # Calculate for each symmetry block the corresponding CI matrix, diagonalize it and append a Multiplet for this block
    multiplets = Multiplet[]
    for  (sym,v) in  symmetries
        # Skip the symmetry block if it not selected
        if  !Basics.selectSymmetry(sym, settings.levelSelectionCI)     continue    end
        matrix = Hamiltonian.setupMatrix(sym, basis, nm, grid, settings; printout=printout)
        eigen  = Basics.diagonalize("matrix: LinearAlgebra", matrix)
        
        # Reassign state vectors to levels
        levels = Level[]
        for  ev = 1:length(eigen.values)
            # Construct the eigenvector with regard to the given basis (not w.r.t the symmetry block)
            evSym = eigen.vectors[ev];    vector = zeros( length(basis.csfs) );   ns = 0
            for  r = 1:length(basis.csfs) 
                if LevelSymmetry(basis.csfs[r].J, basis.csfs[r].parity) == sym    ns = ns + 1;   vector[r] = evSym[ns]   end
            end
            newlevel = Level( sym.J, AngularM64(sym.J.num//sym.J.den), sym.parity, 0, eigen.values[ev], 0., true, basis, vector ) 
            push!( levels, newlevel)
        end
        wa = Multiplet(string(sym) * "+", levels)
        push!( multiplets, wa)
    end
    
    # Merge all multiplets into a single one
    mp = Basics.merge(multiplets)
    mp = Basics.sortByEnergy(mp)
    
    # Determine the level list to be printed out
    levelNos = Int64[]
    for (ilev, level) in  enumerate(mp.levels)
        if  Basics.selectLevel(level, settings.levelSelectionCI)    push!(levelNos, ilev)    end
    end
    
    # Display all level energies and energy splittings
    if  printout
        Basics.tabulate(stdout, mp, levelNos)
    end
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary     
        Basics.tabulate(iostream, mp, levelNos)
    end

    return( mp )
end


"""
`Hamiltonian.performCIwithFrozenOrbitals(configs::Array{Configuration,1}, frozenOrbitals::Dict{Subshell, Orbital}, 
                                         nuclearModel::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings; printout::Bool=true)` 
    ... to generate from the given frozen orbitals a multiplet of single-CSF levels by just using the diagonal 
        part of the Hamiltonian matrix; a multiplet::Multiplet is returned.  
"""
function performCIwithFrozenOrbitals(configs::Array{Configuration,1}, frozenOrbitals::Dict{Subshell, Orbital}, 
                                     nuclearModel::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings; printout::Bool=true)
    if  printout    println("\n... in Hamiltonian.perform...: a mutiplet from orbitals, no CI, CSF diagonal'] ...")    end
    
    if  settings.eeInteractionCI == DiagonalCoulomb()   error("Not yet implemented.")   end
    
    # Generate a list of relativistic configurations and determine an ordered list of subshells for these configurations
    relconfList = ConfigurationR[]
    for  conf in configs
        wa = Basics.generateConfigurations(Basics.RelativisticConfigurations(), conf)
        append!( relconfList, wa)
    end
    if  printout    for  i = 1:length(relconfList)    println(">> include ", relconfList[i])    end   end
    subshellList = Basics.generateSubshellList(relconfList)
    Defaults.setDefaults("relativistic subshell list", subshellList; printout=printout)

    # Generate the relativistic CSF's for the given subshell list
    csfList = CsfR[]
    for  relconf in relconfList
        newCsfs = Basics.generateCsfRs(relconf, subshellList)
        append!( csfList, newCsfs)
    end

    # Determine the number of electrons and the list of coreSubshells
    NoElectrons      = sum( csfList[1].occupation )
    coreSubshellList = Subshell[]
    for  k in 1:length(subshellList)
        mocc = Basics.subshell_2j(subshellList[k]) + 1;    is_filled = true
        for  csf in csfList
            if  csf.occupation[k] != mocc    is_filled = false;    break   end
        end
        if   is_filled    push!( coreSubshellList, subshellList[k])    end
    end
    
    # Set-up a basis for calculating the Hamiltonian matrix
    basis = Basis(true, NoElectrons, subshellList, csfList, coreSubshellList, frozenOrbitals)

    # Calculate the diagonal part of the Hamiltonian matrix and define a multiplet for these approximate ASf;
    # Generate first an effective nuclear charge Z(r) on the given grid
    potential = Nuclear.nuclearPotential(nuclearModel, grid)
    
    keep = true
    InteractionStrength.XL_Coulomb_reset_storage(keep)
    #
    n = length(csfList);    matrix = zeros(Float64, n, n)
    for  r = 1:n
        subshellList = basis.subshells
        opa  = SpinAngular.OneParticleOperator(0, plus, true)
        waG1 = SpinAngular.computeCoefficients(opa, basis.csfs[r], basis.csfs[r], subshellList) 
        opa  = SpinAngular.TwoParticleOperator(0, plus, true)
        waG2 = SpinAngular.computeCoefficients(opa, basis.csfs[r], basis.csfs[r], subshellList)
        wa   = [waG1, waG2]
        #
        me = 0.
        for  coeff in wa[1]
            jj = Basics.subshell_2j(basis.orbitals[coeff.a].subshell)
            me = me + coeff.T * sqrt( jj + 1) * RadialIntegrals.GrantIab(basis.orbitals[coeff.a], basis.orbitals[coeff.b], grid, potential)
        end

        for  coeff in wa[2]
            if  typeof(settings.eeInteractionCI) in [CoulombInteraction, CoulombBreit, CoulombGaunt]
                me = me + coeff.V * InteractionStrength.XL_Coulomb(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                                basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid, keep=keep)  end
            if  typeof(settings.eeInteractionCI) in [BreitInteraction, CoulombBreit]
                me = me + coeff.V * InteractionStrength.XL_Breit(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                           basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid,
                                                                           settings.eeInteractionCI)                                                end
        end
        matrix[r,r] = me
    end 
    #
    eigen  = Basics.diagonalize("matrix: LinearAlgebra", matrix)
    levels = Level[]
    for  ev = 1:length(eigen.values)
        vector = eigen.vectors[ev];   ns = 0
        for  r = 1:length(vector)   
            if      (vector[r])^2 > 0.99    ns = r;   break    
            elseif  (vector[r])^2 > 0.01    error("stop a; unexpected mixing coefficieint; vector = $vector ")    
            end
        end
        J = csfList[ns].J
        newlevel = Level( J, AngularM64(J.num//J.den), csfList[ns].parity, 0, eigen.values[ev], 0., true, basis, vector ) 
        push!( levels, newlevel)
    end
    mp = Multiplet("noName", levels)
    mp = Basics.sortByEnergy(mp)
    
    levelNos = Int64[]
    for (ilev, level) in  enumerate(mp.levels)       push!(levelNos, ilev)   end
    
    # Display all level energies and energy splittings
    if  printout
        Basics.tabulate(stdout, mp, levelNos)
    end
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary  &&  printout    
        Basics.tabulate(iostream, mp, levelNos)
    end

    return( mp )
end



"""
`Hamiltonian.setupMatrix(sym::LevelSymmetry, basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings; printout::Bool=false)
    ... Set-up (computes) the Hamiltonian matrix for all CSF with symmetry sym in the given basis. The individual contributions
        from the Breit or diagonal interaction as well as from QED to this matrix are controlled by the settings.
        A  matrix::Arrays{Float64,2}  is returned.
"""
function setupMatrix(sym::LevelSymmetry, basis::Basis, nm::Nuclear.Model, grid::Radial.Grid, settings::AsfSettings; printout::Bool=false)

    # Determine the dimension of the CI matrix and the indices of the CSF with J^P symmetry in the basis
    idx_csf = Int64[]
    for  idx = 1:length(basis.csfs)
        if  basis.csfs[idx].J ==sym.J   &&   basis.csfs[idx].parity == sym.parity    push!(idx_csf, idx)    end
    end
    n = length(idx_csf)
    if printout    print("> Compute CI matrix of dimension $n x $n for the symmetry $(string(sym.J))^$(string(sym.parity)) ...")    end

    # Generate an effective nuclear charge Z(r) on the given grid to add QED contributions, if requested
    potential = Nuclear.nuclearPotential(nm, grid)
    if  settings.qedModel in [QedPetersburg(), QedSydney()]    
        meanPot = potential
        ## meanPot = compute("radial potential: Dirac-Fock-Slater", grid, basis)
        ## meanPot = Basics.add(potential, meanPot)   
    end   

    # Compute the Coulomb-(Breit-) interaction matrix
    matrix = zeros(Float64, n, n);      keep = true
    InteractionStrength.XL_Coulomb_reset_storage(keep, printout=false)
    InteractionStrength.XL_Breit_reset_storage(keep, printout=false)
    for  r = 1:n
        for  s = 1:n
            if  settings.eeInteractionCI == DiagonalCoulomb()  &&  r != s    continue    end
            # Calculate the spin-angular coefficients
            subshellList = basis.subshells
            opa  = SpinAngular.OneParticleOperator(0, plus, true)
            waG1 = SpinAngular.computeCoefficients(opa, basis.csfs[idx_csf[r]], basis.csfs[idx_csf[s]], subshellList) 
            opa  = SpinAngular.TwoParticleOperator(0, plus, true)
            waG2 = SpinAngular.computeCoefficients(opa, basis.csfs[idx_csf[r]], basis.csfs[idx_csf[s]], subshellList)
            wa   = [waG1, waG2]
            #
            me = 0.
            for  coeff in waG1
                jj = Basics.subshell_2j(basis.orbitals[coeff.a].subshell)
                me = me + coeff.T * sqrt( jj + 1) * RadialIntegrals.GrantIab(basis.orbitals[coeff.a], basis.orbitals[coeff.b], grid, potential)
                if  settings.qedModel != NoneQed()  
                    me = me + InteractionStrengthQED.qedLocal(basis.orbitals[coeff.a], basis.orbitals[coeff.b], nm, settings.qedModel, meanPot, grid)  
                end
            end

            for  coeff in waG2
                if  typeof(settings.eeInteractionCI) in [DiagonalCoulomb, CoulombInteraction, CoulombBreit, CoulombGaunt]
                    me = me + coeff.V * InteractionStrength.XL_Coulomb(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                                 basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid, keep=keep)
                end
                                                                                        
                if      typeof(settings.eeInteractionCI) in [BreitInteraction, CoulombBreit, CoulombGaunt]
                    me = me + coeff.V * InteractionStrength.XL_Breit(coeff.nu, basis.orbitals[coeff.a], basis.orbitals[coeff.b],
                                                                               basis.orbitals[coeff.c], basis.orbitals[coeff.d], grid,
                                                                               settings.eeInteractionCI, keep=keep)     
                end
            end
            matrix[r,s] = me
        end
    end 
    if printout    println("   ... done.")    end

    return( matrix )
end

end # module
