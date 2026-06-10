
using  QuadGK, ..Basics,  ..Defaults, ..Nuclear, ..Radial, ..Math
export Model


"""
`Plasma.computeElectronNumberDensity(grid::Radial.Grid, orbitals::Dict{Subshell, Orbital}, chemMu::Float64, temp::Float64)`  
    ... computes the electron number density ne from the one-particle energies of the electrons, the
        chemical potential chemMu [a.u.] and the temperature [a.u.]. A triple of three Array{Float64,1}'s
        ( totalNe, posNe, negNe ) which refer to the total n_e (r) and the contributions of the positive and negative
        part of the spectrum, respectively, all with regard to the given grid.
"""
function computeElectronNumberDensity(grid::Radial.Grid, orbitals::Dict{Subshell, Orbital}, chemMu::Float64, temp::Float64)
    totalNe = zeros(length(grid.r));   posNe = zeros(length(grid.r));   negNe = zeros(length(grid.r))
    for (k,v) in orbitals
        wa  = Basics.FermiDirac(v.energy, chemMu, temp)
        occ = (Basics.twice(Basics.subshell_j(k)) + 1) 
        for (ir, r)  in  enumerate(grid.r)
            if size(v.P,1) < ir  ||  size(v.Q,1) < ir    break   end
                                totalNe[ir] = totalNe[ir]  +  wa * occ * (v.P[ir]^2 + v.Q[ir]^2)
            if  v.energy < 0.  negNe[ir]   = negNe[ir]    +  wa * occ * (v.P[ir]^2 + v.Q[ir]^2)
            else               posNe[ir]   = posNe[ir]    +  wa * occ * (v.P[ir]^2 + v.Q[ir]^2)
            end
        end
    end
    
    return ( totalNe, posNe, negNe )
end


"""
`Plasma.computeFormFactors(qValues::Array{Float64,1}, orbitals::Dict{Subshell, Orbital}, chemMu::Float64, temp::Float64, 
                            grid::Radial.Grid)`  
    ... computes the (standard) form factor F(q) for the electron density as given by the AA orbitals.
        A list of formFactors::Array{Float64,1} is returned that directly refers to the given q-values.
"""
function computeFormFactors(qValues::Array{Float64,1}, orbitals::Dict{Subshell, Orbital}, chemMu::Float64, temp::Float64, 
                            grid::Radial.Grid)
    formFactors = Float64[];    wa = zeros( grid.NoPoints )
    #
    for  q in qValues
        for (k,v) in orbitals
            occ  = Basics.FermiDirac(v.energy, chemMu, temp) * (Basics.twice(Basics.subshell_j(k)) + 1) 
            nrho = length(v.P)
            for    i = 1:nrho   wa[i] = wa[i] + occ * (v.P[i]^2 + v.Q[i]^2)    end
        end
        # Compute the full integrant; the factor 4pi * r^2 is already in the density
        if     q == 0.    error("q = 0. is not supported for form-factor computations.")
        else   for    i = 2:grid.NoPoints   wa[i] = wa[i] / grid.r[i] * sin(q*grid.r[i]) / q    end
        end
        sF = RadialIntegrals.V0(wa, grid.NoPoints, grid)
        push!(formFactors, sF)
    end
    
    println("\n Form factors: \n" *
            "\n     q [a_o]     form F [a.u.]     " *
            "\n ----------------------------------")
    for  (iq, q)  in  enumerate(qValues)
        sa = "   "     * @sprintf("%.4e", q) * "     "     * @sprintf("%.4e", formFactors[iq])
        println(sa)
    end
    println("  ", TableStrings.hLine(32), "\n")
    
    return ( formFactors )
end


"""
`Plasma.computeMeanCharge(nm::Nuclear.Model, orbitals::Dict{Subshell, Orbital}, chemMu::Float64, temp::Float64)`  
    ... computes the mean charge state of an average-atom ion in the plasma of given density and temperature;
        i.e. Z - sum occ (epsilon > 0); a mean charge state meanCharge::Float64 is returned.
"""
function computeMeanCharge(nm::Nuclear.Model, orbitals::Dict{Subshell, Orbital}, chemMu::Float64, temp::Float64)
    meanCharge = nm.Z
    for (k,v) in orbitals
        if  v.energy > 0. 
            occ        = Basics.FermiDirac(v.energy, chemMu, temp) * (Basics.twice(Basics.subshell_j(k)) + 1) 
            meanCharge = meanCharge - occ
        end
    end
    println(">> Mean charge = $meanCharge ")
    
    return ( meanCharge )
end


"""
`Plasma.determineChemicalPotential(orbitals::Dict{Subshell, Orbital}, temp::Float64, radiusWS::Float64, nm::Nuclear.Model,
                                    grid::Radial.Grid)`  
    ... determines the chemical potential so that Sum_i f(epsilon_i, mu, temp) = Z. 
        The Newton-Raphson methods is used to iterate to the chemical potential; a chemMu::Float64 is returned.
"""
function determineChemicalPotential(orbitals::Dict{Subshell, Orbital}, temp::Float64, radiusWS::Float64, nm::Nuclear.Model,
                                    grid::Radial.Grid)
    function g(mu::Float64, orbitals::Dict{Subshell, Orbital}, temp::Float64, nm::Nuclear.Model)
        wa = - nm.Z
        for  (k,v)  in orbitals
            occ = Basics.twice( Basics.subshell_j(k)) + 1
            ## occ = occ * Plasma.finiteNorm(v, radiusWS, grid)
            wb  = (v.energy - mu) /temp
            if  wb > 300.   wb = 300.   end
            ## @show occ, wb, wa
            wa  = wa + occ / (exp(wb) + 1)    
        end
        return( wa )
    end
    ## for  mu = -18.1:0.005:-18.0  @show  mu, g(mu::Float64, orbitals::Dict{Subshell, Orbital}, temp::Float64, nm::Nuclear.Model)  end
    ## error("stop x")
    #
    function gprime(mu::Float64, orbitals::Dict{Subshell, Orbital}, temp::Float64, nm::Nuclear.Model)
        wa = 0.
        for  (k,v)  in orbitals
            occ = Basics.twice( Basics.subshell_j(k)) + 1
            ## occ = occ # * Plasma.finiteNorm(v, radiusWS, grid)
            wb  = (v.energy - mu) /temp
            if  wb > 300.   wb = 300.   end
            wc  = exp( wb )
            wa  = wa + occ * wc^2 / temp / (wc+1)^2    
            ## @show occ, wb, wc, wa
        end
        return( wa )
    end
    # Iterate for the chemical potential
    chemMu = -0.1;     newMu = 0.;     nx = 0
    while true
        nx = nx + 1
        newMu = chemMu - g(chemMu, orbitals, temp, nm) / gprime(chemMu, orbitals, temp, nm)
        if  abs(newMu - chemMu) < 1.0e-4  break
        else  ## println(">>> Newton-Raphson: $nx)  chemMu = $chemMu  g = $(g(chemMu, orbitals, temp, nm)) ")
                chemMu = newMu
        end
    end
    
    chemMu = chemMu - 0.0011  ## Seems to bring better stability in the SCF computations
    
    println(">>> Newton-Raphson: $nx)  chemMu = $chemMu  g = $(g(chemMu, orbitals, temp, nm)) ")
    return ( chemMu )
end


"""
`Plasma.determineWignerSeitzRadius(rho::Float64, nm::Nuclear.Model)`  
    ... determines the Wigner-Seitz radius R^(WS) from the plasma density rho [g/cm^3] ne and the nuclear charge Z.
"""
function determineWignerSeitzRadius(rho::Float64, nm::Nuclear.Model)
    
    # Convert the density into rho [A_Z * u/a_o^3]
    wa = Defaults.convertUnits("density: from [g/cm^3] to atomic", rho) / nm.mass
    wb = 4pi * wa
    wr = (3.0 / wb)^(1/3)
    
    println(">>> Wigner-Seitz radius R^(WS) = $wr  [a_o] for density $rho [g/cm^3] and atomic mass $(nm.mass).")
    
    return ( wr )
end


"""
`Plasma.displayElectronNumberDensity(grid::Radial.Grid, orbitals::Dict{Subshell, Orbital}, chemMu::Float64, 
                                        temp::Float64, radiusWS::Float64)`  
    ... displays the electron number density 4pi r^2 ne from the one-particle energies of the electrons, the
        chemical potential chemMu [a.u.] and the temperature [a.u.]. It tabulates the radial grid, total density as well
        as the contributions from orbital with negative and positive binding energies separately for a fixed stepsize in r/a_o.
        Nothing is returned from this procedure
"""
function displayElectronNumberDensity(grid::Radial.Grid, orbitals::Dict{Subshell, Orbital}, chemMu::Float64, 
                                        temp::Float64, radiusWS::Float64)
    totalNe, posNe, negNe = Plasma.computeElectronNumberDensity(grid, orbitals, chemMu, temp)
    
    println("\n Radial electron number density: \n" *
            "\n     r [a_o]       total n_e     n_e  (eps < 0.)   n_e  (eps > 0.)        " *
            "\n -------------------------------------------------------------------------")
    for  (ir, r)  in  enumerate(grid.r)
        if  r < 0.1         continue    end
        if  rem(ir ,5) != 1 continue    end
        if  r > radiusWS    break       end
        sa = "   "     * @sprintf("%.4e", grid.r[ir]) * "     "     * @sprintf("%.4e", totalNe[ir]) * 
                "       " * @sprintf("%.4e", negNe[ir])  * "         " * @sprintf("%.4e", posNe[ir])
        println(sa)
    end
    println("  ", TableStrings.hLine(70), "\n")
    
    return ( nothing )
end


"""
`Plasma.finiteNorm(a::Radial.Orbital, radiusWS::Float64, grid::Radial.Grid)`   
    ... computes the (finite norm) integral of two radial orbital functions.
    
        norm = int_0^radiusWS  dr  [P_a^2 + Q_a^2]
"""
function finiteNorm(a::Radial.Orbital, radiusWS::Float64, grid::Radial.Grid)
    mtp = size(a.P, 1)
    # Distinguish the radial integration for different grid definitions
    if  grid.meshType == Radial.MeshGL()
        wa = 0.
        for  i = 2:mtp  
            if  grid.r[i] > radiusWS    break   end
            wa = wa + (a.P[i]^2 + a.Q[i]^2) * grid.wr[i]   
        end
        return( wa )
    else
        error("stop a")
    end
end


    
"""
`Plasma.perform(scheme::Plasma.AverageAtomScheme, computation::Plasma.Computation; output::Bool=true)`  
    ... to perform an average-atom plasma computation for an atoms with nuclear charge Z which generates a self-consistent set 
        of orbitals. For output=true, dictionary is returned from which the relevant results can be can easily accessed by 
        proper keys.
"""
function  perform(scheme::Plasma.AverageAtomScheme, computation::Plasma.Computation; output::Bool=true)
    if  output    results = Dict{String, Any}()    else    results = nothing    end
        
    nm   = computation.nuclearModel
    RWS  = Plasma.determineWignerSeitzRadius(computation.settings.density, nm);   @show RWS
    wa   = Bsplines.generatePrimitives(computation.grid)
    temp = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", computation.settings.temperature)
    println(" ")

    # Generate a subshell list that is taken into account into the average atom computations
    shells    = Basics.generateShellList(1, scheme.nMax, scheme.lMax)
    subshells = Basics.generateSubshellList(shells)
    println(">> Subshells included into the average-atom scheme: \n   $subshells")
    # Generate hydrogenic orbitals for all subshells and display comparison
    orbitals  = Bsplines.generateOrbitalsHydrogenic(wa, nm, subshells; printout=false)
    chemMu    = Plasma.determineChemicalPotential(orbitals, temp, RWS, nm, computation.grid)
    Basics.displayOrbitalProperties(stdout, orbitals, chemMu, temp, RWS, nm, computation.grid)
    # Solve the orbitals and chemical potential self-consistently in the average-atom model
    orbitals  = Bsplines.solveSelfConsistentAverageAtom(wa, nm, orbitals, temp, RWS, scheme.scField, printout=true)
    chemMu    = Plasma.determineChemicalPotential(orbitals, temp, RWS, nm, computation.grid)
    # Diplay orbital properties
    Basics.displayOrbitalProperties(stdout, orbitals, chemMu, temp, RWS, nm, computation.grid)
    Plasma.displayElectronNumberDensity(computation.grid, orbitals, chemMu, temp, RWS)
    # Generate electron number densities and mean charge state
    totalNe, posNe, negNe = Plasma.computeElectronNumberDensity(computation.grid, orbitals, chemMu, temp)
    meanCharge            = Plasma.computeMeanCharge(nm, orbitals, chemMu, temp)
    # Return results if required
    if  output   
        results["chemical mu"]   = chemMu;                  results["mean charge"]   = meanCharge
        results["AA orbitals"]   = orbitals;                results["density n_e"]   = totalNe
        results["negative n_e"]  = negNe;                   results["positive n_e"]  = posNe;                  
    end
    #
    #
    # Calculate photoionization data and cross sections
    if  scheme.calcPhotoionizationCs
        piData = Plasma.computePhotoionizationData(scheme.piSubshells, orbitals, chemMu, temp, computation.grid)
        Plasma.displayPhotoionizationCrossSections(scheme.omegas, piData)
        # Add omegas and cross sections to results ... computePhotoionizationCrossSections(subshell, scheme.omegas, piData)
    end
    #
    # Calculate form factors
    if  scheme.calcFormFactor
        formF = Plasma.computeFormFactors(scheme.qValues, orbitals, chemMu, temp, computation.grid)
        # Add q-values and form factors to results 
    end
    #
    # Calculate form factors
    if  scheme.calcScatteringFactor
        if !scheme.calcPhotoionizationCs    
            error("Scattering factors also require the computation of photoionization cross sections")    end
        scatteringF = Plasma.computeScatteringFactors(scheme.omegas, piData, orbitals, chemMu, temp)
        # Add scattering factors to results 
    end

    
    Defaults.warn(PrintWarnings())
    Defaults.warn(ResetWarnings())
    return( results )
end



