
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
    formFactors = Float64[]
    # The (Fermi-Dirac occupation-weighted) electron density does not depend on q; compute it once here, and take a
    # fresh copy for each q below -- reusing (and mutating) a single array across the q loop previously accumulated
    # the density afresh on top of the already q-transformed values of the *previous* iteration, silently corrupting
    # every form factor but the first.
    density = zeros( grid.NoPoints )
    for (k,v) in orbitals
        occ  = Basics.FermiDirac(v.energy, chemMu, temp) * (Basics.twice(Basics.subshell_j(k)) + 1)
        nrho = length(v.P)
        for    i = 1:nrho   density[i] = density[i] + occ * (v.P[i]^2 + v.Q[i]^2)    end
    end
    #
    for  q in qValues
        # Compute the full integrant; the factor 4pi * r^2 is already in the density
        if     q == 0.    error("q = 0. is not supported for form-factor computations.")
        else   wa = deepcopy(density)
               for    i = 2:grid.NoPoints   wa[i] = wa[i] / grid.r[i] * sin(q*grid.r[i]) / q    end
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
    ... computes the mean charge of an average-atom ion from a given set of orbitals, i.e. the number of electrons
        that are NOT bound to it, `Z* = sum_(eps>0) (2j+1) f(eps,mu,T)`;  a `meanCharge::Float64` is returned.

        **IT USED TO RETURN `Z - sum_(eps>0) occ f`, WHICH IS THE COMPLEMENT -- the number of BOUND electrons --
        under the name of the mean charge.**  Corrected 08-Oct-2026:  for aluminium at solid density and 10 eV it
        printed 9.986 where the mean charge is 3.014, so the error was not small but exactly the complement, and
        at neutrality the two are indistinguishable by eye only when Z/2 electrons happen to be free.

        This form sums over whatever orbitals it is handed, so it is as complete as that set is.
        `Plasma.perform(::AverageAtomScheme, ...)` does NOT use it:  it takes the mean charge from
        `SelfConsistent.solveAverageAtomSpectra`, which counts the FULL per-kappa spectrum of the cell -- about
        seventy-eight states per kappa against the half-dozen a subshell list carries -- and whose count satisfies
        neutrality exactly.  The two agree only where the subshell list happens to hold the whole continuum.
"""
function computeMeanCharge(nm::Nuclear.Model, orbitals::Dict{Subshell, Orbital}, chemMu::Float64, temp::Float64)
    meanCharge = 0.
    for (k,v) in orbitals
        if  v.energy > 0.
            meanCharge = meanCharge + Basics.FermiDirac(v.energy, chemMu, temp) * (Basics.twice(Basics.subshell_j(k)) + 1)
        end
    end
    println(">> Mean charge Z* = $meanCharge  from the given " * string(length(orbitals)) * " subshells.")

    return ( meanCharge )
end


"""
`Plasma.computePhotoionizationData(piSubshells::Array{Subshell,1}, orbitals::Dict{Subshell, Orbital}, chemMu::Float64,
                                   temp::Float64, grid::Radial.Grid, pot::Radial.Potential, omegas::Array{Float64,1})`
    ... computes the bound-free photoionization cross section sigma_bf(omega) of each of the given bound subshells in
        the converged average-atom potential, for every photon energy in omegas [a.u.];  a
        Dict{Subshell, Array{Float64,1}} of cross sections [a.u.] is returned, one entry per subshell and one value
        per omega.

        **THIS IS EQ. (23) OF JOHNSON & NILSEN, HEDP 31 (2019) 92**, `sigma_bf = (8 pi^2/3) alpha omega |D|^2` with
        the radial dipole matrix element `D = int P_eps,l(r) r P_b(r) dr` of Eq. (24).  The continuum electron takes
        the energy `eps = omega + eps_b` left over after the binding energy, so a photon below threshold returns zero.

        **THE FRACTIONAL OCCUPATION IS APPLIED HERE AND IS NOT A DETAIL.**  In a plasma the bound subshell is only
        partly filled, and Johnson multiplies the cross section by the fractional occupation -- occ/2 for the K shell,
        which is just the Fermi-Dirac factor f(eps_b) itself.  Their Fig. 7 shows that factor accounting for a 95 %
        reduction of the opacity between 20 and 500 eV while the cross section proper contributes 18 %, so a
        sigma_bf computed without it is not approximately right but wrong by more than an order of magnitude.

        **THE ANGULAR FACTOR IS THE s -> p ONE.**  The 8 pi^2/3 of Eq. (23) already carries the angular algebra of a
        K-shell (l = 0 -> l = 1) transition, which is what the paper computes and what this reproduces.  For a
        subshell with l > 0 the two channels l-1 and l+1 would each need their own weight, so such a subshell is
        REFUSED here rather than silently given the K-shell factor.  The two relativistic partial waves p_1/2 and
        p_3/2 are combined with their statistical weights, which returns the single non-relativistic l = 1 channel
        in the limit where they coincide.
"""
function computePhotoionizationData(piSubshells::Array{Subshell,1}, orbitals::Dict{Subshell, Orbital}, chemMu::Float64,
                                    temp::Float64, grid::Radial.Grid, pot::Radial.Potential, omegas::Array{Float64,1})
    alpha    = Defaults.getDefaults("alpha")
    piData   = Dict{Subshell, Array{Float64,1}}()
    cSettings = Continuum.Settings(false, grid.NoPoints - 100)
    for  sh  in  piSubshells
        if  !haskey(orbitals, sh)   error("Plasma.computePhotoionizationData(): subshell $sh is not among the " *
                                          "average-atom orbitals;  check scheme.piSubshells against nMax and lMax.")   end
        if  Basics.subshell_l(sh) != 0
            error("Plasma.computePhotoionizationData(): subshell $sh has l > 0, and Eq. (23) of Johnson & Nilsen " *
                  "(2019) carries the angular factor of an s -> p transition only.  Only s subshells are supported.")
        end
        bOrb   = orbitals[sh];      epsB = bOrb.energy
        fOcc   = Basics.FermiDirac(epsB, chemMu, temp)       ## the fractional occupation, occ/(2j+1)
        sigmas = zeros(length(omegas))
        for  (io, omega)  in  enumerate(omegas)
            eps = omega + epsB                               ## epsB < 0, so this is the photoelectron energy
            if  eps <= 0.    continue    end                 ## below threshold
            dSum = 0.;    wSum = 0.
            for  kappa  in  (1, -2)                          ## p_1/2 and p_3/2
                wk = abs(2kappa)                             ## 2j+1
                cOrb, phase = Continuum.generateOrbitalLocalPotential(eps, Subshell(101, kappa), pot, cSettings)
                mtp = min(length(cOrb.P), length(bOrb.P), length(grid.r))
                dip = 0.
                for  ir = 2:mtp    dip = dip + grid.wr[ir] * cOrb.P[ir] * grid.r[ir] * bOrb.P[ir]    end
                dSum = dSum + wk * dip^2;    wSum = wSum + wk
            end
            sigmas[io] = 8pi^2/3 * alpha * omega * (dSum/wSum) * fOcc
        end
        piData[sh] = sigmas
    end

    return( piData )
end


"""
`Plasma.displayPhotoionizationCrossSections(omegas::Array{Float64,1}, piData::Dict{Subshell, Array{Float64,1}})`
    ... displays the bound-free photoionization cross sections of every subshell in piData, in barn and as a
        function of the photon energy; nothing is returned.
"""
function displayPhotoionizationCrossSections(omegas::Array{Float64,1}, piData::Dict{Subshell, Array{Float64,1}})
    println("\n  Bound-free photoionization cross sections of the average atom:\n")
    println("  " * "-"^88)
    print(  "   omega [eV]  ");   for sh in sort(collect(keys(piData)), by=x->string(x))   print("   sigma_bf($sh) [barn]")   end
    println("\n  " * "-"^88)
    for  (io, omega)  in  enumerate(omegas)
        print("  " * @sprintf("%10.2f", Defaults.convertUnits("energy: from atomic to eV", omega)) * "  ")
        for  sh  in  sort(collect(keys(piData)), by=x->string(x))
            print("   " * @sprintf("%18.4e", Defaults.convertUnits("cross section: from atomic to barn", piData[sh][io])))
        end
        println("")
    end
    println("  " * "-"^88)

    return( nothing )
end


"""
`Plasma.computeScatteringFactors(omegas::Array{Float64,1}, piData::Dict{Subshell, Array{Float64,1}}, nm::Nuclear.Model)`
    ... computes the imaginary part f_2 of the atomic scattering factor and the x-ray mass attenuation coefficient
        mu/rho from the bound-free cross sections; a Tuple (f2::Array{Float64,1}, muOverRho::Array{Float64,1}) is
        returned, with f_2 dimensionless and mu/rho in cm^2/g.

        **EQS. (21) AND (27) OF JOHNSON & NILSEN (2019)**: `f_2 = sigma_bf / (2 r_0 lambda)` relates the imaginary
        scattering factor to the cross section through the optical theorem, and `mu/rho = (N_A/A) sigma_bf` is the
        mass attenuation coefficient that the x-ray intensity falls off with, `I = I_0 exp(-mu z)`.  Both are summed
        over the subshells present in piData, since each contributes additively to the absorption.
"""
function computeScatteringFactors(omegas::Array{Float64,1}, piData::Dict{Subshell, Array{Float64,1}}, nm::Nuclear.Model)
    alpha = Defaults.getDefaults("alpha");     r0 = alpha^2        ## classical electron radius in a.u.
    f2 = zeros(length(omegas));     muOverRho = zeros(length(omegas))
    for  (io, omega)  in  enumerate(omegas)
        sigTot = 0.;    for (sh, sg) in piData    sigTot = sigTot + sg[io]    end
        if  omega > 0.
            lambda  = 2pi / (alpha * omega)                         ## lambda = 2 pi c / omega, c = 1/alpha
            f2[io]  = sigTot / (2 * r0 * lambda)
            # mu/rho = (N_A/A) sigma, with sigma converted from a.u. to cm^2 and A in g/mol
            sigCm2        = Defaults.convertUnits("cross section: from atomic to barn", sigTot) * 1.0e-24
            muOverRho[io] = 6.02214076e23 / nm.mass * sigCm2
        end
    end

    return( (f2, muOverRho) )
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
    wa = 0.
    for  i = 2:mtp  
        if  grid.r[i] > radiusWS    break   end
        wa = wa + (a.P[i]^2 + a.Q[i]^2) * grid.wr[i]   
    end
    return( wa )
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
    RWS  = Plasma.determineWignerSeitzRadius(computation.settings.density, nm);
    wa   = Bsplines.generatePrimitives(computation.grid)
    temp = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", computation.settings.temperature)
    println(" ")

    # Generate a subshell list that is taken into account into the average atom computations
    shells    = Basics.generateShellList(1, scheme.nMax, scheme.lMax)
    subshells = Basics.generateSubshellList(shells)
    println(">> Subshells included into the average-atom scheme: \n   $subshells")
    # Generate hydrogenic orbitals for all subshells and display comparison
    orbitals  = Bsplines.generateOrbitalsHydrogenic(subshells, nm, wa; printout=false)
    chemMu    = SelfConsistent.determineChemicalPotential(orbitals, temp, RWS, nm, computation.grid)
    Basics.displayOrbitalProperties(stdout, orbitals, chemMu, temp, RWS, nm, computation.grid)
    # Solve the orbitals, the density and the chemical potential self-consistently in the average-atom model.
    # ONE CALL, AND ITS NUMBERS ARE THE ONES REPORTED.  Until 08-Oct-2026 the chemical potential was solved a
    # second time here, from the subshell list, and that value -- not the converged one -- went into the results
    # and into every property computed below.  The two differ by whatever the subshell list leaves out of the
    # continuum, which is most of it.
    wScf      = SelfConsistent.solveAverageAtomSpectra(orbitals, nm, scheme.scField, temp, RWS, wa; printout=true)
    orbitals  = wScf.orbitals;    chemMu = wScf.chemMu;    meanCharge = wScf.nCont
    # Diplay orbital properties
    Basics.displayOrbitalProperties(stdout, orbitals, chemMu, temp, RWS, nm, computation.grid)
    Plasma.displayElectronNumberDensity(computation.grid, orbitals, chemMu, temp, RWS)
    # Generate electron number densities and mean charge state
    totalNe, posNe, negNe = Plasma.computeElectronNumberDensity(computation.grid, orbitals, chemMu, temp)
    # The table above is built from the SUBSHELL LIST and is therefore as incomplete as that list; say by how
    # much, so that a reader does not take it for the density the field was built from.  That one is in
    # results["radial density"] and integrates to Z.
    wSub      = RadialIntegrals.V0(totalNe, computation.grid.NoPoints, computation.grid)
    println(">> Mean charge Z* = " * @sprintf("%.5f", meanCharge) * " free and " * @sprintf("%.5f", wScf.nBound) *
            " bound electrons at mu = " * @sprintf("%.6f", chemMu) * " Ha, from the full per-kappa spectrum " *
            "(lMax = $(wScf.lMax), free partial waves to $(wScf.lMaxFree)).")
    println("   The subshell-list density tabulated above accounts for " * @sprintf("%.4f", wSub) * " of the " *
            "$(nm.Z) electrons;  the density the field was built from is results[\"radial density\"].")
    # Return results if required
    if  output   
        results["chemical mu"]   = chemMu;                  results["mean charge"]   = meanCharge
        results["AA orbitals"]   = orbitals;                results["density n_e"]   = totalNe
        results["negative n_e"]  = negNe;                   results["positive n_e"]  = posNe;
        results["bound electrons"] = wScf.nBound;           results["radial density"] = wScf.rhot
    end
    #
    #
    # Calculate photoionization data and cross sections
    if  scheme.calcPhotoionizationCs
        # THE CONVERGED AVERAGE-ATOM POTENTIAL IS TAKEN FROM THE SCF, not rebuilt.  The continuum electron must be
        # generated in the SAME field the bound orbital sits in -- a continuum wave from a bare nuclear potential
        # would carry no screening and give a cross section that is wrong wherever the plasma is dense -- and
        # until 08-Oct-2026 it was rebuilt here from the subshell list, which is a DIFFERENT and incomplete field
        # from the one the orbitals were converged in.
        totalPot = wScf.pot
        piData   = Plasma.computePhotoionizationData(scheme.piSubshells, orbitals, chemMu, temp, computation.grid,
                                                     totalPot, scheme.omegas)
        Plasma.displayPhotoionizationCrossSections(scheme.omegas, piData)
        if  output
            results["pi omegas"]       = scheme.omegas
            results["pi cross sections"] = piData
        end
    end
    #
    # Calculate form factors
    if  scheme.calcFormFactor
        formF = Plasma.computeFormFactors(scheme.qValues, orbitals, chemMu, temp, computation.grid)
        if  output    results["ff q-values"] = scheme.qValues;    results["form factors"] = formF    end
    end
    #
    # Calculate the scattering factor f_2 and the mass attenuation coefficient
    if  scheme.calcScatteringFactor
        if !scheme.calcPhotoionizationCs    
            error("Scattering factors also require the computation of photoionization cross sections")    end
        (f2, muOverRho) = Plasma.computeScatteringFactors(scheme.omegas, piData, nm)
        println("\n  Imaginary scattering factor f_2 and mass attenuation coefficient:\n")
        println("  " * "-"^62)
        println("   omega [eV]            f_2        mu/rho [cm^2/g]")
        println("  " * "-"^62)
        for  (io, om)  in  enumerate(scheme.omegas)
            println("  " * @sprintf("%10.2f", Defaults.convertUnits("energy: from atomic to eV", om)) *
                    @sprintf("%17.5e", f2[io]) * @sprintf("%19.5e", muOverRho[io]))
        end
        println("  " * "-"^62)
        if  output    results["scattering factor f2"] = f2;    results["mass attenuation"] = muOverRho    end
    end

    
    Defaults.warn(PrintWarnings())
    Defaults.warn(ResetWarnings())
    return( results isa Dict{String,Any} ? Basics.PerformResults(results) : results )
end



