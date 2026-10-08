

# The finite-temperature AVERAGE-ATOM path. Not reached by SelfConsistent.performSCF -- it is driven
# from Plasma.AverageAtomScheme (module-Plasma-inc-average-atom.jl), and kept here beside the other
# SCF drivers because it iterates orbitals in the same way.

"""
`SelfConsistent.determineChemicalPotential(orbitals::Dict{Subshell, Orbital}, temp::Float64, radiusWS::Float64,
                                           nm::Nuclear.Model, grid::Radial.Grid)`
    ... determines the chemical potential mu from the neutrality condition  Sum_i (2j_i+1) f(eps_i, mu, temp) = Z,
        where f is the Fermi-Dirac factor; a safeguarded Newton iteration inside a bracket is used, and a
        chemMu::Float64 is returned.

        Note: this general finite-temperature Fermi-Dirac root-finding utility was moved here from module Plasma
              (where it originated as `determineChemicalPotential`), since Plasma.perform(::AverageAtomScheme,
              ...) needs SelfConsistent.solveAverageAtomField below, and solveAverageAtomField itself needs this
              function internally at every SCF iteration; keeping it in Plasma would have made the two modules
              depend on each other circularly. Nothing here is Plasma-scheme-specific.
"""
function determineChemicalPotential(orbitals::Dict{Subshell, Orbital}, temp::Float64, radiusWS::Float64, nm::Nuclear.Model,
                                    grid::Radial.Grid)
    # THE FERMI FACTOR IS EVALUATED IN THE STABLE BRANCH, which is also what makes the derivative provable.
    # Writing 1/(exp(w)+1) directly overflows for large positive w, and the previous version avoided that by CLAMPING
    # w at 300 -- which silently changes the function whose root is being sought.  Choosing the branch instead is
    # exact everywhere, and then  df/dmu = f (1-f) / temp  follows in one line and cannot be mistyped.
    fermi(w::Float64)   = w > 0. ? exp(-w) / (1. + exp(-w))  :  1. / (1. + exp(w))
    function g(mu::Float64)
        wa = -nm.Z
        for  (k,v)  in orbitals
            wa = wa + (Basics.twice(Basics.subshell_j(k)) + 1) * fermi( (v.energy - mu) / temp )
        end
        return( wa )
    end
    function gprime(mu::Float64)
        wa = 0.
        for  (k,v)  in orbitals
            wf = fermi( (v.energy - mu) / temp )
            wa = wa + (Basics.twice(Basics.subshell_j(k)) + 1) * wf * (1. - wf) / temp
        end
        return( wa )
    end
    # THE ROOT IS UNIQUE AND BRACKETABLE: every Fermi factor increases with mu, so g is strictly increasing, g -> -Z
    # as mu -> -inf and g -> (sum of all 2j+1) - Z > 0 as mu -> +inf whenever the subshell set can hold Z electrons.
    # A bracket therefore always exists, and finding it first is what stops the iteration running away.
    #
    # WHY THAT MATTERS HERE, measured 08-Oct-2026 on the version this replaces.  Its Newton derivative read
    # `occ * wc^2 / temp / (wc+1)^2` with wc = exp(w) -- exp(2w) where the derivative needs exp(w).  Against a
    # central difference of g itself, at its own starting point mu = -0.1 for Si at T = 10 eV: finite difference
    # 6.267236e-04, correct analytic 6.267246e-04, and that expression 2.406385e-08 -- FOUR ORDERS too small.  The
    # first step went to mu = -2.49e+08 instead of -9.57e+03 and never recovered, ending near -6.4e+24, whereupon
    # every occupation is zero, the electron sum is 0 instead of Z, and the mean charge comes out as exactly Z.
    # The error was worst where it mattered: for a deeply bound orbital exp(w) is tiny and squaring it annihilates
    # the derivative.  A fudge of -0.0011 had been added to the result "for stability", which is what a misbehaving
    # root-finder looks like from the outside.
    epsLo = minimum( v.energy  for (k,v) in orbitals );   epsHi = maximum( v.energy  for (k,v) in orbitals )
    muLo  = epsLo - 60temp - 1.0;    muHi = epsHi + 60temp + 1.0
    nExp  = 0
    while  g(muLo) > 0.  &&  nExp < 200    muLo = muLo - max(1.0, abs(muLo));    nExp = nExp + 1    end
    while  g(muHi) < 0.  &&  nExp < 400    muHi = muHi + max(1.0, abs(muHi));    nExp = nExp + 1    end
    if  g(muLo) > 0.  ||  g(muHi) < 0.
        error("SelfConsistent.determineChemicalPotential(): no bracket for the neutrality condition between " *
              @sprintf("%.3e", muLo) * " and " * @sprintf("%.3e", muHi) * " Ha;  g = " * @sprintf("%.3e", g(muLo)) *
              " and " * @sprintf("%.3e", g(muHi)) * ".  The subshell set holds " *
              string(sum(Basics.twice(Basics.subshell_j(k)) + 1  for (k,v) in orbitals)) * " electrons against " *
              "Z = $(nm.Z);  raise nMax or lMax if that is less than Z.")
    end
    # A SAFEGUARDED NEWTON: the step is taken when it stays inside the bracket and bisection is taken when it does
    # not, so the iteration inherits Newton's speed without its ability to leave the interval.  The exit is on the
    # RESIDUAL of the neutrality condition itself -- electrons, a quantity with a meaning -- and not on the step.
    chemMu = 0.5 * (muLo + muHi);     nx = 0;    gNow = g(chemMu)
    for  it = 1:200
        nx = it
        if  abs(gNow) < 1.0e-10 * max(1.0, nm.Z)    break    end
        if  gNow > 0.    muHi = chemMu    else    muLo = chemMu    end
        gp    = gprime(chemMu)
        newMu = gp > 0. ? chemMu - gNow / gp : 0.5 * (muLo + muHi)
        if  !(muLo < newMu < muHi)    newMu = 0.5 * (muLo + muHi)    end
        if  abs(newMu - chemMu) < 1.0e-14 * max(1.0, abs(chemMu))    chemMu = newMu;    break    end
        chemMu = newMu;    gNow = g(chemMu)
    end
    if  abs(gNow) > 1.0e-6 * max(1.0, nm.Z)
        Defaults.warn(AddWarning(), "SelfConsistent.determineChemicalPotential(): the neutrality condition is " *
                      "satisfied only to " * @sprintf("%.2e", gNow) * " electrons after $nx iterations.")
    end
    println(">>> chemical potential: mu = " * @sprintf("%.8f", chemMu) * " Ha after $nx safeguarded Newton steps;  " *
            "the neutrality residual is " * @sprintf("%.2e", gNow) * " electrons.")

    return( chemMu )
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

    # Defaults.setDefaults("standard grid", primitives.grid; printout=printout)
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
                # Collected as well: in a long run this line scrolls away, and nobody learns afterwards that a
                # field never converged.  The accuracy is rounded so that repeated identical failures collapse
                # into one counted entry; see Defaults.warn.
                Defaults.warn(AddWarning(), "SelfConsistent.solveAverageAtomField(): the SCF did NOT converge -- " *
                              "stopped at accuracy " * @sprintf("%.1e", accuracyScf) * " after 32 iterations.")
            break
        end
        if  printout    println("\nIteration $NoIteration for symmetries ... ")    end
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
                    # println("  $sh  en [a.u.] = $(newOrbital.energy)   self-consistency = $(wcOrbital), $(wcBlock) [kappa=$kappa] ") 
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
