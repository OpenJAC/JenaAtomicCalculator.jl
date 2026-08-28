
"""
`module  JAC.Continuum`  
... a submodel of JAC that contains all functions and methods for computing, either in a given potential or with
    regard to a well-defined bound state.
"""
module Continuum


using  GSL, Printf, SpecialFunctions, DelimitedFiles
using  ..Basics, ..Bsplines, ..Defaults, ..ManyElectron, ..Radial, ..Nuclear


"""
`struct  Continuum.Settings`  ... defines a type for the parameters for computing continuum orbitals.

    + includeExchange         ::Bool    ... True, if the exchange is to be included and false otherwise. SEE THE NOTE BELOW:
                                            the only supported value is `false`.
    + mtp                     ::Int64   ... No of grid points for which the continuum orbital(s) are to be computed.

    DELIBERATELY POSITIONAL, AND DELIBERATELY OUTSIDE THE `Settings` CONVENTION -- decided 28-Aug-2026. Every other
    `Settings` in the package carries a zero-argument constructor and a keyword copy-constructor, because every other one is
    something a USER assembles and then varies. This one is not: it is an internal argument-holder that a process module
    fills in on its way to `Continuum.generateOrbitalForLevel`, and it is built positionally at all 27 sites that construct it in
    `src/` and `examples/` -- 20 of them the identical `Continuum.Settings(false, nrContinuum)`. Adding the convention here
    would give a user a constructor for a struct they have no decision to make about. Do not "fix" this; it is the decision.

    `includeExchange` HAS EXACTLY ONE LEGAL VALUE, `false`, and that is worth knowing before anyone sets it hoping for
    exchange. Measured 28-Aug-2026: it is `false` at all 27 construction sites, and the only code that reads it is
    `generateOrbitalLocalPotential`, which does `settings.includeExchange && error(...)` -- so `true` does not enable
    exchange, it raises. This matches the note on the method below: all continuum orbitals are generated in a LOCAL (DFS)
    potential, exchange with the bound electrons is not treated. The field is therefore a placeholder for a capability that
    does not exist yet, not a switch.
"""
struct Settings 
    includeExchange           ::Bool  
    mtp                       ::Int64  
end 


"""
`Continuum.generateOrbitalForLevel(energy::Float64, sh::Subshell, level::Level, nm::Nuclear.Model, grid::Radial.Grid,
                                    settings::Continuum.Settings; nuclearPot::Union{Nothing,Radial.Potential}=nothing,
                                    primitives::Union{Nothing,Bsplines.Primitives}=nothing)`
    ... to generate a continuum orbital for the (continuum) subshell sh, the energy and the effective charge within the given
        potential. The continuum orbital is generated orthogonal with regard to all subshells of the same symmetry in the basis.
        All further specifications about this generations are made by proper settings. Both keywords carry quantities that are
        CONSTANT for a whole computation -- the nuclear potential depends only on the nuclear model and the grid, the B-spline
        primitives only on the grid -- so that a caller which generates many continuum orbitals need not rebuild them for every
        line and every partial wave; omitting them reproduces the previous behaviour exactly. A tupel of a (continuum)
        (orbital::Orbital, phase::Float64) is returned.

        NOTE ON THE POTENTIAL: all continuum orbitals are generated in a LOCAL (DFS) potential -- exchange with
        the bound electrons is not treated. That is a property of this method, true of every run, and it is
        stated here rather than emitted as a per-orbital warning as it was until 12-Aug-2026 (which produced
        one identical entry in jac-warn.report for every partial wave, and said nothing about the run).
"""
function generateOrbitalForLevel(energy::Float64, sh::Subshell, level::Level, nm::Nuclear.Model, grid::Radial.Grid,
                                 settings::Continuum.Settings; nuclearPot::Union{Nothing,Radial.Potential}=nothing,
                                 primitives::Union{Nothing,Bsplines.Primitives}=nothing)
    #
    # Generate a (local) potential for the given level.  The NUCLEAR part depends only on the nuclear model
    # and the grid, so it is constant for a whole computation -- but this function is called once per line
    # AND per partial wave, so rebuilding it here is pure repeated work.  Profiling an AutoIonization run
    # (10-Aug-2026) put ~19% of the wall clock in the nuclear charge density and a further ~16% in the
    # adaptive quadrature integrating it, of which this call site was responsible for 57.6%.
    # Callers may therefore pass a potential built once; nothing is cached here, and omitting the keyword
    # reproduces the previous behaviour exactly.
    nuclearPotential  = isnothing(nuclearPot) ? Nuclear.nuclearPotential(nm, grid) : nuclearPot
    ## wp1 = compute("radial potential: core-Hartree", grid, wLevel)
    ## wp2 = compute("radial potential: Hartree-Slater", grid, wLevel)
    ## wp3 = compute("radial potential: Kohn-Sham", grid, wLevel)
    ## wp4 = compute("radial potential: Dirac-Fock-Slater", grid, wLevel)
    ## wp  = Basics.computePotential(Basics.DFSField(0.42), grid, level)   
    wp  = Basics.computePotential(Basics.DFSField(1.0), grid, level);
    pot = Basics.add(nuclearPotential, wp)
    
    # Generate a continuum orbital due to the given solution method
    if      Defaults.GBL_CONT_SOLUTION  ==  ContBessel() 
        cOrbital = Continuum.generateOrbitalBessel(energy, sh, grid::Radial.Grid, settings)
    elseif  Defaults.GBL_CONT_SOLUTION  ==  AsymptoticCoulomb() 
        cOrbital = Continuum.generateOrbitalAsymptoticCoulomb(energy, sh, pot, settings)
    elseif  Defaults.GBL_CONT_SOLUTION  ==  NonrelativisticCoulomb()
        cOrbital = Continuum.generateOrbitalNonrelativisticCoulomb(energy, sh, pot.grid, settings)
    elseif  Defaults.GBL_CONT_SOLUTION  ==  BsplineGalerkin()
        cOrbital = Continuum.generateOrbitalGalerkin(energy, sh, pot, settings; primitives=primitives)
    else    error("stop a")
    end
    #
    #
    # Normalize the continuum orbital and determine its phase
    if      Defaults.GBL_CONT_NORMALIZATION  ==  PureSineNorm()
        cOrbital, phase = Continuum.normalizeOrbitalPureSine(cOrbital, pot.grid, settings)
    elseif  Defaults.GBL_CONT_NORMALIZATION  ==  CoulombSineNorm()
        cOrbital, phase = Continuum.normalizeOrbitalCoulombSine(cOrbital, pot, settings)
    elseif  Defaults.GBL_CONT_NORMALIZATION  ==  OngRussekNorm()
        cOrbital, phase = Continuum.normalizeOrbitalOngRussek(cOrbital, pot, settings)
    elseif  Defaults.GBL_CONT_NORMALIZATION  ==  AlokNorm()
        cOrbital, phase = Continuum.normalizeOrbitalAlok(cOrbital, pot, settings)
    else    error("stop b")
    end

    #=open("continuum.txt","w") do file
        write(file,"energy $energy , n $(sh.n), kappa $(sh.kappa), phase $phase ")
        write(file,"potential")
        writedlm(file, pot.Zr)
        write(file,"radial grid")
        writedlm(file, pot.grid.r)
    end=#
    
    return( cOrbital, phase )
end


"""
`Continuum.generateOrbitalForLevel(energy::Float64, sh::Subshell, level::Level, nm::Nuclear.Model, grid::Radial.Grid,
                                    settings::Continuum.Settings, plasmaModel::Basics.AbstractPlasmaModel;
                                    nuclearPot::Union{Nothing,Radial.Potential}=nothing,
                                    primitives::Union{Nothing,Bsplines.Primitives}=nothing)`
    ... to generate a continuum orbital as Continuum.generateOrbitalForLevel(...) above, but with the combined local
        potential (nuclear + bound-electron mean field) Debye-Hueckel-screened according to the given plasma model.
        The WHOLE combined potential is screened, not just its nuclear part: screening only the (attractive) nuclear
        term while leaving the (repulsive) mean-field term unscreened leaves an uncancelled repulsive tail at large r
        once the Debye length becomes comparable to the box radius, giving an unphysical (negative) asymptotic
        effective charge for the continuum electron. For plasmaModel = Basics.NoPlasmaModel(), this delegates
        directly to the field-free method above, so that the computation of the standard (field-free) processes is
        never touched by this plasma-specific code path. A tupel of a (continuum) (orbital::Orbital, phase::Float64)
        is returned.
"""
function generateOrbitalForLevel(energy::Float64, sh::Subshell, level::Level, nm::Nuclear.Model, grid::Radial.Grid,
                                  settings::Continuum.Settings, plasmaModel::Basics.AbstractPlasmaModel;
                                  nuclearPot::Union{Nothing,Radial.Potential}=nothing,
                                  primitives::Union{Nothing,Bsplines.Primitives}=nothing)
    if  typeof(plasmaModel) == Basics.NoPlasmaModel
        return( Continuum.generateOrbitalForLevel(energy, sh, level, nm, grid, settings; nuclearPot=nuclearPot,
                                                  primitives=primitives) )
    end
    typeof(plasmaModel) == Basics.DebyeHueckelModel  ||
        error("Unsupported plasma model = $(plasmaModel)  (only Basics.DebyeHueckelModel is currently supported " *
              "for the plasma-screened continuum orbital).")
    #
    # Generate the (local) potential for the given level exactly as in the field-free case, then Debye-Hueckel
    # screen the combined (nuclear + mean-field) potential as a whole -- see the docstring note above.
    nuclearPotential = isnothing(nuclearPot) ? Nuclear.nuclearPotential(nm, grid) : nuclearPot
    wp  = Basics.computePotential(Basics.DFSField(1.0), grid, level)
    pot = Basics.add(nuclearPotential, wp)
    lambda = 1/plasmaModel.debyeLength
    screenedZr = deepcopy(pot.Zr)
    for  i = 1:length(screenedZr)   screenedZr[i] = screenedZr[i] * exp(-lambda * pot.grid.r[i])    end
    pot = Radial.Potential(pot.name * "+ Debye-Hueckel screening", screenedZr, pot.grid)

    # Generate a continuum orbital due to the given solution method
    if      Defaults.GBL_CONT_SOLUTION  ==  ContBessel()
        cOrbital = Continuum.generateOrbitalBessel(energy, sh, grid::Radial.Grid, settings)
    elseif  Defaults.GBL_CONT_SOLUTION  ==  AsymptoticCoulomb()
        cOrbital = Continuum.generateOrbitalAsymptoticCoulomb(energy, sh, pot, settings)
    elseif  Defaults.GBL_CONT_SOLUTION  ==  NonrelativisticCoulomb()
        cOrbital = Continuum.generateOrbitalNonrelativisticCoulomb(energy, sh, pot.grid, settings)
    elseif  Defaults.GBL_CONT_SOLUTION  ==  BsplineGalerkin()
        cOrbital = Continuum.generateOrbitalGalerkin(energy, sh, pot, settings; primitives=primitives)
    else    error("stop a")
    end
    #
    # Normalize the continuum orbital and determine its phase
    if      Defaults.GBL_CONT_NORMALIZATION  ==  PureSineNorm()
        cOrbital, phase = Continuum.normalizeOrbitalPureSine(cOrbital, pot.grid, settings)
    elseif  Defaults.GBL_CONT_NORMALIZATION  ==  CoulombSineNorm()
        cOrbital, phase = Continuum.normalizeOrbitalCoulombSine(cOrbital, pot, settings)
    elseif  Defaults.GBL_CONT_NORMALIZATION  ==  OngRussekNorm()
        cOrbital, phase = Continuum.normalizeOrbitalOngRussek(cOrbital, pot, settings)
    elseif  Defaults.GBL_CONT_NORMALIZATION  ==  AlokNorm()
        cOrbital, phase = Continuum.normalizeOrbitalAlok(cOrbital, pot, settings)
    else    error("stop b")
    end

    return( cOrbital, phase )
end


"""
`Continuum.generateOrbitalLocalPotential(energy::Float64, sh::Subshell, pot::Radial.Potential, settings::Continuum.Settings)`
    ... to generate a continuum orbital for the (continuum) subshell sh, the energy and the given potential. The effective 
        charge for the normalization of the continuum orbital is derived from the potential. All further specifications about 
        this generations are made by proper settings; however, the function termintates if the settings.includeExchange = true. 
        A tupel of a (continuum) (orbital::Orbital, phase::Float64, normFactor::Float64) is returned.
"""
function generateOrbitalLocalPotential(energy::Float64, sh::Subshell, pot::Radial.Potential, settings::Continuum.Settings)
    settings.includeExchange   &&   error("Continuum orbital for local potential does not allow 'exchange'.")
    
    # Generate a continuum orbital due to the given solution method
    if      Defaults.GBL_CONT_SOLUTION  ==  ContBessel() 
        cOrbital = Continuum.generateOrbitalBessel(energy, sh, pot.grid, settings)
    elseif  Defaults.GBL_CONT_SOLUTION  ==  ContSine()
        cOrbital = Continuum.generateOrbitalPureSine(energy, sh, pot.grid, settings)
    elseif  Defaults.GBL_CONT_SOLUTION  ==  AsymptoticCoulomb()
        cOrbital = Continuum.generateOrbitalAsymptoticCoulomb(energy, sh, pot, settings)
    elseif  Defaults.GBL_CONT_SOLUTION  ==  NonrelativisticCoulomb()
        cOrbital = Continuum.generateOrbitalNonrelativisticCoulomb(energy, sh, pot.grid, settings)
    elseif  Defaults.GBL_CONT_SOLUTION  ==  BsplineGalerkin()
        cOrbital = Continuum.generateOrbitalGalerkin(energy, sh, pot, settings)
    else    error("stop a")
    end
    #
    # Normalize the continuum orbital and determine its phase
    if      Defaults.GBL_CONT_NORMALIZATION  ==  PureSineNorm()
        cOrbital, phase, normFactor = Continuum.normalizeOrbitalPureSine(cOrbital, pot.grid, settings)
    elseif  Defaults.GBL_CONT_NORMALIZATION  ==  CoulombSineNorm()
        cOrbital, phase, normFactor = Continuum.normalizeOrbitalCoulombSine(cOrbital, pot, settings)
    elseif  Defaults.GBL_CONT_NORMALIZATION  ==  OngRussekNorm()
        cOrbital, phase, normFactor = Continuum.normalizeOrbitalOngRussek(cOrbital, pot, settings)
    elseif  Defaults.GBL_CONT_NORMALIZATION  ==  AlokNorm()
        cOrbital, phase             = Continuum.normalizeOrbitalAlok(cOrbital, pot, settings)
        normFactor                  = 1.0
    else    error("stop b")
    end
    #
    return( cOrbital, phase, normFactor )
end


"""
`Continuum.generateOrbitalAsymptoticCoulomb(energy::Float64, sh::Subshell, pot::Radial.Potential, settings::Continuum.Settings)`   
    ... to generate a simple relativistic and non-normalized continuum orbital but by assuming an asymptotic sin/cos behaviour 
        for both components at all mesh points. For this behaviour, both components and their derivatives can be calculated 
        analytically. The effective charge Zbar is determined for the given potential and the non-Coulombic phase shifts is set 
        to zero. An non-normalized cOrbital::Orbital is returned.
        
        Warning: At present, only the real part eta.re is taken into account; check for consistency !!
"""
function generateOrbitalAsymptoticCoulomb(energy::Float64, sh::Subshell, pot::Radial.Potential, settings::Continuum.Settings)  
    P = zeros(settings.mtp);    Q = zeros(settings.mtp);   Pprime = zeros(settings.mtp);    Qprime = zeros(settings.mtp)
    
    phi0 = 0.3
    wc   = Defaults.getDefaults("speed of light: c");    q = sqrt( 2*energy + energy/(wc*wc));    kappa = sh.kappa
    Zbar = Radial.determineZbar(pot)
    y    = Zbar * (energy + wc^2) / (wc^2 * q);   gammaBar = sqrt(kappa^2 - Zbar^2 / wc^2)
    eta  = - (kappa - im*y) / (energy + wc^2) / (gammaBar + im*y) / (2im)
    NP   = sqrt( (energy + 2*wc^2)/(pi*wc^2*q) );   NQ = -sqrt( energy / (pi*wc^2*q) )
    
    ## println("Continuum-asymptotic Coulomb orbital for q = $q,  kappa = $kappa,  Zbar=$Zbar  y = $y,   gammaBar = $gammaBar,  eta = $eta")
    if  Defaults.GBL_PRINT_DEBUG  println("Continuum-asymptotic Coulomb orbital for q = $q,  kappa = $kappa,  Zbar=$Zbar,  eta = $eta")    end

    
    for  i = 2:settings.mtp  
        theta = q * pot.grid.r[i]  +  y * log(2q * pot.grid.r[i])  - angle( SpecialFunctions.gamma(gammaBar + im*y) )  -  pi*gammaBar/2  +  eta.re
        P[i] = NP * cos(theta + phi0);      Pprime[i] = - NP * (q + y/pot.grid.r[i]) * sin(theta + phi0)       
        Q[i] = - NQ * sin(theta + phi0);    Qprime[i] = - NP * (q + y/pot.grid.r[i]) * cos(theta + phi0)       
    end

    cOrbital = Orbital( sh, false, true, energy, P, Q, Pprime, Qprime, Radial.Grid())
    
    return( cOrbital )
end


"""
`Continuum.generateOrbitalBessel(energy::Float64, sh::Subshell, grid::Radial.Grid, settings::Continuum.Settings)`   
    ... to generate a simple, non-relativistic and non-normalized Bessel wave for the large component of the continuum
        orbital and a small component due to the kinetic-balance condition. A (non-normalized) orbital::Orbital is returned.
        While Pprime is obtained from the analytical expression of the large component, Qprime is set to zero here.
"""
function generateOrbitalBessel(energy::Float64, sh::Subshell, grid::Radial.Grid, settings::Continuum.Settings)  
    P = zeros(settings.mtp);    Q = zeros(settings.mtp);   Pprime = zeros(settings.mtp);    Qprime = zeros(settings.mtp)
    
    q = sqrt( 2*energy );    l = Basics.subshell_l(sh);   wc = Defaults.getDefaults("speed of light: c")
    for  i = 2:settings.mtp   P[i] = GSL.sf_bessel_jl(l, q * grid.r[i]) * grid.r[i]    end
    for  i = 2:settings.mtp   qz = q * grid.r[i];   Pprime[i] = (l/qz * GSL.sf_bessel_jl(l, qz) - GSL.sf_bessel_jl(l+1, qz)) + P[i]/grid.r[i]     end
    for  i = 2:settings.mtp   Q[i] = -1/(2 * wc) * Pprime[i]  +  sh.kappa/grid.r[i] * P[i]         end
    
    println("Continuum-Bessel orbital for q = $q,  l = $l")

    cOrbital = Orbital( sh, false, false, energy, P, Q, Pprime, Qprime, grid)
    
    return( cOrbital )
end


"""
`Continuum.generateOrbitalGalerkin(energy::Float64, sh::Subshell, pot::Radial.Potential, settings::Continuum.Settings;
                                    primitives::Union{Nothing,Bsplines.Primitives}=nothing)`
    ... to generate a non-normalized continuum orbital within the given local potential by using the Galerkin method and a
        given B-spline basis. The B-spline primitives depend on the radial grid alone; a caller that generates many
        continuum orbitals on one and the same grid may build them once and pass them in as primitives. Omitting the
        keyword reproduces the previous behaviour exactly. A (non-normalized) orbital::Orbital is returned.
"""
function generateOrbitalGalerkin(energy::Float64, sh::Subshell, pot::Radial.Potential, settings::Continuum.Settings;
                                  primitives::Union{Nothing,Bsplines.Primitives}=nothing)
    P = zeros(settings.mtp);   Q = zeros(settings.mtp);   Pprime = zeros(settings.mtp);    Qprime = zeros(settings.mtp)
    nsL = pot.grid.nsL - 1;    nsS = pot.grid.nsS - 1
    ## The B-spline basis depends only on the grid, but this function is called once per line AND per partial wave;
    ## profiling the Auger rates alone (11-Aug-2026) put 34.3% of them into generatePrimitives.  Nothing is cached
    ## here: a caller that knows its grid is fixed builds the basis once and hands it in.
    if  isnothing(primitives)   wa = Bsplines.generatePrimitives(pot.grid)
    else                        wa = primitives
        wa.grid.nsL == pot.grid.nsL  &&  wa.grid.nsS == pot.grid.nsS  &&  wa.grid.NoPoints == pot.grid.NoPoints  ||
            error("The given primitives belong to a different grid than the potential.")
    end
    wb = Bsplines.generateGalerkinMatrix(sh, energy, pot, wa)
    wc = adjoint(wb) * wb
    
    # Test for 'real-symmetric matrix' ... and symmetrize otherwise
    nx = 0
    for  i = 1:nsL+nsS    
        for  j = i+1:nsL+nsS    
            if  abs(  (wc[i,j] - wc[j,i])/(wc[i,j] + wc[j,i]) ) > 1.0e-7   nx = nx + 1    
            end
        end
    end
    ## Only the LOWEST eigenpair is used below -- the null-space vector of wc = wb'*wb, which is the Galerkin
    ## solution.  Asking LAPACK for that one pair (syevr) instead of decomposing the whole matrix is ~2.7x
    ## faster at this size (671x671, measured 12-Aug-2026); the eigenvector agrees with the full solve to
    ## ~1e-12, and the eigenvalue -- which is ~0 by construction and enters only the printout below -- to the
    ## same absolute accuracy.  Profiling the Auger rates put 18.7% of them in this call.
    wd = Basics.diagonalize(MatrixWithLinearAlgebra(), wc; range=1:1)
    ## println(">>> Galerkin-eigenvalues = $(wd.values[1]), $(wd.values[2]) for  $sh  with  energy = $energy")
    
    cOrbital = Bsplines.generateOrbitalFromPrimitives(sh, energy, settings.mtp, wd.vectors[1], wa)  
    mtp      = size(cOrbital.P,1)
    println(">> Continuum B-spline-Galerkin orbital for energy=" * @sprintf("%.4e",energy) * ",  kappa=$(sh.kappa) " *
            "[mpt=$mtp, r[mtp]=" * @sprintf("%.4e",pot.grid.r[mtp]) * ", smallest eigenvalue=" * @sprintf("%.4e",wd.values[1]) * "].")
    
    return( cOrbital )
end

"""
`Continuum.generateOrbitalNonrelativisticCoulomb(energy::Float64, sh::Subshell, Zeff::Float64, grid::Radial.Grid, 
                                                    settings::Continuum.Settings)`   
    ... to generate a simple, non-relativistic and non-normalized Coulomb wave for the (continuum) subshell sh, 
        the small component is simply obtained by the kinetic-balance condition. A orbital::Orbital is returned.
        *****
        This function does not yet work because there is no Julia implementation of the hypergeometric function with 
        complex arguments available.
"""
function generateOrbitalNonrelativisticCoulomb(energy::Float64, sh::Subshell, Zeff::Float64, grid::Radial.Grid, settings::Continuum.Settings) 
    println("\n error: presently, NO Julia implementation of the hypergeometric function with complex arguments available.")
    
    # Define the derivative of the confluent hypergeometric function
    function hgPrime(a,b,z)    wa = a/z * (GSL.hypergeom(a+1,b,z) - GSL.hypergeom(a,b,z));   return(wa)   end
    
    P = zeros(settings.mtp);    Q = zeros(settings.mtp);   Pprime = zeros(settings.mtp);    Qprime = zeros(settings.mtp)
    
    q  = sqrt( 2*energy );    l = Basics.subshell_l(sh);   wc = Defaults.getDefaults("speed of light: c");   np = Zeff/q;   ws = 1
    wa = 2 * sqrt(Zeff) / sqrt(1 - Base.MathConstants.e^(-2pi*np)) / factorial(2l + 1)
    for  s = 1:l    ws = ws * sqrt(s^2 + np^2)    end
    
    for  i = 2:settings.mtp   P[i] = wa * ws * (2q*grid.r[i])^l * Base.MathConstants.e^(-im*q*grid.r[i]) * 
                                        GSL.hypergeom(l+1+im*np, 2l+2., 2im*q*grid.r[i])                     end
    for  i = 2:settings.mtp   
        Pprime[i] = 2q*l * (2q*grid.r[i])^(l-1) * Base.MathConstants.e^(-im*q*grid.r[i]) * GSL.hypergeom(l+1+im*np, 2l+2., 2im*q*grid.r[i]) -
                    im*q * (2q*grid.r[i])^l * Base.MathConstants.e^(-im*q*grid.r[i]) * GSL.hypergeom(l+1+im*np, 2l+2., 2im*q*grid.r[i])     +
                    (2q*grid.r[i])^l * Base.MathConstants.e^(-im*q*grid.r[i]) * hgPrime(l+1+im*np, 2l+2., 2im*q*grid.r[i])
        Pprime[i] = wa * ws * Pprime[i]          
    end
    for  i = 2:settings.mtp   Q[i] = -1/(2 * wc) * Pprime[i]  +  sh.kappa/grid.r[i] * P[i]                end

    cOrbital = Orbital( sh, false, true, energy, P, Q, Pprime, Qprime, Radial.Grid())
    
    return( cOrbital )
end


"""
+ `(energy::Float64, sh::Subshell, pot::Radial.Potential, settings::Continuum.Settings)`   
    ... to generate a simple, non-relativistic and non-normalized Coulomb wave for the (continuum) subshell sh, 
        the small component is simply obtained by the kinetic-balance condition. A orbital::Orbital is returned.
"""
function generateOrbitalNonrelativisticCoulomb(energy::Float64, sh::Subshell, pot::Radial.Potential, settings::Continuum.Settings) 
    Zbar     = Radial.determineZbar(pot)
    cOrbital = Continuum.generateOrbitalNonrelativisticCoulomb(energy, sh, Zbar, pot.grid, settings)
    
    return( cOrbital )
end


"""
`Continuum.generateOrbitalPureSine(energy::Float64, sh::Subshell, grid::Radial.Grid, settings::Continuum.Settings)`   
    ... to generate a simple, non-relativistic and non-normalized Bessel wave for the large component of the continuum
        orbital and a small component due to the kinetic-balance condition. A (non-normalized) orbital::Orbital is returned.
        While Pprime is obtained from the analytical expression of the large component, Qprime is set to zero here.
"""
function generateOrbitalPureSine(energy::Float64, sh::Subshell, grid::Radial.Grid, settings::Continuum.Settings)  
    P = zeros(settings.mtp);    Q = zeros(settings.mtp);   Pprime = zeros(settings.mtp);    Qprime = zeros(settings.mtp)
    
    phi0 = 0.0
    q = sqrt( 2*energy );    l = Basics.subshell_l(sh);   wc = Defaults.getDefaults("speed of light: c");    N = sqrt(2/(pi*q))
    for  i = 1:settings.mtp   P[i]      = N * sin(q * grid.r[i] - l*pi/2 + phi0)                  end
    for  i = 1:settings.mtp   Pprime[i] = N * q * cos(q * grid.r[i] - l*pi/2 + phi0)              end
    for  i = 2:settings.mtp   Q[i] = -1/(2 * wc) * Pprime[i]  +  sh.kappa/grid.r[i] * P[i]        end
    
    println(">> Continuum-pure sine orbital for q = $q,  l = $l")

    cOrbital = Orbital( sh, false, false, energy, P, Q, Pprime, Qprime, grid)
    
    return( cOrbital )
end


"""
`Continuum.gridConsistency(maxEnergy::Float64, grid::Radial.Grid)`   
    ... to check the consistency of the given grid with the maximum energy of the required continuum electrons; 
        an error message is issued if the grid.hp = 0.   or  15 * grid.hp < wavelength(maxEnergy)   or  if the grid has 
        less than 600 grid points. The function also returns the recommended grid point where the normalization and phase 
        is to be determined. This number is currently set to nrContinuum = grid.NoPoints - 200  ... to correct for the wrong
        'phase behaviour' at large r-values.

        OPEN POINT, stated here rather than emitted as a per-call warning as it was until 12-Aug-2026: the grid
        point at which the normalization and phase are fixed deserves a better choice than the present fixed
        offset. It is a property of the method, not of any particular run.
"""
function gridConsistency(maxEnergy::Float64, grid::Radial.Grid)
    wavenb      = sqrt( 2maxEnergy + maxEnergy * Defaults.getDefaults("alpha")^2 )
    wavelgth    = 2pi / wavenb
    nrContinuum = grid.NoPoints - 200
    
    ptsPerOsc = wavelgth / grid.hp
    #
    if      grid.hp == 0.
        error("\n\nContinuum.gridConsistency():  STOP -- this grid has grid.hp = 0. and therefore no asymptotic region " *
                "with a constant step size.\n" *
                "A continuum orbital oscillates out to the box boundary and cannot be represented on such a grid.\n" *
                ">>> Build the grid with an explicit hp, e.g.  Radial.Grid(Radial.Grid(true), rnt=4.0e-6, h=5.0e-2, " *
                "hp=2.5e-2, rbox=30.0),\n    or let Radial.generateGrid(grid; maximumFreeElectronEnergy = ...) choose hp " *
                "from the shortest de Broglie wavelength.\n")
    elseif  15 * grid.hp > wavelgth
        error("\n\nContinuum.gridConsistency():  STOP -- the radial grid is too coarse for the requested continuum energy.\n" *
                @sprintf("    continuum energy         = %.6e  Hartree\n", maxEnergy) *
                @sprintf("    de Broglie wavelength    = %.6e  a.u.\n",    wavelgth)  *
                @sprintf("    asymptotic step grid.hp  = %.6e  a.u.\n",    grid.hp)   *
                @sprintf("    points per oscillation   = %.2f      (a minimum of 15 is required)\n", ptsPerOsc) *
                "The continuum orbital would be sampled too sparsely to resolve its own oscillation. This does NOT fail " *
                "loudly by itself:\ncross sections first drift by tens of percent while still looking entirely plausible, " *
                "and then diverge by many orders of\nmagnitude once the grid can no longer follow the phase at all. " *
                "The computation is therefore stopped here rather than\nreturning numbers that cannot be trusted.\n" *
                @sprintf(">>> Reduce grid.hp to at most %.3e a.u., or lower the highest continuum energy;\n", wavelgth/15) *
                "    Radial.generateGrid(grid; maximumFreeElectronEnergy = ...) sets hp from this same wavelength.\n" )
    elseif  grid.NoPoints < 600
        error("\n\nContinuum.gridConsistency():  STOP -- too few grid points for continuum processes; " *
                "grid.NoPoints = $(grid.NoPoints) < 600.\n" *
                "The normalization and phase of the continuum orbital are determined near the outer boundary and need a " *
                "well-resolved\nasymptotic region.\n>>> Enlarge the grid (smaller h/hp, or a larger rbox).\n")
    elseif  grid.r[nrContinuum] < 2.0
        error("\n\nContinuum.gridConsistency():  STOP -- the grid does not extend far enough for continuum processes; " *
                @sprintf("grid.r[nrContinuum] = %.4e < 2.\n", grid.r[nrContinuum]) *
                "The normalization point (grid.NoPoints-200) still lies inside the atomic core region, where the orbital " *
                "has not yet\nreached its asymptotic form.\n>>> Enlarge rbox.\n")
    end
    
    return( nrContinuum )
end


"""
`Continuum.normalizeOrbitalPureSine(cOrbital::Orbital, grid::Radial.Grid, settings::Continuum.Settings)`   
    ... to normalize the given continuum orbital with regard to a (asymptotic) pure sine-function. 
        An (on-energy-scale-normalized) orbital::Orbital and its relative phase phi w.r.t.  sin(kr + phi) is returned.
"""
function normalizeOrbitalPureSine(cOrbital::Orbital, grid::Radial.Grid, settings::Continuum.Settings) 
    mtp = size( cOrbital.P, 1);    q = sqrt( 2*cOrbital.energy );        l = Basics.subshell_l(cOrbital.subshell)

    PPprime = cOrbital.P[mtp] / cOrbital.Pprime[mtp];    kr  = q * grid.r[mtp];   at = atan( PPprime * q )
    phi = at - kr + l*pi/2; 
    phi = rem(phi+1000pi, pi)    # to bring phi in the interval 0. <= phi < pi
    # For normalization, determine the maximum of the large component within the last 300 point
    A = maximum( abs.(cOrbital.P[end-300:end]) );    N   = sqrt(2/(pi*q)) / A
    
    println(">> Pure-sine normalized continuum orbital with normalization constant N=" * @sprintf("%.4e",N) *
            " and phase phi=" * @sprintf("%.4e",phi) * " at r=" * @sprintf("%.4e",grid.r[mtp])   * " a.u." )

    if abs(N) > 30.0   error("Pure-sine normalized N=" * @sprintf("%.4e",N) * "==> enlarge box-size")   end
    
    P = N .* cOrbital.P;   Q = N .* cOrbital.Q;   Pprime = N .* cOrbital.Pprime;   Qprime = N .* cOrbital.Qprime;   
    newOrbital = Orbital( cOrbital.subshell, cOrbital.isBound, cOrbital.useStandardGrid, cOrbital.energy, 
                            P, Q, Pprime, Qprime, cOrbital.grid)

    normF = maximum( abs.(P[end-300:end]) )   ## determine the maximum value of the P component as asymptotic measure for the amplitude
                            
    return( newOrbital, phi, normF )
end


"""
`Continuum.normalizeOrbitalCoulombSine(cOrbital::Orbital, pot::Radial.Potential, settings::Continuum.Settings)`   
    ... to normalize the given continuum orbital with regard to a (asymptotic) pure sine-function. 
        An (on-energy-scale-normalized) orbital::Orbital is returned.
"""
function normalizeOrbitalCoulombSine(cOrbital::Orbital, pot::Radial.Potential, settings::Continuum.Settings) 
    mtp = size( cOrbital.P, 1);   energy = cOrbital.energy    
    
    if  false  ## use non-relativistic Coulomb-sine wave representation
        wc   = Defaults.getDefaults("speed of light: c");    q = sqrt( 2*cOrbital.energy );        l = Basics.subshell_l(cOrbital.subshell)
        Zbar = Radial.determineZbar(pot);    
        PPprime = cOrbital.P[mtp] / cOrbital.Pprime[mtp];    kr  = q * pot.grid.r[mtp];   at = atan( PPprime * q )
        sigmac  = SpecialFunctions.gamma(l + 1 - im*Zbar/q)
        phi = at - kr + l*pi/2 - Zbar/q * log(2*q*pot.grid.r[mtp]) - atan(sigmac.im / sigmac.re); 
        phi = rem(phi+1000pi, pi)    # to bring phi in the interval 0. <= phi < pi
        # For normalization, determine the maximum of the large component within the last 300 point
        A = maximum( abs.(cOrbital.P[end-300:end]) );    N   = sqrt(2/(pi*q)) / A
        
        println(">> (Non-relativistic) Coulomb-sine normalized continuum orbital with normalization constant N=" * @sprintf("%.4e",N) *
                " and phase phi=" * @sprintf("%.4e",phi) * " at r=" * @sprintf("%.4e",pot.grid.r[mtp])   * " a.u." )
        #
    else      ## use relativistic Coulomb-sine wave representation
        wc   = Defaults.getDefaults("speed of light: c");    q = sqrt( 2*energy + energy/(wc*wc));      kappa = cOrbital.subshell.kappa
        Zbar = Radial.determineZbar(pot);    y = Zbar * (energy + wc^2) / (wc^2 * q);   gammaBar = sqrt(kappa^2 - Zbar^2 / wc^2)
        eta  = log( - (kappa - im*y) / (energy + wc^2) / (gammaBar + im*y)) / (2*im)
        NP   = sqrt( (energy + 2*wc^2)/(pi*wc^2*q) );   NQ = -sqrt( energy / (pi*wc^2*q) )
        
        # Select the 'ratio' on which the thetaprimephase is determined
        if  false   PprimeP = cOrbital.Pprime[mtp] / cOrbital.P[mtp];        at = atan( - PprimeP / thetaprime );   son = "P'/P"    
        else        PprimeP = cOrbital.Q[mtp] * NP / cOrbital.P[mtp] / NQ;   at = atan( PprimeP );                  son = "Q/P"     end
        #
        argG  = SpecialFunctions.gamma(gammaBar + im*y)
        theta = q * pot.grid.r[mtp]  +  y * log(2q * pot.grid.r[mtp])  - atan(argG.im / argG.re)  -  pi*gammaBar/2 ##  +  eta.re
        phi = at - theta
        phi = rem(phi+1000pi, pi)    # to bring phi in the interval 0. <= phi < pi
        ## A   = cOrbital.P[mtp] / cos(theta + phi);           N   = NP / A
        A = maximum( abs.(cOrbital.P[end-300:end]) );    N   = NP / A
        
        println(">> (Relativistic) Coulomb-sine normalized continuum orbital with normalization constant N=" * @sprintf("%.4e",N) *
                " and phase phi=" * @sprintf("%.4e",phi) * " at r=" * @sprintf("%.4e",pot.grid.r[mtp])   * " a.u." )
    end

    if abs(N) > 30.0   error("Coulomb-sine normalized N=" * @sprintf("%.4e",N) * "==> enlarge box-size")   end
    
    P = N .* cOrbital.P;   Q = N .* cOrbital.Q;   Pprime = N .* cOrbital.Pprime;   Qprime = N .* cOrbital.Qprime;   
    newOrbital = Orbital( cOrbital.subshell, cOrbital.isBound, cOrbital.useStandardGrid, cOrbital.energy, 
                            P, Q, Pprime, Qprime, cOrbital.grid)

    normF = maximum( abs.(P[end-300:end]) )   ## determine the maximum value of the P component as asymptotic measure for the amplitude
    
    return( newOrbital, phi, normF )
end


"""
`Continuum.normalizeOrbitalOngRussek(cOrbital::Orbital, pot::Radial.Potential, settings::Continuum.Settings)`   
    ... to normalize the given continuum orbital with regard to a (asymptotic) pure sine-function. 
        An (on-energy-scale-normalized) orbital::Orbital is returned.
"""
function normalizeOrbitalOngRussek(cOrbital::Orbital, pot::Radial.Potential, settings::Continuum.Settings) 
    mtp = size( cOrbital.P, 1);                         energy = cOrbital.energy;      kappa = cOrbital.subshell.kappa
    wc   = Defaults.getDefaults("speed of light: c");   E = wc^2 + cOrbital.energy;    q = sqrt( 2*energy + energy/(wc*wc))     
    qq = sqrt(E^2/wc^2 - wc^2)
    
    dVdr = pot.Zr[mtp] / (pot.grid.r[mtp]^2);    V  = - pot.Zr[mtp] / pot.grid.r[mtp]
    Pe = cOrbital.P[mtp];   Pep = cOrbital.Pprime[mtp];       r0 = pot.grid.r[mtp]
    # Compute U-function
    U = wc /2. * ((E-V)*dVdr - wc^2 *kappa*(kappa+1) /r0^3) / (( (E-V)^2 -wc^4 - wc^2 *kappa*(kappa+1) / r0^2)^(3/2)) * Pe
    ## U = U - wc * (Pep + dVdr * Pe / (2*(E-V+ 2* wc^2)) ) / (( (E-V)^2 -wc^4 - wc^2 *kappa*(kappa+1)/r0^2)^(1/2))
    U = U - wc * (Pep + dVdr * Pe / (2*(E-V+ 2* wc^2)) ) / (( (E-V)^2 -wc^4 - wc^2 *kappa*(kappa+1)/r0^2)^(1/2))
    # Compute A-function
    A = (( (E-V)^2 -wc^4 - wc^2 *kappa*(kappa+1)/r0^2)^(1/2)) / (E-V+wc^2) * (Pe^2 + U^2)
    A = sqrt(A);              phir0 = atan(U, Pe)
    phi = phir0 - q*r0;       N  = 1 / (A * sqrt(pi*wc))
    phi = rem(phi+1000pi, pi)    ## to bring phi in the interval 0. <= phi < pi
    
    println(">> WKB (Ong-Russek) normalized continuum orbital with normalization constant N=" * @sprintf("%.4e",N) *
            " and phase phi=" * @sprintf("%.4e",phi) *
            " at r=" * @sprintf("%.4e",pot.grid.r[mtp])   * " a.u.  U/Pe = $(U/Pe)" )

    if abs(N) > 100.0   error("WKB (Ong-Russek) normalized N=" * @sprintf("%.4e",N) * "==> enlarge box-size")   end
    
    P = N .* cOrbital.P;   Q = N .* cOrbital.Q;   Pprime = N .* cOrbital.Pprime;   Qprime = N .* cOrbital.Qprime;   
    newOrbital = Orbital( cOrbital.subshell, cOrbital.isBound, cOrbital.useStandardGrid, cOrbital.energy, 
                            P, Q, Pprime, Qprime, cOrbital.grid)

    normF = maximum( abs.(P[end-300:end]) )   ## determine the maximum value of the P component as asymptotic measure for the amplitude

    return( newOrbital, phi, normF )
end


"""
`Continuum.normalizeOrbitalAlok(cOrbital::Orbital, pot::Radial.Potential, settings::Continuum.Settings)`   
    ... to normalize the given continuum orbital with regard to a (asymptotic) wave function as per Salvat Code. 
        The orbitasl are normalized to unit amplitude.
        An ( orbital::Orbital, (δ + Δ)::Float64, N::Float64 ) is returned.
"""
function normalizeOrbitalAlok(cOrbital::Orbital, pot::Radial.Potential, settings::Continuum.Settings)
    mtp = size( cOrbital.P, 1)  ;             energy = cOrbital.energy;      
    kappa = cOrbital.subshell.kappa  ;        l = Basics.subshell_l(cOrbital.subshell) 
    Zbar = -Radial.determineZbar(pot)  ; 
    α = Defaults.getDefaults("alpha")  ;      wc = 1/α
    q  = sqrt( energy * (energy + 2 * wc^2) ) / wc  ;    x  = q * pot.grid.r[mtp]

    ## println("r ", pot.grid.r[mtp])

    if abs(Zbar) > 0.1
        ## println("Normalization with Coulomb functions")
        sa = "Normalization with Coulomb functions"
        λ  = sqrt(kappa^2 - Zbar^2 / wc^2)  ; λm1 = λ - 1.0
        η  = Zbar * α * (energy + wc^2) / sqrt( energy * (energy + 2 * wc^2) ) 

        xTP = η + sqrt(η^2 + λ*(λ+1.0))
        if x < xTP println("The kr is less than the Coulomb turning point") end
    
        Δ = angle(SpecialFunctions.gamma(λm1 + 1 + im * η))
        if Δ >= 0.0 Δ = mod(Δ, 2pi) else Δ = -mod(-Δ, 2pi) end
    
        θ  = x - λm1*pi/2 - η*log(2x) + Δ                               #Eq 7.3
        if θ > 1e4   θ = mod(θ, 2pi) end
    
        hgfλ = twoFzero(im*η - λm1, im*η + λm1 + 1, im*2*x)
        hgfλm1 = twoFzero(im*η - λm1 + 1, im*η + λm1 + 2, im*2*x)
        hgfλm1 = im * hgfλm1 * (im*η - λm1) * (im*η + λm1 + 1) / ( 2.0 * x^2 )
    
        GiFλm1 = hgfλ * exp(im*θ)  ;   GPiFPλm1 = ( hgfλm1 + im *(1.0 - η/x) * hgfλ ) * exp(im*θ)
        gm_1 = GiFλm1.re ; fm_1 = GiFλm1.im   ;   gpm_1 = GPiFPλm1.re ; fpm_1 = GPiFPλm1.im

        if abs(gm_1*fpm_1 - fm_1*gpm_1 - 1.0) > 1e-15 
            ef= 0.0 ; eg = 0.0
            f, fp, g, gp = GSL.sf_coulomb_wave_FG_e(η, x, λm1, 0, ef, eg)
            gm_1 = g.val ; fm_1 = f.val ; gpm_1 = gp.val ; fpm_1 = fp.val
        end
    
        f = λ * ((λ/x + η/λ)*fm_1 - fpm_1) / sqrt(λ^2 + η^2)
        g = λ * ((λ/x + η/λ)*gm_1 - gpm_1) / sqrt(λ^2 + η^2)
    
        N  =  ( α^2 * Zbar^2 * (energy + 2 * wc^2)^2 + (kappa + λ)^2 * wc^2 * q^2 )^(-0.5) / λ

        # Coulomb phase shift for λ
        rnu = angle(Zbar * α * (energy + 2 * wc^2) - im * (kappa+λ) * sqrt(energy*(energy+2*wc^2)))
        Δ =  rnu - (λ-l-1.0)*pi/2 + Δ
        if ( Zbar < 0.0 && kappa < 0)  Δ = Δ - pi ; N = -N end
        if Δ >= 0.0 Δ = mod(Δ, 2pi) else Δ = -mod(-Δ, 2pi) end
    
        fu = N * ( (kappa + λ) * sqrt(λ^2 + η^2) * wc * q * f + α * Zbar * (λ * wc^2 - kappa * (energy + wc^2)) * fm_1 )
        gu = N * ( (kappa + λ) * sqrt(λ^2 + η^2) * wc * q * g + α * Zbar * (λ * wc^2 - kappa * (energy + wc^2)) * gm_1 )
    
        fl = -N * ( α * Zbar * sqrt(λ^2 + η^2) * wc * q * f + (kappa + λ) * (λ * wc^2 - kappa * (energy + wc^2)) * fm_1 )
        gl = -N * ( α * Zbar * sqrt(λ^2 + η^2) * wc * q * g + (kappa + λ) * (λ * wc^2 - kappa * (energy + wc^2)) * gm_1 )
    
        fup = -kappa * fu / pot.grid.r[mtp] + fl * (energy + pot.Zr[mtp]/pot.grid.r[mtp] + 2.0 * wc^2 ) / wc
        gup = -kappa * gu / pot.grid.r[mtp] + gl * (energy + pot.Zr[mtp]/pot.grid.r[mtp] + 2.0 * wc^2 ) / wc

    else
        
        ## println("Normalization with Bessel functions")
        sa = "Normalization with Bessel functions"

        if kappa < 0 ksign = 1 else ksign = -1 end

        fu =  x * SpecialFunctions.sphericalbesselj(l, x)
        gu = -x * SpecialFunctions.sphericalbessely(l, x)
        fl = -sqrt(energy / ( energy + 2*wc^2 )) * ksign * x * SpecialFunctions.sphericalbesselj(l + ksign, x)
        gl =  sqrt(energy / ( energy + 2*wc^2 )) * ksign * x * SpecialFunctions.sphericalbessely(l + ksign, x)

        fup = -kappa * fu / pot.grid.r[mtp] + fl * (energy + pot.Zr[mtp]/pot.grid.r[mtp] + 2.0 * wc^2 ) / wc
        gup = -kappa * gu / pot.grid.r[mtp] + gl * (energy + pot.Zr[mtp]/pot.grid.r[mtp] + 2.0 * wc^2 ) / wc

        Δ = 0.0
    end

    p = cOrbital.P[mtp]    ;   pPrime = cOrbital.Pprime[mtp]

    δ = angle((p*gup - pPrime*gu) + im*(pPrime*fu - p*fup))
    if  ( abs(δ) > pi/2 )   δ = δ * ( 1.0 - pi/abs(δ))     end

    if  ( abs(p) > 1e-15 )  N = ( cos(δ) * fu + sin(δ) * gu ) / p else N = ( cos(δ) * fup + sin(δ) * gup ) / pPrime    end

    # Normalize wrt energy
    N = N * sqrt(q / pi / energy)

    P = N .* cOrbital.P  ;   Q = N .* cOrbital.Q  ;   Pprime = N .* cOrbital.Pprime  ;   Qprime = N .* cOrbital.Qprime  ;   
    newOrbital = Orbital( cOrbital.subshell, cOrbital.isBound, cOrbital.useStandardGrid, cOrbital.energy, 
                            P, Q, Pprime, Qprime, cOrbital.grid)
    
    println(">> $sa:   r = $(pot.grid.r[mtp]),   iPhase = $δ,   cPhase = $Δ ")

    return( newOrbital, δ + Δ, N )
end


"""
`function twoFzero(CA::ComplexF64, CB::ComplexF64, CZ::ComplexF64)`
    ... Calculates the Hypergeometric function 2F0(CA,CB;1/CZ) hypergeometric asymptotic series.
        Taken from Radial package by Salvat et al.
        A ComplexF64 value is returned.  
"""
function twoFzero(CA::ComplexF64, CB::ComplexF64, CZ::ComplexF64)

    EPS=1.0E-16; ACCUR=0.5E-15; NTERM=75

    RRP=1.0
    RRN=0.0
    RIP=0.0
    RIN=0.0
    CDF=1.0 + 0.0im
    ERR2=0.0
    ERR3=1.0
    AR=0.0
    AF=0.0
    CF=0.0 + 0.0im

    for I = 1 : NTERM
        J=I-1
        CDF=CDF*(CA+J)*(CB+J)/(I*CZ)
        ERR1=ERR2
        ERR2=ERR3
        ERR3=abs(CDF)
        if (ERR1 > ERR2 && ERR2 < ERR3) break end
        AR=CDF.re
        if(AR > 0.0) 
        RRP=RRP+AR
        else
        RRN=RRN+AR
        end
        AI=(-im*CDF).re
        if(AI > 0.0) 
        RIP=RIP+AI
        else
        RIN=RIN+AI
        end
        CF=complex(RRP+RRN,RIP+RIN)
        AF=abs(CF)
        if(AF > 1.0e25) 
        CF=0.0 + 0im
        ERR=1.0
        break
        end
        if(ERR3 < 1.0e-25*AF || ERR3 < EPS) 
        ERR=EPS
        break
        end    
    end

    # ****  Round off error.

    TR=abs(RRP+RRN)
    if(TR > 1.0e-25) 
    ERRR=(RRP-RRN)*ACCUR/TR
    else
    ERRR=1.0e0
    end
    TI=abs(RIP+RIN)
    if(TI > 1.0e-25) 
    ERRI=(RIP-RIN)*ACCUR/TI
    else
    ERRI=1.0
    end

    #  ****  ... and truncation error.
    if(AF > 1.0e-25) 
    ERR=max(ERRR,ERRI)+ERR2/AF
    else
    ERR=max(ERRR,ERRI)
    end

    return CF

end


end # module
