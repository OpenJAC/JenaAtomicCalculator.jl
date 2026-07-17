
#################################################################################################################################
### Charge exchange (CX) ########################################################################################################
##
##  Single-electron capture in a slow ion-atom or ion-ion collision,  A^q+ + B --> A^(q-1)+ + B^+,  where B carries a
##  weakly-bound "active" electron of ionization potential I_p -- B's *own*, whatever B's own charge state already
##  is. Charge exchange from a neutral-atom donor and from an already-ionized donor therefore use the *same*
##  formulas below, just with the donor's I_p taken in its actual (possibly ionized) charge state; this is why the
##  Configuration-based convenience methods accept any standard-filling Configuration, atom or ion alike, and simply
##  hand it to Empirical.ionizationPotential.
##
##  Both approximations below are restricted to the *slow-collision* (adiabatic, quasi-molecular) regime, relative
##  collision velocity v << 1 a.u. (25 keV/amu for a proton), where the active electron moves quasi-statically near
##  the saddle point of the two-center potential and the cross section becomes energy-independent to leading order.
##  Outside this regime (fast collisions) neither approximation applies; no velocity-dependent treatment (e.g.
##  Landau-Zener, which needs a coupling matrix element this module does not have access to; or CDW/AOCC) is
##  attempted here.
##  Quantity: a cross section [a.u.] -- a property of the (q, I_p) pair alone, energy-independent within the
##      slow-collision regime; fold with a relative-velocity distribution, or use
##      Empirical.chargeExchangePlasmaAlpha, to obtain a rate coefficient.


"""
`Empirical.chargeExchangeCrossSection(q::Float64, Ip::Float64, approx::Empirical.OverBarrierModel1980;
                                      printout::Bool=false, velocity::Union{Nothing,Float64}=nothing)`
    ... to estimate the total single-electron-capture cross section for a bare (or effectively bare) ion of charge
        q colliding slowly with a donor of ionization potential Ip [a.u.], by the classical (geometric) over-barrier
        model of Bohr & Lindhard, cast into cross-section form by Ryufuku & Watanabe [H. Ryufuku & T. Watanabe,
        Phys. Rev. A 19, 1538 (1979); ibid. 20, 1828 (1979)]:
            R_c = (q + 1) / (2 Ip),        sigma = pi R_c^2 ,
        with R_c the critical (saddle-point) internuclear distance at which the potential barrier between the two
        Coulomb centers drops below the donor's binding energy, allowing the active electron to move classically
        onto the (empty) projectile. A sigma::Float64 [a.u.] is returned.
        Quantity: a cross section [a.u.] -- energy-independent within the slow-collision regime (see the module
            note above); fold with a relative-velocity distribution, or pass to
            Empirical.chargeExchangePlasmaAlpha, to obtain a rate coefficient.

        Note: this is the single-active-electron, energy-independent (adiabatic) limit; its q^2 growth at high
              charge states is known to *overestimate* the cross section there -- cf. Empirical.NiehausScaling1986,
              whose q^1 growth matches measurements better for highly charged ions. Good to a factor of a few for
              low/moderate q; no state (n,l) selectivity or multi-electron capture is provided.
        Note: if velocity (relative collision velocity [a.u.]) is given, a warning is issued once it approaches or
              exceeds the adiabatic regime, v ~ 1 a.u.
"""
function chargeExchangeCrossSection(q::Float64, Ip::Float64, approx::Empirical.OverBarrierModel1980;
                                    printout::Bool=false, velocity::Union{Nothing,Float64}=nothing)
    if  q <= 0.    error("The projectile charge q = $q must be positive.")   end
    if  Ip <= 0.   error("The donor's ionization potential Ip = $Ip [a.u.] must be positive.")   end
    if  !isnothing(velocity)  &&  velocity > 0.5
        sa = "Relative collision velocity v = $velocity [a.u.] approaches or exceeds the slow-collision (adiabatic) " *
             "regime, v << 1 a.u., of the over-barrier model; the returned cross section may be unreliable."
        @warn sa maxlog=5
    end
    Rc    = (q + 1.0) / (2 * Ip)
    sigma = pi * Rc^2

    if  printout
        Empirical.chargeExchangePrintout("the classical over-barrier model (Bohr & Lindhard; Ryufuku & Watanabe " *
            "1979/1980), valid for slow (adiabatic) collisions, v << 1 a.u.; R_c = (q+1)/(2 Ip), sigma = pi R_c^2; " *
            "good to a factor of a few, overestimates at high q", q, Ip, velocity, sigma)
    end

    return( sigma )
end


"""
`Empirical.chargeExchangeCrossSection(q::Float64, targetConf::Configuration, targetZ::Float64,
                                      approx::Empirical.OverBarrierModel1980; printout::Bool=false,
                                      velocity::Union{Nothing,Float64}=nothing)`
    ... convenience wrapper that determines the donor's ionization potential from targetConf via
        Empirical.ionizationPotential(targetZ, targetConf) -- valid for any standard-filling configuration, neutral
        atom or ion alike, so the same call handles both atom-ion and ion-ion charge exchange -- and dispatches to
        the direct (q, Ip) method above. A sigma::Float64 [a.u.] is returned.
"""
function chargeExchangeCrossSection(q::Float64, targetConf::Configuration, targetZ::Float64,
                                    approx::Empirical.OverBarrierModel1980; printout::Bool=false,
                                    velocity::Union{Nothing,Float64}=nothing)
    Ip = Empirical.ionizationPotential(round(Int64, targetZ), targetConf)
    return( Empirical.chargeExchangeCrossSection(q, Ip, approx; printout=printout, velocity=velocity) )
end


"""
`Empirical.chargeExchangeCrossSection(q::Float64, Ip::Float64, approx::Empirical.NiehausScaling1986;
                                      printout::Bool=false, velocity::Union{Nothing,Float64}=nothing)`
    ... to estimate the total single-electron-capture cross section for a bare (or effectively bare) ion of charge
        q colliding slowly with a donor of ionization potential Ip [a.u.], by the empirical high-charge scaling law
            sigma = 2.6e-13 * q / Ip[eV]^2   [cm^2] ,
        established experimentally for slow (v < 1 a.u.) collisions of highly charged ions with atoms [e.g.
        M. Kimura et al., J. Phys. B 28, L643 (1995); K. Hosaka et al., Fus. Eng. Design 34-35, 781 (1997)], and
        shown to be the high-q limit of an extended classical over-barrier model [A. Niehaus, J. Phys. B 19, 2925
        (1986)]. A sigma::Float64 [a.u.] is returned.
        Quantity: a cross section [a.u.] -- energy-independent within the slow-collision regime (see the module
            note above); fold with a relative-velocity distribution, or pass to
            Empirical.chargeExchangePlasmaAlpha, to obtain a rate coefficient.

        Note: this scaling grows only linearly in q, matching experiment better than the plain (q^2-growing)
              Empirical.OverBarrierModel1980 at high charge states; it is a global fit and not expected to be
              accurate for individual, near-resonant systems (where the specific energy defect of the reaction
              matters), nor does it provide state (n,l) selectivity or multi-electron capture.
        Note: if velocity (relative collision velocity [a.u.]) is given, a warning is issued once it approaches or
              exceeds the adiabatic regime, v ~ 1 a.u.
"""
function chargeExchangeCrossSection(q::Float64, Ip::Float64, approx::Empirical.NiehausScaling1986;
                                    printout::Bool=false, velocity::Union{Nothing,Float64}=nothing)
    if  q <= 0.    error("The projectile charge q = $q must be positive.")   end
    if  Ip <= 0.   error("The donor's ionization potential Ip = $Ip [a.u.] must be positive.")   end
    if  !isnothing(velocity)  &&  velocity > 0.5
        sa = "Relative collision velocity v = $velocity [a.u.] approaches or exceeds the slow-collision (adiabatic) " *
             "regime, v << 1 a.u., of this empirical scaling; the returned cross section may be unreliable."
        @warn sa maxlog=5
    end
    ## Constant C = 2.6e-13 cm^2 eV^2, expressed in atomic units [a_o^2 Hartree^2]; cf. the identical unit-conversion
    ## recipe used for Lotz's constant in Empirical.impactIonizationCrossSection.
    Cau   = 2.6e-13 / Defaults.convertUnits("length: from atomic to cm", 1.0)^2 *
            Defaults.convertUnits("energy: from eV to atomic", 1.0)^2
    sigma = Cau * q / Ip^2

    if  printout
        Empirical.chargeExchangePrintout("the empirical high-charge scaling law sigma = 2.6e-13 q/Ip[eV]^2 cm^2 " *
            "(Kimura et al. 1995; Hosaka et al. 1997), the high-q limit of Niehaus' (1986) extended over-barrier " *
            "model; grows linearly in q, matching experiment better than OverBarrierModel1980 at high q", q, Ip, velocity, sigma)
    end

    return( sigma )
end


"""
`Empirical.chargeExchangeCrossSection(q::Float64, targetConf::Configuration, targetZ::Float64,
                                      approx::Empirical.NiehausScaling1986; printout::Bool=false,
                                      velocity::Union{Nothing,Float64}=nothing)`
    ... convenience wrapper that determines the donor's ionization potential from targetConf via
        Empirical.ionizationPotential(targetZ, targetConf) -- valid for any standard-filling configuration, neutral
        atom or ion alike, so the same call handles both atom-ion and ion-ion charge exchange -- and dispatches to
        the direct (q, Ip) method above. A sigma::Float64 [a.u.] is returned.
"""
function chargeExchangeCrossSection(q::Float64, targetConf::Configuration, targetZ::Float64,
                                    approx::Empirical.NiehausScaling1986; printout::Bool=false,
                                    velocity::Union{Nothing,Float64}=nothing)
    Ip = Empirical.ionizationPotential(round(Int64, targetZ), targetConf)
    return( Empirical.chargeExchangeCrossSection(q, Ip, approx; printout=printout, velocity=velocity) )
end


"""
`Empirical.chargeExchangePrintout(sa::String, q::Float64, Ip::Float64, velocity::Union{Nothing,Float64}, sigma::Float64)`
    ... to print a common report about a charge-exchange cross-section estimate: the approximation and its regime
        of validity (sa), the projectile charge and donor ionization potential, the relative velocity if given, and
        the resulting cross section in the user-defined unit. Nothing is returned.
"""
function chargeExchangePrintout(sa::String, q::Float64, Ip::Float64, velocity::Union{Nothing,Float64}, sigma::Float64)
    unCs   = Defaults.getDefaults("unit: cross section");   unEnergy = Defaults.getDefaults("unit: energy")
    Ipx    = Defaults.convertUnits("energy: from atomic to " * unEnergy, Ip)
    sigx   = Defaults.convertUnits("cross section: from atomic to " * unCs, sigma)
    sb = "\n* Estimate empirically the single-electron-capture cross section for a slow ion-atom/ion-ion collision " *
         "with the following assumptions/simplifications: " *
         "\n    + Use " * sa * ". " *
         "\n    + Projectile charge q = $q;  donor ionization potential Ip [$unEnergy] = $Ipx " *
         (isnothing(velocity) ? "" : "\n    + Relative collision velocity [a.u.] = $velocity ") *
         "\n    + Cross section [$unCs] = $sigx " *
         "\n    + Quantity: cross section [$unCs] -- a property of the (q, Ip) pair alone, energy-independent within the slow-collision regime; fold with a relative-velocity distribution, or pass to Empirical.chargeExchangePlasmaAlpha, to obtain a rate. " * "\n"
    println(sb)

    return( nothing )
end


"""
`Empirical.chargeExchangePlasmaAlpha(T::Float64, q::Float64, Ip::Float64, Mproj::Float64, Mtarget::Float64,
                                     approx::Empirical.AbstractEmpiricalApproximation; printout::Bool=false)`
    ... to estimate the charge-exchange plasma rate coefficient alpha^(CX) for a Maxwellian relative-velocity
        distribution at temperature T [a.u.], by folding the (velocity-independent, within the slow-collision
        regime) cross section with the mean relative speed of a Maxwell-Boltzmann distribution,
            alpha^(CX) (T) = <v> sigma = sqrt(8 T / (pi mu)) * sigma ,
        with mu = Mproj Mtarget / (Mproj + Mtarget) [m_e] the reduced mass of the collision partners (in the same
        mass unit as, e.g., Empirical.IonProjectile: for a proton, Defaults.PROTON_MASS_U/Defaults.ELECTRON_MASS_U).
        Unlike the electron-impact plasma rate coefficients elsewhere in this module, no Gauss-Legendre quadrature
        is needed here: sigma is constant within the slow-collision regime (see the module note above), so
        <v sigma> = <v> sigma analytically. approx selects the cross-section formula (OverBarrierModel1980 or
        NiehausScaling1986) via Empirical.chargeExchangeCrossSection; any other approx raises the same
        informative MethodError as an unsupported cross-section call. An alpha::Float64 [a.u.] is returned.
        Quantity: a rate coefficient [cm^3/s] -- multiply by the target (donor) number density n_B [1/cm^3] to
            obtain the rate per (projectile) ion [1/s].

        Note: a warning is issued if the thermal relative speed sqrt(2 T/mu) approaches or exceeds the adiabatic
              regime (v ~< 1 a.u.), since the underlying cross section becomes unreliable for a growing fraction
              of the distribution.
"""
function chargeExchangePlasmaAlpha(T::Float64, q::Float64, Ip::Float64, Mproj::Float64, Mtarget::Float64,
                                   approx::Empirical.AbstractEmpiricalApproximation; printout::Bool=false)
    mu  = Mproj * Mtarget / (Mproj + Mtarget)
    vth = sqrt(2*T / mu)
    sigma = Empirical.chargeExchangeCrossSection(q, Ip, approx; printout=false, velocity=vth)
    vmean = sqrt(8*T / (pi*mu))
    alpha = vmean * sigma

    if  printout
        unEnergy = Defaults.getDefaults("unit: energy");   unCs = Defaults.getDefaults("unit: cross section")
        Tx     = Defaults.convertUnits("energy: from atomic to " * unEnergy, T)
        Ipx    = Defaults.convertUnits("energy: from atomic to " * unEnergy, Ip)
        sigx   = Defaults.convertUnits("cross section: from atomic to " * unCs, sigma)
        factor = Defaults.convertUnits("length: from atomic to cm", 1.0)^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)
        alphax = factor * alpha
        sb = "\n* Estimate empirically the charge-exchange plasma rate coefficient alpha for a Maxwellian relative-" *
             "velocity distribution with the following assumptions/simplifications: " *
             "\n    + Fold a velocity-independent cross section with the Maxwellian mean relative speed " *
             "<v> = sqrt(8 T/(pi mu)); no quadrature is needed since sigma is constant within the slow-collision regime. " *
             "\n    + Temperature T [$unEnergy] = $Tx;  reduced mass mu [m_e] = $mu " *
             "\n    + Projectile charge q = $q;  donor ionization potential Ip [$unEnergy] = $Ipx " *
             "\n    + Cross section [$unCs] = $sigx " *
             "\n    + Plasma rate coefficient alpha^(CX) [cm^3/s] = $alphax " *
             "\n    + Quantity: a rate coefficient [cm^3/s] -- multiply by the target number density n_B [1/cm^3] for the rate per (projectile) ion [1/s]. " * "\n"
        println(sb)
    end

    return( alpha )
end


"""
`Empirical.chargeExchangePlasmaAlpha(T::Float64, q::Float64, targetConf::Configuration, targetZ::Float64,
                                     Mproj::Float64, Mtarget::Float64,
                                     approx::Empirical.AbstractEmpiricalApproximation; printout::Bool=false)`
    ... convenience wrapper that determines the donor's ionization potential from targetConf via
        Empirical.ionizationPotential(targetZ, targetConf) -- valid for any standard-filling configuration, neutral
        atom or ion alike -- and dispatches to the direct (T, q, Ip, Mproj, Mtarget) method above. An alpha::Float64
        [a.u.] is returned.
"""
function chargeExchangePlasmaAlpha(T::Float64, q::Float64, targetConf::Configuration, targetZ::Float64,
                                   Mproj::Float64, Mtarget::Float64,
                                   approx::Empirical.AbstractEmpiricalApproximation; printout::Bool=false)
    Ip = Empirical.ionizationPotential(round(Int64, targetZ), targetConf)
    return( Empirical.chargeExchangePlasmaAlpha(T, q, Ip, Mproj, Mtarget, approx; printout=printout) )
end


#################################################################################################################################
### State-selective (n,l) charge exchange #######################################################################################
##
##  Empirical.chargeExchangeCrossSection() above returns only the *total*, state-summed capture cross section. The
##  functions below additionally estimate *which* shell the electron is most likely captured into -- a considerably
##  cruder, explicitly labeled extension (see the "StateSelective" name and the named-tuple field names below, so
##  the distinction from the total cross section is unambiguous from the call site alone).


"""
`Empirical.chargeExchangeCaptureShell(q::Float64, Ip::Float64)`
    ... to estimate the principal quantum number of the shell into which the active electron is most likely
        captured in a slow-collision charge-exchange event, by the standard hydrogenic level-matching argument
        associated with the over-barrier model: the electron is assumed to be captured near the point where the
        H-like binding energy of shell n in the acceptor ion of charge q, E_n = q^2/(2 n^2), matches the donor's
        own ionization potential Ip [a.u.], i.e.
            n_c = q / sqrt(2 Ip) .
        This level-matching argument is commonly invoked alongside the over-barrier model [e.g. H. Ryufuku,
        K. Sasaki & T. Watanabe, Phys. Rev. A 21, 745 (1980); A. Niehaus, J. Phys. B 19, 2925 (1986)] to identify
        the dominant shell of capture; it is a direct algebraic consequence of the hydrogenic formula already used
        throughout this module (cf. Empirical.scaledBindingEnergy), not an independently fitted quantity. A named
        tuple (nc::Float64=, n::Int64=) is returned: nc is the exact (generally non-integer) level-matching value,
        and n = max(1, round(Int64, nc)) is the nearest physical shell.
        Quantity: n and nc are dimensionless (quantum) labels, not an energy or a cross section.

        Note: this identifies only the single *dominant* shell; the true capture probability is known to spread
              over a range of a few adjacent n (particularly when nc falls near a half-integer), which this
              estimate does not model. Empirical.chargeExchangeStateSelectiveCrossSection() assumes *all* capture
              goes into this one shell -- a further, explicit simplification beyond the total-cross-section
              estimate itself.
"""
function chargeExchangeCaptureShell(q::Float64, Ip::Float64)
    if  q <= 0.    error("The projectile charge q = $q must be positive.")   end
    if  Ip <= 0.   error("The donor's ionization potential Ip = $Ip [a.u.] must be positive.")   end
    nc = q / sqrt(2*Ip)
    n  = max(1, round(Int64, nc))

    return( (nc = nc, n = n) )
end


"""
`Empirical.chargeExchangeStateSelectiveCrossSection(q::Float64, Ip::Float64,
                                                     approx::Empirical.AbstractEmpiricalApproximation;
                                                     printout::Bool=false)`
    ... to split the *total* single-electron-capture cross section (Empirical.chargeExchangeCrossSection) into a
        state-selective (n,l) estimate: all capture is assumed to proceed into the single dominant shell
        n = Empirical.chargeExchangeCaptureShell(q, Ip).n, distributed over its orbital substates l = 0, ..., n-1
        by the purely statistical hydrogenic degeneracy weight,
            sigma_l = sigma_total * (2l+1) / n^2 ,
        with no dynamical (angular-momentum-transfer) preference among them -- the simplest defensible assumption
        absent a detailed (e.g. Landau-Zener or close-coupling) calculation. approx selects the total-cross-section
        formula (OverBarrierModel1980 or NiehausScaling1986) via Empirical.chargeExchangeCrossSection; any other
        approx raises the same informative MethodError as an unsupported cross-section call. A named tuple
        (n::Int64=, nc::Float64=, sigmaTotal::Float64=, states::Array=) is returned, where states is an array of
        named tuples (l::Int64=, sigma::Float64=, fraction::Float64=) with sum(fraction) == 1 and
        sum(sigma) == sigmaTotal by construction.
        Quantity: sigmaTotal and each sigma [a.u.] are cross sections, energy-independent within the slow-collision
            regime, exactly as Empirical.chargeExchangeCrossSection; fraction is dimensionless.

        Note: the function name is deliberately explicit -- "StateSelective" -- to distinguish it from
              Empirical.chargeExchangeCrossSection, which returns only the total, state-summed cross section. The
              n,l-selectivity here rests on the single-dominant-shell and statistical-l simplifications described
              above and in Empirical.chargeExchangeCaptureShell; it is considerably cruder than the total cross
              section itself and should be read as indicative of the capture-shell order and a plausible
              l-spread, not as a precise state-resolved prediction.
"""
function chargeExchangeStateSelectiveCrossSection(q::Float64, Ip::Float64,
                                                   approx::Empirical.AbstractEmpiricalApproximation;
                                                   printout::Bool=false)
    (nc, n)    = Empirical.chargeExchangeCaptureShell(q, Ip)
    sigmaTotal = Empirical.chargeExchangeCrossSection(q, Ip, approx; printout=false)
    states     = [ (l = l, sigma = sigmaTotal * (2l+1)/n^2, fraction = (2l+1)/n^2)   for l = 0:n-1 ]

    if  printout
        unCs = Defaults.getDefaults("unit: cross section");   unEnergy = Defaults.getDefaults("unit: energy")
        Ipx   = Defaults.convertUnits("energy: from atomic to " * unEnergy, Ip)
        sigTx = Defaults.convertUnits("cross section: from atomic to " * unCs, sigmaTotal)
        sb = "\n* Estimate empirically the state-selective (n,l) single-electron-capture cross section for a slow " *
             "ion-atom/ion-ion collision, by splitting the total cross section under the following " *
             "assumptions/simplifications: " *
             "\n    + Total cross section from $approx (cf. Empirical.chargeExchangeCrossSection). " *
             "\n    + ALL capture assumed into the single dominant shell n = $n (level-matching estimate nc = $(round(nc, digits=3))). " *
             "\n    + Statistical (2l+1)/n^2 weighting over l = 0, ..., $(n-1) within that shell; no dynamical " *
             "l-preference is modeled. " *
             "\n    + Projectile charge q = $q;  donor ionization potential Ip [$unEnergy] = $Ipx " *
             "\n    + Total cross section [$unCs] = $sigTx " *
             "\n    + l-resolved cross sections [$unCs]:"
        for  st in states
            sigx = Defaults.convertUnits("cross section: from atomic to " * unCs, st.sigma)
            sb = sb * "\n        l = $(st.l):  fraction = $(round(st.fraction, digits=4)),  sigma [$unCs] = $sigx"
        end
        sb = sb * "\n    + Quantity: sigma [$unCs] per (n,l) state -- a considerably cruder estimate than the " *
             "total cross section; read as indicative of the capture-shell order and l-spread only. " * "\n"
        println(sb)
    end

    return( (n = n, nc = nc, sigmaTotal = sigmaTotal, states = states) )
end


"""
`Empirical.chargeExchangeStateSelectiveCrossSection(q::Float64, targetConf::Configuration, targetZ::Float64,
                                                     approx::Empirical.AbstractEmpiricalApproximation;
                                                     printout::Bool=false)`
    ... convenience wrapper that determines the donor's ionization potential from targetConf via
        Empirical.ionizationPotential(targetZ, targetConf) -- valid for any standard-filling configuration, neutral
        atom or ion alike -- and dispatches to the direct (q, Ip, approx) method above. A named tuple
        (n::Int64=, nc::Float64=, sigmaTotal::Float64=, states::Array=) is returned.
"""
function chargeExchangeStateSelectiveCrossSection(q::Float64, targetConf::Configuration, targetZ::Float64,
                                                   approx::Empirical.AbstractEmpiricalApproximation;
                                                   printout::Bool=false)
    Ip = Empirical.ionizationPotential(round(Int64, targetZ), targetConf)
    return( Empirical.chargeExchangeStateSelectiveCrossSection(q, Ip, approx; printout=printout) )
end
