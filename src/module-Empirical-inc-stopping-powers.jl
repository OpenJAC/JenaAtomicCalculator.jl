
#################################################################################################################################
### Stopping powers (SP) ########################################################################################################
##
##  Energy loss of a fast electron to the free (thermal) electrons of a plasma. All approximations below share the
##  form  - dE/dx = n_e (2 pi e^4 / E) ln Lambda  and differ only in the Coulomb logarithm ln Lambda; the plasma
##  frequency omega_p = sqrt(4 pi n_e e^2/m_e) follows from the electron density and is derived internally.
##  Quantity: a stopping power [Hartree/a_o] -- the energy loss per unit path length; it already contains the electron
##      density n_e. Divide by n_e for the stopping cross section [Hartree a_o^2], or multiply by the electron velocity
##      v(eps) for the energy-loss rate -dE/dt [Hartree per atomic time].


"""
`Empirical.stoppingPower(energy::Float64, ne::Float64;
                         approx::Empirical.AbstractEmpiricalApproximation=KozmaFranson1992(), printout::Bool=false)`
    ... to estimate empirically the stopping power - dE/dx of an electron with kinetic energy [a.u.] due to the free
        electrons of a plasma with electron density ne [a.u.], by applying some simple approximation as determined
        by approx. A sp::Float64 [a.u.] is returned.
"""
function stoppingPower(energy::Float64, ne::Float64;
                       approx::Empirical.AbstractEmpiricalApproximation=KozmaFranson1992(), printout::Bool=false)

    sp = Empirical.stoppingPower([energy], ne, approx, printout=printout)

    return( sp[1] )
end


"""
`Empirical.stoppingPower(energies::Array{Float64,1}, ne::Float64, approx::Empirical.Bohr1913; printout::Bool=false)`
    ... to estimate the stopping power - dE/dx of an electron at the given kinetic energies [a.u.] due to the free
        electrons of a plasma with electron density ne [a.u.], by using Bohr's (1913) classical Coulomb logarithm,
            - dE/dx = n_e (2 pi e^4/E) ln( 1.123 m_e v^3 / (e^2 omega_p) ),
        where the constant 1.123 = 2 exp(-gamma_E) carries Euler's constant gamma_E = 0.5772. The logarithm is
        clamped at zero: for velocities below the validity range no (negative) stopping power is returned.
        A sps::Array{Float64,1} [a.u.] is returned.
        Quantity: a stopping power [Hartree/a_o] -- the energy loss per unit path length; it already contains the
            electron density n_e. Divide by n_e for the stopping cross section [Hartree a_o^2].

        Note: The classical logarithm applies for velocities below v = 1 a.u. (E < 13.6 eV), where the classical
              distance of closest approach exceeds the de Broglie wavelength, but still well above the thermal
              velocity of the plasma electrons. It refers to a heavy (spinless, distinguishable) projectile; the
              electron-projectile variant of Kozma & Fransson (1992) is smaller by ln 2 inside the logarithm.
"""
function stoppingPower(energies::Array{Float64,1}, ne::Float64, approx::Empirical.Bohr1913; printout::Bool=false)
    omegaP = sqrt(4pi * ne)
    sps    = Float64[]
    for  eps in energies
        v = sqrt(2*eps)
        push!(sps, ne * 2pi / eps * max(0., log(1.123 * v^3 / omegaP)) )
    end

    if  printout   Empirical.stoppingPowerPrintout("Bohr's (1913) classical Coulomb logarithm ln(1.123 m v^3/(e^2 omega_p))",
                                                   energies, ne, omegaP, sps)      end
    return( sps )
end


"""
`Empirical.stoppingPower(energies::Array{Float64,1}, ne::Float64, approx::Empirical.Bethe1931; printout::Bool=false)`
    ... to estimate the stopping power - dE/dx of an electron at the given kinetic energies [a.u.] due to the free
        electrons of a plasma with electron density ne [a.u.], by using the quantal (Bethe-type) Coulomb logarithm
        for an electron projectile,
            - dE/dx = n_e (2 pi e^4/E) ln( 2 E / (hbar omega_p) ),
        cf. Kozma & Fransson (1992, Eq. 1). The logarithm is clamped at zero. A sps::Array{Float64,1} [a.u.]
        is returned.
        Quantity: a stopping power [Hartree/a_o] -- the energy loss per unit path length; it already contains the
            electron density n_e. Divide by n_e for the stopping cross section [Hartree a_o^2].

        Note: The quantal logarithm applies for velocities above v = 1 a.u. (E > 13.6 eV), where the de Broglie
              wavelength exceeds the classical distance of closest approach. Bethe's original (heavy-projectile)
              logarithm ln(2 m v^2/(hbar omega_p)) is larger by ln 2; for an electron projectile the maximum energy
              transfer is E/2, which removes this factor.
"""
function stoppingPower(energies::Array{Float64,1}, ne::Float64, approx::Empirical.Bethe1931; printout::Bool=false)
    omegaP = sqrt(4pi * ne)
    sps    = Float64[]
    for  eps in energies
        push!(sps, ne * 2pi / eps * max(0., log(2*eps / omegaP)) )
    end

    if  printout   Empirical.stoppingPowerPrintout("the quantal (Bethe-type) Coulomb logarithm ln(2 E/(hbar omega_p)) " *
                                                   "for an electron projectile", energies, ne, omegaP, sps)      end
    return( sps )
end


"""
`Empirical.stoppingPower(energies::Array{Float64,1}, ne::Float64, approx::Empirical.KozmaFranson1992; printout::Bool=false)`
    ... to estimate the stopping power - dE/dx of an electron at the given kinetic energies [a.u.] due to the free
        electrons of a plasma with electron density ne [a.u.], by using the piecewise electron loss function of
        Kozma & Fransson (1992, Eqs. 1-3),
            - dE/dx = n_e (2 pi e^4/E) ln( 2 E / (hbar omega_p) )                for  E > 14 eV,
            - dE/dx = n_e (2 pi e^4/E) ln( m_e v^3 / (gamma e^2 omega_p) )       for  kT << E < 14 eV,
        with gamma = exp(gamma_E) = 1.7811 (gamma_E = 0.5772, Euler's constant; Schunk & Hays 1971). The boundary
        at 14 eV is physical: it marks E = m e^4/(2 hbar^2) = 1 Ry, where the de Broglie wavelength overtakes the
        classical distance of closest approach. The logarithm is clamped at zero. A sps::Array{Float64,1} [a.u.]
        is returned.
        Quantity: a stopping power [Hartree/a_o] -- the energy loss per unit path length; it already contains the
            electron density n_e. Divide by n_e for the stopping cross section [Hartree a_o^2].

        Note: Both branches refer to an electron projectile and are smaller by ln 2 inside the logarithm than their
              heavy-projectile (Bohr/Bethe) counterparts. At the 14 eV boundary the prescription jumps by
              -gamma_E = -0.577 inside the logarithm (a few percent of ln Lambda), as in the original.
"""
function stoppingPower(energies::Array{Float64,1}, ne::Float64, approx::Empirical.KozmaFranson1992; printout::Bool=false)
    omegaP = sqrt(4pi * ne)
    e14    = Defaults.convertUnits("energy: from eV to atomic", 14.0)
    sps    = Float64[]
    for  eps in energies
        if  eps > e14   push!(sps, ne * 2pi / eps * max(0., log(2*eps / omegaP)) )
        else            v = sqrt(2*eps)
                        push!(sps, ne * 2pi / eps * max(0., log(v^3 / (1.7811 * omegaP))) )
        end
    end

    if  printout   Empirical.stoppingPowerPrintout("the Kozma & Fransson (1992) loss function: quantal above and " *
                                                   "classical below E = 14 eV", energies, ne, omegaP, sps)      end
    return( sps )
end


"""
`Empirical.stoppingPower(energies::Array{Float64,1}, ne::Float64, approx::Empirical.Axelrod1980; printout::Bool=false)`
    ... to estimate the stopping power - dE/dx of an electron at the given kinetic energies [a.u.] due to the free
        electrons of a plasma with electron density ne [a.u.], by using Axelrod's (1980) relativistic plasma energy
        loss, i.e. the Bethe-type loss with the plasma energy hbar omega_p as effective ionization potential,
            - dE/dx = n_e (4 pi e^4/(m_e c^2 beta^2)) [ ln( sqrt(gamma-1) gamma beta m_e c^2 / (hbar omega_p) )
                          + 1/2 ln 2 - (beta^2/12) (23/2 + 7/(gamma+1) + 5/(gamma+1)^2 + 2/(gamma+1)^3) ],
        with the Lorentz factor gamma = 1 + E/(m_e c^2); cf. Milne et al. (1999, Eqs. 1 and 3). The bracket is
        clamped at zero. A sps::Array{Float64,1} [a.u.] is returned.
        Quantity: a stopping power [Hartree/a_o] -- the energy loss per unit path length; it already contains the
            electron density n_e. Divide by n_e for the stopping cross section [Hartree a_o^2].

        Note: In the nonrelativistic limit the bracket reduces exactly to the quantal logarithm ln(2 E/(hbar omega_p))
              of Bethe1931/KozmaFranson1992, since sqrt(gamma-1) gamma beta m_e c^2 -> sqrt(2) E and the 1/2 ln 2
              completes the factor 2. The (beta^2/12) term is the (positron-type) relativistic correction as
              transcribed by Milne et al.; it matters only for E >~ 100 keV.
"""
function stoppingPower(energies::Array{Float64,1}, ne::Float64, approx::Empirical.Axelrod1980; printout::Bool=false)
    omegaP = sqrt(4pi * ne);    c = Defaults.getDefaults("speed of light: c")
    sps    = Float64[]
    for  eps in energies
        gam  = 1.0 + eps / c^2
        bet2 = 1.0 - 1.0 / gam^2
        wa   = log( sqrt(gam - 1.0) * gam * sqrt(bet2) * c^2 / omegaP ) + log(2.0)/2   -
               bet2/12 * ( 23/2 + 7/(gam + 1.0) + 5/(gam + 1.0)^2 + 2/(gam + 1.0)^3 )
        push!(sps, ne * 4pi / (c^2 * bet2) * max(0., wa) )
    end

    if  printout   Empirical.stoppingPowerPrintout("Axelrod's (1980) relativistic plasma energy loss with " *
                                                   "hbar omega_p as effective ionization potential", energies, ne, omegaP, sps)   end
    return( sps )
end


"""
`Empirical.stoppingPowerPrintout(sa::String, energies::Array{Float64,1}, ne::Float64, omegaP::Float64,
                                 sps::Array{Float64,1})`
    ... to print a common report about the given stopping powers: the approximation sa, the electron density and
        plasma frequency, the energies in user-defined units as well as the stopping powers in eV/cm. Nothing is
        returned.
"""
function stoppingPowerPrintout(sa::String, energies::Array{Float64,1}, ne::Float64, omegaP::Float64,
                               sps::Array{Float64,1})
    unEnergy = Defaults.getDefaults("unit: energy")
    aoCm     = Defaults.convertUnits("length: from atomic to cm", 1.0)
    neCm3    = ne / aoCm^3
    omegaPx  = Defaults.convertUnits("energy: from atomic to " * unEnergy, omegaP)
    energiex = [Defaults.convertUnits("energy: from atomic to " * unEnergy, eps)   for eps in energies]
    spsx     = [Defaults.convertUnits("energy: from atomic to eV", sp) / aoCm      for sp  in sps]
    sb = "\n* Estimate empirically the stopping power - dE/dx of an electron due to the free electrons of a plasma " *
         "with the following assumptions/simplifications: " *
         "\n    + Use " * sa * ". " *
         "\n    + Electron density n_e [1/cm^3]  = $neCm3;   plasma energy hbar omega_p [$unEnergy] = $omegaPx " *
         "\n    + Electron energies [$unEnergy]  = $energiex " *
         "\n    + Stopping powers - dE/dx [eV/cm] = $spsx " *
         "\n    + Quantity: a stopping power [eV/cm] -- the energy loss per unit path length; it already contains the electron density n_e. Divide by n_e for the stopping cross section [eV cm^2]. " * "\n"
    println(sb)

    return( nothing )
end
