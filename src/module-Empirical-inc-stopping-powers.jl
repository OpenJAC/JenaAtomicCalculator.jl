
#################################################################################################################################
### Stopping powers (SP) ########################################################################################################
#
#  Energy loss - dE/dx of a projectile (electron, positron or bare ion) to the electrons of a material: the free
#  electrons of a plasma, the bound electrons of a neutral atom gas, or both in a partially ionized gas. All
#  approximations below share the (Bohr/Bethe-type) form
#
#      - dE/dx  =  sum_t  n_t  (4 pi z^2 e^4)/(m_e v^2)  ln Lambda_t ,
#
#  where the sum runs over the electron populations t of the material with densities n_t, v is the projectile
#  velocity, and the Coulomb logarithm ln Lambda_t = ln(b_max/b_min) counts the effective range of impact
#  parameters; each population enters through its density and one characteristic energy (hbar omega_p for free,
#  Ibar for bound electrons). The projectile enters through z^2, its mass (v = sqrt(2E/M)) and its statistics.
#  Quantity: a stopping power [Hartree/a_o] -- the energy loss per unit path length; it already contains the
#      electron densities of the material. Divide by the density for the stopping cross section, or by the mass
#      density for the mass stopping power [MeV cm^2/g].


"""
`Empirical.meanExcitationEnergy(Z::Int64)`
    ... to estimate the mean excitation energy Ibar of the bound electrons of a neutral atom with nuclear charge Z,
            Ibar = 9.1 Z (1 + 1.9 Z^(-2/3))  [eV] ,
        following Segre (1977) and Roy & Reed (1968), as quoted by Milne et al. (1999, Eq. 2). This estimate agrees
        with the ICRU tabulations to typically 5-10% for medium and heavy elements (Al: 159 eV vs. the tabulated
        166 eV) but becomes rough for the lightest atoms (H: 26 eV vs. the tabulated 19.2 eV); since Ibar only
        enters logarithmically, this affects the stopping power at the few-percent level. An energy::Float64 [a.u.]
        is returned.
"""
function meanExcitationEnergy(Z::Int64)
    wa = 9.1 * Z * (1.0 + 1.9 * Z^(-2/3))
    return( Defaults.convertUnits("energy: from eV to atomic", wa) )
end


"""
`Empirical.stoppingProjectileData(projectile::Empirical.AbstractStoppingProjectile)`
    ... to provide the parameters of a projectile as needed in the stopping-power formulas: its squared charge z^2,
        its mass M [m_e] and whether it is a lepton (electron/positron), for which the maximum energy transfer to a
        target electron is limited and the relativistic (Axelrod-type) loss applies.
        A tuple (z2::Float64, M::Float64, lepton::Bool) is returned.
"""
function stoppingProjectileData(projectile::Empirical.AbstractStoppingProjectile)
    if      typeof(projectile) == Empirical.ElectronProjectile   return( (1.0, 1.0, true) )
    elseif  typeof(projectile) == Empirical.PositronProjectile   return( (1.0, 1.0, true) )
    elseif  typeof(projectile) == Empirical.IonProjectile        return( (projectile.z^2, projectile.M, false) )
    else    error("Unsupported projectile $projectile for a stopping-power estimate.")
    end
end


"""
`Empirical.stoppingTargets(material::Empirical.AbstractStoppingMaterial)`
    ... to decompose a material into its electron populations as needed in the stopping-power formulas; each
        population is characterized by its electron density and one characteristic energy: the plasma energy
        hbar omega_p = sqrt(4 pi n_e) [a.u.] for free electrons, and the mean excitation energy Ibar for bound
        electrons. A targets::Array{Tuple{Float64,Float64,String},1} of (density, energy, label) is returned.
"""
function stoppingTargets(material::Empirical.AbstractStoppingMaterial)
    unEnergy = Defaults.getDefaults("unit: energy")
    if      typeof(material) == Empirical.FreeElectronGas
        omegaP = sqrt(4pi * material.ne)
        oPx    = Defaults.convertUnits("energy: from atomic to " * unEnergy, omegaP)
        return( [ (material.ne, omegaP, "free electrons (n_e = $(material.ne) a.u., hbar omega_p = $oPx $unEnergy)") ] )
    elseif  typeof(material) == Empirical.NeutralAtomGas
        nB     = material.Z * material.natom
        iBar   = Empirical.meanExcitationEnergy(material.Z)
        iBx    = Defaults.convertUnits("energy: from atomic to " * unEnergy, iBar)
        return( [ (nB, iBar, "bound electrons (Z n_atom = $nB a.u., mean excitation energy Ibar = $iBx $unEnergy)") ] )
    elseif  typeof(material) == Empirical.PartiallyIonizedGas
        if  material.chie < 0.  ||  material.chie > material.Z
            error("The ionization fraction chie = $(material.chie) must lie between 0 and Z = $(material.Z).")   end
        nB     = (material.Z - material.chie) * material.natom
        nF     = material.chie * material.natom
        iBar   = Empirical.meanExcitationEnergy(material.Z)
        omegaP = sqrt(4pi * nF)
        iBx    = Defaults.convertUnits("energy: from atomic to " * unEnergy, iBar)
        oPx    = Defaults.convertUnits("energy: from atomic to " * unEnergy, omegaP)
        return( [ (nB, iBar,   "bound electrons ((Z - chi_e) n_atom = $nB a.u., Ibar = $iBx $unEnergy)"),
                  (nF, omegaP, "free electrons (chi_e n_atom = $nF a.u., hbar omega_p = $oPx $unEnergy)") ] )
    else    error("Unsupported material $material for a stopping-power estimate.")
    end
end


"""
`Empirical.stoppingPower(energies::Array{Float64,1}, projectile::Empirical.AbstractStoppingProjectile,
                         material::Empirical.AbstractStoppingMaterial, approx::Empirical.Bohr1913;
                         printout::Bool=false)`
    ... to estimate the stopping power - dE/dx of the given projectile at the given kinetic energies [a.u.] in the
        given material, by using Bohr's (1913) classical Coulomb logarithm for each electron population t,
            - dE/dx = sum_t n_t (4 pi z^2 e^4/(m_e v^2)) ln( 1.123 m_e v^3 / (|z| e^2 omega_t) ) ,
        where v = sqrt(2E/M) is the projectile velocity and hbar omega_t the characteristic energy of the population
        (hbar omega_p for free, Ibar for bound electrons). The constant 1.123 = 2 exp(-gamma_E) carries Euler's
        constant gamma_E = 0.5772. Each logarithm is clamped at zero. A sps::Array{Float64,1} [a.u.] is returned.
        Quantity: a stopping power [Hartree/a_o] -- the energy loss per unit path length; it already contains the
            electron densities of the material.

        Note: Bohr's logarithm is the *classical* limit, ln(b_max/b_min) with the adiabatic distance b_max = v/omega_t
              and the classical collision diameter b_min = |z| e^2/(m_e v^2); it applies when the latter exceeds the
              de Broglie wavelength, i.e. for projectile velocities v < |z| e^2/hbar (in a.u.: v < |z|), but still well
              above the thermal/orbital velocities of the target electrons. It refers to a heavy (distinguishable)
              projectile: for electrons the indistinguishability reduces the logarithm by ln 2, cf. KozmaFranson1992.
              Accuracy: a factor ~2 at best; no shell, Barkas (z^3), Bloch or density-effect corrections.
"""
function stoppingPower(energies::Array{Float64,1}, projectile::Empirical.AbstractStoppingProjectile,
                       material::Empirical.AbstractStoppingMaterial, approx::Empirical.Bohr1913;
                       printout::Bool=false)
    (z2, M, lepton) = Empirical.stoppingProjectileData(projectile)
    targets = Empirical.stoppingTargets(material)
    sps     = Float64[]
    for  eps in energies
        v = sqrt(2*eps / M);    sp = 0.
        for  (nt, et, label) in targets
            sp = sp + nt * 4pi * z2 / v^2 * max(0., log(1.123 * v^3 / (sqrt(z2) * et)) )
        end
        push!(sps, sp)
    end

    if  printout
        sa = "Bohr's (1913) classical Coulomb logarithm ln(1.123 m_e v^3/(|z| e^2 omega_t)): the ratio of the " *
             "adiabatic distance b_max = v/omega_t to the classical collision diameter b_min = |z| e^2/(m_e v^2). " *
             "\n    + Regime: classical, i.e. projectile velocities v < |z| a.u. but well above the thermal/orbital " *
             "velocities of the target electrons; heavy (distinguishable) projectile; accuracy a factor ~2"
        Empirical.stoppingPowerPrintout(projectile, material, sa, energies, targets, sps)
    end
    return( sps )
end


"""
`Empirical.stoppingPower(energies::Array{Float64,1}, projectile::Empirical.AbstractStoppingProjectile,
                         material::Empirical.AbstractStoppingMaterial, approx::Empirical.Bethe1931;
                         printout::Bool=false)`
    ... to estimate the stopping power - dE/dx of the given projectile at the given kinetic energies [a.u.] in the
        given material, by using the quantal (Bethe-type) Coulomb logarithm for each electron population t,
            - dE/dx = sum_t n_t (4 pi z^2 e^4/(m_e v^2)) ln( kappa m_e v^2 / (hbar omega_t) ) ,
        where v = sqrt(2E/M) is the projectile velocity and hbar omega_t the characteristic energy of the population
        (hbar omega_p for free, Ibar for bound electrons). For a heavy projectile kappa = 2 (Bethe's original
        logarithm); for an electron or positron projectile kappa = 1, since the maximum energy transfer between
        identical particles is E/2, which reduces the logarithm by ln 2 [ln(m_e v^2/(hbar omega_t)) =
        ln(2 E/(hbar omega_t)), cf. Kozma & Fransson 1992, Eq. 1]. Each logarithm is clamped at zero.
        A sps::Array{Float64,1} [a.u.] is returned.
        Quantity: a stopping power [Hartree/a_o] -- the energy loss per unit path length; it already contains the
            electron densities of the material.

        Note: The quantal logarithm takes b_min from the de Broglie wavelength and applies for projectile velocities
              v > |z| e^2/hbar (in a.u.: v > |z|); for bound electrons additionally v must well exceed the orbital
              velocities, i.e. E >> Ibar. These are bare Born values: no shell, Barkas (z^3), Bloch or density-effect
              corrections and no relativistic terms (use Axelrod1980 for electrons/positrons above ~50 keV).
              Typical accuracy 10-30% well above the characteristic energy.
"""
function stoppingPower(energies::Array{Float64,1}, projectile::Empirical.AbstractStoppingProjectile,
                       material::Empirical.AbstractStoppingMaterial, approx::Empirical.Bethe1931;
                       printout::Bool=false)
    (z2, M, lepton) = Empirical.stoppingProjectileData(projectile)
    kappa   = lepton ? 1.0 : 2.0
    targets = Empirical.stoppingTargets(material)
    sps     = Float64[]
    for  eps in energies
        v = sqrt(2*eps / M);    sp = 0.
        for  (nt, et, label) in targets
            sp = sp + nt * 4pi * z2 / v^2 * max(0., log(kappa * v^2 / et) )
        end
        push!(sps, sp)
    end

    if  printout
        sa = "the quantal (Bethe-type) Coulomb logarithm ln(kappa m_e v^2/(hbar omega_t)) with kappa = $kappa: " *
             "b_min is the de Broglie wavelength" *
             (lepton  ?  "; kappa = 1 since the maximum energy transfer between identical particles is E/2 "  :
                         "; kappa = 2 for a heavy (distinguishable) projectile ") *
             "\n    + Regime: quantal, i.e. projectile velocities v > |z| a.u. and E well above the characteristic " *
             "energies; bare Born values (no shell, Barkas, Bloch, density-effect or relativistic corrections); " *
             "typical accuracy 10-30%"
        Empirical.stoppingPowerPrintout(projectile, material, sa, energies, targets, sps)
    end
    return( sps )
end


"""
`Empirical.stoppingPower(energies::Array{Float64,1}, projectile::Empirical.AbstractStoppingProjectile,
                         material::Empirical.AbstractStoppingMaterial, approx::Empirical.KozmaFranson1992;
                         printout::Bool=false)`
    ... to estimate the stopping power - dE/dx of an *electron* at the given kinetic energies [a.u.] in a
        *free-electron gas*, by using the piecewise electron loss function of Kozma & Fransson (1992, Eqs. 1-3),
            - dE/dx = n_e (2 pi e^4/E) ln( 2 E / (hbar omega_p) )                for  E > 14 eV,
            - dE/dx = n_e (2 pi e^4/E) ln( m_e v^3 / (gamma e^2 omega_p) )       for  kT << E < 14 eV,
        with gamma = exp(gamma_E) = 1.7811 (gamma_E = 0.5772, Euler's constant; Schunk & Hays 1971). The boundary
        at 14 eV is physical: it marks E = m e^4/(2 hbar^2) = 1 Ry, where the de Broglie wavelength overtakes the
        classical distance of closest approach. Both branches carry the same -ln 2 electron-projectile correction
        relative to their heavy-projectile (Bohr/Bethe) counterparts, since 1.123 x 1.7811 = 2. The logarithm is
        clamped at zero. A sps::Array{Float64,1} [a.u.] is returned; any other projectile or material raises an
        informative error, since this loss function is defined just for this one case.
        Quantity: a stopping power [Hartree/a_o] -- the energy loss per unit path length; it already contains the
            electron density n_e.
"""
function stoppingPower(energies::Array{Float64,1}, projectile::Empirical.AbstractStoppingProjectile,
                       material::Empirical.AbstractStoppingMaterial, approx::Empirical.KozmaFranson1992;
                       printout::Bool=false)
    if  typeof(projectile) != Empirical.ElectronProjectile  ||  typeof(material) != Empirical.FreeElectronGas
        error("KozmaFranson1992 is the electron loss function of ApJ 390, 602 and is defined only for an " *
              "ElectronProjectile in a FreeElectronGas; for other projectiles or materials use Bohr1913, " *
              "Bethe1931 or (for relativistic electrons/positrons) Axelrod1980.")
    end
    omegaP  = sqrt(4pi * material.ne)
    e14     = Defaults.convertUnits("energy: from eV to atomic", 14.0)
    targets = Empirical.stoppingTargets(material)
    sps     = Float64[]
    for  eps in energies
        if  eps > e14   push!(sps, material.ne * 2pi / eps * max(0., log(2*eps / omegaP)) )
        else            v = sqrt(2*eps)
                        push!(sps, material.ne * 2pi / eps * max(0., log(v^3 / (1.7811 * omegaP))) )
        end
    end

    if  printout
        sa = "the Kozma & Fransson (1992) electron loss function: quantal logarithm ln(2E/(hbar omega_p)) above " *
             "and classical logarithm ln(m_e v^3/(1.7811 e^2 omega_p)) below E = 14 eV " *
             "\n    + Regime: kT << E; the 14 eV boundary marks 1 Ry, where the de Broglie wavelength overtakes " *
             "the classical closest approach; both branches include the -ln 2 electron-projectile correction"
        Empirical.stoppingPowerPrintout(projectile, material, sa, energies, targets, sps)
    end
    return( sps )
end


"""
`Empirical.stoppingPower(energies::Array{Float64,1}, projectile::Empirical.AbstractStoppingProjectile,
                         material::Empirical.AbstractStoppingMaterial, approx::Empirical.Axelrod1980;
                         printout::Bool=false)`
    ... to estimate the stopping power - dE/dx of an electron or positron at the given kinetic energies [a.u.] in the
        given material, by using Axelrod's (1980) relativistic (Bethe-type) energy loss for each electron
        population t,
            - dE/dx = sum_t n_t (4 pi e^4/(m_e c^2 beta^2)) [ ln( sqrt(gamma-1) gamma beta m_e c^2 / (hbar omega_t) )
                          + 1/2 ln 2 - (beta^2/12) (23/2 + 7/(gamma+1) + 5/(gamma+1)^2 + 2/(gamma+1)^3) ],
        with the Lorentz factor gamma = 1 + E/(m_e c^2) and hbar omega_t the characteristic energy of the population:
        the plasma energy hbar omega_p for free electrons ("plasma losses") and the mean excitation energy Ibar for
        bound electrons ("ionization and excitation losses"); cf. Milne et al. (1999, Eqs. 1-3). The bracket is
        clamped at zero. A sps::Array{Float64,1} [a.u.] is returned; an ion projectile raises an error.
        Quantity: a stopping power [Hartree/a_o] -- the energy loss per unit path length; it already contains the
            electron densities of the material.

        Note: In the nonrelativistic limit the bracket reduces exactly to the quantal electron logarithm
              ln(2 E/(hbar omega_t)) of Bethe1931, since sqrt(gamma-1) gamma beta m_e c^2 -> sqrt(2) E and the
              1/2 ln 2 completes the factor 2; the relativistic terms matter above ~50-100 keV. The (beta^2/12)
              term is the positron-type correction as transcribed by Milne et al., who apply it to electrons and
              positrons alike ("ignoring small differences in the relativistic corrections"); no density-effect
              correction is included, so the values grow logarithmically at highly relativistic energies.
"""
function stoppingPower(energies::Array{Float64,1}, projectile::Empirical.AbstractStoppingProjectile,
                       material::Empirical.AbstractStoppingMaterial, approx::Empirical.Axelrod1980;
                       printout::Bool=false)
    (z2, M, lepton) = Empirical.stoppingProjectileData(projectile)
    if  !lepton
        error("Axelrod1980 provides the relativistic energy loss of electrons and positrons only; no relativistic " *
              "ion stopping is implemented -- for ions with v << c use Bethe1931 (or Bohr1913).")
    end
    c       = Defaults.getDefaults("speed of light: c")
    targets = Empirical.stoppingTargets(material)
    sps     = Float64[]
    for  eps in energies
        gam  = 1.0 + eps / c^2
        bet2 = 1.0 - 1.0 / gam^2
        sp   = 0.
        for  (nt, et, label) in targets
            wa = log( sqrt(gam - 1.0) * gam * sqrt(bet2) * c^2 / et ) + log(2.0)/2   -
                 bet2/12 * ( 23/2 + 7/(gam + 1.0) + 5/(gam + 1.0)^2 + 2/(gam + 1.0)^3 )
            sp = sp + nt * 4pi / (c^2 * bet2) * max(0., wa)
        end
        push!(sps, sp)
    end

    if  printout
        sa = "Axelrod's (1980) relativistic energy loss: the Bethe-type bracket with hbar omega_p (free) or " *
             "Ibar (bound electrons) as effective ionization potential " *
             "\n    + Regime: electrons/positrons up to the MeV range; reduces exactly to Bethe1931 for E << m_e c^2; " *
             "no density-effect correction; positron-type relativistic term applied to e- and e+ alike"
        Empirical.stoppingPowerPrintout(projectile, material, sa, energies, targets, sps)
    end
    return( sps )
end


"""
`Empirical.stoppingPowerPrintout(projectile::Empirical.AbstractStoppingProjectile,
                                 material::Empirical.AbstractStoppingMaterial, sa::String,
                                 energies::Array{Float64,1}, targets::Array{Tuple{Float64,Float64,String},1},
                                 sps::Array{Float64,1})`
    ... to print a common report about the given stopping powers: the projectile and material in words, the
        approximation and its regime of validity (sa), the electron populations that take up the energy, the
        neglected mechanisms, the energies in user-defined units, the stopping powers in eV/cm as well as -- for
        atom gases with known mass number -- the mass stopping power in MeV cm^2/g. Nothing is returned.
"""
function stoppingPowerPrintout(projectile::Empirical.AbstractStoppingProjectile,
                               material::Empirical.AbstractStoppingMaterial, sa::String,
                               energies::Array{Float64,1}, targets::Array{Tuple{Float64,Float64,String},1},
                               sps::Array{Float64,1})
    unEnergy = Defaults.getDefaults("unit: energy")
    aoCm     = Defaults.convertUnits("length: from atomic to cm", 1.0)
    energiex = [Defaults.convertUnits("energy: from atomic to " * unEnergy, eps)   for eps in energies]
    spsx     = [Defaults.convertUnits("energy: from atomic to eV", sp) / aoCm      for sp  in sps]
    sb = "\n* Estimate empirically the stopping power - dE/dx of $projectile in $material " *
         "with the following assumptions/simplifications: " *
         "\n    + Use " * sa * ". " *
         "\n    + All these estimates share - dE/dx = sum_t n_t (4 pi z^2 e^4)/(m_e v^2) ln Lambda_t: the 1/v^2 " *
         "prefactor reflects the interaction time per target electron, and the Coulomb logarithm " *
         "ln Lambda_t = ln(b_max/b_min) counts the effective impact parameters. "
    for  (nt, et, label) in targets
        sb = sb * "\n    + Energy is transferred to " * label * ". "
    end
    sb = sb *
         "\n    + Neglected: nuclear stopping, radiative (bremsstrahlung) losses, electron capture and loss of ions, " *
         "and shell, Barkas (z^3), Bloch and density-effect corrections. " *
         "\n    + Projectile energies [$unEnergy]   = $energiex " *
         "\n    + Stopping powers - dE/dx [eV/cm] = $spsx "
    #  For atom gases with known mass number, also express the mass stopping power - (dE/dx)/rho in MeV cm^2/g.
    if      typeof(material) in [Empirical.NeutralAtomGas, Empirical.PartiallyIonizedGas]
        gPerU   = Defaults.ELECTRON_MASS_IN_G / Defaults.ELECTRON_MASS_U
        rho     = material.natom / aoCm^3 * material.A * gPerU
        spsRho  = [1.0e-6 * sp / rho   for sp in spsx]
        sb = sb *
         "\n    + Mass density rho [g/cm^3] = $rho;  mass stopping power - (dE/dx)/rho [MeV cm^2/g] = $spsRho "
    end
    sb = sb *
         "\n    + Quantity: a stopping power [eV/cm] -- the energy loss per unit path length; it already contains the electron densities of the material. Divide by the electron density for the stopping cross section [eV cm^2], or by rho for the mass stopping power. " * "\n"
    println(sb)

    return( nothing )
end
