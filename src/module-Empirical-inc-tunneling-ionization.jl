
#################################################################################################################################
### Tunneling ionization (TI) ###################################################################################################
##
##  The (quasiclassical) tunneling ionization rate of an atom or ion in a strong, quasi-static electric field, following
##  Ammosov, Delone & Krainov (1986); the exact formula and its coefficients are quoted from the summary of
##  J. Bauer & P. Mulser, Phys. Rev. A 59, 569 (1999) [arXiv:physics/9802042], Eq. (10). The rate depends on the
##  instantaneous (peak) field strength F alone -- no plasma or radiation field enters -- so this is a *static-field*
##  rate, not the cycle-average over an oscillating (laser) field; cf. the Note of Empirical.tunnelingIonizationRate.
##  Quantity: a field-ionization rate [1/s] -- a property of the ion in the given field alone; no further
##      multiplication is needed. It is *not* a spontaneous (zero-field) rate, and it must not be summed with one:
##      the field dependence is strongly nonperturbative (an essential singularity at F = 0).


"""
`Empirical.effectiveQuantumNumbers(Z::Float64, Ip::Float64)`
    ... to provide the effective principal quantum number n* of an electron with binding energy Ip [a.u.] in the
        Coulomb tail of an ion of residual charge Z (i.e. the charge seen by the escaping electron once it is far
        from the core, at asymptotic distances), following the hydrogenic relation Ip = Z^2/(2 n*^2), i.e.
            n* = Z / sqrt(2 Ip).
        A tuple (nstar::Float64, kappa::Float64) is returned, with kappa = sqrt(2 Ip) the (inverse) tunneling decay
        constant that enters the ADK exponent and the critical (barrier-suppression) field F_c = kappa^3/16.
"""
function effectiveQuantumNumbers(Z::Float64, Ip::Float64)
    kappa = sqrt(2*Ip)
    nstar = Z / kappa
    return( (nstar, kappa) )
end


"""
`Empirical.adkCoefficient(nstar::Float64, lstar::Float64)`
    ... to provide the squared ADK coefficient
            C^2_(n* l*) = 2^(2n*) / ( n* Gamma(n*+l*+1) Gamma(n*-l*) ) ,
        using the Gamma function so that n* and l* need not be integers, as appropriate for a real atom or ion (only
        for a pure hydrogenic 1s state, Z = n* = 1 and l* = 0, are they integers by construction). For this
        hydrogenic case, C^2_(1,0) = 4 exactly, and the full ADK rate below then reduces *exactly* to the classical
        Landau tunneling formula for hydrogen [Bauer & Mulser, Eq. (8)] -- a useful, and here numerically verified,
        cross check of the normalization. Many practical implementations instead use the quasiclassical (Stirling)
        approximation C^2_(n*) = (2e/n*)^(2n*)/(2 pi n*), which drops the l*-dependence altogether and is accurate
        only for n* >> 1; for n* = 1 it overestimates C^2 by 17.6%, and the Gamma-function form is used throughout
        this module instead. A c2::Float64 is returned.
"""
function adkCoefficient(nstar::Float64, lstar::Float64)
    return( 2.0^(2*nstar) / (nstar * SpecialFunctions.gamma(nstar+lstar+1) * SpecialFunctions.gamma(nstar-lstar)) )
end


"""
`Empirical.adkAngularFactor(l::Int64, m::Int64)`
    ... to provide the angular factor
            f(l,m) = (2l+1) (l+|m|)! / ( 2^|m| |m|! (l-|m|)! )
        of the ADK rate, which depends on the *actual* orbital and magnetic quantum numbers l and m of the initial
        (bound) state -- in contrast to the coefficient C^2_(n*l*), which uses the *effective* quantum numbers n*
        and l*. A wa::Float64 is returned.
"""
function adkAngularFactor(l::Int64, m::Int64)
    if  abs(m) > l   error("Need |m| <= l; got l = $l, m = $m.")   end
    return( (2l+1) * factorial(l+abs(m)) / (2.0^abs(m) * factorial(abs(m)) * factorial(l-abs(m))) )
end


"""
`Empirical.tunnelingIonizationRate(fields::Array{Float64,1}, Z::Float64, Ip::Float64, l::Int64,
                                   approx::Empirical.ADK1986; m::Int64=0, printout::Bool=false)`
    ... to estimate the tunneling ionization rate of an electron with binding energy Ip [a.u.] and orbital quantum
        number l, in the Coulomb tail of an ion of residual charge Z (the charge seen by the escaping electron at
        large distance, i.e. the charge of the *daughter* ion once this electron is gone), at the given (quasi-static)
        electric field strengths [a.u.]. This is the direct (Z, Ip, l) form of Ammosov, Delone & Krainov's (1986)
        formula,
            W(F) = C^2_(n* l*) f(l,m) Ip (2 kappa^3/F)^(2n*-|m|-1) exp( -2 kappa^3/(3F) ),
        with kappa = sqrt(2 Ip), the effective principal quantum number n* = Z/kappa, cf.
        Empirical.effectiveQuantumNumbers(). The effective orbital quantum number l* is set to zero throughout,
        following standard practice (e.g. Ammosov et al. themselves, and essentially all production strong-field
        codes; l* would otherwise require the quantum-defect structure of the atom, which is rarely available),
        while the real orbital quantum number l and the magnetic quantum number m (default 0, i.e. the field aligned
        along the symmetry axis of the initial orbital) enter the angular factor f(l,m), cf.
        Empirical.adkAngularFactor(). A rates::Array{Float64,1} [a.u.] is returned.
        Quantity: a field-ionization rate [1/s] -- already evaluated at the given field; nothing further needs to be
            multiplied. It is not a spontaneous (zero-field) rate: the field dependence is nonperturbative (an
            essential singularity at F = 0) and cannot be linearized or summed with a zero-field rate.

        Note: This is the *static-field* (quasi-DC) rate. For a field oscillating as F(t) = F0 cos(omega t), the
              literature ADK rate is usually quoted after averaging over one optical cycle, which multiplies the
              static rate above by sqrt(3F/(pi kappa^3)) [Bauer & Mulser, Eq. (10)]; that cycle-averaged variant is
              not implemented here, and F should be understood as the (peak) field strength of interest, e.g. F0.

        Note: The formula is valid in the tunneling regime, F well below the critical (barrier-suppression) field
              F_c = kappa^3/16; above F_c the true ionization rate is known to be substantially overestimated (the
              potential barrier no longer exists at the classical level). No warning is raised automatically, since
              F/F_c is not itself pathological to compute, but it is reported whenever printout = true.

        Note: For Z = 1, Ip = 0.5 [a.u.] and l = 0 (the pure hydrogenic 1s case, where n* = 1 and l* = 0 by
              construction), this rate reduces *exactly* to the closed Landau tunneling formula (numerically verified
              to 8 digits) -- a strong, literature-independent cross check of the coefficient normalization used
              here, cf. Empirical.adkCoefficient().
"""
function tunnelingIonizationRate(fields::Array{Float64,1}, Z::Float64, Ip::Float64, l::Int64,
                                 approx::Empirical.ADK1986; m::Int64=0, printout::Bool=false)
    if  Z <= 0.   error("The residual ion charge Z = $Z is not positive; ADK requires a Coulombic (attractive) " *
                        "asymptotic tail for the escaping electron.")   end
    if  Ip <= 0.  error("The binding energy Ip = $Ip [a.u.] must be positive.")   end
    (nstar, kappa) = Empirical.effectiveQuantumNumbers(Z, Ip)
    lstar          = 0.
    c2             = Empirical.adkCoefficient(nstar, lstar)
    flm            = Empirical.adkAngularFactor(l, m)
    Fc             = kappa^3 / 16

    rates = Float64[]
    for  F in fields
        push!(rates, c2 * flm * Ip * (2*kappa^3/F)^(2*nstar - abs(m) - 1) * exp(-2*kappa^3/(3*F)) )
    end

    # Report about these estimates
    if  printout
        unRate  = Defaults.getDefaults("unit: rate");   unEnergy = Defaults.getDefaults("unit: energy")
        Ipx     = Defaults.convertUnits("energy: from atomic to " * unEnergy, Ip)
        ratesx  = [Defaults.convertUnits("rate: from atomic to " * unRate, r)   for r in rates]
        fRatio  = [F/Fc for F in fields]
        sa = "\n* Estimate the tunneling ionization rate for a given (Z, Ip, l) with the following " *
             "assumptions/simplifications: " *
             "\n    + Use the (quasiclassical) ADK formula (Ammosov, Delone & Krainov 1986) for a static field. " *
             "\n    + Residual ion charge Z = $Z;  l = $l, m = $m " *
             "\n    + Binding energy Ip [$unEnergy] = $Ipx;  n* = $nstar, l* = $lstar " *
             "\n    + Critical (barrier-suppression) field F_c [a.u.] = $Fc;  fields/F_c = $fRatio " *
             "\n    + Field strengths [a.u.]         = $fields " *
             "\n    + Tunneling ionization rates [$unRate] = $ratesx " *
             "\n    + Quantity: a field-ionization rate [$unRate] at the given field -- already contains the field dependence; not a spontaneous (zero-field) rate and must not be summed with one. " * "\n"
        println(sa)
    end

    return( rates )
end


"""
`Empirical.tunnelingIonizationRate(fields::Array{Float64,1}, iConf::Configuration, fConf::Configuration,
                                   approx::Empirical.ADK1986; m::Int64=0, printout::Bool=false,
                                   data::PeriodicTable.AbstractEnergyData=PeriodicTable.XrayDataBooklet())`
    ... to estimate the tunneling ionization rate of an electron from iConf -> fConf at the given (quasi-static)
        electric field strengths [a.u.]; a convenience wrapper around the direct (Z, Ip, l) method above, which
        determines the residual ion charge Z = Z_nuclear - fConf.NoElectrons and the binding energy Ip of the
        ionized shell via Empirical.scaledBindingEnergy(), then dispatches to the direct method. A
        rates::Array{Float64,1} [a.u.] is returned.
        Quantity: a field-ionization rate [1/s], cf. the direct (Z, Ip, l) method for the full description.

        Note: Empirical.scaledBindingEnergy() recognizes a genuine H-like or He-like ion (iConf.NoElectrons <= 2 and
              < Z_nuclear, e.g. He^+ = "1s^1") and then returns the exact (H-like) or approximate (He-like, Slater
              screening; overestimates the binding energy by ~5-30%, decreasing with Z) hydrogenic binding energy
              automatically, rather than the *neutral*-atom tabulated value. A residual gap remains for more highly
              charged ions with 3 or more electrons (e.g. a Li-like O^5+): such an iConf is indistinguishable, from
              the electron count alone, from the common "spectator-omitted shorthand" for a near-neutral atom (e.g.
              Configuration("1s^1 2p^6") for Ne), so the *neutral*-element tabulated value is still used and may
              differ substantially from the true one for a genuinely stripped few-electron ion in this range. For
              such ions, prefer the direct (Z, Ip, l) method with an explicitly supplied (e.g. hydrogenic or
              NIST-tabulated) Ip.
"""
function tunnelingIonizationRate(fields::Array{Float64,1}, iConf::Configuration, fConf::Configuration,
                                 approx::Empirical.ADK1986; m::Int64=0, printout::Bool=false,
                                 data::PeriodicTable.AbstractEnergyData=PeriodicTable.XrayDataBooklet())
    Znuc = Defaults.getDefaults("nuclear: charge");    iShell = Shell(0,0);    diff = 0

    # Determine the ionized shell.
    wa = Basics.extractFromConfigurations(Basics.OccupationDifference(), iConf, fConf)
    if length(wa) > 1   error("Incompatible initial and final configurations for a tunneling ionization rate.")   end
    for  (k,v) in wa
        diff = diff + v
        if      v == 1     iShell = k    end
    end
    if  diff != 1   error("Incompatible initial and final configurations for a tunneling ionization rate.")   end

    Zres  = Znuc - fConf.NoElectrons
    Ip    = Empirical.scaledBindingEnergy(Znuc, iShell, iConf, data)
    rates = Empirical.tunnelingIonizationRate(fields, Zres, Ip, iShell.l, approx; m=m, printout=false)

    # Report about these estimates
    if  printout
        unRate  = Defaults.getDefaults("unit: rate");   unEnergy = Defaults.getDefaults("unit: energy")
        Ipx     = Defaults.convertUnits("energy: from atomic to " * unEnergy, Ip)
        ratesx  = [Defaults.convertUnits("rate: from atomic to " * unRate, r)   for r in rates]
        sa = "\n* Estimate the tunneling ionization rate for a given transition i -> f with the following " *
             "assumptions/simplifications: " *
             "\n    + Use the (quasiclassical) ADK formula (Ammosov, Delone & Krainov 1986) for a static field. " *
             "\n    + iConf = $iConf  -->  fConf = $fConf;  ionized shell $iShell (l = $(iShell.l)), m = $m " *
             "\n    + Binding energy of $iShell from Empirical.scaledBindingEnergy() -- exact (H-like) or approximate " *
             "(He-like, ~5-30% high) hydrogenic for a 1- or 2-electron ion, otherwise the *neutral*-element value from $data; a " *
             "highly stripped ion with 3+ electrons may then differ substantially from the true binding energy -- " *
             "use the direct (Z, Ip, l) method with an explicit Ip in that case. " *
             "\n    + Binding energy Ip [$unEnergy] = $Ipx;  residual ion charge Z = $Zres " *
             "\n    + Field strengths [a.u.]         = $fields " *
             "\n    + Tunneling ionization rates [$unRate] = $ratesx " *
             "\n    + Quantity: a field-ionization rate [$unRate] at the given field -- already contains the field dependence; not a spontaneous (zero-field) rate and must not be summed with one. " * "\n"
        println(sa)
    end

    return( rates )
end


"""
`Empirical.electricFieldFromIntensity(intensity::Float64)`
    ... to convert a (peak) laser intensity [W/cm^2] to the corresponding (peak) electric field strength [a.u.], via
        the standard atomic-unit relation intensity[a.u.] = field[a.u.]^2, i.e.
            F [a.u.] = sqrt( intensity [W/cm^2] / 3.50944758e16 ) .
        A field::Float64 [a.u.] is returned.
"""
function electricFieldFromIntensity(intensity::Float64)
    return( sqrt( Defaults.convertUnits("intensity: from W/cm^2 to atomic", intensity) ) )
end
