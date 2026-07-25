
#################################################################################################################################
### Inelastic ion/atom -- hydrogen collisions (Belyaev-Yakovleva simplified model) ##############################################
##
##  Rate coefficients for neutralization, ion-pair formation, excitation and de-excitation in low-energy collisions
##      A^(Z+1)+ + H^-  <-->  A^Z+(j) + H ,       Z = 0, 1, 2, ... ,
##  important for non-local-thermodynamic-equilibrium (non-LTE) modeling of cool stellar atmospheres, following the
##  simplified model of A. K. Belyaev & S. A. Yakovleva, A&A 606, A147 (2017) [Paper I, Z = 0] and A. K. Belyaev &
##  S. A. Yakovleva, A&A 608, A33 (2017) [the Z = 1, 2 case, whose equation numbers are used throughout below].
##
##  Physics in one paragraph: at low collision energy the reaction proceeds via long-range nonadiabatic transitions
##  where the ionic (Coulomb) A^(Z+1)+ + H^- molecular-state potential crosses a covalent A^Z+(j) + H potential
##  (both taken as simple diabatic curves, Eqs. 6-7); the Landau-Zener model gives the single-pass transition
##  probability at each such crossing from the ionic-covalent coupling matrix element of Olson, Smith & Bauer
##  [Appl. Opt. 10, 1848 (1971)], Eq. (12). This needs only the electronic bound (ionization) energy of the atomic
##  state involved -- no wavefunctions, no CI, no transition matrix elements -- which is the entire point of the
##  model: it "practically does not require any computational effort" (Belyaev & Yakovleva 2017, Sec. 2.4) once the
##  bound energies and statistical weights of the states of interest are known (e.g. from NIST, or from a JAC
##  structure computation where NIST is incomplete).
##
##  Sign/unit convention (matches the paper): all bound energies (Ej, Ei, Ef, EH-) are NEGATIVE, measured from the
##  ionization limit of the A^Z+ ion, in atomic units unless stated otherwise. Z is the (dimensionless) charge of
##  the covalent species A^Z+ (Z=1 for a singly-ionized atom colliding with H, etc.), so the ionic species carries
##  charge Z+1. mu is the ion-H reduced mass in units of m_e (e.g. Defaults.PROTON_MASS_U/Defaults.ELECTRON_MASS_U
##  for atomic hydrogen as the light collision partner).
##
##  Limitations (see the paper for the full discussion): the model targets only rate coefficients with high and
##  moderate values -- those from the *dominant* long-range ionic-covalent crossings -- and deliberately ignores
##  short-range/multichannel effects that matter for small rates. Selecting *which* atomic states A^Z+(j) actually
##  correlate (via a one-electron transition) to the ionic ground molecular-state symmetry, and their statistical
##  probabilities p_stat_j within that symmetry, is a separate, molecular-term-symmetry bookkeeping step (Steps 1-2
##  of the paper) that this module does *not* automate -- the channel list and its p_stat values must be supplied
##  by the caller (as in the paper's own worked example, cf. Empirical.InelasticHChannel below).


"""
`Empirical.hydrogenAnionEnergy()`
    ... returns the electronic bound energy of the hydrogen anion H^-, E_H- = -0.754 eV (Belyaev & Yakovleva 2017,
        Sec. 2.2), converted to atomic units. A value::Float64 is returned.
"""
function hydrogenAnionEnergy()
    return( Defaults.convertUnits("energy: from eV to atomic", -0.754) )
end


"""
`Empirical.nonadiabaticRadius(Ej::Float64, EHminus::Float64, Z::Float64)`
    ... to compute the internuclear distance R_j at which the ionic A^(Z+1)+ + H^- (Coulomb) diabatic potential
        crosses the covalent A^Z+(j) + H (flat) diabatic potential,
            R_j = (Z+1) / (E_H- - E_j) ,
        Eq. (8) of Belyaev & Yakovleva, A&A 608, A33 (2017). Ej and EHminus are electronic bound energies [a.u.,
        negative]. A value::Float64 [a.u.] is returned.
"""
function nonadiabaticRadius(Ej::Float64, EHminus::Float64, Z::Float64)
    if  EHminus <= Ej   error("EHminus = $EHminus must be deeper bound (more negative) than Ej = $Ej.")   end
    return( (Z + 1.0) / (EHminus - Ej) )
end


"""
`Empirical.ionicCovalentCoupling(Ej::Float64, EHminus::Float64, Rj::Float64)`
    ... to compute the (diabatic) ionic-covalent coupling matrix element H_ionic,j(R_j) by the semiempirical
        one-electron-transition formula of Olson, Smith & Bauer [Appl. Opt. 10, 1848 (1971)],
            H_ionic,j = [sqrt(-Ej) + sqrt(-E_H-)] / sqrt(2) * sqrt(Ej E_H-) * R_j
                        * exp( -0.86 R_j [sqrt(-Ej) + sqrt(-E_H-)] / sqrt(2) ) ,
        Eq. (12) of Belyaev & Yakovleva, A&A 608, A33 (2017). A value::Float64 [a.u.] is returned.
"""
function ionicCovalentCoupling(Ej::Float64, EHminus::Float64, Rj::Float64)
    s = sqrt(-Ej) + sqrt(-EHminus)
    return( s / sqrt(2.0) * sqrt(Ej * EHminus) * Rj * exp(-0.86 * Rj * s / sqrt(2.0)) )
end


"""
`Empirical.landauZenerProbability(Ej::Float64, EHminus::Float64, Z::Float64, v::Float64)`
    ... to compute the single-pass Landau-Zener nonadiabatic transition probability at the ionic-covalent crossing
        R_j (Empirical.nonadiabaticRadius) for a radial velocity v [a.u.] there,
            p_j = exp( -2 pi H_ionic,j^2 (Z+1) / [(E_H- - Ej)^2 v] ) ,
        Eq. (11) of Belyaev & Yakovleva, A&A 608, A33 (2017) (the Eq. (10) definition with the diabatic potentials
        of Eqs. (6)-(7) already substituted). A value::Float64 (0 <= p_j <= 1) is returned.
"""
function landauZenerProbability(Ej::Float64, EHminus::Float64, Z::Float64, v::Float64)
    Rj  = Empirical.nonadiabaticRadius(Ej, EHminus, Z)
    Hij = Empirical.ionicCovalentCoupling(Ej, EHminus, Rj)
    return( exp( -2*pi * Hij^2 * (Z + 1.0) / ( (EHminus - Ej)^2 * v ) ) )
end


"""
`Empirical.maxOrbitalAngularMomentum(availableEnergy::Float64, R::Float64, mu::Float64)`
    ... to compute the maximum total angular momentum quantum number J for which the radial velocity at the
        crossing radius R remains real for the given availableEnergy [a.u.] (the collision energy plus/minus any
        energy defect relevant for the crossing under consideration, cf. Eqs. (13), (18), (19)), i.e. the largest
        (non-negative) integer J with J(J+1) <= 2 mu R^2 availableEnergy. Returns -1 if even J=0 is classically
        closed (availableEnergy <= 0), signalling an empty J range. An n::Int64 is returned.
"""
function maxOrbitalAngularMomentum(availableEnergy::Float64, R::Float64, mu::Float64)
    if  availableEnergy <= 0.   return( -1 )   end
    Jmax = floor(Int64, ( -1.0 + sqrt(1.0 + 8.0*mu*R^2*availableEnergy) ) / 2.0)
    return( max(Jmax, 0) )
end


"""
`Empirical.neutralizationCrossSection(E::Float64, Ef::Float64, EHminus::Float64, Z::Float64, mu::Float64;
                                      printout::Bool=false)`
    ... to compute the (two-channel-approximation) cross section for the neutralization process
        A^(Z+1)+ + H^- -> A^Z+(f) + H at collision energy E [a.u.], by summing the Landau-Zener transition
        probability P^N_if(J,E) = 2 p_f (1-p_f) (Eq. 9) over total angular momentum J,
            sigma^N_if(E) = pi / (2 mu E) * sum_J (2J+1) P^N_if(J,E) ,
        Eq. (5) of Belyaev & Yakovleva, A&A 608, A33 (2017), with the radial velocity at the crossing from Eq. (13).
        A sigma::Float64 [a.u.] is returned.
        Quantity: a cross section [a.u.] -- fold with a Maxwellian velocity distribution, or use
            Empirical.neutralizationReducedRate, to obtain a (reduced) rate coefficient.
"""
function neutralizationCrossSection(E::Float64, Ef::Float64, EHminus::Float64, Z::Float64, mu::Float64;
                                    printout::Bool=false)
    if  E <= 0.   return( 0.0 )   end
    Rf   = Empirical.nonadiabaticRadius(Ef, EHminus, Z)
    Jmax = Empirical.maxOrbitalAngularMomentum(E + EHminus - Ef, Rf, mu)
    sigma = 0.0
    for  J = 0:Jmax
        v2 = (2.0/mu) * ( E + EHminus - Ef - J*(J+1)/(2.0*mu*Rf^2) )
        if  v2 <= 0.   continue   end
        pf = Empirical.landauZenerProbability(Ef, EHminus, Z, sqrt(v2))
        sigma = sigma + (2J+1) * 2*pf*(1-pf)
    end
    sigma = pi / (2*mu*E) * sigma

    if  printout
        unCs = Defaults.getDefaults("unit: cross section");   unEnergy = Defaults.getDefaults("unit: energy")
        Ex   = Defaults.convertUnits("energy: from atomic to " * unEnergy, E)
        sigx = Defaults.convertUnits("cross section: from atomic to " * unCs, sigma)
        println("\n* Estimate empirically (Belyaev-Yakovleva 2017 simplified model) the neutralization cross section " *
                "for A^($(Z+1))+ + H^- -> A^($Z)+(f) + H: " *
                "\n    + Collision energy [$unEnergy] = $Ex;  final-state bound energy Ef [a.u.] = $Ef;  Z = $Z " *
                "\n    + Cross section [$unCs] = $sigx " *
                "\n    + Quantity: cross section [$unCs] -- fold with a Maxwellian distribution, or use " *
                "Empirical.neutralizationReducedRate, for a (reduced) rate coefficient. \n")
    end

    return( sigma )
end


"""
`Empirical.neutralizationReducedRate(T::Float64, Ef::Float64, EHminus::Float64, Z::Float64, mu::Float64;
                                     printout::Bool=false, zerosGL::Int64=64)`
    ... to compute the reduced rate coefficient N_if(T; Ef) for the neutralization process A^(Z+1)+ + H^- ->
        A^Z+(f) + H at temperature T [a.u.], by folding Empirical.neutralizationCrossSection with the Maxwell
        distribution,
            N_if(T; Ef) = sqrt(8/(pi mu (kT)^3)) int_0^infty sigma^N_if(E) E exp(-E/kT) dE ,
        Eq. (4) of Belyaev & Yakovleva, A&A 608, A33 (2017), evaluated by Gauss-Legendre quadrature over
        [0, 30 kT] (converged to well below 1% for this smooth, exponentially-cut-off integrand). A value::Float64
        [a.u.] is returned.
        Quantity: a (reduced) rate coefficient [cm^3/s] -- multiply by the statistical probability p_stat of the
            initial (ionic) molecular state, Eq. (1), to obtain the actual rate coefficient K_if(T); cf.
            Empirical.detailedBalanceRate for the inverse (ion-pair formation) process.

        Validated against Belyaev & Yakovleva (2017): N_if(T=6000 K, Ef=-4.29573 eV) for Ba2+ + H- -> Ba+(6d) + H
        gives 7.588e-8 cm^3/s here vs. the paper's quoted 7.59e-8 cm^3/s -- agreement to 3 significant figures.
"""
function neutralizationReducedRate(T::Float64, Ef::Float64, EHminus::Float64, Z::Float64, mu::Float64;
                                   printout::Bool=false, zerosGL::Int64=64)
    gridGL = Radial.GridGL(Radial.GridGaussLegendreFinite(), 0., 30*T, zerosGL; printout=false)
    s = 0.0
    for  n = 1:gridGL.nt
        E = gridGL.t[n]
        s = s + Empirical.neutralizationCrossSection(E, Ef, EHminus, Z, mu) * E * exp(-E/T) * gridGL.wt[n]
    end
    N = sqrt( 8.0 / (pi * mu * T^3) ) * s

    if  printout
        unEnergy = Defaults.getDefaults("unit: energy")
        Tx    = Defaults.convertUnits("temperature: from atomic to Kelvin", T)
        Efx   = Defaults.convertUnits("energy: from atomic to " * unEnergy, Ef)
        factor = Defaults.convertUnits("length: from atomic to cm", 1.0)^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)
        Nx = factor * N
        println("\n* Estimate empirically (Belyaev-Yakovleva 2017 simplified model) the neutralization reduced rate " *
                "coefficient for A^($(Z+1))+ + H^- -> A^($Z)+(f) + H: " *
                "\n    + Temperature T [K] = $Tx;  final-state bound energy Ef [$unEnergy] = $Efx;  Z = $Z " *
                "\n    + Reduced rate coefficient N_if [cm^3/s] = $Nx " *
                "\n    + Quantity: multiply by the statistical probability of the ionic state (Eq. 1) for the actual " *
                "rate coefficient K_if(T); use Empirical.detailedBalanceRate for the inverse (ion-pair formation) rate. \n")
    end

    return( N )
end


"""
`Empirical.deExcitationCrossSection(E::Float64, Ei::Float64, Ef::Float64, EHminus::Float64, Z::Float64, mu::Float64;
                                    printout::Bool=false)`
    ... to compute the (three-channel-approximation) cross section for the de-excitation process
        A^Z+(i) + H -> A^Z+(f) + H at collision energy E [a.u.] (Ei more deeply bound than Ef, i.e. Ei < Ef is NOT
        required by this function -- the caller must ensure Ei corresponds to the higher-lying, i.e. less deeply
        bound, initial state for a genuine de-excitation), by summing the two-crossing Landau-Zener transition
        probability P^D_if(J,E) = 2 p_f (1-p_f) (1-p_i) (Eq. 17) over total angular momentum J,
            sigma^D_if(E) = pi / (2 mu E) * sum_J (2J+1) P^D_if(J,E) ,
        Eq. (16) of Belyaev & Yakovleva, A&A 608, A33 (2017), with the radial velocities at the two crossings from
        Eqs. (18)-(19). A sigma::Float64 [a.u.] is returned.
        Quantity: a cross section [a.u.] -- fold with a Maxwellian velocity distribution, or use
            Empirical.deExcitationReducedRate, to obtain a (reduced) rate coefficient.

        OPEN ISSUE, PAUSED (not resolved -- deliberately set aside, see project_inelastic_h_collisions.md for the
        full trail): checked against Belyaev & Yakovleva (2017), Table 2, for the Ba+(7p) -> Ba+(6d) de-excitation
        rate coefficient at T=6000 K, this formula (after multiplying the reduced rate by p_stat_i) gives
        ~4.6e-10 cm^3/s vs. the paper's quoted 8.94e-10 cm^3/s -- a factor ~1.9 too small, growing to ~2.6 for
        larger energy defects (checked across all 6 de-excitation pairs from the Ba+(7p) row of Table 2). Neither
        dropping the (1-p_i) factor entirely nor removing the Jmax_i restriction resolves this (tested explicitly).
        The likely ROOT CAUSE, identified via Belyaev, Phys. Rev. A 48, 4299 (1993)
        [examples/papers/a93.pra-belyaev-charge-exchange.pdf]: this function treats de-excitation as an isolated
        3-state (i, ionic, f) problem, ignoring that the ionic curve also crosses the OTHER covalent channels of
        the same ion along the way -- that paper's Eq. (3.8) (a general multichannel Landau-Zener branching
        formula) and its accompanying discussion explicitly warn that neglecting such channels can inflate the
        result by a factor of 2 (bounded, for channels beyond the target) or more (unbounded, for channels between
        the initial and target states) -- consistent in both magnitude and trend with what is observed here. A
        first attempt to apply Eq. (3.8) directly (treating the ionic curve as entrance, excluding state i, and
        applying the formula to the remaining channels) made the result *worse* (2 orders of magnitude too small),
        so the correct mapping between the two papers' formalisms is not yet understood. By contrast,
        Empirical.neutralizationReducedRate agrees with the paper to 3 significant figures for the same system.
        Treat de-excitation/excitation rates from this function as order-of-magnitude/qualitative-trend estimates
        only; do not rely on them for a quantitative comparison.
"""
function deExcitationCrossSection(E::Float64, Ei::Float64, Ef::Float64, EHminus::Float64, Z::Float64, mu::Float64;
                                  printout::Bool=false)
    if  E <= 0.   return( 0.0 )   end
    Ri = Empirical.nonadiabaticRadius(Ei, EHminus, Z);   Rf = Empirical.nonadiabaticRadius(Ef, EHminus, Z)
    JmaxI = Empirical.maxOrbitalAngularMomentum(E, Ri, mu)
    JmaxF = Empirical.maxOrbitalAngularMomentum(E + Ei - Ef, Rf, mu)
    Jmax  = min(JmaxI, JmaxF)
    sigma = 0.0
    for  J = 0:Jmax
        v2i = (2.0/mu) * ( E - J*(J+1)/(2.0*mu*Ri^2) )
        v2f = (2.0/mu) * ( E + Ei - Ef - J*(J+1)/(2.0*mu*Rf^2) )
        if  v2i <= 0.  ||  v2f <= 0.   continue   end
        pi_ = Empirical.landauZenerProbability(Ei, EHminus, Z, sqrt(v2i))
        pf  = Empirical.landauZenerProbability(Ef, EHminus, Z, sqrt(v2f))
        sigma = sigma + (2J+1) * 2*pf*(1-pf)*(1-pi_)
    end
    sigma = pi / (2*mu*E) * sigma

    if  printout
        unCs = Defaults.getDefaults("unit: cross section");   unEnergy = Defaults.getDefaults("unit: energy")
        Ex   = Defaults.convertUnits("energy: from atomic to " * unEnergy, E)
        sigx = Defaults.convertUnits("cross section: from atomic to " * unCs, sigma)
        println("\n* Estimate empirically (Belyaev-Yakovleva 2017 simplified model) the de-excitation cross section " *
                "for A^($Z)+(i) + H -> A^($Z)+(f) + H: " *
                "\n    + Collision energy [$unEnergy] = $Ex;  Ei [a.u.] = $Ei;  Ef [a.u.] = $Ef;  Z = $Z " *
                "\n    + Cross section [$unCs] = $sigx " *
                "\n    + Quantity: cross section [$unCs] -- fold with a Maxwellian distribution, or use " *
                "Empirical.deExcitationReducedRate, for a (reduced) rate coefficient. \n")
    end

    return( sigma )
end


"""
`Empirical.deExcitationReducedRate(T::Float64, Ei::Float64, Ef::Float64, EHminus::Float64, Z::Float64, mu::Float64;
                                   printout::Bool=false, zerosGL::Int64=64)`
    ... to compute the reduced rate coefficient D_if(T; Ei, Ef) for the de-excitation process A^Z+(i) + H ->
        A^Z+(f) + H at temperature T [a.u.], by folding Empirical.deExcitationCrossSection with the Maxwell
        distribution,
            D_if(T; Ei, Ef) = sqrt(8/(pi mu (kT)^3)) int_0^infty sigma^D_if(E) E exp(-E/kT) dE ,
        Eq. (15) of Belyaev & Yakovleva, A&A 608, A33 (2017), evaluated by Gauss-Legendre quadrature over
        [0, 30 kT]. A value::Float64 [a.u.] is returned.
        Quantity: a (reduced) rate coefficient [cm^3/s] -- multiply by the statistical probability p_stat of the
            initial state i, Eq. (2), to obtain the actual rate coefficient K_if(T); cf.
            Empirical.detailedBalanceRate for the inverse (excitation) process.
"""
function deExcitationReducedRate(T::Float64, Ei::Float64, Ef::Float64, EHminus::Float64, Z::Float64, mu::Float64;
                                 printout::Bool=false, zerosGL::Int64=64)
    gridGL = Radial.GridGL(Radial.GridGaussLegendreFinite(), 0., 30*T, zerosGL; printout=false)
    s = 0.0
    for  n = 1:gridGL.nt
        E = gridGL.t[n]
        s = s + Empirical.deExcitationCrossSection(E, Ei, Ef, EHminus, Z, mu) * E * exp(-E/T) * gridGL.wt[n]
    end
    D = sqrt( 8.0 / (pi * mu * T^3) ) * s

    if  printout
        unEnergy = Defaults.getDefaults("unit: energy")
        Tx    = Defaults.convertUnits("temperature: from atomic to Kelvin", T)
        Eix   = Defaults.convertUnits("energy: from atomic to " * unEnergy, Ei)
        Efx   = Defaults.convertUnits("energy: from atomic to " * unEnergy, Ef)
        factor = Defaults.convertUnits("length: from atomic to cm", 1.0)^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)
        Dx = factor * D
        println("\n* Estimate empirically (Belyaev-Yakovleva 2017 simplified model) the de-excitation reduced rate " *
                "coefficient for A^($Z)+(i) + H -> A^($Z)+(f) + H: " *
                "\n    + Temperature T [K] = $Tx;  Ei [$unEnergy] = $Eix;  Ef [$unEnergy] = $Efx;  Z = $Z " *
                "\n    + Reduced rate coefficient D_if [cm^3/s] = $Dx " *
                "\n    + Quantity: multiply by the statistical probability of state i (Eq. 2) for the actual rate " *
                "coefficient K_if(T); use Empirical.detailedBalanceRate for the inverse (excitation) rate. \n")
    end

    return( D )
end


"""
`Empirical.detailedBalanceRate(Kjk::Float64, pstatJ::Float64, pstatK::Float64, deltaEjk::Float64, T::Float64)`
    ... to compute the rate coefficient K_kj(T) for the endothermic inverse (ion-pair formation, or excitation)
        process k -> j from the rate coefficient Kjk(T) of the exothermic direct (neutralization, or
        de-excitation) process j -> k, by the detailed balance relation
            K_kj(T) = Kjk(T) * (p_stat_k / p_stat_j) * exp( -deltaEjk / kT ) ,
        Eq. (3) of Belyaev & Yakovleva, A&A 608, A33 (2017), where deltaEjk = Ej - Ek > 0 is the (positive) energy
        defect of the exothermic j -> k process. A value::Float64 [same unit as Kjk] is returned.
"""
function detailedBalanceRate(Kjk::Float64, pstatJ::Float64, pstatK::Float64, deltaEjk::Float64, T::Float64)
    if  deltaEjk < 0.   error("deltaEjk = $deltaEjk must be positive (Ej - Ek for the exothermic j -> k process).")   end
    return( Kjk * (pstatK / pstatJ) * exp(-deltaEjk / T) )
end


"""
`struct  Empirical.InelasticHChannel`
    ... a single atomic/ionic scattering channel A^Z+(j) + H that correlates to the ionic ground molecular state
        A^(Z+1)+ + H^- (Belyaev & Yakovleva 2017, Steps 1-2); this module does not derive the channel list or its
        statistical probabilities automatically -- see the module note above.

    + name        ::String    ... a short label for the channel, e.g. "Ba+(6d 2D)", for display purposes only.
    + E           ::Float64   ... electronic bound energy Ej [a.u.], negative, measured from the ionization limit.
    + pstat       ::Float64   ... statistical probability p_stat_j of this channel within the relevant molecular
                                  symmetry (Eqs. 1-2); e.g. taken from Table 1 of Belyaev & Yakovleva (2017) for a
                                  known case, or derived from the relevant molecular-term-symmetry correlation.
"""
struct  InelasticHChannel
    name          ::String
    E             ::Float64
    pstat         ::Float64
end


# `Base.show(io::IO, ch::InelasticHChannel)`  ... provides a String notation for the variable ch::InelasticHChannel.
function Base.show(io::IO, ch::InelasticHChannel)
    print(io, "$(ch.name):  E = $(ch.E) [a.u.],  pstat = $(ch.pstat)")
end


"""
`Empirical.inelasticHCollisionRateMatrix(T::Float64, channels::Array{Empirical.InelasticHChannel,1},
                                         pstatIonic::Float64, Z::Float64, mu::Float64; printout::Bool=false,
                                         zerosGL::Int64=64)`
    ... to compute the full rate-coefficient matrix K_if(T) for all neutralization, ion-pair formation, excitation
        and de-excitation processes among the given list of atomic/ionic channels and the ionic A^(Z+1)+ + H^-
        channel, following Steps 3-4 of Belyaev & Yakovleva, A&A 608, A33 (2017): for each channel j, the
        neutralization rate K(ionic -> j) = pstatIonic * Empirical.neutralizationReducedRate(T, Ej, ...), Eq. (1),
        and the inverse ion-pair-formation rate K(j -> ionic) from Empirical.detailedBalanceRate, Eq. (3); for each
        pair (i,k) of channels with Ei > Ek (i less deeply bound, i.e. the higher-lying state), the de-excitation
        rate K(i -> k) = pstat_i * Empirical.deExcitationReducedRate(T, Ei, Ek, ...), Eq. (2), and the inverse
        excitation rate K(k -> i) from detailed balance. A named tuple
        (channels=channels, ionicToChannel::Array{Float64,1}=, channelToIonic::Array{Float64,1}=,
        channelToChannel::Array{Float64,2}=) is returned, where channelToChannel[i,k] is the rate for
        channels[i] -> channels[k] (the diagonal is left at 0.0).
        Quantity: all rate coefficients [a.u.]; multiply by
        (Defaults.convertUnits("length: from atomic to cm",1.0)^3 / Defaults.convertUnits("time: from atomic to sec",1.0))
        for the conventional [cm^3/s] used throughout the non-LTE literature (as in the printout below).
"""
function inelasticHCollisionRateMatrix(T::Float64, channels::Array{Empirical.InelasticHChannel,1},
                                       pstatIonic::Float64, Z::Float64, mu::Float64; printout::Bool=false,
                                       zerosGL::Int64=64)
    EHminus = Empirical.hydrogenAnionEnergy()
    nch     = length(channels)

    ionicToChannel = zeros(Float64, nch);   channelToIonic = zeros(Float64, nch)
    for  j = 1:nch
        Nj = Empirical.neutralizationReducedRate(T, channels[j].E, EHminus, Z, mu; zerosGL=zerosGL)
        ionicToChannel[j] = pstatIonic * Nj
        deltaE = EHminus - channels[j].E                          # ionic (j="ionic" in Eq. 3) -> channel[j] is exothermic
        channelToIonic[j] = Empirical.detailedBalanceRate(ionicToChannel[j], pstatIonic, channels[j].pstat, deltaE, T)
    end

    channelToChannel = zeros(Float64, nch, nch)
    for  i = 1:nch,  k = 1:nch
        if  i == k   continue   end
        if  channels[i].E > channels[k].E                                # i less deeply bound -> genuine de-excitation i -> k
            Dik = Empirical.deExcitationReducedRate(T, channels[i].E, channels[k].E, EHminus, Z, mu; zerosGL=zerosGL)
            channelToChannel[i,k] = channels[i].pstat * Dik
        end
    end
    for  i = 1:nch,  k = 1:nch
        if  channels[i].E < channels[k].E                                # fill the endothermic (excitation) half by detailed balance
            deltaE = channels[k].E - channels[i].E
            channelToChannel[i,k] = Empirical.detailedBalanceRate(channelToChannel[k,i], channels[k].pstat, channels[i].pstat, deltaE, T)
        end
    end

    if  printout
        factor = Defaults.convertUnits("length: from atomic to cm", 1.0)^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)
        Tx = Defaults.convertUnits("temperature: from atomic to Kelvin", T)
        println("\n* Estimate empirically (Belyaev-Yakovleva 2017 simplified model) the rate-coefficient matrix for " *
                "$nch atomic/ionic channels + the ionic A^($(Z+1))+ + H^- channel at T [K] = $Tx: " *
                "\n    + ionic->j  = the NEUTRALIZATION rate coefficient A^($(Z+1))+ + H^- -> A^($Z)+(j) + H " *
                "(exothermic; well validated against Belyaev & Yakovleva 2017, agrees to 3 significant figures for " *
                "the systems checked so far). " *
                "\n    + j->ionic  = the (much smaller) inverse ION-PAIR-FORMATION rate coefficient A^($Z)+(j) + H " *
                "-> A^($(Z+1))+ + H^- (endothermic; obtained from ionic->j by detailed balance, Empirical.detailedBalanceRate, " *
                "using channel j's own statistical weight p_stat_j and its energy defect from the ionic state). " *
                "\n    + channel[i]->channel[k] = the DE-EXCITATION rate coefficient A^($Z)+(i) + H -> A^($Z)+(k) + H " *
                "between two of the covalent channels themselves (i more weakly bound than k), with the reverse " *
                "(excitation) entry filled in by detailed balance; 0.0 on the diagonal and where i is NOT more " *
                "weakly bound than k (that direction is filled by its own detailed-balance partner instead). " *
                "\n    + KNOWN OPEN ISSUE: unlike ionic->j/j->ionic, the channel[i]->channel[k] (de-)excitation " *
                "block is only order-of-magnitude/qualitative-trend reliable -- checked against Belyaev & Yakovleva " *
                "(2017) Table 2, this module's values are a factor ~1.9-2.6 too small (root cause diagnosed, not yet " *
                "fixed; see Empirical.deExcitationCrossSection's docstring and project_inelastic_h_collisions.md). " *
                "Do not use this block quantitatively. " *
                "\n    + Quantity: every number above is a RATE COEFFICIENT [cm^3/s], not yet a rate. Multiply " *
                "ionic->j by the number density of H^- [cm^-3] to get the actual neutralization rate [1/s] for one " *
                "A^($(Z+1))+ ion; multiply j->ionic and every channel[i]->channel[k] entry by the number density of " *
                "neutral H [cm^-3] instead, since those processes collide with H, not H^-. These per-ion rates are " *
                "exactly the input a non-LTE statistical-equilibrium (population rate-equation) calculation needs, " *
                "which is the whole motivation behind this model (Belyaev & Yakovleva 2017, Sec. 1). ")
        println("    ionic -> channel[j]  and  channel[j] -> ionic  rate coefficients [cm^3/s]:")
        for  j = 1:nch
            println("      $(channels[j].name):   ionic->j = $(factor*ionicToChannel[j])   j->ionic = $(factor*channelToIonic[j])")
        end
        println("    channel[i] -> channel[k] rate coefficients [cm^3/s] (0.0 on the diagonal and where not applicable; " *
                "NOT quantitatively reliable, see above):")
        for  i = 1:nch
            println("      $(channels[i].name):  " * join([ (@sprintf "%.3e" factor*channelToChannel[i,k]) for k = 1:nch ], "  "))
        end
        println()
    end

    return( (channels = channels, ionicToChannel = ionicToChannel, channelToIonic = channelToIonic,
             channelToChannel = channelToChannel) )
end
