
#################################################################################################################################
### Inelastic ion/atom -- hydrogen collisions (Belyaev-Yakovleva simplified model) ##############################################
#
#  Rate coefficients for neutralization, ion-pair formation, excitation and de-excitation in low-energy collisions
#      A^(Z+1)+ + H^-  <-->  A^Z+(j) + H ,       Z = 0, 1, 2, ... ,
#  important for non-local-thermodynamic-equilibrium (non-LTE) modeling of cool stellar atmospheres, following the
#  simplified model of A. K. Belyaev & S. A. Yakovleva, A&A 606, A147 (2017) [Paper I, Z = 0] and A. K. Belyaev &
#  S. A. Yakovleva, A&A 608, A33 (2017) [the Z = 1, 2 case, whose equation numbers are used throughout below].
#
#  Physics in one paragraph: at low collision energy the reaction proceeds via long-range nonadiabatic transitions
#  where the ionic (Coulomb) A^(Z+1)+ + H^- molecular-state potential crosses a covalent A^Z+(j) + H potential
#  (both taken as simple diabatic curves, Eqs. 6-7); the Landau-Zener model gives the single-pass transition
#  probability at each such crossing from the ionic-covalent coupling matrix element of Olson, Smith & Bauer
#  [Appl. Opt. 10, 1848 (1971)], Eq. (12). This needs only the electronic bound (ionization) energy of the atomic
#  state involved -- no wavefunctions, no CI, no transition matrix elements -- which is the entire point of the
#  model: it "practically does not require any computational effort" (Belyaev & Yakovleva 2017, Sec. 2.4) once the
#  bound energies and statistical weights of the states of interest are known (e.g. from NIST, or from a JAC
#  structure computation where NIST is incomplete).
#
#  Sign/unit convention (matches the paper): all bound energies (Ej, Ei, Ef, EH-) are NEGATIVE, measured from the
#  ionization limit of the A^Z+ ion, in atomic units unless stated otherwise. Z is the (dimensionless) charge of
#  the covalent species A^Z+ (Z=1 for a singly-ionized atom colliding with H, etc.), so the ionic species carries
#  charge Z+1. mu is the ion-H reduced mass in units of m_e (e.g. Defaults.PROTON_MASS_U/Defaults.ELECTRON_MASS_U
#  for atomic hydrogen as the light collision partner).
#
#  Limitations (see the paper for the full discussion): the model targets only rate coefficients with high and
#  moderate values -- those from the *dominant* long-range ionic-covalent crossings -- and deliberately ignores
#  short-range/multichannel effects that matter for small rates. De-excitation/excitation rates specifically are
#  known to be unreliable by a factor ~1.9-2.6 (root cause diagnosed -- neglected multichannel branching among the
#  OTHER covalent channels the ionic curve also crosses -- but not fixed; see Empirical.deExcitationCrossSection's
#  own docstring and project_inelastic_h_collisions.md for the full trail, including a documented, unsuccessful
#  attempt to apply Belyaev's (1993) general multichannel formula). Neutralization/ion-pair-formation rates, by
#  contrast, are validated to 3 significant figures against the paper's own tables.
#
#  This one file covers three layers of interface, low-level to high-level:
#    1. Raw physics (this section): operates directly on bound energies (Ej, Ei, Ef), Z (the COVALENT species' own
#       charge, not the nuclear charge), and mu (reduced mass) -- no knowledge of JAC's Configuration/Nuclear.Model
#       entities is needed here.
#    2. Molecular-symmetry channel correlation (below): automates the Wigner-Witmer-style bookkeeping (Steps 1-2 of
#       the 2017 paper) that decides which atomic levels A^Z+(j) even correlate to the ionic entrance term, and
#       with what statistical weight, given just (L,S) of both sides; verified exactly against Table 1 of the paper.
#    3. Empirical.InelasticHReaction (below): the recommended entry point for a NEW reaction -- specifies
#       everything via JAC's own Nuclear.Model and Configuration, with reduced mass, symmetry, and statistical
#       weights all derived internally. Level ENERGIES are the one thing this file deliberately does not compute
#       itself -- always supplied by the caller (from a real SCF run via Atomic.Computation/perform, from
#       literature, or otherwise), keeping this module free of a dependency on the heavier Atomic/SelfConsistent
#       machinery, and cleanly separating "how were these energies obtained" from "how are rates computed from
#       them." Scope is deliberately restricted to a closed-shell entrance ion and a single-active-electron
#       transfer per final configuration; anything else raises an informative error rather than a wrong number.


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
        properly re-derived attempt at the correct multichannel mapping (25-Jul-2026, see the memory file for the
        full trail) still did NOT succeed -- the correction moved the result FURTHER from the literature value, and
        a sanity check of Eq. (3.8) itself revealed an unresolved inconsistency in even the simplest nontrivial
        case. By contrast, Empirical.neutralizationReducedRate agrees with the paper to 3 significant figures for
        the same system. Treat de-excitation/excitation rates from this function as order-of-magnitude/qualitative-
        trend estimates only; do not rely on them for a quantitative comparison.
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


#################################################################################################################################
### Molecular-symmetry channel correlation #######################################################################################
#
#  Automates Steps 1-2 of Belyaev & Yakovleva (2017): given the ionic entrance term (Lion,Sion) and a candidate channel's own
#  term (Lj,Sj), decides whether that channel correlates to the entrance by molecular symmetry and, if so, with what
#  statistical weight -- the p_stat_j field of Empirical.InelasticHChannel below. The paper states only that "the statistical
#  probabilities p_stat_j can be readily calculated" (Sec. 2.1), without giving an explicit formula; the formula below is this
#  codebase's own reconstruction, verified level-by-level to reproduce all 19 p_stat_j values in Table 1 of the Ba paper
#  exactly (1/4 for every ^2S channel, 1/12 for ^2P, 1/20 for ^2D, 1/28 for ^2F, 1/36 for ^2G).
#
#  Known scope limitation: reflection parity (Sigma+/Sigma-) for Lambda=0 terms is NOT resolved separately -- exact for the
#  Ba+H case (both partners' relevant terms are always Sigma+ there) but unverified where it could matter.


"""
`Empirical.molecularOrbitalFraction(Lj::Int64, Lion::Int64)`
    ... to compute the fraction of channel j's own (2*Lj+1) orbital magnetic sublevels that fall on a molecular Lambda value
        also reachable by the ionic entrance term (Lion combined with a partner of L=0, i.e. Lambda_ion = 0,1,...,Lion). Since
        both Lj and Lion are combined with an L=0 partner (H or H^-), their own Lambda ranges are 0,1,...,Lj and 0,1,...,Lion
        respectively, so the shared range is 0,1,...,min(Lj,Lion) -- always including Lambda=0. A value::Float64 in (0,1] is
        returned.
"""
function molecularOrbitalFraction(Lj::Int64, Lion::Int64)
    sharedMultiplicity = 1.0 + 2.0 * min(Lj, Lion)
    return( sharedMultiplicity / (2*Lj + 1) )
end


"""
`Empirical.molecularSpinFraction(Sj::Float64, Sion::Float64)`
    ... to compute the fraction of channel j's own 2*(2*Sj+1) spin magnetic sublevels (from combining the atomic spin Sj with
        neutral hydrogen's S=1/2) that fall on the molecular spin S_mol = Sion required by the ionic entrance term (combining
        the ionic term's own spin Sion with H^-'s S=0, so the ionic side has the single, unsplit value S_mol = Sion). Returns
        0.0 if Sion is not among {Sj+1/2, Sj-1/2} (no one-electron-transfer path conserves spin between the two terms). A
        value::Float64 in [0,1] is returned.
"""
function molecularSpinFraction(Sj::Float64, Sion::Float64)
    total = 2.0 * (2*Sj + 1)
    if      Sion == Sj + 0.5                  matched = 2*Sion + 1
    elseif  Sj > 0.0   &&   Sion == Sj - 0.5   matched = 2*Sion + 1
    else                                       matched = 0.0
    end
    return( matched / total )
end


"""
`Empirical.statisticalWeight(Lj::Int64, Sj::Float64, Lion::Int64, Sion::Float64)`
    ... to compute the statistical probability p_stat_j (Belyaev & Yakovleva 2017, Eqs. 1-2) that channel j = A^Z+(j) + H
        correlates, by molecular-term symmetry, to the ionic entrance state A^(Z+1)+ + H^- of term (Lion,Sion), as the product
        of Empirical.molecularOrbitalFraction and Empirical.molecularSpinFraction. Returns 0.0 (no correlation, channel j is
        not accessible from this ionic entrance) if the spin fraction is zero; see the section note above for the derivation
        and its verification against Table 1 of the paper. A value::Float64 in [0,1) is returned.

        SCOPE. The counting is carried out in L and S, so it presumes LS coupling on both sides. That is genuinely valid for
        the one-valence-electron systems the model was formulated on -- Ba+ in the paper -- but NOT for an open-shell or
        heavy ion in intermediate coupling, where L and S are not good quantum numbers while J and parity still are. There
        the weights returned here are approximate in a way this function cannot detect or report. A J-primary formulation,
        combining J_j with H's fixed J = 1/2 by the same Omega-projection counting, would remove the assumption; it should
        reproduce these numbers wherever LS coupling is in fact valid.
"""
function statisticalWeight(Lj::Int64, Sj::Float64, Lion::Int64, Sion::Float64)
    return( Empirical.molecularOrbitalFraction(Lj, Lion) * Empirical.molecularSpinFraction(Sj, Sion) )
end


"""
`struct  Empirical.InelasticHChannel`
    ... a single atomic/ionic scattering channel A^Z+(j) + H that correlates to the ionic ground molecular state
        A^(Z+1)+ + H^- (Belyaev & Yakovleva 2017, Steps 1-2); typically built via Empirical.statisticalWeight for the p_stat
        field, or -- for a new reaction specified via Configuration -- automatically by Empirical.InelasticHReaction below.

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

        Also available as Empirical.inelasticHCollisionRateMatrix(T, reaction::Empirical.InelasticHReaction, energies; ...)
        below, which builds channels automatically from Nuclear.Model/Configuration input instead of requiring the caller
        to hand-build the Empirical.InelasticHChannel list.
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


#################################################################################################################################
### Configuration-level reaction interface (Empirical.InelasticHReaction) #######################################################
#
#  Specifies the reaction directly in terms of JAC's own Nuclear.Model and Configuration entities, exactly the way any other
#  JAC structure calculation is specified -- reduced mass, molecular symmetry, and statistical weights are all DERIVED from
#  this input (via the sections above), not supplied by the caller. Level ENERGIES remain the one piece this module does not
#  compute itself (see the module note at the top of this file for why) -- always supplied by the caller as an
#  Array{Pair{Configuration,Float64},1} of total energies (a plain array, not a Dict, because Configuration currently has no
#  matching `hash` method for its working `==` -- a real, separate JAC bug found while building this interface, not fixed here
#  since Configuration is a foundational type used throughout JAC).
#
#  Scope (deliberately restricted, matching the physics already validated for Ba2+ + H- -> Ba+ + H): the entrance ion must be
#  closed-shell (so its own molecular symmetry is the trivial, unique 1S0 -- Lion=0, Sion=0.0, no term-generation needed), and
#  every final ion configuration must differ from the entrance configuration by exactly ONE electron in exactly ONE shell (a
#  genuine single-active-electron transfer, S=1/2 always). Configurations that don't fit this pattern raise an informative
#  error rather than being silently mishandled or guessed at; open-shell entrance ions would need real term generation from
#  the configuration, not yet implemented here.


"""
`struct  Empirical.InelasticHReaction`
    ... specifies an inelastic ion + hydrogen reaction A^(Z+1)+ + H^- <-> A^Z+(f) + H directly in terms of JAC's own Nuclear.Model
        and Configuration entities, for one or several final configurations at once.

    + nm        ::Nuclear.Model            ... the projectile ion's nuclear model; only nm.Z is used here (for the covalent
                                                species' charge and for looking up its standard atomic mass).
    + iConfIon  ::Configuration             ... the (closed-shell) ground configuration of the ionic entrance species A^(Z+1)+.
    + iConfH    ::Configuration             ... either H^-(1s^2) or H(1s^1) -- this module only ever considers hydrogen.
    + fConfIon  ::Array{Configuration,1}    ... the final ion configurations A^Z+(f) of interest; each must differ from iConfIon
                                                by exactly one electron in exactly one shell.
    + fConfH    ::Configuration             ... the complementary final hydrogen configuration; derived automatically from
                                                iConfH by electron-count conservation if not supplied explicitly.
"""
struct  InelasticHReaction
    nm          ::Nuclear.Model
    iConfIon    ::Configuration
    iConfH      ::Configuration
    fConfIon    ::Array{Configuration,1}
    fConfH      ::Configuration
end


"""
`Empirical.InelasticHReaction(nm::Nuclear.Model, iConfIon::Configuration, iConfH::Configuration,
                              fConfIon::Array{Configuration,1}; fConfH::Union{Nothing,Configuration}=nothing)`
    ... outer constructor that validates iConfH/fConfH and derives fConfH from iConfH by electron-count conservation if not
        given explicitly. Raises an informative error if iConfH is not H(1s^1) or H^-(1s^2), or if a given fConfH is
        inconsistent with iConfH's electron count. An Empirical.InelasticHReaction is returned.
"""
function InelasticHReaction(nm::Nuclear.Model, iConfIon::Configuration, iConfH::Configuration,
                            fConfIon::Array{Configuration,1}; fConfH::Union{Nothing,Configuration}=nothing)
    Empirical.validateHydrogenConfiguration(iConfH)
    if  isnothing(fConfH)
        fConfHx = iConfH.NoElectrons == 2 ? Configuration("1s^1") : Configuration("1s^2")
    else
        Empirical.validateHydrogenConfiguration(fConfH)
        fConfHx = fConfH
    end
    if  abs(fConfHx.NoElectrons - iConfH.NoElectrons) != 1
        error("Empirical.InelasticHReaction: iConfH ($iConfH) and fConfH ($fConfHx) must differ by exactly one electron.")
    end
    return( InelasticHReaction(nm, iConfIon, iConfH, fConfIon, fConfHx) )
end


"""
`Empirical.validateHydrogenConfiguration(conf::Configuration)`
    ... to check that conf is either H(1s^1) or H^-(1s^2); raises an informative error otherwise, since
        Empirical.InelasticHReaction is restricted to hydrogen collision partners by design (see the section note above).
"""
function validateHydrogenConfiguration(conf::Configuration)
    if  conf != Configuration("1s^1")  &&  conf != Configuration("1s^2")
        error("Empirical.InelasticHReaction is restricted to hydrogen collision partners: expected H(1s^1) or H^-(1s^2), " *
              "got $conf.")
    end
end


"""
`Empirical.isClosedShell(conf::Configuration)`
    ... to check whether every shell of conf is fully occupied (occupation = 2*(2l+1) for each shell). A Bool is returned.
"""
function isClosedShell(conf::Configuration)
    for  (sh, occ)  in  conf.shells
        if  occ != 2*(2*sh.l + 1)   return( false )   end
    end
    return( true )
end


"""
`Empirical.activeShell(iConf::Configuration, fConf::Configuration)`
    ... to identify the single shell whose occupation differs between iConf and fConf, and check that it differs by exactly
        one electron (a genuine one-electron transfer) -- the scope Empirical.InelasticHReaction is restricted to (see the
        section note above). Raises an informative error if more than one shell differs, or if the occupation change is not
        exactly one electron. A sh::Shell is returned.
"""
function activeShell(iConf::Configuration, fConf::Configuration)
    allShells  = union(keys(iConf.shells), keys(fConf.shells))
    diffShells = Shell[]
    for  sh  in  allShells
        if  get(iConf.shells, sh, 0) != get(fConf.shells, sh, 0)   push!(diffShells, sh)   end
    end
    if  length(diffShells) != 1
        error("Empirical.InelasticHReaction (this version) requires iConfIon and each fConfIon to differ in exactly one " *
              "shell (a single-active-electron transfer); found $(length(diffShells)) differing shells between " *
              "$iConf and $fConf.")
    end
    sh   = diffShells[1]
    dOcc = get(fConf.shells, sh, 0) - get(iConf.shells, sh, 0)
    if  dOcc != 1
        error("Empirical.InelasticHReaction (this version) requires the active shell to GAIN exactly one electron; " *
              "shell $sh changes by $dOcc electrons between $iConf and $fConf.")
    end
    return( sh )
end


"""
`Empirical.reducedMassH(nm::Nuclear.Model; ionMass::Union{Nothing,Float64}=nothing)`
    ... to compute the reduced mass of the projectile ion (nuclear charge nm.Z) with a hydrogen atom, in electron-mass atomic
        units. The projectile's mass is taken from PeriodicTable.getData("mass", nm.Z) (the standard atomic weight, in amu) by
        default; note that nm.mass itself is NOT usable here -- it is a generic A ~ 2Z+0.005Z^2 formula used only to derive the
        nuclear radius, not a real isotope mass. Pass ionMass [amu] explicitly to use a specific isotope instead. A
        value::Float64 [a.u.] is returned.
        Note: this rate model is thermally averaged over energy and summed over many partial waves, which washes out most of
              the mass-dependence of any single collision; even a +/-20% mass error changes typical rate coefficients by well
              under 1% (checked explicitly for Ba+ + H). Getting the mass exactly right is far less important than getting the
              level energies right for this model.
"""
function reducedMassH(nm::Nuclear.Model; ionMass::Union{Nothing,Float64}=nothing)
    Mion_amu = isnothing(ionMass) ? PeriodicTable.getData("mass", round(Int64, nm.Z)) : ionMass
    Mion = Mion_amu / Defaults.ELECTRON_MASS_U
    MH   = PeriodicTable.getData("mass", 1) / Defaults.ELECTRON_MASS_U
    return( Mion * MH / (Mion + MH) )
end


"""
`Empirical.energyOf(energies::Array{Pair{Configuration,Float64},1}, conf::Configuration)`
    ... to look up conf's total energy [a.u.] in energies by VALUE equality (Configuration == is defined and reliable, but
        Configuration currently has no matching `hash` method, so a genuine Dict{Configuration,Float64} silently drops valid
        keys -- a real, separate JAC bug worth fixing centrally at some point, not something to route around by relying on
        Dict here). Raises an informative error if conf is not found. A value::Float64 is returned.
"""
function energyOf(energies::Array{Pair{Configuration,Float64},1}, conf::Configuration)
    idx = findfirst(p -> p.first == conf, energies)
    if  isnothing(idx)   error("energies has no entry for configuration $conf.")   end
    return( energies[idx].second )
end


"""
`Empirical.inelasticHCollisionRateMatrix(T::Float64, reaction::Empirical.InelasticHReaction,
                                         energies::Array{Pair{Configuration,Float64},1};
                                         energyLabel::String="externally supplied", printout::Bool=false, zerosGL::Int64=64)`
    ... convenience method that builds the Array{Empirical.InelasticHChannel,1} directly from reaction -- reduced mass via
        Empirical.reducedMassH, molecular symmetry/statistical weight via Empirical.activeShell + Empirical.statisticalWeight
        (assuming a closed-shell, Lion=0/Sion=0.0, entrance ion) -- and dispatches to
        Empirical.inelasticHCollisionRateMatrix(T, channels, pstatIonic, Z, mu; ...). energies must contain a total-energy
        entry [a.u.] (via Empirical.energyOf, e.g. `[Configuration("[Xe]") => 0.0, Configuration("[Xe] 6s^1") => -0.368]`)
        for reaction.iConfIon and for every reaction.fConfIon; this module never computes these itself (see the module note
        at the top of this file) -- energyLabel is a short caller-supplied description (e.g. "single-configuration Dirac-Fock
        SCF" or "NIST") echoed in the printout so the origin of the numbers is never silently lost. A named tuple, as returned
        by the low-level method, is returned.
"""
function inelasticHCollisionRateMatrix(T::Float64, reaction::Empirical.InelasticHReaction,
                                       energies::Array{Pair{Configuration,Float64},1};
                                       energyLabel::String="externally supplied", printout::Bool=false, zerosGL::Int64=64)
    if  !Empirical.isClosedShell(reaction.iConfIon)
        error("Empirical.InelasticHReaction (this version) requires a closed-shell entrance ion configuration; " *
              "$(reaction.iConfIon) is not closed-shell. An open-shell entrance needs real term generation, not yet " *
              "implemented here.")
    end
    Zcov = reaction.nm.Z - reaction.fConfIon[1].NoElectrons
    for  fConf  in  reaction.fConfIon
        if  reaction.nm.Z - fConf.NoElectrons != Zcov
            error("All entries of reaction.fConfIon must have the same electron count; $(fConf) does not match the others.")
        end
    end

    mu    = Empirical.reducedMassH(reaction.nm)
    E_ion = Empirical.energyOf(energies, reaction.iConfIon)

    channels = Empirical.InelasticHChannel[]
    for  fConf  in  reaction.fConfIon
        Ej    = Empirical.energyOf(energies, fConf) - E_ion
        sh    = Empirical.activeShell(reaction.iConfIon, fConf)
        pstat = Empirical.statisticalWeight(sh.l, 0.5, 0, 0.0)
        push!(channels, Empirical.InelasticHChannel(string(fConf), Ej, pstat))
    end

    if  printout
        println("\n* Empirical.InelasticHReaction:  A^$(round(Int64,Zcov+1))+ ($(reaction.iConfIon)) + $(reaction.iConfH)  " *
                "->  A^$(round(Int64,Zcov))+ (f) + $(reaction.fConfH):" *
                "\n    + Energy source: $energyLabel -- not computed by this module; supplied by the caller via `energies`." *
                "\n    + Reduced mass: $(round(mu,digits=1)) [a.u.], from PeriodicTable's standard atomic weight for Z=" *
                "$(round(Int64,reaction.nm.Z)) (nm.mass itself is a nuclear-radius placeholder, not used here); this rate " *
                "model is insensitive to mass to well under 1% even for a 20% mass error, so this choice is not critical." *
                "\n    + Statistical weights assume a closed-shell (1S0) entrance ion and a single active electron per " *
                "final configuration -- exact for this scope, not a general open-shell treatment." *
                "\n    + Final-state channels and their derived (energy, p_stat):")
        for  (fConf, ch)  in  zip(reaction.fConfIon, channels)
            println("        $fConf:  E = $(round(Defaults.convertUnits("energy: from atomic to eV",ch.E),digits=4)) eV   " *
                    "p_stat = $(round(ch.pstat,digits=6))")
        end
    end

    return( Empirical.inelasticHCollisionRateMatrix(T, channels, 1.0, Zcov, mu; printout=printout, zerosGL=zerosGL) )
end
