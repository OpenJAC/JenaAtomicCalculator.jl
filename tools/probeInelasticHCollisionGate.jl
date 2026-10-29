#
#  [S9-GATE]  A validity gate on the ionic-covalent crossing radius, and a zero statistical weight forcing a
#             zero rate -- the two cheap defects of the Belyaev-Yakovleva simplified model as implemented in
#             `src/module-Empirical-inc-inelastic-h-collisions.jl`.
#
#  THE TWO QUESTIONS THIS PROBE ANSWERS, and they are separate:
#
#    Q1 (CALIBRATION).  The model is a LONG-RANGE model: it assumes a single, isolated, diabatic ionic-covalent
#       crossing in the asymptotic region where the two partners' electron clouds do not yet overlap.  Nothing in
#       the code says so, and `nonadiabaticRadius` answers for R_c = 2 a0 as confidently as for R_c = 26 a0.  At
#       WHICH R_c does the model stop being applicable?  The threshold is MEASURED here, from two independent
#       directions, and not chosen:
#         (a) FROM BELOW, by the published quantum calculations: where both an asymptotic-LCAO/multichannel
#             Landau-Zener treatment and an MRCI/quantum-probability-current treatment give EXACTLY ZERO for a
#             channel that this model answers for, the model is outside its range.
#         (b) FROM ABOVE, by the model's OWN validated worked example: the deepest channel of Belyaev &
#             Yakovleva's own Table 1 fixes an R_c that must NOT be refused, or the paper's own example breaks.
#       A threshold must sit between (a) and (b), and the gap between them is what makes the gate defensible.
#
#    Q2 (BASELINE).  Every public entry point of the module, at full precision, over BOTH systems, so that the
#       change can be verified to move nothing except where a number should become zero.  The six dated example
#       branches of `examples/example-Ni.jl` and `examples/example-Nl.jl` print rounded to 4 significant digits;
#       this prints all of them, which is a strictly stronger baseline.
#
#  RUN:   julia --project=. tools/probeInelasticHCollisionGate.jl
#
using JenaAtomicCalculator, Printf
const JAC = JenaAtomicCalculator

eV      = x -> Defaults.convertUnits("energy: from eV to atomic", x)
toEv    = x -> Defaults.convertUnits("energy: from atomic to eV", x)
toCm3s  = Defaults.convertUnits("length: from atomic to cm", 1.0)^3 / Defaults.convertUnits("time: from atomic to sec", 1.0)
EHminus = Empirical.hydrogenAnionEnergy()
T6000   = Defaults.convertUnits("temperature: from Kelvin to (Hartree) units", 6000.)

##  SYSTEM 1: Ba2+ + H^- -> Ba+(nl) + H, Z = 1.  Belyaev & Yakovleva, A&A 608, A33 (2017), Table 1 -- the model's
##  OWN worked example and the case `neutralizationReducedRate` is validated against to 3 significant figures.
MBa    = 137.327 / Defaults.ELECTRON_MASS_U;   MH = 1.00794 / Defaults.ELECTRON_MASS_U
muBa   = MBa * MH / (MBa + MH)
baTable = [ ("Ba+(6s)",  -10.0080, 1/4 ), ("Ba+(5d)",  -9.34416, 1/20), ("Ba+(6p)",  -7.35615, 1/12),
            ("Ba+(7s)",  -4.75665, 1/4 ), ("Ba+(6d)",  -4.29573, 1/20), ("Ba+(4f)",  -4.00879, 1/28),
            ("Ba+(7p)",  -3.83309, 1/12), ("Ba+(5f)",  -2.87539, 1/28), ("Ba+(8s)",  -2.81381, 1/4 ),
            ("Ba+(7d)",  -2.58669, 1/20), ("Ba+(8p)",  -2.37789, 1/12), ("Ba+(5g)",  -2.19370, 1/36),
            ("Ba+(6f)",  -1.99196, 1/28), ("Ba+(9s)",  -1.86429, 1/4 ), ("Ba+(8d)",  -1.73768, 1/20),
            ("Ba+(9p)",  -1.62363, 1/12), ("Ba+(6g)",  -1.52427, 1/36), ("Ba+(7f)",  -1.42341, 1/28),
            ("Ba+(10s)", -1.32732, 1/4 ) ]

##  SYSTEM 2: O+ (2p3 4So) + H^- -> O(nl) + H, Z = 0.  NIST term energies as tabulated by Barklem's own public
##  archive, so that all three columns rest on the same experimental level structure.  IP(O) = 13.618055 eV.
##  The two published columns: Belyaev et al., MNRAS 487, 5097 (2019) [MRCI + quantum probability current] and
##  Barklem, A&A 610, A57 (2018) [asymptotic LCAO + multichannel Landau-Zener];  -1 = not tabulated there.
IPoxygen = 13.618055;   cm2eV = 1.0 / 8065.543937
muO      = Empirical.reducedMassH(Nuclear.Model(8.))
##        label                 E[cm^-1]    L   S     Belyaev19    Barklem18
oTable = [("O(2p4 3P)   + H",       0.000,  1, 1.0,     0.0,         0.0      ),
          ("O(2p4 1D)   + H",   15790.000,  2, 0.0,     0.0,         0.0      ),
          ("O(2p4 1S)   + H",   33715.000,  0, 0.0,     0.0,         0.0      ),
          ("O(3s 5So)   + H",   73690.000,  0, 2.0,     2.15e-13,    1.550e-14),
          ("O(3s 3So)   + H",   76717.000,  0, 1.0,     2.68e-10,    8.600e-12),
          ("O(3p 5P)    + H",   86551.000,  1, 2.0,     1.73e-09,    3.650e-09),
          ("O(3p 3P)    + H",   88553.000,  1, 1.0,     1.98e-08,    2.340e-08),
          ("O(4s 5So)   + H",   95399.000,  0, 2.0,     2.95e-08,    3.200e-08),
          ("O(4s 3So)   + H",   96147.000,  0, 1.0,     8.48e-09,    1.330e-08),
          ("O(3d 5Do)   + H",   97343.000,  2, 2.0,    -1.0,         6.850e-09),
          ("O(3d 3Do)   + H",   97411.000,  2, 1.0,    -1.0,         3.950e-09),
          ("O(4p 5P)    + H",   99016.000,  1, 2.0,    -1.0,         4.650e-10),
          ("O(4p 3P)    + H",   99603.000,  1, 1.0,    -1.0,         3.650e-11),
          ("O(3s 3Do)   + H",  101065.000,  2, 1.0,    -1.0,         0.0      ),
          ("O(5s 5So)   + H",  102039.000,  0, 2.0,    -1.0,         0.0      ),
          ("O(5s 3So)   + H",  102334.000,  0, 1.0,    -1.0,         0.0      ),
          ("O(3s 1Do)   + H",  102584.000,  2, 0.0,    -1.0,         0.0      )]

println("\n", "="^125)
println("[S9-GATE-Q1]  CALIBRATION.  Where does the crossing stop being a LONG-RANGE crossing?")
println("="^125)

println("\n(a) FROM BELOW -- O+ + H^-, Z = 0, T = 6000 K.  R_c = 1/(E_H- - E_f).  Both published columns are")
println("    EXACTLY ZERO for the three deeply bound terms, where this model returns a large number.")
println("\n    state                E_f[eV]    R_c[a0]    p_stat    JAC[cm^3/s]   Belyaev19     Barklem18")
for  (k, st)  in  enumerate(oTable)
    Ef = eV(st[2]*cm2eV - IPoxygen)
    Rc = Empirical.nonadiabaticRadius(Ef, EHminus, 0.)
    p  = Empirical.statisticalWeight(st[3], st[4], 0, 1.5)
    N  = Empirical.neutralizationReducedRate(T6000, Ef, EHminus, 0., muO) * toCm3s
    @printf("    %-18s %8.4f  %8.2f  %8.5f    %.4e    %.3e     %.3e\n", st[1], toEv(Ef), Rc, p, N, st[5], st[6])
end

println("\n(b) FROM ABOVE -- Ba2+ + H^-, Z = 1, the paper's OWN Table 1.  R_c = 2/(E_H- - E_f).  Every one of")
println("    these channels is tabulated BY THE PAPER, so no gate may refuse them.")
println("\n    state           E_f[eV]    R_c[a0]    N[cm^3/s]")
RcMinBa = Inf;   RcMinBaName = ""
for  (name, EfeV, pst)  in  baTable
    Ef = eV(EfeV);   Rc = Empirical.nonadiabaticRadius(Ef, EHminus, 1.)
    N  = Empirical.neutralizationReducedRate(T6000, Ef, EHminus, 1., muBa) * toCm3s
    global RcMinBa, RcMinBaName
    if  Rc < RcMinBa   RcMinBa = Rc;   RcMinBaName = name   end
    @printf("    %-12s %9.4f  %8.2f    %.4e\n", name, EfeV, Rc, N)
end

RcMaxZero = -Inf;   RcMaxZeroName = ""
RcMinNonzero = Inf; RcMinNonzeroName = ""
for  st  in  oTable
    global RcMaxZero, RcMaxZeroName, RcMinNonzero, RcMinNonzeroName
    Ef = eV(st[2]*cm2eV - IPoxygen);   Rc = Empirical.nonadiabaticRadius(Ef, EHminus, 0.)
    if  st[6] == 0.0  &&  Rc < 10.0  &&  Rc > RcMaxZero      RcMaxZero = Rc;      RcMaxZeroName = st[1]      end
    if  st[6] >  0.0  &&  Rc < RcMinNonzero                  RcMinNonzero = Rc;   RcMinNonzeroName = st[1]   end
end

println("\n", "-"^125)
@printf("    DEEPEST channel the PAPER ITSELF tabulates:                 %-18s  R_c = %6.2f a0\n", RcMinBaName, RcMinBa)
@printf("    DEEPEST O channel either published calculation gives > 0:   %-18s  R_c = %6.2f a0\n", RcMinNonzeroName, RcMinNonzero)
@printf("    SHALLOWEST O channel both published calculations ZERO:      %-18s  R_c = %6.2f a0\n", RcMaxZeroName, RcMaxZero)
println("-"^125)
println("    => a REFUSE threshold must lie strictly between the last line and the two above it;")
println("       a WARN threshold may lie above it, provided it does not refuse anything.")

println("\n", "="^125)
println("[S9-GATE-Q2]  BASELINE.  Every public entry point, full precision, both systems.")
println("="^125)

println("\n-- Empirical.hydrogenAnionEnergy()")
@printf("   %.17e a.u.  = %.17e eV\n", EHminus, toEv(EHminus))

println("\n-- Empirical.nonadiabaticRadius / ionicCovalentCoupling / landauZenerProbability  (Ba, v = 1.0e-3 a.u.)")
for  (name, EfeV, pst)  in  baTable[1:6]
    Ef  = eV(EfeV)
    Rc  = Empirical.nonadiabaticRadius(Ef, EHminus, 1.)
    Hij = Empirical.ionicCovalentCoupling(Ef, EHminus, Rc)
    p   = Empirical.landauZenerProbability(Ef, EHminus, 1., 1.0e-3)
    @printf("   %-12s Rc = %.17e   Hij = %.17e   pLZ = %.17e\n", name, Rc, Hij, p)
end

println("\n-- Empirical.maxOrbitalAngularMomentum")
for  (ae, R)  in  [(1.0e-4, 10.0), (1.0e-3, 26.0), (1.0e-2, 5.0), (-1.0e-5, 10.0)]
    @printf("   availableEnergy = %10.3e   R = %6.2f   ->  Jmax = %d\n", ae, R,
            Empirical.maxOrbitalAngularMomentum(ae, R, muBa))
end

println("\n-- Empirical.neutralizationCrossSection  (Ba+(6d), the paper's largest-rate channel)")
for  E  in  [1.0e-4, 5.0e-4, 1.0e-3, 5.0e-3]
    @printf("   E = %10.3e  ->  sigma = %.17e\n", E,
            Empirical.neutralizationCrossSection(E, eV(-4.29573), EHminus, 1., muBa))
end

println("\n-- Empirical.neutralizationReducedRate  (the WHOLE Ba table; the paper quotes 7.59e-8 for Ba+(6d))")
for  (name, EfeV, pst)  in  baTable
    N = Empirical.neutralizationReducedRate(T6000, eV(EfeV), EHminus, 1., muBa)
    @printf("   %-12s N = %.17e a.u.  = %.6e cm^3/s\n", name, N, N*toCm3s)
end

println("\n-- Empirical.neutralizationReducedRate  (the WHOLE O table)")
for  st  in  oTable
    N = Empirical.neutralizationReducedRate(T6000, eV(st[2]*cm2eV - IPoxygen), EHminus, 0., muO)
    @printf("   %-18s N = %.17e a.u.  = %.6e cm^3/s\n", st[1], N, N*toCm3s)
end

println("\n-- Empirical.deExcitationCrossSection / deExcitationReducedRate  (Ba, 7s -> 6s and 6d -> 6p)")
for  (iN, iE, fN, fE)  in  [("Ba+(7s)", -4.75665, "Ba+(6s)", -10.0080), ("Ba+(6d)", -4.29573, "Ba+(6p)", -7.35615)]
    sig = Empirical.deExcitationCrossSection(1.0e-3, eV(iE), eV(fE), EHminus, 1., muBa)
    D   = Empirical.deExcitationReducedRate(T6000, eV(iE), eV(fE), EHminus, 1., muBa)
    @printf("   %-10s -> %-10s  sigma(1e-3) = %.17e   D = %.17e  = %.6e cm^3/s\n", iN, fN, sig, D, D*toCm3s)
end

println("\n-- Empirical.detailedBalanceRate")
for  (K, pJ, pK, dE)  in  [(1.0e-10, 1.0, 0.25, 0.3), (1.0e-10, 1.0, 1/36, 0.05)]
    @printf("   K = %9.3e  pJ = %.5f  pK = %.5f  dE = %.3f  ->  %.17e\n", K, pJ, pK, dE,
            Empirical.detailedBalanceRate(K, pJ, pK, dE, T6000))
end

println("\n-- Empirical.molecularOrbitalFraction / molecularSpinFraction / statisticalWeight  (O+ 4So: L=0, S=3/2)")
for  st  in  oTable
    fL = Empirical.molecularOrbitalFraction(st[3], 0)
    fS = Empirical.molecularSpinFraction(st[4], 1.5)
    p  = Empirical.statisticalWeight(st[3], st[4], 0, 1.5)
    @printf("   %-18s L=%d S=%.1f   fL = %.17e   fS = %.17e   pstat = %.17e\n", st[1], st[3], st[4], fL, fS, p)
end

println("\n-- Empirical.reducedMassH")
@printf("   Nuclear.Model(8.)  -> mu = %.17e m_e\n", Empirical.reducedMassH(Nuclear.Model(8.)))
@printf("   Nuclear.Model(56.) -> mu = %.17e m_e\n", Empirical.reducedMassH(Nuclear.Model(56.)))

println("\n-- Empirical.inelasticHCollisionRateMatrix  (the full Ba matrix, ionic<->channel legs)")
baChannels = [Empirical.InelasticHChannel(n, eV(e), p)  for (n, e, p) in baTable]
resBa   = Empirical.inelasticHCollisionRateMatrix(T6000, baChannels, 1.0, 1., muBa)
ionToCh = resBa.ionicToChannel;   chToIon = resBa.channelToIonic;   chToCh = resBa.channelToChannel
for  (k, ch)  in  enumerate(baChannels)
    @printf("   %-12s ionic->j = %.17e   j->ionic = %.17e\n", ch.name, ionToCh[k], chToIon[k])
end
@printf("   channel->channel block: sum = %.17e   max = %.17e\n", sum(chToCh), maximum(chToCh))

println("\n-- Empirical.inelasticHCollisionRateMatrix  (O, with THREE channels of pstat = 0 present)")
oChannels = [Empirical.InelasticHChannel(st[1], eV(st[2]*cm2eV - IPoxygen),
                                         Empirical.statisticalWeight(st[3], st[4], 0, 1.5))  for st in oTable]
resO     = Empirical.inelasticHCollisionRateMatrix(T6000, oChannels, 1.0, 0., muO)
ionToChO = resO.ionicToChannel;   chToIonO = resO.channelToIonic;   chToChO = resO.channelToChannel
for  (k, ch)  in  enumerate(oChannels)
    @printf("   %-18s pstat = %.5f   ionic->j = %.17e   j->ionic = %.17e\n",
            ch.name, ch.pstat, ionToChO[k], chToIonO[k])
end
@printf("   channel->channel block: sum = %.17e   max = %.17e\n", sum(chToChO), maximum(chToChO))

println("\n[S9-GATE]  done")
