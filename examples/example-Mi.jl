println("Mi) Test  computeProperties()  for the hyperfine A,B constants and the Lande g_J factors.")
#
# The hyperfine constants scale with the NUCLEAR MOMENTS and not with the charge: A with mu/I and B with Q.
# Neither can be inferred from Z, so every scenario below passes an explicit Nuclear.Model.  A model built as
# Nuclear.Model(Z) carries mu = Q = 0 and would return zeros; scenario 3 shows that on purpose.

setDefaults("unit: energy", "Hz")
setDefaults("print summary: open", "zzz-ForPedestrians.sum")


if  false
    # Last visit:  18Aug2026
    # Last successful:  unknown ... 23Na 3s ground level; A is the textbook case, 885.8 MHz for the 3s_1/2.
    setDefaults("nuclear: charge", 11.0)
    configs = [Configuration("[Ne] 3s")]
    nm      = Nuclear.Model(11., Nuclear.UniformNucleus(), 23., 2.98, AngularJ64(3//2), 2.2176, 0.10, 0.0)
    computeProperties(Basics.HyperfineStructure(AngularJ64(3//2)), configs, nuclearModel=nm)
    #
elseif  true
    # Last successful:  18Aug2026
    # 1H 1s: the hyperfine splitting is the 21 cm line, A = 1420.405751 MHz.  JAC gives 1421.22 MHz, i.e.
    # 0.057 % -- and this is the case that shows why the field matters: with a screened (DFS) potential the
    # same calculation gives 1464.20 MHz, 3.1 % high, because a lone electron would be screening against
    # itself.  computeProperties picks Basics.NuclearField() for a one-electron system on its own.
    # The uniform nuclear model is required at Z = 1: JAC's two-parameter Fermi distribution cannot represent
    # an rms radius below ~1.86 fm, far above hydrogen's 0.88 fm (see examples/example-Cb.jl branch a).
    setDefaults("nuclear: charge", 1.0)
    configs = [Configuration("1s")]
    nm      = Nuclear.Model(1., Nuclear.UniformNucleus(), 1., 0.8797, AngularJ64(1//2), 2.7928, 0.0, 0.0)
    computeProperties(Basics.HyperfineStructure(AngularJ64(1//2)), configs, nuclearModel=nm)
    #
elseif  false
    # Last visit:  18Aug2026
    # Last successful:  unknown ... the same sodium case WITHOUT nuclear moments, to show what the note says:
    # Nuclear.Model(Z) carries mu = Q = 0, so every constant comes out zero by construction rather than by error.
    setDefaults("nuclear: charge", 11.0)
    configs = [Configuration("[Ne] 3s")]
    computeProperties(Basics.HyperfineStructure(AngularJ64(3//2)), configs)
    #
elseif  false
    # Last successful:  18Aug2026
    # Ge II (Z=32) 4s^2 4f: the Lande g_J factors of the (2)F_5/2,7/2 pair, against Andersson & Jonsson,
    # CPC 178 (2008) 156.  JAC gives 1.14318192 and 0.856804697 against their 1.143182 and 0.856804 -- six
    # decimals on both, and with NO grid given: the box is matched by Basics.recommendedGrid.
    # This is the case that shows why that matters.  On Radial.Grid(true), a 614 a.u. box, the same
    # calculation returns g_J(4f_7/2) = -2.263670 against the exact 8/7 -- wrong by a factor and in sign --
    # because the SCF converges onto the WRONG STATE for kappa = -4.  It was recorded for months as an
    # angular-coefficient bug; it was the radial box (see examples/example-Cd.jl branch c).
    setDefaults("nuclear: charge", 32.0)
    configs = [Configuration("[Ar] 3d^10 4s^2 4f")]
    computeProperties(Basics.ZeemanStructure(0.0), configs)
    #
elseif  false
    # Last successful:  18Aug2026
    # He 1s2p, the SECOND independent test case of Andersson & Jonsson (2008), CPC 178, 156 -- a two-electron
    # LS-coupled multiplet rather than a single f electron, so it probes different angular machinery.
    # Neither a grid nor a nuclear model is supplied; both are defaulted by the pedestrian path.
    #     3P_1  g_J = 1.50112146   vs paper 1.5011166
    #     3P_2  g_J = 1.50112284   vs paper 1.5011183
    #     1P_1  g_J = 0.99999320   vs paper 0.9999936      -- six significant figures on all three
    #     3P_0  g_J = 0            trivially, J = 0
    setDefaults("nuclear: charge", 2.0)
    configs = [Configuration("1s 2p")]
    computeProperties(Basics.ZeemanStructure(0.0), configs)
    #
elseif  false
    # Last visit:  18Aug2026
    # Last successful:  unknown ... the same levels IN A FIELD: g_J plus the Zeeman splittings at 1 Tesla.
    # A non-zero field in the theme is what turns the splittings on.
    setDefaults("nuclear: charge", 32.0)
    configs = [Configuration("[Ar] 3d^10 4s^2 4f")]
    computeProperties(Basics.ZeemanStructure(1.0), configs)
    #
end
