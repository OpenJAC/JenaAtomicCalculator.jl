#
# Item 64, second half:  does the known DFS valence underbinding reach an OPEN f SHELL, and does it
# propagate into the quantities an f-shell calculation is actually asked for?
#
# The system is Eu^3+ [Xe] 4f^6, whose ground term is ^7F.  Three yardsticks are used, and NONE of them
# needs a literature value that is not already in this repository:
#
#   (1) THE ORBITAL ENERGY against Koopmans.  `PeriodicTable.ionizationPotentials_Nist2025(63)` carries all
#       63 successive ionization potentials of europium, so the 4th one -- Eu^3+ -> Eu^4+ -- is the energy
#       the 4f eigenvalue should approximate.  Comparing the DFS and the AL field on the SAME grid and the
#       SAME configuration is then the direct measurement of the underbinding for an open f shell.
#   (2) g_J AGAINST AN EXACT FRACTION.  A ^7F term has S = L = 3, so S(S+1) = L(L+1) = 12 and
#       g_J = 1 + [J(J+1) + S(S+1) - L(L+1)] / [2 J(J+1)] = 3/2 for EVERY J.  A flat line at 1.5 is
#       therefore a sharp test both of the angular machinery and of the term identification -- a level that
#       is not ^7F cannot sit on it.
#   (3) THE LANDE INTERVAL RULE.  For a pure LS term split by spin-orbit coupling, E(J) - E(J-1) is
#       proportional to J, so the six ^7F intervals stand in the ratio 1:2:3:4:5:6.  Departures from that
#       measure intermediate coupling, not an error.
#
# RESULTS, 19-Oct-2026, both fields on the same recommendedGrid (rbox 10.7 a.u., 616 points, 94 splines),
# single-configuration 4f^6 CI, 295 CsfR, 50 s for the whole script.  NIST ASD Eu IV supplies the targets.
#
#   ORBITAL                 DFSField              ALField
#   -eps(4f_5/2)/IP4          0.7448               1.0237
#   -eps(4f_7/2)/IP4          0.7290               1.0039        IP4 = 42.94 eV = 1.578016 Ha
#   <r>(4f_5/2)             0.87466 a.u.         0.84593 a.u.    AL is 3.3 % more compact
#   zeta(4f)                 1563.6 cm^-1         1959.1 cm^-1   from the one-electron splitting
#
#   THE TWELVE LEVELS       NIST        DFS      ratio        AL      ratio     [cm^-1 above ^7F_0]
#   ^7F_1                    370       408.3     1.103      381.6     1.031
#   ^7F_2                   1040      1136.7     1.093     1072.5     1.031
#   ^7F_3                   1890      2090.3     1.106     1989.3     1.053
#   ^7F_4                   2860      3195.5     1.117     3064.4     1.072
#   ^7F_5                   3910      4400.3     1.125     4247.9     1.086
#   ^7F_6                   4940      5667.3     1.147     5504.0     1.114
#   ^5D_0                  17270     19220.3     1.113    20433.2     1.183
#   ^5D_1                  19030     20983.0     1.103    22097.3     1.161
#   ^5D_2                  21510     23616.4     1.098    24650.2     1.146
#   ^5D_3                  24390     26798.1     1.099    27778.0     1.139
#   ^5D_4                  27640     30470.9     1.102    31444.3     1.138
#
#   g_J                    exact        DFS                 AL              AL - DFS
#   ^7F_1..^7F_6             1.5    1.498719..1.494278  1.498826..1.495047   <= 7.7e-04
#   ^5D_4  (level 21)        1.5       1.461366            1.410544          -5.1e-02
#
# WHAT THEY SAY.  The underbinding DOES reach an open f shell -- 27 % on DFS -- and the AL field removes it.
# A Lande g_J does not notice: it moves 5e-04 relative while eps(4f) moves 27 %.  The exception is the level
# that is near-degenerate with a same-J partner, ^5D_4 against ^5G_4, whose g_J moves 3.5 %.  The ^7F ladder
# is NOT a Lande ladder here, measured or computed (the measured interval ratios are 1, 1.81, 2.30, 2.62,
# 2.84, 2.78 against 1:2:3:4:5:6) and so cannot be reasoned about as a one-electron zeta observable: AL's
# zeta is 25 % LARGER than DFS's and its ladder comes out NARROWER and closer to experiment.  The ^5D - ^7F
# term separation moves the WRONG WAY on AL because the more compact 4f raises the Slater integrals, which
# a single-configuration calculation has no correlation to bring back down.
#
using JenaAtomicCalculator
const JAC = JenaAtomicCalculator
using Printf, LinearAlgebra

Defaults.setDefaults("print summary: open", "zzz-item64-Eu4f.sum")

const Z       = 63
const nm      = Nuclear.Model(Float64(Z))
const configs = [Configuration("[Xe] 4f^6")]

ipList  = PeriodicTable.ionizationPotentials_Nist2025(Z)
ip4eV   = ipList[4]
ip4Ha   = Defaults.convertUnits("energy: from eV to atomic", ip4eV)

println("\n", "="^125)
println("ITEM 64, SECOND HALF -- Eu^3+ [Xe] 4f^6:  does the DFS underbinding reach an open f shell?")
println("="^125)
@printf("NIST 2025 ionization potentials of Eu, 1st-5th:  %.4f  %.4f  %.4f  %.4f  %.4f eV\n",
        ipList[1], ipList[2], ipList[3], ipList[4], ipList[5])
@printf("The 4th, Eu^3+ -> Eu^4+, is the Koopmans target for the 4f eigenvalue:  %.4f eV = %.6f Ha\n",
        ip4eV, ip4Ha)

grid = Basics.recommendedGrid(configs, nm; printout=true)
@printf("\n>> grid:  r_max = %.3f a.u.,  %d mesh points\n", grid.r[end], grid.NoPoints)

store = Dict{String,Any}()

for (label, field) in [("DFS", Basics.DFSField()), ("AL", Basics.ALField())]
    println("\n", "-"^125)
    println(">>> MEAN FIELD = $label")
    println("-"^125)
    asf = ManyElectron.AsfSettings(ManyElectron.AsfSettings(); scField=field, jjLS=LSjjSettings(true))
    local mp
    try
        mp = SelfConsistent.performSCF(configs, nm, grid, asf; printout=true)
    catch err
        println(">>> $label FAILED:  $err")
        continue
    end
    orbs = mp.levels[1].basis.orbitals

    println("\n>> $label  orbital energies of the outer shells, against the 4th ionization potential:")
    epsStore = Dict{String,Float64}()
    for sh in [Subshell("4d_3/2"), Subshell("4d_5/2"), Subshell("5s_1/2"), Subshell("5p_1/2"),
               Subshell("5p_3/2"), Subshell("4f_5/2"), Subshell("4f_7/2")]
        if haskey(orbs, sh)
            eps = orbs[sh].energy;   epsStore[string(sh)] = eps
            @printf("   %-10s  eps = %+12.6f Ha = %+11.3f eV     -eps/IP4 = %7.4f\n", string(sh), eps,
                    Defaults.convertUnits("energy: from atomic to eV", eps), -eps/ip4Ha)
        end
    end

    # The mean radius and the one-electron spin-orbit splitting of the 4f, so that the many-electron
    # splittings below can be read against the orbital they are built from rather than guessed at.
    if  haskey(orbs, Subshell("4f_5/2"))  &&  haskey(orbs, Subshell("4f_7/2"))
        r5 = RadialIntegrals.rkDiagonal(1, orbs[Subshell("4f_5/2")], orbs[Subshell("4f_5/2")], grid)
        r7 = RadialIntegrals.rkDiagonal(1, orbs[Subshell("4f_7/2")], orbs[Subshell("4f_7/2")], grid)
        dSO   = orbs[Subshell("4f_7/2")].energy - orbs[Subshell("4f_5/2")].energy
        zeta  = 2*dSO/7                                     # E(j=l+1/2) - E(j=l-1/2) = zeta (2l+1)/2
        @printf("\n>> %s  <r>(4f_5/2) = %8.5f a.u.,  <r>(4f_7/2) = %8.5f a.u.\n", label, r5, r7)
        @printf("   one-electron spin-orbit splitting  eps(4f_7/2) - eps(4f_5/2) = %+10.6f Ha = %9.1f cm^-1\n",
                dSO, Defaults.convertUnits("energy: from atomic to Kayser", dSO))
        @printf("   implied zeta(4f) = %9.1f cm^-1,  so a pure ^7F Lande ladder would start at zeta/(2S) = %7.1f cm^-1\n",
                Defaults.convertUnits("energy: from atomic to Kayser", zeta),
                Defaults.convertUnits("energy: from atomic to Kayser", zeta/6))
    end

    try
        LSjj.expandLevelsIntoLS(mp, asf.jjLS)
    catch err
        println("   (jj-LS transformation raised:  $err)")
    end

    nSel = min(30, length(mp.levels))
    lz   = LandeZeeman.Settings(LandeZeeman.Settings(); calcLandeJ=true, includeSchwinger=true,
                                levelSelection=LevelSelection(true; indices=collect(1:nSel)))
    out  = LandeZeeman.computeOutcomes(mp, nm, grid, lz)
    store[label] = (mp, epsStore, out)
end

# ---------------------------------------------------------------------------------------------------------
# The three yardsticks, assembled
# ---------------------------------------------------------------------------------------------------------
cmInv(x) = Defaults.convertUnits("energy: from atomic to Kayser", x)

for label in ["DFS", "AL"]
    haskey(store, label)  ||  continue
    (mp, epsStore, out) = store[label]
    println("\n", "="^125)
    println("YARDSTICKS,  mean field = $label")
    println("="^125)

    @printf("(1) KOOPMANS.  IP4 = %.6f Ha.  -eps(4f_5/2)/IP4 = %6.4f,  -eps(4f_7/2)/IP4 = %6.4f\n",
            ip4Ha, -epsStore["4f_5/2"]/ip4Ha, -epsStore["4f_7/2"]/ip4Ha)

    println("\n(2)+(3)  THE LOWEST LEVELS.  The LS label of each is in the jj-LS table printed above;  levels 1-7")
    println("         are the ^7F ladder and levels 8-11 the first four ^5D.  g_J is compared with the exact")
    println("         LS value 3/2, which holds for BOTH ^7F (S = L = 3) and ^5D (S = L = 2).")
    println("     i   J^P     E - E_1 [cm^-1]        g_J        g_J - 1.5     interval [cm^-1]    ratio to the first")
    e1 = out[1].Jlevel.energy;   prevE = nothing;   firstInterval = nothing
    for (i, oc) in enumerate(out)
        J      = Basics.twice(oc.Jlevel.J) / 2
        dE     = cmInv(oc.Jlevel.energy - e1)
        intStr = "         --        ";  ratStr = "      --"
        if  !isnothing(prevE)
            interval = cmInv(oc.Jlevel.energy - prevE)
            if  isnothing(firstInterval)  &&  interval > 1.0   firstInterval = interval   end
            intStr = @sprintf("%15.1f    ", interval)
            ratStr = isnothing(firstInterval) ? "      --" : @sprintf("%8.3f", interval/firstInterval)
        end
        @printf("    %2d  %4.1f %s  %15.1f   %+10.6f   %+10.6f  %s%s\n",
                i, J, string(oc.Jlevel.parity), dE, oc.LandeJ, oc.LandeJ - 1.5, intStr, ratStr)
        prevE = oc.Jlevel.energy
    end

    # The ^7F ladder against the NIST levels of Eu IV, which are the item's own yardstick
    println("\n     THE ^7F LADDER against NIST ASD Eu IV (the 12 known levels of 4f^6):")
    println("       J    computed [cm^-1]    NIST [cm^-1]      ratio     Lande ratio   computed ratio")
    nist7F = [0.0, 370.0, 1040.0, 1890.0, 2860.0, 3910.0, 4940.0]
    i1 = 1
    for j in 1:7
        j <= length(out)  ||  break
        dE = cmInv(out[j].Jlevel.energy - out[i1].Jlevel.energy)
        rat = nist7F[j] > 0. ? dE/nist7F[j] : NaN
        if  j == 1
            @printf("      %2d  %16.1f  %14.1f        --            --            --\n", j-1, dE, nist7F[j])
        else
            dInt = cmInv(out[j].Jlevel.energy - out[j-1].Jlevel.energy)
            nInt = nist7F[j] - nist7F[j-1]
            @printf("      %2d  %16.1f  %14.1f    %8.4f      %8.3f      %8.3f\n",
                    j-1, dE, nist7F[j], rat, nInt/370.0, dInt/cmInv(out[2].Jlevel.energy - out[1].Jlevel.energy))
        end
    end
    println("       the last two columns are the INTERVAL ratios;  a pure ^7F term split by spin-orbit coupling")
    println("       would give 1, 2, 3, 4, 5, 6 in both -- the measured ladder does NOT, and that is physics.")

    # The ^5D_0 position, i.e. the term separation, which is an electrostatic rather than a spin-orbit quantity
    if  length(out) >= 8
        d5D0 = cmInv(out[8].Jlevel.energy - out[1].Jlevel.energy)
        @printf("\n     ^5D_0 - ^7F_0  =  %10.1f cm^-1   against the NIST 17270 cm^-1   (ratio %6.4f)\n",
                d5D0, d5D0/17270.0)
        println("       This is a TERM separation, set by the f-shell Slater integrals rather than by spin-orbit")
        println("       coupling, so it probes a different part of the same calculation.")
    end
end

println("\n", "="^125)
println("DONE")
println("="^125)
