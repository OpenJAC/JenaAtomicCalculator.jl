
println("Mg) Test  computeResonanceStrength(ForDielectronicRecombination(), ...)  for several ions.")

setDefaults("unit: energy",   "eV")
setDefaults("unit: strength", "cm^2 eV")
setDefaults("print summary: open", "zzz-ForPedestrians.sum")


if  false
    # Last successful:  27Jun2026
    # DR resonance strengths for H-like Ne^9+ with 1s --> 2p core excitation (KL3 series).
    # Initial: 1s^1 (H-like Ne^9+);  core excited 1s --> 2p only;  captured into n=3 only;
    # final (He-like) states decay via radiative stabilization in n=1..2 shells.
    # Restricting to 1s->2p and n=3 gives a compact ~15-20 level table.
    setDefaults("nuclear: charge", 10.0)
    initialConfigs = [Configuration("1s")]
    fromShells     = [Shell("1s")]
    toShells       = [Shell("2p")]
    intoShells     = Basics.generateShellList(3, 3, 3)
    decayShells    = Basics.generateShellList(1, 2, 1)
    theme          = Basics.ForDielectronicRecombination(fromShells, toShells, intoShells, decayShells)
    computeResonanceStrength(theme, initialConfigs)
    #
elseif  true
    # Last successful:  27Jun2026
    # DR resonance strengths for Li-like Sc^18+ with 2s --> 2p core excitation.
    # Initial: 1s^2 2s (Li-like Sc^18+);  core excited 2s --> 2p;  captured into n=12 only.
    # n >= 12 needed: 2s-2p gap (~38 eV) > Rydberg binding at n=12 (~30 eV).
    # At n < 12 the DFS binding exceeds the gap → negative resonance energies → no pathways.
    # rbox = 25 a.u. needed: n=12 Rydberg extends to r_max ~ 16 a.u.
    # Radiative stabilization via 2p --> 2s decay; Rydberg electron is a spectator.
    setDefaults("nuclear: charge", 21.0)
    initialConfigs = [Configuration("1s^2 2s")]
    fromShells     = [Shell("2s")]
    toShells       = [Shell("2p")]
    intoShells     = Basics.generateShellList(12, 12, 3)
    decayShells    = [Shell("2s")]
    grid           = Radial.Grid(Radial.Grid(false), rnt = 1.0e-5, h = 5.0e-2, hp = 2.0e-2, rbox = 25.0)
    theme          = Basics.ForDielectronicRecombination(fromShells, toShells, intoShells, decayShells)
    computeResonanceStrength(theme, initialConfigs, grid=grid)
    #
end


setDefaults("print summary: close", "")
