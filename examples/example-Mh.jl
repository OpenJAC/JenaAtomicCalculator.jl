
println("Mh) Test  computeChargeStateDistribution(ForStepwiseDecay(), ...)  for inner-shell hole configurations.")

setDefaults("unit: energy", "eV")
setDefaults("print summary: open", "zzz-ForPedestrians.sum")


if  false
    # Last successful:  28Jun2026
    # K-hole decay cascade in Ne: charge state distribution after KLL Auger and radiative decay.
    # Initial: 1s vacancy in neutral Ne (9 electrons); 3 electron-loss steps cover KLL and
    # subsequent single-step Auger/radiative channels.  Dominant channel: KLL Auger -> Ne^2+.
    setDefaults("nuclear: charge", 10.0)
    initialConfigs = [Configuration("1s 2s^2 2p^6")]
    computeChargeStateDistribution(Basics.ForStepwiseDecay(3), initialConfigs)
    #
elseif  true
    # Last successful:  28Jun2026
    # L-hole (2p hole) decay cascade in Ar: charge state distribution after LMM Auger and radiative decay.
    # Initial: 2p vacancy in neutral Ar (17 electrons); 2 electron-loss steps cover LMM and
    # subsequent channels.  Expected: Ar^2+ dominant (99.99%), L fluorescence yield ~0.01%.
    setDefaults("nuclear: charge", 18.0)
    initialConfigs = [Configuration("1s^2 2s^2 2p^5 3s^2 3p^6")]
    computeChargeStateDistribution(Basics.ForStepwiseDecay(2), initialConfigs)
    #
end


setDefaults("print summary: close", "")
