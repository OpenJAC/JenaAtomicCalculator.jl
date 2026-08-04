
println("Me) Test  displayCouplings(FineStructure(), ...)  and  displayCouplings(FineStructureLS(), ...)  for argon with two inner-shell holes.")

setDefaults("print summary: open", "zzz-ForPedestrians.sum")


if  false
    # Last successful:  26Jun2026
    # Ar^2+ with two L-shell vacancies: the three KLL Auger final-state configuration groups.
    # Shows jj (FineStructure) and LS (FineStructureLS) coupling of 2p^4, 2s 2p^5, and 2p^6 open shells.
    configs = [Configuration("1s^2 2s^2 2p^4 3s^2 3p^6"),
               Configuration("1s^2 2s^1 2p^5 3s^2 3p^6"),
               Configuration("1s^2 2p^6 3s^2 3p^6")]
    displayCouplings(Basics.FineStructure(),   configs)
    displayCouplings(Basics.FineStructureLS(), configs)
    #
elseif  true
    # Last successful:  26Jun2026
    # Ar^2+ with one K-shell and one L-shell vacancy (KL double holes).
    # Shows jj and LS coupling of the 1s + 2p and 1s + 2s open-shell combinations.
    configs = [Configuration("1s^1 2s^2 2p^5 3s^2 3p^6"),
               Configuration("1s^1 2s^1 2p^6 3s^2 3p^6")]
    displayCouplings(Basics.FineStructure(),   configs)
    displayCouplings(Basics.FineStructureLS(), configs)
    #
end


setDefaults("print summary: close", "")
