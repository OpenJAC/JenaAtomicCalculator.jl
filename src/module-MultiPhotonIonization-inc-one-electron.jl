


"""
`MultiPhotonIonization.oneElectronAmplitude(fOrbital::Orbital, omega2::Float64, mp2::EmMultipole, symx::LevelSymmetry, 
                                                                omega1::Float64, mp1::EmMultipole, iOrbital::Orbital, 
                                                                gauge::UseGauge, orbitals::Dict{Subshell, Orbital}, grid::Radial.Grid)` 
    ... to compute the two-photon ionization amplitudes for the given initial and final orbitals and quantum numbers (symmetries)
        of the intermediate partial waves. An  amp:ComplexF64  is returned.
"""
function  oneElectronAmplitude(fOrbital::Orbital, omega2::Float64, mp2::EmMultipole, symx::LevelSymmetry, 
                                                    omega1::Float64, mp1::EmMultipole, iOrbital::Orbital, 
                                                    gauge::UseGauge, orbitals::Dict{Subshell, Orbital}, grid::Radial.Grid)
    amp = 0.0im
    
    for (k, vOrbital) in orbitals
        if  LevelSymmetry(k) != symx    continue    end
        if      gauge == UseCoulomb     gaugex = Basics.Coulomb     
        elseif  gauge == UseBabushkin   gaugex = Basics.Babushkin
        end
        ## MIGRATED 09-Aug-2026 from MbaAbsorptionCheng to the consolidated InteractionStrength.MabEmission.
        ## TWO THINGS CHANGE HERE, both deliberate and both altering this module's absolute numbers.
        ##  (1) ARGUMENT ORDER. MbaAbsorptionCheng(...,b,a) returned <b||O||a> whereas MabEmission(...,b,a)
        ##      returns <a||O||b>, so the two orbital arguments are exchanged below to keep the same quantity.
        ##  (2) NORMALISATION. The Cheng and Johnson families differ by exactly sqrt(2*pi/(L(2L+1))) -- measured
        ##      1.447203, 0.792665 and 0.546991 for L = 1, 2, 3 against that formula, the L = 3 value being a
        ##      prediction it was not fitted to. Only the Johnson normalisation is validated against data
        ##      (Jitrik & Bunge 2004: E1, M1 and E2 rates in hydrogen to 0.05 %, 0.5 % and 0.05 %), so this
        ##      module now uses it. Its one-electron amplitudes therefore change by sqrt(L(2L+1)/(2pi)), i.e.
        ##      by 0.691 for E1, and are for the first time consistent with the rest of JAC.
        ## Absorption versus emission is a complex conjugation, which is the identity now that the matrix
        ## element is real; see the phase note in AngularMomentum.JohnsonI.
        amp = amp + InteractionStrength.MabEmission(mp2, gaugex, omega2, vOrbital, fOrbital, grid) *
                    InteractionStrength.MabEmission(mp1, gaugex, omega1, iOrbital, vOrbital, grid) /
                    (iOrbital.energy + omega1 - vOrbital.energy)
    end
    
    return( amp )
end


"""
`MultiPhotonIonization.oneElectronComputeTwoPhotonLine(iState::Subshell, multipoles::Array{EmMultipole}, gauges::Array{UseGauge}, 
                                                        omega1::Float64, omega2::Float64, pqnMax::Int64, 
                                                        orbitals::Dict{Subshell, Orbital}, meanPot::Radial.Potential; output=true)` 
    ... to compute the multiphoton transition amplitudes and all properties as requested by the given settings. 
        A list of lines::Array{MultiPhotonIonization.Lines} is returned.
"""
function  oneElectronComputeTwoPhotonLine(iState::Subshell, multipoles::Array{EmMultipole}, gauges::Array{UseGauge},
                                            omega1::Float64, omega2::Float64, pqnMax::Int64, 
                                            orbitals::Dict{Subshell, Orbital}, meanPot::Radial.Potential; output=false)
    println("")
    printstyled("MultiPhotonIonization.oneElectronComputeTwoPhotonLine(): The computation starts now ... \n", color=:light_green)
    printstyled("---------------------------------------------------- ---------------------------------- \n", color=:light_green)
    println("")
    #
    epsilon = omega1 + omega2 + orbitals[iState].energy
    println("Energy i-orbital: $(orbitals[iState].energy) a.u.")
    println("omega1:            $(omega1) a.u.")
    println("omega2:            $(omega2) a.u.")
    println("Excess energy:      $epsilon  a.u.  \n")
    
    symi = LevelSymmetry(iState)
    #
    # Determine and compute all two-photon amplitudes for the given input
    for  mp1 in multipoles
        for  mp2 in multipoles
            symxList = AngularMomentum.allowedMultipoleSymmetries(symi, mp1)
            for  symx in symxList
                symfList = AngularMomentum.allowedMultipoleSymmetries(symx, mp2)
                for symf in symfList
                    for gauge in gauges
                        iOrbital = orbitals[iState]
                        # Generate a proper continuum orbital with symmetry symf and energy  epsilon
                        nrContinuum = Continuum.gridConsistency(epsilon, meanPot.grid)
                        settings    = Continuum.Settings(false, nrContinuum)
                        fOrbitalPhs = Continuum.generateOrbitalLocalPotential(epsilon, Subshell(1001, symf), meanPot, settings)
                        fOrbital    = fOrbitalPhs[1]
                        amp = MultiPhotonIonization.oneElectronAmplitude(fOrbital, omega2, mp2, symx, omega1, mp1, 
                                                                            iOrbital, gauge, orbitals, meanPot.grid)
                        println("$mp1  $mp2  $symi  $symx  $symf  $gauge  $amp")
                    end
                end
            end
        end
    end
    # 
    if    output    return( true )
    else            return( nothing )
    end
end
