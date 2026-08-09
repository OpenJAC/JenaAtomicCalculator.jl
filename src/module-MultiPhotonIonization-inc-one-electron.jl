


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
        amp = amp + InteractionStrength.MbaAbsorptionCheng(mp2, gaugex, omega2, fOrbital, vOrbital, grid) *
                    InteractionStrength.MbaAbsorptionCheng(mp1, gaugex, omega1, vOrbital, iOrbital, grid) /
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
