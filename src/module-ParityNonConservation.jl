
"""
`module  JAC.ParityNonConservation`  
    ... a submodel of JAC that contains all methods for computing P-odd and/or T-odd amplitudes between two 
        bound-state levels.
"""
module ParityNonConservation


using  Printf,  ..Basics, ..Defaults, ..InteractionStrength, ..ManyElectron, ..Nuclear, ..Radial, ..SpinAngular

"""
`ParityNonConservation.schiffMomentAmplitude(finalLevel::Level, initialLevel::Level, nm::Nuclear.Model, 
                                                    grid::Radial.Grid; display::Bool=false)`  
    ... to compute the Schiff moment amplitude  <alpha_f J_f || H^(Schiff) || alpha_i J_i>  for the given final 
        and initial level and for the nuclear density as given by the nuclear model. A value::ComplexF64 is returned.
"""
function schiffMomentAmplitude(finalLevel::Level, initialLevel::Level, nm::Nuclear.Model, grid::Radial.Grid; display::Bool=false)
    nf = length(finalLevel.basis.csfs);    ni = length(initialLevel.basis.csfs)
    printstyled("Compute Schiff moment  matrix of dimension $nf x $ni in the final- and initial-state bases " *
                "for the transition [$(initialLevel.index)- $(finalLevel.index)] ... ", color=:light_green)
    matrix = zeros(ComplexF64, nf, ni)
    #
    for  r = 1:nf
        for  s = 1:ni
            # Calculate the spin-angular coefficients
            if  Defaults.saRatip()
                waR = Basics.compute(AngularCoeffs1pGrasp92(), 0, 1, finalLevel.basis.csfs[r], initialLevel.basis.csfs[s])
                wa  = waR       
            end
            if  Defaults.saGG()
                subshellList = initialLevel.basis.subshells
                opa = SpinAngular.OneParticleOperator(1, plus, true)
                waG = SpinAngular.computeCoefficients(opa, finalLevel.basis.csfs[r], initialLevel.basis.csfs[s], subshellList) 
                wa  = waG
            end
            if  Defaults.saRatip() && Defaults.saGG() && true
                if  length(waR) != 0     println("\n>> Angular coeffients from GRASP/MCT   = $waR ")    end
                if  length(waG) != 0     println(  ">> Angular coeffients from SpinAngular = $waG ")    end
            end
            #
            for  coeff in wa
                tamp  = InteractionStrength.schiffMoment(finalLevel.basis.orbitals[coeff.a], initialLevel.basis.orbitals[coeff.b], nm, grid)
                matrix[r,s] = matrix[r,s] + coeff.T * tamp  
            end
        end
    end
    printstyled("done. \n", color=:light_green)
    amplitude = transpose(finalLevel.mc) * matrix * initialLevel.mc 
    #
    if  display
        sa = @sprintf("%.8e", amplitude.re) * @sprintf("%.8e", amplitude.im)
        println("   Schiff moment amplitude:   "                                                                               *
                "< level=$(finalLevel.index) [J=$(finalLevel.J)$(string(finalLevel.parity))] || H^(Schiff) ($(nm.model)) ||"   *
                " $(initialLevel.index) [$(initialLevel.J)$(string(initialLevel.parity))] >  = " * sa)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        if  printSummary
            println(iostream, "   Schiff moment amplitude:   "                                                                               *
                                "< level=$(finalLevel.index) [J=$(finalLevel.J)$(string(finalLevel.parity))] || H^(Schiff) ($(nm.model)) ||"   *
                                " $(initialLevel.index) [$(initialLevel.J)$(string(initialLevel.parity))] >  = " * sa)
        end
    end
    
    return( amplitude )
end


"""
`ParityNonConservation.anapoleMomentAmplitude(finalLevel::Level, initialLevel::Level, grid::Radial.Grid; display::Bool=false)`  
    ... to compute the anapole moment amplitude <(alpha_f J_f, kappa) J_i || H^(anapole) || alpha_i J_i>  
        for the given final and initial level. A value::ComplexF64 is returned.
"""
function anapoleMomentAmplitude(finalLevel::Level, initialLevel::Level, grid::Radial.Grid; display::Bool=false)
    nf = length(finalLevel.basis.csfs);    ni = length(initialLevel.basis.csfs)
    printstyled("Compute anapole moment matrix of dimension $nf x $ni in the final- and initial-state bases " *
                                    "for the transition [$(initialLevel.index)- $(finalLevel.index)] ... ", color=:light_green)
    matrix = zeros(ComplexF64, nf, ni)
    #
    #
    printstyled("done. \n", color=:light_green)
    #
    amplitude = 3.0 + 3.0*im
    #
    if  display
            sa = @sprintf("%.8e", amplitude.re) * @sprintf("%.8e", amplitude.im)
        println("   Anapole moment amplitude:  "                                                                  *
                "< level=$(finalLevel.index) [J=$(finalLevel.J)$(string(finalLevel.parity))] || H^(anapole) ||"   *
                " $(initialLevel.index) [$(initialLevel.J)$(string(initialLevel.parity))] >  = " * sa)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        if  printSummary
            println(iostream, "   Anapole moment amplitude:  "                                                                  *
                                "< level=$(finalLevel.index) [J=$(finalLevel.J)$(string(finalLevel.parity))] || H^(anapole) ||"   *
                                " $(initialLevel.index) [$(initialLevel.J)$(string(initialLevel.parity))] >  = " * sa)
        end
    end
    
    return( amplitude )
end


"""
`ParityNonConservation.weakChargeAmplitude(finalLevel::Level, initialLevel::Level, nm::Nuclear.Model, 
                                                grid::Radial.Grid; display::Bool=false)`  
    ... to compute the weak-charge amplitude  <alpha_f J_f || H^(weak-charge) || alpha_i J_i>  for the given final and 
        initial level and for the nuclear density as given by the nuclear model. A value::ComplexF64 is returned.
"""
function weakChargeAmplitude(finalLevel::Level, initialLevel::Level, nm::Nuclear.Model, grid::Radial.Grid; display::Bool=false)
    nf = length(finalLevel.basis.csfs);    ni = length(initialLevel.basis.csfs)
    printstyled("Compute weak -charge matrix of dimension $nf x $ni in the final- and initial-state bases " *
                "for the transition [$(initialLevel.index)- $(finalLevel.index)] ... ", color=:light_green)
    matrix = zeros(ComplexF64, nf, ni)
    #
    for  r = 1:nf
        for  s = 1:ni
            # Calculate the spin-angular coefficients
            if  Defaults.saRatip()
                waR = Basics.compute(AngularCoeffs1pGrasp92(), 0, 1, finalLevel.basis.csfs[r], initialLevel.basis.csfs[s])
                wa  = waR       
            end
            if  Defaults.saGG()
                subshellList = initialLevel.basis.subshells
                opa = SpinAngular.OneParticleOperator(1, plus, true)
                waG = SpinAngular.computeCoefficients(opa, finalLevel.basis.csfs[r], initialLevel.basis.csfs[s], subshellList) 
                wa  = waG
            end
            if  Defaults.saRatip() && Defaults.saGG() && true
                if  length(waR) != 0     println("\n>> Angular coeffients from GRASP/MCT   = $waR ")    end
                if  length(waG) != 0     println(  ">> Angular coeffients from SpinAngular = $waG ")    end
            end
            #
            for  coeff in wa
                tamp  = InteractionStrength.weakCharge(finalLevel.basis.orbitals[coeff.a], initialLevel.basis.orbitals[coeff.b], nm, grid)
                matrix[r,s] = matrix[r,s] + coeff.T * tamp  
            end
        end
    end
    printstyled("done. \n", color=:light_green)
    amplitude = transpose(finalLevel.mc) * matrix * initialLevel.mc 
    #
    if  display
        sa = @sprintf("%.5e", amplitude.re) * "  " * @sprintf("%.5e", amplitude.im)
        println("   weak-charge amplitude:     "                                                         *
                "< level=$(finalLevel.index) [J=$(finalLevel.J)$(string(finalLevel.parity))] || H^(weak-charge) ($(nm.model)) ||"   *
                " $(initialLevel.index) [$(initialLevel.J)$(string(initialLevel.parity))] >  = " * sa)
        printSummary, iostream = Defaults.getDefaults("summary flag/stream")
        if  printSummary
            println(iostream, "   weak-charge amplitude:     "                                                         *
                                "< level=$(finalLevel.index) [J=$(finalLevel.J)$(string(finalLevel.parity))] || H^(weak-charge) ($(nm.model)) ||"   *
                                " $(initialLevel.index) [$(initialLevel.J)$(string(initialLevel.parity))] >  = " * sa)
        end
    end
    
    return( amplitude )
end
    
end # module
