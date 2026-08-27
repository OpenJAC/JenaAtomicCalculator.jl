
"""
`module  JAC.ManyElectronAZ`  
... a submodel of JAC that contains methods that support tasks related to spectroscopic computation.
"""
module ManyElectronAZ


using Printf, ..Basics, ..Defaults, ..InteractionStrength, ..ManyElectron, ..Radial, ..SpinAngularNew

export provideSubshellStates


"""
`ManyElectron.matrixElement_Mab(Mp::EmMultipole, gauge::EmGauge, omega::Float64, 
                                    leftBasis::Basis, r::Int64, rightBasis::Basis, s::Int64, grid::Radial.Grid)`  
    ... to compute the electron-photon interaction amplitude  
    
                <leftBasis.csfs[r] || O^(Mp, emission) || rightBasis.csfs[s]>
                
        for two CSF r, s as define in the two given bases. An amplitude::ComplexF64 is returned.
        The procedure assumes that the two bases are defined on the same sequence of subshells and that their orbitals
        are 'orthogonal' in practice.
"""
function ManyElectron.matrixElement_Mab(Mp::EmMultipole, gauge::EmGauge, omega::Float64, 
                                        leftBasis::Basis, r::Int64, rightBasis::Basis, s::Int64, grid::Radial.Grid)
    # Issue a warning if the amplitude is zero due to the given symmetries of the CSFs; this should be tested before
    me  = 0.
        
    subshellList = leftBasis.subshells
    opa = SpinAngularNew.OneParticleOperator(Mp.L, Basics.multipoleParity(Mp))
    wa  = SpinAngularNew.computeCoefficients(opa, leftBasis.csfs[r], rightBasis.csfs[s], subshellList) 
    for  coeff in wa
        MabEm = InteractionStrength.MabEmission(Mp, gauge, omega, leftBasis.orbitals[coeff.a],  
                                                                                rightBasis.orbitals[coeff.b], grid)
        ja = Basics.subshell_2j(leftBasis.orbitals[coeff.a].subshell)
        me = me + coeff.T * MabEm / sqrt( ja + 1) * sqrt( (Basics.twice(leftBasis.csfs[r].J) + 1)) 
    end
    
    return( me )
end


"""
`ManyElectron.matrixElement_Vee(kind::AbstractEeInteraction, leftBasis::Basis, r::Int64, rightBasis::Basis, s::Int64, grid::Radial.Grid)`  
    ... to compute the V^(e-e) electron-electron interaction amplitude  
    
                < leftBasis.csfs[r] || V^(e-e) || rightBasis.csfs[s] > 
                
        for two CSF r, s as define in the two given bases. An amplitude::ComplexF64 is returned.
        The procedure assumes that the two bases are defined on the same sequence of subshells and that their orbitals
        are 'orthogonal' in practice.
"""
function ManyElectron.matrixElement_Vee(kind::AbstractEeInteraction, leftBasis::Basis, r::Int64, rightBasis::Basis, s::Int64, grid::Radial.Grid)
    # Issue a warning if the amplitude is zero due to the given symmetries of the CSFs; this should be tested before
    me   = 0.
    if  kind == DiagonalCoulomb()    return(me)     end
    symr = LevelSymmetry(leftBasis.csfs[r].J, leftBasis.csfs[r].parity)
    syms = LevelSymmetry(rightBasis.csfs[s].J, rightBasis.csfs[s].parity)
    if  symr != syms   @warn("Zero amplitude because of symmetries: $symr != $syms ");   return( me )    end
        
    subshellList = leftBasis.subshells
    opa = SpinAngularNew.TwoParticleOperator(0, plus)
    wa  = SpinAngularNew.computeCoefficients(opa, leftBasis.csfs[r], rightBasis.csfs[s], subshellList) 
    #
    for  coeff in wa
        if   typeof(kind) in [ DiagonalCoulomb, CoulombInteraction, CoulombBreit, CoulombGaunt]    
            me = me + coeff.V * InteractionStrength.XL_Coulomb(coeff.nu, 
                                        leftBasis.orbitals[coeff.a],  leftBasis.orbitals[coeff.b],
                                    rightBasis.orbitals[coeff.c], rightBasis.orbitals[coeff.d], grid)         end
        if   typeof(kind) in [ BreitInteraction, CoulombBreit, CoulombGaunt]    
            me = me + coeff.V * InteractionStrength.XL_Breit(coeff.nu, 
                                        leftBasis.orbitals[coeff.a],  leftBasis.orbitals[coeff.b],
                                    rightBasis.orbitals[coeff.c], rightBasis.orbitals[coeff.d], grid, kind)   end
    end
    
    return( me )
end    


"""
`ManyElectron.provideSubshellStates(sh::Subshell, occ::Int64)`   
    ... to provide all antisymmetric SubshellStates within the relativistic seniority scheme for the given 
        shell and occupation; an Array{SubshellStateR,1} is returned.
"""
function ManyElectron.provideSubshellStates(sh::Subshell, occ::Int64)

    # Use the j-value = wa[4] of the subshell to select the allowed SubshellStateR's
    wb = Basics.SubshellStateR[]

    if       Basics.subshell_2j(sh) == 1
        if       occ == 0  || occ == 2      push!( wb, Basics.SubshellStateR( sh, occ, 0, 0) );       return(wb)
        elseif   occ == 1                   push!( wb, Basics.SubshellStateR( sh, occ, 1, 1) );       return(wb)
        end
    elseif   Basics.subshell_2j(sh) == 3
        if       occ == 0  || occ == 4      push!( wb, Basics.SubshellStateR( sh, occ, 0, 0) );       return(wb)
        elseif   occ == 1  || occ == 3      push!( wb, Basics.SubshellStateR( sh, occ, 1, 3) );       return(wb)
        elseif   occ == 2                   push!( wb, Basics.SubshellStateR( sh, occ, 0, 0) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 2, 4) );       return(wb)
        end
    elseif   Basics.subshell_2j(sh) == 5
        if       occ == 0  || occ == 6      push!( wb, Basics.SubshellStateR( sh, occ, 0, 0) );       return(wb)
        elseif   occ == 1  || occ == 5      push!( wb, Basics.SubshellStateR( sh, occ, 1, 5) );       return(wb)
        elseif   occ == 2  || occ == 4      push!( wb, Basics.SubshellStateR( sh, occ, 0, 0) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 2, 4) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 2, 8) );       return(wb)
        elseif   occ == 3                   push!( wb, Basics.SubshellStateR( sh, occ, 1, 5) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 3, 3) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 3, 9) );       return(wb)
        end
    elseif   Basics.subshell_2j(sh) == 7
        if       occ == 0  || occ == 8      push!( wb, Basics.SubshellStateR( sh, occ, 0, 0) );       return(wb)
        elseif   occ == 1  || occ == 7      push!( wb, Basics.SubshellStateR( sh, occ, 1, 7) );       return(wb)
        elseif   occ == 2  || occ == 6      push!( wb, Basics.SubshellStateR( sh, occ, 0, 0) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 2, 4) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 2, 8) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 2,12) );       return(wb)
        elseif   occ == 3  || occ == 5      push!( wb, Basics.SubshellStateR( sh, occ, 1, 7) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 3, 3) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 3, 5) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 3, 9) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 3,11) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 3,15) );       return(wb)
        elseif   occ == 4                   push!( wb, Basics.SubshellStateR( sh, occ, 0, 0) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 2, 4) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 2, 8) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 2,12) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 4, 4) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 4, 8) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 4,10) )
                                            push!( wb, Basics.SubshellStateR( sh, occ, 4,16) );       return(wb)
        end
    elseif   Basics.subshell_2j(sh) in 9:2:2001
        if       occ == 0                   push!( wb, Basics.SubshellStateR( sh, occ, 0, 0) );                 return(wb)
        elseif   occ == 1                   push!( wb, Basics.SubshellStateR( sh, occ, 1, subshell_2j(sh)) );   return(wb)
        elseif   occ == 2
                                            push!( wb, Basics.SubshellStateR( sh, occ, 0, 0) )
            for N = 2:Basics.subshell_2j(sh)-1
                if iseven(N) == true        push!( wb, Basics.SubshellStateR( sh, occ, 2, 2*N ) )  end
            end
            return(wb)
        else     error("stop b")
        end

    # For an extension to other subshell-occupations, see the Racah-manual; Table ...
    else  error("stop a")
    end
end

end # module
