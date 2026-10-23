
"""
`module  JAC.ManyElectronAZ`  
... a submodel of JAC that contains methods that support tasks related to spectroscopic computation.
"""
module ManyElectronAZ


using Printf, ..Basics, ..Defaults, ..InteractionStrength, ..ManyElectron, ..Radial, ..RadialIntegrals, ..SpinAngular

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
    opa = SpinAngular.OneParticleOperator(Mp.L, Basics.multipoleParity(Mp))
    wa  = SpinAngular.computeCoefficients(opa, leftBasis.csfs[r], rightBasis.csfs[s], subshellList) 
    for  coeff in wa
        MabEm = InteractionStrength.MabEmission(Mp, gauge, omega, rightBasis.orbitals[coeff.b],  
                                                                                leftBasis.orbitals[coeff.a], grid)
        ja = Basics.subshell_2j(leftBasis.orbitals[coeff.a].subshell)
        me = me + coeff.T * MabEm / sqrt( ja + 1) * sqrt( (Basics.twice(leftBasis.csfs[r].J) + 1)) 
    end
    
    return( me )
end


"""
`ManyElectron.matrixElement_HminusE(op::Basics.HminusE, leftBasis::Basis, r::Int64, rightBasis::Basis, s::Int64, grid::Radial.Grid,
                                    nuclearPotential::Radial.Potential)`
    ... computes the transition matrix element  < leftBasis.csfs[r] | H - E | rightBasis.csfs[s] >  for two CSF r, s as defined in the two
        given bases.  This is what a transition amplitude actually requires;  ManyElectron.matrixElement_Vee beside it computes only the
        two-body part, which is the exact answer wherever the two CSFs are built on one orthonormal orbital set and differ in TWO or more
        orbitals -- the first-order Auger case, and the reason that shortcut carries every validated Auger, DR and cascade number in JAC.
        A value::Float64 is returned.

        WHERE THE TWO DIFFER, AND BY HOW MUCH.  For a ONE-orbital difference a -> b the one-body part contributes <b|h|a>, and for orbitals
        that are eigenfunctions of a mean field f = h + U that equals -<b|U|a>, which cancels the two-body spectator direct+exchange terms
        <b|U_HF|a> down to the residual <b|(U_HF - U)|a>.  Keeping only V_ee therefore returns the WHOLE mean field, counted a second time
        after it had already generated both orbitals.  That case is not exotic: it is what every SECOND-ORDER sum produces at its
        Green-multiplet junctions.  Measured 23-Oct-2026 on Ar 2p^-2 -> 3p^-3 against Zitnik et al., PRA 93, 021401(R) (2016) Table I:
        V_ee alone gives 43.5x the measured rate, and dropping the spectator terms by hand gives 12x too little -- a factor 523 apart,
        bracketing the measurement from both sides.

        THE ONE-BODY TERM IS THE BARE h_D, through RadialIntegrals.GrantIab with the NUCLEAR potential, and not -U.  That is deliberate:
        the two are equal whenever the orbitals are eigenfunctions of their mean field, but the bare form needs no knowledge of WHICH mean
        field produced them and stays correct when they are eigenfunctions of nothing -- a frozen set, an imported set, or two different
        SCF solutions.  `nuclearPotential.Zr` must give the effective nuclear charge Z(r) on this grid.

        `op.kind` MUST BE THE INTERACTION OF THE HAMILTONIAN THAT PRODUCED THE LEVELS, because the cancellation is referenced to their
        energy.  Breit on top of Coulomb-only levels is an ADDITIONAL perturbation and carries no one-body counterterm, since it sits in
        neither U nor E;  pairing it here would subtract against the wrong Hamiltonian in silence.

        THE OVERLAP TERM -E <r|s> IS NOT INCLUDED, AND THAT IS A STATED ASSUMPTION RATHER THAN AN OVERSIGHT.  It vanishes identically under
        the orthonormality that ManyElectron.matrixElement_Vee beside it already assumes, and which the whole second-order path assumes
        today.  Lifting it -- genuine non-orthogonality between two SCF orbital sets -- is the biorthogonal-transformation route and is NOT
        addressed here;  until it is, this function is H - E evaluated under that assumption, which is exactly what the one-orbital case
        above needs and no more.
"""
function ManyElectron.matrixElement_HminusE(op::Basics.HminusE, leftBasis::Basis, r::Int64, rightBasis::Basis, s::Int64,
                                             grid::Radial.Grid, nuclearPotential::Radial.Potential)
    me   = 0.
    if  op.kind == DiagonalCoulomb()    return( me )     end
    symr = LevelSymmetry(leftBasis.csfs[r].J, leftBasis.csfs[r].parity)
    syms = LevelSymmetry(rightBasis.csfs[s].J, rightBasis.csfs[s].parity)
    if  symr != syms   @warn("Zero amplitude because of symmetries: $symr != $syms ");   return( me )    end

    subshellList = leftBasis.subshells
    # The ONE-BODY part, which is what distinguishes this function from matrixElement_Vee
    opa = SpinAngular.OneParticleOperator(0, plus)
    wa  = SpinAngular.computeCoefficients(opa, leftBasis.csfs[r], rightBasis.csfs[s], subshellList)
    for  coeff in wa
        me = me + coeff.T * RadialIntegrals.GrantIab(leftBasis.orbitals[coeff.a], rightBasis.orbitals[coeff.b],
                                                      grid, nuclearPotential)
    end
    # The TWO-BODY part, term for term as in ManyElectron.matrixElement_Vee
    opb = SpinAngular.TwoParticleOperator(0, plus)
    wb  = SpinAngular.computeCoefficients(opb, leftBasis.csfs[r], rightBasis.csfs[s], subshellList)
    for  coeff in wb
        if   typeof(op.kind) in [ DiagonalCoulomb, CoulombInteraction, CoulombBreit, CoulombGaunt]
            me = me + coeff.V * InteractionStrength.XL_Coulomb(coeff.nu,
                                        leftBasis.orbitals[coeff.a],  leftBasis.orbitals[coeff.b],
                                    rightBasis.orbitals[coeff.c], rightBasis.orbitals[coeff.d], grid)          end
        if   typeof(op.kind) in [ BreitInteraction, CoulombBreit, CoulombGaunt]
            me = me + coeff.V * InteractionStrength.XL_Breit(coeff.nu,
                                        leftBasis.orbitals[coeff.a],  leftBasis.orbitals[coeff.b],
                                    rightBasis.orbitals[coeff.c], rightBasis.orbitals[coeff.d], grid, op.kind)  end
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
    opa = SpinAngular.TwoParticleOperator(0, plus)
    wa  = SpinAngular.computeCoefficients(opa, leftBasis.csfs[r], rightBasis.csfs[s], subshellList) 
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
