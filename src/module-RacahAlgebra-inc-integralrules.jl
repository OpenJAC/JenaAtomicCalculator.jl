

"""
`RacahAlgebra.integralRulesForOneYlm(rex::RacahAlgebra.RacahExpression)`
    ... attempts to find a simplification of the given Racah expression by using the integration rule for one
        spherical harmonic over the full solid angle (Varshalovich et al. 1988, Eq. 5.9.1). Once a simplification
        is found, no attempt is made to find another simplifcation for this set of rules.
        A (istrue, rex)::Tuple{Bool, RacahExpression} is returned but where rex has no meaning for !istrue.

    AUDITED 21-Nov-2026: numerically verified by quadrature against AngularMomentum.sphericalYlm via
    RacahAlgebra.selectRacahIntegral(1); see project_racah_equivalentform_nondeterminism.md
"""
function integralRulesForOneYlm(rex::RacahAlgebra.RacahExpression)
    # Loop through all spherical harmonics
    for  (iaYlm, aYlm) in enumerate(rex.ylms)
        #
        #   Rule:      2pi          pi
        #   ----      Int    dphi  Int    dtheta  sin(theta)  Y_lm (theta,phi)   =   sqrt(4*pi)  delta(l,0)  delta(m,0)
        #               0            0
        #
        thetaIntegral = Integral(aYlm.theta, Basic(0), SymEngine.PI)
        phiIntegral   = Integral(aYlm.phi,   Basic(0), 2*SymEngine.PI)
        if  thetaIntegral in rex.integrals   &&   phiIntegral in rex.integrals
            testWeight   = rex.weight / sin(aYlm.theta)
            newYlms      = Ylm[];   for (ibYlm, bYlm) in enumerate(rex.ylms)  if  iaYlm != ibYlm   push!(newYlms, bYlm)   end   end
            newIntegrals = RacahAlgebra.removeIntegral(phiIntegral, RacahAlgebra.removeIntegral(thetaIntegral, rex.integrals))
            #
            if  RacahAlgebra.hasNoVars([aYlm.theta, aYlm.phi], rex.phase)        &&    RacahAlgebra.hasNoVars([aYlm.theta, aYlm.phi], testWeight)      &&
                RacahAlgebra.hasNoVars([aYlm.theta, aYlm.phi], rex.deltas)       &&    RacahAlgebra.hasNoVars([aYlm.theta, aYlm.phi], rex.triangles)   &&
                RacahAlgebra.hasNoVars([aYlm.theta, aYlm.phi], rex.w3js)         &&    RacahAlgebra.hasNoVars([aYlm.theta, aYlm.phi], rex.w6js)        &&
                RacahAlgebra.hasNoVars([aYlm.theta, aYlm.phi], rex.w9js)         &&    RacahAlgebra.hasNoVars([aYlm.theta, aYlm.phi], newYlms)         &&
                RacahAlgebra.hasNoVars([aYlm.theta, aYlm.phi], rex.djpqs)        &&    RacahAlgebra.hasNoVars([aYlm.theta, aYlm.phi], newIntegrals)
                newDeltas = rex.deltas;   push!(newDeltas, Kronecker(aYlm.l, 0));   push!(newDeltas, Kronecker(aYlm.m, 0))
                wa = RacahExpression( rex.summations, newIntegrals, rex.phase, testWeight * sqrt(4*SymEngine.PI), newDeltas, rex.triangles,
                                        rex.w3js, rex.w6js, rex.w9js, newYlms, rex.djpqs )
                println(">> Apply integral rule for one Ylm -- Int dOmega Y_lm ...")
                return( (true, wa) )
            end
        end
    end

    return( (false, RacahExpression()) )
end


"""
`RacahAlgebra.integralRulesForThreeYlm(rex::RacahAlgebra.RacahExpression)`
    ... attempts to find a simplification of the given Racah expression by using the Gaunt-type integration rule
        for three spherical harmonics, all taken at the same (theta,phi), over the full solid angle (Varshalovich
        et al. 1988, Eq. 5.9.5). Once a simplification is found, no attempt is made to find another simplifcation
        for this set of rules.
        A (istrue, rex)::Tuple{Bool, RacahExpression} is returned but where rex has no meaning for !istrue.

    AUDITED 21-Nov-2026: numerically verified by quadrature against AngularMomentum.sphericalYlm via
    RacahAlgebra.selectRacahIntegral(2); see project_racah_equivalentform_nondeterminism.md
"""
function integralRulesForThreeYlm(rex::RacahAlgebra.RacahExpression)
    # Make sure that the rule could apply at all
    if  length(rex.ylms) < 3     return( (false, RacahExpression()) )    end
    #
    # Loop through all triples of spherical harmonics taken at the same (theta,phi)
    for  (iaYlm, aYlm) in enumerate(rex.ylms)
        for  (ibYlm, bYlm) in enumerate(rex.ylms)
            if  ibYlm == iaYlm    continue    end
            if  bYlm.theta != aYlm.theta  ||  bYlm.phi != aYlm.phi    continue    end
            for  (icYlm, cYlm) in enumerate(rex.ylms)
                if  icYlm == iaYlm  ||  icYlm == ibYlm    continue    end
                if  cYlm.theta != aYlm.theta  ||  cYlm.phi != aYlm.phi    continue    end
                #
                #   Rule:      2pi          pi
                #   ----      Int    dphi  Int    dtheta  sin(theta)  Y_{l1 m1}(theta,phi)  Y_{l2 m2}(theta,phi)  Y_{l3 m3}(theta,phi)
                #               0            0
                #
                #                 (2l1+1) (2l2+1) (2l3+1)  1/2   ( l1 l2 l3 )   ( l1 l2 l3 )
                #           =    (----------------------- )     (          ) {(          )}
                #                            4 pi                ( 0  0  0  )   ( m1 m2 m3 )
                #
                thetaIntegral = Integral(aYlm.theta, Basic(0), SymEngine.PI)
                phiIntegral   = Integral(aYlm.phi,   Basic(0), 2*SymEngine.PI)
                if  thetaIntegral in rex.integrals   &&   phiIntegral in rex.integrals
                    testWeight   = rex.weight / sin(aYlm.theta)
                    newYlms      = Ylm[];   for (idYlm, dYlm) in enumerate(rex.ylms)
                                                 if  idYlm != iaYlm && idYlm != ibYlm && idYlm != icYlm   push!(newYlms, dYlm)   end   end
                    newIntegrals = RacahAlgebra.removeIntegral(phiIntegral, RacahAlgebra.removeIntegral(thetaIntegral, rex.integrals))
                    #
                    if  RacahAlgebra.hasNoVars([aYlm.theta, aYlm.phi], rex.phase)        &&    RacahAlgebra.hasNoVars([aYlm.theta, aYlm.phi], testWeight)      &&
                        RacahAlgebra.hasNoVars([aYlm.theta, aYlm.phi], rex.deltas)       &&    RacahAlgebra.hasNoVars([aYlm.theta, aYlm.phi], rex.triangles)   &&
                        RacahAlgebra.hasNoVars([aYlm.theta, aYlm.phi], rex.w3js)         &&    RacahAlgebra.hasNoVars([aYlm.theta, aYlm.phi], rex.w6js)        &&
                        RacahAlgebra.hasNoVars([aYlm.theta, aYlm.phi], rex.w9js)         &&    RacahAlgebra.hasNoVars([aYlm.theta, aYlm.phi], newYlms)         &&
                        RacahAlgebra.hasNoVars([aYlm.theta, aYlm.phi], rex.djpqs)        &&    RacahAlgebra.hasNoVars([aYlm.theta, aYlm.phi], newIntegrals)
                        newWeight = testWeight * sqrt( (2*aYlm.l+1)*(2*bYlm.l+1)*(2*cYlm.l+1) / (4*SymEngine.PI) )
                        newW3js   = rex.w3js;   push!(newW3js, W3j(aYlm.l, bYlm.l, cYlm.l, Basic(0), Basic(0), Basic(0)))
                                                 push!(newW3js, W3j(aYlm.l, bYlm.l, cYlm.l, aYlm.m,   bYlm.m,   cYlm.m))
                        wa = RacahExpression( rex.summations, newIntegrals, rex.phase, newWeight, rex.deltas, rex.triangles,
                                                newW3js, rex.w6js, rex.w9js, newYlms, rex.djpqs )
                        println(">> Apply integral rule for three Ylm -- Int dOmega Y_l1m1 Y_l2m2 Y_l3m3 ...")
                        return( (true, wa) )
                    end
                end
            end
        end
    end

    return( (false, RacahExpression()) )
end
