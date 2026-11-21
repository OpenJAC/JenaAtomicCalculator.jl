

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
