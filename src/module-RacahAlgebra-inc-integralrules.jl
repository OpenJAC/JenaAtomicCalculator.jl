

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
`RacahAlgebra.integralRulesForTwoYlm(rex::RacahAlgebra.RacahExpression)`
    ... attempts to find a simplification of the given Racah expression by using the orthonormality integration
        rule for two spherical harmonics, both taken at the same (theta,phi), over the full solid angle
        (Varshalovich et al. 1988, Eq. 5.9.2). Once a simplification is found, no attempt is made to find
        another simplifcation for this set of rules.
        A (istrue, rex)::Tuple{Bool, RacahExpression} is returned but where rex has no meaning for !istrue.

    NOTE ON THE SIGN: `RacahAlgebra.Ylm(l, m, theta, phi, true)` already rewrites a complex-conjugated Y*_lm as
    the PLAIN Y_{l,-m} with an extra (-1)^(-m) folded into the surrounding RacahExpression's phase -- so by the
    time a Ylm reaches rex.ylms it is never itself conjugated, and the identity this rule implements must be
    stated for two PLAIN (unconjugated) harmonics: Int dOmega Y_l1m1 Y_l2m2 = (-1)^m1 delta(l1,l2) delta(m1,-m2),
    which is the textbook Int dOmega Y*_lm Y_l'm' = delta_ll' delta_mm' with m2 -> -m2 and the matching (-1)^m
    from Y_{l2,-m2} = (-1)^m2 Y*_{l2,m2} made explicit instead of hidden in a star flag.

    AUDITED 21-Nov-2026: numerically verified by quadrature against AngularMomentum.sphericalYlm via
    RacahAlgebra.selectRacahIntegral(2); see project_racah_equivalentform_nondeterminism.md
"""
function integralRulesForTwoYlm(rex::RacahAlgebra.RacahExpression)
    # Make sure that the rule could apply at all
    if  length(rex.ylms) < 2     return( (false, RacahExpression()) )    end
    #
    # Loop through all pairs of spherical harmonics taken at the same (theta,phi)
    for  (iaYlm, aYlm) in enumerate(rex.ylms)
        for  (ibYlm, bYlm) in enumerate(rex.ylms)
            if  ibYlm == iaYlm    continue    end
            if  bYlm.theta != aYlm.theta  ||  bYlm.phi != aYlm.phi    continue    end
            #
            #   Rule:      2pi          pi
            #   ----      Int    dphi  Int    dtheta  sin(theta)  Y_{l1 m1}(theta,phi)  Y_{l2 m2}(theta,phi)
            #               0            0
            #
            #                  m1
            #           =  (-1)    delta(l1,l2)  delta(m1,-m2)
            #
            thetaIntegral = Integral(aYlm.theta, Basic(0), SymEngine.PI)
            phiIntegral   = Integral(aYlm.phi,   Basic(0), 2*SymEngine.PI)
            if  thetaIntegral in rex.integrals   &&   phiIntegral in rex.integrals
                testWeight = rex.weight / sin(aYlm.theta)
                newYlms    = Ylm[];   for (icYlm, cYlm) in enumerate(rex.ylms)  if  icYlm != iaYlm && icYlm != ibYlm   push!(newYlms, cYlm)   end   end
                newIntegrals = RacahAlgebra.removeIntegral(phiIntegral, RacahAlgebra.removeIntegral(thetaIntegral, rex.integrals))
                #
                if  RacahAlgebra.hasNoVars([aYlm.theta, aYlm.phi], rex.phase)        &&    RacahAlgebra.hasNoVars([aYlm.theta, aYlm.phi], testWeight)      &&
                    RacahAlgebra.hasNoVars([aYlm.theta, aYlm.phi], rex.deltas)       &&    RacahAlgebra.hasNoVars([aYlm.theta, aYlm.phi], rex.triangles)   &&
                    RacahAlgebra.hasNoVars([aYlm.theta, aYlm.phi], rex.w3js)         &&    RacahAlgebra.hasNoVars([aYlm.theta, aYlm.phi], rex.w6js)        &&
                    RacahAlgebra.hasNoVars([aYlm.theta, aYlm.phi], rex.w9js)         &&    RacahAlgebra.hasNoVars([aYlm.theta, aYlm.phi], newYlms)         &&
                    RacahAlgebra.hasNoVars([aYlm.theta, aYlm.phi], rex.djpqs)        &&    RacahAlgebra.hasNoVars([aYlm.theta, aYlm.phi], newIntegrals)
                    newDeltas = rex.deltas;   push!(newDeltas, Kronecker(aYlm.l, bYlm.l));   push!(newDeltas, Kronecker(aYlm.m, -bYlm.m))
                    wa = RacahExpression( rex.summations, newIntegrals, RacahAlgebra.purifyPhase(rex.phase + aYlm.m), testWeight, newDeltas, rex.triangles,
                                            rex.w3js, rex.w6js, rex.w9js, newYlms, rex.djpqs )
                    println(">> Apply integral rule for two Ylm -- Int dOmega Y_l1m1 Y_l2m2 ...")
                    return( (true, wa) )
                end
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
    RacahAlgebra.selectRacahIntegral(3); see project_racah_equivalentform_nondeterminism.md
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


##### Integration rules for the Wigner D-functions / small-d rotation matrices #####


"""
`RacahAlgebra.integralRulesForOneDjpq(rex::RacahAlgebra.RacahExpression)`
    ... attempts to find a simplification of the given Racah expression by using the normalization integration
        rule for one Wigner D-function D^j_pq(alpha,beta,gamma) = exp(-i(p*alpha+q*gamma)) d^j_pq(beta) over the
        full rotation-group volume (Varshalovich et al. 1988, Eq. 4.11.1). Once a simplification is found, no
        attempt is made to find another simplifcation for this set of rules.
        A (istrue, rex)::Tuple{Bool, RacahExpression} is returned but where rex has no meaning for !istrue.

    A Djpq carries only (j,p,q,beta); the alpha and gamma it is multiplied by exp(-i*p*alpha-i*q*gamma) in (see
    RacahAlgebra.DFunction) are NOT fields of Djpq, so unlike the Ylm rules above, this rule cannot read "its
    own" alpha/gamma off the matched symbol. Instead it searches rex.integrals STRUCTURALLY: for each periodic
    integral candidate (bounds 0 to 2*pi or 0 to 4*pi, for the extended-domain convention), it tests whether
    dividing out exp(-i*p*var) [or exp(-i*q*var)] from the weight removes that candidate's own variable --
    exactly the same division-and-hasNoVars pattern as everywhere else in this file, just applied to find WHICH
    integral plays the alpha/gamma role rather than reading it off a struct field.

    AUDITED 21-Nov-2026: numerically verified by quadrature against AngularMomentum.Wigner_dmatrix via
    RacahAlgebra.selectRacahIntegral(4); see project_racah_equivalentform_nondeterminism.md
"""
function integralRulesForOneDjpq(rex::RacahAlgebra.RacahExpression)
    # Loop through all small Wigner d-matrices
    for  (iaDjpq, aDjpq) in enumerate(rex.djpqs)
        #
        #   Rule:      2pi            2pi            pi
        #   ----      Int    dalpha  Int    dgamma  Int    dbeta  sin(beta)  exp(-i(p*alpha+q*gamma))  d^j_pq(beta)
        #               0              0              0
        #
        #           =   8 pi^2  delta(j,0)  delta(p,0)  delta(q,0)
        #
        betaIntegral = Integral(aDjpq.beta, Basic(0), SymEngine.PI)
        if  !(betaIntegral in rex.integrals)    continue    end
        remaining = RacahAlgebra.removeIntegral(betaIntegral, rex.integrals)
        #
        # Search the remaining integrals structurally for the alpha- and gamma-roles
        testWeight = rex.weight / sin(aDjpq.beta)
        alphaNtg   = nothing;    gammaNtg = nothing
        for  ntg in remaining
            if  !(ntg.low == Basic(0)  &&  (ntg.up == 2*SymEngine.PI  ||  ntg.up == 4*SymEngine.PI))    continue    end
            if  RacahAlgebra.hasNoVars([ntg.var], rex.weight)    continue    end
            if  isnothing(alphaNtg)
                probe = testWeight / exp(-im*aDjpq.p*ntg.var)
                if  RacahAlgebra.hasNoVars([ntg.var], probe)    alphaNtg = ntg;    testWeight = probe;    continue    end
            end
            if  isnothing(gammaNtg)
                probe = testWeight / exp(-im*aDjpq.q*ntg.var)
                if  RacahAlgebra.hasNoVars([ntg.var], probe)    gammaNtg = ntg;    testWeight = probe    end
            end
        end
        if  isnothing(alphaNtg)  ||  isnothing(gammaNtg)    continue    end
        #
        newDjpqs     = Djpq[];   for (ibDjpq, bDjpq) in enumerate(rex.djpqs)  if  iaDjpq != ibDjpq   push!(newDjpqs, bDjpq)   end   end
        newIntegrals = RacahAlgebra.removeIntegral(gammaNtg, RacahAlgebra.removeIntegral(alphaNtg, remaining))
        #
        if  RacahAlgebra.hasNoVars([aDjpq.beta, alphaNtg.var, gammaNtg.var], rex.phase)        &&
            RacahAlgebra.hasNoVars([aDjpq.beta, alphaNtg.var, gammaNtg.var], testWeight)       &&
            RacahAlgebra.hasNoVars([aDjpq.beta, alphaNtg.var, gammaNtg.var], rex.deltas)       &&
            RacahAlgebra.hasNoVars([aDjpq.beta, alphaNtg.var, gammaNtg.var], rex.triangles)    &&
            RacahAlgebra.hasNoVars([aDjpq.beta, alphaNtg.var, gammaNtg.var], rex.w3js)         &&
            RacahAlgebra.hasNoVars([aDjpq.beta, alphaNtg.var, gammaNtg.var], rex.w6js)         &&
            RacahAlgebra.hasNoVars([aDjpq.beta, alphaNtg.var, gammaNtg.var], rex.w9js)         &&
            RacahAlgebra.hasNoVars([aDjpq.beta, alphaNtg.var, gammaNtg.var], rex.ylms)         &&
            RacahAlgebra.hasNoVars([aDjpq.beta, alphaNtg.var, gammaNtg.var], newDjpqs)         &&
            RacahAlgebra.hasNoVars([aDjpq.beta, alphaNtg.var, gammaNtg.var], newIntegrals)
            newDeltas = rex.deltas;   push!(newDeltas, Kronecker(aDjpq.j, 0));   push!(newDeltas, Kronecker(aDjpq.p, 0));   push!(newDeltas, Kronecker(aDjpq.q, 0))
            newWeight = testWeight * 2 * (alphaNtg.up - alphaNtg.low) * (gammaNtg.up - gammaNtg.low)
            wa = RacahExpression( rex.summations, newIntegrals, rex.phase, newWeight, newDeltas, rex.triangles,
                                    rex.w3js, rex.w6js, rex.w9js, rex.ylms, newDjpqs )
            println(">> Apply integral rule for one Djpq -- Int dOmega D^j_pq ...")
            return( (true, wa) )
        end
    end

    return( (false, RacahExpression()) )
end
