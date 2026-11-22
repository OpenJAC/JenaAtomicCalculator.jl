

"""
`RacahAlgebra.differentiateYlm(ylm::RacahAlgebra.Ylm, var::Basic)`
    ... returns the derivative of the spherical harmonic ylm with respect to var, which must be either
        ylm.theta or ylm.phi (Varshalovich et al. 1988, Eq. 5.8.5):

            d/dphi   Y_lm(theta,phi)  =  i m  Y_lm(theta,phi)

            d/dtheta Y_lm(theta,phi)  =  (1/2) sqrt(l(l+1)-m(m+1))  Y_{l,m+1}(theta,phi)  exp(-i phi)
                                        - (1/2) sqrt(l(l+1)-m(m-1))  Y_{l,m-1}(theta,phi)  exp( i phi)

        Unlike the sum and integral rules in the sibling -inc- files, differentiation is not a pattern match
        that can fail -- it always succeeds for a valid (ylm, var) pair -- so there is no (istrue, rex) tuple
        here. Because d/dtheta is a genuine SUM of two terms, this returns an Array{RacahExpression,1} to be
        ADDED (unlike the same container type returned by symmetricForms, which lists EQUIVALENT forms of one
        expression, not terms to sum); a 1-element array for d/dphi, 2-element for d/dtheta.

    AUDITED 22-Nov-2026: numerically verified by finite difference against AngularMomentum.sphericalYlm for
    several (l,m,theta,phi); see project_racah_equivalentform_nondeterminism.md
"""
function differentiateYlm(ylm::Ylm, var::Basic)
    if      var == ylm.phi
            return( [ RacahExpression( Basic[], Integral[], Basic(0), im*ylm.m, Kronecker[], Triangle[], W3j[], W6j[], W9j[], [ylm], Djpq[] ) ] )
    elseif  var == ylm.theta
            term1 = RacahExpression( Basic[], Integral[], Basic(0),  sqrt(ylm.l*(ylm.l+1) - ylm.m*(ylm.m+1)) / 2 * exp(-im*ylm.phi),
                                        Kronecker[], Triangle[], W3j[], W6j[], W9j[], [Ylm(ylm.l, ylm.m+1, ylm.theta, ylm.phi)], Djpq[] )
            term2 = RacahExpression( Basic[], Integral[], Basic(0), -sqrt(ylm.l*(ylm.l+1) - ylm.m*(ylm.m-1)) / 2 * exp( im*ylm.phi),
                                        Kronecker[], Triangle[], W3j[], W6j[], W9j[], [Ylm(ylm.l, ylm.m-1, ylm.theta, ylm.phi)], Djpq[] )

            return( [term1, term2] )
    else    error("RacahAlgebra.differentiateYlm(): var must be ylm.theta or ylm.phi, not $var.")
    end
end
