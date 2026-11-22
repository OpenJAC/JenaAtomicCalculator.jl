"""
`module JAC.SphericalMultipoles`
... a submodule of JAC that builds coupled spherical harmonics -- the Clebsch-Gordan-coupled tensor products
    of ordinary spherical harmonics taken at DIFFERENT solid angles -- as RacahExpressions, ready for
    RacahAlgebra.evaluate() to simplify. These are the natural multipole building blocks for EM multipole
    fields and for nano-optics / quantum-optics radiation-pattern work (Varshalovich et al. 1988); see
    background/Manual-Jac-p.tex, section "Bipolar and tripolar spherical harmonics".

    STATUS, 22-Nov-2026: POSTPONED BY THE MAINTAINER beyond this first, deliberately small step -- this file
    implements only the BIPOLAR harmonic and its own orthogonality integral. Not implemented: the tripolar
    harmonic, the alternate coupling sequence (l1+l2)+l3 the manual mentions as a future option, and
    vector/spinor spherical harmonics. None of this is an oversight; picking any of it up is a new task, not
    a continuation of this one.
"""
module  SphericalMultipoles


using   SymEngine,  ..AngularMomentum, ..Basics,  ..Defaults,  ..RacahAlgebra

export  BipolarHarmonic


"""
`struct  SphericalMultipoles.BipolarHarmonic`
    ... defines the bipolar spherical harmonic {Y_l1(Omega1) x Y_l2(Omega2)}_LM, the irreducible tensor
        product of two ordinary spherical harmonics taken at independent solid angles (Varshalovich et al.
        1988); represented by its underlying RacahExpression, ready for RacahAlgebra.evaluate() to simplify.

    + l1, l2    ::Basic                        ... the two coupled ranks
    + L, M      ::Basic                        ... the total rank and its projection
    + rex       ::RacahAlgebra.RacahExpression ... the Clebsch-Gordan-coupled sum of Y_l1m1(Omega1) Y_l2m2(Omega2)
"""
struct   BipolarHarmonic
    l1          ::Basic
    l2          ::Basic
    L           ::Basic
    M           ::Basic
    rex         ::RacahAlgebra.RacahExpression
end


"""
`SphericalMultipoles.BipolarHarmonic(l1::Basic, l2::Basic, L::Basic, M::Basic,
                                        theta1::Basic, phi1::Basic, theta2::Basic, phi2::Basic)`
    ... constructor for the bipolar spherical harmonic {Y_l1(theta1,phi1) x Y_l2(theta2,phi2)}_LM, built
        entirely from the already-existing RacahAlgebra.ClebschGordanExpansion and RacahAlgebra.Ylm -- nothing
        here is reimplemented. A bh::SphericalMultipoles.BipolarHarmonic is returned.
"""
function BipolarHarmonic(l1::Basic, l2::Basic, L::Basic, M::Basic, theta1::Basic, phi1::Basic, theta2::Basic, phi2::Basic)
    wl1 = RacahAlgebra.AngMomentum(l1);    wl2 = RacahAlgebra.AngMomentum(l2)
    wL  = RacahAlgebra.AngMomentum(L);     wM  = RacahAlgebra.AngMomentum(M)
    # Fresh (randomly-suffixed) magnetic quantum numbers, so that two independently-built BipolarHarmonic's
    # can be multiplied together (as the orthogonality integral needs) without colliding on m1, m2
    m1  = Basic("m1_" * string(floor(Int64, 10000*rand())))
    m2  = Basic("m2_" * string(floor(Int64, 10000*rand())))

    cgRex = RacahAlgebra.ClebschGordanExpansion(wl1, m1, wl2, m2, wL, wM)
    y1Rex = RacahExpression( Basic[], RacahAlgebra.Integral[], Basic(0), Basic(1), RacahAlgebra.Kronecker[], RacahAlgebra.Triangle[],
                                RacahAlgebra.W3j[], RacahAlgebra.W6j[], RacahAlgebra.W9j[], [RacahAlgebra.Ylm(wl1, m1, theta1, phi1)], RacahAlgebra.Djpq[] )
    y2Rex = RacahExpression( Basic[], RacahAlgebra.Integral[], Basic(0), Basic(1), RacahAlgebra.Kronecker[], RacahAlgebra.Triangle[],
                                RacahAlgebra.W3j[], RacahAlgebra.W6j[], RacahAlgebra.W9j[], [RacahAlgebra.Ylm(wl2, m2, theta2, phi2)], RacahAlgebra.Djpq[] )
    rex = cgRex * y1Rex * y2Rex

    return( BipolarHarmonic(wl1, wl2, wL, wM, rex) )
end


function Base.string(bh::BipolarHarmonic)
    return( "{Y_$(bh.l1)(Omega1) x Y_$(bh.l2)(Omega2)}_{$(bh.L),$(bh.M)} = " * string(bh.rex) )
end

function Base.show(io::IO, bh::BipolarHarmonic)
    print(io, Base.string(bh))
end

end # module
