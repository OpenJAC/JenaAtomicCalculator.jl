
"""
`module JAC.RacahAlgebra`  
... a submodel of JAC that contains procedures for defining a Nuclear.Model and for calculating various nuclear
    potentials.
"""
module  RacahAlgebra


using   SymEngine,  ..AngularMomentum, ..Basics,  ..Defaults

export  ClebschGordan, Integral, Kronecker, Triangle, W3j, W6j, W9j, Ylm, Djpq, RacahExpression


"""
`struct RacahAlgebra.AngMomentum` 
    ... defines an (abstract) data types for symbolic angular momenta which accept the types 
        Basic, Symbol, Int64 and Rational{Int64} and check for being consistent with angular momenta.
"""
struct  AngMomentum   end

function AngMomentum(x::Basic)   return(x)                                end
function AngMomentum(x::Int64)   return( Basic(x) )                       end
function AngMomentum(x::Symbol)  return( SymEngine.symbols(string(x)) )   end
function AngMomentum(x::Rational{Int64})  
    if  x.den == 1 || x.den == 2     return( Basic(x.num//x.den) )   
    else    error("Angular momenta must be integer or half-interger.")
    end
end


"""
`struct  RacahAlgebra.Kronecker`  ... defines a type for a Kronecker delta symbol with symbolic arguments.

    + i, k     ::Basic   ... Kronecker indices
"""
struct  Kronecker
    i                ::Basic
    k                ::Basic
end


"""
`RacahAlgebra.Kronecker(i::AngMomentum, k::AngMomentum)`  
    ... constructor for defining a Kronecker delta symbol either by Julia Symbol's or SymEngine Basic variables.
"""
function Kronecker(i::AngMomentum, k::AngMomentum)
    wi = AngMomentum(i);    wk = AngMomentum(k)
    Kronecker(wi,wk)
end


# `Base.show(io::IO, delta::Kronecker)`  ... prepares a proper printout of the variable delta::Kronecker.
    function Base.show(io::IO, delta::Kronecker)
    print(io, "delta($(delta.i), $(delta.k))")
end


"""
`Base.:(==)(wa::Kronecker, wb::Kronecker)`  
    ... compares two (symbolic) Kronecker deltas and return true if all subfields are equal under permutation, and false otherwise.
"""
function  Base.:(==)(wa::Kronecker, wb::Kronecker)
    if  (wa.i == wb.i  &&  wa.k == wb.k)   ||    (wa.i == wb.k  &&  wa.k == wb.i)    return( true )
    else                                                                             return( false )   
    end
end


"""
`struct  RacahAlgebra.Triangle`  ... defines a type for a Triangle delta symbol with symbolic arguments.

    + ja, jb, jc     ::Basic   ... angular momenta in the Triangle delta
"""
struct  Triangle
    ja               ::Basic
    jb               ::Basic
    jc               ::Basic
end


"""
`RacahAlgebra.Triangle(ja::AngMomentum, jb::AngMomentum, jc::AngMomentum)`  
    ... constructor for defining a Triangle delta symbol either by Julia Symbol's or SymEngine Basic variables.
"""
function Triangle(ja::AngMomentum, jb::AngMomentum, jc::AngMomentum)
    wja = AngMomentum(ja);    wjb = AngMomentum(jb);    wjc = AngMomentum(jc)    
    Triangle(wja,wjb,wjc)
end


# `Base.show(io::IO, delta::Triangle)`  ... prepares a proper printout of the variable delta::Triangle.
    function Base.show(io::IO, delta::Triangle)
    print(io, "delta ($(delta.ja), $(delta.jb), $(delta.jc))")
end


"""
`Base.:(==)(wa::Triangle, wb::Triangle)`  
    ... compares two (symbolic) Triangle deltas and return true if all subfields are equal under permutation, 
        and false otherwise.
"""
function  Base.:(==)(wa::Triangle, wb::Triangle)
    wasymbols = [wa.ja, wa.jb, wa.jc];    wbsymbols = [wb.ja, wb.jb, wb.jc]
    for  wasb in wasymbols  if  !(wasb in wbsymbols)    return(false)   end     end
    for  wbsb in wbsymbols  if  !(wbsb in wasymbols)    return(false)   end     end
    return( true )
end


"""
`struct  RacahAlgebra.W3j`  ... defines a type for a Wigner 3-j symbol with symbolic arguments.

    + ja, jb, jc     ::Basic   ... angular momenta
    + ma, mb, mc     ::Basic   ... projections of the angular momenta above
"""
struct  W3j
    ja               ::Basic
    jb               ::Basic
    jc               ::Basic
    ma               ::Basic
    mb               ::Basic
    mc               ::Basic
end


"""
`RacahAlgebra.W3j(ja::AngMomentum, jb::AngMomentum, jc::AngMomentum, ma::AngMomentum, mb::AngMomentum, mc::AngMomentum)`  
    ... constructor for defining the Wigner 3-j symbol either by Julia Symbol's or SymEngine Basic variables.
"""
function W3j(ja::AngMomentum, jb::AngMomentum, jc::AngMomentum, 
                ma::AngMomentum, mb::AngMomentum, mc::AngMomentum)
    wja = AngMomentum(ja);    wjb = AngMomentum(jb);    wjc = AngMomentum(jc)    
    wma = AngMomentum(ma);    wmb = AngMomentum(mb);    wmc = AngMomentum(mc)    
    W3j(wja, wjb, wjc, wma, wmb, wmc)
end


# `Base.show(io::IO, w::W3j)`  ... prepares a proper printout of the variable  w3j::W3j.
    function Base.show(io::IO, w::W3j)
    print(io, "W3j($(w.ja), $(w.jb), $(w.jc); $(w.ma), $(w.mb), $(w.mc))")
end


"""
`Base.:(==)(wa::W3j, wb::W3j)`  
    ... compares two (symbolic) Wigner 3-j symbols and return true if all subfields are equal, and false otherwise.
"""
function  Base.:(==)(wa::W3j, wb::W3j)
    if   wa.ja != wb.ja   ||   wa.jb != wb.jb   ||   wa.jc != wb.jc   ||   
            wa.ma != wb.ma   ||   wa.mb != wb.mb   ||   wa.mc != wb.mc       return( false )    
    else                                                                  return( true )    
    end
end


"""
`struct  RacahAlgebra.W6j`  ... defines a type for a Wigner 6-j symbol with symbolic arguments.

    + a, b, c, d, e, f     ::Basic   ... angular momenta
"""
struct  W6j
    a                      ::Basic
    b                      ::Basic
    c                      ::Basic
    d                      ::Basic
    e                      ::Basic
    f                      ::Basic
end


"""
`RacahAlgebra.W6j(a::AngMomentum, b::AngMomentum, c::AngMomentum, 
                    d::AngMomentum, e::AngMomentum, f::AngMomentum)`  
    ... constructor for defining the Wigner 6-j symbol either by Julia Symbol's or SymEngine Basic variables.
"""
function  W6j(a::AngMomentum, b::AngMomentum, c::AngMomentum, 
                d::AngMomentum, e::AngMomentum, f::AngMomentum)
    wa = AngMomentum(a);    wb = AngMomentum(b);    wc = AngMomentum(c)    
    wd = AngMomentum(d);    we = AngMomentum(e);    wf = AngMomentum(f)    
    println("aaa")
    W6j(wa, wb, wc, wd, we, wf)
end


# `Base.show(io::IO, w::W6j)`  ... prepares a proper printout of the variable  w6j::W6j.
function Base.show(io::IO, w::W6j)
    print(io, "W6j{$(w.a), $(w.b), $(w.c); $(w.d), $(w.e), $(w.f)}")
end


"""
`Base.:(==)(wa::W6j, wb::W6j)`  
    ... compares two (symbolic) Wigner 6j symbols and return true if all subfields are equal, and false otherwise.
"""
function  Base.:(==)(wa::W6j, wb::W6j)
    if   wa.a != wb.a   ||   wa.b != wb.b   ||   wa.c != wb.c   ||   
            wa.d != wb.d   ||   wa.e != wb.e   ||   wa.f != wb.f             return( false )    
    else                                                                  return( true )    
    end
end


"""
`struct  RacahAlgebra.W9j`  ... defines a type for a Wigner 9j symbol with symbolic arguments.

    + a, b, c, d, e, f, g, h, i     ::Basic   ... angular momenta
"""
struct  W9j
    a                      ::Basic
    b                      ::Basic
    c                      ::Basic
    d                      ::Basic
    e                      ::Basic
    f                      ::Basic
    g                      ::Basic
    h                      ::Basic
    i                      ::Basic
end


"""
`RacahAlgebra.W9j(a::AngMomentum, b::AngMomentum, c::AngMomentum, 
                    d::AngMomentum, e::AngMomentum, f::AngMomentum,
                    g::AngMomentum, h::AngMomentum, i::AngMomentum)`  
    ... constructor for defining the Wigner 9j symbol either by Julia Symbol's or SymEngine Basic variables.
"""
function  W9j(a::AngMomentum, b::AngMomentum, c::AngMomentum, 
                d::AngMomentum, e::AngMomentum, f::AngMomentum,
                g::AngMomentum, h::AngMomentum, i::AngMomentum)
    wa = AngMomentum(a);    wb = AngMomentum(b);    wc = AngMomentum(c)    
    wd = AngMomentum(d);    we = AngMomentum(e);    wf = AngMomentum(f)    
    wg = AngMomentum(g);    wh = AngMomentum(h);    wi = AngMomentum(i)    
    W9j(wa, wb, wc, wd, we, wf, wg, wh, wi)
end


# `Base.show(io::IO, w::W9j)`  ... prepares a proper printout of the variable  w::W9j.
function Base.show(io::IO, w::W9j)
    print(io, "W9j{$(w.a), $(w.b), $(w.c); $(w.d), $(w.e), $(w.f); $(w.g), $(w.h), $(w.i)}")
end


"""
`Base.:(==)(wa::W9j, wb::W9j)`  
    ... compares two (symbolic) Wigner 6j symbols and return true if all subfields are equal, and false otherwise.
"""
function  Base.:(==)(wa::W9j, wb::W9j)
    if   wa.a != wb.a   ||   wa.b != wb.b   ||   wa.c != wb.c   ||   
            wa.d != wb.d   ||   wa.e != wb.e   ||   wa.f != wb.f   ||   
            wa.g != wb.g   ||   wa.h != wb.h   ||   wa.i != wb.i             return( false )    
    else                                                                  return( true )    
    end
end


"""
`struct  RacahAlgebra.Ylm`  ... defines a type for a spherical harmonic Y_lm (theta,phi) with symbolic arguments.

    + l, m                 ::Basic   ... angular momenta
    + theta, phi           ::Basic   ... polar and azimuthal angle
"""
struct  Ylm
    l                      ::Basic
    m                      ::Basic
    theta                  ::Basic
    phi                    ::Basic
end


"""
`RacahAlgebra.Ylm(l::AngMomentum, m::AngMomentum, theta::Basic, phi::Basic, star::Bool)`  
    ... constructor for defining the spherical harmonic Y_lm (theta,phi) either by Julia Symbol's or SymEngine Basic variables.
        The variable star indicates a complex-conjugate function, and which case a Racah expression is returned.
        
    + (l::Basic, m::Basic, theta::Basic, phi::Basic, star::Bool)
"""
function  Ylm(l::AngMomentum, m::AngMomentum, theta::Basic, phi::Basic, star::Bool=false)
    wl = AngMomentum(l);    wm = AngMomentum(wm)  
    if  star    rex = RacahExpression( Basic[], Integral[], -wm, Basic(1), Kronecker[], Triangle[], W3j[], W6j[], W9j[], [Ylm(wl, -wm, theta, phi)], Djpq[] )
                return( rex )
    else        return( Ylm(wl, wm, theta, phi) )
    end
end


function  Ylm(wl::Basic, wm::Basic, theta::Basic, phi::Basic, star::Bool)
    if  star    rex = RacahExpression( Basic[], Integral[], -wm, Basic(1), Kronecker[], Triangle[], W3j[], W6j[], W9j[], [Ylm(wl, -wm, theta, phi)], Djpq[] )
                return( rex )
    else        return( Ylm(wl, wm, theta, phi) )
    end
end


# `Base.show(io::IO, y::Ylm)`  ... prepares a proper printout of the variable  ylm::Ylm.
function Base.show(io::IO, y::Ylm)
    print(io, "Y_{$(y.l), $(y.m)} ($(y.theta), $(y.phi))")
end


"""
`Base.:(==)(ya::Ylm, yb::Ylm)`  
    ... compares two (symbolic) spherical harmonics Y_lm (theta,phi) and return true if all subfields are equal, 
        and false otherwise.
"""
function  Base.:(==)(ya::Ylm, yb::Ylm)
    if   ya.l != yb.l   ||   ya.m != yb.m   ||   ya.theta != yb.theta   ||   ya.phi != yb.phi    return( false )    
    else                                                                                         return( true )    
    end
end


"""
`struct  RacahAlgebra.Djpq`  
    ... defines a type for the (small and real !!) Wigner rotation matrix d^(j)_pq (beta) with symbolic arguments.

    + j, p, q              ::Basic   ... angular momenta
    + beta                 ::Basic   ... rotation angle
"""
struct  Djpq
    j                      ::Basic
    p                      ::Basic
    q                      ::Basic
    beta                   ::Basic
end


"""
`RacahAlgebra.Djpq(j::AngMomentum, p::AngMomentum, q::AngMomentum, beta::Basic)`  
    ... constructor for defining the spherical harmonic Y_lm (theta,phi) either by Julia Symbol's or SymEngine Basic variables.
"""
function  Djpq(j::AngMomentum, p::AngMomentum, q::AngMomentum, beta::Basic)
    wj = AngMomentum(j);    wp = AngMomentum(p);    wq = AngMomentum(q)
    Djpq(j, p, q, beta)
end


# `Base.show(io::IO, d::Djpq)`  ... prepares a proper printout of the variable  d::Djpq.
function Base.show(io::IO, d::Djpq)
    print(io, "d^($(d.j))_{$(d.p), $(d.q)} ($(d.beta))")
end


"""
`Base.:(==)(da::Djpq, db::Djpq)`  
    ... compares two (symbolic) small Wigner rotation matrices d^(j)_pq (beta) and return true if all subfields are equal, 
        and false otherwise.
"""
function  Base.:(==)(da::Djpq, db::Djpq)
    if   da.j != db.j   ||   da.p != db.p   ||   da.q != db.q   ||   da.beta != db.beta     return( false )    
    else                                                                                                               return( true )    
    end
end


"""
`struct  RacahAlgebra.Integral`  
    ... defines an (finite) integral, typically over some angle, Int_low^up d var with symbolic arguments.

    + var                  ::Basic   ... integration variable, typically an angle (theta, phi, beta, ...)
    + low, up              ::Basic   ... lower, upper integration bound
"""
struct  Integral
    var                    ::Basic
    low                    ::Basic
    up                     ::Basic
end


# `Base.show(io::IO, ntg::Integral)`  ... prepares a proper printout of the variable  ntg::Integral.
function Base.show(io::IO, ntg::Integral)
    print(io, "Int_$(ntg.low)^$(ntg.up) d$(ntg.var)")
end


"""
`Base.:(==)(ntga::Integral, ntgb::Integral)`  
    ... compares two (symbolic) integrals and return true if all subfields are equal, and false otherwise.
"""
function  Base.:(==)(ntga::Integral, ntgb::Integral)
    if   ntga.var != ntgb.var   ||   ntga.low != ntgb.low   ||   ntga.up != ntgb.up     return( false )    
    else                                                                                return( true )    
    end
end


"""
`struct  RacahAlgebra.RacahExpression`  ... defines a type for a RacahExpression with symbolic arguments.

    + summations    ::Array{Basic,1}      ... Summation indices.
    + integrals     ::Array{Integral,1}   ... Integrals.
    + phase         ::Basic               ... Phase of the Racah expression.
    + weight        ::Basic               ... Weight of the Racah expression.
    + deltas        ::Array{Kronecker,1}  ... List of Kronecker deltas.
    + triangles     ::Array{Triangle,1}   ... List of Triangle deltas.
    + w3js          ::Array{W3j,1}        ... List of Wigner 3-j symbols
    + w6js          ::Array{W6j,1}        ... List of Wigner 6-j symbols
    + w9js          ::Array{W9j,1}        ... List of Wigner 9-j symbols
    + ylms          ::Array{Ylm,1}        ... List of spherical harmonics
    + djpqs         ::Array{Djpq,1}       ... List of (small) Wigner rotation matrices
"""
struct  RacahExpression
    summations             ::Array{Basic,1}
    integrals              ::Array{Integral,1}
    phase                  ::Basic
    weight                 ::Basic
    deltas                 ::Array{Kronecker,1}
    triangles              ::Array{Triangle,1}
    w3js                   ::Array{W3j,1} 
    w6js                   ::Array{W6j,1}
    w9js                   ::Array{W9j,1}
    ylms                   ::Array{Ylm,1}
    djpqs                  ::Array{Djpq,1} 
end


"""
`RacahAlgebra.RacahExpression()`  
    ... constructor for defining an empty RacahExpression.
"""
function  RacahExpression()
    RacahExpression( Basic[], Integral[], Basic(0), Basic(1), Kronecker[], Triangle[], W3j[], W6j[], W9j[], Ylm[], Djpq[] )
end


"""
`RacahAlgebra.RacahExpression(rex::RacahAlgebra.RacahExpression;`

        summations=..,      integrals=..,   phase=..,       weight=..,      deltas=..,     triangles=..,       
        w3js=..,            w6js=..,        w9js=..,        ylms=..,       djpqs=..) 
                    
    ... constructor for modifying a given rex::RacahExpression by 'overwriting' the explicitly selected parts of the 
        expression.
"""
function RacahExpression(rex::RacahAlgebra.RacahExpression;    
                            summations::Union{Nothing,Array{Basic,1}}=nothing,        integrals::Union{Nothing,Array{Integral,1}}=nothing,        
                            phase::Union{Nothing,Basic}=nothing,            
                            weight::Union{Nothing,Basic}=nothing,                     deltas::Union{Nothing,Array{Kronecker,1}}=nothing,     
                            triangles::Union{Nothing,Array{Triangle,1}}=nothing,      w3js::Union{Nothing,Array{W3j,1}}=nothing,             
                            w6js::Union{Nothing,Array{W6j,1}}=nothing,                w9js::Union{Nothing,Array{W9j,1}}=nothing,
                            ylms::Union{Nothing,Array{Ylm,1}}=nothing,                djpqs::Union{Nothing,Array{Djpq,1}}=nothing) 
    
    if  summations == nothing   summationsx = rex.summations    else   summationsx = summations     end 
    if  integrals  == nothing   integralsx  = rex.integrals     else   integralsx  = integrals      end 
    if  phase      == nothing   phasex      = rex.phase         else   phasex = phase               end 
    if  weight     == nothing   weightx     = rex.weight        else   weightx = weight             end 
    if  deltas     == nothing   deltasx     = rex.deltas        else   deltasx = deltas             end 
    if  triangles  == nothing   trianglesx  = rex.triangles     else   trianglesx = triangles       end 
    if  w3js       == nothing   w3jsx       = rex.w3js          else   w3jsx  = w3js                end 
    if  w6js       == nothing   w6jsx       = rex.w6js          else   w6jsx  = w6js                end 
    if  w9js       == nothing   w9jsx       = rex.w9js          else   w9jsx  = w9js                end 
    if  ylms       == nothing   ylmsx       = rex.ylms          else   ylmsx  = ylms                end 
    if  djpqs      == nothing   djpqsx      = rex.djpqs         else   djpqsx = djpqs               end 
    
    RacahExpression( summationsx, integralsx, phasex, weightx, deltasx, trianglesx, w3jsx, w6jsx, w9jsx, ylmsx, djpqsx)
end



# `Base.show(io::IO, rex::RacahExpression)`  ... prepares a proper printout of the variable  rex::RacahExpression.
function Base.show(io::IO, rex::RacahExpression)
    sa = ""
    if  rex.summations != Basic[]    sa = sa * "Sum_[$(rex.summations)]  "      end
    for  ntg      in rex.integrals   sa = sa * "$ntg   "                        end
    if  rex.phase      != Basic(0)   sa = sa * "(-1)^($(rex.phase))  "          end
    if  rex.weight     != Basic(1)   sa = sa * "($(rex.weight))  "              end
    for  delta    in rex.deltas      sa = sa * "$delta  "                       end
    for  triangle in rex.triangles   sa = sa * "$triangle  "                    end
    for  w3j      in rex.w3js        sa = sa * "$w3j  "                         end
    for  w6j      in rex.w6js        sa = sa * "$w6j  "                         end
    for  w9j      in rex.w9js        sa = sa * "$w9j  "                         end
    for  ylm      in rex.ylms        sa = sa * "$ylm  "                         end
    for  djpq     in rex.djpqs       sa = sa * "$djpq  "                        end
    if sa == ""                      sa = "1"                                   end
    print(io, sa)
end


function  Base.:(*)(rexa::RacahExpression, rexb::RacahExpression)
    newSummations = rexa.summations;    for  su in rexb.summations         push!(newSummations, su)         end
    newIntegrals  = rexa.integrals;     for  ntg  in rexb.integrals        push!(newIntegrals, ntg)         end
    newPhase      = rexa.phase  + rexb.phase
    newWeight     = rexa.weight * rexb.weight
    newDeltas     = rexa.deltas;        for  delta in rexb.deltas          push!(newDeltas, delta)          end
    newTriangles  = rexa.triangles;     for  triangle in rexb.triangles    push!(newTriangles, triangle)    end
    newW3js       = rexa.w3js;          for  w3j  in rexb.w3js             push!(newW3js, w3j)              end
    newW6js       = rexa.w6js;          for  w6j  in rexb.w6js             push!(newW6js, w6j)              end
    newW9js       = rexa.w9js;          for  w9j  in rexb.w9js             push!(newW9js, w9j)              end
    newYlms       = rexa.ylms;          for  ylm  in rexb.ylms             push!(newYlms, ylm)              end
    newDjpqs      = rexa.djpqs;         for  djpq in rexb.djpqs            push!(newDjpqs, djpq)            end
    
    return( RacahExpression(newSummations, newIntegrals, newPhase, newWeight, newDeltas, newTriangles, newW3js, newW6js, newW9js, newYlms, newDjpqs) )
end


"""
`Base.:(==)(rexa::RacahExpression, rexb::RacahExpression)`  
    ... compares two (symbolic) Racah expressions and return true if all subfields are equal, and false otherwise.
"""
function  Base.:(==)(rexa::RacahExpression, rexb::RacahExpression)
    if   rexa.phase  != rexb.phase   ||   rexa.weight  != rexb.weight       return( false )     end

    if   length(rexa.summations)     !=   length(rexb.summations)    ||
            length(rexa.integrals)      !=   length(rexb.integrals)     ||
            length(rexa.deltas)         !=   length(rexb.deltas)        ||
            length(rexa.triangles)      !=   length(rexb.triangles)     ||
            length(rexa.w3js)           !=   length(rexb.w3js)          ||
            length(rexa.w6js)           !=   length(rexb.w6js)          ||
            length(rexa.w9js)           !=   length(rexb.w9js)          ||
            length(rexa.ylms)           !=   length(rexb.ylms)          ||
            length(rexa.djpqs)          !=   length(rexb.djpqs)                return( false )     end
    
    for  i = 1:length(rexa.summations)
        if  rexa.summations[i] != rexb.summations[i]                        return( false )     end     
    end
    for  i = 1:length(rexa.integrals)
        if  rexa.integrals[i] != rexb.integrals[i]                          return( false )     end     
    end
    for  i = 1:length(rexa.deltas)
        if  rexa.deltas[i] != rexb.deltas[i]                                return( false )     end     
    end
    for  i = 1:length(rexa.triangles)
        if  rexa.triangles[i] != rexb.triangles[i]                          return( false )     end     
    end
    
    for  i = 1:length(rexa.w3js)
        if  rexa.w3js[i] != rexb.w3js[i]                                    return( false )     end     
    end
    for  i = 1:length(rexa.w6js)
        if  rexa.w6js[i] != rexb.w6js[i]                                    return( false )     end     
    end
    for  i = 1:length(rexa.w9js)
        if  rexa.w9js[i] != rexb.w9js[i]                                    return( false )     end     
    end
    for  i = 1:length(rexa.ylms)
        if  rexa.ylms[i] != rexb.ylms[i]                                    return( false )     end     
    end
    for  i = 1:length(rexa.djpqs)
        if  rexa.djpqs[i] != rexb.djpqs[i]                                  return( false )     end     
    end
    
    return( true )    
end


"""
`struct  RacahAlgebra.Csq`  
    ... defines a type for a coupling sequence of two angular momenta (a,b) c   or  a + b -> c  with symbolic arguments.

    + a    ::Union{Basic,RacahAlgebra.Csq}      ... First angular momentum or coupling sequence.
    + b    ::Union{Basic,RacahAlgebra.Csq}      ... Second angular momentum or coupling sequence.
    + c    ::Basic                              ... Angular momentum to which a + b is coupled.
"""
struct  Csq
    a      ::Union{Basic,RacahAlgebra.Csq} 
    b      ::Union{Basic,RacahAlgebra.Csq} 
    c      ::Basic
end



# `Base.show(io::IO, csq::RacahAlgebra.Csq)`  ... prepares a proper printout of the variable  csq::RacahAlgebra.Csq.
function Base.show(io::IO, csq::RacahAlgebra.Csq)
    sa = "($(csq.a), $(csq.b)) $(csq.c)"
    print(io, sa)
end


"""
`abstract type RacahAlgebra.AbstractRecursionW3j` 
    ... defines an abstract and a number of singleton types for the recursion rules for Wigner 3-j symbols.

    + RecursionW3jMagnetic      ... Recursion wrt. the magnetic quantum numbers.
    + RecursionW3jOneStep       ... Recursion with step-1 of the j-quantum numbers.
    + RecursionW3jHalfStep      ... Recursion with step-1/2 of the j-quantum numbers.
    + RecursionW3jLouck         ... Recursion wrt. j-quantum numbers due to Louck.
"""
abstract type  AbstractRecursionW3j                             end
struct         RecursionW3jMagnetic  <:  AbstractRecursionW3j   end
struct         RecursionW3jOneStep   <:  AbstractRecursionW3j   end
struct         RecursionW3jHalfStep  <:  AbstractRecursionW3j   end
struct         RecursionW3jLouck     <:  AbstractRecursionW3j   end


#############################################################################################################################
#############################################################################################################################
#############################################################################################################################


"""
`RacahAlgebra.ClebschGordan(ja::Basic, ma::Basic, jb::Basic, mb::Basic, jc::Basic, mc::Basic)`  
    ... returns the cg::RacahExpression for a standard Clebsch-Gordan coefficient either by Julia Symbol's or SymEngine 
        Basic variables.
"""
function ClebschGordan(ja::Basic, ma::Basic, jb::Basic, mb::Basic, jc::Basic, mc::Basic)
    wja = AngMomentum(ja);    wjb = AngMomentum(jb);    wjc = AngMomentum(jc)    
    wma = AngMomentum(ma);    wmb = AngMomentum(mb);    wmc = AngMomentum(mc)    
    w3j = W3j(wja, wjb, wjc, wma, wmb, -wmc)
    rex = RacahExpression( Basic[], Integral[], wja - wjb + wmc, sqrt(2*wjc+1), Kronecker[], Triangle[], [w3j], W6j[], W9j[], Ylm[], Djpq[] )
    return( rex )
end


"""
`RacahAlgebra.ClebschGordanExpansion(ja::Basic, ma::Basic, jb::Basic, mb::Basic, jc::Basic, mc::Basic)`  
    ... returns the expansion rex::RacahExpression for a standard Clebsch-Gordan expansion  
        |jc, mc > = Sum(ma, mb) |ja, ma > |jb, mb >  <ja ma, jb mb| jc mc>  either by Julia Symbol's or SymEngine 
        Basic variables.
"""
function ClebschGordanExpansion(ja::Basic, ma::Basic, jb::Basic, mb::Basic, jc::Basic, mc::Basic)
    wja = AngMomentum(ja);    wjb = AngMomentum(jb);    wjc = AngMomentum(jc)    
    wma = AngMomentum(ma);    wmb = AngMomentum(mb);    wmc = AngMomentum(mc)    
    w3j = W3j(wja, wjb, wjc, wma, wmb, -wmc)
    rex = RacahExpression( [wma, wmb], Integral[], wja - wjb + wmc, sqrt(2*wjc+1), Kronecker[], Triangle[], [w3j], W6j[], W9j[], Ylm[], Djpq[] )
    return( rex )
end


"""
`RacahAlgebra.countWignerSymbols(rex::RacahExpression)`  
    ... counts the (total) number of Wigner nj symbols in rex. A nwnjs::Int64 is returned.
"""
function countWignerSymbols(rex::RacahExpression)
    nwnjs = length(rex.w3js) + length(rex.w6js) + length(rex.w9js)
    return( nwnjs )
end


"""
`RacahAlgebra.DFunction(j::Basic, p::Basic, q::Basic, alpha::Basic, beta::Basic, gamma::Basic)`  
    ... returns the rex::RacahExpression for a Wigner D-function D^(j)_pq (alpha, beta, gamma) by SymEngine Basic variables.
"""
function DFunction(j::Basic, p::Basic, q::Basic, alpha::Basic, beta::Basic, gamma::Basic)
    rex = RacahExpression( Basic[], Integral[], Basic(0), exp(-im*p*alpha -im*q*gamma), Kronecker[], Triangle[], W3j[], W6j[], W9j[], Ylm[], [ Djpq(j,p,q,beta) ] )
    return( rex )
end


"""
`RacahAlgebra.equivalentForm(w3j::RacahAlgebra.W3j; regge::Bool=false)`  
    ... generates an (random) equivalent form for the Wigner 3j symbol w3j by using either the classical
        (regge = false) or Regge symmetries (regge = true). A rex:RacahExpression is returned.
"""
function equivalentForm(w3j::RacahAlgebra.W3j; regge::Bool=false)
    wa  = RacahAlgebra.symmetricForms(w3j, regge=regge)
    if  regge   n = rand(1:72)    else     n = rand(1:12)   end
    println(">> Select $(n)th equivalent form for $w3j    ==>   $( wa[n])")
    return( wa[n] )
end


"""
`RacahAlgebra.equivalentForm(w6j::RacahAlgebra.W6j; regge::Bool=false)`  
    ... generates an (random) equivalent form for the Wigner 6j symbol w6j by using either the classical
        (regge = false) or Regge symmetries (regge = true). A rex:RacahExpression is returned.
"""
function equivalentForm(w6j::RacahAlgebra.W6j; regge::Bool=false)
    wa  = RacahAlgebra.symmetricForms(w6j, regge=regge)
    if  regge   n = rand(1:144)    else     n = rand(1:24)   end
    println(">> Select $(n)th equivalent form for $w6j    ==>   $( wa[n])")
    return( wa[n] )
end


"""
`RacahAlgebra.equivalentForm(w9j::RacahAlgebra.W9j; regge::Bool=false)`  
    ... generates an (random) equivalent form for the Wigner 9-j symbol w9j by using either the classical
        (regge = false) or Regge symmetries (regge = true). A rex:RacahExpression is returned.
"""
function equivalentForm(w9j::RacahAlgebra.W9j; regge::Bool=false)
    wa  = RacahAlgebra.symmetricForms(w9j, regge=regge)
    if  regge   n = rand(1:72)    else     n = rand(1:72)   end
    println(">> Select $(n)th equivalent form for $w9j    ==>   $( wa[n])")
    return( wa[n] )
end


"""
`RacahAlgebra.equivalentForm(rex::RacahAlgebra.RacahExpression; regge::Bool=false)`  
    ... generates an (random) equivalent form for the Racah expression rex by using either the classical
        (regge = false) or Regge symmetries (regge = true) for the Wigner 3j, 6j or 9j symbols. 
        A rex:RacahExpression is returned.
"""
function equivalentForm(rex::RacahAlgebra.RacahExpression; regge::Bool=false)
    newRex = rex;   newPhase = newRex.phase;    newW3js  = W3j[] 
    # Generate random equivalent forms of all Wigner 3j symbols
    for (iaW3j, aW3j)  in  enumerate(newRex.w3js)
        xaRex    = RacahAlgebra.equivalentForm(aW3j, regge=regge)
        newPhase = newPhase + xaRex.phase
        push!( newW3js, xaRex.w3js[1])
    end
    newRex = RacahExpression( newRex.summations, newRex.integrals, newPhase, newRex.weight, newRex.deltas, 
                                newRex.triangles,  newW3js, newRex.w6js, newRex.w9js, newRex.ylms, newRex.djpqs )
    
    # Generate random equivalent forms of all Wigner 6j symbols
    newPhase = newRex.phase;    newW6js  = W6j[] 
    for (iaW6j, aW6j)  in  enumerate(newRex.w6js)
        xaRex    = RacahAlgebra.equivalentForm(aW6j, regge=regge)
        newPhase = newPhase + xaRex.phase
        push!( newW6js, xaRex.w6js[1])
    end
    newRex = RacahExpression( newRex.summations, newRex.integrals, newPhase, newRex.weight, newRex.deltas, 
                                newRex.triangles,  newRex.w3js, newW6js, newRex.w9js, newRex.ylms, newRex.djpqs )
                                
    # Generate random equivalent forms of all Wigner 9-j symbols
    newPhase = newRex.phase;    newW9js  = W9j[] 
    for (iaW9j, aW9j)  in  enumerate(newRex.w9js)
        xaRex    = RacahAlgebra.equivalentForm(aW9j, regge=regge)
        newPhase = newPhase + xaRex.phase
        push!( newW9js, xaRex.w9js[1])
    end
    newRex = RacahExpression( newRex.summations, newRex.integrals, newPhase, newRex.weight, newRex.deltas, 
                                newRex.triangles,  newRex.w3js, newRex.w6js, newW9js, newRex.ylms, newRex.djpqs )
    
    
    return( newRex )
end


"""
`RacahAlgebra.evaluate(wj::Union{W3j,W6j,W9j})`  
    ... attempts to evaluate a symbolic Wigner 3j, 6j or 9j symbols by means of special values.
"""
function  evaluate(wj::Union{W3j,W6j,W9j})
    wa = RacahAlgebra.specialValue(wj)
    if    wa[1]   println(">> Special value found for  $wj = $(wa[2]).");     return( wa[2] )
    else          println(">> No special value found for  $(wj).");           return( nothing )
    end
    
end


"""
`RacahAlgebra.evaluate(ylm::Ylm)`  
    ... attempts to evaluate a symbolic spherical harmonic by means of special values.
"""
function  evaluate(ylm::Ylm)
    error("not yet implemented")
    wa = RacahAlgebra.specialValue(ylm)
    if    wa[1]   println(">> Special value found for  $wj = $(wa[2]).");     return( wa[2] )
    else          println(">> No special value found for  $(wj).");           return( nothing )
    end
    
end


"""
`RacahAlgebra.evaluate(djpq::Djpq)`  
    ... attempts to evaluate a symbolic (small) Wigner rotation matrix by means of special values.
"""
function  evaluate(djpq::Djpq)
    error("not yet implemented")
    wa = RacahAlgebra.specialValue(djpq)
    if    wa[1]   println(">> Special value found for  $wj = $(wa[2]).");     return( wa[2] )
    else          println(">> No special value found for  $(wj).");           return( nothing )
    end
    
end


"""
`RacahAlgebra.evaluate(rex::RacahExpression; special::Bool=false)`  
    ... attempts to evaluate and symbolically simplify a Racah expression by means of special values, if
        special=true, or by sum rules. A newrex::RacahExpression is returned once a (single) simplification 
        has been found, and nothing otherwise. No attempt is presently made to find further simplication, 
        once a rule has been applied.
"""
function  evaluate(rx::RacahExpression; special::Bool=false)
    # Discussion(!!)  It's not fully clear why rx/rex is modified by this procedure without such a copy.
    rex = RacahExpression( copy(rx.summations), copy(rx.integrals), copy(rx.phase), copy(rx.weight), 
                            copy(rx.deltas), copy(rx.triangles), copy(rx.w3js), copy(rx.w6js), copy(rx.w9js), copy(rx.ylms), copy(rx.djpqs) )
    if  special
        # Simplify by means of special values if this is requested
        for  (iaW3j, aW3j) in enumerate(rex.w3js)    wa = evaluate(aW3j)
            if wa != nothing    newW3js = W3j[]     
                for  (ibW3j, bW3j) in enumerate(rex.w3js)   if  iaW3j == ibW3j  else  push!(newW3js, ibW3j)   end    end
                rrex = RacahExpression( rex.summations, rex.integrals, rex.phase, rex.weight, rex.deltas, rex.triangles, newW3js, rex.w6js, rex.w9js, rex.ylms, rex.djpqs)
                return( rrex * wa )
            end
        end
        for  (iaW6j, aW6j) in enumerate(rex.w6js)    wa = evaluate(aW6j)
            if wa != nothing    newW6js = W6j[]     
                for  (ibW6j, bW6j) in enumerate(rex.w6js)   if  iaW6j == ibW6j  else  push!(newW6js, ibW6j)   end    end
                rrex = RacahExpression( rex.summations, rex.integrals, rex.phase, rex.weight, rex.deltas, rex.triangles, rex.w3js, newW6js, rex.w9js, rex.ylms, rex.djpqs)
                return( rrex * wa )
            end
        end
        for  (iaW9j, aW9j) in enumerate(rex.w9js)    wa = evaluate(aW9j)
            if wa != nothing    newW9js = W9j[]     
                for  (ibW9j, bW9j) in enumerate(rex.w9js)   if  iaW9j == ibW9j  else  push!(newW9js, ibW9j)   end    end
                rrex = RacahExpression( rex.summations, rex.integrals, rex.phase, rex.weight, rex.deltas, rex.triangles, rex.w3js, rex.w6js, newW9js, rex.ylms, rex.djpqs)
                return( rrex * wa )
            end
        end
    else
        newrex = RacahExpression( rex.summations, rex.integrals, RacahAlgebra.purifyPhase(rex.phase), rex.weight, 
                                    rex.deltas, rex.triangles, rex.w3js, rex.w6js, rex.w9js, rex.ylms, rex.djpqs )
        while true
            cont = false
            wa = RacahAlgebra.sumRulesForTwoW6jOneW9j(newrex);       if    wa[1]  newrex = wa[2];   cont = true  end
            wa = RacahAlgebra.sumRulesForOneW6jTwoW9j(newrex);       if    wa[1]  newrex = wa[2];   cont = true  end
            wa = RacahAlgebra.sumRulesForFourW6j(newrex);            if    wa[1]  newrex = wa[2];   cont = true  end
            #
            wa = RacahAlgebra.sumRulesForOneW3j(newrex);             if    wa[1]  newrex = wa[2];   cont = true  end
            wa = RacahAlgebra.sumRulesForOneW6j(newrex);             if    wa[1]  newrex = wa[2];   cont = true  end
            wa = RacahAlgebra.sumRulesForOneW9j(newrex);             if    wa[1]  newrex = wa[2];   cont = true  end
            wa = RacahAlgebra.sumRulesForTwoW3j(newrex);             if    wa[1]  newrex = wa[2];   cont = true  end
            wa = RacahAlgebra.sumRulesForTwoW6j(newrex);             if    wa[1]  newrex = wa[2];   cont = true  end
            wa = RacahAlgebra.sumRulesForOneW6jOneW9j(newrex);       if    wa[1]  newrex = wa[2];   cont = true  end
            wa = RacahAlgebra.sumRulesForTwoW9j(newrex);             if    wa[1]  newrex = wa[2];   cont = true  end
            wa = RacahAlgebra.sumRulesForThreeW3j(newrex);           if    wa[1]  newrex = wa[2];   cont = true  end
            wa = RacahAlgebra.sumRulesForTwoW3jOneW6j(newrex);       if    wa[1]  newrex = wa[2];   cont = true  end
            wa = RacahAlgebra.sumRulesForThreeW6j(newrex);           if    wa[1]  newrex = wa[2];   cont = true  end
            wa = RacahAlgebra.sumRulesForThreeW9j(newrex);           if    wa[1]  newrex = wa[2];   cont = true  end
            wa = RacahAlgebra.sumRulesForFourW3j(newrex);            if    wa[1]  newrex = wa[2];   cont = true  end
            wa = RacahAlgebra.sumRulesForThreeW6jOneW9j(newrex);     if    wa[1]  newrex = wa[2];   cont = true  end
            wa = RacahAlgebra.sumRulesForTwoW6jTwoW9j(newrex);       if    wa[1]  newrex = wa[2];   cont = true  end
            wa = RacahAlgebra.sumRulesForOneW6jThreeW9j(newrex);     if    wa[1]  newrex = wa[2];   cont = true  end
            wa = RacahAlgebra.sumRulesForFiveW3j(newrex);            if    wa[1]  newrex = wa[2];   cont = true  end
            wa = RacahAlgebra.sumRulesForSixW3j(newrex);             if    wa[1]  newrex = wa[2];   cont = true  end
            
            ## wa = RacahAlgebra.integralRulesForOneYlm(newrex);        if    wa[1]  newrex = wa[2];   cont = true  end
            ## wa = RacahAlgebra.integralRulesForTwoYlm(newrex);        if    wa[1]  newrex = wa[2];   cont = true  end
            
            newrex = RacahAlgebra.simplifyDeltas(newrex)
            newrex = RacahAlgebra.simplifyTriangles(newrex)
            
            if  cont    else     return(newrex)    end
        end
    end
    
end


"""
`RacahAlgebra.evaluate(leftCsq::RacahAlgebra.Csq, rightCsq::RacahAlgebra.Csq)`  
    ... evaluates a general recoupling coefficient that is defined by two coupling sequences. It first translates the coupling
        sequences into Racah expressions, i.e. into sums of Wigner 3-j symbols. A non-zero coefficient is obtained only if
        leftCsq.c == rightCsq.c, since re-coupling coefficients are diagonal in the 'last' angular momentum. 
        A newrex::RacahExpression is returned once a (single) simplification has been found, and nothing otherwise. 
        No attempt is presently made to find further simplications, once a rule has been applied.
"""
function  evaluate(leftCsq::RacahAlgebra.Csq, rightCsq::RacahAlgebra.Csq)
    function prsim(rx::RacahExpression)
        println(">> Simplification found for recoupling coefficient     < $leftCsq | $rightCsq > = $rx ")
    end
    #
    if  leftCsq.c != rightCsq.c     
        println("Recoupling coefficients vanishs because of different total angular momenta $(leftCsq.c) != $(rightCsq.c)")
        return( Basic(0) )
    end
    #
    # Rewrite coupling string as Racah expression and combine(multiply) both sides
    lrex = RacahAlgebra.rewriteCsq(leftCsq, "m");       println(">> Lhs-sequence     = $lrex ")
    rrex = RacahAlgebra.rewriteCsq(rightCsq, "m");      println(">> Rhs-sequence     = $rrex ")
    rex  = lrex * rrex
    summations = unique(rex.summations);    rex = RacahExpression(rex, summations=summations)
    println(">> Total recoupling = $rex")
    #
    wa = RacahAlgebra.evaluate(rex)
    return( wa )
end


"""
`RacahAlgebra.evaluateNumerical(w3j::W3j)`  
    ... attempts to evaluates  a Wigner 3-j symbol numerically; it is supposed that all angular momenta are given numerically.
        A wa::Float64 is returned.
"""
function  evaluateNumerical(w3j::W3j)
    ja = Float64(w3j.ja);   jb = Float64(w3j.jb);   jc = Float64(w3j.jc)
    ma = Float64(w3j.ma);   mb = Float64(w3j.mb);   mc = Float64(w3j.mc)
    wa = AngularMomentum.Wigner_3j(ja, jb, jc, ma, mb, mc)
    return( wa )
end


"""
`RacahAlgebra.evaluateNumerical(w6j::W6j)`  
    ... attempts to evaluates  a Wigner 6-j symbol numerically; it is supposed that all angular momenta are given numerically.
        A wa::Float64 is returned.
"""
function  evaluateNumerical(w6j::W6j)
    a = Float64(w6j.a);   b  = Float64(w6j.b);   c = Float64(w6j.c)
    d = Float64(w6j.e);   ex = Float64(w6j.e);   f = Float64(w6j.f)
    wa = AngularMomentum.Wigner_6j(a,b,c, d,ex,f)
    return( wa )
end


"""
`RacahAlgebra.evaluateNumerical(rex::RacahExpression)`  
    ... attempts to evaluates  a Racah expression numerically; it is supposed that the phase, weight and all angular momenta 
        are given numerically. A newRex::RacahExpression is returned.
"""
function  evaluateNumerical(rex::RacahExpression)
    wa = 1.0
    wa = wa * Float64((-1)^rex.phase) * Float64(rex.weight)
    for  w3j in rex.w3js    wa = wa * RacahAlgebra.evaluateNumerical(w3j)   end
    for  w6j in rex.w6js    wa = wa * RacahAlgebra.evaluateNumerical(w6j)   end
    for  w9j in rex.w9js    wa = wa * RacahAlgebra.evaluateNumerical(w9j)   end
    
    newRex = RacahExpression(rex.summations, Integral[], 0, wa, rex.deltas, rex.triangles, W3j[], W6j[], W9j[], Ylm[], Djpq[])
    return( newRex )
end


"""
`RacahAlgebra.hasAllVars(indexList::Array{SymEngine.Basic,1}, varList::Array{SymEngine.Basic,1})`  
    ... returns true if all indices from indexList are in varList and false otherwise.
"""
function  hasAllVars(indexList::Array{SymEngine.Basic,1}, varList::Array{SymEngine.Basic,1})
    for  index  in  indexList
        if  index in varList   ||  -index in varList    else    return( false )   end
    end
    return( true )
end


"""
`RacahAlgebra.hasNoVars(indexList::Array{SymEngine.Basic,1}, expr::SymEngine.Basic)`  
    ... returns true if no (symbolic) index from indexList occurs in expression and false otherwise.
"""
function hasNoVars(indexList::Array{SymEngine.Basic,1}, expr::SymEngine.Basic)
    sList = SymEngine.free_symbols(expr)
    for  index in indexList
        if  index in sList  || -index in sList   return( false )   end
    end
    
    return( true )    
end


"""
`RacahAlgebra.hasNoVars(indexList::Array{SymEngine.Basic,1}, deltas::Array{RacahAlgebra.Kronecker,1})`  
    ... returns true if no (symbolic) index from indexList occurs in the array deltas and false otherwise.
"""
function hasNoVars(indexList::Array{SymEngine.Basic,1}, deltas::Array{RacahAlgebra.Kronecker,1})
    sList = Basic[]
    for  delta in deltas    push!(sList, delta.i);    push!(sList, delta.k)     end
    for  index in indexList
        if  index in sList   return( false )   end
    end
    
    return( true )    
end


"""
`RacahAlgebra.hasNoVars(indexList::Array{SymEngine.Basic,1}, triangles::Array{RacahAlgebra.Triangle,1})`  
    ... returns true if no (symbolic) index from indexList occurs in the array triangles and false otherwise.
"""
function hasNoVars(indexList::Array{SymEngine.Basic,1}, triangles::Array{RacahAlgebra.Triangle,1})
    sList = Basic[]
    for  triangle in triangles    push!(sList, triangle.ja);    push!(sList, triangle.jb);    push!(sList, triangle.jc)     end
    for  index in indexList
        if  index in sList   return( false )   end
    end
    
    return( true )    
end


"""
`RacahAlgebra.hasNoVars(indexList::Array{SymEngine.Basic,1}, w3js::Array{RacahAlgebra.W3j,1})`  
    ... returns true if no (symbolic) index from indexList occurs in the array w3js and false otherwise.
"""
function hasNoVars(indexList::Array{SymEngine.Basic,1}, w3js::Array{RacahAlgebra.W3j,1})
    sList = Basic[]
    for  w3j in w3js    
        push!(sList, w3j.ja);    push!(sList, w3j.jb);    push!(sList, w3j.jc)
        push!(sList, w3j.ma);    push!(sList, w3j.mb);    push!(sList, w3j.mc)
    end
    for  index in indexList
        if  index in sList   return( false )   end
    end
    
    return( true )    
end


"""
`RacahAlgebra.hasNoVars(indexList::Array{SymEngine.Basic,1}, w6js::Array{RacahAlgebra.W6j,1})`  
    ... returns true if no (symbolic) index from indexList occurs in the array w6js and false otherwise.
"""
function hasNoVars(indexList::Array{SymEngine.Basic,1}, w6js::Array{RacahAlgebra.W6j,1})
    sList = Basic[]
    for  w6j in w6js    
        push!(sList, w6j.a);    push!(sList, w6j.b);    push!(sList, w6j.c)
        push!(sList, w6j.d);    push!(sList, w6j.e);    push!(sList, w6j.f)
    end
    for  index in indexList
        if  index in sList   return( false )   end
    end
    
    return( true )    
end


"""
`RacahAlgebra.hasNoVars(indexList::Array{SymEngine.Basic,1}, w9js::Array{RacahAlgebra.W9j,1})`  
    ... returns true if no (symbolic) index from indexList occurs in the array w3js and false otherwise.
"""
function hasNoVars(indexList::Array{SymEngine.Basic,1}, w9js::Array{RacahAlgebra.W9j,1})
    sList = Basic[]
    for  w9j in w9js    
        push!(sList, w9j.a);    push!(sList, w9j.b);    push!(sList, w9j.c)
        push!(sList, w9j.d);    push!(sList, w9j.e);    push!(sList, w9j.f)
        push!(sList, w9j.g);    push!(sList, w9j.h);    push!(sList, w9j.i)
    end
    for  index in indexList
        if  index in sList   return( false )   end
    end
    
    return( true )    
end


"""
`RacahAlgebra.hasNoVars(indexList::Array{SymEngine.Basic,1}, integrals::Array{RacahAlgebra.Integral,1})`  
    ... returns true if no (symbolic) index from indexList occurs in the array integrals and false otherwise.
"""
function hasNoVars(indexList::Array{SymEngine.Basic,1}, integrals::Array{RacahAlgebra.Integral,1})
    sList = Basic[]
    for  integral in integrals    
        push!(sList, integral.var);    push!(sList, integral.low);    push!(sList, integral.up)
    end
    for  index in indexList
        if  index in sList   return( false )   end
    end
    
    return( true )    
end


"""
`RacahAlgebra.hasNoVars(indexList::Array{SymEngine.Basic,1}, ylms::Array{RacahAlgebra.Ylm,1})`  
    ... returns true if no (symbolic) index from indexList occurs in the array ylms and false otherwise.
"""
function hasNoVars(indexList::Array{SymEngine.Basic,1}, ylms::Array{RacahAlgebra.Ylm,1})
    sList = Basic[]
    for  ylm in ylms    
        push!(sList, ylm.l);    push!(sList, ylm.m);    push!(sList, ylm.theta);    push!(sList, ylm.phi)
    end
    for  index in indexList
        if  index in sList   return( false )   end
    end
    
    return( true )    
end


"""
`RacahAlgebra.hasNoVars(indexList::Array{SymEngine.Basic,1}, djpqs::Array{RacahAlgebra.Djpq,1})`  
    ... returns true if no (symbolic) index from indexList occurs in the array djpqs and false otherwise.
"""
function hasNoVars(indexList::Array{SymEngine.Basic,1}, djpqs::Array{RacahAlgebra.Djpq,1})
    sList = Basic[]
    for  djpq in djpqs    
        push!(sList, djpq.j);    push!(sList, djpq.p);    push!(sList, djpq.q);    push!(sList, djpq.beta)
    end
    for  index in indexList
        if  index in sList   return( false )   end
    end
    
    return( true )    
end


"""
`RacahAlgebra.obsoletePhaseVars(rex::RacahExpression)`  
    ... determines those angular-momentum variables in the phase that should be removed by proper zero-phase terms.
        A list of vars::Array{Basic,1} is returned.
"""
function obsoletePhaseVars(rex::RacahExpression)
    vars  = Basic[]
    sList = SymEngine.free_symbols(rex.phase)
    for  s in sList
        if  !(s in SymEngine.free_symbols(rex.summations))                &&
            !(s in SymEngine.free_symbols(rex.weight))                    &&
            RacahAlgebra.hasNoVars([s], rex.integrals)                    &&
            RacahAlgebra.hasNoVars([s], rex.deltas)                       &&
            RacahAlgebra.hasNoVars([s], rex.triangles)                    &&
            RacahAlgebra.hasNoVars([s], rex.w3js)                         &&
            RacahAlgebra.hasNoVars([s], rex.w6js)                         &&
            RacahAlgebra.hasNoVars([s], rex.w9js)                         &&
            RacahAlgebra.hasNoVars([s], rex.ylms)                         &&
            RacahAlgebra.hasNoVars([s], rex.djpqs)                        
            push!(vars, s)
        end
    end
    
    return( vars )    
end


"""
`RacahAlgebra.purifyPhase(phase::SymEngine.Basic)`  
    ... purifies the phase so that all contributions of 'a' are in the range -a, a, 2a ,3a
        An equivalent newphase::SymEngine.Basic is returned.
"""
function purifyPhase(phase::SymEngine.Basic)
    newphase = Basic(0)
    waList   = SymEngine.get_args(phase + 4*Basic(:Xdummy) + 4*Basic(:Ydummy))  # Add Xdummy, Ydummy to have always more than 1 argument
    for  wa in waList   
        if       length(SymEngine.get_args(wa)) == 0    newphase = newphase + wa
        elseif   length(SymEngine.get_args(wa)) == 2
            wb = SymEngine.get_args(wa);   wc = convert(Int64,wb[1]);   wc = rem(wc+100,4);    if  wc == 3   wc = -1   end                                    
                                                        newphase = newphase + wc * wb[2]
        else     error("stop a")
        end
    end
    
    return( newphase )    
end


"""
`RacahAlgebra.recursionW3j(ww::W3j, rule::AbstractRecursionW3j)`  
    ... applies a given recursion rules to w3j; an rexList::Array{RacahAlgebra.RacahExpression,1} is returned.
"""
function recursionW3j(w3j::W3j, rule::AbstractRecursionW3j)
    #
    rexList = RacahExpression[]
    if      rule == RecursionW3jHalfStep()
        #
        # Half-step rule:   Rotenberg et al. (1959), Eq. (1.45).
        # ---------------
        #                  1/2  ( j1  j2  j3 )                    1/2  ( j1  j2-1/2  j3-1/2 )                    1/2  ( j1  j2-1/2  j3-1/2 )
        #   [(J+1)(J-2*j1)]     (            )  = [(j2+m2)(j3-m3)]     (                    )  - [(j2-m2)(j3+m3)]     (                    )
        #                       ( m1  m2  m3 )                         ( m1  m2-1/2  m3+1/2 )                         ( m1  m2+1/2  m3-1/2 )
        #
        #   with  J = j1+j2+j3
        #
        J = w3j.ja + w3j.jb + w3j.jc;    wa = sqrt( (J+1)*(J-2*w3j.ja) )
        rex = RacahExpression( Basic[], Integral[], 0, sqrt( (w3j.jb+w3j.mb)*(w3j.jc-w3j.mc) )/wa, Kronecker[], Triangle[], 
                                [W3j(w3j.ja, w3j.jb-1//2, w3j.jc-1//2, w3j.ma, w3j.mb-1//2, w3j.mc+1//2)], W6j[], W9j[], Ylm[], Djpq[])
        push!( rexList, rex)
        rex = RacahExpression( Basic[], Integral[], 0, -sqrt( (w3j.jb-w3j.mb)*(w3j.jc+w3j.mc) )/wa, Kronecker[], Triangle[], 
                                [W3j(w3j.ja, w3j.jb-1//2, w3j.jc-1//2, w3j.ma, w3j.mb+1//2, w3j.mc-1//2)], W6j[], W9j[], Ylm[], Djpq[])
        push!( rexList, rex)
        return( rexList )
        #
    elseif  rule == RecursionW3jLouck()
        #
        # Louck's rule:   Rotenberg et al. (1959), Eq. (1.46).
        # -------------
        #          1/2          (  j1     j2   j3 )                                1/2  (   j1    j2-1/2  j3-1/2 )
        #   (j2+m2)   (2*j3+1)  (                 ) == - [(J-2*j1)(J+1)(j3+m3)]     (                        )
        #                       ( m2-m3  -m2   m3 )                                     ( m2-m3  -m2+1/2  m3-1/2 )
        #                                             
        #                                                                                  1/2  (   j1    j2-1/2  j3+1/2 )
        #                                              - [(J-2*j3)(J-2*j2+1)(J+1)(j3-m3+1)]     (                        )
        #                                                                                       ( m2-m3  -m2+1/2  m3-1/2 )
        #
        #  with  J = j1+j2+j3
        #
        J = w3j.ja + w3j.jb + w3j.jc;    wa = sqrt( (w3j.jb+w3j.mb) ) * (2*w3j.jc+1)
        rex = RacahExpression( Basic[], Integral[], 0, -sqrt( (J-2*w3j.ja)* (J+1) * (w3j.jc+w3j.mc) )/wa, Kronecker[], Triangle[], 
                                [W3j(w3j.ja, w3j.jb-1//2, w3j.jc-1//2, w3j.mb-w3j.mc, -w3j.mb+1//2, w3j.mc-1//2)], W6j[], W9j[], Ylm[], Djpq[])
        push!( rexList, rex)
        rex = RacahExpression( Basic[], Integral[], 0, -sqrt( (J-2*w3j.jc) * (J-2*w3j.jb+1) * (J+1) * (w3j.jc-w3j.mc+1) )/wa, Kronecker[], Triangle[], 
                                [W3j(w3j.ja, w3j.jb-1//2, w3j.jc+1//2, w3j.mb-w3j.mc, -w3j.mb+1//2, w3j.mc-1//2)], W6j[], W9j[], Ylm[], Djpq[])
        push!( rexList, rex)
        return( rexList )
        #
    elseif  rule == RecursionW3jOneStep()
        #
        # One-step rule:   Rotenberg et al. (1959), Eq. (1.47).
        # --------------
        #                                    1/2  ( j1  j2  j3 )                                      1/2  ( j1   j2   j3-1 )
        #   [(J+1)(J-2*j1)(J-2*j2)(J-2*j3+1)]     (            ) == [(j2-m2)(j2+m2+1)(j3+m3)(j3+m3-1)]     (                )
        #                                         ( m1  m2  m3 )                                           ( m1  m2+1  m3-1 )
        #
        #                                         1/2  ( j1  j2  j3-1 )                                     1/2  ( j1   j2   j3-1 )
        #                  - 2*m2 [(j3+m3)(j3-m3)]     (              ) - [(j2+m2)(j2-m2+1)(j3-m3)(j3-m3-1)]     (                )
        #                                              ( m1  m2   m3  )                                          ( m1  m2-1  m3+1 )
        #
        #    with  J = j1+j2+j3
        #
        J = w3j.ja + w3j.jb + w3j.jc;    wa = sqrt( (J+1) * (J-2*w3j.ja) * (J-2*w3j.jb) * (J-2*w3j.jc+1) )
        rex = RacahExpression( Basic[], Integral[], 0, sqrt( (w3j.jb-w3j.mb) * (w3j.jb+w3j.mb+1) * (w3j.jc+w3j.mc) * (w3j.jc+w3j.mc-1))/wa, 
                                Kronecker[], Triangle[], [W3j(w3j.ja, w3j.jb, w3j.jc-1, w3j.ma, w3j.mb+1, w3j.mc-1)], W6j[], W9j[], Ylm[], Djpq[])
        push!( rexList, rex)
        rex = RacahExpression( Basic[], Integral[], 0, 2*w3j.jb * sqrt( (w3j.jc+w3j.mc) * (w3j.jc-w3j.mc))/wa, Kronecker[], Triangle[], 
                                [W3j(w3j.ja, w3j.jb, w3j.jc-1, w3j.ma, w3j.mb, w3j.mc)], W6j[], W9j[], Ylm[], Djpq[])
        push!( rexList, rex)
        rex = RacahExpression( Basic[], Integral[], 0, - sqrt( (w3j.jb+w3j.mb) * (w3j.jb-w3j.mb+1) * (w3j.jc-w3j.mc) * (w3j.jc-w3j.mc-1))/wa, 
                                Kronecker[], Triangle[], [W3j(w3j.ja, w3j.jb, w3j.jc-1, w3j.ma, w3j.mb-1, w3j.mc+1)], W6j[], W9j[], Ylm[], Djpq[])
        push!( rexList, rex)
        return( rexList )
        #
    elseif  rule == RecursionW3jMagnetic()
        #
        # Magnetic rule:   Rotenberg et al. (1959), Eq. (1.47).
        # --------------
        #                            1/2   (  j1  j2    j3  )                     1/2   (  j1   j2   j3  )                     1/2   (  j1   j2   j3  )
        #   -[(j3+m1+m2+1)(j3-m1-m2)]      (                ) = [(j1+m1+1)(j1-m1)]      (                ) + [(j2+m2+1)(j2-m2)]      (                )
        #                                  (  m1  m2  -m3+1 )                           (  m1  m2+1  -m3 )                           ( m1+1  m2  -m3  )
        #
        #    with  J = j1+j2+j3
        #
        J = w3j.ja + w3j.jb + w3j.jc;    wa = - sqrt( (w3j.jc+w3j.ma+w3j.mb+1) * (w3j.jc-w3j.ma-w3j.mb) )
        rex = RacahExpression( Basic[], Integral[], 0, sqrt( (w3j.ja+w3j.ma+1) * (w3j.ja-w3j.ma) )/wa, Kronecker[], Triangle[], 
                                [W3j(w3j.ja, w3j.jb, w3j.jc, w3j.ma, w3j.mb+1, w3j.mc-1)], W6j[], W9j[], Ylm[], Djpq[])
        push!( rexList, rex)
        rex = RacahExpression( Basic[], Integral[], 0, sqrt( (w3j.jb+w3j.mb+1) * (w3j.jb-w3j.mb) )/wa, Kronecker[], Triangle[], 
                                [W3j(w3j.ja, w3j.jb, w3j.jc, w3j.ma+1, w3j.mb, w3j.mc-1)], W6j[], W9j[], Ylm[], Djpq[])
        push!( rexList, rex)
        return( rexList )
        #
    else    error("Unrecognized recursion rule of Wigner 3j symbols.")
    end
end


"""
`RacahAlgebra.removeIndex(index::SymEngine.Basic, indexList::Array{SymEngine.Basic,1})`  
    ... removes the index from the given indexList; a newList::Array{SymEngine.Basic,1} is returned.
"""
function  removeIndex(index::SymEngine.Basic, indexList::Array{SymEngine.Basic,1})
    newList = Basic[]
    for  idx in indexList
        if  idx == index  ||  -idx in indices
        else    push!( newList, idx)
        end
    end
    return( newList )
end


"""
`RacahAlgebra.removeIndex(indices::Array{SymEngine.Basic,1}, indexList::Array{SymEngine.Basic,1})`  
    ... removes the indices from the given indexList; a newList::Array{SymEngine.Basic,1} is returned.
"""
function  removeIndex(indices::Array{SymEngine.Basic,1}, indexList::Array{SymEngine.Basic,1})
    newList = Basic[]
    for  idx in indexList
        if  idx in indices  ||  -idx in indices
        else    push!( newList, idx)
        end
    end
    return( newList )
end


"""
`RacahAlgebra.removeIntegral(integral::Integral, integralList::Array{Integral,1})`  
    ... removes the index from the given indexList; a newList::Array{Integral,1} is returned.
"""
function  removeIntegral(integral::Integral, integralList::Array{Integral,1})
    newList = Basic[]
    for  ntg in integralList
        if  ntg == integral
        else    push!( newList, ntg)
        end
    end
    return( newList )
end


"""
`RacahAlgebra.rewriteCsq(csq::RacahAlgebra.Csq, ms::String)`  
    ... attempts to rewrite the coupling sequence csq as Racah expressions, i.e. as (sum of) products of Wigner 3-j symbols. 
        A rex::RacahExpression is returned.
"""
function rewriteCsq(csq::RacahAlgebra.Csq, ms::String)
    function  appendSummation(a::RacahAlgebra.Csq, b::RacahAlgebra.Csq)
        sx = Basic[];   sx = append!( sx, appendSummation(a.a, a.b));   sx = append!( sx, appendSummation(b.a, b.b))
        sx = push!( sx, Basic(ms * string(a.c)));                       sx = push!( sx, Basic(ms * string(b.c))); 
        return( sx )
    end
    function  appendSummation(a::RacahAlgebra.Csq, b::Basic)
        sx = Basic[];   sx = append!( sx, appendSummation(a.a, a.b));   sx = push!( sx, Basic(ms * string(b)))
        sx = push!( sx, Basic(ms * string(a.c)));
        return( sx )
    end
    function  appendSummation(a::Basic, b::RacahAlgebra.Csq)
        sx = Basic[];   sx = push!( sx, Basic(ms * string(a)));         sx = append!( sx, appendSummation(b.a, b.b))
        sx = push!( sx, Basic(ms * string(b.c)));
        return( sx )
    end
    function  appendSummation(a::Basic, b::Basic)
        sx = Basic[];   sx = push!( sx, Basic(ms * string(a)));         sx = push!( sx, Basic(ms * string(b)))
        return( sx )
    end
    #
    function  appendW3j(a::RacahAlgebra.Csq, b::RacahAlgebra.Csq, c::Basic)
        wx = W3j[];   wx = append!( wx, appendW3j(a.a, a.b, a.c));   wx = append!( wx, appendW3j(b.a, b.b, b.c))
        wx = append!( wx, [W3j(a.c,b.c,c, Basic(ms * string(a.c)),  Basic(ms * string(b.c)),  -Basic(ms * string(c))) ] )
        return( wx )
    end
    function  appendW3j(a::RacahAlgebra.Csq, b::Basic, c::Basic)
        wx = W3j[];   wx = append!( wx, appendW3j(a.a, a.b, a.c))
        wx = append!( wx, [W3j(a.c,b,c, Basic(ms * string(a.c)),  Basic(ms * string(b)),  -Basic(ms * string(c))) ] )
        return( wx )
    end
    function  appendW3j(a::Basic, b::RacahAlgebra.Csq, c::Basic)
        wx = W3j[];   wx = append!( wx, appendW3j(b.a, b.b, b.c))
        wx = append!( wx, [W3j(a,b.c,c, Basic(ms * string(a)),  Basic(ms * string(b.c)),  -Basic(ms * string(c))) ] )
        return( wx )
    end
    function  appendW3j(a::Basic, b::Basic, c::Basic)
        wx = W3j[];   wx = append!( wx, [W3j(a,b,c, Basic(ms * string(a)),  Basic(ms * string(b)),  -Basic(ms * string(c)) )] )
        return( wx )
    end
    #
    summations = Basic[];   phase = Basic(0);   weight = Basic(1);   w3js = W3j[]
    summations = appendSummation( csq.a, csq.b)
    w3js       = appendW3j( csq.a, csq.b, csq.c)
    for  w3j in w3js
        phase  = phase  + w3j.ja - w3j.jb - w3j.mc
        weight = weight * sqrt(2*w3j.jc+1)
    end
    
    rex = RacahExpression( summations, Integral[], phase, weight, Kronecker[], Triangle[], w3js, W6j[], W9j[], Ylm[], Djpq[] )    
    ## println("***** rex = $rex")
    return( rex )
end


"""
`RacahAlgebra.rewritePhase(phase::Basic, zeroTerms::Array{SymEngine.Basic,1}, woIndex::Array{SymEngine.Basic,1}; 
                            printout::Bool=false, from::String="Unspecified source")`  
    ... attempts to rewrite the phase by adding one or several 'zero' terms so that it appears without the indices 
        in woIndex. An equivalent newPhase::Basic either 'without' or 'with' the indicated indices is returned.
"""
function rewritePhase(phase::Basic, zeroTerms::Array{SymEngine.Basic,1}, woIndex::Array{SymEngine.Basic,1}; 
                        printout::Bool=false, from::String="Unspecified source")
    # Return the phase if this is OK
    if  RacahAlgebra.hasNoVars(woIndex, phase)    return( phase )   end
    if  printout    println("$from  phase = $phase   without = $woIndex   \n    zeroTerms = $zeroTerms")   end
    # First simply try to add/substract multiples of each zeroTerm and see whether this solves the issue
    purePhase = RacahAlgebra.purifyPhase(phase)
    for  z1Term  in  zeroTerms
        for  k1 in [-2, -1, 0, 1, 2]  
            for  z2Term  in  zeroTerms
                for  k2 in [-2, -1, 1, 2]  
                    newPhase = SymEngine.expand( purePhase + k1 * z1Term + k2 * z2Term )
                    newPhase = RacahAlgebra.purifyPhase(newPhase)
                    if  RacahAlgebra.hasNoVars(woIndex, newPhase)   
                        if  printout    println("$from  newPhase = $newPhase  done !!!!")   end
                        return( RacahAlgebra.purifyPhase(newPhase) )   
                    end
                end
            end
        end
    end
    return( phase )
end


"""
`RacahAlgebra.rewritePhase(rex::RacahExpression, zeroTerms::Array{SymEngine.Basic,1}; printout::Bool=false, from::String="Unspecified source")`  
    ... attempts to rewrite (and shorten) the phase by adding one or several 'zero' terms so that 
        (i) no obsolete indices appear;
        (ii) the number of phase indices becomes smaller
        A modified Racah expression is returned if the re-writing was successful, and nothing is changed otherwise.
"""
function rewritePhase(rex::RacahExpression, zeroTerms::Array{SymEngine.Basic,1}; printout::Bool=false, from::String="Unspecified source")
    #
    function rewritePhaseObsolete(phase::Basic)
        if  vars != Basic[]
            purePhase = RacahAlgebra.purifyPhase(phase)
            for  z1Term  in  zeroTerms
                # Only consider potential zeros with m quantum numbers
                if  length(SymEngine.free_symbols(z1Term)) > 2     continue     end
                for  k1 in [-2, -1, 0, 1, 2]  
                    for  z2Term  in  zeroTerms
                        if  length(SymEngine.free_symbols(z2Term)) > 2     continue     end
                        for  k2 in [-2, -1, 1, 2]  
                            newPhase = SymEngine.expand( purePhase + k1 * z1Term + k2 * z2Term )
                            if  RacahAlgebra.hasNoVars(vars, newPhase)   
                                if  printout    println("$from  newPhase = $newPhase  done !!!!")   end
                                return( true, RacahAlgebra.purifyPhase(newPhase) )   
                            end
                        end
                    end
                end
            end
        end
        return( false, phase )
    end
    #
    function rewritePhaseShorter(phase::Basic)
        Nofs = length(SymEngine.free_symbols(phase))
        if  Nofs > 0
            for  z1Term  in  zeroTerms
                for  k1 in [-2, -1, 0, 1, 2]  
                    for  z2Term  in  zeroTerms
                        for  k2 in [-2, -1, 1, 2]  
                            newPhase = SymEngine.expand( phase + k1 * z1Term + k2 * z2Term )
                            if  length(SymEngine.free_symbols(newPhase)) < Nofs  
                                if  printout    println("$from  newPhase = $newPhase  now shorter !!!!")   end
                                return( true, RacahAlgebra.purifyPhase(newPhase) )   
                            end
                        end
                    end
                end
            end
        end
        return( false, phase )
    end
    #
    hasChanged1 = false;    hasChanged2 = false
    #
    # (i)  Try to remove obsolete indices; simply try to add/substract multiples of each zeroTerm and see whether this solves the issue
    # (ii) Try to shorten/reduce the number of phase variables
    vars = RacahAlgebra.obsoletePhaseVars(rex)
    if  printout    println("$from:::  phase = $(rex.phase)   without = $vars   \n    zeroTerms = $zeroTerms")   end
    hasChanged1, phase = rewritePhaseObsolete(rex.phase)
    hasChanged2, phase = rewritePhaseShorter(phase)
    if  hasChanged1   ||    hasChanged2
        rey = RacahExpression(rex.summations, rex.integrals, phase, rex.weight, rex.deltas, rex.triangles, rex.w3js, rex.w6js, rex.w9js, rex.ylms, rex.djpqs)
    else 
        rey = rex
    end
    #
    return( rey )
end


"""
`RacahAlgebra.selectW3j(n::Int64)`  
    ... selects one of various pre-defined Wigner 3j symbols for which usually special values are known;
        this function has been implemented mainly for test purposes. A w3j::W3j is returned. 
        If n = 99, all pre-defined Wigner 3j symbols are printed to screen and nothing is returned in this
        case.
"""
function selectW3j(n::Int64)
    if  n == 99     for  i = 1:20   println("  $i    $(selectW3j(i))")      end
        return(nothing)
    end
    
    j1 = Basic(:j1);    j2 = Basic(:j2);    j3 = Basic(:j3)
    j  = Basic(:j);     m = Basic(:m);      
    
    if      n ==  1     w3j = W3j(j1, j2, j3, 0, 0, 0)
    elseif  n ==  2     w3j = W3j(j, j, 0, m, -m, 0)
    elseif  n ==  3     w3j = W3j(j, j-1//2, 1//2, m, -m-1//2, 1//2)
    elseif  n ==  4     w3j = W3j(j+1, j, 1, m, -m-1, 1)
    elseif  n ==  5     w3j = W3j(j+1, j, 1, m, -m, 0)
    elseif  n ==  6     w3j = W3j(j, j, 1, m, -m-1, 1)
    elseif  n ==  7     w3j = W3j(j, j, 1, m, -m, 0)
    elseif  n ==  8     w3j = W3j(j+3//2, j, 3//2, m, -m-3//2, 3//2)
    elseif  n ==  9     w3j = W3j(j+3//2, j, 3//2, m, -m-1//2, 1//2)
    elseif  n == 10     w3j = W3j(j+1//2, j, 3//2, m, -m-1//2, 1//2)
    elseif  n == 11     w3j = W3j(j+2, j, 2, m, -m-2, 2)
    elseif  n == 12     w3j = W3j(j+2, j, 2, m, -m, 0)
    elseif  n == 13     w3j = W3j(j+1, j, 2, m, -m-2, 2)
    elseif  n == 14     w3j = W3j(j+1, j, 2, m, -m-1, 1)
    elseif  n == 15     w3j = W3j(j+1, j, 2, m, -m, 0)
    elseif  n == 16     w3j = W3j(j, j, 2, m, -m, 0)

    elseif  n == 20     w3j = FAIL
    else    error("stop a")
    end
    
    return( w3j )
end


"""
`RacahAlgebra.selectW6j(n::Int64)`  
    ... selects one of various pre-defined Wigner 6j symbols for which usually special values are known; this function has been 
        implemented mainly for test purposes. A w6j::W6j is returned. If n = 99, all pre-defined Wigner 6j symbols are printed to 
        screen and nothing is returned in this case.
"""
function selectW6j(n::Int64)
    if  n == 99     for  i = 1:20   println("  $i    $(selectW6j(i))")      end
        return(nothing)
    end
    
    j1 = Basic(:j1);    j2 = Basic(:j2);    j3 = Basic(:j3);    l1 = Basic(:l1);   l2 = Basic(:l2)
    j  = Basic(:j);     m = Basic(:m);      a  = Basic(:a);     b  = Basic(:b);    c  = Basic(:c)
    
    if      n ==  1     w6j = W6j( j1, j2, j3, l1, l2, 0)
    elseif  n ==  2     w6j = W6j( j1, j2, j3, 1//2, j3-1//2, j2+1//2)
    elseif  n ==  3     w6j = W6j( j1, j2, j3, 1//2, j3+1//2, j2+1//2)
    elseif  n ==  4     w6j = W6j( j1, j2, j3, 1, j3-1, j2-1)
    elseif  n ==  5     w6j = W6j( j1, j2, j3, 1, j3-1, j2)
    elseif  n ==  6     w6j = W6j( j1, j2, j3, 1, j3-1, j2+1)
    elseif  n ==  7     w6j = W6j( j1, j2, j3, 1, j3, j2)
    elseif  n ==  8     w6j = W6j( a, b, c, 3//2, c-3//2, b-3//2)
    elseif  n ==  9     w6j = W6j( a, b, c, 3//2, c-3//2, b-1//2)
    elseif  n == 10     w6j = W6j( a, b, c, 3//2, c-3//2, b+1//2)
    elseif  n == 11     w6j = W6j( a, b, c, 3//2, c-3//2, b+3//2)
    elseif  n == 12     w6j = W6j( a, b, c, 3//2, c-1//2, b-1//2)
    elseif  n == 13     w6j = W6j( a, b, c, 3//2, c-1//2, b+1//2)
    elseif  n == 14     w6j = W6j( a, b, c, 2, c-2, b-2)
    elseif  n == 15     w6j = W6j( a, b, c, 2, c-2, b-1)
    elseif  n == 16     w6j = W6j( a, b, c, 2, c-2, b)
    elseif  n == 17     w6j = W6j( a, b, c, 2, c-2, b+1)
    elseif  n == 18     w6j = W6j( a, b, c, 2, c-2, b+2)
    elseif  n == 19     w6j = W6j( a, b, c, 2, c-1, b-1)
    elseif  n == 20     w6j = W6j( a, b, c, 2, c-1, b)
    elseif  n == 21     w6j = W6j( a, b, c, 2, c-1, b+1)
    elseif  n == 22     w6j = W6j( a, b, c, 2, c  , b)
    else    error("stop a")
    end
    
    return( w6j )
end


"""
`RacahAlgebra.selectW9j(n::Int64)`  
    ... selects one of various pre-defined Wigner 9j symbols for which usually special values are known;
        this function has been implemented mainly for test purposes. A w9j::W9j is returned. 
        If n = 99, all pre-defined Wigner 9j symbols are printed to screen and nothing is returned in this
        case.
"""
function selectW9j(n::Int64)
    if  n == 99     for  i = 1:1   println("  $i    $(selectW9j(i))")      end
        return(nothing)
    end
    
    a = Basic(:a);    b = Basic(:b);    c = Basic(:c);    d = Basic(:d);    ee = Basic(:ee)
    f = Basic(:f);    g = Basic(:g);    h = Basic(:h)
    
    if      n ==  1     w9j = W9j(a, b, 0, c, d, 0, ee, f, 0)
    elseif  n ==  2     w9j = W9j(a, b, c, d, ee, f, g, h, 0)
    else    error("stop a")
    end
    
    return( w9j )
end


"""
`RacahAlgebra.selectRacahExpression(n::Int64)`  
    ... selects one of various pre-defined Racah expression as they often occur on the lhs of some sum rule;
        this function has been implemented mainly for test purposes. A rex::RacahExpression is returned. 
        If n = 99, all pre-defined RacahExpression are printed to screen and nothing is returned in this
        case.
"""
function selectRacahExpression(n::Int64)
    if  n == 99     for  i = 1:20   println("  $i    $(selectRacahExpression(i))")      end
        return(nothing)
    end
    
    j = Basic(:j);    J = Basic(:J);    m = Basic(:m);    M  = Basic(:M)
    a = Basic(:a);    b = Basic(:b);    c = Basic(:c);    d  = Basic(:d);    ee  = Basic(:ee);    f  = Basic(:f) 
    g = Basic(:g);    h = Basic(:h);    k = Basic(:k);    l  = Basic(:l);    p   = Basic(:p);     q  = Basic(:q);     
    r = Basic(:r);    s = Basic(:s);    t = Basic(:t);    X = Basic(:X);     Y   = Basic(:Y);     Z  = Basic(:Z);
        
        
    na  = Basic(:na);    np  = Basic(:np);     nq  = Basic(:nq);     ee  = Basic(:ee);    ap  = Basic(:ap) 
    jp = Basic(:jp);     bp  = Basic(:bp);     cp  = Basic(:cp);     dp  = Basic(:dp);    eep = Basic(:eep);   fp  = Basic(:fp);
    j1  = Basic(:j1);    j2  = Basic(:j2);     j3  = Basic(:j3);     m1  = Basic(:m1);    m2  = Basic(:m2);    m3  = Basic(:m3)
    j4  = Basic(:j4);    j5  = Basic(:j5);     j6  = Basic(:j6);     m4  = Basic(:m4);    m5  = Basic(:m5);    m6  = Basic(:m6)
    j7  = Basic(:j7);    j8  = Basic(:j8);     j9  = Basic(:j9);     m7  = Basic(:m7);    m8  = Basic(:m8);    m9  = Basic(:m9)
    j10 = Basic(:j10);   j11 = Basic(:j11);    j12 = Basic(:j12);    m10 = Basic(:m10);   m11 = Basic(:m11);   m12 = Basic(:m12)
    l1  = Basic(:l1);    l2  = Basic(:l2);     l3  = Basic(:l3);     n1  = Basic(:n1);    n2  = Basic(:n2);    n3  = Basic(:n3)
    m1p = Basic(:m1p);   m2p = Basic(:m2p);    j3p = Basic(:j3p);    m3p = Basic(:m3p);   nap = Basic(:nap)
    
    # Sum rules for one Wnj symbol
    if      n ==  1     w3j = W3j(j, j, J, m, -m, M)
                        rex = RacahExpression( [m], Integral[], Basic(-m), Basic(1), Kronecker[], Triangle[], [w3j], W6j[], W9j[], Ylm[], Djpq[] )
    elseif  n ==  2     w6j = W6j(a, b, X, a, b, c)
                        rex = RacahExpression( [X], Integral[], Basic(0), Basic(2*X+1), Kronecker[], Triangle[], W3j[], [w6j], W9j[], Ylm[], Djpq[] )
    elseif  n ==  3     w6j = W6j(a, b, X, b, a, c)
                        rex = RacahExpression( [X], Integral[], Basic(X), Basic(2*X+1), Kronecker[], Triangle[], W3j[], [w6j], W9j[], Ylm[], Djpq[] )
    elseif  n ==  4     w9j = W9j(a, b, ee, c, d, f, ee, f, X)
                        rex = RacahExpression( [X], Integral[], Basic(0), Basic(2*X+1), Kronecker[], Triangle[], W3j[], W6j[], [w9j], Ylm[], Djpq[] )
    elseif  n ==  5     w9j = W9j(a, b, ee, c, d, f, f, ee, X)
                        rex = RacahExpression( [X], Integral[], Basic(-X), Basic(2*X+1), Kronecker[], Triangle[], W3j[], W6j[], [w9j], Ylm[], Djpq[] )
    # Sum rules for two Wnj symbol
    elseif  n ==  6     aw3j = W3j(j1, j2, j3, m1, m2, m3);    bw3j = W3j(j1, j2, j3, m1p, m2p, m3)
                        rex  = RacahExpression( [j3, m3], Integral[], Basic(0), Basic(2*j3+1), Kronecker[], Triangle[], [aw3j, bw3j], W6j[], W9j[], Ylm[], Djpq[] )
    elseif  n ==  7     aw3j = W3j(j1, j2, j3, m1, m2, m3);    bw3j = W3j(j1, j2, j3p, m1, m2, m3p)
                        rex  = RacahExpression( [m1, m2], Integral[], Basic(0), Basic(1), Kronecker[], Triangle[], [aw3j, bw3j], W6j[], W9j[], Ylm[], Djpq[] )
    elseif  n ==  8     aw3j = W3j(a, p, q, -na, np, nq);    bw3j = W3j(p, q, ap, -np, -nq, nap)
                        rex  = RacahExpression( [np, nq], Integral[], Basic(-np-nq), Basic(1), Kronecker[], Triangle[], [aw3j, bw3j], W6j[], W9j[], Ylm[], Djpq[] )
    elseif  n ==  9     aw6j = W6j(X, Y, Z, a, b ,c);    bw6j = W6j(X, Y, Z, a, b ,c)
                        rex  = RacahExpression( [X, Y, Z], Integral[], Basic(0), Basic((2*X+1)*(2*Y+1)*(2*Z+1)), 
                                                Kronecker[], Triangle[], W3j[], W6j[aw6j, bw6j], W9j[], Ylm[], Djpq[] )
    elseif  n ==  10    aw6j = W6j(a, b, X, c, d, p);    bw6j = W6j(c, d, X, b, a, q)
                        rex  = RacahExpression( [X], Integral[], Basic(X), Basic(2*X+1), Kronecker[], Triangle[], W3j[], W6j[aw6j, bw6j], W9j[], Ylm[], Djpq[] )
    elseif  n ==  11    aw6j = W6j(a, b, X, c, d, p);    bw6j = W6j(c, d, X, a, b, q)
                        rex  = RacahExpression( [X], Integral[], Basic(0), Basic(2*X+1), Kronecker[], Triangle[], W3j[], W6j[aw6j, bw6j], W9j[], Ylm[], Djpq[] )
    elseif  n ==  12    aw6j = W6j(X, Y, Z, c, a, b);    bw9j = W9j(X, Y, Z, a, b, c, b, c, a)
                        rex  = RacahExpression( [X,Y,Z], Integral[], Basic(0), Basic((2*X+1)*(2*Y+1)*(2*Z+1)), Kronecker[], Triangle[], W3j[], 
                                                W6j[aw6j], W9j[bw9j], Ylm[], Djpq[] )
    elseif  n ==  13    aw6j = W6j(a, f, X, ee, b, s);   bw9j = W9j(a, f, X, d, q, ee, p, c, b)
                        rex  = RacahExpression( [X], Integral[], Basic(0), Basic(2*X+1), Kronecker[], Triangle[], W3j[], W6j[aw6j], W9j[bw9j], Ylm[], Djpq[] )
    elseif  n ==  14    aw6j = W6j(a, f, X, b, ee, s);   bw9j = W9j(a, f, X, d, q, ee, p, c, b)
                        rex  = RacahExpression( [X], Integral[], X, Basic(2*X+1), Kronecker[], Triangle[], W3j[], W6j[aw6j], W9j[bw9j], Ylm[], Djpq[] )
    elseif  n ==  15    aw9j = W9j(X, Y, Z, a, b, c, d, ee, f);   bw9j = W9j(X, Y, Z, a, b, c, d, ee, f)
                        rex  = RacahExpression( [X,Y,Z], Integral[], Basic(0), Basic((2*X+1)*(2*Y+1)*(2*Z+1)), Kronecker[], Triangle[], W3j[], 
                                                    W6j[], W9j[aw9j, bw9j], Ylm[], Djpq[] )
    elseif  n ==  16    aw9j = W9j(a, b, X, c, d, Y, ee, f, j);   bw9j = W9j(a, b, X, c, d, Y, g, h, j)
                        rex  = RacahExpression( [X,Y], Integral[], Basic(0), Basic((2*X+1)*(2*Y+1)), Kronecker[], Triangle[], W3j[], W6j[], W9j[aw9j, bw9j], Ylm[], Djpq[] )
    elseif  n ==  17    aw9j = W9j(a, b, X, c, d, Y, ee, f, j);   bw9j = W9j(a, b, X, d, c, Y, g, h, j)
                        rex  = RacahExpression( [X,Y], Integral[], Y, (2*X+1)*(2*Y+1), Kronecker[], Triangle[], W3j[], W6j[], W9j[aw9j, bw9j], Ylm[], Djpq[] )
    # Sum rules for three Wnj symbol
    elseif  n ==  18    aw3j = W3j(j5, j1, j6, m5, m1, -m6);    bw3j = W3j(j6, j2, j4, m6, m2, -m4);    cw3j = W3j(j4, j3, j5, m4, m3, -m5)
                        rex  = RacahExpression( [m4, m5, m6], Integral[], -m4 - m5 - m6, Basic(1), Kronecker[], Triangle[], W3j[aw3j, bw3j, cw3j], 
                                                W6j[], W9j[], Ylm[], Djpq[] )
    elseif  n ==  19    aw3j = W3j(l1, j2, l3, n1, m2, n3);    bw3j = W3j(j1, l2, l3, m1, n2, -n3);    cw6j = W6j(j1, j2, j3, l1, l2, l3)
                        rex  = RacahExpression( [l3, n3], Integral[], l3, 2*l3+1, Kronecker[], Triangle[], W3j[aw3j, bw3j], W6j[cw6j], W9j[], Ylm[], Djpq[] )
    elseif  n ==  20    aw6j = W6j(a, b, X, c, d, p);   bw6j = W6j(c, d, X, ee, f, q);   cw6j = W6j(ee, f, X, b, a, r);   
                        rex  = RacahExpression( [X], Integral[], a + b + c + d + ee + f + p + q + r + X, 2*X+1, Kronecker[], Triangle[], 
                                                W3j[], W6j[aw6j, bw6j, cw6j], W9j[], Ylm[], Djpq[] )
    elseif  n ==  21    aw6j = W6j(a, b, X, c, d, p);   bw6j = W6j(c, d, X, ee, f, q);   cw6j = W6j(ee, f, X, a, b, r);   
                        rex  = RacahExpression( [X], Integral[], 2*X, 2*X+1, Kronecker[], Triangle[], W3j[], W6j[aw6j, bw6j, cw6j], W9j[], Ylm[], Djpq[] )
    elseif  n ==  22    aw6j = W6j(a, b, X, c, d, p);   bw6j = W6j(c, d, X, a, b, Y);   cw6j = W6j(a, b, q, c, d, Y);   
                        rex  = RacahExpression( [X, Y], Integral[], Basic(0), (2*X+1)*(2*Y+1), Kronecker[], Triangle[], W3j[], W6j[aw6j, bw6j, cw6j], W9j[], Ylm[], Djpq[] )
    elseif  n ==  23    aw6j = W6j(a, b, X, c, d, p);   bw6j = W6j(d, c, X, a, b, Y);   cw6j = W6j(a, d, q, b, c, Y);   
                        rex  = RacahExpression( [X, Y], Integral[], X + Y, (2*X+1)*(2*Y+1), Kronecker[], Triangle[], W3j[], W6j[aw6j, bw6j, cw6j], W9j[], Ylm[], Djpq[] )
    elseif  n ==  24    aw6j = W6j(p, X, Z, d, ee, c);   bw6j = W6j(q, Y, Z, d, ee, b);   cw9j = W9j(a, b, p, c, d, X, q, Y, Z);   
                        rex  = RacahExpression( [X, Y, Z], Integral[], Basic(0), (2*X+1)*(2*Y+1)*(2*Z+1), Kronecker[], Triangle[], W3j[], 
                                                W6j[aw6j, bw6j], W9j[cw9j], Ylm[], Djpq[] )
    elseif  n ==  25    aw6j = W6j(c, d, X, p, g, s);   bw6j = W6j(b, d, Y, q, g, t);   cw9j = W9j(a, b, p, c, d, X, q, Y, g);   
                        rex  = RacahExpression( [X, Y], Integral[], X + Y, (2*X+1)*(2*Y+1), Kronecker[], Triangle[], W3j[], W6j[aw6j, bw6j], W9j[cw9j], Ylm[], Djpq[] )
    elseif  n ==  26    aw6j = W6j(c, d, X, g, p, s);   bw6j = W6j(b, d, Y, g, q, t);   cw9j = W9j(a, b, p, c, d, X, q, Y, g);   
                        rex  = RacahExpression( [X, Y], Integral[], Basic(0), (2*X+1)*(2*Y+1), Kronecker[], Triangle[], W3j[], W6j[aw6j, bw6j], W9j[cw9j], Ylm[], Djpq[] )
    elseif  n ==  27    aw6j = W6j(a, b, X, Y, g, h);   bw6j = W6j(c, d, Y, b, h, j);   cw9j = W9j(a, b, X, c, d, Y, ee, f, g);   
                        rex  = RacahExpression( [X, Y], Integral[], Basic(0), (2*X+1)*(2*Y+1), Kronecker[], Triangle[], W3j[], W6j[aw6j, bw6j], W9j[cw9j], Ylm[], Djpq[] )
    elseif  n ==  28    aw6j = W6j(a, h, Y, g, b, s);   bw9j = W9j(a, f, X, d, q, ee, p, c, b);   cw9j = W9j(a, f, X, h, r, ee, Y, g, b);   
                        rex  = RacahExpression( [X, Y], Integral[], Basic(0), (2*X+1)*(2*Y+1), Kronecker[], Triangle[], W3j[], W6j[aw6j], W9j[bw9j, cw9j], Ylm[], Djpq[] )
    elseif  n ==  29    aw6j = W6j(a, b, X, g, Z, Y);   bw9j = W9j(a, b, X, c, d, Z, ee, f, g);   cw9j = W9j(g, b, Y, c, d, Z, h, j, a);   
                        rex  = RacahExpression( [X, Y, Z], Integral[], X + Y, (2*X+1)*(2*Y+1)*(2*Z+1), Kronecker[], Triangle[], W3j[], 
                                                W6j[aw6j], W9j[bw9j, cw9j], Ylm[], Djpq[] )
    elseif  n ==  30    aw6j = W6j(p, X, Z, c, f, d);   bw9j = W9j(a, b, p, c, d, X, q, Y, Z);   cw9j = W9j(b, d, Y, c, f, Z, l, k, q);   
                        rex  = RacahExpression( [X, Y, Z], Integral[], X, (2*X+1)*(2*Y+1)*(2*Z+1), Kronecker[], Triangle[], W3j[], W6j[aw6j], W9j[bw9j, cw9j], Ylm[], Djpq[] )
    elseif  n ==  31    aw9j = W9j(a, b, p, c, d, X, q, Y, Z);   bw9j = W9j(c, d, X, ee, f, Z, g, h, p);   cw9j = W9j(b, d, Y, ee, f, Z, j, k, q);   
                        rex  = RacahExpression( [X, Y, Z], Integral[], Basic(0), (2*X+1)*(2*Y+1)*(2*Z+1), Kronecker[], Triangle[], 
                                                W3j[], W6j[], W9j[aw9j, bw9j, cw9j], Ylm[], Djpq[] )
    # Sum rules for four Wnj symbol
    elseif  n ==  32    aw3j = W3j(j1, j5, j2, m1, m5, -m2);    bw3j = W3j(j2, j6, j3, m2, m6, -m3)
                        cw3j = W3j(j3, j7, j4, m3, m7, -m4);    dw3j = W3j(j4, j8, j1, m4, m8, -m1)
                        rex  = RacahExpression( [m1, m2, m3, m4], Integral[], -m1-m2-m3-m4, Basic(1), Kronecker[], Triangle[], W3j[aw3j, bw3j, cw3j, dw3j], 
                                                W6j[], W9j[], Ylm[], Djpq[] )
    elseif  n ==  33    aw6j = W6j(a, b, X, c, d, p);    bw6j = W6j(c, d, X, ee, f, q);    cw6j = W6j(ee, f, X, g, h, r);    dw6j = W6j(g, h, X, b, a, s); 
                        R    = a + b + c + d + ee + f + g + h + p + q + r + s
                        rex  = RacahExpression( [X], Integral[], R - X, (2*X+1), Kronecker[], Triangle[], W3j[], W6j[aw6j, bw6j, cw6j, dw6j], W9j[], Ylm[], Djpq[] )
    elseif  n ==  34    aw6j = W6j(a, b, X, c, d, p);    bw6j = W6j(c, d, X, ee, f, q);    cw6j = W6j(ee, f, X, g, h, r);    dw6j = W6j(g, h, X, a, b, s); 
                        rex  = RacahExpression( [X], Integral[], Basic(0), (2*X+1), Kronecker[], Triangle[], W3j[], W6j[aw6j, bw6j, cw6j, dw6j], W9j[], Ylm[], Djpq[] )
    elseif  n ==  35    aw6j = W6j(a, b, c, Z, X, Y);    bw6j = W6j(a, ee, f, p, X, Y);    cw6j = W6j(c, f, d, p, Z, X);    dw6j = W6j(b, ee, d, p, Z, Y); 
                        rex  = RacahExpression( [X,Y,Z], Integral[], X + Y + Z, (2*X+1)*(2*Y+1)*(2*Z+1), Kronecker[], Triangle[], W3j[], 
                                                W6j[aw6j, bw6j, cw6j, dw6j], W9j[], Ylm[], Djpq[] )
    elseif  n ==  36    aw6j = W6j(a, b, c, d, X, Y);    bw6j = W6j(a, ee, f, g, X, Y);    cw6j = W6j(c, g, p, f, d, X);    dw6j = W6j(b, g, q, ee, d, Y); 
                        rex  = RacahExpression( [X,Y], Integral[], Basic(0), (2*X+1)*(2*Y+1), Kronecker[], Triangle[], W3j[], 
                                                W6j[aw6j, bw6j, cw6j, dw6j], W9j[], Ylm[], Djpq[] )
    elseif  n ==  37    aw6j = W6j(a, b, c, d, X, Y);    bw6j = W6j(a, ee, f, g, X, Y);    cw6j = W6j(c, f, p, g, d, X);    dw6j = W6j(b, ee, q, g, d, Y); 
                        rex  = RacahExpression( [X,Y], Integral[], X + Y, (2*X+1)*(2*Y+1), Kronecker[], Triangle[], W3j[], W6j[aw6j, bw6j, cw6j, dw6j], W9j[], Ylm[], Djpq[] )
    elseif  n ==  38    aw9j = W9j(X, Y, Z, d, ee, f, a, b, c);    bw6j = W6j(X, Y, Z, f, c, g);    cw6j = W6j(X, a, d, bp, g, c);    dw6j = W6j(Y, b, ee, dp, f, g); 
                        rex  = RacahExpression( [X,Y,Z], Integral[], Basic(0), (2*X+1)*(2*Y+1)*(2*Z+1), Kronecker[], Triangle[], 
                                                W3j[], W6j[bw6j, cw6j, dw6j], W9j[aw9j], Ylm[], Djpq[] )
    elseif  n ==  39    aw9j = W9j(a, b, X, c, d, Y, eep, fp, g);    bw6j = W6j(ee, f, g, X, Y, Z);    cw6j = W6j(a, b, X, f, Z, d);    
                        dw6j = W6j(c, d, Y, Z, ee, a); 
                        rex  = RacahExpression( [X,Y,Z], Integral[], 2*Z, (2*X+1)*(2*Y+1)*(2*Z+1), Kronecker[], Triangle[], W3j[], 
                                                W6j[bw6j, cw6j, dw6j], W9j[aw9j], Ylm[], Djpq[] )
    elseif  n ==  40    aw9j = W9j(a, b, X, c, d, Y, p, q, r);    bw6j = W6j(X, Y, r, j, h, g);    cw6j = W6j(a, b, X, h, g, ee);    dw6j = W6j(c, d, Y, j, g, f); 
                        rex  = RacahExpression( [X,Y], Integral[], Y, (2*X+1)*(2*Y+1), Kronecker[], Triangle[], W3j[], W6j[bw6j, cw6j, dw6j], W9j[aw9j], Ylm[], Djpq[] )
    elseif  n ==  41    aw9j = W9j(a, d, X, b, ee, Y, cp, f, g);    bw9j = W9j(l, a, Z, ee, b, Y, j, c, h);    cw6j = W6j(a, d, X, k, Z, l);    
                        dw6j = W6j(X, Y, g, h, k, Z); 
                        rex  = RacahExpression( [X,Y,Z], Integral[], Z, (2*X+1)*(2*Y+1)*(2*Z+1), Kronecker[], Triangle[], W3j[], W6j[cw6j, dw6j], 
                                                W9j[aw9j, bw9j], Ylm[], Djpq[] )
    elseif  n ==  42    aw9j = W9j(a, b, X, c, d, Y, ee, f, Z);    bw9j = W9j(g, h, X, k, l, Y, f, ee, Z);    cw6j = W6j(a, b, X, g, h, j);    
                        dw6j = W6j(c, d, Y, k, l, jp); 
                        rex  = RacahExpression( [X,Y,Z], Integral[], 2*Y - Z, (2*X+1)*(2*Y+1)*(2*Z+1), Kronecker[], Triangle[], 
                                                    W3j[], W6j[cw6j, dw6j], W9j[aw9j, bw9j], Ylm[], Djpq[] )
    elseif  n ==  43    aw9j = W9j(a, b, X, g, c, q, p, Z, Y);    bw9j = W9j(b, d, fp, c, h, j, Z, Y, p);    cw6j = W6j(a, b, X, d, ee, f);    
                        dw6j = W6j(d, h, Y, q, X, ee); 
                        rex  = RacahExpression( [X,Y,Z], Integral[], Y, (2*X+1)*(2*Y+1)*(2*Z+1), Kronecker[], Triangle[], W3j[], W6j[cw6j, dw6j], 
                                                W9j[aw9j, bw9j], Ylm[], Djpq[] )
    elseif  n ==  44    aw9j = W9j(a, b, X, c, d, Y, p, q, s);    bw9j = W9j(ee, f, X, g, h, Y, r, t, s);    cw6j = W6j(a, b, X, f, ee, k);    
                        dw6j = W6j(c, d, Y, h, g, l); 
                        rex  = RacahExpression( [X,Y], Integral[], X + Y, (2*X+1)*(2*Y+1), Kronecker[], Triangle[], W3j[], W6j[cw6j, dw6j], W9j[aw9j, bw9j], Ylm[], Djpq[] )
    elseif  n ==  45    aw9j = W9j(a, b, X, c, d, Y, t, s, r);    bw9j = W9j(a, b, X, h, j, q, ee, f, Z);    cw9j = W9j(k, l, p, c, d, Y, ee, f, Z);    
                        dw6j = W6j(p, q, r, X, Y, Z); 
                        rex  = RacahExpression( [X,Y,Z], Integral[], Basic(0), (2*X+1)*(2*Y+1)*(2*Z+1), Kronecker[], Triangle[], 
                                                W3j[], W6j[dw6j], W9j[aw9j, bw9j, cw9j], Ylm[], Djpq[] )
    # Sum rules for five Wnj symbol
    elseif  n ==  46    aw3j = W3j(j1, j6, j2, m1, m6, -m2);    bw3j = W3j(j2, j7, j3, m2, m7, -m3);    cw3j = W3j(j3, j8, j4, m3, m8, -m4)
                        dw3j = W3j(j4, j9, j5, m4, m9, -m5);    ew3j = W3j(j5, j10, j1, m5, m10, -m1)
                        rex = RacahExpression( [m1, m2, m3, m4, m5], Integral[], -m1-m2-m3-m4-m5, Basic(1), Kronecker[], Triangle[], 
                                                W3j[aw3j, bw3j, cw3j, dw3j, ew3j], W6j[], W9j[], Ylm[], Djpq[] )
    # Sum rules for six Wnj symbol
    elseif  n ==  47    aw3j = W3j(j1, j7, j2, m1, m7, -m2);    bw3j = W3j(j2, j8, j3, m2, m8, -m3);      cw3j = W3j(j3, j9, j4, m3, m9, -m4)
                        dw3j = W3j(j4, j10, j5, m4, m10, -m5);  ew3j = W3j(j5, j11, j6, m5, m11, -m6);    fw3j = W3j(j6, j12, j1, m6, m12, -m1)
                        rex = RacahExpression( [m1, m2, m3, m4, m5, m6], Integral[], -m1-m2-m3-m4-m5-m6, Basic(1), Kronecker[], Triangle[], 
                                                W3j[aw3j, bw3j, cw3j, dw3j, ew3j, fw3j], W6j[], W9j[], Ylm[], Djpq[] )
    else    error("stop a")
    end
    
    return( rex )
end


"""
`RacahAlgebra.selectRacahIntegral(n::Int64)`  
    ... selects one of various pre-defined Racah integrals as they often occur on the lhs of some sum/integration rule;
        this function has been implemented mainly for test purposes. A rex::RacahExpression is returned. 
        If n = 99, all pre-defined RacahExpression are printed to screen and nothing is returned in this
        case.
"""
function selectRacahIntegral(n::Int64)
    if  n == 99     for  i = 1:20   println("  $i    $(selectRacahIntegral(i))")      end
        return(nothing)
    end
    
    j = Basic(:j);    J = Basic(:J);    m = Basic(:m);    M  = Basic(:M)
    a = Basic(:a);    b = Basic(:b);    c = Basic(:c);    d  = Basic(:d);    ee  = Basic(:ee);    f  = Basic(:f) 
    g = Basic(:g);    h = Basic(:h);    k = Basic(:k);    l  = Basic(:l);    p   = Basic(:p);     q  = Basic(:q);     
    r = Basic(:r);    s = Basic(:s);    t = Basic(:t);    X = Basic(:X);     Y   = Basic(:Y);     Z  = Basic(:Z);
    
    # Integral rules for one Ylm symbol   ... to be adapted
    if      n ==  1     w3j = W3j(j, j, J, m, -m, M)
                        rex = RacahExpression( [m], Integral[], Basic(-m), Basic(1), Kronecker[], Triangle[], [w3j], W6j[], W9j[], Ylm[], Djpq[] )
    else    error("stop a")
    end
    
    return( rex )
end
        


"""
`RacahAlgebra.simplifyDeltas(rex::RacahExpression)`  
    ... simplifies the all Kronecker deltas in rex; a new nrex::RacahExpression is returned, typically with less summation
        indices.
"""
function simplifyDeltas(rex::RacahExpression)
    nrex = rex
    #
    for  delta  in  rex.deltas     
        newSummations = nrex.summations;    newPhase = nrex.phase;     newWeight    = nrex.weight;     newDeltas = nrex.deltas;    
        newTriangles  = nrex.triangles;     newW3js  = nrex.w3js;      newW6js      = nrex.w6js;       newW9js   = nrex.w9js
        newYlms       = nrex.ylms;          newDjpqs = nrex.djpqs;     newIntegrals = nrex.djpqs
        #
        if  delta.i == delta.k  &&  (delta.i in newSummations  ||   -delta.i in newSummations)
            newDeltas = Kronecker[]
            for dt in nrex.deltas   if  dt == delta   else   push!(newDeltas, RacahAlgebra.subs(dt, delta.i, delta.k))   end    end
        elseif  delta.i == delta.k 
            newDeltas = Kronecker[]
            for dt in nrex.deltas   
                if  dt == delta   else   push!(newDeltas, RacahAlgebra.subs(dt, delta.i, delta.k))   end    
            end
        elseif  delta.i in newSummations
            newDeltas = Kronecker[]
            for dt in nrex.deltas   if  dt == delta   else   push!(newDeltas, RacahAlgebra.subs(dt, delta.i, delta.k))   end    end
            newSummations = RacahAlgebra.removeIndex([delta.i], newSummations)
            newPhase      = SymEngine.subs(newPhase,        delta.i, delta.k)
            newWeight     = SymEngine.subs(newWeight,       delta.i, delta.k)
            newTriangles  = RacahAlgebra.subs(newTriangles, delta.i, delta.k)
            newW3js       = RacahAlgebra.subs(newW3js,      delta.i, delta.k)
            newW6js       = RacahAlgebra.subs(newW6js,      delta.i, delta.k)
            newW9js       = RacahAlgebra.subs(newW9js,      delta.i, delta.k)
            newYlms       = RacahAlgebra.subs(newYlms,      delta.i, delta.k)
            newDjpqs      = RacahAlgebra.subs(newDjpqs,     delta.i, delta.k)
            newIntegrals  = RacahAlgebra.subs(newIntegrals, delta.i, delta.k)
        elseif  -delta.i in newSummations
            newDeltas = Kronecker[]
            for dt in nrex.deltas   if  dt == delta   else   push!(newDeltas, RacahAlgebra.subs(dt, -delta.i, -delta.k))   end    end
            newSummations = RacahAlgebra.removeIndex([delta.i], newSummations)
            newPhase      = SymEngine.subs(newPhase,        -delta.i, -delta.k)
            newWeight     = SymEngine.subs(newWeight,       -delta.i, -delta.k)
            newTriangles  = RacahAlgebra.subs(newTriangles, -delta.i, -delta.k)
            newW3js       = RacahAlgebra.subs(newW3js,      -delta.i, -delta.k)
            newW6js       = RacahAlgebra.subs(newW6js,      -delta.i, -delta.k)
            newW9js       = RacahAlgebra.subs(newW9js,      -delta.i, -delta.k)
            newYlms       = RacahAlgebra.subs(newYlms,      -delta.i, -delta.k)
            newDjpqs      = RacahAlgebra.subs(newDjpqs,     -delta.i, -delta.k)
            newIntegrals  = RacahAlgebra.subs(newIntegrals, -delta.i, -delta.k)
        elseif  delta.k in newSummations
            newDeltas = Kronecker[]
            for dt in nrex.deltas   if  dt == delta   else   push!(newDeltas, RacahAlgebra.subs(dt, delta.k, delta.i))   end    end
            newSummations = RacahAlgebra.removeIndex([delta.k], newSummations)
            newPhase      = SymEngine.subs(newPhase,        delta.k, delta.i)
            newWeight     = SymEngine.subs(newWeight,       delta.k, delta.i)
            newTriangles  = RacahAlgebra.subs(newTriangles, delta.k, delta.i)
            newW3js       = RacahAlgebra.subs(newW3js,      delta.k, delta.i)
            newW6js       = RacahAlgebra.subs(newW6js,      delta.k, delta.i)
            newW9js       = RacahAlgebra.subs(newW9js,      delta.k, delta.i)
            newYlms       = RacahAlgebra.subs(newYlms,      delta.k, delta.i)
            newDjpqs      = RacahAlgebra.subs(newDjpqs,     delta.k, delta.i)
            newIntegrals  = RacahAlgebra.subs(newIntegrals, delta.k, delta.i)
        elseif  -delta.k in newSummations
            newDeltas = Kronecker[]
            for dt in nrex.deltas   if  dt == delta   else   push!(newDeltas, RacahAlgebra.subs(dt, delta.k, delta.i))   end    end
            newSummations = RacahAlgebra.removeIndex([delta.k], newSummations)
            newPhase      = SymEngine.subs(newPhase,        -delta.k, -delta.i)
            newWeight     = SymEngine.subs(newWeight,       -delta.k, -delta.i)
            newTriangles  = RacahAlgebra.subs(newTriangles, -delta.k, -delta.i)
            newW3js       = RacahAlgebra.subs(newW3js,      -delta.k, -delta.i)
            newW6js       = RacahAlgebra.subs(newW6js,      -delta.k, -delta.i)
            newW9js       = RacahAlgebra.subs(newW9js,      -delta.k, -delta.i)
            newYlms       = RacahAlgebra.subs(newYlms,      -delta.k, -delta.i)
            newDjpqs      = RacahAlgebra.subs(newDjpqs,     -delta.k, -delta.i)
            newIntegrals  = RacahAlgebra.subs(newIntegrals, -delta.k, -delta.i)
        elseif  delta.i != delta.k  && (delta.i in [Basic(0), Basic(1), Basic(2), Basic(3), Basic(4)]  ||
                                        delta.k in [Basic(0), Basic(1), Basic(2), Basic(3), Basic(4)])
        elseif  delta.i != delta.k
            newSummations = RacahAlgebra.subs(newSummations,delta.k, delta.i)
            newPhase      = SymEngine.subs(newPhase,        delta.k, delta.i)
            newWeight     = SymEngine.subs(newWeight,       delta.k, delta.i)
            ## newDeltas     = RacahAlgebra.subs(newDeltas,    delta.k, delta.i)
            newDeltas = Kronecker[]
            for dt in nrex.deltas   
                if     dt == delta    push!(newDeltas, dt)   
                else   push!(newDeltas, RacahAlgebra.subs(dt, delta.i, delta.k))   
                end    
            end
            newTriangles  = RacahAlgebra.subs(newTriangles, delta.k, delta.i)
            newW3js       = RacahAlgebra.subs(newW3js,      delta.k, delta.i)
            newW6js       = RacahAlgebra.subs(newW6js,      delta.k, delta.i)
            newW9js       = RacahAlgebra.subs(newW9js,      delta.k, delta.i)
            newYlms       = RacahAlgebra.subs(newYlms,      delta.k, delta.i)
            newDjpqs      = RacahAlgebra.subs(newDjpqs,     delta.k, delta.i)
            newIntegrals  = RacahAlgebra.subs(newIntegrals, delta.k, delta.i)
        else    error("stop a")
        end
        
        nrex = RacahExpression( newSummations, newIntegrals, newPhase, newWeight, newDeltas, newTriangles, newW3js, newW6js, newW9js, newYlms, newDjpqs )
    end
    return( nrex )
end


"""
`RacahAlgebra.simplifyTriangles(rex::RacahExpression)`  
    ... simplifies the all Tringle deltas in rex; these triangles can be removed if the same condition is places by one of
        the Wigner symbols in rex. A new nrex::RacahExpression is returned.
"""
function simplifyTriangles(rex::RacahExpression)
    newTriangles = Triangle[];      
    for  triangle  in  rex.triangles
        remove = false
        for  w3j in rex.w3js
            if  triangle == Triangle(w3j.ja, w3j.jb, w3j.jc)    remove = true   end
        end
        for  w6j in rex.w6js
            if  triangle == Triangle(w6j.a, w6j.b, w6j.c)  ||  triangle == Triangle(w6j.d, w6j.e, w6j.c)  ||  
                triangle == Triangle(w6j.d, w6j.b, w6j.f)  ||  triangle == Triangle(w6j.a, w6j.e, w6j.f)
                remove = true
            end
        end
        for  w9j in rex.w9js
            if  triangle == Triangle(w9j.a, w9j.b, w9j.c)  ||  triangle == Triangle(w9j.d, w9j.e, w9j.f)  ||  
                triangle == Triangle(w9j.g, w9j.h, w9j.i)  ||  triangle == Triangle(w9j.a, w9j.d, w9j.g)  ||  
                triangle == Triangle(w9j.b, w9j.e, w9j.h)  ||  triangle == Triangle(w9j.c, w9j.f, w9j.i)
                remove = true
            end
        end
        if  remove   else   push!(newTriangles, triangle)   end
    end
    wa = RacahExpression( rex.summations, rex.integrals, rex.phase, rex.weight, rex.deltas, newTriangles, rex.w3js, rex.w6js, rex.w9js, rex.ylms, rex.djpqs )

    return( wa )
end


"""
`RacahAlgebra.subs(list::Array{Basic,1}, wold::Basic, wnew::Basic)`  
    ... substitutes in list all occasions of wold by wnew; a newList::Array{Basic,1} is returned.
"""
function subs(list::Array{Basic,1}, wold::Basic, wnew::Basic)
    newList = Basic[]
    for  wa  in  list
        push!(newList, SymEngine.subs(wa, wold, wnew))
    end
    return( newList )
end


"""
`RacahAlgebra.subs(delta::RacahAlgebra.Kronecker, subList::Array{Pair{Symbol,Rational{Int64}},1})`  
    ... substitutes in delta the symbols by the corresponding rational numbers in subList; a ww:Kronecker is returned.
"""
function subs(delta::RacahAlgebra.Kronecker, subList::Array{Pair{Symbol,Rational{Int64}},1})
    ww = delta
    for  (a,b) in subList
        wx = SymEngine.subs(ww.i, Basic(a) => Basic(b));   ww = Kronecker(wx, ww.k)
        wx = SymEngine.subs(ww.k, Basic(a) => Basic(b));   ww = Kronecker(ww.i, wx)
    end
    return( ww )
end


"""
`RacahAlgebra.subs(delta::RacahAlgebra.Kronecker, wold::Basic, wnew::Basic)`  
    ... substitutes in triangle all occasions of wold by wnew; a newDelta::RacahAlgebra.Kronecker is returned.
"""
function subs(delta::RacahAlgebra.Kronecker, wold::Basic, wnew::Basic)
    ni = SymEngine.subs(delta.i, wold::Basic, wnew::Basic)
    nk = SymEngine.subs(delta.k, wold::Basic, wnew::Basic)
    return( Kronecker(ni, nk) )
end


"""
`RacahAlgebra.subs(deltas::Array{RacahAlgebra.Kronecker,1}, wold::Basic, wnew::Basic)`  
    ... substitutes in triangle all occasions of wold by wnew; a newDeltass::Array{RacahAlgebra.Kronecker,1} is returned.
"""
function subs(deltas::Array{RacahAlgebra.Kronecker,1}, wold::Basic, wnew::Basic)
    newDeltas = Kronecker[]
    for  delta  in  deltas
        push!(newDeltas, RacahAlgebra.subs(delta, wold, wnew))
    end
    return( newDeltas )
end


"""
`RacahAlgebra.subs(triangle::RacahAlgebra.Triangle, subList::Array{Pair{Symbol,Rational{Int64}},1})`  
    ... substitutes in triangle the symbols by the corresponding rational numbers in subList; a ww:Triangle is returned.
"""
function subs(triangle::RacahAlgebra.Triangle, subList::Array{Pair{Symbol,Rational{Int64}},1})
    ww = triangle
    for  (a,b) in subList
        wx = SymEngine.subs(ww.i, Basic(a) => Basic(b));   ww = Triangle(wx, ww.j, ww.k)
        wx = SymEngine.subs(ww.j, Basic(a) => Basic(b));   ww = Triangle(ww.i, wx, ww.k)
        wx = SymEngine.subs(ww.k, Basic(a) => Basic(b));   ww = Triangle(ww.i, ww.j, wx)
    end
    return( ww )
end


"""
`RacahAlgebra.subs(triangles::Array{RacahAlgebra.Triangle,1}, wold::Basic, wnew::Basic)`  
    ... substitutes in triangle all occasions of wold by wnew; a newTriangles::Array{RacahAlgebra.Triangle,1} is returned.
"""
function subs(triangles::Array{RacahAlgebra.Triangle,1}, wold::Basic, wnew::Basic)
    newTriangles = Triangle[]
    for  triangle  in  triangles
        nja = SymEngine.subs(triangle.ja, wold::Basic, wnew::Basic)
        njb = SymEngine.subs(triangle.jb, wold::Basic, wnew::Basic)
        njc = SymEngine.subs(triangle.jc, wold::Basic, wnew::Basic)
        push!(newTriangles, Triangle(nja, njb, njc))
    end
    return( newTriangles )
end


"""
`RacahAlgebra.subs(w3j::RacahAlgebra.W3j, subList::Array{Pair{Symbol,Rational{Int64}},1})`  
    ... substitutes in w3j the symbols by the corresponding rational numbers in subList; a ww:W3j is returned.
"""
function subs(w3j::RacahAlgebra.W3j, subList::Array{Pair{Symbol,Rational{Int64}},1})
    ww = w3j
    for  (a,b) in subList
        wx = SymEngine.subs(ww.ja, Basic(a) => Basic(b));   ww = W3j(wx, ww.jb, ww.jc, ww.ma, ww.mb, ww.mc) 
        wx = SymEngine.subs(ww.jb, Basic(a) => Basic(b));   ww = W3j(ww.ja, wx, ww.jc, ww.ma, ww.mb, ww.mc) 
        wx = SymEngine.subs(ww.jc, Basic(a) => Basic(b));   ww = W3j(ww.ja, ww.jb, wx, ww.ma, ww.mb, ww.mc) 
        wx = SymEngine.subs(ww.ma, Basic(a) => Basic(b));   ww = W3j(ww.ja, ww.jb, ww.jc, wx, ww.mb, ww.mc) 
        wx = SymEngine.subs(ww.mb, Basic(a) => Basic(b));   ww = W3j(ww.ja, ww.jb, ww.jc, ww.ma, wx, ww.mc) 
        wx = SymEngine.subs(ww.mc, Basic(a) => Basic(b));   ww = W3j(ww.ja, ww.jb, ww.jc, ww.ma, ww.mb, wx) 
    end
    return( ww )
end


"""
`RacahAlgebra.subs(w3js::Array{RacahAlgebra.W3j,1}, wold::Basic, wnew::Basic)`  
    ... substitutes in w3js all occasions of wold by wnew; a newW3js::Array{RacahAlgebra.W3j,1} is returned.
"""
function subs(w3js::Array{RacahAlgebra.W3j,1}, wold::Basic, wnew::Basic)
    newW3js = W3j[]
    for  w3j  in  w3js
        nja = SymEngine.subs(w3j.ja, wold::Basic, wnew::Basic)
        njb = SymEngine.subs(w3j.jb, wold::Basic, wnew::Basic)
        njc = SymEngine.subs(w3j.jc, wold::Basic, wnew::Basic)
        nma = SymEngine.subs(w3j.ma, wold::Basic, wnew::Basic)
        nmb = SymEngine.subs(w3j.mb, wold::Basic, wnew::Basic)
        nmc = SymEngine.subs(w3j.mc, wold::Basic, wnew::Basic)
        push!(newW3js, W3j(nja, njb, njc, nma, nmb, nmc))
    end
    return( newW3js )
end


"""
`RacahAlgebra.subs(w6j::RacahAlgebra.W6j, subList::Array{Pair{Symbol,Rational{Int64}},1})`  
    ... substitutes in w6j the symbols by the corresponding rational numbers in subList; a ww:W6j is returned.
"""
function subs(w6j::RacahAlgebra.W6j, subList::Array{Pair{Symbol,Rational{Int64}},1})
    ww = w6j
    for  (a,b) in subList
        wx = SymEngine.subs(ww.a, Basic(a) => Basic(b));   ww = W6j(wx, ww.b, ww.c, ww.d, ww.e, ww.f)
        wx = SymEngine.subs(ww.b, Basic(a) => Basic(b));   ww = W6j(ww.a, wx, ww.c, ww.d, ww.e, ww.f)
        wx = SymEngine.subs(ww.c, Basic(a) => Basic(b));   ww = W6j(ww.a, ww.b, wx, ww.d, ww.e, ww.f)
        wx = SymEngine.subs(ww.d, Basic(a) => Basic(b));   ww = W6j(ww.a, ww.b, ww.c, wx, ww.e, ww.f)
        wx = SymEngine.subs(ww.e, Basic(a) => Basic(b));   ww = W6j(ww.a, ww.b, ww.c, ww.d, wx, ww.f)
        wx = SymEngine.subs(ww.f, Basic(a) => Basic(b));   ww = W6j(ww.a, ww.b, ww.c, ww.d, ww.e, wx)
    end
    return( ww )
end


"""
`RacahAlgebra.subs(w6js::Array{RacahAlgebra.W6j,1}, wold::Basic, wnew::Basic)`  
    ... substitutes in w6js all occasions of wold by wnew; a newW6js::Array{RacahAlgebra.W6j,1} is returned.
"""
function subs(w6js::Array{RacahAlgebra.W6j,1}, wold::Basic, wnew::Basic)
    newW6js = W6j[]
    for  w6j  in  w6js
        na = SymEngine.subs(w6j.a, wold::Basic, wnew::Basic)
        nb = SymEngine.subs(w6j.b, wold::Basic, wnew::Basic)
        nc = SymEngine.subs(w6j.c, wold::Basic, wnew::Basic)
        nd = SymEngine.subs(w6j.d, wold::Basic, wnew::Basic)
        ne = SymEngine.subs(w6j.e, wold::Basic, wnew::Basic)
        nf = SymEngine.subs(w6j.f, wold::Basic, wnew::Basic)
        push!(newW6js, W6j(na, nb, nc, nd, ne, nf))
    end
    return( newW6js )
end


"""
`RacahAlgebra.subs(w9js::Array{RacahAlgebra.W9j,1}, wold::Basic, wnew::Basic)`  
    ... substitutes in w9js all occasions of wold by wnew; a newW9js::Array{RacahAlgebra.W9j,1} is returned.
"""
function subs(w9js::Array{RacahAlgebra.W9j,1}, wold::Basic, wnew::Basic)
    newW9js = W9j[]
    for  w9j  in  w9js
        na = SymEngine.subs(w9j.a, wold::Basic, wnew::Basic)
        nb = SymEngine.subs(w9j.b, wold::Basic, wnew::Basic)
        nc = SymEngine.subs(w9j.c, wold::Basic, wnew::Basic)
        nd = SymEngine.subs(w9j.d, wold::Basic, wnew::Basic)
        ne = SymEngine.subs(w9j.e, wold::Basic, wnew::Basic)
        nf = SymEngine.subs(w9j.f, wold::Basic, wnew::Basic)
        ng = SymEngine.subs(w9j.g, wold::Basic, wnew::Basic)
        nh = SymEngine.subs(w9j.h, wold::Basic, wnew::Basic)
        ni = SymEngine.subs(w9j.i, wold::Basic, wnew::Basic)
        push!(newW9js, W9j(na, nb, nc, nd, ne, nf, ng, nh, ni))
    end
    return( newW9js )
end


"""
`RacahAlgebra.subs(integrals::Array{RacahAlgebra.Integral,1}, wold::Basic, wnew::Basic)`  
    ... substitutes in integrals all occasions of wold by wnew; a newIntegrals::Array{RacahAlgebra.Integral,1} is returned.
"""
function subs(integrals::Array{RacahAlgebra.Integral,1}, wold::Basic, wnew::Basic)
    newIntegrals = Integral[]
    for  integral  in  integrals
        nvar = SymEngine.subs(integral.var,   wold::Basic, wnew::Basic)
        nlow = SymEngine.subs(integral.low,   wold::Basic, wnew::Basic)
        nvar = SymEngine.subs(integral.up,    wold::Basic, wnew::Basic)
        push!(newIntegrals, Integral(nvar, nlow, nup))
    end
    return( newIntegrals )
end


"""
`RacahAlgebra.subs(ylms::Array{RacahAlgebra.Ylm,1}, wold::Basic, wnew::Basic)`  
    ... substitutes in ylms all occasions of wold by wnew; a newYlms::Array{RacahAlgebra.Ylm,1} is returned.
"""
function subs(ylms::Array{RacahAlgebra.Ylm,1}, wold::Basic, wnew::Basic)
    newYlms = Ylm[]
    for  ylm  in  ylms
        nl     = SymEngine.subs(ylm.l,     wold::Basic, wnew::Basic)
        nm     = SymEngine.subs(ylm.m,     wold::Basic, wnew::Basic)
        ntheta = SymEngine.subs(ylm.theta, wold::Basic, wnew::Basic)
        nphi   = SymEngine.subs(ylm.phi,   wold::Basic, wnew::Basic)
        push!(newYlms, Ylm(nl, nm, ntheta, nphi))
    end
    return( newYlms )
end


"""
`RacahAlgebra.subs(djpqs::Array{RacahAlgebra.Djpq,1}, wold::Basic, wnew::Basic)`  
    ... substitutes in djpqs all occasions of wold by wnew; a newYlms::Array{RacahAlgebra.Djpq,1} is returned.
"""
function subs(djpqs::Array{RacahAlgebra.Djpq,1}, wold::Basic, wnew::Basic)
    newDjpqs = Djpq[]
    for  djpq  in  djpqs
        nj     = SymEngine.subs(djpq.j,     wold::Basic, wnew::Basic)
        np     = SymEngine.subs(djpq.p,     wold::Basic, wnew::Basic)
        nq     = SymEngine.subs(djpq.q,     wold::Basic, wnew::Basic)
        nbeta  = SymEngine.subs(djpq.beta,  wold::Basic, wnew::Basic)
        push!(newDjpqs, Djpq(nj, np, nq, nbeta))
    end
    return( newDjpqs )
end


"""
`RacahAlgebra.subs(rex::RacahAlgebra.RacahExpression, subList::Array{Pair{Symbol,Rational{Int64}},1})`  
    ... substitutes in rex the symbols by the corresponding rational numbers in subList; a ww:RacahExpression is returned.
"""
function subs(rex::RacahAlgebra.RacahExpression, subList::Array{Pair{Symbol,Rational{Int64}},1})
    summations = rex.summations;  phase     = rex.phase;    weight = rex.weight
    deltas      = Kronecker[];    triangles = Triangle[];   w3js   = W3j[];          w6js  = W6j[];        
    w9js        = W9j[];          integrals = Integral[];   ylms   = Ylm[];          djpqs = Djpq[]
    for  (a,b) in subList
        if   Basic(a) in summations   summations = SymEngine.subs(summations, Basic(a) => Basic(b))   end
        phase      = SymEngine.subs(phase,      Basic(a) => Basic(b))
        weight     = SymEngine.subs(weight,     Basic(a) => Basic(b))
    end
    for  delta    in rex.deltas        push!(deltas,    RacahAlgebra.subs(delta,     subList))    end
    for  triangle in rex.triangles     push!(triangles, RacahAlgebra.subs(triangle,  subList))    end
    for  w3j      in rex.w3js          push!(w3js,      RacahAlgebra.subs(w3j,       subList))    end
    for  w6j      in rex.w6js          push!(w6js,      RacahAlgebra.subs(w6j,       subList))    end
    for  w9j      in rex.w9js          push!(w9js,      RacahAlgebra.subs(w9j,       subList))    end
    for  ylm      in rex.ylms          push!(ylms,      RacahAlgebra.subs(ylm,       subList))    end
    for  djpq     in rex.djpqs         push!(djpqs,     RacahAlgebra.subs(djpq,      subList))    end
    for  integral in rex.integrals     push!(integrals, RacahAlgebra.subs(integral,  subList))    end
    
    ww = RacahExpression( summations, integrals, phase, weight, deltas, triangles, w3js, w6js, w9js, ylms, djpqs)
    
    return( ww )
end


"""
`RacahAlgebra.symmetricForms(w3j::RacahAlgebra.W3j; regge::Bool=false)`  
    ... generates a list of equivalent symmetric forms of the Wigner 3j symbol w3j. There are 12 basic symmetric forms 
        for a 3j-symbol, including the given one, and 72 symmetries due to Regge, including the 12 classical ones.
        A rexList:Array{RacahExpression,1} is returned.
"""
function symmetricForms(w3j::RacahAlgebra.W3j; regge::Bool=false)
    rexList = RacahExpression[]
    deltas  = Kronecker[];    triangles = Triangle[];   w3js = W3j[];   w6js = W6j[];   w9js = W9j[];   ylms = Ylm[];   djpqs = Djpq[];   ints = Integral[]
    sums    = Basic[];    phase = Basic(0);    weight = Basic(1)
    
    if regge    error("stop a")
    else
        j1 = w3j.ja;  j2 = w3j.jb;  j3 = w3j.jc;  m1 = w3j.ma;  m2 = w3j.mb;  m3 = w3j.mc
        push!( rexList, RacahExpression( sums, ints, phase,            weight, deltas, triangles, [W3j(j1,j2,j3, m1,m2,m3)],     w6js, w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase,            weight, deltas, triangles, [W3j(j2,j3,j1, m2,m3,m1)],     w6js, w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase,            weight, deltas, triangles, [W3j(j3,j1,j2, m3,m1,m2)],     w6js, w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase + j1+j2+j3, weight, deltas, triangles, [W3j(j2,j1,j3, m2,m1,m3)],     w6js, w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase + j1+j2+j3, weight, deltas, triangles, [W3j(j1,j3,j2, m1,m3,m2)],     w6js, w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase + j1+j2+j3, weight, deltas, triangles, [W3j(j3,j2,j1, m3,m2,m1)],     w6js, w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase + j1+j2+j3, weight, deltas, triangles, [W3j(j1,j2,j3, -m1,-m2,-m3)],  w6js, w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase + j1+j2+j3, weight, deltas, triangles, [W3j(j2,j3,j1, -m2,-m3,-m1)],  w6js, w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase + j1+j2+j3, weight, deltas, triangles, [W3j(j3,j1,j2, -m3,-m1,-m2)],  w6js, w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase,            weight, deltas, triangles, [W3j(j2,j1,j3, -m2,-m1,-m3)],  w6js, w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase,            weight, deltas, triangles, [W3j(j1,j3,j2, -m1,-m3,-m2)],  w6js, w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase,            weight, deltas, triangles, [W3j(j3,j2,j1, -m3,-m2,-m1)],  w6js, w9js, ylms, djpqs ) )
    end
    
    return( rexList )
end


"""
`RacahAlgebra.symmetricForms(w6j::RacahAlgebra.W6j; regge::Bool=false)`  
    ... generates a list of equivalent symmetric forms of the Wigner 6j symbol w6j. There are 24 basic symmetric forms 
        for a 6j-symbol, including the given one, and 144 symmetries due to Regge, including the 24 classical ones.
        A rexList:Array{RacahExpression,1} is returned.
"""
function symmetricForms(w6j::RacahAlgebra.W6j; regge::Bool=false)
    rexList = RacahExpression[]
    deltas  = Kronecker[];    triangles = Triangle[];   w3js = W3j[];   w6js = W6j[];   w9js = W9j[];   ylms = Ylm[];   djpqs = Djpq[];   ints = Integral[]
    sums    = Basic[];        phase = Basic(0);    weight = Basic(1)
    
    if regge    error("stop a")
    else
        j1 = w6j.a;  j2 = w6j.b;  j3 = w6j.c;  j4 = w6j.d;  j5 = w6j.e;  j6 = w6j.f
        push!( rexList, RacahExpression( sums, ints, phase, weight, deltas, triangles, w3js, [W6j(j1,j2,j3, j4,j5,j6)], w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase, weight, deltas, triangles, w3js, [W6j(j1,j3,j2, j4,j6,j5)], w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase, weight, deltas, triangles, w3js, [W6j(j2,j1,j3, j5,j4,j6)], w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase, weight, deltas, triangles, w3js, [W6j(j2,j3,j1, j5,j6,j4)], w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase, weight, deltas, triangles, w3js, [W6j(j3,j1,j2, j6,j4,j5)], w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase, weight, deltas, triangles, w3js, [W6j(j3,j2,j1, j6,j5,j4)], w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase, weight, deltas, triangles, w3js, [W6j(j1,j5,j6, j4,j2,j3)], w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase, weight, deltas, triangles, w3js, [W6j(j1,j6,j5, j4,j3,j2)], w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase, weight, deltas, triangles, w3js, [W6j(j5,j1,j6, j2,j4,j3)], w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase, weight, deltas, triangles, w3js, [W6j(j5,j6,j1, j2,j3,j4)], w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase, weight, deltas, triangles, w3js, [W6j(j6,j1,j5, j3,j4,j2)], w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase, weight, deltas, triangles, w3js, [W6j(j6,j5,j1, j3,j2,j4)], w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase, weight, deltas, triangles, w3js, [W6j(j4,j5,j3, j1,j2,j6)], w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase, weight, deltas, triangles, w3js, [W6j(j4,j3,j5, j1,j6,j2)], w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase, weight, deltas, triangles, w3js, [W6j(j5,j4,j3, j2,j1,j6)], w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase, weight, deltas, triangles, w3js, [W6j(j5,j3,j4, j2,j6,j1)], w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase, weight, deltas, triangles, w3js, [W6j(j3,j4,j5, j6,j1,j2)], w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase, weight, deltas, triangles, w3js, [W6j(j3,j5,j4, j6,j2,j1)], w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase, weight, deltas, triangles, w3js, [W6j(j4,j2,j6, j1,j5,j3)], w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase, weight, deltas, triangles, w3js, [W6j(j4,j6,j2, j1,j3,j5)], w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase, weight, deltas, triangles, w3js, [W6j(j2,j4,j6, j5,j1,j3)], w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase, weight, deltas, triangles, w3js, [W6j(j2,j6,j4, j5,j3,j1)], w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase, weight, deltas, triangles, w3js, [W6j(j6,j4,j2, j3,j1,j5)], w9js, ylms, djpqs ) )
        push!( rexList, RacahExpression( sums, ints, phase, weight, deltas, triangles, w3js, [W6j(j6,j2,j4, j3,j5,j1)], w9js, ylms, djpqs ) )
    end
    
    return( rexList )
end


"""
`RacahAlgebra.symmetricForms(w9j::RacahAlgebra.W9j; regge::Bool=false)`  
    ... generates a list of equivalent symmetric forms of the Wigner 9-j symbol w9j. There are 72 basic symmetric forms 
        for a 9-j-symbol, including the given one. The keyword regge has no effect since no additional Regge symmetries 
        are known for the Wigner 9-j symbols. A rexList:Array{RacahExpression,1} is returned.
"""
function symmetricForms(w9j::RacahAlgebra.W9j; regge::Bool=false)
    rexList = RacahExpression[]
    deltas  = Kronecker[];    triangles = Triangle[];   w3js = W3j[];   w6js = W6j[];   w9js = W9j[];   ys = Ylm[];   ds = Djpq[];   ints = Integral[]
    sums    = Basic[];    phase = Basic(0);    weight = Basic(1)
    
    if  true
        j11 = w9j.a;  j12 = w9j.b;  j13 = w9j.c;  j21 = w9j.d;  j22 = w9j.e;  j23 = w9j.f;  j31 = w9j.g;  j32 = w9j.h;  j33 = w9j.i  
        phase0 = Basic(0);          phaseS = j11 + j12 + j13 + j21 + j22 + j23 + j31 + j32 + j33
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j11, j12, j13, j21, j22, j23, j31, j32, j33)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j11, j21, j31, j12, j22, j32, j13, j23, j33)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j11, j13, j12, j21, j23, j22, j31, j33, j32)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j11, j21, j31, j13, j23, j33, j12, j22, j32)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j12, j11, j13, j22, j21, j23, j32, j31, j33)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j12, j22, j32, j11, j21, j31, j13, j23, j33)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j12, j13, j11, j22, j23, j21, j32, j33, j31)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j12, j22, j32, j13, j23, j33, j11, j21, j31)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j13, j11, j12, j23, j21, j22, j33, j31, j32)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j13, j23, j33, j11, j21, j31, j12, j22, j32)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j13, j12, j11, j23, j22, j21, j33, j32, j31)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j13, j23, j33, j12, j22, j32, j11, j21, j31)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j11, j12, j13, j31, j32, j33, j21, j22, j23)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j11, j31, j21, j12, j32, j22, j13, j33, j23)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j11, j13, j12, j31, j33, j32, j21, j23, j22)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j11, j31, j21, j13, j33, j23, j12, j32, j22)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j12, j11, j13, j32, j31, j33, j22, j21, j23)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j12, j32, j22, j11, j31, j21, j13, j33, j23)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j12, j13, j11, j32, j33, j31, j22, j23, j21)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j12, j32, j22, j13, j33, j23, j11, j31, j21)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j13, j11, j12, j33, j31, j32, j23, j21, j22)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j13, j33, j23, j11, j31, j21, j12, j32, j22)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j13, j12, j11, j33, j32, j31, j23, j22, j21)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j13, j33, j23, j12, j32, j22, j11, j31, j21)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j21, j22, j23, j11, j12, j13, j31, j32, j33)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j21, j11, j31, j22, j12, j32, j23, j13, j33)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j21, j23, j22, j11, j13, j12, j31, j33, j32)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j21, j11, j31, j23, j13, j33, j22, j12, j32)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j22, j21, j23, j12, j11, j13, j32, j31, j33)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j22, j12, j32, j21, j11, j31, j23, j13, j33)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j22, j23, j21, j12, j13, j11, j32, j33, j31)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j22, j12, j32, j23, j13, j33, j21, j11, j31)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j23, j21, j22, j13, j11, j12, j33, j31, j32)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j23, j13, j33, j21, j11, j31, j22, j12, j32)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j23, j22, j21, j13, j12, j11, j33, j32, j31)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j23, j13, j33, j22, j12, j32, j21, j11, j31)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j21, j22, j23, j31, j32, j33, j11, j12, j13)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j21, j31, j11, j22, j32, j12, j23, j33, j13)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j21, j23, j22, j31, j33, j32, j11, j13, j12)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j21, j31, j11, j23, j33, j13, j22, j32, j12)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j22, j21, j23, j32, j31, j33, j12, j11, j13)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j22, j32, j12, j21, j31, j11, j23, j33, j13)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j22, j23, j21, j32, j33, j31, j12, j13, j11)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j22, j32, j12, j23, j33, j13, j21, j31, j11)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j23, j21, j22, j33, j31, j32, j13, j11, j12)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j23, j33, j13, j21, j31, j11, j22, j32, j12)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j23, j22, j21, j33, j32, j31, j13, j12, j11)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j23, j33, j13, j22, j32, j12, j21, j31, j11)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j31, j32, j33, j11, j12, j13, j21, j22, j23)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j31, j11, j21, j32, j12, j22, j33, j13, j23)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j31, j33, j32, j11, j13, j12, j21, j23, j22)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j31, j11, j21, j33, j13, j23, j32, j12, j22)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j32, j31, j33, j12, j11, j13, j22, j21, j23)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j32, j12, j22, j31, j11, j21, j33, j13, j23)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j32, j33, j31, j12, j13, j11, j22, j23, j21)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j32, j12, j22, j33, j13, j23, j31, j11, j21)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j33, j31, j32, j13, j11, j12, j23, j21, j22)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j33, j13, j23, j31, j11, j21, j32, j12, j22)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j33, j32, j31, j13, j12, j11, j23, j22, j21)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j33, j13, j23, j32, j12, j22, j31, j11, j21)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j31, j32, j33, j21, j22, j23, j11, j12, j13)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j31, j21, j11, j32, j22, j12, j33, j23, j13)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j31, j33, j32, j21, j23, j22, j11, j13, j12)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j31, j21, j11, j33, j23, j13, j32, j22, j12)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j32, j31, j33, j22, j21, j23, j12, j11, j13)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j32, j22, j12, j31, j21, j11, j33, j23, j13)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j32, j33, j31, j22, j23, j21, j12, j13, j11)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j32, j22, j12, j33, j23, j13, j31, j21, j11)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j33, j31, j32, j23, j21, j22, j13, j11, j12)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phaseS, weight, deltas, triangles, w3js, w6js, [W9j(j33, j23, j13, j31, j21, j11, j32, j22, j12)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j33, j32, j31, j23, j22, j21, j13, j12, j11)], ys, ds ) )
        push!( rexList, RacahExpression( sums, ints, phase0, weight, deltas, triangles, w3js, w6js, [W9j(j33, j23, j13, j32, j22, j12, j31, j21, j11)], ys, ds ) )
    end
    
    return( rexList )
end


"""
`RacahAlgebra.symmetricForms(ylm::RacahAlgebra.Ylm)`  
    ... generates a list of equivalent symmetric forms of the spherical harmonic ylm. There are 8 basic symmetric forms 
        for a spherical harmonic, including the given one. A rexList:Array{RacahExpression,1} is returned.
"""
function symmetricForms(ylm::RacahAlgebra.Ylm)
    rexList = RacahExpression[]
    deltas  = Kronecker[];    triangles = Triangle[];   w3js = W3j[];   w6js = W6j[];   w9js = W9j[];   ylms = Ylm[];   djpqs = Djpq[];   ints = Integral[]
    sums    = Basic[];        phase = Basic(0);         weight = Basic(1)
    
    l = ylm.l;  m = ylm.m;   theta = ylm.theta;    phi = ylm.phi
    push!( rexList, RacahExpression( sums, ints, phase,       weight, deltas, triangles, w3js, w6js, w9js, [Ylm(l,m,theta,phi)],                      djpqs ) )
    push!( rexList, RacahExpression( sums, ints, phase + m,   weight * exp(2*im*m*phi), deltas, triangles, w3js, w6js, w9js, [Ylm(l,-m,theta,phi)],   djpqs ) )
    push!( rexList, RacahExpression( sums, ints, phase + m,   weight,                   deltas, triangles, w3js, w6js, w9js, [Ylm(l,m,-theta,phi)],   djpqs ) )
    push!( rexList, RacahExpression( sums, ints, phase,       weight * exp(2*im*m*phi), deltas, triangles, w3js, w6js, w9js, [Ylm(l,-m,-theta,phi)],  djpqs ) )
    push!( rexList, RacahExpression( sums, ints, phase + m,   weight,                   deltas, triangles, w3js, w6js, w9js, [Ylm(l,-m,theta,-phi)],  djpqs ) )
    push!( rexList, RacahExpression( sums, ints, phase,       weight * exp(2*im*m*phi), deltas, triangles, w3js, w6js, w9js, [Ylm(l,m, theta,-phi)],  djpqs ) )
    push!( rexList, RacahExpression( sums, ints, phase,       weight,                   deltas, triangles, w3js, w6js, w9js, [Ylm(l,-m,-theta,-phi)], djpqs ) )
    push!( rexList, RacahExpression( sums, ints, phase + m,   weight * exp(2*im*m*phi), deltas, triangles, w3js, w6js, w9js, [Ylm(l,m, -theta,-phi)], djpqs ) )

    return( rexList )
end


"""
`RacahAlgebra.symmetricForms(djpq::RacahAlgebra.Djpq)`  
    ... generates a list of equivalent symmetric forms of the small rotation matrix djpq. There are 12 basic symmetric forms 
        for a rotation matix, including the given one. A rexList:Array{RacahExpression,1} is returned.
"""
function symmetricForms(djpq::RacahAlgebra.Djpq)
    rexList = RacahExpression[]
    deltas  = Kronecker[];    triangles = Triangle[];   w3js = W3j[];   w6js = W6j[];   w9js = W9j[];   ylms = Ylm[];   djpqs = Djpq[];   ints = Integral[]
    sums    = Basic[];    phase = Basic(0);    weight = Basic(1);   pi = Basic(:pi)
    
    j = djpq.j;  p = djpq.p;  q = djpq.q;   beta = djpq.beta
    push!( rexList, RacahExpression( sums, ints, phase,             weight, deltas, triangles, w3js, w6js, w9js, ylms, [Djpq( j, p, q, beta)] ) )
    push!( rexList, RacahExpression( sums, ints, phase + p - q,     weight, deltas, triangles, w3js, w6js, w9js, ylms, [Djpq( j,-p,-q, beta)] ) )
    push!( rexList, RacahExpression( sums, ints, phase + p - q,     weight, deltas, triangles, w3js, w6js, w9js, ylms, [Djpq( j, q, p, beta)] ) )
    push!( rexList, RacahExpression( sums, ints, phase,             weight, deltas, triangles, w3js, w6js, w9js, ylms, [Djpq( j,-q,-p, beta)] ) )
    push!( rexList, RacahExpression( sums, ints, phase,             weight, deltas, triangles, w3js, w6js, w9js, ylms, [Djpq( j, q, p,-beta)] ) )
    push!( rexList, RacahExpression( sums, ints, phase + q - p,     weight, deltas, triangles, w3js, w6js, w9js, ylms, [Djpq( j, q, p, beta)] ) )
    push!( rexList, RacahExpression( sums, ints, phase - j - p,     weight, deltas, triangles, w3js, w6js, w9js, ylms, [Djpq( j,-p, q, beta - pi)] ) )
    push!( rexList, RacahExpression( sums, ints, phase - j -2p -q,  weight, deltas, triangles, w3js, w6js, w9js, ylms, [Djpq( j, p,-q, beta - pi)] ) )
    push!( rexList, RacahExpression( sums, ints, phase - 2j,        weight, deltas, triangles, w3js, w6js, w9js, ylms, [Djpq( j, p, q, beta + 2pi)] ) )
    push!( rexList, RacahExpression( sums, ints, phase - 2j,        weight, deltas, triangles, w3js, w6js, w9js, ylms, [Djpq( j, p, q, beta - 2pi)] ) )
    push!( rexList, RacahExpression( sums, ints, phase + j - q,     weight, deltas, triangles, w3js, w6js, w9js, ylms, [Djpq( j, p,-q, beta + pi)] ) )
    push!( rexList, RacahExpression( sums, ints, phase + j -p -2q,  weight, deltas, triangles, w3js, w6js, w9js, ylms, [Djpq( j,-p, q, beta + pi)] ) )
    
    return( rexList )
end


"""
`RacahAlgebra.testRecursions(; short::Bool=true)`  
    ... tests the implemented recursion by just comparing comparing the number of Wigner symbols; this does not include tests on 
        the proper phase nor the algebraic factors of the Racah expression(s). The success::Bool of these tests is returned.
"""
function testRecursions(; short::Bool=true)
    success  = true;   error("nothing tested yet")
    # For each rule in nnList, nwnjList gives to correct number of Wigner symbols after the evaluation
    nnList   =        [  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]
    nwnjList =        [  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  2,  2,  0,  0,  1,  2,  2,  2,  1,  1,  0,  1,  1,  1,  1,  3,  1,  1]
    append!(nnList,   [ 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47])
    append!(nwnjList, [  1,  4,  3,  2,  1,  2,  1,  2,  0,  2,  2,  2,  2,  3,  3,  6,  8])
    #
    for (n, nn) in enumerate(nnList)
        rex = RacahAlgebra.selectRacahExpression(nn);    nrex = RacahAlgebra.evaluate(rex);      nwnjs = RacahAlgebra.countWignerSymbols(nrex)
        if  nwnjs != nwnjList[n]   println(">> Test fails for rule  $nn :");    @show rex;    @show nrex;    @show nwnjs     end
        success = success && (nwnjs == nwnjList[n])
    end
    
    return( success )
end


"""
`RacahAlgebra.testSpecialValuesW3j(; short::Bool=true)`  
    ... tests the implemented special values of the Wigner 3j symbols by just comparing comparing the number of Wigner symbols,
        which must be 0 in all cases. The success::Bool of these tests is returned.
"""
function testSpecialValuesW3j(; short::Bool=true)
    success  = true;   
    #
    for n = 2:16
        w3j = RacahAlgebra.selectW3j(n);;    nrex = RacahAlgebra.evaluate(w3j);      nwnjs = RacahAlgebra.countWignerSymbols(nrex)
        if  nwnjs != 0   println(">> Test fails for special value rule  $n :");      @show w3j;    @show nrex;    @show nwnjs     end
        success = success && (nwnjs == 0)
    end
    
    return( success )
end


"""
`RacahAlgebra.testSpecialValuesW6j(; short::Bool=true)`  
    ... tests the implemented special values of the Wigner 6j symbols by just comparing comparing the number of Wigner symbols,
        which must be 0 in all cases. The success::Bool of these tests is returned.
"""
function testSpecialValuesW6j(; short::Bool=true)
    success  = true;   
    #
    for n = 1:22
        w6j = RacahAlgebra.selectW6j(n);;    nrex = RacahAlgebra.evaluate(w6j);      nwnjs = RacahAlgebra.countWignerSymbols(nrex)
        if  nwnjs != 0   println(">> Test fails for special value rule  $n :");      @show w6j;    @show nrex;    @show nwnjs     end
        success = success && (nwnjs == 0)
    end
    
    return( success )
end


"""
`RacahAlgebra.testSpecialValuesW9j(; short::Bool=true)`  
    ... tests the implemented special values of the Wigner 9j symbols by just comparing comparing the number of Wigner symbols,
        which must be 0 in all cases. The success::Bool of these tests is returned.
"""
function testSpecialValuesW9j(; short::Bool=true)
    success  = true;   
    #
    for n = 1:2
        w9j = RacahAlgebra.selectW6j(n);;    nrex = RacahAlgebra.evaluate(w9j);      nwnjs = RacahAlgebra.countWignerSymbols(nrex)
        if  nwnjs != 0   println(">> Test fails for special value rule  $n :");      @show w6j;    @show nrex;    @show nwnjs     end
        success = success && (nwnjs == 0)
    end
    
    return( success )
end


"""
`RacahAlgebra.testSpecialValuesYlm(; short::Bool=true)`  
    ... tests the implemented special values of the spherical harmonics by just comparing the number of spherical harmonics,
        which must be 0 in all cases. The success::Bool of these tests is returned.
"""
function testSpecialValuesYlm(; short::Bool=true)
    success  = true;   
    #
    for n = 1:2
        error("Not yet implemented.")
        w9j = RacahAlgebra.selectW6j(n);;    nrex = RacahAlgebra.evaluate(w9j);      nwnjs = RacahAlgebra.countWignerSymbols(nrex)
        if  nwnjs != 0   println(">> Test fails for special value rule  $n :");      @show w6j;    @show nrex;    @show nwnjs     end
        success = success && (nwnjs == 0)
    end
    
    return( success )
end


"""
`RacahAlgebra.testSpecialValuesDjpq(; short::Bool=true)`  
    ... tests the implemented special values of the small Wigner matrices by just comparing the number of these matrices,
        which must be 0 in all cases. The success::Bool of these tests is returned.
"""
function testSpecialValuesDjpq(; short::Bool=true)
    success  = true;   
    #
    for n = 1:2
        error("Not yet implemented.")
        w9j = RacahAlgebra.selectW6j(n);;    nrex = RacahAlgebra.evaluate(w9j);      nwnjs = RacahAlgebra.countWignerSymbols(nrex)
        if  nwnjs != 0   println(">> Test fails for special value rule  $n :");      @show w6j;    @show nrex;    @show nwnjs     end
        success = success && (nwnjs == 0)
    end
    
    return( success )
end


"""
`RacahAlgebra.testSumRules(; short::Bool=true)`  
    ... tests the implemented sum rules by just comparing comparing the number of Wigner symbols; this does not include tests on 
        the proper phase nor the algebraic factors of the Racah expression. The success::Bool of these tests is returned.
"""
function testSumRules(; short::Bool=true)
    success  = true
    # For each rule in nnList, nwnjList gives to correct number of Wigner symbols after the evaluation
    nnList   =        [  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]
    nwnjList =        [  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  2,  2,  0,  0,  1,  2,  2,  2,  1,  1,  0,  1,  1,  1,  1,  3,  1,  1]
    append!(nnList,   [ 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47])
    append!(nwnjList, [  1,  4,  3,  2,  1,  2,  1,  0,  0,  2,  2,  2,  2,  3,  3,  6,  8])
    #
    for (n, nn) in enumerate(nnList)
        if  n in [46,47]    continue    end
        rex = RacahAlgebra.selectRacahExpression(nn);    nrex = RacahAlgebra.evaluate(rex);      nwnjs = RacahAlgebra.countWignerSymbols(nrex)
        if  nwnjs != nwnjList[n]   println(">> Test fails for rule  $nn :");    @show rex;    @show nrex;    @show nwnjs     end
        success = success && (nwnjs == nwnjList[n])
        if  !(nwnjs == nwnjList[n])     println("\n**** Test for rule $n  with nwnjs = $(nwnjs)  !=== $(nwnjList[n])")       end
    end
    
    return( success )
end


"""
`RacahAlgebra.testIntegralRules(; short::Bool=true)`  
    ... tests the implemented integral rules by just comparing comparing the number of integrals symbols; this does not include tests on 
        the proper phase nor the algebraic factors of the Racah expression. The success::Bool of these tests is returned.
"""
function testIntegralRules(; short::Bool=true)
    success  = true
    error("Not yet implemented.")
    # For each rule in nnList, nwnjList gives to correct number of Wigner symbols after the evaluation
    nnList   =        [  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]
    nwnjList =        [  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  2,  2,  0,  0,  1,  2,  2,  2,  1,  1,  0,  1,  1,  1,  1,  3,  1,  1]
    append!(nnList,   [ 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47])
    append!(nwnjList, [  1,  4,  3,  2,  1,  2,  1,  0,  0,  2,  2,  2,  2,  3,  3,  6,  8])
    #
    for (n, nn) in enumerate(nnList)
        if  n in [46,47]    continue    end
        rex = RacahAlgebra.selectRacahExpression(nn);    nrex = RacahAlgebra.evaluate(rex);      nwnjs = RacahAlgebra.countWignerSymbols(nrex)
        if  nwnjs != nwnjList[n]   println(">> Test fails for rule  $nn :");    @show rex;    @show nrex;    @show nwnjs     end
        success = success && (nwnjs == nwnjList[n])
        if  !(nwnjs == nwnjList[n])     println("\n**** Test for rule $n  with nwnjs = $(nwnjs)  !=== $(nwnjList[n])")       end
    end
    
    return( success )
end


"""
`RacahAlgebra.zeroPhases(rex::RacahAlgebra.RacahExpression)`  
    ... determines a number of phase terms that do not modify the overall phase of the given Racah expression but may
        allow to elimate certain variables as required for sum-rule evaluation; a phaseList::Basic[] is returned.
"""
function zeroPhases(rex::RacahAlgebra.RacahExpression)
    phaseList = Basic[]
    #
    # Add such 'zero' phase in turn for the triangular deltas, the Wigner 3-j and 6-j symbols and, perhaps,
    # later for others
    for  triangle in rex.triangles    push!(phaseList, 2triangle.ja + 2triangle.jb + 2triangle.jc)        end
    for  w3j      in rex.w3js         push!(phaseList, 2w3j.ja + 2w3j.jb + 2w3j.jc)    
                                        push!(phaseList,  w3j.ma +  w3j.mb +  w3j.mc)     
                                        push!(phaseList, 2w3j.ja + 2w3j.ma) 
                                        push!(phaseList, 2w3j.jb + 2w3j.mb) 
                                        push!(phaseList, 2w3j.jc + 2w3j.mc)                                 end
    for  w6j      in rex.w6js         push!(phaseList, 2w6j.a  + 2w6j.b + 2w6j.c) 
                                        push!(phaseList, 2w6j.a  + 2w6j.e + 2w6j.f) 
                                        push!(phaseList, 2w6j.d  + 2w6j.b + 2w6j.f) 
                                        push!(phaseList, 2w6j.d  + 2w6j.e + 2w6j.c)                         end
    phaseList = unique(phaseList)
    
    return( phaseList )
end

include("module-RacahAlgebra-inc-special.jl")
include("module-RacahAlgebra-inc-sumrules.jl")
include("module-RacahAlgebra-inc-integralrules.jl")

end # module

