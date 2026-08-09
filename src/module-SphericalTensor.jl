
"""
`module JAC.SphericalTensor`  
... a submodel of JAC that help deal with and simplify composed (single- and many-electron) spherical tensor operators,
    matrix elements and spherical amplitudes into a reduced matrix-element representation.
"""
module  SphericalTensor


using   SymEngine,  ..AngularMomentum, ..Basics,  ..Defaults,  ..RacahAlgebra

export  SphericalState, SphericalMatrixElement


"""
`abstract type SphericalTensor.AbstractSphericalTensor` 
    ... defines an abstract type to comprise various spherical (electron) operators and spherical fields.
"""
abstract type  AbstractSphericalTensor                          end


"""
`abstract type SphericalTensor.AbstractSphericalOperator` 
    ... defines an abstract type to comprise various spherical (electron) operators.

    + CkOperator                ... to represent a C^(K) tensor operator.
    + CoulombOperator           ... to represent a V^(0, Coulomb) tensor operator.
    + DipoleOperator            ... to represent a D^(1) dipole operator.
    + MultipoleOperator         ... to represent a O^(L, multipole) operator.
    + TkOperator                ... to represent a general T^(K) tensor operator.
"""
abstract type  AbstractSphericalOperator  <:  AbstractSphericalTensor  end


"""
`struct  SphericalTensor.CkOperator  <: AbstractSphericalOperator`   ... to represent a C^(K) spherical tensor operator.

    + star          ::Bool      ... defines the complex-conjugate form, if true. 
    + k             ::Basic     ... Tensor rank. 
"""
struct   CkOperator  <: AbstractSphericalOperator
    star            ::Bool 
    k               ::Basic
end


"""
`SphericalTensor.CkOperator(k::Basic)`  ... constructor for just setting the rank (star = false).
"""
function CkOperator(k::Basic)
    CkOperator(false, k)
end


function Base.string(op::CkOperator)
    sa = "C^(" * string(op.k);        if op.star   sa = sa * "*"   end;      sa = sa * ")"
    return( sa )
end

function Base.show(io::IO, op::CkOperator)
    sa = string(op) * " spherical tensor operator.";       print(io, sa, "\n")
end


"""
`struct  SphericalTensor.CoulombOperator  <: AbstractSphericalOperator`   ... to represent a V^(0, Coulomb) spherical tensor operator.

    + star          ::Bool      ... defines the complex-conjugate form, if true. 
"""
struct   CoulombOperator  <: AbstractSphericalOperator
    star            ::Bool 
end


"""
`SphericalTensor.CoulombOperator()`  ... constructor for just setting the dipole operator (star = false).
"""
function CoulombOperator()
    CoulombOperator(false)
end


function Base.string(op::CoulombOperator)
    sa = "V^(0, Coulomb";        if op.star   sa = sa * "*"   end;      sa = sa * ")" 
    return( sa )
end

function Base.show(io::IO, op::CoulombOperator)
    sa = string(op) * " rank-0 Coulomb operator.";       print(io, sa, "\n")
end


"""
`struct  SphericalTensor.DipoleOperator  <: AbstractSphericalOperator`   ... to represent a D^(1) spherical tensor operator.

    + star          ::Bool      ... defines the complex-conjugate form, if true. 
"""
struct   DipoleOperator  <: AbstractSphericalOperator
    star            ::Bool 
end


"""
`SphericalTensor.DipoleOperator()`  ... constructor for just setting the dipole operator (star = false).
"""
function DipoleOperator()
    DipoleOperator(false)
end


function Base.string(op::DipoleOperator)
    sa = "D^(1";        if op.star   sa = sa * "*"   end;      sa = sa * ")"
    return( sa )
end

function Base.show(io::IO, op::DipoleOperator)
    sa = string(op) * " spherical dipole operator.";       print(io, sa, "\n")
end


"""
`struct  SphericalTensor.TkOperator  <: AbstractSphericalOperator`   ... to represent a general rank-K T^(K) spherical tensor operator.

    + star          ::Bool      ... defines the complex-conjugate form, if true. 
    + k             ::Basic     ... Tensor rank. 
"""
struct   TkOperator  <: AbstractSphericalOperator
    star            ::Bool 
    k               ::Basic
end


"""
`SphericalTensor.TkOperator(k::Basic)`  ... constructor for just setting the rank (star = false).
"""
function TkOperator(k::Basic)
    TkOperator(false, k)
end


function Base.string(op::TkOperator)
    sa = "T^(" * string(op.k);        if op.star   sa = sa * "*"   end;      sa = sa * ")"
    return( sa )
end

function Base.show(io::IO, op::TkOperator)
    sa = string(op) * " spherical tensor operator.";       print(io, sa, "\n")
end


"""
`abstract type SphericalTensor.AbstractSphericalConstant  <:  AbstractSphericalTensor` 
    ... defines an abstract type to comprise various (parametric-dependent) scalar constants, such as field amplitudes, etc.

    + ScalarConstant            ... to represent a scalar constant to comprise all all factors that do not
                                    affect the tensorial structure of matrix elements and amplitues.
"""
abstract type  AbstractSphericalConstant   <:  AbstractSphericalTensor  end


"""
`struct  SphericalTensor.ScalarConstant  <: AbstractSphericalConstant`   
    ... to represent a scalar constant which does not affect the tensorial evaluation of matrix elements.

    + star          ::Bool      ... defines the complex-conjugate form, if true. 
    + symbol        ::Basic     ... symbol of the constant.
"""
struct   ScalarConstant  <: AbstractSphericalConstant
    star            ::Bool 
    symbol          ::Basic
end


function Base.string(op::ScalarConstant)
    sa = string(op.symbol);        if op.star   sa = sa * "^(*)"   end
    return( sa )
end

function Base.show(io::IO, op::ScalarConstant)
    sa = string(op) * " scalar constant.";       print(io, sa, "\n")
end


"""
`abstract type SphericalTensor.AbstractSphericalField  <:  AbstractSphericalTensor` 
    ... defines an abstract type to comprise various (external and electron-independent) field operators.

    + UVector                   ... to represent a u^(1) polarization (unit) vector.
    + TwistedU                  ... to represent a u^(1, twisted) unit vector.
"""
abstract type  AbstractSphericalField   <:  AbstractSphericalTensor  end


"""
`struct  SphericalTensor.UVector  <: AbstractSphericalField`   ... to represent a u^(1) polarization (unit) vector.

    + star          ::Bool      ... defines the complex-conjugate form, if true. 
"""
struct   UVector  <: AbstractSphericalField
    star            ::Bool 
end


"""
`SphericalTensor.UVector()`  ... constructor for just setting the vector (star = false).
"""
function UVector()
    UVector(false)
end


function Base.string(op::UVector)
    sa = "u^(1";        if op.star   sa = sa * "*"   end;      sa = sa * ")"
    return( sa )
end

function Base.show(io::IO, op::UVector)
    sa = string(op) * " spherical polarization (unit) vector.";       print(io, sa, "\n")
end



"""
`struct  SphericalTensor.TensorProduct  <:  AbstractSphericalTensor` 
    ... defines the coupled tensor product of two spherical tensors [T1^(k1) * T2^(k2)]^(k).

    + star          ::Bool      ... defines the complex-conjugate form, if true. 
    + k             ::Basic     ... (total) tensor rank k.
    + tensor1       ::AbstractSphericalTensor  ... spherical tensor T1.
    + tensor2       ::AbstractSphericalTensor  ... spherical tensor T2.
"""
struct   TensorProduct  <:  AbstractSphericalTensor
    star            ::Bool 
    k               ::Basic
    tensor1         ::AbstractSphericalTensor
    tensor2         ::AbstractSphericalTensor
end


"""
`SphericalTensor.TensorProduct(k::Basic, tensor1::AbstractSphericalTensor, tensor2::AbstractSphericalTensor)`  
    ... constructor for just setting the tensor product (star = false).
"""
function TensorProduct(k::Basic, tensor1::AbstractSphericalTensor, tensor2::AbstractSphericalTensor)
    TensorProduct(false, k, tensor1, tensor2)
end


function Base.string(op::TensorProduct)
    sa = "[" * string(op.tensor1) * " * " * string(op.tensor2) * "]^(" * string(op.k);        if op.star   sa = sa * "*"   end;      sa = sa * ")"
    return( sa )
end

function Base.show(io::IO, op::TensorProduct)
    sa = string(op) * " spherical tensor product.";       print(io, sa, "\n")
end


"""
`abstract type SphericalTensor.AbstractSphericalTensorComp` 
    ... defines an abstract type to comprise various spherical (electron) operators and spherical field tensor components.
"""
abstract type  AbstractSphericalTensorComp     end


"""
`struct  SphericalTensor.SphericalOperatorComp  <: AbstractSphericalTensorComp`   ... to represent a spherical tensor component.

    + tensor        ::AbstractSphericalOperator   ... tensor
    + q             ::Basic                       ... q-component
"""
struct   SphericalOperatorComp  <: AbstractSphericalTensorComp
    tensor          ::AbstractSphericalTensor
    q               ::Basic
end


function Base.string(comp::SphericalOperatorComp)
    sa = string(comp.tensor) * "_" * string(comp.q)
    return( sa )
end

function Base.show(io::IO, comp::SphericalOperatorComp)
    sa = string(comp) * " spherical operator component.";       print(io, sa, "\n")
end


"""
`struct  SphericalTensor.SphericalFieldComp  <: AbstractSphericalTensorComp`   ... to represent a spherical field component.

    + tensor        ::AbstractSphericalField    ... exernal tensor field
    + q             ::Basic                     ... q-component
"""
struct   SphericalFieldComp  <: AbstractSphericalTensorComp
    tensor          ::AbstractSphericalField
    q               ::Basic
end


function Base.string(comp::SphericalFieldComp)
    sa = string(comp.tensor) * "_" * string(comp.q)
    return( sa )
end

function Base.show(io::IO, comp::SphericalFieldComp)
    sa = string(comp) * " spherical field component.";       print(io, sa, "\n")
end



"""
`struct  SphericalTensor.SphericalState`   ... to represent a  |jm > state (spherical representation).

    + star          ::Bool      ... defines the complex-conjugate form, if true. 
    + j             ::Basic     ... total angular momentum
    + m             ::Basic     ... projection of total angular momentum
"""
struct   SphericalState
    star            ::Bool 
    j               ::Basic
    m               ::Basic
end


"""
`SphericalTensor.SphericalState(j::Basic, m::Basic)`  ... constructor for just setting the ket-vector (star = false).
"""
function SphericalState(j::Basic, m::Basic)
    SphericalState(false, j, m)
end


function Base.string(state::SphericalState; plain=false)
    if  plain   sa =       string(state.j) * " " * string(state.m);        if state.star   sa = sa * " *"   end
    else        sa = "|" * string(state.j) * " " * string(state.m);        if state.star   sa = sa * " *"   end;      sa = sa * ">"
    end
    return( sa )
end

function Base.show(io::IO, state::SphericalState)
    sa = string(state) * " spherical state.";       print(io, sa, "\n")
end



"""
`struct  SphericalTensor.SphericalMatrixElement`   
    ... to represent a (full)  <ja ma|| T^(k)_q |jb mb>  matrix element of the spherical operator T^(k).

    + star             ::Bool                         ... defines the complex-conjugate form, if true. 
    + leftState        ::SphericalState               ... (reduced) bra-state.
    + constant         ::ScalarConstant               ... a scalar constant that is treated independently.
    + tensor           ::AbstractSphericalTensor      ... spherical tensor.
    + q                ::Basic                        ... q-quantum number of tensor
    + rightState       ::SphericalState               ... (reduced) ket-state.
"""
struct   SphericalMatrixElement
    star               ::Bool
    leftState          ::SphericalState
    constant           ::ScalarConstant
    tensor             ::AbstractSphericalTensor
    q                  ::Basic 
    rightState         ::SphericalState
end


"""
`SphericalTensor.SphericalMatrixElement(leftState::SphericalState, tensor::AbstractSphericalTensor, q::Basic, rightState::SphericalState)`  
    ... constructor for just setting the matrix element without constant and with star = false.
"""
function SphericalMatrixElement(leftState::SphericalState, tensor::AbstractSphericalTensor, q::Basic, rightState::SphericalState)
    SphericalMatrixElement(false, leftState, ScalarConstant(false, Basic(1)), tensor, q, rightState)
end


function Base.string(me::SphericalMatrixElement)
    sa = "<" * string(me.leftState, plain=true) * " | " 
    if  me.constant.symbol != Basic(1)    sa = sa * string(me.constant)    end
    sa = sa * string(me.tensor) * "_" * string(me.q) * " | " * string(me.rightState, plain=true) * ">" 
    if  me.star                    sa = sa * "^(*)"                 end 
    return( sa )
end

function Base.show(io::IO, me::SphericalMatrixElement)
    sa = string(me) * " reduced matrix element.";       print(io, sa, "\n")
end



"""
`struct  SphericalTensor.ProductMatrixElement`   
    ... to represent a (reduced) <ja|| T^(k) |jb> matrix element of the spherical operator T^(k).

    + leftState        ::SphericalState                    ... (spherical) bra-state.
    + components       ::Array{SphericalOperatorComp,1}    ... (product of) spherical tensors.
    + rightState       ::SphericalState                    ... (spherical) ket-state.
"""
struct   ProductMatrixElement
    leftState          ::SphericalState
    components         ::Array{SphericalOperatorComp,1}
    rightState         ::SphericalState
end


function Base.string(me::ProductMatrixElement)
    sa = "<" * string(me.leftState, plain=true) * " | " 
    for  comp in me.components    sa = sa * "  " * string(comp)      end
    sa = sa * " | " * string(me.rightState, plain=true) * "> " 
    return( sa )
end

function Base.show(io::IO, me::ProductMatrixElement)
    sa = string(me) * " product matrix element.";       print(io, sa, "\n")
end



"""
`struct  SphericalTensor.SphericalAmplitude`   
    ... to represent a (non-reduced, many-electron) spherical amplitude that may include one or several spherical matrix elements, spherical 
        field operator components as well as  various symbols from Racah's algebra. In general, these spherical amplitudes need further to be 
        simplified to (fully) reduced amplitudes prior to the calculations.

    + rex              ::RacahAlgebra.RacahExpression       ... contains all factors from Racah's algebra including summations and delta-factors
    + constants        ::Array{AbstractSphericalConstant,1} ... product (list) of scalar constants.
    + fields           ::Array{SphericalFieldComp,1}        ... product (list) of external field tensors.
    + productMes       ::Array{ProductMatrixElement,1}      ... product (list) of product matrix elements.
"""
struct   SphericalAmplitude
    rex                ::RacahAlgebra.RacahExpression 
    constants          ::Array{AbstractSphericalConstant,1}
    fields             ::Array{SphericalFieldComp,1}
    productMes         ::Array{ProductMatrixElement,1}
end


function Base.string(amp::SphericalAmplitude)
    sa =      "\nRacahExpression:    " * string(amp.rex )
    sa = sa * "\nConstants:        " 
    for cst in amp.constants       sa = sa * "  " * string(cst)      end
    sa = sa * "\nField components: " 
    for field in amp.fields        sa = sa * "  " * string(field)    end
    sa = sa * "\nproductMes:       " 
    for rme in amp.productMes      sa = sa * "  " * string(rme)      end
    return( sa )
end

function Base.show(io::IO, amp::SphericalAmplitude)
    sa = "Spherical amplitude with: " * string(amp);       print(io, sa, "\n")
end



"""
`struct  SphericalTensor.ReducedState`   ... to represent a (reduced) |j> state (spherical representation).

    + j             ::Basic     ... total angular momentum
"""
struct   ReducedState
    j               ::Basic
end


function Base.string(state::ReducedState; plain=false)
    if  plain      sa = string(state.j)    else    sa = "|" * string(state.j) * ">"    end
    return( sa )
end

function Base.show(io::IO, state::ReducedState)
    sa = string(state) * " reduced state.";       print(io, sa, "\n")
end



"""
`struct  SphericalTensor.ReducedMatrixElement`   
    ... to represent a (reduced) <ja|| T^(k) |jb> matrix element of the spherical operator T^(k).

    + leftState        ::ReducedState                      ... (reduced) bra-state.
    + tensors          ::Array{AbstractSphericalTensor,1}  ... (product of) spherical tensors.
    + rightState       ::ReducedState                      ... (reduced) ket-state.
"""
struct   ReducedMatrixElement
    leftState          ::ReducedState
    tensors            ::Array{AbstractSphericalTensor,1}
    rightState         ::ReducedState
end


function Base.string(me::ReducedMatrixElement)
    sa = "<" * string(me.leftState, plain=true) * " || " 
    for  tn in me.tensors   sa = sa * " " * string(tn)      end
    sa = sa * " || " * string(me.rightState, plain=true) * ">" 
    return( sa )
end

function Base.show(io::IO, me::ReducedMatrixElement)
    sa = string(me) * " reduced matrix element.";       print(io, sa, "\n")
end



"""
`struct  SphericalTensor.ReducedAmplitude`   
    ... to represent a (reduced, many-electron) amplitude that may include one or several reduced matrix elements, spherical field operators
        and various symbols from Racah's algebra. In general, these reduced amplitudes need to be simplified prior to the calculations.

    + rex              ::RacahAlgebra.RacahExpression       ... contains all factors from Racah's algebra including summations and delta-factors
    + fields           ::Array{AbstractSphericalField,1}    ... product (list) of external field tensors.
    + reducedMes       ::Array{ReducedMatrixElement,1}      ... product (list) of reduced matrix elements.
"""
struct   ReducedAmplitude
    rex                ::RacahAlgebra.RacahExpression 
    fields             ::Array{AbstractSphericalField,1}
    reducedMes         ::Array{ReducedMatrixElement,1}
end


function Base.string(amp::ReducedAmplitude)
    sa =      "\nRacahExpression:   " * string(amp.rex )
    sa = sa * "\nFields:           " 
    for field in amp.fields    sa = sa * " " * string(field)    end
    sa = sa * "\nreducedMes:       " 
    for rme in amp.reducedMes  sa = sa * " " * string(rme)      end
    return( sa )
end

function Base.show(io::IO, amp::ReducedAmplitude)
    sa = "Reduced amplitude with: " * string(amp);       print(io, sa, "\n")
end




"""
`SphericalTensor.conjugate(me::SphericalTensor.SphericalMatrixElement)`  
    ... conjugates the spherical matrix element by conjugating all parts of the tensor operator and by interchanging the 
        left- and right-hand states. A conjme::SphericalTensor.SphericalMatrixElement is returned.
"""
function conjugate(me::SphericalTensor.SphericalMatrixElement)
    starx       = !me.star
    leftStatex  = me.rightState
    constantx   = ScalarConstant( !me.constant.star,   me.constant.symbol)
    tensorx     = SphericalTensor.conjugate( me.tensor )
    qx          = me.q
    rightStatex = me.leftState
    conjme      = SphericalTensor.SphericalMatrixElement(starx, leftStatex, constantx, tensorx, qx, rightStatex)
    
    return( conjme )    
end


"""
`SphericalTensor.conjugate(tn::SphericalTensor.TensorProduct)`  
    ... conjugates the (abstract) spherical matrix element by conjugating all parts of the tensor operator and by interchanging 
        the individual tensor components of each tensor product.
"""
function conjugate(tn::SphericalTensor.TensorProduct)
    
    starx       = tn.star
    kx          = tn.k
    tensor1x    = SphericalTensor.conjugate( tn.tensor2 )
    tensor2x    = SphericalTensor.conjugate( tn.tensor1 )
    conjtn      = SphericalTensor.TensorProduct(starx, kx, tensor1x, tensor2x)
    
    return( conjtn )    
end

conjugate(tn::SphericalTensor.CkOperator)      = SphericalTensor.CkOperator( !tn.star, tn.k)
conjugate(tn::SphericalTensor.CoulombOperator) = SphericalTensor.CoulombOperator( !tn.star)
conjugate(tn::SphericalTensor.DipoleOperator)  = SphericalTensor.DipoleOperator( !tn.star)
conjugate(tn::SphericalTensor.TkOperator)      = SphericalTensor.TkOperator( !tn.star, tn.k)

conjugate(tn::SphericalTensor.ScalarConstant)  = SphericalTensor.ScalarConstant( !tn.star, tn.symbol)
conjugate(tn::SphericalTensor.UVector)         = SphericalTensor.UVector( !tn.star)


"""
`SphericalTensor.expandSphericalMatrixElements(mes::Array{SphericalTensor.SphericalMatrixElement,1})`  
    ... expand the (product of) spherical matrix elements into a reduced spherical amplitude; it separates the scalar constants from the 
        field operators and the reduced (electron) matrix elements. The Wigner-Eckardt theorem is applied to reduce the spherical
        matrix elements and a Clebsch-Gordan expansion for each product of tensor operators. All additional magnetic quantum numbers,
        that are needed for this expansion, are derived from the corresponding tensor ranks. The Clebsch-Gordan/Wigner symbols are
        'collected' in a Racah expression for further simplification. A sphericalAmp::SphericalTensor.SphericalAmplitude is returned.
"""
function expandSphericalMatrixElements(mes::Array{SphericalTensor.SphericalMatrixElement,1})
    rex         = RacahAlgebra.RacahExpression()
    constants   = AbstractSphericalConstant[]
    fieldComps  = SphericalFieldComp[]
    productMes  = ProductMatrixElement[]
    
    for  mex in mes
        # Always work with the non-conjugated form of spherical ME
        if  mex.star  me = SphericalTensor.conjugate(mex)    else   me = mex     end
        # 'Add' the constants
        if me.constant.symbol != Basic(1)    push!(constants, me.constant)   end
        # Decompose the spherical operator recusively
        wa, wb, wc, wd = SphericalTensor.expandSphericalTensorComponent(me.tensor, me.q)
        rex = rex * wa;    append!(constants, wb);    append!(fieldComps, wc)
        if  me.leftState.star  ||  me.rightState.star   error("SphericalState() assumes non-conjugated spherical states.")     end
        prodMe = SphericalTensor.ProductMatrixElement( me.leftState, wd, me.rightState)
        push!(productMes, prodMe )
    end
    sphericalAmp = SphericalTensor.SphericalAmplitude(rex, constants, fieldComps, productMes)

    return( sphericalAmp )    
end


"""
`SphericalTensor.expandSphericalTensorComponent(tn::SphericalTensor.AbstractSphericalTensor, q::Basic)`  
    ... expand a coupled product of spherical tensors into its (uncoupled) tensor components by means of a Clebsch-Gordan expansion.
        All necessary magnetic quantum numbers are 'derived' from the corresponding tensor ranks: K -> qK, 2 -> q2,  0 -> 0.
        No expansion is made for elementary spherical tensors. The tensor components of this expansion is returned by (the tuple of)
        three variables: (rex::RacahExpression, fields::Array{SphericalTensor.AbstractSphericalField,1}, 
                            operators::Array{SphericalTensor.AbstractSphericalOperator,1}).
"""
function expandSphericalTensorComponent(tn::SphericalTensor.AbstractSphericalTensor, q::Basic)
    rex           = RacahAlgebra.RacahExpression()
    constants     = SphericalTensor.AbstractSphericalConstant[]
    fieldComps    = SphericalTensor.SphericalFieldComp[]
    operatorComps = SphericalTensor.SphericalOperatorComp[]
    
    if  typeof(tn) in [SphericalTensor.CkOperator, SphericalTensor.CoulombOperator, SphericalTensor.DipoleOperator, 
                        SphericalTensor.TkOperator]
        push!(operatorComps, SphericalOperatorComp(tn,q) )
    elseif   typeof(tn) in [SphericalTensor.ScalarConstant]
        push!(constants, tn )
    elseif   typeof(tn) in [SphericalTensor.UVector]
        push!(fieldComps, SphericalFieldComp(tn,q) )
    elseif   typeof(tn) in [SphericalTensor.TensorProduct]
        B0  = Basic(0)
        k1  = getRank(tn.tensor1);   if  k1 == B0 q1 = B0   else   q1 = Basic("q" * string(k1) * "_" *string(floor(10000*rand()))[1:2]) end
        k2  = getRank(tn.tensor2);   if  k2 == B0 q2 = B0   else   q2 = Basic("q" * string(k2) * "_" *string(floor(10000*rand()))[1:2]) end
        kk  = getRank(tn);           if  kk == B0 qq = B0   else   qq = Basic("q" * string(kk) * "_" *string(floor(10000*rand()))[1:2]) end
        
        rx1, cx1, fx1, ox1 = SphericalTensor.expandSphericalTensorComponent(tn.tensor1, q1)
        rx2, cx2, fx2, ox2 = SphericalTensor.expandSphericalTensorComponent(tn.tensor2, q2)
        rex = rex * rx1 * rx2 * RacahAlgebra.ClebschGordanExpansion(k1, q1, k2, q2, kk, qq)
        append!(constants,  cx1);      append!(constants,  cx2)
        append!(fieldComps, fx1);      append!(fieldComps, fx2)
        append!(operatorComps, ox1);   append!(operatorComps, ox2)
    else     error("stop a")
    end
    
    return( rex, constants, fieldComps, operatorComps )    
end


getRank(tn::SphericalTensor.CkOperator)      = tn.k
getRank(tn::SphericalTensor.CoulombOperator) = Basic(0)
getRank(tn::SphericalTensor.DipoleOperator)  = Basic(1)
getRank(tn::SphericalTensor.TkOperator)      = tn.k

getRank(tn::SphericalTensor.ScalarConstant)  = Basic(0)
getRank(tn::SphericalTensor.UVector)         = Basic(1)

getRank(tn::SphericalTensor.TensorProduct)   = tn.k


"""
`SphericalTensor.getWignerEckardtFactor(leftState::SphericalState, k::Basic, q::Basic, rightState::SphericalState)`  
    ... Returns the pre-factor from the Wigner-Eckardt theorem. A rex::RacahAlgebra.RacahExpression is returned.
"""
function getWignerEckardtFactor(leftState::SphericalState, k::Basic, q::Basic, rightState::SphericalState)
    rex = RacahAlgebra.ClebschGordan(rightState.j, rightState.m, k, q, leftState.j, leftState.j)
    return( rex )    
end

end # module

