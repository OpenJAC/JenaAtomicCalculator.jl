


#################################################################################################################################
#################################################################################################################################


"""
`struct  Basics.Eigen`  
    ... defines a simple struct to communicate eigenvalues and eigenvectors if different diagonalization procedures are used.

    + values   ::Array{Float64,1}            ... List of eigenvalues.
    + vectors  ::Array{Vector{Float64},1}    ... List of eigenvectors.
"""
struct  Eigen
    values     ::Array{Float64,1}
    vectors    ::Array{Vector{Float64},1}
end


"""
`Basics.isEqual(a::Eigen, b::Eigen, accuracy::Float64)`  
    ... determines whether a == b with the requested accuracy; a value::Bool is returned. The test is made element-wise on 
        abs(a.number - b.number) > accuracy.
"""
function isEqual(a::Eigen, b::Eigen, accuracy::Float64)
    # Check the dimensions
    if  accuracy <= 0.                          error("Requires accuracy > 0.")                                    end
    if  size(a.values) != size(b.values)   ||   size(a.vectors) != size(b.vectors)                return( false )  end
    # Check the eigenvalues and vectors
    for  i = 1:length(a.values)    
        if  abs(a.values[i] - b.values[i]) > accuracy                                             return( false )  end    
        for  j = 1:length(a.vectors[i])   if  abs(a.vectors[i][j] - b.vectors[i][j]) > accuracy   return( false )  end    end     
    end
    return( true )
end 

#################################################################################################################################
#################################################################################################################################


"""
`struct  Basics.EmStokes`  
    ... defines a type to maintain the (computed) given Stokes parameter for the polarization of emitted photons or electrons.

    + P1              ::EmProperty    ... Stokes P1 parameter in Coulomb and Babushkin gauge
    + P2              ::EmProperty    ... Stokes P2 parameter.
    + P3              ::EmProperty    ... Stokes P3 parameter.
"""
struct  EmStokes 
    P1                ::EmProperty
    P2                ::EmProperty
    P3                ::EmProperty
end 

export EmStokes


# `Base.show(io::IO, ps::EmStokes)`  ... prepares a proper printout of the variable ps::EmStokes.
function Base.show(io::IO, ps::EmStokes) 
    sa = Base.string(ps);                print(io, sa)
end


# `Base.string(ps::EmStokes)`  ... provides a String notation for the variable ps::EmStokes.
function Base.string(ps::EmStokes)
    sa = "Stokes parameter P1 = $(ps.P1),  P2 = $(ps.P2),  P3 = $(ps.P3)  "
    return( sa )
end

#################################################################################################################################
#################################################################################################################################


"""
`struct  Basics.ExpStokes`  
    ... defines a type to maintain the (experimentally) given Stokes parameter for the polarization of incoming 
        photons or electrons.

    + P1              ::Float64    ... Stokes P1 parameter.
    + P2              ::Float64    ... Stokes P2 parameter.
    + P3              ::Float64    ... Stokes P3 parameter.
"""
struct  ExpStokes 
    P1                ::Float64
    P2                ::Float64
    P3                ::Float64
end 

export ExpStokes


"""
`Basics.ExpStokes()`  ... constructor for an (empty or unpolarized) Stokes vector.
"""
function ExpStokes()
    ExpStokes(0., 0., 0.)
end


# `Base.show(io::IO, ps::ExpStokes)`  ... prepares a proper printout of the variable ps::ExpStokes.
function Base.show(io::IO, ps::ExpStokes) 
    sa = Base.string(ps);                print(io, sa)
end


# `Base.string(ps::ExpStokes)`  ... provides a String notation for the variable ps::ExpStokes.
function Base.string(ps::ExpStokes)
    sa = "exp. Stokes parameters P1 = $(ps.P1),  P2 = $(ps.P2),  P3 = $(ps.P3)"
    return( sa )
end

#################################################################################################################################
#################################################################################################################################


"""
`struct  Basics.LevelKey`  ... defines a struct for identifying a level by its symmetry, energy, etc.

    + sym          ::LevelSymmetry    ... Level symmetry.
    + index        ::Int64            ... Index of this level in its original multiplet, or 0.
    + energy       ::Float64          ... Total energy.
    + relativeOcc  ::Float64          ... Relative occupation of this level, if part of some multiplet/cascade.
"""
struct  LevelKey
    sym            ::LevelSymmetry
    index          ::Int64    
    energy         ::Float64 
    relativeOcc    ::Float64 
end

export  LevelKey


"""
`Basics.LevelKey()`  ... constructor for an empty instance of a LevelKey.
"""
function  LevelKey()
    LevelKey( LevelSymmetry(0,Basics.plus), 0, 0., 0.)    
end


# `Base.show(io::IO, key::LevelKey)`  ... prepares a proper printout of the variable key::LevelKey.
function Base.show(io::IO, key::LevelKey) 
    print(io, string(key) )
end


# `Base.string(key::LevelKey)`  ... provides a proper printout of the variable key::LevelKey.
function Base.string(key::LevelKey) 
    return( "LevelKey[$(key.sym) ($(key.index)), en=$(key.energy), occ=$(key.relativeOcc)]")
end

#################################################################################################################################
#################################################################################################################################


"""
`struct  Basics.LineKey`  
    ... defines a struct for identifying a line by the keys of the initial and final level, its transition energy,
        decay strength, etc.

    + iLevelKey    ::LevelKey         ... Key of the initial level.
    + fLevelKey    ::LevelKey         ... Key of the final level.
    + energy       ::Float64          ... Total energy.
    + strength     ::Float64          ... Strength of the given line.
"""
struct  LineKey
    iLevelKey      ::LevelKey  
    fLevelKey      ::LevelKey
    energy         ::Float64   
    strength       ::Float64 
end

export LineKey


"""
`Basics.LineKey()`  ... constructor for an empty instance of a LineKey.
"""
function  LineKey()
    LineKey( LevelKey(), LevelKey(), 0., 0.)    
end


# `Base.show(io::IO, key::LineKey)`  ... prepares a proper printout of the variable key::LineKey.
function Base.show(io::IO, key::LineKey) 
    print(io, string(key) )
end


# `Base.string(key::LineKey)`  ... provides a proper printout of the variable key::LevelKey.
function Base.string(key::LineKey) 
    return( "LineKey[i=$(key.iLevelKey) --> f=$(key.fLevelKey), en=$(key.energy), strength=$(key.strength)]")
end

#################################################################################################################################
#################################################################################################################################


"""
`struct  ScalarProperty{T}`  
    ... defines a type to maintain a scalar function dependence f(x) but for different types input and function values.

    + arg         ::Float64    ... Argument of function f(x).
    + value       ::T          ... Function value f(x).
"""
struct  ScalarProperty{T} 
    arg           ::Float64
    value         ::T
end 


"""
`Basics.ScalarProperty(wa::Float64)`  ... constructor for an `constant` instance of ScalarProperty{T}.
"""
function ScalarProperty(wa::Float64)
    ScalarProperty{T}(wa, 0.)
end


# `Base.show(io::IO, property::ScalarProperty)`  ... prepares a proper printout of the variable property::ScalarProperty.
function Base.show(io::IO, property::ScalarProperty) 
    sa = Base.string(property);                print(io, sa)
end


# `Base.string(property::ScalarProperty)`  ... provides a String notation for the variable property::ScalarProperty.
function Base.string(property::ScalarProperty)
    sa = "f($(property.arg)) = $(property.value)"
    return( sa )
end

#################################################################################################################################
#################################################################################################################################


"""
`struct  Basics.SolidAngle`  ... defines a type for a solid angle Omega = (theta, phi).

    + theta  ::Float64       ... polar angle theta with regard to the z-axis; 0 <= theta <= pi.
    + phi    ::Float64       ... azimuthal angle phi (with regard to the x-z plane in mathematical positive direction); 0 <= phi <= 2 pi.
"""
struct  SolidAngle
    theta    ::Float64
    phi      ::Float64
end


# `Base.show(io::IO, Omega::SolidAngle)`  ... prepares a proper printout of the variable Omega::SolidAngle.
function Base.show(io::IO, Omega::SolidAngle)
    sa = "Omega == (" * string(Omega.theta) * ", " * string(Omega.phi) * ")"
    print(io, sa)
end

export  SolidAngle

#################################################################################################################################
#################################################################################################################################


"""
`struct  Basics.TensorComp`  ... defines a type for a component of the statistical tensor as associated with an atomic or ionic level.

    + k               ::Int64          ... rank of the tensor component
    + q               ::Int64          ... q-component
    + comp            ::ComplexF64     ... value of the tensor component
"""
struct  TensorComp 
    k                 ::Int64
    q                 ::Int64 
    comp              ::ComplexF64
end 

export TensorComp


# `Base.show(io::IO, st::TensorComp)`  ... prepares a proper printout of the variable st::TensorComp.
function Base.show(io::IO, st::TensorComp) 
    sa = Base.string(st);                print(io, sa)
end


# `Base.string(st::TensorComp)`  ... provides a String notation for the variable st::TensorComp.
function Base.string(st::TensorComp)
    sa = "Statistical tensor component rho_($(st.k),$(st.q)) = $(st.comp)  "
    return( sa )
end


"""
`Basics.tensorComp(k::Int64, q::Int64, tcomps::Array{TensorComp,1};  withZeros::Bool = false)`  
    ... returns the value of the statistical tensor component rho_kq if it is contained in tcomps.
"""
function tensorComp(k::Int64, q::Int64, tcomps::Array{TensorComp,1};  withZeros::Bool = false)
    for  tc in  tcomps
        if   k == tc.k    &&    q == tc.q    return( tc.comp )    end
    end
    
    if   withZeros    return( ComplexF64(0.) )    else    error("Statistical tensor component not found for k = $k and q = $q ")    end
end

export  tensorComp

#################################################################################################################################
#################################################################################################################################






