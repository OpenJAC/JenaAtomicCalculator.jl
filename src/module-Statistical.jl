
"""
`module  JenaAtomicCalculator.Statistical`  
    ... a submodel of JAC for the STATISTICAL TENSORS of an atomic ensemble, i.e. for describing an ensemble of atoms or
        ions that is not isotropic because it was made by something with a direction: a light beam, an ion beam, a
        collision.

        WHAT THE OBJECT IS, in plain terms.  Such an ensemble is not described by "how many atoms" alone, since the
        magnetic sublevels M of one level are then populated unevenly.  The statistical tensor rho_kq is a repackaging
        of those sublevel populations and coherences:  rho_00 is the total population;  rho_1q is the ORIENTATION, a net
        direction of spin, which needs circular light or a polarized beam;  rho_2q is the ALIGNMENT, i.e. whether the
        ensemble is cigar- or pancake-shaped about the axis, with the same number of atoms either way.

        WHY REPACKAGE AT ALL.  Turn the axis and the sublevel populations mix into one another messily, whereas a tensor
        of rank k mixes ONLY among its own 2k+1 partners, through a single Wigner D-matrix.  That is the entire reason
        the object exists, and it is why `Statistical.rotate` belongs here rather than beside any one process.

        THE AXIS IS PART OF THE QUANTITY, and this module records it in the type.  A rho_kq or A_kq is meaningless
        without saying which axis it refers to: A_20 of a system whose symmetry axis lies elsewhere is a TRUE statement
        about z and a FALSE statement about the system, and it comes out small rather than wrong-looking.  Carrying
        `axis` makes that a question the code can answer; `Statistical.invariant` gives the frame-INDEPENDENT magnitude
        of a rank, which is what one usually meant to ask.

        CONVENTION, fixed here once for the whole code, following Blum, *Density Matrix Theory and Applications*:

            rho_kq(a,b)  =  SUM_(M,M')  (-1)^(J_b - M')  <J_a M, J_b -M' | k q>  <a J_a M| rho |b J_b M'>

        so that rho_00 = N / sqrt(2J+1) for N atoms, and the NORMALIZED alignment parameters are A_kq = rho_kq / rho_00.
        With that, A_20 = sqrt(2J+1) rho_20 / N -- which is exactly what `CoulombExcitation` computes inline from
        sigma(M_f), and agreement with that existing site is one of this module's tests rather than an assumption.
"""
module Statistical

using  Printf, ..AngularMomentum, ..Basics, ..ManyElectron


"""
`struct  Statistical.ResonanceR`  
    ... defines a type for a resonance state in the continuum with a well-defined bound-ionic core, one or several
        electrons in the continuum, a width as well as a loss rate due to additional decay processes.  It is used where
        an ensemble EVOLVES IN TIME and a width genuinely matters; an ordinary bound level needs no such dressing and is
        labelled by a `Basics.LevelKey` instead.

    + isBound            ::Bool                ... True if it just refers to a bound state with no electron in the continuum.
    + ionLevel           ::Level               ... Level of the bound-state core.
    + widths             ::Float64             ... Widths of the resonance state.
    + lossRate           ::Float64             ... Loss rate due to processes that are not considered explicitly in the time evolution.
"""
struct ResonanceR
    isBound              ::Bool   
    ionLevel             ::Level
    widths               ::Float64
    lossRate             ::Float64
end 


"""
`Statistical.ResonanceR()`  ... constructor for an `empty` instance of Statistical.ResonanceR().
"""
function ResonanceR()
    ResonanceR( true, Level(), 0., 0. )
end


# `Base.show(io::IO, resonance::Statistical.ResonanceR)`  ... prepares a proper printout of the variable resonance::Statistical.ResonanceR.
function Base.show(io::IO, resonance::Statistical.ResonanceR) 
    println(io, "isBound:              $(resonance.isBound)  ")
    println(io, "ionLevel:             $(resonance.ionLevel)  ")
    println(io, "widths:               $(resonance.widths)  ")
    println(io, "lossRate:             $(resonance.lossRate)  ")
end


"""
`struct  Statistical.Tensor`  
    ... defines a type for a statistical tensor of given rank k, projection q and with respect to two levels.  Both
        levels are named because that is the general definition, rho_kq(alpha J, alpha' J'); the diagonal case a == b is
        the usual one and covers alignment, orientation and every angular distribution.

    + k                  ::AngularJ64          ... Rank of the tensor.
    + q                  ::AngularM64          ... Projection of the tensor rank.
    + a                  ::LevelKey            ... Level a with total symmetry alpha J.
    + b                  ::LevelKey            ... Level b with total symmetry alpha' J'.
    + axis               ::Basics.AbstractQuantizationAxis   ... the axis this component refers to; see the module docstring.
    + value              ::ComplexF64          ... value of the statistical tensor.
"""
struct Tensor
    k                    ::AngularJ64
    q                    ::AngularM64
    a                    ::LevelKey
    b                    ::LevelKey
    axis                 ::Basics.AbstractQuantizationAxis
    value                ::ComplexF64
end 


"""
`Statistical.Tensor()`  ... constructor for an `empty` instance of Statistical.Tensor().
"""
function Tensor()
    Tensor( AngularJ64(0), AngularM64(0), LevelKey(), LevelKey(), Basics.DefaultQuantizationAxis(), ComplexF64(0.) )
end


"""
`Statistical.Tensor(k::Int64, q::Int64, a::LevelKey, value::ComplexF64;` 
                    `axis::Basics.AbstractQuantizationAxis=Basics.DefaultQuantizationAxis())`  
    ... constructor for a DIAGONAL statistical tensor component from integer rank and projection, the usual case; a
        `tensor::Statistical.Tensor` is returned.
"""
function Tensor(k::Int64, q::Int64, a::LevelKey, value::ComplexF64;
                axis::Basics.AbstractQuantizationAxis=Basics.DefaultQuantizationAxis())
    Tensor( AngularJ64(k), AngularM64(q), a, a, axis, value )
end


# `Base.show(io::IO, tensor::Statistical.Tensor)`  ... prepares a proper printout of the variable tensor::Statistical.Tensor.
function Base.show(io::IO, tensor::Statistical.Tensor) 
    println(io, "rho_($(tensor.k),$(tensor.q)) = $(tensor.value)   [about $(tensor.axis)]  ")
end


"""
`Statistical.alignmentParameters(tensors::Array{Statistical.Tensor,1})`  
    ... converts a list of statistical tensors into the NORMALIZED alignment parameters A_kq = rho_kq / rho_00, which are
        the quantities an angular distribution or a polarization is usually expressed in.  The rho_00 component must be
        present and non-zero, since an empty ensemble has no alignment to speak of; an
        `aTensors::Array{Statistical.Tensor,1}` is returned, carrying the same ranks and the same axis.
"""
function alignmentParameters(tensors::Array{Statistical.Tensor,1})
    rho00 = ComplexF64(0.);   found = false
    for  tensor in tensors
        if  AngularMomentum.oneJ(tensor.k) == 0.  &&  AngularMomentum.oneM(tensor.q) == 0.   rho00 = tensor.value;   found = true   end
    end
    if  !found            error("Statistical.alignmentParameters(): no rho_00 component in the given list.")            end
    if  abs(rho00) == 0.  error("Statistical.alignmentParameters(): rho_00 = 0; the ensemble is empty, A_kq undefined.")  end

    aTensors = Statistical.Tensor[]
    for  tensor in tensors
        push!( aTensors, Tensor(tensor.k, tensor.q, tensor.a, tensor.b, tensor.axis, tensor.value / rho00) )
    end

    return( aTensors )
end


"""
`Statistical.computeTensors(kmax::Int64, amplitudes::Dict{AngularM64,ComplexF64}, a::LevelKey;`
                            `axis::Basics.AbstractQuantizationAxis=Basics.DefaultQuantizationAxis())`  
    ... computes all statistical tensors rho_kq of ranks k = 0..kmax for an ensemble described by COHERENT amplitudes
        a_M, one per magnetic sublevel, i.e. for the density matrix rho_(M,M') = a_M conj(a_M').  This is the case of a
        single, fully specified initial state; where the initial sublevels are unobserved and must be summed over
        incoherently, use the population form below instead.  A `tensors::Array{Statistical.Tensor,1}` is returned.
"""
function computeTensors(kmax::Int64, amplitudes::Dict{AngularM64,ComplexF64}, a::LevelKey;
                        axis::Basics.AbstractQuantizationAxis=Basics.DefaultQuantizationAxis())
    J   = a.sym.J;    Jx = AngularMomentum.oneJ(J);    MList = AngularMomentum.m_values(J)
    for  M in MList   if !haskey(amplitudes, M)   error("Statistical.computeTensors(): no amplitude given for M = $M.")   end   end
    tensors = Statistical.Tensor[]
    for  k = 0:kmax
        for  q = -k:k
            wa = ComplexF64(0.)
            for  M in MList
                for  Mp in MList
                    if  AngularMomentum.oneM(M) - AngularMomentum.oneM(Mp) != 1.0*q    continue    end
                    Mx = AngularMomentum.oneM(M);    Mpx = AngularMomentum.oneM(Mp)
                    wa = wa + (-1)^Int64(round(Jx - Mpx)) * AngularMomentum.ClebschGordan(Jx, Mx, Jx, -Mpx, 1.0*k, 1.0*q) *
                              amplitudes[M] * conj(amplitudes[Mp])
                end
            end
            if  abs(wa) > 0.   push!( tensors, Tensor(k, q, a, wa; axis=axis) )   end
        end
    end

    return( tensors )
end


"""
`Statistical.computeTensors(kmax::Int64, densityMatrix::Dict{Tuple{AngularM64,AngularM64},ComplexF64}, a::LevelKey;`
                            `axis::Basics.AbstractQuantizationAxis=Basics.DefaultQuantizationAxis())`  
    ... computes the statistical tensors rho_kq of ranks k = 0..kmax from the GENERAL density matrix rho(M,M') of the
        ensemble in the given level, which is the definition itself and of which the population and amplitude forms
        below are special cases.  A pair (M,M') that is absent from the dictionary is taken as zero, so only the
        non-vanishing elements need be supplied.  A `tensors::Array{Statistical.Tensor,1}` is returned.
"""
function computeTensors(kmax::Int64, densityMatrix::Dict{Tuple{AngularM64,AngularM64},ComplexF64}, a::LevelKey;
                        axis::Basics.AbstractQuantizationAxis=Basics.DefaultQuantizationAxis())
    J   = a.sym.J;    Jx = AngularMomentum.oneJ(J);    MList = AngularMomentum.m_values(J)
    tensors = Statistical.Tensor[]
    for  k = 0:kmax
        for  q = -k:k
            wa = ComplexF64(0.)
            for  M in MList
                for  Mp in MList
                    if  AngularMomentum.oneM(M) - AngularMomentum.oneM(Mp) != 1.0*q          continue    end
                    if  !haskey(densityMatrix, (M,Mp))                                       continue    end
                    Mx = AngularMomentum.oneM(M);    Mpx = AngularMomentum.oneM(Mp)
                    wa = wa + (-1)^Int64(round(Jx - Mpx)) * AngularMomentum.ClebschGordan(Jx, Mx, Jx, -Mpx, 1.0*k, 1.0*q) *
                              densityMatrix[(M,Mp)]
                end
            end
            if  abs(wa) > 0.   push!( tensors, Tensor(k, q, a, wa; axis=axis) )   end
        end
    end

    return( tensors )
end


"""
`Statistical.computeTensors(kmax::Int64, populations::Dict{AngularM64,Float64}, a::LevelKey;`
                            `axis::Basics.AbstractQuantizationAxis=Basics.DefaultQuantizationAxis())`  
    ... computes the statistical tensors rho_kq of ranks k = 0..kmax for an ensemble given by its sublevel POPULATIONS
        sigma(M), i.e. for a density matrix that is diagonal in M.  This is the ordinary case in which the initial
        sublevels are unobserved and summed over incoherently, and it yields q = 0 components only, since coherence
        between different M is what a non-zero q describes.  A `tensors::Array{Statistical.Tensor,1}` is returned.
"""
function computeTensors(kmax::Int64, populations::Dict{AngularM64,Float64}, a::LevelKey;
                        axis::Basics.AbstractQuantizationAxis=Basics.DefaultQuantizationAxis())
    J   = a.sym.J;    Jx = AngularMomentum.oneJ(J);    MList = AngularMomentum.m_values(J)
    for  M in MList   if !haskey(populations, M)   error("Statistical.computeTensors(): no population given for M = $M.")   end   end
    tensors = Statistical.Tensor[]
    for  k = 0:kmax
        wa = ComplexF64(0.)
        for  M in MList
            Mx = AngularMomentum.oneM(M)
            wa = wa + (-1)^Int64(round(Jx - Mx)) * AngularMomentum.ClebschGordan(Jx, Mx, Jx, -Mx, 1.0*k, 0.) * populations[M]
        end
        if  abs(wa) > 0.   push!( tensors, Tensor(k, 0, a, wa; axis=axis) )   end
    end

    return( tensors )
end


"""
`Statistical.densityMatrix(tensors::Array{Statistical.Tensor,1})`  
    ... reconstructs the sublevel density matrix rho(M,M') from a COMPLETE set of statistical tensors, i.e. inverts
        `computeTensors`.  The Clebsch-Gordan coefficients connecting the two form an orthogonal matrix, so the inverse
        is the same sum read the other way; it is exact rather than approximate, and a round trip through both returns
        the original matrix.  All ranks k = 0..2J must be present for the result to be complete, since a truncated set
        describes a different ensemble.  A `dm::Dict{Tuple{AngularM64,AngularM64},ComplexF64}` is returned.
"""
function densityMatrix(tensors::Array{Statistical.Tensor,1})
    if  length(tensors) == 0    error("Statistical.densityMatrix(): the given list of tensors is empty.")    end
    J   = tensors[1].a.sym.J;    Jx = AngularMomentum.oneJ(J);    MList = AngularMomentum.m_values(J)
    kmax = Int64(round(2*Jx))
    dm  = Dict{Tuple{AngularM64,AngularM64},ComplexF64}()
    for  M in MList
        for  Mp in MList
            Mx = AngularMomentum.oneM(M);    Mpx = AngularMomentum.oneM(Mp);    q = Int64(round(Mx - Mpx))
            wa = ComplexF64(0.)
            for  k = abs(q):kmax
                value = tensorValue(k, q, tensors; withZeros=true);    if  abs(value) == 0.    continue    end
                wa = wa + (-1)^Int64(round(Jx - Mpx)) * AngularMomentum.ClebschGordan(Jx, Mx, Jx, -Mpx, 1.0*k, 1.0*q) * value
            end
            if  abs(wa) > 0.    dm[(M,Mp)] = wa    end
        end
    end

    return( dm )
end


"""
`Statistical.invariant(tensors::Array{Statistical.Tensor,1}, k::Int64)`  
    ... returns the ROTATIONALLY INVARIANT magnitude sqrt( SUM_q |rho_kq|^2 ) of the rank k, which does not depend on the
        axis the tensors were computed about.  It is the honest quantity to quote when the axis is not the system's own
        symmetry axis: a single component such as rho_20 may be near zero merely because the axis is badly chosen, while
        this magnitude is not.  A `value::Float64` is returned.
"""
function invariant(tensors::Array{Statistical.Tensor,1}, k::Int64)
    wa = 0.
    for  tensor in tensors
        if  AngularMomentum.oneJ(tensor.k) == 1.0*k    wa = wa + abs(tensor.value)^2    end
    end

    return( sqrt(wa) )
end


"""
`Statistical.rotate(tensors::Array{Statistical.Tensor,1}, alpha::Float64, beta::Float64, gamma::Float64;`
                    `axis::Basics.AbstractQuantizationAxis=Basics.DefaultQuantizationAxis())`  
    ... rotates the given statistical tensors into a new frame obtained by the Euler angles (alpha, beta, gamma), using
        rho_kq(new) = SUM_q' rho_kq'(old) * conj( D^k_(q'q)(alpha,beta,gamma) ).  Each rank mixes only among its own
        2k+1 components, which is the property that makes the tensors worth forming at all, and the rank-0 component is
        left untouched.  The new axis is recorded in the returned components, so that a rotated tensor cannot later be
        mistaken for one about the original axis; a `tensors::Array{Statistical.Tensor,1}` is returned.
"""
function rotate(tensors::Array{Statistical.Tensor,1}, alpha::Float64, beta::Float64, gamma::Float64;
                axis::Basics.AbstractQuantizationAxis=Basics.DefaultQuantizationAxis())
    kList = Int64[]
    for  tensor in tensors   k = Int64(AngularMomentum.oneJ(tensor.k));   if  !(k in kList)   push!(kList, k)   end   end
    newTensors = Statistical.Tensor[]
    for  k in kList
        aKey = LevelKey();   bKey = LevelKey()
        for  tensor in tensors   if  Int64(AngularMomentum.oneJ(tensor.k)) == k   aKey = tensor.a;   bKey = tensor.b   end   end
        for  q = -k:k
            wa = ComplexF64(0.)
            for  qp = -k:k
                value = ComplexF64(0.)
                for  tensor in tensors
                    if  Int64(AngularMomentum.oneJ(tensor.k)) == k  &&  AngularMomentum.oneM(tensor.q) == 1.0*qp    value = tensor.value    end
                end
                if  abs(value) == 0.   continue   end
                wa = wa + value * conj( AngularMomentum.Wigner_DFunction(k, qp, q, alpha, beta, gamma) )
            end
            push!( newTensors, Tensor( AngularJ64(k), AngularM64(q), aKey, bKey, axis, wa) )
        end
    end

    return( newTensors )
end


"""
`Statistical.tensorValue(k::Int64, q::Int64, tensors::Array{Statistical.Tensor,1}; withZeros::Bool=false)`  
    ... returns the value of the statistical tensor component rho_kq if it is contained in the given list.  With
        withZeros = true a missing component is reported as zero, which is the right reading for a component that was
        simply not populated; otherwise a missing component raises, so that a typo in k or q is not silently read as a
        vanishing alignment.  A `value::ComplexF64` is returned.
"""
function tensorValue(k::Int64, q::Int64, tensors::Array{Statistical.Tensor,1}; withZeros::Bool=false)
    for  tensor in tensors
        if  AngularMomentum.oneJ(tensor.k) == 1.0*k  &&  AngularMomentum.oneM(tensor.q) == 1.0*q    return( tensor.value )    end
    end

    if   withZeros   return( ComplexF64(0.) )   else   error("Statistical: no tensor component found for k = $k and q = $q.")   end
end

end # module
