
"""
`module  JAC.SpinAngularTables`  
    ... a submodule of JAC that holds Gaigalas's quasispin data -- the coefficients of fractional parentage, the
        completely reduced W tensors and the q-space term list -- together with the small quantum-number helpers
        that index them, and nothing else.

        IT EXISTS SO THAT TWO MODULES CAN SHARE ONE COPY OF THE DATA. `SpinAngular` and `SpinAngular` implement
        the spin-angular algebra independently, but both read these tables; until 28-Aug-2026 the tables lived
        inside `SpinAngular` and `SpinAngular` reached across for them, which made the newer module depend on the
        one it is meant to replace. The algebra is deliberately NOT here: this module is data and indexing.

        THE DATA REACHES j <= 9/2 AND NO FURTHER, which is a property of Gaigalas's tables rather than of this code.
        A caller that needs more must either use a closed form -- an empty or singly occupied subshell needs no
        coefficient of fractional parentage, see `SpinAngular.shellReducedA` and `shellReducedW` -- or say so and
        stop. Returning zero beyond the tables is what hid a 27 % error in the electron-impact excitation cross
        section until 27-Aug-2026.
"""
module SpinAngularTables

using  Printf, ..AngularMomentum, ..Basics, ..ManyElectron

export  QspaceTerm, qshellTermM, qshellTermQ, qspacedelta,
        completlyReducedCfpByIndices, completelyReducedWkk, getTermNumber, qspaceTerms

"""
`struct  SpinAngularTables.QspaceTerm`  
    ... a struct for defining a subshell term/state  |j (nu) alpha Q J> == |j (nu) Q J Nr> for a subshell with well-defined j.

    + j        ::AngularJ64   ... subshell j
    + Q        ::AngularJ64   ... quasi-spin
    + J        ::AngularJ64   ... total J of subshell term
    + Nr       ::Int64        ... Additional quantum number Nr = 0,1,2.
    + min_odd  ::Int64        ... the minimal limits of the subshell terms for odd number operators in second quantization
    + max_odd  ::Int64        ... the maximal limits of the subshell terms for odd number operators in second quantization
    + min_even ::Int64        ... the minimal limits of the subshell terms for even number operators in second quantization
    + max_even ::Int64        ... the maximal limits of the subshell terms for even number operators in second quantization
"""
struct  QspaceTerm
    j          ::AngularJ64
    Q          ::AngularJ64
    J          ::AngularJ64
    Nr         ::Int64
    min_odd    ::Int64
    max_odd    ::Int64
    min_even   ::Int64
    max_even   ::Int64
end


"""
`SpinAngularTables.qshellTermM(j::AngularJ64, N::Int64)`  
    ... computes MQ quantum number; an M::Int64 is returned.

        j - is angular momentum j for the subshell;
        N - is number of electrons in the subshell;
"""
function  qshellTermM(j::AngularJ64, N::Int64)
    M = Int64(N-0.5*(Basics.twice(j)+1));  return( AngularM64(M//2) )
end


"""
`SpinAngularTables.qshellTermQ(j::AngularJ64, nu::Int64)`  
    ... computes Q quantum number; a Q::Int64 is returned.

        j  - is angular momentum j for the subshell;
        nu - is seniority for the subshell;
"""
function  qshellTermQ(j::AngularJ64, nu::Int64)
    Q = Int64((Basics.twice(j)+1)*0.5-nu);  return( AngularJ64(Q//2) )
end


"""
`SpinAngularTables.qspacedelta(q::AngularJ64, mq::AngularM64)`  
    ... computes trivial delta factors for Q space; a value::Int64 = {0,1} is returned.

        Q  - is the subshell total quasispin Q;
        MQ - is the subshell total projection of quasispin MQ;
"""
function  qspacedelta(q::AngularJ64, mq::AngularM64)
    if Basics.twice(q) < abs(Basics.twice(mq))  return( 0 ) end
    if (-1)^Int64(Basics.twice(q) + Basics.twice(mq)) == -1 return( 0 )  end
    return( 1 )
end

include("module-SpinAngularTables-inc-data.jl")

end # module
