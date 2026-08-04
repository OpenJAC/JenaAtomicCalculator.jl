#
println("Al) Apply & test for parallel computating methods/techniques with Julia.")
using Distributed
@everywhere using LinearAlgebra

println("Number of cores = $(nprocs()),  number of workers = $(nworkers())")

@everywhere function inside(n)  
    ni = 0;     for  i in 1:n   x,y = rand(), rand();   if  x^2 + y^2 <= 1      ni = ni + 1     end     end
    return( ni )
end
#
function piParallel(n)  
    p = nworkers();     ni = @distributed (+) for i = 1:p   inside(n/p)     end
    return( 4*ni/n )
end
#
@everywhere function myAction(A,B)
    C = A*B;    return( C[1,1] )
end

nx = 5000
# using Profile
# Profile.clear()
# @profile if  false
@time if false
# if false
    # Last successful:  unknown ...
    # Sequential processing
    for  i in 1:10
        wa = rand(nx,nx)
        wb = LinearAlgebra.eigvals(wa)
        println("wb = $(wb[1:3])")
    end

elseif true
    # Last successful:  unknown ...
    # Broadcasting of values and functions
    A = rand(nx,nx);    B = rand(nx,nx);    Cf = Future[]
    @everywhere A, B = $A, $B
    for  i in 1:15     ii = rem(i,nworkers())+1;     push!(Cf, @spawnat ii (A.^17.3)[1,1])       end
    ## for  i in 1:15     rf = remotecall(*, WorkerPool(workers()), A,B);   push!(Cf, rf)   end
    println("first loop done")
    for  i in 1:15     C = fetch(Cf[i]);   println("$i  $(C[1,1])")    end

elseif false
    # Last successful:  unknown ...
    # Producer-consumer model
    @async foo()

elseif true
    # Last successful:  unknown ...
    # Distributed calculations
    nx = 1000000000
    @time wc = inside(nx);         println("wc = $wc")   
    @time wd = piParallel(nx);     println("wd = $wd")   

elseif false

end
