
#
println("Ak) Apply & test the computation of spin-angular coefficients.")

configs = [Configuration("1s^2 2s"), Configuration("1s^2 2p")]

## configs = [Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 4f^4 6s^1") ] ## , 
           ## Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 4f^4 5d^1"), 
           ## Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 4f^3 5d^1 6p^1"), 
           ## Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 4f^3 6s^1 6p^1"), 
           ## Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 4f^4 6d^1"), 
           ## Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 4f^3 5d^1 7p^1"), 
           ## Configuration("1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 4f^4 7s^1")]
        
# Generate a list of relativistic configurations and determine an ordered list of subshells for these configurations
relconfList = ConfigurationR[]
for  conf in configs    wax = Basics.generate("configuration list: relativistic", conf);     append!( relconfList, wax)   end
for  i = 1:length(relconfList)    println(">> include ", relconfList[i])    end
subshellList = Basics.generate("subshells: ordered list for relativistic configurations", relconfList)
Defaults.setDefaults("relativistic subshell list", subshellList; printout=false)

# Generate the relativistic CSF's for the given subshell list
csfList = CsfR[]
for  relconf in relconfList
    newCsfs = Basics.generate("CSF list: from single ConfigurationR", relconf, subshellList);   append!( csfList, newCsfs)
end

sumSA = 0.;     sumFortran = 0.

# using Profile
# Profile.clear()
# @profile if  false
@time if false
    # Last successful:  unknown ...
    # Calculate angular coefficients for a scalar one- or two-particle operator
    for  leftCsf in csfList
        for rightCsf in csfList
            opa    = SpinAngular.OneParticleOperator(0, plus, true)     ## SpinAngular.TwoParticleOperator(0, plus, true)
            coeffs = SpinAngular.computeCoefficients(opa, leftCsf, rightCsf, subshellList)
            coefft = Basics.compute("angular coefficients: e-e, Ratip2013", leftCsf, rightCsf)
            ## println(">> Angular coeffients for <$leftCsf || O^($(op.rank))  || $rightCsf> = \n $coeffs ")
            if  length(coeffs)    != 0     println("\n>> Angular coeffients from SpinAngular = $coeffs ")         end
            if  length(coefft[1]) != 0     println(  ">> Angular coeffients from Ratip2013   = $(coefft[1]) ")    end
            for coeff in coeffs       global sumSA      = sumSA      + coeff.v     end
            for coeff in coefft[1]    global sumFortran = sumFortran + coeff.T     end
        end
    end
    @show sumSA, sumFortran
    #
elseif  false
    # Last successful:  unknown ...
    # Compute 
    # Calculate angular coefficients for a nonscalar one- particle operator
    for  leftCsf in csfList
        for rightCsf in csfList
            opb    = SpinAngular.OneParticleOperator(1, plus, true)  ## SpinAngular.OneParticleOperator(1, minus, true)
            coeffs = SpinAngular.computeCoefficients(opb, leftCsf, rightCsf, subshellList)
            ## coefft = Basics.compute("angular coefficients: 1-p, Ratip2013", 1, leftCsf, rightCsf)
            coefft = Basics.compute("angular coefficients: 1-p, Grasp92", 0, 1, leftCsf, rightCsf)
            ## println(">> Angular coeffients for <$leftCsf || O^($(op.rank))  || $rightCsf> = \n $coeffs ")
            if  length(coeffs) != 0     println("\n>> Angular coeffients from SpinAngular = $coeffs ")    end
            if  length(coefft) != 0     println(  ">> Angular coeffients from MCT/Grasp92 = $coefft ")    end
            for coeff in coeffs       global sumSA      = sumSA      + coeff.v     end
            for coeff in coefft       global sumFortran = sumFortran + coeff.T     end
        end
    end
    @show sumSA, sumFortran
    #
elseif  false
    # Last successful:  unknown ...
    # Compute 
    # Calculate angular coefficients for a nonscalar one- particle operator
    for  leftCsf in csfList
        for rightCsf in csfList
            opb    = SpinAngular.OneParticleOperator(2, plus, true)  ## SpinAngular.OneParticleOperator(2, plus, true)
            coeffs = SpinAngular.computeCoefficients(opb, leftCsf, rightCsf, subshellList)
            ## coefft = Basics.compute("angular coefficients: 1-p, Ratip2013", 2, leftCsf, rightCsf)
            coefft = Basics.compute("angular coefficients: 1-p, Grasp92", 0, 2, leftCsf, rightCsf)
            ## println(">> Angular coeffients for <$leftCsf || O^($(op.rank))  || $rightCsf> = \n $coeffs ")
            if  length(coeffs) != 0     println("\n>> Angular coeffients from SpinAngular = $coeffs ")    end
            if  length(coefft) != 0     println(  ">> Angular coeffients from MCT/Grasp92 = $coefft ")    end
            for coeff in coeffs       global sumSA      = sumSA      + coeff.v     end
            for coeff in coefft       global sumFortran = sumFortran + coeff.T     end
        end
    end
    @show sumSA, sumFortran
    #
elseif  true
    # Last successful:  unknown ...
    # Compute 
    # Calculate angular coefficients for a scalar two- particle operator
    for  leftCsf in csfList
        for rightCsf in csfList
            opb    = SpinAngular.TwoParticleOperator(0, plus, true)  ## SpinAngular.TwoParticleOperator(0, plus, true)
            coeffs = SpinAngular.computeCoefficients(opb, leftCsf, rightCsf, subshellList)
            coefft = Basics.compute("angular coefficients: e-e, Ratip2013", leftCsf, rightCsf)
            ## println(">> Angular coeffients for <$leftCsf || O^($(op.rank))  || $rightCsf> = \n $coeffs ")
            if  length(coeffs)    != 0     println("\n>> Angular coeffients from SpinAngular = $coeffs ")    end
            if  length(coefft[2]) != 0     println(  ">> Angular coeffients from Ratip2013   = $(coefft[2]) ")    end
            for coeff in coeffs       global sumSA      = sumSA      + coeff.v     end
            for coeff in coefft[2]    global sumFortran = sumFortran + coeff.V     end
        end
    end
    @show sumSA, sumFortran
    #
elseif  true
    # Last successful:  unknown ...
    # Compute 
    N1 = N2 = 0
    # Calculate angular coefficients for a scalar one- or two-particle operator
    op = SpinAngular.TwoParticleOperator(0, plus, true)
    for  leftCsf in csfList
        for rightCsf in csfList
            coeffs = SpinAngular.computeCoefficients(op, leftCsf, rightCsf, subshellList)
            global N2 = N2 + length(coeffs)
        end
    end
    #
    # Calculate angular coefficients for a nonscalar one- particle operator
    op = SpinAngular.OneParticleOperator(0, plus, true)
    for  leftCsf in csfList
        for rightCsf in csfList
            coeffs = SpinAngular.computeCoefficients(op, leftCsf, rightCsf, subshellList)
            global N1 = N1 + length(coeffs)
        end
    end
    @show length(csfList), N1, N2
    println("The Nd II test (Julia) has been completed")
    #
end

# Profile.print(maxdepth=30, mincount=1000)
nothing
