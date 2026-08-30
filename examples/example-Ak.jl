
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
# The string-dispatch form of Basics.generate is retired. The current route is the one src/ itself uses at
# module-BasicsAZ-inc-generate.jl:152-158, verbatim: RelativisticConfigurations() for the configurations and
# Basics.generateSubshellList for the ordered subshells.
for  conf in configs
    wax = Basics.generateConfigurations(Basics.RelativisticConfigurations(), conf);   append!( relconfList, wax)
end
for  i = 1:length(relconfList)    println(">> include ", relconfList[i])    end
subshellList = Basics.generateSubshellList(relconfList)
Defaults.setDefaults("relativistic subshell list", subshellList; printout=false)

# Generate the relativistic CSF's for the given subshell list
csfList = CsfR[]
for  relconf in relconfList
    newCsfs = Basics.generateCsfRs(relconf, subshellList);                             append!( csfList, newCsfs)
end

sumSA = 0.;     sumFortran = 0.

# using Profile
# Profile.clear()
# @profile if  false
# THE FOUR BRANCHES THAT STOOD HERE COMPARED JAC's Julia spin-angular coefficients against the old FORTRAN
# route (AngularCoefficients-Ratip2013, reached through Basics.compute(AngularCoeffs...)). That route was
# RETIRED on 30-Aug-2026: its module source was already absent from src/, its include already commented out
# of JenaAtomicCalculator.jl, and every call site in src/ sat behind Defaults.saRatip(), hardcoded false.
# The comparison therefore could not run and had nothing left to compare against. What remains below is the
# branch that exercises the Julia route on its own.
@time if  true
    # Last successful:  unknown ...
    # Compute 
    N1 = N2 = 0
    # Calculate angular coefficients for a scalar one- or two-particle operator
    op = SpinAngularGaigalas.TwoParticleOperator(0, plus, true)
    for  leftCsf in csfList
        for rightCsf in csfList
            coeffs = SpinAngularGaigalas.computeCoefficients(op, leftCsf, rightCsf, subshellList)
            global N2 = N2 + length(coeffs)
        end
    end
    #
    # Calculate angular coefficients for a nonscalar one- particle operator
    op = SpinAngularGaigalas.OneParticleOperator(0, plus, true)
    for  leftCsf in csfList
        for rightCsf in csfList
            coeffs = SpinAngularGaigalas.computeCoefficients(op, leftCsf, rightCsf, subshellList)
            global N1 = N1 + length(coeffs)
        end
    end
    @show length(csfList), N1, N2
    println("The Nd II test (Julia) has been completed")
    #
end

# Profile.print(maxdepth=30, mincount=1000)
nothing
