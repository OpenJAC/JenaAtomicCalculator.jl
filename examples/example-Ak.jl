
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
    # Last successful:  30-Aug-2026 -- the rank-0 DIAGONAL one-particle coefficients equal the occupation
    #   numbers exactly: 2.0000000000000004 for the doubly-occupied 1s_1/2 and 1.0 for the singly-occupied
    #   subshell of each of the three CSFs, largest departure 4.44e-16. That is an EXACT check needing no
    #   reference data -- under the GRASP convention the sqrt(2j+1) sits inside the coefficient, so the
    #   rank-0 diagonal coefficient IS the occupation (Rule 18), and a wrong convention would be visible
    #   by eye. Dated on the maintainer's decision of 30-Aug-2026.
    #   WHAT THIS DATE DOES NOT COVER, so nobody reads more into it: rank 0 only, on a three-CSF Li-like
    #   system. It says nothing about rank > 0, nothing about the VALUES of the 21 two-particle
    #   coefficients that N2 merely counts, and nothing about any off-diagonal coefficient.
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
    #
    # THE PHYSICS CHECK, and the reason this branch is worth running rather than merely completing.
    # N1 and N2 above are COUNTS: they would be unchanged if every coefficient were wrong by a factor, so they
    # test that the machinery runs and nothing more. The rank-0 DIAGONAL one-particle coefficient, by contrast,
    # has an exact expected value and needs no reference data: under the GRASP convention the sqrt(2j+1) sits
    # INSIDE the coefficient, so that coefficient IS the occupation number of its subshell (Rule 18, and the
    # reason the convention was chosen -- a wrong one is then visible by eye instead of invisible).
    println("\n>> Rank-0 diagonal one-particle coefficients against the occupation numbers.")
    println(">> Each T must EQUAL the occupation of its subshell, exactly.")
    opDiag = SpinAngularGaigalas.OneParticleOperator(0, plus, true)
    worst  = 0.
    for  (i, csf) in enumerate(csfList)
        occ = [csf.occupation[k] for k = 1:length(subshellList)]
        println("     CSF $i, occupations $occ:")
        for  c in SpinAngularGaigalas.computeCoefficients(opDiag, csf, csf, subshellList)
            k        = findfirst(isequal(c.a), subshellList)
            expected = occ[k]
            global worst = max(worst, abs(c.T - expected))
            println("        $(c.a):  T = $(c.T)   expected (occupation) = $expected")
        end
    end
    println(">> Largest departure from the occupation number: $worst   (must be 0 to machine precision)")
    #
    println("The Li-like test (1s^2 2s + 1s^2 2p, Julia route) has been completed")
    #
end

# Profile.print(maxdepth=30, mincount=1000)
nothing
