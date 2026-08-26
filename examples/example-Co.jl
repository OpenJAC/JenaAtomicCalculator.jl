
println("Co) Statistical tensors of an atomic ensemble: alignment and orientation, and the axis they are")
println("    stated about.")

using Printf

# WRITTEN 26-Aug-2026, when module-Statistical.jl was revived.  The module had been written and then commented out of
# JenaAtomicCalculator.jl, so `Statistical.Tensor` had never compiled, while `Basics.TensorComp` held the same object
# live -- two representations of rho_kq, one of them dead, and a dozen modules computing alignment through neither.
#
# WHAT THE OBJECT IS, in plain words.  An ensemble of atoms in one level is not described by "how many atoms" alone.
# If the atoms were made by something with a DIRECTION -- a light beam, an ion beam, a collision -- the magnetic
# sublevels M are populated unevenly, and the statistical tensor rho_kq is a repackaging of those populations:
#   rho_00  the total population, the number of atoms;
#   rho_1q  the ORIENTATION, a net direction of spin, which needs circularly polarized light or a polarized beam;
#   rho_2q  the ALIGNMENT, whether the ensemble is cigar- or pancake-shaped about the axis -- the same number of
#           atoms either way, differently arranged.
#
# WHY REPACKAGE.  Turn the axis and the sublevel populations mix into one another messily, whereas a tensor of rank k
# mixes ONLY among its own 2k+1 partners, through a single Wigner D-matrix.  That is the whole reason the object
# exists, and it is why the rotation belongs with it rather than beside any one process.
#
# THE AXIS IS PART OF THE QUANTITY, which branch b is about.  A_20 of a system whose symmetry axis lies elsewhere is a
# TRUE statement about z and a FALSE statement about the system -- and it comes out SMALL rather than wrong-looking,
# which is the dangerous way to be wrong.


if  true
    # Last visit:      26-Aug-2026
    # Last successful: 26-Aug-2026
    #
    # Branch a: ALIGNMENT IS GEOMETRY, and this branch is dated on a prediction that could have failed.
    #
    #   For J_i = 0 -> J_f = 1 by a single E1 multipole, the normalized alignment and orientation are fixed by the
    #   angular momenta ALONE.  No radial integral enters them, so every atom, every gauge and every line of that
    #   symmetry must return the SAME A_kq -- and the analytic values are
    #
    #       A_20 = 1/sqrt(2)  = 0.707107          for unpolarized, circular OR linear light
    #       A_10 = +-sqrt(3/2) = +-1.224745       for circular light, the sign following the helicity
    #       A_2,+-2 = -sqrt(3)/2 = -0.866025      for light linearly polarized along x
    #
    #   The test is that a REAL computation -- self-consistent field, relativistic amplitudes, the whole chain --
    #   has to reproduce numbers it knows nothing about.
    #
    # REPORT (26-Aug-2026): He-like calcium, 1s^2 -> 1s 2p, which gives two J=1(-) levels: the weak
    #   intercombination-like line 1->2 and the strong resonance line 1->4.  Their statistical tensors differ by a
    #   factor of 42 --
    #        line 1->2   rho_20 = 9.0473e-06        line 1->4   rho_20 = 3.8170e-04
    #   as they must, since rho_kq carries the strength of the transition.  Their ALIGNMENT PARAMETERS are identical:
    #        both lines, both gauges:   A_20    = 7.0711e-01     against 1/sqrt(2)   = 0.707107
    #        both lines, both gauges:   A_10    = 1.2247e+00     against sqrt(3/2)   = 1.224745   [circular]
    #        both lines, both gauges:   A_2,+-2 = -8.6603e-01    against -sqrt(3)/2  = -0.866025  [linear]
    #   Five figures, on two transitions whose oscillator strengths differ by forty, through Coulomb and Babushkin
    #   alike.  What that constrains is the whole chain at once, since the geometry cannot come out right by accident
    #   if the amplitudes, the photon density matrix or the tensor definition is wrong.
    #
    #   THIS TABLE COULD NOT BE PRINTED BEFORE 25-Aug-2026.  PhotoExcitation.computeStatisticalTensor returned a
    #   hard-wired -3.0 without touching the amplitudes, so rho_kq was -3.0 throughout and every A_kq was exactly
    #   1.000000 -- A_kq = tc/tc00 with both from the same constant.  Only the header said "under development".
    setDefaults("print summary: open", "zzz-Statistical-Co-photoexcitation.sum")

    for  (name, stokes)  in  [("unpolarized", ExpStokes(0., 0., 0.)), ("circular   P3=+1", ExpStokes(0., 0., 1.)),
                              ("linear     P1=+1", ExpStokes(1., 0., 0.))]
        println("\n\n########## incident light: $name ##########")
        pSettings = PhotoExcitation.Settings(PhotoExcitation.Settings(), multipoles=[E1],
                                             gauges=[UseCoulomb, UseBabushkin], calcTensors=true, stokes=stokes)
        wa = Atomic.Computation(Atomic.Computation(), name="Co-a", grid=Radial.Grid(true),
                                nuclearModel     = Nuclear.Model(20.),
                                initialConfigs   = [Configuration("1s^2")],
                                finalConfigs     = [Configuration("1s 2p")],
                                processSettings  = pSettings )
        wb = perform(wa)
    end
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      26-Aug-2026
    # Last successful: 26-Aug-2026
    #
    # Branch b: THE AXIS IS PART OF THE QUANTITY.  The same ensembles, stated about three different axes.
    #
    #   Branch a left two definite objects, and the difference between them is the point of this branch:
    #     + UNPOLARIZED light along z leaves an ensemble that is AXIALLY SYMMETRIC about the beam, carrying
    #       A_20 = 1/sqrt(2) and nothing else;
    #     + LINEARLY polarized light leaves one that is NOT, carrying A_2,+-2 = -sqrt(3)/2 as well.
    #   Nothing about the atoms changes when we choose a different axis to describe them about, but every individual
    #   component does.  One number is watched that must NOT move: the rotational invariant sqrt( SUM_q |A_2q|^2 ),
    #   which is what "how aligned is this ensemble" actually means.
    #
    # REPORT (26-Aug-2026):
    #        UNPOLARIZED light -- axially symmetric, invariant 0.707107 throughout
    #             about the beam z                A_20 = +7.071068e-01
    #             about the polarization x        A_20 = -3.535534e-01
    #             about the magic angle           A_20 = +3.925231e-17
    #        LINEAR light -- not axially symmetric, invariant 1.414214 throughout
    #             about the beam z                A_20 = +7.071068e-01
    #             about the polarization x        A_20 = -1.414214e+00
    #             about the magic angle           A_20 = -7.071068e-01
    #
    #   THREE THINGS ARE CHECKED, and each could fail.
    #
    #   FIRST, for an axially symmetric ensemble the rotation must follow A_20(beta) = A_20(0) P_2(cos beta) exactly,
    #   since a q = 0 component of rank 2 transforms by the Legendre polynomial and nothing else.  At beta = pi/2,
    #   P_2(0) = -1/2, and 0.707107 x (-1/2) = -0.353553 is what comes out.
    #
    #   SECOND, A_20 = -sqrt(2) = -1.414214 for the LINEAR case about the polarization axis is a KNOWN ANSWER rather
    #   than an output: linearly polarized light on J_i=0 -> J_f=1 leaves a PURE M = 0 state about the polarization
    #   direction, and a pure M=0 state of J=1 has A_20 = -sqrt(2).  Reproduced to 2e-16, which also PINS the Stokes
    #   convention -- P1 is linear polarization along x -- rather than leaving it assumed.
    #
    #   THIRD, AND THIS IS THE WARNING.  At the magic angle, beta = acos(1/sqrt(3)) = 54.7 degrees, the AXIALLY
    #   SYMMETRIC ensemble reports A_20 = 4e-17.  An ensemble as aligned as it was before reports NO ALIGNMENT,
    #   because P_2(cos beta) vanishes there and the single component was never the physical quantity; the invariant
    #   sits unmoved at 0.707107.  Note that the LINEAR ensemble does NOT vanish at the same angle -- its q = +-2
    #   components feed back into A_20 -- so the trap is not "the magic angle is dangerous" but "a single component
    #   is not a magnitude".
    #
    #   THIS IS NOT HYPOTHETICAL.  On 14-Aug-2026 an A_20 of 2.5e-17 at a C_3v crystal site was reported here as "no
    #   quadrupole interaction".  The site sat on the body diagonal of the cube -- the magic angle, near enough --
    #   and the rotational invariant was 1.4e-03.  A true statement about z, read as a false statement about the
    #   system.  `Statistical.invariant` exists for that reason, and the axis is carried in the type so that the
    #   question can be asked at all.
    setDefaults("print summary: open", "zzz-Statistical-Co-frames.sum")

    fKey = Basics.LevelKey( LevelSymmetry(AngularJ64(1), Basics.minus), 2, 0., 1.)
    unpolarized = [ Statistical.Tensor(0,  0, fKey, ComplexF64( 1.0)),
                    Statistical.Tensor(2,  0, fKey, ComplexF64( 1/sqrt(2))) ]
    linear      = [ Statistical.Tensor(0,  0, fKey, ComplexF64( 1.0)),
                    Statistical.Tensor(2, -2, fKey, ComplexF64(-sqrt(3)/2)),
                    Statistical.Tensor(2,  0, fKey, ComplexF64( 1/sqrt(2))),
                    Statistical.Tensor(2,  2, fKey, ComplexF64(-sqrt(3)/2)) ]
    for  (label, ts)  in  [("UNPOLARIZED light  -- axially symmetric about the beam", unpolarized),
                           ("LINEAR      light  -- NOT axially symmetric",            linear)]
        println("\n  $label\n")
        @printf("     %-30s %18s %14s\n", "described about", "A_20", "invariant")
        for  (name, alpha, beta)  in  [("the beam (z)", 0., 0.), ("the polarization (x)", 0., pi/2),
                                       ("the magic angle, 54.7 deg", 0., acos(1/sqrt(3)))]
            rs = Statistical.rotate(ts, alpha, beta, 0.)
            @printf("     %-30s %+18.6e %14.6f\n", name, real(Statistical.tensorValue(2, 0, rs; withZeros=true)),
                    Statistical.invariant(rs, 2))
        end
    end
    println("\n  A_20 moved with every choice of axis; the invariant moved with none.  Only one of the two is a")
    println("  property of the atoms.")
    setDefaults("print summary: close", "")
    #
elseif  false
    # Last visit:      26-Aug-2026
    # Last successful: 26-Aug-2026
    #
    # Branch c: THE MACHINERY ITSELF -- what rho_kq is, and what it is not allowed to do.
    #
    #   Four properties, each of which the code could fail:
    #     1. rho_kq and the sublevel density matrix rho(M,M') are the SAME information in two arrangements, so the
    #        round trip must return the original matrix exactly;
    #     2. a rank mixes ONLY among its own 2k+1 partners under rotation -- this is the entire reason for forming
    #        the tensors, and it is a statement the code can violate;
    #     3. an ISOTROPIC ensemble has no alignment and no orientation, whatever the rank;
    #     4. a component that was never populated must be REFUSED rather than returned as zero, since a typo in k or
    #        q would otherwise read as a physically vanishing alignment.
    #
    # REPORT (26-Aug-2026):
    #        round trip  dm -> rho_kq -> dm      J=1  2.2e-16     J=2  4.4e-16     J=3  8.9e-16
    #        rank leakage under rotation         J=2  0.0e+00  (a pure rank 2 stays pure rank 2)
    #        isotropic ensemble, A_k for k>0     J=2  1.2e-17
    #        a missing component                 raises, and names k and q
    setDefaults("print summary: open", "zzz-Statistical-Co-machinery.sum")

    println("\n  1. round trip:  rho(M,M')  ->  rho_kq  ->  rho(M,M')\n")
    for  J  in  [1, 2, 3]
        key = Basics.LevelKey( LevelSymmetry(AngularJ64(J), Basics.plus), 1, 0., 1.)
        Ms  = AngularMomentum.m_values(key.sym.J)
        amp = Dict{AngularM64,ComplexF64}( Ms[i] => ComplexF64(0.3*i, 0.15*(J-i))  for i = 1:length(Ms) )
        dm0 = Dict{Tuple{AngularM64,AngularM64},ComplexF64}( (M,Mp) => amp[M]*conj(amp[Mp]) for M in Ms, Mp in Ms )
        dm1 = Statistical.densityMatrix( Statistical.computeTensors(2*J, dm0, key) )
        @printf("     J=%d   max |difference| = %.3e\n", J, maximum(abs(dm1[(M,Mp)] - dm0[(M,Mp)]) for M in Ms, Mp in Ms))
    end

    println("\n  2. a rank mixes only within itself under rotation\n")
    key = Basics.LevelKey( LevelSymmetry(AngularJ64(2), Basics.plus), 1, 0., 1.)
    pure2 = [ Statistical.Tensor(2, q, key, ComplexF64(0.4*q + 0.2, 0.1*q)) for q = -2:2 ]
    rot2  = Statistical.rotate(pure2, 0.6, 1.2, 0.4)
    leak  = maximum( [0.0; [abs(t.value) for t in rot2 if AngularMomentum.oneJ(t.k) != 2.]] )
    @printf("     a pure rank-2 tensor rotated:  largest component outside rank 2 = %.3e\n", leak)

    println("\n  3. an isotropic ensemble has no alignment and no orientation\n")
    Ms  = AngularMomentum.m_values(key.sym.J)
    iso = Statistical.alignmentParameters( Statistical.computeTensors(4, Dict{AngularM64,Float64}(M => 1.0 for M in Ms), key) )
    waMax = maximum( [0.0; [abs(t.value) for t in iso if AngularMomentum.oneJ(t.k) > 0.]] )
    @printf("     largest |A_kq| with k > 0 = %.3e\n", waMax)

    println("\n  4. a component that was never populated is refused, not returned as zero\n")
    try     Statistical.tensorValue(9, 0, iso);   println("     NOT GUARDED")
    catch e println("     ", first(split(sprint(showerror, e), "\n")))
    end
    setDefaults("print summary: close", "")
    #
end
