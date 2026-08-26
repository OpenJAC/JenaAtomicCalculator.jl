
# TWISTED (Bessel) photon beams, in the form-factor approximation, for an atom at FINITE IMPACT PARAMETER from the beam axis.
#
# WHY THIS IS TRACTABLE HERE AND NOWHERE ELSE IN THE MODULE. A Bessel beam is a coherent superposition of plane waves whose
# wave vectors lie on a cone of half-angle theta_k, each carrying the azimuthal phase exp(i m phi_k). For a general scattering
# amplitude that superposition requires the partial-wave/OAM coupling that ParticleScattering explicitly declined to re-derive
# (see module-ParticleScattering-inc-beams.jl, where the Bessel method raises). In the FORM-FACTOR approximation the amplitude
# depends on NOTHING but the momentum transfer q, so each plane-wave component's amplitude is already known and the twisted
# case collapses to a one-dimensional integral over phi_k. That is the whole reason this file can exist while the second-order
# twisted case cannot.
#
#     A_m(b)  =  (1/2pi) Int_0^2pi  d phi_k  exp(i m phi_k)  exp(i kappa b cos(phi_k - phi_b))  F(q(theta_k, phi_k; theta, phi))
#
# with kappa = k sin(theta_k) the transverse momentum and b the impact parameter. The geometry is exact:
#
#     cos(Theta) = sin(theta_k) sin(theta) cos(phi_k - phi) + cos(theta_k) cos(theta) ,      q = 2 k sin(Theta/2) .
#
# COHERENT AND INCOHERENT SCATTERING BEHAVE DIFFERENTLY HERE, and the difference is physical rather than technical.
#   * RAYLEIGH is coherent: the amplitudes from the different cone components ADD before being squared, so the OAM phase and
#     the impact-parameter phase survive and produce the vortex structure. That is what this file computes.
#   * COMPTON, in the closure approximation used by this module, is not an amplitude at all -- S(q) is already a
#     cross-section-level quantity summed over final states. The phases then cancel in the modulus, |exp(i m phi_k)| = 1, and
#     the twisted incoherent cross section degenerates into a plain cone-average of the plane-wave one, with NO m-dependence
#     and NO b-dependence whatever. THAT IS A LIMITATION OF THE CLOSURE, NOT A PHYSICAL PREDICTION: a proper treatment needs
#     the amplitude to each final state separately, which S(q) has already summed away. It is computed and reported, with the
#     degeneracy stated, rather than omitted -- but no OAM physics should be read out of it.
#
# WHAT IS NOT HERE. This is a SCALAR treatment: the vector nature of the Bessel field and its polarization structure are not
# included, so the (1 + cos^2 theta) factor is carried over unchanged from the plane-wave case. A full vector treatment would
# change the angular distribution and is the natural next step; nothing in this file should be quoted as a polarization result.
#
# TWO CHECKS THAT COST NOTHING AND CAN BOTH FAIL.
#   1. theta_k --> 0 must reproduce the PLANE-WAVE result, and example-Pd.jl branch a is DATED, so this is a limit check
#      against an anchored absolute number rather than against another unverified calculation.
#   2. AT b = 0 THE CROSS SECTION MUST VANISH FOR m /= 0. On the axis a vortex beam has a node: with F depending on phi_k only
#      weakly, the integral tends to i^m J_m(kappa b) times a constant, and J_m(0) = 0 for every m /= 0. A non-zero on-axis
#      result for m = 1 would mean the OAM phase is not being applied.


"""
`PhotonScattering.besselAmplitude(process::PhotonScattering.AbstractPhotonProcess, k::Float64, beam::Beam.BesselBeam,
                                  bImpact::Float64, theta::Float64, phi::Float64, level::Level, grid::Radial.Grid)`
    ... to compute the coherent Bessel-beam scattering amplitude of one atom at impact parameter `bImpact` from the beam axis,
        as the cone integral given in the header of this file. An amplitude::ComplexF64 is returned.

        The integral is a 721-point trapezoid over phi_k, which is ample: the integrand is a form factor varying smoothly with
        the momentum transfer, times two phases whose highest frequency is set by max(|m|, kappa b). The azimuthal angle of the
        impact parameter is taken as zero, which costs no generality for a spherically symmetric target -- it merely fixes the
        origin of phi.
"""
function besselAmplitude(process::PhotonScattering.AbstractPhotonProcess, k::Float64, beam::Beam.BesselBeam,
                         bImpact::Float64, theta::Float64, phi::Float64, level::Level, grid::Radial.Grid)
    thetaK = beam.openingAngle;    m = beam.mOAM;    kappa = k * sin(thetaK)
    nphi   = 721;    dphi = 2pi / (nphi - 1);    A = ComplexF64(0.)

    for  i = 1:nphi
        phik  = (i-1) * dphi
        w     = (i == 1 || i == nphi) ? 0.5 : 1.0
        cosTh = sin(thetaK) * sin(theta) * cos(phik - phi) + cos(thetaK) * cos(theta)
        cosTh = max(-1.0, min(1.0, cosTh))
        q     = 2 * k * sqrt( max(0., (1 - cosTh) / 2) )
        F     = FormFactor.standardF(q, level, grid)
        A     = A + w * exp(im * m * phik) * exp(im * kappa * bImpact * cos(phik)) * F * dphi
    end

    return( A / (2pi) )
end


"""
`PhotonScattering.computeTwistedLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model,
                                      grid::Radial.Grid, settings::PhotonScattering.Settings; output=true)`
    ... to compute the scattering of a twisted (Bessel) beam in the form-factor approximation, for every selected level, photon
        energy and impact parameter, and to display the result. An Array{PhotonScattering.Line,1} is returned if output=true,
        and nothing otherwise.

        One Line is produced per (level, photon energy, impact parameter); the impact parameter is carried in the Line's
        `particleEnergy` slot, which is otherwise unused by a photon-in/photon-out process. That is a deliberate reuse rather
        than a new field, and it is documented here because it would otherwise be a trap for the next reader.
"""
function computeTwistedLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                             settings::PhotonScattering.Settings; output=true)
    println("")
    printstyled("PhotonScattering.computeTwistedLines(): The twisted-beam computation starts now ... \n", color=:light_green)
    printstyled("------------------------------------------------------------------------------------ \n", color=:light_green)
    #
    beam = settings.beamType
    println(">> $(beam);  impact parameters $(settings.impactParameters) a.u.")
    if  settings.process == PhotonScattering.ComptonScattering()
        println(">> NOTE: in the closure approximation the INCOHERENT cross section carries no OAM and no impact-parameter " *
                "dependence -- S(q) is already summed over final states, so the phases cancel in the modulus. What follows " *
                "is a cone-average of the plane-wave result. This is a limitation of the closure, not a physical prediction.")
    end

    newLines = PhotonScattering.Line[]
    for  iLevel in initialMultiplet.levels
        for  fLevel in finalMultiplet.levels
            if  !Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)  ||  iLevel.index != fLevel.index    continue    end
            for  omega in settings.photonEnergies
                inOmega = Defaults.convertUnits("energy: to atomic", omega)
                for  bImpact in settings.impactParameters
                    obs, sigma = PhotonScattering.twistedObservables(inOmega, bImpact, iLevel, grid, settings)
                    push!( newLines, PhotonScattering.Line(iLevel, fLevel, inOmega, inOmega, bImpact,
                                                            EmProperty(sigma, sigma), PhotonScattering.Channel[], obs) )
                end
            end
        end
    end

    PhotonScattering.displayTwistedResults(stdout, newLines, grid, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   PhotonScattering.displayTwistedResults(iostream, newLines, grid, settings)     end

    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`PhotonScattering.displayTwistedResults(stream::IO, lines::Array{PhotonScattering.Line,1}, grid::Radial.Grid,
                                        settings::PhotonScattering.Settings)`
    ... to print the twisted-beam cross sections against impact parameter, together with the plane-wave value they must
        approach as the opening angle tends to zero. A neat table is printed but nothing is returned otherwise.
"""
function displayTwistedResults(stream::IO, lines::Array{PhotonScattering.Line,1}, grid::Radial.Grid,
                               settings::PhotonScattering.Settings)
    nx = 128
    alpha = Defaults.getDefaults("alpha")
    println(stream, " ")
    println(stream, "  Twisted-beam scattering cross sections   [process = $(settings.process), m = $(settings.beamType.mOAM), " *
                    "opening angle = $(settings.beamType.openingAngle) rad]:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  " * TableStrings.center(10, "Level"; na=2) * TableStrings.center(16, "Energy"; na=4) *
                TableStrings.center(16, "b [a.u.]"; na=4) * TableStrings.center(18, "Cross section"; na=4) *
                TableStrings.center(20, "sigma / N^2 sigma_T"; na=0)
    println(stream, sa)
    println(stream, "  ", TableStrings.hLine(nx))
    for  line in lines
        N       = FormFactor.standardF(0., line.initialLevel, grid)
        thomson = PhotonScattering.thomsonLimit(settings.process, N, alpha)
        sa  = "  " * TableStrings.center(10, TableStrings.level(line.initialLevel.index); na=2)
        sa  = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.inPhotonEnergy)) * "    "
        sa  = sa * @sprintf("%14.6f", line.particleEnergy) * "    "
        sa  = sa * @sprintf("%.6e", line.crossSection.Coulomb) * "        "
        sa  = sa * @sprintf("%.6f", line.crossSection.Coulomb / thomson)
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    println(stream, "  b is the impact parameter of the atom from the beam axis, carried in the Line's particleEnergy slot.")
    if  settings.beamType.mOAM != 0
        println(stream, "  m /= 0: the b = 0 row MUST vanish -- a vortex beam has a node on its axis, and J_m(0) = 0.")
    end
    println(stream, "  As the opening angle tends to zero every row must approach the PLANE-WAVE value of example-Pd.jl")
    println(stream, "  branch a, which is dated and anchored absolutely against N^2 sigma_Thomson.")

    return( nothing )
end


"""
`PhotonScattering.twistedObservables(inOmega::Float64, bImpact::Float64, level::Level, grid::Radial.Grid,
                                     settings::PhotonScattering.Settings)`
    ... to compute the angle-differential cross sections at the requested angles and the angle-integrated total for a twisted
        beam and one impact parameter. A tuple (observables::Array{PhotonScattering.Observables,1}, sigma::Float64) is returned.

        For RAYLEIGH the coherent cone amplitude of PhotonScattering.besselAmplitude is squared; for COMPTON the incoherent
        S(q) is cone-AVERAGED instead, the closure having already summed over final states so that no phase survives. The
        Stokes parameters are NaN throughout: this is a scalar treatment and carries no polarization information.
"""
function twistedObservables(inOmega::Float64, bImpact::Float64, level::Level, grid::Radial.Grid,
                            settings::PhotonScattering.Settings)
    alpha = Defaults.getDefaults("alpha");    wc = Defaults.getDefaults("speed of light: c")
    reSq  = alpha^4;    k = inOmega / wc;    beam = settings.beamType

    function dcsAt(theta::Float64, phi::Float64)
        if      settings.process == PhotonScattering.RayleighScattering()
            A = PhotonScattering.besselAmplitude(settings.process, k, beam, bImpact, theta, phi, level, grid)
            return( reSq / 2 * (1 + cos(theta)^2) * abs2(A) )
        else
            # incoherent: cone-average the cross-section-level S(q); no phase survives the modulus
            thetaK = beam.openingAngle;    nphi = 361;    dphi = 2pi / (nphi - 1);    acc = 0.
            for  i = 1:nphi
                phik  = (i-1) * dphi;    w = (i == 1 || i == nphi) ? 0.5 : 1.0
                cosTh = sin(thetaK) * sin(theta) * cos(phik - phi) + cos(thetaK) * cos(theta)
                cosTh = max(-1.0, min(1.0, cosTh))
                q     = 2 * k * sqrt( max(0., (1 - cosTh) / 2) )
                acc   = acc + w * PhotonScattering.incoherentScatteringFunction(q, level, grid) * dphi
            end
            return( reSq / 2 * (1 + cos(theta)^2) * acc / (2pi) )
        end
    end

    obs = PhotonScattering.Observables[]
    for  theta in settings.polarThetas
        phi = length(settings.polarPhis) > 0 ? settings.polarPhis[1] : 0.
        d   = dcsAt(theta, phi)
        push!( obs, PhotonScattering.Observables(theta, phi, EmProperty(d, d), NaN, NaN, NaN) )
    end

    # angle-integrated total, over both polar and azimuthal angles: a twisted beam breaks the azimuthal symmetry that the
    # plane-wave case enjoys, so phi cannot be integrated analytically here and is done on a grid like theta.
    nth = 121;   nph = 61;   dth = pi / (nth - 1);   dph = 2pi / (nph - 1);   sigma = 0.
    for  i = 1:nth
        theta = (i-1) * dth;    wt = (i == 1 || i == nth) ? 0.5 : 1.0
        for  j = 1:nph
            phi = (j-1) * dph;   wp = (j == 1 || j == nph) ? 0.5 : 1.0
            sigma = sigma + wt * wp * dcsAt(theta, phi) * sin(theta) * dth * dph
        end
    end

    return( (obs, sigma) )
end
