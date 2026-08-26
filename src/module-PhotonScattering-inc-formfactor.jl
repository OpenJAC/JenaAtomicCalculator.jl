
# The FORM-FACTOR approximation for coherent (Rayleigh) and incoherent (Compton) photon scattering.
#
# THE CHEAP NON-RELATIVISTIC LIMIT, and the point of having it. The atomic electrons are treated as free scatterers whose only
# collective property is the Fourier transform of their charge density -- the atomic form factor F(q), with q = 2 k sin(theta/2)
# the momentum transfer. There is NO sum over intermediate states, NO gauge and NO multipole expansion, which is why it costs
# almost nothing and why it cannot describe a resonance.
#
#     Rayleigh (coherent)    d sigma / d Omega  =  (r_e^2 / 2) (1 + cos^2 theta) |F(q)|^2
#     Compton  (incoherent)  d sigma / d Omega  =  (r_e^2 / 2) (1 + cos^2 theta) S(q)
#
# with r_e = alpha^2 in atomic units.
#
# IT EXISTS TO BE COMPARED AGAINST THE SECOND-ORDER SUM, and the comparison is worth having precisely because the two can
# DISAGREE informatively. Where they agree, the expensive machinery is confirmed by something that shares none of its
# machinery; where they disagree, the disagreement is the physics. A baseline that could only ever agree would be worth much
# less. The two are valid in OPPOSITE regimes: the second-order sum near and below the resonances, where a bound electron
# responds through the atom's polarizability and sigma ~ omega^4; the form factor well ABOVE the binding energies, where the
# electrons do respond freely and sigma tends to a constant.
#
# AND IT CARRIES THE ONLY ABSOLUTE CHECK IN THIS MODULE -- but the two processes approach DIFFERENT limits, from OPPOSITE
# directions in energy, and one number will not serve both:
#
#     RAYLEIGH   as q --> 0     F(q) --> N       sigma --> N^2 sigma_Thomson     electrons scatter IN PHASE, amplitudes add
#     COMPTON    as q --> large S(q) --> N       sigma --> N   sigma_Thomson     electrons resolved, CROSS SECTIONS add
#
# with sigma_Thomson = (8 pi/3) alpha^4. So the coherent ratio must tend to 1 as the photon energy FALLS and the incoherent one
# as it RISES. Both are CLOSED FORMS with no prefactor to derive and nothing taken from a table. Every other test in
# PhotonScattering -- the omega^4 law, the Z^5 scan, detailed balance, the RIXS Lorentzian, low-frequency additivity -- is a
# ratio or a shape, deliberately, because the cross-section prefactors are underived. These are not. If a computed limit misses
# its Thomson value the fault is in this file's own normalization and nowhere else.
#
# MEASURED 22-Aug-2026 on beryllium-like neon: the coherent ratio is 0.999990 at 1 eV, i.e. ONE PART IN 10^5, falling to 0.943
# at 2 keV where q a_0 = 1.07 and the form factor has genuinely left F(0). The incoherent ratio rises from 2e-6 at 10 eV to
# 0.057 at 2 keV, and its low-energy behaviour is sigma ~ omega^2 exactly, which follows from S(q) ~ N q^2 <r^2> / 3 at small q.
# PhotonScattering.thomsonLimit selects the right limit per process; printing the coherent one beside an incoherent cross
# section produced a column of 0.000000 next to perfectly good numbers, which is how that defect was found.
#
# WHAT IS NOT HERE, and is not hidden: the Klein-Nishina factor and the Compton energy shift, so omega_out = omega_in throughout.
# Both matter once the photon energy approaches m c^2 and neither is implemented. S(q) is likewise the closure approximation
# rather than a Hartree-Fock incoherent scattering function; see PhotonScattering.incoherentScatteringFunction.


"""
`PhotonScattering.computeFormFactorLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model,
                                         grid::Radial.Grid, settings::PhotonScattering.Settings; output=true)`
    ... to compute the coherent or incoherent scattering cross sections of all selected lines in the form-factor approximation,
        and to display them together with the Thomson limit they must approach. An Array{PhotonScattering.Line,1} is returned if
        output=true, and nothing otherwise.

        The approximation is ELASTIC by construction -- the form factor describes a target that returns to its initial state --
        so only lines with the final level equal to the initial one are computed, and any other selected pair is skipped with a
        note rather than answered with a number that would mean nothing.
"""
function computeFormFactorLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid,
                                settings::PhotonScattering.Settings; output=true)
    println("")
    printstyled("PhotonScattering.computeFormFactorLines(): The form-factor computation starts now ... \n", color=:light_green)
    printstyled("------------------------------------------------------------------------------------- \n", color=:light_green)
    #
    newLines = PhotonScattering.Line[];    skipped = 0
    for  iLevel in initialMultiplet.levels
        for  fLevel in finalMultiplet.levels
            if  !Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)    continue    end
            if  iLevel.index != fLevel.index    skipped = skipped + 1;   continue    end
            for  omega in settings.photonEnergies
                inOmega = Defaults.convertUnits("energy: to atomic", omega)
                obs, sigma = PhotonScattering.formFactorObservables(settings.process, inOmega, iLevel, grid, settings)
                push!( newLines, PhotonScattering.Line(iLevel, fLevel, inOmega, inOmega, 0., EmProperty(sigma, sigma),
                                                        PhotonScattering.Channel[], obs) )
            end
        end
    end
    if  skipped > 0
        println(">> $skipped inelastic level pair(s) skipped: the form-factor approximation is ELASTIC by construction, the " *
                "form factor describing a target that returns to its initial state. Use SecondOrderGreen() for those.")
    end

    PhotonScattering.displayFormFactorResults(stdout, newLines, grid, settings)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary   PhotonScattering.displayFormFactorResults(iostream, newLines, grid, settings)     end

    if    output    return( newLines )
    else            return( nothing )
    end
end


"""
`PhotonScattering.displayFormFactorResults(stream::IO, lines::Array{PhotonScattering.Line,1}, grid::Radial.Grid,
                                           settings::PhotonScattering.Settings)`
    ... to print the form-factor cross sections together with the two numbers a reader needs in order to judge them: q a_0 at
        backward scattering, which says how far the form factor has fallen from F(0), and the ratio of the computed cross section
        to the THOMSON LIMIT N^2 (8 pi/3) alpha^4, which must tend to 1 as the photon energy falls. A neat table is printed but
        nothing is returned otherwise.
"""
function displayFormFactorResults(stream::IO, lines::Array{PhotonScattering.Line,1}, grid::Radial.Grid,
                                  settings::PhotonScattering.Settings)
    nx = 132
    alpha = Defaults.getDefaults("alpha");   wc = Defaults.getDefaults("speed of light: c")
    println(stream, " ")
    println(stream, "  Form-factor scattering cross sections   [process = $(settings.process)]:")
    println(stream, " ")
    println(stream, "  ", TableStrings.hLine(nx))
    sa = "  " * TableStrings.center(10, "Level"; na=2) * TableStrings.center(16, "Energy"; na=4) *
                TableStrings.center(12, "N"; na=3) * TableStrings.center(14, "q a_0 (180 deg)"; na=3) *
                TableStrings.center(18, "Cross section"; na=3) *
                TableStrings.center(20, settings.process == PhotonScattering.RayleighScattering() ?
                                        "sigma / N^2 sigma_T" : "sigma / N sigma_T"; na=0)
    println(stream, sa)
    println(stream, "  ", TableStrings.hLine(nx))
    for  line in lines
        N        = FormFactor.standardF(0., line.initialLevel, grid)
        thomson  = PhotonScattering.thomsonLimit(settings.process, N, alpha)
        qBack    = 2 * (line.inPhotonEnergy / wc)
        sa  = "  " * TableStrings.center(10, TableStrings.level(line.initialLevel.index); na=2)
        sa  = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.inPhotonEnergy)) * "    "
        sa  = sa * @sprintf("%10.5f", N) * "     " * @sprintf("%12.5f", qBack) * "     "
        sa  = sa * @sprintf("%.6e", line.crossSection.Coulomb) * "      "
        sa  = sa * @sprintf("%.6f", line.crossSection.Coulomb / thomson)
        println(stream, sa)
    end
    println(stream, "  ", TableStrings.hLine(nx))
    if  length(lines) > 0
        N       = FormFactor.standardF(0., lines[1].initialLevel, grid)
        thomson = PhotonScattering.thomsonLimit(settings.process, N, alpha)
        if  settings.process == PhotonScattering.RayleighScattering()
            println(stream, "  Coherent Thomson limit  N^2 (8 pi/3) alpha^4 = " * @sprintf("%.6e", thomson) *
                            " a.u.  for N = " * @sprintf("%.5f", N) * " electrons.")
            println(stream, "  The last column must tend to 1 as the photon energy FALLS -- every electron scattering in phase --")
            println(stream, "  and fall away once q a_0 approaches 1 and the form factor departs from F(0).")
        else
            println(stream, "  Incoherent Thomson limit  N (8 pi/3) alpha^4 = " * @sprintf("%.6e", thomson) *
                            " a.u.  for N = " * @sprintf("%.5f", N) * " electrons.")
            println(stream, "  The last column must tend to 1 as the photon energy RISES -- the opposite of the coherent case --")
            println(stream, "  since S(q) --> N only once the momentum transfer is large enough to resolve individual electrons.")
            println(stream, "  At low energy it must instead vanish as omega^2, S(q) ~ N q^2 <r^2> / 3 for small q.")
        end
    end
    if  length(settings.polarThetas) > 0  &&  length(lines) > 0
        println(stream, " ")
        println(stream, "  Angle-differential cross sections  d sigma / d Omega  [a.u.]:")
        println(stream, "  ", TableStrings.hLine(nx))
        for  line in lines
            sa = "  " * @sprintf("%.4e", Defaults.convertUnits("energy: from atomic", line.inPhotonEnergy)) * " : "
            for  ob in line.observables    sa = sa * @sprintf("  %.4e", ob.dcs.Coulomb)    end
            println(stream, sa)
        end
        println(stream, "  ", TableStrings.hLine(nx))
        println(stream, "  columns at theta = " * string(round.(settings.polarThetas, digits=4)) * " rad")
    end

    return( nothing )
end


"""
`PhotonScattering.formFactorObservables(process::PhotonScattering.AbstractPhotonProcess, inOmega::Float64, level::Level,
                                        grid::Radial.Grid, settings::PhotonScattering.Settings)`
    ... to compute the angle-differential cross sections at the requested angles and the angle-integrated total, both in the
        form-factor approximation and for the given process. A tuple (observables::Array{PhotonScattering.Observables,1},
        sigma::Float64) is returned.

        The integration over the polar angle is a plain 361-point trapezoid, which is ample: the integrand is
        (1 + cos^2 theta) sin(theta) times a form factor that varies smoothly over the same range, with no structure narrower
        than the grid. The Stokes parameters are returned as NaN throughout, following ParticleScattering's convention for an
        observable that does not apply -- this approximation carries no polarization information at all.
"""
function formFactorObservables(process::PhotonScattering.AbstractPhotonProcess, inOmega::Float64, level::Level,
                               grid::Radial.Grid, settings::PhotonScattering.Settings)
    alpha = Defaults.getDefaults("alpha");    wc = Defaults.getDefaults("speed of light: c")
    reSq  = alpha^4                                     # r_e = alpha^2 in atomic units
    k     = inOmega / wc

    function dcsAt(theta::Float64)
        q  = 2 * k * sin(theta/2)
        wa = FormFactor.standardF(q, level, grid)
        if      process == PhotonScattering.RayleighScattering()   weight = wa * wa
        else                                                       weight = PhotonScattering.incoherentScatteringFunction(q, level, grid)
        end
        return( reSq / 2 * (1 + cos(theta)^2) * weight )
    end

    obs = PhotonScattering.Observables[]
    for  theta in settings.polarThetas
        phi = length(settings.polarPhis) > 0 ? settings.polarPhis[1] : 0.
        d   = dcsAt(theta)
        push!( obs, PhotonScattering.Observables(theta, phi, EmProperty(d, d), NaN, NaN, NaN) )
    end

    # angle-integrated total:  sigma = 2 pi Int_0^pi (d sigma / d Omega) sin(theta) d theta
    nth = 361;   dth = pi / (nth - 1);   sigma = 0.
    for  i = 1:nth
        theta = (i-1) * dth
        w     = (i == 1 || i == nth) ? 0.5 : 1.0
        sigma = sigma + w * dcsAt(theta) * sin(theta) * dth
    end
    sigma = 2pi * sigma

    return( (obs, sigma) )
end


"""
`PhotonScattering.incoherentScatteringFunction(q::Float64, level::Level, grid::Radial.Grid)`
    ... to compute the incoherent scattering function S(q) that governs Compton scattering, in the CLOSURE approximation

            S(q)  =  N  -  |F(q)|^2 / N ,

        with N = F(0) the number of electrons. A value::Float64 is returned.

        THIS IS NOT A HARTREE-FOCK INCOHERENT SCATTERING FUNCTION and must not be quoted as one. JAC's FormFactor module supplies
        the coherent F(q) and nothing else, so S(q) is built from it by closure. The approximation is EXACT in both limits that
        matter -- S(0) = N - N^2/N = 0, a soft photon being unable to excite a bound electron incoherently, and S(q -> infinity)
        = N, every electron scattering independently once the momentum transfer is large -- and interpolated in between, which is
        precisely where the tabulations differ from it by some percent and where it is least to be trusted.

        The pairing with the coherent channel is what makes it useful: |F(q)|^2 + S(q) runs from N^2 at q = 0 to N at large q, so
        the two channels trade against one another with nothing lost, which is a check neither can make alone.
"""
function incoherentScatteringFunction(q::Float64, level::Level, grid::Radial.Grid)
    N = FormFactor.standardF(0., level, grid)
    if  N < 1.0e-12    return( 0. )    end
    F = FormFactor.standardF(q, level, grid)

    return( N - F * F / N )
end


"""
`PhotonScattering.thomsonLimit(process::PhotonScattering.AbstractPhotonProcess, N::Float64, alpha::Float64)`
    ... to return the free-electron cross section that the given process must approach, against which the computed one is worth
        measuring. A sigma::Float64 [a.u.] is returned.

        The two processes approach DIFFERENT limits, and from opposite directions, which is why one number will not serve both.
        Coherent (Rayleigh) scattering tends to N^2 sigma_Thomson as q --> 0, the electrons scattering IN PHASE so that the
        amplitudes add before they are squared. Incoherent (Compton) scattering tends to N sigma_Thomson at LARGE q, where the
        momentum transfer resolves the electrons individually and it is the cross sections that add. Printing a coherent
        normalization beside an incoherent cross section yields a number that looks like a failure and means nothing.
"""
function thomsonLimit(process::PhotonScattering.AbstractPhotonProcess, N::Float64, alpha::Float64)
    sigmaT = 8pi / 3 * alpha^4
    if  process == PhotonScattering.RayleighScattering()    return( N * N * sigmaT )
    else                                                    return( N * sigmaT )
    end
end
