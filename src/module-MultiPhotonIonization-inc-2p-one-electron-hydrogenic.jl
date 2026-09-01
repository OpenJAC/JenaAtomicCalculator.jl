
# Two-photon ionization of a ONE-ELECTRON system:   gamma + gamma + |i> --> |eps kappa_f>
#
# WRITTEN FROM SCRATCH 25-Aug-2026.  What stood here before produced partial amplitudes and stopped: it printed one
# number per (multipole, intermediate symmetry, final symmetry) without ever summing the intermediate symmetries
# coherently, without forming a cross section, and without a guard on the energy denominator.  None of it is kept.
#
# WHY A ONE-ELECTRON FILE AT ALL, when the scheme above it handles many-electron ions.  Two reasons, and the second
# is the important one:
#   + IT IS THE ONLY INDEPENDENT CHECK.  Two-photon ionization of hydrogen is a standard benchmark, and here the
#     initial state, the intermediate spectrum and the continuum are all known exactly, so a wrong prefactor or a
#     wrong recoupling has nowhere to hide.
#   + IT IS WHERE SUMMING OVER THE CONTINUUM IS REALISTIC.  The intermediate sum of a second-order process runs
#     over the WHOLE spectrum, bound and free.  A bound sum alone is not a small error: for hydrogen at photon
#     energies above half the ionization potential the continuum carries a substantial part of the amplitude.  In
#     a one-electron system a Green's function can actually be constructed, which is not yet true of the
#     many-electron case.
#
# THE RESONANCE GUARD REFUSES rather than skips.  If a photon energy brings a real intermediate state onto the
# energy shell, the non-resonant second-order expression does not apply at all -- the amplitude diverges and the
# process becomes a two-step sequential one with a genuine lifetime.  Skipping the offending term, as a bound-bound
# module may reasonably do among hundreds of intermediates, would here remove most of the amplitude and return a
# small number instead of no number.  This routine therefore names the resonant state and stops.


"""
`MultiPhotonIonization.reducedAmplitude2pOneElectron(K::AngularJ64, fOrbital::Orbital, omega2::Float64, mp2::EmMultipole,
                                                     symx::LevelSymmetry, omega1::Float64, mp1::EmMultipole,
                                                     iOrbital::Orbital, gauge::EmGauge, orbitals::Dict{Subshell,Orbital},
                                                     grid::Radial.Grid; resonanceTolerance::Float64=1.0e-3)`  
    ... computes the second-order reduced amplitude U^K for one intermediate symmetry symx, i.e. the inner sum over
        all intermediate orbitals of that symmetry.  A U::ComplexF64 is returned.

        The energy denominator is E_i + omega1 - E_v.  If any intermediate comes within `resonanceTolerance` of the
        energy shell the routine RAISES, naming that state: a resonant denominator means the non-resonant
        expression does not describe the process, and returning a number would be worse than returning none.
"""
function  reducedAmplitude2pOneElectron(K::AngularJ64, fOrbital::Orbital, omega2::Float64, mp2::EmMultipole,
                                        symx::LevelSymmetry, omega1::Float64, mp1::EmMultipole, iOrbital::Orbital,
                                        gauge::EmGauge, orbitals::Dict{Subshell,Orbital}, grid::Radial.Grid;
                                        resonanceTolerance::Float64=1.0e-3)
    U = ComplexF64(0.);    found = false
    for  (sh, vOrbital)  in  orbitals
        if  LevelSymmetry(sh) != symx    continue    end
        denom = iOrbital.energy + omega1 - vOrbital.energy
        ## THE TOLERANCE IS RELATIVE TO THE PHOTON ENERGY, and that is not a detail.  An ABSOLUTE 1e-6 a.u. was
        ## tried first and proved useless: asking for hydrogen at omega = 0.375 a.u., which is the 1s -> 2p
        ## resonance, left a denominator of 4.6e-6 -- inflating that one term by a factor 2e5 -- and the guard
        ## passed it, because 4.6e-6 exceeds 1e-6.  A guard that lets the textbook resonance through is not a
        ## guard.  Where "near resonance" begins is physics rather than arithmetic, so the number is a parameter;
        ## the default of 1e-3*omega is about 0.01 eV at optical energies.
        if  abs(denom) < resonanceTolerance * omega1
            error("MultiPhotonIonization: RESONANT INTERMEDIATE STATE.  The intermediate $sh lies within " *
                  "$(resonanceTolerance) x omega of the energy shell at omega1 = $omega1 a.u. (denominator " *
                  "$(denom)).  The non-resonant second-order expression does not apply there: the process becomes " *
                  "a sequential two-step one through a real state with its own lifetime.  Move omega away from " *
                  "this resonance, or treat the sequential process explicitly; this routine will not return a " *
                  "number it cannot stand behind.")
        end
        found = true
        ## MabEmission(mp, gauge, omega, X, Y, grid) returns <Y || O || X>, so the chain reads
        ##     <f || O2 || v>  <v || O1 || i>
        U = U + InteractionStrength.MabEmission(mp2, gauge, omega2, vOrbital, fOrbital, grid) *
                InteractionStrength.MabEmission(mp1, gauge, omega1, iOrbital, vOrbital, grid) / denom
    end
    if  !found    return( ComplexF64(0.) )    end
    ## the two photons are recoupled to the total transferred rank K
    U = U * AngularMomentum.Wigner_6j(AngularMomentum.kappa_j(iOrbital.subshell.kappa),
                                      AngularMomentum.kappa_j(fOrbital.subshell.kappa), K,
                                      AngularJ64(mp2.L), AngularJ64(mp1.L), symx.J)

    return( U )
end


"""
`MultiPhotonIonization.crossSection2pOneElectron(polarization::Symbol, iOrbital::Orbital, omega1::Float64,
                                                 omega2::Float64, multipoles::Array{EmMultipole,1}, gauge::EmGauge,
                                                 orbitals::Dict{Subshell,Orbital}, meanPot::Radial.Potential;
                                                 resonanceTolerance::Float64=1.0e-3)`  
    ... computes the generalized two-photon ionization cross section of a one-electron system.  A tuple
        (sigma::Float64 in ATOMIC UNITS, finalSymmetries::Array{LevelSymmetry,1}) is returned.

        `polarization` is :linear or :unpolarized, and the two differ in exactly the way they do for two-photon
        absorption: for LINEAR light the two helicity components are summed COHERENTLY and then squared, while for
        UNPOLARIZED light each helicity pair is squared on its own and only odd ranks K are dropped.

        THE FINAL CHANNELS ARE SUMMED INCOHERENTLY, which is the one structural difference from the bound-bound
        case.  There, one final level receives everything; here each final partial wave kappa_f is a distinct
        state of the ejected electron, so their cross sections add rather than their amplitudes.
"""
function  crossSection2pOneElectron(polarization::Symbol, iOrbital::Orbital, omega1::Float64, omega2::Float64,
                                    multipoles::Array{EmMultipole,1}, gauge::EmGauge,
                                    orbitals::Dict{Subshell,Orbital}, meanPot::Radial.Potential;
                                    resonanceTolerance::Float64=1.0e-3)
    grid    = meanPot.grid
    epsilon = omega1 + omega2 + iOrbital.energy
    if  epsilon <= 0.
        error("MultiPhotonIonization: the two photons do not reach the continuum -- omega1 + omega2 = " *
              "$(omega1+omega2) a.u. against a binding energy of $(-iOrbital.energy) a.u.  Two-photon ionization " *
              "needs omega1 + omega2 > |E_i|.")
    end
    symi = LevelSymmetry(iOrbital.subshell);    ji = AngularMomentum.kappa_j(iOrbital.subshell.kappa)
    ## collect the reachable final symmetries once, together with the intermediates that reach them
    routes = Tuple{LevelSymmetry,LevelSymmetry,EmMultipole,EmMultipole}[]
    for  mp1 in multipoles,  mp2 in multipoles
        for  symx in AngularMomentum.allowedMultipoleSymmetries(symi, mp1)
            for  symf in AngularMomentum.allowedMultipoleSymmetries(symx, mp2)
                push!(routes, (symf, symx, mp1, mp2))
            end
        end
    end
    finalSyms = unique( r[1] for r in routes )
    #
    sigma = 0.
    for  symf  in  finalSyms
        ## the continuum orbital of this final symmetry, at the energy fixed by energy conservation
        nrContinuum = Continuum.gridConsistency(epsilon, grid)
        cSettings   = Continuum.Settings(false, nrContinuum)
        fOrbital, _ = Continuum.generateOrbitalLocalPotential(epsilon, Subshell(1001, symf), meanPot, cSettings)
        jf          = AngularMomentum.kappa_j(fOrbital.subshell.kappa)
        for  K  in  AngularMomentum.j_values(ji, jf)
            if  polarization == :unpolarized  &&  isodd( Int(Basics.twice(K)/2) )    continue    end
            for  q  in  AngularMomentum.m_values(K)
                ampAll = ComplexF64(0.)
                for  lambda1 in [-1, 1],  lambda2 in [-1, 1]
                    ampPair = ComplexF64(0.)
                    for  (sf, sx, mp1, mp2)  in  routes
                        if  sf != symf    continue    end
                        p1 = mp1.electric ? 1 : 0;    p2 = mp2.electric ? 1 : 0
                        U  = MultiPhotonIonization.reducedAmplitude2pOneElectron(K, fOrbital, omega2, mp2, sx,
                                        omega1, mp1, iOrbital, gauge, orbitals, grid;
                                        resonanceTolerance=resonanceTolerance)
                        if  U == ComplexF64(0.)   continue   end
                        w = (1.0im)^(mp1.L - p1 + mp2.L - p2) * (-lambda1)^p1 * (-lambda2)^p2 *
                            sqrt( (2*mp1.L + 1)*(2*mp2.L + 1) ) * (Basics.twice(K) + 1) *
                            AngularMomentum.Wigner_3j(mp1.L, mp2.L, K, lambda1, lambda2,
                                                      AngularMomentum.oneM(q))
                        ampPair = ampPair + w * U
                    end
                    if  polarization == :unpolarized    sigma  = sigma + abs(ampPair)^2
                    else                                ampAll = ampAll + ampPair
                    end
                end
                if  polarization != :unpolarized    sigma = sigma + abs(ampAll)^2    end
            end
        end
    end
    ## the same prefactor as MultiPhotonTransition's two-photon absorption; sigma is then in atomic units, i.e.
    ## in a0^4 * t_au, and one atomic unit is 1.8966 GM
    sigma = sigma * 2*pi^5 / Defaults.getDefaults("alpha")^2 / (AngularMomentum.oneJ(ji)*2 + 1) / omega1^2

    return( (sigma, finalSyms) )
end


"""
`struct  MultiPhotonIonization.Line2pOneElectron`  
    ... one photon energy of a two-photon one-electron ionization, with the generalized cross sections that belong
        to it.

    + iSubshell       ::Subshell             ... the bound subshell that is ionized.
    + omega           ::Float64              ... photon energy [a.u.]; both photons carry it.
    + electronEnergy  ::Float64              ... energy of the ejected electron, 2*omega + E_i, fixed by conservation.
    + csLinear        ::Basics.EmProperty    ... generalized cross section for linearly polarized light [a.u.].
    + csUnpolarized   ::Basics.EmProperty    ... the same for unpolarized light [a.u.].
    + finalSymmetries ::Array{LevelSymmetry,1}  ... the final partial waves that contribute.
"""
struct  Line2pOneElectron
    iSubshell         ::Subshell
    omega             ::Float64
    electronEnergy    ::Float64
    csLinear          ::Basics.EmProperty
    csUnpolarized     ::Basics.EmProperty
    finalSymmetries   ::Array{LevelSymmetry,1}
end


"""
`MultiPhotonIonization.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model,
                                    grid::Radial.Grid, settings::MultiPhotonIonization.Settings; output::Bool=true)`
    ... the entrance from `Atomic.Computation` and `Basics.perform`, so that a multi-photon study is set up exactly
        like every other process in JAC rather than by calling into this file directly; a
        `lines::Array{Line2pOneElectron,1}` is returned, or `nothing` when `output=false`.

        WHAT THIS DOES AND DOES NOT DO. It is an ADAPTER, not a generalisation. The physics behind it is written
        for ONE electron in a central potential -- that is what `-inc-2p-one-electron-hydrogenic` means -- and this
        function does not change that. It takes the general (initialMultiplet, finalMultiplet) route, checks that
        the system really has a single electron, extracts the occupied subshell, and hands the work to
        `computeLines2pOneElectron`. A many-electron system is REFUSED with a message that says why, rather than
        being silently reduced to something this module can do.

        THE INTERMEDIATE SPECTRUM COMES FROM `settings.intermediateStates`, as in every other second-order module
        of JAC. A second-order amplitude is a sum over a spectrum, and the spectrum is the caller's to choose: a
        mean-field multiplet, or a Green-function expansion. This function therefore does NOT build one of its
        own; if none is supplied it says so instead of quietly summing over whatever orbitals happen to be in the
        initial basis, which would make the result depend on how the initial state was generated.
"""
function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model,
                       grid::Radial.Grid, settings::MultiPhotonIonization.Settings; output::Bool=true)
    # (1) one electron, or nothing
    NoElectrons = length(initialMultiplet.levels) == 0  ?  0  :  initialMultiplet.levels[1].basis.NoElectrons
    if  NoElectrons != 1
        error("\n\nMultiPhotonIonization: this module computes two-photon ionization of a ONE-ELECTRON system " *
              "and the initial state has $NoElectrons.\n"                                                        *
              ">>> The amplitudes are built for a single electron in a central potential -- the file name says\n" *
              ">>> so -- and there is no many-electron version to fall back on. Extending it is real physics,\n"  *
              ">>> not a keyword: the second-order sum would need the full CSF machinery on both steps.\n"        *
              ">>> Use a hydrogen-like initial configuration, or drive a many-electron two-photon process\n"      *
              ">>> through MultiPhotonTransition, which is bound-bound.\n")
    end

    # (2) the intermediate spectrum is supplied, never invented
    if  typeof(settings.intermediateStates) == Multiplet  &&  length(settings.intermediateStates.levels) == 0
        error("\n\nMultiPhotonIonization: no intermediate states were given.\n"                                 *
              ">>> A second-order amplitude is a SUM OVER A SPECTRUM, and which spectrum it is changes the\n"     *
              ">>> answer, so it is the caller's choice rather than this module's. Set\n"                         *
              ">>>     MultiPhotonIonization.Settings(...; intermediateStates = <a Multiplet or GreenChannels>)\n" *
              ">>> exactly as MultiPhotonTransition and PhotonScattering are set up. Summing instead over the\n"  *
              ">>> orbitals that happen to sit in the initial basis would make the result depend on how the\n"    *
              ">>> initial state was generated, which is why it is refused here.\n")
    end

    # (3) the occupied subshell of the one electron, and the orbitals available for the sum
    iSubshell = Subshell[]
    for  (sh, occ)  in  initialMultiplet.levels[1].basis.csfs[1].occupation |> x -> zip(initialMultiplet.levels[1].basis.subshells, x)
        if  occ > 0    push!(iSubshell, sh)    end
    end
    length(iSubshell) == 1  ||  error("MultiPhotonIonization: the one-electron initial state occupies " *
                                      "$(length(iSubshell)) subshells; expected exactly one.")
    orbitals = Dict{Subshell,Orbital}()
    for  (sh, orb)  in  initialMultiplet.levels[1].basis.orbitals      orbitals[sh] = orb      end
    if  typeof(settings.intermediateStates) == Multiplet
        for  lev  in  settings.intermediateStates.levels
            for  (sh, orb)  in  lev.basis.orbitals    orbitals[sh] = orb    end
        end
    end

    scheme = settings.scheme
    return( MultiPhotonIonization.computeLines2pOneElectron(iSubshell[1], scheme.omegas, scheme.multipoles,
                    orbitals, Nuclear.nuclearPotential(nm, grid);
                    resonanceTolerance = scheme.resonanceTolerance, output = output) )
end


"""
`MultiPhotonIonization.computeLines2pOneElectron(iSubshell::Subshell, omegas::Array{Float64,1},
                                                 multipoles::Array{EmMultipole,1}, orbitals::Dict{Subshell,Orbital},
                                                 meanPot::Radial.Potential; resonanceTolerance::Float64=1.0e-3,
                                                 output::Bool=true)`  
    ... computes the generalized two-photon ionization cross sections of a one-electron system, for each of the
        given photon energies and in both gauges.  An Array{Line2pOneElectron,1} is returned if output is true.
"""
function  computeLines2pOneElectron(iSubshell::Subshell, omegas::Array{Float64,1},
                                    multipoles::Array{EmMultipole,1}, orbitals::Dict{Subshell,Orbital},
                                    meanPot::Radial.Potential; resonanceTolerance::Float64=1.0e-3,
                                    output::Bool=true)
    println("")
    printstyled("MultiPhotonIonization: two-photon ionization of a one-electron system starts now ... \n",
                color=:light_green)
    printstyled("----------------------------------------------------------------------------------- \n",
                color=:light_green)
    iOrbital = orbitals[iSubshell];    lines = MultiPhotonIonization.Line2pOneElectron[]
    for  omega  in  omegas
        csL = Basics.EmProperty(0.);   csU = Basics.EmProperty(0.);   fSyms = LevelSymmetry[]
        cL = Float64[];   cU = Float64[]
        for  gauge  in  [Basics.Coulomb, Basics.Babushkin]
            (sL, fSyms) = MultiPhotonIonization.crossSection2pOneElectron(:linear, iOrbital, omega, omega,
                              multipoles, gauge, orbitals, meanPot; resonanceTolerance=resonanceTolerance)
            (sU, _)     = MultiPhotonIonization.crossSection2pOneElectron(:unpolarized, iOrbital, omega, omega,
                              multipoles, gauge, orbitals, meanPot; resonanceTolerance=resonanceTolerance)
            push!(cL, sL);   push!(cU, sU)
        end
        csL = Basics.EmProperty(cL[1], cL[2]);    csU = Basics.EmProperty(cU[1], cU[2])
        push!( lines, MultiPhotonIonization.Line2pOneElectron(iSubshell, omega, 2*omega + iOrbital.energy,
                                                              csL, csU, fSyms) )
    end
    MultiPhotonIonization.displayLines2pOneElectron(stdout, lines)
    printSummary, iostream = Defaults.getDefaults("summary flag/stream")
    if  printSummary    MultiPhotonIonization.displayLines2pOneElectron(iostream, lines)    end

    if  output    return( lines )   else   return( nothing )   end
end


"""
`MultiPhotonIonization.displayLines2pOneElectron(stream::IO, lines::Array{Line2pOneElectron,1})`  
    ... displays the generalized two-photon ionization cross sections, in atomic units and in GM.  Nothing is
        returned.
"""
function  displayLines2pOneElectron(stream::IO, lines::Array{Line2pOneElectron,1})
    nx = 122
    ## one atomic unit of a generalized two-photon cross section is a0^4 * t_au
    auToCm4s = Defaults.convertUnits("length: from atomic to cm", 1.0)^4 *
               Defaults.convertUnits("time: from atomic to sec", 1.0)
    println(stream, " ")
    println(stream, "  Generalized TWO-PHOTON IONIZATION cross sections of a one-electron system:")
    println(stream, " ")
    println(stream, "  ", "-"^nx)
    println(stream, "    omega [a.u.]   eps(e-) [a.u.]     sigma^(2) linear [GM]        sigma^(2) unpolarized [GM]" *
                    "      final waves")
    println(stream, "                                      Coulomb      Babushkin       Coulomb      Babushkin")
    println(stream, "  ", "-"^nx)
    for  line  in  lines
        sa = "   " * @sprintf("%10.5f", line.omega) * "     " * @sprintf("%12.5f", line.electronEnergy) * "   "
        for  cs  in  [line.csLinear, line.csUnpolarized]
            sa = sa * @sprintf("%14.5e", cs.Coulomb * auToCm4s / 1.0e-50) *
                      @sprintf("%14.5e", cs.Babushkin * auToCm4s / 1.0e-50)
        end
        println(stream, sa * "    " * string(length(line.finalSymmetries)))
    end
    println(stream, "  ", "-"^nx)
    println(stream, "    UNITS AND CONVENTION.  A two-photon process needs two photons at once, so its rate goes as")
    println(stream, "    the SQUARE of the photon flux:  W [1/s] = sigma^(2) * F^2,  with F in photons cm^-2 s^-1.")
    println(stream, "    A two-photon cross section is therefore NOT an area: it carries cm^4 s, and is quoted here")
    println(stream, "    in GM (1 GM = 1e-50 cm^4 s, after Maria Goeppert-Mayer, who predicted the process in 1931).")
    println(stream, "    One atomic unit = " * @sprintf("%.4f", auToCm4s/1.0e-50) * " GM.")
    println(stream, "    THE SINGLE-BEAM FACTOR IS A CONVENTION, pinned as in MultiPhotonTransition: for one beam")
    println(stream, "    the two photons are indistinguishable, and the factor is fixed by requiring the")
    println(stream, "    monochromatic result to agree with the bichromatic one, W = sigma^(2) F_1 F_2, as")
    println(stream, "    omega_1 -> omega_2.  Quoting these numbers against another source means checking that the")
    println(stream, "    source pins it the same way.")
    println(stream, "  ", "-"^nx)

    return( nothing )
end
