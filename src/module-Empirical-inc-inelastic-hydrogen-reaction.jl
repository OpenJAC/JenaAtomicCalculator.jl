
#################################################################################################################################
### Inelastic ion -- hydrogen REACTIONS (configuration-level interface) #########################################################
##
##  A higher-level interface to the Belyaev-Yakovleva simplified Landau-Zener model (module-Empirical-inc-inelastic-h-collisions.jl)
##  and the molecular-symmetry channel correlation (module-Empirical-inc-hydrogen-channel-correlation.jl): instead of hand-building
##  Empirical.AtomicLevel/InelasticHChannel objects with manually-assigned L, S, J and energies, an Empirical.InelasticHReaction
##  specifies the physical reaction directly in terms of JAC's own entities -- Nuclear.Model and Configuration -- exactly the way
##  any other JAC structure calculation is specified. Everything else (reduced mass, molecular-term correlation, statistical
##  weights) is then DERIVED, not supplied by the caller. Level ENERGIES are the one piece this module deliberately does NOT
##  compute itself -- they are always supplied by the caller as a Dict{Configuration,Float64} of total energies, whether that
##  dictionary was filled from a real Atomic.Computation/perform SCF run, from literature (NIST-type) values, or from anything
##  else; this keeps the Empirical module free of a dependency on the (much heavier) Atomic/SelfConsistent machinery beyond what
##  it already needs, and cleanly separates "how were these energies obtained" (a caller-level decision, worth reporting) from
##  "how are rates computed from given energies" (this module's actual job).
##
##  Scope (deliberately restricted, matching the physics already validated for Ba2+ + H- -> Ba+ + H): the entrance ion must be
##  closed-shell (so its own molecular symmetry is the trivial, unique 1S0 -- Lion=0, Sion=0.0, no term-generation needed), and
##  every final ion configuration must differ from the entrance configuration by exactly ONE electron in exactly ONE shell (a
##  genuine single-active-electron transfer, S=1/2 always). Configurations that don't fit this pattern raise an informative
##  error rather than being silently mishandled or guessed at; open-shell entrance ions would need real term generation from
##  the configuration, not yet implemented here.


"""
`struct  Empirical.InelasticHReaction`
    ... specifies an inelastic ion + hydrogen reaction A^(Z+1)+ + H^- <-> A^Z+(f) + H directly in terms of JAC's own Nuclear.Model
        and Configuration entities, for one or several final configurations at once.

    + nm        ::Nuclear.Model            ... the projectile ion's nuclear model; only nm.Z is used here (for the covalent
                                                species' charge and for looking up its standard atomic mass).
    + iConfIon  ::Configuration             ... the (closed-shell) ground configuration of the ionic entrance species A^(Z+1)+.
    + iConfH    ::Configuration             ... either H^-(1s^2) or H(1s^1) -- this module only ever considers hydrogen.
    + fConfIon  ::Array{Configuration,1}    ... the final ion configurations A^Z+(f) of interest; each must differ from iConfIon
                                                by exactly one electron in exactly one shell.
    + fConfH    ::Configuration             ... the complementary final hydrogen configuration; derived automatically from
                                                iConfH by electron-count conservation if not supplied explicitly.
"""
struct  InelasticHReaction
    nm          ::Nuclear.Model
    iConfIon    ::Configuration
    iConfH      ::Configuration
    fConfIon    ::Array{Configuration,1}
    fConfH      ::Configuration
end


"""
`Empirical.InelasticHReaction(nm::Nuclear.Model, iConfIon::Configuration, iConfH::Configuration,
                              fConfIon::Array{Configuration,1}; fConfH::Union{Nothing,Configuration}=nothing)`
    ... outer constructor that validates iConfH/fConfH and derives fConfH from iConfH by electron-count conservation if not
        given explicitly. Raises an informative error if iConfH is not H(1s^1) or H^-(1s^2), or if a given fConfH is
        inconsistent with iConfH's electron count. An Empirical.InelasticHReaction is returned.
"""
function InelasticHReaction(nm::Nuclear.Model, iConfIon::Configuration, iConfH::Configuration,
                            fConfIon::Array{Configuration,1}; fConfH::Union{Nothing,Configuration}=nothing)
    Empirical.validateHydrogenConfiguration(iConfH)
    if  isnothing(fConfH)
        fConfHx = iConfH.NoElectrons == 2 ? Configuration("1s^1") : Configuration("1s^2")
    else
        Empirical.validateHydrogenConfiguration(fConfH)
        fConfHx = fConfH
    end
    if  abs(fConfHx.NoElectrons - iConfH.NoElectrons) != 1
        error("Empirical.InelasticHReaction: iConfH ($iConfH) and fConfH ($fConfHx) must differ by exactly one electron.")
    end
    return( InelasticHReaction(nm, iConfIon, iConfH, fConfIon, fConfHx) )
end


"""
`Empirical.validateHydrogenConfiguration(conf::Configuration)`
    ... to check that conf is either H(1s^1) or H^-(1s^2); raises an informative error otherwise, since
        Empirical.InelasticHReaction is restricted to hydrogen collision partners by design (see the module note above).
"""
function validateHydrogenConfiguration(conf::Configuration)
    if  conf != Configuration("1s^1")  &&  conf != Configuration("1s^2")
        error("Empirical.InelasticHReaction is restricted to hydrogen collision partners: expected H(1s^1) or H^-(1s^2), " *
              "got $conf.")
    end
end


"""
`Empirical.isClosedShell(conf::Configuration)`
    ... to check whether every shell of conf is fully occupied (occupation = 2*(2l+1) for each shell). A Bool is returned.
"""
function isClosedShell(conf::Configuration)
    for  (sh, occ)  in  conf.shells
        if  occ != 2*(2*sh.l + 1)   return( false )   end
    end
    return( true )
end


"""
`Empirical.activeShell(iConf::Configuration, fConf::Configuration)`
    ... to identify the single shell whose occupation differs between iConf and fConf, and check that it differs by exactly
        one electron (a genuine one-electron transfer) -- the scope Empirical.InelasticHReaction is restricted to (see the
        module note above). Raises an informative error if more than one shell differs, or if the occupation change is not
        exactly one electron. A sh::Shell is returned.
"""
function activeShell(iConf::Configuration, fConf::Configuration)
    allShells  = union(keys(iConf.shells), keys(fConf.shells))
    diffShells = Shell[]
    for  sh  in  allShells
        if  get(iConf.shells, sh, 0) != get(fConf.shells, sh, 0)   push!(diffShells, sh)   end
    end
    if  length(diffShells) != 1
        error("Empirical.InelasticHReaction (this version) requires iConfIon and each fConfIon to differ in exactly one " *
              "shell (a single-active-electron transfer); found $(length(diffShells)) differing shells between " *
              "$iConf and $fConf.")
    end
    sh   = diffShells[1]
    dOcc = get(fConf.shells, sh, 0) - get(iConf.shells, sh, 0)
    if  dOcc != 1
        error("Empirical.InelasticHReaction (this version) requires the active shell to GAIN exactly one electron; " *
              "shell $sh changes by $dOcc electrons between $iConf and $fConf.")
    end
    return( sh )
end


"""
`Empirical.reducedMassH(nm::Nuclear.Model; ionMass::Union{Nothing,Float64}=nothing)`
    ... to compute the reduced mass of the projectile ion (nuclear charge nm.Z) with a hydrogen atom, in electron-mass atomic
        units. The projectile's mass is taken from PeriodicTable.getData("mass", nm.Z) (the standard atomic weight, in amu) by
        default; note that nm.mass itself is NOT usable here -- it is a generic A ~ 2Z+0.005Z^2 formula used only to derive the
        nuclear radius, not a real isotope mass. Pass ionMass [amu] explicitly to use a specific isotope instead. A
        value::Float64 [a.u.] is returned.
        Note: this rate model is thermally averaged over energy and summed over many partial waves, which washes out most of
              the mass-dependence of any single collision; even a +/-20% mass error changes typical rate coefficients by well
              under 1% (checked explicitly for Ba+ + H). Getting the mass exactly right is far less important than getting the
              level energies right for this model.
"""
function reducedMassH(nm::Nuclear.Model; ionMass::Union{Nothing,Float64}=nothing)
    Mion_amu = isnothing(ionMass) ? PeriodicTable.getData("mass", round(Int64, nm.Z)) : ionMass
    Mion = Mion_amu / Defaults.ELECTRON_MASS_U
    MH   = PeriodicTable.getData("mass", 1) / Defaults.ELECTRON_MASS_U
    return( Mion * MH / (Mion + MH) )
end


"""
`Empirical.energyOf(energies::Array{Pair{Configuration,Float64},1}, conf::Configuration)`
    ... to look up conf's total energy [a.u.] in energies by VALUE equality (Configuration == is defined and reliable, but
        Configuration currently has no matching `hash` method, so a genuine Dict{Configuration,Float64} silently drops valid
        keys -- a real, separate JAC bug worth fixing centrally at some point, not something to route around by relying on
        Dict here). Raises an informative error if conf is not found. A value::Float64 is returned.
"""
function energyOf(energies::Array{Pair{Configuration,Float64},1}, conf::Configuration)
    idx = findfirst(p -> p.first == conf, energies)
    if  isnothing(idx)   error("energies has no entry for configuration $conf.")   end
    return( energies[idx].second )
end


"""
`Empirical.inelasticHCollisionRateMatrix(T::Float64, reaction::Empirical.InelasticHReaction,
                                         energies::Array{Pair{Configuration,Float64},1};
                                         energyLabel::String="externally supplied", printout::Bool=false, zerosGL::Int64=64)`
    ... convenience method that builds the Array{Empirical.InelasticHChannel,1} directly from reaction -- reduced mass via
        Empirical.reducedMassH, molecular symmetry/statistical weight via Empirical.activeShell + Empirical.statisticalWeight
        (assuming a closed-shell, Lion=0/Sion=0.0, entrance ion) -- and dispatches to
        Empirical.inelasticHCollisionRateMatrix(T, channels, pstatIonic, Z, mu; ...). energies must contain a total-energy
        entry [a.u.] (via Empirical.energyOf, e.g. `[Configuration("[Xe]") => 0.0, Configuration("[Xe] 6s^1") => -0.368]`)
        for reaction.iConfIon and for every reaction.fConfIon; this module never computes these itself (see the module note
        above) -- energyLabel is a short caller-supplied description (e.g. "single-configuration Dirac-Fock SCF" or "NIST")
        echoed in the printout so the origin of the numbers is never silently lost. A named tuple, as returned by the
        low-level method, is returned.
"""
function inelasticHCollisionRateMatrix(T::Float64, reaction::Empirical.InelasticHReaction,
                                       energies::Array{Pair{Configuration,Float64},1};
                                       energyLabel::String="externally supplied", printout::Bool=false, zerosGL::Int64=64)
    if  !Empirical.isClosedShell(reaction.iConfIon)
        error("Empirical.InelasticHReaction (this version) requires a closed-shell entrance ion configuration; " *
              "$(reaction.iConfIon) is not closed-shell. An open-shell entrance needs real term generation, not yet " *
              "implemented here.")
    end
    Zcov = reaction.nm.Z - reaction.fConfIon[1].NoElectrons
    for  fConf  in  reaction.fConfIon
        if  reaction.nm.Z - fConf.NoElectrons != Zcov
            error("All entries of reaction.fConfIon must have the same electron count; $(fConf) does not match the others.")
        end
    end

    mu    = Empirical.reducedMassH(reaction.nm)
    E_ion = Empirical.energyOf(energies, reaction.iConfIon)

    channels = Empirical.InelasticHChannel[]
    for  fConf  in  reaction.fConfIon
        Ej    = Empirical.energyOf(energies, fConf) - E_ion
        sh    = Empirical.activeShell(reaction.iConfIon, fConf)
        pstat = Empirical.statisticalWeight(sh.l, 0.5, 0, 0.0)
        push!(channels, Empirical.InelasticHChannel(string(fConf), Ej, pstat))
    end

    if  printout
        println("\n* Empirical.InelasticHReaction:  A^$(round(Int64,Zcov+1))+ ($(reaction.iConfIon)) + $(reaction.iConfH)  " *
                "->  A^$(round(Int64,Zcov))+ (f) + $(reaction.fConfH):" *
                "\n    + Energy source: $energyLabel -- not computed by this module; supplied by the caller via `energies`." *
                "\n    + Reduced mass: $(round(mu,digits=1)) [a.u.], from PeriodicTable's standard atomic weight for Z=" *
                "$(round(Int64,reaction.nm.Z)) (nm.mass itself is a nuclear-radius placeholder, not used here); this rate " *
                "model is insensitive to mass to well under 1% even for a 20% mass error, so this choice is not critical." *
                "\n    + Statistical weights assume a closed-shell (1S0) entrance ion and a single active electron per " *
                "final configuration -- exact for this scope, not a general open-shell treatment." *
                "\n    + Final-state channels and their derived (energy, p_stat):")
        for  (fConf, ch)  in  zip(reaction.fConfIon, channels)
            println("        $fConf:  E = $(round(Defaults.convertUnits("energy: from atomic to eV",ch.E),digits=4)) eV   " *
                    "p_stat = $(round(ch.pstat,digits=6))")
        end
    end

    return( Empirical.inelasticHCollisionRateMatrix(T, channels, 1.0, Zcov, mu; printout=printout, zerosGL=zerosGL) )
end
