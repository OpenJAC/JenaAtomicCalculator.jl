
# How the target is represented to the projectile. This is the axis along which ELSEPA is organised, and it is kept
# separate from the projectile and the beam: a further model is one new subtype of AbstractInteractionModel plus one
# new method of ParticleScattering.scatteringPotential below.
#
# The potentials are taken from JAC's own machinery rather than from tabulated Dirac-Fock densities, which is what
# makes the input "typically for JAC": Basics.CHField() is the static field of the target density without exchange,
# and Basics.DFSField(1.0) is that field plus the local Slater exchange term.


"""
`ParticleScattering.scatteringPotential(projectile::ParticleScattering.Electron, interaction::ParticleScattering.StaticField,
                                        level::Level, nm::Nuclear.Model, grid::Radial.Grid;
                                        nuclearPot::Union{Nothing,Radial.Potential}=nothing)`
    ... to build the local potential in which an ELECTRON partial wave is generated, for the pure static field, i.e. the
        electrostatic interaction with the nuclear and electronic charge density and no exchange term. The nuclear part
        depends only on (nm, grid) and may be passed in by a caller that has built it once. A pot::Radial.Potential is
        returned.
"""
function scatteringPotential(projectile::ParticleScattering.Electron, interaction::ParticleScattering.StaticField,
                             level::Level, nm::Nuclear.Model, grid::Radial.Grid;
                             nuclearPot::Union{Nothing,Radial.Potential}=nothing)
    nucPot = isnothing(nuclearPot) ? Nuclear.nuclearPotential(nm, grid) : nuclearPot
    wp     = Basics.computePotential(Basics.CHField(), grid, level)

    return( Basics.add(nucPot, wp) )
end


"""
`ParticleScattering.scatteringPotential(projectile::ParticleScattering.Electron, interaction::ParticleScattering.StaticFieldExchange,
                                        level::Level, nm::Nuclear.Model, grid::Radial.Grid;
                                        nuclearPot::Union{Nothing,Radial.Potential}=nothing)`
    ... to build the local potential in which an ELECTRON partial wave is generated, for the static field plus a local
        exchange term. The exchange is the Slater term carried by Basics.DFSField(1.0), which is the same field the
        bound-state machinery of JAC uses, so that projectile and target are described consistently.
        A pot::Radial.Potential is returned.
"""
function scatteringPotential(projectile::ParticleScattering.Electron, interaction::ParticleScattering.StaticFieldExchange,
                             level::Level, nm::Nuclear.Model, grid::Radial.Grid;
                             nuclearPot::Union{Nothing,Radial.Potential}=nothing)
    nucPot = isnothing(nuclearPot) ? Nuclear.nuclearPotential(nm, grid) : nuclearPot
    wp     = Basics.computePotential(Basics.DFSField(1.0), grid, level)

    return( Basics.add(nucPot, wp) )
end


"""
`ParticleScattering.scatteringPotential(projectile::ParticleScattering.Positron, interaction::ParticleScattering.StaticField,
                                        level::Level, nm::Nuclear.Model, grid::Radial.Grid;
                                        nuclearPot::Union{Nothing,Radial.Potential}=nothing)`
    ... to build the local potential in which a POSITRON partial wave is generated. The electrostatic interaction changes
        sign throughout -- the nucleus repels and the electrons attract -- so the electron potential is simply negated;
        no exchange term applies, a positron being distinguishable from the target electrons.
        A pot::Radial.Potential is returned.
"""
function scatteringPotential(projectile::ParticleScattering.Positron, interaction::ParticleScattering.StaticField,
                             level::Level, nm::Nuclear.Model, grid::Radial.Grid;
                             nuclearPot::Union{Nothing,Radial.Potential}=nothing)
    ePot = ParticleScattering.scatteringPotential(ParticleScattering.Electron(), interaction, level, nm, grid;
                                                  nuclearPot=nuclearPot)

    return( Radial.Potential("positron: " * ePot.name, -ePot.Zr, ePot.grid) )
end


"""
`ParticleScattering.scatteringPotential(projectile::ParticleScattering.AbstractProjectile,
                                        interaction::ParticleScattering.AbstractInteractionModel,
                                        level::Level, nm::Nuclear.Model, grid::Radial.Grid;
                                        nuclearPot::Union{Nothing,Radial.Potential}=nothing)`
    ... a fallback that raises for every (projectile, interaction) combination which is not implemented above, and says
        which one was asked for, rather than returning a silently wrong potential. Nothing is returned; this method
        always raises.
"""
function scatteringPotential(projectile::ParticleScattering.AbstractProjectile,
                             interaction::ParticleScattering.AbstractInteractionModel,
                             level::Level, nm::Nuclear.Model, grid::Radial.Grid;
                             nuclearPot::Union{Nothing,Radial.Potential}=nothing)
    error("\n\nNo scattering potential is implemented for the combination\n" *
          "    projectile  = $(projectile)\n    interaction = $(interaction)\n\n" *
          "Implemented today: Electron with StaticField or StaticFieldExchange, and Positron with StaticField. " *
          "A Proton needs the reduced mass in the radial equation, which the continuum solver does not carry yet.")
end
