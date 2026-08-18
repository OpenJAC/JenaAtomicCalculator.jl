
# How the target is represented to the projectile. This is the axis along which ELSEPA is organised, and it is kept
# separate from the projectile and the beam: a further model is one new subtype of AbstractInteractionModel plus one
# new method of ParticleScattering.scatteringPotential below.
#
# The potentials are taken from JAC's own machinery rather than from tabulated Dirac-Fock densities, which is what
# makes the input "typically for JAC": Basics.CHField() is the static field of the target density without exchange,
# and Basics.DFSField(1.0) is that field plus the local Slater exchange term.


"""
`ParticleScattering.scatteringPotential(projectile::ParticleScattering.Electron, interaction::ParticleScattering.StaticField,
                                        level::Level, energy::Float64, nm::Nuclear.Model, grid::Radial.Grid;
                                        nuclearPot::Union{Nothing,Radial.Potential}=nothing)`
    ... to build the local potential in which an ELECTRON partial wave is generated, for the pure static field, i.e. the
        electrostatic interaction with the nuclear and electronic charge density and no exchange term. The nuclear part
        depends only on (nm, grid) and may be passed in by a caller that has built it once. A pot::Radial.Potential is
        returned.
"""
function scatteringPotential(projectile::ParticleScattering.Electron, interaction::ParticleScattering.StaticField,
                             level::Level, energy::Float64, nm::Nuclear.Model, grid::Radial.Grid;
                             nuclearPot::Union{Nothing,Radial.Potential}=nothing)
    nucPot = isnothing(nuclearPot) ? Nuclear.nuclearPotential(nm, grid) : nuclearPot
    wp     = Basics.computePotential(Basics.CHField(), grid, level)

    return( Basics.add(nucPot, wp) )
end


"""
`ParticleScattering.scatteringPotential(projectile::ParticleScattering.Electron, interaction::ParticleScattering.StaticFieldSlaterExchange,
                                        level::Level, energy::Float64, nm::Nuclear.Model, grid::Radial.Grid;
                                        nuclearPot::Union{Nothing,Radial.Potential}=nothing)`
    ... to build the local potential in which an ELECTRON partial wave is generated, for the static field plus a local
        exchange term. The exchange is the Slater term carried by Basics.DFSField(1.0), which is the same field the
        bound-state machinery of JAC uses, so that projectile and target are described consistently.
        A pot::Radial.Potential is returned.
"""
function scatteringPotential(projectile::ParticleScattering.Electron, interaction::ParticleScattering.StaticFieldSlaterExchange,
                             level::Level, energy::Float64, nm::Nuclear.Model, grid::Radial.Grid;
                             nuclearPot::Union{Nothing,Radial.Potential}=nothing)
    nucPot = isnothing(nuclearPot) ? Nuclear.nuclearPotential(nm, grid) : nuclearPot
    wp     = Basics.computePotential(Basics.DFSField(1.0), grid, level)

    return( Basics.add(nucPot, wp) )
end


"""
`ParticleScattering.scatteringPotential(projectile::ParticleScattering.Electron,
                                        interaction::ParticleScattering.StaticFieldFurnessMcCarthy,
                                        level::Level, energy::Float64, nm::Nuclear.Model, grid::Radial.Grid;
                                        nuclearPot::Union{Nothing,Radial.Potential}=nothing)`
    ... to build the local potential in which an ELECTRON partial wave is generated, for the static field plus the
        energy-dependent local exchange potential of Furness and McCarthy, J. Phys. B 6, 2280 (1973),

            V_ex(r) = 1/2 [E - V_st(r)]  -  1/2 sqrt( [E - V_st(r)]^2 + 4 pi rho(r) )        [a.u.]

        with rho(r) the electron density of the target. This is the exchange term ELSEPA and the NIST database use
        (NSRDS-64, Eq. (10)), and it is the one appropriate to a CONTINUUM projectile: it depends on the impact energy
        and weakens as the projectile becomes fast, whereas the Slater term of the DFS field is built for a bound
        electron and does not. V_ex is negative everywhere, so it deepens the well. A pot::Radial.Potential is returned.
"""
function scatteringPotential(projectile::ParticleScattering.Electron,
                             interaction::ParticleScattering.StaticFieldFurnessMcCarthy,
                             level::Level, energy::Float64, nm::Nuclear.Model, grid::Radial.Grid;
                             nuclearPot::Union{Nothing,Radial.Potential}=nothing)
    stPot = ParticleScattering.scatteringPotential(ParticleScattering.Electron(), ParticleScattering.StaticField(),
                                                   level, energy, nm, grid; nuclearPot=nuclearPot)
    # Basics.computeDensity returns SUM_a occ_a (P_a^2 + Q_a^2), i.e. the RADIAL density 4 pi r^2 rho(r), so that the
    # combination 4 pi rho(r) needed below is simply that quantity divided by r^2.
    rDensity = Basics.computeDensity(level, grid)
    newZr    = deepcopy(stPot.Zr)

    for  i = 1:length(newZr)
        r = stPot.grid.r[i]
        r <= 0.  &&  continue
        vSt   = -stPot.Zr[i] / r                       # JAC stores Zr with V(r) = -Zr(r)/r
        fourPiRho = (i <= length(rDensity)) ? rDensity[i] / (r*r) : 0.
        vEx   = 0.5 * (energy - vSt) - 0.5 * sqrt( (energy - vSt)^2 + fourPiRho )
        newZr[i] = newZr[i] - r * vEx                  # V -> V + vEx  means  Zr -> Zr - r vEx
    end

    return( Radial.Potential(stPot.name * " + Furness-McCarthy exchange", newZr, stPot.grid) )
end


"""
`ParticleScattering.scatteringPotential(projectile::ParticleScattering.Positron, interaction::ParticleScattering.StaticField,
                                        level::Level, energy::Float64, nm::Nuclear.Model, grid::Radial.Grid;
                                        nuclearPot::Union{Nothing,Radial.Potential}=nothing)`
    ... to build the local potential in which a POSITRON partial wave is generated. The electrostatic interaction changes
        sign throughout -- the nucleus repels and the electrons attract -- so the electron potential is simply negated;
        no exchange term applies, a positron being distinguishable from the target electrons.
        A pot::Radial.Potential is returned.
"""
function scatteringPotential(projectile::ParticleScattering.Positron, interaction::ParticleScattering.StaticField,
                             level::Level, energy::Float64, nm::Nuclear.Model, grid::Radial.Grid;
                             nuclearPot::Union{Nothing,Radial.Potential}=nothing)
    ePot = ParticleScattering.scatteringPotential(ParticleScattering.Electron(), interaction, level, energy, nm, grid;
                                                  nuclearPot=nuclearPot)

    return( Radial.Potential("positron: " * ePot.name, -ePot.Zr, ePot.grid) )
end


"""
`ParticleScattering.scatteringPotential(projectile::ParticleScattering.AbstractProjectile,
                                        interaction::ParticleScattering.AbstractInteractionModel,
                                        level::Level, energy::Float64, nm::Nuclear.Model, grid::Radial.Grid;
                                        nuclearPot::Union{Nothing,Radial.Potential}=nothing)`
    ... a fallback that raises for every (projectile, interaction) combination which is not implemented above, and says
        which one was asked for, rather than returning a silently wrong potential. Nothing is returned; this method
        always raises.
"""
function scatteringPotential(projectile::ParticleScattering.AbstractProjectile,
                             interaction::ParticleScattering.AbstractInteractionModel,
                             level::Level, energy::Float64, nm::Nuclear.Model, grid::Radial.Grid;
                             nuclearPot::Union{Nothing,Radial.Potential}=nothing)
    error("\n\nNo scattering potential is implemented for the combination\n" *
          "    projectile  = $(projectile)\n    interaction = $(interaction)\n\n" *
          "Implemented today: Electron with StaticField, StaticFieldSlaterExchange or StaticFieldFurnessMcCarthy, and Positron with StaticField. " *
          "A Proton needs the reduced mass in the radial equation, which the continuum solver does not carry yet.")
end
